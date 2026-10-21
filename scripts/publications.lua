-- Run by Quarto before rendering, or manually with:
--   quarto run scripts/publications.lua
-- No Python, R, browser-side XML request, or separately installed Lua needed.
local xml = dofile('scripts/vendor/slaxml.lua')
local yaml = dofile('scripts/vendor/tinyyaml.lua')
local input = 'publications.xml' -- Talks and posters only.
local publications_input = 'publications.yml'
local output = 'sections/_publications-generated.qmd'

local function read(path)
  local file, err = io.open(path, 'rb')
  if not file then error(err) end
  local content = file:read('*a')
  file:close()
  return content
end

local function trim(text)
  return (text:gsub('^%s+', ''):gsub('%s+$', ''))
end

local function escape(text)
  return (text:gsub('&', '&amp;'):gsub('<', '&lt;'):gsub('>', '&gt;')
    :gsub('"', '&quot;'):gsub("'", '&#39;'))
end

local function link(href, label)
  return '<a href="' .. escape(href) .. '" target="_blank" rel="noopener">'
    .. label .. '</a>'
end

-- A small DOM over SLAXML's SAX callbacks. Validate nesting here so malformed
-- XML fails the build instead of silently publishing an incomplete list.
local function parse(source)
  local root = { name = '#document', children = {}, text = '' }
  local stack = { root }
  xml:parser({
    startElement = function(name)
      local parent = stack[#stack]
      local node = { name = name, attributes = {}, children = {}, text = '' }
      table.insert(parent.children, node)
      table.insert(stack, node)
    end,
    attribute = function(name, value)
      stack[#stack].attributes[name] = value
    end,
    text = function(text)
      local node = stack[#stack]
      node.text = node.text .. text
    end,
    closeElement = function(name)
      assert(#stack > 1 and stack[#stack].name == name,
        'Mismatched closing XML element: ' .. name)
      table.remove(stack)
    end,
  }):parse(source)
  assert(#stack == 1, 'Unclosed XML element in ' .. input)
  assert(#root.children == 1 and root.children[1].name == 'publications',
    'Expected one <publications> root in ' .. input)
  return root.children[1].children
end

local function field(node, name)
  for _, child in ipairs(node.children) do
    if child.name == name then return trim(child.text) end
  end
  return ''
end

-- Comments and venues in the existing XML intentionally contain escaped HTML.
-- Pandoc also handles their occasional Markdown/math before writing HTML.
local function rich(text)
  if text == '' then return '' end
  local html = pandoc.write(pandoc.read(text, 'markdown'), 'html',
    { html_math_method = 'mathjax', wrap_text = 'none' })
  html = trim(html)
  return (html:gsub('^<p>', ''):gsub('</p>$', ''))
end

local function title(node)
  local text = escape(field(node, 'title'))
  local href = field(node, 'link')
  if href ~= '' then text = link(href, text) end
  return '<h3>' .. text .. ',</h3>'
end

local icons = {
  code = '<svg class="publication-icon" viewBox="0 0 16 16" aria-hidden="true" focusable="false"><path d="m5 4-4 4 4 4m6-8 4 4-4 4M9 2 7 14"/></svg>',
  video = '<svg class="publication-icon" viewBox="0 0 16 16" aria-hidden="true" focusable="false"><path d="m4 2 10 6-10 6Z"/></svg>',
}

local function paper(record)
  local parts = {}
  if record.type ~= 'thesis' then
    local authors = {}
    for i, author in ipairs(record.authors) do authors[i] = escape(author) end
    local text = table.concat(authors, ', ', 1, math.max(1, #authors - 1))
    if #authors > 1 then text = text .. ' and ' .. authors[#authors] end
    table.insert(parts, text .. ',<br>')
  end
  local publication = record.publication
  local href = record.link
    or (record.arxiv and ('https://arxiv.org/abs/' .. record.arxiv))
    or (publication and publication.link)
  local heading = rich(record.title)
  if href then heading = link(href, heading) end
  table.insert(parts, '<h3>' .. heading .. ',</h3><br>')
  if publication then
    local venue = rich(publication.venue)
    if publication.comment then venue = venue .. ' (' .. rich(publication.comment) .. ')' end
    table.insert(parts, venue .. ', ' .. publication.year .. '.')
  else
    local year = 2000 + tonumber(record.arxiv:sub(1, 2))
    table.insert(parts, (record.type == 'note' and 'Notes' or 'Preprint') .. ', ' .. year .. '.')
  end
  if record.comment then
    local html = pandoc.write(pandoc.read(record.comment, 'markdown'), 'html',
      { html_math_method = 'mathjax', wrap_text = 'none' })
    table.insert(parts, '<div class="publication-comment">' .. trim(html) .. '</div>')
  end
  local actions = {}
  if publication and publication.link then
    table.insert(actions, link(publication.link, "publisher's link"))
  end
  if record.code then table.insert(actions, link(record.code, icons.code .. 'code')) end
  if record.video then
    local video = type(record.video) == 'string' and { link = record.video } or record.video
    local text = link(video.link, icons.video .. 'video')
    if video.comment then text = text .. ' (' .. rich(video.comment) .. ')' end
    table.insert(actions, text)
  end
  if record.bibtex then
    table.insert(actions, '<details><summary aria-label="BibTeX for ' .. escape(record.title)
      .. '">bibtex</summary><div class="bibtex-content"><pre><code>' .. escape(record.bibtex)
      .. '</code></pre></div></details>')
  end
  if #actions > 0 then
    -- Keep separators in the same text node, with no whitespace before commas.
    table.insert(parts, '<div class="publication-resources">' .. table.concat(actions, ', ') .. '</div>')
  end
  return table.concat(parts, '\n')
end

local months = { 'Jan.', 'Feb.', 'March', 'April', 'May', 'June',
  'July', 'Aug.', 'Sep.', 'Oct.', 'Nov.', 'Dec.' }

local function presentation(node)
  local parts = { title(node) }
  local coauthors = field(node, 'with')
  if coauthors ~= '' then
    table.insert(parts, '<br><em>(with ' .. escape(coauthors) .. ')</em>')
  end
  table.insert(parts, '<ul class="eventlist">')
  for _, event in ipairs(node.children) do
    if event.name == 'event' then
      local month = tonumber(field(event, 'month'))
      assert(month and months[month], 'Invalid event month in publication ' .. node.attributes.id)
      local day = field(event, 'day')
      local date = months[month] .. (day ~= '' and ' ' .. escape(day) or '')
        .. ', ' .. escape(field(event, 'year'))
      local text = rich(field(event, 'venue')) .. ' in ' .. escape(field(event, 'city'))
        .. ' (' .. escape(field(event, 'country')) .. ') on ' .. date .. '.'
      local video = field(event, 'video')
      if video ~= '' then text = text .. ' ' .. link(video, 'video') end
      table.insert(parts, '<li>' .. text .. '</li>')
    end
  end
  table.insert(parts, '</ul>')
  return table.concat(parts, '\n')
end

local groups = {
  { 'preprint', 'Reports' },
  { 'journal', 'Journal papers' }, { 'conference', 'Conference papers' },
  { 'note', 'Notes' }, { 'thesis', 'Theses' },
  { 'talk', 'Talks' }, { 'poster', 'Posters' },
}
local by_type, seen = {}, {}
for _, group in ipairs(groups) do by_type[group[1]] = {} end

local function register(id, kind, title, object)
  assert(type(id) == 'string' and id:match('^[%w_-]+$'), 'Missing or invalid publication ID: ' .. tostring(id))
  assert(not seen[id], 'Duplicate publication ID: ' .. id)
  assert(by_type[kind], 'Unknown publication type: ' .. tostring(kind))
  assert(type(title) == 'string' and title ~= '', 'Missing publication title: ' .. id)
  seen[id] = true
  object.html_id = id
  table.insert(by_type[kind], object)
end

local function optional_text(record, key, context)
  if record[key] ~= nil then
    assert(type(record[key]) == 'string' and trim(record[key]) ~= '', 'Invalid ' .. context .. '.' .. key)
  end
end

local records = yaml.parse(read(publications_input))
assert(type(records) == 'table' and #records > 0, 'Expected a list in ' .. publications_input)
for _, record in ipairs(records) do
  local id = tostring(record.id)
  assert(record.type ~= 'talk' and record.type ~= 'poster', 'Presentations must remain in ' .. input)
  assert(type(record.authors) == 'table' and #record.authors > 0, 'Missing authors for ' .. id)
  for _, author in ipairs(record.authors) do
    assert(type(author) == 'string' and trim(author) ~= '', 'Invalid author for ' .. id)
  end
  for _, key in ipairs({ 'link', 'arxiv', 'code', 'comment', 'bibtex' }) do
    optional_text(record, key, id)
  end
  if record.arxiv then
    assert(record.arxiv:match('^%d%d%d%d%.%d%d%d%d%d?$'), 'Invalid arXiv ID for ' .. id)
  end
  if record.publication then
    local publication = record.publication
    assert(type(publication) == 'table', 'Invalid publication metadata for ' .. id)
    assert(type(publication.year) == 'number' and publication.year % 1 == 0, 'Missing publication year for ' .. id)
    assert(type(publication.venue) == 'string' and publication.venue ~= '', 'Missing publication venue for ' .. id)
    optional_text(publication, 'link', id .. '.publication')
    optional_text(publication, 'comment', id .. '.publication')
  else
    assert(record.arxiv, 'A publication year or arXiv ID is required for ' .. id)
  end
  if record.video then
    local video = type(record.video) == 'string' and { link = record.video } or record.video
    assert(type(video) == 'table' and type(video.link) == 'string' and video.link ~= '', 'Missing video link for ' .. id)
    optional_text(video, 'comment', id .. '.video')
  end
  register(record.id, record.type, record.title, record)
end

local presentations = parse(read(input))
for _, object in ipairs(presentations) do
  assert(object.name == 'object', 'Expected <object>, got <' .. object.name .. '>')
  local id, kind = object.attributes.id, field(object, 'type')
  assert(kind == 'talk' or kind == 'poster', 'Papers and theses must be maintained in ' .. publications_input)
  register(id, kind, field(object, 'title'), object)
end

local result = { '<!-- Generated from publications.yml and publications.xml by scripts/publications.lua. Do not edit. -->\n' }
for _, group in ipairs(groups) do
  local kind, heading = group[1], group[2]
  table.insert(result, '## ' .. heading .. ' {#publis-' .. kind .. '}\n')
  for _, object in ipairs(by_type[kind]) do
    local body = (kind == 'talk' or kind == 'poster') and presentation(object) or paper(object)
    table.insert(result, '```{=html}\n<article class="publication" id="publication-'
      .. object.html_id .. '" data-type="' .. kind .. '">\n' .. body .. '\n</article>\n```\n')
  end
end

local content = table.concat(result, '\n')
local existing = io.open(output, 'rb')
local previous = existing and existing:read('*a')
if existing then existing:close() end
if previous ~= content then
  local file = assert(io.open(output, 'wb'))
  file:write(content)
  file:close()
end
print(string.format('Publications: %d YAML entries and %d XML presentations -> %s', #records, #presentations, output))
