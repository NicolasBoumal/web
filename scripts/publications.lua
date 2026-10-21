-- Run by Quarto before rendering, or manually with:
--   quarto run scripts/publications.lua
-- No Python, R, browser-side XML request, or separately installed Lua needed.
local yaml = dofile('scripts/vendor/tinyyaml.lua')
local presentations_input = 'presentations.yml'
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

-- Render inline Markdown/math in YAML text fields.
local function rich(text)
  if text == '' then return '' end
  local html = pandoc.write(pandoc.read(text, 'markdown'), 'html',
    { html_math_method = 'mathjax', wrap_text = 'none' })
  html = trim(html)
  return (html:gsub('^<p>', ''):gsub('</p>$', ''))
end

local icons = {
  code = '<svg class="publication-icon" viewBox="0 0 16 16" aria-hidden="true" focusable="false"><path d="m5 4-4 4 4 4m6-8 4 4-4 4M9 2 7 14"/></svg>',
  video = '<svg class="publication-icon" viewBox="0 0 16 16" aria-hidden="true" focusable="false"><path d="m4 2 10 6-10 6Z"/></svg>',
}

local function paper(record)
  local parts = {}
  local publication = record.publication
  local year = publication and publication.year or (2000 + tonumber(record.arxiv:sub(1, 2)))
  if record.type ~= 'thesis' then
    local authors = {}
    for i, author in ipairs(record.authors) do authors[i] = escape(author) end
    local text = table.concat(authors, ', ', 1, math.max(1, #authors - 1))
    if #authors > 1 then text = text .. ' and ' .. authors[#authors] end
    table.insert(parts, text .. (record.type == 'note' and (', ' .. year .. '.') or ',') .. '<br>')
  end
  local href = record.link
    or (record.arxiv and ('https://arxiv.org/abs/' .. record.arxiv))
    or (publication and publication.link)
  local heading = rich(record.title)
  if href then heading = link(href, heading) end
  table.insert(parts, '<h3>' .. heading .. (record.type == 'note' and '.' or ',') .. '</h3><br>')
  if publication then
    local venue = rich(publication.venue)
    if publication.comment then venue = venue .. ' (' .. rich(publication.comment) .. ')' end
    table.insert(parts, venue .. (record.type == 'note' and '' or (', ' .. year)) .. '.')
  elseif record.type ~= 'note' then
    table.insert(parts, 'Preprint, ' .. year .. '.')
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

local presentation = dofile('scripts/presentations.lua')({ rich = rich, link = link, escape = escape })

local groups = {
  { 'preprint', 'Reports', 'reports' },
  { 'journal', 'Journal papers', 'journals' }, { 'conference', 'Conference papers', 'conferences' },
  { 'note', 'Notes', 'notes' }, { 'thesis', 'Theses', 'theses' },
  { 'presentation', 'Presentations', 'talks' },
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
  assert(record.type ~= 'presentation', 'Use ' .. presentations_input .. ' for presentations')
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

local presentations = yaml.parse(read(presentations_input))
assert(type(presentations) == 'table' and #presentations > 0, 'Expected a list in ' .. presentations_input)
for _, record in ipairs(presentations) do
  presentation.validate(record)
  register(record.id, 'presentation', record.title, record)
end
presentation.sort(by_type.presentation)

local result = { '<!-- Generated from publications.yml and presentations.yml by scripts/publications.lua. Do not edit. -->\n' }
local section_links = {}
for _, group in ipairs(groups) do
  table.insert(section_links, '<a href="#pubs-' .. group[3] .. '">' .. escape(group[2]) .. '</a>')
end
table.insert(result, '```{=html}\n<nav class="publication-sections" aria-label="Publication sections">'
  .. table.concat(section_links, ' <span aria-hidden="true">·</span> ') .. '</nav>\n```\n')
for _, group in ipairs(groups) do
  local kind, heading = group[1], group[2]
  table.insert(result, '## ' .. heading .. ' {#pubs-' .. group[3] .. '}\n')
  for _, object in ipairs(by_type[kind]) do
    local body = kind == 'presentation' and presentation.render(object) or paper(object)
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
print(string.format('Publications: %d papers/theses and %d presentation groups -> %s', #records, #presentations, output))
