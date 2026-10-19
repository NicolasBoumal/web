-- Run by Quarto before rendering, or manually with:
--   quarto run scripts/publications.lua
-- No Python, R, browser-side XML request, or separately installed Lua needed.
local xml = dofile('scripts/vendor/slaxml.lua')
local input = 'publications.xml'
local output = 'sections/_publications-generated.qmd'
local site = 'https://www.nicolasboumal.net/'

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

local function url(value)
  if value:match('^[%a][%w+.-]*:') or value:match('^//') or value:match('^#') then
    return value
  end
  if value:match('^downloads/') or value:match('^SIAMOP2023/') then return value end
  return site .. value:gsub('^/', '')
end

local function link(href, label)
  return '<a href="' .. escape(url(href)) .. '" target="_blank" rel="noopener">'
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

local function paper(node, kind)
  local parts = {}
  if kind ~= 'thesis' then
    table.insert(parts, escape(field(node, 'authors')) .. ',<br>')
  end
  table.insert(parts, title(node) .. '<br>')
  table.insert(parts, rich(field(node, 'status')) .. ', ' .. escape(field(node, 'year')) .. '.')
  local comment = field(node, 'comment')
  if comment ~= '' then
    table.insert(parts, '<br><span class="publication-comment">' .. rich(comment) .. '</span>')
  end
  local publisher, bibtex = field(node, 'publisherlink'), field(node, 'bibtex')
  if publisher ~= '' or bibtex ~= '' then table.insert(parts, '<br>') end
  local actions = ''
  if publisher ~= '' then
    actions = '<small>' .. link(publisher, "publisher's link") .. '</small>'
    if bibtex ~= '' then actions = actions .. ', ' end
  end
  if bibtex ~= '' then
    actions = actions .. '<details><summary aria-label="BibTeX for ' .. escape(field(node, 'title'))
      .. '">bibtex</summary><div class="bibtex-content"><pre><code>' .. escape(bibtex)
      .. '</code></pre></div></details>'
  end
  -- Keep inline punctuation together: joining separate fragments with newlines
  -- would introduce a visible space between the publisher link and its comma.
  table.insert(parts, actions)
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
  { 'thesis', 'Theses' }, { 'preprint', 'Reports' },
  { 'journal', 'Journal papers' }, { 'conference', 'Conference papers' },
  { 'talk', 'Talks' }, { 'poster', 'Posters' },
}
local by_type, seen = {}, {}
for _, group in ipairs(groups) do by_type[group[1]] = {} end
local objects = parse(read(input))
for _, object in ipairs(objects) do
  assert(object.name == 'object', 'Expected <object>, got <' .. object.name .. '>')
  local id, kind = object.attributes.id, field(object, 'type')
  assert(id and id:match('^[%w_-]+$'), 'Missing or invalid publication ID: ' .. tostring(id))
  assert(by_type[kind], 'Unknown publication type: ' .. kind)
  assert(field(object, 'title') ~= '', 'Missing publication title: ' .. id)
  seen[id] = (seen[id] or 0) + 1
  object.html_id = id .. (seen[id] > 1 and '-' .. seen[id] or '')
  if seen[id] > 1 then
    io.stderr:write('Warning: duplicate XML ID ' .. id .. '; using publication-' .. object.html_id .. ' in HTML.\n')
  end
  table.insert(by_type[kind], object)
end

local result = { '<!-- Generated from publications.xml by scripts/publications.lua. Do not edit. -->\n' }
for _, group in ipairs(groups) do
  local kind, heading = group[1], group[2]
  table.insert(result, '## ' .. heading .. ' {#publis-' .. kind .. '}\n')
  for _, object in ipairs(by_type[kind]) do
    local body = (kind == 'talk' or kind == 'poster') and presentation(object) or paper(object, kind)
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
print(string.format('Publications: %d entries -> %s', #objects, output))
