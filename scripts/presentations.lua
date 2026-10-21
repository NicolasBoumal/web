-- Presentation validation and HTML generation; called by publications.lua.
return function(helpers)
  local rich, link, escape = helpers.rich, helpers.link, helpers.escape
  local resources = { 'slides', 'poster', 'video', 'notes' }
  local months = { 'Jan.', 'Feb.', 'March', 'April', 'May', 'June',
    'July', 'Aug.', 'Sep.', 'Oct.', 'Nov.', 'Dec.' }
  local function text(value)
    return type(value) == 'string' and value:match('%S') ~= nil
  end
  local function date_parts(date)
    if type(date) ~= 'string' then return end
    local year, month, day = date:match('^(%d%d%d%d)%-(%d%d)%-(%d%d)$')
    if not year then year, month = date:match('^(%d%d%d%d)%-(%d%d)$') end
    year, month, day = tonumber(year), tonumber(month), tonumber(day)
    if not year or year < 1 or not months[month] then return end
    local leap = year % 4 == 0 and (year % 100 ~= 0 or year % 400 == 0)
    local days = {31, leap and 29 or 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}
    if day and (day < 1 or day > days[month]) then return end
    return year, month, day
  end
  -- Stable ordering: missing days compare as the first day of the month.
  local function newest_first(list, key)
    local positions = {}
    for i, value in ipairs(list) do positions[value] = i end
    table.sort(list, function(a, b)
      local x, y = key(a), key(b)
      if #x == 7 then x = x .. '-01' end
      if #y == 7 then y = y .. '-01' end
      if x == y then return positions[a] < positions[b] end
      return x > y
    end)
  end
  local function validate_fields(object, context)
    for _, name in ipairs({ 'comment', 'abstract', 'meta' }) do
      assert(object[name] == nil or text(object[name]), 'Invalid ' .. context .. '.' .. name)
    end
    for _, name in ipairs(resources) do
      local value = object[name]
      if value ~= nil and value ~= false then
        if name == 'poster' and value == true then
          -- Missing poster: display an unlinked label.
        elseif name == 'video' and type(value) == 'table' then
          assert(text(value.link), 'Missing ' .. context .. '.video.link')
          assert(value.comment == nil or text(value.comment), 'Invalid ' .. context .. '.video.comment')
        else
          assert(text(value), 'Invalid ' .. context .. '.' .. name)
        end
      end
    end
  end
  local function validate(record)
    local context = 'presentation ' .. tostring(record.id)
    assert(type(record.events) == 'table' and #record.events > 0, 'Missing events for ' .. context)
    validate_fields(record, context)
    assert(record.link == nil, 'Use notes/resources or event.link, not a title link: ' .. context)
    assert(record.webpage == nil, 'Use event.link instead of webpage: ' .. context)
    if record.coauthors then
      assert(type(record.coauthors) == 'table' and #record.coauthors > 0, 'Invalid coauthors for ' .. context)
      for _, author in ipairs(record.coauthors) do assert(text(author), 'Invalid coauthor for ' .. context) end
    end
    for i, event in ipairs(record.events) do
      local where = context .. '.events[' .. i .. ']'
      assert(type(event) == 'table', 'Invalid ' .. where)
      assert(event.webpage == nil, 'Use event.link instead of webpage: ' .. where)
      assert(date_parts(event.date), 'Invalid date for ' .. where)
      assert(text(event.venue), 'Missing venue for ' .. where)
      for _, name in ipairs({ 'location', 'link' }) do
        assert(event[name] == nil or text(event[name]), 'Invalid ' .. where .. '.' .. name)
      end
      assert(event.online == nil or type(event.online) == 'boolean', 'Invalid ' .. where .. '.online')
      validate_fields(event, where)
    end
    newest_first(record.events, function(event) return event.date end)
  end
  local function actions(object)
    local parts = {}
    for _, name in ipairs(resources) do
      local value = object[name]
      if value then
        local label = name
        local rendered
        if name == 'poster' and value == true then rendered = 'poster'
        elseif name == 'video' and type(value) == 'table' then
          rendered = link(value.link, label)
          if value.comment then rendered = rendered .. ' (' .. rich(value.comment) .. ')' end
        else rendered = link(value, label) end
        table.insert(parts, rendered)
      end
    end
    return table.concat(parts, ', ')
  end
  local function comment(value)
    return '<div class="publication-comment">' .. pandoc.write(pandoc.read(value, 'markdown'), 'html',
      { html_math_method = 'mathjax', wrap_text = 'none' }) .. '</div>'
  end
  local function render(record)
    local parts = { '<h3>' .. rich(record.title) .. '</h3>' }
    if record.coauthors then
      local names = {}
      for _, author in ipairs(record.coauthors) do table.insert(names, escape(author)) end
      local authors = #names == 1 and names[1]
        or table.concat(names, ', ', 1, #names - 1) .. ' and ' .. names[#names]
      table.insert(parts, '<br><em>(with ' .. authors .. ')</em>')
    end
    if record.comment then table.insert(parts, comment(record.comment)) end
    if #record.events > 1 then
      local shared = actions(record)
      if shared ~= '' then table.insert(parts, '<div class="presentation-resources">' .. shared .. '</div>') end
    end
    table.insert(parts, '<ul class="eventlist">')
    for _, event in ipairs(record.events) do
      local year, month, day = date_parts(event.date)
      local date = string.format('%04d, %s', year, months[month]) .. (day and (' ' .. day) or '')
      local venue = rich(event.venue)
      if event.link then venue = link(event.link, venue) end
      if event.location then venue = venue .. ', ' .. escape(event.location) end
      if event.online then venue = venue .. ' (online)' end
      local text = date .. ' · ' .. venue
      local resolved = {}
      for _, name in ipairs(resources) do
        if event[name] ~= nil then resolved[name] = event[name]
        elseif #record.events == 1 then resolved[name] = record[name] end
      end
      local links = actions(resolved)
      if links ~= '' then text = text .. ' · ' .. links end
      if event.comment then text = text .. comment(event.comment) end
      table.insert(parts, '<li>' .. text .. '</li>')
    end
    table.insert(parts, '</ul>')
    -- Abstracts are preserved as metadata; displaying them is a separate choice.
    return table.concat(parts, '\n')
  end
  return { validate = validate, render = render,
    sort = function(records) newest_first(records, function(record) return record.events[1].date end) end }
end
