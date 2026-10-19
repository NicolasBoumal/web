-- Expand after pre-render, rather than using an include shortcode that Quarto
-- resolves during project discovery (before the generated file exists).
function Div(div)
  if not div.classes:includes('publication-list') then return nil end
  local file = assert(io.open('sections/_publications-generated.qmd', 'rb'))
  local markdown = file:read('*a')
  file:close()
  return pandoc.read(markdown, 'markdown').blocks
end
