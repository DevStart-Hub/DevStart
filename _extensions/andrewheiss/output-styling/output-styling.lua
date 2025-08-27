function Div(div)
  -- Check if this div has warning output classes
  -- R/Python warnings: cell-output-stderr
  -- Julia warnings: cell-output-stdout
  if div.classes:includes('cell-output-stderr') or div.classes:includes('cell-output-stdout') then
    -- Check if any CodeBlock or RawBlock in this div contains warning patterns
    for _, block in ipairs(div.content) do
      local should_add_warning = false
      
      if block.tag == "CodeBlock" then
        -- R/Python warnings in CodeBlocks (stderr)
        if div.classes:includes('cell-output-stderr') and 
           string.match(block.text, "^%w*Warning:") then
          should_add_warning = true
        end
      elseif block.tag == "RawBlock" and block.format == "html" then
        -- Julia warnings in HTML RawBlocks (stdout)
        if div.classes:includes('cell-output-stdout') and 
           string.match(block.text, "Warning: ") then
          should_add_warning = true
        end
      elseif block.tag == "Div" and block.classes:includes('ansi-escaped-output') then
        -- Julia warnings are nested in ansi-escaped-output divs
        if div.classes:includes('cell-output-stdout') then
          for _, nested_block in ipairs(block.content) do
            if nested_block.tag == "RawBlock" and nested_block.format == "html" then
              if string.match(nested_block.text, "Warning: ") then
                should_add_warning = true
                break
              end
            end
          end
        end
      end

      if should_add_warning then
        div.classes:insert('cell-output-warning')
        break
      end
    end
  end

  return div
end

-- Add the CSS file to the document
function Pandoc(doc)
  if not (quarto.doc and quarto.doc.is_format and quarto.doc.is_format("html")) then
    return doc
  end

  -- Get value of `appearance` key from YAML fron matter
  local appearance_value = "default"

  if doc.meta['output-styling'] and doc.meta['output-styling'].appearance then
    local raw_appearance = pandoc.utils.stringify(doc.meta['output-styling'].appearance)

    -- Possible values here
    -- MAYBE: someday add others?
    local valid_appearances = {default = true, minimal = true, custom = true}

    if valid_appearances[raw_appearance] then
      appearance_value = raw_appearance
    else
      quarto.log.warning(
        "Invalid value for `appearance` '" .. raw_appearance .. 
        "'. Using 'default'. Valid options are: `default`, `minimal`, or `custom`"
      )
    end
  end

  -- Include CSS file only if not custom
  if appearance_value ~= "custom" then
    local css_file = "output-styling-" .. appearance_value .. ".css"

    quarto.doc.add_html_dependency({
      name = 'output-styling',
      stylesheets = {css_file}
    })
  end

  return doc
end
