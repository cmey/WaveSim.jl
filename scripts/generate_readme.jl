# Script to inject code from files into README.template.md to produce README.md

function process_template()
    template_path = "README.template.md"
    readme_path = "README.md"
    
    if !isfile(template_path)
        error("Template file $template_path not found.")
    end

    content = read(template_path, String)
    
    # Regex to find markers like {{INJECT:path/to/file.jl}}
    pattern = r"\{\{INJECT:(.+?)\}\}"
    
    # Use findall to iteratively replace markers manually, as replace(string, function)
    # in Julia 1.x is tricky for regex markers.
    
    matches = collect(eachmatch(pattern, content))
    processed_content = content
    # Process backwards so offsets don't change
    for m in reverse(matches)
        file_to_inject = String(strip(m.captures[1]))
        if isfile(file_to_inject)
            injected_text = read(file_to_inject, String)
            processed_content = string(
                processed_content[1:prevind(processed_content, m.offset)],
                injected_text,
                processed_content[nextind(processed_content, m.offset + length(m.match) - 1):end]
            )
        else
            @warn "File to inject not found: $file_to_inject"
        end
    end

    write(readme_path, processed_content)
    println("Successfully updated $readme_path from $template_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    process_template()
end
