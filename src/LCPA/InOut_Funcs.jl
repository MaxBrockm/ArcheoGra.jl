"""
    write_Output(vec::Vector, desired_filename::String)

Writes the elements of a vector line by line into a text file.

# Arguments
- `vec::Vector`: The input vector whose elements will be written to the file.
- `desired_filename::String`: The name or path of the target file (e.g., "output.txt").

# Output
- `nothing`: The function writes directly to the disk and returns nothing.
"""
function write_Output(vec::Vector, desired_filename::String)
    open(desired_filename, "w") do file  
        for l in 1:length(vec)
            write(file, "$(vec[l])\n")
        end
    end
end


"""
    tiff_to_matrix(filename::String)

Loads a `.tiff` image file and converts it into a 2D matrix of `Float64` values.

# Arguments
- `filename::String`: The file path to the `.tiff` image.

# Output
- `Matrix{Float64}`: A 2D matrix representing the image data.

# Examples
```julia
img = tiff_to_matrix("example.tiff")
size(img)  # → (1024, 1024)
```
"""
function tiff_to_matrix(filename::String)
    img = load(filename)
    img_array = convert(Array{Float64}, img)
    return img_array
end