"""
    severusdir()

# Description 

This function simply returns the string for the Severus path. This is helpful for instance to load items, such as meshes, from the `assets`` folder. 
"""
function severusdir()
    joinpath(@__DIR__, "..")
end