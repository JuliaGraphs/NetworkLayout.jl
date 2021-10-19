using Documenter
using NetworkLayout
using Graphs
using GraphMakie
using CairoMakie

DocMeta.setdocmeta!(NetworkLayout, :DocTestSetup, :(using NetworkLayout); recursive=true)

makedocs(; modules=[NetworkLayout],
         repo="https://github.com/JuliaGraphs/NetworkLayout.jl/blob/{commit}{path}#{line}",
         sitename="NetworkLayout.jl",
         format=Documenter.HTML(; prettyurls=get(ENV, "CI", "false") == "true",
                                canonical="https://juliagraphs.org/NetworkLayout.jl", assets=String[]),
         pages=["Home" => "index.md",
                "Interface" => "interface.md"])

# if gh_pages branch gets to big, check out
# https://juliadocs.github.io/Documenter.jl/stable/man/hosting/#gh-pages-Branch
deploydocs(;repo="github.com/JuliaGraphs/NetworkLayout.jl",
           push_preview=true)
