using Documenter
using NetworkLayout
using Graphs
using GraphMakie
using CairoMakie
using StableRNGs

NetworkLayout.DEFAULT_RNG[] = StableRNG
DocMeta.setdocmeta!(NetworkLayout, :DocTestSetup, :(using NetworkLayout); recursive=true)

doc = makedocs(; modules=[NetworkLayout],
               repo=Remotes.GitHub("JuliaGraphs", "NetworkLayout.jl"),
               sitename="NetworkLayout.jl",
               build=haskey(ENV, "DOCUMENTER_DRAFT") ? "build_draft" : "build",
               format=Documenter.HTML(; prettyurls=get(ENV, "CI", "false") == "true",
                                      canonical="https://juliagraphs.org/NetworkLayout.jl", assets=String[]),
               pages=["Home" => "index.md",
                      "Interface" => "interface.md"],
               draft=haskey(ENV, "DOCUMENTER_DRAFT"),
               warnonly=true,
               debug=true) # return doc object

# if gh_pages branch gets to big, check out
# https://juliadocs.github.io/Documenter.jl/stable/man/hosting/#gh-pages-Branch

# deploy even if makedocs had errors, so the PR preview is available, but fail CI afterwards
if haskey(ENV, "GITHUB_ACTIONS")
    deploydocs(; repo="github.com/JuliaGraphs/NetworkLayout.jl",
               push_preview=true)
    errors = setdiff(doc.internal.errors, [:missing_docs])
    if !isempty(errors)
        error("makedocs encountered errors: $(errors)")
    end
end
