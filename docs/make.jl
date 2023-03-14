using GapFilling
using Documenter

DocMeta.setdocmeta!(GapFilling, :DocTestSetup, :(using GapFilling); recursive=true)

makedocs(;
    modules=[GapFilling],
    authors="Lokmen Farhat",
    repo="https://github.com/farhatlokmen/GapFilling.jl/blob/{commit}{path}#{line}",
    sitename="GapFilling.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://farhatlokmen.github.io/GapFilling.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/farhatlokmen/GapFilling.jl",
    devbranch="main",
)
