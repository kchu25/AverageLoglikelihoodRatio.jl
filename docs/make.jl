using AverageLoglikelihoodRatio
using Documenter

DocMeta.setdocmeta!(AverageLoglikelihoodRatio, :DocTestSetup, :(using AverageLoglikelihoodRatio); recursive=true)

makedocs(;
    modules=[AverageLoglikelihoodRatio],
    authors="Shane Kuei-Hsien Chu (skchu@wustl.edu)",
    sitename="AverageLoglikelihoodRatio.jl",
    format=Documenter.HTML(;
        canonical="https://kchu25.github.io/AverageLoglikelihoodRatio.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kchu25/AverageLoglikelihoodRatio.jl",
    devbranch="main",
)
