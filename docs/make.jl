using Documenter, ArcheoGra, DocStringExtensions

SUBSECTION_ConnectNetwork = [
    "ConNetw/ConNetw_Landing.md", 
    "ConNetw/Data_Types_and_IO.md", 
    "ConNetw/Network_Generation.md",
    "ConNetw/Feature_Extraction_and_Comparison.md",
    "ConNetw/Clustering_and_GNN.md",
    "ConNetw/Quality_Control.md",
    "ConNetw/Visualization.md"
]

SUBSECTION_LCPA = [
    "LCPA/LCPA_Landing.md",
    "LCPA/Data_Management_and_IO.md",
    "LCPA/Cost_and_Adjacency.md",
    "LCPA/Path_Analysis.md",
    "LCPA/Wavelet_Compression.md",
    "LCPA/Adaptive_Refinement.md"
]

# Placeholder array for future example and common function files
SUBSECTION_ArcheoGra_CommonFuntions = [
    "Examples.md",
    "CommonFunction.md"
]

makedocs(
    sitename = "ArcheoGra",
    format = Documenter.HTML(
    prettyurls = get(ENV, "CI", "false") == "true"
    ),
    remotes = nothing,
    modules = [ArcheoGra],
    checkdocs = :exports,
    pages = [
        "index.md",
        "ArcheoGra - Examples" => SUBSECTION_ArcheoGra_CommonFuntions,
        "Connect Archaeological Network" => SUBSECTION_ConnectNetwork,
        "LCPA and Wavelets" => SUBSECTION_LCPA
        # "Developer API" => "Developer_API.md"
    ]
)

deploydocs(
    repo = "github.com/MaxBrockm/ArcheoGra.jl.git",
    devbranch = "main"
)