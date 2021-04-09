using Documenter, NormalHermiteSplines

makedocs(
    sitename = "NormalHermiteSplines.jl",
	format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
	authors = "Igor Kohanovsky",
    pages = [
				"Home" => "index.md",
				"Public API" => "Public-API.md",
				"Example Usage" => "Usage.md",
				"Selecting a good value of the scaling parameter" => "Parameter-Choice.md",
				"Numerical Tests" => "Numerical-Tests.md",
				"Normal Splines Method" => "Normal-Splines-Method.md",
			]
)

deploydocs(
    repo = "github.com/IgorKohan/NormalHermiteSplines.jl.git",
	devurl = "dev",
    versions = ["stable" => "v^", "v#.#"],
)

