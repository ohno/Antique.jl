using AnalyticalSolutions
using Weave

for file in AnalyticalSolutions.models # [:InfinitePotentialWell]
  weave("./src/jmd/$file.jmd", doctype="github", out_path="./src/", fig_path="./assets/fig/")
  text = AnalyticalSolutions.load("./src/$file.md")
  # remove ```  after include(...s)
  for m in eachmatch(r"\n```julia\n.*?jl\"\)\n```[.\n]*?```", text)
    @show m.match
    @show m.offset
    text = replace(text, m.match => "")
  end
  # remove ``` at the end of file
  for m in eachmatch(r"```\n```[\s\n]*?\z", text)
    @show m.match
    @show m.offset
    text = replace(text, m.match => "\n```\n")
  end
  AnalyticalSolutions.save("./src/$file.md", text)
end
