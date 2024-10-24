

import Pkg 
Pkg.add("Literate")
using Literate

# Convert the Julia script to a Jupyter notebook
Literate.notebook("Partie_2/Directe.jl", "."; name = "Partie_2/example")
