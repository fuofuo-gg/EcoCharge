# EcoCharge

This project simulates various pricing strategies for two competing electric vehicle charging stations. Two scenarios are implemented: one where the stations are of equal size and another where they have different sizes. Key parameters such as electricity costs and operational expenses can be adjusted. The results are analyzed using two plots: the first visualizes game theory dynamics, while the second shows the charging stations' profits.

## Prerequisites
Before running the code, make sure you have the following dependencies installed:

git
Julia (version 1.11 or higher) 

Linux:
```bash
winget install julia -s msstore
```

Windows:
https://julialang.org/

All the packages will be installed automatically during the first run.
- GLMakie.jl
- Plots
- StatsPlots
- MAT
- Makie 
- GLMakie  
- ColorSchemes
- GeometryBasics

## How to Use
1. Clone this repository:
```bash
git clone https://github.com/fuofuo-gg/EcoCharge
cd EcoCharge
```
2. Run the Julia script:
```bash
julia
include("./src/ECOCharge.jl")
```

## Example Output
Below is an example of the tool:
![Example Output](images/affichage_makie.png)

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
