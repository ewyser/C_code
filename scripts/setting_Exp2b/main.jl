# Initialisation
using PyPlot, LaTeXStrings, CSV, DataFrames, DelimitedFiles

@views function main()

#mycommand = `./cpu.out`
#run(mycommand)

## load data
x=Matrix(CSV.read("xp.txt"  ,DataFrame,delim=","))
e=Matrix(CSV.read("epII.txt",DataFrame,delim=","))
u=Matrix(CSV.read("up.txt"  ,DataFrame,delim=","))

fig1 = figure("coordinate",figsize=(4.0,2.0)) # Create a figure and save its handle
ax   = axes([0.12;0.2;0.75;0.7])
h    = scatter(x[:,1],x[:,3], s=1.0, c=e, cmap="jet")
axis("equal")
xlim(minimum(x[:,1]),maximum(x[:,1]))
ylim(minimum(x[:,3]),maximum(x[:,3]))
colorbar(h,orientation="horizontal",extend="max")
show()
savefig("epII.png", dpi=600, format="png")


end
main()