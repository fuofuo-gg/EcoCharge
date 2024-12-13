import Pkg
Pkg.add("Plots")
Pkg.add("StatsPlots")
Pkg.add("MAT")
Pkg.add("Makie")
Pkg.add("GLMakie")
Pkg.add("ColorSchemes")
Pkg.add("GeometryBasics")

using MAT 
using Makie
using GLMakie
using ColorSchemes
using GeometryBasics
using Plots
using StatsPlots

function Neutral(z1, z2)
    z = zeros(length(z1[:,1]), length(z2[:,1]))
    for ii in 1:length(x)
        for jj in 1:length(y) 
            z[ii,jj] = z1[ii,jj] + z2[ii,jj] #- z1[8,8] - z2[8,8]
        end
    end
    return z 
end
function MinimaxA(z1, z2)
    z = zeros(length(z1[:,1]), length(z2[:,1]))
    for ii in 1:length(z1[:,1])
        z[ii,argmin(z1[ii, :])] = minimum(z1[ii, :])
    end
    return z 
end
function MinimaxB(z1, z2)
    z = zeros(length(z1[:,1]), length(z2[:,1]))
    for ii in 1:length(z2[:,1])
        z[argmin(z2[:,ii]),ii] = minimum(z2[:,ii])
    end
    return z 
end
function SharpeA(z1, z2)
    pi_A_star = argmax(MinimaxA(z1,z2))[1]
    pi_B_star = argmax(MinimaxB(z1,z2))[2]
    z = zeros(length(z1[:,1]), length(z2[:,1]))
    for ii in 1:length(z1[:,1])
        xA_bar = 1/length(z2[:,1]) * sum(z1[ii,:])
        sigma = sqrt(1/length(z2[:,1]) * sum((z1[ii,:] .- xA_bar).^2))
        z[ii,:] .= (xA_bar - z1[pi_A_star, pi_B_star] ) / sigma
    end
    return z
end
function SharpeB(z1, z2)
    pi_A_star = argmax(MinimaxA(z1,z2))[1]
    pi_B_star = argmax(MinimaxB(z1,z2))[2]
    z = zeros(length(z1[:,1]), length(z2[:,1]))
    for ii in 1:length(z2[1,:])
        xB_bar = 1/length(z2[:,1]) * sum(z2[:,ii])
        sigma = sqrt(1/length(z1[:,1]) * sum((z2[:,ii] .- xB_bar).^2))
        z[:,ii] .= (xB_bar - z2[pi_A_star, pi_B_star] ) / sigma
    end
    return z
end
function StackelbergA(z1, z2)
    z = zeros(length(z1[:,1]), length(z2[:,1]))
    for ii in 1:length(z1[1,:])
        z[ii,:] .= z1[ii, argmax(z2[ii,:])]
    end
    return z
end
function StackelbergB(z1, z2)
    z = zeros(length(z1[:,1]), length(z2[:,1]))
    for ii in 1:length(z2[1,:])
        z[:,ii] .= z2[argmax(z1[:,ii]), ii]
    end
    return z
end

function Shapley_value(z1,z2)
    z = argmax(MinimaxA(z1, z2))
    posA = z[1]
    z = argmax(MinimaxB(z1, z2))
    posB = z[2]
    JA = z1[posA, posB]
    JB = z2[posA, posB]
    z1z2 = z1 .+ z2
    A2 = (z1z2 .+ JA .- JB)./2
    B2 = (z1z2 .+ JB .- JA)./2
    #z = A2/z1 .+ B2/z2
    z = @. ifelse(z1 > 0 && z2 > 0, A2/JA * B2/JB, 0)
    #min_value = minimum(filter(!isinf, z))
    #z .= [isinf(v) ? 0 : (v < -10 ? 0 : v) for v in z]
    return z
end

function setPriceMatrices(J, Bid, Cost, Elec)
    Price1= zeros(length(J[][1,:,1]), length(J[][1,:,2]))
    Price2= zeros(length(J[][1,:,1]), length(J[][1,:,2]))
    for ii in 1:length(J[][1,:,1])
        Price1[ii,:] .= 0.2 + (0.02 * (ii-1)) .- Cost[] .- Elec[]
        Price2[:,ii] .= 0.2 + (0.02 * (ii-1)) .- Cost[] .- Elec[]
    end
    
    z1 = J[][:,:,1] 
    z2 = J[][:,:,2] 
    z1 = (z1 .* Price1  .+ (Bid[][:,:,1] .* 100)) 
    z2 = (z2 .* Price2  .+ (Bid[][:,:,2] .* 100)) 
    return(z1, z2)
end

function color_map(z, cmap)
    return  cmap[Int(round(z, digits=3)*1000)+1]
end

function plot_triangles(x, y, z_normalized, ax, cmap, reverse=false)
    mesh_matrix = Matrix{Ref{Makie.Mesh}}(undef, length(x), length(x))
    dist = x[2] - x[1]
    pos = dist/2
    for i in 1:length(x)
        for j in 1:length(y)
            # Coordonnées des sommets
            xs = reverse ? [(x[i]+dist-pos), (x[i]-pos), (x[i]-pos)] : [(x[i]-pos), (x[i]+dist-pos), (x[i]+dist-pos)]
            ys = reverse ? [(y[j]+dist-pos), (y[j]+dist-pos), (y[j]-pos)] : [(y[j]-pos), (y[j]-pos), (y[j]+dist-pos)]

            vertices = Point2f0.([xs[k], ys[k]] for k in 1:3)
            
            # Appliquer la couleur
            color_vals = color_map(z_normalized[i,j], cmap)
            if reverse == true
                mesh = Makie.mesh!(ax, vertices, color=color_vals, shading=NoShading)
                mesh_matrix[i,j] = mesh
            else
                mesh = Makie.mesh!(ax, vertices, color=color_vals, shading=NoShading)
                mesh_matrix[i,j] = mesh
            end
        end
    end
    return mesh_matrix
end

function update_triangles(x, y, z_normalized, cmap, reverse=false)
    for i in 1:length(x)
        for j in 1:length(y)
            color_vals = color_map(z_normalized[i,j], cmap)
            if reverse == true
                Save_triangle_plot_true[i,j][].color = color_vals
                Save_triangle_plot_true[i,j][].alpha = 1.0
            else
                Save_triangle_plot_false[i,j][].color = color_vals
                Save_triangle_plot_false[i,j][].alpha = 1.0
            end
        end
    end
end

function add_game_results(z1, z2, Minimax, Sharpe, Stackelberg, Shapley, win_win)
    posAminimax = 0
    posBminimax = 0
    posAsharpe = 0 
    posBsharpe = 0
    posAstackelberg = 0
    posBstackelberg = 0
    posAwinwin = 0
    posBwinwin = 0
    if Minimax == true
        z = argmax(MinimaxA(z1, z2))
        posAminimax = z[1]
        z = argmax(MinimaxB(z1, z2))
        posBminimax = z[2]
        Save_triangle_plot_true[posAminimax,posBminimax][].color = "purple"
        Save_triangle_plot_false[posAminimax,posBminimax][].color = "purple"
        val_label_tpA.text[] = string("Earnings A :", round(Int, z1[posAminimax, posBminimax]))
        val_label_tpB.text[] = string("B :", round(Int, z2[posAminimax, posBminimax]))
    end

    if Sharpe == true
        posAsharpe = argmax(StackelbergA(z1, z2)[:,1])
        posBsharpe = argmax(StackelbergB(z1, z2)[1,:])
        Save_triangle_plot_true[posAsharpe,posBsharpe][].color = "violetred2"
        Save_triangle_plot_false[posAsharpe,posBsharpe][].color = "violetred2"
        val_label_tqA.text[] = string("Earnings A :", round(Int, z1[posAsharpe, posBsharpe]))
        val_label_tqB.text[] = string("B :", round(Int, z2[posAsharpe, posBsharpe]))
    end
    
    if Stackelberg == true
        posAstackelberg = argmax(SharpeA(z1, z2)[:,1])
        posBstackelberg = argmax(SharpeB(z1, z2)[1,:])
        Save_triangle_plot_true[posAstackelberg,posBstackelberg][].color = "navyblue"
        Save_triangle_plot_false[posAstackelberg,posBstackelberg][].color = "navyblue"
        val_label_trA.text[] = string("Earnings A :", round(Int, z1[posAstackelberg, posBstackelberg]))
        val_label_trB.text[] = string("B :", round(Int, z2[posAstackelberg, posBstackelberg]))
    end
    if win_win == true
        posAwinwin = posAminimax
        posBwinwin = posBminimax
        for i in 1:length(z1[:,1])
            for j in 1:length(z1[1,:])
                if i <= posAwinwin && j <= posBwinwin
                    Save_triangle_plot_true[i,j][].alpha = 1.0
                    Save_triangle_plot_false[i,j][].alpha = 1.0
                #elseif i <= posBwinwin && j <= posAwinwin
                #    Save_triangle_plot_true[i,j][].alpha = 1.0
                #    Save_triangle_plot_false[i,j][].alpha = 1.0
                else
                    Save_triangle_plot_true[i,j][].alpha = 0.1
                    Save_triangle_plot_false[i,j][].alpha = 0.1
                end
            end
        end
    end 
    if Shapley == true
        z = argmax(MinimaxA(z1, z2))
        posA = z[1]
        z = argmax(MinimaxB(z1, z2))
        posB = z[2]
        JA = z1[posA, posB]
        JB = z2[posA, posB]
        z1z2 = z1 .+ z2
        A2 = (z1z2 .+ JA .- JB)./2
        B2 = (z1z2 .+ JB .- JA)./2
        z = Shapley_value(z1, z2)
        if win_win == true
            for i in 1:length(z1[:,1])
                for j in 1:length(z1[1,:])
                    if i <= posAwinwin && j <= posBwinwin
                    #elseif i <= posBwinwin && j <= posAwinwin
                    #if i < posAwinwin && j < posBwinwin
                    #elseif i < posBwinwin && j < posAwinwin
                    else
                        z[i,j] = 0
                    end
                end
            end
        end
        posA = argmax(z)
        Save_triangle_plot_true[posA[1],posA[2]][].color = "steelblue"
        Save_triangle_plot_false[posA[1],posA[2]][].color = "steelblue"
        val_label_tsA.text[] = string("Earnings A :", round(Int, A2[posA[1],posA[2]]))
        val_label_tsB.text[] = string("B :", round(Int, B2[posA[1],posA[2]]))
        println(round(Int, B2[posA[1],posA[2]]))
    end
end

function split_heatmap(x, y, z1, z2, Minimax, Sharpe, Stackelberg, Shapley, win_win, color=Symbol("RdYlGn_10"))
    # Normaliser les données
    #z1_normalized = (z1' .- min(minimum(z1), minimum(z2))) ./ (max(maximum(z1), maximum(z2)) - min(minimum(z1), minimum(z2)))
    #z2_normalized = (z2' .- min(minimum(z1), minimum(z2))) ./ (max(maximum(z1), maximum(z2)) - min(minimum(z1), minimum(z2)))

    z1_normalized = ((z1') ./ max(maximum(z1), maximum(z2))) .*0.5 .+ 0.5
    z2_normalized = ((z2') ./ max(maximum(z1), maximum(z2))) .*0.5 .+ 0.5

    z1_normalized = max.(z1_normalized, 0)
    z2_normalized = max.(z2_normalized, 0)
    
    
    # Tracer les triangles inférieurs
    cmap = get(colorschemes[color], range(0, 1, 1001))
    update_triangles(x, y, z1_normalized, cmap)

    # Tracer les triangles supérieurs
    cmap = get(colorschemes[color], range(0, 1, 1001))
    update_triangles(x, y, z2_normalized, cmap, true)

    
    add_game_results(z1, z2, Minimax, Sharpe, Stackelberg, Shapley, win_win)
end

file = matopen("./saveJ_v2_sym.mat")
J_data = read(file, "saveJ")
close(file)

file = matopen("./saveBid_v2_sym.mat")
Bid_data = read(file, "saveBid")
close(file)

J = Observable(J_data)
Bid = Observable(Bid_data)
Cost = Observable(0.0)
Elec = Observable(0.0)
Minimax = Observable(false)
Sharpe = Observable(false)
Stackelberg = Observable(false)
Shapley = Observable(false)
win_win = Observable(false)

x = collect(range(0.20, step=0.02, length=length(J[][:,1,1])))
y = collect(range(0.20, step=0.02, length=length(J[][:,1,2])))
z1,z2 = setPriceMatrices(J, Bid, Cost, Elec)
z = Neutral(z1,z2)

repartition = Makie.Observable(Symbol("Same size"))
strategy = Makie.Observable(Symbol("Earnings A+B"))

screen_resolution = (1900, 1080)
fig = Makie.Figure(size=screen_resolution, fontsize=26)
Label(fig[1, 1], "Comparison of Earnings \nAcross Different \nCharging Price Strategies", fontsize = 24, halign = :left, font = "bold")

ax1 = Makie.Axis3(fig[2:6, 1:4]; aspect=(1,1,1), elevation=π/6,
    perspectiveness=0.5)
ax1.xlabel = "Charging Price A [\$]"
ax1.ylabel = "Charging Price B [\$]"
ax1.zlabeloffset = 110

rectMesh = Rect(Vec(0, 0, 0), Vec(1, 1, 1))
pltobj1 = Makie.meshscatter!(ax1, x, y, 0*z; marker=rectMesh, color=z[:],
markersize=Vec.((x[2] - x[1]), (y[2] - y[1]), z[:]), colormap=:RdYlGn_10,
    shading=NoShading)
Makie.limits!(ax1, minimum(x)-0.01, maximum(x)+0.01, minimum(y)-0.01, maximum(y)+0.01, 
                minimum(z), maximum(z))

ax2 = Makie.Axis(fig[2:6, 5:10])
ax2.xlabel = "Charging Price A [\$]"
ax2.ylabel = "Charging Price B [\$]"

z1_normalized = (z1' .- min(minimum(z1), minimum(z2))) ./ (max(maximum(z1), maximum(z2)) - min(minimum(z1), minimum(z2)))
z2_normalized = (z2' .- min(minimum(z1), minimum(z2))) ./ (max(maximum(z1), maximum(z2)) - min(minimum(z1), minimum(z2)))

cmap = get(colorschemes[Symbol("RdYlGn_10")], range(0, 1, 1001))


Save_triangle_plot_true = plot_triangles(x,y,z2_normalized,ax2,cmap)
Save_triangle_plot_false = plot_triangles(x,y,z1_normalized,ax2,cmap, true)
colorbar = Colorbar(fig[2:6, 11], colormap=cmap, limits=(min(minimum(z1), minimum(z2)), max(maximum(z1), maximum(z2))))

split_heatmap(x,y,z1,z2, Minimax, Sharpe, Stackelberg, Shapley, win_win)

menu = Makie.Menu(fig, options = ["Earnings A+B", "Nominal Non-Collaborative \nNon-Informed for A", 
        "Nominal Non-Collaborative \nNon-Informed for B", "Risk-adjusted Non-Collaborative \nNon-Informed for A", 
        "Risk-adjusted Non-Collaborative \nNon-Informed for B", "Non-Collaborative \nInformed for A", 
        "Non-Collaborative \nInformed for B", "Collaborative Informed"])

fig[7, 1] = vgrid!(Makie.Label(fig, "Strategy", width = nothing, color = :black),
                    menu; tellheight = false, width = 350)

menu_repartition = Makie.Menu(fig, options = ["Same size", "A three times bigger than B"])
fig[1, 2:4] = vgrid!(Makie.Label(fig, "Charging Station size", width = nothing, color = :black),
                    menu_repartition; tellheight = false, width = 350, halign = :left)

sm = GLMakie.Slider(fig, range = 0.05:0.01:0.19, startvalue = 0.1) #0.15
sn = GLMakie.Slider(fig, range = 0.1:0.01:0.30, startvalue = 0.15)  #0.1

val_label_Cost = Makie.Label(fig, width = nothing, color = :black)
val_label_Elec = Makie.Label(fig, width = nothing, color = :black)

Makie.lift(sm.value) do v
    val_label_Cost.text[] = string(v, " \$/kWh")  
end

Makie.lift(sn.value) do v
    val_label_Elec.text[] = string(v, " \$/kWh")  
end

fig[7, 2] = vgrid!(
    Makie.Label(fig, "Operating Cost", width = nothing, color = :black), sm,  val_label_Cost;  
    tellheight = false, width = 200)

fig[7, 3] = vgrid!(
    Makie.Label(fig, "Electricity Cost", width = nothing, color = :black), sn,  val_label_Elec;  
    tellheight = false, width = 200)

gl = GridLayout(fig[7, 6:8], tellwidth = false)

Label(gl[1, 2], "Risk-adjusted Non-Collaborative Non-Informed", halign = :left)
tq = Toggle(gl[1, 1], active = false, buttoncolor = "violetred2")
val_label_tqA = Label(gl[1,4], width = nothing, color = :black, halign = :left)
val_label_tqB = Label(gl[1,17], width = nothing, color = :black, halign = :left)
val_label_tqA.text[] = string("Earnings A : -----")
val_label_tqB.text[] = string("B : -----")

Label(gl[2, 2], "Non-Collaborative Informed", halign = :left)
tr = Toggle(gl[2, 1], active = false, buttoncolor = "navyblue")
val_label_trA = Label(gl[2,4], width = nothing, color = :black, halign = :left)
val_label_trB = Label(gl[2,17], width = nothing, color = :black, halign = :left)
val_label_trA.text[] = string("Earnings A : -----")
val_label_trB.text[] = string("B : -----")

Label(gl[3, 2], "Collaborative Informed", halign = :left)
ts = Toggle(gl[3, 1], active = false, buttoncolor = "steelblue")
val_label_tsA = Label(gl[3,4], width = nothing, color = :black, halign = :left)
val_label_tsB = Label(gl[3,17], width = nothing, color = :black, halign = :left)
val_label_tsA.text[] = string("Earnings A : -----")
val_label_tsB.text[] = string("B : -----")

gm = GridLayout(fig[1, 7:9], tellwidth = false)
Label(gm[2, 2], "Nominal Non-Collaborative Non-Informed", halign = :left)
tp = Toggle(gm[2, 1], active = false, buttoncolor = "purple")
val_label_tpA = Label(gm[2,3], width = nothing, color = :black, halign = :left)
val_label_tpB = Label(gm[2,16], width = nothing, color = :black, halign = :left)
val_label_tpA.text[] = string("Earnings A : -----")
val_label_tpB.text[] = string("B : -----")

Label(gm[1, 2:30], "Activate Restrictive Zone: Access prices under nominal", halign = :left)
tt = Toggle(gm[1, 1], active = false)

connect!(Cost, sm.value)
connect!(Elec, sn.value)

on(menu.selection) do s
    strategy[] = Symbol(s)
end

on(menu_repartition.selection) do s
    repartition[] = Symbol(s)
end

on(tp.active) do is_active
    Minimax[] = is_active 
    if !is_active && tt.active[]
        tt.active[] = false  # Désactive tt si tp est désactivé
    end
    if !is_active
        val_label_tpA.text[] = string("Earnings A : -----")
        val_label_tpB.text[] = string("B : -----") 
    end
end

on(tq.active) do is_active
    Sharpe[] = is_active 
    if !is_active
        val_label_tqA.text[] = string("Earnings A : -----")
        val_label_tqB.text[] = string("B : -----") 
    end
end

on(tr.active) do is_active
    Stackelberg[] = is_active 
    if !is_active
        val_label_trA.text[] = string("Earnings A : -----")
        val_label_trB.text[] = string("B : -----") 
    end
end

on(ts.active) do is_active
    Shapley[] = is_active 
    if !is_active
        val_label_tsA.text[] = string("Earnings A : -----")
        val_label_tsB.text[] = string("B : -----") 
    end
end

on(tt.active) do is_active
    win_win[] = is_active 
    tp.active[] = is_active
    if !is_active
        tp.active[] = false  
    end
end

Makie.lift(strategy, Cost, Elec, Minimax, Sharpe, Stackelberg, Shapley, win_win, repartition) do s,
    Cost, Elec, Minimax, Sharpe, Stackelberg, Shapley, win_win, r
    s = String(s)
    r = String(r)
    if r == "Same size"
        file = matopen("./saveJ_v2_sym.mat")
        J_data = read(file, "saveJ")
        close(file)
        file = matopen("./saveBid_v2_sym.mat")
        Bid_data = read(file, "saveBid")
        close(file) 

        J[] = J_data
        Bid[] =Bid_data
    elseif r == "A three times bigger than B"
        file = matopen("./saveJ_v2_assym.mat")
        J_data = read(file, "saveJ")
        close(file)

        file = matopen("./saveBid_v2_assym.mat")
        Bid_data = read(file, "saveBid")
        close(file) 

        J[] = J_data
        Bid[] =Bid_data
    end
    z1 = zeros(size(J[][:,:,1]))
    z2 = zeros(size(J[][:,:,2]))
    z1, z2 = setPriceMatrices(J, Bid, Cost, Elec)
    z = zeros(size(z1))
    if s == "Nominal Non-Collaborative \nNon-Informed for A"
        z = MinimaxA(z1, z2)
        ax1.zlabel = "Earnings A+B [\$]"
    elseif s == "Nominal Non-Collaborative \nNon-Informed for B"
        z = MinimaxB(z1, z2)
        ax1.zlabel = "Earnings A+B [\$]"
    elseif s == "Risk-adjusted Non-Collaborative \nNon-Informed for A"
        z = SharpeA(z1, z2)
        ax1.zlabel = "Sharpe Value [-]"
    elseif s == "Risk-adjusted Non-Collaborative \nNon-Informed for B"
        z = SharpeB(z1, z2)
        ax1.zlabel = "Sharpe Value [-]"
    elseif s == "Non-Collaborative \nInformed for A"
        z = StackelbergA(z1, z2)
        ax1.zlabel = "Stackelberg Value [-]"
    elseif s == "Non-Collaborative \nInformed for B"
        z = StackelbergB(z1, z2)
        ax1.zlabel = "Stackelberg Value [-]"
    elseif s == "Collaborative Informed"
        z = Shapley_value(z1, z2)
        ax1.zlabel = "Shapley Value [-]"
    elseif s == "Earnings A+B"
        z = Neutral(z1, z2)
        ax1.zlabel = "Earnings [\$]"
    end
    pltobj1.color = z[:]
    if maximum(z)<=0
        pltobj1.colorrange = (minimum(z), maximum(z))
    else
        pltobj1.colorrange = (-maximum(z), maximum(z))
    end
    pltobj1.markersize = Vec.((x[2] - x[1]), (y[2] - y[1]), z[:])
    Makie.limits!(ax1, minimum(x)-0.01, maximum(x)+0.01, minimum(y)-0.01, maximum(y)+0.01, 
                minimum(z), maximum(z))
    split_heatmap(x, y, z1, z2, Minimax, Sharpe, Stackelberg, Shapley, win_win)
    z = Neutral(z1, z2)
    if maximum(z)<=0
        colorbar.limits = (minimum(z), maximum(z))
    else
        colorbar.limits = (-maximum(z), maximum(z))
    end
    
end

display(fig)
