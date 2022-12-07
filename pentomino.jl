"""
All solutions of the pentomino farm problem
"""

using JuMP
using Gurobi
using DelimitedFiles

function rotate(p::Matrix, d::Int64)  # 4 rotations of the pentomino and 4 rotations of the mirrored pentomino
    if (d == 1)
        return p
    elseif (d == 2)
        return rotl90(p)
    elseif (d == 3)
        return rot180(p)
    elseif (d == 4)
        return rotr90(p)
    else
        return reverse(rotate(p, d - 4), dims=2)
    end  
end

function printSol(n::Int64, z, width, height, pieces)  # print n-th solution
    output = fill(0, width, height)

    solution = findall(t->t == 1.0, round.(value.(z; result = n)))
    for i in 1 : length(solution)
        pieceSet = rotate(pieces[solution[i][3]], solution[i][4])

        for j in 1 : size(pieceSet, 1)
            for l in 1 : size(pieceSet, 2)
                if (pieceSet[j, l])
                    output[j + solution[i][2] - 1, l + solution[i][1] - 1] = solution[i][3]
                end
            end
        end
    end

    display(output)
end

function order(n::Int64, z, width, height, pieces)  # determine order of pieces
    output = fill(0, width, height)

    solution = findall(t->t == 1.0, round.(value.(z; result = n)))  # create output matrix
    for i in 1 : length(solution)
        pieceSet = rotate(pieces[solution[i][3]], solution[i][4])

        for j in 1 : size(pieceSet, 1)
            for l in 1 : size(pieceSet, 2)
                if (pieceSet[j, l])
                    output[j + solution[i][2] - 1, l + solution[i][1] - 1] = solution[i][3]
                end
            end
        end
    end

    neighList = Vector{Set{Int}}()  # find two neighbours of every pentomino
    for i in 1 : 12
        push!(neighList, Set{Int}())
    end

    for x in 1 : width
        for y in 1 : height
            for d in ((0,1),(0,-1),(1,0),(-1,0))
                neigh = (x, y) .+ d

                if (neigh[1] > 0 && neigh[2] > 0 && neigh[1] <= width && neigh[2] <= height)
                    if output[x, y] > 0 && output[neigh[1], neigh[2]] > 0 && output[x, y] != output[neigh[1], neigh[2]]
                        push!(neighList[output[x, y]], output[neigh[1], neigh[2]])
                    end
                end
            end
        end
    end

    ord = Vector{Int}()  # bring neighbour list in order
    push!(ord, 1)
    discovered = Set{Int}()
    push!(discovered, 1)
    while !isempty(setdiff(neighList[ord[length(ord)]], discovered))
        n = maximum(setdiff(neighList[ord[length(ord)]], discovered))
        push!(discovered, n)
        push!(ord, n)
    end

    return ord
end

function solvePentomino(area = Pair(20, 20))
    """
    Setup
    """

    # bounding box to search in
    width = area[1]
    height = area[2]

    pieces = Vector{Matrix{Bool}}()

    # Pentominoes as binary arries column by column
    push!(pieces, [false; true; false;; true; true; true;; true; false; false])  # F
    push!(pieces, [true;; true;; true;; true;; true])  # I
    push!(pieces, [true; true; true; true;; false; false; false; true])  # L
    push!(pieces, [false; true; true; true;; true; true; false; false])  # N
    push!(pieces, [true; true; true;; true; true; false])  # P
    push!(pieces, [true; false; false;; true; true; true;; true; false; false])  # T
    push!(pieces, [true; true;; false; true;; true; true])  # U
    push!(pieces, [true; true; true;; false; false; true;; false; false; true])  # V
    push!(pieces, [true; true; false;; false; true; true;; false; false; true])  # W
    push!(pieces, [false; true; false;; true; true; true;; false; true; false])  # X
    push!(pieces, [false; true; false; false;; true; true; true; true])  # Y
    push!(pieces, [true; false; false;; true; true; true;; false; false; true])  # Z
    
    active = Dict{Pair{Int64, Int64}, Set{VariableRef}}()  #  Set of all z[w, h, p, d] that place a tile on (w, h)
    for i in 1 : width
        for j in 1 : height
            active[Pair(i, j)] = Set{VariableRef}()
        end
    end

    # Gurobi Solver
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "PoolSearchMode", 2)
    set_optimizer_attribute(model, "PoolSolutions", typemax(Int32))
    set_optimizer_attribute(model, "PoolGapAbs", 0.01)

    @variable(model, x[1:width, 1:height], Bin)  # if tile is interior
    @variable(model, y[1:width, 1:height], Bin)  # if tile is boundary
    @variable(model, z[1:width, 1:height, 1:length(pieces), 1:8], Bin)  # z[w, h, p, d] the left top corner of polyomino p in orientation d is placed on (w, h) with (1, 1) being top left conrer

    @objective(model, Max, sum(x))  # Maximize area enclosed

    """
    Constraints
    """

    # every polyomino is placed only once in one direction
    for i in 1 : length(pieces)
        @constraint(model, sum(z[:, :, i, :]) == 1)
    end

    # verify pieces dont go out of bounds and fill active set
    for w in 1 : width
        for h in 1 : height
            for p in 1 : length(pieces)
                for d in 1 : 8
                    piece = rotate(pieces[p], d)

                    if (w + size(piece, 2) - 1 > width)  # too far right
                        @constraint(model, z[w, h, p, d] == 0)
                    elseif (h + size(piece, 1) - 1 > height)  # too far down
                        @constraint(model, z[w, h, p, d] == 0)
                    else
                        for i in 1 : size(piece, 2)
                            for j in 1 : size(piece, 1)
                                if (piece[j, i])
                                    push!(active[Pair(w + i - 1, h + j - 1)], z[w, h, p, d])
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    for w in 1 : width
        for h in 1 : height
            @constraint(model, sum(active[Pair(w, h)]) <= 1)  # no polyominoes overlap
            @constraint(model, sum(active[Pair(w, h)]) == y[w, h])  # border definition
        end
    end

    # if neighbour of interior is not boundary, than it is also interior: (not x[l]) or y[l'] or x[l']
    for w in 1 : width
        for h in 1 : height
            if (w != 1)
                @constraint(model, (1 - x[w, h]) + y[w - 1, h] + x[w - 1, h] >= 1)  # left
            end
            if (h != 1)
                @constraint(model, (1 - x[w, h]) + y[w, h - 1] + x[w, h - 1] >= 1)  # above
            end
            if (w != width)
                @constraint(model, (1 - x[w, h]) + y[w + 1, h] + x[w + 1, h] >= 1)  # right
            end
            if (h != height)
                @constraint(model, (1 - x[w, h]) + y[w, h + 1] + x[w, h + 1] >= 1)  # below
            end

            if (w != 1 && h != 1)
                @constraint(model, (1 - x[w, h]) + y[w - 1, h - 1] + x[w - 1, h - 1] >= 1)  # left above
            end
            if (h != 1 && w != width)
                @constraint(model, (1 - x[w, h]) + y[w + 1, h - 1] + x[w + 1, h - 1] >= 1)  # right above
            end
            if (w != width && h != height)
                @constraint(model, (1 - x[w, h]) + y[w + 1, h + 1] + x[w + 1, h + 1] >= 1)  # right below
            end
            if (h != height && w != 1)
                @constraint(model, (1 - x[w, h]) + y[w - 1, h + 1] + x[w - 1, h + 1] >= 1)  # left below
            end
        end
    end
    
    # edges of box are neither boundary nor interior
    for w in 1 : width
        for h in 1 : height
            if (w == 1 || w == width || h == 1 || h == height)
                @constraint(model, x[w, h] == 0)
                @constraint(model, y[w, h] == 0)
            end
        end
    end

    # no tile is boundary and interior
    for w in 1 : width
        for h in 1 : height
            @constraint(model, x[w, h] + y[w, h] <= 1)
        end
    end

    # extra constraints to achieve uniqueness

    # one tile has to be in the second row and one in the second column, prevents duplicate solutions of translation
    @constraint(model, sum(y[2, :]) >= 1)
    @constraint(model, sum(y[:, 2]) >= 1)

    # exclude duplicate orientations caused by symmetries in pentominoes
    @constraint(model, sum(z[:, :, 2, 3:8]) == 0)  # I-Pentomino: 2 axis
    @constraint(model, sum(z[:, :, 6, 5:8]) == 0)  # T-Pentomino: 1 axis
    @constraint(model, sum(z[:, :, 7, 5:8]) == 0)  # U-Pentomino: 1 axis
    @constraint(model, sum(z[:, :, 8, 5:8]) == 0)  # V-Pentomino: 1 axis
    @constraint(model, sum(z[:, :, 9, 5:8]) == 0)  # W-Pentomino: 1 axis
    @constraint(model, sum(z[:, :, 10, 2:8]) == 0)  # X-Pentomino: 4 axis
    @constraint(model, sum(z[:, :, 12, 3:4]) == 0)  # Z-Pentomino: Rotational symmetry
    @constraint(model, sum(z[:, :, 12, 7:8]) == 0)

    # force unsymmetric F-pentomino into specific orientation to prevent rotational-symmetric solutions
    @constraint(model, sum(z[:, :, 1, 1]) == 1)

    """
    Calculate
    """

    optimize!(model)

    # print first solution
    println()
    printSol(1, z, width, height, pieces)
    println()

    println("Number of solutions: " * string(result_count(model)))
    println()

    println("Saving solutions to disk.")

    uniq = Set{Vector}()
    push!(uniq, order(1, z, width, height, pieces))

    output = hcat(Int64(sum(round.(value.(x)))), order(1, z, width, height, pieces)')
    for i in 2 : result_count(model)
        ord = order(i, z, width, height, pieces)

        if !(ord in uniq)
            output = [output; hcat(Int64(sum(round.(value.(x; result = i)))), ord')]
            push!(uniq, ord)
        end
    end
    writedlm("solutions.csv",  output, ',')
    println("Saving completed.")
end

# (c) Muessig
