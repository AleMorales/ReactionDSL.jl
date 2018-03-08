export ReactionODE

# Convert to DiffEqSystem where a first order derivative is created for each
# species. Statements that do not result in Reactions are simply ignored and
# passed along
function ReactionODE(eqns, iv)
    # Separate reactions from the rest of equations
    reactions = Reaction[]
    neqns = Operation[]
    for eq in eqns
        eq isa Reaction ? push!(reactions, eq) : push!(neqns, eq)
    end
    # Create derivatives for species
    spnames, derivs = create_derivatives!(reactions, neqns, iv)
    # Expand derivatives with rate equations
    for reaction in reactions
        process_reaction!(reaction, spnames, derivs)
    end
    # Merge new derivatives with the rest of the system
    append!(neqns, derivs)
    # Convert all Species to DependentVariable
    for eqn in neqns
        convert_species!(eqn)
    end
    DiffEqSystem(neqns, [iv])
end


# Create derivatives for species
function create_derivatives!(reactions, neqns, iv)
    # Get all species and their names
    spnames, species = get_species(reactions)
    # Create derivatives for remaining species
    D = Differential(iv, 1)
    derivs = [D*DependentVariable(sp) == 1 for sp in species]
    return spnames, derivs
end


# Retrieve species used as reactants in reaction formulae
function get_species(reactions::Vector{Reaction})
    spnames = Symbol[]
    species = Species[]
    for reaction in reactions
        form = reaction.formula
        for react in vcat(form.substrates, form.products)
            if !in(react.species.name, spnames)
                push!(spnames, react.species.name)
                push!(species, react.species)
            end
        end
    end
    return spnames, species
end


# Process a reaction and update derivatives in the model
function process_reaction!(reaction, spnames, derivs)
    form = reaction.formula
    arrow = form.arrow
    for react in form.substrates
        expand_derivatives!(react, reaction, spnames, derivs, arrow, -1)
    end
    for react in form.products
        expand_derivatives!(react, reaction, spnames, derivs, arrow, 1)
    end
end


# Process a reaction and transfer expressions to derivatives
# dir indicates whether to substractor add the rate equation
function expand_derivatives!(react, reaction, spnames, derivs, arrow, dir)
    pos = findfirst(spnames, react.species.name)
    derivs[pos].args[2] *= dir*create_rate(react, reaction, arrow)
end


# Create a rate expression depending on the type of arrow and differences in compartments
function create_rate(react, reaction, arrow)
    sc = react.species.compartment
    fc = reaction.formula.compartment
    compcorr = sc == fc ? false : true
    if open[arrow]
        if reversible[arrow]
            out = compcorr ? react.stoich*reaction.rate*fc/sc : react.stoich*reaction.rate
        else
            error("only reversible arrows implemented")
        end
    else
        error("closed arrows not implemented") # Here we would use mass action kinetics
    end
    return out
end

# Convert all Species into DependentVariable
function convert_species!(eqn::Operation)
    for i in 1:length(eqn.args)
        arg = eqn.args[i]
        if arg isa Species
            eqn.args[i] = DependentVariable(arg)
        elseif arg isa Operation
            convert_species!(eqn.args[i])
        end
    end
end
