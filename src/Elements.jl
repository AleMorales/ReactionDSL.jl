export Species, @Species, Environment, ∅

# Chemical species are lowered to DependentVariable
Species(name, value = nothing, value_type = typeof(value)) =
             Variable(name, :Species, value, value_type, nothing)
dvar_species(s) = DependentVariable(s.name, s.value, s.value_type, s.diff)

# Allows to create several Species with or without default values as @Species s1, s2, s3 = 2.5
macro Species(x...)
        esc(SciCompDSL._parse_vars("Species", :Species, x))
end


# To represent sources and sinks in a reaction (creation/degration to/from nothing)
struct Environment end
const ∅ = Environment()

# The basic element of a reaction. Stoichiometry does not have to be a literal
struct Reactant
    species::Variable
    stoich::Variable
end

# The formula of a reaction. The type of arrow is stored as type parameter
struct Formula
    arrow::Symbol
    substrates::Vector{Reactant}
    products::Vector{Reactant}
end
Base.isequal(r1::Formula, r2::Formula) = isequal(r1.substrates, r2.substrates) &&
    isequal(r1.products, r2.products) && isequal(r1.arrow, r2.arrow)

# Reaction holds formula and its rate equation in SciComp IR
struct Reaction
    formula::Formula
    rate::Expression
end
Base.isequal(r1::Reaction, r2::Reaction) = isequal(r1.formula, r2.formula) &&
    isequal(r1.rate, r2.rate)
# Reactions created by pattern (Formula) ~ Operation
Base.:~(f::Formula, o::Expression) = Reaction(f,o)
