export Species, @Species, Environment, ∅

# Default compartment. A compartment is just a SciComp variable (any subtype!)
const _global_compartment_ = Constant(:_global_compartment_,1)

# Chemical species like Variable but they exist inside a compartment and it
# is lowered to a DependentVariable
struct Species <: Expression
    name::Symbol
    compartment::Variable
    value
    value_type::DataType
end
# Facilitate construction of Species
Species(name, value = nothing, value_type = typeof(value)) =
             Species(name, _global_compartment_, value, value_type)
macro Species(x...)
     esc(SciCompDSL._parse_vars("Species", :Species, x))
end
Base.in(s::Species, c::Variable) = Species(s.name, c, s.value, s.value_type)
# Lowering to DependentVariable
SciCompDSL.DependentVariable(s::Species) = SciCompDSL.DependentVariable(s.name, s.value, s.value_type)

# Need to behave like a variable for printing and arithmetic (adapted from SciCompDSL.jl)
Base.Expr(x::Species) = :($(x.name))
function Base.show(io::IO, A::Species)
        str = "Species($(A.name))"
        if A.value != nothing
            str *= ", value = " * string(A.value)
        end
        print(io,str)
end
function Base.:(==)(x::Species,y::Species)
    x.name == y.name && x.compartment == y.compartment && x.value == y.value &&
    x.value_type == y.value_type
end

# To represent sources and sinks in a reaction (creation/degration to/from nothing)
struct Environment end
const ∅ = Environment()


# The basic element of a reaction. Stoichiometry does not have to be a literal
struct Reactant
    species::Species
    stoich::Variable
end


# The formula of a reaction. The type of arrow is stored as type parameter
struct Formula
    arrow::Symbol
    compartment::Variable
    substrates::Vector{Reactant}
    products::Vector{Reactant}
end
Formula(arrow, substrates, products) = Formula(arrow, _global_compartment_, substrates, products)
Base.:(==)(f1::Formula, f2::Formula) = f1.substrates == f2.substrates &&
    f1.products == f2.products && f1.arrow == f2.arrow && f1.compartment == f2.compartment
# Change default compartment of a Formula
Base.in(f::Formula, c::Variable) = Formula(f.arrow, c, f.substrates, f.products)

# Reaction holds formula and its rate equation in SciComp IR
struct Reaction
    formula::Formula
    rate::Expression
end
Base.:(==)(r1::Reaction, r2::Reaction) = r1.formula == r2.formula &&
    r1.rate == r2.rate
# Reactions created by pattern (Formula) ~ Operation
Base.:~(f::Formula, o::Expression) = Reaction(f,o)
