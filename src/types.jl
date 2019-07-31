##### types.jl #####
# Created by Alexander Dowling

##############
##### Functions and Types for setup

mutable struct IrreducibleDegenerateSet
	c::Int64
	nElements::Float64
	lambda::Array{Float64,1}
	y::Array{Float64,1}
end

mutable struct DegenSettings
	includeBounds::Bool
	includeWeaklyActive::Bool
	removeFixedVar::Bool
	onlyCandidateSearch::Bool
	epsiActive::Float64
	epsiLambda::Float64
	lambdaM::Float64
	mySolver
	printSmallestSVD::Float64
end

function DegenSettings(mySolver)
	return DegenSettings(false, true, true, true, 1E-6, 1E-6, 1.0E5, mySolver, 0.0)
end

mutable struct DegenData

# Point to analyze
	x::Array{Float64,1}

# Jacobian expressed as a sparse matrix
	J_active

#	NLP evaluator
	d

# Full Jacobian (pull entries for bounds, if appropriate per settings)
	J::Array{Float64,1}
	iR::Array{Int64,1}
	jC::Array{Int64,1}

# Constraint information
	gLB::Array{Float64,1}
	gUB::Array{Float64,1}
	g::Array{Float64,1}

# Maps from active Jacobian to full Jacobian
	gMap::Array{Int64,1}
	bMap::Array{Int64,1}
	vMap::Array{Int64,1}

# Size info
	nLambda::Int64
	nVarActive::Int64
	nVar::Int64
	nConstr::Int64

end

function DegenData()
	return DegenData(Float64[],
						nothing,
						nothing,
						Float64[],
						Int64[],
						Int64[],
						Float64[],
						Float64[],
						Float64[],
						Int64[],
						Int64[],
						Int64[],
						0,
						0,
						0,
						0)
end

function DegenData(m::Model, f=stdout, status::Symbol=:Unknown)

	dd = DegenData()

# Select point to analyze
	dd.x = get_x(m, status)

	dd.nVar = length(dd.x)

# Evaluate constraint bounds
	(dd.gLB, dd.gUB) = JuMP.constraintbounds(m)

	if(!m.internalModelLoaded)
		JuMP.build(m)
	end
	
	if(typeof(m.internalModel) == MathProgBase.SolverInterface.NonlinearToLPQPBridge)
		processNonlinearModel!(m, dd, f)
	
	# TO DO: Update processLinearModel! to work with quadratic constraints
	elseif(typeof(m.internalModel) <: MathProgBase.SolverInterface.AbstractLinearQuadraticModel)
		# processLinearModel!(m, dd, f)
		processNonlinearModel!(m, dd, f)
		
	elseif(typeof(m.internalModel) <: MathProgBase.SolverInterface.AbstractNonlinearModel)
		processNonlinearModel!(m, dd, f)
	
	else
		println("Warning: Conic models are not yet supported!")
	end

	dd.nConstr = MathProgBase.numconstr(m)

	return dd

end

function get_x(m::Model, status::Symbol = :Unknown)
	if(m.internalModel === nothing || status == :InfeasibleOrUnbounded || status == :Unsolved || status == :Unknown)
		x = m.colVal
	else
	# Note: getsolution will throw an error if model is infeasible or unbounded with Gurobi
		x = MathProgBase.getsolution(m.internalModel)
	end
	
	return x
end

function processNonlinearModel!(m::Model, dd::DegenData, f=stdout)

# Initialize NLP evaluatorSetup

	print(f,"Setting up NLP evaluator... ")
	t0 = time()
	dd.d = JuMP.NLPEvaluator(m)
	tm = time() - t0
	println(f, string(tm, " seconds"))

	print(f,"Initializing... ")
	t0 = time()
	MathProgBase.initialize(dd.d, [:ExprGraph, :Jac])
	tm = time() - t0
	println(f, string(tm," seconds"))

# Evaluate constraints
	n = MathProgBase.numconstr(m)

	dd.g = zeros(n)
	MathProgBase.eval_g(dd.d,dd.g,dd.x)

# Query Jacobian structure

	# Rows: Constraints
	# Columns: Variables

	print(f,"Determining Jacobian structure... ")
	t0 = time()
	(dd.iR, dd.jC) = MathProgBase.jac_structure(dd.d)
	tm = time() - t0
	println(f,string(tm," seconds"))

	dd.J = zeros(length(dd.iR))

# Query Jacobian at point
	print(f,"Evaluating Jacobian... ")
	t0 = time()
	MathProgBase.eval_jac_g(dd.d,dd.J,dd.x)
	tm = time() - t0
	println(f,string(tm, " seconds"))

	return nothing

end

# TODO: Update to work with quadratic models.
function processLinearModel!(m::Model, dd::DegenData, f=stdout)

	# Jacobian as a SparseMatrixCSC
	if(!m.internalModelLoaded)
		JuMP.build(m)
	end
	A = MathProgBase.getconstrmatrix(internalmodel(m))

	# Translate to row, column, element form
	(dd.iR, dd.jC, dd.J) = findnz(A)

	# Get RHS of constraints
	# dd.g = MathProgBase.getconstrsolution(internalmodel(m))
	dd.g = A*dd.x

	return nothing

end
