-- Load utility scripts (e.g. from from ugcore/scripts)
ug_load_script("ug_util.lua")
ug_load_script("util/refinement_util.lua")

-- Parse parameters and print help

ARGS ={
 gridName	= util.GetParam("--grid", "Room.ugx","filename of underlying grid"),
 requiredSubsets= {"Inner", "Wall", "Heater", "North", "South", "East", "West"},
 numRefs		= util.GetParamNumber("--numRefs", 2, "number of refinements"),

 steadyState	= util.GetParamBool("--steadyState",true,"If specified, the steady state of the problem is computed. Else a time-dependent problem is computed."),
}


-- initialize ug with the world dimension 3 and an algebra system with scalar coefficients
InitUG(2, AlgebraType("CPU", 1));

-- Load a domain without initial refinements.
dom = util.CreateDomain(ARGS.gridName, 0, ARGS.requiredSubsets)

-- Refine the domain (redistribution is handled internally for parallel runs)
print("refining...")
util.refinement.CreateRegularHierarchy(dom, ARGS.numRefs, true)


-----------------------------------------
-- A) Modellparameter
--
-----------------------------------------
alpha={ -- achten auf Einheit
--todo
    ["Beton"] = 0.994 , ["Luft"] = 20.0, ["Holz"] = 0.12
}

-----------------------------------------
-- B) Ansatzraum
-----------------------------------------
approxSpaceDesc = { fct = "temp", type = "Lagrange", order = 1 }

approxSpace = ApproximationSpace(dom)
approxSpace:add_fct(approxSpaceDesc.fct, approxSpaceDesc.type, approxSpaceDesc.order)
approxSpace:init_levels()
approxSpace:init_top_surface()

print("Approximation space:")
approxSpace:print_statistic()

-----------------------------------------
-- B) Ansatzraum
-----------------------------------------
approxSpaceDesc = { fct = "temp", type = "Lagrange", order = 1 }

approxSpace = ApproximationSpace(dom)
approxSpace:add_fct(approxSpaceDesc.fct, approxSpaceDesc.type, approxSpaceDesc.order)
approxSpace:init_levels()
approxSpace:init_top_surface()

print("Approximation space:")
approxSpace:print_statistic()


-----------------------------------------
-- C) Elementdiskretisierung (FV)
-----------------------------------------
elemDisc={}
-----------------------------------------
-- C).1 Diffusionstensor
-----------------------------------------
for index, vol in ipairs(requiredSubsets) do
elemDisc[vol]  = ConvectionDiffusion("temp", vol, "fv1")
--todo
elemDisc[vol]: set_mass_scale(1.0)
elemDisc[vol]: set_diffusion(1.0)

end

-----------------------------------------
-- D) Anfangswerte
-----------------------------------------
function InitialValue(x,y,t,si)
  Inner = 4.0
  Wall = 0.0
end

-----------------------------------------
-- E) Randwerte
-----------------------------------------
dirichletBnd = DirichletBoundary()
dirichletBnd:add(30.0, "temp", "Heater")
dirichletBnd:add(4.0, "temp", "North")
dirichletBnd:add(4.0, "temp", "West")

-----------------------------------------
-- F) Diskretisierung auf ganzem Gebiet
-----------------------------------------
domainDisc = DomainDiscretization(approxSpace)
for index, vol in ipairs(requiredSubsets) do
domainDisc:add(elemDisc[vol])
end
domainDisc:add(dirichletBnd)

-----------------------------------------
-- G) Aufsetzen des Lösers
--    (using 'util/solver_util.lua')
-----------------------------------------
solver =LU() 
u = GridFunction(approxSpace)

-----------------------------------------
-- H) Lösen
-----------------------------------------

if ARGS.steadyState then
	local A = AssembledLinearOperator(domainDisc)
	local b = GridFunction(approxSpace)
	domainDisc:adjust_solution(u)
	domainDisc:assemble_linear(A, b)

	solver:init(A, u)
	solver:apply(u, b)

	local solFileName = "sol_temp"
	print("writing solution to '" .. solFileName .. "'...")
	WriteGridFunctionToVTK(u, solFileName)
	SaveVectorForConnectionViewer(u, solFileName .. ".vec")
else

    --* Set initial value  
  Interpolate("InitialValue", u, "temp")
    --Execute time stepping loop w/ fixed time-step 
    print("\nsolving...")
	local startTime = 0.0 --todo
    local endTime=1440 --todo 24hours in minutes
    local dt= 0.6 --todo

	util.SolveLinearTimeProblem(u, domainDisc, solver, VTKOutput(), "sol_temp",
								"ImplEuler", 1, startTime, endTime, dt);

end

print("done")

