function verify_QE_result(results, program)
    if program.exec[1]=="pw.x"
        return occursin("JOB DONE.", results[end-1])
    else
        return true
    end
end


function QE_run(
    scripts::String, 
    program::Cmd;
    workspace = "/tmp/",
    f_in = "tmp.in",
    f_out = "tmp.out"
    )::Tuple{Vector{String},Bool}

    input_file_loc = workspace*f_in
    output_file_loc = workspace*f_out
    @inline split_output(outp) = String[ st for st ∈ split(outp,'\n',keepempty=false) ]

    # write to temporary file 
    open(workspace*f_in,"w") do io
        write(io, scripts)
    end
    
    # run a pipeline :
    # cat tmp.in | pw.x | tmp.out
    try
        program_run = run(pipeline(`cat $input_file_loc`, program, output_file_loc))
    catch
        println("Program $(join(program.exec,' ')) lauching error.")
        println("Please check the scripts.")
        println("Exit with empty results.")
        return String[], false
    end

    # read from temporary file 
    res = nothing
    open(output_file_loc) do io
        res = read(io)
    end
    results = String(Char.(res)) |> split_output
    success = verify_QE_result(results, program)

    # process and return
    return results, success
end


# pw.x
function pw(
    s::String;
    workspace = "",
    fin="pw.x.tmp.in",
    fout="pw.x.tmp.out"
    )
    QE_run( s, `pw.x`, workspace=workspace, f_in=fin, f_out=fout )
end


## -------------------------------------------------------------

# process pw.x results

function extract(res, reg1,reg2)
    pos = findfirst(x->occursin(reg1,x), res)
    m1 = (pos != nothing ? match(reg1,res[pos]).match : "")
    m2 = (m1!="" ? match(reg2,m1).match : "")
    m2
end

energy_line_rstr = r"\s*total\s+energy\s*\=\s*-?\d+\.\d+\s+Ry"

num_f_rstr = r"-?\d+\.\d+"

extract_energy(result_lines) = parse(Float64,extract(result_lines, energy_line_rstr, num_f_rstr))


## -----------------------------------------------------------------

global const PWX_INPUT_FILE = """&control
title = 'TITLE'
pseudo_dir = 'PSEUDO_FD'
calculation = 'scf'
tstress=.false.
/
&system
ibra v = 2
celldm(1) = CELL_DM
nat=  2 
ntyp= NUM_ATOM_TYPES
ecutwfc= 12.0
nbnd= 8
input_dft='hse' 
nqx1=1
nqx2=1
nqx3=1
exxdiv_treatment='gygi-baldereschi'
x_gamma_extrapolation= .true.
/
&electrons
/
ATOMIC_SPECIES
SPECIES_LINES
ATOMIC_POSITIONS (ATOM_POSITION_UNIT)
POSITION_LINES
K_POINTS KP_LIST
KP_LINES"""


function input_str(
    pseudo_fd::String, 
    title::String
    )

end

## -----------------------------------------------------------------

en = []

for d=0.24:0.002:0.26
    res, suc = pw( input_str(d) )
    push!(en, extract_energy(res))
    println(suc)
end

##

using Plots

##

plot(collect(0.24:0.002:0.26), en)

##


## ASE version

ENV["ASE_ESPRESSO_COMMAND"] = "/usr/bin/pw.x -in NaCl.pwi > NaCl.pwo"
ENV["ESPRESSO_PSEUDO"] = "/home/dabajabaza/abinitio/QE/SSSP/efficiency/"

##

using PyCall

ase = pyimport("ase")
build = pyimport("ase.build")
espresso = pyimport("ase.calculators.espresso")
constraints = pyimport("ase.constraints")
optimize = pyimport("ase.optimize")

##

bulk = build.bulk
Espresso = espresso.Espresso
UnitCellFilter = constraints.UnitCellFilter
LBFGS = optimize.LBFGS

##

pseudopotentials = Dict("Na" => "na_pbesol_v1.5.uspp.F.UPF",
                        "Cl" => "cl_pbesol_v1.4.uspp.F.UPF")

##

rocksalt = bulk("NaCl", crystalstructure="rocksalt", a=6.0)

calc = Espresso(pseudopotentials=pseudopotentials,
                tstress=true, tprnfor=true, kpts=(3, 3, 3))

calc.set_label("tmp/NaCl")

rocksalt.calc = calc

##

ucf = UnitCellFilter(rocksalt)

##

opt = LBFGS(ucf)

##

opt.run(fmax=0.005)

##



using PyCall

ase = pyimport("ase")
build = pyimport("ase.build")
espresso = pyimport("ase.calculators.espresso")
constraints = pyimport("ase.constraints")
optimize = pyimport("ase.optimize")
