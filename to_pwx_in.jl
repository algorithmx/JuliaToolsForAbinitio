using Printf

## -------------------------------------------------
## DATA
## -------------------------------------------------

global const _BOHR_ = 0.52917720859

global const AtomSymb = split("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt", " ", keepempty=false)

global const AtomMass = Float64[
1.0079,4.0026,6.941,9.0122,10.811,12.0107,14.0067,15.9994,18.9984,20.1797,22.9897,24.305,26.9815,28.0855,30.9738,32.065,35.453,
39.948,39.0983,40.078,44.9559,47.867,50.9415,51.9961,54.938,55.845,58.9332,58.6934,63.546,65.39,69.723,72.64,74.9216,78.96,79.904,
83.8,85.4678,87.62,88.9059,91.224,92.9064,95.94,98,101.07,102.906,106.42,107.868,112.411,114.818,118.71,121.76,127.6,126.904,131.293,
132.905,137.327,138.905,140.116,140.908,144.24,145,150.36,151.964,157.25,158.925,162.5,164.93,167.259,168.934,173.04,174.967,178.49,
180.948,183.84,186.207,190.23,192.217,195.078,196.966,200.59,204.383,207.2,208.98,209,210,222,223,226,227,232.038,231.036,238.029,
237,244,243,247,247,251,252,257,258,259,262,261,262,266,264,277,268
]

global const Atoms = Dict(AtomSymb[i]=>AtomMass[i] for i=1:length(AtomMass))


# International Tables
# 1-2 Triclinic, 3-15 Monoclinic, 16-74 Orthorhombic, 75-142 Tetragonal, 143-167 Trigonal, 168-194 Hexagonal, 195-230 Cubic

global const Int_Tables = split("P1 P-1 P2 P2(1) C2 Pm Pc Cm Cc P2/m P2(1)/m C2/m P2/c P2(1)/c C2/c P222 P222(1) P2(1)2(1)2 P2(1)2(1)2(1) C222(1) C222 F222 I222 I2(1)2(1)2(1) Pmm2 Pmc2(1) Pcc2 Pma2 Pca2(1) Pnc2 Pmn2(1) Pba2 Pna2(1) Pnn2 Cmm2 Cmc2(1) Ccc2 Amm2 Abm2 Ama2 Aba2 Fmm2 Fdd2 Imm2 Iba2 Ima2 Pmmm Pnnn Pccm Pban Pmma Pnna Pmna Pcca Pbam Pccn Pbcm Pnnm Pmmn Pbcn Pbca Pnma Cmcm Cmca Cmmm Cccm Cmma Ccca Fmmm Fddd Immm Ibam Ibca Imma P4 P4(1) P4(2) P4(3) I4 I4(1) P-4 I-4 P4/m P4(2)/m P4/n P4(2)/n I4/m I4(1)/a P422 P42(1)2 P4(1)22 P4(1)2(1)2 P4(2)22 P4(2)2(1)2 P4(3)22 P4(3)2(1)2 I422 I4(1)22 P4mm P4bm P4(2)cm P4(2)nm P4cc P4nc P4(2)mc P4(2)bc I4mm I4cm I4(1)md I4(1)cd P-42m P-42c P-42(1)m P-42(1)c P-4m2 P-4c2 P-4b2 P-4n2 I-4m2 I-4c2 I-42m I-42d P4/mmm P4/mcc P4/nbm P4/nnc P4/mbm P4/mnc P4/nmm P4/ncc P4(2)/mmc P4(2)/mcm P4(2)/nbc P4(2)/nnm P4(2)/mbc P4(2)/mnm P4(2)/nmc P4(2)/ncm I4/mmm I4/mcm I4(1)/amd I4(1)/acd P3 P3(1) P3(2) R3 P-3 R-3 P312 P321 P3(1)12 P3(1)21 P3(2)12 P3(2)21 R32 P3m1 P31m P3c1 P31c R3m R3c P-31m P-31c P-3m1 P-3c1 R-3m R-3c P6 P6(1) P6(5) P6(2) P6(4) P6(3) P-6 P6/m P6(3)/m P622 P6(1)22 P6(5)22 P6(2)22 P6(4)22 P6(3)22 P6mm P6cc P6(3)cm P6(3)mc P-6m2 P-6c2 P-62m P-62c P6/mmm P6/mcc P6(3)/mcm P6(3)/mmc P23 F23 I23 P2(1)3 I2(1)3 Pm-3 Pn-3 Fm-3 Fd-3 Im-3 Pa-3 Ia-3 P432 P4(2)32 F432 F4(1)32 I432 P4(3)32 P4(1)32 I4(1)32 P-43m F4-3m I-43m P-43n F-43c I-43d Pm-3m Pn-3n Pm-3n Pn-3m Fm-3m Fm-3c Fd-3m Fd-3c Im-3m Ia-3d", " ", keepempty=false) .|> string


## -------------------------------------------------
## FUNCTIONS
## -------------------------------------------------

@inline strtonum(s) = (s isa Number) ? s : (try parse(Int, s) catch _ parse(Float64, s) end)


# find the bravais lattice name from lattice parameters
function Find_Lattice(a,b,c,alpha,beta,gamma)
    reticolo=""
    if alpha≈90.0 && gamma≈90.0
        if beta≈90.0
            if a≈b && a≈c
                reticolo = "cubic"
            elseif a≈b
                reticolo = "tetragonal"
            else
                reticolo = "orthorhombic";
            end
        else
            reticolo = "monoclinic"
        end
    elseif alpha≈90.0 && beta≈90.0 && gamma≈90.0
        reticolo = "hexagonal"
    elseif alpha≈beta && alpha≈gamma && a≈b && a≈c
        reticolo = "rhombohedral"
    else
        reticolo = "triclinic"
    end
    return reticolo
end


function Find_ibrav(spacegroup::String, reticolo::String)
    ibrav=0
    primitive    = occursin("P",spacegroup)
    bodycentered = occursin("I",spacegroup)
    facecentered = occursin("F",spacegroup)
    basecentered = occursin("C",spacegroup)
    if reticolo=="cubic"
        if (primitive) ibrav=1 end
        if (facecentered) ibrav=2 end
        if (bodycentered) ibrav=3 end
    elseif reticolo=="tetragonal"
        if (primitive) ibrav=6 end
        if (bodycentered) ibrav=7 end
    elseif reticolo=="orthorhombic"
        if (primitive) ibrav=8 end
        if (basecentered) ibrav=9 end
        if (facecentered) ibrav=10 end
        if (bodycentered) ibrav=11 end
    elseif reticolo=="monoclinic"
        if (primitive) ibrav=-12 end
        if (basecentered) ibrav=13 end
    elseif reticolo=="triclinic"
        ibrav=14
    elseif reticolo=="hexagonal"
        ibrav=4
    elseif reticolo=="rhombohedral"
        if (primitive) ibrav=4; else ibrav=5; end
    else
        ibrav=0
    end
    return ibrav
end


function celldm(ibrav)
    # https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm200
    if ibrav in [4,6,7]
        return (@sprintf  "                   celldm(1) = %19.14f, celldm(3) = %19.14f" a/_BOHR_ c/a)
    elseif ibrav in [5,-5]
        return (@sprintf  "                   celldm(1) = %19.14f, celldm(4) = %19.14f" a/_BOHR_ cos(alphar))
    elseif ibrav in [14,]
        return (@sprintf  "                   celldm(1) = %19.14f, celldm(2) = %19.14f, celldm(3) = %19.14f, celldm(4) = %19.14f, celldm(5) = %19.14f, celldm(6) = %19.14f"  a/_BOHR_  b/a  c/a  cos(alphar) cos(betar) cos(gammar))
    elseif ibrav in [-12,-13]
        return (@sprintf  "                   celldm(1) = %19.14f, celldm(2) = %19.14f, celldm(3) = %19.14f, celldm(5) = %19.14f" a/_BOHR_  b/a  c/a  cos(betar))
    elseif ibrav in [12,13]
        return (@sprintf  "                   celldm(1) = %19.14f, celldm(2) = %19.14f, celldm(3) = %19.14f, celldm(4) = %19.14f" a/_BOHR_  b/a  c/a  cos(gammar))
    elseif ibrav in [8,9,-9,91,10,11]
        return (@sprintf  "                   celldm(1) = %19.14f, celldm(2) = %19.14f, celldm(3) = %19.14f"  a/_BOHR_  b/a  c/a)
    elseif ibrav in [0,1,2,3]
        return ""
    else
        return ""
    end
end


function generate_kpoints() end


## ----------------------------------------------------------------
## SECTIONS
## ----------------------------------------------------------------


function CELL_PARAMETERS(a, b, c, alphar, betar, gammar)
    cell_px = [a,   b*cos(gammar), c*cos(betar)]
    cell_py = [0.0, b*sin(gammar), c*(cos(alphar)-cos(betar)*cos(gammar))/sin(gammar)]
    cell_pz = [0.0, 0.0, c*sqrt(1.0 - cos(alphar)^2 - cos(betar)^2 - cos(gammar)^2 + 2*cos(alphar)*cos(betar)*cos(gammar))/sin(gammar)]
    lines = ["\nCELL_PARAMETERS {bohr}", ]
    for i=1:3
        s = (@sprintf "  %19.14f   %19.14f   %19.14f"  cell_px[i]/_BOHR_  cell_py[i]/_BOHR_  cell_pz[i]/_BOHR_)
        push!(lines, s)
    end
    push!(lines, "\n\n")
    return lines
end


CONTROL(TITLE, CALC_TYPE, OUT_FD, PSEUDO_FD; FILE_PREFIX="pwscf", TSTRESS=false, DISK_IO="none", F_THR=0.001, EN_THR=0.0001, VERBOSITY="default") = [
"&CONTROL" ,
"                       title = '" * TITLE * "'" ,
"                 calculation = '" * CALC_TYPE * "'" ,
"                restart_mode = 'from_scratch'" ,
"                      outdir = '" * OUT_FD * "'" ,
"                  pseudo_dir = '" * PSEUDO_FD * "'" ,
"                      prefix = '" * FILE_PREFIX * "'" ,
"                     disk_io = '$DISK_IO'" ,
"                   verbosity = '$VERBOSITY'" ,
"               etot_conv_thr = $EN_THR" ,     ## Convergence threshold on total energy (a.u) for ionic minimization
"               forc_conv_thr = $F_THR" ,      ## Convergence threshold on forces (a.u) for ionic minimization
"                     tstress = ." * (TSTRESS ? "true" : "false") * "." ,
" /",
]


ELECTRONS(MIX_BETA; DIAG_METHOD="david", MIX_NDIM=8, ELECTRON_STEP=200, CONVERGENCE=1e-7, DIAG_THR=1e-4) = [
" &ELECTRONS",
"            electron_maxstep = $ELECTRON_STEP",
"                    conv_thr = $CONVERGENCE",
"              diago_thr_init = $DIAG_THR",
"                 startingpot = 'atomic'",
"                 startingwfc = 'atomic'",
"                 mixing_mode = 'plain'",
"                 mixing_beta = $MIX_BETA",
"                 mixing_ndim = $MIX_NDIM",
"             diagonalization = '$DIAG_METHOD'",
" /",
]


SYSTEM(ibrav, SG::Int, NATOM::Int, NTYPE::Int, ecutwfc::Int, ecutrho::Int; VdW_CORR_MODE="none", xdm_a1=1.2153, xdm_a2=2.3704) = [
"&SYSTEM",
"                       ibrav = $ibrav",
celldm(ibrav),
"                 space_group = $SG",
"                         nat = $NATOM",    ## number of atoms in the unit cell (ALL atoms, except if space_group is set, in which case, INEQUIVALENT atoms),
"                        ntyp = $NTYPE",    ## number of types of atoms in the unit cell,
"                     ecutwfc = $ecutwfc",  # kinetic energy cutoff (Ry) for wavefunctions,
"                     ecutrho = $ecutrho",  # Kinetic energy cutoff (Ry) for charge density and potential,
"                    vdw_corr = '$VdW_CORR_MODE'", ## xdm => Exchange-hole dipole-moment model.,
"                      xdm_a1 = $xdm_a1",   ## 1.2153,
"                      xdm_a2 = $xdm_a2",   ## 2.3704",
" /",
]


IONS(M) = [
"&IONS",
(length(M)>0) ? "                ion_dynamics = '$M'" : "",
" /",
]


CELL(M) = [
"&CELL",
" /",
]


function ATOMIC_SPECIES(pseudo_dic::Dict)
    lines = ["\nATOMIC_SPECIES",]
    for (atm, pseudo)  in  pseudo_dic
        @assert (atm in AtomSymb)
        s = (@sprintf "  %3s  %14.10f  %s"  atm  Atoms[atm]  pseudo)
        push!(lines, s)
    end
    return lines
end


function ATOMIC_POSITIONS(pos)
    lines = ["\nATOMIC_POSITIONS crystal"]
    for (atm,x,y,z) in pos
        s = (@sprintf "%3s   %19.14f   %19.14f   %19.14f"  atm  x  y  z)
        push!(lines, s)
    end
    return lines
end


function ATOMIC_POSITIONS_SG(pos_sg)
    lines = ["\nATOMIC_POSITIONS crystal_sg"]
    for (atm,x,y,z) in pos_sg
        s = (@sprintf "%3s   %19.14f   %19.14f   %19.14f"  atm  x  y  z)
        push!(lines, s)
    end
    return lines
end


function K_POINTS(mode="gamma")
    lines = ["\nK_POINTS $mode",]
    return lines
end


## =======================================================================

function pwx_inupt_compute_energy(TITLE, CIF, PSEUDO_FILES, POSITIONS, OUT_FD, PSEUDO_FD)
    a, b, c  = strtonum.([CIF["_cell_length_a"], CIF["_cell_length_b"], CIF["_cell_length_c"]])
    alpha, beta, gamma  = strtonum.([CIF["_cell_angle_alpha"], CIF["_cell_angle_beta"], CIF["_cell_angle_gamma"]])
    alphar   = π * alpha / 180.0
    betar    = π * beta  / 180.0
    gammar   = π * gamma / 180.0
    IT_num   = CIF["_symmetry_Int_Tables_number"]
    @assert IT_num > 0 
    SG       = Int_Tables[IT_num]
    latt     = Find_Lattice(a,b,c,alpha,beta,gamma)
    ibrav    = Find_ibrav(SG, latt)
    natom    = length(POSITIONS)
    ntype    = length(PSEUDO_FILES) 

    return join(vcat(   CONTROL(TITLE, "scf", OUT_FD, PSEUDO_FD),
                        SYSTEM(ibrav, IT_num, natom, ntype, 120, 600),
                        ELECTRONS(0.5),
                        IONS(""),
                        CELL(""),
                        ATOMIC_SPECIES(PSEUDO_FILES),
                        ATOMIC_POSITIONS(POSITIONS),
                        K_POINTS("gamma"),
                        ibrav==0 ? CELL_PARAMETERS(a, b, c, alphar, betar, gammar) : String[]
            ), "\n" )
end

## =======================================================================

to_cif(a::Float64,b::Float64,c::Float64,alpha::Float64,beta::Float64,gamma::Float64,IT::Int) = Dict(
    "_cell_length_a"=>a,
    "_cell_length_b"=>b,
    "_cell_length_c"=>c,
    "_cell_angle_alpha"=>alpha,
    "_cell_angle_beta"=>beta,
    "_cell_angle_gamma"=>gamma,
    "_symmetry_Int_Tables_number"=>IT
)


## =======================================================================
pwx_inupt_compute_energy(   "TITLE", 
                            to_cif(1.0,1.0,1.0,90.0,90.0,90.0,220), 
                            Dict("Si"=>"si.udf"),
                            [("Si", 1.0,1.0,1.0)], 
                            "./", 
                            "./" ) |> print


##

