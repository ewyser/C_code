# Initialisation
using Printf, LinearAlgebra, DelimitedFiles

typeD = Float64  # Precision (double=Float64 or single=Float32)

include("../src/functions.jl")
include("../src/types.jl")

@views function main()
	# ---------------------------------------------------------------------------
	nel = 40
	@printf("\n------------------------")
	@printf("\nRun: nel = %.0f",nel)
	@printf("\n------------------------")
	# ---------------------------------------------------------------------------
	# non-dimensional constant 
	# ---------------------------------------------------------------------------
	ν       = 0.3                                                        
    ni      = 2
    nstr    = 6                                                                                                                      
	# ---------------------------------------------------------------------------
	# independant physical constant
	# ---------------------------------------------------------------------------
    g       = 9.81                                                        # gravitationnal acceleration [m/s^2]
    E       = 1.0e6                                                       # young's modulus             [Pa]
    Gc      = E/(2.0*(1.0+ν))                                                # shear modulus               [Pa]
    Kc      = E/(3.0*(1.0-2.0*ν))                                              # bulk modulus                [Pa]
    ρ0      = 2700.0                                                        # density                     [kg/m^3]
    yd      = sqrt((Kc+4.0/3.0*Gc)/ρ0)                                      # elastic wave speed          [m/s]
    c0      = 20.0e3                                                      # cohesion                    [Pa]
    ϕ0      = 20.0*pi/180                                                 # friction angle              [Rad]
    ψ0      = 0.0                                                          # dilantancy angle
    H       = -60.0e3                                                       # softening modulus           [Pa]
    cr      =   4.0e3                                                      # residual cohesion           [Pa]
    ϕr      = 7.0*pi/180                                                  # residual friction angle     [Rad]
    t       = 15.0                                                        # simulation time             [s]
    te      = 10.0                                                         # elastic loading             [s]
    tg      = te/1.5                                                      # gravity increase 
	# ---------------------------------------------------------------------------
    # mesh & mp setup
    # ---------------------------------------------------------------------------
    lx      = 64.1584
    ly      = lx/4
    lz      = 12.80
    meD,bc  = meshSetup(nel,lx,ly,lz)	
    mpD     = pointSetup(meD,ni,lz,c0,cr,ϕ0,ϕr,ρ0,nstr,typeD)
    # isotropic elastic matrix
    Del     = [ Kc+4/3*Gc Kc-2/3*Gc Kc-2/3*Gc 0.0 0.0 0.0;
                Kc-2/3*Gc Kc+4/3*Gc Kc-2/3*Gc 0.0 0.0 0.0;
                Kc-2/3*Gc Kc-2/3*Gc Kc+4/3*Gc 0.0 0.0 0.0;
                0.0       0.0       0.0       Gc  0.0 0.0;                             
                0.0       0.0       0.0       0.0 Gc  0.0;
                0.0       0.0       0.0       0.0 0.0 Gc ] 
    Hp      = H*meD.h[1]                                                  
    # ---------------------------------------------------------------------------
    # display parameters & runtime
    # ---------------------------------------------------------------------------
    C       = 0.5                                                            #
    dt      = C*meD.h[1]/yd                                                  # unconditionally stable timestep
    nit     = ceil(t/dt)                                                     # maximum number of interation
    nf      = max(2,ceil(round(1/dt)/25))                                    # number of frame interval
    # runtime parameters
    it      = 1                                                              # initialize iteration
    tw      = 0.0                                                            # initialize time while statement
    
    char    = save2txt(meD,mpD,bc)
    p       = [g;ρ0;ψ0;ν;E;Kc;Gc;cr;Hp;t;te;tg]
    writedlm("phys.txt" ,p)

    # ---------------------------------------------------------------------------
    # explicit solver for finite deformation mechanics
    # ---------------------------------------------------------------------------
    @printf("\nmpm solver: %.0f elements"       ,meD.nel[4])
    @printf("\n          : %.0f nodes"          ,meD.nno[4])
    @printf("\n          : %.0f material points\n", mpD.nmp                     )
end
main()