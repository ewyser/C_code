#----------------------------------------------------------------------------------------------------------
function meshSetup(nel,Lx,Ly,Lz)
    # geometry
    nel = [nel,nel]                                               
    L   = [Lx,Ly,ceil(Lz)]
    h   = [L[1]/nel[1],L[1]/nel[1],L[1]/nel[1]]
    # mesh 
    xn  = (0.0-2*h[1]):h[1]:(L[1]+2.0*h[1])
    yn  = (0.0-2*h[2]):h[2]:(L[2]+2.0*h[2])
    zn  = (0.0-2*h[3]):h[3]:(L[3]+2.0*h[3])
    zn  = reverse(zn)
    nno = [length(xn),length(yn),length(zn),length(xn)*length(yn)*length(zn)] 
    nel = [nno[1]-1,nno[2]-1,nno[3]-1,(nno[1]-1)*(nno[2]-1)*(nno[3]-1)]
    nn  = 64
    xn  = (xn'.*ones(typeD,nno[3],1     ))     .*ones(typeD,1,1,nno[2])
    zn  = (     ones(typeD,nno[1],1     )'.*zn).*ones(typeD,1,1,nno[2])
    yn  = (     ones(typeD,nno[3],nno[1]))     .*reshape(yn,1,1,nno[2])
    xn  = vec(xn)
    yn  = vec(yn)
    zn  = vec(zn)
    # nodal quantities
    mn  = zeros(typeD,nno[4],1) 
    fen = zeros(typeD,nno[4],3) 
    fin = zeros(typeD,nno[4],3)
    fn  = zeros(typeD,nno[4],3)
    an  = zeros(typeD,nno[4],3)
    pn  = zeros(typeD,nno[4],3)
    vn  = zeros(typeD,nno[4],3)
    un  = zeros(typeD,nno[4],3)
    # mesh-to-node topology
    e2n = e2N(nno,nel,nn)
    # boundary conditions
    xB  = [minimum(xn)+2*h[1],maximum(xn)-2*h[1],minimum(yn)+2*h[1],maximum(yn)-2*h[1],0.0,Inf]                                    
    bcx = vcat(findall(x->x<=xB[1], xn),findall(x->x>=xB[2], xn))
    bcy = vcat(findall(x->x<=xB[3], yn),findall(x->x>=xB[4], yn))
    bcz = findall(x->x<=xB[5], zn)
    # push to struct
    meD = mesh(nel,nno,nn,L,h,xn,yn,zn,mn,fen,fin,fn,an,pn,vn,un,e2n,xB)
    bc  = boundary(bcx,bcy,bcz)
    return(meD,bc)
end
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
function pointSetup(meD,ni,lz,coh0,cohr,phi0,phir,rho0,nstr,typeD)
    # mpm initialization
    xL  = meD.xB[1]+(0.5*meD.h[1]/ni):meD.h[1]/ni:meD.xB[2]
    yL  = meD.xB[3]+(0.5*meD.h[2]/ni):meD.h[2]/ni:meD.xB[4]
    zL  = meD.xB[5]+(0.5*meD.h[3]/ni):meD.h[3]/ni:lz-0.5*meD.h[3]/ni
    npx = length(xL)
    npy = length(yL)
    npz = length(zL)
    xp  = vec((xL'.*ones(npz,1  )      ).*ones(1,1,npy))
    yp  = vec((     ones(npz,npx)      ).*reshape(yL,1,1,npy))
    zp  = vec((     ones(npx,1  )'.*zL ).*ones(1,1,npy))    
    wl  = 0.15*lz
    x   = LinRange(minimum(xp),maximum(xp),200)
    a   = -1.25
    z   = a.*x
    x   = x.+0.5.*meD.L[1]

    xlt = xp
    ylt = yp
    zlt = yp
    
    xlt= Float64[]
    ylt= Float64[]
    zlt= Float64[]
    pos= Float64 
    for mp in 1:length(xp)
        for p in 1:length(z)
            Δx = xp[mp]-x[p]
            Δz = zp[mp]-z[p]
            nx = a
            nz = -1.0
            s  = Δx*nx+Δz*nz        
            if(s>0)
                pos = 1
            else
                pos = 0
            end
            if(zp[mp]<wl) 
                pos = 1
            end
        end
        if(pos==1)
            push!(xlt, xp[mp]) # push!(inArray, What), incremental construction of an array of arbitrary size
            push!(ylt, yp[mp]) # push!(inArray, What), incremental construction of an array of arbitrary size
            push!(zlt, zp[mp]) # push!(inArray, What), incremental construction of an array of arbitrary size
        end
    end
    # material point's quantities
    # scalars & vectors
    nmp  = length(xlt)
    l0   =  ones(typeD,nmp,3).*0.5.*(meD.h[1]./ni)
    l    =  ones(typeD,nmp,3).*0.5.*(meD.h[1]./ni)
    v0   =  ones(typeD,nmp,1).*(2.0.*l0[:,1].*2.0.*l0[:,2].*2.0.*l0[:,3])
    v    =  ones(typeD,nmp,1).*(2.0.*l[:,1].*2.0.*l[:,2].*2.0.*l[:,3])
    m    = rho0.*v0
    xp   = hcat(xlt,ylt,zlt)
    up   = zeros(typeD,nmp,3) 
    vp   = zeros(typeD,nmp,3)
    pp   = zeros(typeD,nmp,3)
    coh  =  ones(typeD,nmp,1).*coh0
    cohr =  ones(typeD,nmp,1).*cohr
    phi  =  ones(typeD,nmp,1).*phi0
    p    = findall(x->x<=2*wl, xp[:,3])
    phi[p] .= phir

    epII = zeros(typeD,nmp,1)
    # tensors
    dF   = zeros(typeD,nmp,9) 
    F    = zeros(typeD,nmp,9)
    b    = zeros(typeD,nmp,9)
    bT   = zeros(typeD,nmp,9)
    e    = zeros(typeD,nstr,nmp)
    s    = zeros(typeD,nstr,nmp)
    ep   = zeros(typeD,nstr,nmp)
    # additional quantities
    nn   = convert(UInt64,meD.nn)
    S    = zeros(typeD,nmp ,nn       )
    dSx  = zeros(typeD,nmp ,nn       )
    dSy  = zeros(typeD,nmp ,nn       )
    dSz  = zeros(typeD,nmp ,nn       )
    B    = zeros(typeD,nstr,nn.*2,nmp)
    # connectivity
    p2e  = zeros(UInt64,nmp,1)
    p2n  = zeros(UInt64,nmp,nn)
    # push to struct
    mpD  = point(nmp,l0,l,v0,v,m,xp,up,vp,pp,coh,cohr,phi,epII,
                 dF,F,b,bT,e,s,ep,S,dSx,dSy,dSz,B,p2e,p2n)
    return(mpD)
end
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
function e2N(nno,nel,nn)
	gnum = reverse(reshape(1:(nno[4]),nno[3],nno[1],nno[2]),dims=1)
	e2n  = zeros(nel[4],nn)
	iel  = 1
	for k in 1:nel[2]#nely
		for i in 1:nel[1]#nelx
            for j in 1:nel[3]#nelz
			    if(i>1 && i<nel[1] && j>1 && j<nel[3] && k>1 && k<nel[2])
                    e2n[iel,1 ] = gnum[j-1,i-1,k-1]
                    e2n[iel,2 ] = gnum[j-0,i-1,k-1]
                    e2n[iel,3 ] = gnum[j+1,i-1,k-1]
                    e2n[iel,4 ] = gnum[j+2,i-1,k-1]
                    e2n[iel,5 ] = gnum[j-1,i  ,k-1]
                    e2n[iel,6 ] = gnum[j-0,i  ,k-1]
                    e2n[iel,7 ] = gnum[j+1,i  ,k-1]
                    e2n[iel,8 ] = gnum[j+2,i  ,k-1]
                    e2n[iel,9 ] = gnum[j-1,i+1,k-1]
                    e2n[iel,10] = gnum[j-0,i+1,k-1]
                    e2n[iel,11] = gnum[j+1,i+1,k-1]
                    e2n[iel,12] = gnum[j+2,i+1,k-1]
                    e2n[iel,13] = gnum[j-1,i+2,k-1]
                    e2n[iel,14] = gnum[j-0,i+2,k-1]
                    e2n[iel,15] = gnum[j+1,i+2,k-1]
                    e2n[iel,16] = gnum[j+2,i+2,k-1]
                    
                    e2n[iel,17] = gnum[j-1,i-1,k  ]
                    e2n[iel,18] = gnum[j-0,i-1,k  ]
                    e2n[iel,19] = gnum[j+1,i-1,k  ]
                    e2n[iel,20] = gnum[j+2,i-1,k  ]
                    e2n[iel,21] = gnum[j-1,i  ,k  ]
                    e2n[iel,22] = gnum[j-0,i  ,k  ]
                    e2n[iel,23] = gnum[j+1,i  ,k  ]
                    e2n[iel,24] = gnum[j+2,i  ,k  ]
                    e2n[iel,25] = gnum[j-1,i+1,k  ]
                    e2n[iel,26] = gnum[j-0,i+1,k  ]
                    e2n[iel,27] = gnum[j+1,i+1,k  ]
                    e2n[iel,28] = gnum[j+2,i+1,k  ]
                    e2n[iel,29] = gnum[j-1,i+2,k  ]
                    e2n[iel,30] = gnum[j-0,i+2,k  ]
                    e2n[iel,31] = gnum[j+1,i+2,k  ]
                    e2n[iel,32] = gnum[j+2,i+2,k  ]
                    
                    e2n[iel,33] = gnum[j-1,i-1,k+1]
                    e2n[iel,34] = gnum[j-0,i-1,k+1]
                    e2n[iel,35] = gnum[j+1,i-1,k+1]
                    e2n[iel,36] = gnum[j+2,i-1,k+1]
                    e2n[iel,37] = gnum[j-1,i  ,k+1]
                    e2n[iel,38] = gnum[j-0,i  ,k+1]
                    e2n[iel,39] = gnum[j+1,i  ,k+1]
                    e2n[iel,40] = gnum[j+2,i  ,k+1]
                    e2n[iel,41] = gnum[j-1,i+1,k+1]
                    e2n[iel,42] = gnum[j-0,i+1,k+1]
                    e2n[iel,43] = gnum[j+1,i+1,k+1]
                    e2n[iel,44] = gnum[j+2,i+1,k+1]
                    e2n[iel,45] = gnum[j-1,i+2,k+1]
                    e2n[iel,46] = gnum[j-0,i+2,k+1]
                    e2n[iel,47] = gnum[j+1,i+2,k+1]
                    e2n[iel,48] = gnum[j+2,i+2,k+1]
                     
                    e2n[iel,49] = gnum[j-1,i-1,k+2]
                    e2n[iel,50] = gnum[j-0,i-1,k+2]
                    e2n[iel,51] = gnum[j+1,i-1,k+2]
                    e2n[iel,52] = gnum[j+2,i-1,k+2]
                    e2n[iel,53] = gnum[j-1,i  ,k+2]
                    e2n[iel,54] = gnum[j-0,i  ,k+2]
                    e2n[iel,55] = gnum[j+1,i  ,k+2]
                    e2n[iel,56] = gnum[j+2,i  ,k+2]
                    e2n[iel,57] = gnum[j-1,i+1,k+2]
                    e2n[iel,58] = gnum[j-0,i+1,k+2]
                    e2n[iel,59] = gnum[j+1,i+1,k+2]
                    e2n[iel,60] = gnum[j+2,i+1,k+2]
                    e2n[iel,61] = gnum[j-1,i+2,k+2]
                    e2n[iel,62] = gnum[j-0,i+2,k+2]
                    e2n[iel,63] = gnum[j+1,i+2,k+2]
                    e2n[iel,64] = gnum[j+2,i+2,k+2]
			    end
                iel = iel+1;
            end

		end
	end
	return(convert(Array{Int64},e2n))
end
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
function save2txt(meD,mpD,bc)
    xn = [meD.xn meD.yn meD.zn]
    writedlm("xn.txt",vec(xn))
    writedlm("e2n.txt",vec(meD.e2n.-1))
    
    writedlm("cohp.txt",vec(mpD.coh))
    writedlm("phip.txt",vec(mpD.phi))

    p = [mpD.nmp meD.nn meD.nno[4] meD.h[1] meD.h[2] meD.h[3] minimum(meD.xn) minimum(meD.yn) minimum(meD.zn) meD.nno[1] meD.nno[2] meD.nno[3]]
    writedlm("param.txt",vec(p))    

    writedlm("mp.txt" ,vec(mpD.m) )    
    writedlm("xp.txt" ,vec(mpD.xp))    
    writedlm("vol.txt",vec(mpD.v) )    
    writedlm("lp.txt" ,vec(mpD.l)) 


    bcx = bc.x.+(0*meD.nno[4])
    bcy = bc.y.+(1*meD.nno[4])
    bcz = bc.z.+(2*meD.nno[4])
    BC  = ones(Int64,meD.nno[4]*3,1)
    BC[vcat(bcx,bcy,bcz)].= 0
    writedlm("bcs.txt" ,vec(BC)) 

    return("done")
end