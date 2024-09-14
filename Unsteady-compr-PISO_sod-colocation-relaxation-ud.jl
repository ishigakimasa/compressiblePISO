# 
# - i=1,2の間が境界となるように設定。 
#

using Parameters
using HDF5

@with_kw struct Param{T1,T2}
    dt::T1 = 1e-5 # s

    # mesh size
    nx::T2 = 1000
    ny::T2 = 20
    
    # length
    Lx::T1 = 10.0
    Ly::T1 = 1.0
    dx::T1 = Lx/nx
    dy::T1 = Ly/ny

    # viscosity
    μ0::T1 = 1e-5 # Pa s
    k0::T1 = 0.0241 # W/Km
    M::T1 = 29.0e-3 # molucular mass 29x10^-3 kg/mol
    R0::T1 = 8.31 # ideal gas constant J/K mol
    R::T1 = R0/M
    γ::T1 = 1.4 # specific heat ratio
    Cv::T1 = R/(γ -1.0) # J/kgK

    nstep::T2 = 1000
    niter::T2 = 25
    ipter::T2 = 200
    alphaua::T1 = 0.5
    alphava::T1 = 0.5
    alphapa::T1 = 0.5
    alphaea::T1 = 0.5
    
    ndata::T2 = 200

    outfile::String = "Data-ud-01/uvTp" 

end

function main(para)
    @unpack nx, ny, nstep, niter, ndata, dx, dy = para

    μ = zeros(Float64, nx+2, ny+2)
    k = zeros(Float64, nx+2, ny+2)

    # auxiliary velocity
    ua = zeros(Float64, nx+2, ny+2)
    utemp = zeros(Float64, nx+2, ny+2)
    utemp2 = zeros(Float64, nx+2, ny+2)
    uold = zeros(Float64, nx+2, ny+2)
    apua = zeros(Float64, nx+2, ny+2)
    aeua = zeros(Float64, nx+2, ny+2)
    awua  = zeros(Float64, nx+2, ny+2)
    anua = zeros(Float64, nx+2, ny+2)
    asua = zeros(Float64, nx+2, ny+2)
    bua  = zeros(Float64, nx+2, ny+2)
    buae = zeros(Float64, nx+2, ny+2)
    buaw = zeros(Float64, nx+2, ny+2) 
    buan = zeros(Float64, nx+2, ny+2) 
    buas = zeros(Float64, nx+2, ny+2)
    
    va = zeros(Float64, nx+2, ny+2)
    vtemp = zeros(Float64, nx+2, ny+2)
    vtemp2 = zeros(Float64, nx+2, ny+2)
    vold = zeros(Float64, nx+2, ny+2)
    apva = zeros(Float64, nx+2, ny+2)
    aeva = zeros(Float64, nx+2, ny+2)
    awva  = zeros(Float64, nx+2, ny+2)
    anva = zeros(Float64, nx+2, ny+2)
    asva = zeros(Float64, nx+2, ny+2)
    bva = zeros(Float64, nx+2, ny+2)
    bvae = zeros(Float64, nx+2, ny+2)
    bvaw = zeros(Float64, nx+2, ny+2)
    bvan = zeros(Float64, nx+2, ny+2)
    bvas = zeros(Float64, nx+2, ny+2)
    
    Ueq = zeros(Float64, nx+2, ny+2)
    Veq = zeros(Float64, nx+2, ny+2)

    uf = zeros(Float64, nx+2, ny+2)
    vf = zeros(Float64, nx+2, ny+2)

    ea = zeros(Float64, nx+2, ny+2)
    eold = zeros(Float64, nx+2, ny+2)
    etemp = zeros(Float64, nx+2, ny+2)
    apea = zeros(Float64, nx+2, ny+2)
    aeea = zeros(Float64, nx+2, ny+2)
    awea  = zeros(Float64, nx+2, ny+2)
    anea = zeros(Float64, nx+2, ny+2)
    asea = zeros(Float64, nx+2, ny+2)
    bea = zeros(Float64, nx+2, ny+2)

    ρa = zeros(Float64, nx+2, ny+2)
    ρold = zeros(Float64, nx+2, ny+2)
    Ta = zeros(Float64, nx+2, ny+2)
    Told = zeros(Float64, nx+2, ny+2)


    # pressure, corrected pressure
    pa = zeros(Float64, nx+2, ny+2)
    pc = zeros(Float64, nx+2, ny+2)
    ptemp = zeros(Float64, nx+2, ny+2)
    pold = zeros(Float64, nx+2, ny+2)
    appc = zeros(Float64, nx+2, ny+2)
    aepc = zeros(Float64, nx+2, ny+2)
    awpc = zeros(Float64, nx+2, ny+2)
    anpc = zeros(Float64, nx+2, ny+2)
    aspc = zeros(Float64, nx+2, ny+2)
    bpc = zeros(Float64, nx+2, ny+2)

    # coordinate
    x = zeros(Float64, nx+1, ny+1)
    y = zeros(Float64, nx+1, ny+1)
    for j = 2:ny+1
        for i = 2:nx+1
            x[i,j] = (i-2 + 0.5)*dx
        end
    end
    for j = 2:ny+1
        for i = 2:nx+1
            y[i,j] = (j-2 + 0.5)*dy
        end
    end

    rest = 0.0
    init(para, x, y, pa, ea, Ta, ρa, ua, va)

    for l = 1:nstep
        #println(l, "===========================================")

        for j = 1:ny+2
            for i = 1:nx+2
                uold[i,j] = ua[i,j]
                vold[i,j] = va[i,j]
                eold[i,j] = ea[i,j]
                ρold[i,j] = ρa[i,j]
                Told[i,j] = Ta[i,j]
                pold[i,j] = pa[i,j]
            end
        end

        BC(para, ua, va, pa, pc, ea, ρa, l)
        updateTransport(para, k, μ)
        
        for n = 1:niter
            BC(para, ua, va, pa, pc, ea, ρa, l)

            makeae(para, apea, aeea, awea, anea, asea, ρa, ua, va, k, ρold)
            solveEA(para, ua, va, pa, ea, eold, ρa, ρold, etemp, bea, buae, buaw, buan, buas, bvae, bvaw, bvan, bvas, apea, aeea, awea, anea, asea, k, μ)

            makeau(para, apua, aeua, awua, anua, asua, ρa, ua, va, μ, ρold)
            makeav(para, apva, aeva, awva, anva, asva, ρa, ua, va, μ, ρold)
            makebuv(para, ua, va, pa, bua, buae, buaw, buan, buas, uold, bva, bvae, bvaw, bvan, bvas, vold, ρold, μ)
            makeUVeq(para, ua, va, pa, apua, aeua, awua, anua, asua, bua,  apva, aeva, awva, anva, asva, bva, Ueq, Veq)
            solveUA(para, Ueq, ua, pa, utemp, apua)            
            solveVA(para, Veq, va, pa, vtemp, apva)
            RhieChow(para, ua, va, pa, uf, vf, apua, apva, l)

            convT(para, ua, va, ea, Ta)
            makeap(para, appc, aepc, awpc, anpc, aspc, apua, apva, ρa, Ta)
            solvePC(para, ρa, pc, ptemp, appc, aepc, awpc, anpc, aspc, bpc, uf, vf, Ta, ρold, pa)
            update(para, ua, va, pc, apua, apva, ρa, pa, Ta, )

            makeae(para, apea, aeea, awea, anea, asea, ρa, ua, va, k, ρold)
            solveEA(para, ua, va, pa, ea, eold, ρa, ρold, etemp, bea, buae, buaw, buan, buas, bvae, bvaw, bvan, bvas, apea, aeea, awea, anea, asea, k, μ)

            makeau(para, apua, aeua, awua, anua, asua, ρa, ua, va, μ, ρold)
            makeav(para, apva, aeva, awva, anva, asva, ρa, ua, va, μ, ρold)
            makebuv(para, ua, va, pa, bua, buae, buaw, buan, buas, uold, bva, bvae, bvaw, bvan, bvas, vold, ρold, μ)
            convT(para, ua, va, ea, Ta)
            makeap(para, appc, aepc, awpc, anpc, aspc, apua, apva, ρa, Ta)
            makeUVeq(para, ua, va, pa, apua, aeua, awua, anua, asua, bua,  apva, aeva, awva, anva, asva, bva, Ueq, Veq)
            solvePA2(para, pa, apua, apva, appc, aepc, awpc, anpc, aspc, bpc, ptemp, ρold, Ueq, Veq, ρa, l)
            update2(para, ua, va, utemp, vtemp, pa, apua, aeua, awua, anua, asua, apva, aeva, awva, anva, asva, bua, bva, ρa, Ta)
            makeae(para, apea, aeea, awea, anea, asea, ρa, ua, va, k, ρold)
            solveEA(para, ua, va, pa, ea, eold, ρa, ρold, etemp, bea, buae, buaw, buan, buas, bvae, bvaw, bvan, bvas, apea, aeea, awea, anea, asea, k, μ)

            rest = converg(para, ua, va, utemp2, vtemp2)
            #println("iter, redisual ", n,", ",rest)
            for j=2:ny+1
                for i=2:nx+1
                    utemp2[i,j] = ua[i,j]
                    vtemp2[i,j] = va[i,j]
                end
            end
            
            if rest <=1e-7 && n>=2
                if l%10 == 0
                    println(l, " th step, ", n, " th iteration, residual=", rest)
                end
                break
            end

            BC(para, ua, va, pa, pc, ea, ρa, l)
        end
        
        if l%ndata == 0
            println("residual=", rest)
            convT(para, ua, va, ea, Ta)
            output(para, ua, va, pa, Ta, ea, ρa,  x, y, l)
        end
    end
    #streamfunction(para, u, v, psi)
    
end

function init(para, x, y, pa, ea, Ta, ρa, ua, va)
    @unpack nx, ny, dx, dy, R, Cv, Lx = para

    #for j=1:ny+2
    #    for i=1:nx+2
    #        pa[i,j] = 1e5 + 1e4*cos(2π/(nx+1)*(i-1))
    #    end
    #end
    #T0 = 348.4
    T0 = 348.4 
    for j = 1:ny+2
        for i = 1:div(nx+2, 2)
            pa[i,j] = 1e5
            Ta[i,j] = T0
            ea[i,j] = Cv*Ta[i,j] #+ 0.5*(ua[i,j]^2 + va[i,j]^2)
            ρa[i,j] = pa[i,j]/(R*Ta[i,j])
        end
    end
    T0 = 278.7 
    for j = 1:ny+2
        for i = div(nx+2, 2)+1:nx+2
            pa[i,j] = 1e4
            Ta[i,j] = T0
            ea[i,j] = Cv*Ta[i,j] #+ 0.5*(ua[i,j]^2 + va[i,j]^2)
            ρa[i,j] = pa[i,j]/(R*Ta[i,j])
        end
    end

end

function BC(para, ua, va, pa, pc, ea, ρa, l)
    @unpack nx, ny, dt = para
    for i = 1:nx+2
        ua[i,ny+2] = ua[i,ny+1]#0.0
        va[i,ny+2] = -va[i,ny+1]#0.0
        ua[i,1] = ua[i,2]#0.0
        va[i,1] = -va[i,2]#0.0
        
        pa[i,1] = pa[i,2]
        pa[i,ny+2] = pa[i,ny+1]
        pc[i,1] = pc[i,2]
        pc[i,ny+2] = pc[i,ny+1]
        ea[i,1] = ea[i,2]
        ea[i,ny+2] = ea[i,ny+1]
        ρa[i,1] = ρa[i,2]
        ρa[i,ny+2] = ρa[i,ny+1]
        
    end
    
    for j = 1:ny+2
        ua[nx+2,j] = -ua[nx+1,j]#0.0
        va[nx+2,j] = va[nx+1,j]#0.0
        ua[1,j] = -ua[2,j]
        va[1,j] = va[2,j]#0.0

        pa[1,j] = pa[2,j]
        pa[nx+2,j] = pa[nx+1,j]
        pc[nx+2,j] = pc[nx+1,j]
        pc[1,j] = pc[2,j]
        ea[1,j] = ea[2,j] 
        ea[nx+2,j] = ea[nx+1,j]
        ρa[1,j] = ρa[2,j]
        ρa[nx+2,j] = ρa[nx+1,j]
    end
end

function updateTransport(para, k, μ)
    @unpack nx, ny, k0, μ0 = para

    for j = 1:ny+2
        for i = 1:nx+2
            k[i,j] = k0
            μ[i,j] = μ0
        end
    end
end

function makeau(para, apua, aeua, awua, anua, asua, ρ, u, v, μ, ρold)
    @unpack nx, ny, dx, dy, dt = para
    for j = 2:ny+1
        for i = 2:nx+1
            aeua[i,j] = 4.0/3.0*0.5*(μ[i+1,j] + μ[i,j])*dy/dx + max(-0.5dy*(ρ[i+1,j]*u[i+1,j] + ρ[i,j]*u[i,j]), 0)
            awua[i,j] = 4.0/3.0*0.5*(μ[i,j] + μ[i-1,j])*dy/dx + max(0.5*dy*(ρ[i,j]*u[i,j] + ρ[i-1,j]*u[i-1,j]), 0)
            asua[i,j] = 0.5(μ[i,j] + μ[i,j-1])*dx/dy + max(0.5dx*(ρ[i,j]*v[i,j] + ρ[i,j-1]*v[i,j-1]), 0)
            anua[i,j] = 0.5(μ[i,j] + μ[i,j+1])*dx/dy + max(-0.5dx*(ρ[i,j+1]*v[i,j+1] + ρ[i,j]*v[i,j]), 0)
            
            apua[i,j] = aeua[i,j] + awua[i,j] + anua[i,j] + asua[i,j] + ρold[i,j]*dx*dy/dt
        end
    end
end

function makeav(para, apva, aeva, awva, anva, asva, ρ, u, v, μ, ρold)
    @unpack nx, ny, dx, dy, dt, = para
    for j=2:ny+1
        for i=2:nx+1
            aeva[i,j] = 0.5*(μ[i+1,j] + μ[i,j])*dy/dx + max(-0.5*dy*(ρ[i+1,j]*u[i+1,j] + ρ[i,j]*u[i,j]), 0)
            awva[i,j] = 0.5*(μ[i,j] + μ[i-1,j])*dy/dx + max(0.5*dy*(ρ[i,j]*u[i,j] + ρ[i-1,j]*u[i-1,j]), 0)
            asva[i,j] = 4.0/3.0*0.5(μ[i,j] + μ[i,j-1])*dx/dy + max(0.5dx*(ρ[i,j]*v[i,j] + ρ[i,j-1]*v[i,j-1]), 0)
            anva[i,j] = 4.0/3.0*0.5(μ[i,j] + μ[i,j+1])*dx/dy + max(-0.5dx*(ρ[i,j+1]*v[i,j+1] + ρ[i,j]*v[i,j]), 0)

            apva[i,j] = aeva[i,j] + awva[i,j] + anva[i,j] + asva[i,j] + ρold[i,j]*dx*dy/dt
        end
    end
end

function makebuv(para, ua, va, pa, bua, buae, buaw, buan, buas, uold, bva, bvae, bvaw, bvan, bvas, vold, ρold, μ, )
    @unpack nx, ny, dx, dy, dt, = para

    for j=2:ny+1
        for i=2:nx+1
            buae[i,j] = -0.25*(2.0/3.0*0.5*(μ[i+1,j] + μ[i,j])*(va[i+1,j+1] - va[i+1,j-1] + va[i,j+1] - va[i,j-1]))
            buaw[i,j] = 0.25*(2.0/3.0*0.5*(μ[i,j] + μ[i-1,j])*(va[i-1,j+1] - va[i-1,j-1] + va[i,j+1] - va[i,j-1]))
            buan[i,j] = 0.25*(0.5*(μ[i,j+1] + μ[i,j])*(va[i+1,j] - va[i-1,j] + va[i+1,j+1] - va[i-1,j+1]))
            buas[i,j] = -0.25*(0.5*(μ[i,j] + μ[i,j-1])*(va[i+1,j] - va[i-1,j] + va[i+1,j-1] - va[i-1,j-1]))
            
            bua[i,j] = dx*dy/dt*ρold[i,j]*uold[i,j] + buae[i,j] + buaw[i,j] + buan[i,j] + buas[i,j]

            bvae[i,j] = 0.25*(0.5*(μ[i,j] + μ[i+1,j])*(ua[i,j+1] - ua[i,j-1] + ua[i+1,j+1] - ua[i+1,j-1]))
            bvaw[i,j] = -0.25*(0.5*(μ[i-1,j] + μ[i,j])*(ua[i,j+1] - ua[i,j-1] + ua[i-1,j+1] - ua[i-1,j-1]))
            bvan[i,j] = -0.25*(2.0/3.0*0.5*(μ[i,j] + μ[i,j+1])*(ua[i+1,j+1] - ua[i-1,j+1] + ua[i+1,j] - ua[i-1,j]))
            bvas[i,j] = 0.25*(2.0/3.0*0.5*(μ[i,j+1] + μ[i,j])*(ua[i+1,j] - ua[i-1,j] + ua[i+1,j-1] - ua[i-1,j-1]))

            bva[i,j] = dx*dy/dt*ρold[i,j]*vold[i,j] + bvae[i,j] + bvaw[i,j] + bvan[i,j] + bvas[i,j]
        end
    end
end

function makeUVeq(para, ua, va, pa, apua, aeua, awua, anua, asua, bua,  apva, aeva, awva, anva, asva, bva, Ueq, Veq)
    @unpack nx, ny,  = para

    for j=2:ny+1
        for i=2:nx+1
            Ueq[i,j] = (aeua[i,j]*ua[i+1,j] + awua[i,j]*ua[i-1,j] + anua[i,j]*ua[i,j+1] + asua[i,j]*ua[i,j-1] + bua[i,j])/apua[i,j]
            Veq[i,j] = (aeva[i,j]*va[i+1,j] + awva[i,j]*va[i-1,j] + anva[i,j]*va[i,j+1] + asva[i,j]*va[i,j-1] + bva[i,j])/apva[i,j]
        end
    end
end

function solveUA(para, Ueq, ua, pa, utemp, apua)
    @unpack nx, ny, dx, dy, dt, alphaua = para

    for j=2:ny+1
        for i=2:nx+1            
            utemp[i,j] = ua[i,j]*(1.0-alphaua) + alphaua*(Ueq[i,j] - 0.5dy*(pa[i+1,j] - pa[i-1,j])/apua[i,j])
        end
    end

    for j=2:ny+1
        for i=2:nx+1
            ua[i,j] = utemp[i,j] # u
        end
    end

    #println("ua=",ua[2:3,2:3])
end

function solveVA(para, Veq, va, pa, vtemp, apva)
    @unpack nx, ny, dx, dy, dt, alphava = para
    
    for j=2:ny+1
        for i=2:nx+1 
           vtemp[i,j] = va[i,j] * (1.0 - alphava) + alphava *(Veq[i,j] - 0.5dx*(pa[i,j+1] - pa[i,j-1])/apva[i,j])
        end
    end

    for j=2:ny+1
        for i=2:nx+1
            va[i,j] = vtemp[i,j]
        end
    end

    #println("va=",va[2:3,2:3])
end

function RhieChow(para, ua, va, pa, uf, vf, apua, apva, l)
    @unpack nx, ny, dx, dy, dt = para

    for j=2:ny+1
        for i=2:nx
            uf[i,j] = 0.5(ua[i+1,j] + ua[i,j]) + 0.5(dy/apua[i,j] + dy/apua[i+1,j])*(pa[i,j] - pa[i+1,j]) -
                0.25dy/apua[i,j]*(pa[i-1,j] - pa[i+1,j]) -0.25dy/apua[i+1,j]*(pa[i,j] - pa[i+2,j])   
        end
        uf[1,j] = 0.0 
        uf[nx+1,j] = 0.0
    end

    for j=2:ny
        for i=2:nx+1
            vf[i,j] = 0.5(va[i,j+1] + va[i,j]) + 0.5(dx/apva[i,j] + dx/apva[i,j+1])*(pa[i,j] - pa[i,j+1]) -
                0.25dx/apva[i,j]*(pa[i,j-1] - pa[i,j+1]) -0.25dx/apva[i,j+1]*(pa[i,j] - pa[i,j+2])   
        end
    end
    for i=2:nx+1
        vf[i,1] = 0.0
        vf[i,ny+1] = 0.0
    end
end

function makeae(para, apea, aeea, awea, anea, asea, ρ, u, v, k, ρold)
    @unpack nx, ny, dx, dy, dt, Cv, γ= para
    for j=2:ny+1
        for i=2:nx+1
            aeea[i,j] = 0.5*(k[i+1,j] + k[i,j])/Cv*dy/dx  + max(-0.5γ*dy*(ρ[i+1,j]*u[i+1,j] + ρ[i,j]*u[i,j]), 0)
            awea[i,j] = 0.5*(k[i,j] + k[i-1,j])/Cv*dy/dx + max(0.5γ*dy*(ρ[i,j]*u[i,j] + ρ[i-1,j]*u[i-1,j]), 0)
            asea[i,j] = 0.5(k[i,j] + k[i,j-1])/Cv*dx/dy + max(0.5γ*dx*(ρ[i,j]*v[i,j] + ρ[i,j-1]*v[i,j-1]), 0)
            anea[i,j] = 0.5(k[i,j] + k[i,j+1])/Cv*dx/dy + max(-0.5γ*dx*(ρ[i,j+1]*v[i,j+1] + ρ[i,j]*v[i,j]), 0)

            apea[i,j] = aeea[i,j] + awea[i,j] + anea[i,j] + asea[i,j] + (γ*ρold[i,j] + (1.0 - γ)ρ[i,j])*dx*dy/dt
        end
    end
end

function solveEA(para, ua, va, pa, ea, eold, ρa, ρold, etemp, bea, buae, buaw, buan, buas, bvae, bvaw, bvan, bvas, apea, aeea, awea, anea, asea, k, μ)
    @unpack nx, ny, dx, dy, dt, alphaea, Cv, γ = para
    for j=2:ny+1
        for i=2:nx+1 
            bea[i,j] = dx*dy/dt*ρold[i,j]*eold[i,j] -
            0.125dy*(γ-1)*(-(ρa[i+1,j]*ua[i+1,j] + ρa[i,j]*ua[i,j])*(ua[i,j]^2 + va[i,j]^2 + ua[i+1,j]^2 + va[i+1,j]^2) + (ρa[i,j]*ua[i,j] + ρa[i-1,j]*ua[i-1,j])*(ua[i,j]^2 + va[i,j]^2 + ua[i-1,j]^2 + va[i-1,j]^2) ) -
            0.125dx*(γ-1)*(-(ρa[i,j+1]*va[i,j+1] + ρa[i,j]*va[i,j])*(ua[i,j]^2 + va[i,j]^2 + ua[i,j+1]^2 + va[i,j+1]^2) + (ρa[i,j]*va[i,j] + ρa[i,j-1]*va[i,j-1])*(ua[i,j]^2 + va[i,j]^2 + ua[i,j-1]^2 + va[i,j-1]^2) )  +
            0.5dy/dx/Cv*(0.5(k[i+1,j] + k[i,j])*((ua[i,j]^2 + va[i,j]^2) -(ua[i+1,j]^2 + va[i+1,j]^2)) - 0.5(k[i,j] + k[i-1,j])*(-(ua[i,j]^2 + va[i,j]^2) + (ua[i-1,j]^2 + va[i-1,j]^2))) + 
            0.5dx/dy/Cv*(0.5(k[i,j+1] + k[i,j])*((ua[i,j]^2 + va[i,j]^2) -(ua[i,j+1]^2 + va[i,j+1]^2)) - 0.5(k[i,j] + k[i,j-1])*(-(ua[i,j]^2 + va[i,j]^2) + (ua[i,j-1]^2 + va[i,j-1]^2))) +
            0.5*(ua[i+1,j] + ua[i,j])*(buae[i,j] + 4.0/3.0*0.5*(μ[i,j] + μ[i+1,j])*(ua[i+1,j] - ua[i,j]))*dy +
            0.5*(ua[i-1,j] + ua[i,j])*(buaw[i,j] - 4.0/3.0*0.5*(μ[i,j] + μ[i-1,j])*(ua[i-1,j] - ua[i,j]))*dy +
            0.5*(ua[i,j+1] + ua[i,j])*(buan[i,j] + 0.5*(μ[i,j] + μ[i,j+1])*(ua[i,j+1] - ua[i,j]))*dx +
            0.5*(ua[i,j-1] + ua[i,j])*(buas[i,j] - 0.5*(μ[i,j] + μ[i,j-1])*(ua[i,j]   - ua[i,j-1]))*dx +
            0.5*(va[i+1,j] + va[i,j])*(bvae[i,j] + 0.5*(μ[i,j] + μ[i+1,j])*(va[i+1,j] - va[i,j]))*dy +
            0.5*(va[i-1,j] + va[i,j])*(bvaw[i,j] - 0.5*(μ[i,j] + μ[i-1,j])*(va[i,j]   - va[i-1,j]))*dy +
            0.5*(va[i,j+1] + va[i,j])*(bvan[i,j] + 4.0/3.0*0.5*(μ[i,j] + μ[i,j+1])*(va[i,j+1] - va[i,j]))*dx +
            0.5*(va[i,j-1] + va[i,j])*(bvas[i,j] - 4.0/3.0*0.5*(μ[i,j] + μ[i,j-1])*(va[i,j]   - va[i,j-1]))*dx
            
            etemp[i,j] = ea[i,j] * (1.0 - alphaea) + alphaea *(aeea[i,j]*ea[i+1,j] + awea[i,j]*ea[i-1,j] + anea[i,j]*ea[i,j+1] + asea[i,j]*ea[i,j-1] + bea[i,j])/apea[i,j]
        end
    end
    for j=2:ny+1
        for i=2:nx+1
            ea[i,j] = etemp[i,j]
        end
    end

    #println("ea,=",ea[2:3,2:3])
    #println("EA ua, va, rhoa, pa,=", ua[502:503,2:3], va[2:3,2:3], ρa[502:503,2:3], pa[502:503,2:3])
end

function makeap(para, appc, aepc, awpc, anpc, aspc, apua, apva, ρ, T)
    @unpack nx, ny, dx, dy, dt, R = para

    for j=2:ny+1
        for i=2:nx+1
            if i == 2
                awpc[i,j] = 0.0
            else
                awpc[i,j] = 0.5dy*(dy*ρ[i-1,j]/apua[i-1,j] + dy*ρ[i,j]/apua[i,j])
            end
            
            if i == nx+1
                aepc[i,j] = 0.0
            else          
                aepc[i,j] = 0.5dy*(dy*ρ[i,j]/apua[i,j] + dy*ρ[i+1,j]/apua[i+1,j])
            end
            
            if j == 2
                aspc[i,j] = 0.0
            else
                aspc[i,j] = 0.5dx*(dx*ρ[i,j-1]/apva[i,j-1] + dx*ρ[i,j]/apva[i,j])
            end
                    
            if j == ny+1
                anpc[i,j] = 0.0
            else
                anpc[i,j] = 0.5dx*(dx*ρ[i,j]/apva[i,j] + dx*ρ[i,j+1]/apva[i,j+1])
            end
                        
            appc[i,j] = aepc[i,j] + awpc[i,j] + anpc[i,j] + aspc[i,j] + dx*dy/R/T[i,j]/dt
        end
    end
end

function solvePC(para, ρa, pc, ptemp, appc, aepc, awpc, anpc, aspc, bpc, uf, vf, Ta, ρold, pa,)
    @unpack nx, ny, dx, dy, dt, ipter, alphapa, R = para
    
    for j=2:ny+1
        for i=2:nx+1
            pc[i,j] = 0.0
        end
    end

    for j = 2:ny+1
        for i = 2:nx+1
            bpc[i,j] = -dy*0.5*((ρa[i+1,j] + ρa[i,j])*uf[i,j] - (ρa[i,j] + ρa[i-1,j])*uf[i-1,j]) - dx*0.5*((ρa[i,j+1] + ρa[i,j])*vf[i,j] - (ρa[i,j] + ρa[i,j-1])*vf[i,j-1]) - pa[i,j]*dx*dy/R/Ta[i,j]/dt + ρold[i,j]*dx*dy/dt
        end
    end
    
    
    # SOR法で圧力の解を求めてる
    for ip=1:ipter
        for j=2:ny+1
            for i=2:nx+1
                ptemp[i,j] = pc[i,j]
            end
        end
        for j=2:ny+1
            for i=2:nx+1
                pc[i,j] = pc[i,j]*(1.0 - alphapa) + alphapa*(aepc[i,j]*pc[i+1,j] + awpc[i,j]*pc[i-1,j] + anpc[i,j]*pc[i,j+1] + aspc[i,j]*pc[i,j-1] + bpc[i,j])/appc[i,j]
            end
        end
        for i=1:nx+2
            pc[i, 1] = pc[i, 2]
            pc[i, ny+2] = pc[i, ny+1]
        end
        for j=1:ny+2
            pc[nx+2, j] = pc[nx+1, j]
            pc[1, j] = pc[2, j]
        end
        rest = convergP(para, pc, ptemp)
        #println("rest PC ", rest, ", ip=", ip)
        if rest < 1e-5 && ip >= 2
            #println("rest PC ", rest, ", ip=", ip)
            break
        end
        
    end
    #println("pc=",pc[2:3,2:3])
    #println("max pc=",findmax(abs.(pc)))
end    

function update(para, ua, va, pc, apua, apva, ρa, pa, Ta, )
    @unpack nx, ny, dx, dy, R = para
    
    for j = 2:ny+1
        for i = 2:nx+1
            pa[i,j] += pc[i,j]
        end
    end
    # u** = u* + u'
    for j=2:ny+1
        for i=2:nx+1
            ua[i,j] = ρa[i,j]*(ua[i,j] -0.5dy*(pc[i+1,j] - pc[i-1,j])/apua[i,j])/(pa[i,j]/R/Ta[i,j]) 
        end
    end
    
    for j=2:ny+1
        for i=2:nx+1
            va[i,j] = ρa[i,j]*(va[i,j] -0.5dx*(pc[i,j+1] - pc[i,j-1])/apva[i,j])/(pa[i,j]/R/Ta[i,j]) 
        end
    end

    for j = 2:ny+1
        for i = 2:nx+1
            ρa[i,j] = pa[i,j]/R/Ta[i,j]
        end
    end


    #println("1st up ua, va, rhoa, pa, Ta=", ua[502:503,2:3], va[2:3,2:3], ρa[502:503,2:3], pa[502:503,2:3], Ta[502:503,2:3])
end

function solvePA2(para, pa, apua, apva, appc, aepc, awpc, anpc, aspc, bpc, ptemp, ρold, Ueq, Veq, ρa, l)
    @unpack nx, ny, dx, dy, dt, ipter, alphapa, = para

    for j=2:ny+1
        for i=2:nx+1
            bpc[i,j] = ρold[i,j]*dx*dy/dt
        end
    end

    for j=3:ny
        for i=3:nx
            bpc[i,j] += -0.5dy*(ρa[i+1,j]Ueq[i+1,j] - ρa[i-1,j]Ueq[i-1,j]) - 0.5dx*(ρa[i,j+1]Veq[i,j+1] - ρa[i,j-1]Veq[i,j-1]) 
        end
    end

    j=2
    for i=2:nx+1
        if i==2
            bpc[i,j] += -0.5dy*(ρa[i+1,j]Ueq[i+1,j] + ρa[i,j]Ueq[i,j]) - 0.5dx*(ρa[i,j+1]Veq[i,j+1] + ρa[i,j]Veq[i,j])  
        elseif i==nx+1
            bpc[i,j] += -0.5dy*(-ρa[i,j]Ueq[i,j] - ρa[i-1,j]Ueq[i-1,j]) - 0.5dx*(ρa[i,j+1]Veq[i,j+1] + ρa[i,j]Veq[i,j]) 
        else
            bpc[i,j] += -0.5dy*(ρa[i+1,j]Ueq[i+1,j] - ρa[i-1,j]Ueq[i-1,j]) - 0.5dx*(ρa[i,j+1]Veq[i,j+1] + ρa[i,j]Veq[i,j]) 
        end
    end
    j=ny+1
    for i=2:nx+1
        if i==2
            bpc[i,j] += -0.5dy*(ρa[i+1,j]Ueq[i+1,j] + ρa[i,j]Ueq[i,j]) - 0.5dx*(-ρa[i,j]Veq[i,j] - ρa[i,j-1]Veq[i,j-1]) 
        elseif i==nx+1
            bpc[i,j] += -0.5dy*(-ρa[i,j]Ueq[i,j] - ρa[i-1,j]Ueq[i-1,j]) - 0.5dx*(-ρa[i,j]Veq[i,j] - ρa[i,j-1]Veq[i,j-1]) 
        else
            bpc[i,j] += -0.5dy*(ρa[i+1,j]Ueq[i+1,j] - ρa[i-1,j]Ueq[i-1,j]) - 0.5dx*(-ρa[i,j]Veq[i,j] - ρa[i,j-1]Veq[i,j-1]) 
        end
    end
    i=2
    #println(0.5dy*(Ueq[i+1,j] + Ueq[i,j]))
    #println(dy*0.5(ρa[i,j] + ρa[i-1,j])*sin(2π*8.5*(l-1)*dt))
    for j=3:ny
        bpc[i,j] += -0.5dy*(ρa[i+1,j]Ueq[i+1,j] + ρa[i,j]Ueq[i,j]) - 0.5dx*(ρa[i,j+1]Veq[i,j+1] - ρa[i,j-1]Veq[i,j-1]) 
    end
    i=nx+1
    for j=3:ny
        bpc[i,j] += -0.5dy*(-ρa[i,j]Ueq[i,j] - ρa[i-1,j]Ueq[i-1,j]) - 0.5dx*(ρa[i,j+1]Veq[i,j+1] - ρa[i,j-1]Veq[i,j-1]) 
    end

    #println("ue, uw, vn, vs ", findmax(abs.(uE)), ",", findmax(abs.(uW)), ", ", findmax(abs.(vN)), ",", findmax(abs.(vS)))
    #println("apua ", findmax(apua))
    #println("appc, aepc ", findmax(appc), findmax(aepc))
    #println("bpc ",findmax(abs.(bpc)))
    
    # SOR法で圧力の解を求めてる
    for ip=1:ipter
        for j=2:ny+1
            for i=2:nx+1
                ptemp[i,j] = pa[i,j]
            end
        end
        for j=2:ny+1
            for i=2:nx+1
                pa[i,j] = pa[i,j]*(1.0 - alphapa) + alphapa*(aepc[i,j]*pa[i+1,j] + awpc[i,j]*pa[i-1,j] + anpc[i,j]*pa[i,j+1] + aspc[i,j]*pa[i,j-1] + bpc[i,j])/appc[i,j]
            end
        end

        rest=convergP(para, pa, ptemp)

        #println("rest PC2 ", rest, ", ip=", ip)
        if rest < 1e-5 && ip>=2
            #println("rest PC2 ", rest, ", ip=", ip)
            break
        end
    end
    #println("pa=",pa[2:3,2:3])
    #println("max pc2=",findmax(abs.(pc2)))
end

function update2(para, ua, va, utemp, vtemp, pa, apua, aeua, awua, anua, asua, apva, aeva, awva, anva, asva, bua, bva, ρa, Ta, )
    @unpack nx, ny, dx, dy, R = para
    
    #println("max ua2-ua=",findmax(abs.(ua2[2:nx,2:ny+1]-ua[2:nx,2:ny+1])))
    #println("max va2-va=",findmax(abs.(va2[2:nx+1,2:ny]-va[2:nx+1,2:ny])))
    
    #println("ua2[2,ny+2], ua[2,ny+2]", ua2[2,ny+2], ",",ua[2,ny+2])

    # u** = u* + u'
    for j=2:ny+1
        for i=2:nx+1
            utemp[i,j] = ρa[i,j]*(aeua[i,j]*ua[i+1,j] + awua[i,j]*ua[i-1,j] + anua[i,j]*ua[i,j+1] + asua[i,j]*ua[i,j-1] -
            0.5dy*(pa[i+1,j] - pa[i-1,j]) + bua[i,j])/apua[i,j]/(pa[i,j]/R/Ta[i,j]) 
        end
    end

    for j=2:ny+1
        for i=2:nx+1
            vtemp[i,j] =  ρa[i,j]*(aeva[i,j]*va[i+1,j] + awva[i,j]*va[i-1,j] + anva[i,j]*va[i,j+1] + asva[i,j]*va[i,j-1] -
            0.5dx*(pa[i,j+1] - pa[i,j-1]) + bva[i,j])/apva[i,j]/(pa[i,j]/R/Ta[i,j]) 
        end
    end

    for j=2:ny+1
        for i=2:nx+1
            ua[i,j] = utemp[i,j]
            va[i,j] = vtemp[i,j]
        end
    end

    for j = 2:ny+1
        for i = 2:nx+1
            ρa[i,j] = pa[i,j]/R/Ta[i,j]
        end
    end
    #println("2nd up ua, va, rhoa, pa, Ta=", ua[502:503,2:3], va[2:3,2:3], ρa[502:503,2:3], pa[502:503,2:3], Ta[502:503,2:3])
end

function convT(para, ua, va, ea, Ta)
    @unpack nx, ny, Cv = para

    for j = 1:ny+2
        for i = 1:nx+2
            Ta[i,j] = (ea[i,j] - 0.5*(ua[i,j]^2 + va[i,j]^2)) / Cv
        end
    end
end

function converg(para, ua, va, utemp, vtemp)
    @unpack nx, ny = para
    
    restua = 0.0
    restva = 0.0
    
    for j=2:ny+1
        for i=2:nx+1
            restua = restua + abs(ua[i,j] - utemp[i,j])
            
        end
    end
    for j=2:ny+1
        for i=2:nx+1
            restva = restva + abs(va[i,j] - vtemp[i,j])
        end
    end
    #rest = max(restua, restva)
    rest = (restua + restva)/2.0
    return rest/(nx*ny)
end

function convergP(para, pc, ptemp)
    @unpack nx, ny = para
    
    rest = 0.0
    
    for j=2:ny+1
        for i=2:nx+1
            rest = rest + abs(pc[i,j] - ptemp[i,j])
        end
    end

    
    return rest/(nx*ny)
end

function output(para, ua, va, pa, Ta, ea, ρa,  x, y, l)
    @unpack outfile, dt, nx, ny = para
    #println("time=", time)
    filename = outfile*string(l)*".h5"
    h5open(filename, "w") do file
        write(file, "U", ua[2:nx+1, 2:ny+1])
        write(file, "V", va[2:nx+1, 2:ny+1])
        write(file, "P", pa[2:nx+1, 2:ny+1])
        write(file, "T", Ta[2:nx+1, 2:ny+1])
        write(file, "E", ea[2:nx+1, 2:ny+1])
        write(file, "rho", ρa[2:nx+1, 2:ny+1])
        write(file, "X", x[2:nx+1, 2:ny+1])
        write(file, "Y", y[2:nx+1, 2:ny+1])
        write(file, "ua_latest", ua)
        write(file, "va_latest", va)
        write(file, "pa_latest", pa)
        write(file, "ea_latest", ea)
        write(file, "ρa_latest", ρa)
        write(file, "l", l)
    end
    
    #open( outfile*time*".csv", "w")　 do f
    #    println(f, "x, y, u, v, pa")
    #    for j=2:ny+1
    #        for i=2:nx+1
    #            println(f, x[i,j], ",", y[i,j], ",",(u[i,j]+u[i-1,j])/2., ",", (v[i,j]+v[i,j-1])/2.0, ", ", pa[i,j] )
    #        end
    #    end
    #end
end

function main(para, restartFile)
    @unpack nx, ny, nstep, niter, ndata, dx, dy = para

    μ = zeros(Float64, nx+2, ny+2)
    k = zeros(Float64, nx+2, ny+2)

    # auxiliary velocity
    ua = zeros(Float64, nx+2, ny+2)
    utemp = zeros(Float64, nx+2, ny+2)
    utemp2 = zeros(Float64, nx+2, ny+2)
    uold = zeros(Float64, nx+2, ny+2)
    apua = zeros(Float64, nx+2, ny+2)
    aeua = zeros(Float64, nx+2, ny+2)
    awua  = zeros(Float64, nx+2, ny+2)
    anua = zeros(Float64, nx+2, ny+2)
    asua = zeros(Float64, nx+2, ny+2)
    bua  = zeros(Float64, nx+2, ny+2)
    buae = zeros(Float64, nx+2, ny+2)
    buaw = zeros(Float64, nx+2, ny+2) 
    buan = zeros(Float64, nx+2, ny+2) 
    buas = zeros(Float64, nx+2, ny+2)
    
    va = zeros(Float64, nx+2, ny+2)
    vtemp = zeros(Float64, nx+2, ny+2)
    vtemp2 = zeros(Float64, nx+2, ny+2)
    vold = zeros(Float64, nx+2, ny+2)
    apva = zeros(Float64, nx+2, ny+2)
    aeva = zeros(Float64, nx+2, ny+2)
    awva  = zeros(Float64, nx+2, ny+2)
    anva = zeros(Float64, nx+2, ny+2)
    asva = zeros(Float64, nx+2, ny+2)
    bva = zeros(Float64, nx+2, ny+2)
    bvae = zeros(Float64, nx+2, ny+2)
    bvaw = zeros(Float64, nx+2, ny+2)
    bvan = zeros(Float64, nx+2, ny+2)
    bvas = zeros(Float64, nx+2, ny+2)
    
    Ueq = zeros(Float64, nx+2, ny+2)
    Veq = zeros(Float64, nx+2, ny+2)

    uf = zeros(Float64, nx+2, ny+2)
    vf = zeros(Float64, nx+2, ny+2)

    ea = zeros(Float64, nx+2, ny+2)
    eold = zeros(Float64, nx+2, ny+2)
    etemp = zeros(Float64, nx+2, ny+2)
    apea = zeros(Float64, nx+2, ny+2)
    aeea = zeros(Float64, nx+2, ny+2)
    awea  = zeros(Float64, nx+2, ny+2)
    anea = zeros(Float64, nx+2, ny+2)
    asea = zeros(Float64, nx+2, ny+2)
    bea = zeros(Float64, nx+2, ny+2)

    ρa = zeros(Float64, nx+2, ny+2)
    ρold = zeros(Float64, nx+2, ny+2)
    Ta = zeros(Float64, nx+2, ny+2)
    Told = zeros(Float64, nx+2, ny+2)

    phix = zeros(Float64, nx+2, ny+2)
    phiy = zeros(Float64, nx+2, ny+2)

    # pressure, corrected pressure
    pa = zeros(Float64, nx+2, ny+2)
    pc = zeros(Float64, nx+2, ny+2)
    ptemp = zeros(Float64, nx+2, ny+2)
    pold = zeros(Float64, nx+2, ny+2)
    appc = zeros(Float64, nx+2, ny+2)
    aepc = zeros(Float64, nx+2, ny+2)
    awpc = zeros(Float64, nx+2, ny+2)
    anpc = zeros(Float64, nx+2, ny+2)
    aspc = zeros(Float64, nx+2, ny+2)
    bpc = zeros(Float64, nx+2, ny+2)

    # coordinate
    x = zeros(Float64, nx+1, ny+1)
    y = zeros(Float64, nx+1, ny+1)
    for j = 2:ny+1
        for i = 2:nx+1
            x[i,j] = (i-2 + 0.5)*dx
        end
    end
    for j = 2:ny+1
        for i = 2:nx+1
            y[i,j] = (j-2 + 0.5)*dy
        end
    end

    rest = 0.0
    
    nstart = restart(para, restartFile, ua, va, pa, ea, ρa )
    convT(para, ua, va, ea, Ta)
    println("nstart=", nstart)


    for l = nstart + 1:nstart + nstep
        #println(l, "===========================================")

        for j = 1:ny+2
            for i = 1:nx+2
                uold[i,j] = ua[i,j]
                vold[i,j] = va[i,j]
                eold[i,j] = ea[i,j]
                ρold[i,j] = ρa[i,j]
                Told[i,j] = Ta[i,j]
                pold[i,j] = pa[i,j]
            end
        end

        BC(para, ua, va, pa, pc, ea, ρa, l)
        updateTransport(para, k, μ)
        
        for n = 1:niter
            BC(para, ua, va, pa, pc, ea, ρa, l)

            makeae(para, apea, aeea, awea, anea, asea, ρa, ua, va, k, ρold)
            solveEA(para, ua, va, pa, ea, eold, ρa, ρold, etemp, bea, buae, buaw, buan, buas, bvae, bvaw, bvan, bvas, apea, aeea, awea, anea, asea, k, μ)

            makeau(para, apua, aeua, awua, anua, asua, ρa, ua, va, μ, ρold)
            makeav(para, apva, aeva, awva, anva, asva, ρa, ua, va, μ, ρold)
            makebuv(para, ua, va, pa, bua, buae, buaw, buan, buas, uold, bva, bvae, bvaw, bvan, bvas, vold, ρold, μ)
            makeUVeq(para, ua, va, pa, apua, aeua, awua, anua, asua, bua,  apva, aeva, awva, anva, asva, bva, Ueq, Veq)
            solveUA(para, Ueq, ua, pa, utemp, apua)            
            solveVA(para, Veq, va, pa, vtemp, apva)
            RhieChow(para, ua, va, pa, uf, vf, apua, apva, l)

            convT(para, ua, va, ea, Ta)
            makeap(para, appc, aepc, awpc, anpc, aspc, apua, apva, ρa, Ta)
            solvePC(para, ρa, pc, ptemp, appc, aepc, awpc, anpc, aspc, bpc, uf, vf, Ta, ρold, pa)
            update(para, ua, va, pc, apua, apva, ρa, pa, Ta, )

            makeae(para, apea, aeea, awea, anea, asea, ρa, ua, va, k, ρold)
            solveEA(para, ua, va, pa, ea, eold, ρa, ρold, etemp, bea, buae, buaw, buan, buas, bvae, bvaw, bvan, bvas, apea, aeea, awea, anea, asea, k, μ)

            makeau(para, apua, aeua, awua, anua, asua, ρa, ua, va, μ, ρold)
            makeav(para, apva, aeva, awva, anva, asva, ρa, ua, va, μ, ρold)
            makebuv(para, ua, va, pa, bua, buae, buaw, buan, buas, uold, bva, bvae, bvaw, bvan, bvas, vold, ρold, μ)
            convT(para, ua, va, ea, Ta)
            makeap(para, appc, aepc, awpc, anpc, aspc, apua, apva, ρa, Ta)
            makeUVeq(para, ua, va, pa, apua, aeua, awua, anua, asua, bua,  apva, aeva, awva, anva, asva, bva, Ueq, Veq)
            solvePA2(para, pa, apua, apva, appc, aepc, awpc, anpc, aspc, bpc, ptemp, ρold, Ueq, Veq, ρa, l)
            update2(para, ua, va, utemp, vtemp, pa, apua, aeua, awua, anua, asua, apva, aeva, awva, anva, asva, bua, bva, ρa, Ta)
            makeae(para, apea, aeea, awea, anea, asea, ρa, ua, va, k, ρold)
            solveEA(para, ua, va, pa, ea, eold, ρa, ρold, etemp, bea, buae, buaw, buan, buas, bvae, bvaw, bvan, bvas, apea, aeea, awea, anea, asea, k, μ)

            rest = converg(para, ua, va, utemp2, vtemp2)
            #println("iter, redisual ", n,", ",rest)
            for j=2:ny+1
                for i=2:nx+1
                    utemp2[i,j] = ua[i,j]
                    vtemp2[i,j] = va[i,j]
                end
            end
            
            if rest <=1e-7 && n>=2
                if l%10 == 0
                    println(l, " th step, ", n, " th iteration, residual=", rest)
                end
                break
            end

            BC(para, ua, va, pa, pc, ea, ρa, l)
        end
        
        if l%ndata == 0
            println("residual=", rest)
            convT(para, ua, va, ea, Ta)
            output(para, ua, va, pa, Ta, ea, ρa,  x, y, l)
        end
    end
    #streamfunction(para, u, v, psi)
    
end

function restart(para, restartFile, ua, va, pa, ea, ρa )
    @unpack nx, ny,  = para
    
    file = h5open(restartFile, "r") 
    ua .= read(file, "ua_latest")
    va .= read(file, "va_latest")
    pa .= read(file, "pa_latest")
    ea .= read(file, "ea_latest")
    ρa .= read(file, "ρa_latest")
    nstart = read(file, "l")
    
    return nstart
end

#para = Param(nu=0.02, nstep=500, outfile="Data-simple/uv", nx=30, ny=30, ndata=50, alphapa=0.5)

#@time main(para)

