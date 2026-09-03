using ERT3D
using BenchmarkTools

function profile_scheme(scheme)

    for N in (32, 64, 128)

        println()
        println("N = $N")

        grid = Grid(N)
        D = Central4(grid)

        X = reshape(grid.x, :, 1, 1)
        Y = reshape(grid.y, 1, :, 1)
        Z = reshape(grid.z, 1, 1, :)

        rho = 1.0 .+ 0.1 .* sin.(X .+ 0.0 .* Y .+ 0.0 .* Z)
        u   = 1.0 .+ 0.2 .* cos.(X .+ 0.0 .* Y .+ 0.0 .* Z)
        v   = 1.0 .+ 0.2 .* sin.(0.0 .* X .+ Y .+ 0.0 .* Z)
        w   = 1.0 .+ 0.2 .* cos.(0.0 .* X .+ 0.0 .* Y .+ Z)

        phi = 1.0 .+ 0.1 .* sin.(X .+ Y .+ Z)

        out = similar(rho)

        @btime $scheme(
            $out,
            $rho,
            $u,
            $v,
            $w,
            $phi,
            $D,
            $grid,
        )
    end
end

profile_scheme(ERT3D.KennedyGruber!)
profile_scheme(ERT3D.feiereisen!)