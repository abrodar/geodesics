using DifferentialEquations
using LinearAlgebra
using GeometryBasics
using GLMakie
using Makie
using ColorTypes: RGBAf0

# Euler's method
function euler(f, Y0, tspan, h, p)
    # Obtain the range of steps from
    # the tuple tspan.
    t0, tk = tspan
    # Create an array of fixed steps 
    # of size h stored in t.
    t = t0:h:tk
    # Y is an array containing the points
    # on a surface forming a geodesic.
    Y = [Y0]
    # Creates a new uninitialized array that 
    # has the same shape and type as Y0.
    # The function f (e.g., geodesic_sphere!) 
    # writes the computed derivatives into du,
    # modifying it in-place.
    du = similar(Y0) 
    for ti in t[1:(end-1)]
        # Perform a numerical step of 
        # the Euler's method.
        f(du, Y[end], p, ti)
        # Add the new result to Y
        # using the computed du
        # from f.
        push!(Y, Y[end] + h * du)
    end
    # Return the steps array and a 
    # geodesic value for each step.
    return t, Y
end

#Runge-Kutta 4 method
function rk4(f, Y0, tspan, h, p)
    # Obtain the range of steps from
    # the tuple tspan.
    t0, tk = tspan
    # Create an array of fixed steps 
    # of size h stored in t.
    t = t0:h:tk
    # Y is an array containing the points
    # on a surface forming a geodesic.
    Y = [Y0]
    # Creates a new uninitialized array that 
    # has the same shape and type as Y0.
    du = similar(Y0) 
    for ti in t[1:(end-1)]
        k1 = f(du, Y[end], p, ti); k1 = h * du
        k2 = f(du, Y[end] + k1/2, p, ti + h/2); k2 = h * du
        k3 = f(du, Y[end] + k2/2, p, ti + h/2); k3 = h * du
        k4 = f(du, Y[end] + k3, p, ti + h); k4 = h * du
        # Add the next approximation to Y
        push!(Y, Y[end] + (k1 + 2*k2 + 2*k3 + k4)/6)
    end    
    # Return the steps array and a 
    # geodesic value for each step.
    return t, Y
end

# Geodesic 1st order ODE system for a sphere
function geodesic_sphere!(du, u, p, t)
    y1, y2, y3, y4 = u
    du[1] = y2
    du[2] = cos(y1) * sin(y1) * y4^2
    du[3] = y4
    du[4] = -2 * (cos(y1) / sin(y1)) * y2 * y4
    return du
end

# Geodesic 1st order ODE system for a torus
function geodesic_torus!(du, u, p, t)
    R, r = p  # Unpack parameters
    y1, y2, y3, y4 = u
    du[1] = y2
    du[2] = -((R/r) * sin(y1) + sin(y1) * cos(y1)) * y4^2
    du[3] = y4
    du[4] = ((2 * r * sin(y1)) / (R + r * cos(y1))) * y2 * y4
    return du
end

# Convert sphere into 3D cartesian coordinates
function sphere_to_cartesian(u, v)
    x = sin(u) * cos(v)
    y = sin(u) * sin(v)
    z = cos(u)
    return Point3f0(x, y, z)  # Use Point3f0 for GLMakie
end

# Ploting functionality using GLMakie
function plot_geodesics_on_sphere(list_of_us_vs)
    # Sphere mesh setup
    θ = range(0, π, length=50)
    φ = range(0, 2π, length=50)
    points = [sphere_to_cartesian(θi, φj) for θi in θ, φj in φ]
    faces = TriangleFace{Int}[]
    nθ, nφ = length(θ), length(φ)
    for i in 1:nθ-1, j in 1:nφ-1
        a = (i-1)*nφ + j
        b = a + 1
        c = a + nφ
        d = c + 1
        push!(faces, TriangleFace(a, b, c))
        push!(faces, TriangleFace(b, d, c))
    end
    sphere_mesh = GeometryBasics.Mesh(vec(points), faces)
    fig = Figure()
    ax = LScene(fig[1, 1], show_axis=true)
    mesh!(ax, sphere_mesh, color=RGBAf(0.3, 0.6, 1.0, 0.3), shading=Makie.automatic,  transparency=true, alpha=0.8)
    colors = [:red, :orange, :limegreen] 
    labels = ["Euler", "RK4", "DOPRI5"]
    for (i, (us, vs)) in enumerate(list_of_us_vs)
        geo_points = [sphere_to_cartesian(u, v) for (u, v) in zip(us, vs)]
        lines!(ax, geo_points, color=colors[i], linewidth=7, label=labels[i])
    end
    axislegend(ax, position=:rb)
    display(fig)
    #save("geodesics_on_sphere.png", fig)
    return fig
end

# Convert torus into 3D cartesian coordinates
function torus_to_cartesian(u, v, R, r)
    x = (R + r * cos(u)) * cos(v)
    y = (R + r * cos(u)) * sin(v)
    z = r * sin(u)
    return Point3f0(x, y, z)
end

# Ploting functionality using GLMakie
function plot_geodesics_on_torus(list_of_us_vs, p)
    R, r = p
    u_range = range(0, 2π, length=60)
    v_range = range(0, 2π, length=60)
    points = [torus_to_cartesian(ui, vi, R, r) for ui in u_range, vi in v_range]
    faces = TriangleFace{Int}[]
    nu, nv = length(u_range), length(v_range)
    for i in 1:nu-1, j in 1:nv-1
        a = (i-1)*nv + j
        b = a + 1
        c = a + nv
        d = c + 1
        push!(faces, TriangleFace(a, b, c))
        push!(faces, TriangleFace(b, d, c))
    end
    torus_mesh = GeometryBasics.Mesh(vec(points), faces)
    fig = Figure()
    ax = LScene(fig[1, 1], show_axis=true)
    mesh!(ax, torus_mesh, color=RGBAf(0.3, 0.6, 1.0, 0.3), shading=Makie.automatic,  transparency=true, alpha=0.8)
    colors = [:red, :orange, :limegreen] 
    labels = ["Euler", "RK4", "DOPRI5"]
    for (i, (us, vs)) in enumerate(list_of_us_vs)
        geo_points = [torus_to_cartesian(u, v, R, r) for (u, v) in zip(us, vs)]
        lines!(ax, geo_points, color=colors[i], linewidth=7, label=labels[i])
    end
    axislegend(ax, position=:rb)
    display(fig)
    #save("geodesics_on_torus.png", fig)
    return fig
end

# Main function
function main()
    # [u, du/dt, v, dv/dt], 
    # (u, v) is the initial point on the sphere,
    # while (d1, d2) are the derivatives
    Y0 = [pi/2, 0.5, 0.0, 0.5]  
    # step interval
    tspan = (0.0, 15.0)
    # Step size
    h = 0.0001
    # No extra parameters for a sphere
    p = -1

    # 1. Sphere
    # Euler's method
    # List of computed geodesics 
    geodesics = []
    # Solve with Euler
    t, Y = euler(geodesic_sphere!, Y0, tspan, h, p)
    # Convert Y to a matrix for 
    # easier handling.
    # With "..." expand a collection 
    # of vectors into separate arguments.
    # hcat horizontally concatenates these vectors 
    # into 4xN then transpose them into Nx4
    # now each row is the state 
    # at a time step, more natural for plotting
    Y_matrix = hcat(Y...)'
    # since y1 = u, and y3 = v 
    # extract these from the matrix
    # for plotting
    push!(geodesics, (Y_matrix[:, 1], Y_matrix[:, 3]))

    #RK4 
    t, Y = rk4(geodesic_sphere!, Y0, tspan, h, p)
    Y_matrix = hcat(Y...)'
    push!(geodesics, (Y_matrix[:, 1], Y_matrix[:, 3]))

    #DOPRI5
    # Define the ODEProblem
    p = 0.0
    prob = ODEProblem(geodesic_sphere!, Y0, tspan, p)
    # Solve with DOPRI5
    t_eval = 0:0.0001:15.0
    sol = solve(prob, DP5(), abstol=1e-9, reltol=1e-9, saveat=t_eval)
    # Extract us, vs from sol and plot
    t = sol.t
    Y_matrix = hcat(sol.u...)'
    push!(geodesics, (Y_matrix[:, 1], Y_matrix[:, 3]))
    # Plot all geodesics on the same sphere
    plot_geodesics_on_sphere(geodesics)

    # 2. Torus
    # Initial values
    Y0 = [0.0, 0.485, 0.0, 0.1] 
    tspan = (0.0, 76)
    h = 0.0001
    # Big radius of a torus
    R = 3
    # Small radius of a torus
    r = 1
    # Make a tuple p, for parameters
    p = (R, r)
    # List of computed geodesics 
    geodesics = []

    # Euler's method
    t, Y = euler(geodesic_torus!, Y0, tspan, h, p)
    Y_matrix = hcat(Y...)'
    push!(geodesics, (Y_matrix[:, 1], Y_matrix[:, 3]))

    #RK4 
    t, Y = rk4(geodesic_torus!, Y0, tspan, h, p)
    Y_matrix = hcat(Y...)'
    push!(geodesics, (Y_matrix[:, 1], Y_matrix[:, 3]))

    #DOPRI5
    # Define the ODE problem
    prob = ODEProblem(geodesic_torus!, Y0, tspan, p)
    # Set evaluation points for solution saving
    t_eval = 0:0.0001:76
    # Solve using Dormand-Prince 5 method with strict tolerances
    sol = solve(prob, DP5(), abstol=1e-12, reltol=1e-12, saveat=t_eval)	
    # Extract time steps and geodesic values
    t = sol.t
    Y_matrix = hcat(sol.u...)'
    push!(geodesics, (Y_matrix[:, 1], Y_matrix[:, 3]))

    # Plot all geodesics on the same torus
    plot_geodesics_on_torus(geodesics, p)
end

# Run
main()


