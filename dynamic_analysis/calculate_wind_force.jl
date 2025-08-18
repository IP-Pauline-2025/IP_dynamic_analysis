function constant_sign(phi::Float64)::Float64
    Xi_U2 = 2 / phi^2 * (phi - 1 + exp(-phi))
    return Xi_U2
end

function calculate_wind_force(category::Int, v_b0::Float64, lat::Float64, z::Float64, geometry::Int, A::Float64, A_encl::Float64, A_f::Float64, A_c::Float64, A_csup::Float64)
    # Calculate the wind force on elongated structures (lattice towers) based on Euro Code
    
    #= 
    input values:
        category = terrain category (0-4)
        v_b0 = basic wind velocity in m/s
        lat = latitude of the site in degrees
        z = height in m
        geometry = (4 = square, 3 = triangular)
        𝐴  is the sum of the projected area of the members and gusset plates of the face projected normal to the face
        A_encl is the area enclosed by the boundaries of the face projected normal to the face
        𝐴_f is the total projected area when viewed normal to the face of the flat-sided section members in the face
        𝐴_c is the total projected area when viewed normal to the face of the circular-section members in the face in sub critical regimes
        𝐴_csup is the total projected area when viewed normal to the face, of the circular-section members in the face in supercritical regimes
    
    
        output value:
        F_wp = wind force in kN
    =#

    # Define a dictionary with parameters for each terrain category according to Table B.3 and B.4
    terrain_parameters = Dict(
    0 => (z_0 = 0.003, z_min = 1, c_x = 1.28, alpha_1 = 2.10, alpha_2 = 70),
    1 => (z_0 = 0.01, z_min = 1, c_x = 1.17, alpha_1 = 2.15, alpha_2 = 55),
    2 => (z_0 = 0.05, z_min = 2, c_x = 1.00, alpha_1 = 2.30, alpha_2 = 35),
    3 => (z_0 = 0.3, z_min = 5, c_x = 0.75, alpha_1 = 2.45, alpha_2 = 25),
    4 => (z_0 = 1.0, z_min = 10, c_x = 0.55, alpha_1 = 2.50, alpha_2 = 18),
    )
    # Get the parameters for the selected terrain category
    params = terrain_parameters[category]
    z_0 = params.z_0
    z_min = params.z_min
    c_x = params.c_x
    alpha_1 = params.alpha_1
    alpha_2 = params.alpha_2

    # Calculate the basic wind velocity 
        # v_b0 = basic wind velocity in m/s 
        # irrespective of wind direction and time of year at hight of 10 m above ground level in flat open county
        c_dir = 1.0 # directional factor 1,0 unless the National Annex gives a different value
        c_season = 1.0 # seasonal factor 1,0 unless the National Annex gives a different value
        c_alt = 1.0 # 1,0 altitude factor 1,0 unless the National Annex gives a different value
        c_prop = 1.0 # probability factor is 1,0 because the data already contains the probability
    v_b = c_prop * c_dir * c_season * c_alt * v_b0



    #calculate mean wind velocity according to 6.4
        # orography factor is 1,0 unless otherwise specified in 6.3.3
            c_0 = 1.0 
            
        # Calculate roughness factor according to B.15
            # Calculate the gradient height based on the given parameters
                #Coriolis Parameter depented of latidude of the site (only applies for storm winds)
                f = 1.454 * 10^(-4) * sin(lat * pi / 180)
                if f == 0
                    error("Coriolis parameter f is zero.")
                end
                #tarrain factor k_r according to 6.7
                z_0II = 0.05 # roughness length for terrain category II in m
                k_r = 0.19 * (z_0 / z_0II)^(0.07)
                z_g = k_r * abs(v_b) / (15 * f)
                if z_g == 0 || z_0 <= 0
                    error("Invalid gradient height or roughness length.")
                end
            a = z / z_g
            z_max = 300 # maximum height for the roughness factor is to be taken as 300m
            if z >= z_min && z <= z_max
            c_r = c_x * (log(z / z_0) + 5.75 *a) / log(10/z_0)
            elseif z <= z_min
            c_r = c_x * (log(z_min / z_0) + 5.75 *a) / log(10/z_0)
            else 
                error("z >= z_max")
            end
    v_m = c_r * c_0 * v_b

    # Calculate the peak wind velocity according to 6.12
        # Calculate the turbulence intensity I_u according to B.17
        if z >= z_min && z <= z_max
        I_u = (alpha_1 * ( 1- z/z_g) * (0.538 + 0.090 *log(z/z_0))^((1-z/z_g)^16)) / (c_0 * (log(z/z_0) + alpha_2 * (z/z_g)) * (1 + 0.156 *  log((6 * z_g)/z_0)))
        elseif z <= z_min
        I_u = (alpha_1 * ( 1- z_min/z_g) * (0.538 + 0.090 *log(z_min/z_0))^((1-z_min/z_g)^16))/ (c_0 * (log(z_min/z_0) + alpha_2 * (z_min/z_g)) * (1 + 0.156 * log((6 * z_g)/z_0)))
        else 
            error("z >= z_max")
        end
        k_u = 2.8 #peak factor for turbulence is 2,8 unless the National Annex gives different values
    v_p = v_m * (1 + k_u * I_u)  # Calculate the peak wind velocity m/s
   
    # Calculate the peak velocity pressure according to 6.11
        roh = 1.25 #air density is 1,25 kg/m^3 unless the National Annex gives different values
    q_p = 0.5 * roh * v_p^2


    # Calculate the force coefficient c_f according to Annex E.1
        # overall longitudinal force coefficient of a section without end-effects (E.22)
            # Define a dictionary with parameters for each structure form
            geometry_parameter = Dict(
            4 => (C_1 = 2.25, C_2 = 1.5),
            3 => (C_1 = 1.90, C_2 = 1.4)
            )
            # Get the parameters for the selected structur form
            param = geometry_parameter[geometry]
            C_1 = param.C_1
            C_2 = param.C_2

            # phi is the solidity ratio
            phi = A / A_encl #(E.21) 

            # force coefficients for sections composed of flat-sided, sub critial circular and supercritical circula-section members
            c_f0f = 1.76 * C_1 * (1 - C_1 * phi + phi^2) #(E.23)
            c_f0c = C_1 *(1 - C_2 *phi) + (C_1 + 0.875) * phi^2 #(E.24)
            c_f0csup = 1.9 - sqrt((1 - phi) * (2.8 - 1.14 * C_1 + phi)) #(E.25)

            #The reference area shall be taken as 𝐴_ref = 𝐴_f + 𝐴_c + 𝐴_csup
            A_ref = A_f + A_c + A_csup
        c_f0 = c_f0f * A_f / A_ref + c_f0c * A_c / A_ref + c_f0csup * A_csup / A_ref
    
        # psi_r reduction factor for square sections with rounded corners
        psi_r = 1.0

        # psi_lamda end-effect factor for elements with free-end flow
        psi_lamda = 1.0

    c_f = c_f0 * psi_r * psi_lamda
    

    # calculate wind force F_wp in kN
    F_wp = c_f * q_p * A_ref * 0.001 #0,001 is the conversion factor from N to kN
    return F_wp
end
