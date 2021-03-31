function gvnvaq(
    given_vaq,
    fio2,
    hb,
    be,
    po2_mix,
    pco2_mix,
    last_r;
    minimum_r = 0.01,
    maximum_r = 100,
)

    if last_r >= 5
        last_r = last_r - 3
    end
    #include("fndvaq.jl")
    #include("physiology_functions.jl")

    array_r = log_distribution(minimum_r, maximum_r, 1000)

    for r in collect(Iterators.rest(array_r, last_r))
        global vaq_calculated, cont_o2, cont_co2, ph, _, _ =
            fndvaq_modified(fio2, hb, be, po2_mix, pco2_mix, r)
        gradient_vaq = abs((given_vaq - vaq_calculated) / given_vaq)

        #println("Gradient VAQ = ", gradient_vaq)    # Debug only

        if given_vaq <= 0.1 && given_vaq >= 0
            tolerance = 0.25
        elseif given_vaq > 0.1 && given_vaq <= 1
            tolerance = 0.1
        elseif given_vaq > 1 && given_vaq <= 10
            tolerance = 0.1
        else
            tolerance = 0.05
        end

        print(".")

        if gradient_vaq < tolerance
            global last_index = findfirst(x -> x == r, array_r)
            break
        else
            continue
        end

    end

    return vaq_calculated, cont_o2, cont_co2, ph, last_index

end
