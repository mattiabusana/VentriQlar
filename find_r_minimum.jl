function find_r_minimum(fio2, hb, be, po2, pco2)



    println("")
    println("Finding mimimum VA/Q compartment...")

    for r = 0.01:0.005:2
        vaq, _, _, _, _ = fndvaq_modified(fio2, hb, be, po2, pco2, r)
        print("-")

        if vaq >= 0
            vaq = ceil(vaq, digits = 2)
            r = floor(r, digits = 1)
            println("")
            println("Found!")
            println("")
            return vaq, r
            break
        end

    end

end
