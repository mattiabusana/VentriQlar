function find_vq_min(fio2, hb, be, po2, pco2)



    println("")
    println("Finding mimimum VA/Q compartment...")

    for r = 0.01:0.005:2
        vaq, _, _, _, _ = fndvaq_modified(fio2, hb, be, po2, pco2, r)
        print("-")

        if vaq >= 0
            vaq = ceil(vaq, digits = 2)
            r = floor(r, digits = 2)
            println("")
            println("Found!")
            println("Minimum R = ", r)
            println("Minimum VAQ compartment = ", vaq)
            println("")
            return vaq, r
            break
        end

    end

end


function find_vq_max(fio2, hb, be, po2, pco2)

    println("")
    println("Finding maximum VA/Q compartment...")

    global last_r = nothing
    global last_vaq = nothing

    for r = 1:0.1:20
        global vaq, _, co2_compart, _, _, pco2_compart =
            fndvaq_modified(fio2, hb, be, po2, pco2, r)
        print("-")

        #println("VAQ = ", vaq)          #For debug
        #println("PCO2 = ", pco2_compart)
        #println("CO2 content = ", co2_compart)

        if co2_compart <= 0 || vaq >= 105
            println("")
            println("Found!")
            println("Maximum R = ", last_r)
            println("Maximum VAQ compartment = ", last_vaq)
            println("")
            break
        else
            last_r = ceil(r, digits = 2)
            last_vaq = floor(vaq, digits = 2)
        end

    end

    return last_vaq, last_r

end
