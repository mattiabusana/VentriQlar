function fndvaq_modified(fio2, hb, be, po2_mix, pco2_mix, r_gas)

    #include("physiology_functions.jl")

    global sat_mix = severinghaus(po2_mix) * 100

    for i = 1:3
        global ph_mix = ph_from_be_zander(be, pco2_mix, hb, sat_mix)
        global sat_mix = zander_saturation(po2_mix, hb, pco2_mix, ph_mix)
        global o2_cont_mix = oxygen_content(sat_mix, hb, po2_mix)
        global co2_cont_mix = douglas_to_content(ph_mix, pco2_mix, hb, sat_mix)

    end

    global pco2_compart = 2
    r_blood = 100000000000000
    gradient_r = 1000000000000
    global index = 0

    while gradient_r >= 0.05

        index += 1

        global po2_compart =
            fio2 * 713 - (pco2_compart * (1 - fio2 * (1 - r_gas))) / r_gas
        global sat_compart = severinghaus(po2_compart) * 100
        global ph_compart = ph_from_be_zander(be, pco2_compart, hb, sat_compart)

        if po2_compart < 0
            global o2_cont_compart = 0
        else
            global sat_compart =
                zander_saturation(po2_compart, hb, pco2_compart, ph_compart)
            global o2_cont_compart =
                oxygen_content(sat_compart, hb, po2_compart)
        end

        global co2_cont_compart =
            douglas_to_content(ph_compart, pco2_compart, hb, sat_compart)


        if co2_cont_compart < 0
            co2_cont_compart = 0
        end

        #println("PO2 = ", po2_compart)   # Debug only
        #println("PCO2 = ", pco2_compart)


        diff_alv_ven_o2 = o2_cont_compart - o2_cont_mix
        diff_ven_alv_co2 = co2_cont_mix - co2_cont_compart

        r_blood = diff_ven_alv_co2 / diff_alv_ven_o2

        gradient_r = abs(r_gas - r_blood) / r_gas

        #println("R = ", r_blood)
        #println("Gradient = ", gradient_r)
        #println("")

        if po2_compart <= 0
            vaq = -500
            break
        end

        if gradient_r >= 1
            pco2_compart += 1
        elseif gradient_r < 1 && gradient_r >= 0.5
            pco2_compart += 0.1
        elseif gradient_r < 0.5 && gradient_r >= 0.05
            pco2_compart += 0.1
        else
            #Debug only
            #println("Gradient = ", gradient_r)
            #println("Found R blood = R gas")
            #println("PaCO2 = ", pco2_compart)
            #println("PaO2 = ", po2_compart)
            #println("R blood = ", r_blood)
            #println("Diff alv ven O2 = ", diff_alv_ven_o2)
            #println("Diff alv ven CO2 = ", diff_ven_alv_co2)
            #println("Cont CO2 compart = ", co2_cont_compart)
            #println("Cont O2 compart = ", o2_cont_compart)
            #println("pH compart = ", ph_compart)
            #println("Index = ", index)

            global vaq = (8.63 * diff_ven_alv_co2) / pco2_compart

            #println("VAQ = ", vaq)
            break
        end
    end

    return vaq,
    o2_cont_compart,
    co2_cont_compart,
    ph_compart,
    po2_compart,
    pco2_compart
end
