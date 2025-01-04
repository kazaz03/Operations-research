using LinearAlgebra; 
#=Radili skupa 
Amina Kazazović 19364
Daris Mujkić 19413=#

# Algoritam za topološko sortiranje 
function cpm(A, P, T)
    pocetnoVrijeme = zeros(Float64, length(A))
    
    # najraniji pocetak za sve ove zadatke
    for i in 1:length(A)
        prethodnik = P[i]
        if prethodnik == "-"
            pocetnoVrijeme[i] = 0
        else
            prethodnici = split(prethodnik, ",")
            maxVrijeme = 0.0
            for pred in prethodnici
                indeksPrethodnika = findfirst(x -> x == pred, A)
                if indeksPrethodnika !== nothing
                    maxVrijeme = max(maxVrijeme, pocetnoVrijeme[indeksPrethodnika] + T[indeksPrethodnika])
                end
            end
            pocetnoVrijeme[i] = maxVrijeme
        end
    end

    # najkasniji zavrsetak
    maksimalnoVrijeme = 0.0
    posljednjiZadatak = ""
    for i in eachindex(A)
        zavrsetak = pocetnoVrijeme[i] + T[i]
        if zavrsetak > maksimalnoVrijeme
            maksimalnoVrijeme = zavrsetak
            posljednjiZadatak = A[i]
        end
    end

    ukupnoVrijeme = maksimalnoVrijeme - T[findfirst(x -> x == posljednjiZadatak, A)]
    kriticniPut = [posljednjiZadatak]
    trenutniZadatak = posljednjiZadatak
    
    while true
        i = findfirst(x -> x == trenutniZadatak, A)
        prethodnik = P[i]
        if prethodnik == "-"
            break
        else
            prethodnici = split(prethodnik, ",")
            for pred in prethodnici
                indeksPrethodnika = findfirst(x -> x == pred, A)
                if pocetnoVrijeme[indeksPrethodnika] + T[indeksPrethodnika] == ukupnoVrijeme
                    ukupnoVrijeme -= T[indeksPrethodnika]
                    trenutniZadatak = pred
                    push!(kriticniPut, trenutniZadatak)
                    break
                end
            end
        end
    end
    
    return maksimalnoVrijeme, reverse(kriticniPut)
end

#=Primjer iz laba  
Z=12
put= C -D - G - I=#
A = ["A", "B", "C", "D", "E", "F", "G", "H", "I"]
P = ["-", "-", "-","C","A","A","B,D","E","F,G"]
T = [3, 3, 2, 2, 4, 1, 4, 1, 4]
Z,put = cpm(A, P, T)
println("Z = ",Z)
print("put = ")
for i in eachindex(put)
    if i!=length(put)
        print(put[i])
        print(" - ")
    else 
        println(put[i])
    end
end

#=Primjer iz predavanja 8 strana 10
Z=121
put = A - B - C -G=#

A = ["A", "B", "C", "D", "E", "F", "G"]
P=["-","A","B","A","D","E","C,F"]
T=[25,30,60,1,50,4,6]
Z,put = cpm(A, P, T)
println("Z = ",Z)
print("put = ")
for i in eachindex(put)
    if i!=length(put)
        print(put[i])
        print(" - ")
    else 
        println(put[i])
    end
end

#=primjer iz predavanja 8 strana 17
Z=11
put = B - E - G=#

A = ["A", "B", "C", "D", "E", "F", "G", "H", "I"]
P = ["-", "-", "-", "A", "A,B", "C", "D,E,F", "C", "H"]
T = [2,5,4,4,2,2,4,1,2]
Z, put = cpm(A, P, T);
println("Z = ",Z)
print("put = ")
for i in eachindex(put)
    if i!=length(put)
        print(put[i])
        print(" - ")
    else 
        println(put[i])
    end
end

#=primjer iz zsr strana 1
Z=38
put= A-E-G-J-K=# 
A = ["A", "B", "C", "D", "E", "F", "G", "H", "I","J","K"]
P = ["-", "-", "-", "A", "A", "C", "B,E", "D", "F","F,G,H","I,J"]
T = [4,10,12,6,8,8,10,10,8,10,6]
Z, put = cpm(A, P, T);
println("Z = ",Z)
print("put = ")
for i in eachindex(put)
    if i!=length(put)
        print(put[i])
        print(" - ")
    else 
        println(put[i])
    end
end

#=primjer iz zsr strana 2
Z=14
put= A -C -H=#

A = ["A", "B", "C", "D", "E", "F", "G", "H", "I","J"]
P = ["-", "-", "A", "A", "B", "C", "C,D", "C", "F,G","D,E"]
T = [4,5,3,5,3,4,3,7,2,5]
Z, put = cpm(A, P, T);
println("Z = ",Z)
print("put = ")
for i in eachindex(put)
    if i!=length(put)
        print(put[i])
        print(" - ")
    else 
        println(put[i])
    end
end








