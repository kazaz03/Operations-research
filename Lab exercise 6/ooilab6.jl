#=Radili skupa
Amina Kazazovic 19364
Daris Mujkic 19413=#

using LinearAlgebra 

function najkraci_put(M)
    broj_redovaM,broj_kolonaM=size(M)
    putevi=zeros(broj_redovaM,3)
    for i in 1:broj_redovaM
        putevi[i,1]=i
    end
    #u prvoj koloni je cvor do kojeg se trazi najkraci_put
    #u drugoj koloni se nalaze duzine puteva 
    #u trecoj je preko kojeg je dobiven
    putevi[1,2]=0 
    for i in 2:broj_redovaM
        putevi[i,2]=Inf
    end
    putevi[1,3]=1
    for i in 2:broj_redovaM
        #trebam trazit u matrici one koje na tom indeksu za taj cvor imaju razlciito od nula 
        #znaci mogu imat vise suma i uzimam onu najmanju 
        if i==1 continue end 
        #vanjska petlja da znam za koji cvor gledam
        sume=Float64[]
        for j in 1:broj_redovaM
            #ova petlja za prolazak kroz matricu
            #j je cvor taj neki drugi 
            suma=0.0
            if i==j suma=Inf end
            if M[j,i]!=0 
                suma=putevi[j,2]+M[j,i]
            else 
                suma=Inf #u slucaju kada nema nikakvog puta od tog cvora do drugog 
            end
            push!(sume,suma)
        end      
        vrijednost,indeks=findmin(sume)
        putevi[i,2]=vrijednost 
        putevi[i,3]=indeks
    end
    return Int.(putevi)
end

M = [0 1 3 0 0 0; 0 0 2 3 0 0; 0 0 0 -4 9 0; 0 0 0 0 1 2; 0 0 0 0 0 2; 0 0 0 0 0 0];
putevi=najkraci_put(M)
for row in eachrow(putevi)
    println(row)
end

#dodatni primjeri 

#primjer1 
M=[0 2 4 0; 0 0 6 1; 0 0 0 1; 0 0 0 0]
putevi=najkraci_put(M)
for row in eachrow(putevi)
    println(row)
end

#primjer2 
M=[0 3 0 0 0 0; 0 0 4 2 -3 0; 0 0 0 0 7 0; 0 0 0 0 -3 -1; 0 0 0 0 0 1; 0 0 0 0 0 0]
putevi=najkraci_put(M)
for row in eachrow(putevi)
    println(row)
end

#primjer3
M=[0 -2 -1 3 0 0 0; 0 0 1 0 8 0 0; 0 0 0 0 0 0 0; 0 0 0 0 0 4 0; 0 0 0 0 0 7 -5; 0 0 0 0 0 0 0; 0 0 0 0 0 0 0]
putevi=najkraci_put(M)
for row in eachrow(putevi)
    println(row)
end