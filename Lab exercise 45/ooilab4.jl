using LinearAlgebra
#=
radili skupa
Amina Kazazović 19364
Daris Mujkić 19413
=#

function Validacija(A,b,c,csigns,vsigns)
    broj_redova, broj_kolona= size(A) 
    if length(b)!=broj_redova || length(c)!=broj_kolona
        error("Neadekvatne dimenzije matrice/vektora")
    end
    #csigns oznacava kakvo je koje ogranicenje jel <=, >= ili =
    #zbog toga mora imat iste dimenzije kao i b 
    if(length(csigns)!=length(b)) error("Neadekvatne dimenzije vektora csigns ili b") end 
    #vsigns oznacava ogranicenje na xi, da li su >=,<= ili =
    #zbog toga ih mora imat isto kao i kolona od A 
    if(broj_kolona!=length(vsigns)) error("Neadekvatne dimenzije vektora vsigns ili matrice A") end
end

function postojiTaVjestacka(i,vektorBaze)
    for (indeks,tip) in vektorBaze
        if(i-1==indeks && tip==1) return true end
    end
    return false
end

function rijesi_simplex(goal,A,b,c,csigns,vsigns)
    Validacija(A,b,c,csigns,vsigns)
    broj_redova,broj_kolona=size(A)
    #sad slijedi namjestanje u slucaju da je neki od x-eva negativan ili =0
    neograniceneVarijable=Tuple[]
    for i in 1:length(vsigns)
        #ako je negativan onda zamijenit sa drugom varijablom x'=-x 
        if vsigns[i] == -1
			A[:, i] *= -1
			c[i] *= -1
        elseif vsigns[i]==0 
            #ako je =0 onda imamo x=x''-x' 
            #zato dodamo dvije 
            c = [c -c[i]]
			A = [A -A[:, i]]
			push!( neograniceneVarijable, (i, size(A, 2))) #prvi je indeks varijable te a drugo je indeks te dodane 
            #originalna se dobije tako što oduzmem te dvije 
        end
    end
    #namjestanje u slucaju da je neki bi negativan
    for i in 1:broj_redova
        if b[i]<0 
            b[i]=-b[i]
            A[i,:]=-A[i,:]
            csigns[i]=-csigns[i]
        end
    end
    #namjestanje sada za csigns dodavanje vjestackih i dopunskih
    vektorDopunskih=Tuple[] #prva vrijednost je indeks dopunske, indeks ogranicenja u kojem je  i njen koeficijent
    vektorVjestackih=Tuple[] #isto ko za gore
    #namjestanje ogranicenja u slucaju da je >= onda i dopunska i vjestacka idu
    #u slucaju da je <= onda samo dopunska
    #u slucaju da je = onda samo vjestacka
    #csigns je constraintsigns
    indeksDopunskih=broj_kolona;
    #prvo dodavanje dopunskih
    for i in 1:length(csigns)
        if(csigns[i]==+1) 
            indeksDopunskih+=1
            push!(vektorDopunskih,(indeksDopunskih,i,-1))
        elseif (csigns[i]==-1)
            #ako je <=
            indeksDopunskih+=1
            push!(vektorDopunskih,(indeksDopunskih,i,1))
        end
    end
    indeksVjestackih=indeksDopunskih
    for i in 1:length(csigns)
        if(csigns[i]==+1) 
            indeksVjestackih+=1
            push!(vektorVjestackih,(indeksVjestackih,i,1))
        elseif(csigns[i]==0)
            indeksVjestackih+=1
            push!(vektorVjestackih,(indeksVjestackih,i,1))
        end
    end
    minimizacija=false
    if(goal=="Min" || goal=="min")
        c= -1 .*c #pomnozit je cijelu s minus 1 da postane maksimizacija
        minimizacija=true
    end

    #pravljenje simpleks matrice
    simpleks_matrica=hcat(b,A)
    #dodat s desne strane ove od dopunskih tako što ću napravit matricu od vektora pojedinacnih
    for (indeks, i, koef) in vektorDopunskih
        simpleks_matrica = hcat(simpleks_matrica, zeros(broj_redova))  # Dodaj novu kolonu
        for j in 1:broj_redova
            if i == j
                simpleks_matrica[j, end] = koef  # Postavi koeficijent u odgovarajući red
            end
        end
    end

    # dodavanje vjestackih
    for (indeks, i, koef) in vektorVjestackih
        simpleks_matrica = hcat(simpleks_matrica, zeros(broj_redova))  # Dodaj novu kolonu
        for j in 1:broj_redova
            if i == j
                simpleks_matrica[j, end] = koef  # Postavi koeficijent u odgovarajući red
            end
        end
    end
    #sada treba da dodamo vektor koji ce oznacavat ko nam je u bazi
    vektorBaze=Tuple[] #prva vrijednost je indekspromjenljive, a druga jel dopunska il vjestacka il obicna (1 za vjestacku, za dopunsku 0 i -1 za obicnu)
    for i in 1:broj_redova
        # Provjera da li je red u vektoru dopunskih ili vještackih
        postojiRedDopunske = false
        postojiRedVjestacke = false
        indeksDopunske=0
        indeksVjestacke=0
        for (indeks, i_dopunski, _) in vektorDopunskih
            if i_dopunski == i
                postojiRedDopunske = true
                indeksDopunske=indeks
                break
            end
        end
        for (indeks, i_vjestacki, _) in vektorVjestackih
            if i_vjestacki == i
                postojiRedVjestacke = true
                indeksVjestacke=indeks
                break
            end
        end

        if postojiRedDopunske && postojiRedVjestacke
            push!(vektorBaze, (indeksVjestacke, 1))  # Vjestacka
        elseif postojiRedDopunske && !postojiRedVjestacke
            push!(vektorBaze, (indeksDopunske, 0))  # Dopunska
        elseif postojiRedVjestacke
            push!(vektorBaze,(indeksVjestacke,1)) #Vjestacka
        end
    end
    vektorVjestackih2=[]
    vektorDopunskih2=[]
    vektorBaze2=[]
    for i in 1:length(vektorBaze) 
        indeks,tip=vektorBaze[i]
        push!(vektorBaze2,indeks)
    end

    for i in 1:length(vektorVjestackih) 
        indeks,tip=vektorVjestackih[i]
        push!(vektorVjestackih2,indeks)
    end

    for i in 1:length(vektorDopunskih)
        indeks,tip=vektorDopunskih[i]
        push!(vektorDopunskih2,indeks)
    end
    #sad imamo formiranu simpleks matricu i trebamo zapocet algoritam 
    imaVjestackih=true
    if isempty(vektorVjestackih2)
        imaVjestackih=false 
    end

    broj_redovaSim,broj_kolonaSim=size(simpleks_matrica)
    vektorM=Float64[]
    if imaVjestackih
    for j in 1:broj_kolonaSim
        zbir=0
        if(postojiTaVjestacka(j,vektorBaze))
            zbir=0
        else
        for i in 1:broj_redovaSim 
            indeksVar,tip=vektorBaze[i]
            if(tip==1)
                zbir+=simpleks_matrica[i,j]
            end
        end
    end
        push!(vektorM,zbir)
    end
    simpleks_matrica=vcat(simpleks_matrica,vektorM')
    end

    if imaVjestackih==false
        prazna=collect(0 for i in 1:size(simpleks_matrica,2))
        simpleks_matrica=vcat(simpleks_matrica,prazna')
    end
    vel1=length(vektorDopunskih2)
    vel2=length(vektorVjestackih2)
    vel=vel1+vel2
    c_nova=vcat(0,c',collect(0 for i in 1:vel))
    simpleks_matrica=vcat(simpleks_matrica,c_nova')
    brojac=0;

	redM = deepcopy(simpleks_matrica[end-1, :])
	popfirst!(redM)
	cMaxM, indeks_koloneM = findmax(redM)
	indeks_koloneM =indeks_koloneM+1

	indeks_koloneZadnjiRed = 0
	predzadnjiRed = simpleks_matrica[end-1, :]
	zadnjiRed = simpleks_matrica[end, :]
	cMax = -Inf
	for i in 2:lastindex(zadnjiRed)
		if zadnjiRed[i] > cMax && (predzadnjiRed[i] >= 0 || predzadnjiRed[i] == -0)
			cMax = zadnjiRed[i]
			indeks_koloneZadnjiRed = i
		end
	end

	while cMax > 0 || cMaxM > 0
		if cMaxM > 0
			indeks_kolone = indeks_koloneM
		else
			indeks_kolone = indeks_koloneZadnjiRed
		end

		t_max, indeks_reda = Inf,-1

        for i in 1:(size(simpleks_matrica, 1) - 2)
            kolona_vrijednost = simpleks_matrica[i, indeks_kolone]
            if kolona_vrijednost > 0
                t_temp = simpleks_matrica[i, 1] / kolona_vrijednost
                if t_temp < t_max || (t_temp == t_max && rand() > 0.5)
                    t_max, indeks_reda = t_temp, i
                end
            end
        end

		if t_max == Inf
			throw("Rjesenje je neograniceno");
		end
		vektorBaze2[indeks_reda] = indeks_kolone - 1

		simpleks_matrica[indeks_reda, :] ./= simpleks_matrica[indeks_reda,indeks_kolone] #dijeljenje reda s pivotom

		for i in 1:size(simpleks_matrica, 1)
			if i != indeks_reda
				faktor = simpleks_matrica[i, indeks_kolone]
				for j in 1:size(simpleks_matrica, 2)
					simpleks_matrica[i, j] -= simpleks_matrica[indeks_reda, j] * faktor
				end
			end
		end

		redM = deepcopy(simpleks_matrica[end-1, :])
		popfirst!(redM)
		cMaxM, indeks_koloneM = findmax(redM)
		indeks_koloneM += 1

		if cMaxM <= 1e-9
			cMaxM = 0
		end

		if cMaxM <= 0
			predzadnjiRed = simpleks_matrica[end-1, :]
			zadnjiRed = simpleks_matrica[end, :]
			cMax = -Inf
			for i in 2:lastindex(zadnjiRed)
				if zadnjiRed[i] > cMax && (predzadnjiRed[i] >= 0 || predzadnjiRed[i] == -0)
					cMax = zadnjiRed[i]
					indeks_koloneZadnjiRed = i
				end
			end
		end
	end

    #provjera imal vjestackih u bazi
    if imaVjestackih
        for i in 1:length(vektorVjestackih2)
            if vektorVjestackih2[i] in vektorBaze2
                error("Dopustiva oblast ne postoji jer imamo vještackih u bazi!")
            end
        end
    end

    x = zeros(1, size(b, 1) + size(c, 2))
    for i in 1:lastindex(vektorBaze2)
		x[vektorBaze2[i]] = simpleks_matrica[i, 1]
	end

	for i in 1:(lastindex(simpleks_matrica[:, 1])-2)
		if simpleks_matrica[i, 1] == 0
			println("Rjesenje je degenerirano")
            break
		end
	end

    println("Rjesenje nije degenerirano")

    jedinstveno = all(x[i-1] != 0 || simpleks_matrica[end, i] != 0 for i in 2:(lastindex(simpleks_matrica[end, :])-lastindex(vektorVjestackih2)))
    
    if jedinstveno
        println("Rjesenje je jedinstveno")
    else println("Rjesenje nije jedinstveno")
    end 

    if !isempty(neograniceneVarijable)
        for varijabla in neograniceneVarijable
            prviEl, drugiEl = findall(==(varijabla[1]), x), findall(==(varijabla[2]), x)
            
            if isempty(prviEl) && !isempty(drugiEl)
                x[varijabla[1]] = -x[drugiEl[1]]
                deleteat!(x, drugiEl[1])
            end
        end
    end

    Z=simpleks_matrica[end,1]
    if goal=="max"
        Z=-Z
    end
    
    return Z,x
end

rijesi_simplex("max", [[0.3 0.1]; [0.5 0.5];[0.6 0.4]], [2.7,6,6], [0.4 0.5], [-1,0,1],[1,1])

#primjer1 iz predavanja strana 38
#=
argmaxZ = 3x1+x2
p.o
0.5x1+0.3x2<=150
0.1x1+0.2x2<=60
x1>=0,x2>=0  strana 34 
Rjesenje: 
Z = 900
X= (300, 0, 0, 30)
=#

Z,x=rijesi_simplex("max", [[0.5 0.3]; [0.1 0.2]], [150, 60], [3 1], [-1,-1],[1,1])
println("Z = ",Z)
println("X = ",x)


#=primjer2 iz predavanja strana 39
argmaxZ=800X1+1000X2
P.O 
30x1+16x2<=22800
14x1+19x2<=14100
11x1+26x2<=15950
x2<=550
x1>=0 i x2>=0
Rjesenje:
Z = 780000
X = (600, 300, 0, 0, 1550, 250)
=# 

Z,x=rijesi_simplex("max", [[30 16]; [14 19]; [11 26]; [0 1]], [22800, 14100, 15950, 550], [800 1000], [-1, -1, -1, -1], [1,1])
println("Z = ",Z)
println("X = ",x)


#zadatak sa stranice 53 
#=
argminZ=40x1+30x2
0.1x1>=0.2
0.1x2>=0.3
0.5x1+0.3x2>=1.2
0.1x1+0.2x2>=1.2
x1>=0 i x2>=0
Rjesenje: 
Z = 265.714
X = (3.43, 4.3, 0.14, 0.12, 0, 0)
=#

Z,x=rijesi_simplex("min", [0.1 0; 0 0.1; 0.5 0.3; 0.1 0.2], [0.2,0.3,3,1.2], [40 30], [1,1,1,1], [1,1])
println("Z = ",Z)
println("X = ",x)


#=zadatak iz predavanja stranica 57
argminZ = 32x1+56x2+50x3+60x4
x1+x2+x3+x4=1
250x1+150x2+400x3+200x4>=300
x4<=0.3
x2+x3<=0.5
x1>=0 x2>=0 x3>=0 x4>=0 
Rjesenje:
Z = 38
X = (0.67, 0, 0.33, 0, 0, 0.3, 0.167, 0)
=#

Z,x=rijesi_simplex("min", [[1 1 1 1]; [250 150 400 200]; [0 0 0 1]; [0 1 1 0]],[1, 300, 0.3, 0.5], [32 56 50 60],[0, 1, -1, -1], [1, 1, 1, 1])
println("Z = ",Z)
println("X = ",x)


#=Zadatak iz Aminine zadaće 4
argminZ=3.5x1+5x2 
5.8x1+3.5x2>=422.3
3.5x1+0.5x2<=2095
4.5x1+4.5x2=422.3
x1>=0 x2>=0
=#

Z,x= rijesi_simplex("min",[5.8 3.5; 3.5 0.5; 4.5 4.5],[422.3, 2095, 422.3],[3.5 5],[1 -1 0],[1 1]);
println("Z = ",Z)
println("X = ",x)

#Dual od proslog zadatka 
#=
argmaxZ=422.3y1+2095y2+422.3y3 
5.8y1+3.5y2+4.5y3<=3.5
3.5y1+0.5y2+4.5y3<=5
y1>=0, y2<=0 i y3 neogr po znaku
=#

Z,x=rijesi_simplex("max",[5.8 3.5 4.5; 3.5 0.5 4.5],[3.5, 5],[422.3 2095 422.3],[-1 -1],[1 -1 0])
println("Z = ",Z)
println("X = ",x)


