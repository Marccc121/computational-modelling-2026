U blok1 dijelu zadataka proučavao sam numeričku analizu kroz pad tijela uz otpor zraka, numeričko deriviranje i gubitak značajnih znamenki. Cilj mi je bio vidjeti kako dt ili odsječci vremenske crte i preciznost mojig računica utječu na pitanje da li će moji rezultati biti manje točni, i kada će biti više.


Korištene formule:
Diferencijalna jednadžba pada:  dv/dt = g - (c_d / m) * v^2
Numerička derivacija (Forward):  f'(x) approx frac{f(x+h) - f(x)}{h}
Numerička derivacija (Central):  f'(x) approx frac{f(x+h) -f(x-h)}{2h}
Stabilna formula za 1C:  1 -cos(x) = 2sin^2(x/2)

- Rezultati iz 1. grafa nam pokazuju kako gušći raspored podataka bolje prati krivulju stepeničastim odstupanjem, dok iz 2. vidimo linije različitih nagiba ukazujući na opadanje pogrešaka, a iz tablica iz 1C učimo o gubljenju preciznosti direktnim oduzimanjem te gubitkom samog broja..


### Zadatak 1A


Zašto se rješenje poboljšava kad je dt manji?

Rješenje se poboljšava kad je dt manji zato što se vremenska vrpca koja nam je fiksno zadana (20 sekundi) dijeli na više dijelova, točnije x2 dijelova ili 40 ako prelazimo s dt=1 na dt=0.5. To dijeljenje dozvoljava našoj y-duljini da se "osvježi" duplo brže nego što bi se osvježavala da je dt veći, u drugim riječima linija postane manje isprekidana zbog veće gustoće točaka što pokazuje varijabla n_points preko n_points =int(t_max/dt +1). Dijeljenje s manjim brojem rezultira većim rezultatom koji je u ovom slučaju broj točaka na vremenskoj skali.


Zašto Eulerova metoda odstupa od analitičkog rješenja?

Eulerova metoda odstupa od analitičkog rješenja zato sto je u prirodi taj odnos (kao i svaki drugi) kontinuiran i u kontekstu ovog zadatka ima beskonačno mnogo "točaka" dok ja u ovom primjeru imam od 11 do 41 i ništa više od toga. Po tome je dobiveni rezultat logičan jer u slici grafa vidimo kako su sve 3 linije koje su derivirane eulerovom metodom zaostaju u sklapanju s crnom analitičkom linijom koja ne koristi pretpostavku koju ja koristim, a to je da brzina promjene v u odnosu na t ostaje ista dok ne dođe do druge točke mjerenja ovisno koji je dt.


Koja je veza između te pogreške i truncation error iz 4. poglavlja?

Veza je ta da recimo pogreška koju ja vidim na svom grafu kojega sam generirao sadrži zbroj svih lokalnih "truncation errora" direktno rečeno.


### Zadatak 1B


Koja metoda je točnija i koji red točnosti ima?

Metoda koja je točnija je središnja razlika zato što ima red točnosti O(H*2) a prednja ima O(h). To direktno znači da je pogreška proporcionalna kvadratu tog jednog pomaka, što bi značilo da ako ja sad smanjim h za 10 (dijelim s 10) greška bi pala za 100 puta ili 10x10 (dijeljeno sa 100). To se događa zbog truncation-a i saznajemo s grafa u kojemu O(h na kvadrat) središnje razlike ima omjer penjanja 1:2 dok druga ima 1:1.


Što se događa kad h postane jako malen? Može li premalen h pogoršati rezultat?

Kompjutori kao ovaj na kojemu su ovi grafovi napravljeni imaju preciznost floata64 kojeg nakon pretraživanja sam saznao sa ima preciznost oko 15-16 decimala tako da ako uzmem jako malen h koji rubi tim decimalama kad dođe do oduzimanja ono gubi većinu "significant digits" što se zove catastrophic cancellation. Može pogoršati rezultat jer tokom dijeljenja se taj strgani broj onda proširi i postane još lošiji.


### Zadatak 1C

  P.legend()
x          #2 * sin^2(x/2)       # 1 - cos(x)           # Razlika             
############################################################################
0.004995834721974179 # 0.004995834721974234 # 0.1        # 5.464378949326942347e-17
0.000049999583334737 # 0.000049999583334722 # 0.01       # 1.442666515763524338e-17
0.000000004999999970 # 0.000000004999999996 # 0.0001     # 2.622068904947756727e-17
0.000000000000500044 # 0.000000000000500000 # 1e-06      # 4.445029121221820637e-17
0.000000000000000000 # 0.000000000000000050 # 1e-08      # 5.000000000000000512e-17
0.000000000000000000 # 0.000000000000000000 # 1e-10      # 5.000000000000000478e-21


Zašto dva matematički jednaka izraza numerički ne daju isti rezultat?

Kao što sam već rekao, računala spremaju informacije u binarnom obliku s ograničenom memorijom za svaku tu informaciju (15-16 znač. znamenki) i kada je cos(x) izračunat s malim x koji je vrlo blizu 0 ali nije, dobije se nešto kao 0 s velikim brojem jedinica na kraju (što opet nije cijeli broj 1).


Što je cancellation error?

Cancellation error je kada se brojevi koji su jako bliski jedan drugome (konkretno u A-b, A je veći od b) oduzimaju te se ostane samo s jednim značajnim brojem u rezultatu upravo zato što su se sve druge znamenke pretvore u 0. To se dešava kada imamo puno zaredanih dekatskih prenašanja.