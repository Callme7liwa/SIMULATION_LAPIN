
#include <stdio.h>
#include <stdlib.h>

    Function CalculerNouveauxNes(nombreLapinesAdultes)
        nombreNouveauxNes = 0
        Pour chaque lapine adulte
            nombrePortees = ObtenirNombrePorteesParAn(lapine)
            Pour chaque portée
                nombreNouveauxNes = nombreNouveauxNes + GenererNombreLapereaux()
            Fin Pour
        Fin Pour
        Retourner nombreNouveauxNes
    finfunction




