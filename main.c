#include <stdio.h>
#include "./headers/lapins.h"
#include "./headers/mersenne_twister.h"

int main()
{
    /* Initialisation du générateur de nombres aléatoires */
    unsigned long init[4] = {0x123, 0x234, 0x345, 0x456};
    int length = 4;
    init_by_array(init, length);

    /* Initialisation de la structure Rabbit */
    static struct Rabbit Rabbit;

    /* Mesure du temps de début de simulation */
    time_t start_time = 0;
    time_t end_time = 0;
    start_time = clock();

    /* Affichage du début de la simulation */
    printf("-- Début de la simulation !\n");

    /* Simulation de dix expériences avec un seul lapin initial pendant une durée prédéfinie */
    simulateMultipleExperiments(2, 0, &Rabbit, DURATION, 20);

    /* Mesure du temps de fin de simulation */
    end_time = clock();

    /* Affichage du temps écoulé de la simulation */
    printf("\n -- Le temps ecoule de la simulation : %fs\n", (double)(end_time - start_time) / CLOCKS_PER_SEC);

    /* Retour du programme avec succès */
    return EXIT_SUCCESS;

}