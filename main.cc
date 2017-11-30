#include <cstdlib>
#include <unistd.h>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <algorithm>    // std::sort
#include <mpi/mpi.h>


using namespace std;

/**
 * Fait la reunion entre le tableau inférieur et supérieur dans le global
 * @param dataInf tableau des valeurs inférieur
 * @param tailleInf taille du tableau inférieur
 * @param dataSup tableau des valeurs supérieur
 * @param tailleSup taille du tableau supérieur
 * @param data tableau global
 * @param taille taille du tableau global
 */
void reunion(int *dataInf, int tailleInf, int *dataSup, int tailleSup, int *&data, int &taille) {
    taille = tailleInf + tailleSup;

    delete[] data;
    data = new int[taille];

    int inf = 0;
    int sup = 0;

    for (int i = 0; i < taille; i++) {
        if (sup >= tailleSup || (inf < tailleInf && dataInf[inf] < dataSup[sup])) {
            data[i] = dataInf[inf];
            inf++;
        } else {
            data[i] = dataSup[sup];
            sup++;
        }

    }

}

/**
 * Partitionne le tableau global dans les tableau supérieur ou inférieur en fonction du pivot
 * @param pivot
 * @param data tableau global
 * @param taille taille du tableau global
 * @param dataInf tableau des valeurs inférieur
 * @param tailleInf taille du tableau inférieur
 * @param dataSup tableau des valeurs supérieur
 * @param tailleSup taille du tableau supérieur
 */
void partition(int pivot, int *data, int taille, int *&dataInf, int &tailleInf, int *&dataSup, int &tailleSup) {
    tailleInf = 0;
    tailleSup = 0;

    for (int i = 0; i < taille; i++) {
        if (data[i] < pivot) {
            tailleInf++;
        } else {
            tailleSup++;
        }
    }

    delete[] dataInf;
    delete[] dataSup;
    dataInf = new int[tailleInf];
    dataSup = new int[tailleSup];

    int inf = 0;
    int sup = 0;

    for (int i = 0; i < taille; i++) {
        if (data[i] < pivot) {
            dataInf[inf] = data[i];
            inf++;
        } else {
            dataSup[sup] = data[i];
            sup++;
        }
    }
}

/**
 * Echange le tableau supérieur ou inférieur au noeud correspondant
 * @param data tableau inférieur ou supérieur dépendant du myPE
 * @param taille taille du tableau data
 * @param etape numéro de l'étape
 */
void exchange(int *&data, int &taille, int etape) {
    int p, myPE;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &myPE);

    //echange le bit de l'étape pour savoir à qui envoyer
    int numEnvoyeur = myPE ^(0x1 << etape - 1);

    MPI_Send(&taille, 1, MPI_INT, numEnvoyeur, 666, MPI_COMM_WORLD);
    if (taille > 0) {
        MPI_Send(data, taille, MPI_INT, numEnvoyeur, 666, MPI_COMM_WORLD);
    }

    MPI_Recv(&taille, 1, MPI_INT, numEnvoyeur, 666, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (taille > 0) {
        delete[] data;
        data = new int[taille];
        MPI_Recv(data, taille, MPI_INT, numEnvoyeur, 666, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

/**
 * Diffuse le pivot a tout le monde en fonction de l'étape
 * @param pivot
 * @param etape numéro de l'étape
 * @param data tableau global
 * @param taille taille du tableau global
 */
void diffusion(int &pivot, int etape, int *data, int &taille) {
    int p, myPE;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &myPE);

    int diffuseur;
    if (etape == log2(p)) {
        diffuseur = 0;
    } else {
        diffuseur = ((myPE >> etape) << etape);
    }

    if(myPE == diffuseur){
        pivot = data[taille / 2];
    }

    for (int i = 0; i < etape; i++) {
        if ((myPE - diffuseur) < (0x1 << i)) {
            MPI_Send(&pivot, 1, MPI_INT, myPE + (0x1 << i), 666, MPI_COMM_WORLD);
        } else if ((myPE - diffuseur) < (0x1 << (i + 1))) {
            MPI_Recv(&pivot, 1, MPI_INT, myPE - (0x1 << i), 666, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
}

/**
 * Realise le quickSort en appelant toutes les fonctions
 * @param data tableau global de données
 * @param taille taille du tableau global
 */
void quickSort(int *&data, int &taille) {
    int p, myPE;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &myPE);

    //tri local
    sort(data, data + taille);

    int dimension = log2(p);
    int pivot = 0;
    int *dataInf = new int;
    int *dataSup = new int;
    int tailleInf;
    int tailleSup;
    for (int i = dimension; i > 0; i--) {
        diffusion(pivot, i, data, taille);
        partition(pivot, data, taille, dataInf, tailleInf, dataSup, tailleSup);
        //Si le bit est de poid faible
        if (!(myPE >> (i - 1) & 0x1)) {
            exchange(dataSup, tailleSup, i);
        } else {
            exchange(dataInf, tailleInf, i);
        }

        reunion(dataInf, tailleInf, dataSup, tailleSup, data, taille);

    }
}

/**
 * Affiche le resultat du quickSort
 * @param data tableau global de données
 * @param taille taille du tableau global
 */
void printAll(int *data, int taille) {
    int p, myPE;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &myPE);
    int *recvcounts, *displs, *recvbuf;
    if (myPE == 0) recvcounts = new int[p];
    MPI_Gather(&taille, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (myPE == 0) {
        displs = new int[p];
        displs[0] = 0;
        for (int pe = 1; pe < p; pe++) displs[pe] = displs[pe - 1] + recvcounts[pe - 1];
        recvbuf = new int[displs[p - 1] + recvcounts[p - 1]];
    }
    MPI_Gatherv(data, taille, MPI_INT, recvbuf, recvcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);
    if (myPE == 0)
        for (int k = 0; k < displs[p - 1] + recvcounts[p - 1]; k++) cout << recvbuf[k] << endl;
    if (myPE == 0) delete recvbuf, recvcounts, displs;
}

int random(int min, int max){
    return rand()%(max-min + 1) + min;
}

void test_reunion(){
    int size = random(10,1000);
    auto *array1 = new int[size];
    auto *array2 = new int[size];
    for (int i = 0; i < size; ++i) {
        array1[i] = i;
        array2[i] = i+5;
    }
    auto *result = new int[size*2];
    int size2 = size*2;
    reunion(array1, size, array2, size, result, size2);
    auto *resultMerge = new int[size*2];
    merge (array1,array1+size,array2,array2+size,resultMerge);
    if(equal(result, result+size2, resultMerge)){
        cout << "reunion is ok" << endl;
    }else {
        cout << "reunion is not ok" << endl;
    }
}

void test_partition(){
    int pivot = random(10,1000);
    int size = random(10,1000);
    auto *array1 = new int[size];
    auto *array2 = new int[size];
    for (int i = 0; i < size; ++i) {
        array1[i] = random(0, pivot-1);
        array2[i] = random(pivot, 1000);
    }
    sort(array1, array1+size);
    sort(array2, array2+size);
    auto *resultMerge = new int[size*2];
    merge (array1,array1+size,array2,array2+size,resultMerge);
    auto *dataInf = new int;
    auto *dataSup = new int;
    int tailleInf;
    int tailleSup;
    partition(pivot, resultMerge, size*2, dataInf, tailleInf, dataSup, tailleSup);
    cout << "partition ";
    if(equal(array1, array1+size, dataInf)){
        cout << "dataInf is ok" << endl;
    }else {
        cout << "dataInf is not ok" << endl;
        for (int j = 0; j < tailleInf; ++j) {
            cout << dataInf[j]  << "==" << array1[j] << endl;
        }
    }
    cout << "partition ";
    if(equal(array2, array2+size, dataSup)){
        cout << "dataSup is ok" << endl;
    }else {
        cout << "dataSup is not ok" << endl;
        for (int j = 0; j < tailleSup; ++j) {
            cout << dataSup[j]  << "==" << array2[j] << endl;
        }
    }



}

void tests(){
    /*test_reunion();
    test_partition();*/
}


int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int myPE;
    MPI_Comm_rank(MPI_COMM_WORLD, &myPE);
    int seed = atoi(argv[1]) + myPE;
    srand(seed);
    int tailleLoc = atoi(argv[2]); // tailleLoc vaut n/p
    int *dataLoc = new int[tailleLoc];
    for (int k = 0; k < tailleLoc; k++) dataLoc[k] = rand() % 1000;
    quickSort(dataLoc, tailleLoc);
    printAll(dataLoc, tailleLoc);
    MPI_Finalize();
    if(myPE == 0){
        tests();
    }
    return 0;
}


