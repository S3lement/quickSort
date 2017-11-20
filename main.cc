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

void reunion(int* dataInf, int tailleInf, int* dataSup, int tailleSup, int*& data, int& taille) {
    taille = tailleInf + tailleSup;

    delete[] data;
    data = new int[taille];

    int inf = 0;
	 int sup = 0;

    for (int i = 0; i < taille; i++) {
		  if(sup >= tailleSup || (inf < tailleInf && dataInf[inf] < dataSup[sup])){
		     data[i] = dataInf[inf];
			  inf++;
		  }else{
			  data[i] = dataSup[sup];
			  sup++;
		  }
        
    }

}

void partition(int pivot, int* data, int taille, int*& dataInf, int& tailleInf, int*& dataSup, int& tailleSup) {
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

void exchange(int*& data, int& taille, int etape) {
    int p, myPE;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &myPE);

    int numEnvoyeur = myPE ^ (0x1 << etape - 1);

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

void diffusion(int& pivot, int etape, int* data, int& taille) {
    /*J'ai un bug dans ma fonction
	 int p, myPE;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &myPE);
    int nbSender = -1;
    int dec;
    int dimension = log2(p);

    if ((myPE % (int) pow(2.0, (double) etape)) == 0) {
        pivot = data[taille / 2];
        nbSender = myPE;
    }

    if (nbSender == -1) {
        MPI_Recv(&pivot, 1, MPI_INT, MPI_ANY_SOURCE, 666, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //cout << "recv " << pivot << endl;
    }

    int nbSousHyperCube = pow(2, dimension - etape);
    for (int i = 0; i < etape; i++) {
        for (int j = 0; j < nbSousHyperCube; j++) {
            int dec = pow(2, i);
            if (myPE < dec + (pow(2, etape) * j)) {
                MPI_Send(&pivot, 1, MPI_INT, myPE + dec, 666, MPI_COMM_WORLD);
                //cout << "send to " << myPE + dec << " by " << myPE << endl;
            }
        }
    }*/
	 // Send pivot depending on the etape
    int p, myPE;
    MPI_Comm_size(MPI_COMM_WORLD,&p);
    MPI_Comm_rank(MPI_COMM_WORLD,&myPE);

    // Root is the process of broadcast
    int root;
    if (etape == log2(p))
        root = 0;
    else
        root = ((myPE >> etape) << etape);

    // myPE - root = relative index in sub-hypercube of dimension i
    for (int k = 0; k < etape; k++) {
        if ((myPE - root) < (0x1 << k)) {
            // Send the pivot
            MPI_Send(&pivot, 1, MPI_INT, myPE + (0x1 << k), 666, MPI_COMM_WORLD);
        }
        else if ((myPE - root) < (0x1 << (k + 1))) {
            // Receive the pivot
            MPI_Recv(&pivot, 1, MPI_INT, myPE - (0x1 << k), 666, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
}

void quickSort(int*& data, int& taille) {
    int p, myPE;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &myPE);

    //tri local
    sort(data, data + taille);

    int dimension = log2(p);
    int pivot = 0;
    int* dataInf = new int;
    int* dataSup = new int;
    int tailleInf;
    int tailleSup;
    for (int i = dimension; i > 0; i--) {
		  if (taille != 0) {
            pivot = data[taille / 2];
        }
        
		  diffusion(pivot, i, data, taille);
        
		  partition(pivot, data, taille, dataInf, tailleInf, dataSup, tailleSup);

        if (!(myPE >> (i - 1) & 0x1)) {
            exchange(dataSup, tailleSup, i);
        } else {
            exchange(dataInf, tailleInf, i);
        }

        reunion(dataInf, tailleInf, dataSup, tailleSup, data, taille);
    }
}

void printAll(int* data,int taille) {
    int p,myPE;
    MPI_Comm_size(MPI_COMM_WORLD,&p);
    MPI_Comm_rank(MPI_COMM_WORLD,&myPE);
    int *recvcounts, *displs, *recvbuf;
    if (myPE == 0) recvcounts = new int[p];
    MPI_Gather(&taille,1,MPI_INT,recvcounts,1,MPI_INT,0,MPI_COMM_WORLD);
    if (myPE == 0) {
        displs = new int[p];
        displs[0] = 0;
        for (int pe=1;pe<p;pe++) displs[pe] = displs[pe-1]+recvcounts[pe-1];
        recvbuf = new int[displs[p-1]+recvcounts[p-1]];
    }
    MPI_Gatherv(data,taille,MPI_INT,recvbuf,recvcounts,displs,MPI_INT,0,MPI_COMM_WORLD);
    if (myPE == 0)
        for (int k=0;k<displs[p-1]+recvcounts[p-1];k++) cout << recvbuf[k] << endl;
    if (myPE == 0) delete recvbuf,recvcounts,displs;
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
    return 0;
}
