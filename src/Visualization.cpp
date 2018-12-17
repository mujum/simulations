//
// Created by lad on 12/17/18.
//

#include "Visualization.h"
#include <iostream>
#include <fstream>


Visualization::Visualization(double ***PP, int size1, int size2, int size3) {
    std::ofstream fout;
    fout.open("data3d.dat");
    for(int index = 0; index < size2; index++) {
        fout    << PP[0][index][size3/2] << " "
                << PP[1][index][size3/2] << " "
                << PP[2][index][size3/2] << " "
                << PP[3][index][size3/2] << " "
                << PP[4][index][size3/2] << " "
                << PP[5][index][size3/2] << "\n";
    }
    fout.close();

    FILE *pipe = popen("gnuplot -persisit", "w");
    fprintf(pipe, "splot \"./data3d.dat\" with vectors\n");
    fclose(pipe);
}