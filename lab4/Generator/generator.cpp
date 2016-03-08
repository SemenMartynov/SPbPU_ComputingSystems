#include <iostream>
#include <random>
#include <fstream>

const size_t matrix_size = 512;

int main(int agrc, char* argv[]) {
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 100.0);

    std::ofstream matrix("matrix.txt");

    for (size_t i = 0; i != matrix_size; ++i) {
        for(size_t j = 0; j != matrix_size; ++j) {
            if (j) {
                matrix << "\t" << dist(mt);
            } else {
                matrix << dist(mt);
            }
        }
        matrix << std::endl;
    }

    matrix.flush();
    return 0;
}
