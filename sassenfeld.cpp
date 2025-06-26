#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

// Função para calcular os coeficientes beta do critério de Sassenfeld
bool criterioSassenfeld(const vector<vector<double>>& A, vector<double>& beta) {
    int n = A.size();
    beta.assign(n, 0.0);

    for (int i = 0; i < n; ++i) {
        double soma = 0.0;

        for (int j = 0; j < n; ++j) {
            if (j == i) continue;
            soma += fabs(A[i][j]);
            //soma += fabs(A[i][j]) * (j < i ? beta[j] : 1.0);
        }

        if (fabs(A[i][i]) < 1e-12) return false; // divisão por zero
        beta[i] = soma / fabs(A[i][i]);
    }

    for (double b : beta) {
        if (b >= 1.0) return false;
    }

    return true;
}

// Tenta permutar linhas para satisfazer o critério de Sassenfeld
bool tentaPermutarLinhas(vector<vector<double>>& A) {
    int n = A.size();
    vector<int> indices(n);
    for (int i = 0; i < n; ++i) indices[i] = i;

    do {
        vector<vector<double>> permA(n);
        for (int i = 0; i < n; ++i)
            permA[i] = A[indices[i]];

        vector<double> beta;
        if (criterioSassenfeld(permA, beta)) {
            A = permA;
            return true;
        }

    } while (next_permutation(indices.begin(), indices.end()));

    return false;
}

int main() {
    int n;
    cout << "Digite a ordem da matriz: ";
    cin >> n;

    vector<vector<double>> A(n, vector<double>(n));
    cout << "Digite a matriz dos coeficientes A (" << n << "x" << n << "):\n";
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            cin >> A[i][j];

    vector<double> beta;
    if (criterioSassenfeld(A, beta)) {
        cout << "\nCritério de Sassenfeld satisfeito.\n";
    } else {
        cout << "\nCritério de Sassenfeld NÃO satisfeito.\n";
        cout << "Tentando permutar linhas...\n";

        if (tentaPermutarLinhas(A)) {
            cout << "Permutação encontrada que satisfaz o critério.\n";
            cout << "Nova matriz A:\n";
            for (const auto& linha : A) {
                for (double val : linha)
                    cout << val << "\t";
                cout << "\n";
            }
        } else {
            cout << "Não foi possível satisfazer o critério com permutação de linhas.\n";
        }
    }

    return 0;
}
