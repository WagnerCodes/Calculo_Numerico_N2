#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

void gaussElimination(vector<vector<double>>& A, vector<double>& b, vector<double>& x) {
    int n = A.size();

    // Eliminação de Gauss (triangularização)
    for (int k = 0; k < n - 1; ++k) {
        for (int i = k + 1; i < n; ++i) {
            double m = A[i][k] / A[k][k];
            for (int j = k; j < n; ++j) {
                A[i][j] -= m * A[k][j];
            }
            b[i] -= m * b[k];
        }
    }

    // Substituição regressiva (sistema triangular superior)
    x.resize(n);
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) {
            sum += A[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / A[i][i];
    }
}

int main() {
    int n;
    cout << "Digite o número de equações (n): ";
    cin >> n;

    vector<vector<double>> A(n, vector<double>(n));
    vector<double> b(n), x;

    cout << "Digite a matriz A (" << n << "x" << n << "):\n";
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            cin >> A[i][j];

    cout << "Digite o vetor b:\n";
    for (int i = 0; i < n; ++i)
        cin >> b[i];

    gaussElimination(A, b, x);

    cout << "\nSolução aproximada do sistema:\n";
    for (int i = 0; i < n; ++i)
        cout << "x[" << i << "] = " << fixed << setprecision(6) << x[i] << endl;

    return 0;
}
