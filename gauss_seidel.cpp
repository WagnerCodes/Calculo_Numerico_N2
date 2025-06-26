#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

const int MAX_ITER = 100;

// Função para resolver sistema n x n pelo método de Gauss-Seidel
void gaussSeidel(const vector<vector<double>> &A, const vector<double> &b, vector<double> &x, double TOLERANCIA)
{
    int n = A.size();
    int iter = 0;
    bool convergiu = false;

    while (iter < MAX_ITER && !convergiu)
    {
        convergiu = true;
        vector<double> x_ant = x;

        for (int i = 0; i < n; i++)
        {
            double soma = 0.0;
            for (int j = 0; j < n; j++)
            {
                if (i != j)
                {
                    soma += A[i][j] * x[j]; // usa valores já atualizados
                }
            }

            double x_novo = (b[i] - soma) / A[i][i];

            if (fabs(x_novo - x[i]) > TOLERANCIA)
            {
                convergiu = false;
            }

            x[i] = x_novo;
        }

        iter++;
    }

    cout << endl << endl;
    cout << "===========================================================" << endl;
    cout << " Método de Gauss-Seidel" << endl;
    cout << "===========================================================" << endl;

    if (convergiu)
    {
        cout << "Convergiu em " << iter << " iteracoes." << endl;
    }
    else
    {
        cout << "Nao convergiu apos " << MAX_ITER << " iteracoes." << endl;
    }

    cout << "Solucoes aproximadas:" << endl;
    for (int i = 0; i < n; i++)
    {
        cout << "x" << i + 1 << " = " << x[i] << endl;
    }
}

int main()
{
    cout << "Metodo de Gauss-Seidel para resolver sistemas lineares" << endl;

    int n;
    cout << "Digite o numero de variaveis: ";
    cin >> n;

    vector<vector<double>> A(n, vector<double>(n));
    cout << "Digite a matriz A (linha por linha):" << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cin >> A[i][j];
        }
    }

    vector<double> b(n);
    cout << "Digite o vetor b:" << endl;
    for (int i = 0; i < n; i++)
    {
        cin >> b[i];
    }

    vector<double> x(n);
    cout << "Digite o vetor chute inicial:" << endl;
    for (int i = 0; i < n; i++)
    {
        cin >> x[i];
    }

    double TOLERANCIA;
    cout << "Digite a tolerancia: ";
    cin >> TOLERANCIA;

    gaussSeidel(A, b, x, TOLERANCIA);

    return 0;
}
