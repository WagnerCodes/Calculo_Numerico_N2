#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

const int MAX_ITER = 100;

// Função para resolver sistema n x n pelo método de Gauss-Jacobi
void gaussJacobi(const vector<vector<double>> &A, const vector<double> &b, vector<double> &x, double TOLERANCIA)
{
    int n = A.size();
    vector<double> x_ant= x;

    int iter = 0;
    bool convergiu = false;

    while (iter < MAX_ITER && !convergiu)
    {
        convergiu = true;

        vector<double> x_novo(n, 0.0);

        for (int i = 0; i < n; i++)
        {
            double soma = 0.0;
            for (int j = 0; j < n; j++)
            {
                if (i != j)
                {
                    soma += A[i][j] * x_ant[j];
                }
            }

            x_novo[i] = (b[i] - soma) / A[i][i];

            if (fabs(x_novo[i] - x_ant[i]) > TOLERANCIA)
            {
                convergiu = false;
            }
        }

        x_ant = x_novo;  // Atualiza valores antigos
        iter++;
    }
    cout<<endl;
    cout<<endl;
    cout<<endl;
    cout<<"==========================================================="<<endl;
    cout<< " Método de Gauss-Jacobi"<<endl;
    cout<<"==========================================================="<<endl;
    // Resultado
    if (convergiu)
    {
        cout << "Convergiu em " << iter -1 << " iteracoes." << endl;
    }
    else
    {
        cout << "Nao convergiu apos " << MAX_ITER << " iteracoes." << endl;
    }

    cout << "Solucoes aproximadas:" << endl;
    for (int i = 0; i < n; i++)
    {
        cout << "x" << i + 1 << " = " << x_ant[i] << endl;
    }

    // Copia os resultados finais para o vetor x
    x = x_ant;
}

int main()
{
    cout << "Metodo de Gauss-Jacobi para resolver sistemas lineares" << endl;
    
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

    gaussJacobi(A, b, x, TOLERANCIA);

    return 0;
}
