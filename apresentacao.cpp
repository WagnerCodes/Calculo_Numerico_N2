#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>
using namespace std;

/* RETURN 1 - Apenas encerra a função main(). Se usado em outra função,
 encerra aquela função, não o programa todo.

 EXIT(1) - Encerra imediatamente todo o programa, independente de onde for chamado.
 Libera memória, fecha arquivos abertos, e executa funções registradas com atexit().

 */
const int MAX_ITERACOES = 1000;
double f(double x)
{
    // return pow(x,3) - 9*x +3;
    return exp(-pow(x, 2)) - cos(x);
}

double metodoSecanteRec(double x0, double x1, double epsilon, int k = 0)
{
    cout << fixed << setprecision(8);

    if (k > MAX_ITERACOES)
    {
        cerr << "Erro: número máximo de iteracoes excedido. O método pode não estar convergindo.\n";
        exit(1);
    }

    double fx0 = f(x0);
    double fx1 = f(x1);

    if (fabs(fx1 - fx0) < 1e-12) // risco divisão por zero ou overflow (1 × 10^−12)
    {
        cerr << "Erro: divisão por zero na fórmula da secante.\n";
        exit(1);
    }

    double x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0);
    k++;

    cout << "Iteracao " << k << ": x = " << x2 << ", f(x) = " << f(x2) << endl;

    if (fabs(f(x2)) < epsilon || fabs(x2 - x1) < epsilon) // fabs(x2 - x1) < epsilon) parou de variar significativamente
    {
        cout << "\nTotal de iteracoes: " << k << endl;
        return x2;
    }

    return metodoSecanteRec(x1, x2, epsilon, k);
}

const int MAX_ITER = 100;

// Função para resolver sistema n x n pelo método de Gauss-Jacobi
void gaussJacobi(const vector<vector<double>> &A, const vector<double> &b, vector<double> &x, double TOLERANCIA)
{
    int n = A.size();
    // vector<double> x_ant(n, 0.0);
    vector<double> x_ant = x;

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
            // if (fabs(x_novo[i] - x_ant[i])/ max(1.0,fabs(x_novo[i])) > TOLERANCIA)
            {
                convergiu = false;
            }
        }

        x_ant = x_novo; // Atualiza valores antigos
        iter++;
    }
    cout << endl;
    cout << endl;
    cout << endl;

    if (convergiu)
    {
        cout << "Convergiu em " << iter - 1 << " iteracoes." << endl;
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

    cout << endl
         << endl;

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

bool criterioSassenfeld(const vector<vector<double>> &A, vector<double> &beta)
{
    int n = A.size();
    beta.assign(n, 0.0);

    for (int i = 0; i < n; ++i)
    {
        double soma = 0.0;

        for (int j = 0; j < n; ++j)
        {
            if (j == i)
                continue;
           // soma += fabs(A[i][j]);
             soma += fabs(A[i][j]) * (j < i ? beta[j] : 1.0);
        }

        if (fabs(A[i][i]) < 1e-12)
            return false; // divisão por zero
        beta[i] = soma / fabs(A[i][i]);
    }




    for (double b : beta)
    {
        if (b >= 1.0)
            return false;
    }

    return true;
}

// Tenta permutar linhas para satisfazer o critério de Sassenfeld
bool tentaPermutarLinhas(vector<vector<double>> &A)
{
    int n = A.size();
    vector<int> indices(n);
    for (int i = 0; i < n; ++i)
        indices[i] = i;

    do
    {
        vector<vector<double>> permA(n);
        for (int i = 0; i < n; ++i)
            permA[i] = A[indices[i]];

        vector<double> beta;
        if (criterioSassenfeld(permA, beta))
        {
            A = permA;
            return true;
        }

    } while (next_permutation(indices.begin(), indices.end()));

    return false;
}

int main()
{
    double x0, x1, epsilon;
    cout << "\n\nMetodo da Secante (Recursivo, com limite de iteracoes)\n\n";

    double raiz = metodoSecanteRec(1.5, 2, 0.0001);
    cout << "\nRaiz aproximada: " << raiz << endl;
    cout << endl
         << endl;
    cout << "===========================================================" << endl;
    cout << " Método de Gauss-Jacobi" << endl;
    cout << "===========================================================" << endl;

    cout << "Metodo de Gauss-Jacobi para resolver sistemas lineares" << endl;

    int n = 3; // numero de variaveis

    vector<vector<double>> A(n, vector<double>(n));
    vector<double> values = {10, 2, 3, 1, 5, 1, 2, 3, 10};

    int k = 0;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            A[i][j] = values[k++];
        }
    }

    vector<double> b = {7, 8, 6};
    vector<double> x = {0.7, -1.6, 0.6}; // chute inicial

    double TOLERANCIA = 0.5;
    gaussJacobi(A, b, x, TOLERANCIA);
    cout << endl
         << endl;
    cout << "===========================================================" << endl;
    cout << " Método de Gauss-Seidel" << endl;
    cout << "===========================================================" << endl;

    n = 3; // numero de variaveis

    vector<vector<double>> As(n, vector<double>(n));
    vector<double> valuess = {5, 1, 1, 3, 4, 1, 3, 3, 6};

    int ks = 0;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            As[i][j] = valuess[ks++];
        }
    }

    vector<double> bs = {5, 6, 0};
    vector<double> xs = {0, 0, 0}; // chute inicial

    double TOLERANCIAs = 0.1;
    gaussSeidel(As, bs, xs, TOLERANCIAs);

    cout << endl
         << endl;
    cout << "===========================================================" << endl;
    cout << " Critério Sassenfeld" << endl;
    cout << "===========================================================" << endl;

    n = 3; // numero de variaveis

    vector<vector<double>> As2(n, vector<double>(n));
    vector<double> valuess2 = {2, 1, 3, 0, -1, 1, 1, 0, 3}; // questão da prova  --> n=3
    // vector<double> valuess2 = {1,0.5,-0.1,0.1,0.2,1,-0.2,-0.1,-0.1,-0.2,1,0.2,0.1,0.3,0.2,1};// usar n=4
    // vector<double> valuess2 = {1,3,1,2,1,1,2,3,10};
    //vector<double> valuess2 = {1,3,1,2,1,1,2,3,10};
    //vector<double> valuess2 = {1,3,1,2,1,1,2,3,10};
    //vector<double> valuess2 = {3,0,1,1,-1,0,3,1,2};
    

    int ks2 = 0;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            As2[i][j] = valuess2[ks2++];
        }
    }

    vector<double> beta;
    if (criterioSassenfeld(As2, beta))
    {
        cout << "\nCritério de Sassenfeld satisfeito.\n";
        cout << endl<< endl;
    }
    else
    {
        cout << "\nCritério de Sassenfeld NÃO satisfeito.\n";
        cout << "Tentando permutar linhas...\n";

        if (tentaPermutarLinhas(As2))
        {
            cout << "Permutação encontrada que satisfaz o critério.\n";
            cout << "Nova matriz A:\n";
            for (const auto &linha : As2)
            {
                for (double val : linha)
                    cout << val << "\t";
                cout << "\n";
            }
        }
        else
        {
            cout << "Não foi possível satisfazer o critério com permutação de linhas.\n";
            cout << endl
                 << endl;
        }
    }

    return 0;
}