#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

/* RETURN 1 - Apenas encerra a função main(). Se usado em outra função,
 encerra aquela função, não o programa todo.

 EXIT(1) - Encerra imediatamente todo o programa, independente de onde for chamado.
 Libera memória, fecha arquivos abertos, e executa funções registradas com atexit().
 
 */
const int MAX_ITERACOES = 1000;

double f(double x)
{
    // x3 - 9x +3 -- Exemplo de função
     return exp(- pow(x,2)) - cos(x); // Exemplo de função
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

    if (fabs(f(x2)) < epsilon || fabs(x2 - x1) < epsilon)//fabs(x2 - x1) < epsilon) parou de variar significativamente 
    {
        cout << "\nTotal de iteracoes: " << k << endl;
        return x2;
    }

    return metodoSecanteRec(x1, x2, epsilon, k);
}

int main()
{
    double x0, x1, epsilon;
    cout << "\n\nMetodo da Secante (Recursivo, com limite de iteracoes)\n\n";
    cout << "Digite dois valores iniciais:\nDigite x0: ";
    cin >> x0;
    cout << "Digite x1: ";
    cin >> x1;
    cout << "Digite a precisao epsilon: ";
    cin >> epsilon;
    cout << endl;

    double raiz = metodoSecanteRec(x0, x1, epsilon);
    cout << "\nRaiz aproximada: " << raiz << endl;

    return 0;
}
