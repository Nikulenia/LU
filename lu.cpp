#include <iostream>
#include <vector>

using namespace std;


// Функция для нахождения LU-разложения матрицы
void LU_decomposition(vector<vector<double>>& A, vector<vector<double>>& L,vector<vector<double>>& U) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        for (int k = i; k < n; k++) {
            double sum = 0;
            for (int j = 0; j < i; j++) {
                sum += (L[i][j] * U[j][k]);
            }
            U[i][k] = A[i][k] - sum;
        }
        for (int k = i; k < n; k++) {
            if (i == k) {
                L[i][i] = 1;
            }
            else {
                double sum = 0;
                for (int j = 0; j < i; j++) {
                    sum += (L[k][j] * U[j][i]);
                }
                L[k][i] = (A[k][i] - sum) / U[i][i];
            }
        }
    }
}

// Функция для решения СЛАУ методом LU-разложения
vector<double> solve_LU(vector<vector<double>>& L, vector<vector<double>>& U,vector<double>& b) {
    int n = L.size();
    vector<double> y(n, 0);
    vector<double> x(n, 0);

    // Решаем Ly = b
    for (int i = 0; i < n; i++) {
        double sum = 0;
        for (int j = 0; j < i; j++) {
            sum += L[i][j] * y[j];
        }
        y[i] = b[i] - sum;
    }

    // Решаем Ux = y
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += U[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / U[i][i];
    }

    return x;
}

int main() {
    setlocale(LC_ALL, "ru");

    int n;
    cout << "Введите размерность матрицы: ";
    cin >> n;

    cout << "Введите элементы матрицы A:" << endl;
    vector<vector<double>> A(n, vector<double>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> A[i][j];
        }
    }

    cout << "Введите элементы вектора b:" << endl;
    vector<double> b(n);
    for (int i = 0; i < n; i++) {
        cin >> b[i];
    }

    vector<vector<double>> L(n, vector<double>(n, 0));
    vector<vector<double>> U(n, vector<double>(n, 0));

    LU_decomposition(A, L, U);

    vector<double> x = solve_LU(L, U, b);

    cout << "Решение:" << endl;
    for (int i = 0; i < n; i++) {
        cout << "x_" << i + 1 << " = " << x[i] << endl;
    }

    return 0;
}
