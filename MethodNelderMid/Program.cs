using System;

class Program
{
    static void Main()
    {
        // Начальная точка
        double x = 1;
        double y = 2;
        double[] initialPoint = new double[] { x, y };

        double penaltyParam = 100;

        // Вызываем метод Нелдера-Мида
        var result = NelderMeadMinimize((point) => PenaltyFunction(point, penaltyParam), initialPoint);

        // Выводим результат
        Console.WriteLine($"Минимум функции: ({result.MinimizedPoint[0]}, {result.MinimizedPoint[1]})");
        Console.WriteLine($"Минимальное значение функции: {result.MinValue}");
        Console.WriteLine($"Количество итераций: {result.Iterations}");
    }

    // Определяем функцию, которую будем минимизировать
    static double F(double[] point)
    {
        double x = point[0];
        double y = point[1];
        return Math.Sqrt(x * x + y * y + 1) + x / 2 - y / 2;
    }

    static double Plane(double[] point)
    {
        double x = point[0];
        double y = point[1];
        return x - y + F(point) - 1;
    }

    static double PenaltyFunction(double[] point, double penaltyParam)
    {
        return F(point) + penaltyParam * Math.Pow(Plane(point), 2);
    }

    // Метод Нелдера-Мида для поиска минимума функции
    static NelderMeadResult NelderMeadMinimize(Func<double[], double> func, double[] initialGuess)
    {
        int n = initialGuess.Length;
        int maxIterations = 1000;
        double tolerance = 1e-6;

        // Инициализация начального симплекса
        double[][] simplex = new double[n + 1][];
        simplex[0] = initialGuess;
        for (int i = 0; i < n; i++)
        {
            double[] point = new double[n];
            Array.Copy(initialGuess, point, n);
            point[i] += 1.0;
            simplex[i + 1] = point;
        }

        // Основной цикл метода
        int iterations = 0;
        while (iterations < maxIterations)
        {
            // Вычисление значений функции для точек симплекса
            double[] values = new double[n + 1];
            for (int i = 0; i <= n; i++)
            {
                values[i] = func(simplex[i]);
            }

            // Сортировка точек симплекса по значениям функции
            Array.Sort(values, simplex, Comparer<double>.Create((a, b) => func(simplex[Array.IndexOf(values, a)]).CompareTo(func(simplex[Array.IndexOf(values, b)]))));

            // Проверка условия сходимости
            double maxDiff = 0;
            for (int i = 0; i < n; i++)
            {
                double diff = Math.Abs(values[i] - values[n]);
                if (diff > maxDiff)
                {
                    maxDiff = diff;
                }
            }

            if (maxDiff < tolerance)
            {
                break;
            }

            // Вычисление центроида симплекса
            double[] centroid = new double[n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    centroid[j] += simplex[i][j];
                }
            }

            for (int j = 0; j < n; j++)
            {
                centroid[j] /= n;
            }

            // Отражение
            double[] reflected = Reflect(centroid, simplex[n]);
            double reflectedValue = func(reflected);

            if (reflectedValue < values[0] && reflectedValue >= values[n - 1])
            {
                Array.Copy(reflected, simplex[n], n);
                iterations++;
                continue;
            }

            // Экспансия
            if (reflectedValue < values[n - 1])
            {
                double[] expanded = Expand(centroid, reflected);
                double expandedValue = func(expanded);

                if (expandedValue < reflectedValue)
                {
                    Array.Copy(expanded, simplex[n], n);
                    iterations++;
                }
                else
                {
                    Array.Copy(reflected, simplex[n], n);
                    iterations++;
                }

                continue;
            }

            // Контракция
            double[] contracted = Contract(centroid, simplex[n]);
            double contractedValue = func(contracted);

            if (contractedValue < values[n])
            {
                Array.Copy(contracted, simplex[n], n);
                iterations++;
                continue;
            }

            // Уменьшение
            for (int i = 1; i <= n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    simplex[i][j] = (simplex[i][j] + simplex[0][j]) / 2;
                    iterations++;
                }
            }
        }

        // Возвращаем результат
        return new NelderMeadResult
        {
            MinimizedPoint = simplex[0],
            MinValue = func(simplex[0]),
            Iterations = iterations
        };
    }

    // Операции Нелдера-Мида: Отражение, Экспансия, Контракция
    static double[] Reflect(double[] centroid, double[] point)
    {
        int n = centroid.Length;
        double[] reflected = new double[n];
        for (int i = 0; i < n; i++)
        {
            reflected[i] = 2 * centroid[i] - point[i];
        }
        return reflected;
    }

    static double[] Expand(double[] centroid, double[] reflected)
    {
        int n = centroid.Length;
        double[] expanded = new double[n];
        for (int i = 0; i < n; i++)
        {
            expanded[i] = 3 * centroid[i] - 2 * reflected[i];
        }
        return expanded;
    }

    static double[] Contract(double[] centroid, double[] point)
    {
        int n = centroid.Length;
        double[] contracted = new double[n];
        for (int i = 0; i < n; i++)
        {
            contracted[i] = 0.5 * centroid[i] + 0.5 * point[i];
        }
        return contracted;
    }
}

// Класс для хранения результата метода Нелдера-Мида
class NelderMeadResult
{
    public double[] MinimizedPoint { get; set; }
    public double MinValue { get; set; }
    public int Iterations { get; set; }
}
