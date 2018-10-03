using System;
using System.Runtime.Remoting.Messaging;

//[Д]+[И]*x+[М]*x^2+[А]*x^3+[П]+[А]*y+[В]*y^2+[Л]*y^3+[О]*y^4+[В]*y^5+[Л]+[Ь]*xy+[В]*(xy)^2+[О]*(xy)^3

//f(x)=5+10*x+14*x^2+x^3
//f(y)=17+y+3*y^2+13*y^3+16*y^4+3*y^5
//f(x,y)=13+30*xy+3*(xy)^2+16*(xy)^3
//z=f(x)+f(y)+f(x,y)=5+10*x+14*x^2+x^3+17+y+3*y^2+13*y^3+16*y^4+3*y^5+13+30*xy+3*(xy)^2+16*(xy)^3
//g(z)=x^2+y^2-1<=0
//h(z)=x+y<=0
namespace Задача_3
{
    class Point
    {
        public double X;
        public double Y;

        public Point(double x, double y)//конструктор
        {
            X = x;
            Y = y;
        }
        public static Point operator +(Point a, Point b)
        {
            return new Point(a.X + b.X, a.Y + b.Y);
        }
        public static Point operator -(Point a, Point b)
        {
            return new Point(a.X - b.X, a.Y - b.Y);
        }
        public static Point operator *(double a, Point b)
        {
            return new Point(a * b.X, a * b.Y);
        }

        public double Norm_2() //норма в квадрате
        {
            return X*X+Y*Y;
        }
    }

    class Poly//весь полином P
    {
        public int x_power, y_power;
        public double[,] coef;

        public Poly(double[,] Coef)//конструктор полинома
        {
            coef = Coef;
            x_power = coef.GetLength(0);//максимальная степень х
            y_power = coef.GetLength(1);//максимальная степень у
        }

        public double Value(Point a)//вычисление значения полинома в данной точке
        {
            double c = 0;
            for (int i = 0; i < x_power; i++)
            {
                for (int j = 0; j < y_power; j++)
                {
                    c = c + coef[i, j] * Math.Pow(a.X, i) * Math.Pow(a.Y, j);
                }
            }
            return c;
        }
        public Point Grad(Point a)//градиент от полинома
        {
            Point c = new Point(0, 0);
            for (int i = 0; i < x_power; i++)
            {
                for (int j = 0; j < y_power; j++)
                {
                    if (i!=x_power-1)
                    {
                        c.X = c.X + coef[i + 1, j] * (i + 1) * Math.Pow(a.X, i) * Math.Pow(a.Y, j);
                    }
                    if (j!=y_power-1)
                    {
                        c.Y = c.Y + coef[i, j + 1] * (j + 1) * Math.Pow(a.X, i) * Math.Pow(a.Y, j);
                    }
                }
            }
            return c;
        }
        
    }

    class Program
    {
        private const double Eps = 1e-8;
        static double[,] PolyCoef =//задание коеффициентов 
        {   /*x*/
       /*y*/{35,  1, 18,  1, 12, 16, 3, 1},//по горизонтали - фамилия
            {10, 30,  0,  0,  0,  0, 0, 0},//по вертикали - имя
            {14,  0,  3,  0,  0,  0, 0, 0},//по диагонали - отчество
            { 1,  0,  0, 16,  0,  0, 0, 0}
        };
        static Poly Polynom=new Poly(PolyCoef);//создание экземпляра класса для работы с ним
        static double[] Condition(Point a)//1 ограничение, 2 ограничение
        {
            return new[] {a.X*a.X+a.Y*a.Y-1,a.X+a.Y};
        }
        static Point[] GradCondition(Point a)//градиент 1 ограничения, градиент 2 ограничения
        {
            return new[] {new Point(2 * a.X, 2 * a.Y), new Point(1, 1)};
        }

        static double Penalty(Point a)//вычисление H
        {
            double c = 0;
            double penalty_sum=0;//Значение H
            for (int i = 0; i < 2; i++)
            {
                double check=Condition(a)[i]; //проверка условий g,h
                if (check > 0)
                {
                    penalty_sum += Math.Pow(check, 2);
                }
            }
            return penalty_sum;
        }

        static double GetPointValue(Point a,double r)//Значение f(x,y)+r*H
        {
            return Polynom.Value(a)+r * Penalty(a);
        }

        static Point PenaltyCalculate()//метод штрафов
        {
            Point currentPoint = new Point(1, 1);//выбор начальной точки
            double r = 1e6;//выбор r
            double value_CP = GetPointValue(currentPoint, r);//подсчет значения в 1 точке
            double value_PP = value_CP;//предыдущая точка = текущей
            int super_iter = 1;
            do//вычисления с одним r(k)
            {
                Console.WriteLine("Супер итерация: {0}.\n", super_iter);
                super_iter++;
                double dx = 0.5;//шаг по х
                double dy = 0.5;//шаг по у
                double epsilon = Eps*100;//точность для безусловной минимизации
                int iter = 1;
                double[] constrains = Condition(currentPoint);//вывод ограничений
                Console.WriteLine("{0})  X={1:0.00000}, Y={2:0.00000} | {3:0.00000} | {4:E2}  {5:E2} | {6:0.00000} | {7:E2}  {8:E2}", iter, currentPoint.X,currentPoint.Y, value_CP, dx, dy, Polynom.Value(currentPoint), constrains[0], constrains[1]);
                iter++;
                //безусловная минимизация методом покоординатного спуска
                do
                {
                    value_PP = value_CP;//предыдущая точка = текущей
                    Point nextPoint;//следующая точка
                    double value_NP;//значение сл.т.
                    //по х
                    do
                    {
                        nextPoint = currentPoint + new Point(dx, 0);//х + dx
                        value_NP = GetPointValue(nextPoint, r);//Значение в х + dx
                        if (value_NP > value_CP)
                        {
                            nextPoint = currentPoint - new Point(dx, 0);
                            value_NP = GetPointValue(nextPoint, r);
                            if (value_NP > value_CP)//если значение следующей точки меньше предыдущего изменяем направление движения
                                dx *= 0.5;
                            else
                                dx *= -1;
                        }
                    } while (value_NP > value_CP);//пока значение следующей точки больше предыдущей
                    while (value_NP < value_CP)//пока значение текущей точки меньше следующей
                    {
                        currentPoint = nextPoint;
                        value_CP = value_NP;
                        nextPoint = currentPoint + new Point(dx, 0);
                        value_NP = GetPointValue(nextPoint, r);
                    }
                    //по у
                    do
                    {
                        nextPoint = currentPoint + new Point(0, dy);
                        value_NP = GetPointValue(nextPoint, r);
                        if (value_NP > value_CP)
                        {
                            nextPoint = currentPoint - new Point(0, dy);
                            value_NP = GetPointValue(nextPoint, r);
                            if (value_NP > value_CP)
                                dy *= 0.5;
                            else
                                dy *= -1;
                        }
                    } while (value_NP > value_CP);
                    while (value_NP < value_CP)
                    {
                        currentPoint = nextPoint;
                        value_CP = value_NP;
                        nextPoint = currentPoint + new Point(0, dy);
                        value_NP = GetPointValue(nextPoint, r);
                    }
                    Console.WriteLine("{0})  X={1:0.00000}, Y={2:0.00000} | {3:0.00000} | {4:E2}  {5:E2} | {6:0.00000} | {7:E2}  {8:E2}", iter, currentPoint.X, currentPoint.Y, value_CP, dx, dy, Polynom.Value(currentPoint), Condition(currentPoint)[0], Condition(currentPoint)[1]);
                    iter++;
                } while ((dx >= epsilon) || (dy >= epsilon) || (value_PP - value_CP >= epsilon));
                //конец безусловной минимизации

                r *= 10;//r(k+1)=10*r(k)
            } while (Penalty(currentPoint) > Eps);//проверка EPS
            Console.WriteLine("\nPoint X={0:0.00000000},Y={1:0.00000000}, f={2}", currentPoint.X, currentPoint.Y,Polynom.Value(currentPoint));
            return currentPoint;
        }

        static double GetScalarProduct(Point point, Point antigradient)//получение скаляра для проверок
        {
            double result = Math.Abs(point.X * antigradient.X + point.Y * antigradient.Y) / (Math.Sqrt(point.Norm_2()) * Math.Sqrt(antigradient.Norm_2()));

            return result;
        }
        static void CheckPoint(Point point)
        {
            Point antigradient = -1 * Polynom.Grad(point);
            double epsilon = 1e-4;
            bool isOutOfRange = false;
            double[] Constraints = Condition(point);
            Point[] constraintsGradients = GradCondition(point);
            //вывод значения точки, антиградиента, норма антиградиента, значение в первом ограничении, значение во 2
            Console.WriteLine($"Value, antigradient of initial polynomial f(x) and antigradient norm: {Polynom.Value(point):0.00000} " +
                $"X={antigradient.X},Y={antigradient.Y} {Math.Sqrt(antigradient.Norm_2()):0.00000}.");
            Console.WriteLine($"Value: first boundary = {Constraints[0]:E2} and second boundary = {Constraints[1]:E2}.");
            //проверка, соответствует ли точка области ограничений
            if (Constraints[0] < 0 && Constraints[1] < 0)
                Console.WriteLine("\nPoint is within the allowable range.");
            else { Console.WriteLine("\nPoint is out-of-range point."); isOutOfRange = true; }
            
            if (Math.Abs(Constraints[0]) < epsilon && Math.Abs(Constraints[1]) < epsilon)
            {
                Console.WriteLine("Point is in a neighbourhood of the intersection of boundaries.");

                Point normal = new Point(1 / Math.Sqrt(2), 1 / Math.Sqrt(2));
                if (GetScalarProduct(normal, antigradient) >= 0 && GetScalarProduct(point, antigradient) >= 0)//принадлежит ли градиент конусу направлений
                    Console.WriteLine("Gradient belongs to cone of directions.");
                else Console.WriteLine("Gradient doesn't belong to cone of directions.");

                if (isOutOfRange)
                {
                    point = point - Constraints[0] / constraintsGradients[0].Norm_2() * constraintsGradients[0];
                    point = point - Constraints[1] / constraintsGradients[1].Norm_2() * constraintsGradients[1];
                }
            }

            if (Math.Abs(Constraints[0]) < epsilon && Math.Abs(Constraints[1]) > epsilon)
            {
                Console.WriteLine("Point is in a neighbourhood of the first boundary.");
                //проверка на коллинеарность
                if (GetScalarProduct(point, antigradient) - 1 < epsilon && Math.Sign(point.X) == Math.Sign(antigradient.X) && Math.Sign(point.Y) == Math.Sign(antigradient.Y))
                    Console.WriteLine("Gradient and point are collinear.");
                else Console.WriteLine("Gradient and point are not collinear.");

                if (isOutOfRange)
                    point = point - Constraints[0] / constraintsGradients[0].Norm_2() * constraintsGradients[0];
            }

            if (Math.Abs(Constraints[0]) > epsilon && Math.Abs(Constraints[1]) < epsilon)
            {
                Console.WriteLine("Point is in a neighbourhood of the second boundary.");

                Point normal = new Point(1 / Math.Sqrt(2), 1 / Math.Sqrt(2));

                if (GetScalarProduct(normal, antigradient) - 1 < epsilon && Math.Sign(normal.X) == Math.Sign(antigradient.X) && Math.Sign(normal.Y) == Math.Sign(antigradient.Y))
                    Console.WriteLine("Gradient and point are collinear.");
                else Console.WriteLine("Gradient and point are not collinear.");

                if (isOutOfRange)
                    point = point - Constraints[1] / constraintsGradients[1].Norm_2() * constraintsGradients[1];
            }

            if (isOutOfRange)
            {
                antigradient = -1 * Polynom.Grad(point);
                Constraints = Condition(point);

                Console.WriteLine($"\nChanged point: X={point.X:0.00000000} Y={point.Y:0.00000000} f={Polynom.Value(point)}.\n");
                Console.WriteLine($"Value, antigradient of initial polynomial f(x) and antigradient norm: {Polynom.Value(point):0.00000}" +
                    $" X={antigradient.X},Y={antigradient.Y} {Math.Sqrt(antigradient.Norm_2()):0.00000}.");
                Console.WriteLine($"Value: first boundary = {Constraints[0]:E2} and second boundary = {Constraints[1]:E2}.");
            }
            return;
        }
        static void Main(string[] args)
        {
            Console.WriteLine("Iter) Optimum | f(x,y)+r*H(x,y) | step_X  step_Y | f(x) | g1(x)  g2(x)\n");
            Point a=PenaltyCalculate();
            CheckPoint(a);
            Console.ReadKey();
        }
    }
}
