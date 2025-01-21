using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CurveFit;

namespace Executer
{
    public class Program
    {
        static void Main(string[] args)
        {
            double[] xData = new double[] { 1, 2, 3, 4, 5, 6, 7 };
            double[] yData = new double[] { 7.3, 3.5, 3.5, 3, 2.1, 2.5, 3 };

            LMResult result = LevenbergMarquardt.Fit(PolyFunction, 3, xData, yData, derivative: PolyDerivative);

        }

        private static double PolyFunction(double x, double[] parameters)
        {
            return parameters[0] * x * x + parameters[1] * x + parameters[2];
        }

        private static double[] PolyDerivative(double x, double[] parameters)
        {
            double[] derivatives = new double[parameters.Length];
            derivatives[0] = x * x;
            derivatives[1] = x;
            derivatives[2] = 1;
            return derivatives;
        }
    }
}
