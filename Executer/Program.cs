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
            double[] xData = new double[] { 1, 2, 3, 4, 5, 6 };
            double[] yData = new double[] { 7.3, 3.5, 3.0, 2.1, 2.5, 1.6 };

            LMResult result = LevenbergMarquardt.Fit(LinearFunction, 2, xData, yData);

        }

        private static double LinearFunction(double x, double[] parameters)
        {
            return parameters[0] * x + parameters[1];
        }
    }
}
