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
            double[] xData = new double[] { 0.1, 1.2, 2.1, 2.9, 4.2, 5.2 };
            double[] yData = new double[] { 7.3, 3.5, 3.0, 2.1, 2.5, 1.6 };

            LMResult result = LevenbergMarquardt.Fit(ExponentialFunction, 3, xData, yData);

        }

        private static double ExponentialFunction(double x, double[] parameters)
        {
            return parameters[0] * Math.Exp(-parameters[1] * x) + parameters[2];
        }
    }
}
