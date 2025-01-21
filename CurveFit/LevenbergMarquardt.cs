using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CurveFit
{
    public static class LevenbergMarquardt
    {
        public delegate double FitFunction(double x, double[] parameters); // Function to fit
        public delegate double FitDerivative(double x, double[] parameters, int i); // Derivative of function wrt i'th parameter

        public static LMResult Fit(FitFunction model, int numParameters, double[] xData, double[] yData,
            double[] initialParameters = null, FitDerivative derivative = null, int maxIterations = 0,
            double tolerance = 1E-8, double lambda = 0.001, double eps = 1E-8)
        {   
            double[] parameters;
            if (initialParameters == null)
            {
                parameters = Enumerable.Repeat(1.0, numParameters).ToArray();
            }
            else
            {
                parameters = (double[])initialParameters.Clone();
            }

            if (maxIterations == 0)
            {
                if (derivative == null)
                {
                    maxIterations = 200 * (numParameters + 1);
                }
                else
                {
                    maxIterations = 100 * (numParameters + 1);
                }
            }

            if (derivative == null)
            {
                derivative = CreateDerivative(model, eps);
            }

            

            LMResult result = new LMResult();
            bool success = false;

            int n = xData.Length;
            int p = numParameters;

            double[] residuals = new double[n];
            double[,] jacobian = new double[n, p];

            int iteration;
            for (iteration = 0; iteration < maxIterations; iteration++)
            {
                // Compute residuals and Jacobian
                ComputeResiduals(model, parameters, xData, yData, residuals);
                ComputeJacobian(derivative, parameters, xData, jacobian);

                // Calculate J^T * J (Jacobian transposed times Jacobian) and J^T * r (Jacobian transposed times residuals)
                double[,] jTj = new double[p, p];
                double[] jTr = new double[p];
                for (int i = 0; i < n; i++)
                {
                    for (int row = 0; row < p; row++)
                    {
                        for (int col = 0; col < p; col++)
                            jTj[row, col] += jacobian[i, row] * jacobian[i, col];

                        jTr[row] -= jacobian[i, row] * residuals[i];
                    }
                }

                // Add damping factor to diagonal
                for (int i = 0; i < p; i++)
                {
                    jTj[i, i] += lambda;
                }
                    
                // Solve for parameter update using J^T * J * dp = -J^T * r
                double[] dp = SolveLinearSystem(jTj, jTr);

                // Update parameters
                for (int i = 0; i < p; i++)
                {
                    parameters[i] += dp[i];
                }
                    
                // Check for convergence
                double maxUpdate = 0;
                foreach (double d in dp)
                {
                    maxUpdate = Math.Max(maxUpdate, Math.Abs(d));
                }

                if (maxUpdate < tolerance)
                {
                    success = true;
                    break;
                }
            }

            result.Success = success;
            result.Iterations = iteration + 1;
            result.SumResiduals = SumResiduals(residuals);
            result.Residuals = residuals;
            result.Parameters = parameters;

            return result;
        }

        private static double SumResiduals(double[] residuals)
        {
            double sum = 0;
            foreach (double r in residuals)
            {
                sum += Math.Pow(r, 2);
            }
            return sum;
        }

        private static FitDerivative CreateDerivative(FitFunction model, double eps)
        {
            return (x, parameters, i) =>
            {
                double[] parametersCloned = (double[])parameters.Clone();
                double y1 = model(x, parametersCloned);
                parametersCloned[i] += eps;
                double y2 = model(x, parametersCloned);
                return (y2 - y1) / eps;
            };
        }

        private static void ComputeResiduals(FitFunction model, double[] parameters, double[] xData, double[] yData, double[] residuals)
        {
            for (int i = 0; i < xData.Length; i++)
            {
                double x = xData[i];
                double yPredicted = model(x, parameters);
                residuals[i] = yPredicted - yData[i];
            }
        }

        private static void ComputeJacobian(FitDerivative derivative, double[] parameters, double[] xData, double[,] jacobian)
        {
            for (int i = 0; i < parameters.Length; i++)
            {
                for (int j = 0; j < xData.Length; j++)
                {
                    double x = xData[j];
                    jacobian[j, i] = derivative(x, parameters, i);
                }
            }
        }

        private static double[] SolveLinearSystem(double[,] A, double[] b)
        {
            int n = b.Length;
            double[] x = new double[n];
            double[,] lu = (double[,])A.Clone();

            // LU Decomposition
            for (int k = 0; k < n; k++)
            {
                for (int i = k + 1; i < n; i++)
                {
                    lu[i, k] /= lu[k, k];
                    for (int j = k + 1; j < n; j++)
                        lu[i, j] -= lu[i, k] * lu[k, j];
                }
            }

            // Forward substitution
            double[] y = new double[n];
            for (int i = 0; i < n; i++)
            {
                y[i] = b[i];
                for (int j = 0; j < i; j++)
                    y[i] -= lu[i, j] * y[j];
            }

            // Backward substitution
            for (int i = n - 1; i >= 0; i--)
            {
                x[i] = y[i];
                for (int j = i + 1; j < n; j++)
                    x[i] -= lu[i, j] * x[j];
                x[i] /= lu[i, i];
            }

            return x;
        }
    }
}
