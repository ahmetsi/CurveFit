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
        public delegate double[] FitDerivative(double x, double[] parameters); // Partial derivatives of the function 

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

            double minLambda = 1E-7; // Minimum allowed value for lambda
            double maxLambda = 1E+6; // Maximum allowed value for lambda

            double[] residuals = new double[n];
            double[,] jacobian = new double[n, p];

            double lastError = double.MaxValue;
            int iteration;
            for (iteration = 0; iteration < maxIterations; iteration++)
            {
                // Compute residuals and Jacobian
                ComputeResiduals(model, parameters, xData, yData, residuals);
                ComputeJacobian(derivative, parameters, xData, jacobian);

                double currentError = SumResiduals(residuals);

                // Adjust the damping factor dynamically
                if (currentError < lastError)
                {
                    // Error reduced: Decrease lambda to move towards Gauss-Newton method
                    lambda = Math.Max(lambda / 10, minLambda);
                }
                else
                {
                    // Error increased: Increase lambda to move towards gradient descent
                    lambda = Math.Min(lambda * 10, maxLambda);
                }

                // Check if progress is bad
                if (Math.Abs(lastError - currentError) < 1E-8)
                {
                    result.Message = "Iterations did not make good progress.";
                }

                lastError = currentError;

                // Calculate J^T * J (Jacobian transposed times Jacobian) and J^T * r (Jacobian transposed times residuals)
                double[,] jTj = new double[p, p];
                double[] jTr = new double[p];
                for (int i = 0; i < n; i++)
                {
                    for (int row = 0; row < p; row++)
                    {
                        for (int col = 0; col < p; col++)
                            jTj[row, col] += jacobian[i, row] * jacobian[i, col];

                        jTr[row] += jacobian[i, row] * residuals[i];
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
                    parameters[i] -= dp[i];
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
                    result.Message = "Solution successfully converged.";
                    break;
                }
            }

            if (string.IsNullOrEmpty(result.Message))
            {
                result.Message = "Maximum number of iterations have been reached.";
            }

            result.Success = success;
            result.Iterations = iteration;
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
                sum += r * r;
            }
            return sum;
        }

        private static FitDerivative CreateDerivative(FitFunction model, double eps)
        {
            return (x, parameters) =>
            {
                double[] parametersCloned = (double[])parameters.Clone();
                double y1 = model(x, parameters);

                double[] derivatives = new double[parameters.Length];
                for (int i = 0; i < parameters.Length; i++)
                {
                    parametersCloned = (double[])parameters.Clone();
                    parametersCloned[i] += eps;
                    double y2 = model(x, parametersCloned);
                    derivatives[i] = (y2 - y1) / eps;
                }
                return derivatives;
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
            for (int i = 0; i < xData.Length; i++)
            {
                double x = xData[i];
                double[] derivatives = derivative(x, parameters);
                for (int j = 0; j < parameters.Length; j++)
                {
                    jacobian[i, j] = derivatives[j];
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
