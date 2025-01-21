using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CurveFit
{
    public class LMResult
    {
        public bool Success { get; set; }
        public int Iterations { get; set; }
        public double SumResiduals { get; set; }
        public double[] Residuals { get; set; }
        public double[] Parameters { get; set; }
       
        public LMResult()
        {
            
        }

        public LMResult(bool success, int iterations, double sumResiduals, double[] residuals, double[] parameters)
        {
            Success = success;
            Iterations = iterations;
            SumResiduals = sumResiduals;
            Residuals = residuals;
            Parameters = parameters;
        }
    }
}
