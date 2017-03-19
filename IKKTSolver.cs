using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MathNet.Numerics.LinearAlgebra;

namespace QP.solvers
{
   public interface IKKTSolver
   {
     Vector<double>  Solve();
     Vector<double> residual(Vector<double> x);
     Matrix<double> AMatrix { get; set; }
     Matrix<double> QMatrix { get; set; }
     Vector<double> GVector { get; set; }
     Vector<double> HVector { get; set; }
     double Tolerance { get; set; }
   }
}
