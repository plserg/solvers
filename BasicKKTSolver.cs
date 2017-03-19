using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;

namespace solvers.solvers
{
   public class BasicKKTSolver : KKTSolver
   {
      public override Vector<double> Solve()
      {

         int n = this.QMatrix.ColumnCount;
         int m = this.AMatrix.RowCount;
         var A = this.AMatrix;
         var AT = this.AMatrix.Transpose();
         var h = this.HVector;
         var g = this.GVector;

         var Qinv = this.QMatrix.Inverse();
         var S = A.Multiply(Qinv.Multiply(AT));//-Schur complement 
         var q = A.Multiply(Qinv.Multiply(g)) - h;
         var y = S.Solve(q);
         var rhs = -g - AT.Multiply(y);
         var x = Qinv.Multiply(rhs);

         var  sol = CreateVector.Dense<double>(n + m);
         for (int i = 0; i < n; ++i)
         {
            sol[i] = x[i];
         }
         for (int i = 0; i < m; ++i)
         {
            sol[n + i] = y[i];
         }
         return sol;
      }
   }
}
