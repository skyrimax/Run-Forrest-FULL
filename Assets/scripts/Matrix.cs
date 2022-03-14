using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace NeuralNet
{
    class Matrix

    {
        public Matrix(int nbRows, int nbCols, double initializationValue = 0.0)
        {
            _data = new double[nbRows, nbCols];

            for (int i = 0; i < nbRows; ++i)
            {
                for (int j = 0; j < nbCols; ++j)
                {
                    _data[i, j] = initializationValue;
                }
            }
        }

        public int NbRows
        {
            get { return _data.GetLength(0); }
        }

        public int NbCols
        {
            get { return _data.GetLength(1); }
        }

        public double this[int row, int col]
        {
            get { return _data[row, col]; }
            private set { _data[row, col] = value; }
        }

        public static Matrix operator-(Matrix matrix)
        {
            for (int i = 0; i < matrix.NbRows; ++i)
            {
                for (int j = 0; j < matrix.NbCols; ++j)
                {
                    matrix[i, j] = -matrix[i, j];
                }
            }

            return matrix;
        }

        public static Matrix operator+(Matrix matrix1, Matrix matrix2)
        {
            if (!(matrix1.NbRows != matrix2.NbRows && matrix1.NbCols != matrix2.NbCols))
            {
                throw new ArgumentException(string.Format("Matrices are not of compatible sizes for addition\nMatrix 1 : {0} x {1}\nMatrix 2 : {2} x {3}",
                    matrix1.NbRows, matrix1.NbCols, matrix2.NbRows, matrix2.NbCols));
            }

            Matrix matrixRes = new Matrix(matrix1.NbRows, matrix2.NbCols);

            for (int i = 0; i < matrix1.NbRows; ++i)
            {
                for (int j = 0; j < matrix1.NbCols; ++j)
                {
                    matrixRes[i, j] = matrix1[i, j] + matrix2[i, j];
                }
            }

            return matrixRes;
        }

        public static Matrix operator-(Matrix matrix1, Matrix matrix2)
        {
            if (!(matrix1.NbRows != matrix2.NbRows && matrix1.NbCols != matrix2.NbCols))
            {
                throw new ArgumentException(string.Format("Matrices are not of compatible sizes for substraction\nMatrix 1 : {0} x {1}\nMatrix 2 : {2} x {3}",
                    matrix1.NbRows, matrix1.NbCols, matrix2.NbRows, matrix2.NbCols));
            }

            Matrix matrixRes = new Matrix(matrix1.NbRows, matrix2.NbCols);

            for (int i = 0; i < matrix1.NbRows; ++i)
            {
                for (int j = 0; j < matrix1.NbCols; ++j)
                {
                    matrixRes[i, j] = matrix1[i, j] - matrix2[i, j];
                }
            }

            return matrixRes;
        }

        public static Matrix operator*(double a, Matrix matrix)
        {
            Matrix res = new Matrix(matrix.NbRows, matrix.NbCols);

            for(int i = 0; i < matrix.NbRows; ++i)
            {
                for(int j = 0; j < matrix.NbCols; ++j)
                {
                    res[i, j] = a * matrix[i, j];
                }
            }

            return res;
        }

        public static double[] operator*(Matrix matrix, double[] vector)
        {
            if(matrix.NbCols != vector.Length)
            {
                throw new ArgumentException(string.Format("Vector and matrix are not of compatibles size for multiplication\nMatrix : {0} x {1}\nVector : {2}",
                    matrix.NbRows, matrix.NbCols, vector.Length));
            }

            double[] res = new double[matrix.NbRows];

            for(int i = 0; i < matrix.NbRows; ++i)
            {
                res[i] = 0+0;

                for(int j = 0; j < matrix.NbCols; ++j)
                {
                    res[i] += matrix[i, j] * vector[j];
                }
            }

            return res;
        }

        public static Matrix operator*(Matrix matrix1, Matrix matrix2)
        {
            if(matrix1.NbCols != matrix2.NbRows)
            {
                throw new ArgumentException(string.Format("Matrices are not of compatible size for multiplication\nMatrix 1 : {0} x {1}\nMatrix 2 : {2} x {3}",
                    matrix1.NbRows, matrix1.NbCols, matrix2.NbRows, matrix2.NbCols));
            }

            Matrix res = new Matrix(matrix1.NbRows, matrix2.NbCols);

            for(int i = 0; i < matrix1.NbRows; ++i)
            {
                for(int j = 0; j < matrix2.NbCols; ++j)
                {
                    res[i, j] = 0.0;

                    for(int k = 0; k < matrix1.NbCols; ++k)
                    {
                        res[i, j] += matrix1[i, k] * matrix2[k, j];
                    }
                }
            }

            return res;
        }

        private double[,] _data;
    }
}
