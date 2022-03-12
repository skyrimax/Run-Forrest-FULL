using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace NeuralNet
{
    class Matrix<T> where T : struct

    {
        public Matrix(int nbRows, int nbCols, T initializationValue = default(T))
        {
            _data = new T[nbRows, nbCols];

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

        public T this[int row, int col]
        {
            get { return _data[row, col]; }
            private set { _data[row, col] = value; }
        }

        public static Matrix<T> operator-(Matrix<T> matrix)
        {
            for (int i = 0; i < matrix.NbRows; ++i)
            {
                for (int j = 0; j < matrix.NbCols; ++j)
                {
                    matrix[i, j] = Inverse(matrix[i, j]);
                }
            }

            return matrix;
        }

        public static Matrix<T> operator+(Matrix<T> matrix1, Matrix<T> matrix2)
        {
            if (!(matrix1.NbRows != matrix2.NbRows && matrix1.NbCols != matrix2.NbCols))
            {
                throw new ArgumentException(string.Format("Matrices are not of compatible sizes for addition\nMatrix 1 : {0} x {1}\nMatrix 2 : {2} x {3}", matrix1.NbRows, matrix1.NbCols, matrix2.NbRows, matrix2.NbCols));
            }

            Matrix<T> matrixRes = new Matrix<T>(matrix1.NbRows, matrix2.NbCols);

            for (int i = 0; i < matrix1.NbRows; ++i)
            {
                for (int j = 0; j < matrix1.NbCols; ++j)
                {
                    matrixRes[i, j] = Add(matrix1[i, j], matrix2[i, j]);
                }
            }

            return matrixRes;
        }

        public static Matrix<T> operator-(Matrix<T> matrix1, Matrix<T> matrix2)
        {
            if (!(matrix1.NbRows != matrix2.NbRows && matrix1.NbCols != matrix2.NbCols))
            {
                throw new ArgumentException(string.Format("Matrices are not of compatible sizes for substraction\nMatrix 1 : {0} x {1}\nMatrix 2 : {2} x {3}", matrix1.NbRows, matrix1.NbCols, matrix2.NbRows, matrix2.NbCols));
            }

            Matrix<T> matrixRes = new Matrix<T>(matrix1.NbRows, matrix2.NbCols);

            for (int i = 0; i < matrix1.NbRows; ++i)
            {
                for (int j = 0; j < matrix1.NbCols; ++j)
                {
                    matrixRes[i, j] = Substract(matrix1[i, j], matrix2[i, j]);
                }
            }

            return matrixRes;
        }

        public static Matrix<T> operator*(T a, Matrix<T> matrix)
        {
            Matrix<T> res = new Matrix<T>(matrix.NbRows, matrix.NbCols);

            for(int i = 0; i < matrix.NbRows; ++i)
            {
                for(int j = 0; j < matrix.NbCols; ++j)
                {
                    res[i, j] = Multiply(a, matrix[i, j]);
                }
            }

            return res;
        }

        public static T[] operator*(Matrix<T> matrix, T[] vector)
        {
            if(matrix.NbCols != vector.Length)
            {
                throw new ArgumentException(string.Format("Vector and matrix are not of compatibles size for multiplication\nMatrix : {0} x {1}\nVector : {2}", matrix.NbRows, matrix.NbCols, vector.Length));
            }

            T[] res = new T[matrix.NbRows];

            for(int i =0; i < matrix.NbRows; ++i)
            {
                res[i] = default(T);

                for(int j = 0; j < matrix.NbCols; ++j)
                {
                    res[i] = Add(res[i], Multiply(matrix[i, j], vector[j]));
                }
            }

            return res;
        }

        public static Matrix<T> operator*(Matrix<T> matrix1, Matrix<T> matrix2)
        {
            if(matrix1.NbCols != matrix2.NbRows)
            {
                throw new ArgumentException(string.Format("Matrices are not of compatible size for multiplication\nMatrix 1 : {0} x {1}\nMatrix 2 : {2} x {3}", matrix1.NbRows, matrix1.NbCols, matrix2.NbRows, matrix2.NbCols));
            }

            Matrix<T> res = new Matrix<T>(matrix1.NbRows, matrix2.NbCols);

            for(int i = 0; i < matrix1.NbRows; ++i)
            {
                for(int j = 0; j < matrix2.NbCols; ++j)
                {
                    res[i, j] = default(T);

                    for(int k = 0; k < matrix1.NbCols; ++k)
                    {
                        res[i, j] = Add(res[i, j], Multiply(matrix1[i, k], matrix2[k, j]));
                    }
                }
            }

            return res;
        }

        private static T Inverse(T x)
        {
            dynamic dx = x;
            return -dx;
        }

        private static T Add<T>(T x, T y)
        {
            dynamic dx = x, dy = y;
            return dx + dy;
        }

        private static T Substract<T>(T x, T y)
        {
            dynamic dx = x, dy = y;
            return dx - dy;
        }

        private static T Multiply<T>(T x, T y)
        {
            dynamic dx = x, dy = y;
            return dx * dy;
        }

        private T[,] _data;
    }
}
