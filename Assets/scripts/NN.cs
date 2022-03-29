using System;
using System.Collections;
using System.Collections.Generic;

namespace NeuralNet
{
    public class NN
    {
        string _name = "Neural Network";
        public string Name { get { return "Forrest " + _name; } }

        public int NbInputs { get { return _nbInputs; } }
        int _nbInputs;

        public int NbHiddenLayers { get { return _nbHiddenLayers; } }
        int _nbHiddenLayers;

        public int[] NbNodesHiddenLayers { get { return _nbNodesHiddenLayers; } }
        int[] _nbNodesHiddenLayers;

        Matrix[] _hiddenLayersWeights;
        double[][] _hiddenLayersBiases;

        public int WeightCount { get { return TotalWeightCount(); } }

        public double Fitness { get { return _fitness; } }
        double _fitness;

        //float[] _hL;
        //public int hL { get { return _hL.Length; } }
        //int _inputs;
        //public int inputs { get { return _inputs; } }
        //float[][] _hLw;
        //float[] _hLb;
        //float[] Ow;
        //float Ob;

        /// <summary>
        /// Sometimes you just want access to variables like fitness & don't need to pass in node counts
        /// </summary>
        public NN()
        {

        }

        /// <summary>
        /// In most cases its important to assign input & Hidden Layer node count
        /// </summary>
        /// <param name="inpsNumb"></param>
        /// <param name="_hLNodes"></param>
        //public NN(int inpsNumb, int _hLNodes)
        //{
        //    _hL = new float[_hLNodes];
        //    IniWeights(inpsNumb);
        //    _inputs = inpsNumb;
        //}
        public NN(int nbInputs, int[] nbNodesHiddenLayers)
        {
            // Comminting parameters to memory for futur use
            _nbInputs = nbInputs;
            _nbHiddenLayers = nbNodesHiddenLayers.Length;
            _nbNodesHiddenLayers = nbNodesHiddenLayers;

            // Creating the arrays containing the weights and biases of each hidden layer plus output layer
            _hiddenLayersWeights = new Matrix[_nbHiddenLayers + 1];
            _hiddenLayersBiases = new double[_nbHiddenLayers + 1][];

            for(int i = 0; i <= _nbHiddenLayers; ++i)
            {
                // Creating weight matrix for the hidden layer i
                // Size is number of nodes of current layer x number of nodes in previous layer
                _hiddenLayersWeights[i] = new Matrix(
                    i == _nbHiddenLayers ? 1 : _nbNodesHiddenLayers[i],
                    i == 0 ? _nbInputs : _nbNodesHiddenLayers[i - 1]);

                // Creating bias vector for hidden layer i
                _hiddenLayersBiases[i] = new double[i == _nbHiddenLayers ? 1 :_nbNodesHiddenLayers[i]];
            }
        }

        /// <summary>
        /// Add a fitness value to this network
        /// </summary>
        /// <param name="fit"></param>
        public void SetFitness(double fit)
        {
            _fitness = fit;
        }

        /// <summary>
        /// Inputing an interager (number of _inputs) will initialize all weights with zeros
        /// </summary>
        /// <param name="inps">number of _inputs</param>
        //void IniWeights(int inps)
        //{
        //    // 
        //    _inputs = inps;

        //    // 
        //    Ow = new float[_hL.Length];

        //    // Set a double nested array, [hidden layer length][input length]
        //    _hLw = new float[_hL.Length][];

        //    // Set the array for biases that are = to the length of nodes
        //    _hLb = new float[_hLw.Length];

        //    // Next we loop through each input array nest & set the length of each to the number of _inputs
        //    for (int i = 0; i < _hLw.GetLength(0); i++)
        //    {
        //        _hLw[i] = new float[inps];
        //        _hLb[i] = 0;// r.Next(-cap, cap);
        //        Ow[i] = 0;// r.Next(-cap, cap);

        //        // 
        //        for (int j = 0; j < _hLw[i].Length; j++)
        //        {
        //            _hLw[i][j] = 0;// r.Next(-cap, cap);
        //        }
        //    }
        //}

        /// <summary>
        /// This function is used for loading weight matricies (brains/DNA/whatever you want to call it) Inputting a float array will set the weights to this array
        /// </summary>
        /// <param name="w"></param>
        //public void IniWeights(float[] w)
        //{
        //    // Okay don't freak out from this code let's walk through it step by step.

        //    // first were going to run a for loop the length of the number of nodes in our Hidden Layer
        //    for (int i = 0; i < _hLw.Length; i++)
        //    {
        //        // then we're going to map all of the hidden layer weights from our float vector input w
        //        // & because I encode the input vector w as all _hL[0] weights then _hL[0] bias, then _hL[1] weight -> _hL[1] bias
        //        // all the way to it's length, then do all Output weights, then Output bias, I am going to map my Neural Net accordingly 
        //        for (int j = 0; j < _hLw[i].Length; j++)
        //        {
        //            // this is exactly what I described above but in math form
        //            _hLw[i][j] = w[j + (i * 6)];
        //        }

        //        // then we're going to map all the biases for ALL nodes in the Hidden Layer, as described above but in math form
        //        _hLb[i] = w[5 + (i * 6)];

        //        // Then we're going to map all of the Output weights, as described above but in math form
        //        Ow[i] = w[i + 24];

        //        // & lastly the Output bias is the last on the input w sequence & we're done!
        //        Ob = w[w.Length - 1];
        //    }
        //}
        public void IniWeights(double[] dna)
        {
            // Making sure the DNA sequence is the right size for the neural network
            if(dna.Length != WeightCount)
            {
                throw new ArgumentException(string.Format("DNA is not the right size\nDNA lenght : {0}\nNb of weights and biases : {1}", dna.Length, WeightCount));
            }

            // Interator used to read throught the DNA sequence
            int dnaIterator = 0;

            // First for loops throught all hidden layers plus output layer
            for(int i = 0; i < _hiddenLayersWeights.Length; ++i)
            {
                // These 2 for loops fill the matrix of the layer i, filling row by row from left to right
                for(int j = 0; j < _hiddenLayersWeights[i].NbRows; ++j)
                {
                    for(int k = 0; k < _hiddenLayersWeights[i].NbCols; ++k)
                    {
                        // Assign weight of current position in DNA sequence to matrix element [j, k] of layer [i]
                        _hiddenLayersWeights[i][j, k] = dna[dnaIterator++];
                    }

                    // Assigns the biase for row j of the layer i
                    _hiddenLayersBiases[i][j] = dna[dnaIterator++];
                }
            }
        }

        /// <summary>
        /// This function is used for loading weight matricies (brains/DNA/whatever you want to call it) Inputting a string with values separeted with ","s will convert the string into a float array & use that as weights
        /// </summary>
        /// <param name="w"></param>
        //public void IniWeights(string inp)
        //{
        //    float[] w = new float[WeightCount];

        //    for (int i = 0; i < w.Length; i++)
        //    {
        //        w[i] = float.Parse(inp.Split(',')[i]);
        //    }

        //    // Okay don't freak out from this code let's walk through it step by step.

        //    // first were going to run a for loop the length of the number of nodes in our Hidden Layer
        //    for (int i = 0; i < _hLw.Length; i++)
        //    {
        //        // then we're going to map all of the hidden layer weights from our float vector input w
        //        // & because I encode the input vector w as all _hL[0] weights then _hL[0] bias, then _hL[1] weight -> _hL[1] bias
        //        // all the way to it's length, then do all Output weights, then Output bias, I am going to map my Neural Net accordingly 
        //        for (int j = 0; j < _hLw[i].Length; j++)
        //        {
        //            // this is exactly what I described above but in math form
        //            _hLw[i][j] = w[j + (i * 6)];
        //        }

        //        // then we're going to map all the biases for ALL nodes in the Hidden Layer, as described above but in math form
        //        _hLb[i] = w[5 + (i * 6)];

        //        // Then we're going to map all of the Output weights, as described above but in math form
        //        Ow[i] = w[i + 24];

        //        // & lastly the Output bias is the last on the input w sequence & we're done!
        //        Ob = w[w.Length - 1];
        //    }
        //}
        public void IniWeights(string dna)
        {
            // Splits DNA sequence in an array of genes
            string[] splitedDNA = dna.Split(',');

            // Making sure the DNA sequence is the right size for the neural network
            if (splitedDNA.Length != WeightCount)
            {
                throw new ArgumentException(string.Format("DNA is not the right size\nDNA lenght : {0}\nNb of weights and biases : {1}", dna.Length, WeightCount));
            }

            // Array of number converted genes
            double[] numericSplitedDNA = new double[splitedDNA.Length];

            // Converts each gene into numbers
            for(int i = 0; i < splitedDNA.Length; ++i)
            {
                numericSplitedDNA[i] = double.Parse(splitedDNA[i]);
            }

            // Initializes weights based on array of numeric genes
            IniWeights(numericSplitedDNA);
        }

        /// <summary>
        /// This will return the output aka guess for the joystick movement
        /// </summary>
        /// <param name="inps"></param>
        /// <returns></returns>
        //public float CalculateNN(float[] inps)
        //{
        //    // Set the value of all hidden layer outputs
        //    for (int i = 0; i < _hL.Length; i++)
        //    {
        //        _hL[i] = ReLU(Sum(inps, _hLw[i]) + _hLb[i]);
        //    }

        //    float O = SoftSign(Sum(Ow, _hL) + Ob);
        //    //float O = ACTIVATION(Sum(Ow, _hL) + Ob);

        //    return O;
        //}
        public double CalculateNN(double[] inputs)
        {
            if(inputs.Length != _nbInputs)
            {
                throw new ArgumentException(string.Format("Number of input provided does not match with expectation from neural network\n" +
                    "Number of inputs provided : {0}\nNumber of inputs expected : {1}", inputs.Length, _nbInputs));
            }

            for(int i = 0; i < _nbHiddenLayers; ++i)
            {
                inputs = ReLU(Sum(_hiddenLayersWeights[i] * inputs, _hiddenLayersBiases[i]));
            }

            return SoftSign(Sum(_hiddenLayersWeights[_nbHiddenLayers] * inputs, _hiddenLayersBiases[_nbHiddenLayers])[0]);
        }

        /// <summary>
        /// The summation of a[0-n] * b[0-n]
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        //float Sum(float[] a, float[] b)
        //{
        //    float ret = 0;
        //    for (int i = 0; i < a.Length; i++)
        //    {
        //        ret += a[i] * b[i];
        //    }
        //    return ret;
        //}

        /// <summary>
        /// The element-wise addition of two vectors
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        static double[] Sum(double[] a, double[] b)
        {
            if(a.Length != b.Length)
            {
                throw new ArgumentException(string.Format("Arrays are not of compatiple sizes\nArray a : {0}\nArray b : {1}", a.Length, b.Length));
            }

            double[] res = new double[a.Length];

            for(int i = 0; i < a.Length; ++i)
            {
                res[i] = a[i] + b[i];
            }

            return res;
        }

        /// <summary>
        /// ReLU activation function returns the max between 0 & your inp
        /// </summary>
        /// <param name="inp"></param>
        /// <returns></returns>
        //float ReLU(float inp)
        //{
        //    return Math.Max(0, inp);
        //}
        double ReLU(double x)
        {
            return Math.Max(0, x);
        }

        /// <summary>
        /// Applies ReLU activation function to each element of array
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        double[] ReLU(double[] x)
        {
            for(int i = 0; i < x.Length; ++i)
            {
                x[i] = ReLU(x[i]);
            }

            return x;
        }

        /// <summary>
        /// Softsign is an activation function that maps your input in between -1 & 1
        /// </summary>
        //float SoftSign(float inp)
        //{
        //    return inp / (1 + Math.Abs(inp));
        //}
        double SoftSign(double x)
        {
            return x / (1 + Math.Abs(x));
        }


        /// <summary>
        /// This returns the Genetic Code Sequence for the attempt in an easy readable string format
        /// </summary>
        /// <returns></returns>
        //public string ReadBrain()
        //{
        //    string dna = "";

        //    for (int i = 0; i < _hLw.Length; i++)
        //    {
        //        for (int j = 0; j < _hLw[i].Length; j++)
        //        {
        //            dna += _hLw[i][j] + ",";
        //        }
        //        dna += _hLb[i] + ",";
        //    }

        //    for (int i = 0; i < Ow.Length; i++)
        //    {
        //        dna += Ow[i] + ",";
        //    }

        //    dna += Ob;

        //    return dna;
        //}
        public string ReadBrain()
        {
            string dna = "";

            for(int i = 0; i < _nbHiddenLayers + 1; ++i)
            {
                for(int j = 0; j < _hiddenLayersWeights[i].NbRows; ++j)
                {
                    for(int k = 0; k < _hiddenLayersWeights[i].NbCols; ++k)
                    {
                        dna += _hiddenLayersWeights[i][j, k] + ",";
                    }

                    dna += _hiddenLayersBiases[i][j] + ",";
                }
            }

            dna = dna.Remove(dna.Length - 1);

            return dna;
        }

        /// <summary>
        /// This returns the Genetic Code Sequence for the attempt in a float array
        /// </summary>
        /// <returns></returns>
        //public float[] GetBrain()
        //{
        //    string[] dna = ReadBrain().Split(',');

        //    float[] ret = new float[dna.Length];

        //    for (int i = 0; i < dna.Length; i++)
        //    {
        //        ret[i] = float.Parse(dna[i]);
        //    }

        //    return ret;
        //}
        public double[] GetBrain()
        {
            string[] dna = ReadBrain().Split(',');

            double[] ret = new double[dna.Length];

            for (int i = 0; i < dna.Length; i++)
            {
                ret[i] = double.Parse(dna[i]);
            }

            return ret;
        }


        /// <summary>
        /// Returns the total weight count
        /// </summary>
        /// <returns></returns>
        //int TotalWeightCount()
        //{
        //    return (_hLw[0].Length * _inputs) + Ow.Length;
        //}
        int TotalWeightCount()
        {
            int weightCount = 0;

            for(int i = 0; i < _nbHiddenLayers + 1; ++i)
            {
                weightCount += _hiddenLayersWeights[i].NbRows * _hiddenLayersWeights[i].NbCols + _hiddenLayersBiases[i].Length;
            }

            return weightCount;
        }

        /// <summary>
        /// Sets the name of the NN
        /// </summary>
        /// <param name="n"></param>
        public void SetName(string n)
        {
            _name = n;
        }
    }
}