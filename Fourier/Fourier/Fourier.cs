using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Drawing;
using System.Drawing.Imaging;

namespace Fourier
{
    public class Fourier
    {

        public enum Axis { x, y, u, v }
        public struct number
        {
            public double real, image;
            public number(double real, double image) { this.real = real; this.image = image; }
            public void FromAbsAngle(double abs, double angle)
            {
                this.real = (double)(Math.Cos(angle) * abs);
                this.image = (double)(Math.Sin(angle) * abs);
            }
            public void Add(double real, double image) { this.real += real; this.image += image; }
            public void Add(number Num) { this.real += Num.real; this.image += Num.image; }
            public void AddAbsAngle(double abs, double angle)
            {
                this.real += (double)(Math.Cos(angle) * abs);
                this.image += (double)(Math.Sin(angle) * abs);
            }
            public double abs()
            {
                double abs;
                abs = (double)Math.Sqrt(Math.Pow(this.real,2.0)+Math.Pow(this.image,2.0));
                return abs;
            }
            public double angle()
            {
                double angle;
                angle = (double)Math.Atan2(image, real);
                return angle;
            }
            public override string ToString()
            {
                string str; if (image >= 0) str = " +"; else str = " ";
                return (real.ToString() + str + image.ToString() + "j");
            }

        }
        public number Multiplication(number num1, number num2)
        {
            number result = new number();
            result.real = num1.real * num2.real - num1.image * num2.image;
            result.image = num1.real * num2.image + num2.real * num1.image;
            return result;
        }
        public number Summation(number num1, number num2)
        {
            number result = new number();
            result.real = num1.real + num2.real;
            result.image = num1.image + num2.image;
            return result;
        }
        public number Subtraction(number num1, number num2)
        {
            number result = new number();
            result.real = num1.real - num2.real;
            result.image = num1.image - num2.image;
            return result;
        }
        public number Division(number num1, number num2)
        {
            number result = new number();
            result.real = (num1.real * num2.real + num1.image * num2.image) / (num2.real * num2.real + num2.image * num2.image);
            result.image = (num1.image * num2.real - num1.real * num2.image) / (num2.real * num2.real + num2.image * num2.image);
            return result;
        }
        public number Conjucate(number num)
        {
            number result = new number();
            result.real = num.real;
            result.image = -num.image;
            return result;
        }


        public number[] DFT1(number[] ft)
        {
            number buf = new number();
            number[] Ff = new number[ft.Count()];
            for (int f = 0; f < ft.Count(); f++)
            {
                Ff[f] = new number(0, 0);
                for (int t = 0; t < ft.Count(); t++)
                {
                    buf.FromAbsAngle(ft[t].abs(), ft[t].angle() - 2 * Math.PI * f * t / ft.Count());
                    Ff[f].Add(buf);
                }
            }
            return Ff;
        }
        public void DFT1(ref number[,] fxy, ref number[,] Fuv, Axis axis, int Row_Column)
        {
            number buf = new number();
            switch (axis)
            {
                case Axis.x:
                    for (int u = 0; u < fxy.GetLength(0); u++)
                    {
                        Fuv[u, Row_Column] = new number(0, 0);
                        for (int x = 0; x < fxy.GetLength(0); x++)
                        {
                            buf.FromAbsAngle(fxy[x, Row_Column].abs(), fxy[x, Row_Column].angle() - 2 * Math.PI * x * u / fxy.GetLength(0));
                            Fuv[u, Row_Column].Add(buf);
                        }
                    }
                    break;
                case Axis.y:
                    for (int v = 0; v < fxy.GetLength(1); v++)
                    {
                        Fuv[Row_Column, v] = new number(0, 0);
                        for (int y = 0; y < fxy.GetLength(1); y++)
                        {
                            buf.FromAbsAngle(fxy[Row_Column, y].abs(), fxy[Row_Column, y].angle() - 2 * Math.PI * y * v / fxy.GetLength(1));
                            Fuv[Row_Column, v].Add(buf);
                        }
                    }
                    break;
            }
        }
        public number[,] DFT2(number[,] fxy)
        {
            number[,] buffer = new number[fxy.GetLength(0), fxy.GetLength(1)];
            number[,] Fuv = new number[fxy.GetLength(0), fxy.GetLength(1)];
            for (int x = 0; x < fxy.GetLength(0); x++)
            {
                DFT1(ref fxy, ref buffer, Axis.y, x);
            }
            for (int y = 0; y < fxy.GetLength(1); y++)
            {
                DFT1(ref buffer, ref Fuv, Axis.x, y);
            }
            return Fuv;
        }
        public number[] iDFT1(number[] Ff)
        {
            number buf = new number();
            number[] ft = new number[Ff.Count()];
            for (int t = 0; t < Ff.Count(); t++)
            {
                ft[t] = new number(0, 0);
                for (int f = 0; f < Ff.Count(); f++)
                {
                    buf.FromAbsAngle(Ff[f].abs() / Ff.Count(), Ff[f].angle() + 2 * Math.PI * f * t / Ff.Count());
                    ft[t].Add(buf);
                }
            }
            return ft;
        }
        public void iDFT1(ref number[,] Fuv, ref number[,] fxy, Axis axis, int Row_Column)
        {
            number buf = new number();
            switch (axis)
            {
                case Axis.u:
                    for (int x = 0; x < Fuv.GetLength(0); x++)
                    {
                        fxy[x, Row_Column] = new number(0, 0);
                        for (int u = 0; u < Fuv.GetLength(0); u++)
                        {
                            buf.FromAbsAngle(Fuv[u, Row_Column].abs() / Fuv.GetLength(0), Fuv[u, Row_Column].angle() + 2 * Math.PI * u * x / Fuv.GetLength(0));
                            fxy[x, Row_Column].Add(buf);
                        }
                    }
                    break;
                case Axis.v:
                    for (int y = 0; y < Fuv.GetLength(1); y++)
                    {
                        fxy[Row_Column, y] = new number(0, 0);
                        for (int v = 0; v < Fuv.GetLength(1); v++)
                        {
                            buf.FromAbsAngle(Fuv[Row_Column, v].abs() / Fuv.GetLength(1), Fuv[Row_Column, v].angle() + 2 * Math.PI * v * y / Fuv.GetLength(1));
                            fxy[Row_Column, y].Add(buf);
                        }
                    }
                    break;
            }
        }
        public number[,] iDFT2(number[,] Fuv)
        {
            number[,] buffer = new number[Fuv.GetLength(0), Fuv.GetLength(1)];
            number[,] fxy = new number[Fuv.GetLength(0), Fuv.GetLength(1)];
            for (int u = 0; u < Fuv.GetLength(0); u++)
            {
                iDFT1(ref Fuv, ref buffer, Axis.v, u);
            }
            for (int v = 0; v < Fuv.GetLength(1); v++)
            {
                iDFT1(ref buffer, ref fxy, Axis.u, v);
            }
            return fxy;
        }


        public number[] FFT1(number[] ft)
        {
            number[] Ff = new number[ft.Length];
            int N = ft.Length;
            if (N == 1) return ft;
            number Wn = new number();
            Wn.FromAbsAngle(1, -2 * Math.PI / N);
            number W = new number(1, 0);
            number[] A0 = new number[N / 2];
            number[] A1 = new number[N / 2];
            for (int n = 0; n < N/2; n++)
            {
                A0[n] = ft[2 * n];
                A1[n] = ft[2 * n + 1];
            }
            number[] Y0 = FFT1(A0);
            number[] Y1 = FFT1(A1);
            for (int k = 0; k < N / 2; k++)
            {
                Ff[k] = Summation(Y0[k], Multiplication(W, Y1[k]));
                Ff[k+N/2] = Subtraction(Y0[k], Multiplication(W, Y1[k]));
                W = Multiplication(W, Wn);
            }
            return Ff;
        }
        public number[] iFFT1(number[] Ff)
        {
            number[] ft = new number[Ff.Length];
            int N = Ff.Length;
            if (N == 1) return Ff;
            number Wn = new number();
            Wn.FromAbsAngle(1 / Ff.Length, 2 * Math.PI / N);
            number W = new number(1, 0);
            number[] A0 = new number[N / 2];
            number[] A1 = new number[N / 2];
            for (int n = 0; n < N / 2; n++)
            {
                A0[n] = Ff[2 * n];
                A1[n] = Ff[2 * n + 1];
            }
            number[] Y0 = iFFT1(A0);
            number[] Y1 = iFFT1(A1);
            for (int k = 0; k < N / 2; k++)
            {
                ft[k] = Summation(Y0[k], Multiplication(W, Y1[k]));
                ft[k + N / 2] = Subtraction(Y0[k], Multiplication(W, Y1[k]));
                W = Multiplication(W, Wn);
            }
            return ft;
        }

        public number[,] FFT2(number[,] fxy)
        {
            number[,] buf = new number[fxy.GetLength(0), fxy.GetLength(1)];
            for (int x = 0; x < fxy.GetLength(0); x++)
            {
                number[] ft = new number[fxy.GetLength(1)];
                for (int y = 0; y < fxy.GetLength(1); y++)
                {
                    ft[y] = fxy[x, y];
                }
                number[] Ff = FFT1(ft);
                int u = x;
                for (int v = 0; v < fxy.GetLength(1); v++)
                {
                    buf[u, v] = Ff[v];
                }
            }

            number[,] Fuv = new number[fxy.GetLength(0), fxy.GetLength(1)];
            for (int y = 0; y < fxy.GetLength(1); y++)
            {
                number[] ft = new number[fxy.GetLength(0)];
                for (int x = 0; x < fxy.GetLength(0); x++)
                {
                    ft[x] = buf[x, y];
                }
                number[] Ff = FFT1(ft);
                int v = y;
                for (int u = 0; u < fxy.GetLength(0); u++)
                {
                    Fuv[u, v] = Ff[u];
                }
            }

            return Fuv;
        }
        public number[,] iFFT2(number[,] Fuv)
        {
            number[,] buf = new number[Fuv.GetLength(0), Fuv.GetLength(1)];
            number MN = new number(Fuv.GetLength(0)*Fuv.GetLength(1),0);
            for (int u = 0; u < Fuv.GetLength(0); u++)
            {
                number[] ft = new number[Fuv.GetLength(1)];
                for (int v = 0; v < Fuv.GetLength(1); v++)
                {
                    Fuv[u, v].image *= -1;
                    Fuv[u, v] = Division(Fuv[u, v], MN);
                    ft[v] = Fuv[u, v];
                }
                number[] Ff = FFT1(ft);
                int x = u;
                for (int y = 0; y < Fuv.GetLength(1); y++)
                {
                    buf[u, y] = Ff[y];
                }
            }

            number[,] fxy = new number[Fuv.GetLength(0), Fuv.GetLength(1)];
            for (int v = 0; v < Fuv.GetLength(1); v++)
            {
                number[] ft = new number[Fuv.GetLength(0)];
                for (int u = 0; u < Fuv.GetLength(0); u++)
                {
                    ft[u] = buf[u, v];
                }
                number[] Ff = FFT1(ft);
                int y = v;
                for (int x = 0; x < Fuv.GetLength(0); x++)
                {
                    fxy[x, y] =Conjucate(Ff[x]);
                }
            }

            return fxy;
        }
        
        public number[] complex(double[] array)
        {
            number[] complex = new number[array.Length];
            for (int x = 0; x < array.Length; x++)
            {
                complex[x] = new number(array[x], 0);
            }
            return complex;
        }
        public number[,] complex(double[,] array)
        {
            number[,] complex = new number[array.GetLength(0), array.GetLength(1)];
            for (int x = 0; x < array.GetLength(0); x++)
            {
                for (int y = 0; y < array.GetLength(1); y++)
                {
                    complex[x, y] = new number(array[x, y], 0);
                }
            }
            return complex;
        }
        public number[,] complex(double[,] abs,double [,]angle)
        {
            number[,] complex = new number[abs.GetLength(0), abs.GetLength(1)];
            number buf = new number(0, 0);
            for (int x = 0; x < abs.GetLength(0); x++)
            {
                for (int y = 0; y < abs.GetLength(1); y++)
                {
                    buf.FromAbsAngle(abs[x,y],angle[x,y]);
                    complex[x, y] = buf;
                }
            }
            return complex;
        }
        public number[,] complex(Bitmap image)
        {
            int w=0,h=0;
            while (image.Width > Math.Pow(2, w)) w++;
            while (image.Height > Math.Pow(2, h)) h++;
            w = (int)Math.Pow(2, w);
            h = (int)Math.Pow(2, h);

            number[,] complex = new number[w, h];
            unsafe
            {
                BitmapData Data = image.LockBits(new Rectangle(0, 0, image.Width, image.Height), ImageLockMode.ReadOnly, PixelFormat.Format24bppRgb);
                int pixelSize = 3;
                for (int y = 0; y < h; y++)
                {
                    for (int x = 0; x < w; x++)
                    {
                        if (x < image.Width && y < image.Height)
                        {
                            double grayScale = ((((byte*)Data.Scan0 + (y * Data.Stride))[x * pixelSize] * 0.11) + (((byte*)Data.Scan0 + (y * Data.Stride))[x * pixelSize + 1] * 0.59) + (((byte*)Data.Scan0 + (y * Data.Stride))[x * pixelSize + 2] * 0.3));
                            complex[x, y] = new number(grayScale, 0);
                        }
                        else complex[x, y] = new number(-1, 0);
                    }
                }
                image.UnlockBits(Data);
            }
            return complex;
        }
        public Bitmap NormalizeImage(double[,] numbers, int width, int height)
        {
            double Min, Max;
            Min = Max = numbers[0, 0];
            for (int x = 0; x < width; x++)
            {
                for (int y = 0; y < height; y++)
                {
                    if (numbers[x, y] > Max) Max = numbers[x, y];
                    else if (numbers[x, y] < Min) Min = numbers[x, y];
                }
            }

            System.Drawing.Bitmap Pic = new System.Drawing.Bitmap(width, height);
            unsafe
            {
                BitmapData newData = Pic.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.WriteOnly, PixelFormat.Format24bppRgb);
                int pixelSize = 3;
                for (int y = 0; y < height; y++)
                {
                    byte* nRow = (byte*)newData.Scan0 + (y * newData.Stride);
                    for (int x = 0; x < width; x++)
                    {
                        double grayScale = (numbers[x, y] - Min) * 255 / (Max - Min);
                        nRow[x * pixelSize] = (byte)grayScale ;
                        nRow[x * pixelSize + 1] = (byte)grayScale;
                        nRow[x * pixelSize + 2] = (byte)grayScale;
                    }
                }
                Pic.UnlockBits(newData);
            }
            return Pic;
        }

        public double[,] abs(number[,] array)
        {
            double[,] abs = new double[array.GetLength(0), array.GetLength(1)];
            for (int x = 0; x < array.GetLength(0); x++)
            {
                for (int y = 0; y < array.GetLength(1); y++)
                {
                    abs[x, y] = array[x, y].abs();
                    if (abs[x, y] < 0)
                    {
                        System.Windows.Forms.MessageBox.Show("g");
                    }
                }
            }
            return abs;
        }
        public double[,] angle(number[,] array)
        {
            double[,] angle = new double[array.GetLength(0), array.GetLength(1)];
            for (int x = 0; x < array.GetLength(0); x++)
            {
                for (int y = 0; y < array.GetLength(1); y++)
                {
                    angle[x, y] = array[x, y].angle();
                }
            }
            return angle;
        }
        public double[,] real(number[,] array)
        {
            double[,] real = new double[array.GetLength(0), array.GetLength(1)];
            for (int x = 0; x < array.GetLength(0); x++)
            {
                for (int y = 0; y < array.GetLength(1); y++)
                {
                    real[x, y] = array[x, y].real;
                }
            }
            return real;
        }
        public double[,] image(number[,] array)
        {
            double[,] image = new double[array.GetLength(0), array.GetLength(1)];
            for (int x = 0; x < array.GetLength(0); x++)
            {
                for (int y = 0; y < array.GetLength(1); y++)
                {
                    image[x, y] = array[x, y].image;
                }
            }
            return image;
        }
        public double[] abs(number[] array)
        {
            double[] abs = new double[array.Count()];
            for (int x = 0; x < array.Count(); x++)
            {
                abs[x] = array[x].abs();
            }
            return abs;
        }
        public double[] angle(number[] array)
        {
            double[] angle = new double[array.GetLength(0)];
            for (int x = 0; x < array.GetLength(0); x++)
            {
                angle[x] = array[x].angle();
            }
            return angle;
        }
        public double[] real(number[] array)
        {
            double[] real = new double[array.Count()];
            for (int x = 0; x < array.Count(); x++)
            {
                real[x] = array[x].real;
            }
            return real;
        }
        public double[] image(number[] array)
        {
            double[] image = new double[array.Count()];
            for (int x = 0; x < array.Count(); x++)
            {
                image[x] = array[x].image;
            }
            return image;
        }
        public double[,] pow(double[,] array)
        {
            double[,] pow = new double[array.GetLength(0), array.GetLength(1)];
            for (int x = 0; x < array.GetLength(0); x++)
            {
                for (int y = 0; y < array.GetLength(1); y++)
                {
                    pow[x, y] = Math.Pow(array[x, y],0.2);
                }
            }
            return pow;
        }

        public double[,] log(double[,] array)
        {
            double[,] log = new double[array.GetLength(0), array.GetLength(1)];
            for (int x = 0; x < array.GetLength(0); x++)
            {
                for (int y = 0; y < array.GetLength(1); y++)
                {
                    log[x, y] = Math.Log(array[x, y]+0.001);
                }
            }
            return log;
        }
        public double[] log(double[] array)
        {
            double[] log = new double[array.Count()];
            for (int x = 0; x < array.Count(); x++)
            {
                log[x] = Math.Log(array[x]+0.001);
            }
            return log;
        }

        public number[,] fftshift(number[,] Fuv)
        {
            number[,] buf = new number[Fuv.GetLength(0), Fuv.GetLength(1)];
            int M = Fuv.GetLength(0) / 2, N = Fuv.GetLength(1) / 2;
            for (int u = 0; u < M; u++)
            {
                for (int v = 0; v < N; v++)
                {
                    buf[u, v] = Fuv[u + M, v + N];
                    buf[u + M, v] = Fuv[u, v + N];
                    buf[u, v + N] = Fuv[u + M, v];
                    buf[u + M, v + N] = Fuv[u, v];
                }
            }
            return buf;
        }
        public number[,] ifftshift(number[,] Fuv)
        {
            return fftshift(Fuv);
        }




    }
}
