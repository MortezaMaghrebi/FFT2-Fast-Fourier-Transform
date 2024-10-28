using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Imaging;
using System.Linq;
using System.Text;

namespace Fourier
{
    public class ImageProcessingFourier
    {

        public class FFT2
        {
            public Fourier fourier = new Fourier();
            Fourier.number[,] Fuv,fxy;
            int width, height;

            public void fft2(Bitmap picture)
            {
                width = picture.Width;
                height = picture.Height;
                this.fxy =  fourier.complex(picture);
                Fuv = fourier.FFT2(fxy);
            }
            public Bitmap ifft2()
            {
                fxy = fourier.iFFT2(Fuv);
                return fourier.NormalizeImage(fourier.abs(fxy), width, height);
            }

            public Bitmap abs()
            {
                return fourier.NormalizeImage(fourier.pow(fourier.abs(fftshift(Fuv))),Fuv.GetLength(0),Fuv.GetLength(1));
            }
            public Bitmap angle()
            {
                return fourier.NormalizeImage(fourier.angle(fftshift(Fuv)), Fuv.GetLength(0), Fuv.GetLength(1));
            }
            
            public void Filter(double[,] Hf)
            {
                Fourier.number[,] shifted = fftshift(Fuv);
                double abs, angle;
                int Dx;
                int Dy;

                for (int x = 0; x < Fuv.GetLength(0); x++)
                {
                    for (int y = 0; y < Fuv.GetLength(1); y++)
                    {
                        Dx = x * Hf.GetLength(0) / Fuv.GetLength(0);
                        Dy = y * Hf.GetLength(1) / Fuv.GetLength(1);
                        abs = shifted[x, y].abs();
                        angle = shifted[x, y].angle();
                        if (Hf[Dx, Dy] < 0) Hf[Dx, Dy] = 0;
                        else if (Hf[Dx, Dy] > 1) Hf[Dx, Dy] = 1;
                        abs *= Hf[Dx, Dy];
                        shifted[x, y].FromAbsAngle(abs, angle);
                    }
                }
                Fuv = ifftshift(shifted);
            }
            
            
            Fourier.number[,] fftshift(Fourier.number[,] Fuv)
            {
                Fourier.number[,] buf = new Fourier.number[Fuv.GetLength(0), Fuv.GetLength(1)];
                int M = Fuv.GetLength(0)/2,N = Fuv.GetLength(1)/2;
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
            Fourier.number[,] ifftshift(Fourier.number[,] Fuv)
            {
                return fftshift(Fuv);
            }

        }

 

    }
}
