using UnityEngine;

namespace FluidFlow
{
    #region Simulator
    public class FluidFlow2DSimulator
    {
        private int Size;
        private float Dt; // time step
        private float Diff; // Diffuse
        private float Visc; // Viscosity

        float[] S; // Previous density
        float[] Density;

        float[] Vx;
        float[] Vy;

        float[] Vx0;
        float[] Vy0;

        public FluidFlow2DSimulator(int size, float diffusion, float viscosity, float deltaTime)
        {
            Size = size;
            Diff = diffusion;
            Visc = viscosity;
            Dt = deltaTime;

            int totalSize = size * size;

            S = new float[totalSize];
            Density = new float[totalSize];

            Vx = new float[totalSize];
            Vy = new float[totalSize];

            Vx0 = new float[totalSize];
            Vy0 = new float[totalSize];
        }
        
        private int IndexOf(int x, int y)
        {
            return x + y * Size;
        }

        /// <summary>
        /// Update Fluid flow system per step
        /// </summary>
        public void Step()
        {
            int size          = this.Size;
            float visc     = this.Visc;
            float diff     = this.Diff;
            float dt       = this.Dt;
            float[] vx      = this.Vx;
            float[] vy      = this.Vy;
            float[] vx0     = this.Vx0;
            float[] vy0     = this.Vy0;
            float[] s       = this.S;
            float[] density = this.Density;
    
            Diffuse(1, vx0, vx, visc, dt, 4, size);
            Diffuse(2, vy0, vy, visc, dt, 4, size);

            Project(vx0, vy0, vx, vy, 4, size);
    
            Advect(1, vx, vx0, vx0, vy0, dt, size);
            Advect(2, vy, vy0, vx0, vy0, dt, size);

            Project(vx, vy, Vx0, Vy0, 4, size);
    
            Diffuse(0, s, density, diff, dt, 4, size);
            Advect(0, density, s, vx, vy, dt, size);
        }
        
        private void AddDensity(int x, int y, float amount)
        {
            int index = IndexOf(x, y);
            Density[index] += amount;
        }

        private void AddVelocity(int x, int y, int amountX, int amountY)
        {
            int index = IndexOf(x, y);

            Vx[index] += amountX;
            Vy[index] += amountY;
        }

        /// <summary>
        /// Diffuse
        /// </summary>
        /// <param name="b">for calc boundary</param>
        /// <param name="x">Density or velocity array</param>
        /// <param name="x0">Original density or original velocity array</param>
        /// <param name="diff">Diffusion</param>
        /// <param name="dt">Delta time</param>
        /// <param name="iter">Iterations</param>
        /// <param name="size">Size</param>
        public void Diffuse(int b, float[] x, float[] x0, float diff, float dt, int iter, int size)
        {
            float a = dt * diff * (size - 2) * (size - 2);
            Lin_Solve(b, x, x0, a, 1 + 6 * a, iter, size);
        }

        public void Project(float[] velocX, float[] velocY, float[] p, float[] div, int iter, int size)
        {
            for (int j = 1; j < size - 1; j++)
            {
                for (int i = 1; i < size - 1; i++)
                {
                    div[IndexOf(i, j)] = -0.5f * (
                        velocX[IndexOf(i + 1, j)]
                        - velocX[IndexOf(i - 1, j)]
                        + velocY[IndexOf(i, j + 1)]
                        - velocY[IndexOf(i, j - 1)]
                    ) / size;
                    p[IndexOf(i, j)] = 0;
                }
            }

            Set_Bnd(0, div, size);
            Set_Bnd(0, p, size);
            Lin_Solve(0, p, div, 1, 6, iter, size);

            for (int j = 1; j < size - 1; j++)
            {
                for (int i = 1; i < size - 1; i++)
                {
                    velocX[IndexOf(i, j)] -= 0.5f * (p[IndexOf(i + 1, j)]
                                                     - p[IndexOf(i - 1, j)]) * size;
                    velocY[IndexOf(i, j)] -= 0.5f * (p[IndexOf(i, j + 1)]
                                                     - p[IndexOf(i, j - 1)]) * size;
                }
            }

            Set_Bnd(1, velocX, size);
            Set_Bnd(2, velocY, size);
        }

        public void Advect(int b, float[] d, float[] d0, float[] velocX, float[] velocY, float dt, int size)
        {
            float i0, i1, j0, j1;

            float dtx = dt * (size - 2);
            float dty = dt * (size - 2);

            float s0, s1, t0, t1;
            float tmp1, tmp2, x, y;

            float Nfloat = size;
            float ifloat, jfloat;
            int i, j;

            for (j = 1, jfloat = 1; j < size - 1; j++, jfloat++)
            {
                for (i = 1, ifloat = 1; i < size - 1; i++, ifloat++)
                {
                    tmp1 = dtx * velocX[IndexOf(i, j)];
                    tmp2 = dty * velocY[IndexOf(i, j)];

                    x = ifloat - tmp1;
                    y = jfloat - tmp2;

                    if (x < 0.5f)
                    {
                        x = 0.5f;
                    }

                    if (x > Nfloat + 0.5f)
                    {
                        x = Nfloat + 0.5f;
                    }

                    i0 = Mathf.Floor(x);
                    i1 = i0 + 1.0f;

                    if (y < 0.5f)
                    {
                        y = 0.5f;
                    }

                    if (y > Nfloat + 0.5f)
                    {
                        y = Nfloat + 0.5f;
                    }

                    j0 = Mathf.Floor(y);
                    j1 = j0 + 1.0f;

                    s1 = x - i0;
                    s0 = 1.0f - s1;
                    t1 = y - j0;
                    t0 = 1.0f - t1;

                    int i0i = (int) i0;
                    int i1i = (int) i1;
                    int j0i = (int) j0;
                    int j1i = (int) j1;

                    d[IndexOf(i, j)] =
                        s0 * (t0 * d0[IndexOf(i0i, j0i)] + t1 * d0[IndexOf(i0i, j1i)])
                        + s1 * (t0 * d0[IndexOf(i1i, j0i)] + t1 * d0[IndexOf(i1i, j1i)]);
                }
            }

            Set_Bnd(b, d, size);
        }

        public void Lin_Solve(int b, float[] x, float[] x0, float a, float c, int iter, int size)
        {
            float cRecip = 1.0f / c;
            for (int k = 0; k < iter; k++) // This is Solve. The more iterations the more accurate the results 
            {
                for (int j = 1; j < size - 1; j++) // Begin 1 and end at (size - 2) is this for boundary? 
                {
                    for (int i = 1; i < size - 1; i++) // Begin 1 and end at (size - 2) is this for boundary?
                    {
                        x[IndexOf(i, j)] =
                            (x0[IndexOf(i, j)]
                             + a * (x[IndexOf(i - 1, j)]
                                    + x[IndexOf(i + 1, j)]
                                    + x[IndexOf(i, j - 1)]
                                    + x[IndexOf(i, j + 1)])
                            ) * cRecip;
                    }
                }
            }

            Set_Bnd(b, x, size);
        }

        /// <summary>
        /// Set boundary
        /// </summary>
        /// <param name="b">Left-right-top-down</param>
        /// <param name="x"></param>
        /// <param name="size"></param>
        private void Set_Bnd(int b, float[] x, int size)
        {
            for (int k = 1; k < size - 1; k++)
            {
                for (int i = 1; i < size - 1; i++)
                {
                    x[IndexOf(i, 0)] = b == 2 ? -x[IndexOf(i, 1)] : x[IndexOf(i, 1)];
                    x[IndexOf(i, size - 1)] = b == 2 ? -x[IndexOf(i, size - 2)] : x[IndexOf(i, size - 2)];
                }
            }

            for (int k = 1; k < size - 1; k++)
            {
                for (int j = 1; j < size - 1; j++)
                {
                    x[IndexOf(0, j)] = b == 1 ? -x[IndexOf(1, j)] : x[IndexOf(1, j)];
                    x[IndexOf(size - 1, j)] = b == 1 ? -x[IndexOf(size - 2, j)] : x[IndexOf(size - 2, j)];
                }
            }

            x[IndexOf(0, 0)] = 0.5f * (x[IndexOf(1, 0)] + x[IndexOf(0, 1)]);
            
            x[IndexOf(0, size- 1)] = 0.5f * (x[IndexOf(1, size - 1)] + x[IndexOf(0, size - 2)]);

            x[IndexOf(size - 1, 0)] = 0.5f * (x[IndexOf(size - 2, 0)] + x[IndexOf(size - 1, 1)]);
            
            x[IndexOf(size - 1, size - 1)] = 0.5f * (x[IndexOf(size - 2, size - 1)] + x[IndexOf(size - 1, size - 2)]);
        }
    }
    #endregion

    #region Class
    public class Fluid2DPoint
    {
        
    }
    #endregion
}
