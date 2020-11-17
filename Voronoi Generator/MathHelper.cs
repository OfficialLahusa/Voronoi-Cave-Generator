using SFML.System;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Voronoi_Generator
{
    public static class MathHelper
    {
        // Get 2D position of sample point on subdivided circle with optional scaling 
        public static Vector2f getSamplePoint(uint pos, uint max, float xscale = 1.0f, float yscale = 1.0f)
        {
            float radiansPerStep = 2.0f * (float)Math.PI / max;
            return new Vector2f(xscale * (0.5f + (float)Math.Cos(pos * radiansPerStep) / 2.0f), yscale * (0.5f + (float)Math.Sin(pos * radiansPerStep) / 2.0f));
        }
    }

    public static class DecimalExtensions
    {
        public static decimal Map(this decimal value, decimal fromSource, decimal toSource, decimal fromTarget, decimal toTarget)
        {
            return (value - fromSource) / (toSource - fromSource) * (toTarget - fromTarget) + fromTarget;
        }
    }

    public static class IntExtensions
    {
        public static int Map(this int value, int fromSource, int toSource, int fromTarget, int toTarget)
        {
            return (value - fromSource) / (toSource - fromSource) * (toTarget - fromTarget) + fromTarget;
        }
    }

    public static class FloatExtensions
    {
        public static float Map(this float value, float fromSource, float toSource, float fromTarget, float toTarget)
        {
            return (value - fromSource) / (toSource - fromSource) * (toTarget - fromTarget) + fromTarget;
        }
    }

    public static class DoubleExtensions
    {
        public static double Map(this double value, double fromSource, double toSource, double fromTarget, double toTarget)
        {
            return (value - fromSource) / (toSource - fromSource) * (toTarget - fromTarget) + fromTarget;
        }
    }
}
