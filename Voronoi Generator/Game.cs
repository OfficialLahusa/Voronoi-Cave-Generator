using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;
using SFML.System;
using SFML.Window;
using SFML.Graphics;

namespace Voronoi_Generator
{
    public class Game
    {
        private const uint POINT_COUNT = 1200; // 400
        private const uint WIDTH = 1200;
        private const uint HEIGHT = 1200;
        private const double SUBDIV_MIN_LEN = 1.0;
        private const double SUBDIVISION_DISTANCE = 1.0f; // 6.0f
        private const float EDGE_DISTORTION_SCALE = 0.5f; // 0.5f
        private const uint SMOOTHING_ITERATIONS = 10;

        private RenderWindow window = null;
        private Clock clock;
        private float deltaTime;
        private Random random;

        private List<VoronoiLib.Structures.FortuneSite> points;
        private LinkedList<VoronoiLib.Structures.VEdge> voronoi;
        private Dictionary<VoronoiLib.Structures.FortuneSite, VoronoiCell> cellDictionary;
        private FastNoise cellNoise;
        private FastNoise lineDistortionNoise;

        /*private Image image;
        private Texture imageTexture;
        private RectangleShape imageRect;
        private const uint DEBUG_POLYCOUNT = 200;
        private VertexArray smoothingCircle;
        private VertexArray distortedLine;*/

        private VertexArray walls;
        private VertexArray centers;
        private VertexArray cells;
        private VertexArray subdividedCells;


        public Game()
        {

        }

        // http://www-cs-students.stanford.edu/~amitp/game-programming/polygon-map-generation/
        // https://leatherbee.org/index.php/2018/10/28/terrain-generation-4-plates-continents-coasts/
        public void Run()
        {
            window = new RenderWindow(new VideoMode(WIDTH, HEIGHT), "Voronoi Cave Generation", Styles.Default, new ContextSettings(0, 0, 4));
            window.Closed += OnClose;

            /*image = new Image(WIDTH, HEIGHT);
            imageRect = new RectangleShape(new Vector2f(WIDTH, HEIGHT));
            imageRect.Position = new Vector2f(0, 0);

            smoothingCircle = new VertexArray(PrimitiveType.LineStrip);
            for(uint i = 0; i <= DEBUG_POLYCOUNT; i++)
            {
                smoothingCircle.Append(new Vertex(MathHelper.getSamplePoint(i, DEBUG_POLYCOUNT, WIDTH, HEIGHT), Color.Blue));
            }*/

            random = new Random();

            cellNoise = new FastNoise(random.Next());
            cellNoise.SetNoiseType(FastNoise.NoiseType.PerlinFractal);
            cellNoise.SetFrequency(0.005f);

            lineDistortionNoise = new FastNoise(random.Next());
            lineDistortionNoise.SetNoiseType(FastNoise.NoiseType.PerlinFractal);
            lineDistortionNoise.SetFrequency(1.0f);

            RegenerateVoronoi();

            //Game Loop
            clock = new Clock();
            while (window.IsOpen)
            {
                deltaTime = clock.Restart().AsSeconds();
                Console.Write("\r" + Math.Round(1.0f / deltaTime) + " FPS - Frametime: " + deltaTime + "s                           ");

                //Event handling
                window.DispatchEvents();

                //Input handling
                #region InputHandling
                if (Keyboard.IsKeyPressed(Keyboard.Key.Escape) && window.HasFocus())
                {
                    window.Close();
                }
                if(Keyboard.IsKeyPressed(Keyboard.Key.Space) && window.HasFocus())
                {
                    RegenerateVoronoi();
                }
                #endregion

                // Rendering
                #region Rendering
                window.Clear(new Color(0, 0, 0));

                if(Keyboard.IsKeyPressed(Keyboard.Key.Tab) && window.HasFocus())
                {
                    window.Draw(cells);
                } else
                {
                    window.Draw(subdividedCells);
                }

                //window.Draw(walls);

                //window.Draw(centers);

                /*window.Draw(imageRect);

                window.Draw(smoothingCircle);
                window.Draw(distortedLine);*/

                window.Display();
                #endregion
            }
        }

        private void RegenerateVoronoi()
        {
            cellNoise.SetSeed(random.Next());
            lineDistortionNoise.SetSeed(random.Next());

            #region DistortionDebugCircle
            /*for (uint x = 0; x < WIDTH; x++)
            {
                for(uint y = 0; y < HEIGHT; y++)
                {
                    float noiseVal = (1 + lineDistortionNoise.GetPerlinFractal(x / (float)WIDTH, y / (float)HEIGHT)) / 2.0f;
                    byte grayscaleVal = (byte)(255.0f * noiseVal);
                    image.SetPixel(x, y, new Color(grayscaleVal, grayscaleVal, grayscaleVal));
                }
            }

            imageTexture = new Texture(image);
            imageRect.Texture = imageTexture;

            // Precalc highest noise value for distorted line as described in https://leatherbee.org/index.php/2018/10/28/terrain-generation-4-plates-continents-coasts/
            float cutoffValue = lineDistortionNoise.GetPerlinFractal(MathHelper.getSamplePoint(0, DEBUG_POLYCOUNT).X, MathHelper.getSamplePoint(0, DEBUG_POLYCOUNT).Y);
            float highestAbs = 1.0f;
            for (uint i = 0; i <= DEBUG_POLYCOUNT; i++)
            {
                float noiseVal = Math.Abs(lineDistortionNoise.GetPerlinFractal(MathHelper.getSamplePoint(i, DEBUG_POLYCOUNT).X, MathHelper.getSamplePoint(i, DEBUG_POLYCOUNT).Y) - cutoffValue);
                if (noiseVal > highestAbs) {
                    highestAbs = noiseVal;
                }
            }

            distortedLine = new VertexArray(PrimitiveType.LineStrip);
            for (uint i = 0; i <= DEBUG_POLYCOUNT; i++)
            {
                distortedLine.Append(new Vertex(new Vector2f(WIDTH * i * (1.0f / DEBUG_POLYCOUNT), HEIGHT * (0.5f + 0.5f * (lineDistortionNoise.GetPerlinFractal(MathHelper.getSamplePoint(i, DEBUG_POLYCOUNT).X, MathHelper.getSamplePoint(i, DEBUG_POLYCOUNT).Y) - cutoffValue) / highestAbs)), Color.Red));
            }

            Console.WriteLine("\r" + lineDistortionNoise.GetPerlinFractal(MathHelper.getSamplePoint(0, DEBUG_POLYCOUNT).X, MathHelper.getSamplePoint(0, DEBUG_POLYCOUNT).Y) + "                                                 ");*/
            #endregion

            // Original randomized point cloud
            #region PointCloud
            points = new List<VoronoiLib.Structures.FortuneSite>();

            for (uint i = 0; i < POINT_COUNT; i++)
            {
                points.Add(new VoronoiLib.Structures.FortuneSite(random.NextDouble() * WIDTH, random.NextDouble() * HEIGHT));
            }
            #endregion

            // Multiple iterations of Lloyd's algorithm to relax the point cloud
            #region Relaxations
            for (uint j = 0; j < SMOOTHING_ITERATIONS; j++)
            {
                // Pre-Gen to identify points that are too close to eachother
                voronoi = VoronoiLib.FortunesAlgorithm.Run(points, 0, 0, WIDTH, HEIGHT);

                cellDictionary = new Dictionary<VoronoiLib.Structures.FortuneSite, VoronoiCell>();

                foreach (var line in voronoi)
                {
                    if (cellDictionary.ContainsKey(line.Left))
                    {
                        if (cellDictionary[line.Left] == null)
                        {
                            cellDictionary[line.Left] = new VoronoiCell();
                        }
                        cellDictionary[line.Left].cellWalls.Add(line);
                    }
                    else
                    {
                        cellDictionary.Add(line.Left, new VoronoiCell());
                        cellDictionary[line.Left].cellWalls.Add(line);
                    }

                    if (cellDictionary.ContainsKey(line.Right))
                    {
                        if (cellDictionary[line.Right] == null)
                        {
                            cellDictionary[line.Right] = new VoronoiCell();
                        }
                        cellDictionary[line.Right].cellWalls.Add(line);
                    }
                    else
                    {
                        cellDictionary.Add(line.Right, new VoronoiCell());
                        cellDictionary[line.Right].cellWalls.Add(line);
                    }
                }

                points = new List<VoronoiLib.Structures.FortuneSite>();

                foreach (var cell in cellDictionary.Values)
                {
                    Vector2f addedPositions = new Vector2f(0, 0);

                    for (int i = 0; i < cell.cellWalls.Count; i++)
                    {
                        addedPositions.X += (float)cell.cellWalls[i].Start.X + (float)cell.cellWalls[i].End.X;
                        addedPositions.Y += (float)cell.cellWalls[i].Start.Y + (float)cell.cellWalls[i].End.Y;
                    }

                    addedPositions /= cell.cellWalls.Count * 2;
                    points.Add(new VoronoiLib.Structures.FortuneSite(addedPositions.X, addedPositions.Y));
                }
            }
            #endregion

            // Generate final voronoi graph
            voronoi = VoronoiLib.FortunesAlgorithm.Run(points, 0, 0, WIDTH, HEIGHT);

            // Generate outline VertexArray for the voronoi graph
            #region Outlines
            walls = new VertexArray(PrimitiveType.Lines);
            foreach (var line in voronoi)
            {
                walls.Append(new Vertex(new Vector2f((float)line.Start.X, (float)line.Start.Y), Color.Black));
                walls.Append(new Vertex(new Vector2f((float)line.End.X, (float)line.End.Y), Color.Black));
            }
            #endregion

            // Generate VertexArray for center points
            #region CenterPoints
            centers = new VertexArray(PrimitiveType.Points);
            foreach (var center in points)
            {
                centers.Append(new Vertex(new Vector2f((float)center.X, (float)center.Y), new Color((byte)(Math.Min(center.Neighbors.Count * (255.0 / 10.0), 255)), 20, 20)));
            }
            #endregion

            // Organize cells into dictionary structure
            #region OrganizeCells
            cellDictionary = new Dictionary<VoronoiLib.Structures.FortuneSite, VoronoiCell>();

            foreach(var line in voronoi)
            {
                if (cellDictionary.ContainsKey(line.Left))
                {
                    if (cellDictionary[line.Left] == null) {
                        cellDictionary[line.Left] = new VoronoiCell();
                    }
                    cellDictionary[line.Left].cellWalls.Add(line);
                } else
                {
                    cellDictionary.Add(line.Left, new VoronoiCell());
                    cellDictionary[line.Left].cellWalls.Add(line);
                }

                if (cellDictionary.ContainsKey(line.Right))
                {
                    if (cellDictionary[line.Right] == null)
                    {
                        cellDictionary[line.Right] = new VoronoiCell();
                    }
                    cellDictionary[line.Right].cellWalls.Add(line);
                }
                else
                {
                    cellDictionary.Add(line.Right, new VoronoiCell());
                    cellDictionary[line.Right].cellWalls.Add(line);
                }
            }
            #endregion

            // Register missing margin triangles (also accounting for corners)
            #region MarginStitching
            foreach (var cell in cellDictionary.Keys)
            {
                // Counters for how many points intersect with which boundary
                uint top = 0, left = 0, bottom = 0, right = 0;

                foreach (var line in cellDictionary[cell].cellWalls)
                {
                    if(line.Start.X == 0 || line.End.X == 0)
                    {
                        left++;
                    }
                    if (line.Start.Y == 0 || line.End.Y == 0)
                    {
                        top++;
                    }
                    if (line.Start.X == WIDTH || line.End.X == WIDTH)
                    {
                        right++;
                    }
                    if (line.Start.Y == HEIGHT || line.End.Y == HEIGHT)
                    {
                        bottom++;
                    }
                }

                if(top != 0 || left != 0 || bottom != 0 || right != 0)
                {
                    cellDictionary[cell].isBorder = true;
                }

                if(top == 2)
                {
                    VoronoiLib.Structures.VEdge firstEdge = cellDictionary[cell].cellWalls.Find(x => x.Start.Y == 0 || x.End.Y == 0);
                    VoronoiLib.Structures.VPoint firstPoint;
                    if (firstEdge.Start.Y == 0) firstPoint = firstEdge.Start;
                    else firstPoint = firstEdge.End;

                    VoronoiLib.Structures.VEdge secondEdge = cellDictionary[cell].cellWalls.Find(x => x != firstEdge && (x.Start.Y == 0 || x.End.Y == 0));
                    VoronoiLib.Structures.VPoint secondPoint;
                    if (secondEdge.Start.Y == 0) secondPoint = secondEdge.Start;
                    else secondPoint = secondEdge.End;

                    cellDictionary[cell].additionalEdges.Add(new KeyValuePair<Vector2f, Vector2f>(new Vector2f((float)firstPoint.X, (float)firstPoint.Y), new Vector2f((float)secondPoint.X, (float)secondPoint.Y)));
                }
                if (left == 2)
                {
                    VoronoiLib.Structures.VEdge firstEdge = cellDictionary[cell].cellWalls.Find(x => x.Start.X == 0 || x.End.X == 0);
                    VoronoiLib.Structures.VPoint firstPoint;
                    if (firstEdge.Start.X == 0) firstPoint = firstEdge.Start;
                    else firstPoint = firstEdge.End;

                    VoronoiLib.Structures.VEdge secondEdge = cellDictionary[cell].cellWalls.Find(x => x != firstEdge && (x.Start.X == 0 || x.End.X == 0));
                    VoronoiLib.Structures.VPoint secondPoint;
                    if (secondEdge.Start.X == 0) secondPoint = secondEdge.Start;
                    else secondPoint = secondEdge.End;

                    cellDictionary[cell].additionalEdges.Add(new KeyValuePair<Vector2f, Vector2f>(new Vector2f((float)firstPoint.X, (float)firstPoint.Y), new Vector2f((float)secondPoint.X, (float)secondPoint.Y)));
                }
                if (bottom == 2)
                {
                    VoronoiLib.Structures.VEdge firstEdge = cellDictionary[cell].cellWalls.Find(x => x.Start.Y == HEIGHT || x.End.Y == HEIGHT);
                    VoronoiLib.Structures.VPoint firstPoint;
                    if (firstEdge.Start.Y == HEIGHT) firstPoint = firstEdge.Start;
                    else firstPoint = firstEdge.End;

                    VoronoiLib.Structures.VEdge secondEdge = cellDictionary[cell].cellWalls.Find(x => x != firstEdge && (x.Start.Y == HEIGHT || x.End.Y == HEIGHT));
                    VoronoiLib.Structures.VPoint secondPoint;
                    if (secondEdge.Start.Y == HEIGHT) secondPoint = secondEdge.Start;
                    else secondPoint = secondEdge.End;

                    cellDictionary[cell].additionalEdges.Add(new KeyValuePair<Vector2f, Vector2f>(new Vector2f((float)firstPoint.X, (float)firstPoint.Y), new Vector2f((float)secondPoint.X, (float)secondPoint.Y)));
                }
                if (right == 2)
                {
                    VoronoiLib.Structures.VEdge firstEdge = cellDictionary[cell].cellWalls.Find(x => x.Start.X == WIDTH || x.End.X == WIDTH);
                    VoronoiLib.Structures.VPoint firstPoint;
                    if (firstEdge.Start.X == WIDTH) firstPoint = firstEdge.Start;
                    else firstPoint = firstEdge.End;

                    VoronoiLib.Structures.VEdge secondEdge = cellDictionary[cell].cellWalls.Find(x => x != firstEdge && (x.Start.X == WIDTH || x.End.X == WIDTH));
                    VoronoiLib.Structures.VPoint secondPoint;
                    if (secondEdge.Start.X == WIDTH) secondPoint = secondEdge.Start;
                    else secondPoint = secondEdge.End;

                    cellDictionary[cell].additionalEdges.Add(new KeyValuePair<Vector2f, Vector2f>(new Vector2f((float)firstPoint.X, (float)firstPoint.Y), new Vector2f((float)secondPoint.X, (float)secondPoint.Y)));
                }
                if (top == 1 && left == 1)
                {
                    VoronoiLib.Structures.VEdge firstEdge = cellDictionary[cell].cellWalls.Find(x => x.Start.Y == 0 || x.End.Y == 0);
                    VoronoiLib.Structures.VPoint firstPoint;
                    if (firstEdge.Start.Y == 0) firstPoint = firstEdge.Start;
                    else firstPoint = firstEdge.End;

                    VoronoiLib.Structures.VEdge secondEdge = cellDictionary[cell].cellWalls.Find(x =>x.Start.X == 0 || x.End.X == 0);
                    VoronoiLib.Structures.VPoint secondPoint;
                    if (secondEdge.Start.X == 0) secondPoint = secondEdge.Start;
                    else secondPoint = secondEdge.End;

                    cellDictionary[cell].additionalEdges.Add(new KeyValuePair<Vector2f, Vector2f>(new Vector2f((float)firstPoint.X, (float)firstPoint.Y), new Vector2f(0, 0)));
                    cellDictionary[cell].additionalEdges.Add(new KeyValuePair<Vector2f, Vector2f>(new Vector2f(0, 0), new Vector2f((float)secondPoint.X, (float)secondPoint.Y)));
                }
                if (top == 1 && right == 1)
                {
                    VoronoiLib.Structures.VEdge firstEdge = cellDictionary[cell].cellWalls.Find(x => x.Start.Y == 0 || x.End.Y == 0);
                    VoronoiLib.Structures.VPoint firstPoint;
                    if (firstEdge.Start.Y == 0) firstPoint = firstEdge.Start;
                    else firstPoint = firstEdge.End;

                    VoronoiLib.Structures.VEdge secondEdge = cellDictionary[cell].cellWalls.Find(x => x.Start.X == WIDTH || x.End.X == WIDTH);
                    VoronoiLib.Structures.VPoint secondPoint;
                    if (secondEdge.Start.X == WIDTH) secondPoint = secondEdge.Start;
                    else secondPoint = secondEdge.End;

                    cellDictionary[cell].additionalEdges.Add(new KeyValuePair<Vector2f, Vector2f>(new Vector2f((float)firstPoint.X, (float)firstPoint.Y), new Vector2f(WIDTH, 0)));
                    cellDictionary[cell].additionalEdges.Add(new KeyValuePair<Vector2f, Vector2f>(new Vector2f(WIDTH, 0), new Vector2f((float)secondPoint.X, (float)secondPoint.Y)));
                }
                if (bottom == 1 && left == 1)
                {
                    VoronoiLib.Structures.VEdge firstEdge = cellDictionary[cell].cellWalls.Find(x => x.Start.Y == HEIGHT || x.End.Y == HEIGHT);
                    VoronoiLib.Structures.VPoint firstPoint;
                    if (firstEdge.Start.Y == HEIGHT) firstPoint = firstEdge.Start;
                    else firstPoint = firstEdge.End;

                    VoronoiLib.Structures.VEdge secondEdge = cellDictionary[cell].cellWalls.Find(x => x.Start.X == 0 || x.End.X == 0);
                    VoronoiLib.Structures.VPoint secondPoint;
                    if (secondEdge.Start.X == 0) secondPoint = secondEdge.Start;
                    else secondPoint = secondEdge.End;

                    cellDictionary[cell].additionalEdges.Add(new KeyValuePair<Vector2f, Vector2f>(new Vector2f((float)firstPoint.X, (float)firstPoint.Y), new Vector2f(0, HEIGHT)));
                    cellDictionary[cell].additionalEdges.Add(new KeyValuePair<Vector2f, Vector2f>(new Vector2f(0, HEIGHT), new Vector2f((float)secondPoint.X, (float)secondPoint.Y)));
                }
                if (bottom == 1 && right == 1)
                {
                    VoronoiLib.Structures.VEdge firstEdge = cellDictionary[cell].cellWalls.Find(x => x.Start.Y == HEIGHT || x.End.Y == HEIGHT);
                    VoronoiLib.Structures.VPoint firstPoint;
                    if (firstEdge.Start.Y == HEIGHT) firstPoint = firstEdge.Start;
                    else firstPoint = firstEdge.End;

                    VoronoiLib.Structures.VEdge secondEdge = cellDictionary[cell].cellWalls.Find(x => x.Start.X == WIDTH || x.End.X == WIDTH);
                    VoronoiLib.Structures.VPoint secondPoint;
                    if (secondEdge.Start.X == WIDTH) secondPoint = secondEdge.Start;
                    else secondPoint = secondEdge.End;

                    cellDictionary[cell].additionalEdges.Add(new KeyValuePair<Vector2f, Vector2f>(new Vector2f((float)firstPoint.X, (float)firstPoint.Y), new Vector2f(WIDTH, HEIGHT)));
                    cellDictionary[cell].additionalEdges.Add(new KeyValuePair<Vector2f, Vector2f>(new Vector2f(WIDTH, HEIGHT), new Vector2f((float)secondPoint.X, (float)secondPoint.Y)));
                }
            }
            #endregion

            // Precalculate Noise
            #region PrecalcNoise

            foreach(var cell in cellDictionary.Keys)
            {
                byte noiseVal = (byte)(255.0 * cellNoise.GetPerlinFractal((float)cell.X, (float)cell.Y));
                cellDictionary[cell].isSolid = cellDictionary[cell].isBorder || noiseVal > 127;
            }
            #endregion

            // Generate VertexArray for Cell background
            #region Cells
            cells = new VertexArray(PrimitiveType.Triangles);
            subdividedCells = new VertexArray(PrimitiveType.Triangles);
            foreach (var cell in cellDictionary.Keys)
            {
                if(!cellDictionary[cell].isSolid)
                {
                    continue;
                }

                byte colorVal = (cellDictionary[cell].isSolid) ? (byte)255 : (byte)0;
                Color cellColor = new Color(colorVal, colorVal, colorVal);

                foreach (var line in cellDictionary[cell].cellWalls)
                {
                    cells.Append(new Vertex(new Vector2f((float)cell.X, (float)cell.Y), cellColor));
                    cells.Append(new Vertex(new Vector2f((float)line.Start.X, (float)line.Start.Y), cellColor));
                    cells.Append(new Vertex(new Vector2f((float)line.End.X, (float)line.End.Y), cellColor));

                    float length = (float)Math.Sqrt(Math.Pow(line.End.X - line.Start.X, 2) + Math.Pow(line.End.Y - line.Start.Y, 2));

                    // Subdivide the lines ONLY if it's between a solid and a passable cell AND if it's longer than the minimum required length for subdivision
                    if (cellDictionary[line.Left].isSolid != cellDictionary[line.Right].isSolid && length > SUBDIV_MIN_LEN)
                    {
                        uint subdivisionCount = (uint)Math.Abs(Math.Floor(Math.Sqrt(Math.Pow(line.Start.X - line.End.X, 2) + Math.Pow(line.Start.Y - line.End.Y, 2)) / SUBDIVISION_DISTANCE));

                        Vector2f v = new Vector2f((float)line.End.X - (float)line.Start.X, (float)line.End.Y - (float)line.Start.Y);
                        Vector2f inc = new Vector2f(v.X / (subdivisionCount + 1), v.Y / (subdivisionCount + 1));
                        Vector2f normal = new Vector2f(-v.Y, v.X);
                        Vector2f mid = new Vector2f((float)line.Start.X + v.X / 2, (float)line.Start.Y + v.Y / 2);
                        float len = (float)Math.Sqrt(Math.Pow(normal.X, 2) + Math.Pow(normal.Y, 2));
                        normal = new Vector2f(normal.X / len, normal.Y / len);

                        if (Math.Sqrt(Math.Pow(cell.X - (line.Start.X + v.X / 2 + normal.X), 2) + Math.Pow(cell.Y - (line.Start.Y + v.Y / 2 + normal.Y), 2)) >= Math.Sqrt(Math.Pow(cell.X - (line.Start.X + v.X / 2), 2) + Math.Pow(cell.Y - (line.Start.Y + v.Y / 2), 2)))
                        {
                            normal *= -1;
                        }

                        float cutoffValue = lineDistortionNoise.GetPerlinFractal(MathHelper.getSamplePoint(0, subdivisionCount).X + mid.X, MathHelper.getSamplePoint(0, subdivisionCount).Y + mid.Y);
                        float highestAbs = 1.0f; // Default to 1.0f to avoid overstretching of smaller noise peaks
                        for (uint i = 0; i <= subdivisionCount; i++)
                        {
                            float noiseVal = Math.Abs(lineDistortionNoise.GetPerlinFractal(MathHelper.getSamplePoint(i, subdivisionCount + 1).X, MathHelper.getSamplePoint(i, subdivisionCount + 1).Y) - cutoffValue);
                            if (noiseVal > highestAbs)
                            {
                                highestAbs = noiseVal;
                            }
                        }

                        float offx1 = 0, offy1 = 0, offx2 = 0, offy2 = 0;
                        float distNeeded = (subdivisionCount + 1.0f) / 2.0f;

                        for(uint i = 0; i <= subdivisionCount; i++)
                        {
                            // Compute offsets of current and next point
                            offx1 = normal.X * EDGE_DISTORTION_SCALE * len / 2.0f * (lineDistortionNoise.GetPerlinFractal(MathHelper.getSamplePoint(i, subdivisionCount + 1).X + mid.X, MathHelper.getSamplePoint(i, subdivisionCount + 1).Y + mid.Y) - cutoffValue) / highestAbs;
                            offy1 = normal.Y * EDGE_DISTORTION_SCALE * len / 2.0f * (lineDistortionNoise.GetPerlinFractal(MathHelper.getSamplePoint(i, subdivisionCount + 1).X + mid.X, MathHelper.getSamplePoint(i, subdivisionCount + 1).Y + mid.Y) - cutoffValue) / highestAbs;
                            offx2 = normal.X * EDGE_DISTORTION_SCALE * len / 2.0f * (lineDistortionNoise.GetPerlinFractal(MathHelper.getSamplePoint(i + 1, subdivisionCount + 1).X + mid.X, MathHelper.getSamplePoint(i + 1, subdivisionCount + 1).Y + mid.Y) - cutoffValue) / highestAbs;
                            offy2 = normal.Y * EDGE_DISTORTION_SCALE * len / 2.0f * (lineDistortionNoise.GetPerlinFractal(MathHelper.getSamplePoint(i + 1, subdivisionCount + 1).X + mid.X, MathHelper.getSamplePoint(i + 1, subdivisionCount + 1).Y + mid.Y) - cutoffValue) / highestAbs;

                            // Place current, next and midpoint into vertex
                            subdividedCells.Append(new Vertex(new Vector2f((float)cell.X, (float)cell.Y), cellColor));
                            subdividedCells.Append(new Vertex(new Vector2f((float)line.Start.X + i * inc.X + offx1, (float)line.Start.Y + i * inc.Y + offy1), cellColor));
                            subdividedCells.Append(new Vertex(new Vector2f((float)line.Start.X + (i + 1) * inc.X + offx2, (float)line.Start.Y + (i + 1) * inc.Y + offy2), cellColor));
                        }

                    } else
                    {
                        subdividedCells.Append(new Vertex(new Vector2f((float)cell.X, (float)cell.Y), cellColor));
                        subdividedCells.Append(new Vertex(new Vector2f((float)line.Start.X, (float)line.Start.Y), cellColor));
                        subdividedCells.Append(new Vertex(new Vector2f((float)line.End.X, (float)line.End.Y), cellColor));
                    }
                }

                if (cellDictionary[cell].additionalEdges.Count > 0)
                {
                    foreach (var edge in cellDictionary[cell].additionalEdges)
                    {
                        cells.Append(new Vertex(new Vector2f((float)cell.X, (float)cell.Y), cellColor));
                        cells.Append(new Vertex(edge.Key, cellColor));
                        cells.Append(new Vertex(edge.Value, cellColor));

                        subdividedCells.Append(new Vertex(new Vector2f((float)cell.X, (float)cell.Y), cellColor));
                        subdividedCells.Append(new Vertex(edge.Key, cellColor));
                        subdividedCells.Append(new Vertex(edge.Value, cellColor));
                    }
                }

                // 
                //  TODO:
                //  Subdivide solid->passable transition edges and store them in another vertex array with copies of all non bordering edges for rendering comparison
                //
            }
            #endregion
        }

        private void OnClose(object sender, EventArgs e)
        {
            window.Close();
        }
    }
}
