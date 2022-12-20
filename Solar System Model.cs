using Rhino;
using Rhino.Geometry;
using Rhino.DocObjects;
using Rhino.Commands;
using Rhino.Input;
using Rhino.Geometry.Intersect;
using System.Collections.Generic;
using Rhino.UI;
using System;
using Rhino.Display;
using System.IO;
using System.Runtime.Remoting;
using ObjRef = Rhino.DocObjects.ObjRef;

/*
 *Things to be added before our model can become actually intresting :
There are a few ways to make the calculations in the provided code more accurate:

Use a higher precision for the gravitational constant: The value used for the gravitational constant (6.67 * 10^-11) 
is only accurate to two significant figures. A more accurate value would be 6.674184 * 10^-11.
Use a more accurate value for the mass of the sun: The mass of the sun used in the code (1.989 * 10^30 kilograms) 
is only accurate to two significant figures. A more accurate value would be 1.989 x 10^30 kilograms.
Use a more accurate value for the distance from the sun for each planet: The distance from the sun for each planet is given in meters, 
but the actual distances are much larger (on the order of 10^11 meters). This can lead to significant errors in the calculations. To improve the 
accuracy, you should use the average distance from the sun for each planet in astronomical units (AU), which is a unit of distance specifically 
defined to represent the distance from the sun to Earth (1 AU is equal to approximately 149.6 million kilometers).
Take into account the elliptical shape of planetary orbits: Planetary orbits are not perfectly circular, but are actually elliptical. 
This means that the distance from the sun changes over the course of the orbit. To account for this, you could use the average distance 
from the sun over the course of the orbit in your calculations. This average distance is known as the semi-major axis of the elliptical orbit.
Use more accurate values for the mass of each planet: The masses of the planets are only given to two significant figures in the provided code.
Using more accurate values for the masses of the planets would improve the accuracy of the calculations. You can find more accurate values 
for the masses of the planets in kilograms on various online resources such as NASA’s Planetary Fact Sheets (Planetary Fact Sheet).
Another set of tasks to accomplish (that you can commit on the repo)

The influence of other celestial bodies: The gravitational force on a planet is not only due to the sun, but also due to the presence of 
other celestial bodies such as other planets and moons. To account for these additional forces, you would need to consider the masses and 
distances of these bodies as well.
The general theory of relativity: The theory of relativity, developed by Albert Einstein, describes how the presence of matter and energy 
can affect the curvature of space-time. This theory can cause small deviations from the predictions of classical mechanics, such as the 
laws of gravitation and motion. To incorporate the effects of relativity, you would need to use more complex mathematical models such as 
the Schwarzschild metric or the Kerr metric.
The oblateness of celestial bodies: Celestial bodies are not perfectly spherical, but are slightly flattened at the poles due to their rotation.
This means that the gravitational force at the surface of a celestial body depends on the distance from the body’s center of mass and the latitude. 
To account for this effect, you could use a more accurate model of the gravitational field, such as the gravitational potential of a rotating, 
oblate spheroid.
The non-uniform distribution of mass within celestial bodies: The mass of celestial bodies is not evenly distributed, but can be concentrated in 
certain regions such as the core or mantle. This can affect the gravitational force experienced by objects at the surface. To account for this 
effect, you could use a more accurate model of the gravitational field, such as the gravitational potential of a body with a non-uniform mass 
distribution.
 *
 */
namespace My_Rhino_Commands
{
    public class SectionAllObjects : Command
    {
        public override string EnglishName { get { return "SectionAllObjects"; } }
        public class Planet
        {
            public string Name { get; set; }
            public Sphere Shape { get; set; }
            public double DistanceFromSun { get; set; }  // in meters
            public double Mass { get; set; }  // in kilograms
            public double Eccentricity { get; set; } //needed for accurate orbital calculations
            public double OrbitPeriod { get; private set; }
            public Guid Guid { get; set; }
            public Planet(string name, Sphere shape, double distanceFromSun, double mass, Guid guid, double eccentricity)
            {
                Name = name;
                Shape = shape;
                DistanceFromSun = distanceFromSun;
                Mass = mass;
                Guid = guid;
                Eccentricity = eccentricity;
                OrbitPeriod = 2 * Math.PI * Math.Sqrt(Math.Pow(DistanceFromSun, 3) / (Mass * 1.989 * Math.Pow(10, 30)));

            }
        }
        protected override Result RunCommand(RhinoDoc doc, RunMode mode)
        {
            //all uranus values are wrong to be fixxed with the manual  perigee apogee and SemiMajorAxes
            List<string> planetNames = new List<string>() { "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune" };
            List<double> distancesFromSun = new List<double>() { 5.79 * Math.Pow(10, 7), 1.08 * Math.Pow(10, 8), 1.50 * Math.Pow(10, 8), 2.28 * Math.Pow(10, 8), 7.78 * Math.Pow(10, 8), 1.43 * Math.Pow(10, 9), 2.87 * Math.Pow(10, 9), 4.50 * Math.Pow(10, 9) };
            List<double> masses = new List<double>() { 3.30 * Math.Pow(10, 23), 4.87 * Math.Pow(10, 24), 5.97 * Math.Pow(10, 24), 6.42 * Math.Pow(10, 23), 1.90 * Math.Pow(10, 27), 5.68 * Math.Pow(10, 26), 8.68 * Math.Pow(10, 25), 1.02 * Math.Pow(10, 26) };
            List<double> perigees = new List<double>() { 4.092e10, 7.232e10, 1.471e11, 1.383e11, 4.220e11, 9.050e10, 2.735e11, 4.506e11 };
            List<double> apogees = new List<double>() { 4.794e10, 7.850e10, 1.521e11, 1.666e11, 8.157e11, 1.352e11, 3.004e11, 4.600e11 }; 
            List<double> semiMajorAxes = new List<double>() { 5.792e10, 1.082e11, 1.496e11, 2.279e11, 7.785e11, 1.429e12, 2.870e12, 4.498e12 };
            List<double> semiMinorAxes = new List<double>() { 0.38709893,0.72333199,0.9932373,1.523666231,5.20336301,9.53707032,19.19126393,30.06896348 };
            List<Planet> planets = DrawPlanets(planetNames, distancesFromSun, masses,  perigees, apogees, semiMajorAxes,semiMinorAxes);
            double sunMass = 1.989 * Math.Pow(10, 30);  // in kilograms
            CalculateOrbits(planets, sunMass);
            DrawOrbits(planets, sunMass,semiMajorAxes,semiMinorAxes);
            AddMaterialsToPlanets(planets);

            return Result.Success;
        }
        public static double CalculateOrbitPeriod(double distanceFromSun, double mass, double sunMass)
        {
            return 2 * Math.PI * Math.Sqrt(Math.Pow(distanceFromSun, 3) / (mass * sunMass));
        }
        public static List<Planet> DrawPlanets(List<string> planetNames, List<double> distancesFromSun,
            List<double> masses,List<double> perigee, List<double> apogee,List<double> semiMajorAxes,List<double> semiMinorAxes)
        {
            // Set the display mode to wireframe
            RhinoDoc.ActiveDoc.Views.RedrawEnabled = false;
            RhinoDoc.ActiveDoc.Layers.SetCurrentLayerIndex(0, true);
            RhinoDoc.ActiveDoc.Views.RedrawEnabled = true;

            // Create a list of planets
            List<Planet> planets = new List<Planet>();

            // Draw each planet
            for (int i = 0; i < planetNames.Count; i++)
            {
                // Calculate the position of the planet in its orbit
                double r = distancesFromSun[i];
                double e = GetEccentricity(semiMajorAxes[i], perigee[i], apogee[i]);
                double a = r / (1 - e);
                double b = Math.Sqrt(a * a * (1 - e * e));
                double T = CalculateOrbitPeriod(distancesFromSun[i], masses[i], 1.989e30);
                double M = 2 * Math.PI * 25000 / T;
                double E = M + e * Math.Sin(M);
                double x = a * Math.Cos(E);
                double y = b * Math.Sin(E);
                double z = 0;  // Initialize the Z coordinate to 0

                // Modify the Z coordinate to take into account the orbit's inclination (angle of orbit relative to the XY plane)
                double inclination = 45;  // Replace this with the actual inclination of the orbit
                z = Math.Sin(inclination) * y;  // Modify the Z coordinate based on the inclination of the orbit

                // Create a sphere for the planet
                double radius = masses[i] / (4.0 * Math.PI * Math.Pow(r, 3) / 3.0);
                Sphere sphere = new Sphere(new Point3d(x, y, z), radius);

                // Add the sphere to the document and get its Guid
                Guid guid = RhinoDoc.ActiveDoc.Objects.AddSphere(sphere);

                // Create a planet object and add it to the list
                Planet planet = new Planet(planetNames[i], sphere, r, masses[i], guid, e);
                planets.Add(planet);

                // Add a text dot for the planet
                Point3d textDotLocation = new Point3d(x, y, z);
                TextDot textDot = new TextDot(planetNames[i], textDotLocation);
                RhinoDoc.ActiveDoc.Objects.AddTextDot(textDot);
            }

            // Redraw the view
            RhinoDoc.ActiveDoc.Views.Redraw();

            return planets;
        }

        private static double GetEccentricity(double semiMajorAxes, double perigee, double apogee)
        {
            return (apogee - perigee) / (2 * semiMajorAxes);
        }

        public static void CalculateOrbits(List<Planet> planets, double sunMass)
        {
            /*Eccentricity is a measure of how much an orbit deviates from a perfect circle. 
            For the orbital velocity:

            v = sqrt(GM(2/r - 1/a))
            v = sqrt(GM(1 + ecos(E)) / (r * (1 - e^2)))
            For the orbital period:

            T = 2 * Pi * sqrt(a^3 / GM)
            T = 2 * Pi * sqrt((a^3 / GM) / (1 - e^2))
            Where:

            v: orbital velocity
            T: orbital period
            G: gravitational constant (6.67 * 10^-11 N * m^2 / kg^2)
            M: mass of the sun (in kilograms)
            r: distance of the planet from the sun (in meters)
            a: semi-major axis of the orbit (in meters)
            e: eccentricity of the orbit
            E: eccentric anomaly
            */
            // Calculate the orbits of the planets


            //ok let's integrate the eliptical shape first
            //we need the semi-major axis, semi minor axis, the focus(F), the eccentric anomaly (E) now we
            foreach (Planet planet in planets)
            {
                // Calculate the gravitational force between the sun and the planet
                double force = (6.67 * Math.Pow(10, -11) * sunMass * planet.Mass) / Math.Pow(planet.DistanceFromSun, 2);

                // Calculate the orbital velocity of the planet
                double velocity = Math.Sqrt(force / planet.Mass);

                // Calculate the orbital period of the planet
                double period = 2 * Math.PI * planet.DistanceFromSun / velocity;

                // Print the results
                RhinoApp.WriteLine($"Planet: {planet.Name}");
                RhinoApp.WriteLine($"  Distance from sun: {planet.DistanceFromSun} meters");
                RhinoApp.WriteLine($"  Mass: {planet.Mass} kilograms");
                RhinoApp.WriteLine($"  Orbital velocity: {velocity} meters/second");
                RhinoApp.WriteLine($"  Orbital period: {period} seconds");
            }
        }




        //getEllipseParameters won't work for long distances over 10ly, to fix the formulas,
        //it doesnt account for the gravitational field of close celestial bodies
        private static void GetEllipseParameters(Planet planet, out double a, out double b, out double F)
        {
            // Calculate the semi-major axis of the ellipse
            a = planet.DistanceFromSun / (1 - planet.Eccentricity);
                                                                    

            // Calculate the semi-minor axis of the ellipse
            b = Math.Sqrt(a * a * (1 - planet.Eccentricity * planet.Eccentricity));
            // Calculate the focus of the ellipse
            F = a * planet.Eccentricity;
        }

        public static void DrawOrbits(List<Planet> planets, double sunMass,List<double> semiMajorAxis,List<double>semiMinorAxis)
        {
            int i = 0;
            // Calculate the orbits of the planets
            foreach (Planet planet in planets)
            {
                // Calculate the gravitational force between the sun and the planet
                double force = (6.67 * Math.Pow(10, -11) * sunMass * planet.Mass) / Math.Pow(planet.DistanceFromSun, 2);

                // Calculate the orbital velocity of the planet
                double velocity = Math.Sqrt(force / planet.Mass);

                // Calculate the orbital period of the planet
                double period = 2 * Math.PI * planet.DistanceFromSun / velocity;

                // Calculate the radius of the orbit
                double orbitRadius = planet.DistanceFromSun;

                //We must get A,B,F values, which rappresent the properties of an eliptical shape
                double a, b, f;
                f = GetFocus(semiMajorAxis[i], semiMinorAxis[i]);
                GetEllipseParameters(planet, out a, out b, out f); //F = √(a^2 - b^2)

                Curve orbit = Curve.CreateInterpolatedCurve(GetOrbitPoints(a, b, f), 3);
                //make a method to determine the degree of the curve based on the eleptical value


                RhinoDoc.ActiveDoc.Objects.AddCurve(orbit);
                i++;
            }

            // Redraw the view
            RhinoDoc.ActiveDoc.Views.Redraw();
        }

        private static double GetFocus(double semiMajorAxis, double semiMinorAxis)
        {
            // Calculate the eccentricity of the ellipse
            double eccentricity = Math.Sqrt(semiMajorAxis * semiMajorAxis - semiMinorAxis * semiMinorAxis) / semiMajorAxis;

            // Calculate the focus of the ellipse
            double focus = Math.Sqrt(semiMajorAxis * semiMajorAxis - semiMinorAxis * semiMinorAxis);

            return focus;
        }

       
        private static Point3d[] GetOrbitPoints(double a, double b, double F)
        {
            /*Determine the semi-major axis (a) and the semi-minor axis (b) of the ellipse. These are the lengths of the major and
             minor axes of the ellipse, respectively, and are related to the distance from the sun and the eccentricity of the orbit.
            Determine the focus (F) of the ellipse. This is a point within the ellipse that is closer to the sun than the center of the ellipse.
            Determine the eccentric anomaly (E) at each point along the orbit. This is a measure of the angle between the focus and 
            the position of the planet at a given point in time.
            */
            // Set the number of points in the orbit
            int pointCount = 360;

            // Create an array of points
            Point3d[] points = new Point3d[pointCount];

            // Generate points along the orbit
            for (int i = 0; i < pointCount; i++)
            {
                double E = i * 2 * Math.PI / pointCount;
                double x = a * Math.Cos(E) - F;
                double y = b * Math.Sin(E);
                points[i] = new Point3d(x, y, 0);
            }

            return points;
        }
        public static void AddMaterialsToPlanets(List<Planet> planets) 
        {
            // Set up a dictionary of material names and colors for the planets
            Dictionary<string, Color4f> colors = new Dictionary<string, Color4f>()
            {
                { "Mercury", new Color4f(0.75f, 0.75f, 0.75f, 1.0f) },
                { "Venus", new Color4f(1.0f, 1.0f, 0.88f, 1.0f) },
                { "Earth", new Color4f(0.125f, 0.5f, 1.0f, 1.0f) },
                { "Mars", new Color4f(1.0f, 0.5f, 0.25f, 1.0f) },
                { "Jupiter", new Color4f(1.0f, 0.84f, 0.0f, 1.0f) },
                { "Saturn", new Color4f(1.0f, 0.75f, 0.0f, 1.0f) },
                { "Uranus", new Color4f(0.0f, 0.5f, 1.0f, 1.0f) },
                { "Neptune", new Color4f(0.0f, 0.5f, 1.0f, 1.0f) },
                { "Pluto", new Color4f(0.75f, 0.75f, 0.75f, 1.0f) }
            };


            // Set the material of each planet
            // Set the color of each planet
            foreach (Planet planet in planets)
            {
                if (colors.ContainsKey(planet.Name))
                {
                    // Get the object representing the planet
                    RhinoObject obj = RhinoDoc.ActiveDoc.Objects.Find(planet.Guid);
                    Color4f thecolor = colors[planet.Name];
                    // Set the object color
                    obj.Attributes.ObjectColor = thecolor.AsSystemColor();
                }
            }
        }
        public  void  CalculateCurrentPositions(Planet[] planets, double sunMass) //to account for eccentricity and orbit elasticity 
        {
            // Loop through the list of planets
            foreach (Planet planet in planets)
            {
                // Get the current date and time
                DateTime currentTime = DateTime.Now;

                // Calculate the gravitational force between the sun and the planet
                double force = (6.67 * Math.Pow(10, -11) * sunMass * planet.Mass) / Math.Pow(planet.DistanceFromSun, 2);

                // Calculate the orbital velocity of the planet
                double velocity = Math.Sqrt(force / planet.Mass);

                // Calculate the time elapsed since the start of the year
                TimeSpan elapsedTime = currentTime - new DateTime(currentTime.Year, 1, 1);

                // Calculate the current angle of the planet around its orbit
                double angle = velocity * elapsedTime.TotalSeconds / (2 * Math.PI * planet.DistanceFromSun);

                // Calculate the current position of the planet
                Point3d currentPosition = new Point3d(planet.DistanceFromSun * Math.Cos(angle), planet.DistanceFromSun * Math.Sin(angle), 0);

                // Find the RhinoObject associated with the planet's Guid
                RhinoObject obj = RhinoDoc.ActiveDoc.Objects.Find(planet.Guid);

                // Modify the geometry of the object
                Sphere sphere = new Sphere(currentPosition, planet.Shape.Radius);
                Brep brepSphere = Brep.CreateFromSphere(sphere);
    
                // Get a reference to the underlying geometry object
                ObjRef planetRef = new ObjRef(planet.Guid);
                Guid sphereGuid=RhinoDoc.ActiveDoc.Objects.Add(brepSphere);
                RhinoDoc.ActiveDoc.Objects.Replace(planetRef, RhinoDoc.ActiveDoc.Objects.Find(sphereGuid));

                // Redraw the view
                RhinoDoc.ActiveDoc.Views.Redraw();

            }
        }
    } 
}
