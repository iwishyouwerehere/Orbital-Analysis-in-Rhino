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
            public Guid Guid { get; set; }
            public Planet(string name, Sphere shape, double distanceFromSun, double mass, Guid guid)
            {
                Name = name;
                Shape = shape;
                DistanceFromSun = distanceFromSun;
                Mass = mass;
                Guid = guid;
            }
        }
        protected override Result RunCommand(RhinoDoc doc, RunMode mode)
        {
            List<string> planetNames = new List<string>() { "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune" };
            List<double> distancesFromSun = new List<double>() { 5.79 * Math.Pow(10, 7), 1.08 * Math.Pow(10, 8), 1.50 * Math.Pow(10, 8), 2.28 * Math.Pow(10, 8), 7.78 * Math.Pow(10, 8), 1.43 * Math.Pow(10, 9), 2.87 * Math.Pow(10, 9), 4.50 * Math.Pow(10, 9) };
            List<double> masses = new List<double>() { 3.30 * Math.Pow(10, 23), 4.87 * Math.Pow(10, 24), 5.97 * Math.Pow(10, 24), 6.42 * Math.Pow(10, 23), 1.90 * Math.Pow(10, 27), 5.68 * Math.Pow(10, 26), 8.68 * Math.Pow(10, 25), 1.02 * Math.Pow(10, 26) };

            List<Planet> planets = DrawPlanets(planetNames, distancesFromSun, masses);
            double sunMass = 1.989 * Math.Pow(10, 30);  // in kilograms
            CalculateOrbits(planets, sunMass);
            DrawOrbits(planets, sunMass);
            AddMaterialsToPlanets(planets);

            return Result.Success;
        }

        public static List<Planet> DrawPlanets(List<string> planetNames, List<double> distancesFromSun,
            List<double> masses)
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
                // Create a sphere for the planet
                double radius = masses[i] / (4.0 * Math.PI * Math.Pow(distancesFromSun[i], 3) / 3.0);
                Sphere sphere = new Sphere(Point3d.Origin, radius);

                // Add the sphere to the document and get its Guid
                Guid guid=RhinoDoc.ActiveDoc.Objects.AddSphere(sphere);

                // Create a planet object and add it to the list
                Planet planet = new Planet(planetNames[i], sphere, distancesFromSun[i], masses[i], guid);
                planets.Add(planet);
            }

            // Redraw the view
            RhinoDoc.ActiveDoc.Views.Redraw();

            return planets;
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
        public static void DrawOrbits(List<Planet> planets, double sunMass)
        {
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

                // Create a curve for the orbit
                Curve orbit = Curve.CreateInterpolatedCurve(GetOrbitPoints(orbitRadius), 3);

                // Add the curve to the document
                RhinoDoc.ActiveDoc.Objects.AddCurve(orbit);
            }

            // Redraw the view
            RhinoDoc.ActiveDoc.Views.Redraw();
        }
        private static Point3d[] GetOrbitPoints(double orbitRadius)
        {
            // Set the number of points in the orbit
            int pointCount = 360;

            // Create an array of points
            Point3d[] points = new Point3d[pointCount];

            // Generate points along the orbit
            for (int i = 0; i < pointCount; i++)
            {
                double angle = i * 2 * Math.PI / pointCount;
                double x = orbitRadius * Math.Cos(angle);
                double y = orbitRadius * Math.Sin(angle);
                points[i] = new Point3d(x, y, 0);
            }

            return points;
        }
        public static void AddMaterialsToPlanets(List<Planet> planets) //add planets here
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
                // Get a reference to the underlying geometry object


                // Redraw the view
                RhinoDoc.ActiveDoc.Views.Redraw();
            }
        }
    }
}