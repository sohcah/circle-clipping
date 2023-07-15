import './App.css'
import MapGL, {Layer, Source} from 'react-map-gl/maplibre';
import {
    circle,
    destination,
    Feature, point,

} from '@turf/turf';
import seedrandom from 'seed-random';
import KDBush from "kdbush";
import * as geokdbush from "geokdbush-tk";

const pointsPerCircle = 32;

const data: [number, number, number][] = [
    [-122.444768, 37.782551, 0.03],
    [-122.444586, 37.782745, 0.035],
    [-122.443688, 37.782842, 0.04],
    [-122.442715, 37.782932, 0.04],
];

const r = seedrandom("hi");
//
for (let i = 0; i < 0; i++) {
    data.push([data[0][0] + r() * 0.05, data[0][1] + r() * 0.05, 0.03]);
}
console.log(r());

interface IntersectionResult {
    points: [number, number][];
    midpoint: [number, number];
}

type CircleArray = [number, number, number];
const earthRadius = 6371; // Radius of the Earth in kilometers

function calculateIntersection(circle1: CircleArray, circle2: CircleArray): IntersectionResult | null {

    // Convert the coordinates to radians
    const lat1 = toRadians(circle1[1]);
    const lon1 = toRadians(circle1[0]);
    const lat2 = toRadians(circle2[1]);
    const lon2 = toRadians(circle2[0]);

    // Calculate the distance between the circle centers
    const distance = haversine(lat1, lon1, lat2, lon2, earthRadius);

    // Check if the circles intersect
    if (distance > circle1[2] + circle2[2] || distance < Math.abs(circle1[2] - circle2[2])) {
        // No intersection points or circles are contained within each other
        return null;
    }

    // Calculate the intersection angles
    const intersectionAngle = Math.acos(
        (circle1[2] * circle1[2] + distance * distance - circle2[2] * circle2[2]) /
        (2 * circle1[2] * distance)
    );

    // Calculate the coordinates of the intersection points
    const angle1 = Math.atan2(lon2 - lon1, Math.log(Math.tan(lat2 / 2 + Math.PI / 4) / Math.tan(lat1 / 2 + Math.PI / 4)));
    const angle2 = angle1 + intersectionAngle;
    const angle3 = angle1 - intersectionAngle;

    const p1lat = toDegrees(Math.asin(Math.sin(lat1) * Math.cos(circle1[2] / earthRadius) +
        Math.cos(lat1) * Math.sin(circle1[2] / earthRadius) * Math.cos(angle2)));
    const p1lng = toDegrees(lon1 + Math.atan2(Math.sin(angle2) * Math.sin(circle1[2] / earthRadius) * Math.cos(lat1),
        Math.cos(circle1[2] / earthRadius) - Math.sin(lat1) * Math.sin(toRadians(p1lat))))

    const p2lat = toDegrees(Math.asin(Math.sin(lat1) * Math.cos(circle1[2] / earthRadius) +
        Math.cos(lat1) * Math.sin(circle1[2] / earthRadius) * Math.cos(angle3)));
    const p2lng = toDegrees(lon1 + Math.atan2(Math.sin(angle3) * Math.sin(circle1[2] / earthRadius) * Math.cos(lat1),
        Math.cos(circle1[2] / earthRadius) - Math.sin(lat1) * Math.sin(toRadians(p2lat))));

    return {
        points: [
            [p1lng, p1lat],
            [p2lng, p2lat]
        ], midpoint: [
            (p1lng + p2lng) / 2,
            (p1lat + p2lat) / 2
        ]
    };
}

// Helper function to convert degrees to radians
function toRadians(degrees: number): number {
    return degrees * (Math.PI / 180);
}

// Helper function to convert radians to degrees
function toDegrees(radians: number): number {
    return radians * (180 / Math.PI);
}

// Haversine formula to calculate the great-circle distance between two points
function haversine(lat1: number, lon1: number, lat2: number, lon2: number, radius: number): number {
    const dlat = lat2 - lat1;
    const dlon = lon2 - lon1;
    const a =
        Math.sin(dlat / 2) * Math.sin(dlat / 2) +
        Math.cos(lat1) * Math.cos(lat2) * Math.sin(dlon / 2) * Math.sin(dlon / 2);
    const c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
    return radius * c;
}

function generateFlattened(circles: CircleArray[]): Feature[] {
    const features: Feature[] = [];

// Benchmark the performance of calculating the intersection points
    const start = performance.now();

    const kdbush = new KDBush(circles.length);

    for (const circle of circles) {
        kdbush.add(circle[0], circle[1]);
    }

    kdbush.finish();

    const intersectionPointsByCircle: [number, number][][] = new Array(circles.length).fill(null!).map(() => [])

    const maximumRadius = circles.reduce((max, c) => c[2] > max ? c[2] : max, 0);

    for (let i = 0; i < circles.length; i++) {
        const circle = circles[i];
        const potentialIntersections = geokdbush.around(kdbush, circle[0], circle[1], Infinity, maximumRadius + circle[2], c => c > i);
        for (const potentialIntersection of potentialIntersections) {
            const otherCircle = circles[potentialIntersection];
            const intersection = calculateIntersection(
                circle,
                otherCircle,
            );

            if (intersection) {
                for (const pt of intersection.points) {
                    features.push(point(pt, {
                        radius: 10,
                    }));
                    intersectionPointsByCircle[i].push(pt);
                    intersectionPointsByCircle[potentialIntersection].push(pt);
                }
            }
        }
    }

    const circleArcs: [number, number][][] = new Array(circles.length).fill(null!).map(() => []);

    for (let i = 0; i < circles.length; i++) {
        const circle = circles[i];
        const intersectionPoints = intersectionPointsByCircle[i];

        if (intersectionPoints.length === 0) {
            circleArcs[i].push([0, 2 * Math.PI]);
            continue;
        }

        const sortedIntersectionPoints = intersectionPoints.sort((a, b) => {
            const angleA = Math.atan2(a[1] - circle[1], a[0] - circle[0]);
            const angleB = Math.atan2(b[1] - circle[1], b[0] - circle[0]);
            return angleA - angleB;
        });

        let prevAngle = Math.atan2(sortedIntersectionPoints[sortedIntersectionPoints.length - 1][1] - circle[1], sortedIntersectionPoints[sortedIntersectionPoints.length - 1][0] - circle[0]);

        for (const pt of sortedIntersectionPoints) {
            const angle = Math.atan2(pt[1] - circle[1], pt[0] - circle[0]);
            circleArcs[i].push([(prevAngle + 2 * Math.PI) % (2 * Math.PI), (angle + 2 * Math.PI) % (2 * Math.PI)]);
            prevAngle = angle;
        }
    }

    console.log(circleArcs);

    // Create line segments for each circle arc with 10 points per arc
    for (let i = 0; i < circles.length; i++) {
        const circle = circles[i];
        const arcs = circleArcs[i];

        for (const arc of arcs) {
            let [startAngle, endAngle] = arc;
            if (startAngle > endAngle) {
                endAngle += 2 * Math.PI;
            }
            const noPoints = Math.floor(Math.max(10, Math.abs(endAngle - startAngle) * 10));
            const deltaAngle = (endAngle - startAngle) / noPoints;

            const points: [number, number][] = [];
            for (let p = 0; p <= noPoints; p++) {
                const angle = startAngle + p * deltaAngle;
                const pt = destination([circle[0], circle[1]], circle[2], angle * 180 / Math.PI, {units: "kilometers"});
                points.push(pt.geometry.coordinates);
            }

            features.push({
                type: "Feature",
                properties: {
                    color: `#${Math.floor(Math.random() * 0x1000000).toString(16).padStart(6, "0")}`
                },
                geometry: {
                    type: "LineString",
                    coordinates: points
                }
            });
        }
    }

    const end = performance.now();
    console.log(`Time taken: ${end - start}ms`);

    return features;
}


function App() {
    return (
        <>
            <div style={{
                display: "flex",
                height: "100vh",
                width: "100vw",
            }}>
                {/*<div style={{flex: 1}}>*/}
                {/*    <MapGL*/}
                {/*        initialViewState={{*/}
                {/*            longitude: -122.4,*/}
                {/*            latitude: 37.8,*/}
                {/*            zoom: 14*/}
                {/*        }}*/}
                {/*        mapStyle="https://api.maptiler.com/maps/streets/style.json?key=k02Kl04VqEmUthdcLti3"*/}
                {/*    >*/}
                {/*        <Source type="geojson" data={{*/}
                {/*            type: "FeatureCollection",*/}
                {/*            features: (() => {*/}
                {/*                const circles = data.map((d) => circle([d[0], d[1]], d[2], {*/}
                {/*                    steps: pointsPerCircle,*/}
                {/*                }))*/}

                {/*                return union(circles.map(i => i.geometry.coordinates as any)).map(i => ({*/}
                {/*                    type: "Feature" as const,*/}
                {/*                    geometry: {*/}
                {/*                        type: "Polygon" as const,*/}
                {/*                        properties: {},*/}
                {/*                        coordinates: i,*/}
                {/*                    },*/}
                {/*                }));*/}
                {/*            })()*/}
                {/*        }}>*/}
                {/*            <Layer type="fill" paint={{*/}
                {/*                "fill-color": "#007cbf",*/}
                {/*                "fill-opacity": 1*/}
                {/*            }}/>*/}
                {/*            <Layer type="line" paint={{*/}
                {/*                "line-color": "#ad0000",*/}
                {/*                "line-width": 2*/}
                {/*            }}/>*/}

                {/*            /!*<Layer type="circle" paint={{*!/*/}
                {/*            /!*    "circle-radius": 4,*!/*/}
                {/*            /!*    "circle-color": "#007cbf"*!/*/}
                {/*            /!*}}/>*!/*/}
                {/*        </Source>*/}
                {/*    </MapGL>*/}
                {/*</div>*/}
                <div style={{flex: 1}}>
                    <MapGL
                        initialViewState={{
                            longitude: -122.4,
                            latitude: 37.8,
                            zoom: 14
                        }}
                        mapStyle="https://api.maptiler.com/maps/streets/style.json?key=k02Kl04VqEmUthdcLti3"
                    >
                        <Source type="geojson" data={{
                            type: "FeatureCollection",
                            features: data.map((d) => circle([d[0], d[1]], d[2], {
                                steps: pointsPerCircle,
                                units: "kilometers",
                            }))
                        }}>
                            <Layer type="fill" paint={{
                                "fill-color": "#007cbf",
                                "fill-opacity": 0.1
                            }}/>
                        </Source>
                        {/*<Source type="geojson" data={{*/}
                        {/*    type: "FeatureCollection",*/}
                        {/*    features: flatten(data)*/}
                        {/*    // .map(d => d.length < 2 ? point(d[0]) : d.length < 4 ? lineString(d) : polygon(d))*/}
                        {/*}}>*/}
                        {/*    /!*<Layer type="line" paint={{*!/*/}
                        {/*    /!*    "line-color": "#ad0000",*!/*/}
                        {/*    /!*    "line-width": 2*!/*/}
                        {/*    /!*}}/>*!/*/}

                        {/*    /!*<Layer type="circle" paint={{*!/*/}
                        {/*    /!*    "circle-radius": 4,*!/*/}
                        {/*    /!*    "circle-color": "#007cbf"*!/*/}
                        {/*    /!*}}/>*!/*/}
                        {/*    <Layer type="circle" paint={{*/}
                        {/*        "circle-radius": 3,*/}
                        {/*        "circle-color": ["get", "color"],*/}
                        {/*        "circle-opacity": ["get", "opacity"]*/}
                        {/*    }}/>*/}
                        {/*    <Layer type="line" paint={{*/}
                        {/*        "line-color": ["get", "color"],*/}
                        {/*        "line-width": 4,*/}
                        {/*        "line-opacity": ["get", "opacity"]*/}
                        {/*    }}/>*/}

                        {/*    /!*<Layer type="fill" paint={{*!/*/}
                        {/*    /!*    "fill-color": ["get", "color"],*!/*/}
                        {/*    /!*    "fill-opacity": 0.5*!/*/}
                        {/*    /!*}}/>*!/*/}
                        {/*</Source>*/}
                        <Source type="geojson" data={{
                            type: "FeatureCollection",
                            features: generateFlattened(data)
                            // .map(d => d.length < 2 ? point(d[0]) : d.length < 4 ? lineString(d) : polygon(d))
                        }}>
                            <Layer type="line"
                                   paint={{
                                       "line-color": ["get", "color"],
                                       "line-width": 2
                                   }}/>
                            <Layer type="circle"
                                   paint={{
                                       "circle-color": ["get", "color"],
                                       "circle-radius": ["get", "radius"],
                                   }}/>
                        </Source>
                    </MapGL>
                </div>
            </div>
        </>
    )
}

export default App
