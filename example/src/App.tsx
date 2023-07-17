import './App.css'
import MapGL, {Layer, Source} from 'react-map-gl';
import {
    destination,
    Feature,
    bearing,
    // distance,
    polygon, area, circle,
} from "@turf/turf";
import seedrandom from 'seed-random';
import KDBush from "kdbush";
import * as geokdbush from "geokdbush-tk";
import CheapRuler from "cheap-ruler";
import Clipping from "polygon-clipping";
const {union} = Clipping;

const rulers = new Map<number, CheapRuler>();

function getRuler(lat: number): CheapRuler {
    const key = Math.floor(lat * 10);
    let ruler = rulers.get(key);
    if (!ruler) {
        ruler = new CheapRuler(lat, "kilometers");
        rulers.set(key, ruler);
    }
    return ruler;
}

function distance(c1: [number, number], c2: [number, number]): number {
    return getRuler(c1[1]).distance(c1, c2);
    // return haversine(toRadians(lat1), toRadians(lon1), toRadians(lat2), toRadians(lon2), earthRadius);
}

const data: [number, number, number][] = [
    [-122.444768, 37.782551, 0.03],
    [-122.444586, 37.782745, 0.035],
    [-122.443688, 37.782842, 0.04],
    [-122.442715, 37.782932, 0.04],
];

const r = seedrandom("hi");

// 9720: kinda broken
// 9750: broken
// for (let i = 0; i < 5000; i++) {
//     data.push([data[0][0] + r() * 0.05, data[0][1] + r() * 0.05, r() * 0.1]);
// }
for (let i = 0; i < 5000; i++) {
    data.push([data[0][0] + r() * 0.05, data[0][1] + r() * 0.05, 0.03]);
}

interface IntersectionResult {
    points: [number, number][];
    midpoint: [number, number];
}

type CircleArray = [number, number, number];
const earthRadius = 6371.0088; // Radius of the Earth in kilometers

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
    const startClipping = performance.now();

    const circlePolys = circles.map(c => circle(c, c[2], {units: "kilometers"}));


    union(circlePolys.map(i => i.geometry.coordinates as any)).map(i => ({
        type: "Feature" as const,
        properties: {},
        geometry: {
            type: "Polygon" as const,
            coordinates: i,
        },
    }));

    const endClipping = performance.now();

    console.log(`Clipping took ${endClipping - startClipping}ms`);


    const features: Feature[] = [];

    // Benchmark the performance of calculating the intersection point
    const start = performance.now();

    const kdbush = new KDBush(circles.length);

    for (const circle of circles) {
        kdbush.add(circle[0], circle[1]);
    }

    kdbush.finish();

    const intersectionPointsByCircle: [number, number][][] = new Array(circles.length).fill(null!).map(() => [])

    const maximumRadius = circles.reduce((max, c) => c[2] > max ? c[2] : max, 0);

    const circleGroups: number[] = new Array(circles.length).fill(0).map((_, i) => i);

    for (let i = 0; i < circles.length; i++) {
        const circle = circles[i];
        const potentialIntersections = geokdbush.around(kdbush, circle[0], circle[1], Infinity, maximumRadius + circle[2], c => c !== i);
        for (const potentialIntersection of potentialIntersections) {
            const otherCircle = circles[potentialIntersection];

            if (circle[2] + otherCircle[2] > distance([circle[0], circle[1]], [otherCircle[0], otherCircle[1]])) {
                if (circleGroups[potentialIntersection] === undefined) {
                    circleGroups[potentialIntersection] = circleGroups[i];
                } else {
                    const groupToUpdate = circleGroups[potentialIntersection];
                    for (let j = 0; j < circles.length; j++) {
                        if (circleGroups[j] === groupToUpdate) {
                            circleGroups[j] = i;
                        }
                    }
                }
            }

            if (potentialIntersection < i) continue;

            const intersection = calculateIntersection(
                circle,
                otherCircle,
            );

            if (intersection) {
                for (const pt of intersection.points) {
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

        const angles = intersectionPoints.map(pt => bearing(pt, circle));

        angles.sort((a, b) => a - b);

        for (let j = 0; j < angles.length; j++) {
            const bearing1 = angles[j]
            const bearing2 = angles[(j + 1) % angles.length];

            let startAngle = bearing1 * Math.PI / 180;
            let endAngle = bearing2 * Math.PI / 180;

            if (startAngle < 0) {
                startAngle += 2 * Math.PI;
            }

            if (endAngle < 0) {
                endAngle += 2 * Math.PI;
            }

            if (endAngle < startAngle) {
                endAngle += 2 * Math.PI;
            }

            circleArcs[i].push([startAngle, endAngle]);

        }
    }

    const arcSegments: ([number, [number, number][]] | null)[] = [];

    // Create line segments for each circle arc with 10 points per arc
    for (let i = 0; i < circles.length; i++) {
        const circle = circles[i];
        const arcs = circleArcs[i];

        for (const arc of arcs) {
            let [startAngle, endAngle] = arc;
            if (startAngle > endAngle) {
                endAngle += 2 * Math.PI;
            }
            const noPoints = Math.floor(Math.max(5, Math.abs(endAngle - startAngle) * 10));
            const deltaAngle = (endAngle - startAngle) / noPoints;

            const points: [number, number][] = [];
            for (let p = 0; p <= noPoints; p++) {
                const angle = startAngle + p * deltaAngle;
                const pt = destination([circle[0], circle[1]], circle[2], angle * 180 / Math.PI + 180, {units: "kilometers"});
                points.push(pt.geometry.coordinates);
            }

            arcSegments.push([circleGroups[i], points]);
        }
    }

    const arcSegmentsKDBush = new KDBush(arcSegments.length);

    for (let i = 0; i < arcSegments.length; i++) {
        const [, segment] = arcSegments[i]!;
        const middleCoordinate = segment[Math.floor(segment.length / 2)];
        arcSegmentsKDBush.add(middleCoordinate[0], middleCoordinate[1]);
    }

    arcSegmentsKDBush.finish();

    for (const circle of circles) {
        const overlappingArcs = geokdbush.around(arcSegmentsKDBush, circle[0], circle[1], Infinity, circle[2] - 0.0001);
        for (const overlappingArc of overlappingArcs) {
            arcSegments[overlappingArc] = null;
        }
    }


    const filteredSegmentsByGroup = arcSegments.filter(s => s !== null)
        .reduce((acc, segment) => {
            const [group, arc] = segment!;
            acc[group] ??= [];
            acc[group].push(arc!);
            return acc;
        }, {} as {
            [key: number]: [number, number][][]
        });

    // Final step: Stitch together groups of line segments that are connected to make polygons
    const polygons: [number, number][][][] = [];

    for (const group in filteredSegmentsByGroup) {
        const segments = filteredSegmentsByGroup[group];
        const groupPolygon: [number, number][][] = [];
        let currentPolygon: [number, number][] = [];

        let i = 0;
        while (segments.length) {
            if (i++ > 1_000_000) {
                throw new Error("Infinite loop");
            }

            if (currentPolygon.length === 0) {
                currentPolygon.push(...segments[0]);
                segments.splice(0, 1);
                continue;
            }

            let found = false;
            for (const segment of segments) {
                const start = segment[0];
                const end = segment[segment.length - 1];

                const startDistance = distance(currentPolygon[currentPolygon.length - 1], start);
                const endDistance = distance(currentPolygon[0], end);

                if (startDistance < 0.0005) {
                    currentPolygon.push(...segment);
                    found = true;
                    segments.splice(segments.indexOf(segment), 1);
                    break;
                } else if (endDistance < 0.0005) {
                    currentPolygon.unshift(...segment);
                    found = true;
                    segments.splice(segments.indexOf(segment), 1);
                    break;
                }
            }
            if (currentPolygon.length > 2 && distance(currentPolygon[0], currentPolygon[currentPolygon.length - 1]) < 0.0005) {
                groupPolygon.push([...currentPolygon, currentPolygon[0]]);
                currentPolygon = [];
            } else if (!found) {
                if (currentPolygon.length > 2 && distance(currentPolygon[0], currentPolygon[currentPolygon.length - 1]) < 0.05) {
                    groupPolygon.push([...currentPolygon, currentPolygon[0]]);
                    currentPolygon = [];
                } else if (currentPolygon.length) {
                    console.error("Could not find a segment to connect to", currentPolygon, segments, distance(currentPolygon[0], currentPolygon[currentPolygon.length - 1]));
                    throw new Error("Could not find a segment to connect to");
                }
            }
        }

        if (currentPolygon.length > 2 && distance(currentPolygon[0], currentPolygon[currentPolygon.length - 1]) < 0.05) {
            groupPolygon.push([...currentPolygon, currentPolygon[0]]);
            currentPolygon = [];
        } else if (currentPolygon.length) {
            console.error("Could not find a segment to connect to", currentPolygon, segments, distance(currentPolygon[0], currentPolygon[currentPolygon.length - 1]));
            throw new Error("Could not find a segment to connect to");
        }

        // Sort so that the largest polygon is first
        const areasCache = new Map<unknown, number>();
        const getArea = (poly: [number, number][]) => {
            if (areasCache.has(poly)) {
                return areasCache.get(poly)!;
            }
            const ar = area(polygon([poly]));
            areasCache.set(poly, ar);
            return ar;
        }
        groupPolygon.sort((a, b) => {
            const areaA = getArea(a);
            const areaB = getArea(b);
            return areaB - areaA;
        });

        polygons.push(groupPolygon);
    }

    // let i = 0;
    for (const poly of polygons) {
        if (poly.length === 0) continue;
        features.push(polygon(poly, {
            color: `#${Math.floor(Math.random() * 0x1000000).toString(16).padStart(6, "0")}`,
        }));
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
                <div style={{flex: 1}}>
                    <MapGL
                        initialViewState={{
                            longitude: -122.4,
                            latitude: 37.8,
                            zoom: 14
                        }}
                        mapboxAccessToken="pk.eyJ1Ijoic29oY2FoIiwiYSI6ImNqeWVqcm8wdTAxc2MzaXFpa282Yzd2aHEifQ.afYbt2sVMZ-kbwdx5_PekQ"
                        mapStyle="mapbox://styles/mapbox/streets-v11"
                        fog={{}}
                        light={{}}
                    >
                        <Source type="geojson" data={{
                            type: "FeatureCollection",
                            features: generateFlattened(data)
                        }}>
                            <Layer type="line"
                                   paint={{
                                       "line-color": ["get", "color"],
                                       "line-width": 2
                                   }}/>
                            <Layer type="fill"
                                   paint={{
                                       "fill-color": ["get", "color"],
                                       "fill-opacity": 0.1
                                   }}/>
                        </Source>
                    </MapGL>
                </div>
            </div>
        </>
    )
}

export default App
