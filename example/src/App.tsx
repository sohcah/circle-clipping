import './App.css'
import MapGL, {Layer, Source} from 'react-map-gl';
import {
    Feature,
    polygon,
    area,
} from "@turf/turf";
import seedrandom from 'seed-random';
import KDBush from "kdbush";
import * as geokdbush from "geokdbush-tk";
import CheapRuler from "cheap-ruler";

const rulers = new Array<CheapRuler>(1800);

function getRuler(lat: number): CheapRuler {
    const key = Math.floor(lat * 10);
    let ruler = rulers[key];
    if (!ruler) {
        ruler = new CheapRuler(lat, "kilometers");
        rulers[key] = ruler;
    }
    return ruler;
}

function distance(c1: [number, number, ...any[]], c2: [number, number, ...any[]]): number {
    return getRuler(c1[1]).distance(c1 as [number, number], c2 as [number, number]);
}

function destination(c1: [number, number, ...any[]], distance: number, bearing: number): [number, number] {
    return getRuler(c1[1]).destination(c1 as [number, number], distance, bearing);
}

const data: [number, number, number][] = [
    [-122.444768, 37.782551, 0.03],
    [-122.444586, 37.782745, 0.035],
    [-122.443688, 37.782842, 0.04],
    [-122.442715, 37.782932, 0.04],
];

const r = seedrandom("hi");

for (let i = 0; i < 5000; i++) {
    data.push([data[0][0] + r() * 0.01, data[0][1] + r() * 0.01, 0.03]);
}

interface IntersectionResult {
    angles: [number, number][];
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
    const dist = distance(circle1, circle2);

    // Check if the circles intersect
    if (dist > circle1[2] + circle2[2] || dist < Math.abs(circle1[2] - circle2[2])) {
        // No intersection points or circles are contained within each other
        return null;
    }

    const intersectionAngle = Math.acos(
        (circle1[2] * circle1[2] + dist * dist - circle2[2] * circle2[2]) /
        (2 * circle1[2] * dist)
    );

    const angle1 = Math.atan2(
        lon2 - lon1,
        Math.log(
            Math.tan(lat2 / 2 + Math.PI / 4) / Math.tan(lat1 / 2 + Math.PI / 4)
        )
    );

    const radiusAngle = circle1[2] / earthRadius;
    const cosRadiusAngle = Math.cos(radiusAngle);
    const sinRadiusAngle = Math.sin(radiusAngle);

    const angle2 = angle1 + intersectionAngle;
    const angle3 = angle1 - intersectionAngle;

    const sinLat1 = Math.sin(lat1);
    const cosLat1 = Math.cos(lat1);

    const sinAngle2 = Math.sin(angle2);
    const cosAngle2 = Math.cos(angle2);
    const p1latSinPart = sinLat1 * cosRadiusAngle +
        cosLat1 * sinRadiusAngle * cosAngle2;
    const p1latCosPart = Math.sqrt(1 - p1latSinPart * p1latSinPart);
    const p1lngSinPart = sinAngle2 * sinRadiusAngle * cosLat1;
    const p1lngCosPart = cosRadiusAngle - sinLat1 * p1latSinPart;
    const p1lngrad = lon1 + Math.atan2(p1lngSinPart, p1lngCosPart);

    const sinAngle3 = Math.sin(angle3);
    const cosAngle3 = Math.cos(angle3);
    const p2latSinPart = sinLat1 * cosRadiusAngle +
        cosLat1 * sinRadiusAngle * cosAngle3;
    const p2latCosPart = Math.sqrt(1 - p2latSinPart * p2latSinPart);

    const p2lngSinPart = sinAngle3 * sinRadiusAngle * cosLat1;
    const p2lngCosPart = cosRadiusAngle - sinLat1 * p2latSinPart;
    const p2lngrad = lon1 + Math.atan2(p2lngSinPart, p2lngCosPart);


    const bearingC1P1A = Math.sin(p1lngrad - lon1) * p1latCosPart;
    const bearingC1P1B = Math.cos(lat1) * p1latSinPart -
        Math.sin(lat1) * p1latCosPart * Math.cos(p1lngrad - lon1);
    const bearingC1P1 = Math.atan2(bearingC1P1A, bearingC1P1B);

    const bearingC1P2A = Math.sin(p2lngrad - lon1) * p2latCosPart;
    const bearingC1P2B = Math.cos(lat1) * p2latSinPart -
        Math.sin(lat1) * p2latCosPart * Math.cos(p2lngrad - lon1);
    const bearingC1P2 = Math.atan2(bearingC1P2A, bearingC1P2B);

    const bearingC2P1A = Math.sin(p1lngrad - lon2) * p1latCosPart;
    const bearingC2P1B = Math.cos(lat2) * p1latSinPart -
        Math.sin(lat2) * p1latCosPart * Math.cos(p1lngrad - lon2);
    const bearingC2P1 = Math.atan2(bearingC2P1A, bearingC2P1B);

    const bearingC2P2A = Math.sin(p2lngrad - lon2) * p2latCosPart;
    const bearingC2P2B = Math.cos(lat2) * p2latSinPart -
        Math.sin(lat2) * p2latCosPart * Math.cos(p2lngrad - lon2);
    const bearingC2P2 = Math.atan2(bearingC2P2A, bearingC2P2B);

    return {
        angles: [
            [toDegrees(bearingC1P1), toDegrees(bearingC2P1)],
            [toDegrees(bearingC1P2), toDegrees(bearingC2P2)],
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

function generateFlattened(circles: CircleArray[]): Feature[] {
    // const startClipping = performance.now();
    //
    // const circlePolys = circles.map(c => circle(c, c[2], {units: "kilometers"}));
    //
    //
    // union(circlePolys.map(i => i.geometry.coordinates as any)).map(i => ({
    //     type: "Feature" as const,
    //     properties: {},
    //     geometry: {
    //         type: "Polygon" as const,
    //         coordinates: i,
    //     },
    // }));
    //
    // const endClipping = performance.now();
    //
    // console.log(`Clipping took ${endClipping - startClipping}ms`);

    const features: Feature[] = [];

    // Benchmark the performance of calculating the intersection point
    const start = performance.now();
    let segmentStart = performance.now();

    const kdbush = new KDBush(circles.length);

    for (const circle of circles) {
        kdbush.add(circle[0], circle[1]);
    }

    kdbush.finish();

    const intersectingCircles: number[][] = new Array(circles.length).fill(null!).map(() => []);
    const intersectionAnglesByCircle: [number, number][][] = new Array(circles.length).fill(null!).map(() => [])

    const maximumRadius = circles.reduce((max, c) => c[2] > max ? c[2] : max, 0);

    const circleGroups: number[] = new Array(circles.length).fill(-1);

    const unprocessedCircles = new Set<number>();
    for (let i = 0; i < circles.length; i++) {
        unprocessedCircles.add(i);
    }

    const queue = new Set([0]);

    let currentGroup = 0;

    let attempts = 0;
    while (queue.size || unprocessedCircles.size) {
        if (!queue.size) {
            const nextCircle = unprocessedCircles.values().next().value;
            queue.add(nextCircle);
            currentGroup++;
        }
        if (attempts++ > circles.length) {
            throw new Error("Too many grouping attempts");
        }
        const i = queue.values().next().value;
        queue.delete(i);
        if (!unprocessedCircles.has(i)) {
            continue;
        }
        unprocessedCircles.delete(i);
        const circle = circles[i];
        circleGroups[i] = currentGroup;
        const potentialIntersections = geokdbush.around(kdbush, circle[0], circle[1], Infinity, maximumRadius + circle[2], c => c !== i && unprocessedCircles.has(c));
        for (const potentialIntersection of potentialIntersections) {
            const otherCircle = circles[potentialIntersection];

            if (circle[2] + otherCircle[2] > distance(circle, otherCircle)) {
                queue.add(potentialIntersection);
            }

            const intersection = calculateIntersection(
                circle,
                otherCircle,
            );

            if (intersection) {
                for (const [a1, a2] of intersection.angles) {
                    intersectionAnglesByCircle[i].push([a1, potentialIntersection]);
                    intersectionAnglesByCircle[potentialIntersection].push([a2, i]);
                    intersectingCircles[i].push(potentialIntersection);
                    intersectingCircles[potentialIntersection].push(i);
                }
            }
        }
    }

    console.log("number of intersection points", intersectionAnglesByCircle.reduce((sum, i) => sum + i.length, 0) / 2);

    console.log(`Grouping took ${performance.now() - segmentStart}ms`);
    segmentStart = performance.now();

    const circleArcs: [number, number, number, number][][] = new Array(circles.length).fill(null!).map(() => []);

    for (let i = 0; i < circles.length; i++) {
        const intersectionAngles = intersectionAnglesByCircle[i];

        if (intersectionAngles.length === 0) {
            circleArcs[i].push([0, 360, -1, -1]);
            continue;
        }

        const sortedAngles = intersectionAngles.slice().sort((a, b) => a[0] - b[0]);

        for (let j = 0; j < sortedAngles.length; j++) {
            let [startAngle, startIntersect] = sortedAngles[j];
            let [endAngle, endIntersect] = sortedAngles[(j + 1) % sortedAngles.length];

            if (endAngle < startAngle) {
                endAngle += 360;
            }

            circleArcs[i].push([startAngle, endAngle, startIntersect, endIntersect]);
        }
    }

    console.log(`Arcing took ${performance.now() - segmentStart}ms`);
    segmentStart = performance.now();

    const arcSegments: ([number, [number, number][]] | null)[] = [];

    // Create line segments for each circle arc with 10 points per arc
    for (let i = 0; i < circles.length; i++) {
        const circle = circles[i];
        const arcs = circleArcs[i];

        for (const arc of arcs) {
            let [startAngle, endAngle] = arc;
            if (startAngle > endAngle) {
                endAngle += 360;
            }
            const noPoints = Math.floor(Math.max(5, Math.abs(endAngle - startAngle) * 0.2));
            const deltaAngle = (endAngle - startAngle) / noPoints;

            const midpoint = destination(circle, circle[2], (startAngle + endAngle) / 2);
            if (intersectingCircles[i].some(c => distance(circles[c], midpoint) < circles[c][2])) {
                continue;
            }

            const points: [number, number][] = [];
            for (let p = 0; p <= noPoints; p++) {
                const angle = startAngle + p * deltaAngle;
                const pt = destination(circle, circle[2], angle);
                points.push(pt);
            }

            arcSegments.push([circleGroups[i], points]);
        }
    }

    console.log(`Found ${arcSegments.length} arc segments`);

    console.log(`Segmenting took ${performance.now() - segmentStart}ms`);
    segmentStart = performance.now();

    const arcSegmentsKDBush = new KDBush(arcSegments.length);

    for (let i = 0; i < arcSegments.length; i++) {
        const [, segment] = arcSegments[i]!;
        const middleCoordinate = segment[Math.floor(segment.length / 2)];
        arcSegmentsKDBush.add(middleCoordinate[0], middleCoordinate[1]);
    }

    arcSegmentsKDBush.finish();

    console.log(`Segment indexing took ${performance.now() - segmentStart}ms`);
    segmentStart = performance.now();

    const removedArcs = new Array(arcSegments.length).fill(false);
    for (const circle of circles) {
        const overlappingArcs = geokdbush.around(arcSegmentsKDBush, circle[0], circle[1], Infinity, circle[2] - 0.0001);
        for (const overlappingArc of overlappingArcs) {
            removedArcs[overlappingArc] = true;
        }
    }

    console.log(removedArcs.filter(a => a).length, "arcs removed");

    const filteredSegments = arcSegments.filter((_, n) => !removedArcs[n]);

    console.log(`Filtering took ${performance.now() - segmentStart}ms`);
    segmentStart = performance.now();


    const filteredSegmentsByGroup = filteredSegments
        .reduce((acc, segment) => {
            const [group, arc] = segment!;
            acc[group] ??= [];
            acc[group].push(arc!);
            return acc;
        }, {} as {
            [key: number]: [number, number][][]
        });

    console.log(`Grouping took ${performance.now() - segmentStart}ms`);
    segmentStart = performance.now();

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

    console.log(`Merging took ${performance.now() - segmentStart}ms`);
    segmentStart = performance.now();

    // let i = 0;
    for (const poly of polygons) {
        if (poly.length === 0) continue;
        features.push(polygon(poly, {
            color: `#${Math.floor(Math.random() * 0x1000000).toString(16).padStart(6, "0")}`,
        }));
    }

    console.log(`Featuring took ${performance.now() - segmentStart}ms`);
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
