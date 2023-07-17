import './App.css'
import MapGL, {Layer, Source} from 'react-map-gl';
import {
    Feature,
    polygon,
    area, bearing,
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

for (let i = 0; i < 500; i++) {
    data.push([data[0][0] + r() * 0.01, data[0][1] + r() * 0.01, 0.03]);
}

function bearingRadians(c1: [number, number, ...any[]], c2: [number, number, ...any[]]): number {
    const a = Math.sin(c2[0] - c1[0]) * Math.cos(c2[1]);
    const b = Math.cos(c1[1]) * Math.sin(c2[1]) -
        Math.sin(c1[1]) * Math.cos(c2[1]) * Math.cos(c2[0] - c1[0]);
    return Math.atan2(a, b);
}

interface IntersectionResult {
    points: [number, number][];
    angles: [number, number][];
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

// Calculate the coordinates of the intersection points
    const angle1 = Math.atan2(
        lon2 - lon1,
        Math.log(
            Math.tan(lat2 / 2 + Math.PI / 4) / Math.tan(lat1 / 2 + Math.PI / 4)
        )
    );

    const angle2 = angle1 + intersectionAngle;
    const angle3 = angle1 - intersectionAngle;

    const bearingToIntersection1FromCircle1o = angle2;
    const bearingToIntersection2FromCircle1o = angle3 < angle2 ? angle3 + Math.PI * 2 : angle3;

// Calculate intersection point 1 coordinates
    const sinLat1 = Math.sin(lat1);
    const cosLat1 = Math.cos(lat1);
    const cosAngle2 = Math.cos(angle2);
    const p1latSinPart = sinLat1 * Math.cos(circle1[2] / earthRadius) +
        cosLat1 * Math.sin(circle1[2] / earthRadius) * cosAngle2;
    const p1latrad = Math.asin(p1latSinPart);
    const p1lat = toDegrees(Math.asin(p1latSinPart));
    const p1lngSinPart = Math.sin(angle2) * Math.sin(circle1[2] / earthRadius) * cosLat1;
    const p1lngCosPart = Math.cos(circle1[2] / earthRadius) - sinLat1 * Math.sin(toRadians(p1lat));
    const p1lngrad = lon1 + Math.atan2(p1lngSinPart, p1lngCosPart);
    const p1lng = toDegrees(lon1 + Math.atan2(p1lngSinPart, p1lngCosPart));

// Calculate intersection point 2 coordinates
    const cosAngle3 = Math.cos(angle3);
    const p2latSinPart = sinLat1 * Math.cos(circle1[2] / earthRadius) +
        cosLat1 * Math.sin(circle1[2] / earthRadius) * cosAngle3;
    const p2latrad = Math.asin(p2latSinPart);
    const p2lat = toDegrees(Math.asin(p2latSinPart));
    const p2lngSinPart = Math.sin(angle3) * Math.sin(circle1[2] / earthRadius) * cosLat1;
    const p2lngCosPart = Math.cos(circle1[2] / earthRadius) - sinLat1 * Math.sin(toRadians(p2lat));
    const p2lngrad = lon1 + Math.atan2(p2lngSinPart, p2lngCosPart);
    const p2lng = toDegrees(lon1 + Math.atan2(p2lngSinPart, p2lngCosPart));

    // Using p1 and p2, find bearings from circle2
    // const bearingToIntersection1FromCircle2 = Math.atan2(
    //     Math.sin(lon1 + Math.atan2(p1lngSinPart, p1lngCosPart) - lon2) * Math.cos(Math.asin(p1latSinPart)),
    //     Math.cos(lat2) * Math.sin(Math.asin(p1latSinPart)) - Math.sin(lat2) * Math.cos(Math.asin(p1latSinPart)) * Math.cos(lon1 + Math.atan2(p1lngSinPart, p1lngCosPart) - lon2)
    // );
    //
    // const bearingToIntersection2FromCircle2 = Math.atan2(
    //     Math.sin(lon1 + Math.atan2(p2lngSinPart, p2lngCosPart) - lon2) * Math.cos(Math.asin(p2latSinPart)),
    //     Math.cos(lat2) * Math.sin(Math.asin(p2latSinPart)) - Math.sin(lat2) * Math.cos(Math.asin(p2latSinPart)) * Math.cos(lon1 + Math.atan2(p2lngSinPart, p2lngCosPart) - lon2)
    // );
    //
    // const bearingToIntersection1FromCircle1 = Math.atan2(
    //     Math.sin(lon2 + Math.atan2(p1lngSinPart, p1lngCosPart) - lon1) * Math.cos(Math.asin(p1latSinPart)),
    //     Math.cos(lat1) * Math.sin(Math.asin(p1latSinPart)) - Math.sin(lat1) * Math.cos(Math.asin(p1latSinPart)) * Math.cos(lon2 + Math.atan2(p1lngSinPart, p1lngCosPart) - lon1)
    // );
    //
    // const bearingToIntersection2FromCircle1 = Math.atan2(
    //     Math.sin(lon2 + Math.atan2(p2lngSinPart, p2lngCosPart) - lon1) * Math.cos(Math.asin(p2latSinPart)),
    //     Math.cos(lat1) * Math.sin(Math.asin(p2latSinPart)) - Math.sin(lat1) * Math.cos(Math.asin(p2latSinPart)) * Math.cos(lon2 + Math.atan2(p2lngSinPart, p2lngCosPart) - lon1)
    // );

    // console.log(toDegrees(bearingToIntersection1FromCircle2), toDegrees(bearingToIntersection2FromCircle2));

    return {
        points: [
            [p1lng, p1lat],
            [p2lng, p2lat]
        ],
        angles: [
            // [toDegrees(bearingToIntersection1FromCircle1), toDegrees(bearingToIntersection1FromCircle2)],
            // [toDegrees(bearingToIntersection2FromCircle1), toDegrees(bearingToIntersection2FromCircle2)],
            // [bearing(circle1, [p1lng, p1lat]), bearing(circle2, [p1lng, p1lat])],
            // [bearing(circle1, [p2lng, p2lat]), bearing(circle2, [p2lng, p2lat])],
            [toDegrees(bearingRadians([lon1, lat1], [p1lngrad, p1latrad])), toDegrees(bearingRadians([lon2, lat2], [p1lngrad, p1latrad]))],
            [toDegrees(bearingRadians([lon1, lat1], [p2lngrad, p2latrad])), toDegrees(bearingRadians([lon2, lat2], [p2lngrad, p2latrad]))],
        ],
        midpoint: [
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

    const intersectionAnglesByCircle: number[][] = new Array(circles.length).fill(null!).map(() => [])

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
                    intersectionAnglesByCircle[i].push(a1);
                    intersectionAnglesByCircle[potentialIntersection].push(a2);
                }
            }
        }
    }

    console.log(intersectionAnglesByCircle);

    console.log("number of intersection points", intersectionAnglesByCircle.reduce((sum, i) => sum + i.length, 0) / 2);

    console.log(`Grouping took ${performance.now() - segmentStart}ms`);
    segmentStart = performance.now();

    const circleArcs: [number, number][][] = new Array(circles.length).fill(null!).map(() => []);

    for (let i = 0; i < circles.length; i++) {
        const intersectionAngles = intersectionAnglesByCircle[i];

        if (intersectionAngles.length === 0) {
            circleArcs[i].push([0, 360]);
            continue;
        }

        const sortedAngles = intersectionAngles.slice().sort((a, b) => a - b);

        for (let j = 0; j < sortedAngles.length; j++) {
            let startAngle = sortedAngles[j]
            let endAngle = sortedAngles[(j + 1) % sortedAngles.length];

            if (endAngle < startAngle) {
                endAngle += 360;
            }

            circleArcs[i].push([startAngle, endAngle]);

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

            const points: [number, number][] = [];
            for (let p = 0; p <= noPoints; p++) {
                const angle = startAngle + p * deltaAngle;
                const pt = destination([circle[0], circle[1]], circle[2], angle);
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
