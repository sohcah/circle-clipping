import './App.css'
import MapGL, {Layer, Source} from 'react-map-gl/maplibre';
import {
    circle,
    distance,
    // union,
    destination,
    featureCollection,
    polygon,
    lineString,
    point,
    Feature,
    lineIntersect
} from '@turf/turf';
import seedrandom from 'seed-random';
import KDBush from "kdbush";
import * as geokdbush from "geokdbush-tk";
import Polygon from "polygon-clipping";
import {Delaunay} from "d3-delaunay";

const {union} = Polygon;

import delaunate from "delaunator";

const data: [number, number, number][] = [
    [-122.444768, 37.782551, 0.03],
    [-122.444586, 37.782745, 0.035],
    [-122.443688, 37.782842, 0.04],
    [-122.442715, 37.782932, 0.04],
];

const r = seedrandom("hi");
//
for (let i = 0; i < 100; i++) {
    data.push([data[0][0] + r() * 0.05, data[0][1] + r() * 0.05, 0.03]);
}
console.log(r());

interface Point {
    lat: number;
    lng: number;
}

interface Circle {
    center: Point;
    radius: number;
}

function calculateIntersectionPoints(circle1: Circle, circle2: Circle): Point[] {
    // Convert latitude and longitude to radians
    const lat1 = toRadians(circle1.center.lat);
    const lng1 = toRadians(circle1.center.lng);
    const lat2 = toRadians(circle2.center.lat);
    const lng2 = toRadians(circle2.center.lng);

    // Calculate angular distance between the centers of the circles
    const centralAngle = getCentralAngle(lat1, lng1, lat2, lng2);

    // Check if the circles are disjoint or identical
    if (centralAngle > circle1.radius + circle2.radius || centralAngle < Math.abs(circle1.radius - circle2.radius)) {
        // No intersection points
        return [];
    }

    // Calculate the azimuth angle from circle1 to circle2
    const azimuth12 = Math.atan2(
        Math.sin(lng2 - lng1) * Math.cos(lat2),
        Math.cos(lat1) * Math.sin(lat2) - Math.sin(lat1) * Math.cos(lat2) * Math.cos(lng2 - lng1)
    );

    // Calculate the azimuth angle from circle2 to circle1
    const azimuth21 = Math.atan2(
        Math.sin(lng1 - lng2) * Math.cos(lat1),
        Math.cos(lat2) * Math.sin(lat1) - Math.sin(lat2) * Math.cos(lat1) * Math.cos(lng1 - lng2)
    );

    // Calculate the angular distance from the center of circle1 to the intersection points
    const angularDistanceFromCenter1 = Math.acos(
        (Math.cos(circle2.radius) - Math.cos(circle1.radius) * Math.cos(centralAngle - azimuth21)) /
        (Math.sin(circle1.radius) * Math.sin(centralAngle - azimuth21))
    );

    // Calculate the angular distance from the center of circle2 to the intersection points
    const angularDistanceFromCenter2 = Math.acos(
        (Math.cos(circle1.radius) - Math.cos(circle2.radius) * Math.cos(centralAngle - azimuth12)) /
        (Math.sin(circle2.radius) * Math.sin(centralAngle - azimuth12))
    );

    // Calculate the coordinates of the intersection points
    const intersectionPoint1: Point = {
        lat: toDegrees(lat1 + Math.atan2(Math.cos(circle1.radius) * Math.sin(angularDistanceFromCenter1), Math.sin(circle1.radius) * Math.cos(angularDistanceFromCenter1) * Math.cos(azimuth12))),
        lng: toDegrees(lng1 + Math.atan2(Math.sin(azimuth12) * Math.sin(angularDistanceFromCenter1), Math.cos(angularDistanceFromCenter1) - Math.sin(lat1) * Math.sin(lat1)))
    };

    const intersectionPoint2: Point = {
        lat: toDegrees(lat2 + Math.atan2(Math.cos(circle2.radius) * Math.sin(angularDistanceFromCenter2), Math.sin(circle2.radius) * Math.cos(angularDistanceFromCenter2) * Math.cos(azimuth21))),
        lng: toDegrees(lng2 + Math.atan2(Math.sin(azimuth21) * Math.sin(angularDistanceFromCenter2), Math.cos(angularDistanceFromCenter2) - Math.sin(lat2) * Math.sin(lat2)))
    };

    return [intersectionPoint1, intersectionPoint2];
}

function toRadians(degrees: number): number {
    return (degrees * Math.PI) / 180;
}

function toDegrees(radians: number): number {
    return (radians * 180) / Math.PI;
}

function getCentralAngle(lat1: number, lng1: number, lat2: number, lng2: number): number {
    const deltaLng = Math.abs(lng2 - lng1);
    const numerator = Math.sqrt(
        Math.pow(Math.cos(lat2) * Math.sin(deltaLng), 2) +
        Math.pow(Math.cos(lat1) * Math.sin(lat2) - Math.sin(lat1) * Math.cos(lat2) * Math.cos(deltaLng), 2)
    );
    const denominator = Math.sin(lat1) * Math.sin(lat2) + Math.cos(lat1) * Math.cos(lat2) * Math.cos(deltaLng);
    return Math.atan2(numerator, denominator);
}

const triangleFeatures: Feature[] = [];

// const intPoints = findIntersectionPoints(data);
// for (const [lon, lat] of intPoints) {
//     triangleFeatures.push(point([lon, lat]));
// }

const intersectionOfFirstTwo = calculateIntersectionPoints({
    center: {lat: data[0][1], lng: data[0][0]},
    radius: data[0][2] * 1000 / 6378137,
}, {
    center: {lat: data[1][1], lng: data[1][0]},
    radius: data[1][2] * 1000 / 6378137,
});
console.log(intersectionOfFirstTwo, data[0], data[1]);
triangleFeatures.push(point([intersectionOfFirstTwo[0].lng, intersectionOfFirstTwo[0].lat]));
triangleFeatures.push(point([intersectionOfFirstTwo[1].lng, intersectionOfFirstTwo[1].lat]));


type CircleArray = [number, number, number];


const pointsPerCircle = 32;

function getCircleIndex(pointIndex: number) {
    return Math.floor(pointIndex / pointsPerCircle);
}

function getLeftNeighbor(pointIndex: number) {
    return ((pointIndex + pointsPerCircle - 1) % pointsPerCircle) + getCircleIndex(pointIndex) * pointsPerCircle
}

function getRightNeighbor(pointIndex: number) {
    return ((pointIndex + 1) % pointsPerCircle) + getCircleIndex(pointIndex) * pointsPerCircle;
}

function group(circles: CircleArray[]): CircleArray[][] {

    const maximumRadius = circles.reduce((max, c) => c[2] > max ? c[2] : max, 0);

    const kdbush = new KDBush(circles.length);

    for (const circle of circles) {
        kdbush.add(circle[0], circle[1]);
    }

    kdbush.finish();

    const circleGroups = new Array<number>(circles.length);

    for (let i = 0; i < circles.length; i++) {
        const potentialIntersections = geokdbush.around(
            kdbush,
            circles[i][0],
            circles[i][1],
            Infinity,
            circles[i][2] + maximumRadius,
            (c) => c !== i
        );

        const intersections = potentialIntersections.filter((j) => {
            const circle = circles[j];
            const dist = distance([circles[i][0], circles[i][1]], [circle[0], circle[1]]);
            return dist < circles[i][2] + circle[2];
        });

        for (const intersection of intersections) {
            if (circleGroups[intersection] === undefined) {
                circleGroups[intersection] = i;
            } else {
                const groupToUpdate = circleGroups[intersection];
                for (let j = 0; j < circles.length; j++) {
                    if (circleGroups[j] === groupToUpdate) {
                        circleGroups[j] = i;
                    }
                }
            }
        }

        if (circleGroups[i] === undefined) {
            circleGroups[i] = i;
        }
    }

    const groups = new Map<number, CircleArray[]>();

    for (let i = 0; i < circles.length; i++) {
        const group = groups.get(circleGroups[i]);
        if (group === undefined) {
            groups.set(circleGroups[i], [circles[i]]);
        } else {
            group.push(circles[i]);
        }
    }

    return Array.from(groups.values());
}

// Union a set of circles into a set of polygons
function flatten(circles: CircleArray[]): Feature[] {
    let log: ((...v: any[]) => void) | null = null;
    if (Math.PI === 3) {
        log = console.log;
    }
    const start = performance.now();
    /*
    First, identify circles that belong to the same group of overlapping circles. Obviously, two circles overlap if the absolute distance of their centres is less than their combined radii. For each circle, store the circles it itersecs with in a hashmap/dictionary. (You can extend this to entire groups of circles using the union-find algorithm, or disjoint sets, but this is not really needed.)
When creating the points belonging to one circle, memorize the left and right neighbor of each point. You get this for free by just storing them in an ordered array.
Remove the "inner" points, i.e. those that are closer than one radius to any of the centres of the circles intersecting the first circle.
Mark the neighbors to those inner points that were not removed the "loose ends" of the circles.
Connect the points that were not removed, including the loose ends, to their original left and right neighbors.
For each loose end point, find the closest loose end of a different circle in the same group and connect them.
     */
    const colors = new Array<string>(circles.length).fill(null!).map(() => `#${Math.floor(Math.random() * 0x1000000).toString(16).padStart(6, "0")}`);

    const maximumRadius = circles.reduce((max, c) => c[2] > max ? c[2] : max, 0);

    const kdbush = new KDBush(circles.length);

    for (const circle of circles) {
        kdbush.add(circle[0], circle[1]);
    }

    kdbush.finish();

    // const circleGroups = new Array<number>(circles.length);
    // for (let i = 0; i < circles.length; i++) {
    //     const potentialIntersections = geokdbush.around(kdbush, circles[i][0], circles[i][1], Infinity, circles[i][2] + maximumRadius, c => c !== i);
    //     const intersections = potentialIntersections.filter((j) => {
    //         const circle = circles[j];
    //         const dist = distance([circles[i][0], circles[i][1]], [circle[0], circle[1]]);
    //         return dist < circles[i][2] + circle[2];
    //     });
    //     console.log(i, intersections);
    //     for (const intersection of intersections) {
    //         if (circleGroups[intersection] === undefined) continue;
    //         if (circleGroups[i] === undefined) {
    //             circleGroups[i] = intersection;
    //         } else {
    //             const groupToUpdate = circleGroups[intersection];
    //             for (let j = 0; j < circles.length; j++) {
    //                 if (circleGroups[j] === groupToUpdate) {
    //                     circleGroups[j] = circleGroups[i];
    //                 }
    //             }
    //         }
    //     }
    //     if (circleGroups[i] === undefined) circleGroups[i] = i;
    // }
    const circleGroups = new Array<number>(circles.length);

    for (let i = 0; i < circles.length; i++) {
        const potentialIntersections = geokdbush.around(
            kdbush,
            circles[i][0],
            circles[i][1],
            Infinity,
            circles[i][2] + maximumRadius,
            (c) => c !== i
        );

        const intersections = potentialIntersections.filter((j) => {
            const circle = circles[j];
            const dist = distance([circles[i][0], circles[i][1]], [circle[0], circle[1]]);
            return dist < circles[i][2] + circle[2];
        });

        for (const intersection of intersections) {
            if (circleGroups[intersection] === undefined) {
                circleGroups[intersection] = i;
            } else {
                const groupToUpdate = circleGroups[intersection];
                for (let j = 0; j < circles.length; j++) {
                    if (circleGroups[j] === groupToUpdate) {
                        circleGroups[j] = i;
                    }
                }
            }
        }

        if (circleGroups[i] === undefined) {
            circleGroups[i] = i;
        }
    }

    // const circlesByCenter = new Array<number[]>(circles.length);
    // for (let i = 0; i < circles.length; i++) {
    //     const potentialIntersections = geokdbush.around(kdbush, circles[i][0], circles[i][1], Infinity, circles[i][2] + maximumRadius, c => c !== i);
    //     circlesByCenter[i] = potentialIntersections.filter((j) => {
    //         const circle = circles[j];
    //         const dist = distance([circles[i][0], circles[i][1]], [circle[0], circle[1]]);
    //         return dist < circles[i][2] + circle[2];
    //     });
    // }
    // console.log(circlesByCenter);
    // // Make a circleIndex -> sectionId array
    // const circleSectionIds = new Array<number>(circles.length).fill(null!);
    // let id = 0;
    //
    // let count = 0;
    // let changes = true;
    // while (circleSectionIds.includes(null!) || changes) {
    //     changes = false;
    //     if (count++ > circles.length) throw new Error("Infinite loop");
    //     let hasSetIdThisIteration = false;
    //     for (let i = 0; i < circles.length; i++) {
    //         if (circleSectionIds[i] !== null) continue;
    //         if (!hasSetIdThisIteration) {
    //             changes = true;
    //             circleSectionIds[i] = id++;
    //             hasSetIdThisIteration = true;
    //         }
    //         const neighbors = circlesByCenter[i];
    //         for (const neighbor of neighbors) {
    //             if (circleSectionIds[neighbor] !== null) {
    //                 if (circleSectionIds[i] !== null && circleSectionIds[neighbor] !== circleSectionIds[i]) {
    //                     console.warn("Circle is in two sections");
    //                     for (let j = 0; j < circles.length; j++) {
    //                         if (circleSectionIds[j] === circleSectionIds[i]) {
    //                             circleSectionIds[j] = circleSectionIds[neighbor];
    //                             changes = true;
    //                         }
    //                     }
    //                 }
    //             } else {
    //                 circleSectionIds[neighbor] = circleSectionIds[i];
    //                 changes = true;
    //             }
    //         }
    //     }
    // }

    const circleSectionIds = circleGroups;

    // console.log("circlesByCenter", circlesByCenter);

    const ringPoints = new Float64Array(circles.length * pointsPerCircle * 2);

    for (let i = 0; i < circles.length; i++) {
        for (let j = 0; j < pointsPerCircle; j++) {
            const dest = destination([circles[i][0], circles[i][1]], circles[i][2], j * 360 / pointsPerCircle);
            ringPoints[i * pointsPerCircle * 2 + j * 2] = dest.geometry.coordinates[0];
            ringPoints[i * pointsPerCircle * 2 + j * 2 + 1] = dest.geometry.coordinates[1];
        }
    }

    function getIntersectionLineStrings(a: number, b: number) {
        const a1 = getLeftNeighbor(a);
        const a2 = a;
        const a3 = getRightNeighbor(a);
        const a4 = getRightNeighbor(a3);
        const b1 = getLeftNeighbor(b);
        const b2 = b;
        const b3 = getRightNeighbor(b);
        const b4 = getRightNeighbor(b3);
        return [
            lineString([
                [ringPoints[a1 * 2], ringPoints[a1 * 2 + 1]],
                [ringPoints[a2 * 2], ringPoints[a2 * 2 + 1]],
                [ringPoints[a3 * 2], ringPoints[a3 * 2 + 1]],
                [ringPoints[a4 * 2], ringPoints[a4 * 2 + 1]],
            ]),
            lineString([
                [ringPoints[b1 * 2], ringPoints[b1 * 2 + 1]],
                [ringPoints[b2 * 2], ringPoints[b2 * 2 + 1]],
                [ringPoints[b3 * 2], ringPoints[b3 * 2 + 1]],
                [ringPoints[b4 * 2], ringPoints[b4 * 2 + 1]],
            ])
        ]
    }

    function getIntersect(a: number, b: number) {

        return lineIntersect(
            ...getIntersectionLineStrings(a, b)
        ).features[0]?.geometry.coordinates
    }

    const ringPointsKdbush = new KDBush(ringPoints.length / 2);

    for (let i = 0; i < ringPoints.length; i += 2) {
        ringPointsKdbush.add(ringPoints[i], ringPoints[i + 1]);
    }

    ringPointsKdbush.finish();

    const overlaps = new Array<number>(ringPoints.length / 2).fill(0);

    let overlapId = 0;
    for (const circle of circles) {
        overlapId++;
        const overlappingPoints = geokdbush.around(ringPointsKdbush, circle[0], circle[1], Infinity, circle[2] - 0.0001);
        for (const point of overlappingPoints) {
            overlaps[point] = overlapId;
        }
    }

    const polygons: number[][] = [];
    let currentPolygon: number[] = [];
    const hasBeenVisited = overlaps.map(i => !!i);
    const extraPoints: [number, number][] = [[0, 0]];

    let nextPointIndex: number | null = 0;

    const extraFeatures: Feature[] = [];

    while (nextPointIndex !== null) {
        if (nextPointIndex >= ringPoints.length / 2) {
            const potentialNextPoint = hasBeenVisited.indexOf(false);
            if (potentialNextPoint === -1) {
                nextPointIndex = null;
            } else {
                nextPointIndex = potentialNextPoint;
            }
            continue;
        }
        const pointIndex: number = nextPointIndex;
        nextPointIndex = null;
        if (hasBeenVisited[pointIndex]) {
            // Complete the circle
            if (currentPolygon.length) {
                log?.("complete circle visited", pointIndex);
                polygons.push(currentPolygon);
                currentPolygon = [];
            }
            nextPointIndex = pointIndex + 1;
            continue;
        }
        hasBeenVisited[pointIndex] = true;
        const circleIndex = getCircleIndex(pointIndex);
        const rightIndex = getRightNeighbor(pointIndex);
        const isLastInCircle = pointIndex % pointsPerCircle === pointsPerCircle - 1;

        currentPolygon.push(pointIndex);

        if (overlaps[rightIndex]) {
            // We need to merge into the next circle
            const mergeToIndex = geokdbush.around(ringPointsKdbush, ringPoints[pointIndex * 2], ringPoints[pointIndex * 2 + 1], 1, circles[circleIndex][2] / 2, nextPointIndex => {
                const nextCircleIndex = getCircleIndex(nextPointIndex);
                const nextRightIndex = getRightNeighbor(nextPointIndex);
                return !!overlaps[nextPointIndex] && !overlaps[nextRightIndex] && nextCircleIndex !== circleIndex && circleSectionIds[circleIndex] === circleSectionIds[nextCircleIndex];
            })[0];
            //     ,
            //     // ...geokdbush.around(ringPointsKdbush, ringPoints[pointIndex * 2], ringPoints[pointIndex * 2 + 1], 1, circles[circleIndex][2] / 2, nextPointIndex => {
            //     //     const nextCircleIndex = getCircleIndex(nextPointIndex);
            //     //     const nextRightIndex = getRightNeighbor(nextPointIndex);
            //     //     const nextLeftIndex = getLeftNeighbor(nextPointIndex);
            //     //     return !!overlaps[nextPointIndex] && !!overlaps[nextRightIndex] && (overlaps[nextPointIndex] !== overlaps[nextRightIndex]) && nextCircleIndex !== circleIndex && circleSectionIds[circleIndex] === circleSectionIds[nextCircleIndex];
            //     // })
            // ];

            log?.("mergeToIndex", mergeToIndex);
            if (mergeToIndex !== undefined) {
                const intersection = getIntersect(mergeToIndex, pointIndex);
                log?.("intersection", intersection);
                if (intersection) {
                    extraPoints.push([intersection[0], intersection[1]]);

                    currentPolygon.push(-(extraPoints.length - 1));
                    nextPointIndex = getRightNeighbor(mergeToIndex);
                } else {
                    currentPolygon.push(getRightNeighbor(mergeToIndex));
                    nextPointIndex = getRightNeighbor(getRightNeighbor(mergeToIndex))
                    extraFeatures.push(
                        ...getIntersectionLineStrings(mergeToIndex, pointIndex)
                    );
                }
            } else {
                // Complete the circle
                log?.("complete circle no merge", pointIndex);
                polygons.push(currentPolygon);
                currentPolygon = [];
                nextPointIndex = pointIndex + 1;
            }

        } else if (isLastInCircle) {
            nextPointIndex = rightIndex;
        } else {
            // Continue the circle
            nextPointIndex = rightIndex;
        }
    }
    if (currentPolygon.length) polygons.push(currentPolygon);

    log?.(polygons);


    const result = [...polygons.flatMap(polygon => {
        if (!polygon.length) return [];
        if (polygon.length < 2) return [point([ringPoints[polygon[0] * 2], ringPoints[polygon[0] * 2 + 1]])];
        return [
            lineString([...polygon, polygon[0]].map(i => i < 0 ? extraPoints[-i] : [ringPoints[i * 2], ringPoints[i * 2 + 1]]), {
                color: colors[circleSectionIds[getCircleIndex(polygon[0])]],
                opacity: 0.5,
            }),
        ]
    }), ...extraFeatures];

    console.log(performance.now() - start);

    return result;


    //
    const ringPointsFeatures = [];

    for (let circleIndex = 0; circleIndex < circles.length; circleIndex++) {
        const points = [];
        for (let i = circleIndex * pointsPerCircle * 2; i < (circleIndex * pointsPerCircle + pointsPerCircle) * 2; i += 2) {
            const pointIndex = i / 2;
            const leftIndex = getLeftNeighbor(pointIndex)
            const rightIndex = getRightNeighbor(pointIndex)
            if (!(overlaps[pointIndex] && overlaps[leftIndex] && overlaps[rightIndex])) {

                points.push([ringPoints[i], ringPoints[i + 1]]);
            }
            if (overlaps[pointIndex]) {
                if (!overlaps[leftIndex] || !overlaps[rightIndex]) {
                    ringPointsFeatures.push(
                        point([ringPoints[i], ringPoints[i + 1]], {
                            color: "#ff0000",
                        })
                    );
                } else {
                    ringPointsFeatures.push(
                        point([ringPoints[i], ringPoints[i + 1]], {
                            color: "#000000",
                        })
                    );
                }
            } else {
                if (overlaps[leftIndex] || overlaps[rightIndex]) {
                    ringPointsFeatures.push(
                        point([ringPoints[i], ringPoints[i + 1]], {
                            color: "#0000ff",
                        })
                    );
                } else {
                    ringPointsFeatures.push(
                        point([ringPoints[i], ringPoints[i + 1]], {
                            color: "#00ff00",
                        })
                    );
                }
            }

        }
        points.push(points[0])
        if (points.length >= 2) {
            ringPointsFeatures.unshift(lineString(points, {
                opacity: 0,
            }))
        }
    }

    return ringPointsFeatures;

    // console.log(ringPoints);
    //
    //
    // const circlesByCenter: number[][] = [];
    // for (let i = 0; i < circles.length; i++) {
    //     const potentialIntersections = geokdbush.around(kdbush, circles[i][0], circles[i][1], maximumRadius);
    //     circlesByCenter.push(potentialIntersections.filter((j) => {
    //         const circle = circles[j];
    //         return geokdbush.distance(circles[i][0], circles[i][1], circle[0], circle[1]) < circles[i][2] + circle[2];
    //     }));
    // }


    // const remainingCirclesSet = new Set(circles);
    // if (circles.length === 0) return [];
    //
    // const polygons: [number, number][][] = [];
    //
    // let currentPolygon: [number, number][] | null = null;
    //
    // const queue = [circles[0]];
    // remainingCirclesSet.delete(circles[0]);
    //
    // let i = 0;
    // while (queue.length > 0) {
    //     if(i++ > 100_000) {
    //         break;
    //     }
    //     const circleEntry = queue.pop()!;
    //     if (currentPolygon) {
    //         currentPolygon = union(currentPolygon, circle(circleEntry.center, circleEntry.radius));
    //     } else {
    //         currentPolygon = circle(circleEntry.center, circleEntry.radius);
    //     }
    //     const intersections = [...remainingCirclesSet].filter((c) => {
    //         return distance(circleEntry.center, c.center) < circleEntry.radius + c.radius;
    //     });
    //     intersections.forEach((c) => remainingCirclesSet.delete(c));
    //     queue.push(...intersections);
    //     if (!queue.length) {
    //         polygons.push(currentPolygon!);
    //         currentPolygon = null;
    //         const next = remainingCirclesSet.values().next().value;
    //         if(next) {
    //             queue.push(next);
    //         }
    //     }
    // }
    //
    // return polygons;
}

let flattenTime = 0;
let unionTime = 0;
// let groupUnionTime = 0;
// for(let i = 0; i < 1;i++) {
//     const start = performance.now();
//     flatten(data);
//     flattenTime += performance.now() - start;
//
//     const start2 = performance.now();
//     const circles = data.map((d) => circle([d[0], d[1]], d[2], {
//         steps: pointsPerCircle,
//     }))
//
//     union(circles.map(i => i.geometry.coordinates as any)).map(i => ({
//         type: "Feature" as const,
//         geometry: {
//             type: "Polygon" as const,
//             properties: {},
//             coordinates: i,
//         },
//     }))
//     unionTime += performance.now() - start2;
//
//     // const start3 = performance.now();
//     // const groupedCircles = group(data);
//     // for(const group of groupedCircles) {
//     //     const circles2 = group.map((d) => circle([d[0], d[1]], d[2], {
//     //         steps: pointsPerCircle,
//     //     }))
//     //
//     //     union(circles2.map(i => i.geometry.coordinates as any)).map(i => ({
//     //         type: "Feature" as const,
//     //         geometry: {
//     //             type: "Polygon" as const,
//     //             properties: {},
//     //             coordinates: i,
//     //         },
//     //     }))
//     // }
//     // console.log(groupedCircles);
//     // groupUnionTime += performance.now() - start3;
//
// }

console.log(flattenTime, unionTime);//, groupUnionTime);

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
                            features: triangleFeatures
                            // .map(d => d.length < 2 ? point(d[0]) : d.length < 4 ? lineString(d) : polygon(d))
                        }}>
                            <Layer type="line"
                                   paint={{
                                       "line-color": "#000000",
                                       "line-width": 2
                                   }}/>
                            <Layer type="circle"
                                   paint={{
                                       "circle-color": "#000000",
                                       "circle-radius": 2
                                   }}/>
                        </Source>
                    </MapGL>
                </div>
            </div>
        </>
    )
}

export default App
