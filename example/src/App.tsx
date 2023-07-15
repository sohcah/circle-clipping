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
for (let i = 0; i < 10000; i++) {
    data.push([data[0][0] + r() * 0.05, data[0][1] + r() * 0.05, 0.03]);
}
console.log(r());

// const delaunay = delaunate.from(data.map((d) => [d[0], d[1]]));
// const triangleFeatures: Feature[] = [];
//
// // const triangles = delaunay.triangles;
// // for (let i = 0; i < triangles.length; i += 3) {
// //     triangleFeatures.push(lineString([
// //         [data[triangles[i]][0], data[triangles[i]][1]],
// //         [data[triangles[i + 1]][0], data[triangles[i + 1]][1]],
// //         [data[triangles[i + 2]][0], data[triangles[i + 2]][1]],
// //         [data[triangles[i]][0], data[triangles[i]][1]],
// //     ]));
// // }
// const halfedges = delaunay.halfedges;
// console.log(halfedges);
// for (let i = 0; i < halfedges.length; i++) {
//     const p = data[halfedges[i]];
//     if (!p) continue;
//     triangleFeatures.push(point(
//         [p[0], p[1]],
//     ));
// }

const delaunay = Delaunay.from(data.map((d) => [d[0], d[1]]));

// Compute the Laguerre–Voronoi diagram
const laguerreDiagram = delaunay.voronoi([
    Math.min(...data.map((d) => d[0])),
    Math.min(...data.map((d) => d[1])),
    Math.max(...data.map((d) => d[0])),
    Math.max(...data.map((d) => d[1])),
]);

// Calculate the weighted circumcenters (Laguerre–Voronoi vertices)
const circumcenters = [...laguerreDiagram.cellPolygons()].map((polygon, i) => {
    const weight = data[i][2] * 1000;
    const center = polygon.reduce((sum, point) => [sum[0] + point[0], sum[1] + point[1]], [0, 0]);
    return [center[0] / weight, center[1] / weight];
});

// type Point = [number, number]; // [latitude, longitude]

// function findCircleIntersections(circles: [longitude: number, latitude: number, radiusKm: number][]): [longitude: number, latitude: number][] {
//     const intersections: [longitude: number, latitude: number][] = [];
//
//     // Loop through each pair of circles
//     for (let i = 0; i < circles.length - 1; i++) {
//         const [lon1, lat1, r1] = circles[i];
//
//         for (let j = i + 1; j < circles.length; j++) {
//             const [lon2, lat2, r2] = circles[j];
//
//             // Calculate the distance between the centers of the circles
//             const dx = lon2 - lon1;
//             const dy = lat2 - lat1;
//             const distance = Math.sqrt(dx * dx + dy * dy);
//
//             // Check if the circles intersect
//             if (distance < r1 + r2) {
//                 // Calculate the intersection points
//                 const a = (r1 * r1 - r2 * r2 + distance * distance) / (2 * distance);
//                 const h = Math.sqrt(r1 * r1 - a * a);
//                 const cx = lon1 + (a * dx) / distance;
//                 const cy = lat1 + (a * dy) / distance;
//
//                 // Calculate the coordinates of the intersection points
//                 const x1 = cx + (h * dy) / distance;
//                 const y1 = cy - (h * dx) / distance;
//                 const x2 = cx - (h * dy) / distance;
//                 const y2 = cy + (h * dx) / distance;
//
//                 intersections.push([x1, y1]);
//                 intersections.push([x2, y2]);
//             }
//         }
//     }
//
//     return intersections;
// }

// console.log(findIntersectionPoints(data));


//
// const circumcenters: [number, number][] = [];
// for (let i = 0; i < delaunay.triangles.length; i += 3) {
//     const t0 = delaunay.triangles[i];
//     const t1 = delaunay.triangles[i + 1];
//     const t2 = delaunay.triangles[i + 2];
//
//     const x0 = data[t0][0];
//     const y0 = data[t0][1];
//     const x1 = data[t1][0];
//     const y1 = data[t1][1];
//     const x2 = data[t2][0];
//     const y2 = data[t2][1];
//
//     const D = 2 * (x0 * (y1 - y2) + x1 * (y2 - y0) + x2 * (y0 - y1));
//     const Ux = ((x0 * x0 + y0 * y0) * (y1 - y2) + (x1 * x1 + y1 * y1) * (y2 - y0) + (x2 * x2 + y2 * y2) * (y0 - y1)) / D;
//     const Uy = ((x0 * x0 + y0 * y0) * (x2 - x1) + (x1 * x1 + y1 * y1) * (x0 - x2) + (x2 * x2 + y2 * y2) * (x1 - x0)) / D;
//
//     circumcenters.push([Ux, Uy]);
// }
//
// const dualPoints: [number, number][] = [];
// for (let i = 0; i < data.length; i++) {
//     const region = delaunay.neighbors(i);
//     const vertices = [...region].map(index => circumcenters[index]);
//     const center = vertices.reduce((sum, vertex) => [sum[0] + vertex[0], sum[1] + vertex[1]], [0, 0]);
//     dualPoints.push([center[0] / vertices.length, center[1] / vertices.length]);
// }
//
// console.log(dualPoints);

const triangleFeatures: Feature[] = [];

// for (const dp of circumcenters) {
//     triangleFeatures.push(point(
//         [dp[0], dp[1]],
//     ));
// }

// const ints = findCircleIntersections(data);
// for (const intersection of ints) {
//     triangleFeatures.push(point(
//         [intersection[0], intersection[1]],
//     ));
// }
//
// console.log(ints);

type Circle = [number, number, number];


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

function group(circles: Circle[]): Circle[][] {

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

    const groups = new Map<number, Circle[]>();

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
function flatten(circles: Circle[]): Feature[] {
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
