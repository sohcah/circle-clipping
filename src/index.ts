import {
    type Feature,
    polygon,
} from "@turf/helpers";
import area from "@turf/area";
import turfDistance from "@turf/distance";
import turfDestination from "@turf/destination";
import KDBush from "kdbush";
import { around } from "geokdbush-tk";
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

function fastDistance(c1: [number, number, ...any[]], c2: [number, number, ...any[]]): number {
  return getRuler(c1[1]).distance(c1 as [number, number], c2 as [number, number]);
}
function slowDistance(c1: [number, number, ...any[]], c2: [number, number, ...any[]]): number {
    return turfDistance(c1, c2, {units: "kilometers"});
}

function fastDestination(c1: [number, number, ...any[]], distance: number, bearing: number): [number, number] {
  return getRuler(c1[1]).destination(c1 as [number, number], distance, bearing);
}
function slowDestination(c1: [number, number, ...any[]], distance: number, bearing: number): [number, number] {
  const res = turfDestination(c1, distance, bearing, {units: "kilometers"});
  return res.geometry.coordinates as [number, number];
}


type CircleArray = [number, number, number];
const earthRadius = 6371.0088; // Radius of the Earth in kilometers

function calculateIntersection(circle1: CircleArray, circle2: CircleArray, useFastMath: boolean): [number, number][] | null {
    // Calculate the distance between the circle centers
    const dist = (useFastMath ? fastDistance : slowDistance)(circle1, circle2);

    // Check if the circles intersect
    if (dist > circle1[2] + circle2[2] || dist < Math.abs(circle1[2] - circle2[2])) {
        // No intersection points or circles are contained within each other
        return null;
    }

    const intersectionAngle = Math.acos(
        (circle1[2] * circle1[2] + dist * dist - circle2[2] * circle2[2]) /
        (2 * circle1[2] * dist)
    );

  // Convert the coordinates to radians
  const lat1 = toRadians(circle1[1]);
  const lon1 = toRadians(circle1[0]);
  const lat2 = toRadians(circle2[1]);
  const lon2 = toRadians(circle2[0]);


  const angle1 = Math.atan2(
    Math.sin(lon2 - lon1) * Math.cos(lat2),
    Math.cos(lat1) * Math.sin(lat2) - Math.sin(lat1) * Math.cos(lat2) * Math.cos(lon2 - lon1)
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

  return [
        [toDegrees(bearingC1P1), toDegrees(bearingC2P1)],
        [toDegrees(bearingC1P2), toDegrees(bearingC2P2)],
    ];
}

// Helper function to convert degrees to radians
function toRadians(degrees: number): number {
    return degrees * (Math.PI / 180);
}

// Helper function to convert radians to degrees
function toDegrees(radians: number): number {
    return radians * (180 / Math.PI);
}

export function circleUnion<T extends {
    [key: string]: any
}>(circles: CircleArray[], properties: T, { logPerf = false, useFastMath = true }: {
    logPerf?: boolean
    useFastMath?: boolean
} = {}): Feature[] {
    const distance = useFastMath ? fastDistance : slowDistance;
    const destination = useFastMath ? fastDestination : slowDestination;

    // This is a workaround to deal with the function sometimes
    // fails to merge segments when it receives longitudes > 180
    // so this function will convert all longitudes to be between
    // -180 and 180
    // At some point we should fix the underlying issue, however
    // there isn't any obvious reason why this is happening
    circles = circles.reduce((array, circle) => {
      if (Math.abs(circle[1]) > 85) {
        console.error(`Received invalid latitude, ${circle[1]}, skipping`);
        return array;
      }
      if (circle[2] <= 0) {
        console.error(`Received invalid radius, ${circle[2]}, skipping`);
        return array;
      }
      const boundedLng = circle[0] % 360;
      array.push([
        boundedLng >= 180 ? boundedLng - 360 : (boundedLng < -180 ? boundedLng + 360 : boundedLng),
        circle[1],
        circle[2],
      ]);
      return array;
    }, [] as CircleArray[]);

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
    // if (logPerf) console.log(`Clipping took ${endClipping - startClipping}ms`);

    // Deduplicate circles that are exactly the same
    const circlesSet = new Set<string>();

    circles = circles.filter(c => {
        const key = c.join(",");
        if (circlesSet.has(key)) {
            return false;
        } else {
            circlesSet.add(key);
            return true;
        }
    }).sort((a, b) => b[2] - a[2]);

    const features: Feature[] = [];

    // Benchmark the performance of calculating the intersection point
    const start = performance.now();
    let segmentStart = performance.now();

    const intersectingCircles: number[][] = new Array(circles.length).fill(null!).map(() => []);
    const intersectionAnglesByCircle: number[][] = new Array(circles.length).fill(null!).map(() => [])

    const maximumRadius = circles.reduce((max, c) => c[2] > max ? c[2] : max, 0);

    const circleGroups: number[] = new Array(circles.length).fill(-1);

    const unprocessedCircles = new Set<number>();
    for (let i = 0; i < circles.length; i++) {
        unprocessedCircles.add(i);
    }

    let kdbush: KDBush = null!;

    let kdbushIndex: number[] = null!;
    function generateKDBush() {
      kdbush = new KDBush(unprocessedCircles.size)
      kdbushIndex = new Array(unprocessedCircles.size).fill(0);

      let n = 0;
      for (const i of unprocessedCircles) {
        const circle = circles[i];
        kdbush.add(circle[0], circle[1]);
        kdbushIndex[n] = i;
        n++;
      }

      kdbush.finish();
    }

    generateKDBush();

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
        const potentialIntersections = around(kdbush, circle[0], circle[1], Infinity, maximumRadius + circle[2], c => kdbushIndex[c] !== i && unprocessedCircles.has(kdbushIndex[c]));
        let overlaps = 0;
        for (const potentialIntersectionKdbush of potentialIntersections) {
            const potentialIntersection = kdbushIndex[potentialIntersectionKdbush];
            const otherCircle = circles[potentialIntersection];

            if (circle[2] + otherCircle[2] > distance(circle, otherCircle)) {

                const intersectionAngles = calculateIntersection(
                    circle,
                    otherCircle,
                    useFastMath
                );

                if (intersectionAngles) {
                    queue.add(potentialIntersection);
                    for (const [a1, a2] of intersectionAngles) {
                        intersectionAnglesByCircle[i].push(a1);
                        intersectionAnglesByCircle[potentialIntersection].push(a2);
                    }
                } else if (otherCircle[2] < circle[2]) {
                  // Other Circle is entire contained within Circle, so we can skip it
                  unprocessedCircles.delete(potentialIntersection);
                  overlaps++;
                }
                intersectingCircles[i].push(potentialIntersection);
                intersectingCircles[potentialIntersection].push(i);
            }
        }
        if (overlaps > 50) {
          // If we've skipped lots of circles, we should regenerate the KDBush index.
          if (logPerf) console.log(`Regenerating KDBush with ${unprocessedCircles.size} circles (${overlaps} overlaps)`);
          generateKDBush();
        }
    }

    if (logPerf) console.log(`Grouping took ${performance.now() - segmentStart}ms`);
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
            let startAngle = sortedAngles[j];
            let endAngle = sortedAngles[(j + 1) % sortedAngles.length];

            if (endAngle < startAngle) {
                endAngle += 360;
            }

            circleArcs[i].push([startAngle, endAngle]);
        }
    }

    if (logPerf) console.log(`Arcing took ${performance.now() - segmentStart}ms`);
    segmentStart = performance.now();

    const arcSegments: ([number, [number, number][], number] | null)[] = [];

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

            arcSegments.push([circleGroups[i], points, i]);
        }
    }

    if (logPerf) console.log(`Found ${arcSegments.length} arc segments`);

    if (logPerf) console.log(`Segmenting took ${performance.now() - segmentStart}ms`);
    segmentStart = performance.now();

    const filteredSegmentsByGroup = arcSegments
        .reduce((acc, segment) => {
            const [group, arc] = segment!;
            acc[group] ??= [];
            acc[group].push(arc!);
            return acc;
        }, {} as {
            [key: number]: [number, number][][]
        });

    if (logPerf) console.log(`Grouping took ${performance.now() - segmentStart}ms`);
    segmentStart = performance.now();

    // Final step: Stitch together groups of line segments that are connected to make polygons
    const polygons: [number, number][][][] = [];

    for (const group in filteredSegmentsByGroup) {
        const segments = filteredSegmentsByGroup[group].filter(s => s.length > 1);
        const groupPolygon: [number, number][][] = [];
        let currentPolygon: [number, number][] = [];

        let i = 0;
        while (segments.length) {
            if (i++ > 1_000_000) {
                throw new Error("Infinite loop");
            }

            if (currentPolygon.length === 0) {
              const fixedSection = segments[0].slice();

              // If > 180 degrees latitude away from group polygon start, adjust longitudes in the segment by 360
              const firstPoint = groupPolygon[0]?.[0];
              if (firstPoint) {
                if (fixedSection[0][0] - firstPoint[0] > 180) {
                  for (let i = 0; i < fixedSection.length; i++) {
                    fixedSection[i] = [fixedSection[i][0] - 360, fixedSection[i][1]];
                  }
                } else if (fixedSection[0][0] - firstPoint[0] < -180) {
                  for (let i = 0; i < fixedSection.length; i++) {
                    fixedSection[i] = [fixedSection[i][0] + 360, fixedSection[i][1]];
                  }
                }
              }
                currentPolygon.push(...fixedSection);
                segments.splice(0, 1);
                continue;
            }

            let closestSegment: [number, number][] | null = null;
            let closestDistance = Infinity;
            for (const segment of segments) {
                const startDistance = distance(currentPolygon[currentPolygon.length - 1], segment[0]);

                if (startDistance < closestDistance) {
                    closestDistance = startDistance;
                    closestSegment = segment;
                }
            }
            if (closestSegment && closestDistance < maximumRadius * 0.2) {
                if (distance(currentPolygon[0], currentPolygon[currentPolygon.length - 1]) < closestDistance) {
                    groupPolygon.push([...currentPolygon, currentPolygon[0]]);
                    currentPolygon = [];
                } else {
                    const fixedSection = closestSegment.slice();

                    // If > 180 degrees latitude away from group/current polygon start, adjust longitudes in the segment by 360
                    const firstPoint = groupPolygon[0]?.[0] ?? currentPolygon[0];
                    if (fixedSection[0][0] - firstPoint[0] > 180) {
                      for (let i = 0; i < fixedSection.length; i++) {
                        fixedSection[i] = [fixedSection[i][0] - 360, fixedSection[i][1]];
                      }
                    } else if (fixedSection[0][0] - firstPoint[0] < -180) {
                      for (let i = 0; i < fixedSection.length; i++) {
                        fixedSection[i] = [fixedSection[i][0] + 360, fixedSection[i][1]];
                      }
                    }
                  currentPolygon.push(...fixedSection);
                  segments.splice(segments.indexOf(closestSegment), 1);
                }
            } else {
                if (currentPolygon.length > 2 && distance(currentPolygon[0], currentPolygon[currentPolygon.length - 1]) < maximumRadius * 0.5) {
                    groupPolygon.push([...currentPolygon, currentPolygon[0]]);
                    currentPolygon = [];
                } else if (currentPolygon.length) {
                    console.error("Could not find a segment to connect to", {currentPolygon, segments, dist: distance(currentPolygon[0], currentPolygon[currentPolygon.length - 1])});
                    throw new Error("Could not find a segment to connect to");
                }
            }
        }

        if (currentPolygon.length > 2 && distance(currentPolygon[0], currentPolygon[currentPolygon.length - 1]) < maximumRadius * 0.5) {
            groupPolygon.push([...currentPolygon, currentPolygon[0]]);
            currentPolygon = [];
        } else if (currentPolygon.length) {
            console.error("Could not find a segment to connect to", {circles, currentPolygon, segments, dist: distance(currentPolygon[0], currentPolygon[currentPolygon.length - 1])});
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

    if (logPerf) console.log(`Merging took ${performance.now() - segmentStart}ms`);
    segmentStart = performance.now();

    for (const poly of polygons) {
        if (poly.length === 0) continue;
        features.push(polygon(poly, properties));
    }

    if (logPerf) console.log(`Featuring took ${performance.now() - segmentStart}ms`);
    const end = performance.now();
    if (logPerf) console.log(`Time taken: ${end - start}ms`);

    return features;
}
