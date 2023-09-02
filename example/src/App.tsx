import './App.css'
import MapGL, {Layer, Source} from 'react-map-gl';
import seedrandom from 'seed-random';
import {circleUnion} from "../../src";
import {circle} from "@turf/turf";

let data: [number, number, number][] = [];

const r = seedrandom("hi");

for (let i = 0; i < 439; i++) {
    data.push([(r() * 2 + 179) % 360, r() * 10, r() * 40]);
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
                          // longitude: 180,
                          // latitude: 5,
                          // zoom: 5
                          // longitude: 0,
                          // latitude: 1.6,
                          // zoom: 10
                          longitude: 0,
                          latitude: 0,
                          zoom: 1
                        }}
                        mapboxAccessToken="pk.eyJ1Ijoic29oY2FoIiwiYSI6ImNqeWVqcm8wdTAxc2MzaXFpa282Yzd2aHEifQ.afYbt2sVMZ-kbwdx5_PekQ"
                        mapStyle="mapbox://styles/mapbox/streets-v11"
                        fog={{}}
                        light={{}}
                    >
                        <Source type="geojson" data={{
                            type: "FeatureCollection",
                            features: circleUnion(data, {}, true).map(i => ({
                                ...i,
                                properties: {
                                    color: `#${Math.floor(Math.random() * 0x1000000).toString(16).padStart(6, "0")}`,
                                },
                            }))
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
                            {/*<Layer type="circle"*/}
                            {/*       paint={{*/}
                            {/*           "circle-color": ["get", "color"],*/}
                            {/*            "circle-radius": 5*/}
                            {/*       }}/>*/}
                        </Source>
                      <Source type="geojson" data={{
                        type: "FeatureCollection",
                        features: data.map(i => circle([i[0], i[1]], i[2], {
                          units: "kilometers"
                        }))
                      }}>
                        <Layer type="fill"
                               paint={{
                                 "fill-color": "#f00",
                                 "fill-opacity": 0.1
                               }}
                        />
                      </Source>
                    </MapGL>
                </div>
            </div>
        </>
    )
}

export default App
