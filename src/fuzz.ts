import {circleUnion} from './index';
import seedrandom from 'seed-random';

function fuzz360Edge() {
  let success = 0;
  let fail = 0;
  for (let i = 0; i < 1000; i++) {
    const r = seedrandom(i.toString() + "seed4");
    let data: [number, number, number][] = [];
    for (let n = 0; n < 500; n++) {
      data.push([(r() * 0.8 + 360), r() * 2, r() * 20]);
    }
    try {
      circleUnion(data, {});
      success++;
    } catch (err) {
      console.info("Failed on line: ", i)
      fail++;
    }
  }
  console.log("Success: ", success);
  console.log("Fail: ", fail);
}

function fuzzRange() {
  let success = 0;
  let fail = 0;
  for (let i = 0; i < 1000; i++) {
    const r = seedrandom(i.toString() + "seed4");
    const lngConst = r() * 720;
    const latConst = r() * 2;
    const radiusConst = r() * 20;
    let data: [number, number, number][] = [];
    for (let n = 0; n < 1000; n++) {
      data.push([(r() - 0.5) * lngConst, r() * latConst, r() * radiusConst]);
    }
    try {
      circleUnion(data, {});
      success++;
    } catch (err) {
      console.info("Failed on line: ", i)
      fail++;
    }
  }
  console.log("Success: ", success);
  console.log("Fail: ", fail);
}

fuzz360Edge();
fuzzRange();
