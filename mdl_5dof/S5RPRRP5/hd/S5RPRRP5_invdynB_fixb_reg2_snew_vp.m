% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RPRRP5_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP5_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP5_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP5_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:41:01
% EndTime: 2019-12-31 18:41:09
% DurationCPUTime: 6.67s
% Computational Cost: add. (16629->390), mult. (25393->519), div. (0->0), fcn. (14266->8), ass. (0->281)
t453 = qJD(4) ^ 2;
t440 = qJD(1) + qJD(3);
t438 = t440 ^ 2;
t447 = sin(qJ(4));
t442 = t447 ^ 2;
t515 = t442 * t438;
t419 = t453 + t515;
t450 = cos(qJ(4));
t424 = t450 * t438 * t447;
t414 = qJDD(4) - t424;
t496 = t450 * t414;
t368 = -t419 * t447 + t496;
t486 = qJD(4) * t440;
t481 = t450 * t486;
t439 = qJDD(1) + qJDD(3);
t505 = t447 * t439;
t395 = 0.2e1 * t481 + t505;
t448 = sin(qJ(3));
t451 = cos(qJ(3));
t319 = t368 * t448 + t395 * t451;
t322 = t368 * t451 - t395 * t448;
t445 = sin(pkin(8));
t446 = cos(pkin(8));
t268 = t319 * t446 + t322 * t445;
t272 = t319 * t445 - t322 * t446;
t449 = sin(qJ(1));
t452 = cos(qJ(1));
t220 = t268 * t452 - t272 * t449;
t571 = pkin(5) * t220;
t224 = t268 * t449 + t272 * t452;
t570 = pkin(5) * t224;
t569 = qJ(2) * t268;
t568 = pkin(1) * t268 + pkin(2) * t319 + pkin(7) * t368;
t508 = t447 * t414;
t362 = t419 * t450 + t508;
t567 = -pkin(1) * t362 - qJ(2) * t272;
t482 = t447 * t486;
t495 = t450 * t439;
t396 = -0.2e1 * t482 + t495;
t353 = t450 * t396;
t354 = t447 * t395;
t335 = -t353 + t354;
t443 = t450 ^ 2;
t409 = (t442 - t443) * t438;
t311 = t335 * t448 + t409 * t451;
t313 = t335 * t451 - t409 * t448;
t256 = t311 * t446 + t313 * t445;
t257 = t311 * t445 - t313 * t446;
t566 = t256 * t452 - t257 * t449;
t565 = t256 * t449 + t257 * t452;
t514 = t443 * t438;
t421 = -t453 + t514;
t366 = -t421 * t450 + t508;
t492 = t451 * t439;
t329 = t366 * t448 + t450 * t492;
t332 = t366 * t451 - t448 * t495;
t285 = t329 * t446 + t332 * t445;
t287 = t329 * t445 - t332 * t446;
t564 = t285 * t452 - t287 * t449;
t563 = t285 * t449 + t287 * t452;
t502 = t448 * t439;
t404 = t438 * t451 + t502;
t407 = t438 * t448 - t492;
t341 = t404 * t446 - t407 * t445;
t444 = g(3) - qJDD(2);
t379 = pkin(6) * t404 - t444 * t451;
t540 = pkin(6) * t407 - t444 * t448;
t283 = qJ(2) * t341 + t379 * t446 - t445 * t540;
t345 = t404 * t445 + t407 * t446;
t298 = t341 * t449 + t345 * t452;
t552 = qJ(2) * t345 + t379 * t445 + t446 * t540;
t562 = pkin(5) * t298 + t283 * t449 + t452 * t552;
t539 = t341 * t452 - t345 * t449;
t561 = pkin(5) * t539 + t283 * t452 - t449 * t552;
t560 = pkin(6) * t319;
t556 = -pkin(2) * t362 + pkin(6) * t322;
t426 = g(1) * t452 + g(2) * t449;
t454 = qJD(1) ^ 2;
t460 = -pkin(1) * t454 - t426;
t425 = g(1) * t449 - g(2) * t452;
t463 = qJDD(1) * pkin(1) + t425;
t352 = t445 * t463 + t446 * t460;
t350 = -pkin(2) * t454 + t352;
t456 = -t445 * t460 + t446 * t463;
t455 = qJDD(1) * pkin(2) + t456;
t304 = t448 * t350 - t451 * t455;
t305 = t350 * t451 + t448 * t455;
t478 = t304 * t448 + t305 * t451;
t246 = t304 * t451 - t305 * t448;
t512 = t446 * t246;
t207 = -t445 * t478 + t512;
t513 = t445 * t246;
t547 = t446 * t478 + t513;
t188 = t207 * t449 + t452 * t547;
t555 = t207 * t452 - t449 * t547;
t477 = t352 * t446 - t445 * t456;
t308 = -t352 * t445 - t446 * t456;
t490 = t452 * t308;
t548 = -t449 * t477 + t490;
t500 = t449 * t308;
t249 = t452 * t477 + t500;
t415 = qJDD(1) * t445 + t446 * t454;
t416 = qJDD(1) * t446 - t445 * t454;
t356 = -t415 * t449 + t416 * t452;
t385 = qJ(2) * t415 - t444 * t446;
t465 = -qJ(2) * t416 - t444 * t445;
t546 = -pkin(5) * t356 + t385 * t449 + t452 * t465;
t545 = 2 * qJD(5);
t542 = pkin(3) * t362;
t541 = pkin(7) * t362;
t529 = t415 * t452 + t416 * t449;
t537 = pkin(5) * t529 + t385 * t452 - t449 * t465;
t524 = pkin(4) * t450;
t468 = -qJ(5) * t447 - t524;
t394 = t468 * t440;
t293 = -pkin(3) * t438 + pkin(7) * t439 + t305;
t488 = -t293 * t450 + t444 * t447;
t466 = t394 * t440 * t450 + qJDD(4) * qJ(5) + (qJD(4) * t545) - t488;
t528 = t421 * t447 + t496;
t516 = t440 * t447;
t527 = t394 * t516 + qJDD(5);
t525 = pkin(1) * t444;
t422 = -t453 - t514;
t413 = qJDD(4) + t424;
t509 = t447 * t413;
t365 = t422 * t450 - t509;
t318 = t365 * t448 + t396 * t451;
t321 = t365 * t451 - t396 * t448;
t267 = t318 * t446 + t321 * t445;
t270 = -t318 * t445 + t321 * t446;
t219 = t267 * t452 + t270 * t449;
t523 = pkin(5) * t219;
t487 = t442 + t443;
t400 = t487 * t439;
t408 = t487 * t438;
t347 = t400 * t448 + t408 * t451;
t348 = t400 * t451 - t408 * t448;
t300 = t347 * t446 + t348 * t445;
t301 = -t347 * t445 + t348 * t446;
t240 = t300 * t452 + t301 * t449;
t522 = pkin(5) * t240;
t521 = pkin(6) * t318;
t520 = pkin(6) * t347;
t401 = t450 * t413;
t360 = t422 * t447 + t401;
t519 = pkin(7) * t360;
t518 = qJ(2) * t267;
t517 = qJ(2) * t300;
t292 = -pkin(3) * t439 - pkin(7) * t438 + t304;
t511 = t447 * t292;
t510 = t447 * t396;
t497 = t450 * t292;
t274 = t293 * t447 + t444 * t450;
t489 = t408 - t453;
t480 = -pkin(1) * t360 + qJ(2) * t270;
t479 = -pkin(2) * t360 + pkin(6) * t321;
t229 = t274 * t447 - t450 * t488;
t375 = -t425 * t449 - t426 * t452;
t474 = t448 * t424;
t473 = t451 * t424;
t472 = pkin(1) * t267 + pkin(2) * t318 + pkin(3) * t396 + pkin(7) * t365;
t471 = pkin(1) * t300 + pkin(2) * t347 + pkin(3) * t408 + pkin(7) * t400;
t259 = -pkin(3) * t360 + t274;
t418 = qJDD(1) * t452 - t449 * t454;
t469 = -pkin(5) * t418 - g(3) * t449;
t467 = pkin(4) * t447 - qJ(5) * t450;
t228 = t274 * t450 + t447 * t488;
t464 = t395 * t450 + t510;
t374 = t425 * t452 - t426 * t449;
t462 = t481 + t505;
t461 = -t482 + t495;
t459 = -qJDD(4) * pkin(4) + t274 + t527;
t458 = -t461 * pkin(4) + t292 + (-t462 - t481) * qJ(5);
t457 = t516 * t545 - t458;
t420 = t453 - t515;
t417 = qJDD(1) * t449 + t452 * t454;
t392 = -pkin(5) * t417 + g(3) * t452;
t391 = t467 * t439;
t390 = t487 * t486;
t373 = qJDD(4) * t448 + t390 * t451;
t372 = -qJDD(4) * t451 + t390 * t448;
t371 = -t442 * t486 + t450 * t462;
t370 = -t443 * t486 - t447 * t461;
t367 = -t420 * t447 + t401;
t361 = t420 * t450 + t509;
t339 = pkin(6) * t348;
t333 = t367 * t451 + t447 * t502;
t330 = t367 * t448 - t447 * t492;
t327 = t371 * t451 - t474;
t326 = t370 * t451 + t474;
t325 = t371 * t448 + t473;
t324 = t370 * t448 - t473;
t315 = -t372 * t445 + t373 * t446;
t314 = t372 * t446 + t373 * t445;
t303 = qJ(2) * t477 + t525;
t294 = qJ(2) * t301;
t289 = -t330 * t445 + t333 * t446;
t286 = t330 * t446 + t333 * t445;
t279 = -t325 * t445 + t327 * t446;
t278 = -t324 * t445 + t326 * t446;
t277 = t325 * t446 + t327 * t445;
t276 = t324 * t446 + t326 * t445;
t264 = t497 + t541;
t263 = t511 - t519;
t262 = -t314 * t449 + t315 * t452;
t261 = t314 * t452 + t315 * t449;
t260 = -t488 + t542;
t254 = t453 * qJ(5) - t459;
t253 = -pkin(4) * t453 + t466;
t252 = qJ(5) * t489 + t459;
t251 = pkin(4) * t489 + t466;
t250 = (pkin(4) * qJD(4) - (2 * qJD(5))) * t516 + t458;
t243 = (t396 - t482) * pkin(4) + t457;
t242 = -pkin(4) * t482 + qJ(5) * t395 + t457;
t241 = -t300 * t449 + t301 * t452;
t239 = pkin(2) * t444 + pkin(6) * t478;
t238 = pkin(5) * t241;
t237 = (-t422 - t453) * qJ(5) + (-qJDD(4) - t413) * pkin(4) + t259 + t527;
t236 = -t542 - qJ(5) * t414 + (-t419 + t453) * pkin(4) - t466;
t235 = -t286 * t449 + t289 * t452;
t234 = t286 * t452 + t289 * t449;
t233 = -t277 * t449 + t279 * t452;
t232 = -t276 * t449 + t278 * t452;
t231 = t277 * t452 + t279 * t449;
t230 = t276 * t452 + t278 * t449;
t226 = -pkin(4) * t354 + t242 * t450 - t541;
t225 = qJ(5) * t353 - t243 * t447 - t519;
t222 = -t267 * t449 + t270 * t452;
t218 = pkin(5) * t222;
t217 = t228 * t451 - t520;
t216 = t228 * t448 + t339;
t215 = t253 * t450 - t254 * t447;
t214 = t253 * t447 + t254 * t450;
t213 = t229 * t451 + t292 * t448;
t212 = t229 * t448 - t292 * t451;
t211 = -t251 * t447 + t252 * t450;
t210 = -t260 * t448 + t264 * t451 + t560;
t209 = -t259 * t448 + t263 * t451 - t521;
t204 = t260 * t451 + t264 * t448 - t556;
t203 = t259 * t451 + t263 * t448 + t479;
t202 = t211 * t451 - t391 * t448 - t520;
t201 = t211 * t448 + t391 * t451 + t339;
t200 = t215 * t451 + t250 * t448;
t199 = t215 * t448 - t250 * t451;
t198 = t225 * t451 - t237 * t448 - t521;
t197 = t226 * t451 - t236 * t448 - t560;
t196 = t225 * t448 + t237 * t451 + t479;
t195 = t226 * t448 + t236 * t451 + t556;
t194 = -pkin(3) * t214 - pkin(4) * t254 - qJ(5) * t253;
t193 = -pkin(7) * t214 + t250 * t467;
t192 = -t216 * t445 + t217 * t446 - t517;
t191 = t216 * t446 + t217 * t445 + t294;
t190 = -t212 * t445 + t213 * t446;
t189 = t212 * t446 + t213 * t445;
t186 = pkin(6) * t512 + qJ(2) * t207 - t239 * t445;
t185 = pkin(6) * t513 + qJ(2) * t547 + t239 * t446 + t525;
t184 = -pkin(6) * t212 - (pkin(3) * t448 - pkin(7) * t451) * t228;
t183 = -t204 * t445 + t210 * t446 + t569;
t182 = -t203 * t445 + t209 * t446 - t518;
t181 = t204 * t446 + t210 * t445 - t567;
t180 = t203 * t446 + t209 * t445 + t480;
t179 = -t201 * t445 + t202 * t446 - t517;
t178 = t201 * t446 + t202 * t445 + t294;
t177 = -t199 * t445 + t200 * t446;
t176 = t199 * t446 + t200 * t445;
t175 = pkin(6) * t213 - (-pkin(3) * t451 - pkin(7) * t448 - pkin(2)) * t228;
t174 = -t196 * t445 + t198 * t446 - t518;
t173 = -t195 * t445 + t197 * t446 - t569;
t172 = t196 * t446 + t198 * t445 + t480;
t171 = t195 * t446 + t197 * t445 + t567;
t170 = -t189 * t449 + t190 * t452;
t169 = t189 * t452 + t190 * t449;
t168 = -pkin(6) * t199 + t193 * t451 - t194 * t448;
t167 = -t176 * t449 + t177 * t452;
t166 = t176 * t452 + t177 * t449;
t165 = -pkin(2) * t214 + pkin(6) * t200 + t193 * t448 + t194 * t451;
t164 = -qJ(2) * t189 - t175 * t445 + t184 * t446;
t163 = pkin(1) * t228 + qJ(2) * t190 + t175 * t446 + t184 * t445;
t162 = -qJ(2) * t176 - t165 * t445 + t168 * t446;
t161 = -pkin(1) * t214 + qJ(2) * t177 + t165 * t446 + t168 * t445;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t417, -t418, 0, t375, 0, 0, 0, 0, 0, 0, -t529, -t356, 0, t249, 0, 0, 0, 0, 0, 0, -t539, t298, 0, t188, 0, 0, 0, 0, 0, 0, t222, t224, t241, t170, 0, 0, 0, 0, 0, 0, t222, t241, -t224, t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t418, -t417, 0, t374, 0, 0, 0, 0, 0, 0, t356, -t529, 0, -t548, 0, 0, 0, 0, 0, 0, -t298, -t539, 0, -t555, 0, 0, 0, 0, 0, 0, t219, -t220, t240, t169, 0, 0, 0, 0, 0, 0, t219, t240, t220, t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t444, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t444, 0, 0, 0, 0, 0, 0, t360, -t362, 0, -t228, 0, 0, 0, 0, 0, 0, t360, 0, t362, t214; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t418, 0, -t417, 0, t469, -t392, -t374, -pkin(5) * t374, 0, 0, t356, 0, -t529, 0, t546, t537, t548, pkin(5) * t548 + qJ(2) * t490 - t303 * t449, 0, 0, -t298, 0, -t539, 0, t562, t561, t555, pkin(5) * t555 - t185 * t449 + t186 * t452, t233, t565, t235, t232, t563, t262, -t180 * t449 + t182 * t452 - t523, -t181 * t449 + t183 * t452 + t571, -t191 * t449 + t192 * t452 - t522, -pkin(5) * t169 - t163 * t449 + t164 * t452, t233, t235, -t565, t262, -t563, t232, -t172 * t449 + t174 * t452 - t523, -t178 * t449 + t179 * t452 - t522, -t171 * t449 + t173 * t452 - t571, -pkin(5) * t166 - t161 * t449 + t162 * t452; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t417, 0, t418, 0, t392, t469, t375, pkin(5) * t375, 0, 0, t529, 0, t356, 0, -t537, t546, t249, pkin(5) * t249 + qJ(2) * t500 + t303 * t452, 0, 0, t539, 0, -t298, 0, -t561, t562, t188, pkin(5) * t188 + t185 * t452 + t186 * t449, t231, -t566, t234, t230, -t564, t261, t180 * t452 + t182 * t449 + t218, t181 * t452 + t183 * t449 + t570, t191 * t452 + t192 * t449 + t238, pkin(5) * t170 + t163 * t452 + t164 * t449, t231, t234, t566, t261, t564, t230, t172 * t452 + t174 * t449 + t218, t178 * t452 + t179 * t449 + t238, t171 * t452 + t173 * t449 - t570, pkin(5) * t167 + t161 * t452 + t162 * t449; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t425, t426, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * t416 + t456, -pkin(1) * t415 - t352, 0, -pkin(1) * t308, 0, 0, 0, 0, 0, t439, -pkin(1) * t345 - pkin(2) * t407 - t304, -pkin(1) * t341 - pkin(2) * t404 - t305, 0, -pkin(1) * t207 - pkin(2) * t246, t354, t464, t361, t353, t528, 0, t472 - t497, -pkin(3) * t395 + t511 - t568, t229 + t471, pkin(1) * t189 + pkin(2) * t212 - pkin(3) * t292 + pkin(7) * t229, t354, t361, -t464, 0, -t528, t353, qJ(5) * t510 + t243 * t450 + t472, t251 * t450 + t252 * t447 + t471, t447 * t242 + (pkin(3) + t524) * t395 + t568, pkin(1) * t176 + pkin(2) * t199 + pkin(7) * t215 + (-pkin(3) + t468) * t250;];
tauB_reg = t1;
