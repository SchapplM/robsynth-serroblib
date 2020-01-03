% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RRRR3_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR3_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR3_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR3_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:42
% EndTime: 2019-12-31 17:24:48
% DurationCPUTime: 5.39s
% Computational Cost: add. (24012->441), mult. (52723->668), div. (0->0), fcn. (36596->8), ass. (0->311)
t481 = sin(qJ(4));
t482 = sin(qJ(3));
t486 = cos(qJ(3));
t487 = cos(qJ(2));
t483 = sin(qJ(2));
t514 = qJD(1) * t483;
t439 = -t486 * t487 * qJD(1) + t482 * t514;
t441 = (t482 * t487 + t483 * t486) * qJD(1);
t485 = cos(qJ(4));
t404 = t485 * t439 + t481 * t441;
t406 = -t481 * t439 + t485 * t441;
t354 = t406 * t404;
t477 = qJDD(2) + qJDD(3);
t499 = qJDD(4) + t477;
t543 = -t354 + t499;
t547 = t481 * t543;
t412 = t441 * t439;
t541 = -t412 + t477;
t546 = t482 * t541;
t545 = t485 * t543;
t544 = t486 * t541;
t478 = qJD(2) + qJD(3);
t433 = t478 * t439;
t511 = qJD(1) * qJD(2);
t502 = t487 * t511;
t510 = t483 * qJDD(1);
t450 = t502 + t510;
t472 = t487 * qJDD(1);
t503 = t483 * t511;
t451 = t472 - t503;
t492 = t439 * qJD(3) - t486 * t450 - t482 * t451;
t366 = -t433 + t492;
t471 = qJD(4) + t478;
t395 = t471 * t404;
t500 = t482 * t450 - t486 * t451;
t385 = -t441 * qJD(3) - t500;
t491 = t404 * qJD(4) - t481 * t385 + t485 * t492;
t542 = -t395 - t491;
t540 = -t433 - t492;
t501 = -t485 * t385 - t481 * t492;
t310 = (qJD(4) - t471) * t406 + t501;
t362 = (qJD(3) - t478) * t441 + t500;
t402 = t404 ^ 2;
t403 = t406 ^ 2;
t437 = t439 ^ 2;
t438 = t441 ^ 2;
t469 = t471 ^ 2;
t476 = t478 ^ 2;
t539 = t471 * t481;
t538 = t471 * t485;
t537 = t478 * t482;
t536 = t478 * t486;
t479 = t483 ^ 2;
t490 = qJD(1) ^ 2;
t535 = t479 * t490;
t480 = t487 ^ 2;
t474 = t480 * t490;
t484 = sin(qJ(1));
t488 = cos(qJ(1));
t460 = t484 * g(1) - t488 * g(2);
t493 = qJDD(1) * pkin(1) + t460;
t494 = qJD(2) * pkin(2) - pkin(6) * t514;
t388 = t451 * pkin(2) - t494 * t514 + (pkin(6) * t480 + pkin(5)) * t490 + t493;
t496 = t478 * pkin(3) - t441 * pkin(7);
t321 = t385 * pkin(3) + t437 * pkin(7) - t441 * t496 + t388;
t534 = t481 * t321;
t349 = t354 + t499;
t533 = t481 * t349;
t524 = t483 * t490;
t461 = t488 * g(1) + t484 * g(2);
t443 = -t490 * pkin(1) + qJDD(1) * pkin(5) - t461;
t527 = t483 * t443;
t376 = qJDD(2) * pkin(2) - t450 * pkin(6) - t527 + (pkin(2) * t524 + pkin(6) * t511 - g(3)) * t487;
t427 = -t483 * g(3) + t487 * t443;
t380 = -pkin(2) * t474 + t451 * pkin(6) - qJD(2) * t494 + t427;
t342 = -t486 * t376 + t482 * t380;
t296 = t541 * pkin(3) + t366 * pkin(7) - t342;
t343 = t482 * t376 + t486 * t380;
t300 = -t437 * pkin(3) + t385 * pkin(7) - t478 * t496 + t343;
t258 = -t485 * t296 + t481 * t300;
t259 = t481 * t296 + t485 * t300;
t232 = -t485 * t258 + t481 * t259;
t532 = t482 * t232;
t531 = t482 * t388;
t409 = t412 + t477;
t530 = t482 * t409;
t292 = -t486 * t342 + t482 * t343;
t529 = t483 * t292;
t442 = t490 * pkin(5) + t493;
t528 = t483 * t442;
t467 = t487 * t524;
t458 = qJDD(2) + t467;
t526 = t483 * t458;
t459 = qJDD(2) - t467;
t525 = t483 * t459;
t523 = t485 * t321;
t522 = t485 * t349;
t521 = t486 * t232;
t520 = t486 * t388;
t519 = t486 * t409;
t518 = t487 * t292;
t517 = t487 * t442;
t516 = t487 * t459;
t515 = t479 + t480;
t509 = t484 * qJDD(1);
t508 = t488 * qJDD(1);
t507 = t484 * t354;
t506 = t484 * t412;
t505 = t488 * t354;
t504 = t488 * t412;
t233 = t481 * t258 + t485 * t259;
t293 = t482 * t342 + t486 * t343;
t426 = t487 * g(3) + t527;
t379 = t483 * t426 + t487 * t427;
t418 = -t484 * t460 - t488 * t461;
t498 = t484 * t467;
t497 = t488 * t467;
t455 = -t484 * t490 + t508;
t495 = -pkin(4) * t455 - t484 * g(3);
t378 = t487 * t426 - t483 * t427;
t417 = t488 * t460 - t484 * t461;
t489 = qJD(2) ^ 2;
t465 = -t474 - t489;
t464 = t474 - t489;
t463 = -t489 - t535;
t462 = t489 - t535;
t457 = t474 - t535;
t456 = t474 + t535;
t454 = t488 * t490 + t509;
t453 = t515 * qJDD(1);
t452 = t472 - 0.2e1 * t503;
t449 = 0.2e1 * t502 + t510;
t447 = t487 * t458;
t446 = t515 * t511;
t436 = -pkin(4) * t454 + t488 * g(3);
t431 = -t438 + t476;
t430 = t437 - t476;
t429 = t487 * t450 - t479 * t511;
t428 = -t483 * t451 - t480 * t511;
t425 = -t438 - t476;
t424 = -t483 * t463 - t516;
t423 = -t483 * t462 + t447;
t422 = t487 * t465 - t526;
t421 = t487 * t464 - t525;
t420 = t487 * t463 - t525;
t419 = t483 * t465 + t447;
t415 = t488 * t453 - t484 * t456;
t414 = t484 * t453 + t488 * t456;
t413 = -t483 * t449 + t487 * t452;
t411 = -t438 + t437;
t407 = -t476 - t437;
t400 = t488 * t424 + t484 * t449;
t399 = t488 * t422 - t484 * t452;
t398 = t484 * t424 - t488 * t449;
t397 = t484 * t422 + t488 * t452;
t394 = -t403 + t469;
t393 = t402 - t469;
t392 = (-t439 * t486 + t441 * t482) * t478;
t391 = (-t439 * t482 - t441 * t486) * t478;
t390 = -pkin(5) * t420 - t517;
t389 = -pkin(5) * t419 - t528;
t387 = -t437 - t438;
t384 = -pkin(1) * t420 + t427;
t383 = -pkin(1) * t419 + t426;
t382 = -t403 - t469;
t372 = t486 * t430 - t530;
t371 = -t482 * t431 + t544;
t370 = t482 * t430 + t519;
t369 = t486 * t431 + t546;
t368 = -t482 * t425 - t519;
t367 = t486 * t425 - t530;
t361 = (qJD(3) + t478) * t441 + t500;
t360 = -t441 * t537 - t486 * t492;
t359 = t441 * t536 - t482 * t492;
t358 = -t482 * t385 + t439 * t536;
t357 = t486 * t385 + t439 * t537;
t356 = t488 * t379 - t484 * t442;
t355 = t484 * t379 + t488 * t442;
t353 = -t403 + t402;
t352 = t486 * t407 - t546;
t351 = t482 * t407 + t544;
t347 = -t469 - t402;
t346 = (-t404 * t485 + t406 * t481) * t471;
t345 = (-t404 * t481 - t406 * t485) * t471;
t344 = -t483 * t391 + t487 * t392;
t341 = -t402 - t403;
t339 = -pkin(6) * t367 - t520;
t338 = -pkin(6) * t351 - t531;
t337 = -t483 * t370 + t487 * t372;
t336 = -t483 * t369 + t487 * t371;
t335 = t485 * t393 - t533;
t334 = -t481 * t394 + t545;
t333 = t481 * t393 + t522;
t332 = t485 * t394 + t547;
t331 = -t481 * t382 - t522;
t330 = t485 * t382 - t533;
t329 = -t483 * t367 + t487 * t368;
t328 = t487 * t367 + t483 * t368;
t327 = -t362 * t486 - t482 * t366;
t326 = -t486 * t361 - t482 * t540;
t325 = -t362 * t482 + t486 * t366;
t324 = -t482 * t361 + t486 * t540;
t322 = -t406 * qJD(4) - t501;
t320 = -t483 * t359 + t487 * t360;
t319 = -t483 * t357 + t487 * t358;
t318 = -t483 * t351 + t487 * t352;
t317 = t487 * t351 + t483 * t352;
t316 = t485 * t347 - t547;
t315 = t481 * t347 + t545;
t314 = -t395 + t491;
t309 = (qJD(4) + t471) * t406 + t501;
t308 = -pkin(2) * t540 + pkin(6) * t368 - t531;
t307 = -t482 * t345 + t486 * t346;
t306 = t486 * t345 + t482 * t346;
t305 = -t406 * t539 - t485 * t491;
t304 = t406 * t538 - t481 * t491;
t303 = -t481 * t322 + t404 * t538;
t302 = t485 * t322 + t404 * t539;
t301 = -pkin(2) * t361 + pkin(6) * t352 + t520;
t298 = t488 * t329 + t484 * t540;
t297 = t484 * t329 - t488 * t540;
t291 = t488 * t318 + t484 * t361;
t290 = t484 * t318 - t488 * t361;
t289 = -t482 * t333 + t486 * t335;
t288 = -t482 * t332 + t486 * t334;
t287 = t486 * t333 + t482 * t335;
t286 = t486 * t332 + t482 * t334;
t285 = -t482 * t330 + t486 * t331;
t284 = t486 * t330 + t482 * t331;
t283 = pkin(2) * t388 + pkin(6) * t293;
t282 = -t483 * t325 + t487 * t327;
t281 = -t483 * t324 + t487 * t326;
t280 = t487 * t325 + t483 * t327;
t279 = -pkin(7) * t330 - t523;
t278 = -pkin(1) * t328 - pkin(2) * t367 + t343;
t277 = -pkin(7) * t315 - t534;
t276 = t488 * t282 + t484 * t387;
t275 = t484 * t282 - t488 * t387;
t274 = -t482 * t315 + t486 * t316;
t273 = t486 * t315 + t482 * t316;
t272 = -pkin(1) * t317 - pkin(2) * t351 + t342;
t271 = -t310 * t485 - t481 * t314;
t270 = -t485 * t309 - t481 * t542;
t269 = -t310 * t481 + t485 * t314;
t268 = -t481 * t309 + t485 * t542;
t267 = -pkin(6) * t325 - t292;
t266 = -t483 * t306 + t487 * t307;
t265 = -t482 * t304 + t486 * t305;
t264 = -t482 * t302 + t486 * t303;
t263 = t486 * t304 + t482 * t305;
t262 = t486 * t302 + t482 * t303;
t261 = -pkin(2) * t387 + pkin(6) * t327 + t293;
t260 = -pkin(1) * t280 - pkin(2) * t325;
t256 = -pkin(5) * t328 - t483 * t308 + t487 * t339;
t255 = t487 * t293 - t529;
t254 = t483 * t293 + t518;
t253 = -pkin(3) * t542 + pkin(7) * t331 - t534;
t252 = -pkin(5) * t317 - t483 * t301 + t487 * t338;
t251 = -pkin(3) * t309 + pkin(7) * t316 + t523;
t250 = t488 * t255 - t484 * t388;
t249 = t484 * t255 + t488 * t388;
t248 = -t483 * t287 + t487 * t289;
t247 = -t483 * t286 + t487 * t288;
t246 = -t483 * t284 + t487 * t285;
t245 = t487 * t284 + t483 * t285;
t244 = -t483 * t273 + t487 * t274;
t243 = t487 * t273 + t483 * t274;
t242 = -pkin(1) * t254 - pkin(2) * t292;
t241 = -t482 * t269 + t486 * t271;
t240 = -t482 * t268 + t486 * t270;
t239 = t486 * t269 + t482 * t271;
t238 = t486 * t268 + t482 * t270;
t237 = t488 * t246 + t484 * t542;
t236 = t484 * t246 - t488 * t542;
t235 = -t483 * t263 + t487 * t265;
t234 = -t483 * t262 + t487 * t264;
t231 = t488 * t244 + t484 * t309;
t230 = t484 * t244 - t488 * t309;
t229 = -pkin(5) * t254 - pkin(6) * t518 - t483 * t283;
t228 = pkin(3) * t321 + pkin(7) * t233;
t227 = -pkin(6) * t284 - t482 * t253 + t486 * t279;
t226 = -pkin(5) * t280 - t483 * t261 + t487 * t267;
t225 = -pkin(6) * t273 - t482 * t251 + t486 * t277;
t224 = -pkin(2) * t542 + pkin(6) * t285 + t486 * t253 + t482 * t279;
t223 = -pkin(7) * t269 - t232;
t222 = -pkin(3) * t341 + pkin(7) * t271 + t233;
t221 = -pkin(2) * t309 + pkin(6) * t274 + t486 * t251 + t482 * t277;
t220 = -pkin(1) * t245 - pkin(2) * t284 - pkin(3) * t330 + t259;
t219 = -t483 * t239 + t487 * t241;
t218 = -t483 * t238 + t487 * t240;
t217 = t487 * t239 + t483 * t241;
t216 = t488 * t219 + t484 * t341;
t215 = t484 * t219 - t488 * t341;
t214 = -pkin(1) * t243 - pkin(2) * t273 - pkin(3) * t315 + t258;
t213 = t486 * t233 - t532;
t212 = t482 * t233 + t521;
t211 = -pkin(1) * t217 - pkin(2) * t239 - pkin(3) * t269;
t210 = -pkin(5) * t245 - t483 * t224 + t487 * t227;
t209 = -pkin(5) * t243 - t483 * t221 + t487 * t225;
t208 = -pkin(6) * t239 - t482 * t222 + t486 * t223;
t207 = -pkin(2) * t341 + pkin(6) * t241 + t486 * t222 + t482 * t223;
t206 = -t483 * t212 + t487 * t213;
t205 = t487 * t212 + t483 * t213;
t204 = -pkin(6) * t212 - pkin(7) * t521 - t482 * t228;
t203 = t488 * t206 - t484 * t321;
t202 = t484 * t206 + t488 * t321;
t201 = pkin(2) * t321 + pkin(6) * t213 - pkin(7) * t532 + t486 * t228;
t200 = -pkin(1) * t205 - pkin(2) * t212 - pkin(3) * t232;
t199 = -pkin(5) * t217 - t483 * t207 + t487 * t208;
t198 = -pkin(5) * t205 - t483 * t201 + t487 * t204;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t454, -t455, 0, t418, 0, 0, 0, 0, 0, 0, t399, t400, t415, t356, 0, 0, 0, 0, 0, 0, t291, t298, t276, t250, 0, 0, 0, 0, 0, 0, t231, t237, t216, t203; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t455, -t454, 0, t417, 0, 0, 0, 0, 0, 0, t397, t398, t414, t355, 0, 0, 0, 0, 0, 0, t290, t297, t275, t249, 0, 0, 0, 0, 0, 0, t230, t236, t215, t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t419, t420, 0, -t378, 0, 0, 0, 0, 0, 0, t317, t328, t280, t254, 0, 0, 0, 0, 0, 0, t243, t245, t217, t205; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t455, 0, -t454, 0, t495, -t436, -t417, -pkin(4) * t417, t488 * t429 - t498, t488 * t413 - t484 * t457, t488 * t423 + t483 * t509, t488 * t428 + t498, t488 * t421 + t472 * t484, t484 * qJDD(2) + t488 * t446, -pkin(4) * t397 - t484 * t383 + t488 * t389, -pkin(4) * t398 - t484 * t384 + t488 * t390, -pkin(4) * t414 + t488 * t378, -pkin(4) * t355 - (pkin(1) * t484 - pkin(5) * t488) * t378, t488 * t320 + t506, t488 * t281 - t484 * t411, t488 * t336 - t484 * t366, t488 * t319 - t506, t488 * t337 - t484 * t362, t488 * t344 + t484 * t477, -pkin(4) * t290 + t488 * t252 - t484 * t272, -pkin(4) * t297 + t488 * t256 - t484 * t278, -pkin(4) * t275 + t488 * t226 - t484 * t260, -pkin(4) * t249 + t488 * t229 - t484 * t242, t488 * t235 + t507, t488 * t218 - t484 * t353, t488 * t247 - t484 * t314, t488 * t234 - t507, t488 * t248 - t484 * t310, t488 * t266 + t484 * t499, -pkin(4) * t230 + t488 * t209 - t484 * t214, -pkin(4) * t236 + t488 * t210 - t484 * t220, -pkin(4) * t215 + t488 * t199 - t484 * t211, -pkin(4) * t202 + t488 * t198 - t484 * t200; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t454, 0, t455, 0, t436, t495, t418, pkin(4) * t418, t484 * t429 + t497, t484 * t413 + t488 * t457, t484 * t423 - t483 * t508, t484 * t428 - t497, t484 * t421 - t487 * t508, -t488 * qJDD(2) + t484 * t446, pkin(4) * t399 + t488 * t383 + t484 * t389, pkin(4) * t400 + t488 * t384 + t484 * t390, pkin(4) * t415 + t484 * t378, pkin(4) * t356 - (-pkin(1) * t488 - pkin(5) * t484) * t378, t484 * t320 - t504, t484 * t281 + t488 * t411, t484 * t336 + t488 * t366, t484 * t319 + t504, t484 * t337 + t488 * t362, t484 * t344 - t488 * t477, pkin(4) * t291 + t484 * t252 + t488 * t272, pkin(4) * t298 + t484 * t256 + t488 * t278, pkin(4) * t276 + t484 * t226 + t488 * t260, pkin(4) * t250 + t484 * t229 + t488 * t242, t484 * t235 - t505, t484 * t218 + t488 * t353, t484 * t247 + t488 * t314, t484 * t234 + t505, t484 * t248 + t488 * t310, t484 * t266 - t488 * t499, pkin(4) * t231 + t484 * t209 + t488 * t214, pkin(4) * t237 + t484 * t210 + t488 * t220, pkin(4) * t216 + t484 * t199 + t488 * t211, pkin(4) * t203 + t484 * t198 + t488 * t200; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t460, t461, 0, 0, (t450 + t502) * t483, t487 * t449 + t483 * t452, t487 * t462 + t526, (t451 - t503) * t487, t483 * t464 + t516, 0, pkin(1) * t452 + pkin(5) * t422 + t517, -pkin(1) * t449 + pkin(5) * t424 - t528, pkin(1) * t456 + pkin(5) * t453 + t379, pkin(1) * t442 + pkin(5) * t379, t487 * t359 + t483 * t360, t487 * t324 + t483 * t326, t487 * t369 + t483 * t371, t487 * t357 + t483 * t358, t487 * t370 + t483 * t372, t487 * t391 + t483 * t392, -pkin(1) * t361 + pkin(5) * t318 + t487 * t301 + t483 * t338, -pkin(1) * t540 + pkin(5) * t329 + t487 * t308 + t483 * t339, -pkin(1) * t387 + pkin(5) * t282 + t487 * t261 + t483 * t267, pkin(1) * t388 + pkin(5) * t255 - pkin(6) * t529 + t487 * t283, t487 * t263 + t483 * t265, t487 * t238 + t483 * t240, t487 * t286 + t483 * t288, t487 * t262 + t483 * t264, t487 * t287 + t483 * t289, t487 * t306 + t483 * t307, -pkin(1) * t309 + pkin(5) * t244 + t487 * t221 + t483 * t225, -pkin(1) * t542 + pkin(5) * t246 + t487 * t224 + t483 * t227, -pkin(1) * t341 + pkin(5) * t219 + t487 * t207 + t483 * t208, pkin(1) * t321 + pkin(5) * t206 + t487 * t201 + t483 * t204;];
tauB_reg = t1;