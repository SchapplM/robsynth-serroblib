% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR15
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR15_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR15_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:38:58
% EndTime: 2019-02-26 22:39:02
% DurationCPUTime: 4.50s
% Computational Cost: add. (16749->210), mult. (50114->407), div. (948->12), fcn. (63586->15), ass. (0->176)
t480 = cos(pkin(6));
t481 = sin(qJ(2));
t420 = t480 * t481;
t379 = cos(qJ(2));
t482 = sin(qJ(1));
t435 = t482 * t379;
t483 = cos(qJ(1));
t362 = t420 * t483 + t435;
t393 = t435 * t480 + t481 * t483;
t352 = qJD(1) * t393 + qJD(2) * t362;
t371 = t482 * t481;
t406 = t482 * t420;
t421 = t480 * t483;
t432 = t483 * qJD(1);
t353 = -qJD(1) * t406 - qJD(2) * t371 + (qJD(2) * t421 + t432) * t379;
t374 = cos(pkin(7));
t376 = sin(qJ(3));
t378 = cos(qJ(3));
t372 = sin(pkin(7));
t373 = sin(pkin(6));
t436 = t373 * t482;
t422 = qJD(1) * t436;
t414 = t372 * t422;
t405 = -t379 * t421 + t371;
t394 = t405 * t374;
t437 = t373 * t483;
t426 = t372 * t437;
t388 = t394 + t426;
t486 = t362 * t376 + t378 * t388;
t299 = t486 * qJD(3) - (-t352 * t374 + t414) * t376 - t353 * t378;
t461 = t362 * t378;
t341 = t376 * t388 - t461;
t375 = sin(qJ(4));
t377 = cos(qJ(4));
t390 = t372 * t405 - t374 * t437;
t319 = t341 * t375 + t377 * t390;
t397 = t352 * t372 + t374 * t422;
t293 = qJD(4) * t319 - t299 * t377 + t375 * t397;
t320 = t341 * t377 - t375 * t390;
t496 = qJD(4) * t320 + t299 * t375 + t377 * t397;
t434 = t481 * t376;
t453 = t378 * t379;
t402 = -t374 * t434 + t453;
t403 = t374 * t453 - t434;
t431 = t480 * t372;
t419 = qJD(3) * t431;
t327 = t378 * t419 + (qJD(2) * t402 + qJD(3) * t403) * t373;
t433 = t481 * t378;
t454 = t376 * t379;
t404 = t374 * t454 + t433;
t357 = t373 * t404 + t376 * t431;
t457 = t372 * t379;
t361 = -t373 * t457 + t374 * t480;
t332 = t357 * t377 + t361 * t375;
t438 = t372 * t481;
t427 = t373 * t438;
t413 = qJD(2) * t427;
t309 = qJD(4) * t332 + t327 * t375 - t377 * t413;
t331 = t357 * t375 - t361 * t377;
t329 = 0.1e1 / t331 ^ 2;
t491 = t309 * t329;
t328 = 0.1e1 / t331;
t356 = t373 * t403 + t378 * t431;
t467 = t319 * t329;
t408 = t328 * t486 - t356 * t467;
t490 = t375 * t408;
t424 = t372 * t436;
t386 = -t374 * t393 + t424;
t392 = -t379 * t483 + t406;
t343 = t376 * t386 - t378 * t392;
t350 = t405 * qJD(1) + t392 * qJD(2);
t351 = qJD(1) * t362 + qJD(2) * t393;
t423 = t373 * t432;
t415 = t372 * t423;
t455 = t374 * t378;
t295 = qJD(3) * t343 - t350 * t455 - t351 * t376 - t378 * t415;
t391 = t393 * t378;
t460 = t392 * t376;
t342 = t374 * t391 - t378 * t424 - t460;
t335 = 0.1e1 / t342;
t296 = -t351 * t378 + (t350 * t374 + t415) * t376 + (t378 * t386 + t460) * qJD(3);
t387 = t372 * t393 + t374 * t436;
t321 = t343 * t375 - t377 * t387;
t398 = -t350 * t372 + t374 * t423;
t291 = -qJD(4) * t321 + t296 * t377 + t375 * t398;
t322 = t343 * t377 + t375 * t387;
t315 = t322 ^ 2;
t336 = 0.1e1 / t342 ^ 2;
t308 = t315 * t336 + 0.1e1;
t337 = t335 * t336;
t466 = t322 * t336;
t477 = (-t295 * t315 * t337 + t291 * t466) / t308 ^ 2;
t445 = t335 * t477;
t306 = 0.1e1 / t308;
t470 = t306 * t336;
t489 = -t295 * t470 - 0.2e1 * t445;
t446 = 0.2e1 * t322 * t337;
t449 = 0.2e1 * t477;
t488 = t295 * t306 * t446 - t291 * t470 + t449 * t466;
t305 = atan2(t319, t331);
t300 = sin(t305);
t301 = cos(t305);
t289 = t300 * t319 + t301 * t331;
t286 = 0.1e1 / t289;
t287 = 0.1e1 / t289 ^ 2;
t485 = 0.2e1 * t319;
t484 = 0.2e1 * t321;
t314 = t321 ^ 2;
t285 = t287 * t314 + 0.1e1;
t290 = qJD(4) * t322 + t296 * t375 - t377 * t398;
t474 = t290 * t287;
t313 = t319 ^ 2;
t304 = t313 * t329 + 0.1e1;
t302 = 0.1e1 / t304;
t411 = -t309 * t467 + t328 * t496;
t277 = t411 * t302;
t416 = -t300 * t331 + t301 * t319;
t272 = t277 * t416 + t300 * t496 + t301 * t309;
t288 = t286 * t287;
t478 = t272 * t288;
t479 = (-t314 * t478 + t321 * t474) / t285 ^ 2;
t469 = t328 * t491;
t476 = (-t313 * t469 + t467 * t496) / t304 ^ 2;
t475 = t287 * t321;
t473 = t300 * t321;
t472 = t301 * t321;
t471 = t306 * t335;
t468 = t319 * t328;
t465 = t335 * t342;
t464 = t342 * t375;
t459 = t372 * t375;
t458 = t372 * t377;
t456 = t374 * t376;
t452 = qJD(4) * t375;
t451 = qJD(4) * t377;
t450 = 0.2e1 * t479;
t448 = -0.2e1 * t476;
t447 = t288 * t484;
t444 = t328 * t476;
t443 = t287 * t473;
t442 = t287 * t472;
t439 = t306 * t466;
t430 = t272 * t447;
t428 = t469 * t485;
t410 = t320 * t328 - t332 * t467;
t346 = -t362 * t456 - t378 * t405;
t323 = t346 * t375 - t362 * t458;
t360 = t402 * t373;
t349 = t360 * t375 - t377 * t427;
t409 = -t323 * t328 - t349 * t467;
t348 = t392 * t456 - t391;
t407 = -t348 * t375 - t392 * t458;
t325 = t348 * t377 - t392 * t459;
t401 = -t374 * t433 - t454;
t400 = -t300 + (-t301 * t468 + t300) * t302;
t399 = -t352 * t455 + t378 * t414 + (qJD(3) * t426 - t353) * t376;
t395 = t405 * t376;
t347 = -t376 * t393 - t392 * t455;
t326 = -t376 * t419 + (qJD(2) * t401 - qJD(3) * t404) * t373;
t312 = t360 * t451 + ((qJD(3) * t401 + qJD(4) * t438) * t375 + (-t375 * t404 - t377 * t457) * qJD(2)) * t373;
t311 = -qJD(3) * t347 + t350 * t378 + t351 * t456;
t310 = -qJD(4) * t331 + t327 * t377 + t375 * t413;
t297 = (t374 * t395 - t461) * qJD(3) + t399;
t294 = (-t353 * t456 - t352 * t378 + (-t362 * t455 + t395) * qJD(3)) * t375 + t346 * t451 - t353 * t458 + t362 * t372 * t452;
t283 = 0.1e1 / t285;
t282 = t302 * t490;
t281 = t409 * t302;
t280 = t410 * t302;
t276 = t400 * t321;
t275 = (t300 * t486 + t301 * t356) * t375 + t416 * t282;
t274 = t281 * t416 - t300 * t323 + t301 * t349;
t273 = t280 * t416 + t300 * t320 + t301 * t332;
t271 = t409 * t448 + (t349 * t428 - t294 * t328 + (t309 * t323 - t312 * t319 - t349 * t496) * t329) * t302;
t269 = t410 * t448 + (t332 * t428 - t293 * t328 + (-t309 * t320 - t310 * t319 - t332 * t496) * t329) * t302;
t268 = t448 * t490 + (t408 * t451 + (t356 * t428 - t297 * t328 + (-t309 * t486 - t319 * t326 - t356 * t496) * t329) * t375) * t302;
t1 = [t444 * t484 + (-t290 * t328 + t321 * t491) * t302, t271, t268, t269, 0, 0; -0.2e1 * t319 * t286 * t479 + (t496 * t286 + (-t319 * t272 - t276 * t290) * t287) * t283 + (t276 * t287 * t450 + (0.2e1 * t276 * t478 - (t277 * t302 * t468 + t448) * t443 - (t444 * t485 - t277 + (t277 - t411) * t302) * t442 - t400 * t474) * t283) * t321 (t274 * t475 + t286 * t407) * t450 + ((qJD(4) * t325 + t311 * t375 + t351 * t458) * t286 + t274 * t430 + (t407 * t272 - t274 * t290 - (t271 * t319 + t281 * t496 + t312 + (-t281 * t331 - t323) * t277) * t472 - (-t271 * t331 - t281 * t309 - t294 + (-t281 * t319 - t349) * t277) * t473) * t287) * t283 (t275 * t475 + t286 * t464) * t450 + (-t275 * t474 + (-t295 * t375 - t342 * t451) * t286 + (t275 * t447 + t287 * t464) * t272 - (t356 * t451 + t268 * t319 + t282 * t496 + t326 * t375 + (-t282 * t331 + t375 * t486) * t277) * t442 - (t486 * t451 - t268 * t331 - t282 * t309 - t297 * t375 + (-t282 * t319 - t356 * t375) * t277) * t443) * t283 (t273 * t475 - t286 * t322) * t450 + (t273 * t430 + t291 * t286 + (-t322 * t272 - t273 * t290 - (t269 * t319 + t280 * t496 + t310 + (-t280 * t331 + t320) * t277) * t472 - (-t269 * t331 - t280 * t309 - t293 + (-t280 * t319 - t332) * t277) * t473) * t287) * t283, 0, 0; -t293 * t471 - ((t376 * t394 - t461) * qJD(3) + t399) * t439 + t489 * t320 - t488 * t486 (qJD(4) * t407 + t311 * t377 - t351 * t459) * t471 - (qJD(3) * t348 + t350 * t376 - t351 * t455) * t439 + t489 * t325 + t488 * t347 (t343 * t466 + t377 * t465) * t449 + (t452 * t465 + (-t291 * t343 - t296 * t322) * t336 + (t343 * t446 + (t336 * t342 - t335) * t377) * t295) * t306, t445 * t484 + (t295 * t321 * t336 - t290 * t335) * t306, 0, 0;];
JaD_rot  = t1;
