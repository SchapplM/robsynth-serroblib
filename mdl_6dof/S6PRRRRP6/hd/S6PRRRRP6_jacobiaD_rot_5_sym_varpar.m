% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRRP6_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP6_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:18:03
% EndTime: 2019-02-26 20:18:05
% DurationCPUTime: 2.72s
% Computational Cost: add. (14691->202), mult. (44088->385), div. (816->12), fcn. (56595->17), ass. (0->160)
t391 = sin(qJ(3));
t395 = cos(qJ(3));
t386 = sin(pkin(7));
t388 = cos(pkin(7));
t392 = sin(qJ(2));
t396 = cos(qJ(2));
t471 = cos(pkin(12));
t472 = cos(pkin(6));
t425 = t472 * t471;
t470 = sin(pkin(12));
t411 = t470 * t392 - t396 * t425;
t387 = sin(pkin(6));
t435 = t387 * t471;
t403 = -t386 * t435 - t411 * t388;
t410 = -t392 * t425 - t470 * t396;
t359 = t391 * t410 + t403 * t395;
t376 = t411 * qJD(2);
t377 = t410 * qJD(2);
t451 = t388 * t391;
t337 = t359 * qJD(3) - t376 * t395 + t377 * t451;
t360 = t403 * t391 - t395 * t410;
t390 = sin(qJ(4));
t394 = cos(qJ(4));
t404 = t411 * t386 - t388 * t435;
t347 = t360 * t394 + t404 * t390;
t453 = t386 * t394;
t312 = t347 * qJD(4) + t337 * t390 + t377 * t453;
t345 = t360 * t390 - t404 * t394;
t343 = t345 ^ 2;
t447 = t392 * t395;
t448 = t391 * t396;
t414 = t388 * t448 + t447;
t436 = t386 * t472;
t372 = t414 * t387 + t391 * t436;
t452 = t386 * t396;
t380 = -t387 * t452 + t472 * t388;
t363 = t372 * t390 - t380 * t394;
t357 = 0.1e1 / t363 ^ 2;
t328 = t343 * t357 + 0.1e1;
t326 = 0.1e1 / t328;
t446 = t395 * t396;
t449 = t391 * t392;
t413 = -t388 * t449 + t446;
t416 = t388 * t446 - t449;
t427 = qJD(3) * t436;
t354 = t395 * t427 + (t413 * qJD(2) + t416 * qJD(3)) * t387;
t364 = t372 * t394 + t380 * t390;
t454 = t386 * t392;
t437 = t387 * t454;
t429 = qJD(2) * t437;
t330 = t364 * qJD(4) + t354 * t390 - t394 * t429;
t356 = 0.1e1 / t363;
t460 = t345 * t357;
t295 = (-t312 * t356 + t330 * t460) * t326;
t329 = atan2(-t345, t363);
t324 = sin(t329);
t325 = cos(t329);
t423 = -t324 * t363 - t325 * t345;
t290 = t423 * t295 - t324 * t312 + t325 * t330;
t308 = -t324 * t345 + t325 * t363;
t305 = 0.1e1 / t308;
t306 = 0.1e1 / t308 ^ 2;
t475 = t290 * t305 * t306;
t424 = t472 * t470;
t408 = t471 * t392 + t396 * t424;
t434 = t387 * t470;
t428 = t386 * t434;
t405 = -t408 * t388 + t428;
t409 = t392 * t424 - t471 * t396;
t362 = t405 * t391 - t395 * t409;
t406 = t408 * t386 + t388 * t434;
t348 = t362 * t390 - t406 * t394;
t432 = 0.2e1 * t348 * t475;
t371 = t416 * t387 + t395 * t436;
t418 = -t356 * t359 + t371 * t460;
t474 = t390 * t418;
t461 = t330 * t356 * t357;
t473 = -0.2e1 * (t312 * t460 - t343 * t461) / t328 ^ 2;
t349 = t362 * t394 + t406 * t390;
t407 = t408 * t395;
t456 = t409 * t391;
t361 = t388 * t407 - t395 * t428 - t456;
t389 = sin(qJ(5));
t393 = cos(qJ(5));
t323 = t349 * t393 + t361 * t389;
t319 = 0.1e1 / t323;
t320 = 0.1e1 / t323 ^ 2;
t378 = t408 * qJD(2);
t379 = t409 * qJD(2);
t339 = t379 * t451 - t378 * t395 + (t405 * t395 + t456) * qJD(3);
t455 = t386 * t390;
t315 = -t348 * qJD(4) + t339 * t394 - t379 * t455;
t450 = t388 * t395;
t338 = t362 * qJD(3) - t378 * t391 - t379 * t450;
t322 = t349 * t389 - t361 * t393;
t443 = qJD(5) * t322;
t304 = t315 * t393 + t338 * t389 - t443;
t469 = t304 * t319 * t320;
t468 = t306 * t348;
t303 = t323 * qJD(5) + t315 * t389 - t338 * t393;
t318 = t322 ^ 2;
t311 = t318 * t320 + 0.1e1;
t465 = t320 * t322;
t467 = 0.1e1 / t311 ^ 2 * (t303 * t465 - t318 * t469);
t466 = t319 * t389;
t464 = t322 * t393;
t463 = t324 * t348;
t462 = t325 * t348;
t459 = t361 * t390;
t458 = t361 * t394;
t445 = qJD(4) * t390;
t444 = qJD(4) * t394;
t344 = t348 ^ 2;
t302 = t344 * t306 + 0.1e1;
t314 = t349 * qJD(4) + t339 * t390 + t379 * t453;
t442 = 0.2e1 * (t314 * t468 - t344 * t475) / t302 ^ 2;
t440 = -0.2e1 * t467;
t439 = 0.2e1 * t467;
t438 = t322 * t469;
t431 = 0.2e1 * t438;
t430 = -0.2e1 * t345 * t461;
t426 = qJD(5) * t458 + t339;
t367 = t409 * t451 - t407;
t352 = t367 * t394 - t409 * t455;
t366 = -t408 * t391 - t409 * t450;
t335 = t352 * t393 + t366 * t389;
t334 = t352 * t389 - t366 * t393;
t421 = t320 * t464 - t466;
t420 = -t347 * t356 + t364 * t460;
t365 = -t411 * t395 + t410 * t451;
t350 = t365 * t390 + t410 * t453;
t375 = t413 * t387;
t368 = t375 * t390 - t394 * t437;
t419 = -t350 * t356 + t368 * t460;
t417 = -t367 * t390 - t409 * t453;
t415 = -t388 * t447 - t448;
t412 = qJD(5) * t362 - t338 * t394 + t361 * t445;
t353 = -t391 * t427 + (t415 * qJD(2) - t414 * qJD(3)) * t387;
t342 = -t366 * qJD(3) + t378 * t451 + t379 * t395;
t341 = t367 * qJD(3) - t378 * t450 + t379 * t391;
t340 = t375 * t444 + ((t415 * qJD(3) + qJD(4) * t454) * t390 + (-t414 * t390 - t394 * t452) * qJD(2)) * t387;
t336 = -t360 * qJD(3) + t376 * t391 + t377 * t450;
t333 = t362 * t389 - t393 * t458;
t332 = -t362 * t393 - t389 * t458;
t331 = -t363 * qJD(4) + t354 * t394 + t390 * t429;
t317 = t417 * qJD(4) + t342 * t394 - t378 * t455;
t316 = (t376 * t451 + t377 * t395 + (t411 * t391 + t410 * t450) * qJD(3)) * t390 + t365 * t444 + t376 * t453 - t410 * t386 * t445;
t313 = -t345 * qJD(4) + t337 * t394 - t377 * t455;
t309 = 0.1e1 / t311;
t300 = 0.1e1 / t302;
t299 = t326 * t474;
t298 = t419 * t326;
t297 = t420 * t326;
t293 = (-t324 * t359 + t325 * t371) * t390 + t423 * t299;
t292 = t423 * t298 - t324 * t350 + t325 * t368;
t291 = t423 * t297 - t324 * t347 + t325 * t364;
t289 = t419 * t473 + (t368 * t430 - t316 * t356 + (t312 * t368 + t330 * t350 + t340 * t345) * t357) * t326;
t287 = t420 * t473 + (t364 * t430 - t313 * t356 + (t312 * t364 + t330 * t347 + t331 * t345) * t357) * t326;
t286 = t473 * t474 + (t418 * t444 + (t371 * t430 - t336 * t356 + (t312 * t371 + t330 * t359 + t345 * t353) * t357) * t390) * t326;
t1 = [0, t289, t286, t287, 0, 0; 0 (t292 * t468 + t305 * t417) * t442 + ((t352 * qJD(4) + t342 * t390 + t378 * t453) * t305 + t292 * t432 + (t417 * t290 - t292 * t314 - (-t289 * t345 - t298 * t312 + t340 + (-t298 * t363 - t350) * t295) * t462 - (-t289 * t363 - t298 * t330 - t316 + (t298 * t345 - t368) * t295) * t463) * t306) * t300 (t293 * t468 + t305 * t459) * t442 + ((-t338 * t390 - t361 * t444) * t305 + t293 * t432 + (-t293 * t314 + t459 * t290 - (t371 * t444 - t286 * t345 - t299 * t312 + t353 * t390 + (-t299 * t363 - t359 * t390) * t295) * t462 - (-t359 * t444 - t286 * t363 - t299 * t330 - t336 * t390 + (t299 * t345 - t371 * t390) * t295) * t463) * t306) * t300 (t291 * t468 - t305 * t349) * t442 + (t291 * t432 + t315 * t305 + (-t349 * t290 - t291 * t314 - (-t287 * t345 - t297 * t312 + t331 + (-t297 * t363 - t347) * t295) * t462 - (-t287 * t363 - t297 * t330 - t313 + (t297 * t345 - t364) * t295) * t463) * t306) * t300, 0, 0; 0 (-t319 * t334 + t335 * t465) * t439 + ((t335 * qJD(5) + t317 * t389 - t341 * t393) * t319 + t335 * t431 + (-t334 * t304 - (-t334 * qJD(5) + t317 * t393 + t341 * t389) * t322 - t335 * t303) * t320) * t309 (-t319 * t332 + t333 * t465) * t439 + (t333 * t431 - t426 * t319 * t393 + t412 * t466 + (-t426 * t322 * t389 - t333 * t303 - t332 * t304 - t412 * t464) * t320) * t309, t421 * t348 * t440 + (t421 * t314 + ((-qJD(5) * t319 - 0.2e1 * t438) * t393 + (t303 * t393 + (t304 - t443) * t389) * t320) * t348) * t309, t440 + 0.2e1 * (t303 * t320 * t309 + (-t309 * t469 - t320 * t467) * t322) * t322, 0;];
JaD_rot  = t1;
