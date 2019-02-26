% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR6
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRR6_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR6_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_jacobiaD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:07:07
% EndTime: 2019-02-26 20:07:10
% DurationCPUTime: 2.86s
% Computational Cost: add. (19059->203), mult. (44088->383), div. (816->12), fcn. (56595->17), ass. (0->161)
t404 = sin(qJ(3));
t407 = cos(qJ(3));
t400 = sin(pkin(7));
t402 = cos(pkin(7));
t405 = sin(qJ(2));
t408 = cos(qJ(2));
t483 = cos(pkin(12));
t484 = cos(pkin(6));
t437 = t484 * t483;
t482 = sin(pkin(12));
t423 = t482 * t405 - t408 * t437;
t401 = sin(pkin(6));
t447 = t401 * t483;
t415 = -t400 * t447 - t423 * t402;
t422 = -t405 * t437 - t482 * t408;
t372 = t404 * t422 + t415 * t407;
t387 = t423 * qJD(2);
t388 = t422 * qJD(2);
t463 = t402 * t404;
t349 = t372 * qJD(3) - t387 * t407 + t388 * t463;
t373 = t415 * t404 - t407 * t422;
t399 = pkin(13) + qJ(5);
t397 = sin(t399);
t398 = cos(t399);
t416 = t423 * t400 - t402 * t447;
t356 = t373 * t398 + t416 * t397;
t465 = t398 * t400;
t323 = t356 * qJD(5) + t349 * t397 + t388 * t465;
t354 = t373 * t397 - t416 * t398;
t352 = t354 ^ 2;
t459 = t405 * t407;
t460 = t404 * t408;
t426 = t402 * t460 + t459;
t448 = t400 * t484;
t383 = t426 * t401 + t404 * t448;
t464 = t400 * t408;
t391 = -t401 * t464 + t484 * t402;
t367 = t383 * t397 - t391 * t398;
t365 = 0.1e1 / t367 ^ 2;
t339 = t352 * t365 + 0.1e1;
t337 = 0.1e1 / t339;
t368 = t383 * t398 + t391 * t397;
t458 = t407 * t408;
t461 = t404 * t405;
t425 = -t402 * t461 + t458;
t428 = t402 * t458 - t461;
t439 = qJD(3) * t448;
t370 = t407 * t439 + (t425 * qJD(2) + t428 * qJD(3)) * t401;
t449 = t400 * t401 * t405;
t441 = qJD(2) * t449;
t341 = t368 * qJD(5) + t370 * t397 - t398 * t441;
t364 = 0.1e1 / t367;
t472 = t354 * t365;
t306 = (-t323 * t364 + t341 * t472) * t337;
t340 = atan2(-t354, t367);
t333 = sin(t340);
t334 = cos(t340);
t435 = -t333 * t367 - t334 * t354;
t301 = t435 * t306 - t333 * t323 + t334 * t341;
t319 = -t333 * t354 + t334 * t367;
t316 = 0.1e1 / t319;
t317 = 0.1e1 / t319 ^ 2;
t487 = t301 * t316 * t317;
t436 = t484 * t482;
t420 = t483 * t405 + t408 * t436;
t446 = t401 * t482;
t440 = t400 * t446;
t417 = -t420 * t402 + t440;
t421 = t405 * t436 - t483 * t408;
t375 = t417 * t404 - t407 * t421;
t418 = t420 * t400 + t402 * t446;
t357 = t375 * t397 - t418 * t398;
t444 = 0.2e1 * t357 * t487;
t382 = t428 * t401 + t407 * t448;
t430 = -t364 * t372 + t382 * t472;
t486 = t397 * t430;
t473 = t341 * t364 * t365;
t485 = -0.2e1 * (t323 * t472 - t352 * t473) / t339 ^ 2;
t358 = t375 * t398 + t418 * t397;
t406 = cos(qJ(6));
t419 = t420 * t407;
t467 = t421 * t404;
t374 = t402 * t419 - t407 * t440 - t467;
t403 = sin(qJ(6));
t470 = t374 * t403;
t336 = t358 * t406 + t470;
t330 = 0.1e1 / t336;
t331 = 0.1e1 / t336 ^ 2;
t389 = t420 * qJD(2);
t390 = t421 * qJD(2);
t351 = t390 * t463 - t389 * t407 + (t417 * t407 + t467) * qJD(3);
t466 = t397 * t400;
t326 = -t357 * qJD(5) + t351 * t398 - t390 * t466;
t462 = t402 * t407;
t350 = t375 * qJD(3) - t389 * t404 - t390 * t462;
t469 = t374 * t406;
t335 = t358 * t403 - t469;
t455 = qJD(6) * t335;
t315 = t326 * t406 + t350 * t403 - t455;
t481 = t315 * t330 * t331;
t480 = t317 * t357;
t314 = t336 * qJD(6) + t326 * t403 - t350 * t406;
t329 = t335 ^ 2;
t322 = t329 * t331 + 0.1e1;
t477 = t331 * t335;
t479 = 0.1e1 / t322 ^ 2 * (t314 * t477 - t329 * t481);
t478 = t330 * t403;
t476 = t333 * t357;
t475 = t334 * t357;
t474 = t335 * t406;
t471 = t374 * t397;
t457 = qJD(5) * t398;
t456 = qJD(5) * t400;
t353 = t357 ^ 2;
t313 = t317 * t353 + 0.1e1;
t325 = t358 * qJD(5) + t351 * t397 + t390 * t465;
t454 = 0.2e1 * (t325 * t480 - t353 * t487) / t313 ^ 2;
t452 = -0.2e1 * t479;
t451 = 0.2e1 * t479;
t450 = t335 * t481;
t443 = 0.2e1 * t450;
t442 = -0.2e1 * t354 * t473;
t438 = qJD(6) * t374 * t398 + t351;
t378 = t421 * t463 - t419;
t363 = t378 * t398 - t421 * t466;
t377 = -t420 * t404 - t421 * t462;
t344 = t363 * t406 + t377 * t403;
t343 = t363 * t403 - t377 * t406;
t433 = t331 * t474 - t478;
t432 = -t356 * t364 + t368 * t472;
t376 = -t423 * t407 + t422 * t463;
t361 = t376 * t397 + t422 * t465;
t386 = t425 * t401;
t379 = t386 * t397 - t398 * t449;
t431 = -t361 * t364 + t379 * t472;
t429 = -t378 * t397 - t421 * t465;
t427 = -t402 * t459 - t460;
t424 = qJD(5) * t471 + qJD(6) * t375 - t350 * t398;
t369 = -t404 * t439 + (t427 * qJD(2) - t426 * qJD(3)) * t401;
t360 = -t377 * qJD(3) + t389 * t463 + t390 * t407;
t359 = t378 * qJD(3) - t389 * t462 + t390 * t404;
t348 = -t373 * qJD(3) + t387 * t404 + t388 * t462;
t347 = t386 * t457 + ((t427 * qJD(3) + t405 * t456) * t397 + (-t426 * t397 - t398 * t464) * qJD(2)) * t401;
t346 = t375 * t403 - t398 * t469;
t345 = -t375 * t406 - t398 * t470;
t342 = -t367 * qJD(5) + t370 * t398 + t397 * t441;
t328 = t429 * qJD(5) + t360 * t398 - t389 * t466;
t327 = t376 * t457 + t387 * t465 + (t387 * t463 + t388 * t407 + (t423 * t404 + t422 * t462) * qJD(3) - t422 * t456) * t397;
t324 = -t354 * qJD(5) + t349 * t398 - t388 * t466;
t320 = 0.1e1 / t322;
t311 = 0.1e1 / t313;
t310 = t337 * t486;
t309 = t431 * t337;
t308 = t432 * t337;
t304 = (-t333 * t372 + t334 * t382) * t397 + t435 * t310;
t303 = t435 * t309 - t333 * t361 + t334 * t379;
t302 = t435 * t308 - t333 * t356 + t334 * t368;
t300 = t431 * t485 + (t379 * t442 - t327 * t364 + (t323 * t379 + t341 * t361 + t347 * t354) * t365) * t337;
t298 = t432 * t485 + (t368 * t442 - t324 * t364 + (t323 * t368 + t341 * t356 + t342 * t354) * t365) * t337;
t297 = t485 * t486 + (t430 * t457 + (t382 * t442 - t348 * t364 + (t323 * t382 + t341 * t372 + t354 * t369) * t365) * t397) * t337;
t1 = [0, t300, t297, 0, t298, 0; 0 (t303 * t480 + t316 * t429) * t454 + ((t363 * qJD(5) + t360 * t397 + t389 * t465) * t316 + t303 * t444 + (t429 * t301 - t303 * t325 - (-t300 * t354 - t309 * t323 + t347 + (-t309 * t367 - t361) * t306) * t475 - (-t300 * t367 - t309 * t341 - t327 + (t309 * t354 - t379) * t306) * t476) * t317) * t311 (t304 * t480 + t316 * t471) * t454 + ((-t350 * t397 - t374 * t457) * t316 + t304 * t444 + (-t304 * t325 + t471 * t301 - (t382 * t457 - t297 * t354 - t310 * t323 + t369 * t397 + (-t310 * t367 - t372 * t397) * t306) * t475 - (-t372 * t457 - t297 * t367 - t310 * t341 - t348 * t397 + (t310 * t354 - t382 * t397) * t306) * t476) * t317) * t311, 0 (t302 * t480 - t316 * t358) * t454 + (t302 * t444 + t326 * t316 + (-t358 * t301 - t302 * t325 - (-t298 * t354 - t308 * t323 + t342 + (-t308 * t367 - t356) * t306) * t475 - (-t298 * t367 - t308 * t341 - t324 + (t308 * t354 - t368) * t306) * t476) * t317) * t311, 0; 0 (-t330 * t343 + t344 * t477) * t451 + ((t344 * qJD(6) + t328 * t403 - t359 * t406) * t330 + t344 * t443 + (-t343 * t315 - (-t343 * qJD(6) + t328 * t406 + t359 * t403) * t335 - t344 * t314) * t331) * t320 (-t330 * t345 + t346 * t477) * t451 + (t346 * t443 - t438 * t330 * t406 + t424 * t478 + (-t438 * t335 * t403 - t346 * t314 - t345 * t315 - t424 * t474) * t331) * t320, 0, t433 * t357 * t452 + (t433 * t325 + ((-qJD(6) * t330 - 0.2e1 * t450) * t406 + (t314 * t406 + (t315 - t455) * t403) * t331) * t357) * t320, t452 + 0.2e1 * (t314 * t331 * t320 + (-t320 * t481 - t331 * t479) * t335) * t335;];
JaD_rot  = t1;
