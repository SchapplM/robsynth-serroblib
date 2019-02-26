% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRPR5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_jacobiaD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:12:55
% EndTime: 2019-02-26 20:12:58
% DurationCPUTime: 2.92s
% Computational Cost: add. (19059->203), mult. (44088->383), div. (816->12), fcn. (56595->17), ass. (0->161)
t405 = sin(qJ(3));
t408 = cos(qJ(3));
t401 = sin(pkin(7));
t403 = cos(pkin(7));
t406 = sin(qJ(2));
t409 = cos(qJ(2));
t484 = cos(pkin(12));
t485 = cos(pkin(6));
t438 = t485 * t484;
t483 = sin(pkin(12));
t424 = t483 * t406 - t409 * t438;
t402 = sin(pkin(6));
t448 = t402 * t484;
t416 = -t401 * t448 - t424 * t403;
t423 = -t406 * t438 - t483 * t409;
t373 = t405 * t423 + t416 * t408;
t388 = t424 * qJD(2);
t389 = t423 * qJD(2);
t464 = t403 * t405;
t350 = t373 * qJD(3) - t388 * t408 + t389 * t464;
t374 = t416 * t405 - t408 * t423;
t400 = qJ(4) + pkin(13);
t398 = sin(t400);
t399 = cos(t400);
t417 = t424 * t401 - t403 * t448;
t357 = t374 * t399 + t417 * t398;
t466 = t399 * t401;
t324 = t357 * qJD(4) + t350 * t398 + t389 * t466;
t355 = t374 * t398 - t417 * t399;
t353 = t355 ^ 2;
t460 = t406 * t408;
t461 = t405 * t409;
t427 = t403 * t461 + t460;
t449 = t401 * t485;
t384 = t427 * t402 + t405 * t449;
t465 = t401 * t409;
t392 = -t402 * t465 + t485 * t403;
t368 = t384 * t398 - t392 * t399;
t366 = 0.1e1 / t368 ^ 2;
t340 = t353 * t366 + 0.1e1;
t338 = 0.1e1 / t340;
t369 = t384 * t399 + t392 * t398;
t459 = t408 * t409;
t462 = t405 * t406;
t426 = -t403 * t462 + t459;
t429 = t403 * t459 - t462;
t440 = qJD(3) * t449;
t371 = t408 * t440 + (t426 * qJD(2) + t429 * qJD(3)) * t402;
t450 = t401 * t402 * t406;
t442 = qJD(2) * t450;
t342 = t369 * qJD(4) + t371 * t398 - t399 * t442;
t365 = 0.1e1 / t368;
t473 = t355 * t366;
t307 = (-t324 * t365 + t342 * t473) * t338;
t341 = atan2(-t355, t368);
t334 = sin(t341);
t335 = cos(t341);
t436 = -t334 * t368 - t335 * t355;
t302 = t436 * t307 - t324 * t334 + t335 * t342;
t320 = -t334 * t355 + t335 * t368;
t317 = 0.1e1 / t320;
t318 = 0.1e1 / t320 ^ 2;
t488 = t302 * t317 * t318;
t437 = t485 * t483;
t421 = t484 * t406 + t409 * t437;
t447 = t402 * t483;
t441 = t401 * t447;
t418 = -t421 * t403 + t441;
t422 = t406 * t437 - t484 * t409;
t376 = t418 * t405 - t408 * t422;
t419 = t421 * t401 + t403 * t447;
t358 = t376 * t398 - t419 * t399;
t445 = 0.2e1 * t358 * t488;
t383 = t429 * t402 + t408 * t449;
t431 = -t365 * t373 + t383 * t473;
t487 = t398 * t431;
t474 = t342 * t365 * t366;
t486 = -0.2e1 * (t324 * t473 - t353 * t474) / t340 ^ 2;
t359 = t376 * t399 + t419 * t398;
t407 = cos(qJ(6));
t420 = t421 * t408;
t468 = t422 * t405;
t375 = t403 * t420 - t408 * t441 - t468;
t404 = sin(qJ(6));
t471 = t375 * t404;
t337 = t359 * t407 + t471;
t331 = 0.1e1 / t337;
t332 = 0.1e1 / t337 ^ 2;
t390 = t421 * qJD(2);
t391 = t422 * qJD(2);
t352 = t391 * t464 - t390 * t408 + (t418 * t408 + t468) * qJD(3);
t467 = t398 * t401;
t327 = -t358 * qJD(4) + t352 * t399 - t391 * t467;
t463 = t403 * t408;
t351 = t376 * qJD(3) - t390 * t405 - t391 * t463;
t315 = t337 * qJD(6) + t327 * t404 - t351 * t407;
t470 = t375 * t407;
t336 = t359 * t404 - t470;
t330 = t336 ^ 2;
t323 = t330 * t332 + 0.1e1;
t478 = t332 * t336;
t456 = qJD(6) * t336;
t316 = t327 * t407 + t351 * t404 - t456;
t481 = t316 * t331 * t332;
t482 = (t315 * t478 - t330 * t481) / t323 ^ 2;
t480 = t318 * t358;
t479 = t331 * t404;
t477 = t334 * t358;
t476 = t335 * t358;
t475 = t336 * t407;
t472 = t375 * t398;
t458 = qJD(4) * t399;
t457 = qJD(4) * t401;
t354 = t358 ^ 2;
t314 = t318 * t354 + 0.1e1;
t326 = t359 * qJD(4) + t352 * t398 + t391 * t466;
t455 = 0.2e1 * (t326 * t480 - t354 * t488) / t314 ^ 2;
t454 = -0.2e1 * t482;
t453 = 0.2e1 * t482;
t451 = t336 * t481;
t444 = 0.2e1 * t451;
t443 = -0.2e1 * t355 * t474;
t439 = qJD(6) * t375 * t399 + t352;
t379 = t422 * t464 - t420;
t364 = t379 * t399 - t422 * t467;
t378 = -t421 * t405 - t422 * t463;
t345 = t364 * t407 + t378 * t404;
t344 = t364 * t404 - t378 * t407;
t434 = t332 * t475 - t479;
t433 = -t357 * t365 + t369 * t473;
t377 = -t424 * t408 + t423 * t464;
t362 = t377 * t398 + t423 * t466;
t387 = t426 * t402;
t380 = t387 * t398 - t399 * t450;
t432 = -t362 * t365 + t380 * t473;
t430 = -t379 * t398 - t422 * t466;
t428 = -t403 * t460 - t461;
t425 = qJD(4) * t472 + qJD(6) * t376 - t351 * t399;
t370 = -t405 * t440 + (t428 * qJD(2) - t427 * qJD(3)) * t402;
t361 = -t378 * qJD(3) + t390 * t464 + t391 * t408;
t360 = t379 * qJD(3) - t390 * t463 + t391 * t405;
t349 = -qJD(3) * t374 + t388 * t405 + t389 * t463;
t348 = t387 * t458 + ((t428 * qJD(3) + t406 * t457) * t398 + (-t427 * t398 - t399 * t465) * qJD(2)) * t402;
t347 = t376 * t404 - t399 * t470;
t346 = -t376 * t407 - t399 * t471;
t343 = -t368 * qJD(4) + t371 * t399 + t398 * t442;
t329 = t430 * qJD(4) + t361 * t399 - t390 * t467;
t328 = t377 * t458 + t388 * t466 + (t388 * t464 + t389 * t408 + (t424 * t405 + t423 * t463) * qJD(3) - t423 * t457) * t398;
t325 = -t355 * qJD(4) + t350 * t399 - t389 * t467;
t321 = 0.1e1 / t323;
t312 = 0.1e1 / t314;
t311 = t338 * t487;
t310 = t432 * t338;
t309 = t433 * t338;
t305 = (-t334 * t373 + t335 * t383) * t398 + t436 * t311;
t304 = t436 * t310 - t334 * t362 + t335 * t380;
t303 = t436 * t309 - t334 * t357 + t335 * t369;
t301 = t432 * t486 + (t380 * t443 - t328 * t365 + (t324 * t380 + t342 * t362 + t348 * t355) * t366) * t338;
t299 = t433 * t486 + (t369 * t443 - t325 * t365 + (t324 * t369 + t342 * t357 + t343 * t355) * t366) * t338;
t298 = t486 * t487 + (t431 * t458 + (t383 * t443 - t349 * t365 + (t324 * t383 + t342 * t373 + t355 * t370) * t366) * t398) * t338;
t1 = [0, t301, t298, t299, 0, 0; 0 (t304 * t480 + t317 * t430) * t455 + ((t364 * qJD(4) + t361 * t398 + t390 * t466) * t317 + t304 * t445 + (t430 * t302 - t304 * t326 - (-t301 * t355 - t310 * t324 + t348 + (-t310 * t368 - t362) * t307) * t476 - (-t301 * t368 - t310 * t342 - t328 + (t310 * t355 - t380) * t307) * t477) * t318) * t312 (t305 * t480 + t317 * t472) * t455 + ((-t351 * t398 - t375 * t458) * t317 + t305 * t445 + (-t305 * t326 + t472 * t302 - (t383 * t458 - t298 * t355 - t311 * t324 + t370 * t398 + (-t311 * t368 - t373 * t398) * t307) * t476 - (-t373 * t458 - t298 * t368 - t311 * t342 - t349 * t398 + (t311 * t355 - t383 * t398) * t307) * t477) * t318) * t312 (t303 * t480 - t317 * t359) * t455 + (t303 * t445 + t327 * t317 + (-t359 * t302 - t303 * t326 - (-t299 * t355 - t309 * t324 + t343 + (-t309 * t368 - t357) * t307) * t476 - (-t299 * t368 - t309 * t342 - t325 + (t309 * t355 - t369) * t307) * t477) * t318) * t312, 0, 0; 0 (-t331 * t344 + t345 * t478) * t453 + ((t345 * qJD(6) + t329 * t404 - t360 * t407) * t331 + t345 * t444 + (-t344 * t316 - (-t344 * qJD(6) + t329 * t407 + t360 * t404) * t336 - t345 * t315) * t332) * t321 (-t331 * t346 + t347 * t478) * t453 + (t347 * t444 - t439 * t331 * t407 + t425 * t479 + (-t439 * t336 * t404 - t347 * t315 - t346 * t316 - t425 * t475) * t332) * t321, t434 * t358 * t454 + (t434 * t326 + ((-qJD(6) * t331 - 0.2e1 * t451) * t407 + (t315 * t407 + (t316 - t456) * t404) * t332) * t358) * t321, 0, t454 + 0.2e1 * (t315 * t332 * t321 + (-t321 * t481 - t332 * t482) * t336) * t336;];
JaD_rot  = t1;
