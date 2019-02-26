% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR7
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
% Datum: 2019-02-26 20:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRPR7_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR7_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_jacobiaD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:14:05
% EndTime: 2019-02-26 20:14:07
% DurationCPUTime: 2.73s
% Computational Cost: add. (15112->203), mult. (44088->385), div. (816->12), fcn. (56595->17), ass. (0->162)
t402 = sin(qJ(3));
t405 = cos(qJ(3));
t398 = sin(pkin(7));
t400 = cos(pkin(7));
t403 = sin(qJ(2));
t406 = cos(qJ(2));
t482 = cos(pkin(12));
t483 = cos(pkin(6));
t435 = t483 * t482;
t481 = sin(pkin(12));
t421 = t481 * t403 - t406 * t435;
t399 = sin(pkin(6));
t445 = t399 * t482;
t413 = -t398 * t445 - t421 * t400;
t420 = -t403 * t435 - t481 * t406;
t368 = t402 * t420 + t413 * t405;
t385 = t421 * qJD(2);
t386 = t420 * qJD(2);
t461 = t400 * t402;
t346 = t368 * qJD(3) - t385 * t405 + t386 * t461;
t369 = t413 * t402 - t405 * t420;
t401 = sin(qJ(4));
t404 = cos(qJ(4));
t414 = t421 * t398 - t400 * t445;
t356 = t369 * t404 + t414 * t401;
t463 = t398 * t404;
t321 = t356 * qJD(4) + t346 * t401 + t386 * t463;
t354 = t369 * t401 - t414 * t404;
t352 = t354 ^ 2;
t457 = t403 * t405;
t458 = t402 * t406;
t424 = t400 * t458 + t457;
t446 = t398 * t483;
t381 = t424 * t399 + t402 * t446;
t462 = t398 * t406;
t389 = -t399 * t462 + t483 * t400;
t372 = t381 * t401 - t389 * t404;
t366 = 0.1e1 / t372 ^ 2;
t337 = t352 * t366 + 0.1e1;
t335 = 0.1e1 / t337;
t456 = t405 * t406;
t459 = t402 * t403;
t423 = -t400 * t459 + t456;
t426 = t400 * t456 - t459;
t437 = qJD(3) * t446;
t364 = t405 * t437 + (t423 * qJD(2) + t426 * qJD(3)) * t399;
t373 = t381 * t404 + t389 * t401;
t464 = t398 * t403;
t447 = t399 * t464;
t439 = qJD(2) * t447;
t339 = t373 * qJD(4) + t364 * t401 - t404 * t439;
t365 = 0.1e1 / t372;
t472 = t354 * t366;
t304 = (-t321 * t365 + t339 * t472) * t335;
t338 = atan2(-t354, t372);
t333 = sin(t338);
t334 = cos(t338);
t433 = -t333 * t372 - t334 * t354;
t299 = t433 * t304 - t333 * t321 + t334 * t339;
t317 = -t333 * t354 + t334 * t372;
t314 = 0.1e1 / t317;
t315 = 0.1e1 / t317 ^ 2;
t486 = t299 * t314 * t315;
t434 = t483 * t481;
t418 = t482 * t403 + t406 * t434;
t444 = t399 * t481;
t438 = t398 * t444;
t415 = -t418 * t400 + t438;
t419 = t403 * t434 - t482 * t406;
t371 = t415 * t402 - t405 * t419;
t416 = t418 * t398 + t400 * t444;
t357 = t371 * t401 - t416 * t404;
t442 = 0.2e1 * t357 * t486;
t380 = t426 * t399 + t405 * t446;
t428 = -t365 * t368 + t380 * t472;
t485 = t401 * t428;
t473 = t339 * t365 * t366;
t484 = -0.2e1 * (t321 * t472 - t352 * t473) / t337 ^ 2;
t358 = t371 * t404 + t416 * t401;
t417 = t418 * t405;
t468 = t419 * t402;
t370 = t400 * t417 - t405 * t438 - t468;
t397 = pkin(13) + qJ(6);
t395 = sin(t397);
t396 = cos(t397);
t332 = t358 * t396 + t370 * t395;
t328 = 0.1e1 / t332;
t329 = 0.1e1 / t332 ^ 2;
t387 = t418 * qJD(2);
t388 = t419 * qJD(2);
t348 = t388 * t461 - t387 * t405 + (t415 * t405 + t468) * qJD(3);
t465 = t398 * t401;
t324 = -t357 * qJD(4) + t348 * t404 - t388 * t465;
t460 = t400 * t405;
t347 = t371 * qJD(3) - t387 * t402 - t388 * t460;
t312 = t332 * qJD(6) + t324 * t395 - t347 * t396;
t331 = t358 * t395 - t370 * t396;
t327 = t331 ^ 2;
t320 = t327 * t329 + 0.1e1;
t476 = t329 * t331;
t453 = qJD(6) * t331;
t313 = t324 * t396 + t347 * t395 - t453;
t479 = t313 * t328 * t329;
t480 = (t312 * t476 - t327 * t479) / t320 ^ 2;
t478 = t315 * t357;
t323 = t358 * qJD(4) + t348 * t401 + t388 * t463;
t477 = t323 * t315;
t475 = t333 * t357;
t474 = t334 * t357;
t471 = t370 * t401;
t470 = t370 * t404;
t467 = t395 * t328;
t466 = t396 * t331;
t455 = qJD(4) * t401;
t454 = qJD(4) * t404;
t353 = t357 ^ 2;
t311 = t353 * t315 + 0.1e1;
t452 = 0.2e1 * (-t353 * t486 + t357 * t477) / t311 ^ 2;
t451 = -0.2e1 * t480;
t450 = 0.2e1 * t480;
t448 = t331 * t479;
t441 = 0.2e1 * t448;
t440 = -0.2e1 * t354 * t473;
t436 = qJD(6) * t470 + t348;
t376 = t419 * t461 - t417;
t361 = t376 * t404 - t419 * t465;
t375 = -t418 * t402 - t419 * t460;
t344 = t361 * t396 + t375 * t395;
t343 = t361 * t395 - t375 * t396;
t431 = t329 * t466 - t467;
t430 = -t356 * t365 + t373 * t472;
t374 = -t421 * t405 + t420 * t461;
t359 = t374 * t401 + t420 * t463;
t384 = t423 * t399;
t377 = t384 * t401 - t404 * t447;
t429 = -t359 * t365 + t377 * t472;
t427 = -t376 * t401 - t419 * t463;
t425 = -t400 * t457 - t458;
t422 = qJD(6) * t371 - t347 * t404 + t370 * t455;
t363 = -t402 * t437 + (t425 * qJD(2) - t424 * qJD(3)) * t399;
t351 = -t375 * qJD(3) + t387 * t461 + t388 * t405;
t350 = t376 * qJD(3) - t387 * t460 + t388 * t402;
t349 = t384 * t454 + ((t425 * qJD(3) + qJD(4) * t464) * t401 + (-t424 * t401 - t404 * t462) * qJD(2)) * t399;
t345 = -qJD(3) * t369 + t385 * t402 + t386 * t460;
t342 = t371 * t395 - t396 * t470;
t341 = -t371 * t396 - t395 * t470;
t340 = -t372 * qJD(4) + t364 * t404 + t401 * t439;
t326 = t427 * qJD(4) + t351 * t404 - t387 * t465;
t325 = (t385 * t461 + t386 * t405 + (t421 * t402 + t420 * t460) * qJD(3)) * t401 + t374 * t454 + t385 * t463 - t420 * t398 * t455;
t322 = -t354 * qJD(4) + t346 * t404 - t386 * t465;
t318 = 0.1e1 / t320;
t309 = 0.1e1 / t311;
t308 = t335 * t485;
t307 = t429 * t335;
t306 = t430 * t335;
t302 = (-t333 * t368 + t334 * t380) * t401 + t433 * t308;
t301 = t433 * t307 - t333 * t359 + t334 * t377;
t300 = t433 * t306 - t333 * t356 + t334 * t373;
t298 = t429 * t484 + (t377 * t440 - t325 * t365 + (t321 * t377 + t339 * t359 + t349 * t354) * t366) * t335;
t296 = t430 * t484 + (t373 * t440 - t322 * t365 + (t321 * t373 + t339 * t356 + t340 * t354) * t366) * t335;
t295 = t484 * t485 + (t428 * t454 + (t380 * t440 - t345 * t365 + (t321 * t380 + t339 * t368 + t354 * t363) * t366) * t401) * t335;
t1 = [0, t298, t295, t296, 0, 0; 0 (t301 * t478 + t314 * t427) * t452 + ((t361 * qJD(4) + t351 * t401 + t387 * t463) * t314 + t301 * t442 + (t427 * t299 - t301 * t323 - (-t298 * t354 - t307 * t321 + t349 + (-t307 * t372 - t359) * t304) * t474 - (-t298 * t372 - t307 * t339 - t325 + (t307 * t354 - t377) * t304) * t475) * t315) * t309 (t302 * t478 + t314 * t471) * t452 + ((-t347 * t401 - t370 * t454) * t314 + (-t477 + t442) * t302 + (t471 * t299 - (t380 * t454 - t295 * t354 - t308 * t321 + t363 * t401 + (-t308 * t372 - t368 * t401) * t304) * t474 - (-t368 * t454 - t295 * t372 - t308 * t339 - t345 * t401 + (t308 * t354 - t380 * t401) * t304) * t475) * t315) * t309 (t300 * t478 - t314 * t358) * t452 + (t300 * t442 + t324 * t314 + (-t358 * t299 - t300 * t323 - (-t296 * t354 - t306 * t321 + t340 + (-t306 * t372 - t356) * t304) * t474 - (-t296 * t372 - t306 * t339 - t322 + (t306 * t354 - t373) * t304) * t475) * t315) * t309, 0, 0; 0 (-t328 * t343 + t344 * t476) * t450 + ((t344 * qJD(6) + t326 * t395 - t350 * t396) * t328 + t344 * t441 + (-t343 * t313 - (-t343 * qJD(6) + t326 * t396 + t350 * t395) * t331 - t344 * t312) * t329) * t318 (-t328 * t341 + t342 * t476) * t450 + (t342 * t441 - t436 * t328 * t396 + t422 * t467 + (-t436 * t331 * t395 - t342 * t312 - t341 * t313 - t422 * t466) * t329) * t318, t431 * t357 * t451 + (t431 * t323 + ((-qJD(6) * t328 - 0.2e1 * t448) * t396 + (t312 * t396 + (t313 - t453) * t395) * t329) * t357) * t318, 0, t451 + 0.2e1 * (t312 * t329 * t318 + (-t318 * t479 - t329 * t480) * t331) * t331;];
JaD_rot  = t1;
