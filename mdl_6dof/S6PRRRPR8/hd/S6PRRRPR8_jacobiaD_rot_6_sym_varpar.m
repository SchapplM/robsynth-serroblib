% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRPR8_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR8_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:14:38
% EndTime: 2019-02-26 20:14:40
% DurationCPUTime: 2.74s
% Computational Cost: add. (14691->204), mult. (44088->386), div. (816->12), fcn. (56595->17), ass. (0->168)
t400 = sin(qJ(4));
t404 = cos(qJ(4));
t394 = sin(pkin(7));
t401 = sin(qJ(3));
t405 = cos(qJ(3));
t398 = cos(pkin(6));
t406 = cos(qJ(2));
t489 = sin(pkin(12));
t444 = t489 * t406;
t396 = cos(pkin(12));
t402 = sin(qJ(2));
t469 = t396 * t402;
t423 = -t398 * t469 - t444;
t397 = cos(pkin(7));
t395 = sin(pkin(6));
t470 = t395 * t396;
t452 = t394 * t470;
t445 = t489 * t402;
t468 = t396 * t406;
t491 = -t398 * t468 + t445;
t428 = -t397 * t491 - t452;
t415 = t428 * t401 - t405 * t423;
t418 = t423 * qJD(2);
t493 = qJD(4) * t415 + t394 * t418;
t377 = t394 * t491 - t397 * t470;
t466 = t397 * t405;
t447 = qJD(3) * t466;
t456 = t423 * qJD(3);
t467 = t397 * t401;
t384 = t491 * qJD(2);
t494 = qJD(3) * t452 + t384;
t495 = t377 * qJD(4) + t401 * t456 - t494 * t405 + t418 * t467 - t447 * t491;
t317 = -t493 * t400 + t495 * t404;
t351 = t377 * t400 + t415 * t404;
t348 = t351 ^ 2;
t462 = t402 * t405;
t463 = t401 * t406;
t425 = t397 * t463 + t462;
t475 = t394 * t398;
t376 = t425 * t395 + t401 * t475;
t471 = t394 * t406;
t387 = -t395 * t471 + t398 * t397;
t368 = t376 * t404 + t387 * t400;
t362 = 0.1e1 / t368 ^ 2;
t332 = t348 * t362 + 0.1e1;
t330 = 0.1e1 / t332;
t460 = t405 * t406;
t464 = t401 * t402;
t424 = -t397 * t464 + t460;
t427 = t397 * t460 - t464;
t449 = qJD(3) * t475;
t360 = t405 * t449 + (t424 * qJD(2) + t427 * qJD(3)) * t395;
t367 = -t376 * t400 + t387 * t404;
t473 = t394 * t402;
t451 = t395 * t473;
t440 = qJD(2) * t451;
t335 = t367 * qJD(4) + t360 * t404 + t400 * t440;
t361 = 0.1e1 / t368;
t480 = t351 * t362;
t299 = (-t317 * t361 + t335 * t480) * t330;
t333 = atan2(-t351, t368);
t328 = sin(t333);
t329 = cos(t333);
t436 = -t328 * t368 - t329 * t351;
t294 = t436 * t299 - t328 * t317 + t329 * t335;
t312 = -t328 * t351 + t329 * t368;
t309 = 0.1e1 / t312;
t310 = 0.1e1 / t312 ^ 2;
t498 = t294 * t309 * t310;
t421 = t398 * t444 + t469;
t446 = t395 * t489;
t438 = t394 * t446;
t413 = -t421 * t397 + t438;
t422 = t398 * t445 - t468;
t366 = t413 * t401 - t405 * t422;
t378 = t421 * t394 + t397 * t446;
t353 = t366 * t400 - t378 * t404;
t419 = t421 * t405;
t476 = t422 * t401;
t365 = t397 * t419 - t405 * t438 - t476;
t399 = sin(qJ(6));
t403 = cos(qJ(6));
t435 = t353 * t403 - t365 * t399;
t497 = t435 * qJD(6);
t354 = t366 * t404 + t378 * t400;
t441 = 0.2e1 * t354 * t498;
t364 = t401 * t423 + t428 * t405;
t375 = t427 * t395 + t405 * t475;
t430 = -t361 * t364 + t375 * t480;
t496 = t404 * t430;
t492 = qJD(3) * t401 * t491 + t405 * t418;
t481 = t335 * t361 * t362;
t490 = -0.2e1 * (t317 * t480 - t348 * t481) / t332 ^ 2;
t478 = t365 * t403;
t327 = t353 * t399 + t478;
t323 = 0.1e1 / t327;
t324 = 0.1e1 / t327 ^ 2;
t385 = t421 * qJD(2);
t386 = t422 * qJD(2);
t343 = t386 * t467 - t385 * t405 + (t413 * t405 + t476) * qJD(3);
t472 = t394 * t404;
t318 = t354 * qJD(4) + t343 * t400 + t386 * t472;
t342 = t366 * qJD(3) - t385 * t401 - t386 * t466;
t307 = t327 * qJD(6) - t318 * t403 + t342 * t399;
t322 = t435 ^ 2;
t315 = t322 * t324 + 0.1e1;
t484 = t324 * t435;
t308 = t318 * t399 + t342 * t403 + t497;
t487 = t308 * t323 * t324;
t488 = (-t307 * t484 - t322 * t487) / t315 ^ 2;
t486 = t310 * t354;
t474 = t394 * t400;
t319 = -t353 * qJD(4) + t343 * t404 - t386 * t474;
t485 = t319 * t310;
t483 = t328 * t354;
t482 = t329 * t354;
t479 = t365 * t400;
t477 = t365 * t404;
t465 = t399 * t435;
t461 = t403 * t323;
t459 = qJD(4) * t400;
t458 = qJD(4) * t404;
t349 = t354 ^ 2;
t306 = t349 * t310 + 0.1e1;
t455 = 0.2e1 * (-t349 * t498 + t354 * t485) / t306 ^ 2;
t454 = 0.2e1 * t488;
t443 = -0.2e1 * t435 * t487;
t442 = -0.2e1 * t351 * t481;
t437 = -qJD(6) * t479 + t343;
t370 = -t421 * t401 - t422 * t466;
t371 = t422 * t467 - t419;
t429 = -t371 * t400 - t422 * t472;
t434 = -t370 * t399 - t403 * t429;
t339 = t370 * t403 - t399 * t429;
t433 = -t324 * t465 + t461;
t350 = -t377 * t404 + t415 * t400;
t432 = t350 * t361 + t367 * t480;
t369 = -t405 * t491 + t423 * t467;
t355 = t369 * t404 - t423 * t474;
t383 = t424 * t395;
t372 = t383 * t404 + t400 * t451;
t431 = -t355 * t361 + t372 * t480;
t357 = t371 * t404 - t422 * t474;
t426 = -t397 * t462 - t463;
t420 = qJD(6) * t366 + t342 * t400 + t365 * t458;
t359 = -t401 * t449 + (t426 * qJD(2) - t425 * qJD(3)) * t395;
t347 = -t370 * qJD(3) + t385 * t467 + t386 * t405;
t346 = t371 * qJD(3) - t385 * t466 + t386 * t401;
t344 = -t383 * t459 + ((t426 * qJD(3) + qJD(4) * t473) * t404 + (t400 * t471 - t425 * t404) * qJD(2)) * t395;
t340 = t492 * t397 + t494 * t401 + t405 * t456;
t337 = t366 * t403 - t399 * t479;
t336 = t366 * t399 + t400 * t478;
t334 = -t368 * qJD(4) - t360 * t400 + t404 * t440;
t321 = t357 * qJD(4) + t347 * t400 + t385 * t472;
t320 = (t384 * t467 + t423 * t447 + t492) * t404 - t369 * t459 - t384 * t474 - t423 * t394 * t458;
t316 = t495 * t400 + t493 * t404;
t313 = 0.1e1 / t315;
t304 = 0.1e1 / t306;
t303 = t330 * t496;
t302 = t431 * t330;
t301 = t432 * t330;
t297 = (-t328 * t364 + t329 * t375) * t404 + t436 * t303;
t296 = t436 * t302 - t328 * t355 + t329 * t372;
t295 = t436 * t301 + t328 * t350 + t329 * t367;
t293 = t431 * t490 + (t372 * t442 - t320 * t361 + (t317 * t372 + t335 * t355 + t344 * t351) * t362) * t330;
t291 = t432 * t490 + (t367 * t442 + t316 * t361 + (t317 * t367 + t334 * t351 - t335 * t350) * t362) * t330;
t290 = t490 * t496 + (-t430 * t459 + (t375 * t442 - t340 * t361 + (t317 * t375 + t335 * t364 + t351 * t359) * t362) * t404) * t330;
t1 = [0, t293, t290, t291, 0, 0; 0 (t296 * t486 - t309 * t357) * t455 + ((t429 * qJD(4) + t347 * t404 - t385 * t474) * t309 + t296 * t441 + (-t357 * t294 - t296 * t319 - (-t293 * t351 - t302 * t317 + t344 + (-t302 * t368 - t355) * t299) * t482 - (-t293 * t368 - t302 * t335 - t320 + (t302 * t351 - t372) * t299) * t483) * t310) * t304 (t297 * t486 + t309 * t477) * t455 + ((-t342 * t404 + t365 * t459) * t309 + (-t485 + t441) * t297 + (t477 * t294 - (-t375 * t459 - t290 * t351 - t303 * t317 + t359 * t404 + (-t303 * t368 - t364 * t404) * t299) * t482 - (t364 * t459 - t290 * t368 - t303 * t335 - t340 * t404 + (t303 * t351 - t375 * t404) * t299) * t483) * t310) * t304 (t295 * t486 + t309 * t353) * t455 + (t295 * t441 - t318 * t309 + (t353 * t294 - t295 * t319 - (-t291 * t351 - t301 * t317 + t334 + (-t301 * t368 + t350) * t299) * t482 - (-t291 * t368 - t301 * t335 + t316 + (t301 * t351 - t367) * t299) * t483) * t310) * t304, 0, 0; 0 (t323 * t434 - t339 * t484) * t454 + ((t339 * qJD(6) - t321 * t403 + t346 * t399) * t323 + t339 * t443 + (t434 * t308 + (t434 * qJD(6) + t321 * t399 + t346 * t403) * t435 - t339 * t307) * t324) * t313 (-t323 * t336 - t337 * t484) * t454 + (t337 * t443 + t437 * t323 * t399 + t420 * t461 + (t403 * t435 * t437 - t337 * t307 - t336 * t308 - t420 * t465) * t324) * t313, t433 * t354 * t454 + (-t433 * t319 + ((qJD(6) * t323 + t443) * t399 + (-t307 * t399 + (t308 + t497) * t403) * t324) * t354) * t313, 0, -0.2e1 * t488 - 0.2e1 * (t307 * t324 * t313 - (-t313 * t487 - t324 * t488) * t435) * t435;];
JaD_rot  = t1;
