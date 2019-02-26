% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRR8_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:08:08
% EndTime: 2019-02-26 20:08:11
% DurationCPUTime: 2.73s
% Computational Cost: add. (14691->203), mult. (44088->385), div. (816->12), fcn. (56595->17), ass. (0->166)
t395 = sin(qJ(6));
t484 = sin(pkin(12));
t486 = cos(pkin(6));
t434 = t486 * t484;
t485 = cos(pkin(12));
t487 = sin(qJ(2));
t489 = cos(qJ(2));
t388 = -t434 * t487 + t485 * t489;
t394 = cos(pkin(7));
t397 = sin(qJ(3));
t416 = t434 * t489 + t485 * t487;
t392 = sin(pkin(7));
t393 = sin(pkin(6));
t449 = t393 * t484;
t439 = t392 * t449;
t488 = cos(qJ(3));
t368 = t388 * t488 + (-t394 * t416 + t439) * t397;
t398 = cos(qJ(6));
t472 = t368 * t398;
t396 = sin(qJ(5));
t399 = cos(qJ(5));
t422 = -t392 * t416 - t394 * t449;
t454 = t394 * t488;
t490 = -t388 * t397 - t416 * t454 + t488 * t439;
t492 = t396 * t490 + t422 * t399;
t326 = -t395 * t492 - t472;
t496 = 0.2e1 * t326;
t435 = t486 * t485;
t414 = -t489 * t435 + t484 * t487;
t415 = t435 * t487 + t484 * t489;
t450 = t393 * t485;
t440 = t392 * t450;
t367 = t415 * t488 + (-t394 * t414 - t440) * t397;
t382 = t414 * qJD(2);
t383 = t415 * qJD(2);
t340 = qJD(3) * t367 - t382 * t397 + t383 * t454;
t410 = t414 * t488;
t407 = t394 * t410 + t397 * t415 + t440 * t488;
t409 = t392 * t414 - t394 * t450;
t351 = t407 * t396 + t409 * t399;
t468 = t392 * t396;
t316 = -t351 * qJD(5) + t340 * t399 - t383 * t468;
t349 = t396 * t409 - t399 * t407;
t347 = t349 ^ 2;
t443 = t489 * t488;
t452 = t487 * t397;
t419 = t394 * t443 - t452;
t451 = t392 * t486;
t430 = t488 * t451;
t377 = -t393 * t419 - t430;
t469 = t392 * t393;
t444 = t489 * t469;
t386 = t394 * t486 - t444;
t432 = t377 * t399 - t386 * t396;
t365 = 0.1e1 / t432 ^ 2;
t332 = t347 * t365 + 0.1e1;
t330 = 0.1e1 / t332;
t442 = t487 * t488;
t453 = t489 * t397;
t417 = t394 * t442 + t453;
t420 = t394 * t453 + t442;
t441 = t397 * t451;
t359 = qJD(3) * t441 + (qJD(2) * t417 + qJD(3) * t420) * t393;
t370 = t377 * t396 + t386 * t399;
t445 = t487 * t469;
t431 = qJD(2) * t445;
t334 = qJD(5) * t370 - t359 * t399 + t396 * t431;
t364 = 0.1e1 / t432;
t474 = t349 * t365;
t299 = (-t316 * t364 + t334 * t474) * t330;
t333 = atan2(-t349, -t432);
t328 = sin(t333);
t329 = cos(t333);
t433 = t328 * t432 - t329 * t349;
t294 = t299 * t433 + t328 * t316 + t329 * t334;
t312 = -t328 * t349 - t329 * t432;
t309 = 0.1e1 / t312;
t310 = 0.1e1 / t312 ^ 2;
t495 = t294 * t309 * t310;
t352 = -t396 * t422 + t399 * t490;
t446 = 0.2e1 * t352 * t495;
t475 = t334 * t364 * t365;
t494 = (-t316 * t474 + t347 * t475) / t332 ^ 2;
t378 = t393 * t420 + t441;
t425 = t364 * t367 + t378 * t474;
t493 = t399 * t425;
t491 = -0.2e1 * t494;
t327 = t368 * t395 - t398 * t492;
t323 = 0.1e1 / t327;
t324 = 0.1e1 / t327 ^ 2;
t384 = t416 * qJD(2);
t385 = t388 * qJD(2);
t342 = qJD(3) * t368 - t384 * t397 + t385 * t454;
t467 = t392 * t399;
t319 = -qJD(5) * t352 + t342 * t396 + t385 * t467;
t466 = t394 * t397;
t343 = t490 * qJD(3) - t384 * t488 - t385 * t466;
t307 = qJD(6) * t327 + t319 * t395 - t343 * t398;
t322 = t326 ^ 2;
t315 = t322 * t324 + 0.1e1;
t478 = t324 * t326;
t463 = qJD(6) * t326;
t308 = t319 * t398 + t343 * t395 - t463;
t482 = t308 * t323 * t324;
t483 = (t307 * t478 - t322 * t482) / t315 ^ 2;
t481 = t310 * t352;
t313 = 0.1e1 / t315;
t480 = t313 * t324;
t318 = t492 * qJD(5) + t342 * t399 - t385 * t468;
t479 = t318 * t310;
t477 = t328 * t352;
t476 = t329 * t352;
t473 = t368 * t396;
t471 = t368 * t399;
t465 = qJD(5) * t396;
t464 = qJD(5) * t399;
t348 = t352 ^ 2;
t306 = t310 * t348 + 0.1e1;
t462 = 0.2e1 * (-t348 * t495 - t352 * t479) / t306 ^ 2;
t461 = -0.2e1 * t483;
t460 = 0.2e1 * t483;
t458 = t324 * t483;
t457 = t307 * t480;
t456 = t326 * t482;
t455 = t349 * t475;
t448 = 0.2e1 * t456;
t447 = 0.2e1 * t455;
t436 = qJD(6) * t473 + t342;
t423 = -t388 * t454 + t397 * t416;
t356 = t388 * t467 - t396 * t423;
t373 = -t388 * t466 - t416 * t488;
t339 = t356 * t398 + t373 * t395;
t338 = t356 * t395 - t373 * t398;
t428 = t323 * t395 - t398 * t478;
t427 = t351 * t364 + t370 * t474;
t411 = t394 * t415;
t371 = -t397 * t414 + t411 * t488;
t412 = t392 * t415;
t354 = t371 * t399 - t396 * t412;
t381 = t417 * t393;
t374 = -t381 * t399 + t396 * t445;
t426 = -t354 * t364 + t374 * t474;
t355 = t388 * t468 + t399 * t423;
t418 = -t394 * t452 + t443;
t413 = qJD(6) * t490 + t343 * t396 + t368 * t464;
t360 = qJD(3) * t430 + (qJD(2) * t418 + qJD(3) * t419) * t393;
t346 = qJD(3) * t423 + t384 * t466 - t385 * t488;
t345 = qJD(3) * t373 - t384 * t454 - t385 * t397;
t344 = qJD(2) * t396 * t444 + t445 * t464 - (qJD(2) * t419 + qJD(3) * t418) * t393 * t399 + t381 * t465;
t341 = -qJD(3) * t407 - t382 * t488 - t383 * t466;
t337 = t395 * t490 + t396 * t472;
t335 = qJD(5) * t432 + t359 * t396 + t399 * t431;
t321 = -qJD(5) * t355 + t345 * t396 - t384 * t467;
t320 = t382 * t468 - t412 * t464 + (-t382 * t454 - t383 * t397 + (-t397 * t411 - t410) * qJD(3)) * t399 - t371 * t465;
t317 = qJD(5) * t349 - t340 * t396 - t383 * t467;
t304 = 0.1e1 / t306;
t303 = t330 * t493;
t302 = t426 * t330;
t301 = t427 * t330;
t297 = (t328 * t367 - t329 * t378) * t399 - t433 * t303;
t296 = t302 * t433 + t328 * t354 + t329 * t374;
t295 = t301 * t433 - t328 * t351 + t329 * t370;
t293 = t426 * t491 + (t374 * t447 - t320 * t364 + (-t316 * t374 - t334 * t354 + t344 * t349) * t365) * t330;
t291 = t427 * t491 + (t370 * t447 - t317 * t364 + (-t316 * t370 + t334 * t351 + t335 * t349) * t365) * t330;
t290 = 0.2e1 * t493 * t494 + (t425 * t465 + (-0.2e1 * t378 * t455 - t341 * t364 + (t316 * t378 - t334 * t367 - t349 * t360) * t365) * t399) * t330;
t1 = [0, t293, t290, 0, t291, 0; 0 (t296 * t481 - t309 * t355) * t462 + ((qJD(5) * t356 - t345 * t399 - t384 * t468) * t309 + t296 * t446 + (-t355 * t294 + t296 * t318 - (-t293 * t349 + t302 * t316 + t344 + (t302 * t432 + t354) * t299) * t476 - (t293 * t432 - t302 * t334 + t320 + (t302 * t349 - t374) * t299) * t477) * t310) * t304 (t297 * t481 + t309 * t471) * t462 + ((-t343 * t399 + t368 * t465) * t309 + (t479 + t446) * t297 + (t471 * t294 - (t378 * t465 - t290 * t349 - t303 * t316 - t360 * t399 + (-t303 * t432 + t367 * t399) * t299) * t476 - (-t367 * t465 + t290 * t432 + t303 * t334 + t341 * t399 + (-t303 * t349 + t378 * t399) * t299) * t477) * t310) * t304, 0 (t295 * t481 + t309 * t492) * t462 + (t295 * t446 + t319 * t309 + (t492 * t294 + t295 * t318 - (-t291 * t349 + t301 * t316 + t335 + (t301 * t432 - t351) * t299) * t476 - (t291 * t432 - t301 * t334 + t317 + (t301 * t349 - t370) * t299) * t477) * t310) * t304, 0; 0 (-t323 * t338 + t339 * t478) * t460 + ((qJD(6) * t339 + t321 * t395 - t346 * t398) * t323 + t339 * t448 + (-t338 * t308 - (-qJD(6) * t338 + t321 * t398 + t346 * t395) * t326 - t339 * t307) * t324) * t313 (t458 * t496 - t457) * t337 + (-t308 * t480 + t323 * t461) * (t395 * t473 - t398 * t490) + ((t395 * t413 + t398 * t436) * t323 - (-t395 * t436 + t398 * t413) * t478 + t337 * t448) * t313, 0, t428 * t352 * t460 + (t428 * t318 + ((-qJD(6) * t323 - 0.2e1 * t456) * t398 + (t307 * t398 + (t308 - t463) * t395) * t324) * t352) * t313, t461 + (t457 + (-t313 * t482 - t458) * t326) * t496;];
JaD_rot  = t1;
