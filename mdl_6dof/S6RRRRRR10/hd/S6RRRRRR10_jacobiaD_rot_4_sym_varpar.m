% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRRR10
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRR10_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobiaD_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:52:57
% EndTime: 2019-02-26 22:53:01
% DurationCPUTime: 3.13s
% Computational Cost: add. (13117->203), mult. (39618->394), div. (699->12), fcn. (49954->17), ass. (0->172)
t400 = cos(pkin(6));
t406 = cos(qJ(2));
t504 = sin(qJ(1));
t452 = t504 * t406;
t403 = sin(qJ(2));
t407 = cos(qJ(1));
t466 = t407 * t403;
t387 = t400 * t466 + t452;
t402 = sin(qJ(3));
t405 = cos(qJ(3));
t396 = sin(pkin(7));
t397 = sin(pkin(6));
t476 = t397 * t407;
t455 = t396 * t476;
t399 = cos(pkin(7));
t394 = t504 * t403;
t467 = t406 * t407;
t505 = -t400 * t467 + t394;
t506 = t505 * t399;
t422 = t506 + t455;
t369 = t387 * t402 + t422 * t405;
t385 = -t396 * t505 + t399 * t476;
t395 = sin(pkin(8));
t398 = cos(pkin(8));
t512 = -t369 * t395 + t385 * t398;
t484 = t387 * t405;
t511 = t422 * t402 - t484;
t428 = t400 * t452 + t466;
t379 = t428 * qJD(1) + t387 * qJD(2);
t453 = t397 * t504;
t444 = qJD(1) * t453;
t424 = t379 * t396 + t399 * t444;
t510 = t424 * t398;
t445 = t400 * t394;
t380 = -qJD(1) * t445 - qJD(2) * t394 + (qJD(2) * t400 + qJD(1)) * t467;
t440 = t396 * t444;
t423 = t379 * t399 - t440;
t478 = t396 * t405;
t450 = qJD(3) * t478;
t509 = (qJD(3) * t506 - t380) * t405 + t450 * t476 + (t387 * qJD(3) + t423) * t402;
t465 = qJD(3) * t395;
t472 = t399 * t405;
t486 = t380 * t402;
t323 = (-t379 * t472 + t405 * t440 - t486) * t395 - t510 + t511 * t465;
t350 = t512 ^ 2;
t468 = t405 * t406;
t471 = t402 * t403;
t432 = t399 * t468 - t471;
t477 = t396 * t406;
t365 = -(t432 * t397 + t400 * t478) * t395 + (-t397 * t477 + t399 * t400) * t398;
t363 = 0.1e1 / t365 ^ 2;
t342 = t350 * t363 + 0.1e1;
t362 = 0.1e1 / t365;
t469 = t403 * t405;
t470 = t402 * t406;
t479 = t396 * t398;
t420 = -(-t399 * t469 - t470) * t395 + t403 * t479;
t431 = -t399 * t470 - t469;
t457 = t396 * t400 * t402;
t356 = t457 * t465 + (t420 * qJD(2) - t431 * t465) * t397;
t491 = t356 * t363;
t490 = t362 * t491;
t492 = t512 * t363;
t500 = (t323 * t492 - t350 * t490) / t342 ^ 2;
t508 = -0.2e1 * t500;
t417 = t396 * t453 - t428 * t399;
t429 = t445 - t467;
t372 = t417 * t402 - t405 * t429;
t343 = atan2(t512, t365);
t338 = sin(t343);
t339 = cos(t343);
t318 = t338 * t512 + t339 * t365;
t315 = 0.1e1 / t318;
t401 = sin(qJ(4));
t404 = cos(qJ(4));
t483 = t429 * t402;
t371 = t417 * t405 + t483;
t418 = -t428 * t396 - t399 * t453;
t415 = t418 * t395;
t414 = t371 * t398 - t415;
t337 = t372 * t404 + t414 * t401;
t326 = 0.1e1 / t337;
t316 = 0.1e1 / t318 ^ 2;
t327 = 0.1e1 / t337 ^ 2;
t354 = t371 * t395 + t418 * t398;
t351 = t354 ^ 2;
t312 = t316 * t351 + 0.1e1;
t378 = t387 * qJD(1) + t428 * qJD(2);
t377 = t505 * qJD(1) + t429 * qJD(2);
t451 = qJD(1) * t476;
t426 = t377 * t399 + t396 * t451;
t329 = -t372 * qJD(3) + t378 * t402 + t426 * t405;
t427 = -t377 * t396 + t399 * t451;
t322 = t329 * t395 - t427 * t398;
t497 = t322 * t316;
t340 = 0.1e1 / t342;
t436 = t323 * t362 - t491 * t512;
t306 = t436 * t340;
t442 = -t338 * t365 + t339 * t512;
t301 = t442 * t306 + t323 * t338 + t339 * t356;
t502 = t301 * t315 * t316;
t503 = (-t351 * t502 + t354 * t497) / t312 ^ 2;
t330 = qJD(3) * t483 + t426 * t402 + (t417 * qJD(3) - t378) * t405;
t419 = t329 * t398 + t427 * t395;
t313 = t337 * qJD(4) + t330 * t401 - t419 * t404;
t474 = t398 * t404;
t488 = t372 * t401;
t336 = -t371 * t474 + t404 * t415 + t488;
t325 = t336 ^ 2;
t321 = t325 * t327 + 0.1e1;
t496 = t327 * t336;
t314 = t330 * t404 + t419 * t401 + (t414 * t404 - t488) * qJD(4);
t499 = t314 * t326 * t327;
t501 = (t313 * t496 - t325 * t499) / t321 ^ 2;
t498 = t316 * t354;
t495 = t338 * t354;
t494 = t339 * t354;
t493 = t512 * t362;
t446 = t402 * t455 - t484;
t370 = t402 * t506 + t446;
t489 = t370 * t401;
t473 = t399 * t402;
t375 = -t428 * t405 + t429 * t473;
t487 = t375 * t401;
t481 = t395 * t401;
t480 = t395 * t404;
t475 = t398 * t401;
t464 = -0.2e1 * t503;
t463 = 0.2e1 * t503;
t462 = -0.2e1 * t502;
t461 = 0.2e1 * t501;
t460 = 0.2e1 * t500;
t459 = t512 * t490;
t458 = t396 * t480;
t449 = t362 * t508;
t448 = 0.2e1 * t336 * t499;
t447 = t354 * t462;
t441 = t369 * t398 + t385 * t395;
t437 = t505 * t402;
t357 = (-t387 * t472 + t437) * t395 - t387 * t479;
t376 = t420 * t397;
t435 = -t357 * t362 + t376 * t492;
t368 = t399 * t437 + t446;
t384 = t431 * t397 - t457;
t434 = t362 * t368 + t384 * t492;
t348 = t371 * t401 + t372 * t474;
t349 = t371 * t404 - t372 * t475;
t374 = t428 * t402 + t429 * t472;
t433 = -t395 * t396 * t429 + t374 * t398;
t430 = t399 * t471 - t468;
t425 = t338 + (t339 * t493 - t338) * t340;
t335 = t370 * t404 + t441 * t401;
t347 = t375 * t404 + t433 * t401;
t361 = -t400 * t450 + (t430 * qJD(2) - t432 * qJD(3)) * t397;
t359 = (-t430 * t465 + (t432 * t395 + t398 * t477) * qJD(2)) * t397;
t358 = -t374 * t395 - t429 * t479;
t346 = -t433 * t404 + t487;
t345 = t374 * qJD(3) + t377 * t405 + t378 * t473;
t344 = -t375 * qJD(3) - t377 * t402 + t378 * t472;
t334 = -t441 * t404 + t489;
t332 = -t511 * qJD(3) + t423 * t405 + t486;
t324 = (-t380 * t472 + t379 * t402 + (t387 * t473 + t405 * t505) * qJD(3)) * t395 - t380 * t479;
t319 = 0.1e1 / t321;
t310 = 0.1e1 / t312;
t308 = t434 * t395 * t340;
t307 = t435 * t340;
t305 = t425 * t354;
t303 = (t338 * t368 - t339 * t384) * t395 + t442 * t308;
t302 = -t442 * t307 + t338 * t357 + t339 * t376;
t300 = (t434 * t508 + (-0.2e1 * t384 * t459 + t509 * t362 + (t323 * t384 - t356 * t368 + t361 * t512) * t363) * t340) * t395;
t299 = t435 * t460 + (0.2e1 * t376 * t459 + t324 * t362 + (-t323 * t376 - t356 * t357 - t359 * t512) * t363) * t340;
t1 = [t354 * t449 + (t322 * t362 - t354 * t491) * t340, t299, t300, 0, 0, 0; t512 * t315 * t464 + ((-t332 * t395 - t510) * t315 + (-t301 * t512 + t305 * t322) * t316) * t310 + ((t305 * t462 + t425 * t497) * t310 + (t305 * t464 + ((-t306 * t340 * t493 + t460) * t495 + (t512 * t449 + t306 + (-t306 + t436) * t340) * t494) * t310) * t316) * t354 (-t302 * t498 - t315 * t358) * t463 + ((-t344 * t395 - t378 * t479) * t315 + t302 * t447 + (-t358 * t301 + t302 * t322 + (t299 * t512 - t307 * t323 + t359 + (t307 * t365 + t357) * t306) * t494 + (-t299 * t365 + t307 * t356 + t324 + (t307 * t512 - t376) * t306) * t495) * t316) * t310 (-t315 * t372 * t395 - t303 * t498) * t463 + ((t442 * t300 + (-t318 * t306 + t323 * t339 - t338 * t356) * t308) * t498 + (t447 + t497) * t303 + (t330 * t315 + (-t372 * t301 + (t509 * t338 - t339 * t361 + (t338 * t384 + t339 * t368) * t306) * t354) * t316) * t395) * t310, 0, 0, 0; (-t326 * t334 + t335 * t496) * t461 + ((-t332 * t474 + t401 * t509 + t424 * t480) * t326 + t335 * t448 + (-t334 * t314 - (t332 * t475 + t404 * t509 - t424 * t481) * t336 - t335 * t313) * t327 + (t335 * t326 - (t369 * t474 + t385 * t480 - t489) * t496) * qJD(4)) * t319 (-t326 * t346 + t347 * t496) * t461 + ((-t344 * t474 + t345 * t401 + t378 * t458) * t326 + t347 * t448 + (-t346 * t314 - (-t378 * t396 * t481 + t344 * t475 + t345 * t404) * t336 - t347 * t313) * t327 + (t347 * t326 - (t374 * t474 - t429 * t458 - t487) * t496) * qJD(4)) * t319 (-t326 * t348 + t349 * t496) * t461 + ((t349 * qJD(4) + t329 * t401 + t330 * t474) * t326 + t349 * t448 + (-t348 * t314 - (-t348 * qJD(4) + t329 * t404 - t330 * t475) * t336 - t349 * t313) * t327) * t319, -0.2e1 * t501 + 0.2e1 * (t313 * t327 * t319 + (-t319 * t499 - t327 * t501) * t336) * t336, 0, 0;];
JaD_rot  = t1;
