% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR14_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR14_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_jacobiaD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:38:12
% EndTime: 2019-02-26 22:38:18
% DurationCPUTime: 5.57s
% Computational Cost: add. (18619->234), mult. (55569->452), div. (959->12), fcn. (70289->17), ass. (0->190)
t523 = cos(pkin(6));
t524 = sin(qJ(2));
t469 = t523 * t524;
t422 = cos(qJ(2));
t525 = sin(qJ(1));
t482 = t525 * t422;
t526 = cos(qJ(1));
t403 = t469 * t526 + t482;
t442 = t482 * t523 + t524 * t526;
t393 = qJD(1) * t442 + qJD(2) * t403;
t412 = t525 * t524;
t455 = t525 * t469;
t470 = t523 * t526;
t394 = -qJD(1) * t455 - qJD(2) * t412 + (qJD(1) * t526 + qJD(2) * t470) * t422;
t417 = cos(pkin(7));
t419 = sin(qJ(3));
t421 = cos(qJ(3));
t414 = sin(pkin(7));
t415 = sin(pkin(6));
t483 = t415 * t525;
t471 = qJD(1) * t483;
t464 = t414 * t471;
t532 = -t422 * t470 + t412;
t444 = t532 * t417;
t484 = t415 * t526;
t474 = t414 * t484;
t435 = t444 + t474;
t529 = t403 * t419 + t421 * t435;
t335 = qJD(3) * t529 - (-t393 * t417 + t464) * t419 - t394 * t421;
t506 = t403 * t421;
t383 = t419 * t435 - t506;
t418 = sin(qJ(4));
t420 = cos(qJ(4));
t437 = t414 * t532 - t417 * t484;
t364 = t383 * t418 + t420 * t437;
t447 = t393 * t414 + t417 * t471;
t326 = qJD(4) * t364 - t335 * t420 + t418 * t447;
t365 = t383 * t420 - t418 * t437;
t539 = qJD(4) * t365 + t335 * t418 + t420 * t447;
t481 = t524 * t419;
t497 = t421 * t422;
t451 = -t417 * t481 + t497;
t452 = t417 * t497 - t481;
t479 = t414 * t523;
t468 = qJD(3) * t479;
t372 = t421 * t468 + (qJD(2) * t451 + qJD(3) * t452) * t415;
t480 = t524 * t421;
t498 = t419 * t422;
t453 = t417 * t498 + t480;
t398 = t415 * t453 + t419 * t479;
t502 = t414 * t422;
t402 = -t415 * t502 + t417 * t523;
t377 = t398 * t420 + t402 * t418;
t485 = t414 * t524;
t475 = t415 * t485;
t463 = qJD(2) * t475;
t350 = qJD(4) * t377 + t372 * t418 - t420 * t463;
t376 = t398 * t418 - t402 * t420;
t374 = 0.1e1 / t376 ^ 2;
t534 = t350 * t374;
t373 = 0.1e1 / t376;
t397 = t415 * t452 + t421 * t479;
t511 = t364 * t374;
t458 = t373 * t529 - t397 * t511;
t533 = t418 * t458;
t440 = t442 * t421;
t472 = t414 * t483;
t531 = -t417 * t440 + t421 * t472;
t349 = atan2(t364, t376);
t344 = sin(t349);
t345 = cos(t349);
t319 = t344 * t364 + t345 * t376;
t316 = 0.1e1 / t319;
t441 = -t526 * t422 + t455;
t385 = -t441 * t421 + (-t417 * t442 + t472) * t419;
t434 = t414 * t442 + t417 * t483;
t367 = t385 * t420 + t418 * t434;
t384 = -t419 * t441 - t531;
t413 = sin(pkin(13));
t416 = cos(pkin(13));
t343 = t367 * t416 + t384 * t413;
t337 = 0.1e1 / t343;
t317 = 0.1e1 / t319 ^ 2;
t338 = 0.1e1 / t343 ^ 2;
t528 = 0.2e1 * t364;
t433 = t434 * t420;
t366 = t385 * t418 - t433;
t527 = 0.2e1 * t366;
t360 = t366 ^ 2;
t313 = t317 * t360 + 0.1e1;
t392 = qJD(1) * t403 + qJD(2) * t442;
t439 = t441 * qJD(2);
t432 = qJD(1) * t532 + t439;
t431 = t432 * t419;
t465 = qJD(1) * t474;
t496 = qJD(3) * t419;
t331 = qJD(3) * t531 - t392 * t421 + t417 * t431 + t419 * t465 + t441 * t496;
t429 = -t437 * qJD(1) - t414 * t439;
t323 = qJD(4) * t367 + t331 * t418 - t420 * t429;
t517 = t323 * t317;
t359 = t364 ^ 2;
t348 = t359 * t374 + 0.1e1;
t346 = 0.1e1 / t348;
t461 = -t350 * t511 + t373 * t539;
t306 = t461 * t346;
t466 = -t344 * t376 + t345 * t364;
t300 = t306 * t466 + t344 * t539 + t345 * t350;
t318 = t316 * t317;
t521 = t300 * t318;
t522 = (-t360 * t521 + t366 * t517) / t313 ^ 2;
t513 = t373 * t534;
t520 = (-t359 * t513 + t511 * t539) / t348 ^ 2;
t495 = qJD(4) * t418;
t324 = qJD(4) * t433 + t331 * t420 - t385 * t495 + t418 * t429;
t430 = t432 * t421;
t330 = qJD(3) * t385 - t392 * t419 - t417 * t430 - t421 * t465;
t315 = t324 * t416 + t330 * t413;
t519 = t315 * t337 * t338;
t518 = t317 * t366;
t342 = t367 * t413 - t384 * t416;
t516 = t338 * t342;
t515 = t344 * t366;
t514 = t345 * t366;
t512 = t364 * t373;
t510 = t384 * t418;
t509 = t384 * t420;
t505 = t413 * t337;
t504 = t414 * t418;
t503 = t414 * t420;
t501 = t416 * t342;
t500 = t417 * t419;
t499 = t417 * t421;
t494 = qJD(4) * t420;
t493 = 0.2e1 * t522;
t314 = t324 * t413 - t330 * t416;
t336 = t342 ^ 2;
t322 = t336 * t338 + 0.1e1;
t492 = 0.2e1 * (t314 * t516 - t336 * t519) / t322 ^ 2;
t491 = -0.2e1 * t520;
t490 = t318 * t527;
t489 = t373 * t520;
t488 = t342 * t519;
t487 = t317 * t515;
t486 = t317 * t514;
t478 = t300 * t490;
t477 = 0.2e1 * t488;
t476 = t513 * t528;
t460 = t365 * t373 - t377 * t511;
t388 = -t403 * t500 - t421 * t532;
t368 = t388 * t418 - t403 * t503;
t401 = t451 * t415;
t391 = t401 * t418 - t420 * t475;
t459 = -t368 * t373 - t391 * t511;
t390 = t441 * t500 - t440;
t457 = -t390 * t418 - t441 * t503;
t370 = t390 * t420 - t441 * t504;
t456 = -t330 * t420 + t384 * t495;
t450 = -t417 * t480 - t498;
t449 = -t344 + (-t345 * t512 + t344) * t346;
t448 = -t393 * t499 - t394 * t419 + t421 * t464 + t474 * t496;
t445 = t532 * t419;
t389 = -t419 * t442 - t441 * t499;
t371 = -t419 * t468 + (qJD(2) * t450 - qJD(3) * t453) * t415;
t358 = t401 * t494 + ((qJD(3) * t450 + qJD(4) * t485) * t418 + (-t418 * t453 - t420 * t502) * qJD(2)) * t415;
t357 = t370 * t416 + t389 * t413;
t356 = t370 * t413 - t389 * t416;
t355 = t385 * t413 - t416 * t509;
t354 = -t385 * t416 - t413 * t509;
t353 = -qJD(3) * t389 + t392 * t500 + t430;
t352 = qJD(3) * t390 - t392 * t499 + t431;
t351 = -qJD(4) * t376 + t372 * t420 + t418 * t463;
t341 = t365 * t416 - t413 * t529;
t340 = t365 * t413 + t416 * t529;
t334 = (t419 * t444 - t506) * qJD(3) + t448;
t332 = (t417 * t445 - t506) * qJD(3) + t448;
t329 = (-t394 * t500 - t393 * t421 + (-t403 * t499 + t445) * qJD(3)) * t418 + t388 * t494 - t394 * t503 + t403 * t414 * t495;
t328 = qJD(4) * t457 + t353 * t420 - t392 * t504;
t320 = 0.1e1 / t322;
t311 = 0.1e1 / t313;
t310 = t346 * t533;
t309 = t459 * t346;
t308 = t460 * t346;
t305 = t449 * t366;
t303 = (t344 * t529 + t345 * t397) * t418 + t466 * t310;
t302 = t309 * t466 - t344 * t368 + t345 * t391;
t301 = t308 * t466 + t344 * t365 + t345 * t377;
t299 = t459 * t491 + (t391 * t476 - t329 * t373 + (t350 * t368 - t358 * t364 - t391 * t539) * t374) * t346;
t297 = t460 * t491 + (t377 * t476 - t326 * t373 + (-t350 * t365 - t351 * t364 - t377 * t539) * t374) * t346;
t296 = t491 * t533 + (t458 * t494 + (t397 * t476 - t332 * t373 + (-t350 * t529 - t364 * t371 - t397 * t539) * t374) * t418) * t346;
t1 = [t489 * t527 + (-t323 * t373 + t366 * t534) * t346, t299, t296, t297, 0, 0; -0.2e1 * t364 * t316 * t522 + (t539 * t316 + (-t364 * t300 - t305 * t323) * t317) * t311 + (t305 * t317 * t493 + (0.2e1 * t305 * t521 - (t306 * t346 * t512 + t491) * t487 - (t489 * t528 - t306 + (t306 - t461) * t346) * t486 - t449 * t517) * t311) * t366 (t302 * t518 + t316 * t457) * t493 + ((qJD(4) * t370 + t353 * t418 + t392 * t503) * t316 + t302 * t478 + (t457 * t300 - t302 * t323 - (t299 * t364 + t309 * t539 + t358 + (-t309 * t376 - t368) * t306) * t514 - (-t299 * t376 - t309 * t350 - t329 + (-t309 * t364 - t391) * t306) * t515) * t317) * t311 (t303 * t518 + t316 * t510) * t493 + (-t303 * t517 + (-t330 * t418 - t384 * t494) * t316 + (t303 * t490 + t317 * t510) * t300 - (t397 * t494 + t296 * t364 + t310 * t539 + t371 * t418 + (-t310 * t376 + t418 * t529) * t306) * t486 - (t529 * t494 - t296 * t376 - t310 * t350 - t332 * t418 + (-t310 * t364 - t397 * t418) * t306) * t487) * t311 (t301 * t518 - t316 * t367) * t493 + (t301 * t478 + t324 * t316 + (-t367 * t300 - t301 * t323 - (t297 * t364 + t308 * t539 + t351 + (-t308 * t376 + t365) * t306) * t514 - (-t297 * t376 - t308 * t350 - t326 + (-t308 * t364 - t377) * t306) * t515) * t317) * t311, 0, 0; (-t337 * t340 + t341 * t516) * t492 + ((-t326 * t413 - t334 * t416) * t337 + t341 * t477 + (-t340 * t315 - (-t326 * t416 + t334 * t413) * t342 - t341 * t314) * t338) * t320 (-t337 * t356 + t357 * t516) * t492 + ((t328 * t413 - t352 * t416) * t337 + t357 * t477 + (-t356 * t315 - (t328 * t416 + t352 * t413) * t342 - t357 * t314) * t338) * t320 (-t337 * t354 + t355 * t516) * t492 + ((-t331 * t416 + t413 * t456) * t337 + t355 * t477 + (-t354 * t315 - (t331 * t413 + t416 * t456) * t342 - t355 * t314) * t338) * t320 (-t338 * t501 + t505) * t366 * t492 + (-0.2e1 * t366 * t416 * t488 - t323 * t505 + (t323 * t501 + (t314 * t416 + t315 * t413) * t366) * t338) * t320, 0, 0;];
JaD_rot  = t1;
