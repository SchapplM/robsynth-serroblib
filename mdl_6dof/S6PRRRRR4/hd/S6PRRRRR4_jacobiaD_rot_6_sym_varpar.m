% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRRR4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_jacobiaD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:20:37
% EndTime: 2019-02-26 20:20:40
% DurationCPUTime: 3.34s
% Computational Cost: add. (25849->204), mult. (58272->387), div. (1079->12), fcn. (74841->17), ass. (0->168)
t432 = qJ(4) + qJ(5);
t429 = sin(t432);
t430 = cos(t432);
t437 = sin(qJ(3));
t440 = cos(qJ(3));
t433 = sin(pkin(7));
t435 = cos(pkin(7));
t438 = sin(qJ(2));
t441 = cos(qJ(2));
t522 = cos(pkin(13));
t523 = cos(pkin(6));
t475 = t523 * t522;
t521 = sin(pkin(13));
t458 = t438 * t521 - t441 * t475;
t434 = sin(pkin(6));
t486 = t434 * t522;
t450 = -t433 * t486 - t435 * t458;
t457 = -t438 * t475 - t441 * t521;
t404 = t437 * t457 + t440 * t450;
t419 = t458 * qJD(2);
t420 = t457 * qJD(2);
t431 = qJD(4) + qJD(5);
t451 = t433 * t458 - t435 * t486;
t500 = t435 * t437;
t449 = qJD(3) * t404 - t419 * t440 + t420 * t500 + t431 * t451;
t405 = t437 * t450 - t440 * t457;
t472 = t405 * t431 + t420 * t433;
t355 = t429 * t449 + t430 * t472;
t386 = t405 * t429 - t430 * t451;
t384 = t386 ^ 2;
t496 = t438 * t440;
t497 = t437 * t441;
t462 = t435 * t497 + t496;
t487 = t433 * t523;
t415 = t434 * t462 + t437 * t487;
t501 = t433 * t441;
t423 = -t434 * t501 + t435 * t523;
t401 = t415 * t429 - t423 * t430;
t397 = 0.1e1 / t401 ^ 2;
t371 = t384 * t397 + 0.1e1;
t369 = 0.1e1 / t371;
t502 = t433 * t438;
t488 = t434 * t502;
t460 = qJD(2) * t488 - t415 * t431;
t495 = t440 * t441;
t498 = t437 * t438;
t461 = -t435 * t498 + t495;
t464 = t435 * t495 - t498;
t478 = qJD(3) * t487;
t480 = t423 * t431 + t440 * t478 + (qJD(2) * t461 + qJD(3) * t464) * t434;
t373 = t429 * t480 - t430 * t460;
t396 = 0.1e1 / t401;
t511 = t386 * t397;
t338 = (-t355 * t396 + t373 * t511) * t369;
t372 = atan2(-t386, t401);
t367 = sin(t372);
t368 = cos(t372);
t473 = -t367 * t401 - t368 * t386;
t333 = t338 * t473 - t355 * t367 + t368 * t373;
t351 = -t367 * t386 + t368 * t401;
t348 = 0.1e1 / t351;
t349 = 0.1e1 / t351 ^ 2;
t526 = t333 * t348 * t349;
t474 = t523 * t521;
t455 = t438 * t522 + t441 * t474;
t485 = t434 * t521;
t479 = t433 * t485;
t452 = -t435 * t455 + t479;
t456 = t438 * t474 - t441 * t522;
t407 = t437 * t452 - t440 * t456;
t453 = t433 * t455 + t435 * t485;
t389 = t407 * t429 - t430 * t453;
t483 = 0.2e1 * t389 * t526;
t414 = t434 * t464 + t440 * t487;
t465 = -t396 * t404 + t414 * t511;
t525 = t429 * t465;
t512 = t373 * t396 * t397;
t524 = -0.2e1 * (t355 * t511 - t384 * t512) / t371 ^ 2;
t390 = t407 * t430 + t429 * t453;
t439 = cos(qJ(6));
t454 = t455 * t440;
t506 = t456 * t437;
t406 = t435 * t454 - t440 * t479 - t506;
t436 = sin(qJ(6));
t509 = t406 * t436;
t366 = t390 * t439 + t509;
t362 = 0.1e1 / t366;
t363 = 0.1e1 / t366 ^ 2;
t421 = t455 * qJD(2);
t422 = t456 * qJD(2);
t383 = t422 * t500 - t421 * t440 + (t440 * t452 + t506) * qJD(3);
t448 = t431 * t453 + t383;
t471 = t407 * t431 + t422 * t433;
t358 = -t429 * t471 + t430 * t448;
t499 = t435 * t440;
t382 = qJD(3) * t407 - t421 * t437 - t422 * t499;
t508 = t406 * t439;
t365 = t390 * t436 - t508;
t494 = qJD(6) * t365;
t347 = t358 * t439 + t382 * t436 - t494;
t520 = t347 * t362 * t363;
t519 = t349 * t389;
t346 = qJD(6) * t366 + t358 * t436 - t382 * t439;
t361 = t365 ^ 2;
t354 = t361 * t363 + 0.1e1;
t516 = t363 * t365;
t518 = 0.1e1 / t354 ^ 2 * (t346 * t516 - t361 * t520);
t517 = t362 * t436;
t515 = t365 * t439;
t514 = t367 * t389;
t513 = t368 * t389;
t510 = t406 * t429;
t505 = t430 * t431;
t504 = t430 * t433;
t503 = t431 * t433;
t385 = t389 ^ 2;
t345 = t349 * t385 + 0.1e1;
t357 = t429 * t448 + t430 * t471;
t493 = 0.2e1 * (t357 * t519 - t385 * t526) / t345 ^ 2;
t491 = -0.2e1 * t518;
t490 = 0.2e1 * t518;
t489 = t365 * t520;
t482 = 0.2e1 * t489;
t481 = -0.2e1 * t386 * t512;
t409 = -t437 * t455 - t456 * t499;
t477 = -qJD(3) * t409 + t421 * t500 + t422 * t440 - t456 * t503;
t476 = qJD(6) * t406 * t430 + t383;
t410 = t456 * t500 - t454;
t395 = -t429 * t433 * t456 + t410 * t430;
t376 = t395 * t439 + t409 * t436;
t375 = t395 * t436 - t409 * t439;
t470 = t410 * t431 + t421 * t433;
t468 = t363 * t515 - t517;
t388 = t405 * t430 + t429 * t451;
t402 = t415 * t430 + t423 * t429;
t467 = -t388 * t396 + t402 * t511;
t408 = -t440 * t458 + t457 * t500;
t393 = t408 * t429 + t457 * t504;
t418 = t461 * t434;
t411 = t418 * t429 - t430 * t488;
t466 = -t393 * t396 + t411 * t511;
t463 = -t435 * t496 - t497;
t459 = qJD(6) * t407 - t382 * t430 + t431 * t510;
t399 = -t437 * t478 + (qJD(2) * t463 - qJD(3) * t462) * t434;
t394 = t410 * t429 + t456 * t504;
t391 = qJD(3) * t410 - t421 * t499 + t422 * t437;
t380 = -qJD(3) * t405 + t419 * t437 + t420 * t499;
t379 = t418 * t505 + ((qJD(3) * t463 + t431 * t502) * t429 + (-t429 * t462 - t430 * t501) * qJD(2)) * t434;
t378 = t407 * t436 - t430 * t508;
t377 = -t407 * t439 - t430 * t509;
t374 = t429 * t460 + t430 * t480;
t360 = -t429 * t470 + t430 * t477;
t359 = t408 * t505 + t419 * t504 + (t419 * t500 + t420 * t440 + (t437 * t458 + t457 * t499) * qJD(3) - t457 * t503) * t429;
t356 = -t429 * t472 + t430 * t449;
t352 = 0.1e1 / t354;
t343 = 0.1e1 / t345;
t342 = t369 * t525;
t341 = t466 * t369;
t340 = t467 * t369;
t336 = (-t367 * t404 + t368 * t414) * t429 + t473 * t342;
t335 = t341 * t473 - t367 * t393 + t368 * t411;
t334 = t340 * t473 - t367 * t388 + t368 * t402;
t332 = t466 * t524 + (t411 * t481 - t359 * t396 + (t355 * t411 + t373 * t393 + t379 * t386) * t397) * t369;
t330 = t467 * t524 + (t402 * t481 - t356 * t396 + (t355 * t402 + t373 * t388 + t374 * t386) * t397) * t369;
t329 = t524 * t525 + (t465 * t505 + (t414 * t481 - t380 * t396 + (t355 * t414 + t373 * t404 + t386 * t399) * t397) * t429) * t369;
t328 = t468 * t389 * t491 + (t468 * t357 + ((-qJD(6) * t362 - 0.2e1 * t489) * t439 + (t346 * t439 + (t347 - t494) * t436) * t363) * t389) * t352;
t327 = (t334 * t519 - t348 * t390) * t493 + (t334 * t483 + t358 * t348 + (-t390 * t333 - t334 * t357 - (-t330 * t386 - t340 * t355 + t374 + (-t340 * t401 - t388) * t338) * t513 - (-t330 * t401 - t340 * t373 - t356 + (t340 * t386 - t402) * t338) * t514) * t349) * t343;
t1 = [0, t332, t329, t330, t330, 0; 0 (t335 * t519 - t348 * t394) * t493 + ((t429 * t477 + t430 * t470) * t348 + t335 * t483 + (-t394 * t333 - t335 * t357 - (-t332 * t386 - t341 * t355 + t379 + (-t341 * t401 - t393) * t338) * t513 - (-t332 * t401 - t341 * t373 - t359 + (t341 * t386 - t411) * t338) * t514) * t349) * t343 (t336 * t519 + t348 * t510) * t493 + ((-t382 * t429 - t406 * t505) * t348 + t336 * t483 + (-t336 * t357 + t510 * t333 - (t414 * t505 - t329 * t386 - t342 * t355 + t399 * t429 + (-t342 * t401 - t404 * t429) * t338) * t513 - (-t404 * t505 - t329 * t401 - t342 * t373 - t380 * t429 + (t342 * t386 - t414 * t429) * t338) * t514) * t349) * t343, t327, t327, 0; 0 (-t362 * t375 + t376 * t516) * t490 + ((qJD(6) * t376 + t360 * t436 - t391 * t439) * t362 + t376 * t482 + (-t375 * t347 - (-qJD(6) * t375 + t360 * t439 + t391 * t436) * t365 - t376 * t346) * t363) * t352 (-t362 * t377 + t378 * t516) * t490 + (t378 * t482 - t476 * t362 * t439 + t459 * t517 + (-t365 * t436 * t476 - t378 * t346 - t377 * t347 - t459 * t515) * t363) * t352, t328, t328, t491 + 0.2e1 * (t346 * t363 * t352 + (-t352 * t520 - t363 * t518) * t365) * t365;];
JaD_rot  = t1;
