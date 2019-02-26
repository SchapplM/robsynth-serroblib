% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR5
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
% Datum: 2019-02-26 20:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRRR5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_jacobiaD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:21:11
% EndTime: 2019-02-26 20:21:14
% DurationCPUTime: 2.94s
% Computational Cost: add. (15903->204), mult. (45948->386), div. (834->12), fcn. (58926->17), ass. (0->164)
t423 = sin(qJ(3));
t426 = cos(qJ(3));
t419 = sin(pkin(7));
t421 = cos(pkin(7));
t424 = sin(qJ(2));
t427 = cos(qJ(2));
t505 = cos(pkin(13));
t506 = cos(pkin(6));
t456 = t506 * t505;
t504 = sin(pkin(13));
t442 = t504 * t424 - t427 * t456;
t420 = sin(pkin(6));
t470 = t420 * t505;
t434 = -t419 * t470 - t442 * t421;
t441 = -t424 * t456 - t504 * t427;
t388 = t423 * t441 + t434 * t426;
t405 = t442 * qJD(2);
t406 = t441 * qJD(2);
t485 = t421 * t423;
t366 = t388 * qJD(3) - t405 * t426 + t406 * t485;
t389 = t434 * t423 - t426 * t441;
t422 = sin(qJ(4));
t425 = cos(qJ(4));
t435 = t442 * t419 - t421 * t470;
t376 = t389 * t425 + t435 * t422;
t487 = t419 * t425;
t341 = t376 * qJD(4) + t366 * t422 + t406 * t487;
t374 = t389 * t422 - t435 * t425;
t372 = t374 ^ 2;
t481 = t424 * t426;
t482 = t423 * t427;
t445 = t421 * t482 + t481;
t471 = t419 * t506;
t401 = t445 * t420 + t423 * t471;
t486 = t419 * t427;
t409 = -t420 * t486 + t506 * t421;
t392 = t401 * t422 - t409 * t425;
t386 = 0.1e1 / t392 ^ 2;
t357 = t372 * t386 + 0.1e1;
t355 = 0.1e1 / t357;
t480 = t426 * t427;
t483 = t423 * t424;
t444 = -t421 * t483 + t480;
t447 = t421 * t480 - t483;
t458 = qJD(3) * t471;
t384 = t426 * t458 + (t444 * qJD(2) + t447 * qJD(3)) * t420;
t393 = t401 * t425 + t409 * t422;
t488 = t419 * t424;
t472 = t420 * t488;
t460 = qJD(2) * t472;
t359 = t393 * qJD(4) + t384 * t422 - t425 * t460;
t385 = 0.1e1 / t392;
t494 = t374 * t386;
t324 = (-t341 * t385 + t359 * t494) * t355;
t358 = atan2(-t374, t392);
t353 = sin(t358);
t354 = cos(t358);
t454 = -t353 * t392 - t354 * t374;
t319 = t454 * t324 - t353 * t341 + t354 * t359;
t337 = -t353 * t374 + t354 * t392;
t334 = 0.1e1 / t337;
t335 = 0.1e1 / t337 ^ 2;
t509 = t319 * t334 * t335;
t455 = t506 * t504;
t439 = t505 * t424 + t427 * t455;
t469 = t420 * t504;
t459 = t419 * t469;
t436 = -t439 * t421 + t459;
t440 = t424 * t455 - t505 * t427;
t391 = t436 * t423 - t426 * t440;
t437 = t439 * t419 + t421 * t469;
t377 = t391 * t422 - t437 * t425;
t467 = 0.2e1 * t377 * t509;
t400 = t447 * t420 + t426 * t471;
t449 = -t385 * t388 + t400 * t494;
t508 = t422 * t449;
t495 = t359 * t385 * t386;
t507 = -0.2e1 * (t341 * t494 - t372 * t495) / t357 ^ 2;
t378 = t391 * t425 + t437 * t422;
t438 = t439 * t426;
t490 = t440 * t423;
t390 = t421 * t438 - t426 * t459 - t490;
t418 = qJ(5) + qJ(6);
t415 = sin(t418);
t416 = cos(t418);
t352 = t378 * t416 + t390 * t415;
t348 = 0.1e1 / t352;
t349 = 0.1e1 / t352 ^ 2;
t407 = t439 * qJD(2);
t408 = t440 * qJD(2);
t484 = t421 * t426;
t367 = t391 * qJD(3) - t407 * t423 - t408 * t484;
t417 = qJD(5) + qJD(6);
t462 = t378 * t417 - t367;
t368 = t408 * t485 - t407 * t426 + (t436 * t426 + t490) * qJD(3);
t489 = t419 * t422;
t344 = -t377 * qJD(4) + t368 * t425 - t408 * t489;
t464 = t390 * t417 + t344;
t332 = t464 * t415 + t462 * t416;
t351 = t378 * t415 - t390 * t416;
t347 = t351 ^ 2;
t340 = t347 * t349 + 0.1e1;
t499 = t349 * t351;
t333 = -t462 * t415 + t464 * t416;
t502 = t333 * t348 * t349;
t503 = (t332 * t499 - t347 * t502) / t340 ^ 2;
t501 = t335 * t377;
t500 = t348 * t415;
t498 = t351 * t416;
t497 = t353 * t377;
t496 = t354 * t377;
t493 = t390 * t422;
t492 = t390 * t425;
t479 = qJD(4) * t422;
t478 = qJD(4) * t425;
t373 = t377 ^ 2;
t331 = t335 * t373 + 0.1e1;
t343 = t378 * qJD(4) + t368 * t422 + t408 * t487;
t477 = 0.2e1 * (t343 * t501 - t373 * t509) / t331 ^ 2;
t476 = -0.2e1 * t503;
t475 = 0.2e1 * t503;
t473 = t351 * t502;
t466 = 0.2e1 * t473;
t465 = -0.2e1 * t374 * t495;
t395 = -t439 * t423 - t440 * t484;
t371 = -t395 * qJD(3) + t407 * t485 + t408 * t426;
t396 = t440 * t485 - t438;
t448 = -t396 * t422 - t440 * t487;
t463 = t448 * qJD(4) + t371 * t425 + t395 * t417 - t407 * t489;
t381 = t396 * t425 - t440 * t489;
t461 = -t396 * qJD(3) + t381 * t417 + t407 * t484 - t408 * t423;
t457 = t417 * t492 + t368;
t452 = t349 * t498 - t500;
t451 = -t376 * t385 + t393 * t494;
t394 = -t442 * t426 + t441 * t485;
t379 = t394 * t422 + t441 * t487;
t404 = t444 * t420;
t397 = t404 * t422 - t425 * t472;
t450 = -t379 * t385 + t397 * t494;
t446 = -t421 * t481 - t482;
t443 = -t367 * t425 + t390 * t479 + t391 * t417;
t383 = -t423 * t458 + (t446 * qJD(2) - t445 * qJD(3)) * t420;
t369 = t404 * t478 + ((t446 * qJD(3) + qJD(4) * t488) * t422 + (-t445 * t422 - t425 * t486) * qJD(2)) * t420;
t365 = -t389 * qJD(3) + t405 * t423 + t406 * t484;
t364 = t381 * t416 + t395 * t415;
t363 = t381 * t415 - t395 * t416;
t362 = t391 * t415 - t416 * t492;
t361 = -t391 * t416 - t415 * t492;
t360 = -t392 * qJD(4) + t384 * t425 + t422 * t460;
t345 = (t405 * t485 + t406 * t426 + (t442 * t423 + t441 * t484) * qJD(3)) * t422 + t394 * t478 + t405 * t487 - t441 * t419 * t479;
t342 = -t374 * qJD(4) + t366 * t425 - t406 * t489;
t338 = 0.1e1 / t340;
t329 = 0.1e1 / t331;
t328 = t355 * t508;
t327 = t450 * t355;
t326 = t451 * t355;
t322 = (-t353 * t388 + t354 * t400) * t422 + t454 * t328;
t321 = t454 * t327 - t353 * t379 + t354 * t397;
t320 = t454 * t326 - t353 * t376 + t354 * t393;
t318 = t450 * t507 + (t397 * t465 - t345 * t385 + (t341 * t397 + t359 * t379 + t369 * t374) * t386) * t355;
t316 = t451 * t507 + (t393 * t465 - t342 * t385 + (t341 * t393 + t359 * t376 + t360 * t374) * t386) * t355;
t315 = t507 * t508 + (t449 * t478 + (t400 * t465 - t365 * t385 + (t341 * t400 + t359 * t388 + t374 * t383) * t386) * t422) * t355;
t314 = t476 + 0.2e1 * (t332 * t338 * t349 + (-t338 * t502 - t349 * t503) * t351) * t351;
t1 = [0, t318, t315, t316, 0, 0; 0 (t321 * t501 + t334 * t448) * t477 + ((t381 * qJD(4) + t371 * t422 + t407 * t487) * t334 + t321 * t467 + (t448 * t319 - t321 * t343 - (-t318 * t374 - t327 * t341 + t369 + (-t327 * t392 - t379) * t324) * t496 - (-t318 * t392 - t327 * t359 - t345 + (t327 * t374 - t397) * t324) * t497) * t335) * t329 (t322 * t501 + t334 * t493) * t477 + ((-t367 * t422 - t390 * t478) * t334 + t322 * t467 + (-t322 * t343 + t493 * t319 - (t400 * t478 - t315 * t374 - t328 * t341 + t383 * t422 + (-t328 * t392 - t388 * t422) * t324) * t496 - (-t388 * t478 - t315 * t392 - t328 * t359 - t365 * t422 + (t328 * t374 - t400 * t422) * t324) * t497) * t335) * t329 (t320 * t501 - t334 * t378) * t477 + (t320 * t467 + t344 * t334 + (-t378 * t319 - t320 * t343 - (-t316 * t374 - t326 * t341 + t360 + (-t326 * t392 - t376) * t324) * t496 - (-t316 * t392 - t326 * t359 - t342 + (t326 * t374 - t393) * t324) * t497) * t335) * t329, 0, 0; 0 (-t348 * t363 + t364 * t499) * t475 + ((t463 * t415 + t461 * t416) * t348 + t364 * t466 + (-t363 * t333 - (-t461 * t415 + t463 * t416) * t351 - t364 * t332) * t349) * t338 (-t348 * t361 + t362 * t499) * t475 + (t362 * t466 - t457 * t348 * t416 + t443 * t500 + (-t457 * t351 * t415 - t362 * t332 - t361 * t333 - t443 * t498) * t349) * t338, t452 * t377 * t476 + (t452 * t343 + ((-t348 * t417 - 0.2e1 * t473) * t416 + (t332 * t416 + (-t351 * t417 + t333) * t415) * t349) * t377) * t338, t314, t314;];
JaD_rot  = t1;
