% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PPPRRR1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_jacobiaD_rot_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:38:45
% EndTime: 2019-02-26 19:38:48
% DurationCPUTime: 3.30s
% Computational Cost: add. (25805->134), mult. (75366->274), div. (559->12), fcn. (100504->21), ass. (0->137)
t475 = sin(pkin(12));
t480 = cos(pkin(13));
t442 = t475 * t480;
t474 = sin(pkin(13));
t481 = cos(pkin(12));
t447 = t481 * t474;
t484 = cos(pkin(6));
t430 = -t484 * t442 - t447;
t483 = cos(pkin(7));
t427 = t430 * t483;
t441 = t475 * t474;
t448 = t481 * t480;
t431 = -t484 * t441 + t448;
t473 = sin(pkin(14));
t477 = sin(pkin(7));
t443 = t477 * t473;
t478 = sin(pkin(6));
t434 = t478 * t443;
t479 = cos(pkin(14));
t398 = t473 * t427 + t431 * t479 + t475 * t434;
t406 = sin(qJ(4));
t409 = cos(qJ(4));
t444 = t477 * t479;
t435 = t478 * t444;
t422 = -t479 * t427 + t431 * t473 - t475 * t435;
t449 = t483 * t478;
t425 = -t430 * t477 + t475 * t449;
t476 = sin(pkin(8));
t482 = cos(pkin(8));
t486 = t422 * t482 - t425 * t476;
t388 = t398 * t406 + t409 * t486;
t428 = t484 * t448 - t441;
t426 = t428 * t483;
t429 = t484 * t447 + t442;
t397 = t473 * t426 + t429 * t479 - t481 * t434;
t421 = -t479 * t426 + t429 * t473 + t481 * t435;
t424 = -t428 * t477 - t481 * t449;
t417 = -t421 * t482 + t424 * t476;
t387 = t397 * t409 + t417 * t406;
t405 = sin(qJ(5));
t408 = cos(qJ(5));
t416 = t421 * t476 + t424 * t482;
t373 = t387 * t408 + t416 * t405;
t386 = -t397 * t406 + t417 * t409;
t376 = t386 * qJD(4);
t349 = t373 * qJD(5) + t376 * t405;
t371 = t387 * t405 - t416 * t408;
t369 = t371 ^ 2;
t446 = t478 * t480;
t436 = t483 * t446;
t445 = t478 * t474;
t402 = t473 * t436 + t484 * t443 + t479 * t445;
t401 = t479 * t436 + t484 * t444 - t473 * t445;
t403 = -t477 * t446 + t484 * t483;
t433 = t482 * t401 + t476 * t403;
t395 = t402 * t409 + t433 * t406;
t399 = -t401 * t476 + t403 * t482;
t384 = t395 * t405 - t399 * t408;
t382 = 0.1e1 / t384 ^ 2;
t363 = t369 * t382 + 0.1e1;
t361 = 0.1e1 / t363;
t385 = t395 * t408 + t399 * t405;
t394 = -t402 * t406 + t433 * t409;
t390 = t394 * qJD(4);
t367 = t385 * qJD(5) + t390 * t405;
t381 = 0.1e1 / t384;
t463 = t371 * t382;
t333 = (-t349 * t381 + t367 * t463) * t361;
t364 = atan2(-t371, t384);
t359 = sin(t364);
t360 = cos(t364);
t440 = -t359 * t384 - t360 * t371;
t329 = t440 * t333 - t359 * t349 + t360 * t367;
t343 = -t359 * t371 + t360 * t384;
t340 = 0.1e1 / t343;
t341 = 0.1e1 / t343 ^ 2;
t489 = t329 * t340 * t341;
t389 = t398 * t409 - t406 * t486;
t418 = t422 * t476 + t425 * t482;
t374 = t389 * t405 - t418 * t408;
t488 = 0.2e1 * t374 * t489;
t437 = -t381 * t386 + t394 * t463;
t487 = t405 * t437;
t464 = t367 * t381 * t382;
t485 = -0.2e1 * (t349 * t463 - t369 * t464) / t363 ^ 2;
t375 = t389 * t408 + t418 * t405;
t404 = sin(qJ(6));
t407 = cos(qJ(6));
t358 = t375 * t407 + t388 * t404;
t354 = 0.1e1 / t358;
t355 = 0.1e1 / t358 ^ 2;
t472 = t341 * t374;
t378 = t388 * qJD(4);
t352 = -t374 * qJD(5) - t378 * t408;
t379 = t389 * qJD(4);
t357 = t375 * t404 - t388 * t407;
t456 = qJD(6) * t357;
t345 = t352 * t407 + t379 * t404 - t456;
t471 = t345 * t354 * t355;
t344 = t358 * qJD(6) + t352 * t404 - t379 * t407;
t353 = t357 ^ 2;
t348 = t353 * t355 + 0.1e1;
t468 = t355 * t357;
t470 = 0.1e1 / t348 ^ 2 * (t344 * t468 - t353 * t471);
t469 = t354 * t404;
t467 = t357 * t407;
t466 = t359 * t374;
t465 = t360 * t374;
t462 = t388 * t405;
t461 = t388 * t408;
t457 = qJD(5) * t408;
t370 = t374 ^ 2;
t339 = t341 * t370 + 0.1e1;
t351 = t375 * qJD(5) - t378 * t405;
t455 = 0.2e1 * (t351 * t472 - t370 * t489) / t339 ^ 2;
t453 = -0.2e1 * t470;
t452 = t357 * t471;
t451 = -0.2e1 * t371 * t464;
t450 = qJD(6) * t461 - t378;
t439 = t355 * t467 - t469;
t438 = -t373 * t381 + t385 * t463;
t432 = qJD(5) * t462 + qJD(6) * t389 - t379 * t408;
t391 = t395 * qJD(4);
t377 = t387 * qJD(4);
t368 = -t384 * qJD(5) + t390 * t408;
t366 = t389 * t404 - t407 * t461;
t365 = -t389 * t407 - t404 * t461;
t350 = -t371 * qJD(5) + t376 * t408;
t346 = 0.1e1 / t348;
t337 = 0.1e1 / t339;
t335 = t361 * t487;
t334 = t438 * t361;
t331 = (-t359 * t386 + t360 * t394) * t405 + t440 * t335;
t330 = t440 * t334 - t359 * t373 + t360 * t385;
t327 = t438 * t485 + (t385 * t451 - t350 * t381 + (t349 * t385 + t367 * t373 + t368 * t371) * t382) * t361;
t326 = t485 * t487 + (t437 * t457 + (t394 * t451 + t377 * t381 + (t349 * t394 + t367 * t386 - t371 * t391) * t382) * t405) * t361;
t1 = [0, 0, 0, t326, t327, 0; 0, 0, 0 (t331 * t472 + t340 * t462) * t455 + ((-t379 * t405 - t388 * t457) * t340 + t331 * t488 + (-t331 * t351 + t462 * t329 - (t394 * t457 - t326 * t371 - t335 * t349 - t391 * t405 + (-t335 * t384 - t386 * t405) * t333) * t465 - (-t386 * t457 - t326 * t384 - t335 * t367 + t377 * t405 + (t335 * t371 - t394 * t405) * t333) * t466) * t341) * t337 (t330 * t472 - t340 * t375) * t455 + (t330 * t488 + t352 * t340 + (-t375 * t329 - t330 * t351 - (-t327 * t371 - t334 * t349 + t368 + (-t334 * t384 - t373) * t333) * t465 - (-t327 * t384 - t334 * t367 - t350 + (t334 * t371 - t385) * t333) * t466) * t341) * t337, 0; 0, 0, 0, 0.2e1 * (-t354 * t365 + t366 * t468) * t470 + (0.2e1 * t366 * t452 - t450 * t354 * t407 + t432 * t469 + (-t450 * t357 * t404 - t366 * t344 - t365 * t345 - t432 * t467) * t355) * t346, t439 * t374 * t453 + (t439 * t351 + ((-qJD(6) * t354 - 0.2e1 * t452) * t407 + (t344 * t407 + (t345 - t456) * t404) * t355) * t374) * t346, t453 + 0.2e1 * (t344 * t355 * t346 + (-t346 * t471 - t355 * t470) * t357) * t357;];
JaD_rot  = t1;
