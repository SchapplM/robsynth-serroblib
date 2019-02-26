% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR3
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
% Datum: 2019-02-26 20:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRR3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_jacobiaD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:15
% EndTime: 2019-02-26 20:05:18
% DurationCPUTime: 3.45s
% Computational Cost: add. (22420->219), mult. (65370->407), div. (816->12), fcn. (84818->19), ass. (0->160)
t431 = sin(pkin(7));
t430 = sin(pkin(13));
t435 = sin(qJ(3));
t509 = cos(pkin(13));
t513 = cos(qJ(3));
t456 = -t435 * t430 + t513 * t509;
t409 = t456 * t431;
t405 = qJD(3) * t409;
t511 = cos(pkin(7));
t468 = t511 * t509;
t481 = t430 * t511;
t411 = -t435 * t481 + t513 * t468;
t407 = t411 * qJD(3);
t436 = sin(qJ(2));
t439 = cos(qJ(2));
t510 = cos(pkin(12));
t512 = cos(pkin(6));
t470 = t512 * t510;
t508 = sin(pkin(12));
t452 = t508 * t436 - t439 * t470;
t413 = t452 * qJD(2);
t451 = -t436 * t470 - t439 * t508;
t414 = t451 * qJD(2);
t457 = t430 * t513 + t435 * t509;
t422 = t457 * qJD(3);
t448 = -t435 * t468 - t481 * t513;
t432 = sin(pkin(6));
t480 = t432 * t510;
t365 = -t405 * t480 - t407 * t452 - t413 * t456 - t414 * t448 + t422 * t451;
t410 = t457 * t431;
t390 = -t410 * t480 + t448 * t452 - t451 * t456;
t434 = sin(qJ(5));
t438 = cos(qJ(5));
t446 = t431 * t452 - t480 * t511;
t377 = t390 * t438 + t434 * t446;
t493 = t431 * t438;
t341 = qJD(5) * t377 + t365 * t434 + t414 * t493;
t465 = t436 * t448 + t439 * t456;
t381 = t512 * t405 + (qJD(2) * t465 + t407 * t439 - t422 * t436) * t432;
t466 = t436 * t456 - t439 * t448;
t397 = t410 * t512 + t432 * t466;
t492 = t431 * t439;
t417 = -t432 * t492 + t511 * t512;
t395 = t397 * t438 + t417 * t434;
t482 = t431 * t432 * t436;
t474 = qJD(2) * t482;
t347 = qJD(5) * t395 + t381 * t434 - t438 * t474;
t375 = t390 * t434 - t438 * t446;
t373 = t375 ^ 2;
t394 = t397 * t434 - t417 * t438;
t387 = 0.1e1 / t394 ^ 2;
t359 = t373 * t387 + 0.1e1;
t357 = 0.1e1 / t359;
t386 = 0.1e1 / t394;
t497 = t375 * t387;
t324 = (-t341 * t386 + t347 * t497) * t357;
t360 = atan2(-t375, t394);
t355 = sin(t360);
t356 = cos(t360);
t467 = -t355 * t394 - t356 * t375;
t319 = t324 * t467 - t355 * t341 + t356 * t347;
t337 = -t355 * t375 + t356 * t394;
t334 = 0.1e1 / t337;
t335 = 0.1e1 / t337 ^ 2;
t516 = t319 * t334 * t335;
t469 = t512 * t508;
t449 = t436 * t510 + t439 * t469;
t450 = t436 * t469 - t439 * t510;
t479 = t432 * t508;
t453 = t410 * t479 + t448 * t449 - t450 * t456;
t455 = t431 * t449 + t479 * t511;
t378 = t434 * t453 - t438 * t455;
t475 = 0.2e1 * t378 * t516;
t389 = -t409 * t480 - t411 * t452 + t451 * t457;
t396 = t512 * t409 + (t411 * t439 - t436 * t457) * t432;
t460 = -t386 * t389 + t396 * t497;
t515 = t434 * t460;
t503 = t347 * t386 * t387;
t514 = -0.2e1 * (t341 * t497 - t373 * t503) / t359 ^ 2;
t379 = t434 * t455 + t438 * t453;
t392 = t409 * t479 - t411 * t449 + t450 * t457;
t433 = sin(qJ(6));
t437 = cos(qJ(6));
t354 = t379 * t437 - t392 * t433;
t350 = 0.1e1 / t354;
t351 = 0.1e1 / t354 ^ 2;
t416 = t450 * qJD(2);
t415 = t449 * qJD(2);
t447 = t405 * t479 - t407 * t449 - t415 * t456 - t416 * t448 + t422 * t450;
t494 = t431 * t434;
t344 = -qJD(5) * t378 - t416 * t494 + t438 * t447;
t406 = t431 * t422;
t408 = t448 * qJD(3);
t421 = t456 * qJD(3);
t367 = -t406 * t479 - t408 * t449 + t416 * t411 + t415 * t457 + t421 * t450;
t332 = qJD(6) * t354 + t344 * t433 + t367 * t437;
t353 = t379 * t433 + t392 * t437;
t349 = t353 ^ 2;
t340 = t349 * t351 + 0.1e1;
t501 = t351 * t353;
t488 = qJD(6) * t353;
t333 = t344 * t437 - t367 * t433 - t488;
t506 = t333 * t350 * t351;
t507 = (t332 * t501 - t349 * t506) / t340 ^ 2;
t505 = t335 * t378;
t343 = qJD(5) * t379 + t416 * t493 + t434 * t447;
t504 = t343 * t335;
t502 = t350 * t433;
t500 = t353 * t437;
t499 = t355 * t378;
t498 = t356 * t378;
t496 = t392 * t434;
t495 = t392 * t438;
t489 = qJD(5) * t438;
t374 = t378 ^ 2;
t331 = t335 * t374 + 0.1e1;
t487 = 0.2e1 * (-t374 * t516 + t378 * t504) / t331 ^ 2;
t486 = -0.2e1 * t507;
t485 = 0.2e1 * t507;
t483 = t353 * t506;
t477 = -0.2e1 * t375 * t503;
t476 = 0.2e1 * t483;
t471 = qJD(6) * t495 - t447;
t400 = -t448 * t450 - t449 * t456;
t385 = t400 * t438 - t450 * t494;
t399 = -t411 * t450 - t449 * t457;
t370 = t385 * t437 + t399 * t433;
t369 = t385 * t433 - t399 * t437;
t463 = t500 * t351 - t502;
t462 = -t377 * t386 + t395 * t497;
t398 = -t448 * t451 - t452 * t456;
t383 = t398 * t434 + t451 * t493;
t402 = t465 * t432;
t401 = t402 * t434 - t438 * t482;
t461 = -t383 * t386 + t401 * t497;
t459 = -t400 * t434 - t450 * t493;
t454 = -qJD(5) * t496 + qJD(6) * t453 + t367 * t438;
t380 = -t512 * t406 + (t408 * t439 - t421 * t436 + (-t411 * t436 - t439 * t457) * qJD(2)) * t432;
t372 = t407 * t450 - t415 * t448 + t416 * t456 + t422 * t449;
t371 = -t408 * t450 - t411 * t415 + t416 * t457 - t421 * t449;
t364 = t406 * t480 - t408 * t452 + t414 * t411 + t413 * t457 + t421 * t451;
t363 = t402 * t489 + ((-t422 * t439 + (qJD(5) * t431 - t407) * t436) * t434 + (-t434 * t466 - t438 * t492) * qJD(2)) * t432;
t362 = t433 * t453 + t437 * t495;
t361 = t433 * t495 - t437 * t453;
t348 = -qJD(5) * t394 + t381 * t438 + t434 * t474;
t346 = qJD(5) * t459 + t372 * t438 - t415 * t494;
t345 = (t407 * t451 - t413 * t448 + t414 * t456 + t422 * t452) * t434 + t413 * t493 + (t398 * t438 - t451 * t494) * qJD(5);
t342 = -qJD(5) * t375 + t365 * t438 - t414 * t494;
t338 = 0.1e1 / t340;
t329 = 0.1e1 / t331;
t328 = t461 * t357;
t327 = t357 * t515;
t326 = t462 * t357;
t322 = t328 * t467 - t355 * t383 + t356 * t401;
t321 = (-t355 * t389 + t356 * t396) * t434 + t467 * t327;
t320 = t326 * t467 - t355 * t377 + t356 * t395;
t318 = t461 * t514 + (t401 * t477 - t345 * t386 + (t341 * t401 + t347 * t383 + t363 * t375) * t387) * t357;
t316 = t462 * t514 + (t395 * t477 - t342 * t386 + (t341 * t395 + t347 * t377 + t348 * t375) * t387) * t357;
t315 = t514 * t515 + (t460 * t489 + (t396 * t477 - t364 * t386 + (t341 * t396 + t347 * t389 + t375 * t380) * t387) * t434) * t357;
t1 = [0, t318, t315, 0, t316, 0; 0 (t322 * t505 + t334 * t459) * t487 + ((qJD(5) * t385 + t372 * t434 + t415 * t493) * t334 + t322 * t475 + (t459 * t319 - t322 * t343 - (-t318 * t375 - t328 * t341 + t363 + (-t328 * t394 - t383) * t324) * t498 - (-t318 * t394 - t328 * t347 - t345 + (t328 * t375 - t401) * t324) * t499) * t335) * t329 (t321 * t505 - t334 * t496) * t487 + ((t367 * t434 + t392 * t489) * t334 + (-t504 + t475) * t321 + (-t496 * t319 - (t396 * t489 - t315 * t375 - t327 * t341 + t380 * t434 + (-t327 * t394 - t389 * t434) * t324) * t498 - (-t389 * t489 - t315 * t394 - t327 * t347 - t364 * t434 + (t327 * t375 - t396 * t434) * t324) * t499) * t335) * t329, 0 (t320 * t505 - t334 * t379) * t487 + (t320 * t475 + t344 * t334 + (-t379 * t319 - t320 * t343 - (-t316 * t375 - t326 * t341 + t348 + (-t326 * t394 - t377) * t324) * t498 - (-t316 * t394 - t326 * t347 - t342 + (t326 * t375 - t395) * t324) * t499) * t335) * t329, 0; 0 (-t350 * t369 + t370 * t501) * t485 + ((qJD(6) * t370 + t346 * t433 - t371 * t437) * t350 + t370 * t476 + (-t369 * t333 - (-qJD(6) * t369 + t346 * t437 + t371 * t433) * t353 - t370 * t332) * t351) * t338 (-t350 * t361 + t362 * t501) * t485 + (t362 * t476 + t471 * t350 * t437 + t454 * t502 + (t353 * t433 * t471 - t362 * t332 - t361 * t333 - t454 * t500) * t351) * t338, 0, t463 * t378 * t486 + (t463 * t343 + ((-qJD(6) * t350 - 0.2e1 * t483) * t437 + (t332 * t437 + (t333 - t488) * t433) * t351) * t378) * t338, t486 + 0.2e1 * (t332 * t351 * t338 + (-t338 * t506 - t351 * t507) * t353) * t353;];
JaD_rot  = t1;
