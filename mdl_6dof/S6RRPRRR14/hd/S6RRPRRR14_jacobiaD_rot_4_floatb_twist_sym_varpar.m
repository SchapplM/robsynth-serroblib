% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR14_jacobiaD_rot_4_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_rot_4_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_rot_4_floatb_twist_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobiaD_rot_4_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:38:22
% EndTime: 2018-12-10 18:38:25
% DurationCPUTime: 2.79s
% Computational Cost: add. (29002->190), mult. (32096->339), div. (449->12), fcn. (31350->29), ass. (0->163)
t534 = pkin(6) - qJ(2);
t518 = sin(t534);
t472 = pkin(6) + qJ(2);
t554 = sin(t472);
t526 = t554 / 0.2e1;
t448 = t526 - t518 / 0.2e1;
t479 = sin(qJ(1));
t481 = cos(qJ(2));
t482 = cos(qJ(1));
t428 = t448 * t482 + t479 * t481;
t470 = pkin(7) + pkin(14);
t459 = sin(t470) / 0.2e1;
t471 = pkin(7) - pkin(14);
t468 = sin(t471);
t441 = t459 + t468 / 0.2e1;
t460 = cos(t471) / 0.2e1;
t469 = cos(t470);
t443 = t460 + t469 / 0.2e1;
t473 = sin(pkin(14));
t520 = cos(t534);
t527 = cos(t472) / 0.2e1;
t498 = t520 / 0.2e1 + t527;
t555 = sin(qJ(2));
t492 = t479 * t555 - t482 * t498;
t551 = sin(pkin(6));
t521 = t482 * t551;
t406 = t428 * t473 + t441 * t521 + t492 * t443;
t474 = sin(pkin(8));
t552 = cos(pkin(8));
t475 = sin(pkin(7));
t553 = cos(pkin(7));
t504 = t553 * t551;
t556 = t492 * t475 - t482 * t504;
t559 = -t406 * t474 - t552 * t556;
t507 = t518 / 0.2e1;
t438 = (t507 - t554 / 0.2e1) * qJD(2);
t491 = t479 * t498 + t482 * t555;
t537 = qJD(2) * t481;
t421 = t491 * qJD(1) - t482 * t438 + t479 * t537;
t439 = t498 * qJD(2);
t524 = qJD(2) * t555;
t458 = t479 * t524;
t538 = t479 * t448;
t422 = -qJD(1) * t538 - t458 + (qJD(1) * t481 + t439) * t482;
t515 = qJD(1) * t551;
t512 = t479 * t515;
t388 = t421 * t443 + t422 * t473 - t441 * t512;
t501 = qJD(1) * t504;
t493 = t421 * t475 + t479 * t501;
t377 = -t388 * t474 - t493 * t552;
t431 = -t481 * t482 + t538;
t522 = t479 * t551;
t408 = t431 * t473 + t441 * t522 - t491 * t443;
t442 = t459 - t468 / 0.2e1;
t444 = t460 - t469 / 0.2e1;
t476 = cos(pkin(14));
t409 = -t431 * t476 - t491 * t442 + t444 * t522;
t532 = pkin(8) + qJ(4);
t516 = sin(t532);
t505 = t516 / 0.2e1;
t533 = pkin(8) - qJ(4);
t517 = sin(t533);
t506 = t517 / 0.2e1;
t445 = t505 + t506;
t508 = cos(t532) / 0.2e1;
t519 = cos(t533);
t450 = t519 / 0.2e1 + t508;
t478 = sin(qJ(4));
t490 = -t491 * t475 - t479 * t504;
t371 = -t408 * t450 + t409 * t478 + t490 * t445;
t365 = t371 ^ 2;
t446 = t505 - t517 / 0.2e1;
t449 = t508 - t519 / 0.2e1;
t480 = cos(qJ(4));
t489 = t408 * t446 + t409 * t480 + t490 * t449;
t367 = 0.1e1 / t489 ^ 2;
t558 = t365 * t367;
t447 = t526 + t507;
t451 = t527 - t520 / 0.2e1;
t477 = cos(pkin(6));
t404 = -(t441 * t477 + t443 * t447 + t451 * t473) * t474 + (-t447 * t475 + t477 * t553) * t552;
t402 = 0.1e1 / t404 ^ 2;
t437 = t447 * qJD(2);
t440 = t451 * qJD(2);
t523 = t475 * t552;
t410 = -(-t437 * t473 + t440 * t443) * t474 - t440 * t523;
t541 = t402 * t410;
t419 = t428 * qJD(1) + t479 * t439 + t482 * t524;
t381 = atan2(t559, t404);
t374 = sin(t381);
t375 = cos(t381);
t359 = t374 * t559 + t375 * t404;
t356 = 0.1e1 / t359;
t366 = 0.1e1 / t489;
t401 = 0.1e1 / t404;
t357 = 0.1e1 / t359 ^ 2;
t394 = t408 * t474 + t490 * t552;
t391 = t394 ^ 2;
t354 = t357 * t391 + 0.1e1;
t418 = qJD(1) * t492 - t479 * t438 - t482 * t537;
t511 = t482 * t515;
t386 = t418 * t443 + t419 * t473 + t441 * t511;
t494 = t418 * t475 - t482 * t501;
t376 = t386 * t474 + t494 * t552;
t547 = t357 * t394;
t390 = t559 ^ 2;
t380 = t390 * t402 + 0.1e1;
t378 = 0.1e1 / t380;
t500 = t377 * t401 - t541 * t559;
t350 = t500 * t378;
t503 = -t374 * t404 + t375 * t559;
t346 = t503 * t350 + t374 * t377 + t375 * t410;
t549 = t346 * t356 * t357;
t550 = (t376 * t547 - t391 * t549) / t354 ^ 2;
t540 = t401 * t541;
t548 = (t377 * t402 * t559 - t390 * t540) / t380 ^ 2;
t387 = t418 * t442 - t419 * t476 + t444 * t511;
t433 = t445 * qJD(4);
t435 = t450 * qJD(4);
t536 = qJD(4) * t478;
t361 = t386 * t446 + t387 * t480 + t408 * t435 - t409 * t536 - t490 * t433 + t494 * t449;
t368 = t366 * t367;
t546 = t361 * t368;
t545 = t367 * t371;
t543 = t559 * t401;
t415 = -(t443 * t451 - t447 * t473) * t474 - t451 * t523;
t542 = t559 * t415;
t539 = t431 * t475;
t535 = qJD(4) * t480;
t531 = -0.2e1 * t550;
t530 = -0.2e1 * t549;
t364 = 0.1e1 + t558;
t434 = (t506 - t516 / 0.2e1) * qJD(4);
t436 = t449 * qJD(4);
t360 = -t386 * t450 + t387 * t478 - t408 * t434 + t409 * t535 + t490 * t436 + t494 * t445;
t525 = t360 * t545;
t529 = 0.2e1 * (-t365 * t546 + t525) / t364 ^ 2;
t528 = 0.2e1 * t548;
t514 = -0.2e1 * t401 * t548;
t513 = 0.2e1 * t371 * t546;
t398 = (-t428 * t443 + t492 * t473) * t474 - t428 * t523;
t499 = -t398 * t401 + t402 * t542;
t423 = t431 * qJD(1) - t439 * t482 + t458;
t417 = t431 * t442 - t491 * t476;
t416 = t431 * t443 + t491 * t473;
t411 = -(-t437 * t443 - t440 * t473) * t474 + t437 * t523;
t407 = -t428 * t476 + t492 * t442 + t444 * t521;
t399 = -t416 * t474 - t431 * t523;
t397 = t418 * t476 + t419 * t442;
t396 = -t418 * t473 + t419 * t443;
t389 = t421 * t442 - t422 * t476 - t444 * t512;
t384 = t416 * t446 + t417 * t480 + t449 * t539;
t383 = -t416 * t450 + t417 * t478 + t445 * t539;
t382 = (t421 * t473 + t423 * t443) * t474 + t423 * t523;
t370 = t406 * t446 + t407 * t480 + t449 * t556;
t369 = -t406 * t450 + t407 * t478 + t445 * t556;
t362 = 0.1e1 / t364;
t352 = 0.1e1 / t354;
t351 = t499 * t378;
t349 = (t374 + (t375 * t543 - t374) * t378) * t394;
t347 = -t503 * t351 + t374 * t398 + t375 * t415;
t345 = t499 * t528 + (0.2e1 * t540 * t542 + t382 * t401 + (-t377 * t415 - t398 * t410 - t411 * t559) * t402) * t378;
t1 = [t394 * t514 + (t376 * t401 - t394 * t541) * t378, t345, 0, 0, 0, 0; t559 * t356 * t531 + (t377 * t356 + (-t346 * t559 + t349 * t376) * t357) * t352 + (t349 * t530 * t352 + (t349 * t531 + ((-t376 * t378 + t376 + (-t350 * t378 * t543 + t528) * t394) * t374 + ((t514 * t559 + t350) * t394 + (t376 * t543 + (-t350 + t500) * t394) * t378) * t375) * t352) * t357) * t394, 0.2e1 * (-t347 * t547 - t356 * t399) * t550 + ((-t396 * t474 - t419 * t523) * t356 + (-t399 * t346 + t347 * t376) * t357 + (t347 * t530 + ((t345 * t559 - t351 * t377 + t411 + (t351 * t404 + t398) * t350) * t375 + (-t345 * t404 + t351 * t410 + t382 + (t351 * t559 - t415) * t350) * t374) * t357) * t394) * t352, 0, 0, 0, 0; (-t366 * t369 + t370 * t545) * t529 + ((-t388 * t450 + t389 * t478 - t406 * t434 + t407 * t535 + t436 * t556 + t445 * t493) * t366 + t370 * t513 + (-t369 * t361 - (t388 * t446 + t389 * t480 + t406 * t435 - t407 * t536 - t433 * t556 + t449 * t493) * t371 - t370 * t360) * t367) * t362 (-t366 * t383 + t384 * t545) * t529 + ((t417 * t535 - t396 * t450 + t397 * t478 - t416 * t434 + (t419 * t445 + t431 * t436) * t475) * t366 + t384 * t513 + (-t383 * t361 - (-t417 * t536 + t396 * t446 + t397 * t480 + t416 * t435 + (t419 * t449 - t431 * t433) * t475) * t371 - t384 * t360) * t367) * t362, 0 (-t366 * t489 - t558) * t529 + (0.2e1 * t525 + (-0.2e1 * t368 * t365 - t367 * t489 + t366) * t361) * t362, 0, 0;];
JaD_rot  = t1;
