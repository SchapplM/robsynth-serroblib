% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRRR6_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR6_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_jacobiaD_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:22:02
% EndTime: 2019-02-26 20:22:06
% DurationCPUTime: 3.90s
% Computational Cost: add. (24618->237), mult. (75913->438), div. (795->12), fcn. (96601->19), ass. (0->184)
t574 = cos(pkin(14));
t575 = cos(pkin(6));
t525 = t575 * t574;
t573 = sin(pkin(14));
t576 = sin(qJ(2));
t577 = cos(qJ(2));
t502 = -t525 * t577 + t573 * t576;
t463 = t502 * qJD(2);
t503 = t525 * t576 + t573 * t577;
t464 = t503 * qJD(2);
t481 = cos(qJ(3));
t478 = sin(qJ(3));
t472 = sin(pkin(7));
t475 = cos(pkin(7));
t473 = sin(pkin(6));
t532 = t473 * t574;
t498 = -t472 * t532 - t475 * t502;
t493 = -t478 * t503 + t481 * t498;
t549 = t475 * t478;
t429 = qJD(3) * t493 - t463 * t481 - t464 * t549;
t471 = sin(pkin(8));
t497 = (t472 * t502 - t475 * t532) * t471;
t585 = qJD(4) * t497 + t429;
t477 = sin(qJ(4));
t480 = cos(qJ(4));
t449 = t478 * t498 + t503 * t481;
t548 = t475 * t481;
t489 = -qJD(3) * t449 + t463 * t478 - t464 * t548;
t491 = qJD(4) * t493;
t584 = t489 * t477 + t480 * t491;
t474 = cos(pkin(8));
t583 = t493 * t474 + t497;
t501 = t475 * t503;
t495 = t478 * t502 - t481 * t501;
t557 = t471 * t472;
t582 = t474 * t495 + t503 * t557;
t535 = t576 * t481;
t538 = t577 * t478;
t508 = -t475 * t535 - t538;
t553 = t472 * t473;
t523 = t576 * t471 * t553;
t581 = t508 * t473 * t474 + t523;
t554 = t471 * t480;
t539 = t472 * t554;
t544 = t449 * qJD(4);
t550 = t474 * t480;
t551 = t474 * t477;
t384 = -t464 * t539 + t477 * t585 + t480 * t544 - t489 * t550 + t491 * t551;
t414 = t449 * t477 - t480 * t583;
t412 = t414 ^ 2;
t511 = t475 * t538 + t535;
t533 = t472 * t575;
t459 = t473 * t511 + t478 * t533;
t536 = t576 * t478;
t537 = t577 * t481;
t510 = t475 * t537 - t536;
t458 = t473 * t510 + t481 * t533;
t527 = t577 * t553;
t521 = t458 * t474 + (t475 * t575 - t527) * t471;
t578 = -t459 * t477 + t480 * t521;
t434 = 0.1e1 / t578 ^ 2;
t402 = t412 * t434 + 0.1e1;
t564 = t414 * t434;
t437 = t459 * t480 + t477 * t521;
t509 = t475 * t536 - t537;
t526 = qJD(3) * t533;
t448 = t481 * t526 + (-qJD(2) * t509 + qJD(3) * t510) * t473;
t447 = -t478 * t526 + (qJD(2) * t508 - qJD(3) * t511) * t473;
t506 = qJD(2) * t523 + t447 * t474;
t404 = qJD(4) * t437 + t448 * t477 - t480 * t506;
t433 = 0.1e1 / t578;
t565 = t404 * t433 * t434;
t580 = -0.2e1 * (t384 * t564 + t412 * t565) / t402 ^ 2;
t524 = t575 * t573;
t469 = -t524 * t576 + t574 * t577;
t504 = t524 * t577 + t574 * t576;
t531 = t473 * t573;
t512 = t472 * t531 - t475 * t504;
t451 = t469 * t481 + t478 * t512;
t514 = t469 * t549 + t481 * t504;
t454 = -t469 * t548 + t478 * t504;
t515 = t454 * t474 + t469 * t557;
t579 = t477 * t514 + t480 * t515;
t403 = atan2(-t414, -t578);
t398 = sin(t403);
t399 = cos(t403);
t380 = -t398 * t414 - t399 * t578;
t377 = 0.1e1 / t380;
t559 = t469 * t478;
t450 = t481 * t512 - t559;
t513 = t472 * t504 + t475 * t531;
t507 = t513 * t471;
t505 = t450 * t474 + t507;
t418 = t451 * t480 + t477 * t505;
t440 = -t450 * t471 + t474 * t513;
t476 = sin(qJ(5));
t479 = cos(qJ(5));
t397 = t418 * t479 + t440 * t476;
t393 = 0.1e1 / t397;
t378 = 0.1e1 / t380 ^ 2;
t394 = 0.1e1 / t397 ^ 2;
t400 = 0.1e1 / t402;
t367 = (t384 * t433 + t404 * t564) * t400;
t522 = t398 * t578 - t399 * t414;
t362 = t367 * t522 - t398 * t384 + t399 * t404;
t572 = t362 * t377 * t378;
t465 = t504 * qJD(2);
t466 = t469 * qJD(2);
t430 = -qJD(3) * t451 + t465 * t478 - t466 * t548;
t431 = -t466 * t549 - qJD(3) * t559 + (qJD(3) * t512 - t465) * t481;
t563 = t451 * t477;
t387 = t431 * t480 + (t430 * t474 + t466 * t557) * t477 + (t480 * t505 - t563) * qJD(4);
t552 = t472 * t474;
t419 = -t430 * t471 + t466 * t552;
t396 = t418 * t476 - t440 * t479;
t545 = qJD(5) * t396;
t376 = t387 * t479 + t419 * t476 - t545;
t571 = t376 * t393 * t394;
t417 = -t450 * t550 - t480 * t507 + t563;
t570 = t378 * t417;
t375 = qJD(5) * t397 + t387 * t476 - t419 * t479;
t392 = t396 ^ 2;
t383 = t392 * t394 + 0.1e1;
t568 = t394 * t396;
t569 = 0.1e1 / t383 ^ 2 * (t375 * t568 - t392 * t571);
t567 = t398 * t417;
t566 = t399 * t417;
t556 = t471 * t476;
t555 = t471 * t479;
t547 = qJD(4) * t477;
t546 = qJD(4) * t480;
t413 = t417 ^ 2;
t374 = t378 * t413 + 0.1e1;
t386 = qJD(4) * t418 - t430 * t550 + t431 * t477 - t466 * t539;
t543 = 0.2e1 * (t386 * t570 - t413 * t572) / t374 ^ 2;
t542 = -0.2e1 * t569;
t541 = 0.2e1 * t569;
t540 = t396 * t571;
t534 = t477 * t544;
t530 = 0.2e1 * t417 * t572;
t529 = 0.2e1 * t540;
t528 = 0.2e1 * t414 * t565;
t423 = t477 * t515 - t480 * t514;
t442 = -t454 * t471 + t469 * t552;
t408 = t423 * t479 + t442 * t476;
t407 = t423 * t476 - t442 * t479;
t520 = -t393 * t476 + t479 * t568;
t416 = t449 * t480 + t477 * t583;
t519 = t416 * t433 + t437 * t564;
t453 = -t478 * t501 - t481 * t502;
t421 = t453 * t477 - t480 * t582;
t461 = t509 * t473;
t443 = -t461 * t477 - t480 * t581;
t518 = t421 * t433 + t443 * t564;
t424 = t449 * t550 + t477 * t493;
t441 = t458 * t477 + t459 * t550;
t517 = t424 * t433 + t441 * t564;
t426 = t450 * t480 - t451 * t551;
t516 = -t426 * t476 + t451 * t555;
t411 = t426 * t479 + t451 * t556;
t425 = t450 * t477 + t451 * t550;
t439 = qJD(3) * t454 + t465 * t549 - t466 * t481;
t438 = qJD(3) * t514 + t465 * t548 + t466 * t478;
t420 = -t438 * t471 - t465 * t552;
t409 = -qJD(2) * t527 * t554 - t461 * t546 + t581 * t547 + ((-qJD(2) * t511 + qJD(3) * t508) * t477 - (-qJD(2) * t510 + qJD(3) * t509) * t550) * t473;
t406 = t448 * t550 + t447 * t477 + (t458 * t480 - t459 * t551) * qJD(4);
t405 = qJD(4) * t578 + t448 * t480 + t506 * t477;
t391 = -qJD(4) * t425 + t430 * t480 - t431 * t551;
t390 = t429 * t550 - t474 * t534 + t584;
t389 = t439 * t480 + (t438 * t474 - t465 * t557) * t477 + t579 * qJD(4);
t388 = (qJD(3) * t495 + t463 * t549 - t464 * t481) * t477 + t453 * t546 - (-qJD(3) * t453 + t463 * t548 + t464 * t478) * t550 + t463 * t539 + t582 * t547;
t385 = t464 * t477 * t557 + t584 * t474 + t480 * t585 - t534;
t381 = 0.1e1 / t383;
t372 = 0.1e1 / t374;
t371 = t518 * t400;
t370 = t517 * t400;
t369 = t519 * t400;
t365 = t371 * t522 - t398 * t421 + t399 * t443;
t364 = t370 * t522 - t398 * t424 + t399 * t441;
t363 = t369 * t522 - t398 * t416 + t399 * t437;
t361 = t518 * t580 + (t443 * t528 + t388 * t433 + (t384 * t443 + t404 * t421 + t409 * t414) * t434) * t400;
t360 = t517 * t580 + (t441 * t528 + t390 * t433 + (t384 * t441 + t404 * t424 + t406 * t414) * t434) * t400;
t358 = t519 * t580 + (t437 * t528 + t385 * t433 + (t384 * t437 + t404 * t416 + t405 * t414) * t434) * t400;
t1 = [0, t361, t360, t358, 0, 0; 0 (t365 * t570 + t377 * t579) * t543 + (t365 * t530 + (t579 * t362 - t365 * t386 - (-t361 * t414 - t371 * t384 + t409 + (t371 * t578 - t421) * t367) * t566 - (t361 * t578 - t371 * t404 - t388 + (t371 * t414 - t443) * t367) * t567) * t378 + (qJD(4) * t423 - t438 * t550 + t439 * t477 + t465 * t539) * t377) * t372 (t364 * t570 - t377 * t425) * t543 + ((qJD(4) * t426 + t430 * t477 + t431 * t550) * t377 + t364 * t530 + (-t425 * t362 - t364 * t386 - (-t360 * t414 - t370 * t384 + t406 + (t370 * t578 - t424) * t367) * t566 - (t360 * t578 - t370 * t404 - t390 + (t370 * t414 - t441) * t367) * t567) * t378) * t372 (t363 * t570 - t377 * t418) * t543 + (t363 * t530 + t387 * t377 + (-t418 * t362 - t363 * t386 - (-t358 * t414 - t369 * t384 + t405 + (t369 * t578 - t416) * t367) * t566 - (t358 * t578 - t369 * t404 - t385 + (t369 * t414 - t437) * t367) * t567) * t378) * t372, 0, 0; 0 (-t393 * t407 + t408 * t568) * t541 + ((qJD(5) * t408 + t389 * t476 - t420 * t479) * t393 + t408 * t529 + (-t407 * t376 - (-qJD(5) * t407 + t389 * t479 + t420 * t476) * t396 - t408 * t375) * t394) * t381 (t393 * t516 + t411 * t568) * t541 + ((qJD(5) * t411 + t391 * t476 - t431 * t555) * t393 + t411 * t529 + (t516 * t376 - (qJD(5) * t516 + t391 * t479 + t431 * t556) * t396 - t411 * t375) * t394) * t381, t520 * t417 * t542 + (t520 * t386 + ((-qJD(5) * t393 - 0.2e1 * t540) * t479 + (t375 * t479 + (t376 - t545) * t476) * t394) * t417) * t381, t542 + 0.2e1 * (t375 * t394 * t381 + (-t381 * t571 - t394 * t569) * t396) * t396, 0;];
JaD_rot  = t1;
