% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR9_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR9_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_jacobiaD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:20:39
% EndTime: 2019-02-26 22:20:44
% DurationCPUTime: 5.71s
% Computational Cost: add. (29654->261), mult. (86169->486), div. (983->12), fcn. (110623->19), ass. (0->189)
t490 = cos(pkin(6));
t497 = cos(qJ(2));
t593 = sin(qJ(1));
t554 = t593 * t497;
t494 = sin(qJ(2));
t498 = cos(qJ(1));
t570 = t498 * t494;
t477 = t490 * t570 + t554;
t517 = t490 * t554 + t570;
t460 = qJD(1) * t517 + qJD(2) * t477;
t487 = t593 * t494;
t541 = t490 * t487;
t569 = t498 * t497;
t461 = -qJD(1) * t541 - qJD(2) * t487 + (qJD(2) * t490 + qJD(1)) * t569;
t493 = sin(qJ(3));
t590 = sin(pkin(7));
t591 = cos(pkin(13));
t533 = t591 * t590;
t488 = sin(pkin(13));
t552 = t488 * t590;
t594 = cos(qJ(3));
t470 = -t493 * t552 + t594 * t533;
t466 = t470 * qJD(3);
t592 = cos(pkin(7));
t534 = t592 * t591;
t553 = t488 * t592;
t472 = -t493 * t553 + t594 * t534;
t468 = t472 * qJD(3);
t471 = t493 * t533 + t552 * t594;
t473 = t493 * t534 + t553 * t594;
t480 = -t488 * t594 - t493 * t591;
t476 = t480 * qJD(3);
t479 = t493 * t488 - t591 * t594;
t489 = sin(pkin(6));
t599 = -t490 * t569 + t487;
t399 = t460 * t473 + t461 * t479 + t468 * t599 - t477 * t476 - t489 * (qJD(1) * t471 * t593 - t498 * t466);
t572 = t489 * t498;
t442 = -t471 * t572 - t473 * t599 - t477 * t479;
t551 = t489 * t592;
t539 = t498 * t551;
t464 = -t590 * t599 + t539;
t492 = sin(qJ(5));
t496 = cos(qJ(5));
t423 = t442 * t492 + t464 * t496;
t529 = t593 * t551;
t510 = qJD(1) * t529 + t460 * t590;
t388 = qJD(5) * t423 + t399 * t496 - t492 * t510;
t425 = t442 * t496 - t464 * t492;
t386 = qJD(5) * t425 - t399 * t492 - t496 * t510;
t530 = -t473 * t494 - t479 * t497;
t431 = t490 * t466 + (qJD(2) * t530 + t468 * t497 + t476 * t494) * t489;
t531 = t473 * t497 - t479 * t494;
t449 = t490 * t471 + t489 * t531;
t548 = t497 * t590;
t474 = -t489 * t548 + t490 * t592;
t437 = t449 * t496 + t474 * t492;
t540 = t489 * t494 * t590;
t528 = qJD(2) * t540;
t391 = qJD(5) * t437 + t431 * t492 - t496 * t528;
t436 = t449 * t492 - t474 * t496;
t434 = 0.1e1 / t436 ^ 2;
t601 = t391 * t434;
t433 = 0.1e1 / t436;
t448 = t490 * t470 + (t472 * t497 + t480 * t494) * t489;
t575 = t423 * t434;
t597 = -t470 * t572 - t472 * t599 + t477 * t480;
t523 = -t433 * t597 + t448 * t575;
t600 = t492 * t523;
t467 = t471 * qJD(3);
t469 = t473 * qJD(3);
t475 = t479 * qJD(3);
t555 = t489 * t593;
t542 = t470 * t555;
t598 = qJD(1) * t542 - t460 * t472 + t461 * t480 + t467 * t572 + t469 * t599 + t477 * t475;
t413 = atan2(-t423, t436);
t408 = sin(t413);
t409 = cos(t413);
t380 = -t408 * t423 + t409 * t436;
t377 = 0.1e1 / t380;
t506 = t517 * t590 + t529;
t518 = t541 - t569;
t507 = t471 * t555 - t473 * t517 + t479 * t518;
t429 = t492 * t506 + t496 * t507;
t446 = -t472 * t517 - t480 * t518 + t542;
t491 = sin(qJ(6));
t495 = cos(qJ(6));
t407 = t429 * t495 - t446 * t491;
t401 = 0.1e1 / t407;
t378 = 0.1e1 / t380 ^ 2;
t402 = 0.1e1 / t407 ^ 2;
t596 = -0.2e1 * t423;
t428 = t492 * t507 - t496 * t506;
t595 = 0.2e1 * t428;
t422 = t428 ^ 2;
t376 = t422 * t378 + 0.1e1;
t458 = t599 * qJD(1) + t518 * qJD(2);
t459 = qJD(1) * t477 + qJD(2) * t517;
t568 = qJD(1) * t498;
t505 = t459 * t479 - t518 * t476 + t458 * t473 - t517 * t468 + (t466 * t593 + t471 * t568) * t489;
t511 = qJD(1) * t539 - t458 * t590;
t384 = qJD(5) * t429 + t492 * t505 - t496 * t511;
t582 = t378 * t428;
t421 = t423 ^ 2;
t412 = t421 * t434 + 0.1e1;
t410 = 0.1e1 / t412;
t527 = -t386 * t433 + t391 * t575;
t367 = t527 * t410;
t532 = -t408 * t436 - t409 * t423;
t361 = t367 * t532 - t408 * t386 + t409 * t391;
t379 = t377 * t378;
t588 = t361 * t379;
t589 = (t384 * t582 - t422 * t588) / t376 ^ 2;
t385 = -qJD(5) * t428 + t492 * t511 + t496 * t505;
t394 = -t459 * t480 - t518 * t475 + t458 * t472 + t517 * t469 + (-t467 * t593 + t470 * t568) * t489;
t372 = qJD(6) * t407 + t385 * t491 + t394 * t495;
t406 = t429 * t491 + t446 * t495;
t400 = t406 ^ 2;
t383 = t400 * t402 + 0.1e1;
t580 = t402 * t406;
t566 = qJD(6) * t406;
t373 = t385 * t495 - t394 * t491 - t566;
t584 = t373 * t401 * t402;
t586 = (t372 * t580 - t400 * t584) / t383 ^ 2;
t581 = t433 * t601;
t585 = (t386 * t575 - t421 * t581) / t412 ^ 2;
t583 = t378 * t384;
t579 = t406 * t495;
t578 = t408 * t428;
t577 = t409 * t428;
t576 = t423 * t433;
t574 = t446 * t492;
t573 = t446 * t496;
t571 = t491 * t401;
t567 = qJD(5) * t496;
t565 = 0.2e1 * t589;
t564 = -0.2e1 * t586;
t563 = 0.2e1 * t586;
t562 = -0.2e1 * t585;
t561 = t379 * t595;
t560 = t433 * t585;
t559 = t378 * t578;
t558 = t378 * t577;
t557 = t406 * t584;
t550 = t492 * t590;
t549 = t496 * t590;
t547 = -0.2e1 * t377 * t589;
t546 = t378 * t565;
t545 = t361 * t561;
t544 = t581 * t596;
t543 = 0.2e1 * t557;
t535 = qJD(6) * t573 - t505;
t405 = -t425 * t495 + t491 * t597;
t404 = -t425 * t491 - t495 * t597;
t452 = t473 * t518 + t479 * t517;
t440 = t452 * t496 - t518 * t550;
t451 = -t472 * t518 + t480 * t517;
t420 = t440 * t495 + t451 * t491;
t419 = t440 * t491 - t451 * t495;
t526 = t579 * t402 - t571;
t525 = -t425 * t433 + t437 * t575;
t450 = -t477 * t473 + t479 * t599;
t438 = t450 * t492 - t477 * t549;
t457 = t530 * t489;
t453 = t457 * t492 - t496 * t540;
t524 = -t433 * t438 + t453 * t575;
t515 = -t452 * t492 - t518 * t549;
t514 = -t408 + (t409 * t576 + t408) * t410;
t512 = -qJD(5) * t574 + qJD(6) * t507 + t394 * t496;
t430 = -t490 * t467 + (-t469 * t497 + t475 * t494 + (-t472 * t494 + t480 * t497) * qJD(2)) * t489;
t418 = -t458 * t479 + t459 * t473 + t468 * t518 - t476 * t517;
t417 = -t458 * t480 - t459 * t472 + t469 * t518 + t475 * t517;
t416 = t457 * t567 + ((t476 * t497 + (qJD(5) * t590 - t468) * t494) * t492 + (-t492 * t531 - t496 * t548) * qJD(2)) * t489;
t415 = t491 * t507 + t495 * t573;
t414 = t491 * t573 - t495 * t507;
t392 = -qJD(5) * t436 + t431 * t496 + t492 * t528;
t390 = (t460 * t479 - t461 * t473 - t477 * t468 - t476 * t599) * t492 - t461 * t549 + (t450 * t496 + t477 * t550) * qJD(5);
t389 = qJD(5) * t515 + t418 * t496 - t459 * t550;
t381 = 0.1e1 / t383;
t374 = 0.1e1 / t376;
t371 = t524 * t410;
t370 = t410 * t600;
t369 = t525 * t410;
t366 = t514 * t428;
t363 = (-t408 * t597 + t409 * t448) * t492 + t532 * t370;
t362 = t369 * t532 - t408 * t425 + t409 * t437;
t360 = t524 * t562 + (t453 * t544 - t390 * t433 + (t386 * t453 + t391 * t438 + t416 * t423) * t434) * t410;
t358 = t525 * t562 + (t437 * t544 + t388 * t433 + (t386 * t437 + t391 * t425 + t392 * t423) * t434) * t410;
t357 = t562 * t600 + (t523 * t567 + (t448 * t544 - t598 * t433 + (t386 * t448 + t391 * t597 + t423 * t430) * t434) * t492) * t410;
t1 = [t560 * t595 + (-t384 * t433 + t428 * t601) * t410, t360, t357, 0, t358, 0; -t423 * t547 + (-t386 * t377 + (t361 * t423 - t366 * t384) * t378) * t374 + (t366 * t546 + (0.2e1 * t366 * t588 - (-t367 * t410 * t576 + t562) * t559 - (t560 * t596 - t367 + (t367 - t527) * t410) * t558 - t514 * t583) * t374) * t428, -t515 * t547 + ((qJD(5) * t440 + t418 * t492 + t459 * t549) * t377 + t515 * t378 * t361 - ((-t360 * t423 - t371 * t386 + t416 + (-t371 * t436 - t438) * t367) * t409 + (-t360 * t436 - t371 * t391 - t390 + (t371 * t423 - t453) * t367) * t408) * t582) * t374 + (t428 * t546 + (-t583 + t545) * t374) * (t371 * t532 - t408 * t438 + t409 * t453) (t363 * t582 - t377 * t574) * t565 + (-t363 * t583 + (t394 * t492 + t446 * t567) * t377 + (t363 * t561 - t378 * t574) * t361 - (t448 * t567 - t357 * t423 - t370 * t386 + t430 * t492 + (-t370 * t436 - t492 * t597) * t367) * t558 - (-t597 * t567 - t357 * t436 - t370 * t391 - t598 * t492 + (t370 * t423 - t448 * t492) * t367) * t559) * t374, 0 (t362 * t582 - t377 * t429) * t565 + (t362 * t545 + t385 * t377 + (-t429 * t361 - t362 * t384 - (-t358 * t423 - t369 * t386 + t392 + (-t369 * t436 - t425) * t367) * t577 - (-t358 * t436 - t369 * t391 + t388 + (t369 * t423 - t437) * t367) * t578) * t378) * t374, 0; (-t401 * t404 + t405 * t580) * t563 + ((qJD(6) * t405 + t388 * t491 - t495 * t598) * t401 + t405 * t543 + (-t404 * t373 - (-qJD(6) * t404 + t388 * t495 + t491 * t598) * t406 - t405 * t372) * t402) * t381 (-t401 * t419 + t420 * t580) * t563 + ((qJD(6) * t420 + t389 * t491 - t417 * t495) * t401 + t420 * t543 + (-t419 * t373 - (-qJD(6) * t419 + t389 * t495 + t417 * t491) * t406 - t420 * t372) * t402) * t381 (-t401 * t414 + t415 * t580) * t563 + (t415 * t543 + t535 * t401 * t495 + t512 * t571 + (t406 * t491 * t535 - t415 * t372 - t414 * t373 - t512 * t579) * t402) * t381, 0, t526 * t428 * t564 + (t526 * t384 + ((-qJD(6) * t401 - 0.2e1 * t557) * t495 + (t372 * t495 + (t373 - t566) * t491) * t402) * t428) * t381, t564 + 0.2e1 * (t372 * t402 * t381 + (-t381 * t584 - t402 * t586) * t406) * t406;];
JaD_rot  = t1;
