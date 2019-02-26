% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRR9_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR9_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_jacobiaD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:52:09
% EndTime: 2019-02-26 22:52:14
% DurationCPUTime: 5.88s
% Computational Cost: add. (21300->243), mult. (61275->462), div. (1001->12), fcn. (77398->17), ass. (0->193)
t600 = cos(pkin(6));
t601 = sin(qJ(2));
t544 = t600 * t601;
t491 = cos(qJ(2));
t602 = sin(qJ(1));
t565 = t602 * t491;
t603 = cos(qJ(1));
t472 = t544 * t603 + t565;
t510 = t565 * t600 + t601 * t603;
t462 = qJD(1) * t510 + qJD(2) * t472;
t481 = t602 * t601;
t525 = t602 * t544;
t545 = t600 * t603;
t463 = -qJD(1) * t525 - qJD(2) * t481 + (qJD(1) * t603 + qJD(2) * t545) * t491;
t488 = sin(qJ(3));
t490 = cos(qJ(3));
t486 = sin(pkin(6));
t598 = sin(pkin(7));
t563 = t486 * t598;
t533 = t602 * t563;
t523 = qJD(1) * t533;
t599 = cos(pkin(7));
t536 = t603 * t563;
t521 = -t491 * t545 + t481;
t610 = t521 * t599;
t502 = t610 + t536;
t606 = t472 * t488 + t490 * t502;
t411 = qJD(3) * t606 - (-t462 * t599 + t523) * t488 - t463 * t490;
t577 = t472 * t490;
t452 = t488 * t502 - t577;
t487 = sin(qJ(4));
t489 = cos(qJ(4));
t564 = t486 * t599;
t537 = t603 * t564;
t503 = t521 * t598 - t537;
t432 = t452 * t487 + t489 * t503;
t534 = t602 * t564;
t508 = qJD(1) * t534 + t462 * t598;
t394 = qJD(4) * t432 - t411 * t489 + t487 * t508;
t433 = t452 * t489 - t487 * t503;
t618 = qJD(4) * t433 + t411 * t487 + t489 * t508;
t615 = -t510 * t599 + t533;
t559 = t490 * t599;
t576 = qJD(3) * t488;
t408 = (t488 * t610 - t577) * qJD(3) - t462 * t559 - t463 * t488 + t490 * t523 + t536 * t576;
t543 = t599 * t601;
t516 = -t488 * t543 + t490 * t491;
t558 = t491 * t599;
t517 = -t488 * t601 + t490 * t558;
t540 = t600 * t598;
t532 = qJD(3) * t540;
t440 = t490 * t532 + (qJD(2) * t516 + qJD(3) * t517) * t486;
t518 = t488 * t558 + t490 * t601;
t467 = t486 * t518 + t488 * t540;
t557 = t491 * t598;
t471 = -t486 * t557 + t599 * t600;
t446 = t467 * t489 + t471 * t487;
t542 = t598 * t601;
t538 = t486 * t542;
t522 = qJD(2) * t538;
t418 = qJD(4) * t446 + t440 * t487 - t489 * t522;
t445 = t467 * t487 - t471 * t489;
t443 = 0.1e1 / t445 ^ 2;
t612 = t418 * t443;
t442 = 0.1e1 / t445;
t466 = t486 * t517 + t490 * t540;
t582 = t432 * t443;
t527 = t442 * t606 - t466 * t582;
t611 = t487 * t527;
t509 = -t491 * t603 + t525;
t500 = qJD(1) * t521 + t509 * qJD(2);
t609 = qJD(1) * t536 + t500 * t599;
t608 = t615 * t490;
t417 = atan2(t432, t445);
t412 = sin(t417);
t413 = cos(t417);
t387 = t412 * t432 + t413 * t445;
t384 = 0.1e1 / t387;
t454 = t615 * t488 - t509 * t490;
t501 = t510 * t598 + t534;
t435 = t454 * t489 + t487 * t501;
t453 = -t488 * t509 - t608;
t485 = qJ(5) + qJ(6);
t482 = sin(t485);
t483 = cos(t485);
t405 = t435 * t483 + t453 * t482;
t399 = 0.1e1 / t405;
t385 = 0.1e1 / t387 ^ 2;
t400 = 0.1e1 / t405 ^ 2;
t605 = 0.2e1 * t432;
t434 = t454 * t487 - t489 * t501;
t604 = 0.2e1 * t434;
t428 = t434 ^ 2;
t383 = t385 * t428 + 0.1e1;
t461 = qJD(1) * t472 + qJD(2) * t510;
t407 = qJD(3) * t608 - t461 * t490 + t488 * t609 + t509 * t576;
t498 = qJD(1) * t537 - t500 * t598;
t391 = qJD(4) * t435 + t407 * t487 - t489 * t498;
t590 = t385 * t434;
t427 = t432 ^ 2;
t416 = t427 * t443 + 0.1e1;
t414 = 0.1e1 / t416;
t531 = -t418 * t582 + t442 * t618;
t374 = t531 * t414;
t539 = -t412 * t445 + t413 * t432;
t368 = t374 * t539 + t412 * t618 + t413 * t418;
t386 = t384 * t385;
t596 = t368 * t386;
t597 = (t391 * t590 - t428 * t596) / t383 ^ 2;
t406 = t454 * qJD(3) - t461 * t488 - t490 * t609;
t484 = qJD(5) + qJD(6);
t548 = t435 * t484 - t406;
t392 = -qJD(4) * t434 + t407 * t489 + t487 * t498;
t551 = t453 * t484 + t392;
t377 = t482 * t551 + t483 * t548;
t404 = t435 * t482 - t453 * t483;
t398 = t404 ^ 2;
t390 = t398 * t400 + 0.1e1;
t588 = t400 * t404;
t378 = -t482 * t548 + t483 * t551;
t592 = t378 * t399 * t400;
t594 = (t377 * t588 - t398 * t592) / t390 ^ 2;
t584 = t442 * t612;
t593 = (-t427 * t584 + t582 * t618) / t416 ^ 2;
t591 = t385 * t391;
t589 = t399 * t482;
t587 = t404 * t483;
t586 = t412 * t434;
t585 = t413 * t434;
t583 = t432 * t442;
t581 = t453 * t487;
t580 = t453 * t489;
t575 = qJD(4) * t489;
t574 = 0.2e1 * t597;
t573 = -0.2e1 * t594;
t572 = 0.2e1 * t594;
t571 = -0.2e1 * t593;
t570 = t386 * t604;
t569 = t442 * t593;
t568 = t404 * t592;
t567 = t385 * t586;
t566 = t385 * t585;
t562 = t487 * t598;
t561 = t488 * t599;
t560 = t489 * t598;
t556 = -0.2e1 * t384 * t597;
t555 = t385 * t574;
t554 = t368 * t570;
t553 = 0.2e1 * t568;
t552 = t584 * t605;
t550 = -t484 * t606 - t394;
t458 = -t488 * t510 - t509 * t559;
t421 = -qJD(3) * t458 + t461 * t561 + t490 * t500;
t459 = -t490 * t510 + t509 * t561;
t520 = -t459 * t487 - t509 * t560;
t549 = qJD(4) * t520 + t421 * t489 + t458 * t484 - t461 * t562;
t547 = t433 * t484 - t408;
t438 = t459 * t489 - t509 * t562;
t546 = -qJD(3) * t459 + t438 * t484 + t461 * t559 - t488 * t500;
t541 = t484 * t580 + t407;
t530 = t400 * t587 - t589;
t529 = t433 * t442 - t446 * t582;
t457 = -t472 * t561 - t490 * t521;
t436 = t457 * t487 - t472 * t560;
t470 = t516 * t486;
t460 = t470 * t487 - t489 * t538;
t528 = -t436 * t442 - t460 * t582;
t519 = -t412 + (-t413 * t583 + t412) * t414;
t515 = -t488 * t491 - t490 * t543;
t514 = qJD(4) * t581 - t406 * t489 + t454 * t484;
t439 = -t488 * t532 + (qJD(2) * t515 - qJD(3) * t518) * t486;
t426 = t470 * t575 + ((qJD(3) * t515 + qJD(4) * t542) * t487 + (-t487 * t518 - t489 * t557) * qJD(2)) * t486;
t425 = t438 * t483 + t458 * t482;
t424 = t438 * t482 - t458 * t483;
t423 = t454 * t482 - t483 * t580;
t422 = -t454 * t483 - t482 * t580;
t419 = -qJD(4) * t445 + t440 * t489 + t487 * t522;
t403 = t433 * t483 - t482 * t606;
t402 = t433 * t482 + t483 * t606;
t397 = (-t463 * t561 - t462 * t490 + (-t472 * t559 + t488 * t521) * qJD(3)) * t487 + t457 * t575 - t463 * t560 + t472 * qJD(4) * t562;
t388 = 0.1e1 / t390;
t381 = 0.1e1 / t383;
t380 = t414 * t611;
t379 = t528 * t414;
t376 = t529 * t414;
t373 = t519 * t434;
t371 = (t412 * t606 + t413 * t466) * t487 + t539 * t380;
t369 = t376 * t539 + t412 * t433 + t413 * t446;
t367 = t528 * t571 + (t460 * t552 - t397 * t442 + (t418 * t436 - t426 * t432 - t460 * t618) * t443) * t414;
t365 = t529 * t571 + (t446 * t552 - t394 * t442 + (-t418 * t433 - t419 * t432 - t446 * t618) * t443) * t414;
t364 = t571 * t611 + (t527 * t575 + (t466 * t552 - t408 * t442 + (-t418 * t606 - t432 * t439 - t466 * t618) * t443) * t487) * t414;
t363 = t573 + 0.2e1 * (t377 * t388 * t400 + (-t388 * t592 - t400 * t594) * t404) * t404;
t1 = [t569 * t604 + (-t391 * t442 + t434 * t612) * t414, t367, t364, t365, 0, 0; t432 * t556 + (t618 * t384 + (-t368 * t432 - t373 * t391) * t385) * t381 + (t373 * t555 + (0.2e1 * t373 * t596 - (t374 * t414 * t583 + t571) * t567 - (t569 * t605 - t374 + (t374 - t531) * t414) * t566 - t519 * t591) * t381) * t434, -t520 * t556 + ((qJD(4) * t438 + t421 * t487 + t461 * t560) * t384 + t520 * t385 * t368 - ((t367 * t432 + t379 * t618 + t426 + (-t379 * t445 - t436) * t374) * t413 + (-t367 * t445 - t379 * t418 - t397 + (-t379 * t432 - t460) * t374) * t412) * t590) * t381 + (t434 * t555 + (-t591 + t554) * t381) * (t379 * t539 - t412 * t436 + t413 * t460) (t371 * t590 + t384 * t581) * t574 + (-t371 * t591 + (-t406 * t487 - t453 * t575) * t384 + (t371 * t570 + t385 * t581) * t368 - (t466 * t575 + t364 * t432 + t380 * t618 + t439 * t487 + (-t380 * t445 + t487 * t606) * t374) * t566 - (t606 * t575 - t364 * t445 - t380 * t418 - t408 * t487 + (-t380 * t432 - t466 * t487) * t374) * t567) * t381 (t369 * t590 - t384 * t435) * t574 + (t369 * t554 + t392 * t384 + (-t435 * t368 - t369 * t391 - (t365 * t432 + t376 * t618 + t419 + (-t376 * t445 + t433) * t374) * t585 - (-t365 * t445 - t376 * t418 - t394 + (-t376 * t432 - t446) * t374) * t586) * t385) * t381, 0, 0; (-t399 * t402 + t403 * t588) * t572 + ((t482 * t550 + t483 * t547) * t399 + t403 * t553 + (-t402 * t378 - (-t482 * t547 + t483 * t550) * t404 - t403 * t377) * t400) * t388 (-t399 * t424 + t425 * t588) * t572 + ((t482 * t549 + t483 * t546) * t399 + t425 * t553 + (-t424 * t378 - (-t482 * t546 + t483 * t549) * t404 - t425 * t377) * t400) * t388 (-t399 * t422 + t423 * t588) * t572 + (t423 * t553 - t541 * t399 * t483 + t514 * t589 + (-t404 * t482 * t541 - t423 * t377 - t422 * t378 - t514 * t587) * t400) * t388, t530 * t434 * t573 + (t530 * t391 + ((-t399 * t484 - 0.2e1 * t568) * t483 + (t377 * t483 + (-t404 * t484 + t378) * t482) * t400) * t434) * t388, t363, t363;];
JaD_rot  = t1;
