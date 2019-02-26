% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR8
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
% Datum: 2019-02-26 22:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRR8_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR8_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_jacobiaD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:51:24
% EndTime: 2019-02-26 22:51:31
% DurationCPUTime: 6.36s
% Computational Cost: add. (32556->246), mult. (74391->465), div. (1246->12), fcn. (94105->17), ass. (0->196)
t503 = qJ(4) + qJ(5);
t500 = sin(t503);
t620 = cos(pkin(6));
t621 = sin(qJ(2));
t565 = t620 * t621;
t509 = cos(qJ(2));
t622 = sin(qJ(1));
t583 = t622 * t509;
t623 = cos(qJ(1));
t490 = t565 * t623 + t583;
t530 = t583 * t620 + t621 * t623;
t480 = qJD(1) * t530 + qJD(2) * t490;
t499 = t622 * t621;
t544 = t622 * t565;
t566 = t620 * t623;
t481 = -qJD(1) * t544 - qJD(2) * t499 + (qJD(1) * t623 + qJD(2) * t566) * t509;
t502 = qJD(4) + qJD(5);
t506 = sin(qJ(3));
t508 = cos(qJ(3));
t541 = -t509 * t566 + t499;
t504 = sin(pkin(6));
t619 = cos(pkin(7));
t580 = t504 * t619;
t558 = t623 * t580;
t618 = sin(pkin(7));
t523 = t541 * t618 - t558;
t579 = t504 * t618;
t554 = t622 * t579;
t542 = qJD(1) * t554;
t557 = t623 * t579;
t630 = t541 * t619;
t522 = t630 + t557;
t626 = t490 * t506 + t508 * t522;
t569 = t626 * qJD(3) - t502 * t523 - (-t480 * t619 + t542) * t506 - t481 * t508;
t639 = t500 * t569;
t501 = cos(t503);
t598 = t490 * t508;
t470 = t506 * t522 - t598;
t555 = t622 * t580;
t528 = qJD(1) * t555 + t480 * t618;
t567 = t470 * t502 + t528;
t410 = -t500 * t567 + t501 * t569;
t450 = t470 * t500 + t501 * t523;
t451 = t470 * t501 - t500 * t523;
t635 = -t530 * t619 + t554;
t577 = t508 * t619;
t594 = qJD(3) * t506;
t633 = (t506 * t630 - t598) * qJD(3) - t480 * t577 - t481 * t506 + t508 * t542 + t557 * t594;
t576 = t509 * t619;
t539 = t506 * t576 + t508 * t621;
t561 = t620 * t618;
t485 = t504 * t539 + t506 * t561;
t563 = t618 * t621;
t559 = t504 * t563;
t533 = qJD(2) * t559 - t485 * t502;
t575 = t509 * t618;
t489 = -t504 * t575 + t619 * t620;
t564 = t619 * t621;
t537 = -t506 * t564 + t508 * t509;
t538 = -t506 * t621 + t508 * t576;
t553 = qJD(3) * t561;
t568 = t489 * t502 + t508 * t553 + (qJD(2) * t537 + qJD(3) * t538) * t504;
t436 = t500 * t568 - t501 * t533;
t462 = t485 * t500 - t489 * t501;
t458 = 0.1e1 / t462 ^ 2;
t632 = t436 * t458;
t457 = 0.1e1 / t462;
t484 = t504 * t538 + t508 * t561;
t604 = t450 * t458;
t547 = t457 * t626 - t484 * t604;
t631 = t500 * t547;
t529 = -t509 * t623 + t544;
t518 = t541 * qJD(1) + t529 * qJD(2);
t629 = qJD(1) * t557 + t518 * t619;
t628 = t635 * t508;
t429 = atan2(t450, t462);
t416 = sin(t429);
t417 = cos(t429);
t405 = t416 * t450 + t417 * t462;
t402 = 0.1e1 / t405;
t472 = t635 * t506 - t529 * t508;
t520 = t530 * t618 + t555;
t453 = t472 * t501 + t500 * t520;
t507 = cos(qJ(6));
t471 = -t506 * t529 - t628;
t505 = sin(qJ(6));
t602 = t471 * t505;
t427 = t453 * t507 + t602;
t419 = 0.1e1 / t427;
t403 = 0.1e1 / t405 ^ 2;
t420 = 0.1e1 / t427 ^ 2;
t625 = 0.2e1 * t450;
t452 = t472 * t500 - t501 * t520;
t624 = 0.2e1 * t452;
t446 = t452 ^ 2;
t401 = t403 * t446 + 0.1e1;
t516 = qJD(1) * t558 - t518 * t618;
t479 = qJD(1) * t490 + qJD(2) * t530;
t431 = t628 * qJD(3) - t479 * t508 + t629 * t506 + t529 * t594;
t519 = t502 * t520 + t431;
t597 = t501 * t502;
t406 = t472 * t597 + t500 * t519 - t501 * t516;
t610 = t406 * t403;
t445 = t450 ^ 2;
t428 = t445 * t458 + 0.1e1;
t422 = 0.1e1 / t428;
t408 = -t470 * t597 - t501 * t528 - t639;
t551 = -t408 * t457 - t436 * t604;
t392 = t551 * t422;
t560 = -t416 * t462 + t417 * t450;
t386 = t392 * t560 - t408 * t416 + t417 * t436;
t404 = t402 * t403;
t616 = t386 * t404;
t617 = (-t446 * t616 + t452 * t610) / t401 ^ 2;
t407 = t519 * t501 + (-t472 * t502 + t516) * t500;
t430 = t472 * qJD(3) - t479 * t506 - t629 * t508;
t396 = qJD(6) * t427 + t407 * t505 - t430 * t507;
t601 = t471 * t507;
t426 = t453 * t505 - t601;
t418 = t426 ^ 2;
t413 = t418 * t420 + 0.1e1;
t607 = t420 * t426;
t593 = qJD(6) * t426;
t397 = t407 * t507 + t430 * t505 - t593;
t612 = t397 * t419 * t420;
t614 = (t396 * t607 - t418 * t612) / t413 ^ 2;
t606 = t457 * t632;
t613 = (-t408 * t604 - t445 * t606) / t428 ^ 2;
t611 = t403 * t452;
t609 = t416 * t452;
t608 = t417 * t452;
t605 = t450 * t457;
t603 = t471 * t500;
t596 = t505 * t419;
t595 = t507 * t426;
t592 = 0.2e1 * t617;
t591 = -0.2e1 * t614;
t590 = 0.2e1 * t614;
t589 = -0.2e1 * t613;
t588 = t404 * t624;
t587 = t457 * t613;
t586 = t403 * t609;
t585 = t403 * t608;
t584 = t426 * t612;
t582 = t501 * t618;
t581 = t502 * t618;
t578 = t506 * t619;
t574 = -0.2e1 * t402 * t617;
t573 = t403 * t592;
t572 = t386 * t588;
t571 = 0.2e1 * t584;
t570 = t606 * t625;
t562 = qJD(6) * t471 * t501 + t431;
t425 = t451 * t507 - t505 * t626;
t424 = t451 * t505 + t507 * t626;
t477 = -t508 * t530 + t529 * t578;
t456 = -t500 * t529 * t618 + t477 * t501;
t476 = -t506 * t530 - t529 * t577;
t441 = t456 * t507 + t476 * t505;
t440 = t456 * t505 - t476 * t507;
t552 = -qJD(3) * t476 + t479 * t578 + t508 * t518 - t529 * t581;
t550 = t420 * t595 - t596;
t463 = t485 * t501 + t489 * t500;
t549 = t451 * t457 - t463 * t604;
t475 = -t490 * t578 - t508 * t541;
t454 = t475 * t500 - t490 * t582;
t488 = t537 * t504;
t478 = t488 * t500 - t501 * t559;
t548 = -t454 * t457 - t478 * t604;
t545 = t477 * t502 + t479 * t618;
t540 = -t416 + (-t417 * t605 + t416) * t422;
t536 = -t506 * t509 - t508 * t564;
t535 = qJD(6) * t472 - t430 * t501 + t502 * t603;
t460 = -t506 * t553 + (qJD(2) * t536 - qJD(3) * t539) * t504;
t455 = t477 * t500 + t529 * t582;
t444 = t488 * t597 + ((qJD(3) * t536 + t502 * t563) * t500 + (-t500 * t539 - t501 * t575) * qJD(2)) * t504;
t443 = t472 * t505 - t501 * t601;
t442 = -t472 * t507 - t501 * t602;
t438 = qJD(3) * t477 - t479 * t577 + t506 * t518;
t437 = t500 * t533 + t501 * t568;
t415 = t475 * t597 - t481 * t582 + (-t481 * t578 - t480 * t508 + (-t490 * t577 + t506 * t541) * qJD(3) + t490 * t581) * t500;
t414 = -t500 * t545 + t501 * t552;
t411 = 0.1e1 / t413;
t399 = 0.1e1 / t401;
t398 = t422 * t631;
t395 = t548 * t422;
t394 = t549 * t422;
t391 = t540 * t452;
t389 = (t416 * t626 + t417 * t484) * t500 + t560 * t398;
t387 = t394 * t560 + t416 * t451 + t417 * t463;
t385 = t548 * t589 + (t478 * t570 - t415 * t457 + (t408 * t478 + t436 * t454 - t444 * t450) * t458) * t422;
t383 = t549 * t589 + (t463 * t570 + t410 * t457 + (t408 * t463 - t436 * t451 - t437 * t450) * t458) * t422;
t382 = t589 * t631 + (t547 * t597 + (t484 * t570 - t633 * t457 + (t408 * t484 - t436 * t626 - t450 * t460) * t458) * t500) * t422;
t381 = t550 * t452 * t591 + (t550 * t406 + ((-qJD(6) * t419 - 0.2e1 * t584) * t507 + (t396 * t507 + (t397 - t593) * t505) * t420) * t452) * t411;
t380 = (t387 * t611 - t402 * t453) * t592 + (t387 * t572 + t407 * t402 + (-t453 * t386 - t387 * t406 - (t383 * t450 - t394 * t408 + t437 + (-t394 * t462 + t451) * t392) * t608 - (-t383 * t462 - t394 * t436 + t410 + (-t394 * t450 - t463) * t392) * t609) * t403) * t399;
t1 = [t587 * t624 + (-t406 * t457 + t452 * t632) * t422, t385, t382, t383, t383, 0; t450 * t574 + ((t501 * t567 + t639) * t402 + (-t450 * t386 - t391 * t406) * t403) * t399 + (t391 * t573 + (0.2e1 * t391 * t616 - (t392 * t422 * t605 + t589) * t586 - (t587 * t625 - t392 + (t392 - t551) * t422) * t585 - t540 * t610) * t399) * t452, t455 * t574 + ((t500 * t552 + t501 * t545) * t402 - t455 * t403 * t386 - ((t385 * t450 - t395 * t408 + t444 + (-t395 * t462 - t454) * t392) * t417 + (-t385 * t462 - t395 * t436 - t415 + (-t395 * t450 - t478) * t392) * t416) * t611) * t399 + (t452 * t573 + (-t610 + t572) * t399) * (t395 * t560 - t416 * t454 + t417 * t478) (t389 * t611 + t402 * t603) * t592 + (-t389 * t610 + (-t430 * t500 - t471 * t597) * t402 + (t389 * t588 + t403 * t603) * t386 - (t484 * t597 + t382 * t450 - t398 * t408 + t460 * t500 + (-t398 * t462 + t500 * t626) * t392) * t585 - (t626 * t597 - t382 * t462 - t398 * t436 - t633 * t500 + (-t398 * t450 - t484 * t500) * t392) * t586) * t399, t380, t380, 0; (-t419 * t424 + t425 * t607) * t590 + ((qJD(6) * t425 + t410 * t505 - t507 * t633) * t419 + t425 * t571 + (-t424 * t397 - (-qJD(6) * t424 + t410 * t507 + t505 * t633) * t426 - t425 * t396) * t420) * t411 (-t419 * t440 + t441 * t607) * t590 + ((qJD(6) * t441 + t414 * t505 - t438 * t507) * t419 + t441 * t571 + (-t440 * t397 - (-qJD(6) * t440 + t414 * t507 + t438 * t505) * t426 - t441 * t396) * t420) * t411 (-t419 * t442 + t443 * t607) * t590 + (t443 * t571 - t562 * t419 * t507 + t535 * t596 + (-t426 * t505 * t562 - t443 * t396 - t442 * t397 - t535 * t595) * t420) * t411, t381, t381, t591 + 0.2e1 * (t396 * t420 * t411 + (-t411 * t612 - t420 * t614) * t426) * t426;];
JaD_rot  = t1;
