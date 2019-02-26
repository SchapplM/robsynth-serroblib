% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP12
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRP12_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP12_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:46:24
% EndTime: 2019-02-26 22:46:33
% DurationCPUTime: 8.10s
% Computational Cost: add. (38608->306), mult. (114623->564), div. (1208->12), fcn. (145320->17), ass. (0->230)
t522 = sin(qJ(3));
t519 = cos(pkin(7));
t639 = sin(qJ(2));
t640 = sin(qJ(1));
t516 = t640 * t639;
t525 = cos(qJ(2));
t638 = cos(pkin(6));
t642 = cos(qJ(1));
t575 = t638 * t642;
t558 = -t525 * t575 + t516;
t550 = t558 * t519;
t517 = sin(pkin(7));
t518 = sin(pkin(6));
t598 = t518 * t642;
t583 = t517 * t598;
t574 = t638 * t639;
t594 = t640 * t525;
t546 = t642 * t574 + t594;
t641 = cos(qJ(3));
t647 = t546 * t641;
t481 = (t550 + t583) * t522 - t647;
t499 = t558 * t517 - t519 * t598;
t521 = sin(qJ(4));
t524 = cos(qJ(4));
t464 = t481 * t524 - t499 * t521;
t520 = sin(qJ(5));
t661 = t464 * t520;
t523 = cos(qJ(5));
t660 = t464 * t523;
t461 = t481 * t521 + t499 * t524;
t545 = t638 * t594 + t642 * t639;
t494 = t545 * qJD(1) + t546 * qJD(2);
t597 = t518 * t640;
t577 = qJD(1) * t597;
t485 = t494 * t517 + t519 * t577;
t659 = t461 * qJD(4) + t485 * t521;
t658 = t464 * qJD(4) + t485 * t524;
t542 = t546 * t522;
t560 = t641 * t583;
t652 = -t560 - t542;
t559 = t640 * t574;
t544 = -t642 * t525 + t559;
t580 = t517 * t597;
t482 = -t544 * t641 + (-t545 * t519 + t580) * t522;
t500 = t545 * t517 + t519 * t597;
t466 = t482 * t524 + t500 * t521;
t510 = t641 * t580;
t540 = t545 * t641;
t573 = t519 * t540 - t522 * t544 - t510;
t443 = t466 * t520 - t573 * t523;
t651 = -0.2e1 * t443;
t492 = t558 * qJD(1) + t544 * qJD(2);
t493 = t546 * qJD(1) + t545 * qJD(2);
t592 = t642 * qJD(1);
t578 = t518 * t592;
t429 = -t493 * t641 + (t492 * t519 + t517 * t578) * t522 - t573 * qJD(3);
t465 = -t482 * t521 + t500 * t524;
t483 = -t492 * t517 + t519 * t578;
t408 = t465 * qJD(4) + t429 * t524 + t483 * t521;
t596 = t519 * t641;
t428 = -qJD(1) * t560 + t482 * qJD(3) - t492 * t596 - t493 * t522;
t402 = -t443 * qJD(5) + t408 * t523 + t428 * t520;
t444 = t466 * t523 + t573 * t520;
t436 = 0.1e1 / t444 ^ 2;
t650 = t402 * t436;
t593 = t639 * t522;
t595 = t641 * t525;
t553 = -t519 * t593 + t595;
t555 = t519 * t595 - t593;
t591 = t517 * t638;
t570 = t641 * t591;
t473 = qJD(3) * t570 + (t553 * qJD(2) + t555 * qJD(3)) * t518;
t579 = t641 * t639;
t614 = t522 * t525;
t556 = t519 * t614 + t579;
t576 = t522 * t591;
t497 = t556 * t518 + t576;
t617 = t517 * t525;
t505 = -t518 * t617 + t638 * t519;
t477 = -t497 * t521 + t505 * t524;
t599 = t517 * t639;
t584 = t518 * t599;
t571 = qJD(2) * t584;
t446 = t477 * qJD(4) + t473 * t524 + t521 * t571;
t478 = t497 * t524 + t505 * t521;
t496 = -t555 * t518 - t570;
t459 = t478 * t523 + t496 * t520;
t554 = t519 * t579 + t614;
t472 = qJD(3) * t576 + (t554 * qJD(2) + t556 * qJD(3)) * t518;
t412 = t459 * qJD(5) + t446 * t520 - t472 * t523;
t458 = t478 * t520 - t496 * t523;
t456 = 0.1e1 / t458 ^ 2;
t649 = t412 * t456;
t455 = 0.1e1 / t458;
t547 = t558 * t641;
t646 = -t519 * t547 + t652;
t438 = t523 * t646 - t661;
t621 = t438 * t456;
t562 = -t455 * t461 + t477 * t621;
t648 = t520 * t562;
t495 = -qJD(1) * t559 - qJD(2) * t516 + (qJD(2) * t575 + t592) * t525;
t645 = (-t494 * t519 + t517 * t577) * t522 + t495 * t641;
t548 = qJD(1) * t510 - t494 * t596 + (qJD(3) * t583 - t495) * t522;
t551 = t558 * t522;
t532 = (-t519 * t551 + t647) * qJD(3) - t548;
t419 = atan2(-t438, t458);
t414 = sin(t419);
t415 = cos(t419);
t400 = -t414 * t438 + t415 * t458;
t397 = 0.1e1 / t400;
t435 = 0.1e1 / t444;
t398 = 0.1e1 / t400 ^ 2;
t644 = -0.2e1 * t438;
t643 = 0.2e1 * t443;
t434 = t443 ^ 2;
t396 = t398 * t434 + 0.1e1;
t401 = t444 * qJD(5) + t408 * t520 - t428 * t523;
t631 = t398 * t443;
t433 = t438 ^ 2;
t418 = t433 * t456 + 0.1e1;
t416 = 0.1e1 / t418;
t430 = t646 * qJD(3) + t645;
t410 = t430 * t524 + t659;
t440 = -t520 * t646 - t660;
t403 = t440 * qJD(5) + t410 * t520 - t532 * t523;
t566 = -t403 * t455 + t412 * t621;
t387 = t566 * t416;
t572 = -t414 * t458 - t415 * t438;
t381 = t572 * t387 - t403 * t414 + t412 * t415;
t399 = t397 * t398;
t636 = t381 * t399;
t637 = (t401 * t631 - t434 * t636) / t396 ^ 2;
t407 = -t466 * qJD(4) - t429 * t521 + t483 * t524;
t460 = t465 ^ 2;
t624 = t436 * t460;
t422 = 0.1e1 + t624;
t630 = t435 * t650;
t602 = t460 * t630;
t623 = t436 * t465;
t634 = (t407 * t623 - t602) / t422 ^ 2;
t629 = t455 * t649;
t633 = (t403 * t621 - t433 * t629) / t418 ^ 2;
t632 = t398 * t401;
t628 = t414 * t443;
t627 = t415 * t443;
t622 = t438 * t455;
t620 = t465 * t520;
t619 = t517 * t521;
t618 = t517 * t524;
t616 = t519 * t522;
t615 = t520 * t524;
t613 = qJD(4) * t521;
t612 = qJD(4) * t524;
t611 = qJD(5) * t520;
t610 = qJD(5) * t523;
t609 = 0.2e1 * t637;
t608 = 0.2e1 * t634;
t607 = -0.2e1 * t633;
t606 = t399 * t643;
t605 = t455 * t633;
t604 = t398 * t628;
t603 = t398 * t627;
t601 = t465 * t630;
t590 = -0.2e1 * t397 * t637;
t589 = t398 * t609;
t588 = t381 * t606;
t587 = 0.2e1 * t601;
t586 = t629 * t644;
t585 = t623 * t634;
t480 = -t641 * t550 + t652;
t442 = t480 * t520 + t660;
t441 = -t480 * t523 + t661;
t490 = t544 * t616 - t540;
t470 = t490 * t524 - t544 * t619;
t489 = -t545 * t522 - t544 * t596;
t454 = t470 * t523 + t489 * t520;
t453 = t470 * t520 - t489 * t523;
t568 = t521 * t573;
t567 = t524 * t573;
t565 = -t440 * t455 + t459 * t621;
t535 = t524 * t646;
t449 = t481 * t523 + t520 * t535;
t467 = -t496 * t615 - t497 * t523;
t564 = -t449 * t455 + t467 * t621;
t488 = -t519 * t542 - t547;
t543 = t517 * t546;
t468 = t488 * t524 + t521 * t543;
t534 = t519 * t647 - t551;
t452 = t468 * t520 - t534 * t523;
t504 = t553 * t518;
t491 = t504 * t524 + t521 * t584;
t503 = t554 * t518;
t471 = t491 * t520 - t503 * t523;
t563 = -t452 * t455 + t471 * t621;
t469 = -t490 * t521 - t544 * t618;
t561 = qJD(4) * t573;
t557 = -t414 + (t415 * t622 + t414) * t416;
t549 = -qJD(5) * t567 - t429;
t541 = t482 * qJD(5) - t428 * t524 + t521 * t561;
t451 = t482 * t520 - t523 * t567;
t450 = -t482 * t523 - t520 * t567;
t448 = -t489 * qJD(3) + t492 * t641 + t493 * t616;
t447 = t490 * qJD(3) + t492 * t522 - t493 * t596;
t445 = -t478 * qJD(4) - t473 * t521 + t524 * t571;
t432 = -t480 * qJD(3) - t645;
t431 = (t522 * t550 - t647) * qJD(3) + t548;
t425 = -t504 * t520 * t613 + (t491 * t523 + t503 * t520) * qJD(5) + (t520 * t599 * t612 + (-t553 * t523 - t554 * t615) * qJD(3) + ((t521 * t617 - t556 * t524) * t520 - t555 * t523) * qJD(2)) * t518;
t424 = (-qJD(5) * t496 * t524 - t473) * t523 + (qJD(5) * t497 - t472 * t524 + t496 * t613) * t520;
t423 = t469 * qJD(4) + t448 * t524 - t493 * t619;
t420 = 0.1e1 / t422;
t413 = -t458 * qJD(5) + t446 * t523 + t472 * t520;
t411 = t432 * t524 - t659;
t409 = -t430 * t521 + t658;
t406 = ((-t534 * qJD(3) - t494 * t641 - t495 * t616) * t524 - t488 * t613 + t495 * t619 + t543 * t612) * t520 + t468 * t610 - (t488 * qJD(3) - t494 * t522 + t495 * t596) * t523 + t534 * t611;
t405 = (qJD(5) * t535 - t430) * t523 + (-qJD(5) * t481 - t532 * t524 - t613 * t646) * t520;
t404 = -t438 * qJD(5) + t410 * t523 + t532 * t520;
t394 = 0.1e1 / t396;
t393 = t416 * t648;
t392 = t563 * t416;
t391 = t564 * t416;
t390 = t565 * t416;
t386 = t557 * t443;
t385 = (-t414 * t461 + t415 * t477) * t520 + t572 * t393;
t384 = t572 * t392 - t414 * t452 + t415 * t471;
t382 = t572 * t390 - t414 * t440 + t415 * t459;
t380 = t563 * t607 + (t471 * t586 - t406 * t455 + (t403 * t471 + t412 * t452 + t425 * t438) * t456) * t416;
t379 = t564 * t607 + (t467 * t586 - t405 * t455 + (t403 * t467 + t412 * t449 + t424 * t438) * t456) * t416;
t377 = t565 * t607 + (t459 * t586 - t404 * t455 + (t403 * t459 + t412 * t440 + t413 * t438) * t456) * t416;
t376 = t607 * t648 + (t562 * t610 + (t477 * t586 - t409 * t455 + (t403 * t477 + t412 * t461 + t438 * t445) * t456) * t520) * t416;
t1 = [t605 * t643 + (-t401 * t455 + t443 * t649) * t416, t380, t379, t376, t377, 0; t441 * t590 + ((t442 * qJD(5) + t411 * t520 - t431 * t523) * t397 + (-t441 * t381 - t386 * t401) * t398) * t394 + (t386 * t589 + (0.2e1 * t386 * t636 - (-t387 * t416 * t622 + t607) * t604 - (t605 * t644 - t387 + (t387 - t566) * t416) * t603 - t557 * t632) * t394) * t443 (t384 * t631 - t397 * t453) * t609 + ((t454 * qJD(5) + t423 * t520 - t447 * t523) * t397 + t384 * t588 + (-t453 * t381 - t384 * t401 - (-t380 * t438 - t392 * t403 + t425 + (-t392 * t458 - t452) * t387) * t627 - (-t380 * t458 - t392 * t412 - t406 + (t392 * t438 - t471) * t387) * t628) * t398) * t394, t450 * t590 + ((t541 * t520 + t549 * t523) * t397 - t450 * t398 * t381 - ((-t379 * t438 - t391 * t403 + t424 + (-t391 * t458 - t449) * t387) * t415 + (-t379 * t458 - t391 * t412 - t405 + (t391 * t438 - t467) * t387) * t414) * t631) * t394 + (t443 * t589 + (-t632 + t588) * t394) * (t572 * t391 - t414 * t449 + t415 * t467) (t385 * t631 - t397 * t620) * t609 + (-t385 * t632 + (t407 * t520 + t465 * t610) * t397 + (t385 * t606 - t398 * t620) * t381 - (t477 * t610 - t376 * t438 - t393 * t403 + t445 * t520 + (-t393 * t458 - t461 * t520) * t387) * t603 - (-t461 * t610 - t376 * t458 - t393 * t412 - t409 * t520 + (t393 * t438 - t477 * t520) * t387) * t604) * t394 (t382 * t631 - t397 * t444) * t609 + (t382 * t588 + t402 * t397 + (-t444 * t381 - t382 * t401 - (-t377 * t438 - t390 * t403 + t413 + (-t390 * t458 - t440) * t387) * t627 - (-t377 * t458 - t390 * t412 - t404 + (t390 * t438 - t459) * t387) * t628) * t398) * t394, 0; (t435 * t461 + t442 * t623) * t608 + ((-t432 * t521 - t658) * t435 + t442 * t587 + (t461 * t402 - (-t441 * qJD(5) + t411 * t523 + t431 * t520) * t465 - t442 * t407) * t436) * t420 (-t435 * t469 + t454 * t623) * t608 + ((-t470 * qJD(4) - t448 * t521 - t493 * t618) * t435 + t454 * t587 + (-t469 * t402 - (-t453 * qJD(5) + t423 * t523 + t447 * t520) * t465 - t454 * t407) * t436) * t420, -t435 * t568 * t608 + 0.2e1 * t451 * t585 + ((t428 * t521 + t524 * t561) * t435 + (-t407 * t436 + t587) * t451 - t568 * t650 - (-t549 * t520 + t541 * t523) * t623) * t420 (t435 * t466 + t523 * t624) * t608 + (0.2e1 * t523 * t602 - t408 * t435 + (-0.2e1 * t407 * t465 * t523 + t402 * t466 + t460 * t611) * t436) * t420, t585 * t651 + (t601 * t651 + (t401 * t465 + t407 * t443) * t436) * t420, 0;];
JaD_rot  = t1;
