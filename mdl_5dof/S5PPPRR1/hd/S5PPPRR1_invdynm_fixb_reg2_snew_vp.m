% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5PPPRR1
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5PPPRR1_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPPRR1_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR1_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_invdynm_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:58:16
% EndTime: 2019-12-05 14:58:22
% DurationCPUTime: 6.28s
% Computational Cost: add. (19854->392), mult. (32360->559), div. (0->0), fcn. (26492->10), ass. (0->265)
t580 = sin(pkin(9));
t582 = sin(pkin(7));
t585 = cos(pkin(7));
t562 = g(1) * t585 + g(2) * t582;
t581 = sin(pkin(8));
t584 = cos(pkin(8));
t631 = g(3) - qJDD(1);
t541 = -t562 * t584 - t581 * t631;
t561 = g(1) * t582 - g(2) * t585;
t555 = -qJDD(2) + t561;
t583 = cos(pkin(9));
t492 = t541 * t580 + t555 * t583;
t493 = t541 * t583 - t555 * t580;
t588 = sin(qJ(4));
t590 = cos(qJ(4));
t457 = t492 * t590 + t493 * t588;
t458 = -t492 * t588 + t493 * t590;
t614 = t457 * t588 + t458 * t590;
t408 = t457 * t590 - t458 * t588;
t644 = t583 * t408;
t382 = -t580 * t614 + t644;
t649 = t580 * t408;
t661 = t583 * t614 + t649;
t592 = qJD(4) ^ 2;
t624 = t588 * qJDD(4);
t557 = t590 * t592 + t624;
t623 = t590 * qJDD(4);
t558 = -t588 * t592 + t623;
t510 = t557 * t583 + t558 * t580;
t660 = t585 * t510;
t540 = -t562 * t581 + t584 * t631;
t535 = qJDD(3) + t540;
t496 = pkin(5) * t557 + t535 * t590;
t606 = -pkin(5) * t558 + t535 * t588;
t602 = t496 * t583 - t580 * t606;
t659 = -qJ(3) * t510 - t602;
t514 = t557 * t580 - t558 * t583;
t603 = -t580 * t496 - t583 * t606;
t430 = qJ(3) * t514 - t603;
t637 = t584 * t585;
t658 = t510 * t582 - t514 * t637;
t655 = t582 * t631;
t652 = t585 * t631;
t450 = -pkin(4) * t592 + qJDD(4) * pkin(6) + t458;
t587 = sin(qJ(5));
t589 = cos(qJ(5));
t426 = t450 * t587 - t535 * t589;
t427 = t450 * t589 + t535 * t587;
t403 = t426 * t587 + t427 * t589;
t651 = pkin(3) * t408;
t577 = t587 ^ 2;
t650 = t577 * t592;
t648 = t580 * t535;
t647 = t581 * t535;
t542 = t581 * t555;
t646 = t582 * t555;
t645 = t582 * t584;
t643 = t583 * t535;
t642 = t584 * t382;
t455 = t492 * t583 - t493 * t580;
t641 = t584 * t455;
t628 = qJD(4) * qJD(5);
t578 = t589 ^ 2;
t629 = t577 + t578;
t548 = t629 * t628;
t537 = -qJDD(5) * t590 + t548 * t588;
t539 = qJDD(5) * t588 + t548 * t590;
t478 = -t537 * t580 + t539 * t583;
t640 = t584 * t478;
t639 = t584 * t514;
t638 = t584 * t510;
t518 = t584 * t535;
t543 = t584 * t555;
t636 = t585 * t555;
t449 = -qJDD(4) * pkin(4) - pkin(6) * t592 + t457;
t442 = t587 * t449;
t569 = t587 * t592 * t589;
t563 = qJDD(5) + t569;
t635 = t587 * t563;
t564 = qJDD(5) - t569;
t634 = t587 * t564;
t443 = t589 * t449;
t633 = t589 * t563;
t632 = t589 * t564;
t630 = -pkin(4) * t449 + pkin(6) * t403;
t627 = t581 * qJDD(4);
t626 = t584 * qJDD(4);
t625 = t587 * qJDD(4);
t573 = t589 * qJDD(4);
t591 = qJD(5) ^ 2;
t566 = -t591 - t650;
t533 = -t566 * t587 - t632;
t572 = t589 * t628;
t550 = 0.2e1 * t572 + t625;
t622 = -pkin(4) * t550 + pkin(6) * t533 + t442;
t575 = t578 * t592;
t568 = -t575 - t591;
t531 = t568 * t589 - t635;
t620 = t587 * t628;
t553 = t573 - 0.2e1 * t620;
t621 = pkin(4) * t553 + pkin(6) * t531 - t443;
t619 = -pkin(1) * t581 - qJ(3);
t556 = t629 * qJDD(4);
t559 = t575 + t650;
t516 = t556 * t588 + t559 * t590;
t517 = t556 * t590 - t559 * t588;
t474 = t516 * t583 + t517 * t580;
t609 = pkin(4) * t559 + pkin(6) * t556 + t403;
t597 = pkin(3) * t516 + t609;
t384 = -pkin(2) * t474 - t597;
t475 = -t516 * t580 + t517 * t583;
t618 = qJ(2) * t475 + t384;
t607 = -pkin(3) * t557 - t458;
t416 = pkin(2) * t510 - t607;
t617 = qJ(2) * t514 + t416;
t598 = pkin(3) * t558 - t457;
t417 = pkin(2) * t514 - t598;
t616 = -qJ(2) * t510 + t417;
t534 = pkin(2) * t535;
t604 = t492 * t580 + t493 * t583;
t615 = qJ(3) * t604 - t534;
t613 = t540 * t581 + t541 * t584;
t612 = -t561 * t582 - t562 * t585;
t611 = t588 * t569;
t610 = t590 * t569;
t385 = t403 * t588 - t449 * t590;
t608 = pkin(3) * t385 + t630;
t402 = t426 * t589 - t427 * t587;
t394 = pkin(5) * t517 + t402 * t588;
t395 = -pkin(5) * t516 + t402 * t590;
t605 = -t583 * t394 - t580 * t395;
t481 = t540 * t584 - t541 * t581;
t601 = t561 * t585 - t562 * t582;
t489 = t533 * t588 - t550 * t590;
t600 = pkin(3) * t489 + t622;
t488 = t531 * t588 + t553 * t590;
t599 = pkin(3) * t488 + t621;
t405 = -pkin(3) * t535 + pkin(5) * t614;
t596 = pkin(5) * t649 + qJ(3) * t661 + t405 * t583 - t534;
t386 = t403 * t590 + t449 * t588;
t361 = pkin(5) * t386 - (-pkin(4) * t590 - pkin(6) * t588 - pkin(3)) * t402;
t367 = -pkin(5) * t385 - (pkin(4) * t588 - pkin(6) * t590) * t402;
t372 = -t385 * t580 + t386 * t583;
t595 = pkin(2) * t402 + qJ(3) * t372 + t361 * t583 + t367 * t580;
t527 = t568 * t587 + t633;
t410 = -pkin(4) * t527 + t426;
t424 = -pkin(6) * t527 + t442;
t490 = t531 * t590 - t553 * t588;
t387 = -pkin(3) * t527 + pkin(5) * t490 + t410 * t590 + t424 * t588;
t389 = -pkin(5) * t488 - t410 * t588 + t424 * t590;
t447 = -t488 * t580 + t490 * t583;
t594 = -pkin(2) * t527 + qJ(3) * t447 + t387 * t583 + t389 * t580;
t529 = t566 * t589 - t634;
t411 = -pkin(4) * t529 + t427;
t425 = -pkin(6) * t529 + t443;
t491 = t533 * t590 + t550 * t588;
t388 = -pkin(3) * t529 + pkin(5) * t491 + t411 * t590 + t425 * t588;
t390 = -pkin(5) * t489 - t411 * t588 + t425 * t590;
t448 = -t489 * t580 + t491 * t583;
t593 = -pkin(2) * t529 + qJ(3) * t448 + t388 * t583 + t390 * t580;
t567 = t575 - t591;
t565 = t591 - t650;
t560 = -t575 + t650;
t552 = t573 - t620;
t551 = t572 + t625;
t538 = t551 * t589 - t577 * t628;
t536 = -t552 * t587 - t578 * t628;
t532 = -t565 * t587 + t633;
t530 = t567 * t589 - t634;
t528 = t565 * t589 + t635;
t526 = t567 * t587 + t632;
t525 = (t551 + t572) * t587;
t524 = -t552 * t589 + t572 * t587;
t509 = -t550 * t587 + t553 * t589;
t508 = t550 * t589 + t553 * t587;
t507 = t581 * t510;
t506 = t581 * t514;
t505 = t538 * t590 - t611;
t504 = t536 * t590 + t611;
t503 = t538 * t588 + t610;
t502 = t536 * t588 - t610;
t501 = t532 * t590 + t587 * t624;
t500 = t530 * t590 + t573 * t588;
t499 = t532 * t588 - t587 * t623;
t498 = t530 * t588 - t589 * t623;
t484 = t509 * t590 + t560 * t588;
t483 = t509 * t588 - t560 * t590;
t477 = t537 * t583 + t539 * t580;
t476 = t581 * t478;
t473 = pkin(1) * t555 + qJ(2) * t613;
t472 = -t510 * t637 - t514 * t582;
t471 = t514 * t585 - t582 * t638;
t470 = -t493 * t581 + t518 * t583;
t469 = -t492 * t581 + t518 * t580;
t468 = t493 * t584 + t581 * t643;
t467 = t492 * t584 + t580 * t647;
t466 = -t503 * t580 + t505 * t583;
t465 = -t502 * t580 + t504 * t583;
t464 = t503 * t583 + t505 * t580;
t463 = t502 * t583 + t504 * t580;
t462 = -t499 * t580 + t501 * t583;
t461 = -t498 * t580 + t500 * t583;
t460 = t499 * t583 + t501 * t580;
t459 = t498 * t583 + t500 * t580;
t451 = t581 * t455;
t446 = t489 * t583 + t491 * t580;
t445 = t488 * t583 + t490 * t580;
t441 = -t483 * t580 + t484 * t583;
t440 = t483 * t583 + t484 * t580;
t439 = t466 * t584 + t525 * t581;
t438 = t465 * t584 - t524 * t581;
t437 = t466 * t581 - t525 * t584;
t436 = t465 * t581 + t524 * t584;
t435 = t462 * t584 + t528 * t581;
t434 = t461 * t584 + t526 * t581;
t433 = t462 * t581 - t528 * t584;
t432 = t461 * t581 - t526 * t584;
t429 = t584 * t604 + t647;
t428 = t581 * t604 - t518;
t421 = t448 * t584 + t529 * t581;
t420 = t447 * t584 + t527 * t581;
t419 = t448 * t581 - t529 * t584;
t418 = t447 * t581 - t527 * t584;
t415 = t441 * t584 + t508 * t581;
t414 = t441 * t581 - t508 * t584;
t413 = t514 * t619 + t603;
t412 = -t510 * t619 + t602;
t404 = -pkin(1) * t428 - t615;
t399 = -pkin(2) * t446 - t600;
t398 = -pkin(2) * t445 - t599;
t397 = -t581 * t617 - t584 * t659;
t396 = t584 * t430 - t581 * t616;
t393 = -qJ(2) * t428 - (pkin(2) * t581 - qJ(3) * t584) * t455;
t392 = pkin(1) * t510 - t581 * t659 + t584 * t617;
t391 = pkin(1) * t514 + t581 * t430 + t584 * t616;
t379 = t581 * t382;
t378 = qJ(2) * t429 - (-pkin(2) * t584 - qJ(3) * t581 - pkin(1)) * t455;
t377 = t584 * t661 + t647;
t376 = t581 * t661 - t518;
t375 = pkin(2) * t382 + t651;
t374 = -qJ(3) * t474 - t394 * t580 + t395 * t583;
t373 = t475 * t619 + t605;
t371 = t385 * t583 + t386 * t580;
t370 = pkin(5) * t644 + qJ(3) * t382 - t405 * t580;
t369 = -qJ(3) * t446 - t388 * t580 + t390 * t583;
t368 = -qJ(3) * t445 - t387 * t580 + t389 * t583;
t366 = -pkin(1) * t419 - t593;
t365 = -pkin(1) * t418 - t594;
t364 = t372 * t584 - t402 * t581;
t363 = t372 * t581 + t402 * t584;
t362 = t584 * t374 - t581 * t618;
t360 = -pkin(1) * t474 + t581 * t374 + t584 * t618;
t359 = -qJ(2) * t419 + t369 * t584 - t399 * t581;
t358 = -qJ(2) * t418 + t368 * t584 - t398 * t581;
t357 = -pkin(1) * t376 - t596;
t356 = -pkin(1) * t446 + qJ(2) * t421 + t369 * t581 + t399 * t584;
t355 = -pkin(1) * t445 + qJ(2) * t420 + t368 * t581 + t398 * t584;
t354 = -pkin(2) * t371 - t608;
t353 = -qJ(2) * t376 + t370 * t584 - t375 * t581;
t352 = pkin(1) * t382 + qJ(2) * t377 + t370 * t581 + t375 * t584;
t351 = -qJ(3) * t371 - t361 * t580 + t367 * t583;
t350 = -pkin(1) * t363 - t595;
t349 = -qJ(2) * t363 + t351 * t584 - t354 * t581;
t348 = -pkin(1) * t371 + qJ(2) * t364 + t351 * t581 + t354 * t584;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t655, -t652, -t601, -qJ(1) * t601, 0, 0, 0, 0, 0, 0, -t540 * t582 - t581 * t636, -t541 * t582 - t584 * t636, t585 * t481, -qJ(1) * (t582 * t613 + t636) - (pkin(1) * t582 - qJ(2) * t585) * t481, 0, 0, 0, 0, 0, 0, t469 * t585 - t582 * t643, t470 * t585 + t582 * t648, t455 * t637 + t582 * t604, t585 * t393 - t582 * t404 - qJ(1) * (t429 * t582 + t455 * t585), 0, 0, t658, 0, t472, t585 * t627, -qJ(1) * t471 + t396 * t585 - t412 * t582, t585 * t397 - t582 * t413 - qJ(1) * (t514 * t645 + t660), t382 * t637 + t582 * t661, t585 * t353 - t582 * t357 - qJ(1) * (t377 * t582 + t382 * t585), t439 * t585 + t464 * t582, t415 * t585 + t440 * t582, t435 * t585 + t460 * t582, t438 * t585 + t463 * t582, t434 * t585 + t459 * t582, t477 * t582 + t478 * t637, t585 * t358 - t582 * t365 - qJ(1) * (t420 * t582 - t445 * t585), t585 * t359 - t582 * t366 - qJ(1) * (t421 * t582 - t446 * t585), t585 * t362 - t582 * t373 - qJ(1) * (-t474 * t585 + t475 * t645), t585 * t349 - t582 * t350 - qJ(1) * (t364 * t582 - t371 * t585); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t652, -t655, t612, qJ(1) * t612, 0, 0, 0, 0, 0, 0, t540 * t585 - t581 * t646, t541 * t585 - t543 * t582, t582 * t481, qJ(1) * (t585 * t613 - t646) - (-pkin(1) * t585 - qJ(2) * t582) * t481, 0, 0, 0, 0, 0, 0, t469 * t582 + t585 * t643, t470 * t582 - t585 * t648, t582 * t641 - t585 * t604, t582 * t393 + t585 * t404 + qJ(1) * (t429 * t585 - t455 * t582), 0, 0, -t582 * t639 - t660, 0, t471, t582 * t627, qJ(1) * t472 + t396 * t582 + t412 * t585, -qJ(1) * t658 + t582 * t397 + t585 * t413, t582 * t642 - t585 * t661, t582 * t353 + t585 * t357 + qJ(1) * (t377 * t585 - t382 * t582), t439 * t582 - t464 * t585, t415 * t582 - t440 * t585, t435 * t582 - t460 * t585, t438 * t582 - t463 * t585, t434 * t582 - t459 * t585, -t477 * t585 + t582 * t640, t582 * t358 + t585 * t365 + qJ(1) * (t420 * t585 + t445 * t582), t582 * t359 + t585 * t366 + qJ(1) * (t421 * t585 + t446 * t582), t582 * t362 + t585 * t373 + qJ(1) * (t474 * t582 + t475 * t637), t582 * t349 + t585 * t350 + qJ(1) * (t364 * t585 + t371 * t582); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t561, t562, 0, 0, 0, 0, 0, 0, 0, 0, t543, -t542, t613, t473, 0, 0, 0, 0, 0, 0, t467, t468, t451, t378, 0, 0, -t506, 0, -t507, -t626, t391, t392, t379, t352, t437, t414, t433, t436, t432, t476, t355, t356, t360, t348; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t631, -t561, 0, 0, 0, 0, 0, 0, 0, -t542, -t543, t481, qJ(2) * t481, 0, 0, 0, 0, 0, 0, t469, t470, t641, t393, 0, 0, -t639, 0, -t638, t627, t396, t397, t642, t353, t439, t415, t435, t438, t434, t640, t358, t359, t362, t349; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t631, 0, -t562, 0, 0, 0, 0, 0, 0, 0, t540, t541, 0, pkin(1) * t481, 0, 0, 0, 0, 0, 0, t643, -t648, -t604, t404, 0, 0, -t510, 0, t514, 0, t412, t413, -t661, t357, -t464, -t440, -t460, -t463, -t459, -t477, t365, t366, t373, t350; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t561, t562, 0, 0, 0, 0, 0, 0, 0, 0, t543, -t542, t613, t473, 0, 0, 0, 0, 0, 0, t467, t468, t451, t378, 0, 0, -t506, 0, -t507, -t626, t391, t392, t379, t352, t437, t414, t433, t436, t432, t476, t355, t356, t360, t348; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t555, t540, 0, 0, 0, 0, 0, 0, 0, t648, t643, t455, qJ(3) * t455, 0, 0, -t514, 0, -t510, 0, t430, -t659, t382, t370, t466, t441, t462, t465, t461, t478, t368, t369, t374, t351; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t555, 0, t541, 0, 0, 0, 0, 0, 0, 0, t492, t493, 0, pkin(2) * t455, 0, 0, 0, 0, 0, -qJDD(4), t417, t416, 0, t375, -t525, -t508, -t528, t524, -t526, 0, t398, t399, t384, t354; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t540, -t541, 0, 0, 0, 0, 0, 0, 0, 0, -t643, t648, t604, t615, 0, 0, t510, 0, -t514, 0, t659, t430, t661, t596, t464, t440, t460, t463, t459, t477, t594, t593, qJ(3) * t475 - t605, t595; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t535, t492, 0, 0, 0, t558, 0, -t557, 0, t606, t496, t408, pkin(5) * t408, t505, t484, t501, t504, t500, t539, t389, t390, t395, t367; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t535, 0, t493, 0, 0, 0, t557, 0, t558, 0, -t496, t606, t614, t405, t503, t483, t499, t502, t498, t537, t387, t388, t394, t361; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t492, -t493, 0, 0, 0, 0, 0, 0, 0, qJDD(4), t598, t607, 0, -t651, t525, t508, t528, -t524, t526, 0, t599, t600, t597, t608; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4), 0, -t592, 0, 0, t535, t457, 0, t538, t509, t532, t536, t530, t548, t424, t425, t402, pkin(6) * t402; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t592, 0, qJDD(4), 0, -t535, 0, t458, 0, t569, -t560, -t625, -t569, -t573, -qJDD(5), t410, t411, 0, pkin(4) * t402; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4), -t457, -t458, 0, 0, t525, t508, t528, -t524, t526, 0, t621, t622, t609, t630; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t551, t553, t563, -t572, t567, t572, 0, t449, t426, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t620, t550, t565, t552, t564, -t620, -t449, 0, t427, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t569, t560, t625, t569, t573, qJDD(5), -t426, -t427, 0, 0;];
m_new_reg = t1;