% Calculate vector of inverse dynamics joint torques for
% S5RRPRR14
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR14_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR14_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR14_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR14_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR14_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR14_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:38:56
% EndTime: 2019-12-31 20:39:11
% DurationCPUTime: 9.34s
% Computational Cost: add. (6195->562), mult. (15648->785), div. (0->0), fcn. (12844->14), ass. (0->248)
t548 = cos(pkin(5));
t646 = qJD(1) * t548;
t532 = qJD(2) + t646;
t545 = sin(pkin(10));
t547 = cos(pkin(10));
t551 = sin(qJ(2));
t546 = sin(pkin(5));
t647 = qJD(1) * t546;
t619 = t551 * t647;
t475 = t532 * t547 - t545 * t619;
t476 = t532 * t545 + t547 * t619;
t550 = sin(qJ(4));
t554 = cos(qJ(4));
t423 = -t554 * t475 + t476 * t550;
t420 = qJD(5) + t423;
t695 = t420 ^ 2;
t555 = cos(qJ(2));
t645 = qJD(1) * t555;
t618 = t546 * t645;
t517 = -qJD(4) + t618;
t694 = t423 * t517;
t501 = t545 * t550 - t554 * t547;
t663 = t546 * t555;
t568 = t501 * t663;
t693 = -qJD(1) * t568 + t501 * qJD(4);
t502 = t545 * t554 + t547 * t550;
t569 = t502 * t663;
t650 = -qJD(1) * t569 + t502 * qJD(4);
t552 = sin(qJ(1));
t657 = t552 * t555;
t556 = cos(qJ(1));
t658 = t551 * t556;
t497 = t548 * t658 + t657;
t542 = pkin(10) + qJ(4);
t539 = sin(t542);
t540 = cos(t542);
t662 = t546 * t556;
t447 = t497 * t540 - t539 * t662;
t655 = t555 * t556;
t659 = t551 * t552;
t496 = -t548 * t655 + t659;
t549 = sin(qJ(5));
t553 = cos(qJ(5));
t692 = t447 * t549 - t496 * t553;
t691 = t447 * t553 + t496 * t549;
t684 = pkin(1) * t551;
t649 = pkin(7) * t663 + t548 * t684;
t483 = qJ(3) * t548 + t649;
t579 = -pkin(2) * t555 - qJ(3) * t551 - pkin(1);
t484 = t579 * t546;
t428 = -t483 * t545 + t547 * t484;
t665 = t546 * t551;
t493 = t545 * t548 + t547 * t665;
t393 = -pkin(3) * t663 - pkin(8) * t493 + t428;
t429 = t547 * t483 + t545 * t484;
t492 = t545 * t665 - t548 * t547;
t403 = -pkin(8) * t492 + t429;
t690 = t550 * t393 + t554 * t403;
t591 = pkin(2) * t551 - qJ(3) * t555;
t488 = t591 * t647;
t627 = pkin(1) * t646;
t489 = -pkin(7) * t619 + t555 * t627;
t432 = t547 * t488 - t489 * t545;
t572 = (-pkin(8) * t547 * t555 + pkin(3) * t551) * t546;
t404 = qJD(1) * t572 + t432;
t433 = t545 * t488 + t547 * t489;
t601 = t545 * t618;
t416 = -pkin(8) * t601 + t433;
t680 = pkin(8) + qJ(3);
t512 = t680 * t545;
t513 = t680 * t547;
t581 = -t512 * t554 - t513 * t550;
t689 = -qJD(3) * t501 + qJD(4) * t581 - t550 * t404 - t554 * t416;
t460 = -t512 * t550 + t513 * t554;
t688 = qJD(3) * t502 + qJD(4) * t460 + t404 * t554 - t416 * t550;
t632 = qJDD(1) * t548;
t531 = qJDD(2) + t632;
t687 = -pkin(2) * t531 + qJDD(3);
t630 = qJDD(1) * t555;
t530 = t546 * t630;
t633 = qJD(1) * qJD(2);
t615 = t551 * t633;
t599 = t546 * t615;
t487 = qJDD(4) - t530 + t599;
t525 = pkin(7) * t618;
t490 = t551 * t627 + t525;
t463 = qJ(3) * t532 + t490;
t468 = qJD(1) * t484;
t406 = -t463 * t545 + t547 * t468;
t385 = -pkin(3) * t618 - pkin(8) * t476 + t406;
t407 = t547 * t463 + t545 * t468;
t388 = pkin(8) * t475 + t407;
t361 = t385 * t550 + t388 * t554;
t626 = pkin(1) * qJD(2) * t548;
t603 = qJD(1) * t626;
t625 = pkin(1) * t632;
t620 = -pkin(7) * t530 - t551 * t625 - t555 * t603;
t562 = -pkin(7) * t599 - t620;
t411 = qJ(3) * t531 + qJD(3) * t532 + t562;
t559 = qJD(2) * t591 - qJD(3) * t551;
t417 = (qJD(1) * t559 + qJDD(1) * t579) * t546;
t381 = -t411 * t545 + t547 * t417;
t614 = t555 * t633;
t631 = qJDD(1) * t551;
t577 = t614 + t631;
t563 = t577 * t546;
t670 = t531 * t545;
t445 = t547 * t563 + t670;
t576 = t615 - t630;
t364 = pkin(3) * t546 * t576 - pkin(8) * t445 + t381;
t382 = t547 * t411 + t545 * t417;
t507 = t547 * t531;
t444 = t545 * t563 - t507;
t368 = -pkin(8) * t444 + t382;
t558 = -qJD(4) * t361 + t554 * t364 - t368 * t550;
t349 = -pkin(4) * t487 - t558;
t499 = -t548 * t659 + t655;
t664 = t546 * t552;
t449 = -t499 * t539 + t540 * t664;
t672 = t497 * t539;
t573 = g(1) * t449 + g(2) * (-t540 * t662 - t672) + g(3) * (-t539 * t665 + t540 * t548);
t582 = t475 * t550 + t476 * t554;
t686 = t420 * (pkin(4) * t582 + t420 * pkin(9)) + t349 + t573;
t379 = qJD(4) * t582 + t554 * t444 + t445 * t550;
t466 = t559 * t546;
t642 = qJD(2) * t551;
t617 = t546 * t642;
t580 = -pkin(7) * t617 + t555 * t626;
t472 = qJD(3) * t548 + t580;
t412 = t547 * t466 - t472 * t545;
t391 = qJD(2) * t572 + t412;
t413 = t545 * t466 + t547 * t472;
t641 = qJD(2) * t555;
t616 = t546 * t641;
t600 = t545 * t616;
t401 = -pkin(8) * t600 + t413;
t685 = -qJD(4) * t690 + t391 * t554 - t401 * t550;
t682 = g(1) * t556;
t681 = g(3) * t546;
t679 = MDP(6) * t546;
t639 = qJD(4) * t554;
t640 = qJD(4) * t550;
t378 = -t550 * t444 + t554 * t445 + t475 * t639 - t476 * t640;
t637 = qJD(5) * t553;
t621 = t553 * t378 + t549 * t487 - t517 * t637;
t638 = qJD(5) * t549;
t356 = -t582 * t638 + t621;
t678 = t356 * t549;
t675 = t582 * t549;
t398 = t553 * t517 + t675;
t677 = t398 * t420;
t400 = -t517 * t549 + t553 * t582;
t676 = t400 * t420;
t671 = t502 * t553;
t669 = t531 * t548;
t668 = t540 * t549;
t667 = t540 * t553;
t541 = t546 ^ 2;
t666 = t541 * qJD(1) ^ 2;
t377 = qJDD(5) + t379;
t661 = t549 * t377;
t660 = t549 * t555;
t375 = t553 * t377;
t656 = t553 * t555;
t652 = pkin(4) * t619 + t688;
t491 = pkin(7) * t616 + t551 * t626;
t543 = t551 ^ 2;
t648 = -t555 ^ 2 + t543;
t644 = qJD(2) * t545;
t643 = qJD(2) * t547;
t636 = qJD(2) - t532;
t634 = qJ(3) * qJDD(1);
t629 = 0.2e1 * t541;
t628 = g(3) * t663;
t533 = pkin(7) * t665;
t624 = t555 * t666;
t623 = t546 * t660;
t622 = t546 * t656;
t454 = pkin(3) * t601 + t490;
t455 = pkin(3) * t600 + t491;
t538 = -pkin(3) * t547 - pkin(2);
t575 = t550 * t364 + t554 * t368 + t385 * t639 - t388 * t640;
t348 = pkin(9) * t487 + t575;
t602 = qJD(2) * t525 + qJDD(1) * t533 + t551 * t603 - t555 * t625;
t427 = t602 + t687;
t387 = pkin(3) * t444 + t427;
t353 = pkin(4) * t379 - pkin(9) * t378 + t387;
t612 = -t549 * t348 + t553 * t353;
t610 = t378 * t549 - t553 * t487;
t609 = t693 * t549 - t553 * t619;
t608 = t549 * t619 + t693 * t553;
t606 = t420 * t553;
t605 = t532 + t646;
t604 = t531 + t632;
t498 = t548 * t657 + t658;
t598 = g(1) * t498 + g(2) * t496;
t597 = -g(1) * t496 + g(2) * t498;
t596 = g(1) * t499 + g(2) * t497;
t595 = g(1) * t497 - g(2) * t499;
t594 = g(2) * t552 + t682;
t443 = pkin(4) * t501 - pkin(9) * t502 + t538;
t593 = pkin(9) * t619 - qJD(5) * t443 - t689;
t592 = -t650 * pkin(4) - t693 * pkin(9) + qJD(5) * t460 + t454;
t590 = t553 * t348 + t549 * t353;
t359 = -pkin(9) * t517 + t361;
t456 = -pkin(2) * t532 + qJD(3) - t489;
t421 = -pkin(3) * t475 + t456;
t367 = pkin(4) * t423 - pkin(9) * t582 + t421;
t351 = t359 * t553 + t367 * t549;
t589 = t359 * t549 - t367 * t553;
t370 = -pkin(9) * t663 + t690;
t436 = t554 * t492 + t493 * t550;
t437 = -t492 * t550 + t493 * t554;
t486 = t533 + (-pkin(1) * t555 - pkin(2)) * t548;
t441 = pkin(3) * t492 + t486;
t380 = pkin(4) * t436 - pkin(9) * t437 + t441;
t588 = t370 * t553 + t380 * t549;
t587 = -t370 * t549 + t380 * t553;
t360 = t385 * t554 - t388 * t550;
t585 = t393 * t554 - t403 * t550;
t578 = -t427 + t598;
t418 = t437 * t549 + t622;
t574 = t550 * t391 + t393 * t639 + t554 * t401 - t403 * t640;
t571 = t502 * t637 - t609;
t570 = -t502 * t638 - t608;
t567 = -t598 + t628;
t566 = -g(3) * t665 - t596;
t564 = -qJ(3) * t642 + (qJD(3) - t456) * t555;
t561 = t598 - t602;
t358 = pkin(4) * t517 - t360;
t560 = -pkin(9) * t377 + (t358 + t360) * t420;
t482 = t539 * t548 + t540 * t665;
t450 = t499 * t540 + t539 * t664;
t419 = t437 * t553 - t623;
t415 = t450 * t553 + t498 * t549;
t414 = -t450 * t549 + t498 * t553;
t396 = qJD(2) * t569 + qJD(4) * t437;
t395 = -qJD(2) * t568 - qJD(4) * t436;
t372 = -qJD(5) * t623 + t395 * t549 + t437 * t637 - t553 * t617;
t371 = -qJD(5) * t418 + t395 * t553 + t549 * t617;
t369 = pkin(4) * t663 - t585;
t365 = pkin(4) * t396 - pkin(9) * t395 + t455;
t357 = qJD(5) * t400 + t610;
t355 = -pkin(4) * t617 - t685;
t354 = pkin(9) * t617 + t574;
t347 = -t351 * qJD(5) + t612;
t346 = -t589 * qJD(5) + t590;
t1 = [(-t356 * t418 - t357 * t419 - t371 * t398 - t372 * t400) * MDP(23) + (t356 * t419 + t371 * t400) * MDP(22) + (-t491 * t532 - t533 * t531 - t602 * t548 + (t555 * t669 - t576 * t629) * pkin(1) + t595) * MDP(9) + (t382 * t429 + t407 * t413 + t381 * t428 + t406 * t412 + t427 * t486 + t456 * t491 - g(1) * (-pkin(1) * t552 - pkin(2) * t497 + pkin(7) * t662 - qJ(3) * t496) - g(2) * (pkin(1) * t556 + pkin(2) * t499 + pkin(7) * t664 + qJ(3) * t498)) * MDP(14) + (-pkin(1) * t577 * t629 - t531 * t649 - t532 * t580 - t548 * t562 + t597) * MDP(10) + (t396 * t517 - t436 * t487) * MDP(18) + (-t395 * t517 + t437 * t487) * MDP(17) + (0.2e1 * (t551 * t630 - t633 * t648) * MDP(5) + (qJDD(1) * t543 + 0.2e1 * t551 * t614) * MDP(4)) * t541 + t594 * MDP(3) + (t427 * t493 + t486 * t445 + t491 * t476 - t595 * t545) * MDP(12) + (t427 * t492 + t486 * t444 - t491 * t475 + t595 * t547) * MDP(11) + (g(1) * t447 - g(2) * t450 + t441 * t379 + t387 * t436 + t421 * t396 + t455 * t423 + t585 * t487 - t517 * t685) * MDP(20) + (-g(1) * t672 - g(2) * t449 + t441 * t378 + t387 * t437 + t421 * t395 + t455 * t582 - t487 * t690 + t574 * t517) * MDP(21) + ((-t361 * t642 - t540 * t682 + t555 * t575) * MDP(21) + (t360 * t642 - t555 * t558) * MDP(20) + (-t378 * t555 + t582 * t642) * MDP(17) + (-t594 * t547 + (-qJD(1) * t429 - t407) * t642 + (qJD(1) * t413 + qJDD(1) * t429 + t456 * t643 + t382) * t555) * MDP(12) + (-t594 * t545 + (qJD(1) * t428 + t406) * t642 + (-qJD(1) * t412 - qJDD(1) * t428 + t456 * t644 - t381) * t555) * MDP(11) + (t379 * t555 - t423 * t642) * MDP(18) + (t555 * t604 - t605 * t642) * MDP(7) + (-t487 * t555 - t517 * t642) * MDP(19)) * t546 + (t378 * t437 + t395 * t582) * MDP(15) + (-t378 * t436 - t379 * t437 - t395 * t423 - t396 * t582) * MDP(16) + (-t381 * t493 - t382 * t492 - t412 * t476 + t413 * t475 - t428 * t445 - t429 * t444 + (-t406 * t547 - t407 * t545) * t616 - t597) * MDP(13) + ((-qJD(5) * t588 - t354 * t549 + t365 * t553) * t420 + t587 * t377 + t347 * t436 - t589 * t396 + t355 * t398 + t369 * t357 + t349 * t418 + t358 * t372 + g(1) * t691 - g(2) * t415) * MDP(27) + (-(qJD(5) * t587 + t354 * t553 + t365 * t549) * t420 - t588 * t377 - t346 * t436 - t351 * t396 + t355 * t400 + t369 * t356 + t349 * t419 + t358 * t371 - g(1) * t692 - g(2) * t414) * MDP(28) + (-t357 * t436 - t372 * t420 - t377 * t418 - t396 * t398) * MDP(25) + (t356 * t436 + t371 * t420 + t377 * t419 + t396 * t400) * MDP(24) + (t377 * t436 + t396 * t420) * MDP(26) + qJDD(1) * MDP(1) + (g(1) * t552 - g(2) * t556) * MDP(2) + MDP(8) * t669 + (t551 * t604 + t605 * t641) * t679; (t490 * t532 + t666 * t684 + t561 - t628) * MDP(9) + (t609 * t400 + t608 * t398 + (-t678 - t357 * t553 + (t398 * t549 - t400 * t553) * qJD(5)) * t502) * MDP(23) + (t356 * t671 + t400 * t570) * MDP(22) + (pkin(1) * t624 + t489 * t532 + (pkin(7) * t633 + g(3)) * t665 + t596 + t620) * MDP(10) + (-t357 * t501 - t398 * t650 - t420 * t571 - t502 * t661) * MDP(25) + t517 * MDP(19) * t619 + (t356 * t501 + t375 * t502 + t400 * t650 + t420 * t570) * MDP(24) + (t377 * t501 + t420 * t650) * MDP(26) + (t423 * t619 - t487 * t501 + t517 * t650) * MDP(18) - t551 * MDP(4) * t624 + (-t360 * t619 + t538 * t379 + t387 * t501 + t650 * t421 - t454 * t423 + t487 * t581 + t517 * t688 - t567 * t540) * MDP(20) + ((t443 * t553 - t460 * t549) * t377 + t347 * t501 - t581 * t357 + t349 * t549 * t502 - g(1) * (-t498 * t667 + t499 * t549) - g(2) * (-t496 * t667 + t497 * t549) - (t540 * t656 + t549 * t551) * t681 + (t549 * t593 - t553 * t592) * t420 + t652 * t398 - t650 * t589 + t571 * t358) * MDP(27) + (-(t443 * t549 + t460 * t553) * t377 - t346 * t501 - t581 * t356 + t349 * t671 - g(1) * (t498 * t668 + t499 * t553) - g(2) * (t496 * t668 + t497 * t553) - (-t540 * t660 + t551 * t553) * t681 + (t549 * t592 + t553 * t593) * t420 + t652 * t400 - t650 * t351 + t570 * t358) * MDP(28) + (-t619 * t636 + t530) * MDP(7) + (t361 * t619 + t538 * t378 + t387 * t502 - t421 * t693 - t454 * t582 - t460 * t487 + t517 * t689 + t567 * t539) * MDP(21) + (t378 * t502 - t582 * t693) * MDP(15) + (-t378 * t501 - t379 * t502 + t423 * t693 - t582 * t650) * MDP(16) + (t487 * t502 + t517 * t693 - t582 * t619) * MDP(17) + (t432 * t476 - t433 * t475 + (-qJ(3) * t444 + qJD(3) * t475 + t406 * t618 + t382) * t547 + (qJ(3) * t445 + qJD(3) * t476 + t407 * t618 - t381) * t545 + t566) * MDP(13) + (-t406 * t432 - t407 * t433 - t456 * t490 + (-t406 * t545 + t407 * t547) * qJD(3) + (-t427 - t567) * pkin(2) + (-t381 * t545 + t382 * t547 + t566) * qJ(3)) * MDP(14) + t648 * MDP(5) * t666 + t531 * MDP(8) + (-pkin(2) * t445 - t476 * t490 - t578 * t545 + ((g(3) * t545 + t547 * t634) * t555 + (t407 * t551 - t433 * t555 + t547 * t564) * qJD(1)) * t546) * MDP(12) + (-pkin(2) * t444 + t475 * t490 + t578 * t547 + ((-g(3) * t547 + t545 * t634) * t555 + (-t406 * t551 + t432 * t555 + t545 * t564) * qJD(1)) * t546) * MDP(11) + (t636 * t645 + t631) * t679; -t507 * MDP(11) + MDP(12) * t670 + (-t475 ^ 2 - t476 ^ 2) * MDP(13) + (t406 * t476 - t407 * t475 - t561 + t687) * MDP(14) + (-t517 * t582 + t379) * MDP(20) + (t378 + t694) * MDP(21) + (-t398 * t582 + t375) * MDP(27) + (-t400 * t582 - t661) * MDP(28) + ((t545 * MDP(11) + t547 * MDP(12)) * t631 + (g(3) * MDP(14) + ((-t476 + t644) * MDP(11) + (-t475 + t643) * MDP(12)) * qJD(1)) * t555) * t546 - (MDP(27) * t549 + MDP(28) * t553) * t695; -t423 ^ 2 * MDP(16) + (t378 - t694) * MDP(17) - t379 * MDP(18) + t487 * MDP(19) + (-t361 * t517 + t558 - t573) * MDP(20) + (g(1) * t450 + g(2) * t447 + g(3) * t482 - t360 * t517 + t421 * t423 - t575) * MDP(21) + (t400 * t606 + t678) * MDP(22) + ((t356 - t677) * t553 + (-t357 - t676) * t549) * MDP(23) + (t420 * t606 + t661) * MDP(24) + (-t695 * t549 + t375) * MDP(25) + (-pkin(4) * t357 - t361 * t398 + t560 * t549 - t553 * t686) * MDP(27) + (-pkin(4) * t356 - t361 * t400 + t549 * t686 + t560 * t553) * MDP(28) + (MDP(15) * t423 + MDP(16) * t582 - t517 * MDP(18) - t421 * MDP(20) - t400 * MDP(24) + t398 * MDP(25) - t420 * MDP(26) + MDP(27) * t589 + t351 * MDP(28)) * t582; t400 * t398 * MDP(22) + (-t398 ^ 2 + t400 ^ 2) * MDP(23) + (t621 + t677) * MDP(24) + (-t610 + t676) * MDP(25) + t377 * MDP(26) + (t351 * t420 - t358 * t400 - g(1) * t414 + g(2) * t692 - g(3) * (-t482 * t549 - t622) + t612) * MDP(27) + (-t589 * t420 + t358 * t398 + g(1) * t415 + g(2) * t691 - g(3) * (-t482 * t553 + t623) - t590) * MDP(28) + (-MDP(24) * t675 - MDP(25) * t400 - MDP(27) * t351 + MDP(28) * t589) * qJD(5);];
tau = t1;
