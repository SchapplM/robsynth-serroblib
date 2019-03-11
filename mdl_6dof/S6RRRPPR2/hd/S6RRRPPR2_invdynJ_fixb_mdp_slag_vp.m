% Calculate vector of inverse dynamics joint torques for
% S6RRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPPR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPPR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:26:55
% EndTime: 2019-03-09 15:27:04
% DurationCPUTime: 6.87s
% Computational Cost: add. (7374->532), mult. (17529->646), div. (0->0), fcn. (13038->14), ass. (0->248)
t545 = qJ(2) + qJ(3);
t535 = sin(t545);
t536 = cos(t545);
t551 = sin(qJ(1));
t555 = cos(qJ(1));
t597 = g(1) * t555 + g(2) * t551;
t693 = -g(3) * t536 + t535 * t597;
t553 = cos(qJ(3));
t554 = cos(qJ(2));
t636 = qJD(1) * t554;
t616 = t553 * t636;
t549 = sin(qJ(3));
t550 = sin(qJ(2));
t637 = qJD(1) * t550;
t617 = t549 * t637;
t472 = -t616 + t617;
t474 = -t549 * t636 - t553 * t637;
t546 = sin(pkin(10));
t547 = cos(pkin(10));
t441 = t472 * t547 - t474 * t546;
t670 = pkin(5) * t441;
t541 = qJD(2) + qJD(3);
t548 = sin(qJ(6));
t552 = cos(qJ(6));
t428 = -t552 * t441 + t541 * t548;
t586 = -t472 * t546 - t547 * t474;
t684 = qJD(6) + t586;
t692 = t428 * t684;
t430 = t441 * t548 + t541 * t552;
t691 = t430 * t684;
t690 = t441 * t541;
t678 = pkin(7) + pkin(8);
t499 = t678 * t554;
t492 = qJD(1) * t499;
t479 = t553 * t492;
t498 = t678 * t550;
t490 = qJD(1) * t498;
t606 = t490 * t549 - t479;
t663 = qJ(4) * t472;
t431 = t606 + t663;
t468 = t474 * qJ(4);
t475 = t549 * t492;
t641 = -t553 * t490 - t475;
t432 = t468 + t641;
t516 = t546 * t549 * pkin(2);
t634 = qJD(3) * t553;
t642 = -qJD(3) * t516 - t431 * t546 + (pkin(2) * t634 - t432) * t547;
t624 = qJD(1) * qJD(2);
t614 = t554 * t624;
t623 = qJDD(1) * t550;
t688 = t614 + t623;
t687 = t586 ^ 2;
t671 = pkin(5) * t586;
t537 = t554 * pkin(2);
t665 = pkin(1) + t537;
t603 = t684 * t548;
t645 = -qJD(5) - t642;
t651 = t547 * t549;
t643 = t547 * t431 - t432 * t546 + (t546 * t553 + t651) * qJD(3) * pkin(2);
t664 = qJD(2) * pkin(2);
t482 = -t490 + t664;
t607 = t553 * t482 - t475;
t423 = t468 + t607;
t640 = -t549 * t498 + t553 * t499;
t534 = pkin(10) + t545;
t521 = sin(t534);
t522 = cos(t534);
t594 = t522 * pkin(4) + t521 * qJ(5);
t686 = g(1) * t551 - g(2) * t555;
t622 = qJDD(1) * t554;
t425 = qJD(3) * t616 - t541 * t617 + t549 * t622 + t688 * t553;
t485 = t549 * t554 + t550 * t553;
t450 = t541 * t485;
t593 = t549 * t623 - t553 * t622;
t426 = qJD(1) * t450 + t593;
t400 = t425 * t547 - t426 * t546;
t685 = -qJ(5) * t400 - qJD(5) * t586;
t585 = -t482 * t549 - t479;
t424 = -t585 - t663;
t418 = t546 * t424;
t398 = t423 * t547 - t418;
t627 = qJD(5) - t398;
t683 = qJDD(1) * t665;
t682 = pkin(3) * t426 + qJDD(4);
t497 = t665 * qJD(1);
t452 = pkin(3) * t472 + qJD(4) - t497;
t568 = -qJ(5) * t586 + t452;
t402 = pkin(4) * t441 + t568;
t509 = g(3) * t522;
t539 = qJDD(2) + qJDD(3);
t451 = qJDD(2) * pkin(2) - t678 * t688;
t615 = t550 * t624;
t453 = t678 * (-t615 + t622);
t570 = qJD(3) * t585 + t553 * t451 - t549 * t453;
t370 = pkin(3) * t539 - qJ(4) * t425 + qJD(4) * t474 + t570;
t635 = qJD(3) * t549;
t680 = (qJD(3) * t482 + t453) * t553 + t549 * t451 - t492 * t635;
t375 = -qJ(4) * t426 - qJD(4) * t472 + t680;
t354 = t370 * t547 - t546 * t375;
t592 = qJDD(5) - t354;
t562 = t402 * t586 - t521 * t597 + t509 + t592;
t396 = qJDD(6) + t400;
t528 = pkin(2) * t553 + pkin(3);
t466 = t528 * t547 - t516;
t462 = -pkin(4) - t466;
t455 = -pkin(9) + t462;
t681 = -t684 * (-t670 - t643) + t455 * t396;
t679 = pkin(4) + pkin(9);
t676 = pkin(3) * t474;
t675 = pkin(3) * t535;
t399 = t425 * t546 + t547 * t426;
t674 = pkin(4) * t399;
t673 = pkin(4) * t521;
t672 = pkin(4) * t539;
t661 = qJ(5) * t522;
t629 = qJD(6) * t552;
t621 = t548 * t399 + t441 * t629 + t552 * t539;
t630 = qJD(6) * t548;
t373 = -t541 * t630 + t621;
t660 = t373 * t552;
t484 = t549 * t550 - t553 * t554;
t446 = t547 * t484 + t485 * t546;
t447 = -t484 * t546 + t485 * t547;
t599 = pkin(3) * t484 - t665;
t581 = -qJ(5) * t447 + t599;
t386 = t446 * t679 + t581;
t659 = t386 * t396;
t658 = t396 * t548;
t652 = t547 * t424;
t397 = t423 * t546 + t652;
t657 = t397 * t586;
t656 = t428 * t441;
t655 = t430 * t441;
t654 = t446 * t548;
t650 = t548 * t551;
t649 = t548 * t555;
t648 = t551 * t552;
t393 = t552 * t396;
t647 = t552 * t555;
t355 = t546 * t370 + t547 * t375;
t417 = pkin(3) * t541 + t423;
t392 = t546 * t417 + t652;
t644 = t671 - t645;
t525 = pkin(3) * t536;
t639 = t525 + t537;
t543 = t550 ^ 2;
t638 = -t554 ^ 2 + t543;
t381 = t441 * t679 + t568;
t632 = qJD(6) * t381;
t631 = qJD(6) * t541;
t628 = t671 + t627;
t532 = t550 * t664;
t620 = t525 + t594;
t523 = -pkin(3) * t547 - pkin(4);
t618 = qJD(2) * t678;
t613 = pkin(3) * t450 + t532;
t494 = -pkin(2) * t550 - t675;
t612 = t494 - t673;
t610 = t643 * t586;
t352 = -t539 * qJ(5) - t541 * qJD(5) - t355;
t350 = -pkin(5) * t399 - t352;
t391 = t417 * t547 - t418;
t595 = qJD(5) - t391;
t376 = -t541 * t679 + t595 + t671;
t361 = t376 * t548 + t381 * t552;
t609 = t350 * t552 - t361 * t441;
t491 = t550 * t618;
t493 = t554 * t618;
t578 = -t553 * t491 - t549 * t493 - t498 * t634 - t499 * t635;
t403 = -qJ(4) * t450 - qJD(4) * t484 + t578;
t449 = t541 * t484;
t569 = -qJD(3) * t640 + t491 * t549 - t553 * t493;
t404 = qJ(4) * t449 - qJD(4) * t485 + t569;
t365 = t403 * t546 - t547 * t404;
t605 = -t553 * t498 - t499 * t549;
t438 = -qJ(4) * t485 + t605;
t439 = -qJ(4) * t484 + t640;
t410 = -t547 * t438 + t439 * t546;
t390 = -qJ(5) * t541 - t392;
t378 = -t390 - t670;
t604 = t684 * t378;
t602 = t684 * t552;
t519 = pkin(2) * t615;
t469 = t519 - t683;
t571 = t469 + t682;
t561 = t571 + t685;
t351 = t399 * t679 + t561;
t600 = qJD(6) * t376 + t351;
t598 = -t673 - t675;
t591 = -t632 + t509;
t518 = -pkin(9) + t523;
t590 = -(t397 - t670) * t684 + t518 * t396;
t389 = -pkin(4) * t541 + t595;
t589 = t389 * t441 - t390 * t586;
t588 = -t391 * t441 + t392 * t586;
t366 = t403 * t547 + t404 * t546;
t411 = t438 * t546 + t439 * t547;
t360 = t376 * t552 - t381 * t548;
t584 = t350 * t548 + t360 * t441 + (t552 * t586 + t629) * t378;
t467 = pkin(2) * t651 + t528 * t546;
t583 = -0.2e1 * pkin(1) * t624 - pkin(7) * qJDD(2);
t414 = -t449 * t546 + t547 * t450;
t582 = t414 * t548 + t446 * t629;
t408 = pkin(4) * t586 + qJ(5) * t441 - t676;
t387 = pkin(5) * t447 + t410;
t577 = t350 * t446 + t378 * t414 - t387 * t396;
t576 = t519 - t686 + t682;
t531 = pkin(2) * t637;
t407 = t408 + t531;
t415 = -t449 * t547 - t450 * t546;
t575 = -qJ(5) * t415 - qJD(5) * t447 + t613;
t574 = -g(3) * t521 - t522 * t597;
t556 = qJD(2) ^ 2;
t573 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t556 + t686;
t557 = qJD(1) ^ 2;
t572 = pkin(1) * t557 - pkin(7) * qJDD(1) + t597;
t436 = t586 * pkin(9);
t566 = (-qJD(6) * t455 + t407 + t436) * t684 + t574;
t565 = (-qJD(6) * t518 + t408 + t436) * t684 + t574;
t564 = t365 * t586 - t366 * t441 - t399 * t411 + t400 * t410 - t597;
t563 = g(3) * t535 - t497 * t472 + t597 * t536 - t680;
t395 = t552 * t399;
t374 = qJD(6) * t430 + t539 * t548 - t395;
t560 = -t474 * t472 * MDP(11) + t684 * t441 * MDP(28) + ((-t374 - t691) * t552 + (-t373 + t692) * t548) * MDP(25) + (-t603 * t684 + t393 + t655) * MDP(26) + (-t602 * t684 - t656 - t658) * MDP(27) + (-t430 * t603 + t660) * MDP(24) + (t472 * t541 + t425) * MDP(13) + (-t593 + (-qJD(1) * t485 - t474) * t541) * MDP(14) + (-t472 ^ 2 + t474 ^ 2) * MDP(12) + t539 * MDP(15);
t559 = -t402 * t441 - t352 + t574;
t558 = -t497 * t474 + t570 + t693;
t540 = -qJ(4) - t678;
t520 = pkin(3) * t546 + qJ(5);
t496 = t555 * t661;
t495 = t551 * t661;
t489 = pkin(1) + t639;
t481 = t555 * t489;
t461 = qJ(5) + t467;
t460 = -t521 * t650 + t647;
t459 = t521 * t648 + t649;
t458 = t521 * t649 + t648;
t457 = t521 * t647 - t650;
t409 = pkin(4) * t446 + t581;
t388 = -pkin(5) * t446 + t411;
t367 = pkin(4) * t414 + t575;
t364 = t414 * t679 + t575;
t363 = -pkin(5) * t414 + t366;
t362 = pkin(5) * t415 + t365;
t356 = t561 + t674;
t353 = t592 - t672;
t349 = pkin(5) * t400 - t539 * t679 + t592;
t346 = t552 * t349;
t1 = [(t396 * t447 + t415 * t684) * MDP(28) + (t446 * t393 - t374 * t447 - t415 * t428 + (t414 * t552 - t446 * t630) * t684) * MDP(27) + (t373 * t447 + t396 * t654 + t415 * t430 + t582 * t684) * MDP(26) + (g(1) * t459 - g(2) * t457 - t361 * t415 + t363 * t430 + t388 * t373 + (-(qJD(6) * t387 + t364) * t684 - t659 - t600 * t447 + t378 * qJD(6) * t446) * t552 + (-(-qJD(6) * t386 + t362) * t684 - (t349 - t632) * t447 + t577) * t548) * MDP(30) + (-g(1) * t460 - g(2) * t458 + t346 * t447 + t360 * t415 + t363 * t428 + t388 * t374 + (-t351 * t447 - t364 * t684 - t659) * t548 + (t362 * t684 - t577) * t552 + ((-t386 * t552 - t387 * t548) * t684 - t361 * t447 + t378 * t654) * qJD(6)) * MDP(29) + (-t356 * t447 + t366 * t541 - t367 * t586 - t400 * t409 - t402 * t415 + t411 * t539 + t521 * t686) * MDP(22) + (-t356 * t446 + t365 * t541 - t367 * t441 - t399 * t409 - t402 * t414 + t410 * t539 - t522 * t686) * MDP(21) + t686 * MDP(2) + (-t449 * t541 + t485 * t539) * MDP(13) + (-t425 * t484 - t426 * t485 + t449 * t472 + t450 * t474) * MDP(12) + (t425 * t485 + t449 * t474) * MDP(11) + (t550 * t583 + t554 * t573) * MDP(9) + (-t550 * t573 + t554 * t583) * MDP(10) + (-t354 * t447 - t355 * t446 - t391 * t415 - t392 * t414 + t564) * MDP(18) + (t352 * t446 + t353 * t447 + t389 * t415 + t390 * t414 + t564) * MDP(20) + (-t426 * t665 - t497 * t450 + t469 * t484 + t472 * t532 + t536 * t686 + t539 * t605 + t541 * t569) * MDP(16) + (-t425 * t665 + t497 * t449 + t469 * t485 - t474 * t532 - t535 * t686 - t539 * t640 - t541 * t578) * MDP(17) + qJDD(1) * MDP(1) + (qJDD(2) * t550 + t554 * t556) * MDP(6) + (qJDD(2) * t554 - t550 * t556) * MDP(7) + (-t450 * t541 - t484 * t539) * MDP(14) + (-g(2) * t481 - t352 * t411 + t353 * t410 + t356 * t409 + t389 * t365 - t390 * t366 + t402 * t367 + (g(1) * t540 - g(2) * t594) * t555 + (-g(1) * (-t489 - t594) + g(2) * t540) * t551) * MDP(23) + t597 * MDP(3) + (t355 * t411 + t392 * t366 - t354 * t410 - t391 * t365 + t571 * t599 + t452 * t613 - g(1) * (-t489 * t551 - t540 * t555) - g(2) * (-t540 * t551 + t481)) * MDP(19) + (qJDD(1) * t543 + 0.2e1 * t550 * t614) * MDP(4) + 0.2e1 * (t550 * t622 - t624 * t638) * MDP(5) + (t373 * t654 + t430 * t582) * MDP(24) + ((-t428 * t548 + t430 * t552) * t414 + (t660 - t374 * t548 + (-t428 * t552 - t430 * t548) * qJD(6)) * t446) * MDP(25); (t641 * t541 + (t474 * t637 - t549 * t539 - t541 * t634) * pkin(2) + t563) * MDP(17) + (t407 * t441 + t643 * t541 + (-pkin(4) + t462) * t539 + t562) * MDP(21) + (-t399 * t467 - t400 * t466 - t441 * t642 + t588 + t610) * MDP(18) + (-t399 * t461 + t400 * t462 + t441 * t645 + t589 + t610) * MDP(20) + (-g(3) * t554 + t550 * t572) * MDP(9) + (g(3) * t550 + t554 * t572) * MDP(10) + (-t352 * t461 + t353 * t462 - t402 * t407 - g(1) * (t555 * t612 + t496) - g(2) * (t551 * t612 + t495) - g(3) * (t537 + t620) + t645 * t390 + t643 * t389) * MDP(23) + t560 + (t407 * t586 + t461 * t539 - t541 * t645 + t559) * MDP(22) + (-t606 * t541 + (-t472 * t637 + t539 * t553 - t541 * t635) * pkin(2) + t558) * MDP(16) + (t461 * t373 + t644 * t430 + t566 * t552 + (-t604 - t681) * t548 + t609) * MDP(30) + qJDD(2) * MDP(8) + (t355 * t467 + t354 * t466 - t452 * (t531 - t676) - g(3) * t639 - t597 * t494 + t642 * t392 - t643 * t391) * MDP(19) + (t461 * t374 + t644 * t428 + t566 * t548 + t552 * t681 + t584) * MDP(29) + MDP(6) * t623 + MDP(7) * t622 + (-MDP(4) * t550 * t554 + MDP(5) * t638) * t557; (-t397 * t541 + t408 * t441 + (-pkin(4) + t523) * t539 + t562) * MDP(21) + (t391 * t397 - t392 * t398 + (t354 * t547 + t355 * t546 + t452 * t474 + t693) * pkin(3)) * MDP(19) + (-t399 * t520 + t400 * t523 - t441 * t627 + t589 - t657) * MDP(20) + (-t657 + t398 * t441 + (-t399 * t546 - t400 * t547) * pkin(3) + t588) * MDP(18) + (-t352 * t520 + t353 * t523 - t402 * t408 - t389 * t397 - g(1) * (t555 * t598 + t496) - g(2) * (t551 * t598 + t495) - g(3) * t620 - t627 * t390) * MDP(23) + t560 + (t408 * t586 + t520 * t539 + t541 * t627 + t559) * MDP(22) + (-t541 * t585 + t558) * MDP(16) + (t541 * t607 + t563) * MDP(17) + (t520 * t374 + t428 * t628 + t548 * t565 + t552 * t590 + t584) * MDP(29) + (t520 * t373 + t628 * t430 + (-t604 - t590) * t548 + t565 * t552 + t609) * MDP(30); (t391 * t586 + t392 * t441 + t576) * MDP(19) + (-t541 * t586 - t399) * MDP(21) + (-t400 + t690) * MDP(22) + (-t389 * t586 - t390 * t441 + t576 + t674 + t685) * MDP(23) + (t656 - t658) * MDP(29) + (-t393 + t655) * MDP(30) + (-MDP(29) * t602 + MDP(30) * t603) * t684 - (MDP(19) + MDP(23)) * t683 + (MDP(18) + MDP(20)) * (-t441 ^ 2 - t687); (t400 + t690) * MDP(20) + (-t441 * t586 + t539) * MDP(21) + (-t541 ^ 2 - t687) * MDP(22) + (t390 * t541 + t562 - t672) * MDP(23) + (-t428 * t541 + t393) * MDP(29) + (-t430 * t541 - t658) * MDP(30) - (MDP(29) * t548 + MDP(30) * t552) * t684 ^ 2; t430 * t428 * MDP(24) + (-t428 ^ 2 + t430 ^ 2) * MDP(25) + (t621 + t692) * MDP(26) + (t395 + t691) * MDP(27) + t396 * MDP(28) + (-g(1) * t457 - g(2) * t459 + t361 * t684 - t378 * t430 + t346) * MDP(29) + (g(1) * t458 - g(2) * t460 + t360 * t684 + t378 * t428) * MDP(30) + (-MDP(27) * t631 + MDP(29) * t591 - MDP(30) * t600) * t552 + (-MDP(26) * t631 + (-qJD(6) * t441 - t539) * MDP(27) - t600 * MDP(29) + (-t349 - t591) * MDP(30)) * t548;];
tau  = t1;
