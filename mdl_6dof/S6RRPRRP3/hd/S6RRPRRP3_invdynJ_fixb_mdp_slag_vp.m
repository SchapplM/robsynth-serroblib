% Calculate vector of inverse dynamics joint torques for
% S6RRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RRPRRP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:50:51
% EndTime: 2019-03-09 11:51:03
% DurationCPUTime: 8.80s
% Computational Cost: add. (8840->550), mult. (20548->701), div. (0->0), fcn. (15484->14), ass. (0->250)
t601 = sin(pkin(10));
t606 = sin(qJ(2));
t677 = qJD(1) * t606;
t602 = cos(pkin(10));
t609 = cos(qJ(2));
t698 = t602 * t609;
t546 = qJD(1) * t698 - t601 * t677;
t558 = t601 * t609 + t602 * t606;
t548 = t558 * qJD(1);
t476 = pkin(2) * t677 + pkin(3) * t548 - pkin(8) * t546;
t608 = cos(qJ(4));
t469 = t608 * t476;
t720 = qJ(3) + pkin(7);
t572 = t720 * t609;
t563 = qJD(1) * t572;
t551 = t601 * t563;
t571 = t720 * t606;
t562 = qJD(1) * t571;
t503 = -t562 * t602 - t551;
t605 = sin(qJ(4));
t581 = pkin(2) * t601 + pkin(8);
t721 = pkin(9) + t581;
t654 = qJD(4) * t721;
t757 = pkin(4) * t548 - t503 * t605 + t469 + (-pkin(9) * t546 + t654) * t608;
t683 = t605 * t476 + t608 * t503;
t709 = t546 * t605;
t756 = -pkin(9) * t709 + t605 * t654 + t683;
t676 = qJD(4) * t605;
t755 = t676 - t709;
t673 = qJD(1) * qJD(2);
t658 = t606 * t673;
t573 = t601 * t658;
t657 = t609 * t673;
t636 = t602 * t657 - t573;
t614 = qJDD(1) * t558 + t636;
t752 = qJD(2) * qJD(4) + t614;
t449 = t605 * qJDD(2) - t548 * t676 + t752 * t608;
t512 = qJD(2) * t608 - t605 * t548;
t513 = qJD(2) * t605 + t548 * t608;
t604 = sin(qJ(5));
t675 = qJD(4) * t608;
t666 = t548 * t675 + t752 * t605;
t635 = t608 * qJDD(2) - t666;
t733 = cos(qJ(5));
t659 = t733 * qJD(5);
t674 = qJD(5) * t604;
t404 = -t733 * t449 - t512 * t659 + t513 * t674 - t604 * t635;
t633 = t604 * t512 + t513 * t733;
t405 = qJD(5) * t633 + t604 * t449 - t733 * t635;
t453 = -t733 * t512 + t513 * t604;
t451 = t453 ^ 2;
t547 = t558 * qJD(2);
t670 = qJDD(1) * t609;
t671 = qJDD(1) * t606;
t639 = -t601 * t671 + t602 * t670;
t495 = qJD(1) * t547 + qJDD(4) - t639;
t490 = qJDD(5) + t495;
t536 = qJD(4) - t546;
t530 = qJD(5) + t536;
t734 = t633 ^ 2;
t754 = t490 * MDP(24) + (t530 * t633 - t405) * MDP(23) + t453 * MDP(20) * t633 + (t453 * t530 - t404) * MDP(22) + (-t451 + t734) * MDP(21);
t753 = t453 * qJ(6);
t665 = t733 * t605;
t561 = t604 * t608 + t665;
t741 = qJD(4) + qJD(5);
t506 = t741 * t561;
t685 = t561 * t546 - t506;
t697 = t604 * t605;
t632 = t733 * t608 - t697;
t742 = t733 * qJD(4) + t659;
t684 = t632 * t546 - t608 * t742 + t697 * t741;
t557 = t601 * t606 - t698;
t550 = t557 * qJD(2);
t662 = t558 * t675;
t751 = -t550 * t605 + t662;
t719 = qJD(2) * pkin(2);
t554 = -t562 + t719;
t498 = t554 * t602 - t551;
t486 = -qJD(2) * pkin(3) - t498;
t450 = -pkin(4) * t512 + t486;
t597 = qJ(2) + pkin(10);
t590 = cos(t597);
t600 = qJ(4) + qJ(5);
t593 = cos(t600);
t607 = sin(qJ(1));
t701 = t593 * t607;
t592 = sin(t600);
t610 = cos(qJ(1));
t702 = t592 * t610;
t516 = -t590 * t701 + t702;
t700 = t593 * t610;
t703 = t592 * t607;
t518 = t590 * t700 + t703;
t589 = sin(t597);
t500 = -qJD(2) * t548 + t639;
t669 = pkin(2) * t658 + qJDD(3);
t722 = t609 * pkin(2);
t588 = pkin(1) + t722;
t746 = -pkin(8) * t558 - t588;
t436 = -t500 * pkin(3) - t636 * pkin(8) + qJDD(1) * t746 + t669;
t433 = t608 * t436;
t569 = -qJD(1) * t588 + qJD(3);
t465 = -pkin(3) * t546 - pkin(8) * t548 + t569;
t699 = t602 * t563;
t499 = t601 * t554 + t699;
t487 = qJD(2) * pkin(8) + t499;
t435 = t465 * t605 + t487 * t608;
t655 = qJD(2) * t720;
t545 = -qJD(3) * t606 - t609 * t655;
t494 = qJDD(2) * pkin(2) + qJD(1) * t545 - qJDD(1) * t571;
t544 = qJD(3) * t609 - t606 * t655;
t501 = qJD(1) * t544 + qJDD(1) * t572;
t445 = t601 * t494 + t602 * t501;
t443 = qJDD(2) * pkin(8) + t445;
t385 = pkin(4) * t495 - pkin(9) * t449 - qJD(4) * t435 - t443 * t605 + t433;
t629 = t605 * t436 + t608 * t443 + t465 * t675 - t487 * t676;
t388 = pkin(9) * t635 + t629;
t434 = t608 * t465 - t487 * t605;
t424 = -pkin(9) * t513 + t434;
t411 = pkin(4) * t536 + t424;
t425 = pkin(9) * t512 + t435;
t646 = -t604 * t385 - t733 * t388 - t411 * t659 + t425 * t674;
t724 = g(3) * t593;
t750 = g(1) * t518 - g(2) * t516 + t450 * t453 + t589 * t724 + t646;
t748 = qJ(6) * t633;
t416 = pkin(5) * t453 + qJD(6) + t450;
t747 = t416 * t633;
t497 = pkin(3) * t557 + t746;
t485 = t608 * t497;
t511 = -t571 * t601 + t572 * t602;
t705 = t558 * t608;
t430 = pkin(4) * t557 - pkin(9) * t705 - t511 * t605 + t485;
t504 = t608 * t511;
t682 = t605 * t497 + t504;
t706 = t558 * t605;
t439 = -pkin(9) * t706 + t682;
t686 = t604 * t430 + t733 * t439;
t502 = -t601 * t562 + t699;
t645 = pkin(4) * t755 - t502;
t555 = t721 * t605;
t556 = t721 * t608;
t681 = -t604 * t555 + t733 * t556;
t745 = g(1) * t607 - g(2) * t610;
t744 = -qJD(5) * t681 + t756 * t604 - t733 * t757;
t743 = t555 * t659 + t556 * t674 + t604 * t757 + t756 * t733;
t644 = g(1) * t610 + g(2) * t607;
t444 = t602 * t494 - t601 * t501;
t442 = -qJDD(2) * pkin(3) - t444;
t622 = g(3) * t590 - t589 * t644;
t740 = qJD(4) * t581 * t536 + t442 + t622;
t515 = t590 * t703 + t700;
t517 = -t590 * t702 + t701;
t726 = g(3) * t589;
t739 = -g(1) * t517 + g(2) * t515 + t592 * t726;
t422 = t733 * t425;
t397 = t604 * t411 + t422;
t617 = -qJD(5) * t397 + t733 * t385 - t604 * t388;
t738 = -t450 * t633 + t617 + t739;
t737 = t490 * t561 - t530 * t684;
t736 = t404 * t632 - t633 * t685;
t723 = g(3) * t609;
t717 = t449 * t605;
t716 = t453 * t548;
t715 = t633 * t548;
t713 = t512 * t536;
t712 = t512 * t548;
t711 = t513 * t536;
t710 = t513 * t548;
t707 = t550 * t608;
t565 = pkin(4) * t605 + pkin(5) * t592;
t704 = t565 * t590;
t420 = t604 * t425;
t696 = t605 * t495;
t695 = t605 * t607;
t694 = t605 * t610;
t693 = t607 * t608;
t481 = t608 * t495;
t692 = t608 * t610;
t396 = t733 * t411 - t420;
t390 = t396 - t748;
t389 = pkin(5) * t530 + t390;
t691 = -t390 + t389;
t690 = t685 * qJ(6) + qJD(6) * t632 - t743;
t689 = -pkin(5) * t548 + t684 * qJ(6) - t561 * qJD(6) + t744;
t688 = t733 * t424 - t420;
t680 = t565 + t720;
t594 = t608 * pkin(4);
t566 = pkin(5) * t593 + t594;
t598 = t606 ^ 2;
t679 = -t609 ^ 2 + t598;
t668 = t606 * t719;
t582 = -pkin(2) * t602 - pkin(3);
t663 = t558 * t676;
t652 = -t424 * t604 - t422;
t651 = t733 * t430 - t439 * t604;
t474 = t544 * t601 - t602 * t545;
t649 = -t733 * t555 - t556 * t604;
t510 = t602 * t571 + t572 * t601;
t648 = t536 * t608;
t647 = -qJD(4) * t465 - t443;
t570 = t582 - t594;
t642 = -t561 * t405 + t684 * t453;
t641 = t632 * t490 + t685 * t530;
t473 = pkin(4) * t706 + t510;
t640 = -t487 * t675 + t433;
t564 = pkin(3) + t566;
t596 = -qJ(6) - pkin(9) - pkin(8);
t638 = t564 * t590 - t589 * t596;
t441 = t751 * pkin(4) + t474;
t637 = -t536 * t755 + t481;
t634 = -0.2e1 * pkin(1) * t673 - pkin(7) * qJDD(2);
t630 = -t663 - t707;
t475 = t544 * t602 + t545 * t601;
t477 = pkin(3) * t547 + pkin(8) * t550 + t668;
t628 = t608 * t475 + t605 * t477 + t497 * t675 - t511 * t676;
t470 = t608 * t477;
t403 = pkin(9) * t707 + pkin(4) * t547 - t475 * t605 + t470 + (-t504 + (pkin(9) * t558 - t497) * t605) * qJD(4);
t407 = -pkin(9) * t751 + t628;
t627 = t604 * t403 + t733 * t407 + t430 * t659 - t439 * t674;
t625 = t486 * t536 - t581 * t495;
t623 = -qJDD(1) * t588 + t669;
t611 = qJD(2) ^ 2;
t621 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t611 + t745;
t612 = qJD(1) ^ 2;
t620 = pkin(1) * t612 - pkin(7) * qJDD(1) + t644;
t616 = -qJD(5) * t686 + t733 * t403 - t604 * t407;
t408 = -pkin(4) * t635 + t442;
t386 = t405 * pkin(5) + qJDD(6) + t408;
t587 = pkin(4) * t733 + pkin(5);
t574 = t610 * t588;
t540 = t590 * t692 + t695;
t539 = -t590 * t694 + t693;
t538 = -t590 * t693 + t694;
t537 = t590 * t695 + t692;
t483 = t632 * t558;
t482 = t561 * t558;
t462 = qJ(6) * t632 + t681;
t461 = -qJ(6) * t561 + t649;
t413 = -t550 * t665 - t604 * t663 - t674 * t706 + (-t550 * t604 + t558 * t742) * t608;
t412 = t506 * t558 + t550 * t632;
t399 = -qJ(6) * t482 + t686;
t398 = pkin(5) * t557 - qJ(6) * t483 + t651;
t393 = t688 - t748;
t392 = t652 + t753;
t391 = t397 - t753;
t382 = -qJ(6) * t413 - qJD(6) * t482 + t627;
t381 = t547 * pkin(5) + t412 * qJ(6) - t483 * qJD(6) + t616;
t380 = -qJ(6) * t405 - qJD(6) * t453 - t646;
t379 = t490 * pkin(5) + t404 * qJ(6) - qJD(6) * t633 + t617;
t1 = [(t512 * t547 - t536 * t751 + t557 * t635 - t558 * t696) * MDP(16) + (-g(1) * t516 - g(2) * t518 + t396 * t547 + t473 * t405 + t408 * t482 + t450 * t413 + t441 * t453 + t490 * t651 + t530 * t616 + t557 * t617) * MDP(25) + ((-t511 * t675 + t470) * t536 + t485 * t495 + t640 * t557 + t434 * t547 - t474 * t512 - t510 * t635 + t486 * t662 - g(1) * t538 - g(2) * t540 + ((-qJD(4) * t497 - t475) * t536 - t511 * t495 + t647 * t557 + t442 * t558 - t486 * t550) * t605) * MDP(18) + (-(t512 * t608 - t513 * t605) * t550 + (t608 * t635 - t717 + (-t512 * t605 - t513 * t608) * qJD(4)) * t558) * MDP(14) + t644 * MDP(3) + (t449 * t705 + t513 * t630) * MDP(13) + (-g(1) * t537 - g(2) * t539 - t435 * t547 + t442 * t705 + t510 * t449 + t474 * t513 + t486 * t630 - t495 * t682 - t536 * t628 - t557 * t629) * MDP(19) + (t449 * t557 + t481 * t558 + t513 * t547 + t536 * t630) * MDP(15) + 0.2e1 * (t606 * t670 - t673 * t679) * MDP(5) + (t380 * t399 + t391 * t382 + t379 * t398 + t389 * t381 + t386 * (pkin(5) * t482 + t473) + t416 * (pkin(5) * t413 + t441) - g(2) * t574 + (-g(1) * t680 - g(2) * t638) * t610 + (-g(1) * (-t588 - t638) - g(2) * t680) * t607) * MDP(28) + (qJDD(1) * t598 + 0.2e1 * t606 * t657) * MDP(4) + (t606 * t634 + t609 * t621) * MDP(9) + (-t606 * t621 + t609 * t634) * MDP(10) + (-g(1) * t515 - g(2) * t517 - t397 * t547 - t473 * t404 + t408 * t483 - t450 * t412 + t441 * t633 - t490 * t686 - t530 * t627 + t557 * t646) * MDP(26) + (-t404 * t557 - t412 * t530 + t483 * t490 + t547 * t633) * MDP(22) + (-t404 * t483 - t412 * t633) * MDP(20) + (t404 * t482 - t405 * t483 + t412 * t453 - t413 * t633) * MDP(21) + (-t444 * t558 - t445 * t557 + t474 * t548 + t475 * t546 + t498 * t550 - t499 * t547 + t511 * t500 + t510 * t614 - t644) * MDP(11) + qJDD(1) * MDP(1) + (qJDD(2) * t609 - t606 * t611) * MDP(7) + (qJDD(2) * t606 + t609 * t611) * MDP(6) + (t445 * t511 + t499 * t475 - t444 * t510 - t498 * t474 - t623 * t588 + t569 * t668 - g(1) * (-t588 * t607 + t610 * t720) - g(2) * (t607 * t720 + t574)) * MDP(12) + (-t379 * t483 - t380 * t482 - t381 * t633 - t382 * t453 + t389 * t412 - t391 * t413 + t398 * t404 - t399 * t405 + t589 * t745) * MDP(27) + t745 * MDP(2) + (t495 * t557 + t536 * t547) * MDP(17) + (t490 * t557 + t530 * t547) * MDP(24) + (-t405 * t557 - t413 * t530 - t453 * t547 - t482 * t490) * MDP(23); (g(3) * t606 + t609 * t620) * MDP(10) + (-t379 * t561 + t380 * t632 + t389 * t684 + t391 * t685 + t404 * t461 - t405 * t462 - t453 * t690 - t590 * t644 - t633 * t689 - t726) * MDP(27) + (t582 * t666 - t469 * t536 + t502 * t512 + (t503 * t536 + t625) * t605 + (-t582 * qJDD(2) - t740) * t608) * MDP(18) + (t582 * t449 - t502 * t513 + t683 * t536 + t605 * t740 + t625 * t608) * MDP(19) + (-t715 + t737) * MDP(22) + (t570 * t405 - t408 * t632 + t490 * t649 - t590 * t724 + (g(1) * t700 + g(2) * t701) * t589 + t744 * t530 + t645 * t453 - t685 * t450) * MDP(25) + (-t570 * t404 + t408 * t561 - t684 * t450 - t681 * t490 + t530 * t743 + t622 * t592 + t645 * t633) * MDP(26) + (t380 * t462 + t379 * t461 + t386 * (-pkin(5) * t632 + t570) - g(3) * (t638 + t722) + (-pkin(5) * t685 + t645) * t416 + t690 * t391 + t689 * t389 + t644 * (pkin(2) * t606 + t564 * t589 + t590 * t596)) * MDP(28) + (t536 * t648 + t696 - t710) * MDP(15) + (t637 - t712) * MDP(16) + (t641 + t716) * MDP(23) + (-t404 * t561 - t633 * t684) * MDP(20) + ((-t503 + t498) * t546 + (t601 * t500 + (-t601 * t670 + t573 + (-t657 - t671) * t602) * t602) * pkin(2)) * MDP(11) + (t642 - t736) * MDP(21) + MDP(7) * t670 + MDP(6) * t671 + qJDD(2) * MDP(8) + ((t449 + t713) * t608 + (t635 - t711) * t605) * MDP(14) + (t513 * t648 + t717) * MDP(13) + (t606 * t620 - t723) * MDP(9) + (t498 * t502 - t499 * t503 + (-t723 + t444 * t602 + t445 * t601 + (-qJD(1) * t569 + t644) * t606) * pkin(2)) * MDP(12) + (-MDP(4) * t606 * t609 + MDP(5) * t679) * t612 + (-t434 * MDP(18) + t435 * MDP(19) - t396 * MDP(25) + t397 * MDP(26) - t536 * MDP(17) - t530 * MDP(24) + (t499 - t502) * MDP(11)) * t548; (-t546 ^ 2 - t548 ^ 2) * MDP(11) + (t498 * t548 - t499 * t546 + t623 - t745) * MDP(12) + (t637 + t712) * MDP(18) + (-t536 ^ 2 * t608 - t696 - t710) * MDP(19) + (t641 - t716) * MDP(25) + (-t715 - t737) * MDP(26) + (t642 + t736) * MDP(27) + (t379 * t632 + t380 * t561 + t389 * t685 - t391 * t684 - t416 * t548 - t745) * MDP(28); -t513 * t512 * MDP(13) + (-t512 ^ 2 + t513 ^ 2) * MDP(14) + (t449 - t713) * MDP(15) + (t635 + t711) * MDP(16) + t495 * MDP(17) + (-g(1) * t539 + g(2) * t537 + t435 * t536 - t486 * t513 + (t647 + t726) * t605 + t640) * MDP(18) + (g(1) * t540 - g(2) * t538 + t434 * t536 - t486 * t512 + t608 * t726 - t629) * MDP(19) + (-t652 * t530 + (-t513 * t453 + t490 * t733 - t530 * t674) * pkin(4) + t738) * MDP(25) + (t688 * t530 + (-t604 * t490 - t513 * t633 - t530 * t659) * pkin(4) + t750) * MDP(26) + (-t389 * t453 + t391 * t633 + t392 * t633 + t393 * t453 + t587 * t404 + (-t405 * t604 + (-t453 * t733 + t604 * t633) * qJD(5)) * pkin(4)) * MDP(27) + (t379 * t587 - t391 * t393 - t389 * t392 - pkin(5) * t747 - g(1) * (t566 * t607 - t610 * t704) - g(2) * (-t566 * t610 - t607 * t704) + t565 * t726 + (t380 * t604 - t416 * t513 + (-t389 * t604 + t391 * t733) * qJD(5)) * pkin(4)) * MDP(28) + t754; (t397 * t530 + t738) * MDP(25) + (t396 * t530 + t750) * MDP(26) + (pkin(5) * t404 - t453 * t691) * MDP(27) + (t691 * t391 + (t379 + t739 - t747) * pkin(5)) * MDP(28) + t754; (-t451 - t734) * MDP(27) + (t389 * t633 + t391 * t453 + t386 + t622) * MDP(28);];
tau  = t1;
