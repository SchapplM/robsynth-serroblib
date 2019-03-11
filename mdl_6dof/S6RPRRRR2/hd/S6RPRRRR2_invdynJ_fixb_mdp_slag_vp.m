% Calculate vector of inverse dynamics joint torques for
% S6RPRRRR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRRR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:58:47
% EndTime: 2019-03-09 06:58:58
% DurationCPUTime: 8.02s
% Computational Cost: add. (6813->521), mult. (14193->679), div. (0->0), fcn. (10544->18), ass. (0->247)
t599 = qJD(3) + qJD(4);
t609 = sin(qJ(4));
t614 = cos(qJ(3));
t737 = cos(qJ(4));
t677 = qJD(1) * t737;
t610 = sin(qJ(3));
t698 = qJD(1) * t610;
t752 = -t609 * t698 + t614 * t677;
t755 = t599 * t752;
t604 = qJ(3) + qJ(4);
t594 = sin(t604);
t600 = qJ(1) + pkin(11);
t589 = sin(t600);
t590 = cos(t600);
t745 = g(1) * t590 + g(2) * t589;
t754 = t745 * t594;
t739 = qJD(5) + qJD(6);
t753 = -t752 + t739;
t605 = sin(pkin(11));
t580 = pkin(1) * t605 + pkin(7);
t732 = pkin(8) + t580;
t670 = t732 * qJD(1);
t516 = qJD(2) * t610 + t614 * t670;
t503 = t609 * t516;
t515 = t614 * qJD(2) - t670 * t610;
t462 = t515 * t737 - t503;
t676 = qJD(4) * t737;
t746 = -pkin(3) * t676 + t462;
t711 = t609 * t614;
t534 = -qJD(1) * t711 - t610 * t677;
t491 = -pkin(4) * t534 - pkin(9) * t752;
t482 = pkin(3) * t698 + t491;
t608 = sin(qJ(5));
t613 = cos(qJ(5));
t751 = -t613 * t482 + t608 * t746;
t506 = -t534 * t608 - t613 * t599;
t612 = cos(qJ(6));
t607 = sin(qJ(6));
t639 = t534 * t613 - t599 * t608;
t721 = t639 * t607;
t449 = t612 * t506 - t721;
t529 = qJD(5) - t752;
t525 = qJD(6) + t529;
t750 = t449 * t525;
t640 = t506 * t607 + t612 * t639;
t749 = t525 * t640;
t713 = t607 * t613;
t544 = t608 * t612 + t713;
t748 = t753 * t544;
t542 = t607 * t608 - t612 * t613;
t703 = t753 * t542;
t504 = t737 * t516;
t730 = qJD(3) * pkin(3);
t505 = t515 + t730;
t456 = t609 * t505 + t504;
t446 = pkin(9) * t599 + t456;
t606 = cos(pkin(11));
t581 = -pkin(1) * t606 - pkin(2);
t556 = -pkin(3) * t614 + t581;
t535 = t556 * qJD(1);
t474 = -pkin(4) * t752 + pkin(9) * t534 + t535;
t421 = t446 * t613 + t474 * t608;
t409 = -pkin(10) * t506 + t421;
t693 = qJD(6) * t607;
t407 = t409 * t693;
t455 = t505 * t737 - t503;
t445 = -t599 * pkin(4) - t455;
t429 = t506 * pkin(5) + t445;
t603 = qJ(5) + qJ(6);
t593 = sin(t603);
t595 = cos(t603);
t596 = cos(t604);
t716 = t595 * t596;
t510 = -t589 * t716 + t590 * t593;
t512 = t589 * t593 + t590 * t716;
t582 = g(3) * t594;
t744 = g(1) * t512 - g(2) * t510 + t429 * t449 + t595 * t582 + t407;
t717 = t593 * t596;
t509 = t589 * t717 + t590 * t595;
t511 = t589 * t595 - t590 * t717;
t598 = qJDD(3) + qJDD(4);
t591 = t614 * qJDD(2);
t560 = t580 * qJDD(1);
t669 = pkin(8) * qJDD(1) + t560;
t460 = qJDD(3) * pkin(3) - qJD(3) * t516 - t669 * t610 + t591;
t469 = qJD(3) * t515 + t610 * qJDD(2) + t669 * t614;
t696 = qJD(4) * t609;
t623 = t609 * t460 + t469 * t737 + t505 * t676 - t516 * t696;
t398 = t598 * pkin(9) + t623;
t671 = qJDD(1) * t737;
t689 = qJDD(1) * t614;
t471 = t609 * t689 + t610 * t671 + t755;
t545 = t610 * t737 + t711;
t498 = t599 * t545;
t690 = qJDD(1) * t610;
t648 = t609 * t690 - t614 * t671;
t472 = qJD(1) * t498 + t648;
t691 = qJD(1) * qJD(3);
t675 = t610 * t691;
t517 = pkin(3) * t675 + qJDD(1) * t556;
t415 = pkin(4) * t472 - pkin(9) * t471 + t517;
t414 = t613 * t415;
t694 = qJD(5) * t613;
t695 = qJD(5) * t608;
t431 = t613 * t471 + t534 * t695 + t608 * t598 + t599 * t694;
t468 = qJDD(5) + t472;
t380 = pkin(5) * t468 - pkin(10) * t431 - qJD(5) * t421 - t398 * t608 + t414;
t432 = -qJD(5) * t639 + t471 * t608 - t613 * t598;
t630 = t613 * t398 + t608 * t415 - t446 * t695 + t474 * t694;
t381 = -pkin(10) * t432 + t630;
t668 = t612 * t380 - t607 * t381;
t743 = -g(1) * t511 + g(2) * t509 + t429 * t640 + t593 * t582 + t668;
t464 = qJDD(6) + t468;
t742 = t464 * MDP(30) + (-t449 ^ 2 + t640 ^ 2) * MDP(27) - t449 * MDP(26) * t640;
t487 = t544 * t545;
t461 = t609 * t515 + t504;
t659 = pkin(3) * t696 - t461;
t720 = t752 * t608;
t741 = (t695 - t720) * pkin(5);
t740 = t608 * t482 + t613 * t746;
t667 = t431 * t607 + t612 * t432;
t395 = -qJD(6) * t640 + t667;
t738 = -pkin(9) - pkin(10);
t734 = g(3) * t596;
t733 = t613 * pkin(5);
t597 = t613 * pkin(10);
t586 = pkin(3) * t609 + pkin(9);
t731 = -pkin(10) - t586;
t420 = -t446 * t608 + t613 * t474;
t408 = pkin(10) * t639 + t420;
t402 = pkin(5) * t529 + t408;
t729 = t402 * t612;
t728 = t409 * t612;
t727 = t431 * t608;
t726 = t445 * t752;
t634 = -t609 * t610 + t614 * t737;
t497 = t599 * t634;
t725 = t497 * t608;
t724 = t497 * t613;
t723 = t506 * t529;
t722 = t639 * t529;
t719 = t545 * t608;
t718 = t545 * t613;
t715 = t596 * t608;
t714 = t596 * t613;
t712 = t608 * t468;
t536 = t732 * t610;
t537 = t732 * t614;
t490 = -t609 * t536 + t537 * t737;
t483 = t613 * t490;
t710 = qJDD(2) - g(3);
t692 = qJD(6) * t612;
t681 = t612 * t431 - t607 * t432 - t506 * t692;
t394 = t639 * t693 + t681;
t709 = -t394 * t634 - t498 * t640;
t678 = t545 * t695;
t406 = t497 * t713 - t607 * t678 - t693 * t719 + (t718 * t739 + t725) * t612;
t708 = -t406 * t525 - t487 * t464;
t707 = -t431 * t634 - t498 * t639;
t706 = t613 * t455 + t608 * t491;
t484 = -pkin(4) * t634 - pkin(9) * t545 + t556;
t702 = t608 * t484 + t483;
t701 = t741 + t659;
t601 = t610 ^ 2;
t700 = -t614 ^ 2 + t601;
t563 = qJD(1) * t581;
t688 = pkin(10) * t720;
t687 = t610 * t730;
t684 = qJD(5) * pkin(9) * t529;
t683 = t545 * t712;
t682 = t468 * t718;
t680 = qJD(5) * t738;
t444 = t445 * t694;
t661 = -t737 * t460 + t609 * t469 + t505 * t696 + t516 * t676;
t399 = -pkin(4) * t598 + t661;
t674 = -t399 - t734;
t673 = qJD(3) * t732;
t672 = qJD(5) * t731;
t665 = t529 * t613;
t664 = -qJD(5) * t474 - t398;
t663 = qJD(6) * t402 + t381;
t587 = -pkin(3) * t737 - pkin(4);
t660 = -t534 * pkin(5) - t597 * t752;
t658 = -t456 + t741;
t656 = g(1) * t589 - g(2) * t590;
t611 = sin(qJ(1));
t615 = cos(qJ(1));
t655 = g(1) * t611 - g(2) * t615;
t654 = -t446 * t694 + t414;
t539 = t586 * t613 + t597;
t653 = qJD(6) * t539 - t613 * t672 + t660 - t751;
t486 = t613 * t491;
t565 = pkin(9) * t613 + t597;
t652 = qJD(6) * t565 - t455 * t608 - t613 * t680 + t486 + t660;
t538 = t731 * t608;
t651 = -qJD(6) * t538 - t608 * t672 - t688 + t740;
t564 = t738 * t608;
t650 = -qJD(6) * t564 - t608 * t680 - t688 + t706;
t649 = -pkin(9) * t468 - t726;
t645 = t395 * t634 - t449 * t498;
t386 = t402 * t607 + t728;
t405 = -t487 * t739 - t542 * t497;
t488 = t542 * t545;
t644 = -t405 * t525 + t464 * t488;
t643 = t432 * t634 - t498 * t506;
t642 = -t468 * t586 - t726;
t641 = t497 * t599 + t545 * t598;
t638 = g(3) * t715 + t399 * t608 - t421 * t534 + t444;
t637 = t420 * t534 + t445 * t695 + t613 * t754;
t635 = -t536 * t737 - t609 * t537;
t633 = t545 * t694 + t725;
t632 = t678 - t724;
t526 = t610 * t673;
t527 = t614 * t673;
t434 = qJD(4) * t635 - t526 * t737 - t609 * t527;
t442 = pkin(4) * t498 - pkin(9) * t497 + t687;
t629 = t613 * t434 + t608 * t442 + t484 * t694 - t490 * t695;
t628 = -qJD(1) * t563 - t560 + t745;
t627 = 0.2e1 * qJD(3) * t563 - qJDD(3) * t580;
t625 = t535 * t534 - t661 - t734 + t754;
t616 = qJD(3) ^ 2;
t624 = -0.2e1 * qJDD(1) * t581 - t580 * t616 + t656;
t385 = -t409 * t607 + t729;
t390 = pkin(5) * t432 + t399;
t622 = -g(3) * t716 + t385 * t534 + t390 * t542 + t429 * t748 + t595 * t754;
t435 = qJD(4) * t490 - t609 * t526 + t527 * t737;
t621 = g(3) * t717 - t386 * t534 + t390 * t544 - t429 * t703 - t593 * t754;
t620 = -t535 * t752 + t596 * t745 + t582 - t623;
t619 = (-t394 * t542 - t395 * t544 + t449 * t703 + t640 * t748) * MDP(27) + (t394 * t544 + t640 * t703) * MDP(26) + ((t431 - t723) * t613 + (-t432 + t722) * t608) * MDP(20) + (t464 * t544 - t525 * t703 - t534 * t640) * MDP(28) + (-t449 * t534 - t464 * t542 - t525 * t748) * MDP(29) + (-t639 * t665 + t727) * MDP(19) + (-t529 ^ 2 * t608 + t468 * t613 - t506 * t534) * MDP(22) + (t529 * t665 - t534 * t639 + t712) * MDP(21) + (t471 - t755) * MDP(14) + (-t648 + (-qJD(1) * t545 - t534) * t599) * MDP(15) + (t534 ^ 2 - t752 ^ 2) * MDP(13) + t598 * MDP(16) + (MDP(12) * t752 + MDP(23) * t529 + MDP(30) * t525) * t534;
t588 = -pkin(4) - t733;
t559 = t587 - t733;
t558 = qJDD(3) * t614 - t610 * t616;
t557 = qJDD(3) * t610 + t614 * t616;
t521 = t589 * t608 + t590 * t714;
t520 = t589 * t613 - t590 * t715;
t519 = -t589 * t714 + t590 * t608;
t518 = t589 * t715 + t590 * t613;
t481 = t613 * t484;
t470 = -t498 * t599 + t598 * t634;
t466 = pkin(5) * t719 - t635;
t439 = t613 * t442;
t424 = -pkin(10) * t719 + t702;
t422 = -pkin(5) * t634 - pkin(10) * t718 - t490 * t608 + t481;
t412 = pkin(5) * t633 + t435;
t387 = -pkin(10) * t633 + t629;
t384 = -pkin(10) * t724 + pkin(5) * t498 - t434 * t608 + t439 + (-t483 + (pkin(10) * t545 - t484) * t608) * qJD(5);
t1 = [(t471 * t634 - t472 * t545 + t497 * t752 + t498 * t534) * MDP(13) + (-t435 * t599 + t472 * t556 + t498 * t535 - t517 * t634 + t596 * t656 + t598 * t635 - t687 * t752) * MDP(17) + t655 * MDP(2) + (-t394 * t487 + t395 * t488 - t405 * t449 + t406 * t640) * MDP(27) + (-t394 * t488 - t405 * t640) * MDP(26) + ((-t506 * t613 + t608 * t639) * t497 + (-t727 - t432 * t613 + (t506 * t608 + t613 * t639) * qJD(5)) * t545) * MDP(20) + (t431 * t718 + t632 * t639) * MDP(19) + (g(1) * t615 + g(2) * t611) * MDP(3) + ((t384 * t612 - t387 * t607) * t525 + (t422 * t612 - t424 * t607) * t464 - t668 * t634 + t385 * t498 + t412 * t449 + t466 * t395 + t390 * t487 + t429 * t406 - g(1) * t510 - g(2) * t512 + ((-t422 * t607 - t424 * t612) * t525 + t386 * t634) * qJD(6)) * MDP(31) + (-t468 * t634 + t498 * t529) * MDP(23) + (-t464 * t634 + t498 * t525) * MDP(30) + (-g(1) * t509 - g(2) * t511 - t386 * t498 - t390 * t488 + t466 * t394 + t429 * t405 - t407 * t634 - t412 * t640 + (-(-qJD(6) * t424 + t384) * t525 - t422 * t464 + t380 * t634) * t607 + (-(qJD(6) * t422 + t387) * t525 - t424 * t464 + t663 * t634) * t612) * MDP(32) + ((-t490 * t694 + t439) * t529 + t481 * t468 - t654 * t634 + t420 * t498 + t435 * t506 - t635 * t432 + t545 * t444 - g(1) * t519 - g(2) * t521 + ((-qJD(5) * t484 - t434) * t529 - t490 * t468 - t664 * t634 + t399 * t545 + t445 * t497) * t608) * MDP(24) + (-g(1) * t518 - g(2) * t520 + t399 * t718 - t421 * t498 - t431 * t635 - t435 * t639 - t445 * t632 - t468 * t702 - t529 * t629 + t630 * t634) * MDP(25) + (-t529 * t633 + t643 - t683) * MDP(22) + (-t529 * t632 + t682 + t707) * MDP(21) + (t645 + t708) * MDP(29) + (-t644 + t709) * MDP(28) + qJDD(1) * MDP(1) + t470 * MDP(15) + (-t434 * t599 + t471 * t556 - t490 * t598 + t497 * t535 + t517 * t545 - t534 * t687 - t594 * t656) * MDP(18) + t557 * MDP(7) + t558 * MDP(8) + (t471 * t545 - t497 * t534) * MDP(12) + (t610 * t627 + t614 * t624) * MDP(10) + (-t610 * t624 + t614 * t627) * MDP(11) + 0.2e1 * (t610 * t689 - t691 * t700) * MDP(6) + t641 * MDP(14) + (t655 + (t605 ^ 2 + t606 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (qJDD(1) * t601 + 0.2e1 * t614 * t675) * MDP(5); t710 * MDP(4) + t558 * MDP(10) - t557 * MDP(11) + t470 * MDP(17) - t641 * MDP(18) + (-t643 - t683) * MDP(24) + (-t682 + t707) * MDP(25) + (-t645 + t708) * MDP(31) + (t644 + t709) * MDP(32) + (-MDP(24) * t633 + MDP(25) * t632) * t529; t619 + (-(t538 * t607 + t539 * t612) * t464 + t559 * t394 + (t607 * t653 + t612 * t651) * t525 - t701 * t640 + t621) * MDP(32) + MDP(7) * t690 + MDP(8) * t689 + ((t538 * t612 - t539 * t607) * t464 + t559 * t395 + (t607 * t651 - t612 * t653) * t525 + t701 * t449 + t622) * MDP(31) + (t462 * t599 + (t534 * t698 - t598 * t609 - t599 * t676) * pkin(3) + t620) * MDP(18) + qJDD(3) * MDP(9) + (t461 * t599 + (t598 * t737 - t599 * t696 + t698 * t752) * pkin(3) + t625) * MDP(17) + (t587 * t432 + t674 * t613 + t642 * t608 + t659 * t506 + (-t586 * t694 + t751) * t529 + t637) * MDP(24) + (t587 * t431 + t642 * t613 - t608 * t754 - t659 * t639 + (t586 * t695 + t740) * t529 + t638) * MDP(25) + (-g(3) * t614 + t610 * t628 + t591) * MDP(10) + (-t610 * t710 + t614 * t628) * MDP(11) + (-MDP(5) * t610 * t614 + MDP(6) * t700) * qJD(1) ^ 2; t619 + (t455 * t599 + t620) * MDP(18) + ((t564 * t612 - t565 * t607) * t464 + t588 * t395 + (t607 * t650 - t612 * t652) * t525 + t658 * t449 + t622) * MDP(31) + (-pkin(4) * t432 - t456 * t506 - t486 * t529 + (t455 * t529 + t649) * t608 + (t674 - t684) * t613 + t637) * MDP(24) + (t456 * t599 + t625) * MDP(17) + (-pkin(4) * t431 + t706 * t529 + t456 * t639 + t649 * t613 + (-t754 + t684) * t608 + t638) * MDP(25) + (-(t564 * t607 + t565 * t612) * t464 + t588 * t394 + (t607 * t652 + t612 * t650) * t525 - t658 * t640 + t621) * MDP(32); -t639 * t506 * MDP(19) + (-t506 ^ 2 + t639 ^ 2) * MDP(20) + (t431 + t723) * MDP(21) + (-t432 - t722) * MDP(22) + t468 * MDP(23) + (-g(1) * t520 + g(2) * t518 + t421 * t529 + t445 * t639 + (t664 + t582) * t608 + t654) * MDP(24) + (g(1) * t521 - g(2) * t519 + t420 * t529 + t445 * t506 + t582 * t613 - t630) * MDP(25) + (t394 + t750) * MDP(28) + (-t395 - t749) * MDP(29) + (-(-t408 * t607 - t728) * t525 - t386 * qJD(6) + (t449 * t639 + t464 * t612 - t525 * t693) * pkin(5) + t743) * MDP(31) + ((-t409 * t525 - t380) * t607 + (t408 * t525 - t663) * t612 + (-t464 * t607 - t525 * t692 - t639 * t640) * pkin(5) + t744) * MDP(32) + t742; (t681 + t750) * MDP(28) + (-t667 - t749) * MDP(29) + (t386 * t525 + t743) * MDP(31) + (-t607 * t380 - t612 * t381 + t385 * t525 + t744) * MDP(32) + (MDP(28) * t721 + MDP(29) * t640 - MDP(31) * t386 - MDP(32) * t729) * qJD(6) + t742;];
tau  = t1;
