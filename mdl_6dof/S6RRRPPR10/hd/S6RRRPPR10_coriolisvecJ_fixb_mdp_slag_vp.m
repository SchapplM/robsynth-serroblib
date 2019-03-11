% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR10_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPPR10_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR10_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRPPR10_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:28:14
% EndTime: 2019-03-09 16:28:34
% DurationCPUTime: 12.38s
% Computational Cost: add. (8400->644), mult. (22026->859), div. (0->0), fcn. (16808->10), ass. (0->250)
t633 = sin(pkin(6));
t641 = cos(qJ(2));
t724 = qJD(1) * t641;
t701 = t633 * t724;
t611 = -qJD(3) + t701;
t638 = sin(qJ(2));
t759 = cos(pkin(6));
t693 = t759 * qJD(1);
t679 = pkin(1) * t693;
t618 = t638 * t679;
t579 = pkin(8) * t701 + t618;
t637 = sin(qJ(3));
t682 = t637 * t701;
t720 = qJD(3) * t637;
t778 = qJD(4) * t637 + t579 + (t682 - t720) * pkin(3);
t640 = cos(qJ(3));
t659 = t693 + qJD(2);
t648 = t640 * t659;
t725 = qJD(1) * t633;
t702 = t638 * t725;
t555 = t637 * t702 - t648;
t632 = sin(pkin(11));
t634 = cos(pkin(11));
t512 = t634 * t555 + t611 * t632;
t639 = cos(qJ(6));
t513 = t555 * t632 - t611 * t634;
t636 = sin(qJ(6));
t752 = t513 * t636;
t450 = -t512 * t639 + t752;
t557 = t637 * t659 + t640 * t702;
t550 = qJD(6) + t557;
t777 = t450 * t550;
t576 = -pkin(8) * t702 + t641 * t679;
t655 = t633 * (pkin(2) * t638 - pkin(9) * t641);
t577 = qJD(1) * t655;
t689 = -t637 * t576 + t577 * t640;
t761 = pkin(3) + qJ(5);
t697 = t761 * t638;
t719 = qJD(3) * t640;
t739 = t640 * t641;
t763 = pkin(4) + pkin(9);
t776 = (pkin(4) * t739 - t697) * t725 - t689 - t763 * t719;
t775 = qJD(5) * t640 + t778 + t611 * (-qJ(4) * t640 + qJ(5) * t637);
t591 = t632 * t639 + t634 * t636;
t584 = t591 * qJD(6);
t728 = -t591 * t557 - t584;
t712 = qJD(1) * qJD(2);
t695 = t633 * t712;
t677 = t641 * t695;
t744 = t633 * t638;
t707 = t637 * t744;
t680 = qJD(3) * t707;
t516 = qJD(1) * t680 - qJD(3) * t648 - t640 * t677;
t661 = t632 * t636 - t634 * t639;
t717 = qJD(6) * t639;
t718 = qJD(6) * t636;
t727 = -t557 * t661 - t632 * t718 + t634 * t717;
t774 = t591 * t516 - t550 * t727;
t663 = -t512 * t636 - t513 * t639;
t772 = t550 * t663;
t540 = t632 * t702 - t634 * t682;
t770 = -t634 * t720 - t540;
t538 = pkin(9) * t659 + t579;
t571 = (-pkin(2) * t641 - pkin(9) * t638 - pkin(1)) * t633;
t549 = qJD(1) * t571;
t475 = t538 * t637 - t640 * t549;
t716 = -qJD(4) - t475;
t629 = t633 ^ 2;
t769 = -0.2e1 * t629 * t712;
t768 = MDP(5) * (t638 ^ 2 - t641 ^ 2);
t738 = -t632 * t775 + t634 * t776;
t737 = t632 * t776 + t634 * t775;
t678 = t638 * t695;
t671 = pkin(3) * t678;
t578 = qJD(2) * t655;
t565 = qJD(1) * t578;
t703 = pkin(1) * t759;
t765 = -pkin(8) * t744 + t641 * t703;
t580 = t765 * qJD(2);
t566 = qJD(1) * t580;
t685 = t538 * t719 + t549 * t720 - t640 * t565 + t637 * t566;
t437 = -t671 + t685;
t476 = t640 * t538 + t637 * t549;
t603 = t611 * qJ(4);
t470 = t603 - t476;
t766 = -t470 * t611 + t437;
t721 = qJD(3) * t557;
t517 = t637 * t677 + t721;
t764 = t516 * t661 + t550 * t728;
t714 = pkin(4) * t557 - t716;
t554 = t557 ^ 2;
t642 = qJD(1) ^ 2;
t762 = pkin(10) * t557;
t760 = -pkin(10) - t761;
t758 = qJ(4) * t555;
t607 = qJ(4) * t678;
t686 = t538 * t720 - t549 * t719 - t637 * t565 - t640 * t566;
t670 = -qJD(4) * t611 - t686;
t435 = -t607 - t670;
t418 = -pkin(4) * t517 - t435;
t757 = t418 * t632;
t756 = t418 * t634;
t537 = -pkin(2) * t659 - t576;
t643 = -t557 * qJ(4) + t537;
t464 = t555 * pkin(3) + t643;
t755 = t464 * t557;
t751 = t516 * t761;
t750 = t555 * t557;
t749 = t555 * t611;
t748 = t557 * t611;
t653 = t611 * t640;
t746 = t629 * t642;
t743 = t633 * t641;
t741 = t634 * t640;
t740 = t637 * t641;
t567 = pkin(8) * t677 + qJD(2) * t618;
t647 = qJ(4) * t516 - qJD(4) * t557 + t567;
t415 = qJD(5) * t555 + t517 * t761 + t647;
t676 = qJD(2) * t697;
t416 = -pkin(4) * t516 + qJD(5) * t611 - t676 * t725 + t685;
t396 = t634 * t415 + t632 * t416;
t698 = qJD(2) * t743;
t529 = -t680 + (qJD(3) * t759 + t698) * t640;
t684 = t638 * t703;
t570 = pkin(8) * t743 + pkin(9) * t759 + t684;
t657 = -t570 * t719 - t571 * t720 + t578 * t640 - t637 * t580;
t429 = pkin(4) * t529 + (qJD(5) * t641 - t676) * t633 - t657;
t587 = t637 * t759 + t640 * t744;
t528 = qJD(3) * t587 + t637 * t698;
t586 = -t640 * t759 + t707;
t581 = pkin(8) * t698 + qJD(2) * t684;
t646 = -qJ(4) * t529 - qJD(4) * t587 + t581;
t430 = qJD(5) * t586 + t528 * t761 + t646;
t399 = t632 * t429 + t634 * t430;
t442 = t611 * t761 + t714;
t447 = t555 * t761 + t643;
t404 = t632 * t442 + t634 * t447;
t690 = -t637 * t570 + t571 * t640;
t488 = pkin(3) * t743 - t690;
t459 = pkin(4) * t587 + qJ(5) * t743 + t488;
t569 = -pkin(2) * t759 - t765;
t644 = -t587 * qJ(4) + t569;
t465 = t586 * t761 + t644;
t423 = t632 * t459 + t634 * t465;
t458 = -pkin(4) * t555 + t476;
t468 = t557 * t761 + t758;
t424 = t632 * t458 + t634 * t468;
t730 = t640 * t576 + t637 * t577;
t493 = -qJ(4) * t702 - t730;
t473 = -pkin(4) * t682 - t493;
t736 = -t763 * t720 - t473;
t541 = (t632 * t740 + t634 * t638) * t725;
t573 = t661 * t640;
t735 = qJD(6) * t573 + t540 * t636 - t541 * t639 + t591 * t720;
t734 = -t639 * t540 - t541 * t636 - t584 * t640 + t661 * t720;
t681 = t640 * t701;
t733 = (-t681 + t719) * qJ(4) + t778;
t731 = t640 * t570 + t637 * t571;
t704 = -pkin(5) * t634 - pkin(4);
t729 = -pkin(5) * t540 - t473 + (-pkin(9) + t704) * t720;
t694 = -qJ(4) * t637 - pkin(2);
t590 = -t640 * t761 + t694;
t612 = t763 * t637;
t521 = t634 * t590 + t632 * t612;
t613 = t763 * t640;
t723 = qJD(2) * t638;
t722 = qJD(2) * t640;
t715 = -t557 * t704 - t716;
t448 = qJD(5) + t458 - t603;
t713 = qJD(5) - t448;
t711 = pkin(9) * t611 * t637;
t710 = pkin(9) * t653;
t709 = pkin(9) * t723;
t708 = pkin(9) * t722;
t490 = -t634 * t517 + t632 * t678;
t491 = t517 * t632 + t634 * t678;
t706 = -t636 * t490 + t639 * t491 + t512 * t717;
t700 = t633 * t723;
t395 = -t415 * t632 + t634 * t416;
t392 = -pkin(5) * t516 - pkin(10) * t491 + t395;
t393 = -pkin(10) * t490 + t396;
t692 = t639 * t392 - t393 * t636;
t398 = t634 * t429 - t430 * t632;
t403 = t634 * t442 - t447 * t632;
t421 = t634 * t459 - t465 * t632;
t691 = t639 * t490 + t636 * t491;
t687 = t629 * t638 * t641 * MDP(4);
t674 = pkin(1) * t769;
t515 = -pkin(10) * t741 + t521;
t673 = pkin(5) * t681 - pkin(10) * t541 + qJD(6) * t515 - (-pkin(10) * t632 * t637 + pkin(5) * t640) * qJD(3) + t738;
t595 = t634 * t612;
t503 = pkin(5) * t637 + t595 + (pkin(10) * t640 - t590) * t632;
t672 = pkin(10) * t770 - qJD(6) * t503 + t737;
t668 = t392 * t636 + t393 * t639;
t667 = t395 * t634 + t396 * t632;
t400 = pkin(5) * t557 - pkin(10) * t513 + t403;
t401 = pkin(10) * t512 + t404;
t389 = t400 * t639 - t401 * t636;
t390 = t400 * t636 + t401 * t639;
t666 = -t403 * t632 + t404 * t634;
t527 = -t586 * t632 + t634 * t743;
t406 = pkin(5) * t587 + pkin(10) * t527 + t421;
t526 = t586 * t634 + t632 * t743;
t409 = pkin(10) * t526 + t423;
t665 = t406 * t639 - t409 * t636;
t664 = t406 * t636 + t409 * t639;
t662 = t639 * t526 + t527 * t636;
t472 = t526 * t636 - t527 * t639;
t487 = qJ(4) * t743 - t731;
t454 = t634 * t458;
t600 = t760 * t632;
t652 = qJD(5) * t634 + qJD(6) * t600 - pkin(5) * t555 + t454 + (-t468 - t762) * t632;
t601 = t760 * t634;
t651 = qJD(5) * t632 - qJD(6) * t601 + t634 * t762 + t424;
t650 = -t476 * t611 - t685;
t649 = -t570 * t720 + t571 * t719 + t637 * t578 + t640 * t580;
t407 = -t513 * t718 + t706;
t467 = -pkin(4) * t586 - t487;
t645 = -t516 - t749;
t408 = -qJD(6) * t663 + t691;
t440 = -qJ(4) * t700 + qJD(4) * t743 - t649;
t431 = -pkin(4) * t528 - t440;
t624 = pkin(5) * t632 + qJ(4);
t606 = -pkin(3) * t640 + t694;
t583 = pkin(5) * t741 + t613;
t574 = t591 * t640;
t520 = -t590 * t632 + t595;
t514 = t516 * t637;
t501 = t528 * t632 + t634 * t700;
t500 = -t528 * t634 + t632 * t700;
t499 = pkin(3) * t557 + t758;
t494 = -pkin(3) * t702 - t689;
t489 = t516 * t587;
t486 = t586 * pkin(3) + t644;
t469 = pkin(3) * t611 - t716;
t446 = pkin(3) * t528 + t646;
t444 = -pkin(3) * t700 - t657;
t443 = -pkin(5) * t526 + t467;
t436 = pkin(3) * t517 + t647;
t432 = -pkin(5) * t512 + t448;
t422 = -t468 * t632 + t454;
t420 = qJD(6) * t472 + t639 * t500 + t501 * t636;
t419 = qJD(6) * t662 - t500 * t636 + t501 * t639;
t412 = pkin(5) * t500 + t431;
t402 = pkin(5) * t490 + t418;
t397 = -pkin(10) * t500 + t399;
t394 = pkin(5) * t529 - pkin(10) * t501 + t398;
t388 = -qJD(6) * t390 + t692;
t387 = qJD(6) * t389 + t668;
t1 = [(t649 * t611 + t581 * t557 - t569 * t516 + t567 * t587 + t537 * t529 + (-t686 * t641 + (-qJD(1) * t731 - t476) * t723) * t633) * MDP(17) + (MDP(6) * t698 - MDP(7) * t700) * (0.2e1 * t693 + qJD(2)) + (-(qJD(6) * t665 + t394 * t636 + t397 * t639) * t550 + t664 * t516 - t387 * t587 - t390 * t529 - t412 * t663 + t443 * t407 + t402 * t472 + t432 * t419) * MDP(32) + (t407 * t587 + t419 * t550 - t472 * t516 - t529 * t663) * MDP(28) + (t407 * t472 - t419 * t663) * MDP(26) + (t407 * t662 - t408 * t472 - t419 * t450 + t420 * t663) * MDP(27) + ((-qJD(6) * t664 + t394 * t639 - t397 * t636) * t550 - t665 * t516 + t388 * t587 + t389 * t529 + t412 * t450 + t443 * t408 - t402 * t662 + t432 * t420) * MDP(31) + (-t408 * t587 - t420 * t550 - t450 * t529 - t516 * t662) * MDP(29) + (t395 * t587 + t398 * t557 + t403 * t529 - t418 * t526 - t421 * t516 - t431 * t512 + t448 * t500 + t467 * t490) * MDP(22) + (t395 * t527 + t396 * t526 - t398 * t513 + t399 * t512 - t403 * t501 - t404 * t500 - t421 * t491 - t423 * t490) * MDP(24) + (t529 * t557 - t489) * MDP(11) + (t529 * t550 - t489) * MDP(30) + (-t566 * t759 - t580 * t659 + t641 * t674) * MDP(10) + (-t567 * t759 - t581 * t659 + t638 * t674) * MDP(9) + t768 * t769 + (-t396 * t587 - t399 * t557 - t404 * t529 - t418 * t527 + t423 * t516 + t431 * t513 + t448 * t501 + t467 * t491) * MDP(23) + (t435 * t586 + t437 * t587 + t440 * t555 + t444 * t557 + t469 * t529 + t470 * t528 + t487 * t517 - t488 * t516) * MDP(18) + (t516 * t586 - t517 * t587 - t528 * t557 - t529 * t555) * MDP(12) + (t435 * t487 + t436 * t486 + t437 * t488 + t440 * t470 + t444 * t469 + t446 * t464) * MDP(21) + (t395 * t421 + t396 * t423 + t398 * t403 + t399 * t404 + t418 * t467 + t431 * t448) * MDP(25) + (-t529 * t611 + (t516 * t641 + (qJD(1) * t587 + t557) * t723) * t633) * MDP(13) + (t528 * t611 + (t517 * t641 + (-qJD(1) * t586 - t555) * t723) * t633) * MDP(14) + (-t657 * t611 + t581 * t555 + t569 * t517 + t567 * t586 + t537 * t528 + (t685 * t641 + (qJD(1) * t690 - t475) * t723) * t633) * MDP(16) + (-t436 * t587 + t440 * t611 - t446 * t557 - t464 * t529 + t486 * t516 + (t435 * t641 + (-qJD(1) * t487 - t470) * t723) * t633) * MDP(20) + (-t436 * t586 - t444 * t611 - t446 * t555 - t464 * t528 - t486 * t517 + (-t437 * t641 + (qJD(1) * t488 + t469) * t723) * t633) * MDP(19) + (-t611 * t633 - t629 * t724) * MDP(15) * t723 + 0.2e1 * t687 * t712; ((t503 * t636 + t515 * t639) * t516 - t387 * t637 + t583 * t407 - t402 * t574 + (t636 * t673 + t639 * t672) * t550 - t729 * t663 + t735 * t432 + t390 * t653) * MDP(32) + (-t407 * t574 - t663 * t735) * MDP(26) + (t407 * t573 + t408 * t574 - t450 * t735 + t663 * t734) * MDP(27) + (t403 * t541 + t404 * t540 - t490 * t521 - t491 * t520 + (t395 * t632 - t396 * t634) * t640 + t738 * t513 - t737 * t512 + t666 * t720) * MDP(24) + (t407 * t637 + t516 * t574 + t550 * t735 + t653 * t663) * MDP(28) + (-t557 * t653 - t514) * MDP(11) + (-t550 * t653 - t514) * MDP(30) + (-(t503 * t639 - t515 * t636) * t516 + t388 * t637 + t583 * t408 - t402 * t573 + (t636 * t672 - t639 * t673) * t550 + t729 * t450 + t734 * t432 - t389 * t653) * MDP(31) + (pkin(1) * t638 * t746 + t579 * t659 - t567) * MDP(9) + (pkin(8) * t678 + t576 * t659 + (-qJD(2) * t693 + t746) * t641 * pkin(1)) * MDP(10) + ((-t516 + t749) * t640 + (-t517 + t748) * t637) * MDP(12) + (-t687 + (-MDP(6) * t641 + MDP(7) * t638) * t633 * t759) * t642 + (t395 * t637 + t490 * t613 - t516 * t520 - t738 * t557 - t736 * t512 + t770 * t448 + (-t403 * t611 + t756) * t640) * MDP(22) + (t436 * t640 + t494 * t611 - t517 * t606 + t733 * t555 + (-t464 * t637 - t710) * qJD(3) + (-t469 * t638 + (t464 * t641 + t709) * t637) * t725) * MDP(19) + (t436 * t606 - t469 * t494 - t470 * t493 - t733 * t464 + (-t435 * t640 + t437 * t637 + (t469 * t640 + t470 * t637) * qJD(3)) * pkin(9)) * MDP(21) + (-t408 * t637 + t450 * t653 - t516 * t573 - t550 * t734) * MDP(29) + t611 * MDP(15) * t702 + (-t396 * t637 + t491 * t613 + t516 * t521 + t737 * t557 + t736 * t513 + (t632 * t720 - t541) * t448 + (t404 * t611 - t757) * t640) * MDP(23) + (-t493 * t555 - t494 * t557 + (-t435 - t611 * t469 + (-t517 + t721) * pkin(9)) * t640 + ((qJD(3) * t555 - t516) * pkin(9) + t766) * t637) * MDP(18) + (t395 * t520 + t396 * t521 - t403 * t738 - t404 * t737 + t418 * t613 + t448 * t736) * MDP(25) + (pkin(2) * t516 + t567 * t637 - t730 * t611 - t579 * t557 + (t537 * t640 - t711) * qJD(3) + (-t537 * t739 + (t476 - t708) * t638) * t725) * MDP(17) + (-t611 * t719 + (t611 * t739 + (qJD(2) * t637 - t557) * t638) * t725) * MDP(13) + (-t436 * t637 - t493 * t611 + t516 * t606 + t733 * t557 + (-t464 * t640 + t711) * qJD(3) + (t464 * t739 + (t470 + t708) * t638) * t725) * MDP(20) + (t611 * t720 + (-t611 * t740 + (t555 + t722) * t638) * t725) * MDP(14) + t746 * t768 + (-pkin(2) * t517 - t567 * t640 + t689 * t611 - t579 * t555 + (t537 * t637 + t710) * qJD(3) + (t475 * t638 + (-t537 * t641 - t709) * t637) * t725) * MDP(16); MDP(11) * t750 + (-t555 ^ 2 + t554) * MDP(12) + t645 * MDP(13) + (-t748 - t517) * MDP(14) + MDP(15) * t678 + (-t537 * t557 + t650) * MDP(16) + (t475 * t611 + t537 * t555 + t686) * MDP(17) + (pkin(3) * t516 - qJ(4) * t517 + (-t470 - t476) * t557 + (t469 + t716) * t555) * MDP(18) + (t499 * t555 - t650 - 0.2e1 * t671 + t755) * MDP(19) + (-t464 * t555 + t499 * t557 + t611 * t716 + 0.2e1 * t607 + t670) * MDP(20) + (-pkin(3) * t437 - qJ(4) * t435 - t464 * t499 - t469 * t476 + t470 * t716) * MDP(21) + (t634 * t751 + qJ(4) * t490 + t403 * t555 + t757 - t714 * t512 + (-t634 * t713 - t422) * t557) * MDP(22) + (-t632 * t751 + qJ(4) * t491 - t404 * t555 + t756 + t714 * t513 + (t632 * t713 + t424) * t557) * MDP(23) + (t422 * t513 - t424 * t512 + (qJD(5) * t513 - t404 * t557 + t491 * t761 - t395) * t634 + (-qJD(5) * t512 + t403 * t557 + t490 * t761 - t396) * t632) * MDP(24) + (qJ(4) * t418 - t403 * t422 - t404 * t424 - t667 * t761 + t714 * t448 + (-t403 * t634 - t404 * t632) * qJD(5)) * MDP(25) + (-t407 * t661 - t663 * t728) * MDP(26) + (-t407 * t591 + t408 * t661 - t450 * t728 + t663 * t727) * MDP(27) + (-t555 * t663 + t764) * MDP(28) + (-t450 * t555 + t774) * MDP(29) + t550 * t555 * MDP(30) + (-(-t600 * t636 + t601 * t639) * t516 + t624 * t408 + t402 * t591 + t389 * t555 + (t636 * t651 - t639 * t652) * t550 + t715 * t450 + t727 * t432) * MDP(31) + ((t600 * t639 + t601 * t636) * t516 + t624 * t407 - t402 * t661 - t390 * t555 + (t636 * t652 + t639 * t651) * t550 - t715 * t663 + t728 * t432) * MDP(32); t645 * MDP(18) + (t678 - t750) * MDP(19) + (-t611 ^ 2 - t554) * MDP(20) + (t755 + t766) * MDP(21) + (-t512 * t611 - t516 * t634 - t554 * t632) * MDP(22) + (t513 * t611 + t516 * t632 - t554 * t634) * MDP(23) + (-t490 * t632 - t491 * t634 + (t512 * t634 + t513 * t632) * t557) * MDP(24) + (t448 * t611 + t557 * t666 + t667) * MDP(25) + (t611 * t450 + t764) * MDP(31) + (-t611 * t663 + t774) * MDP(32); (t513 * t557 + t490) * MDP(22) + (t512 * t557 + t491) * MDP(23) + (-t512 ^ 2 - t513 ^ 2) * MDP(24) + (t403 * t513 - t404 * t512 + t418) * MDP(25) + (t408 - t772) * MDP(31) + (t407 - t777) * MDP(32); -t663 * t450 * MDP(26) + (-t450 ^ 2 + t663 ^ 2) * MDP(27) + (t706 + t777) * MDP(28) + (-t691 - t772) * MDP(29) - t516 * MDP(30) + (t390 * t550 + t432 * t663 + t692) * MDP(31) + (t389 * t550 + t432 * t450 - t668) * MDP(32) + (-MDP(28) * t752 + MDP(29) * t663 - MDP(31) * t390 - MDP(32) * t389) * qJD(6);];
tauc  = t1;
