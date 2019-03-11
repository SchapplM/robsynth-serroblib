% Calculate vector of inverse dynamics joint torques for
% S6RRPRPP3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RRPRPP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:57:26
% EndTime: 2019-03-09 09:57:39
% DurationCPUTime: 10.63s
% Computational Cost: add. (7106->645), mult. (15650->755), div. (0->0), fcn. (11050->10), ass. (0->262)
t777 = MDP(23) - MDP(28);
t606 = sin(qJ(4));
t609 = cos(qJ(2));
t603 = cos(pkin(9));
t758 = cos(qJ(4));
t693 = t758 * t603;
t673 = t609 * t693;
t602 = sin(pkin(9));
t708 = qJD(1) * t609;
t692 = t602 * t708;
t687 = qJD(4) * t758;
t703 = qJD(4) * t606;
t765 = -t602 * t703 + t603 * t687;
t716 = -qJD(1) * t673 + t606 * t692 + t765;
t607 = sin(qJ(2));
t709 = qJD(1) * t607;
t688 = t602 * t709;
t706 = qJD(2) * t603;
t534 = -t688 + t706;
t691 = t603 * t709;
t707 = qJD(2) * t602;
t535 = t691 + t707;
t476 = -t758 * t534 + t535 * t606;
t760 = t476 ^ 2;
t646 = -t606 * t534 - t535 * t758;
t474 = t646 ^ 2;
t581 = -qJD(4) + t708;
t775 = t476 * t581;
t739 = t646 * t581;
t592 = t609 * qJDD(1);
t699 = qJD(1) * qJD(2);
t641 = t607 * t699 - t592;
t537 = qJDD(4) + t641;
t746 = pkin(4) + qJ(6);
t685 = t746 * t537;
t660 = pkin(2) * t607 - qJ(3) * t609;
t541 = t660 * qJD(1);
t497 = pkin(7) * t688 + t603 * t541;
t730 = t603 * t609;
t652 = pkin(3) * t607 - pkin(8) * t730;
t464 = qJD(1) * t652 + t497;
t525 = t602 * t541;
t731 = t603 * t607;
t732 = t602 * t609;
t645 = -pkin(7) * t731 - pkin(8) * t732;
t481 = qJD(1) * t645 + t525;
t744 = pkin(8) + qJ(3);
t553 = t744 * t602;
t554 = t744 * t603;
t774 = qJD(3) * t693 - t758 * t481 - t553 * t687 + (-qJD(3) * t602 - qJD(4) * t554 - t464) * t606;
t539 = t602 * t758 + t606 * t603;
t524 = t539 * qJD(4);
t633 = t609 * t539;
t715 = -qJD(1) * t633 + t524;
t682 = t609 * t699;
t698 = qJDD(1) * t607;
t642 = t682 + t698;
t610 = cos(qJ(1));
t728 = t607 * t610;
t608 = sin(qJ(1));
t729 = t607 * t608;
t773 = g(1) * t728 + g(2) * t729;
t661 = pkin(2) * t609 + qJ(3) * t607;
t549 = -pkin(1) - t661;
t527 = t549 * qJD(1);
t588 = pkin(7) * t708;
t555 = qJD(2) * qJ(3) + t588;
t484 = t603 * t527 - t555 * t602;
t442 = -pkin(3) * t708 - pkin(8) * t535 + t484;
t485 = t602 * t527 + t603 * t555;
t449 = pkin(8) * t534 + t485;
t411 = -t758 * t442 + t449 * t606;
t649 = -pkin(5) * t646 + t411;
t701 = qJD(5) + t649;
t772 = MDP(24) + MDP(27);
t749 = g(2) * t608;
t667 = g(1) * t610 + t749;
t771 = t609 * t667;
t648 = t667 * t607;
t697 = qJDD(2) * t602;
t622 = t642 * t603 + t697;
t713 = t642 * t602;
t659 = qJDD(2) * t603 - t713;
t616 = qJD(4) * t646 - t606 * t622 + t758 * t659;
t770 = -qJ(6) * t616 + t476 * qJD(6);
t720 = qJ(5) * t709 - t774;
t492 = -t606 * t553 + t554 * t758;
t769 = qJD(3) * t539 + qJD(4) * t492 + t464 * t758 - t606 * t481;
t528 = t537 * qJ(5);
t563 = qJD(5) * t581;
t767 = t563 - t528;
t529 = pkin(3) * t692 + t588;
t766 = -qJ(5) * t716 - qJD(5) * t539 - t529;
t764 = MDP(22) + MDP(26);
t763 = pkin(7) * t682 + qJDD(3);
t762 = pkin(5) * t616 + qJDD(6);
t761 = -0.2e1 * pkin(1);
t575 = t581 ^ 2;
t759 = 0.2e1 * t528;
t757 = pkin(3) * t602;
t756 = pkin(4) * t537;
t422 = -t534 * t687 + t535 * t703 - t606 * t659 - t758 * t622;
t755 = pkin(5) * t422;
t753 = pkin(5) * t476;
t752 = pkin(7) * t534;
t751 = g(1) * t608;
t748 = g(2) * t610;
t747 = g(3) * t607;
t598 = g(3) * t609;
t745 = pkin(5) + t744;
t743 = qJ(5) * t616;
t742 = qJ(5) * t476;
t741 = qJDD(2) * pkin(2);
t412 = t606 * t442 + t758 * t449;
t740 = t412 * t581;
t738 = t476 * t646;
t599 = pkin(9) + qJ(4);
t590 = sin(t599);
t737 = t590 * t607;
t736 = t590 * t609;
t591 = cos(t599);
t735 = t591 * t607;
t734 = t591 * t609;
t733 = t602 * t607;
t727 = t608 * t609;
t585 = pkin(3) * t603 + pkin(2);
t558 = t609 * t585;
t726 = t609 * t610;
t725 = t610 * t590;
t538 = t602 * t606 - t693;
t724 = qJD(6) * t538 + t715 * t746 + t766;
t723 = -t715 * pkin(5) - t720;
t686 = t607 * t746;
t722 = t716 * pkin(5) + qJD(1) * t686 + t769;
t721 = t715 * pkin(4) + t766;
t719 = pkin(4) * t709 + t769;
t533 = t603 * t549;
t483 = -pkin(8) * t731 + t533 + (-pkin(7) * t602 - pkin(3)) * t609;
t501 = pkin(7) * t730 + t602 * t549;
t490 = -pkin(8) * t733 + t501;
t717 = t606 * t483 + t758 * t490;
t521 = qJD(2) * t660 - qJD(3) * t607;
t472 = qJD(1) * t521 + qJDD(1) * t549;
t510 = -pkin(7) * t641 + qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t434 = t602 * t472 + t603 * t510;
t705 = qJD(2) * t607;
t694 = pkin(7) * t705;
t488 = t603 * t521 + t602 * t694;
t586 = pkin(7) * t698;
t712 = -t586 - t598;
t704 = qJD(2) * t609;
t690 = t602 * t704;
t530 = pkin(3) * t690 + pkin(7) * t704;
t542 = pkin(3) * t733 + t607 * pkin(7);
t711 = t610 * pkin(1) + t608 * pkin(7);
t600 = t607 ^ 2;
t710 = -t609 ^ 2 + t600;
t702 = -qJD(5) - t411;
t402 = t412 - t753;
t700 = -qJD(6) - t402;
t684 = t745 * t610;
t681 = -pkin(1) - t558;
t519 = t586 - t741 + t763;
t680 = -t519 - t598;
t679 = -qJD(2) * pkin(2) + qJD(3);
t511 = t590 * t727 + t591 * t610;
t512 = t591 * t727 - t725;
t678 = -t511 * pkin(4) + qJ(5) * t512;
t513 = -t608 * t591 + t609 * t725;
t514 = t590 * t608 + t591 * t726;
t677 = -t513 * pkin(4) + qJ(5) * t514;
t676 = -qJ(5) * t590 - t585;
t433 = t603 * t472 - t602 * t510;
t415 = pkin(3) * t641 - pkin(8) * t622 + t433;
t424 = pkin(8) * t659 + t434;
t675 = -t606 * t415 - t758 * t424 - t442 * t687 + t449 * t703;
t674 = -t758 * t415 + t606 * t424 + t442 * t703 + t449 * t687;
t672 = g(3) * (pkin(4) * t734 + qJ(5) * t736 + t558);
t583 = g(1) * t729;
t670 = -g(2) * t728 + t583;
t669 = g(1) * t511 - g(2) * t513;
t668 = g(1) * t512 - g(2) * t514;
t666 = -t748 + t751;
t548 = pkin(7) * t709 + t679;
t407 = qJ(5) * t581 - t412;
t429 = qJ(5) * t609 - t717;
t518 = -t606 * t733 + t607 * t693;
t664 = -qJ(5) * t518 + t542;
t662 = t483 * t758 - t606 * t490;
t657 = -qJDD(5) - t674;
t655 = -g(3) * t736 + t773 * t590;
t654 = -g(3) * t734 + t773 * t591;
t651 = -qJ(5) * t539 - t585;
t430 = t609 * pkin(4) - t662;
t391 = t675 + t767;
t647 = -pkin(7) * qJDD(2) + t699 * t761;
t491 = t553 * t758 + t606 * t554;
t596 = t610 * pkin(7);
t644 = -t512 * pkin(4) - qJ(5) * t511 + t610 * t757 - t729 * t744 + t596;
t453 = qJD(2) * t652 + t488;
t508 = t602 * t521;
t465 = qJD(2) * t645 + t508;
t643 = -t453 * t758 + t606 * t465 + t483 * t703 + t490 * t687;
t640 = t603 * t698 + t697;
t639 = t606 * t453 + t758 * t465 + t483 * t687 - t490 * t703;
t612 = qJD(1) ^ 2;
t638 = pkin(1) * t612 + t667;
t496 = -pkin(3) * t534 + t548;
t611 = qJD(2) ^ 2;
t637 = pkin(7) * t611 + qJDD(1) * t761 + t748;
t458 = -qJD(2) * t673 + t524 * t607 + t606 * t690;
t636 = qJ(5) * t458 - qJD(5) * t518 + t530;
t634 = t514 * pkin(4) + qJ(5) * t513 + t585 * t726 + t608 * t757 + t711;
t632 = -t491 * t537 + t654;
t631 = t492 * t537 + t655;
t629 = -t747 - t771;
t628 = -t648 - t741;
t404 = -t422 - t775;
t626 = qJ(5) * t646 + t496;
t625 = g(1) * t513 + g(2) * t511 + g(3) * t737 - t674;
t624 = t680 + t648;
t623 = -qJDD(5) + t625;
t621 = g(1) * t514 + g(2) * t512 + g(3) * t735 + t675;
t454 = -pkin(3) * t659 + t519;
t397 = -qJ(5) * t705 + qJD(5) * t609 - t639;
t417 = pkin(4) * t476 + t626;
t620 = -t417 * t646 - t623;
t403 = t476 * t746 + t626;
t618 = -t403 * t646 - t623 - t755;
t617 = -t403 * t476 - t621 + t762;
t393 = -pkin(4) * t616 + t422 * qJ(5) + qJD(5) * t646 + t454;
t615 = t616 + t739;
t613 = t598 - t648 + t393;
t556 = qJ(5) * t735;
t517 = t539 * t607;
t500 = -pkin(7) * t732 + t533;
t498 = -t691 * pkin(7) + t525;
t489 = -t603 * t694 + t508;
t473 = pkin(4) * t538 + t651;
t459 = qJD(2) * t633 + t607 * t765;
t457 = -t538 * pkin(5) + t492;
t456 = t539 * pkin(5) + t491;
t447 = t538 * t746 + t651;
t443 = pkin(4) * t517 + t664;
t432 = -pkin(4) * t646 + t742;
t431 = t517 * t746 + t664;
t418 = -pkin(5) * t517 - t429;
t416 = t518 * pkin(5) + t609 * qJ(6) + t430;
t410 = -t646 * t746 + t742;
t406 = pkin(4) * t581 - t702;
t405 = pkin(4) * t459 + t636;
t400 = qJD(6) - t407 - t753;
t399 = t581 * t746 + t701;
t398 = -pkin(4) * t705 + t643;
t396 = qJD(6) * t517 + t459 * t746 + t636;
t395 = -pkin(5) * t459 - t397;
t394 = -t458 * pkin(5) - qJD(2) * t686 + t609 * qJD(6) + t643;
t392 = -t657 - t756;
t390 = t393 + t770;
t389 = -t391 + t762;
t388 = qJD(6) * t581 - t657 - t685 - t755;
t1 = [(t391 * t517 + t392 * t518 + t397 * t476 - t398 * t646 - t406 * t458 + t407 * t459 - t422 * t430 - t429 * t616 + t670) * MDP(22) + (t388 * t518 - t389 * t517 - t394 * t646 - t395 * t476 - t399 * t458 - t400 * t459 - t416 * t422 + t418 * t616 + t670) * MDP(26) + (t422 * t517 + t458 * t476 + t459 * t646 + t518 * t616) * MDP(16) + (-t411 * t705 + t454 * t517 + t496 * t459 + t530 * t476 + t537 * t662 - t542 * t616 + t581 * t643 + t609 * t674 + t668) * MDP(20) + (t459 * t581 - t476 * t705 - t517 * t537 - t609 * t616) * MDP(18) + (t388 * t609 + t390 * t517 + t394 * t581 + t396 * t476 - t399 * t705 + t403 * t459 - t416 * t537 - t431 * t616 + t668) * MDP(28) + (-t392 * t609 - t393 * t517 - t398 * t581 - t405 * t476 + t406 * t705 - t417 * t459 + t430 * t537 + t443 * t616 - t668) * MDP(23) + (qJDD(1) * t600 + 0.2e1 * t607 * t682) * MDP(4) + (t607 * t637 + t609 * t647 - t583) * MDP(10) + t666 * MDP(2) + t667 * MDP(3) + (-t422 * t518 + t458 * t646) * MDP(15) + (t391 * t609 - t393 * t518 + t397 * t581 + t405 * t646 - t407 * t705 + t417 * t458 + t422 * t443 - t429 * t537 + t669) * MDP(24) + (t422 * t609 + t458 * t581 + t518 * t537 - t646 * t705) * MDP(17) + (-t389 * t609 - t390 * t518 - t395 * t581 + t396 * t646 + t400 * t705 + t403 * t458 + t418 * t537 + t422 * t431 + t669) * MDP(27) + (-t412 * t705 - t542 * t422 + t454 * t518 - t496 * t458 - t530 * t646 - t537 * t717 + t581 * t639 - t609 * t675 - t669) * MDP(21) + (t393 * t443 + t417 * t405 + t391 * t429 + t407 * t397 + t392 * t430 + t406 * t398 - g(1) * (t608 * t681 + t644) - g(2) * (t728 * t744 + t634)) * MDP(25) + qJDD(1) * MDP(1) + (-t667 * t603 + (t519 * t603 + (-qJD(1) * t501 - t485) * qJD(2) + t640 * pkin(7)) * t607 + (t489 * qJD(1) + t501 * qJDD(1) + t434 - t666 * t602 + (t548 * t603 + (t535 + t691) * pkin(7)) * qJD(2)) * t609) * MDP(12) + (qJDD(2) * t607 + t609 * t611) * MDP(6) + (qJDD(2) * t609 - t607 * t611) * MDP(7) + (-t537 * t609 - t581 * t705) * MDP(19) + 0.2e1 * (t592 * t607 - t699 * t710) * MDP(5) + (t489 * t534 - t501 * t713 - t488 * t535 + (-qJDD(2) * t500 - t434 * t607 - t485 * t704) * t602 + (t501 * qJDD(2) - t433 * t607 - t484 * t704 - t500 * t642) * t603 + t670) * MDP(13) + (t390 * t431 + t403 * t396 + t388 * t416 + t399 * t394 + t389 * t418 + t400 * t395 - g(1) * (-qJ(6) * t512 + t644) - g(2) * (qJ(6) * t514 + t607 * t684 + t634) - (-pkin(5) * t607 + t681) * t751) * MDP(29) + (t434 * t501 + t485 * t489 + t433 * t500 + t484 * t488 - g(1) * t596 - g(2) * (t610 * t661 + t711) - t549 * t751 + (t519 * t607 + t548 * t704) * pkin(7)) * MDP(14) + (t647 * t607 + (-t637 + t751) * t609) * MDP(9) + (-t667 * t602 + (-pkin(7) * t659 + t519 * t602 + (qJD(1) * t500 + t484) * qJD(2)) * t607 + (-t488 * qJD(1) - t500 * qJDD(1) - t433 + t666 * t603 + (t548 * t602 - t752) * qJD(2)) * t609) * MDP(11); (t412 * t709 + t585 * t422 + t454 * t539 + t716 * t496 + t529 * t646 + t774 * t581 - t631) * MDP(21) + (t411 * t709 + t454 * t538 - t529 * t476 + t715 * t496 + t581 * t769 + t585 * t616 + t632) * MDP(20) + (t422 * t538 - t476 * t716 + t539 * t616 + t646 * t715) * MDP(16) + (t391 * t538 + t392 * t539 + t406 * t716 + t407 * t715 - t422 * t491 + t476 * t720 + t492 * t616 - t646 * t719 + t629) * MDP(22) + (t388 * t539 - t389 * t538 + t399 * t716 - t400 * t715 - t422 * t456 + t457 * t616 - t476 * t723 - t646 * t722 + t629) * MDP(26) + (-t393 * t538 - t406 * t709 - t417 * t715 + t473 * t616 - t476 * t721 - t581 * t719 - t632) * MDP(23) + (t390 * t538 + t399 * t709 + t403 * t715 - t447 * t616 - t456 * t537 + t476 * t724 + t581 * t722 + t654) * MDP(28) + MDP(7) * t592 + MDP(6) * t698 + (t390 * t447 + t388 * t456 + t389 * t457 - t672 + t724 * t403 + t723 * t400 + t722 * t399 + (-g(3) * qJ(6) * t591 - g(1) * t684 - t745 * t749) * t609 + (-g(3) * t745 + t667 * (t591 * t746 - t676)) * t607) * MDP(29) + (-t422 * t539 - t646 * t716) * MDP(15) + (t537 * t539 - t581 * t716 + t646 * t709) * MDP(17) + (-t393 * t539 + t407 * t709 - t417 * t716 + t422 * t473 + t581 * t720 + t646 * t721 + t631) * MDP(24) + (-t390 * t539 - t400 * t709 - t403 * t716 + t422 * t447 + t457 * t537 - t581 * t723 + t646 * t724 + t655) * MDP(27) + (t393 * t473 - t391 * t492 + t392 * t491 - t672 - t744 * t771 + t721 * t417 + t720 * t407 + t719 * t406 + (-g(3) * t744 + t667 * (pkin(4) * t591 - t676)) * t607) * MDP(25) + (-t548 * t588 - t484 * t497 - t485 * t498 + (-t484 * t602 + t485 * t603) * qJD(3) + t624 * pkin(2) + (-t433 * t602 + t434 * t603 + t629) * qJ(3)) * MDP(14) + qJDD(2) * MDP(8) + (-t660 * t603 * qJDD(1) + (t628 - t680) * t602 + ((-qJ(3) * t706 + t485) * t607 + (-pkin(7) * t535 - t498 + (-t548 + t679) * t603) * t609) * qJD(1)) * MDP(12) + (t497 * t535 - t498 * t534 + (qJ(3) * t659 + qJD(3) * t534 + t484 * t708 + t434) * t603 + (qJ(3) * t622 + qJD(3) * t535 + t485 * t708 - t433) * t602 + t629) * MDP(13) + (t607 * t638 + t712) * MDP(9) + (t476 * t709 - t537 * t538 + t581 * t715) * MDP(18) + (-MDP(4) * t607 * t609 + MDP(5) * t710) * t612 + t581 * MDP(19) * t709 + (t747 + (-pkin(7) * qJDD(1) + t638) * t609) * MDP(10) + (t602 * qJ(3) * t592 - pkin(2) * t713 + (t624 + t741) * t603 + ((-qJ(3) * t707 - t484) * t607 + (t752 + t497 + (qJD(3) - t548) * t602) * t609) * qJD(1)) * MDP(11); (-t535 * t708 - t659) * MDP(11) + ((-t534 + t706) * t708 + t640) * MDP(12) + (-t534 ^ 2 - t535 ^ 2) * MDP(13) + (t484 * t535 - t485 * t534 + t628 - t712 + t763) * MDP(14) + (t406 * t646 - t407 * t476 + t613) * MDP(25) + (t399 * t646 + t400 * t476 + t613 + t770) * MDP(29) + t764 * (-t474 - t760) + (-MDP(21) + t772) * (t422 - t775) + (-MDP(20) + t777) * (t616 - t739); -MDP(15) * t738 + (t474 - t760) * MDP(16) + t404 * MDP(17) + t615 * MDP(18) + t537 * MDP(19) + (t496 * t646 + t625 - t740) * MDP(20) + (t411 * t581 + t476 * t496 + t621) * MDP(21) + (pkin(4) * t422 + t743 - (-t407 - t412) * t646 + (t406 + t702) * t476) * MDP(22) + (t432 * t476 + t620 + t740 - 0.2e1 * t756) * MDP(23) + (-t417 * t476 - t432 * t646 + t581 * t702 - t563 - t621 + t759) * MDP(24) + (-t391 * qJ(5) - t392 * pkin(4) - t417 * t432 - t406 * t412 - g(1) * t677 - g(2) * t678 - g(3) * (-pkin(4) * t737 + t556) + t702 * t407) * MDP(25) + (t743 + t422 * t746 - (t400 + t700) * t646 + (t399 - t701) * t476) * MDP(26) + (-t410 * t646 - t581 * t649 - 0.2e1 * t563 + t617 + t759) * MDP(27) + (-t410 * t476 + (-0.2e1 * qJD(6) - t402) * t581 + 0.2e1 * t685 - t618) * MDP(28) + (-t388 * t746 + t389 * qJ(5) - t403 * t410 - g(1) * (-qJ(6) * t513 + t677) - g(2) * (-qJ(6) * t511 + t678) - g(3) * (-t590 * t686 + t556) + t701 * t400 + t700 * t399) * MDP(29); (-t407 * t581 + t620 - t756) * MDP(25) + ((qJD(6) + t400) * t581 - t685 + t618) * MDP(29) + t777 * (t537 + t738) + t764 * t404 + t772 * (-t474 - t575); t615 * MDP(26) + (t537 - t738) * MDP(27) + (-t575 - t760) * MDP(28) + (-t399 * t581 + t617 - t767) * MDP(29);];
tau  = t1;
