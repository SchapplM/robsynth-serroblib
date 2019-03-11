% Calculate vector of inverse dynamics joint torques for
% S6RRRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRPRR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:13:43
% EndTime: 2019-03-09 18:13:55
% DurationCPUTime: 9.72s
% Computational Cost: add. (7478->551), mult. (16448->679), div. (0->0), fcn. (12254->12), ass. (0->244)
t628 = qJD(2) + qJD(3);
t636 = sin(qJ(3));
t637 = sin(qJ(2));
t746 = t636 * t637;
t685 = t628 * t746;
t641 = cos(qJ(2));
t772 = cos(qJ(3));
t713 = t772 * t641;
t691 = qJD(1) * t713;
t707 = qJDD(1) * t772;
t719 = qJDD(1) * t641;
t692 = -t628 * t691 - t636 * t719 - t637 * t707;
t473 = qJD(1) * t685 + t692;
t565 = t636 * t641 + t637 * t772;
t512 = t628 * t565;
t720 = qJDD(1) * t637;
t683 = t636 * t720 - t641 * t707;
t474 = qJD(1) * t512 + t683;
t733 = qJD(1) * t637;
t712 = t636 * t733;
t549 = -t691 + t712;
t551 = t565 * qJD(1);
t635 = sin(qJ(5));
t640 = cos(qJ(5));
t730 = qJD(5) * t640;
t731 = qJD(5) * t635;
t427 = -t640 * t473 + t635 * t474 + t549 * t730 - t551 * t731;
t627 = qJDD(2) + qJDD(3);
t619 = -qJDD(5) + t627;
t634 = sin(qJ(6));
t639 = cos(qJ(6));
t723 = qJD(5) - t628;
t726 = qJD(6) * t639;
t715 = t639 * t427 - t634 * t619 + t723 * t726;
t727 = qJD(6) * t634;
t780 = t549 * t635 + t640 * t551;
t414 = -t727 * t780 + t715;
t477 = t634 * t723 + t639 * t780;
t592 = t639 * t619;
t415 = qJD(6) * t477 + t427 * t634 + t592;
t428 = qJD(5) * t780 - t473 * t635 - t640 * t474;
t475 = t634 * t780 - t639 * t723;
t426 = qJDD(6) + t428;
t745 = t639 * t426;
t748 = t634 * t426;
t782 = -t640 * t549 + t551 * t635;
t791 = qJD(6) + t782;
t759 = t477 * t791;
t760 = t475 * t791;
t762 = t414 * t634;
t800 = t791 * t639;
t808 = -((t415 + t759) * t634 - (t414 - t760) * t639) * MDP(30) + (t477 * t800 + t762) * MDP(29) + (-t477 * t780 + t791 * t800 + t748) * MDP(31) - (t634 * t791 ^ 2 - t475 * t780 - t745) * MDP(32) + (t780 ^ 2 - t782 ^ 2) * MDP(23) + (t723 * t782 + t427) * MDP(24) + (t723 * t780 - t428) * MDP(25) - t619 * MDP(26) + (MDP(22) * t782 - MDP(33) * t791) * t780;
t541 = t551 * pkin(9);
t773 = pkin(8) + pkin(7);
t580 = t773 * t637;
t570 = qJD(1) * t580;
t763 = qJD(2) * pkin(2);
t559 = -t570 + t763;
t581 = t773 * t641;
t572 = qJD(1) * t581;
t747 = t636 * t572;
t498 = t772 * t559 - t747;
t779 = qJD(4) - t498;
t792 = -t541 + t779;
t558 = t772 * t572;
t509 = -t636 * t570 + t558;
t732 = qJD(3) * t636;
t690 = pkin(2) * t732 - t509;
t510 = -t772 * t570 - t747;
t710 = qJD(3) * t772;
t738 = pkin(2) * t710 + qJD(4) - t510;
t797 = pkin(5) * t780;
t626 = t641 * pkin(2);
t613 = t626 + pkin(1);
t579 = t613 * qJD(1);
t482 = t549 * pkin(3) - t551 * qJ(4) - t579;
t460 = -pkin(4) * t549 - t482;
t429 = pkin(5) * t782 - pkin(10) * t780 + t460;
t644 = -pkin(3) - pkin(4);
t458 = t628 * t644 + t792;
t499 = t636 * t559 + t558;
t764 = t549 * pkin(9);
t472 = t499 + t764;
t621 = t628 * qJ(4);
t464 = t472 + t621;
t436 = t458 * t635 + t464 * t640;
t433 = pkin(10) * t723 + t436;
t411 = t429 * t634 + t433 * t639;
t796 = t411 * t780;
t795 = t764 - t690;
t794 = -t541 + t738;
t721 = qJD(1) * qJD(2);
t793 = t637 * t721 - t719;
t638 = sin(qJ(1));
t631 = qJ(2) + qJ(3);
t624 = cos(t631);
t750 = t624 * t638;
t623 = sin(t631);
t752 = t623 * t640;
t521 = t635 * t750 - t638 * t752;
t642 = cos(qJ(1));
t749 = t624 * t642;
t751 = t623 * t642;
t523 = t635 * t749 - t640 * t751;
t675 = t623 * t635 + t624 * t640;
t666 = g(1) * t523 + g(2) * t521 + g(3) * t675;
t708 = t641 * t721;
t515 = qJDD(2) * pkin(2) + t773 * (-t708 - t720);
t516 = t773 * t793;
t694 = -t772 * t515 - t636 * t516 + t559 * t732 + t572 * t710;
t620 = t627 * pkin(3);
t778 = qJDD(4) - t620;
t442 = t694 + t778;
t422 = -pkin(4) * t627 + pkin(9) * t473 + t442;
t693 = -t636 * t515 + t772 * t516 - t559 * t710 + t572 * t732;
t616 = t627 * qJ(4);
t618 = t628 * qJD(4);
t781 = t616 + t618;
t439 = -t693 + t781;
t423 = pkin(9) * t474 + t439;
t673 = -t640 * t422 + t635 * t423 + t458 * t731 + t464 * t730;
t405 = pkin(5) * t619 + t673;
t410 = t429 * t639 - t433 * t634;
t705 = t405 * t639 + t410 * t780;
t524 = t675 * t642;
t667 = t635 * t422 + t640 * t423 + t458 * t730 - t464 * t731;
t522 = t675 * t638;
t686 = g(2) * t522 + g(3) * (-t624 * t635 + t752);
t648 = -g(1) * t524 - t460 * t782 + t667 - t686;
t649 = t460 * t780 - t666 + t673;
t496 = t551 * pkin(3) + t549 * qJ(4);
t771 = pkin(4) * t551;
t470 = -t496 - t771;
t434 = -pkin(10) * t782 + t470 - t797;
t736 = t640 * qJ(4) + t635 * t644;
t569 = -pkin(10) + t736;
t784 = t791 * (qJD(6) * t569 + t434);
t612 = -pkin(2) * t772 - pkin(3);
t605 = -pkin(4) + t612;
t607 = pkin(2) * t636 + qJ(4);
t739 = t635 * t605 + t640 * t607;
t526 = -pkin(10) + t739;
t718 = pkin(2) * t733;
t783 = t791 * (qJD(6) * t526 + t434 - t718);
t520 = -t636 * t580 + t772 * t581;
t737 = t624 * pkin(3) + t623 * qJ(4);
t777 = -t791 * (pkin(10) * t791 + t797) + t666;
t714 = qJD(2) * t773;
t571 = t637 * t714;
t573 = t641 * t714;
t461 = -t772 * t571 - t636 * t573 - t580 * t710 - t581 * t732;
t766 = g(2) * t642;
t687 = g(1) * t638 - t766;
t776 = t461 * t628 + t520 * t627 + t623 * t687;
t447 = pkin(9) * t512 + t461;
t462 = qJD(3) * t520 - t636 * t571 + t772 * t573;
t511 = -qJD(2) * t713 - t641 * t710 + t685;
t448 = t511 * pkin(9) + t462;
t519 = t580 * t772 + t636 * t581;
t487 = -t565 * pkin(9) + t519;
t564 = -t713 + t746;
t488 = pkin(9) * t564 + t520;
t679 = t487 * t640 - t488 * t635;
t412 = qJD(5) * t679 + t447 * t640 + t448 * t635;
t435 = t458 * t640 - t464 * t635;
t432 = -pkin(5) * t723 - t435;
t497 = t564 * pkin(3) - t565 * qJ(4) - t613;
t479 = -pkin(4) * t564 - t497;
t504 = t564 * t635 + t565 * t640;
t677 = t640 * t564 - t565 * t635;
t437 = -pkin(5) * t677 - pkin(10) * t504 + t479;
t445 = qJD(5) * t677 - t511 * t640 + t512 * t635;
t450 = t487 * t635 + t488 * t640;
t698 = -pkin(10) * t619 + qJD(6) * t429 + t667;
t768 = g(1) * t642;
t775 = t405 * t504 - t450 * t426 + t432 * t445 - (qJD(6) * t437 + t412) * t791 + t677 * t698 + t768;
t774 = t551 ^ 2;
t770 = g(1) * t522;
t632 = qJDD(1) * pkin(1);
t761 = t432 * t504;
t758 = t498 * t628;
t757 = t499 * t628;
t756 = t549 * t551;
t753 = t623 * t638;
t676 = t605 * t640 - t607 * t635;
t743 = qJD(5) * t676 - t635 * t795 + t640 * t794;
t742 = t739 * qJD(5) + t635 * t794 + t640 * t795;
t681 = -qJ(4) * t635 + t640 * t644;
t741 = qJD(5) * t681 - t472 * t635 + t640 * t792;
t740 = t736 * qJD(5) + t472 * t640 + t635 * t792;
t629 = t637 ^ 2;
t735 = -t641 ^ 2 + t629;
t729 = qJD(6) * t433;
t728 = qJD(6) * t780;
t717 = t637 * t763;
t543 = pkin(2) * t793 - t632;
t703 = t723 ^ 2;
t702 = t640 * t723;
t701 = t791 * t432;
t689 = -pkin(2) * t637 - pkin(3) * t623;
t688 = g(2) * t638 + t768;
t684 = t437 * t426 + t770;
t682 = -t729 - t766;
t486 = t718 + t496;
t674 = t634 * t666 - t796;
t672 = t613 + t737;
t430 = t474 * pkin(3) + t473 * qJ(4) - t551 * qJD(4) + t543;
t671 = -0.2e1 * pkin(1) * t721 - pkin(7) * qJDD(2);
t670 = t445 * t639 - t504 * t727;
t451 = t512 * pkin(3) + t511 * qJ(4) - t565 * qJD(4) + t717;
t665 = -g(1) * t749 - g(2) * t750 - g(3) * t623 - t693;
t664 = t686 - t698;
t418 = -pkin(4) * t474 - t430;
t663 = g(1) * t751 + g(2) * t753 - g(3) * t624 - t694;
t438 = -pkin(4) * t512 - t451;
t662 = -pkin(10) * t426 + (t432 + t435) * t791;
t645 = qJD(2) ^ 2;
t661 = -pkin(7) * t645 + 0.2e1 * t632 + t687;
t646 = qJD(1) ^ 2;
t660 = pkin(1) * t646 - pkin(7) * qJDD(1) + t688;
t659 = g(1) * t750 - g(2) * t749 - t462 * t628 - t519 * t627;
t658 = -t482 * t549 + t665;
t657 = -t579 * t549 - t665;
t655 = t579 * t551 + t663;
t654 = -t526 * t426 - t743 * t791 - t701;
t653 = -t569 * t426 - t741 * t791 - t701;
t651 = t482 * t551 - t663 + t778;
t455 = -t692 + (t549 - t712) * t628;
t647 = MDP(11) * t756 + t455 * MDP(13) - t683 * MDP(14) + (-t549 ^ 2 + t774) * MDP(12) + t627 * MDP(15) - t808;
t583 = qJ(4) * t749;
t582 = qJ(4) * t750;
t568 = pkin(5) - t681;
t525 = pkin(5) - t676;
t507 = t524 * t639 - t634 * t638;
t506 = -t524 * t634 - t638 * t639;
t490 = t621 + t499;
t485 = -pkin(3) * t628 + t779;
t463 = -t486 - t771;
t446 = qJD(5) * t504 - t511 * t635 - t640 * t512;
t413 = qJD(5) * t450 + t447 * t635 - t448 * t640;
t409 = pkin(5) * t446 - pkin(10) * t445 + t438;
t402 = pkin(5) * t428 - pkin(10) * t427 + t418;
t401 = t639 * t402;
t1 = [0.2e1 * (t637 * t719 - t721 * t735) * MDP(5) + (t427 * t504 + t445 * t780) * MDP(22) + (t427 * t677 - t428 * t504 - t445 * t782 - t446 * t780) * MDP(23) + (-g(1) * t521 + g(2) * t523 - t412 * t723 + t418 * t504 + t427 * t479 + t438 * t780 + t445 * t460 + t450 * t619) * MDP(28) + (-g(2) * t524 - t413 * t723 - t418 * t677 + t428 * t479 + t438 * t782 + t446 * t460 - t619 * t679 + t770) * MDP(27) + (-t446 * t723 - t619 * t677) * MDP(25) + (t445 * t723 - t504 * t619) * MDP(24) + (t430 * t497 + t439 * t520 + t442 * t519 + t482 * t451 + t490 * t461 + t485 * t462 + (-g(1) * t773 - g(2) * t672) * t642 + (g(1) * t672 - g(2) * t773) * t638) * MDP(21) + ((-t475 * t639 - t477 * t634) * t445 + (-t762 - t415 * t639 + (t475 * t634 - t477 * t639) * qJD(6)) * t504) * MDP(30) + (t473 * t564 - t474 * t565 + t511 * t549 - t512 * t551) * MDP(12) + (-t473 * t565 - t511 * t551) * MDP(11) + t687 * MDP(2) + (-g(2) * t506 - t411 * t446 + t413 * t477 - t679 * t414 + (-(-qJD(6) * t450 + t409) * t791 + (t402 - t729) * t677 - qJD(6) * t761 - t684) * t634 + t775 * t639) * MDP(35) + (-g(2) * t507 - t401 * t677 + t410 * t446 + t413 * t475 - t679 * t415 + (t409 * t791 + (t433 * t677 - t450 * t791 + t761) * qJD(6) + t684) * t639 + t775 * t634) * MDP(34) + (-t414 * t677 + t446 * t477 + t504 * t745 + t670 * t791) * MDP(31) + (-t504 * t748 + t415 * t677 - t446 * t475 + (-t445 * t634 - t504 * t726) * t791) * MDP(32) + (-t426 * t677 + t446 * t791) * MDP(33) + (-t474 * t613 - t512 * t579 + t543 * t564 + t549 * t717 + t659) * MDP(16) + (t473 * t613 + t511 * t579 + t543 * t565 + t551 * t717 - t776) * MDP(17) + (qJDD(1) * t629 + 0.2e1 * t637 * t708) * MDP(4) + (-t430 * t565 - t451 * t551 + t473 * t497 + t482 * t511 + t776) * MDP(20) + (qJDD(2) * t637 + t641 * t645) * MDP(6) + (qJDD(2) * t641 - t637 * t645) * MDP(7) + (t430 * t564 + t451 * t549 + t474 * t497 + t482 * t512 + t659) * MDP(18) + (-t512 * t628 - t564 * t627) * MDP(14) + (-t511 * t628 + t565 * t627) * MDP(13) + t688 * MDP(3) + (t637 * t671 + t641 * t661) * MDP(9) + (-t637 * t661 + t641 * t671) * MDP(10) + qJDD(1) * MDP(1) + (t414 * t504 * t639 + t477 * t670) * MDP(29) + (-t439 * t564 + t442 * t565 - t461 * t549 + t462 * t551 - t473 * t519 - t474 * t520 - t485 * t511 - t490 * t512 - t688) * MDP(19); MDP(6) * t720 + (t525 * t415 + t742 * t475 + (-t666 - t783) * t639 + t654 * t634 + t705) * MDP(34) + t647 + MDP(7) * t719 + (-t463 * t782 - t619 * t676 - t723 * t742 + t649) * MDP(27) + (-t463 * t780 + t619 * t739 - t723 * t743 + t648) * MDP(28) + qJDD(2) * MDP(8) + (-t473 * t612 - t474 * t607 + (t490 + t690) * t551 + (t485 - t738) * t549) * MDP(19) + (t510 * t628 + (-t551 * t733 - t627 * t636 - t628 * t710) * pkin(2) + t657) * MDP(17) + (t509 * t628 + (-t549 * t733 + t627 * t772 - t628 * t732) * pkin(2) + t655) * MDP(16) + (t486 * t551 + t607 * t627 + t628 * t738 + t658 + t781) * MDP(20) + (-t486 * t549 - t612 * t627 - t628 * t690 - t651) * MDP(18) + (t439 * t607 + t442 * t612 - t482 * t486 - g(1) * (t642 * t689 + t583) - g(2) * (t638 * t689 + t582) - g(3) * (t626 + t737) + t738 * t490 + t690 * t485) * MDP(21) + (-g(3) * t641 + t637 * t660) * MDP(9) + (g(3) * t637 + t641 * t660) * MDP(10) + (t525 * t414 + t742 * t477 + (-t405 + t783) * t634 + t654 * t639 + t674) * MDP(35) + (-MDP(4) * t637 * t641 + MDP(5) * t735) * t646; t647 + (t657 + t758) * MDP(17) + (-t470 * t782 - t619 * t681 - t723 * t740 + t649) * MDP(27) + (t439 * qJ(4) - t442 * pkin(3) - t482 * t496 - t485 * t499 - g(1) * (-pkin(3) * t751 + t583) - g(2) * (-pkin(3) * t753 + t582) - g(3) * t737 + t779 * t490) * MDP(21) + (-t470 * t780 + t619 * t736 - t723 * t741 + t648) * MDP(28) + (t568 * t415 + t740 * t475 + (-t666 - t784) * t639 + t653 * t634 + t705) * MDP(34) + (-t496 * t549 + t620 - t651 + t757) * MDP(18) + (t496 * t551 + 0.2e1 * t616 + 0.2e1 * t618 + t658 - t758) * MDP(20) + (t655 + t757) * MDP(16) + (pkin(3) * t473 - qJ(4) * t474 + (t490 - t499) * t551 + (t485 - t779) * t549) * MDP(19) + (t568 * t414 + t740 * t477 + (-t405 + t784) * t634 + t653 * t639 + t674) * MDP(35); (-t627 + t756) * MDP(18) + t455 * MDP(19) + (-t628 ^ 2 - t774) * MDP(20) + (-t490 * t628 + t651) * MDP(21) + (-t551 * t782 - t619 * t640 - t635 * t703) * MDP(27) + (-t551 * t780 + t619 * t635 - t640 * t703) * MDP(28) + (-t640 * t415 + (-t551 * t639 - t634 * t702) * t791 + (t475 * t723 - t726 * t791 - t748) * t635) * MDP(34) + (-t640 * t414 + (t551 * t634 - t639 * t702) * t791 + (t477 * t723 + t727 * t791 - t745) * t635) * MDP(35); (t436 * t723 - t649) * MDP(27) + (t435 * t723 - t648) * MDP(28) + (-pkin(5) * t415 - t436 * t475 + t662 * t634 + t639 * t777 - t705) * MDP(34) + (-pkin(5) * t414 + t796 - t436 * t477 + t662 * t639 + (t405 - t777) * t634) * MDP(35) + t808; t477 * t475 * MDP(29) + (-t475 ^ 2 + t477 ^ 2) * MDP(30) + (t715 + t760) * MDP(31) + (-t592 + t759) * MDP(32) + t426 * MDP(33) + (-g(1) * t506 + t411 * t791 - t432 * t477 + t401) * MDP(34) + (g(1) * t507 + t410 * t791 + t432 * t475) * MDP(35) + (-MDP(32) * t728 + MDP(34) * t682 + MDP(35) * t664) * t639 + (-MDP(31) * t728 + (-qJD(6) * t723 - t427) * MDP(32) + t664 * MDP(34) + (-t402 - t682) * MDP(35)) * t634;];
tau  = t1;
