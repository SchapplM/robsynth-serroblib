% Calculate vector of inverse dynamics joint torques for
% S6RRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRPRRR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:53:43
% EndTime: 2019-03-09 13:53:58
% DurationCPUTime: 10.09s
% Computational Cost: add. (6980->540), mult. (15034->683), div. (0->0), fcn. (11056->12), ass. (0->251)
t645 = sin(qJ(4));
t650 = cos(qJ(4));
t651 = cos(qJ(2));
t746 = qJD(1) * t651;
t646 = sin(qJ(2));
t747 = qJD(1) * t646;
t550 = t645 * t747 + t650 * t746;
t553 = -t645 * t746 + t650 * t747;
t644 = sin(qJ(5));
t649 = cos(qJ(5));
t494 = t550 * t644 - t649 * t553;
t634 = qJDD(2) - qJDD(4);
t627 = -qJDD(5) + t634;
t643 = sin(qJ(6));
t648 = cos(qJ(6));
t567 = t645 * t646 + t650 * t651;
t673 = t567 * qJD(4);
t733 = qJD(1) * qJD(2);
t723 = t646 * t733;
t731 = qJDD(1) * t651;
t816 = t723 - t731;
t722 = t651 * t733;
t732 = qJDD(1) * t646;
t817 = t722 + t732;
t482 = -qJD(1) * t673 + t645 * t816 + t650 * t817;
t742 = qJD(4) * t650;
t743 = qJD(4) * t645;
t744 = qJD(2) * t651;
t835 = t645 * t744 + t646 * t742 - t651 * t743;
t483 = qJD(1) * t835 + qJDD(1) * t567 - t650 * t723;
t740 = qJD(5) * t649;
t741 = qJD(5) * t644;
t675 = -t649 * t482 + t644 * t483 + t550 * t740 + t553 * t741;
t736 = qJD(6) * t648;
t635 = qJD(2) - qJD(4);
t832 = -qJD(5) + t635;
t728 = -t643 * t627 - t648 * t675 - t736 * t832;
t737 = qJD(6) * t643;
t431 = t494 * t737 + t728;
t429 = t431 * t643;
t430 = t431 * t648;
t486 = -t494 * t648 - t643 * t832;
t593 = t648 * t627;
t432 = qJD(6) * t486 - t643 * t675 + t593;
t660 = qJD(5) * t494 - t482 * t644 - t649 * t483;
t442 = qJDD(6) - t660;
t439 = t643 * t442;
t440 = t648 * t442;
t484 = -t494 * t643 + t648 * t832;
t799 = t649 * t550 + t553 * t644;
t813 = qJD(6) + t799;
t842 = t813 * t643;
t822 = t799 * t648;
t847 = t736 + t822;
t848 = -t627 * MDP(26) + (t494 * t832 + t660) * MDP(25) + (t486 * t847 + t429) * MDP(29) + (t486 * t494 + t813 * t847 + t439) * MDP(31) + (t494 ^ 2 - t799 ^ 2) * MDP(23) - (t799 * t832 + t675) * MDP(24) + (-MDP(22) * t799 + MDP(33) * t813) * t494 - (t643 * t432 + t484 * t847 + t486 * t842 - t430) * MDP(30) - (t484 * t494 + t813 * t842 - t440) * MDP(32);
t623 = pkin(7) * t747;
t844 = -pkin(8) * t747 + qJD(3) + t623;
t843 = (-t550 * t635 + t482) * MDP(17) + t553 * t550 * MDP(15) - (t550 ^ 2 - t553 ^ 2) * MDP(16) - t634 * MDP(19) - (t553 * t635 + t483) * MDP(18) + t848;
t624 = pkin(7) * t746;
t577 = -pkin(8) * t746 + t624;
t653 = -pkin(2) - pkin(3);
t719 = -qJ(3) * t645 + t650 * t653;
t838 = qJD(4) * t719 - t645 * t577 + t650 * t844;
t579 = qJ(3) * t650 + t645 * t653;
t837 = qJD(4) * t579 + t650 * t577 + t645 * t844;
t607 = pkin(7) * t722;
t618 = pkin(7) * t732;
t721 = qJDD(3) + t607 + t618;
t507 = -pkin(8) * t817 + t653 * qJDD(2) + t721;
t619 = pkin(7) * t731;
t636 = qJDD(2) * qJ(3);
t637 = qJD(2) * qJD(3);
t532 = -pkin(7) * t723 + t619 + t636 + t637;
t509 = pkin(8) * t816 + t532;
t727 = t653 * qJD(2);
t534 = t727 + t844;
t638 = qJD(2) * qJ(3);
t554 = t577 + t638;
t688 = t534 * t645 + t554 * t650;
t834 = -t688 * qJD(4) + t650 * t507 - t645 * t509;
t435 = -pkin(4) * t634 - pkin(9) * t482 + t834;
t677 = -t645 * t507 - t650 * t509 - t534 * t742 + t554 * t743;
t438 = -pkin(9) * t483 - t677;
t716 = t650 * t534 - t554 * t645;
t788 = pkin(9) * t553;
t474 = t716 - t788;
t471 = -pkin(4) * t635 + t474;
t789 = pkin(9) * t550;
t475 = t688 - t789;
t704 = -t649 * t435 + t644 * t438 + t471 * t741 + t475 * t740;
t420 = pkin(5) * t627 + t704;
t774 = t475 * t649;
t450 = t471 * t644 + t774;
t448 = -pkin(10) * t832 + t450;
t555 = -qJD(1) * pkin(1) - pkin(2) * t746 - qJ(3) * t747;
t530 = pkin(3) * t746 - t555;
t501 = pkin(4) * t550 + t530;
t453 = pkin(5) * t799 + pkin(10) * t494 + t501;
t427 = t448 * t648 + t453 * t643;
t775 = t475 * t644;
t449 = t471 * t649 - t775;
t447 = pkin(5) * t832 - t449;
t833 = -t420 * t643 + t427 * t494 - t447 * t736;
t824 = t447 * t799;
t820 = -t789 + t837;
t819 = -t788 + t838;
t647 = sin(qJ(1));
t685 = t645 * t651 - t646 * t650;
t541 = t685 * t647;
t652 = cos(qJ(1));
t764 = t651 * t652;
t767 = t646 * t652;
t543 = t645 * t764 - t650 * t767;
t818 = -g(1) * t543 - g(2) * t541 - g(3) * t567 + t530 * t553 - t834;
t763 = qJ(4) + qJ(5);
t630 = sin(t763);
t724 = cos(t763);
t696 = t646 * t724;
t798 = -t651 * t630 + t696;
t526 = t798 * t647;
t528 = t630 * t764 - t652 * t696;
t668 = t646 * t630 + t651 * t724;
t672 = g(1) * t528 - g(2) * t526 + g(3) * t668;
t426 = -t448 * t643 + t453 * t648;
t811 = t426 * t494 + t447 * t737;
t810 = -t643 * t672 - t833;
t807 = t501 * t494 + t672 - t704;
t529 = t668 * t652;
t527 = t668 * t647;
t697 = g(2) * t527 + g(3) * t798;
t703 = -t644 * t435 - t649 * t438 - t471 * t740 + t475 * t741;
t806 = g(1) * t529 + t799 * t501 + t697 + t703;
t804 = pkin(5) * t494 - pkin(10) * t799;
t750 = t651 * pkin(2) + t646 * qJ(3);
t583 = -pkin(1) - t750;
t614 = qJ(3) * t746;
t540 = t653 * t747 + t614;
t790 = pkin(4) * t553;
t506 = t540 - t790;
t574 = -pkin(4) + t719;
t754 = t644 * t574 + t649 * t579;
t515 = -pkin(10) + t754;
t801 = t813 * (qJD(6) * t515 + t506 + t804);
t791 = pkin(7) - pkin(8);
t586 = t791 * t646;
t587 = t791 * t651;
t753 = t645 * t586 + t650 * t587;
t781 = g(2) * t652;
t785 = g(1) * t647;
t797 = -t781 + t785;
t782 = g(2) * t647;
t784 = g(1) * t652;
t698 = t782 + t784;
t779 = pkin(7) * qJDD(2);
t794 = (qJD(1) * t583 + t555) * qJD(2) - t779;
t745 = qJD(2) * t646;
t518 = -t650 * t745 + t835;
t576 = t791 * t745;
t578 = qJD(2) * t587;
t676 = -t650 * t576 + t645 * t578 + t586 * t742 - t587 * t743;
t462 = -pkin(9) * t518 + t676;
t519 = qJD(2) * t567 - t673;
t661 = -qJD(4) * t753 + t576 * t645 + t650 * t578;
t463 = -pkin(9) * t519 + t661;
t714 = t650 * t586 - t587 * t645;
t491 = pkin(9) * t685 + t714;
t492 = -pkin(9) * t567 + t753;
t689 = t491 * t649 - t492 * t644;
t423 = qJD(5) * t689 + t462 * t649 + t463 * t644;
t510 = t649 * t567 - t644 * t685;
t457 = -qJD(5) * t510 - t518 * t644 + t519 * t649;
t511 = -t567 * t644 - t649 * t685;
t563 = t651 * pkin(3) - t583;
t520 = pkin(4) * t567 + t563;
t461 = pkin(5) * t510 - pkin(10) * t511 + t520;
t465 = t491 * t644 + t492 * t649;
t709 = -pkin(10) * t627 + qJD(6) * t453 - t703;
t792 = t420 * t511 - t465 * t442 + t447 * t457 - (qJD(6) * t461 + t423) * t813 - t510 * t709 + t784;
t787 = g(1) * t527;
t641 = qJDD(1) * pkin(1);
t777 = qJDD(2) * pkin(2);
t776 = t447 * t511;
t655 = qJD(1) ^ 2;
t766 = t646 * t655;
t687 = t574 * t649 - t579 * t644;
t759 = qJD(5) * t687 - t644 * t820 + t649 * t819;
t758 = qJD(5) * t754 + t644 * t819 + t649 * t820;
t566 = t644 * t645 - t649 * t650;
t757 = t832 * t566;
t569 = t644 * t650 + t645 * t649;
t756 = t832 * t569;
t629 = t646 * qJD(3);
t751 = qJ(3) * t744 + t629;
t639 = t646 ^ 2;
t640 = t651 ^ 2;
t749 = t639 - t640;
t739 = qJD(6) * t448;
t738 = qJD(6) * t494;
t730 = t651 * t766;
t720 = -qJD(2) * pkin(2) + qJD(3);
t616 = pkin(4) * t644 + pkin(10);
t706 = qJD(6) * t616 + t790 - t804;
t705 = t635 ^ 2;
t702 = t646 * t727;
t451 = t474 * t644 + t774;
t701 = pkin(4) * t741 - t451;
t452 = t474 * t649 - t775;
t700 = -pkin(4) * t740 + t452;
t654 = qJD(2) ^ 2;
t699 = pkin(7) * t654 + t781;
t695 = pkin(2) * t646 - qJ(3) * t651;
t694 = t461 * t442 + t787;
t692 = pkin(2) * t731 + qJ(3) * t817 + qJD(1) * t629 + t641;
t691 = -t739 - t781;
t690 = -t616 * t442 + t824;
t580 = t623 + t720;
t584 = t624 + t638;
t686 = t580 * t651 - t584 * t646;
t684 = g(1) * t767 - g(3) * t651 + t646 * t782 - t618;
t683 = qJD(6) * t569 + t747;
t681 = -0.2e1 * pkin(1) * t733 - t779;
t680 = t457 * t648 - t511 * t737;
t678 = -qJDD(3) + t684;
t524 = t702 + t751;
t670 = -t699 + 0.2e1 * t641;
t667 = -t420 + t672;
t666 = t697 - t709;
t478 = pkin(4) * t518 + t524;
t665 = -t515 * t442 - t759 * t813 - t824;
t505 = pkin(2) * t723 - t692;
t545 = pkin(2) * t745 - t751;
t659 = -qJD(1) * t545 - qJDD(1) * t583 - t505 - t699;
t487 = pkin(3) * t731 + qJD(1) * t702 + t692;
t459 = pkin(4) * t483 + t487;
t542 = t567 * t647;
t544 = t567 * t652;
t657 = g(1) * t544 + g(2) * t542 - g(3) * t685 + t530 * t550 + t677;
t539 = t721 - t777;
t656 = qJD(2) * t686 + t532 * t651 + t539 * t646 - t698;
t617 = -pkin(4) * t649 - pkin(5);
t609 = t651 * t785;
t573 = pkin(2) * t747 - t614;
t514 = pkin(5) - t687;
t513 = t529 * t648 - t643 * t647;
t512 = -t529 * t643 - t647 * t648;
t458 = qJD(5) * t511 + t649 * t518 + t519 * t644;
t425 = pkin(5) * t458 - pkin(10) * t457 + t478;
t424 = qJD(5) * t465 + t462 * t644 - t463 * t649;
t422 = -pkin(5) * t660 + pkin(10) * t675 + t459;
t421 = t648 * t422;
t1 = [(-t457 * t799 + t458 * t494 + t510 * t675 + t511 * t660) * MDP(23) + (-g(2) * t529 + t424 * t832 + t458 * t501 + t459 * t510 + t478 * t799 - t520 * t660 - t627 * t689 + t787) * MDP(27) + (g(1) * t526 + g(2) * t528 + t423 * t832 + t457 * t501 + t459 * t511 + t465 * t627 - t478 * t494 - t520 * t675) * MDP(28) + (t458 * t832 + t510 * t627) * MDP(25) + (-t457 * t832 - t511 * t627) * MDP(24) + (-g(2) * t513 + t421 * t510 + t424 * t484 + t426 * t458 - t689 * t432 + (t425 * t813 + (-t448 * t510 - t465 * t813 + t776) * qJD(6) + t694) * t648 + t792 * t643) * MDP(34) + (t431 * t510 + t440 * t511 + t458 * t486 + t680 * t813) * MDP(31) + (-t511 * t439 - t432 * t510 - t458 * t484 + (-t457 * t643 - t511 * t736) * t813) * MDP(32) + (t442 * t510 + t458 * t813) * MDP(33) + (-g(2) * t512 + t424 * t486 - t427 * t458 - t689 * t431 + (-(-qJD(6) * t465 + t425) * t813 - (t422 - t739) * t510 - qJD(6) * t776 - t694) * t643 + t792 * t648) * MDP(35) + (g(1) * t542 - g(2) * t544 + t563 * t483 + t487 * t567 + t530 * t518 + t524 * t550 - t634 * t714 - t635 * t661) * MDP(20) + (-t457 * t494 - t511 * t675) * MDP(22) + (-g(1) * t541 + g(2) * t543 + t563 * t482 - t487 * t685 + t530 * t519 + t524 * t553 + t634 * t753 + t635 * t676) * MDP(21) + (-t482 * t685 + t519 * t553) * MDP(15) + qJDD(1) * MDP(1) + t698 * MDP(3) + (t656 * pkin(7) + t555 * t545 + (t505 - t797) * t583) * MDP(14) + (t430 * t511 + t486 * t680) * MDP(29) + ((-t484 * t648 - t486 * t643) * t457 + (-t429 - t432 * t648 + (t484 * t643 - t486 * t648) * qJD(6)) * t511) * MDP(30) + (t646 * t681 + t651 * t670 + t609) * MDP(9) + (-t482 * t567 + t483 * t685 - t518 * t553 - t519 * t550) * MDP(16) + (-t519 * t635 + t634 * t685) * MDP(17) + (qJDD(1) * t639 + 0.2e1 * t646 * t722) * MDP(4) + 0.2e1 * (t646 * t731 - t733 * t749) * MDP(5) + (t681 * t651 + (-t670 - t785) * t646) * MDP(10) + ((t639 + t640) * qJDD(1) * pkin(7) + t656) * MDP(12) + (qJDD(2) * t651 - t646 * t654) * MDP(7) + (qJDD(2) * t646 + t651 * t654) * MDP(6) + (t646 * t794 + t659 * t651 + t609) * MDP(11) + (-t794 * t651 + (t659 + t785) * t646) * MDP(13) + t797 * MDP(2) + (t518 * t635 + t567 * t634) * MDP(18); (t514 * t431 + t758 * t486 + t665 * t648 + (t672 + t801) * t643 + t833) * MDP(35) + (-t540 * t550 - t719 * t634 + t635 * t837 + t818) * MDP(20) + (-t540 * t553 + t579 * t634 + t635 * t838 - t657) * MDP(21) + (t514 * t432 + t758 * t484 + t665 * t643 + (-t667 - t801) * t648 - t811) * MDP(34) + (-t506 * t799 - t627 * t687 + t758 * t832 - t807) * MDP(27) + (t494 * t506 + t627 * t754 + t759 * t832 - t806) * MDP(28) + qJDD(2) * MDP(8) + (t619 + 0.2e1 * t636 + 0.2e1 * t637 + (qJD(1) * t573 - g(3)) * t646 + (qJD(1) * t555 - t698) * t651) * MDP(13) + (g(3) * t646 - t619 + (pkin(1) * t655 + t698) * t651) * MDP(10) - t843 + (0.2e1 * t777 + (-t555 * t646 + t573 * t651) * qJD(1) + t678) * MDP(11) + (pkin(1) * t766 + t684) * MDP(9) + MDP(7) * t731 + MDP(6) * t732 - MDP(4) * t730 + (-t695 * qJDD(1) + ((t584 - t638) * t646 + (-t580 + t720) * t651) * qJD(1)) * MDP(12) + t749 * MDP(5) * t655 + (-pkin(7) * qJD(1) * t686 - t539 * pkin(2) - g(3) * t750 + t532 * qJ(3) + t584 * qJD(3) - t555 * t573 + t695 * t698) * MDP(14); (-qJDD(2) - t730) * MDP(11) + MDP(12) * t732 + (-t639 * t655 - t654) * MDP(13) + (-qJD(2) * t584 + t555 * t747 + t607 - t678 - t777) * MDP(14) + (-t550 * t747 - t634 * t650 - t645 * t705) * MDP(20) + (-t553 * t747 + t634 * t645 - t650 * t705) * MDP(21) + (t566 * t627 - t747 * t799 - t756 * t832) * MDP(27) + (t494 * t747 + t569 * t627 + t757 * t832) * MDP(28) + (-t569 * t439 + t566 * t432 - t756 * t484 + (-t643 * t757 - t648 * t683) * t813) * MDP(34) + (-t569 * t440 + t566 * t431 - t756 * t486 + (t643 * t683 - t648 * t757) * t813) * MDP(35); (-t635 * t688 - t818) * MDP(20) + (-t635 * t716 + t657) * MDP(21) + (-t451 * t832 + (-t553 * t799 - t627 * t649 + t741 * t832) * pkin(4) + t807) * MDP(27) + (-t452 * t832 + (t494 * t553 + t627 * t644 + t740 * t832) * pkin(4) + t806) * MDP(28) + (t617 * t432 + t701 * t484 + (t700 * t813 + t690) * t643 + (-t706 * t813 + t667) * t648 + t811) * MDP(34) + (t617 * t431 + t690 * t648 + t701 * t486 + (t643 * t706 + t648 * t700) * t813 + t810) * MDP(35) + t843; (-t450 * t832 + t807) * MDP(27) + (-t449 * t832 + t806) * MDP(28) + (-pkin(5) * t432 - t450 * t484 + (-pkin(10) * t442 + t449 * t813 + t824) * t643 + ((-pkin(10) * qJD(6) + t804) * t813 + t667) * t648 + t811) * MDP(34) + (-pkin(5) * t431 + (t449 * t648 - t643 * t804) * t813 - t450 * t486 + t447 * t822 + (t737 * t813 - t440) * pkin(10) + t810) * MDP(35) + t848; t486 * t484 * MDP(29) + (-t484 ^ 2 + t486 ^ 2) * MDP(30) + (t484 * t813 + t728) * MDP(31) + (t486 * t813 - t593) * MDP(32) + t442 * MDP(33) + (-g(1) * t512 + t427 * t813 - t447 * t486 + t421) * MDP(34) + (g(1) * t513 + t426 * t813 + t447 * t484) * MDP(35) + (MDP(32) * t738 + MDP(34) * t691 + MDP(35) * t666) * t648 + (MDP(31) * t738 + (qJD(6) * t832 + t675) * MDP(32) + t666 * MDP(34) + (-t422 - t691) * MDP(35)) * t643;];
tau  = t1;
