% Calculate vector of inverse dynamics joint torques for
% S6RRPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPR8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR8_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPRPR8_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:53:39
% EndTime: 2019-03-09 10:53:56
% DurationCPUTime: 13.21s
% Computational Cost: add. (7213->665), mult. (16111->834), div. (0->0), fcn. (11912->12), ass. (0->264)
t663 = sin(pkin(10));
t804 = pkin(8) + qJ(3);
t619 = t804 * t663;
t664 = cos(pkin(10));
t620 = t804 * t664;
t667 = sin(qJ(4));
t815 = cos(qJ(4));
t556 = -t667 * t619 + t815 * t620;
t671 = cos(qJ(2));
t653 = t671 * qJDD(1);
t668 = sin(qJ(2));
t759 = qJD(1) * qJD(2);
t698 = t668 * t759 - t653;
t607 = qJDD(4) + t698;
t660 = pkin(10) + qJ(4);
t651 = sin(t660);
t659 = g(3) * t671;
t669 = sin(qJ(1));
t672 = cos(qJ(1));
t735 = g(1) * t672 + g(2) * t669;
t705 = t735 * t668;
t834 = t705 - t659;
t843 = t556 * t607 + t651 * t834;
t752 = t815 * t664;
t702 = -t667 * t663 + t752;
t771 = qJD(1) * t671;
t747 = qJD(4) * t815;
t766 = qJD(4) * t667;
t825 = -t663 * t766 + t664 * t747;
t779 = -t702 * t771 + t825;
t772 = qJD(1) * t668;
t748 = t663 * t772;
t769 = qJD(2) * t664;
t602 = -t748 + t769;
t751 = t664 * t772;
t770 = qJD(2) * t663;
t603 = t751 + t770;
t538 = -t815 * t602 + t603 * t667;
t744 = qJD(2) * pkin(2) - qJD(3);
t733 = -pkin(7) * t772 + t744;
t693 = pkin(3) * t602 + t733;
t703 = -t667 * t602 - t603 * t815;
t684 = -qJ(5) * t703 + t693;
t816 = pkin(4) + pkin(5);
t448 = -t538 * t816 + t684;
t666 = sin(qJ(6));
t670 = cos(qJ(6));
t479 = t538 * t666 - t670 * t703;
t788 = t672 * t651;
t652 = cos(t660);
t791 = t669 * t652;
t578 = t671 * t788 - t791;
t789 = t671 * t672;
t579 = t651 * t669 + t652 * t789;
t507 = t578 * t670 - t579 * t666;
t715 = t651 * t670 - t652 * t666;
t745 = t671 * t759;
t758 = qJDD(1) * t668;
t699 = t745 + t758;
t757 = t663 * qJDD(2);
t681 = t664 * t699 + t757;
t777 = t699 * t663;
t724 = qJDD(2) * t664 - t777;
t469 = -t602 * t747 + t603 * t766 - t667 * t724 - t815 * t681;
t725 = pkin(2) * t668 - qJ(3) * t671;
t587 = qJD(2) * t725 - qJD(3) * t668;
t726 = pkin(2) * t671 + qJ(3) * t668;
t616 = -pkin(1) - t726;
t535 = qJD(1) * t587 + qJDD(1) * t616;
t575 = -pkin(7) * t698 + qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t487 = t664 * t535 - t663 * t575;
t460 = pkin(3) * t698 - pkin(8) * t681 + t487;
t488 = t663 * t535 + t664 * t575;
t471 = pkin(8) * t724 + t488;
t594 = t616 * qJD(1);
t649 = pkin(7) * t771;
t621 = qJD(2) * qJ(3) + t649;
t548 = t664 * t594 - t621 * t663;
t756 = pkin(3) * t771;
t496 = -pkin(8) * t603 + t548 - t756;
t549 = t663 * t594 + t664 * t621;
t505 = pkin(8) * t602 + t549;
t740 = -t815 * t460 + t667 * t471 + t496 * t766 + t505 * t747;
t713 = qJDD(5) + t740;
t425 = pkin(9) * t469 - t607 * t816 + t713;
t595 = t607 * qJ(5);
t640 = -qJD(4) + t771;
t624 = t640 * qJD(5);
t696 = t667 * t460 + t815 * t471 + t496 * t747 - t505 * t766;
t429 = t595 - t624 + t696;
t470 = -qJD(4) * t703 + t667 * t681 - t815 * t724;
t427 = pkin(9) * t470 + t429;
t743 = t670 * t425 - t666 * t427;
t805 = g(3) * t668;
t790 = t669 * t671;
t576 = t651 * t790 + t652 * t672;
t577 = t652 * t790 - t788;
t824 = t576 * t670 - t577 * t666;
t842 = g(1) * t507 + g(2) * t824 + t448 * t479 + t715 * t805 - t743;
t764 = qJD(6) * t670;
t765 = qJD(6) * t666;
t433 = -t670 * t469 + t666 * t470 + t538 * t764 + t703 * t765;
t596 = -qJDD(6) + t607;
t719 = -t670 * t538 - t666 * t703;
t760 = -qJD(6) - t640;
t836 = t719 * t760;
t841 = (t479 ^ 2 - t719 ^ 2) * MDP(27) + MDP(26) * t479 * t719 - (-t433 + t836) * MDP(28) - t596 * MDP(30);
t840 = t538 ^ 2;
t817 = t703 ^ 2;
t839 = t479 * t760;
t838 = t538 * t640;
t837 = t640 * t703;
t611 = t725 * qJD(1);
t561 = pkin(7) * t748 + t664 * t611;
t794 = t664 * t671;
t710 = pkin(3) * t668 - pkin(8) * t794;
t525 = qJD(1) * t710 + t561;
t592 = t663 * t611;
t795 = t664 * t668;
t796 = t663 * t671;
t701 = -pkin(7) * t795 - pkin(8) * t796;
t545 = qJD(1) * t701 + t592;
t835 = qJD(3) * t752 - t815 * t545 - t619 * t747 + (-qJD(3) * t663 - qJD(4) * t620 - t525) * t667;
t609 = t663 * t815 + t667 * t664;
t591 = t609 * qJD(4);
t690 = t671 * t609;
t778 = -qJD(1) * t690 + t591;
t783 = -qJ(5) * t772 + t835;
t829 = t609 * qJD(3) + qJD(4) * t556 + t525 * t815 - t667 * t545;
t597 = t663 * t756 + t649;
t826 = -qJ(5) * t779 - qJD(5) * t609 - t597;
t455 = t815 * t496 - t667 * t505;
t761 = qJD(5) - t455;
t823 = -pkin(7) * t745 - qJDD(3);
t742 = t666 * t469 + t670 * t470;
t434 = t479 * qJD(6) - t742;
t508 = t578 * t666 + t579 * t670;
t714 = t651 * t666 + t652 * t670;
t718 = t576 * t666 + t577 * t670;
t762 = -pkin(9) * t703 - t761;
t442 = t640 * t816 - t762;
t753 = t666 * t425 + t670 * t427 + t442 * t764;
t819 = -g(1) * t508 - g(2) * t718 - t448 * t719 - t714 * t805 + t753;
t818 = -0.2e1 * pkin(1);
t814 = pkin(3) * t663;
t813 = pkin(4) * t607;
t812 = pkin(7) * t602;
t809 = g(1) * t669;
t806 = g(2) * t672;
t803 = qJ(5) * t538;
t802 = qJDD(2) * pkin(2);
t456 = t667 * t496 + t815 * t505;
t801 = t456 * t640;
t800 = t538 * t703;
t797 = t663 * t668;
t445 = pkin(9) * t538 + t456;
t625 = t640 * qJ(5);
t443 = t445 - t625;
t793 = t666 * t443;
t792 = t668 * t672;
t716 = -t609 * t666 - t670 * t702;
t787 = qJD(6) * t716 + t666 * t778 + t670 * t779;
t544 = t609 * t670 - t666 * t702;
t786 = qJD(6) * t544 + t666 * t779 - t670 * t778;
t785 = -t778 * t816 - t826;
t784 = pkin(4) * t778 + t826;
t782 = pkin(4) * t772 + t829;
t601 = t664 * t616;
t547 = -pkin(8) * t795 + t601 + (-pkin(7) * t663 - pkin(3)) * t671;
t569 = pkin(7) * t794 + t663 * t616;
t554 = -pkin(8) * t797 + t569;
t780 = t667 * t547 + t815 * t554;
t768 = qJD(2) * t668;
t755 = pkin(7) * t768;
t552 = t664 * t587 + t663 * t755;
t647 = pkin(7) * t758;
t776 = -t647 - t659;
t767 = qJD(2) * t671;
t750 = t663 * t767;
t598 = pkin(3) * t750 + pkin(7) * t767;
t612 = pkin(3) * t797 + t668 * pkin(7);
t775 = t672 * pkin(1) + t669 * pkin(7);
t661 = t668 ^ 2;
t774 = -t671 ^ 2 + t661;
t754 = t668 * t816;
t644 = pkin(3) * t664 + pkin(2);
t642 = t668 * t809;
t738 = -g(2) * t792 + t642;
t737 = -g(1) * t576 + g(2) * t578;
t736 = g(1) * t577 - g(2) * t579;
t734 = -t806 + t809;
t482 = -qJ(5) * t671 + t780;
t584 = t702 * t668;
t731 = qJ(5) * t584 - t612;
t729 = t547 * t815 - t667 * t554;
t518 = -pkin(9) * t702 + t556;
t728 = pkin(9) * t779 - qJD(1) * t754 + qJD(6) * t518 - t829;
t555 = t619 * t815 + t667 * t620;
t517 = -t609 * pkin(9) + t555;
t727 = -pkin(9) * t778 - qJD(6) * t517 - t783;
t723 = qJ(5) * t670 - t666 * t816;
t722 = qJ(5) * t666 + t670 * t816;
t432 = t666 * t442 + t670 * t443;
t483 = t671 * pkin(4) - t729;
t461 = t671 * pkin(5) - t584 * pkin(9) + t483;
t583 = t609 * t668;
t463 = pkin(9) * t583 + t482;
t721 = t461 * t670 - t463 * t666;
t720 = t461 * t666 + t463 * t670;
t717 = t670 * t583 - t584 * t666;
t513 = t583 * t666 + t584 * t670;
t709 = qJ(5) * t609 + t644;
t708 = pkin(4) * t652 + qJ(5) * t651 + t644;
t706 = -t647 + t802 + t823;
t704 = -pkin(7) * qJDD(2) + t759 * t818;
t514 = qJD(2) * t710 + t552;
t573 = t663 * t587;
t528 = qJD(2) * t701 + t573;
t700 = t514 * t815 - t667 * t528 - t547 * t766 - t554 * t747;
t697 = t664 * t758 + t757;
t695 = t667 * t514 + t815 * t528 + t547 * t747 - t554 * t766;
t675 = qJD(1) ^ 2;
t694 = pkin(1) * t675 + t735;
t674 = qJD(2) ^ 2;
t692 = pkin(7) * t674 + qJDD(1) * t818 + t806;
t519 = t591 * t668 + t667 * t750 - t752 * t767;
t691 = -qJ(5) * t519 + qJD(5) * t584 - t598;
t689 = g(2) * t668 * t791 - t555 * t607 + (g(1) * t792 - t659) * t652;
t687 = -t671 * t735 - t805;
t686 = -t705 - t802;
t683 = g(1) * t578 + g(2) * t576 + t651 * t805 - t740;
t682 = t706 + t834;
t438 = qJ(5) * t768 - qJD(5) * t671 + t695;
t515 = -pkin(3) * t724 - t706;
t462 = pkin(4) * t538 - t684;
t680 = -t462 * t703 + qJDD(5) - t683;
t678 = g(1) * t579 + g(2) * t577 - t455 * t640 + t652 * t805 - t696;
t677 = -t469 * qJ(5) - qJD(5) * t703 - t515;
t435 = t470 * pkin(4) - t677;
t657 = t672 * pkin(7);
t568 = -pkin(7) * t796 + t601;
t562 = -pkin(7) * t751 + t592;
t553 = -t664 * t755 + t573;
t536 = -pkin(4) * t702 - t709;
t520 = qJD(2) * t690 + t668 * t825;
t502 = t702 * t816 + t709;
t497 = pkin(4) * t583 - t731;
t485 = -pkin(4) * t703 + t803;
t484 = -t583 * t816 + t731;
t457 = t703 * t816 - t803;
t452 = -t625 + t456;
t451 = pkin(4) * t640 + t761;
t450 = pkin(4) * t520 - t691;
t449 = -t469 - t838;
t447 = qJD(6) * t513 - t519 * t666 - t670 * t520;
t446 = qJD(6) * t717 - t519 * t670 + t520 * t666;
t441 = -t520 * t816 + t691;
t440 = -pkin(4) * t768 - t700;
t437 = pkin(9) * t520 + t438;
t436 = t519 * pkin(9) - qJD(2) * t754 - t700;
t431 = t442 * t670 - t793;
t430 = t713 - t813;
t428 = -t470 * t816 + t677;
t1 = [((qJD(6) * t721 + t436 * t666 + t437 * t670) * t760 + t720 * t596 - (-t443 * t765 + t753) * t671 + t432 * t768 + t441 * t479 + t484 * t433 + t428 * t513 + t448 * t446 + g(1) * t824 - g(2) * t507) * MDP(32) + (-t429 * t583 + t430 * t584 - t438 * t538 - t440 * t703 - t451 * t519 - t452 * t520 - t469 * t483 - t470 * t482 + t738) * MDP(23) + (t469 * t671 + t519 * t640 + t584 * t607 - t703 * t768) * MDP(17) + (-t429 * t671 - t435 * t584 - t438 * t640 + t450 * t703 + t452 * t768 + t462 * t519 + t469 * t497 + t482 * t607 - t737) * MDP(24) + (-t469 * t584 + t519 * t703) * MDP(15) + (t469 * t583 - t470 * t584 + t519 * t538 + t520 * t703) * MDP(16) + (-t456 * t768 - t612 * t469 + t515 * t584 + t519 * t693 - t598 * t703 - t607 * t780 + t640 * t695 + t671 * t696 + t737) * MDP(21) + (t455 * t768 + t612 * t470 + t515 * t583 - t520 * t693 + t598 * t538 + t607 * t729 - t640 * t700 + t671 * t740 + t736) * MDP(20) + (t488 * t569 + t549 * t553 + t487 * t568 + t548 * t552 - g(1) * t657 - g(2) * (t672 * t726 + t775) - t616 * t809 + (-t668 * t706 - t733 * t767) * pkin(7)) * MDP(14) + (t433 * t717 - t434 * t513 - t446 * t719 - t447 * t479) * MDP(27) + (-t735 * t663 + (-pkin(7) * t724 - t706 * t663 + (qJD(1) * t568 + t548) * qJD(2)) * t668 + (-t552 * qJD(1) - t568 * qJDD(1) - t487 + t734 * t664 + (-t663 * t733 - t812) * qJD(2)) * t671) * MDP(11) + (-t735 * t664 + (-t706 * t664 + (-qJD(1) * t569 - t549) * qJD(2) + t697 * pkin(7)) * t668 + (t553 * qJD(1) + t569 * qJDD(1) + t488 - t734 * t663 + (-t733 * t664 + (t603 + t751) * pkin(7)) * qJD(2)) * t671) * MDP(12) + (-(t436 * t670 - t437 * t666) * t760 - t721 * t596 + t743 * t671 - t431 * t768 + t441 * t719 + t484 * t434 - t428 * t717 + t448 * t447 + g(1) * t718 - g(2) * t508 + (-t432 * t671 + t720 * t760) * qJD(6)) * MDP(31) + (-t434 * t671 + t447 * t760 - t596 * t717 + t719 * t768) * MDP(29) + (-t596 * t671 + t760 * t768) * MDP(30) + (t433 * t671 - t446 * t760 - t479 * t768 - t513 * t596) * MDP(28) + (t429 * t482 + t452 * t438 + t435 * t497 + t462 * t450 + t430 * t483 + t451 * t440 - g(1) * (-pkin(4) * t577 - qJ(5) * t576 + t672 * t814 + t657) - g(2) * (pkin(4) * t579 + qJ(5) * t578 + t644 * t789 + t792 * t804 + t775) + (-g(1) * (-t644 * t671 - t668 * t804 - pkin(1)) - g(2) * t814) * t669) * MDP(25) + qJDD(1) * MDP(1) + (t704 * t668 + (-t692 + t809) * t671) * MDP(9) + (t668 * t692 + t671 * t704 - t642) * MDP(10) + t734 * MDP(2) + t735 * MDP(3) + 0.2e1 * (t653 * t668 - t759 * t774) * MDP(5) + (t553 * t602 - t569 * t777 - t552 * t603 + (-qJDD(2) * t568 - t488 * t668 - t549 * t767) * t663 + (t569 * qJDD(2) - t487 * t668 - t548 * t767 - t568 * t699) * t664 + t738) * MDP(13) + (-t607 * t671 - t640 * t768) * MDP(19) + (t470 * t671 + t520 * t640 - t538 * t768 - t583 * t607) * MDP(18) + (t430 * t671 + t435 * t583 + t440 * t640 + t450 * t538 - t451 * t768 + t462 * t520 + t470 * t497 - t483 * t607 + t736) * MDP(22) + (qJDD(2) * t668 + t671 * t674) * MDP(6) + (qJDD(2) * t671 - t668 * t674) * MDP(7) + (qJDD(1) * t661 + 0.2e1 * t668 * t745) * MDP(4) + (t433 * t513 + t446 * t479) * MDP(26); (-MDP(4) * t668 * t671 + MDP(5) * t774) * t675 + (-(t517 * t670 - t518 * t666) * t596 + t502 * t434 - t428 * t716 - (t666 * t727 - t670 * t728) * t760 + t785 * t719 + t786 * t448 + t834 * t714) * MDP(31) + ((t517 * t666 + t518 * t670) * t596 + t502 * t433 + t428 * t544 - (t666 * t728 + t670 * t727) * t760 + t785 * t479 + t787 * t448 + t834 * t715) * MDP(32) + (-t469 * t609 - t703 * t779) * MDP(15) + (t433 * t716 - t434 * t544 - t479 * t786 - t719 * t787) * MDP(27) + (t429 * t702 + t430 * t609 + t451 * t779 - t452 * t778 - t469 * t555 - t470 * t556 - t538 * t783 - t703 * t782 + t687) * MDP(23) + (-t469 * t702 - t470 * t609 - t538 * t779 + t703 * t778) * MDP(16) + (-t644 * t470 - t515 * t702 - t597 * t538 + t640 * t829 - t693 * t778 + t689) * MDP(20) + (-t435 * t702 + t462 * t778 + t470 * t536 + t538 * t784 + t640 * t782 + t689) * MDP(22) + (t607 * t702 + t640 * t778) * MDP(18) + (-t725 * t664 * qJDD(1) + (-t706 + t686 + t659) * t663 + ((-qJ(3) * t769 + t549) * t668 + (-pkin(7) * t603 - t562 + (t733 - t744) * t664) * t671) * qJD(1)) * MDP(12) + (t663 * qJ(3) * t653 - pkin(2) * t777 + (t682 + t802) * t664 + ((-qJ(3) * t770 - t548) * t668 + (t812 + t561 + (qJD(3) + t733) * t663) * t671) * qJD(1)) * MDP(11) + (t733 * t649 - t548 * t561 - t549 * t562 + (-t548 * t663 + t549 * t664) * qJD(3) + t682 * pkin(2) + (-t487 * t663 + t488 * t664 + t687) * qJ(3)) * MDP(14) + (MDP(17) * t703 + t538 * MDP(18) + t640 * MDP(19) - t455 * MDP(20) + t456 * MDP(21) + t451 * MDP(22) - t452 * MDP(24) + t479 * MDP(28) - MDP(29) * t719 - MDP(30) * t760 + t431 * MDP(31) - t432 * MDP(32)) * t772 + (-t596 * t716 + t760 * t786) * MDP(29) + (-t544 * t596 - t760 * t787) * MDP(28) + (t429 * t556 + t430 * t555 + t435 * t536 + t784 * t462 + t783 * t452 + t782 * t451 + (-g(3) * t708 - t735 * t804) * t671 + (-g(3) * t804 + t708 * t735) * t668) * MDP(25) + qJDD(2) * MDP(8) + MDP(7) * t653 + MDP(6) * t758 + (t433 * t544 + t479 * t787) * MDP(26) + (t607 * t609 - t640 * t779) * MDP(17) + (t805 + (-pkin(7) * qJDD(1) + t694) * t671) * MDP(10) + (t668 * t694 + t776) * MDP(9) + (t561 * t603 - t562 * t602 + (qJ(3) * t724 + qJD(3) * t602 + t548 * t771 + t488) * t664 + (qJ(3) * t681 + qJD(3) * t603 + t549 * t771 - t487) * t663 + t687) * MDP(13) + (t644 * t469 + t515 * t609 + t597 * t703 + t640 * t835 - t779 * t693 - t843) * MDP(21) + (-t435 * t609 - t462 * t779 + t469 * t536 - t640 * t783 + t703 * t784 + t843) * MDP(24); (-t603 * t771 - t724) * MDP(11) + ((-t602 + t769) * t771 + t697) * MDP(12) + (-t602 ^ 2 - t603 ^ 2) * MDP(13) + (t548 * t603 - t549 * t602 + t686 - t776 - t823) * MDP(14) + (-t817 - t840) * MDP(23) + (t451 * t703 + t452 * t538 + t435 - t834) * MDP(25) + (-t434 + t839) * MDP(31) + (-t433 - t836) * MDP(32) + (-MDP(21) + MDP(24)) * (t469 - t838) + (MDP(20) + MDP(22)) * (t470 + t837); -MDP(15) * t800 + (t817 - t840) * MDP(16) + t449 * MDP(17) + (-t470 + t837) * MDP(18) + t607 * MDP(19) + (-t693 * t703 + t683 - t801) * MDP(20) + (-t538 * t693 + t678) * MDP(21) + (-t485 * t538 - t680 - t801 + 0.2e1 * t813) * MDP(22) + (pkin(4) * t469 - qJ(5) * t470 - (t452 - t456) * t703 + (t451 - t761) * t538) * MDP(23) + (-t462 * t538 - t485 * t703 + 0.2e1 * t595 - 0.2e1 * t624 - t678) * MDP(24) + (t429 * qJ(5) - t430 * pkin(4) - t462 * t485 - t451 * t456 - g(1) * (-pkin(4) * t578 + qJ(5) * t579) - g(2) * (-pkin(4) * t576 + qJ(5) * t577) - (-pkin(4) * t651 + qJ(5) * t652) * t805 + t761 * t452) * MDP(25) + (t434 + t839) * MDP(29) + (t722 * t596 - t457 * t719 - (-t670 * t445 + t666 * t762) * t760 + (t723 * t760 + t432) * qJD(6) + t842) * MDP(31) + (t723 * t596 - t457 * t479 - (t666 * t445 + t670 * t762) * t760 + (-t722 * t760 - t793) * qJD(6) + t819) * MDP(32) - t841; (-t607 - t800) * MDP(22) + t449 * MDP(23) + (-t640 ^ 2 - t817) * MDP(24) + (t452 * t640 + t680 - t813) * MDP(25) + (-t596 * t670 + t703 * t719) * MDP(31) + (t479 * t703 + t596 * t666) * MDP(32) - (MDP(31) * t666 + MDP(32) * t670) * t760 ^ 2; (t742 - t839) * MDP(29) + (-t432 * t760 - t842) * MDP(31) + (-t431 * t760 - t819) * MDP(32) + (-MDP(29) * t479 - MDP(31) * t432 + MDP(32) * t793) * qJD(6) + t841;];
tau  = t1;
