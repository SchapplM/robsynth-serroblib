% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRPRR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:39:51
% EndTime: 2019-03-09 18:40:16
% DurationCPUTime: 16.59s
% Computational Cost: add. (16484->601), mult. (42947->826), div. (0->0), fcn. (35447->12), ass. (0->264)
t705 = cos(qJ(6));
t782 = qJD(6) * t705;
t706 = cos(qJ(5));
t700 = cos(pkin(6));
t792 = qJD(1) * t700;
t688 = qJD(2) + t792;
t707 = cos(qJ(3));
t703 = sin(qJ(3));
t704 = sin(qJ(2));
t698 = sin(pkin(6));
t793 = qJD(1) * t698;
t769 = t704 * t793;
t745 = t703 * t769;
t622 = -t688 * t707 + t745;
t624 = t688 * t703 + t707 * t769;
t697 = sin(pkin(12));
t699 = cos(pkin(12));
t725 = t622 * t699 + t624 * t697;
t586 = t706 * t725;
t588 = t622 * t697 - t624 * t699;
t702 = sin(qJ(5));
t532 = -t588 * t702 + t586;
t860 = t532 * t705;
t869 = t782 + t860;
t708 = cos(qJ(2));
t791 = qJD(1) * t708;
t768 = t698 * t791;
t868 = qJD(3) - t768;
t777 = qJD(1) * qJD(2);
t761 = t698 * t777;
t742 = t708 * t761;
t786 = qJD(3) * t707;
t598 = -qJD(3) * t745 + t688 * t786 + t707 * t742;
t788 = qJD(2) * t708;
t765 = t703 * t788;
t787 = qJD(3) * t703;
t599 = (t704 * t786 + t765) * t793 + t688 * t787;
t556 = t598 * t699 - t599 * t697;
t726 = -t598 * t697 - t599 * t699;
t785 = qJD(5) * t702;
t478 = -qJD(5) * t586 + t556 * t706 + t588 * t785 + t702 * t726;
t670 = -qJD(5) - t868;
t701 = sin(qJ(6));
t743 = t704 * t761;
t773 = t478 * t705 - t670 * t782 + t701 * t743;
t783 = qJD(6) * t701;
t840 = -t588 * t706 - t702 * t725;
t457 = -t783 * t840 + t773;
t455 = t457 * t701;
t525 = -t670 * t701 + t705 * t840;
t756 = t478 * t701 - t705 * t743;
t458 = qJD(6) * t525 + t756;
t479 = qJD(5) * t840 + t702 * t556 - t706 * t726;
t821 = t840 * t701;
t523 = t670 * t705 + t821;
t859 = -qJD(6) - t532;
t750 = t859 * t701;
t475 = t701 * t479;
t808 = -t782 * t859 + t475;
t823 = t532 * t670;
t825 = t840 * t670;
t826 = t525 * t840;
t867 = MDP(24) * t743 + (-t479 - t825) * MDP(23) - t532 ^ 2 * MDP(21) + (MDP(20) * t532 + MDP(21) * t840 + MDP(31) * t859) * t840 + (t478 - t823) * MDP(22) + (t525 * t869 + t455) * MDP(27) + (-t859 * t860 + t808 - t826) * MDP(29) + (t457 * t705 - t458 * t701 - t523 * t869 + t525 * t750) * MDP(28);
t776 = pkin(1) * t792;
t646 = -pkin(8) * t769 + t708 * t776;
t721 = (pkin(2) * t704 - pkin(9) * t708) * t698;
t647 = qJD(1) * t721;
t751 = -t646 * t703 + t647 * t707;
t828 = -qJ(4) - pkin(9);
t760 = qJD(3) * t828;
t811 = t707 * t708;
t865 = -(pkin(3) * t704 - qJ(4) * t811) * t793 - t751 - qJD(4) * t703 + t707 * t760;
t744 = t703 * t768;
t799 = t646 * t707 + t647 * t703;
t864 = -qJ(4) * t744 - qJD(4) * t707 - t703 * t760 + t799;
t685 = t704 * t776;
t649 = pkin(8) * t768 + t685;
t612 = pkin(9) * t688 + t649;
t642 = (-pkin(2) * t708 - pkin(9) * t704 - pkin(1)) * t698;
t618 = qJD(1) * t642;
t581 = -t612 * t703 + t618 * t707;
t554 = -qJ(4) * t624 + t581;
t541 = pkin(3) * t868 + t554;
t582 = t612 * t707 + t618 * t703;
t555 = -qJ(4) * t622 + t582;
t547 = t697 * t555;
t501 = t541 * t699 - t547;
t844 = pkin(10) * t588;
t482 = pkin(4) * t868 + t501 + t844;
t812 = t699 * t555;
t502 = t541 * t697 + t812;
t843 = pkin(10) * t725;
t484 = t502 - t843;
t445 = t482 * t706 - t484 * t702;
t443 = pkin(5) * t670 - t445;
t861 = t443 * t532;
t662 = t697 * t707 + t699 * t703;
t849 = t868 * t662;
t661 = -t697 * t703 + t699 * t707;
t856 = t868 * t661;
t490 = pkin(5) * t840 + pkin(11) * t532;
t611 = -t688 * pkin(2) - t646;
t587 = pkin(3) * t622 + qJD(4) + t611;
t536 = pkin(4) * t725 + t587;
t648 = qJD(2) * t721;
t632 = qJD(1) * t648;
t814 = t698 * t704;
t689 = pkin(8) * t814;
t830 = pkin(1) * t708;
t650 = (t700 * t830 - t689) * qJD(2);
t633 = qJD(1) * t650;
t711 = -qJD(3) * t582 + t632 * t707 - t633 * t703;
t503 = pkin(3) * t743 - qJ(4) * t598 - qJD(4) * t624 + t711;
t720 = -t612 * t787 + t618 * t786 + t632 * t703 + t633 * t707;
t508 = -qJ(4) * t599 - qJD(4) * t622 + t720;
t465 = t503 * t699 - t508 * t697;
t451 = pkin(4) * t743 - pkin(10) * t556 + t465;
t466 = t503 * t697 + t508 * t699;
t453 = pkin(10) * t726 + t466;
t784 = qJD(5) * t706;
t746 = -t451 * t702 - t453 * t706 - t482 * t784 + t484 * t785;
t858 = t532 * t536 + t746;
t802 = t697 * t864 + t699 * t865;
t834 = t697 * t865 - t699 * t864;
t827 = t523 * t840;
t851 = pkin(4) * t769 + pkin(10) * t856 - t802;
t850 = -pkin(10) * t849 + t834;
t446 = t482 * t702 + t484 * t706;
t712 = -qJD(5) * t446 + t451 * t706 - t453 * t702;
t433 = -pkin(5) * t743 - t712;
t431 = t433 * t701;
t444 = -pkin(11) * t670 + t446;
t471 = pkin(5) * t532 - pkin(11) * t840 + t536;
t437 = t444 * t705 + t471 * t701;
t848 = t437 * t840 + t443 * t782 + t431;
t847 = -t536 * t840 + t712;
t735 = t444 * t701 - t471 * t705;
t845 = -t433 * t705 + t443 * t783 + t735 * t840;
t724 = t661 * t706 - t662 * t702;
t804 = qJD(5) * t724 - t702 * t849 + t706 * t856;
t601 = t661 * t702 + t662 * t706;
t803 = qJD(5) * t601 + t702 * t856 + t706 * t849;
t841 = -t649 + (-t744 + t787) * pkin(3);
t694 = t698 ^ 2;
t839 = -0.2e1 * t694 * t777;
t838 = (t704 ^ 2 - t708 ^ 2) * MDP(5);
t680 = t828 * t703;
t681 = t828 * t707;
t609 = t680 * t699 + t681 * t697;
t594 = -pkin(10) * t662 + t609;
t610 = t680 * t697 - t681 * t699;
t595 = pkin(10) * t661 + t610;
t728 = t594 * t706 - t595 * t702;
t837 = qJD(5) * t728 - t702 * t851 + t706 * t850;
t544 = t594 * t702 + t595 * t706;
t836 = qJD(5) * t544 + t702 * t850 + t706 * t851;
t658 = t700 * t703 + t707 * t814;
t813 = t698 * t708;
t831 = pkin(1) * t704;
t641 = pkin(8) * t813 + (pkin(9) + t831) * t700;
t752 = -t641 * t703 + t642 * t707;
t565 = -pkin(3) * t813 - qJ(4) * t658 + t752;
t657 = -t700 * t707 + t703 * t814;
t800 = t641 * t707 + t642 * t703;
t572 = -qJ(4) * t657 + t800;
t517 = t565 * t699 - t572 * t697;
t597 = -t657 * t697 + t658 * t699;
t495 = -pkin(4) * t813 - pkin(10) * t597 + t517;
t518 = t565 * t697 + t572 * t699;
t596 = -t657 * t699 - t658 * t697;
t497 = pkin(10) * t596 + t518;
t835 = t495 * t702 + t497 * t706;
t801 = pkin(4) * t849 + t841;
t692 = pkin(3) * t699 + pkin(4);
t829 = pkin(3) * t697;
t796 = t692 * t702 + t706 * t829;
t477 = t705 * t479;
t833 = -t783 * t859 - t477;
t766 = t698 * t788;
t604 = -qJD(3) * t657 + t707 * t766;
t713 = -qJD(3) * t800 + t648 * t707 - t650 * t703;
t790 = qJD(2) * t704;
t767 = t698 * t790;
t515 = pkin(3) * t767 - qJ(4) * t604 - qJD(4) * t658 + t713;
t603 = qJD(3) * t658 + t698 * t765;
t719 = -t641 * t787 + t642 * t786 + t648 * t703 + t650 * t707;
t519 = -qJ(4) * t603 - qJD(4) * t657 + t719;
t473 = t515 * t699 - t519 * t697;
t571 = -t603 * t697 + t604 * t699;
t462 = pkin(4) * t767 - pkin(10) * t571 + t473;
t474 = t515 * t697 + t519 * t699;
t570 = -t603 * t699 - t604 * t697;
t464 = pkin(10) * t570 + t474;
t832 = -qJD(5) * t835 + t462 * t706 - t464 * t702;
t820 = t601 * t705;
t819 = t622 * t868;
t818 = t624 * t868;
t817 = t868 * t703;
t816 = t868 * t707;
t709 = qJD(1) ^ 2;
t815 = t694 * t709;
t807 = pkin(5) * t769 + t836;
t507 = t554 * t699 - t547;
t506 = -t554 * t697 - t812;
t488 = t506 + t843;
t489 = t507 + t844;
t722 = t692 * t706 - t702 * t829;
t798 = qJD(5) * t722 - t488 * t702 - t489 * t706;
t797 = qJD(5) * t796 + t488 * t706 - t489 * t702;
t634 = pkin(8) * t742 + qJD(2) * t685;
t651 = pkin(1) * t700 * t790 + pkin(8) * t766;
t789 = qJD(2) * t707;
t780 = qJD(2) - t688;
t775 = t701 * t813;
t771 = -pkin(3) * t707 - pkin(2);
t770 = t694 * t791;
t432 = pkin(11) * t743 - t746;
t578 = pkin(3) * t599 + t634;
t520 = -pkin(4) * t726 + t578;
t439 = pkin(5) * t479 - pkin(11) * t478 + t520;
t759 = -t432 * t701 + t439 * t705;
t754 = -t701 * t804 - t705 * t769;
t753 = t701 * t769 - t705 * t804;
t558 = pkin(3) * t624 - pkin(4) * t588;
t645 = pkin(11) + t796;
t748 = qJD(6) * t645 + t490 + t558;
t747 = t694 * t704 * t708 * MDP(4);
t741 = pkin(3) * t603 + t651;
t740 = pkin(1) * t839;
t631 = -pkin(4) * t661 + t771;
t542 = -pkin(5) * t724 - pkin(11) * t601 + t631;
t739 = pkin(11) * t769 - qJD(6) * t542 - t837;
t738 = -pkin(5) * t803 + pkin(11) * t804 + qJD(6) * t544 - t801;
t737 = t432 * t705 + t439 * t701;
t736 = -t479 * t645 + t861;
t460 = -pkin(11) * t813 + t835;
t546 = t596 * t702 + t597 * t706;
t640 = t689 + (-pkin(2) - t830) * t700;
t714 = pkin(3) * t657 + t640;
t557 = -pkin(4) * t596 + t714;
t727 = t596 * t706 - t597 * t702;
t480 = -pkin(5) * t727 - pkin(11) * t546 + t557;
t734 = t460 * t705 + t480 * t701;
t733 = -t460 * t701 + t480 * t705;
t731 = t495 * t706 - t497 * t702;
t723 = t532 * t750 - t833;
t537 = t546 * t701 + t705 * t813;
t528 = -pkin(4) * t570 + t741;
t718 = t462 * t702 + t464 * t706 + t495 * t784 - t497 * t785;
t717 = t601 * t782 - t754;
t716 = -t601 * t783 - t753;
t644 = -pkin(5) - t722;
t538 = t546 * t705 - t775;
t486 = qJD(5) * t546 - t570 * t706 + t571 * t702;
t485 = qJD(5) * t727 + t570 * t702 + t571 * t706;
t468 = -qJD(6) * t775 + t485 * t701 + t546 * t782 - t705 * t767;
t467 = -qJD(6) * t537 + t485 * t705 + t701 * t767;
t459 = pkin(5) * t813 - t731;
t440 = pkin(5) * t486 - pkin(11) * t485 + t528;
t435 = -pkin(5) * t767 - t832;
t434 = pkin(11) * t767 + t718;
t430 = -qJD(6) * t437 + t759;
t429 = -qJD(6) * t735 + t737;
t1 = [(t478 * t727 - t479 * t546 - t485 * t532 - t486 * t840) * MDP(21) + (t718 * t670 + t528 * t840 + t557 * t478 + t520 * t546 + t536 * t485 + (-t746 * t708 + (-qJD(1) * t835 - t446) * t790) * t698) * MDP(26) + (-t485 * t670 + (-t478 * t708 + (qJD(1) * t546 + t840) * t790) * t698) * MDP(22) + (t478 * t546 + t485 * t840) * MDP(20) + (-t465 * t597 + t466 * t596 + t473 * t588 - t474 * t725 - t501 * t571 + t502 * t570 - t517 * t556 + t518 * t726) * MDP(18) + (t698 * t868 - t770) * MDP(15) * t790 + (-t719 * t868 + t651 * t624 + t640 * t598 + t634 * t658 + t611 * t604 + (t720 * t708 + (-qJD(1) * t800 - t582) * t790) * t698) * MDP(17) + (t604 * t868 + (-t598 * t708 + (qJD(1) * t658 + t624) * t790) * t698) * MDP(13) + (-t603 * t868 + (t599 * t708 + (-qJD(1) * t657 - t622) * t790) * t698) * MDP(14) + (t713 * t868 + t651 * t622 + t640 * t599 + t634 * t657 + t611 * t603 + (-t711 * t708 + (qJD(1) * t752 + t581) * t790) * t698) * MDP(16) + (MDP(6) * t766 - MDP(7) * t767) * (t688 + t792) + (-t457 * t537 - t458 * t538 - t467 * t523 - t468 * t525) * MDP(28) + (t457 * t538 + t467 * t525) * MDP(27) + (-t598 * t657 - t599 * t658 - t603 * t624 - t604 * t622) * MDP(12) + (t598 * t658 + t604 * t624) * MDP(11) + (-t832 * t670 + t528 * t532 + t557 * t479 - t520 * t727 + t536 * t486 + (-t712 * t708 + (qJD(1) * t731 + t445) * t790) * t698) * MDP(25) + (t486 * t670 + (t479 * t708 + (qJD(1) * t727 - t532) * t790) * t698) * MDP(23) + ((qJD(6) * t733 + t434 * t705 + t440 * t701) * t859 - t734 * t479 + t429 * t727 - t437 * t486 + t435 * t525 + t459 * t457 + t433 * t538 + t443 * t467) * MDP(33) + (t458 * t727 + t468 * t859 - t479 * t537 - t486 * t523) * MDP(30) + (-t457 * t727 - t467 * t859 + t479 * t538 + t486 * t525) * MDP(29) + (-t479 * t727 - t486 * t859) * MDP(31) + (-(-qJD(6) * t734 - t434 * t701 + t440 * t705) * t859 + t733 * t479 - t430 * t727 - t735 * t486 + t435 * t523 + t459 * t458 + t433 * t537 + t443 * t468) * MDP(32) + 0.2e1 * t747 * t777 + (-t670 * t698 - t770) * MDP(24) * t790 + (-t633 * t700 - t650 * t688 + t708 * t740) * MDP(10) + (-t634 * t700 - t651 * t688 + t704 * t740) * MDP(9) + (t465 * t517 + t466 * t518 + t501 * t473 + t502 * t474 + t578 * t714 + t587 * t741) * MDP(19) + t838 * t839; (t465 * t609 + t466 * t610 + t802 * t501 + t834 * t502 + t578 * t771 + t587 * t841) * MDP(19) + ((t598 - t819) * t707 + (-t599 - t818) * t703) * MDP(12) + (-t465 * t662 + t466 * t661 - t501 * t856 - t849 * t502 - t609 * t556 + t802 * t588 + t610 * t726 - t834 * t725) * MDP(18) + (t598 * t703 + t624 * t816) * MDP(11) + (t478 * t724 - t479 * t601 - t532 * t804 - t803 * t840) * MDP(21) + (t631 * t478 + t520 * t601 + t804 * t536 + t801 * t840) * MDP(26) + (t478 * t601 + t804 * t840) * MDP(20) + (-MDP(22) * t804 + MDP(23) * t803 + MDP(25) * t836 + MDP(26) * t837) * t670 + (t631 * t479 - t520 * t724 + t801 * t532 + t803 * t536) * MDP(25) + (-pkin(2) * t598 + t634 * t703 + t799 * t868 - t649 * t624 + (pkin(9) * t817 + t611 * t707) * qJD(3) + (-t611 * t811 + (-pkin(9) * t789 + t582) * t704) * t793) * MDP(17) + ((qJD(2) * t724 + t532) * MDP(23) + (qJD(2) * t728 - t445) * MDP(25) + (-qJD(2) * t544 + t446) * MDP(26) - t868 * MDP(15) + t670 * MDP(24) + (qJD(2) * t601 - t840) * MDP(22) - t780 * MDP(7)) * t769 + (t868 * t786 + (-t868 * t811 + (qJD(2) * t703 - t624) * t704) * t793) * MDP(13) + (-pkin(2) * t599 - t634 * t707 - t751 * t868 - t649 * t622 + (-pkin(9) * t816 + t611 * t703) * qJD(3) + (-t581 * t704 + (-pkin(9) * t790 - t611 * t708) * t703) * t793) * MDP(16) + (-t868 * t787 + (t708 * t817 + (t622 + t789) * t704) * t793) * MDP(14) + (pkin(8) * t743 + t646 * t688 + (-t700 * t777 + t815) * t830) * MDP(10) + (t649 * t688 + t815 * t831 - t634) * MDP(9) + (t754 * t525 + t753 * t523 + (-t455 - t458 * t705 + (t523 * t701 - t525 * t705) * qJD(6)) * t601) * MDP(28) + (t457 * t820 + t525 * t716) * MDP(27) + (t458 * t724 - t475 * t601 - t523 * t803 + t717 * t859) * MDP(30) + (-t457 * t724 + t477 * t601 + t525 * t803 - t716 * t859) * MDP(29) + (-t479 * t724 - t803 * t859) * MDP(31) + (-(t542 * t701 + t544 * t705) * t479 + t429 * t724 - t728 * t457 + t433 * t820 - (t701 * t738 + t705 * t739) * t859 + t807 * t525 - t803 * t437 + t716 * t443) * MDP(33) + ((t542 * t705 - t544 * t701) * t479 - t430 * t724 - t728 * t458 + t601 * t431 - (t701 * t739 - t705 * t738) * t859 + t807 * t523 - t803 * t735 + t717 * t443) * MDP(32) - t709 * t747 + t815 * t838 + t780 * MDP(6) * t768; (t644 * t458 + t736 * t701 + t797 * t523 - (-t701 * t798 - t705 * t748) * t859 + t845) * MDP(32) + (-t532 * t558 + t670 * t797 + t722 * t743 + t847) * MDP(25) + (t644 * t457 + t736 * t705 + t797 * t525 - (t701 * t748 - t705 * t798) * t859 + t848) * MDP(33) + (-t599 + t818) * MDP(14) + (-t501 * t506 - t502 * t507 + (t465 * t699 + t466 * t697 - t587 * t624) * pkin(3)) * MDP(19) + (-t558 * t840 + t670 * t798 - t743 * t796 + t858) * MDP(26) + ((-t699 * t556 + t697 * t726) * pkin(3) + (t507 - t501) * t725 + (-t502 - t506) * t588) * MDP(18) + (-t622 ^ 2 + t624 ^ 2) * MDP(12) + (t723 + t827) * MDP(30) + (t598 + t819) * MDP(13) + (t582 * t868 - t611 * t624 + t711) * MDP(16) + t624 * t622 * MDP(11) + MDP(15) * t743 + (t581 * t868 + t611 * t622 - t720) * MDP(17) + t867; (-t588 ^ 2 - t725 ^ 2) * MDP(18) + (-t501 * t588 + t502 * t725 + t578) * MDP(19) + (t479 - t825) * MDP(25) + (t478 + t823) * MDP(26) + (t723 - t827) * MDP(32) + (-t705 * t859 ^ 2 - t475 - t826) * MDP(33); (-t446 * t670 + t847) * MDP(25) + (-t445 * t670 + t858) * MDP(26) + (-t750 * t859 + t477 + t827) * MDP(30) + (-pkin(5) * t458 + (-t445 * t701 + t490 * t705) * t859 - t446 * t523 + t701 * t861 - t808 * pkin(11) + t845) * MDP(32) + (-pkin(5) * t457 - (t445 * t705 + t490 * t701) * t859 - t446 * t525 + t443 * t860 + t833 * pkin(11) + t848) * MDP(33) + t867; t525 * t523 * MDP(27) + (-t523 ^ 2 + t525 ^ 2) * MDP(28) + (-t523 * t859 + t773) * MDP(29) + (-t525 * t859 - t756) * MDP(30) + t479 * MDP(31) + (-t437 * t859 - t443 * t525 + t759) * MDP(32) + (t443 * t523 + t735 * t859 - t737) * MDP(33) + (-MDP(29) * t821 - MDP(30) * t525 - MDP(32) * t437 + MDP(33) * t735) * qJD(6);];
tauc  = t1;
