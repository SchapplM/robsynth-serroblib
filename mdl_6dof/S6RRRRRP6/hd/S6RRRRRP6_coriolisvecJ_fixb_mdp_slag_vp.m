% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRP6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRRP6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:31:50
% EndTime: 2019-03-10 01:32:09
% DurationCPUTime: 13.69s
% Computational Cost: add. (13981->642), mult. (33524->811), div. (0->0), fcn. (23894->8), ass. (0->264)
t659 = cos(qJ(3));
t657 = sin(qJ(2));
t660 = cos(qJ(2));
t781 = t659 * t660;
t696 = pkin(3) * t657 - pkin(9) * t781;
t798 = pkin(8) + pkin(9);
t726 = qJD(3) * t798;
t700 = pkin(2) * t657 - pkin(8) * t660;
t609 = t700 * qJD(1);
t656 = sin(qJ(3));
t750 = qJD(1) * t657;
t723 = t656 * t750;
t757 = pkin(7) * t723 + t659 * t609;
t841 = qJD(1) * t696 + t659 * t726 + t757;
t588 = t656 * t609;
t783 = t657 * t659;
t840 = -t588 - (-pkin(9) * t656 * t660 - pkin(7) * t783) * qJD(1) - t656 * t726;
t654 = sin(qJ(5));
t745 = qJD(2) * t656;
t603 = t659 * t750 + t745;
t740 = qJD(3) * t659;
t720 = t657 * t740;
t743 = qJD(2) * t660;
t722 = t656 * t743;
t676 = t720 + t722;
t731 = qJD(2) * qJD(3);
t563 = qJD(1) * t676 + t656 * t731;
t734 = t659 * qJD(2);
t601 = t723 - t734;
t655 = sin(qJ(4));
t658 = cos(qJ(4));
t741 = qJD(3) * t657;
t721 = t656 * t741;
t677 = t660 * t734 - t721;
t668 = qJD(1) * t677 + t659 * t731;
t738 = qJD(4) * t658;
t728 = -t655 * t563 - t601 * t738 + t658 * t668;
t739 = qJD(4) * t655;
t683 = t603 * t739 - t728;
t695 = t658 * t563 - t601 * t739 + t603 * t738 + t655 * t668;
t698 = -t601 * t655 + t658 * t603;
t699 = -t601 * t658 - t655 * t603;
t797 = cos(qJ(5));
t719 = qJD(5) * t797;
t737 = qJD(5) * t654;
t442 = t654 * t695 + t797 * t683 + t698 * t737 - t699 * t719;
t495 = t654 * t698 - t699 * t797;
t749 = qJD(1) * t660;
t638 = -qJD(3) + t749;
t625 = -qJD(4) + t638;
t617 = -qJD(5) + t625;
t426 = -t495 * t617 - t442;
t802 = t654 * t699 + t698 * t797;
t443 = qJD(5) * t802 - t654 * t683 + t797 * t695;
t732 = qJD(1) * qJD(2);
t718 = t657 * t732;
t799 = t802 ^ 2;
t837 = t495 * t802;
t839 = MDP(29) * t718 + (-t617 * t802 - t443) * MDP(28) + MDP(25) * t837 + (-t495 ^ 2 + t799) * MDP(26) + t426 * MDP(27);
t604 = t655 * t656 - t658 * t659;
t684 = t604 * t660;
t808 = qJD(3) + qJD(4);
t824 = t808 * t604;
t763 = qJD(1) * t684 - t824;
t622 = t798 * t656;
t623 = t798 * t659;
t758 = -t655 * t622 + t658 * t623;
t838 = qJD(4) * t758 + t840 * t655 + t658 * t841;
t614 = t617 * qJD(6);
t635 = qJ(6) * t718;
t616 = -pkin(2) * t660 - pkin(8) * t657 - pkin(1);
t594 = t616 * qJD(1);
t647 = pkin(7) * t749;
t621 = qJD(2) * pkin(8) + t647;
t554 = t659 * t594 - t621 * t656;
t521 = -pkin(9) * t603 + t554;
t513 = -pkin(3) * t638 + t521;
t785 = t656 * t594;
t555 = t659 * t621 + t785;
t522 = -pkin(9) * t601 + t555;
t518 = t658 * t522;
t678 = t696 * qJD(2);
t612 = t700 * qJD(2);
t595 = qJD(1) * t612;
t705 = pkin(7) * t718;
t761 = t659 * t595 + t656 * t705;
t477 = qJD(1) * t678 - qJD(3) * t522 + t761;
t742 = qJD(3) * t656;
t762 = t594 * t740 + t656 * t595;
t669 = -t621 * t742 - t659 * t705 + t762;
t488 = -pkin(9) * t563 + t669;
t713 = -t658 * t477 + t655 * t488;
t422 = -t728 * pkin(10) + pkin(4) * t718 + (-t518 + (pkin(10) * t603 - t513) * t655) * qJD(4) - t713;
t702 = -t655 * t477 - t658 * t488 - t513 * t738 + t522 * t739;
t425 = -pkin(10) * t695 - t702;
t516 = t655 * t522;
t467 = t658 * t513 - t516;
t817 = pkin(10) * t698;
t459 = t467 - t817;
t455 = -pkin(4) * t625 + t459;
t468 = t655 * t513 + t518;
t816 = pkin(10) * t699;
t460 = t468 + t816;
t704 = -t654 * t422 - t797 * t425 - t455 * t719 + t460 * t737;
t414 = t635 - t614 - t704;
t791 = qJD(2) * pkin(2);
t620 = pkin(7) * t750 - t791;
t565 = pkin(3) * t601 + t620;
t511 = -pkin(4) * t699 + t565;
t446 = pkin(5) * t495 - qJ(6) * t802 + t511;
t827 = t446 * t495;
t836 = t414 - t827;
t605 = t655 * t659 + t656 * t658;
t569 = t605 * t749;
t680 = t605 * qJD(3);
t667 = -qJD(4) * t605 - t680;
t829 = t569 + t667;
t826 = t495 * t511;
t834 = t704 + t826;
t833 = -t622 * t738 - t655 * t841 + t840 * t658;
t832 = MDP(22) * t718 + (t698 ^ 2 - t699 ^ 2) * MDP(19) + (t625 * t699 - t683) * MDP(20) + (-t625 * t698 - t695) * MDP(21) - t699 * MDP(18) * t698 + t839;
t831 = pkin(4) * t750 + pkin(10) * t763 + t838;
t786 = t655 * t623;
t792 = pkin(10) * t605;
t830 = (-t786 - t792) * qJD(4) + t833 + (t569 - t680) * pkin(10);
t456 = pkin(5) * t802 + qJ(6) * t495;
t813 = t446 * t802;
t825 = t511 * t802;
t771 = t604 * t719 + t605 * t737 - t654 * t829 - t763 * t797;
t551 = -t654 * t604 + t605 * t797;
t770 = qJD(5) * t551 + t654 * t763 - t797 * t829;
t796 = pkin(3) * t656;
t597 = t749 * t796 + t647;
t823 = pkin(3) * t742 - t597;
t820 = -0.2e1 * t732;
t819 = pkin(4) * t698;
t815 = MDP(4) * t657;
t652 = t657 ^ 2;
t814 = MDP(5) * (-t660 ^ 2 + t652);
t710 = -t658 * t622 - t786;
t527 = t710 - t792;
t528 = -pkin(10) * t604 + t758;
t690 = t527 * t797 - t654 * t528;
t812 = qJD(5) * t690 - t654 * t831 + t797 * t830;
t490 = t654 * t527 + t528 * t797;
t811 = qJD(5) * t490 + t654 * t830 + t797 * t831;
t712 = -t521 * t655 - t518;
t462 = t712 - t816;
t768 = t658 * t521 - t516;
t463 = t768 - t817;
t644 = pkin(3) * t658 + pkin(4);
t730 = t654 * t655 * pkin(3);
t810 = -t654 * t462 + t644 * t719 + (-qJD(4) - qJD(5)) * t730 + (pkin(3) * t738 - t463) * t797;
t578 = t604 * t657;
t600 = t659 * t616;
t795 = pkin(7) * t656;
t553 = -pkin(9) * t783 + t600 + (-pkin(3) - t795) * t660;
t641 = pkin(7) * t781;
t754 = t656 * t616 + t641;
t784 = t656 * t657;
t559 = -pkin(9) * t784 + t754;
t711 = t658 * t553 - t559 * t655;
t485 = -pkin(4) * t660 + pkin(10) * t578 + t711;
t577 = t605 * t657;
t764 = t655 * t553 + t658 * t559;
t491 = -pkin(10) * t577 + t764;
t809 = t654 * t485 + t797 * t491;
t767 = -pkin(4) * t829 + t823;
t746 = qJD(2) * t605;
t663 = t657 * t824 - t660 * t746;
t736 = t652 * qJD(1);
t703 = -t797 * t422 + t654 * t425 + t455 * t737 + t460 * t719;
t807 = -t703 - t825;
t670 = -qJD(4) * t468 - t713;
t806 = -t565 * t698 + t670;
t805 = -t565 * t699 + t702;
t509 = -qJD(2) * t684 + t657 * t667;
t744 = qJD(2) * t657;
t759 = t659 * t612 + t744 * t795;
t503 = t678 + (-t641 + (pkin(9) * t657 - t616) * t656) * qJD(3) + t759;
t760 = t656 * t612 + t616 * t740;
t508 = -t676 * pkin(9) + (-t657 * t734 - t660 * t742) * pkin(7) + t760;
t671 = -qJD(4) * t764 + t658 * t503 - t508 * t655;
t437 = pkin(4) * t744 - pkin(10) * t509 + t671;
t729 = t655 * t503 + t658 * t508 + t553 * t738;
t441 = (-t783 * t808 - t722) * pkin(10) * t658 + (-qJD(4) * t559 + (qJD(4) * t784 - t677) * pkin(10)) * t655 + t729;
t803 = -qJD(5) * t809 + t437 * t797 - t654 * t441;
t794 = pkin(8) * t638;
t725 = t797 * t460;
t430 = t654 * t455 + t725;
t790 = t430 * t617;
t789 = t603 * t638;
t788 = t620 * t656;
t787 = t654 * t460;
t661 = qJD(2) ^ 2;
t782 = t657 * t661;
t780 = t660 * t638;
t779 = t660 * t661;
t662 = qJD(1) ^ 2;
t778 = t660 * t662;
t432 = t459 * t797 - t787;
t777 = -pkin(4) * t719 - qJD(6) + t432;
t776 = pkin(5) * t770 + qJ(6) * t771 - t551 * qJD(6) + t767;
t775 = -qJ(6) * t750 + t812;
t774 = pkin(5) * t750 + t811;
t773 = -qJD(6) - t810;
t724 = t797 * t655;
t766 = t462 * t797 - t654 * t463 + t644 * t737 + (t655 * t719 + (t654 * t658 + t724) * qJD(4)) * pkin(3);
t755 = pkin(3) * t724 + t654 * t644;
t613 = pkin(3) * t784 + t657 * pkin(7);
t748 = qJD(2) * t690;
t747 = qJD(2) * t490;
t735 = t657 * MDP(15);
t429 = t455 * t797 - t787;
t733 = qJD(6) - t429;
t566 = pkin(3) * t676 + pkin(7) * t743;
t645 = -pkin(3) * t659 - pkin(2);
t717 = t660 * t732;
t715 = qJD(2) * t735;
t542 = t563 * pkin(3) + pkin(7) * t717;
t714 = pkin(1) * t820;
t709 = t638 + t749;
t708 = t601 + t734;
t707 = -t603 + t745;
t706 = qJD(3) + t749;
t556 = pkin(4) * t577 + t613;
t431 = t654 * t459 + t725;
t701 = pkin(4) * t737 - t431;
t515 = pkin(3) * t603 + t819;
t637 = pkin(5) * t718;
t415 = -t637 + t703;
t571 = pkin(4) * t604 + t645;
t693 = t485 * t797 - t654 * t491;
t524 = -t654 * t577 - t578 * t797;
t688 = t706 * t745;
t687 = -t429 * t617 + t704;
t686 = t703 + t813;
t685 = t565 * t605;
t682 = t654 * t437 + t797 * t441 + t485 * t719 - t491 * t737;
t679 = t644 * t797 - t730;
t673 = t617 * t766 - t703;
t464 = pkin(4) * t695 + t542;
t492 = -pkin(4) * t663 + t566;
t643 = -pkin(4) * t797 - pkin(5);
t642 = pkin(4) * t654 + qJ(6);
t579 = -pkin(5) - t679;
t576 = qJ(6) + t755;
t550 = t604 * t797 + t605 * t654;
t523 = t577 * t797 - t578 * t654;
t487 = pkin(5) * t550 - qJ(6) * t551 + t571;
t465 = pkin(5) * t523 - qJ(6) * t524 + t556;
t452 = qJD(5) * t524 + t654 * t509 - t663 * t797;
t451 = -t509 * t797 + t577 * t719 - t578 * t737 - t654 * t663;
t450 = t456 + t819;
t449 = t660 * pkin(5) - t693;
t448 = -qJ(6) * t660 + t809;
t447 = t456 + t515;
t428 = -t617 * qJ(6) + t430;
t427 = t617 * pkin(5) + t733;
t419 = t452 * pkin(5) + t451 * qJ(6) - t524 * qJD(6) + t492;
t418 = t443 * pkin(5) + t442 * qJ(6) - qJD(6) * t802 + t464;
t417 = -pkin(5) * t744 - t803;
t416 = qJ(6) * t744 - qJD(6) * t660 + t682;
t1 = [(-(-t616 * t742 + t759) * t638 + (t620 * t740 + pkin(7) * t563 + (t600 * qJD(1) + t554) * qJD(2)) * t657 + ((pkin(7) * t601 + t788) * qJD(2) + (t785 + (pkin(7) * t638 + t621) * t659) * qJD(3) - t761) * t660) * MDP(16) + ((-t601 * t659 - t603 * t656) * t743 + ((t601 + t723) * t742 + (-t603 * qJD(3) - t563 - t688) * t659) * t657) * MDP(12) + (-t671 * t625 - t566 * t699 + t613 * t695 + t542 * t577 + (qJD(2) * t685 - t670) * t660 + (-t565 * t824 + (qJD(1) * t711 + t467) * qJD(2)) * t657) * MDP(23) + ((t625 * t746 + t695) * t660 + (-t824 * t625 + (-t577 * qJD(1) + t699) * qJD(2)) * t657) * MDP(21) - MDP(7) * t782 + (pkin(7) * t782 + t660 * t714) * MDP(10) + (t603 * t677 + t668 * t783) * MDP(11) + (-pkin(7) * t779 + t657 * t714) * MDP(9) + (t638 * t720 + t563 * t660 + (-t601 * t657 + (-t736 + t780) * t656) * qJD(2)) * MDP(14) + (-t414 * t660 - t416 * t617 - t418 * t524 - t419 * t802 + t442 * t465 + t446 * t451) * MDP(34) + (-t556 * t442 - t511 * t451 + t464 * t524 + t492 * t802 + t682 * t617 - t704 * t660) * MDP(31) + (-t442 * t524 - t451 * t802) * MDP(25) + (t442 * t523 - t443 * t524 + t451 * t495 - t452 * t802) * MDP(26) + (-t414 * t523 + t415 * t524 - t416 * t495 + t417 * t802 - t427 * t451 - t428 * t452 - t442 * t449 - t443 * t448) * MDP(33) + ((-qJD(1) * t578 + t698) * MDP(20) + (qJD(1) * t524 + t802) * MDP(27) + (-qJD(1) * t523 - t495) * MDP(28) + (qJD(1) * t448 + t428) * MDP(34) + (-qJD(1) * t449 - t427) * MDP(32) + t429 * MDP(30) + (-qJD(1) * t809 - t430) * MDP(31) + (-qJD(1) * t764 - t468) * MDP(24) + (-t625 - t749) * MDP(22) + (-t617 - t749) * MDP(29)) * t744 + (t415 * t660 + t417 * t617 + t418 * t523 + t419 * t495 + t443 * t465 + t446 * t452) * MDP(32) + (t443 * t660 + t452 * t617) * MDP(28) + (t442 * t660 + t451 * t617) * MDP(27) + (-t509 * t625 + t683 * t660) * MDP(20) + (t509 * t698 + t578 * t683) * MDP(18) + ((-t559 * t739 + t729) * t625 - t702 * t660 + t566 * t698 - t613 * t683 - t542 * t578 + t565 * t509) * MDP(24) + (t709 * t721 + (t603 * t657 + (t736 + (-t638 - t706) * t660) * t659) * qJD(2)) * MDP(13) + t814 * t820 + MDP(6) * t779 - t709 * t715 + 0.2e1 * t717 * t815 + (t509 * t699 + t577 * t683 + t578 * t695 + t663 * t698) * MDP(19) + (t760 * t638 + t762 * t660 + (-t620 * t657 - t621 * t660 + (-t780 - t736) * pkin(7)) * t742 + ((pkin(7) * t603 + t620 * t659) * t660 + (-t754 * qJD(1) - t555 + (-t638 + t706) * pkin(7) * t659) * t657) * qJD(2)) * MDP(17) + (t556 * t443 + t511 * t452 + t464 * t523 + t492 * t495 - t617 * t803 + t703 * t660 + t693 * t718) * MDP(30) + (t414 * t448 + t415 * t449 + t416 * t428 + t417 * t427 + t418 * t465 + t419 * t446) * MDP(35); (-t656 ^ 2 * qJD(1) * t741 + (t688 - t789) * t659) * MDP(11) + ((-t563 + t789) * t656 + ((-t601 + t734) * qJD(3) + (t660 * t708 - t721) * qJD(1)) * t659) * MDP(12) + (t638 * t742 + (-t656 * t780 + t657 * t708) * qJD(1)) * MDP(14) + (-t638 * t740 + (t657 * t707 + t659 * t780) * qJD(1)) * MDP(13) + (t418 * t550 + t443 * t487 + t770 * t446 + t776 * t495) * MDP(32) + (t571 * t443 + t464 * t550 + t767 * t495 + t770 * t511) * MDP(30) + (t414 * t490 - t415 * t690 + t418 * t487 + t427 * t774 + t428 * t775 + t446 * t776) * MDP(35) + (-t414 * t550 + t415 * t551 - t427 * t771 - t428 * t770 + t442 * t690 - t443 * t490 - t495 * t775 + t774 * t802) * MDP(33) + (t625 * MDP(22) + t617 * MDP(29) + (-qJD(2) * t604 - t699) * MDP(21) + (t430 - t747) * MDP(31) + (t427 + t748) * MDP(32) + (-t428 + t747) * MDP(34) + (-qJD(2) * t758 + t468) * MDP(24) + (-qJD(2) * t550 + t495) * MDP(28) + (-t429 + t748) * MDP(30) + (qJD(2) * t551 - t802) * MDP(27) + (-t698 + t746) * MDP(20) + (qJD(2) * t710 - t467) * MDP(23)) * t750 + (-t571 * t442 + t464 * t551 - t771 * t511 + t767 * t802) * MDP(31) + (-t418 * t551 + t442 * t487 + t771 * t446 - t776 * t802) * MDP(34) + (t442 * t550 - t443 * t551 + t495 * t771 - t770 * t802) * MDP(26) + (-t442 * t551 - t771 * t802) * MDP(25) + t638 * qJD(1) * t735 + (-t588 * t638 + (-t656 * t794 + (t620 - t791) * t659) * qJD(3) + ((-t620 - t791) * t781 + (pkin(2) * t742 - pkin(8) * t734 + t555) * t657 + (t638 * t783 + t660 * t707) * pkin(7)) * qJD(1)) * MDP(17) + (-pkin(2) * t563 + t757 * t638 + (t659 * t794 + t788) * qJD(3) + ((-pkin(8) * t745 - t554) * t657 + (-pkin(7) * t708 - t788) * t660) * qJD(1)) * MDP(16) + (t683 * t604 - t605 * t695 + t698 * t829 + t763 * t699) * MDP(19) + (t645 * t695 + t542 * t604 + t597 * t699 - t565 * t569 + t685 * qJD(4) + (-t699 * t796 + t685) * qJD(3)) * MDP(23) + (-t605 * t683 + t698 * t763) * MDP(18) + t662 * t814 + ((-t623 * t739 + t833) * MDP(24) + t838 * MDP(23) + (t605 * t808 - t569) * MDP(21) - t763 * MDP(20)) * t625 + (MDP(9) * t657 * t662 + MDP(10) * t778) * pkin(1) - t778 * t815 + (t771 * MDP(27) + t770 * MDP(28) + MDP(30) * t811 + MDP(31) * t812 + t774 * MDP(32) - t775 * MDP(34)) * t617 + (t542 * t605 + t763 * t565 - t645 * t683 + t823 * t698) * MDP(24); (-t563 - t789) * MDP(14) + (-t603 * t620 + t761 + (-qJD(3) - t638) * t555) * MDP(16) + (t414 * t576 + t415 * t579 + t427 * t766 - t428 * t773 - t446 * t447) * MDP(35) + (-t447 * t495 - t579 * t718 + t637 + t673 - t813) * MDP(32) + (-t515 * t495 + t679 * t718 + t673 - t825) * MDP(30) + (-t554 * t638 + t601 * t620 - t669) * MDP(17) + (-t601 * t638 + t668) * MDP(13) + (-t515 * t802 + t617 * t810 - t755 * t718 + t834) * MDP(31) + (t447 * t802 + t576 * t718 + t617 * t773 + t836) * MDP(34) + (-t442 * t579 - t443 * t576 + (t428 + t766) * t802 + (t427 + t773) * t495) * MDP(33) + (t712 * t625 + (t603 * t699 + t625 * t739 + t658 * t718) * pkin(3) + t806) * MDP(23) + (-t768 * t625 + (-t603 * t698 + t625 * t738 - t655 * t718) * pkin(3) + t805) * MDP(24) + t603 * t601 * MDP(11) + qJD(1) * t715 + (-t601 ^ 2 + t603 ^ 2) * MDP(12) + t832; (-t468 * t625 + t806) * MDP(23) + (-t467 * t625 + t805) * MDP(24) + (-t431 * t617 + (-t495 * t698 + t617 * t737 + t718 * t797) * pkin(4) + t807) * MDP(30) + (-t432 * t617 + (t617 * t719 - t654 * t718 - t698 * t802) * pkin(4) + t834) * MDP(31) + (-t450 * t495 + t617 * t701 - t643 * t718 - t415 - t813) * MDP(32) + (-t442 * t643 - t443 * t642 + (t428 + t701) * t802 + (t427 + t777) * t495) * MDP(33) + (t450 * t802 + t617 * t777 + t642 * t718 + t836) * MDP(34) + (t414 * t642 + t415 * t643 + t427 * t701 - t428 * t777 - t446 * t450) * MDP(35) + t832; (-t790 + t807) * MDP(30) + (t687 + t826) * MDP(31) + (-t456 * t495 + 0.2e1 * t637 - t686 - t790) * MDP(32) + (pkin(5) * t442 - qJ(6) * t443 + (t428 - t430) * t802 + (t427 - t733) * t495) * MDP(33) + (t456 * t802 - 0.2e1 * t614 + 0.2e1 * t635 - t687 - t827) * MDP(34) + (-pkin(5) * t415 + qJ(6) * t414 - t427 * t430 + t428 * t733 - t446 * t456) * MDP(35) + t839; (-t718 + t837) * MDP(32) + t426 * MDP(33) + (-t617 ^ 2 - t799) * MDP(34) + (t428 * t617 - t637 + t686) * MDP(35);];
tauc  = t1;
