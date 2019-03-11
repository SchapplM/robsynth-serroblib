% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRPRR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:18:03
% EndTime: 2019-03-09 18:18:20
% DurationCPUTime: 10.39s
% Computational Cost: add. (10878->553), mult. (26128->727), div. (0->0), fcn. (20320->10), ass. (0->246)
t690 = sin(qJ(3));
t694 = cos(qJ(2));
t799 = cos(qJ(3));
t744 = qJD(1) * t799;
t691 = sin(qJ(2));
t759 = qJD(1) * t691;
t817 = -t690 * t759 + t694 * t744;
t627 = qJD(5) - t817;
t619 = qJD(6) + t627;
t778 = t690 * t694;
t635 = -qJD(1) * t778 - t691 * t744;
t683 = qJD(2) + qJD(3);
t686 = sin(pkin(11));
t687 = cos(pkin(11));
t603 = t635 * t687 - t683 * t686;
t605 = t635 * t686 + t687 * t683;
t689 = sin(qJ(5));
t693 = cos(qJ(5));
t547 = t603 * t689 + t605 * t693;
t692 = cos(qJ(6));
t688 = sin(qJ(6));
t811 = t603 * t693 - t605 * t689;
t788 = t811 * t688;
t821 = t547 * t692 + t788;
t816 = t619 * t821;
t822 = t547 * t627;
t716 = -t547 * t688 + t692 * t811;
t815 = t619 * t716;
t654 = t686 * t693 + t687 * t689;
t632 = t654 * qJD(5);
t814 = -t654 * t817 + t632;
t776 = t693 * t687;
t653 = t686 * t689 - t776;
t767 = t627 * t653;
t820 = t627 * t811;
t587 = -pkin(3) * t635 - qJ(4) * t817;
t800 = -pkin(8) - pkin(7);
t667 = t800 * t694;
t660 = qJD(1) * t667;
t638 = t690 * t660;
t666 = t800 * t691;
t658 = qJD(1) * t666;
t792 = qJD(2) * pkin(2);
t642 = t658 + t792;
t591 = t642 * t799 + t638;
t529 = t686 * t587 + t687 * t591;
t784 = t817 * t686;
t751 = pkin(9) * t784;
t505 = -t751 + t529;
t818 = -qJD(4) * t687 + t505;
t594 = -t653 * t688 + t654 * t692;
t771 = -qJD(6) * t594 + t688 * t767 - t692 * t814;
t593 = t692 * t653 + t654 * t688;
t770 = qJD(6) * t593 + t688 * t814 + t692 * t767;
t639 = t799 * t660;
t596 = t690 * t658 - t639;
t758 = qJD(3) * t690;
t724 = pkin(2) * t758 - t596;
t681 = -pkin(2) * t694 - pkin(1);
t665 = t681 * qJD(1);
t566 = -pkin(3) * t817 + qJ(4) * t635 + t665;
t592 = t690 * t642 - t639;
t576 = qJ(4) * t683 + t592;
t510 = t687 * t566 - t576 * t686;
t472 = -pkin(4) * t817 + pkin(9) * t603 + t510;
t511 = t686 * t566 + t687 * t576;
t483 = pkin(9) * t605 + t511;
t444 = t472 * t689 + t483 * t693;
t431 = pkin(10) * t547 + t444;
t754 = qJD(6) * t688;
t429 = t431 * t754;
t573 = -t683 * pkin(3) + qJD(4) - t591;
t534 = -pkin(4) * t605 + t573;
t473 = -pkin(5) * t547 + t534;
t813 = -t473 * t821 + t429;
t585 = t817 * t683;
t755 = qJD(5) * t693;
t787 = t585 * t686;
t468 = t585 * t776 + t605 * t755 + (qJD(5) * t603 - t787) * t689;
t656 = t691 * t799 + t778;
t600 = t683 * t656;
t586 = t600 * qJD(1);
t752 = qJD(1) * qJD(2);
t742 = t691 * t752;
t501 = pkin(2) * t742 + pkin(3) * t586 - qJ(4) * t585 + qJD(4) * t635;
t747 = qJD(2) * t800;
t727 = qJD(1) * t747;
t645 = t691 * t727;
t646 = t694 * t727;
t743 = qJD(3) * t799;
t699 = t642 * t743 + t645 * t799 + t690 * t646 + t660 * t758;
t527 = t683 * qJD(4) + t699;
t460 = t687 * t501 - t527 * t686;
t786 = t585 * t687;
t441 = pkin(4) * t586 - pkin(9) * t786 + t460;
t461 = t686 * t501 + t687 * t527;
t449 = -pkin(9) * t787 + t461;
t737 = t693 * t441 - t449 * t689;
t700 = -qJD(5) * t444 + t737;
t415 = pkin(5) * t586 - pkin(10) * t468 + t700;
t469 = -qJD(5) * t811 + t585 * t654;
t756 = qJD(5) * t689;
t709 = t689 * t441 + t693 * t449 + t472 * t755 - t483 * t756;
t416 = -pkin(10) * t469 + t709;
t738 = t692 * t415 - t688 * t416;
t812 = t473 * t716 + t738;
t810 = t586 * MDP(33) + (t716 ^ 2 - t821 ^ 2) * MDP(30) + t821 * MDP(29) * t716;
t809 = -0.2e1 * t752;
t808 = MDP(5) * (t691 ^ 2 - t694 ^ 2);
t719 = -t460 * t686 + t461 * t687;
t711 = -t690 * t691 + t694 * t799;
t590 = -pkin(3) * t711 - qJ(4) * t656 + t681;
t607 = t690 * t666 - t667 * t799;
t535 = t687 * t590 - t607 * t686;
t781 = t656 * t687;
t504 = -pkin(4) * t711 - pkin(9) * t781 + t535;
t536 = t686 * t590 + t687 * t607;
t782 = t656 * t686;
t519 = -pkin(9) * t782 + t536;
t769 = t689 * t504 + t693 * t519;
t612 = pkin(4) * t784;
t725 = -t612 + t724;
t677 = pkin(2) * t690 + qJ(4);
t640 = (-pkin(9) - t677) * t686;
t682 = t687 * pkin(9);
t641 = t677 * t687 + t682;
t764 = t689 * t640 + t693 * t641;
t807 = t814 * pkin(10);
t663 = (-pkin(9) - qJ(4)) * t686;
t664 = qJ(4) * t687 + t682;
t762 = t689 * t663 + t693 * t664;
t806 = -t635 * pkin(5) - t767 * pkin(10);
t805 = t799 * t666 + t690 * t667;
t804 = t814 * pkin(5);
t803 = qJD(1) * t656;
t735 = t688 * t468 + t692 * t469;
t424 = -qJD(6) * t716 + t735;
t528 = t687 * t587 - t591 * t686;
t783 = t817 * t687;
t726 = -t635 * pkin(4) - pkin(9) * t783;
t490 = t528 + t726;
t713 = qJD(4) * t686 + qJD(5) * t664;
t802 = -t663 * t755 + t818 * t693 + (t490 + t713) * t689;
t574 = pkin(2) * t759 + t587;
t597 = t658 * t799 + t638;
t525 = t687 * t574 - t597 * t686;
t487 = t525 + t726;
t526 = t686 * t574 + t687 * t597;
t503 = -t751 + t526;
t670 = pkin(2) * t743 + qJD(4);
t780 = t670 * t686;
t714 = qJD(5) * t641 + t780;
t801 = t693 * t503 - t640 * t755 - t670 * t776 + (t487 + t714) * t689;
t796 = pkin(10) * t654;
t794 = t653 * pkin(5);
t793 = t687 * pkin(4);
t443 = t693 * t472 - t483 * t689;
t430 = pkin(10) * t811 + t443;
t428 = pkin(5) * t627 + t430;
t791 = t428 * t692;
t790 = t431 * t692;
t599 = t683 * t711;
t785 = t599 * t686;
t779 = t670 * t687;
t695 = qJD(2) ^ 2;
t777 = t691 * t695;
t775 = t694 * t695;
t696 = qJD(1) ^ 2;
t774 = t694 * t696;
t750 = t691 * t792;
t520 = pkin(3) * t600 - qJ(4) * t599 - qJD(4) * t656 + t750;
t659 = t691 * t747;
t661 = t694 * t747;
t548 = qJD(3) * t805 + t799 * t659 + t690 * t661;
t464 = t686 * t520 + t687 * t548;
t766 = t804 + t725;
t753 = qJD(6) * t692;
t748 = t692 * t468 - t688 * t469 + t547 * t753;
t552 = t612 + t592;
t678 = -pkin(3) - t793;
t741 = -t552 + t804;
t739 = pkin(1) * t809;
t463 = t687 * t520 - t548 * t686;
t454 = pkin(4) * t600 - t599 * t682 + t463;
t462 = -pkin(9) * t785 + t464;
t736 = t693 * t454 - t462 * t689;
t734 = t693 * t504 - t519 * t689;
t533 = t642 * t758 + t690 * t645 - t799 * t646 - t660 * t743;
t733 = -t511 * t635 + t533 * t686;
t732 = t525 + t780;
t731 = -t526 + t779;
t730 = t693 * t640 - t641 * t689;
t729 = t693 * t663 - t664 * t689;
t728 = qJD(6) * t428 + t416;
t680 = -pkin(2) * t799 - pkin(3);
t482 = t693 * t487;
t644 = t653 * pkin(10);
t554 = -t644 + t764;
t723 = qJD(5) * t764 + qJD(6) * t554 - t503 * t689 + t654 * t670 + t482 + t806;
t486 = t693 * t490;
t570 = -t644 + t762;
t722 = t654 * qJD(4) + qJD(5) * t762 + qJD(6) * t570 - t505 * t689 + t486 + t806;
t553 = t730 - t796;
t721 = -qJD(6) * t553 + t801 + t807;
t569 = t729 - t796;
t720 = -qJD(6) * t569 + t802 + t807;
t491 = pkin(4) * t787 + t533;
t422 = t428 * t688 + t790;
t717 = t510 * t635 - t533 * t687;
t579 = t654 * t656;
t580 = t653 * t656;
t521 = t692 * t579 - t580 * t688;
t522 = -t579 * t688 - t580 * t692;
t712 = t510 * t783 + t511 * t784 + t719;
t568 = pkin(4) * t782 - t805;
t662 = t680 - t793;
t710 = t635 * t665 - t533;
t549 = t690 * t659 - t661 * t799 + t666 * t758 - t667 * t743;
t708 = t689 * t454 + t693 * t462 + t504 * t755 - t519 * t756;
t423 = t754 * t811 + t748;
t421 = -t431 * t688 + t791;
t445 = pkin(5) * t469 + t491;
t707 = t421 * t635 + t445 * t593 - t473 * t771;
t706 = -t422 * t635 + t445 * t594 - t473 * t770;
t705 = t443 * t635 + t491 * t653 + t534 * t814;
t704 = -t444 * t635 + t491 * t654 - t534 * t767;
t509 = pkin(4) * t785 + t549;
t703 = t533 * t656 + t573 * t599 - t585 * t805;
t702 = -pkin(3) * t585 - qJ(4) * t586 - (-qJD(4) + t573) * t817;
t701 = t585 * t680 - t586 * t677 - (t573 - t670) * t817;
t698 = -t665 * t817 - t699;
t697 = (-t423 * t593 - t424 * t594 - t716 * t771 - t770 * t821) * MDP(30) + (t423 * t594 + t716 * t770) * MDP(29) + (-t468 * t653 - t469 * t654 - t547 * t767 + t811 * t814) * MDP(23) + (t586 * t594 - t619 * t770 - t635 * t716) * MDP(31) + (-t586 * t593 + t619 * t771 + t635 * t821) * MDP(32) + (t468 * t654 + t767 * t811) * MDP(22) + (t586 * t654 - t627 * t767 - t635 * t811) * MDP(24) + (t547 * t635 - t586 * t653 - t627 * t814) * MDP(25) + t585 * MDP(13) + (t635 ^ 2 - t817 ^ 2) * MDP(12) + (MDP(11) * t817 + MDP(26) * t627 + MDP(33) * t619) * t635 + (-t817 * MDP(13) + (-t635 - t803) * MDP(14)) * t683;
t610 = t678 + t794;
t608 = t662 + t794;
t557 = t586 * t711;
t523 = t579 * pkin(5) + t568;
t489 = t599 * t654 + t755 * t781 - t756 * t782;
t488 = -t599 * t653 - t632 * t656;
t455 = t489 * pkin(5) + t509;
t447 = -pkin(10) * t579 + t769;
t446 = -pkin(5) * t711 + pkin(10) * t580 + t734;
t433 = qJD(6) * t522 + t488 * t688 + t692 * t489;
t432 = -qJD(6) * t521 + t488 * t692 - t489 * t688;
t419 = -pkin(10) * t489 + t708;
t418 = pkin(5) * t600 - pkin(10) * t488 - qJD(5) * t769 + t736;
t1 = [(t460 * t535 + t461 * t536 + t463 * t510 + t464 * t511 - t533 * t805 + t549 * t573) * MDP(21) + (-t468 * t580 - t488 * t811) * MDP(22) + (-t444 * t600 + t568 * t468 + t534 * t488 - t491 * t580 - t509 * t811 - t586 * t769 - t627 * t708 + t709 * t711) * MDP(28) + (MDP(13) * t599 - MDP(14) * t600 - MDP(16) * t549 - MDP(17) * t548) * t683 + (t423 * t522 - t432 * t716) * MDP(29) + (t463 * t603 + t464 * t605 + (-t460 * t656 - t510 * t599 - t535 * t585) * t687 + (-t461 * t656 - t511 * t599 - t536 * t585) * t686) * MDP(20) + (-t422 * t600 + t523 * t423 - t429 * t711 + t473 * t432 + t445 * t522 - t455 * t716 + (-(-qJD(6) * t447 + t418) * t619 - t446 * t586 + t415 * t711) * t688 + (-(qJD(6) * t446 + t419) * t619 - t447 * t586 + t728 * t711) * t692) * MDP(35) + (-t423 * t711 + t432 * t619 + t522 * t586 - t600 * t716) * MDP(31) + (-t423 * t521 - t424 * t522 + t432 * t821 + t433 * t716) * MDP(30) + ((t418 * t692 - t419 * t688) * t619 + (t446 * t692 - t447 * t688) * t586 - t738 * t711 + t421 * t600 - t455 * t821 + t523 * t424 + t445 * t521 + t473 * t433 + ((-t446 * t688 - t447 * t692) * t619 + t422 * t711) * qJD(6)) * MDP(34) + (t424 * t711 - t433 * t619 - t521 * t586 + t600 * t821) * MDP(32) + (-t468 * t711 + t488 * t627 - t580 * t586 - t600 * t811) * MDP(24) + (-pkin(7) * t775 + t691 * t739) * MDP(9) - MDP(7) * t777 + (pkin(7) * t777 + t694 * t739) * MDP(10) + (-t460 * t711 - t463 * t817 + t510 * t600 + t535 * t586 - t549 * t605 + t686 * t703) * MDP(18) + (t461 * t711 + t464 * t817 - t511 * t600 - t536 * t586 - t549 * t603 + t687 * t703) * MDP(19) + (t586 * t681 + t600 * t665 + (-qJD(1) * t711 - t817) * t750) * MDP(16) + (t585 * t711 - t586 * t656 + t599 * t817 + t600 * t635) * MDP(12) + (t585 * t681 + t599 * t665 + (-t635 + t803) * t750) * MDP(17) + 0.2e1 * t694 * MDP(4) * t742 + (-t468 * t579 + t469 * t580 + t488 * t547 + t489 * t811) * MDP(23) + (t469 * t711 - t489 * t627 + t547 * t600 - t579 * t586) * MDP(25) + (t736 * t627 + t734 * t586 - t737 * t711 + t443 * t600 - t509 * t547 + t568 * t469 + t491 * t579 + t534 * t489 + (t444 * t711 - t627 * t769) * qJD(5)) * MDP(27) + (t600 * t627 - t557) * MDP(26) + (t600 * t619 - t557) * MDP(33) + (t585 * t656 - t599 * t635) * MDP(11) + t808 * t809 + MDP(6) * t775; t697 + (t596 * t683 + (-t683 * t758 + t759 * t817) * pkin(2) + t710) * MDP(16) + t696 * t808 - t691 * MDP(4) * t774 + (t525 * t817 - t605 * t724 + t686 * t701 + t717) * MDP(18) + (-t603 * t732 + t605 * t731 + t712) * MDP(20) + (-t526 * t817 - t603 * t724 + t687 * t701 + t733) * MDP(19) + (t662 * t468 - t764 * t586 + t627 * t801 - t725 * t811 + t704) * MDP(28) + (-(t553 * t688 + t554 * t692) * t586 + t608 * t423 + (t688 * t723 + t692 * t721) * t619 - t766 * t716 + t706) * MDP(35) + ((t553 * t692 - t554 * t688) * t586 + t608 * t424 + (t688 * t721 - t692 * t723) * t619 - t766 * t821 + t707) * MDP(34) + (-t510 * t732 + t511 * t731 + t533 * t680 + t573 * t724 + t677 * t719) * MDP(21) + (t730 * t586 + t662 * t469 + (-t482 - t714 * t693 + (-qJD(5) * t640 + t503 - t779) * t689) * t627 - t725 * t547 + t705) * MDP(27) + (t597 * t683 + (t635 * t759 - t683 * t743) * pkin(2) + t698) * MDP(17) + (MDP(9) * t691 * t696 + MDP(10) * t774) * pkin(1); t697 + (t528 * t817 + t592 * t605 + t686 * t702 + t717) * MDP(18) + (-t529 * t817 + t592 * t603 + t687 * t702 + t733) * MDP(19) + (-pkin(3) * t533 - t510 * t528 - t511 * t529 - t573 * t592 + (-t510 * t686 + t511 * t687) * qJD(4) + t719 * qJ(4)) * MDP(21) + (-t528 * t603 - t529 * t605 + (-t603 * t686 + t605 * t687) * qJD(4) + t712) * MDP(20) + ((t569 * t692 - t570 * t688) * t586 + t610 * t424 + (t688 * t720 - t692 * t722) * t619 - t741 * t821 + t707) * MDP(34) + (-(t569 * t688 + t570 * t692) * t586 + t610 * t423 + (t688 * t722 + t692 * t720) * t619 - t741 * t716 + t706) * MDP(35) + (t678 * t468 + t552 * t811 - t762 * t586 + t627 * t802 + t704) * MDP(28) + (t729 * t586 + t678 * t469 + t552 * t547 + (-t486 - t713 * t693 + (-qJD(5) * t663 + t818) * t689) * t627 + t705) * MDP(27) + (t591 * t683 + t698) * MDP(17) + (t592 * t683 + t710) * MDP(16); (t603 * t817 + t787) * MDP(18) + (-t605 * t817 + t786) * MDP(19) + (-t603 ^ 2 - t605 ^ 2) * MDP(20) + (-t510 * t603 - t511 * t605 + t533) * MDP(21) + (t469 - t820) * MDP(27) + (t468 + t822) * MDP(28) + (t424 - t815) * MDP(34) + (t423 + t816) * MDP(35); t811 * t547 * MDP(22) + (-t547 ^ 2 + t811 ^ 2) * MDP(23) + (t468 - t822) * MDP(24) + (-t469 - t820) * MDP(25) + t586 * MDP(26) + (t444 * t627 + t534 * t811 + t700) * MDP(27) + (t443 * t627 - t534 * t547 - t709) * MDP(28) + (t423 - t816) * MDP(31) + (-t424 - t815) * MDP(32) + (-(-t430 * t688 - t790) * t619 - t422 * qJD(6) + (t586 * t692 - t619 * t754 - t811 * t821) * pkin(5) + t812) * MDP(34) + ((-t431 * t619 - t415) * t688 + (t430 * t619 - t728) * t692 + (-t586 * t688 - t619 * t753 - t716 * t811) * pkin(5) + t813) * MDP(35) + t810; (t748 - t816) * MDP(31) + (-t735 - t815) * MDP(32) + (t422 * t619 + t812) * MDP(34) + (-t688 * t415 - t692 * t416 + t421 * t619 + t813) * MDP(35) + (MDP(31) * t788 + MDP(32) * t716 - MDP(34) * t422 - MDP(35) * t791) * qJD(6) + t810;];
tauc  = t1;
