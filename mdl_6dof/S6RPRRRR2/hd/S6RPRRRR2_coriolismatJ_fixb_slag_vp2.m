% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRRR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:57:14
% EndTime: 2019-03-09 06:57:36
% DurationCPUTime: 14.12s
% Computational Cost: add. (37413->715), mult. (72227->911), div. (0->0), fcn. (78892->10), ass. (0->424)
t839 = qJD(3) + qJD(4);
t477 = sin(qJ(5));
t473 = t477 ^ 2;
t478 = cos(qJ(5));
t474 = t478 ^ 2;
t641 = t473 + t474;
t476 = sin(qJ(6));
t744 = cos(qJ(6));
t745 = cos(qJ(4));
t561 = t745 * t744;
t605 = t478 * t745;
t394 = (-t476 * t605 - t477 * t561) * pkin(3);
t606 = t477 * t745;
t395 = (-t476 * t606 + t478 * t561) * pkin(3);
t634 = t745 * pkin(3);
t567 = -t634 / 0.2e1;
t550 = t478 * t567;
t794 = m(7) * pkin(5);
t639 = t794 / 0.2e1;
t694 = t477 * mrSges(6,1);
t876 = -t395 / 0.2e1;
t878 = mrSges(7,1) / 0.2e1;
t826 = mrSges(7,2) * t876 + t394 * t878;
t891 = (t394 * t744 + t395 * t476) * t639 + t567 * t694 + mrSges(6,2) * t550 + t826;
t742 = sin(qJ(4));
t743 = sin(qJ(3));
t746 = cos(qJ(3));
t443 = -t742 * t746 - t743 * t745;
t735 = t443 * pkin(4);
t442 = t742 * t743 - t745 * t746;
t736 = t442 * pkin(9);
t367 = -t735 + t736;
t559 = sin(pkin(11)) * pkin(1) + pkin(7);
t529 = t746 * t559;
t423 = pkin(8) * t746 + t529;
t528 = t743 * t559;
t511 = -pkin(8) * t743 - t528;
t345 = t423 * t742 - t745 * t511;
t848 = t477 * t345;
t206 = t478 * t367 + t848;
t666 = t442 * t478;
t560 = -t443 * pkin(5) + pkin(10) * t666;
t132 = t206 + t560;
t847 = t478 * t345;
t207 = t477 * t367 - t847;
t667 = t442 * t477;
t637 = pkin(10) * t667;
t163 = t637 + t207;
t103 = t132 * t744 - t476 * t163;
t104 = t476 * t132 + t163 * t744;
t792 = -mrSges(6,2) / 0.2e1;
t793 = mrSges(6,1) / 0.2e1;
t877 = -mrSges(7,2) / 0.2e1;
t832 = t103 * t878 + t104 * t877;
t890 = -t206 * t793 - t207 * t792 - (t103 * t744 + t104 * t476) * t639 - t832;
t633 = t743 * pkin(3);
t353 = t633 + t367;
t196 = t478 * t353 + t848;
t197 = t477 * t353 - t847;
t131 = t196 + t560;
t153 = t637 + t197;
t91 = t131 * t744 - t476 * t153;
t92 = t476 * t131 + t153 * t744;
t833 = t92 * t877 + t878 * t91;
t889 = -t196 * t793 - t197 * t792 - (t476 * t92 + t744 * t91) * t639 - t833;
t533 = t476 * t478 + t477 * t744;
t332 = t533 * t442;
t231 = mrSges(7,2) * t443 + mrSges(7,3) * t332;
t817 = -t476 * t477 + t744 * t478;
t334 = t817 * t442;
t233 = -mrSges(7,1) * t443 + mrSges(7,3) * t334;
t635 = pkin(5) * t744;
t738 = pkin(5) * t476;
t778 = -t334 / 0.2e1;
t781 = t332 / 0.2e1;
t539 = Ifges(7,5) * t778 + Ifges(7,6) * t781;
t875 = t443 / 0.2e1;
t815 = Ifges(7,3) * t875 - t539;
t859 = -t443 / 0.2e1;
t888 = t815 - t233 * t635 / 0.2e1 - t231 * t738 / 0.2e1 - Ifges(6,3) * t859;
t632 = t742 * pkin(3);
t466 = t632 + pkin(9);
t727 = pkin(10) + t466;
t437 = t727 * t477;
t438 = t727 * t478;
t351 = -t476 * t437 + t438 * t744;
t573 = -t744 * t437 - t438 * t476;
t887 = -t351 * mrSges(7,1) - t573 * mrSges(7,2);
t790 = -pkin(10) - pkin(9);
t459 = t790 * t477;
t460 = t790 * t478;
t382 = t476 * t459 - t460 * t744;
t572 = t744 * t459 + t460 * t476;
t886 = -t382 * mrSges(7,1) - t572 * mrSges(7,2);
t885 = t832 - t815;
t884 = t833 - t815;
t149 = -Ifges(7,4) * t334 + Ifges(7,2) * t332 - Ifges(7,6) * t443;
t151 = -Ifges(7,1) * t334 + Ifges(7,4) * t332 - t443 * Ifges(7,5);
t470 = Ifges(6,4) * t478;
t552 = Ifges(6,2) * t477 - t470;
t266 = -Ifges(6,6) * t443 + t442 * t552;
t723 = Ifges(6,4) * t477;
t455 = Ifges(6,1) * t478 - t723;
t268 = -Ifges(6,5) * t443 - t442 * t455;
t719 = Ifges(6,2) * t478;
t452 = t719 + t723;
t551 = Ifges(6,5) * t477 + Ifges(6,6) * t478;
t591 = -t666 / 0.2e1;
t592 = t667 / 0.2e1;
t747 = t478 / 0.2e1;
t749 = t477 / 0.2e1;
t760 = t533 / 0.2e1;
t762 = t817 / 0.2e1;
t429 = Ifges(7,4) * t817;
t366 = Ifges(7,1) * t533 + t429;
t768 = t366 / 0.2e1;
t721 = Ifges(7,4) * t533;
t364 = Ifges(7,2) * t817 + t721;
t771 = t364 / 0.2e1;
t724 = Ifges(6,1) * t477;
t824 = t470 + t724;
t508 = -Ifges(5,5) * t442 + Ifges(5,6) * t443 + t149 * t762 + t151 * t760 + t266 * t747 + t268 * t749 + t332 * t771 - t334 * t768 + t452 * t592 + t824 * t591 + (Ifges(7,5) * t533 + Ifges(7,6) * t817 + t551) * t859;
t813 = -mrSges(6,1) * t478 + t477 * mrSges(6,2);
t823 = t745 * t423 + t742 * t511;
t846 = t823 * t813;
t852 = t823 * mrSges(5,1);
t855 = t345 * mrSges(5,2);
t362 = -mrSges(7,1) * t817 + mrSges(7,2) * t533;
t417 = pkin(5) * t667;
t844 = -t417 + t823;
t867 = t844 * t362;
t883 = t508 + t846 + t855 - t852 + t867;
t590 = t666 / 0.2e1;
t711 = t334 * mrSges(7,2);
t712 = t332 * mrSges(7,1);
t650 = t712 / 0.2e1 + t711 / 0.2e1;
t497 = (t332 * t744 - t334 * t476) * t639 + mrSges(6,1) * t592 + mrSges(6,2) * t590 + t650;
t360 = mrSges(7,1) * t533 + mrSges(7,2) * t817;
t671 = t442 * t360;
t331 = t817 * t443;
t695 = t533 * mrSges(7,3);
t613 = -t695 / 0.2e1;
t614 = t695 / 0.2e1;
t821 = (t613 + t614) * t331;
t857 = m(7) * t417;
t692 = t478 * mrSges(6,2);
t451 = t692 + t694;
t670 = t442 * t451;
t858 = -t670 / 0.2e1;
t865 = -t857 / 0.2e1 - t671 / 0.2e1 + t858 + t497 + t821;
t882 = qJD(2) * t865;
t881 = t855 / 0.2e1 + t867 / 0.2e1;
t609 = t671 / 0.2e1 + t821;
t549 = t609 + t650;
t333 = t533 * t443;
t378 = t442 * t443;
t60 = m(7) * (t331 * t334 + t332 * t333 - t378) + m(6) * (t641 - 0.1e1) * t378;
t651 = t60 * qJD(2);
t864 = t857 / 0.2e1 + t670 / 0.2e1 + t497 + t609;
t880 = t864 * qJD(5) + t549 * qJD(6) + t651;
t548 = t609 - t650;
t879 = -qJD(5) * t865 + qJD(6) * t548 - t651;
t607 = -cos(pkin(11)) * pkin(1) - pkin(2);
t447 = -pkin(3) * t746 + t607;
t326 = t442 * pkin(4) + t443 * pkin(9) + t447;
t184 = t478 * t326 - t477 * t823;
t664 = t443 * t478;
t134 = pkin(10) * t664 + t184;
t124 = t442 * pkin(5) + t134;
t185 = t477 * t326 + t478 * t823;
t665 = t443 * t477;
t135 = pkin(10) * t665 + t185;
t657 = t476 * t135;
t84 = t124 * t744 - t657;
t96 = t134 * t744 - t657;
t874 = t84 - t96;
t873 = mrSges(6,3) * t641;
t246 = -pkin(5) * t665 + t345;
t872 = t246 * t844;
t467 = -t634 - pkin(4);
t734 = t478 * pkin(5);
t448 = t467 - t734;
t871 = t448 * t844;
t468 = -pkin(4) - t734;
t870 = t468 * t844;
t642 = Ifges(7,5) * t817 - Ifges(7,6) * t533;
t64 = t642 + t887;
t869 = t64 * qJD(6);
t74 = t642 + t886;
t868 = t74 * qJD(6);
t186 = t331 * mrSges(7,1) - t333 * mrSges(7,2);
t866 = qJD(6) * t186;
t863 = t846 / 0.2e1 - t852 / 0.2e1;
t722 = Ifges(7,4) * t331;
t150 = Ifges(7,2) * t333 + t442 * Ifges(7,6) - t722;
t316 = Ifges(7,4) * t333;
t152 = -Ifges(7,1) * t331 + t442 * Ifges(7,5) + t316;
t187 = -t711 - t712;
t188 = -mrSges(7,1) * t333 - mrSges(7,2) * t331;
t356 = mrSges(6,2) * t443 + mrSges(6,3) * t667;
t358 = -t443 * mrSges(6,1) + mrSges(6,3) * t666;
t361 = -t443 * mrSges(5,1) - t442 * mrSges(5,2);
t779 = t333 / 0.2e1;
t783 = -t331 / 0.2e1;
t604 = t744 * t135;
t85 = t476 * t124 + t604;
t862 = t149 * t779 + t150 * t781 + t151 * t783 + t152 * t778 + t184 * t358 + t185 * t356 + t246 * t187 + t844 * t188 + t85 * t231 + t84 * t233 - t345 * t670 + t447 * t361;
t861 = -t246 * t443 + t332 * t84 - t334 * t85 + t442 * t844;
t774 = t351 / 0.2e1;
t763 = t382 / 0.2e1;
t856 = pkin(4) * t823;
t95 = -t476 * t134 - t604;
t834 = t85 + t95;
t469 = Ifges(6,5) * t478;
t718 = Ifges(6,6) * t477;
t851 = Ifges(5,4) - t469 / 0.2e1 + t718 / 0.2e1;
t850 = t345 * t742;
t677 = t345 * t823;
t849 = t467 * t823;
t845 = -t362 - t813;
t841 = t641 * t442;
t652 = t478 * t356;
t654 = t477 * t358;
t535 = t652 / 0.2e1 - t654 / 0.2e1;
t684 = t197 * t478;
t685 = t196 * t477;
t542 = t684 - t685;
t659 = t468 * t187;
t672 = t382 * t231;
t673 = t572 * t233;
t690 = t92 * t817;
t691 = t91 * t533;
t740 = pkin(4) * t670;
t797 = m(7) / 0.2e1;
t799 = m(6) / 0.2e1;
t840 = (pkin(9) * t542 - t856) * t799 + (t382 * t92 + t572 * t91 + t870) * t797 + t740 / 0.2e1 + t673 / 0.2e1 + t672 / 0.2e1 + t659 / 0.2e1 + (t684 / 0.2e1 - t685 / 0.2e1) * mrSges(6,3) + (-t691 / 0.2e1 + t690 / 0.2e1) * mrSges(7,3) + t535 * pkin(9) + t881;
t757 = t442 / 0.4e1;
t521 = t246 * t360 / 0.2e1 + t642 * t757;
t574 = t469 - t718;
t365 = Ifges(7,1) * t817 - t721;
t577 = t364 / 0.4e1 - t365 / 0.4e1;
t363 = -Ifges(7,2) * t533 + t429;
t578 = t363 / 0.4e1 + t366 / 0.4e1;
t269 = Ifges(6,5) * t442 - t443 * t455;
t653 = t478 * t269;
t267 = Ifges(6,6) * t442 + t443 * t552;
t655 = t477 * t267;
t737 = pkin(5) * t477;
t754 = t451 / 0.2e1;
t838 = t345 * t754 + t574 * t757 - t655 / 0.4e1 + t653 / 0.4e1 + t188 * t737 / 0.2e1 + t521 + t578 * t333 + t577 * t331;
t686 = t185 * t478;
t837 = (t184 * t477 - t686 + t823) * t442;
t748 = -t478 / 0.2e1;
t780 = -t333 / 0.2e1;
t782 = t331 / 0.2e1;
t836 = -t451 * t823 - t851 * t443 + (-Ifges(6,3) - Ifges(7,3) - Ifges(5,2) + Ifges(5,1)) * t442 + Ifges(7,5) * t782 + Ifges(7,6) * t780 + t266 * t749 + t268 * t748;
t531 = t813 * t443;
t189 = Ifges(7,2) * t331 + t316;
t827 = t152 + t189;
t825 = t824 - t552;
t190 = Ifges(7,1) * t333 + t722;
t557 = (t150 / 0.4e1 - t190 / 0.4e1) * t533;
t822 = (t189 / 0.4e1 + t152 / 0.4e1) * t817 - t557;
t631 = mrSges(6,3) * t665;
t357 = -mrSges(6,2) * t442 + t631;
t359 = t442 * mrSges(6,1) + mrSges(6,3) * t664;
t750 = -t477 / 0.2e1;
t534 = t357 * t750 + t359 * t748;
t766 = -t572 / 0.2e1;
t819 = t331 * t763 + t333 * t766;
t777 = -t573 / 0.2e1;
t818 = t331 * t774 + t333 * t777;
t812 = qJD(2) * t548;
t725 = mrSges(7,3) * t333;
t232 = -mrSges(7,2) * t442 + t725;
t751 = -t468 / 0.2e1;
t587 = t186 * t751;
t765 = t572 / 0.2e1;
t713 = t331 * mrSges(7,3);
t234 = mrSges(7,1) * t442 + t713;
t784 = t234 / 0.2e1;
t810 = t232 * t765 - t382 * t784 + t587;
t755 = -t448 / 0.2e1;
t588 = t186 * t755;
t776 = t573 / 0.2e1;
t809 = t232 * t776 - t351 * t784 + t588;
t806 = t873 * t875 + t534;
t803 = t331 ^ 2;
t802 = t333 ^ 2;
t801 = 2 * qJD(4);
t798 = -m(7) / 0.2e1;
t796 = -pkin(4) / 0.2e1;
t795 = m(6) * pkin(3);
t788 = -t188 / 0.2e1;
t786 = t232 / 0.2e1;
t785 = -t234 / 0.2e1;
t775 = -t351 / 0.2e1;
t773 = -t362 / 0.2e1;
t764 = -t382 / 0.2e1;
t759 = -t442 / 0.2e1;
t758 = t442 / 0.2e1;
t753 = -t455 / 0.4e1;
t752 = t467 / 0.2e1;
t739 = pkin(4) * t451;
t733 = t84 * mrSges(7,2);
t732 = t85 * mrSges(7,1);
t729 = t95 * mrSges(7,1);
t728 = t96 * mrSges(7,2);
t696 = t817 * mrSges(7,3);
t584 = t359 * t749;
t484 = (t788 + (t692 / 0.2e1 + t694 / 0.2e1) * t443 - t535) * t443 + (t357 * t748 + t584 + t187 / 0.2e1 + t858) * t442 + t231 * t783 + t234 * t781 + t233 * t779 + t232 * t778;
t15 = (-t331 * t92 + t333 * t91 + t861) * t797 + ((-t345 - t542) * t443 + t837) * t799 + t484;
t689 = t15 * qJD(1);
t682 = t207 * t478;
t683 = t206 * t477;
t541 = t682 - t683;
t17 = (t103 * t333 - t104 * t331 + t861) * t797 + ((-t345 - t541) * t443 + t837) * t799 + t484;
t688 = t17 * qJD(1);
t518 = t186 * t759 + t232 * t779 + t234 * t782;
t558 = mrSges(6,3) * (t474 / 0.2e1 + t473 / 0.2e1);
t18 = (t331 * t874 + t834 * t333) * t798 + (t802 / 0.2e1 + t803 / 0.2e1) * mrSges(7,3) + (t443 * t558 + t590 * t794 - t758 * t813 + t534) * t443 - t518;
t687 = t18 * qJD(1);
t21 = t518 - (t802 + t803) * mrSges(7,3) / 0.2e1;
t681 = t21 * qJD(1);
t675 = t573 * t233;
t674 = t351 * t231;
t663 = t448 * t187;
t662 = t448 * t360;
t661 = t467 * t670;
t660 = t467 * t451;
t658 = t468 * t360;
t649 = Ifges(7,5) * t333 + Ifges(7,6) * t331;
t640 = mrSges(7,3) * t738;
t638 = pkin(5) * t664;
t630 = mrSges(6,3) * t683;
t629 = mrSges(6,3) * t682;
t624 = -t84 / 0.2e1 + t96 / 0.2e1;
t623 = t95 / 0.2e1 + t85 / 0.2e1;
t622 = t466 * t654;
t621 = t466 * t652;
t616 = -t696 / 0.2e1;
t615 = t696 / 0.2e1;
t610 = t332 * t573 - t351 * t334 - t448 * t443;
t608 = -t467 * t443 - t466 * t841;
t583 = t452 * t750;
t580 = t777 + t776;
t579 = t775 + t774;
t576 = t766 + t765;
t575 = t764 + t763;
t570 = mrSges(7,3) * t635;
t569 = t443 * t632;
t419 = t442 * t632;
t555 = (t448 / 0.2e1 + t468 / 0.2e1) * t360;
t494 = -t653 / 0.2e1 + t655 / 0.2e1 + t539 + t851 * t442;
t527 = mrSges(4,1) * t743 + mrSges(4,2) * t746;
t7 = (-mrSges(5,2) * t633 + t836) * t443 + m(6) * (t184 * t196 + t185 * t197 + t677) + m(7) * (t84 * t91 + t85 * t92 + t872) + m(5) * t447 * t633 + (mrSges(5,1) * t633 + t494) * t442 + t607 * t527 + t197 * t357 + t196 * t359 + t92 * t232 + t91 * t234 + (-Ifges(4,2) + Ifges(4,1)) * t746 * t743 + (-t743 ^ 2 + t746 ^ 2) * Ifges(4,4) + t862;
t547 = t7 * qJD(1) + t15 * qJD(2);
t9 = t836 * t443 + m(6) * (t184 * t206 + t185 * t207 + t677) + m(7) * (t103 * t84 + t104 * t85 + t872) + t494 * t442 + t207 * t357 + t206 * t359 + t103 * t234 + t104 * t232 + t862;
t546 = t9 * qJD(1) + t17 * qJD(2);
t530 = -t246 * t186 + t649 * t758 + t85 * t713;
t10 = -t95 * t234 + t188 * t638 - t96 * t232 - m(7) * (t84 * t95 + t85 * t96) - t345 * t531 + t185 * t359 + (-t150 / 0.2e1 + t190 / 0.2e1) * t331 + (t84 * mrSges(7,3) - t189 / 0.2e1 - t152 / 0.2e1) * t333 + (t269 * t750 + t267 * t748 + m(7) * t246 * t734 + t551 * t759 - mrSges(6,3) * t686 + (t747 * t824 + t583) * t443) * t443 - t530 + (t631 - t357) * t184;
t545 = -t10 * qJD(1) - t18 * qJD(2);
t16 = -t234 * t85 + t150 * t782 + t190 * t783 + (t232 - t725) * t84 + t530 + t827 * t779;
t544 = t16 * qJD(1) + t21 * qJD(2);
t537 = -t552 / 0.4e1 + t824 / 0.4e1 + t724 / 0.4e1;
t526 = t641 * t745;
t479 = -Ifges(6,6) * t667 / 0.2e1 + t96 * t615 + t84 * t616 + Ifges(6,5) * t590 + t638 * t773 + t827 * t817 / 0.4e1 + (t452 / 0.2e1 + t753) * t664 + t834 * t613 - t557 + (t824 + t825) * t665 / 0.4e1 + t838 + t888;
t230 = t246 * t737;
t513 = -t351 * t874 + t834 * t573 + t230;
t3 = t479 + t531 * t752 + (-t448 * t638 + t513) * t797 + t818 * mrSges(7,3) + t806 * t466 + t809 + t889;
t435 = t448 * t737;
t514 = -(-t365 / 0.2e1 + t771) * t533 + (t768 + t363 / 0.2e1) * t817;
t493 = (t824 / 0.2e1 - t552 / 0.2e1) * t478 + (pkin(5) * t362 + t455 / 0.2e1 - t452 / 0.2e1) * t477 + t514;
t39 = m(7) * t435 + t493 + t660 + t662;
t524 = t3 * qJD(1) + t39 * qJD(3) - t882;
t486 = (-pkin(3) * t443 * t526 + t419 + t608) * t799 + (-t331 * t395 + t333 * t394 + t419 + t610) * t797;
t496 = (t332 * t572 - t334 * t382 - t443 * t468) * t797 + (-t641 * t736 + t735) * t799;
t40 = t486 - t496;
t487 = -t394 * t695 + t395 * t696 + (-mrSges(5,1) - t845) * t632 + (-mrSges(5,2) + t873) * t634;
t62 = m(7) * (t351 * t395 + t394 * t573 + t448 * t632) + (t466 * t526 + t467 * t742) * t795 + t487;
t480 = t630 / 0.2e1 - t629 / 0.2e1 - t621 / 0.2e1 + t622 / 0.2e1 - t863 + t104 * t616 - m(6) * (t849 + t541 * t466 + (-t184 * t606 + t185 * t605 + t850) * pkin(3)) / 0.2e1 + t103 * t614 + t661 / 0.2e1 - t674 / 0.2e1 - t675 / 0.2e1 + t394 * t785 + t632 * t788 + t569 * t754 - t663 / 0.2e1 + t357 * t550 + t232 * t876 + (t103 * t573 + t351 * t104 + t246 * t632 + t394 * t84 + t395 * t85 + t871) * t798 + t584 * t634 - t881;
t8 = t480 + t840 + t863;
t523 = -t8 * qJD(1) + t40 * qJD(2) + t62 * qJD(3);
t498 = t521 + t822;
t483 = (mrSges(7,3) * t774 + t577) * t331 + (mrSges(7,3) * t777 + t578) * t333 + t573 * t786 + t234 * t775 + t588 + t498;
t12 = t483 - t884;
t51 = t514 + t662;
t522 = t12 * qJD(1) + t51 * qJD(3) + t812;
t519 = t452 / 0.4e1 + t753 + t723 / 0.2e1 + t719 / 0.4e1;
t517 = -t533 * t623 + t624 * t817;
t449 = t468 * t737;
t503 = (t435 + t449) * t797;
t24 = t503 + t362 * t737 + t455 * t749 + t583 + t660 / 0.2e1 + t658 / 0.2e1 + t662 / 0.2e1 - t739 / 0.2e1 - t533 * t771 + t365 * t760 + (t366 + t363) * t762 + t825 * t747 + (t616 + t615) * (t572 + t573) - t891;
t43 = m(7) * t449 + t493 + t658 - t739;
t512 = -t382 * t874 + t834 * t572 + t230;
t5 = t479 + (-t468 * t638 + t512) * t797 + t531 * t796 + t819 * mrSges(7,3) + t806 * pkin(9) + t810 + t890;
t516 = t5 * qJD(1) + t24 * qJD(3) + t43 * qJD(4) - t882;
t482 = (mrSges(7,3) * t766 + t578) * t333 + (mrSges(7,3) * t763 + t577) * t331 + t572 * t786 + t234 * t764 + t587 + t498;
t14 = t482 - t885;
t500 = t555 + t514;
t35 = t500 - t826;
t56 = t514 + t658;
t515 = t14 * qJD(1) + t35 * qJD(3) + t56 * qJD(4) + t812;
t510 = -t533 * t640 - t570 * t817 + t574 + t642;
t509 = -mrSges(6,3) * t841 - t332 * t695 - t334 * t696 + t443 * t845 - t361;
t499 = (t476 * t785 + t744 * t786 + (t476 * t782 + t744 * t780) * mrSges(7,3)) * pkin(5);
t20 = -mrSges(7,1) * t623 + mrSges(7,2) * t624 + t499;
t446 = (mrSges(7,1) * t476 + mrSges(7,2) * t744) * pkin(5);
t68 = mrSges(7,1) * t579 + mrSges(7,2) * t580;
t88 = mrSges(7,1) * t575 + mrSges(7,2) * t576;
t504 = -t20 * qJD(1) - t68 * qJD(3) - t88 * qJD(4) + t446 * qJD(5);
t481 = Ifges(6,5) * t591 + Ifges(6,6) * t592 + t822 + t838 - t888;
t439 = t446 * qJD(6);
t36 = t500 + t826;
t34 = t486 + t496 + t509;
t25 = t503 + (t752 + t796) * t451 + t555 + (-(-t575 - t579) * t533 + (t576 + t580) * t817) * mrSges(7,3) + t493 + t891;
t19 = -t733 / 0.2e1 - t732 / 0.2e1 - t728 / 0.2e1 + t729 / 0.2e1 + t499 + t649;
t13 = t482 + t885;
t11 = t483 + t884;
t6 = t481 + t512 * t797 + t534 * pkin(9) + (t517 + t819) * mrSges(7,3) + (pkin(9) * t558 + (pkin(4) * t792 + t537) * t477 + (pkin(4) * t793 + (m(7) * t751 + t773) * pkin(5) + t519) * t478) * t443 + t810 - t890;
t4 = t481 + t513 * t797 + (t466 * t558 + (mrSges(6,2) * t752 + t537) * t477 + (-t467 * mrSges(6,1) / 0.2e1 + (m(7) * t755 + t773) * pkin(5) + t519) * t478) * t443 + (t517 + t818) * mrSges(7,3) + t534 * t466 + t809 - t889;
t2 = -t480 + t508 + (t813 / 0.2e1 - mrSges(5,1) / 0.2e1) * t823 + t840;
t1 = qJD(3) * t15 + qJD(4) * t17 - qJD(5) * t18 + qJD(6) * t21;
t22 = [qJD(3) * t7 + qJD(4) * t9 - qJD(5) * t10 + qJD(6) * t16, t1 (Ifges(4,5) * t746 - Ifges(4,6) * t743 + t621 - t622 + m(5) * (-t745 * t823 - t850) * pkin(3) - t661 + mrSges(4,2) * t528 + t674 + t675 + t663 + m(6) * (t466 * t542 + t849) + (t442 * t634 + t569) * mrSges(5,3) + t542 * mrSges(6,3) - mrSges(4,1) * t529 + m(7) * (t351 * t92 + t573 * t91 + t871) + (-t691 + t690) * mrSges(7,3) + t883) * qJD(3) + t2 * qJD(4) + t4 * qJD(5) + t11 * qJD(6) + t547, t2 * qJD(3) + (-t103 * t695 + t104 * t696 + t629 - t630 + t659 + t672 + t673 + t740 + (t652 - t654) * pkin(9) + t883) * qJD(4) + t6 * qJD(5) + t13 * qJD(6) + ((pkin(9) * t541 - t856) * t799 + (t103 * t572 + t104 * t382 + t870) * t797) * t801 + t546, t4 * qJD(3) + t6 * qJD(4) + (Ifges(6,5) * t665 + Ifges(6,6) * t664 - t333 * t570 - t728 + t729 + (t476 * t96 + t744 * t95) * t794 + t331 * t640 - t185 * mrSges(6,1) - t184 * mrSges(6,2) + t649) * qJD(5) + t19 * qJD(6) + t545, t11 * qJD(3) + t13 * qJD(4) + t19 * qJD(5) + (t649 - t732 - t733) * qJD(6) + t544; t1, t839 * t60, t34 * qJD(4) + t689 + (t509 - t527 + 0.2e1 * t610 * t797 + 0.2e1 * t608 * t799 + (t443 * t634 - t419) * m(5)) * qJD(3) + t880, t34 * qJD(3) + qJD(4) * t509 + t496 * t801 + t688 + t880, -t687 + (m(7) * (t331 * t635 + t333 * t738) - t531 + t186) * qJD(5) + t866 + t839 * t864, qJD(5) * t186 + t549 * t839 + t681 + t866; -qJD(4) * t8 + qJD(5) * t3 + qJD(6) * t12 - t547, qJD(4) * t40 - t689 + t879, qJD(4) * t62 + qJD(5) * t39 + qJD(6) * t51 (m(7) * (t382 * t395 + t394 * t572 + t468 * t632) + (-pkin(4) * t742 + pkin(9) * t526) * t795 + t487) * qJD(4) + t25 * qJD(5) + t36 * qJD(6) + t523, t25 * qJD(4) + ((-t351 * t744 + t476 * t573) * t794 + t510 + t813 * t466 + t887) * qJD(5) + t869 + t524, t36 * qJD(4) + t64 * qJD(5) + t522 + t869; qJD(3) * t8 + qJD(5) * t5 + qJD(6) * t14 - t546, -qJD(3) * t40 - t688 + t879, qJD(5) * t24 + qJD(6) * t35 - t523, qJD(5) * t43 + qJD(6) * t56 ((-t382 * t744 + t476 * t572) * t794 + t510 + t813 * pkin(9) + t886) * qJD(5) + t868 + t516, t74 * qJD(5) + t515 + t868; -qJD(3) * t3 - qJD(4) * t5 + qJD(6) * t20 - t545, t839 * t865 + t687, -qJD(4) * t24 + qJD(6) * t68 - t524, qJD(6) * t88 - t516, -t439, -t439 - t504; -qJD(3) * t12 - qJD(4) * t14 - qJD(5) * t20 - t544, -t548 * t839 - t681, -qJD(4) * t35 - qJD(5) * t68 - t522, -qJD(5) * t88 - t515, t504, 0;];
Cq  = t22;
