% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRRPRR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:11:05
% EndTime: 2019-03-09 18:11:38
% DurationCPUTime: 20.00s
% Computational Cost: add. (33580->733), mult. (62854->906), div. (0->0), fcn. (69811->8), ass. (0->431)
t460 = sin(qJ(5));
t760 = sin(qJ(2));
t658 = t760 * pkin(7);
t527 = -pkin(8) * t760 - t658;
t763 = cos(qJ(2));
t661 = t763 * pkin(7);
t528 = pkin(8) * t763 + t661;
t759 = sin(qJ(3));
t762 = cos(qJ(3));
t336 = -t527 * t762 + t528 * t759;
t422 = -t759 * t763 - t760 * t762;
t490 = pkin(9) * t422 + t336;
t761 = cos(qJ(5));
t421 = t759 * t760 - t762 * t763;
t804 = t527 * t759 + t528 * t762;
t847 = pkin(9) * t421 + t804;
t174 = t460 * t490 + t761 * t847;
t300 = t421 * t460 - t422 * t761;
t459 = sin(qJ(6));
t449 = t459 * mrSges(7,1);
t461 = cos(qJ(6));
t451 = t461 * mrSges(7,2);
t819 = t451 + t449;
t203 = t819 * t300;
t535 = -t421 * t761 - t422 * t460;
t828 = t535 * t459;
t205 = mrSges(7,2) * t300 - mrSges(7,3) * t828;
t830 = t461 * t535;
t207 = mrSges(7,1) * t300 + mrSges(7,3) * t830;
t384 = t422 * qJ(4);
t447 = -pkin(2) * t763 - pkin(1);
t293 = t421 * pkin(3) + t384 + t447;
t234 = -pkin(4) * t421 - t293;
t105 = pkin(5) * t535 - pkin(10) * t300 + t234;
t62 = t105 * t461 - t174 * t459;
t63 = t105 * t459 + t174 * t461;
t823 = t819 * t535;
t897 = t460 * t847 - t490 * t761;
t915 = -t174 * t203 + t63 * t205 - t62 * t207 + t897 * t823;
t899 = m(6) * t234 + mrSges(6,1) * t535 + mrSges(6,2) * t300;
t748 = Ifges(7,4) * t461;
t555 = -Ifges(7,2) * t459 + t748;
t744 = Ifges(7,6) * t300;
t117 = -t535 * t555 + t744;
t496 = t461 * t117;
t435 = Ifges(7,1) * t459 + t748;
t671 = t461 * t435;
t749 = Ifges(7,4) * t459;
t434 = Ifges(7,2) * t461 + t749;
t678 = t459 * t434;
t536 = t678 / 0.2e1 - t671 / 0.2e1;
t710 = t461 * Ifges(7,6);
t713 = t459 * Ifges(7,5);
t433 = t710 + t713;
t692 = t300 * t433;
t747 = Ifges(6,5) * t535;
t556 = Ifges(7,1) * t461 - t749;
t746 = Ifges(7,5) * t300;
t119 = -t535 * t556 + t746;
t497 = t459 * t119;
t745 = Ifges(6,6) * t300;
t898 = t497 / 0.2e1 - t745;
t482 = t496 / 0.2e1 + t692 / 0.2e1 - t747 + t536 * t535 + (-Ifges(5,6) + Ifges(4,6)) * t422 + (-Ifges(4,5) - Ifges(5,4)) * t421 + t898;
t857 = t804 * mrSges(5,1);
t858 = t804 * mrSges(4,1);
t881 = t336 * mrSges(5,3);
t882 = t336 * mrSges(4,2);
t901 = t897 * mrSges(6,2);
t711 = t461 * mrSges(7,1);
t752 = mrSges(7,2) * t459;
t558 = t711 - t752;
t907 = t174 * t558;
t908 = t174 * mrSges(6,1);
t892 = t901 - t907 - t908;
t913 = t482 - t857 - t858 + t882 - t881 + t892;
t912 = t907 / 0.2e1;
t677 = t460 * t897;
t798 = m(6) / 0.2e1;
t905 = t761 * t174;
t890 = (t905 + t677) * t798;
t454 = Ifges(7,5) * t461;
t743 = Ifges(7,6) * t459;
t554 = -t743 + t454;
t793 = Ifges(7,3) / 0.2e1;
t640 = t793 + Ifges(6,2) / 0.2e1;
t764 = t461 / 0.2e1;
t766 = -t459 / 0.2e1;
t841 = t535 / 0.2e1;
t860 = Ifges(6,4) * t300;
t734 = t535 * Ifges(7,5);
t120 = t300 * t556 + t734;
t675 = t461 * t120;
t733 = t535 * Ifges(7,6);
t118 = t300 * t555 + t733;
t681 = t459 * t118;
t776 = t300 / 0.2e1;
t866 = -t234 * mrSges(6,2) - Ifges(6,1) * t776 + Ifges(6,4) * t535 + t681 / 0.2e1 - t675 / 0.2e1;
t839 = -t860 / 0.2e1;
t867 = -t234 * mrSges(6,1) - t554 * t776 - t839 - (Ifges(6,2) + Ifges(7,3)) * t841;
t910 = t293 * (-mrSges(5,1) * t422 + mrSges(5,3) * t421) + t447 * (-mrSges(4,1) * t422 - mrSges(4,2) * t421) + (t554 * t841 - t866) * t535 + (-t535 * t640 + (Ifges(7,1) * t830 - t746) * t764 + (-Ifges(7,2) * t828 - t744) * t766 + Ifges(6,1) * t841 + t860 / 0.2e1 - t830 * t749 + t867) * t300 + t915;
t909 = pkin(5) * t174;
t732 = t535 * mrSges(6,3);
t660 = t762 * pkin(2);
t446 = -t660 - pkin(3);
t440 = -pkin(4) + t446;
t657 = t759 * pkin(2);
t442 = t657 + qJ(4);
t352 = t440 * t761 - t442 * t460;
t350 = pkin(5) - t352;
t879 = t174 * t350;
t462 = -pkin(3) - pkin(4);
t429 = -qJ(4) * t460 + t462 * t761;
t427 = pkin(5) - t429;
t878 = t174 * t427;
t906 = t174 * t897;
t833 = t350 * t819;
t863 = -t833 / 0.2e1;
t832 = t427 * t819;
t864 = -t832 / 0.2e1;
t829 = t461 * t555;
t831 = t459 * t556;
t872 = t536 - t829 / 0.2e1 - t831 / 0.2e1;
t904 = t872 - t863 - t864;
t353 = t440 * t460 + t442 * t761;
t854 = t353 * t897;
t903 = t174 * t352 + t854;
t430 = qJ(4) * t761 + t460 * t462;
t853 = t430 * t897;
t902 = t174 * t429 + t853;
t852 = t459 * t897;
t851 = t461 * t897;
t679 = t207 * t459;
t894 = t461 * t205;
t900 = -t679 - t894;
t846 = m(5) * t293 + mrSges(5,1) * t421 + mrSges(5,3) * t422;
t457 = t459 ^ 2;
t458 = t461 ^ 2;
t665 = t457 + t458;
t419 = t460 * t558;
t525 = t665 * t761;
t848 = -mrSges(6,1) * t460 - mrSges(6,2) * t761 + mrSges(7,3) * t525;
t893 = t848 - t419;
t765 = -t460 / 0.2e1;
t602 = t459 * t765;
t784 = -t823 / 0.2e1;
t889 = -t207 * t602 + t761 * t784;
t633 = -t901 / 0.2e1;
t888 = t633 + t912 + t908 / 0.2e1;
t887 = t857 / 0.2e1 + t858 / 0.2e1 + t881 / 0.2e1 - t882 / 0.2e1;
t795 = -pkin(10) / 0.2e1;
t859 = t300 * mrSges(6,3);
t885 = t859 / 0.2e1;
t756 = m(7) * (-pkin(5) * t460 + pkin(10) * t525);
t771 = -t419 / 0.2e1;
t539 = t756 / 0.2e1 + t771;
t871 = t539 + t771 + t848;
t487 = -mrSges(5,3) + t893;
t623 = t461 * t761;
t565 = -t623 / 0.2e1;
t624 = t459 * t761;
t566 = -t624 / 0.2e1;
t501 = mrSges(7,1) * t566 + mrSges(7,2) * t565;
t822 = t761 * t819;
t844 = t822 / 0.2e1 + t501;
t870 = qJD(6) * t844;
t307 = t844 * qJD(4);
t845 = -t822 / 0.2e1 + t501;
t305 = t845 * qJD(6);
t865 = -pkin(3) * t804 - qJ(4) * t336;
t796 = m(7) / 0.2e1;
t781 = -t535 / 0.2e1;
t861 = -Ifges(5,5) + Ifges(4,4);
t601 = t894 / 0.2e1;
t849 = t679 / 0.2e1 + t601;
t843 = t901 / 0.2e1 - t908 / 0.2e1;
t820 = -t831 / 0.4e1 - t829 / 0.4e1;
t214 = pkin(5) * t300 + pkin(10) * t535;
t589 = t897 * t819;
t840 = t589 / 0.2e1;
t592 = pkin(5) * t819;
t639 = mrSges(6,3) * t761;
t825 = t639 * t781;
t824 = (t454 / 0.2e1 - t743 / 0.2e1) * t535;
t386 = (t460 * t759 + t761 * t762) * pkin(2);
t583 = t665 * t386;
t310 = t353 * t558;
t348 = t352 * mrSges(6,2);
t584 = t665 * t352;
t721 = t353 * mrSges(6,1);
t821 = mrSges(7,3) * t584 - t310 - t348 - t721;
t317 = -pkin(3) * t422 + qJ(4) * t421;
t753 = t422 * pkin(4);
t267 = -t317 + t753;
t115 = -t214 + t267;
t70 = t115 * t461 - t852;
t71 = t115 * t459 + t851;
t551 = -t70 * t459 + t71 * t461;
t818 = (t460 * t551 + t905) * t796 + t890;
t816 = t430 * t203 / 0.2e1 + t427 * t784;
t773 = t353 / 0.2e1;
t815 = t203 * t773 + t350 * t784;
t791 = -t71 / 0.2e1;
t792 = t70 / 0.2e1;
t814 = mrSges(7,1) * t792 + mrSges(7,2) * t791;
t659 = t760 * pkin(2);
t106 = -t659 + t115;
t65 = t106 * t461 - t852;
t66 = t106 * t459 + t851;
t813 = -t66 * mrSges(7,2) / 0.2e1 + t65 * mrSges(7,1) / 0.2e1;
t87 = t214 * t461 + t852;
t88 = t214 * t459 - t851;
t550 = -t87 * t459 + t88 * t461;
t709 = t461 * t66;
t712 = t459 * t65;
t552 = t709 - t712;
t800 = m(5) / 0.2e1;
t811 = (t460 * t552 + t905) * t796 + t890 + t804 * t800;
t687 = t300 * t461;
t209 = mrSges(7,1) * t535 - mrSges(7,3) * t687;
t673 = t461 * t209;
t751 = mrSges(7,3) * t300;
t810 = -t673 / 0.2e1 - t665 * t751 / 0.2e1;
t506 = -t62 * t624 + t623 * t63;
t809 = (t506 + t677) * t796 + t890;
t522 = (-t711 / 0.2e1 + t752 / 0.2e1) * t422;
t593 = t458 / 0.2e1 + t457 / 0.2e1;
t541 = t593 * t751;
t200 = t558 * t300;
t619 = t761 * t200;
t688 = t300 * t459;
t206 = -mrSges(7,2) * t535 - mrSges(7,3) * t688;
t680 = t459 * t206;
t28 = t619 / 0.2e1 + t522 + (t680 / 0.2e1 + t673 / 0.2e1 + t541) * t460;
t806 = (qJD(2) + qJD(3)) * t845 + t28 * qJD(1);
t628 = mrSges(7,3) * t764;
t714 = t459 * mrSges(7,3);
t630 = -t714 / 0.2e1;
t803 = t628 * t66 + t630 * t65 - t843;
t629 = t714 / 0.2e1;
t663 = mrSges(7,3) * t709;
t802 = t65 * t629 - t663 / 0.2e1 + t843;
t607 = -t680 / 0.2e1;
t801 = t607 + t810;
t799 = -m(6) / 0.2e1;
t797 = -m(7) / 0.2e1;
t794 = m(5) * pkin(2);
t790 = t87 / 0.2e1;
t789 = -t88 / 0.2e1;
t788 = -t117 / 0.4e1;
t787 = -t118 / 0.4e1;
t786 = -t119 / 0.4e1;
t785 = t120 / 0.4e1;
t783 = t206 / 0.2e1;
t782 = -t207 / 0.2e1;
t780 = -t535 / 0.4e1;
t778 = t535 / 0.4e1;
t351 = -pkin(10) + t353;
t775 = t351 / 0.2e1;
t774 = -t352 / 0.2e1;
t385 = t460 * t660 - t657 * t761;
t772 = t385 / 0.2e1;
t428 = -pkin(10) + t430;
t770 = t428 / 0.2e1;
t769 = -t429 / 0.2e1;
t768 = -t434 / 0.4e1;
t767 = t435 / 0.4e1;
t758 = m(5) * t804;
t757 = m(7) * (-t761 + t525) * t460;
t755 = pkin(5) * t823;
t296 = t659 + t317;
t508 = (Ifges(4,1) + Ifges(5,1) - Ifges(4,2) - Ifges(5,3)) * t422 + t861 * t421;
t591 = t861 * t422;
t1 = (-mrSges(4,2) * t659 - t591) * t422 + (mrSges(4,1) * t659 + t508) * t421 + m(7) * (t62 * t65 + t63 * t66 - t906) + m(4) * t447 * t659 - pkin(1) * (mrSges(3,1) * t760 + mrSges(3,2) * t763) + t66 * t206 + t65 * t209 + (-Ifges(3,2) + Ifges(3,1)) * t763 * t760 + (-t760 ^ 2 + t763 ^ 2) * Ifges(3,4) + t846 * t296 + t910 + t899 * (-t296 + t753);
t742 = t1 * qJD(1);
t3 = -t591 * t422 + t508 * t421 + m(7) * (t62 * t70 + t63 * t71 - t906) + t71 * t206 + t70 * t209 + t846 * t317 + t899 * t267 + t910;
t731 = t3 * qJD(1);
t720 = t385 * mrSges(6,1);
t719 = t421 * mrSges(5,2);
t718 = t422 * mrSges(5,2);
t717 = t429 * mrSges(6,2);
t716 = t457 * mrSges(7,3);
t715 = t458 * mrSges(7,3);
t27 = (-t673 - t680 + m(7) * (-t459 * t63 - t461 * t62) + t846 - t899) * t422;
t704 = qJD(1) * t27;
t553 = t459 * t62 - t461 * t63;
t10 = t897 * t200 + t62 * t206 - t63 * t209 + (t433 * t781 + t120 * t766 - t461 * t118 / 0.2e1 + t536 * t300 + t553 * mrSges(7,3)) * t300;
t703 = t10 * qJD(1);
t698 = t897 * t385;
t686 = t350 * t200;
t685 = t350 * t823;
t339 = t385 * t558;
t684 = t427 * t200;
t683 = t427 * t823;
t682 = t430 * t558;
t676 = t460 * t203;
t311 = qJD(4) * t757;
t670 = t311 + t870;
t656 = mrSges(5,2) * t384;
t655 = t352 * t732;
t654 = t429 * t732;
t653 = t430 * t859;
t652 = t442 * t718;
t650 = t205 * t795;
t649 = t207 * t795;
t648 = t209 * t795;
t643 = mrSges(7,3) * t790;
t642 = mrSges(7,3) * t789;
t641 = -Ifges(7,2) / 0.4e1 + Ifges(7,1) / 0.4e1;
t638 = t428 * t679;
t637 = t428 * t894;
t634 = t300 * t793;
t378 = -t719 / 0.2e1;
t627 = -t433 / 0.4e1 + Ifges(6,6) / 0.2e1;
t626 = t427 * t460 + t428 * t525;
t625 = t353 * t761;
t387 = t430 * t761;
t617 = t761 * t385;
t614 = t535 * t768;
t613 = t535 * t767;
t612 = t454 * t778;
t600 = t206 * t764;
t582 = t665 * t428;
t581 = t665 * t429;
t580 = t665 * t460;
t577 = t421 * t660;
t576 = -t657 / 0.2e1;
t573 = t859 * t765;
t572 = t209 * t566 + t761 * t600 + t676 / 0.2e1 + t460 * t885;
t564 = -t592 / 0.2e1;
t563 = t627 * t300;
t562 = t435 / 0.2e1 + t748 / 0.2e1;
t561 = t593 * t352;
t478 = t628 * t71 + t630 * t70 + t888;
t463 = (t698 + t903) * t798 + t652 / 0.2e1 + t446 * t378 - t655 / 0.2e1 - mrSges(5,2) * t577 / 0.2e1 + t576 * t718 + t203 * t772 + t685 / 0.2e1 + ((t657 - t442) * t336 + (t446 + t660) * t804) * t800 + (t698 - t879) * t796 - t478 + (t551 * t796 + t849) * t351 + (mrSges(6,3) * t781 + t174 * t798 + t209 * t766 - t553 * t796 + t600) * t386 - t887 + (t385 + t353) * t885;
t465 = -m(5) * t865 / 0.2e1 + t902 * t799 + (t428 * t552 - t878) * t797 + t912 - t683 / 0.2e1 + pkin(3) * t378 - t656 / 0.2e1 + t654 / 0.2e1 - t653 / 0.2e1 - t638 / 0.2e1 - t637 / 0.2e1 + t887;
t4 = t465 + t463 + t803;
t375 = t386 * mrSges(6,2);
t472 = t339 + t375 + t720 + (-mrSges(5,1) - mrSges(4,1)) * t657 + (-t715 - t716) * t386 + (-mrSges(4,2) + mrSges(5,3)) * t660;
t48 = -m(7) * (t350 * t385 + t351 * t583) - m(6) * (-t352 * t385 + t353 * t386) - (t442 * t762 + t446 * t759) * t794 - t472;
t549 = qJD(1) * t4 - qJD(2) * t48;
t47 = -m(7) * (t350 * t353 + t351 * t584) + t821;
t467 = t563 - (-t713 / 0.4e1 - t710 / 0.4e1 + t627) * t300 + t820 * t535 + t755 / 0.2e1;
t481 = (pkin(10) * t552 + t909) * t797 - t907 / 0.2e1;
t499 = t209 * t774 + t351 * t782 + t643 + t786;
t500 = -t205 * t775 + t352 * t783 + t642 + t788;
t545 = t854 + t879;
t6 = t633 + (t558 / 0.2e1 + mrSges(6,1) / 0.2e1) * t174 + (t650 + t500) * t461 + (t649 + t499) * t459 + (t351 * t550 - t352 * t553 + t545) * t796 + t467 + t481 + t802 + t815;
t548 = qJD(1) * t6 - qJD(2) * t47;
t22 = t573 + (-t905 + (t897 + t550) * t460 + t506) * t796 + t894 * t765 + t572 - t889;
t9 = t87 * t209 + t88 * t206 + m(7) * (t62 * t87 + t63 * t88 + t906) + (t117 * t766 + t119 * t764 + t839 - t867) * t300 + (-t824 + (-Ifges(6,1) / 0.2e1 + t640) * t300 + t866) * t535 - t915;
t547 = qJD(1) * t9 + qJD(4) * t22;
t476 = t825 + (t601 + t885) * t460 + t889;
t466 = -t676 / 0.2e1 + t209 * t624 / 0.2e1 + t573 + t206 * t565 + t639 * t841 + t476 - t809;
t17 = t466 - t758 / 0.2e1 + t811;
t509 = t525 * t351;
t491 = t350 * t460 + t509;
t80 = -m(7) * t491 - m(6) * (-t352 * t460 + t625) - m(5) * t442 + t487;
t546 = qJD(1) * t17 + qJD(2) * t80;
t544 = t853 + t878;
t543 = t593 * pkin(10) * mrSges(7,3);
t530 = t461 * t641 + t768;
t489 = (t530 - t749) * t300 + t785;
t479 = Ifges(7,5) * t841 + t489;
t529 = -t435 / 0.4e1 - t641 * t459;
t480 = (t781 + t780) * Ifges(7,6) + t529 * t300 + t787;
t498 = -t634 + t840 + t612;
t11 = -t686 / 0.2e1 + t351 * t541 + (t206 * t775 + t480) * t459 + (t209 * t775 + t479) * t461 + t498 + t813;
t175 = -t833 - t872;
t542 = -qJD(1) * t11 + qJD(2) * t175;
t540 = -t87 * mrSges(7,1) / 0.2e1 + t88 * mrSges(7,2) / 0.2e1;
t537 = -t451 / 0.2e1 - t449 / 0.2e1;
t526 = t419 / 0.2e1 + t539;
t523 = t537 * t386;
t521 = -pkin(5) * t200 / 0.2e1 + t840;
t380 = t429 * t716;
t381 = t429 * t715;
t420 = t430 * mrSges(6,1);
t470 = t310 / 0.2e1 + t348 / 0.2e1 - t380 / 0.2e1 - t381 / 0.2e1 + t420 / 0.2e1 + (t350 * t430 + t351 * t581 + t352 * t582 + t353 * t427) * t796 + t717 / 0.2e1;
t495 = t375 / 0.2e1 + (-pkin(5) * t385 + pkin(10) * t583) * t797;
t31 = -(-t430 / 0.2e1 - t385 / 0.2e1) * t558 + (t773 + t772) * mrSges(6,1) + (-t386 * t593 - t561) * mrSges(7,3) + t470 + t495;
t513 = t380 + t381 - t420 - t682 - t717;
t67 = m(7) * (t427 * t430 + t428 * t581) - t513;
t502 = m(7) * (pkin(10) * t551 + t909);
t517 = t209 * t769 + t428 * t782 + t786;
t518 = -t205 * t770 + t429 * t783 + t788;
t8 = (t650 + t518) * t461 + (t649 + t517) * t459 + (t428 * t550 - t429 * t553 + t544) * t796 - t502 / 0.2e1 + ((t789 + t791) * t461 + (t790 + t792) * t459) * mrSges(7,3) + t467 + t816;
t520 = qJD(1) * t8 + qJD(2) * t31 + qJD(3) * t67;
t101 = m(5) * qJ(4) + m(7) * t626 + m(6) * (-t429 * t460 + t387) - t487;
t20 = t466 + t818;
t468 = m(5) * (0.2e1 * t657 + 0.4e1 * qJ(4)) / 0.4e1 + (t625 + t387 + (-t352 - t429) * t460) * t798 + (t491 + t626) * t796 - t487;
t473 = (t386 * t460 - t617) * t799 + (t386 * t580 - t617) * t797 + m(5) * t576;
t46 = t468 + t473;
t515 = qJD(1) * t20 - qJD(2) * t46 - qJD(3) * t101;
t13 = -t684 / 0.2e1 + t428 * t541 + (t206 * t770 + t480) * t459 + (t209 * t770 + t479) * t461 + t498 + t814;
t198 = -t832 - t872;
t90 = t523 + t904;
t514 = -qJD(1) * t13 - qJD(2) * t90 + qJD(3) * t198;
t507 = t592 / 0.2e1 + t536 + 0.2e1 * t820;
t503 = (-Ifges(7,2) / 0.2e1 + Ifges(7,1) / 0.2e1) * t461 - t434 / 0.2e1 - t749 / 0.2e1;
t475 = m(7) * (-t625 + t509 + (t350 + t584) * t460);
t49 = -t475 / 0.2e1 + t871;
t477 = (t429 * t580 - t387 + t626) * t796;
t64 = -t477 + t871;
t493 = qJD(1) * t22 - qJD(2) * t49 - qJD(3) * t64 + t311;
t16 = t612 + (-Ifges(7,3) / 0.2e1 - t543) * t300 + (t648 + t785 + t734 / 0.2e1 + t530 * t300) * t461 + (-0.3e1 / 0.4e1 * t733 + t206 * t795 + t787 + (t529 - t748) * t300) * t459 + t521 + t540;
t215 = t592 + t872;
t73 = t863 + t564 + (mrSges(7,2) * t774 + t562) * t461 + (mrSges(7,1) * t774 + t503) * t459;
t81 = t564 + t864 + (mrSges(7,2) * t769 + t562) * t461 + (mrSges(7,1) * t769 + t503) * t459;
t484 = qJD(1) * t16 - qJD(2) * t73 - qJD(3) * t81 - qJD(5) * t215 - t307;
t474 = t824 + t681 / 0.4e1 - t675 / 0.4e1 - t589 / 0.2e1 - t554 * t778 - t634 + (0.2e1 * t767 + t555 / 0.4e1) * t688 + (t434 / 0.2e1 - t556 / 0.4e1) * t687 + (t629 + t630) * t63;
t471 = -t497 / 0.4e1 - t496 / 0.4e1 - t755 / 0.2e1 - t692 / 0.4e1 + t745 / 0.2e1 + t747 / 0.2e1 + t563 + Ifges(6,5) * t841 + (-t678 / 0.4e1 + t671 / 0.4e1) * t535 + t849 * pkin(10) + t888;
t469 = 0.2e1 * t378 + t476 + t572 + t809 + t825;
t304 = t845 * qJD(4);
t91 = t523 - t904;
t82 = t832 / 0.2e1 + t537 * t429 + t507;
t75 = t477 + t526;
t74 = t833 / 0.2e1 + t537 * t352 + t507;
t50 = t475 / 0.2e1 + t526;
t45 = t468 - t473;
t30 = t721 / 0.2e1 + t682 / 0.2e1 - t339 / 0.2e1 - t720 / 0.2e1 + t470 - t495 + (-t561 + t583 / 0.2e1) * mrSges(7,3);
t29 = t206 * t602 - t619 / 0.2e1 + t522 + t810 * t460;
t21 = t22 * qJD(5);
t19 = t469 + t758 + t818;
t18 = t469 + t758 / 0.2e1 + t811;
t15 = -t554 * t780 + pkin(10) * t607 - t681 / 0.4e1 - Ifges(7,5) * t830 / 0.2e1 + Ifges(7,6) * t828 / 0.2e1 + Ifges(7,3) * t776 + (t459 * t529 - t543) * t300 + (t648 + t489) * t461 + t521 - t540;
t14 = t684 / 0.2e1 + t474 + t801 * t428 + t814;
t12 = t686 / 0.2e1 + t474 + t801 * t351 + t813;
t7 = t471 + t502 / 0.2e1 + t544 * t796 + ((t428 * t88 + t429 * t63) * t796 + t613 + t642 + t518) * t461 + ((-t428 * t87 - t429 * t62) * t796 + t614 + t643 + t517) * t459 + t478 + t816;
t5 = t471 + t545 * t796 + ((t351 * t88 + t352 * t63) * t796 + t613 + t500) * t461 + ((-t351 * t87 - t352 * t62) * t796 + t614 + t499) * t459 - t481 + t803 + t815;
t2 = t463 - t465 + t482 + t802;
t23 = [qJD(2) * t1 + qJD(3) * t3 + qJD(4) * t27 + qJD(5) * t9 + qJD(6) * t10, t742 + (m(4) * (-t336 * t759 - t762 * t804) * pkin(2) + t685 - t663 - t655 + t652 + mrSges(3,2) * t658 + Ifges(3,5) * t763 - Ifges(3,6) * t760 - t446 * t719 - mrSges(3,1) * t661 + (m(7) * t552 - t900) * t351 + mrSges(7,3) * t712 - m(7) * t879 + m(6) * t903 + (t422 * t657 + t577) * mrSges(4,3) + m(5) * (-t336 * t442 + t446 * t804) + t353 * t859 + t913) * qJD(2) + t2 * qJD(3) + t18 * qJD(4) + t5 * qJD(5) + t12 * qJD(6), t731 + t2 * qJD(2) + (-mrSges(7,3) * t551 + pkin(3) * t719 + t637 + t638 + t653 - t654 + t656 + t683 + t913) * qJD(3) + t19 * qJD(4) + t7 * qJD(5) + t14 * qJD(6) + 0.2e1 * ((t428 * t551 - t878) * t796 + t902 * t798 + t865 * t800) * qJD(3), qJD(2) * t18 + qJD(3) * t19 + qJD(6) * t29 + t21 + t704, t5 * qJD(2) + t7 * qJD(3) + t15 * qJD(6) + t547 + (t433 * t776 + t117 * t764 + (-Ifges(6,5) + t536) * t535 + (-m(7) * t174 + t823) * pkin(5) + (m(7) * t550 + t900) * pkin(10) + t550 * mrSges(7,3) + t892 + t898) * qJD(5), t703 + t12 * qJD(2) + t14 * qJD(3) + t29 * qJD(4) + t15 * qJD(5) + (-mrSges(7,1) * t63 - mrSges(7,2) * t62 - t692) * qJD(6); qJD(3) * t4 - qJD(4) * t17 + qJD(5) * t6 - qJD(6) * t11 - t742, -qJD(3) * t48 - qJD(4) * t80 - qJD(5) * t47 + qJD(6) * t175 (m(7) * (t385 * t427 + t386 * t582) + m(6) * (-t385 * t429 + t386 * t430) + (-pkin(3) * t759 + qJ(4) * t762) * t794 + t472) * qJD(3) + t45 * qJD(4) + t30 * qJD(5) + t91 * qJD(6) + t549, qJD(3) * t45 + qJD(5) * t50 - t546 + t670, t30 * qJD(3) + t50 * qJD(4) + (m(7) * (-pkin(5) * t353 + pkin(10) * t584) + t821) * qJD(5) + t74 * qJD(6) + t548, t91 * qJD(3) + t307 + t74 * qJD(5) + (-t351 * t558 - t554) * qJD(6) + t542; -qJD(2) * t4 - qJD(4) * t20 + qJD(5) * t8 - qJD(6) * t13 - t731, qJD(4) * t46 + qJD(5) * t31 - qJD(6) * t90 - t549, qJD(4) * t101 + qJD(5) * t67 + qJD(6) * t198, qJD(5) * t75 - t515 + t670, t75 * qJD(4) + (m(7) * (-pkin(5) * t430 + pkin(10) * t581) + t513) * qJD(5) + t82 * qJD(6) + t520, t307 + t82 * qJD(5) + (-t428 * t558 - t554) * qJD(6) + t514; qJD(2) * t17 + qJD(3) * t20 - qJD(6) * t28 + t21 - t704, -qJD(3) * t46 - qJD(5) * t49 - t305 + t546, -qJD(5) * t64 - t305 + t515, qJD(5) * t757 (t756 + t893) * qJD(5) + t305 + t493, qJD(5) * t845 - qJD(6) * t419 - t806; -qJD(2) * t6 - qJD(3) * t8 + qJD(6) * t16 - t547, -qJD(3) * t31 + qJD(4) * t49 - qJD(6) * t73 - t548, qJD(4) * t64 - qJD(6) * t81 - t520, -t493 - t870, -t215 * qJD(6) (-pkin(10) * t558 + t554) * qJD(6) + t484; qJD(2) * t11 + qJD(3) * t13 + qJD(4) * t28 - qJD(5) * t16 - t703, qJD(3) * t90 + qJD(5) * t73 + t304 - t542, qJD(5) * t81 + t304 - t514, qJD(5) * t844 + t806, -t484, 0;];
Cq  = t23;
