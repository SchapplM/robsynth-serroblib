% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPRRR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:16:52
% EndTime: 2019-03-09 13:17:32
% DurationCPUTime: 25.26s
% Computational Cost: add. (64885->795), mult. (125002->1034), div. (0->0), fcn. (151205->10), ass. (0->427)
t485 = sin(qJ(4));
t482 = sin(pkin(11));
t483 = cos(pkin(11));
t486 = sin(qJ(2));
t722 = -qJ(3) - pkin(7);
t572 = t722 * t486;
t488 = cos(qJ(2));
t841 = t722 * t488;
t411 = -t482 * t572 + t483 * t841;
t448 = -t482 * t486 + t483 * t488;
t524 = -t448 * pkin(8) + t411;
t745 = cos(qJ(4));
t449 = -t482 * t488 - t483 * t486;
t867 = t482 * t841 + t483 * t572;
t892 = t449 * pkin(8) + t867;
t253 = t485 * t524 + t745 * t892;
t484 = sin(qJ(5));
t487 = cos(qJ(5));
t743 = sin(qJ(6));
t744 = cos(qJ(6));
t528 = t484 * t744 + t487 * t743;
t832 = t743 * t484 - t744 * t487;
t401 = mrSges(7,1) * t832 + mrSges(7,2) * t528;
t714 = Ifges(7,4) * t528;
t404 = -Ifges(7,2) * t832 + t714;
t443 = Ifges(7,4) * t832;
t406 = Ifges(7,1) * t528 - t443;
t759 = t528 / 0.2e1;
t761 = -t832 / 0.2e1;
t396 = t745 * t448 + t449 * t485;
t879 = t832 * t396;
t790 = -t879 / 0.2e1;
t875 = t528 * t396;
t792 = -t875 / 0.2e1;
t536 = t485 * t448 - t449 * t745;
t908 = -Ifges(7,4) * t879 - Ifges(7,2) * t875 + Ifges(7,6) * t536;
t909 = -Ifges(7,1) * t879 - Ifges(7,4) * t875 + Ifges(7,5) * t536;
t877 = t484 * t396;
t932 = t485 * t892 - t745 * t524;
t936 = pkin(5) * t877 + t932;
t948 = -t253 * mrSges(5,2) + t936 * t401 + t404 * t792 + t406 * t790 + t759 * t909 + t908 * t761;
t661 = t536 * t484;
t183 = pkin(5) * t661 - t253;
t946 = t183 * t936;
t741 = pkin(2) * t483;
t471 = pkin(3) + t741;
t742 = pkin(2) * t482;
t434 = t485 * t471 + t742 * t745;
t945 = t253 * t434;
t944 = t253 * t484;
t669 = t253 * t932;
t433 = t471 * t745 - t485 * t742;
t430 = -pkin(4) - t433;
t733 = t487 * pkin(5);
t425 = t430 - t733;
t943 = t425 * t936;
t474 = -pkin(4) - t733;
t942 = t474 * t936;
t941 = t487 * t253;
t475 = -pkin(2) * t488 - pkin(1);
t426 = -t448 * pkin(3) + t475;
t250 = -pkin(4) * t396 - pkin(9) * t536 + t426;
t125 = t487 * t250 - t484 * t932;
t126 = t250 * t484 + t487 * t932;
t269 = t832 * t536;
t272 = t528 * t536;
t156 = mrSges(7,1) * t272 - mrSges(7,2) * t269;
t683 = t487 * mrSges(6,2);
t684 = t484 * mrSges(6,1);
t461 = t683 + t684;
t290 = t461 * t536;
t101 = -pkin(10) * t661 + t126;
t602 = t743 * t101;
t660 = t536 * t487;
t100 = -pkin(10) * t660 + t125;
t93 = -pkin(5) * t396 + t100;
t58 = t744 * t93 - t602;
t604 = t744 * t101;
t59 = t743 * t93 + t604;
t791 = -t272 / 0.2e1;
t793 = -t269 / 0.2e1;
t878 = t461 * t396;
t888 = -mrSges(7,2) * t536 - mrSges(7,3) * t875;
t876 = t487 * t396;
t889 = mrSges(6,1) * t536 - mrSges(6,3) * t876;
t890 = -mrSges(6,2) * t536 - mrSges(6,3) * t877;
t910 = mrSges(7,1) * t536 + mrSges(7,3) * t879;
t700 = t879 * mrSges(7,2);
t702 = t875 * mrSges(7,1);
t926 = -t700 + t702;
t939 = t125 * t889 + t126 * t890 + t936 * t156 + t183 * t926 - t253 * t878 + t932 * t290 + t58 * t910 + t59 * t888 + t908 * t791 + t909 * t793;
t938 = pkin(4) * t932;
t937 = t430 * t932;
t934 = t474 * t926;
t740 = pkin(4) * t536;
t310 = -pkin(9) * t396 + t740;
t148 = t487 * t310 - t944;
t930 = -t148 / 0.2e1;
t149 = t484 * t310 + t941;
t929 = t149 / 0.2e1;
t924 = -t449 * mrSges(4,1) + t448 * mrSges(4,2);
t630 = pkin(5) * t744;
t565 = t630 / 0.2e1;
t629 = pkin(5) * t743;
t805 = Ifges(6,3) / 0.2e1;
t923 = -t910 * t565 - t888 * t629 / 0.2e1 - t536 * t805;
t736 = pkin(9) * t889;
t737 = pkin(9) * t890;
t431 = pkin(9) + t434;
t920 = t431 * t890;
t692 = Ifges(6,5) * t396;
t609 = t692 / 0.2e1;
t711 = Ifges(6,6) * t484;
t613 = t711 / 0.2e1;
t806 = m(7) * pkin(5);
t633 = t806 / 0.2e1;
t898 = pkin(10) * t877;
t110 = t149 - t898;
t891 = pkin(5) * t536 - pkin(10) * t876;
t98 = t148 + t891;
t73 = -t110 * t743 + t744 * t98;
t74 = t110 * t744 + t743 * t98;
t540 = Ifges(7,5) * t790 + Ifges(7,6) * t792;
t851 = Ifges(7,3) * t536;
t903 = t73 / 0.2e1;
t905 = -mrSges(7,2) / 0.2e1;
t820 = t851 / 0.2e1 + mrSges(7,1) * t903 + t74 * t905 + t540;
t916 = mrSges(6,1) * t930 + mrSges(6,2) * t929 - (t73 * t744 + t74 * t743) * t633 - t487 * t609 + t396 * t613 - t820 + t923;
t344 = t528 * t433;
t345 = t832 * t433;
t539 = t344 * mrSges(7,1) / 0.2e1 + t345 * t905;
t605 = -t683 / 0.2e1;
t763 = -t433 / 0.2e1;
t915 = (-t344 * t744 - t345 * t743) * t633 + t684 * t763 + t433 * t605 - t539;
t585 = -t877 / 0.2e1;
t479 = t486 * pkin(2);
t427 = -pkin(3) * t449 + t479;
t255 = t310 + t427;
t130 = t484 * t255 + t941;
t102 = t130 - t898;
t129 = t487 * t255 - t944;
t94 = t129 + t891;
t60 = -t102 * t743 + t744 * t94;
t61 = t102 * t744 + t743 * t94;
t799 = -t130 / 0.2e1;
t800 = t129 / 0.2e1;
t710 = Ifges(7,6) * t875;
t713 = Ifges(7,5) * t879;
t904 = -t60 / 0.2e1;
t818 = -t851 / 0.2e1 + t710 / 0.2e1 + t713 / 0.2e1 + t61 * mrSges(7,2) / 0.2e1 + mrSges(7,1) * t904;
t914 = t818 - mrSges(6,1) * t800 - mrSges(6,2) * t799 - (t60 * t744 + t61 * t743) * t633 - Ifges(6,5) * t876 / 0.2e1 - Ifges(6,6) * t585 + t923;
t747 = t487 / 0.2e1;
t751 = -t484 / 0.2e1;
t476 = Ifges(6,5) * t487;
t837 = t476 - t711;
t864 = t536 / 0.2e1;
t865 = -t536 / 0.2e1;
t477 = Ifges(6,4) * t487;
t836 = -Ifges(6,2) * t484 + t477;
t868 = Ifges(6,6) * t536 + t396 * t836;
t716 = Ifges(6,4) * t484;
t466 = Ifges(6,1) * t487 - t716;
t869 = Ifges(6,5) * t536 + t396 * t466;
t884 = -t396 / 0.2e1;
t913 = 0.2e1 * Ifges(5,4) * t865 + t869 * t747 + t868 * t751 + Ifges(7,5) * t793 + Ifges(7,6) * t791 + t837 * t864 + (Ifges(5,2) + Ifges(6,3) + Ifges(7,3)) * t884;
t804 = -pkin(10) - pkin(9);
t467 = t804 * t484;
t468 = t804 * t487;
t418 = -t467 * t743 + t468 * t744;
t550 = t744 * t467 + t468 * t743;
t912 = t418 * mrSges(7,1) - t550 * mrSges(7,2);
t723 = pkin(10) + t431;
t416 = t723 * t484;
t417 = t723 * t487;
t322 = t416 * t743 - t417 * t744;
t551 = -t744 * t416 - t417 * t743;
t911 = t322 * mrSges(7,1) - t551 * mrSges(7,2);
t902 = t868 / 0.2e1;
t901 = t869 / 0.2e1;
t900 = pkin(4) * t878;
t70 = t100 * t744 - t602;
t897 = t58 - t70;
t480 = t484 ^ 2;
t481 = t487 ^ 2;
t635 = t480 + t481;
t882 = mrSges(6,3) * t635;
t636 = -Ifges(7,5) * t832 - Ifges(7,6) * t528;
t75 = t636 + t911;
t895 = t75 * qJD(6);
t99 = t636 + t912;
t893 = t99 * qJD(6);
t785 = -t322 / 0.2e1;
t781 = -t396 / 0.4e1;
t883 = t396 / 0.2e1;
t767 = -t418 / 0.2e1;
t758 = -t528 / 0.2e1;
t69 = -t100 * t743 - t604;
t860 = t59 + t69;
t400 = mrSges(7,1) * t528 - mrSges(7,2) * t832;
t872 = t400 * qJD(6);
t405 = -Ifges(7,1) * t832 - t714;
t577 = t404 / 0.4e1 - t405 / 0.4e1;
t403 = -Ifges(7,2) * t528 - t443;
t578 = t403 / 0.4e1 + t406 / 0.4e1;
t794 = -t253 / 0.2e1;
t871 = t577 * t269 - t578 * t272 + t461 * t794;
t153 = -mrSges(7,1) * t269 - mrSges(7,2) * t272;
t701 = t272 * mrSges(7,3);
t195 = mrSges(7,2) * t396 - t701;
t766 = t425 / 0.2e1;
t863 = t551 / 0.2e1;
t870 = t153 * t766 + t195 * t863;
t685 = t528 * mrSges(7,3);
t831 = -mrSges(6,1) * t487 + mrSges(6,2) * t484;
t848 = t831 - mrSges(5,1);
t402 = Ifges(7,5) * t528 - Ifges(7,6) * t832;
t462 = Ifges(6,5) * t484 + Ifges(6,6) * t487;
t830 = t402 / 0.4e1 + t462 / 0.4e1;
t844 = (-Ifges(5,6) / 0.2e1 + t830) * t536;
t465 = Ifges(6,1) * t484 + t477;
t770 = -t550 / 0.2e1;
t834 = t269 * t767 - t272 * t770;
t789 = -t551 / 0.2e1;
t833 = t269 * t785 - t272 * t789;
t463 = Ifges(6,2) * t487 + t716;
t574 = -t463 / 0.4e1 + t466 / 0.4e1;
t829 = t884 + t883;
t212 = t466 * t536 - t692;
t646 = t487 * t212;
t691 = Ifges(6,6) * t396;
t209 = t536 * t836 - t691;
t650 = t484 * t209;
t828 = Ifges(5,4) * t396 + t426 * mrSges(5,2) + Ifges(5,1) * t864 - t650 / 0.2e1 + t646 / 0.2e1;
t827 = (t344 * t528 + t345 * t832) * mrSges(7,3) + (t401 + t848) * t434 + (-mrSges(5,2) + t882) * t433;
t198 = -mrSges(7,1) * t396 + t269 * mrSges(7,3);
t826 = t195 * t761 + t198 * t758;
t640 = -t702 / 0.2e1 + t700 / 0.2e1;
t686 = t832 * mrSges(7,3);
t608 = t686 / 0.2e1;
t560 = -t272 * t608 - t685 * t793 + t826;
t22 = t560 - t640;
t825 = t22 * qJD(1);
t302 = mrSges(6,2) * t396 - mrSges(6,3) * t661;
t305 = -mrSges(6,1) * t396 - mrSges(6,3) * t660;
t748 = -t487 / 0.2e1;
t824 = t302 * t751 + t305 * t748 + t865 * t882;
t287 = t831 * t536;
t752 = t474 / 0.2e1;
t587 = t153 * t752;
t769 = t550 / 0.2e1;
t795 = t198 / 0.2e1;
t809 = -pkin(4) / 0.2e1;
t823 = t195 * t769 - t287 * t809 + t418 * t795 + t587;
t765 = t430 / 0.2e1;
t822 = -t287 * t765 + t322 * t795 + t870;
t715 = Ifges(7,4) * t269;
t115 = -Ifges(7,2) * t272 - Ifges(7,6) * t396 - t715;
t158 = -Ifges(7,1) * t272 + t715;
t817 = t183 * t400 / 0.2e1 + t636 * t781 - (-t158 / 0.4e1 + t115 / 0.4e1) * t528;
t815 = 2 * qJD(4);
t814 = m(5) / 0.2e1;
t813 = -m(6) / 0.2e1;
t812 = m(6) / 0.2e1;
t811 = -m(7) / 0.2e1;
t810 = m(7) / 0.2e1;
t808 = -pkin(9) / 0.2e1;
t807 = m(4) * pkin(2);
t803 = t115 / 0.2e1;
t266 = Ifges(7,4) * t272;
t118 = -Ifges(7,1) * t269 - Ifges(7,5) * t396 - t266;
t801 = t118 / 0.2e1;
t798 = t156 / 0.2e1;
t796 = t195 / 0.2e1;
t786 = t322 / 0.2e1;
t778 = t396 / 0.4e1;
t776 = t401 / 0.2e1;
t768 = t418 / 0.2e1;
t764 = -t431 / 0.2e1;
t762 = t832 / 0.2e1;
t757 = t831 / 0.2e1;
t754 = -t465 / 0.4e1;
t750 = t484 / 0.2e1;
t749 = t484 / 0.4e1;
t746 = t487 / 0.4e1;
t739 = pkin(4) * t461;
t738 = pkin(5) * t484;
t732 = t58 * mrSges(7,2);
t731 = t59 * mrSges(7,1);
t728 = t69 * mrSges(7,1);
t727 = t70 * mrSges(7,2);
t721 = m(7) * qJD(3);
t707 = t148 * mrSges(6,3);
t706 = t149 * mrSges(6,3);
t705 = t932 * mrSges(5,1);
t391 = t536 * mrSges(5,1);
t558 = t805 + Ifges(5,2) / 0.2e1 + Ifges(7,3) / 0.2e1;
t3 = (-mrSges(4,2) * t479 - Ifges(4,4) * t449) * t449 + m(7) * (t58 * t60 + t59 * t61 + t946) + (m(4) * t479 + t924) * t475 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t488 + (Ifges(3,1) - Ifges(3,2)) * t486) * t488 - t879 * t801 + (t427 * mrSges(5,2) + Ifges(5,1) * t883 + t913) * t536 + (-t710 - t713) * t884 + (-t427 * mrSges(5,1) - t558 * t536 + t884 * t837 + t828) * t396 - t875 * t803 + t130 * t302 + t129 * t305 + (-mrSges(4,1) * t479 + Ifges(4,4) * t448 + (-Ifges(4,1) + Ifges(4,2)) * t449) * t448 + t61 * t195 + t60 * t198 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t486) * t486 + m(6) * (t125 * t129 + t126 * t130 - t669) + t939 + (m(5) * t427 + t391) * t426;
t699 = t3 * qJD(1);
t682 = t73 * t528;
t681 = t74 * t832;
t8 = (t426 * mrSges(5,1) - (-Ifges(5,1) / 0.2e1 + t558) * t396 + t913) * t536 + m(6) * (t125 * t148 + t126 * t149 - t669) + m(7) * (t58 * t73 + t59 * t74 + t946) - (-(-t476 / 0.2e1 + t613) * t396 + t540 - t828) * t396 + t118 * t790 + t115 * t792 + t149 * t302 + t148 * t305 + t74 * t195 + t73 * t198 + t939;
t680 = t8 * qJD(1);
t291 = t536 * t463;
t292 = t465 * t536;
t157 = Ifges(7,2) * t269 - t266;
t639 = -Ifges(7,5) * t272 + Ifges(7,6) * t269;
t502 = t183 * t153 - (t801 + t157 / 0.2e1) * t272 + (-t158 / 0.2e1 + t803 + t59 * mrSges(7,3)) * t269 + t58 * t701 + t639 * t884;
t9 = t69 * t198 + t70 * t195 + m(7) * (t58 * t69 + t59 * t70) + t253 * t287 + t125 * t302 - t126 * t305 + ((t609 + t125 * mrSges(6,3) - t212 / 0.2e1 + t291 / 0.2e1) * t484 + (t691 / 0.2e1 - t126 * mrSges(6,3) - t292 / 0.2e1 - t209 / 0.2e1 + (m(7) * t183 + t156) * pkin(5)) * t487) * t536 + t502;
t679 = t9 * qJD(1);
t569 = t635 * t396;
t495 = -m(5) * (t396 * t434 - t433 * t536) / 0.2e1 + (t430 * t536 + t431 * t569) * t813 + (t322 * t879 + t425 * t536 - t551 * t875) * t811 - (t448 * t482 + t449 * t483) * t807 / 0.2e1;
t634 = t807 / 0.2e1;
t497 = t427 * t814 + (t129 * t487 + t130 * t484) * t812 + (t528 * t61 - t60 * t832) * t810 + t910 * t761 + t888 * t759 + t890 * t750 + t889 * t747 + t486 * t634;
t556 = (-t481 / 0.2e1 - t480 / 0.2e1) * mrSges(6,3);
t20 = t391 + (-t401 / 0.2e1 - t831 / 0.2e1) * t536 + (-t759 * t875 - t762 * t879) * mrSges(7,3) + (mrSges(5,2) + t556) * t396 + t495 + t497 + t924;
t678 = qJD(1) * t20;
t677 = t129 * t484;
t676 = t130 * t487;
t14 = t58 * t195 - t59 * t198 + t502;
t675 = t14 * qJD(1);
t503 = (-t743 * t879 - t744 * t875) * t633 + mrSges(6,1) * t585 + t396 * t605 + t640;
t645 = t487 * t302;
t649 = t484 * t305;
t506 = (-t528 * t897 - t860 * t832) * t810 - t649 / 0.2e1 + t645 / 0.2e1;
t606 = -t685 / 0.2e1;
t607 = -t686 / 0.2e1;
t16 = t269 * t606 - t272 * t607 + t503 - t506 - t826;
t674 = t16 * qJD(1);
t545 = -t125 * t484 + t126 * t487;
t667 = t253 * t536;
t18 = -t879 * t195 - t875 * t198 + (t448 ^ 2 + t449 ^ 2) * mrSges(4,3) + (mrSges(5,3) * t396 + t645 - t649) * t396 + (mrSges(5,3) * t536 + t156 + t290) * t536 + m(7) * (t183 * t536 - t58 * t875 - t59 * t879) + m(6) * (t396 * t545 - t667) + m(5) * (t396 * t932 - t667) + m(4) * (-t411 * t448 + t449 * t867);
t673 = t18 * qJD(1);
t668 = t932 * t831;
t659 = t425 * t400;
t658 = t430 * t461;
t651 = t474 * t400;
t647 = t484 * t463;
t642 = t487 * t465;
t632 = pkin(5) * t660;
t192 = t344 * t832 - t345 * t528;
t628 = t192 * t810;
t631 = qJD(4) * t628;
t626 = pkin(5) * t798;
t621 = t69 / 0.2e1 + t59 / 0.2e1;
t620 = t70 / 0.2e1 - t58 / 0.2e1;
t601 = t743 * t198;
t583 = t212 / 0.4e1 - t291 / 0.4e1;
t581 = t789 + t863;
t580 = t786 + t785;
t579 = t778 + t781;
t576 = t770 + t769;
t575 = t768 + t767;
t568 = t635 * t433;
t567 = mrSges(7,3) * t630;
t566 = mrSges(7,3) * t629;
t561 = t462 / 0.2e1 + t402 / 0.2e1 - Ifges(5,6);
t554 = (t766 + t752) * t400;
t553 = t272 * t567;
t552 = t269 * t566;
t544 = t676 - t677;
t543 = -t148 * t484 + t149 * t487;
t542 = t937 - t945;
t541 = t536 * t556;
t538 = -t875 * t606 - t879 * t607 + t883 * t882 + (t757 + t776) * t536;
t537 = t647 / 0.2e1 - t642 / 0.2e1;
t499 = (t148 * t487 + t149 * t484) * t813 + (t528 * t74 - t73 * t832) * t811 + mrSges(5,1) * t865 + t910 * t762 + t888 * t758 + t890 * t751 + t889 * t748;
t500 = -t391 / 0.2e1 + (pkin(9) * t569 - t740) * t812 + (t418 * t879 + t474 * t536 - t550 * t875) * t810 + t538;
t24 = 0.2e1 * t884 * mrSges(5,2) + t499 + t500;
t531 = -t24 * qJD(1) + qJD(2) * t628;
t530 = -Ifges(5,5) + t537;
t526 = t868 / 0.4e1 + t920 / 0.2e1 + t433 * t302 / 0.2e1;
t525 = t869 / 0.4e1 + t889 * t764 + t305 * t763;
t415 = t425 * t738;
t517 = -(t403 / 0.2e1 + t406 / 0.2e1) * t832 - (t404 / 0.2e1 - t405 / 0.2e1) * t528;
t501 = (t465 / 0.2e1 + t836 / 0.2e1) * t487 + (pkin(5) * t401 + t466 / 0.2e1 - t463 / 0.2e1) * t484 + t517;
t44 = m(7) * t415 + t501 + t658 + t659;
t489 = -t650 / 0.4e1 + t646 / 0.4e1 + t632 * t776 + t837 * t781 - t291 * t746 - t292 * t749 + t70 * t607 + t58 * t608 + t484 * t626 - (t157 + t118) * t832 / 0.4e1 + t574 * t660 + t860 * t606 - (t465 + t836) * t661 / 0.4e1 + t817 + t871;
t181 = t183 * t738;
t516 = t322 * t897 + t860 * t551 + t181;
t5 = t489 + (t425 * t632 + t516) * t810 + t833 * mrSges(7,3) + t824 * t431 + t822 + t914;
t523 = -t5 * qJD(1) - t44 * qJD(2);
t505 = -(t118 / 0.4e1 + t157 / 0.4e1) * t832 + t817;
t492 = (mrSges(7,3) * t785 + t577) * t269 - (mrSges(7,3) * t789 + t578) * t272 + t198 * t786 + t505 + t870;
t11 = t492 + t818;
t46 = t517 + t659;
t522 = -t11 * qJD(1) - t46 * qJD(2);
t490 = (t290 / 0.2e1 + t798) * t434 + t844 + (t183 * t434 - t322 * t74 - t344 * t58 - t345 * t59 + t551 * t73 + t943) * t810 + t910 * t863 + t888 * t785 - t344 * t795 - t345 * t796 + t926 * t766 + t878 * t765;
t498 = (pkin(9) * t544 - t938) * t813 + (-t418 * t61 + t550 * t60 + t942) * t811 + t900 / 0.2e1 + t910 * t770 + t888 * t768 - t934 / 0.2e1;
t2 = t829 * Ifges(5,5) + (t794 + t253 / 0.2e1) * mrSges(5,2) + t498 + t490 - t844 + (t431 * t543 + t433 * t545 + t542) * t812 + (t736 / 0.2e1 - t869 / 0.4e1 + t579 * t463 + (t930 + t800) * mrSges(6,3) + t525) * t484 + (-t737 / 0.2e1 - t868 / 0.4e1 - t579 * t465 + (t929 + t799) * mrSges(6,3) + t526) * t487 + (-(t903 + t904) * t528 - (t74 / 0.2e1 - t61 / 0.2e1) * t832) * mrSges(7,3);
t45 = m(7) * (t322 * t345 - t344 * t551 + t425 * t434) + m(6) * (t430 * t434 + t431 * t568) + t827;
t521 = -t2 * qJD(1) - t45 * qJD(2) - t192 * t721 / 0.2e1;
t519 = -t528 * t621 - t620 * t832;
t515 = t418 * t897 + t860 * t550 + t181;
t459 = t474 * t738;
t508 = (t415 + t459) * t810;
t33 = t401 * t738 + t508 + t836 * t747 + t466 * t750 + t651 / 0.2e1 - t739 / 0.2e1 + t658 / 0.2e1 + t404 * t758 + t405 * t759 + t659 / 0.2e1 - t537 + (t406 + t403) * t761 + (t608 + t607) * (t550 + t551) - t915;
t47 = m(7) * t459 + t501 + t651 - t739;
t7 = t489 + (t474 * t632 + t515) * t810 + t834 * mrSges(7,3) + t824 * pkin(9) + t823 + t916;
t514 = -t7 * qJD(1) - t33 * qJD(2) - t47 * qJD(4);
t491 = (mrSges(7,3) * t767 + t577) * t269 - (mrSges(7,3) * t770 + t578) * t272 + t550 * t796 + t198 * t768 + t587 + t505;
t13 = t491 - t820;
t507 = t554 + t517;
t41 = t507 + t539;
t79 = t517 + t651;
t513 = -t13 * qJD(1) - t41 * qJD(2) - t79 * qJD(4);
t512 = -t528 * t566 + t567 * t832 + t636 + t837;
t510 = (t754 - t836 / 0.4e1) * t536 - t209 / 0.4e1 - t292 / 0.4e1 + t626 + t691 / 0.4e1;
t109 = mrSges(7,1) * t575 + mrSges(7,2) * t576;
t19 = -t620 * mrSges(7,2) + t621 * mrSges(7,1) + (-t744 * t195 / 0.2e1 + t601 / 0.2e1 + (t743 * t793 + t744 * t791) * mrSges(7,3)) * pkin(5);
t458 = (mrSges(7,1) * t743 + mrSges(7,2) * t744) * pkin(5);
t78 = mrSges(7,1) * t580 + mrSges(7,2) * t581;
t509 = t19 * qJD(1) - t78 * qJD(2) - t109 * qJD(4) + t458 * qJD(5);
t496 = t476 * t781 + t505 + t871;
t450 = t458 * qJD(6);
t42 = t507 - t539;
t34 = t508 + (t809 + t765) * t461 + t554 + (-(-t575 - t580) * t528 - (t576 + t581) * t832) * mrSges(7,3) + t501 + t915;
t25 = mrSges(5,2) * t829 - t499 + t500;
t23 = t560 + t640;
t21 = -t495 + t497 + t538;
t17 = t503 + t506 + t560;
t15 = -t731 / 0.2e1 - t732 / 0.2e1 + t195 * t565 + t552 / 0.2e1 - pkin(5) * t601 / 0.2e1 + t553 / 0.2e1 + t728 / 0.2e1 - t727 / 0.2e1 + t639;
t12 = t491 + t820;
t10 = t492 - t818;
t6 = t496 + (t519 + t834) * mrSges(7,3) + (t302 * t808 + t510) * t484 + pkin(9) * t541 + t515 * t810 + (t305 * t808 + ((m(7) * t752 + t776) * pkin(5) + t574) * t536 + t583) * t487 + t823 - t916;
t4 = (t519 + t833) * mrSges(7,3) + (t302 * t764 + t510) * t484 + t496 + (t305 * t764 + ((m(7) * t766 + t776) * pkin(5) + t574) * t536 + t583) * t487 + t431 * t541 + t516 * t810 + t822 - t914;
t1 = (t676 / 0.2e1 - t677 / 0.2e1) * mrSges(6,3) + t948 + t668 / 0.2e1 - t498 - t705 / 0.2e1 + t490 + (-t682 / 0.2e1 - t681 / 0.2e1) * mrSges(7,3) + Ifges(5,6) * t865 + (-t396 * t754 + t706 / 0.2e1 + (t126 * t433 + t149 * t431) * t812 + t526) * t487 + t542 * t812 + t642 * t778 + t737 * t747 + t736 * t751 + 0.2e1 * t883 * Ifges(5,5) + t868 * t746 + t869 * t749 + (t757 - mrSges(5,1) / 0.2e1) * t932 + ((-t125 * t433 - t148 * t431) * t812 - t707 / 0.2e1 + t463 * t781 + t525) * t484 + t60 * t606 + t61 * t607 + t647 * t781 + t830 * t536;
t26 = [qJD(2) * t3 + qJD(3) * t18 + qJD(4) * t8 + qJD(5) * t9 + qJD(6) * t14, t21 * qJD(3) + t1 * qJD(4) + t4 * qJD(5) + t10 * qJD(6) + t699 + (t948 + t668 + (-t434 * mrSges(5,3) + t561) * t536 - t705 + (-t129 * mrSges(6,3) - t431 * t889 + t901) * t484 - t322 * t888 + (t130 * mrSges(6,3) + t902 + t920) * t487 + 0.2e1 * (-t433 * t932 + t945) * t814 + 0.2e1 * (-t322 * t61 + t551 * t60 + t943) * t810 + t430 * t878 - t867 * mrSges(4,2) + 0.2e1 * (t411 * t483 + t482 * t867) * t634 - t60 * t685 - t61 * t686 + (-mrSges(3,1) * t488 + mrSges(3,2) * t486) * pkin(7) + (-t448 * t741 + t449 * t742) * mrSges(4,3) + (-t433 * mrSges(5,3) - t530) * t396 + 0.2e1 * (t431 * t544 + t937) * t812 + Ifges(3,5) * t488 - Ifges(3,6) * t486 + Ifges(4,5) * t448 + Ifges(4,6) * t449 + t411 * mrSges(4,1) + t425 * t926 + t551 * t910) * qJD(2), t673 + t21 * qJD(2) + (-t528 * t879 + t832 * t875) * t721 + t25 * qJD(4) + t17 * qJD(5) + t23 * qJD(6), t680 + t1 * qJD(2) + t25 * qJD(3) + t6 * qJD(5) + t12 * qJD(6) + ((pkin(9) * t543 - t938) * t812 + (-t418 * t74 + t550 * t73 + t942) * t810) * t815 + (t934 - t418 * t888 + t550 * t910 - t900 + (t706 + t902 + t737) * t487 + (-t707 + t901 - t736) * t484 + t561 * t536 - t530 * t396 + t848 * t932 + (-t682 - t681) * mrSges(7,3) + t948) * qJD(4), t679 + t4 * qJD(2) + t17 * qJD(3) + t6 * qJD(4) + (-Ifges(6,5) * t661 - Ifges(6,6) * t660 + t728 - t727 + (t69 * t744 + t70 * t743) * t806 + t552 + t553 - t125 * mrSges(6,2) - t126 * mrSges(6,1) + t639) * qJD(5) + t15 * qJD(6), t675 + t10 * qJD(2) + t23 * qJD(3) + t12 * qJD(4) + t15 * qJD(5) + (t639 - t731 - t732) * qJD(6); -qJD(3) * t20 + qJD(4) * t2 + qJD(5) * t5 + qJD(6) * t11 - t699, qJD(4) * t45 + qJD(5) * t44 + qJD(6) * t46, t631 - t678, t34 * qJD(5) + t42 * qJD(6) + ((-t344 * t550 + t345 * t418 + t434 * t474) * t810 + (-pkin(4) * t434 + pkin(9) * t568) * t812) * t815 - t521 + t827 * qJD(4), t34 * qJD(4) + ((t322 * t744 + t551 * t743) * t806 + t512 + t831 * t431 + t911) * qJD(5) + t895 - t523, t42 * qJD(4) + t75 * qJD(5) - t522 + t895; qJD(2) * t20 - qJD(4) * t24 - qJD(5) * t16 + qJD(6) * t22 - t673, t631 + t678, 0, t531, -t674 - t872 + (-t400 - t461 + (-t528 * t630 - t629 * t832) * m(7)) * qJD(5), -qJD(5) * t400 + t825 - t872; -qJD(2) * t2 + qJD(3) * t24 + qJD(5) * t7 + qJD(6) * t13 - t680, t33 * qJD(5) + t41 * qJD(6) + t521, -t531, qJD(5) * t47 + qJD(6) * t79 ((t418 * t744 + t550 * t743) * t806 + t512 + t831 * pkin(9) + t912) * qJD(5) + t893 - t514, t99 * qJD(5) - t513 + t893; -qJD(2) * t5 + qJD(3) * t16 - qJD(4) * t7 - qJD(6) * t19 - t679, -t33 * qJD(4) + t78 * qJD(6) + t523, t674, t109 * qJD(6) + t514, -t450, -t450 - t509; -qJD(2) * t11 - qJD(3) * t22 - qJD(4) * t13 + qJD(5) * t19 - t675, -t41 * qJD(4) - t78 * qJD(5) + t522, -t825, -t109 * qJD(5) + t513, t509, 0;];
Cq  = t26;
