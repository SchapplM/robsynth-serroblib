% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2019-03-09 15:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRRPPR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR5_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_coriolismatJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR5_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:38:35
% EndTime: 2019-03-09 15:39:32
% DurationCPUTime: 31.04s
% Computational Cost: add. (68478->1120), mult. (162343->1583), div. (0->0), fcn. (187873->12), ass. (0->540)
t660 = sin(pkin(12));
t662 = cos(pkin(12));
t663 = cos(pkin(6));
t664 = sin(qJ(3));
t661 = sin(pkin(6));
t665 = sin(qJ(2));
t798 = t661 * t665;
t888 = cos(qJ(3));
t611 = t663 * t888 - t664 * t798;
t612 = t663 * t664 + t798 * t888;
t825 = sin(pkin(11));
t826 = cos(pkin(11));
t693 = t611 * t825 + t612 * t826;
t666 = cos(qJ(2));
t797 = t661 * t666;
t440 = t660 * t797 - t662 * t693;
t703 = -t660 * t693 - t662 * t797;
t886 = sin(qJ(6));
t887 = cos(qJ(6));
t283 = t440 * t887 - t703 * t886;
t1024 = t283 / 0.2e1;
t1025 = mrSges(7,1) * t1024;
t742 = t826 * t664;
t626 = -t825 * t888 - t742;
t983 = t886 * t660 - t887 * t662;
t1015 = t983 * t626;
t1023 = t1015 / 0.2e1;
t1007 = t440 * t886 + t703 * t887;
t1016 = -t1007 / 0.2e1;
t1022 = mrSges(7,2) * t1016;
t696 = t660 * t887 + t662 * t886;
t895 = -t696 / 0.2e1;
t1021 = t283 * t895;
t741 = t825 * t664;
t624 = -t826 * t888 + t741;
t807 = t624 * t660;
t762 = t888 * qJ(4);
t778 = t888 * pkin(9);
t639 = t778 + t762;
t878 = -qJ(4) - pkin(9);
t988 = t826 * t639 + t878 * t741;
t474 = -pkin(5) * t807 + t988;
t913 = -t988 / 0.2e1;
t968 = -m(7) / 0.2e1;
t510 = t983 * t624;
t843 = t510 * mrSges(7,2);
t507 = t696 * t624;
t845 = t507 * mrSges(7,1);
t990 = t845 / 0.2e1 - t843 / 0.2e1;
t1020 = m(6) * t913 + t474 * t968 + t990;
t781 = pkin(8) * t797;
t884 = pkin(1) * t665;
t590 = t781 + (pkin(9) + t884) * t663;
t591 = (-pkin(2) * t666 - pkin(9) * t665 - pkin(1)) * t661;
t466 = -t664 * t590 + t888 * t591;
t417 = -t612 * qJ(4) + t466;
t467 = t590 * t888 + t664 * t591;
t418 = t611 * qJ(4) + t467;
t743 = t826 * t418;
t247 = t417 * t825 + t743;
t513 = -t826 * t611 + t612 * t825;
t812 = t513 * t660;
t190 = -pkin(5) * t812 + t247;
t970 = -m(6) / 0.2e1;
t342 = t983 * t513;
t849 = t342 * mrSges(7,2);
t341 = t696 * t513;
t850 = t341 * mrSges(7,1);
t993 = t850 / 0.2e1 - t849 / 0.2e1;
t1019 = t190 * t968 + t247 * t970 + t993;
t1009 = Ifges(5,4) * t513;
t347 = Ifges(5,1) * t693 - Ifges(5,5) * t797 - t1009;
t1018 = -Ifges(5,2) * t693 - t1009 + t347;
t1017 = t1015 * t895;
t965 = -mrSges(7,1) / 0.2e1;
t1000 = Ifges(5,4) * t693;
t995 = -Ifges(6,5) * t440 - Ifges(7,5) * t283 + Ifges(6,6) * t703 + Ifges(7,6) * t1007 + (Ifges(6,3) + Ifges(7,3)) * t513;
t1014 = -Ifges(5,1) * t513 - t1000 + t995;
t739 = -t626 * mrSges(5,1) - t624 * mrSges(5,2);
t1013 = -t440 / 0.2e1;
t1012 = -t513 / 0.2e1;
t1010 = t513 / 0.2e1;
t928 = t513 / 0.4e1;
t757 = t825 * pkin(3);
t648 = t757 + qJ(5);
t1008 = m(6) * t648 + mrSges(6,3);
t870 = Ifges(7,4) * t283;
t119 = Ifges(7,2) * t1007 + t513 * Ifges(7,6) - t870;
t1006 = -t119 / 0.2e1;
t209 = -mrSges(7,2) * t513 + mrSges(7,3) * t1007;
t958 = t209 / 0.2e1;
t691 = t696 * t626;
t869 = Ifges(7,4) * t1015;
t322 = Ifges(7,2) * t691 + t624 * Ifges(7,6) + t869;
t1005 = -t322 / 0.4e1;
t934 = -t691 / 0.2e1;
t935 = t691 / 0.2e1;
t1004 = t693 / 0.2e1;
t1001 = Ifges(5,3) + Ifges(4,3);
t999 = Ifges(3,6) * t798;
t646 = pkin(8) * t798;
t883 = pkin(1) * t666;
t615 = t663 * t883 - t646;
t998 = t615 * mrSges(3,2);
t841 = t693 * mrSges(5,3);
t785 = t660 ^ 2 + t662 ^ 2;
t976 = t785 * mrSges(6,3);
t736 = t785 * t648;
t583 = t624 * t797;
t527 = t583 * t660 + t662 * t798;
t529 = -t583 * t662 + t660 * t798;
t582 = t626 * t797;
t855 = Ifges(7,3) * t582;
t368 = t527 * t887 - t529 * t886;
t858 = Ifges(7,6) * t368;
t369 = t527 * t886 + t529 * t887;
t864 = Ifges(7,5) * t369;
t994 = Ifges(6,5) * t529 + Ifges(6,6) * t527 - Ifges(6,3) * t582 - t855 + t858 + t864;
t320 = Ifges(7,5) * t1015 + Ifges(7,6) * t691 + t624 * Ifges(7,3);
t860 = Ifges(6,6) * t660;
t866 = Ifges(6,5) * t662;
t722 = -t860 + t866;
t450 = t624 * Ifges(6,3) - t626 * t722;
t991 = t450 + t320;
t987 = Ifges(6,5) * t660 + Ifges(7,5) * t696 + Ifges(6,6) * t662 - Ifges(7,6) * t983;
t657 = Ifges(4,4) * t888;
t737 = -Ifges(4,2) * t664 + t657;
t856 = Ifges(7,3) * t693;
t859 = Ifges(7,6) * t341;
t865 = Ifges(7,5) * t342;
t986 = Ifges(6,3) * t693 - t513 * t722 + t856 + t859 + t865;
t641 = Ifges(4,1) * t664 + t657;
t985 = t737 + t641;
t608 = Ifges(4,4) * t611;
t479 = Ifges(4,1) * t612 - Ifges(4,5) * t797 + t608;
t984 = -Ifges(4,2) * t612 + t479 + t608;
t829 = t662 * mrSges(6,2);
t830 = t660 * mrSges(6,1);
t729 = t829 + t830;
t535 = t729 * t626;
t833 = t626 * mrSges(5,3);
t981 = -t833 / 0.2e1 - t535 / 0.2e1;
t980 = -t888 * mrSges(4,1) + t664 * mrSges(4,2);
t979 = Ifges(4,5) * t611 - Ifges(5,5) * t513 - Ifges(4,6) * t612 - Ifges(5,6) * t693;
t966 = m(5) * pkin(3);
t975 = t825 * t966 - mrSges(5,2);
t634 = -mrSges(6,1) * t662 + mrSges(6,2) * t660;
t759 = t826 * pkin(3);
t655 = -t759 - pkin(4);
t974 = m(6) * t655 - t826 * t966 - mrSges(5,1) + t634;
t972 = 0.2e1 * t626;
t971 = m(5) / 0.2e1;
t969 = m(6) / 0.2e1;
t967 = m(7) / 0.2e1;
t964 = mrSges(7,2) / 0.2e1;
t963 = mrSges(6,3) / 0.2e1;
t278 = Ifges(7,4) * t1007;
t120 = -Ifges(7,1) * t283 + t513 * Ifges(7,5) + t278;
t962 = t120 / 0.2e1;
t961 = t120 / 0.4e1;
t139 = -mrSges(7,1) * t1007 - mrSges(7,2) * t283;
t960 = t139 / 0.2e1;
t402 = -pkin(3) * t797 + t417;
t410 = t825 * t418;
t238 = t402 * t826 - t410;
t230 = pkin(4) * t797 - t238;
t159 = -pkin(5) * t703 + t230;
t959 = t159 / 0.2e1;
t210 = mrSges(7,1) * t513 + mrSges(7,3) * t283;
t957 = -t210 / 0.2e1;
t954 = t1007 / 0.2e1;
t951 = -t283 / 0.2e1;
t296 = -mrSges(6,1) * t703 - mrSges(6,2) * t440;
t950 = t296 / 0.2e1;
t496 = Ifges(7,4) * t691;
t324 = Ifges(7,1) * t1015 + t624 * Ifges(7,5) + t496;
t949 = t324 / 0.2e1;
t948 = t324 / 0.4e1;
t325 = -mrSges(6,2) * t513 + mrSges(6,3) * t703;
t947 = t325 / 0.2e1;
t946 = t368 / 0.2e1;
t945 = t369 / 0.2e1;
t426 = -mrSges(7,2) * t624 + mrSges(7,3) * t691;
t944 = -t426 / 0.2e1;
t943 = t426 / 0.2e1;
t428 = mrSges(7,1) * t624 - mrSges(7,3) * t1015;
t942 = -t428 / 0.2e1;
t940 = t703 / 0.2e1;
t459 = mrSges(5,2) * t797 - t513 * mrSges(5,3);
t938 = t459 / 0.2e1;
t575 = t639 * t825 - t878 * t742;
t805 = t626 * t660;
t475 = -pkin(5) * t805 + t575;
t937 = t475 / 0.2e1;
t923 = t527 / 0.2e1;
t879 = pkin(10) + t648;
t613 = t879 * t660;
t614 = t879 * t662;
t528 = -t613 * t887 - t614 * t886;
t922 = t528 / 0.2e1;
t921 = t529 / 0.2e1;
t530 = -t613 * t886 + t614 * t887;
t920 = -t530 / 0.2e1;
t549 = mrSges(7,1) * t696 - mrSges(7,2) * t983;
t918 = t549 / 0.2e1;
t551 = mrSges(7,1) * t983 + mrSges(7,2) * t696;
t917 = t551 / 0.2e1;
t868 = Ifges(7,4) * t696;
t555 = -Ifges(7,2) * t983 + t868;
t916 = t555 / 0.2e1;
t622 = Ifges(7,4) * t983;
t558 = Ifges(7,1) * t696 - t622;
t915 = t558 / 0.2e1;
t914 = t558 / 0.4e1;
t912 = t582 / 0.2e1;
t911 = -t582 / 0.2e1;
t909 = -t583 / 0.2e1;
t908 = t611 / 0.2e1;
t906 = t612 / 0.2e1;
t905 = -t624 / 0.2e1;
t903 = t624 / 0.2e1;
t902 = t624 / 0.4e1;
t901 = -t626 / 0.2e1;
t899 = t626 / 0.2e1;
t898 = t983 / 0.2e1;
t897 = -t983 / 0.2e1;
t896 = t696 / 0.2e1;
t633 = -t662 * pkin(5) + t655;
t894 = t633 / 0.2e1;
t893 = t634 / 0.2e1;
t892 = t660 / 0.2e1;
t891 = t662 / 0.2e1;
t890 = t664 / 0.2e1;
t889 = -t666 / 0.2e1;
t882 = pkin(3) * t612;
t881 = pkin(3) * t664;
t880 = pkin(9) * t664;
t877 = m(7) * qJD(4);
t876 = Ifges(3,4) * t665;
t875 = Ifges(4,4) * t612;
t874 = Ifges(4,4) * t664;
t873 = Ifges(5,4) * t626;
t872 = Ifges(6,4) * t660;
t871 = Ifges(6,4) * t662;
t867 = Ifges(5,5) * t583;
t863 = Ifges(7,5) * t510;
t861 = Ifges(5,6) * t582;
t857 = Ifges(7,6) * t507;
t854 = Ifges(7,3) * t626;
t853 = t1007 * mrSges(7,2);
t852 = t283 * mrSges(7,1);
t239 = t825 * t402 + t743;
t229 = -qJ(5) * t797 + t239;
t589 = t646 + (-pkin(2) - t883) * t663;
t533 = -t611 * pkin(3) + t589;
t265 = t513 * pkin(4) - qJ(5) * t693 + t533;
t112 = -t229 * t660 + t662 * t265;
t113 = t662 * t229 + t660 * t265;
t248 = t417 * t826 - t410;
t315 = pkin(4) * t693 + qJ(5) * t513 + t882;
t131 = -t248 * t660 + t662 * t315;
t132 = t662 * t248 + t660 * t315;
t143 = Ifges(7,4) * t342 + Ifges(7,2) * t341 + Ifges(7,6) * t693;
t144 = Ifges(7,1) * t342 + Ifges(7,4) * t341 + Ifges(7,5) * t693;
t181 = t849 - t850;
t219 = -Ifges(6,4) * t440 + Ifges(6,2) * t703 + t513 * Ifges(6,6);
t220 = -Ifges(6,1) * t440 + Ifges(6,4) * t703 + t513 * Ifges(6,5);
t245 = -mrSges(7,2) * t693 + mrSges(7,3) * t341;
t246 = mrSges(7,1) * t693 - mrSges(7,3) * t342;
t724 = -Ifges(6,2) * t660 + t871;
t249 = Ifges(6,6) * t693 - t513 * t724;
t726 = Ifges(6,1) * t662 - t872;
t250 = Ifges(6,5) * t693 - t513 * t726;
t326 = mrSges(6,1) * t513 + mrSges(6,3) * t440;
t346 = -Ifges(5,2) * t513 - Ifges(5,6) * t797 + t1000;
t348 = t729 * t513;
t355 = -mrSges(6,2) * t693 + mrSges(6,3) * t812;
t811 = t513 * t662;
t356 = mrSges(6,1) * t693 + mrSges(6,3) * t811;
t364 = mrSges(5,1) * t513 + mrSges(5,2) * t693;
t478 = Ifges(4,2) * t611 - Ifges(4,6) * t797 + t875;
t836 = t611 * mrSges(4,3);
t578 = mrSges(4,2) * t797 + t836;
t835 = t612 * mrSges(4,3);
t579 = -mrSges(4,1) * t797 - t835;
t86 = pkin(5) * t513 + pkin(10) * t440 + t112;
t93 = pkin(10) * t703 + t113;
t60 = t86 * t887 - t886 * t93;
t61 = t86 * t886 + t887 * t93;
t106 = pkin(10) * t812 + t132;
t94 = pkin(5) * t693 + pkin(10) * t811 + t131;
t67 = -t106 * t886 + t887 * t94;
t68 = t106 * t887 + t886 * t94;
t728 = Ifges(4,1) * t611 - t875;
t730 = mrSges(4,1) * t612 + mrSges(4,2) * t611;
t740 = mrSges(5,1) * t693 - mrSges(5,2) * t513;
t755 = -t811 / 0.2e1;
t756 = t812 / 0.2e1;
t460 = -mrSges(5,1) * t797 - t841;
t788 = -t460 + t296;
t3 = m(6) * (t112 * t131 + t113 * t132 + t230 * t247) + m(7) * (t159 * t190 + t60 * t67 + t61 * t68) - t693 * t346 / 0.2e1 + (-t466 * t611 - t467 * t612) * mrSges(4,3) - t341 * t1006 + (t238 * t513 - t239 * t693) * mrSges(5,3) + t979 * t661 * t889 + t220 * t755 + t219 * t756 + m(5) * (-t238 * t247 + t239 * t248 + t533 * t882) + t250 * t1013 + t986 * t1010 + t788 * t247 - t612 * t478 / 0.2e1 + t466 * t578 - t467 * t579 + t248 * t459 + t113 * t355 + t112 * t356 - t230 * t348 + t132 * t325 + t131 * t326 + t1018 * t1012 + t61 * t245 + t60 * t246 + t67 * t210 + t68 * t209 + t190 * t139 + t159 * t181 + t1014 * t1004 + t533 * t740 + t143 * t954 + t342 * t962 + t589 * t730 + t364 * t882 + t728 * t906 + t249 * t940 + t144 * t951 + t984 * t908;
t851 = t3 * qJD(1);
t848 = t368 * mrSges(7,1);
t847 = t369 * mrSges(7,2);
t616 = (pkin(2) * t665 - pkin(9) * t666) * t661;
t531 = -t664 * t615 + t888 * t616;
t439 = (pkin(3) * t665 - t666 * t762) * t661 + t531;
t532 = t888 * t615 + t664 * t616;
t792 = t664 * t666;
t771 = t661 * t792;
t465 = -qJ(4) * t771 + t532;
t298 = t825 * t439 + t826 * t465;
t268 = qJ(5) * t798 + t298;
t617 = t663 * t884 + t781;
t581 = pkin(3) * t771 + t617;
t383 = -t582 * pkin(4) + t583 * qJ(5) + t581;
t160 = -t268 * t660 + t662 * t383;
t161 = t662 * t268 + t660 * t383;
t177 = Ifges(7,4) * t369 + Ifges(7,2) * t368 - Ifges(7,6) * t582;
t178 = Ifges(7,1) * t369 + Ifges(7,4) * t368 - Ifges(7,5) * t582;
t196 = t847 - t848;
t297 = t439 * t826 - t825 * t465;
t269 = -pkin(4) * t798 - t297;
t211 = -t527 * pkin(5) + t269;
t292 = mrSges(7,2) * t582 + mrSges(7,3) * t368;
t293 = -mrSges(7,1) * t582 - mrSges(7,3) * t369;
t305 = Ifges(6,4) * t529 + Ifges(6,2) * t527 - Ifges(6,6) * t582;
t306 = Ifges(6,1) * t529 + Ifges(6,4) * t527 - Ifges(6,5) * t582;
t839 = t529 * mrSges(6,2);
t840 = t527 * mrSges(6,1);
t374 = t839 - t840;
t413 = mrSges(6,2) * t582 + mrSges(6,3) * t527;
t414 = -mrSges(6,1) * t582 - mrSges(6,3) * t529;
t420 = -Ifges(5,4) * t583 + Ifges(5,2) * t582 + Ifges(5,6) * t798;
t421 = -Ifges(5,1) * t583 + Ifges(5,4) * t582 + Ifges(5,5) * t798;
t837 = t583 * mrSges(5,2);
t838 = t582 * mrSges(5,1);
t444 = -t837 - t838;
t539 = -mrSges(5,2) * t798 + mrSges(5,3) * t582;
t540 = mrSges(5,1) * t798 + mrSges(5,3) * t583;
t563 = (Ifges(4,6) * t665 + t666 * t737) * t661;
t709 = Ifges(4,1) * t888 - t874;
t564 = (Ifges(4,5) * t665 + t666 * t709) * t661;
t706 = t664 * mrSges(4,1) + mrSges(4,2) * t888;
t587 = t706 * t797;
t604 = (-mrSges(4,2) * t665 - mrSges(4,3) * t792) * t661;
t605 = (-mrSges(4,3) * t666 * t888 + mrSges(4,1) * t665) * t661;
t708 = Ifges(4,5) * t888 - Ifges(4,6) * t664;
t695 = t708 * t666;
t749 = t797 / 0.2e1;
t714 = t888 * t749;
t731 = Ifges(3,5) * t797 - t999;
t750 = -t797 / 0.2e1;
t732 = t664 * t750;
t123 = -pkin(5) * t582 - pkin(10) * t529 + t160;
t134 = pkin(10) * t527 + t161;
t75 = t123 * t887 - t134 * t886;
t751 = t798 / 0.2e1;
t76 = t123 * t886 + t134 * t887;
t4 = t421 * t1004 + (0.2e1 * Ifges(3,4) * t797 + (Ifges(3,1) - Ifges(3,2)) * t798) * t749 + t479 * t714 + t478 * t732 + (-(Ifges(3,2) * t666 + t876) * t798 / 0.2e1 + (-pkin(1) * (mrSges(3,1) * t665 + mrSges(3,2) * t666) + t665 * (Ifges(3,1) * t666 - t876) / 0.2e1 + (Ifges(4,3) * t665 + t695) * t889) * t661) * t661 + t306 * t1013 + t994 * t1010 + t420 * t1012 + t617 * (-mrSges(4,1) * t611 + mrSges(4,2) * t612) + t467 * t604 + t466 * t605 + t589 * t587 + t581 * t364 + t532 * t578 + t531 * t579 + t239 * t539 + t238 * t540 + t533 * t444 + t298 * t459 + t297 * t460 + t113 * t413 + t112 * t414 + t230 * t374 + t161 * t325 + t160 * t326 + t269 * t296 + t61 * t292 + t60 * t293 + t75 * t210 + t211 * t139 + t76 * t209 + t159 * t196 + (Ifges(5,3) * t798 + t861 - t867) * t750 + t177 * t954 + m(6) * (t112 * t160 + t113 * t161 + t230 * t269) + m(7) * (t159 * t211 + t60 * t75 + t61 * t76) + m(5) * (t238 * t297 + t239 * t298 + t533 * t581) + m(4) * (t466 * t531 + t467 * t532 + t589 * t617) + t564 * t906 + t563 * t908 + t347 * t909 + t346 * t912 + t220 * t921 + t219 * t923 + t305 * t940 + t120 * t945 + t119 * t946 + t178 * t951 + t995 * t911 + (-t999 / 0.2e1 + Ifges(3,5) * t749 - t617 * mrSges(3,1) - t998 + t731 / 0.2e1) * t663 + (Ifges(4,5) * t612 + Ifges(5,5) * t693 + Ifges(4,6) * t611 - Ifges(5,6) * t513 - t1001 * t797) * t751;
t846 = t4 * qJD(1);
t844 = t691 * mrSges(7,2);
t842 = t1015 * mrSges(7,1);
t834 = t624 * mrSges(5,3);
t832 = t983 * mrSges(7,3);
t831 = t696 * mrSges(7,3);
t138 = -t852 + t853;
t140 = Ifges(7,5) * t1007 + Ifges(7,6) * t283;
t141 = Ifges(7,2) * t283 + t278;
t142 = Ifges(7,1) * t1007 + t870;
t9 = t60 * t209 - t61 * t210 + t159 * t138 + t140 * t1010 - (-t61 * mrSges(7,3) + t142 / 0.2e1 + t1006) * t283 + (-t60 * mrSges(7,3) + t962 + t141 / 0.2e1) * t1007;
t827 = t9 * qJD(1);
t24 = t1007 * t209 + t283 * t210 + t703 * t325 + t440 * t326 + m(7) * (t1007 * t61 + t283 * t60) + m(6) * (t112 * t440 + t113 * t703);
t824 = qJD(1) * t24;
t716 = t112 * t660 - t113 * t662;
t796 = t662 * t325;
t801 = t660 * t326;
t15 = t342 * t209 + t341 * t210 - (t459 + t796 - t801) * t513 + (t139 + t788) * t693 + m(7) * (t159 * t693 + t341 * t60 + t342 * t61) + m(6) * (t230 * t693 + t513 * t716) + m(5) * (-t238 * t693 - t239 * t513);
t821 = t15 * qJD(1);
t820 = t160 * t660;
t819 = t161 * t662;
t818 = t247 * t575;
t815 = t440 * t660;
t814 = t703 * t662;
t813 = t693 * t575;
t808 = t575 * t626;
t806 = t624 * t662;
t804 = t626 * t662;
t803 = t648 * t660;
t802 = t648 * t662;
t452 = Ifges(6,6) * t624 - t626 * t724;
t800 = t660 * t452;
t544 = mrSges(6,1) * t624 + mrSges(6,3) * t804;
t799 = t660 * t544;
t795 = t662 * t326;
t454 = Ifges(6,5) * t624 - t626 * t726;
t794 = t662 * t454;
t542 = -mrSges(6,2) * t624 + mrSges(6,3) * t805;
t793 = t662 * t542;
t365 = Ifges(7,5) * t691 - Ifges(7,6) * t1015;
t538 = -pkin(4) * t626 + qJ(5) * t624 + t881;
t406 = t660 * t538 - t575 * t662;
t656 = -pkin(3) * t888 - pkin(2);
t537 = t624 * pkin(4) + t626 * qJ(5) + t656;
t404 = t660 * t537 + t662 * t988;
t552 = -Ifges(7,5) * t983 - Ifges(7,6) * t696;
t784 = t549 * qJD(6);
t783 = t966 / 0.2e1;
t782 = -m(7) / 0.4e1 - m(6) / 0.4e1;
t780 = t888 / 0.2e1;
t774 = t882 / 0.2e1;
t773 = t881 / 0.2e1;
t768 = mrSges(5,3) * t1012;
t767 = -t834 / 0.2e1;
t765 = -t832 / 0.2e1;
t764 = -t831 / 0.2e1;
t763 = t830 / 0.2e1;
t754 = t807 / 0.2e1;
t753 = -t806 / 0.2e1;
t752 = t805 / 0.2e1;
t748 = t961 + t141 / 0.4e1;
t747 = t142 / 0.4e1 - t119 / 0.4e1;
t366 = -Ifges(7,2) * t1015 + t496;
t746 = t948 + t366 / 0.4e1;
t367 = Ifges(7,1) * t691 - t869;
t745 = t367 / 0.4e1 + t1005;
t744 = m(6) * t785;
t403 = t662 * t537 - t660 * t988;
t405 = t662 * t538 + t575 * t660;
t735 = mrSges(5,3) * t759;
t734 = mrSges(5,3) * t757;
t723 = -Ifges(5,5) * t624 + Ifges(5,6) * t626;
t311 = pkin(5) * t624 + pkin(10) * t804 + t403;
t352 = pkin(10) * t805 + t404;
t170 = t311 * t887 - t352 * t886;
t171 = t311 * t886 + t352 * t887;
t312 = -pkin(5) * t626 + pkin(10) * t806 + t405;
t354 = pkin(10) * t807 + t406;
t174 = t312 * t887 - t354 * t886;
t175 = t312 * t886 + t354 * t887;
t321 = Ifges(7,4) * t510 + Ifges(7,2) * t507 - Ifges(7,6) * t626;
t323 = Ifges(7,1) * t510 + Ifges(7,4) * t507 - Ifges(7,5) * t626;
t362 = t843 - t845;
t363 = -mrSges(7,1) * t691 + mrSges(7,2) * t1015;
t425 = mrSges(7,2) * t626 + mrSges(7,3) * t507;
t427 = -mrSges(7,1) * t626 - mrSges(7,3) * t510;
t451 = -Ifges(6,6) * t626 - t624 * t724;
t453 = -Ifges(6,5) * t626 - t624 * t726;
t534 = t729 * t624;
t541 = mrSges(6,2) * t626 + mrSges(6,3) * t807;
t543 = -mrSges(6,1) * t626 + mrSges(6,3) * t806;
t550 = mrSges(5,1) * t624 - mrSges(5,2) * t626;
t556 = -t624 * Ifges(5,2) - t873;
t623 = Ifges(5,4) * t624;
t559 = -t626 * Ifges(5,1) - t623;
t640 = Ifges(4,2) * t888 + t874;
t707 = t863 / 0.2e1 + t857 / 0.2e1;
t12 = -t988 * t535 + m(7) * (t170 * t174 + t171 * t175 + t474 * t475) + t323 * t1023 + (-t450 / 0.2e1 - t320 / 0.2e1 + t556 / 0.2e1 + t451 * t892 - t662 * t453 / 0.2e1 - t873 / 0.2e1 + (-Ifges(6,3) / 0.2e1 + Ifges(5,1) / 0.2e1 - Ifges(7,3) / 0.2e1 - Ifges(5,2) / 0.2e1) * t624) * t626 + m(6) * (t403 * t405 + t404 * t406 + t988 * t575) - t664 * t640 / 0.2e1 - t575 * t534 + t404 * t541 + t406 * t542 + t403 * t543 + t405 * t544 + t507 * t322 / 0.2e1 + t474 * t363 + t475 * t362 + t171 * t425 + t175 * t426 + t170 * t427 + t174 * t428 + (m(5) * t881 + t739) * t656 - pkin(2) * t706 + (-t559 / 0.2e1 + t623 / 0.2e1 + t800 / 0.2e1 - t794 / 0.2e1 + (-t866 / 0.2e1 + t860 / 0.2e1) * t624 + t707) * t624 + t550 * t881 + t709 * t890 + t321 * t935 + t510 * t949 + t985 * t780;
t667 = (t841 + t460) * t913 + t664 * t728 / 0.4e1 - t693 * t556 / 0.4e1 - t661 * t695 / 0.4e1 - t341 * t1005 + (t112 * t405 + t113 * t406 + t131 * t403 + t132 * t404 + t230 * t988 + t818) * t969 + t988 * t950 + t364 * t773 + t550 * t774 + t248 * t767 - (t579 + t835) * t778 / 0.2e1 + t703 * t451 / 0.4e1 + t691 * t143 / 0.4e1 - t283 * t323 / 0.4e1 + t1015 * t144 / 0.4e1 - (Ifges(5,2) * t626 + t559 - t623 + t794) * t513 / 0.4e1 + (t836 / 0.2e1 - t578 / 0.2e1) * t880 + t239 * t833 / 0.2e1 + t238 * t834 / 0.2e1 + t1007 * t321 / 0.4e1 + t249 * t805 / 0.4e1 - t220 * t806 / 0.4e1 + t219 * t807 / 0.4e1 - t250 * t804 / 0.4e1 - t723 * t797 / 0.4e1 - t664 * t478 / 0.4e1 - t612 * t640 / 0.4e1 + t113 * t541 / 0.2e1 + t132 * t542 / 0.2e1 + t112 * t543 / 0.2e1 + t131 * t544 / 0.2e1 - t230 * t534 / 0.2e1 + t507 * t119 / 0.4e1 + t61 * t425 / 0.2e1 + t60 * t427 / 0.2e1 + t67 * t428 / 0.2e1 + t403 * t356 / 0.2e1 + t404 * t355 / 0.2e1 + t405 * t326 / 0.2e1 + t190 * t363 / 0.2e1 + (-t624 * t722 + t800 - t854 + t857 + t863) * t928 - t1018 * t624 / 0.4e1 + t171 * t245 / 0.2e1 + t170 * t246 / 0.2e1 + t174 * t210 / 0.2e1 + t589 * t706 / 0.2e1 + t612 * t709 / 0.4e1 + t533 * t739 / 0.2e1 + t656 * t740 / 0.2e1 + (t768 - t348 / 0.2e1 - t938 - t239 * t971) * t575 + t175 * t958 + t362 * t959 + t474 * t960 + t510 * t961 + (t159 * t474 + t170 * t67 + t171 * t68 + t174 * t60 + t175 * t61 + t190 * t475) * t967 - pkin(2) * t730 / 0.2e1 + (t346 / 0.4e1 - t1014 / 0.4e1 - Ifges(6,3) * t928) * t626 + (-Ifges(5,1) * t624 + t873 + t991) * t693 / 0.4e1 - t440 * t453 / 0.4e1 + t181 * t937 + t68 * t943 + t406 * t947 + t342 * t948 + t981 * t247 + t984 * t888 / 0.4e1 + t985 * t611 / 0.4e1 + t986 * t902 + (t818 + (t533 * t664 + t612 * t656) * pkin(3) + (-t238 + t248) * t988) * t971;
t636 = Ifges(6,2) * t662 + t872;
t637 = Ifges(6,1) * t660 + t871;
t668 = t696 * t178 / 0.4e1 + t75 * t764 + t76 * t765 + Ifges(4,5) * t714 + Ifges(4,6) * t732 - t867 / 0.2e1 + t861 / 0.2e1 + (t297 * t826 + t298 * t825) * t783 - mrSges(6,3) * t820 / 0.2e1 - t983 * t177 / 0.4e1 - t414 * t803 / 0.2e1 + t413 * t802 / 0.2e1 + t540 * t759 / 0.2e1 + t539 * t757 / 0.2e1 + t662 * t305 / 0.4e1 + t660 * t306 / 0.4e1 + t655 * t374 / 0.2e1 + t527 * t636 / 0.4e1 + t529 * t637 / 0.4e1 + t368 * t555 / 0.4e1 + t530 * t292 / 0.2e1 + t531 * mrSges(4,1) / 0.2e1 - t532 * mrSges(4,2) / 0.2e1 + t297 * mrSges(5,1) / 0.2e1 - t298 * mrSges(5,2) / 0.2e1 + (t269 * t655 + (t819 - t820) * t648) * t969 + t819 * t963 + (t211 * t633 + t528 * t75 + t530 * t76) * t967 + t269 * t893 + t196 * t894 + t369 * t914 + t211 * t917 + t293 * t922 - t987 * t582 / 0.4e1 + t1001 * t751;
t2 = t668 - t667;
t721 = -t2 * qJD(1) + t12 * qJD(2);
t361 = t842 + t844;
t29 = t365 * t903 + t475 * t361 - t171 * t428 + t170 * t426 + (t367 / 0.2e1 - t322 / 0.2e1 - t171 * mrSges(7,3)) * t1015 + (t949 + t366 / 0.2e1 - t170 * mrSges(7,3)) * t691;
t670 = t748 * t691 + t747 * t1015 + t746 * t1007 - t745 * t283 + (t170 * t1016 + t171 * t1024 + t60 * t934 - t61 * t1015 / 0.2e1) * mrSges(7,3) + t361 * t959 + t170 * t958 + t171 * t957 + t138 * t937 + t365 * t928 + t60 * t943 + t61 * t942 + t140 * t902;
t686 = -t864 / 0.2e1 - t858 / 0.2e1 + t855 / 0.2e1 + t75 * t965 + t76 * t964;
t6 = t670 + t686;
t720 = t6 * qJD(1) + t29 * qJD(2);
t715 = t403 * t660 - t404 * t662;
t669 = -m(5) * (t238 * t626 - t239 * t624 - t513 * t988 + t813) / 0.2e1 + (-t230 * t626 + t513 * t715 + t624 * t716 + t813) * t970 + (-t159 * t626 + t170 * t341 + t171 * t342 + t475 * t693 + t507 * t60 + t510 * t61) * t968 + t341 * t942 + t342 * t944 + t507 * t957 - t510 * t958;
t678 = t581 * t971 + (t160 * t662 + t161 * t660) * t969 + (t696 * t76 - t75 * t983) * t967 - t838 / 0.2e1 - t837 / 0.2e1 + t293 * t897 + t292 * t896 + t413 * t892 + t414 * t891;
t699 = t799 / 0.2e1 - t793 / 0.2e1;
t700 = -t801 / 0.2e1 + t796 / 0.2e1;
t10 = t669 + (t535 / 0.2e1 - t363 / 0.2e1) * t693 - t699 * t513 + (t950 + t960 - t460 / 0.2e1 + t841 / 0.2e1) * t626 + (t938 + t768 + t700) * t624 + t678;
t35 = t510 * t426 + t507 * t428 + (-t363 + t535 + t833) * t626 + (-t793 + t799 + t834) * t624 + m(7) * (t170 * t507 + t171 * t510 - t475 * t626) + m(6) * (t624 * t715 - t808) + m(5) * (-t624 * t988 - t808);
t719 = -t10 * qJD(1) + t35 * qJD(2);
t673 = (t403 * t440 + t404 * t703 + (t112 * t662 + t113 * t660) * t626) * t970 + (t1007 * t171 - t1015 * t60 + t170 * t283 + t61 * t691) * t968 + t283 * t942 + t1007 * t944 + t544 * t1013 - t703 * t542 / 0.2e1 + t210 * t1023 + t209 * t934;
t683 = t269 * t969 + t211 * t967 - t848 / 0.2e1 + t847 / 0.2e1 - t840 / 0.2e1 + t839 / 0.2e1;
t13 = (-t795 / 0.2e1 - t660 * t325 / 0.2e1) * t626 + t673 + t683;
t45 = m(7) * (-t1015 * t170 + t171 * t691) - t1015 * t428 + t691 * t426 + (m(6) * (t403 * t662 + t404 * t660) + t660 * t542 + t662 * t544) * t626;
t718 = qJD(1) * t13 - qJD(2) * t45;
t672 = (t341 * t528 + t342 * t530) * t967 + t341 * t764 + t342 * t765 + t1012 * t976 - (t736 * t969 + t783 * t825) * t513 + (t633 * t967 + t655 * t969 - t783 * t826 + t893 + t917) * t693;
t680 = (t131 * t662 + t132 * t660) * t969 + (-t67 * t983 + t68 * t696) * t967 + t246 * t897 + t245 * t896 + t355 * t892 + t356 * t891 + m(5) * t774;
t22 = t672 - t680 - t740;
t671 = (-t624 * t736 - t626 * t655) * t969 + (t507 * t528 + t510 * t530 - t626 * t633) * t967 + (-t624 * t825 + t626 * t826) * t783 + t507 * t764 + t510 * t765 + (t551 + t634) * t901 + t905 * t976;
t679 = (t405 * t662 + t406 * t660) * t969 + (-t174 * t983 + t175 * t696) * t967 + t427 * t897 + t425 * t896 + t541 * t892 + t543 * t891 + m(5) * t773;
t34 = -t671 + t679 + t739;
t717 = -qJD(1) * t22 + qJD(2) * t34;
t698 = m(7) * (t1015 * t983 + t691 * t696);
t186 = -t698 / 0.2e1 + (-t744 / 0.4e1 + t782) * t972;
t684 = (t1007 * t696 - t283 * t983) * t967 + (t440 * t662 + t660 * t703) * t969;
t80 = 0.2e1 * t693 * t782 + t684;
t713 = t80 * qJD(1) - t186 * qJD(2);
t557 = -Ifges(7,1) * t983 - t868;
t705 = t557 / 0.4e1 - t555 / 0.4e1 + mrSges(7,3) * t920;
t554 = -Ifges(7,2) * t696 - t622;
t704 = t914 + t554 / 0.4e1 - t528 * mrSges(7,3) / 0.2e1;
t702 = t209 * t897 + t210 * t895;
t701 = t426 * t897 + t428 * t895;
t676 = (t1007 * t897 + t1021) * mrSges(7,3) + ((t814 - t815) * t648 - t716) * t969 + (t1007 * t530 + t283 * t528 - t60 * t696 - t61 * t983) * t967 + t702;
t17 = (mrSges(6,2) * t1010 + t703 * t963 + t947) * t662 + (-t326 / 0.2e1 + mrSges(6,3) * t1013 + mrSges(6,1) * t1010) * t660 + t676 + t1019;
t184 = (t696 ^ 2 + t983 ^ 2) * mrSges(7,3) + t976 + m(7) * (-t528 * t696 - t530 * t983) + m(6) * t736;
t677 = (t691 * t897 - t1017) * mrSges(7,3) + t715 * t970 + (-t1015 * t528 - t170 * t696 - t171 * t983 + t530 * t691) * t967 - t699 + t701;
t40 = (t829 / 0.2e1 + t763) * t624 + t677 + t1020;
t694 = t17 * qJD(1) + t40 * qJD(2) + t184 * qJD(3);
t199 = 0.2e1 * mrSges(7,1) * t1023 + 0.2e1 * t935 * mrSges(7,2);
t81 = 0.2e1 * t1022 + 0.2e1 * t1025;
t692 = qJD(1) * t81 - qJD(2) * t199 - qJD(3) * t549;
t682 = (t1007 * t898 - t1021) * mrSges(7,3) + t702;
t20 = t682 - t993;
t681 = (t691 * t898 + t1017) * mrSges(7,3) + t701;
t41 = t681 - t990;
t690 = -t20 * qJD(1) - t41 * qJD(2);
t674 = t1015 * t705 + t361 * t894 + t426 * t922 + t428 * t920 + t475 * t918 + t552 * t902 + t691 * t704 + t696 * t745 - t746 * t983;
t685 = t854 / 0.2e1 + t174 * t965 + t175 * t964 - t707;
t26 = t674 + t685;
t78 = t633 * t549 - (-t557 / 0.2e1 + t916) * t696 - (t915 + t554 / 0.2e1) * t983;
t675 = t1007 * t704 + t138 * t894 + t159 * t918 + t209 * t922 + t210 * t920 - t283 * t705 + t552 * t928 + t696 * t747 - t748 * t983;
t687 = -t865 / 0.2e1 - t859 / 0.2e1 - t856 / 0.2e1 + t67 * t965 + t68 * t964;
t8 = t675 + t687;
t688 = -t8 * qJD(1) - t26 * qJD(2) - t78 * qJD(3);
t200 = t844 / 0.2e1 + t842 / 0.2e1 + mrSges(7,2) * t934 + t1015 * t965;
t187 = t698 / 0.2e1 + t744 * t899 + t782 * t972;
t82 = t853 / 0.2e1 - t852 / 0.2e1 + t1025 + t1022;
t79 = t684 + (m(6) + m(7)) * t1004;
t42 = t681 + t990;
t39 = mrSges(6,2) * t753 - mrSges(6,1) * t807 / 0.2e1 + t677 - t1020;
t36 = t671 + t679;
t25 = t674 - t685;
t23 = t672 + t680;
t19 = t682 + t993;
t16 = mrSges(6,2) * t755 - t513 * t763 + (t814 / 0.2e1 - t815 / 0.2e1) * mrSges(6,3) + t676 + t700 - t1019;
t14 = t325 * t752 + t795 * t899 - t673 + t683;
t11 = t325 * t753 + t326 * t754 + t363 * t1004 + t459 * t905 + t460 * t899 - t513 * t767 + t542 * t755 + t544 * t756 - t669 + t678 + (t139 + t296) * t901 + t981 * t693;
t7 = t675 - t687;
t5 = t670 - t686;
t1 = t668 + t667;
t18 = [qJD(2) * t4 + qJD(3) * t3 + qJD(4) * t15 + qJD(5) * t24 + qJD(6) * t9, t846 + (t731 + m(5) * (t298 * t988 + t581 * t656) + t988 * t539 + m(6) * (t160 * t403 + t161 * t404) + t305 * t752 + (Ifges(4,5) * t664 - Ifges(5,5) * t626 + Ifges(4,6) * t888 - Ifges(5,6) * t624) * t751 + t297 * t833 + t604 * t778 + t563 * t780 + (-m(5) * t297 + m(6) * t269 + t374 - t540) * t575 - t998 + t178 * t1023 + t641 * t714 + t640 * t732 - t605 * t880 - t298 * t834 - t306 * t804 / 0.2e1 + t656 * t444 - pkin(2) * t587 + t581 * t550 - t269 * t535 + t161 * t542 + t160 * t544 + t475 * t196 + t76 * t426 + t75 * t428 + t404 * t413 + t403 * t414 + t211 * t363 + t171 * t292 + t170 * t293 + (m(4) * pkin(9) + mrSges(4,3)) * (-t531 * t664 + t888 * t532) + t564 * t890 + t421 * t901 + t420 * t905 + t559 * t909 + t556 * t912 + t454 * t921 + t452 * t923 + t177 * t935 + t324 * t945 + t322 * t946 + m(7) * (t170 * t75 + t171 * t76 + t211 * t475) + (-m(4) * pkin(2) - mrSges(3,1) + t980) * t617 + t991 * t911 + t994 * t903) * qJD(2) + t1 * qJD(3) + t11 * qJD(4) + t14 * qJD(5) + t5 * qJD(6), t851 + t1 * qJD(2) + (-t693 * t734 + t513 * t735 + t987 * t1004 + t637 * t755 + t636 * t756 - t67 * t831 - t68 * t832 + t355 * t802 + t979 - t356 * t803 - t655 * t348 + t633 * t181 + t190 * t551 + t530 * t245 + t528 * t246 - t466 * mrSges(4,2) - t467 * mrSges(4,1) + t1008 * (-t131 * t660 + t132 * t662) + t974 * t247 + m(7) * (t190 * t633 + t528 * t67 + t530 * t68) + t975 * t248 + t249 * t891 + t250 * t892 + t144 * t896 + t143 * t897 + t342 * t915 + t341 * t916) * qJD(3) + t23 * qJD(4) + t16 * qJD(5) + t7 * qJD(6), t821 + t11 * qJD(2) + t23 * qJD(3) + (-t341 * t983 + t342 * t696) * t877 + t79 * qJD(5) + t19 * qJD(6), qJD(2) * t14 + qJD(3) * t16 + qJD(4) * t79 + qJD(6) * t82 + t824, t827 + t5 * qJD(2) + t7 * qJD(3) + t19 * qJD(4) + t82 * qJD(5) + (-mrSges(7,1) * t61 - mrSges(7,2) * t60 + t140) * qJD(6); -qJD(3) * t2 - qJD(4) * t10 - qJD(5) * t13 + qJD(6) * t6 - t846, qJD(3) * t12 + qJD(4) * t35 + qJD(5) * t45 + qJD(6) * t29 (t626 * t734 + t624 * t735 + t974 * t988 - t975 * t575 + t637 * t753 + t636 * t754 + t723 + t708 - t174 * t831 - t175 * t832 + t541 * t802 - t543 * t803 - t655 * t534 + t633 * t362 + t474 * t551 + t530 * t425 + t528 * t427 + t1008 * (-t405 * t660 + t406 * t662) + m(7) * (t174 * t528 + t175 * t530 + t474 * t633) + t451 * t891 + t453 * t892 + t323 * t896 + t321 * t897 + t510 * t915 + t507 * t916 + t980 * pkin(9) + t987 * t901) * qJD(3) + t36 * qJD(4) + t39 * qJD(5) + t25 * qJD(6) + t721, t36 * qJD(3) + (-t507 * t983 + t510 * t696) * t877 + t187 * qJD(5) + t42 * qJD(6) + t719, qJD(3) * t39 + qJD(4) * t187 + qJD(6) * t200 - t718, t25 * qJD(3) + t42 * qJD(4) + t200 * qJD(5) + (-mrSges(7,1) * t171 - mrSges(7,2) * t170 + t365) * qJD(6) + t720; qJD(2) * t2 + qJD(4) * t22 + qJD(5) * t17 + qJD(6) * t8 - t851, -qJD(4) * t34 + qJD(5) * t40 + qJD(6) * t26 - t721, qJD(5) * t184 + qJD(6) * t78, -t717, t694 (-mrSges(7,1) * t530 - mrSges(7,2) * t528 + t552) * qJD(6) - t688; qJD(2) * t10 - qJD(3) * t22 + qJD(5) * t80 + qJD(6) * t20 - t821, qJD(3) * t34 - qJD(5) * t186 + qJD(6) * t41 - t719, t717, 0, t713, -t690 - t784; qJD(2) * t13 - qJD(3) * t17 - qJD(4) * t80 - qJD(6) * t81 - t824, -qJD(3) * t40 + qJD(4) * t186 + qJD(6) * t199 + t718, -t694 + t784, -t713, 0, -t692; -qJD(2) * t6 - qJD(3) * t8 - qJD(4) * t20 + qJD(5) * t81 - t827, -qJD(3) * t26 - qJD(4) * t41 - qJD(5) * t199 - t720, -t549 * qJD(5) + t688, t690, t692, 0;];
Cq  = t18;
