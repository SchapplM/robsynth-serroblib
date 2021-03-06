% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRPR9_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR9_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_coriolismatJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR9_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR9_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR9_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:27:21
% EndTime: 2019-03-09 05:28:03
% DurationCPUTime: 27.58s
% Computational Cost: add. (93762->1042), mult. (254608->1492), div. (0->0), fcn. (297970->14), ass. (0->566)
t619 = sin(qJ(3));
t616 = sin(pkin(7));
t847 = sin(pkin(6));
t849 = cos(pkin(12));
t690 = t849 * t847;
t850 = cos(pkin(7));
t851 = cos(pkin(6));
t637 = t851 * t616 + t850 * t690;
t615 = sin(pkin(12));
t724 = t615 * t847;
t906 = cos(qJ(3));
t526 = t619 * t637 + t724 * t906;
t567 = -t616 * t690 + t850 * t851;
t618 = sin(qJ(4));
t621 = cos(qJ(4));
t447 = -t526 * t618 + t567 * t621;
t448 = t526 * t621 + t567 * t618;
t846 = sin(pkin(13));
t848 = cos(pkin(13));
t332 = -t848 * t447 + t448 * t846;
t1038 = t332 * Ifges(6,4);
t645 = t447 * t846 + t448 * t848;
t1027 = t645 * Ifges(6,1);
t525 = t619 * t724 - t637 * t906;
t187 = t525 * Ifges(6,5) + t1027 - t1038;
t1047 = -Ifges(6,2) * t645 - t1038 + t187;
t617 = sin(qJ(6));
t913 = t617 / 0.2e1;
t721 = t848 * t618;
t606 = t621 * qJ(5);
t610 = t621 * pkin(10);
t794 = t610 + t606;
t897 = -qJ(5) - pkin(10);
t540 = -t897 * t721 + t794 * t846;
t811 = t616 * t619;
t569 = -t618 * t811 + t621 * t850;
t570 = t618 * t850 + t621 * t811;
t644 = t569 * t846 + t570 * t848;
t1046 = t540 * t644;
t1045 = t618 ^ 2 + t621 ^ 2;
t1028 = Ifges(6,4) * t645;
t620 = cos(qJ(6));
t274 = t525 * t620 - t617 * t645;
t275 = t525 * t617 + t620 * t645;
t110 = Ifges(7,5) * t275 + Ifges(7,6) * t274 + t332 * Ifges(7,3);
t1044 = -Ifges(6,1) * t332 - t1028 + t110;
t988 = m(7) / 0.4e1;
t1043 = 0.2e1 * t988;
t281 = -mrSges(6,2) * t525 - mrSges(6,3) * t332;
t1042 = t281 / 0.2e1;
t1041 = -t332 / 0.2e1;
t965 = -t332 / 0.4e1;
t1040 = t332 / 0.2e1;
t963 = t332 / 0.4e1;
t719 = t846 * t618;
t574 = -t621 * t848 + t719;
t398 = t525 * t574;
t344 = -t398 * t617 + t526 * t620;
t958 = t344 / 0.2e1;
t345 = t398 * t620 + t526 * t617;
t957 = t345 / 0.2e1;
t575 = -t621 * t846 - t721;
t397 = t525 * t575;
t1032 = t397 / 0.2e1;
t488 = -t848 * t569 + t570 * t846;
t1039 = t488 / 0.2e1;
t860 = t574 * mrSges(6,3);
t583 = -mrSges(5,1) * t621 + t618 * mrSges(5,2);
t1037 = t583 - mrSges(4,1);
t759 = t846 * pkin(4);
t603 = t759 + pkin(11);
t787 = 0.2e1 * t603;
t989 = m(7) / 0.2e1;
t1036 = t787 * t989 + mrSges(7,3);
t1007 = t897 * t719 + t794 * t848;
t904 = m(6) * t1007;
t1035 = -t860 + t904;
t152 = mrSges(7,1) * t275 + mrSges(7,2) * t274;
t853 = t620 * mrSges(7,2);
t856 = t617 * mrSges(7,1);
t584 = t853 + t856;
t760 = t848 * pkin(4);
t604 = -t760 - pkin(5);
t608 = Ifges(7,4) * t620;
t1021 = -Ifges(7,2) * t617 + t608;
t590 = Ifges(7,1) * t617 + t608;
t916 = t590 / 0.4e1;
t734 = t916 + t1021 / 0.4e1;
t766 = pkin(1) * t851;
t598 = t849 * t766;
t691 = t850 * t847;
t527 = t851 * pkin(2) + t598 + (-pkin(9) * t691 - qJ(2) * t847) * t615;
t698 = t616 * t724;
t553 = -pkin(1) * t847 - pkin(2) * t690 - pkin(9) * t698;
t428 = -t527 * t616 + t850 * t553;
t308 = pkin(3) * t525 - pkin(10) * t526 + t428;
t795 = qJ(2) * t690 + t615 * t766;
t516 = pkin(9) * t637 + t795;
t722 = t850 * t527;
t359 = t906 * t516 + (t553 * t616 + t722) * t619;
t319 = t567 * pkin(10) + t359;
t193 = t621 * t308 - t319 * t618;
t161 = -qJ(5) * t448 + t193;
t133 = pkin(4) * t525 + t161;
t194 = t618 * t308 + t319 * t621;
t162 = t447 * qJ(5) + t194;
t720 = t846 * t162;
t82 = t133 * t848 - t720;
t77 = -t525 * pkin(5) - t82;
t978 = -t77 / 0.2e1;
t1034 = -t604 * t152 / 0.2e1 + t584 * t978 - t734 * t274;
t898 = t618 * pkin(4);
t513 = -pkin(5) * t575 + pkin(11) * t574 + t898;
t419 = t513 * t620 + t540 * t617;
t420 = t513 * t617 - t540 * t620;
t765 = t616 * t906;
t454 = -t617 * t644 - t620 * t765;
t455 = -t617 * t765 + t620 * t644;
t676 = t1007 * t488 + t1046;
t1033 = t419 * t454 + t420 * t455 + t676;
t814 = t575 * t620;
t521 = mrSges(7,1) * t574 + mrSges(7,3) * t814;
t815 = t575 * t617;
t782 = mrSges(7,3) * t815;
t519 = -mrSges(7,2) * t574 + t782;
t800 = t620 * t519;
t660 = t521 * t913 - t800 / 0.2e1;
t1031 = t645 / 0.2e1;
t1030 = -t526 * Ifges(6,6) / 0.2e1;
t1029 = Ifges(5,3) + Ifges(6,3);
t611 = t617 ^ 2;
t613 = t620 ^ 2;
t793 = t611 + t613;
t680 = t793 * t787;
t507 = t584 * t575;
t977 = m(6) + m(7);
t1024 = t977 * t540 - t507;
t824 = t525 * t618;
t300 = -pkin(4) * t824 + t359;
t529 = mrSges(6,1) * t574 - mrSges(6,2) * t575;
t605 = -pkin(4) * t621 - pkin(3);
t903 = m(6) * t605;
t1023 = t529 + t903;
t607 = Ifges(7,5) * t620;
t879 = Ifges(7,6) * t617;
t1022 = t607 - t879;
t881 = Ifges(5,6) * t618;
t885 = Ifges(5,5) * t621;
t1020 = -Ifges(6,5) * t574 + Ifges(6,6) * t575 - t881 + t885;
t528 = -t575 * mrSges(6,1) - t574 * mrSges(6,2);
t895 = mrSges(5,2) * t621;
t585 = t618 * mrSges(5,1) + t895;
t1019 = -(t528 + t585) * t906 / 0.2e1;
t609 = Ifges(5,4) * t621;
t1012 = Ifges(5,2) * t618 - t609;
t592 = Ifges(5,1) * t618 + t609;
t1018 = -t1012 + t592;
t887 = Ifges(7,4) * t275;
t111 = Ifges(7,2) * t274 + t332 * Ifges(7,6) + t887;
t273 = Ifges(7,4) * t274;
t112 = Ifges(7,1) * t275 + t332 * Ifges(7,5) + t273;
t908 = t620 / 0.4e1;
t1017 = t112 * t908 - t617 * t111 / 0.4e1;
t1016 = Ifges(5,5) * t447 - Ifges(6,5) * t332 - Ifges(5,6) * t448 - Ifges(6,6) * t645;
t857 = t575 * mrSges(6,3);
t770 = -t857 / 0.2e1;
t938 = -t507 / 0.2e1;
t1011 = t938 + t770;
t1010 = t583 / 0.2e1 + t529 / 0.2e1 - mrSges(4,1) / 0.2e1;
t1009 = mrSges(7,3) * t793;
t792 = 0.2e1 * m(6);
t718 = t792 / 0.2e1;
t706 = pkin(4) * t718;
t1008 = t706 * t846 - mrSges(6,2);
t1006 = Ifges(7,5) * t957 + Ifges(7,6) * t958 + Ifges(7,3) * t1032;
t582 = -mrSges(7,1) * t620 + mrSges(7,2) * t617;
t995 = 0.2e1 * t604;
t1005 = -t706 * t848 + t989 * t995 - mrSges(6,1) + t582;
t669 = t615 * t691;
t555 = -t619 * t669 + t690 * t906;
t790 = 0.2e1 * pkin(10);
t493 = t621 * t555 + t618 * t698;
t825 = t493 * t621;
t492 = -t618 * t555 + t621 * t698;
t826 = t492 * t618;
t993 = m(5) / 0.4e1;
t1004 = -t555 * mrSges(4,2) / 0.2e1 + (t825 - t826) * t790 * t993 + (t825 / 0.2e1 - t826 / 0.2e1) * mrSges(5,3);
t512 = t574 * pkin(5) + t575 * pkin(11) + t605;
t414 = -t1007 * t617 + t512 * t620;
t717 = t792 / 0.4e1;
t771 = -t860 / 0.2e1;
t415 = t1007 * t620 + t512 * t617;
t835 = t415 * t620;
t1003 = t1007 * t717 + (-t414 * t617 + t835) * t1043 + t771 - t660;
t1002 = -0.2e1 * pkin(3);
t1001 = -0.2e1 * t82;
t157 = t848 * t162;
t83 = t846 * t133 + t157;
t1000 = 0.2e1 * t83;
t999 = 0.2e1 * t644;
t998 = 0.2e1 * t569;
t997 = 0.2e1 * t570;
t996 = -0.2e1 * t575;
t994 = m(5) / 0.2e1;
t992 = m(6) / 0.2e1;
t991 = m(6) / 0.4e1;
t990 = -m(7) / 0.4e1;
t987 = m(5) * pkin(3);
t986 = m(6) * pkin(4);
t984 = mrSges(5,1) / 0.2e1;
t983 = mrSges(7,1) / 0.2e1;
t982 = -mrSges(5,2) / 0.2e1;
t981 = -mrSges(6,2) / 0.2e1;
t980 = -mrSges(7,2) / 0.2e1;
t979 = mrSges(7,2) / 0.2e1;
t899 = t448 * pkin(4);
t170 = pkin(5) * t645 + pkin(11) * t332 + t899;
t87 = t161 * t848 - t720;
t66 = t170 * t620 - t617 * t87;
t976 = m(7) * t66;
t67 = t170 * t617 + t620 * t87;
t975 = m(7) * t67;
t156 = Ifges(7,1) * t274 - t887;
t974 = -t156 / 0.4e1;
t875 = t274 * mrSges(7,3);
t183 = -mrSges(7,2) * t332 + t875;
t973 = t183 / 0.2e1;
t186 = -t332 * Ifges(6,2) + t525 * Ifges(6,6) + t1028;
t972 = -t186 / 0.2e1;
t971 = t274 / 0.2e1;
t970 = -t275 / 0.4e1;
t969 = t275 / 0.2e1;
t282 = mrSges(6,1) * t525 - mrSges(6,3) * t645;
t968 = t282 / 0.2e1;
t956 = -t397 / 0.2e1;
t955 = t398 / 0.2e1;
t954 = -t414 / 0.2e1;
t953 = -t415 / 0.2e1;
t952 = -t419 / 0.2e1;
t950 = t447 / 0.2e1;
t949 = t448 / 0.2e1;
t948 = -t454 / 0.2e1;
t947 = t454 / 0.2e1;
t946 = -t455 / 0.2e1;
t945 = t455 / 0.2e1;
t944 = -t644 / 0.2e1;
t943 = -t488 / 0.2e1;
t940 = t644 / 0.2e1;
t694 = t906 * t846;
t695 = t906 * t848;
t552 = (-t618 * t694 + t621 * t695) * t616;
t500 = -t552 * t617 + t620 * t811;
t939 = t500 / 0.2e1;
t937 = -t519 / 0.2e1;
t936 = t519 / 0.2e1;
t935 = -t521 / 0.2e1;
t934 = t521 / 0.2e1;
t931 = t526 / 0.2e1;
t929 = t540 / 0.2e1;
t551 = (t618 * t695 + t621 * t694) * t616;
t928 = -t551 / 0.2e1;
t927 = -t567 / 0.2e1;
t926 = t569 / 0.2e1;
t925 = t570 / 0.2e1;
t924 = -t574 / 0.2e1;
t923 = t574 / 0.4e1;
t922 = -t575 / 0.2e1;
t920 = t582 / 0.2e1;
t886 = Ifges(7,4) * t617;
t587 = Ifges(7,2) * t620 + t886;
t918 = t587 / 0.4e1;
t915 = -t603 / 0.2e1;
t914 = t604 / 0.2e1;
t912 = t617 / 0.4e1;
t911 = t618 / 0.2e1;
t910 = -t620 / 0.2e1;
t909 = t620 / 0.2e1;
t907 = t621 / 0.2e1;
t905 = m(6) * t644;
t902 = m(7) * t414;
t901 = m(7) * t415;
t894 = mrSges(4,3) * t525;
t893 = mrSges(4,3) * t526;
t890 = Ifges(5,4) * t448;
t889 = Ifges(5,4) * t618;
t888 = Ifges(6,4) * t575;
t878 = Ifges(7,3) * t645;
t876 = Ifges(7,3) * t575;
t874 = t275 * mrSges(7,3);
t358 = -t619 * t516 + t553 * t765 + t722 * t906;
t425 = pkin(3) * t526 + pkin(10) * t525;
t251 = -t618 * t358 + t621 * t425;
t209 = t526 * pkin(4) + t525 * t606 + t251;
t252 = t621 * t358 + t618 * t425;
t221 = qJ(5) * t824 + t252;
t100 = t209 * t848 - t221 * t846;
t101 = t846 * t209 + t848 * t221;
t153 = -mrSges(7,1) * t274 + mrSges(7,2) * t275;
t166 = Ifges(7,4) * t345 + Ifges(7,2) * t344 + Ifges(7,6) * t397;
t167 = Ifges(7,1) * t345 + Ifges(7,4) * t344 + Ifges(7,5) * t397;
t184 = mrSges(7,1) * t332 - t874;
t212 = mrSges(6,1) * t332 + mrSges(6,2) * t645;
t239 = -mrSges(7,2) * t397 + mrSges(7,3) * t344;
t863 = t526 * Ifges(6,5);
t242 = Ifges(6,1) * t398 - Ifges(6,4) * t397 + t863;
t318 = -t567 * pkin(3) - t358;
t250 = -t447 * pkin(4) + t318;
t871 = t398 * mrSges(6,2);
t872 = t397 * mrSges(6,1);
t280 = t871 + t872;
t350 = -mrSges(6,2) * t526 - mrSges(6,3) * t397;
t351 = mrSges(6,1) * t526 - mrSges(6,3) * t398;
t356 = Ifges(5,6) * t526 + t1012 * t525;
t688 = Ifges(5,1) * t621 - t889;
t357 = Ifges(5,5) * t526 - t525 * t688;
t868 = t447 * mrSges(5,3);
t369 = -mrSges(5,2) * t525 + t868;
t867 = t448 * mrSges(5,3);
t370 = mrSges(5,1) * t525 - t867;
t403 = t585 * t525;
t456 = -mrSges(4,2) * t567 - t894;
t114 = t332 * pkin(5) - pkin(11) * t645 + t250;
t78 = pkin(11) * t525 + t83;
t51 = t114 * t620 - t617 * t78;
t52 = t114 * t617 + t620 * t78;
t668 = Ifges(6,5) * t955 + Ifges(6,6) * t956;
t343 = -mrSges(5,1) * t447 + mrSges(5,2) * t448;
t457 = mrSges(4,1) * t567 - t893;
t699 = m(5) * t318 + t343 - t457;
t173 = t397 * pkin(5) - t398 * pkin(11) + t300;
t94 = pkin(11) * t526 + t101;
t72 = t173 * t620 - t617 * t94;
t412 = -mrSges(5,2) * t526 + mrSges(5,3) * t824;
t726 = m(5) * t252 + t412;
t823 = t525 * t621;
t413 = t526 * mrSges(5,1) + mrSges(5,3) * t823;
t727 = m(5) * t251 + t413;
t240 = mrSges(7,1) * t397 - mrSges(7,3) * t345;
t728 = m(7) * t72 + t240;
t215 = -mrSges(7,1) * t344 + mrSges(7,2) * t345;
t93 = -t526 * pkin(5) - t100;
t729 = m(7) * t93 + t215;
t73 = t173 * t617 + t620 * t94;
t732 = m(6) * t82 + t282;
t733 = m(6) * t83 + t281;
t738 = -Ifges(6,4) * t398 / 0.2e1 + Ifges(6,2) * t1032 + t1030 + t1006;
t796 = -Ifges(4,5) * t525 - Ifges(4,6) * t526;
t442 = Ifges(5,4) * t447;
t297 = t448 * Ifges(5,1) + Ifges(5,5) * t525 + t442;
t799 = t621 * t297;
t296 = t447 * Ifges(5,2) + Ifges(5,6) * t525 + t890;
t805 = t618 * t296;
t3 = (t428 * mrSges(4,1) - t359 * mrSges(4,3) - Ifges(4,4) * t526 + Ifges(5,5) * t949 + Ifges(6,5) * t1031 + Ifges(4,6) * t927 + Ifges(5,6) * t950 + Ifges(6,6) * t1041) * t526 + (m(7) * t73 + t239) * t52 + (m(6) * t300 + t280) * t250 + t358 * t456 - t318 * t403 + t252 * t369 + t251 * t370 + t83 * t350 + t82 * t351 + t242 * t1031 + t300 * t212 + t72 * t184 + t73 * t183 + t93 * t153 + t699 * t359 + t738 * t332 + t728 * t51 + t729 * t77 + t732 * t100 + t733 * t101 + t726 * t194 + t727 * t193 + (t972 + t110 / 0.2e1) * t397 + t567 * t796 / 0.2e1 + t357 * t949 + t356 * t950 + t187 * t955 + t112 * t957 + t111 * t958 + t167 * t969 + t166 * t971 + (-t799 / 0.2e1 + t805 / 0.2e1 + t358 * mrSges(4,3) + Ifges(4,5) * t927 - t428 * mrSges(4,2) + (Ifges(4,2) - Ifges(4,1) + t1029) * t526 + t668 + (Ifges(4,4) - t885 / 0.2e1 + t881 / 0.2e1) * t525) * t525;
t873 = t3 * qJD(1);
t124 = Ifges(7,6) * t645 - t1021 * t332;
t591 = Ifges(7,1) * t620 - t886;
t125 = Ifges(7,5) * t645 - t332 * t591;
t208 = t584 * t332;
t839 = t332 * t617;
t210 = -mrSges(7,2) * t645 + mrSges(7,3) * t839;
t838 = t332 * t620;
t211 = mrSges(7,1) * t645 + mrSges(7,3) * t838;
t641 = -t1022 * t332 + t878;
t657 = m(7) * t77 + t153 - t732;
t714 = mrSges(6,1) * t645 - mrSges(6,2) * t332;
t725 = m(6) * t250 + t212;
t757 = -t838 / 0.2e1;
t758 = t839 / 0.2e1;
t86 = t161 * t846 + t157;
t4 = (t332 * t82 - t645 * t83) * mrSges(6,3) + t112 * t757 + t111 * t758 + t193 * t369 - t194 * t370 - t77 * t208 + t66 * t184 + t67 * t183 + (t442 / 0.2e1 - t193 * mrSges(5,3) + t297 / 0.2e1 + t318 * mrSges(5,2)) * t447 + t733 * t87 + t250 * t714 + t657 * t86 + (t210 + t975) * t52 + (t211 + t976) * t51 + (-t194 * mrSges(5,3) + t318 * mrSges(5,1) - t296 / 0.2e1 - t890 / 0.2e1 + (-Ifges(5,2) / 0.2e1 + Ifges(5,1) / 0.2e1) * t447 + t725 * pkin(4)) * t448 + t641 * t1040 + t125 * t969 + t124 * t971 + t645 * t972 + t1016 * t525 / 0.2e1 + t1044 * t1031 + t1047 * t1041;
t870 = t4 * qJD(1);
t869 = t415 * mrSges(7,3);
t866 = t488 * t86;
t865 = t51 * t617;
t864 = t52 * t620;
t861 = t540 * t86;
t859 = t574 * Ifges(7,5);
t858 = t574 * Ifges(7,6);
t730 = m(7) * t51 + t184;
t731 = m(7) * t52 + t183;
t17 = t657 * t645 - (-t617 * t730 + t620 * t731 + t733) * t332;
t845 = qJD(1) * t17;
t646 = t906 * t669;
t829 = t455 * t620;
t830 = t454 * t617;
t661 = -t829 / 0.2e1 + t830 / 0.2e1;
t670 = t619 * t690;
t828 = t488 * t645;
t756 = t828 / 0.2e1;
t385 = t492 * t846 + t493 * t848;
t554 = t670 + t646;
t346 = -t385 * t617 + t554 * t620;
t803 = t620 * t346;
t347 = t385 * t620 + t554 * t617;
t809 = t617 * t347;
t55 = (t756 - t803 / 0.2e1 - t809 / 0.2e1 + t661 * t332) * m(7) + (t756 - t332 * t940 - t646 / 0.2e1 - t670 / 0.2e1) * m(6);
t844 = qJD(1) * t55;
t154 = Ifges(7,5) * t274 - Ifges(7,6) * t275;
t155 = -Ifges(7,2) * t275 + t273;
t11 = t77 * t152 + t154 * t1040 + t51 * t183 - t52 * t184 + (-t52 * mrSges(7,3) - t111 / 0.2e1 + t156 / 0.2e1) * t275 + (-t51 * mrSges(7,3) + t155 / 0.2e1 + t112 / 0.2e1) * t274;
t843 = t11 * qJD(1);
t384 = -t492 * t848 + t493 * t846;
t16 = (m(4) * t359 + t456) * t555 + (m(5) * t194 + t369) * t493 + t733 * t385 + t731 * t347 + (m(5) * t193 + t370) * t492 + t730 * t346 + (-m(3) * t598 - t851 * mrSges(3,1) + (m(4) * t428 + mrSges(4,1) * t525 + mrSges(4,2) * t526) * t616 + (m(3) * qJ(2) + mrSges(3,3)) * t724) * t724 + t657 * t384 + (-m(4) * t358 + t699 + t725) * t554 + (m(3) * t795 - mrSges(3,2) * t851 + mrSges(3,3) * t690) * t690;
t842 = t16 * qJD(1);
t634 = (t274 * t948 + t275 * t946) * mrSges(7,3) + t183 * t947 + t184 * t946 + t152 * t1039;
t666 = -t346 * mrSges(7,1) / 0.2e1 + t347 * t979;
t19 = t634 + t666;
t841 = t19 * qJD(1);
t840 = t645 * t540;
t837 = t384 * t488;
t836 = t384 * t540;
t827 = t488 * t584;
t820 = t540 * t551;
t819 = t574 * t644;
t818 = t574 * t617;
t817 = t574 * t620;
t816 = t575 * t488;
t813 = t603 * t617;
t812 = t603 * t620;
t810 = t617 * t184;
t462 = -t1021 * t575 + t858;
t808 = t617 * t462;
t501 = t552 * t620 + t617 * t811;
t807 = t617 * t501;
t804 = t620 * t183;
t464 = -t575 * t591 + t859;
t802 = t620 * t464;
t801 = t620 * t500;
t791 = 0.2e1 * m(7);
t788 = 0.2e1 * t574;
t786 = 0.2e1 * t616;
t785 = t986 / 0.2e1;
t784 = t989 + t992;
t781 = t904 / 0.2e1;
t780 = t488 * t977;
t776 = -t875 / 0.2e1;
t775 = mrSges(6,3) * t1031;
t774 = mrSges(6,3) * t1040;
t773 = t867 / 0.2e1;
t769 = t857 / 0.2e1;
t768 = -t617 * mrSges(7,3) / 0.2e1;
t767 = mrSges(7,3) * t909;
t764 = t618 * t906;
t763 = t621 * t906;
t761 = t906 * t554;
t753 = t824 / 0.2e1;
t752 = -t823 / 0.2e1;
t751 = t818 / 0.2e1;
t750 = -t817 / 0.2e1;
t749 = t816 / 0.2e1;
t748 = -t813 / 0.2e1;
t745 = t811 / 0.2e1;
t743 = t810 / 0.2e1;
t740 = -t804 / 0.2e1;
t509 = t590 * t575;
t737 = -t462 / 0.4e1 + t509 / 0.4e1;
t508 = t575 * t587;
t736 = t464 / 0.4e1 + t508 / 0.4e1;
t735 = t918 - t591 / 0.4e1;
t716 = t791 / 0.2e1;
t715 = t791 / 0.4e1;
t711 = t448 * t785;
t710 = pkin(4) * t765;
t708 = mrSges(6,3) * t760;
t707 = mrSges(6,3) * t759;
t705 = pkin(4) * t717;
t704 = -t765 / 0.2e1;
t703 = t765 / 0.2e1;
t702 = -t764 / 0.2e1;
t697 = t936 + t901 / 0.2e1;
t696 = t902 / 0.2e1 + t934;
t693 = t973 + t776;
t692 = t874 / 0.2e1 + t184 / 0.2e1;
t689 = mrSges(5,1) * t448 + mrSges(5,2) * t447;
t586 = Ifges(7,5) * t617 + Ifges(7,6) * t620;
t681 = -t864 + t865;
t624 = 0.2e1 * (t454 * t66 + t455 * t67 + t488 * t681 + t644 * t77 + t866) * t990 + t211 * t948 + t210 * t946 + t153 * t944 + t644 * t968 - t208 * t943 + t281 * t1039 - t569 * t369 / 0.2e1 + t370 * t925 - (-t448 * t710 - t488 * t83 + t866 + (-t82 + t87) * t644) * t792 / 0.4e1 + t488 * t774 + t644 * t775 + t868 * t926 + t570 * t773 - t488 * t743 - t488 * t740 + (t714 + t689) * t703;
t629 = (-t346 * t617 + t347 * t620) * t787 * t988 + t492 * t984 + t493 * t982 + t346 * t768 + t347 * t767 + (t705 * t846 + t981) * t385 + (t995 * t988 - mrSges(6,1) / 0.2e1 + t920 - t848 * t705) * t384;
t10 = t624 + t629;
t677 = -t829 + t830;
t652 = t677 * t990;
t92 = t644 * t780 - (0.4e1 * t652 + t905) * t488;
t679 = -t10 * qJD(1) + t92 * qJD(2);
t581 = -0.4e1 * t616 ^ 2 * t906 * t619;
t107 = (0.4e1 * t552 * t644 + t581) * t991 + t551 * t780 + 0.4e1 * (t454 * t500 + t455 * t501) * t988 + (t581 + 0.4e1 * (-t569 * t764 + t570 * t763) * t616) * t993;
t675 = t616 * t702;
t623 = t551 * t153 / 0.2e1 + t552 * t1042 + (t454 * t72 + t455 * t73 + t488 * t93 + t500 * t51 + t501 * t52 + t551 * t77) * t715 + (-0.2e1 * t488 * t100 + t101 * t999 + t551 * t1001 + t552 * t1000 + (t250 * t619 - t300 * t906) * t786) * t991 + (t251 * t998 + t252 * t997 + (-t193 * t764 + t194 * t763 + t318 * t619 - t359 * t906) * t786) * t993 + t850 * (mrSges(4,1) * t526 - mrSges(4,2) * t525) / 0.2e1 + t370 * t675 + t412 * t925 + t413 * t926 + t282 * t928 + t184 * t939 + t350 * t940 + t215 * t1039 + t351 * t943 + t239 * t945 + t240 * t947 + t501 * t973 - (t457 + t893) * t811 / 0.2e1 + (t343 + t212) * t745 + (-t403 + t280) * t704 + (t369 * t621 + t456 + t894) * t703;
t8 = -t623 + (t346 * t414 + t347 * t415 + t836) * t715 + (t1007 * t385 + t836) * t717 + t385 * t771 + t347 * t936 + t346 * t934 + t1011 * t384 + (t1002 * t993 + t605 * t717 + t1010) * t554 + t1004;
t678 = -t8 * qJD(1) + t107 * qJD(2);
t631 = (-t574 * t680 + t604 * t996) * t988 + t582 * t922 + (-t574 * t846 + t575 * t848) * t705 + t924 * t1009;
t518 = mrSges(7,2) * t575 + mrSges(7,3) * t818;
t520 = -mrSges(7,1) * t575 + mrSges(7,3) * t817;
t633 = t518 * t913 + t520 * t909 + t618 * t785 + (t419 * t620 + t420 * t617) * t715;
t113 = -t528 + t631 - t633;
t632 = (-t332 * t680 + t645 * t995) * t988 + t645 * t920 + (-t332 * t846 - t645 * t848) * t705 + t1041 * t1009;
t25 = t210 * t913 + t211 * t909 + (t617 * t67 + t620 * t66) * t715 + t711 - t632 + t714;
t674 = qJD(1) * t25 - qJD(3) * t113;
t663 = -t853 / 0.2e1 - t856 / 0.2e1;
t654 = t663 * t574;
t191 = -t654 + t660;
t28 = (mrSges(7,2) * t1040 - t693) * t620 + (mrSges(7,1) * t1040 + t692) * t617;
t673 = t28 * qJD(1) + t191 * qJD(3);
t672 = (t613 / 0.2e1 + t611 / 0.2e1) * t603 * mrSges(7,3);
t667 = t607 / 0.2e1 - t879 / 0.2e1;
t665 = mrSges(7,1) * t952 + t420 * t979;
t664 = mrSges(7,1) * t939 + t501 * t980;
t659 = t519 * t915 + t737;
t658 = t521 * t915 + t736;
t653 = t663 * t488;
t460 = t574 * Ifges(7,3) - t1022 * t575;
t461 = -Ifges(7,6) * t575 - t1021 * t574;
t463 = -Ifges(7,5) * t575 - t574 * t591;
t506 = t584 * t574;
t530 = -t574 * Ifges(6,2) - t888;
t573 = Ifges(6,4) * t574;
t531 = -t575 * Ifges(6,1) - t573;
t589 = t621 * Ifges(5,2) + t889;
t622 = t540 * t1042 - m(6) * (0.2e1 * t861 - t540 * t1000 + 0.2e1 * (t250 * t618 + t448 * t605) * pkin(4) + (t1001 + 0.2e1 * t87) * t1007) / 0.4e1 + (-t868 / 0.2e1 + t369 / 0.2e1) * pkin(10) * t618 + t805 / 0.4e1 - t799 / 0.4e1 + (-t1022 * t574 + t808 - t876) * t965 - t621 * (-Ifges(5,2) * t448 + t442) / 0.4e1 - t318 * t585 / 0.2e1 + t448 * t589 / 0.4e1 - t250 * t528 / 0.2e1 - t52 * t518 / 0.2e1 - t51 * t520 / 0.2e1 - t274 * t461 / 0.4e1 - t420 * t183 / 0.2e1 + (t208 / 0.2e1 + t774) * t540 - t1007 * t153 / 0.2e1 - (t1007 * t77 + t414 * t66 + t415 * t67 + t419 * t51 + t420 * t52 + t861) * t791 / 0.4e1 + t1007 * t968 + t1007 * t775 + pkin(3) * t689 / 0.2e1 - t448 * t688 / 0.4e1 + (t1027 / 0.4e1 - t641 / 0.4e1 + t1017) * t574 + (t530 / 0.4e1 - t460 / 0.4e1 - t888 / 0.4e1) * t645 + (t802 + t531 - t573) * t963 - t605 * t714 / 0.2e1 + t83 * t770 + t82 * t771 + (-t186 / 0.4e1 + Ifges(6,2) * t963 + t1044 / 0.4e1) * t575 + t1047 * t923 - t529 * t899 / 0.2e1 - t212 * t898 / 0.2e1 - t618 * (Ifges(5,1) * t447 - t890) / 0.4e1 + t87 * t860 / 0.2e1 + t125 * t814 / 0.4e1 - t124 * t815 / 0.4e1 - t1018 * t447 / 0.4e1 - t1020 * t525 / 0.4e1 + t66 * t935 + t67 * t937 + t184 * t952 + t210 * t953 + t211 * t954 + t463 * t970 - t506 * t978 + (t507 / 0.2e1 + t769) * t86 + (t370 / 0.2e1 + t773) * t610;
t625 = t668 + Ifges(5,5) * t752 + t239 * t812 / 0.2e1 + t73 * t767 + Ifges(5,6) * t753 + t240 * t748 + t72 * t768 + t350 * t759 / 0.2e1 + t351 * t760 / 0.2e1 + t166 * t908 + t167 * t912 + t215 * t914 + t397 * t586 / 0.4e1 + t344 * t918 + t345 * t916 + t93 * t920 + (t100 * t848 + t101 * t846) * t705 + t251 * t984 + t252 * t982 + t100 * mrSges(6,1) / 0.2e1 + t101 * t981 + (t93 * t995 + (-t72 * t617 + t73 * t620) * t787) * t988 + t1029 * t931;
t2 = t625 + t622;
t628 = (t551 * t995 + (-t500 * t617 + t501 * t620) * t787) * t988 + mrSges(6,1) * t928 + t551 * t920 + t552 * t981 + t500 * t768 + t501 * t767 + (-t551 * t848 + t552 * t846) * t705 + mrSges(5,1) * t675 + t704 * t895;
t636 = -t1039 * t506 + t518 * t945 + t520 * t947 + t644 * t938;
t649 = t676 - t1046;
t37 = (-t618 * t710 + t649) * t717 + t1033 * t1043 - t628 + t636 + t1019 * t616 + (t770 + t769) * t644 + (t771 - t1003) * t488;
t39 = t688 * t911 - t618 * t589 / 0.2e1 + t605 * t528 - pkin(3) * t585 + t420 * t519 + t419 * t521 + (t463 * t910 + t461 * t913 - t460 / 0.2e1 + t530 / 0.2e1 - t888 / 0.2e1) * t575 + (-t802 / 0.2e1 + t808 / 0.2e1 - t531 / 0.2e1 + t573 / 0.2e1 - t667 * t574 + (-Ifges(7,3) / 0.2e1 + Ifges(6,1) / 0.2e1 - Ifges(6,2) / 0.2e1) * t575) * t574 + (m(7) * t420 + t518) * t415 + (m(7) * t419 + t520) * t414 + t1023 * t898 + t1024 * t1007 + t1018 * t907 + (-t904 - t506) * t540;
t651 = -t2 * qJD(1) + t37 * qJD(2) + t39 * qJD(3);
t505 = t575 * t582;
t58 = t415 * t521 - t540 * t505 + ((-t464 / 0.2e1 - t508 / 0.2e1 - t859 / 0.2e1) * t617 + (t509 / 0.2e1 - t462 / 0.2e1 - t869 - t858 / 0.2e1) * t620) * t575 + (-t519 + t782) * t414;
t627 = (mrSges(7,3) * t954 + t736) * t274 + (-t869 / 0.2e1 + t737) * t275 + (t620 * t974 + t111 * t908 + t586 * t963 + (-t865 / 0.2e1 + t864 / 0.2e1) * mrSges(7,3) + (t155 + t112) * t912) * t575 + t414 * t973 + t184 * t953 + t51 * t936 + t52 * t935 + t152 * t929 + t154 * t923 + t77 * t505 / 0.2e1;
t635 = t72 * t983 + t73 * t980 + t1006;
t6 = t627 - t635;
t640 = t454 * t937 + t455 * t934 + t505 * t943;
t61 = mrSges(7,3) * t575 * t661 + t640 + t664;
t650 = t6 * qJD(1) - t61 * qJD(2) - t58 * qJD(3);
t115 = (t745 + t749 + t819 / 0.2e1) * m(6) + (t801 / 0.2e1 + t807 / 0.2e1 + t749 - t661 * t574) * m(7);
t626 = (t681 * t788 + t77 * t996 + 0.2e1 * t840) * t988 + t281 * t924 + t153 * t922 + t575 * t968 + (t575 * t82 + t840) * t717 + (-t717 * t83 + t740 + t743) * t574 + t1011 * t645 - t1003 * t332;
t630 = 0.2e1 * t300 * t991 + t872 / 0.2e1 + t871 / 0.2e1 + t239 * t913 + t240 * t909 + (t617 * t73 + t620 * t72) * t715;
t14 = -t626 + t630;
t96 = (-t857 + t1024) * t575 + (t800 + m(7) * t835 + (-t521 - t902) * t617 + t1035) * t574;
t648 = qJD(1) * t14 + qJD(2) * t115 + qJD(3) * t96;
t642 = t878 / 0.2e1 + t66 * t983 + t67 * t980;
t639 = t505 * t914 + t584 * t929 + t607 * t923;
t12 = t607 * t965 + t735 * t275 + (t974 + t111 / 0.4e1 + t693 * t603 + (t1040 + t963) * Ifges(7,6)) * t617 + (Ifges(7,5) * t1041 - t112 / 0.4e1 - t155 / 0.4e1 + t692 * t603) * t620 + t642 + t1034;
t122 = -t827 / 0.2e1 - t653;
t409 = t604 * t584 + (t590 / 0.2e1 + t1021 / 0.2e1) * t620 + (t591 / 0.2e1 - t587 / 0.2e1) * t617;
t57 = (Ifges(7,3) / 0.2e1 + t672) * t575 + (t859 / 0.2e1 + t735 * t575 + t658) * t620 + (-0.3e1 / 0.4e1 * t858 + t734 * t575 + t659) * t617 + t639 + t665;
t638 = t12 * qJD(1) + t122 * qJD(2) - t57 * qJD(3) - t409 * qJD(4);
t532 = -0.2e1 * t616 * t761;
t192 = -t654 - t660;
t123 = t827 / 0.2e1 - t653;
t117 = t631 + t633;
t116 = (-t816 - t819) * t717 + (t677 * t788 - 0.2e1 * t816) * t988 + m(6) * t745 + (t801 + t807) * t715;
t62 = -t640 + t664 + (t454 * t768 + t455 * t767) * t575;
t56 = -t876 / 0.2e1 + Ifges(7,5) * t750 + Ifges(7,6) * t751 + (-t858 / 0.4e1 + t659) * t617 + t658 * t620 + (t617 * t734 + t620 * t735 + t672) * t575 + t639 - t665;
t54 = (t803 + t809) * t715 + t554 * t717 + t784 * t828 - (0.2e1 * t652 + t905 / 0.2e1) * t332;
t36 = t1033 * t715 + t649 * t717 + ((t940 + t944) * t575 + (t943 + t1039) * t574) * mrSges(6,3) - (-t617 * t696 + t620 * t697 + t781) * t488 + (t702 * t986 + t1019) * t616 + t628 + t636;
t29 = t804 / 0.2e1 + t275 * t768 - t810 / 0.2e1 + t620 * t776 - t663 * t332;
t27 = t711 + (t211 / 0.2e1 + t976 / 0.2e1) * t620 + (t210 / 0.2e1 + t975 / 0.2e1) * t617 + t632;
t18 = t634 - t666;
t15 = t626 + t630;
t13 = t587 * t970 + t1022 * t963 + t275 * t591 / 0.4e1 + t156 * t912 + t155 * t908 + t183 * t748 + t813 * t875 / 0.2e1 - t667 * t332 + t642 - (t874 + t184) * t812 / 0.2e1 + t1017 - t1034;
t9 = -t624 + t629;
t7 = t623 + (t771 + t781) * t385 + t696 * t346 + (-t987 / 0.2e1 + t903 / 0.2e1 + t1010) * t554 + t697 * t347 + (t540 * t784 + t1011) * t384 + t1004;
t5 = t627 + t635;
t1 = t625 - t622;
t20 = [qJD(2) * t16 + qJD(3) * t3 + qJD(4) * t4 + qJD(5) * t17 + qJD(6) * t11, t842 + ((t346 * t454 + t347 * t455 + t837) * t716 + (t385 * t999 + t532 + 0.2e1 * t837) * t992 + (t492 * t998 + t493 * t997 + t532) * t994 + m(4) * (t555 * t619 + t669 - t761) * t786 / 0.2e1) * qJD(2) + t7 * qJD(3) + t9 * qJD(4) + t54 * qJD(5) + t18 * qJD(6), t7 * qJD(2) + t1 * qJD(4) + t15 * qJD(5) + t5 * qJD(6) + t873 + (t796 + t605 * t280 + t72 * t521 - t93 * t507 + t415 * t239 + pkin(3) * t403 + (t519 + t901) * t73 - t358 * mrSges(4,2) + t592 * t752 + t589 * t753 + t1007 * t350 + (t167 * t910 + t166 * t913 + t100 * mrSges(6,3) - t863 / 0.2e1 - t242 / 0.2e1) * t575 + t728 * t414 + (-m(6) * t100 - t351 + t729) * t540 + (-t618 * t727 + t621 * t726) * pkin(10) + (-t987 + t1037) * t359 + t1035 * t101 + (-t251 * t618 + t252 * t621) * mrSges(5,3) + t1023 * t300 + t356 * t907 + t357 * t911 + (Ifges(5,5) * t618 + Ifges(5,6) * t621) * t931 + t531 * t955 + t530 * t956 + t464 * t957 + t462 * t958 + t460 * t1032 + (t1030 + t738) * t574) * qJD(3), t870 + t9 * qJD(2) + t1 * qJD(3) + (-t194 * mrSges(5,1) - t193 * mrSges(5,2) + t1005 * t86 + t1008 * t87 + t124 * t909 + t125 * t913 - t604 * t208 + t210 * t812 - t211 * t813 + t332 * t708 + t586 * t1031 + t587 * t758 + t590 * t757 - t645 * t707 + t1016 + t1036 * (-t617 * t66 + t620 * t67)) * qJD(4) + t27 * qJD(5) + t13 * qJD(6), qJD(2) * t54 + qJD(3) * t15 + qJD(4) * t27 + qJD(6) * t29 + t845, t843 + t18 * qJD(2) + t5 * qJD(3) + t13 * qJD(4) + t29 * qJD(5) + (-mrSges(7,1) * t52 - mrSges(7,2) * t51 + t154) * qJD(6); -qJD(3) * t8 - qJD(4) * t10 + qJD(5) * t55 + qJD(6) * t19 - t842, t107 * qJD(3) + t92 * qJD(4), t36 * qJD(4) + t116 * qJD(5) + t62 * qJD(6) + t678 + ((t1007 * t552 + t820) * t718 + (t414 * t500 + t415 * t501 + t820) * t716 + t501 * t519 + t500 * t521 - t552 * t860 + (t1045 * t790 * t906 + t619 * t1002) * t616 * t994 + (t605 * t718 + t1037 + t529) * t811 + (-t507 - t857) * t551 + (mrSges(5,3) * t1045 - mrSges(4,2)) * t765) * qJD(3), t36 * qJD(3) + (-t570 * mrSges(5,1) - t569 * mrSges(5,2) + t1005 * t644 + (-t680 * t989 - t1008 - t1009) * t488) * qJD(4) + t123 * qJD(6) + t679, qJD(3) * t116 + t844, t841 + t62 * qJD(3) + t123 * qJD(4) + (-mrSges(7,1) * t455 - mrSges(7,2) * t454) * qJD(6); qJD(2) * t8 - qJD(4) * t2 - qJD(5) * t14 + qJD(6) * t6 - t873, qJD(4) * t37 - qJD(5) * t115 - qJD(6) * t61 - t678, qJD(4) * t39 - qJD(5) * t96 - qJD(6) * t58 (pkin(10) * t583 + t1005 * t1007 - t1008 * t540 + t461 * t909 + t463 * t913 - t604 * t506 + t518 * t812 - t520 * t813 + t574 * t708 + t575 * t707 + t586 * t922 + t587 * t751 + t590 * t750 + t1020 + t1036 * (-t419 * t617 + t420 * t620)) * qJD(4) + t117 * qJD(5) + t56 * qJD(6) + t651, qJD(4) * t117 + qJD(6) * t192 - t648, t56 * qJD(4) + t192 * qJD(5) + (-mrSges(7,1) * t415 - mrSges(7,2) * t414 + t575 * t586) * qJD(6) + t650; qJD(2) * t10 + qJD(3) * t2 - qJD(5) * t25 - qJD(6) * t12 - t870, -qJD(3) * t37 - qJD(6) * t122 - t679, qJD(5) * t113 + qJD(6) * t57 - t651, t409 * qJD(6), -t674 (t582 * t603 + t1022) * qJD(6) - t638; -qJD(2) * t55 + qJD(3) * t14 + qJD(4) * t25 - qJD(6) * t28 - t845, qJD(3) * t115 - t844, -qJD(4) * t113 - qJD(6) * t191 + t648, t674, 0, -t584 * qJD(6) - t673; -qJD(2) * t19 - qJD(3) * t6 + qJD(4) * t12 + qJD(5) * t28 - t843, t61 * qJD(3) + t122 * qJD(4) - t841, -qJD(4) * t57 + qJD(5) * t191 - t650, t638, t673, 0;];
Cq  = t20;
