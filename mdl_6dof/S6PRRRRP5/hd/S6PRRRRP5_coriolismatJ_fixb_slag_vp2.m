% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRRRRP5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP5_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_coriolismatJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP5_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP5_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP5_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:20:31
% EndTime: 2019-03-09 00:21:11
% DurationCPUTime: 23.18s
% Computational Cost: add. (32223->1205), mult. (86035->1646), div. (0->0), fcn. (93239->12), ass. (0->577)
t1082 = -m(7) / 0.2e1;
t1084 = -m(6) / 0.2e1;
t1110 = mrSges(7,3) + mrSges(6,3);
t731 = sin(pkin(7));
t927 = sin(pkin(6));
t928 = cos(pkin(7));
t810 = t928 * t927;
t929 = cos(pkin(6));
t996 = cos(qJ(2));
t1097 = t929 * t731 + t996 * t810;
t734 = sin(qJ(3));
t735 = sin(qJ(2));
t831 = t735 * t927;
t995 = cos(qJ(3));
t500 = -t1097 * t995 + t734 * t831;
t813 = t995 * t927;
t501 = t1097 * t734 + t735 * t813;
t736 = cos(qJ(5));
t732 = sin(qJ(5));
t737 = cos(qJ(4));
t903 = t732 * t737;
t250 = t500 * t903 + t501 * t736;
t900 = t736 * t737;
t251 = -t500 * t900 + t501 * t732;
t981 = -qJ(6) - pkin(11);
t670 = t981 * t732;
t674 = t981 * t736;
t715 = -pkin(5) * t736 - pkin(4);
t733 = sin(qJ(4));
t909 = t500 * t733;
t917 = t251 * t736;
t918 = t250 * t732;
t1120 = -t1110 * (t917 / 0.2e1 - t918 / 0.2e1) + (pkin(4) * t909 + (t917 - t918) * pkin(11)) * t1084 + (t250 * t670 - t251 * t674 - t715 * t909) * t1082;
t1076 = mrSges(5,2) / 0.2e1;
t790 = t735 * t810;
t556 = -t734 * t790 + t813 * t996;
t816 = t731 * t831;
t458 = t556 * t737 + t733 * t816;
t814 = t996 * t927;
t555 = t734 * t814 + t790 * t995;
t292 = -t458 * t732 + t555 * t736;
t293 = t458 * t736 + t555 * t732;
t457 = t556 * t733 - t737 * t816;
t914 = t293 * t736;
t915 = t292 * t732;
t1119 = t1110 * (t914 / 0.2e1 - t915 / 0.2e1) - (-pkin(4) * t457 + (t914 - t915) * pkin(11)) * t1084 - (t292 * t670 - t293 * t674 + t457 * t715) * t1082 - t458 * t1076;
t1083 = m(6) / 0.2e1;
t1081 = m(7) / 0.2e1;
t873 = pkin(2) * t928;
t906 = t731 * t734;
t619 = -pkin(9) * t906 + t873 * t995;
t568 = -pkin(3) * t928 - t619;
t615 = t733 * t906 - t737 * t928;
t616 = t733 * t928 + t737 * t906;
t356 = t615 * pkin(4) - t616 * pkin(11) + t568;
t825 = t734 * t873;
t872 = t731 * t995;
t621 = pkin(9) * t872 + t825;
t569 = pkin(10) * t928 + t621;
t570 = (-pkin(3) * t995 - pkin(10) * t734 - pkin(2)) * t731;
t392 = t737 * t569 + t733 * t570;
t363 = -pkin(11) * t872 + t392;
t139 = t356 * t732 + t363 * t736;
t510 = -t732 * t616 - t736 * t872;
t95 = qJ(6) * t510 + t139;
t1118 = m(7) * t95;
t511 = t616 * t736 - t732 * t872;
t308 = -mrSges(7,1) * t510 + mrSges(7,2) * t511;
t391 = -t733 * t569 + t570 * t737;
t362 = pkin(4) * t872 - t391;
t207 = -t510 * pkin(5) + t362;
t994 = m(7) * t207;
t1117 = t308 + t994;
t468 = mrSges(5,1) * t615 + mrSges(5,2) * t616;
t887 = mrSges(4,3) * t906;
t649 = mrSges(4,1) * t928 - t887;
t1114 = m(5) * t568 + t468 - t649;
t669 = -pkin(4) * t737 - t733 * pkin(11) - pkin(3);
t648 = t736 * t669;
t902 = t733 * t736;
t811 = -qJ(6) * t902 + t648;
t495 = (-pkin(10) * t732 - pkin(5)) * t737 + t811;
t560 = pkin(10) * t900 + t732 * t669;
t904 = t732 * t733;
t515 = -qJ(6) * t904 + t560;
t891 = pkin(10) * t903;
t559 = t648 - t891;
t990 = pkin(5) * t732;
t874 = pkin(10) + t990;
t663 = t874 * t733;
t934 = t733 * mrSges(5,2);
t673 = -mrSges(5,1) * t737 + t934;
t658 = -mrSges(6,1) * t737 - mrSges(6,3) * t902;
t1017 = t658 / 0.2e1;
t933 = t737 * mrSges(7,1);
t657 = -mrSges(7,3) * t902 - t933;
t1019 = t657 / 0.2e1;
t838 = t1017 + t1019;
t911 = t458 * t737;
t912 = t457 * t733;
t1113 = -(t673 / 0.2e1 - mrSges(4,1) / 0.2e1) * t555 - t838 * t292 - m(5) * (-pkin(3) * t555 + (t911 + t912) * pkin(10)) / 0.2e1 + (pkin(10) * t912 + t292 * t559 + t293 * t560) * t1084 + (t292 * t495 + t293 * t515 + t457 * t663) * t1082 + t556 * mrSges(4,2) / 0.2e1 - mrSges(5,3) * t911 / 0.2e1;
t1112 = m(7) + m(6);
t1111 = mrSges(6,2) + mrSges(7,2);
t1109 = Ifges(6,5) + Ifges(7,5);
t1108 = Ifges(6,6) + Ifges(7,6);
t1107 = Ifges(6,3) + Ifges(7,3);
t1106 = -mrSges(4,1) + t673;
t1080 = m(7) * pkin(5);
t878 = -mrSges(7,1) - t1080;
t870 = t737 * t995;
t557 = (-t732 * t870 + t734 * t736) * t731;
t558 = (t732 * t734 + t736 * t870) * t731;
t871 = t733 * t995;
t826 = t731 * t871;
t344 = Ifges(7,4) * t558 + Ifges(7,2) * t557 + Ifges(7,6) * t826;
t345 = Ifges(6,4) * t558 + Ifges(6,2) * t557 + Ifges(6,6) * t826;
t1104 = t345 + t344;
t346 = Ifges(7,1) * t558 + Ifges(7,4) * t557 + Ifges(7,5) * t826;
t347 = Ifges(6,1) * t558 + Ifges(6,4) * t557 + Ifges(6,5) * t826;
t1103 = t347 + t346;
t309 = -mrSges(6,1) * t510 + mrSges(6,2) * t511;
t936 = t616 * mrSges(5,3);
t526 = -mrSges(5,1) * t872 - t936;
t1102 = -t526 + t309;
t678 = Ifges(7,5) * t732 + Ifges(7,6) * t736;
t680 = Ifges(6,5) * t732 + Ifges(6,6) * t736;
t1101 = t680 + t678;
t719 = t732 * mrSges(7,1);
t720 = t736 * mrSges(7,2);
t1100 = t720 + t719;
t724 = Ifges(7,4) * t736;
t1099 = -Ifges(7,2) * t732 + t724;
t688 = Ifges(7,1) * t732 + t724;
t725 = Ifges(6,4) * t736;
t1098 = -Ifges(6,2) * t732 + t725;
t690 = Ifges(6,1) * t732 + t725;
t1093 = t678 / 0.2e1 + t680 / 0.2e1;
t654 = mrSges(6,2) * t737 - mrSges(6,3) * t904;
t1023 = t654 / 0.2e1;
t931 = t737 * mrSges(7,2);
t653 = -mrSges(7,3) * t904 + t931;
t1025 = t653 / 0.2e1;
t1092 = t1025 + t1023;
t676 = t732 * mrSges(6,1) + t736 * mrSges(6,2);
t627 = t676 * t733;
t1029 = t627 / 0.2e1;
t626 = t1100 * t733;
t1030 = t626 / 0.2e1;
t840 = t1030 + t1029;
t970 = Ifges(7,4) * t732;
t689 = Ifges(7,1) * t736 - t970;
t583 = -Ifges(7,5) * t737 + t689 * t733;
t972 = Ifges(6,4) * t732;
t691 = Ifges(6,1) * t736 - t972;
t585 = -Ifges(6,5) * t737 + t691 * t733;
t841 = t585 / 0.2e1 + t583 / 0.2e1;
t1001 = -t732 / 0.2e1;
t1020 = -t657 / 0.2e1;
t939 = t615 * mrSges(7,1);
t948 = t511 * mrSges(7,3);
t383 = t939 - t948;
t1061 = -t383 / 0.2e1;
t938 = t615 * mrSges(7,2);
t950 = t510 * mrSges(7,3);
t381 = -t938 + t950;
t138 = t736 * t356 - t363 * t732;
t94 = -qJ(6) * t511 + t138;
t86 = pkin(5) * t615 + t94;
t1090 = (-t495 * t511 + t510 * t515 + (-t732 * t95 - t736 * t86) * t733) * t1081 + t510 * t1025 + t511 * t1020 + (t381 * t1001 + t736 * t1061) * t733;
t660 = t733 * mrSges(6,1) - mrSges(6,3) * t900;
t1015 = t660 / 0.2e1;
t659 = t733 * mrSges(7,1) - mrSges(7,3) * t900;
t1016 = t659 / 0.2e1;
t656 = -t733 * mrSges(6,2) - mrSges(6,3) * t903;
t1021 = t656 / 0.2e1;
t655 = -t733 * mrSges(7,2) - mrSges(7,3) * t903;
t1022 = t655 / 0.2e1;
t629 = t676 * t737;
t1027 = t629 / 0.2e1;
t628 = t1100 * t737;
t1028 = t628 / 0.2e1;
t768 = -t731 * t814 + t928 * t929;
t354 = t501 * t737 + t733 * t768;
t181 = -t354 * t732 + t500 * t736;
t182 = t354 * t736 + t500 * t732;
t353 = t501 * t733 - t737 * t768;
t697 = t733 * pkin(4) - pkin(11) * t737;
t563 = pkin(10) * t904 + t736 * t697;
t496 = t733 * pkin(5) - qJ(6) * t900 + t563;
t564 = -pkin(10) * t902 + t732 * t697;
t522 = -qJ(6) * t903 + t564;
t664 = t874 * t737;
t1024 = -t654 / 0.2e1;
t1026 = -t653 / 0.2e1;
t839 = t1026 + t1024;
t987 = pkin(10) * t737;
t988 = pkin(10) * t733;
t1089 = t840 * t354 + (t1022 + t1021) * t182 + (t1016 + t1015) * t181 + (t563 * t181 + t564 * t182 + t354 * t988) * t1083 + (t181 * t496 + t182 * t522 + t354 * t663) * t1081 + ((t559 * t732 - t560 * t736 + t987) * t1083 + (t495 * t732 - t515 * t736 + t664) * t1081 + t732 * t838 + t736 * t839 + t1027 + t1028) * t353;
t1088 = 0.2e1 * m(7);
t1087 = t731 ^ 2;
t1086 = 2 * qJD(3);
t1085 = m(5) / 0.2e1;
t1079 = -mrSges(5,1) / 0.2e1;
t1078 = mrSges(5,1) / 0.2e1;
t1077 = -mrSges(6,1) / 0.2e1;
t1075 = mrSges(6,2) / 0.2e1;
t1074 = mrSges(7,2) / 0.2e1;
t1073 = -mrSges(6,3) / 0.2e1;
t1072 = mrSges(6,3) / 0.2e1;
t1071 = mrSges(7,3) / 0.2e1;
t471 = pkin(4) * t616 + pkin(11) * t615;
t189 = -t391 * t732 + t736 * t471;
t907 = t615 * t736;
t115 = pkin(5) * t616 + qJ(6) * t907 + t189;
t1070 = -t115 / 0.2e1;
t620 = (pkin(3) * t734 - pkin(10) * t995) * t731;
t456 = t737 * t619 + t733 * t620;
t407 = pkin(11) * t906 + t456;
t476 = t825 + (pkin(4) * t871 + pkin(9) * t995 - pkin(11) * t870) * t731;
t199 = -t732 * t407 + t736 * t476;
t133 = pkin(5) * t826 - t558 * qJ(6) + t199;
t1069 = -t133 / 0.2e1;
t307 = mrSges(6,1) * t511 + mrSges(6,2) * t510;
t1066 = t307 / 0.2e1;
t1065 = t308 / 0.2e1;
t1064 = t354 / 0.2e1;
t1063 = -t381 / 0.2e1;
t951 = t510 * mrSges(6,3);
t382 = -mrSges(6,2) * t615 + t951;
t1062 = -t382 / 0.2e1;
t1060 = t383 / 0.2e1;
t949 = t511 * mrSges(6,3);
t384 = mrSges(6,1) * t615 - t949;
t1059 = -t384 / 0.2e1;
t1058 = t384 / 0.2e1;
t394 = -mrSges(6,1) * t557 + mrSges(6,2) * t558;
t1057 = t394 / 0.2e1;
t908 = t615 * t732;
t448 = -mrSges(7,2) * t616 + mrSges(7,3) * t908;
t1056 = t448 / 0.2e1;
t449 = -mrSges(6,2) * t616 + mrSges(6,3) * t908;
t1055 = t449 / 0.2e1;
t450 = mrSges(7,1) * t616 + mrSges(7,3) * t907;
t1054 = t450 / 0.2e1;
t1053 = t457 / 0.2e1;
t486 = mrSges(7,1) * t826 - t558 * mrSges(7,3);
t1052 = -t486 / 0.2e1;
t487 = mrSges(6,1) * t826 - t558 * mrSges(6,3);
t1051 = t487 / 0.2e1;
t1050 = -t496 / 0.2e1;
t1049 = t500 / 0.2e1;
t1047 = t511 / 0.2e1;
t969 = Ifges(5,5) * t734;
t1046 = (Ifges(5,1) * t870 - Ifges(5,4) * t871 + t969) * t731 / 0.2e1;
t514 = t811 - t891;
t1045 = t514 / 0.2e1;
t1044 = t522 / 0.2e1;
t937 = t615 * mrSges(5,3);
t525 = mrSges(5,2) * t872 - t937;
t1043 = -t525 / 0.2e1;
t1042 = t557 / 0.2e1;
t1038 = t559 / 0.2e1;
t1037 = t564 / 0.2e1;
t624 = mrSges(7,1) * t902 - mrSges(7,2) * t904;
t1033 = t624 / 0.2e1;
t980 = mrSges(6,1) * t736;
t809 = -mrSges(6,2) * t732 + t980;
t625 = t809 * t733;
t1032 = -t625 / 0.2e1;
t1031 = t625 / 0.2e1;
t1018 = -t658 / 0.2e1;
t1014 = t663 / 0.2e1;
t1013 = t664 / 0.2e1;
t1012 = t670 / 0.2e1;
t671 = -mrSges(7,1) * t736 + mrSges(7,2) * t732;
t1011 = -t671 / 0.2e1;
t1010 = t671 / 0.2e1;
t1009 = t809 / 0.2e1;
t1008 = t674 / 0.2e1;
t1007 = t1100 / 0.2e1;
t1006 = t676 / 0.2e1;
t932 = t737 * mrSges(5,2);
t677 = t733 * mrSges(5,1) + t932;
t1005 = t677 / 0.2e1;
t1002 = t715 / 0.2e1;
t999 = -t733 / 0.2e1;
t998 = t736 / 0.2e1;
t993 = m(7) * t663;
t992 = m(7) * t715;
t991 = m(7) * t733;
t989 = pkin(10) * t500;
t485 = -mrSges(6,2) * t826 + t557 * mrSges(6,3);
t986 = pkin(11) * t485;
t985 = pkin(11) * t732;
t984 = t86 * mrSges(7,3);
t983 = t95 * mrSges(7,3);
t979 = mrSges(7,3) * t736;
t976 = Ifges(4,4) * t734;
t975 = Ifges(5,4) * t616;
t974 = Ifges(5,4) * t733;
t973 = Ifges(6,4) * t511;
t971 = Ifges(7,4) * t511;
t968 = Ifges(6,5) * t558;
t722 = Ifges(6,5) * t736;
t967 = Ifges(7,5) * t558;
t721 = Ifges(7,5) * t736;
t964 = Ifges(6,6) * t557;
t963 = Ifges(6,6) * t732;
t962 = Ifges(7,6) * t557;
t961 = Ifges(7,6) * t732;
t960 = Ifges(5,3) * t734;
t959 = Ifges(6,3) * t616;
t958 = Ifges(6,3) * t733;
t957 = Ifges(7,3) * t616;
t956 = Ifges(7,3) * t733;
t955 = t138 * mrSges(6,3);
t954 = t139 * mrSges(6,3);
t953 = t495 * mrSges(7,3);
t952 = t495 * t95;
t947 = t514 * t95;
t946 = t515 * mrSges(7,3);
t945 = t515 * t86;
t944 = t515 * t94;
t943 = t557 * mrSges(7,1);
t942 = t558 * mrSges(7,2);
t941 = t559 * mrSges(6,3);
t940 = t560 * mrSges(6,3);
t935 = t732 * mrSges(7,3);
t930 = -t809 - mrSges(5,1);
t926 = t181 * t511;
t925 = t181 * t736;
t924 = t182 * t510;
t923 = t182 * t732;
t922 = t199 * t732;
t200 = t736 * t407 + t732 * t476;
t921 = t200 * t736;
t22 = (t181 * t250 + t182 * t251 - t353 * t909) * t1112 + m(5) * (-t353 * t733 - t354 * t737 + t501) * t500;
t920 = t22 * qJD(1);
t910 = t500 * t555;
t913 = t353 * t457;
t23 = m(5) * (t354 * t458 + t910 + t913) + m(4) * (t501 * t556 + t768 * t816 + t910) + (t181 * t292 + t182 * t293 + t913) * t1112;
t919 = t23 * qJD(1);
t901 = t736 * t182;
t905 = t732 * t181;
t798 = t901 - t905;
t29 = (t354 - t798) * t353 * t1112;
t916 = t29 * qJD(1);
t190 = t736 * t391 + t732 * t471;
t899 = -t495 + t514;
t310 = Ifges(7,5) * t510 - Ifges(7,6) * t511;
t311 = Ifges(6,5) * t510 - Ifges(6,6) * t511;
t898 = -Ifges(5,5) * t615 - Ifges(5,6) * t616;
t727 = t732 ^ 2;
t729 = t736 ^ 2;
t897 = t727 + t729;
t896 = qJD(4) * t732;
t895 = qJD(4) * t736;
t894 = t1080 / 0.2e1;
t892 = pkin(5) * t902;
t890 = -t995 / 0.2e1;
t889 = -t995 / 0.4e1;
t285 = -pkin(5) * t908 + t392;
t886 = t285 * t1081;
t885 = m(7) * t1013;
t884 = mrSges(7,1) / 0.2e1 + mrSges(6,1) / 0.2e1;
t883 = t1074 + t1075;
t882 = t1071 + t1072;
t881 = Ifges(6,5) / 0.2e1 + Ifges(7,5) / 0.2e1;
t880 = Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1;
t879 = -Ifges(6,3) / 0.2e1 - Ifges(7,3) / 0.2e1;
t877 = m(7) * (-t86 + t94);
t876 = t979 / 0.2e1;
t875 = mrSges(5,3) * t999;
t865 = -t909 / 0.2e1;
t862 = -t906 / 0.2e1;
t861 = t906 / 0.2e1;
t218 = Ifges(7,2) * t510 + Ifges(7,6) * t615 + t971;
t219 = Ifges(6,2) * t510 + Ifges(6,6) * t615 + t973;
t856 = t218 / 0.2e1 + t219 / 0.2e1;
t855 = -t219 / 0.4e1 - t218 / 0.4e1;
t507 = Ifges(7,4) * t510;
t220 = Ifges(7,1) * t511 + Ifges(7,5) * t615 + t507;
t508 = Ifges(6,4) * t510;
t221 = Ifges(6,1) * t511 + Ifges(6,5) * t615 + t508;
t854 = -t220 / 0.2e1 - t221 / 0.2e1;
t853 = t221 / 0.4e1 + t220 / 0.4e1;
t852 = -t309 / 0.2e1 + t526 / 0.2e1;
t337 = Ifges(7,6) * t616 - t1099 * t615;
t338 = Ifges(6,6) * t616 - t1098 * t615;
t851 = t337 / 0.2e1 + t338 / 0.2e1;
t339 = Ifges(7,5) * t616 - t615 * t689;
t340 = Ifges(6,5) * t616 - t615 * t691;
t850 = t339 / 0.2e1 + t340 / 0.2e1;
t849 = t1062 + t1063;
t848 = t1058 + t1060;
t420 = t676 * t615;
t847 = t1043 - t420 / 0.2e1;
t579 = -Ifges(7,6) * t737 + t1099 * t733;
t581 = -Ifges(6,6) * t737 + t1098 * t733;
t846 = -t579 / 0.2e1 - t581 / 0.2e1;
t580 = Ifges(7,6) * t733 + t1099 * t737;
t582 = Ifges(6,6) * t733 + t1098 * t737;
t845 = -t580 / 0.2e1 - t582 / 0.2e1;
t844 = t581 / 0.4e1 + t579 / 0.4e1;
t584 = Ifges(7,5) * t733 + t689 * t737;
t586 = Ifges(6,5) * t733 + t691 * t737;
t843 = t584 / 0.2e1 + t586 / 0.2e1;
t842 = -t585 / 0.4e1 - t583 / 0.4e1;
t837 = t722 / 0.4e1 - t963 / 0.4e1 + t721 / 0.4e1 - t961 / 0.4e1;
t682 = Ifges(7,2) * t736 + t970;
t684 = Ifges(6,2) * t736 + t972;
t836 = -t682 / 0.4e1 - t684 / 0.4e1;
t835 = t682 / 0.2e1 + t684 / 0.2e1;
t834 = -t688 / 0.4e1 - t690 / 0.4e1;
t833 = t688 / 0.2e1 + t690 / 0.2e1;
t832 = m(7) * t899;
t829 = m(7) * t892;
t455 = -t733 * t619 + t620 * t737;
t828 = -m(6) * pkin(4) + t930;
t827 = mrSges(4,3) * t872;
t822 = -t872 / 0.2e1;
t821 = t872 / 0.2e1;
t819 = t870 / 0.2e1;
t723 = Ifges(5,5) * t737;
t818 = t723 * t889;
t817 = t883 * t736;
t815 = t1087 * t831;
t812 = t1009 + t1011 + t1078;
t805 = t722 - t963;
t804 = t721 - t961;
t146 = qJ(6) * t557 + t200;
t406 = -pkin(4) * t906 - t455;
t284 = -t557 * pkin(5) + t406;
t393 = t942 - t943;
t484 = -mrSges(7,2) * t826 + t557 * mrSges(7,3);
t561 = (mrSges(5,1) * t871 + mrSges(5,2) * t870) * t731;
t592 = (-mrSges(5,2) * t734 - mrSges(5,3) * t871) * t731;
t593 = (mrSges(5,1) * t734 - mrSges(5,3) * t870) * t731;
t618 = (mrSges(4,1) * t734 + mrSges(4,2) * t995) * t731;
t650 = -mrSges(4,2) * t928 + t827;
t738 = -t849 * t251 + t848 * t250 + (t484 / 0.2e1 + t485 / 0.2e1) * t182 + (t486 / 0.2e1 + t1051) * t181 + (-t593 / 0.2e1 + t1057 + t393 / 0.2e1) * t353 + (-t649 / 0.2e1 + t468 / 0.2e1 + mrSges(4,3) * t862) * t501 + (-t650 / 0.2e1 + t561 / 0.2e1 + mrSges(4,3) * t821 + t737 * t1043 + (-t308 / 0.2e1 + t852) * t733) * t500 + (-t455 * t353 + t456 * t354 + t568 * t501 + (t391 * t733 - t392 * t737 + t621) * t500) * t1085 + (t138 * t250 + t139 * t251 + t181 * t199 + t182 * t200 + t353 * t406 - t362 * t909) * t1083 + (t133 * t181 + t146 * t182 - t207 * t909 + t250 * t86 + t251 * t95 + t284 * t353) * t1081 + t592 * t1064 + t768 * t618 / 0.2e1;
t4 = t738 + (t875 - t840) * t457 + t839 * t293 + t1113;
t216 = Ifges(7,5) * t511 + Ifges(7,6) * t510 + Ifges(7,3) * t615;
t217 = Ifges(6,5) * t511 + Ifges(6,6) * t510 + Ifges(6,3) * t615;
t342 = Ifges(7,3) * t826 + t962 + t967;
t343 = Ifges(6,3) * t826 + t964 + t968;
t408 = -Ifges(5,2) * t615 - Ifges(5,6) * t872 + t975;
t602 = Ifges(5,4) * t615;
t409 = Ifges(5,1) * t616 - Ifges(5,5) * t872 - t602;
t512 = (Ifges(5,4) * t870 - Ifges(5,2) * t871 + Ifges(5,6) * t734) * t731;
t707 = Ifges(4,5) * t872;
t793 = t731 * t819;
t795 = t733 * t821;
t796 = t733 * t822;
t7 = t409 * t793 + t408 * t796 + (Ifges(4,6) * t928 + (Ifges(4,2) * t995 + t976) * t731) * t862 + (t217 + t216) * t795 + t928 * (-Ifges(4,6) * t906 + t707) / 0.2e1 + (-t827 + t650) * t619 + m(5) * (t391 * t455 + t392 * t456) + (Ifges(5,5) * t616 - Ifges(5,6) * t615 - Ifges(5,3) * t872) * t861 + (t734 * (Ifges(4,1) * t995 - t976) / 0.2e1 + (Ifges(5,5) * t870 - Ifges(5,6) * t871 + t960) * t890) * t1087 + (t219 + t218) * t1042 + t1103 * t1047 + t1104 * t510 / 0.2e1 - t615 * t512 / 0.2e1 + t392 * t592 + t391 * t593 + t568 * t561 + t456 * t525 + t455 * t526 + t86 * t486 + t138 * t487 + t95 * t484 + t139 * t485 + t406 * t309 + t362 * t394 + t207 * t393 + t146 * t381 + t200 * t382 + t133 * t383 + t199 * t384 + (Ifges(4,5) * t928 + 0.2e1 * Ifges(4,4) * t872 + (Ifges(4,1) - Ifges(4,2)) * t906) * t821 + (t343 + t342) * t615 / 0.2e1 + (t221 + t220) * t558 / 0.2e1 + t284 * t308 + (t1114 - t887) * t621 + m(7) * (t133 * t86 + t146 * t95 + t207 * t284) + m(6) * (t138 * t199 + t139 * t200 + t362 * t406) + t616 * t1046 - t731 * pkin(2) * t618;
t803 = t4 * qJD(1) + t7 * qJD(2);
t802 = t86 * t732 - t95 * t736;
t145 = qJ(6) * t908 + t190;
t335 = -t615 * t804 + t957;
t336 = -t615 * t805 + t959;
t419 = t1100 * t615;
t451 = mrSges(6,1) * t616 + mrSges(6,3) * t907;
t467 = mrSges(5,1) * t616 - mrSges(5,2) * t615;
t469 = -Ifges(5,2) * t616 - t602;
t470 = -Ifges(5,1) * t615 - t975;
t10 = t568 * t467 + t391 * t525 + t138 * t451 + t95 * t448 + t139 * t449 + t86 * t450 - t207 * t419 - t362 * t420 + t145 * t381 + t190 * t382 + t115 * t383 + t189 * t384 + t285 * t308 + t898 * t822 + t850 * t511 + t851 * t510 + t1102 * t392 + m(7) * (t115 * t86 + t145 * t95 + t207 * t285) + m(6) * (t138 * t189 + t139 * t190 + t362 * t392) + (t470 / 0.2e1 + t216 / 0.2e1 + t217 / 0.2e1 - t408 / 0.2e1 - t392 * mrSges(5,3)) * t616 + (-t469 / 0.2e1 - t409 / 0.2e1 + t335 / 0.2e1 + t336 / 0.2e1 + t391 * mrSges(5,3) + t854 * t736 + t856 * t732) * t615;
t786 = -t937 / 0.2e1 + t847;
t787 = -t936 / 0.2e1 - t852;
t740 = (t1056 + t1055) * t182 + (t1054 + t451 / 0.2e1) * t181 + (t1065 + t787) * t354 + (-t419 / 0.2e1 + t849 * t736 + t848 * t732 + (t138 * t732 - t139 * t736 + t392) * t1083 + (t285 + t802) * t1081 + t786) * t353 + (t181 * t189 + t182 * t190 + t354 * t362) * t1083 + (t115 * t181 + t145 * t182 + t207 * t354) * t1081 + t467 * t1049;
t9 = t457 * t812 - t1119 + t740;
t801 = t9 * qJD(1) + t10 * qJD(2);
t502 = t510 * mrSges(7,2);
t306 = mrSges(7,1) * t511 + t502;
t312 = -Ifges(7,2) * t511 + t507;
t313 = -Ifges(6,2) * t511 + t508;
t314 = Ifges(7,1) * t510 - t971;
t315 = Ifges(6,1) * t510 - t973;
t13 = t138 * t382 - t139 * t384 + t207 * t306 + t362 * t307 + t94 * t381 + (t310 / 0.2e1 + t311 / 0.2e1) * t615 + (t877 - t383) * t95 + (-t954 + t314 / 0.2e1 + t315 / 0.2e1 - t983 + t1117 * pkin(5) - t856) * t511 + (t313 / 0.2e1 + t312 / 0.2e1 - t984 - t955 - t854) * t510;
t791 = t1061 + t877 / 0.2e1;
t748 = (-t510 * t882 - t849) * t181 + (t511 * t894 + t306 / 0.2e1 + t1066) * t353 + (-t511 * t882 + t1059 + t791) * t182;
t792 = t894 + t884;
t14 = -t292 * t792 + t293 * t883 + t748;
t800 = t14 * qJD(1) + t13 * qJD(2);
t47 = m(7) * (t510 * t95 - t511 * t86) - t511 * t383 + t510 * t381;
t60 = (t457 / 0.4e1 + t926 / 0.4e1 - t924 / 0.4e1) * t1088;
t799 = -qJD(1) * t60 + qJD(2) * t47;
t797 = t670 * t732 + t674 * t736;
t789 = t1020 + t832 / 0.2e1;
t788 = t720 / 0.2e1 + t719 / 0.2e1;
t783 = t312 / 0.4e1 + t313 / 0.4e1 + t853;
t632 = t682 * t733;
t633 = t684 * t733;
t782 = -t632 / 0.4e1 - t633 / 0.4e1 - t842;
t17 = (t1005 - t932 / 0.2e1 - t812 * t733) * t500 + t1089 + t1120;
t575 = -Ifges(7,3) * t737 + t733 * t804;
t577 = -Ifges(6,3) * t737 + t733 * t805;
t686 = Ifges(5,2) * t737 + t974;
t693 = Ifges(5,1) * t737 - t974;
t739 = t308 * t1013 - t419 * t1014 + t138 * t1015 + t86 * t1016 + t189 * t1017 + t115 * t1019 + t139 * t1021 + t95 * t1022 + t190 * t1023 + t568 * t1005 + (t586 / 0.4e1 + t584 / 0.4e1) * t511 + (t580 / 0.4e1 + t582 / 0.4e1) * t510 - pkin(3) * t467 / 0.2e1 + (-t686 / 0.4e1 + t693 / 0.4e1 + t577 / 0.4e1 + t575 / 0.4e1) * t616 + t145 * t1025 + t362 * t1027 + t207 * t1028 + t392 * t1029 + t285 * t1030 + t382 * t1037 + t451 * t1038 + t381 * t1044 + t495 * t1054 + t560 * t1055 + t515 * t1056 + t563 * t1058 + t496 * t1060 + (t115 * t495 + t145 * t515 + t207 * t664 + t285 * t663 + t496 * t86 + t522 * t95) * t1081;
t742 = (-pkin(4) * t406 + (t921 - t922) * pkin(11)) * t1084 + (t133 * t670 - t146 * t674 + t284 * t715) * t1082 + pkin(4) * t1057 + t284 * t1011 + t406 * t1009 + t455 * t1079 + t456 * t1076 + t670 * t1052 + t484 * t1008 - t715 * t393 / 0.2e1;
t749 = t469 / 0.4e1 + t409 / 0.4e1 - t336 / 0.4e1 - t335 / 0.4e1 + t853 * t736 + t855 * t732 + (t1083 * t362 + t787) * pkin(10);
t755 = (-t338 / 0.4e1 - t337 / 0.4e1) * t732 + (t340 / 0.4e1 + t339 / 0.4e1) * t736 + t216 / 0.4e1 + t217 / 0.4e1 - t408 / 0.4e1 + t470 / 0.4e1;
t576 = t737 * t804 + t956;
t578 = t737 * t805 + t958;
t726 = Ifges(5,4) * t737;
t687 = -Ifges(5,2) * t733 + t726;
t692 = Ifges(5,1) * t733 + t726;
t756 = -t687 / 0.4e1 - t692 / 0.4e1 + t578 / 0.4e1 + t576 / 0.4e1 + t842 * t736 + t844 * t732;
t771 = t563 * t138 + t564 * t139 + t559 * t189 + t560 * t190;
t2 = t739 + (-t345 / 0.4e1 - t344 / 0.4e1 - t986 / 0.2e1 + t200 * t1073 - t146 * mrSges(7,3) / 0.2e1) * t736 + (-t347 / 0.4e1 - t346 / 0.4e1 + pkin(11) * t1051 + t199 * t1072 + t133 * t1071) * t732 + t742 + t834 * t558 + t836 * t557 + t771 * t1083 + (t818 - t960 / 0.2e1) * t731 + ((0.3e1 / 0.4e1 * t995 * Ifges(5,6) + t1101 * t889) * t731 + (t1083 * t392 + t786) * pkin(10) + t755) * t733 + (Ifges(5,5) * t822 + t749) * t737 + t756 * t615;
t30 = -pkin(3) * t677 + t495 * t659 + t496 * t657 + t515 * t655 + t522 * t653 + t559 * t660 + t560 * t656 + t563 * t658 + t564 * t654 + t664 * t626 + t663 * t628 + m(6) * (t559 * t563 + t560 * t564) + m(7) * (t495 * t496 + t515 * t522 + t663 * t664) + (t687 / 0.2e1 + t692 / 0.2e1 - t576 / 0.2e1 - t578 / 0.2e1 + pkin(10) * t627 + t841 * t736 + t846 * t732) * t737 + (-t686 / 0.2e1 + t693 / 0.2e1 + t575 / 0.2e1 + t577 / 0.2e1 + t843 * t736 + t845 * t732 + (m(6) * t987 + t629) * pkin(10)) * t733;
t779 = t17 * qJD(1) + t2 * qJD(2) + t30 * qJD(3);
t765 = t250 * t792 - t251 * t883;
t18 = (-t829 / 0.2e1 - t624 / 0.2e1 + t1032) * t353 + (t882 * t902 + t1017 - t789) * t182 + (-t882 * t904 + t839) * t181 + t765;
t630 = t733 * t678;
t631 = t733 * t680;
t634 = t688 * t733;
t635 = t690 * t733;
t36 = t514 * t653 + t559 * t654 - t560 * t658 + t663 * t624 + (t631 / 0.2e1 + t630 / 0.2e1) * t737 + (t832 - t657) * t515 + (pkin(10) * t625 + (t941 + t953 + t632 / 0.2e1 + t633 / 0.2e1 - t841) * t732 + (-t940 - t946 - t635 / 0.2e1 - t634 / 0.2e1 + (t626 + t993) * pkin(5) + t846) * t736) * t733;
t763 = (t993 / 0.2e1 + t1030) * pkin(5) - t634 / 0.4e1 - t635 / 0.4e1 - t844;
t741 = (-t310 / 0.4e1 - t311 / 0.4e1) * t737 + (-t630 / 0.4e1 - t631 / 0.4e1) * t615 + (-t953 / 0.2e1 - t941 / 0.2e1 + t782) * t510 + (-t946 / 0.2e1 - t940 / 0.2e1 + t763) * t511 + t138 * t1023 + t139 * t1018 + t207 * t1033 + t362 * t1031 + t381 * t1045 + t515 * t1061 + t382 * t1038 + t560 * t1059 + t306 * t1014 + t94 * t1025 + t95 * t1020;
t764 = (t1065 + t994 / 0.2e1) * pkin(5) + t314 / 0.4e1 + t315 / 0.4e1 + t855;
t746 = pkin(10) * t1066 + (t955 / 0.2e1 + t984 / 0.2e1 - t783) * t732 + (-t954 / 0.2e1 - t983 / 0.2e1 + t764) * t736;
t759 = mrSges(7,1) * t1069 + pkin(5) * t1052 + t1074 * t146 + t1075 * t200 + t1077 * t199;
t6 = t741 - t880 * t557 + (pkin(5) * t1069 - t952 / 0.2e1 + t947 / 0.2e1 - t945 / 0.2e1 + t944 / 0.2e1) * m(7) + (t1107 * t731 * t890 + t746) * t733 - t881 * t558 + t759;
t778 = -t18 * qJD(1) + t6 * qJD(2) + t36 * qJD(3);
t141 = (t732 * t653 - m(7) * (-t495 * t736 - t515 * t732) + t736 * t657) * t733;
t774 = t284 * t1082 + t943 / 0.2e1 - t942 / 0.2e1;
t39 = t774 + t1090;
t67 = (-t500 / 0.2e1 + t925 / 0.2e1 + t923 / 0.2e1) * t991;
t777 = -qJD(1) * t67 + qJD(2) * t39 - qJD(3) * t141;
t180 = t511 * t878 - t502;
t548 = -t624 - t829;
t636 = -m(7) * t990 - t1100;
t776 = qJD(2) * t180 + qJD(3) * t548 + qJD(4) * t636;
t773 = mrSges(7,3) * t1012 - t1098 / 0.4e1 - t1099 / 0.4e1 + t834;
t757 = t689 / 0.4e1 + t691 / 0.4e1 + mrSges(7,3) * t1008 + (t1010 + t992 / 0.2e1) * pkin(5) + t836;
t743 = -t791 * t674 + t837 * t615 + t757 * t511 - pkin(4) * t307 / 0.2e1 + t207 * t1007 + t362 * t1006 + t381 * t1012 + t306 * t1002;
t754 = (-t949 / 0.2e1 + t1059) * pkin(11) + (t94 / 0.2e1 - t86 / 0.2e1) * mrSges(7,3) + t783;
t767 = mrSges(7,1) * t1070 + t1074 * t145 + t1075 * t190 + t1077 * t189;
t12 = t743 + (t615 * t881 + t754) * t736 + (m(7) * t1070 - t450 / 0.2e1) * pkin(5) - t773 * t510 + t879 * t616 + (-t880 * t615 + (t951 / 0.2e1 + t1062) * pkin(11) + t764) * t732 + t767;
t20 = (-t1100 / 0.2e1 - t676 / 0.2e1 + t817 + t884 * t732) * t353;
t745 = pkin(10) * t1006 + (-t729 / 0.2e1 - t727 / 0.2e1) * pkin(11) * mrSges(6,3) + t773 * t732 + t757 * t736;
t750 = pkin(4) * t1032 + t1002 * t624 + t1007 * t663 + t1012 * t653 - t674 * t789 - t737 * t837;
t758 = pkin(11) * t1024 + t763;
t760 = pkin(11) * t1018 + (t1045 - t495 / 0.2e1) * mrSges(7,3) + t782;
t766 = mrSges(7,1) * t1050 + mrSges(6,2) * t1037 + mrSges(7,2) * t1044 + t1077 * t563;
t25 = (-t659 / 0.2e1 + m(7) * t1050) * pkin(5) + (t737 * t880 + t758) * t732 + (-t737 * t881 + t760) * t736 + (t745 + t879) * t733 + t750 + t766;
t70 = pkin(4) * t676 - t715 * t1100 - t671 * t990 + (-t1098 / 0.2e1 - t1099 / 0.2e1 - t833) * t736 + (-pkin(5) * t992 - t689 / 0.2e1 - t691 / 0.2e1 + t835) * t732;
t772 = -t20 * qJD(1) + t12 * qJD(2) + t25 * qJD(3) - t70 * qJD(4);
t769 = m(7) * (-t510 * t674 - t511 * t670 - t802);
t41 = (-t938 / 0.2e1 + t1063 - t950 / 0.2e1) * t736 + (-t939 / 0.2e1 + t1060 - t948 / 0.2e1) * t732 + t886 - t769 / 0.2e1;
t436 = -m(7) * t797 + mrSges(7,3) * t897;
t63 = (-t905 / 0.4e1 + t901 / 0.4e1 - t354 / 0.4e1) * t1088;
t762 = m(7) * ((-t670 * t733 + t515) * t736 + (t674 * t733 - t495) * t732);
t84 = (t931 / 0.2e1 + t1026) * t736 + (t933 / 0.2e1 + t1019) * t732 + t885 - t762 / 0.2e1;
t770 = qJD(1) * t63 - qJD(2) * t41 - qJD(3) * t84 + qJD(4) * t436;
t730 = t737 ^ 2;
t728 = t733 ^ 2;
t475 = t728 * t989;
t85 = t653 * t998 + t762 / 0.2e1 + t657 * t1001 + t885 + t788 * t737;
t68 = (-t923 - t925) * t991 / 0.2e1 + m(7) * t865;
t62 = m(7) * t1064 + t1081 * t798;
t61 = (t924 - t926) * t1081 + m(7) * t1053;
t42 = t769 / 0.2e1 + t383 * t1001 + t381 * t998 + t510 * t876 + t935 * t1047 + t886 - t788 * t615;
t38 = -t774 + t1090;
t24 = t958 / 0.2e1 + t956 / 0.2e1 + pkin(5) * t1016 + t750 + t745 * t733 + t760 * t736 + t758 * t732 + t496 * t894 - t766 - t1108 * t903 / 0.2e1 + t1109 * t900 / 0.2e1;
t21 = (t1006 + t1007 + t817 + (t792 + t894) * t732) * t353;
t19 = t765 + (t1081 * t892 + t1031 + t1033) * t353 + (t1081 * t899 + t1018 + t1020) * t182 + t1092 * t181 + t1110 * (t901 * t999 + t181 * t904 / 0.2e1);
t16 = t909 * t1078 + t932 * t1049 + t500 * t1005 + (t671 - t809) * t865 + t1089 - t1120;
t15 = t292 * t894 + t748 + (mrSges(6,1) + mrSges(7,1)) * t292 / 0.2e1 - t1111 * t293 / 0.2e1;
t11 = t743 + t115 * t894 + t959 / 0.2e1 + t957 / 0.2e1 + pkin(5) * t1054 + t754 * t736 + (pkin(11) * t1062 + t764) * t732 + (t1072 * t985 - t773) * t510 - t767 + t1108 * t908 / 0.2e1 - t1109 * t907 / 0.2e1;
t8 = -t1053 * t809 + t740 + (t1010 + t1079) * t457 + t1119;
t5 = t741 + t133 * t894 + t968 / 0.2e1 + t967 / 0.2e1 + t964 / 0.2e1 + t962 / 0.2e1 + t746 * t733 + (t944 - t945 + t947 - t952) * t1081 - t759 + t1107 * t795;
t3 = t738 + mrSges(5,3) * t912 / 0.2e1 + t840 * t457 + t1092 * t293 - t1113;
t1 = t739 + t986 * t998 + Ifges(5,5) * t793 + Ifges(5,6) * t796 + t146 * t876 + Ifges(5,3) * t861 + t731 * t818 - t742 - t487 * t985 / 0.2e1 + t749 * t737 + (pkin(10) * t875 + t756) * t615 + t935 * t1069 + t921 * t1072 + t922 * t1073 + (t392 * t988 + t771) * t1083 + t1104 * t736 / 0.4e1 + (t684 + t682) * t557 / 0.4e1 + (t690 + t688) * t558 / 0.4e1 + t1103 * t732 / 0.4e1 + (pkin(10) * t847 + t755 + (Ifges(5,6) + t1101) * t872 / 0.4e1) * t733;
t26 = [qJD(2) * t23 + qJD(3) * t22 + qJD(4) * t29, t3 * qJD(3) + t8 * qJD(4) + t15 * qJD(5) + t61 * qJD(6) + t919 + (-mrSges(3,2) * t814 - mrSges(3,1) * t831 + m(4) * (-pkin(2) * t815 + t621 * t556) + t556 * t650 + (-mrSges(4,1) * t995 + mrSges(4,2) * t734) * t815 + (-m(5) * t391 + m(6) * t362 + t1102 + t1117) * t457 + (m(5) * t392 + t525) * t458 + (-m(4) * t619 + t1114) * t555 + (m(6) * t139 + t1118 + t381 + t382) * t293 + (m(6) * t138 + m(7) * t86 + t383 + t384) * t292) * qJD(2), t920 + t3 * qJD(2) + t16 * qJD(4) + t19 * qJD(5) + t68 * qJD(6) + ((t250 * t495 + t251 * t515 - t663 * t909) * t1081 + (t250 * t559 + t251 * t560 - t475) * t1083 + (-pkin(3) * t501 - t730 * t989 - t475) * t1085) * t1086 + ((mrSges(4,2) + (-t626 - t627) * t733 + (-t728 - t730) * mrSges(5,3)) * t500 + t1106 * t501 + (t653 + t654) * t251 + (t657 + t658) * t250) * qJD(3), t8 * qJD(2) + t16 * qJD(3) + t21 * qJD(5) + t62 * qJD(6) + t916 + ((t671 + t930) * t354 + (-t1110 * t897 + mrSges(5,2)) * t353 + 0.2e1 * (t353 * t797 + t354 * t715) * t1081 + 0.2e1 * (-pkin(11) * t353 * t897 - pkin(4) * t354) * t1083) * qJD(4), t15 * qJD(2) + t19 * qJD(3) + t21 * qJD(4) + (-t1111 * t181 + (-mrSges(6,1) + t878) * t182) * qJD(5), qJD(2) * t61 + qJD(3) * t68 + qJD(4) * t62; qJD(3) * t4 + qJD(4) * t9 + qJD(5) * t14 - qJD(6) * t60 - t919, qJD(3) * t7 + qJD(4) * t10 + qJD(5) * t13 + qJD(6) * t47, t1 * qJD(4) + t5 * qJD(5) + t38 * qJD(6) + ((-pkin(3) * t621 + (-t455 * t733 + t456 * t737) * pkin(10)) * t1085 + (t199 * t559 + t200 * t560 + t406 * t988) * t1083 + (t133 * t495 + t146 * t515 + t284 * t663) * t1081) * t1086 + t803 + (-t619 * mrSges(4,2) - pkin(3) * t561 + t133 * t657 + t146 * t653 + t199 * t658 + t200 * t654 + t284 * t626 + t663 * t393 + t406 * t627 + t515 * t484 + t560 * t485 + t495 * t486 + t559 * t487 + t707 + t1106 * t621 + t841 * t558 + (t579 + t581) * t1042 + (t512 / 0.2e1 - t342 / 0.2e1 - t343 / 0.2e1 + t456 * mrSges(5,3) + pkin(10) * t592) * t737 + (t1046 - t455 * mrSges(5,3) + (t346 / 0.2e1 + t347 / 0.2e1) * t736 + (-t344 / 0.2e1 - t345 / 0.2e1) * t732 + (-t593 + t394) * pkin(10)) * t733 + (t692 * t819 + (Ifges(5,6) * t737 / 0.2e1 - Ifges(4,6)) * t734 + (t686 * t890 + t969 / 0.2e1 + (t575 + t577) * t995 / 0.2e1) * t733) * t731) * qJD(3), t1 * qJD(3) + t11 * qJD(5) + t42 * qJD(6) + (t190 * mrSges(6,3) + t145 * mrSges(7,3) - t833 * t615 + (m(6) * t190 + t449) * pkin(11) + t851) * t895 + (-t189 * mrSges(6,3) - t115 * mrSges(7,3) + t835 * t615 + (-m(6) * t189 - t451) * pkin(11) + t850) * t896 + t801 + (-t715 * t419 - t674 * t448 + t670 * t450 + t285 * t671 + pkin(4) * t420 - t391 * mrSges(5,2) + m(7) * (t115 * t670 - t145 * t674 + t285 * t715) + t898 + t828 * t392 + t1093 * t616) * qJD(4), t5 * qJD(3) + t11 * qJD(4) + (-mrSges(6,1) * t139 - mrSges(7,1) * t95 - mrSges(6,2) * t138 - mrSges(7,2) * t94 + (-t950 - t1118) * pkin(5) + t311 + t310) * qJD(5) + t800, qJD(3) * t38 + qJD(4) * t42 + t799; -qJD(2) * t4 + qJD(4) * t17 - qJD(5) * t18 - qJD(6) * t67 - t920, qJD(4) * t2 + qJD(5) * t6 + qJD(6) * t39 - t803, qJD(4) * t30 + qJD(5) * t36 - qJD(6) * t141, t24 * qJD(5) + t85 * qJD(6) + (t522 * mrSges(7,3) + t564 * mrSges(6,3) + t833 * t737 + (m(6) * t564 + t656) * pkin(11) - t845) * t895 + (-t496 * mrSges(7,3) - t563 * mrSges(6,3) - t835 * t737 + (-m(6) * t563 - t660) * pkin(11) + t843) * t896 + t779 + (t723 + t715 * t628 - t674 * t655 + t670 * t659 + t664 * t671 - pkin(4) * t629 + m(7) * (t496 * t670 - t522 * t674 + t664 * t715) + (-Ifges(5,6) + t1093) * t733 + (t737 * t828 + t934) * pkin(10)) * qJD(4), t24 * qJD(4) + t778 + (-mrSges(6,1) * t560 - mrSges(6,2) * t559 - mrSges(7,2) * t514 + (-t1108 * t736 + (mrSges(7,3) * pkin(5) - t1109) * t732) * t733 + t878 * t515) * qJD(5), qJD(4) * t85 + t777; -qJD(2) * t9 - qJD(3) * t17 - qJD(5) * t20 + qJD(6) * t63 - t916, -qJD(3) * t2 + qJD(5) * t12 - qJD(6) * t41 - t801, qJD(5) * t25 - qJD(6) * t84 - t779, -qJD(5) * t70 + qJD(6) * t436, t772 + (-mrSges(7,2) * t670 - pkin(5) * t979 - pkin(11) * t980 + t721 + t722 + (mrSges(6,2) * pkin(11) - t1108) * t732 - t878 * t674) * qJD(5), t770; -qJD(2) * t14 + qJD(3) * t18 + qJD(4) * t20, -qJD(3) * t6 - qJD(4) * t12 + qJD(6) * t180 - t800, -qJD(4) * t25 + qJD(6) * t548 - t778, qJD(6) * t636 - t772, 0, t776; qJD(2) * t60 + qJD(3) * t67 - qJD(4) * t63, -qJD(3) * t39 + qJD(4) * t41 - qJD(5) * t180 - t799, qJD(4) * t84 - qJD(5) * t548 - t777, -qJD(5) * t636 - t770, -t776, 0;];
Cq  = t26;
