% Calculate time derivative of joint inertia matrix for
% S6RRRRRP3
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
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP3_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP3_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP3_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP3_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:05:54
% EndTime: 2019-03-10 01:07:19
% DurationCPUTime: 60.12s
% Computational Cost: add. (109304->1284), mult. (108425->1663), div. (0->0), fcn. (103483->10), ass. (0->638)
t1028 = Icges(6,4) + Icges(7,4);
t1027 = Icges(6,1) + Icges(7,1);
t1023 = Icges(7,5) + Icges(6,5);
t1026 = Icges(6,2) + Icges(7,2);
t1022 = -Icges(6,6) - Icges(7,6);
t622 = qJ(4) + qJ(5);
t611 = cos(t622);
t1025 = t1028 * t611;
t609 = sin(t622);
t1024 = t1028 * t609;
t1021 = Icges(7,3) + Icges(6,3);
t629 = cos(qJ(1));
t623 = qJ(2) + qJ(3);
t612 = cos(t623);
t626 = sin(qJ(1));
t878 = t612 * t626;
t526 = -t609 * t878 - t611 * t629;
t527 = -t609 * t629 + t611 * t878;
t610 = sin(t623);
t883 = t610 * t626;
t392 = Icges(7,5) * t527 + Icges(7,6) * t526 + Icges(7,3) * t883;
t394 = Icges(6,5) * t527 + Icges(6,6) * t526 + Icges(6,3) * t883;
t1009 = t392 + t394;
t877 = t612 * t629;
t528 = -t609 * t877 + t611 * t626;
t529 = t609 * t626 + t611 * t877;
t882 = t610 * t629;
t393 = Icges(7,5) * t529 + Icges(7,6) * t528 + Icges(7,3) * t882;
t395 = Icges(6,5) * t529 + Icges(6,6) * t528 + Icges(6,3) * t882;
t1008 = t393 + t395;
t396 = Icges(7,4) * t527 + Icges(7,2) * t526 + Icges(7,6) * t883;
t398 = Icges(6,4) * t527 + Icges(6,2) * t526 + Icges(6,6) * t883;
t1020 = t396 + t398;
t397 = Icges(7,4) * t529 + Icges(7,2) * t528 + Icges(7,6) * t882;
t399 = Icges(6,4) * t529 + Icges(6,2) * t528 + Icges(6,6) * t882;
t1019 = t397 + t399;
t400 = Icges(7,1) * t527 + Icges(7,4) * t526 + Icges(7,5) * t883;
t402 = Icges(6,1) * t527 + Icges(6,4) * t526 + Icges(6,5) * t883;
t1018 = t400 + t402;
t401 = Icges(7,1) * t529 + Icges(7,4) * t528 + Icges(7,5) * t882;
t403 = Icges(6,1) * t529 + Icges(6,4) * t528 + Icges(6,5) * t882;
t1017 = t401 + t403;
t1016 = t1022 * t609 + t1023 * t611;
t1015 = t1026 * t609 - t1025;
t1014 = t1027 * t611 - t1024;
t1013 = t1009 * t882 + t1018 * t529 + t1020 * t528;
t1012 = t1008 * t882 + t1017 * t529 + t1019 * t528;
t618 = qJD(4) + qJD(5);
t619 = qJD(2) + qJD(3);
t879 = t612 * t619;
t1011 = t1016 * t879 + (t1021 * t619 + (t1022 * t611 - t1023 * t609) * t618) * t610;
t998 = t1015 * t879 + (t1022 * t619 + (t1026 * t611 + t1024) * t618) * t610;
t1010 = t1014 * t879 + (t1023 * t619 + (-t1027 * t609 - t1025) * t618) * t610;
t1006 = t1016 * t610 - t1021 * t612;
t997 = -t1015 * t610 + t1022 * t612;
t996 = -t1014 * t610 + t1023 * t612;
t1007 = t609 * t997;
t1005 = t996 * t611;
t767 = -t612 * t618 + qJD(1);
t682 = t626 * t767;
t832 = qJD(1) * t612;
t766 = -t618 + t832;
t876 = t619 * t626;
t815 = t610 * t876;
t946 = t629 * t766 - t815;
t359 = -t609 * t946 + t611 * t682;
t360 = t609 * t682 + t611 * t946;
t812 = t612 * t876;
t830 = qJD(1) * t629;
t653 = t610 * t830 + t812;
t237 = Icges(7,5) * t360 + Icges(7,6) * t359 + Icges(7,3) * t653;
t239 = Icges(6,5) * t360 + Icges(6,6) * t359 + Icges(6,3) * t653;
t241 = Icges(7,4) * t360 + Icges(7,2) * t359 + Icges(7,6) * t653;
t243 = Icges(6,4) * t360 + Icges(6,2) * t359 + Icges(6,6) * t653;
t245 = Icges(7,1) * t360 + Icges(7,4) * t359 + Icges(7,5) * t653;
t247 = Icges(6,1) * t360 + Icges(6,4) * t359 + Icges(6,5) * t653;
t681 = t629 * t767;
t874 = t619 * t629;
t814 = t610 * t874;
t947 = t626 * t766 + t814;
t357 = t609 * t947 + t611 * t681;
t358 = t609 * t681 - t611 * t947;
t831 = qJD(1) * t626;
t792 = t610 * t831;
t811 = t612 * t874;
t652 = -t792 + t811;
t1004 = (t237 + t239) * t882 + t1009 * t652 + (t245 + t247) * t529 + (t241 + t243) * t528 + t1018 * t358 + t1020 * t357;
t236 = Icges(7,5) * t358 + Icges(7,6) * t357 + Icges(7,3) * t652;
t238 = Icges(6,5) * t358 + Icges(6,6) * t357 + Icges(6,3) * t652;
t240 = Icges(7,4) * t358 + Icges(7,2) * t357 + Icges(7,6) * t652;
t242 = Icges(6,4) * t358 + Icges(6,2) * t357 + Icges(6,6) * t652;
t244 = Icges(7,1) * t358 + Icges(7,4) * t357 + Icges(7,5) * t652;
t246 = Icges(6,1) * t358 + Icges(6,4) * t357 + Icges(6,5) * t652;
t1003 = (t236 + t238) * t882 + t1008 * t652 + (t244 + t246) * t529 + (t240 + t242) * t528 + t1017 * t358 + t1019 * t357;
t978 = t1006 * t652 + t1010 * t529 + t1011 * t882 + t357 * t997 - t358 * t996 - t528 * t998;
t977 = t1006 * t653 + t1010 * t527 + t1011 * t883 + t359 * t997 - t360 * t996 - t526 * t998;
t694 = -t398 * t609 + t402 * t611;
t696 = -t396 * t609 + t400 * t611;
t1002 = -t1009 * t612 + (t694 + t696) * t610;
t693 = -t399 * t609 + t403 * t611;
t695 = -t397 * t609 + t401 * t611;
t1001 = -t1008 * t612 + (t693 + t695) * t610;
t1000 = t1006 * t883 + t526 * t997 - t527 * t996;
t999 = t1006 * t882 + t528 * t997 - t529 * t996;
t733 = rSges(7,1) * t611 - rSges(7,2) * t609;
t630 = -pkin(10) - pkin(9);
t617 = -qJ(6) + t630;
t836 = t617 - t630;
t627 = cos(qJ(4));
t605 = t627 * pkin(4) + pkin(3);
t577 = pkin(5) * t611 + t605;
t841 = t577 - t605;
t970 = (-rSges(7,3) + t836) * t612 + (t733 + t841) * t610;
t995 = t1012 * t629 + t1013 * t626;
t989 = -t1012 * t626 + t1013 * t629;
t186 = t394 * t883 + t398 * t526 + t402 * t527;
t187 = t395 * t883 + t399 * t526 + t403 * t527;
t710 = t186 * t626 + t187 * t629;
t184 = t392 * t883 + t396 * t526 + t400 * t527;
t185 = t393 * t883 + t397 * t526 + t401 * t527;
t712 = t184 * t626 + t185 * t629;
t994 = t710 + t712;
t988 = (t184 + t186) * t629 + (-t185 - t187) * t626;
t624 = sin(qJ(4));
t872 = t624 * t626;
t603 = pkin(4) * t872;
t842 = -t605 * t877 - t603;
t931 = pkin(4) * t624;
t580 = pkin(5) * t609 + t931;
t991 = t529 * rSges(7,1) + t528 * rSges(7,2) + rSges(7,3) * t882 + t577 * t877 + t626 * t580;
t972 = -t836 * t882 + t842 + t991;
t659 = -t610 * t617 + t612 * t841;
t734 = -t527 * rSges(7,1) - t526 * rSges(7,2);
t871 = t624 * t629;
t604 = pkin(4) * t871;
t881 = t610 * t630;
t839 = t626 * t881 + t604;
t973 = rSges(7,3) * t883 - t629 * t580 + t626 * t659 - t734 + t839;
t993 = -t626 * t973 - t629 * t972;
t992 = t1006 * t612 + (t1005 + t1007) * t610;
t880 = t611 * t618;
t884 = t610 * t619;
t990 = -t1006 * t884 + t879 * t1005 + t1011 * t612 + (-t1010 * t611 + t997 * t880) * t610;
t826 = qJD(4) * t624;
t821 = pkin(4) * t826;
t886 = t609 * t618;
t555 = -pkin(5) * t886 - t821;
t825 = qJD(4) * t627;
t820 = pkin(4) * t825;
t556 = pkin(5) * t880 + t820;
t987 = t358 * rSges(7,1) + t357 * rSges(7,2) + rSges(7,3) * t811 + qJD(6) * t882 + t555 * t877 + t626 * t556 + t580 * t830 + t617 * t792;
t785 = t879 / 0.2e1;
t967 = -t629 * t785 + t792 / 0.2e1;
t780 = t830 / 0.2e1;
t966 = -t610 * t780 - t626 * t785;
t986 = (t619 * t995 - t978) * t612 + (t989 * qJD(1) + t1003 * t629 + t1004 * t626 + t999 * t619) * t610;
t52 = t237 * t883 + t241 * t526 + t245 * t527 + t359 * t396 + t360 * t400 + t392 * t653;
t53 = t236 * t883 + t240 * t526 + t244 * t527 + t359 * t397 + t360 * t401 + t393 * t653;
t54 = t239 * t883 + t243 * t526 + t247 * t527 + t359 * t398 + t360 * t402 + t394 * t653;
t55 = t238 * t883 + t242 * t526 + t246 * t527 + t359 * t399 + t360 * t403 + t395 * t653;
t985 = (t619 * t994 - t977) * t612 + ((t53 + t55) * t629 + (t52 + t54) * t626 + t1000 * t619 + t988 * qJD(1)) * t610;
t984 = t995 * qJD(1) + t1003 * t626 - t1004 * t629;
t27 = qJD(1) * t712 - t52 * t629 + t53 * t626;
t28 = qJD(1) * t710 - t54 * t629 + t55 * t626;
t983 = t27 + t28;
t66 = (t619 * t696 - t237) * t612 + (t392 * t619 + (-t396 * t618 + t245) * t611 + (-t400 * t618 - t241) * t609) * t610;
t68 = (t619 * t694 - t239) * t612 + (t394 * t619 + (-t398 * t618 + t247) * t611 + (-t402 * t618 - t243) * t609) * t610;
t982 = -t66 - t68;
t67 = (t619 * t695 - t236) * t612 + (t393 * t619 + (-t397 * t618 + t244) * t611 + (-t401 * t618 - t240) * t609) * t610;
t69 = (t619 * t693 - t238) * t612 + (t395 * t619 + (-t399 * t618 + t246) * t611 + (-t403 * t618 - t242) * t609) * t610;
t981 = -t67 - t69;
t980 = -t1000 * t612 + t610 * t994;
t979 = t610 * t995 - t612 * t999;
t867 = t990 + (t997 * t879 + (-t618 * t996 - t998) * t610) * t609;
t649 = -t619 * t836 + t821;
t965 = qJD(1) * t881 + t820;
t793 = qJD(1) * t604 + t626 * t965;
t976 = -t841 * t814 + (t629 * t649 - t831 * t841) * t612 - t793 - rSges(7,3) * t792 + t987;
t810 = t612 * t872;
t731 = qJD(4) * pkin(4) * t810 + t605 * t815 + t629 * t965 + t630 * t812;
t735 = rSges(7,1) * t360 + rSges(7,2) * t359;
t954 = (t577 * t619 - qJD(6)) * t610;
t975 = (qJD(1) * t659 - t556) * t629 + ((-t617 * t619 + t555) * t612 - t954 + (t580 - t931) * qJD(1)) * t626 + t731 + rSges(7,3) * t653 + t735;
t974 = (t619 * t841 - qJD(6)) * t612 + t733 * t879 + (t555 + t649 + rSges(7,3) * t619 + (-rSges(7,1) * t609 - rSges(7,2) * t611) * t618) * t610;
t971 = t970 * t831;
t969 = t1001 * t629 + t1002 * t626;
t968 = -t1001 * t626 + t1002 * t629;
t631 = -pkin(8) - pkin(7);
t625 = sin(qJ(2));
t828 = qJD(2) * t625;
t823 = pkin(2) * t828;
t964 = qJD(1) * t631 + t823;
t963 = t1000 + t1002;
t962 = t999 + t1001;
t628 = cos(qJ(2));
t586 = rSges(3,1) * t625 + rSges(3,2) * t628;
t663 = qJD(2) * t586;
t960 = t626 * t663;
t910 = Icges(3,4) * t628;
t724 = -Icges(3,2) * t625 + t910;
t539 = Icges(3,6) * t626 + t629 * t724;
t911 = Icges(3,4) * t625;
t730 = Icges(3,1) * t628 - t911;
t541 = Icges(3,5) * t626 + t629 * t730;
t684 = t539 * t625 - t541 * t628;
t959 = t626 * t684;
t908 = Icges(4,4) * t612;
t722 = -Icges(4,2) * t610 + t908;
t514 = Icges(4,6) * t626 + t629 * t722;
t909 = Icges(4,4) * t610;
t728 = Icges(4,1) * t612 - t909;
t516 = Icges(4,5) * t626 + t629 * t728;
t686 = t514 * t610 - t516 * t612;
t958 = t626 * t686;
t606 = pkin(2) * t628 + pkin(1);
t928 = pkin(1) - t606;
t957 = t626 * t928;
t538 = -Icges(3,6) * t629 + t626 * t724;
t540 = -Icges(3,5) * t629 + t626 * t730;
t685 = t538 * t625 - t540 * t628;
t956 = t629 * t685;
t513 = -Icges(4,6) * t629 + t626 * t722;
t515 = -Icges(4,5) * t629 + t626 * t728;
t687 = t513 * t610 - t515 * t612;
t955 = t629 * t687;
t953 = t992 * t884;
t737 = -t527 * rSges(6,1) - t526 * rSges(6,2);
t407 = rSges(6,3) * t883 - t737;
t409 = t529 * rSges(6,1) + t528 * rSges(6,2) + rSges(6,3) * t882;
t952 = -t626 * t407 - t629 * t409;
t716 = Icges(5,5) * t627 - Icges(5,6) * t624;
t385 = t716 * t879 + (Icges(5,3) * t619 + (-Icges(5,5) * t624 - Icges(5,6) * t627) * qJD(4)) * t610;
t906 = Icges(5,4) * t627;
t721 = -Icges(5,2) * t624 + t906;
t503 = -Icges(5,6) * t612 + t610 * t721;
t895 = t503 * t624;
t951 = -t619 * t895 - t385;
t869 = t627 * t629;
t550 = -t810 - t869;
t870 = t626 * t627;
t551 = t612 * t870 - t871;
t740 = -t551 * rSges(5,1) - t550 * rSges(5,2);
t448 = rSges(5,3) * t883 - t740;
t552 = -t612 * t871 + t870;
t553 = t612 * t869 + t872;
t449 = t553 * rSges(5,1) + t552 * rSges(5,2) + rSges(5,3) * t882;
t950 = -t626 * t448 - t629 * t449;
t717 = Icges(4,5) * t612 - Icges(4,6) * t610;
t511 = -Icges(4,3) * t629 + t626 * t717;
t949 = qJD(1) * t511;
t718 = Icges(3,5) * t628 - Icges(3,6) * t625;
t536 = -Icges(3,3) * t629 + t626 * t718;
t763 = -qJD(4) + t832;
t948 = t626 * t763 + t814;
t561 = Icges(4,2) * t612 + t909;
t562 = Icges(4,1) * t610 + t908;
t683 = t561 * t610 - t562 * t612;
t945 = qJD(1) * t683 + t717 * t619;
t944 = 2 * m(3);
t943 = 2 * m(4);
t942 = 2 * m(5);
t941 = 2 * m(6);
t940 = 2 * m(7);
t620 = t626 ^ 2;
t621 = t629 ^ 2;
t939 = -t612 / 0.2e1;
t938 = t626 / 0.2e1;
t937 = -t629 / 0.2e1;
t936 = -rSges(5,3) - pkin(9);
t935 = m(3) * t586;
t564 = rSges(4,1) * t610 + rSges(4,2) * t612;
t934 = m(4) * t564;
t933 = pkin(2) * t625;
t932 = pkin(3) * t612;
t930 = pkin(9) * t610;
t929 = t626 * pkin(7);
t616 = t629 * pkin(7);
t927 = pkin(3) - t605;
t926 = -pkin(7) - t631;
t925 = pkin(9) + t630;
t924 = t953 + (-t619 * t969 - t867) * t612 + (-qJD(1) * t968 + t626 * t982 + t629 * t981) * t610;
t922 = rSges(3,1) * t628;
t921 = rSges(4,1) * t612;
t920 = rSges(3,2) * t625;
t919 = rSges(3,3) * t629;
t614 = t626 * rSges(3,3);
t613 = t626 * rSges(4,3);
t918 = t66 * t629;
t917 = t67 * t626;
t916 = t68 * t629;
t915 = t69 * t626;
t764 = -qJD(4) * t612 + qJD(1);
t679 = t627 * t764;
t419 = t626 * t679 + (-t629 * t763 + t815) * t624;
t678 = t764 * t624;
t875 = t619 * t627;
t420 = t763 * t869 + (-t610 * t875 + t678) * t626;
t280 = Icges(5,5) * t420 + Icges(5,6) * t419 + Icges(5,3) * t653;
t282 = Icges(5,4) * t420 + Icges(5,2) * t419 + Icges(5,6) * t653;
t284 = Icges(5,1) * t420 + Icges(5,4) * t419 + Icges(5,5) * t653;
t442 = Icges(5,5) * t551 + Icges(5,6) * t550 + Icges(5,3) * t883;
t444 = Icges(5,4) * t551 + Icges(5,2) * t550 + Icges(5,6) * t883;
t446 = Icges(5,1) * t551 + Icges(5,4) * t550 + Icges(5,5) * t883;
t692 = -t444 * t624 + t446 * t627;
t78 = (t619 * t692 - t280) * t612 + (-t282 * t624 + t284 * t627 + t442 * t619 + (-t444 * t627 - t446 * t624) * qJD(4)) * t610;
t914 = t78 * t629;
t417 = t624 * t948 + t629 * t679;
t418 = -t627 * t948 + t629 * t678;
t279 = Icges(5,5) * t418 + Icges(5,6) * t417 + Icges(5,3) * t652;
t281 = Icges(5,4) * t418 + Icges(5,2) * t417 + Icges(5,6) * t652;
t283 = Icges(5,1) * t418 + Icges(5,4) * t417 + Icges(5,5) * t652;
t443 = Icges(5,5) * t553 + Icges(5,6) * t552 + Icges(5,3) * t882;
t445 = Icges(5,4) * t553 + Icges(5,2) * t552 + Icges(5,6) * t882;
t447 = Icges(5,1) * t553 + Icges(5,4) * t552 + Icges(5,5) * t882;
t691 = -t445 * t624 + t447 * t627;
t79 = (t619 * t691 - t279) * t612 + (-t281 * t624 + t283 * t627 + t443 * t619 + (-t445 * t627 - t447 * t624) * qJD(4)) * t610;
t913 = t79 * t626;
t912 = -rSges(7,3) + t617;
t907 = Icges(5,4) * t624;
t742 = -rSges(4,2) * t610 + t921;
t533 = t742 * t619;
t894 = t533 * t626;
t893 = t538 * t628;
t892 = t539 * t628;
t891 = t540 * t625;
t890 = t541 * t625;
t889 = t561 * t619;
t888 = t562 * t619;
t386 = t721 * t879 + (Icges(5,6) * t619 + (-Icges(5,2) * t627 - t907) * qJD(4)) * t610;
t873 = t624 * t386;
t868 = t629 * t631;
t738 = rSges(6,1) * t360 + rSges(6,2) * t359;
t251 = rSges(6,3) * t653 + t738;
t865 = t251 * t882 + t407 * t811;
t804 = t358 * rSges(6,1) + t357 * rSges(6,2) + rSges(6,3) * t811;
t249 = -rSges(6,3) * t792 + t804;
t576 = pkin(9) * t811;
t670 = -t619 * t630 - t821;
t289 = -t576 + (pkin(9) * t831 + t874 * t927) * t610 + (t629 * t670 + t831 * t927) * t612 + t793;
t864 = -t249 - t289;
t575 = pkin(3) * t815;
t782 = t927 * t612;
t668 = -t782 - t930;
t290 = -pkin(9) * t812 + t575 + (t629 * t668 + t603) * qJD(1) - t731;
t438 = t626 * t668 - t839;
t863 = t290 * t882 + t438 * t811;
t861 = t973 * t882;
t736 = rSges(6,1) * t611 - rSges(6,2) * t609;
t356 = t736 * t879 + (rSges(6,3) * t619 + (-rSges(6,1) * t609 - rSges(6,2) * t611) * t618) * t610;
t788 = t610 * t826;
t391 = -pkin(4) * t788 + (-t610 * t925 - t782) * t619;
t858 = -t356 - t391;
t493 = -rSges(6,3) * t612 + t610 * t736;
t483 = t493 * t831;
t857 = t409 * t884 + t610 * t483;
t594 = pkin(3) * t877;
t548 = pkin(9) * t882 + t594;
t674 = -t629 * t881 - t842;
t439 = t674 - t548;
t484 = -t610 * t927 + t612 * t925;
t480 = t484 * t831;
t856 = t439 * t884 + t610 * t480;
t739 = rSges(5,1) * t627 - rSges(5,2) * t624;
t390 = t739 * t879 + (rSges(5,3) * t619 + (-rSges(5,1) * t624 - rSges(5,2) * t627) * qJD(4)) * t610;
t744 = t930 + t932;
t534 = t744 * t619;
t855 = -t390 - t534;
t854 = -t407 - t438;
t853 = -t409 - t439;
t852 = t612 * t438 + t484 * t883;
t851 = -t449 - t548;
t311 = t612 * t407 + t493 * t883;
t565 = pkin(3) * t610 - pkin(9) * t612;
t549 = t565 * t831;
t849 = t480 + t549;
t848 = -t484 - t493;
t505 = -rSges(5,3) * t612 + t610 * t739;
t485 = t505 * t831;
t847 = t485 + t549;
t508 = t616 + t868 - t957;
t588 = t629 * t606;
t509 = -t629 * pkin(1) + t626 * t926 + t588;
t846 = t626 * t508 + t629 * t509;
t517 = -t629 * rSges(4,3) + t626 * t742;
t518 = rSges(4,1) * t877 - rSges(4,2) * t882 + t613;
t421 = t626 * t517 + t629 * t518;
t845 = -t505 - t565;
t547 = t744 * t626;
t844 = t626 * t547 + t629 * t548;
t840 = rSges(4,2) * t792 + rSges(4,3) * t830;
t838 = t964 * t626;
t837 = t629 * t922 + t614;
t835 = t620 + t621;
t512 = Icges(4,3) * t626 + t629 * t717;
t834 = qJD(1) * t512;
t537 = Icges(3,3) * t626 + t629 * t718;
t833 = qJD(1) * t537;
t827 = qJD(2) * t628;
t824 = t629 * t920;
t822 = pkin(2) * t827;
t809 = -t289 - t976;
t808 = -t391 - t974;
t807 = -t438 - t973;
t806 = -t439 - t972;
t803 = -t534 + t858;
t727 = Icges(5,1) * t627 - t907;
t387 = t727 * t879 + (Icges(5,5) * t619 + (-Icges(5,1) * t624 - t906) * qJD(4)) * t610;
t502 = -Icges(5,3) * t612 + t610 * t716;
t504 = -Icges(5,5) * t612 + t610 * t727;
t802 = t610 * t627 * t387 + t612 * t504 * t875 + t502 * t884;
t654 = -t612 * t831 - t814;
t665 = t564 * t619;
t801 = t626 * (-t626 * t665 + (t629 * t742 + t613) * qJD(1)) + t629 * (rSges(4,1) * t654 - rSges(4,2) * t811 + t840) + t517 * t830;
t800 = -t548 + t853;
t799 = t418 * rSges(5,1) + t417 * rSges(5,2) + rSges(5,3) * t811;
t798 = t626 * (pkin(9) * t653 + qJD(1) * t594 - t575) + t629 * (pkin(3) * t654 - pkin(9) * t792 + t576) + t547 * t830;
t797 = -t484 - t970;
t796 = t626 * ((-t629 * t928 - t929) * qJD(1) - t838) + t629 * (-t629 * t823 + (t629 * t926 + t957) * qJD(1)) + t508 * t830;
t795 = t483 + t849;
t794 = -t565 + t848;
t789 = t625 * t831;
t787 = t883 / 0.2e1;
t786 = t882 / 0.2e1;
t230 = -t442 * t612 + t610 * t692;
t298 = t502 * t883 + t503 * t550 + t504 * t551;
t784 = t230 / 0.2e1 + t298 / 0.2e1;
t231 = -t443 * t612 + t610 * t691;
t299 = t502 * t882 + t503 * t552 + t504 * t553;
t783 = t231 / 0.2e1 + t299 / 0.2e1;
t781 = t831 / 0.2e1;
t779 = -t564 - t933;
t778 = -t565 - t933;
t777 = t629 * t970;
t776 = t629 * t848;
t775 = t972 * t612;
t774 = t853 * t612;
t452 = t845 * t629;
t424 = -qJD(1) * t513 - t629 * t889;
t773 = t516 * t619 + t424;
t425 = qJD(1) * t514 - t626 * t889;
t772 = t515 * t619 + t425;
t426 = -qJD(1) * t515 - t629 * t888;
t771 = -t514 * t619 + t426;
t427 = qJD(1) * t516 - t626 * t888;
t770 = t513 * t619 - t427;
t769 = -t577 * t612 - t606;
t768 = -t626 * t631 + t588;
t762 = t811 * t973 + t882 * t975;
t761 = t612 * t251 + t356 * t883 + t493 * t653;
t760 = t612 * t290 + t391 * t883 + t484 * t653;
t759 = -t534 + t808;
t758 = t610 * t971 + t884 * t972;
t757 = -t548 + t806;
t177 = t612 * t973 + t883 * t970;
t756 = t626 * t438 + t629 * t439 + t844;
t755 = t849 + t971;
t754 = -t565 + t797;
t266 = t844 - t950;
t749 = -t484 + t778;
t748 = -t505 + t778;
t747 = -t534 - t822;
t746 = t629 * t797;
t745 = t806 * t612;
t326 = t794 * t629;
t743 = -t920 + t922;
t741 = rSges(5,1) * t420 + rSges(5,2) * t419;
t729 = Icges(3,1) * t625 + t910;
t723 = Icges(3,2) * t628 + t911;
t560 = Icges(4,5) * t610 + Icges(4,6) * t612;
t217 = t442 * t883 + t444 * t550 + t446 * t551;
t218 = t443 * t883 + t445 * t550 + t447 * t551;
t701 = t217 * t629 - t218 * t626;
t700 = t217 * t626 + t218 * t629;
t219 = t442 * t882 + t444 * t552 + t446 * t553;
t220 = t443 * t882 + t445 * t552 + t447 * t553;
t699 = t219 * t629 - t220 * t626;
t698 = t219 * t626 + t220 * t629;
t697 = t230 * t626 + t231 * t629;
t690 = t448 * t629 - t449 * t626;
t680 = -t493 + t749;
t677 = -t390 + t747;
t676 = -t391 + t747;
t288 = t754 * t629;
t675 = -pkin(1) - t743;
t432 = t748 * t629;
t673 = -t606 - t742;
t285 = -rSges(5,3) * t792 + t799;
t286 = rSges(5,3) * t653 + t741;
t672 = t629 * t285 + t626 * t286 + t448 * t830 + t798;
t671 = t629 * t289 + t626 * t290 + t438 * t830 + t798;
t170 = t756 - t952;
t669 = -rSges(6,3) * t610 - t605 * t612 - t606;
t666 = t749 - t970;
t664 = -t356 + t676;
t660 = t619 * t560;
t658 = qJD(2) * t729;
t657 = qJD(2) * t723;
t656 = qJD(2) * (-Icges(3,5) * t625 - Icges(3,6) * t628);
t315 = t680 * t629;
t655 = t610 * t936 - t606 - t932;
t651 = t676 - t974;
t650 = t669 * t626;
t648 = t610 * t912 + t769;
t269 = t666 * t629;
t647 = t975 * t612 + t653 * t970 + t974 * t883;
t139 = t756 - t993;
t646 = t629 * t249 + t626 * t251 + t407 * t830 + t671;
t645 = (t969 * t610 + t612 * t992) * t884 + t985 * t883 + t986 * t882 + t979 * t811 + t653 * t980;
t644 = rSges(3,2) * t789 + rSges(3,3) * t830 - t629 * t663;
t316 = -t511 * t629 - t626 * t687;
t317 = -t512 * t629 - t958;
t318 = t511 * t626 - t955;
t319 = t512 * t626 - t629 * t686;
t70 = t280 * t882 + t282 * t552 + t284 * t553 + t417 * t444 + t418 * t446 + t442 * t652;
t71 = t279 * t882 + t281 * t552 + t283 * t553 + t417 * t445 + t418 * t447 + t443 * t652;
t35 = qJD(1) * t698 + t626 * t71 - t629 * t70;
t422 = -t629 * t660 - t949;
t423 = -t626 * t660 + t834;
t641 = (-t316 * t629 - t701 - t988) * t831 + (-t318 * t629 - t699 - t989) * t830 + (t35 + t317 * t831 + t319 * t830 + (t319 * qJD(1) + (t425 * t610 - t427 * t612 + t513 * t879 + t515 * t884 - t949) * t629) * t629 + ((t318 + t958) * qJD(1) + (-t423 + t771 * t612 - t773 * t610 + (t512 - t687) * qJD(1)) * t629 + t626 * t422) * t626 + t984) * t626;
t640 = t626 * t655 - t868;
t639 = t626 * t975 + t629 * t976 + t830 * t973 + t671;
t531 = t722 * t619;
t532 = t728 * t619;
t638 = qJD(1) * t560 + (t532 - t889) * t612 + (-t531 - t888) * t610;
t72 = t280 * t883 + t282 * t550 + t284 * t551 + t419 * t444 + t420 * t446 + t442 * t653;
t73 = t279 * t883 + t281 * t550 + t283 * t551 + t419 * t445 + t420 * t447 + t443 * t653;
t36 = qJD(1) * t700 + t626 * t73 - t629 * t72;
t57 = (t629 * t423 + (t317 + t955) * qJD(1)) * t629 + (t316 * qJD(1) + (-t424 * t610 + t426 * t612 - t514 * t879 - t516 * t884 + t834) * t626 + (-t422 + t770 * t612 + t772 * t610 + (-t511 - t686) * qJD(1)) * t629) * t626;
t637 = (-t36 - t57 - t983) * t629 + t641;
t636 = t612 * t924 - t792 * t979 + t645;
t635 = -t968 * t884 / 0.2e1 + (qJD(1) * t969 + t915 - t916 + t917 - t918) * t939 + t986 * t938 + t985 * t937 + t983 * t787 + t984 * t786 + t980 * t781 + t979 * t780 + t989 * t967 + t988 * t966;
t634 = -t953 + (t977 - t982) * t787 + (t978 - t981) * t786 - t966 * t963 - t967 * t962;
t135 = t385 * t882 + t386 * t552 + t387 * t553 + t417 * t503 + t418 * t504 + t502 * t652;
t136 = t385 * t883 + t386 * t550 + t387 * t551 + t419 * t503 + t420 * t504 + t502 * t653;
t633 = -t918 / 0.2e1 + t913 / 0.2e1 + t915 / 0.2e1 - t914 / 0.2e1 - t916 / 0.2e1 + t917 / 0.2e1 + (t610 * t771 + t612 * t773 + t626 * t945 + t638 * t629 + t135 + t978) * t938 + (-t610 * t770 + t612 * t772 + t638 * t626 - t629 * t945 + t136 + t977) * t937 + (t513 * t612 + t515 * t610 - t560 * t629 - t626 * t683 + t230 + t298 + t963) * t781 + (t514 * t612 + t516 * t610 + t560 * t626 - t629 * t683 + t231 + t299 + t962) * t780;
t113 = -t298 * t612 + t610 * t700;
t114 = -t299 * t612 + t610 * t698;
t17 = (t619 * t698 - t135) * t612 + (qJD(1) * t699 + t299 * t619 + t626 * t70 + t629 * t71) * t610;
t18 = (t619 * t700 - t136) * t612 + (qJD(1) * t701 + t298 * t619 + t626 * t72 + t629 * t73) * t610;
t632 = t17 * t938 + t18 * t937 + t35 * t786 + t36 * t787 + t113 * t781 + t114 * t780 + (-t230 * t629 + t231 * t626) * t884 / 0.2e1 + (qJD(1) * t697 + t913 - t914) * t939 + t635 + t967 * t699 + t966 * t701;
t599 = pkin(2) * t789;
t574 = t743 * qJD(2);
t544 = -t824 + t837;
t543 = t626 * t743 - t919;
t507 = t779 * t629;
t506 = t779 * t626;
t497 = t929 + (pkin(1) - t920) * t629 + t837;
t496 = t626 * t675 + t616 + t919;
t476 = t518 + t768;
t475 = (rSges(4,3) - t631) * t629 + t673 * t626;
t459 = t626 * t656 + t833;
t458 = -qJD(1) * t536 + t629 * t656;
t451 = t845 * t626;
t434 = t960 + ((-rSges(3,3) - pkin(7)) * t626 + t675 * t629) * qJD(1);
t433 = (t616 + (-pkin(1) - t922) * t626) * qJD(1) + t644;
t431 = t748 * t626;
t389 = t438 * t882;
t379 = -t564 * t830 - t894 + (-t625 * t830 - t626 * t827) * pkin(2);
t378 = t564 * t831 + t599 + (-t533 - t822) * t629;
t372 = t407 * t882;
t344 = t537 * t626 - t629 * t684;
t343 = t536 * t626 - t956;
t342 = -t537 * t629 - t959;
t341 = -t536 * t629 - t626 * t685;
t336 = t564 * t876 + (t629 * t673 - t613) * qJD(1) + t838;
t335 = (-t606 - t921) * t831 + (-t665 - t964) * t629 + t840;
t332 = t768 - t851;
t331 = t640 + t740;
t325 = t794 * t626;
t322 = -t449 * t612 - t505 * t882;
t321 = t448 * t612 + t505 * t883;
t314 = t680 * t626;
t312 = -t409 * t612 - t493 * t882;
t310 = -t502 * t612 + (t504 * t627 - t895) * t610;
t309 = t674 + t768 + t409;
t308 = t650 + t737 + t839 - t868;
t306 = t310 * t884;
t305 = t421 + t846;
t304 = t690 * t610;
t301 = -t617 * t882 + t768 + t991;
t300 = (t580 - t631) * t629 + t648 * t626 + t734;
t293 = -t409 * t883 + t372;
t292 = qJD(1) * t452 + t626 * t855;
t291 = t629 * t855 + t847;
t287 = t754 * t626;
t268 = t666 * t626;
t261 = qJD(1) * t432 + t626 * t677;
t260 = t629 * t677 + t599 + t847;
t227 = -t518 * t831 + t801;
t216 = t610 * t776 + t774;
t215 = t311 + t852;
t204 = t266 + t846;
t197 = t655 * t830 + t812 * t936 + t575 - t741 + t838;
t196 = t576 + (-pkin(3) * t884 - t823) * t629 + t640 * qJD(1) + t799;
t178 = -t610 * t777 - t775;
t175 = t853 * t883 + t372 + t389;
t174 = qJD(1) * t326 + t626 * t803;
t173 = t629 * t803 + t795;
t172 = qJD(1) * t315 + t626 * t664;
t171 = t629 * t664 + t599 + t795;
t169 = -rSges(6,3) * t812 + (t629 * t669 - t603) * qJD(1) + t731 - t738 + t838;
t168 = qJD(1) * t650 + (-t605 * t884 + t612 * t670 - t964) * t629 + t793 + t804;
t167 = -t883 * t972 + t861;
t166 = (-t509 - t518) * t831 + t796 + t801;
t165 = t170 + t846;
t164 = t610 * t746 + t745;
t163 = t177 + t852;
t162 = (qJD(1) * t648 + t556) * t629 + (-qJD(1) * t580 + t954 + (t619 * t912 - t555) * t612) * t626 - t735 + t838;
t161 = (-t823 + (-t577 * t610 - t612 * t617) * t619) * t629 + (-t868 + (-rSges(7,3) * t610 + t769) * t626) * qJD(1) + t987;
t160 = (t505 * t876 + t286) * t612 + (t390 * t626 - t448 * t619 + t505 * t830) * t610;
t159 = (-t505 * t874 - t285) * t612 + (-t390 * t629 + t449 * t619 + t485) * t610;
t157 = t951 * t612 + (-t873 + (-t503 * t627 - t504 * t624) * qJD(4)) * t610 + t802;
t154 = -t407 * t884 + t761;
t153 = -t356 * t882 + (-t493 * t874 - t249) * t612 + t857;
t150 = qJD(1) * t288 + t626 * t759;
t149 = t629 * t759 + t755;
t142 = qJD(1) * t269 + t626 * t651;
t141 = t629 * t651 + t599 + t755;
t140 = t806 * t883 + t389 + t861;
t116 = t139 + t846;
t115 = t690 * t879 + (qJD(1) * t950 - t285 * t626 + t286 * t629) * t610;
t112 = t831 * t851 + t672;
t91 = -t409 * t812 + (qJD(1) * t952 - t249 * t626) * t610 + t865;
t80 = (-t509 + t851) * t831 + t672 + t796;
t75 = t854 * t884 + t760 + t761;
t74 = t858 * t882 + (t619 * t776 + t864) * t612 + t856 + t857;
t47 = -t884 * t973 + t647;
t46 = -t974 * t882 + (-t619 * t777 - t976) * t612 + t758;
t45 = t800 * t831 + t646;
t44 = (-t509 + t800) * t831 + t646 + t796;
t43 = t774 * t876 + (t864 * t626 + (t626 * t854 + t629 * t853) * qJD(1)) * t610 + t863 + t865;
t42 = t807 * t884 + t647 + t760;
t41 = t808 * t882 + (t619 * t746 + t809) * t612 + t758 + t856;
t40 = -t775 * t876 + (qJD(1) * t993 - t976 * t626) * t610 + t762;
t39 = t757 * t831 + t639;
t38 = t639 + (-t509 + t757) * t831 + t796;
t31 = t745 * t876 + (t809 * t626 + (t626 * t807 + t629 * t806) * qJD(1)) * t610 + t762 + t863;
t1 = [-t504 * t788 + t802 + (t161 * t301 + t162 * t300) * t940 + (t168 * t309 + t169 * t308) * t941 + (t196 * t332 + t197 * t331) * t942 + (t335 * t476 + t336 * t475) * t943 + (t433 * t497 + t434 * t496) * t944 - t561 * t884 + (-t723 + t730) * t828 + (t724 + t729) * t827 + (t531 + t951) * t612 + (t562 - t1007) * t879 + (-t503 * t825 + t609 * t998 + t886 * t996 + t532 - t873) * t610 - t990; (-qJD(2) * t685 + (qJD(1) * t539 - t626 * t657) * t628 + (qJD(1) * t541 - t626 * t658) * t625) * t937 + m(3) * ((-t433 * t626 - t434 * t629) * t586 + (-t496 * t629 - t497 * t626) * t574) + (t621 / 0.2e1 + t620 / 0.2e1) * t718 * qJD(2) + t633 + m(7) * (t141 * t300 + t142 * t301 + t161 * t268 + t162 * t269) + m(6) * (t168 * t314 + t169 * t315 + t171 * t308 + t172 * t309) + m(5) * (t196 * t431 + t197 * t432 + t260 * t331 + t261 * t332) + m(4) * (t335 * t506 + t336 * t507 + t378 * t475 + t379 * t476) + (-qJD(2) * t684 + (-qJD(1) * t538 - t629 * t657) * t628 + (-qJD(1) * t540 - t629 * t658) * t625) * t938 + ((t892 / 0.2e1 + t890 / 0.2e1 - t497 * t935) * t629 + (t496 * t935 + t893 / 0.2e1 + t891 / 0.2e1) * t626) * qJD(1); t626 * ((t626 * t458 + (t343 + t959) * qJD(1)) * t626 + (t344 * qJD(1) + (t538 * t827 + t540 * t828) * t629 + (-t459 + (-t890 - t892) * qJD(2) + (t537 - t685) * qJD(1)) * t626) * t629) + t641 - t629 * ((t629 * t459 + (t342 + t956) * qJD(1)) * t629 + (t341 * qJD(1) + (-t539 * t827 - t541 * t828 + t833) * t626 + (-t458 + (t891 + t893) * qJD(2) - t684 * qJD(1)) * t629) * t626) - t629 * t27 - t629 * t28 - t629 * t36 - t629 * t57 + (-t341 * t629 + t342 * t626) * t831 + (-t343 * t629 + t344 * t626) * t830 + ((t543 * t626 + t544 * t629) * ((qJD(1) * t543 + t644) * t629 + (-t960 + (-t544 - t824 + t614) * qJD(1)) * t626) + t835 * t586 * t574) * t944 + (t116 * t38 + t141 * t269 + t142 * t268) * t940 + (t165 * t44 + t171 * t315 + t172 * t314) * t941 + (t204 * t80 + t260 * t432 + t261 * t431) * t942 + (t166 * t305 + t378 * t507 + t379 * t506) * t943; m(4) * (-t475 * t629 - t476 * t626) * t533 + (-t335 * t626 - t336 * t629 + (t475 * t626 - t476 * t629) * qJD(1)) * t934 + m(7) * (t149 * t300 + t150 * t301 + t161 * t287 + t162 * t288) + m(6) * (t168 * t325 + t169 * t326 + t173 * t308 + t174 * t309) + m(5) * (t196 * t451 + t197 * t452 + t291 * t331 + t292 * t332) + t633; (-t378 * t629 - t379 * t626 + (-t506 * t629 + t507 * t626) * qJD(1)) * t934 + t637 + m(7) * (t116 * t39 + t139 * t38 + t141 * t288 + t142 * t287 + t149 * t269 + t150 * t268) + m(6) * (t165 * t45 + t170 * t44 + t171 * t326 + t172 * t325 + t173 * t315 + t174 * t314) + m(5) * (t112 * t204 + t260 * t452 + t261 * t451 + t266 * t80 + t291 * t432 + t292 * t431) + m(4) * (-t507 * t533 * t629 + t166 * t421 + t227 * t305 - t506 * t894); (t139 * t39 + t149 * t288 + t150 * t287) * t940 + (t170 * t45 + t173 * t326 + t174 * t325) * t941 + (t112 * t266 + t291 * t452 + t292 * t451) * t942 + t637 + (t533 * t564 * t835 + t227 * t421) * t943; t634 + m(7) * (t161 * t164 + t162 * t163 + t300 * t42 + t301 * t41) + m(6) * (t168 * t216 + t169 * t215 + t308 * t75 + t309 * t74) + m(5) * (t159 * t332 + t160 * t331 + t196 * t322 + t197 * t321) + ((t79 / 0.2e1 + t135 / 0.2e1) * t629 + (t78 / 0.2e1 + t136 / 0.2e1) * t626 + (-t626 * t783 + t629 * t784) * qJD(1)) * t610 + t306 + (-t157 + (t626 * t784 + t629 * t783) * t619 + t867) * t612; m(7) * (t116 * t31 + t140 * t38 + t141 * t163 + t142 * t164 + t268 * t41 + t269 * t42) + m(6) * (t165 * t43 + t171 * t215 + t172 * t216 + t175 * t44 + t314 * t74 + t315 * t75) + m(5) * (t115 * t204 + t159 * t431 + t160 * t432 + t260 * t321 + t261 * t322 + t304 * t80) + t632; t632 + m(7) * (t139 * t31 + t140 * t39 + t149 * t163 + t150 * t164 + t287 * t41 + t288 * t42) + m(6) * (t170 * t43 + t173 * t215 + t174 * t216 + t175 * t45 + t325 * t74 + t326 * t75) + m(5) * (t112 * t304 + t115 * t266 + t159 * t451 + t160 * t452 + t291 * t321 + t292 * t322); (t157 * t612 - t306 + (t626 * t113 + t629 * t114 - t612 * t697) * t619 + t924) * t612 + t645 + (t629 * t17 + t626 * t18 + t697 * t884 + (-t310 * t619 - t626 * t78 - t629 * t79) * t612 + ((-t230 * t612 + t113) * t629 + (t231 * t612 - t114 - t979) * t626) * qJD(1)) * t610 + (t140 * t31 + t163 * t42 + t164 * t41) * t940 + (t175 * t43 + t215 * t75 + t216 * t74) * t941 + (t115 * t304 + t159 * t322 + t160 * t321) * t942; t634 + m(7) * (t161 * t178 + t162 * t177 + t300 * t47 + t301 * t46) + m(6) * (t153 * t309 + t154 * t308 + t168 * t312 + t169 * t311) + t867 * t612; m(7) * (t116 * t40 + t141 * t177 + t142 * t178 + t167 * t38 + t268 * t46 + t269 * t47) + m(6) * (t153 * t314 + t154 * t315 + t165 * t91 + t171 * t311 + t172 * t312 + t293 * t44) + t635; m(7) * (t139 * t40 + t149 * t177 + t150 * t178 + t167 * t39 + t287 * t46 + t288 * t47) + m(6) * (t153 * t325 + t154 * t326 + t170 * t91 + t173 * t311 + t174 * t312 + t293 * t45) + t635; m(7) * (t140 * t40 + t163 * t47 + t164 * t46 + t167 * t31 + t177 * t42 + t178 * t41) + m(6) * (t153 * t216 + t154 * t215 + t175 * t91 + t293 * t43 + t311 * t75 + t312 * t74) + t636; t636 + (t167 * t40 + t177 * t47 + t178 * t46) * t940 + (t153 * t312 + t154 * t311 + t293 * t91) * t941; m(7) * ((t300 * t629 + t301 * t626) * t879 + (t161 * t626 + t162 * t629 + (-t300 * t626 + t301 * t629) * qJD(1)) * t610); m(7) * ((-t38 + (t268 * t626 + t269 * t629) * t619) * t612 + (t116 * t619 + t141 * t629 + t142 * t626 + (t268 * t629 - t269 * t626) * qJD(1)) * t610); m(7) * ((-t39 + (t287 * t626 + t288 * t629) * t619) * t612 + (t139 * t619 + t149 * t629 + t150 * t626 + (t287 * t629 - t288 * t626) * qJD(1)) * t610); m(7) * ((-t31 + (t163 * t629 + t164 * t626) * t619) * t612 + (t140 * t619 + t41 * t626 + t42 * t629 + (-t163 * t626 + t164 * t629) * qJD(1)) * t610); m(7) * ((-t40 + (t177 * t629 + t178 * t626) * t619) * t612 + (t167 * t619 + t46 * t626 + t47 * t629 + (-t177 * t626 + t178 * t629) * qJD(1)) * t610); (-0.1e1 + t835) * t610 * t879 * t940;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
