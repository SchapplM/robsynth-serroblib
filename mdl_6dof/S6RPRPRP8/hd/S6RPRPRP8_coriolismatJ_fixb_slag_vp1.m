% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPRP8_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_coriolismatJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP8_coriolismatJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP8_coriolismatJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP8_coriolismatJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP8_coriolismatJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:24:28
% EndTime: 2019-03-09 03:25:22
% DurationCPUTime: 45.33s
% Computational Cost: add. (71655->1097), mult. (97181->1493), div. (0->0), fcn. (105450->8), ass. (0->612)
t714 = qJ(3) + pkin(9);
t702 = sin(t714);
t703 = cos(t714);
t658 = pkin(4) * t703 + pkin(8) * t702;
t719 = sin(qJ(1));
t721 = cos(qJ(3));
t917 = t719 * t721;
t695 = pkin(3) * t917;
t859 = t719 * t658 + t695;
t717 = sin(qJ(5));
t720 = cos(qJ(5));
t801 = rSges(6,1) * t720 - rSges(6,2) * t717;
t560 = rSges(6,3) * t702 + t703 * t801;
t933 = t560 * t719;
t433 = t859 + t933;
t722 = cos(qJ(1));
t972 = pkin(3) * t721;
t812 = t658 + t972;
t435 = (t560 + t812) * t722;
t1110 = t433 * t719 + t722 * t435;
t970 = -qJ(4) - pkin(7);
t698 = t722 * t970;
t718 = sin(qJ(3));
t923 = t718 * t719;
t635 = pkin(3) * t923 - pkin(7) * t722 - t698;
t926 = t703 * t719;
t687 = pkin(8) * t926;
t929 = t702 * t719;
t860 = -pkin(4) * t929 - t635 + t687;
t973 = pkin(3) * t718;
t696 = t722 * t973;
t856 = t719 * t970 + t696;
t585 = t722 * (t719 * pkin(7) + t856);
t925 = t703 * t722;
t928 = t702 * t722;
t858 = pkin(4) * t928 - pkin(8) * t925;
t866 = -t722 * t858 - t585;
t1031 = rSges(7,1) + pkin(5);
t918 = t719 * t720;
t924 = t717 * t722;
t638 = t702 * t924 + t918;
t916 = t720 * t722;
t919 = t719 * t717;
t639 = t702 * t916 - t919;
t963 = rSges(7,3) + qJ(6);
t883 = -rSges(7,2) * t925 + t1031 * t639 + t638 * t963;
t636 = t702 * t919 - t916;
t637 = t702 * t918 + t924;
t884 = -rSges(7,2) * t926 + t1031 * t637 + t636 * t963;
t196 = -t883 * t722 + (t860 - t884) * t719 + t866;
t459 = t637 * rSges(6,1) - t636 * rSges(6,2) - rSges(6,3) * t926;
t463 = -t639 * rSges(6,1) + t638 * rSges(6,2) + rSges(6,3) * t925;
t242 = t463 * t722 + (-t459 + t860) * t719 + t866;
t384 = -t1031 * t636 + t637 * t963;
t880 = t1031 * t638 - t639 * t963;
t262 = -t719 * t384 + t722 * t880;
t501 = -rSges(6,1) * t636 - rSges(6,2) * t637;
t505 = rSges(6,1) * t638 + rSges(6,2) * t639;
t376 = t719 * t501 - t505 * t722;
t800 = rSges(7,1) * t720 + rSges(7,3) * t717;
t965 = rSges(7,2) * t702;
t559 = t703 * t800 + t965;
t799 = pkin(5) * t720 + qJ(6) * t717;
t871 = t799 * t703 + t559;
t809 = t871 * t719;
t387 = t809 + t859;
t389 = (t812 + t871) * t722;
t861 = (-t1031 * t717 + t720 * t963) * t703;
t472 = t861 * t719;
t473 = t861 * t722;
t618 = (-rSges(6,1) * t717 - rSges(6,2) * t720) * t703;
t1199 = -m(6) * (t1110 * t618 - t242 * t376) - m(7) * (t196 * t262 + t387 * t472 + t389 * t473);
t1067 = m(7) / 0.2e1;
t1071 = m(6) / 0.2e1;
t1182 = -t384 * t722 - t719 * t880;
t374 = -t501 * t722 - t719 * t505;
t901 = t1067 * t1182 + t374 * t1071;
t1163 = -t463 * t702 + t560 * t925;
t1173 = t702 * t883 + t871 * t925;
t282 = t702 * t884 + t703 * t809;
t367 = t459 * t702 + t560 * t926;
t905 = t1067 * (t1173 * t722 + t282 * t719) + t1071 * (t1163 * t722 + t367 * t719);
t21 = t905 - t901;
t1198 = qJD(1) * t21;
t900 = -t1067 * t262 + t376 * t1071;
t1072 = -m(6) / 0.2e1;
t1104 = t1163 * t719 - t367 * t722;
t276 = t282 * t722;
t906 = (-t1173 * t719 + t276) * t1067 + t1104 * t1072;
t20 = t906 - t900;
t1197 = qJD(1) * t20;
t704 = t722 * qJ(2);
t1131 = -t719 * pkin(1) + t704;
t764 = t1131 + t856;
t740 = t764 + t858;
t1193 = t740 + t883;
t813 = qJ(2) + t973;
t971 = pkin(4) * t702;
t974 = pkin(1) * t722;
t731 = t974 + (t813 + t971) * t719 - t687 - t698;
t322 = t731 + t884;
t194 = t1193 * t638 + t322 * t636;
t1196 = t194 * m(7) * qJD(1);
t1195 = Icges(5,3) + Icges(4,3);
t787 = Icges(5,5) * t702 + Icges(5,6) * t703;
t789 = Icges(4,5) * t718 + Icges(4,6) * t721;
t1194 = t787 + t789;
t609 = Icges(7,5) * t637;
t438 = -Icges(7,6) * t926 + Icges(7,3) * t636 + t609;
t444 = Icges(7,4) * t637 - Icges(7,2) * t926 + Icges(7,6) * t636;
t951 = Icges(7,5) * t636;
t450 = Icges(7,1) * t637 - Icges(7,4) * t926 + t951;
t237 = -t638 * t438 + t444 * t925 - t639 * t450;
t1192 = t237 * t722;
t441 = Icges(6,5) * t637 - Icges(6,6) * t636 - Icges(6,3) * t926;
t956 = Icges(6,4) * t637;
t447 = -Icges(6,2) * t636 - Icges(6,6) * t926 + t956;
t612 = Icges(6,4) * t636;
t453 = Icges(6,1) * t637 - Icges(6,5) * t926 - t612;
t239 = t441 * t925 + t638 * t447 - t639 * t453;
t1191 = t239 * t722;
t1190 = t719 * t237;
t1189 = t719 * t239;
t489 = -Icges(6,5) * t636 - Icges(6,6) * t637;
t491 = -Icges(7,4) * t636 + Icges(7,6) * t637;
t1188 = t489 + t491;
t490 = Icges(6,5) * t638 + Icges(6,6) * t639;
t492 = Icges(7,4) * t638 - Icges(7,6) * t639;
t1187 = t490 + t492;
t613 = Icges(6,4) * t638;
t454 = Icges(6,1) * t639 - Icges(6,5) * t925 - t613;
t885 = Icges(6,2) * t639 - t454 + t613;
t610 = Icges(7,5) * t638;
t452 = -Icges(7,1) * t639 + Icges(7,4) * t925 - t610;
t887 = Icges(7,3) * t639 + t452 - t610;
t1186 = t885 + t887;
t886 = -Icges(6,2) * t637 + t453 - t612;
t888 = -Icges(7,3) * t637 + t450 + t951;
t1185 = t886 + t888;
t614 = Icges(6,4) * t639;
t449 = Icges(6,2) * t638 + Icges(6,6) * t925 - t614;
t889 = -Icges(6,1) * t638 + t449 - t614;
t611 = Icges(7,5) * t639;
t439 = -Icges(7,6) * t925 + Icges(7,3) * t638 + t611;
t891 = Icges(7,1) * t638 - t439 - t611;
t1184 = t889 - t891;
t890 = Icges(6,1) * t636 + t447 + t956;
t892 = -Icges(7,1) * t636 + t438 + t609;
t1183 = t890 - t892;
t1161 = t449 * t717 + t454 * t720;
t443 = -Icges(6,5) * t639 + Icges(6,6) * t638 + Icges(6,3) * t925;
t939 = t443 * t702;
t267 = t1161 * t703 - t939;
t1162 = t439 * t717 - t452 * t720;
t446 = -Icges(7,4) * t639 + Icges(7,2) * t925 - Icges(7,6) * t638;
t936 = t446 * t702;
t264 = t1162 * t703 - t936;
t915 = t722 * t638;
t1106 = t636 * t719 + t915;
t1172 = m(7) * t1106;
t471 = -t1172 / 0.2e1;
t466 = t1172 / 0.2e1;
t1073 = m(5) / 0.2e1;
t656 = rSges(5,1) * t703 - rSges(5,2) * t702;
t804 = (t656 + t972) * t722;
t1115 = t719 * t804;
t630 = pkin(4) * t926 + pkin(8) * t929;
t834 = t695 + t630;
t840 = t703 * t918;
t842 = t703 * t919;
t878 = -rSges(7,2) * t929 - t1031 * t840 - t842 * t963;
t378 = t834 - t878;
t857 = -pkin(4) * t925 - pkin(8) * t928;
t379 = (t965 + t972 + (t1031 * t720 + t717 * t963) * t703) * t722 - t857;
t527 = rSges(6,1) * t840 - rSges(6,2) * t842 + rSges(6,3) * t929;
t422 = t527 + t834;
t423 = (t560 + t972) * t722 - t857;
t805 = rSges(5,1) * t926 - rSges(5,2) * t929;
t555 = t695 + t805;
t677 = rSges(4,1) * t721 - rSges(4,2) * t718;
t648 = t677 * t719;
t649 = t677 * t722;
t808 = (-t378 * t722 + t719 * t379) * t1067 + (-t422 * t722 + t719 * t423) * t1071 + (-t555 * t722 + t1115) * t1073 + m(4) * (-t648 * t722 + t719 * t649) / 0.2e1;
t566 = t656 * t719 + t695;
t838 = (t387 * t722 - t389 * t719) * t1067 + (t433 * t722 - t435 * t719) * t1071 + (t566 * t722 - t1115) * t1073;
t59 = t838 - t808;
t1181 = qJD(1) * t59;
t1114 = t722 * t804;
t383 = t722 * t389;
t836 = (-t387 * t719 - t383) * t1067 + t1110 * t1072 + (-t566 * t719 - t1114) * t1073;
t839 = (t719 * t378 + t379 * t722) * t1067 + (t719 * t422 + t423 * t722) * t1071 + (t719 * t555 + t1114) * t1073;
t63 = t839 - t836;
t1180 = t63 * qJD(1);
t777 = t438 * t717 + t450 * t720;
t938 = t444 * t702;
t263 = t703 * t777 + t938;
t774 = -t447 * t717 + t453 * t720;
t941 = t441 * t702;
t266 = t703 * t774 + t941;
t1179 = -t263 - t266;
t1178 = -t264 - t267;
t1140 = t1194 * t722 - t1195 * t719;
t958 = Icges(5,4) * t702;
t793 = Icges(5,2) * t703 + t958;
t573 = Icges(5,6) * t722 + t719 * t793;
t680 = Icges(5,4) * t926;
t575 = Icges(5,1) * t929 + Icges(5,5) * t722 + t680;
t960 = Icges(4,4) * t718;
t794 = Icges(4,2) * t721 + t960;
t592 = Icges(4,6) * t722 + t719 * t794;
t694 = Icges(4,4) * t917;
t594 = Icges(4,1) * t923 + Icges(4,5) * t722 + t694;
t1137 = t573 * t703 + t575 * t702 + t592 * t721 + t594 * t718;
t1177 = Icges(4,5) * t721 + Icges(5,5) * t703 - Icges(4,6) * t718 - Icges(5,6) * t702;
t776 = -t638 * t439 + t639 * t452;
t894 = t636 * t438 + t637 * t450;
t1176 = t894 + (-t444 * t719 - t446 * t722) * t703 + t776;
t773 = -t638 * t449 - t639 * t454;
t893 = -t636 * t447 + t637 * t453;
t1175 = t893 + (-t441 * t719 - t443 * t722) * t703 + t773;
t1105 = -t636 * t722 + t719 * t638;
t1159 = m(7) * t1105;
t468 = -t1159 / 0.2e1;
t464 = t1159 / 0.2e1;
t1171 = -t1183 * t637 - t1185 * t636 - t1188 * t926;
t1170 = -t1184 * t637 - t1186 * t636 - t1187 * t926;
t1169 = t1183 * t639 + t1185 * t638 + t1188 * t925;
t1168 = t1184 * t639 + t1186 * t638 + t1187 * t925;
t600 = (-Icges(7,4) * t717 + Icges(7,6) * t720) * t703;
t950 = Icges(7,5) * t717;
t795 = Icges(7,1) * t720 + t950;
t552 = Icges(7,4) * t702 + t703 * t795;
t596 = (Icges(7,3) * t720 - t950) * t703;
t874 = t552 - t596;
t949 = Icges(7,5) * t720;
t679 = t703 * t949;
t927 = t703 * t717;
t946 = Icges(7,6) * t702;
t544 = Icges(7,3) * t927 + t679 + t946;
t604 = -Icges(7,1) * t927 + t679;
t876 = t544 + t604;
t227 = -t600 * t926 - t636 * t874 + t637 * t876;
t597 = (-Icges(6,5) * t717 - Icges(6,6) * t720) * t703;
t955 = Icges(6,4) * t717;
t796 = Icges(6,1) * t720 - t955;
t554 = Icges(6,5) * t702 + t703 * t796;
t601 = (-Icges(6,2) * t720 - t955) * t703;
t873 = t554 + t601;
t954 = Icges(6,4) * t720;
t792 = -Icges(6,2) * t717 + t954;
t550 = Icges(6,6) * t702 + t703 * t792;
t605 = (-Icges(6,1) * t717 - t954) * t703;
t875 = t550 - t605;
t228 = -t597 * t926 - t636 * t873 - t637 * t875;
t1167 = t228 + t227;
t229 = t600 * t925 + t638 * t874 - t639 * t876;
t230 = t597 * t925 + t638 * t873 + t639 * t875;
t1166 = t230 + t229;
t574 = -Icges(5,6) * t719 + t722 * t793;
t957 = Icges(5,4) * t703;
t797 = Icges(5,1) * t702 + t957;
t576 = -Icges(5,5) * t719 + t722 * t797;
t593 = -Icges(4,6) * t719 + t722 * t794;
t959 = Icges(4,4) * t721;
t798 = Icges(4,1) * t718 + t959;
t595 = -Icges(4,5) * t719 + t722 * t798;
t1164 = (-t574 * t703 - t576 * t702 - t593 * t721 - t595 * t718) * t722;
t234 = -t439 * t636 - t446 * t926 + t637 * t452;
t236 = -t443 * t926 - t636 * t449 - t454 * t637;
t1069 = -m(7) / 0.2e1;
t1035 = t719 / 0.2e1;
t791 = Icges(7,4) * t720 + Icges(7,6) * t717;
t548 = Icges(7,2) * t702 + t703 * t791;
t330 = -t638 * t544 + t548 * t925 - t639 * t552;
t1158 = t330 * t702;
t786 = Icges(6,5) * t720 - Icges(6,6) * t717;
t546 = Icges(6,3) * t702 + t703 * t786;
t331 = t546 * t925 + t638 * t550 - t639 * t554;
t1157 = t331 * t702;
t802 = rSges(5,1) * t702 + rSges(5,2) * t703;
t1154 = t722 * t802;
t785 = Icges(7,3) * t717 + t949;
t1089 = t703 * t785 + t946;
t512 = t1089 * t719;
t520 = t552 * t719;
t752 = -t548 * t719 + t777;
t188 = (t512 * t717 + t520 * t720 + t444) * t703 - t752 * t702;
t518 = t550 * t719;
t522 = t554 * t719;
t750 = t546 * t719 - t774;
t190 = (-t518 * t717 + t522 * t720 + t441) * t703 + t750 * t702;
t1153 = -t188 - t190;
t513 = t1089 * t722;
t521 = t552 * t722;
t751 = t548 * t722 - t1162;
t189 = (-t513 * t717 - t521 * t720 + t446) * t703 - t751 * t702;
t519 = t550 * t722;
t523 = t554 * t722;
t749 = -t546 * t722 + t1161;
t191 = (t519 * t717 - t523 * t720 + t443) * t703 + t749 * t702;
t1152 = t189 + t191;
t204 = t491 * t702 + (-t717 * t888 + t720 * t892) * t703;
t206 = t489 * t702 + (-t717 * t886 - t720 * t890) * t703;
t1151 = -t204 - t206;
t205 = t492 * t702 + (-t717 * t887 + t720 * t891) * t703;
t207 = t490 * t702 + (-t717 * t885 - t720 * t889) * t703;
t1150 = t205 + t207;
t547 = Icges(7,2) * t703 - t702 * t791;
t770 = t717 * t544 + t720 * t552;
t748 = -t547 + t770;
t934 = t548 * t702;
t1094 = t703 * t748 + t934;
t543 = Icges(7,6) * t703 - t702 * t785;
t551 = Icges(7,4) * t703 - t702 * t795;
t212 = t1094 * t719 + t543 * t636 + t551 * t637;
t545 = Icges(6,3) * t703 - t702 * t786;
t769 = -t717 * t550 + t720 * t554;
t747 = t545 - t769;
t935 = t546 * t702;
t1095 = t703 * t747 - t935;
t549 = Icges(6,6) * t703 - t702 * t792;
t553 = Icges(6,5) * t703 - t702 * t796;
t213 = -t1095 * t719 - t549 * t636 + t553 * t637;
t1149 = t212 + t213;
t214 = -t1094 * t722 - t638 * t543 - t639 * t551;
t215 = t1095 * t722 + t638 * t549 - t639 * t553;
t1148 = t214 + t215;
t220 = (t543 * t717 + t551 * t720 + t548) * t703 - t748 * t702;
t221 = (-t549 * t717 + t553 * t720 + t546) * t703 + t747 * t702;
t1147 = t220 + t221;
t326 = t544 * t636 - t548 * t926 + t552 * t637;
t327 = -t546 * t926 - t550 * t636 + t554 * t637;
t1146 = t326 + t327;
t1145 = t330 + t331;
t346 = t703 * t770 + t934;
t347 = t703 * t769 + t935;
t1144 = t346 + t347;
t1143 = t573 * t926 + t575 * t929 + t592 * t917 + t594 * t923 + (t1194 * t719 + t1195 * t722) * t722;
t1142 = -t1140 * t722 - t574 * t926 - t576 * t929 - t593 * t917 - t595 * t923;
t1141 = -t1140 * t719 - t1164;
t1139 = t1177 * t719;
t1138 = t1177 * t722;
t1136 = t1178 * t722 + t1179 * t719;
t672 = -Icges(4,2) * t718 + t959;
t1132 = (t798 / 0.2e1 + t672 / 0.2e1) * t721;
t1130 = t1137 * t722;
t602 = -Icges(5,2) * t929 + t680;
t652 = -Icges(5,2) * t702 + t957;
t603 = t652 * t722;
t654 = Icges(5,1) * t703 - t958;
t606 = t654 * t719;
t607 = t654 * t722;
t644 = -Icges(4,2) * t923 + t694;
t645 = t672 * t722;
t674 = Icges(4,1) * t721 - t960;
t646 = t674 * t719;
t647 = t674 * t722;
t1127 = -(t719 * (t574 - t607) - t722 * (t573 - t606)) * t702 + (t719 * (t576 + t603) - t722 * (t575 + t602)) * t703 - (t719 * (t593 - t647) - t722 * (t592 - t646)) * t718 + (t719 * (t595 + t645) - t722 * (t594 + t644)) * t721;
t1126 = -t702 / 0.2e1;
t1041 = t702 / 0.2e1;
t1125 = -t703 / 0.2e1;
t1039 = t703 / 0.2e1;
t1037 = -t719 / 0.2e1;
t1033 = t722 / 0.2e1;
t1122 = (t1170 * t722 - t1171 * t719) * t703 + t1167 * t702;
t1121 = (t1168 * t722 - t1169 * t719) * t703 + t1166 * t702;
t1120 = t234 * t719;
t1119 = t234 * t722;
t1118 = t236 * t719;
t1117 = t236 * t722;
t1116 = t560 * t722;
t1100 = t740 - t463;
t372 = t459 + t731;
t1111 = -t1100 * t719 + t722 * t372;
t1103 = t191 / 0.2e1 + t189 / 0.2e1;
t1102 = t190 / 0.2e1 + t188 / 0.2e1;
t1098 = (-t1111 * t618 - t433 * t505 - t435 * t501) * t1072 + (t1193 * t472 - t322 * t473 - t384 * t389 - t387 * t880) * t1069;
t872 = rSges(7,2) * t703 + (-t799 - t800) * t702;
t192 = (t719 * t872 + t884) * t703 + (-t878 - t809) * t702;
t877 = -t559 * t722 - t799 * t925;
t193 = (t722 * t872 + t883) * t703 + (-t722 * t871 - t877) * t702;
t558 = rSges(6,3) * t703 - t702 * t801;
t280 = (t558 * t719 + t459) * t703 + (t527 - t933) * t702;
t281 = (t558 * t722 - t463) * t703;
t1097 = (t1173 * t379 + t1193 * t193 + t192 * t322 + t282 * t378) * t1069 + (t1100 * t281 + t1163 * t423 + t280 * t372 + t367 * t422) * t1072;
t1093 = t703 * t749 - t939;
t1092 = t703 * t750 - t941;
t1091 = t703 * t751 + t936;
t1090 = t703 * t752 + t938;
t816 = t554 / 0.2e1 + t552 / 0.2e1;
t818 = t550 / 0.2e1 - t544 / 0.2e1;
t1088 = -t717 * (t601 / 0.2e1 - t596 / 0.2e1 + t816) + t720 * (t605 / 0.2e1 + t604 / 0.2e1 - t818);
t1083 = -t719 * t883 + t722 * t884;
t1082 = t717 * (t549 / 0.2e1 - t543 / 0.2e1) - t720 * (t553 / 0.2e1 + t551 / 0.2e1) - t546 / 0.2e1 - t548 / 0.2e1 + t652 / 0.2e1 + t797 / 0.2e1;
t1081 = t717 * t818 - t720 * t816 + t545 / 0.2e1 + t547 / 0.2e1 + t793 / 0.2e1 - t654 / 0.2e1;
t715 = t719 ^ 2;
t716 = t722 ^ 2;
t693 = t715 + t716;
t1080 = 0.2e1 * t693;
t1078 = 0.4e1 * qJD(1);
t1077 = 2 * qJD(3);
t1075 = 2 * qJD(5);
t1074 = 4 * qJD(5);
t1066 = -pkin(1) - pkin(7);
t771 = t459 * t722 + t463 * t719;
t222 = (t1116 * t719 - t527 * t722) * t703 + t771 * t702;
t340 = t771 * t703;
t1064 = m(6) * (t1163 * t281 - t222 * t340 + t280 * t367);
t153 = (-t719 * t877 + t722 * t878) * t703 + t1083 * t702;
t218 = t1083 * t703;
t261 = t1173 * t842;
t1061 = m(7) * (-t638 * t192 + t636 * t193 + t261 + (t218 * t702 + (t153 - t276) * t703) * t717);
t1060 = m(7) * (t1173 * t193 - t153 * t218 + t192 * t282);
t436 = t1105 * t703;
t701 = t703 ^ 2;
t510 = t636 * t702 + t701 * t919;
t511 = t702 * t638 + t701 * t924;
t841 = t703 * t924;
t1058 = m(7) * (t1106 * t218 + t436 * t196 - t282 * t841 + t511 * t387 - t510 * t389 + t261);
t1054 = m(7) * (t387 * t637 + t389 * t639 + t472 * t636 + t473 * t638 + (t196 * t720 + t262 * t717) * t703);
t1052 = m(7) * (t1173 * t638 + t1193 * t511 + t282 * t636 + t322 * t510);
t1043 = t693 / 0.2e1;
t1034 = t719 / 0.4e1;
t1032 = t722 / 0.4e1;
t1030 = m(3) * ((rSges(3,3) * t722 + t1131) * t722 + (t974 + (rSges(3,3) + qJ(2)) * t719) * t719);
t803 = rSges(4,1) * t718 + rSges(4,2) * t721;
t735 = -t719 * rSges(4,3) + t722 * t803;
t534 = t1066 * t719 + t704 + t735;
t535 = (rSges(4,3) - t1066) * t722 + (qJ(2) + t803) * t719;
t1029 = m(4) * (t534 * t649 + t535 * t648);
t1028 = m(4) * (t534 * t722 + t535 * t719);
t736 = -t719 * rSges(5,3) + t1154;
t476 = t736 + t764;
t477 = -t698 + (rSges(5,3) + pkin(1)) * t722 + (t802 + t813) * t719;
t1026 = m(5) * (t476 * t804 + t477 * t555);
t1025 = m(5) * (-t476 * t719 + t722 * t477);
t1024 = m(5) * (t476 * t722 + t719 * t477);
t1011 = m(6) * (t1100 * t423 + t372 * t422);
t1010 = m(6) * (-t1100 * t505 + t372 * t501);
t1007 = m(6) * t1111;
t1006 = m(6) * (t1100 * t722 + t719 * t372);
t997 = m(7) * (t1193 * t637 - t322 * t639 - t384 * t638 - t636 * t880);
t762 = t1193 * t842 - t322 * t841;
t996 = m(7) * (-t638 * t378 + t636 * t379 + t762);
t995 = m(7) * (t387 * t638 - t636 * t389 + t762);
t994 = m(7) * (t1193 * t379 + t322 * t378);
t993 = m(7) * (-t1193 * t880 + t322 * t384);
t990 = m(7) * (-t1193 * t719 + t722 * t322);
t989 = m(7) * (t1193 * t722 + t719 * t322);
t968 = m(7) * qJD(3);
t967 = m(7) * qJD(5);
t966 = m(7) * qJD(6);
t238 = t446 * t925 - t776;
t914 = t238 + t1176;
t240 = t443 * t925 - t773;
t912 = t240 + t1175;
t233 = -t444 * t926 + t894;
t909 = -t233 + t1176;
t235 = -t441 * t926 + t893;
t907 = -t235 + t1175;
t855 = qJD(1) * t702;
t854 = qJD(1) * t703;
t853 = qJD(3) * t719;
t852 = qJD(3) * t722;
t851 = qJD(5) * t702;
t850 = qJD(5) * t703;
t400 = (t1069 + t1072 - m(5) / 0.2e1) * t1080;
t849 = t400 * qJD(1);
t166 = t1090 * t719 + t512 * t636 + t520 * t637;
t167 = t1091 * t719 - t513 * t636 - t521 * t637;
t168 = -t1092 * t719 - t518 * t636 + t522 * t637;
t169 = -t1093 * t719 + t519 * t636 - t523 * t637;
t782 = -t235 * t719 + t1117;
t783 = -t233 * t719 + t1119;
t846 = (-t782 - t783 + t1149) * t1126 + ((t167 + t169) * t722 + (-t166 - t168) * t719 + t1146) * t1125;
t170 = -t1090 * t722 - t638 * t512 - t639 * t520;
t171 = -t1091 * t722 + t638 * t513 + t639 * t521;
t172 = t1092 * t722 + t638 * t518 - t639 * t522;
t173 = t1093 * t722 - t638 * t519 + t639 * t523;
t780 = t240 * t722 - t1189;
t781 = t238 * t722 - t1190;
t845 = (-t780 - t781 + t1148) * t1041 + ((t171 + t173) * t722 + (-t170 - t172) * t719 + t1145) * t1039;
t844 = (-t1136 + t1147) * t1126 + (t1152 * t722 + t1153 * t719 + t1144) * t1125;
t843 = t1033 * t1171 + t1035 * t1170;
t831 = -t926 / 0.4e1;
t829 = t925 / 0.4e1;
t828 = t1033 * t1169 + t1035 * t1168;
t313 = t326 * t702;
t113 = t703 * t783 + t313;
t314 = t327 * t702;
t114 = t703 * t782 + t314;
t827 = t114 / 0.2e1 + t113 / 0.2e1;
t115 = t703 * t781 + t1158;
t116 = t703 * t780 + t1157;
t826 = -t116 / 0.2e1 - t115 / 0.2e1;
t825 = t1125 * t1136 + t1126 * t1144;
t137 = t233 * t722 + t1120;
t138 = t235 * t722 + t1118;
t824 = t137 / 0.2e1 + t138 / 0.2e1;
t139 = t238 * t719 + t1192;
t140 = t240 * t719 + t1191;
t823 = -t139 / 0.2e1 - t140 / 0.2e1;
t822 = t204 / 0.2e1 + t206 / 0.2e1;
t821 = t205 / 0.2e1 + t207 / 0.2e1;
t815 = t597 / 0.2e1 + t600 / 0.2e1;
t814 = -t787 / 0.2e1 - t789 / 0.2e1;
t810 = t803 * t693;
t806 = t693 * t972;
t726 = (-t280 * t722 + t281 * t719) * t1071 + (-t192 * t722 + t193 * t719) * t1067;
t727 = (t472 * t719 + t473 * t722) * t1067 + m(6) * t618 * t1043;
t61 = -t726 + t727;
t725 = (t280 * t719 + t281 * t722) * t1072 + (t192 * t719 + t193 * t722) * t1069;
t739 = (t472 * t722 - t719 * t473) * t1067;
t67 = t739 + t725;
t784 = t61 * qJD(2) + t67 * qJD(4);
t737 = (t637 * t719 + t639 * t722) * t1067;
t742 = m(7) * (-t510 * t722 + t511 * t719);
t317 = t737 - t742 / 0.2e1;
t738 = (t637 * t722 - t639 * t719) * t1067;
t741 = m(7) * (t510 * t719 + t511 * t722);
t318 = t738 - t741 / 0.2e1;
t763 = t317 * qJD(2) + t318 * qJD(4);
t761 = -t843 - t846;
t760 = t722 * t857 - t806;
t759 = t828 - t845;
t41 = -t1158 + (t722 * t909 + t1190) * t703;
t42 = -t1157 + (t722 * t907 + t1189) * t703;
t758 = -t41 / 0.2e1 - t42 / 0.2e1 + t826;
t39 = t313 + (-t719 * t914 + t1119) * t703;
t40 = t314 + (-t719 * t912 + t1117) * t703;
t757 = t40 / 0.2e1 + t39 / 0.2e1 - t827;
t734 = -t1098 + (t1150 + t1166) * t1034 + (-t1151 + t1167) * t1032;
t53 = t719 * t909 - t1192;
t54 = t719 * t907 - t1191;
t730 = t53 / 0.2e1 + t54 / 0.2e1 - t823 + t1140 * t715 / 0.2e1 + t1141 * t1035 + ((t1137 + t1140) * t722 - t1130 + t1142) * t1033;
t51 = t722 * t914 + t1120;
t52 = t722 * t912 + t1118;
t729 = t51 / 0.2e1 + t52 / 0.2e1 - t824 - t1143 * t722 / 0.2e1 + t1142 * t1037 + (t1130 + t1142) * t1035 + ((-t1137 + t1140) * t719 + t1141 + t1143 + t1164) * t1033;
t724 = -(t137 + t138) * t925 / 0.4e1 + (t51 + t52) * t829 - (t113 + t114) * t719 / 0.4e1 + (t39 + t40) * t1034 + (t139 + t140 + t53 + t54) * t831 + (t115 + t116 + t41 + t42) * t1032;
t723 = -t1097 + t1147 * t1041 + t1144 * t1039 + (t1146 - t1179) * t929 / 0.4e1 - (t1145 + t1178) * t928 / 0.4e1 + (t1149 - t1153) * t831 + (t1148 + t1152) * t829;
t657 = pkin(8) * t703 - t971;
t640 = t719 * t657;
t581 = t693 * t927;
t567 = t696 + t1154;
t565 = (-t802 - t973) * t719;
t536 = t636 * t842;
t434 = t696 + (-t558 - t657) * t722;
t432 = t640 + (t558 - t973) * t719;
t399 = (-m(7) / 0.4e1 - m(6) / 0.4e1 - m(5) / 0.4e1) * t1080 + (m(5) + m(6) + m(7)) * t1043;
t398 = t536 + (-t702 * t717 + t915) * t927;
t395 = -t702 * t505 + t618 * t925;
t394 = t501 * t702 + t618 * t926;
t388 = t696 + (-t657 - t872) * t722;
t386 = t640 + (t872 - t973) * t719;
t377 = t701 * t717 * t720 + t636 * t637 + t638 * t639;
t373 = t387 * t842;
t351 = t374 * t703;
t320 = t741 / 0.2e1 + t738;
t319 = t742 / 0.2e1 + t737;
t309 = -t1116 * t722 + (-t527 - t630) * t719 + t760;
t303 = 0.2e1 * t471;
t302 = 0.2e1 * t468;
t301 = t466 + t471;
t300 = 0.2e1 * t466;
t299 = t464 + t468;
t298 = 0.2e1 * t464;
t294 = -t702 * t880 + t861 * t925;
t293 = t384 * t702 + t472 * t703;
t292 = (t702 * t597 + (-t717 * t873 - t720 * t875) * t703) * t702;
t291 = (t702 * t600 + (-t717 * t874 + t720 * t876) * t703) * t702;
t246 = t877 * t722 + (-t630 + t878) * t719 + t760;
t243 = t1182 * t703;
t122 = t995 / 0.2e1;
t120 = t996 / 0.2e1;
t118 = t997 / 0.2e1;
t117 = -t1106 * t196 + t389 * t841 + t373;
t108 = t1007 + t990 + t1025;
t97 = t1173 * t511 - t218 * t436 + t282 * t510;
t95 = t172 * t722 + t173 * t719;
t94 = t170 * t722 + t171 * t719;
t93 = t168 * t722 + t169 * t719;
t92 = t166 * t722 + t167 * t719;
t84 = t1006 + t1024 + t1028 + t989 + t1030;
t82 = t1052 / 0.2e1;
t79 = t1054 / 0.2e1;
t68 = t1088 * t703 + t815 * t702 + t1010 + t993;
t66 = t739 - t725;
t64 = t836 + t839;
t60 = t726 + t727;
t57 = t808 + t838;
t55 = t1058 / 0.2e1;
t27 = -t1132 + (-t674 / 0.2e1 + t794 / 0.2e1) * t718 + t1029 + t1026 + t1011 + t994 - t1082 * t703 + t1081 * t702;
t25 = t1061 / 0.2e1;
t24 = t120 - t995 / 0.2e1;
t23 = t122 + t120;
t22 = t122 - t996 / 0.2e1;
t18 = t900 + t906;
t17 = t901 + t905;
t13 = t82 + t118;
t12 = t82 - t997 / 0.2e1;
t11 = t118 - t1052 / 0.2e1;
t10 = t55 + t79 - t1061 / 0.2e1;
t9 = t55 + t25 - t1054 / 0.2e1;
t8 = t79 + t25 - t1058 / 0.2e1;
t7 = t719 * t828 + t722 * t843 - t1199;
t6 = (t719 * t758 + t722 * t757) * t703;
t5 = t719 * t729 + t722 * t730;
t4 = t1064 + t1060 + (t719 * t846 + t722 * t845 - t825) * t703 + (t719 * t827 + t722 * t826 - t844) * t702;
t3 = t723 + t734 + (t114 / 0.4e1 + t113 / 0.4e1 - t40 / 0.4e1 - t39 / 0.4e1) * t719 + (-t116 / 0.4e1 - t115 / 0.4e1 - t42 / 0.4e1 - t41 / 0.4e1) * t722 + ((-t52 / 0.4e1 - t51 / 0.4e1 + t138 / 0.4e1 + t137 / 0.4e1) * t722 + (t54 / 0.4e1 + t53 / 0.4e1 + t140 / 0.4e1 + t139 / 0.4e1) * t719) * t703;
t2 = t724 + t723 + (-t230 / 0.4e1 - t229 / 0.4e1 - t207 / 0.4e1 - t205 / 0.4e1) * t719 + (-t228 / 0.4e1 - t227 / 0.4e1 - t206 / 0.4e1 - t204 / 0.4e1) * t722 + t1098;
t1 = t724 + (-t347 / 0.2e1 - t346 / 0.2e1 + (-t215 / 0.4e1 - t214 / 0.4e1 - t191 / 0.4e1 - t189 / 0.4e1) * t722 + (t213 / 0.4e1 + t212 / 0.4e1 + t190 / 0.4e1 + t188 / 0.4e1) * t719) * t703 + t734 + (-t221 / 0.2e1 - t220 / 0.2e1 + (t331 / 0.4e1 + t330 / 0.4e1 - t267 / 0.4e1 - t264 / 0.4e1) * t722 + (-t327 / 0.4e1 - t326 / 0.4e1 - t266 / 0.4e1 - t263 / 0.4e1) * t719) * t702 + t1097;
t14 = [t84 * qJD(2) + t27 * qJD(3) + t108 * qJD(4) + t68 * qJD(5) + t194 * t966, qJD(1) * t84 + qJD(3) * t57 + qJD(4) * t399 + qJD(5) * t17 + qJD(6) * t299, t27 * qJD(1) + t57 * qJD(2) + t64 * qJD(4) + t3 * qJD(5) + t23 * qJD(6) + ((t476 * t565 + t477 * t567 + (-t555 + t566) * t804) * t1073 + (t1100 * t432 + t372 * t434 - t422 * t435 + t423 * t433) * t1071 + (t1193 * t386 + t322 * t388 - t378 * t389 + t379 * t387) * t1067) * t1077 + (t213 / 0.2e1 + t212 / 0.2e1 + (-t592 / 0.2e1 + t646 / 0.2e1) * t721 + (-t594 / 0.2e1 - t644 / 0.2e1) * t718 + (-t573 / 0.2e1 + t606 / 0.2e1) * t703 + (-t575 / 0.2e1 - t602 / 0.2e1) * t702 + t814 * t722 + m(4) * (t535 * t803 - t648 * t677) - t730 + t1102) * t852 + ((-t607 / 0.2e1 + t574 / 0.2e1) * t703 + (t603 / 0.2e1 + t576 / 0.2e1) * t702 + t814 * t719 + t215 / 0.2e1 + t214 / 0.2e1 + (-t647 / 0.2e1 + t593 / 0.2e1) * t721 + (t645 / 0.2e1 + t595 / 0.2e1) * t718 + m(4) * (-t534 * t803 + t649 * t677) - t729 + t1103) * t853, qJD(1) * t108 + qJD(2) * t399 + qJD(3) * t64 + qJD(5) * t18 + qJD(6) * t301, t68 * qJD(1) + t17 * qJD(2) + t3 * qJD(3) + t18 * qJD(4) + (t292 + t291) * qJD(5) + t13 * qJD(6) + ((-t1173 * t880 + t1193 * t294 + t282 * t384 + t293 * t322) * t1067 + (t1100 * t395 - t1163 * t505 + t367 * t501 + t372 * t394) * t1071) * t1075 + ((t229 / 0.2e1 + t230 / 0.2e1 - t757 + t821) * t722 + (-t228 / 0.2e1 - t227 / 0.2e1 - t758 - t822) * t719) * t850, t299 * qJD(2) + t23 * qJD(3) + t301 * qJD(4) + t13 * qJD(5) + t1196; -t59 * qJD(3) + t400 * qJD(4) - t21 * qJD(5) + t298 * qJD(6) + (-t989 / 0.4e1 - t1006 / 0.4e1 - t1024 / 0.4e1 - t1028 / 0.4e1 - t1030 / 0.4e1) * t1078, 0, -t1181 + t60 * qJD(5) + (-m(4) * t810 / 0.2e1 + (t565 * t719 - t567 * t722) * t1073 + (t432 * t719 - t434 * t722) * t1071 + (t386 * t719 - t388 * t722) * t1067) * t1077 + t581 * t966, t849, -t1198 + t60 * qJD(3) + ((-t394 * t722 + t395 * t719) * t1071 + (-t293 * t722 + t294 * t719) * t1067) * t1075 + t319 * qJD(6), t298 * qJD(1) + t319 * qJD(5) + t581 * t968; t59 * qJD(2) + t5 * qJD(3) - t63 * qJD(4) + t1 * qJD(5) + t22 * qJD(6) + (-t1029 / 0.4e1 - t1026 / 0.4e1 - t1011 / 0.4e1 - t994 / 0.4e1) * t1078 + t1082 * t854 - t1081 * t855 + (t1132 + (-t794 + t674) * t718 / 0.2e1) * qJD(1), qJD(5) * t61 + t1181, t5 * qJD(1) + t7 * qJD(5) + t117 * t966 + (m(7) * (t196 * t246 + t386 * t387 - t388 * t389) + m(6) * (t242 * t309 + t432 * t433 - t434 * t435) + m(5) * (t566 * t565 - t804 * t567 + (-t722 * t736 - t585 + (-t722 * rSges(5,3) - t719 * t802 - t635) * t719) * (-t656 * t716 - t719 * t805 - t806)) + m(4) * ((-t722 * t735 + (-t722 * rSges(4,3) - t719 * t803) * t719) * (-t719 * t648 - t649 * t722) - t677 * t810) + (t94 + t95 - t1138 * t715 + (t1139 * t719 + t1127) * t722) * t1035 + (t92 + t93 + t1139 * t716 + (-t1138 * t722 - t1127) * t719) * t1033) * qJD(3), qJD(5) * t67 - t1180, t1 * qJD(1) + t7 * qJD(3) + (t1033 * t1122 + t1035 * t1121) * qJD(5) + t10 * qJD(6) + (-t1064 / 0.4e1 - t1060 / 0.4e1) * t1074 + ((t1104 * t618 + t351 * t242 + t340 * t376 - t394 * t435 + t395 * t433) * t1071 + (t1173 * t472 + t196 * t243 - t218 * t262 - t282 * t473 - t293 * t389 + t294 * t387) * t1067) * t1075 + ((t822 - t826) * t722 + (t821 - t827) * t719 + t844) * t851 + (t719 * t761 + t722 * t759 + t825) * t850 + t784, t22 * qJD(1) + t117 * t968 + t10 * qJD(5) + (t536 + (-t1106 + t915) * t927 - t398) * t966; -t400 * qJD(2) + t63 * qJD(3) - t20 * qJD(5) + t300 * qJD(6) + (-t990 / 0.4e1 - t1007 / 0.4e1 - t1025 / 0.4e1) * t1078, -t849, t1180 + t66 * qJD(5) + ((t386 * t722 + t719 * t388) * t1067 + (t432 * t722 + t719 * t434) * t1071 + (t565 * t722 + t719 * t567) * t1073) * t1077, 0, -t1197 + t66 * qJD(3) + ((t394 * t719 + t395 * t722) * t1071 + (t293 * t719 + t294 * t722) * t1067) * t1075 + t320 * qJD(6), qJD(1) * t300 + qJD(5) * t320; t21 * qJD(2) + t2 * qJD(3) + t20 * qJD(4) + t6 * qJD(5) + t12 * qJD(6) - t815 * t855 + (-t993 / 0.4e1 - t1010 / 0.4e1) * t1078 - t1088 * t854, -qJD(3) * t61 - qJD(6) * t317 + t1198, t2 * qJD(1) + t4 * qJD(5) + t9 * qJD(6) + ((t1173 * t386 + t153 * t196 - t192 * t389 + t193 * t387 - t218 * t246 + t282 * t388) * t1067 + (t1163 * t432 + t222 * t242 - t280 * t435 + t281 * t433 - t309 * t340 + t367 * t434) * t1071) * t1077 + t761 * t852 - t759 * t853 - t784 + (((t266 / 0.2e1 + t263 / 0.2e1 + t94 / 0.2e1 + t95 / 0.2e1) * t722 + (-t93 / 0.2e1 - t92 / 0.2e1 - t264 / 0.2e1 - t267 / 0.2e1) * t719) * t703 + ((t823 + t1102) * t722 + (t824 + t1103) * t719) * t702 + t1199) * qJD(3), -qJD(3) * t67 - qJD(6) * t318 + t1197, t6 * qJD(1) + t4 * qJD(3) + (t292 / 0.2e1 + t291 / 0.2e1) * t851 + (m(7) * (t1173 * t294 - t218 * t243 + t282 * t293) / 0.4e1 + (t1163 * t395 - t340 * t351 + t367 * t394) * m(6) / 0.4e1) * t1074 + t97 * t966 + ((t1150 * t722 + t1151 * t719) * t1041 + t1122 * t1037 + t1121 * t1033) * t850, t12 * qJD(1) + t9 * qJD(3) + t97 * t967 + (t436 * t927 - t510 * t638 + t511 * t636 - t377) * t966 - t763; t302 * qJD(2) + t24 * qJD(3) + t303 * qJD(4) + t11 * qJD(5) - t1196, qJD(1) * t302 + qJD(5) * t317, t24 * qJD(1) + (t636 * t386 - t638 * t388 + t373 + (-t196 * t702 + (t246 + t383) * t703) * t717 - t117) * t968 + t8 * qJD(5) + t398 * t966, qJD(1) * t303 + qJD(5) * t318, t11 * qJD(1) + t8 * qJD(3) + (-t282 * t639 + t1173 * t637 - t293 * t638 + t294 * t636 + (-t218 * t720 + t243 * t717) * t703 - t97) * t967 + t377 * t966 + t763, 0.4e1 * (t398 * qJD(3) / 0.4e1 + t377 * qJD(5) / 0.4e1) * m(7);];
Cq  = t14;
