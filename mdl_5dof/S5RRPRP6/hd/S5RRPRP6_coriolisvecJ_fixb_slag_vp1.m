% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP6_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP6_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP6_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP6_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:56:46
% EndTime: 2019-12-31 19:58:15
% DurationCPUTime: 76.45s
% Computational Cost: add. (27275->1081), mult. (40162->1420), div. (0->0), fcn. (37163->8), ass. (0->532)
t1048 = -Icges(5,4) - Icges(6,4);
t983 = Icges(5,1) + Icges(6,1);
t1044 = Icges(5,5) + Icges(6,5);
t1005 = -Icges(5,2) - Icges(6,2);
t1004 = Icges(5,6) + Icges(6,6);
t1047 = Icges(3,3) + Icges(4,3);
t486 = cos(qJ(4));
t1046 = t1048 * t486;
t483 = sin(qJ(4));
t1045 = t1048 * t483;
t1043 = Icges(5,3) + Icges(6,3);
t478 = qJ(2) + pkin(8);
t460 = sin(t478);
t461 = cos(t478);
t484 = sin(qJ(2));
t487 = cos(qJ(2));
t1013 = Icges(3,5) * t487 + Icges(4,5) * t461 - Icges(3,6) * t484 - Icges(4,6) * t460;
t1042 = -t1004 * t483 + t1044 * t486;
t1041 = t1005 * t483 - t1046;
t1040 = t983 * t486 + t1045;
t488 = cos(qJ(1));
t1039 = t1047 * t488;
t1038 = t1042 * t461 + t1043 * t460;
t1015 = -t1004 * t461 + t1041 * t460;
t930 = t1004 * t460 + t1041 * t461;
t991 = t1040 * t460 - t1044 * t461;
t929 = t1040 * t461 + t1044 * t460;
t485 = sin(qJ(1));
t773 = t485 * t487;
t776 = t484 * t485;
t779 = t461 * t485;
t781 = t460 * t485;
t1001 = -Icges(3,5) * t773 - Icges(4,5) * t779 + Icges(3,6) * t776 + Icges(4,6) * t781 + t1039;
t1014 = t1013 * t488 + t1047 * t485;
t1037 = (-t1004 * t486 - t1044 * t483) * t460;
t1036 = (t1005 * t486 + t1045) * t460;
t1035 = (-t983 * t483 + t1046) * t460;
t1029 = Icges(3,5) * t484 + Icges(4,5) * t460 + Icges(3,6) * t487 + Icges(4,6) * t461;
t807 = Icges(4,6) * t488;
t304 = Icges(4,4) * t779 - Icges(4,2) * t781 - t807;
t808 = Icges(3,6) * t488;
t326 = Icges(3,4) * t773 - Icges(3,2) * t776 - t808;
t1034 = t304 * t460 + t326 * t484;
t426 = Icges(4,4) * t781;
t813 = Icges(4,5) * t488;
t306 = Icges(4,1) * t779 - t426 - t813;
t446 = Icges(3,4) * t776;
t814 = Icges(3,5) * t488;
t328 = Icges(3,1) * t773 - t446 - t814;
t989 = -t306 * t461 - t328 * t487 + t1034;
t982 = t1001 * t488 - t485 * t989;
t715 = qJD(4) * t460;
t717 = qJD(2) * t485;
t381 = t488 * t715 + t717;
t716 = qJD(2) * t488;
t382 = -t485 * t715 + t716;
t714 = qJD(4) * t461;
t433 = qJD(1) - t714;
t770 = t488 * t483;
t774 = t485 * t486;
t368 = -t461 * t770 + t774;
t772 = t486 * t488;
t777 = t483 * t485;
t369 = t461 * t772 + t777;
t780 = t460 * t488;
t993 = t1042 * t460 - t1043 * t461;
t908 = t1015 * t368 + t991 * t369 + t993 * t780;
t344 = Icges(6,4) * t368;
t188 = Icges(6,1) * t369 + Icges(6,5) * t780 + t344;
t347 = Icges(5,4) * t368;
t191 = Icges(5,1) * t369 + Icges(5,5) * t780 + t347;
t1016 = t188 + t191;
t817 = Icges(6,4) * t369;
t182 = Icges(6,2) * t368 + Icges(6,6) * t780 + t817;
t820 = Icges(5,4) * t369;
t185 = Icges(5,2) * t368 + Icges(5,6) * t780 + t820;
t1017 = t182 + t185;
t176 = Icges(6,5) * t369 + Icges(6,6) * t368 + Icges(6,3) * t780;
t179 = Icges(5,5) * t369 + Icges(5,6) * t368 + Icges(5,3) * t780;
t994 = t176 + t179;
t944 = t1016 * t369 + t1017 * t368 + t994 * t780;
t367 = t461 * t774 - t770;
t343 = Icges(6,4) * t367;
t366 = t461 * t777 + t772;
t180 = -Icges(6,2) * t366 + Icges(6,6) * t781 + t343;
t346 = Icges(5,4) * t367;
t183 = -Icges(5,2) * t366 + Icges(5,6) * t781 + t346;
t1018 = t180 + t183;
t174 = Icges(6,5) * t367 - Icges(6,6) * t366 + Icges(6,3) * t781;
t177 = Icges(5,5) * t367 - Icges(5,6) * t366 + Icges(5,3) * t781;
t1019 = t174 + t177;
t342 = Icges(6,4) * t366;
t187 = -Icges(6,1) * t367 - Icges(6,5) * t781 + t342;
t345 = Icges(5,4) * t366;
t190 = -Icges(5,1) * t367 - Icges(5,5) * t781 + t345;
t904 = t187 + t190;
t945 = t1018 * t368 + t1019 * t780 - t369 * t904;
t951 = t381 * t944 - t945 * t382 + t908 * t433;
t909 = -t1015 * t366 + t991 * t367 + t993 * t781;
t946 = t1016 * t367 - t1017 * t366 + t994 * t781;
t947 = -t1018 * t366 + t1019 * t781 - t367 * t904;
t952 = t381 * t946 - t947 * t382 + t909 * t433;
t821 = Icges(4,4) * t460;
t390 = Icges(4,1) * t461 - t821;
t307 = Icges(4,5) * t485 + t390 * t488;
t822 = Icges(3,4) * t484;
t416 = Icges(3,1) * t487 - t822;
t329 = Icges(3,5) * t485 + t416 * t488;
t1030 = -t307 * t779 - t329 * t773;
t387 = Icges(4,2) * t461 + t821;
t450 = Icges(4,4) * t461;
t389 = Icges(4,1) * t460 + t450;
t413 = Icges(3,2) * t487 + t822;
t470 = Icges(3,4) * t487;
t415 = Icges(3,1) * t484 + t470;
t1012 = t387 * t460 - t389 * t461 + t413 * t484 - t415 * t487;
t578 = t486 * t433;
t646 = qJD(1) * t461 - qJD(4);
t685 = t460 * t716;
t886 = t485 * t646 + t685;
t158 = t483 * t886 + t488 * t578;
t577 = t433 * t483;
t159 = -t486 * t886 + t488 * t577;
t684 = t461 * t716;
t720 = qJD(1) * t485;
t691 = t460 * t720;
t540 = t684 - t691;
t85 = Icges(6,5) * t159 + Icges(6,6) * t158 + Icges(6,3) * t540;
t87 = Icges(5,5) * t159 + Icges(5,6) * t158 + Icges(5,3) * t540;
t1028 = t85 + t87;
t686 = t460 * t717;
t160 = t485 * t578 + (-t488 * t646 + t686) * t483;
t718 = qJD(2) * t460;
t161 = t646 * t772 + (-t486 * t718 + t577) * t485;
t719 = qJD(1) * t488;
t541 = t460 * t719 + t461 * t717;
t86 = Icges(6,5) * t161 + Icges(6,6) * t160 + Icges(6,3) * t541;
t88 = Icges(5,5) * t161 + Icges(5,6) * t160 + Icges(5,3) * t541;
t1027 = t86 + t88;
t89 = Icges(6,4) * t159 + Icges(6,2) * t158 + Icges(6,6) * t540;
t91 = Icges(5,4) * t159 + Icges(5,2) * t158 + Icges(5,6) * t540;
t1026 = t89 + t91;
t90 = Icges(6,4) * t161 + Icges(6,2) * t160 + Icges(6,6) * t541;
t92 = Icges(5,4) * t161 + Icges(5,2) * t160 + Icges(5,6) * t541;
t1025 = t90 + t92;
t93 = Icges(6,1) * t159 + Icges(6,4) * t158 + Icges(6,5) * t540;
t95 = Icges(5,1) * t159 + Icges(5,4) * t158 + Icges(5,5) * t540;
t1024 = t93 + t95;
t94 = Icges(6,1) * t161 + Icges(6,4) * t160 + Icges(6,5) * t541;
t96 = Icges(5,1) * t161 + Icges(5,4) * t160 + Icges(5,5) * t541;
t1023 = t94 + t96;
t1022 = t1038 * qJD(2) + t1037 * qJD(4);
t1021 = t930 * qJD(2) + t1036 * qJD(4);
t1020 = t929 * qJD(2) + t1035 * qJD(4);
t1011 = -t1015 * t483 + t991 * t486;
t1010 = t1014 * t488 + t1030;
t771 = t487 * t488;
t778 = t461 * t488;
t891 = -t1014 * t485 - t307 * t778 - t329 * t771;
t1009 = t1001 * t485 - t306 * t778 - t328 * t771;
t918 = t1029 * t488;
t917 = t1029 * t485;
t599 = -Icges(4,2) * t460 + t450;
t305 = Icges(4,6) * t485 + t488 * t599;
t600 = -Icges(3,2) * t484 + t470;
t327 = Icges(3,6) * t485 + t488 * t600;
t1008 = t305 * t460 + t327 * t484;
t981 = -t305 * t781 - t327 * t776 - t1010;
t775 = t484 * t488;
t980 = -t304 * t780 - t326 * t775 - t1009;
t979 = -t305 * t780 - t327 * t775 - t891;
t1003 = -t1012 * t485 - t918;
t1002 = -t1012 * t488 + t917;
t937 = t304 * t461 + t306 * t460 + t326 * t487 + t328 * t484;
t936 = t305 * t461 + t307 * t460 + t327 * t487 + t329 * t484;
t999 = t1029 * qJD(2);
t998 = t307 * t461 + t329 * t487 - t1008;
t958 = t1018 * t158 + t1019 * t540 + t1023 * t369 + t1025 * t368 + t1027 * t780 - t904 * t159;
t957 = t1016 * t159 + t1017 * t158 + t1024 * t369 + t1026 * t368 + t1028 * t780 + t994 * t540;
t956 = t1018 * t160 + t1019 * t541 + t1023 * t367 - t1025 * t366 + t1027 * t781 - t904 * t161;
t955 = t1016 * t161 + t1017 * t160 + t1024 * t367 - t1026 * t366 + t1028 * t781 + t994 * t541;
t915 = t1015 * t158 + t1020 * t369 + t1021 * t368 + t1022 * t780 + t991 * t159 + t993 * t540;
t997 = t1015 * t160 + t1020 * t367 - t1021 * t366 + t1022 * t781 + t991 * t161 + t993 * t541;
t594 = -t180 * t483 - t187 * t486;
t73 = -t174 * t461 + t460 * t594;
t592 = -t183 * t483 - t190 * t486;
t75 = -t177 * t461 + t460 * t592;
t996 = t73 + t75;
t593 = -t182 * t483 + t188 * t486;
t74 = -t176 * t461 + t460 * t593;
t591 = -t185 * t483 + t191 * t486;
t76 = -t179 * t461 + t460 * t591;
t995 = t74 + t76;
t907 = t1011 * t460 - t993 * t461;
t371 = t599 * qJD(2);
t372 = t390 * qJD(2);
t399 = t600 * qJD(2);
t400 = t416 * qJD(2);
t990 = -t371 * t460 + t372 * t461 - t399 * t484 + t400 * t487 + (-t387 * t461 - t389 * t460 - t413 * t487 - t415 * t484) * qJD(2) + t1029 * qJD(1);
t988 = t1014 * qJD(1);
t987 = (-t1011 + t1038) * t433 + (t993 * t485 + t592 + t594) * t382 + (-t993 * t488 - t591 - t593) * t381;
t986 = t1012 * qJD(1) + t1013 * qJD(2);
t985 = t381 * (t1005 * t369 + t1016 + t344 + t347) - t382 * (t1005 * t367 - t342 - t345 - t904) + t433 * (t991 + t1036);
t984 = (t1011 * qJD(2) - t1022) * t461 + (t1020 * t486 - t1021 * t483 + (-t1015 * t486 - t991 * t483) * qJD(4) + t993 * qJD(2)) * t460;
t978 = t1002 * qJD(1);
t455 = pkin(4) * t486 + pkin(3);
t977 = -rSges(6,1) * t367 + rSges(6,2) * t366 - t455 * t779;
t710 = qJD(5) * t488;
t418 = t460 * t710;
t702 = pkin(4) * t770;
t858 = pkin(3) * t461;
t481 = -qJ(5) - pkin(7);
t852 = pkin(7) + t481;
t911 = t460 * t852;
t768 = t702 + (t858 + t911) * t485 - rSges(6,3) * t781 + t977;
t976 = t433 * t768 + t418;
t552 = qJD(2) * t387;
t168 = -t488 * t552 + (-t485 * t599 + t807) * qJD(1);
t554 = qJD(2) * t389;
t170 = -t488 * t554 + (-t390 * t485 + t813) * qJD(1);
t553 = qJD(2) * t413;
t224 = -t488 * t553 + (-t485 * t600 + t808) * qJD(1);
t555 = qJD(2) * t415;
t226 = -t488 * t555 + (-t416 * t485 + t814) * qJD(1);
t975 = -t936 * qJD(2) - t168 * t460 + t170 * t461 - t224 * t484 + t226 * t487 + t988;
t169 = qJD(1) * t305 - t485 * t552;
t171 = qJD(1) * t307 - t485 * t554;
t225 = qJD(1) * t327 - t485 * t553;
t227 = qJD(1) * t329 - t485 * t555;
t974 = t1001 * qJD(1) + t937 * qJD(2) + t169 * t460 - t171 * t461 + t225 * t484 - t227 * t487;
t973 = -t1037 * t433 + (-t1004 * t367 - t1044 * t366) * t382 + (t1004 * t369 - t1044 * t368) * t381;
t972 = (t979 * t485 - t980 * t488) * qJD(2);
t971 = (t981 * t485 - t982 * t488) * qJD(2);
t970 = t1003 * qJD(1);
t966 = t989 * qJD(1) - t999 * t485 + t988;
t965 = -t999 * t488 + (-t1013 * t485 + t1039 - t998) * qJD(1);
t962 = 0.2e1 * qJD(2);
t705 = qJD(2) * qJD(4);
t671 = t461 * t705;
t273 = qJD(1) * t381 + t485 * t671;
t274 = qJD(1) * t382 + t488 * t671;
t672 = t460 * t705;
t961 = t945 * t273 + t944 * t274 + t957 * t381 - t958 * t382 + t915 * t433 + t908 * t672;
t960 = t947 * t273 + t946 * t274 + t955 * t381 - t956 * t382 + t997 * t433 + t909 * t672;
t959 = rSges(6,1) + pkin(4);
t17 = (qJD(2) * t594 - t86) * t461 + (qJD(2) * t174 - t483 * t90 + t486 * t94 + (-t180 * t486 + t187 * t483) * qJD(4)) * t460;
t19 = (qJD(2) * t592 - t88) * t461 + (qJD(2) * t177 - t483 * t92 + t486 * t96 + (-t183 * t486 + t190 * t483) * qJD(4)) * t460;
t954 = t17 + t19;
t18 = (qJD(2) * t593 - t85) * t461 + (qJD(2) * t176 - t483 * t89 + t486 * t93 + (-t182 * t486 - t188 * t483) * qJD(4)) * t460;
t20 = (qJD(2) * t591 - t87) * t461 + (qJD(2) * t179 - t483 * t91 + t486 * t95 + (-t185 * t486 - t191 * t483) * qJD(4)) * t460;
t953 = t18 + t20;
t950 = t381 * t995 - t382 * t996 + t433 * t907;
t949 = t970 + t971;
t948 = t972 + t978;
t943 = t485 * t986 + t488 * t990;
t942 = t485 * t990 - t488 * t986;
t941 = qJD(2) * t989 - t169 * t461 - t171 * t460 - t225 * t487 - t227 * t484;
t940 = t998 * qJD(2) + t168 * t461 + t170 * t460 + t224 * t487 + t226 * t484;
t308 = rSges(4,1) * t779 - rSges(4,2) * t781 - t488 * rSges(4,3);
t471 = t485 * rSges(4,3);
t309 = rSges(4,1) * t778 - rSges(4,2) * t780 + t471;
t584 = t308 * t485 + t309 * t488;
t476 = t488 * pkin(6);
t431 = pkin(1) * t485 - t476;
t482 = -qJ(3) - pkin(6);
t453 = t488 * t482;
t859 = pkin(2) * t487;
t456 = pkin(1) + t859;
t728 = -t485 * t456 - t453;
t300 = t431 + t728;
t475 = t485 * pkin(6);
t432 = t488 * pkin(1) + t475;
t437 = t488 * t456;
t649 = -t482 * t485 + t437;
t301 = t649 - t432;
t752 = -t300 * t717 + t301 * t716;
t105 = qJD(2) * t584 + t752;
t751 = -t485 * t300 + t488 * t301;
t544 = t584 + t751;
t938 = qJD(2) * t544 + t105;
t935 = t1015 * t485;
t934 = t1015 * t488;
t933 = t991 * t485;
t932 = t991 * t488;
t853 = pkin(3) - t455;
t263 = -t461 * t853 - t911;
t619 = rSges(6,1) * t486 - rSges(6,2) * t483;
t845 = rSges(6,3) * t460;
t931 = t461 * t619 + t263 + t845;
t928 = t987 * t460;
t528 = t326 * t488 - t327 * t485;
t529 = t304 * t488 - t305 * t485;
t877 = t485 * (-t387 * t488 + t307) - t488 * (-Icges(4,2) * t779 + t306 - t426);
t878 = t485 * (-t413 * t488 + t329) - t488 * (-Icges(3,2) * t773 + t328 - t446);
t927 = -t460 * t877 + t529 * t461 - t484 * t878 + t528 * t487;
t730 = t415 + t600;
t731 = -t413 + t416;
t734 = t389 + t599;
t735 = -t387 + t390;
t926 = (-t460 * t734 + t461 * t735 - t484 * t730 + t487 * t731) * qJD(1);
t674 = t852 * t461;
t925 = rSges(6,3) * t461 - t674 + (-t619 + t853) * t460;
t924 = -t1019 * t382 + t381 * t994 + t433 * t993;
t898 = t485 * t996 + t488 * t995;
t923 = t485 * t995 - t488 * t996;
t897 = t485 * t945 + t488 * t944;
t922 = t485 * t944 - t488 * t945;
t896 = t485 * t947 + t946 * t488;
t921 = t946 * t485 - t488 * t947;
t435 = pkin(3) * t778;
t363 = pkin(7) * t780 + t435;
t396 = pkin(3) * t460 - pkin(7) * t461;
t364 = t396 * t717;
t860 = pkin(2) * t484;
t440 = t717 * t860;
t726 = qJD(3) * t488 + t440;
t743 = t301 + t432;
t527 = (t363 + t743) * qJD(1) - t364 - t726;
t711 = qJD(5) * t485;
t679 = t460 * t711;
t903 = t369 * rSges(6,1) + t368 * rSges(6,2) + rSges(6,3) * t780 + pkin(4) * t777 + t455 * t778;
t767 = -t481 * t780 - t363 + t903;
t45 = t381 * t925 + t433 * t767 + t527 + t679;
t656 = t45 * t925;
t397 = pkin(7) * t460 + t858;
t361 = t397 * t485;
t644 = t361 * t717 + t363 * t716 + t752;
t712 = qJD(5) * t461;
t33 = -t381 * t768 + t382 * t767 + t644 - t712;
t658 = t33 * t768;
t920 = -t656 + t658;
t919 = t1013 * qJD(1);
t623 = rSges(5,1) * t367 - rSges(5,2) * t366;
t193 = rSges(5,3) * t781 + t623;
t622 = rSges(5,1) * t486 - rSges(5,2) * t483;
t293 = -rSges(5,3) * t461 + t460 * t622;
t916 = -t193 * t433 - t293 * t382;
t405 = qJD(1) * t431;
t902 = qJD(1) * t300 - t405;
t900 = (-t1015 + t1035) * t433 + (t366 * t983 + t1018 + t343 + t346) * t382 + (t368 * t983 - t1017 - t817 - t820) * t381;
t899 = t973 * t460;
t843 = pkin(4) * qJD(4);
t699 = t486 * t843;
t893 = t159 * rSges(6,1) + t158 * rSges(6,2) + rSges(6,3) * t684 + qJD(1) * t702 + t481 * t691 + t485 * t699 + t418;
t892 = t1001 + t1008;
t890 = t433 * t984 + t907 * t672;
t410 = pkin(7) * t684;
t542 = -t461 * t720 - t685;
t220 = pkin(3) * t542 - pkin(7) * t691 + t410;
t409 = pkin(3) * t686;
t221 = pkin(7) * t541 + qJD(1) * t435 - t409;
t707 = qJD(1) * qJD(2);
t673 = t488 * t707;
t459 = pkin(6) * t719;
t462 = qJD(3) * t485;
t683 = t484 * t716;
t575 = -pkin(2) * t683 + t462;
t854 = pkin(1) - t456;
t218 = -t459 + (t485 * t854 - t453) * qJD(1) + t575;
t692 = t482 * t720 + t726;
t219 = (-t488 * t854 - t475) * qJD(1) - t692;
t696 = t218 * t716 + t219 * t717 - t300 * t673;
t571 = t220 * t716 + t221 * t717 + t361 * t673 + t696;
t744 = t301 + t363;
t629 = t744 * t720;
t713 = qJD(5) * t460;
t621 = rSges(6,1) * t161 + rSges(6,2) * t160;
t783 = t455 * t460;
t850 = t409 + (qJD(1) * t263 - t699) * t488 + (t713 + pkin(4) * t577 + (-t674 - t783) * qJD(2)) * t485 + rSges(6,3) * t541 + t621;
t700 = t483 * t843;
t851 = -t410 + (pkin(7) * t720 + t716 * t853) * t460 + ((-qJD(2) * t481 - t700) * t488 + t853 * t720) * t461 - rSges(6,3) * t691 + t893;
t5 = t851 * t382 + t850 * t381 - t768 * t274 - t767 * t273 + (-t629 + t713) * qJD(2) + t571;
t879 = t33 * t851 + t5 * t767;
t874 = -m(6) / 0.2e1;
t873 = m(6) / 0.2e1;
t872 = t273 / 0.2e1;
t871 = t274 / 0.2e1;
t870 = -t381 / 0.2e1;
t869 = t381 / 0.2e1;
t868 = -t382 / 0.2e1;
t867 = t382 / 0.2e1;
t866 = -t433 / 0.2e1;
t865 = t433 / 0.2e1;
t863 = t485 / 0.2e1;
t862 = -t488 / 0.2e1;
t861 = -rSges(5,3) - pkin(7);
t857 = pkin(4) * t483;
t849 = rSges(3,1) * t487;
t848 = rSges(4,1) * t461;
t847 = rSges(5,3) * t460;
t842 = t17 * t382;
t841 = t18 * t381;
t840 = t19 * t382;
t839 = t20 * t381;
t295 = t461 * t622 + t847;
t351 = (-rSges(5,1) * t483 - rSges(5,2) * t486) * t460;
t163 = qJD(2) * t295 + qJD(4) * t351;
t197 = t369 * rSges(5,1) + t368 * rSges(5,2) + rSges(5,3) * t780;
t374 = t397 * qJD(2);
t395 = qJD(1) * (-pkin(1) * t720 + t459);
t706 = qJD(1) * qJD(3);
t694 = qJD(1) * t218 + t485 * t706 + t395;
t703 = qJD(2) ^ 2 * t859;
t526 = qJD(1) * t220 - t485 * t703 + t694;
t659 = -t396 - t860;
t627 = t659 * t488;
t572 = qJD(1) * t627;
t697 = t159 * rSges(5,1) + t158 * rSges(5,2) + rSges(5,3) * t684;
t98 = -rSges(5,3) * t691 + t697;
t25 = -t163 * t381 - t274 * t293 + t433 * t98 + (t197 * t715 - t374 * t485 + t572) * qJD(2) + t526;
t838 = t25 * t488;
t624 = rSges(5,1) * t161 + rSges(5,2) * t160;
t100 = rSges(5,3) * t541 + t624;
t729 = qJD(1) * t440 + t488 * t706;
t543 = qJD(1) * t364 - t488 * t703 + t729;
t404 = t432 * qJD(1);
t760 = -t219 - t404;
t693 = -t221 + t760;
t26 = -t100 * t433 - t163 * t382 + t273 * t293 + (-t193 * t715 - t374 * t488) * qJD(2) + t693 * qJD(1) + t543;
t837 = t26 * t485;
t472 = t485 * rSges(3,3);
t66 = t197 * t433 - t293 * t381 + t527;
t832 = t485 * t66;
t831 = t73 * t273;
t830 = t74 * t274;
t829 = t75 * t273;
t828 = t76 * t274;
t827 = -rSges(6,3) + t481;
t391 = rSges(4,1) * t460 + rSges(4,2) * t461;
t661 = -t391 - t860;
t628 = t488 * t661;
t573 = qJD(2) * t628;
t539 = t462 + t573;
t746 = t300 - t431;
t112 = (-t308 + t746) * qJD(1) + t539;
t798 = t112 * t391;
t725 = rSges(3,2) * t776 + t488 * rSges(3,3);
t358 = rSges(3,1) * t773 - t725;
t420 = rSges(3,1) * t484 + rSges(3,2) * t487;
t687 = t420 * t716;
t204 = -t687 + (-t358 - t431) * qJD(1);
t797 = t204 * t485;
t796 = t204 * t488;
t688 = t420 * t717;
t359 = rSges(3,1) * t771 - rSges(3,2) * t775 + t472;
t737 = t359 + t432;
t205 = qJD(1) * t737 - t688;
t384 = t420 * t488;
t795 = t205 * t384;
t782 = t460 * t481;
t350 = (-rSges(6,1) * t483 - rSges(6,2) * t486) * t460;
t769 = qJD(2) * t931 + qJD(4) * t350 - t460 * t700 - t712;
t762 = t925 * t485;
t759 = -rSges(6,2) * t367 - t366 * t959;
t758 = -rSges(6,2) * t369 + t368 * t959;
t360 = t396 * t485;
t362 = t396 * t488;
t740 = -t360 * t717 - t362 * t716;
t690 = t484 * t720;
t441 = pkin(2) * t690;
t736 = t396 * t720 + t441;
t732 = rSges(4,2) * t691 + rSges(4,3) * t719;
t727 = rSges(3,2) * t690 + rSges(3,3) * t719;
t704 = pkin(2) * t775;
t701 = qJD(2) * t859;
t695 = t488 * t218 + t485 * t219 - t300 * t719;
t689 = t391 * t717;
t682 = t487 * t716;
t675 = -pkin(1) - t849;
t669 = t719 / 0.2e1;
t668 = t718 / 0.2e1;
t667 = -t717 / 0.2e1;
t666 = t717 / 0.2e1;
t664 = t716 / 0.2e1;
t392 = -rSges(4,2) * t460 + t848;
t660 = -t392 - t859;
t538 = qJD(2) * t627 + t462;
t507 = (-t361 + t746) * qJD(1) + t538;
t44 = t382 * t925 + t507 + t976;
t657 = t44 * t925;
t654 = (-t485 ^ 2 - t488 ^ 2) * t484;
t648 = qJD(1) * t360 + t441;
t647 = -t374 + t712;
t643 = t485 * t361 + t488 * t363 + t751;
t635 = t460 * t857 - t350;
t634 = -t293 + t659;
t631 = qJD(4) * t668;
t630 = -t374 - t701;
t625 = -rSges(3,2) * t484 + t849;
t618 = -t44 * t488 - t45 * t485;
t65 = t507 + t916;
t609 = t488 * t65 + t832;
t590 = t193 * t488 - t197 * t485;
t589 = -t205 * t485 - t796;
t579 = t659 + t925;
t576 = (-t362 - t704) * qJD(1);
t574 = -t163 + t630;
t570 = t488 * t220 + t485 * t221 + t361 * t719 + t695;
t373 = t392 * qJD(2);
t569 = -qJD(2) * t373 - t703;
t568 = -t455 * t461 - t456 - t845;
t383 = t420 * t485;
t352 = t391 * t485;
t557 = t630 - t769;
t198 = (t358 * t485 + t359 * t488) * qJD(2);
t545 = -t397 - t847;
t537 = t33 * t850 - t5 * t768;
t531 = t44 * t768 + t45 * t767;
t530 = -t33 * t767 - t657;
t513 = -qJD(1) * t361 + t538 + t902;
t52 = t193 * t381 + t197 * t382 + t644;
t502 = t52 * t590 + (t485 * t65 - t488 * t66) * t293;
t495 = t530 * t485 - t488 * t920;
t402 = t625 * qJD(2);
t353 = t391 * t488;
t261 = t293 * t488;
t259 = t293 * t485;
t245 = rSges(5,1) * t368 - rSges(5,2) * t369;
t243 = -rSges(5,1) * t366 - rSges(5,2) * t367;
t229 = -qJD(2) * t383 + (t488 * t625 + t472) * qJD(1);
t228 = -rSges(3,2) * t682 + (-t487 * t720 - t683) * rSges(3,1) + t727;
t173 = -qJD(2) * t352 + (t392 * t488 + t471) * qJD(1);
t172 = rSges(4,1) * t542 - rSges(4,2) * t684 + t732;
t119 = -t402 * t716 + (-t229 - t404 + t688) * qJD(1);
t118 = -t402 * t717 + t395 + (t228 - t687) * qJD(1);
t113 = -t689 + (t309 + t743) * qJD(1) - t726;
t70 = t569 * t488 + (-t173 + t689 + t760) * qJD(1) + t729;
t69 = t569 * t485 + (t172 + t573) * qJD(1) + t694;
t16 = -qJD(2) * t629 + t100 * t381 + t193 * t274 - t197 * t273 + t382 * t98 + t571;
t7 = -t850 * t433 - t769 * t382 - t925 * t273 + (t488 * t647 + t715 * t768) * qJD(2) + (-t679 + t693) * qJD(1) + t543;
t6 = qJD(1) * t418 + t851 * t433 - t769 * t381 + t925 * t274 + (t485 * t647 + t715 * t767 + t572) * qJD(2) + t526;
t1 = [t951 * t867 + (t7 * (t728 + t977) + t44 * (-t621 + t692) + t6 * (t437 + t903) + t45 * (t462 + t893) + (-t6 * t782 + t45 * (-t461 * t481 - t783 - t860) * qJD(2) + (t7 * t483 + (-t45 * t461 * t483 + t44 * t486) * qJD(4)) * pkin(4) + (t44 * (t568 + t782) - t45 * t482) * qJD(1)) * t488 + (-t6 * t482 + t44 * (qJD(2) * t827 + t700) * t461 + (t7 * t827 + t44 * (qJD(2) * t455 - qJD(5))) * t460 + (-t44 * t857 + t45 * t568) * qJD(1)) * t485 - (-t44 + t513 + t976) * t45 - t656 * t382) * m(6) + (-t1012 * qJD(2) + t371 * t461 + t372 * t460 + t399 * t487 + t400 * t484) * qJD(1) + (-(-qJD(1) * t358 - t204 - t405 - t687) * t205 + t119 * (t485 * t675 + t476 + t725) + t118 * t737 + t205 * (t459 + t727) + (t420 * t797 - t795) * qJD(2) + ((-pkin(1) - t625) * t796 + (t204 * (-rSges(3,3) - pkin(6)) + t205 * t675) * t485) * qJD(1)) * m(3) + t890 + (t1019 * t461 + (t1018 * t483 + t486 * t904) * t460 + t996) * t381 * t866 + t915 * t869 + ((t936 + t1002) * t488 + (t937 + t1003) * t485) * t707 / 0.2e1 + (((t488 * t892 + t891 + t979) * t488 + (t485 * t892 + t1010 + t980) * t485) * qJD(2) + t949 - t970) * t667 + t908 * t871 + t909 * t872 + (-(-qJD(1) * t308 - t112 + t539 + t902) * t113 + t70 * (-t308 + t728) + t112 * t692 + t69 * (t309 + t649) + t113 * (t462 + t732) + (t113 * t628 + t485 * t798) * qJD(2) + ((-t112 * rSges(4,3) + t113 * (-t456 - t848)) * t485 + (t112 * (-t392 - t456) - t113 * t482) * t488) * qJD(1)) * m(4) - t842 / 0.2e1 - t840 / 0.2e1 + t841 / 0.2e1 + t839 / 0.2e1 + t831 / 0.2e1 + t829 / 0.2e1 + t830 / 0.2e1 + (t997 + t951) * t868 + t828 / 0.2e1 + (t940 + t943) * t666 - (-t941 + t942 + t948) * t716 / 0.2e1 + (t26 * (-t623 + t728) + t65 * (t409 - t624 + t692) + t25 * (t437 + t197 + t363) + t66 * (-pkin(3) * t685 + t410 + t575 + t697) + (qJD(2) * t461 * t65 * t861 - t25 * t482 + t26 * t545) * t485 + ((t460 * t861 - t456 - t858) * t832 + (t65 * (-t456 + t545) - t66 * t482) * t488) * qJD(1) - (t513 - t65 + t916) * t66) * m(5) + ((((t1014 + t1034) * t488 + t981 + t1009 + t1030) * t488 - t891 * t485) * qJD(2) + t978) * t664; t921 * t872 + t922 * t871 + ((t460 * t908 + t779 * t945) * qJD(4) + ((qJD(4) * t944 + t924) * t461 + t928) * t488 + (t368 * t930 + t369 * t929) * t433 + (t368 * t935 + t369 * t933) * t382 + (-t368 * t934 - t369 * t932) * t381) * t870 + (qJD(1) * t897 + t485 * t957 - t488 * t958) * t869 + (qJD(1) * t896 + t485 * t955 - t488 * t956) * t868 + ((t460 * t909 + t778 * t946) * qJD(4) + ((qJD(4) * t947 + t924) * t461 + t928) * t485 + (-t366 * t930 + t367 * t929) * t433 + (-t366 * t935 + t367 * t933) * t382 + (t366 * t934 - t367 * t932) * t381) * t867 + ((qJD(4) * t898 - t987) * t461 + ((-t483 * t930 + t486 * t929 + t993) * t433 + (-t483 * t935 + t486 * t933 - t1019) * t382 + (t483 * t934 - t486 * t932 + t994) * t381 + t907 * qJD(4)) * t460) * t866 + (qJD(1) * t898 + t485 * t953 - t488 * t954) * t865 - ((t529 * t460 + t461 * t877 + t528 * t484 + t487 * t878) * qJD(2) + (t460 * t735 + t461 * t734 + t484 * t731 + t487 * t730) * qJD(1)) * qJD(1) / 0.2e1 + (t941 * t488 + t940 * t485 + (t937 * t485 + t936 * t488) * qJD(1)) * qJD(1) / 0.2e1 + ((-t717 * t918 + t919) * t485 + ((t485 * t917 + t927) * qJD(2) + t926) * t488) * t667 + ((-t716 * t917 - t919) * t488 + ((t488 * t918 + t927) * qJD(2) + t926) * t485) * t664 - t950 * t715 / 0.2e1 + t923 * t631 + (t5 * t643 + t33 * t570 + (t6 * t579 + t45 * t557 + (-t657 + t33 * (-t744 - t767)) * qJD(1) + t537) * t485 - t45 * (t461 * t711 + t576) - t33 * (t713 + t740) - (t33 * t762 - t45 * t931) * t381 - (t618 * t397 + (t33 * t654 + t487 * t618) * pkin(2)) * qJD(2) - (t460 * t531 + t461 * t495) * qJD(4) + (t382 * t931 + t762 * t433 - t461 * t710 - t648 + t736) * t44 + (t7 * t579 + (t45 * t579 - t658) * qJD(1) + t879 + (-t33 * t382 - t45 * t433) * t925 + t44 * t557) * t488) * m(6) + (t65 * t736 + t16 * t643 + t52 * t570 + (t26 * t634 + t65 * t574 + t16 * t197 + t52 * t98 + (t52 * t193 + t634 * t66) * qJD(1)) * t488 + (t25 * t634 + t66 * t574 + t16 * t193 + t52 * t100 + (t65 * t293 + t52 * (-t197 - t744)) * qJD(1)) * t485 - t65 * (t259 * t433 - t295 * t382 + t648) - t66 * (-t261 * t433 - t295 * t381 + t576) - t52 * (-t259 * t381 - t261 * t382 + t740) - ((-t193 * t65 + t197 * t66) * t460 + t502 * t461) * qJD(4) - (-t609 * t397 + (-t487 * t609 + t52 * t654) * pkin(2)) * qJD(2)) * m(5) + (-(t112 * t352 + t113 * (-t353 - t704)) * qJD(1) - (t105 * pkin(2) * t654 + (-t105 * t353 + t112 * t660) * t488 + (-t105 * t352 + t113 * t660) * t485) * qJD(2) + t70 * t628 - t112 * pkin(2) * t682 + (t172 * t716 + t173 * t717 + t696) * t544 + t105 * t695 + (-t112 * t373 + t105 * t172 + (t113 * t661 + t308 * t938) * qJD(1)) * t488 + (t69 * t661 + t113 * (-t373 - t701) + t105 * t173 + (t798 + t938 * (-t301 - t309)) * qJD(1)) * t485) * m(4) + (-(t204 * t383 - t795) * qJD(1) - (t198 * (-t383 * t485 - t384 * t488) + t589 * t625) * qJD(2) + 0.2e1 * t198 * (t228 * t488 + t229 * t485 + (t358 * t488 - t359 * t485) * qJD(1)) + t589 * t402 + (-t118 * t485 - t119 * t488 + (-t205 * t488 + t797) * qJD(1)) * t420) * m(3) + (t943 * qJD(1) + t961 + ((t979 * qJD(1) + t974 * t488) * t488 + (t965 * t485 + t980 * qJD(1) + (-t966 + t975) * t488) * t485) * t962) * t863 + (t942 * qJD(1) + t960 + ((t981 * qJD(1) + t966 * t488) * t488 + (t975 * t485 + t982 * qJD(1) + (-t965 + t974) * t488) * t485) * t962) * t862 + (t949 + t952 + t971) * t720 / 0.2e1 + (t948 + t951 + t972) * t669 - (t952 * t485 + t951 * t488) * t714 / 0.2e1; 0.2e1 * (t6 * t862 + t7 * t863) * m(6) + 0.2e1 * (-t838 / 0.2e1 + t837 / 0.2e1) * m(5) + 0.2e1 * (t69 * t862 + t70 * t863) * m(4); (t460 * t896 - t461 * t909) * t872 + (t460 * t897 - t461 * t908) * t871 + (t368 * t985 + t900 * t369 - t899 * t488) * t870 + ((qJD(2) * t897 - t915) * t461 + (-qJD(1) * t922 + t908 * qJD(2) + t485 * t958 + t488 * t957) * t460) * t869 + ((qJD(2) * t896 - t997) * t461 + (-qJD(1) * t921 + t909 * qJD(2) + t485 * t956 + t488 * t955) * t460) * t868 + (-t366 * t985 + t367 * t900 - t485 * t899) * t867 + (t973 * t461 + (-t483 * t985 + t486 * t900) * t460) * t866 + ((qJD(2) * t898 - t984) * t461 + (-qJD(1) * t923 + t907 * qJD(2) + t485 * t954 + t488 * t953) * t460) * t865 - (t830 + t831 + t841 - t842 + t828 + t829 + t839 - t840 + t890) * t461 / 0.2e1 + t960 * t781 / 0.2e1 + t961 * t780 / 0.2e1 + t950 * t668 + (t460 * t898 - t461 * t907) * t631 + (-(-t44 * t759 + t45 * t758) * t433 - (t33 * t758 + t44 * t635) * t382 - (t33 * t759 + t45 * t635) * t381 + (qJD(2) * t495 + t44 * t850 - t45 * t851 - t6 * t767 - t7 * t768) * t461 + (t531 * qJD(2) + (qJD(1) * t530 - t45 * t769 + t6 * t925 + t537) * t488 + (qJD(1) * t920 + t44 * t769 - t7 * t925 - t879) * t485) * t460) * m(6) + ((qJD(2) * t502 + t65 * t100 + t26 * t193 - t25 * t197 - t66 * t98) * t461 + (t65 * (-qJD(2) * t193 + t163 * t485) + t66 * (qJD(2) * t197 - t163 * t488) + t16 * t590 + t52 * (t100 * t488 - t193 * t720 - t197 * t719 - t485 * t98) + (qJD(1) * t609 + t837 - t838) * t293) * t460 - t65 * (-t243 * t433 - t351 * t382) - t66 * (t245 * t433 - t351 * t381) - t52 * (t243 * t381 + t245 * t382)) * m(5) + t952 * (t460 * t669 + t461 * t666) + t951 * (t461 * t664 - t691 / 0.2e1); 0.2e1 * ((t44 * t716 + t45 * t717 - t5) * t873 + (t381 * t45 + t382 * t44) * t874) * t461 + 0.2e1 * ((qJD(2) * t33 - t44 * t720 + t45 * t719 + t485 * t6 + t488 * t7) * t873 + (t33 * (t381 * t485 + t382 * t488) + (-t44 * t485 + t45 * t488) * t433) * t874) * t460;];
tauc = t1(:);
