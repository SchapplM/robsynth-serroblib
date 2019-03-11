% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR1_coriolisvecJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:57:32
% EndTime: 2019-03-09 04:58:38
% DurationCPUTime: 58.27s
% Computational Cost: add. (51066->1197), mult. (38106->1540), div. (0->0), fcn. (33590->12), ass. (0->633)
t554 = qJ(3) + qJ(4);
t545 = pkin(11) + t554;
t538 = sin(t545);
t539 = cos(t545);
t546 = sin(t554);
t547 = cos(t554);
t1049 = Icges(5,5) * t546 + Icges(6,5) * t538 + Icges(5,6) * t547 + Icges(6,6) * t539;
t553 = qJ(1) + pkin(10);
t543 = sin(t553);
t862 = t539 * t543;
t865 = t538 * t543;
t544 = cos(t553);
t889 = Icges(6,3) * t544;
t308 = Icges(6,5) * t862 - Icges(6,6) * t865 - t889;
t858 = t543 * t547;
t859 = t543 * t546;
t890 = Icges(5,3) * t544;
t325 = Icges(5,5) * t858 - Icges(5,6) * t859 - t890;
t1048 = -t308 - t325;
t425 = Icges(6,5) * t539 - Icges(6,6) * t538;
t309 = Icges(6,3) * t543 + t425 * t544;
t446 = Icges(5,5) * t547 - Icges(5,6) * t546;
t326 = Icges(5,3) * t543 + t446 * t544;
t1043 = t309 + t326;
t906 = Icges(6,4) * t538;
t426 = Icges(6,2) * t539 + t906;
t512 = Icges(6,4) * t539;
t428 = Icges(6,1) * t538 + t512;
t907 = Icges(5,4) * t546;
t447 = Icges(5,2) * t547 + t907;
t537 = Icges(5,4) * t547;
t449 = Icges(5,1) * t546 + t537;
t1047 = t426 * t538 - t428 * t539 + t447 * t546 - t449 * t547;
t429 = Icges(6,1) * t539 - t906;
t313 = Icges(6,5) * t543 + t429 * t544;
t450 = Icges(5,1) * t547 - t907;
t330 = Icges(5,5) * t543 + t450 * t544;
t1046 = -t313 * t862 - t330 * t858;
t1039 = t425 + t446;
t1035 = t1049 * t544;
t1034 = t1049 * t543;
t1045 = t1047 * t543 + t1035;
t1044 = -t1047 * t544 + t1034;
t669 = -Icges(5,2) * t546 + t537;
t328 = Icges(5,6) * t543 + t544 * t669;
t1042 = t309 * t544 + t328 * t859 + t1046;
t851 = t544 * t547;
t861 = t539 * t544;
t1041 = -t1043 * t543 - t313 * t861 - t330 * t851;
t456 = Icges(6,4) * t865;
t900 = Icges(6,5) * t544;
t312 = Icges(6,1) * t862 - t456 - t900;
t473 = Icges(5,4) * t859;
t901 = Icges(5,5) * t544;
t329 = Icges(5,1) * t858 - t473 - t901;
t1040 = t1048 * t543 - t312 * t861 - t329 * t851;
t894 = Icges(5,6) * t544;
t327 = Icges(5,4) * t858 - Icges(5,2) * t859 - t894;
t878 = t327 * t546;
t656 = -t329 * t547 + t878;
t624 = t656 * t543;
t893 = Icges(6,6) * t544;
t310 = Icges(6,4) * t862 - Icges(6,2) * t865 - t893;
t658 = t310 * t538 - t312 * t539;
t625 = t658 * t543;
t880 = t308 * t544;
t1007 = -t325 * t544 - t624 - t625 - t880;
t668 = -Icges(6,2) * t538 + t512;
t311 = Icges(6,6) * t543 + t544 * t668;
t1006 = -t311 * t865 - t326 * t544 - t1042;
t852 = t544 * t546;
t864 = t538 * t544;
t1005 = -t310 * t864 - t327 * t852 - t1040;
t1004 = -t311 * t864 - t328 * t852 - t1041;
t552 = qJD(3) + qJD(4);
t992 = -t447 + t450;
t709 = t992 * t552;
t991 = t449 + t669;
t710 = t991 * t552;
t994 = -t426 + t429;
t711 = t994 * t552;
t993 = t428 + t668;
t712 = t993 * t552;
t1038 = qJD(1) * t1049 - t538 * t712 + t539 * t711 - t546 * t710 + t547 * t709;
t775 = rSges(5,1) * t858;
t561 = -pkin(8) - pkin(7);
t511 = t544 * t561;
t559 = cos(qJ(3));
t540 = t559 * pkin(3) + pkin(2);
t810 = t543 * t540 + t511;
t557 = sin(qJ(1));
t938 = pkin(1) * t557;
t1037 = -t775 - t938 - t810;
t1036 = t656 + t658;
t1033 = qJD(1) * t1047 + t1039 * t552;
t1032 = t1044 * qJD(1);
t1031 = t1045 * qJD(1);
t1030 = t1043 * qJD(1);
t555 = sin(qJ(6));
t850 = t544 * t555;
t558 = cos(qJ(6));
t854 = t543 * t558;
t390 = t539 * t854 - t850;
t372 = Icges(7,4) * t390;
t848 = t544 * t558;
t856 = t543 * t555;
t389 = t539 * t856 + t848;
t211 = -Icges(7,2) * t389 + Icges(7,6) * t865 + t372;
t371 = Icges(7,4) * t389;
t215 = -Icges(7,1) * t390 - Icges(7,5) * t865 + t371;
t1015 = t211 * t555 + t215 * t558;
t208 = Icges(7,5) * t390 - Icges(7,6) * t389 + Icges(7,3) * t865;
t97 = -t1015 * t538 - t208 * t539;
t517 = qJD(3) * t543;
t439 = qJD(4) * t543 + t517;
t440 = t544 * t552;
t1029 = t1006 * t439 - t1007 * t440 - t1031;
t1028 = t1004 * t439 - t1005 * t440 + t1032;
t998 = -t449 * t543 - t327;
t721 = qJD(1) * t330 + t552 * t998;
t857 = t543 * t552;
t723 = qJD(1) * t328 + t329 * t552 - t447 * t857;
t1001 = -t428 * t543 - t310;
t725 = qJD(1) * t313 + t1001 * t552;
t727 = qJD(1) * t311 + t312 * t552 - t426 * t857;
t1027 = -t538 * t725 - t539 * t727 - t546 * t721 - t547 * t723;
t997 = -t449 * t544 - t328;
t722 = (-t450 * t543 + t901) * qJD(1) + t997 * t552;
t996 = -t447 * t544 + t330;
t724 = (-t543 * t669 + t894) * qJD(1) + t996 * t552;
t1000 = -t428 * t544 - t311;
t726 = (-t429 * t543 + t900) * qJD(1) + t1000 * t552;
t999 = -t426 * t544 + t313;
t728 = (-t543 * t668 + t893) * qJD(1) + t999 * t552;
t1026 = t538 * t726 + t539 * t728 + t546 * t722 + t547 * t724;
t1025 = t1033 * t543 + t1038 * t544;
t1024 = -t1033 * t544 + t1038 * t543;
t1023 = t310 * t539 + t312 * t538 + t327 * t547 + t329 * t546;
t1022 = t311 * t539 + t313 * t538 + t328 * t547 + t330 * t546;
t1021 = -t538 * t728 + t539 * t726 - t546 * t724 + t547 * t722 + t1030;
t1020 = qJD(1) * t1048 + t538 * t727 - t539 * t725 + t546 * t723 - t547 * t721;
t877 = t328 * t546;
t715 = -t325 + t877;
t879 = t311 * t538;
t1019 = t715 - t308 + t879;
t1018 = qJD(1) * t1036 - t1034 * t552 + t1030;
t1017 = -t1035 * t552 + (-t1039 * t543 - t313 * t539 - t330 * t547 + t877 + t879 + t889 + t890) * qJD(1);
t1016 = -t624 + t1041;
t935 = pkin(4) * t547;
t469 = t540 + t935;
t551 = -qJ(5) + t561;
t817 = -t543 * t469 - t544 * t551;
t665 = Icges(7,5) * t558 - Icges(7,6) * t555;
t331 = -Icges(7,3) * t539 + t538 * t665;
t903 = Icges(7,4) * t558;
t667 = -Icges(7,2) * t555 + t903;
t333 = -Icges(7,6) * t539 + t538 * t667;
t904 = Icges(7,4) * t555;
t672 = Icges(7,1) * t558 - t904;
t335 = -Icges(7,5) * t539 + t538 * t672;
t124 = t331 * t865 - t333 * t389 + t335 * t390;
t785 = qJD(6) * t544;
t345 = t538 * t785 + t439;
t786 = qJD(6) * t543;
t346 = -t538 * t786 + t440;
t787 = qJD(6) * t539;
t483 = qJD(1) - t787;
t87 = t208 * t865 - t211 * t389 - t215 * t390;
t391 = -t539 * t850 + t854;
t392 = t539 * t848 + t856;
t210 = Icges(7,5) * t392 + Icges(7,6) * t391 + Icges(7,3) * t864;
t905 = Icges(7,4) * t392;
t213 = Icges(7,2) * t391 + Icges(7,6) * t864 + t905;
t373 = Icges(7,4) * t391;
t216 = Icges(7,1) * t392 + Icges(7,5) * t864 + t373;
t88 = t210 * t865 - t389 * t213 + t390 * t216;
t33 = t124 * t483 + t345 * t88 - t346 * t87;
t125 = t331 * t864 + t333 * t391 + t335 * t392;
t89 = t208 * t864 + t391 * t211 - t215 * t392;
t90 = t210 * t864 + t391 * t213 + t392 * t216;
t34 = t125 * t483 + t345 * t90 - t346 * t89;
t563 = qJD(1) ^ 2;
t1014 = -t1018 * t543 + t1020 * t544;
t1013 = t1017 * t543 + t1021 * t544;
t1012 = t1018 * t544 + t1020 * t543;
t1011 = -t1017 * t544 + t1021 * t543;
t685 = rSges(7,1) * t390 - rSges(7,2) * t389;
t218 = rSges(7,3) * t865 + t685;
t933 = pkin(5) * t539;
t435 = pkin(9) * t538 + t933;
t366 = t435 * t543;
t1003 = t218 + t366;
t220 = t392 * rSges(7,1) + t391 * rSges(7,2) + rSges(7,3) * t864;
t368 = pkin(5) * t861 + pkin(9) * t864;
t1002 = t220 + t368;
t684 = rSges(7,1) * t558 - rSges(7,2) * t555;
t339 = -rSges(7,3) * t539 + t538 * t684;
t934 = pkin(5) * t538;
t434 = -pkin(9) * t539 + t934;
t995 = t339 + t434;
t583 = qJD(1) * t992 + t439 * t997 - t440 * t998;
t584 = qJD(1) * t994 + t1000 * t439 - t1001 * t440;
t957 = qJD(1) * t993 + t439 * t999 - t440 * (-Icges(6,2) * t862 + t312 - t456);
t958 = qJD(1) * t991 + t439 * t996 - t440 * (-Icges(5,2) * t858 + t329 - t473);
t990 = -t538 * t957 + t584 * t539 - t546 * t958 + t583 * t547;
t989 = qJD(1) * t1039 + t1034 * t440 - t1035 * t439;
t795 = qJD(1) * t543;
t759 = t538 * t795;
t768 = t539 * t440;
t987 = t759 - t768;
t986 = t33 * t543 + t34 * t544;
t985 = -t218 * t483 - t339 * t346;
t984 = 0.2e1 * qJD(3);
t946 = t439 / 0.2e1;
t945 = -t440 / 0.2e1;
t527 = t543 * rSges(6,3);
t317 = rSges(6,1) * t861 - rSges(6,2) * t864 + t527;
t431 = rSges(6,1) * t538 + rSges(6,2) * t539;
t423 = t544 * t469;
t475 = t544 * t540;
t802 = -t551 + t561;
t264 = t543 * t802 + t423 - t475;
t534 = t543 * pkin(7);
t444 = t544 * pkin(2) + t534;
t708 = -t543 * t561 + t475;
t315 = t708 - t444;
t560 = cos(qJ(1));
t550 = t560 * pkin(1);
t736 = t444 + t550;
t697 = t315 + t736;
t645 = t264 + t697;
t936 = pkin(4) * t546;
t404 = t439 * t936;
t556 = sin(qJ(3));
t790 = qJD(3) * t556;
t771 = pkin(3) * t790;
t481 = t543 * t771;
t515 = qJD(5) * t544;
t760 = -t404 - t481 - t515;
t106 = -t431 * t439 + (t317 + t645) * qJD(1) + t760;
t361 = t431 * t543;
t928 = rSges(6,1) * t539;
t432 = -rSges(6,2) * t538 + t928;
t181 = -t552 * t361 + (t432 * t544 + t527) * qJD(1);
t381 = t432 * t552;
t781 = qJD(1) * qJD(3);
t505 = t543 * t781;
t780 = qJD(1) * qJD(4);
t420 = t543 * t780 + t505;
t406 = t440 * t935;
t776 = t563 * t550;
t844 = t559 * qJD(3) ^ 2;
t601 = -pkin(3) * t544 * t844 + qJD(1) * t481 - t776;
t779 = qJD(1) * qJD(5);
t588 = -t552 * t406 + t420 * t936 + t544 * t779 + t601;
t846 = t546 * t552;
t437 = -pkin(4) * t846 - t771;
t794 = qJD(1) * t544;
t806 = t561 * t795 + t481;
t807 = t551 * t795 + t515;
t811 = t469 - t540;
t149 = t437 * t543 + t794 * t811 + t806 - t807;
t930 = pkin(2) - t540;
t259 = (-t544 * t930 - t534) * qJD(1) - t806;
t430 = t444 * qJD(1);
t840 = -t259 - t430;
t766 = -t149 + t840;
t62 = -t381 * t440 + t420 * t431 + (-t181 + t766) * qJD(1) + t588;
t983 = qJD(1) * t106 + t62;
t633 = t684 * t539;
t683 = -rSges(7,1) * t555 - rSges(7,2) * t558;
t217 = t552 * t633 + (rSges(7,3) * t552 + qJD(6) * t683) * t538;
t395 = t435 * t552;
t981 = -t217 - t395;
t535 = t544 * pkin(7);
t443 = pkin(2) * t543 - t535;
t314 = t443 - t810;
t438 = qJD(1) * t443;
t980 = qJD(1) * t314 - t438;
t529 = t543 * rSges(4,3);
t847 = t544 * t559;
t849 = t544 * t556;
t364 = rSges(4,1) * t847 - rSges(4,2) * t849 + t529;
t979 = t364 + t736;
t978 = -rSges(5,2) * t859 - t544 * rSges(5,3);
t730 = t544 * rSges(3,1) - rSges(3,2) * t543;
t548 = Icges(4,4) * t559;
t670 = -Icges(4,2) * t556 + t548;
t488 = Icges(4,1) * t556 + t548;
t977 = t550 + t730;
t274 = t339 * t543;
t925 = rSges(7,3) * t538;
t340 = t633 + t925;
t365 = t434 * t543;
t793 = qJD(1) * t546;
t758 = t543 * t793;
t466 = pkin(4) * t758;
t750 = t539 * t786;
t788 = qJD(6) * t538;
t976 = -qJD(1) * t365 + t218 * t788 - t274 * t483 - t339 * t750 + t346 * t340 + t440 * t435 + t795 * t995 + t406 + t466;
t845 = t547 * t552;
t772 = rSges(5,2) * t845;
t975 = rSges(5,1) * t846 + t772;
t337 = t775 + t978;
t528 = t543 * rSges(5,3);
t809 = rSges(5,1) * t851 + t528;
t338 = -rSges(5,2) * t852 + t809;
t791 = qJD(3) * t544;
t746 = -t314 * t517 + t315 * t791 + qJD(2);
t120 = t337 * t439 + t338 * t440 + t746;
t451 = rSges(5,1) * t546 + rSges(5,2) * t547;
t752 = t544 * t790;
t706 = pkin(3) * t752;
t617 = -t440 * t451 - t706;
t737 = -t443 - t938;
t698 = t314 + t737;
t131 = (-t337 + t698) * qJD(1) + t617;
t132 = -t439 * t451 - t481 + (t338 + t697) * qJD(1);
t393 = t451 * t543;
t394 = t451 * t544;
t926 = rSges(5,2) * t546;
t452 = rSges(5,1) * t547 - t926;
t972 = -t131 * (qJD(1) * t393 - t440 * t452) - t120 * (-t439 * t393 - t394 * t440) - t132 * (-qJD(1) * t394 - t439 * t452);
t853 = t543 * t559;
t855 = t543 * t556;
t891 = Icges(4,3) * t544;
t353 = Icges(4,5) * t853 - Icges(4,6) * t855 - t891;
t497 = Icges(4,4) * t855;
t902 = Icges(4,5) * t544;
t357 = Icges(4,1) * t853 - t497 - t902;
t895 = Icges(4,6) * t544;
t355 = Icges(4,4) * t853 - Icges(4,2) * t855 - t895;
t874 = t355 * t556;
t653 = -t357 * t559 + t874;
t141 = -t353 * t544 - t543 * t653;
t485 = Icges(4,5) * t559 - Icges(4,6) * t556;
t484 = Icges(4,5) * t556 + Icges(4,6) * t559;
t620 = qJD(3) * t484;
t908 = Icges(4,4) * t556;
t489 = Icges(4,1) * t559 - t908;
t358 = Icges(4,5) * t543 + t489 * t544;
t356 = Icges(4,6) * t543 + t544 * t670;
t873 = t356 * t556;
t652 = -t358 * t559 + t873;
t967 = -t544 * t620 + (-t485 * t543 + t652 + t891) * qJD(1);
t354 = Icges(4,3) * t543 + t485 * t544;
t798 = qJD(1) * t354;
t966 = qJD(1) * t653 - t543 * t620 + t798;
t486 = Icges(4,2) * t559 + t908;
t648 = t486 * t556 - t488 * t559;
t963 = t648 * qJD(1) + t485 * qJD(3);
t821 = -Icges(4,2) * t853 + t357 - t497;
t823 = t488 * t543 + t355;
t962 = -t556 * t821 - t559 * t823;
t863 = t538 * t552;
t596 = t483 * t558 + t555 * t863;
t705 = qJD(1) * t539 - qJD(6);
t183 = t544 * t596 + t705 * t856;
t595 = t483 * t555 - t558 * t863;
t184 = t544 * t595 - t705 * t854;
t765 = t184 * rSges(7,1) + t183 * rSges(7,2) + rSges(7,3) * t768;
t118 = -rSges(7,3) * t759 + t765;
t436 = pkin(9) * t768;
t769 = t538 * t440;
t614 = -t539 * t795 - t769;
t225 = pkin(5) * t614 - pkin(9) * t759 + t436;
t506 = t544 * t781;
t421 = t544 * t780 + t506;
t262 = -qJD(6) * t987 + t421;
t514 = qJD(5) * t543;
t818 = t544 * t437 + t514;
t148 = t706 + (-t543 * t811 + t544 * t802) * qJD(1) + t818;
t513 = pkin(7) * t794;
t258 = -t706 - t513 + (t543 * t930 - t511) * qJD(1);
t777 = t563 * t938;
t701 = qJD(1) * (-pkin(2) * t795 + t513) - t777;
t569 = qJD(1) * t258 + (-t506 * t556 - t543 * t844) * pkin(3) + t701;
t568 = qJD(1) * t148 + t543 * t779 + (-t421 * t546 - t439 * t845) * pkin(4) + t569;
t751 = t552 * t788;
t26 = qJD(1) * t225 + t483 * t118 - t345 * t217 + t220 * t751 - t262 * t339 - t439 * t395 - t421 * t434 + t568;
t613 = t538 * t794 + t539 * t857;
t185 = t543 * t596 - t705 * t850;
t186 = t543 * t595 + t705 * t848;
t686 = rSges(7,1) * t186 + rSges(7,2) * t185;
t119 = rSges(7,3) * t613 + t686;
t226 = t613 * pkin(9) + (-t538 * t857 + t539 * t794) * pkin(5);
t261 = qJD(6) * t613 + t420;
t27 = -t218 * t751 - t119 * t483 - t217 * t346 + t261 * t339 - t395 * t440 + t420 * t434 + (-t226 + t766) * qJD(1) + t588;
t77 = t220 * t483 - t339 * t345 - t434 * t439 + (t368 + t645) * qJD(1) + t760;
t961 = (qJD(1) * t77 + t27) * t544 + t26 * t543;
t626 = t665 * t539;
t654 = -t333 * t555 + t335 * t558;
t660 = -t213 * t555 + t216 * t558;
t960 = t345 * (-t331 * t544 - t660) - t346 * (-t331 * t543 + t1015) + t483 * (Icges(7,3) * t538 + t626 - t654);
t666 = -Icges(7,2) * t558 - t904;
t959 = t345 * (-Icges(7,2) * t392 + t216 + t373) - t346 * (-Icges(7,2) * t390 - t215 - t371) + t483 * (t666 * t538 + t335);
t956 = t261 / 0.2e1;
t955 = t262 / 0.2e1;
t954 = t543 * t945 + t544 * t946;
t953 = -t345 / 0.2e1;
t952 = t345 / 0.2e1;
t951 = -t346 / 0.2e1;
t950 = t346 / 0.2e1;
t949 = t420 / 0.2e1;
t948 = t421 / 0.2e1;
t947 = -t439 / 0.2e1;
t944 = t440 / 0.2e1;
t943 = -t483 / 0.2e1;
t942 = t483 / 0.2e1;
t941 = t543 / 0.2e1;
t940 = -t544 / 0.2e1;
t939 = -rSges(7,3) - pkin(9);
t937 = pkin(3) * t556;
t932 = -qJD(1) / 0.2e1;
t931 = qJD(1) / 0.2e1;
t929 = rSges(4,1) * t559;
t113 = Icges(7,5) * t186 + Icges(7,6) * t185 + Icges(7,3) * t613;
t115 = Icges(7,4) * t186 + Icges(7,2) * t185 + Icges(7,6) * t613;
t117 = Icges(7,1) * t186 + Icges(7,4) * t185 + Icges(7,5) * t613;
t24 = (-t1015 * t552 - t113) * t539 + (-t115 * t555 + t117 * t558 + t208 * t552 + (-t211 * t558 + t215 * t555) * qJD(6)) * t538;
t923 = t24 * t346;
t112 = Icges(7,5) * t184 + Icges(7,6) * t183 - Icges(7,3) * t987;
t114 = Icges(7,4) * t184 + Icges(7,2) * t183 - Icges(7,6) * t987;
t116 = Icges(7,1) * t184 + Icges(7,4) * t183 - Icges(7,5) * t987;
t25 = (t552 * t660 - t112) * t539 + (-t114 * t555 + t116 * t558 + t210 * t552 + (-t213 * t558 - t216 * t555) * qJD(6)) * t538;
t922 = t25 * t345;
t921 = t26 * t544;
t920 = t27 * t543;
t915 = t543 * t77;
t914 = t97 * t261;
t98 = -t210 * t539 + t538 * t660;
t913 = t98 * t262;
t137 = -t331 * t539 + t538 * t654;
t664 = -Icges(7,5) * t555 - Icges(7,6) * t558;
t199 = t552 * t626 + (Icges(7,3) * t552 + qJD(6) * t664) * t538;
t627 = t667 * t539;
t200 = t552 * t627 + (Icges(7,6) * t552 + qJD(6) * t666) * t538;
t628 = t672 * t539;
t671 = -Icges(7,1) * t555 - t903;
t201 = t552 * t628 + (Icges(7,5) * t552 + qJD(6) * t671) * t538;
t55 = (t552 * t654 - t199) * t539 + (-t200 * t555 + t201 * t558 + t331 * t552 + (-t333 * t558 - t335 * t555) * qJD(6)) * t538;
t912 = t137 * t751 + t55 * t483;
t263 = t810 + t817;
t703 = -t439 * t263 + t746;
t838 = t264 + t368;
t47 = t218 * t345 + t220 * t346 + t366 * t439 + t440 * t838 + t703;
t887 = qJD(1) * t47;
t773 = rSges(6,1) * t862;
t813 = -rSges(6,2) * t865 - t544 * rSges(6,3);
t316 = t773 + t813;
t839 = t264 + t317;
t80 = t316 * t439 + t440 * t839 + t703;
t885 = qJD(1) * t80;
t803 = rSges(4,2) * t855 + t544 * rSges(4,3);
t363 = rSges(4,1) * t853 - t803;
t490 = rSges(4,1) * t556 + rSges(4,2) * t559;
t753 = t490 * t791;
t194 = -t753 + (-t363 + t737) * qJD(1);
t883 = t194 * t543;
t882 = t194 * t544;
t755 = t490 * t517;
t195 = qJD(1) * t979 - t755;
t418 = t490 * t544;
t881 = t195 * t418;
t870 = t439 * t547;
t867 = t484 * t543;
t866 = t484 * t544;
t860 = t539 * t552;
t841 = -t543 * t263 + t544 * t264;
t832 = -t543 * t314 + t544 * t315;
t831 = t543 * t337 + t544 * t338;
t830 = -t543 * t353 - t357 * t847;
t829 = t543 * t354 + t358 * t847;
t822 = -t488 * t544 - t356;
t820 = -t486 * t544 + t358;
t819 = t431 * t795 + t466;
t816 = t423 + t550;
t812 = rSges(5,2) * t758 + rSges(5,3) * t794;
t792 = qJD(1) * t556;
t808 = rSges(4,2) * t543 * t792 + rSges(4,3) * t794;
t805 = -t486 + t489;
t804 = t488 + t670;
t478 = -t936 - t937;
t734 = t478 + t937;
t359 = t734 * t543;
t797 = qJD(1) * t359;
t796 = qJD(1) * t485;
t789 = qJD(3) * t559;
t245 = -t543 * t648 - t866;
t782 = t245 * qJD(1);
t778 = pkin(4) * t845;
t405 = t440 * t936;
t770 = pkin(3) * t789;
t767 = t544 * t148 + t543 * t149 - t263 * t794;
t221 = -t544 * t772 + (-t544 * t846 - t547 * t795) * rSges(5,1) + t812;
t222 = -t552 * t393 + (t452 * t544 + t528) * qJD(1);
t764 = t544 * t221 + t543 * t222 + t337 * t794;
t763 = -t220 - t838;
t762 = t544 * t258 + t543 * t259 - t314 * t794;
t756 = t544 * t792;
t749 = t539 * t785;
t748 = t863 / 0.2e1;
t745 = -pkin(2) - t929;
t744 = t795 / 0.2e1;
t743 = t794 / 0.2e1;
t742 = -t517 / 0.2e1;
t739 = t791 / 0.2e1;
t735 = -t451 - t937;
t733 = -t431 - t936;
t732 = -t434 - t936;
t729 = t556 * (-t543 ^ 2 - t544 ^ 2);
t360 = t734 * t544;
t719 = t439 * t359 + t360 * t440;
t362 = t431 * t544;
t718 = -t439 * t361 - t362 * t440;
t305 = t358 * t853;
t717 = t354 * t544 - t305;
t714 = -qJD(1) * t362 - t432 * t439;
t713 = -t353 + t873;
t704 = t543 * t316 + t544 * t317 + t841;
t702 = qJD(6) * t748;
t700 = -pkin(4) * t870 + qJD(1) * t360;
t699 = -t381 - t778;
t696 = t817 - t938;
t410 = t452 * t552;
t692 = -t410 - t770;
t690 = qJD(1) * t361 - t440 * t432 - t406;
t442 = rSges(3,1) * t543 + rSges(3,2) * t544;
t687 = -rSges(4,2) * t556 + t929;
t644 = t514 - t706;
t631 = -t405 + t644;
t646 = t263 + t698;
t76 = -t434 * t440 + (-t366 + t646) * qJD(1) + t631 + t985;
t682 = t544 * t76 + t915;
t681 = t543 * t88 - t544 * t87;
t680 = t543 * t87 + t544 * t88;
t679 = t543 * t90 - t544 * t89;
t678 = t543 * t89 + t544 * t90;
t677 = t543 * t98 - t544 * t97;
t676 = t543 * t97 + t544 * t98;
t663 = -t131 * t544 - t132 * t543;
t662 = -t195 * t543 - t882;
t659 = t218 * t544 - t220 * t543;
t227 = t355 * t559 + t357 * t556;
t228 = t356 * t559 + t358 * t556;
t651 = t363 * t543 + t364 * t544;
t647 = -t778 + t981;
t643 = -t431 + t478;
t632 = rSges(6,2) * t987 + rSges(6,3) * t794;
t180 = rSges(6,1) * t614 + t632;
t640 = t544 * t180 + t543 * t181 + t316 * t794 + t767;
t639 = t1002 * t544 + t1003 * t543 + t841;
t417 = t490 * t543;
t623 = -t770 - t778;
t622 = qJD(3) * t488;
t621 = qJD(3) * t486;
t142 = -t356 * t855 - t717;
t619 = (-t141 * t544 + t142 * t543) * qJD(3);
t143 = -t355 * t849 - t830;
t144 = -t356 * t849 + t829;
t618 = (-t143 * t544 + t144 * t543) * qJD(3);
t616 = -t435 - t925;
t611 = (-t544 * t793 - t870) * pkin(4);
t610 = -t208 * t346 + t210 * t345 + t331 * t483;
t609 = (-Icges(7,5) * t389 - Icges(7,6) * t390) * t346 - (Icges(7,5) * t391 - Icges(7,6) * t392) * t345 - t664 * t538 * t483;
t608 = -t381 + t623;
t605 = t258 * t791 + t259 * t517 - t314 * t506 - t315 * t505;
t275 = t339 * t544;
t367 = t434 * t544;
t604 = t218 * t749 - t345 * t274 - t275 * t346 - t439 * t365 - t367 * t440;
t603 = -t556 * t820 + t559 * t822;
t602 = qJD(1) * t263 + t644 + t980;
t600 = t623 + t981;
t599 = t767 + t1003 * t794 + (t118 + t225) * t544 + (t119 + t226) * t543;
t598 = t538 * t609;
t591 = t439 * t149 - t421 * t263 + t605;
t590 = (-t556 * t804 + t559 * t805) * qJD(1);
t587 = (Icges(7,1) * t391 - t213 - t905) * t345 - (-Icges(7,1) * t389 - t211 - t372) * t346 + (t671 * t538 - t333) * t483;
t247 = -rSges(4,2) * t544 * t789 + (-t559 * t795 - t752) * rSges(4,1) + t808;
t248 = -qJD(3) * t417 + (t544 * t687 + t529) * qJD(1);
t585 = t247 * t544 + t248 * t543 + (t363 * t544 - t364 * t543) * qJD(1);
t582 = -qJD(1) * t367 + t220 * t788 - t483 * t275 - t339 * t749 - t340 * t345 - t435 * t439;
t242 = qJD(1) * t356 - t543 * t621;
t244 = qJD(1) * t358 - t543 * t622;
t572 = qJD(1) * t353 - qJD(3) * t227 - t242 * t556 + t244 * t559;
t241 = -t544 * t621 + (-t543 * t670 + t895) * qJD(1);
t243 = -t544 * t622 + (-t489 * t543 + t902) * qJD(1);
t571 = -qJD(3) * t228 - t241 * t556 + t243 * t559 + t798;
t461 = t670 * qJD(3);
t462 = t489 * qJD(3);
t570 = qJD(1) * t484 - t461 * t556 + t462 * t559 + (-t486 * t559 - t488 * t556) * qJD(3);
t567 = t960 * t538;
t20 = t113 * t864 + t115 * t391 + t117 * t392 + t183 * t211 - t184 * t215 - t208 * t987;
t21 = t112 * t864 + t114 * t391 + t116 * t392 + t183 * t213 + t184 * t216 - t210 * t987;
t22 = t113 * t865 - t115 * t389 + t117 * t390 + t185 * t211 - t186 * t215 + t208 * t613;
t23 = t112 * t865 - t114 * t389 + t116 * t390 + t185 * t213 + t186 * t216 + t210 * t613;
t270 = t333 * t543;
t271 = t333 * t544;
t272 = t335 * t543;
t273 = t335 * t544;
t43 = t183 * t333 + t184 * t335 + t199 * t864 + t200 * t391 + t201 * t392 - t331 * t987;
t3 = t125 * t751 - t20 * t346 + t21 * t345 + t261 * t89 + t262 * t90 + t43 * t483;
t334 = Icges(7,6) * t538 + t627;
t336 = Icges(7,5) * t538 + t628;
t37 = t137 * t483 + t345 * t98 - t346 * t97;
t44 = t185 * t333 + t186 * t335 + t199 * t865 - t200 * t389 + t201 * t390 + t331 * t613;
t4 = t124 * t751 - t22 * t346 + t23 * t345 + t261 * t87 + t262 * t88 + t44 * t483;
t564 = (((t271 * t555 - t273 * t558 + t210) * t345 - (t270 * t555 - t272 * t558 + t208) * t346 + (-t334 * t555 + t336 * t558 + t331) * t483 + t137 * qJD(6)) * t538 + (qJD(6) * t676 - t960) * t539) * t943 + t677 * t702 + (qJD(1) * t676 - t24 * t544 + t25 * t543) * t942 + ((t271 * t389 - t273 * t390) * t345 - (t270 * t389 - t272 * t390) * t346 + (-t334 * t389 + t336 * t390) * t483 + (t124 * t538 + t861 * t88) * qJD(6) + ((qJD(6) * t87 + t610) * t539 + t567) * t543) * t950 + (qJD(1) * t680 - t22 * t544 + t23 * t543) * t951 + (qJD(1) * t678 - t20 * t544 + t21 * t543) * t952 + ((-t271 * t391 - t273 * t392) * t345 - (-t270 * t391 - t272 * t392) * t346 + (t334 * t391 + t336 * t392) * t483 + (t125 * t538 + t862 * t89) * qJD(6) + ((qJD(6) * t90 + t610) * t539 + t567) * t544) * t953 + t679 * t955 + t681 * t956 - t37 * t788 / 0.2e1 + (t1006 * t543 - t1007 * t544) * t949 + (t1004 * t543 - t1005 * t544) * t948 + (t543 * t989 + t544 * t990) * t947 + (t1014 * t544 + t1013 * t543 + (t1004 * t544 + t1005 * t543) * qJD(1)) * t946 + (t1012 * t544 + t1011 * t543 + (t1006 * t544 + t1007 * t543) * qJD(1)) * t945 + (t543 * t990 - t544 * t989) * t944 + (t584 * t538 + t539 * t957 + t583 * t546 + t547 * t958) * t932 + (t1027 * t544 + t1026 * t543 + (t1022 * t544 + t1023 * t543) * qJD(1)) * t931 - t986 * t787 / 0.2e1 + (qJD(1) * t1025 + t1004 * t421 + t1005 * t420 + t1013 * t439 + t1014 * t440 + t3) * t941 + (qJD(1) * t1024 + t1006 * t421 + t1007 * t420 + t1011 * t439 + t1012 * t440 + t4) * t940 + (t33 + t1029) * t744 + (t34 + t1028) * t743;
t467 = t687 * qJD(3);
t402 = t683 * t538;
t369 = t544 * t405;
t257 = rSges(7,1) * t391 - rSges(7,2) * t392;
t256 = -rSges(7,1) * t389 - rSges(7,2) * t390;
t246 = -t544 * t648 + t867;
t235 = t246 * qJD(1);
t191 = qJD(3) * t651 + qJD(2);
t139 = -t776 - t467 * t791 + (-t248 - t430 + t755) * qJD(1);
t138 = -t467 * t517 + (t247 - t753) * qJD(1) + t701;
t123 = t570 * t543 - t544 * t963;
t122 = t543 * t963 + t570 * t544;
t111 = -qJD(3) * t652 + t241 * t559 + t243 * t556;
t110 = -qJD(3) * t653 + t242 * t559 + t244 * t556;
t107 = t585 * qJD(3);
t105 = -t431 * t440 + (-t316 + t646) * qJD(1) + t631;
t96 = -t410 * t440 + t420 * t451 + (-t222 + t840) * qJD(1) + t601;
t95 = qJD(1) * t221 - t410 * t439 - t421 * t451 + t569;
t84 = t235 + t618;
t83 = t619 + t782;
t61 = qJD(1) * t180 - t381 * t439 - t421 * t431 + t568;
t53 = t221 * t440 + t222 * t439 + t337 * t421 - t338 * t420 + t605;
t28 = t181 * t439 + t316 * t421 + (t148 + t180) * t440 - t839 * t420 + t591;
t13 = t118 * t346 + t119 * t345 + t218 * t262 - t220 * t261 + t226 * t439 + t366 * t421 + (t148 + t225) * t440 - t838 * t420 + t591;
t1 = [((t228 + t246) * t544 + (t227 + t245) * t543) * t781 / 0.2e1 + (t235 + ((t142 - t305 + (t354 + t874) * t544 + t830) * t544 + t829 * t543) * qJD(3)) * t739 + (t83 - t782 + ((t544 * t713 + t144 - t829) * t544 + (t543 * t713 + t143 + t717) * t543) * qJD(3)) * t742 + m(3) * ((-t442 * t563 - t777) * t977 + (-t776 + (-0.2e1 * t730 - t550 + t977) * t563) * (-t442 - t938)) + (t27 * (-t685 + t696) + t76 * (-t686 + t807) + t26 * (t816 + t1002) + t77 * (-pkin(5) * t769 + t436 + t765 + t818) + (-t26 * t551 + t27 * t616 + (-t437 + (t539 * t939 + t934) * t552) * t76) * t543 + ((-t557 * t77 - t560 * t76) * pkin(1) + (t76 * (-t469 + t616) - t77 * t551) * t544 + (t538 * t939 - t469 - t933) * t915) * qJD(1) - (-t76 + t732 * t440 + (-t366 - t938) * qJD(1) + t602 + t985) * t77) * m(7) + (t1022 + t1044) * t948 + (t1023 - t1045) * t949 + (t139 * (t543 * t745 + t535 + t803 - t938) + t138 * t979 + t195 * (t513 + t808) + (t490 * t883 - t881) * qJD(3) + ((-t194 * t560 - t195 * t557) * pkin(1) + (-pkin(2) - t687) * t882 + (t194 * (-rSges(4,3) - pkin(7)) + t195 * t745) * t543) * qJD(1) - (-t753 - t194 - t438 + (-t363 - t938) * qJD(1)) * t195) * m(4) + (t62 * (t696 - t813) + t61 * (t317 + t816) + (-t61 * t551 - t62 * t928) * t543 + (t807 + (rSges(6,1) * t863 + rSges(6,2) * t860 - t437) * t543 + (-t527 - t550 + (-t432 - t469) * t544) * qJD(1)) * t105 + (-rSges(6,1) * t769 - t733 * t440 + t105 - t602 + t632 + t818 + (t316 - t773 + t817) * qJD(1)) * t106) * m(6) + ((t1019 * t544 + t1004 + t1016 - t625) * t440 + ((t1036 + t1043) * t544 + t1019 * t543 + t1005 + t1046) * t439 + t1029 + t1031) * t947 - (t110 + t123 + t84) * t791 / 0.2e1 + (t111 + t122) * t517 / 0.2e1 + (-qJD(3) * t648 + t461 * t559 + t462 * t556 + t538 * t711 + t539 * t712 + t546 * t709 + t547 * t710) * qJD(1) + t922 / 0.2e1 + (t1025 + t1026) * t946 + (((t310 * t544 + t311 * t543) * t538 + (t326 + t878) * t544 + t1006 + t1040 + t1042) * t440 + (-t544 * t715 - t312 * t862 + t880 + (t310 * t543 - t311 * t544) * t538 - t1016 + t1007) * t439 + t1032) * t944 + t912 + t913 / 0.2e1 + t914 / 0.2e1 + (t96 * (-t978 + t1037) + (-t544 * t926 + t550 + t708 + t809) * t95 + (t975 * t543 + t806 + (-t528 - t550 + (-t452 - t540) * t544) * qJD(1)) * t131 + (t812 + (-t771 - t975) * t544 + t131 - t617 - t980 + (t337 + t938 + t1037) * qJD(1)) * t132) * m(5) - t923 / 0.2e1 + (t950 + t951) * t34 + (t1024 - t1027 + t1028) * t945 + t44 * t951 + t43 * t952 + t125 * t955 + t124 * t956; m(4) * t107 + m(5) * t53 + m(6) * t28 + m(7) * t13; t564 + ((t556 * t805 + t559 * t804) * qJD(1) + ((t543 * t820 - t544 * t821) * t559 + (t543 * t822 + t544 * t823) * t556) * qJD(3)) * t932 + (-t110 * t544 + t111 * t543 + (t227 * t543 + t228 * t544) * qJD(1)) * t931 + ((-t517 * t866 + t796) * t543 + (t590 + (-t962 * t544 + (t867 + t603) * t543) * qJD(3)) * t544) * t742 + ((-t791 * t867 - t796) * t544 + (t590 + (t603 * t543 + (t866 - t962) * t544) * qJD(3)) * t543) * t739 + (qJD(1) * t122 + (t543 * (t543 * t967 + t571 * t544) - t544 * (t543 * t966 + t572 * t544) + (t143 * t543 + t144 * t544) * qJD(1)) * t984) * t941 + (qJD(1) * t123 + (t543 * (t571 * t543 - t544 * t967) - t544 * (t572 * t543 - t544 * t966) + (t141 * t543 + t142 * t544) * qJD(1)) * t984) * t940 + (t619 + t83) * t744 + (t618 + t84) * t743 + (t13 * (t639 + t832) + t47 * (t599 + t762) + (t77 * t600 + (-t315 + t763) * t887) * t543 - t77 * (t582 + t700) - t47 * (-t220 * t750 + t604 + t719) - (-t77 * t756 + (t47 * t729 - t559 * t682) * qJD(3)) * pkin(3) + t961 * (t478 - t995) + (t544 * t600 + t797 + t976) * t76) * m(7) + (t105 * t819 + t28 * (t704 + t832) + t80 * (t640 + t762) + (t105 * t608 + t643 * t983) * t544 + (t61 * t643 + t106 * t608 + (-t315 - t839) * t885) * t543 - t105 * (t690 - t797) - t106 * (t700 + t714) - t80 * (t718 + t719) - (-t106 * t756 + ((-t105 * t544 - t106 * t543) * t559 + t80 * t729) * qJD(3)) * pkin(3)) * m(6) + (t53 * (t831 + t832) + t120 * (t762 + t764) + (t131 * t692 + (qJD(1) * t132 + t96) * t735) * t544 + (t95 * t735 + t132 * t692 + (t131 * t451 + t120 * (-t315 - t338)) * qJD(1)) * t543 - (-t132 * t756 + (t120 * t729 + t559 * t663) * qJD(3)) * pkin(3) + t972) * m(5) + (t107 * t651 + t191 * t585 + t662 * t467 + (-t138 * t543 - t139 * t544 + (-t195 * t544 + t883) * qJD(1)) * t490 - (t194 * t417 - t881) * qJD(1) - (t191 * (-t417 * t543 - t418 * t544) + t662 * t687) * qJD(3)) * m(4); t564 + (t13 * t639 + (t647 * t77 + t763 * t887) * t543 - t77 * (t611 + t582) + t961 * (-t339 + t732) + (t544 * t647 - t466 + t976) * t76 + (t599 + t369 - (-t220 * t787 - t404) * t543 - t604) * t47) * m(7) + (t28 * t704 + (t106 * t699 - t839 * t885) * t543 - t106 * (t611 + t714) + (t404 * t543 + t369 + t640 - t718) * t80 + (t61 * t543 + t544 * t983) * t733 + (t544 * t699 - t466 - t690 + t819) * t105) * m(6) + (t53 * t831 + t120 * (-t338 * t795 + t764) + t663 * t410 + (-t95 * t543 - t96 * t544 + (t131 * t543 - t132 * t544) * qJD(1)) * t451 + t972) * m(5); 0.2e1 * (-t921 / 0.2e1 + t920 / 0.2e1 + t47 * t954) * m(7) + 0.2e1 * (t61 * t940 + t62 * t941 + t80 * t954) * m(6); -t34 * t759 / 0.2e1 + t3 * t864 / 0.2e1 + (-t125 * t539 + t538 * t678) * t955 + ((t552 * t678 - t43) * t539 + (-qJD(1) * t679 + t125 * t552 + t20 * t543 + t21 * t544) * t538) * t952 + t538 * t33 * t743 + t4 * t865 / 0.2e1 + (-t124 * t539 + t538 * t680) * t956 + ((t552 * t680 - t44) * t539 + (-qJD(1) * t681 + t124 * t552 + t22 * t543 + t23 * t544) * t538) * t951 + t37 * t748 - t539 * (t912 + t913 + t914 + t922 - t923) / 0.2e1 + (-t137 * t539 + t538 * t676) * t702 + ((t552 * t676 - t55) * t539 + (-qJD(1) * t677 + t137 * t552 + t24 * t543 + t25 * t544) * t538) * t942 + (t391 * t959 + t587 * t392 - t544 * t598) * t953 + (-t389 * t959 + t390 * t587 - t543 * t598) * t950 + (t609 * t539 + (-t555 * t959 + t558 * t587) * t538) * t943 + t986 * t860 / 0.2e1 + ((-t77 * t118 + t76 * t119 + t27 * t218 - t26 * t220 + (t47 * t659 + (t543 * t76 - t544 * t77) * t339) * t552) * t539 + (t76 * (t217 * t543 - t218 * t552) + t77 * (-t217 * t544 + t220 * t552) + t13 * t659 + t47 * (-t118 * t543 + t119 * t544 - t218 * t795 - t220 * t794) + (qJD(1) * t682 + t920 - t921) * t339) * t538 - t76 * (-t256 * t483 - t346 * t402) - t77 * (t257 * t483 - t345 * t402) - t47 * (t256 * t345 + t257 * t346)) * m(7);];
tauc  = t1(:);
