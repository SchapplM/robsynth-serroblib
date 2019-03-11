% Calculate vector of inverse dynamics joint torques for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PPRRRP1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_invdynJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_invdynJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRP1_invdynJ_fixb_slag_vp1: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_invdynJ_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP1_invdynJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRRP1_invdynJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPRRRP1_invdynJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:52:22
% EndTime: 2019-03-08 18:54:23
% DurationCPUTime: 112.56s
% Computational Cost: add. (146583->1526), mult. (421760->2073), div. (0->0), fcn. (531773->14), ass. (0->531)
t955 = Icges(6,4) + Icges(7,4);
t947 = Icges(6,1) + Icges(7,1);
t946 = Icges(6,5) + Icges(7,5);
t954 = Icges(6,2) + Icges(7,2);
t953 = Icges(6,6) + Icges(7,6);
t952 = Icges(6,3) + Icges(7,3);
t590 = sin(pkin(11));
t592 = cos(pkin(11));
t809 = sin(pkin(12));
t813 = cos(pkin(6));
t676 = t813 * t809;
t811 = cos(pkin(12));
t576 = t590 * t811 + t592 * t676;
t678 = t813 * t811;
t624 = t590 * t809 - t592 * t678;
t812 = cos(pkin(7));
t607 = t624 * t812;
t591 = sin(pkin(6));
t810 = sin(pkin(7));
t694 = t591 * t810;
t821 = sin(qJ(3));
t661 = t821 * t694;
t823 = cos(qJ(3));
t537 = t576 * t823 - t592 * t661 - t607 * t821;
t695 = t591 * t812;
t567 = -t592 * t695 + t624 * t810;
t595 = sin(qJ(4));
t822 = cos(qJ(4));
t465 = t537 * t822 + t567 * t595;
t662 = t823 * t694;
t536 = t576 * t821 + t592 * t662 + t607 * t823;
t594 = sin(qJ(5));
t596 = cos(qJ(5));
t353 = -t465 * t594 + t536 * t596;
t971 = t955 * t353;
t577 = -t590 * t676 + t592 * t811;
t623 = t590 * t678 + t592 * t809;
t606 = t623 * t812;
t539 = t577 * t823 + t590 * t661 - t606 * t821;
t568 = t590 * t695 + t623 * t810;
t467 = t539 * t822 + t568 * t595;
t538 = t577 * t821 - t590 * t662 + t606 * t823;
t355 = -t467 * t594 + t538 * t596;
t970 = t955 * t355;
t675 = t812 * t811;
t677 = t813 * t810;
t566 = t821 * t677 + (t675 * t821 + t809 * t823) * t591;
t575 = -t694 * t811 + t812 * t813;
t541 = t566 * t822 + t575 * t595;
t893 = t591 * (t823 * t675 - t821 * t809) + t823 * t677;
t468 = -t541 * t594 - t596 * t893;
t969 = t955 * t468;
t770 = t893 * t594;
t469 = t541 * t596 - t770;
t968 = t955 * t469;
t777 = t538 * t594;
t356 = t467 * t596 + t777;
t967 = t955 * t356;
t781 = t536 * t594;
t354 = t465 * t596 + t781;
t966 = t955 * t354;
t507 = t536 * qJD(3);
t655 = -t537 * t595 + t567 * t822;
t338 = qJD(4) * t655 - t507 * t822;
t508 = t537 * qJD(3);
t163 = -qJD(5) * t354 - t338 * t594 + t508 * t596;
t628 = qJD(5) * t353 + t508 * t594;
t164 = t338 * t596 + t628;
t337 = qJD(4) * t465 - t507 * t595;
t965 = -t163 * t953 - t164 * t946 - t337 * t952;
t509 = t538 * qJD(3);
t654 = -t539 * t595 + t568 * t822;
t340 = qJD(4) * t654 - t509 * t822;
t510 = t539 * qJD(3);
t165 = -qJD(5) * t356 - t340 * t594 + t510 * t596;
t627 = qJD(5) * t355 + t510 * t594;
t166 = t340 * t596 + t627;
t339 = qJD(4) * t467 - t509 * t595;
t964 = -t165 * t953 - t166 * t946 - t339 * t952;
t963 = t163 * t954 + t164 * t955 + t337 * t953;
t962 = t165 * t954 + t166 * t955 + t339 * t953;
t961 = t955 * t163 + t164 * t947 + t946 * t337;
t960 = t955 * t165 + t166 * t947 + t946 * t339;
t548 = t893 * qJD(3);
t653 = -t566 * t595 + t575 * t822;
t450 = qJD(4) * t653 + t548 * t822;
t549 = t566 * qJD(3);
t268 = -qJD(5) * t469 - t450 * t594 + t549 * t596;
t626 = qJD(5) * t468 + t549 * t594;
t269 = t450 * t596 + t626;
t449 = qJD(4) * t541 + t548 * t595;
t959 = -t268 * t953 - t269 * t946 - t449 * t952;
t958 = t268 * t954 + t269 * t955 + t449 * t953;
t957 = t955 * t268 + t269 * t947 + t946 * t449;
t925 = t353 * t953 + t354 * t946 - t655 * t952;
t924 = t355 * t953 + t356 * t946 - t654 * t952;
t923 = t353 * t954 - t655 * t953 + t966;
t922 = t355 * t954 - t654 * t953 + t967;
t921 = t354 * t947 - t946 * t655 + t971;
t920 = t356 * t947 - t946 * t654 + t970;
t907 = t468 * t953 + t469 * t946 - t653 * t952;
t906 = t468 * t954 - t653 * t953 + t968;
t905 = t469 * t947 - t946 * t653 + t969;
t938 = t923 * t163 + t921 * t164 + t925 * t337 + t963 * t353 + t961 * t354 + t655 * t965;
t937 = t163 * t922 + t164 * t920 + t337 * t924 + t353 * t962 + t354 * t960 + t655 * t964;
t936 = t923 * t165 + t921 * t166 + t925 * t339 + t963 * t355 + t961 * t356 + t654 * t965;
t935 = t165 * t922 + t166 * t920 + t339 * t924 + t355 * t962 + t356 * t960 + t654 * t964;
t932 = t923 * t268 + t921 * t269 + t925 * t449 + t963 * t468 + t961 * t469 + t653 * t965;
t931 = t268 * t922 + t269 * t920 + t449 * t924 + t468 * t962 + t469 * t960 + t653 * t964;
t929 = t163 * t906 + t164 * t905 + t337 * t907 + t353 * t958 + t354 * t957 + t655 * t959;
t928 = t165 * t906 + t166 * t905 + t339 * t907 + t355 * t958 + t356 * t957 + t654 * t959;
t927 = t268 * t906 + t269 * t905 + t449 * t907 + t468 * t958 + t469 * t957 + t653 * t959;
t892 = t353 * t923 + t354 * t921 - t655 * t925;
t891 = t353 * t922 + t354 * t920 - t655 * t924;
t890 = t355 * t923 + t356 * t921 - t654 * t925;
t889 = t355 * t922 + t356 * t920 - t654 * t924;
t888 = t468 * t923 + t469 * t921 - t653 * t925;
t887 = t468 * t922 + t469 * t920 - t653 * t924;
t886 = t353 * t906 + t354 * t905 - t655 * t907;
t885 = t355 * t906 + t356 * t905 - t654 * t907;
t884 = t468 * t906 + t469 * t905 - t653 * t907;
t593 = -qJ(6) - pkin(10);
t956 = rSges(7,3) - t593;
t951 = t594 * t953 - t596 * t946;
t950 = -t594 * t954 + t596 * t955;
t949 = -t594 * t955 + t596 * t947;
t948 = rSges(7,1) + pkin(5);
t926 = -pkin(10) + t956;
t943 = 0.2e1 * qJD(3);
t942 = 2 * qJDD(3);
t560 = qJDD(3) * t567;
t359 = qJD(4) * t508 + qJDD(4) * t536 + t560;
t157 = qJD(5) * t337 - qJDD(5) * t655 + t359;
t561 = qJDD(3) * t568;
t360 = qJD(4) * t510 + qJDD(4) * t538 + t561;
t158 = qJD(5) * t339 - qJDD(5) * t654 + t360;
t573 = qJDD(3) * t575;
t472 = qJD(4) * t549 - qJDD(4) * t893 + t573;
t265 = qJD(5) * t449 - qJDD(5) * t653 + t472;
t563 = qJD(3) * t567;
t470 = qJD(4) * t536 + t563;
t305 = -qJD(5) * t655 + t470;
t564 = qJD(3) * t568;
t471 = qJD(4) * t538 + t564;
t306 = -qJD(5) * t654 + t471;
t574 = qJD(3) * t575;
t542 = -qJD(4) * t893 + t574;
t430 = -qJD(5) * t653 + t542;
t941 = t892 * t157 + t891 * t158 + t886 * t265 + t305 * t938 + t937 * t306 + t929 * t430;
t940 = t157 * t890 + t158 * t889 + t265 * t885 + t305 * t936 + t306 * t935 + t430 * t928;
t939 = t157 * t888 + t158 * t887 + t265 * t884 + t305 * t932 + t306 * t931 + t430 * t927;
t934 = t305 * t892 + t306 * t891 + t430 * t886;
t933 = t305 * t890 + t306 * t889 + t430 * t885;
t930 = t305 * t888 + t306 * t887 + t430 * t884;
t703 = t594 * t822;
t404 = t536 * t703 + t537 * t596;
t702 = t596 * t822;
t778 = t537 * t594;
t405 = -t536 * t702 + t778;
t780 = t536 * t595;
t919 = -t404 * t953 - t405 * t946 + t780 * t952;
t406 = t538 * t703 + t539 * t596;
t774 = t539 * t594;
t407 = -t538 * t702 + t774;
t776 = t538 * t595;
t918 = -t406 * t953 - t407 * t946 + t776 * t952;
t917 = t404 * t954 + t405 * t955 - t780 * t953;
t916 = t406 * t954 + t407 * t955 - t776 * t953;
t915 = t404 * t955 + t405 * t947 - t780 * t946;
t914 = t406 * t955 + t407 * t947 - t776 * t946;
t913 = t465 * t953 + t655 * t950;
t912 = t467 * t953 + t654 * t950;
t911 = t465 * t946 + t655 * t949;
t910 = t467 * t946 + t654 * t949;
t586 = pkin(5) * t596 + pkin(4);
t787 = t655 * t596;
t788 = t655 * t594;
t909 = rSges(7,1) * t787 - rSges(7,2) * t788 + t465 * t956 + t586 * t655;
t784 = t654 * t596;
t785 = t654 * t594;
t908 = rSges(7,1) * t784 - rSges(7,2) * t785 + t467 * t956 + t586 * t654;
t300 = pkin(4) * t465 - pkin(10) * t655;
t819 = -pkin(4) + t586;
t762 = rSges(7,1) * t354 + rSges(7,2) * t353 + pkin(5) * t781 + t465 * t819 - t655 * t926;
t710 = -t300 - t762;
t302 = pkin(4) * t467 - pkin(10) * t654;
t761 = rSges(7,1) * t356 + rSges(7,2) * t355 + pkin(5) * t777 + t467 * t819 - t654 * t926;
t709 = -t302 - t761;
t486 = t566 * t596 - t703 * t893;
t767 = t566 * t594;
t487 = t702 * t893 + t767;
t769 = t893 * t595;
t904 = -t486 * t953 - t487 * t946 - t769 * t952;
t903 = t486 * t954 + t487 * t955 + t769 * t953;
t902 = t486 * t955 + t487 * t947 + t769 * t946;
t901 = t541 * t953 + t653 * t950;
t900 = t541 * t946 + t653 * t949;
t772 = t653 * t596;
t773 = t653 * t594;
t899 = rSges(7,1) * t772 - rSges(7,2) * t773 + t541 * t956 + t586 * t653;
t429 = pkin(4) * t541 - pkin(10) * t653;
t750 = rSges(7,1) * t469 + rSges(7,2) * t468 - pkin(5) * t770 + t541 * t819 - t653 * t926;
t708 = -t429 - t750;
t898 = (-t541 * t952 - t594 * t906 + t596 * t905 + t653 * t951) * t430 + (-t467 * t952 - t594 * t922 + t596 * t920 + t654 * t951) * t306 + (-t465 * t952 - t594 * t923 + t596 * t921 + t655 * t951) * t305;
t897 = -t586 * t822 + t593 * t595;
t374 = Icges(5,5) * t541 + Icges(5,6) * t653 - Icges(5,3) * t893;
t803 = Icges(5,4) * t541;
t375 = Icges(5,2) * t653 - Icges(5,6) * t893 + t803;
t532 = Icges(5,4) * t653;
t376 = Icges(5,1) * t541 - Icges(5,5) * t893 + t532;
t112 = t374 * t536 + t375 * t655 + t376 * t465;
t167 = Icges(5,5) * t338 - Icges(5,6) * t337 + Icges(5,3) * t508;
t169 = Icges(5,4) * t338 - Icges(5,2) * t337 + Icges(5,6) * t508;
t171 = Icges(5,1) * t338 - Icges(5,4) * t337 + Icges(5,5) * t508;
t247 = Icges(5,5) * t465 + Icges(5,6) * t655 + Icges(5,3) * t536;
t805 = Icges(5,4) * t465;
t249 = Icges(5,2) * t655 + Icges(5,6) * t536 + t805;
t451 = Icges(5,4) * t655;
t251 = Icges(5,1) * t465 + Icges(5,5) * t536 + t451;
t46 = t167 * t536 + t169 * t655 + t171 * t465 + t247 * t508 - t249 * t337 + t251 * t338;
t168 = Icges(5,5) * t340 - Icges(5,6) * t339 + Icges(5,3) * t510;
t170 = Icges(5,4) * t340 - Icges(5,2) * t339 + Icges(5,6) * t510;
t172 = Icges(5,1) * t340 - Icges(5,4) * t339 + Icges(5,5) * t510;
t248 = Icges(5,5) * t467 + Icges(5,6) * t654 + Icges(5,3) * t538;
t804 = Icges(5,4) * t467;
t250 = Icges(5,2) * t654 + Icges(5,6) * t538 + t804;
t452 = Icges(5,4) * t654;
t252 = Icges(5,1) * t467 + Icges(5,5) * t538 + t452;
t47 = t168 * t536 + t170 * t655 + t172 * t465 + t248 * t508 - t250 * t337 + t252 * t338;
t270 = Icges(5,5) * t450 - Icges(5,6) * t449 + Icges(5,3) * t549;
t271 = Icges(5,4) * t450 - Icges(5,2) * t449 + Icges(5,6) * t549;
t272 = Icges(5,1) * t450 - Icges(5,4) * t449 + Icges(5,5) * t549;
t56 = t270 * t536 + t271 * t655 + t272 * t465 - t337 * t375 + t338 * t376 + t374 * t508;
t77 = t247 * t536 + t249 * t655 + t251 * t465;
t78 = t248 * t536 + t250 * t655 + t252 * t465;
t896 = t112 * t472 + t359 * t77 + t360 * t78 + t46 * t470 + t47 * t471 + t542 * t56 + t941;
t113 = t374 * t538 + t375 * t654 + t376 * t467;
t48 = t167 * t538 + t169 * t654 + t171 * t467 + t247 * t510 - t249 * t339 + t251 * t340;
t49 = t168 * t538 + t170 * t654 + t172 * t467 + t248 * t510 - t250 * t339 + t252 * t340;
t57 = t270 * t538 + t271 * t654 + t272 * t467 - t339 * t375 + t340 * t376 + t374 * t510;
t79 = t247 * t538 + t249 * t654 + t251 * t467;
t80 = t248 * t538 + t250 * t654 + t252 * t467;
t895 = t113 * t472 + t359 * t79 + t360 * t80 + t470 * t48 + t471 * t49 + t542 * t57 + t940;
t106 = -t247 * t893 + t249 * t653 + t251 * t541;
t107 = -t248 * t893 + t250 * t653 + t252 * t541;
t122 = -t374 * t893 + t375 * t653 + t376 * t541;
t51 = -t167 * t893 + t169 * t653 + t171 * t541 + t247 * t549 - t249 * t449 + t251 * t450;
t52 = -t168 * t893 + t170 * t653 + t172 * t541 + t248 * t549 - t250 * t449 + t252 * t450;
t68 = -t270 * t893 + t271 * t653 + t272 * t541 + t374 * t549 - t375 * t449 + t376 * t450;
t894 = t106 * t359 + t107 * t360 + t122 * t472 + t470 * t51 + t471 * t52 + t542 * t68 + t939;
t299 = pkin(4) * t655 + pkin(10) * t465;
t883 = -t299 + t909;
t301 = pkin(4) * t654 + pkin(10) * t467;
t882 = -t301 + t908;
t428 = pkin(4) * t653 + pkin(10) * t541;
t881 = -t428 + t899;
t880 = (-t469 * t954 + t905 + t969) * t430 + (-t356 * t954 + t920 + t970) * t306 + (-t354 * t954 + t921 + t971) * t305;
t879 = (t468 * t947 - t906 - t968) * t430 + (t355 * t947 - t922 - t967) * t306 + (t353 * t947 - t923 - t966) * t305;
t878 = (t468 * t946 - t469 * t953) * t430 + (t355 * t946 - t356 * t953) * t306 + (t353 * t946 - t354 * t953) * t305;
t392 = -rSges(4,1) * t507 - rSges(4,2) * t508;
t393 = -rSges(4,1) * t509 - rSges(4,2) * t510;
t665 = t392 * t568 - t393 * t567;
t370 = rSges(4,1) * t537 - rSges(4,2) * t536 + rSges(4,3) * t567;
t371 = rSges(4,1) * t539 - rSges(4,2) * t538 + rSges(4,3) * t568;
t668 = t370 * t568 - t371 * t567;
t876 = t665 * qJD(3) + t668 * qJDD(3);
t875 = t487 * rSges(7,1) + t486 * rSges(7,2) + rSges(7,3) * t769 + pkin(5) * t767 - t893 * t897;
t874 = t407 * rSges(7,1) + t406 * rSges(7,2) - rSges(7,3) * t776 + pkin(5) * t774 + t538 * t897;
t873 = t405 * rSges(7,1) + t404 * rSges(7,2) - rSges(7,3) * t780 + pkin(5) * t778 + t536 * t897;
t716 = qJDD(2) * t813 + qJDD(1);
t872 = -m(3) - m(4) - m(5) - m(6) - m(7);
t727 = qJD(6) * t654;
t815 = rSges(7,1) * t166 + rSges(7,2) * t165 + pkin(5) * t627 + t339 * t926 + t340 * t819 - t727;
t728 = qJD(6) * t655;
t817 = rSges(7,1) * t164 + rSges(7,2) * t163 + pkin(5) * t628 + t337 * t926 + t338 * t819 - t728;
t871 = -t761 * t157 + t762 * t158 - t815 * t305 + t817 * t306;
t185 = pkin(4) * t338 + pkin(10) * t337;
t140 = t471 * t185;
t182 = t360 * t300;
t394 = -pkin(3) * t507 + pkin(9) * t508;
t357 = t394 * t564;
t421 = pkin(3) * t537 + pkin(9) * t536;
t369 = t421 * t561;
t395 = -pkin(3) * t509 + pkin(9) * t510;
t423 = pkin(3) * t539 + pkin(9) * t538;
t620 = t357 + t369 + (-qJD(3) * t395 - qJDD(3) * t423) * t567 + t716;
t186 = pkin(4) * t340 + pkin(10) * t339;
t782 = t470 * t186;
t789 = t359 * t302;
t603 = t140 + t182 + t620 - t782 - t789;
t718 = qJDD(6) * t653;
t729 = qJD(6) * t449;
t7 = t603 - t718 + t729 + t871;
t870 = 0.2e1 * t7;
t284 = pkin(4) * t450 + pkin(10) * t449;
t481 = pkin(3) * t548 + pkin(9) * t549;
t495 = pkin(3) * t566 - pkin(9) * t893;
t719 = qJDD(2) * t591;
t583 = t590 * t719;
t621 = t481 * t563 + t495 * t560 + t583 + (-qJD(3) * t394 - qJDD(3) * t421) * t575;
t604 = -t185 * t542 + t470 * t284 - t300 * t472 + t359 * t429 + t621;
t726 = qJD(6) * t653;
t764 = rSges(7,1) * t269 + rSges(7,2) * t268 + pkin(5) * t626 + t449 * t926 + t450 * t819 - t726;
t11 = qJD(6) * t339 - qJDD(6) * t654 + t157 * t750 - t265 * t762 + t305 * t764 - t430 * t817 + t604;
t869 = 0.2e1 * t11;
t699 = t592 * t719;
t605 = t395 * t574 + t423 * t573 + (-qJD(3) * t481 - qJDD(3) * t495) * t568 - t699;
t601 = t542 * t186 - t471 * t284 + t472 * t302 - t360 * t429 + t605;
t12 = qJD(6) * t337 - qJDD(6) * t655 - t158 * t750 + t265 * t761 - t306 * t764 + t430 * t815 + t601;
t868 = 0.2e1 * t12;
t724 = qJD(2) * t813 + qJD(1);
t647 = t421 * t564 - t423 * t563 + t724;
t622 = t471 * t300 - t302 * t470 + t647;
t38 = -t305 * t761 + t306 * t762 + t622 - t726;
t867 = 0.2e1 * t38;
t734 = qJD(2) * t591;
t584 = t590 * t734;
t657 = -t421 * t574 + t495 * t563 + t584;
t635 = -t300 * t542 + t470 * t429 + t657;
t43 = t305 * t750 - t430 * t762 + t635 - t727;
t866 = 0.2e1 * t43;
t701 = t592 * t734;
t636 = t423 * t574 - t495 * t564 - t701;
t619 = t542 * t302 - t429 * t471 + t636;
t44 = -t306 * t750 + t430 * t761 + t619 - t728;
t865 = 0.2e1 * t44;
t864 = m(3) / 0.2e1;
t863 = m(4) / 0.2e1;
t862 = m(5) / 0.2e1;
t861 = m(6) / 0.2e1;
t860 = m(7) / 0.2e1;
t859 = t157 / 0.2e1;
t858 = t158 / 0.2e1;
t857 = t265 / 0.2e1;
t856 = -t305 / 0.2e1;
t855 = t305 / 0.2e1;
t854 = -t306 / 0.2e1;
t853 = t306 / 0.2e1;
t850 = t359 / 0.2e1;
t849 = t360 / 0.2e1;
t846 = -t430 / 0.2e1;
t845 = t430 / 0.2e1;
t841 = -t470 / 0.2e1;
t840 = t470 / 0.2e1;
t839 = -t471 / 0.2e1;
t838 = t471 / 0.2e1;
t837 = t472 / 0.2e1;
t830 = -t542 / 0.2e1;
t829 = t542 / 0.2e1;
t99 = rSges(6,1) * t164 + rSges(6,2) * t163 + rSges(6,3) * t337;
t816 = t306 * t99;
t814 = -t185 - t99;
t808 = Icges(4,4) * t537;
t807 = Icges(4,4) * t539;
t806 = Icges(4,4) * t566;
t101 = rSges(6,1) * t166 + rSges(6,2) * t165 + rSges(6,3) * t339;
t796 = t101 * t305;
t154 = rSges(6,1) * t354 + rSges(6,2) * t353 - rSges(6,3) * t655;
t795 = t154 * t158;
t156 = rSges(6,1) * t356 + rSges(6,2) * t355 - rSges(6,3) * t654;
t794 = t156 * t157;
t173 = rSges(5,1) * t338 - rSges(5,2) * t337 + rSges(5,3) * t508;
t793 = t173 * t471;
t174 = rSges(5,1) * t340 - rSges(5,2) * t339 + rSges(5,3) * t510;
t792 = t174 * t470;
t253 = rSges(5,1) * t465 + rSges(5,2) * t655 + rSges(5,3) * t536;
t791 = t253 * t360;
t254 = rSges(5,1) * t467 + rSges(5,2) * t654 + rSges(5,3) * t538;
t790 = t254 * t359;
t765 = -t101 - t186;
t130 = rSges(6,1) * t269 + rSges(6,2) * t268 + rSges(6,3) * t449;
t763 = -t130 - t284;
t760 = -t154 - t300;
t759 = -t156 - t302;
t758 = t538 * t185 + t510 * t300;
t757 = -t186 * t893 + t549 * t302;
t358 = t568 * t394;
t756 = t568 * t185 + t358;
t381 = t575 * t395;
t755 = t575 * t186 + t381;
t754 = -t354 * rSges(7,2) + t353 * t948;
t753 = -t356 * rSges(7,2) + t355 * t948;
t707 = t536 * t822;
t738 = pkin(4) * t707 + pkin(10) * t780;
t752 = t738 + t873;
t706 = t538 * t822;
t737 = pkin(4) * t706 + pkin(10) * t776;
t751 = t737 + t874;
t749 = t536 * t284 + t508 * t429;
t263 = rSges(6,1) * t469 + rSges(6,2) * t468 - rSges(6,3) * t653;
t748 = -t263 - t429;
t432 = t567 * t481;
t747 = t567 * t284 + t432;
t378 = t568 * t421;
t746 = t568 * t300 + t378;
t385 = t575 * t423;
t745 = t575 * t302 + t385;
t744 = -t469 * rSges(7,2) + t468 * t948;
t705 = t893 * t822;
t736 = -pkin(4) * t705 - pkin(10) * t769;
t743 = t736 + t875;
t446 = t567 * t495;
t742 = t567 * t429 + t446;
t735 = 0.2e1 * t716;
t733 = qJD(4) * t537;
t732 = qJD(4) * t539;
t731 = qJD(4) * t566;
t730 = qJD(5) * t595;
t725 = qJD(6) * t595;
t723 = 0.2e1 * m(4);
t722 = 0.2e1 * m(5);
t721 = 0.2e1 * m(6);
t720 = 0.2e1 * m(7);
t715 = 0.2e1 * t813;
t714 = 0.2e1 * t591;
t713 = -t185 - t817;
t712 = -t186 - t815;
t711 = -t284 - t764;
t218 = t405 * rSges(6,1) + t404 * rSges(6,2) - rSges(6,3) * t780;
t220 = t407 * rSges(6,1) + t406 * rSges(6,2) - rSges(6,3) * t776;
t315 = t487 * rSges(6,1) + t486 * rSges(6,2) + rSges(6,3) * t769;
t420 = -t536 * pkin(3) + pkin(9) * t537;
t422 = -t538 * pkin(3) + pkin(9) * t539;
t494 = pkin(3) * t893 + pkin(9) * t566;
t692 = -t722 / 0.2e1;
t691 = t722 / 0.2e1;
t690 = -t721 / 0.2e1;
t689 = t721 / 0.2e1;
t688 = -t720 / 0.2e1;
t685 = t471 * t299 - t301 * t470;
t684 = t542 * t301 - t428 * t471;
t683 = -t299 * t542 + t470 * t428;
t681 = t420 * t564 - t422 * t563;
t680 = t422 * t574 - t494 * t564;
t679 = -t420 * t574 + t494 * t563;
t476 = rSges(4,1) * t566 + rSges(4,2) * t893 + rSges(4,3) * t575;
t667 = -t370 * t575 + t476 * t567;
t666 = t371 * t575 - t476 * t568;
t480 = rSges(4,1) * t548 - rSges(4,2) * t549;
t664 = -t392 * t575 + t480 * t567;
t663 = t393 * t575 - t480 * t568;
t322 = -rSges(5,1) * t707 + rSges(5,2) * t780 + t537 * rSges(5,3);
t323 = -rSges(5,1) * t706 + rSges(5,2) * t776 + t539 * rSges(5,3);
t436 = rSges(5,1) * t705 - rSges(5,2) * t769 + t566 * rSges(5,3);
t236 = rSges(6,1) * t787 - rSges(6,2) * t788 + t465 * rSges(6,3);
t238 = rSges(6,1) * t784 - rSges(6,2) * t785 + t467 * rSges(6,3);
t331 = rSges(6,1) * t772 - rSges(6,2) * t773 + t541 * rSges(6,3);
t656 = -0.2e1 * t395 * t563 - 0.2e1 * t423 * t560 + 0.2e1 * t357 + 0.2e1 * t369 + t735;
t652 = -Icges(5,1) * t822 + Icges(5,4) * t595;
t651 = -Icges(5,4) * t822 + Icges(5,2) * t595;
t650 = -Icges(5,5) * t822 + Icges(5,6) * t595;
t643 = t38 * t817 + t7 * t762;
t640 = (Icges(5,5) * t655 - Icges(5,6) * t465) * t470 + (Icges(5,5) * t654 - Icges(5,6) * t467) * t471 + (Icges(5,5) * t653 - Icges(5,6) * t541) * t542;
t639 = (-Icges(4,5) * t536 - Icges(4,6) * t537) * t567 + (-Icges(4,5) * t538 - Icges(4,6) * t539) * t568 + (Icges(4,5) * t893 - Icges(4,6) * t566) * t575;
t638 = t12 * t761 + t44 * t815;
t637 = t11 * t750 + t43 * t764;
t625 = -0.2e1 * t782 + 0.2e1 * t140 - 0.2e1 * t789 + 0.2e1 * t182 + t656;
t614 = (Icges(5,1) * t654 - t250 - t804) * t471 + (Icges(5,1) * t655 - t249 - t805) * t470 + (Icges(5,1) * t653 - t375 - t803) * t542;
t613 = (Icges(5,2) * t467 - t252 - t452) * t471 + (Icges(5,2) * t465 - t251 - t451) * t470 + (Icges(5,2) * t541 - t376 - t532) * t542;
t365 = -Icges(4,2) * t536 + Icges(4,6) * t567 + t808;
t366 = -Icges(4,2) * t538 + Icges(4,6) * t568 + t807;
t474 = Icges(4,2) * t893 + Icges(4,6) * t575 + t806;
t612 = (-Icges(4,1) * t538 - t366 - t807) * t568 + (-Icges(4,1) * t536 - t365 - t808) * t567 + (Icges(4,1) * t893 - t474 - t806) * t575;
t526 = Icges(4,4) * t536;
t367 = Icges(4,1) * t537 + Icges(4,5) * t567 - t526;
t527 = Icges(4,4) * t538;
t368 = Icges(4,1) * t539 + Icges(4,5) * t568 - t527;
t558 = Icges(4,4) * t893;
t475 = Icges(4,1) * t566 + Icges(4,5) * t575 + t558;
t611 = (Icges(4,2) * t539 - t368 + t527) * t568 + (Icges(4,2) * t537 - t367 + t526) * t567 + (Icges(4,2) * t566 - t475 - t558) * t575;
t610 = t300 * t732 - t302 * t733 + t470 * t737 - t471 * t738 + t681;
t609 = t302 * t731 - t429 * t732 + t471 * t736 - t542 * t737 + t680;
t608 = -t300 * t731 + t429 * t733 - t470 * t736 + t542 * t738 + t679;
t598 = (Icges(5,3) * t539 + t250 * t595 - t252 * t822 + t538 * t650) * t471 + (Icges(5,3) * t537 + t249 * t595 - t251 * t822 + t536 * t650) * t470 + (Icges(5,3) * t566 + t375 * t595 - t376 * t822 - t650 * t893) * t542;
t493 = rSges(4,1) * t893 - rSges(4,2) * t566;
t489 = t730 * t893 + t731;
t479 = Icges(4,1) * t548 - Icges(4,4) * t549;
t478 = Icges(4,4) * t548 - Icges(4,2) * t549;
t477 = Icges(4,5) * t548 - Icges(4,6) * t549;
t473 = Icges(4,5) * t566 + Icges(4,6) * t893 + Icges(4,3) * t575;
t435 = Icges(5,5) * t566 - t652 * t893;
t434 = Icges(5,6) * t566 - t651 * t893;
t427 = rSges(5,1) * t653 - rSges(5,2) * t541;
t419 = -rSges(4,1) * t538 - rSges(4,2) * t539;
t418 = -rSges(4,1) * t536 - rSges(4,2) * t537;
t411 = -t538 * t730 + t732;
t410 = -t536 * t730 + t733;
t391 = -Icges(4,1) * t509 - Icges(4,4) * t510;
t390 = -Icges(4,1) * t507 - Icges(4,4) * t508;
t389 = -Icges(4,4) * t509 - Icges(4,2) * t510;
t388 = -Icges(4,4) * t507 - Icges(4,2) * t508;
t387 = -Icges(4,5) * t509 - Icges(4,6) * t510;
t386 = -Icges(4,5) * t507 - Icges(4,6) * t508;
t377 = rSges(5,1) * t541 + rSges(5,2) * t653 - rSges(5,3) * t893;
t364 = Icges(4,5) * t539 - Icges(4,6) * t538 + Icges(4,3) * t568;
t363 = Icges(4,5) * t537 - Icges(4,6) * t536 + Icges(4,3) * t567;
t335 = t536 * t429;
t321 = Icges(5,5) * t539 + t538 * t652;
t320 = Icges(5,5) * t537 + t536 * t652;
t319 = Icges(5,6) * t539 + t538 * t651;
t318 = Icges(5,6) * t537 + t536 * t651;
t304 = rSges(6,1) * t468 - rSges(6,2) * t469;
t292 = rSges(5,1) * t654 - rSges(5,2) * t467;
t291 = rSges(5,1) * t655 - rSges(5,2) * t465;
t278 = t893 * t302;
t274 = rSges(5,1) * t450 - rSges(5,2) * t449 + rSges(5,3) * t549;
t264 = t538 * t300;
t222 = qJD(3) * t666 - t701;
t221 = qJD(3) * t667 + t584;
t202 = rSges(6,1) * t355 - rSges(6,2) * t356;
t200 = rSges(6,1) * t353 - rSges(6,2) * t354;
t177 = qJD(3) * t668 + t724;
t133 = qJD(3) * t663 + qJDD(3) * t666 - t699;
t132 = qJD(3) * t664 + qJDD(3) * t667 + t583;
t114 = t716 + t876;
t111 = t254 * t542 - t377 * t471 + t636;
t110 = -t253 * t542 + t377 * t470 + t657;
t83 = t253 * t471 - t254 * t470 + t647;
t59 = t156 * t430 - t263 * t306 + t619;
t58 = -t154 * t430 + t263 * t305 + t635;
t55 = t174 * t542 + t254 * t472 - t274 * t471 - t360 * t377 + t605;
t54 = -t173 * t542 - t253 * t472 + t274 * t470 + t359 * t377 + t621;
t53 = t154 * t306 - t156 * t305 + t622;
t50 = t620 - t790 + t791 - t792 + t793;
t45 = t106 * t470 + t107 * t471 + t122 * t542;
t42 = t113 * t542 + t470 * t79 + t471 * t80;
t41 = t112 * t542 + t470 * t77 + t471 * t78;
t31 = t101 * t430 - t130 * t306 + t156 * t265 - t158 * t263 + t601;
t30 = t130 * t305 - t154 * t265 + t157 * t263 - t430 * t99 + t604;
t13 = t603 - t794 + t795 - t796 + t816;
t1 = [m(2) * qJDD(1) + (t656 - 0.2e1 * t790 + 0.2e1 * t791 - 0.2e1 * t792 + 0.2e1 * t793) * t862 + (t625 - 0.2e1 * t794 + 0.2e1 * t795 - 0.2e1 * t796 + 0.2e1 * t816) * t861 + (-m(2) + t872) * g(3) + (t864 + t863) * t735 + 0.2e1 * t876 * t863 + (t625 - 0.2e1 * t718 + 0.2e1 * t729 + 0.2e1 * t871) * t860; (t716 * t715 + 0.2e1 * (t590 ^ 2 + t592 ^ 2) * t591 ^ 2 * qJDD(2)) * t864 + (t114 * t715 + (t132 * t590 - t133 * t592) * t714) * t863 + (t50 * t715 + (t54 * t590 - t55 * t592) * t714) * t862 + (t13 * t715 + (t30 * t590 - t31 * t592) * t714) * t861 + (t7 * t715 + (t11 * t590 - t12 * t592) * t714) * t860 + t872 * (g(3) * t813 + (g(1) * t590 - g(2) * t592) * t591); (t884 * t489 + (t468 * t903 + t469 * t902 + t486 * t906 + t487 * t905 + t653 * t904 + t769 * t907) * t430 + t887 * t411 + t888 * t410 + (t468 * t916 + t469 * t914 + t486 * t922 + t487 * t920 + t653 * t918 + t769 * t924) * t306 + (t468 * t917 + t469 * t915 + t486 * t923 + t487 * t921 + t653 * t919 + t769 * t925) * t305) * t846 - m(7) * (g(1) * (t422 + t874) + g(2) * (t420 + t873) + g(3) * (t494 + t875)) + ((t539 * t248 + t319 * t654 + t467 * t321) * t471 + (t539 * t247 + t318 * t654 + t467 * t320) * t470 + (t539 * t374 + t434 * t654 + t467 * t435) * t542 + (t113 * t566 + t537 * t79 + t539 * t80) * qJD(4) + t598 * t538) * t839 + (t114 * t668 + t132 * t667 + t133 * t666 + t177 * t665 + t221 * t664 + t222 * t663) * t723 / 0.2e1 - (t221 * (-t418 * t575 + t493 * t567) + t222 * (t419 * t575 - t493 * t568) + t177 * (t418 * t568 - t419 * t567)) * qJD(3) * t723 / 0.2e1 - t930 * t489 / 0.2e1 + (t567 * t932 + t568 * t931 + t575 * t927) * t845 - t933 * t411 / 0.2e1 - t934 * t410 / 0.2e1 + (t567 * t936 + t568 * t935 + t575 * t928) * t853 + (t567 * t938 + t568 * t937 + t575 * t929) * t855 + (t30 * (t263 * t567 + (-t421 + t760) * t575 + t742) + t58 * (t130 * t567 + (-t394 + t814) * t575 + t747) + t31 * (t156 * t575 + (-t495 + t748) * t568 + t745) + t59 * (t101 * t575 + (-t481 + t763) * t568 + t755) + t13 * (t154 * t568 + (-t423 + t759) * t567 + t746) + t53 * (t568 * t99 + (-t395 + t765) * t567 + t756)) * t689 + ((t537 * t248 + t319 * t655 + t465 * t321) * t471 + (t537 * t247 + t318 * t655 + t465 * t320) * t470 + (t537 * t374 + t434 * t655 + t465 * t435) * t542 + (t112 * t566 + t537 * t77 + t539 * t78) * qJD(4) + t598 * t536) * t841 - t45 * t731 / 0.2e1 - t42 * t732 / 0.2e1 + (((t363 * t568 - t365 * t538 + t367 * t539) * t567 + (t364 * t568 - t366 * t538 + t368 * t539) * t568 + (t473 * t568 - t474 * t538 + t475 * t539) * t575) * t942 + ((-t365 * t510 - t367 * t509 + t386 * t568 - t388 * t538 + t390 * t539) * t567 + (-t366 * t510 - t368 * t509 + t387 * t568 - t389 * t538 + t391 * t539) * t568 + (-t474 * t510 - t475 * t509 + t477 * t568 - t478 * t538 + t479 * t539) * t575) * t943 + t895) * t568 / 0.2e1 + (((-t365 * t508 - t367 * t507 + t386 * t567 - t388 * t536 + t390 * t537) * t567 + (-t366 * t508 - t368 * t507 + t387 * t567 - t389 * t536 + t391 * t537) * t568 + (-t474 * t508 - t475 * t507 + t477 * t567 - t478 * t536 + t479 * t537) * t575) * t943 + ((t363 * t567 - t365 * t536 + t367 * t537) * t567 + (t364 * t567 - t366 * t536 + t368 * t537) * t568 + (t473 * t567 - t474 * t536 + t475 * t537) * t575) * t942 + t896) * t567 / 0.2e1 + (t58 * (-t154 * t489 - t218 * t430 + t263 * t410 + t305 * t315 + t608) + t59 * (t156 * t489 + t220 * t430 - t263 * t411 - t306 * t315 + t609) + t53 * (t154 * t411 - t156 * t410 + t218 * t306 - t220 * t305 + t610)) * t690 - m(5) * (g(1) * (t323 + t422) + g(2) * (t322 + t420) + g(3) * (t436 + t494)) - (t575 * (t566 * t612 + t575 * t639 - t611 * t893) + t568 * (t538 * t611 + t539 * t612 + t568 * t639) + t567 * (t536 * t611 + t537 * t612 + t567 * t639)) * qJD(3) ^ 2 / 0.2e1 + ((t566 * t248 + t319 * t653 + t541 * t321) * t471 + (t566 * t247 + t318 * t653 + t541 * t320) * t470 + (t566 * t374 + t434 * t653 + t541 * t435) * t542 + (t106 * t537 + t107 * t539 + t122 * t566) * qJD(4) - t598 * t893) * t830 + (((t363 * t575 + t365 * t893 + t367 * t566) * t567 + (t364 * t575 + t366 * t893 + t368 * t566) * t568 + (t473 * t575 + t474 * t893 + t475 * t566) * t575) * t942 + ((-t365 * t549 + t367 * t548 + t386 * t575 + t388 * t893 + t390 * t566) * t567 + (-t366 * t549 + t368 * t548 + t387 * t575 + t389 * t893 + t391 * t566) * t568 + (-t474 * t549 + t475 * t548 + t477 * t575 + t478 * t893 + t479 * t566) * t575) * t943 + t894) * t575 / 0.2e1 + (t43 * (t305 * t743 + t410 * t750 - t430 * t752 - t489 * t762 - t538 * t725 + t608) + t44 * (-t306 * t743 - t411 * t750 + t430 * t751 + t489 * t761 - t536 * t725 + t609) + t38 * (-t305 * t751 + t306 * t752 - t410 * t761 + t411 * t762 + t725 * t893 + t610)) * t688 + (t46 * t567 + t47 * t568 + t56 * t575) * t840 + (t113 * t575 + t567 * t79 + t568 * t80) * t849 + (t112 * t575 + t567 * t77 + t568 * t78) * t850 + (t51 * t567 + t52 * t568 + t575 * t68) * t829 + (t106 * t567 + t107 * t568 + t122 * t575) * t837 + (t48 * t567 + t49 * t568 + t57 * t575) * t838 - m(6) * (g(1) * (t422 + t220 - t737) + g(2) * (t420 + t218 - t738) + g(3) * (t494 + t315 - t736)) - t41 * t733 / 0.2e1 + (t54 * (t377 * t567 + t446 + (-t253 - t421) * t575) + t110 * (t274 * t567 + t432 + (-t173 - t394) * t575) + t55 * (t254 * t575 + t385 + (-t377 - t495) * t568) + t111 * (t174 * t575 + t381 + (-t274 - t481) * t568) + t50 * (t253 * t568 + t378 + (-t254 - t423) * t567) + t83 * (t173 * t568 + t358 + (-t174 - t395) * t567)) * t691 + (t110 * (-t322 * t542 + t436 * t470 + (-t253 * t566 + t377 * t537) * qJD(4) + t679) + t111 * (t323 * t542 - t436 * t471 + (t254 * t566 - t377 * t539) * qJD(4) + t680) + t83 * (t322 * t471 - t323 * t470 + (t253 * t539 - t254 * t537) * qJD(4) + t681)) * t692 + (t567 * t888 + t568 * t887 + t575 * t884) * t857 + (t567 * t890 + t568 * t889 + t575 * t885) * t858 + (t885 * t489 + (t355 * t903 + t356 * t902 + t406 * t906 + t407 * t905 + t654 * t904 - t776 * t907) * t430 + t889 * t411 + t890 * t410 + (t355 * t916 + t356 * t914 + t406 * t922 + t407 * t920 + t654 * t918 - t776 * t924) * t306 + (t355 * t917 + t356 * t915 + t406 * t923 + t407 * t921 + t654 * t919 - t776 * t925) * t305) * t854 + (t567 * t892 + t568 * t891 + t575 * t886) * t859 + (t886 * t489 + (t353 * t903 + t354 * t902 + t404 * t906 + t405 * t905 + t655 * t904 - t780 * t907) * t430 + t891 * t411 + t892 * t410 + (t353 * t916 + t354 * t914 + t404 * t922 + t405 * t920 + t655 * t918 - t780 * t924) * t306 + (t353 * t917 + t354 * t915 + t404 * t923 + t405 * t921 + t655 * t919 - t780 * t925) * t305) * t856 - m(4) * (g(1) * t419 + g(2) * t418 + g(3) * t493) + (t742 * t869 + t745 * t868 + t746 * t870 + t747 * t866 + t755 * t865 + t756 * t867 + 0.2e1 * (t11 * (-t421 + t710) + t43 * (-t394 + t713) + t638) * t575 + 0.2e1 * (t12 * (-t495 + t708) + t44 * (-t481 + t711) + t643) * t568 + 0.2e1 * (t7 * (-t423 + t709) + t38 * (-t395 + t712) + t637) * t567) * t860; -(t465 * t934 + t467 * t933 + t541 * t930) * qJD(5) / 0.2e1 + (t467 * t614 + t538 * t640 - t613 * t654) * t839 + (t45 + t930) * t549 / 0.2e1 + (t42 + t933) * t510 / 0.2e1 + (t41 + t934) * t508 / 0.2e1 + (t465 * t614 + t536 * t640 - t613 * t655) * t841 + t895 * t538 / 0.2e1 + t896 * t536 / 0.2e1 + (t54 * (t253 * t893 + t377 * t536) + t110 * (t173 * t893 - t253 * t549 + t274 * t536 + t377 * t508) + t55 * (-t254 * t893 - t377 * t538) + t111 * (-t174 * t893 + t254 * t549 - t274 * t538 - t377 * t510) + t50 * (t253 * t538 - t254 * t536) + t83 * (t173 * t538 - t174 * t536 + t253 * t510 - t254 * t508)) * t691 + (t536 * t888 + t538 * t887 - t884 * t893) * t857 + (t536 * t890 + t538 * t889 - t885 * t893) * t858 + (t536 * t892 + t538 * t891 - t886 * t893) * t859 - t894 * t893 / 0.2e1 + (t541 * t614 - t613 * t653 - t640 * t893) * t830 + (t508 * t888 + t510 * t887 + t536 * t932 + t538 * t931 + t549 * t884 - t893 * t927) * t845 + (t508 * t890 + t510 * t889 + t536 * t936 + t538 * t935 + t549 * t885 - t893 * t928) * t853 + (t508 * t892 + t510 * t891 + t536 * t938 + t538 * t937 + t549 * t886 - t893 * t929) * t855 + (t112 * t549 + t46 * t536 + t47 * t538 + t508 * t77 + t510 * t78 - t56 * t893) * t840 + (-t113 * t893 + t536 * t79 + t538 * t80) * t849 + (-t112 * t893 + t536 * t77 + t538 * t78) * t850 - m(6) * (g(1) * (t238 + t301) + g(2) * (t236 + t299) + g(3) * (t331 + t428)) + (t58 * (-t236 * t430 + t305 * t331 + (-t154 * t541 + t263 * t465) * qJD(5) + t683) + t59 * (t238 * t430 - t306 * t331 + (t156 * t541 - t263 * t467) * qJD(5) + t684) + t53 * (t236 * t306 - t238 * t305 + (t154 * t467 - t156 * t465) * qJD(5) + t685)) * t690 + (t106 * t508 + t107 * t510 + t122 * t549 + t51 * t536 + t52 * t538 - t68 * t893) * t829 + (t106 * t536 + t107 * t538 - t122 * t893) * t837 + (t113 * t549 + t48 * t536 + t49 * t538 + t508 * t79 + t510 * t80 - t57 * t893) * t838 + (t264 * t870 - t278 * t868 + t335 * t869 + t749 * t866 + t757 * t865 + t758 * t867 + (t710 * t866 + t761 * t865) * t549 + (t708 * t865 + t762 * t867) * t510 + (t709 * t867 + t750 * t866) * t508 - 0.2e1 * (t11 * t710 + t43 * t713 + t638) * t893 + 0.2e1 * (t12 * t708 + t44 * t711 + t643) * t538 + 0.2e1 * (t38 * t712 + t7 * t709 + t637) * t536) * t860 + (t30 * (t263 * t536 - t760 * t893 + t335) + t58 * (t130 * t536 + t263 * t508 + t549 * t760 - t814 * t893 + t749) + t31 * (-t156 * t893 + t538 * t748 - t278) + t59 * (-t101 * t893 + t156 * t549 + t510 * t748 + t538 * t763 + t757) + t13 * (t154 * t538 + t536 * t759 + t264) + t53 * (t154 * t510 + t508 * t759 + t536 * t765 + t538 * t99 + t758)) * t689 - m(7) * (g(1) * t908 + g(2) * t909 + g(3) * t899) - m(7) * ((qJD(6) * t467 + t305 * t881 - t430 * t883 + t683) * t866 + (qJD(6) * t465 - t306 * t881 + t430 * t882 + t684) * t865 + (qJD(6) * t541 - t305 * t882 + t306 * t883 + t685) * t867 + 0.2e1 * ((-t43 * t762 + t44 * t761) * t541 + (t38 * t762 - t44 * t750) * t467 + (-t38 * t761 + t43 * t750) * t465) * qJD(5)) / 0.2e1 + (t110 * (-t291 * t542 + t427 * t470) + t111 * (t292 * t542 - t427 * t471) + t83 * (t291 * t471 - t292 * t470)) * t692 + (t898 * t653 + (t468 * t901 + t469 * t900 + t541 * t907) * t430 + (t468 * t912 + t469 * t910 + t541 * t924) * t306 + (t468 * t913 + t469 * t911 + t541 * t925) * t305 + (t465 * t888 + t467 * t887 + t541 * t884) * qJD(5)) * t846 + (t898 * t654 + (t355 * t901 + t356 * t900 + t467 * t907) * t430 + (t355 * t912 + t356 * t910 + t467 * t924) * t306 + (t355 * t913 + t356 * t911 + t467 * t925) * t305 + (t465 * t890 + t467 * t889 + t541 * t885) * qJD(5)) * t854 + (t898 * t655 + (t353 * t901 + t354 * t900 + t465 * t907) * t430 + (t353 * t912 + t354 * t910 + t465 * t924) * t306 + (t353 * t913 + t354 * t911 + t465 * t925) * t305 + (t465 * t892 + t467 * t891 + t541 * t886) * qJD(5)) * t856 - m(5) * (g(1) * t292 + g(2) * t291 + g(3) * t427); -m(7) * (g(1) * t753 + g(2) * t754 + g(3) * t744) + (t11 * (t653 * t762 - t655 * t750) + t43 * (t337 * t750 - t449 * t762 + t653 * t817 - t655 * t764) + t12 * (-t653 * t761 + t654 * t750) + t44 * (-t339 * t750 + t449 * t761 - t653 * t815 + t654 * t764) + t7 * (-t654 * t762 + t655 * t761) + t38 * (-t337 * t761 + t339 * t762 - t654 * t817 + t655 * t815)) * t720 / 0.2e1 + (t30 * (t154 * t653 - t263 * t655) + t58 * (-t130 * t655 - t154 * t449 + t263 * t337 + t653 * t99) + t31 * (-t156 * t653 + t263 * t654) + t59 * (-t101 * t653 + t130 * t654 + t156 * t449 - t263 * t339) + t13 * (-t154 * t654 + t156 * t655) + t53 * (t101 * t655 + t154 * t339 - t156 * t337 - t654 * t99)) * t689 + (t58 * (-t200 * t430 + t304 * t305) + t59 * (t202 * t430 - t304 * t306) + t53 * (t200 * t306 - t202 * t305)) * t690 - m(6) * (g(1) * t202 + g(2) * t200 + g(3) * t304) + (t43 * (t305 * t744 - t430 * t754) + t44 * (-t306 * t744 + t430 * t753) + t38 * (-t305 * t753 + t306 * t754)) * t688 + (-t653 * t886 - t654 * t891 - t655 * t892) * t859 + (-t653 * t885 - t654 * t889 - t655 * t890) * t858 + (-t653 * t884 - t654 * t887 - t655 * t888) * t857 + (t353 * t880 + t354 * t879 - t655 * t878) * t856 + (t892 * t337 + t891 * t339 + t886 * t449 - t653 * t929 - t654 * t937 - t655 * t938) * t855 + (t355 * t880 + t356 * t879 - t654 * t878) * t854 + (t890 * t337 + t889 * t339 + t885 * t449 - t653 * t928 - t654 * t935 - t655 * t936) * t853 + t934 * t337 / 0.2e1 + t933 * t339 / 0.2e1 + (t468 * t880 + t469 * t879 - t653 * t878) * t846 + (t888 * t337 + t887 * t339 + t884 * t449 - t653 * t927 - t654 * t931 - t655 * t932) * t845 + t930 * t449 / 0.2e1 - t941 * t655 / 0.2e1 - t940 * t654 / 0.2e1 - t939 * t653 / 0.2e1; (-(-g(3) + t7) * t653 - (-g(1) + t11) * t654 - (t12 - g(2)) * t655 + (-t306 * t653 + t430 * t654 + t337) * t44 + (t305 * t653 - t430 * t655 + t339) * t43 + (-t305 * t654 + t306 * t655 + t449) * t38) * m(7);];
tau  = t1;
