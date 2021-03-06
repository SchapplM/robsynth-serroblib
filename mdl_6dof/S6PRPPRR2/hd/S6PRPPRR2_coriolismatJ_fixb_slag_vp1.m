% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRPPRR2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_coriolismatJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_coriolismatJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_coriolismatJ_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR2_coriolismatJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPPRR2_coriolismatJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPPRR2_coriolismatJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:17:49
% EndTime: 2019-03-08 19:18:47
% DurationCPUTime: 56.48s
% Computational Cost: add. (196223->1307), mult. (513557->1969), div. (0->0), fcn. (678370->12), ass. (0->608)
t906 = Icges(4,2) + Icges(5,3) - Icges(4,1) - Icges(5,2);
t899 = (Icges(5,6) + Icges(4,4));
t916 = 2 * t899;
t568 = cos(pkin(6));
t566 = sin(pkin(6));
t567 = cos(pkin(10));
t727 = t566 * t567;
t565 = sin(pkin(10));
t728 = t565 * t566;
t772 = sin(pkin(11));
t773 = cos(pkin(11));
t801 = sin(qJ(2));
t802 = cos(qJ(2));
t556 = t801 * t772 - t802 * t773;
t580 = t568 * t556;
t583 = t772 * t802 + t773 * t801;
t505 = -t565 * t583 - t567 * t580;
t860 = t568 * t583;
t506 = -t565 * t556 + t567 * t860;
t901 = Icges(5,5) - Icges(4,6);
t907 = t505 * t906 + t506 * t916 + t727 * t901;
t507 = t556 * t567 + t565 * t860;
t508 = t565 * t580 - t567 * t583;
t908 = t507 * t916 - t508 * t906 + t728 * t901;
t541 = t556 * t566;
t542 = t583 * t566;
t909 = t541 * t906 - 0.2e1 * t899 * t542 + t568 * t901;
t915 = t909 * t568 + t727 * t907 + t728 * t908;
t572 = cos(qJ(5));
t570 = sin(qJ(5));
t726 = t566 * t570;
t461 = -t505 * t572 + t567 * t726;
t725 = t566 * t572;
t462 = t505 * t570 + t567 * t725;
t767 = Icges(6,4) * t462;
t334 = Icges(6,2) * t461 + Icges(6,6) * t506 - t767;
t457 = Icges(6,4) * t461;
t336 = -Icges(6,1) * t462 + Icges(6,5) * t506 + t457;
t613 = Icges(6,5) * t570 + Icges(6,6) * t572;
t912 = Icges(6,3) * t505 + t334 * t572 + t336 * t570 + t506 * t613;
t459 = t508 * t572 + t565 * t726;
t460 = -t508 * t570 + t565 * t725;
t768 = Icges(6,4) * t460;
t333 = -Icges(6,2) * t459 - Icges(6,6) * t507 + t768;
t456 = Icges(6,4) * t459;
t335 = Icges(6,1) * t460 - Icges(6,5) * t507 - t456;
t911 = Icges(6,3) * t508 + t333 * t572 + t335 * t570 - t507 * t613;
t519 = -t541 * t572 + t568 * t570;
t520 = t541 * t570 + t568 * t572;
t766 = Icges(6,4) * t520;
t436 = -Icges(6,2) * t519 + Icges(6,6) * t542 + t766;
t518 = Icges(6,4) * t519;
t437 = Icges(6,1) * t520 + Icges(6,5) * t542 - t518;
t910 = -Icges(6,3) * t541 + t436 * t572 + t437 * t570 + t542 * t613;
t905 = t506 / 0.4e1;
t904 = -t507 / 0.4e1;
t902 = Icges(5,4) - Icges(4,5);
t898 = 2 * Icges(5,6) + 2 * Icges(4,4);
t865 = ((Icges(3,5) * t802 - Icges(3,6) * t801) * t566 + t901 * t542 + t902 * t541) * t568;
t656 = t568 * t802;
t548 = -t565 * t801 + t567 * t656;
t655 = t568 * t801;
t549 = t565 * t802 + t567 * t655;
t867 = Icges(3,5) * t548 - Icges(3,6) * t549 - t505 * t902 + t506 * t901;
t550 = -t565 * t656 - t567 * t801;
t551 = -t565 * t655 + t567 * t802;
t868 = Icges(3,5) * t550 - Icges(3,6) * t551 - t507 * t901 - t508 * t902;
t897 = t727 * t867 - t728 * t868 - t865;
t569 = sin(qJ(6));
t571 = cos(qJ(6));
t376 = -t460 * t569 - t507 * t571;
t377 = t460 * t571 - t507 * t569;
t285 = Icges(7,5) * t377 + Icges(7,6) * t376 + Icges(7,3) * t459;
t765 = Icges(7,4) * t377;
t287 = Icges(7,2) * t376 + Icges(7,6) * t459 + t765;
t374 = Icges(7,4) * t376;
t289 = Icges(7,1) * t377 + Icges(7,5) * t459 + t374;
t191 = t285 * t459 + t287 * t376 + t289 * t377;
t378 = t462 * t569 + t506 * t571;
t379 = -t462 * t571 + t506 * t569;
t286 = Icges(7,5) * t379 + Icges(7,6) * t378 - Icges(7,3) * t461;
t764 = Icges(7,4) * t379;
t288 = Icges(7,2) * t378 - Icges(7,6) * t461 + t764;
t375 = Icges(7,4) * t378;
t290 = Icges(7,1) * t379 - Icges(7,5) * t461 + t375;
t192 = t286 * t459 + t288 * t376 + t290 * t377;
t454 = -t520 * t569 + t542 * t571;
t455 = t520 * t571 + t542 * t569;
t339 = Icges(7,5) * t455 + Icges(7,6) * t454 + Icges(7,3) * t519;
t763 = Icges(7,4) * t455;
t340 = Icges(7,2) * t454 + Icges(7,6) * t519 + t763;
t451 = Icges(7,4) * t454;
t341 = Icges(7,1) * t455 + Icges(7,5) * t519 + t451;
t227 = t339 * t459 + t340 * t376 + t341 * t377;
t896 = t191 * t507 - t192 * t506 - t227 * t542;
t193 = -t285 * t461 + t287 * t378 + t289 * t379;
t194 = -t286 * t461 + t288 * t378 + t290 * t379;
t228 = -t339 * t461 + t340 * t378 + t341 * t379;
t895 = t193 * t507 - t194 * t506 - t228 * t542;
t199 = t285 * t519 + t287 * t454 + t289 * t455;
t200 = t286 * t519 + t288 * t454 + t290 * t455;
t242 = t339 * t519 + t340 * t454 + t341 * t455;
t894 = t199 * t507 - t200 * t506 - t242 * t542;
t889 = -2 * Icges(3,4);
t892 = Icges(3,1) - Icges(3,2);
t888 = 2 * Icges(3,4);
t673 = 0.2e1 * t566;
t815 = -t541 / 0.2e1;
t811 = t542 / 0.4e1;
t887 = t568 / 0.4e1;
t675 = 2 * m(7);
t886 = -t675 / 0.4e1;
t643 = t675 / 0.4e1;
t676 = 2 * m(6);
t885 = t676 / 0.4e1;
t606 = t505 * t565 + t508 * t567;
t884 = pkin(8) * t606;
t331 = Icges(6,5) * t460 - Icges(6,6) * t459 - Icges(6,3) * t507;
t615 = Icges(6,4) * t570 + Icges(6,2) * t572;
t346 = Icges(6,6) * t508 - t507 * t615;
t617 = Icges(6,1) * t570 + Icges(6,4) * t572;
t348 = Icges(6,5) * t508 - t507 * t617;
t179 = t505 * t331 + t461 * t346 - t462 * t348 + t506 * t911;
t332 = -Icges(6,5) * t462 + Icges(6,6) * t461 + Icges(6,3) * t506;
t345 = Icges(6,6) * t505 + t506 * t615;
t347 = Icges(6,5) * t505 + t506 * t617;
t180 = t505 * t332 + t461 * t345 - t462 * t347 + t506 * t912;
t435 = Icges(6,5) * t520 - Icges(6,6) * t519 + Icges(6,3) * t542;
t442 = -Icges(6,6) * t541 + t542 * t615;
t443 = -Icges(6,5) * t541 + t542 * t617;
t218 = t505 * t435 + t461 * t442 - t462 * t443 + t506 * t910;
t724 = t569 * t570;
t409 = t507 * t724 + t508 * t571;
t723 = t570 * t571;
t410 = -t507 * t723 + t508 * t569;
t739 = t507 * t572;
t303 = Icges(7,5) * t410 + Icges(7,6) * t409 + Icges(7,3) * t739;
t305 = Icges(7,4) * t410 + Icges(7,2) * t409 + Icges(7,6) * t739;
t307 = Icges(7,1) * t410 + Icges(7,4) * t409 + Icges(7,5) * t739;
t407 = t505 * t571 - t506 * t724;
t408 = t505 * t569 + t506 * t723;
t745 = t506 * t572;
t130 = -t285 * t745 + t407 * t287 + t408 * t289 - t461 * t303 + t378 * t305 + t379 * t307;
t302 = Icges(7,5) * t408 + Icges(7,6) * t407 - Icges(7,3) * t745;
t304 = Icges(7,4) * t408 + Icges(7,2) * t407 - Icges(7,6) * t745;
t306 = Icges(7,1) * t408 + Icges(7,4) * t407 - Icges(7,5) * t745;
t131 = -t286 * t745 + t407 * t288 + t408 * t290 - t461 * t302 + t378 * t304 + t379 * t306;
t469 = -t541 * t571 - t542 * t724;
t470 = -t541 * t569 + t542 * t723;
t731 = t542 * t572;
t351 = Icges(7,5) * t470 + Icges(7,6) * t469 - Icges(7,3) * t731;
t352 = Icges(7,4) * t470 + Icges(7,2) * t469 - Icges(7,6) * t731;
t353 = Icges(7,1) * t470 + Icges(7,4) * t469 - Icges(7,5) * t731;
t167 = -t339 * t745 + t407 * t340 + t408 * t341 - t461 * t351 + t378 * t352 + t379 * t353;
t46 = t167 * t568 + (t130 * t565 - t131 * t567) * t566;
t883 = t46 + t218 * t568 + (t179 * t565 - t180 * t567) * t566;
t181 = t508 * t331 - t459 * t346 + t460 * t348 - t507 * t911;
t182 = t508 * t332 - t459 * t345 + t460 * t347 - t507 * t912;
t219 = t508 * t435 - t459 * t442 + t460 * t443 - t507 * t910;
t132 = t285 * t739 + t409 * t287 + t410 * t289 + t459 * t303 + t376 * t305 + t377 * t307;
t133 = t286 * t739 + t409 * t288 + t410 * t290 + t459 * t302 + t376 * t304 + t377 * t306;
t168 = t339 * t739 + t409 * t340 + t410 * t341 + t459 * t351 + t376 * t352 + t377 * t353;
t47 = t168 * t568 + (t132 * t565 - t133 * t567) * t566;
t882 = t47 + t219 * t568 + (t181 * t565 - t182 * t567) * t566;
t195 = -t541 * t331 - t519 * t346 + t520 * t348 + t542 * t911;
t196 = -t541 * t332 - t519 * t345 + t520 * t347 + t542 * t912;
t239 = -t541 * t435 - t519 * t442 + t520 * t443 + t542 * t910;
t137 = -t285 * t731 + t469 * t287 + t470 * t289 + t519 * t303 + t454 * t305 + t455 * t307;
t138 = -t286 * t731 + t469 * t288 + t470 * t290 + t519 * t302 + t454 * t304 + t455 * t306;
t190 = -t339 * t731 + t469 * t340 + t470 * t341 + t519 * t351 + t454 * t352 + t455 * t353;
t59 = t190 * t568 + (t137 * t565 - t138 * t567) * t566;
t881 = t239 * t568 + (t195 * t565 - t196 * t567) * t566 + t59;
t618 = rSges(7,1) * t571 - rSges(7,2) * t569;
t327 = rSges(7,3) * t460 - t459 * t618;
t712 = -pkin(5) * t459 + pkin(9) * t460 + t327;
t638 = t712 * t567;
t328 = -rSges(7,3) * t462 + t461 * t618;
t711 = pkin(5) * t461 - pkin(9) * t462 + t328;
t874 = t565 * t711 + t638;
t873 = -t505 * t568 + t541 * t727;
t872 = t541 * t916 + t542 * t906 + t902 * t568;
t871 = t507 * t906 + t898 * t508 - t902 * t728;
t870 = t898 * t505 - t506 * t906 + t902 * t727;
t292 = rSges(7,1) * t379 + rSges(7,2) * t378 - rSges(7,3) * t461;
t715 = -pkin(5) * t462 - pkin(9) * t461 + t292;
t291 = rSges(7,1) * t377 + rSges(7,2) * t376 + rSges(7,3) * t459;
t716 = pkin(5) * t460 + pkin(9) * t459 + t291;
t203 = -t506 * t716 - t507 * t715;
t342 = rSges(7,1) * t455 + rSges(7,2) * t454 + rSges(7,3) * t519;
t704 = pkin(5) * t520 + pkin(9) * t519 + t342;
t225 = t506 * t704 - t542 * t715;
t226 = t507 * t704 + t542 * t716;
t337 = rSges(6,1) * t460 - rSges(6,2) * t459 - rSges(6,3) * t507;
t338 = -rSges(6,1) * t462 + rSges(6,2) * t461 + rSges(6,3) * t506;
t269 = -t337 * t506 - t338 * t507;
t438 = rSges(6,1) * t520 - rSges(6,2) * t519 + rSges(6,3) * t542;
t283 = -t338 * t542 + t438 * t506;
t284 = t337 * t542 + t438 * t507;
t388 = t606 * t566;
t668 = t541 * t728;
t439 = -t508 * t568 - t668;
t722 = (-t203 * t388 - t225 * t873 + t226 * t439) * t643 + (-t269 * t388 - t283 * t873 + t284 * t439) * t885;
t308 = t408 * rSges(7,1) + t407 * rSges(7,2) - rSges(7,3) * t745;
t620 = pkin(5) * t570 - pkin(9) * t572;
t714 = t506 * t620 + t308;
t309 = t410 * rSges(7,1) + t409 * rSges(7,2) + rSges(7,3) * t739;
t713 = -t620 * t507 + t309;
t854 = -t716 * t505 - t713 * t506 + t715 * t508;
t146 = -t507 * t714 + t854;
t354 = t470 * rSges(7,1) + t469 * rSges(7,2) - rSges(7,3) * t731;
t703 = t620 * t542 + t354;
t165 = t505 * t704 + t506 * t703 + t541 * t715 - t542 * t714;
t166 = t507 * t703 - t508 * t704 - t541 * t716 + t542 * t713;
t619 = rSges(6,1) * t570 + rSges(6,2) * t572;
t349 = t505 * rSges(6,3) + t506 * t619;
t350 = t508 * rSges(6,3) - t507 * t619;
t220 = -t337 * t505 + t338 * t508 - t349 * t507 - t350 * t506;
t444 = -t541 * rSges(6,3) + t542 * t619;
t245 = t338 * t541 - t349 * t542 + t438 * t505 + t444 * t506;
t246 = -t337 * t541 + t350 * t542 - t438 * t508 + t444 * t507;
t792 = (t146 * t541 - t165 * t508 - t166 * t505 + t203 * t542 - t225 * t507 + t226 * t506) * t643 + (t220 * t541 - t245 * t508 - t246 * t505 + t269 * t542 - t283 * t507 + t284 * t506) * t885;
t869 = t722 - t792;
t823 = t505 / 0.2e1;
t822 = t506 / 0.2e1;
t819 = t508 / 0.2e1;
t818 = -t507 / 0.2e1;
t812 = t542 / 0.2e1;
t803 = t568 / 0.2e1;
t277 = t435 * t542 - t436 * t519 + t437 * t520;
t862 = t277 * t815;
t843 = m(7) / 0.4e1;
t845 = m(6) / 0.4e1;
t847 = m(5) / 0.4e1;
t634 = t845 + t843 + t847;
t233 = 0.4e1 * t634 * (-t505 * t506 + t507 * t508 + t541 * t542);
t383 = rSges(7,3) * t520 - t519 * t618;
t211 = t291 * t520 + t327 * t519 - t342 * t460 - t383 * t459;
t316 = rSges(7,1) * t376 - rSges(7,2) * t377;
t361 = rSges(7,1) * t454 - rSges(7,2) * t455;
t271 = t316 * t542 + t361 * t507;
t859 = -t211 - t271;
t212 = -t292 * t520 - t328 * t519 - t342 * t462 - t383 * t461;
t317 = rSges(7,1) * t378 - rSges(7,2) * t379;
t270 = -t317 * t542 + t361 * t506;
t858 = t212 + t270;
t751 = t328 * t459;
t752 = t327 * t461;
t756 = t292 * t460;
t758 = t291 * t462;
t189 = t751 + t752 + t756 + t758;
t736 = t507 * t317;
t741 = t506 * t316;
t254 = -t736 - t741;
t857 = t254 + t189;
t688 = Icges(3,6) * t727 + t548 * t892 + t549 * t889;
t687 = -Icges(3,6) * t728 + t550 * t892 + t551 * t889;
t686 = -Icges(3,5) * t727 + t548 * t888 + t549 * t892;
t685 = Icges(3,5) * t728 + t550 * t888 + t551 * t892;
t545 = t548 * pkin(2);
t546 = t550 * pkin(2);
t681 = t545 * t728 + t546 * t727;
t680 = -Icges(3,6) * t568 + (t801 * t889 + t802 * t892) * t566;
t657 = t566 * t801;
t679 = Icges(3,4) * t802 * t673 + Icges(3,5) * t568 + t657 * t892;
t844 = m(7) / 0.2e1;
t856 = m(6) / 0.2e1 + t844;
t855 = m(5) + m(6) + m(7);
t748 = t461 * t309;
t750 = t459 * t308;
t755 = t292 * t507;
t757 = t291 * t506;
t175 = (t755 + t757) * t572 + t748 + t750;
t241 = t874 * t566;
t702 = -pkin(5) * t519 + pkin(9) * t520 + t383;
t637 = t702 * t566;
t258 = -t565 * t637 + t568 * t712;
t259 = -t567 * t637 - t568 * t711;
t368 = -rSges(6,1) * t459 - rSges(6,2) * t460;
t369 = rSges(6,1) * t461 + rSges(6,2) * t462;
t608 = t368 * t567 + t369 * t565;
t298 = t608 * t566;
t448 = -rSges(6,1) * t519 - rSges(6,2) * t520;
t319 = t368 * t568 - t448 * t728;
t320 = -t369 * t568 - t448 * t727;
t721 = (t241 * t541 - t258 * t505 - t259 * t508) * t643 + (t298 * t541 - t319 * t505 - t320 * t508) * t885;
t102 = t191 * t459 - t192 * t461 + t227 * t519;
t103 = t193 * t459 - t194 * t461 + t228 * t519;
t115 = t199 * t459 - t200 * t461 + t242 * t519;
t23 = t130 * t459 - t131 * t461 + t167 * t519 + t572 * t895;
t24 = t132 * t459 - t133 * t461 + t168 * t519 + t572 * t896;
t37 = t137 * t459 - t138 * t461 + t190 * t519 + t572 * t894;
t813 = -t542 / 0.2e1;
t817 = t519 / 0.2e1;
t821 = t507 / 0.2e1;
t825 = -t506 / 0.2e1;
t827 = -t461 / 0.2e1;
t829 = t459 / 0.2e1;
t853 = t572 * (t102 * t821 + t103 * t825 + t115 * t813) + t23 * t827 + t24 * t829 + t37 * t817;
t463 = pkin(4) * t728 - pkin(8) * t507;
t464 = -pkin(4) * t727 + pkin(8) * t506;
t430 = pkin(3) * t506 - qJ(4) * t505;
t434 = -pkin(3) * t507 - qJ(4) * t508;
t646 = pkin(2) * t655 - qJ(3) * t566;
t672 = t802 * pkin(2);
t488 = t565 * t672 + t567 * t646;
t489 = -t565 * t646 + t567 * t672;
t691 = t488 * t728 + t489 * t727;
t631 = t430 * t728 + t434 * t727 + t691;
t600 = t463 * t727 + t464 * t728 + t631;
t178 = (t565 * t715 + t567 * t716) * t566 + t600;
t558 = pkin(2) * t657 + t568 * qJ(3);
t684 = -pkin(3) * t542 - qJ(4) * t541 - t558;
t658 = -pkin(4) * t568 - pkin(8) * t542 + t684;
t595 = (t658 - t704) * t566;
t478 = t568 * t489;
t696 = t568 * t434 + t478;
t660 = t568 * t463 + t696;
t209 = t565 * t595 + t568 * t716 + t660;
t694 = -t430 - t488;
t659 = -t464 + t694;
t210 = (t659 - t715) * t568 + t567 * t595;
t753 = t317 * t565;
t754 = t316 * t567;
t263 = (t753 + t754) * t566;
t278 = t316 * t568 - t361 * t728;
t279 = -t317 * t568 - t361 * t727;
t310 = Icges(7,5) * t376 - Icges(7,6) * t377;
t718 = -Icges(7,2) * t377 + t289 + t374;
t720 = Icges(7,1) * t376 - t287 - t765;
t151 = t310 * t459 + t376 * t718 + t377 * t720;
t311 = Icges(7,5) * t378 - Icges(7,6) * t379;
t717 = -Icges(7,2) * t379 + t290 + t375;
t719 = Icges(7,1) * t378 - t288 - t764;
t152 = t311 * t459 + t376 * t717 + t377 * t719;
t358 = Icges(7,5) * t454 - Icges(7,6) * t455;
t705 = -Icges(7,2) * t455 + t341 + t451;
t706 = Icges(7,1) * t454 - t340 - t763;
t186 = t358 * t459 + t376 * t705 + t377 * t706;
t71 = t186 * t568 + (t151 * t565 - t152 * t567) * t566;
t153 = -t310 * t461 + t378 * t718 + t379 * t720;
t154 = -t311 * t461 + t378 * t717 + t379 * t719;
t187 = -t358 * t461 + t378 * t705 + t379 * t706;
t72 = t187 * t568 + (t153 * t565 - t154 * t567) * t566;
t162 = t310 * t519 + t454 * t718 + t455 * t720;
t163 = t311 * t519 + t454 * t717 + t455 * t719;
t208 = t358 * t519 + t454 * t705 + t455 * t706;
t86 = -t162 * t507 + t163 * t506 + t208 * t542;
t88 = t208 * t568 + (t162 * t565 - t163 * t567) * t566;
t852 = t86 * t887 + t88 * t811 + t71 * t904 + t72 * t905 + (t178 * t254 + t203 * t263 + t209 * t271 + t210 * t270 + t225 * t279 + t226 * t278) * t643;
t110 = t227 * t568 + (t191 * t565 - t192 * t567) * t566;
t111 = t228 * t568 + (t193 * t565 - t194 * t567) * t566;
t123 = t242 * t568 + (t199 * t565 - t200 * t567) * t566;
t229 = t291 * t461 + t292 * t459;
t255 = t291 * t519 - t342 * t459;
t256 = -t292 * t519 - t342 * t461;
t614 = Icges(7,4) * t571 - Icges(7,2) * t569;
t323 = Icges(7,6) * t460 - t459 * t614;
t616 = Icges(7,1) * t571 - Icges(7,4) * t569;
t325 = Icges(7,5) * t460 - t459 * t616;
t612 = Icges(7,5) * t571 - Icges(7,6) * t569;
t598 = Icges(7,3) * t460 + t287 * t569 - t289 * t571 - t459 * t612;
t158 = t285 * t520 + t323 * t454 + t325 * t455 + t519 * t598;
t324 = -Icges(7,6) * t462 + t461 * t614;
t326 = -Icges(7,5) * t462 + t461 * t616;
t597 = -Icges(7,3) * t462 + t288 * t569 - t290 * t571 + t461 * t612;
t159 = t286 * t520 + t324 * t454 + t326 * t455 + t519 * t597;
t381 = Icges(7,6) * t520 - t519 * t614;
t382 = Icges(7,5) * t520 - t519 * t616;
t596 = Icges(7,3) * t520 + t340 * t569 - t341 * t571 - t519 * t612;
t198 = t339 * t520 + t381 * t454 + t382 * t455 + t519 * t596;
t43 = t158 * t459 - t159 * t461 + t198 * t519 + t199 * t460 - t200 * t462 + t242 * t520;
t139 = t285 * t460 + t323 * t376 + t325 * t377 + t459 * t598;
t140 = t286 * t460 + t324 * t376 + t326 * t377 + t459 * t597;
t176 = t339 * t460 + t376 * t381 + t377 * t382 + t459 * t596;
t56 = t176 * t568 + (t139 * t565 - t140 * t567) * t566;
t141 = -t285 * t462 + t323 * t378 + t325 * t379 - t461 * t598;
t142 = -t286 * t462 + t324 * t378 + t326 * t379 - t461 * t597;
t177 = -t339 * t462 + t378 * t381 + t379 * t382 - t461 * t596;
t57 = t177 * t568 + (t141 * t565 - t142 * t567) * t566;
t82 = t198 * t568 + (t158 * t565 - t159 * t567) * t566;
t851 = t520 * t123 / 0.4e1 - t462 * t111 / 0.4e1 + t460 * t110 / 0.4e1 + t43 * t887 + t519 * t82 / 0.4e1 - t461 * t57 / 0.4e1 + t459 * t56 / 0.4e1 + (t178 * t189 + t209 * t211 + t210 * t212 + t229 * t241 + t255 * t258 + t256 * t259) * t643;
t201 = -t519 * t308 - t461 * t354 + (t292 * t542 - t342 * t506) * t572;
t202 = t519 * t309 - t459 * t354 + (-t291 * t542 - t342 * t507) * t572;
t25 = -t130 * t507 + t131 * t506 + t167 * t542 + t193 * t508 + t194 * t505 - t228 * t541;
t26 = -t132 * t507 + t133 * t506 + t168 * t542 + t191 * t508 + t192 * t505 - t227 * t541;
t39 = -t137 * t507 + t138 * t506 + t190 * t542 + t199 * t508 + t200 * t505 - t242 * t541;
t850 = t541 * t115 / 0.4e1 - t508 * t102 / 0.4e1 - t505 * t103 / 0.4e1 - t542 * t37 / 0.4e1 - t519 * t39 / 0.4e1 + t507 * t24 / 0.4e1 - t506 * t23 / 0.4e1 + t461 * t25 / 0.4e1 - t459 * t26 / 0.4e1 + (t146 * t229 + t165 * t256 + t166 * t255 + t175 * t203 + t201 * t225 + t202 * t226) * t886;
t564 = t566 ^ 2;
t849 = 0.2e1 * t568;
t848 = m(5) / 0.2e1;
t837 = (t175 * t541 - t201 * t508 - t202 * t505 + t229 * t542 + t255 * t506 - t256 * t507) * t675;
t836 = t241 / 0.2e1;
t835 = -t258 / 0.2e1;
t834 = t259 / 0.2e1;
t833 = t298 / 0.2e1;
t832 = -t319 / 0.2e1;
t831 = t320 / 0.2e1;
t828 = t460 / 0.2e1;
t826 = -t462 / 0.2e1;
t816 = t520 / 0.2e1;
t809 = -t565 / 0.4e1;
t808 = t565 / 0.2e1;
t806 = -t567 / 0.2e1;
t804 = t567 / 0.4e1;
t797 = (-t229 * t388 + t255 * t439 - t256 * t873) * t675;
t795 = (t263 * t541 - t278 * t505 - t279 * t508) * t675;
t794 = pkin(8) * t508;
t791 = m(7) * qJD(2);
t790 = m(7) * qJD(5);
t789 = m(7) * qJD(6);
t774 = t568 * t88;
t746 = t506 * t567;
t740 = t507 * t565;
t732 = t542 * t568;
t710 = Icges(6,1) * t459 + t333 + t768;
t709 = -Icges(6,1) * t461 + t334 - t767;
t708 = -Icges(6,2) * t460 + t335 - t456;
t707 = Icges(6,2) * t462 + t336 + t457;
t431 = pkin(3) * t508 - qJ(4) * t507;
t535 = t568 * t546;
t697 = t568 * t431 + t535;
t427 = pkin(3) * t505 + qJ(4) * t506;
t695 = -t427 - t545;
t693 = Icges(6,1) * t519 + t436 + t766;
t692 = -Icges(6,2) * t520 + t437 - t518;
t683 = t873 * pkin(8);
t682 = 0.2e1 * t681;
t678 = qJD(2) * t566;
t236 = t331 * t506 + t333 * t461 - t335 * t462;
t237 = t332 * t506 + t334 * t461 - t336 * t462;
t268 = t435 * t506 + t436 * t461 - t437 * t462;
t48 = -t179 * t507 + t180 * t506 + t218 * t542 + t236 * t508 + t237 * t505 - t268 * t541;
t671 = t25 / 0.2e1 + t48 / 0.2e1;
t234 = -t331 * t507 - t333 * t459 + t335 * t460;
t235 = -t332 * t507 - t334 * t459 + t336 * t460;
t267 = -t435 * t507 - t436 * t459 + t437 * t460;
t49 = -t181 * t507 + t182 * t506 + t219 * t542 + t234 * t508 + t235 * t505 - t267 * t541;
t670 = t26 / 0.2e1 + t49 / 0.2e1;
t247 = t331 * t542 - t333 * t519 + t335 * t520;
t248 = t332 * t542 - t334 * t519 + t336 * t520;
t669 = t39 / 0.2e1 + t195 * t818 + t196 * t822 + t239 * t812 + t247 * t819 + t248 * t823 + t862;
t362 = -Icges(6,5) * t459 - Icges(6,6) * t460;
t204 = -t362 * t507 - t459 * t708 - t460 * t710;
t363 = Icges(6,5) * t461 + Icges(6,6) * t462;
t205 = -t363 * t507 - t459 * t707 - t460 * t709;
t445 = -Icges(6,5) * t519 - Icges(6,6) * t520;
t231 = -t445 * t507 - t459 * t692 - t460 * t693;
t124 = t231 * t568 + (t204 * t565 - t205 * t567) * t566;
t666 = t124 / 0.2e1 + t56 / 0.2e1;
t52 = -t139 * t507 + t140 * t506 + t176 * t542;
t665 = t52 / 0.2e1 + t204 * t818 + t205 * t822 + t231 * t812;
t206 = t362 * t506 + t461 * t708 + t462 * t710;
t207 = t363 * t506 + t461 * t707 + t462 * t709;
t232 = t445 * t506 + t461 * t692 + t462 * t693;
t53 = -t141 * t507 + t142 * t506 + t177 * t542;
t664 = t53 / 0.2e1 + t206 * t818 + t207 * t822 + t232 * t812;
t125 = t232 * t568 + (t206 * t565 - t207 * t567) * t566;
t663 = -t57 / 0.2e1 - t125 / 0.2e1;
t213 = t362 * t542 - t519 * t708 - t520 * t710;
t214 = t363 * t542 - t519 * t707 - t520 * t709;
t252 = t445 * t542 - t519 * t692 - t520 * t693;
t662 = t82 / 0.2e1 + t252 * t803 + (t213 * t565 - t214 * t567) * t566 / 0.2e1;
t661 = pkin(8) * t668 + t697;
t653 = t728 / 0.2e1;
t652 = t728 / 0.4e1;
t651 = -t727 / 0.2e1;
t650 = -t727 / 0.4e1;
t649 = t189 / 0.2e1 - t254 / 0.2e1;
t648 = -t211 / 0.2e1 + t271 / 0.2e1;
t647 = t212 / 0.2e1 - t270 / 0.2e1;
t644 = t675 / 0.2e1;
t636 = (-rSges(4,1) * t542 + rSges(4,2) * t541 - rSges(4,3) * t568 - t558) * t566;
t635 = t564 * t672;
t403 = t427 * t728;
t405 = t431 * t727;
t633 = 0.2e1 * t403 + 0.2e1 * t405 + t682;
t632 = t403 + t405 + t681;
t630 = 0.2e1 * t714;
t629 = 0.2e1 * t711;
t628 = t713 + t794;
t627 = (-rSges(5,1) * t568 + rSges(5,2) * t542 - rSges(5,3) * t541 + t684) * t566;
t624 = qJD(6) * t644;
t221 = -t506 * t712 - t507 * t711;
t282 = -t368 * t506 - t369 * t507;
t623 = m(6) * t282 + m(7) * t221;
t243 = t506 * t702 - t542 * t711;
t299 = -t369 * t542 + t448 * t506;
t622 = m(6) * t299 + m(7) * t243;
t244 = t507 * t702 + t542 * t712;
t300 = t368 * t542 + t448 * t507;
t621 = -m(6) * t300 - m(7) * t244;
t18 = ((t836 - t146 / 0.2e1) * t568 + ((t835 + t166 / 0.2e1) * t567 + (t834 - t165 / 0.2e1) * t565) * t566) * m(7) + ((t833 - t220 / 0.2e1) * t568 + ((t832 + t246 / 0.2e1) * t567 + (t831 - t245 / 0.2e1) * t565) * t566) * m(6);
t574 = t220 * t885 + (-t507 * t630 + 0.2e1 * t854) * t843;
t591 = t608 * t676;
t90 = (-t591 / 0.4e1 + t874 * t886) * t566 + t574;
t611 = t90 * qJD(1) - t18 * qJD(3);
t251 = t316 * t461 + t317 * t459;
t149 = (t750 / 0.2e1 + t748 / 0.2e1 + (t757 / 0.2e1 + t755 / 0.2e1) * t572 + (-t754 / 0.2e1 - t753 / 0.2e1) * t566) * m(7);
t78 = ((t263 / 0.2e1 - t175 / 0.2e1) * t568 + ((-t278 / 0.2e1 + t202 / 0.2e1) * t567 + (t279 / 0.2e1 - t201 / 0.2e1) * t565) * t566) * m(7);
t605 = -t149 * qJD(1) + t78 * qJD(3);
t296 = -t388 * t849 + (-t439 * t567 - t565 * t873) * t673;
t184 = t855 * (t296 / 0.4e1 - t732 / 0.2e1 + (t746 / 0.2e1 + t740 / 0.2e1) * t566);
t387 = -0.2e1 * t388;
t573 = 0.2e1 * t855 * t811;
t273 = -t573 + t855 * t387 / 0.4e1;
t604 = -t273 * qJD(1) - t184 * qJD(3);
t603 = (-t438 + t658) * t566;
t602 = t662 - t669;
t65 = -t151 * t507 + t152 * t506 + t186 * t542;
t66 = -t153 * t507 + t154 * t506 + t187 * t542;
t599 = t65 * t652 + t66 * t650 + t852;
t593 = -(-rSges(4,1) * t541 - rSges(4,2) * t542) * t566 - t635;
t155 = (t758 / 0.2e1 + t756 / 0.2e1 + t752 / 0.2e1 + t751 / 0.2e1 + t741 / 0.2e1 + t736 / 0.2e1) * m(7);
t73 = (t505 * t648 - t508 * t647 + t541 * t649) * m(7);
t80 = (t649 * t568 + (t565 * t647 + t567 * t648) * t566) * m(7);
t592 = -t155 * qJD(1) - t80 * qJD(3) - t73 * qJD(4);
t497 = -pkin(3) * t541 + qJ(4) * t542;
t590 = -t635 + (-t444 - t497) * t566;
t589 = -t635 + (-rSges(5,2) * t541 - rSges(5,3) * t542 - t497) * t566;
t587 = t65 * t818 + t66 * t822 + t812 * t86;
t32 = t139 * t459 - t140 * t461 + t176 * t519 + t191 * t460 - t192 * t462 + t227 * t520;
t33 = t141 * t459 - t142 * t461 + t177 * t519 + t193 * t460 - t194 * t462 + t228 * t520;
t586 = t32 * t652 + t33 * t650 + t851;
t584 = -t635 + (-t497 - t703) * t566;
t581 = t895 * t745 / 0.4e1 - t896 * t739 / 0.4e1 + t894 * t731 / 0.4e1 - t850;
t575 = t102 * t828 + t103 * t826 + t115 * t816 + t32 * t829 + t33 * t827 + t43 * t817;
t555 = (rSges(3,1) * t802 - rSges(3,2) * t801) * t566;
t540 = t568 * rSges(3,3) + (rSges(3,1) * t801 + rSges(3,2) * t802) * t566;
t517 = rSges(3,1) * t550 - rSges(3,2) * t551;
t516 = rSges(3,1) * t548 - rSges(3,2) * t549;
t485 = rSges(3,1) * t551 + rSges(3,2) * t550 + rSges(3,3) * t728;
t484 = rSges(3,1) * t549 + rSges(3,2) * t548 - rSges(3,3) * t727;
t433 = -rSges(5,2) * t508 - rSges(5,3) * t507;
t432 = rSges(4,1) * t508 + rSges(4,2) * t507;
t429 = -rSges(5,2) * t505 + rSges(5,3) * t506;
t428 = rSges(4,1) * t505 - rSges(4,2) * t506;
t400 = -rSges(4,1) * t507 + rSges(4,2) * t508 + rSges(4,3) * t728;
t399 = rSges(4,1) * t506 + rSges(4,2) * t505 - rSges(4,3) * t727;
t398 = -rSges(5,1) * t727 - rSges(5,2) * t506 - rSges(5,3) * t505;
t397 = rSges(5,1) * t728 + rSges(5,2) * t507 - rSges(5,3) * t508;
t330 = (-t428 - t545) * t568 + t593 * t567;
t329 = t568 * t432 + t565 * t593 + t535;
t301 = (t428 * t565 + t432 * t567) * t566 + t681;
t281 = (-t429 + t695) * t568 + t589 * t567;
t280 = t568 * t433 + t565 * t589 + t697;
t275 = (-t398 + t694) * t568 + t567 * t627;
t274 = t397 * t568 + t565 * t627 + t696;
t272 = t387 * t634 + t573;
t266 = -t317 * t519 - t361 * t461;
t265 = t316 * t519 - t361 * t459;
t264 = (t429 * t565 + t433 * t567) * t566 + t632;
t261 = (-t349 + t695) * t568 + t590 * t567 + t683;
t260 = (t350 + t794) * t568 + t590 * t565 + t661;
t257 = (t397 * t567 + t398 * t565) * t566 + t631;
t250 = (-t338 + t659) * t568 + t567 * t603;
t249 = t337 * t568 + t565 * t603 + t660;
t238 = (t349 * t565 + t350 * t567 + t884) * t566 + t632;
t222 = (t337 * t567 + t338 * t565) * t566 + t600;
t216 = (t695 - t714) * t568 + t584 * t567 + t683;
t215 = t565 * t584 + t568 * t628 + t661;
t197 = (t628 * t567 + (pkin(8) * t505 + t714) * t565) * t566 + t632;
t183 = (t296 + 0.2e1 * t732 + (-t740 - t746) * t673) * t634;
t172 = t795 / 0.4e1;
t161 = -0.4e1 * t257 * t388 + 0.4e1 * t274 * t439 - 0.4e1 * t275 * t873;
t160 = -t247 * t507 + t248 * t506 + t277 * t542;
t156 = t857 * t643;
t150 = 0.2e1 * t175 * t843 + t263 * t643;
t144 = -t236 * t507 + t237 * t506 + t268 * t542;
t143 = -t234 * t507 + t235 * t506 + t267 * t542;
t135 = t797 / 0.4e1;
t129 = -0.4e1 * t222 * t388 + 0.4e1 * t249 * t439 - 0.4e1 * t250 * t873;
t127 = -t213 * t507 + t214 * t506 + t252 * t542;
t126 = 0.4e1 * t222 * t298 + 0.4e1 * t249 * t319 + 0.4e1 * t250 * t320;
t109 = 0.4e1 * t220 * t269 + 0.4e1 * t245 * t283 + 0.4e1 * t246 * t284;
t92 = -0.4e1 * t178 * t388 + 0.4e1 * t209 * t439 - 0.4e1 * t210 * t873;
t91 = t574 + (t591 + m(7) * (t565 * t629 + 0.2e1 * t638)) * t566 / 0.4e1;
t89 = 0.4e1 * t203 * t254 + 0.4e1 * t225 * t270 + 0.4e1 * t226 * t271;
t84 = t162 * t459 - t163 * t461 + t208 * t519;
t83 = 0.4e1 * t178 * t263 + 0.4e1 * t209 * t278 + 0.4e1 * t210 * t279;
t81 = (t857 * t849 + (t565 * t858 + t567 * t859) * t673) * t843;
t77 = ((t263 + t175) * t849 + ((-t202 - t278) * t567 + (t201 + t279) * t565) * t673) * t843;
t76 = -t158 * t507 + t159 * t506 + t198 * t542;
t74 = (t505 * t859 - t508 * t858 + t541 * t857) * t643;
t67 = 0.4e1 * t189 * t229 + 0.4e1 * t211 * t255 + 0.4e1 * t212 * t256;
t62 = t153 * t459 - t154 * t461 + t187 * t519;
t61 = t151 * t459 - t152 * t461 + t186 * t519;
t60 = 0.4e1 * t178 * t241 + 0.4e1 * t209 * t258 + 0.4e1 * t210 * t259;
t58 = 0.4e1 * t175 * t229 + 0.4e1 * t201 * t256 + 0.4e1 * t202 * t255;
t50 = t837 / 0.4e1;
t45 = t129 * t845 + t161 * t847 + t843 * t92;
t44 = 0.4e1 * t146 * t203 + 0.4e1 * t165 * t225 + 0.4e1 * t166 * t226;
t29 = t135 + t50 - t795 / 0.4e1;
t28 = t172 + t135 - t837 / 0.4e1;
t27 = t172 + t50 - t797 / 0.4e1;
t17 = (t220 * t849 + (t245 * t565 - t246 * t567) * t673) * t845 + (t146 * t849 + (t165 * t565 - t166 * t567) * t673) * t843 + (m(6) * t833 + m(7) * t836) * t568 + ((m(6) * t832 + m(7) * t835) * t567 + (m(6) * t831 + m(7) * t834) * t565) * t566;
t12 = t721 - t869;
t11 = t721 + t869;
t10 = t722 + t792 - t721;
t9 = t843 * t89 + t587;
t8 = t83 * t843 + t774 / 0.2e1 + (t71 * t808 + t72 * t806) * t566;
t7 = t60 * t843 + t126 * t845 + t662 * t568 + (t565 * t666 + t567 * t663) * t566;
t6 = t67 * t843 + t575;
t5 = t58 * t843 + t853;
t4 = t44 * t843 + t109 * t845 + t669 * t542 + (t894 / 0.2e1 - t160 / 0.2e1) * t541 - t670 * t507 + (-t896 / 0.2e1 + t143 / 0.2e1) * t508 + t671 * t506 + (-t895 / 0.2e1 + t144 / 0.2e1) * t505;
t3 = t581 + t586 + (t65 * t809 + t66 * t804) * t566 - t852;
t2 = t581 + (t32 * t809 + t33 * t804) * t566 + t599 - t851;
t1 = t586 + (-t811 * t894 - t895 * t905 - t896 * t904) * t572 + t599 + t850;
t13 = [0 (m(4) * t682 / 0.2e1 + t633 * t848 + t856 * (t673 * t884 + t633)) * qJD(2) + t272 * qJD(4) + t91 * qJD(5) + t150 * qJD(6) + ((m(3) * t517 + m(4) * t432 + m(5) * t433 + m(6) * t350 + t644 * t713) * t567 + (m(3) * t516 + m(4) * t428 + m(5) * t429 + m(6) * t349 + t630 * t844) * t565) * t678, 0, t272 * qJD(2), t91 * qJD(2) + (-(m(6) * t369 + t629 * t844) * t507 + (-m(6) * t368 - t644 * t712) * t506) * qJD(5) + t156 * qJD(6), t150 * qJD(2) + t156 * qJD(5) + t251 * t624; qJD(4) * t273 - qJD(5) * t90 - qJD(6) * t149, t45 * qJD(4) + t7 * qJD(5) + t8 * qJD(6) + ((t178 * t197 + t209 * t215 + t210 * t216) * m(7) + (t222 * t238 + t249 * t260 + t250 * t261) * m(6) + m(3) * ((-t484 * t568 - t540 * t727) * (-t516 * t568 - t555 * t727) + (t485 * t568 - t540 * t728) * (t517 * t568 - t555 * t728) + (t484 * t565 + t485 * t567) * t564 * (t516 * t565 + t517 * t567)) + m(5) * (t257 * t264 + t274 * t280 + t275 * t281) + m(4) * (((-t399 - t488) * t568 + t567 * t636) * t330 + (t400 * t568 + t565 * t636 + t478) * t329 + ((t399 * t565 + t400 * t567) * t566 + t691) * t301) + ((((t565 * t685 - t567 * t686) * t802 + (t565 * t687 - t567 * t688) * t801) * t566 + (t565 * t908 + t567 * t907) * t542 + (-t565 * t871 + t567 * t870) * t541) * t566 + (t865 + t909 * t542 + t872 * t541 + (t565 * t868 - t567 * t867 + t679 * t802 + t680 * t801) * t566) * t568 + t881) * t803 + ((-t508 * t870 - t550 * t686 - t551 * t688) * t727 + (-t508 * t872 + t550 * t679 + t551 * t680) * t568 + (t508 * t871 + t550 * t685 + t551 * t687 - t897) * t728 + t882 - t915 * t507) * t653 + ((t505 * t871 + t548 * t685 + t549 * t687) * t728 + (-t505 * t872 + t548 * t679 + t549 * t680) * t568 + (-t505 * t870 - t548 * t686 - t549 * t688 + t897) * t727 + t883 + t915 * t506) * t651) * qJD(2), qJD(4) * t184 + qJD(5) * t18 + qJD(6) * t78, t45 * qJD(2) + (-t233 + 0.2e1 * (t848 + t856) * (-t388 * t541 - t439 * t505 + t508 * t873)) * qJD(4) + t11 * qJD(5) + t28 * qJD(6) - t604, t7 * qJD(2) + t11 * qJD(4) + t1 * qJD(6) + (-t44 / 0.4e1 + t178 * t221 + t203 * t241 + t209 * t244 + t210 * t243 + t225 * t259 + t226 * t258) * t790 - t611 + ((t565 * t665 - t567 * t664) * t566 + t602 * t542 - (t666 - t670) * t507 + (-t663 - t671) * t506 + (t222 * t282 + t249 * t300 + t250 * t299 + t269 * t298 + t283 * t320 + t284 * t319 - t109 / 0.4e1) * m(6) - (-t895 + t144) * t505 / 0.2e1 - (-t896 + t143) * t508 / 0.2e1 + (-t894 + t160) * t541 / 0.2e1 + (t127 + t76) * t803) * qJD(5), t8 * qJD(2) + t28 * qJD(4) + t1 * qJD(5) + (t178 * t251 + t209 * t265 + t210 * t266 + t229 * t263 + t255 * t278 + t256 * t279 - t58 / 0.4e1) * t789 + t605 + (t61 * t653 + t62 * t651 + t71 * t829 + t72 * t827 + t803 * t84 + t817 * t88 - t853) * qJD(6); 0 ((m(4) * t301 + m(5) * t264 + m(6) * t238 + m(7) * t197) * t568 + ((-m(4) * t329 - m(5) * t280 - m(6) * t260 - m(7) * t215) * t567 + (m(4) * t330 + m(5) * t281 + m(6) * t261 + m(7) * t216) * t565) * t566) * qJD(2) + t183 * qJD(4) + t17 * qJD(5) + t77 * qJD(6), 0, t183 * qJD(2), t17 * qJD(2) + (t623 * t568 + (t565 * t622 + t567 * t621) * t566) * qJD(5) + t81 * qJD(6), t77 * qJD(2) + t81 * qJD(5) + (t251 * t849 + (-t265 * t567 + t266 * t565) * t673) * t789 / 0.2e1; -t273 * qJD(2) ((t178 * t542 + t197 * t541 + t209 * t506 - t210 * t507 - t215 * t505 - t216 * t508 - t92 / 0.4e1) * m(7) + (t222 * t542 + t238 * t541 + t249 * t506 - t250 * t507 - t260 * t505 - t261 * t508 - t129 / 0.4e1) * m(6) + (t257 * t542 + t264 * t541 + t274 * t506 - t275 * t507 - t280 * t505 - t281 * t508 - t161 / 0.4e1) * m(5)) * qJD(2) + t233 * qJD(4) + t12 * qJD(5) + t27 * qJD(6) + t604, -t184 * qJD(2), t233 * qJD(2), t12 * qJD(2) + (t505 * t621 - t508 * t622 + t541 * t623) * qJD(5) + t74 * qJD(6), t27 * qJD(2) + t74 * qJD(5) + (t251 * t541 - t265 * t505 - t266 * t508) * t624; qJD(2) * t90 - qJD(6) * t155, t10 * qJD(4) + t4 * qJD(5) + t2 * qJD(6) + ((t247 * t565 - t248 * t567) * t815 + (t234 * t565 - t235 * t567) * t819 + (t236 * t565 - t237 * t567) * t823 - (t124 + t56) * t565 / 0.2e1 + (t49 + t26) * t808 + (t48 + t25) * t806 + (t57 + t125) * t567 / 0.2e1) * t678 + (-t60 / 0.4e1 + t146 * t178 + t165 * t210 + t166 * t209 + t197 * t203 + t215 * t226 + t216 * t225) * t791 + t611 + (t110 * t819 + t111 * t823 + t123 * t815 + (t267 * t819 + t268 * t823 - t602 + t862) * t568 + (t220 * t222 + t238 * t269 + t245 * t250 + t246 * t249 + t260 * t284 + t261 * t283 - t126 / 0.4e1) * m(6) + t883 * t822 + t882 * t818 + t881 * t812) * qJD(2), -qJD(2) * t18 - qJD(6) * t80, qJD(2) * t10 - qJD(6) * t73, t4 * qJD(2) + ((t203 * t221 + t225 * t243 + t226 * t244) * m(7) + (t269 * t282 + t283 * t299 + t284 * t300) * m(6) + (t76 / 0.2e1 + t127 / 0.2e1) * t542 - t665 * t507 + t664 * t506) * qJD(5) + t9 * qJD(6), t2 * qJD(2) + t9 * qJD(5) + (t61 * t818 + t62 * t822 + t65 * t829 + t66 * t827 + t812 * t84 + t817 * t86 - t575) * qJD(6) + (t203 * t251 + t225 * t266 + t226 * t265 + t229 * t254 + t255 * t271 + t256 * t270 - t67 / 0.4e1) * t789 + t592; qJD(2) * t149 + qJD(5) * t155, t29 * qJD(4) + t3 * qJD(5) + t5 * qJD(6) + ((-t23 / 0.2e1 + t72 / 0.2e1) * t567 + (t24 / 0.2e1 - t71 / 0.2e1) * t565) * t678 + (t175 * t178 + t197 * t229 + t201 * t210 + t202 * t209 + t215 * t255 + t216 * t256 - t83 / 0.4e1) * t791 - t605 + (t46 * t827 + t47 * t829 + t59 * t817 + t37 * t803 - t774 / 0.2e1 + (t110 * t821 + t111 * t825 + t123 * t813) * t572) * qJD(2), -qJD(2) * t78 + qJD(5) * t80, qJD(2) * t29 + qJD(5) * t73, t3 * qJD(2) + (t32 * t818 + t33 * t822 + t43 * t812 + t52 * t829 + t53 * t827 + t76 * t817 - t816 * t894 - t826 * t895 - t828 * t896 - t587) * qJD(5) + t6 * qJD(6) + (t189 * t203 + t211 * t226 + t212 * t225 + t221 * t229 + t243 * t256 + t244 * t255 - t89 / 0.4e1) * t790 - t592, t5 * qJD(2) + t6 * qJD(5) + ((t229 * t251 + t255 * t265 + t256 * t266) * m(7) + t61 * t829 + t62 * t827 + t84 * t817) * qJD(6);];
Cq  = t13;
