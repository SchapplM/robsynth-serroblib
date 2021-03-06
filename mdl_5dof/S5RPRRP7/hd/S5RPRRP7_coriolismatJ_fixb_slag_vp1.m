% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRP7_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP7_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP7_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP7_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP7_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:44:35
% EndTime: 2019-12-31 18:45:16
% DurationCPUTime: 34.06s
% Computational Cost: add. (67036->801), mult. (82675->1129), div. (0->0), fcn. (91490->8), ass. (0->468)
t510 = qJ(1) + pkin(8);
t507 = sin(t510);
t508 = cos(t510);
t512 = sin(qJ(4));
t515 = cos(qJ(4));
t516 = cos(qJ(3));
t697 = t515 * t516;
t435 = t507 * t697 - t508 * t512;
t421 = Icges(6,5) * t435;
t692 = t516 * t512;
t434 = t507 * t692 + t508 * t515;
t513 = sin(qJ(3));
t710 = t507 * t513;
t321 = -Icges(6,6) * t710 - Icges(6,3) * t434 - t421;
t323 = Icges(5,5) * t435 - Icges(5,6) * t434 + Icges(5,3) * t710;
t326 = Icges(6,4) * t435 + Icges(6,2) * t710 + Icges(6,6) * t434;
t424 = Icges(5,4) * t435;
t329 = -Icges(5,2) * t434 + Icges(5,6) * t710 + t424;
t420 = Icges(6,5) * t434;
t332 = Icges(6,1) * t435 + Icges(6,4) * t710 + t420;
t423 = Icges(5,4) * t434;
t336 = -Icges(5,1) * t435 - Icges(5,5) * t710 + t423;
t889 = (t323 + t326) * t710 + (t332 - t336) * t435 + (-t321 - t329) * t434;
t437 = t507 * t512 + t508 * t697;
t422 = Icges(6,5) * t437;
t436 = -t507 * t515 + t508 * t692;
t706 = t508 * t513;
t322 = Icges(6,6) * t706 + Icges(6,3) * t436 + t422;
t325 = Icges(5,5) * t437 - Icges(5,6) * t436 + Icges(5,3) * t706;
t328 = Icges(6,4) * t437 + Icges(6,2) * t706 + Icges(6,6) * t436;
t731 = Icges(5,4) * t437;
t331 = -Icges(5,2) * t436 + Icges(5,6) * t706 + t731;
t726 = Icges(6,5) * t436;
t334 = Icges(6,1) * t437 + Icges(6,4) * t706 + t726;
t425 = Icges(5,4) * t436;
t337 = Icges(5,1) * t437 + Icges(5,5) * t706 - t425;
t888 = (t325 + t328) * t710 + (t334 + t337) * t435 + (t322 - t331) * t434;
t698 = t513 * t515;
t501 = Icges(6,5) * t698;
t703 = t512 * t513;
t721 = Icges(6,6) * t516;
t438 = Icges(6,3) * t703 + t501 - t721;
t725 = Icges(6,5) * t512;
t584 = Icges(6,1) * t515 + t725;
t446 = -Icges(6,4) * t516 + t513 * t584;
t582 = Icges(6,4) * t515 + Icges(6,6) * t512;
t442 = -Icges(6,2) * t516 + t513 * t582;
t699 = t513 * t442;
t254 = t434 * t438 + t435 * t446 + t507 * t699;
t729 = Icges(5,4) * t515;
t583 = -Icges(5,2) * t512 + t729;
t444 = -Icges(5,6) * t516 + t513 * t583;
t730 = Icges(5,4) * t512;
t585 = Icges(5,1) * t515 - t730;
t448 = -Icges(5,5) * t516 + t513 * t585;
t580 = Icges(5,5) * t515 - Icges(5,6) * t512;
t440 = -Icges(5,3) * t516 + t513 * t580;
t700 = t513 * t440;
t255 = -t434 * t444 + t435 * t448 + t507 * t700;
t831 = t254 + t255;
t887 = t889 * t507 + t888 * t508;
t760 = t513 / 0.2e1;
t759 = -t516 / 0.2e1;
t734 = rSges(6,3) + qJ(5);
t757 = rSges(6,1) + pkin(4);
t886 = t512 * t734 + t515 * t757;
t825 = t886 * t513;
t766 = t507 / 0.2e1;
t764 = -t508 / 0.2e1;
t622 = t831 * t759 + t887 * t760;
t745 = cos(qJ(1)) * pkin(1);
t592 = t507 * pkin(6) + t745;
t744 = pkin(3) * t516;
t618 = pkin(2) + t744;
t756 = rSges(6,2) + pkin(7);
t802 = t513 * t756 + t618;
t854 = t734 * t436 + t757 * t437;
t241 = t508 * t802 + t592 + t854;
t597 = -sin(qJ(1)) * pkin(1) + t508 * pkin(6);
t812 = -t434 * t734 - t435 * t757;
t866 = -t507 * t802 + t597 + t812;
t871 = -t241 * t436 + t434 * t866;
t885 = t871 * m(6) * qJD(1);
t258 = t436 * t438 + t437 * t446 + t508 * t699;
t187 = -t321 * t436 + t326 * t706 + t437 * t332;
t188 = t436 * t322 + t328 * t706 + t437 * t334;
t576 = t187 * t507 + t188 * t508;
t873 = -t516 * t258 + t513 * t576;
t259 = -t436 * t444 + t437 * t448 + t508 * t700;
t189 = t323 * t706 - t436 * t329 - t336 * t437;
t190 = t325 * t706 - t436 * t331 + t437 * t337;
t575 = t189 * t507 + t190 * t508;
t874 = -t516 * t259 + t513 * t575;
t621 = t873 / 0.2e1 + t874 / 0.2e1;
t351 = -Icges(5,5) * t434 - Icges(5,6) * t435;
t353 = -Icges(6,4) * t434 + Icges(6,6) * t435;
t881 = t351 + t353;
t352 = -Icges(5,5) * t436 - Icges(5,6) * t437;
t354 = -Icges(6,4) * t436 + Icges(6,6) * t437;
t880 = t352 + t354;
t738 = rSges(5,1) * t515;
t589 = -rSges(5,2) * t512 + t738;
t543 = t589 * t513;
t459 = -rSges(5,3) * t516 + t543;
t810 = -t435 * rSges(5,1) + t434 * rSges(5,2);
t339 = rSges(5,3) * t710 - t810;
t714 = t339 * t516;
t273 = t459 * t710 + t714;
t666 = Icges(5,2) * t437 - t337 + t425;
t668 = Icges(6,3) * t437 - t334 - t726;
t878 = t666 + t668;
t667 = Icges(5,2) * t435 + t336 + t423;
t669 = Icges(6,3) * t435 - t332 - t420;
t877 = t667 + t669;
t670 = -Icges(5,1) * t436 - t331 - t731;
t672 = -Icges(6,1) * t436 + t322 + t422;
t876 = t670 + t672;
t671 = -Icges(5,1) * t434 - t329 - t424;
t673 = -Icges(6,1) * t434 - t321 + t421;
t875 = t671 + t673;
t716 = t326 * t516;
t853 = t321 * t512 - t332 * t515;
t213 = t853 * t513 + t716;
t718 = t323 * t516;
t852 = t329 * t512 + t336 * t515;
t216 = t852 * t513 + t718;
t868 = -t213 - t216;
t569 = t322 * t512 + t334 * t515;
t715 = t328 * t516;
t214 = t513 * t569 - t715;
t567 = -t331 * t512 + t337 * t515;
t717 = t325 * t516;
t217 = t513 * t567 - t717;
t867 = t214 + t217;
t826 = rSges(6,2) * t710 - t812;
t596 = t826 * t516;
t856 = rSges(6,2) * t516 - t825;
t865 = -t710 * t856 + t596;
t864 = -t187 * t508 + t188 * t507;
t863 = -t189 * t508 + t190 * t507;
t862 = t877 * t434 + t875 * t435 + t881 * t710;
t861 = t878 * t434 + t876 * t435 + t880 * t710;
t860 = t877 * t436 + t875 * t437 + t881 * t706;
t859 = t878 * t436 + t876 * t437 + t880 * t706;
t467 = (-Icges(6,4) * t512 + Icges(6,6) * t515) * t513;
t465 = (Icges(6,3) * t515 - t725) * t513;
t648 = -t446 + t465;
t469 = -Icges(6,1) * t703 + t501;
t650 = t438 + t469;
t218 = t434 * t648 + t435 * t650 + t467 * t710;
t466 = (-Icges(5,5) * t512 - Icges(5,6) * t515) * t513;
t468 = (-Icges(5,2) * t515 - t730) * t513;
t647 = -t448 - t468;
t470 = (-Icges(5,1) * t512 - t729) * t513;
t649 = -t444 + t470;
t219 = t434 * t647 + t435 * t649 + t466 * t710;
t858 = -t218 - t219;
t220 = t436 * t648 + t437 * t650 + t467 * t706;
t221 = t436 * t647 + t437 * t649 + t466 * t706;
t857 = -t220 - t221;
t795 = t508 ^ 2;
t796 = t507 ^ 2;
t855 = t795 + t796;
t271 = t434 * t757 - t435 * t734;
t851 = -m(5) / 0.2e1;
t850 = -m(6) / 0.2e1;
t579 = Icges(6,5) * t515 + Icges(6,3) * t512;
t531 = -t513 * t579 + t721;
t383 = t531 * t507;
t391 = t446 * t507;
t553 = t442 * t507 - t853;
t158 = t553 * t516 + (t383 * t512 - t391 * t515 + t326) * t513;
t389 = t444 * t507;
t393 = t448 * t507;
t551 = -t440 * t507 + t852;
t160 = -t551 * t516 + (t389 * t512 - t393 * t515 + t323) * t513;
t838 = t158 + t160;
t384 = t531 * t508;
t392 = t446 * t508;
t552 = t442 * t508 + t569;
t159 = t552 * t516 + (t384 * t512 - t392 * t515 + t328) * t513;
t390 = t444 * t508;
t394 = t448 * t508;
t550 = -t440 * t508 - t567;
t161 = -t550 * t516 + (t390 * t512 - t394 * t515 + t325) * t513;
t837 = t159 + t161;
t170 = -t353 * t516 + (t512 * t669 + t515 * t673) * t513;
t172 = -t351 * t516 + (t512 * t667 + t515 * t671) * t513;
t836 = t170 + t172;
t171 = -t354 * t516 + (t512 * t668 + t515 * t672) * t513;
t173 = -t352 * t516 + (t512 * t666 + t515 * t670) * t513;
t835 = t171 + t173;
t439 = Icges(6,6) * t513 + t516 * t579;
t447 = Icges(6,4) * t513 + t516 * t584;
t443 = Icges(6,2) * t513 + t516 * t582;
t565 = t512 * t438 + t515 * t446;
t549 = -t443 + t565;
t711 = t442 * t516;
t523 = -t513 * t549 + t711;
t199 = t434 * t439 + t435 * t447 + t507 * t523;
t445 = Icges(5,6) * t513 + t516 * t583;
t449 = Icges(5,5) * t513 + t516 * t585;
t441 = Icges(5,3) * t513 + t516 * t580;
t564 = -t512 * t444 + t515 * t448;
t548 = t441 - t564;
t712 = t440 * t516;
t524 = t513 * t548 + t712;
t200 = -t434 * t445 + t435 * t449 + t507 * t524;
t834 = t199 + t200;
t201 = t436 * t439 + t437 * t447 + t508 * t523;
t202 = -t436 * t445 + t437 * t449 + t508 * t524;
t833 = t201 + t202;
t234 = t549 * t516 + (t512 * t439 + t515 * t447 + t442) * t513;
t235 = -t548 * t516 + (-t512 * t445 + t515 * t449 + t440) * t513;
t832 = t234 + t235;
t830 = t258 + t259;
t732 = Icges(4,4) * t513;
t481 = Icges(4,1) * t516 - t732;
t408 = Icges(4,5) * t507 + t481 * t508;
t693 = t516 * t408;
t380 = t507 * t693;
t477 = Icges(4,5) * t516 - Icges(4,6) * t513;
t404 = Icges(4,3) * t507 + t477 * t508;
t595 = t508 * t404 - t380;
t709 = t507 * t516;
t403 = Icges(4,5) * t709 - Icges(4,6) * t710 - Icges(4,3) * t508;
t488 = Icges(4,4) * t710;
t407 = Icges(4,1) * t709 - Icges(4,5) * t508 - t488;
t694 = t516 * t407;
t660 = -t507 * t403 - t508 * t694;
t509 = Icges(4,4) * t516;
t724 = Icges(4,2) * t513;
t406 = Icges(4,6) * t507 + (t509 - t724) * t508;
t701 = t513 * t406;
t405 = Icges(4,4) * t709 - Icges(4,2) * t710 - Icges(4,6) * t508;
t702 = t513 * t405;
t829 = -t507 * t701 - t508 * t702 - t595 - t660;
t283 = t513 * t565 - t711;
t284 = t513 * t564 - t712;
t828 = t283 + t284;
t755 = rSges(5,3) + pkin(7);
t803 = t513 * t755 + t618;
t269 = -t507 * t803 + t597 + t810;
t827 = t868 * t507 + t867 * t508;
t824 = -t507 / 0.2e1;
t762 = t508 / 0.2e1;
t823 = -t513 / 0.2e1;
t822 = t516 / 0.2e1;
t820 = t858 * t516 + (t862 * t507 + t861 * t508) * t513;
t819 = t857 * t516 + (t860 * t507 + t859 * t508) * t513;
t662 = rSges(6,2) * t706 + t854;
t814 = t516 * t662;
t347 = t434 * t507 + t436 * t508;
t299 = (-t347 + t692) * t703;
t611 = t692 / 0.2e1;
t303 = (t611 - t347 / 0.2e1) * m(6);
t740 = m(6) * qJD(5);
t809 = t303 * qJD(2) + t299 * t740;
t808 = t158 / 0.2e1 + t160 / 0.2e1;
t807 = -t159 / 0.2e1 - t161 / 0.2e1;
t590 = t437 * rSges(5,1) - t436 * rSges(5,2);
t270 = t508 * t803 + t590 + t592;
t272 = -t757 * t436 + t734 * t437;
t493 = pkin(3) * t513 - pkin(7) * t516;
t620 = t493 - t856;
t314 = t620 * t507;
t316 = t620 * t508;
t363 = -rSges(5,1) * t434 - rSges(5,2) * t435;
t368 = -rSges(5,1) * t436 - rSges(5,2) * t437;
t642 = (-t757 * t512 + t734 * t515) * t513;
t372 = t642 * t507;
t373 = t642 * t508;
t645 = t459 + t493;
t374 = t645 * t507;
t376 = t645 * t508;
t473 = (-rSges(5,1) * t512 - rSges(5,2) * t515) * t513;
t806 = (t363 * t376 - t368 * t374 + (-t269 * t508 - t270 * t507) * t473) * t851 + (-t241 * t372 - t271 * t316 - t272 * t314 - t373 * t866) * t850;
t658 = t856 * t507;
t644 = rSges(6,2) * t513 + t516 * t886;
t799 = t513 * t644 - t516 * t856;
t168 = t507 * t799 - t513 * t826 + t658 * t516;
t705 = t508 * t516;
t492 = rSges(6,2) * t705;
t657 = -t825 * t508 + t492;
t169 = -t508 * t799 + t662 * t513 - t657 * t516;
t229 = -t706 * t856 + t814;
t396 = t459 * t507;
t461 = rSges(5,3) * t513 + t516 * t589;
t563 = t459 * t516 + t461 * t513;
t236 = -t339 * t513 - t396 * t516 + t507 * t563;
t343 = rSges(5,3) * t706 + t590;
t641 = t508 * rSges(5,2) * t703 + rSges(5,3) * t705;
t398 = -rSges(5,1) * t508 * t698 + t641;
t237 = t343 * t513 - t398 * t516 - t508 * t563;
t713 = t343 * t516;
t275 = t459 * t706 + t713;
t496 = pkin(7) * t705;
t281 = t492 + t496 + (-pkin(3) - t886) * t706;
t495 = pkin(3) * t710;
t282 = t495 + (-t756 * t516 + t825) * t507;
t345 = t495 + (-t516 * t755 + t543) * t507;
t346 = t496 + (-pkin(3) - t738) * t706 + t641;
t805 = (t168 * t866 + t169 * t241 - t229 * t281 + t282 * t865) * t850 + (t236 * t269 + t237 * t270 + t273 * t345 - t275 * t346) * t851;
t601 = t448 / 0.2e1 + t446 / 0.2e1;
t603 = t444 / 0.2e1 - t438 / 0.2e1;
t801 = -t512 * (t468 / 0.2e1 - t465 / 0.2e1 + t601) + t515 * (t470 / 0.2e1 + t469 / 0.2e1 - t603);
t733 = Icges(4,1) * t513;
t798 = t512 * t603 - t515 * t601 + t441 / 0.2e1 + t443 / 0.2e1 - t509 + t724 / 0.2e1 - t733 / 0.2e1;
t478 = Icges(4,2) * t516 + t732;
t797 = t512 * (t445 / 0.2e1 - t439 / 0.2e1) - t515 * (t449 / 0.2e1 + t447 / 0.2e1) - t440 / 0.2e1 - t442 / 0.2e1 + t478 / 0.2e1 - t481 / 0.2e1;
t793 = 0.4e1 * qJD(1);
t792 = 2 * qJD(3);
t790 = 2 * qJD(4);
t789 = 4 * qJD(4);
t788 = m(5) / 0.2e1;
t786 = m(6) / 0.2e1;
t192 = (-t396 * t513 + t714) * t508 + (-t398 * t513 - t713) * t507;
t247 = (t339 * t508 - t343 * t507) * t513;
t784 = m(5) * (t192 * t247 + t236 * t273 - t237 * t275);
t128 = (t513 * t658 + t596) * t508 + (-t513 * t657 - t814) * t507;
t179 = (-t507 * t662 + t508 * t826) * t513;
t572 = t229 * t507 - t508 * t865;
t779 = m(6) * (t168 * t436 + t169 * t434 + (t179 * t516 + (t128 + t572) * t513) * t512);
t778 = m(6) * (t128 * t179 + t168 * t865 - t169 * t229);
t494 = pkin(7) * t513 + t744;
t654 = t855 * t494;
t175 = t507 * t826 + t508 * t662 + t654;
t318 = (t434 * t508 - t436 * t507) * t513;
t704 = t513 ^ 2 * t512;
t378 = t434 * t516 + t507 * t704;
t379 = -t436 * t516 - t508 * t704;
t776 = m(6) * (t175 * t318 + t179 * t347 - t314 * t379 - t316 * t378 + t572 * t703);
t198 = -t271 * t507 + t272 * t508;
t772 = m(6) * (-t314 * t435 - t316 * t437 - t372 * t434 - t373 * t436 + (t175 * t515 + t198 * t512) * t513);
t771 = m(6) * (-t436 * t229 + t241 * t379 + t378 * t866 - t434 * t865);
t770 = m(6) * (t175 * t198 + t314 * t372 + t316 * t373);
t767 = m(6) * (t241 * t435 + t271 * t436 + t272 * t434 + t437 * t866);
t763 = -t508 / 0.4e1;
t739 = rSges(4,1) * t516;
t598 = pkin(2) + t739;
t640 = rSges(4,2) * t710 + t508 * rSges(4,3);
t370 = -t507 * t598 + t597 + t640;
t490 = rSges(4,2) * t706;
t371 = t745 - t490 + t598 * t508 + (rSges(4,3) + pkin(6)) * t507;
t482 = rSges(4,1) * t513 + rSges(4,2) * t516;
t456 = t482 * t507;
t457 = t482 * t508;
t754 = m(4) * (t370 * t456 - t371 * t457);
t225 = t339 * t507 + t343 * t508 + t654;
t266 = t363 * t507 + t368 * t508;
t753 = m(5) * (t225 * t266 + (t374 * t507 + t376 * t508) * t473);
t751 = m(5) * (t269 * t345 + t270 * t346);
t750 = m(5) * (-t269 * t363 + t270 * t368);
t540 = (-t241 * t507 - t508 * t866) * t703;
t749 = m(6) * (t281 * t434 + t282 * t436 + t540);
t748 = m(6) * (-t436 * t314 + t434 * t316 + t540);
t747 = m(6) * (t241 * t272 + t271 * t866);
t746 = m(6) * (t241 * t281 + t282 * t866);
t742 = m(6) * qJD(3);
t741 = m(6) * qJD(4);
t55 = 0.2e1 * (t128 / 0.4e1 - t198 / 0.4e1) * m(6) + 0.2e1 * (t192 / 0.4e1 - t266 / 0.4e1) * m(5);
t691 = t55 * qJD(2);
t659 = t507 * t404 + t508 * t693;
t655 = t507 * (pkin(7) * t709 - t495) + t508 * (-pkin(3) * t706 + t496);
t643 = -t461 - t494;
t639 = qJD(1) * t513;
t638 = qJD(1) * t516;
t637 = qJD(3) * t507;
t636 = qJD(3) * t508;
t635 = qJD(4) * t513;
t634 = qJD(4) * t516;
t612 = t698 / 0.2e1;
t300 = (t612 - t318 / 0.2e1) * m(6);
t632 = t300 * qJD(2);
t526 = -t513 * t553 + t716;
t137 = t383 * t434 - t391 * t435 + t507 * t526;
t525 = -t513 * t552 + t715;
t138 = t384 * t434 - t392 * t435 + t507 * t525;
t528 = t513 * t551 + t718;
t139 = t389 * t434 - t393 * t435 + t507 * t528;
t527 = t513 * t550 + t717;
t140 = t390 * t434 - t394 * t435 + t507 * t527;
t627 = (-t834 + t887) * t822 + ((t138 + t140) * t508 + (t137 + t139) * t507 + t831) * t760;
t141 = t383 * t436 - t391 * t437 + t508 * t526;
t142 = t384 * t436 - t392 * t437 + t508 * t525;
t143 = t389 * t436 - t393 * t437 + t508 * t528;
t144 = t390 * t436 - t394 * t437 + t508 * t527;
t626 = (t576 + t575 - t833) * t822 + ((t142 + t144) * t508 + (t141 + t143) * t507 + t830) * t760;
t625 = (t838 * t507 + t837 * t508 + t828) * t823 + (t827 - t832) * t759;
t624 = t862 * t762 + t861 * t824;
t623 = t860 * t764 + t859 * t766;
t619 = -t494 - t644;
t616 = t710 / 0.4e1;
t610 = t828 * t822 + t827 * t823;
t606 = t170 / 0.2e1 + t172 / 0.2e1;
t605 = -t171 / 0.2e1 - t173 / 0.2e1;
t599 = t466 / 0.2e1 + t467 / 0.2e1;
t594 = -t403 + t701;
t586 = -t509 - t733;
t581 = -Icges(4,5) * t513 - Icges(4,6) * t516;
t571 = t314 * t507 + t316 * t508;
t560 = -t624 - t627;
t559 = t623 - t626;
t452 = -Icges(4,2) * t709 - t488;
t454 = t586 * t507;
t539 = (t405 - t454) * t516 + (t407 + t452) * t513;
t453 = t478 * t508;
t455 = t586 * t508;
t538 = (-t406 + t455) * t516 + (-t408 + t453) * t513;
t265 = -t508 * t701 + t659;
t537 = (t508 * t594 + t265 - t659) * t762 + (-t507 * (-t694 + t702) - t508 * t403) * t764 + (t507 * t594 + t595 + t829) * t766;
t536 = t265 * t766 + t659 * t824 + (-t380 + (t404 + t702) * t508 + t660 + t829) * t764;
t529 = -t806 + (t835 - t857) * t507 / 0.4e1 + (t836 - t858) * t763;
t519 = (-t710 / 0.4e1 + t616) * (t864 + t863) + (t763 + t508 / 0.4e1) * (t873 + t874);
t518 = -t805 + t828 * t760 + t832 * t759 + (t834 + t838) * t616 + (t831 + t868) * t709 / 0.4e1 + (t833 + t837) * t706 / 0.4e1 + (t830 + t867) * t705 / 0.4e1;
t484 = -rSges(4,2) * t513 + t739;
t451 = t581 * t508;
t450 = t581 * t507;
t377 = t643 * t508;
t375 = t643 * t507;
t348 = -t456 * t507 - t457 * t508;
t317 = t619 * t508;
t315 = t619 * t507;
t304 = m(6) * t611 + t347 * t786;
t301 = m(6) * t612 + t318 * t786;
t298 = -t368 * t516 - t473 * t706;
t297 = t363 * t516 + t473 * t710;
t267 = t434 * t435 + t436 * t437 + t515 * t704;
t260 = (t363 * t508 - t368 * t507) * t513;
t253 = -t516 * t466 + (t512 * t647 + t515 * t649) * t513;
t252 = -t516 * t467 + (t512 * t648 + t515 * t650) * t513;
t243 = -t396 * t507 + t398 * t508 + t655;
t239 = -t272 * t516 - t373 * t513;
t238 = -t271 * t516 + t642 * t710;
t203 = t507 * t658 + t508 * t657 + t655;
t191 = (-t271 * t508 - t272 * t507) * t513;
t102 = t748 / 0.2e1;
t101 = t175 * t347 + t571 * t703;
t99 = t749 / 0.2e1;
t94 = t767 / 0.2e1;
t78 = t179 * t318 - t229 * t379 + t378 * t865;
t71 = t771 / 0.2e1;
t70 = -t143 * t508 + t144 * t507;
t69 = -t141 * t508 + t142 * t507;
t68 = -t139 * t508 + t140 * t507;
t67 = -t137 * t508 + t138 * t507;
t63 = t772 / 0.2e1;
t56 = (t192 + t266) * t788 + (t128 + t198) * t786;
t52 = t513 * t801 - t599 * t516 + t747 + t750;
t47 = t776 / 0.2e1;
t46 = -t513 * t797 - t516 * t798 + t746 + t751 + t754;
t19 = t779 / 0.2e1;
t18 = t99 - t748 / 0.2e1;
t17 = t102 + t99;
t16 = t102 - t749 / 0.2e1;
t13 = t94 - t771 / 0.2e1;
t12 = t71 + t94;
t11 = t71 - t767 / 0.2e1;
t10 = t63 + t19 - t776 / 0.2e1;
t9 = t47 + t63 - t779 / 0.2e1;
t8 = t47 + t19 - t772 / 0.2e1;
t7 = t507 * t623 + t508 * t624 + t753 + t770;
t6 = t507 * t537 + t508 * t536;
t4 = t784 + t778 + (t507 * t622 + t508 * t621 + t625) * t516 + (t507 * t627 + t508 * t626 - t610) * t513;
t3 = t529 + t518;
t2 = t519 + (-t221 / 0.4e1 - t220 / 0.4e1 - t173 / 0.4e1 - t171 / 0.4e1) * t507 + (t219 / 0.4e1 + t218 / 0.4e1 + t172 / 0.4e1 + t170 / 0.4e1) * t508 + t518 + t806;
t1 = t529 + t519 + (t235 / 0.2e1 + t234 / 0.2e1 + (-t217 / 0.4e1 - t214 / 0.4e1 - t259 / 0.4e1 - t258 / 0.4e1) * t508 + (t216 / 0.4e1 + t213 / 0.4e1 - t255 / 0.4e1 - t254 / 0.4e1) * t507) * t516 + (-t284 / 0.2e1 - t283 / 0.2e1 + (-t161 / 0.4e1 - t159 / 0.4e1 - t202 / 0.4e1 - t201 / 0.4e1) * t508 + (-t160 / 0.4e1 - t158 / 0.4e1 - t200 / 0.4e1 - t199 / 0.4e1) * t507) * t513 + t805;
t5 = [t46 * qJD(3) + t52 * qJD(4) - t740 * t871, 0, t46 * qJD(1) + t3 * qJD(4) + t17 * qJD(5) + ((t269 * t377 + t270 * t375 - t345 * t376 - t346 * t374) * t788 + (t241 * t315 - t281 * t314 - t282 * t316 + t317 * t866) * t786) * t792 + (m(4) * (-t370 * t484 - t456 * t482) - t200 / 0.2e1 - t199 / 0.2e1 + t477 * t762 + (-t407 / 0.2e1 - t452 / 0.2e1) * t516 + (t405 / 0.2e1 - t454 / 0.2e1) * t513 - t536 - t808) * t636 + (m(4) * (-t371 * t484 + t457 * t482) + t201 / 0.2e1 + t202 / 0.2e1 + t477 * t766 + (t408 / 0.2e1 - t453 / 0.2e1) * t516 + (-t406 / 0.2e1 + t455 / 0.2e1) * t513 - t537 - t807) * t637, t52 * qJD(1) + t3 * qJD(3) + t12 * qJD(5) + ((-t229 * t272 + t238 * t866 + t239 * t241 + t271 * t865) * t786 + (t269 * t297 + t270 * t298 - t273 * t363 - t275 * t368) * t788) * t790 + (-t253 - t252) * t634 + ((t221 / 0.2e1 + t220 / 0.2e1 - t605) * t508 + (t218 / 0.2e1 + t219 / 0.2e1 + t606) * t507) * t635, t17 * qJD(3) + t12 * qJD(4) - t885; 0, 0, (m(4) * t348 / 0.2e1 + t243 * t788 + t203 * t786) * t792 + t56 * qJD(4) + t304 * qJD(5), t56 * qJD(3) + (t191 * t786 + t260 * t788) * t790 + t301 * qJD(5), qJD(3) * t304 + qJD(4) * t301; t6 * qJD(3) + t1 * qJD(4) + t16 * qJD(5) + (-t746 / 0.4e1 - t751 / 0.4e1 - t754 / 0.4e1) * t793 + t798 * t638 + t797 * t639, -qJD(4) * t55 - qJD(5) * t303, t6 * qJD(1) + t7 * qJD(4) + t101 * t740 + (m(6) * (t175 * t203 - t314 * t315 - t316 * t317) + m(5) * (t225 * t243 - t374 * t375 - t376 * t377) + m(4) * ((t507 * (rSges(4,1) * t709 - t640) + t508 * (rSges(4,1) * t705 + t507 * rSges(4,3) - t490)) * t348 + t855 * t484 * t482) + (t70 + t69 + t796 * t451 + (t539 * t508 + (-t450 + t538) * t507) * t508) * t766 + (t68 + t67 + t795 * t450 + (t538 * t507 + (-t451 + t539) * t508) * t507) * t764) * qJD(3), t1 * qJD(1) - t691 + t7 * qJD(3) + (t820 * t764 + t819 * t766) * qJD(4) + t9 * qJD(5) + (-t778 / 0.4e1 - t784 / 0.4e1) * t789 + ((t225 * t260 + t247 * t266 - t297 * t376 - t298 * t374 + (-t273 * t508 + t275 * t507) * t473) * t788 + (t175 * t191 + t179 * t198 + t229 * t372 - t238 * t316 - t239 * t314 - t373 * t865) * t786) * t790 + ((t606 - t621) * t508 + (t605 - t622) * t507 - t625) * t634 + (t507 * t560 + t508 * t559 + t610) * t635, t16 * qJD(1) + t9 * qJD(4) + t101 * t742 - t809; t2 * qJD(3) + t11 * qJD(5) + t599 * t638 + (-t747 / 0.4e1 - t750 / 0.4e1) * t793 - t801 * t639, qJD(3) * t55 - qJD(5) * t300, t2 * qJD(1) + t691 + t4 * qJD(4) + t8 * qJD(5) + 0.4e1 * (-t753 / 0.4e1 - t770 / 0.4e1) * qJD(3) + ((t128 * t175 - t168 * t316 - t169 * t314 + t179 * t203 - t229 * t315 + t317 * t865) * t786 + (t192 * t225 - t236 * t376 - t237 * t374 + t243 * t247 + t273 * t377 - t275 * t375) * t788) * t792 + ((t864 / 0.2e1 + t863 / 0.2e1 + t808) * t516 + (t70 / 0.2e1 + t69 / 0.2e1 + t213 / 0.2e1 + t216 / 0.2e1) * t513 + t560) * t636 + ((t889 * t764 + t888 * t766 + t807) * t516 + (t214 / 0.2e1 + t217 / 0.2e1 + t67 / 0.2e1 + t68 / 0.2e1) * t513 - t559) * t637, t4 * qJD(3) + (t253 / 0.2e1 + t252 / 0.2e1) * qJD(4) * t516 ^ 2 + (m(6) * (t179 * t191 - t229 * t239 + t238 * t865) / 0.4e1 + (t247 * t260 + t273 * t297 - t275 * t298) * m(5) / 0.4e1) * t789 + t78 * t740 + (t820 * t766 + t819 * t762 + (t507 * t836 + t508 * t835) * t759) * t635, t11 * qJD(1) - t632 + t8 * qJD(3) + t78 * t741 + (t318 * t703 + t378 * t436 + t379 * t434 - t267) * t740; t18 * qJD(3) + t13 * qJD(4) + t885, qJD(3) * t303 + qJD(4) * t300, t18 * qJD(1) + (t315 * t434 + t317 * t436 - t101 + (t175 * t516 + (t203 + t571) * t513) * t512) * t742 + t10 * qJD(4) + t809, t13 * qJD(1) + t632 + t10 * qJD(3) + (t437 * t865 - t435 * t229 + t436 * t238 + t434 * t239 + (t179 * t515 + t191 * t512) * t513 - t78) * t741 + t267 * t740, 0.4e1 * (t299 * qJD(3) / 0.4e1 + t267 * qJD(4) / 0.4e1) * m(6);];
Cq = t5;
