% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2020-01-03 12:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRP2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP2_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP2_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:11:06
% EndTime: 2020-01-03 12:11:26
% DurationCPUTime: 14.11s
% Computational Cost: add. (61343->500), mult. (45768->622), div. (0->0), fcn. (41936->8), ass. (0->340)
t844 = Icges(5,4) + Icges(6,4);
t843 = Icges(5,1) + Icges(6,1);
t837 = Icges(5,5) + Icges(6,5);
t841 = Icges(5,2) + Icges(6,2);
t836 = Icges(5,6) + Icges(6,6);
t513 = qJ(3) + qJ(4);
t505 = sin(t513);
t842 = t844 * t505;
t514 = qJ(1) + qJ(2);
t506 = sin(t514);
t507 = cos(t513);
t508 = cos(t514);
t624 = t507 * t508;
t638 = t505 * t508;
t829 = t836 * t506 + t844 * t624 - t841 * t638;
t840 = -t836 * t505 + t837 * t507;
t822 = t843 * t507 - t842;
t839 = t844 * t638;
t838 = t844 * t507;
t835 = Icges(5,3) + Icges(6,3);
t827 = t837 * t506 + t843 * t624 - t839;
t834 = -t841 * t505 + t838;
t828 = t822 * t506 - t837 * t508;
t823 = t841 * t507 + t842;
t833 = t840 * t506;
t832 = t829 * t505;
t831 = t843 * t505 + t838;
t808 = -t508 * t835 + t833;
t807 = t835 * t506 + t624 * t837 - t836 * t638;
t830 = t506 * t834 - t508 * t836;
t826 = -t507 * t827 + t832;
t631 = t506 * t507;
t825 = t828 * t631;
t824 = t827 * t631;
t821 = t505 * t837 + t836 * t507;
t820 = t841 * t624 - t827 + t839;
t819 = -t823 * t506 + t828;
t818 = t831 * t508 + t829;
t817 = -t831 * t506 - t830;
t816 = t831 + t834;
t688 = rSges(6,2) * t507;
t458 = rSges(6,1) * t505 + t688;
t698 = pkin(4) * t505;
t515 = sin(qJ(3));
t699 = pkin(3) * t515;
t465 = -t698 - t699;
t539 = t458 - t465;
t333 = t539 * t508;
t815 = t333 * t508;
t459 = rSges(5,1) * t505 + rSges(5,2) * t507;
t410 = t459 * t508;
t561 = t459 + t699;
t793 = t561 * t506;
t317 = t793 * t410;
t639 = t505 * t506;
t814 = t808 * t508 + t830 * t639 - t825;
t813 = -t807 * t508 - t829 * t639 + t824;
t517 = cos(qJ(3));
t510 = t517 * pkin(3);
t502 = t510 + pkin(2);
t697 = pkin(4) * t507;
t464 = t502 + t697;
t519 = -pkin(8) - pkin(7);
t512 = -qJ(5) + t519;
t576 = -rSges(6,2) * t639 - t508 * rSges(6,3);
t689 = rSges(6,1) * t507;
t307 = t508 * t512 + (t464 + t689) * t506 + t576;
t696 = sin(qJ(1)) * pkin(1);
t299 = t307 + t696;
t812 = -t299 + t307;
t789 = rSges(6,1) * t624 - rSges(6,2) * t638 + t508 * t464;
t308 = (rSges(6,3) - t512) * t506 + t789;
t511 = cos(qJ(1)) * pkin(1);
t300 = t511 + t308;
t811 = -t300 + t308;
t575 = -rSges(5,2) * t639 - t508 * rSges(5,3);
t616 = t508 * t519;
t690 = rSges(5,1) * t507;
t318 = t616 + (t502 + t690) * t506 + t575;
t309 = t318 + t696;
t810 = -t309 + t318;
t472 = t508 * t502;
t555 = rSges(5,1) * t624 - rSges(5,2) * t638;
t319 = t472 + (rSges(5,3) - t519) * t506 + t555;
t310 = t511 + t319;
t809 = -t310 + t319;
t509 = Icges(4,4) * t517;
t480 = -Icges(4,2) * t515 + t509;
t778 = Icges(4,1) * t515 + t509;
t804 = t480 + t778;
t803 = t828 * t507;
t802 = t830 * t505;
t801 = -t807 * t506 + t826 * t508 - t825;
t800 = t808 * t506 + t828 * t624;
t757 = m(4) / 0.2e1;
t756 = m(5) / 0.2e1;
t755 = m(6) / 0.2e1;
t734 = -t506 / 0.2e1;
t797 = t506 / 0.2e1;
t733 = -t508 / 0.2e1;
t730 = m(3) * (-t511 * (rSges(3,1) * t506 + rSges(3,2) * t508) + t696 * (t508 * rSges(3,1) - rSges(3,2) * t506));
t691 = rSges(4,1) * t517;
t485 = -rSges(4,2) * t515 + t691;
t794 = t485 * t757;
t792 = t561 * t508;
t791 = t821 * t506;
t790 = t821 * t508;
t503 = t506 ^ 2;
t504 = t508 ^ 2;
t570 = t503 + t504;
t462 = -rSges(5,2) * t505 + t690;
t645 = t462 * t508;
t283 = t318 * t645;
t439 = t506 * t465;
t335 = -t458 * t506 + t439;
t560 = t458 + t698;
t362 = t560 * t506;
t364 = t560 * t508;
t600 = t333 * t362 + t364 * t335;
t461 = -rSges(6,2) * t505 + t689;
t559 = -t461 - t697;
t363 = t559 * t506;
t365 = pkin(4) * t624 + t508 * t461;
t607 = t365 * t307 + t363 * t308;
t661 = t792 * t459;
t663 = t319 * t462;
t679 = (t283 - t317 + (t661 - t663) * t506) * t756 + (t600 + t607) * t755;
t408 = t459 * t506;
t596 = -t408 * t792 + t317;
t331 = t539 * t506;
t538 = -t688 + (-rSges(6,1) - pkin(4)) * t505;
t352 = t538 * t506;
t353 = t538 * t508;
t602 = -t331 * t353 + t333 * t352;
t646 = t462 * t506;
t680 = (-t319 * t646 + t283 + t596) * t756 + (t602 + t607) * t755;
t788 = t679 - t680;
t274 = t309 * t645;
t608 = t365 * t299 + t363 * t300;
t665 = t310 * t462;
t681 = (t274 - t317 + (t661 - t665) * t506) * t756 + (t600 + t608) * t755;
t682 = (-t310 * t646 + t274 + t596) * t756 + (t602 + t608) * t755;
t787 = t681 - t682;
t732 = t508 / 0.2e1;
t786 = t732 + t733;
t785 = (-t822 + t823) * t507 + t816 * t505;
t784 = (-t818 * t506 - t817 * t508) * t507 + (t820 * t506 + t819 * t508) * t505;
t483 = rSges(4,1) * t515 + rSges(4,2) * t517;
t433 = t483 * t506;
t434 = t483 * t508;
t138 = t299 * t308 - t307 * t300;
t146 = t309 * t319 - t318 * t310;
t500 = t508 * pkin(7);
t630 = t506 * t515;
t572 = -rSges(4,2) * t630 - t508 * rSges(4,3);
t337 = -t500 + (pkin(2) + t691) * t506 + t572;
t322 = t337 + t696;
t617 = t508 * t517;
t618 = t508 * t515;
t537 = rSges(4,1) * t617 - rSges(4,2) * t618 + rSges(4,3) * t506;
t571 = t508 * pkin(2) + t506 * pkin(7);
t338 = t537 + t571;
t323 = t338 + t511;
t159 = t322 * t338 - t337 * t323;
t339 = t506 * (t616 + t500 + (-pkin(2) + t502) * t506);
t356 = t506 * t519 - t472 + t571;
t569 = t512 - t519;
t594 = -t472 + t789 + (rSges(6,3) - t569) * t506;
t595 = (rSges(6,1) * t631 + t569 * t508 + t576 + (t464 - t502) * t506) * t506;
t118 = t339 + (-t356 + t594) * t508 + t595;
t384 = t503 * t458;
t409 = t458 * t508;
t247 = -t409 * t508 - t570 * t698 - t384;
t601 = -t331 * t363 + t333 * t365;
t55 = t118 * t247 + t601;
t343 = t506 * (rSges(5,1) * t631 + t575);
t379 = t506 * rSges(5,3) + t555;
t163 = t339 + t343 + (-t356 + t379) * t508;
t305 = -t506 * t408 - t410 * t508;
t91 = t163 * t305 + t645 * t792 + t646 * t793;
t777 = -m(5) * t91 - m(6) * t55;
t568 = qJD(1) + qJD(2);
t683 = (t810 * t506 + t809 * t508) * t459 * t756 + (t812 * t362 + t811 * t364) * t755;
t153 = t299 * t352 + t300 * t353;
t154 = t307 * t352 + t308 * t353;
t161 = -t309 * t408 - t310 * t410;
t164 = -t318 * t408 - t319 * t410;
t775 = (t154 + t153) * t755 + (t164 + t161) * t756;
t566 = ((-t323 + t338) * t508 + (-t322 + t337) * t506) * t483 * t757 + (t812 * t331 + t811 * t333) * t755 + (t809 * t792 + t810 * t793) * t756;
t148 = t299 * t335 - t333 * t300;
t152 = t307 * t335 - t333 * t308;
t157 = -t309 * t793 - t310 * t792;
t158 = -t318 * t793 - t319 * t792;
t175 = -t322 * t433 - t323 * t434;
t185 = -t337 * t433 - t338 * t434;
t773 = (t185 + t175) * t757 + (t152 + t148) * t755 + (t158 + t157) * t756;
t672 = Icges(4,4) * t515;
t479 = Icges(4,2) * t517 + t672;
t482 = Icges(4,1) * t517 - t672;
t768 = t804 * t515 + (t479 - t482) * t517;
t489 = Icges(4,4) * t618;
t394 = Icges(4,1) * t617 + Icges(4,5) * t506 - t489;
t581 = -Icges(4,2) * t617 + t394 - t489;
t392 = Icges(4,4) * t617 - Icges(4,2) * t618 + Icges(4,6) * t506;
t583 = t508 * t778 + t392;
t767 = -t515 * t581 - t517 * t583;
t393 = -Icges(4,5) * t508 + t482 * t506;
t582 = -t479 * t506 + t393;
t391 = -Icges(4,6) * t508 + t480 * t506;
t584 = t506 * t778 + t391;
t766 = -t515 * t582 - t517 * t584;
t765 = t804 * t517 / 0.2e1 + (t482 / 0.2e1 - t479 / 0.2e1) * t515;
t535 = (t813 * t506 + t814 * t508) * t797 + (((t808 - t826) * t508 + t801) * t508 + ((t808 - t832) * t506 + (t802 + t803) * t508 - t800 + t824) * t506) * t734 + (((-t803 + t807) * t508 + t800 + t813) * t508 + ((t802 + t807) * t506 + t801 - t814) * t506) * t733;
t534 = t816 * t507 / 0.2e1 + (-t823 / 0.2e1 + t822 / 0.2e1) * t505;
t763 = 4 * qJD(1);
t761 = 4 * qJD(2);
t760 = 2 * qJD(3);
t759 = 2 * qJD(4);
t156 = t594 * t508 + t595;
t182 = t506 * (pkin(3) * t630 + t439) - t384 + (-(-t465 - t699) * t508 - t409) * t508;
t598 = -t362 * t363 + t364 * t365;
t744 = m(6) * (t156 * t182 + t598);
t629 = t506 * t517;
t348 = t393 * t629;
t478 = Icges(4,5) * t517 - Icges(4,6) * t515;
t644 = t478 * t506;
t389 = -Icges(4,3) * t508 + t644;
t206 = -t389 * t508 - t391 * t630 + t348;
t349 = t394 * t629;
t390 = Icges(4,5) * t617 - Icges(4,6) * t618 + Icges(4,3) * t506;
t207 = t390 * t508 + t392 * t630 - t349;
t144 = -t206 * t508 - t207 * t506;
t350 = t391 * t618;
t208 = -t389 * t506 - t393 * t617 + t350;
t653 = t392 * t515;
t541 = t394 * t517 - t653;
t209 = t390 * t506 + t541 * t508;
t145 = -t208 * t508 - t209 * t506;
t652 = t393 * t517;
t654 = t391 * t515;
t46 = (t208 + t349 - t350 + (t389 - t653) * t506) * t506 + (-t348 - t209 + (t389 + t541) * t508 + (t652 + t654) * t506) * t508;
t47 = (-t207 + t350 + (t390 - t652) * t508) * t508 + (t206 - t348 + (t390 + t654) * t506) * t506;
t4 = (-t145 / 0.2e1 - t47 / 0.2e1) * t508 + (-t46 / 0.2e1 + t144 / 0.2e1) * t506 + t535;
t702 = m(6) * (-t331 * t506 - t815);
t218 = t702 / 0.2e1;
t703 = m(6) * (-t335 * t506 + t815);
t98 = t218 - t703 / 0.2e1;
t731 = t4 * qJD(3) + t98 * qJD(5);
t726 = m(4) * t159;
t724 = m(4) * t175;
t723 = m(4) * t185;
t717 = m(5) * t146;
t280 = t379 * t508 + t343;
t151 = t570 * t459 * t462 + t280 * t305;
t149 = m(5) * t151;
t715 = m(5) * t157;
t714 = m(5) * t158;
t713 = m(5) * t161;
t712 = m(5) * t164;
t286 = t300 * t506;
t296 = t308 * t506;
t711 = m(6) * (t812 * t508 + t286 - t296);
t710 = m(6) * (t286 + t296 + (-t299 - t307) * t508);
t709 = m(6) * t138;
t707 = m(6) * t148;
t706 = m(6) * t152;
t705 = m(6) * t153;
t704 = m(6) * t154;
t701 = m(6) * (-t352 * t506 - t353 * t508);
t700 = m(6) * (-t362 * t506 - t364 * t508);
t257 = t700 / 0.2e1;
t136 = t257 - t701 / 0.2e1;
t695 = qJD(4) * t535 + t136 * qJD(5);
t693 = m(6) * qJD(3);
t692 = m(6) * qJD(5);
t238 = t701 / 0.2e1;
t134 = t238 - t700 / 0.2e1;
t217 = t703 / 0.2e1;
t97 = t217 - t702 / 0.2e1;
t678 = t97 * qJD(3) + t134 * qJD(4);
t135 = t257 + t238;
t99 = t218 + t217;
t677 = t99 * qJD(3) + t135 * qJD(4);
t651 = t410 * t459;
t604 = (t331 + t335) * t333;
t599 = t364 * t352 - t362 * t353;
t172 = -t299 * t508 + t286;
t564 = m(6) * t172 * qJD(1);
t177 = -t307 * t508 + t296;
t563 = t177 * m(6) * qJD(2);
t562 = t156 * t247 + t598;
t558 = ((-t791 * t506 - t784) * t508 + t790 * t503) * t734 + ((t790 * t508 + t784) * t506 - t791 * t504) * t733;
t557 = t570 * t699;
t556 = t149 + t558;
t549 = Icges(4,5) * t515 + Icges(4,6) * t517;
t528 = t534 + t775;
t527 = t534 + t765;
t526 = t527 + t773;
t525 = -t535 + (t818 * t505 + t820 * t507 + t785 * t508 - t833) * t734 + (t817 * t505 - t785 * t506 + t819 * t507 - t840 * t508) * t733;
t524 = -t534 + (t827 * t505 + t829 * t507) * t786;
t522 = t525 * qJD(4) + t135 * qJD(5);
t521 = t524 - t765 + t786 * (t517 * t392 + t515 * t394);
t520 = t99 * qJD(5) + (t46 * t797 + t525 + (t768 * t508 + t583 * t515 - t581 * t517 + t144 - t644) * t734 + (-t478 * t508 - t768 * t506 - t584 * t515 + t582 * t517) * t733 + (t47 + t145) * t732) * qJD(3);
t492 = pkin(3) * t617;
t428 = t549 * t508;
t427 = t506 * t549;
t383 = t492 + t645;
t381 = (-t462 - t510) * t506;
t334 = t492 + t365;
t332 = (t559 - t510) * t506;
t324 = t408 * t410;
t269 = -t363 * t508 - t365 * t506;
t262 = t305 - t557;
t249 = m(6) * t269 * qJD(4);
t166 = -t557 + t182;
t103 = t710 / 0.2e1;
t102 = t711 / 0.2e1;
t52 = t534 + t704 + t712;
t51 = t534 + t705 + t713;
t43 = t103 - t711 / 0.2e1;
t42 = t103 + t102;
t41 = t102 - t710 / 0.2e1;
t36 = t709 + t717 + t726 + t730;
t35 = t527 + t706 + t714 + t723;
t22 = t527 + t707 + t715 + t724;
t18 = t528 + t683;
t17 = t528 - t683;
t16 = t556 + t744;
t15 = t558 - t777;
t14 = t524 + t683 - t775;
t13 = t526 + t566;
t12 = t526 - t566;
t11 = t521 + t566 - t773;
t8 = t535 + t788;
t7 = t535 - t788;
t6 = t535 + t787;
t5 = t535 - t787;
t2 = t525 + t679 + t680;
t1 = t525 + t681 + t682;
t3 = [t36 * qJD(2) + t22 * qJD(3) + t51 * qJD(4) + t172 * t692, t36 * qJD(1) + t13 * qJD(3) + t18 * qJD(4) + t42 * qJD(5) + 0.2e1 * (t730 / 0.2e1 + t138 * t755 + t146 * t756 + t159 * t757) * qJD(2), t22 * qJD(1) + t13 * qJD(2) + t1 * qJD(4) + ((t322 * t508 - t323 * t506) * t794 + (t309 * t383 + t310 * t381) * t756 + (t299 * t334 + t300 * t332 + t604) * t755) * t760 + t520, t51 * qJD(1) + t18 * qJD(2) + t1 * qJD(3) + ((t274 - t324 + (t651 - t665) * t506) * t756 + (t599 + t608) * t755) * t759 + t522, t42 * qJD(2) + t564 + t677; t12 * qJD(3) + t17 * qJD(4) + t43 * qJD(5) + (-t730 / 0.4e1 - t726 / 0.4e1 - t717 / 0.4e1 - t709 / 0.4e1) * t763, t35 * qJD(3) + t52 * qJD(4) + t177 * t692, t12 * qJD(1) + t35 * qJD(2) + t2 * qJD(4) + ((t337 * t508 - t338 * t506) * t794 + (t318 * t383 + t319 * t381) * t756 + (t307 * t334 + t308 * t332 + t604) * t755) * t760 + t520, t17 * qJD(1) + t52 * qJD(2) + t2 * qJD(3) + ((t599 + t607) * t755 + (t283 - t324 + (t651 - t663) * t506) * t756) * t759 + t522, t43 * qJD(1) + t563 + t677; t521 * qJD(1) + t11 * qJD(2) + t5 * qJD(4) + (-t724 / 0.4e1 - t715 / 0.4e1 - t707 / 0.4e1) * t763 + t731, t11 * qJD(1) + t521 * qJD(2) + t7 * qJD(4) + (-t723 / 0.4e1 - t706 / 0.4e1 - t714 / 0.4e1) * t761 + t731, (m(6) * (t118 * t166 - t331 * t332 + t333 * t334) + m(5) * (t163 * t262 - t381 * t793 + t383 * t792) + (t503 * t428 + (t766 * t508 + (-t427 - t767) * t506) * t508) * t734 + (-t504 * t427 + (t767 * t506 + (t428 - t766) * t508) * t506) * t733 + m(4) * ((t508 * t537 + t506 * (rSges(4,1) * t629 + t572)) * (-t433 * t506 - t434 * t508) + t570 * t485 * t483) + t558) * qJD(3) + t15 * qJD(4) + t568 * t4, t5 * qJD(1) + t7 * qJD(2) + t15 * qJD(3) + ((t562 + t55) * t755 + (t91 + t151) * t756) * t759 + (t558 - t149 - t744) * qJD(4), t568 * t98; t524 * qJD(1) + t14 * qJD(2) + t6 * qJD(3) + (-t713 / 0.4e1 - t705 / 0.4e1) * t763 + t695, t14 * qJD(1) + t524 * qJD(2) + t8 * qJD(3) + (-t704 / 0.4e1 - t712 / 0.4e1) * t761 + t695, t6 * qJD(1) + t8 * qJD(2) + t16 * qJD(4) + ((t118 * t182 + t156 * t166 - t332 * t362 + t334 * t364 + t601) * t755 + (t262 * t280 + (-t381 * t506 + t383 * t508) * t459 + t91) * t756) * t760 + (t558 + t777) * qJD(3), t16 * qJD(3) + (m(6) * t562 + t556) * qJD(4) + t568 * t535, t568 * t136; t41 * qJD(2) - t564 + t678, t41 * qJD(1) - t563 + t678, (-t332 * t508 - t334 * t506) * t693 + t249 + t568 * t97, t134 * t568 + t269 * t693 + t249, 0;];
Cq = t3;
