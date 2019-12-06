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
% Datum: 2019-12-05 18:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:47:34
% EndTime: 2019-12-05 18:47:55
% DurationCPUTime: 13.97s
% Computational Cost: add. (61343->509), mult. (45768->623), div. (0->0), fcn. (41936->8), ass. (0->338)
t859 = Icges(5,4) + Icges(6,4);
t856 = Icges(5,5) + Icges(6,5);
t852 = Icges(5,2) + Icges(6,2);
t855 = Icges(5,6) + Icges(6,6);
t533 = qJ(3) + qJ(4);
t528 = cos(t533);
t860 = t859 * t528;
t857 = Icges(5,1) + Icges(6,1);
t526 = sin(t533);
t858 = -t855 * t526 + t856 * t528;
t834 = t852 * t526 - t860;
t854 = Icges(5,3) + Icges(6,3);
t534 = qJ(1) + qJ(2);
t529 = cos(t534);
t653 = t526 * t529;
t853 = t859 * t653;
t527 = sin(t534);
t848 = t527 * t834 + t855 * t529;
t643 = t528 * t529;
t847 = t527 * t855 + t643 * t859 - t653 * t852;
t845 = t856 * t527 + t643 * t857 - t853;
t654 = t526 * t527;
t851 = t859 * t654;
t850 = t858 * t527;
t833 = t526 * t857 + t860;
t849 = t854 * t527 + t856 * t643 - t855 * t653;
t650 = t527 * t528;
t846 = -t856 * t529 + t857 * t650 - t851;
t844 = t859 * t526;
t842 = t528 * t852 + t844;
t841 = t857 * t528 - t844;
t840 = t856 * t526 + t855 * t528;
t839 = -t643 * t852 + t845 - t853;
t838 = t650 * t852 - t846 + t851;
t837 = t529 * t833 + t847;
t836 = -t527 * t833 + t848;
t835 = t848 * t654 + (t854 * t529 - t850) * t529;
t832 = (t526 * t847 - t528 * t845) * t529;
t831 = t846 * t650 + t835;
t830 = t849 * t529 - t845 * t650 + t847 * t654;
t537 = cos(qJ(3));
t531 = t537 * pkin(3);
t523 = t531 + pkin(2);
t710 = pkin(4) * t528;
t476 = t523 + t710;
t539 = -pkin(8) - pkin(7);
t532 = -qJ(5) + t539;
t589 = -rSges(6,2) * t654 - t529 * rSges(6,3);
t700 = rSges(6,1) * t528;
t311 = t529 * t532 + (t476 + t700) * t527 + t589;
t709 = sin(qJ(1)) * pkin(1);
t302 = t311 + t709;
t829 = -t302 + t311;
t312 = -rSges(6,1) * t643 + rSges(6,2) * t653 - t529 * t476 - (rSges(6,3) - t532) * t527;
t713 = cos(qJ(1)) * pkin(1);
t303 = t312 - t713;
t828 = t303 - t312;
t701 = rSges(5,1) * t528;
t575 = t523 + t701;
t588 = -rSges(5,2) * t654 - t529 * rSges(5,3);
t634 = t529 * t539;
t321 = t527 * t575 + t588 + t634;
t314 = t321 + t709;
t827 = -t314 + t321;
t826 = t846 * t528 - t849;
t825 = t848 * t526;
t824 = t849 * t527 - t832 + t835;
t823 = t848 * t653;
t822 = -(t527 * t839 + t529 * t838) * t526 - (t527 * t837 + t529 * t836) * t528;
t773 = m(4) / 0.2e1;
t772 = m(5) / 0.2e1;
t771 = m(6) / 0.2e1;
t819 = -t527 / 0.2e1;
t751 = t527 / 0.2e1;
t749 = t529 / 0.2e1;
t747 = m(3) * (-t713 * (rSges(3,1) * t527 + rSges(3,2) * t529) - (-rSges(3,1) * t529 + t527 * rSges(3,2)) * t709);
t638 = t529 * t311;
t177 = -t312 * t527 - t638;
t815 = t177 * m(6) * qJD(2);
t814 = t840 * t527;
t813 = t840 * t529;
t524 = t527 ^ 2;
t525 = t529 ^ 2;
t583 = t524 + t525;
t515 = t527 * t539;
t574 = rSges(5,2) * t653 - rSges(5,3) * t527;
t322 = -t529 * t575 + t515 + t574;
t474 = -rSges(5,2) * t526 + t701;
t663 = t474 * t527;
t285 = t322 * t663;
t471 = rSges(5,1) * t526 + rSges(5,2) * t528;
t535 = sin(qJ(3));
t712 = pkin(3) * t535;
t566 = (t471 + t712) * t529;
t664 = t471 * t527;
t320 = t566 * t664;
t420 = rSges(6,1) * t654 + rSges(6,2) * t650;
t711 = pkin(4) * t526;
t477 = -t711 - t712;
t340 = -t477 * t527 + t420;
t587 = rSges(6,1) * t653 + rSges(6,2) * t643;
t341 = -t477 * t529 + t587;
t470 = rSges(6,1) * t526 + rSges(6,2) * t528;
t496 = pkin(4) * t654;
t372 = t527 * t470 + t496;
t374 = (t470 + t711) * t529;
t618 = -t374 * t340 + t372 * t341;
t473 = -rSges(6,2) * t526 + t700;
t373 = pkin(4) * t650 + t527 * t473;
t577 = -t473 - t710;
t375 = t577 * t529;
t624 = -t375 * t311 + t373 * t312;
t421 = rSges(5,1) * t654 + rSges(5,2) * t650;
t649 = t527 * t535;
t513 = pkin(3) * t649;
t366 = t513 + t421;
t672 = t366 * t471;
t677 = t321 * t474;
t691 = (t285 + t320 + (-t672 + t677) * t529) * t772 + (t618 + t624) * t771;
t392 = t513 + t664;
t422 = t471 * t529;
t614 = t392 * t422 - t421 * t566;
t336 = t513 + t372;
t338 = (t470 - t477) * t529;
t358 = t496 + t420;
t359 = pkin(4) * t653 + t587;
t620 = t336 * t359 - t338 * t358;
t662 = t474 * t529;
t692 = (t321 * t662 + t285 + t614) * t772 + (t620 + t624) * t771;
t811 = t691 - t692;
t315 = t322 - t713;
t275 = t315 * t663;
t627 = -t375 * t302 + t373 * t303;
t678 = t314 * t474;
t693 = (t275 + t320 + (-t672 + t678) * t529) * t772 + (t618 + t627) * t771;
t694 = (t314 * t662 + t275 + t614) * t772 + (t620 + t627) * t771;
t810 = t693 - t694;
t750 = -t529 / 0.2e1;
t809 = t749 + t750;
t806 = (-t841 + t842) * t528 + (t833 - t834) * t526;
t288 = t529 * t302;
t172 = -t303 * t527 - t288;
t804 = t172 * m(6) * qJD(1);
t702 = rSges(4,1) * t537;
t505 = -rSges(4,2) * t535 + t702;
t803 = t505 * t773;
t635 = t529 * t537;
t636 = t529 * t535;
t404 = Icges(4,4) * t635 - Icges(4,2) * t636 + Icges(4,6) * t527;
t510 = Icges(4,4) * t636;
t406 = Icges(4,1) * t635 + Icges(4,5) * t527 - t510;
t802 = (t404 * t535 - t406 * t537) * t529;
t240 = t566 * t315;
t255 = t566 * t322;
t138 = -t312 * t302 + t303 * t311;
t146 = -t322 * t314 + t315 * t321;
t522 = t529 * pkin(7);
t578 = pkin(2) + t702;
t584 = -rSges(4,2) * t649 - t529 * rSges(4,3);
t342 = t527 * t578 - t522 + t584;
t327 = t342 + t709;
t512 = rSges(4,2) * t636;
t343 = t512 - t578 * t529 + (-rSges(4,3) - pkin(7)) * t527;
t328 = t343 - t713;
t159 = -t343 * t327 + t328 * t342;
t530 = Icges(4,4) * t537;
t797 = Icges(4,1) * t535 + t530;
t794 = Icges(4,2) * t535 - t530;
t708 = pkin(2) - t523;
t344 = t529 * (pkin(7) * t527 + t708 * t529 + t515);
t362 = t708 * t527 - t522 - t634;
t608 = (-t523 * t529 - t312 + t515) * t529;
t613 = -(-t532 + t539) * t529 - (-t476 + t523) * t527 + rSges(6,1) * t650 + t589;
t118 = -t344 + (-t362 + t613) * t527 + t608;
t396 = t529 * t587;
t247 = -t420 * t527 - t583 * t711 - t396;
t619 = t336 * t373 - t338 * t375;
t55 = t118 * t247 + t619;
t353 = t529 * (rSges(5,1) * t643 - t574);
t389 = -rSges(5,1) * t650 - t588;
t163 = -t344 + t353 + (-t362 - t389) * t527;
t308 = -t421 * t527 - t529 * t422;
t91 = t163 * t308 + t392 * t663 + t566 * t662;
t793 = -m(5) * t91 - m(6) * t55;
t582 = qJD(1) + qJD(2);
t695 = ((t315 - t322) * t422 + t827 * t664) * t772 + (t829 * t372 + t828 * t374) * t771;
t153 = -t302 * t358 + t303 * t359;
t154 = -t311 * t358 + t312 * t359;
t161 = -t314 * t421 + t315 * t422;
t164 = -t321 * t421 + t322 * t422;
t792 = (t154 + t153) * t771 + (t164 + t161) * t772;
t503 = rSges(4,1) * t535 + rSges(4,2) * t537;
t446 = t503 * t529;
t445 = t503 * t527;
t674 = t342 * t445;
t675 = t327 * t445;
t581 = (-t675 + t674 + (t328 - t343) * t446) * t773 + (t829 * t336 + t828 * t338) * t771 + (t827 * t392 + t240 - t255) * t772;
t148 = -t302 * t340 + t303 * t341;
t152 = -t311 * t340 + t312 * t341;
t157 = -t314 * t366 + t240;
t158 = -t321 * t366 + t255;
t175 = t328 * t446 - t675;
t185 = t343 * t446 - t674;
t791 = (t185 + t175) * t773 + (t152 + t148) * t771 + (t158 + t157) * t772;
t685 = Icges(4,4) * t535;
t499 = Icges(4,2) * t537 + t685;
t502 = Icges(4,1) * t537 - t685;
t784 = (-t794 + t797) * t535 + (t499 - t502) * t537;
t594 = -Icges(4,2) * t635 + t406 - t510;
t596 = t529 * t797 + t404;
t783 = t535 * t594 + t537 * t596;
t509 = Icges(4,4) * t649;
t648 = t527 * t537;
t405 = -Icges(4,1) * t648 + Icges(4,5) * t529 + t509;
t595 = Icges(4,2) * t648 + t405 + t509;
t403 = Icges(4,6) * t529 + t794 * t527;
t597 = -t527 * t797 + t403;
t782 = -t535 * t595 - t537 * t597;
t781 = (t797 / 0.2e1 - t794 / 0.2e1) * t537 + (t502 / 0.2e1 - t499 / 0.2e1) * t535;
t552 = (t830 * t527 + t831 * t529) * t819 + ((t824 + t832) * t529 + ((-t825 + t826) * t529 + t823 + t830) * t527) * t751 + (((t825 - t849) * t529 - t823 + t830) * t529 + (t826 * t527 + t824 - t831) * t527) * t749;
t551 = (-t834 / 0.2e1 + t833 / 0.2e1) * t528 + (-t842 / 0.2e1 + t841 / 0.2e1) * t526;
t779 = 0.4e1 * qJD(1);
t777 = 0.4e1 * qJD(2);
t776 = 2 * qJD(3);
t775 = 2 * qJD(4);
t156 = t527 * t613 + t608;
t182 = t525 * (t477 + t712) - t396 + (t513 - t340) * t527;
t616 = t372 * t373 - t374 * t375;
t760 = m(6) * (t156 * t182 + t616);
t498 = Icges(4,5) * t537 - Icges(4,6) * t535;
t659 = t498 * t527;
t401 = Icges(4,3) * t529 - t659;
t607 = t529 * t401 + t403 * t649;
t206 = -t405 * t648 + t607;
t402 = Icges(4,5) * t635 - Icges(4,6) * t636 + Icges(4,3) * t527;
t207 = t529 * t402 + t404 * t649 - t406 * t648;
t144 = t206 * t529 + t207 * t527;
t606 = t527 * t401 + t405 * t635;
t208 = -t403 * t636 + t606;
t209 = t527 * t402 - t802;
t145 = t208 * t529 + t209 * t527;
t571 = -t405 * t537 - t402;
t669 = t403 * t535;
t46 = (t209 + t607 + t802) * t529 + (-t208 + (t571 - t669) * t529 + t207 + t606) * t527;
t47 = (t207 + (-t402 + t669) * t529 - t606) * t529 + (t527 * t571 - t206 + t607) * t527;
t4 = (t145 / 0.2e1 + t47 / 0.2e1) * t529 + (t46 / 0.2e1 - t144 / 0.2e1) * t527 + t552;
t717 = m(6) * (-t336 * t527 - t529 * t338);
t218 = t717 / 0.2e1;
t718 = m(6) * (t340 * t527 + t341 * t529);
t98 = t218 - t718 / 0.2e1;
t748 = t4 * qJD(3) + t98 * qJD(5);
t743 = m(4) * t159;
t741 = m(4) * t175;
t740 = m(4) * t185;
t734 = m(5) * t146;
t281 = -t389 * t527 + t353;
t151 = t583 * t471 * t474 + t281 * t308;
t149 = m(5) * t151;
t732 = m(5) * t157;
t731 = m(5) * t158;
t730 = m(5) * t161;
t729 = m(5) * t164;
t728 = m(6) * (-t828 * t527 - t288 + t638);
t727 = m(6) * (-t288 - t638 + (-t303 - t312) * t527);
t726 = m(6) * t138;
t724 = m(6) * t148;
t723 = m(6) * t152;
t722 = m(6) * t153;
t721 = m(6) * t154;
t716 = m(6) * (t358 * t527 + t359 * t529);
t714 = m(6) * (-t372 * t527 - t529 * t374);
t257 = t714 / 0.2e1;
t136 = t257 - t716 / 0.2e1;
t707 = qJD(4) * t552 + t136 * qJD(5);
t704 = m(6) * qJD(3);
t703 = m(6) * qJD(5);
t238 = t716 / 0.2e1;
t134 = t238 - t714 / 0.2e1;
t217 = t718 / 0.2e1;
t97 = t217 - t717 / 0.2e1;
t690 = t97 * qJD(3) + t134 * qJD(4);
t135 = t257 + t238;
t99 = t218 + t217;
t689 = t99 * qJD(3) + t135 * qJD(4);
t668 = t421 * t471;
t622 = t336 * t341 - t338 * t340;
t617 = -t374 * t358 + t372 * t359;
t615 = (-t366 + t392) * t566;
t579 = t156 * t247 + t616;
t576 = ((t814 * t527 + t822) * t529 - t813 * t524) * t751 + ((-t813 * t529 - t822) * t527 + t814 * t525) * t749;
t570 = t583 * t712;
t569 = t149 + t576;
t562 = Icges(4,5) * t535 + Icges(4,6) * t537;
t439 = t527 * t562;
t548 = t551 + t792;
t547 = t551 + t781;
t546 = t547 + t791;
t545 = -t552 + (-t526 * t837 + t528 * t839 - t806 * t529 + t850) * t751 + (-t836 * t526 + t806 * t527 + t838 * t528 + t858 * t529) * t749;
t544 = -t551 + (t845 * t526 + t847 * t528) * t809;
t542 = t545 * qJD(4) + t135 * qJD(5);
t541 = t544 - t781 + t809 * (t537 * t404 + t535 * t406);
t540 = t99 * qJD(5) + (t46 * t819 + t545 + (-t784 * t529 - t535 * t596 + t537 * t594 + t144 + t659) * t751 + (t145 + t47) * t750 + (t498 * t529 + t784 * t527 - t535 * t597 + t537 * t595) * t749) * qJD(3);
t514 = pkin(3) * t648;
t440 = t562 * t529;
t395 = (-t474 - t531) * t529;
t393 = t514 + t663;
t339 = (t577 - t531) * t529;
t337 = t514 + t373;
t329 = t422 * t664;
t269 = t373 * t529 + t375 * t527;
t262 = t308 - t570;
t249 = m(6) * t269 * qJD(4);
t166 = -t570 + t182;
t103 = t727 / 0.2e1;
t102 = t728 / 0.2e1;
t52 = t551 + t721 + t729;
t51 = t551 + t722 + t730;
t43 = t103 - t728 / 0.2e1;
t42 = t103 + t102;
t41 = t102 - t727 / 0.2e1;
t36 = t726 + t734 + t743 + t747;
t35 = t547 + t723 + t731 + t740;
t22 = t547 + t724 + t732 + t741;
t18 = t548 + t695;
t17 = t548 - t695;
t16 = t569 + t760;
t15 = t576 - t793;
t14 = t544 + t695 - t792;
t13 = t546 + t581;
t12 = t546 - t581;
t11 = t541 + t581 - t791;
t8 = t552 + t811;
t7 = t552 - t811;
t6 = t552 + t810;
t5 = t552 - t810;
t2 = t545 + t691 + t692;
t1 = t545 + t693 + t694;
t3 = [t36 * qJD(2) + t22 * qJD(3) + t51 * qJD(4) + t172 * t703, t36 * qJD(1) + t13 * qJD(3) + t18 * qJD(4) + t42 * qJD(5) + 0.2e1 * (t747 / 0.2e1 + t138 * t771 + t146 * t772 + t159 * t773) * qJD(2), t22 * qJD(1) + t13 * qJD(2) + t1 * qJD(4) + ((t327 * t529 + t328 * t527) * t803 + (-t314 * t395 + t315 * t393 + t615) * t772 + (-t302 * t339 + t303 * t337 + t622) * t771) * t776 + t540, t51 * qJD(1) + t18 * qJD(2) + t1 * qJD(3) + ((t275 + t329 + (-t668 + t678) * t529) * t772 + (t617 + t627) * t771) * t775 + t542, t42 * qJD(2) + t689 + t804; t12 * qJD(3) + t17 * qJD(4) + t43 * qJD(5) + (-t747 / 0.4e1 - t743 / 0.4e1 - t734 / 0.4e1 - t726 / 0.4e1) * t779, t35 * qJD(3) + t52 * qJD(4) + t177 * t703, t12 * qJD(1) + t35 * qJD(2) + t2 * qJD(4) + ((t342 * t529 + t343 * t527) * t803 + (-t321 * t395 + t322 * t393 + t615) * t772 + (-t311 * t339 + t312 * t337 + t622) * t771) * t776 + t540, t17 * qJD(1) + t52 * qJD(2) + t2 * qJD(3) + ((t617 + t624) * t771 + (t285 + t329 + (-t668 + t677) * t529) * t772) * t775 + t542, t43 * qJD(1) + t689 + t815; t541 * qJD(1) + t11 * qJD(2) + t5 * qJD(4) + (-t741 / 0.4e1 - t724 / 0.4e1 - t732 / 0.4e1) * t779 + t748, t11 * qJD(1) + t541 * qJD(2) + t7 * qJD(4) + (-t740 / 0.4e1 - t723 / 0.4e1 - t731 / 0.4e1) * t777 + t748, (m(6) * (t118 * t166 + t336 * t337 - t338 * t339) + m(5) * (t163 * t262 + t392 * t393 - t395 * t566) + m(4) * ((t529 * (rSges(4,1) * t635 + rSges(4,3) * t527 - t512) - t527 * (-rSges(4,1) * t648 - t584)) * (-t445 * t527 - t446 * t529) + t583 * t505 * t503) + (-t524 * t440 + (t782 * t529 + (t439 - t783) * t527) * t529) * t751 + (t525 * t439 + (t783 * t527 + (-t440 - t782) * t529) * t527) * t749 + t576) * qJD(3) + t15 * qJD(4) + t582 * t4, t5 * qJD(1) + t7 * qJD(2) + t15 * qJD(3) + ((t579 + t55) * t771 + (t91 + t151) * t772) * t775 + (t576 - t149 - t760) * qJD(4), t582 * t98; t544 * qJD(1) + t14 * qJD(2) + t6 * qJD(3) + (-t730 / 0.4e1 - t722 / 0.4e1) * t779 + t707, t14 * qJD(1) + t544 * qJD(2) + t8 * qJD(3) + (-t721 / 0.4e1 - t729 / 0.4e1) * t777 + t707, t6 * qJD(1) + t8 * qJD(2) + t16 * qJD(4) + ((t118 * t182 + t156 * t166 + t337 * t372 - t339 * t374 + t619) * t771 + (t262 * t281 + (t393 * t527 - t395 * t529) * t471 + t91) * t772) * t776 + (t576 + t793) * qJD(3), t16 * qJD(3) + (m(6) * t579 + t569) * qJD(4) + t582 * t552, t582 * t136; t41 * qJD(2) + t690 - t804, t41 * qJD(1) + t690 - t815, (t337 * t529 + t339 * t527) * t704 + t249 + t582 * t97, t134 * t582 + t269 * t704 + t249, 0;];
Cq = t3;
