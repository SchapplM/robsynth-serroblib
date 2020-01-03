% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRP8_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP8_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP8_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP8_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:09
% EndTime: 2019-12-31 18:47:24
% DurationCPUTime: 11.71s
% Computational Cost: add. (18837->500), mult. (39318->634), div. (0->0), fcn. (46307->6), ass. (0->339)
t409 = sin(qJ(4));
t410 = cos(qJ(4));
t755 = (Icges(5,6) - Icges(6,6)) * t410 + (Icges(6,4) + Icges(5,5)) * t409;
t606 = sin(qJ(1));
t608 = cos(qJ(1));
t457 = -t606 * pkin(1) + t608 * qJ(2);
t433 = -t606 * pkin(2) + t457;
t605 = sin(qJ(3));
t607 = cos(qJ(3));
t381 = t608 * t605 - t606 * t607;
t380 = -t606 * t605 - t608 * t607;
t585 = t380 * pkin(7);
t491 = t381 * pkin(3) + t585;
t420 = t433 + t491;
t373 = t380 * rSges(6,2);
t535 = t409 * qJ(5);
t536 = t381 * t410;
t537 = t381 * t409;
t609 = rSges(6,1) + pkin(4);
t453 = rSges(6,3) * t537 + t381 * t535 + t536 * t609 + t373;
t141 = t420 + t453;
t423 = t585 + t433;
t708 = rSges(6,3) + qJ(5);
t749 = t609 * t410;
t434 = t409 * t708 + t749;
t424 = pkin(3) + t434;
t715 = t381 * t424 + t373;
t142 = t423 + t715;
t521 = t141 - t142;
t754 = m(6) * t521;
t445 = Icges(5,5) * t410 - Icges(5,6) * t409;
t263 = -Icges(5,3) * t380 - t445 * t381;
t447 = Icges(6,4) * t410 + Icges(6,6) * t409;
t266 = Icges(6,2) * t380 + t447 * t381;
t753 = t263 - t266;
t568 = Icges(6,5) * t409;
t450 = Icges(6,1) * t410 + t568;
t273 = -Icges(6,4) * t381 + t450 * t380;
t574 = Icges(5,4) * t409;
t452 = Icges(5,1) * t410 - t574;
t277 = -Icges(5,5) * t381 + t452 * t380;
t722 = t273 + t277;
t274 = Icges(6,4) * t380 + t450 * t381;
t278 = Icges(5,5) * t380 + t452 * t381;
t721 = t274 + t278;
t566 = Icges(5,2) * t410;
t448 = t566 + t574;
t437 = t409 * t448;
t573 = Icges(5,4) * t410;
t451 = Icges(5,1) * t409 + t573;
t523 = t410 * t451;
t567 = Icges(6,5) * t410;
t392 = -Icges(6,1) * t409 + t567;
t524 = t410 * t392;
t442 = -Icges(6,3) * t410 + t568;
t530 = t409 * t442;
t752 = t524 - t530 + t437 - t523;
t647 = m(6) / 0.2e1;
t648 = m(5) / 0.2e1;
t650 = m(4) / 0.2e1;
t165 = -t585 - t715;
t161 = t608 * t165;
t578 = t381 * rSges(6,2);
t584 = t381 * pkin(7);
t163 = t380 * t424 - t578 - t584;
t667 = t163 * t606 + t161;
t462 = -rSges(5,2) * t537 + t380 * rSges(5,3);
t579 = rSges(5,1) * t410;
t463 = pkin(3) + t579;
t659 = t463 * t381 + t462;
t210 = -t585 - t659;
t206 = t608 * t210;
t543 = t380 * t409;
t700 = rSges(5,2) * t543 + t381 * rSges(5,3);
t208 = t463 * t380 - t584 - t700;
t668 = t208 * t606 + t206;
t328 = t380 * rSges(4,1) + t381 * rSges(4,2);
t454 = -t381 * rSges(4,1) + t380 * rSges(4,2);
t669 = t328 * t606 + t608 * t454;
t477 = t667 * t647 + t668 * t648 + t669 * t650;
t470 = t381 * t606;
t338 = t409 * t470;
t471 = t380 * t608;
t458 = t409 * t471;
t307 = t338 + t458;
t200 = 0.2e1 * (-t458 / 0.4e1 - t338 / 0.4e1 + t307 / 0.4e1) * m(6);
t130 = t165 * t537;
t131 = t163 * t543;
t89 = -t130 + t131;
t751 = t89 * m(6) * qJD(3) - t200 * qJD(2);
t747 = t442 - t448;
t674 = t755 * t381;
t673 = t755 * t380;
t436 = t608 * pkin(1) + t606 * qJ(2);
t419 = t608 * pkin(2) + t436;
t143 = -t163 + t419;
t746 = -t142 * t163 + t143 * t165;
t703 = -t409 * rSges(6,3) - t535 - t749;
t292 = t380 * t703;
t289 = t381 * t703;
t744 = t721 * t410;
t743 = t722 * t410;
t264 = Icges(5,3) * t381 - t445 * t380;
t265 = -Icges(6,2) * t381 + t447 * t380;
t542 = t380 * t410;
t346 = Icges(6,5) * t542;
t560 = Icges(6,6) * t381;
t260 = -Icges(6,3) * t543 - t346 + t560;
t533 = t409 * t260;
t742 = -t381 * t533 + t722 * t536 + (-t264 + t265) * t380;
t345 = Icges(6,5) * t536;
t561 = Icges(6,6) * t380;
t259 = -Icges(6,3) * t537 - t345 - t561;
t534 = t409 * t259;
t741 = -t380 * t534 + t753 * t381 + t721 * t542;
t740 = t752 * t380 + t674;
t739 = t752 * t381 - t673;
t738 = t450 + t452;
t709 = m(6) * t410;
t449 = -Icges(5,2) * t409 + t573;
t269 = -Icges(5,6) * t381 + t449 * t380;
t270 = Icges(5,6) * t380 + t449 * t381;
t734 = (-t269 * t381 - t270 * t380) * t409 + t741 + t742;
t733 = m(6) * qJD(1);
t581 = m(6) * qJD(4);
t731 = t380 * t266;
t730 = t381 * t264;
t729 = t409 * t273;
t728 = t409 * t274;
t681 = t409 * t609;
t531 = t269 * t409;
t725 = -t381 * t531 + t742;
t532 = t270 * t409;
t724 = -t380 * t532 + t741;
t675 = -t264 - t532;
t723 = -t265 - t675;
t720 = -t747 * t380 - t722;
t719 = t747 * t381 + t721;
t718 = -Icges(6,1) * t543 - t451 * t380 - t260 - t269 + t346;
t174 = t269 * t410 + t409 * t277;
t553 = t174 * t380;
t170 = t260 * t410 + t729;
t555 = t170 * t380;
t677 = t555 / 0.2e1 + t553 / 0.2e1;
t173 = t270 * t410 + t409 * t278;
t554 = t173 * t381;
t169 = t259 * t410 + t728;
t556 = t169 * t381;
t678 = t556 / 0.2e1 + t554 / 0.2e1;
t717 = t677 - t678;
t181 = t423 + t659;
t182 = t419 - t208;
t716 = -t181 * t208 + t182 * t210;
t714 = -m(6) / 0.2e1;
t499 = t578 + t292;
t500 = -t373 + t289;
t98 = t499 * t380 + t500 * t381;
t713 = m(6) * t98;
t692 = -t380 / 0.2e1;
t691 = t380 / 0.2e1;
t712 = -t381 / 0.2e1;
t689 = t381 / 0.2e1;
t711 = t409 / 0.2e1;
t685 = t410 / 0.2e1;
t398 = -t409 * rSges(5,1) - rSges(5,2) * t410;
t591 = m(5) * t398;
t242 = t381 * t265;
t106 = t273 * t542 - t380 * t533 - t242;
t475 = -t274 * t536 + t381 * t534 - t731;
t704 = t106 - t475;
t443 = Icges(6,3) * t409 + t567;
t702 = (-t442 + t566 - t738) * t410 + (-t443 - t392 + t449 + t451 + t573) * t409;
t439 = rSges(5,1) * t536 + t462;
t478 = -m(5) * t668 / 0.2e1 - m(4) * t669 / 0.2e1 + t667 * t714;
t164 = -t453 - t491;
t209 = -t439 - t491;
t482 = (-t164 * t608 + t161) * t647 + (-t209 * t608 + t206) * t648;
t701 = t478 - t482;
t699 = t447 / 0.2e1 + t445 / 0.2e1;
t698 = t739 * t381 / 0.4e1;
t697 = t740 * t380 / 0.4e1;
t322 = t381 * t398;
t324 = t398 * t380;
t693 = t381 * t681 - t708 * t536;
t694 = t380 * t681 - t708 * t542;
t682 = (-(t143 - t163) * t694 - (-t142 + t165) * t693) * t647 + ((t182 - t208) * t324 + (-t181 + t210) * t322) * t648;
t490 = t410 * t708 - t681;
t461 = t490 * t380;
t679 = t490 * t381;
t696 = -t461 * t608 - t606 * t679;
t306 = -t328 + t419;
t415 = t433 - t454;
t695 = t306 * t454 - t415 * t328;
t587 = m(6) * t307;
t300 = -t587 / 0.2e1;
t201 = 0.2e1 * t300;
t686 = -t409 / 0.2e1;
t514 = t209 - t210;
t684 = m(5) * t514;
t518 = t164 - t165;
t683 = m(6) * t518;
t676 = -t263 + t531;
t672 = t720 * t409 + t718 * t410;
t506 = -t451 * t381 - t270;
t508 = Icges(6,1) * t537 + t259 - t345;
t671 = (-t506 + t508) * t410 + t719 * t409;
t513 = t380 * t263 - t278 * t536;
t103 = -t381 * t532 - t513;
t257 = t443 * t380 - t560;
t665 = -t242 * t689 + (t675 * t381 - t103 - t513) * t712 + (((t257 + t260) * t409 + t753) * t381 + (-t534 - t723 + t744) * t380 - t725) * t692;
t510 = t277 * t542 + t730;
t108 = -t380 * t531 + t510;
t258 = t443 * t381 + t561;
t379 = t380 ^ 2;
t664 = t379 * t266 / 0.2e1 + (t676 * t380 + t108 - t510 + t730) * t691 + ((-t265 + (-t258 - t259) * t409) * t380 + (t266 + t533 + t676 - t743) * t381 + t724) * t689;
t663 = (-t108 + t513 - t704 + t723 * t381 + (t743 + (t257 - t269) * t409) * t380) * t691 + (-t724 + t734) * t689;
t662 = (-t731 + t103 + t510 + t704 + (-t744 + (-t258 + t270) * t409) * t381) * t712 + (-t381 * t676 + t725 - t734) * t692;
t411 = (-t606 * t693 - t608 * t694) * t647 + (t606 * t322 + t324 * t608) * t648;
t661 = t679 * t683 / 0.2e1 + t322 * t684 / 0.2e1;
t520 = t143 + t163;
t660 = (t461 * t520 + (-t142 - t165) * t679) * t714 - ((-t181 - t210) * t381 + (t182 + t208) * t380) * t591 / 0.2e1;
t658 = t697 + t698;
t431 = t530 / 0.2e1 - t437 / 0.2e1 - t410 * t443 / 0.2e1 + t449 * t685 - t524 / 0.2e1 + t523 / 0.2e1 + t738 * t711;
t656 = t381 ^ 2;
t655 = 0.2e1 * qJD(1);
t654 = 0.4e1 * qJD(1);
t653 = 0.2e1 * qJD(3);
t652 = 0.4e1 * qJD(3);
t642 = m(5) * t716;
t638 = t208 * t684;
t637 = m(5) * (t181 * t322 - t182 * t324);
t636 = m(5) * (t208 * t324 - t210 * t322);
t96 = (t453 + t500) * t380;
t632 = t96 * t713;
t97 = (t434 * t380 + t499 - t578) * t381;
t631 = t97 * t713;
t123 = t142 * t537;
t124 = t143 * t543;
t87 = t123 - t124;
t630 = m(6) * (-t89 - t87);
t629 = m(6) * (t520 * t543 - t123 - t130);
t628 = m(6) * (-t131 + (t163 * t380 - t518 * t381) * t409);
t627 = t96 * t709;
t626 = t97 * t709;
t625 = m(6) * t746;
t125 = t142 * t542;
t476 = t693 * t543;
t551 = t694 * t409;
t623 = m(6) * (t476 - t125 + (-t143 * t410 - t551) * t381);
t622 = t163 * t683;
t132 = t165 * t542;
t621 = m(6) * (-t476 - t132 + (-t163 * t410 + t551) * t381);
t509 = -t461 * t537 + t543 * t679;
t620 = m(6) * (-t143 * t536 - t125 + t509);
t619 = m(6) * (-t142 * t693 + t143 * t694);
t618 = m(6) * (-t163 * t536 - t132 - t509);
t617 = m(6) * (-t163 * t694 + t165 * t693);
t128 = t142 * t608;
t616 = m(6) * (t143 * t606 + t128);
t604 = m(3) * ((t608 * rSges(3,3) + t457) * t608 + (t606 * rSges(3,3) + t436) * t606);
t603 = m(4) * t695;
t600 = m(4) * (t306 * t606 + t415 * t608);
t178 = t181 * t608;
t597 = m(5) * (t182 * t606 + t178);
t592 = (-t470 - t471) * t591;
t588 = m(6) * t696;
t580 = m(6) * qJD(5);
t549 = t322 * t398;
t548 = t324 * t398;
t527 = t409 * t410;
t254 = t592 / 0.2e1;
t412 = t647 * t696 + t254;
t73 = t411 + t412;
t522 = t73 * qJD(2);
t517 = t588 / 0.2e1 + t254;
t183 = t420 + t439;
t516 = t181 - t183;
t488 = qJD(4) * t380;
t487 = qJD(4) * t381;
t483 = t379 + t656;
t299 = t587 / 0.2e1;
t480 = t103 * t692 + t475 * t691 + t725 * t689;
t479 = t724 * t692 + (t106 + t108) * t689;
t441 = t461 * t693 - t679 * t694;
t435 = -t606 * t380 + t608 * t381;
t422 = t699 * t380 + t506 * t711 + t508 * t686 + t719 * t685 + t702 * t712 + t479;
t421 = t699 * t381 + t720 * t685 + t718 * t686 + t702 * t691 + t480;
t418 = -t431 - t717;
t417 = -t431 + t717;
t167 = t257 * t410 - t729;
t414 = -t556 / 0.4e1 - t554 / 0.4e1 + t658 - t660 - (t167 - t174 + t740) * t380 / 0.4e1 - t698;
t168 = t258 * t410 - t728;
t413 = -t555 / 0.4e1 - t553 / 0.4e1 + t658 - t661 - (t168 - t173 + t739) * t381 / 0.4e1 - t697;
t402 = t409 * rSges(5,2) - t579;
t339 = t379 * t527;
t296 = t339 + (-0.1e1 + t656) * t527;
t295 = t483 * t409;
t249 = t461 * t542;
t202 = t435 * t709;
t138 = 0.2e1 * t299;
t137 = t299 + t300;
t109 = t380 * t694 + t381 * t693;
t78 = t618 / 0.2e1;
t74 = -t411 + t412;
t72 = -t98 * t295 + t536 * t679 + t249;
t68 = t620 / 0.2e1;
t61 = t621 / 0.2e1;
t57 = t623 / 0.2e1;
t50 = t626 / 0.2e1;
t49 = t627 / 0.2e1;
t46 = t628 / 0.2e1;
t43 = t629 / 0.2e1;
t42 = t630 / 0.2e1;
t41 = t597 + t600 + t604 + t616;
t36 = t517 - t411;
t35 = -t588 / 0.2e1 - t592 / 0.2e1 + t411;
t34 = t517 + t411;
t31 = t431 + t617 + t636;
t28 = t431 + t619 + t637;
t27 = t622 + t638;
t26 = t603 + t625 + t642;
t17 = t78 + t61 - t626 / 0.2e1;
t16 = t78 + t50 - t621 / 0.2e1;
t15 = t61 + t50 - t618 / 0.2e1;
t14 = t68 + t57 - t627 / 0.2e1;
t13 = t68 + t49 - t623 / 0.2e1;
t12 = t57 + t49 - t620 / 0.2e1;
t11 = t43 + t46 - t630 / 0.2e1;
t10 = t42 + t46 - t629 / 0.2e1;
t9 = t42 + t43 - t628 / 0.2e1;
t8 = t477 - t701;
t7 = t478 + t482 - t477;
t6 = t477 + t701;
t5 = t413 + t414 + t431 - t682;
t4 = t414 + (t174 / 0.4e1 + t170 / 0.4e1) * t380 + (-t173 / 0.4e1 + t168 / 0.4e1) * t381 + t418 + t661 + t682;
t3 = t413 + (t173 / 0.4e1 + t169 / 0.4e1) * t381 + (-t174 / 0.4e1 + t167 / 0.4e1) * t380 + t417 + t660 + t682;
t2 = t631 + (-t480 + t663) * t381 + (-t479 + t665) * t380;
t1 = t632 + (t480 + t664) * t381 + (t479 + t662) * t380;
t18 = [(-m(5) * t516 * t182 / 0.4e1 + t143 * t754 / 0.4e1) * t654 + t41 * qJD(2) + t26 * qJD(3) + t28 * qJD(4) + t87 * t580, qJD(1) * t41 + qJD(3) * t6 + qJD(4) * t34 + qJD(5) * t137, t26 * qJD(1) + t6 * qJD(2) + t4 * qJD(4) + t9 * qJD(5) + (-t622 / 0.4e1 - t638 / 0.4e1) * t652 + (-t746 * t647 - t716 * t648 - t695 * t650) * t653, t28 * qJD(1) + t34 * qJD(2) + t4 * qJD(3) + (-t142 * t292 - t143 * t289 + t441) * t581 + t14 * qJD(5) - t632 * qJD(4) + (m(5) * (-t182 * t402 + t548) - t421 - t664) * t487 + (m(5) * (-t181 * t402 - t549) - t422 - t662) * t488, t137 * qJD(2) + t9 * qJD(3) + t14 * qJD(4) + t87 * t733; t7 * qJD(3) + t35 * qJD(4) + t138 * qJD(5) + (-t616 / 0.4e1 - t597 / 0.4e1 - t600 / 0.4e1 - t604 / 0.4e1) * t654 + ((-t608 * t141 + t128) * t647 + (-t608 * t183 + t178) * t648) * t655, 0, t7 * qJD(1) + t74 * qJD(4) + t201 * qJD(5) + t477 * t653, t35 * qJD(1) + t74 * qJD(3) + 0.2e1 * ((t289 * t608 - t292 * t606) * t647 + t435 * t402 * t648) * qJD(4) + t202 * qJD(5), qJD(1) * t138 + qJD(3) * t201 + qJD(4) * t202; t8 * qJD(2) + t27 * qJD(3) + t3 * qJD(4) + t10 * qJD(5) + (-t625 / 0.4e1 - t642 / 0.4e1 - t603 / 0.4e1) * t654 + ((t518 * t143 + t163 * t521) * t647 + (t514 * t182 - t208 * t516) * t648) * t655, qJD(1) * t8 - qJD(4) * t73 - qJD(5) * t200, t27 * qJD(1) + t31 * qJD(4) + t89 * t580, t3 * qJD(1) - t522 + t31 * qJD(3) + (-t163 * t289 - t165 * t292 - t441) * t581 + t17 * qJD(5) - t631 * qJD(4) + (m(5) * (-t208 * t402 - t548) + t421 - t663) * t487 + (m(5) * (-t210 * t402 + t549) + t422 - t665) * t488, t10 * qJD(1) + t17 * qJD(4) + t751; t36 * qJD(2) + t5 * qJD(3) + t1 * qJD(4) + t13 * qJD(5) + (-t619 / 0.4e1 - t637 / 0.4e1) * t654 + (-t679 * t754 + t418 + (t168 / 0.2e1 - t173 / 0.2e1 + t516 * t591) * t381 + t677) * qJD(1), qJD(1) * t36 + qJD(3) * t73, t5 * qJD(1) + t522 + t2 * qJD(4) + t16 * qJD(5) + (-t617 / 0.4e1 - t636 / 0.4e1) * t652 + (t417 + (-t174 / 0.2e1 + t167 / 0.2e1) * t380 + t678) * qJD(3), t1 * qJD(1) + t2 * qJD(3) + (m(5) * ((-t439 * t381 + (-rSges(5,1) * t542 + t700) * t380) * (-t322 * t381 - t324 * t380) + t483 * t402 * t398) + m(6) * (t109 * t98 + t289 * t679 + t292 * t461) + (t674 * t379 + (t672 * t381 + (t671 - t673) * t380) * t381) * t692 + (t673 * t656 + (t671 * t380 + (t672 - t674) * t381) * t380) * t689) * qJD(4) + t72 * t580, t13 * qJD(1) + t16 * qJD(3) + t72 * t581 + (t339 + (t656 * t409 - t295) * t410 - t296) * t580; t201 * qJD(2) + t11 * qJD(3) + t12 * qJD(4) + (-t124 + (t143 * t380 - t521 * t381) * t409 - t87) * t733, qJD(1) * t201 + qJD(3) * t200, t11 * qJD(1) + t15 * qJD(4) - t751, t12 * qJD(1) + t15 * qJD(3) + (t249 + (t381 * t679 + t109) * t410 + (t289 * t381 + t292 * t380 - t98) * t409 - t72) * t581 + t296 * t580, t296 * t581;];
Cq = t18;
