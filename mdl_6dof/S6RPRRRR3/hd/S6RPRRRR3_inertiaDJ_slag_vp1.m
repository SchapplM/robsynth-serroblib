% Calculate time derivative of joint inertia matrix for
% S6RPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR3_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR3_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR3_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR3_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:00:23
% EndTime: 2019-03-09 07:01:02
% DurationCPUTime: 22.34s
% Computational Cost: add. (98981->1146), mult. (82756->1518), div. (0->0), fcn. (79682->12), ass. (0->582)
t486 = sin(qJ(3));
t660 = qJD(1) * t486;
t779 = -t660 / 0.2e1;
t482 = qJ(1) + pkin(11);
t474 = cos(t482);
t489 = cos(qJ(3));
t656 = qJD(3) * t489;
t615 = t656 / 0.2e1;
t473 = sin(t482);
t634 = t473 * t660;
t778 = t474 * t615 - t634 / 0.2e1;
t661 = qJD(1) * t474;
t616 = t661 / 0.2e1;
t777 = t473 * t615 + t486 * t616;
t488 = cos(qJ(4));
t472 = t488 * pkin(4) + pkin(3);
t749 = pkin(3) - t472;
t618 = t749 * t489;
t751 = pkin(8) * t486;
t776 = t618 + t751;
t626 = t473 * t656;
t632 = t474 * t660;
t511 = t626 + t632;
t491 = -pkin(9) - pkin(8);
t654 = qJD(4) * t488;
t651 = pkin(4) * t654;
t775 = t491 * t660 + t651;
t761 = t473 / 0.2e1;
t657 = qJD(3) * t486;
t774 = -t657 / 0.2e1;
t745 = rSges(4,1) * t489;
t581 = -rSges(4,2) * t486 + t745;
t744 = rSges(4,3) * t474;
t386 = t473 * t581 - t744;
t703 = t474 * t489;
t705 = t474 * t486;
t770 = -rSges(4,2) * t705 + t473 * rSges(4,3);
t387 = rSges(4,1) * t703 + t770;
t454 = rSges(4,1) * t486 + rSges(4,2) * t489;
t520 = qJD(3) * t454;
t502 = rSges(4,2) * t634 + rSges(4,3) * t661 - t474 * t520;
t178 = (qJD(1) * t386 + t502) * t474 + (-t473 * t520 + (-t387 + t770) * qJD(1)) * t473;
t765 = 2 * m(4);
t773 = t178 * t765;
t740 = Icges(4,4) * t486;
t570 = Icges(4,1) * t489 - t740;
t381 = Icges(4,5) * t473 + t474 * t570;
t721 = t381 * t489;
t739 = Icges(4,4) * t489;
t565 = -Icges(4,2) * t486 + t739;
t379 = Icges(4,6) * t473 + t474 * t565;
t726 = t379 * t486;
t528 = -t721 + t726;
t772 = t474 * t528;
t771 = t486 * t749;
t480 = cos(qJ(1)) * pkin(1);
t769 = t473 * pkin(7) + t480;
t434 = t581 * qJD(3);
t768 = t454 * t434 * t765;
t560 = Icges(4,5) * t489 - Icges(4,6) * t486;
t376 = -Icges(4,3) * t474 + t473 * t560;
t378 = -Icges(4,6) * t474 + t473 * t565;
t380 = -Icges(4,5) * t474 + t473 * t570;
t481 = qJD(4) + qJD(5);
t475 = qJD(6) + t481;
t659 = qJD(1) * t489;
t605 = -t475 + t659;
t625 = t474 * t657;
t767 = t473 * t605 + t625;
t627 = t473 * t657;
t766 = t474 * t605 - t627;
t483 = -pkin(10) + t491;
t664 = t483 - t491;
t484 = qJ(4) + qJ(5);
t477 = cos(t484);
t438 = pkin(5) * t477 + t472;
t666 = t438 - t472;
t355 = t486 * t666 + t489 * t664;
t764 = 2 * m(5);
t763 = 2 * m(6);
t762 = 2 * m(7);
t760 = -t474 / 0.2e1;
t759 = t474 / 0.2e1;
t758 = -t489 / 0.2e1;
t757 = -rSges(5,3) - pkin(8);
t756 = m(4) * t454;
t755 = sin(qJ(1)) * pkin(1);
t754 = pkin(3) * t489;
t485 = sin(qJ(4));
t753 = pkin(4) * t485;
t752 = pkin(5) * t481;
t750 = qJD(1) / 0.2e1;
t748 = pkin(8) + t491;
t476 = sin(t484);
t735 = Icges(6,4) * t477;
t562 = -Icges(6,2) * t476 + t735;
t699 = t481 * t486;
t736 = Icges(6,4) * t476;
t327 = (-Icges(6,2) * t477 - t736) * t699 + (Icges(6,6) * t486 + t489 * t562) * qJD(3);
t567 = Icges(6,1) * t477 - t736;
t406 = -Icges(6,5) * t489 + t486 * t567;
t558 = Icges(6,5) * t477 - Icges(6,6) * t476;
t326 = (-Icges(6,5) * t476 - Icges(6,6) * t477) * t699 + (Icges(6,3) * t486 + t489 * t558) * qJD(3);
t328 = (-Icges(6,1) * t476 - t735) * t699 + (Icges(6,5) * t486 + t489 * t567) * qJD(3);
t404 = -Icges(6,3) * t489 + t486 * t558;
t405 = -Icges(6,6) * t489 + t486 * t562;
t719 = t406 * t477;
t499 = -t489 * t326 + t404 * t657 + t656 * t719 + (t328 * t486 - t405 * t699) * t477;
t629 = t405 * t656;
t117 = (-t629 + (-t406 * t481 - t327) * t486) * t476 + t499;
t478 = qJ(6) + t484;
t471 = cos(t478);
t470 = sin(t478);
t733 = Icges(7,4) * t471;
t561 = -Icges(7,2) * t470 + t733;
t702 = t475 * t486;
t734 = Icges(7,4) * t470;
t303 = (-Icges(7,2) * t471 - t734) * t702 + (Icges(7,6) * t486 + t489 * t561) * qJD(3);
t566 = Icges(7,1) * t471 - t734;
t390 = -Icges(7,5) * t489 + t486 * t566;
t557 = Icges(7,5) * t471 - Icges(7,6) * t470;
t302 = (-Icges(7,5) * t470 - Icges(7,6) * t471) * t702 + (Icges(7,3) * t486 + t489 * t557) * qJD(3);
t304 = (-Icges(7,1) * t470 - t733) * t702 + (Icges(7,5) * t486 + t489 * t566) * qJD(3);
t388 = -Icges(7,3) * t489 + t486 * t557;
t389 = -Icges(7,6) * t489 + t486 * t561;
t720 = t390 * t471;
t500 = -t489 * t302 + t388 * t657 + t656 * t720 + (t304 * t486 - t389 * t702) * t471;
t630 = t389 * t656;
t112 = (-t630 + (-t390 * t475 - t303) * t486) * t470 + t500;
t238 = -t388 * t489 + (-t389 * t470 + t720) * t486;
t231 = t238 * t657;
t606 = -t475 * t489 + qJD(1);
t527 = t473 * t606;
t257 = -t470 * t766 + t471 * t527;
t258 = t470 * t527 + t471 * t766;
t169 = Icges(7,5) * t258 + Icges(7,6) * t257 + Icges(7,3) * t511;
t171 = Icges(7,4) * t258 + Icges(7,2) * t257 + Icges(7,6) * t511;
t173 = Icges(7,1) * t258 + Icges(7,4) * t257 + Icges(7,5) * t511;
t709 = t473 * t489;
t382 = -t470 * t709 - t471 * t474;
t383 = -t470 * t474 + t471 * t709;
t711 = t473 * t486;
t281 = Icges(7,5) * t383 + Icges(7,6) * t382 + Icges(7,3) * t711;
t283 = Icges(7,4) * t383 + Icges(7,2) * t382 + Icges(7,6) * t711;
t285 = Icges(7,1) * t383 + Icges(7,4) * t382 + Icges(7,5) * t711;
t539 = -t283 * t470 + t285 * t471;
t46 = (qJD(3) * t539 - t169) * t489 + (qJD(3) * t281 + (-t283 * t475 + t173) * t471 + (-t285 * t475 - t171) * t470) * t486;
t526 = t474 * t606;
t255 = t470 * t767 + t471 * t526;
t256 = t470 * t526 - t471 * t767;
t624 = t474 * t656;
t510 = t624 - t634;
t168 = Icges(7,5) * t256 + Icges(7,6) * t255 + Icges(7,3) * t510;
t170 = Icges(7,4) * t256 + Icges(7,2) * t255 + Icges(7,6) * t510;
t172 = Icges(7,1) * t256 + Icges(7,4) * t255 + Icges(7,5) * t510;
t384 = -t470 * t703 + t471 * t473;
t385 = t470 * t473 + t471 * t703;
t282 = Icges(7,5) * t385 + Icges(7,6) * t384 + Icges(7,3) * t705;
t284 = Icges(7,4) * t385 + Icges(7,2) * t384 + Icges(7,6) * t705;
t286 = Icges(7,1) * t385 + Icges(7,4) * t384 + Icges(7,5) * t705;
t538 = -t284 * t470 + t286 * t471;
t47 = (qJD(3) * t538 - t168) * t489 + (qJD(3) * t282 + (-t284 * t475 + t172) * t471 + (-t286 * t475 - t170) * t470) * t486;
t150 = -t281 * t489 + t486 * t539;
t151 = -t282 * t489 + t486 * t538;
t547 = t150 * t473 + t151 * t474;
t548 = t150 * t474 - t151 * t473;
t13 = t231 + (qJD(3) * t547 - t112) * t489 + (qJD(1) * t548 + t46 * t473 + t47 * t474) * t486;
t244 = -t404 * t489 + (-t405 * t476 + t719) * t486;
t241 = t244 * t657;
t604 = -t481 * t489 + qJD(1);
t508 = t476 * t657 + t477 * t604;
t603 = -t481 + t659;
t708 = t474 * t476;
t273 = t473 * t508 - t603 * t708;
t507 = t476 * t604 - t477 * t657;
t707 = t474 * t477;
t274 = t473 * t507 + t603 * t707;
t184 = Icges(6,5) * t274 + Icges(6,6) * t273 + Icges(6,3) * t511;
t186 = Icges(6,4) * t274 + Icges(6,2) * t273 + Icges(6,6) * t511;
t188 = Icges(6,1) * t274 + Icges(6,4) * t273 + Icges(6,5) * t511;
t701 = t476 * t489;
t400 = -t473 * t701 - t707;
t700 = t477 * t489;
t401 = t473 * t700 - t708;
t292 = Icges(6,5) * t401 + Icges(6,6) * t400 + Icges(6,3) * t711;
t294 = Icges(6,4) * t401 + Icges(6,2) * t400 + Icges(6,6) * t711;
t296 = Icges(6,1) * t401 + Icges(6,4) * t400 + Icges(6,5) * t711;
t537 = -t294 * t476 + t296 * t477;
t52 = (qJD(3) * t537 - t184) * t489 + (qJD(3) * t292 + (-t294 * t481 + t188) * t477 + (-t296 * t481 - t186) * t476) * t486;
t714 = t473 * t476;
t271 = t474 * t508 + t603 * t714;
t713 = t473 * t477;
t272 = t474 * t507 - t603 * t713;
t183 = Icges(6,5) * t272 + Icges(6,6) * t271 + Icges(6,3) * t510;
t185 = Icges(6,4) * t272 + Icges(6,2) * t271 + Icges(6,6) * t510;
t187 = Icges(6,1) * t272 + Icges(6,4) * t271 + Icges(6,5) * t510;
t402 = -t474 * t701 + t713;
t403 = t474 * t700 + t714;
t293 = Icges(6,5) * t403 + Icges(6,6) * t402 + Icges(6,3) * t705;
t295 = Icges(6,4) * t403 + Icges(6,2) * t402 + Icges(6,6) * t705;
t297 = Icges(6,1) * t403 + Icges(6,4) * t402 + Icges(6,5) * t705;
t536 = -t295 * t476 + t297 * t477;
t53 = (qJD(3) * t536 - t183) * t489 + (qJD(3) * t293 + (-t295 * t481 + t187) * t477 + (-t297 * t481 - t185) * t476) * t486;
t160 = -t292 * t489 + t486 * t537;
t161 = -t293 * t489 + t486 * t536;
t541 = t160 * t473 + t161 * t474;
t542 = t160 * t474 - t161 * t473;
t747 = -t13 - t241 - (qJD(3) * t541 - t117) * t489 - (qJD(1) * t542 + t473 * t52 + t474 * t53) * t486;
t198 = t384 * t389 + t385 * t390 + t388 * t705;
t133 = t281 * t705 + t283 * t384 + t285 * t385;
t134 = t282 * t705 + t284 * t384 + t286 * t385;
t553 = t133 * t473 + t134 * t474;
t67 = -t198 * t489 + t486 * t553;
t221 = t402 * t405 + t403 * t406 + t404 * t705;
t138 = t292 * t705 + t294 * t402 + t296 * t403;
t139 = t293 * t705 + t295 * t402 + t297 * t403;
t549 = t138 * t473 + t139 * t474;
t77 = -t221 * t489 + t486 * t549;
t746 = -t67 - t77;
t743 = rSges(6,3) * t486;
t742 = rSges(7,3) * t486;
t741 = -rSges(7,3) + t483;
t738 = Icges(5,4) * t485;
t737 = Icges(5,4) * t488;
t697 = t485 * t489;
t704 = t474 * t488;
t412 = -t473 * t697 - t704;
t694 = t488 * t489;
t706 = t474 * t485;
t413 = t473 * t694 - t706;
t579 = -rSges(5,1) * t413 - rSges(5,2) * t412;
t324 = rSges(5,3) * t711 - t579;
t729 = t324 * t474;
t728 = t378 * t486;
t727 = t378 * t489;
t725 = t379 * t489;
t724 = t380 * t486;
t723 = t380 * t489;
t722 = t381 * t486;
t568 = Icges(5,1) * t488 - t738;
t418 = -Icges(5,5) * t489 + t486 * t568;
t718 = t418 * t488;
t652 = qJD(4) * t753;
t425 = -t476 * t752 - t652;
t717 = t425 * t489;
t426 = t477 * t752 + t651;
t716 = t426 * t474;
t440 = pkin(5) * t476 + t753;
t715 = t440 * t474;
t712 = t473 * t485;
t710 = t473 * t488;
t698 = t483 * t486;
t696 = t486 * t438;
t695 = t486 * t491;
t693 = t489 * t491;
t692 = -t117 - t112;
t645 = t256 * rSges(7,1) + t255 * rSges(7,2) + rSges(7,3) * t624;
t174 = -rSges(7,3) * t634 + t645;
t591 = t425 * t703 + t473 * t426 + t440 * t661 + t483 * t634;
t602 = t489 * t652;
t633 = t473 * t659;
t461 = pkin(4) * t706;
t637 = qJD(1) * t461 + t473 * t775;
t691 = -t666 * t633 + (-qJD(3) * t355 + t602) * t474 + t591 - t637 + t174;
t574 = rSges(7,1) * t258 + rSges(7,2) * t257;
t175 = rSges(7,3) * t511 + t574;
t573 = -rSges(7,1) * t383 - rSges(7,2) * t382;
t287 = rSges(7,3) * t711 - t573;
t690 = t175 * t705 + t287 * t624;
t577 = rSges(6,1) * t274 + rSges(6,2) * t273;
t191 = rSges(6,3) * t511 + t577;
t576 = -rSges(6,1) * t401 - rSges(6,2) * t400;
t299 = rSges(6,3) * t711 - t576;
t689 = t191 * t705 + t299 * t624;
t644 = t272 * rSges(6,1) + t271 * rSges(6,2) + rSges(6,3) * t624;
t190 = -rSges(6,3) * t634 + t644;
t449 = pkin(8) * t624;
t662 = qJD(1) * t473;
t216 = -t449 + t776 * t662 + (-t602 + (-t693 + t771) * qJD(3)) * t474 + t637;
t688 = -t190 - t216;
t447 = pkin(3) * t627;
t460 = pkin(4) * t712;
t571 = t472 * t627 + t473 * t602 + t474 * t775 + t491 * t626;
t217 = -pkin(8) * t626 + t447 + (-t474 * t776 + t460) * qJD(1) - t571;
t665 = t473 * t695 + t461;
t329 = -t473 * t776 - t665;
t687 = t217 * t705 + t329 * t624;
t608 = t666 * t489;
t519 = t608 - t698;
t251 = t473 * t519 + t665 - t715;
t263 = t287 * t705;
t686 = t251 * t705 + t263;
t685 = t251 + t287;
t607 = t664 * t486;
t667 = -t472 * t703 - t460;
t669 = t438 * t703 + t473 * t440;
t252 = -t474 * t607 + t667 + t669;
t288 = t385 * rSges(7,1) + t384 * rSges(7,2) + rSges(7,3) * t705;
t684 = t252 + t288;
t270 = (t425 + t652) * t486 + (t608 - t607) * qJD(3);
t572 = rSges(7,1) * t471 - rSges(7,2) * t470;
t306 = (-rSges(7,1) * t470 - rSges(7,2) * t471) * t702 + (t489 * t572 + t742) * qJD(3);
t683 = -t270 - t306;
t392 = -rSges(7,3) * t489 + t486 * t572;
t682 = t288 * t657 + t392 * t634;
t300 = t403 * rSges(6,1) + t402 * rSges(6,2) + rSges(6,3) * t705;
t575 = rSges(6,1) * t477 - rSges(6,2) * t476;
t407 = -rSges(6,3) * t489 + t486 * t575;
t681 = t300 * t657 + t407 * t634;
t680 = -t299 - t329;
t463 = pkin(3) * t703;
t421 = pkin(8) * t705 + t463;
t524 = -t474 * t695 - t667;
t330 = t524 - t421;
t679 = -t300 - t330;
t399 = t489 * t748 - t771;
t678 = t330 * t657 + t399 * t634;
t677 = t489 * t329 + t399 * t711;
t414 = -t474 * t697 + t710;
t415 = t474 * t694 + t712;
t325 = t415 * rSges(5,1) + t414 * rSges(5,2) + rSges(5,3) * t705;
t676 = -t325 - t421;
t337 = (-rSges(6,1) * t476 - rSges(6,2) * t477) * t699 + (t489 * t575 + t743) * qJD(3);
t655 = qJD(4) * t486;
t623 = t485 * t655;
t354 = -pkin(4) * t623 + (-t486 * t748 - t618) * qJD(3);
t675 = -t337 - t354;
t578 = rSges(5,1) * t488 - rSges(5,2) * t485;
t351 = (-rSges(5,1) * t485 - rSges(5,2) * t488) * t655 + (rSges(5,3) * t486 + t489 * t578) * qJD(3);
t584 = t751 + t754;
t439 = t584 * qJD(3);
t674 = -t351 - t439;
t673 = t355 + t392;
t232 = t489 * t287 + t392 * t711;
t239 = t489 * t299 + t407 * t711;
t459 = t486 * pkin(3) - t489 * pkin(8);
t424 = t459 * t662;
t672 = t399 * t662 + t424;
t420 = t584 * t473;
t671 = t473 * t420 + t474 * t421;
t670 = -t399 - t407;
t419 = -rSges(5,3) * t489 + t486 * t578;
t668 = -t419 - t459;
t377 = Icges(4,3) * t473 + t474 * t560;
t663 = qJD(1) * t377;
t658 = qJD(3) * t473;
t601 = -qJD(4) * t489 + qJD(1);
t506 = t485 * t657 + t488 * t601;
t600 = -qJD(4) + t659;
t313 = t474 * t506 + t600 * t712;
t505 = t485 * t601 - t488 * t657;
t314 = t474 * t505 - t600 * t710;
t559 = Icges(5,5) * t488 - Icges(5,6) * t485;
t348 = (-Icges(5,5) * t485 - Icges(5,6) * t488) * t655 + (Icges(5,3) * t486 + t489 * t559) * qJD(3);
t563 = -Icges(5,2) * t485 + t737;
t349 = (-Icges(5,2) * t488 - t738) * t655 + (Icges(5,6) * t486 + t489 * t563) * qJD(3);
t350 = (-Icges(5,1) * t485 - t737) * t655 + (Icges(5,5) * t486 + t489 * t568) * qJD(3);
t416 = -Icges(5,3) * t489 + t486 * t559;
t417 = -Icges(5,6) * t489 + t486 * t563;
t113 = t313 * t417 + t314 * t418 + t348 * t705 + t349 * t414 + t350 * t415 + t416 * t510;
t200 = Icges(5,5) * t314 + Icges(5,6) * t313 + Icges(5,3) * t510;
t202 = Icges(5,4) * t314 + Icges(5,2) * t313 + Icges(5,6) * t510;
t204 = Icges(5,1) * t314 + Icges(5,4) * t313 + Icges(5,5) * t510;
t319 = Icges(5,5) * t415 + Icges(5,6) * t414 + Icges(5,3) * t705;
t321 = Icges(5,4) * t415 + Icges(5,2) * t414 + Icges(5,6) * t705;
t323 = Icges(5,1) * t415 + Icges(5,4) * t414 + Icges(5,5) * t705;
t534 = -t321 * t485 + t323 * t488;
t59 = (qJD(3) * t534 - t200) * t489 + (qJD(3) * t319 - t202 * t485 + t204 * t488 + (-t321 * t488 - t323 * t485) * qJD(4)) * t486;
t650 = t113 / 0.2e1 + t59 / 0.2e1;
t315 = t473 * t506 - t600 * t706;
t316 = t473 * t505 + t600 * t704;
t114 = t315 * t417 + t316 * t418 + t348 * t711 + t349 * t412 + t350 * t413 + t416 * t511;
t201 = Icges(5,5) * t316 + Icges(5,6) * t315 + Icges(5,3) * t511;
t203 = Icges(5,4) * t316 + Icges(5,2) * t315 + Icges(5,6) * t511;
t205 = Icges(5,1) * t316 + Icges(5,4) * t315 + Icges(5,5) * t511;
t318 = Icges(5,5) * t413 + Icges(5,6) * t412 + Icges(5,3) * t711;
t320 = Icges(5,4) * t413 + Icges(5,2) * t412 + Icges(5,6) * t711;
t322 = Icges(5,1) * t413 + Icges(5,4) * t412 + Icges(5,5) * t711;
t535 = -t320 * t485 + t322 * t488;
t58 = (qJD(3) * t535 - t201) * t489 + (qJD(3) * t318 - t203 * t485 + t205 * t488 + (-t320 * t488 - t322 * t485) * qJD(4)) * t486;
t649 = t114 / 0.2e1 + t58 / 0.2e1;
t648 = -t216 - t691;
t647 = -t329 - t685;
t646 = -t330 - t684;
t643 = -t354 + t683;
t642 = t314 * rSges(5,1) + t313 * rSges(5,2) + rSges(5,3) * t624;
t641 = -t439 + t675;
t640 = t473 * (pkin(8) * t511 + qJD(1) * t463 - t447) + t474 * (-pkin(8) * t634 + t449 + (-t625 - t633) * pkin(3)) + t420 * t661;
t639 = -t399 - t673;
t638 = -t459 + t670;
t636 = t474 * pkin(2) + t769;
t635 = t419 * t662;
t628 = t417 * t656;
t622 = t711 / 0.2e1;
t621 = t705 / 0.2e1;
t179 = -t318 * t489 + t486 * t535;
t227 = t412 * t417 + t413 * t418 + t416 * t711;
t620 = t179 / 0.2e1 + t227 / 0.2e1;
t180 = -t319 * t489 + t486 * t534;
t228 = t414 * t417 + t415 * t418 + t416 * t705;
t619 = t180 / 0.2e1 + t228 / 0.2e1;
t617 = t662 / 0.2e1;
t468 = t474 * pkin(7);
t614 = t468 - t755;
t613 = -t438 * t489 - pkin(2);
t612 = t684 * t486;
t611 = t684 * t489;
t610 = t679 * t486;
t609 = t679 * t489;
t353 = t668 * t474;
t514 = (-t489 * t483 - t696) * qJD(3);
t149 = -t716 + (t514 + t717) * t473 + ((t440 - t753) * t473 + t519 * t474) * qJD(1) + t571;
t599 = t149 * t705 + t251 * t624 + t690;
t598 = t489 * t175 + t306 * t711 + t392 * t511;
t597 = t489 * t191 + t337 * t711 + t407 * t511;
t596 = t489 * t217 + t354 * t711 + t399 * t511;
t595 = t252 * t657 + t355 * t634 + t682;
t594 = -t439 + t643;
t593 = t473 * t329 + t474 * t330 + t671;
t140 = t489 * t251 + t355 * t711 + t232;
t592 = -t459 + t639;
t586 = t486 * t646;
t585 = t646 * t489;
t265 = t638 * t474;
t35 = t169 * t705 + t171 * t384 + t173 * t385 + t255 * t283 + t256 * t285 + t281 * t510;
t36 = t168 * t705 + t170 * t384 + t172 * t385 + t255 * t284 + t256 * t286 + t282 * t510;
t554 = t133 * t474 - t134 * t473;
t91 = t255 * t389 + t256 * t390 + t302 * t705 + t303 * t384 + t304 * t385 + t388 * t510;
t5 = (qJD(3) * t553 - t91) * t489 + (qJD(1) * t554 + qJD(3) * t198 + t35 * t473 + t36 * t474) * t486;
t197 = t382 * t389 + t383 * t390 + t388 * t711;
t37 = t169 * t711 + t171 * t382 + t173 * t383 + t257 * t283 + t258 * t285 + t281 * t511;
t38 = t168 * t711 + t170 * t382 + t172 * t383 + t257 * t284 + t258 * t286 + t282 * t511;
t131 = t281 * t711 + t283 * t382 + t285 * t383;
t132 = t282 * t711 + t284 * t382 + t286 * t383;
t555 = t131 * t473 + t132 * t474;
t556 = t131 * t474 - t132 * t473;
t92 = t257 * t389 + t258 * t390 + t302 * t711 + t303 * t382 + t304 * t383 + t388 * t511;
t6 = (qJD(3) * t555 - t92) * t489 + (qJD(1) * t556 + qJD(3) * t197 + t37 * t473 + t38 * t474) * t486;
t66 = -t197 * t489 + t486 * t555;
t583 = t5 * t705 + t6 * t711 + t67 * t624 + (-t238 * t489 + t486 * t547) * t657 + t511 * t66;
t580 = rSges(5,1) * t316 + rSges(5,2) * t315;
t564 = Icges(4,2) * t489 + t740;
t136 = t292 * t711 + t294 * t400 + t296 * t401;
t137 = t293 * t711 + t295 * t400 + t297 * t401;
t552 = t136 * t474 - t137 * t473;
t551 = t136 * t473 + t137 * t474;
t550 = t138 * t474 - t139 * t473;
t156 = t318 * t711 + t320 * t412 + t322 * t413;
t157 = t319 * t711 + t321 * t412 + t323 * t413;
t546 = t156 * t474 - t157 * t473;
t545 = t156 * t473 + t157 * t474;
t158 = t318 * t705 + t320 * t414 + t322 * t415;
t159 = t319 * t705 + t321 * t414 + t323 * t415;
t544 = t158 * t474 - t159 * t473;
t543 = t158 * t473 + t159 * t474;
t540 = t179 * t473 + t180 * t474;
t533 = -t325 * t473 + t729;
t532 = -t324 * t473 - t325 * t474;
t529 = -t723 + t728;
t224 = t592 * t474;
t525 = -pkin(2) - t581;
t523 = -t472 * t489 - pkin(2) - t743;
t522 = t474 * t216 + t473 * t217 + t329 * t661 + t640;
t518 = t486 * t757 - pkin(2) - t754;
t516 = qJD(3) * t564;
t515 = qJD(3) * (-Icges(4,5) * t486 - Icges(4,6) * t489);
t512 = t486 * t741 + t613;
t509 = t489 * t149 + t270 * t711 + t355 * t511 + t598;
t103 = t271 * t405 + t272 * t406 + t326 * t705 + t327 * t402 + t328 * t403 + t404 * t510;
t42 = t184 * t705 + t186 * t402 + t188 * t403 + t271 * t294 + t272 * t296 + t292 * t510;
t43 = t183 * t705 + t185 * t402 + t187 * t403 + t271 * t295 + t272 * t297 + t293 * t510;
t11 = (qJD(3) * t549 - t103) * t489 + (qJD(1) * t550 + qJD(3) * t221 + t42 * t473 + t43 * t474) * t486;
t104 = t273 * t405 + t274 * t406 + t326 * t711 + t327 * t400 + t328 * t401 + t404 * t511;
t220 = t400 * t405 + t401 * t406 + t404 * t711;
t44 = t184 * t711 + t186 * t400 + t188 * t401 + t273 * t294 + t274 * t296 + t292 * t511;
t45 = t183 * t711 + t185 * t400 + t187 * t401 + t273 * t295 + t274 * t297 + t293 * t511;
t12 = (qJD(3) * t551 - t104) * t489 + (qJD(1) * t552 + qJD(3) * t220 + t44 * t473 + t45 * t474) * t486;
t76 = -t220 * t489 + t486 * t551;
t504 = t11 * t705 + t12 * t711 + t77 * t624 + t583 + (-t244 * t489 + t486 * t541) * t657 + t511 * t76;
t503 = t473 * t523 - t755;
t19 = qJD(1) * t553 - t35 * t474 + t36 * t473;
t20 = qJD(1) * t555 - t37 * t474 + t38 * t473;
t501 = t19 * t621 + t20 * t622 + t5 * t761 + t6 * t760 + t548 * t774 + (qJD(1) * t547 - t46 * t474 + t47 * t473) * t758 + t66 * t617 + t67 * t616 - t778 * t554 - t777 * t556;
t498 = -t489 * t348 + t416 * t657 + t656 * t718 + (t350 * t488 - t417 * t654) * t486;
t497 = t473 * t518 - t755;
t496 = -t489 * t13 - t634 * t67 + t583;
t495 = t231 + (t46 + t92) * t622 + (t47 + t91) * t621 + (t151 + t198) * t778 + (t150 + t197) * t777;
t494 = t489 * t747 + t634 * t746 + t504;
t23 = qJD(1) * t549 - t42 * t474 + t43 * t473;
t24 = qJD(1) * t551 - t44 * t474 + t45 * t473;
t493 = t11 * t761 + t12 * t760 + t23 * t621 + t24 * t622 + t542 * t774 + (qJD(1) * t541 + t53 * t473 - t52 * t474) * t758 + t501 + t76 * t617 + t77 * t616 - t778 * t550 - t777 * t552;
t492 = t241 + t495 + (t104 + t52) * t622 + (t103 + t53) * t621 + (t161 + t221) * t778 + (t160 + t220) * t777;
t465 = pkin(7) * t661;
t431 = t560 * qJD(3);
t352 = t668 * t473;
t347 = t387 + t636;
t346 = t473 * t525 + t614 + t744;
t332 = t473 * t515 + t663;
t331 = -qJD(1) * t376 + t474 * t515;
t298 = t329 * t705;
t278 = t299 * t705;
t267 = t454 * t658 + (-t480 + (-rSges(4,3) - pkin(7)) * t473 + t525 * t474) * qJD(1);
t266 = t465 + (-t755 + (-pkin(2) - t745) * t473) * qJD(1) + t502;
t264 = t638 * t473;
t262 = -t416 * t489 + (-t417 * t485 + t718) * t486;
t259 = t262 * t657;
t250 = -t325 * t489 - t419 * t705;
t249 = t324 * t489 + t419 * t711;
t246 = t636 - t676;
t245 = t468 + t497 + t579;
t240 = -t300 * t489 - t407 * t705;
t237 = t377 * t473 - t772;
t236 = t376 * t473 - t474 * t529;
t235 = -t377 * t474 - t473 * t528;
t234 = -t376 * t474 - t473 * t529;
t233 = -t288 * t489 - t392 * t705;
t230 = qJD(1) * t353 + t473 * t674;
t229 = t474 * t674 + t424 + t635;
t226 = t524 + t636 + t300;
t225 = t468 + t503 + t576 + t665;
t223 = t592 * t473;
t222 = t533 * t486;
t219 = -t474 * t698 + t288 + t636 + t669;
t218 = t473 * t512 + t573 + t614 + t715;
t212 = -t300 * t711 + t278;
t207 = rSges(5,3) * t511 + t580;
t206 = -rSges(5,3) * t634 + t642;
t196 = -t288 * t711 + t263;
t189 = -t532 + t671;
t177 = t670 * t705 + t609;
t176 = t239 + t677;
t165 = qJD(1) * t265 + t473 * t641;
t164 = t407 * t662 + t474 * t641 + t672;
t163 = t447 + t757 * t626 + (t474 * t518 - t769) * qJD(1) - t580;
t162 = -pkin(3) * t625 + qJD(1) * t497 + t449 + t465 + t642;
t141 = -t673 * t705 - t611;
t135 = t473 * t610 + t278 + t298;
t130 = (-t628 + (-qJD(4) * t418 - t349) * t486) * t485 + t498;
t129 = -rSges(6,3) * t626 + (-t480 + (-pkin(7) - t753) * t473 + t523 * t474) * qJD(1) + t571 - t577;
t128 = t465 + (-t602 + (-t486 * t472 - t693) * qJD(3)) * t474 + t503 * qJD(1) + t637 + t644;
t127 = t299 * t473 + t300 * t474 + t593;
t126 = (t419 * t658 + t207) * t489 + (-qJD(3) * t324 + t351 * t473 + t419 * t661) * t486;
t125 = (-qJD(3) * t419 * t474 - t206) * t489 + (qJD(3) * t325 - t351 * t474 + t635) * t486;
t124 = t639 * t705 + t585;
t123 = t140 + t677;
t122 = qJD(1) * t224 + t473 * t594;
t121 = t474 * t594 + t662 * t673 + t672;
t120 = -t473 * t612 + t686;
t119 = t716 + (-t717 + (t489 * t741 + t696) * qJD(3)) * t473 + (-t480 + (-pkin(7) - t440) * t473 + t512 * t474) * qJD(1) - t574;
t118 = t465 + t474 * t514 + (-t755 + (t613 - t742) * t473) * qJD(1) + t591 + t645;
t116 = -t299 * t657 + t597;
t115 = -t190 * t489 + (-t337 * t486 - t407 * t656) * t474 + t681;
t111 = -t287 * t657 + t598;
t110 = -t174 * t489 + (-t306 * t486 - t392 * t656) * t474 + t682;
t106 = t473 * t586 + t298 + t686;
t100 = t473 * t685 + t474 * t684 + t593;
t82 = t533 * t656 + (qJD(1) * t532 - t206 * t473 + t207 * t474) * t486;
t81 = -t228 * t489 + t486 * t543;
t80 = -t227 * t489 + t486 * t545;
t79 = t206 * t474 + t207 * t473 + (t473 * t676 + t729) * qJD(1) + t640;
t75 = -t300 * t632 + (-t300 * t656 + (-qJD(1) * t299 - t190) * t486) * t473 + t689;
t69 = t657 * t680 + t596 + t597;
t68 = t688 * t489 + (t486 * t675 + t656 * t670) * t474 + t678 + t681;
t65 = -t288 * t632 + (-t288 * t656 + (-qJD(1) * t287 - t174) * t486) * t473 + t690;
t57 = t200 * t711 + t202 * t412 + t204 * t413 + t315 * t321 + t316 * t323 + t319 * t511;
t56 = t201 * t711 + t203 * t412 + t205 * t413 + t315 * t320 + t316 * t322 + t318 * t511;
t55 = t200 * t705 + t202 * t414 + t204 * t415 + t313 * t321 + t314 * t323 + t319 * t510;
t54 = t201 * t705 + t203 * t414 + t205 * t415 + t313 * t320 + t314 * t322 + t318 * t510;
t51 = -t657 * t685 + t509;
t50 = -t691 * t489 + (t486 * t683 - t656 * t673) * t474 + t595;
t39 = t190 * t474 + t191 * t473 + (t299 * t474 + (-t421 + t679) * t473) * qJD(1) + t522;
t34 = t610 * t661 + (qJD(3) * t609 + (qJD(1) * t680 + t688) * t486) * t473 + t687 + t689;
t33 = t647 * t657 + t509 + t596;
t32 = t648 * t489 + (t486 * t643 + t639 * t656) * t474 + t595 + t678;
t31 = -t612 * t661 + (-qJD(3) * t611 + (-qJD(1) * t685 - t691) * t486) * t473 + t599;
t30 = t691 * t474 + (t149 + t175) * t473 + (t685 * t474 + (-t421 + t646) * t473) * qJD(1) + t522;
t29 = t586 * t661 + (qJD(3) * t585 + (qJD(1) * t647 + t648) * t486) * t473 + t599 + t687;
t28 = qJD(1) * t545 + t473 * t57 - t474 * t56;
t27 = qJD(1) * t543 + t473 * t55 - t474 * t54;
t16 = (qJD(3) * t545 - t114) * t489 + (qJD(1) * t546 + qJD(3) * t227 + t473 * t56 + t474 * t57) * t486;
t15 = (qJD(3) * t543 - t113) * t489 + (qJD(1) * t544 + qJD(3) * t228 + t473 * t54 + t474 * t55) * t486;
t1 = [t499 + t500 + t498 + (t266 * t347 + t267 * t346) * t765 + (t162 * t246 + t163 * t245) * t764 + (t118 * t219 + t119 * t218) * t762 + (t128 * t226 + t129 * t225) * t763 - t418 * t623 + (-t564 + t570) * t657 + (Icges(4,1) * t486 + t565 + t739) * t656 + (-t349 * t486 - t628) * t485 + (-t327 * t486 - t406 * t699 - t629) * t476 + (-t303 * t486 - t390 * t702 - t630) * t470; 0; 0; m(5) * (t162 * t352 + t163 * t353 + t229 * t245 + t230 * t246) + m(6) * (t128 * t264 + t129 * t265 + t164 * t225 + t165 * t226) + m(7) * (t118 * t223 + t119 * t224 + t121 * t218 + t122 * t219) + ((qJD(1) * t379 - t473 * t516) * t758 + t381 * t779 - t52 / 0.2e1 - t46 / 0.2e1 - t104 / 0.2e1 - t92 / 0.2e1 + m(4) * (-t267 * t454 - t346 * t434) + t431 * t759 + (t728 / 0.2e1 - t723 / 0.2e1) * qJD(3) - t649) * t474 + ((-qJD(1) * t378 - t474 * t516) * t489 / 0.2e1 + t380 * t779 + t53 / 0.2e1 + t47 / 0.2e1 + t103 / 0.2e1 + t91 / 0.2e1 + m(4) * (-t266 * t454 - t347 * t434) + t431 * t761 + (-t726 / 0.2e1 + t721 / 0.2e1) * qJD(3) + t650) * t473 + ((t725 / 0.2e1 + t722 / 0.2e1 + t221 / 0.2e1 + t198 / 0.2e1 + t161 / 0.2e1 + t151 / 0.2e1 - t347 * t756 + t619) * t474 + (t727 / 0.2e1 + t724 / 0.2e1 + t346 * t756 + t220 / 0.2e1 + t197 / 0.2e1 + t160 / 0.2e1 + t150 / 0.2e1 + t620) * t473) * qJD(1); m(4) * t178 + m(5) * t79 + m(6) * t39 + m(7) * t30; (t100 * t30 + t121 * t224 + t122 * t223) * t762 + (t127 * t39 + t164 * t265 + t165 * t264) * t763 + (t189 * t79 + t229 * t353 + t230 * t352) * t764 + (t387 * t773 - t20 - t24 - t28 + (t768 - t235 * qJD(1) + (-qJD(1) * t529 - t332) * t474) * t474) * t474 + (t19 + t23 + t27 + t386 * t773 + (t236 * qJD(1) + t768 + (t528 * qJD(1) + t331) * t473) * t473 + ((-t332 + (-t722 - t725) * qJD(3) + t379 * t656 + t381 * t657 - t663) * t473 + (t378 * t656 + t380 * t657 + t331 - (t724 + t727) * qJD(3)) * t474 + (t237 - t234 + (t377 - t529) * t473 + t772) * qJD(1)) * t474) * t473 + (-t234 * t474 + t235 * t473 - t546 - t552 - t556) * t662 + (-t236 * t474 + t237 * t473 - t544 - t550 - t554) * t661; t492 + (t650 * t474 + t649 * t473 + (-t473 * t619 + t474 * t620) * qJD(1)) * t486 + m(5) * (t125 * t246 + t126 * t245 + t162 * t250 + t163 * t249) + m(7) * (t118 * t124 + t119 * t123 + t218 * t33 + t219 * t32) + t259 + (-t130 + (t473 * t620 + t474 * t619) * qJD(3) + t692) * t489 + m(6) * (t128 * t177 + t129 * t176 + t225 * t69 + t226 * t68); m(5) * t82 + m(6) * t34 + m(7) * t29; ((qJD(1) * t180 - t58) * t758 - t16 / 0.2e1 - t544 * t615 + t81 * t750) * t474 + ((qJD(1) * t179 + t59) * t758 + t15 / 0.2e1 - t546 * t615 + t80 * t750) * t473 + m(5) * (t125 * t352 + t126 * t353 + t189 * t82 + t222 * t79 + t229 * t249 + t230 * t250) + m(6) * (t127 * t34 + t135 * t39 + t164 * t176 + t165 * t177 + t264 * t68 + t265 * t69) + m(7) * (t100 * t29 + t106 * t30 + t121 * t123 + t122 * t124 + t223 * t32 + t224 * t33) + t493 + (qJD(3) * (-t179 * t474 + t180 * t473) / 0.2e1 + t27 * t759 + t28 * t761 + (t544 * t761 - t546 * t759) * qJD(1)) * t486; t504 + (t125 * t250 + t126 * t249 + t222 * t82) * t764 + (t135 * t34 + t176 * t69 + t177 * t68) * t763 + (t106 * t29 + t123 * t33 + t124 * t32) * t762 + (t130 * t489 - t259 + (t473 * t80 + t474 * t81 - t489 * t540) * qJD(3) + t747) * t489 + (-t489 * (t473 * t58 + t474 * t59) + t474 * t15 + t473 * t16 + (-t262 * t489 + t486 * t540) * qJD(3) + ((-t179 * t489 + t80) * t474 + (t180 * t489 + t746 - t81) * t473) * qJD(1)) * t486; t492 + m(7) * (t118 * t141 + t119 * t140 + t218 * t51 + t219 * t50) + t692 * t489 + m(6) * (t115 * t226 + t116 * t225 + t128 * t240 + t129 * t239); m(6) * t75 + m(7) * t31; t493 + m(7) * (t100 * t31 + t120 * t30 + t121 * t140 + t122 * t141 + t223 * t50 + t224 * t51) + m(6) * (t115 * t264 + t116 * t265 + t127 * t75 + t164 * t239 + t165 * t240 + t212 * t39); t494 + m(7) * (t106 * t31 + t120 * t29 + t123 * t51 + t124 * t50 + t140 * t33 + t141 * t32) + m(6) * (t115 * t177 + t116 * t176 + t135 * t75 + t212 * t34 + t239 * t69 + t240 * t68); t494 + (t120 * t31 + t140 * t51 + t141 * t50) * t762 + (t115 * t240 + t116 * t239 + t212 * t75) * t763; -t112 * t489 + t495 + m(7) * (t110 * t219 + t111 * t218 + t118 * t233 + t119 * t232); m(7) * t65; t501 + m(7) * (t100 * t65 + t110 * t223 + t111 * t224 + t121 * t232 + t122 * t233 + t196 * t30); m(7) * (t106 * t65 + t110 * t124 + t111 * t123 + t196 * t29 + t232 * t33 + t233 * t32) + t496; m(7) * (t110 * t141 + t111 * t140 + t120 * t65 + t196 * t31 + t232 * t51 + t233 * t50) + t496; (t110 * t233 + t111 * t232 + t196 * t65) * t762 + t496;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
