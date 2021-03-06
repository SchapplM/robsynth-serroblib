% Calculate time derivative of joint inertia matrix for
% S6RRPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR10_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR10_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR10_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR10_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR10_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:34:11
% EndTime: 2019-03-09 09:34:57
% DurationCPUTime: 30.36s
% Computational Cost: add. (40785->1146), mult. (49520->1574), div. (0->0), fcn. (46678->10), ass. (0->566)
t431 = cos(qJ(2));
t429 = sin(qJ(2));
t675 = Icges(4,6) * t429;
t685 = Icges(3,4) * t429;
t753 = t675 + t685 + (Icges(3,2) + Icges(4,3)) * t431;
t674 = Icges(4,6) * t431;
t684 = Icges(3,4) * t431;
t752 = t674 + t684 + (Icges(3,1) + Icges(4,2)) * t429;
t751 = t753 * qJD(2);
t750 = t752 * qJD(2);
t749 = -Icges(5,3) / 0.2e1;
t427 = cos(pkin(10));
t426 = sin(pkin(10));
t430 = sin(qJ(1));
t654 = t430 * t426;
t432 = cos(qJ(1));
t659 = t429 * t432;
t339 = t427 * t659 - t654;
t617 = qJD(2) * t431;
t585 = t430 * t617;
t272 = qJD(1) * t339 + t427 * t585;
t748 = -t272 / 0.2e1;
t598 = t426 * t659;
t653 = t430 * t427;
t340 = t598 + t653;
t273 = qJD(1) * t340 + t426 * t585;
t747 = -t273 / 0.2e1;
t703 = t430 / 0.2e1;
t746 = -t432 / 0.2e1;
t745 = -qJD(1) / 0.2e1;
t696 = qJD(1) / 0.2e1;
t619 = qJD(2) * t429;
t577 = -t619 / 0.2e1;
t621 = qJD(1) * t430;
t589 = t431 * t621;
t744 = t432 * t577 - t589 / 0.2e1;
t620 = qJD(1) * t432;
t578 = t620 / 0.2e1;
t743 = t430 * t577 + t431 * t578;
t616 = qJD(2) * t432;
t586 = t429 * t616;
t465 = t586 + t589;
t420 = pkin(10) + qJ(5);
t409 = cos(t420);
t408 = sin(t420);
t682 = Icges(6,4) * t408;
t514 = Icges(6,2) * t409 + t682;
t290 = Icges(6,6) * t429 - t431 * t514;
t681 = Icges(6,4) * t409;
t520 = Icges(6,1) * t408 + t681;
t291 = Icges(6,5) * t429 - t431 * t520;
t742 = t290 * t409 + t291 * t408;
t410 = qJ(6) + t420;
t402 = cos(t410);
t401 = sin(t410);
t680 = Icges(7,4) * t401;
t513 = Icges(7,2) * t402 + t680;
t281 = Icges(7,6) * t429 - t431 * t513;
t679 = Icges(7,4) * t402;
t519 = Icges(7,1) * t401 + t679;
t282 = Icges(7,5) * t429 - t431 * t519;
t741 = t281 * t402 + t282 * t401;
t625 = t430 ^ 2 + t432 ^ 2;
t740 = -0.1e1 + t625;
t739 = qJD(2) / 0.2e1;
t506 = -Icges(4,3) * t429 + t674;
t318 = Icges(4,5) * t430 - t432 * t506;
t508 = Icges(4,2) * t431 - t675;
t320 = Icges(4,4) * t430 - t432 * t508;
t488 = t318 * t429 - t320 * t431;
t738 = t430 * t488;
t518 = -Icges(3,2) * t429 + t684;
t315 = Icges(3,6) * t430 + t432 * t518;
t523 = Icges(3,1) * t431 - t685;
t317 = Icges(3,5) * t430 + t432 * t523;
t489 = t315 * t429 - t317 * t431;
t737 = t430 * t489;
t728 = Icges(4,4) * t432 + t430 * t508;
t729 = Icges(4,5) * t432 + t430 * t506;
t487 = -t429 * t729 + t431 * t728;
t736 = t432 * t487;
t314 = -Icges(3,6) * t432 + t430 * t518;
t316 = -Icges(3,5) * t432 + t430 * t523;
t490 = t314 * t429 - t316 * t431;
t735 = t432 * t490;
t428 = -pkin(8) - qJ(4);
t419 = -pkin(9) + t428;
t697 = pkin(4) * t426;
t361 = pkin(5) * t408 + t697;
t575 = -t361 + t697;
t547 = t575 * t429;
t462 = -t419 * t431 - t547;
t403 = t427 * pkin(4) + pkin(3);
t652 = t430 * t431;
t630 = t432 * t403 + t428 * t652;
t350 = pkin(5) * t409 + t403;
t663 = t350 * t432;
t173 = t430 * t462 + t630 - t663;
t657 = t430 * t402;
t300 = t401 * t432 + t429 * t657;
t658 = t430 * t401;
t301 = -t402 * t432 + t429 * t658;
t537 = -t301 * rSges(7,1) - t300 * rSges(7,2);
t207 = rSges(7,3) * t652 - t537;
t644 = t173 + t207;
t734 = t644 * t432;
t651 = t431 * t432;
t733 = t430 * rSges(4,1) - rSges(4,2) * t651;
t359 = t430 * pkin(3) + qJ(4) * t651;
t732 = -rSges(3,2) * t659 + t430 * rSges(3,3);
t407 = pkin(3) * t620;
t731 = -qJ(4) * t465 + t407;
t536 = rSges(7,1) * t401 + rSges(7,2) * t402;
t285 = rSges(7,3) * t429 - t431 * t536;
t268 = t285 * t652;
t149 = -t207 * t429 + t268;
t298 = t402 * t659 - t658;
t299 = t401 * t659 + t657;
t206 = t299 * rSges(7,1) + t298 * rSges(7,2) + rSges(7,3) * t651;
t192 = t429 * t206;
t150 = -t285 * t651 + t192;
t421 = qJD(5) + qJD(6);
t622 = qJD(1) * t429;
t570 = t421 + t622;
t584 = t431 * t616;
t446 = -t430 * t570 + t584;
t571 = t421 * t429 + qJD(1);
t486 = t401 * t571;
t180 = t402 * t446 - t432 * t486;
t485 = t402 * t571;
t181 = t401 * t446 + t432 * t485;
t643 = t181 * rSges(7,1) + t180 * rSges(7,2);
t116 = -rSges(7,3) * t465 + t643;
t618 = qJD(2) * t430;
t587 = t429 * t618;
t588 = t431 * t620;
t466 = -t587 + t588;
t447 = t432 * t570 + t585;
t178 = t402 * t447 - t430 * t486;
t179 = t401 * t447 + t430 * t485;
t538 = t179 * rSges(7,1) + t178 * rSges(7,2);
t115 = rSges(7,3) * t466 + t538;
t648 = t115 * t651 + t206 * t587;
t41 = -t207 * t586 + (-t116 * t430 + (-t206 * t432 - t430 * t207) * qJD(1)) * t431 + t648;
t730 = qJD(2) * (t149 * t432 + t150 * t430) - t41;
t512 = Icges(3,5) * t431 - Icges(3,6) * t429;
t312 = -Icges(3,3) * t432 + t430 * t512;
t516 = Icges(4,4) * t431 - Icges(4,5) * t429;
t727 = Icges(4,1) * t432 + t430 * t516;
t626 = t419 - t428;
t726 = t429 * t626 - t431 * t575;
t573 = t626 * t431;
t632 = -pkin(4) * t598 - t430 * t403;
t636 = t430 * t350 + t361 * t659;
t172 = -t432 * t573 + t632 + t636;
t645 = t172 + t206;
t725 = -t430 * t644 - t432 * t645;
t724 = 2 * m(3);
t723 = 2 * m(4);
t722 = 2 * m(5);
t721 = 2 * m(6);
t720 = 2 * m(7);
t719 = m(4) / 0.2e1;
t718 = -m(5) / 0.2e1;
t717 = m(5) / 0.2e1;
t716 = -m(6) / 0.2e1;
t715 = m(6) / 0.2e1;
t714 = -m(7) / 0.2e1;
t713 = m(7) / 0.2e1;
t515 = Icges(5,4) * t426 + Icges(5,2) * t427;
t305 = Icges(5,6) * t429 - t431 * t515;
t712 = t305 / 0.2e1;
t521 = Icges(5,1) * t426 + Icges(5,4) * t427;
t306 = Icges(5,5) * t429 - t431 * t521;
t711 = t306 / 0.2e1;
t710 = t339 / 0.2e1;
t709 = t340 / 0.2e1;
t708 = -t426 / 0.2e1;
t707 = t426 / 0.2e1;
t706 = -t427 / 0.2e1;
t705 = t427 / 0.2e1;
t704 = t429 / 0.2e1;
t702 = t432 / 0.2e1;
t701 = -rSges(6,3) - pkin(2);
t700 = -rSges(7,3) - pkin(2);
t378 = rSges(3,1) * t429 + rSges(3,2) * t431;
t699 = m(3) * t378;
t698 = pkin(2) * t431;
t695 = rSges(4,1) * t432;
t694 = rSges(4,2) * t429;
t693 = rSges(3,3) * t432;
t692 = rSges(5,3) * t429;
t691 = pkin(5) * qJD(5);
t690 = -rSges(4,3) - qJ(3);
t689 = rSges(7,3) - t419;
t509 = Icges(7,5) * t401 + Icges(7,6) * t402;
t280 = Icges(7,3) * t429 - t431 * t509;
t142 = t280 * t429 - t431 * t741;
t662 = t421 * t431;
t190 = (Icges(7,2) * t401 - t679) * t662 + (Icges(7,6) * t431 + t429 * t513) * qJD(2);
t189 = (-Icges(7,5) * t402 + Icges(7,6) * t401) * t662 + (Icges(7,3) * t431 + t429 * t509) * qJD(2);
t526 = t401 * t281 * t662 + t429 * t189 + t280 * t617 + t619 * t741;
t191 = (-Icges(7,1) * t402 + t680) * t662 + (Icges(7,5) * t431 + t429 * t519) * qJD(2);
t670 = t191 * t401;
t688 = t142 * t617 + ((-t670 + (-t282 * t421 - t190) * t402) * t431 + t526) * t429;
t510 = Icges(6,5) * t408 + Icges(6,6) * t409;
t289 = Icges(6,3) * t429 - t431 * t510;
t146 = t289 * t429 - t431 * t742;
t613 = qJD(5) * t431;
t221 = (Icges(6,2) * t408 - t681) * t613 + (Icges(6,6) * t431 + t429 * t514) * qJD(2);
t220 = (-Icges(6,5) * t409 + Icges(6,6) * t408) * t613 + (Icges(6,3) * t431 + t429 * t510) * qJD(2);
t525 = t408 * t290 * t613 + t429 * t220 + t289 * t617 + t619 * t742;
t222 = (-Icges(6,1) * t409 + t682) * t613 + (Icges(6,5) * t431 + t429 * t520) * qJD(2);
t668 = t222 * t408;
t687 = t146 * t617 + ((-t668 + (-qJD(5) * t291 - t221) * t409) * t431 + t525) * t429;
t672 = qJ(3) * t429;
t671 = qJ(3) * t431;
t655 = t430 * t409;
t310 = t408 * t432 + t429 * t655;
t656 = t430 * t408;
t311 = -t409 * t432 + t429 * t656;
t540 = -rSges(6,1) * t311 - rSges(6,2) * t310;
t219 = rSges(6,3) * t652 - t540;
t669 = t219 * t432;
t661 = t426 * t431;
t660 = t429 * t430;
t650 = -qJ(3) - t361;
t649 = -qJ(4) - t428;
t368 = t428 * t588;
t605 = t409 * t691;
t569 = t429 * t605;
t606 = t408 * t691;
t113 = t368 + (qJD(1) * t462 + t606) * t432 + (t569 + (t350 - t403) * qJD(1) + t726 * qJD(2)) * t430;
t647 = t113 + t115;
t524 = t350 * t620 + t361 * t584 + t419 * t465 + t432 * t569;
t558 = t426 * t584;
t562 = pkin(4) * t558 + t403 * t620 + t428 * t465;
t114 = (t575 * t622 - t606) * t430 + t524 - t562;
t646 = t114 + t116;
t195 = (-rSges(7,1) * t402 + rSges(7,2) * t401) * t662 + (rSges(7,3) * t431 + t429 * t536) * qJD(2);
t642 = t195 * t652 + t285 * t588;
t565 = qJD(5) + t622;
t443 = -t430 * t565 + t584;
t566 = qJD(5) * t429 + qJD(1);
t484 = t408 * t566;
t198 = t409 * t443 - t432 * t484;
t483 = t409 * t566;
t199 = t408 * t443 + t432 * t483;
t641 = t199 * rSges(6,1) + t198 * rSges(6,2);
t583 = t409 * t613;
t209 = -pkin(5) * t583 + (-t573 - t547) * qJD(2);
t640 = -t195 - t209;
t639 = -t726 + t285;
t341 = t426 * t432 + t429 * t653;
t274 = -qJD(1) * t341 + t427 * t584;
t342 = -t427 * t432 + t429 * t654;
t275 = -qJD(1) * t342 + t558;
t638 = t275 * rSges(5,1) + t274 * rSges(5,2);
t533 = t672 + t698;
t346 = t533 * t430;
t347 = pkin(2) * t651 + qJ(3) * t659;
t637 = t430 * t346 + t432 * t347;
t337 = qJD(2) * t533 - qJD(3) * t431;
t535 = -rSges(4,2) * t431 + rSges(4,3) * t429;
t635 = -t535 * qJD(2) - t337;
t634 = -t347 - t359;
t376 = pkin(2) * t429 - t671;
t348 = t376 * t621;
t591 = t429 * t621;
t633 = qJ(4) * t591 + t348;
t534 = rSges(4,3) * t431 + t694;
t631 = -t376 + t534;
t615 = qJD(3) * t429;
t629 = qJ(3) * t584 + t432 * t615;
t628 = rSges(3,2) * t591 + rSges(3,3) * t620;
t627 = t432 * pkin(1) + t430 * pkin(7);
t313 = Icges(3,3) * t430 + t432 * t512;
t624 = qJD(1) * t313;
t322 = Icges(4,1) * t430 - t432 * t516;
t623 = qJD(1) * t322;
t614 = qJD(4) * t431;
t612 = -rSges(5,3) - pkin(2) - qJ(4);
t610 = t429 * t697;
t609 = pkin(4) * t661;
t122 = Icges(6,5) * t199 + Icges(6,6) * t198 - Icges(6,3) * t465;
t124 = Icges(6,4) * t199 + Icges(6,2) * t198 - Icges(6,6) * t465;
t126 = Icges(6,1) * t199 + Icges(6,4) * t198 - Icges(6,5) * t465;
t308 = t409 * t659 - t656;
t309 = t408 * t659 + t655;
t212 = Icges(6,5) * t309 + Icges(6,6) * t308 + Icges(6,3) * t651;
t214 = Icges(6,4) * t309 + Icges(6,2) * t308 + Icges(6,6) * t651;
t216 = Icges(6,1) * t309 + Icges(6,4) * t308 + Icges(6,5) * t651;
t498 = t214 * t409 + t216 * t408;
t34 = (qJD(2) * t498 + t122) * t429 + (qJD(2) * t212 - t124 * t409 - t126 * t408 + (t214 * t408 - t216 * t409) * qJD(5)) * t431;
t58 = t198 * t290 + t199 * t291 + t220 * t651 + t308 * t221 + t309 * t222 - t289 * t465;
t604 = t34 / 0.2e1 + t58 / 0.2e1;
t444 = t432 * t565 + t585;
t196 = t409 * t444 - t430 * t484;
t197 = t408 * t444 + t430 * t483;
t121 = Icges(6,5) * t197 + Icges(6,6) * t196 + Icges(6,3) * t466;
t123 = Icges(6,4) * t197 + Icges(6,2) * t196 + Icges(6,6) * t466;
t125 = Icges(6,1) * t197 + Icges(6,4) * t196 + Icges(6,5) * t466;
t213 = Icges(6,5) * t311 + Icges(6,6) * t310 + Icges(6,3) * t652;
t215 = Icges(6,4) * t311 + Icges(6,2) * t310 + Icges(6,6) * t652;
t217 = Icges(6,1) * t311 + Icges(6,4) * t310 + Icges(6,5) * t652;
t497 = t215 * t409 + t217 * t408;
t35 = (qJD(2) * t497 + t121) * t429 + (qJD(2) * t213 - t123 * t409 - t125 * t408 + (t215 * t408 - t217 * t409) * qJD(5)) * t431;
t57 = t196 * t290 + t197 * t291 + t220 * t652 + t310 * t221 + t311 * t222 + t289 * t466;
t603 = t35 / 0.2e1 + t57 / 0.2e1;
t229 = Icges(5,5) * t340 + Icges(5,6) * t339 + Icges(5,3) * t651;
t602 = t229 * t652;
t601 = t229 * t651;
t230 = Icges(5,5) * t342 + Icges(5,6) * t341 + Icges(5,3) * t652;
t600 = t230 * t652;
t599 = t230 * t651;
t134 = t289 * t651 + t308 * t290 + t309 * t291;
t96 = t212 * t429 - t431 * t498;
t597 = -t134 / 0.2e1 - t96 / 0.2e1;
t135 = t289 * t652 + t290 * t310 + t291 * t311;
t97 = t213 * t429 - t431 * t497;
t596 = -t97 / 0.2e1 - t135 / 0.2e1;
t389 = pkin(2) * t587;
t595 = t430 * (pkin(2) * t588 + t430 * t615 - t389 + (t429 * t620 + t585) * qJ(3)) + t432 * (-pkin(2) * t465 - qJ(3) * t591 + t629) + t346 * t620;
t480 = -t428 * t651 - t632;
t253 = t480 - t359;
t594 = -t253 + t634;
t330 = t429 * t649 - t609;
t593 = t330 * t621 + t633;
t218 = t309 * rSges(6,1) + t308 * rSges(6,2) + rSges(6,3) * t651;
t237 = t340 * rSges(5,1) + t339 * rSges(5,2) + rSges(5,3) * t651;
t406 = pkin(7) * t620;
t592 = t406 + t629;
t539 = rSges(6,1) * t408 + rSges(6,2) * t409;
t295 = rSges(6,3) * t429 - t431 * t539;
t590 = t295 * t621;
t581 = t652 / 0.2e1;
t580 = t651 / 0.2e1;
t579 = -t516 * qJD(2) / 0.2e1 + t512 * t739;
t576 = -qJ(3) - t697;
t86 = t212 * t652 + t214 * t310 + t216 * t311;
t87 = t213 * t652 + t215 * t310 + t217 * t311;
t529 = t430 * t87 + t432 * t86;
t43 = t135 * t429 + t431 * t529;
t574 = t429 * t97 + t43;
t284 = t631 * t432;
t572 = -qJ(4) * t429 - t376;
t568 = t717 + t715 + t713;
t108 = Icges(7,5) * t181 + Icges(7,6) * t180 - Icges(7,3) * t465;
t110 = Icges(7,4) * t181 + Icges(7,2) * t180 - Icges(7,6) * t465;
t112 = Icges(7,1) * t181 + Icges(7,4) * t180 - Icges(7,5) * t465;
t200 = Icges(7,5) * t299 + Icges(7,6) * t298 + Icges(7,3) * t651;
t202 = Icges(7,4) * t299 + Icges(7,2) * t298 + Icges(7,6) * t651;
t204 = Icges(7,1) * t299 + Icges(7,4) * t298 + Icges(7,5) * t651;
t500 = t202 * t402 + t204 * t401;
t26 = (qJD(2) * t500 + t108) * t429 + (qJD(2) * t200 + (-t204 * t421 - t110) * t402 + (t202 * t421 - t112) * t401) * t431;
t107 = Icges(7,5) * t179 + Icges(7,6) * t178 + Icges(7,3) * t466;
t109 = Icges(7,4) * t179 + Icges(7,2) * t178 + Icges(7,6) * t466;
t111 = Icges(7,1) * t179 + Icges(7,4) * t178 + Icges(7,5) * t466;
t201 = Icges(7,5) * t301 + Icges(7,6) * t300 + Icges(7,3) * t652;
t203 = Icges(7,4) * t301 + Icges(7,2) * t300 + Icges(7,6) * t652;
t205 = Icges(7,1) * t301 + Icges(7,4) * t300 + Icges(7,5) * t652;
t499 = t203 * t402 + t205 * t401;
t27 = (qJD(2) * t499 + t107) * t429 + (qJD(2) * t201 + (-t205 * t421 - t109) * t402 + (t203 * t421 - t111) * t401) * t431;
t130 = t280 * t652 + t281 * t300 + t282 * t301;
t80 = t200 * t652 + t202 * t300 + t204 * t301;
t81 = t201 * t652 + t203 * t300 + t205 * t301;
t531 = t430 * t81 + t432 * t80;
t40 = t130 * t429 + t431 * t531;
t19 = t108 * t652 + t300 * t110 + t301 * t112 + t178 * t202 + t179 * t204 + t200 * t466;
t20 = t107 * t652 + t300 * t109 + t301 * t111 + t178 * t203 + t179 * t205 + t201 * t466;
t49 = t178 * t281 + t179 * t282 + t189 * t652 + t300 * t190 + t301 * t191 + t280 * t466;
t56 = t80 * t430 - t432 * t81;
t5 = (-qJD(2) * t531 + t49) * t429 + (-qJD(1) * t56 + qJD(2) * t130 + t19 * t432 + t20 * t430) * t431;
t92 = t200 * t429 - t431 * t500;
t93 = t201 * t429 - t431 * t499;
t527 = t92 * t430 - t432 * t93;
t528 = t93 * t430 + t92 * t432;
t129 = t280 * t651 + t298 * t281 + t299 * t282;
t21 = t108 * t651 + t298 * t110 + t299 * t112 + t180 * t202 + t181 * t204 - t200 * t465;
t22 = t107 * t651 + t298 * t109 + t299 * t111 + t180 * t203 + t181 * t205 - t201 * t465;
t50 = t180 * t281 + t181 * t282 + t189 * t651 + t298 * t190 + t299 * t191 - t280 * t465;
t78 = t200 * t651 + t298 * t202 + t299 * t204;
t79 = t201 * t651 + t298 * t203 + t299 * t205;
t532 = t430 * t79 + t432 * t78;
t55 = t78 * t430 - t432 * t79;
t6 = (-qJD(2) * t532 + t50) * t429 + (-qJD(1) * t55 + qJD(2) * t129 + t21 * t432 + t22 * t430) * t431;
t567 = t5 * t652 + t6 * t651 + t40 * t588 + (t142 * t429 + t431 * t528) * t617 + t429 * (-t528 * t619 + (-qJD(1) * t527 + t26 * t432 + t27 * t430) * t431 + t688);
t564 = t429 * t116 + t206 * t617 + t285 * t465;
t417 = t432 * pkin(3);
t360 = qJ(4) * t652 - t417;
t563 = t432 * t359 + t430 * t360 + t637;
t391 = t432 * t614;
t561 = t391 + t592;
t560 = rSges(4,1) * t620 + rSges(4,2) * t465 + rSges(4,3) * t584;
t559 = t627 + t347;
t39 = t129 * t429 + t431 * t532;
t84 = t212 * t651 + t308 * t214 + t309 * t216;
t85 = t213 * t651 + t308 * t215 + t309 * t217;
t530 = t430 * t85 + t432 * t84;
t42 = t134 * t429 + t431 * t530;
t557 = -t429 * t96 - t39 - t42;
t542 = rSges(5,1) * t426 + rSges(5,2) * t427;
t307 = -t431 * t542 + t692;
t552 = -t307 + t572;
t551 = -t330 + t572;
t550 = t429 * t690 - pkin(1);
t549 = t229 / 0.2e1 + t317 / 0.2e1 - t320 / 0.2e1;
t548 = t728 / 0.2e1 + t316 / 0.2e1 + t230 / 0.2e1;
t546 = t429 * t650 - pkin(1);
t545 = rSges(3,1) * t431 - rSges(3,2) * t429;
t544 = -t273 * rSges(5,1) - t272 * rSges(5,2);
t543 = -rSges(5,1) * t342 - rSges(5,2) * t341;
t541 = t197 * rSges(6,1) + t196 * rSges(6,2);
t60 = t84 * t430 - t432 * t85;
t61 = t86 * t430 - t432 * t87;
t511 = Icges(5,5) * t426 + Icges(5,6) * t427;
t416 = t432 * pkin(7);
t445 = (-pkin(2) - t689) * t431 + t546;
t132 = t430 * t445 + t416 + t537 + t663;
t133 = -t419 * t651 + t206 + t559 + t636;
t504 = t132 * t432 + t133 * t430;
t442 = t429 * t576 + t431 * t701 - pkin(1);
t143 = t430 * t442 + t416 + t540 + t630;
t144 = t480 + t559 + t218;
t503 = t143 * t432 + t144 * t430;
t461 = t431 * t612 - pkin(1) - t672;
t448 = t461 * t430;
t162 = t416 + t417 + t448 + t543;
t163 = t559 + t237 + t359;
t501 = t162 * t432 + t163 * t430;
t496 = -t218 * t432 - t430 * t219;
t495 = t218 * t430 - t669;
t482 = -t614 - t615;
t326 = rSges(3,1) * t651 + t732;
t327 = rSges(4,3) * t659 + t733;
t481 = -pkin(1) - t545;
t479 = -t295 + t551;
t478 = -qJ(4) * t431 + t610;
t382 = qJ(4) * t587;
t477 = t430 * (qJD(1) * t359 + t430 * t614 - t382) + t432 * (t391 + t731) + t360 * t620 + t595;
t254 = t430 * t478 + t417 - t630;
t476 = t432 * t253 + t430 * t254 + t563;
t228 = t552 * t432;
t474 = qJD(2) * t378;
t473 = t551 - t639;
t470 = qJD(2) * (Icges(4,4) * t429 + Icges(4,5) * t431);
t469 = qJD(2) * (-Icges(3,5) * t429 - Icges(3,6) * t431);
t175 = t479 * t432;
t464 = -qJ(4) * t617 - qJD(4) * t429 - t337;
t18 = (t172 * t430 - t734) * t619 + (qJD(1) * t725 + t113 * t432 - t646 * t430) * t431 + t648;
t82 = -t429 * t644 - t652 * t726 + t268;
t83 = t429 * t172 - t639 * t651 + t192;
t460 = t616 * t82 + t618 * t83 - t18;
t140 = t473 * t432;
t139 = t473 * t430;
t452 = t430 * (-t368 + t382 + (t428 * t429 + t609) * t618 + (t478 * t432 + (-pkin(3) + t403) * t430) * qJD(1)) + t432 * (-t591 * t697 + t562 - t731) + t254 * t620 + t477;
t17 = t646 * t432 + t647 * t430 + (t734 + (t594 - t645) * t430) * qJD(1) + t452;
t459 = t139 * t618 + t140 * t616 - t17;
t151 = -t219 * t429 + t295 * t652;
t152 = t429 * t218 - t295 * t651;
t127 = rSges(6,3) * t466 + t541;
t128 = -rSges(6,3) * t465 + t641;
t44 = t495 * t619 + (qJD(1) * t496 + t127 * t432 - t128 * t430) * t431;
t458 = t151 * t616 + t152 * t618 - t44;
t174 = t479 * t430;
t23 = t430 * t127 + t128 * t432 + (t669 + (-t218 + t594) * t430) * qJD(1) + t452;
t457 = t174 * t618 + t175 * t616 - t23;
t227 = t552 * t430;
t238 = rSges(5,3) * t652 - t543;
t48 = t430 * (-rSges(5,3) * t587 - t544) + t432 * (-rSges(5,3) * t586 + t638) + (t432 * t238 + (-t237 + t634) * t430) * qJD(1) + t477;
t456 = t227 * t618 + t228 * t616 - t48;
t455 = (rSges(4,2) - pkin(2)) * t431 + t550;
t454 = -(rSges(5,3) * t431 + t429 * t542) * qJD(2) + t464;
t453 = -(t431 * t649 + t610) * qJD(2) + t464;
t231 = Icges(5,4) * t340 + Icges(5,2) * t339 + Icges(5,6) * t651;
t233 = Icges(5,1) * t340 + Icges(5,4) * t339 + Icges(5,5) * t651;
t451 = t231 * t705 + t233 * t707 - t315 / 0.2e1 + t318 / 0.2e1;
t232 = Icges(5,4) * t342 + Icges(5,2) * t341 + Icges(5,6) * t652;
t234 = Icges(5,1) * t342 + Icges(5,4) * t341 + Icges(5,5) * t652;
t450 = t232 * t706 + t234 * t708 + t314 / 0.2e1 + t729 / 0.2e1;
t223 = (-rSges(6,1) * t409 + rSges(6,2) * t408) * t613 + (rSges(6,3) * t431 + t429 * t539) * qJD(2);
t449 = -t223 + t453;
t12 = qJD(1) * t531 + t19 * t430 - t20 * t432;
t13 = qJD(1) * t532 + t21 * t430 - t22 * t432;
t441 = t6 * t703 + t12 * t581 + t13 * t580 + (qJD(1) * t528 + t26 * t430 - t27 * t432) * t704 + t40 * t621 / 0.2e1 + t39 * t578 + t5 * t746 + t527 * t617 / 0.2e1 + t743 * t56 + t744 * t55;
t440 = t453 + t640;
t439 = qJD(1) * t442;
t438 = t688 + (t27 + t49) * t581 + (t26 + t50) * t580 + (t129 + t92) * t744 + (t130 + t93) * t743;
t187 = t207 * t651;
t131 = -t206 * t652 + t187;
t64 = -t195 * t651 + t564;
t65 = -t429 * t115 + (-t207 * t431 - t285 * t660) * qJD(2) + t642;
t437 = qJD(2) * t131 + t430 * t64 + t432 * t65 + (-t149 * t430 + t150 * t432) * qJD(1);
t436 = (-t39 * t432 - t40 * t430) * t619 - t39 * t589 + t567;
t136 = t495 * t431;
t32 = (-t616 * t726 + t114) * t429 + (qJD(2) * t172 + t432 * t640 - t621 * t726) * t431 + t564;
t33 = (t430 * t209 - t620 * t726) * t431 - t647 * t429 + (-t431 * t644 - t639 * t660) * qJD(2) + t642;
t69 = (t295 * t616 + t128) * t429 + (qJD(2) * t218 - t223 * t432 + t590) * t431;
t70 = (-t295 * t618 - t127) * t429 + (-qJD(2) * t219 + t430 * t223 + t295 * t620) * t431;
t73 = t187 + (t173 * t432 - t430 * t645) * t431;
t435 = (-qJD(2) * t136 - t151 * t621 + t152 * t620 + t430 * t69 + t432 * t70) * t715 + (qJD(2) * t73 + t32 * t430 + t33 * t432 + t620 * t83 - t621 * t82) * t713;
t102 = qJD(1) * t448 + t586 * t612 + t407 + t561 + t638;
t103 = t382 + t389 + ((-t671 + t692) * qJD(2) + t482) * t430 + ((-pkin(3) - pkin(7)) * t430 + t461 * t432) * qJD(1) + t544;
t67 = t700 * t586 + (-t606 + (t431 * t700 + t546) * qJD(1)) * t430 + t524 + t561 + t643;
t68 = t389 + (qJD(1) * t445 - t606) * t432 + (-t614 + (-qJD(3) - t605) * t429 + (-pkin(7) - t350) * qJD(1) + (t429 * t689 + t431 * t650) * qJD(2)) * t430 - t538;
t76 = t430 * t439 + t586 * t701 + t561 + t562 + t641;
t77 = t368 + t389 + t432 * t439 + ((-pkin(7) - t403) * qJD(1) + (t576 * t431 + (rSges(6,3) - t428) * t429) * qJD(2) + t482) * t430 - t541;
t434 = (-t132 * t621 + t133 * t620 + t430 * t67 + t432 * t68) * t713 + (-t143 * t621 + t144 * t620 + t430 * t76 + t432 * t77) * t715 + (t102 * t430 + t103 * t432 - t162 * t621 + t163 * t620) * t717;
t105 = t237 * t432 + t430 * t238 + t563;
t137 = t307 * t621 + t432 * t454 + t633;
t138 = qJD(1) * t228 + t430 * t454;
t59 = t476 - t725;
t71 = t432 * t440 + t621 * t639 + t593;
t72 = qJD(1) * t140 + t430 * t440;
t75 = t476 - t496;
t94 = t432 * t449 + t590 + t593;
t95 = qJD(1) * t175 + t430 * t449;
t433 = (qJD(2) * t59 + t139 * t620 - t140 * t621 + t430 * t72 + t432 * t71) * t713 + (qJD(2) * t75 + t174 * t620 - t175 * t621 + t430 * t95 + t432 * t94) * t715 + (qJD(2) * t105 + t137 * t432 + t138 * t430 + t227 * t620 - t228 * t621) * t717;
t358 = t545 * qJD(2);
t328 = t430 * t535 - t695;
t325 = t430 * t545 - t693;
t288 = (Icges(5,5) * t431 + t429 * t521) * qJD(2);
t287 = (Icges(5,6) * t431 + t429 * t515) * qJD(2);
t283 = t631 * t430;
t279 = t326 + t627;
t278 = t430 * t481 + t416 + t693;
t277 = t740 * t429 * t617;
t252 = qJD(1) * t727 + t432 * t470;
t251 = t430 * t470 + t623;
t242 = t430 * t469 + t624;
t241 = -qJD(1) * t312 + t432 * t469;
t240 = t327 + t559;
t239 = t430 * t455 + t416 + t695;
t211 = t378 * t618 + ((-rSges(3,3) - pkin(7)) * t430 + t481 * t432) * qJD(1);
t210 = -rSges(3,1) * t465 - rSges(3,2) * t584 - pkin(1) * t621 + t406 + t628;
t183 = qJD(1) * t284 + t430 * t635;
t182 = t432 * t635 - t534 * t621 + t348;
t171 = t430 * t487 + t432 * t727;
t170 = -t322 * t432 + t738;
t169 = -t430 * t727 + t736;
t168 = t430 * t322 + t432 * t488;
t167 = t430 * t313 - t432 * t489;
t166 = t430 * t312 - t735;
t165 = -t313 * t432 - t737;
t164 = -t312 * t432 - t430 * t490;
t161 = Icges(5,1) * t275 + Icges(5,4) * t274 - Icges(5,5) * t465;
t160 = Icges(5,1) * t273 + Icges(5,4) * t272 + Icges(5,5) * t466;
t159 = Icges(5,4) * t275 + Icges(5,2) * t274 - Icges(5,6) * t465;
t158 = Icges(5,4) * t273 + Icges(5,2) * t272 + Icges(5,6) * t466;
t153 = t327 * t432 + t430 * t328 + t637;
t148 = t389 + (-t615 + (t431 * t690 - t694) * qJD(2)) * t430 + ((-rSges(4,1) - pkin(7)) * t430 + t455 * t432) * qJD(1);
t147 = -pkin(2) * t586 + (t550 - t698) * t621 + t560 + t592;
t101 = t232 * t341 + t234 * t342 + t600;
t100 = t231 * t341 + t233 * t342 + t602;
t99 = t339 * t232 + t340 * t234 + t599;
t98 = t339 * t231 + t340 * t233 + t601;
t74 = (qJD(1) * t328 + t560) * t432 + (t534 * t618 + (-t327 - t347 + t733) * qJD(1)) * t430 + t595;
t31 = t121 * t651 + t308 * t123 + t309 * t125 + t198 * t215 + t199 * t217 - t213 * t465;
t30 = t122 * t651 + t308 * t124 + t309 * t126 + t198 * t214 + t199 * t216 - t212 * t465;
t29 = t121 * t652 + t310 * t123 + t311 * t125 + t196 * t215 + t197 * t217 + t213 * t466;
t28 = t122 * t652 + t310 * t124 + t311 * t126 + t196 * t214 + t197 * t216 + t212 * t466;
t16 = qJD(1) * t530 + t30 * t430 - t31 * t432;
t15 = qJD(1) * t529 + t28 * t430 - t29 * t432;
t8 = (-qJD(2) * t530 + t58) * t429 + (-qJD(1) * t60 + qJD(2) * t134 + t30 * t432 + t31 * t430) * t431;
t7 = (-qJD(2) * t529 + t57) * t429 + (-qJD(1) * t61 + qJD(2) * t135 + t28 * t432 + t29 * t430) * t431;
t1 = [t525 + t526 - t402 * t282 * t662 - t288 * t661 - t291 * t583 + (t132 * t68 + t133 * t67) * t720 + (t143 * t77 + t144 * t76) * t721 + (t102 * t163 + t103 * t162) * t722 + (t147 * t240 + t148 * t239) * t723 + (t210 * t279 + t211 * t278) * t724 + (Icges(5,3) * t429 + t506 + t518 + t752) * t617 + (-t190 * t402 - t221 * t409 - t287 * t427 - t511 * t617 - t668 - t670) * t431 + (Icges(5,3) * t431 + t305 * t427 + t306 * t426 + t429 * t511 + t508 + t523 - t753) * t619; m(4) * (t147 * t283 + t148 * t284 + t182 * t239 + t183 * t240) + m(5) * (t102 * t227 + t103 * t228 + t137 * t162 + t138 * t163) + m(6) * (t143 * t94 + t144 * t95 + t174 * t76 + t175 * t77) + m(7) * (t132 * t71 + t133 * t72 + t139 * t67 + t140 * t68) + (t305 * t748 + t306 * t747 - t341 * t287 / 0.2e1 - t342 * t288 / 0.2e1 + m(3) * (-t211 * t378 - t278 * t358) - t27 / 0.2e1 - t49 / 0.2e1 + t579 * t432 + (Icges(5,5) * t747 + Icges(5,6) * t748 + t317 * t745 + t320 * t696 + t466 * t749 + t703 * t750) * t429 - t603) * t432 + (m(3) * (-t210 * t378 - t279 * t358) + t50 / 0.2e1 + t26 / 0.2e1 + t274 * t712 + t275 * t711 + t287 * t710 + t288 * t709 + t579 * t430 + (Icges(5,5) * t275 / 0.2e1 + Icges(5,6) * t274 / 0.2e1 + t465 * t749 + t750 * t746 + (t316 + t728) * t745) * t429 + t604) * t430 + ((t158 * t705 + t160 * t707 + t315 * t745 + t318 * t696 + t703 * t751) * t432 + (t159 * t706 + t161 * t708 + t751 * t746 + (t314 + t729) * t745) * t430) * t431 + ((t430 * t549 - t432 * t548) * t431 + (t430 * t451 + t432 * t450) * t429) * qJD(2) + ((-t279 * t699 + t129 / 0.2e1 + t305 * t710 + t306 * t709 + t92 / 0.2e1 + t549 * t429 - t451 * t431 - t597) * t432 + (t278 * t699 + t93 / 0.2e1 + t130 / 0.2e1 + t341 * t712 + t342 * t711 + t548 * t429 + t450 * t431 - t596) * t430) * qJD(1); -t432 * t15 - t432 * t12 + t430 * t16 + t430 * t13 - t432 * ((t251 * t432 + (t170 - t736) * qJD(1)) * t432 + (t171 * qJD(1) + (t318 * t617 + t320 * t619 + t623) * t430 + (-t252 + (t429 * t728 + t431 * t729) * qJD(2) + t488 * qJD(1)) * t432) * t430) - t432 * ((t242 * t432 + (t165 + t735) * qJD(1)) * t432 + (t164 * qJD(1) + (-t315 * t617 - t317 * t619 + t624) * t430 + (-t241 + (t314 * t431 + t316 * t429) * qJD(2) - t489 * qJD(1)) * t432) * t430) + t430 * ((t430 * t241 + (t166 + t737) * qJD(1)) * t430 + (t167 * qJD(1) + (t314 * t617 + t316 * t619) * t432 + (-t242 + (-t315 * t431 - t317 * t429) * qJD(2) + (t313 - t490) * qJD(1)) * t430) * t432) + t430 * ((t430 * t252 + (t169 - t738) * qJD(1)) * t430 + (t168 * qJD(1) + (t617 * t729 + t619 * t728) * t432 + (-t251 + (t318 * t431 + t320 * t429) * qJD(2) + (t322 + t487) * qJD(1)) * t430) * t432) - t432 * ((-t341 * t158 - t342 * t160 - t272 * t232 - t273 * t234 + (t100 - t599) * qJD(1)) * t432 + (t341 * t159 + t342 * t161 + t272 * t231 + t273 * t233 + (t101 + t601) * qJD(1)) * t430) + t430 * ((t339 * t159 + t340 * t161 + t274 * t231 + t275 * t233 + (t99 - t602) * qJD(1)) * t430 + (-t339 * t158 - t340 * t160 - t274 * t232 - t275 * t234 + (t98 + t600) * qJD(1)) * t432) + (t139 * t72 + t140 * t71 + t59 * t17) * t720 + (t174 * t95 + t175 * t94 + t23 * t75) * t721 + (t105 * t48 + t137 * t228 + t138 * t227) * t722 + (t153 * t74 + t182 * t284 + t183 * t283) * t723 + ((t430 * t325 + t326 * t432) * ((qJD(1) * t325 - t432 * t474 + t628) * t432 + (-t430 * t474 + (-t326 + t732) * qJD(1)) * t430) + t625 * t378 * t358) * t724 + (t56 + t61 + (-t101 - t164 - t171) * t432 + (t100 + t165 + t170) * t430) * t621 + (t55 + t60 + (-t166 - t169 - t99) * t432 + (t167 + t168 + t98) * t430) * t620; 0.2e1 * (t504 * t713 + t503 * t715 + t501 * t717 + (t239 * t432 + t240 * t430) * t719) * t617 + 0.2e1 * ((t147 * t430 + t148 * t432 - t239 * t621 + t240 * t620) * t719 + t434) * t429; 0.2e1 * (t459 * t713 + t457 * t715 + t456 * t717 + (t283 * t618 + t284 * t616 - t74) * t719) * t431 + 0.2e1 * ((qJD(2) * t153 + t182 * t432 + t183 * t430 + t283 * t620 - t284 * t621) * t719 + t433) * t429; 0.4e1 * (t719 + t568) * t277; 0.2e1 * (t501 * t718 + t503 * t716 + t504 * t714) * t619 + 0.2e1 * t434 * t431; 0.2e1 * (t456 * t718 + t457 * t716 + t459 * t714) * t429 + 0.2e1 * t433 * t431; 0.2e1 * t568 * (-t429 ^ 2 + t431 ^ 2) * t740 * qJD(2); -0.4e1 * t568 * t277; t438 + (t604 * t432 + t603 * t430 + (t430 * t597 - t432 * t596) * qJD(1)) * t431 + m(6) * (t143 * t70 + t144 * t69 + t151 * t77 + t152 * t76) + m(7) * (t132 * t33 + t133 * t32 + t67 * t83 + t68 * t82) + (t430 * t596 + t432 * t597) * t619 + t687; t441 + ((qJD(1) * t96 - t35) * t704 - t7 / 0.2e1 + t60 * t577 + t42 * t696) * t432 + ((qJD(1) * t97 + t34) * t704 + t8 / 0.2e1 + t61 * t577 + t43 * t696) * t430 + (t16 * t702 + (t96 * t430 - t432 * t97) * t739 + t15 * t703 + (t61 * t702 - t430 * t60 / 0.2e1) * qJD(1)) * t431 + m(6) * (-t136 * t23 + t151 * t94 + t152 * t95 + t174 * t69 + t175 * t70 + t44 * t75) + m(7) * (t139 * t32 + t140 * t33 + t73 * t17 + t18 * t59 + t71 * t82 + t72 * t83); 0.2e1 * (t458 * t715 + t460 * t713) * t431 + 0.2e1 * t435 * t429; 0.2e1 * (t458 * t716 + t460 * t714) * t429 + 0.2e1 * t435 * t431; (t73 * t18 + t32 * t83 + t33 * t82) * t720 + (-t136 * t44 + t151 * t70 + t152 * t69) * t721 + ((t557 * t432 + (-t40 - t574) * t430) * qJD(2) + t687) * t429 + (t432 * t8 + t430 * t7 + t429 * (t34 * t432 + t35 * t430) + (t146 * t429 + (t430 * t97 + t432 * t96) * t431) * qJD(2) + (t430 * t557 + t432 * t574) * qJD(1)) * t431 + t567; t438 + m(7) * (t132 * t65 + t133 * t64 + t149 * t68 + t150 * t67); m(7) * (t131 * t17 + t139 * t64 + t140 * t65 + t149 * t71 + t150 * t72 + t41 * t59) + t441; m(7) * (t437 * t429 + t431 * t730); m(7) * (-t429 * t730 + t437 * t431); m(7) * (t131 * t18 + t149 * t33 + t150 * t32 + t41 * t73 + t64 * t83 + t65 * t82) + t436; (t131 * t41 + t149 * t65 + t150 * t64) * t720 + t436;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
