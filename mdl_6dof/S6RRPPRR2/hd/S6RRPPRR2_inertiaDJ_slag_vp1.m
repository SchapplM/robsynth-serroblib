% Calculate time derivative of joint inertia matrix for
% S6RRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 08:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR2_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR2_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR2_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR2_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:49:52
% EndTime: 2019-03-09 08:50:35
% DurationCPUTime: 26.36s
% Computational Cost: add. (53639->1159), mult. (46869->1595), div. (0->0), fcn. (44354->12), ass. (0->553)
t422 = sin(qJ(1));
t687 = t422 / 0.2e1;
t424 = cos(qJ(1));
t686 = -t424 / 0.2e1;
t715 = -qJD(1) / 0.2e1;
t414 = qJ(2) + pkin(10);
t403 = cos(t414);
t598 = qJD(2) * t403;
t551 = t598 / 0.2e1;
t401 = sin(t414);
t601 = qJD(1) * t422;
t567 = t401 * t601;
t714 = t424 * t551 - t567 / 0.2e1;
t600 = qJD(1) * t424;
t552 = t600 / 0.2e1;
t713 = t401 * t552 + t422 * t551;
t418 = cos(pkin(11));
t396 = t418 * pkin(4) + pkin(3);
t675 = pkin(3) - t396;
t553 = t675 * t403;
t651 = qJ(4) * t401;
t712 = t553 + t651;
t596 = qJD(2) * t422;
t559 = t403 * t596;
t438 = t401 * t600 + t559;
t711 = qJD(2) / 0.2e1;
t710 = t401 * t675;
t423 = cos(qJ(2));
t397 = pkin(2) * t423 + pkin(1);
t676 = pkin(1) - t397;
t709 = t422 * t676;
t421 = sin(qJ(2));
t664 = Icges(3,4) * t423;
t502 = -Icges(3,2) * t421 + t664;
t331 = Icges(3,6) * t422 + t424 * t502;
t665 = Icges(3,4) * t421;
t509 = Icges(3,1) * t423 - t665;
t333 = Icges(3,5) * t422 + t424 * t509;
t473 = t331 * t421 - t333 * t423;
t448 = t473 * t422;
t330 = -Icges(3,6) * t424 + t422 * t502;
t332 = -Icges(3,5) * t424 + t422 * t509;
t474 = t330 * t421 - t332 * t423;
t449 = t474 * t424;
t662 = Icges(4,4) * t403;
t500 = -Icges(4,2) * t401 + t662;
t310 = Icges(4,6) * t422 + t424 * t500;
t663 = Icges(4,4) * t401;
t507 = Icges(4,1) * t403 - t663;
t312 = Icges(4,5) * t422 + t424 * t507;
t475 = t310 * t401 - t312 * t403;
t450 = t475 * t422;
t309 = -Icges(4,6) * t424 + t422 * t500;
t311 = -Icges(4,5) * t424 + t422 * t507;
t476 = t309 * t401 - t311 * t403;
t451 = t476 * t424;
t408 = t422 * rSges(4,3);
t640 = t401 * t424;
t708 = -rSges(4,2) * t640 + t408;
t494 = Icges(4,5) * t403 - Icges(4,6) * t401;
t307 = -Icges(4,3) * t424 + t422 * t494;
t495 = Icges(3,5) * t423 - Icges(3,6) * t421;
t328 = -Icges(3,3) * t424 + t422 * t495;
t603 = qJD(1) * t403;
t539 = -qJD(5) + t603;
t594 = qJD(2) * t424;
t560 = t401 * t594;
t707 = t422 * t539 + t560;
t413 = qJD(5) + qJD(6);
t542 = -t413 + t603;
t706 = t422 * t542 + t560;
t561 = t401 * t596;
t705 = t424 * t539 - t561;
t704 = t424 * t542 - t561;
t420 = -pkin(8) - qJ(4);
t411 = -pkin(9) + t420;
t607 = t411 - t420;
t546 = t607 * t401;
t417 = sin(pkin(11));
t628 = t422 * t417;
t391 = pkin(4) * t628;
t636 = t403 * t424;
t614 = -t396 * t636 - t391;
t412 = pkin(11) + qJ(5);
t402 = cos(t412);
t357 = pkin(5) * t402 + t396;
t400 = sin(t412);
t679 = pkin(5) * t400;
t680 = pkin(4) * t417;
t365 = t679 + t680;
t616 = t357 * t636 + t422 * t365;
t170 = -t424 * t546 + t614 + t616;
t404 = qJ(6) + t412;
t394 = sin(t404);
t395 = cos(t404);
t631 = t422 * t395;
t303 = -t394 * t636 + t631;
t632 = t422 * t394;
t304 = t395 * t636 + t632;
t207 = t304 * rSges(7,1) + t303 * rSges(7,2) + rSges(7,3) * t640;
t622 = t170 + t207;
t613 = t357 - t396;
t547 = t613 * t403;
t446 = -t401 * t411 + t547;
t635 = t417 * t424;
t392 = pkin(4) * t635;
t641 = t401 * t422;
t611 = t420 * t641 + t392;
t169 = -t365 * t424 + t422 * t446 + t611;
t301 = -t395 * t424 - t403 * t632;
t302 = -t394 * t424 + t403 * t631;
t519 = -t302 * rSges(7,1) - t301 * rSges(7,2);
t206 = rSges(7,3) * t641 - t519;
t623 = t169 + t206;
t703 = -t422 * t623 - t424 * t622;
t702 = 2 * m(3);
t701 = 2 * m(4);
t700 = 2 * m(5);
t699 = 2 * m(6);
t698 = 2 * m(7);
t697 = m(5) / 0.2e1;
t696 = m(6) / 0.2e1;
t695 = m(7) / 0.2e1;
t498 = Icges(5,4) * t418 - Icges(5,2) * t417;
t294 = -Icges(5,6) * t403 + t401 * t498;
t694 = t294 / 0.2e1;
t505 = Icges(5,1) * t418 - Icges(5,4) * t417;
t295 = -Icges(5,5) * t403 + t401 * t505;
t693 = t295 / 0.2e1;
t627 = t422 * t418;
t341 = -t403 * t635 + t627;
t692 = t341 / 0.2e1;
t634 = t418 * t424;
t342 = t403 * t634 + t628;
t691 = t342 / 0.2e1;
t690 = -t403 / 0.2e1;
t689 = -t417 / 0.2e1;
t688 = t418 / 0.2e1;
t685 = t424 / 0.2e1;
t378 = rSges(3,1) * t421 + rSges(3,2) * t423;
t684 = m(3) * t378;
t683 = pkin(2) * t421;
t682 = pkin(3) * t401;
t681 = pkin(3) * t403;
t678 = t422 * pkin(7);
t410 = t424 * pkin(7);
t677 = qJD(1) / 0.2e1;
t419 = -qJ(3) - pkin(7);
t674 = -pkin(7) - t419;
t673 = rSges(3,1) * t423;
t672 = rSges(4,1) * t403;
t671 = rSges(3,2) * t421;
t670 = rSges(3,3) * t424;
t669 = rSges(6,3) * t401;
t668 = rSges(7,3) * t401;
t409 = t422 * rSges(3,3);
t667 = -rSges(5,3) - qJ(4);
t666 = -rSges(7,3) + t411;
t661 = Icges(6,4) * t400;
t660 = Icges(6,4) * t402;
t659 = Icges(7,4) * t394;
t658 = Icges(7,4) * t395;
t630 = t422 * t400;
t321 = -t402 * t424 - t403 * t630;
t629 = t422 * t402;
t322 = -t400 * t424 + t403 * t629;
t522 = -t322 * rSges(6,1) - t321 * rSges(6,2);
t218 = rSges(6,3) * t641 - t522;
t650 = t218 * t424;
t503 = Icges(7,1) * t395 - t659;
t278 = -Icges(7,5) * t403 + t401 * t503;
t649 = t278 * t395;
t504 = Icges(6,1) * t402 - t661;
t285 = -Icges(6,5) * t403 + t401 * t504;
t648 = t285 * t402;
t643 = t357 * t401;
t642 = t401 * t413;
t639 = t403 * t417;
t638 = t403 * t418;
t637 = t403 * t420;
t633 = t419 * t424;
t626 = qJ(4) + t420;
t543 = -t403 * t413 + qJD(1);
t471 = t422 * t543;
t182 = -t394 * t704 + t395 * t471;
t183 = t394 * t471 + t395 * t704;
t520 = t183 * rSges(7,1) + t182 * rSges(7,2);
t117 = rSges(7,3) * t438 + t520;
t558 = t403 * t594;
t625 = t117 * t640 + t206 * t558;
t470 = t424 * t543;
t180 = t394 * t706 + t395 * t470;
t181 = t394 * t470 - t395 * t706;
t576 = t181 * rSges(7,1) + t180 * rSges(7,2) + rSges(7,3) * t558;
t116 = -rSges(7,3) * t567 + t576;
t591 = qJD(5) * t402;
t585 = pkin(5) * t591;
t570 = t365 * t600 + t411 * t567 + t422 * t585;
t586 = qJD(5) * t679;
t612 = qJD(1) * t392 + t420 * t567;
t624 = -t613 * t560 + (-t613 * t601 + (-qJD(2) * t607 - t586) * t424) * t403 + t570 - t612 + t116;
t518 = rSges(7,1) * t395 - rSges(7,2) * t394;
t190 = (-rSges(7,1) * t394 - rSges(7,2) * t395) * t642 + (t403 * t518 + t668) * qJD(2);
t592 = qJD(5) * t401;
t557 = t400 * t592;
t193 = -pkin(5) * t557 + (t547 - t546) * qJD(2);
t621 = -t190 - t193;
t279 = -rSges(7,3) * t403 + t401 * t518;
t599 = qJD(2) * t401;
t620 = t207 * t599 + t279 * t567;
t247 = t401 * t613 + t403 * t607;
t619 = t247 + t279;
t148 = t403 * t206 + t279 * t641;
t305 = t410 + t633 - t709;
t385 = t424 * t397;
t306 = -pkin(1) * t424 + t422 * t674 + t385;
t618 = t422 * t305 + t424 * t306;
t384 = pkin(3) * t636;
t335 = qJ(4) * t640 + t384;
t617 = -t306 - t335;
t354 = -qJ(4) * t403 + t682;
t565 = t421 * t601;
t389 = pkin(2) * t565;
t615 = t354 * t601 + t389;
t610 = rSges(4,2) * t567 + rSges(4,3) * t600;
t593 = qJD(4) * t401;
t377 = t424 * t593;
t405 = qJD(3) * t422;
t609 = t377 + t405;
t608 = t424 * t673 + t409;
t606 = t422 ^ 2 + t424 ^ 2;
t308 = Icges(4,3) * t422 + t424 * t494;
t605 = qJD(1) * t308;
t329 = Icges(3,3) * t422 + t424 * t495;
t604 = qJD(1) * t329;
t602 = qJD(1) * t419;
t597 = qJD(2) * t421;
t595 = qJD(2) * t423;
t590 = t424 * t671;
t588 = pkin(2) * t597;
t587 = pkin(2) * t595;
t540 = -qJD(5) * t403 + qJD(1);
t468 = t422 * t540;
t198 = -t400 * t705 + t402 * t468;
t199 = t400 * t468 + t402 * t705;
t125 = Icges(6,5) * t199 + Icges(6,6) * t198 + Icges(6,3) * t438;
t127 = Icges(6,4) * t199 + Icges(6,2) * t198 + Icges(6,6) * t438;
t129 = Icges(6,1) * t199 + Icges(6,4) * t198 + Icges(6,5) * t438;
t212 = Icges(6,5) * t322 + Icges(6,6) * t321 + Icges(6,3) * t641;
t214 = Icges(6,4) * t322 + Icges(6,2) * t321 + Icges(6,6) * t641;
t216 = Icges(6,1) * t322 + Icges(6,4) * t321 + Icges(6,5) * t641;
t484 = -t214 * t400 + t216 * t402;
t34 = (qJD(2) * t484 - t125) * t403 + (qJD(2) * t212 - t127 * t400 + t129 * t402 + (-t214 * t402 - t216 * t400) * qJD(5)) * t401;
t492 = Icges(6,5) * t402 - Icges(6,6) * t400;
t208 = (-Icges(6,5) * t400 - Icges(6,6) * t402) * t592 + (Icges(6,3) * t401 + t403 * t492) * qJD(2);
t497 = -Icges(6,2) * t400 + t660;
t209 = (-Icges(6,2) * t402 - t661) * t592 + (Icges(6,6) * t401 + t403 * t497) * qJD(2);
t210 = (-Icges(6,1) * t400 - t660) * t592 + (Icges(6,5) * t401 + t403 * t504) * qJD(2);
t283 = -Icges(6,3) * t403 + t401 * t492;
t284 = -Icges(6,6) * t403 + t401 * t497;
t59 = t198 * t284 + t199 * t285 + t208 * t641 + t321 * t209 + t322 * t210 + t283 * t438;
t584 = t34 / 0.2e1 + t59 / 0.2e1;
t467 = t424 * t540;
t196 = t400 * t707 + t402 * t467;
t197 = t400 * t467 - t402 * t707;
t437 = t558 - t567;
t124 = Icges(6,5) * t197 + Icges(6,6) * t196 + Icges(6,3) * t437;
t126 = Icges(6,4) * t197 + Icges(6,2) * t196 + Icges(6,6) * t437;
t128 = Icges(6,1) * t197 + Icges(6,4) * t196 + Icges(6,5) * t437;
t323 = -t400 * t636 + t629;
t324 = t402 * t636 + t630;
t213 = Icges(6,5) * t324 + Icges(6,6) * t323 + Icges(6,3) * t640;
t215 = Icges(6,4) * t324 + Icges(6,2) * t323 + Icges(6,6) * t640;
t217 = Icges(6,1) * t324 + Icges(6,4) * t323 + Icges(6,5) * t640;
t483 = -t215 * t400 + t217 * t402;
t35 = (qJD(2) * t483 - t124) * t403 + (qJD(2) * t213 - t126 * t400 + t128 * t402 + (-t215 * t402 - t217 * t400) * qJD(5)) * t401;
t58 = t196 * t284 + t197 * t285 + t208 * t640 + t323 * t209 + t324 * t210 + t283 * t437;
t583 = t35 / 0.2e1 + t58 / 0.2e1;
t340 = t403 * t627 - t635;
t458 = t403 * t628 + t634;
t233 = Icges(5,5) * t340 - Icges(5,6) * t458 + Icges(5,3) * t641;
t582 = t233 * t641;
t581 = t233 * t640;
t234 = Icges(5,5) * t342 + Icges(5,6) * t341 + Icges(5,3) * t640;
t580 = t234 * t641;
t579 = t234 * t640;
t134 = t283 * t641 + t284 * t321 + t285 * t322;
t98 = -t212 * t403 + t401 * t484;
t578 = t98 / 0.2e1 + t134 / 0.2e1;
t135 = t283 * t640 + t323 * t284 + t324 * t285;
t99 = -t213 * t403 + t401 * t483;
t577 = t99 / 0.2e1 + t135 / 0.2e1;
t575 = t197 * rSges(6,1) + t196 * rSges(6,2) + rSges(6,3) * t558;
t463 = -t420 * t640 - t614;
t230 = t463 - t335;
t574 = -t230 + t617;
t568 = qJD(3) * t424 + t419 * t601 + t422 * t588;
t573 = t422 * ((-t424 * t676 - t678) * qJD(1) - t568) + t424 * (-t424 * t588 + t405 + (t424 * t674 + t709) * qJD(1)) + t305 * t600;
t271 = t403 * t626 - t710;
t572 = t271 * t601 + t615;
t272 = qJD(1) * t458 + t417 * t560;
t273 = -qJD(1) * t340 - t418 * t560;
t571 = t273 * rSges(5,1) + t272 * rSges(5,2) + rSges(5,3) * t558;
t219 = t324 * rSges(6,1) + t323 * rSges(6,2) + rSges(6,3) * t640;
t246 = t342 * rSges(5,1) + t341 * rSges(5,2) + rSges(5,3) * t640;
t569 = t396 * t561 + t420 * t438;
t521 = rSges(6,1) * t402 - rSges(6,2) * t400;
t287 = -rSges(6,3) * t403 + t401 * t521;
t564 = t287 * t601;
t496 = -Icges(7,2) * t394 + t658;
t277 = -Icges(7,6) * t403 + t401 * t496;
t563 = t277 * t598;
t562 = t284 * t598;
t556 = t641 / 0.2e1;
t555 = t640 / 0.2e1;
t554 = (t494 + t495) * t711;
t550 = -t354 - t683;
t355 = rSges(4,1) * t401 + rSges(4,2) * t403;
t549 = -t355 - t683;
t548 = t619 * t424;
t545 = -t357 * t403 - t397;
t544 = -t422 * t419 + t385;
t541 = t403 * t586;
t538 = t403 * t117 + t190 * t641 + t279 * t438;
t517 = t651 + t681;
t334 = t517 * t422;
t537 = t422 * t334 + t424 * t335 + t618;
t532 = -t271 + t550;
t524 = rSges(5,1) * t418 - rSges(5,2) * t417;
t296 = -rSges(5,3) * t403 + t401 * t524;
t531 = -t296 + t550;
t530 = -qJD(2) * t517 + qJD(4) * t403 - t587;
t491 = Icges(7,5) * t395 - Icges(7,6) * t394;
t276 = -Icges(7,3) * t403 + t401 * t491;
t139 = -t276 * t403 + (-t277 * t394 + t649) * t401;
t122 = t276 * t641 + t277 * t301 + t278 * t302;
t200 = Icges(7,5) * t302 + Icges(7,6) * t301 + Icges(7,3) * t641;
t202 = Icges(7,4) * t302 + Icges(7,2) * t301 + Icges(7,6) * t641;
t204 = Icges(7,1) * t302 + Icges(7,4) * t301 + Icges(7,5) * t641;
t81 = t200 * t641 + t202 * t301 + t204 * t302;
t201 = Icges(7,5) * t304 + Icges(7,6) * t303 + Icges(7,3) * t640;
t203 = Icges(7,4) * t304 + Icges(7,2) * t303 + Icges(7,6) * t640;
t205 = Icges(7,1) * t304 + Icges(7,4) * t303 + Icges(7,5) * t640;
t82 = t201 * t641 + t203 * t301 + t205 * t302;
t516 = t422 * t81 + t424 * t82;
t41 = -t122 * t403 + t401 * t516;
t123 = t276 * t640 + t303 * t277 + t304 * t278;
t83 = t200 * t640 + t303 * t202 + t304 * t204;
t84 = t201 * t640 + t303 * t203 + t304 * t205;
t515 = t422 * t83 + t424 * t84;
t42 = -t123 * t403 + t401 * t515;
t111 = Icges(7,5) * t183 + Icges(7,6) * t182 + Icges(7,3) * t438;
t113 = Icges(7,4) * t183 + Icges(7,2) * t182 + Icges(7,6) * t438;
t115 = Icges(7,1) * t183 + Icges(7,4) * t182 + Icges(7,5) * t438;
t20 = t111 * t640 + t303 * t113 + t304 * t115 + t180 * t202 + t181 * t204 + t200 * t437;
t110 = Icges(7,5) * t181 + Icges(7,6) * t180 + Icges(7,3) * t437;
t112 = Icges(7,4) * t181 + Icges(7,2) * t180 + Icges(7,6) * t437;
t114 = Icges(7,1) * t181 + Icges(7,4) * t180 + Icges(7,5) * t437;
t21 = t110 * t640 + t303 * t112 + t304 * t114 + t180 * t203 + t181 * t205 + t201 * t437;
t187 = (-Icges(7,5) * t394 - Icges(7,6) * t395) * t642 + (Icges(7,3) * t401 + t403 * t491) * qJD(2);
t188 = (-Icges(7,2) * t395 - t659) * t642 + (Icges(7,6) * t401 + t403 * t496) * qJD(2);
t189 = (-Icges(7,1) * t394 - t658) * t642 + (Icges(7,5) * t401 + t403 * t503) * qJD(2);
t51 = t180 * t277 + t181 * t278 + t187 * t640 + t303 * t188 + t304 * t189 + t276 * t437;
t61 = t84 * t422 - t424 * t83;
t5 = (qJD(2) * t515 - t51) * t403 + (-qJD(1) * t61 + qJD(2) * t123 + t20 * t422 + t21 * t424) * t401;
t486 = -t202 * t394 + t204 * t395;
t95 = -t200 * t403 + t401 * t486;
t485 = -t203 * t394 + t205 * t395;
t96 = -t201 * t403 + t401 * t485;
t511 = t422 * t95 + t424 * t96;
t22 = t111 * t641 + t301 * t113 + t302 * t115 + t182 * t202 + t183 * t204 + t200 * t438;
t23 = t110 * t641 + t301 * t112 + t302 * t114 + t182 * t203 + t183 * t205 + t201 * t438;
t52 = t182 * t277 + t183 * t278 + t187 * t641 + t301 * t188 + t302 * t189 + t276 * t438;
t60 = t82 * t422 - t424 * t81;
t6 = (qJD(2) * t516 - t52) * t403 + (-qJD(1) * t60 + qJD(2) * t122 + t22 * t422 + t23 * t424) * t401;
t529 = t42 * t558 + t5 * t640 + t6 * t641 + (-t139 * t403 + t401 * t511) * t599 + t438 * t41;
t528 = -t671 + t673;
t527 = -rSges(4,2) * t401 + t672;
t274 = qJD(1) * t341 + t417 * t561;
t275 = qJD(1) * t342 - t418 * t561;
t526 = -t275 * rSges(5,1) - t274 * rSges(5,2);
t525 = -t340 * rSges(5,1) + rSges(5,2) * t458;
t523 = t199 * rSges(6,1) + t198 * rSges(6,2);
t91 = t212 * t641 + t214 * t321 + t216 * t322;
t92 = t213 * t641 + t215 * t321 + t217 * t322;
t63 = t92 * t422 - t424 * t91;
t514 = t422 * t91 + t424 * t92;
t93 = t212 * t640 + t323 * t214 + t324 * t216;
t94 = t213 * t640 + t323 * t215 + t324 * t217;
t64 = t94 * t422 - t424 * t93;
t513 = t422 * t93 + t424 * t94;
t512 = t96 * t422 - t424 * t95;
t510 = t422 * t98 + t424 * t99;
t506 = Icges(4,1) * t401 + t662;
t499 = Icges(4,2) * t403 + t663;
t493 = Icges(5,5) * t418 - Icges(5,6) * t417;
t435 = t401 * t666 + t545;
t143 = (t365 - t419) * t424 + t435 * t422 + t519;
t144 = -t411 * t640 + t207 + t544 + t616;
t490 = t143 * t424 + t144 * t422;
t149 = -t403 * t207 - t279 * t640;
t489 = t148 * t424 + t149 * t422;
t459 = -t396 * t403 - t397 - t669;
t428 = t422 * t459 - t633;
t150 = t428 + t522 + t611;
t151 = t463 + t544 + t219;
t488 = t150 * t424 + t151 * t422;
t439 = t401 * t667 - t397 - t681;
t425 = t422 * t439 - t633;
t166 = t425 + t525;
t167 = t544 + t246 + t335;
t487 = t166 * t424 + t167 * t422;
t482 = -t219 * t422 + t650;
t481 = -t422 * t218 - t219 * t424;
t472 = -t403 * t411 - t643;
t469 = -t287 + t532;
t318 = rSges(4,1) * t636 + t708;
t466 = -(-t401 * t626 - t553) * qJD(2) + t530;
t465 = -(rSges(5,3) * t401 + t403 * t524) * qJD(2) + t530;
t464 = -pkin(1) - t528;
t228 = t531 * t424;
t462 = -t397 - t527;
t367 = qJ(4) * t558;
t372 = pkin(3) * t561;
t461 = t422 * (qJ(4) * t438 + qJD(1) * t384 + t422 * t593 - t372) + t424 * (-qJ(4) * t567 + t367 + t377 + (-t403 * t601 - t560) * pkin(3)) + t334 * t600 + t573;
t229 = -t422 * t712 - t611;
t460 = t422 * t229 + t424 * t230 + t537;
t456 = t532 - t619;
t211 = (-rSges(6,1) * t400 - rSges(6,2) * t402) * t592 + (t403 * t521 + t669) * qJD(2);
t454 = -t211 + t466;
t453 = qJD(2) * t378;
t452 = qJD(2) * t355;
t447 = t169 * t424 - t422 * t622;
t444 = qJD(2) * t506;
t442 = qJD(2) * t499;
t441 = qJD(2) * (-Icges(3,5) * t421 - Icges(3,6) * t423);
t440 = qJD(2) * (-Icges(4,5) * t401 - Icges(4,6) * t403);
t155 = t469 * t424;
t436 = t466 + t621;
t133 = t456 * t424;
t433 = t422 * (-qJ(4) * t559 + t372 + (-t424 * t712 + t391) * qJD(1) - t569) + t424 * (-t367 + (-t637 + t710) * t594 + t712 * t601 + t612) + t229 * t600 + t461;
t12 = qJD(1) * t515 - t20 * t424 + t21 * t422;
t13 = qJD(1) * t516 - t22 * t424 + t23 * t422;
t26 = (qJD(2) * t486 - t111) * t403 + (qJD(2) * t200 + (-t202 * t413 + t115) * t395 + (-t204 * t413 - t113) * t394) * t401;
t27 = (qJD(2) * t485 - t110) * t403 + (qJD(2) * t201 + (-t203 * t413 + t114) * t395 + (-t205 * t413 - t112) * t394) * t401;
t432 = t12 * t555 + t13 * t556 + t5 * t687 + t6 * t686 + (qJD(1) * t511 - t26 * t424 + t27 * t422) * t690 + t41 * t601 / 0.2e1 + t42 * t552 + t512 * t599 / 0.2e1 + t714 * t61 + t713 * t60;
t431 = rSges(3,2) * t565 + rSges(3,3) * t600 - t424 * t453;
t430 = -t403 * t187 + t276 * t599 + t598 * t649 + (t189 * t401 - t277 * t642) * t395;
t429 = -t403 * t208 + t283 * t599 + t598 * t648 + (t210 * t402 - t284 * t591) * t401;
t137 = t139 * t599;
t62 = (-t563 + (-t278 * t413 - t188) * t401) * t394 + t430;
t9 = t137 + (qJD(2) * t511 - t62) * t403 + (-qJD(1) * t512 + t26 * t422 + t27 * t424) * t401;
t427 = -t403 * t9 - t42 * t567 + t529;
t426 = t137 + (t26 + t52) * t556 + (t27 + t51) * t555 + (t123 + t96) * t714 + (t122 + t95) * t713;
t364 = t528 * qJD(2);
t346 = t527 * qJD(2);
t337 = -t590 + t608;
t336 = t422 * t528 - t670;
t317 = -rSges(4,3) * t424 + t422 * t527;
t300 = t549 * t424;
t299 = t549 * t422;
t292 = t678 + (pkin(1) - t671) * t424 + t608;
t291 = t422 * t464 + t410 + t670;
t282 = (Icges(5,5) * t401 + t403 * t505) * qJD(2);
t281 = (Icges(5,6) * t401 + t403 * t498) * qJD(2);
t266 = t318 + t544;
t265 = (rSges(4,3) - t419) * t424 + t462 * t422;
t258 = t422 * t441 + t604;
t257 = -qJD(1) * t328 + t424 * t441;
t245 = rSges(5,3) * t641 - t525;
t240 = t422 * t440 + t605;
t239 = -qJD(1) * t307 + t424 * t440;
t238 = Icges(5,1) * t342 + Icges(5,4) * t341 + Icges(5,5) * t640;
t237 = Icges(5,1) * t340 - Icges(5,4) * t458 + Icges(5,5) * t641;
t236 = Icges(5,4) * t342 + Icges(5,2) * t341 + Icges(5,6) * t640;
t235 = Icges(5,4) * t340 - Icges(5,2) * t458 + Icges(5,6) * t641;
t232 = t378 * t596 + ((-rSges(3,3) - pkin(7)) * t422 + t464 * t424) * qJD(1);
t231 = (t410 + (-pkin(1) - t673) * t422) * qJD(1) + t431;
t227 = t531 * t422;
t224 = -t355 * t600 - t422 * t346 + (-t421 * t600 - t422 * t595) * pkin(2);
t223 = t355 * t601 + t389 + (-t346 - t587) * t424;
t186 = t206 * t640;
t179 = t422 * t329 - t424 * t473;
t178 = t422 * t328 - t449;
t177 = -t329 * t424 - t448;
t176 = -t328 * t424 - t422 * t474;
t172 = t355 * t596 + (t424 * t462 - t408) * qJD(1) + t568;
t171 = t405 + (-t397 - t672) * t601 + (qJD(2) * t549 - t602) * t424 + t610;
t165 = t422 * t308 - t424 * t475;
t164 = t422 * t307 - t451;
t163 = -t308 * t424 - t450;
t162 = -t307 * t424 - t422 * t476;
t161 = Icges(5,1) * t275 + Icges(5,4) * t274 + Icges(5,5) * t438;
t160 = Icges(5,1) * t273 + Icges(5,4) * t272 + Icges(5,5) * t437;
t159 = Icges(5,4) * t275 + Icges(5,2) * t274 + Icges(5,6) * t438;
t158 = Icges(5,4) * t273 + Icges(5,2) * t272 + Icges(5,6) * t437;
t154 = t469 * t422;
t153 = -t403 * t219 - t287 * t640;
t152 = t218 * t403 + t287 * t641;
t145 = -t283 * t403 + (-t284 * t400 + t648) * t401;
t142 = t145 * t599;
t141 = qJD(1) * t228 + t422 * t465;
t140 = t296 * t601 + t424 * t465 + t615;
t138 = t482 * t401;
t136 = -t207 * t641 + t186;
t132 = t456 * t422;
t131 = rSges(6,3) * t438 + t523;
t130 = -rSges(6,3) * t567 + t575;
t109 = t372 + (t598 * t667 - t593) * t422 + t439 * t600 + t526 + t568;
t108 = t367 + (-t682 - t683) * t594 + t425 * qJD(1) + t571 + t609;
t106 = -t424 * t585 + (qJD(2) * t472 - t541) * t422 + ((t365 - t680) * t422 + t446 * t424) * qJD(1) + t569;
t103 = t341 * t236 + t342 * t238 + t579;
t102 = t341 * t235 + t342 * t237 + t581;
t101 = -t236 * t458 + t238 * t340 + t580;
t100 = -t235 * t458 + t237 * t340 + t582;
t97 = t422 * t245 + t246 * t424 + t537;
t86 = qJD(1) * t155 + t422 * t454;
t85 = t424 * t454 + t564 + t572;
t80 = -t401 * t548 - t403 * t622;
t79 = t169 * t403 + t247 * t641 + t148;
t78 = (-rSges(6,3) * t598 - t593) * t422 + (t424 * t459 - t391) * qJD(1) - t523 + t568 + t569;
t77 = (-t396 * t401 - t637 - t683) * t594 + t428 * qJD(1) + t575 + t609 + t612;
t76 = t401 * t447 + t186;
t75 = (qJD(1) * t435 + t585) * t424 + (t541 - qJD(1) * t365 - t593 + (t403 * t666 + t643) * qJD(2)) * t422 - t520 + t568;
t74 = (t545 - t668) * t601 + (-t541 - t602 + (t472 - t683) * qJD(2)) * t424 + t570 + t576 + t609;
t73 = t460 - t481;
t72 = (t287 * t596 + t131) * t403 + (-qJD(2) * t218 + t422 * t211 + t287 * t600) * t401;
t71 = (-t287 * t594 - t130) * t403 + (qJD(2) * t219 - t211 * t424 + t564) * t401;
t70 = qJD(1) * t133 + t422 * t436;
t69 = t424 * t436 + t601 * t619 + t572;
t68 = -t206 * t599 + t538;
t67 = -t190 * t640 + (-t279 * t594 - t116) * t403 + t620;
t66 = (-t562 + (-qJD(5) * t285 - t209) * t401) * t400 + t429;
t57 = t460 - t703;
t48 = t482 * t598 + (qJD(1) * t481 - t130 * t422 + t131 * t424) * t401;
t46 = t422 * (rSges(5,3) * t559 - t526) + t424 * t571 + (t424 * t245 + (-t246 + t617) * t422) * qJD(1) + t461;
t45 = -t135 * t403 + t401 * t513;
t44 = -t134 * t403 + t401 * t514;
t43 = -t207 * t559 + (-t116 * t422 + (-t422 * t206 - t207 * t424) * qJD(1)) * t401 + t625;
t33 = (t247 * t596 + t106) * t403 + (-qJD(2) * t623 + t422 * t193 + t247 * t600) * t401 + t538;
t32 = (-qJD(2) * t548 - t624) * t403 + (qJD(2) * t170 + t247 * t601 + t424 * t621) * t401 + t620;
t31 = t124 * t641 + t321 * t126 + t322 * t128 + t198 * t215 + t199 * t217 + t213 * t438;
t30 = t125 * t641 + t321 * t127 + t322 * t129 + t198 * t214 + t199 * t216 + t212 * t438;
t29 = t124 * t640 + t323 * t126 + t324 * t128 + t196 * t215 + t197 * t217 + t213 * t437;
t28 = t125 * t640 + t323 * t127 + t324 * t129 + t196 * t214 + t197 * t216 + t212 * t437;
t19 = t130 * t424 + t422 * t131 + (t650 + (-t219 + t574) * t422) * qJD(1) + t433;
t18 = t447 * t598 + (qJD(1) * t703 + t106 * t424 - t624 * t422) * t401 + t625;
t17 = t624 * t424 + (t106 + t117) * t422 + (t623 * t424 + (t574 - t622) * t422) * qJD(1) + t433;
t16 = qJD(1) * t514 - t30 * t424 + t31 * t422;
t15 = qJD(1) * t513 - t28 * t424 + t29 * t422;
t8 = (qJD(2) * t514 - t59) * t403 + (-qJD(1) * t63 + qJD(2) * t134 + t30 * t422 + t31 * t424) * t401;
t7 = (qJD(2) * t513 - t58) * t403 + (-qJD(1) * t64 + qJD(2) * t135 + t28 * t422 + t29 * t424) * t401;
t1 = [(t143 * t75 + t144 * t74) * t698 + (t150 * t78 + t151 * t77) * t699 + (t108 * t167 + t109 * t166) * t700 + (t171 * t266 + t172 * t265) * t701 + (t231 * t292 + t232 * t291) * t702 + t429 + t430 - t400 * t562 - t285 * t557 + (-Icges(3,2) * t423 + t509 - t665) * t597 + (Icges(3,1) * t421 + t502 + t664) * t595 + (-t278 * t642 - t563) * t394 + (-Icges(5,3) * t403 - t499 + t507) * t599 + (-t188 * t394 - t209 * t400 - t281 * t417 + t282 * t418 + t493 * t599) * t401 + (-Icges(5,3) * t401 - t294 * t417 + t295 * t418 - t403 * t493 + t500 + t506) * t598; (-t274 * t294 / 0.2e1 - t275 * t295 / 0.2e1 + t458 * t281 / 0.2e1 - t340 * t282 / 0.2e1 - t52 / 0.2e1 - t26 / 0.2e1 + t554 * t424 - t584) * t424 + (t51 / 0.2e1 + t27 / 0.2e1 + t272 * t694 + t273 * t693 + t281 * t692 + t282 * t691 + t554 * t422 + t583) * t422 + m(3) * ((-t231 * t422 - t232 * t424) * t378 + (-t291 * t424 - t292 * t422) * t364) + m(4) * (t171 * t299 + t172 * t300 + t223 * t265 + t224 * t266) + m(5) * (t108 * t227 + t109 * t228 + t140 * t166 + t141 * t167) + m(6) * (t150 * t85 + t151 * t86 + t154 * t77 + t155 * t78) + m(7) * (t132 * t74 + t133 * t75 + t143 * t69 + t144 * t70) + ((t310 * t715 + t442 * t687 + Icges(5,5) * t275 / 0.2e1 + Icges(5,6) * t274 / 0.2e1 + Icges(5,3) * t438 / 0.2e1) * t424 + (t309 * t715 + t442 * t686 - Icges(5,5) * t273 / 0.2e1 - Icges(5,6) * t272 / 0.2e1 - Icges(5,3) * t437 / 0.2e1) * t422) * t403 + ((-t158 * t417 + t160 * t418 - t424 * t444) * t687 + (-t159 * t417 + t161 * t418 - t422 * t444) * t686) * t401 + (t449 / 0.2e1 + t451 / 0.2e1 + (t233 * t401 - t235 * t639 + t237 * t638) * t686 - t450 / 0.2e1 + (t234 * t401 - t236 * t639 + t238 * t638) * t687 - t448 / 0.2e1) * qJD(2) + ((-t311 * t687 + t312 * t686) * t401 + (t123 / 0.2e1 - t292 * t684 + t294 * t692 + t295 * t691 + t96 / 0.2e1 + (-t234 / 0.2e1 + t310 / 0.2e1) * t403 + (t236 * t689 + t238 * t688 + t312 / 0.2e1) * t401 + t577) * t424 + (t291 * t684 - t458 * t694 + t340 * t693 + t122 / 0.2e1 + t95 / 0.2e1 + (-t233 / 0.2e1 + t309 / 0.2e1) * t403 + (t235 * t689 + t237 * t688 + t311 / 0.2e1) * t401 + t578) * t422) * qJD(1); (t132 * t70 + t133 * t69 + t17 * t57) * t698 + (t154 * t86 + t155 * t85 + t19 * t73) * t699 + (t140 * t228 + t141 * t227 + t46 * t97) * t700 + (t300 * t223 + t299 * t224 + (t422 * t317 + t318 * t424 + t618) * ((qJD(1) * t317 - t424 * t452 + t610) * t424 + (-t422 * t452 + (-t306 - t318 + t708) * qJD(1)) * t422 + t573)) * t701 + ((t422 * t336 + t337 * t424) * ((qJD(1) * t336 + t431) * t424 + (-t422 * t453 + (-t337 - t590 + t409) * qJD(1)) * t422) + t606 * t378 * t364) * t702 + t422 * ((t422 * t257 + (t178 + t448) * qJD(1)) * t422 + (t179 * qJD(1) + (t330 * t595 + t332 * t597) * t424 + (-t258 + (-t331 * t423 - t333 * t421) * qJD(2) + (t329 - t474) * qJD(1)) * t422) * t424) - t424 * ((t424 * t258 + (t177 + t449) * qJD(1)) * t424 + (t176 * qJD(1) + (-t331 * t595 - t333 * t597 + t604) * t422 + (-t257 + (t330 * t423 + t332 * t421) * qJD(2) - t473 * qJD(1)) * t424) * t422) - t424 * ((t424 * t240 + (t163 + t451) * qJD(1)) * t424 + (t162 * qJD(1) + (-t310 * t598 - t312 * t599 + t605) * t422 + (-t239 + (t309 * t403 + t311 * t401) * qJD(2) - t475 * qJD(1)) * t424) * t422) + t422 * ((t422 * t239 + (t164 + t450) * qJD(1)) * t422 + (t165 * qJD(1) + (t309 * t598 + t311 * t599) * t424 + (-t240 + (-t310 * t403 - t312 * t401) * qJD(2) + (t308 - t476) * qJD(1)) * t422) * t424) - t424 * ((t458 * t159 - t340 * t161 - t274 * t235 - t275 * t237 + (t101 - t581) * qJD(1)) * t424 + (-t458 * t158 + t340 * t160 + t274 * t236 + t275 * t238 + (t100 + t579) * qJD(1)) * t422) + t422 * ((t341 * t158 + t342 * t160 + t272 * t236 + t273 * t238 + (t102 - t580) * qJD(1)) * t422 + (-t341 * t159 - t342 * t161 - t272 * t235 - t273 * t237 + (t103 + t582) * qJD(1)) * t424) - t424 * t16 - t424 * t13 + t422 * t15 + t422 * t12 + (t60 + t63 + (-t100 - t162 - t176) * t424 + (t101 + t163 + t177) * t422) * t601 + (t61 + t64 + (-t102 - t164 - t178) * t424 + (t103 + t165 + t179) * t422) * t600; m(7) * (qJD(1) * t490 + t422 * t75 - t424 * t74) + m(6) * (qJD(1) * t488 + t422 * t78 - t424 * t77) + m(5) * (qJD(1) * t487 - t108 * t424 + t422 * t109) + m(4) * (-t171 * t424 + t422 * t172 + (t265 * t424 + t266 * t422) * qJD(1)); m(7) * (t422 * t69 - t424 * t70 + (t132 * t422 + t133 * t424) * qJD(1)) + m(6) * (t422 * t85 - t424 * t86 + (t154 * t422 + t155 * t424) * qJD(1)) + m(5) * (t422 * t140 - t141 * t424 + (t227 * t422 + t228 * t424) * qJD(1)) + m(4) * (t422 * t223 - t224 * t424 + (t299 * t422 + t300 * t424) * qJD(1)); 0; 0.2e1 * (t487 * t697 + t488 * t696 + t490 * t695) * t598 + 0.2e1 * ((-t143 * t601 + t144 * t600 + t422 * t74 + t424 * t75) * t695 + (-t150 * t601 + t151 * t600 + t422 * t77 + t424 * t78) * t696 + (t108 * t422 + t109 * t424 - t166 * t601 + t167 * t600) * t697) * t401; 0.2e1 * ((t132 * t596 + t133 * t594 - t17) * t695 + (t154 * t596 + t155 * t594 - t19) * t696 + (t227 * t596 + t228 * t594 - t46) * t697) * t403 + 0.2e1 * ((qJD(2) * t57 + t132 * t600 - t133 * t601 + t422 * t70 + t424 * t69) * t695 + (qJD(2) * t73 + t154 * t600 - t155 * t601 + t422 * t86 + t424 * t85) * t696 + (qJD(2) * t97 + t140 * t424 + t141 * t422 + t227 * t600 - t228 * t601) * t697) * t401; 0; 0.4e1 * (t697 + t696 + t695) * (-0.1e1 + t606) * t401 * t598; t426 + t142 + (t583 * t424 + t584 * t422 + (-t422 * t577 + t424 * t578) * qJD(1)) * t401 + m(7) * (t143 * t33 + t144 * t32 + t74 * t80 + t75 * t79) + m(6) * (t150 * t72 + t151 * t71 + t152 * t78 + t153 * t77) + (-t62 - t66 + (t422 * t578 + t424 * t577) * qJD(2)) * t403; t432 + ((t99 * t422 - t424 * t98) * t711 + t15 * t685 + t16 * t687 + (t63 * t685 - t422 * t64 / 0.2e1) * qJD(1)) * t401 + m(6) * (t138 * t19 + t152 * t85 + t153 * t86 + t154 * t71 + t155 * t72 + t48 * t73) + m(7) * (t132 * t32 + t133 * t33 + t17 * t76 + t18 * t57 + t69 * t79 + t70 * t80) + (t45 * t677 + (qJD(1) * t99 - t34) * t690 - t8 / 0.2e1 + t64 * t551) * t424 + ((qJD(1) * t98 + t35) * t690 + t44 * t677 + t7 / 0.2e1 + t63 * t551) * t422; m(6) * (t72 * t422 - t424 * t71 + (t152 * t424 + t153 * t422) * qJD(1)) + m(7) * (-t32 * t424 + t33 * t422 + (t422 * t80 + t424 * t79) * qJD(1)); 0.2e1 * ((t152 * t594 + t153 * t596 - t48) * t696 + (t594 * t79 + t596 * t80 - t18) * t695) * t403 + 0.2e1 * ((qJD(2) * t138 - t152 * t601 + t153 * t600 + t422 * t71 + t424 * t72) * t696 + (qJD(2) * t76 + t32 * t422 + t33 * t424 + t600 * t80 - t601 * t79) * t695) * t401; (t18 * t76 + t32 * t80 + t33 * t79) * t698 + (t138 * t48 + t152 * t72 + t153 * t71) * t699 + (t66 * t403 - t142 - t9 + (-t403 * t510 + t422 * t44 + t424 * t45) * qJD(2)) * t403 + (t424 * t7 + t422 * t8 - t403 * (t34 * t422 + t35 * t424) + (-t145 * t403 + t401 * t510) * qJD(2) + ((-t403 * t98 + t44) * t424 + (t403 * t99 - t42 - t45) * t422) * qJD(1)) * t401 + t529; t426 - t62 * t403 + m(7) * (t143 * t68 + t144 * t67 + t148 * t75 + t149 * t74); t432 + m(7) * (t132 * t67 + t133 * t68 + t136 * t17 + t148 * t69 + t149 * t70 + t43 * t57); m(7) * (qJD(1) * t489 + t68 * t422 - t424 * t67); m(7) * ((qJD(2) * t489 - t43) * t403 + (qJD(2) * t136 + t422 * t67 + t424 * t68 + (-t148 * t422 + t149 * t424) * qJD(1)) * t401); m(7) * (t136 * t18 + t148 * t33 + t149 * t32 + t43 * t76 + t67 * t80 + t68 * t79) + t427; (t136 * t43 + t148 * t68 + t149 * t67) * t698 + t427;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
