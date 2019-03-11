% Calculate time derivative of joint inertia matrix for
% S6PRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRP4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP4_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP4_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP4_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRP4_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:13:11
% EndTime: 2019-03-09 00:14:04
% DurationCPUTime: 32.65s
% Computational Cost: add. (138316->1377), mult. (301627->1851), div. (0->0), fcn. (346867->12), ass. (0->567)
t759 = rSges(7,1) + pkin(5);
t758 = rSges(7,3) + qJ(6);
t566 = sin(pkin(11));
t568 = cos(pkin(11));
t574 = cos(qJ(2));
t569 = cos(pkin(6));
t572 = sin(qJ(2));
t698 = t569 * t572;
t550 = t566 * t574 + t568 * t698;
t567 = sin(pkin(6));
t571 = sin(qJ(3));
t700 = t567 * t571;
t715 = cos(qJ(3));
t526 = t550 * t715 - t568 * t700;
t696 = qJ(4) + qJ(5);
t564 = sin(t696);
t697 = t569 * t574;
t590 = -t566 * t572 + t568 * t697;
t630 = cos(t696);
t478 = t526 * t564 + t590 * t630;
t479 = t526 * t630 - t564 * t590;
t636 = t567 * t715;
t585 = -t550 * t571 - t568 * t636;
t367 = Icges(7,5) * t479 - Icges(7,6) * t585 + Icges(7,3) * t478;
t371 = Icges(7,4) * t479 - Icges(7,2) * t585 + Icges(7,6) * t478;
t375 = Icges(7,1) * t479 - Icges(7,4) * t585 + Icges(7,5) * t478;
t198 = t367 * t478 - t371 * t585 + t375 * t479;
t369 = Icges(6,5) * t479 - Icges(6,6) * t478 - Icges(6,3) * t585;
t373 = Icges(6,4) * t479 - Icges(6,2) * t478 - Icges(6,6) * t585;
t377 = Icges(6,1) * t479 - Icges(6,4) * t478 - Icges(6,5) * t585;
t200 = -t369 * t585 - t373 * t478 + t377 * t479;
t757 = t198 + t200;
t655 = t566 * t698;
t552 = t568 * t574 - t655;
t528 = t552 * t715 + t566 * t700;
t551 = t566 * t697 + t568 * t572;
t480 = t528 * t564 - t551 * t630;
t481 = t528 * t630 + t551 * t564;
t586 = -t552 * t571 + t566 * t636;
t368 = Icges(7,5) * t481 - Icges(7,6) * t586 + Icges(7,3) * t480;
t372 = Icges(7,4) * t481 - Icges(7,2) * t586 + Icges(7,6) * t480;
t376 = Icges(7,1) * t481 - Icges(7,4) * t586 + Icges(7,5) * t480;
t199 = t368 * t478 - t372 * t585 + t376 * t479;
t370 = Icges(6,5) * t481 - Icges(6,6) * t480 - Icges(6,3) * t586;
t374 = Icges(6,4) * t481 - Icges(6,2) * t480 - Icges(6,6) * t586;
t378 = Icges(6,1) * t481 - Icges(6,4) * t480 - Icges(6,5) * t586;
t201 = -t370 * t585 - t374 * t478 + t378 * t479;
t756 = t199 + t201;
t202 = t367 * t480 - t371 * t586 + t375 * t481;
t204 = -t369 * t586 - t373 * t480 + t377 * t481;
t755 = t202 + t204;
t203 = t368 * t480 - t372 * t586 + t376 * t481;
t205 = -t370 * t586 - t374 * t480 + t378 * t481;
t754 = t203 + t205;
t554 = t569 * t571 + t572 * t636;
t599 = t567 * t630;
t596 = t574 * t599;
t523 = t554 * t564 + t596;
t699 = t567 * t574;
t653 = t564 * t699;
t524 = t554 * t630 - t653;
t553 = -t569 * t715 + t572 * t700;
t214 = t367 * t523 + t371 * t553 + t375 * t524;
t216 = t369 * t553 - t373 * t523 + t377 * t524;
t753 = t214 + t216;
t215 = t368 * t523 + t372 * t553 + t376 * t524;
t217 = t370 * t553 - t374 * t523 + t378 * t524;
t752 = t215 + t217;
t435 = Icges(7,5) * t524 + Icges(7,6) * t553 + Icges(7,3) * t523;
t437 = Icges(7,4) * t524 + Icges(7,2) * t553 + Icges(7,6) * t523;
t439 = Icges(7,1) * t524 + Icges(7,4) * t553 + Icges(7,5) * t523;
t226 = t435 * t478 - t437 * t585 + t439 * t479;
t436 = Icges(6,5) * t524 - Icges(6,6) * t523 + Icges(6,3) * t553;
t438 = Icges(6,4) * t524 - Icges(6,2) * t523 + Icges(6,6) * t553;
t440 = Icges(6,1) * t524 - Icges(6,4) * t523 + Icges(6,5) * t553;
t227 = -t436 * t585 - t438 * t478 + t440 * t479;
t751 = t226 + t227;
t228 = t435 * t480 - t437 * t586 + t439 * t481;
t229 = -t436 * t586 - t438 * t480 + t440 * t481;
t750 = t228 + t229;
t270 = t435 * t523 + t437 * t553 + t439 * t524;
t271 = t436 * t553 - t438 * t523 + t440 * t524;
t749 = t270 + t271;
t658 = qJD(2) * t574;
t634 = t567 * t658;
t530 = -qJD(3) * t553 + t634 * t715;
t565 = qJD(4) + qJD(5);
t600 = t565 * t630;
t659 = qJD(2) * t572;
t444 = t530 * t564 + t554 * t600 - t565 * t653 - t599 * t659;
t635 = t567 * t659;
t445 = -t565 * t596 + t530 * t630 + (-t554 * t565 + t635) * t564;
t529 = qJD(3) * t554 + t571 * t634;
t318 = Icges(7,5) * t445 + Icges(7,6) * t529 + Icges(7,3) * t444;
t320 = Icges(7,4) * t445 + Icges(7,2) * t529 + Icges(7,6) * t444;
t322 = Icges(7,1) * t445 + Icges(7,4) * t529 + Icges(7,5) * t444;
t540 = t590 * qJD(2);
t484 = qJD(3) * t585 + t540 * t715;
t541 = t550 * qJD(2);
t386 = t526 * t600 - t541 * t630 + (-t565 * t590 + t484) * t564;
t387 = -t590 * t600 + t484 * t630 + (-t526 * t565 + t541) * t564;
t483 = qJD(3) * t526 + t540 * t571;
t151 = t318 * t478 - t320 * t585 + t322 * t479 + t386 * t435 + t387 * t439 + t437 * t483;
t319 = Icges(6,5) * t445 - Icges(6,6) * t444 + Icges(6,3) * t529;
t321 = Icges(6,4) * t445 - Icges(6,2) * t444 + Icges(6,6) * t529;
t323 = Icges(6,1) * t445 - Icges(6,4) * t444 + Icges(6,5) * t529;
t152 = -t319 * t585 - t321 * t478 + t323 * t479 - t386 * t438 + t387 * t440 + t436 * t483;
t542 = t551 * qJD(2);
t485 = qJD(3) * t528 - t542 * t571;
t254 = Icges(7,5) * t387 + Icges(7,6) * t483 + Icges(7,3) * t386;
t258 = Icges(7,4) * t387 + Icges(7,2) * t483 + Icges(7,6) * t386;
t262 = Icges(7,1) * t387 + Icges(7,4) * t483 + Icges(7,5) * t386;
t73 = t254 * t478 - t258 * t585 + t262 * t479 + t367 * t386 + t371 * t483 + t375 * t387;
t486 = qJD(3) * t586 - t542 * t715;
t543 = -qJD(2) * t655 + t568 * t658;
t388 = t528 * t600 - t543 * t630 + (t551 * t565 + t486) * t564;
t389 = t551 * t600 + t486 * t630 + (-t528 * t565 + t543) * t564;
t255 = Icges(7,5) * t389 + Icges(7,6) * t485 + Icges(7,3) * t388;
t259 = Icges(7,4) * t389 + Icges(7,2) * t485 + Icges(7,6) * t388;
t263 = Icges(7,1) * t389 + Icges(7,4) * t485 + Icges(7,5) * t388;
t74 = t255 * t478 - t259 * t585 + t263 * t479 + t368 * t386 + t372 * t483 + t376 * t387;
t256 = Icges(6,5) * t387 - Icges(6,6) * t386 + Icges(6,3) * t483;
t260 = Icges(6,4) * t387 - Icges(6,2) * t386 + Icges(6,6) * t483;
t264 = Icges(6,1) * t387 - Icges(6,4) * t386 + Icges(6,5) * t483;
t75 = -t256 * t585 - t260 * t478 + t264 * t479 + t369 * t483 - t373 * t386 + t377 * t387;
t257 = Icges(6,5) * t389 - Icges(6,6) * t388 + Icges(6,3) * t485;
t261 = Icges(6,4) * t389 - Icges(6,2) * t388 + Icges(6,6) * t485;
t265 = Icges(6,1) * t389 - Icges(6,4) * t388 + Icges(6,5) * t485;
t76 = -t257 * t585 - t261 * t478 + t265 * t479 + t370 * t483 - t374 * t386 + t378 * t387;
t748 = (-t74 - t76) * t586 + (-t73 - t75) * t585 + (t151 + t152) * t553 + t751 * t529 + t756 * t485 + t757 * t483;
t153 = t318 * t480 - t320 * t586 + t322 * t481 + t388 * t435 + t389 * t439 + t437 * t485;
t154 = -t319 * t586 - t321 * t480 + t323 * t481 - t388 * t438 + t389 * t440 + t436 * t485;
t77 = t254 * t480 - t258 * t586 + t262 * t481 + t367 * t388 + t371 * t485 + t375 * t389;
t78 = t255 * t480 - t259 * t586 + t263 * t481 + t368 * t388 + t372 * t485 + t376 * t389;
t79 = -t256 * t586 - t260 * t480 + t264 * t481 + t369 * t485 - t373 * t388 + t377 * t389;
t80 = -t257 * t586 - t261 * t480 + t265 * t481 + t370 * t485 - t374 * t388 + t378 * t389;
t747 = (-t78 - t80) * t586 + (-t77 - t79) * t585 + (t153 + t154) * t553 + t750 * t529 + t754 * t485 + t755 * t483;
t21 = t198 * t541 + t199 * t543 - t590 * t73 + t551 * t74 + (-t151 * t574 + t226 * t659) * t567;
t22 = t200 * t541 + t201 * t543 - t590 * t75 + t551 * t76 + (-t152 * t574 + t227 * t659) * t567;
t746 = t21 + t22;
t23 = t202 * t541 + t203 * t543 - t590 * t77 + t551 * t78 + (-t153 * t574 + t228 * t659) * t567;
t24 = t204 * t541 + t205 * t543 - t590 * t79 + t551 * t80 + (-t154 * t574 + t229 * t659) * t567;
t745 = t23 + t24;
t162 = t318 * t523 + t320 * t553 + t322 * t524 + t435 * t444 + t437 * t529 + t439 * t445;
t163 = t319 * t553 - t321 * t523 + t323 * t524 + t436 * t529 - t438 * t444 + t440 * t445;
t92 = t254 * t523 + t258 * t553 + t262 * t524 + t367 * t444 + t371 * t529 + t375 * t445;
t93 = t255 * t523 + t259 * t553 + t263 * t524 + t368 * t444 + t372 * t529 + t376 * t445;
t94 = t256 * t553 - t260 * t523 + t264 * t524 + t369 * t529 - t373 * t444 + t377 * t445;
t95 = t257 * t553 - t261 * t523 + t265 * t524 + t370 * t529 - t374 * t444 + t378 * t445;
t744 = (-t93 - t95) * t586 + (-t92 - t94) * t585 + (t162 + t163) * t553 + t749 * t529 + t752 * t485 + t753 * t483;
t39 = t214 * t541 + t215 * t543 - t590 * t92 + t551 * t93 + (-t162 * t574 + t270 * t659) * t567;
t40 = t216 * t541 + t217 * t543 - t590 * t94 + t551 * t95 + (-t163 * t574 + t271 * t659) * t567;
t743 = t39 + t40;
t742 = t553 * t751 - t585 * t757 - t586 * t756;
t741 = t553 * t750 - t585 * t755 - t586 * t754;
t740 = t551 * t756 - t590 * t757 - t699 * t751;
t739 = t551 * t754 - t590 * t755 - t699 * t750;
t738 = t553 * t749 - t585 * t753 - t586 * t752;
t737 = t551 * t752 - t590 * t753 - t699 * t749;
t692 = rSges(7,2) * t485 + qJD(6) * t480 + t758 * t388 + t389 * t759;
t686 = rSges(7,2) * t529 + qJD(6) * t523 + t758 * t444 + t445 * t759;
t379 = rSges(7,1) * t479 - rSges(7,2) * t585 + rSges(7,3) * t478;
t424 = pkin(5) * t479 + qJ(6) * t478;
t676 = t379 + t424;
t675 = -rSges(7,2) * t586 + t758 * t480 + t481 * t759;
t668 = rSges(7,2) * t553 + t758 * t523 + t524 * t759;
t269 = rSges(6,1) * t389 - rSges(6,2) * t388 + rSges(6,3) * t485;
t382 = rSges(6,1) * t481 - rSges(6,2) * t480 - rSges(6,3) * t586;
t267 = rSges(6,1) * t387 - rSges(6,2) * t386 + rSges(6,3) * t483;
t380 = rSges(6,1) * t479 - rSges(6,2) * t478 - rSges(6,3) * t585;
t695 = -t267 * t586 + t485 * t380;
t161 = t269 * t585 - t483 * t382 + t695;
t517 = pkin(2) * t540 + pkin(8) * t541;
t518 = -pkin(2) * t542 + pkin(8) * t543;
t701 = t567 * t568;
t702 = t566 * t567;
t662 = t517 * t702 + t518 * t701;
t727 = m(7) / 0.2e1;
t728 = m(6) / 0.2e1;
t735 = t728 + t727;
t573 = cos(qJ(4));
t570 = sin(qJ(4));
t703 = t551 * t570;
t490 = t528 * t573 + t703;
t404 = -qJD(4) * t490 - t486 * t570 + t543 * t573;
t489 = -t528 * t570 + t551 * t573;
t581 = qJD(4) * t489 + t543 * t570;
t405 = t486 * t573 + t581;
t283 = rSges(5,1) * t405 + rSges(5,2) * t404 + rSges(5,3) * t485;
t414 = rSges(4,1) * t486 - rSges(4,2) * t485 + rSges(4,3) * t543;
t734 = -m(4) * t414 - m(5) * t283 - m(6) * t269 - m(7) * t692;
t733 = -0.2e1 * t483;
t732 = 0.2e1 * t585;
t731 = -0.2e1 * t541;
t730 = 0.2e1 * t590;
t729 = m(5) / 0.2e1;
t726 = t483 / 0.2e1;
t725 = t485 / 0.2e1;
t724 = -t585 / 0.2e1;
t723 = -t586 / 0.2e1;
t722 = t529 / 0.2e1;
t721 = t541 / 0.2e1;
t720 = t543 / 0.2e1;
t719 = -t590 / 0.2e1;
t718 = t551 / 0.2e1;
t717 = t553 / 0.2e1;
t716 = t569 / 0.2e1;
t714 = pkin(4) * t573;
t712 = Icges(3,4) * t572;
t711 = Icges(3,4) * t574;
t449 = Icges(4,5) * t526 + Icges(4,6) * t585 - Icges(4,3) * t590;
t451 = Icges(4,4) * t526 + Icges(4,2) * t585 - Icges(4,6) * t590;
t453 = Icges(4,1) * t526 + Icges(4,4) * t585 - Icges(4,5) * t590;
t288 = -t449 * t590 + t451 * t585 + t453 * t526;
t708 = t541 * t288;
t503 = Icges(4,5) * t554 - Icges(4,6) * t553 - Icges(4,3) * t699;
t504 = Icges(4,4) * t554 - Icges(4,2) * t553 - Icges(4,6) * t699;
t505 = Icges(4,1) * t554 - Icges(4,4) * t553 - Icges(4,5) * t699;
t314 = -t503 * t590 + t504 * t585 + t505 * t526;
t707 = t541 * t314;
t450 = Icges(4,5) * t528 + Icges(4,6) * t586 + Icges(4,3) * t551;
t452 = Icges(4,4) * t528 + Icges(4,2) * t586 + Icges(4,6) * t551;
t454 = Icges(4,1) * t528 + Icges(4,4) * t586 + Icges(4,5) * t551;
t291 = t450 * t551 + t452 * t586 + t454 * t528;
t706 = t543 * t291;
t315 = t503 * t551 + t504 * t586 + t505 * t528;
t705 = t543 * t315;
t704 = t590 * t570;
t694 = t553 * t269 + t529 * t382;
t252 = pkin(5) * t387 + qJ(6) * t386 + qJD(6) * t478;
t266 = rSges(7,1) * t387 + rSges(7,2) * t483 + rSges(7,3) * t386;
t693 = t252 + t266;
t300 = pkin(4) * t581 + pkin(10) * t485 + t486 * t714;
t691 = -t269 - t300;
t487 = -t526 * t570 - t573 * t590;
t582 = qJD(4) * t487 + t541 * t570;
t299 = pkin(4) * t582 + pkin(10) * t483 + t484 * t714;
t275 = t586 * t299;
t365 = -pkin(4) * t704 - pkin(10) * t585 + t526 * t714;
t324 = t485 * t365;
t690 = -t275 + t324;
t431 = t486 * pkin(3) + t485 * pkin(9);
t689 = -t283 - t431;
t366 = pkin(4) * t703 - pkin(10) * t586 + t528 * t714;
t688 = t553 * t300 + t529 * t366;
t326 = rSges(6,1) * t445 - rSges(6,2) * t444 + rSges(6,3) * t529;
t442 = rSges(6,1) * t524 - rSges(6,2) * t523 + rSges(6,3) * t553;
t687 = -t326 * t585 + t483 * t442;
t531 = -t554 * t570 - t573 * t699;
t578 = qJD(4) * t531 + t570 * t635;
t351 = pkin(4) * t578 + pkin(10) * t529 + t530 * t714;
t685 = -t326 - t351;
t654 = t570 * t699;
t443 = -pkin(4) * t654 + pkin(10) * t553 + t554 * t714;
t684 = -t351 * t585 + t483 * t443;
t591 = -t554 * t573 + t654;
t463 = qJD(4) * t591 - t530 * t570 + t573 * t635;
t464 = t530 * t573 + t578;
t341 = rSges(5,1) * t464 + rSges(5,2) * t463 + rSges(5,3) * t529;
t477 = t530 * pkin(3) + t529 * pkin(9);
t683 = -t341 - t477;
t682 = t676 * t586;
t475 = t526 * pkin(3) - pkin(9) * t585;
t460 = t551 * t475;
t681 = t551 * t365 + t460;
t680 = t675 * t553;
t476 = t528 * pkin(3) - pkin(9) * t586;
t469 = t476 * t635;
t679 = t366 * t635 + t469;
t678 = -t365 - t380;
t677 = -t366 - t382;
t488 = t526 * t573 - t704;
t400 = rSges(5,1) * t488 + rSges(5,2) * t487 - rSges(5,3) * t585;
t674 = -t400 - t475;
t401 = rSges(5,1) * t490 + rSges(5,2) * t489 - rSges(5,3) * t586;
t673 = -t401 - t476;
t430 = t484 * pkin(3) + t483 * pkin(9);
t417 = t551 * t430;
t448 = t543 * t475;
t672 = t417 + t448;
t671 = t668 * t585;
t496 = t569 * t518;
t670 = t569 * t431 + t496;
t669 = -t430 - t517;
t667 = -t442 - t443;
t461 = -rSges(5,1) * t591 + rSges(5,2) * t531 + rSges(5,3) * t553;
t522 = t554 * pkin(3) + t553 * pkin(9);
t666 = -t461 - t522;
t665 = t475 * t699 - t522 * t590;
t521 = pkin(2) * t552 + pkin(8) * t551;
t519 = t569 * t521;
t664 = t569 * t476 + t519;
t663 = 0.2e1 * t662;
t520 = pkin(2) * t550 - pkin(8) * t590;
t661 = t520 * t702 + t521 * t701;
t660 = qJD(2) * t567;
t651 = -t300 - t692;
t650 = -t431 + t691;
t649 = t569 * t300 + t670;
t648 = -t299 + t669;
t647 = -t351 - t686;
t646 = -t477 + t685;
t645 = t569 * t366 + t664;
t644 = -t365 - t676;
t643 = -t475 + t678;
t642 = -t366 - t675;
t641 = -t476 + t677;
t640 = t430 * t699 - t477 * t590 + t541 * t522;
t639 = -t443 - t668;
t638 = -t522 + t667;
t629 = t659 / 0.2e1;
t628 = 0.2e1 * m(4);
t626 = 0.2e1 * m(5);
t624 = 0.2e1 * m(6);
t622 = 0.2e1 * m(7);
t468 = rSges(4,1) * t530 - rSges(4,2) * t529 + rSges(4,3) * t635;
t548 = (pkin(2) * t574 + pkin(8) * t572) * t660;
t621 = t567 * (-t468 - t548);
t506 = rSges(4,1) * t554 - rSges(4,2) * t553 - rSges(4,3) * t699;
t555 = (pkin(2) * t572 - pkin(8) * t574) * t567;
t620 = (-t506 - t555) * t567;
t234 = t586 * t252;
t239 = t586 * t266;
t331 = t485 * t379;
t350 = t485 * t424;
t619 = -t234 - t239 + t331 + t350;
t618 = t529 * t675 + t553 * t692;
t617 = -t431 + t651;
t616 = t300 * t732 + t366 * t733 - 0.2e1 * t275 + 0.2e1 * t324;
t286 = t551 * t299;
t355 = t543 * t365;
t615 = t286 + t355 + t672;
t614 = t483 * t668 - t585 * t686;
t613 = -t477 + t647;
t612 = 0.2e1 * t161;
t611 = t365 * t699 - t443 * t590 + t665;
t610 = -t475 + t644;
t609 = -t476 + t642;
t608 = t431 * t730 + t476 * t731 + 0.2e1 * t417 + 0.2e1 * t448;
t421 = t430 * t702;
t422 = t431 * t701;
t607 = 0.2e1 * t421 + 0.2e1 * t422 + t663;
t606 = t421 + t422 + t662;
t605 = -t522 + t639;
t604 = t475 * t702 + t476 * t701 + t661;
t602 = t567 * (-t548 + t683);
t601 = (-t555 + t666) * t567;
t598 = (-t548 + t646) * t567;
t597 = (-t555 + t638) * t567;
t294 = t299 * t702;
t295 = t300 * t701;
t594 = t294 + t295 + t606;
t593 = t299 * t699 - t351 * t590 + t541 * t443 + t640;
t592 = t365 * t702 + t366 * t701 + t604;
t589 = (-t548 + t613) * t567;
t588 = (-t555 + t605) * t567;
t583 = t675 * t733 + t692 * t732 - 0.2e1 * t234 - 0.2e1 * t239 + 0.2e1 * t331 + 0.2e1 * t350;
t580 = t742 * t483 + t741 * t485 + t738 * t529 + t744 * t553 - t585 * t748 - t747 * t586;
t402 = -qJD(4) * t488 - t484 * t570 + t541 * t573;
t403 = t484 * t573 + t582;
t282 = rSges(5,1) * t403 + rSges(5,2) * t402 + rSges(5,3) * t483;
t164 = -t282 * t586 + t283 * t585 + t400 * t485 - t401 * t483;
t413 = rSges(4,1) * t484 - rSges(4,2) * t483 + rSges(4,3) * t541;
t579 = m(4) * t413 + m(5) * t282 + m(6) * t267 + m(7) * t693;
t125 = t226 * t569 + (-t198 * t568 + t199 * t566) * t567;
t126 = t227 * t569 + (-t200 * t568 + t201 * t566) * t567;
t127 = t228 * t569 + (-t202 * t568 + t203 * t566) * t567;
t128 = t229 * t569 + (-t204 * t568 + t205 * t566) * t567;
t148 = t270 * t569 + (-t214 * t568 + t215 * t566) * t567;
t149 = t271 * t569 + (-t216 * t568 + t217 * t566) * t567;
t47 = t151 * t569 + (t566 * t74 - t568 * t73) * t567;
t48 = t152 * t569 + (t566 * t76 - t568 * t75) * t567;
t49 = t153 * t569 + (t566 * t78 - t568 * t77) * t567;
t50 = t154 * t569 + (t566 * t80 - t568 * t79) * t567;
t55 = t162 * t569 + (t566 * t93 - t568 * t92) * t567;
t56 = t163 * t569 + (t566 * t95 - t568 * t94) * t567;
t577 = (t125 + t126) * t726 + (t127 + t128) * t725 + (t47 + t48) * t724 + (t49 + t50) * t723 + (t148 + t149) * t722 + (t55 + t56) * t717 + t744 * t716 + t747 * t702 / 0.2e1 - t748 * t701 / 0.2e1;
t576 = t740 * t726 + t739 * t725 + t746 * t724 + t745 * t723 + t737 * t722 + t742 * t721 + t741 * t720 + t748 * t719 + t747 * t718 + t743 * t717 - t744 * t699 / 0.2e1 + t738 * t567 * t629;
t547 = (rSges(3,1) * t574 - rSges(3,2) * t572) * t660;
t546 = (Icges(3,1) * t574 - t712) * t660;
t545 = (-Icges(3,2) * t572 + t711) * t660;
t544 = (Icges(3,5) * t574 - Icges(3,6) * t572) * t660;
t539 = rSges(3,3) * t569 + (rSges(3,1) * t572 + rSges(3,2) * t574) * t567;
t538 = Icges(3,5) * t569 + (Icges(3,1) * t572 + t711) * t567;
t537 = Icges(3,6) * t569 + (Icges(3,2) * t574 + t712) * t567;
t516 = -rSges(3,1) * t542 - rSges(3,2) * t543;
t515 = rSges(3,1) * t540 - rSges(3,2) * t541;
t514 = -Icges(3,1) * t542 - Icges(3,4) * t543;
t513 = Icges(3,1) * t540 - Icges(3,4) * t541;
t512 = -Icges(3,4) * t542 - Icges(3,2) * t543;
t511 = Icges(3,4) * t540 - Icges(3,2) * t541;
t510 = -Icges(3,5) * t542 - Icges(3,6) * t543;
t509 = Icges(3,5) * t540 - Icges(3,6) * t541;
t502 = rSges(3,1) * t552 - rSges(3,2) * t551 + rSges(3,3) * t702;
t501 = rSges(3,1) * t550 + rSges(3,2) * t590 - rSges(3,3) * t701;
t500 = Icges(3,1) * t552 - Icges(3,4) * t551 + Icges(3,5) * t702;
t499 = Icges(3,1) * t550 + Icges(3,4) * t590 - Icges(3,5) * t701;
t498 = Icges(3,4) * t552 - Icges(3,2) * t551 + Icges(3,6) * t702;
t497 = Icges(3,4) * t550 + Icges(3,2) * t590 - Icges(3,6) * t701;
t467 = Icges(4,1) * t530 - Icges(4,4) * t529 + Icges(4,5) * t635;
t466 = Icges(4,4) * t530 - Icges(4,2) * t529 + Icges(4,6) * t635;
t465 = Icges(4,5) * t530 - Icges(4,6) * t529 + Icges(4,3) * t635;
t459 = -Icges(5,1) * t591 + Icges(5,4) * t531 + Icges(5,5) * t553;
t458 = -Icges(5,4) * t591 + Icges(5,2) * t531 + Icges(5,6) * t553;
t457 = -Icges(5,5) * t591 + Icges(5,6) * t531 + Icges(5,3) * t553;
t456 = rSges(4,1) * t528 + rSges(4,2) * t586 + rSges(4,3) * t551;
t455 = rSges(4,1) * t526 + rSges(4,2) * t585 - rSges(4,3) * t590;
t428 = t585 * t443;
t427 = t585 * t442;
t411 = Icges(4,1) * t486 - Icges(4,4) * t485 + Icges(4,5) * t543;
t410 = Icges(4,1) * t484 - Icges(4,4) * t483 + Icges(4,5) * t541;
t409 = Icges(4,4) * t486 - Icges(4,2) * t485 + Icges(4,6) * t543;
t408 = Icges(4,4) * t484 - Icges(4,2) * t483 + Icges(4,6) * t541;
t407 = Icges(4,5) * t486 - Icges(4,6) * t485 + Icges(4,3) * t543;
t406 = Icges(4,5) * t484 - Icges(4,6) * t483 + Icges(4,3) * t541;
t398 = Icges(5,1) * t490 + Icges(5,4) * t489 - Icges(5,5) * t586;
t397 = Icges(5,1) * t488 + Icges(5,4) * t487 - Icges(5,5) * t585;
t396 = Icges(5,4) * t490 + Icges(5,2) * t489 - Icges(5,6) * t586;
t395 = Icges(5,4) * t488 + Icges(5,2) * t487 - Icges(5,6) * t585;
t394 = Icges(5,5) * t490 + Icges(5,6) * t489 - Icges(5,3) * t586;
t393 = Icges(5,5) * t488 + Icges(5,6) * t487 - Icges(5,3) * t585;
t384 = -t456 * t699 - t506 * t551;
t383 = t455 * t699 - t506 * t590;
t359 = t553 * t382;
t357 = t553 * t366;
t354 = -t503 * t699 - t504 * t553 + t505 * t554;
t345 = t586 * t380;
t342 = t586 * t365;
t340 = Icges(5,1) * t464 + Icges(5,4) * t463 + Icges(5,5) * t529;
t339 = Icges(5,4) * t464 + Icges(5,2) * t463 + Icges(5,6) * t529;
t338 = Icges(5,5) * t464 + Icges(5,6) * t463 + Icges(5,3) * t529;
t336 = t455 * t551 + t456 * t590;
t335 = (-t455 - t520) * t569 + t568 * t620;
t334 = t456 * t569 + t566 * t620 + t519;
t309 = (t455 * t566 + t456 * t568) * t567 + t661;
t308 = t401 * t553 + t461 * t586;
t307 = -t400 * t553 - t461 * t585;
t306 = t442 * t586 + t359;
t305 = -t380 * t553 - t427;
t304 = -t450 * t699 - t452 * t553 + t454 * t554;
t303 = -t449 * t699 - t451 * t553 + t453 * t554;
t302 = (-t413 - t517) * t569 + t568 * t621;
t301 = t414 * t569 + t566 * t621 + t496;
t297 = t457 * t553 + t458 * t531 - t459 * t591;
t290 = t449 * t551 + t451 * t586 + t453 * t528;
t289 = -t450 * t590 + t452 * t585 + t454 * t526;
t281 = Icges(5,1) * t405 + Icges(5,4) * t404 + Icges(5,5) * t485;
t280 = Icges(5,1) * t403 + Icges(5,4) * t402 + Icges(5,5) * t483;
t279 = Icges(5,4) * t405 + Icges(5,2) * t404 + Icges(5,6) * t485;
t278 = Icges(5,4) * t403 + Icges(5,2) * t402 + Icges(5,6) * t483;
t277 = Icges(5,5) * t405 + Icges(5,6) * t404 + Icges(5,3) * t485;
t276 = Icges(5,5) * t403 + Icges(5,6) * t402 + Icges(5,3) * t483;
t272 = -t400 * t586 + t401 * t585;
t251 = t382 * t585 - t345;
t250 = t551 * t666 + t673 * t699;
t249 = t400 * t699 - t461 * t590 + t665;
t245 = -t457 * t586 + t458 * t489 + t459 * t490;
t244 = -t457 * t585 + t458 * t487 + t459 * t488;
t243 = (t413 * t566 + t414 * t568) * t567 + t662;
t242 = -t468 * t551 - t506 * t543 + (-t414 * t574 + t456 * t659) * t567;
t241 = -t468 * t590 + t506 * t541 + (t413 * t574 - t455 * t659) * t567;
t231 = (-t520 + t674) * t569 + t568 * t601;
t230 = t401 * t569 + t566 * t601 + t664;
t225 = -t466 * t553 + t467 * t554 - t504 * t529 + t505 * t530 + (-t465 * t574 + t503 * t659) * t567;
t224 = t400 * t551 - t590 * t673 + t460;
t223 = t394 * t553 + t396 * t531 - t398 * t591;
t222 = t393 * t553 + t395 * t531 - t397 * t591;
t221 = t413 * t551 + t414 * t590 + t455 * t543 - t456 * t541;
t220 = t586 * t668 + t680;
t219 = -t553 * t676 - t671;
t218 = (t400 * t566 + t401 * t568) * t567 + t604;
t213 = t465 * t551 + t466 * t586 + t467 * t528 - t485 * t504 + t486 * t505 + t503 * t543;
t212 = -t465 * t590 + t466 * t585 + t467 * t526 - t483 * t504 + t484 * t505 + t503 * t541;
t211 = -t394 * t586 + t396 * t489 + t398 * t490;
t210 = -t393 * t586 + t395 * t489 + t397 * t490;
t209 = -t394 * t585 + t396 * t487 + t398 * t488;
t208 = -t393 * t585 + t395 * t487 + t397 * t488;
t207 = -t586 * t667 + t357 + t359;
t206 = t553 * t678 - t427 - t428;
t197 = t551 * t638 + t641 * t699;
t196 = t380 * t699 - t442 * t590 + t611;
t195 = t585 * t675 - t682;
t194 = (-t520 + t643) * t569 + t568 * t597;
t193 = t382 * t569 + t566 * t597 + t645;
t192 = -t585 * t677 - t342 - t345;
t191 = -t409 * t553 + t411 * t554 - t452 * t529 + t454 * t530 + (-t407 * t574 + t450 * t659) * t567;
t190 = -t408 * t553 + t410 * t554 - t451 * t529 + t453 * t530 + (-t406 * t574 + t449 * t659) * t567;
t189 = (-t282 + t669) * t569 + t568 * t602;
t188 = t283 * t569 + t566 * t602 + t670;
t187 = -t586 * t639 + t357 + t680;
t186 = t553 * t644 - t428 - t671;
t185 = t380 * t551 - t590 * t641 + t681;
t184 = (t380 * t566 + t382 * t568) * t567 + t592;
t183 = t551 * t605 + t609 * t699;
t182 = -t590 * t668 + t676 * t699 + t611;
t181 = t407 * t551 + t409 * t586 + t411 * t528 + t450 * t543 - t452 * t485 + t454 * t486;
t180 = t406 * t551 + t408 * t586 + t410 * t528 + t449 * t543 - t451 * t485 + t453 * t486;
t179 = -t407 * t590 + t409 * t585 + t411 * t526 + t450 * t541 - t452 * t483 + t454 * t484;
t178 = -t406 * t590 + t408 * t585 + t410 * t526 + t449 * t541 - t451 * t483 + t453 * t484;
t177 = (-t520 + t610) * t569 + t568 * t588;
t176 = t566 * t588 + t569 * t675 + t645;
t175 = t283 * t553 + t341 * t586 + t401 * t529 - t461 * t485;
t174 = -t282 * t553 - t341 * t585 - t400 * t529 + t461 * t483;
t173 = t326 * t586 - t442 * t485 + t694;
t172 = -t267 * t553 - t380 * t529 + t687;
t171 = (t282 * t566 + t283 * t568) * t567 + t606;
t170 = -t585 * t642 - t342 - t682;
t169 = t338 * t553 + t339 * t531 - t340 * t591 + t457 * t529 + t458 * t463 + t459 * t464;
t168 = t551 * t676 - t590 * t609 + t681;
t167 = (t566 * t676 + t568 * t675) * t567 + t592;
t166 = t469 + t683 * t551 + t666 * t543 + (t401 * t659 + t574 * t689) * t567;
t165 = -t341 * t590 + t461 * t541 + (t282 * t574 + t659 * t674) * t567 + t640;
t160 = -t338 * t586 + t339 * t489 + t340 * t490 + t404 * t458 + t405 * t459 + t457 * t485;
t159 = -t338 * t585 + t339 * t487 + t340 * t488 + t402 * t458 + t403 * t459 + t457 * t483;
t158 = (-t267 + t648) * t569 + t568 * t598;
t157 = t269 * t569 + t566 * t598 + t649;
t156 = t297 * t569 + (-t222 * t568 + t223 * t566) * t567;
t155 = -t222 * t590 + t223 * t551 - t297 * t699;
t150 = -t222 * t585 - t223 * t586 + t297 * t553;
t147 = t282 * t551 + t400 * t543 + t541 * t673 - t590 * t689 + t672;
t134 = t245 * t569 + (-t210 * t568 + t211 * t566) * t567;
t133 = t244 * t569 + (-t208 * t568 + t209 * t566) * t567;
t132 = -t210 * t590 + t211 * t551 - t245 * t699;
t131 = -t208 * t590 + t209 * t551 - t244 * t699;
t130 = -t210 * t585 - t211 * t586 + t245 * t553;
t129 = -t208 * t585 - t209 * t586 + t244 * t553;
t112 = t277 * t553 + t279 * t531 - t281 * t591 + t394 * t529 + t396 * t463 + t398 * t464;
t111 = t276 * t553 + t278 * t531 - t280 * t591 + t393 * t529 + t395 * t463 + t397 * t464;
t102 = (t267 * t566 + t269 * t568) * t567 + t594;
t97 = t485 * t667 - t586 * t685 + t688 + t694;
t96 = (-t267 - t299) * t553 + t678 * t529 + t684 + t687;
t91 = -t485 * t668 + t586 * t686 + t618;
t90 = -t529 * t676 - t553 * t693 + t614;
t89 = -t277 * t586 + t279 * t489 + t281 * t490 + t394 * t485 + t396 * t404 + t398 * t405;
t88 = -t276 * t586 + t278 * t489 + t280 * t490 + t393 * t485 + t395 * t404 + t397 * t405;
t87 = -t277 * t585 + t279 * t487 + t281 * t488 + t394 * t483 + t396 * t402 + t398 * t403;
t86 = -t276 * t585 + t278 * t487 + t280 * t488 + t393 * t483 + t395 * t402 + t397 * t403;
t85 = (t648 - t693) * t569 + t568 * t589;
t84 = t566 * t589 + t569 * t692 + t649;
t83 = t646 * t551 + t638 * t543 + (t382 * t659 + t574 * t650) * t567 + t679;
t82 = -t326 * t590 + t442 * t541 + (t267 * t574 + t643 * t659) * t567 + t593;
t81 = t225 * t569 + (-t190 * t568 + t191 * t566) * t567;
t72 = t483 * t677 - t585 * t691 + t690 + t695;
t71 = t213 * t569 + (-t180 * t568 + t181 * t566) * t567;
t70 = t212 * t569 + (-t178 * t568 + t179 * t566) * t567;
t69 = -t483 * t675 + t585 * t692 + t619;
t68 = (t566 * t693 + t568 * t692) * t567 + t594;
t67 = t267 * t551 + t380 * t543 + t541 * t641 - t590 * t650 + t615;
t66 = t485 * t639 - t586 * t647 + t618 + t688;
t65 = (-t299 - t693) * t553 + t644 * t529 + t614 + t684;
t64 = t613 * t551 + t605 * t543 + (t574 * t617 + t659 * t675) * t567 + t679;
t63 = -t686 * t590 + t668 * t541 + (t574 * t693 + t610 * t659) * t567 + t593;
t62 = -t190 * t590 + t191 * t551 + t303 * t541 + t304 * t543 + (-t225 * t574 + t354 * t659) * t567;
t61 = -t180 * t590 + t181 * t551 + t290 * t541 + t706 + (-t213 * t574 + t315 * t659) * t567;
t60 = -t178 * t590 + t179 * t551 + t708 + t289 * t543 + (-t212 * t574 + t314 * t659) * t567;
t59 = t483 * t642 - t585 * t651 + t619 + t690;
t58 = t541 * t609 + t543 * t676 + t551 * t693 - t590 * t617 + t615;
t57 = t169 * t569 + (-t111 * t568 + t112 * t566) * t567;
t52 = t160 * t569 + (t566 * t89 - t568 * t88) * t567;
t51 = t159 * t569 + (t566 * t87 - t568 * t86) * t567;
t42 = -t111 * t590 + t112 * t551 + t222 * t541 + t223 * t543 + (-t169 * t574 + t297 * t659) * t567;
t41 = -t111 * t585 - t112 * t586 + t169 * t553 + t222 * t483 + t223 * t485 + t297 * t529;
t28 = t210 * t541 + t211 * t543 - t590 * t88 + t551 * t89 + (-t160 * t574 + t245 * t659) * t567;
t27 = t208 * t541 + t209 * t543 - t590 * t86 + t551 * t87 + (-t159 * t574 + t244 * t659) * t567;
t26 = t160 * t553 + t210 * t483 + t211 * t485 + t245 * t529 - t585 * t88 - t586 * t89;
t25 = t159 * t553 + t208 * t483 + t209 * t485 + t244 * t529 - t585 * t86 - t586 * t87;
t1 = [0; m(4) * t663 / 0.2e1 + t607 * t729 + (m(3) * t516 - t734) * t701 + (m(3) * t515 + t579) * t702 + t735 * (0.2e1 * t294 + 0.2e1 * t295 + t607); t50 * t702 - t48 * t701 + t49 * t702 - t47 * t701 + t52 * t702 - t51 * t701 + t71 * t702 - t70 * t701 + ((-t498 * t543 - t500 * t542 + t510 * t702 - t512 * t551 + t514 * t552) * t702 - (-t497 * t543 - t499 * t542 + t509 * t702 - t511 * t551 + t513 * t552) * t701 + (-t537 * t543 - t538 * t542 + t544 * t702 - t545 * t551 + t546 * t552) * t569) * t702 - ((-t498 * t541 + t500 * t540 - t510 * t701 + t512 * t590 + t514 * t550) * t702 - (-t497 * t541 + t499 * t540 - t509 * t701 + t511 * t590 + t513 * t550) * t701 + (-t537 * t541 + t538 * t540 - t544 * t701 + t545 * t590 + t546 * t550) * t569) * t701 + (t167 * t68 + t176 * t84 + t177 * t85) * t622 + (t102 * t184 + t157 * t193 + t158 * t194) * t624 + t569 * t55 + t569 * t56 + (t171 * t218 + t188 * t230 + t189 * t231) * t626 + t569 * t57 + t569 * t81 + (t243 * t309 + t301 * t334 + t302 * t335) * t628 + t569 * (t569 ^ 2 * t544 + (((t512 * t574 + t514 * t572) * t566 - (t511 * t574 + t513 * t572) * t568 + ((-t498 * t572 + t500 * t574) * t566 - (-t497 * t572 + t499 * t574) * t568) * qJD(2)) * t567 + (-t509 * t568 + t510 * t566 + t545 * t574 + t546 * t572 + (-t537 * t572 + t538 * t574) * qJD(2)) * t569) * t567) + 0.2e1 * m(3) * ((-t501 * t569 - t539 * t701) * (-t515 * t569 - t547 * t701) + (t502 * t569 - t539 * t702) * (t516 * t569 - t547 * t702) + (t501 * t566 + t502 * t568) * t567 ^ 2 * (t515 * t566 + t516 * t568)); t608 * t729 + t579 * t551 - t734 * t590 + (m(4) * t455 + m(5) * t400 + m(6) * t380 + m(7) * t676) * t543 + (-m(4) * t456 - m(5) * t401 - m(6) * t382 - m(7) * t675) * t541 + t735 * (t300 * t730 + t366 * t731 + 0.2e1 * t286 + 0.2e1 * t355 + t608); m(4) * (t221 * t309 + t241 * t335 + t242 * t334 + t243 * t336 + t301 * t384 + t302 * t383) + (t147 * t218 + t165 * t231 + t166 * t230 + t171 * t224 + t188 * t250 + t189 * t249) * m(5) + (t102 * t185 + t157 * t197 + t158 * t196 + t184 * t67 + t193 * t83 + t194 * t82) * m(6) + (t167 * t58 + t168 * t68 + t176 * t64 + t177 * t63 + t182 * t85 + t183 * t84) * m(7) + (t134 / 0.2e1 + t128 / 0.2e1 + t127 / 0.2e1) * t543 + (t133 / 0.2e1 + t125 / 0.2e1 + t126 / 0.2e1) * t541 + (t71 / 0.2e1 + t52 / 0.2e1 + t50 / 0.2e1 + t49 / 0.2e1) * t551 - (t70 / 0.2e1 + t51 / 0.2e1 + t48 / 0.2e1 + t47 / 0.2e1) * t590 + (t705 / 0.2e1 + t707 / 0.2e1 + t62 / 0.2e1 + t42 / 0.2e1 + t40 / 0.2e1 + t39 / 0.2e1) * t569 + ((-t543 * t290 / 0.2e1 - t708 / 0.2e1 - t60 / 0.2e1 - t27 / 0.2e1 - t21 / 0.2e1 - t22 / 0.2e1) * t568 + (t706 / 0.2e1 + t289 * t721 + t61 / 0.2e1 + t28 / 0.2e1 + t23 / 0.2e1 + t24 / 0.2e1) * t566 + (-t57 / 0.2e1 - t55 / 0.2e1 - t56 / 0.2e1 - t81 / 0.2e1) * t574 + (t148 / 0.2e1 + t149 / 0.2e1 + t156 / 0.2e1 + t354 * t716 + (-t303 * t568 + t304 * t566) * t567 / 0.2e1) * t659) * t567; (t221 * t336 + t241 * t383 + t242 * t384) * t628 + (t147 * t224 + t165 * t249 + t166 * t250) * t626 + (t185 * t67 + t196 * t82 + t197 * t83) * t624 + (t168 * t58 + t182 * t63 + t183 * t64) * t622 + (t61 + t28 + t745) * t551 - (t60 + t27 + t746) * t590 + (-t290 * t590 + t291 * t551 + t132 + t739) * t543 + (-t288 * t590 + t289 * t551 + t131 + t740) * t541 + ((-t303 * t590 + t304 * t551 + t155 + t737) * t659 + (-t354 * t635 - t42 - t62 - t705 - t707 - t743) * t574) * t567; t164 * m(5) + (t612 + t616) * t728 + (t583 + t616) * t727; t577 + t41 * t716 + t57 * t717 + t156 * t722 + t52 * t723 + t51 * t724 + (t102 * t192 + t157 * t207 + t158 * t206 + t184 * t72 + t193 * t97 + t194 * t96) * m(6) + t134 * t725 + t133 * t726 + (t164 * t218 + t171 * t272 + t174 * t231 + t175 * t230 + t188 * t308 + t189 * t307) * m(5) + (t167 * t59 + t170 * t68 + t176 * t66 + t177 * t65 + t186 * t85 + t187 * t84) * m(7) + (-t568 * t25 / 0.2e1 + t566 * t26 / 0.2e1) * t567; t42 * t717 + t26 * t718 + t25 * t719 + t130 * t720 + t129 * t721 + t155 * t722 + t28 * t723 + t27 * t724 + (t147 * t272 + t164 * t224 + t165 * t307 + t166 * t308 + t174 * t249 + t175 * t250) * m(5) + t132 * t725 + t131 * t726 + t576 + (t168 * t59 + t170 * t58 + t182 * t65 + t183 * t66 + t186 * t63 + t187 * t64) * m(7) + (t185 * t72 + t192 * t67 + t196 * t96 + t197 * t97 + t206 * t82 + t207 * t83) * m(6) + (-t574 * t41 / 0.2e1 + t150 * t629) * t567; t553 * t41 + t529 * t150 - t586 * t26 - t585 * t25 + (t170 * t59 + t186 * t65 + t187 * t66) * t622 + (t192 * t72 + t206 * t96 + t207 * t97) * t624 + t485 * t130 + t483 * t129 + t580 + (t164 * t272 + t174 * t307 + t175 * t308) * t626; t583 * t727 + t612 * t728; t577 + (t102 * t251 + t157 * t306 + t158 * t305 + t161 * t184 + t172 * t194 + t173 * t193) * m(6) + (t167 * t69 + t176 * t91 + t177 * t90 + t195 * t68 + t219 * t85 + t220 * t84) * m(7); t576 + (t161 * t185 + t172 * t196 + t173 * t197 + t251 * t67 + t305 * t82 + t306 * t83) * m(6) + (t168 * t69 + t182 * t90 + t183 * t91 + t195 * t58 + t219 * t63 + t220 * t64) * m(7); (t161 * t192 + t172 * t206 + t173 * t207 + t251 * t72 + t305 * t96 + t306 * t97) * m(6) + (t170 * t69 + t186 * t90 + t187 * t91 + t195 * t59 + t219 * t65 + t220 * t66) * m(7) + t580; (t195 * t69 + t219 * t90 + t220 * t91) * t622 + t580 + (t161 * t251 + t172 * t305 + t173 * t306) * t624; t444 * m(7); (t167 * t444 + t176 * t386 + t177 * t388 + t478 * t84 + t480 * t85 + t523 * t68) * m(7); (t168 * t444 + t182 * t388 + t183 * t386 + t478 * t64 + t480 * t63 + t523 * t58) * m(7); (t170 * t444 + t186 * t388 + t187 * t386 + t478 * t66 + t480 * t65 + t523 * t59) * m(7); (t195 * t444 + t219 * t388 + t220 * t386 + t478 * t91 + t480 * t90 + t523 * t69) * m(7); (t386 * t478 + t388 * t480 + t444 * t523) * t622;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
