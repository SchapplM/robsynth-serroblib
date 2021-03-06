% Calculate time derivative of joint inertia matrix for
% S6PRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 00:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_inertiaDJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR1_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRR1_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRR1_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:37:36
% EndTime: 2019-03-09 00:38:23
% DurationCPUTime: 27.86s
% Computational Cost: add. (189453->1323), mult. (246299->1820), div. (0->0), fcn. (275861->14), ass. (0->583)
t564 = sin(pkin(12));
t566 = cos(pkin(12));
t573 = cos(qJ(2));
t567 = cos(pkin(6));
t570 = sin(qJ(2));
t712 = t567 * t570;
t541 = t564 * t573 + t566 * t712;
t563 = qJ(3) + qJ(4);
t663 = qJ(5) + t563;
t613 = cos(t663);
t555 = sin(t663);
t565 = sin(pkin(6));
t718 = t565 * t566;
t664 = t555 * t718;
t508 = t541 * t613 - t664;
t568 = sin(qJ(6));
t571 = cos(qJ(6));
t711 = t567 * t573;
t594 = -t564 * t570 + t566 * t711;
t455 = -t508 * t568 - t571 * t594;
t456 = t508 * t571 - t568 * t594;
t604 = t565 * t613;
t581 = -t541 * t555 - t566 * t604;
t297 = Icges(7,5) * t456 + Icges(7,6) * t455 - Icges(7,3) * t581;
t299 = Icges(7,4) * t456 + Icges(7,2) * t455 - Icges(7,6) * t581;
t301 = Icges(7,1) * t456 + Icges(7,4) * t455 - Icges(7,5) * t581;
t161 = -t297 * t581 + t299 * t455 + t301 * t456;
t666 = t564 * t712;
t543 = t566 * t573 - t666;
t719 = t564 * t565;
t510 = t543 * t613 + t555 * t719;
t542 = t564 * t711 + t566 * t570;
t457 = -t510 * t568 + t542 * t571;
t458 = t510 * t571 + t542 * t568;
t582 = -t543 * t555 + t564 * t604;
t298 = Icges(7,5) * t458 + Icges(7,6) * t457 - Icges(7,3) * t582;
t300 = Icges(7,4) * t458 + Icges(7,2) * t457 - Icges(7,6) * t582;
t302 = Icges(7,1) * t458 + Icges(7,4) * t457 - Icges(7,5) * t582;
t162 = -t298 * t581 + t300 * t455 + t302 * t456;
t590 = t570 * t604;
t524 = t567 * t555 + t590;
t714 = t565 * t573;
t511 = -t524 * t568 - t571 * t714;
t716 = t565 * t570;
t523 = t555 * t716 - t567 * t613;
t597 = -t524 * t571 + t568 * t714;
t378 = -Icges(7,5) * t597 + Icges(7,6) * t511 + Icges(7,3) * t523;
t379 = -Icges(7,4) * t597 + Icges(7,2) * t511 + Icges(7,6) * t523;
t380 = -Icges(7,1) * t597 + Icges(7,4) * t511 + Icges(7,5) * t523;
t185 = -t378 * t581 + t379 * t455 + t380 * t456;
t531 = t541 * qJD(2);
t674 = qJD(2) * t573;
t533 = -qJD(2) * t666 + t566 * t674;
t530 = t594 * qJD(2);
t561 = qJD(3) + qJD(4);
t557 = qJD(5) + t561;
t440 = t530 * t613 + t557 * t581;
t306 = -qJD(6) * t456 - t440 * t568 + t531 * t571;
t307 = qJD(6) * t455 + t440 * t571 + t531 * t568;
t605 = t557 * t613;
t439 = t530 * t555 + t541 * t605 - t557 * t664;
t197 = Icges(7,5) * t307 + Icges(7,6) * t306 + Icges(7,3) * t439;
t199 = Icges(7,4) * t307 + Icges(7,2) * t306 + Icges(7,6) * t439;
t201 = Icges(7,1) * t307 + Icges(7,4) * t306 + Icges(7,5) * t439;
t61 = -t197 * t581 + t199 * t455 + t201 * t456 + t297 * t439 + t299 * t306 + t301 * t307;
t532 = t542 * qJD(2);
t442 = -t532 * t613 + t557 * t582;
t308 = -qJD(6) * t458 - t442 * t568 + t533 * t571;
t309 = qJD(6) * t457 + t442 * t571 + t533 * t568;
t441 = t543 * t605 + (t557 * t719 - t532) * t555;
t198 = Icges(7,5) * t309 + Icges(7,6) * t308 + Icges(7,3) * t441;
t200 = Icges(7,4) * t309 + Icges(7,2) * t308 + Icges(7,6) * t441;
t202 = Icges(7,1) * t309 + Icges(7,4) * t308 + Icges(7,5) * t441;
t62 = -t198 * t581 + t200 * t455 + t202 * t456 + t298 * t439 + t300 * t306 + t302 * t307;
t675 = qJD(2) * t570;
t478 = -t523 * t557 + t604 * t674;
t643 = t565 * t675;
t382 = qJD(6) * t597 - t478 * t568 + t571 * t643;
t383 = qJD(6) * t511 + t478 * t571 + t568 * t643;
t642 = t565 * t674;
t477 = t557 * t590 + (t557 * t567 + t642) * t555;
t242 = Icges(7,5) * t383 + Icges(7,6) * t382 + Icges(7,3) * t477;
t243 = Icges(7,4) * t383 + Icges(7,2) * t382 + Icges(7,6) * t477;
t244 = Icges(7,1) * t383 + Icges(7,4) * t382 + Icges(7,5) * t477;
t92 = -t242 * t581 + t243 * t455 + t244 * t456 + t306 * t379 + t307 * t380 + t378 * t439;
t11 = t161 * t531 + t162 * t533 - t594 * t61 + t542 * t62 + (t185 * t675 - t573 * t92) * t565;
t283 = Icges(6,5) * t440 - Icges(6,6) * t439 + Icges(6,3) * t531;
t285 = Icges(6,4) * t440 - Icges(6,2) * t439 + Icges(6,6) * t531;
t287 = Icges(6,1) * t440 - Icges(6,4) * t439 + Icges(6,5) * t531;
t394 = Icges(6,5) * t508 + Icges(6,6) * t581 - Icges(6,3) * t594;
t396 = Icges(6,4) * t508 + Icges(6,2) * t581 - Icges(6,6) * t594;
t398 = Icges(6,1) * t508 + Icges(6,4) * t581 - Icges(6,5) * t594;
t113 = -t283 * t594 + t285 * t581 + t287 * t508 + t394 * t531 - t396 * t439 + t398 * t440;
t284 = Icges(6,5) * t442 - Icges(6,6) * t441 + Icges(6,3) * t533;
t286 = Icges(6,4) * t442 - Icges(6,2) * t441 + Icges(6,6) * t533;
t288 = Icges(6,1) * t442 - Icges(6,4) * t441 + Icges(6,5) * t533;
t395 = Icges(6,5) * t510 + Icges(6,6) * t582 + Icges(6,3) * t542;
t397 = Icges(6,4) * t510 + Icges(6,2) * t582 + Icges(6,6) * t542;
t399 = Icges(6,1) * t510 + Icges(6,4) * t582 + Icges(6,5) * t542;
t114 = -t284 * t594 + t286 * t581 + t288 * t508 + t395 * t531 - t397 * t439 + t399 * t440;
t372 = Icges(6,5) * t478 - Icges(6,6) * t477 + Icges(6,3) * t643;
t373 = Icges(6,4) * t478 - Icges(6,2) * t477 + Icges(6,6) * t643;
t374 = Icges(6,1) * t478 - Icges(6,4) * t477 + Icges(6,5) * t643;
t459 = Icges(6,5) * t524 - Icges(6,6) * t523 - Icges(6,3) * t714;
t460 = Icges(6,4) * t524 - Icges(6,2) * t523 - Icges(6,6) * t714;
t461 = Icges(6,1) * t524 - Icges(6,4) * t523 - Icges(6,5) * t714;
t152 = -t372 * t594 + t373 * t581 + t374 * t508 - t439 * t460 + t440 * t461 + t459 * t531;
t212 = -t394 * t594 + t396 * t581 + t398 * t508;
t213 = -t395 * t594 + t397 * t581 + t399 * t508;
t240 = -t459 * t594 + t460 * t581 + t461 * t508;
t758 = t11 - t113 * t594 + t114 * t542 + t212 * t531 + t213 * t533 + (-t152 * t573 + t240 * t675) * t565;
t115 = t283 * t542 + t285 * t582 + t287 * t510 + t394 * t533 - t396 * t441 + t398 * t442;
t116 = t284 * t542 + t286 * t582 + t288 * t510 + t395 * t533 - t397 * t441 + t399 * t442;
t163 = -t297 * t582 + t299 * t457 + t301 * t458;
t164 = -t298 * t582 + t300 * t457 + t302 * t458;
t186 = -t378 * t582 + t379 * t457 + t380 * t458;
t63 = -t197 * t582 + t199 * t457 + t201 * t458 + t297 * t441 + t299 * t308 + t301 * t309;
t64 = -t198 * t582 + t200 * t457 + t202 * t458 + t298 * t441 + t300 * t308 + t302 * t309;
t93 = -t242 * t582 + t243 * t457 + t244 * t458 + t308 * t379 + t309 * t380 + t378 * t441;
t12 = t163 * t531 + t164 * t533 - t594 * t63 + t542 * t64 + (t186 * t675 - t573 * t93) * t565;
t153 = t372 * t542 + t373 * t582 + t374 * t510 - t441 * t460 + t442 * t461 + t459 * t533;
t214 = t394 * t542 + t396 * t582 + t398 * t510;
t215 = t395 * t542 + t397 * t582 + t399 * t510;
t241 = t459 * t542 + t460 * t582 + t461 * t510;
t757 = t12 - t115 * t594 + t116 * t542 + t214 * t531 + t215 * t533 + (-t153 * t573 + t241 * t675) * t565;
t130 = -t285 * t523 + t287 * t524 - t396 * t477 + t398 * t478 + (-t283 * t573 + t394 * t675) * t565;
t131 = -t286 * t523 + t288 * t524 - t397 * t477 + t399 * t478 + (-t284 * t573 + t395 * t675) * t565;
t165 = -t373 * t523 + t374 * t524 - t460 * t477 + t461 * t478 + (-t372 * t573 + t459 * t675) * t565;
t100 = t242 * t523 + t243 * t511 - t244 * t597 + t378 * t477 + t379 * t382 + t380 * t383;
t167 = t297 * t523 + t299 * t511 - t301 * t597;
t168 = t298 * t523 + t300 * t511 - t302 * t597;
t206 = t378 * t523 + t379 * t511 - t380 * t597;
t65 = t197 * t523 + t199 * t511 - t201 * t597 + t297 * t477 + t299 * t382 + t301 * t383;
t66 = t198 * t523 + t200 * t511 - t202 * t597 + t298 * t477 + t300 * t382 + t302 * t383;
t17 = t167 * t531 + t168 * t533 - t594 * t65 + t542 * t66 + (-t100 * t573 + t206 * t675) * t565;
t222 = -t394 * t714 - t396 * t523 + t398 * t524;
t223 = -t395 * t714 - t397 * t523 + t399 * t524;
t249 = -t459 * t714 - t460 * t523 + t461 * t524;
t756 = t17 - t130 * t594 + t131 * t542 + t222 * t531 + t223 * t533 + (-t165 * t573 + t249 * t675) * t565;
t203 = rSges(7,1) * t307 + rSges(7,2) * t306 + rSges(7,3) * t439;
t321 = pkin(5) * t440 + pkin(11) * t439;
t710 = t203 + t321;
t204 = rSges(7,1) * t309 + rSges(7,2) * t308 + rSges(7,3) * t441;
t709 = pkin(5) * t442 + pkin(11) * t441 + t204;
t245 = rSges(7,1) * t383 + rSges(7,2) * t382 + rSges(7,3) * t477;
t708 = -pkin(5) * t478 - pkin(11) * t477 - t245;
t303 = rSges(7,1) * t456 + rSges(7,2) * t455 - rSges(7,3) * t581;
t432 = pkin(5) * t508 - pkin(11) * t581;
t701 = t303 + t432;
t304 = rSges(7,1) * t458 + rSges(7,2) * t457 - rSges(7,3) * t582;
t700 = pkin(5) * t510 - pkin(11) * t582 + t304;
t381 = -rSges(7,1) * t597 + rSges(7,2) * t511 + rSges(7,3) * t523;
t755 = pkin(5) * t524 + pkin(11) * t523 + t381;
t558 = sin(t563);
t608 = t561 * t719 - t532;
t559 = cos(t563);
t720 = t559 * t561;
t450 = -t543 * t720 - t558 * t608;
t721 = t558 * t561;
t451 = -t543 * t721 + t559 * t608;
t317 = rSges(5,1) * t451 + rSges(5,2) * t450 + rSges(5,3) * t533;
t515 = -t543 * t558 + t559 * t719;
t516 = t543 * t559 + t558 * t719;
t416 = rSges(5,1) * t516 + rSges(5,2) * t515 + rSges(5,3) * t542;
t609 = t561 * t718 - t530;
t448 = -t541 * t720 + t558 * t609;
t449 = -t541 * t721 - t559 * t609;
t316 = rSges(5,1) * t449 + rSges(5,2) * t448 + rSges(5,3) * t531;
t513 = -t541 * t558 - t559 * t718;
t514 = t541 * t559 - t558 * t718;
t415 = rSges(5,1) * t514 + rSges(5,2) * t513 - rSges(5,3) * t594;
t703 = t542 * t316 + t533 * t415;
t172 = t317 * t594 - t531 * t416 + t703;
t290 = rSges(6,1) * t442 - rSges(6,2) * t441 + rSges(6,3) * t533;
t401 = rSges(6,1) * t510 + rSges(6,2) * t582 + rSges(6,3) * t542;
t289 = rSges(6,1) * t440 - rSges(6,2) * t439 + rSges(6,3) * t531;
t400 = rSges(6,1) * t508 + rSges(6,2) * t581 - rSges(6,3) * t594;
t705 = t542 * t289 + t533 * t400;
t166 = t290 * t594 - t531 * t401 + t705;
t500 = t530 * pkin(2) + t531 * pkin(8);
t501 = -t532 * pkin(2) + t533 * pkin(8);
t680 = t500 * t719 + t501 * t718;
t752 = -0.2e1 * t531;
t751 = 0.2e1 * t594;
t750 = m(5) / 0.2e1;
t749 = m(6) / 0.2e1;
t748 = m(7) / 0.2e1;
t747 = t439 / 0.2e1;
t746 = t441 / 0.2e1;
t745 = t477 / 0.2e1;
t744 = -t581 / 0.2e1;
t743 = -t582 / 0.2e1;
t742 = t523 / 0.2e1;
t741 = t531 / 0.2e1;
t740 = t533 / 0.2e1;
t739 = -t594 / 0.2e1;
t738 = t542 / 0.2e1;
t737 = t564 / 0.2e1;
t736 = -t566 / 0.2e1;
t735 = t567 / 0.2e1;
t572 = cos(qJ(3));
t733 = t572 * pkin(3);
t730 = pkin(3) * qJD(3);
t729 = Icges(3,4) * t570;
t728 = Icges(3,4) * t573;
t569 = sin(qJ(3));
t544 = t567 * t572 - t569 * t716;
t713 = t567 * t569;
t715 = t565 * t572;
t545 = t570 * t715 + t713;
t485 = Icges(4,5) * t545 + Icges(4,6) * t544 - Icges(4,3) * t714;
t486 = Icges(4,4) * t545 + Icges(4,2) * t544 - Icges(4,6) * t714;
t487 = Icges(4,1) * t545 + Icges(4,4) * t544 - Icges(4,5) * t714;
t717 = t565 * t569;
t665 = t566 * t717;
t595 = -t541 * t572 + t665;
t596 = t541 * t569 + t566 * t715;
t251 = -t485 * t594 - t486 * t596 - t487 * t595;
t727 = t531 * t251;
t519 = -t543 * t569 + t564 * t715;
t667 = t564 * t717;
t520 = t543 * t572 + t667;
t252 = t485 * t542 + t486 * t519 + t487 * t520;
t724 = t533 * t252;
t670 = t569 * t730;
t546 = -pkin(4) * t721 - t670;
t669 = t572 * t730;
t547 = pkin(4) * t720 + t669;
t678 = pkin(4) * t559;
t266 = pkin(10) * t531 + t530 * t678 + t541 * t546 - t547 * t718 + t596 * t730;
t255 = t542 * t266;
t638 = pkin(4) * t558;
t353 = -pkin(10) * t594 + t541 * t678 - t638 * t718;
t328 = t533 * t353;
t707 = t255 + t328;
t267 = pkin(10) * t533 - t519 * t730 - t532 * t678 + t543 * t546 + t547 * t719;
t706 = -t267 - t290;
t704 = t701 * t542;
t702 = t700 * t643;
t588 = t519 * qJD(3);
t351 = pkin(3) * t588 + pkin(9) * t533 - t532 * t733;
t699 = -t317 - t351;
t589 = t596 * qJD(3);
t350 = -pkin(3) * t589 + pkin(9) * t531 + t530 * t733;
t327 = t542 * t350;
t417 = -pkin(3) * t665 - pkin(9) * t594 + t541 * t733;
t368 = t533 * t417;
t698 = t327 + t368;
t329 = t542 * t353;
t362 = t542 * t400;
t697 = t329 + t362;
t354 = pkin(10) * t542 + t543 * t678 + t638 * t719;
t344 = t354 * t643;
t376 = t401 * t643;
t696 = t344 + t376;
t476 = t567 * t501;
t695 = t567 * t351 + t476;
t438 = t638 * t567 + (-pkin(10) * t573 + t570 * t678) * t565;
t694 = t353 * t714 - t438 * t594;
t693 = -t350 - t500;
t692 = -t353 - t400;
t691 = -t354 - t401;
t356 = (t547 - t669) * t567 + ((t546 + t670) * t570 + (pkin(10) * t570 + t573 * t678) * qJD(2)) * t565;
t375 = rSges(6,1) * t478 - rSges(6,2) * t477 + rSges(6,3) * t643;
t690 = -t356 - t375;
t489 = pkin(3) * t713 + (-pkin(9) * t573 + t570 * t733) * t565;
t688 = t417 * t714 - t489 * t594;
t418 = pkin(3) * t667 + pkin(9) * t542 + t543 * t733;
t506 = t543 * pkin(2) + t542 * pkin(8);
t502 = t567 * t506;
t687 = t567 * t418 + t502;
t593 = t561 * t567 + t642;
t668 = t561 * t716;
t503 = -t558 * t593 - t559 * t668;
t504 = -t558 * t668 + t559 * t593;
t407 = rSges(5,1) * t504 + rSges(5,2) * t503 + rSges(5,3) * t643;
t453 = t567 * t669 + (-t570 * t670 + (pkin(9) * t570 + t573 * t733) * qJD(2)) * t565;
t686 = -t407 - t453;
t685 = -t415 - t417;
t684 = -t416 - t418;
t462 = rSges(6,1) * t524 - rSges(6,2) * t523 - rSges(6,3) * t714;
t275 = t400 * t714 - t462 * t594;
t683 = -t438 - t462;
t528 = -t558 * t716 + t559 * t567;
t529 = t558 * t567 + t559 * t716;
t470 = rSges(5,1) * t529 + rSges(5,2) * t528 - rSges(5,3) * t714;
t292 = t415 * t714 - t470 * t594;
t682 = -t470 - t489;
t681 = 0.2e1 * t680;
t505 = t541 * pkin(2) - pkin(8) * t594;
t679 = t505 * t719 + t506 * t718;
t676 = qJD(2) * t565;
t310 = Icges(5,5) * t449 + Icges(5,6) * t448 + Icges(5,3) * t531;
t312 = Icges(5,4) * t449 + Icges(5,2) * t448 + Icges(5,6) * t531;
t314 = Icges(5,1) * t449 + Icges(5,4) * t448 + Icges(5,5) * t531;
t409 = Icges(5,5) * t514 + Icges(5,6) * t513 - Icges(5,3) * t594;
t411 = Icges(5,4) * t514 + Icges(5,2) * t513 - Icges(5,6) * t594;
t413 = Icges(5,1) * t514 + Icges(5,4) * t513 - Icges(5,5) * t594;
t134 = t312 * t528 + t314 * t529 + t411 * t503 + t413 * t504 + (-t310 * t573 + t409 * t675) * t565;
t311 = Icges(5,5) * t451 + Icges(5,6) * t450 + Icges(5,3) * t533;
t313 = Icges(5,4) * t451 + Icges(5,2) * t450 + Icges(5,6) * t533;
t315 = Icges(5,1) * t451 + Icges(5,4) * t450 + Icges(5,5) * t533;
t410 = Icges(5,5) * t516 + Icges(5,6) * t515 + Icges(5,3) * t542;
t412 = Icges(5,4) * t516 + Icges(5,2) * t515 + Icges(5,6) * t542;
t414 = Icges(5,1) * t516 + Icges(5,4) * t515 + Icges(5,5) * t542;
t135 = t313 * t528 + t315 * t529 + t412 * t503 + t414 * t504 + (-t311 * t573 + t410 * t675) * t565;
t404 = Icges(5,5) * t504 + Icges(5,6) * t503 + Icges(5,3) * t643;
t405 = Icges(5,4) * t504 + Icges(5,2) * t503 + Icges(5,6) * t643;
t406 = Icges(5,1) * t504 + Icges(5,4) * t503 + Icges(5,5) * t643;
t467 = Icges(5,5) * t529 + Icges(5,6) * t528 - Icges(5,3) * t714;
t468 = Icges(5,4) * t529 + Icges(5,2) * t528 - Icges(5,6) * t714;
t469 = Icges(5,1) * t529 + Icges(5,4) * t528 - Icges(5,5) * t714;
t171 = t405 * t528 + t406 * t529 + t468 * t503 + t469 * t504 + (-t404 * t573 + t467 * t675) * t565;
t232 = -t409 * t714 + t411 * t528 + t413 * t529;
t233 = -t410 * t714 + t412 * t528 + t414 * t529;
t256 = -t467 * t714 + t468 * t528 + t469 * t529;
t39 = -t134 * t594 + t135 * t542 + t232 * t531 + t233 * t533 + (-t171 * t573 + t256 * t675) * t565;
t671 = -t39 - t756;
t662 = -t267 - t709;
t661 = -t356 + t708;
t660 = t266 * t714 - t356 * t594 + t531 * t438;
t659 = t567 * t267 + t695;
t658 = -t266 + t693;
t657 = -t351 + t706;
t656 = t329 + t704;
t655 = t289 * t714 - t375 * t594 + t531 * t462;
t654 = t344 + t702;
t653 = -t353 - t701;
t652 = -t354 - t700;
t651 = t316 * t714 - t407 * t594 + t531 * t470;
t650 = t350 * t714 - t453 * t594 + t531 * t489;
t649 = t567 * t354 + t687;
t648 = -t417 + t692;
t647 = -t418 + t691;
t646 = -t453 + t690;
t645 = -t438 - t755;
t644 = -t489 + t683;
t641 = t719 / 0.2e1;
t640 = -t718 / 0.2e1;
t639 = -t714 / 0.2e1;
t636 = 2 * m(4);
t634 = 0.2e1 * m(5);
t632 = 0.2e1 * m(6);
t630 = 0.2e1 * m(7);
t521 = -qJD(3) * t545 - t569 * t642;
t522 = qJD(3) * t544 + t572 * t642;
t443 = rSges(4,1) * t522 + rSges(4,2) * t521 + rSges(4,3) * t643;
t538 = (pkin(2) * t573 + pkin(8) * t570) * t676;
t629 = (-t443 - t538) * t565;
t488 = rSges(4,1) * t545 + rSges(4,2) * t544 - rSges(4,3) * t714;
t548 = (pkin(2) * t570 - pkin(8) * t573) * t565;
t628 = (-t488 - t548) * t565;
t195 = t542 * t203;
t273 = t533 * t303;
t296 = t542 * t321;
t391 = t533 * t432;
t627 = t195 + t273 + t296 + t391;
t626 = -t351 + t662;
t625 = -t453 + t661;
t624 = t267 * t751 + t354 * t752 + 0.2e1 * t255 + 0.2e1 * t328;
t623 = t705 + t707;
t622 = -t417 + t653;
t621 = -t418 + t652;
t620 = t351 * t751 + t418 * t752 + 0.2e1 * t327 + 0.2e1 * t368;
t341 = t350 * t719;
t342 = t351 * t718;
t619 = 0.2e1 * t341 + 0.2e1 * t342 + t681;
t618 = t341 + t342 + t680;
t617 = 0.2e1 * t166;
t616 = 0.2e1 * t172;
t615 = -t489 + t645;
t614 = t417 * t719 + t418 * t718 + t679;
t207 = t275 + t694;
t190 = -t594 * t755 + t701 * t714;
t612 = t643 / 0.2e1;
t611 = (-t538 + t686) * t565;
t610 = (-t548 + t682) * t565;
t607 = (-t538 + t646) * t565;
t606 = (-t548 + t644) * t565;
t603 = t627 + t707;
t602 = t531 * t755 + t594 * t708 + t710 * t714;
t259 = t266 * t719;
t260 = t267 * t718;
t600 = t259 + t260 + t618;
t599 = t655 + t660;
t598 = t353 * t719 + t354 * t718 + t614;
t159 = t190 + t694;
t592 = (-t538 + t625) * t565;
t591 = (-t548 + t615) * t565;
t587 = t617 + t624;
t586 = t700 * t752 + t709 * t751 + 0.2e1 * t195 + 0.2e1 * t273 + 0.2e1 * t296 + 0.2e1 * t391;
t80 = -t161 * t594 + t162 * t542 - t185 * t714;
t81 = -t163 * t594 + t164 * t542 - t186 * t714;
t91 = -t167 * t594 + t168 * t542 - t206 * t714;
t585 = (-t249 * t714 + t91) * t643 + (t223 * t643 + t757) * t542 + (-t222 * t643 - t758) * t594 + (-t214 * t594 + t215 * t542 - t241 * t714 + t81) * t533 + (-t212 * t594 + t213 * t542 - t240 * t714 + t80) * t531;
t584 = t602 + t660;
t14 = t100 * t523 + t167 * t439 + t168 * t441 + t206 * t477 - t581 * t65 - t582 * t66;
t3 = t161 * t439 + t162 * t441 + t185 * t477 + t523 * t92 - t581 * t61 - t582 * t62;
t4 = t163 * t439 + t164 * t441 + t186 * t477 + t523 * t93 - t581 * t63 - t582 * t64;
t75 = -t161 * t581 - t162 * t582 + t185 * t523;
t76 = -t163 * t581 - t164 * t582 + t186 * t523;
t88 = -t167 * t581 - t168 * t582 + t206 * t523;
t583 = t11 * t744 + t12 * t743 + t14 * t639 + t17 * t742 + t3 * t739 + t4 * t738 + t88 * t612 + t76 * t740 + t75 * t741 + t91 * t745 + t81 * t746 + t80 * t747;
t103 = -t203 * t582 + t204 * t581 + t303 * t441 - t304 * t439;
t463 = qJD(3) * t595 - t530 * t569;
t464 = t530 * t572 - t589;
t336 = rSges(4,1) * t464 + rSges(4,2) * t463 + rSges(4,3) * t531;
t465 = -qJD(3) * t520 + t532 * t569;
t466 = -t532 * t572 + t588;
t337 = rSges(4,1) * t466 + rSges(4,2) * t465 + rSges(4,3) * t533;
t428 = -rSges(4,1) * t595 - rSges(4,2) * t596 - rSges(4,3) * t594;
t429 = rSges(4,1) * t520 + rSges(4,2) * t519 + rSges(4,3) * t542;
t180 = t336 * t542 + t337 * t594 + t428 * t533 - t429 * t531;
t580 = t586 + t624;
t218 = -t409 * t594 + t411 * t513 + t413 * t514;
t219 = -t410 * t594 + t412 * t513 + t414 * t514;
t220 = t409 * t542 + t411 * t515 + t413 * t516;
t221 = t410 * t542 + t412 * t515 + t414 * t516;
t246 = -t467 * t594 + t468 * t513 + t469 * t514;
t247 = t467 * t542 + t468 * t515 + t469 * t516;
t125 = -t310 * t594 + t312 * t513 + t314 * t514 + t409 * t531 + t411 * t448 + t413 * t449;
t126 = -t311 * t594 + t313 * t513 + t315 * t514 + t410 * t531 + t412 * t448 + t414 * t449;
t157 = -t404 * t594 + t405 * t513 + t406 * t514 + t448 * t468 + t449 * t469 + t467 * t531;
t34 = -t125 * t594 + t126 * t542 + t218 * t531 + t219 * t533 + (-t157 * t573 + t246 * t675) * t565;
t127 = t310 * t542 + t312 * t515 + t314 * t516 + t409 * t533 + t411 * t450 + t413 * t451;
t128 = t311 * t542 + t313 * t515 + t315 * t516 + t410 * t533 + t412 * t450 + t414 * t451;
t158 = t404 * t542 + t405 * t515 + t406 * t516 + t450 * t468 + t451 * t469 + t467 * t533;
t35 = -t127 * t594 + t128 * t542 + t220 * t531 + t221 * t533 + (-t158 * t573 + t247 * t675) * t565;
t579 = -t594 * t34 + t542 * t35 + t531 * (-t218 * t594 + t219 * t542 - t246 * t714) + t533 * (-t220 * t594 + t221 * t542 - t247 * t714) + (-t232 * t594 + t233 * t542 - t256 * t714) * t643 + t585;
t578 = -t714 * t756 + t585;
t20 = t567 * t92 + (t564 * t62 - t566 * t61) * t565;
t21 = t567 * t93 + (t564 * t64 - t566 * t63) * t565;
t23 = t100 * t567 + (t564 * t66 - t566 * t65) * t565;
t45 = t152 * t567 + (-t113 * t566 + t114 * t564) * t565;
t46 = t153 * t567 + (-t115 * t566 + t116 * t564) * t565;
t53 = t165 * t567 + (-t130 * t566 + t131 * t564) * t565;
t84 = t185 * t567 + (-t161 * t566 + t162 * t564) * t565;
t85 = t186 * t567 + (-t163 * t566 + t164 * t564) * t565;
t95 = t206 * t567 + (-t167 * t566 + t168 * t564) * t565;
t577 = (t84 + t240 * t567 + (-t212 * t566 + t213 * t564) * t565) * t741 + (t85 + t241 * t567 + (-t214 * t566 + t215 * t564) * t565) * t740 + (t20 + t45) * t739 + (t21 + t46) * t738 + t756 * t735 + t757 * t641 + t758 * t640 + (t23 + t53) * t639 + (t95 + t249 * t567 + (-t222 * t566 + t223 * t564) * t565) * t612;
t576 = t671 * t714 + t579;
t50 = t157 * t567 + (-t125 * t566 + t126 * t564) * t565;
t51 = t158 * t567 + (-t127 * t566 + t128 * t564) * t565;
t56 = t171 * t567 + (-t134 * t566 + t135 * t564) * t565;
t575 = t34 * t640 + t35 * t641 + t39 * t735 + t50 * t739 + t51 * t738 + t56 * t639 + (t246 * t567 + (-t218 * t566 + t219 * t564) * t565) * t741 + (t247 * t567 + (-t220 * t566 + t221 * t564) * t565) * t740 + (t256 * t567 + (-t232 * t566 + t233 * t564) * t565) * t612 + t577;
t537 = (rSges(3,1) * t573 - rSges(3,2) * t570) * t676;
t536 = (Icges(3,1) * t573 - t729) * t676;
t535 = (-Icges(3,2) * t570 + t728) * t676;
t534 = (Icges(3,5) * t573 - Icges(3,6) * t570) * t676;
t527 = rSges(3,3) * t567 + (rSges(3,1) * t570 + rSges(3,2) * t573) * t565;
t526 = Icges(3,5) * t567 + (Icges(3,1) * t570 + t728) * t565;
t525 = Icges(3,6) * t567 + (Icges(3,2) * t573 + t729) * t565;
t499 = -rSges(3,1) * t532 - rSges(3,2) * t533;
t498 = rSges(3,1) * t530 - rSges(3,2) * t531;
t497 = -Icges(3,1) * t532 - Icges(3,4) * t533;
t496 = Icges(3,1) * t530 - Icges(3,4) * t531;
t495 = -Icges(3,4) * t532 - Icges(3,2) * t533;
t494 = Icges(3,4) * t530 - Icges(3,2) * t531;
t493 = -Icges(3,5) * t532 - Icges(3,6) * t533;
t492 = Icges(3,5) * t530 - Icges(3,6) * t531;
t484 = rSges(3,1) * t543 - rSges(3,2) * t542 + rSges(3,3) * t719;
t483 = rSges(3,1) * t541 + rSges(3,2) * t594 - rSges(3,3) * t718;
t482 = Icges(3,1) * t543 - Icges(3,4) * t542 + Icges(3,5) * t719;
t481 = Icges(3,1) * t541 + Icges(3,4) * t594 - Icges(3,5) * t718;
t480 = Icges(3,4) * t543 - Icges(3,2) * t542 + Icges(3,6) * t719;
t479 = Icges(3,4) * t541 + Icges(3,2) * t594 - Icges(3,6) * t718;
t437 = Icges(4,1) * t522 + Icges(4,4) * t521 + Icges(4,5) * t643;
t436 = Icges(4,4) * t522 + Icges(4,2) * t521 + Icges(4,6) * t643;
t435 = Icges(4,5) * t522 + Icges(4,6) * t521 + Icges(4,3) * t643;
t426 = Icges(4,1) * t520 + Icges(4,4) * t519 + Icges(4,5) * t542;
t425 = -Icges(4,1) * t595 - Icges(4,4) * t596 - Icges(4,5) * t594;
t424 = Icges(4,4) * t520 + Icges(4,2) * t519 + Icges(4,6) * t542;
t423 = -Icges(4,4) * t595 - Icges(4,2) * t596 - Icges(4,6) * t594;
t422 = Icges(4,5) * t520 + Icges(4,6) * t519 + Icges(4,3) * t542;
t421 = -Icges(4,5) * t595 - Icges(4,6) * t596 - Icges(4,3) * t594;
t385 = t418 * t643;
t384 = t416 * t643;
t371 = t542 * t417;
t370 = t542 * t415;
t335 = Icges(4,1) * t466 + Icges(4,4) * t465 + Icges(4,5) * t533;
t334 = Icges(4,1) * t464 + Icges(4,4) * t463 + Icges(4,5) * t531;
t333 = Icges(4,4) * t466 + Icges(4,2) * t465 + Icges(4,6) * t533;
t332 = Icges(4,4) * t464 + Icges(4,2) * t463 + Icges(4,6) * t531;
t331 = Icges(4,5) * t466 + Icges(4,6) * t465 + Icges(4,3) * t533;
t330 = Icges(4,5) * t464 + Icges(4,6) * t463 + Icges(4,3) * t531;
t320 = -t429 * t714 - t488 * t542;
t319 = t428 * t714 - t488 * t594;
t293 = -t416 * t714 - t470 * t542;
t280 = -t485 * t714 + t486 * t544 + t487 * t545;
t276 = -t401 * t714 - t462 * t542;
t264 = t428 * t542 + t429 * t594;
t263 = (-t428 - t505) * t567 + t566 * t628;
t262 = t429 * t567 + t564 * t628 + t502;
t250 = t416 * t594 + t370;
t248 = t401 * t594 + t362;
t238 = (t428 * t564 + t429 * t566) * t565 + t679;
t237 = -t422 * t714 + t424 * t544 + t426 * t545;
t236 = -t421 * t714 + t423 * t544 + t425 * t545;
t235 = (-t336 - t500) * t567 + t566 * t629;
t234 = t337 * t567 + t564 * t629 + t476;
t231 = t542 * t682 + t684 * t714;
t230 = t292 + t688;
t229 = t304 * t523 + t381 * t582;
t228 = -t303 * t523 - t381 * t581;
t227 = t422 * t542 + t424 * t519 + t426 * t520;
t226 = t421 * t542 + t423 * t519 + t425 * t520;
t225 = -t422 * t594 - t424 * t596 - t426 * t595;
t224 = -t421 * t594 - t423 * t596 - t425 * t595;
t217 = (-t505 + t685) * t567 + t566 * t610;
t216 = t416 * t567 + t564 * t610 + t687;
t211 = (t336 * t564 + t337 * t566) * t565 + t680;
t210 = -t443 * t542 - t488 * t533 + (-t337 * t573 + t429 * t675) * t565;
t209 = -t443 * t594 + t488 * t531 + (t336 * t573 - t428 * t675) * t565;
t208 = t542 * t683 + t691 * t714;
t205 = -t303 * t582 + t304 * t581;
t192 = -t594 * t684 + t370 + t371;
t191 = -t542 * t755 - t700 * t714;
t189 = -t317 * t714 - t407 * t542 - t470 * t533 + t384;
t188 = -t415 * t643 + t651;
t187 = t436 * t544 + t437 * t545 + t486 * t521 + t487 * t522 + (-t435 * t573 + t485 * t675) * t565;
t184 = (t415 * t564 + t416 * t566) * t565 + t614;
t183 = -t290 * t714 - t375 * t542 - t462 * t533 + t376;
t182 = -t400 * t643 + t655;
t181 = -t594 * t691 + t697;
t179 = (-t316 + t693) * t567 + t566 * t611;
t178 = t317 * t567 + t564 * t611 + t695;
t177 = t594 * t700 + t704;
t176 = t542 * t644 + t647 * t714;
t175 = t207 + t688;
t174 = t435 * t542 + t436 * t519 + t437 * t520 + t465 * t486 + t466 * t487 + t485 * t533;
t173 = -t435 * t594 - t436 * t596 - t437 * t595 + t463 * t486 + t464 * t487 + t485 * t531;
t170 = (-t505 + t648) * t567 + t566 * t606;
t169 = t401 * t567 + t564 * t606 + t649;
t160 = t542 * t645 + t652 * t714;
t156 = (t316 * t564 + t317 * t566) * t565 + t618;
t155 = -t594 * t647 + t371 + t697;
t154 = (t400 * t564 + t401 * t566) * t565 + t598;
t151 = t333 * t544 + t335 * t545 + t424 * t521 + t426 * t522 + (-t331 * t573 + t422 * t675) * t565;
t150 = t332 * t544 + t334 * t545 + t423 * t521 + t425 * t522 + (-t330 * t573 + t421 * t675) * t565;
t149 = t542 * t615 + t621 * t714;
t148 = t159 + t688;
t147 = t533 * t682 + t542 * t686 + t699 * t714 + t384 + t385;
t146 = t643 * t685 + t650 + t651;
t145 = (-t505 + t622) * t567 + t566 * t591;
t144 = t564 * t591 + t567 * t700 + t649;
t143 = -t594 * t652 + t656;
t142 = t331 * t542 + t333 * t519 + t335 * t520 + t422 * t533 + t424 * t465 + t426 * t466;
t141 = t330 * t542 + t332 * t519 + t334 * t520 + t421 * t533 + t423 * t465 + t425 * t466;
t140 = -t331 * t594 - t333 * t596 - t335 * t595 + t422 * t531 + t424 * t463 + t426 * t464;
t139 = -t330 * t594 - t332 * t596 - t334 * t595 + t421 * t531 + t423 * t463 + t425 * t464;
t137 = (-t289 + t658) * t567 + t566 * t607;
t136 = t290 * t567 + t564 * t607 + t659;
t120 = t533 * t683 + t542 * t690 + t706 * t714 + t696;
t119 = t643 * t692 + t599;
t118 = -t594 * t621 + t371 + t656;
t117 = (t564 * t701 + t566 * t700) * t565 + t598;
t110 = t531 * t684 - t594 * t699 + t698 + t703;
t107 = t204 * t523 + t245 * t582 + t304 * t477 - t381 * t441;
t106 = -t203 * t523 - t245 * t581 - t303 * t477 + t381 * t439;
t105 = (t289 * t564 + t290 * t566) * t565 + t600;
t104 = t531 * t691 - t594 * t706 + t623;
t102 = -t533 * t755 + t542 * t708 - t709 * t714 + t702;
t101 = -t643 * t701 + t602;
t99 = t533 * t644 + t542 * t646 + t657 * t714 + t385 + t696;
t98 = t643 * t648 + t599 + t650;
t97 = (t658 - t710) * t567 + t566 * t592;
t96 = t564 * t592 + t567 * t709 + t659;
t86 = -t531 * t700 + t594 * t709 + t627;
t79 = t531 * t647 - t594 * t657 + t623 + t698;
t70 = (t564 * t710 + t566 * t709) * t565 + t600;
t69 = t533 * t645 + t542 * t661 + t662 * t714 + t654;
t68 = t643 * t653 + t584;
t67 = t187 * t567 + (-t150 * t566 + t151 * t564) * t565;
t60 = t533 * t615 + t542 * t625 + t626 * t714 + t385 + t654;
t59 = t622 * t643 + t584 + t650;
t58 = t174 * t567 + (-t141 * t566 + t142 * t564) * t565;
t57 = t173 * t567 + (-t139 * t566 + t140 * t564) * t565;
t54 = t531 * t652 - t594 * t662 + t603;
t47 = -t150 * t594 + t151 * t542 + t236 * t531 + t237 * t533 + (-t187 * t573 + t280 * t675) * t565;
t42 = t531 * t621 - t594 * t626 + t603 + t698;
t41 = -t141 * t594 + t142 * t542 + t226 * t531 + t533 * t227 + (-t174 * t573 + t252 * t675) * t565;
t40 = -t139 * t594 + t140 * t542 + t531 * t224 + t225 * t533 + (-t173 * t573 + t251 * t675) * t565;
t1 = [0; m(4) * t681 / 0.2e1 + t619 * t750 + (m(3) * t499 + m(4) * t337 + m(5) * t317 + m(6) * t290 + m(7) * t709) * t718 + (m(3) * t498 + m(4) * t336 + m(5) * t316 + m(6) * t289 + m(7) * t710) * t719 + (t749 + t748) * (0.2e1 * t259 + 0.2e1 * t260 + t619); -t20 * t718 + t21 * t719 + t46 * t719 - t45 * t718 - t50 * t718 + t51 * t719 - t57 * t718 + t58 * t719 - ((-t480 * t531 + t482 * t530 - t493 * t718 + t495 * t594 + t497 * t541) * t719 - (-t479 * t531 + t481 * t530 - t492 * t718 + t494 * t594 + t496 * t541) * t718 + (-t525 * t531 + t526 * t530 - t534 * t718 + t535 * t594 + t536 * t541) * t567) * t718 + ((-t480 * t533 - t482 * t532 + t493 * t719 - t495 * t542 + t497 * t543) * t719 - (-t479 * t533 - t481 * t532 + t492 * t719 - t494 * t542 + t496 * t543) * t718 + (-t525 * t533 - t526 * t532 + t534 * t719 - t535 * t542 + t536 * t543) * t567) * t719 + (t117 * t70 + t144 * t96 + t145 * t97) * t630 + t567 * t23 + t567 * t53 + (t105 * t154 + t136 * t169 + t137 * t170) * t632 + t567 * t56 + (t156 * t184 + t178 * t216 + t179 * t217) * t634 + t567 * t67 + (t211 * t238 + t234 * t262 + t235 * t263) * t636 + 0.2e1 * m(3) * ((-t483 * t567 - t527 * t718) * (-t498 * t567 - t537 * t718) + (t484 * t567 - t527 * t719) * (t499 * t567 - t537 * t719) + (t483 * t564 + t484 * t566) * t565 ^ 2 * (t498 * t564 + t499 * t566)) + t567 * (t567 ^ 2 * t534 + (((t495 * t573 + t497 * t570) * t564 - (t494 * t573 + t496 * t570) * t566 + ((-t480 * t570 + t482 * t573) * t564 - (-t479 * t570 + t481 * t573) * t566) * qJD(2)) * t565 + (-t492 * t566 + t493 * t564 + t535 * t573 + t536 * t570 + (-t525 * t570 + t526 * t573) * qJD(2)) * t567) * t565); t180 * m(4) + (t616 + t620) * t750 + (t587 + t620) * t749 + (t580 + t620) * t748; t58 * t738 + t57 * t739 + ((-t226 * t566 + t227 * t564) * t740 + (-t224 * t566 + t225 * t564) * t741 - t573 * t67 / 0.2e1 + t40 * t736 + t41 * t737 + (t280 * t735 + (-t236 * t566 + t237 * t564) * t565 / 0.2e1) * t675) * t565 + (t724 / 0.2e1 + t727 / 0.2e1 + t47 / 0.2e1) * t567 + (t105 * t155 + t136 * t176 + t137 * t175 + t154 * t79 + t169 * t99 + t170 * t98) * m(6) + (t110 * t184 + t146 * t217 + t147 * t216 + t156 * t192 + t178 * t231 + t179 * t230) * m(5) + (t117 * t42 + t118 * t70 + t144 * t60 + t145 * t59 + t148 * t97 + t149 * t96) * m(7) + (t180 * t238 + t209 * t263 + t210 * t262 + t211 * t264 + t234 * t320 + t235 * t319) * m(4) + t575; t542 * t41 - t594 * t40 + ((-t236 * t594 + t237 * t542) * t675 + (-t280 * t643 - t47 + t671 - t724 - t727) * t573) * t565 + t533 * (-t226 * t594 + t227 * t542) + t531 * (-t224 * t594 + t225 * t542) + (t180 * t264 + t209 * t319 + t210 * t320) * t636 + (t110 * t192 + t146 * t230 + t147 * t231) * t634 + (t118 * t42 + t148 * t59 + t149 * t60) * t630 + (t155 * t79 + t175 * t98 + t176 * t99) * t632 + t579; t580 * t748 + t587 * t749 + t616 * t750; (t104 * t154 + t105 * t181 + t119 * t170 + t120 * t169 + t136 * t208 + t137 * t207) * m(6) + (t117 * t54 + t143 * t70 + t144 * t69 + t145 * t68 + t159 * t97 + t160 * t96) * m(7) + (t156 * t250 + t172 * t184 + t178 * t293 + t179 * t292 + t188 * t217 + t189 * t216) * m(5) + t575; (t104 * t155 + t119 * t175 + t120 * t176 + t181 * t79 + t207 * t98 + t208 * t99) * m(6) + (t110 * t250 + t146 * t292 + t147 * t293 + t172 * t192 + t188 * t230 + t189 * t231) * m(5) + (t118 * t54 + t143 * t42 + t148 * t68 + t149 * t69 + t159 * t59 + t160 * t60) * m(7) + t576; (t172 * t250 + t188 * t292 + t189 * t293) * t634 + (t104 * t181 + t119 * t207 + t120 * t208) * t632 + (t143 * t54 + t159 * t68 + t160 * t69) * t630 + t576; t586 * t748 + t617 * t749; t577 + (t105 * t248 + t136 * t276 + t137 * t275 + t154 * t166 + t169 * t183 + t170 * t182) * m(6) + (t101 * t145 + t102 * t144 + t117 * t86 + t177 * t70 + t190 * t97 + t191 * t96) * m(7); (t155 * t166 + t175 * t182 + t176 * t183 + t248 * t79 + t275 * t98 + t276 * t99) * m(6) + (t101 * t148 + t102 * t149 + t118 * t86 + t177 * t42 + t190 * t59 + t191 * t60) * m(7) + t578; (t104 * t248 + t119 * t275 + t120 * t276 + t166 * t181 + t182 * t207 + t183 * t208) * m(6) + (t101 * t159 + t102 * t160 + t143 * t86 + t177 * t54 + t190 * t68 + t191 * t69) * m(7) + t578; (t166 * t248 + t182 * t275 + t183 * t276) * t632 + (t101 * t190 + t102 * t191 + t177 * t86) * t630 + t578; t103 * m(7); t84 * t747 + t20 * t744 + t85 * t746 + t21 * t743 + (t103 * t117 + t106 * t145 + t107 * t144 + t205 * t70 + t228 * t97 + t229 * t96) * m(7) + t14 * t735 + t95 * t745 + t23 * t742 + (t3 * t736 + t4 * t737) * t565; (t103 * t118 + t106 * t148 + t107 * t149 + t205 * t42 + t228 * t59 + t229 * t60) * m(7) + t583; (t103 * t143 + t106 * t159 + t107 * t160 + t205 * t54 + t228 * t68 + t229 * t69) * m(7) + t583; (t101 * t228 + t102 * t229 + t103 * t177 + t106 * t190 + t107 * t191 + t205 * t86) * m(7) + t583; t441 * t76 - t582 * t4 + t439 * t75 - t581 * t3 + t477 * t88 + t523 * t14 + (t103 * t205 + t106 * t228 + t107 * t229) * t630;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
