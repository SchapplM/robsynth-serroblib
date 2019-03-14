% Calculate time derivative of joint inertia matrix for
% S6PRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR2_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_inertiaDJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR2_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR2_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPR2_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:05:25
% EndTime: 2019-03-08 23:06:15
% DurationCPUTime: 32.53s
% Computational Cost: add. (130334->1261), mult. (204180->1730), div. (0->0), fcn. (230541->14), ass. (0->544)
t549 = sin(pkin(11));
t552 = cos(pkin(11));
t558 = cos(qJ(2));
t553 = cos(pkin(6));
t556 = sin(qJ(2));
t682 = t553 * t556;
t530 = t549 * t558 + t552 * t682;
t678 = qJ(3) + qJ(4);
t614 = cos(t678);
t545 = sin(t678);
t550 = sin(pkin(6));
t687 = t550 * t552;
t638 = t545 * t687;
t506 = t530 * t614 - t638;
t548 = sin(pkin(12));
t551 = cos(pkin(12));
t680 = t553 * t558;
t575 = -t549 * t556 + t552 * t680;
t452 = -t506 * t548 - t551 * t575;
t690 = t575 * t548;
t453 = t506 * t551 - t690;
t586 = t550 * t614;
t566 = -t530 * t545 - t552 * t586;
t319 = Icges(6,5) * t453 + Icges(6,6) * t452 - Icges(6,3) * t566;
t408 = Icges(5,4) * t506 + Icges(5,2) * t566 - Icges(5,6) * t575;
t741 = t319 - t408;
t640 = t549 * t682;
t532 = t552 * t558 - t640;
t688 = t549 * t550;
t508 = t532 * t614 + t545 * t688;
t531 = t549 * t680 + t552 * t556;
t454 = -t508 * t548 + t531 * t551;
t689 = t531 * t548;
t455 = t508 * t551 + t689;
t567 = -t532 * t545 + t549 * t586;
t320 = Icges(6,5) * t455 + Icges(6,6) * t454 - Icges(6,3) * t567;
t409 = Icges(5,4) * t508 + Icges(5,2) * t567 + Icges(5,6) * t531;
t740 = t320 - t409;
t583 = t556 * t586;
t519 = t553 * t545 + t583;
t684 = t550 * t558;
t503 = -t519 * t548 - t551 * t684;
t642 = t548 * t684;
t504 = t519 * t551 - t642;
t518 = t550 * t556 * t545 - t553 * t614;
t390 = Icges(6,5) * t504 + Icges(6,6) * t503 + Icges(6,3) * t518;
t463 = Icges(5,4) * t519 - Icges(5,2) * t518 - Icges(5,6) * t684;
t739 = t390 - t463;
t322 = Icges(6,4) * t455 + Icges(6,2) * t454 - Icges(6,6) * t567;
t324 = Icges(6,1) * t455 + Icges(6,4) * t454 - Icges(6,5) * t567;
t407 = Icges(5,5) * t508 + Icges(5,6) * t567 + Icges(5,3) * t531;
t411 = Icges(5,1) * t508 + Icges(5,4) * t567 + Icges(5,5) * t531;
t734 = t322 * t452 + t324 * t453 - t407 * t575 + t411 * t506 - t566 * t740;
t321 = Icges(6,4) * t453 + Icges(6,2) * t452 - Icges(6,6) * t566;
t323 = Icges(6,1) * t453 + Icges(6,4) * t452 - Icges(6,5) * t566;
t406 = Icges(5,5) * t506 + Icges(5,6) * t566 - Icges(5,3) * t575;
t410 = Icges(5,1) * t506 + Icges(5,4) * t566 - Icges(5,5) * t575;
t738 = t321 * t454 + t323 * t455 + t406 * t531 + t410 * t508 - t567 * t741;
t732 = t322 * t454 + t324 * t455 + t407 * t531 + t411 * t508 - t567 * t740;
t730 = t322 * t503 + t324 * t504 - t407 * t684 + t411 * t519 + t518 * t740;
t737 = t321 * t452 + t323 * t453 - t406 * t575 + t410 * t506 - t566 * t741;
t736 = t321 * t503 + t323 * t504 - t406 * t684 + t410 * t519 + t518 * t741;
t391 = Icges(6,4) * t504 + Icges(6,2) * t503 + Icges(6,6) * t518;
t392 = Icges(6,1) * t504 + Icges(6,4) * t503 + Icges(6,5) * t518;
t462 = Icges(5,5) * t519 - Icges(5,6) * t518 - Icges(5,3) * t684;
t464 = Icges(5,1) * t519 - Icges(5,4) * t518 - Icges(5,5) * t684;
t729 = t391 * t452 + t392 * t453 - t462 * t575 + t464 * t506 - t566 * t739;
t728 = t391 * t454 + t392 * t455 + t462 * t531 + t464 * t508 - t567 * t739;
t726 = t391 * t503 + t392 * t504 - t462 * t684 + t464 * t519 + t518 * t739;
t520 = t575 * qJD(2);
t547 = qJD(3) + qJD(4);
t440 = t520 * t614 + t547 * t566;
t546 = pkin(12) + qJ(6);
t543 = sin(t546);
t544 = cos(t546);
t446 = t506 * t544 - t543 * t575;
t521 = t530 * qJD(2);
t315 = -qJD(6) * t446 - t440 * t543 + t521 * t544;
t445 = -t506 * t543 - t544 * t575;
t316 = qJD(6) * t445 + t440 * t544 + t521 * t543;
t587 = t547 * t614;
t439 = t520 * t545 + t530 * t587 - t547 * t638;
t206 = rSges(7,1) * t316 + rSges(7,2) * t315 + rSges(7,3) * t439;
t694 = t521 * t548;
t701 = pkin(5) * t551;
t232 = pkin(5) * t694 + pkin(10) * t439 + t440 * t701;
t676 = t206 + t232;
t522 = t531 * qJD(2);
t442 = -t522 * t614 + t547 * t567;
t448 = t508 * t544 + t531 * t543;
t647 = qJD(2) * t558;
t523 = -qJD(2) * t640 + t552 * t647;
t317 = -qJD(6) * t448 - t442 * t543 + t523 * t544;
t447 = -t508 * t543 + t531 * t544;
t318 = qJD(6) * t447 + t442 * t544 + t523 * t543;
t441 = t532 * t587 + (t547 * t688 - t522) * t545;
t207 = rSges(7,1) * t318 + rSges(7,2) * t317 + rSges(7,3) * t441;
t692 = t523 * t548;
t675 = pkin(5) * t692 + pkin(10) * t441 + t442 * t701 + t207;
t498 = -t518 * t547 + t586 * t647;
t577 = -t519 * t544 + t543 * t684;
t648 = qJD(2) * t556;
t619 = t550 * t648;
t382 = qJD(6) * t577 - t498 * t543 + t544 * t619;
t499 = -t519 * t543 - t544 * t684;
t383 = qJD(6) * t499 + t498 * t544 + t543 * t619;
t618 = t550 * t647;
t497 = t547 * t583 + (t547 * t553 + t618) * t545;
t258 = rSges(7,1) * t383 + rSges(7,2) * t382 + rSges(7,3) * t497;
t591 = t548 * t619;
t727 = -pkin(5) * t591 - pkin(10) * t497 - t498 * t701 - t258;
t299 = -pkin(5) * t690 - pkin(10) * t566 + t506 * t701;
t310 = rSges(7,1) * t446 + rSges(7,2) * t445 - rSges(7,3) * t566;
t670 = t299 + t310;
t311 = rSges(7,1) * t448 + rSges(7,2) * t447 - rSges(7,3) * t567;
t669 = pkin(5) * t689 - pkin(10) * t567 + t508 * t701 + t311;
t381 = -rSges(7,1) * t577 + rSges(7,2) * t499 + rSges(7,3) * t518;
t725 = -pkin(5) * t642 + pkin(10) * t518 + t519 * t701 + t381;
t334 = rSges(5,1) * t442 - rSges(5,2) * t441 + rSges(5,3) * t523;
t413 = rSges(5,1) * t508 + rSges(5,2) * t567 + rSges(5,3) * t531;
t333 = rSges(5,1) * t440 - rSges(5,2) * t439 + rSges(5,3) * t521;
t412 = rSges(5,1) * t506 + rSges(5,2) * t566 - rSges(5,3) * t575;
t671 = t531 * t333 + t523 * t412;
t174 = t334 * t575 - t521 * t413 + t671;
t460 = -t498 * t548 + t551 * t619;
t461 = t498 * t551 + t591;
t335 = Icges(6,5) * t461 + Icges(6,6) * t460 + Icges(6,3) * t497;
t336 = Icges(6,4) * t461 + Icges(6,2) * t460 + Icges(6,6) * t497;
t337 = Icges(6,1) * t461 + Icges(6,4) * t460 + Icges(6,5) * t497;
t397 = -t442 * t548 + t523 * t551;
t398 = t442 * t551 + t692;
t118 = -t335 * t567 + t336 * t454 + t337 * t455 + t390 * t441 + t391 * t397 + t392 * t398;
t255 = Icges(7,5) * t383 + Icges(7,6) * t382 + Icges(7,3) * t497;
t256 = Icges(7,4) * t383 + Icges(7,2) * t382 + Icges(7,6) * t497;
t257 = Icges(7,1) * t383 + Icges(7,4) * t382 + Icges(7,5) * t497;
t378 = -Icges(7,5) * t577 + Icges(7,6) * t499 + Icges(7,3) * t518;
t379 = -Icges(7,4) * t577 + Icges(7,2) * t499 + Icges(7,6) * t518;
t380 = -Icges(7,1) * t577 + Icges(7,4) * t499 + Icges(7,5) * t518;
t104 = -t255 * t567 + t256 * t447 + t257 * t448 + t317 * t379 + t318 * t380 + t378 * t441;
t304 = Icges(7,5) * t446 + Icges(7,6) * t445 - Icges(7,3) * t566;
t306 = Icges(7,4) * t446 + Icges(7,2) * t445 - Icges(7,6) * t566;
t308 = Icges(7,1) * t446 + Icges(7,4) * t445 - Icges(7,5) * t566;
t161 = -t304 * t567 + t306 * t447 + t308 * t448;
t305 = Icges(7,5) * t448 + Icges(7,6) * t447 - Icges(7,3) * t567;
t307 = Icges(7,4) * t448 + Icges(7,2) * t447 - Icges(7,6) * t567;
t309 = Icges(7,1) * t448 + Icges(7,4) * t447 - Icges(7,5) * t567;
t162 = -t305 * t567 + t307 * t447 + t309 * t448;
t185 = -t378 * t567 + t379 * t447 + t380 * t448;
t200 = Icges(7,5) * t316 + Icges(7,6) * t315 + Icges(7,3) * t439;
t202 = Icges(7,4) * t316 + Icges(7,2) * t315 + Icges(7,6) * t439;
t204 = Icges(7,1) * t316 + Icges(7,4) * t315 + Icges(7,5) * t439;
t63 = -t200 * t567 + t202 * t447 + t204 * t448 + t304 * t441 + t306 * t317 + t308 * t318;
t201 = Icges(7,5) * t318 + Icges(7,6) * t317 + Icges(7,3) * t441;
t203 = Icges(7,4) * t318 + Icges(7,2) * t317 + Icges(7,6) * t441;
t205 = Icges(7,1) * t318 + Icges(7,4) * t317 + Icges(7,5) * t441;
t64 = -t201 * t567 + t203 * t447 + t205 * t448 + t305 * t441 + t307 * t317 + t309 * t318;
t12 = t161 * t521 + t162 * t523 - t575 * t63 + t531 * t64 + (-t104 * t558 + t185 * t648) * t550;
t325 = Icges(5,5) * t440 - Icges(5,6) * t439 + Icges(5,3) * t521;
t327 = Icges(5,4) * t440 - Icges(5,2) * t439 + Icges(5,6) * t521;
t329 = Icges(5,1) * t440 - Icges(5,4) * t439 + Icges(5,5) * t521;
t131 = t325 * t531 + t327 * t567 + t329 * t508 + t406 * t523 - t408 * t441 + t410 * t442;
t326 = Icges(5,5) * t442 - Icges(5,6) * t441 + Icges(5,3) * t523;
t328 = Icges(5,4) * t442 - Icges(5,2) * t441 + Icges(5,6) * t523;
t330 = Icges(5,1) * t442 - Icges(5,4) * t441 + Icges(5,5) * t523;
t132 = t326 * t531 + t328 * t567 + t330 * t508 + t407 * t523 - t409 * t441 + t411 * t442;
t399 = Icges(5,5) * t498 - Icges(5,6) * t497 + Icges(5,3) * t619;
t400 = Icges(5,4) * t498 - Icges(5,2) * t497 + Icges(5,6) * t619;
t401 = Icges(5,1) * t498 - Icges(5,4) * t497 + Icges(5,5) * t619;
t158 = t399 * t531 + t400 * t567 + t401 * t508 - t441 * t463 + t442 * t464 + t462 * t523;
t395 = -t440 * t548 + t521 * t551;
t396 = t440 * t551 + t694;
t245 = Icges(6,5) * t396 + Icges(6,6) * t395 + Icges(6,3) * t439;
t247 = Icges(6,4) * t396 + Icges(6,2) * t395 + Icges(6,6) * t439;
t249 = Icges(6,1) * t396 + Icges(6,4) * t395 + Icges(6,5) * t439;
t97 = -t245 * t567 + t247 * t454 + t249 * t455 + t319 * t441 + t321 * t397 + t323 * t398;
t246 = Icges(6,5) * t398 + Icges(6,6) * t397 + Icges(6,3) * t441;
t248 = Icges(6,4) * t398 + Icges(6,2) * t397 + Icges(6,6) * t441;
t250 = Icges(6,1) * t398 + Icges(6,4) * t397 + Icges(6,5) * t441;
t98 = -t246 * t567 + t248 * t454 + t250 * t455 + t320 * t441 + t322 * t397 + t324 * t398;
t724 = t12 + (-t97 - t131) * t575 + (t728 * t648 + (-t118 - t158) * t558) * t550 + (t98 + t132) * t531 + t732 * t523 + t738 * t521;
t103 = -t255 * t566 + t256 * t445 + t257 * t446 + t315 * t379 + t316 * t380 + t378 * t439;
t159 = -t304 * t566 + t306 * t445 + t308 * t446;
t160 = -t305 * t566 + t307 * t445 + t309 * t446;
t184 = -t378 * t566 + t379 * t445 + t380 * t446;
t61 = -t200 * t566 + t202 * t445 + t204 * t446 + t304 * t439 + t306 * t315 + t308 * t316;
t62 = -t201 * t566 + t203 * t445 + t205 * t446 + t305 * t439 + t307 * t315 + t309 * t316;
t11 = t159 * t521 + t160 * t523 - t575 * t61 + t531 * t62 + (-t103 * t558 + t184 * t648) * t550;
t117 = -t335 * t566 + t336 * t452 + t337 * t453 + t390 * t439 + t391 * t395 + t392 * t396;
t129 = -t325 * t575 + t327 * t566 + t329 * t506 + t406 * t521 - t408 * t439 + t410 * t440;
t130 = -t326 * t575 + t328 * t566 + t330 * t506 + t407 * t521 - t409 * t439 + t411 * t440;
t157 = -t399 * t575 + t400 * t566 + t401 * t506 - t439 * t463 + t440 * t464 + t462 * t521;
t95 = -t245 * t566 + t247 * t452 + t249 * t453 + t319 * t439 + t321 * t395 + t323 * t396;
t96 = -t246 * t566 + t248 * t452 + t250 * t453 + t320 * t439 + t322 * t395 + t324 * t396;
t723 = t11 + (-t129 - t95) * t575 + (t729 * t648 + (-t117 - t157) * t558) * t550 + (t130 + t96) * t531 + t734 * t523 + t737 * t521;
t109 = t245 * t518 + t247 * t503 + t249 * t504 + t319 * t497 + t321 * t460 + t323 * t461;
t110 = t246 * t518 + t248 * t503 + t250 * t504 + t320 * t497 + t322 * t460 + t324 * t461;
t128 = t335 * t518 + t336 * t503 + t337 * t504 + t390 * t497 + t391 * t460 + t392 * t461;
t141 = -t327 * t518 + t329 * t519 - t408 * t497 + t410 * t498 + (-t325 * t558 + t406 * t648) * t550;
t142 = -t328 * t518 + t330 * t519 - t409 * t497 + t411 * t498 + (-t326 * t558 + t407 * t648) * t550;
t112 = t255 * t518 + t256 * t499 - t257 * t577 + t378 * t497 + t379 * t382 + t380 * t383;
t169 = t304 * t518 + t306 * t499 - t308 * t577;
t170 = t305 * t518 + t307 * t499 - t309 * t577;
t197 = t378 * t518 + t379 * t499 - t380 * t577;
t68 = t200 * t518 + t202 * t499 - t204 * t577 + t304 * t497 + t306 * t382 + t308 * t383;
t69 = t201 * t518 + t203 * t499 - t205 * t577 + t305 * t497 + t307 * t382 + t309 * t383;
t17 = t169 * t521 + t170 * t523 - t575 * t68 + t531 * t69 + (-t112 * t558 + t197 * t648) * t550;
t171 = -t400 * t518 + t401 * t519 - t463 * t497 + t464 * t498 + (-t399 * t558 + t462 * t648) * t550;
t722 = t17 + (-t141 - t109) * t575 + (t726 * t648 + (-t128 - t171) * t558) * t550 + (t142 + t110) * t531 + t730 * t523 + t736 * t521;
t493 = t520 * pkin(2) + t521 * pkin(8);
t494 = -t522 * pkin(2) + t523 * pkin(8);
t651 = t493 * t688 + t494 * t687;
t716 = m(7) / 0.2e1;
t717 = m(6) / 0.2e1;
t644 = t717 + t716;
t720 = -0.2e1 * t521;
t719 = 0.2e1 * t575;
t718 = m(5) / 0.2e1;
t715 = t439 / 0.2e1;
t714 = t441 / 0.2e1;
t713 = t497 / 0.2e1;
t712 = -t566 / 0.2e1;
t711 = -t567 / 0.2e1;
t710 = t518 / 0.2e1;
t709 = t521 / 0.2e1;
t708 = t523 / 0.2e1;
t707 = -t575 / 0.2e1;
t706 = t531 / 0.2e1;
t705 = t549 / 0.2e1;
t704 = -t552 / 0.2e1;
t703 = t553 / 0.2e1;
t557 = cos(qJ(3));
t702 = pkin(3) * t557;
t699 = pkin(3) * qJD(3);
t698 = Icges(3,4) * t556;
t697 = Icges(3,4) * t558;
t555 = sin(qJ(3));
t679 = t555 * t556;
t681 = t553 * t557;
t533 = -t550 * t679 + t681;
t683 = t553 * t555;
t685 = t550 * t557;
t534 = t556 * t685 + t683;
t477 = Icges(4,5) * t534 + Icges(4,6) * t533 - Icges(4,3) * t684;
t478 = Icges(4,4) * t534 + Icges(4,2) * t533 - Icges(4,6) * t684;
t479 = Icges(4,1) * t534 + Icges(4,4) * t533 - Icges(4,5) * t684;
t509 = -t530 * t555 - t552 * t685;
t686 = t550 * t555;
t639 = t552 * t686;
t576 = -t530 * t557 + t639;
t260 = -t477 * t575 + t478 * t509 - t479 * t576;
t696 = t521 * t260;
t511 = -t532 * t555 + t549 * t685;
t641 = t549 * t686;
t512 = t532 * t557 + t641;
t261 = t477 * t531 + t478 * t511 + t479 * t512;
t693 = t523 * t261;
t252 = rSges(6,1) * t398 + rSges(6,2) * t397 + rSges(6,3) * t441;
t298 = pkin(4) * t442 + qJ(5) * t441 - qJD(5) * t567;
t674 = -t252 - t298;
t297 = pkin(4) * t440 + qJ(5) * t439 - qJD(5) * t566;
t272 = t531 * t297;
t435 = pkin(4) * t506 - qJ(5) * t566;
t405 = t523 * t435;
t673 = t272 + t405;
t331 = rSges(6,1) * t453 + rSges(6,2) * t452 - rSges(6,3) * t566;
t416 = t531 * t435;
t672 = t531 * t331 + t416;
t332 = rSges(6,1) * t455 + rSges(6,2) * t454 - rSges(6,3) * t567;
t436 = pkin(4) * t508 - qJ(5) * t567;
t425 = t436 * t619;
t668 = t332 * t619 + t425;
t667 = -t331 - t435;
t666 = -t332 - t436;
t569 = t511 * qJD(3);
t360 = pkin(3) * t569 + pkin(9) * t523 - t522 * t702;
t665 = -t334 - t360;
t338 = rSges(6,1) * t461 + rSges(6,2) * t460 + rSges(6,3) * t497;
t376 = pkin(4) * t498 + qJ(5) * t497 + qJD(5) * t518;
t664 = -t338 - t376;
t570 = t509 * qJD(3);
t359 = pkin(3) * t570 + pkin(9) * t521 + t520 * t702;
t343 = t531 * t359;
t414 = -pkin(3) * t639 - pkin(9) * t575 + t530 * t702;
t372 = t523 * t414;
t663 = t343 + t372;
t470 = t553 * t494;
t662 = t553 * t360 + t470;
t661 = -t359 - t493;
t482 = pkin(3) * t683 + (-pkin(9) * t558 + t556 * t702) * t550;
t660 = t414 * t684 - t482 * t575;
t393 = rSges(6,1) * t504 + rSges(6,2) * t503 + rSges(6,3) * t518;
t480 = pkin(4) * t519 + qJ(5) * t518;
t659 = -t393 - t480;
t415 = pkin(3) * t641 + pkin(9) * t531 + t532 * t702;
t502 = t532 * pkin(2) + t531 * pkin(8);
t495 = t553 * t502;
t658 = t553 * t415 + t495;
t404 = rSges(5,1) * t498 - rSges(5,2) * t497 + rSges(5,3) * t619;
t449 = t681 * t699 + (-t679 * t699 + (pkin(9) * t556 + t558 * t702) * qJD(2)) * t550;
t657 = -t404 - t449;
t656 = -t412 - t414;
t655 = -t413 - t415;
t654 = t435 * t684 - t480 * t575;
t465 = rSges(5,1) * t519 - rSges(5,2) * t518 - rSges(5,3) * t684;
t302 = t412 * t684 - t465 * t575;
t653 = -t465 - t482;
t652 = 0.2e1 * t651;
t501 = t530 * pkin(2) - pkin(8) * t575;
t650 = t501 * t688 + t502 * t687;
t649 = qJD(2) * t550;
t637 = -t298 - t675;
t636 = -t360 + t674;
t635 = -t376 + t727;
t634 = t531 * t670 + t416;
t633 = t619 * t669 + t425;
t632 = t297 * t684 - t376 * t575 + t521 * t480;
t631 = t553 * t298 + t662;
t630 = -t297 + t661;
t629 = -t435 - t670;
t628 = -t436 - t669;
t627 = t333 * t684 - t404 * t575 + t521 * t465;
t626 = -t414 + t667;
t625 = -t415 + t666;
t624 = -t449 + t664;
t623 = t359 * t684 - t449 * t575 + t521 * t482;
t622 = -t480 - t725;
t621 = -t482 + t659;
t620 = t553 * t436 + t658;
t615 = -t684 / 0.2e1;
t612 = 2 * m(4);
t610 = 0.2e1 * m(5);
t608 = 0.2e1 * m(6);
t606 = 0.2e1 * m(7);
t513 = -qJD(3) * t534 - t555 * t618;
t514 = qJD(3) * t533 + t557 * t618;
t434 = rSges(4,1) * t514 + rSges(4,2) * t513 + rSges(4,3) * t619;
t528 = (pkin(2) * t558 + pkin(8) * t556) * t649;
t605 = (-t434 - t528) * t550;
t481 = rSges(4,1) * t534 + rSges(4,2) * t533 - rSges(4,3) * t684;
t535 = (pkin(2) * t556 - pkin(8) * t558) * t550;
t604 = (-t481 - t535) * t550;
t603 = -t360 + t637;
t251 = rSges(6,1) * t396 + rSges(6,2) * t395 + rSges(6,3) * t439;
t240 = t531 * t251;
t288 = t523 * t331;
t602 = t240 + t288 + t673;
t601 = -t449 + t635;
t600 = t298 * t719 + t436 * t720 + 0.2e1 * t272 + 0.2e1 * t405;
t599 = -t414 + t629;
t598 = -t415 + t628;
t597 = t360 * t719 + t415 * t720 + 0.2e1 * t343 + 0.2e1 * t372;
t354 = t359 * t688;
t355 = t360 * t687;
t596 = 0.2e1 * t354 + 0.2e1 * t355 + t652;
t595 = t354 + t355 + t651;
t594 = 0.2e1 * t174;
t593 = -t482 + t622;
t592 = t414 * t688 + t415 * t687 + t650;
t192 = t331 * t684 - t393 * t575 + t654;
t590 = t619 / 0.2e1;
t589 = (-t528 + t657) * t550;
t588 = (-t535 + t653) * t550;
t585 = (-t528 + t624) * t550;
t584 = (-t535 + t621) * t550;
t196 = t531 * t206;
t230 = t531 * t232;
t269 = t523 * t299;
t276 = t523 * t310;
t582 = t196 + t230 + t269 + t276 + t673;
t581 = t251 * t684 - t338 * t575 + t521 * t393 + t632;
t282 = t297 * t688;
t283 = t298 * t687;
t579 = t282 + t283 + t595;
t578 = t435 * t688 + t436 * t687 + t592;
t154 = -t575 * t725 + t670 * t684 + t654;
t574 = (-t528 + t601) * t550;
t573 = (-t535 + t593) * t550;
t568 = t252 * t719 + t332 * t720 + 0.2e1 * t240 + 0.2e1 * t288 + t600;
t565 = t521 * t725 + t575 * t727 + t676 * t684 + t632;
t100 = -t169 * t575 + t170 * t531 - t197 * t684;
t14 = t112 * t518 + t169 * t439 + t170 * t441 + t197 * t497 - t566 * t68 - t567 * t69;
t3 = t103 * t518 + t159 * t439 + t160 * t441 + t184 * t497 - t566 * t61 - t567 * t62;
t4 = t104 * t518 + t161 * t439 + t162 * t441 + t185 * t497 - t566 * t63 - t567 * t64;
t76 = -t159 * t566 - t160 * t567 + t184 * t518;
t77 = -t161 * t566 - t162 * t567 + t185 * t518;
t80 = -t159 * t575 + t160 * t531 - t184 * t684;
t81 = -t161 * t575 + t162 * t531 - t185 * t684;
t93 = -t169 * t566 - t170 * t567 + t197 * t518;
t564 = t100 * t713 + t11 * t712 + t12 * t711 + t14 * t615 + t17 * t710 + t3 * t707 + t4 * t706 + t93 * t590 + t77 * t708 + t76 * t709 + t81 * t714 + t80 * t715;
t113 = -t206 * t567 + t207 * t566 + t310 * t441 - t311 * t439;
t456 = qJD(3) * t576 - t520 * t555;
t457 = t520 * t557 + t570;
t350 = rSges(4,1) * t457 + rSges(4,2) * t456 + rSges(4,3) * t521;
t458 = -qJD(3) * t512 + t522 * t555;
t459 = -t522 * t557 + t569;
t351 = rSges(4,1) * t459 + rSges(4,2) * t458 + rSges(4,3) * t523;
t423 = -rSges(4,1) * t576 + rSges(4,2) * t509 - rSges(4,3) * t575;
t424 = rSges(4,1) * t512 + rSges(4,2) * t511 + rSges(4,3) * t531;
t182 = t350 * t531 + t351 * t575 + t423 * t523 - t424 * t521;
t563 = t669 * t720 + t675 * t719 + 0.2e1 * t196 + 0.2e1 * t230 + 0.2e1 * t269 + 0.2e1 * t276 + t600;
t562 = (-t684 * t726 + t100) * t619 + (t619 * t730 + t724) * t531 + (-t619 * t736 - t723) * t575 + (t531 * t732 - t575 * t738 - t684 * t728 + t81) * t523 + (t734 * t531 - t575 * t737 - t729 * t684 + t80) * t521;
t561 = -t684 * t722 + t562;
t102 = t197 * t553 + (-t169 * t552 + t170 * t549) * t550;
t26 = t103 * t553 + (t549 * t62 - t552 * t61) * t550;
t27 = t104 * t553 + (t549 * t64 - t552 * t63) * t550;
t31 = t112 * t553 + (t549 * t69 - t552 * t68) * t550;
t34 = t117 * t553 + (t549 * t96 - t552 * t95) * t550;
t35 = t118 * t553 + (t549 * t98 - t552 * t97) * t550;
t37 = t128 * t553 + (-t109 * t552 + t110 * t549) * t550;
t53 = t157 * t553 + (-t129 * t552 + t130 * t549) * t550;
t54 = t158 * t553 + (-t131 * t552 + t132 * t549) * t550;
t58 = t171 * t553 + (-t141 * t552 + t142 * t549) * t550;
t84 = t184 * t553 + (-t159 * t552 + t160 * t549) * t550;
t85 = t185 * t553 + (-t161 * t552 + t162 * t549) * t550;
t560 = (t84 + t729 * t553 + (t734 * t549 - t552 * t737) * t550) * t709 + (t85 + t728 * t553 + (t549 * t732 - t552 * t738) * t550) * t708 + (t53 + t34 + t26) * t707 + (t54 + t35 + t27) * t706 + t722 * t703 + t724 * t688 / 0.2e1 - t723 * t687 / 0.2e1 + (t37 + t31 + t58) * t615 + (t102 + t726 * t553 + (t549 * t730 - t552 * t736) * t550) * t590;
t527 = (rSges(3,1) * t558 - rSges(3,2) * t556) * t649;
t526 = (Icges(3,1) * t558 - t698) * t649;
t525 = (-Icges(3,2) * t556 + t697) * t649;
t524 = (Icges(3,5) * t558 - Icges(3,6) * t556) * t649;
t517 = rSges(3,3) * t553 + (rSges(3,1) * t556 + rSges(3,2) * t558) * t550;
t516 = Icges(3,5) * t553 + (Icges(3,1) * t556 + t697) * t550;
t515 = Icges(3,6) * t553 + (Icges(3,2) * t558 + t698) * t550;
t492 = -rSges(3,1) * t522 - rSges(3,2) * t523;
t491 = rSges(3,1) * t520 - rSges(3,2) * t521;
t490 = -Icges(3,1) * t522 - Icges(3,4) * t523;
t489 = Icges(3,1) * t520 - Icges(3,4) * t521;
t488 = -Icges(3,4) * t522 - Icges(3,2) * t523;
t487 = Icges(3,4) * t520 - Icges(3,2) * t521;
t486 = -Icges(3,5) * t522 - Icges(3,6) * t523;
t485 = Icges(3,5) * t520 - Icges(3,6) * t521;
t476 = rSges(3,1) * t532 - rSges(3,2) * t531 + rSges(3,3) * t688;
t475 = rSges(3,1) * t530 + rSges(3,2) * t575 - rSges(3,3) * t687;
t474 = Icges(3,1) * t532 - Icges(3,4) * t531 + Icges(3,5) * t688;
t473 = Icges(3,1) * t530 + Icges(3,4) * t575 - Icges(3,5) * t687;
t472 = Icges(3,4) * t532 - Icges(3,2) * t531 + Icges(3,6) * t688;
t471 = Icges(3,4) * t530 + Icges(3,2) * t575 - Icges(3,6) * t687;
t433 = Icges(4,1) * t514 + Icges(4,4) * t513 + Icges(4,5) * t619;
t432 = Icges(4,4) * t514 + Icges(4,2) * t513 + Icges(4,6) * t619;
t431 = Icges(4,5) * t514 + Icges(4,6) * t513 + Icges(4,3) * t619;
t422 = Icges(4,1) * t512 + Icges(4,4) * t511 + Icges(4,5) * t531;
t421 = -Icges(4,1) * t576 + Icges(4,4) * t509 - Icges(4,5) * t575;
t420 = Icges(4,4) * t512 + Icges(4,2) * t511 + Icges(4,6) * t531;
t419 = -Icges(4,4) * t576 + Icges(4,2) * t509 - Icges(4,6) * t575;
t418 = Icges(4,5) * t512 + Icges(4,6) * t511 + Icges(4,3) * t531;
t417 = -Icges(4,5) * t576 + Icges(4,6) * t509 - Icges(4,3) * t575;
t385 = t415 * t619;
t384 = t413 * t619;
t375 = t531 * t414;
t374 = t531 * t412;
t349 = Icges(4,1) * t459 + Icges(4,4) * t458 + Icges(4,5) * t523;
t348 = Icges(4,1) * t457 + Icges(4,4) * t456 + Icges(4,5) * t521;
t347 = Icges(4,4) * t459 + Icges(4,2) * t458 + Icges(4,6) * t523;
t346 = Icges(4,4) * t457 + Icges(4,2) * t456 + Icges(4,6) * t521;
t345 = Icges(4,5) * t459 + Icges(4,6) * t458 + Icges(4,3) * t523;
t344 = Icges(4,5) * t457 + Icges(4,6) * t456 + Icges(4,3) * t521;
t340 = -t424 * t684 - t481 * t531;
t339 = t423 * t684 - t481 * t575;
t303 = -t413 * t684 - t465 * t531;
t292 = -t477 * t684 + t478 * t533 + t479 * t534;
t266 = t423 * t531 + t424 * t575;
t265 = (-t423 - t501) * t553 + t552 * t604;
t264 = t424 * t553 + t549 * t604 + t495;
t259 = t413 * t575 + t374;
t241 = (t423 * t549 + t424 * t552) * t550 + t650;
t237 = -t418 * t684 + t420 * t533 + t422 * t534;
t236 = -t417 * t684 + t419 * t533 + t421 * t534;
t235 = (-t350 - t493) * t553 + t552 * t605;
t234 = t351 * t553 + t549 * t605 + t470;
t225 = t311 * t518 + t381 * t567;
t224 = -t310 * t518 - t381 * t566;
t223 = t531 * t653 + t655 * t684;
t222 = t302 + t660;
t221 = t418 * t531 + t420 * t511 + t422 * t512;
t220 = t417 * t531 + t419 * t511 + t421 * t512;
t219 = -t418 * t575 + t420 * t509 - t422 * t576;
t218 = -t417 * t575 + t419 * t509 - t421 * t576;
t213 = (-t501 + t656) * t553 + t552 * t588;
t212 = t413 * t553 + t549 * t588 + t658;
t210 = (t350 * t549 + t351 * t552) * t550 + t651;
t209 = -t434 * t531 - t481 * t523 + (-t351 * t558 + t424 * t648) * t550;
t208 = -t434 * t575 + t481 * t521 + (t350 * t558 - t423 * t648) * t550;
t199 = -t310 * t567 + t311 * t566;
t193 = t531 * t659 + t666 * t684;
t191 = -t575 * t655 + t374 + t375;
t188 = -t334 * t684 - t404 * t531 - t465 * t523 + t384;
t187 = -t412 * t619 + t627;
t186 = t432 * t533 + t433 * t534 + t478 * t513 + t479 * t514 + (-t431 * t558 + t477 * t648) * t550;
t183 = (t412 * t549 + t413 * t552) * t550 + t592;
t181 = -t575 * t666 + t672;
t180 = (-t333 + t661) * t553 + t552 * t589;
t179 = t334 * t553 + t549 * t589 + t662;
t178 = t431 * t531 + t432 * t511 + t433 * t512 + t458 * t478 + t459 * t479 + t477 * t523;
t177 = -t431 * t575 + t432 * t509 - t433 * t576 + t456 * t478 + t457 * t479 + t477 * t521;
t173 = t531 * t621 + t625 * t684;
t172 = t192 + t660;
t168 = (-t501 + t626) * t553 + t552 * t584;
t167 = t332 * t553 + t549 * t584 + t620;
t156 = (t333 * t549 + t334 * t552) * t550 + t595;
t155 = t531 * t622 + t628 * t684;
t153 = -t575 * t625 + t375 + t672;
t152 = (t331 * t549 + t332 * t552) * t550 + t578;
t151 = t347 * t533 + t349 * t534 + t420 * t513 + t422 * t514 + (-t345 * t558 + t418 * t648) * t550;
t150 = t346 * t533 + t348 * t534 + t419 * t513 + t421 * t514 + (-t344 * t558 + t417 * t648) * t550;
t149 = t523 * t653 + t531 * t657 + t665 * t684 + t384 + t385;
t148 = t619 * t656 + t623 + t627;
t147 = t345 * t531 + t347 * t511 + t349 * t512 + t418 * t523 + t420 * t458 + t422 * t459;
t146 = t344 * t531 + t346 * t511 + t348 * t512 + t417 * t523 + t419 * t458 + t421 * t459;
t145 = -t345 * t575 + t347 * t509 - t349 * t576 + t418 * t521 + t420 * t456 + t422 * t457;
t144 = -t344 * t575 + t346 * t509 - t348 * t576 + t417 * t521 + t419 * t456 + t421 * t457;
t140 = t531 * t593 + t598 * t684;
t139 = t154 + t660;
t137 = (-t501 + t599) * t553 + t552 * t573;
t136 = t549 * t573 + t553 * t669 + t620;
t135 = -t575 * t628 + t634;
t134 = (-t251 + t630) * t553 + t552 * t585;
t133 = t252 * t553 + t549 * t585 + t631;
t123 = t207 * t518 + t258 * t567 + t311 * t497 - t381 * t441;
t122 = -t206 * t518 - t258 * t566 - t310 * t497 + t381 * t439;
t121 = t521 * t655 - t575 * t665 + t663 + t671;
t120 = -t575 * t598 + t375 + t634;
t119 = (t549 * t670 + t552 * t669) * t550 + t578;
t116 = t523 * t659 + t531 * t664 + t674 * t684 + t668;
t115 = t619 * t667 + t581;
t114 = (t251 * t549 + t252 * t552) * t550 + t579;
t111 = t521 * t666 - t575 * t674 + t602;
t107 = t523 * t621 + t531 * t624 + t636 * t684 + t385 + t668;
t106 = t619 * t626 + t581 + t623;
t91 = (t630 - t676) * t553 + t552 * t574;
t90 = t549 * t574 + t553 * t675 + t631;
t71 = t521 * t625 - t575 * t636 + t602 + t663;
t70 = t186 * t553 + (-t150 * t552 + t151 * t549) * t550;
t67 = (t549 * t676 + t552 * t675) * t550 + t579;
t66 = t523 * t622 + t531 * t635 + t637 * t684 + t633;
t65 = t619 * t629 + t565;
t60 = t178 * t553 + (-t146 * t552 + t147 * t549) * t550;
t59 = t177 * t553 + (-t144 * t552 + t145 * t549) * t550;
t56 = t523 * t593 + t531 * t601 + t603 * t684 + t385 + t633;
t55 = t599 * t619 + t565 + t623;
t50 = t521 * t628 - t575 * t637 + t582;
t49 = -t150 * t575 + t151 * t531 + t236 * t521 + t237 * t523 + (-t186 * t558 + t292 * t648) * t550;
t48 = -t146 * t575 + t147 * t531 + t220 * t521 + t523 * t221 + (-t178 * t558 + t261 * t648) * t550;
t47 = -t144 * t575 + t145 * t531 + t521 * t218 + t219 * t523 + (-t177 * t558 + t260 * t648) * t550;
t44 = t521 * t598 - t575 * t603 + t582 + t663;
t1 = [0; m(4) * t652 / 0.2e1 + t596 * t718 + (m(3) * t492 + m(4) * t351 + m(5) * t334 + m(6) * t252 + m(7) * t675) * t687 + (m(3) * t491 + m(4) * t350 + m(5) * t333 + m(6) * t251 + m(7) * t676) * t688 + t644 * (0.2e1 * t282 + 0.2e1 * t283 + t596); t27 * t688 - t26 * t687 + t35 * t688 + t54 * t688 - t53 * t687 - t34 * t687 - t59 * t687 + t60 * t688 - ((-t472 * t521 + t474 * t520 - t486 * t687 + t488 * t575 + t490 * t530) * t688 - (-t471 * t521 + t473 * t520 - t485 * t687 + t487 * t575 + t489 * t530) * t687 + (-t515 * t521 + t516 * t520 - t524 * t687 + t525 * t575 + t526 * t530) * t553) * t687 + ((-t472 * t523 - t474 * t522 + t486 * t688 - t488 * t531 + t490 * t532) * t688 - (-t471 * t523 - t473 * t522 + t485 * t688 - t487 * t531 + t489 * t532) * t687 + (-t515 * t523 - t516 * t522 + t524 * t688 - t525 * t531 + t526 * t532) * t553) * t688 + t553 * t31 + (t119 * t67 + t136 * t90 + t137 * t91) * t606 + (t114 * t152 + t133 * t167 + t134 * t168) * t608 + t553 * t58 + t553 * t37 + (t156 * t183 + t179 * t212 + t180 * t213) * t610 + (t210 * t241 + t234 * t264 + t235 * t265) * t612 + t553 * t70 + t553 * (t553 ^ 2 * t524 + (((t488 * t558 + t490 * t556) * t549 - (t487 * t558 + t489 * t556) * t552 + ((-t472 * t556 + t474 * t558) * t549 - (-t471 * t556 + t473 * t558) * t552) * qJD(2)) * t550 + (-t485 * t552 + t486 * t549 + t525 * t558 + t526 * t556 + (-t515 * t556 + t516 * t558) * qJD(2)) * t553) * t550) + 0.2e1 * m(3) * ((-t475 * t553 - t517 * t687) * (-t491 * t553 - t527 * t687) + (t476 * t553 - t517 * t688) * (t492 * t553 - t527 * t688) + (t475 * t549 + t476 * t552) * t550 ^ 2 * (t491 * t549 + t492 * t552)); t182 * m(4) + (t594 + t597) * t718 + (t568 + t597) * t717 + (t563 + t597) * t716; t60 * t706 + t59 * t707 + ((-t220 * t552 + t221 * t549) * t708 + (-t218 * t552 + t219 * t549) * t709 - t558 * t70 / 0.2e1 + t47 * t704 + t48 * t705 + (t292 * t703 + (-t236 * t552 + t237 * t549) * t550 / 0.2e1) * t648) * t550 + (t49 / 0.2e1 + t693 / 0.2e1 + t696 / 0.2e1) * t553 + t560 + (t106 * t168 + t107 * t167 + t114 * t153 + t133 * t173 + t134 * t172 + t152 * t71) * m(6) + (t119 * t44 + t120 * t67 + t136 * t56 + t137 * t55 + t139 * t91 + t140 * t90) * m(7) + (t121 * t183 + t148 * t213 + t149 * t212 + t156 * t191 + t179 * t223 + t180 * t222) * m(5) + (t182 * t241 + t208 * t265 + t209 * t264 + t210 * t266 + t234 * t340 + t235 * t339) * m(4); t562 + t531 * t48 - t575 * t47 + ((-t236 * t575 + t237 * t531) * t648 + (-t292 * t619 - t49 - t693 - t696 - t722) * t558) * t550 + t523 * (-t220 * t575 + t221 * t531) + t521 * (-t218 * t575 + t219 * t531) + (t106 * t172 + t107 * t173 + t153 * t71) * t608 + (t120 * t44 + t139 * t55 + t140 * t56) * t606 + (t121 * t191 + t148 * t222 + t149 * t223) * t610 + (t182 * t266 + t208 * t339 + t209 * t340) * t612; t563 * t716 + t568 * t717 + t594 * t718; t560 + (t111 * t152 + t114 * t181 + t115 * t168 + t116 * t167 + t133 * t193 + t134 * t192) * m(6) + (t119 * t50 + t135 * t67 + t136 * t66 + t137 * t65 + t154 * t91 + t155 * t90) * m(7) + (t156 * t259 + t174 * t183 + t179 * t303 + t180 * t302 + t187 * t213 + t188 * t212) * m(5); t561 + (t106 * t192 + t107 * t193 + t111 * t153 + t115 * t172 + t116 * t173 + t181 * t71) * m(6) + (t120 * t50 + t135 * t44 + t139 * t65 + t140 * t66 + t154 * t55 + t155 * t56) * m(7) + (t121 * t259 + t148 * t302 + t149 * t303 + t174 * t191 + t187 * t222 + t188 * t223) * m(5); t561 + (t111 * t181 + t115 * t192 + t116 * t193) * t608 + (t135 * t50 + t154 * t65 + t155 * t66) * t606 + (t174 * t259 + t187 * t302 + t188 * t303) * t610; 0.2e1 * t644 * t497; (m(6) * t114 + m(7) * t67) * t518 - (m(6) * t134 + m(7) * t91) * t567 - (m(6) * t133 + m(7) * t90) * t566 + (m(6) * t152 + m(7) * t119) * t497 + (m(6) * t168 + m(7) * t137) * t441 + (m(6) * t167 + m(7) * t136) * t439; (m(6) * t71 + m(7) * t44) * t518 - (m(6) * t106 + m(7) * t55) * t567 - (m(6) * t107 + m(7) * t56) * t566 + (m(6) * t153 + m(7) * t120) * t497 + (m(6) * t172 + m(7) * t139) * t441 + (m(6) * t173 + m(7) * t140) * t439; (m(6) * t111 + m(7) * t50) * t518 - (m(6) * t115 + m(7) * t65) * t567 - (m(6) * t116 + m(7) * t66) * t566 + (m(6) * t181 + m(7) * t135) * t497 + (m(6) * t192 + m(7) * t154) * t441 + (m(6) * t193 + m(7) * t155) * t439; 0.4e1 * t644 * (-t439 * t566 - t441 * t567 + t497 * t518); t113 * m(7); t85 * t714 + t27 * t711 + t14 * t703 + (t113 * t119 + t122 * t137 + t123 * t136 + t199 * t67 + t224 * t91 + t225 * t90) * m(7) + t102 * t713 + t31 * t710 + t84 * t715 + t26 * t712 + (t3 * t704 + t4 * t705) * t550; t564 + (t113 * t120 + t122 * t139 + t123 * t140 + t199 * t44 + t224 * t55 + t225 * t56) * m(7); t564 + (t113 * t135 + t122 * t154 + t123 * t155 + t199 * t50 + t224 * t65 + t225 * t66) * m(7); (t113 * t518 - t122 * t567 - t123 * t566 + t199 * t497 + t224 * t441 + t225 * t439) * m(7); t441 * t77 - t567 * t4 + t439 * t76 - t566 * t3 + t497 * t93 + t518 * t14 + (t113 * t199 + t122 * t224 + t123 * t225) * t606;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;