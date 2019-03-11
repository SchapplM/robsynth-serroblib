% Calculate vector of inverse dynamics joint torques for
% S6PRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRRP5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:20:15
% EndTime: 2019-03-09 00:21:27
% DurationCPUTime: 42.86s
% Computational Cost: add. (14253->908), mult. (35921->1247), div. (0->0), fcn. (30238->14), ass. (0->406)
t627 = Ifges(6,4) + Ifges(7,4);
t628 = Ifges(6,1) + Ifges(7,1);
t626 = -Ifges(7,5) - Ifges(6,5);
t625 = Ifges(6,2) + Ifges(7,2);
t624 = Ifges(7,6) + Ifges(6,6);
t512 = cos(pkin(7));
t324 = qJD(2) * t512 + qJD(3);
t338 = sin(qJ(4));
t342 = cos(qJ(4));
t339 = sin(qJ(3));
t333 = sin(pkin(7));
t486 = qJD(2) * t333;
t463 = t339 * t486;
t246 = t324 * t342 - t338 * t463;
t343 = cos(qJ(3));
t476 = qJD(2) * qJD(3);
t278 = (qJDD(2) * t339 + t343 * t476) * t333;
t323 = qJDD(2) * t512 + qJDD(3);
t162 = qJD(4) * t246 + t278 * t342 + t323 * t338;
t561 = t162 / 0.2e1;
t658 = Ifges(5,4) * t561;
t402 = pkin(4) * t338 - pkin(11) * t342;
t334 = sin(pkin(6));
t540 = cos(qJ(2));
t466 = t334 * t540;
t420 = qJD(1) * t466;
t306 = qJD(2) * pkin(2) + t420;
t430 = t512 * t306;
t485 = qJD(2) * t343;
t335 = cos(pkin(6));
t487 = qJD(1) * t335;
t340 = sin(qJ(2));
t498 = t334 * t340;
t464 = qJD(1) * t498;
t293 = pkin(9) * t486 + t464;
t491 = t343 * t293;
t657 = -t339 * t430 - t491 - (t339 * t487 + t402 * t485) * t333 + t402 * qJD(4);
t358 = t333 * t487 + t430;
t247 = t324 * t338 + t342 * t463;
t163 = -qJD(4) * t247 - t278 * t338 + t323 * t342;
t158 = qJDD(5) - t163;
t462 = t333 * t485;
t310 = qJD(4) - t462;
t337 = sin(qJ(5));
t341 = cos(qJ(5));
t202 = t247 * t341 + t310 * t337;
t483 = qJD(3) * t343;
t316 = qJD(2) * t420;
t280 = qJDD(1) * t498 + t316;
t640 = pkin(9) * qJDD(2) * t333 + qJD(3) * t358 + t280;
t461 = qJD(2) * t498;
t419 = qJD(1) * t461;
t279 = qJDD(1) * t466 - t419;
t255 = qJDD(2) * pkin(2) + t279;
t644 = qJDD(1) * t333 * t335 + t255 * t512;
t79 = -t293 * t483 - t339 * t640 + t343 * t644;
t66 = -pkin(3) * t323 - t79;
t30 = -pkin(4) * t163 - pkin(11) * t162 + t66;
t184 = t339 * t358 + t491;
t165 = pkin(10) * t324 + t184;
t434 = t335 * t512;
t320 = qJD(1) * t434;
t404 = -pkin(3) * t343 - pkin(10) * t339;
t197 = t320 + (qJD(2) * t404 - t306) * t333;
t95 = t165 * t342 + t197 * t338;
t83 = pkin(11) * t310 + t95;
t183 = -t293 * t339 + t343 * t358;
t164 = -pkin(3) * t324 - t183;
t97 = -pkin(4) * t246 - pkin(11) * t247 + t164;
t32 = t337 * t97 + t341 * t83;
t209 = qJDD(1) * t434 - t255 * t333;
t277 = (-qJDD(2) * t343 + t339 * t476) * t333;
t127 = pkin(3) * t277 - pkin(10) * t278 + t209;
t481 = qJD(4) * t342;
t482 = qJD(4) * t338;
t484 = qJD(3) * t339;
t78 = -t293 * t484 + t339 * t644 + t343 * t640;
t65 = pkin(10) * t323 + t78;
t17 = t127 * t338 - t165 * t482 + t197 * t481 + t342 * t65;
t264 = qJDD(4) + t277;
t9 = pkin(11) * t264 + t17;
t4 = -qJD(5) * t32 + t30 * t341 - t337 * t9;
t201 = -t247 * t337 + t310 * t341;
t67 = qJD(5) * t201 + t162 * t341 + t264 * t337;
t1 = pkin(5) * t158 - qJ(6) * t67 - qJD(6) * t202 + t4;
t545 = t264 / 0.2e1;
t560 = t163 / 0.2e1;
t562 = t158 / 0.2e1;
t68 = -qJD(5) * t202 - t162 * t337 + t264 * t341;
t566 = t68 / 0.2e1;
t567 = t67 / 0.2e1;
t623 = Ifges(7,3) + Ifges(6,3);
t478 = qJD(5) * t341;
t480 = qJD(5) * t337;
t3 = t30 * t337 + t341 * t9 + t478 * t97 - t480 * t83;
t2 = qJ(6) * t68 + qJD(6) * t201 + t3;
t636 = t4 * mrSges(6,1) - t3 * mrSges(6,2) - t2 * mrSges(7,2);
t656 = t1 * mrSges(7,1) - 0.2e1 * Ifges(5,2) * t560 - 0.2e1 * Ifges(5,6) * t545 + t623 * t562 + t624 * t566 - t626 * t567 + t636 - t658;
t332 = pkin(5) * t341 + pkin(4);
t655 = m(7) * t332;
t654 = t627 * t201;
t492 = t342 * t343;
t231 = (-t337 * t492 + t339 * t341) * t486;
t458 = t337 * t481;
t601 = t338 * t478 + t231 + t458;
t506 = t246 * t337;
t653 = t480 - t506;
t652 = t627 * t202;
t237 = qJD(5) - t246;
t31 = -t337 * t83 + t341 * t97;
t24 = -qJ(6) * t202 + t31;
t23 = pkin(5) * t237 + t24;
t25 = qJ(6) * t201 + t32;
t615 = t310 * Ifges(5,6);
t651 = -t164 * mrSges(5,1) - t31 * mrSges(6,1) - t23 * mrSges(7,1) + t32 * mrSges(6,2) + t25 * mrSges(7,2) + t615 / 0.2e1;
t336 = -qJ(6) - pkin(11);
t650 = -m(7) * t336 + mrSges(6,3) + mrSges(7,3);
t620 = -t158 * t626 + t627 * t68 + t628 * t67;
t649 = t620 / 0.2e1;
t621 = t158 * t624 + t625 * t68 + t627 * t67;
t648 = t621 / 0.2e1;
t618 = t201 * t625 + t237 * t624 + t652;
t617 = t202 * t628 - t237 * t626 + t654;
t647 = t246 * Ifges(5,2);
t375 = (pkin(3) * t339 - pkin(10) * t343) * t333;
t271 = qJD(2) * t375;
t129 = t183 * t342 + t271 * t338;
t115 = pkin(11) * t463 + t129;
t474 = pkin(10) * t482;
t646 = t657 * t341 + (t115 + t474) * t337;
t403 = -pkin(4) * t342 - pkin(11) * t338;
t309 = -pkin(3) + t403;
t645 = -t341 * t115 + t309 * t478 + t337 * t657;
t128 = -t183 * t338 + t271 * t342;
t114 = -pkin(4) * t463 - t128;
t473 = pkin(10) * t481;
t643 = t473 - t114;
t642 = t627 * t341;
t641 = t627 * t337;
t427 = mrSges(4,3) * t463;
t574 = -m(5) * t164 + mrSges(4,1) * t324 + mrSges(5,1) * t246 - mrSges(5,2) * t247 - t427;
t639 = -m(4) * t183 - t574;
t524 = Ifges(5,4) * t247;
t149 = t524 + t615 + t647;
t619 = t201 * t624 - t202 * t626 + t237 * t623;
t638 = -t95 * mrSges(5,3) - t149 / 0.2e1 + t619 / 0.2e1;
t637 = Ifges(5,1) * t561 + Ifges(5,5) * t545;
t569 = m(7) * pkin(5);
t633 = t78 * mrSges(4,2);
t632 = t79 * mrSges(4,1);
t631 = -mrSges(7,1) - mrSges(6,1);
t630 = mrSges(4,2) - mrSges(5,3);
t629 = mrSges(6,2) + mrSges(7,2);
t622 = t158 * t623 + t624 * t68 - t626 * t67;
t616 = t310 * Ifges(5,5);
t116 = mrSges(5,1) * t264 - mrSges(5,3) * t162;
t27 = -mrSges(6,1) * t68 + mrSges(6,2) * t67;
t614 = -t116 + t27;
t232 = (t337 * t339 + t341 * t492) * t486;
t493 = t341 * t342;
t330 = pkin(10) * t493;
t421 = t338 * t462;
t477 = qJD(6) * t341;
t613 = -pkin(5) * t421 + qJ(6) * t232 - t338 * t477 + (pkin(5) * t338 - qJ(6) * t493) * qJD(4) + (-t330 + (qJ(6) * t338 - t309) * t337) * qJD(5) + t646;
t495 = t338 * t341;
t612 = -qJ(6) * t231 + (-pkin(10) * qJD(4) - qJ(6) * qJD(5)) * t495 + (-qJD(6) * t338 + (-pkin(10) * qJD(5) - qJ(6) * qJD(4)) * t342) * t337 + t645;
t611 = (-t341 * t482 - t342 * t480) * pkin(10) + t645;
t257 = t309 * t337 + t330;
t610 = -qJD(5) * t257 + t646;
t609 = pkin(5) * t601 + t643;
t437 = qJD(5) * t336;
t182 = pkin(4) * t247 - pkin(11) * t246;
t94 = -t165 * t338 + t197 * t342;
t44 = t182 * t337 + t341 * t94;
t608 = qJ(6) * t506 + t337 * t437 - t44 + t477;
t43 = t182 * t341 - t337 * t94;
t505 = t246 * t341;
t607 = -pkin(5) * t247 + qJ(6) * t505 - qJD(6) * t337 + t341 * t437 - t43;
t606 = pkin(5) * t653 - t95;
t605 = t569 + mrSges(7,1);
t528 = mrSges(7,3) * t201;
t140 = -mrSges(7,2) * t237 + t528;
t530 = mrSges(6,3) * t201;
t141 = -mrSges(6,2) * t237 + t530;
t604 = -t140 - t141;
t527 = mrSges(7,3) * t202;
t142 = mrSges(7,1) * t237 - t527;
t529 = mrSges(6,3) * t202;
t143 = mrSges(6,1) * t237 - t529;
t603 = -t142 - t143;
t432 = t343 * t512;
t501 = t333 * t339;
t290 = pkin(2) * t432 - pkin(9) * t501;
t259 = -pkin(3) * t512 - t290;
t286 = t338 * t501 - t342 * t512;
t500 = t333 * t342;
t287 = t338 * t512 + t339 * t500;
t172 = pkin(4) * t286 - pkin(11) * t287 + t259;
t433 = t339 * t512;
t499 = t333 * t343;
t291 = pkin(2) * t433 + pkin(9) * t499;
t260 = pkin(10) * t512 + t291;
t261 = (-pkin(2) + t404) * t333;
t186 = t260 * t342 + t261 * t338;
t174 = -pkin(11) * t499 + t186;
t81 = t172 * t337 + t174 * t341;
t272 = qJD(3) * t375;
t273 = t290 * qJD(3);
t106 = -t260 * t482 + t261 * t481 + t272 * t338 + t273 * t342;
t355 = -t340 * t433 + t343 * t540;
t249 = t355 * t334;
t227 = qJD(1) * t249;
t423 = t333 * t464;
t196 = t227 * t342 + t338 * t423;
t602 = -t196 + t106;
t479 = qJD(5) * t338;
t600 = t337 * t479 - t341 * t481 + t232;
t599 = t273 - t227;
t511 = cos(pkin(12));
t407 = t511 * t540;
t510 = sin(pkin(12));
t428 = t510 * t340;
t349 = -t335 * t407 + t428;
t436 = t334 * t511;
t597 = t333 * t436 + t349 * t512;
t406 = t510 * t540;
t429 = t511 * t340;
t350 = t335 * t406 + t429;
t435 = t334 * t510;
t596 = -t333 * t435 + t350 * t512;
t397 = mrSges(7,1) * t337 + mrSges(7,2) * t341;
t399 = mrSges(6,1) * t337 + mrSges(6,2) * t341;
t82 = -pkin(4) * t310 - t94;
t47 = -pkin(5) * t201 + qJD(6) + t82;
t595 = t397 * t47 + t399 * t82;
t408 = t512 * t540;
t594 = -t339 * t340 + t343 * t408;
t593 = -t337 * t626 + t341 * t624;
t592 = -t337 * t624 - t341 * t626;
t591 = t341 * t625 + t641;
t590 = -t337 * t625 + t642;
t589 = t337 * t628 + t642;
t588 = t341 * t628 - t641;
t587 = t421 - t482;
t586 = -t478 + t505;
t584 = t3 * t341 - t337 * t4;
t583 = -m(7) - m(6) - m(5);
t234 = -t306 * t333 + t320;
t582 = (-t324 * (Ifges(4,5) * t343 - Ifges(4,6) * t339) / 0.2e1 - t234 * (mrSges(4,1) * t339 + mrSges(4,2) * t343)) * t333;
t581 = -mrSges(6,1) - t605;
t398 = -mrSges(7,1) * t341 + mrSges(7,2) * t337;
t400 = -mrSges(6,1) * t341 + mrSges(6,2) * t337;
t580 = m(6) * pkin(4) + mrSges(5,1) - t398 - t400 + t655;
t572 = -m(6) * pkin(11) + mrSges(5,2) - t650;
t401 = -mrSges(5,1) * t342 + mrSges(5,2) * t338;
t577 = -m(6) * t403 + t338 * t650 + t342 * t655 + mrSges(4,1) - t401;
t467 = pkin(5) * t337 + pkin(10);
t576 = -m(7) * t467 + t630;
t113 = -mrSges(6,1) * t201 + mrSges(6,2) * t202;
t531 = mrSges(5,3) * t247;
t206 = mrSges(5,1) * t310 - t531;
t575 = -m(6) * t82 - t113 + t206;
t112 = -mrSges(7,1) * t201 + mrSges(7,2) * t202;
t571 = m(5) * t94 - m(7) * t47 - t112 + t575;
t570 = t333 ^ 2;
t344 = qJD(2) ^ 2;
t568 = Ifges(5,4) * t560 + t637;
t559 = -t201 / 0.2e1;
t558 = t201 / 0.2e1;
t557 = -t202 / 0.2e1;
t556 = t202 / 0.2e1;
t550 = -t237 / 0.2e1;
t549 = t237 / 0.2e1;
t548 = -t246 / 0.2e1;
t547 = -t247 / 0.2e1;
t546 = t247 / 0.2e1;
t539 = pkin(2) * t333;
t538 = pkin(5) * t202;
t532 = mrSges(5,3) * t246;
t526 = Ifges(4,4) * t339;
t525 = Ifges(4,4) * t343;
t523 = Ifges(5,4) * t338;
t522 = Ifges(5,4) * t342;
t18 = t127 * t342 - t165 * t481 - t197 * t482 - t338 * t65;
t10 = -pkin(4) * t264 - t18;
t517 = t10 * t338;
t516 = t17 * t342;
t513 = t338 * mrSges(5,3);
t284 = t335 * t429 + t406;
t167 = t284 * t343 - t339 * t597;
t509 = t167 * t337;
t285 = -t335 * t428 + t407;
t169 = t285 * t343 - t339 * t596;
t508 = t169 * t337;
t353 = t339 * t408 + t340 * t343;
t213 = t334 * t353 + t335 * t501;
t507 = t213 * t337;
t504 = t284 * t333;
t503 = t285 * t333;
t502 = t333 * t338;
t497 = t337 * t338;
t496 = t337 * t342;
t471 = t333 * t498;
t488 = pkin(2) * t466 + pkin(9) * t471;
t470 = Ifges(5,5) * t162 + Ifges(5,6) * t163 + Ifges(5,3) * t264;
t469 = pkin(3) * t249 + t488;
t468 = Ifges(4,5) * t278 - Ifges(4,6) * t277 + Ifges(4,3) * t323;
t460 = t333 * t484;
t459 = t333 * t483;
t455 = t501 / 0.2e1;
t26 = -mrSges(7,1) * t68 + mrSges(7,2) * t67;
t447 = t481 / 0.2e1;
t442 = -t479 / 0.2e1;
t80 = t172 * t341 - t174 * t337;
t185 = -t260 * t338 + t261 * t342;
t426 = mrSges(4,3) * t462;
t422 = t333 * t461;
t412 = -pkin(2) * t349 + pkin(9) * t504;
t411 = -pkin(2) * t350 + pkin(9) * t503;
t173 = pkin(4) * t499 - t185;
t396 = Ifges(5,1) * t342 - t523;
t391 = -Ifges(5,2) * t338 + t522;
t386 = Ifges(5,5) * t342 - Ifges(5,6) * t338;
t356 = -t333 * t466 + t434;
t171 = t213 * t342 + t338 * t356;
t212 = -t334 * t594 - t335 * t499;
t101 = t171 * t341 + t212 * t337;
t100 = -t171 * t337 + t212 * t341;
t192 = -t284 * t433 - t343 * t349;
t379 = pkin(3) * t192 + t412;
t194 = -t285 * t433 - t343 * t350;
t378 = pkin(3) * t194 + t411;
t354 = t339 * t540 + t340 * t432;
t248 = t354 * t334;
t377 = pkin(10) * t248 + t469;
t107 = -t260 * t481 - t261 * t482 + t272 * t342 - t273 * t338;
t218 = -t287 * t337 - t341 * t499;
t374 = -t287 * t341 + t337 * t499;
t370 = (-mrSges(4,1) * t343 + mrSges(4,2) * t339) * t333;
t369 = (Ifges(4,2) * t343 + t526) * t333;
t216 = -qJD(4) * t286 + t342 * t459;
t217 = qJD(4) * t287 + t338 * t459;
t274 = t291 * qJD(3);
t118 = pkin(4) * t217 - pkin(11) * t216 + t274;
t98 = pkin(11) * t460 + t106;
t19 = t118 * t337 + t172 * t478 - t174 * t480 + t341 * t98;
t363 = t339 * t570 * (Ifges(4,1) * t343 - t526);
t191 = t284 * t432 - t339 * t349;
t360 = pkin(10) * t191 + t379;
t193 = t285 * t432 - t339 * t350;
t359 = pkin(10) * t193 + t378;
t99 = -pkin(4) * t460 - t107;
t20 = -qJD(5) * t81 + t118 * t341 - t337 * t98;
t170 = t213 * t338 - t342 * t356;
t346 = t333 * t350 + t435 * t512;
t345 = t333 * t349 - t436 * t512;
t317 = Ifges(4,4) * t462;
t312 = t336 * t341;
t311 = t336 * t337;
t305 = t467 * t338;
t297 = t341 * t309;
t270 = qJD(2) * t370;
t269 = -mrSges(4,2) * t324 + t426;
t256 = -pkin(10) * t496 + t297;
t235 = Ifges(5,4) * t246;
t226 = qJD(1) * t248;
t222 = Ifges(4,1) * t463 + Ifges(4,5) * t324 + t317;
t221 = Ifges(4,6) * t324 + qJD(2) * t369;
t220 = -qJ(6) * t497 + t257;
t215 = mrSges(4,1) * t323 - mrSges(4,3) * t278;
t214 = -mrSges(4,2) * t323 - mrSges(4,3) * t277;
t211 = -qJ(6) * t495 + t297 + (-pkin(10) * t337 - pkin(5)) * t342;
t205 = -mrSges(5,2) * t310 + t532;
t204 = t249 * t342 + t338 * t471;
t200 = mrSges(4,1) * t277 + mrSges(4,2) * t278;
t168 = t285 * t339 + t343 * t596;
t166 = t284 * t339 + t343 * t597;
t153 = t335 * t459 + (qJD(2) * t355 + qJD(3) * t594) * t334;
t152 = t335 * t460 + (qJD(2) * t354 + qJD(3) * t353) * t334;
t150 = t247 * Ifges(5,1) + t235 + t616;
t148 = Ifges(5,5) * t247 + t246 * Ifges(5,6) + t310 * Ifges(5,3);
t136 = t194 * t342 + t285 * t502;
t134 = t192 * t342 + t284 * t502;
t124 = t196 * t341 + t226 * t337;
t123 = -t196 * t337 + t226 * t341;
t122 = qJD(5) * t374 - t216 * t337 + t341 * t460;
t121 = qJD(5) * t218 + t216 * t341 + t337 * t460;
t117 = -mrSges(5,2) * t264 + mrSges(5,3) * t163;
t108 = -pkin(5) * t218 + t173;
t105 = t169 * t342 + t338 * t346;
t104 = t169 * t338 - t342 * t346;
t103 = t167 * t342 + t338 * t345;
t102 = t167 * t338 - t342 * t345;
t77 = -mrSges(5,1) * t163 + mrSges(5,2) * t162;
t56 = -qJD(4) * t170 + t153 * t342 + t338 * t422;
t52 = qJ(6) * t218 + t81;
t42 = pkin(5) * t286 + qJ(6) * t374 + t80;
t41 = -pkin(5) * t122 + t99;
t38 = -mrSges(6,2) * t158 + mrSges(6,3) * t68;
t37 = -mrSges(7,2) * t158 + mrSges(7,3) * t68;
t36 = mrSges(6,1) * t158 - mrSges(6,3) * t67;
t35 = mrSges(7,1) * t158 - mrSges(7,3) * t67;
t7 = qJ(6) * t122 + qJD(6) * t218 + t19;
t6 = pkin(5) * t217 - qJ(6) * t121 + qJD(6) * t374 + t20;
t5 = -pkin(5) * t68 + qJDD(6) + t10;
t8 = [m(4) * (t153 * t184 + t209 * t356 + t213 * t78 + t234 * t422) + t356 * t200 + m(5) * (t17 * t171 + t56 * t95) + m(3) * (t335 ^ 2 * qJDD(1) + (t279 * t540 + t280 * t340) * t334) + t270 * t422 + t153 * t269 + t213 * t214 + t56 * t205 + t171 * t117 + m(2) * qJDD(1) + (-m(4) * t79 + m(5) * t66 - t215 + t77) * t212 + (m(6) * t32 + m(7) * t25 - t604) * (qJD(5) * t100 + t152 * t337 + t341 * t56) + (m(6) * t31 + m(7) * t23 - t603) * (-qJD(5) * t101 + t152 * t341 - t337 * t56) + t639 * t152 + (m(6) * t3 + m(7) * t2 + t37 + t38) * t101 + (m(6) * t4 + m(7) * t1 + t35 + t36) * t100 + (-qJDD(2) * t498 - t344 * t466) * mrSges(3,2) + (qJDD(2) * t466 - t344 * t498) * mrSges(3,1) - t571 * (qJD(4) * t171 + t153 * t338 - t342 * t422) + (-m(5) * t18 + m(6) * t10 + m(7) * t5 + t26 + t614) * t170 + (-m(2) - m(3) - m(4) + t583) * g(3); (Ifges(4,4) * t278 - Ifges(4,2) * t277 + Ifges(4,6) * t323) * t499 / 0.2e1 + (Ifges(4,1) * t278 - Ifges(4,4) * t277 + Ifges(4,5) * t323) * t455 + (t10 * t173 + t3 * t81 + t4 * t80 + t82 * t99 + (-t124 + t19) * t32 + (-t123 + t20) * t31) * m(6) + (-t183 * t459 - t184 * t460 + t499 * t78 - t501 * t79) * mrSges(4,3) + (mrSges(6,3) * t32 + mrSges(7,3) * t25 - mrSges(6,1) * t82 - mrSges(7,1) * t47 + t618 / 0.2e1 + t624 * t549 + t625 * t558 + t627 * t556) * t122 + t512 * t468 / 0.2e1 - t277 * (Ifges(4,6) * t512 + t369) / 0.2e1 + t323 * (Ifges(4,3) * t512 + (Ifges(4,5) * t339 + Ifges(4,6) * t343) * t333) / 0.2e1 + (t107 * t94 + t17 * t186 + t18 * t185 + t259 * t66 + t602 * t95) * m(5) + (t4 * mrSges(6,3) + t1 * mrSges(7,3) - t10 * mrSges(6,2) - t5 * mrSges(7,2) - t620 / 0.2e1 + t626 * t562 - t627 * t566 - t628 * t567) * t374 + (t419 + t279) * mrSges(3,1) + (t164 * t216 + t17 * t499 + t287 * t66 - t460 * t95) * mrSges(5,2) - t200 * t539 + t310 * (Ifges(5,5) * t216 + Ifges(5,3) * t460) / 0.2e1 + (Ifges(5,5) * t287 - Ifges(5,3) * t499) * t545 - t221 * t460 / 0.2e1 + t94 * (mrSges(5,1) * t460 - mrSges(5,3) * t216) - t270 * t423 + t278 * (Ifges(4,5) * t512 + (Ifges(4,1) * t339 + t525) * t333) / 0.2e1 + t363 * t476 / 0.2e1 + (t570 * qJD(2) * (-Ifges(4,2) * t339 + t525) + t333 * t222) * t483 / 0.2e1 + t571 * (t227 * t338 - t342 * t423) + (t316 - t280) * mrSges(3,2) + (Ifges(5,1) * t287 - Ifges(5,5) * t499) * t561 + (Ifges(5,1) * t216 + Ifges(5,5) * t460) * t546 + (-t647 / 0.2e1 - Ifges(5,4) * t546 + t623 * t549 + t624 * t558 - t626 * t556 + t638 - t651) * t217 + (-mrSges(6,3) * t31 - mrSges(7,3) * t23 + mrSges(6,2) * t82 + mrSges(7,2) * t47 + t617 / 0.2e1 - t626 * t549 + t627 * t558 + t628 * t556) * t121 - t470 * t499 / 0.2e1 + t18 * (-mrSges(5,1) * t499 - mrSges(5,3) * t287) + (-mrSges(6,1) * t10 - mrSges(7,1) * t5 + t3 * mrSges(6,3) + t2 * mrSges(7,3) + t562 * t624 + t566 * t625 + t567 * t627 + t648) * t218 + t209 * t370 + t639 * (t274 - t226) + (t184 * t599 - t209 * t539 - t234 * t423 + t290 * t79 + t291 * t78) * m(4) + (t1 * t42 + t108 * t5 + t2 * t52 + t41 * t47 + (-t124 + t7) * t25 + (-t123 + t6) * t23) * m(7) + t246 * (Ifges(5,4) * t216 + Ifges(5,6) * t460) / 0.2e1 + (Ifges(5,4) * t287 - Ifges(5,6) * t499) * t560 + Ifges(3,3) * qJDD(2) + t290 * t215 + t291 * t214 + t259 * t77 + t216 * t150 / 0.2e1 + t107 * t206 + t186 * t117 + t185 * t116 + t173 * t27 + t6 * t142 + t20 * t143 + t7 * t140 + t19 * t141 + t41 * t112 + t99 * t113 + t108 * t26 + t80 * t36 + t81 * t38 + t52 * t37 + t42 * t35 + t287 * t568 + t512 * t632 + (t148 * t455 - t582) * qJD(3) + t599 * t269 + t602 * t205 + t603 * t123 + t604 * t124 + (-t17 * mrSges(5,3) + t66 * mrSges(5,1) - t658 + t622 / 0.2e1 + t656) * t286 + (-m(7) * (t134 * t332 + t379) - m(4) * t412 - t192 * mrSges(4,1) - mrSges(4,3) * t504 + mrSges(3,1) * t349 + t284 * mrSges(3,2) - m(6) * (pkin(4) * t134 + t360) - m(5) * t360 - t134 * mrSges(5,1) + t631 * (t134 * t341 + t191 * t337) - t629 * (-t134 * t337 + t191 * t341) + t576 * t191 + t572 * (t192 * t338 - t284 * t500)) * g(2) + (-m(7) * (t136 * t332 + t378) - m(4) * t411 - t194 * mrSges(4,1) - mrSges(4,3) * t503 + mrSges(3,1) * t350 + t285 * mrSges(3,2) - m(6) * (pkin(4) * t136 + t359) - m(5) * t359 - t136 * mrSges(5,1) + t631 * (t136 * t341 + t193 * t337) - t629 * (-t136 * t337 + t193 * t341) + t576 * t193 + t572 * (t194 * t338 - t285 * t500)) * g(1) + (-(mrSges(3,1) * t540 - mrSges(3,2) * t340) * t334 - m(7) * (t204 * t332 + t469) - m(4) * t488 - t249 * mrSges(4,1) - mrSges(4,3) * t471 - m(6) * (pkin(4) * t204 + t377) - m(5) * t377 - t204 * mrSges(5,1) + t576 * t248 + t631 * (t204 * t341 + t248 * t337) - t629 * (-t204 * t337 + t248 * t341) + t572 * (t249 * t338 - t342 * t471)) * g(3) - t512 * t633; (-t474 - t129) * t205 + (t516 + (t486 * t492 - t481) * t94) * mrSges(5,3) - (t310 * (Ifges(5,3) * t339 + t343 * t386) + t247 * (Ifges(5,5) * t339 + t343 * t396) + t246 * (Ifges(5,6) * t339 + t343 * t391) + t339 * t148 + (-Ifges(4,2) * t463 + t150 * t342 + t338 * t619 + t222 + t317) * t343) * t486 / 0.2e1 + (t427 + t574) * t184 - t18 * t513 + t614 * pkin(10) * t338 + (-t1 * t495 - t2 * t497) * mrSges(7,3) + t495 * t649 + (t5 * t397 + t562 * t592 + t566 * t590 + t588 * t567 + t568 + t637) * t338 + t638 * t482 + t643 * t113 + t310 * t164 * (mrSges(5,1) * t338 + mrSges(5,2) * t342) + t617 * (t337 * t442 + t341 * t447 - t232 / 0.2e1) + t618 * (t341 * t442 - t458 / 0.2e1 - t231 / 0.2e1) + (-t128 * t94 - t129 * t95 - pkin(3) * t66 + (t516 - t18 * t338 + (-t338 * t95 - t342 * t94) * qJD(4)) * pkin(10)) * m(5) + (-t3 * t497 - t4 * t495) * mrSges(6,3) + (t246 * t391 + t247 * t396 + t310 * t386) * qJD(4) / 0.2e1 + (t426 - t269) * t183 + t66 * t401 + t523 * t560 + (-t95 * (-mrSges(5,2) * t339 - t343 * t513) - t94 * mrSges(5,1) * t339) * t486 - t344 * t363 / 0.2e1 - t633 + t632 + t149 * t421 / 0.2e1 + (-t473 - t128) * t206 + t468 + t522 * t561 + t305 * t26 + t256 * t36 + t257 * t38 + t220 * t37 + t211 * t35 - pkin(3) * t77 + (-t509 * t569 + t631 * (-t166 * t493 + t509) - t629 * (t166 * t496 + t167 * t341) + t583 * (-pkin(3) * t166 + pkin(10) * t167) + t630 * t167 + t577 * t166) * g(2) + (-t508 * t569 + t631 * (-t168 * t493 + t508) - t629 * (t168 * t496 + t169 * t341) + t583 * (-pkin(3) * t168 + pkin(10) * t169) + t630 * t169 + t577 * t168) * g(1) + (-t507 * t569 + t583 * (-pkin(3) * t212 + pkin(10) * t213) + t630 * t213 + t631 * (-t212 * t493 + t507) - t629 * (t212 * t496 + t213 * t341) + t577 * t212) * g(3) + t399 * t517 + (pkin(10) * t117 - t656) * t342 + (t221 * t455 + t582) * qJD(2) + t150 * t447 + (-mrSges(7,1) * t587 + mrSges(7,3) * t600) * t23 + (-mrSges(6,1) * t587 + mrSges(6,3) * t600) * t31 + (mrSges(6,1) * t601 - mrSges(6,2) * t600) * t82 + (mrSges(7,1) * t601 - mrSges(7,2) * t600) * t47 + (mrSges(6,2) * t587 - mrSges(6,3) * t601) * t32 + (mrSges(7,2) * t587 - mrSges(7,3) * t601) * t25 + t609 * t112 + t610 * t143 + t611 * t141 + (t256 * t4 + t257 * t3 + (t481 * t82 + t517) * pkin(10) - t114 * t82 + t611 * t32 + t610 * t31) * m(6) + t612 * t140 + t613 * t142 + (t1 * t211 + t2 * t220 + t23 * t613 + t25 * t612 + t305 * t5 + t47 * t609) * m(7) - t621 * t497 / 0.2e1 - t622 * t342 / 0.2e1 + (-t593 * t479 + (t338 * t623 + t342 * t592) * qJD(4)) * t549 + (-t591 * t479 + (t338 * t624 + t342 * t590) * qJD(4)) * t558 + (-t589 * t479 + (-t338 * t626 + t342 * t588) * qJD(4)) * t556 + (t231 * t624 - t232 * t626 + t421 * t623) * t550 + (t231 * t625 + t232 * t627 + t421 * t624) * t559 + (t231 * t627 + t232 * t628 - t421 * t626) * t557; (-t205 + t532) * t94 + (-t480 / 0.2e1 + t506 / 0.2e1) * t618 + (t531 + t575) * t95 + t341 * t648 + t337 * t649 + (t235 + t150) * t548 + (t201 * t590 + t202 * t588 + t237 * t592) * qJD(5) / 0.2e1 + (t478 / 0.2e1 - t505 / 0.2e1) * t617 + t5 * t398 + t10 * t400 + (-Ifges(5,2) * t548 + t623 * t550 - t626 * t557 + t624 * t559 + t651) * t247 + (-t1 * t337 + t2 * t341 + t23 * t586 - t25 * t653) * mrSges(7,3) + (t31 * t586 - t32 * t653 + t584) * mrSges(6,3) + t470 - t332 * t26 + t311 * t35 - t312 * t37 - t43 * t143 - t44 * t141 - pkin(4) * t27 + t18 * mrSges(5,1) - t17 * mrSges(5,2) + (-pkin(4) * t10 - t31 * t43 - t32 * t44) * m(6) + t149 * t546 + (t104 * t580 + t105 * t572) * g(1) + (t102 * t580 + t103 * t572) * g(2) + (t170 * t580 + t171 * t572) * g(3) + (m(6) * ((-t31 * t341 - t32 * t337) * qJD(5) + t584) - t141 * t480 - t143 * t478 - t337 * t36 + t341 * t38) * pkin(11) + t589 * t567 + t591 * t566 + t593 * t562 + t595 * qJD(5) + t606 * t112 + t607 * t142 + t608 * t140 + (t1 * t311 - t2 * t312 + t23 * t607 + t25 * t608 - t332 * t5 + t47 * t606) * m(7) + (-t616 / 0.2e1 - t164 * mrSges(5,2) + Ifges(5,1) * t547 + t590 * t559 + t588 * t557 + t592 * t550 - t595) * t246 + (-t524 + t619) * t547; (-t629 * (-t103 * t341 - t166 * t337) + t581 * (-t103 * t337 + t166 * t341)) * g(2) + (-t629 * (-t105 * t341 - t168 * t337) + t581 * (-t105 * t337 + t168 * t341)) * g(1) + (-t141 + t530) * t31 + t636 + t605 * t1 - t112 * t538 + (-t202 * t625 + t617 + t654) * t559 + (t100 * t581 + t101 * t629) * g(3) + t618 * t556 + (-m(7) * (-t23 + t24) + t142 + t527) * t25 + (t201 * t628 - t652) * t557 + (-t201 * t626 - t202 * t624) * t550 + (-m(7) * t538 - mrSges(7,1) * t202 - mrSges(7,2) * t201) * t47 + (t143 + t529) * t32 - t82 * (mrSges(6,1) * t202 + mrSges(6,2) * t201) - t24 * t140 + pkin(5) * t35 + t23 * t528 + t622; -t201 * t140 + t202 * t142 + (-g(1) * t104 - g(2) * t102 - g(3) * t170 - t201 * t25 + t202 * t23 + t5) * m(7) + t26;];
tau  = t8;
