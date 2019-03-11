% Calculate vector of inverse dynamics joint torques for
% S6RRPRPR13
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 11:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPR13_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR13_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR13_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR13_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR13_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR13_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR13_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:24:11
% EndTime: 2019-03-09 11:25:16
% DurationCPUTime: 43.31s
% Computational Cost: add. (16258->1009), mult. (38660->1335), div. (0->0), fcn. (29723->14), ass. (0->441)
t366 = cos(qJ(2));
t358 = cos(pkin(6));
t474 = qJD(1) * t358;
t463 = pkin(1) * t474;
t332 = t366 * t463;
t643 = qJD(3) - t332;
t365 = cos(qJ(4));
t361 = sin(qJ(4));
t399 = pkin(4) * t365 + qJ(5) * t361;
t362 = sin(qJ(2));
t356 = sin(pkin(6));
t475 = qJD(1) * t356;
t452 = t362 * t475;
t551 = pkin(3) + pkin(8);
t642 = -(-t399 - t551) * t452 + qJD(4) * t399 - qJD(5) * t365 + t643;
t357 = cos(pkin(11));
t348 = pkin(5) * t357 + pkin(4);
t355 = sin(pkin(11));
t409 = -mrSges(6,1) * t357 + mrSges(6,2) * t355;
t641 = -m(6) * pkin(4) - m(7) * t348 + t409;
t490 = t356 * t362;
t450 = qJD(2) * t490;
t488 = t356 * t366;
t277 = qJD(1) * t450 - qJDD(1) * t488;
t464 = qJDD(1) * t358;
t335 = qJDD(2) + t464;
t337 = qJD(2) + t474;
t451 = t366 * t475;
t422 = t361 * t451;
t467 = qJD(4) * t365;
t152 = -qJD(4) * t422 - t365 * t277 + t335 * t361 + t337 * t467;
t542 = t152 / 0.2e1;
t247 = t337 * t361 + t365 * t451;
t469 = qJD(4) * t247;
t151 = t277 * t361 + t335 * t365 - t469;
t471 = qJD(2) * t366;
t278 = (qJD(1) * t471 + qJDD(1) * t362) * t356;
t258 = qJDD(4) + t278;
t107 = t151 * t357 + t258 * t355;
t547 = t107 / 0.2e1;
t106 = -t151 * t355 + t258 * t357;
t548 = t106 / 0.2e1;
t640 = Ifges(6,5) * t547 + Ifges(6,6) * t548 + Ifges(6,3) * t542;
t327 = pkin(2) * t452;
t398 = pkin(9) * t362 - qJ(3) * t366;
t230 = t398 * t475 + t327;
t330 = pkin(8) * t451;
t274 = t362 * t463 + t330;
t232 = pkin(3) * t451 + t274;
t136 = t365 * t230 + t361 * t232;
t124 = qJ(5) * t451 + t136;
t552 = pkin(2) + pkin(9);
t466 = qJD(4) * t552;
t446 = t365 * t466;
t601 = (-t124 - t446) * t357 + t642 * t355;
t639 = t124 * t355 + t357 * t642;
t354 = pkin(11) + qJ(6);
t350 = sin(t354);
t351 = cos(t354);
t580 = -mrSges(7,1) * t351 + mrSges(7,2) * t350;
t587 = -t580 - t641;
t635 = -m(5) - m(7);
t465 = m(6) - t635;
t631 = -mrSges(3,1) + mrSges(4,2);
t638 = -pkin(9) * t465 - mrSges(5,3) + t631;
t318 = t337 * qJ(3);
t189 = t318 + t232;
t312 = qJD(4) + t452;
t607 = t312 * Ifges(5,5);
t637 = t189 * mrSges(5,2) + t607 / 0.2e1;
t248 = t337 * t365 - t422;
t180 = t248 * t357 + t312 * t355;
t109 = pkin(4) * t247 - qJ(5) * t248 + t189;
t169 = -t337 * t552 + t452 * t551 + t643;
t437 = -qJ(3) * t362 - pkin(1);
t197 = (-t366 * t552 + t437) * t475;
t104 = t169 * t361 + t197 * t365;
t93 = qJ(5) * t312 + t104;
t52 = t357 * t109 - t355 * t93;
t36 = pkin(5) * t247 - pkin(10) * t180 + t52;
t360 = sin(qJ(6));
t364 = cos(qJ(6));
t432 = -t248 * t355 + t357 * t312;
t53 = t355 * t109 + t357 * t93;
t42 = pkin(10) * t432 + t53;
t13 = t36 * t364 - t360 * t42;
t14 = t36 * t360 + t364 * t42;
t508 = t104 * mrSges(5,3);
t606 = t312 * Ifges(5,6);
t636 = -t508 - t606 / 0.2e1 + t189 * mrSges(5,1) + t52 * mrSges(6,1) + t13 * mrSges(7,1) - t53 * mrSges(6,2) - t14 * mrSges(7,2);
t143 = qJDD(6) + t152;
t614 = -t180 * t360 + t364 * t432;
t34 = qJD(6) * t614 + t106 * t360 + t107 * t364;
t100 = t180 * t364 + t360 * t432;
t35 = -qJD(6) * t100 + t106 * t364 - t107 * t360;
t6 = Ifges(7,5) * t34 + Ifges(7,6) * t35 + Ifges(7,3) * t143;
t633 = t6 / 0.2e1 + t640;
t630 = Ifges(6,6) * t432;
t103 = t169 * t365 - t361 * t197;
t629 = t103 * mrSges(5,3);
t628 = t180 * Ifges(6,5);
t484 = t361 * t362;
t238 = (t355 * t366 + t357 * t484) * t475;
t423 = t365 * t452;
t435 = t355 * t552 + pkin(5);
t486 = t357 * t361;
t626 = pkin(5) * t423 + pkin(10) * t238 + (pkin(10) * t486 + t365 * t435) * qJD(4) + t639;
t237 = (-t355 * t484 + t357 * t366) * t475;
t468 = qJD(4) * t361;
t448 = t355 * t468;
t625 = t601 + (-t237 + t448) * pkin(10);
t407 = t350 * mrSges(7,1) + t351 * mrSges(7,2);
t499 = t357 * mrSges(6,2);
t408 = t355 * mrSges(6,1) + t499;
t525 = pkin(5) * t355;
t624 = m(7) * t525 + t407 + t408;
t520 = pkin(10) + qJ(5);
t623 = m(7) * t520 + mrSges(6,3) + mrSges(7,3);
t338 = pkin(8) * t490;
t527 = pkin(1) * t358;
t462 = qJD(2) * t527;
t426 = qJD(1) * t462;
t457 = pkin(1) * t464;
t182 = -qJD(2) * t330 - qJDD(1) * t338 - t362 * t426 + t366 * t457;
t376 = qJDD(3) - t182;
t115 = pkin(3) * t278 - t335 * t552 + t376;
t470 = qJD(3) * t362;
t372 = -qJ(3) * t278 + (-pkin(1) * qJDD(1) - qJD(1) * t470) * t356;
t118 = t277 * t552 + t372;
t37 = t361 * t115 + t365 * t118 + t169 * t467 - t197 * t468;
t27 = qJ(5) * t258 + qJD(5) * t312 + t37;
t181 = -pkin(8) * t277 + t362 * t457 + t366 * t426;
t146 = -t335 * qJ(3) - t337 * qJD(3) - t181;
t117 = -pkin(3) * t277 - t146;
t45 = pkin(4) * t152 - qJ(5) * t151 - qJD(5) * t248 + t117;
t11 = -t27 * t355 + t357 * t45;
t5 = pkin(5) * t152 - pkin(10) * t107 + t11;
t12 = t357 * t27 + t355 * t45;
t9 = pkin(10) * t106 + t12;
t1 = qJD(6) * t13 + t360 * t5 + t364 * t9;
t2 = -qJD(6) * t14 - t360 * t9 + t364 * t5;
t622 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t532 = t258 / 0.2e1;
t543 = -t152 / 0.2e1;
t544 = t151 / 0.2e1;
t558 = Ifges(5,1) * t544 + Ifges(5,4) * t543 + Ifges(5,5) * t532;
t504 = t248 * Ifges(5,4);
t139 = -t247 * Ifges(5,2) + t504 + t606;
t244 = qJD(6) + t247;
t537 = t244 / 0.2e1;
t549 = t100 / 0.2e1;
t553 = t614 / 0.2e1;
t608 = t100 * Ifges(7,5) + Ifges(7,6) * t614 + t247 * Ifges(6,3) + t244 * Ifges(7,3) + t628 + t630;
t621 = Ifges(7,5) * t549 + Ifges(7,6) * t553 + Ifges(7,3) * t537 - t139 / 0.2e1 + t608 / 0.2e1;
t459 = m(4) + t465;
t616 = pkin(2) * t459 + t624 - t638;
t410 = t361 * mrSges(5,1) + t365 * mrSges(5,2);
t615 = -t361 * t587 - t410;
t567 = t34 / 0.2e1;
t566 = t35 / 0.2e1;
t545 = t143 / 0.2e1;
t522 = mrSges(4,3) - mrSges(3,2);
t612 = Ifges(3,5) - Ifges(4,4);
t611 = Ifges(4,5) - Ifges(3,6);
t433 = -qJ(5) * t365 + qJ(3);
t302 = pkin(4) * t361 + t433;
t296 = t357 * t302;
t485 = t357 * t365;
t192 = -pkin(10) * t485 + t361 * t435 + t296;
t483 = t361 * t552;
t241 = t355 * t302 - t357 * t483;
t491 = t355 * t365;
t201 = -pkin(10) * t491 + t241;
t119 = t192 * t364 - t201 * t360;
t610 = qJD(6) * t119 + t360 * t626 + t364 * t625;
t120 = t192 * t360 + t201 * t364;
t609 = -qJD(6) * t120 - t360 * t625 + t364 * t626;
t111 = mrSges(5,1) * t258 - mrSges(5,3) * t151;
t58 = -t106 * mrSges(6,1) + t107 * mrSges(6,2);
t605 = t111 - t58;
t306 = t520 * t355;
t307 = t520 * t357;
t221 = -t306 * t364 - t307 * t360;
t393 = t355 * t360 - t357 * t364;
t492 = t247 * t357;
t164 = pkin(4) * t248 + qJ(5) * t247;
t71 = -t103 * t355 + t357 * t164;
t51 = pkin(5) * t248 + pkin(10) * t492 + t71;
t493 = t247 * t355;
t72 = t357 * t103 + t355 * t164;
t59 = pkin(10) * t493 + t72;
t604 = -qJD(5) * t393 + qJD(6) * t221 - t360 * t51 - t364 * t59;
t222 = -t306 * t360 + t307 * t364;
t300 = t355 * t364 + t357 * t360;
t603 = -qJD(5) * t300 - qJD(6) * t222 + t360 * t59 - t364 * t51;
t602 = t355 * t446 + t639;
t135 = -t361 * t230 + t232 * t365;
t125 = -pkin(4) * t451 - t135;
t434 = t552 + t525;
t600 = pkin(5) * t237 - t434 * t468 - t125;
t599 = t189 * (mrSges(5,1) * t365 - mrSges(5,2) * t361);
t147 = t237 * t364 - t238 * t360;
t579 = qJD(6) * t393;
t175 = t300 * t468 + t365 * t579;
t596 = t147 - t175;
t148 = t237 * t360 + t238 * t364;
t285 = t300 * qJD(6);
t173 = -t285 * t365 + t393 * t468;
t595 = t148 - t173;
t149 = t300 * t247;
t594 = t149 + t285;
t150 = t393 * t247;
t593 = t150 + t579;
t235 = t337 * t357 + t355 * t423;
t236 = t337 * t355 - t357 * t423;
t265 = t393 * t365;
t592 = -qJD(4) * t265 - t235 * t360 - t236 * t364 - t285 * t361;
t263 = t300 * t365;
t591 = -qJD(4) * t263 - t235 * t364 + t236 * t360 + t361 * t579;
t108 = -mrSges(6,1) * t432 + mrSges(6,2) * t180;
t191 = mrSges(5,1) * t312 - mrSges(5,3) * t248;
t590 = t191 - t108;
t165 = mrSges(5,1) * t247 + mrSges(5,2) * t248;
t429 = mrSges(4,1) * t451;
t269 = -mrSges(4,3) * t337 - t429;
t589 = -t269 + t165;
t428 = mrSges(3,3) * t452;
t430 = mrSges(4,1) * t452;
t588 = t337 * t631 + t428 + t430;
t92 = -pkin(4) * t312 + qJD(5) - t103;
t586 = -t92 * t408 + t629;
t585 = t423 + t467;
t584 = -m(6) * qJ(5) - t623;
t132 = -mrSges(6,2) * t247 + mrSges(6,3) * t432;
t133 = mrSges(6,1) * t247 - mrSges(6,3) * t180;
t583 = t132 * t357 - t133 * t355;
t38 = t115 * t365 - t361 * t118 - t169 * t468 - t197 * t467;
t582 = -t361 * t37 - t365 * t38;
t69 = -mrSges(6,2) * t152 + mrSges(6,3) * t106;
t70 = mrSges(6,1) * t152 - mrSges(6,3) * t107;
t581 = -t355 * t70 + t357 * t69;
t397 = -t11 * t355 + t12 * t357;
t517 = Ifges(3,4) * t362;
t578 = pkin(1) * (mrSges(3,1) * t362 + mrSges(3,2) * t366) - t362 * (Ifges(3,1) * t366 - t517) / 0.2e1;
t576 = mrSges(5,1) + t587;
t379 = mrSges(5,2) + t584;
t575 = qJ(3) * t459 + t522;
t574 = t38 * mrSges(5,1) - t37 * mrSges(5,2) + Ifges(5,5) * t151 - Ifges(5,6) * t152 + Ifges(5,3) * t258;
t572 = -t355 * (m(7) * pkin(5) + mrSges(6,1)) - t499 + t638;
t570 = -m(6) * t433 - t522 + t615 + t623 * t365 + (-m(4) + t635) * qJ(3);
t569 = Ifges(7,4) * t567 + Ifges(7,2) * t566 + Ifges(7,6) * t545;
t568 = Ifges(7,1) * t567 + Ifges(7,4) * t566 + Ifges(7,5) * t545;
t40 = t107 * Ifges(6,4) + t106 * Ifges(6,2) + t152 * Ifges(6,6);
t565 = t40 / 0.2e1;
t564 = Ifges(6,1) * t547 + Ifges(6,4) * t548 + Ifges(6,5) * t542;
t512 = Ifges(7,4) * t100;
t47 = Ifges(7,2) * t614 + Ifges(7,6) * t244 + t512;
t563 = -t47 / 0.2e1;
t562 = t47 / 0.2e1;
t96 = Ifges(7,4) * t614;
t48 = Ifges(7,1) * t100 + Ifges(7,5) * t244 + t96;
t561 = -t48 / 0.2e1;
t560 = t48 / 0.2e1;
t559 = -t151 * Ifges(5,4) / 0.2e1 + Ifges(5,2) * t542 - t258 * Ifges(5,6) / 0.2e1;
t88 = t180 * Ifges(6,4) + Ifges(6,2) * t432 + Ifges(6,6) * t247;
t557 = -t88 / 0.2e1;
t556 = t88 / 0.2e1;
t89 = t180 * Ifges(6,1) + Ifges(6,4) * t432 + Ifges(6,5) * t247;
t555 = t89 / 0.2e1;
t554 = -t614 / 0.2e1;
t550 = -t100 / 0.2e1;
t541 = -t432 / 0.2e1;
t540 = -t180 / 0.2e1;
t538 = -t244 / 0.2e1;
t536 = -t247 / 0.2e1;
t535 = t247 / 0.2e1;
t533 = t248 / 0.2e1;
t528 = pkin(1) * t356;
t526 = pkin(1) * t366;
t521 = Ifges(3,4) + Ifges(4,6);
t329 = pkin(2) * t450;
t194 = t329 + (qJD(2) * t398 - t470) * t356;
t453 = -pkin(2) - t526;
t200 = pkin(3) * t490 + t338 + (-pkin(9) + t453) * t358;
t477 = pkin(2) * t488 + qJ(3) * t490;
t252 = -t477 - t528;
t341 = pkin(9) * t488;
t224 = t252 - t341;
t347 = t362 * t527;
t233 = (t488 * t551 + t347) * qJD(2);
t76 = t365 * t194 + t200 * t467 - t224 * t468 + t361 * t233;
t65 = (qJ(5) * t471 + qJD(5) * t362) * t356 + t76;
t333 = t366 * t462;
t349 = t358 * qJD(3);
t431 = t551 * t490;
t199 = -qJD(2) * t431 + t333 + t349;
t286 = t358 * t361 + t365 * t488;
t212 = -qJD(4) * t286 + t361 * t450;
t456 = t361 * t488;
t213 = -qJD(4) * t456 + t358 * t467 - t365 * t450;
t287 = t358 * t365 - t456;
t86 = pkin(4) * t213 - qJ(5) * t212 - qJD(5) * t287 + t199;
t29 = t355 * t86 + t357 * t65;
t516 = Ifges(5,4) * t361;
t515 = Ifges(5,4) * t365;
t514 = Ifges(6,4) * t355;
t513 = Ifges(6,4) * t357;
t511 = Ifges(4,6) * t362;
t510 = Ifges(4,6) * t366;
t505 = t247 * Ifges(5,6);
t503 = t248 * Ifges(5,5);
t30 = -pkin(4) * t258 + qJDD(5) - t38;
t502 = t30 * t365;
t501 = t312 * Ifges(5,3);
t363 = sin(qJ(1));
t489 = t356 * t363;
t367 = cos(qJ(1));
t487 = t356 * t367;
t482 = t362 * t363;
t481 = t362 * t367;
t480 = t363 * t366;
t478 = t366 * t367;
t131 = t361 * t200 + t365 * t224;
t121 = qJ(5) * t490 + t131;
t294 = pkin(8) * t488 + t347;
t251 = -t358 * qJ(3) - t294;
t223 = pkin(3) * t488 - t251;
t134 = pkin(4) * t286 - qJ(5) * t287 + t223;
t68 = t357 * t121 + t355 * t134;
t476 = t367 * pkin(1) + pkin(8) * t489;
t472 = qJD(1) ^ 2 * t356 ^ 2;
t461 = Ifges(3,5) / 0.2e1 - Ifges(4,4) / 0.2e1;
t460 = -Ifges(3,6) / 0.2e1 + Ifges(4,5) / 0.2e1;
t291 = -t358 * t482 + t478;
t455 = t291 * pkin(2) + t476;
t449 = t356 * t471;
t447 = t361 * t466;
t10 = -t35 * mrSges(7,1) + t34 * mrSges(7,2);
t436 = -pkin(1) * t363 + pkin(8) * t487;
t28 = -t355 * t65 + t357 * t86;
t211 = t278 * mrSges(4,1) + t335 * mrSges(4,2);
t67 = -t121 * t355 + t357 * t134;
t130 = t200 * t365 - t361 * t224;
t427 = mrSges(3,3) * t451;
t424 = t361 * t452;
t289 = t358 * t481 + t480;
t416 = -t289 * pkin(2) + t436;
t54 = -mrSges(7,1) * t614 + mrSges(7,2) * t100;
t415 = (t54 - t590) * t361;
t412 = mrSges(5,1) * t286 + mrSges(5,2) * t287;
t406 = mrSges(4,2) * t366 - mrSges(4,3) * t362;
t405 = Ifges(5,1) * t361 + t515;
t404 = Ifges(6,1) * t357 - t514;
t403 = Ifges(5,2) * t365 + t516;
t402 = -Ifges(6,2) * t355 + t513;
t401 = Ifges(5,5) * t361 + Ifges(5,6) * t365;
t400 = Ifges(6,5) * t357 - Ifges(6,6) * t355;
t396 = -t355 * t52 + t357 * t53;
t209 = t287 * t357 + t355 * t490;
t49 = pkin(5) * t286 - pkin(10) * t209 + t67;
t208 = -t287 * t355 + t357 * t490;
t55 = pkin(10) * t208 + t68;
t15 = -t360 * t55 + t364 * t49;
t16 = t360 * t49 + t364 * t55;
t394 = t103 * t361 - t104 * t365;
t126 = t208 * t364 - t209 * t360;
t127 = t208 * t360 + t209 * t364;
t273 = pkin(8) * t452 - t332;
t275 = -pkin(8) * t450 + t333;
t77 = -t361 * t194 - t200 * t468 - t224 * t467 + t233 * t365;
t288 = -t358 * t478 + t482;
t218 = -t288 * t361 + t365 * t487;
t216 = t288 * t365 + t361 * t487;
t387 = t362 * (-Ifges(4,2) * t366 + t511);
t386 = t366 * (Ifges(4,3) * t362 - t510);
t123 = -pkin(4) * t490 - t130;
t290 = t358 * t480 + t481;
t214 = -t290 * t365 + t361 * t489;
t382 = -g(1) * t214 + g(2) * t216 - g(3) * t286;
t163 = -pkin(2) * t335 + t376;
t374 = t182 * mrSges(3,1) - t181 * mrSges(3,2) + t163 * mrSges(4,2) - t146 * mrSges(4,3);
t373 = -mrSges(5,1) + t641;
t66 = -pkin(4) * t449 - t77;
t371 = -qJD(4) * t394 - t582;
t325 = Ifges(3,4) * t451;
t317 = Ifges(4,1) * t335;
t316 = Ifges(3,3) * t335;
t297 = t434 * t365;
t293 = t358 * t526 - t338;
t292 = (-mrSges(3,1) * t366 + mrSges(3,2) * t362) * t356;
t276 = t294 * qJD(2);
t272 = -qJ(3) * t451 + t327;
t271 = t406 * t475;
t268 = -mrSges(3,2) * t337 + t427;
t264 = t393 * t361;
t262 = t300 * t361;
t257 = Ifges(4,4) * t278;
t256 = Ifges(3,5) * t278;
t255 = Ifges(4,5) * t277;
t254 = Ifges(3,6) * t277;
t253 = t358 * t453 + t338;
t243 = -t275 - t349;
t242 = Ifges(5,4) * t247;
t240 = t355 * t483 + t296;
t239 = (-pkin(2) * t366 + t437) * t475;
t234 = t329 + (-qJ(3) * t471 - t470) * t356;
t231 = -qJD(1) * t431 + t332;
t229 = -t318 - t274;
t228 = t337 * Ifges(4,4) + (-t362 * Ifges(4,2) - t510) * t475;
t227 = t337 * Ifges(4,5) + (-t366 * Ifges(4,3) - t511) * t475;
t226 = Ifges(3,1) * t452 + t337 * Ifges(3,5) + t325;
t225 = t337 * Ifges(3,6) + (t366 * Ifges(3,2) + t517) * t475;
t220 = -pkin(2) * t337 + qJD(3) + t273;
t215 = t290 * t361 + t365 * t489;
t210 = mrSges(4,1) * t277 - mrSges(4,3) * t335;
t190 = -mrSges(5,2) * t312 - mrSges(5,3) * t247;
t171 = t212 * t357 + t355 * t449;
t170 = -t212 * t355 + t357 * t449;
t155 = pkin(2) * t277 + t372;
t154 = t215 * t351 + t291 * t350;
t153 = -t215 * t350 + t291 * t351;
t140 = t248 * Ifges(5,1) - t242 + t607;
t138 = t501 + t503 - t505;
t112 = -mrSges(5,2) * t258 - mrSges(5,3) * t152;
t90 = -pkin(5) * t208 + t123;
t85 = -pkin(5) * t493 + t104;
t84 = mrSges(7,1) * t244 - mrSges(7,3) * t100;
t83 = -mrSges(7,2) * t244 + mrSges(7,3) * t614;
t80 = mrSges(5,1) * t152 + mrSges(5,2) * t151;
t73 = -pkin(5) * t432 + t92;
t61 = -qJD(6) * t127 + t170 * t364 - t171 * t360;
t60 = qJD(6) * t126 + t170 * t360 + t171 * t364;
t50 = -pkin(5) * t170 + t66;
t25 = pkin(10) * t170 + t29;
t24 = -mrSges(7,2) * t143 + mrSges(7,3) * t35;
t23 = mrSges(7,1) * t143 - mrSges(7,3) * t34;
t20 = pkin(5) * t213 - pkin(10) * t171 + t28;
t19 = -pkin(5) * t106 + t30;
t4 = -qJD(6) * t16 + t20 * t364 - t25 * t360;
t3 = qJD(6) * t15 + t20 * t360 + t25 * t364;
t7 = [(t11 * mrSges(6,1) - t12 * mrSges(6,2) - t37 * mrSges(5,3) - Ifges(5,4) * t544 + Ifges(7,5) * t567 - Ifges(5,2) * t543 - Ifges(5,6) * t532 + Ifges(7,6) * t566 + Ifges(7,3) * t545 + t559 + t622 + t633 + t640) * t286 + (Ifges(7,1) * t60 + Ifges(7,4) * t61) * t549 + (Ifges(6,1) * t209 + Ifges(6,4) * t208) * t547 + (Ifges(6,5) * t171 + Ifges(6,6) * t170) * t535 + (Ifges(7,5) * t60 + Ifges(7,6) * t61) * t537 + (-m(3) * t476 - m(4) * t455 - t154 * mrSges(7,1) - t153 * mrSges(7,2) - mrSges(2,1) * t367 + mrSges(2,2) * t363 - t575 * t290 + t373 * t215 + t379 * t214 + t572 * t291 - t465 * (pkin(3) * t489 + t455)) * g(2) + (-m(3) * t436 - m(4) * t416 + mrSges(2,1) * t363 + mrSges(2,2) * t367 + t575 * t288 + t379 * t216 + (t373 + t580) * t218 + (t407 - t572) * t289 + t465 * (-pkin(3) * t487 - t416)) * g(1) + (Ifges(7,4) * t60 + Ifges(7,2) * t61) * t553 + t171 * t555 + t170 * t556 + (Ifges(6,4) * t209 + Ifges(6,2) * t208) * t548 + t432 * (Ifges(6,4) * t171 + Ifges(6,2) * t170) / 0.2e1 + (t1 * t126 - t127 * t2 - t13 * t60 + t14 * t61) * mrSges(7,3) + (Ifges(7,5) * t127 + Ifges(7,6) * t126) * t545 + t117 * t412 + (Ifges(7,4) * t127 + Ifges(7,2) * t126) * t566 + (t317 / 0.2e1 - t257 / 0.2e1 + t255 / 0.2e1 + t256 / 0.2e1 - t254 / 0.2e1 + t316 / 0.2e1 + (Ifges(4,1) / 0.2e1 + Ifges(3,3) / 0.2e1) * t335 + t461 * t278 + t460 * t277 + t374) * t358 + t180 * (Ifges(6,1) * t171 + Ifges(6,4) * t170) / 0.2e1 + m(3) * (t181 * t294 + t182 * t293 + t273 * t276 + t274 * t275) + (Ifges(6,5) * t209 + Ifges(6,6) * t208) * t542 + t60 * t560 + t61 * t562 + t209 * t564 + t208 * t565 + t127 * t568 + t126 * t569 + (-t11 * t209 + t12 * t208 + t170 * t53 - t171 * t52) * mrSges(6,3) + m(4) * (t146 * t251 + t155 * t252 + t163 * t253 + t220 * t276 + t229 * t243 + t234 * t239) + m(5) * (t103 * t77 + t104 * t76 + t117 * t223 + t130 * t38 + t131 * t37 + t189 * t199) + m(7) * (t1 * t16 + t13 * t4 + t14 * t3 + t15 * t2 + t19 * t90 + t50 * t73) + m(6) * (t11 * t67 + t12 * t68 + t123 * t30 + t28 * t52 + t29 * t53 + t66 * t92) + (-mrSges(5,3) * t38 + 0.2e1 * t558) * t287 + t293 * (mrSges(3,1) * t335 - mrSges(3,3) * t278) + t294 * (-mrSges(3,2) * t335 - mrSges(3,3) * t277) + Ifges(2,3) * qJDD(1) + t252 * (-mrSges(4,2) * t277 - mrSges(4,3) * t278) + t275 * t268 + t243 * t269 + t234 * t271 + t251 * t210 + t253 * t211 + t223 * t80 + t30 * (-mrSges(6,1) * t208 + mrSges(6,2) * t209) + t199 * t165 + t76 * t190 + t77 * t191 + t92 * (-mrSges(6,1) * t170 + mrSges(6,2) * t171) + t19 * (-mrSges(7,1) * t126 + mrSges(7,2) * t127) + t130 * t111 + t131 * t112 + t29 * t132 + t28 * t133 + t123 * t58 + t66 * t108 + t90 * t10 + t3 * t83 + t4 * t84 + t73 * (-mrSges(7,1) * t61 + mrSges(7,2) * t60) + t68 * t69 + t67 * t70 + t50 * t54 + (Ifges(5,1) * t533 + Ifges(5,4) * t536 + t140 / 0.2e1 - t629 + t637) * t212 + (-Ifges(5,4) * t533 + Ifges(6,3) * t535 - Ifges(5,2) * t536 + t630 / 0.2e1 + t628 / 0.2e1 + t621 + t636) * t213 + t15 * t23 + t16 * t24 + t588 * t276 + ((-mrSges(3,1) * t277 - mrSges(3,2) * t278 + (m(3) * t528 - t292) * qJDD(1)) * pkin(1) + (-t146 * mrSges(4,1) + t155 * mrSges(4,2) + t181 * mrSges(3,3) - t611 * t335 + t521 * t278 + (-Ifges(3,2) - Ifges(4,3)) * t277) * t366 + (-t182 * mrSges(3,3) + t163 * mrSges(4,1) - t155 * mrSges(4,3) + t612 * t335 + (Ifges(3,1) + Ifges(4,2)) * t278 - t521 * t277 + t574) * t362 + ((-t239 * mrSges(4,2) - t225 / 0.2e1 + t227 / 0.2e1 - t274 * mrSges(3,3) + t229 * mrSges(4,1) + t460 * t337) * t362 + (-t239 * mrSges(4,3) + t226 / 0.2e1 - t228 / 0.2e1 + t138 / 0.2e1 - t104 * mrSges(5,2) + t103 * mrSges(5,1) + t220 * mrSges(4,1) + t273 * mrSges(3,3) + t501 / 0.2e1 - t505 / 0.2e1 + t503 / 0.2e1 + t461 * t337) * t366 + (-t387 / 0.2e1 - t386 / 0.2e1 + t366 * (Ifges(3,4) * t366 - Ifges(3,2) * t362) / 0.2e1 - t578) * t475) * qJD(2) + (mrSges(3,3) + mrSges(4,1)) * (-g(1) * t367 - g(2) * t363)) * t356 + (Ifges(7,1) * t127 + Ifges(7,4) * t126) * t567; t578 * t472 + ((t365 * t608 + t225) * t362 + t247 * (Ifges(5,6) * t366 + t362 * t403) + t366 * t228) * t475 / 0.2e1 + (-pkin(2) * t163 - qJ(3) * t146 - qJD(3) * t229 - t239 * t272) * m(4) + t361 * t633 + (-m(4) * t477 + t292 - t465 * (t341 + t477) + (t406 + (-mrSges(5,3) - t624) * t366 + (-t365 * t584 + t615) * t362) * t356) * g(3) + t582 * mrSges(5,3) + t586 * t468 + (-t210 + t80) * qJ(3) - t53 * (mrSges(6,2) * t423 + mrSges(6,3) * t237) + t256 - t257 - t254 + t255 + (Ifges(5,1) * t365 - t516) * t544 + (Ifges(6,5) * t361 + t365 * t404) * t547 + (Ifges(6,6) * t361 + t365 * t402) * t548 + (Ifges(7,1) * t148 + Ifges(7,4) * t147 - Ifges(7,5) * t423) * t550 + (Ifges(7,4) * t148 + Ifges(7,2) * t147 - Ifges(7,6) * t423) * t554 + t448 * t556 + t237 * t557 + t365 * t558 + t361 * t559 + t11 * (mrSges(6,1) * t361 - mrSges(6,3) * t485) + t117 * t410 + (Ifges(5,5) * t365 - Ifges(5,6) * t361) * t532 + (Ifges(6,5) * t238 + Ifges(6,6) * t237 - Ifges(6,3) * t423) * t536 + (Ifges(7,5) * t148 + Ifges(7,6) * t147 - Ifges(7,3) * t423) * t538 + (Ifges(6,1) * t238 + Ifges(6,4) * t237 - Ifges(6,5) * t423) * t540 + (Ifges(6,4) * t238 + Ifges(6,2) * t237 - Ifges(6,6) * t423) * t541 + (Ifges(6,3) * t361 + t365 * t400) * t542 + (-Ifges(5,2) * t361 + t515) * t543 - t605 * t365 * t552 + (Ifges(7,5) * t173 + Ifges(7,6) * t175) * t537 - t52 * (-mrSges(6,1) * t423 - mrSges(6,3) * t238) + (Ifges(7,4) * t173 + Ifges(7,2) * t175) * t553 + t316 + t317 + t408 * t502 - t220 * t429 - t229 * t430 + t600 * t54 + t601 * t132 + t602 * t133 + (-m(4) * t229 + t268 - t269 - t427) * t273 + (t53 * (mrSges(6,3) * t355 * t361 - mrSges(6,2) * t365) + t52 * (mrSges(6,1) * t365 + mrSges(6,3) * t486) + t599) * qJD(4) + (-t508 + t621) * t467 - t40 * t491 / 0.2e1 + t12 * (-mrSges(6,2) * t361 - mrSges(6,3) * t491) + t173 * t560 + t148 * t561 + t175 * t562 + t147 * t563 + t485 * t564 - t265 * t568 - t263 * t569 + (-t103 * (mrSges(5,1) * t366 - mrSges(5,3) * t484) - t239 * (-mrSges(4,2) * t362 - mrSges(4,3) * t366) - t104 * (mrSges(5,3) * t362 * t365 - mrSges(5,2) * t366)) * t475 + t374 - ((t362 * t611 + t366 * t612) * t337 + (t365 * t139 + t361 * t140 + t227) * t362 + (-Ifges(3,2) * t452 + t138 + t226 + t325) * t366 + (Ifges(5,3) * t366 + t362 * t401) * t312 + (Ifges(5,5) * t366 + t362 * t405) * t248) * t475 / 0.2e1 - (t248 * t405 + t312 * t401) * qJD(4) / 0.2e1 + (Ifges(6,3) * t365 - t361 * t400 + t403) * t469 / 0.2e1 + (t11 * t240 + t12 * t241 - (t468 * t92 - t502) * t552 - t125 * t92 + t601 * t53 + t602 * t52) * m(6) + (t117 * qJ(3) - t103 * t135 - t104 * t136 - t552 * t371 + (qJD(3) - t231) * t189) * m(5) + (-t446 - t136) * t190 + (-t447 - t125) * t108 + (t447 - t135) * t191 - t112 * t483 + t297 * t10 + t2 * (mrSges(7,1) * t361 + mrSges(7,3) * t265) + (-Ifges(7,5) * t265 - Ifges(7,6) * t263 + Ifges(7,3) * t361) * t545 + (-Ifges(7,4) * t265 - Ifges(7,2) * t263 + Ifges(7,6) * t361) * t566 + (-Ifges(7,1) * t265 - Ifges(7,4) * t263 + Ifges(7,5) * t361) * t567 + t19 * (mrSges(7,1) * t263 - mrSges(7,2) * t265) + t1 * (-mrSges(7,2) * t361 - mrSges(7,3) * t263) - (t357 * t89 + t140) * t468 / 0.2e1 + (t387 + t386) * t472 / 0.2e1 + (t180 * (Ifges(6,5) * t365 - t361 * t404) + t432 * (Ifges(6,6) * t365 - t361 * t402)) * qJD(4) / 0.2e1 - t272 * t271 + t240 * t70 + t241 * t69 - t92 * (-mrSges(6,1) * t237 + mrSges(6,2) * t238) - t238 * t89 / 0.2e1 - t231 * t165 - pkin(2) * t211 + (Ifges(7,1) * t173 + Ifges(7,4) * t175) * t549 + t120 * t24 + t119 * t23 + (-m(4) * t220 + t428 - t588) * t274 + t589 * qJD(3) + (mrSges(7,1) * t585 + mrSges(7,3) * t595) * t13 + (mrSges(7,1) * t596 - mrSges(7,2) * t595) * t73 + (-mrSges(7,2) * t585 - mrSges(7,3) * t596) * t14 + t609 * t84 + (t1 * t120 + t119 * t2 + t13 * t609 + t14 * t610 + t19 * t297 + t600 * t73) * m(7) + t610 * t83 + (t288 * t616 + t289 * t570) * g(2) + (t290 * t616 + t291 * t570) * g(1) + t452 * t599; t591 * t84 + t592 * t83 + t211 + (t112 + t581) * t361 + (-t10 + t605) * t365 - t589 * t337 + ((t190 + t583) * t365 + t415) * qJD(4) + (t365 * t190 + t271 + t415) * t452 - t264 * t24 - t262 * t23 - t235 * t133 - t236 * t132 + (-g(1) * t290 - g(2) * t288 + g(3) * t488) * t459 + (-t1 * t264 - t19 * t365 - t2 * t262 + (t424 + t468) * t73 + t592 * t14 + t591 * t13) * m(7) + (-t502 + t397 * t361 + (t361 * t92 + t365 * t396) * qJD(4) - t235 * t52 - t236 * t53 + t92 * t424) * m(6) + (-t189 * t337 - t394 * t452 + t371) * m(5) + (t229 * t337 + t239 * t452 + t163) * m(4); (t214 * t576 + t215 * t379) * g(1) + (-t216 * t576 - t218 * t379) * g(2) + t574 + (-t492 * t52 - t493 * t53 + t397) * mrSges(6,3) + t581 * qJ(5) + t583 * qJD(5) + (t286 * t587 + t287 * t584 + t412) * g(3) + (Ifges(6,1) * t355 + t513) * t547 + (Ifges(6,2) * t357 + t514) * t548 + t492 * t555 + t493 * t557 + (Ifges(7,4) * t150 + Ifges(7,2) * t149) * t554 + t30 * t409 + (Ifges(7,1) * t150 + Ifges(7,4) * t149) * t550 + t139 * t533 + (Ifges(6,5) * t355 + Ifges(6,6) * t357) * t542 + (t140 - t242) * t535 + t603 * t84 + (t1 * t222 + t13 * t603 + t14 * t604 - t19 * t348 + t2 * t221 - t73 * t85) * m(7) + t604 * t83 + (-pkin(4) * t30 + qJ(5) * t397 + qJD(5) * t396 - t104 * t92 - t52 * t71 - t53 * t72) * m(6) + t150 * t561 - t285 * t562 + t149 * t563 + t355 * t564 + t357 * t565 + t300 * t568 + (Ifges(7,5) * t150 + Ifges(7,6) * t149) * t538 - t579 * t560 + (-Ifges(7,1) * t579 - Ifges(7,4) * t285) * t549 + (-Ifges(7,4) * t579 - Ifges(7,2) * t285) * t553 + (-Ifges(7,5) * t579 - Ifges(7,6) * t285) * t537 - t393 * t569 + (Ifges(7,5) * t300 - Ifges(7,6) * t393) * t545 + (Ifges(7,4) * t300 - Ifges(7,2) * t393) * t566 + (Ifges(7,1) * t300 - Ifges(7,4) * t393) * t567 + t19 * (mrSges(7,1) * t393 + mrSges(7,2) * t300) + (-t1 * t393 + t13 * t593 - t14 * t594 - t2 * t300) * mrSges(7,3) - t348 * t10 + t221 * t23 + t222 * t24 - t103 * t190 - t72 * t132 - t71 * t133 - t85 * t54 - pkin(4) * t58 + (-t400 * t536 - t402 * t541 - t404 * t540 - t586 + t637) * t247 + (Ifges(6,5) * t540 + Ifges(7,5) * t550 - Ifges(5,2) * t535 + Ifges(6,6) * t541 + Ifges(7,6) * t554 + Ifges(6,3) * t536 + Ifges(7,3) * t538 - t636) * t248 + t590 * t104 + (mrSges(7,1) * t594 - mrSges(7,2) * t593) * t73 - (-Ifges(5,1) * t247 - t504 + t608) * t248 / 0.2e1; -t432 * t132 + t180 * t133 - t614 * t83 + t100 * t84 + t10 + t58 + (t100 * t13 - t14 * t614 + t19 + t382) * m(7) + (t180 * t52 - t432 * t53 + t30 + t382) * m(6); -t73 * (mrSges(7,1) * t100 + mrSges(7,2) * t614) + (Ifges(7,1) * t614 - t512) * t550 + t47 * t549 + (Ifges(7,5) * t614 - Ifges(7,6) * t100) * t538 - t13 * t83 + t14 * t84 - g(1) * (mrSges(7,1) * t153 - mrSges(7,2) * t154) - g(2) * ((t218 * t350 + t289 * t351) * mrSges(7,1) + (t218 * t351 - t289 * t350) * mrSges(7,2)) - g(3) * ((-t287 * t350 + t351 * t490) * mrSges(7,1) + (-t287 * t351 - t350 * t490) * mrSges(7,2)) + (t100 * t14 + t13 * t614) * mrSges(7,3) + t6 + (-Ifges(7,2) * t100 + t48 + t96) * t554 + t622;];
tau  = t7;
