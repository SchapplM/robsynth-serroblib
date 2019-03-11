% Calculate vector of inverse dynamics joint torques for
% S6RRPRPR11
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 11:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPR11_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR11_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR11_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR11_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR11_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR11_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR11_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:11:45
% EndTime: 2019-03-09 11:12:42
% DurationCPUTime: 36.26s
% Computational Cost: add. (13166->903), mult. (27281->1196), div. (0->0), fcn. (18056->14), ass. (0->399)
t339 = sin(qJ(2));
t319 = t339 * qJD(1);
t310 = pkin(2) * t319;
t343 = cos(qJ(2));
t469 = qJ(3) * t343;
t376 = pkin(8) * t339 - t469;
t216 = qJD(1) * t376 + t310;
t439 = qJD(1) * t343;
t311 = pkin(7) * t439;
t266 = pkin(3) * t439 + t311;
t338 = sin(qJ(4));
t342 = cos(qJ(4));
t154 = -t216 * t338 + t342 * t266;
t431 = t342 * qJD(5);
t435 = qJD(4) * t338;
t519 = pkin(2) + pkin(8);
t444 = qJ(5) + t519;
t455 = t338 * t339;
t611 = -(pkin(4) * t343 - qJ(5) * t455) * qJD(1) - t154 + t435 * t444 - t431;
t155 = t342 * t216 + t338 * t266;
t274 = t444 * t342;
t416 = t342 * t319;
t610 = qJ(5) * t416 + qJD(4) * t274 + t338 * qJD(5) + t155;
t335 = cos(pkin(10));
t334 = sin(pkin(10));
t461 = t334 * t338;
t196 = -t319 * t461 + t335 * t416;
t434 = qJD(4) * t342;
t226 = t334 * t435 - t335 * t434;
t567 = t196 - t226;
t370 = t334 * t342 + t335 * t338;
t360 = t370 * t339;
t197 = qJD(1) * t360;
t227 = t370 * qJD(4);
t608 = t197 + t227;
t609 = -mrSges(6,3) - mrSges(7,3);
t574 = t610 * t334 + t335 * t611;
t573 = t334 * t611 - t610 * t335;
t337 = sin(qJ(6));
t341 = cos(qJ(6));
t295 = t319 + qJD(4);
t320 = t339 * qJ(3);
t406 = -pkin(1) - t320;
t200 = (-t343 * t519 + t406) * qJD(1);
t309 = pkin(7) * t319;
t265 = -pkin(3) * t319 - t309;
t552 = qJD(3) - t265;
t210 = -qJD(2) * t519 + t552;
t128 = -t200 * t338 + t342 * t210;
t437 = qJD(2) * t342;
t361 = t338 * t439 - t437;
t113 = qJ(5) * t361 + t128;
t105 = pkin(4) * t295 + t113;
t129 = t200 * t342 + t210 * t338;
t258 = -qJD(2) * t338 - t342 * t439;
t114 = qJ(5) * t258 + t129;
t107 = t334 * t114;
t53 = t335 * t105 - t107;
t169 = t258 * t334 - t335 * t361;
t595 = pkin(9) * t169;
t39 = pkin(5) * t295 + t53 - t595;
t460 = t335 * t114;
t54 = t334 * t105 + t460;
t400 = t335 * t258 + t334 * t361;
t585 = pkin(9) * t400;
t42 = t54 + t585;
t13 = -t337 * t42 + t341 * t39;
t430 = qJD(1) * qJD(2);
t269 = -t343 * qJDD(1) + t339 * t430;
t158 = qJD(4) * t258 + qJDD(2) * t342 + t269 * t338;
t270 = qJDD(1) * t339 + t343 * t430;
t255 = qJDD(4) + t270;
t468 = qJDD(1) * pkin(1);
t354 = -qJ(3) * t270 - qJD(3) * t319 - t468;
t121 = t269 * t519 + t354;
t254 = t270 * pkin(7);
t407 = qJDD(3) + t254;
t160 = pkin(3) * t270 - qJDD(2) * t519 + t407;
t52 = -qJD(4) * t129 - t121 * t338 + t342 * t160;
t36 = pkin(4) * t255 - qJ(5) * t158 + qJD(5) * t361 + t52;
t159 = qJD(4) * t361 - qJDD(2) * t338 + t269 * t342;
t51 = t342 * t121 + t338 * t160 - t200 * t435 + t210 * t434;
t38 = qJ(5) * t159 + qJD(5) * t258 + t51;
t11 = -t334 * t38 + t335 * t36;
t87 = t158 * t335 + t159 * t334;
t6 = pkin(5) * t255 - pkin(9) * t87 + t11;
t12 = t334 * t36 + t335 * t38;
t86 = -t158 * t334 + t159 * t335;
t9 = pkin(9) * t86 + t12;
t2 = qJD(6) * t13 + t337 * t6 + t341 * t9;
t607 = t2 * mrSges(7,2);
t14 = t337 * t39 + t341 * t42;
t3 = -qJD(6) * t14 - t337 * t9 + t341 * t6;
t606 = t3 * mrSges(7,1);
t605 = t11 * mrSges(6,1);
t604 = t12 * mrSges(6,2);
t583 = -Ifges(4,4) + Ifges(3,5);
t582 = Ifges(4,5) - Ifges(3,6);
t603 = Ifges(5,3) + Ifges(6,3);
t480 = t128 * mrSges(5,3);
t479 = t129 * mrSges(5,3);
t602 = -pkin(5) * t439 + pkin(9) * t608 + t574;
t601 = -pkin(9) * t567 + t573;
t123 = t196 * t341 - t197 * t337;
t369 = -t335 * t342 + t461;
t564 = -t337 * t370 - t341 * t369;
t89 = qJD(6) * t564 - t226 * t341 - t227 * t337;
t592 = t89 + t123;
t388 = t343 * mrSges(4,2) - t339 * mrSges(4,3);
t392 = mrSges(3,1) * t343 - mrSges(3,2) * t339;
t600 = t388 - t392;
t599 = qJD(3) + t309;
t329 = qJ(4) + pkin(10);
t315 = sin(t329);
t496 = pkin(4) * t338;
t271 = pkin(5) * t315 + t496;
t318 = qJ(6) + t329;
t301 = sin(t318);
t302 = cos(t318);
t316 = cos(t329);
t389 = mrSges(5,1) * t338 + mrSges(5,2) * t342;
t598 = -m(6) * t496 - m(7) * t271 - t315 * mrSges(6,1) - t301 * mrSges(7,1) - t316 * mrSges(6,2) - t302 * mrSges(7,2) - t389;
t336 = -qJ(5) - pkin(8);
t328 = -pkin(9) + t336;
t597 = -m(6) * (-pkin(2) + t336) + m(5) * t519 + mrSges(5,3) - m(7) * (-pkin(2) + t328) - t609;
t483 = Ifges(4,6) * t343;
t378 = -t339 * Ifges(4,2) - t483;
t596 = t14 * mrSges(7,2) + t54 * mrSges(6,2) + Ifges(4,4) * qJD(2) / 0.2e1 + qJD(1) * t378 / 0.2e1 - t13 * mrSges(7,1) - t53 * mrSges(6,1);
t340 = sin(qJ(1));
t594 = g(2) * t340;
t593 = t169 * Ifges(6,4);
t322 = t342 * pkin(4);
t305 = t322 + pkin(3);
t565 = pkin(4) * t434 + t305 * t319 + t599;
t520 = -m(6) - m(7);
t551 = t520 - m(4) - m(5);
t591 = t11 * t369 - t12 * t370 + t53 * t608 - t54 * t567;
t590 = t343 * t609 + t600;
t589 = -t169 * t337 + t341 * t400;
t99 = t169 * t341 + t337 * t400;
t118 = mrSges(5,1) * t255 - mrSges(5,3) * t158;
t119 = -mrSges(5,2) * t255 + mrSges(5,3) * t159;
t185 = -mrSges(5,2) * t295 + mrSges(5,3) * t258;
t186 = mrSges(5,1) * t295 + mrSges(5,3) * t361;
t374 = t342 * t185 - t338 * t186;
t588 = t374 * qJD(4) + t342 * t118 + t338 * t119;
t162 = t337 * t369 - t341 * t370;
t402 = t196 * t337 + t341 * t197;
t587 = -t13 * t402 - t162 * t2 + t3 * t564;
t26 = qJD(6) * t589 + t337 * t86 + t341 * t87;
t539 = t26 / 0.2e1;
t27 = -qJD(6) * t99 - t337 * t87 + t341 * t86;
t538 = t27 / 0.2e1;
t526 = t86 / 0.2e1;
t525 = t87 / 0.2e1;
t245 = qJDD(6) + t255;
t511 = t245 / 0.2e1;
t510 = t255 / 0.2e1;
t273 = t444 * t338;
t173 = t273 * t334 - t335 * t274;
t135 = pkin(9) * t369 + t173;
t174 = -t335 * t273 - t334 * t274;
t136 = -pkin(9) * t370 + t174;
t70 = t135 * t341 - t136 * t337;
t581 = qJD(6) * t70 + t337 * t602 + t341 * t601;
t71 = t135 * t337 + t136 * t341;
t580 = -qJD(6) * t71 - t337 * t601 + t341 * t602;
t579 = Ifges(6,4) * t400;
t497 = pkin(4) * t335;
t303 = pkin(5) + t497;
t498 = pkin(4) * t334;
t221 = t303 * t337 + t341 * t498;
t58 = -t113 * t334 - t460;
t43 = t58 - t585;
t59 = t335 * t113 - t107;
t44 = t59 - t595;
t577 = -t221 * qJD(6) + t337 * t44 - t341 * t43;
t220 = t303 * t341 - t337 * t498;
t576 = t220 * qJD(6) - t337 * t43 - t341 * t44;
t540 = m(6) * pkin(4);
t575 = -t540 - mrSges(5,1);
t344 = cos(qJ(1));
t553 = g(1) * t344 + t594;
t570 = t343 * t553;
t569 = pkin(5) * t567 + t565;
t175 = -mrSges(5,1) * t258 - mrSges(5,2) * t361;
t419 = mrSges(4,1) * t439;
t282 = -qJD(2) * mrSges(4,3) - t419;
t568 = t175 - t282;
t325 = t343 * pkin(2);
t441 = t325 + t320;
t396 = pkin(8) * t343 + t441;
t250 = -pkin(1) - t396;
t518 = pkin(3) + pkin(7);
t285 = t518 * t339;
t172 = t342 * t250 + t338 * t285;
t563 = qJD(2) * mrSges(3,2) - mrSges(3,3) * t439 + t282;
t420 = mrSges(4,1) * t319;
t562 = mrSges(3,3) * t319 + t420 + (-mrSges(3,1) + mrSges(4,2)) * qJD(2);
t350 = qJD(6) * t162 + t226 * t337 - t341 * t227;
t561 = t350 - t402;
t484 = Ifges(4,6) * t339;
t560 = t339 * (-Ifges(4,2) * t343 + t484) + t343 * (Ifges(4,3) * t339 - t483);
t559 = t339 * t582 + t343 * t583;
t558 = Ifges(5,5) * t158 + Ifges(6,5) * t87 + Ifges(5,6) * t159 + Ifges(6,6) * t86 + t255 * t603;
t253 = t269 * pkin(7);
t556 = -t253 * t343 + t254 * t339;
t201 = -qJDD(2) * qJ(3) - qJD(2) * qJD(3) + t253;
t215 = -qJDD(2) * pkin(2) + t407;
t555 = -t201 * t343 + t215 * t339;
t554 = -t338 * t51 - t342 * t52;
t476 = t361 * Ifges(5,4);
t147 = t258 * Ifges(5,2) + t295 * Ifges(5,6) - t476;
t251 = Ifges(5,4) * t258;
t148 = -Ifges(5,1) * t361 + t295 * Ifges(5,5) + t251;
t377 = -Ifges(4,3) * t343 - t484;
t549 = Ifges(4,5) * qJD(2) + qJD(1) * t377 + t342 * t147 + t338 * t148;
t287 = qJD(6) + t295;
t308 = Ifges(3,4) * t439;
t548 = Ifges(3,1) * t319 + Ifges(3,5) * qJD(2) - Ifges(5,5) * t361 + t169 * Ifges(6,5) + t99 * Ifges(7,5) + t258 * Ifges(5,6) + Ifges(6,6) * t400 + t589 * Ifges(7,6) + t287 * Ifges(7,3) + t295 * t603 + t308;
t473 = t343 * mrSges(4,3);
t547 = -t473 + t598 * t343 + (m(4) * pkin(2) - mrSges(4,2) + t597) * t339;
t272 = pkin(5) * t316 + t322;
t544 = -m(5) * pkin(3) - m(6) * t305 - m(7) * (pkin(3) + t272) - mrSges(4,1) + mrSges(2,2) - mrSges(3,3);
t542 = Ifges(7,4) * t539 + Ifges(7,2) * t538 + Ifges(7,6) * t511;
t541 = Ifges(7,1) * t539 + Ifges(7,4) * t538 + Ifges(7,5) * t511;
t537 = Ifges(6,4) * t525 + Ifges(6,2) * t526 + Ifges(6,6) * t510;
t536 = Ifges(6,1) * t525 + Ifges(6,4) * t526 + Ifges(6,5) * t510;
t500 = Ifges(7,4) * t99;
t46 = Ifges(7,2) * t589 + Ifges(7,6) * t287 + t500;
t535 = -t46 / 0.2e1;
t534 = t46 / 0.2e1;
t88 = Ifges(7,4) * t589;
t47 = Ifges(7,1) * t99 + Ifges(7,5) * t287 + t88;
t533 = -t47 / 0.2e1;
t532 = t47 / 0.2e1;
t531 = -t158 * Ifges(5,4) / 0.2e1 - t159 * Ifges(5,2) / 0.2e1 - t255 * Ifges(5,6) / 0.2e1;
t84 = Ifges(6,2) * t400 + t295 * Ifges(6,6) + t593;
t530 = -t84 / 0.2e1;
t529 = t84 / 0.2e1;
t85 = t169 * Ifges(6,1) + t295 * Ifges(6,5) + t579;
t528 = -t85 / 0.2e1;
t527 = t85 / 0.2e1;
t524 = -t589 / 0.2e1;
t523 = t589 / 0.2e1;
t522 = -t99 / 0.2e1;
t521 = t99 / 0.2e1;
t517 = t158 / 0.2e1;
t516 = t159 / 0.2e1;
t515 = -t400 / 0.2e1;
t514 = t400 / 0.2e1;
t513 = -t169 / 0.2e1;
t512 = t169 / 0.2e1;
t508 = -t361 / 0.2e1;
t507 = -t287 / 0.2e1;
t506 = t287 / 0.2e1;
t505 = -t295 / 0.2e1;
t504 = t295 / 0.2e1;
t502 = mrSges(7,3) * t13;
t501 = mrSges(7,3) * t14;
t499 = pkin(4) * t361;
t495 = pkin(7) * t339;
t492 = g(3) * t343;
t323 = t343 * pkin(7);
t438 = qJD(2) * t339;
t399 = pkin(2) * t438 - qJD(3) * t339;
t187 = qJD(2) * t376 + t399;
t436 = qJD(2) * t343;
t268 = t518 * t436;
t242 = t342 * t268;
t403 = qJ(5) * t343 - t250;
t65 = pkin(4) * t436 + t242 + t403 * t434 + (-qJ(5) * t438 - qJD(4) * t285 + qJD(5) * t343 - t187) * t338;
t433 = qJD(4) * t343;
t415 = t338 * t433;
t357 = t339 * t437 + t415;
t94 = t342 * t187 - t250 * t435 + t338 * t268 + t285 * t434;
t69 = qJ(5) * t357 - t343 * t431 + t94;
t35 = t334 * t65 + t335 * t69;
t489 = mrSges(7,1) * t302;
t488 = Ifges(3,4) * t339;
t487 = Ifges(3,4) * t343;
t486 = Ifges(5,4) * t338;
t485 = Ifges(5,4) * t342;
t462 = t328 * t343;
t459 = t336 * t343;
t454 = t338 * t343;
t453 = t339 * t340;
t452 = t339 * t344;
t451 = t340 * t342;
t447 = t342 * t343;
t446 = t342 * t344;
t445 = t343 * t344;
t261 = t342 * t285;
t142 = pkin(4) * t339 + t338 * t403 + t261;
t149 = -qJ(5) * t447 + t172;
t77 = t334 * t142 + t335 * t149;
t189 = -t301 * t340 + t302 * t452;
t190 = t301 * t452 + t302 * t340;
t443 = t189 * mrSges(7,1) - t190 * mrSges(7,2);
t191 = t301 * t344 + t302 * t453;
t192 = -t301 * t453 + t302 * t344;
t442 = t191 * mrSges(7,1) + t192 * mrSges(7,2);
t286 = t343 * pkin(3) + t323;
t440 = t344 * pkin(1) + t340 * pkin(7);
t425 = Ifges(7,5) * t26 + Ifges(7,6) * t27 + Ifges(7,3) * t245;
t333 = qJD(2) * qJ(3);
t233 = t333 + t266;
t228 = pkin(4) * t447 + t286;
t48 = -t86 * mrSges(6,1) + t87 * mrSges(6,2);
t10 = -t27 * mrSges(7,1) + t26 * mrSges(7,2);
t409 = -t434 / 0.2e1;
t304 = qJ(3) + t496;
t34 = -t334 * t69 + t335 * t65;
t405 = -t430 / 0.2e1;
t219 = t270 * mrSges(4,1) + qJDD(2) * mrSges(4,2);
t76 = t335 * t142 - t149 * t334;
t391 = mrSges(3,1) * t339 + mrSges(3,2) * t343;
t390 = mrSges(5,1) * t342 - mrSges(5,2) * t338;
t387 = Ifges(5,1) * t342 - t486;
t386 = Ifges(5,1) * t338 + t485;
t385 = Ifges(3,2) * t343 + t488;
t383 = -Ifges(5,2) * t338 + t485;
t382 = Ifges(5,2) * t342 + t486;
t380 = Ifges(5,5) * t342 - Ifges(5,6) * t338;
t379 = Ifges(5,5) * t338 + Ifges(5,6) * t342;
t214 = t370 * t343;
t61 = pkin(5) * t339 + pkin(9) * t214 + t76;
t213 = t369 * t343;
t62 = pkin(9) * t213 + t77;
t30 = -t337 * t62 + t341 * t61;
t31 = t337 * t61 + t341 * t62;
t375 = t128 * t338 - t129 * t342;
t139 = t213 * t341 + t214 * t337;
t140 = t213 * t337 - t214 * t341;
t275 = -qJD(2) * pkin(2) + t599;
t280 = -t311 - t333;
t371 = t275 * t343 + t280 * t339;
t368 = t406 - t325;
t367 = t425 + t606 - t607;
t366 = pkin(1) * t391;
t235 = -t338 * t340 + t339 * t446;
t237 = t338 * t344 + t339 * t451;
t234 = t368 * qJD(1);
t365 = t234 * (-mrSges(4,2) * t339 - t473);
t364 = t339 * (Ifges(3,1) * t343 - t488);
t170 = -pkin(4) * t258 + qJD(5) + t233;
t356 = t338 * t438 - t342 * t433;
t353 = Ifges(5,5) * t343 + t339 * t386;
t352 = Ifges(5,6) * t343 + t339 * t382;
t351 = Ifges(5,3) * t343 + t339 * t379;
t165 = -pkin(3) * t269 - t201;
t176 = -pkin(4) * t415 + (-pkin(7) - t305) * t438;
t349 = -qJD(4) * t375 - t554;
t102 = -pkin(4) * t159 + qJDD(5) + t165;
t277 = -pkin(1) - t441;
t276 = t343 * t301 * mrSges(7,2);
t267 = t518 * t438;
t264 = -qJ(3) * t439 + t310;
t263 = t388 * qJD(1);
t246 = t390 * t343;
t238 = -t338 * t453 + t446;
t236 = t338 * t452 + t451;
t229 = Ifges(3,6) * qJD(2) + qJD(1) * t385;
t222 = -qJ(3) * t436 + t399;
t218 = mrSges(4,1) * t269 - qJDD(2) * mrSges(4,3);
t208 = -t315 * t453 + t316 * t344;
t207 = t315 * t344 + t316 * t453;
t206 = t315 * t452 + t316 * t340;
t205 = -t315 * t340 + t316 * t452;
t199 = pkin(5) * t370 + t304;
t171 = -t250 * t338 + t261;
t157 = -pkin(5) * t213 + t228;
t156 = pkin(2) * t269 + t354;
t144 = qJD(2) * t360 + qJD(4) * t213;
t143 = -t369 * t438 + t370 * t433;
t138 = mrSges(6,1) * t295 - mrSges(6,3) * t169;
t137 = -mrSges(6,2) * t295 + mrSges(6,3) * t400;
t122 = pkin(5) * t169 - t499;
t106 = -pkin(5) * t400 + t170;
t104 = -mrSges(6,1) * t400 + mrSges(6,2) * t169;
t103 = -pkin(5) * t143 + t176;
t95 = -qJD(4) * t172 - t187 * t338 + t242;
t93 = -mrSges(5,1) * t159 + mrSges(5,2) * t158;
t79 = mrSges(7,1) * t287 - mrSges(7,3) * t99;
t78 = -mrSges(7,2) * t287 + mrSges(7,3) * t589;
t75 = t158 * Ifges(5,1) + t159 * Ifges(5,4) + t255 * Ifges(5,5);
t73 = mrSges(6,1) * t255 - mrSges(6,3) * t87;
t72 = -mrSges(6,2) * t255 + mrSges(6,3) * t86;
t56 = -qJD(6) * t140 + t143 * t341 - t144 * t337;
t55 = qJD(6) * t139 + t143 * t337 + t144 * t341;
t50 = -pkin(5) * t86 + t102;
t49 = -mrSges(7,1) * t589 + mrSges(7,2) * t99;
t22 = pkin(9) * t143 + t35;
t21 = -mrSges(7,2) * t245 + mrSges(7,3) * t27;
t20 = mrSges(7,1) * t245 - mrSges(7,3) * t26;
t19 = pkin(5) * t436 - pkin(9) * t144 + t34;
t5 = -qJD(6) * t31 + t19 * t341 - t22 * t337;
t4 = qJD(6) * t30 + t19 * t337 + t22 * t341;
t1 = [(Ifges(6,4) * t144 + Ifges(6,2) * t143) * t514 + (Ifges(6,1) * t144 + Ifges(6,4) * t143) * t512 + t392 * t468 + (t343 * (-Ifges(3,2) * t339 + t487) + t364) * t430 / 0.2e1 + t555 * mrSges(4,1) + (t563 * pkin(7) - t229 / 0.2e1 + t549 / 0.2e1 + t280 * mrSges(4,1)) * t438 + (Ifges(7,1) * t55 + Ifges(7,4) * t56) * t521 + (Ifges(7,5) * t55 + Ifges(7,6) * t56) * t506 + (Ifges(7,4) * t55 + Ifges(7,2) * t56) * t523 + t258 * (qJD(2) * t352 - t383 * t433) / 0.2e1 + (qJD(2) * t353 - t387 * t433) * t508 + (Ifges(7,5) * t140 + Ifges(7,6) * t139 + Ifges(7,3) * t339) * t511 + t270 * t487 / 0.2e1 + (-Ifges(3,4) * t269 + Ifges(3,5) * qJDD(2) + t425 + t558) * t339 / 0.2e1 + (-m(3) * t440 - t236 * mrSges(5,1) - t206 * mrSges(6,1) - t190 * mrSges(7,1) - t235 * mrSges(5,2) - t205 * mrSges(6,2) - t189 * mrSges(7,2) + (-m(5) * pkin(8) - mrSges(5,3)) * t445 + t551 * (pkin(2) * t445 + qJ(3) * t452 + t440) + t544 * t340 + (-m(7) * (t271 * t339 - t462) - m(6) * (pkin(4) * t455 - t459) - mrSges(2,1) + t590) * t344) * g(2) - t269 * t385 / 0.2e1 + t156 * t388 + t269 * t377 / 0.2e1 - t270 * t378 / 0.2e1 + t106 * (-mrSges(7,1) * t56 + mrSges(7,2) * t55) + t103 * t49 + (-qJDD(2) * mrSges(3,2) - t218) * t323 + t51 * (-mrSges(5,2) * t339 - mrSges(5,3) * t447) + t357 * t479 + t339 * t605 + t339 * t606 + (Ifges(6,5) * t144 + Ifges(6,6) * t143 + qJD(2) * t351 - t380 * t433) * t504 + t270 * t339 * Ifges(3,1) - t366 * t430 - t75 * t454 / 0.2e1 + t52 * (mrSges(5,1) * t339 + mrSges(5,3) * t454) + m(7) * (t103 * t106 + t13 * t5 + t14 * t4 + t157 * t50 + t2 * t31 + t3 * t30) + m(6) * (t102 * t228 + t11 * t76 + t12 * t77 + t170 * t176 + t34 * t53 + t35 * t54) + t147 * t415 / 0.2e1 + t76 * t73 + t77 * t72 + t4 * t78 + t5 * t79 + t30 * t20 + Ifges(2,3) * qJDD(1) + t31 * t21 + (-t13 * t55 + t139 * t2 + t14 * t56 - t140 * t3) * mrSges(7,3) + (t11 * t214 + t12 * t213 + t143 * t54 - t144 * t53) * mrSges(6,3) - t339 * t607 - t339 * t604 + (-Ifges(6,5) * t214 + Ifges(6,6) * t213 + t339 * t603 - t343 * t379) * t510 - t356 * t480 + t343 * (Ifges(3,4) * t270 - Ifges(3,2) * t269 + Ifges(3,6) * qJDD(2)) / 0.2e1 - t343 * (Ifges(4,5) * qJDD(2) - Ifges(4,6) * t270 + Ifges(4,3) * t269) / 0.2e1 + (-qJDD(2) * mrSges(3,1) + t219) * t495 - t339 * (Ifges(4,4) * qJDD(2) - Ifges(4,2) * t270 + Ifges(4,6) * t269) / 0.2e1 + (t339 * t583 - t343 * t582) * qJDD(2) / 0.2e1 + (Ifges(5,6) * t339 - t343 * t382) * t516 + (Ifges(5,5) * t339 - t343 * t386) * t517 + (-Ifges(6,1) * t214 + Ifges(6,4) * t213 + Ifges(6,5) * t339) * t525 + (-Ifges(6,4) * t214 + Ifges(6,2) * t213 + Ifges(6,6) * t339) * t526 + t102 * (-mrSges(6,1) * t213 - mrSges(6,2) * t214) + t35 * t137 + t34 * t138 + t50 * (-mrSges(7,1) * t139 + mrSges(7,2) * t140) + m(5) * (t128 * t95 + t129 * t94 + t165 * t286 + t171 * t52 + t172 * t51 - t233 * t267) + t157 * t10 + t343 * t148 * t409 + t170 * (-mrSges(6,1) * t143 + mrSges(6,2) * t144) + t171 * t118 + t172 * t119 + t176 * t104 + t94 * t185 + t95 * t186 + t233 * (-mrSges(5,1) * t357 + mrSges(5,2) * t356) + t144 * t527 + t143 * t529 + t447 * t531 + t55 * t532 + t56 * t534 - t214 * t536 + t213 * t537 + (Ifges(7,4) * t140 + Ifges(7,2) * t139 + Ifges(7,6) * t339) * t538 + (Ifges(7,1) * t140 + Ifges(7,4) * t139 + Ifges(7,5) * t339) * t539 + t140 * t541 + t139 * t542 + qJD(2) * t365 + (t562 * pkin(7) + Ifges(6,5) * t512 + Ifges(6,3) * t504 + t128 * mrSges(5,1) - t129 * mrSges(5,2) + Ifges(6,6) * t514 + Ifges(7,5) * t521 + Ifges(7,6) * t523 + t548 / 0.2e1 + t275 * mrSges(4,1) + Ifges(7,3) * t506 - t596) * t436 + t228 * t48 + t165 * t246 + t222 * t263 - t267 * t175 + m(4) * (t156 * t277 + t222 * t234 + (qJD(2) * t371 + t555) * pkin(7)) + (-t269 * t323 + t270 * t495 + t556) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t556) + t559 * qJD(2) ^ 2 / 0.2e1 + t560 * t405 - pkin(1) * (mrSges(3,1) * t269 + mrSges(3,2) * t270) + t277 * (-mrSges(4,2) * t269 - mrSges(4,3) * t270) + t286 * t93 + (-t238 * mrSges(5,1) - t208 * mrSges(6,1) - t192 * mrSges(7,1) + t237 * mrSges(5,2) + t207 * mrSges(6,2) + t191 * mrSges(7,2) + (m(3) * pkin(1) - m(6) * (-t304 * t339 - pkin(1)) - m(4) * t368 + mrSges(2,1) - m(5) * t406 - m(7) * (-pkin(1) + (-qJ(3) - t271) * t339) + t597 * t343 - t600) * t340 + ((-m(3) + t551) * pkin(7) + t544) * t344) * g(1); t295 * t390 * t233 - t370 * t537 + (-Ifges(6,1) * t369 - Ifges(6,4) * t370) * t525 + (-Ifges(6,4) * t369 - Ifges(6,2) * t370) * t526 + (-Ifges(6,5) * t369 - Ifges(6,6) * t370 + t380) * t510 + t102 * (mrSges(6,1) * t370 - mrSges(6,2) * t369) - t369 * t536 - (t258 * t382 + t295 * t379 - t361 * t386) * qJD(4) / 0.2e1 - (t258 * t352 + t295 * t351 - t353 * t361) * qJD(1) / 0.2e1 + (Ifges(7,5) * t564 + Ifges(7,6) * t162) * t511 + t50 * (-mrSges(7,1) * t162 + mrSges(7,2) * t564) + (Ifges(7,4) * t564 + Ifges(7,2) * t162) * t538 + (Ifges(7,1) * t564 + Ifges(7,4) * t162) * t539 + t564 * t541 + (t165 * qJ(3) - t128 * t154 - t129 * t155 + t233 * t552) * m(5) + (Ifges(7,1) * t521 + Ifges(7,4) * t523 + Ifges(7,5) * t506 - t502 + t532) * t350 + (Ifges(7,5) * t402 + Ifges(7,6) * t123) * t507 + (Ifges(7,1) * t402 + Ifges(7,4) * t123) * t522 + (Ifges(7,4) * t402 + Ifges(7,2) * t123) * t524 + t402 * t533 - (Ifges(7,4) * t521 + Ifges(7,2) * t523 + Ifges(7,6) * t506 + t501 + t534) * t89 + t165 * t389 + (Ifges(3,3) + Ifges(4,1)) * qJDD(2) - t275 * t419 - t280 * t420 + (mrSges(6,1) * t567 - mrSges(6,2) * t608) * t170 + t591 * mrSges(6,3) + (t480 - t148 / 0.2e1) * t435 + (Ifges(6,4) * t197 + Ifges(6,2) * t196) * t515 + t229 * t319 / 0.2e1 + (Ifges(6,1) * t197 + Ifges(6,4) * t196) * t513 + (Ifges(6,5) * t197 + Ifges(6,6) * t196) * t505 + (-m(4) * t371 * pkin(7) - t365 - t128 * (mrSges(5,1) * t343 - mrSges(5,3) * t455) - t129 * (mrSges(5,3) * t339 * t342 - mrSges(5,2) * t343)) * qJD(1) + (-t364 / 0.2e1 + t366 + t560 / 0.2e1) * qJD(1) ^ 2 + (t93 - t218) * qJ(3) + (-pkin(2) * t215 - qJ(3) * t201 - qJD(3) * t280 - t234 * t264) * m(4) + (-t123 * t14 - t587) * mrSges(7,3) + (-m(5) * t349 - t588) * t519 - t434 * t479 + t70 * t20 + t71 * t21 + t342 * t75 / 0.2e1 + (mrSges(7,1) * t592 + t561 * mrSges(7,2)) * t106 + t580 * t79 + t581 * t78 + (t569 * t106 + t580 * t13 + t581 * t14 + t199 * t50 + t2 * t71 + t3 * t70) * m(7) + t582 * t269 + t583 * t270 + t383 * t516 + t387 * t517 + (-Ifges(6,1) * t227 + Ifges(6,4) * t226) * t512 + (-Ifges(6,4) * t227 + Ifges(6,2) * t226) * t514 + (-Ifges(6,5) * t227 + Ifges(6,6) * t226) * t504 + t173 * t73 + t174 * t72 - t155 * t185 - t154 * t186 + t199 * t10 - t201 * mrSges(4,3) + t147 * t409 - t227 * t527 + t197 * t528 + t226 * t529 + t196 * t530 + t338 * t531 + t123 * t535 + t162 * t542 + (t469 * t551 + t547) * t594 + (Ifges(6,5) * t513 + Ifges(7,5) * t522 + Ifges(6,6) * t515 + Ifges(7,6) * t524 + Ifges(6,3) * t505 + Ifges(7,3) * t507 + t596) * t439 + t215 * mrSges(4,2) - (-Ifges(3,2) * t319 + t308 + t548) * t439 / 0.2e1 - pkin(2) * t219 - t549 * t319 / 0.2e1 + t253 * mrSges(3,2) - t254 * mrSges(3,1) + (qJ(3) * t445 * t551 + t344 * t547) * g(1) - t264 * t263 - t265 * t175 + t553 * t391 + t554 * mrSges(5,3) + t559 * t405 - t562 * t311 - t563 * t309 + t565 * t104 + t568 * qJD(3) + t569 * t49 + t304 * t48 + (-m(5) * t396 - t343 * mrSges(5,3) - m(4) * t441 - m(7) * (t441 - t462) - m(6) * (t441 - t459) + t598 * t339 + t590) * g(3) + t573 * t137 + t574 * t138 + (t102 * t304 + t11 * t173 + t12 * t174 + t170 * t565 + t53 * t574 + t54 * t573) * m(6); t561 * t79 + t564 * t20 + t592 * t78 + t219 + t567 * t137 + (-t104 - t49 - t568) * qJD(2) - t551 * t492 - t608 * t138 - t162 * t21 + t370 * t72 - t369 * t73 + ((t263 + t374) * qJD(1) + t553 * t551) * t339 + (-qJD(2) * t106 + t13 * t350 + t14 * t592 + t587) * m(7) + (-qJD(2) * t170 - t591) * m(6) + (-qJD(2) * t233 - t319 * t375 + t349) * m(5) + (qJD(2) * t280 + t234 * t319 + t215) * m(4) + t588; (Ifges(6,1) * t400 - t593) * t513 + (Ifges(5,5) * t258 + Ifges(6,5) * t400 + Ifges(5,6) * t361 - Ifges(6,6) * t169) * t505 - t361 * t479 - (Ifges(5,2) * t361 + t148 + t251) * t258 / 0.2e1 - t233 * (-mrSges(5,1) * t361 + mrSges(5,2) * t258) + t361 * (Ifges(5,1) * t258 + t476) / 0.2e1 + t400 * t528 - (mrSges(7,1) * t106 + Ifges(7,4) * t522 + Ifges(7,2) * t524 + Ifges(7,6) * t507 - t501 + t535) * t99 + (-mrSges(7,2) * t106 + Ifges(7,1) * t522 + Ifges(7,4) * t524 + Ifges(7,5) * t507 + t502 + t533) * t589 + t367 + t605 + t104 * t499 + t147 * t508 - (-m(6) * t322 - mrSges(6,1) * t316 + mrSges(6,2) * t315) * t492 + t558 - t604 + (-Ifges(6,2) * t169 + t579) * t515 + (t169 * t54 + t400 * t53) * mrSges(6,3) - t170 * (mrSges(6,1) * t169 + mrSges(6,2) * t400) - t169 * t530 + t73 * t497 + t258 * t480 + t72 * t498 + (-t276 - (-m(7) * t272 - t489) * t343 + t246) * g(3) + t52 * mrSges(5,1) - t51 * mrSges(5,2) - t122 * t49 - t59 * t137 - t58 * t138 - m(6) * (-t170 * t499 + t53 * t58 + t54 * t59) - t128 * t185 + t129 * t186 + (t11 * t335 + t12 * t334) * t540 + t220 * t20 + t221 * t21 + (-t205 * mrSges(6,1) + t206 * mrSges(6,2) - m(7) * (-t271 * t340 + t272 * t452) - t443 + mrSges(5,2) * t236 + t575 * t235) * g(1) + (-t207 * mrSges(6,1) - t208 * mrSges(6,2) - m(7) * (t271 * t344 + t272 * t453) - t442 - mrSges(5,2) * t238 + t575 * t237) * g(2) + t576 * t78 + t577 * t79 + (-t106 * t122 + t13 * t577 + t14 * t576 + t2 * t221 + t220 * t3) * m(7); t520 * t339 * g(3) - t400 * t137 + t169 * t138 - t589 * t78 + t99 * t79 + t10 + t48 + (t13 * t99 - t14 * t589 + t50 - t570) * m(7) + (t169 * t53 - t400 * t54 + t102 - t570) * m(6); -t106 * (mrSges(7,1) * t99 + mrSges(7,2) * t589) + (Ifges(7,1) * t589 - t500) * t522 + t46 * t521 + (Ifges(7,5) * t589 - Ifges(7,6) * t99) * t507 - t13 * t78 + t14 * t79 - g(1) * t443 - g(2) * t442 - g(3) * (-t343 * t489 + t276) + (t13 * t589 + t14 * t99) * mrSges(7,3) + t367 + (-Ifges(7,2) * t99 + t47 + t88) * t524;];
tau  = t1;
