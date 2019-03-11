% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPRRR8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR8_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR8_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR8_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR8_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:35:37
% EndTime: 2019-03-09 02:35:49
% DurationCPUTime: 7.93s
% Computational Cost: add. (22887->551), mult. (42951->743), div. (0->0), fcn. (48286->8), ass. (0->309)
t341 = sin(pkin(10));
t518 = pkin(1) + qJ(3);
t435 = -pkin(7) - t518;
t317 = t435 * t341;
t345 = sin(qJ(4));
t342 = cos(pkin(10));
t391 = t342 * t435;
t532 = cos(qJ(4));
t252 = t317 * t345 - t532 * t391;
t309 = t341 * t532 + t345 * t342;
t308 = t341 * t345 - t342 * t532;
t530 = pkin(4) * t308;
t254 = pkin(8) * t309 - t530;
t344 = sin(qJ(5));
t346 = cos(qJ(5));
t152 = t252 * t344 + t346 * t254;
t153 = -t346 * t252 + t344 * t254;
t343 = sin(qJ(6));
t531 = cos(qJ(6));
t583 = -t343 * t344 + t531 * t346;
t608 = t583 * t309;
t167 = -mrSges(7,1) * t308 + mrSges(7,3) * t608;
t460 = t309 * t344;
t403 = t460 / 0.2e1;
t335 = Ifges(6,5) * t346;
t422 = -t335 / 0.2e1;
t430 = pkin(5) * t531;
t568 = m(7) * pkin(5);
t433 = t568 / 0.2e1;
t528 = pkin(5) * t343;
t543 = -t308 / 0.2e1;
t567 = mrSges(6,1) / 0.2e1;
t556 = -t608 / 0.2e1;
t378 = t343 * t346 + t344 * t531;
t590 = t378 * t309;
t558 = t590 / 0.2e1;
t380 = Ifges(7,5) * t556 + Ifges(7,6) * t558;
t542 = t308 / 0.2e1;
t124 = pkin(9) * t460 + t153;
t459 = t309 * t346;
t99 = -pkin(5) * t308 + pkin(9) * t459 + t152;
t69 = -t343 * t124 + t531 * t99;
t70 = t124 * t531 + t343 * t99;
t578 = Ifges(7,3) * t542 + t70 * mrSges(7,2) / 0.2e1 - t69 * mrSges(7,1) / 0.2e1 - t380;
t165 = mrSges(7,2) * t308 + mrSges(7,3) * t590;
t607 = t165 / 0.2e1;
t613 = Ifges(6,3) * t543 + t152 * t567 - t153 * mrSges(6,2) / 0.2e1 + (t343 * t70 + t531 * t69) * t433 + t309 * t422 + Ifges(6,6) * t403 + t528 * t607 + t167 * t430 / 0.2e1 - t578;
t488 = t346 * mrSges(6,1);
t489 = t344 * mrSges(6,2);
t491 = t378 * mrSges(7,2);
t493 = t583 * mrSges(7,1);
t586 = t493 / 0.2e1 - t491 / 0.2e1;
t612 = t489 / 0.2e1 - t488 / 0.2e1 - (t343 * t378 + t531 * t583) * t433 - t586;
t611 = t378 * t608 - t583 * t590;
t327 = t341 * pkin(3) + qJ(2);
t525 = pkin(8) * t308;
t529 = pkin(4) * t309;
t235 = t327 + t525 + t529;
t253 = t317 * t532 + t345 * t391;
t141 = t235 * t344 + t253 * t346;
t463 = t308 * t344;
t120 = pkin(9) * t463 + t141;
t456 = t343 * t120;
t140 = t346 * t235 - t253 * t344;
t462 = t308 * t346;
t119 = pkin(9) * t462 + t140;
t97 = pkin(5) * t309 + t119;
t65 = t531 * t97 - t456;
t74 = t119 * t531 - t456;
t602 = t65 - t74;
t339 = t344 ^ 2;
t340 = t346 ^ 2;
t436 = t339 + t340;
t610 = mrSges(6,3) * t436;
t566 = -pkin(9) - pkin(8);
t325 = t566 * t344;
t326 = t566 * t346;
t276 = t343 * t325 - t326 * t531;
t393 = t531 * t325 + t326 * t343;
t438 = Ifges(7,5) * t583 - Ifges(7,6) * t378;
t49 = -t276 * mrSges(7,1) - t393 * mrSges(7,2) + t438;
t609 = t49 * qJD(6);
t220 = t583 * t308;
t223 = t378 * t308;
t130 = -mrSges(7,1) * t220 + mrSges(7,2) * t223;
t166 = -mrSges(7,2) * t309 + t223 * mrSges(7,3);
t470 = t590 * t223;
t472 = t608 * t220;
t559 = -t590 / 0.2e1;
t168 = mrSges(7,1) * t309 + t220 * mrSges(7,3);
t561 = t168 / 0.2e1;
t355 = t130 * t542 + t166 * t559 - t608 * t561 + (t472 / 0.2e1 + t470 / 0.2e1) * mrSges(7,3);
t606 = -t168 / 0.2e1;
t605 = t276 / 0.2e1;
t540 = t583 / 0.2e1;
t487 = t346 * mrSges(6,2);
t320 = t344 * mrSges(6,1) + t487;
t464 = t308 * t320;
t604 = -t464 / 0.2e1;
t527 = pkin(5) * t344;
t411 = t531 * t120;
t66 = t343 * t97 + t411;
t73 = -t343 * t119 - t411;
t593 = t66 + t73;
t319 = -t488 + t489;
t236 = t308 * t319;
t272 = t308 ^ 2;
t526 = pkin(5) * t346;
t569 = m(7) / 0.2e1;
t599 = (-t272 * t526 - t590 * t593 - t602 * t608) * t569 + t236 * t542 + t355;
t396 = -mrSges(7,1) * t608 + mrSges(7,2) * t590;
t597 = qJD(6) * t396;
t257 = mrSges(7,1) * t378 + mrSges(7,2) * t583;
t596 = t257 * qJD(6);
t431 = pkin(5) * t463;
t192 = t252 - t431;
t336 = Ifges(6,4) * t346;
t323 = Ifges(6,1) * t344 + t336;
t240 = t308 * t323;
t331 = -pkin(4) - t526;
t541 = t309 / 0.4e1;
t368 = t192 * t257 / 0.2e1 + t438 * t541 + t331 * t130 / 0.2e1;
t506 = Ifges(6,6) * t344;
t395 = t335 - t506;
t508 = Ifges(7,4) * t378;
t260 = Ifges(7,2) * t583 + t508;
t261 = Ifges(7,1) * t583 - t508;
t398 = t260 / 0.4e1 - t261 / 0.4e1;
t307 = Ifges(7,4) * t583;
t259 = -Ifges(7,2) * t378 + t307;
t262 = Ifges(7,1) * t378 + t307;
t399 = t259 / 0.4e1 + t262 / 0.4e1;
t246 = -mrSges(6,2) * t309 + mrSges(6,3) * t463;
t453 = t344 * t246;
t402 = -t453 / 0.2e1;
t384 = Ifges(6,2) * t344 - t336;
t494 = t309 * Ifges(6,6);
t189 = t308 * t384 + t494;
t454 = t344 * t189;
t545 = t393 / 0.2e1;
t502 = t223 * mrSges(7,1);
t504 = t220 * mrSges(7,2);
t132 = -t502 - t504;
t564 = t132 / 0.2e1;
t595 = -pkin(4) * t236 / 0.2e1 + t252 * t320 / 0.2e1 + t166 * t545 + t395 * t541 - t454 / 0.4e1 + t344 * t240 / 0.4e1 + t527 * t564 + pkin(8) * t402 + t368 + (t192 * t527 + t593 * t393) * t569 + (-t569 * t602 - t561) * t276 + t399 * t223 + t398 * t220;
t387 = mrSges(6,3) * (t340 / 0.2e1 + t339 / 0.2e1);
t490 = t378 * mrSges(7,3);
t591 = Ifges(5,4) + t422 + t506 / 0.2e1;
t437 = t341 ^ 2 + t342 ^ 2;
t215 = Ifges(7,4) * t223;
t113 = -Ifges(7,1) * t220 + t309 * Ifges(7,5) + t215;
t133 = Ifges(7,2) * t220 + t215;
t509 = Ifges(7,4) * t220;
t111 = Ifges(7,2) * t223 + t309 * Ifges(7,6) - t509;
t134 = Ifges(7,1) * t223 + t509;
t386 = (t111 / 0.4e1 - t134 / 0.4e1) * t378;
t585 = (t133 / 0.4e1 + t113 / 0.4e1) * t583 - t386;
t467 = t276 * t220;
t468 = t393 * t223;
t584 = t467 / 0.2e1 - t468 / 0.2e1;
t382 = -t152 * t344 + t153 * t346;
t510 = Ifges(6,4) * t344;
t321 = Ifges(6,2) * t346 + t510;
t324 = Ifges(6,1) * t346 - t510;
t581 = -t324 / 0.4e1 + t321 / 0.4e1;
t248 = t309 * mrSges(6,1) + mrSges(6,3) * t462;
t446 = t346 * t248;
t401 = -t446 / 0.2e1;
t580 = t542 * t610 + t401;
t579 = t166 * t540 + t378 * t606;
t443 = t502 / 0.2e1 + t504 / 0.2e1;
t501 = t608 * mrSges(7,2);
t503 = t590 * mrSges(7,1);
t444 = t503 / 0.2e1 + t501 / 0.2e1;
t574 = t309 ^ 2;
t573 = 2 * qJD(4);
t572 = -m(4) / 0.2e1;
t571 = m(6) / 0.2e1;
t570 = -m(7) / 0.2e1;
t562 = t166 / 0.2e1;
t560 = -t220 / 0.2e1;
t557 = t223 / 0.2e1;
t239 = t308 * t321;
t555 = t239 / 0.4e1;
t258 = t491 - t493;
t553 = -t258 / 0.2e1;
t551 = t260 / 0.2e1;
t548 = t262 / 0.2e1;
t546 = -t393 / 0.2e1;
t544 = -t276 / 0.2e1;
t538 = t378 / 0.2e1;
t535 = t344 / 0.2e1;
t534 = -t346 / 0.2e1;
t533 = t346 / 0.2e1;
t524 = t65 * mrSges(7,2);
t523 = t66 * mrSges(7,1);
t520 = t73 * mrSges(7,1);
t519 = t74 * mrSges(7,2);
t514 = m(7) * qJD(2);
t513 = m(7) * qJD(3);
t110 = -Ifges(7,4) * t608 + Ifges(7,2) * t590 - Ifges(7,6) * t308;
t112 = -Ifges(7,1) * t608 + Ifges(7,4) * t590 - Ifges(7,5) * t308;
t131 = -t501 - t503;
t188 = -Ifges(6,6) * t308 + t309 * t384;
t190 = -Ifges(6,5) * t308 - t309 * t324;
t193 = -pkin(5) * t460 + t253;
t237 = t320 * t309;
t245 = mrSges(6,2) * t308 + mrSges(6,3) * t460;
t247 = -mrSges(6,1) * t308 + mrSges(6,3) * t459;
t301 = t309 * mrSges(5,2);
t495 = t309 * Ifges(6,5);
t191 = -t308 * t324 + t495;
t449 = t346 * t191;
t3 = t153 * t246 + t140 * t247 + t152 * t248 - t252 * t237 - t253 * t464 + t141 * t245 + t113 * t556 + t111 * t558 + t110 * t557 + t112 * t560 + t192 * t131 + t193 * t132 + t70 * t166 + t65 * t167 + t69 * t168 + t66 * t165 - t327 * t301 + m(6) * (t140 * t152 + t141 * t153 + t252 * t253) + m(7) * (t192 * t193 + t65 * t69 + t66 * t70) + (-t449 / 0.2e1 + t454 / 0.2e1 + t380 + t591 * t309) * t309 + (-t327 * mrSges(5,1) + Ifges(7,5) * t220 / 0.2e1 - Ifges(7,6) * t223 / 0.2e1 + t190 * t534 + t188 * t535 - t591 * t308 + (-Ifges(6,3) - Ifges(7,3) - Ifges(5,2) + Ifges(5,1)) * t309) * t308;
t496 = t3 * qJD(1);
t492 = t583 * mrSges(7,3);
t442 = Ifges(7,5) * t223 + Ifges(7,6) * t220;
t485 = t65 * t223;
t358 = t192 * t130 + (t133 / 0.2e1 + t113 / 0.2e1) * t223 + (t66 * mrSges(7,3) + t111 / 0.2e1 - t134 / 0.2e1) * t220 - mrSges(7,3) * t485 + t309 * t442 / 0.2e1;
t4 = t74 * t166 + m(7) * (t65 * t73 + t66 * t74) + t73 * t168 + t140 * t246 - t141 * t248 + t252 * t236 + ((-t140 * mrSges(6,3) + t495 / 0.2e1 + t191 / 0.2e1 + t239 / 0.2e1) * t344 + (t141 * mrSges(6,3) + t494 / 0.2e1 - t240 / 0.2e1 + t189 / 0.2e1 + (-m(7) * t192 - t132) * pkin(5)) * t346) * t308 + t358;
t486 = t4 * qJD(1);
t9 = t65 * t166 - t66 * t168 + t358;
t484 = t9 * qJD(1);
t356 = (-t590 ^ 2 - t608 ^ 2 - t272) * t569 + (-t436 * t574 - t272) * t571 + m(5) * (-t272 - t574) / 0.2e1 + t437 * t572;
t365 = t572 - m(5) / 0.2e1 - m(6) * t436 / 0.2e1 + (t378 ^ 2 + t583 ^ 2) * t570;
t41 = t356 + t365;
t481 = qJD(1) * t41;
t10 = (t402 + t580) * t309 + t599 + t612;
t480 = t10 * qJD(1);
t413 = t487 / 0.2e1;
t360 = (-t343 * t608 + t531 * t590) * t433 + mrSges(6,1) * t403 + t309 * t413 + t444;
t447 = t346 * t246;
t451 = t344 * t248;
t379 = t451 / 0.2e1 - t447 / 0.2e1;
t363 = (-t378 * t602 + t583 * t593) * t569 - t379;
t414 = -t490 / 0.2e1;
t417 = t492 / 0.2e1;
t12 = t220 * t414 + t223 * t417 + t360 - t363 - t579;
t479 = t12 * qJD(1);
t15 = t355 - t586;
t478 = t15 * qJD(1);
t383 = t140 * t344 - t141 * t346;
t469 = t252 * t308;
t18 = -t608 * t166 + t590 * t168 + (mrSges(5,3) * t309 - t447 + t451) * t309 + (mrSges(5,3) * t308 - t132 + t464) * t308 + m(7) * (-t192 * t308 + t590 * t65 - t608 * t66) + m(6) * (t309 * t383 - t469) + m(5) * (-t253 * t309 - t469) + (m(4) * t518 + mrSges(4,3)) * t437;
t475 = t18 * qJD(1);
t415 = t490 / 0.2e1;
t418 = -t492 / 0.2e1;
t388 = t220 * t415 + t223 * t418 + t579;
t20 = t388 - t444;
t473 = t20 * qJD(1);
t471 = t220 * t343;
t302 = t309 * mrSges(5,1);
t29 = t341 * mrSges(4,1) + t342 * mrSges(4,2) - t308 * mrSges(5,2) + t378 * t166 + t583 * t168 + t453 + t446 + mrSges(3,3) + t302 + (m(4) + m(3)) * qJ(2) + m(7) * (t378 * t66 + t583 * t65) + m(6) * (t140 * t346 + t141 * t344) + m(5) * t327;
t466 = t29 * qJD(1);
t465 = t308 * t257;
t461 = t309 * t319;
t452 = t344 * t247;
t450 = t344 * t321;
t448 = t346 * t245;
t434 = mrSges(7,3) * t528;
t98 = -t220 * t378 + t223 * t583;
t432 = t98 * qJD(4) * t569;
t427 = -t65 / 0.2e1 + t74 / 0.2e1;
t426 = t66 / 0.2e1 + t73 / 0.2e1;
t412 = t465 / 0.2e1 - (t414 + t415) * t608;
t410 = t531 * t223;
t400 = t553 - t319 / 0.2e1;
t394 = t436 * t309;
t392 = mrSges(7,3) * t430;
t376 = m(7) * (t276 * t378 + t393 * t583);
t353 = (-t378 * t558 - t540 * t608) * mrSges(7,3) - t309 * t387 + (-pkin(8) * t394 + t530) * t571 + (-t276 * t608 - t308 * t331 + t393 * t590) * t569;
t354 = (t152 * t346 + t153 * t344) * t571 + (t378 * t70 + t583 * t69) * t569 + t167 * t540 + t165 * t538 + t245 * t535 + t247 * t533;
t19 = t301 + (mrSges(5,1) + t400) * t308 + t353 - t354;
t375 = -t19 * qJD(1) + t98 * t514 / 0.2e1;
t361 = (t531 * t562 + t343 * t606 + (t471 / 0.2e1 - t410 / 0.2e1) * mrSges(7,3)) * pkin(5);
t17 = -mrSges(7,1) * t426 + mrSges(7,2) * t427 + t361;
t318 = (mrSges(7,1) * t343 + mrSges(7,2) * t531) * pkin(5);
t59 = (t546 + t545) * mrSges(7,2) + (t544 + t605) * mrSges(7,1);
t372 = -t17 * qJD(1) - t59 * qJD(4) + t318 * qJD(5);
t271 = t308 * t309;
t43 = m(7) * (t271 - t470 - t472) + m(6) * (-t308 * t394 + t271);
t349 = (t131 / 0.2e1 - t237 / 0.2e1 + t379) * t308 + (t564 + t448 / 0.2e1 - t452 / 0.2e1 + t604) * t309 + ((t252 + t382) * t309 + (t253 + t383) * t308) * t571 + (t192 * t309 + t308 * t193 - t66 * t220 - t590 * t69 + t608 * t70 + t485) * t569 + t608 * t607 + t166 * t560 + t167 * t559 + t168 * t557;
t8 = -t376 / 0.2e1 + t349;
t371 = -t8 * qJD(1) - t43 * qJD(2) - t513 * t98 / 0.2e1;
t1 = t346 * t555 + t449 / 0.4e1 + t74 * t417 + t65 * t418 + (t133 + t113) * t583 / 0.4e1 + ((-t331 * t569 + t553) * pkin(5) + t581) * t462 + t593 * t414 + (t323 - t384) * t463 / 0.4e1 - t386 + t584 * mrSges(7,3) + t580 * pkin(8) + t595 - t613;
t359 = (t410 - t471) * t433 + t463 * t567 + t308 * t413 + t443;
t362 = t431 * t570 + t604;
t31 = -t465 / 0.2e1 + t359 + t362;
t42 = t331 * t257 - (-t261 / 0.2e1 + t551) * t378 + (t548 + t259 / 0.2e1) * t583;
t33 = -t450 / 0.2e1 + t324 * t535 - pkin(4) * t320 + (-t384 / 0.2e1 + t323 / 0.2e1) * t346 + t42 + (m(7) * t331 + t258) * t527;
t367 = -t1 * qJD(1) + t31 * qJD(2) - t33 * qJD(4);
t34 = t412 - t443;
t350 = (mrSges(7,3) * t546 + t399) * t223 + (mrSges(7,3) * t605 + t398) * t220 + t393 * t562 + t168 * t544 + t368 + t585;
t5 = t350 + t578;
t366 = -t5 * qJD(1) - t34 * qJD(2) - t42 * qJD(4);
t311 = t318 * qJD(6);
t40 = t356 - t365;
t35 = t412 + t443;
t32 = t359 - t362 + t412;
t22 = t308 * t400 + t353 + t354;
t21 = t388 + t444;
t16 = t355 + t586;
t14 = -t524 / 0.2e1 - t523 / 0.2e1 - t519 / 0.2e1 + t520 / 0.2e1 + t361 + t442;
t13 = t360 + t363 + t388;
t11 = (t308 * t387 + t401 + t402) * t309 + t599 - t612;
t7 = t376 / 0.2e1 + t349;
t6 = t350 - t578;
t2 = (-t378 * t426 + t427 * t583 + t584) * mrSges(7,3) + (t555 + t191 / 0.4e1 - pkin(8) * t248 / 0.2e1) * t346 + ((t323 / 0.4e1 - t384 / 0.4e1) * t344 + pkin(8) * t387 + ((t331 * t570 + t553) * pkin(5) + t581) * t346) * t308 + t585 + t595 + t613;
t23 = [qJD(2) * t29 + qJD(3) * t18 + qJD(4) * t3 + qJD(5) * t4 + qJD(6) * t9, t40 * qJD(3) + t7 * qJD(4) + t11 * qJD(5) + t16 * qJD(6) + t611 * t514 + t466, t40 * qJD(2) + t22 * qJD(4) + t13 * qJD(5) + t21 * qJD(6) - t611 * t513 + t475, t496 + t7 * qJD(2) + t22 * qJD(3) + t2 * qJD(5) + t6 * qJD(6) + ((-pkin(4) * t253 + pkin(8) * t382) * t571 + (t193 * t331 + t276 * t70 + t393 * t69) * t569) * t573 + (t188 * t533 + t190 * t535 + t331 * t131 + t112 * t538 + t110 * t540 + Ifges(5,6) * t308 + t393 * t167 + t276 * t165 + t590 * t551 - t608 * t548 + t193 * t258 + t252 * mrSges(5,2) + pkin(4) * t237 - t69 * t490 + t70 * t492 + (-Ifges(5,5) + t323 * t534 + t450 / 0.2e1) * t309 + (Ifges(6,5) * t344 + Ifges(7,5) * t378 + Ifges(6,6) * t346 + Ifges(7,6) * t583) * t543 + (t319 - mrSges(5,1)) * t253 + (t448 - t452) * pkin(8) + t382 * mrSges(6,3)) * qJD(4), t486 + t11 * qJD(2) + t13 * qJD(3) + t2 * qJD(4) + (Ifges(6,5) * t463 + Ifges(6,6) * t462 - t519 + (t343 * t74 + t531 * t73) * t568 + t220 * t434 - t223 * t392 + t520 - t140 * mrSges(6,2) - t141 * mrSges(6,1) + t442) * qJD(5) + t14 * qJD(6), t484 + t16 * qJD(2) + t21 * qJD(3) + t6 * qJD(4) + t14 * qJD(5) + (t442 - t523 - t524) * qJD(6); qJD(3) * t41 + qJD(4) * t8 + qJD(5) * t10 + qJD(6) * t15 - t466, t43 * qJD(4), t432 + t481, t32 * qJD(5) + t35 * qJD(6) + ((t331 * t309 - t467 + t468) * t569 + (-t436 * t525 - t529) * t571) * t573 - t371 + (-t220 * t492 - t223 * t490 + t309 * t258 - t302 + t461 + (mrSges(5,2) - t610) * t308) * qJD(4), t480 + t32 * qJD(4) + (m(7) * (-t430 * t608 - t528 * t590) + t461 + t396) * qJD(5) + t597, t35 * qJD(4) + qJD(5) * t396 + t478 + t597; -qJD(2) * t41 - qJD(4) * t19 - qJD(5) * t12 + qJD(6) * t20 - t475, t432 - t481, 0, t375, -t479 - t596 + (-t257 - t320 + (-t378 * t430 + t528 * t583) * m(7)) * qJD(5), -qJD(5) * t257 + t473 - t596; -qJD(2) * t8 + qJD(3) * t19 + qJD(5) * t1 + qJD(6) * t5 - t496, -t31 * qJD(5) + t34 * qJD(6) + t371, -t375, qJD(5) * t33 + qJD(6) * t42 (-t378 * t434 - t583 * t392 + (-t276 * t531 + t343 * t393) * t568 + t395 + t319 * pkin(8) + t49) * qJD(5) + t609 - t367, t49 * qJD(5) - t366 + t609; -qJD(2) * t10 + qJD(3) * t12 - qJD(4) * t1 + qJD(6) * t17 - t486, t31 * qJD(4) - t480, t479, t59 * qJD(6) + t367, -t311, -t311 - t372; -qJD(2) * t15 - qJD(3) * t20 - qJD(4) * t5 - qJD(5) * t17 - t484, -t34 * qJD(4) - t478, -t473, -t59 * qJD(5) + t366, t372, 0;];
Cq  = t23;
