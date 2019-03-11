% Calculate vector of inverse dynamics joint torques for
% S6PRRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-03-08 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRPP3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPP3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:54:28
% EndTime: 2019-03-08 22:55:04
% DurationCPUTime: 24.23s
% Computational Cost: add. (5185->738), mult. (11951->922), div. (0->0), fcn. (8594->10), ass. (0->345)
t522 = -Ifges(7,5) - Ifges(5,5);
t541 = Ifges(6,4) + t522;
t523 = -Ifges(6,5) - Ifges(7,4);
t548 = Ifges(5,6) + t523;
t545 = Ifges(5,1) + Ifges(7,3);
t544 = Ifges(7,2) + Ifges(6,3);
t549 = -Ifges(5,4) + Ifges(7,6);
t250 = sin(qJ(4));
t253 = cos(qJ(4));
t251 = sin(qJ(3));
t254 = cos(qJ(3));
t390 = qJD(2) * qJD(3);
t278 = qJDD(2) * t251 + t254 * t390;
t393 = qJD(4) * t251;
t366 = t250 * t393;
t391 = t253 * qJD(3);
t90 = qJD(2) * t366 - qJD(4) * t391 - t250 * qJDD(3) - t253 * t278;
t474 = -t90 / 0.2e1;
t401 = qJD(2) * t251;
t192 = qJD(3) * t250 + t253 * t401;
t395 = qJD(4) * t192;
t91 = -t253 * qJDD(3) + t250 * t278 + t395;
t472 = -t91 / 0.2e1;
t471 = t91 / 0.2e1;
t198 = qJDD(2) * t254 - t251 * t390;
t189 = qJDD(4) - t198;
t468 = t189 / 0.2e1;
t191 = t250 * t401 - t391;
t467 = -t191 / 0.2e1;
t464 = t192 / 0.2e1;
t399 = qJD(2) * t254;
t232 = -qJD(4) + t399;
t462 = t232 / 0.2e1;
t546 = qJD(3) / 0.2e1;
t438 = Ifges(7,6) * t253;
t440 = Ifges(6,6) * t253;
t496 = t544 * t250 + t438 - t440;
t439 = Ifges(7,6) * t250;
t443 = Ifges(5,4) * t250;
t495 = t545 * t253 + t439 - t443;
t426 = qJ(5) * t250;
t457 = pkin(4) * t253;
t543 = t426 + t457;
t542 = Ifges(6,1) + Ifges(7,1) + Ifges(5,3);
t485 = -t548 * t250 - t253 * t541;
t458 = pkin(3) * t254;
t206 = -pkin(9) * t251 - pkin(2) - t458;
t277 = t206 * qJD(2);
t255 = cos(qJ(2));
t247 = sin(pkin(6));
t406 = qJD(1) * t247;
t371 = t255 * t406;
t259 = t277 - t371;
t252 = sin(qJ(2));
t372 = t252 * t406;
t199 = qJD(2) * pkin(8) + t372;
t248 = cos(pkin(6));
t388 = qJDD(1) * t248;
t398 = qJD(3) * t251;
t402 = qJD(2) * t247;
t359 = qJD(1) * t402;
t217 = t255 * t359;
t389 = qJDD(1) * t247;
t159 = t252 * t389 + t217;
t405 = qJD(1) * t248;
t538 = qJDD(2) * pkin(8) + qJD(3) * t405 + t159;
t35 = -t199 * t398 + t251 * t388 + t254 * t538;
t539 = qJDD(3) * pkin(9) + qJD(4) * t259 + t35;
t315 = mrSges(6,2) * t253 - mrSges(6,3) * t250;
t317 = mrSges(5,1) * t253 - mrSges(5,2) * t250;
t362 = m(7) * qJ(6) + mrSges(7,3);
t434 = t250 * mrSges(7,2);
t537 = t253 * t362 + mrSges(4,1) - t315 + t317 + t434;
t136 = t254 * t199 + t251 * t405;
t394 = qJD(4) * t250;
t536 = -qJD(5) * t250 - t136 + (-t250 * t399 + t394) * pkin(4);
t392 = qJD(4) * t253;
t397 = qJD(3) * t254;
t535 = t250 * t397 + t251 * t392;
t121 = qJD(3) * pkin(9) + t136;
t39 = t121 * t250 - t253 * t259;
t288 = pkin(5) * t192 + t39;
t534 = qJD(5) + t288;
t533 = t191 * pkin(5) - qJD(6);
t238 = Ifges(4,4) * t399;
t188 = Ifges(5,4) * t191;
t437 = t191 * Ifges(7,6);
t518 = t192 * t545 + t232 * t522 - t188 + t437;
t186 = Ifges(7,6) * t192;
t435 = t192 * Ifges(6,6);
t519 = t191 * t544 + t232 * t523 + t186 - t435;
t532 = Ifges(4,1) * t401 + Ifges(4,5) * qJD(3) + t250 * t519 + t253 * t518 + t238;
t445 = Ifges(4,4) * t251;
t515 = t254 * Ifges(4,2);
t310 = t445 + t515;
t531 = Ifges(4,6) * t546 + qJD(2) * t310 / 0.2e1 + t542 * t462 + t541 * t464 - t548 * t467;
t216 = t252 * t359;
t158 = t255 * t389 - t216;
t140 = -qJDD(2) * pkin(2) - t158;
t59 = -t198 * pkin(3) - pkin(9) * t278 + t140;
t7 = -t121 * t392 - t250 * t539 + t253 * t59;
t281 = qJDD(5) - t7;
t449 = pkin(4) + qJ(6);
t1 = -pkin(5) * t90 + qJD(6) * t232 - t189 * t449 + t281;
t530 = t549 * t471 + mrSges(7,1) * t1 + Ifges(6,6) * t472 - t541 * t189 / 0.2e1 + (t545 + Ifges(6,2)) * t474;
t6 = -t121 * t394 + t250 * t59 + t253 * t539;
t4 = -qJ(5) * t189 + qJD(5) * t232 - t6;
t2 = -pkin(5) * t91 + qJDD(6) - t4;
t529 = t2 * mrSges(7,1) + (Ifges(6,6) - t549) * t474 + (Ifges(5,2) + t544) * t472 + t548 * t468;
t470 = -m(6) - m(7);
t528 = t198 / 0.2e1;
t527 = t278 / 0.2e1;
t525 = pkin(8) * (t251 * t391 + t254 * t394);
t524 = mrSges(3,2) - mrSges(4,3);
t50 = -t90 * mrSges(7,1) - t189 * mrSges(7,3);
t53 = -t90 * mrSges(6,1) + t189 * mrSges(6,2);
t520 = t50 + t53;
t428 = cos(pkin(10));
t336 = t428 * t252;
t246 = sin(pkin(10));
t421 = t246 * t255;
t165 = t248 * t336 + t421;
t337 = t247 * t428;
t105 = -t165 * t251 - t254 * t337;
t517 = t543 * t105;
t335 = t428 * t255;
t422 = t246 * t252;
t167 = -t248 * t422 + t335;
t419 = t247 * t254;
t107 = -t167 * t251 + t246 * t419;
t516 = t543 * t107;
t19 = mrSges(5,1) * t91 - mrSges(5,2) * t90;
t514 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t278 + t19;
t318 = mrSges(4,1) * t254 - mrSges(4,2) * t251;
t513 = t318 + mrSges(3,1);
t469 = pkin(5) + pkin(9);
t215 = t469 * t253;
t414 = t253 * t254;
t385 = pkin(5) * t414;
t135 = -t251 * t199 + t254 * t405;
t321 = pkin(3) * t251 - pkin(9) * t254;
t194 = t321 * qJD(2);
t70 = -t250 * t135 + t194 * t253;
t512 = -(-t251 * t449 + t385) * qJD(2) + t70 + qJD(4) * t215;
t511 = -m(7) * t469 - mrSges(7,1) + mrSges(4,2);
t416 = t250 * t254;
t386 = pkin(5) * t416;
t71 = t253 * t135 + t250 * t194;
t510 = -(qJ(5) * t251 - t386) * qJD(2) - t71 - t469 * t394;
t425 = qJ(5) * t253;
t291 = qJ(6) * t250 - t425;
t279 = t291 * t254;
t509 = -qJD(2) * t279 + qJD(4) * t291 - qJD(6) * t253 + t536;
t508 = -qJ(5) * t392 + t399 * t425 + t536;
t507 = -qJD(5) - t39;
t329 = t250 * t371;
t40 = t253 * t121 + t250 * t277 - t329;
t506 = -t40 + t533;
t101 = -mrSges(7,2) * t192 + mrSges(7,3) * t191;
t104 = -mrSges(6,2) * t191 - mrSges(6,3) * t192;
t412 = t101 + t104;
t420 = t247 * t252;
t169 = -t248 * t254 + t251 * t420;
t369 = t255 * t402;
t110 = -qJD(3) * t169 + t254 * t369;
t418 = t247 * t255;
t505 = -qJD(4) * t418 + t110;
t447 = mrSges(5,3) * t191;
t129 = mrSges(5,2) * t232 - t447;
t132 = mrSges(6,1) * t191 + mrSges(6,3) * t232;
t133 = -mrSges(7,1) * t191 - mrSges(7,2) * t232;
t411 = -t132 + t133;
t504 = t129 + t411;
t131 = mrSges(7,1) * t192 + mrSges(7,3) * t232;
t446 = mrSges(5,3) * t192;
t130 = -mrSges(5,1) * t232 - t446;
t134 = mrSges(6,1) * t192 - mrSges(6,2) * t232;
t410 = t134 - t130;
t503 = t131 + t410;
t502 = t543 * t169;
t378 = mrSges(4,3) * t401;
t501 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t191 + mrSges(5,2) * t192 + t378;
t500 = -t251 * t522 + t254 * t495;
t499 = -t251 * t523 + t254 * t496;
t292 = mrSges(5,1) - mrSges(6,2) + t362;
t498 = -pkin(4) * t470 + t292;
t441 = Ifges(6,6) * t250;
t497 = t253 * t544 - t439 + t441;
t36 = -t199 * t397 - t251 * t538 + t254 * t388;
t494 = -t251 * t36 + t254 * t35;
t187 = Ifges(6,6) * t191;
t77 = -t232 * Ifges(6,4) - t192 * Ifges(6,2) + t187;
t436 = t192 * Ifges(5,4);
t78 = -t191 * Ifges(5,2) - t232 * Ifges(5,6) + t436;
t493 = t250 * t78 + t253 * t77;
t492 = -t250 * t7 + t253 * t6;
t5 = -pkin(4) * t189 + t281;
t491 = t250 * t5 - t253 * t4;
t490 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t487 = t251 * t542 + t254 * t485;
t120 = -qJD(3) * pkin(3) - t135;
t260 = -qJ(5) * t192 + t120;
t27 = t191 * t449 + t260;
t313 = -mrSges(7,2) * t253 + mrSges(7,3) * t250;
t314 = -mrSges(6,2) * t250 - mrSges(6,3) * t253;
t316 = mrSges(5,1) * t250 + mrSges(5,2) * t253;
t41 = pkin(4) * t191 + t260;
t486 = -t120 * t316 - t27 * t313 - t41 * t314;
t484 = t189 * t542 + t541 * t90 - t548 * t91;
t480 = -m(5) * t120 - t501;
t479 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1) + mrSges(5,3);
t106 = t165 * t254 - t251 * t337;
t108 = t246 * t247 * t251 + t167 * t254;
t170 = t248 * t251 + t252 * t419;
t477 = -g(1) * t108 - g(2) * t106 - g(3) * t170;
t476 = -t7 * mrSges(5,1) + t6 * mrSges(5,2) - t5 * mrSges(6,2) - t2 * mrSges(7,2) + t4 * mrSges(6,3) + t1 * mrSges(7,3);
t256 = qJD(2) ^ 2;
t473 = t90 / 0.2e1;
t466 = t191 / 0.2e1;
t465 = -t192 / 0.2e1;
t463 = -t232 / 0.2e1;
t244 = t251 * pkin(8);
t51 = mrSges(6,1) * t91 - mrSges(6,3) * t189;
t52 = -t91 * mrSges(7,1) + t189 * mrSges(7,2);
t448 = -t51 + t52;
t444 = Ifges(4,4) * t254;
t442 = Ifges(5,4) * t253;
t34 = -qJDD(3) * pkin(3) - t36;
t432 = t251 * t34;
t427 = qJ(5) * t191;
t164 = -t248 * t335 + t422;
t424 = t164 * t251;
t166 = t248 * t421 + t336;
t423 = t166 * t251;
t417 = t250 * t251;
t415 = t251 * t253;
t413 = t254 * t255;
t197 = t321 * qJD(3);
t409 = t250 * t197 + t206 * t392;
t408 = pkin(2) * t418 + pkin(8) * t420;
t407 = pkin(4) * t417 + t244;
t235 = pkin(8) * t414;
t143 = t250 * t206 + t235;
t400 = qJD(2) * t252;
t234 = pkin(8) * t416;
t381 = pkin(8) * t398;
t380 = pkin(9) * t394;
t379 = pkin(9) * t392;
t240 = pkin(8) * t397;
t377 = mrSges(4,3) * t399;
t376 = t251 * t418;
t375 = t250 * t418;
t373 = -pkin(8) * t250 - pkin(4);
t370 = t247 * t400;
t97 = t105 * pkin(3);
t361 = pkin(9) * t106 + t97;
t98 = t107 * pkin(3);
t360 = pkin(9) * t108 + t98;
t345 = -t394 / 0.2e1;
t344 = t394 / 0.2e1;
t343 = -t392 / 0.2e1;
t342 = t392 / 0.2e1;
t341 = -pkin(3) - t426;
t340 = -t164 * pkin(2) + pkin(8) * t165;
t339 = -t166 * pkin(2) + pkin(8) * t167;
t162 = t169 * pkin(3);
t338 = pkin(9) * t170 - t162;
t334 = t390 / 0.2e1;
t142 = t206 * t253 - t234;
t332 = t247 * pkin(3) * t413 + pkin(9) * t376 + t408;
t331 = pkin(4) * t535 + qJ(5) * t366 + t240;
t123 = qJ(5) * t254 - t143;
t320 = qJD(4) * t235 - t197 * t253 + t206 * t394;
t319 = -qJD(5) * t254 + t409;
t311 = Ifges(5,1) * t250 + t442;
t309 = -Ifges(5,2) * t250 + t442;
t308 = Ifges(5,2) * t253 + t443;
t306 = Ifges(6,4) * t250 + Ifges(6,5) * t253;
t305 = Ifges(7,4) * t253 - Ifges(7,5) * t250;
t303 = Ifges(4,5) * t254 - Ifges(4,6) * t251;
t301 = Ifges(5,5) * t250 + Ifges(5,6) * t253;
t300 = -Ifges(6,2) * t253 + t441;
t299 = Ifges(6,2) * t250 + t440;
t294 = -Ifges(7,3) * t250 + t438;
t290 = -pkin(9) * t424 - t164 * t458 + t340;
t289 = -pkin(9) * t423 - t166 * t458 + t339;
t287 = qJ(5) * t470 + t490;
t200 = -qJD(2) * pkin(2) - t371;
t283 = t200 * (mrSges(4,1) * t251 + mrSges(4,2) * t254);
t282 = t251 * (Ifges(4,1) * t254 - t445);
t111 = t170 * t250 + t253 * t418;
t44 = t106 * t250 - t164 * t253;
t46 = t108 * t250 - t166 * t253;
t280 = -g(1) * t46 - g(2) * t44 - g(3) * t111;
t139 = (t250 * t252 + t253 * t413) * t247;
t267 = Ifges(6,4) * t251 + t254 * t300;
t262 = Ifges(5,6) * t251 + t254 * t309;
t257 = qJ(5) * t90 - qJD(5) * t192 + t34;
t29 = t232 * qJ(5) - t40;
t245 = t254 * pkin(4);
t214 = t469 * t250;
t210 = -qJD(3) * mrSges(4,2) + t377;
t201 = t341 - t457;
t193 = t318 * qJD(2);
t181 = -t253 * t449 + t341;
t160 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t198;
t155 = -qJ(5) * t415 + t407;
t138 = -t253 * t420 + t254 * t375;
t124 = -t142 + t245;
t119 = t251 * t291 + t407;
t114 = -t198 * mrSges(4,1) + mrSges(4,2) * t278;
t112 = t170 * t253 - t375;
t109 = qJD(3) * t170 + t251 * t369;
t102 = pkin(4) * t192 + t427;
t96 = -pkin(5) * t417 - t123;
t89 = qJ(6) * t254 + t234 + t245 + (pkin(5) * t251 - t206) * t253;
t69 = -t166 * t414 + t167 * t250;
t68 = -t166 * t416 - t167 * t253;
t67 = -t164 * t414 + t165 * t250;
t66 = -t164 * t416 - t165 * t253;
t63 = t192 * t449 + t427;
t62 = t250 * t381 - t320;
t61 = t409 - t525;
t60 = (-qJ(5) * t397 - qJD(5) * t251) * t253 + t331;
t58 = -pkin(4) * t401 - t70;
t57 = -qJ(5) * t401 - t71;
t54 = t373 * t398 + t320;
t49 = -mrSges(5,2) * t189 - mrSges(5,3) * t91;
t48 = mrSges(5,1) * t189 + mrSges(5,3) * t90;
t47 = t108 * t253 + t166 * t250;
t45 = t106 * t253 + t164 * t250;
t37 = -qJ(5) * t398 - t319 + t525;
t28 = pkin(4) * t232 - t507;
t24 = qJD(3) * t279 + (qJD(6) * t250 + (qJ(6) * qJD(4) - qJD(5)) * t253) * t251 + t331;
t23 = (-pkin(5) * t415 - t234) * qJD(4) + (-t386 + (-pkin(8) * t253 + qJ(5)) * t251) * qJD(3) + t319;
t22 = -t170 * t394 + t250 * t370 + t253 * t505;
t21 = t170 * t392 + t250 * t505 - t253 * t370;
t20 = mrSges(7,2) * t90 + mrSges(7,3) * t91;
t18 = -mrSges(6,2) * t91 + mrSges(6,3) * t90;
t17 = -pkin(5) * t366 + qJD(6) * t254 + (t385 + (-qJ(6) + t373) * t251) * qJD(3) + t320;
t16 = -t29 - t533;
t15 = t232 * t449 + t534;
t8 = pkin(4) * t91 + t257;
t3 = qJD(6) * t191 + t449 * t91 + t257;
t9 = [m(2) * qJDD(1) + t110 * t210 + t170 * t160 + t504 * t22 + t503 * t21 + (t49 + t448) * t112 + (-t48 + t520) * t111 + (t18 + t20 + t514) * t169 + (t412 + t501) * t109 + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t256 - t114) * t255 + (-mrSges(3,1) * t256 - mrSges(3,2) * qJDD(2) - qJD(2) * t193) * t252) * t247 + (-m(2) - m(3) - m(4) - m(5) + t470) * g(3) + m(7) * (t1 * t111 + t109 * t27 + t112 * t2 + t15 * t21 + t16 * t22 + t169 * t3) + m(6) * (t109 * t41 + t111 * t5 - t112 * t4 + t169 * t8 + t21 * t28 - t22 * t29) + m(5) * (t109 * t120 - t111 * t7 + t112 * t6 + t169 * t34 + t21 * t39 + t22 * t40) + m(4) * (-t109 * t135 + t110 * t136 - t169 * t36 + t170 * t35 + (-t140 * t255 + t200 * t400) * t247) + m(3) * (qJDD(1) * t248 ^ 2 + (t158 * t255 + t159 * t252) * t247); (t216 + t158) * mrSges(3,1) + (Ifges(4,4) * t527 - Ifges(6,4) * t473 + Ifges(4,2) * t528 - Ifges(5,6) * t472 + pkin(8) * t160 - t210 * t371 + t444 * t334 + t471 * t523 + t474 * t522 + t476) * t254 + (Ifges(4,1) * t278 + Ifges(4,4) * t528 + t3 * t313 + t300 * t473 + t309 * t472 + t8 * t314 - t334 * t515 + t78 * t343 + t77 * t344 + t471 * t496 + t474 * t495) * t251 + (-m(4) * t408 - m(5) * t332 + t470 * (t139 * pkin(4) + qJ(5) * t138 + t332) + (t252 * t524 - t255 * t513) * t247 - t479 * t376 - t292 * t139 + t490 * t138) * g(3) + (-m(4) * t339 - m(5) * t289 + t470 * (t69 * pkin(4) + qJ(5) * t68 + t289) + t524 * t167 + t513 * t166 - t292 * t69 + t490 * t68 + t479 * t423) * g(1) + (-m(4) * t340 - m(5) * t290 + t470 * (t67 * pkin(4) + qJ(5) * t66 + t290) + t524 * t165 + t513 * t164 - t292 * t67 + t490 * t66 + t479 * t424) * g(2) + (t4 * mrSges(6,1) - t6 * mrSges(5,3) - t529) * t417 + t518 * t251 * t345 + t519 * t251 * t342 + (-t135 * t397 + t494) * mrSges(4,3) + (mrSges(6,1) * t28 + mrSges(7,1) * t15 + mrSges(5,2) * t120 - mrSges(7,2) * t27 + mrSges(5,3) * t39 - mrSges(6,3) * t41) * (t254 * t391 - t366) + (-m(5) * t40 + m(6) * t29 - m(7) * t16 - t504) * qJD(1) * t139 + (-m(6) * t41 - m(7) * t27 - t412 + t480) * t251 * t371 + (-t39 * mrSges(5,1) - t40 * mrSges(5,2) + t28 * mrSges(6,2) + t16 * mrSges(7,2) - t136 * mrSges(4,3) - t29 * mrSges(6,3) - t15 * mrSges(7,3) - t531) * t398 + t532 * t397 / 0.2e1 + m(5) * (t142 * t7 + t143 * t6 - t39 * t62 + t40 * t61 + (t120 * t397 + t432) * pkin(8)) + (t485 * t251 - t254 * t542) * t468 - t140 * t318 + t514 * t244 + (t217 - t159) * mrSges(3,2) + (t262 * t467 + t267 * t465 + t303 * t546 + t500 * t464 + t499 * t466 + t283) * qJD(3) + (mrSges(6,1) * t5 - mrSges(5,3) * t7 + t530) * t415 + (-m(5) * t39 - m(6) * t28 - m(7) * t15 - t503) * (-t253 * t372 + t254 * t329) + t193 * t372 - t493 * t397 / 0.2e1 + (-pkin(2) * t140 + ((-t135 * t254 - t136 * t251) * qJD(3) + t494) * pkin(8) - (t200 * t252 + (-t135 * t251 + t136 * t254) * t255) * t406) * m(4) + (t299 * t465 - t308 * t467 + t497 * t466 + (t294 - t311) * t464) * t393 + t501 * t240 + t282 * t334 - t484 * t254 / 0.2e1 + ((-t301 + t305 + t306) * t393 + t487 * qJD(3)) * t463 + m(7) * (t1 * t89 + t119 * t3 + t15 * t17 + t16 * t23 + t2 * t96 + t24 * t27) + m(6) * (t123 * t4 + t124 * t5 + t155 * t8 + t28 * t54 + t29 * t37 + t41 * t60) + Ifges(3,3) * qJDD(2) + (mrSges(5,1) * t120 + mrSges(6,1) * t29 - mrSges(7,1) * t16 - mrSges(6,2) * t41 - mrSges(5,3) * t40 + mrSges(7,3) * t27) * t535 + t89 * t50 + t96 * t52 + t24 * t101 + t60 * t104 - pkin(2) * t114 + t119 * t20 + t123 * t51 + t124 * t53 + t61 * t129 + t62 * t130 + t17 * t131 + t37 * t132 + t23 * t133 + t54 * t134 - t210 * t381 + t142 * t48 + t143 * t49 + t155 * t18 + t444 * t527 + t310 * t528 + t316 * t432 + qJDD(3) * (Ifges(4,5) * t251 + Ifges(4,6) * t254); (t463 * t485 - t486) * qJD(4) + (t378 + t480) * t136 + (t15 * t392 - t16 * t394) * mrSges(7,1) + (t496 / 0.2e1 - t309 / 0.2e1) * qJD(4) * t191 + (-t3 * mrSges(7,3) + t529) * t253 + t531 * t401 - (-Ifges(4,2) * t401 + t238 + t532) * t399 / 0.2e1 + (t377 - t210) * t135 - t34 * t317 + t8 * t315 + t518 * t342 + t519 * t344 + (t487 * t462 - t40 * (-mrSges(5,2) * t251 - mrSges(5,3) * t416) - t16 * (-mrSges(7,1) * t416 + mrSges(7,2) * t251) - t29 * (mrSges(6,1) * t416 - mrSges(6,3) * t251) - t283 + t39 * (mrSges(5,1) * t251 - mrSges(5,3) * t414) - t15 * (mrSges(7,1) * t414 - mrSges(7,3) * t251) - t28 * (mrSges(6,1) * t414 + mrSges(6,2) * t251) + (t262 / 0.2e1 - t499 / 0.2e1) * t191 + (t267 / 0.2e1 - t500 / 0.2e1) * t192) * qJD(2) + (-t379 - t70) * t130 + (-m(5) * t360 - m(7) * (t98 + t516) - m(6) * (t360 + t516) + t511 * t108 - t537 * t107) * g(1) + (-m(5) * t361 - m(7) * (t97 + t517) - m(6) * (t361 + t517) + t511 * t106 - t537 * t105) * g(2) + (-m(5) * t338 - m(6) * (t338 - t502) - m(7) * (-t162 - t502) + t511 * t170 + t537 * t169) * g(3) + t77 * t343 + t78 * t345 + (t201 * t8 - t28 * t58 - t29 * t57 + t508 * t41) * m(6) + (-pkin(3) * t34 + t39 * t70 - t40 * t71) * m(5) + ((-t48 + t53) * t250 + (-t51 + t49) * t253 + ((t250 * t29 + t253 * t28) * qJD(4) + t491) * m(6) + ((-t250 * t40 + t253 * t39) * qJD(4) + t492) * m(5)) * pkin(9) + (t28 * t392 + t29 * t394 + t477 + t491) * mrSges(6,1) - t3 * t434 - t256 * t282 / 0.2e1 + (t493 / 0.2e1 + t486) * t399 - t303 * t390 / 0.2e1 + (t380 - t57) * t132 + t508 * t104 + t509 * t101 + t510 * t133 + t512 * t131 + (t1 * t214 + t512 * t15 + t510 * t16 + t181 * t3 + t2 * t215 + t509 * t27) * m(7) + (t311 + t299) * t474 + (-t71 - t380) * t129 + (t379 - t58) * t134 + t530 * t250 + (-t306 / 0.2e1 - t305 / 0.2e1) * t189 + (t39 * t392 - t394 * t40 + t477 + t492) * mrSges(5,3) + (t308 + t497) * t472 + Ifges(4,3) * qJDD(3) + Ifges(4,5) * t278 + (t495 / 0.2e1 - t300 / 0.2e1) * t395 - pkin(3) * t19 - t35 * mrSges(4,2) + t36 * mrSges(4,1) + t181 * t20 + Ifges(4,6) * t198 + t201 * t18 + t214 * t50 + t215 * t52 + t301 * t468 + t294 * t473; -t476 + t484 + t288 * t133 - t449 * t50 + (qJ(5) * t2 - t1 * t449 + t506 * t15 + t16 * t534 - t27 * t63) * m(7) + (Ifges(6,2) * t191 + t435 + t78) * t464 + t448 * qJ(5) + (-Ifges(5,2) * t192 - t188 + t518) * t466 + (-t410 + t446) * t40 + (t129 - t132 + t447) * t39 + t411 * qJD(5) + (t191 * t541 - t192 * t548) * t462 + (t192 * t544 + t187 - t437 + t77) * t467 + (-t191 * t545 + t186 - t436 + t519) * t465 + t506 * t131 + (-pkin(4) * t5 - qJ(5) * t4 - t102 * t41 - t28 * t40 + t29 * t507) * m(6) + (t287 * t47 + t46 * t498) * g(1) + (t287 * t45 + t44 * t498) * g(2) + (t111 * t498 + t112 * t287) * g(3) + (t191 * t28 - t192 * t29) * mrSges(6,1) - pkin(4) * t53 + (t15 * t191 + t16 * t192) * mrSges(7,1) - t63 * t101 - t102 * t104 - t41 * (-mrSges(6,2) * t192 + mrSges(6,3) * t191) - t27 * (mrSges(7,2) * t191 + mrSges(7,3) * t192) - t120 * (mrSges(5,1) * t192 - mrSges(5,2) * t191); t412 * t192 + t411 * t232 + (t16 * t232 + t192 * t27 + t1 + t280) * m(7) + (t192 * t41 - t232 * t29 + t280 + t5) * m(6) + t520; -t191 * t101 - t232 * t131 + (-g(1) * t47 - g(2) * t45 - g(3) * t112 - t15 * t232 - t27 * t191 + t2) * m(7) + t52;];
tau  = t9;
