% Calculate vector of inverse dynamics joint torques for
% S6RPRRRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRP7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP7_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP7_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP7_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP7_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP7_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:18:40
% EndTime: 2019-03-09 06:19:31
% DurationCPUTime: 32.35s
% Computational Cost: add. (14883->746), mult. (35215->930), div. (0->0), fcn. (26710->14), ass. (0->336)
t281 = cos(qJ(4));
t263 = pkin(4) * t281 + pkin(3);
t271 = pkin(10) + qJ(3);
t266 = cos(t271);
t237 = t266 * t263;
t265 = sin(t271);
t284 = -pkin(9) - pkin(8);
t383 = t265 * t284;
t548 = mrSges(6,3) + mrSges(7,2);
t533 = t548 * t265;
t539 = m(6) + m(7);
t550 = -t533 - t539 * (t237 - t383);
t489 = -Ifges(6,4) + Ifges(7,5);
t549 = t489 + Ifges(7,5);
t273 = sin(pkin(10));
t274 = cos(pkin(10));
t278 = sin(qJ(3));
t282 = cos(qJ(3));
t511 = -t273 * t278 + t282 * t274;
t222 = t511 * qJD(3);
t228 = t273 * t282 + t274 * t278;
t164 = qJD(1) * t222 + qJDD(1) * t228;
t277 = sin(qJ(4));
t221 = t228 * qJD(1);
t301 = qJD(3) * t281 - t221 * t277;
t293 = t301 * qJD(4);
t117 = qJDD(3) * t277 + t164 * t281 + t293;
t179 = qJD(3) * t277 + t221 * t281;
t118 = -qJD(4) * t179 + qJDD(3) * t281 - t164 * t277;
t276 = sin(qJ(5));
t280 = cos(qJ(5));
t125 = t276 * t179 - t280 * t301;
t48 = -qJD(5) * t125 + t280 * t117 + t276 * t118;
t460 = t48 / 0.2e1;
t290 = t280 * t179 + t276 * t301;
t49 = qJD(5) * t290 + t276 * t117 - t280 * t118;
t458 = t49 / 0.2e1;
t223 = t228 * qJD(3);
t349 = qJDD(1) * t274;
t350 = qJDD(1) * t273;
t165 = -qJD(1) * t223 - t278 * t350 + t282 * t349;
t158 = qJDD(4) - t165;
t155 = qJDD(5) + t158;
t442 = t155 / 0.2e1;
t538 = mrSges(6,1) + mrSges(7,1);
t537 = mrSges(6,2) - mrSges(7,3);
t220 = t511 * qJD(1);
t211 = qJD(4) - t220;
t203 = qJD(5) + t211;
t436 = t203 / 0.2e1;
t445 = t290 / 0.2e1;
t448 = t125 / 0.2e1;
t449 = -t125 / 0.2e1;
t488 = Ifges(7,4) + Ifges(6,5);
t490 = Ifges(6,1) + Ifges(7,1);
t414 = pkin(7) + qJ(2);
t241 = t414 * t273;
t229 = qJD(1) * t241;
t242 = t414 * t274;
t230 = qJD(1) * t242;
t167 = -t282 * t229 - t278 * t230;
t159 = -qJD(3) * pkin(3) - t167;
t119 = -pkin(4) * t301 + t159;
t260 = pkin(2) * t274 + pkin(1);
t239 = -qJD(1) * t260 + qJD(2);
t135 = -pkin(3) * t220 - pkin(8) * t221 + t239;
t168 = -t278 * t229 + t282 * t230;
t160 = qJD(3) * pkin(8) + t168;
t91 = t277 * t135 + t281 * t160;
t80 = pkin(9) * t301 + t91;
t394 = t276 * t80;
t90 = t281 * t135 - t160 * t277;
t79 = -pkin(9) * t179 + t90;
t70 = pkin(4) * t211 + t79;
t27 = t280 * t70 - t394;
t534 = qJD(6) - t27;
t23 = -pkin(5) * t203 + t534;
t123 = Ifges(6,4) * t125;
t400 = Ifges(7,5) * t125;
t484 = t488 * t203 + t290 * t490 - t123 + t400;
t524 = t484 / 0.2e1;
t56 = t125 * pkin(5) - qJ(6) * t290 + t119;
t509 = t119 * mrSges(6,2) + mrSges(7,2) * t23 - mrSges(6,3) * t27 - t56 * mrSges(7,3) + t524;
t547 = Ifges(6,4) * t449 + Ifges(7,5) * t448 + t436 * t488 + t445 * t490 + t509;
t393 = t280 * t80;
t28 = t276 * t70 + t393;
t24 = qJ(6) * t203 + t28;
t122 = Ifges(7,5) * t290;
t64 = t203 * Ifges(7,6) + t125 * Ifges(7,3) + t122;
t401 = Ifges(6,4) * t290;
t67 = -t125 * Ifges(6,2) + t203 * Ifges(6,6) + t401;
t546 = -mrSges(7,2) * t24 - mrSges(6,3) * t28 + mrSges(6,1) * t119 + t56 * mrSges(7,1) + t64 / 0.2e1 - t67 / 0.2e1;
t487 = -Ifges(6,6) + Ifges(7,6);
t520 = Ifges(6,3) + Ifges(7,2);
t161 = pkin(3) * t221 - pkin(8) * t220;
t108 = t277 * t161 + t281 * t167;
t340 = qJD(4) * t284;
t392 = t220 * t277;
t545 = pkin(9) * t392 + t277 * t340 - t108;
t107 = t281 * t161 - t167 * t277;
t391 = t220 * t281;
t544 = -pkin(4) * t221 + pkin(9) * t391 + t281 * t340 - t107;
t360 = t273 ^ 2 + t274 ^ 2;
t352 = qJD(1) * qJD(2);
t251 = qJ(2) * qJDD(1) + t352;
t357 = qJD(4) * t277;
t543 = t357 - t392;
t542 = -Ifges(6,2) * t449 + Ifges(7,3) * t448 + t436 * t487 + t445 * t489 + t546;
t235 = -qJDD(1) * t260 + qJDD(2);
t98 = -pkin(3) * t165 - pkin(8) * t164 + t235;
t329 = pkin(7) * qJDD(1) + t251;
t205 = t329 * t273;
t206 = t329 * t274;
t358 = qJD(3) * t282;
t359 = qJD(3) * t278;
t105 = -t278 * t205 + t282 * t206 - t229 * t358 - t230 * t359;
t99 = qJDD(3) * pkin(8) + t105;
t26 = -t91 * qJD(4) - t277 * t99 + t281 * t98;
t20 = pkin(4) * t158 - pkin(9) * t117 + t26;
t356 = qJD(4) * t281;
t25 = t135 * t356 - t160 * t357 + t277 * t98 + t281 * t99;
t22 = pkin(9) * t118 + t25;
t6 = -qJD(5) * t28 + t20 * t280 - t22 * t276;
t3 = -pkin(5) * t155 + qJDD(6) - t6;
t459 = -t49 / 0.2e1;
t106 = -t205 * t282 - t278 * t206 + t229 * t359 - t230 * t358;
t100 = -qJDD(3) * pkin(3) - t106;
t59 = -pkin(4) * t118 + t100;
t9 = pkin(5) * t49 - qJ(6) * t48 - qJD(6) * t290 + t59;
t540 = mrSges(6,2) * t59 + mrSges(7,2) * t3 - mrSges(6,3) * t6 - mrSges(7,3) * t9 + Ifges(6,4) * t459 + 0.2e1 * t442 * t488 + t458 * t549 + 0.2e1 * t460 * t490;
t231 = t276 * t277 - t280 * t281;
t140 = t231 * t220;
t467 = qJD(4) + qJD(5);
t170 = t467 * t231;
t474 = -t170 + t140;
t376 = t276 * t281;
t232 = t277 * t280 + t376;
t139 = t232 * t220;
t171 = t467 * t232;
t473 = t171 - t139;
t409 = mrSges(4,3) * t221;
t472 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t301 - t179 * mrSges(5,2) - t409;
t532 = -m(5) * t159 + t472;
t272 = qJ(4) + qJ(5);
t267 = sin(t272);
t268 = cos(t272);
t315 = -mrSges(5,1) * t281 + mrSges(5,2) * t277;
t530 = m(5) * pkin(3) - t267 * t537 + t268 * t538 - t315;
t354 = qJD(5) * t280;
t355 = qJD(5) * t276;
t5 = t276 * t20 + t280 * t22 + t70 * t354 - t355 * t80;
t2 = qJ(6) * t155 + qJD(6) * t203 + t5;
t528 = mrSges(6,1) * t59 + mrSges(7,1) * t9 - mrSges(7,2) * t2 - mrSges(6,3) * t5 - t48 * Ifges(6,4) / 0.2e1 - t155 * Ifges(6,6) / 0.2e1 + 0.2e1 * Ifges(7,3) * t458 + (-t459 + t458) * Ifges(6,2) + t549 * t460 + (t487 + Ifges(7,6)) * t442;
t525 = -m(5) - m(4);
t522 = t90 * mrSges(5,1);
t521 = t91 * mrSges(5,2);
t518 = Ifges(5,3) * t211;
t516 = t301 * Ifges(5,6);
t244 = t284 * t277;
t245 = t284 * t281;
t302 = t280 * t244 + t245 * t276;
t481 = qJD(5) * t302 + t276 * t544 + t280 * t545;
t178 = t244 * t276 - t245 * t280;
t480 = -qJD(5) * t178 - t276 * t545 + t280 * t544;
t515 = Ifges(4,5) * qJD(3);
t514 = Ifges(4,6) * qJD(3);
t406 = mrSges(6,3) * t290;
t103 = mrSges(6,1) * t203 - t406;
t104 = -mrSges(7,1) * t203 + mrSges(7,2) * t290;
t476 = t104 - t103;
t475 = pkin(4) * t543 - t168;
t336 = t228 * t356;
t390 = t222 * t277;
t299 = t336 + t390;
t279 = sin(qJ(1));
t371 = t279 * t281;
t283 = cos(qJ(1));
t373 = t277 * t283;
t214 = -t266 * t373 + t371;
t512 = t155 * t520 + t48 * t488 + t487 * t49;
t508 = Ifges(5,5) * t179 + t125 * t487 + t203 * t520 + t290 * t488 + t516 + t518;
t334 = m(3) * qJ(2) + mrSges(3,3);
t507 = -t334 - mrSges(4,3) + mrSges(2,2);
t437 = -t203 / 0.2e1;
t446 = -t290 / 0.2e1;
t505 = -Ifges(6,2) * t448 + Ifges(7,3) * t449 + t437 * t487 + t446 * t489 - t546;
t504 = -m(7) * qJ(6) - mrSges(7,3);
t75 = pkin(5) * t290 + qJ(6) * t125;
t503 = m(7) * pkin(5) + t538;
t317 = t266 * mrSges(4,1) - t265 * mrSges(4,2);
t318 = -mrSges(3,1) * t274 + mrSges(3,2) * t273;
t502 = m(3) * pkin(1) + t265 * mrSges(5,3) + mrSges(2,1) + t317 - t318;
t501 = -mrSges(6,2) - t504;
t500 = t26 * mrSges(5,1) - t25 * mrSges(5,2);
t497 = t6 * mrSges(6,1) - t3 * mrSges(7,1) - t5 * mrSges(6,2) + t2 * mrSges(7,3);
t451 = t117 / 0.2e1;
t450 = t118 / 0.2e1;
t441 = t158 / 0.2e1;
t492 = -t301 / 0.2e1;
t29 = -mrSges(7,2) * t49 + mrSges(7,3) * t155;
t32 = -mrSges(6,2) * t155 - mrSges(6,3) * t49;
t486 = t29 + t32;
t30 = mrSges(6,1) * t155 - mrSges(6,3) * t48;
t31 = -t155 * mrSges(7,1) + t48 * mrSges(7,2);
t485 = t31 - t30;
t483 = pkin(5) * t473 - qJ(6) * t474 - qJD(6) * t232 + t475;
t482 = -qJ(6) * t221 + t481;
t479 = pkin(5) * t221 - t480;
t150 = t232 * t228;
t101 = -mrSges(7,2) * t125 + mrSges(7,3) * t203;
t407 = mrSges(6,3) * t125;
t102 = -mrSges(6,2) * t203 - t407;
t477 = t102 + t101;
t163 = -pkin(3) * t511 - pkin(8) * t228 - t260;
t175 = -t241 * t278 + t242 * t282;
t169 = t281 * t175;
t116 = t277 * t163 + t169;
t471 = -t282 * t241 - t242 * t278;
t367 = t283 * t267;
t372 = t279 * t268;
t197 = t266 * t367 - t372;
t381 = t268 * t283;
t382 = t267 * t279;
t198 = t266 * t381 + t382;
t470 = t197 * t538 + t198 * t537;
t195 = t266 * t382 + t381;
t196 = t266 * t372 - t367;
t469 = t195 * t538 + t196 * t537;
t468 = t25 * t281 - t26 * t277;
t136 = qJD(2) * t511 + qJD(3) * t471;
t162 = pkin(3) * t223 - pkin(8) * t222;
t326 = -t136 * t277 + t281 * t162;
t389 = t222 * t281;
t38 = -pkin(9) * t389 + pkin(4) * t223 + (-t169 + (pkin(9) * t228 - t163) * t277) * qJD(4) + t326;
t115 = t281 * t163 - t175 * t277;
t387 = t228 * t281;
t86 = -pkin(4) * t511 - pkin(9) * t387 + t115;
t388 = t228 * t277;
t92 = -pkin(9) * t388 + t116;
t413 = t276 * t86 + t280 * t92;
t57 = t281 * t136 + t277 * t162 + t163 * t356 - t175 * t357;
t51 = -pkin(9) * t299 + t57;
t11 = -qJD(5) * t413 - t276 * t51 + t280 * t38;
t457 = Ifges(5,1) * t451 + Ifges(5,4) * t450 + Ifges(5,5) * t441;
t439 = -t179 / 0.2e1;
t438 = t179 / 0.2e1;
t435 = -t211 / 0.2e1;
t433 = t220 / 0.2e1;
t431 = t221 / 0.2e1;
t422 = pkin(4) * t179;
t421 = pkin(4) * t276;
t420 = pkin(4) * t277;
t419 = pkin(4) * t280;
t416 = g(3) * t265;
t412 = mrSges(6,2) * t268;
t410 = mrSges(4,3) * t220;
t408 = mrSges(5,3) * t179;
t405 = Ifges(4,4) * t221;
t404 = Ifges(5,4) * t179;
t403 = Ifges(5,4) * t277;
t402 = Ifges(5,4) * t281;
t377 = t414 * t283;
t111 = Ifges(5,2) * t301 + Ifges(5,6) * t211 + t404;
t375 = t277 * t111;
t374 = t277 * t279;
t176 = Ifges(5,4) * t301;
t112 = t179 * Ifges(5,1) + t211 * Ifges(5,5) + t176;
t370 = t281 * t112;
t369 = t281 * t283;
t345 = pkin(4) * t355;
t344 = pkin(4) * t354;
t343 = m(5) * pkin(8) + mrSges(5,3);
t341 = Ifges(5,5) * t117 + Ifges(5,6) * t118 + Ifges(5,3) * t158;
t337 = t228 * t357;
t335 = t370 / 0.2e1;
t333 = -t357 / 0.2e1;
t332 = t504 * t265 * t268;
t331 = -t165 * mrSges(4,1) + t164 * mrSges(4,2);
t328 = -t195 * pkin(5) + qJ(6) * t196;
t327 = -t197 * pkin(5) + qJ(6) * t198;
t324 = t283 * t260 + t279 * t414;
t321 = pkin(3) * t266 + pkin(8) * t265;
t138 = pkin(4) * t388 - t471;
t319 = -mrSges(3,1) * t349 + mrSges(3,2) * t350;
t314 = mrSges(5,1) * t277 + t281 * mrSges(5,2);
t311 = Ifges(5,1) * t281 - t403;
t310 = -Ifges(5,2) * t277 + t402;
t309 = Ifges(5,5) * t281 - Ifges(5,6) * t277;
t308 = pkin(5) * t268 + qJ(6) * t267;
t52 = -t276 * t92 + t280 * t86;
t296 = t301 * mrSges(5,3);
t133 = -t211 * mrSges(5,2) + t296;
t134 = mrSges(5,1) * t211 - t408;
t303 = t133 * t281 - t134 * t277;
t300 = t214 * pkin(4);
t212 = t266 * t374 + t369;
t10 = t276 * t38 + t280 * t51 + t86 * t354 - t355 * t92;
t298 = t337 - t389;
t297 = t159 * t314;
t292 = t212 * pkin(4);
t289 = t497 + t512;
t137 = qJD(2) * t228 + qJD(3) * t175;
t96 = pkin(4) * t299 + t137;
t264 = -qJDD(1) * pkin(1) + qJDD(2);
t262 = -pkin(5) - t419;
t259 = qJ(6) + t421;
t252 = qJD(6) + t344;
t215 = t266 * t369 + t374;
t213 = -t266 * t371 + t373;
t207 = Ifges(4,4) * t220;
t190 = -qJD(3) * mrSges(4,2) + t410;
t166 = pkin(5) * t231 - qJ(6) * t232 - t263;
t151 = t231 * t228;
t143 = t221 * Ifges(4,1) + t207 + t515;
t142 = t220 * Ifges(4,2) + t405 + t514;
t83 = -mrSges(5,2) * t158 + mrSges(5,3) * t118;
t82 = mrSges(5,1) * t158 - mrSges(5,3) * t117;
t77 = mrSges(6,1) * t125 + mrSges(6,2) * t290;
t76 = mrSges(7,1) * t125 - mrSges(7,3) * t290;
t73 = pkin(5) * t150 + qJ(6) * t151 + t138;
t72 = t222 * t376 - t276 * t337 - t355 * t388 + (t387 * t467 + t390) * t280;
t71 = -t150 * t467 - t231 * t222;
t62 = t422 + t75;
t61 = -mrSges(5,1) * t118 + mrSges(5,2) * t117;
t58 = -qJD(4) * t116 + t326;
t54 = t117 * Ifges(5,4) + t118 * Ifges(5,2) + t158 * Ifges(5,6);
t42 = pkin(5) * t511 - t52;
t41 = -qJ(6) * t511 + t413;
t34 = t280 * t79 - t394;
t33 = t276 * t79 + t393;
t18 = pkin(5) * t72 - qJ(6) * t71 + qJD(6) * t151 + t96;
t17 = mrSges(6,1) * t49 + mrSges(6,2) * t48;
t16 = mrSges(7,1) * t49 - mrSges(7,3) * t48;
t8 = -pkin(5) * t223 - t11;
t7 = qJ(6) * t223 - qJD(6) * t511 + t10;
t1 = [t542 * t72 - (-m(4) * t106 + m(5) * t100 - qJDD(3) * mrSges(4,1) + mrSges(4,3) * t164 + t61) * t471 + (-m(4) * t167 - t532) * t137 + 0.2e1 * t360 * t251 * mrSges(3,3) - (-t239 * mrSges(4,2) - t515 / 0.2e1 - Ifges(4,1) * t431 - Ifges(4,4) * t433 + t375 / 0.2e1 - t143 / 0.2e1 - t335 + t167 * mrSges(4,3)) * t222 + (-Ifges(5,1) * t298 - Ifges(5,4) * t299) * t438 + t528 * t150 - (t235 * mrSges(4,1) - t105 * mrSges(4,3) - Ifges(4,4) * t164 + Ifges(5,5) * t451 - Ifges(4,2) * t165 - Ifges(4,6) * qJDD(3) + Ifges(5,6) * t450 + Ifges(6,6) * t459 + Ifges(7,6) * t458 + Ifges(5,3) * t441 + t520 * t442 + t488 * t460 + t497 + t500) * t511 - (t341 + t512) * t511 / 0.2e1 + m(6) * (t10 * t28 + t11 * t27 + t119 * t96 + t138 * t59 + t413 * t5 + t52 * t6) + t413 * t32 + t547 * t71 - t54 * t388 / 0.2e1 - t540 * t151 + m(5) * (t115 * t26 + t116 * t25 + t57 * t91 + t58 * t90) + m(4) * (t105 * t175 + t136 * t168 - t235 * t260) + (t508 / 0.2e1 + t239 * mrSges(4,1) - t514 / 0.2e1 - t142 / 0.2e1 + t27 * mrSges(6,1) + t24 * mrSges(7,3) - t28 * mrSges(6,2) - t23 * mrSges(7,1) + t516 / 0.2e1 + t522 - t521 + t518 / 0.2e1 - mrSges(4,3) * t168 - Ifges(4,4) * t431 - Ifges(4,2) * t433 + Ifges(5,5) * t438 + Ifges(7,6) * t448 + Ifges(6,6) * t449 + t488 * t445 + t520 * t436) * t223 + (-t25 * t388 - t26 * t387 + t298 * t90 - t299 * t91) * mrSges(5,3) + t211 * (-Ifges(5,5) * t298 - Ifges(5,6) * t299) / 0.2e1 + m(7) * (t18 * t56 + t2 * t41 + t23 * t8 + t24 * t7 + t3 * t42 + t73 * t9) + t301 * (-Ifges(5,4) * t298 - Ifges(5,2) * t299) / 0.2e1 + (-t215 * mrSges(5,1) - t214 * mrSges(5,2) + t525 * t324 - t539 * (pkin(4) * t374 + t324) - t503 * t198 - t501 * t197 + t507 * t279 + (-m(5) * t321 - t502 + t550) * t283) * g(2) + (-m(5) * t377 - t213 * mrSges(5,1) - t212 * mrSges(5,2) - t539 * (pkin(4) * t373 + t279 * t383 + t377) + t503 * t196 + t501 * t195 + (-m(4) * t414 + t507) * t283 + (-m(5) * (-t260 - t321) + m(4) * t260 - t539 * (-t260 - t237) + t502 + t533) * t279) * g(1) + Ifges(2,3) * qJDD(1) - t111 * t336 / 0.2e1 - t260 * t331 - pkin(1) * t319 + t264 * t318 + t136 * t190 + t175 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t165) + (t235 * mrSges(4,2) - t106 * mrSges(4,3) + Ifges(4,1) * t164 + Ifges(4,4) * t165 + Ifges(4,5) * qJDD(3) + t100 * t314 + t112 * t333 + t309 * t441 + t310 * t450 + t311 * t451) * t228 + t138 * t17 + t57 * t133 + t58 * t134 + t159 * (mrSges(5,1) * t299 - mrSges(5,2) * t298) + t115 * t82 + t116 * t83 + t7 * t101 + t10 * t102 + t11 * t103 + t8 * t104 + t96 * t77 + t18 * t76 + t73 * t16 + t387 * t457 + m(3) * (-pkin(1) * t264 + (t251 + t352) * qJ(2) * t360) + t41 * t29 + t42 * t31 + t52 * t30 + (Ifges(3,4) * t273 + Ifges(3,2) * t274) * t349 + (Ifges(3,1) * t273 + Ifges(3,4) * t274) * t350; t319 + t486 * t232 + t485 * t231 + (-t76 - t77 + t472) * t221 + m(3) * t264 + t331 + t281 * t82 + t277 * t83 + (-t190 - t303) * t220 + t303 * qJD(4) + t477 * t474 + t476 * t473 + (-g(1) * t279 + g(2) * t283) * (m(3) + t539 - t525) - t334 * t360 * qJD(1) ^ 2 + (t2 * t232 - t221 * t56 + t23 * t473 + t231 * t3 + t24 * t474) * m(7) + (-t119 * t221 - t231 * t6 + t232 * t5 - t27 * t473 + t28 * t474) * m(6) + (-t159 * t221 + t25 * t277 + t26 * t281 + t211 * (-t277 * t90 + t281 * t91)) * m(5) + (t167 * t221 - t168 * t220 + t235) * m(4); t542 * t171 + (-t543 * t91 + (-t356 + t391) * t90 + t468) * mrSges(5,3) + (-t488 * t140 + t221 * t520) * t437 + (t409 + t532) * t168 + t475 * t77 + (t119 * t475 + t178 * t5 - t263 * t59 + t27 * t480 + t28 * t481 + t302 * t6) * m(6) + (t166 * t9 + t178 * t2 + t23 * t479 + t24 * t482 - t3 * t302 + t483 * t56) * m(7) - t485 * t302 + t528 * t231 + (t335 + t297) * qJD(4) + t505 * t139 - (Ifges(4,1) * t220 - t405 + t508) * t221 / 0.2e1 + (-Ifges(6,4) * t140 + Ifges(6,6) * t221) * t448 + (-Ifges(7,5) * t140 + Ifges(7,6) * t221) * t449 - t23 * (-mrSges(7,1) * t221 - mrSges(7,2) * t140) - t27 * (mrSges(6,1) * t221 + mrSges(6,3) * t140) + (t179 * t311 + t211 * t309) * qJD(4) / 0.2e1 - (-Ifges(4,2) * t221 + t143 + t207 + t370) * t220 / 0.2e1 + (t119 * t140 + t221 * t28) * mrSges(6,2) + (-t140 * t56 - t221 * t24) * mrSges(7,3) + (-pkin(3) * t100 - t107 * t90 - t108 * t91) * m(5) + (m(5) * ((-t277 * t91 - t281 * t90) * qJD(4) + t468) + t281 * t83 - t277 * t82 - t134 * t356 - t133 * t357) * pkin(8) - t547 * t170 + (Ifges(5,2) * t281 + t403) * t450 + (Ifges(5,1) * t277 + t402) * t451 + (Ifges(5,6) * t221 + t220 * t310) * t492 + t540 * t232 + (-t265 * t343 - t317 + (-m(7) * t308 - t530) * t266 + t550) * g(3) - t221 * t522 + (g(1) * t283 + g(2) * t279) * ((t284 * t539 + mrSges(4,2) - t343 - t548) * t266 + (m(6) * t263 - m(7) * (-t263 - t308) + mrSges(4,1) + t530) * t265) + t281 * t54 / 0.2e1 + (-t190 + t410) * t167 + Ifges(4,3) * qJDD(3) - t263 * t17 - t239 * (mrSges(4,1) * t221 + mrSges(4,2) * t220) - qJD(3) * (Ifges(4,5) * t220 - Ifges(4,6) * t221) / 0.2e1 + t100 * t315 + t310 * t293 / 0.2e1 + Ifges(4,5) * t164 + Ifges(4,6) * t165 + t166 * t16 - t108 * t133 - t107 * t134 - t220 * t297 - t105 * mrSges(4,2) + t106 * mrSges(4,1) + t277 * t457 + t479 * t104 + t480 * t103 + t481 * t102 + t482 * t101 + t111 * t333 + t483 * t76 + t486 * t178 + (-t140 * t490 + t221 * t488) * t446 - pkin(3) * t61 + t221 * t521 + t140 * t524 + t142 * t431 + t375 * t433 + (Ifges(5,3) * t221 + t220 * t309) * t435 + (Ifges(5,5) * t221 + t220 * t311) * t439 + (Ifges(5,5) * t277 + Ifges(5,6) * t281) * t441; (-t133 + t296) * t90 + (-Ifges(6,4) * t448 - Ifges(7,5) * t449 - t488 * t437 - t490 * t446 + t509) * t125 + t505 * t290 - g(3) * ((m(7) * (-pkin(5) * t267 - t420) - t267 * mrSges(7,1)) * t265 - t332) - t77 * t422 + (-m(7) * (-t292 + t328) + m(6) * t292 + mrSges(5,1) * t212 - mrSges(5,2) * t213 + t469) * g(2) + (-m(6) * t300 - m(7) * (t300 + t327) - mrSges(5,1) * t214 + mrSges(5,2) * t215 + t470) * g(1) + ((t276 * t5 + t280 * t6 + (-t27 * t276 + t28 * t280) * qJD(5)) * pkin(4) + t420 * t416 - t119 * t422 + t27 * t33 - t28 * t34) * m(6) - m(7) * (t23 * t33 + t24 * t34 + t56 * t62) + t341 + (mrSges(6,1) * t267 + t314 + t412) * t416 + (-t34 + t344) * t102 + (-Ifges(5,2) * t179 + t112 + t176) * t492 + t289 + (t252 - t34) * t101 + m(7) * (t2 * t259 + t23 * t345 + t24 * t252 + t262 * t3) + t259 * t29 + t262 * t31 - t159 * (t179 * mrSges(5,1) + mrSges(5,2) * t301) - t62 * t76 - t476 * (-t345 + t33) + t500 + (t134 + t408) * t91 + t30 * t419 + t32 * t421 + (Ifges(5,5) * t301 - Ifges(5,6) * t179) * t435 + t111 * t438 + (Ifges(5,1) * t301 - t404) * t439; (Ifges(7,3) * t290 - t400) * t449 + ((t267 * t503 + t412) * t265 + t332) * g(3) + (t406 - t476) * t28 + (-t407 - t477) * t27 + t470 * g(1) + t469 * g(2) + t289 - t56 * (mrSges(7,1) * t290 + mrSges(7,3) * t125) - t119 * (mrSges(6,1) * t290 - mrSges(6,2) * t125) + qJD(6) * t101 - t75 * t76 + (t125 * t23 + t24 * t290) * mrSges(7,2) + qJ(6) * t29 - pkin(5) * t31 + t67 * t445 + (-t125 * t488 + t290 * t487) * t437 + (-Ifges(6,2) * t290 - t123 + t484) * t448 + (-t125 * t490 + t122 - t401 + t64) * t446 + (-pkin(5) * t3 - t327 * g(1) - t328 * g(2) + qJ(6) * t2 - t23 * t28 + t24 * t534 - t56 * t75) * m(7); -t203 * t101 + t290 * t76 + (-g(1) * t197 - g(2) * t195 - t24 * t203 - t267 * t416 + t290 * t56 + t3) * m(7) + t31;];
tau  = t1;
