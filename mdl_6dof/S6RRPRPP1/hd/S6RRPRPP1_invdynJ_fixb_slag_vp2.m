% Calculate vector of inverse dynamics joint torques for
% S6RRPRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
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
% Datum: 2019-03-09 09:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:45:09
% EndTime: 2019-03-09 09:45:54
% DurationCPUTime: 27.82s
% Computational Cost: add. (12126->745), mult. (27694->932), div. (0->0), fcn. (20287->14), ass. (0->349)
t530 = mrSges(6,1) + mrSges(7,1);
t529 = -mrSges(6,2) + mrSges(7,3);
t275 = sin(qJ(2));
t278 = cos(qJ(2));
t386 = sin(pkin(9));
t388 = cos(pkin(9));
t231 = t388 * t275 + t386 * t278;
t214 = t231 * qJD(1);
t274 = sin(qJ(4));
t277 = cos(qJ(4));
t184 = qJD(2) * t277 - t214 * t274;
t185 = qJD(2) * t274 + t214 * t277;
t271 = sin(pkin(10));
t387 = cos(pkin(10));
t119 = -t387 * t184 + t185 * t271;
t442 = -t119 / 0.2e1;
t229 = t275 * t386 - t278 * t388;
t213 = t229 * qJD(1);
t204 = t213 + qJD(4);
t431 = -t204 / 0.2e1;
t430 = t204 / 0.2e1;
t290 = t271 * t184 + t185 * t387;
t439 = -t290 / 0.2e1;
t438 = t290 / 0.2e1;
t336 = t386 * pkin(2);
t257 = t336 + pkin(8);
t363 = qJ(5) + t257;
t317 = qJD(4) * t363;
t384 = t213 * t274;
t360 = qJD(1) * t275;
t347 = pkin(2) * t360;
t142 = pkin(3) * t214 + pkin(8) * t213 + t347;
t273 = -qJ(3) - pkin(7);
t245 = t273 * t278;
t234 = qJD(1) * t245;
t218 = t386 * t234;
t243 = t273 * t275;
t233 = qJD(1) * t243;
t174 = t233 * t388 + t218;
t97 = t274 * t142 + t277 * t174;
t528 = -qJ(5) * t384 + qJD(5) * t277 - t274 * t317 - t97;
t383 = t213 * t277;
t96 = t277 * t142 - t174 * t274;
t527 = -pkin(4) * t214 - qJ(5) * t383 - qJD(5) * t274 - t277 * t317 - t96;
t354 = qJD(4) * t274;
t526 = t354 + t384;
t223 = qJD(2) * pkin(2) + t233;
t165 = t223 * t388 + t218;
t153 = -qJD(2) * pkin(3) - t165;
t117 = -t184 * pkin(4) + qJD(5) + t153;
t268 = t278 * pkin(2);
t261 = t268 + pkin(1);
t239 = -qJD(1) * t261 + qJD(3);
t131 = pkin(3) * t213 - pkin(8) * t214 + t239;
t328 = t388 * t234;
t166 = t386 * t223 - t328;
t154 = qJD(2) * pkin(8) + t166;
t86 = t277 * t131 - t154 * t274;
t78 = -qJ(5) * t185 + t86;
t70 = pkin(4) * t204 + t78;
t87 = t131 * t274 + t154 * t277;
t79 = qJ(5) * t184 + t87;
t73 = t387 * t79;
t26 = t271 * t70 + t73;
t22 = qJ(6) * t204 + t26;
t41 = t119 * pkin(5) - qJ(6) * t290 + t117;
t525 = Ifges(6,4) * t438 + Ifges(7,5) * t439 + Ifges(6,6) * t430 + Ifges(7,6) * t431 + (Ifges(6,2) + Ifges(7,3)) * t442 - t117 * mrSges(6,1) - t41 * mrSges(7,1) + mrSges(7,2) * t22 + mrSges(6,3) * t26;
t352 = qJD(1) * qJD(2);
t332 = t275 * t352;
t351 = qJDD(1) * t278;
t235 = -t332 + t351;
t236 = qJDD(1) * t275 + t278 * t352;
t176 = t235 * t386 + t236 * t388;
t115 = qJD(4) * t184 + qJDD(2) * t274 + t176 * t277;
t116 = -qJD(4) * t185 + qJDD(2) * t277 - t176 * t274;
t60 = t115 * t271 - t116 * t387;
t451 = -t60 / 0.2e1;
t61 = t115 * t387 + t271 * t116;
t449 = t61 / 0.2e1;
t524 = m(5) + m(4);
t445 = m(6) + m(7);
t175 = t235 * t388 - t236 * t386;
t172 = qJDD(4) - t175;
t435 = t172 / 0.2e1;
t523 = mrSges(4,2) - mrSges(5,3);
t522 = -mrSges(6,3) - mrSges(7,2);
t502 = Ifges(6,1) + Ifges(7,1);
t501 = Ifges(6,4) - Ifges(7,5);
t500 = Ifges(7,4) + Ifges(6,5);
t499 = Ifges(6,6) - Ifges(7,6);
t515 = pkin(4) * t445;
t521 = mrSges(5,1) + t515;
t230 = t271 * t277 + t387 * t274;
t132 = t230 * t213;
t212 = t230 * qJD(4);
t480 = t212 + t132;
t325 = t387 * t277;
t289 = -t271 * t274 + t325;
t133 = t289 * t213;
t216 = t289 * qJD(4);
t479 = t216 + t133;
t390 = t271 * t79;
t25 = t387 * t70 - t390;
t21 = -t204 * pkin(5) + qJD(6) - t25;
t520 = -t117 * mrSges(6,2) - mrSges(7,2) * t21 + mrSges(6,3) * t25 + t41 * mrSges(7,3);
t441 = t119 / 0.2e1;
t518 = Ifges(6,2) * t442 - Ifges(7,3) * t441 + t499 * t430 + t501 * t438 + t525;
t475 = Ifges(5,3) + Ifges(6,3) + Ifges(7,2);
t269 = qJ(4) + pkin(10);
t264 = sin(t269);
t266 = cos(t269);
t310 = -mrSges(5,1) * t277 + mrSges(5,2) * t274;
t517 = m(5) * pkin(3) + t529 * t264 + t266 * t530 - t310;
t495 = -t119 * t501 + t204 * t500 + t290 * t502;
t516 = t495 / 0.2e1;
t244 = -mrSges(3,1) * t278 + mrSges(3,2) * t275;
t514 = m(3) * pkin(1) + mrSges(2,1) - t244;
t492 = -t271 * t528 + t387 * t527;
t489 = t271 * t527 + t387 * t528;
t173 = t233 * t386 - t328;
t482 = pkin(4) * t526 - t173;
t217 = t229 * qJD(2);
t353 = qJD(4) * t277;
t334 = t231 * t353;
t294 = -t217 * t274 + t334;
t270 = qJ(2) + pkin(9);
t267 = cos(t270);
t276 = sin(qJ(1));
t367 = t276 * t277;
t279 = cos(qJ(1));
t371 = t274 * t279;
t207 = -t267 * t371 + t367;
t513 = g(1) * t279 + g(2) * t276;
t227 = t236 * pkin(7);
t357 = qJD(3) * t275;
t161 = qJDD(2) * pkin(2) - qJ(3) * t236 - qJD(1) * t357 - t227;
t262 = pkin(7) * t351;
t358 = qJD(2) * t275;
t344 = pkin(7) * t358;
t356 = qJD(3) * t278;
t170 = qJ(3) * t235 + t262 + (-t344 + t356) * qJD(1);
t106 = t386 * t161 + t388 * t170;
t100 = qJDD(2) * pkin(8) + t106;
t385 = qJDD(1) * pkin(1);
t201 = -pkin(2) * t235 + qJDD(3) - t385;
t94 = -pkin(3) * t175 - pkin(8) * t176 + t201;
t24 = -qJD(4) * t87 - t100 * t274 + t277 * t94;
t11 = pkin(4) * t172 - qJ(5) * t115 - qJD(5) * t185 + t24;
t23 = t277 * t100 + t131 * t353 - t154 * t354 + t274 * t94;
t13 = qJ(5) * t116 + qJD(5) * t184 + t23;
t3 = t11 * t387 - t271 * t13;
t2 = -t172 * pkin(5) + qJDD(6) - t3;
t450 = t60 / 0.2e1;
t105 = t161 * t388 - t386 * t170;
t99 = -qJDD(2) * pkin(3) - t105;
t50 = -t116 * pkin(4) + qJDD(5) + t99;
t6 = t60 * pkin(5) - t61 * qJ(6) - qJD(6) * t290 + t50;
t509 = mrSges(6,2) * t50 + mrSges(7,2) * t2 - mrSges(6,3) * t3 - mrSges(7,3) * t6 + Ifges(7,5) * t450 + 0.2e1 * t435 * t500 + 0.2e1 * t449 * t502 + (Ifges(6,4) + t501) * t451;
t508 = -Ifges(6,2) * t441 + Ifges(7,3) * t442 - t439 * t501 + t525;
t507 = Ifges(6,4) * t442 + Ifges(7,5) * t441 + t430 * t500 + t438 * t502 + t516 - t520;
t506 = Ifges(6,4) * t441 + Ifges(7,5) * t442 + t500 * t431 + t439 * t502 + t520;
t444 = t115 / 0.2e1;
t443 = t116 / 0.2e1;
t503 = t235 / 0.2e1;
t407 = qJD(2) / 0.2e1;
t42 = -mrSges(7,2) * t60 + mrSges(7,3) * t172;
t43 = -mrSges(6,2) * t172 - mrSges(6,3) * t60;
t497 = t42 + t43;
t44 = mrSges(6,1) * t172 - mrSges(6,3) * t61;
t45 = -t172 * mrSges(7,1) + t61 * mrSges(7,2);
t496 = t45 - t44;
t494 = t184 * Ifges(5,6);
t389 = qJDD(2) / 0.2e1;
t493 = pkin(5) * t480 - qJ(6) * t479 - qJD(6) * t230 + t482;
t491 = t214 * pkin(5) - t492;
t490 = -qJ(6) * t214 + t489;
t488 = Ifges(4,5) * qJD(2);
t487 = Ifges(4,6) * qJD(2);
t265 = sin(t270);
t486 = t265 * t513;
t30 = t387 * t78 - t390;
t485 = qJD(6) - t30;
t101 = -mrSges(7,2) * t119 + mrSges(7,3) * t204;
t102 = -mrSges(6,2) * t204 - mrSges(6,3) * t119;
t484 = t102 + t101;
t103 = mrSges(6,1) * t204 - mrSges(6,3) * t290;
t104 = -mrSges(7,1) * t204 + mrSges(7,2) * t290;
t483 = t104 - t103;
t164 = pkin(3) * t229 - pkin(8) * t231 - t261;
t181 = t243 * t386 - t245 * t388;
t177 = t277 * t181;
t112 = t274 * t164 + t177;
t481 = qJD(2) * mrSges(4,1) + mrSges(5,1) * t184 - mrSges(5,2) * t185 - mrSges(4,3) * t214;
t478 = -t267 * mrSges(4,1) + t265 * t523;
t226 = -pkin(7) * t332 + t262;
t477 = t226 * t278 + t227 * t275;
t476 = t23 * t277 - t24 * t274;
t474 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t473 = t185 * Ifges(5,5) - t119 * t499 + t204 * t475 + t290 * t500 + t494;
t472 = 0.2e1 * t389;
t471 = Ifges(5,5) * t115 + Ifges(5,6) * t116 + t172 * t475 - t499 * t60 + t500 * t61;
t470 = -t239 * mrSges(4,2) + t165 * mrSges(4,3);
t469 = m(7) * pkin(5) + t530;
t468 = -t265 * t522 - t478;
t467 = m(7) * qJ(6) + t529;
t4 = t271 * t11 + t387 * t13;
t1 = qJ(6) * t172 + qJD(6) * t204 + t4;
t462 = t24 * mrSges(5,1) + t3 * mrSges(6,1) - t2 * mrSges(7,1) - t23 * mrSges(5,2) - t4 * mrSges(6,2) + t1 * mrSges(7,3);
t459 = t239 * mrSges(4,1) + t86 * mrSges(5,1) + t25 * mrSges(6,1) - t21 * mrSges(7,1) - t87 * mrSges(5,2) - t26 * mrSges(6,2) - t166 * mrSges(4,3) + t22 * mrSges(7,3);
t456 = mrSges(6,1) * t50 + mrSges(7,1) * t6 - mrSges(7,2) * t1 - mrSges(6,3) * t4 + 0.2e1 * Ifges(7,3) * t450 - t61 * Ifges(6,4) / 0.2e1 - t172 * Ifges(6,6) / 0.2e1 + (-t501 + Ifges(7,5)) * t449 + (-t499 + Ifges(7,6)) * t435 + (-t451 + t450) * Ifges(6,2);
t452 = Ifges(5,1) * t444 + Ifges(5,4) * t443 + Ifges(5,5) * t435;
t434 = -t184 / 0.2e1;
t433 = -t185 / 0.2e1;
t432 = t185 / 0.2e1;
t428 = -t213 / 0.2e1;
t427 = t214 / 0.2e1;
t417 = pkin(2) * t275;
t416 = pkin(4) * t185;
t415 = pkin(4) * t271;
t413 = pkin(7) * t278;
t412 = pkin(8) * t265;
t409 = g(3) * t265;
t408 = t277 * pkin(4);
t215 = t231 * qJD(2);
t297 = qJ(5) * t217 - qJD(5) * t231;
t329 = qJD(2) * t273;
t210 = t275 * t329 + t356;
t211 = t278 * t329 - t357;
t141 = t210 * t388 + t211 * t386;
t346 = pkin(2) * t358;
t143 = pkin(3) * t215 + pkin(8) * t217 + t346;
t320 = -t141 * t274 + t277 * t143;
t32 = pkin(4) * t215 + t297 * t277 + (-t177 + (qJ(5) * t231 - t164) * t274) * qJD(4) + t320;
t339 = t277 * t141 + t274 * t143 + t164 * t353;
t36 = -qJ(5) * t334 + (-qJD(4) * t181 + t297) * t274 + t339;
t9 = t271 * t32 + t387 * t36;
t111 = t277 * t164 - t181 * t274;
t378 = t231 * t277;
t84 = pkin(4) * t229 - qJ(5) * t378 + t111;
t379 = t231 * t274;
t91 = -qJ(5) * t379 + t112;
t40 = t271 * t84 + t387 * t91;
t405 = mrSges(5,3) * t184;
t404 = mrSges(5,3) * t185;
t403 = Ifges(3,4) * t275;
t402 = Ifges(3,4) * t278;
t401 = Ifges(5,4) * t274;
t400 = Ifges(5,4) * t277;
t397 = t185 * Ifges(5,4);
t396 = t214 * Ifges(4,4);
t272 = -qJ(5) - pkin(8);
t377 = t265 * t272;
t376 = t265 * t279;
t260 = pkin(3) + t408;
t237 = t267 * t260;
t375 = t267 * t279;
t373 = t273 * t279;
t108 = t184 * Ifges(5,2) + t204 * Ifges(5,6) + t397;
t372 = t274 * t108;
t370 = t276 * t264;
t369 = t276 * t266;
t368 = t276 * t274;
t182 = Ifges(5,4) * t184;
t109 = t185 * Ifges(5,1) + t204 * Ifges(5,5) + t182;
t366 = t277 * t109;
t365 = t277 * t279;
t364 = t279 * t264;
t359 = qJD(1) * t278;
t355 = qJD(4) * t231;
t350 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t360) * t413;
t338 = t388 * pkin(2);
t337 = t387 * pkin(4);
t335 = t231 * t354;
t333 = t366 / 0.2e1;
t20 = t60 * mrSges(6,1) + t61 * mrSges(6,2);
t19 = t60 * mrSges(7,1) - t61 * mrSges(7,3);
t330 = -t354 / 0.2e1;
t322 = -t175 * mrSges(4,1) + t176 * mrSges(4,2);
t321 = t363 * t274;
t319 = t279 * t261 - t276 * t273;
t259 = -t338 - pkin(3);
t316 = pkin(3) * t267 + t412;
t312 = mrSges(3,1) * t275 + mrSges(3,2) * t278;
t309 = mrSges(5,1) * t274 + mrSges(5,2) * t277;
t306 = Ifges(5,1) * t277 - t401;
t305 = t278 * Ifges(3,2) + t403;
t304 = -Ifges(5,2) * t274 + t400;
t303 = Ifges(3,5) * t278 - Ifges(3,6) * t275;
t302 = Ifges(5,5) * t277 - Ifges(5,6) * t274;
t300 = pkin(5) * t266 + qJ(6) * t264;
t140 = t210 * t386 - t388 * t211;
t180 = -t388 * t243 - t245 * t386;
t129 = -mrSges(5,2) * t204 + t405;
t130 = mrSges(5,1) * t204 - t404;
t298 = t129 * t277 - t130 * t274;
t139 = pkin(4) * t379 + t180;
t295 = pkin(1) * t312;
t205 = t267 * t368 + t365;
t8 = -t271 * t36 + t32 * t387;
t39 = -t271 * t91 + t387 * t84;
t293 = t217 * t277 + t335;
t292 = t153 * t309;
t291 = t275 * (Ifges(3,1) * t278 - t403);
t240 = t259 - t408;
t98 = pkin(4) * t294 + t140;
t263 = Ifges(3,4) * t359;
t258 = -t337 - pkin(5);
t254 = qJ(6) + t415;
t242 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t359;
t224 = t363 * t277;
t222 = Ifges(3,1) * t360 + Ifges(3,5) * qJD(2) + t263;
t221 = Ifges(3,6) * qJD(2) + qJD(1) * t305;
t208 = t267 * t365 + t368;
t206 = -t267 * t367 + t371;
t203 = Ifges(4,4) * t213;
t193 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t213;
t192 = t266 * t375 + t370;
t191 = t267 * t364 - t369;
t190 = t267 * t369 - t364;
t189 = t266 * t279 + t267 * t370;
t162 = mrSges(4,1) * t213 + mrSges(4,2) * t214;
t158 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t176;
t157 = -qJDD(2) * mrSges(4,2) + mrSges(4,3) * t175;
t156 = t224 * t387 - t271 * t321;
t155 = t224 * t271 + t321 * t387;
t150 = t214 * Ifges(4,1) - t203 + t488;
t149 = -t213 * Ifges(4,2) + t396 + t487;
t148 = t289 * t231;
t147 = t230 * t231;
t144 = -pkin(5) * t289 - t230 * qJ(6) + t240;
t93 = -t217 * t289 - t230 * t355;
t92 = t217 * t230 + t271 * t335 - t325 * t355;
t83 = -mrSges(5,2) * t172 + mrSges(5,3) * t116;
t82 = mrSges(5,1) * t172 - mrSges(5,3) * t115;
t77 = mrSges(6,1) * t119 + mrSges(6,2) * t290;
t76 = mrSges(7,1) * t119 - mrSges(7,3) * t290;
t71 = t147 * pkin(5) - t148 * qJ(6) + t139;
t69 = -mrSges(5,1) * t116 + mrSges(5,2) * t115;
t52 = pkin(5) * t290 + qJ(6) * t119 + t416;
t49 = -qJD(4) * t112 + t320;
t48 = -t181 * t354 + t339;
t46 = t115 * Ifges(5,4) + t116 * Ifges(5,2) + t172 * Ifges(5,6);
t38 = -t229 * pkin(5) - t39;
t37 = qJ(6) * t229 + t40;
t29 = t271 * t78 + t73;
t18 = -t92 * pkin(5) - t93 * qJ(6) - t148 * qJD(6) + t98;
t7 = -t215 * pkin(5) - t8;
t5 = qJ(6) * t215 + qJD(6) * t229 + t9;
t10 = [t278 * t222 * t407 - qJDD(2) * mrSges(3,2) * t413 + Ifges(3,6) * t278 * t389 + (t471 / 0.2e1 + t201 * mrSges(4,1) - t106 * mrSges(4,3) - Ifges(4,4) * t176 + Ifges(5,5) * t444 - Ifges(4,2) * t175 - t472 * Ifges(4,6) + Ifges(5,6) * t443 + Ifges(6,6) * t451 + Ifges(7,6) * t450 + t475 * t435 + t500 * t449 + t462) * t229 + t236 * t402 / 0.2e1 + (t201 * mrSges(4,2) - t105 * mrSges(4,3) + Ifges(4,1) * t176 + Ifges(4,4) * t175 + Ifges(4,5) * t472 + t109 * t330 + t302 * t435 + t304 * t443 + t306 * t444 + t99 * t309) * t231 + (t291 + t278 * (-Ifges(3,2) * t275 + t402)) * t352 / 0.2e1 + (-Ifges(5,1) * t293 - Ifges(5,4) * t294) * t432 + (t473 / 0.2e1 + t494 / 0.2e1 + Ifges(6,6) * t442 + Ifges(5,5) * t432 + Ifges(7,6) * t441 - Ifges(4,4) * t427 - Ifges(4,2) * t428 - Ifges(4,6) * t407 - t149 / 0.2e1 + t459 + t500 * t438 + t475 * t430) * t215 + t378 * t452 + t518 * t92 + t456 * t147 + (m(5) * t373 - t206 * mrSges(5,1) - t205 * mrSges(5,2) - t445 * (pkin(4) * t371 + t276 * t377 - t373) + t469 * t190 + t467 * t189 + (m(4) * t273 + t474) * t279 + (-m(5) * (-t261 - t316) + m(4) * t261 - t445 * (-t261 - t237) + t468 + t514) * t276) * g(1) + m(7) * (t1 * t37 + t18 * t41 + t2 * t38 + t21 * t7 + t22 * t5 + t6 * t71) + m(6) * (t117 * t98 + t139 * t50 + t25 * t8 + t26 * t9 + t3 * t39 + t4 * t40) + t509 * t148 + t184 * (-Ifges(5,4) * t293 - Ifges(5,2) * t294) / 0.2e1 + (-t208 * mrSges(5,1) - t207 * mrSges(5,2) + t522 * t376 - t524 * t319 - t445 * (pkin(4) * t368 + t260 * t375 - t272 * t376 + t319) - t469 * t192 - t467 * t191 + t474 * t276 + (-m(5) * t316 + t478 - t514) * t279) * g(2) - (t333 + Ifges(4,1) * t427 + Ifges(4,4) * t428 + Ifges(4,5) * t407 + t150 / 0.2e1 - t372 / 0.2e1 - t470) * t217 + t507 * t93 + (-Ifges(5,5) * t293 - Ifges(5,6) * t294) * t430 + (t235 * t413 + t477) * mrSges(3,3) + (-m(4) * t105 + m(5) * t99 - t158 + t69) * t180 + (-t23 * t379 - t24 * t378 + t293 * t86 - t294 * t87) * mrSges(5,3) + m(5) * (t111 * t24 + t112 * t23 + t48 * t87 + t49 * t86) + m(4) * (t106 * t181 + t141 * t166 - t201 * t261 + t239 * t346) + t153 * (mrSges(5,1) * t294 - mrSges(5,2) * t293) + Ifges(2,3) * qJDD(1) + t278 * (Ifges(3,4) * t236 + Ifges(3,2) * t235 + Ifges(3,6) * qJDD(2)) / 0.2e1 - pkin(1) * (-mrSges(3,1) * t235 + mrSges(3,2) * t236) - t244 * t385 - t46 * t379 / 0.2e1 + t141 * t193 + t181 * t157 - t221 * t358 / 0.2e1 - t295 * t352 + t139 * t20 + t48 * t129 + t49 * t130 - t242 * t344 + t7 * t104 + t111 * t82 + t112 * t83 + t98 * t77 + t5 * t101 + t9 * t102 + t8 * t103 + t18 * t76 - t108 * t334 / 0.2e1 + t71 * t19 + (t303 * t407 - t350) * qJD(2) + t162 * t346 + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t477) + (-m(4) * t165 + m(5) * t153 - t481) * t140 + (Ifges(3,1) * t236 + Ifges(3,4) * t503 - pkin(7) * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t236) + t472 * Ifges(3,5)) * t275 + t305 * t503 - t261 * t322 + t37 * t42 + t40 * t43 + t39 * t44 + t38 * t45; t513 * (t312 + t524 * t417 - t445 * (-t267 * t272 - t417) + (-m(5) * pkin(8) + t522 + t523) * t267 + (mrSges(4,1) + m(6) * t260 - m(7) * (-t260 - t300) + t517) * t265) + (-m(5) * (t268 + t412) - m(4) * t268 + t244 - t445 * (t237 + t268 - t377) + (-m(7) * t300 - t517) * t267 - t468) * g(3) + (t350 + (t295 - t291 / 0.2e1) * qJD(1)) * qJD(1) + (t333 + t292) * qJD(4) - t456 * t289 - (t306 * t433 + t304 * t434 - t292 - t488 / 0.2e1 + t470) * t213 + t489 * t102 + t490 * t101 + t491 * t104 + (t117 * t482 - t155 * t3 + t156 * t4 + t240 * t50 + t25 * t492 + t26 * t489) * m(6) + t492 * t103 + (t1 * t156 + t144 * t6 + t155 * t2 + t21 * t491 + t22 * t490 + t41 * t493) * m(7) + t493 * t76 - (-Ifges(3,2) * t360 + t222 + t263) * t359 / 0.2e1 + (t184 * t304 + t185 * t306 + t204 * t302) * qJD(4) / 0.2e1 + (Ifges(5,2) * t277 + t401) * t443 + (Ifges(5,1) * t274 + t400) * t444 + t274 * t452 + (Ifges(5,5) * t274 + Ifges(5,6) * t277) * t435 + t149 * t427 + t372 * t428 + t157 * t336 + t158 * t338 - t518 * t212 - (-t495 / 0.2e1 + t506) * t133 - t508 * t132 + t509 * t230 + t482 * t77 + t507 * t216 + (-t130 * t353 - t129 * t354 + m(5) * ((-t87 * t274 - t86 * t277) * qJD(4) + t476) - t274 * t82 + t277 * t83) * t257 + (-t526 * t87 + (-t353 - t383) * t86 + t476) * mrSges(5,3) + (t221 / 0.2e1 + pkin(7) * t242) * t360 + t108 * t330 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + (t165 * t173 - t166 * t174 - t239 * t347 + (t105 * t388 + t106 * t386) * pkin(2)) * m(4) + t277 * t46 / 0.2e1 + t259 * t69 + Ifges(3,6) * t235 + Ifges(3,5) * t236 + t240 * t20 - t226 * mrSges(3,2) - t227 * mrSges(3,1) - t174 * t193 + Ifges(4,6) * t175 + Ifges(4,5) * t176 - t303 * t352 / 0.2e1 + t144 * t19 - t97 * t129 - t162 * t347 - t96 * t130 + t105 * mrSges(4,1) - t106 * mrSges(4,2) + t496 * t155 + t497 * t156 + t99 * t310 - (-Ifges(5,5) * t433 - Ifges(5,6) * t434 - Ifges(6,6) * t441 - Ifges(7,6) * t442 - t487 / 0.2e1 - t500 * t439 + t459) * t214 + t481 * t173 + (-t153 * t173 + t259 * t99 - t86 * t96 - t87 * t97) * m(5) + (t132 * t499 - t213 * t302 + t214 * t475) * t431 + (-Ifges(4,2) * t214 + t150 - t203 + t366) * t213 / 0.2e1 - (-Ifges(4,1) * t213 - t396 + t473) * t214 / 0.2e1; -(-t193 - t298) * t213 + t298 * qJD(4) + t322 + t497 * t230 - t496 * t289 - (t76 + t77 - t481) * t214 + t277 * t82 + t274 * t83 + t483 * t480 + t484 * t479 + (-g(1) * t276 + g(2) * t279) * (t524 + t445) + (t1 * t230 - t2 * t289 + t21 * t480 - t214 * t41 + t22 * t479) * m(7) + (-t117 * t214 + t230 * t4 - t25 * t480 + t26 * t479 + t289 * t3) * m(6) + (-t153 * t214 + t23 * t274 + t24 * t277 + t204 * (-t274 * t86 + t277 * t87)) * m(5) + (t165 * t214 + t166 * t213 + t201) * m(4); (Ifges(5,5) * t184 - Ifges(5,6) * t185) * t431 + t108 * t432 + (Ifges(5,1) * t184 - t397) * t433 + t43 * t415 + t44 * t337 + (-t431 * t499 + t508) * t290 - t483 * t29 + t485 * t101 + (t1 * t254 + t2 * t258 - t21 * t29 + t22 * t485 - t41 * t52) * m(7) + (-t117 * t416 + t25 * t29 - t26 * t30 + (t271 * t4 + t3 * t387) * pkin(4)) * m(6) + (-mrSges(5,2) * t206 + t189 * t469 - t190 * t467 + t205 * t521) * g(2) + (mrSges(5,2) * t208 + t469 * t191 - t467 * t192 - t207 * t521) * g(1) + t471 - t77 * t416 + t254 * t42 + t258 * t45 + (-t506 + t516) * t119 - t153 * (mrSges(5,1) * t185 + mrSges(5,2) * t184) + (t404 + t130) * t87 - t30 * t102 - t52 * t76 + (t405 - t129) * t86 + (t264 * t469 - t266 * t467 + t274 * t515 + t309) * t409 + t462 + (-Ifges(5,2) * t185 + t109 + t182) * t434; t484 * t119 - t483 * t290 + t445 * t267 * g(3) + t19 + t20 + (t119 * t22 - t21 * t290 - t486 + t6) * m(7) + (t119 * t26 + t25 * t290 - t486 + t50) * m(6); -t204 * t101 + t290 * t76 + (-g(1) * t191 - g(2) * t189 - t22 * t204 - t264 * t409 + t290 * t41 + t2) * m(7) + t45;];
tau  = t10;
