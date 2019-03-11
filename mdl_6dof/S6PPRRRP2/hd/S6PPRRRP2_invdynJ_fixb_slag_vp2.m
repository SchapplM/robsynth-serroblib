% Calculate vector of inverse dynamics joint torques for
% S6PPRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-03-08 18:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PPRRRP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRP2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRP2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRP2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:56:04
% EndTime: 2019-03-08 18:56:29
% DurationCPUTime: 14.26s
% Computational Cost: add. (6276->591), mult. (15963->813), div. (0->0), fcn. (14134->14), ass. (0->282)
t452 = Ifges(6,1) + Ifges(7,1);
t451 = Ifges(7,4) + Ifges(6,5);
t450 = Ifges(6,6) - Ifges(7,6);
t209 = sin(qJ(5));
t456 = pkin(9) * t209;
t355 = cos(pkin(6));
t197 = qJD(1) * t355 + qJD(2);
t205 = sin(pkin(12));
t207 = sin(pkin(6));
t211 = sin(qJ(3));
t208 = cos(pkin(12));
t354 = cos(pkin(7));
t287 = t208 * t354;
t274 = t211 * t287;
t381 = cos(qJ(3));
t221 = (t205 * t381 + t274) * t207;
t206 = sin(pkin(7));
t343 = t206 * t211;
t101 = qJD(1) * t221 + t197 * t343;
t213 = cos(qJ(4));
t210 = sin(qJ(4));
t272 = pkin(4) * t210 - pkin(10) * t213;
t455 = -qJD(5) * t213 * pkin(9) + t272 * qJD(4) - t101;
t322 = qJD(3) * qJD(4);
t185 = qJDD(3) * t210 + t213 * t322;
t212 = cos(qJ(5));
t324 = t212 * qJD(4);
t325 = t210 * qJD(3);
t175 = t209 * t325 - t324;
t330 = qJD(5) * t175;
t117 = qJDD(4) * t209 + t185 * t212 - t330;
t393 = t117 / 0.2e1;
t176 = qJD(4) * t209 + t212 * t325;
t118 = qJD(5) * t176 - t212 * qJDD(4) + t185 * t209;
t391 = t118 / 0.2e1;
t184 = qJDD(3) * t213 - t210 * t322;
t173 = qJDD(5) - t184;
t390 = t173 / 0.2e1;
t388 = t175 / 0.2e1;
t454 = -t176 / 0.2e1;
t323 = t213 * qJD(3);
t198 = qJD(5) - t323;
t384 = t198 / 0.2e1;
t453 = qJD(4) / 0.2e1;
t449 = -Ifges(6,3) - Ifges(7,2);
t418 = -t209 * t450 + t212 * t451;
t362 = Ifges(7,5) * t209;
t365 = Ifges(6,4) * t209;
t416 = t212 * t452 + t362 - t365;
t392 = -t118 / 0.2e1;
t448 = t451 * t390 + (-Ifges(6,4) + Ifges(7,5)) * t391 + t452 * t393;
t380 = pkin(4) * t213;
t192 = -pkin(10) * t210 - pkin(3) - t380;
t241 = t381 * t287;
t234 = t207 * t241;
t229 = qJD(1) * t234;
t345 = t205 * t207;
t306 = qJD(1) * t345;
t307 = t206 * t381;
t215 = -t197 * t307 + t211 * t306 - t229;
t327 = qJD(5) * t212;
t337 = t212 * t213;
t435 = -t210 * t324 * pkin(9) + t192 * t327 + t209 * t455 + t215 * t337;
t329 = qJD(5) * t209;
t332 = qJD(4) * t210;
t339 = t209 * t213;
t434 = -t192 * t329 + t212 * t455 - t215 * t339 + t332 * t456;
t170 = Ifges(6,4) * t175;
t363 = Ifges(7,5) * t175;
t447 = t176 * t452 + t198 * t451 - t170 + t363;
t368 = Ifges(5,4) * t210;
t433 = Ifges(5,2) * t213;
t257 = t368 + t433;
t446 = Ifges(5,6) * t453 + qJD(3) * t257 / 0.2e1 + t449 * t384 + t451 * t454 + t450 * t388;
t445 = -m(7) - m(6);
t444 = t184 / 0.2e1;
t443 = t185 / 0.2e1;
t442 = mrSges(6,3) + mrSges(7,2);
t439 = -pkin(5) * t332 - t434;
t438 = qJ(6) * t332 - qJD(6) * t213 + t435;
t85 = mrSges(6,1) * t173 - mrSges(6,3) * t117;
t86 = -t173 * mrSges(7,1) + t117 * mrSges(7,2);
t437 = t86 - t85;
t87 = -mrSges(6,2) * t173 - mrSges(6,3) * t118;
t88 = -mrSges(7,2) * t118 + mrSges(7,3) * t173;
t436 = t88 + t87;
t246 = pkin(5) * t209 - qJ(6) * t212;
t316 = t207 * t208 * t206;
t146 = -qJD(1) * t316 + t197 * t354;
t347 = t146 * t210;
t98 = qJD(3) * pkin(9) + t101;
t432 = qJD(5) * t246 - qJD(6) * t209 - t347 - (qJD(3) * t246 + t98) * t213;
t266 = -mrSges(5,1) * t213 + mrSges(5,2) * t210;
t431 = -t266 + mrSges(4,1);
t56 = mrSges(6,1) * t118 + mrSges(6,2) * t117;
t430 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t185 + t56;
t318 = mrSges(5,3) * t325;
t428 = -qJD(4) * mrSges(5,1) + mrSges(6,1) * t175 + mrSges(6,2) * t176 + t318;
t369 = mrSges(6,3) * t176;
t140 = mrSges(6,1) * t198 - t369;
t141 = -mrSges(7,1) * t198 + mrSges(7,2) * t176;
t427 = t140 - t141;
t370 = mrSges(6,3) * t175;
t139 = -mrSges(6,2) * t198 - t370;
t142 = -mrSges(7,2) * t175 + mrSges(7,3) * t198;
t426 = -t142 - t139;
t65 = t213 * t98 + t347;
t61 = qJD(4) * pkin(10) + t65;
t89 = qJD(3) * t192 + t215;
t19 = -t209 * t61 + t212 * t89;
t425 = -t19 + qJD(6);
t353 = cos(pkin(11));
t268 = t355 * t353;
t352 = sin(pkin(11));
t219 = t205 * t352 - t208 * t268;
t289 = t207 * t353;
t424 = t206 * t289 + t219 * t354;
t267 = t355 * t352;
t220 = t205 * t353 + t208 * t267;
t288 = t207 * t352;
t423 = -t206 * t288 + t220 * t354;
t422 = -t449 * t210 + t418 * t213;
t421 = t451 * t210 + t416 * t213;
t262 = mrSges(7,1) * t209 - mrSges(7,3) * t212;
t264 = mrSges(6,1) * t209 + mrSges(6,2) * t212;
t64 = t146 * t213 - t210 * t98;
t60 = -qJD(4) * pkin(4) - t64;
t27 = pkin(5) * t175 - qJ(6) * t176 + t60;
t420 = t27 * t262 + t60 * t264;
t419 = t451 * t209 + t450 * t212;
t361 = Ifges(7,5) * t212;
t364 = Ifges(6,4) * t212;
t417 = t452 * t209 - t361 + t364;
t413 = t451 * t117 - t450 * t118 - t449 * t173;
t133 = -mrSges(5,1) * t184 + mrSges(5,2) * t185;
t412 = mrSges(4,1) * qJDD(3) - t133;
t169 = Ifges(7,5) * t176;
t102 = t198 * Ifges(7,6) + t175 * Ifges(7,3) + t169;
t201 = Ifges(5,4) * t323;
t411 = Ifges(5,1) * t325 + Ifges(5,5) * qJD(4) + t209 * t102 + t201;
t177 = t266 * qJD(3);
t410 = qJD(3) * mrSges(4,1) - t177;
t196 = qJDD(1) * t355 + qJDD(2);
t145 = -qJDD(1) * t316 + t196 * t354;
t331 = qJD(4) * t213;
t243 = t207 * t274;
t279 = qJD(3) * t306;
t280 = qJD(3) * t307;
t302 = qJDD(1) * t345;
t53 = qJD(3) * t229 + qJDD(1) * t243 + t196 * t343 + t197 * t280 - t211 * t279 + t381 * t302;
t51 = qJDD(3) * pkin(9) + t53;
t11 = t210 * t145 + t146 * t331 + t213 * t51 - t332 * t98;
t12 = t145 * t213 - t146 * t332 - t210 * t51 - t98 * t331;
t409 = t11 * t213 - t12 * t210;
t333 = qJD(3) * t211;
t305 = t206 * t333;
t54 = -qJD(3) * qJD(1) * t243 + qJDD(1) * t234 + t196 * t307 - t197 * t305 - t211 * t302 - t381 * t279;
t52 = -qJDD(3) * pkin(3) - t54;
t26 = -t184 * pkin(4) - t185 * pkin(10) + t52;
t9 = qJDD(4) * pkin(10) + t11;
t3 = t209 * t26 + t212 * t9 + t89 * t327 - t329 * t61;
t20 = t209 * t89 + t212 * t61;
t4 = -qJD(5) * t20 - t209 * t9 + t212 * t26;
t408 = -t209 * t4 + t212 * t3;
t1 = qJ(6) * t173 + qJD(6) * t198 + t3;
t2 = -pkin(5) * t173 + qJDD(6) - t4;
t407 = t1 * t212 + t2 * t209;
t406 = -t117 * Ifges(7,5) / 0.2e1 - t173 * Ifges(7,6) / 0.2e1 + Ifges(6,4) * t393 + Ifges(6,6) * t390 + (Ifges(7,3) + Ifges(6,2)) * t392;
t247 = pkin(5) * t212 + qJ(6) * t209;
t263 = -t212 * mrSges(7,1) - t209 * mrSges(7,3);
t265 = mrSges(6,1) * t212 - mrSges(6,2) * t209;
t405 = -m(7) * t247 - mrSges(5,1) + t263 - t265;
t404 = m(5) + m(4) + m(3) - t445;
t285 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t284 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t403 = -m(6) * t60 - t428;
t402 = m(6) * t3 + m(7) * t1 + t436;
t401 = m(6) * t4 - m(7) * t2 - t437;
t18 = qJ(6) * t198 + t20;
t400 = m(6) * t20 + m(7) * t18 - t426;
t17 = -pkin(5) * t198 + t425;
t399 = -m(6) * t19 + m(7) * t17 - t427;
t124 = t206 * t219 - t289 * t354;
t154 = t205 * t268 + t208 * t352;
t80 = t154 * t381 - t211 * t424;
t35 = t124 * t210 + t213 * t80;
t125 = t206 * t220 + t288 * t354;
t155 = -t205 * t267 + t208 * t353;
t82 = t155 * t381 - t211 * t423;
t37 = t125 * t210 + t213 * t82;
t290 = t206 * t355;
t123 = t211 * t290 + t221;
t228 = t354 * t355 - t316;
t84 = t123 * t213 + t210 * t228;
t398 = -g(1) * t37 - g(2) * t35 - g(3) * t84;
t127 = mrSges(7,1) * t175 - mrSges(7,3) * t176;
t397 = -m(7) * t27 - t127 + t403;
t396 = -t4 * mrSges(6,1) + t2 * mrSges(7,1) + t3 * mrSges(6,2) - t1 * mrSges(7,3);
t10 = -qJDD(4) * pkin(4) - t12;
t5 = pkin(5) * t118 - qJ(6) * t117 - qJD(6) * t176 + t10;
t55 = mrSges(7,1) * t118 - mrSges(7,3) * t117;
t395 = -m(5) * t12 + m(6) * t10 + m(7) * t5 + t430 + t55;
t394 = -m(5) * t64 - t397;
t214 = qJD(3) ^ 2;
t389 = -t175 / 0.2e1;
t386 = t176 / 0.2e1;
t372 = -qJD(3) / 0.2e1;
t371 = qJD(5) / 0.2e1;
t367 = Ifges(5,4) * t213;
t366 = Ifges(6,4) * t176;
t360 = t10 * t210;
t79 = t154 * t211 + t381 * t424;
t357 = t210 * t79;
t81 = t155 * t211 + t381 * t423;
t356 = t210 * t81;
t182 = t272 * qJD(3);
t41 = t209 * t182 + t212 * t64;
t351 = mrSges(4,2) * qJD(3);
t242 = t381 * t290;
t344 = t205 * t211;
t122 = t207 * t344 - t234 - t242;
t348 = t122 * t210;
t346 = t192 * t212;
t105 = -t175 * Ifges(6,2) + t198 * Ifges(6,6) + t366;
t341 = t209 * t105;
t338 = t210 * t212;
t148 = pkin(9) * t337 + t209 * t192;
t335 = mrSges(4,2) * qJDD(3);
t328 = qJD(5) * t210;
t320 = pkin(10) * t329;
t319 = pkin(10) * t327;
t317 = mrSges(5,3) * t323;
t313 = -t79 * pkin(3) + pkin(9) * t80;
t312 = -t81 * pkin(3) + pkin(9) * t82;
t304 = -t341 / 0.2e1;
t298 = t331 / 0.2e1;
t295 = -t328 / 0.2e1;
t294 = t327 / 0.2e1;
t292 = -t323 / 0.2e1;
t291 = -t122 * pkin(3) + pkin(9) * t123;
t286 = t322 / 0.2e1;
t281 = t209 * t307;
t256 = -Ifges(6,2) * t209 + t364;
t255 = Ifges(6,2) * t212 + t365;
t252 = Ifges(5,5) * t213 - Ifges(5,6) * t210;
t249 = Ifges(7,3) * t209 + t361;
t248 = -Ifges(7,3) * t212 + t362;
t245 = t122 * t212 - t209 * t84;
t39 = t122 * t209 + t212 * t84;
t40 = t182 * t212 - t209 * t64;
t240 = pkin(9) + t246;
t97 = -qJD(3) * pkin(3) + t215;
t237 = t97 * (mrSges(5,1) * t210 + mrSges(5,2) * t213);
t235 = t210 * (Ifges(5,1) * t213 - t368);
t157 = t210 * t354 + t213 * t343;
t233 = -t209 * t157 - t212 * t307;
t232 = -t209 * t328 + t213 * t324;
t231 = t209 * t331 + t210 * t327;
t156 = t210 * t343 - t213 * t354;
t224 = Ifges(6,6) * t210 + t213 * t256;
t223 = Ifges(7,6) * t210 + t213 * t249;
t83 = t123 * t210 - t213 * t228;
t194 = -qJD(4) * mrSges(5,2) + t317;
t189 = -pkin(4) - t247;
t151 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t184;
t150 = t240 * t210;
t147 = -pkin(9) * t339 + t346;
t138 = -t346 + (pkin(5) + t456) * t213;
t137 = -qJ(6) * t213 + t148;
t129 = -qJD(4) * t156 + t213 * t280;
t126 = pkin(5) * t176 + qJ(6) * t175;
t114 = t123 * qJD(3);
t113 = (t242 + (t241 - t344) * t207) * qJD(3);
t91 = (qJD(5) * t247 - qJD(6) * t212) * t210 + t240 * t331;
t36 = t125 * t213 - t210 * t82;
t34 = t124 * t213 - t210 * t80;
t33 = -pkin(5) * t325 - t40;
t30 = qJ(6) * t325 + t41;
t29 = -qJD(4) * t83 + t113 * t213;
t15 = t209 * t37 - t81 * t212;
t13 = t209 * t35 - t79 * t212;
t6 = [m(5) * (t11 * t84 + t29 * t65) + m(4) * (t101 * t113 + t53 * t123 + t145 * t228) + m(3) * (t196 * t355 + (t205 ^ 2 + t208 ^ 2) * t207 ^ 2 * qJDD(1)) - t123 * t335 + t29 * t194 + t84 * t151 - t113 * t351 + m(2) * qJDD(1) + t400 * (qJD(5) * t245 + t114 * t209 + t212 * t29) + t399 * (qJD(5) * t39 - t114 * t212 + t209 * t29) + t402 * t39 + t401 * t245 + (-m(4) * t54 + m(5) * t52 - t412) * t122 + (m(4) * t215 + m(5) * t97 - t410) * t114 + t395 * t83 + t394 * (qJD(4) * t84 + t113 * t210) + (-m(2) - t404) * g(3); m(3) * t196 + t177 * t305 + m(4) * (t145 * t354 + (t381 * t54 + t211 * t53 + (t381 * t101 + t211 * t215) * qJD(3)) * t206) + m(5) * (t11 * t157 + t65 * t129 + (t333 * t97 - t381 * t52) * t206) + t129 * t194 + t157 * t151 + t399 * (-qJD(5) * t281 + t129 * t209 + t157 * t327 - t212 * t305) + t400 * (qJD(5) * t233 + t212 * t129 + t209 * t305) + (-mrSges(4,1) * t214 - t335) * t343 + t402 * (t212 * t157 - t281) + t401 * t233 + (-mrSges(4,2) * t214 + t412) * t307 + t395 * t156 + t394 * (qJD(4) * t157 + t210 * t280) + (-g(1) * t288 + g(2) * t289 - g(3) * t355) * t404; t338 * t448 + t367 * t443 + t257 * t444 + (-g(1) * t82 - g(2) * t80 - g(3) * t123 - t331 * t64 + t409) * mrSges(5,3) + ((-t210 * t64 + t213 * t65) * m(5) - t351 + t194 * t213 - t397 * t210) * t215 + (t17 * t232 - t18 * t231 + t2 * t338) * mrSges(7,2) + (-m(5) * t291 + mrSges(4,2) * t123 + t445 * (-pkin(10) * t348 - t122 * t380 + t291) - t285 * (-t122 * t337 + t123 * t209) + t284 * (-t122 * t339 - t123 * t212) + t442 * t348 + t431 * t122) * g(3) + (-m(5) * t312 + mrSges(4,2) * t82 + t431 * t81 + t442 * t356 + t445 * (-pkin(10) * t356 - t81 * t380 + t312) - t285 * (t209 * t82 - t337 * t81) + t284 * (-t82 * t212 - t339 * t81)) * g(1) + (-m(5) * t313 + mrSges(4,2) * t80 + t431 * t79 + t442 * t357 + t445 * (-pkin(10) * t357 - t79 * t380 + t313) - t285 * (t209 * t80 - t337 * t79) + t284 * (-t80 * t212 - t339 * t79)) * g(2) + (Ifges(5,1) * t185 + Ifges(5,4) * t444 + t102 * t294 + t249 * t391 + t256 * t392 + t5 * t262 - t286 * t433 + t390 * t418 + t393 * t416 + (-t1 * mrSges(7,2) - t3 * mrSges(6,3) - t406) * t209 + Ifges(5,5) * qJDD(4)) * t210 + (-t248 * t388 - t255 * t389 - t384 * t419 - t386 * t417) * t328 + t410 * t101 + t411 * t298 + t264 * t360 + (t147 * t4 + t148 * t3 + t434 * t19 + t435 * t20) * m(6) + (((-t210 * t65 - t213 * t64) * qJD(4) + t409) * m(5) + t151 * t213 + t428 * t331 + t430 * t210 - t194 * t332 + (t331 * t60 + t360) * m(6)) * pkin(9) + (-t19 * t232 - t20 * t231 - t4 * t338) * mrSges(6,3) + t52 * t266 + (t19 * mrSges(6,1) - t17 * mrSges(7,1) - t20 * mrSges(6,2) - t65 * mrSges(5,3) + t18 * mrSges(7,3) - t446) * t332 + t304 * t331 + t438 * t142 + t439 * t141 + (t1 * t137 + t138 * t2 + t150 * t5 + t17 * t439 + t18 * t438 + t27 * t91) * m(7) + t434 * t140 + t435 * t139 + t60 * (mrSges(6,1) * t231 + mrSges(6,2) * t232) + t27 * (mrSges(7,1) * t231 - mrSges(7,3) * t232) + (-t413 / 0.2e1 - Ifges(7,6) * t391 - Ifges(6,6) * t392 + t367 * t286 + Ifges(5,4) * t443 + Ifges(5,2) * t444 - t451 * t393 + t449 * t390 + t396 + Ifges(5,6) * qJDD(4)) * t213 + t212 * t105 * t295 + t447 * (t209 * t295 + t212 * t298) + t235 * t286 + (-pkin(3) * t52 - t101 * t97) * m(5) + Ifges(4,3) * qJDD(3) + t147 * t85 + t148 * t87 + t150 * t55 - pkin(3) * t133 + t137 * t88 + t138 * t86 + t91 * t127 + (t223 * t388 + t224 * t389 + t252 * t453 + t422 * t384 + t421 * t386 + t237) * qJD(4) - t53 * mrSges(4,2) + t54 * mrSges(4,1); t209 * t448 + (-t237 - t19 * (mrSges(6,1) * t210 - mrSges(6,3) * t337) - t17 * (-mrSges(7,1) * t210 + mrSges(7,2) * t337) - t20 * (-mrSges(6,2) * t210 - mrSges(6,3) * t339) - t18 * (-mrSges(7,2) * t339 + mrSges(7,3) * t210) + (t224 / 0.2e1 - t223 / 0.2e1) * t175) * qJD(3) + (t201 + t411) * t292 + (t371 * t418 + t372 * t422) * t198 + (t371 * t416 + t372 * t421) * t176 + (mrSges(5,2) * t84 + t445 * (-t83 * pkin(4) + pkin(10) * t84) - t405 * t83) * g(3) + (mrSges(5,2) * t37 + t445 * (t36 * pkin(4) + pkin(10) * t37) + t405 * t36) * g(1) + (mrSges(5,2) * t35 + t445 * (t34 * pkin(4) + pkin(10) * t35) + t405 * t34) * g(2) + t419 * t390 + (t341 / 0.2e1 - t420) * t323 + (t304 + t420) * qJD(5) + (t319 - t33) * t141 + (t17 * t327 - t18 * t329 + t398 + t407) * mrSges(7,2) + (-pkin(4) * t10 + ((-t19 * t212 - t20 * t209) * qJD(5) + t408) * pkin(10) - t19 * t40 - t20 * t41) * m(6) + (-t19 * t327 - t20 * t329 + t398 + t408) * mrSges(6,3) + t417 * t393 + t406 * t212 + (t318 + t403) * t65 - t10 * t265 + t5 * t263 + (-Ifges(5,2) * t292 + t446) * t325 + t436 * pkin(10) * t212 + t437 * pkin(10) * t209 + (-t320 - t41) * t139 + (-t319 - t40) * t140 + (-t256 / 0.2e1 + t249 / 0.2e1) * t330 - t214 * t235 / 0.2e1 + (-t17 * t33 - t18 * t30 + t189 * t5 + ((t17 * t212 - t18 * t209) * qJD(5) + t407) * pkin(10) + t432 * t27) * m(7) + t432 * t127 + t447 * (t212 * t292 + t294) + (-t320 - t30) * t142 + t102 * t329 / 0.2e1 - t252 * t322 / 0.2e1 + Ifges(5,3) * qJDD(4) + Ifges(5,6) * t184 + Ifges(5,5) * t185 + t189 * t55 + t248 * t391 + t255 * t392 + (t317 - t194) * t64 - t11 * mrSges(5,2) + t12 * mrSges(5,1) - pkin(4) * t56; (-t452 * t175 + t102 + t169 - t366) * t454 + (-t245 * t285 + t284 * t39) * g(3) - (-t175 * t451 - t450 * t176) * t198 / 0.2e1 + (-Ifges(6,2) * t176 - t170 + t447) * t388 - t396 + (t369 + t427) * t20 + (-t370 + t426) * t19 + (t17 * t175 + t176 * t18) * mrSges(7,2) + (t284 * (t209 * t79 + t212 * t35) + t285 * t13) * g(2) + (t284 * (t209 * t81 + t212 * t37) + t285 * t15) * g(1) - t27 * (mrSges(7,1) * t176 + mrSges(7,3) * t175) - t60 * (mrSges(6,1) * t176 - mrSges(6,2) * t175) + qJD(6) * t142 - t126 * t127 + t105 * t386 + (Ifges(7,3) * t176 - t363) * t389 + (-pkin(5) * t2 + qJ(6) * t1 - t126 * t27 - t17 * t20 + t18 * t425) * m(7) - pkin(5) * t86 + qJ(6) * t88 + t413; t176 * t127 - t198 * t142 + (-g(1) * t15 - g(2) * t13 + g(3) * t245 + t27 * t176 - t18 * t198 + t2) * m(7) + t86;];
tau  = t6;
