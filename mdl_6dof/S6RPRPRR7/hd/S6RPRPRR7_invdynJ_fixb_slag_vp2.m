% Calculate vector of inverse dynamics joint torques for
% S6RPRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 03:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR7_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR7_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR7_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR7_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR7_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:54:58
% EndTime: 2019-03-09 03:55:21
% DurationCPUTime: 14.77s
% Computational Cost: add. (11657->627), mult. (23373->831), div. (0->0), fcn. (16538->14), ass. (0->294)
t242 = qJD(3) + qJD(5);
t247 = sin(pkin(10));
t255 = cos(qJ(3));
t355 = cos(pkin(10));
t306 = t355 * t255;
t251 = sin(qJ(3));
t332 = qJD(1) * t251;
t173 = -qJD(1) * t306 + t247 * t332;
t184 = -t247 * t255 - t251 * t355;
t174 = t184 * qJD(1);
t250 = sin(qJ(5));
t254 = cos(qJ(5));
t280 = t254 * t173 - t174 * t250;
t376 = mrSges(6,3) * t280;
t249 = sin(qJ(6));
t253 = cos(qJ(6));
t95 = t242 * t253 + t249 * t280;
t96 = t242 * t249 - t253 * t280;
t356 = mrSges(6,1) * t242 + mrSges(7,1) * t95 - mrSges(7,2) * t96 + t376;
t257 = -pkin(1) - pkin(7);
t211 = qJD(1) * t257 + qJD(2);
t344 = t211 * t251;
t164 = -qJ(4) * t332 + t344;
t144 = t247 * t164;
t196 = t255 * t211;
t331 = qJD(1) * t255;
t165 = -qJ(4) * t331 + t196;
t148 = qJD(3) * pkin(3) + t165;
t102 = t355 * t148 - t144;
t391 = pkin(8) * t173;
t82 = qJD(3) * pkin(4) + t102 + t391;
t307 = t355 * t164;
t103 = t247 * t148 + t307;
t390 = pkin(8) * t174;
t85 = t103 + t390;
t46 = -t250 * t85 + t254 * t82;
t43 = -pkin(5) * t242 - t46;
t420 = -m(7) * t43 + t356;
t448 = -m(6) * t46 - t420;
t445 = t255 / 0.2e1;
t243 = qJ(3) + pkin(10);
t232 = sin(t243);
t233 = cos(t243);
t238 = t251 * pkin(3);
t297 = t251 * mrSges(4,1) + t255 * mrSges(4,2);
t444 = m(5) * t238 + mrSges(5,1) * t232 + mrSges(5,2) * t233 + t297;
t234 = qJ(5) + t243;
t222 = sin(t234);
t223 = cos(t234);
t443 = mrSges(6,1) * t222 + (mrSges(6,2) - mrSges(7,3)) * t223;
t380 = mrSges(7,2) * t249;
t382 = mrSges(7,1) * t253;
t442 = t380 - t382;
t47 = t250 * t82 + t254 * t85;
t323 = qJD(1) * qJD(3);
t192 = qJDD(1) * t255 - t251 * t323;
t193 = -qJDD(1) * t251 - t255 * t323;
t128 = t192 * t355 + t247 * t193;
t210 = qJDD(1) * t257 + qJDD(2);
t329 = qJD(3) * t255;
t134 = t251 * t210 + t211 * t329;
t322 = qJD(1) * qJD(4);
t104 = qJ(4) * t193 - t251 * t322 + t134;
t330 = qJD(3) * t251;
t133 = t255 * t210 - t211 * t330;
t98 = qJDD(3) * pkin(3) - qJ(4) * t192 - t255 * t322 + t133;
t66 = -t104 * t247 + t355 * t98;
t49 = qJDD(3) * pkin(4) - pkin(8) * t128 + t66;
t127 = -t247 * t192 + t193 * t355;
t67 = t355 * t104 + t247 * t98;
t50 = pkin(8) * t127 + t67;
t15 = -qJD(5) * t47 - t250 * t50 + t254 * t49;
t240 = qJDD(3) + qJDD(5);
t301 = t173 * t250 + t254 * t174;
t61 = qJD(5) * t301 + t127 * t250 + t128 * t254;
t54 = mrSges(6,1) * t240 - mrSges(6,3) * t61;
t441 = m(6) * t15 + t54;
t76 = -pkin(5) * t280 - pkin(9) * t301;
t440 = -m(5) - m(4);
t439 = -m(6) - m(7);
t308 = -pkin(5) * t222 + t223 * pkin(9);
t438 = m(7) * t308;
t44 = pkin(9) * t242 + t47;
t197 = pkin(3) * t332 + qJD(1) * qJ(2) + qJD(4);
t131 = -pkin(4) * t174 + t197;
t63 = -pkin(5) * t301 + pkin(9) * t280 + t131;
t20 = -t249 * t44 + t253 * t63;
t436 = t20 * mrSges(7,1);
t21 = t249 * t63 + t253 * t44;
t435 = t21 * mrSges(7,2);
t377 = mrSges(6,3) * t301;
t100 = -mrSges(6,2) * t242 + t377;
t112 = qJD(6) - t301;
t70 = -mrSges(7,2) * t112 + mrSges(7,3) * t95;
t71 = mrSges(7,1) * t112 - mrSges(7,3) * t96;
t286 = -t249 * t71 + t253 * t70;
t272 = -t100 - t286;
t185 = -t247 * t251 + t306;
t125 = t184 * t250 + t254 * t185;
t256 = cos(qJ(1));
t320 = t223 * t380;
t381 = mrSges(6,2) * t222;
t430 = (-t320 - t381) * t256;
t252 = sin(qJ(1));
t319 = t223 * t382;
t342 = t223 * t252;
t343 = t222 * t252;
t429 = -mrSges(6,1) * t342 - mrSges(7,3) * t343 - (t319 - t320) * t252;
t313 = t355 * pkin(3);
t224 = t313 + pkin(4);
t394 = pkin(3) * t247;
t170 = t250 * t224 + t254 * t394;
t245 = qJD(1) * qJD(2);
t213 = qJDD(1) * qJ(2) + t245;
t428 = -t222 * t442 + t443;
t326 = qJD(6) * t253;
t175 = -qJD(3) * t306 + t247 * t330;
t176 = t184 * qJD(3);
t279 = t254 * t184 - t185 * t250;
t78 = qJD(5) * t279 + t175 * t250 + t254 * t176;
t271 = t125 * t326 + t249 * t78;
t281 = t133 * t255 + t134 * t251;
t35 = qJD(6) * t95 + t240 * t249 + t253 * t61;
t62 = qJD(5) * t280 + t127 * t254 - t128 * t250;
t60 = qJDD(6) - t62;
t17 = mrSges(7,1) * t60 - mrSges(7,3) * t35;
t36 = -qJD(6) * t96 + t240 * t253 - t249 * t61;
t18 = -mrSges(7,2) * t60 + mrSges(7,3) * t36;
t427 = -t249 * t17 + t253 * t18;
t426 = mrSges(6,1) * t223 + t222 * mrSges(7,3) + t319;
t425 = -g(1) * t252 + g(2) * t256;
t298 = mrSges(4,1) * t255 - mrSges(4,2) * t251;
t372 = Ifges(4,4) * t255;
t424 = (-Ifges(4,1) * t251 - t372) * t445 + qJ(2) * t298;
t12 = -pkin(5) * t240 - t15;
t16 = -mrSges(7,1) * t36 + mrSges(7,2) * t35;
t423 = -m(7) * t12 - t16;
t206 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t332;
t207 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t331;
t422 = (t206 * t255 - t207 * t251) * qJD(3);
t421 = -t102 * t176 + t103 * t175 + t184 * t67 - t185 * t66;
t419 = mrSges(2,1) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3) - mrSges(3,2);
t14 = qJD(5) * t46 + t250 * t49 + t254 * t50;
t11 = pkin(9) * t240 + t14;
t138 = -pkin(3) * t193 + qJDD(4) + t213;
t90 = -pkin(4) * t127 + t138;
t19 = -pkin(5) * t62 - pkin(9) * t61 + t90;
t2 = qJD(6) * t20 + t11 * t253 + t19 * t249;
t3 = -qJD(6) * t21 - t11 * t249 + t19 * t253;
t418 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t417 = m(4) * t281 + t255 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t192) + t251 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t193);
t416 = -mrSges(3,3) + t438 + mrSges(2,2) - t443 - t444;
t288 = t20 * t253 + t21 * t249;
t261 = -qJD(6) * t288 + t2 * t253 - t249 * t3;
t327 = qJD(6) * t249;
t415 = m(7) * t261 - t71 * t326 - t70 * t327 + t427;
t397 = -t242 / 0.2e1;
t402 = -t301 / 0.2e1;
t403 = -t112 / 0.2e1;
t406 = -t96 / 0.2e1;
t407 = -t95 / 0.2e1;
t414 = -t131 * mrSges(6,1) + Ifges(7,5) * t406 - Ifges(6,2) * t402 - Ifges(6,6) * t397 + Ifges(7,6) * t407 + Ifges(7,3) * t403 + t435 - t436;
t295 = mrSges(7,1) * t249 + mrSges(7,2) * t253;
t269 = t43 * t295;
t289 = Ifges(7,5) * t253 - Ifges(7,6) * t249;
t368 = Ifges(7,4) * t253;
t291 = -Ifges(7,2) * t249 + t368;
t369 = Ifges(7,4) * t249;
t293 = Ifges(7,1) * t253 - t369;
t92 = Ifges(7,4) * t95;
t42 = t96 * Ifges(7,1) + t112 * Ifges(7,5) + t92;
t360 = t253 * t42;
t374 = mrSges(7,3) * t253;
t375 = mrSges(7,3) * t249;
t396 = t249 / 0.2e1;
t401 = t280 / 0.2e1;
t385 = t96 * Ifges(7,4);
t41 = t95 * Ifges(7,2) + t112 * Ifges(7,6) + t385;
t413 = -t131 * mrSges(6,2) + Ifges(6,1) * t401 + Ifges(6,5) * t397 + t20 * t374 + t21 * t375 + t289 * t403 + t291 * t407 + t293 * t406 - t269 - t360 / 0.2e1 + t41 * t396;
t412 = qJD(1) ^ 2;
t411 = t35 / 0.2e1;
t410 = t36 / 0.2e1;
t408 = t60 / 0.2e1;
t405 = t96 / 0.2e1;
t404 = -m(3) - m(4);
t400 = -t280 / 0.2e1;
t399 = -t173 / 0.2e1;
t393 = pkin(3) * t255;
t195 = pkin(4) * t233 + t393;
t395 = m(6) * t195;
t387 = t46 * mrSges(6,3);
t386 = t47 * mrSges(6,3);
t248 = -qJ(4) - pkin(7);
t379 = mrSges(5,3) * t173;
t378 = mrSges(5,3) * t174;
t373 = Ifges(4,4) * t251;
t371 = Ifges(5,4) * t173;
t370 = Ifges(6,4) * t280;
t352 = t125 * t249;
t351 = t125 * t253;
t341 = t249 * t252;
t340 = t249 * t256;
t338 = t252 * t253;
t337 = t253 * t256;
t225 = qJ(2) + t238;
t335 = qJ(4) - t257;
t160 = -qJD(4) * t255 + t330 * t335;
t199 = t335 * t255;
t161 = -qJD(3) * t199 - qJD(4) * t251;
t106 = t247 * t160 + t355 * t161;
t110 = t355 * t165 - t144;
t198 = t335 * t251;
t130 = -t355 * t198 - t247 * t199;
t334 = pkin(5) * t342 + pkin(9) * t343;
t333 = t256 * pkin(1) + t252 * qJ(2);
t325 = qJDD(1) * mrSges(3,2);
t214 = pkin(3) * t329 + qJD(2);
t324 = -m(5) + t439;
t321 = Ifges(7,5) * t35 + Ifges(7,6) * t36 + Ifges(7,3) * t60;
t228 = pkin(3) * t331;
t314 = t360 / 0.2e1;
t311 = -t62 * mrSges(6,1) + t61 * mrSges(6,2);
t310 = -t327 / 0.2e1;
t237 = t256 * qJ(2);
t309 = -pkin(1) * t252 + t237;
t137 = -pkin(4) * t173 + t228;
t194 = pkin(4) * t232 + t238;
t305 = -t323 / 0.2e1;
t304 = -t127 * mrSges(5,1) + t128 * mrSges(5,2);
t302 = (t213 + t245) * qJ(2);
t105 = t355 * t160 - t161 * t247;
t109 = -t165 * t247 - t307;
t129 = t198 * t247 - t355 * t199;
t149 = -pkin(4) * t184 + t225;
t132 = -pkin(4) * t175 + t214;
t299 = -pkin(5) * t223 - pkin(9) * t222;
t294 = t255 * Ifges(4,1) - t373;
t292 = -t251 * Ifges(4,2) + t372;
t290 = -Ifges(4,5) * t251 - Ifges(4,6) * t255;
t287 = -t20 * t249 + t21 * t253;
t107 = -pkin(8) * t185 + t129;
t108 = pkin(8) * t184 + t130;
t69 = t107 * t250 + t108 * t254;
t72 = -pkin(5) * t279 - pkin(9) * t125 + t149;
t29 = t249 * t72 + t253 * t69;
t28 = -t249 * t69 + t253 * t72;
t285 = -t249 * t70 - t253 * t71;
t283 = t254 * t107 - t108 * t250;
t169 = t224 * t254 - t250 * t394;
t274 = -pkin(8) * t176 + t105;
t273 = t109 - t390;
t270 = t125 * t327 - t253 * t78;
t267 = t251 * (-Ifges(4,2) * t255 - t373);
t260 = qJD(5) * t125 - t254 * t175 + t176 * t250;
t8 = t35 * Ifges(7,4) + t36 * Ifges(7,2) + t60 * Ifges(7,6);
t9 = t35 * Ifges(7,1) + t36 * Ifges(7,4) + t60 * Ifges(7,5);
t259 = -t14 * mrSges(6,2) - t3 * t375 + t2 * t374 + t12 * t442 + t15 * mrSges(6,1) + Ifges(6,3) * t240 + (Ifges(7,1) * t249 + t368) * t411 + (Ifges(7,2) * t253 + t369) * t410 + t41 * t310 + (Ifges(7,5) * t249 + Ifges(7,6) * t253) * t408 + Ifges(6,6) * t62 + Ifges(6,5) * t61 + t9 * t396 + t253 * t8 / 0.2e1 + (-t20 * t326 - t21 * t327) * mrSges(7,3) + (t269 + t314) * qJD(6) + (t112 * t289 + t291 * t95 + t293 * t96) * qJD(6) / 0.2e1;
t241 = -pkin(8) + t248;
t231 = -pkin(1) * qJDD(1) + qJDD(2);
t189 = t297 * qJD(1);
t178 = Ifges(4,5) * qJD(3) + qJD(1) * t294;
t177 = Ifges(4,6) * qJD(3) + qJD(1) * t292;
t166 = Ifges(5,4) * t174;
t162 = -pkin(5) - t169;
t159 = t222 * t337 - t341;
t158 = t222 * t340 + t338;
t157 = t222 * t338 + t340;
t156 = -t222 * t341 + t337;
t143 = qJD(3) * mrSges(5,1) + t379;
t142 = -qJD(3) * mrSges(5,2) + t378;
t121 = -mrSges(5,1) * t174 - mrSges(5,2) * t173;
t120 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t128;
t119 = -qJDD(3) * mrSges(5,2) + mrSges(5,3) * t127;
t114 = -t173 * Ifges(5,1) + Ifges(5,5) * qJD(3) + t166;
t113 = t174 * Ifges(5,2) + Ifges(5,6) * qJD(3) - t371;
t111 = Ifges(6,4) * t301;
t89 = t110 + t391;
t88 = pkin(8) * t175 + t106;
t75 = -mrSges(6,1) * t301 - mrSges(6,2) * t280;
t74 = -Ifges(6,1) * t280 + t242 * Ifges(6,5) + t111;
t73 = Ifges(6,2) * t301 + t242 * Ifges(6,6) - t370;
t64 = t137 + t76;
t55 = -mrSges(6,2) * t240 + mrSges(6,3) * t62;
t52 = t250 * t273 + t254 * t89;
t40 = t96 * Ifges(7,5) + t95 * Ifges(7,6) + t112 * Ifges(7,3);
t32 = pkin(5) * t260 - pkin(9) * t78 + t132;
t27 = t249 * t76 + t253 * t46;
t26 = -t249 * t46 + t253 * t76;
t24 = qJD(5) * t283 + t250 * t274 + t254 * t88;
t23 = t249 * t64 + t253 * t52;
t22 = -t249 * t52 + t253 * t64;
t5 = -qJD(6) * t29 - t24 * t249 + t253 * t32;
t4 = qJD(6) * t28 + t24 * t253 + t249 * t32;
t1 = [-t260 * t435 + t260 * t436 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-t2 * t352 + t20 * t270 - t21 * t271 - t3 * t351) * mrSges(7,3) + (Ifges(4,1) * t192 + Ifges(4,4) * t193) * t445 + t448 * (qJD(5) * t69 + t250 * t88 - t254 * t274) + (t90 * mrSges(6,2) - t15 * mrSges(6,3) + Ifges(6,1) * t61 + Ifges(6,4) * t62 + Ifges(6,5) * t240 + t12 * t295 + t289 * t408 + t291 * t410 + t293 * t411 + t310 * t42) * t125 + m(5) * (t102 * t105 + t103 * t106 + t129 * t66 + t130 * t67 + t138 * t225 + t197 * t214) - t251 * (Ifges(4,4) * t192 + Ifges(4,2) * t193) / 0.2e1 - t281 * mrSges(4,3) - t271 * t41 / 0.2e1 + t424 * t323 - t177 * t329 / 0.2e1 - t178 * t330 / 0.2e1 - pkin(1) * t325 + m(3) * (-pkin(1) * t231 + t302) + t225 * t304 - (-t423 - t441) * t283 + t193 * t292 / 0.2e1 + t192 * t294 / 0.2e1 + qJD(3) ^ 2 * t290 / 0.2e1 + t421 * mrSges(5,3) + t417 * t257 + m(4) * t302 + (-t157 * mrSges(7,1) - t156 * mrSges(7,2) + t439 * (t252 * t194 - t241 * t256 + t333) + (-m(3) + t440) * t333 + (-m(4) * pkin(7) + m(5) * t248 - t419) * t256 + t416 * t252) * g(2) + (-m(3) * t309 - t159 * mrSges(7,1) + t158 * mrSges(7,2) + t439 * (t256 * t194 + t252 * t241 + t309) + t440 * t237 + (-m(5) * (-pkin(1) + t248) - m(4) * t257 + t419) * t252 + t416 * t256) * g(1) + (Ifges(5,1) * t176 + Ifges(5,4) * t175) * t399 + t149 * t311 + t257 * t422 + (-Ifges(7,1) * t270 - Ifges(7,4) * t271 + Ifges(7,5) * t260) * t405 + (Ifges(6,1) * t78 - Ifges(6,4) * t260) * t400 + t95 * (-Ifges(7,4) * t270 - Ifges(7,2) * t271 + Ifges(7,6) * t260) / 0.2e1 + t112 * (-Ifges(7,5) * t270 - Ifges(7,6) * t271 + Ifges(7,3) * t260) / 0.2e1 + t242 * (Ifges(6,5) * t78 - Ifges(6,6) * t260) / 0.2e1 + t131 * (mrSges(6,1) * t260 + mrSges(6,2) * t78) + t260 * t40 / 0.2e1 - t260 * t73 / 0.2e1 - t260 * t386 + t301 * (Ifges(6,4) * t78 - Ifges(6,2) * t260) / 0.2e1 + m(7) * (t2 * t29 + t20 * t5 + t21 * t4 + t28 * t3) + m(6) * (t131 * t132 + t14 * t69 + t149 * t90 + t24 * t47) - (t321 / 0.2e1 + Ifges(7,3) * t408 + Ifges(7,6) * t410 + Ifges(7,5) * t411 - Ifges(6,6) * t240 - Ifges(6,4) * t61 - Ifges(6,2) * t62 + t90 * mrSges(6,1) - t14 * mrSges(6,3) + t418) * t279 + (t297 + 0.2e1 * mrSges(3,3)) * t213 + t43 * (mrSges(7,1) * t271 - mrSges(7,2) * t270) + t267 * t305 + t231 * mrSges(3,2) + t214 * t121 + qJ(2) * (-mrSges(4,1) * t193 + mrSges(4,2) * t192) + t197 * (-mrSges(5,1) * t175 + mrSges(5,2) * t176) + qJD(2) * t189 + t138 * (-mrSges(5,1) * t184 + mrSges(5,2) * t185) + t127 * (Ifges(5,4) * t185 + Ifges(5,2) * t184) + t128 * (Ifges(5,1) * t185 + Ifges(5,4) * t184) + t174 * (Ifges(5,4) * t176 + Ifges(5,2) * t175) / 0.2e1 + qJD(3) * (Ifges(5,5) * t176 + Ifges(5,6) * t175) / 0.2e1 + t176 * t114 / 0.2e1 + t175 * t113 / 0.2e1 + t106 * t142 + t105 * t143 + t129 * t120 + t130 * t119 + t132 * t75 + t24 * t100 + t78 * t74 / 0.2e1 + t4 * t70 + t5 * t71 + t69 * t55 + t28 * t17 + t29 * t18 + (0.2e1 * Ifges(4,5) * t445 + Ifges(5,5) * t185 - Ifges(4,6) * t251 + Ifges(5,6) * t184) * qJDD(3) + t9 * t351 / 0.2e1 - t8 * t352 / 0.2e1 + t78 * t314 - t78 * t387; t325 - t184 * t119 + t185 * t120 - t175 * t142 + t176 * t143 + t356 * t78 - (t16 - t54) * t125 + t422 - t272 * t260 + (qJ(2) * t404 - mrSges(3,3)) * t412 - (qJD(6) * t285 + t427 + t55) * t279 + m(6) * (t125 * t15 - t14 * t279 + t260 * t47 + t46 * t78) - m(5) * t421 + m(3) * t231 + m(7) * (-t12 * t125 + t260 * t287 - t261 * t279 - t43 * t78) + (-m(5) * t197 - m(6) * t131 - m(7) * t288 - t121 - t189 + t285 - t75) * qJD(1) + t425 * (-t324 - t404) + t417; (-t102 * t109 - t103 * t110 - t197 * t228 + (t247 * t67 + t355 * t66) * pkin(3)) * m(5) + (t12 * t162 - t20 * t22 - t21 * t23) * m(7) + t448 * (t170 * qJD(5) - t250 * t89 + t254 * t273) + (m(5) * t393 + mrSges(5,1) * t233 - mrSges(5,2) * t232 + t298) * t425 - t121 * t228 + (t267 / 0.2e1 - t424) * t412 + (Ifges(5,3) + Ifges(4,3)) * qJDD(3) + t177 * t331 / 0.2e1 + t178 * t332 / 0.2e1 + (Ifges(6,4) * t402 + t413 + t387 - t74 / 0.2e1) * t301 + ((m(6) * t47 + m(7) * t287 - t272) * qJD(5) + t441) * t169 + t259 + (-t131 * t137 + t14 * t170 - t47 * t52) * m(6) + t415 * (pkin(9) + t170) + (-m(7) * (t195 * t252 + t334) - (-t381 + t395) * t252 + t429) * g(1) + ((-m(7) * (-t195 + t299) + t395 + t426) * t256 + t430) * g(2) + t119 * t394 + t113 * t399 - (Ifges(5,2) * t173 + t114 + t166) * t174 / 0.2e1 + (m(6) * t194 - m(7) * (-t194 + t308) + t428 + t444) * g(3) - (-Ifges(6,4) * t401 + t386 + t73 / 0.2e1 - t40 / 0.2e1 + t414) * t280 + t102 * t378 + t207 * t344 + t120 * t313 + t290 * t305 + Ifges(4,6) * t193 - t197 * (-mrSges(5,1) * t173 + mrSges(5,2) * t174) + Ifges(4,5) * t192 - qJD(3) * (Ifges(5,5) * t174 + Ifges(5,6) * t173) / 0.2e1 + t170 * t55 + t162 * t16 - t110 * t142 - t109 * t143 + t133 * mrSges(4,1) - t134 * mrSges(4,2) - t137 * t75 + Ifges(5,6) * t127 + Ifges(5,5) * t128 - t52 * t100 - t23 * t70 - t22 * t71 + t66 * mrSges(5,1) - t67 * mrSges(5,2) - t206 * t196 + t173 * (Ifges(5,1) * t174 + t371) / 0.2e1 - t103 * t379; t286 * qJD(6) - t356 * t280 + t272 * t301 - t174 * t142 - t173 * t143 + t253 * t17 + t249 * t18 + t304 + t311 + (g(1) * t256 + g(2) * t252) * t324 + (t112 * t287 + t2 * t249 + t3 * t253 + t280 * t43) * m(7) + (-t280 * t46 - t301 * t47 + t90) * m(6) + (-t102 * t173 - t103 * t174 + t138) * m(5); t259 - m(7) * (t20 * t26 + t21 * t27) + t73 * t400 - t27 * t70 - t26 * t71 + (t377 - t100) * t46 + (t111 + t74) * t402 + (t40 + t370) * t401 + (t428 - t438) * g(3) + ((-m(7) * t299 + t426) * t256 + t430) * g(2) + (-m(7) * t334 + mrSges(6,2) * t343 + t429) * g(1) + t423 * pkin(5) + (-t376 + t420) * t47 - t414 * t280 + t413 * t301 + t415 * pkin(9); -t43 * (mrSges(7,1) * t96 + mrSges(7,2) * t95) + (Ifges(7,1) * t95 - t385) * t406 + t41 * t405 + (Ifges(7,5) * t95 - Ifges(7,6) * t96) * t403 - t20 * t70 + t21 * t71 - g(1) * (mrSges(7,1) * t156 - mrSges(7,2) * t157) - g(2) * (mrSges(7,1) * t158 + mrSges(7,2) * t159) + g(3) * t295 * t223 + (t20 * t95 + t21 * t96) * mrSges(7,3) + t321 + (-Ifges(7,2) * t96 + t42 + t92) * t407 + t418;];
tau  = t1;
