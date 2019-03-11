% Calculate vector of inverse dynamics joint torques for
% S6RRPPPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
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
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPPR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:13:53
% EndTime: 2019-03-09 08:14:23
% DurationCPUTime: 20.44s
% Computational Cost: add. (5214->694), mult. (10431->870), div. (0->0), fcn. (6056->10), ass. (0->312)
t226 = sin(pkin(9));
t227 = cos(pkin(9));
t235 = cos(qJ(6));
t329 = qJD(6) * t235;
t232 = sin(qJ(6));
t330 = qJD(6) * t232;
t126 = -t226 * t329 - t227 * t330;
t149 = t226 * t235 + t227 * t232;
t233 = sin(qJ(2));
t337 = qJD(1) * t233;
t99 = t149 * t337;
t356 = -t99 + t126;
t308 = t226 * t337;
t343 = t235 * t227;
t100 = -t232 * t308 + t337 * t343;
t125 = t226 * t330 - t227 * t329;
t341 = t100 - t125;
t445 = Ifges(3,1) + Ifges(6,3);
t464 = Ifges(5,1) + Ifges(4,3);
t463 = Ifges(4,4) + Ifges(3,5);
t443 = Ifges(3,6) - Ifges(4,6);
t462 = Ifges(4,6) - Ifges(5,5);
t461 = -Ifges(5,5) - t443;
t379 = -pkin(3) - qJ(5);
t218 = -pkin(2) + t379;
t231 = -pkin(8) - qJ(5);
t238 = -pkin(2) - pkin(3);
t460 = -m(7) * (t231 + t238) + mrSges(7,3) - m(6) * t218 + mrSges(6,3) - m(5) * t238 - mrSges(5,2);
t325 = qJD(1) * qJD(2);
t306 = t233 * t325;
t236 = cos(qJ(2));
t323 = t236 * qJDD(1);
t159 = t306 - t323;
t101 = -qJDD(2) * t227 - t159 * t226;
t160 = qJDD(1) * t233 + t236 * t325;
t206 = t233 * qJD(3);
t225 = qJDD(1) * pkin(1);
t71 = t159 * pkin(2) - t160 * qJ(3) - qJD(1) * t206 - t225;
t264 = qJDD(4) - t71;
t331 = qJD(5) * t236;
t28 = pkin(4) * t160 + qJD(1) * t331 + t159 * t379 + t264;
t141 = t160 * pkin(7);
t298 = qJDD(3) + t141;
t333 = qJD(4) * t233;
t244 = -qJ(4) * t160 - qJD(1) * t333 + t298;
t52 = -qJD(2) * qJD(5) + qJDD(2) * t218 + t244;
t13 = t226 * t28 + t227 * t52;
t11 = pkin(8) * t101 + t13;
t102 = -qJDD(2) * t226 + t159 * t227;
t12 = -t226 * t52 + t227 * t28;
t7 = pkin(5) * t160 - pkin(8) * t102 + t12;
t328 = t226 * qJD(2);
t336 = qJD(1) * t236;
t143 = -t227 * t336 - t328;
t435 = t233 * pkin(4) + t236 * qJ(5);
t135 = -qJD(1) * pkin(1) - pkin(2) * t336 - qJ(3) * t337;
t96 = pkin(3) * t336 + qJD(4) - t135;
t81 = qJD(1) * t435 + t96;
t189 = qJ(4) * t337;
t198 = pkin(7) * t337;
t327 = qJD(3) + t198;
t307 = -t189 + t327;
t94 = qJD(2) * t218 + t307;
t37 = -t226 * t94 + t227 * t81;
t27 = pkin(5) * t337 - pkin(8) * t143 + t37;
t142 = -qJD(2) * t227 + t226 * t336;
t38 = t226 * t81 + t227 * t94;
t32 = pkin(8) * t142 + t38;
t9 = -t232 * t32 + t235 * t27;
t1 = qJD(6) * t9 + t11 * t235 + t232 * t7;
t459 = t1 * mrSges(7,2);
t10 = t232 * t27 + t235 * t32;
t2 = -qJD(6) * t10 - t11 * t232 + t235 * t7;
t458 = t2 * mrSges(7,1);
t446 = t236 / 0.2e1;
t234 = sin(qJ(1));
t457 = g(2) * t234;
t456 = Ifges(5,6) + t463;
t284 = t236 * mrSges(4,1) + t233 * mrSges(4,3);
t286 = mrSges(3,1) * t236 - mrSges(3,2) * t233;
t455 = -t284 - t286;
t211 = t233 * mrSges(5,1);
t434 = -t236 * mrSges(5,2) + t211;
t454 = -t236 * mrSges(7,3) - t434;
t148 = t226 * t232 - t343;
t453 = t1 * t148 - t10 * t356 + t149 * t2 + t341 * t9;
t199 = pkin(7) * t336;
t157 = -qJ(4) * t336 + t199;
t223 = qJD(2) * qJ(3);
t128 = -t157 - t223;
t106 = qJD(2) * pkin(4) + qJD(5) - t128;
t161 = -qJD(2) * pkin(2) + t327;
t167 = t199 + t223;
t371 = Ifges(6,4) * t227;
t274 = -Ifges(6,2) * t226 + t371;
t372 = Ifges(6,4) * t226;
t278 = Ifges(6,1) * t227 - t372;
t281 = t226 * mrSges(6,1) + t227 * mrSges(6,2);
t348 = t227 * t233;
t350 = t226 * t233;
t359 = t236 * mrSges(4,3);
t361 = t236 * mrSges(5,1);
t452 = m(4) * pkin(7) * (t161 * t236 - t167 * t233) + t106 * t233 * t281 + t38 * (-mrSges(6,2) * t236 - mrSges(6,3) * t350) + t37 * (mrSges(6,1) * t236 - mrSges(6,3) * t348) + t135 * (mrSges(4,1) * t233 - t359) + t96 * (t233 * mrSges(5,2) + t361) + t142 * (Ifges(6,6) * t236 + t233 * t274) / 0.2e1 + t143 * (Ifges(6,5) * t236 + t233 * t278) / 0.2e1;
t373 = Ifges(5,4) * t236;
t275 = -t233 * Ifges(5,2) - t373;
t451 = t10 * mrSges(7,2) - Ifges(5,6) * qJD(2) / 0.2e1 + qJD(1) * t275 / 0.2e1 - t9 * mrSges(7,1);
t291 = t235 * t142 - t143 * t232;
t25 = qJD(6) * t291 + t101 * t232 + t102 * t235;
t410 = t25 / 0.2e1;
t80 = t142 * t232 + t143 * t235;
t26 = -qJD(6) * t80 + t101 * t235 - t102 * t232;
t409 = t26 / 0.2e1;
t450 = -m(4) - m(6);
t449 = m(6) + m(7);
t448 = -m(7) - m(5);
t147 = qJDD(6) + t160;
t398 = t147 / 0.2e1;
t447 = -t160 / 0.2e1;
t263 = pkin(5) * t236 - pkin(8) * t348;
t191 = qJ(3) * t336;
t254 = pkin(4) * t236 + t218 * t233;
t87 = qJD(1) * t254 + t191;
t50 = -t157 * t226 + t227 * t87;
t39 = qJD(1) * t263 + t50;
t51 = t227 * t157 + t226 * t87;
t44 = -pkin(8) * t308 + t51;
t377 = pkin(8) - t218;
t152 = t377 * t226;
t153 = t377 * t227;
t84 = t152 * t235 + t153 * t232;
t442 = qJD(5) * t148 + qJD(6) * t84 - t232 * t39 - t235 * t44;
t85 = t152 * t232 - t153 * t235;
t441 = qJD(5) * t149 - qJD(6) * t85 + t232 * t44 - t235 * t39;
t312 = mrSges(5,3) * t336;
t440 = qJD(2) * mrSges(5,1) - mrSges(6,1) * t142 + mrSges(6,2) * t143 - t312;
t113 = t148 * t236;
t237 = cos(qJ(1));
t426 = g(1) * t237 + t457;
t439 = t236 * t426;
t187 = t227 * pkin(5) + pkin(4);
t178 = qJ(3) + t187;
t315 = mrSges(4,2) * t337;
t438 = mrSges(3,3) * t337 + t315 + (-mrSges(3,1) - mrSges(4,1)) * qJD(2);
t314 = mrSges(4,2) * t336;
t170 = qJD(2) * mrSges(4,3) + t314;
t34 = -mrSges(7,1) * t291 + mrSges(7,2) * t80;
t320 = t34 + t440;
t437 = -t170 - t320;
t436 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t336 + t170;
t433 = qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t431 = -t233 * t443 + t236 * t463;
t195 = pkin(7) * t323;
t140 = -pkin(7) * t306 + t195;
t430 = t140 * t236 + t141 * t233;
t105 = t140 + t433;
t114 = -qJDD(2) * pkin(2) + t298;
t429 = t105 * t236 + t114 * t233;
t69 = -mrSges(6,2) * t160 + mrSges(6,3) * t101;
t70 = mrSges(6,1) * t160 - mrSges(6,3) * t102;
t428 = -t226 * t70 + t227 * t69;
t97 = -mrSges(6,2) * t337 + mrSges(6,3) * t142;
t98 = mrSges(6,1) * t337 - mrSges(6,3) * t143;
t427 = t226 * t98 - t227 * t97;
t270 = -t12 * t226 + t13 * t227;
t425 = -mrSges(2,1) + t455;
t369 = Ifges(4,5) * t236;
t374 = Ifges(5,4) * t233;
t424 = t233 * t374 + (t369 - t373 + (-Ifges(5,2) + t464) * t233) * t236;
t183 = qJD(6) + t337;
t197 = Ifges(3,4) * t336;
t280 = t233 * Ifges(4,1) - t369;
t423 = t143 * Ifges(6,5) + t80 * Ifges(7,5) + t142 * Ifges(6,6) + Ifges(7,6) * t291 + t183 * Ifges(7,3) + qJD(1) * t280 + t463 * qJD(2) + t337 * t445 + t197;
t335 = qJD(2) * t233;
t318 = pkin(7) * t335;
t332 = qJD(4) * t236;
t261 = t318 + t332;
t65 = -t159 * qJ(4) + qJD(1) * t261 - t195 - t433;
t219 = pkin(9) + qJ(6);
t203 = sin(t219);
t204 = cos(t219);
t251 = m(7) * t187 + mrSges(7,1) * t204 - mrSges(7,2) * t203;
t282 = t227 * mrSges(6,1) - t226 * mrSges(6,2);
t253 = m(6) * pkin(4) + t282;
t422 = -t359 - t361 + (-t251 - t253) * t236 + (m(4) * pkin(2) + mrSges(4,1) + t460) * t233;
t420 = qJ(3) * (t448 + t450);
t196 = Ifges(4,5) * t337;
t279 = -t236 * Ifges(5,1) - t374;
t419 = t227 * (t143 * Ifges(6,1) + t142 * Ifges(6,4) + Ifges(6,5) * t337) - Ifges(4,3) * t336 + t196 + qJD(1) * t279 + t462 * qJD(2);
t418 = m(5) * t128 - m(6) * t106;
t76 = -pkin(5) * t142 + t106;
t416 = m(7) * t76 - t418;
t390 = pkin(5) * t226;
t415 = m(6) * qJ(4) + m(7) * t390 + mrSges(2,2) - mrSges(4,2) - mrSges(3,3) + mrSges(5,3) + t281;
t411 = Ifges(7,4) * t410 + Ifges(7,2) * t409 + Ifges(7,6) * t398;
t391 = Ifges(7,4) * t80;
t30 = Ifges(7,2) * t291 + Ifges(7,6) * t183 + t391;
t408 = t30 / 0.2e1;
t75 = Ifges(7,4) * t291;
t31 = Ifges(7,1) * t80 + Ifges(7,5) * t183 + t75;
t407 = t31 / 0.2e1;
t406 = -t102 * Ifges(6,1) / 0.2e1 - t101 * Ifges(6,4) / 0.2e1 + Ifges(6,5) * t447;
t405 = -t291 / 0.2e1;
t404 = t291 / 0.2e1;
t403 = -t80 / 0.2e1;
t402 = t80 / 0.2e1;
t400 = t101 / 0.2e1;
t399 = t102 / 0.2e1;
t396 = t160 / 0.2e1;
t395 = -t183 / 0.2e1;
t394 = t183 / 0.2e1;
t392 = -t236 / 0.2e1;
t389 = pkin(7) * t233;
t388 = pkin(7) * t236;
t385 = g(3) * t236;
t215 = t236 * pkin(2);
t378 = pkin(7) - qJ(4);
t230 = qJ(3) + pkin(4);
t376 = Ifges(3,4) * t233;
t375 = Ifges(3,4) * t236;
t370 = Ifges(4,5) * t233;
t334 = qJD(2) * t236;
t122 = t334 * t378 - t333;
t340 = qJ(3) * t334 + t206;
t66 = qJD(2) * t254 + t331 + t340;
t41 = t227 * t122 + t226 * t66;
t171 = t378 * t233;
t207 = t233 * qJ(3);
t339 = t215 + t207;
t311 = t236 * pkin(3) + t339;
t268 = t311 + t435;
t92 = pkin(1) + t268;
t64 = t227 * t171 + t226 * t92;
t349 = t226 * t236;
t347 = t227 * t236;
t346 = t231 * t236;
t345 = t233 * t234;
t344 = t233 * t237;
t342 = t236 * t237;
t118 = qJDD(2) * mrSges(5,2) - t160 * mrSges(5,3);
t154 = t198 - t189;
t338 = t237 * pkin(1) + t234 * pkin(7);
t326 = m(5) + t449;
t319 = Ifges(7,5) * t25 + Ifges(7,6) * t26 + Ifges(7,3) * t147;
t316 = m(4) + t326;
t313 = mrSges(5,3) * t337;
t310 = -pkin(7) + t390;
t309 = t238 * qJD(2);
t8 = -t26 * mrSges(7,1) + t25 * mrSges(7,2);
t303 = t337 / 0.2e1;
t54 = -t101 * mrSges(6,1) + t102 * mrSges(6,2);
t297 = -pkin(1) - t207;
t296 = -t325 / 0.2e1;
t295 = t325 / 0.2e1;
t294 = t160 * mrSges(5,1) + t159 * mrSges(5,2);
t40 = -t122 * t226 + t227 * t66;
t63 = -t171 * t226 + t227 * t92;
t117 = -qJDD(2) * mrSges(4,1) + t160 * mrSges(4,2);
t289 = pkin(2) * t342 + qJ(3) * t344 + t338;
t285 = mrSges(3,1) * t233 + mrSges(3,2) * t236;
t277 = t236 * Ifges(3,2) + t376;
t272 = Ifges(5,5) * t233 - Ifges(5,6) * t236;
t271 = Ifges(6,5) * t227 - Ifges(6,6) * t226;
t45 = pkin(5) * t233 + pkin(8) * t347 + t63;
t53 = pkin(8) * t349 + t64;
t16 = -t232 * t53 + t235 * t45;
t17 = t232 * t45 + t235 * t53;
t269 = pkin(3) * t342 + t289;
t262 = pkin(1) * t285;
t258 = t233 * (Ifges(3,1) * t236 - t376);
t112 = t149 * t236;
t241 = t233 * (Ifges(6,3) * t236 + t233 * t271);
t57 = qJDD(2) * pkin(4) + qJDD(5) - t65;
t216 = t237 * pkin(7);
t210 = t236 * qJ(4);
t192 = qJ(4) * t335;
t164 = qJD(2) * mrSges(5,2) - t313;
t162 = -pkin(1) - t339;
t158 = t434 * qJD(1);
t156 = pkin(2) * t337 - t191;
t155 = t284 * qJD(1);
t137 = pkin(1) + t311;
t131 = Ifges(3,6) * qJD(2) + qJD(1) * t277;
t124 = -t236 * t310 - t210;
t121 = pkin(2) * t335 - t340;
t119 = -mrSges(4,2) * t159 + qJDD(2) * mrSges(4,3);
t116 = -qJDD(2) * mrSges(5,1) - mrSges(5,3) * t159;
t115 = t238 * t337 + t191;
t111 = t309 + t307;
t110 = -t203 * t234 + t204 * t344;
t109 = -t203 * t344 - t204 * t234;
t108 = -t203 * t237 - t204 * t345;
t107 = t203 * t345 - t204 * t237;
t103 = pkin(5) * t308 - t154;
t95 = t233 * t309 + t340;
t91 = t310 * t335 + t192 - t332;
t73 = t143 * Ifges(6,4) + t142 * Ifges(6,2) + Ifges(6,6) * t337;
t62 = qJDD(2) * t238 + t244;
t61 = -qJD(6) * t113 - t149 * t335;
t60 = qJD(6) * t112 - t148 * t335;
t56 = mrSges(7,1) * t183 - mrSges(7,3) * t80;
t55 = -mrSges(7,2) * t183 + mrSges(7,3) * t291;
t46 = -pkin(3) * t159 + t264;
t42 = t102 * Ifges(6,4) + t101 * Ifges(6,2) + t160 * Ifges(6,6);
t36 = -pkin(5) * t101 + t57;
t35 = -pkin(8) * t233 * t328 + t41;
t33 = qJD(2) * t263 + t40;
t19 = -mrSges(7,2) * t147 + mrSges(7,3) * t26;
t18 = mrSges(7,1) * t147 - mrSges(7,3) * t25;
t6 = t25 * Ifges(7,1) + t26 * Ifges(7,4) + t147 * Ifges(7,5);
t4 = -qJD(6) * t17 - t232 * t35 + t235 * t33;
t3 = qJD(6) * t16 + t232 * t33 + t235 * t35;
t5 = [(t452 + (-t272 / 0.2e1 + t431 / 0.2e1) * qJD(2)) * qJD(2) + (-t128 * t335 - t233 * t62 + t236 * t65) * mrSges(5,3) + m(7) * (t1 * t17 + t10 * t3 + t124 * t36 + t16 * t2 + t4 * t9 + t76 * t91) + t3 * t55 + t4 * t56 + t16 * t18 + t17 * t19 + (Ifges(7,4) * t60 + Ifges(7,2) * t61) * t404 + (-m(3) * t338 - m(4) * t289 - m(6) * t269 - t110 * mrSges(7,1) - t109 * mrSges(7,2) + t448 * (-qJ(4) * t234 + t269) + t415 * t234 + (-m(7) * (t187 * t233 - t346) - (m(6) * qJ(5) + mrSges(6,3)) * t236 - t253 * t233 + t425 + t454) * t237) * g(2) + t286 * t225 + t46 * t434 + (t233 * (Ifges(4,1) * t236 + t370) + t236 * (-Ifges(3,2) * t233 + t375) + t258 + t241) * t295 + (t233 * t445 - t236 * t271 + t280 + t375) * t396 - (t226 * t73 + t131) * t335 / 0.2e1 + (-t236 * Ifges(4,3) + t279 + t370) * t159 / 0.2e1 + (Ifges(6,5) * t102 + Ifges(6,6) * t101 + t319 + (-Ifges(3,4) + Ifges(4,5)) * t159 + (Ifges(4,1) + t445) * t160) * t233 / 0.2e1 + t233 * t458 + (-t108 * mrSges(7,1) - t107 * mrSges(7,2) + t448 * (-qJ(4) * t237 + t216) + (-m(3) + t450) * t216 + t415 * t237 + (-m(4) * (t297 - t215) - m(5) * t297 + t211 + (m(6) * t230 + m(7) * t178 + t282) * t233 + (m(3) + t449) * pkin(1) + t460 * t236 - t425) * t234) * g(1) - t233 * t459 + t419 * t335 / 0.2e1 + (t438 * pkin(7) + t423 / 0.2e1 + t161 * mrSges(4,2) - t111 * mrSges(5,3) + Ifges(7,3) * t394 + Ifges(7,5) * t402 + Ifges(7,6) * t404 - t451) * t334 + t424 * t296 + (-t167 * t335 + t429) * mrSges(4,2) + (t418 - t440) * (-t192 + t261) + (-m(5) * t65 + m(6) * t57 - t116 + t54) * (-t210 + t388) - t436 * t318 - t233 * (Ifges(5,4) * t159 - Ifges(5,2) * t160) / 0.2e1 + m(4) * (pkin(7) * t429 + t121 * t135 + t162 * t71) + (-t159 * t388 + t160 * t389 + t430) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t430) + (Ifges(7,1) * t60 + Ifges(7,4) * t61) * t402 + t42 * t349 / 0.2e1 + t13 * (-mrSges(6,2) * t233 + mrSges(6,3) * t349) + t12 * (mrSges(6,1) * t233 + mrSges(6,3) * t347) + m(6) * (t12 * t63 + t13 * t64 + t37 * t40 + t38 * t41) + m(5) * (t111 * t122 + t137 * t46 + t171 * t62 + t95 * t96) + Ifges(2,3) * qJDD(1) - t262 * t325 + t137 * t294 + (Ifges(3,4) * t160 - Ifges(3,2) * t159) * t446 - t159 * t277 / 0.2e1 + t64 * t69 + t63 * t70 + t76 * (-mrSges(7,1) * t61 + mrSges(7,2) * t60) - t71 * t284 + t91 * t34 + t41 * t97 + t40 * t98 + t113 * t6 / 0.2e1 + t36 * (-mrSges(7,1) * t112 + mrSges(7,2) * t113) + t124 * t8 + (t1 * t112 + t10 * t61 - t113 * t2 - t60 * t9) * mrSges(7,3) - t121 * t155 + t95 * t158 - t57 * t281 * t236 + t275 * t447 + t119 * t388 + t117 * t389 + (Ifges(7,5) * t113 + Ifges(7,6) * t112 + Ifges(7,3) * t233) * t398 + (Ifges(6,5) * t233 - t236 * t278) * t399 + (Ifges(6,6) * t233 - t236 * t274) * t400 + t347 * t406 + t60 * t407 + t61 * t408 + (Ifges(7,4) * t113 + Ifges(7,2) * t112 + Ifges(7,6) * t233) * t409 + (Ifges(7,1) * t113 + Ifges(7,4) * t112 + Ifges(7,5) * t233) * t410 + t112 * t411 - pkin(1) * (mrSges(3,1) * t159 + mrSges(3,2) * t160) + t162 * (mrSges(4,1) * t159 - mrSges(4,3) * t160) + t122 * t164 + t171 * t118 + (Ifges(7,5) * t60 + Ifges(7,6) * t61) * t394 + (-mrSges(3,1) * t389 - mrSges(3,2) * t388 + t456 * t233 + t462 * t392 + (Ifges(3,6) - t461) * t446) * qJDD(2) + ((-Ifges(5,4) + Ifges(4,5)) * t160 + t464 * t159) * t392; t426 * t285 + t62 * mrSges(5,2) - t65 * mrSges(5,1) + t453 * mrSges(7,3) + (-t116 + t119) * qJ(3) + (Ifges(7,4) * t100 - Ifges(7,2) * t99) * t405 + (-Ifges(6,6) * t396 - Ifges(6,2) * t400 - t42 / 0.2e1) * t227 - t372 * t400 + (Ifges(7,5) * t100 - Ifges(7,6) * t99) * t395 + (-Ifges(6,1) * t399 - Ifges(6,5) * t396 + t303 * t73 + t406) * t226 + (t236 * t420 + t422) * t457 + (Ifges(7,1) * t100 - Ifges(7,4) * t99) * t403 + t111 * t312 + t128 * t313 + t167 * t315 - (Ifges(4,1) * t336 + t196 + t419) * t337 / 0.2e1 + (t237 * t422 + t342 * t420) * g(1) - (-Ifges(3,2) * t337 + t197 + t423) * t336 / 0.2e1 + (Ifges(7,5) * t403 + Ifges(7,6) * t405 + Ifges(7,3) * t395 + t451) * t336 + (mrSges(7,1) * t356 - mrSges(7,2) * t341) * t76 + t440 * t154 + t436 * t198 + (m(4) * t167 + t416 - t437) * qJD(3) - t438 * t199 + t441 * t56 + t442 * t55 + (t1 * t85 + t10 * t442 - t103 * t76 + t178 * t36 + t2 * t84 + t441 * t9) * m(7) + t456 * t160 - t270 * mrSges(6,3) + t427 * qJD(5) + t428 * t218 + t431 * t296 + (-t452 + (t424 / 0.2e1 - t258 / 0.2e1 - t241 / 0.2e1 + t262) * qJD(1)) * qJD(1) - t161 * t314 + (t106 * t154 - t37 * t50 - t38 * t51 + t230 * t57 + t270 * t218 + (t226 * t37 - t227 * t38) * qJD(5)) * m(6) + (-t65 * qJ(3) - t111 * t157 - t115 * t96 - t128 * t154 + t238 * t62) * m(5) + (-pkin(2) * t114 + qJ(3) * t105 - t135 * t156) * m(4) + t238 * t118 - t371 * t399 + t272 * t295 + t461 * t159 + t57 * t282 + t84 * t18 + t85 * t19 - t51 * t97 - t50 * t98 - t100 * t31 / 0.2e1 - t103 * t34 + t105 * mrSges(4,3) - t114 * mrSges(4,1) - pkin(2) * t117 - t126 * t30 / 0.2e1 - t140 * mrSges(3,2) - t141 * mrSges(3,1) - t149 * t6 / 0.2e1 + t36 * (-mrSges(7,1) * t148 - mrSges(7,2) * t149) + t156 * t155 - t115 * t158 + t131 * t303 + (Ifges(7,5) * t125 - Ifges(7,6) * t126) * t394 + (-Ifges(7,5) * t149 + Ifges(7,6) * t148) * t398 + (Ifges(7,1) * t125 - Ifges(7,4) * t126) * t402 + (Ifges(7,4) * t125 - Ifges(7,2) * t126) * t404 + t125 * t407 + t99 * t408 + (-Ifges(7,4) * t149 + Ifges(7,2) * t148) * t409 + (-Ifges(7,1) * t149 + Ifges(7,4) * t148) * t410 + t148 * t411 - t157 * t164 + t178 * t8 + (Ifges(3,3) + Ifges(5,3) + Ifges(4,2)) * qJDD(2) + t230 * t54 + (-m(7) * (t311 - t346) - m(4) * t339 - m(5) * t311 - m(6) * t268 - t236 * mrSges(6,3) + (-t251 - t282) * t233 + t454 + t455) * g(3); -t148 * t19 - t149 * t18 - t341 * t56 + t356 * t55 + t437 * qJD(2) + t316 * t385 + ((-t226 * t97 - t227 * t98 - t155 - t158) * qJD(1) - t426 * t316) * t233 + t117 + t118 + t428 + (-qJD(2) * t76 - t453) * m(7) + (t270 - qJD(2) * t106 - (t226 * t38 + t227 * t37) * t337) * m(6) + (qJD(2) * t128 - t337 * t96 + t62) * m(5) + (-qJD(2) * t167 + t135 * t337 + t114) * m(4); -t148 * t18 + t149 * t19 + t226 * t69 + t227 * t70 + t356 * t56 + t341 * t55 + m(5) * t46 + m(6) * (t12 * t227 + t13 * t226) + (-m(6) * (-t348 * t38 + t350 * t37) + (m(5) * t111 + t164 - t427) * t233 + (t320 + t416) * t236) * qJD(1) + t294 + (g(1) * t234 - g(2) * t237) * t326 + (t1 * t149 + t10 * t341 - t148 * t2 + t356 * t9) * m(7); -t449 * t233 * g(3) - t142 * t97 + t143 * t98 - t291 * t55 + t80 * t56 + t54 + t8 + (-t10 * t291 + t80 * t9 + t36 - t439) * m(7) + (-t142 * t38 + t143 * t37 - t439 + t57) * m(6); -t459 + t458 - t76 * (mrSges(7,1) * t80 + mrSges(7,2) * t291) + (Ifges(7,1) * t291 - t391) * t403 + t30 * t402 + (Ifges(7,5) * t291 - Ifges(7,6) * t80) * t395 - t9 * t55 + t10 * t56 - g(1) * (mrSges(7,1) * t109 - mrSges(7,2) * t110) - g(2) * (-mrSges(7,1) * t107 + mrSges(7,2) * t108) - (mrSges(7,1) * t203 + mrSges(7,2) * t204) * t385 + (t10 * t80 + t291 * t9) * mrSges(7,3) + t319 + (-Ifges(7,2) * t80 + t31 + t75) * t405;];
tau  = t5;
