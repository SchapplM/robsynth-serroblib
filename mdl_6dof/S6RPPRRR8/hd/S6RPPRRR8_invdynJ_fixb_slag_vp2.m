% Calculate vector of inverse dynamics joint torques for
% S6RPPRRR8
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRR8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR8_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR8_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR8_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR8_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR8_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:35:25
% EndTime: 2019-03-09 02:35:51
% DurationCPUTime: 16.98s
% Computational Cost: add. (10445->648), mult. (21799->838), div. (0->0), fcn. (15554->14), ass. (0->299)
t239 = -pkin(9) - pkin(8);
t440 = -m(6) * pkin(8) + m(7) * t239 + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t223 = pkin(10) + qJ(4);
t213 = sin(t223);
t214 = cos(t223);
t439 = mrSges(5,1) * t213 + t440 * t214;
t228 = sin(pkin(10));
t229 = cos(pkin(10));
t305 = t228 ^ 2 + t229 ^ 2;
t278 = t305 * mrSges(4,3);
t233 = sin(qJ(5));
t290 = qJD(5) * t239;
t357 = cos(qJ(4));
t282 = qJD(1) * t357;
t234 = sin(qJ(4));
t304 = qJD(1) * t234;
t167 = -t228 * t282 - t229 * t304;
t324 = t167 * t233;
t231 = -pkin(1) - qJ(3);
t190 = qJD(1) * t231 + qJD(2);
t276 = -pkin(7) * qJD(1) + t190;
t159 = t276 * t228;
t160 = t276 * t229;
t100 = -t234 * t159 + t160 * t357;
t168 = -t228 * t304 + t229 * t282;
t125 = pkin(4) * t168 - pkin(8) * t167;
t237 = cos(qJ(5));
t63 = t237 * t100 + t233 * t125;
t438 = pkin(9) * t324 + t233 * t290 - t63;
t323 = t167 * t237;
t62 = -t100 * t233 + t237 * t125;
t437 = -pkin(5) * t168 + pkin(9) * t323 + t237 * t290 - t62;
t302 = qJD(5) * t233;
t436 = t302 - t324;
t235 = sin(qJ(1));
t238 = cos(qJ(1));
t435 = g(1) * t235 - g(2) * t238;
t434 = -qJD(1) * qJD(3) + qJDD(1) * t231;
t211 = qJD(1) * qJ(2) + qJD(3);
t217 = t228 * pkin(3);
t185 = qJD(1) * t217 + t211;
t407 = Ifges(5,5) * qJD(4);
t433 = -t407 / 0.2e1 - t185 * mrSges(5,2);
t236 = cos(qJ(6));
t232 = sin(qJ(6));
t140 = qJD(4) * t237 - t168 * t233;
t101 = t159 * t357 + t234 * t160;
t94 = qJD(4) * pkin(8) + t101;
t95 = -pkin(4) * t167 - pkin(8) * t168 + t185;
t55 = t233 * t95 + t237 * t94;
t40 = pkin(9) * t140 + t55;
t332 = t232 * t40;
t162 = qJD(5) - t167;
t141 = qJD(4) * t233 + t168 * t237;
t54 = -t233 * t94 + t237 * t95;
t39 = -pkin(9) * t141 + t54;
t34 = pkin(5) * t162 + t39;
t13 = t236 * t34 - t332;
t330 = t236 * t40;
t14 = t232 * t34 + t330;
t406 = Ifges(5,6) * qJD(4);
t432 = -t185 * mrSges(5,1) - t54 * mrSges(6,1) - t13 * mrSges(7,1) + t55 * mrSges(6,2) + t14 * mrSges(7,2) + t406 / 0.2e1;
t210 = pkin(5) * t237 + pkin(4);
t431 = -m(6) * pkin(4) - m(7) * t210;
t248 = -t234 * t228 + t357 * t229;
t289 = t357 * t228;
t249 = -t234 * t229 - t289;
t126 = qJD(1) * qJD(4) * t249 + qJDD(1) * t248;
t75 = qJD(5) * t140 + qJDD(4) * t233 + t126 * t237;
t376 = t75 / 0.2e1;
t76 = -qJD(5) * t141 + qJDD(4) * t237 - t126 * t233;
t375 = t76 / 0.2e1;
t298 = qJDD(1) * t229;
t127 = -qJD(4) * t168 - qJDD(1) * t289 - t234 * t298;
t123 = qJDD(5) - t127;
t369 = t123 / 0.2e1;
t161 = Ifges(5,4) * t167;
t429 = t140 * Ifges(6,6);
t428 = t162 * Ifges(6,3);
t427 = t167 * Ifges(5,2);
t383 = m(7) * pkin(5);
t270 = t236 * t140 - t141 * t232;
t25 = qJD(6) * t270 + t232 * t76 + t236 * t75;
t382 = t25 / 0.2e1;
t84 = t140 * t232 + t141 * t236;
t26 = -qJD(6) * t84 - t232 * t75 + t236 * t76;
t381 = t26 / 0.2e1;
t423 = Ifges(6,1) * t376 + Ifges(6,4) * t375 + Ifges(6,5) * t369;
t119 = qJDD(6) + t123;
t370 = t119 / 0.2e1;
t158 = qJD(6) + t162;
t418 = t141 * Ifges(6,5) + t84 * Ifges(7,5) + Ifges(7,6) * t270 + t158 * Ifges(7,3) + t428 + t429;
t187 = t239 * t233;
t188 = t239 * t237;
t138 = t187 * t236 + t188 * t232;
t417 = qJD(6) * t138 + t232 * t437 + t236 * t438;
t139 = t187 * t232 - t188 * t236;
t416 = -qJD(6) * t139 - t232 * t438 + t236 * t437;
t415 = -mrSges(6,1) - t383;
t255 = t232 * t233 - t236 * t237;
t396 = qJD(5) + qJD(6);
t131 = t396 * t255;
t99 = t255 * t167;
t414 = -t131 + t99;
t179 = t232 * t237 + t233 * t236;
t132 = t396 * t179;
t98 = t179 * t167;
t413 = -t132 + t98;
t346 = mrSges(5,3) * t168;
t412 = qJD(4) * mrSges(5,1) + mrSges(6,1) * t140 - mrSges(6,2) * t141 - t346;
t281 = qJD(4) * t357;
t303 = qJD(4) * t234;
t170 = -t228 * t303 + t229 * t281;
t404 = t255 * t249;
t411 = t255 * qJD(1) - t179 * t170 - t396 * t404;
t410 = -t179 * qJD(1) + t132 * t249 - t170 * t255;
t35 = -mrSges(6,1) * t76 + mrSges(6,2) * t75;
t409 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t126 + t35;
t408 = pkin(5) * t436 - t101;
t405 = t179 * t248;
t203 = qJ(2) + t217;
t124 = -pkin(4) * t249 - pkin(8) * t248 + t203;
t348 = -pkin(7) + t231;
t183 = t348 * t228;
t184 = t348 * t229;
t130 = t183 * t357 + t234 * t184;
t128 = t237 * t130;
t68 = t233 * t124 + t128;
t403 = -t234 * t183 + t357 * t184;
t225 = qJD(1) * qJD(2);
t193 = -qJDD(1) * qJ(2) - t225;
t227 = qJ(5) + qJ(6);
t215 = sin(t227);
t216 = cos(t227);
t265 = -mrSges(6,1) * t237 + mrSges(6,2) * t233;
t402 = mrSges(7,1) * t216 - mrSges(7,2) * t215 - t265 - t431;
t47 = mrSges(6,1) * t123 - mrSges(6,3) * t75;
t48 = -mrSges(6,2) * t123 + mrSges(6,3) * t76;
t400 = -t233 * t47 + t237 * t48;
t301 = qJD(5) * t237;
t180 = qJDD(2) + t434;
t271 = -pkin(7) * qJDD(1) + t180;
t149 = t271 * t228;
t150 = t271 * t229;
t56 = t357 * t149 + t234 * t150 - t159 * t303 + t160 * t281;
t52 = qJDD(4) * pkin(8) + t56;
t186 = qJDD(3) - t193;
t299 = qJDD(1) * t228;
t171 = pkin(3) * t299 + t186;
t64 = -pkin(4) * t127 - pkin(8) * t126 + t171;
t11 = t233 * t64 + t237 * t52 + t95 * t301 - t302 * t94;
t12 = -qJD(5) * t55 - t233 * t52 + t237 * t64;
t399 = t11 * t237 - t12 * t233;
t398 = -m(6) - m(5) - m(7);
t93 = -qJD(4) * pkin(4) - t100;
t69 = -t140 * pkin(5) + t93;
t394 = -t69 * mrSges(7,1) + t14 * mrSges(7,3);
t393 = t69 * mrSges(7,2) - t13 * mrSges(7,3);
t392 = mrSges(3,2) - mrSges(2,1) - mrSges(4,3) - mrSges(5,3);
t169 = -t228 * t281 - t229 * t303;
t57 = -t234 * t149 + t150 * t357 - t159 * t281 - t160 * t303;
t391 = -t100 * t169 - t101 * t170 - t248 * t57 + t249 * t56;
t390 = -m(6) * t93 + t412;
t10 = pkin(9) * t76 + t11;
t9 = pkin(5) * t123 - pkin(9) * t75 + t12;
t2 = qJD(6) * t13 + t10 * t236 + t232 * t9;
t3 = -qJD(6) * t14 - t10 * t232 + t236 * t9;
t389 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t388 = t12 * mrSges(6,1) - t11 * mrSges(6,2);
t268 = mrSges(4,1) * t228 + mrSges(4,2) * t229;
t387 = t213 * t431 + mrSges(2,2) - mrSges(3,3) - t268 - t439;
t386 = m(4) * t211 + m(5) * t185 - mrSges(5,1) * t167 + mrSges(5,2) * t168 + t268 * qJD(1);
t385 = Ifges(7,4) * t382 + Ifges(7,2) * t381 + Ifges(7,6) * t370;
t384 = Ifges(7,1) * t382 + Ifges(7,4) * t381 + Ifges(7,5) * t370;
t356 = Ifges(7,4) * t84;
t37 = Ifges(7,2) * t270 + Ifges(7,6) * t158 + t356;
t380 = -t37 / 0.2e1;
t379 = t37 / 0.2e1;
t80 = Ifges(7,4) * t270;
t38 = Ifges(7,1) * t84 + Ifges(7,5) * t158 + t80;
t378 = -t38 / 0.2e1;
t377 = t38 / 0.2e1;
t374 = -t270 / 0.2e1;
t373 = t270 / 0.2e1;
t372 = -t84 / 0.2e1;
t371 = t84 / 0.2e1;
t368 = -t140 / 0.2e1;
t367 = -t141 / 0.2e1;
t366 = t141 / 0.2e1;
t365 = -t158 / 0.2e1;
t364 = t158 / 0.2e1;
t363 = -t162 / 0.2e1;
t362 = -t167 / 0.2e1;
t360 = t168 / 0.2e1;
t355 = pkin(5) * t141;
t352 = g(3) * t214;
t347 = mrSges(5,3) * t167;
t345 = mrSges(6,3) * t140;
t344 = mrSges(6,3) * t141;
t343 = Ifges(5,4) * t168;
t342 = Ifges(6,4) * t141;
t341 = Ifges(6,4) * t233;
t340 = Ifges(6,4) * t237;
t53 = -qJDD(4) * pkin(4) - t57;
t336 = t248 * t53;
t328 = pkin(1) * qJDD(1);
t322 = t169 * t237;
t321 = t248 * t233;
t320 = t248 * t237;
t319 = t215 * t235;
t318 = t215 * t238;
t317 = t216 * t235;
t316 = t216 * t238;
t315 = t233 * t235;
t314 = t233 * t238;
t311 = t235 * t237;
t310 = t237 * t238;
t154 = -t213 * t319 + t316;
t155 = t213 * t317 + t318;
t309 = t154 * mrSges(7,1) - t155 * mrSges(7,2);
t156 = t213 * t318 + t317;
t157 = t213 * t316 - t319;
t308 = t156 * mrSges(7,1) + t157 * mrSges(7,2);
t307 = mrSges(4,1) * t299 + mrSges(4,2) * t298;
t306 = t238 * pkin(1) + t235 * qJ(2);
t46 = -mrSges(7,1) * t270 + mrSges(7,2) * t84;
t297 = -t46 + t412;
t296 = Ifges(7,5) * t25 + Ifges(7,6) * t26 + Ifges(7,3) * t119;
t295 = Ifges(6,5) * t75 + Ifges(6,6) * t76 + Ifges(6,3) * t123;
t293 = -m(4) + t398;
t71 = t140 * Ifges(6,2) + t162 * Ifges(6,6) + t342;
t291 = -t233 * t71 / 0.2e1;
t280 = t301 / 0.2e1;
t219 = t238 * qJ(2);
t279 = -pkin(1) * t235 + t219;
t277 = -t127 * mrSges(5,1) + t126 * mrSges(5,2);
t120 = pkin(4) * t170 - pkin(8) * t169 + qJD(2);
t89 = t249 * qJD(3) + qJD(4) * t403;
t275 = t237 * t120 - t233 * t89;
t274 = t180 * t305;
t273 = t305 * t190;
t67 = t237 * t124 - t130 * t233;
t264 = mrSges(6,1) * t233 + mrSges(6,2) * t237;
t263 = -mrSges(7,1) * t215 - mrSges(7,2) * t216;
t262 = Ifges(6,1) * t237 - t341;
t261 = -Ifges(6,2) * t233 + t340;
t260 = Ifges(6,5) * t237 - Ifges(6,6) * t233;
t50 = -pkin(5) * t249 - pkin(9) * t320 + t67;
t58 = -pkin(9) * t321 + t68;
t27 = -t232 * t58 + t236 * t50;
t28 = t232 * t50 + t236 * t58;
t259 = t233 * t55 + t237 * t54;
t258 = -t233 * t54 + t237 * t55;
t96 = -mrSges(6,2) * t162 + t345;
t97 = mrSges(6,1) * t162 - t344;
t257 = -t233 * t97 + t237 * t96;
t256 = -t233 * t96 - t237 * t97;
t252 = t296 + t389;
t152 = -qJD(4) * mrSges(5,2) + t347;
t251 = -t152 - t257;
t165 = t213 * t314 + t311;
t163 = -t213 * t315 + t310;
t250 = t93 * t264;
t247 = -t169 * t233 - t248 * t301;
t246 = -t248 * t302 + t322;
t32 = t233 * t120 + t124 * t301 - t130 * t302 + t237 * t89;
t241 = -qJD(5) * t259 + t399;
t90 = qJD(3) * t248 + qJD(4) * t130;
t240 = qJD(1) ^ 2;
t230 = -pkin(7) - qJ(3);
t212 = qJDD(2) - t328;
t166 = t213 * t310 - t315;
t164 = t213 * t311 + t314;
t134 = Ifges(6,4) * t140;
t117 = t255 * t248;
t114 = t179 * t249;
t108 = t168 * Ifges(5,1) + t161 + t407;
t107 = t343 + t406 + t427;
t106 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t127;
t91 = pkin(5) * t321 - t403;
t72 = t141 * Ifges(6,1) + t162 * Ifges(6,5) + t134;
t66 = mrSges(7,1) * t158 - mrSges(7,3) * t84;
t65 = -mrSges(7,2) * t158 + mrSges(7,3) * t270;
t61 = -pkin(5) * t247 + t90;
t44 = t131 * t248 - t169 * t179;
t42 = -t255 * t169 - t396 * t405;
t33 = -qJD(5) * t68 + t275;
t31 = -t76 * pkin(5) + t53;
t29 = t75 * Ifges(6,4) + t76 * Ifges(6,2) + t123 * Ifges(6,6);
t22 = pkin(9) * t247 + t32;
t21 = -pkin(9) * t322 + pkin(5) * t170 + (-t128 + (pkin(9) * t248 - t124) * t233) * qJD(5) + t275;
t18 = -mrSges(7,2) * t119 + mrSges(7,3) * t26;
t17 = mrSges(7,1) * t119 - mrSges(7,3) * t25;
t16 = t236 * t39 - t332;
t15 = -t232 * t39 - t330;
t8 = -mrSges(7,1) * t26 + mrSges(7,2) * t25;
t5 = -qJD(6) * t28 + t21 * t236 - t22 * t232;
t4 = qJD(6) * t27 + t21 * t232 + t22 * t236;
t1 = [t320 * t423 + (-t180 - t434) * t278 - (-t171 * mrSges(5,2) - Ifges(5,1) * t126 - Ifges(5,4) * t127 - Ifges(5,5) * qJDD(4) - t260 * t369 - t261 * t375 - t262 * t376 + t280 * t71) * t248 - (t296 + t295) * t249 / 0.2e1 - (mrSges(5,1) * t171 - Ifges(5,4) * t126 + Ifges(6,5) * t376 + Ifges(7,5) * t382 - Ifges(5,2) * t127 - Ifges(5,6) * qJDD(4) + Ifges(6,6) * t375 + Ifges(7,6) * t381 + Ifges(6,3) * t369 + Ifges(7,3) * t370 + t388 + t389) * t249 + (Ifges(7,1) * t42 + Ifges(7,4) * t44) * t371 + (Ifges(6,1) * t246 + Ifges(6,4) * t247) * t366 + t162 * (Ifges(6,5) * t246 + Ifges(6,6) * t247) / 0.2e1 + (t212 - t328) * mrSges(3,2) + (Ifges(2,3) + Ifges(3,1)) * qJDD(1) + (t418 / 0.2e1 + t428 / 0.2e1 + t429 / 0.2e1 - t427 / 0.2e1 - t107 / 0.2e1 - Ifges(5,4) * t360 + Ifges(7,3) * t364 + Ifges(6,5) * t366 + Ifges(7,5) * t371 + Ifges(7,6) * t373 - t432) * t170 + (t161 / 0.2e1 + t108 / 0.2e1 + Ifges(5,1) * t360 + t291 - t433) * t169 + (Ifges(7,5) * t42 + Ifges(7,6) * t44) * t364 + (-t314 * t383 - t164 * mrSges(6,1) - t155 * mrSges(7,1) - t163 * mrSges(6,2) - t154 * mrSges(7,2) + (-m(4) - m(3)) * t306 + t398 * (t235 * t217 - t230 * t238 + t306) + (-m(4) * qJ(3) + t392) * t238 + t387 * t235) * g(2) + (t315 * t383 - m(3) * t279 - m(4) * t219 - t166 * mrSges(6,1) - t157 * mrSges(7,1) + t165 * mrSges(6,2) + t156 * mrSges(7,2) + t398 * (t238 * t217 + t235 * t230 + t279) + (-m(4) * t231 - t392) * t235 + t387 * t238) * g(1) + m(7) * (t13 * t5 + t14 * t4 + t2 * t28 + t27 * t3 + t31 * t91 + t61 * t69) - (qJD(5) * t72 + t29) * t321 / 0.2e1 + (-Ifges(7,1) * t117 - Ifges(7,4) * t405) * t382 + (t117 * t3 - t13 * t42 + t14 * t44 - t2 * t405) * mrSges(7,3) + (-Ifges(7,5) * t117 - Ifges(7,6) * t405) * t370 + (-Ifges(7,4) * t117 - Ifges(7,2) * t405) * t381 + t31 * (mrSges(7,1) * t405 - mrSges(7,2) * t117) - t405 * t385 - (-m(5) * t57 + m(6) * t53 + t409) * t403 + t140 * (Ifges(6,4) * t246 + Ifges(6,2) * t247) / 0.2e1 + (Ifges(4,1) * t229 - Ifges(4,4) * t228) * t298 + m(5) * (t101 * t89 + t130 * t56 + t171 * t203) + m(6) * (t11 * t68 + t12 * t67 + t32 * t55 + t33 * t54) + t391 * mrSges(5,3) + t93 * (-mrSges(6,1) * t247 + mrSges(6,2) * t246) + m(4) * (qJ(2) * t186 - qJD(3) * t273 + t231 * t274) + (-t11 * t321 - t12 * t320 - t246 * t54 + t247 * t55) * mrSges(6,3) + t386 * qJD(2) + (Ifges(7,4) * t42 + Ifges(7,2) * t44) * t373 - 0.2e1 * t193 * mrSges(3,3) + t89 * t152 + t130 * t106 + t264 * t336 + t32 * t96 + t33 * t97 + t91 * t8 + t68 * t48 + t69 * (-mrSges(7,1) * t44 + mrSges(7,2) * t42) + t72 * t322 / 0.2e1 + t4 * t65 + t5 * t66 + t67 * t47 + t61 * t46 + t203 * t277 + t42 * t377 + t44 * t379 - t117 * t384 + (-m(5) * t100 - t390) * t90 - (Ifges(4,4) * t229 - Ifges(4,2) * t228) * t299 + t27 * t17 + t28 * t18 + m(3) * (-pkin(1) * t212 + (-t193 + t225) * qJ(2)) + qJ(2) * t307 + t186 * t268; t114 * t17 + t404 * t18 + t411 * t66 + t410 * t65 + (-m(3) * qJ(2) - mrSges(3,3)) * t240 - (t8 + t409) * t248 - t251 * t170 + t297 * t169 + (mrSges(3,2) - t278) * qJDD(1) - (qJD(5) * t256 + t106 + t400) * t249 + (t256 - t386) * qJD(1) - m(5) * t391 + m(3) * t212 + m(4) * t274 - t435 * (m(3) - t293) + (t114 * t3 + t13 * t411 + t14 * t410 - t169 * t69 + t2 * t404 - t248 * t31) * m(7) + (-qJD(1) * t259 - t169 * t93 + t170 * t258 - t241 * t249 - t336) * m(6); t413 * t66 + t307 + t414 * t65 - t240 * t278 + t251 * t167 + t257 * qJD(5) + t237 * t47 + t233 * t48 - t255 * t17 + t179 * t18 + t277 + t297 * t168 + (g(1) * t238 + g(2) * t235) * t293 + (t13 * t413 + t14 * t414 - t168 * t69 + t179 * t2 - t255 * t3) * m(7) + (t11 * t233 + t12 * t237 + t162 * t258 - t168 * t93) * m(6) + (t100 * t168 - t101 * t167 + t171) * m(5) + (qJD(1) * t273 + t186) * m(4); t233 * t423 + (-t13 * t99 + t14 * t98 - t179 * t3 - t2 * t255) * mrSges(7,3) + t31 * (mrSges(7,1) * t255 + mrSges(7,2) * t179) + (Ifges(7,5) * t179 - Ifges(7,6) * t255) * t370 + (Ifges(7,4) * t179 - Ifges(7,2) * t255) * t381 + (Ifges(7,1) * t179 - Ifges(7,4) * t255) * t382 - t255 * t385 + (-t436 * t55 + (-t301 + t323) * t54 + t399) * mrSges(6,3) + t416 * t66 + (t13 * t416 + t138 * t3 + t139 * t2 + t14 * t417 - t210 * t31 + t408 * t69) * m(7) + t417 * t65 - (Ifges(5,1) * t167 - t343 + t418) * t168 / 0.2e1 + t408 * t46 + (m(6) * t241 - t301 * t97 - t302 * t96 + t400) * pkin(8) - (Ifges(7,1) * t371 + Ifges(7,4) * t373 + Ifges(7,5) * t364 + t377 + t393) * t131 - (Ifges(7,4) * t371 + Ifges(7,2) * t373 + Ifges(7,6) * t364 + t379 + t394) * t132 + (-t323 / 0.2e1 + t280) * t72 + (-pkin(4) * t53 - t54 * t62 - t55 * t63) * m(6) + (Ifges(6,5) * t367 + Ifges(7,5) * t372 - Ifges(5,2) * t362 + Ifges(6,6) * t368 + Ifges(7,6) * t374 + Ifges(6,3) * t363 + Ifges(7,3) * t365 + t432) * t168 + (t260 * t363 + t261 * t368 + t262 * t367 - t250 + t433) * t167 + (t213 * t402 + t439) * g(3) + (t140 * t261 + t141 * t262 + t162 * t260) * qJD(5) / 0.2e1 + (t250 + t291) * qJD(5) + (t161 + t108) * t362 + (-Ifges(7,5) * t99 - Ifges(7,6) * t98) * t365 + (-Ifges(7,4) * t99 - Ifges(7,2) * t98) * t374 + (-Ifges(7,1) * t99 - Ifges(7,4) * t98) * t372 - t69 * (mrSges(7,1) * t98 - mrSges(7,2) * t99) + (t347 - t152) * t100 + t53 * t265 + t237 * t29 / 0.2e1 - t210 * t8 + Ifges(5,3) * qJDD(4) + t138 * t17 + t139 * t18 + Ifges(5,5) * t126 + Ifges(5,6) * t127 - t63 * t96 - t62 * t97 + t71 * t324 / 0.2e1 + t57 * mrSges(5,1) - t56 * mrSges(5,2) + t107 * t360 + (Ifges(6,5) * t233 + Ifges(6,6) * t237) * t369 + (Ifges(6,2) * t237 + t341) * t375 + (Ifges(6,1) * t233 + t340) * t376 - t99 * t378 - t98 * t380 + t179 * t384 + (t346 + t390) * t101 + t435 * ((-mrSges(5,1) - t402) * t214 + t440 * t213) - pkin(4) * t35; t388 + (t344 + t97) * t55 + (t345 - t96) * t54 + (t233 * t383 - t263 + t264) * t352 - (Ifges(7,4) * t372 + Ifges(7,2) * t374 + Ifges(7,6) * t365 + t380 - t394) * t84 + (Ifges(7,1) * t372 + Ifges(7,4) * t374 + Ifges(7,5) * t365 + t378 - t393) * t270 + (-Ifges(6,2) * t141 + t134 + t72) * t368 + t295 + t252 + (mrSges(6,2) * t164 + t163 * t415 - t309) * g(1) + (-mrSges(6,2) * t166 + t165 * t415 - t308) * g(2) - t46 * t355 - m(7) * (t13 * t15 + t14 * t16 + t355 * t69) - t93 * (mrSges(6,1) * t141 + mrSges(6,2) * t140) - t16 * t65 - t15 * t66 + (Ifges(6,5) * t140 - Ifges(6,6) * t141) * t363 + t71 * t366 + (Ifges(6,1) * t140 - t342) * t367 + (t2 * t232 + t236 * t3 + (-t13 * t232 + t14 * t236) * qJD(6)) * t383 + ((-t232 * t66 + t236 * t65) * qJD(6) + t236 * t17 + t232 * t18) * pkin(5); -t69 * (mrSges(7,1) * t84 + mrSges(7,2) * t270) + (Ifges(7,1) * t270 - t356) * t372 + t37 * t371 + (Ifges(7,5) * t270 - Ifges(7,6) * t84) * t365 - t13 * t65 + t14 * t66 - g(1) * t309 - g(2) * t308 - t263 * t352 + (t13 * t270 + t14 * t84) * mrSges(7,3) + t252 + (-Ifges(7,2) * t84 + t38 + t80) * t374;];
tau  = t1;
