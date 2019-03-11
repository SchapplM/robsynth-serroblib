% Calculate vector of inverse dynamics joint torques for
% S6RPRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRP5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:14:33
% EndTime: 2019-03-09 03:15:10
% DurationCPUTime: 24.83s
% Computational Cost: add. (10650->653), mult. (25729->803), div. (0->0), fcn. (19765->14), ass. (0->283)
t398 = -Ifges(6,4) + Ifges(7,5);
t443 = t398 + Ifges(7,5);
t240 = sin(pkin(9));
t242 = cos(pkin(9));
t246 = sin(qJ(3));
t348 = cos(qJ(3));
t206 = t240 * t348 + t246 * t242;
t195 = t206 * qJD(1);
t239 = sin(pkin(10));
t241 = cos(pkin(10));
t169 = qJD(3) * t239 + t195 * t241;
t170 = qJD(3) * t241 - t195 * t239;
t245 = sin(qJ(5));
t248 = cos(qJ(5));
t111 = t245 * t169 - t170 * t248;
t291 = t348 * t242;
t259 = -t246 * t240 + t291;
t198 = t259 * qJD(3);
t157 = qJD(1) * t198 + qJDD(1) * t206;
t136 = qJDD(3) * t241 - t157 * t239;
t137 = qJDD(3) * t239 + t157 * t241;
t47 = -qJD(5) * t111 + t245 * t136 + t248 * t137;
t376 = t47 / 0.2e1;
t417 = t248 * t169 + t170 * t245;
t48 = qJD(5) * t417 - t248 * t136 + t245 * t137;
t374 = t48 / 0.2e1;
t199 = t206 * qJD(3);
t296 = qJDD(1) * t240;
t158 = qJD(1) * t199 - qJDD(1) * t291 + t246 * t296;
t154 = qJDD(5) + t158;
t361 = t154 / 0.2e1;
t401 = mrSges(6,1) + mrSges(7,1);
t400 = mrSges(6,2) - mrSges(7,3);
t439 = -mrSges(6,3) - mrSges(7,2);
t399 = Ifges(6,1) + Ifges(7,1);
t397 = Ifges(7,4) + Ifges(6,5);
t396 = -Ifges(6,6) + Ifges(7,6);
t426 = Ifges(6,3) + Ifges(7,2);
t297 = qJD(1) * qJD(2);
t221 = qJ(2) * qJDD(1) + t297;
t342 = pkin(7) + qJ(2);
t216 = t342 * t242;
t208 = qJD(1) * t216;
t192 = t246 * t208;
t214 = t342 * t240;
t207 = qJD(1) * t214;
t159 = -t207 * t348 - t192;
t149 = -qJD(3) * pkin(3) + qJD(4) - t159;
t108 = -pkin(4) * t170 + t149;
t39 = t111 * pkin(5) - qJ(6) * t417 + t108;
t109 = Ifges(7,5) * t417;
t194 = t259 * qJD(1);
t190 = qJD(5) - t194;
t53 = t190 * Ifges(7,6) + t111 * Ifges(7,3) + t109;
t328 = Ifges(6,4) * t417;
t56 = -t111 * Ifges(6,2) + t190 * Ifges(6,6) + t328;
t442 = t39 * mrSges(7,1) + t53 / 0.2e1 - t56 / 0.2e1;
t358 = t190 / 0.2e1;
t366 = t417 / 0.2e1;
t368 = t111 / 0.2e1;
t369 = -t111 / 0.2e1;
t110 = Ifges(6,4) * t111;
t327 = Ifges(7,5) * t111;
t393 = t190 * t397 + t399 * t417 - t110 + t327;
t421 = t39 * mrSges(7,3) - t393 / 0.2e1;
t441 = -Ifges(6,4) * t369 - Ifges(7,5) * t368 - t397 * t358 - t399 * t366 + t421;
t229 = pkin(2) * t242 + pkin(1);
t210 = -qJDD(1) * t229 + qJDD(2);
t79 = pkin(3) * t158 - qJ(4) * t157 - qJD(4) * t195 + t210;
t281 = pkin(7) * qJDD(1) + t221;
t183 = t281 * t240;
t184 = t281 * t242;
t284 = qJD(3) * t348;
t292 = -t246 * t183 + t348 * t184 - t207 * t284;
t90 = qJDD(3) * qJ(4) + (qJD(4) - t192) * qJD(3) + t292;
t37 = -t239 * t90 + t241 * t79;
t20 = pkin(4) * t158 - pkin(8) * t137 + t37;
t211 = -qJD(1) * t229 + qJD(2);
t125 = -pkin(3) * t194 - qJ(4) * t195 + t211;
t160 = -t246 * t207 + t348 * t208;
t151 = qJD(3) * qJ(4) + t160;
t84 = t241 * t125 - t151 * t239;
t59 = -pkin(4) * t194 - pkin(8) * t169 + t84;
t85 = t239 * t125 + t241 * t151;
t70 = pkin(8) * t170 + t85;
t22 = t245 * t59 + t248 * t70;
t38 = t239 * t79 + t241 * t90;
t24 = pkin(8) * t136 + t38;
t4 = -qJD(5) * t22 + t20 * t248 - t24 * t245;
t2 = -pkin(5) * t154 + qJDD(6) - t4;
t375 = -t48 / 0.2e1;
t302 = qJD(3) * t246;
t99 = -t183 * t348 - t246 * t184 + t207 * t302 - t208 * t284;
t91 = -qJDD(3) * pkin(3) + qJDD(4) - t99;
t60 = -t136 * pkin(4) + t91;
t6 = t48 * pkin(5) - t47 * qJ(6) - qJD(6) * t417 + t60;
t440 = mrSges(6,2) * t60 + mrSges(7,2) * t2 - mrSges(6,3) * t4 - mrSges(7,3) * t6 + Ifges(6,4) * t375 + 0.2e1 * t361 * t397 + t443 * t374 + 0.2e1 * t376 * t399;
t406 = -m(6) - m(7);
t275 = -mrSges(3,1) * t242 + mrSges(3,2) * t240;
t438 = -m(3) * pkin(1) - mrSges(2,1) + t275;
t205 = t239 * t248 + t241 * t245;
t134 = t205 * t194;
t197 = t205 * qJD(5);
t385 = -t134 + t197;
t203 = t239 * t245 - t248 * t241;
t135 = t203 * t194;
t196 = t203 * qJD(5);
t384 = -t135 + t196;
t303 = t240 ^ 2 + t242 ^ 2;
t237 = pkin(10) + qJ(5);
t231 = sin(t237);
t233 = cos(t237);
t437 = -t231 * t400 + t233 * t401;
t334 = mrSges(4,3) * t195;
t383 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t170 - t169 * mrSges(5,2) - t334;
t436 = -m(5) * t149 + t383;
t329 = Ifges(5,4) * t241;
t267 = -Ifges(5,2) * t239 + t329;
t330 = Ifges(5,4) * t239;
t268 = Ifges(5,1) * t241 - t330;
t271 = mrSges(5,1) * t239 + mrSges(5,2) * t241;
t349 = t241 / 0.2e1;
t389 = Ifges(4,5) * qJD(3);
t404 = t170 / 0.2e1;
t405 = t169 / 0.2e1;
t435 = t211 * mrSges(4,2) + t267 * t404 + t268 * t405 + (Ifges(5,1) * t169 + Ifges(5,4) * t170 - Ifges(5,5) * t194) * t349 - t239 * (Ifges(5,4) * t169 + Ifges(5,2) * t170 - Ifges(5,6) * t194) / 0.2e1 - t159 * mrSges(4,3) + t149 * t271 + t389 / 0.2e1;
t21 = -t245 * t70 + t248 * t59;
t17 = -pkin(5) * t190 + qJD(6) - t21;
t19 = qJ(6) * t190 + t22;
t388 = Ifges(4,6) * qJD(3);
t434 = -t211 * mrSges(4,1) - t84 * mrSges(5,1) - t21 * mrSges(6,1) + t17 * mrSges(7,1) + t85 * mrSges(5,2) + t22 * mrSges(6,2) - t19 * mrSges(7,3) + t388 / 0.2e1;
t431 = mrSges(6,1) * t60 + mrSges(7,1) * t6 - t47 * Ifges(6,4) / 0.2e1 - t154 * Ifges(6,6) / 0.2e1 + 0.2e1 * Ifges(7,3) * t374 + (-t375 + t374) * Ifges(6,2) + t443 * t376 + (t396 + Ifges(7,6)) * t361;
t429 = -Ifges(6,2) * t369 + Ifges(7,3) * t368 + t358 * t396 + t366 * t398 + t442;
t423 = t154 * t426 + t396 * t48 + t397 * t47;
t247 = sin(qJ(1));
t249 = cos(qJ(1));
t422 = g(1) * t249 + g(2) * t247;
t391 = t170 * Ifges(5,6);
t392 = Ifges(5,5) * t169;
t420 = -Ifges(5,3) * t194 + t111 * t396 + t190 * t426 + t397 * t417 + t391 + t392;
t238 = pkin(9) + qJ(3);
t232 = sin(t238);
t234 = cos(t238);
t274 = t234 * mrSges(4,1) - t232 * mrSges(4,2);
t419 = t232 * t439 - t274;
t286 = m(3) * qJ(2) + mrSges(3,3);
t416 = -m(5) * t342 + mrSges(2,2) - mrSges(4,3) - t271 - t286;
t415 = m(7) * pkin(5) + t401;
t414 = m(7) * qJ(6) - t400;
t300 = qJD(5) * t248;
t301 = qJD(5) * t245;
t3 = t245 * t20 + t248 * t24 + t59 * t300 - t301 * t70;
t1 = qJ(6) * t154 + qJD(6) * t190 + t3;
t411 = -t4 * mrSges(6,1) + t2 * mrSges(7,1) + t3 * mrSges(6,2) - t1 * mrSges(7,3);
t364 = t136 / 0.2e1;
t363 = t137 / 0.2e1;
t360 = t158 / 0.2e1;
t33 = -mrSges(7,2) * t48 + mrSges(7,3) * t154;
t36 = -mrSges(6,2) * t154 - mrSges(6,3) * t48;
t395 = t33 + t36;
t34 = mrSges(6,1) * t154 - mrSges(6,3) * t47;
t35 = -t154 * mrSges(7,1) + t47 * mrSges(7,2);
t394 = -t34 + t35;
t337 = mrSges(7,2) * t111;
t94 = mrSges(7,3) * t190 - t337;
t333 = mrSges(6,3) * t111;
t95 = -mrSges(6,2) * t190 - t333;
t339 = t94 + t95;
t332 = mrSges(6,3) * t417;
t96 = mrSges(6,1) * t190 - t332;
t336 = mrSges(7,2) * t417;
t97 = -mrSges(7,1) * t190 + t336;
t338 = t97 - t96;
t317 = t194 * t239;
t118 = pkin(4) * t317 + t160;
t390 = pkin(5) * t385 + qJ(6) * t384 - qJD(6) * t205 - t118;
t272 = -t241 * mrSges(5,1) + t239 * mrSges(5,2);
t257 = m(5) * pkin(3) - t272;
t386 = t257 * t234;
t382 = -t348 * t214 - t246 * t216;
t131 = t194 * mrSges(5,2) + mrSges(5,3) * t170;
t132 = -mrSges(5,1) * t194 - mrSges(5,3) * t169;
t381 = t131 * t241 - t132 * t239;
t380 = -t239 * t37 + t241 * t38;
t156 = -pkin(3) * t259 - qJ(4) * t206 - t229;
t167 = -t246 * t214 + t216 * t348;
t105 = t241 * t156 - t167 * t239;
t312 = t206 * t241;
t78 = -pkin(4) * t259 - pkin(8) * t312 + t105;
t106 = t239 * t156 + t241 * t167;
t313 = t206 * t239;
t86 = -pkin(8) * t313 + t106;
t340 = t245 * t78 + t248 * t86;
t314 = t198 * t241;
t120 = pkin(3) * t199 - qJ(4) * t198 - qJD(4) * t206;
t129 = t259 * qJD(2) + qJD(3) * t382;
t76 = t241 * t120 - t129 * t239;
t51 = pkin(4) * t199 - pkin(8) * t314 + t76;
t315 = t198 * t239;
t77 = t239 * t120 + t241 * t129;
t65 = -pkin(8) * t315 + t77;
t9 = -qJD(5) * t340 - t245 * t65 + t248 * t51;
t371 = Ifges(5,1) * t363 + Ifges(5,4) * t364 + Ifges(5,5) * t360;
t367 = -t417 / 0.2e1;
t359 = -t190 / 0.2e1;
t357 = t194 / 0.2e1;
t356 = -t194 / 0.2e1;
t354 = t195 / 0.2e1;
t346 = pkin(4) * t239;
t343 = g(3) * t232;
t341 = pkin(8) + qJ(4);
t152 = pkin(3) * t195 - qJ(4) * t194;
t100 = t241 * t152 - t159 * t239;
t316 = t194 * t241;
t71 = pkin(4) * t195 - pkin(8) * t316 + t100;
t101 = t239 * t152 + t241 * t159;
t81 = -pkin(8) * t317 + t101;
t28 = t245 * t71 + t248 * t81;
t331 = Ifges(4,4) * t195;
t311 = t231 * t247;
t310 = t232 * t341;
t309 = t232 * t249;
t228 = pkin(4) * t241 + pkin(3);
t209 = t234 * t228;
t308 = t234 * t249;
t305 = t247 * t233;
t304 = t249 * t231;
t299 = m(5) - t406;
t295 = qJDD(1) * t242;
t285 = m(5) * qJ(4) + mrSges(5,3);
t15 = t48 * mrSges(6,1) + t47 * mrSges(6,2);
t14 = t48 * mrSges(7,1) - t47 * mrSges(7,3);
t283 = t158 * mrSges(4,1) + t157 * mrSges(4,2);
t83 = -t136 * mrSges(5,1) + t137 * mrSges(5,2);
t217 = t249 * t229;
t279 = t247 * t342 + t217;
t133 = pkin(4) * t313 - t382;
t276 = -mrSges(3,1) * t295 + mrSges(3,2) * t296;
t266 = Ifges(5,5) * t241 - Ifges(5,6) * t239;
t265 = pkin(5) * t233 + qJ(6) * t231;
t263 = t239 * t84 - t241 * t85;
t27 = -t245 * t81 + t248 * t71;
t31 = -t245 * t86 + t248 * t78;
t213 = t341 * t239;
t215 = t341 * t241;
t260 = -t248 * t213 - t215 * t245;
t166 = -t213 * t245 + t215 * t248;
t8 = t245 * t51 + t248 * t65 + t78 * t300 - t301 * t86;
t254 = t232 * t285 + t386;
t130 = qJD(2) * t206 + qJD(3) * t167;
t107 = pkin(4) * t315 + t130;
t230 = -qJDD(1) * pkin(1) + qJDD(2);
t186 = Ifges(4,4) * t194;
t175 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t194;
t174 = t233 * t308 + t311;
t173 = t234 * t304 - t305;
t172 = t234 * t305 - t304;
t171 = t233 * t249 + t234 * t311;
t155 = pkin(5) * t203 - qJ(6) * t205 - t228;
t146 = t195 * Ifges(4,1) + t186 + t389;
t145 = t194 * Ifges(4,2) + t331 + t388;
t142 = t203 * t206;
t141 = t205 * t206;
t128 = qJD(4) * t205 + qJD(5) * t166;
t127 = -qJD(4) * t203 + qJD(5) * t260;
t98 = -t208 * t302 + t292;
t93 = mrSges(5,1) * t158 - mrSges(5,3) * t137;
t92 = -mrSges(5,2) * t158 + mrSges(5,3) * t136;
t89 = t198 * t205 + t300 * t312 - t301 * t313;
t88 = -t197 * t206 - t198 * t203;
t68 = mrSges(6,1) * t111 + mrSges(6,2) * t417;
t67 = mrSges(7,1) * t111 - mrSges(7,3) * t417;
t66 = pkin(5) * t417 + qJ(6) * t111;
t64 = pkin(5) * t141 + qJ(6) * t142 + t133;
t62 = t137 * Ifges(5,4) + t136 * Ifges(5,2) + t158 * Ifges(5,6);
t30 = pkin(5) * t259 - t31;
t29 = -qJ(6) * t259 + t340;
t26 = -pkin(5) * t195 - t27;
t25 = qJ(6) * t195 + t28;
t16 = t89 * pkin(5) - t88 * qJ(6) + t142 * qJD(6) + t107;
t7 = -pkin(5) * t199 - t9;
t5 = qJ(6) * t199 - qJD(6) * t259 + t8;
t10 = [(t210 * mrSges(4,2) - t99 * mrSges(4,3) + Ifges(4,1) * t157 - Ifges(4,4) * t158 + Ifges(4,5) * qJDD(3) + t266 * t360 + t267 * t364 + t268 * t363 + t271 * t91) * t206 + (t406 * (-t247 * t310 + (t342 + t346) * t249) + t415 * t172 + t414 * t171 + (-m(4) * t342 + t416) * t249 + (m(4) * t229 - m(5) * (-qJ(4) * t232 - t229) + t232 * mrSges(5,3) + t386 + t406 * (-t229 - t209) - t419 - t438) * t247) * g(1) - (Ifges(5,5) * t137 + Ifges(5,6) * t136 + Ifges(5,3) * t158 + t423) * t259 / 0.2e1 + (mrSges(6,1) * t108 - mrSges(7,2) * t19 - mrSges(6,3) * t22 + t429) * t89 + (-mrSges(7,2) * t1 - mrSges(6,3) * t3 + t431) * t141 + (-t312 * t37 - t313 * t38 - t314 * t84 - t315 * t85) * mrSges(5,3) - (-m(4) * t99 + m(5) * t91 - qJDD(3) * mrSges(4,1) + mrSges(4,3) * t157 + t83) * t382 + (t266 * t356 + Ifges(4,4) * t357 + Ifges(4,1) * t354 + t146 / 0.2e1 + t435) * t198 + (-m(4) * t159 - t436) * t130 + (t420 / 0.2e1 + Ifges(5,6) * t404 + Ifges(5,5) * t405 + t426 * t358 + Ifges(5,3) * t356 - Ifges(4,2) * t357 - Ifges(4,4) * t354 + Ifges(7,6) * t368 + Ifges(6,6) * t369 - t145 / 0.2e1 + t397 * t366 - mrSges(4,3) * t160 - t434) * t199 + (-m(4) * t279 - m(5) * t217 + t439 * t309 + t406 * (t228 * t308 + t247 * t346 + t309 * t341 + t279) - t415 * t174 - t414 * t173 + (-t274 - t254 + t438) * t249 + t416 * t247) * g(2) + (mrSges(6,2) * t108 + mrSges(7,2) * t17 - mrSges(6,3) * t21 - t441) * t88 + t340 * t36 + m(6) * (t107 * t108 + t133 * t60 + t21 * t9 + t22 * t8 + t3 * t340 + t31 * t4) - t440 * t142 + Ifges(2,3) * qJDD(1) + m(3) * (-pkin(1) * t230 + (t221 + t297) * qJ(2) * t303) + t129 * t175 + t167 * (-qJDD(3) * mrSges(4,2) - mrSges(4,3) * t158) + t133 * t15 + t77 * t131 + t76 * t132 + 0.2e1 * t303 * t221 * mrSges(3,3) + (Ifges(3,4) * t240 + Ifges(3,2) * t242) * t295 + (Ifges(3,1) * t240 + Ifges(3,4) * t242) * t296 + t105 * t93 + t106 * t92 + t107 * t68 + t5 * t94 + t8 * t95 + t9 * t96 + t7 * t97 - t62 * t313 / 0.2e1 + t16 * t67 + t64 * t14 + t29 * t33 + t31 * t34 + t30 * t35 + m(5) * (t105 * t37 + t106 * t38 + t76 * t84 + t77 * t85) + m(4) * (t129 * t160 + t167 * t98 - t210 * t229) - t229 * t283 - pkin(1) * t276 + t230 * t275 + (-t210 * mrSges(4,1) - t37 * mrSges(5,1) + t38 * mrSges(5,2) + t98 * mrSges(4,3) + Ifges(4,4) * t157 - Ifges(5,5) * t363 - Ifges(4,2) * t158 + Ifges(4,6) * qJDD(3) - Ifges(5,6) * t364 - Ifges(6,6) * t375 - Ifges(7,6) * t374 - Ifges(5,3) * t360 - t426 * t361 - t397 * t376 + t411) * t259 + t312 * t371 + m(7) * (t1 * t29 + t16 * t39 + t17 * t7 + t19 * t5 + t2 * t30 + t6 * t64); t283 + t276 + m(3) * t230 + t241 * t93 + t239 * t92 - (t175 + t381) * t194 + t395 * t205 + t394 * t203 + (-t67 - t68 + t383) * t195 - t384 * t339 + t385 * t338 + (-g(1) * t247 + g(2) * t249) * (m(3) + m(4) + t299) - t286 * t303 * qJD(1) ^ 2 + (t1 * t205 + t17 * t385 - t19 * t384 - t195 * t39 + t2 * t203) * m(7) + (-t108 * t195 - t203 * t4 + t205 * t3 - t21 * t385 - t22 * t384) * m(6) + (-t149 * t195 + t194 * t263 + t239 * t38 + t241 * t37) * m(5) + (t159 * t195 - t160 * t194 + t210) * m(4); (-t1 * t203 - t17 * t384 - t19 * t385) * mrSges(7,2) + t422 * ((t341 * t406 + mrSges(4,2) - t285 + t439) * t234 + (mrSges(4,1) - m(7) * (-t228 - t265) + t257 + m(6) * t228 + t437) * t232) + (-t239 * t93 + t241 * t92) * qJ(4) + (-t203 * t3 + t21 * t384 - t22 * t385) * mrSges(6,3) + (-pkin(3) * t91 + qJ(4) * t380 - t263 * qJD(4) - t100 * t84 - t101 * t85) * m(5) + t381 * qJD(4) + (mrSges(6,1) * t385 - mrSges(6,2) * t384) * t108 - (Ifges(6,4) * t368 + Ifges(7,5) * t369 + t397 * t359 + t399 * t367 + t421) * t135 + (t1 * t166 + t155 * t6 - t260 * t2 + t390 * t39 + (t127 - t25) * t19 + (t128 - t26) * t17) * m(7) - t394 * t260 + (-t108 * t118 + t260 * t4 + t166 * t3 - t228 * t60 + (t127 - t28) * t22 + (-t128 - t27) * t21) * m(6) + t429 * t197 + t431 * t203 + (Ifges(5,5) * t239 + Ifges(5,6) * t241) * t360 + (Ifges(5,1) * t239 + t329) * t363 + (Ifges(5,2) * t241 + t330) * t364 + t62 * t349 + t145 * t354 - (-t266 * t357 + t435) * t194 + (t334 + t436) * t160 + (-t254 + t406 * (t209 + t310) + (-m(7) * t265 - t437) * t234 + t419) * g(3) + (-Ifges(4,2) * t356 + Ifges(5,3) * t357 + Ifges(6,6) * t368 + Ifges(7,6) * t369 - t391 / 0.2e1 - t392 / 0.2e1 + t397 * t367 + t426 * t359 + t434) * t195 - (Ifges(4,1) * t194 - t331 + t420) * t195 / 0.2e1 + (t316 * t84 + t317 * t85 + t380) * mrSges(5,3) + (t186 + t146) * t356 + t441 * t196 + (-Ifges(6,2) * t368 + Ifges(7,3) * t369 + t396 * t359 + t398 * t367 - t442) * t134 + t440 * t205 - t228 * t15 - t159 * t175 - Ifges(4,6) * t158 + Ifges(4,5) * t157 + t155 * t14 - t101 * t131 - t100 * t132 - t118 * t68 - t25 * t94 - t28 * t95 - t27 * t96 - t26 * t97 - t98 * mrSges(4,2) + t99 * mrSges(4,1) - pkin(3) * t83 + Ifges(4,3) * qJDD(3) + t390 * t67 + t338 * t128 + t339 * t127 + t395 * t166 + t239 * t371 + t91 * t272; -t338 * t417 + t339 * t111 - t170 * t131 + t169 * t132 + t14 + t15 + t83 + (t111 * t19 - t17 * t417 + t6) * m(7) + (t111 * t22 + t21 * t417 + t60) * m(6) + (t169 * t84 - t170 * t85 + t91) * m(5) + (t234 * g(3) - t232 * t422) * t299; (-t39 * t66 - pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t19 - g(1) * (-pkin(5) * t173 + qJ(6) * t174) - g(2) * (-pkin(5) * t171 + qJ(6) * t172) - (-pkin(5) * t231 + qJ(6) * t233) * t343) * m(7) + t423 - t411 + t19 * t336 + t17 * t337 - t108 * (mrSges(6,1) * t417 - mrSges(6,2) * t111) - t39 * (mrSges(7,1) * t417 + mrSges(7,3) * t111) + qJD(6) * t94 - t66 * t67 + qJ(6) * t33 - pkin(5) * t35 + t56 * t366 + (Ifges(7,3) * t417 - t327) * t369 + (-m(7) * t17 + t332 - t338) * t22 + (-m(7) * t19 - t333 - t339) * t21 + (-Ifges(6,2) * t417 - t110 + t393) * t368 + (-t111 * t397 + t396 * t417) * t359 + (-t111 * t399 + t109 - t328 + t53) * t367 + (t173 * t401 + t174 * t400) * g(1) + (t171 * t401 + t172 * t400) * g(2) + (t231 * t401 + t233 * t400) * t343; t417 * t67 - t190 * t94 + (-g(1) * t173 - g(2) * t171 - t19 * t190 - t231 * t343 + t39 * t417 + t2) * m(7) + t35;];
tau  = t10;
