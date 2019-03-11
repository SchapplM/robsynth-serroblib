% Calculate vector of inverse dynamics joint torques for
% S6RPRPRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRP6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP6_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP6_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP6_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP6_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP6_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:18:04
% EndTime: 2019-03-09 03:18:34
% DurationCPUTime: 21.57s
% Computational Cost: add. (6539->605), mult. (15354->704), div. (0->0), fcn. (10941->10), ass. (0->274)
t381 = Ifges(5,4) - Ifges(4,5);
t400 = Ifges(6,4) + Ifges(7,4);
t415 = Ifges(5,6) + Ifges(4,4);
t192 = sin(pkin(9));
t193 = cos(pkin(9));
t336 = cos(qJ(3));
t262 = t336 * t193;
t241 = qJD(1) * t262;
t248 = qJDD(1) * t336;
t197 = sin(qJ(3));
t277 = qJD(3) * t197;
t259 = t192 * t277;
t267 = qJDD(1) * t197;
t111 = qJD(1) * t259 - qJD(3) * t241 - t192 * t248 - t193 * t267;
t109 = qJDD(5) - t111;
t351 = t109 / 0.2e1;
t152 = t192 * t336 + t197 * t193;
t148 = t152 * qJD(3);
t112 = qJD(1) * t148 + t192 * t267 - t193 * t248;
t286 = t192 * t197;
t145 = qJD(1) * t286 - t241;
t196 = sin(qJ(5));
t199 = cos(qJ(5));
t119 = -qJD(3) * t196 + t145 * t199;
t59 = qJD(5) * t119 + qJDD(3) * t199 + t112 * t196;
t354 = t59 / 0.2e1;
t379 = Ifges(6,5) + Ifges(7,5);
t382 = Ifges(6,1) + Ifges(7,1);
t414 = t351 * t379 + t354 * t382;
t120 = qJD(3) * t199 + t145 * t196;
t60 = -qJD(5) * t120 - qJDD(3) * t196 + t112 * t199;
t353 = t60 / 0.2e1;
t413 = -mrSges(5,1) - mrSges(4,3);
t412 = -mrSges(5,2) + mrSges(4,1);
t411 = Ifges(6,2) + Ifges(7,2);
t399 = Ifges(7,6) + Ifges(6,6);
t398 = Ifges(7,3) + Ifges(6,3);
t278 = t192 ^ 2 + t193 ^ 2;
t271 = qJD(1) * qJD(2);
t170 = qJDD(1) * qJ(2) + t271;
t410 = t400 * t353 + t414;
t380 = Ifges(5,5) - Ifges(4,6);
t409 = t400 * t119;
t232 = mrSges(6,1) * t196 + mrSges(6,2) * t199;
t263 = m(7) * pkin(5) + mrSges(7,1);
t297 = t199 * mrSges(7,2);
t408 = -t196 * t263 - t232 - t297;
t194 = -qJ(6) - pkin(8);
t352 = pkin(3) + pkin(8);
t407 = m(6) * t352 + mrSges(6,3) - m(7) * (-pkin(3) + t194) + mrSges(7,3);
t406 = t400 * t120;
t405 = t400 * t199;
t404 = t400 * t196;
t403 = -m(5) - m(7);
t198 = sin(qJ(1));
t402 = g(2) * t198;
t401 = -mrSges(6,1) - mrSges(7,1);
t383 = mrSges(6,2) + mrSges(7,2);
t397 = t109 * t399 + t400 * t59 + t411 * t60;
t146 = t152 * qJD(1);
t137 = qJD(5) + t146;
t376 = t119 * t411 + t137 * t399 + t406;
t375 = t120 * t382 + t137 * t379 + t409;
t121 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t145;
t123 = mrSges(5,1) * t145 - qJD(3) * mrSges(5,3);
t395 = t121 - t123;
t370 = qJD(3) * t412 + t146 * t413;
t368 = t196 * t379 + t199 * t399;
t367 = t199 * t411 + t404;
t366 = t196 * t382 + t405;
t187 = pkin(9) + qJ(3);
t183 = sin(t187);
t184 = cos(t187);
t394 = -t412 * t184 + (mrSges(4,2) - mrSges(5,3)) * t183;
t393 = t109 * t398 + t379 * t59 + t399 * t60;
t321 = pkin(7) + qJ(2);
t159 = t321 * t192;
t153 = qJD(1) * t159;
t160 = t321 * t193;
t154 = qJD(1) * t160;
t113 = t336 * t153 + t154 * t197;
t215 = pkin(4) * t146 + t113;
t392 = qJD(4) + t215;
t360 = -t113 - qJD(4);
t133 = Ifges(4,4) * t145;
t391 = t146 * Ifges(4,1) + Ifges(4,5) * qJD(3) + t119 * t399 + t120 * t379 + t137 * t398 - t133;
t236 = -mrSges(3,1) * t193 + mrSges(3,2) * t192;
t388 = -m(3) * pkin(1) - mrSges(2,1) + t236 + t394;
t179 = pkin(2) * t193 + pkin(1);
t155 = -qJDD(1) * t179 + qJDD(2);
t205 = qJ(4) * t111 - qJD(4) * t146 + t155;
t21 = t112 * t352 + t205;
t275 = qJD(5) * t199;
t276 = qJD(5) * t196;
t245 = pkin(7) * qJDD(1) + t170;
t130 = t245 * t192;
t131 = t245 * t193;
t251 = qJD(3) * t336;
t45 = -t130 * t336 - t197 * t131 + t153 * t277 - t154 * t251;
t208 = qJDD(4) - t45;
t28 = -t111 * pkin(4) - qJDD(3) * t352 + t208;
t156 = -qJD(1) * t179 + qJD(2);
t207 = -qJ(4) * t146 + t156;
t64 = t145 * t352 + t207;
t70 = -qJD(3) * t352 + t392;
t3 = t196 * t28 + t199 * t21 + t70 * t275 - t276 * t64;
t2 = qJ(6) * t60 + qJD(6) * t119 + t3;
t25 = t196 * t70 + t199 * t64;
t4 = -qJD(5) * t25 - t196 * t21 + t199 * t28;
t387 = t4 * mrSges(6,1) - t3 * mrSges(6,2) - t2 * mrSges(7,2);
t114 = -t197 * t153 + t336 * t154;
t191 = qJD(3) * qJ(4);
t104 = -t191 - t114;
t84 = pkin(3) * t145 + t207;
t386 = t156 * mrSges(4,1) + t104 * mrSges(5,1) - t84 * mrSges(5,2) - t114 * mrSges(4,3);
t181 = pkin(5) * t199 + pkin(4);
t252 = m(3) * qJ(2) + mrSges(3,3);
t385 = -m(6) * (pkin(4) + t321) - m(7) * t181 + mrSges(2,2) - t252 + t413;
t103 = -qJD(3) * pkin(3) - t360;
t24 = -t196 * t64 + t199 * t70;
t15 = -qJ(6) * t120 + t24;
t14 = pkin(5) * t137 + t15;
t16 = qJ(6) * t119 + t25;
t384 = t103 * mrSges(5,1) + t24 * mrSges(6,1) + t14 * mrSges(7,1) + t156 * mrSges(4,2) - t25 * mrSges(6,2) - t16 * mrSges(7,2) + t113 * mrSges(4,3) - t84 * mrSges(5,3);
t330 = g(3) * t184;
t78 = -mrSges(6,1) * t119 + mrSges(6,2) * t120;
t374 = -t123 + t78;
t273 = qJD(6) * t199;
t281 = qJ(6) + t352;
t82 = -pkin(4) * t145 + t114;
t75 = t199 * t82;
t295 = qJ(4) * t145;
t76 = t146 * t352 + t295;
t373 = t276 * t281 - t273 + pkin(5) * t145 - t75 - (-qJ(6) * t146 - t76) * t196;
t158 = t281 * t199;
t293 = t146 * t199;
t32 = t196 * t82 + t199 * t76;
t372 = -qJ(6) * t293 - qJD(5) * t158 - qJD(6) * t196 - t32;
t371 = pkin(5) * t275 + t146 * t181 - t360;
t175 = t183 * qJ(4);
t200 = cos(qJ(1));
t287 = t184 * t200;
t369 = pkin(3) * t287 + t200 * t175;
t365 = -t275 - t293;
t294 = t146 * t196;
t364 = t276 + t294;
t363 = -t196 * t3 - t199 * t4;
t1 = pkin(5) * t109 - qJ(6) * t59 - qJD(6) * t120 + t4;
t362 = -t1 * t199 - t196 * t2;
t361 = g(1) * t200 + t402;
t272 = m(6) - t403;
t359 = -mrSges(6,1) - t263;
t358 = (-mrSges(5,3) + t408) * t184 + (m(5) * pkin(3) - mrSges(5,2) + t407) * t183;
t89 = t197 * (qJD(2) * t192 + qJD(3) * t160) - qJD(2) * t262 + t159 * t251;
t357 = qJ(4) * t272;
t316 = mrSges(7,2) * t196;
t231 = mrSges(7,1) * t199 - t316;
t233 = mrSges(6,1) * t199 - mrSges(6,2) * t196;
t73 = t191 + t82;
t43 = -pkin(5) * t119 + qJD(6) + t73;
t356 = t231 * t43 + t233 * t73;
t350 = -t119 / 0.2e1;
t349 = t119 / 0.2e1;
t348 = -t120 / 0.2e1;
t347 = t120 / 0.2e1;
t346 = -t137 / 0.2e1;
t345 = t137 / 0.2e1;
t344 = -t145 / 0.2e1;
t343 = t145 / 0.2e1;
t342 = -t146 / 0.2e1;
t341 = t146 / 0.2e1;
t333 = pkin(5) * t196;
t176 = t184 * pkin(3);
t325 = -qJD(3) / 0.2e1;
t324 = qJD(3) / 0.2e1;
t33 = mrSges(7,1) * t109 - mrSges(7,3) * t59;
t34 = mrSges(6,1) * t109 - mrSges(6,3) * t59;
t320 = -t33 - t34;
t35 = -mrSges(7,2) * t109 + mrSges(7,3) * t60;
t36 = -mrSges(6,2) * t109 + mrSges(6,3) * t60;
t319 = t35 + t36;
t151 = -t262 + t286;
t217 = -qJ(4) * t152 - t179;
t80 = t151 * t352 + t217;
t115 = t336 * t159 + t160 * t197;
t93 = pkin(4) * t152 + t115;
t40 = t196 * t93 + t199 * t80;
t313 = mrSges(7,3) * t119;
t85 = -mrSges(7,2) * t137 + t313;
t315 = mrSges(6,3) * t119;
t86 = -mrSges(6,2) * t137 + t315;
t318 = t85 + t86;
t312 = mrSges(7,3) * t120;
t87 = mrSges(7,1) * t137 - t312;
t314 = mrSges(6,3) * t120;
t88 = mrSges(6,1) * t137 - t314;
t317 = -t87 - t88;
t301 = t146 * Ifges(4,4);
t300 = t146 * Ifges(5,6);
t298 = t184 * mrSges(7,3);
t292 = t148 * t196;
t291 = t148 * t199;
t290 = t151 * t196;
t289 = t151 * t199;
t288 = t184 * t194;
t285 = t196 * t198;
t284 = t196 * t200;
t283 = t198 * t199;
t282 = t199 * t200;
t279 = t176 + t175;
t274 = qJD(5) * t352;
t269 = qJDD(1) * t192;
t268 = qJDD(1) * t193;
t77 = -mrSges(7,1) * t119 + mrSges(7,2) * t120;
t266 = -t77 - t374;
t18 = -t60 * mrSges(7,1) + t59 * mrSges(7,2);
t250 = -t276 / 0.2e1;
t180 = qJ(4) + t333;
t96 = -t111 * mrSges(5,1) + qJDD(3) * mrSges(5,2);
t247 = -qJ(6) * t151 - t80;
t244 = -t179 - t175;
t164 = t200 * t179;
t243 = t198 * t321 + t164;
t237 = -mrSges(3,1) * t268 + mrSges(3,2) * t269;
t220 = t196 * t24 - t199 * t25;
t147 = -t193 * t251 + t259;
t219 = qJ(4) * t147 - qJD(4) * t152;
t138 = t183 * t282 - t285;
t140 = t183 * t283 + t284;
t61 = t148 * t352 + t219;
t116 = -t197 * t159 + t160 * t336;
t90 = qJD(2) * t152 + qJD(3) * t116;
t66 = -t147 * pkin(4) + t90;
t7 = t196 * t66 + t199 * t61 + t93 * t275 - t276 * t80;
t214 = t151 * t275 + t292;
t213 = t151 * t276 - t291;
t44 = -t197 * t130 + t336 * t131 - t153 * t251 - t154 * t277;
t209 = -t196 * t318 + t199 * t317;
t41 = -qJDD(3) * qJ(4) - qJD(3) * qJD(4) - t44;
t29 = -pkin(4) * t112 - t41;
t206 = -qJD(5) * t220 - t363;
t65 = -pkin(4) * t148 - t89;
t182 = -qJDD(1) * pkin(1) + qJDD(2);
t157 = t281 * t196;
t141 = -t183 * t285 + t282;
t139 = t183 * t284 + t283;
t132 = Ifges(5,6) * t145;
t110 = pkin(3) * t151 + t217;
t108 = t111 * mrSges(4,2);
t107 = t111 * mrSges(5,3);
t106 = -mrSges(5,2) * t145 - mrSges(5,3) * t146;
t105 = pkin(3) * t146 + t295;
t99 = -t145 * Ifges(4,2) + Ifges(4,6) * qJD(3) + t301;
t98 = Ifges(5,4) * qJD(3) - t146 * Ifges(5,2) + t132;
t97 = Ifges(5,5) * qJD(3) + t145 * Ifges(5,3) - t300;
t95 = mrSges(5,1) * t112 - qJDD(3) * mrSges(5,3);
t94 = -t151 * pkin(4) + t116;
t92 = t199 * t93;
t79 = pkin(3) * t148 + t219;
t69 = -t151 * t181 + t116;
t63 = t199 * t66;
t42 = -qJDD(3) * pkin(3) + t208;
t39 = -t196 * t80 + t92;
t38 = pkin(3) * t112 + t205;
t37 = pkin(5) * t213 + t65;
t31 = -t196 * t76 + t75;
t30 = qJ(6) * t289 + t40;
t22 = pkin(5) * t152 + t196 * t247 + t92;
t19 = -mrSges(6,1) * t60 + mrSges(6,2) * t59;
t9 = -pkin(5) * t60 + qJDD(6) + t29;
t8 = -qJD(5) * t40 - t196 * t61 + t63;
t6 = -qJ(6) * t213 + t151 * t273 + t7;
t5 = -pkin(5) * t147 + t63 + t247 * t275 + (-qJ(6) * t148 - qJD(5) * t93 - qJD(6) * t151 - t61) * t196;
t10 = [0.2e1 * t278 * t170 * mrSges(3,3) + (t401 * t141 + t383 * t140 + ((-m(4) + t403) * t321 + t385) * t200 + (m(4) * t179 - m(5) * (t244 - t176) - m(6) * t244 - m(7) * (-t180 * t183 - t179) + t407 * t184 - t388) * t198) * g(1) + (qJD(5) * t375 + t397) * t289 / 0.2e1 + (-t213 * t399 + t214 * t379) * t345 + (t97 / 0.2e1 - t99 / 0.2e1 - Ifges(4,4) * t341 + Ifges(5,6) * t342 + Ifges(5,3) * t343 - Ifges(4,2) * t344 + t380 * t324 + t386) * t148 + (-t213 * t400 + t214 * t382) * t347 + (-t399 * t349 - t379 * t347 - t384 + t98 / 0.2e1 + t381 * t324 - Ifges(4,1) * t341 + Ifges(5,2) * t342 + Ifges(5,6) * t343 - Ifges(4,4) * t344 - t398 * t345 - t391 / 0.2e1) * t147 + (m(4) * t113 + m(5) * t103 - t370) * t90 + (-m(4) * t114 + m(5) * t104 - t395) * t89 + (m(4) * t44 - m(5) * t41 - qJDD(3) * mrSges(4,2) - mrSges(4,3) * t112 - t95) * t116 + t182 * t236 - pkin(1) * t237 + m(7) * (t1 * t22 + t14 * t5 + t16 * t6 + t2 * t30 + t37 * t43 + t69 * t9) + m(6) * (t24 * t8 + t25 * t7 + t29 * t94 + t3 * t40 + t39 * t4 + t65 * t73) + (-t213 * t411 + t214 * t400) * t349 + (t155 * mrSges(4,1) + t41 * mrSges(5,1) - t38 * mrSges(5,2) - t44 * mrSges(4,3) - t9 * t231 - t29 * t233 + t366 * t354 + t367 * t353 + t368 * t351 + (Ifges(5,3) + Ifges(4,2)) * t112 + t415 * t111 + t380 * qJDD(3) + t376 * t250) * t151 + t375 * t292 / 0.2e1 + t376 * t291 / 0.2e1 + m(3) * (-pkin(1) * t182 + (t170 + t271) * qJ(2) * t278) + t290 * t410 + (Ifges(3,4) * t192 + Ifges(3,2) * t193) * t268 + (Ifges(3,1) * t192 + Ifges(3,4) * t193) * t269 + (t42 * mrSges(5,1) + t1 * mrSges(7,1) + t155 * mrSges(4,2) - t45 * mrSges(4,3) - t38 * mrSges(5,3) + t398 * t351 + t399 * t353 + t379 * t354 + t387 + t393 / 0.2e1 - t415 * t112 + (-Ifges(4,1) - Ifges(5,2)) * t111 - t381 * qJDD(3)) * t152 + (-m(6) * (pkin(8) * t287 + t164 + t369) - mrSges(6,3) * t287 - m(4) * t243 + t403 * (t243 + t369) + t401 * t139 - t383 * t138 + t385 * t198 + (-m(7) * (t183 * t333 - t288) - t298 + t388) * t200) * g(2) + t43 * (mrSges(7,1) * t213 + mrSges(7,2) * t214) + t73 * (mrSges(6,1) * t213 + mrSges(6,2) * t214) + (-m(4) * t155 - t112 * mrSges(4,1) + t108) * t179 + Ifges(2,3) * qJDD(1) + (-m(4) * t45 + m(5) * t42 - qJDD(3) * mrSges(4,1) - mrSges(4,3) * t111 + t96) * t115 + t110 * (-t112 * mrSges(5,2) + t107) + t79 * t106 + t94 * t19 + t6 * t85 + t7 * t86 + t5 * t87 + t8 * t88 + t37 * t77 + t65 * t78 + t69 * t18 + t39 * t34 + t40 * t36 + t22 * t33 + t30 * t35 + m(5) * (t110 * t38 + t79 * t84) + (-t1 * t290 - t14 * t214 - t16 * t213 + t2 * t289) * mrSges(7,3) + (-t213 * t25 - t214 * t24 + t289 * t3 - t290 * t4) * mrSges(6,3); t237 + t412 * t112 + t107 - t108 + t319 * t199 + (t121 - t266) * t145 + m(3) * t182 + t320 * t196 + (t209 + t370) * t146 + t209 * qJD(5) + (-g(1) * t198 + g(2) * t200) * (m(3) + m(4) + t272) - t252 * t278 * qJD(1) ^ 2 + (-t1 * t196 + t145 * t43 + t199 * t2 - t137 * (t14 * t199 + t16 * t196)) * m(7) + (t145 * t73 - t196 * t4 + t199 * t3 - t137 * (t196 * t25 + t199 * t24)) * m(6) + (-t103 * t146 - t104 * t145 + t38) * m(5) + (-t113 * t146 + t114 * t145 + t155) * m(4); (t250 - t294 / 0.2e1) * t375 + (qJ(4) * t29 - t352 * t206 - t24 * t31 - t25 * t32 + t392 * t73) * m(6) - t404 * t354 + t405 * t353 + (-m(6) * (pkin(8) * t184 + t279) - m(7) * (t279 - t288) - t298 - m(5) * t279 + t408 * t183 + t394) * g(3) + (-t275 / 0.2e1 - t293 / 0.2e1) * t376 + (-Ifges(4,1) * t342 + Ifges(5,2) * t341 + t381 * t325 - t346 * t398 - t379 * t348 - t350 * t399 + t384) * t145 + t395 * t113 + (-Ifges(4,2) * t343 + Ifges(5,3) * t344 + t325 * t380 + t346 * t368 + t348 * t366 + t350 * t367 + t356 - t386) * t146 + (-t274 * t86 - t34 * t352 + t410 + t414) * t199 + t29 * t232 + (-t133 + t391) * t343 + (-t397 / 0.2e1 + mrSges(7,1) * t9 + t274 * t88 - t351 * t399 - t352 * t36 - t353 * t411) * t196 + t361 * (mrSges(4,1) * t183 + mrSges(4,2) * t184) + (-pkin(3) * t42 - qJ(4) * t41 - t103 * t114 + t104 * t360 - t105 * t84) * m(5) + (-t301 + t97) * t342 + t356 * qJD(5) + t371 * t77 + t372 * t85 + (-t1 * t158 + t14 * t373 - t157 * t2 + t16 * t372 + t180 * t9 + t371 * t43) * m(7) + t373 * t87 + t374 * qJD(4) + (t132 + t98) * t344 + (t24 * t364 + t25 * t365 - t330 + t363) * mrSges(6,3) + t380 * t112 + t381 * t111 + t370 * t114 + (t14 * t364 + t16 * t365 + t362) * mrSges(7,3) - (t119 * t367 + t120 * t366 + t137 * t368) * qJD(5) / 0.2e1 + (-t95 + t19) * qJ(4) + (Ifges(4,3) + Ifges(5,1)) * qJDD(3) + t180 * t18 - t157 * t35 - t158 * t33 + t9 * t297 - t105 * t106 - pkin(3) * t96 - t32 * t86 - t31 * t88 - t44 * mrSges(4,2) + t45 * mrSges(4,1) - t41 * mrSges(5,3) + t42 * mrSges(5,2) + (-t184 * t357 + t358) * t402 + t215 * t78 + (t200 * t358 - t287 * t357) * g(1) + (t300 + t99) * t341; t146 * t106 + t266 * qJD(3) + (t137 * t318 - t320) * t199 + (t137 * t317 + t319) * t196 + t96 + (-qJD(3) * t43 - t137 * (t14 * t196 - t16 * t199) - t362) * m(7) + (-qJD(3) * t73 - t146 * t220 + t206) * m(6) + (qJD(3) * t104 + t146 * t84 + t42) * m(5) + (-t183 * t361 + t330) * t272; (t88 + t314) * t25 + (-m(7) * (-t14 + t15) + t87 + t312) * t16 + (t33 + (-m(7) * t43 - t77) * t120) * pkin(5) + t263 * t1 + (t119 * t382 - t406) * t348 + (t138 * t359 + t139 * t383) * g(1) + (t140 * t359 - t141 * t383) * g(2) + t376 * t347 + (-t120 * t411 + t375 + t409) * t350 + (t119 * t379 - t120 * t399) * t346 + (-t86 + t315) * t24 + t387 - t73 * (mrSges(6,1) * t120 + mrSges(6,2) * t119) - t43 * (mrSges(7,1) * t120 + mrSges(7,2) * t119) - t15 * t85 + t14 * t313 + (t199 * t263 + t233 - t316) * t330 + t393; -t119 * t85 + t120 * t87 + (-g(3) * t183 - t16 * t119 + t14 * t120 - t184 * t361 + t9) * m(7) + t18;];
tau  = t10;
