% Calculate vector of inverse dynamics joint torques for
% S6RRPPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% Datum: 2019-03-09 08:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPRP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:29:41
% EndTime: 2019-03-09 08:30:18
% DurationCPUTime: 22.53s
% Computational Cost: add. (6908->670), mult. (15373->796), div. (0->0), fcn. (10566->10), ass. (0->312)
t415 = Ifges(4,5) - Ifges(5,4);
t313 = qJDD(2) / 0.2e1;
t385 = 0.2e1 * t313;
t421 = -m(7) - m(5);
t287 = m(6) - t421;
t416 = Ifges(6,4) + Ifges(7,4);
t435 = -Ifges(5,6) - Ifges(4,4);
t210 = sin(qJ(2));
t286 = qJD(1) * qJD(2);
t268 = t210 * t286;
t213 = cos(qJ(2));
t285 = qJDD(1) * t213;
t169 = -t268 + t285;
t170 = qJDD(1) * t210 + t213 * t286;
t206 = sin(pkin(9));
t312 = cos(pkin(9));
t123 = -t169 * t312 + t170 * t206;
t261 = t312 * t213;
t295 = qJD(1) * t210;
t149 = -qJD(1) * t261 + t206 * t295;
t209 = sin(qJ(5));
t212 = cos(qJ(5));
t129 = -qJD(2) * t209 + t149 * t212;
t60 = qJD(5) * t129 + qJDD(2) * t212 + t123 * t209;
t374 = t60 / 0.2e1;
t130 = qJD(2) * t212 + t149 * t209;
t61 = -qJD(5) * t130 - qJDD(2) * t209 + t123 * t212;
t373 = t61 / 0.2e1;
t124 = t206 * t169 + t170 * t312;
t120 = qJDD(5) + t124;
t371 = t120 / 0.2e1;
t434 = -mrSges(5,1) - mrSges(4,3);
t433 = -mrSges(5,2) + mrSges(4,1);
t417 = Ifges(6,1) + Ifges(7,1);
t414 = Ifges(6,5) + Ifges(7,5);
t432 = Ifges(6,2) + Ifges(7,2);
t411 = Ifges(7,6) + Ifges(6,6);
t410 = Ifges(7,3) + Ifges(6,3);
t247 = mrSges(6,1) * t209 + mrSges(6,2) * t212;
t277 = m(7) * pkin(5) + mrSges(7,1);
t314 = t212 * mrSges(7,2);
t431 = -t209 * t277 - t247 - t314;
t207 = -qJ(6) - pkin(8);
t372 = pkin(3) + pkin(8);
t430 = m(6) * t372 + mrSges(6,3) - m(7) * (-pkin(3) + t207) + mrSges(7,3);
t429 = t414 * t371 + t416 * t373 + t417 * t374;
t211 = sin(qJ(1));
t428 = g(2) * t211;
t412 = -Ifges(4,6) + Ifges(5,5);
t427 = t416 * t129;
t176 = -mrSges(3,1) * t213 + mrSges(3,2) * t210;
t201 = qJ(2) + pkin(9);
t198 = sin(t201);
t199 = cos(t201);
t426 = t176 - t433 * t199 - (-mrSges(4,2) + mrSges(5,3)) * t198;
t425 = t416 * t130;
t424 = t416 * t212;
t423 = t416 * t209;
t208 = -qJ(3) - pkin(7);
t177 = t208 * t213;
t168 = qJD(1) * t177;
t153 = t206 * t168;
t175 = t208 * t210;
t167 = qJD(1) * t175;
t122 = t167 * t312 + t153;
t422 = -qJD(4) + t122;
t420 = t169 / 0.2e1;
t341 = qJD(2) / 0.2e1;
t419 = -mrSges(7,1) - mrSges(6,1);
t418 = mrSges(6,2) + mrSges(7,2);
t409 = t120 * t411 + t416 * t60 + t432 * t61;
t200 = t213 * pkin(2);
t193 = t200 + pkin(1);
t165 = t206 * t213 + t210 * t312;
t151 = t165 * qJD(1);
t142 = qJD(5) + t151;
t407 = t129 * t432 + t142 * t411 + t425;
t406 = t130 * t417 + t142 * t414 + t427;
t288 = qJD(6) * t212;
t290 = qJD(5) * t209;
t275 = t312 * pkin(2);
t191 = -t275 - pkin(3);
t186 = -pkin(8) + t191;
t296 = qJ(6) - t186;
t196 = pkin(2) * t295;
t258 = qJ(4) * t149 + t196;
t69 = t151 * t372 + t258;
t262 = t312 * t168;
t121 = t167 * t206 - t262;
t354 = pkin(4) * t149;
t81 = t121 - t354;
t79 = t212 * t81;
t405 = t290 * t296 - t288 + pkin(5) * t149 - t79 - (-qJ(6) * t151 - t69) * t209;
t157 = t296 * t212;
t307 = t151 * t212;
t32 = t209 * t81 + t212 * t69;
t404 = -qJ(6) * t307 - qJD(5) * t157 - qJD(6) * t209 - t32;
t334 = mrSges(5,1) * t149;
t135 = -qJD(2) * mrSges(5,3) + t334;
t77 = -mrSges(6,1) * t129 + mrSges(6,2) * t130;
t403 = -t135 + t77;
t192 = pkin(5) * t212 + pkin(4);
t289 = qJD(5) * t212;
t402 = pkin(5) * t289 + t151 * t192 - t422;
t401 = qJD(2) * t433 + t151 * t434;
t318 = t149 * mrSges(4,3);
t133 = -qJD(2) * mrSges(4,2) - t318;
t400 = t135 - t133;
t187 = t198 * qJ(4);
t214 = cos(qJ(1));
t301 = t199 * t214;
t399 = pkin(3) * t301 + t214 * t187;
t398 = t209 * t414 + t212 * t411;
t397 = t212 * t432 + t423;
t396 = t209 * t417 + t424;
t394 = t120 * t410 + t411 * t61 + t414 * t60;
t393 = -t289 - t307;
t308 = t151 * t209;
t392 = t290 + t308;
t194 = pkin(7) * t285;
t162 = -pkin(7) * t268 + t194;
t163 = t170 * pkin(7);
t391 = t162 * t213 + t163 * t210;
t311 = qJDD(1) * pkin(1);
t139 = -pkin(2) * t169 + qJDD(3) - t311;
t218 = -qJ(4) * t124 - qJD(4) * t151 + t139;
t22 = t123 * t372 + t218;
t292 = qJD(3) * t210;
t109 = qJDD(2) * pkin(2) - qJ(3) * t170 - qJD(1) * t292 - t163;
t293 = qJD(2) * t210;
t280 = pkin(7) * t293;
t291 = qJD(3) * t213;
t117 = qJ(3) * t169 + t194 + (-t280 + t291) * qJD(1);
t45 = t109 * t312 - t206 * t117;
t229 = qJDD(4) - t45;
t28 = t124 * pkin(4) - qJDD(2) * t372 + t229;
t172 = -qJD(1) * t193 + qJD(3);
t220 = -qJ(4) * t151 + t172;
t63 = t149 * t372 + t220;
t161 = qJD(2) * pkin(2) + t167;
t113 = t161 * t312 + t153;
t230 = qJD(4) - t113;
t346 = t151 * pkin(4);
t68 = -qJD(2) * t372 + t230 + t346;
t3 = t209 * t28 + t212 * t22 + t68 * t289 - t290 * t63;
t25 = t209 * t68 + t212 * t63;
t4 = -qJD(5) * t25 - t209 * t22 + t212 * t28;
t390 = -t209 * t3 - t212 * t4;
t1 = pkin(5) * t120 - qJ(6) * t60 - qJD(6) * t130 + t4;
t2 = qJ(6) * t61 + qJD(6) * t129 + t3;
t389 = -t1 * t212 - t2 * t209;
t388 = g(1) * t214 + t428;
t387 = -mrSges(6,1) - t277;
t141 = Ifges(4,4) * t149;
t386 = t151 * Ifges(4,1) + Ifges(4,5) * qJD(2) + t129 * t411 + t130 * t414 + t142 * t410 - t141;
t355 = pkin(2) * t210;
t384 = t287 * t355 + (-mrSges(5,3) + t431) * t199 + (m(5) * pkin(3) - mrSges(5,2) + t430) * t198;
t383 = qJ(4) * t287;
t382 = -m(3) * pkin(1) - mrSges(2,1) + t426;
t333 = mrSges(7,2) * t209;
t246 = mrSges(7,1) * t212 - t333;
t248 = mrSges(6,1) * t212 - mrSges(6,2) * t209;
t114 = t206 * t161 - t262;
t102 = -qJD(2) * qJ(4) - t114;
t74 = -t102 - t354;
t42 = -pkin(5) * t129 + qJD(6) + t74;
t381 = t246 * t42 + t248 * t74;
t380 = t4 * mrSges(6,1) - t3 * mrSges(6,2) - t2 * mrSges(7,2);
t88 = pkin(3) * t149 + t220;
t379 = t172 * mrSges(4,1) + t102 * mrSges(5,1) - t88 * mrSges(5,2) - t114 * mrSges(4,3);
t378 = -m(3) * pkin(7) - m(6) * (pkin(4) - t208) - m(7) * t192 + mrSges(2,2) - mrSges(3,3) + t434;
t24 = -t209 * t63 + t212 * t68;
t15 = -qJ(6) * t130 + t24;
t10 = pkin(5) * t142 + t15;
t16 = qJ(6) * t129 + t25;
t377 = t24 * mrSges(6,1) + t10 * mrSges(7,1) + t172 * mrSges(4,2) - t25 * mrSges(6,2) - t16 * mrSges(7,2) - t88 * mrSges(5,3);
t370 = -t129 / 0.2e1;
t369 = t129 / 0.2e1;
t368 = -t130 / 0.2e1;
t367 = t130 / 0.2e1;
t366 = -t142 / 0.2e1;
t365 = t142 / 0.2e1;
t364 = -t149 / 0.2e1;
t363 = t149 / 0.2e1;
t362 = -t151 / 0.2e1;
t361 = t151 / 0.2e1;
t356 = pkin(2) * t206;
t352 = pkin(5) * t209;
t351 = pkin(7) * t213;
t348 = g(3) * t199;
t189 = t199 * pkin(3);
t342 = -qJD(2) / 0.2e1;
t33 = mrSges(7,1) * t120 - mrSges(7,3) * t60;
t34 = mrSges(6,1) * t120 - mrSges(6,3) * t60;
t338 = -t33 - t34;
t35 = -mrSges(7,2) * t120 + mrSges(7,3) * t61;
t36 = -mrSges(6,2) * t120 + mrSges(6,3) * t61;
t337 = t35 + t36;
t164 = t206 * t210 - t261;
t252 = -qJ(4) * t165 - t193;
t80 = t164 * t372 + t252;
t125 = -t312 * t175 - t177 * t206;
t94 = pkin(4) * t165 + t125;
t39 = t209 * t94 + t212 * t80;
t330 = mrSges(7,3) * t129;
t83 = -mrSges(7,2) * t142 + t330;
t332 = mrSges(6,3) * t129;
t84 = -mrSges(6,2) * t142 + t332;
t336 = t83 + t84;
t329 = mrSges(7,3) * t130;
t85 = mrSges(7,1) * t142 - t329;
t331 = mrSges(6,3) * t130;
t86 = mrSges(6,1) * t142 - t331;
t335 = -t85 - t86;
t328 = Ifges(3,4) * t210;
t327 = Ifges(3,4) * t213;
t317 = t151 * Ifges(4,4);
t316 = t151 * Ifges(5,6);
t315 = t199 * mrSges(7,3);
t46 = t206 * t109 + t312 * t117;
t150 = t165 * qJD(2);
t310 = t150 * t209;
t309 = t150 * t212;
t304 = t164 * t209;
t303 = t164 * t212;
t302 = t199 * t207;
t300 = t209 * t214;
t299 = t211 * t209;
t298 = t211 * t212;
t297 = t212 * t214;
t294 = qJD(1) * t213;
t284 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t295) * t351;
t76 = -mrSges(7,1) * t129 + mrSges(7,2) * t130;
t283 = -t76 - t403;
t197 = pkin(2) * t293;
t276 = t189 + t187 + t200;
t19 = -t61 * mrSges(7,1) + t60 * mrSges(7,2);
t266 = -t290 / 0.2e1;
t264 = -qJ(4) - t352;
t263 = qJD(2) * t208;
t106 = t124 * mrSges(5,1) + qJDD(2) * mrSges(5,2);
t259 = -qJ(6) * t164 - t80;
t257 = -t193 - t187;
t147 = t210 * t263 + t291;
t148 = t213 * t263 - t292;
t92 = t147 * t206 - t312 * t148;
t182 = t214 * t193;
t256 = -t211 * t208 + t182;
t41 = -qJDD(2) * qJ(4) - qJD(2) * qJD(4) - t46;
t250 = mrSges(3,1) * t210 + mrSges(3,2) * t213;
t242 = t213 * Ifges(3,2) + t328;
t239 = Ifges(3,5) * t213 - Ifges(3,6) * t210;
t233 = t209 * t24 - t212 * t25;
t231 = pkin(1) * t250;
t143 = t198 * t297 - t299;
t145 = t198 * t298 + t300;
t152 = qJD(2) * t261 - t206 * t293;
t226 = -qJ(4) * t152 - qJD(4) * t165 + t197;
t47 = t150 * t372 + t226;
t70 = pkin(4) * t152 + t92;
t7 = t209 * t70 + t212 * t47 + t94 * t289 - t290 * t80;
t29 = -pkin(4) * t123 - t41;
t228 = t164 * t289 + t310;
t227 = t164 * t290 - t309;
t225 = t210 * (Ifges(3,1) * t213 - t328);
t93 = t147 * t312 + t206 * t148;
t126 = t206 * t175 - t177 * t312;
t222 = -t209 * t336 + t212 * t335;
t71 = -t150 * pkin(4) + t93;
t217 = -qJD(5) * t233 - t390;
t195 = Ifges(3,4) * t294;
t188 = qJ(4) + t356;
t174 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t294;
t171 = -t264 + t356;
t159 = Ifges(3,1) * t295 + Ifges(3,5) * qJD(2) + t195;
t158 = Ifges(3,6) * qJD(2) + qJD(1) * t242;
t156 = t296 * t209;
t146 = -t198 * t299 + t297;
t144 = t198 * t300 + t298;
t140 = Ifges(5,6) * t149;
t119 = t124 * mrSges(5,3);
t118 = t124 * mrSges(4,2);
t112 = pkin(3) * t164 + t252;
t111 = -mrSges(5,2) * t149 - mrSges(5,3) * t151;
t110 = mrSges(4,1) * t149 + mrSges(4,2) * t151;
t105 = mrSges(5,1) * t123 - qJDD(2) * mrSges(5,3);
t104 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t124;
t103 = -qJDD(2) * mrSges(4,2) - mrSges(4,3) * t123;
t100 = -t149 * Ifges(4,2) + Ifges(4,6) * qJD(2) + t317;
t99 = Ifges(5,4) * qJD(2) - t151 * Ifges(5,2) + t140;
t98 = Ifges(5,5) * qJD(2) + t149 * Ifges(5,3) - t316;
t97 = -qJD(2) * pkin(3) + t230;
t95 = -t164 * pkin(4) + t126;
t91 = pkin(3) * t151 + t258;
t90 = t212 * t94;
t82 = t122 - t346;
t73 = pkin(3) * t150 + t226;
t72 = -t164 * t192 + t126;
t67 = t212 * t70;
t43 = -qJDD(2) * pkin(3) + t229;
t40 = pkin(5) * t227 + t71;
t38 = -t209 * t80 + t90;
t37 = pkin(3) * t123 + t218;
t31 = -t209 * t69 + t79;
t30 = qJ(6) * t303 + t39;
t23 = pkin(5) * t165 + t209 * t259 + t90;
t20 = -mrSges(6,1) * t61 + mrSges(6,2) * t60;
t9 = -pkin(5) * t61 + qJDD(6) + t29;
t8 = -qJD(5) * t39 - t209 * t47 + t67;
t6 = -qJ(6) * t227 + t164 * t288 + t7;
t5 = pkin(5) * t152 + t67 + t259 * t289 + (-qJ(6) * t150 - qJD(5) * t94 - qJD(6) * t164 - t47) * t209;
t11 = [(-t227 * t416 + t228 * t417) * t367 + (t239 * t341 - t284) * qJD(2) + Ifges(3,6) * t213 * t313 + t170 * t327 / 0.2e1 + (-t1 * t304 - t10 * t228 - t16 * t227 + t2 * t303) * mrSges(7,3) + (-m(6) * (pkin(8) * t301 + t182 + t399) - mrSges(6,3) * t301 - m(4) * t256 + t421 * (t256 + t399) + t419 * t144 - t418 * t143 + t378 * t211 + (-m(7) * (t198 * t352 - t302) - t315 + t382) * t214) * g(2) + (m(4) * t114 - m(5) * t102 - t400) * t93 + (-m(4) * t113 + m(5) * t97 - t401) * t92 + (t169 * t351 + t391) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t391) - qJDD(2) * mrSges(3,2) * t351 + (t386 / 0.2e1 - t99 / 0.2e1 + Ifges(4,1) * t361 - Ifges(5,2) * t362 - Ifges(5,6) * t363 + Ifges(4,4) * t364 - t113 * mrSges(4,3) + t97 * mrSges(5,1) + t411 * t369 + t414 * t367 + t410 * t365 + t415 * t341 + t377) * t152 + (t139 * mrSges(4,1) + t41 * mrSges(5,1) - t37 * mrSges(5,2) - t46 * mrSges(4,3) - t9 * t246 - t29 * t248 + t396 * t374 + t397 * t373 + t398 * t371 + t435 * t124 + (Ifges(5,3) + Ifges(4,2)) * t123 + t412 * t385 + t407 * t266) * t164 + (t419 * t146 + t418 * t145 + ((m(4) - t421) * t208 + t378) * t214 + (-m(5) * (t257 - t189) - m(6) * t257 - m(7) * (t198 * t264 - t193) + m(4) * t193 + t430 * t199 - t382) * t211) * g(1) + (Ifges(3,1) * t170 + Ifges(3,4) * t420 - pkin(7) * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t170) + t385 * Ifges(3,5)) * t210 + (t225 + t213 * (-Ifges(3,2) * t210 + t327)) * t286 / 0.2e1 + (-t100 / 0.2e1 + t98 / 0.2e1 - Ifges(4,4) * t361 + Ifges(5,6) * t362 + Ifges(5,3) * t363 - Ifges(4,2) * t364 + t412 * t341 + t379) * t150 + (-t227 * t411 + t228 * t414) * t365 + m(4) * (-t139 * t193 + t172 * t197) + m(5) * (t112 * t37 + t73 * t88) + (-t227 * t25 - t228 * t24 + t3 * t303 - t304 * t4) * mrSges(6,3) - t176 * t311 - t158 * t293 / 0.2e1 - t231 * t286 - t174 * t280 + (t394 / 0.2e1 + t43 * mrSges(5,1) + t1 * mrSges(7,1) + t139 * mrSges(4,2) - t45 * mrSges(4,3) - t37 * mrSges(5,3) + t410 * t371 + t411 * t373 + t414 * t374 + t380 + (Ifges(4,1) + Ifges(5,2)) * t124 + t435 * t123 + t415 * t385) * t165 + t406 * t310 / 0.2e1 + t407 * t309 / 0.2e1 + (m(4) * t46 - m(5) * t41 + t103 - t105) * t126 + (qJD(5) * t406 + t409) * t303 / 0.2e1 + t242 * t420 + t112 * (-t123 * mrSges(5,2) - t119) + (-m(4) * t45 + m(5) * t43 - t104 + t106) * t125 + m(7) * (t1 * t23 + t10 * t5 + t16 * t6 + t2 * t30 + t40 * t42 + t72 * t9) + m(6) * (t24 * t8 + t25 * t7 + t29 * t95 + t3 * t39 + t38 * t4 + t71 * t74) + Ifges(2,3) * qJDD(1) + t213 * (Ifges(3,4) * t170 + Ifges(3,2) * t169 + Ifges(3,6) * qJDD(2)) / 0.2e1 - t193 * (t123 * mrSges(4,1) + t118) - pkin(1) * (-mrSges(3,1) * t169 + mrSges(3,2) * t170) + t73 * t111 + t95 * t20 + t6 * t83 + t7 * t84 + t5 * t85 + t8 * t86 + t40 * t76 + t71 * t77 + t72 * t19 + t30 * t35 + t38 * t34 + t39 * t36 + t23 * t33 + t110 * t197 + t213 * t159 * t341 + t304 * t429 + t74 * (mrSges(6,1) * t227 + mrSges(6,2) * t228) + t42 * (mrSges(7,1) * t227 + mrSges(7,2) * t228) + (-t227 * t432 + t228 * t416) * t369; (t140 + t99) * t364 + (t284 + (t231 - t225 / 0.2e1) * qJD(1)) * qJD(1) + (t10 * t392 + t16 * t393 + t389) * mrSges(7,3) + (t24 * t392 + t25 * t393 + t390) * mrSges(6,3) - (t129 * t397 + t130 * t396 + t142 * t398) * qJD(5) / 0.2e1 + t400 * t122 + t401 * t121 + (-t141 + t386) * t363 + (t102 * t422 - t121 * t97 - t188 * t41 + t191 * t43 - t88 * t91) * m(5) + (-m(4) * t200 - m(7) * (t276 - t302) - t315 - m(6) * (pkin(8) * t199 + t276) - t199 * mrSges(6,3) - m(5) * t276 + t431 * t198 + t426) * g(3) - (-Ifges(3,2) * t295 + t159 + t195) * t294 / 0.2e1 + (m(6) * t217 + t209 * t36 + t212 * t34 + t289 * t84 - t290 * t86) * t186 + (t316 + t100) * t361 + (-Ifges(4,2) * t363 + Ifges(5,3) * t364 + t342 * t412 + t366 * t398 + t368 * t396 + t370 * t397 - t379 + t381) * t151 + t412 * t123 + (-t209 * t411 + t212 * t414) * t371 + t415 * t124 + (-Ifges(4,1) * t362 + Ifges(5,2) * t361 - t342 * t415 - t366 * t410 - t368 * t414 - t370 * t411 + t377) * t149 - t113 * t318 + t9 * (mrSges(7,1) * t209 + t314) + t381 * qJD(5) + (Ifges(4,3) + Ifges(3,3) + Ifges(5,1)) * qJDD(2) + (-t317 + t98) * t362 - t239 * t286 / 0.2e1 - t110 * t196 + (t113 * t121 - t114 * t122 - t172 * t196 + (t206 * t46 + t312 * t45) * pkin(2)) * m(4) + (t214 * t384 - t301 * t383) * g(1) - t409 * t209 / 0.2e1 + t402 * t76 + t403 * qJD(4) + t404 * t83 + (-t1 * t157 + t10 * t405 - t156 * t2 + t16 * t404 + t171 * t9 + t402 * t42) * m(7) + t405 * t85 + (t188 * t29 - t24 * t31 - t25 * t32 + (-t82 + qJD(4)) * t74) * m(6) + (t20 - t105) * t188 + (-t308 / 0.2e1 + t266) * t406 + (m(4) * t355 + mrSges(4,1) * t198 + mrSges(4,2) * t199 + t250) * t388 + t191 * t106 + t171 * t19 + Ifges(3,6) * t169 + Ifges(3,5) * t170 - t162 * mrSges(3,2) - t163 * mrSges(3,1) - t156 * t35 - t157 * t33 - t91 * t111 - t82 * t77 - t32 * t84 - t31 * t86 - t46 * mrSges(4,2) - t41 * mrSges(5,3) + t43 * mrSges(5,2) + t45 * mrSges(4,1) + (t158 / 0.2e1 + pkin(7) * t174) * t295 + t29 * t247 + t104 * t275 + (-t199 * t383 + t384) * t428 + t212 * t429 + t97 * t334 + t103 * t356 + (-t289 / 0.2e1 - t307 / 0.2e1) * t407 + (-t209 * t432 + t424) * t373 + (t212 * t417 - t423) * t374; t118 - t119 + t337 * t212 + t338 * t209 + t433 * t123 + (t133 - t283) * t149 + t222 * qJD(5) + (t222 + t401) * t151 + (-g(1) * t211 + g(2) * t214) * (m(4) + t287) + (-t1 * t209 + t149 * t42 + t2 * t212 - t142 * (t10 * t212 + t16 * t209)) * m(7) + (t149 * t74 - t209 * t4 + t212 * t3 - t142 * (t209 * t25 + t212 * t24)) * m(6) + (-t102 * t149 - t151 * t97 + t37) * m(5) + (t113 * t151 + t114 * t149 + t139) * m(4); t151 * t111 + t283 * qJD(2) + (t142 * t336 - t338) * t212 + (t142 * t335 + t337) * t209 + t106 + (-qJD(2) * t42 - t142 * (t10 * t209 - t16 * t212) - t389) * m(7) + (-qJD(2) * t74 - t151 * t233 + t217) * m(6) + (qJD(2) * t102 + t151 * t88 + t43) * m(5) + (-t198 * t388 + t348) * t287; (-t130 * t432 + t406 + t427) * t370 + (t143 * t387 + t144 * t418) * g(1) + (t145 * t387 - t146 * t418) * g(2) + (t129 * t414 - t130 * t411) * t366 + (t129 * t417 - t425) * t368 + t277 * t1 + t380 + t407 * t367 + (-m(7) * (-t10 + t15) + t85 + t329) * t16 + (t212 * t277 + t248 - t333) * t348 + (t86 + t331) * t25 + (-t84 + t332) * t24 - t74 * (mrSges(6,1) * t130 + mrSges(6,2) * t129) - t42 * (mrSges(7,1) * t130 + mrSges(7,2) * t129) - t15 * t83 + t10 * t330 + t394 + ((-m(7) * t42 - t76) * t130 + t33) * pkin(5); -t129 * t83 + t130 * t85 + (-g(3) * t198 + t10 * t130 - t129 * t16 - t199 * t388 + t9) * m(7) + t19;];
tau  = t11;
