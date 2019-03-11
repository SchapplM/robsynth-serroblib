% Calculate vector of inverse dynamics joint torques for
% S6RPPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRP4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:05:10
% EndTime: 2019-03-09 02:05:32
% DurationCPUTime: 13.45s
% Computational Cost: add. (4696->559), mult. (8581->719), div. (0->0), fcn. (4865->8), ass. (0->256)
t419 = Ifges(6,1) + Ifges(7,1);
t401 = Ifges(7,4) + Ifges(6,5);
t418 = Ifges(6,6) - Ifges(7,6);
t172 = sin(pkin(9));
t177 = cos(qJ(4));
t315 = cos(pkin(9));
t243 = qJD(1) * t315;
t175 = sin(qJ(4));
t294 = qJD(4) * t175;
t420 = -t172 * t294 - t177 * t243;
t174 = sin(qJ(5));
t176 = cos(qJ(5));
t293 = qJD(4) * t176;
t299 = qJD(1) * t175;
t127 = t174 * t299 + t293;
t351 = t127 / 0.2e1;
t285 = t174 * qJD(4);
t298 = qJD(1) * t176;
t128 = t175 * t298 - t285;
t350 = -t128 / 0.2e1;
t297 = qJD(1) * t177;
t152 = qJD(5) + t297;
t347 = t152 / 0.2e1;
t402 = mrSges(6,3) + mrSges(7,2);
t417 = Ifges(6,3) + Ifges(7,2);
t208 = -pkin(5) * t176 - qJ(6) * t174;
t205 = pkin(4) - t208;
t225 = -t176 * mrSges(7,1) - t174 * mrSges(7,3);
t227 = -mrSges(6,1) * t176 + mrSges(6,2) * t174;
t416 = m(6) * pkin(4) + m(7) * t205 - t225 - t227;
t384 = -t174 * t418 + t176 * t401;
t322 = Ifges(7,5) * t174;
t325 = Ifges(6,4) * t174;
t382 = t176 * t419 + t322 - t325;
t307 = t174 * t177;
t111 = t172 * t307 + t176 * t315;
t245 = t177 * t315;
t260 = t175 * t293;
t394 = -qJD(5) * t111 - t172 * t260 - (t172 * t174 + t176 * t245) * qJD(1);
t305 = t176 * t177;
t112 = t172 * t305 - t174 * t315;
t393 = qJD(5) * t112 + t172 * t298 + t420 * t174;
t178 = -pkin(1) - pkin(2);
t147 = qJDD(1) * t178 + qJDD(2);
t284 = qJD(1) * qJD(2);
t149 = qJDD(1) * qJ(2) + t284;
t83 = t172 * t147 + t315 * t149;
t415 = -qJDD(1) * pkin(7) + qJD(3) * qJD(4) + t83;
t272 = mrSges(5,3) * t299;
t389 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t127 - mrSges(6,2) * t128 - t272;
t295 = qJD(3) * t175;
t148 = qJD(1) * t178 + qJD(2);
t96 = qJ(2) * t243 + t172 * t148;
t90 = -qJD(1) * pkin(7) + t96;
t70 = t177 * t90 + t295;
t64 = qJD(4) * pkin(8) + t70;
t233 = t177 * pkin(4) + t175 * pkin(8);
t300 = qJD(1) * t172;
t95 = -qJ(2) * t300 + t148 * t315;
t89 = qJD(1) * pkin(3) - t95;
t65 = qJD(1) * t233 + t89;
t19 = -t174 * t64 + t176 * t65;
t20 = t174 * t65 + t176 * t64;
t23 = t175 * qJDD(3) + t177 * t415 - t294 * t90;
t21 = qJDD(4) * pkin(8) + t23;
t288 = qJD(5) * t176;
t290 = qJD(5) * t174;
t283 = qJD(1) * qJD(4);
t134 = -qJDD(1) * t177 + t175 * t283;
t135 = -qJDD(1) * t175 - t177 * t283;
t82 = t147 * t315 - t172 * t149;
t77 = qJDD(1) * pkin(3) - t82;
t30 = -t134 * pkin(4) - t135 * pkin(8) + t77;
t3 = t174 * t30 + t176 * t21 + t65 * t288 - t290 * t64;
t4 = -qJD(5) * t20 - t174 * t21 + t176 * t30;
t230 = -t174 * t4 + t176 * t3;
t414 = -t19 * t288 - t20 * t290 + t230;
t14 = -pkin(5) * t152 + qJD(6) - t19;
t15 = qJ(6) * t152 + t20;
t124 = qJDD(5) - t134;
t1 = qJ(6) * t124 + qJD(6) * t152 + t3;
t2 = -pkin(5) * t124 + qJDD(6) - t4;
t231 = t1 * t176 + t174 * t2;
t413 = t14 * t288 - t15 * t290 + t231;
t289 = qJD(5) * t175;
t61 = -t174 * t135 - qJD(5) * t285 + (qJD(1) * t289 + qJDD(4)) * t176;
t33 = -mrSges(6,2) * t124 + mrSges(6,3) * t61;
t34 = mrSges(7,2) * t61 + mrSges(7,3) * t124;
t397 = t33 + t34;
t291 = qJD(5) * t127;
t60 = qJDD(4) * t174 + t135 * t176 + t291;
t31 = mrSges(6,1) * t124 - mrSges(6,3) * t60;
t32 = -t124 * mrSges(7,1) + t60 * mrSges(7,2);
t398 = -t31 + t32;
t412 = t398 * t174 + t397 * t176;
t327 = Ifges(5,4) * t177;
t223 = -t175 * Ifges(5,1) - t327;
t123 = Ifges(6,4) * t127;
t323 = Ifges(7,5) * t127;
t395 = -t128 * t419 + t152 * t401 + t123 - t323;
t122 = Ifges(7,5) * t128;
t42 = Ifges(7,6) * t152 - Ifges(7,3) * t127 - t122;
t411 = Ifges(5,5) * qJD(4) + qJD(1) * t223 + t174 * t42 + t176 * t395;
t328 = Ifges(5,4) * t175;
t218 = -t177 * Ifges(5,2) - t328;
t410 = t417 * t347 + t401 * t350 + t418 * t351 - Ifges(5,6) * qJD(4) / 0.2e1 - qJD(1) * t218 / 0.2e1;
t356 = t60 / 0.2e1;
t354 = t61 / 0.2e1;
t409 = m(6) + m(7);
t353 = t124 / 0.2e1;
t408 = -t134 / 0.2e1;
t407 = -t135 / 0.2e1;
t406 = g(3) * t175;
t405 = mrSges(3,1) + mrSges(2,1);
t404 = -mrSges(3,3) + mrSges(2,2);
t403 = -mrSges(5,3) + mrSges(4,2);
t400 = (Ifges(6,4) - Ifges(7,5)) * t61 + t419 * t60 + t401 * t124;
t18 = -mrSges(6,1) * t61 + mrSges(6,2) * t60;
t399 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t135 + t18;
t330 = mrSges(6,3) * t127;
t85 = -mrSges(6,2) * t152 + t330;
t332 = mrSges(7,2) * t127;
t88 = mrSges(7,3) * t152 + t332;
t335 = t85 + t88;
t329 = mrSges(6,3) * t128;
t86 = mrSges(6,1) * t152 + t329;
t331 = mrSges(7,2) * t128;
t87 = -mrSges(7,1) * t152 - t331;
t392 = t86 - t87;
t146 = mrSges(5,1) * t177 - mrSges(5,2) * t175;
t391 = t146 + mrSges(4,1);
t207 = pkin(5) * t174 - qJ(6) * t176;
t390 = -t295 - (-qJD(1) * t207 + t90) * t177 + qJD(5) * t207 - qJD(6) * t174;
t388 = -t175 * t417 - t177 * t384;
t387 = -t175 * t401 - t177 * t382;
t386 = t177 * (Ifges(5,2) * t175 - t327) + t175 * (-Ifges(5,1) * t177 + t328);
t385 = t174 * t401 + t176 * t418;
t321 = Ifges(7,5) * t176;
t324 = Ifges(6,4) * t176;
t383 = t174 * t419 - t321 + t324;
t137 = t315 * qJ(2) + t172 * t178;
t130 = -pkin(7) + t137;
t242 = t315 * qJD(2);
t234 = t177 * t242;
t380 = -t130 * t294 + t234;
t377 = t124 * t417 + t401 * t60 + t418 * t61;
t292 = qJD(4) * t177;
t24 = qJDD(3) * t177 - t175 * t415 - t90 * t292;
t376 = -t175 * t24 + t177 * t23;
t374 = Ifges(6,4) * t356 + Ifges(6,6) * t353 - t60 * Ifges(7,5) / 0.2e1 - t124 * Ifges(7,6) / 0.2e1 + (Ifges(6,2) + Ifges(7,3)) * t354;
t17 = -mrSges(7,1) * t61 - mrSges(7,3) * t60;
t373 = t17 + t399;
t372 = pkin(8) * t409;
t72 = -mrSges(7,1) * t127 + mrSges(7,3) * t128;
t371 = t72 + t389;
t370 = m(4) + m(5) + t409;
t228 = mrSges(5,1) * t175 + mrSges(5,2) * t177;
t369 = -t175 * t416 + t177 * t402 - t228;
t363 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t69 = qJD(3) * t177 - t175 * t90;
t63 = -qJD(4) * pkin(4) - t69;
t362 = -m(6) * t63 - t389;
t361 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t232 = -pkin(4) * t175 + pkin(8) * t177;
t296 = qJD(2) * t172;
t107 = qJD(4) * t232 + t296;
t287 = qJD(5) * t177;
t136 = -t172 * qJ(2) + t178 * t315;
t129 = pkin(3) - t136;
t91 = t129 + t233;
t9 = -t174 * (qJD(5) * t91 + t380) - t176 * (t130 * t287 - t107);
t224 = t174 * mrSges(7,1) - t176 * mrSges(7,3);
t226 = mrSges(6,1) * t174 + mrSges(6,2) * t176;
t26 = -pkin(5) * t127 + qJ(6) * t128 + t63;
t326 = Ifges(6,4) * t128;
t45 = Ifges(6,2) * t127 + Ifges(6,6) * t152 - t326;
t360 = t224 * t26 + t226 * t63 - t174 * t45 / 0.2e1;
t359 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t179 = qJD(1) ^ 2;
t355 = -t61 / 0.2e1;
t352 = -t127 / 0.2e1;
t349 = t128 / 0.2e1;
t348 = -t152 / 0.2e1;
t344 = cos(qJ(1));
t343 = sin(qJ(1));
t22 = -qJDD(4) * pkin(4) - t24;
t5 = -pkin(5) * t61 - qJ(6) * t60 + qJD(6) * t128 + t22;
t337 = t175 * t5;
t41 = t130 * t305 + t174 * t91;
t333 = mrSges(4,1) * t172;
t318 = t176 * t91;
t98 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t134;
t316 = t177 * t98;
t133 = t232 * qJD(1);
t36 = t174 * t133 + t176 * t69;
t125 = -t172 * t343 - t315 * t344;
t314 = t125 * t175;
t313 = t125 * t177;
t126 = t172 * t344 - t315 * t343;
t312 = t126 * t175;
t311 = t126 * t177;
t308 = t174 * t175;
t306 = t175 * t176;
t304 = t344 * pkin(1) + t343 * qJ(2);
t303 = mrSges(4,1) * qJDD(1);
t302 = mrSges(4,2) * qJDD(1);
t286 = qJDD(1) * mrSges(3,1);
t279 = t174 * t107 + t176 * t234 + t91 * t288;
t271 = mrSges(5,3) * t297;
t264 = t344 * pkin(2) + t304;
t258 = t177 * t285;
t248 = t290 / 0.2e1;
t247 = t288 / 0.2e1;
t244 = -t283 / 0.2e1;
t238 = -pkin(1) * t343 + t344 * qJ(2);
t235 = t175 * t242;
t217 = -Ifges(6,2) * t174 + t324;
t216 = Ifges(6,2) * t176 + t325;
t213 = -Ifges(5,5) * t177 + Ifges(5,6) * t175;
t210 = Ifges(7,3) * t174 + t321;
t209 = -Ifges(7,3) * t176 + t322;
t35 = t133 * t176 - t174 * t69;
t204 = t130 - t207;
t203 = t175 * t22 + t292 * t63;
t55 = t125 * t174 + t126 * t305;
t54 = -t125 * t176 + t126 * t307;
t202 = t89 * t228;
t201 = -t172 * t95 + t315 * t96;
t200 = -t125 * pkin(3) + pkin(7) * t126 + t264;
t196 = -t174 * t289 + t176 * t292;
t195 = t175 * t288 + t258;
t194 = -pkin(2) * t343 + t238;
t189 = -Ifges(6,6) * t175 - t177 * t217;
t188 = -Ifges(7,6) * t175 - t177 * t210;
t186 = t126 * pkin(3) + t125 * pkin(7) + t194;
t183 = (-t175 * t70 - t177 * t69) * qJD(4) + t376;
t182 = -t175 * t315 * t69 + t172 * t89 + t245 * t70;
t162 = -qJDD(1) * pkin(1) + qJDD(2);
t145 = -qJD(4) * mrSges(5,2) - t271;
t131 = t146 * qJD(1);
t121 = t226 * t175;
t76 = -mrSges(5,1) * t134 + mrSges(5,2) * t135;
t71 = -pkin(5) * t128 - qJ(6) * t127;
t68 = t204 * t175;
t59 = -t125 * t305 + t126 * t174;
t58 = -t125 * t307 - t126 * t176;
t40 = -t130 * t307 + t318;
t38 = -t318 + (t130 * t174 - pkin(5)) * t177;
t37 = qJ(6) * t177 + t41;
t29 = pkin(5) * t299 - t35;
t28 = -qJ(6) * t299 + t36;
t25 = t204 * t292 + (qJD(5) * t208 + qJD(6) * t176 + t242) * t175;
t8 = (-t174 * t287 - t260) * t130 + t279;
t7 = pkin(5) * t294 - t9;
t6 = (-t130 * t290 + qJD(6)) * t177 + (-t130 * t176 - qJ(6)) * t294 + t279;
t10 = [t130 * t316 + t284 * t333 + t131 * t296 + t137 * t302 + pkin(1) * t286 + t389 * (t130 * t292 + t235) - (qJD(5) * t42 + t400) * t306 / 0.2e1 + m(6) * (t130 * t203 + t19 * t9 + t20 * t8 + t235 * t63 + t3 * t41 + t4 * t40) + t135 * t223 / 0.2e1 + t134 * t218 / 0.2e1 + t374 * t308 + m(4) * (qJD(2) * t201 + t82 * t136 + t83 * t137) + (t19 * t196 + t195 * t20 + t3 * t308 + t306 * t4) * mrSges(6,3) + (qJD(1) * t242 + t83) * mrSges(4,2) + t380 * t145 + (Ifges(2,3) + Ifges(4,3) + Ifges(3,2)) * qJDD(1) + (Ifges(5,1) * t407 + Ifges(5,4) * t408 - Ifges(5,5) * qJDD(4) + t130 * t399 - t210 * t355 - t217 * t354 + t45 * t247 + t248 * t395 - t353 * t384 - t356 * t382) * t175 + (t209 * t352 + t216 * t351 + t347 * t385 + t350 * t383) * t289 + t386 * t244 + (t188 * t352 + t189 * t351 - t202 + t213 * qJD(4) / 0.2e1 + t387 * t350 + t388 * t347) * qJD(4) + t45 * t258 / 0.2e1 + t26 * (-mrSges(7,1) * t195 + mrSges(7,3) * t196) + t63 * (-mrSges(6,1) * t195 - mrSges(6,2) * t196) + t8 * t85 + t9 * t86 + t7 * t87 + t6 * t88 + t25 * t72 - t82 * mrSges(4,1) + t68 * t17 + t40 * t31 + t41 * t33 + t37 * t34 + t38 * t32 + (t1 * t308 - t14 * t196 + t15 * t195 - t2 * t306) * mrSges(7,2) - t224 * t337 + (t377 / 0.2e1 + Ifges(6,6) * t354 + Ifges(7,6) * t355 + Ifges(5,4) * t407 + Ifges(5,2) * t408 + t401 * t356 + t417 * t353 + t359 - Ifges(5,6) * qJDD(4)) * t177 + m(3) * (-pkin(1) * t162 + (t149 + t284) * qJ(2)) + m(5) * (qJD(2) * t182 + t77 * t129 + t130 * t183) + (t292 * t69 - t376) * mrSges(5,3) - t136 * t303 - t411 * t292 / 0.2e1 + (-t19 * mrSges(6,1) + t14 * mrSges(7,1) + t20 * mrSges(6,2) + t70 * mrSges(5,3) - t15 * mrSges(7,3) - t410) * t294 + (-m(3) * t238 - m(4) * t194 - m(5) * t186 - t409 * (pkin(4) * t311 + pkin(8) * t312 + t186) - t363 * t55 + t361 * t54 + t404 * t344 + t405 * t343 - t402 * t312 - t391 * t126 + t403 * t125) * g(1) + (-m(3) * t304 - m(4) * t264 - m(5) * t200 - t409 * (-pkin(4) * t313 - pkin(8) * t314 + t200) - t363 * t59 + t361 * t58 - t405 * t344 + t404 * t343 + t402 * t314 + t403 * t126 + t391 * t125) * g(2) - t22 * t121 + t129 * t76 + t77 * t146 + 0.2e1 * t149 * mrSges(3,3) - t162 * mrSges(3,1) + m(7) * (t1 * t37 + t14 * t7 + t15 * t6 + t2 * t38 + t25 * t26 + t5 * t68); m(3) * t162 - t131 * t300 - t286 + (-t76 - t303) * t315 + t420 * t145 + t397 * t112 + t398 * t111 + (t1 * t112 + t111 * t2 + t393 * t14 + t394 * t15) * m(7) + (-t111 * t4 + t112 * t3 - t19 * t393 + t20 * t394) * m(6) + (-qJD(1) * t182 - t315 * t77) * m(5) + (-qJD(1) * t201 + t315 * t82) * m(4) + (-m(7) * t26 + t362 - t72) * t175 * t243 + (-m(3) * qJ(2) - mrSges(4,2) * t315 - mrSges(3,3) - t333) * t179 - t393 * t392 + (-t343 * g(1) + t344 * g(2)) * (m(3) + t370) + t394 * t335 + (t302 + t316 + (t26 * t292 + t337) * m(7) + t203 * m(6) + t183 * m(5) + t83 * m(4) + t373 * t175 + t371 * t292) * t172; m(4) * qJDD(3) + t370 * g(3) + ((-t174 * t392 + t176 * t335 + t145) * qJD(4) + m(5) * (t70 * qJD(4) + t24) + m(7) * (t14 * t285 + t15 * t293 - t5) + m(6) * (-t19 * t285 + t20 * t293 - t22) - t373) * t177 + (t98 + t371 * qJD(4) + (-t174 * t335 - t176 * t392) * qJD(5) + m(5) * (-qJD(4) * t69 + t23) + m(7) * (qJD(4) * t26 + t413) + m(6) * (qJD(4) * t63 + t414) + t412) * t175; t216 * t354 + t209 * t355 + t390 * t72 + t395 * t247 + t400 * t174 / 0.2e1 + t5 * t225 + t22 * t227 + (-pkin(4) * t22 - t19 * t35 - t20 * t36) * m(6) + (t175 * t372 + t177 * t416 + t146) * g(3) + (t360 + t411 / 0.2e1) * t297 + (-t272 + t362) * t70 + t374 * t176 + (t125 * t369 + t313 * t372) * g(1) + (t126 * t369 + t311 * t372) * g(2) + (t217 / 0.2e1 - t210 / 0.2e1) * t291 + t383 * t356 + t385 * t353 + t386 * t179 / 0.2e1 + ((-t189 / 0.2e1 + t188 / 0.2e1) * t127 + t202 - t14 * (mrSges(7,1) * t175 - mrSges(7,2) * t305) - t19 * (-mrSges(6,1) * t175 + mrSges(6,3) * t305) - t20 * (mrSges(6,2) * t175 + mrSges(6,3) * t307) - t15 * (mrSges(7,2) * t307 - mrSges(7,3) * t175) + t387 * t349 + t388 * t348) * qJD(1) - t36 * t85 - t35 * t86 - t29 * t87 - t28 * t88 + (-t14 * t29 - t15 * t28 - t205 * t5 + t390 * t26) * m(7) - t205 * t17 - pkin(4) * t18 - t23 * mrSges(5,2) + t24 * mrSges(5,1) + (-t271 - t145) * t69 + (t347 * t384 + t350 * t382 + t360) * qJD(5) + (((t14 * t176 - t15 * t174) * qJD(5) + t231) * m(7) - t392 * t288 - t335 * t290 + ((-t20 * t174 - t19 * t176) * qJD(5) + t230) * m(6) + t412) * pkin(8) + (t406 + t413) * mrSges(7,2) + (t406 + t414) * mrSges(6,3) + t410 * t299 + Ifges(5,6) * t134 + Ifges(5,5) * t135 + t213 * t244 + t42 * t248 + Ifges(5,3) * qJDD(4); t45 * t350 + (-Ifges(7,3) * t128 + t323) * t351 + t377 + (t401 * t127 + t128 * t418) * t348 + (t361 * t59 + t363 * t58) * g(1) + (-t361 * t55 - t363 * t54) * g(2) - t14 * t332 - t15 * t331 + qJD(6) * t88 - t71 * t72 - pkin(5) * t32 + qJ(6) * t34 + (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t15 - t26 * t71) * m(7) + t359 + (t127 * t419 - t122 + t326 + t42) * t349 + (-m(7) * t14 - t329 + t392) * t20 + (-m(7) * t15 + t330 - t335) * t19 + (Ifges(6,2) * t128 + t123 + t395) * t352 + (-t121 - (m(7) * t207 + t224) * t175) * g(3) - t63 * (-mrSges(6,1) * t128 + mrSges(6,2) * t127) - t26 * (-mrSges(7,1) * t128 - mrSges(7,3) * t127); -t128 * t72 - t152 * t88 + (-g(1) * t58 + g(2) * t54 + g(3) * t308 - t26 * t128 - t15 * t152 + t2) * m(7) + t32;];
tau  = t10;
