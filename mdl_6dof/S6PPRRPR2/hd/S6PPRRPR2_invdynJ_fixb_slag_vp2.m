% Calculate vector of inverse dynamics joint torques for
% S6PPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
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
% Datum: 2019-03-08 18:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PPRRPR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRPR2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRPR2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:49:08
% EndTime: 2019-03-08 18:49:29
% DurationCPUTime: 11.79s
% Computational Cost: add. (4847->550), mult. (12192->773), div. (0->0), fcn. (10568->14), ass. (0->273)
t396 = -mrSges(5,1) + mrSges(6,2);
t387 = m(7) + m(6);
t385 = Ifges(5,5) - Ifges(6,4);
t384 = Ifges(6,5) - Ifges(5,6);
t177 = sin(qJ(4));
t180 = cos(qJ(4));
t217 = pkin(10) * t177 - qJ(5) * t180;
t284 = qJD(5) * t177;
t287 = qJD(4) * t177;
t242 = pkin(4) * t287 - t284;
t308 = cos(pkin(6));
t162 = qJD(1) * t308 + qJD(2);
t172 = sin(pkin(12));
t174 = sin(pkin(6));
t178 = sin(qJ(3));
t175 = cos(pkin(12));
t307 = cos(pkin(7));
t247 = t175 * t307;
t236 = t178 * t247;
t337 = cos(qJ(3));
t192 = (t172 * t337 + t236) * t174;
t173 = sin(pkin(7));
t300 = t173 * t178;
t71 = qJD(1) * t192 + t162 * t300;
t395 = qJD(4) * t217 + t242 - t71;
t176 = sin(qJ(6));
t179 = cos(qJ(6));
t288 = qJD(3) * t180;
t134 = -qJD(4) * t176 - t179 * t288;
t290 = qJD(3) * t177;
t163 = qJD(6) + t290;
t102 = -mrSges(7,2) * t163 + mrSges(7,3) * t134;
t286 = qJD(4) * t179;
t205 = t176 * t288 - t286;
t103 = mrSges(7,1) * t163 + mrSges(7,3) * t205;
t215 = t179 * t102 - t176 * t103;
t279 = qJD(3) * qJD(4);
t146 = qJDD(3) * t177 + t180 * t279;
t133 = qJDD(6) + t146;
t145 = -t180 * qJDD(3) + t177 * t279;
t82 = qJD(6) * t134 + qJDD(4) * t179 + t145 * t176;
t65 = mrSges(7,1) * t133 - mrSges(7,3) * t82;
t83 = qJD(6) * t205 - qJDD(4) * t176 + t145 * t179;
t66 = -mrSges(7,2) * t133 + mrSges(7,3) * t83;
t394 = t215 * qJD(6) + t176 * t66 + t179 * t65;
t160 = qJDD(1) * t308 + qJDD(2);
t212 = t337 * t247;
t202 = t174 * t212;
t214 = t174 * t236;
t302 = t172 * t174;
t261 = qJD(1) * t302;
t238 = qJD(3) * t261;
t258 = qJDD(1) * t302;
t289 = qJD(3) * t178;
t260 = t173 * t289;
t262 = t173 * t337;
t37 = -qJD(3) * qJD(1) * t214 + qJDD(1) * t202 + t160 * t262 - t162 * t260 - t178 * t258 - t337 * t238;
t34 = -qJDD(3) * pkin(3) - t37;
t393 = m(5) * t34;
t343 = pkin(4) + pkin(10);
t267 = t174 * t175 * t173;
t108 = -qJD(1) * t267 + t162 * t307;
t101 = t180 * t108;
t69 = qJD(3) * pkin(9) + t71;
t39 = -t177 * (pkin(5) * qJD(3) + t69) + t101;
t376 = qJD(5) - t39;
t35 = -qJD(4) * t343 + t376;
t296 = t177 * qJ(5);
t252 = -pkin(3) - t296;
t129 = -t180 * t343 + t252;
t197 = qJD(1) * t202;
t186 = -t162 * t262 + t178 * t261 - t197;
t54 = qJD(3) * t129 + t186;
t11 = -t176 * t54 + t179 * t35;
t183 = -t146 * qJ(5) - qJD(3) * t284 + t34;
t13 = t145 * t343 + t183;
t106 = -qJDD(1) * t267 + t160 * t307;
t285 = qJD(4) * t180;
t239 = qJD(3) * t262;
t36 = qJD(3) * t197 + qJDD(1) * t214 + t160 * t300 + t162 * t239 - t178 * t238 + t337 * t258;
t33 = qJDD(3) * pkin(9) + t36;
t10 = t106 * t180 - t108 * t287 - t177 * t33 - t69 * t285;
t203 = qJDD(5) - t10;
t5 = pkin(5) * t146 - qJDD(4) * t343 + t203;
t1 = qJD(6) * t11 + t13 * t179 + t176 * t5;
t392 = t1 * mrSges(7,2);
t12 = t176 * t35 + t179 * t54;
t2 = -qJD(6) * t12 - t13 * t176 + t179 * t5;
t391 = t2 * mrSges(7,1);
t270 = mrSges(5,3) * t290;
t272 = mrSges(6,1) * t290;
t390 = qJD(4) * t396 + t270 + t272;
t389 = m(7) * pkin(10) + pkin(4) * t387 - t396;
t49 = t108 * t177 + t180 * t69;
t45 = -qJD(4) * qJ(5) - t49;
t386 = m(6) * t45;
t342 = pkin(5) + pkin(9);
t158 = t342 * t180;
t144 = qJD(4) * t158;
t295 = t177 * t179;
t157 = t342 * t177;
t90 = t129 * t179 + t157 * t176;
t383 = -qJD(6) * t90 + t144 * t179 - t176 * t395 + t186 * t295;
t298 = t176 * t177;
t89 = -t129 * t176 + t157 * t179;
t382 = qJD(6) * t89 + t144 * t176 + t179 * t395 + t186 * t298;
t306 = cos(pkin(11));
t234 = t308 * t306;
t305 = sin(pkin(11));
t116 = t172 * t234 + t175 * t305;
t190 = t172 * t305 - t175 * t234;
t249 = t174 * t306;
t371 = t173 * t249 + t190 * t307;
t58 = t116 * t178 + t337 * t371;
t313 = t180 * t58;
t381 = -pkin(4) * t313 - t58 * t296;
t233 = t308 * t305;
t117 = -t172 * t233 + t175 * t306;
t191 = t172 * t306 + t175 * t233;
t248 = t174 * t305;
t370 = -t173 * t248 + t191 * t307;
t60 = t117 * t178 + t337 * t370;
t312 = t180 * t60;
t380 = -pkin(4) * t312 - t60 * t296;
t250 = t173 * t308;
t213 = t337 * t250;
t301 = t172 * t178;
t87 = t174 * t301 - t202 - t213;
t310 = t180 * t87;
t379 = -pkin(4) * t310 - t87 * t296;
t165 = Ifges(5,4) * t288;
t378 = Ifges(5,1) * t290 + Ifges(5,5) * qJD(4) - Ifges(7,5) * t205 + t134 * Ifges(7,6) + t163 * Ifges(7,3) + t165;
t271 = mrSges(6,1) * t288;
t154 = -qJD(4) * mrSges(6,3) - t271;
t91 = -mrSges(7,1) * t134 - mrSges(7,2) * t205;
t377 = t91 - t154;
t114 = t146 * mrSges(6,1) + qJDD(4) * mrSges(6,2);
t375 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t146 + t114;
t113 = mrSges(6,1) * t145 - qJDD(4) * mrSges(6,3);
t374 = -qJDD(4) * mrSges(5,2) - mrSges(5,3) * t145 - t113;
t229 = mrSges(6,2) * t180 - mrSges(6,3) * t177;
t136 = t229 * qJD(3);
t232 = -mrSges(5,1) * t180 + mrSges(5,2) * t177;
t373 = t232 * qJD(3) + t136;
t269 = mrSges(5,3) * t288;
t372 = -qJD(4) * mrSges(5,2) - t154 + t269;
t322 = Ifges(6,6) * t180;
t323 = Ifges(6,6) * t177;
t369 = t177 * (-Ifges(6,2) * t180 + t323) + t180 * (Ifges(6,3) * t177 - t322);
t368 = t177 * t384 + t180 * t385;
t118 = t177 * t300 - t180 * t307;
t276 = t177 * t106 + t108 * t285 + t180 * t33;
t9 = -t287 * t69 + t276;
t365 = -t10 * t177 + t180 * t9;
t364 = t1 * t176 + t179 * t2;
t317 = t177 * t69;
t7 = -qJDD(4) * qJ(5) + qJD(4) * (-qJD(5) + t317) - t276;
t8 = -qJDD(4) * pkin(4) + t203;
t363 = t177 * t8 - t180 * t7;
t96 = mrSges(5,1) * t145 + mrSges(5,2) * t146;
t97 = -mrSges(6,2) * t145 - mrSges(6,3) * t146;
t361 = mrSges(4,1) * qJDD(3) - t96 - t97;
t218 = -t180 * Ifges(6,3) - t323;
t320 = t205 * Ifges(7,4);
t73 = t134 * Ifges(7,2) + t163 * Ifges(7,6) - t320;
t130 = Ifges(7,4) * t134;
t74 = -Ifges(7,1) * t205 + Ifges(7,5) * t163 + t130;
t360 = Ifges(6,5) * qJD(4) + qJD(3) * t218 + t176 * t74 + t179 * t73;
t359 = m(5) + m(4) + m(3) + t387;
t48 = -t101 + t317;
t41 = -qJD(4) * pkin(4) + qJD(5) + t48;
t357 = -m(6) * t41 - t390;
t231 = t179 * mrSges(7,1) - t176 * mrSges(7,2);
t230 = mrSges(7,1) * t176 + mrSges(7,2) * t179;
t356 = -qJ(5) * t387 + mrSges(5,2) - mrSges(6,3) - t230;
t167 = pkin(5) * t288;
t38 = t167 - t45;
t355 = -m(7) * t38 - t372 - t91;
t354 = -m(5) * t10 + m(6) * t8 + t375;
t353 = m(5) * t48 - t357;
t184 = t173 * t190 - t249 * t307;
t59 = t116 * t337 - t178 * t371;
t21 = t177 * t59 - t180 * t184;
t185 = t173 * t191 + t248 * t307;
t61 = t117 * t337 - t178 * t370;
t23 = t177 * t61 - t180 * t185;
t196 = t307 * t308 - t267;
t88 = t178 * t250 + t192;
t62 = t177 * t88 - t180 * t196;
t204 = -g(1) * t23 - g(2) * t21 - g(3) * t62;
t352 = -g(1) * t61 - g(2) * t59 - g(3) * t88;
t150 = -pkin(4) * t180 + t252;
t64 = qJD(3) * t150 + t186;
t68 = -qJD(3) * pkin(3) + t186;
t351 = -m(5) * t68 - m(6) * t64 + mrSges(4,1) * qJD(3) - t373;
t350 = mrSges(7,1) * t298 + mrSges(7,2) * t295 + mrSges(4,1) - t229 - t232;
t349 = -m(7) * t342 + mrSges(4,2) - t231;
t42 = -mrSges(7,1) * t83 + mrSges(7,2) * t82;
t6 = -pkin(5) * t145 - t7;
t348 = m(5) * t9 - m(6) * t7 + m(7) * t6 + t374 + t42;
t347 = m(5) * t49 - t355 - t386;
t182 = qJD(3) ^ 2;
t346 = -t82 * Ifges(7,4) / 0.2e1 - t83 * Ifges(7,2) / 0.2e1 - t133 * Ifges(7,6) / 0.2e1;
t345 = t82 / 0.2e1;
t344 = t83 / 0.2e1;
t341 = t133 / 0.2e1;
t339 = -t205 / 0.2e1;
t327 = Ifges(5,4) * t177;
t326 = Ifges(5,4) * t180;
t325 = Ifges(7,4) * t176;
t324 = Ifges(7,4) * t179;
t303 = qJD(3) * mrSges(4,2);
t297 = t176 * t180;
t293 = t179 * t180;
t291 = mrSges(4,2) * qJDD(3);
t283 = qJD(6) * t176;
t282 = qJD(6) * t179;
t281 = qJD(6) * t180;
t275 = Ifges(7,5) * t82 + Ifges(7,6) * t83 + Ifges(7,3) * t133;
t55 = t58 * pkin(3);
t266 = pkin(9) * t59 - t55;
t56 = t60 * pkin(3);
t265 = pkin(9) * t61 - t56;
t86 = t87 * pkin(3);
t264 = pkin(9) * t88 - t86;
t259 = t176 * t281;
t253 = -t282 / 0.2e1;
t245 = -t279 / 0.2e1;
t228 = Ifges(7,1) * t179 - t325;
t227 = Ifges(7,1) * t176 + t324;
t226 = t180 * Ifges(5,2) + t327;
t224 = -Ifges(7,2) * t176 + t324;
t223 = Ifges(7,2) * t179 + t325;
t221 = Ifges(7,5) * t179 - Ifges(7,6) * t176;
t220 = Ifges(7,5) * t176 + Ifges(7,6) * t179;
t219 = -t177 * Ifges(6,2) - t322;
t216 = t11 * t176 - t12 * t179;
t25 = -t176 * t87 + t179 * t62;
t26 = t176 * t62 + t179 * t87;
t210 = t64 * (-mrSges(6,2) * t177 - mrSges(6,3) * t180);
t209 = t68 * (mrSges(5,1) * t177 + mrSges(5,2) * t180);
t208 = t177 * (Ifges(5,1) * t180 - t327);
t94 = t179 * t118 + t176 * t262;
t200 = -t176 * t118 + t179 * t262;
t199 = t177 * t286 + t259;
t198 = t176 * t287 - t179 * t281;
t119 = t177 * t307 + t180 * t300;
t195 = Ifges(7,5) * t180 + t177 * t227;
t194 = Ifges(7,6) * t180 + t177 * t223;
t193 = Ifges(7,3) * t180 + t177 * t220;
t189 = -qJD(6) * t216 + t364;
t63 = t177 * t196 + t88 * t180;
t166 = pkin(4) * t290;
t143 = t342 * t287;
t138 = -qJ(5) * t288 + t166;
t123 = Ifges(6,4) * qJD(4) + qJD(3) * t219;
t120 = Ifges(5,6) * qJD(4) + qJD(3) * t226;
t115 = -qJ(5) * t285 + t242;
t110 = qJD(3) * t217 + t166;
t93 = qJD(4) * t119 + t177 * t239;
t80 = t88 * qJD(3);
t79 = (t213 + (t212 - t301) * t174) * qJD(3);
t47 = qJD(6) * t200 - t176 * t260 + t179 * t93;
t46 = qJD(6) * t94 + t176 * t93 + t179 * t260;
t40 = t167 + t49;
t28 = t82 * Ifges(7,1) + t83 * Ifges(7,4) + t133 * Ifges(7,5);
t17 = qJD(4) * t63 + t79 * t177;
t16 = t110 * t179 + t176 * t40;
t15 = -t110 * t176 + t179 * t40;
t14 = t145 * pkin(4) + t183;
t4 = qJD(6) * t25 + t17 * t176 + t179 * t80;
t3 = -qJD(6) * t26 + t17 * t179 - t176 * t80;
t18 = [-t88 * t291 + m(4) * (t106 * t196 + t36 * t88 + t71 * t79) + m(3) * (t160 * t308 + (t172 ^ 2 + t175 ^ 2) * t174 ^ 2 * qJDD(1)) - t79 * t303 + m(7) * (t1 * t26 + t11 * t3 + t12 * t4 + t2 * t25) + t25 * t65 + t26 * t66 + t4 * t102 + t3 * t103 + m(2) * qJDD(1) + t354 * t62 + t353 * t17 + (-m(4) * t37 + m(6) * t14 - t361 + t393) * t87 + (m(4) * t186 - t351) * t80 + t348 * t63 + t347 * (-t88 * t287 + (qJD(4) * t196 + t79) * t180) + (-m(2) - t359) * g(3); m(4) * t106 * t307 + m(7) * (-t1 * t200 + t11 * t47 + t12 * t46 + t2 * t94) + m(3) * t160 + t94 * t65 - t200 * t66 + t46 * t102 + t47 * t103 + (m(5) * (t289 * t68 - t337 * t34) + m(6) * (-t14 * t337 + t289 * t64) + m(4) * (t337 * t37 + t178 * t36 + (t178 * t186 + t337 * t71) * qJD(3))) * t173 + t353 * t93 + (-mrSges(4,1) * t182 - t291) * t300 + t373 * t260 + t354 * t118 - t347 * (qJD(4) * t118 - t180 * t239) + t348 * t119 + (-mrSges(4,2) * t182 + t361) * t262 + (-g(1) * t248 + g(2) * t249 - g(3) * t308) * t359; t390 * (pkin(9) * t285 + t177 * t186) + t6 * t231 * t180 + (t180 * (-Ifges(5,2) * t177 + t326) + t208) * t279 / 0.2e1 + (-t123 / 0.2e1 - t12 * mrSges(7,2) + t11 * mrSges(7,1) + t378 / 0.2e1 + t48 * mrSges(5,3) + t41 * mrSges(6,1)) * t285 + t134 * (qJD(4) * t194 - t224 * t281) / 0.2e1 + t163 * (qJD(4) * t193 - t221 * t281) / 0.2e1 + t146 * t326 / 0.2e1 + (t177 * t385 - t180 * t384) * qJDD(4) / 0.2e1 + t382 * t102 + (t1 * t90 + t11 * t383 + t12 * t382 - t143 * t38 + t158 * t6 + t2 * t89) * m(7) + t383 * t103 + t368 * qJD(4) ^ 2 / 0.2e1 + t369 * t245 + (-t120 / 0.2e1 - t49 * mrSges(5,3) + t45 * mrSges(6,1) + t360 / 0.2e1) * t287 - t28 * t297 / 0.2e1 + (m(6) * (t177 * t41 - t180 * t45) + m(5) * (t177 * t48 + t180 * t49) - t303 - t355 * t180) * t186 + m(6) * (t115 * t64 + t14 * t150) + t38 * (-mrSges(7,1) * t199 + mrSges(7,2) * t198) + t177 * t391 + t351 * t71 + (m(5) * ((-t177 * t49 + t180 * t48) * qJD(4) + t365) + m(6) * ((t177 * t45 + t180 * t41) * qJD(4) + t363) - t372 * t287 + t374 * t180 + t375 * t177) * pkin(9) + t180 * t74 * t253 + t14 * t229 + t34 * t232 + t145 * t218 / 0.2e1 - t146 * t219 / 0.2e1 - t145 * t226 / 0.2e1 + Ifges(4,3) * qJDD(3) + t73 * t259 / 0.2e1 + (-Ifges(5,4) * t145 + Ifges(5,5) * qJDD(4) + t275) * t177 / 0.2e1 + (-m(6) * (t264 + t379) - m(7) * (-pkin(10) * t310 + t379 - t86) - m(5) * t264 + t349 * t88 + t350 * t87) * g(3) + (-m(6) * (t265 + t380) - m(7) * (-pkin(10) * t312 + t380 - t56) - m(5) * t265 + t349 * t61 + t350 * t60) * g(1) + (-m(6) * (t266 + t381) - m(7) * (-pkin(10) * t313 + t381 - t55) - m(5) * t266 + t349 * t59 + t350 * t58) * g(2) + (g(1) * t312 + g(2) * t313 + g(3) * t310 - t1 * t293 - t11 * t198 + t12 * t199 + t2 * t297) * mrSges(7,3) + t180 * (Ifges(5,4) * t146 - Ifges(5,2) * t145 + Ifges(5,6) * qJDD(4)) / 0.2e1 - t180 * (Ifges(6,5) * qJDD(4) - Ifges(6,6) * t146 + Ifges(6,3) * t145) / 0.2e1 - t177 * (Ifges(6,4) * qJDD(4) - Ifges(6,2) * t146 + Ifges(6,6) * t145) / 0.2e1 - t177 * t392 + (-t96 - t393) * pkin(3) + (t352 + t365) * mrSges(5,3) + (t352 + t363) * mrSges(6,1) - t36 * mrSges(4,2) + t37 * mrSges(4,1) + t89 * t65 + t90 * t66 + t115 * t136 + t146 * t177 * Ifges(5,1) - t143 * t91 + t150 * t97 + t158 * t42 + qJD(4) * t209 + qJD(4) * t210 + (qJD(4) * t195 - t228 * t281) * t339 + (Ifges(7,3) * t177 - t180 * t220) * t341 + (Ifges(7,6) * t177 - t180 * t223) * t344 + (Ifges(7,5) * t177 - t180 * t227) * t345 + t293 * t346; (t270 + t357) * t49 + t120 * t290 / 0.2e1 + (-t209 - t210 - t12 * (-mrSges(7,2) * t180 + mrSges(7,3) * t295) - t11 * (mrSges(7,1) * t180 - mrSges(7,3) * t298)) * qJD(3) + (Ifges(5,3) + Ifges(6,1)) * qJDD(4) + (-m(7) * t189 - t394) * t343 + (-pkin(4) * t8 - qJ(5) * t7 - qJD(5) * t45 - t138 * t64) * m(6) + t123 * t288 / 0.2e1 - t74 * t283 / 0.2e1 + (-t269 + t372 - t386) * t48 + t384 * t145 + t385 * t146 + t377 * qJD(5) - (-Ifges(5,2) * t290 + t165 + t378) * t288 / 0.2e1 + t368 * t245 + (t11 * t283 - t12 * t282 - t204 - t364) * mrSges(7,3) - t360 * t290 / 0.2e1 + (-t208 / 0.2e1 + t369 / 0.2e1) * t182 - (t134 * t223 + t163 * t220 - t205 * t227) * qJD(6) / 0.2e1 - (t134 * t194 + t163 * t193 - t195 * t205) * qJD(3) / 0.2e1 + (qJ(5) * t6 - t11 * t15 - t12 * t16 + t376 * t38) * m(7) + t163 * t38 * t231 + (t356 * (t177 * t185 + t61 * t180) + t389 * t23) * g(1) + (t356 * (t177 * t184 + t59 * t180) + t389 * t21) * g(2) + (t356 * t63 + t389 * t62) * g(3) + t6 * t230 - t41 * t271 - t45 * t272 + t179 * t28 / 0.2e1 - t7 * mrSges(6,3) + t8 * mrSges(6,2) - t9 * mrSges(5,2) + t10 * mrSges(5,1) + (t42 - t113) * qJ(5) + t73 * t253 - t39 * t91 - t16 * t102 - t15 * t103 - pkin(4) * t114 - t138 * t136 + t221 * t341 + t224 * t344 + t228 * t345 + t176 * t346; -t377 * qJD(4) + (t136 + t215) * t290 + t114 + (-qJD(4) * t38 - t216 * t290 + t189 + t204) * m(7) + (qJD(4) * t45 + t290 * t64 + t204 + t8) * m(6) + t394; -t392 + t391 - t38 * (-mrSges(7,1) * t205 + mrSges(7,2) * t134) + t205 * (Ifges(7,1) * t134 + t320) / 0.2e1 + t73 * t339 - t163 * (Ifges(7,5) * t134 + Ifges(7,6) * t205) / 0.2e1 - t11 * t102 + t12 * t103 - g(1) * ((-t176 * t60 + t179 * t23) * mrSges(7,1) + (-t176 * t23 - t179 * t60) * mrSges(7,2)) - g(2) * ((-t176 * t58 + t179 * t21) * mrSges(7,1) + (-t176 * t21 - t179 * t58) * mrSges(7,2)) - g(3) * (mrSges(7,1) * t25 - mrSges(7,2) * t26) + (t11 * t134 - t12 * t205) * mrSges(7,3) + t275 - (Ifges(7,2) * t205 + t130 + t74) * t134 / 0.2e1;];
tau  = t18;
