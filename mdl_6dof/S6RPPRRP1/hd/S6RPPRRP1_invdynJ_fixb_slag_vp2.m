% Calculate vector of inverse dynamics joint torques for
% S6RPPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-03-09 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:57:37
% EndTime: 2019-03-09 01:57:56
% DurationCPUTime: 12.64s
% Computational Cost: add. (6224->523), mult. (13787->646), div. (0->0), fcn. (9717->14), ass. (0->233)
t373 = mrSges(6,3) + mrSges(7,3);
t372 = Ifges(6,4) + Ifges(7,4);
t166 = sin(pkin(10));
t168 = cos(pkin(10));
t173 = sin(qJ(4));
t176 = cos(qJ(4));
t140 = t166 * t176 + t168 * t173;
t132 = t140 * qJD(1);
t172 = sin(qJ(5));
t175 = cos(qJ(5));
t111 = qJD(4) * t175 - t132 * t172;
t139 = t166 * t173 - t176 * t168;
t133 = t139 * qJD(4);
t96 = -qJD(1) * t133 + qJDD(1) * t140;
t55 = qJD(5) * t111 + qJDD(4) * t172 + t175 * t96;
t312 = t55 / 0.2e1;
t112 = qJD(4) * t172 + t132 * t175;
t56 = -qJD(5) * t112 + qJDD(4) * t175 - t172 * t96;
t311 = t56 / 0.2e1;
t134 = t140 * qJD(4);
t97 = -qJD(1) * t134 - qJDD(1) * t139;
t93 = qJDD(5) - t97;
t310 = t93 / 0.2e1;
t357 = Ifges(6,1) + Ifges(7,1);
t355 = Ifges(6,5) + Ifges(7,5);
t354 = Ifges(6,2) + Ifges(7,2);
t353 = Ifges(6,6) + Ifges(7,6);
t352 = Ifges(7,3) + Ifges(6,3);
t152 = pkin(5) * t175 + pkin(4);
t202 = -mrSges(7,1) * t175 + mrSges(7,2) * t172;
t204 = -mrSges(6,1) * t175 + mrSges(6,2) * t172;
t371 = -m(6) * pkin(4) - m(7) * t152 + t202 + t204;
t170 = -qJ(6) - pkin(8);
t370 = -m(6) * pkin(8) + m(7) * t170 - t373;
t167 = sin(pkin(9));
t147 = pkin(1) * t167 + qJ(3);
t137 = qJD(1) * qJD(3) + qJDD(1) * t147;
t369 = t355 * t310 + t372 * t311 + t357 * t312;
t368 = t372 * t111;
t367 = t166 ^ 2 + t168 ^ 2;
t154 = t168 * qJDD(2);
t104 = t154 + (-pkin(7) * qJDD(1) - t137) * t166;
t114 = t166 * qJDD(2) + t168 * t137;
t243 = qJDD(1) * t168;
t105 = pkin(7) * t243 + t114;
t144 = t147 * qJD(1);
t156 = t168 * qJD(2);
t268 = pkin(7) * qJD(1);
t109 = t156 + (-t144 - t268) * t166;
t122 = t166 * qJD(2) + t168 * t144;
t110 = t168 * t268 + t122;
t248 = qJD(4) * t176;
t249 = qJD(4) * t173;
t23 = t173 * t104 + t176 * t105 + t109 * t248 - t110 * t249;
t21 = qJDD(4) * pkin(8) + t23;
t246 = qJD(5) * t175;
t247 = qJD(5) * t172;
t150 = pkin(3) * t168 + pkin(2);
t169 = cos(pkin(9));
t296 = pkin(1) * t169;
t143 = -t150 - t296;
t125 = qJDD(1) * t143 + qJDD(3);
t39 = -pkin(4) * t97 - pkin(8) * t96 + t125;
t62 = t109 * t173 + t110 * t176;
t59 = qJD(4) * pkin(8) + t62;
t128 = qJD(1) * t143 + qJD(3);
t131 = t139 * qJD(1);
t70 = pkin(4) * t131 - pkin(8) * t132 + t128;
t3 = t172 * t39 + t175 * t21 + t70 * t246 - t247 * t59;
t27 = t172 * t70 + t175 * t59;
t4 = -qJD(5) * t27 - t172 * t21 + t175 * t39;
t210 = -t172 * t4 + t175 * t3;
t26 = -t172 * t59 + t175 * t70;
t366 = -t26 * t246 - t27 * t247 + t210;
t365 = t372 * t112;
t364 = t372 * t175;
t363 = t372 * t172;
t165 = qJ(1) + pkin(9);
t158 = sin(t165);
t160 = cos(t165);
t362 = g(1) * t160 + g(2) * t158;
t361 = qJD(5) + t131;
t313 = m(7) * pkin(5);
t360 = -m(4) - m(3);
t359 = -mrSges(7,1) - mrSges(6,1);
t358 = mrSges(6,2) + mrSges(7,2);
t351 = t353 * t93 + t354 * t56 + t372 * t55;
t18 = -mrSges(6,1) * t56 + mrSges(6,2) * t55;
t349 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t96 + t18;
t348 = t353 * t111 + t355 * t112 + t352 * t361;
t347 = t354 * t111 + t353 * t361 + t365;
t346 = t357 * t112 + t355 * t361 + t368;
t345 = Ifges(5,4) * t131;
t344 = t131 * Ifges(5,2);
t278 = mrSges(5,3) * t132;
t343 = qJD(4) * mrSges(5,1) + mrSges(6,1) * t111 - mrSges(6,2) * t112 - t278;
t214 = qJD(5) * t170;
t259 = t131 * t175;
t61 = t109 * t176 - t173 * t110;
t94 = pkin(4) * t132 + pkin(8) * t131;
t34 = -t172 * t61 + t175 * t94;
t342 = -pkin(5) * t132 - qJ(6) * t259 - qJD(6) * t172 + t175 * t214 - t34;
t260 = t131 * t172;
t341 = -t62 + (t247 + t260) * pkin(5);
t35 = t172 * t94 + t175 * t61;
t340 = -qJ(6) * t260 + qJD(6) * t175 + t172 * t214 - t35;
t339 = t313 + mrSges(7,1);
t338 = Ifges(5,5) * qJD(4);
t337 = Ifges(5,6) * qJD(4);
t285 = pkin(7) + t147;
t135 = t285 * t166;
t136 = t285 * t168;
t336 = -t176 * t135 - t136 * t173;
t280 = mrSges(7,2) * t175;
t201 = mrSges(7,1) * t172 + t280;
t203 = mrSges(6,1) * t172 + mrSges(6,2) * t175;
t58 = -qJD(4) * pkin(4) - t61;
t40 = -pkin(5) * t111 + qJD(6) + t58;
t335 = t40 * t201 + t58 * t203;
t334 = -t353 * t172 + t355 * t175;
t333 = -t354 * t172 + t364;
t332 = t357 * t175 - t363;
t329 = t352 * t93 + t353 * t56 + t355 * t55;
t113 = -t137 * t166 + t154;
t328 = -t113 * t166 + t114 * t168;
t325 = -m(5) - m(7) - m(6);
t322 = mrSges(6,1) + t339;
t1 = pkin(5) * t93 - qJ(6) * t55 - qJD(6) * t112 + t4;
t15 = -qJ(6) * t112 + t26;
t14 = pkin(5) * t361 + t15;
t16 = qJ(6) * t111 + t27;
t2 = qJ(6) * t56 + qJD(6) * t111 + t3;
t320 = -t1 * t172 - t14 * t246 - t16 * t247 + t175 * t2;
t319 = -m(6) * t58 + t343;
t164 = pkin(10) + qJ(4);
t157 = sin(t164);
t159 = cos(t164);
t206 = mrSges(5,1) * t159 - mrSges(5,2) * t157;
t207 = -mrSges(4,1) * t168 + mrSges(4,2) * t166;
t318 = m(4) * pkin(2) + t373 * t157 + mrSges(3,1) + t206 - t207;
t317 = -m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t316 = t128 * mrSges(5,2) + t338 / 0.2e1;
t315 = t4 * mrSges(6,1) - t3 * mrSges(6,2) - t2 * mrSges(7,2);
t314 = t16 * mrSges(7,2) + t27 * mrSges(6,2) + t337 / 0.2e1 - t128 * mrSges(5,1) - t14 * mrSges(7,1) - t26 * mrSges(6,1);
t309 = -t111 / 0.2e1;
t308 = t111 / 0.2e1;
t307 = -t112 / 0.2e1;
t306 = t112 / 0.2e1;
t305 = -t361 / 0.2e1;
t304 = t361 / 0.2e1;
t303 = t131 / 0.2e1;
t302 = -t132 / 0.2e1;
t301 = t132 / 0.2e1;
t174 = sin(qJ(1));
t295 = pkin(1) * t174;
t177 = cos(qJ(1));
t161 = t177 * pkin(1);
t29 = mrSges(7,1) * t93 - mrSges(7,3) * t55;
t30 = mrSges(6,1) * t93 - mrSges(6,3) * t55;
t284 = -t29 - t30;
t31 = -mrSges(7,2) * t93 + mrSges(7,3) * t56;
t32 = -mrSges(6,2) * t93 + mrSges(6,3) * t56;
t283 = t31 + t32;
t275 = mrSges(7,3) * t111;
t71 = -mrSges(7,2) * t361 + t275;
t277 = mrSges(6,3) * t111;
t72 = -mrSges(6,2) * t361 + t277;
t282 = t71 + t72;
t274 = mrSges(7,3) * t112;
t73 = mrSges(7,1) * t361 - t274;
t276 = mrSges(6,3) * t112;
t74 = mrSges(6,1) * t361 - t276;
t281 = t73 + t74;
t89 = -t135 * t173 + t136 * t176;
t78 = t175 * t89;
t79 = pkin(4) * t139 - pkin(8) * t140 + t143;
t42 = t172 * t79 + t78;
t279 = mrSges(5,3) * t131;
t273 = Ifges(5,4) * t132;
t258 = t133 * t172;
t257 = t133 * t175;
t255 = t140 * t172;
t254 = t140 * t175;
t253 = t158 * t172;
t252 = t158 * t175;
t251 = t160 * t172;
t250 = t160 * t175;
t244 = qJDD(1) * t166;
t64 = -t139 * qJD(3) + qJD(4) * t336;
t95 = pkin(4) * t134 + pkin(8) * t133;
t240 = t172 * t95 + t175 * t64 + t79 * t246;
t68 = -mrSges(7,1) * t111 + mrSges(7,2) * t112;
t239 = -t68 + t343;
t237 = m(4) - t325;
t151 = -pkin(2) - t296;
t229 = t140 * t246;
t221 = -t97 * mrSges(5,1) + t96 * mrSges(5,2);
t17 = -t56 * mrSges(7,1) + t55 * mrSges(7,2);
t217 = -t247 / 0.2e1;
t215 = -t172 * t64 + t175 * t95;
t41 = -t172 * t89 + t175 * t79;
t211 = pkin(4) * t159 + pkin(8) * t157;
t208 = -mrSges(4,1) * t243 + mrSges(4,2) * t244;
t191 = -(-t144 * t166 + t156) * t166 + t122 * t168;
t190 = t152 * t159 - t157 * t170;
t189 = qJ(6) * t133 - qJD(6) * t140;
t24 = t104 * t176 - t173 * t105 - t109 * t249 - t110 * t248;
t117 = -t159 * t251 + t252;
t115 = t159 * t253 + t250;
t186 = t229 - t258;
t185 = t140 * t247 + t257;
t22 = -qJDD(4) * pkin(4) - t24;
t65 = qJD(3) * t140 + qJD(4) * t89;
t171 = -pkin(7) - qJ(3);
t146 = t170 * t175;
t145 = t170 * t172;
t142 = qJDD(1) * t151 + qJDD(3);
t119 = -qJD(4) * mrSges(5,2) - t279;
t118 = t159 * t250 + t253;
t116 = -t159 * t252 + t251;
t83 = t132 * Ifges(5,1) + t338 - t345;
t82 = t273 + t337 - t344;
t81 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t97;
t66 = pkin(5) * t255 - t336;
t36 = pkin(5) * t186 + t65;
t33 = -qJ(6) * t255 + t42;
t28 = pkin(5) * t139 - qJ(6) * t254 + t41;
t9 = -qJD(5) * t42 + t215;
t8 = -t247 * t89 + t240;
t7 = -pkin(5) * t56 + qJDD(6) + t22;
t6 = -qJ(6) * t229 + (-qJD(5) * t89 + t189) * t172 + t240;
t5 = pkin(5) * t134 + t189 * t175 + (-t78 + (qJ(6) * t140 - t79) * t172) * qJD(5) + t215;
t10 = [m(4) * (t191 * qJD(3) + t142 * t151 + t147 * t328) + (-t1 * t254 + t14 * t185 - t16 * t186 - t2 * t255) * mrSges(7,3) - (-m(5) * t24 + m(6) * t22 + t349) * t336 + (t185 * t26 - t186 * t27 - t254 * t4 - t255 * t3) * mrSges(6,3) + m(7) * (t1 * t28 + t14 * t5 + t16 * t6 + t2 * t33 + t36 * t40 + t66 * t7) + (-t185 * t372 - t186 * t354) * t308 + (-t185 * t357 - t186 * t372) * t306 + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t169 - 0.2e1 * mrSges(3,2) * t167 + m(3) * (t167 ^ 2 + t169 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + (-m(5) * t61 - t319) * t65 + t142 * t207 + t151 * t208 + t58 * (mrSges(6,1) * t186 - mrSges(6,2) * t185) + t40 * (mrSges(7,1) * t186 - mrSges(7,2) * t185) + (-t229 / 0.2e1 + t258 / 0.2e1) * t347 + (t125 * mrSges(5,2) - t24 * mrSges(5,3) + Ifges(5,1) * t96 + Ifges(5,4) * t97 + Ifges(5,5) * qJDD(4) + t201 * t7 + t203 * t22 + t217 * t346 + t310 * t334 + t311 * t333 + t312 * t332) * t140 + (-t62 * mrSges(5,3) - Ifges(5,4) * t301 - t82 / 0.2e1 + t344 / 0.2e1 + t353 * t308 + t355 * t306 + t352 * t304 - t314 + t348 / 0.2e1) * t134 + (t125 * mrSges(5,1) + t1 * mrSges(7,1) - t23 * mrSges(5,3) - Ifges(5,4) * t96 - Ifges(5,2) * t97 - Ifges(5,6) * qJDD(4) + t352 * t310 + t353 * t311 + t355 * t312 + t315 + t329 / 0.2e1) * t139 + (t367 * t137 + t328) * mrSges(4,3) + (t61 * mrSges(5,3) - Ifges(5,1) * t301 - t83 / 0.2e1 + t345 / 0.2e1 - t316) * t133 + (-t253 * t313 - mrSges(2,1) * t177 + mrSges(2,2) * t174 + t325 * (t160 * t150 - t158 * t171 + t161) + t360 * t161 + t359 * t118 - t358 * t117 + t317 * t158 + (-m(6) * t211 - m(7) * t190 - t318) * t160) * g(2) + (-t251 * t313 + mrSges(2,1) * t174 + mrSges(2,2) * t177 - t360 * t295 + t325 * (-t160 * t171 - t295) + t359 * t116 - t358 * t115 + t317 * t160 + (m(5) * t150 - m(7) * (-t150 - t190) - m(6) * (-t150 - t211) + t318) * t158) * g(1) + (-t185 * t355 - t186 * t353) * t304 - t346 * t257 / 0.2e1 - t351 * t255 / 0.2e1 + t254 * t369 + m(5) * (t125 * t143 + t23 * t89 + t62 * t64) + m(6) * (t26 * t9 + t27 * t8 + t3 * t42 + t4 * t41) + (Ifges(4,4) * t166 + Ifges(4,2) * t168) * t243 + (Ifges(4,1) * t166 + Ifges(4,4) * t168) * t244 + t28 * t29 + t33 * t31 + t41 * t30 + t42 * t32 + t66 * t17 + t36 * t68 + t6 * t71 + t8 * t72 + t5 * t73 + t9 * t74 + t89 * t81 + t143 * t221 + t64 * t119; m(3) * qJDD(2) + (t17 + t349) * t139 - t239 * t134 - (-t172 * t281 + t175 * t282 + t119) * t133 + (-m(3) - t237) * g(3) + m(5) * (-t133 * t62 - t134 * t61 - t139 * t24) + m(4) * (t113 * t168 + t114 * t166) + m(7) * (t134 * t40 + t139 * t7 + t14 * t258 - t16 * t257) + m(6) * (t134 * t58 + t139 * t22 - t257 * t27 + t258 * t26) + (t81 + t283 * t175 + t284 * t172 + (-t172 * t282 - t175 * t281) * qJD(5) + m(5) * t23 + m(7) * t320 + m(6) * t366) * t140; t131 * t119 - t367 * qJD(1) ^ 2 * mrSges(4,3) + t239 * t132 + (t282 * t361 - t284) * t175 + (-t281 * t361 + t283) * t172 + t208 + t221 + (-g(1) * t158 + g(2) * t160) * t237 + (t1 * t175 - t132 * t40 + t2 * t172 + t361 * (-t14 * t172 + t16 * t175)) * m(7) + (-t132 * t58 + t3 * t172 + t4 * t175 + t361 * (-t26 * t172 + t27 * t175)) * m(6) + (t131 * t62 + t132 * t61 + t125) * m(5) + (-qJD(1) * t191 + t142) * m(4); (-t345 + t83) * t303 + (t111 * t333 + t112 * t332 + t334 * t361) * qJD(5) / 0.2e1 + (-t74 * t246 - t72 * t247 + m(6) * ((-t27 * t172 - t26 * t175) * qJD(5) + t210) - t172 * t30 + t175 * t32) * pkin(8) + (-pkin(4) * t22 - t26 * t34 - t27 * t35) * m(6) + (t278 + t319) * t62 + (-t14 * t259 - t16 * t260 + t320) * mrSges(7,3) + (-t279 - t119) * t61 + t7 * t202 + t22 * t204 - g(3) * t206 + (t175 * t354 + t363) * t311 + (t172 * t357 + t364) * t312 + (-t259 * t26 - t260 * t27 + t366) * mrSges(6,3) + (t370 * g(3) + t362 * (mrSges(5,1) - t371)) * t157 + (t371 * g(3) + t362 * (mrSges(5,2) + t370)) * t159 + (t217 - t260 / 0.2e1) * t347 + (t246 / 0.2e1 + t259 / 0.2e1) * t346 + t340 * t71 + (t172 * t355 + t175 * t353) * t310 + (-Ifges(5,2) * t303 + t305 * t352 + t307 * t355 + t309 * t353 + t314) * t132 + (-t273 + t348) * t302 + t351 * t175 / 0.2e1 + t341 * t68 + t342 * t73 + (t1 * t145 + t14 * t342 - t146 * t2 - t152 * t7 + t16 * t340 + t341 * t40) * m(7) - (Ifges(5,1) * t302 + t305 * t334 + t307 * t332 + t309 * t333 - t316 - t335) * t131 + t335 * qJD(5) + t172 * t369 + t82 * t301 - pkin(4) * t18 - t23 * mrSges(5,2) + t24 * mrSges(5,1) - t35 * t72 - t34 * t74 + Ifges(5,5) * t96 + Ifges(5,6) * t97 + t145 * t29 - t146 * t31 - t152 * t17 + Ifges(5,3) * qJDD(4); (t111 * t355 - t112 * t353) * t305 + (t111 * t357 - t365) * t307 + t339 * t1 + t347 * t306 + (t276 + t74) * t27 + (t274 + t73 - m(7) * (-t14 + t15)) * t16 + t315 + (t277 - t72) * t26 + (t172 * t339 + t203 + t280) * g(3) * t157 + (-t112 * t354 + t346 + t368) * t309 + (-t117 * t322 + t118 * t358) * g(1) + (t115 * t322 - t116 * t358) * g(2) + t14 * t275 - t15 * t71 - t40 * (mrSges(7,1) * t112 + mrSges(7,2) * t111) - t58 * (mrSges(6,1) * t112 + mrSges(6,2) * t111) + t329 + ((-m(7) * t40 - t68) * t112 + t29) * pkin(5); -t111 * t71 + t112 * t73 + (g(3) * t159 - t16 * t111 + t14 * t112 - t157 * t362 + t7) * m(7) + t17;];
tau  = t10;
