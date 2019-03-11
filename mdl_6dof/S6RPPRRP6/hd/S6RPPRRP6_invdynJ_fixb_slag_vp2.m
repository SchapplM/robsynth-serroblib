% Calculate vector of inverse dynamics joint torques for
% S6RPPRRP6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRP6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP6_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP6_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP6_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP6_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP6_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:10:12
% EndTime: 2019-03-09 02:10:31
% DurationCPUTime: 12.22s
% Computational Cost: add. (3332->518), mult. (6073->651), div. (0->0), fcn. (3127->6), ass. (0->233)
t362 = Ifges(6,1) + Ifges(7,1);
t343 = Ifges(6,5) + Ifges(7,4);
t361 = Ifges(6,6) - Ifges(7,6);
t137 = sin(qJ(5));
t140 = cos(qJ(5));
t234 = t140 * qJD(4);
t141 = cos(qJ(4));
t243 = qJD(1) * t141;
t91 = t137 * t243 - t234;
t263 = qJD(5) * t91;
t138 = sin(qJ(4));
t232 = qJD(1) * qJD(4);
t97 = qJDD(1) * t141 - t138 * t232;
t42 = qJDD(4) * t137 + t140 * t97 - t263;
t310 = t42 / 0.2e1;
t240 = qJD(4) * t137;
t92 = t140 * t243 + t240;
t43 = qJD(5) * t92 - t140 * qJDD(4) + t137 * t97;
t308 = t43 / 0.2e1;
t98 = -qJDD(1) * t138 - t141 * t232;
t87 = qJDD(5) - t98;
t307 = t87 / 0.2e1;
t305 = t91 / 0.2e1;
t366 = -t92 / 0.2e1;
t244 = qJD(1) * t138;
t112 = qJD(5) + t244;
t365 = -t112 / 0.2e1;
t364 = qJD(4) / 0.2e1;
t363 = mrSges(6,3) + mrSges(7,2);
t360 = Ifges(6,3) + Ifges(7,2);
t330 = -t137 * t361 + t140 * t343;
t272 = Ifges(7,5) * t137;
t274 = Ifges(6,4) * t137;
t328 = t140 * t362 + t272 - t274;
t185 = t140 * mrSges(7,1) + t137 * mrSges(7,3);
t187 = mrSges(6,1) * t140 - mrSges(6,2) * t137;
t359 = -t185 - t187;
t309 = -t43 / 0.2e1;
t358 = t343 * t307 + (-Ifges(6,4) + Ifges(7,5)) * t308 + t362 * t310;
t357 = -qJD(4) * mrSges(5,1) + mrSges(6,1) * t91 + mrSges(6,2) * t92 + mrSges(5,3) * t243;
t236 = qJD(5) * t140;
t237 = qJD(5) * t137;
t136 = pkin(1) + qJ(3);
t233 = qJD(1) * qJD(3);
t94 = qJDD(1) * t136 - qJDD(2) + t233;
t27 = -pkin(4) * t98 - pkin(8) * t97 + t94;
t121 = qJD(1) * qJ(2) + qJD(3);
t108 = -qJD(1) * pkin(7) + t121;
t238 = qJD(4) * t141;
t131 = qJD(1) * qJD(2);
t335 = qJDD(1) * qJ(2) + t131;
t106 = qJDD(3) + t335;
t93 = -qJDD(1) * pkin(7) + t106;
t52 = t108 * t238 + t138 * t93;
t45 = qJDD(4) * pkin(8) + t52;
t128 = t141 * pkin(8);
t199 = pkin(4) * t138 + t136;
t100 = -t128 + t199;
t64 = qJD(1) * t100 - qJD(2);
t99 = t138 * t108;
t73 = qJD(4) * pkin(8) + t99;
t3 = t137 * t27 + t140 * t45 + t64 * t236 - t237 * t73;
t25 = t137 * t64 + t140 * t73;
t4 = -qJD(5) * t25 - t137 * t45 + t140 * t27;
t190 = -t137 * t4 + t140 * t3;
t24 = -t137 * t73 + t140 * t64;
t356 = -t24 * t236 - t25 * t237 + t190;
t1 = qJ(6) * t87 + qJD(6) * t112 + t3;
t2 = -pkin(5) * t87 + qJDD(6) - t4;
t191 = t1 * t140 + t137 * t2;
t21 = -pkin(5) * t112 + qJD(6) - t24;
t22 = qJ(6) * t112 + t25;
t355 = t21 * t236 - t22 * t237 + t191;
t135 = qJ(2) - pkin(7);
t327 = qJD(2) * t138 + t135 * t238;
t239 = qJD(4) * t138;
t354 = -qJD(2) * t141 + t135 * t239;
t139 = sin(qJ(1));
t142 = cos(qJ(1));
t320 = g(1) * t142 + g(2) * t139;
t276 = Ifges(5,4) * t138;
t183 = t141 * Ifges(5,1) - t276;
t84 = Ifges(7,5) * t92;
t28 = Ifges(7,6) * t112 + Ifges(7,3) * t91 + t84;
t292 = Ifges(7,5) * t91;
t85 = Ifges(6,4) * t91;
t340 = t112 * t343 + t362 * t92 + t292 - t85;
t353 = Ifges(5,5) * qJD(4) + qJD(1) * t183 + t137 * t28 + t140 * t340;
t275 = Ifges(5,4) * t141;
t178 = -t138 * Ifges(5,2) + t275;
t352 = Ifges(5,6) * t364 + qJD(1) * t178 / 0.2e1 + t343 * t366 + t361 * t305 + t360 * t365;
t19 = -mrSges(6,2) * t87 - mrSges(6,3) * t43;
t20 = -mrSges(7,2) * t43 + mrSges(7,3) * t87;
t279 = t19 + t20;
t17 = mrSges(6,1) * t87 - mrSges(6,3) * t42;
t18 = -t87 * mrSges(7,1) + t42 * mrSges(7,2);
t280 = -t17 + t18;
t351 = t280 * t137 + t279 * t140;
t350 = t97 / 0.2e1;
t349 = t98 / 0.2e1;
t348 = m(5) + m(4);
t347 = m(6) + m(7);
t345 = -mrSges(4,2) - mrSges(3,3);
t344 = -mrSges(4,3) + mrSges(3,2);
t11 = mrSges(6,1) * t43 + mrSges(6,2) * t42;
t342 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t97 - t11;
t295 = mrSges(6,3) * t91;
t59 = -mrSges(6,2) * t112 - t295;
t297 = mrSges(7,2) * t91;
t62 = mrSges(7,3) * t112 - t297;
t278 = t59 + t62;
t294 = mrSges(6,3) * t92;
t60 = mrSges(6,1) * t112 - t294;
t296 = mrSges(7,2) * t92;
t61 = -mrSges(7,1) * t112 + t296;
t277 = t60 - t61;
t167 = pkin(5) * t137 - qJ(6) * t140;
t339 = -qJD(6) * t137 + t112 * t167 - t99;
t103 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t244;
t155 = t137 * t277 - t140 * t278;
t336 = -t103 + t155;
t333 = -t138 * t330 + t141 * t360;
t332 = -t138 * t328 + t141 * t343;
t331 = t137 * t343 + t140 * t361;
t271 = Ifges(7,5) * t140;
t273 = Ifges(6,4) * t140;
t329 = t137 * t362 - t271 + t273;
t324 = t343 * t42 + t360 * t87 - t361 * t43;
t322 = t363 * t141;
t51 = -t108 * t239 + t141 * t93;
t321 = -t138 * t52 - t141 * t51;
t319 = -t42 * Ifges(7,5) / 0.2e1 - t87 * Ifges(7,6) / 0.2e1 + Ifges(6,4) * t310 + Ifges(6,6) * t307 + (Ifges(7,3) + Ifges(6,2)) * t309;
t317 = mrSges(2,2) + mrSges(5,3) + t345;
t188 = t138 * mrSges(5,1) + t141 * mrSges(5,2);
t316 = -mrSges(2,1) - t188 + t344;
t315 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t314 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t249 = t141 * t142;
t251 = t139 * t141;
t254 = t138 * t142;
t290 = pkin(8) * t138;
t313 = -g(1) * (pkin(4) * t249 + pkin(8) * t254) - g(2) * (pkin(4) * t251 + t139 * t290);
t259 = t135 * t138;
t217 = qJD(5) * t259;
t193 = pkin(4) * t141 + t290;
t88 = qJD(4) * t193 + qJD(3);
t15 = -t137 * (qJD(5) * t100 + t327) - t140 * (-t88 + t217);
t311 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t143 = qJD(1) ^ 2;
t306 = -t91 / 0.2e1;
t303 = t92 / 0.2e1;
t300 = t112 / 0.2e1;
t293 = Ifges(6,4) * t92;
t282 = -qJD(1) / 0.2e1;
t281 = qJD(5) / 0.2e1;
t250 = t140 * t141;
t96 = t193 * qJD(1);
t47 = t108 * t250 + t137 * t96;
t268 = t140 * t96;
t255 = t138 * t140;
t58 = t137 * t100 + t135 * t255;
t262 = qJDD(1) * pkin(1);
t261 = t100 * t140;
t260 = t108 * t141;
t257 = t137 * t138;
t256 = t137 * t141;
t253 = t139 * t137;
t252 = t139 * t140;
t248 = t142 * t140;
t245 = t142 * pkin(1) + t139 * qJ(2);
t235 = qJD(5) * t141;
t49 = mrSges(7,1) * t91 - mrSges(7,3) * t92;
t229 = t49 + t357;
t226 = -t348 - t347;
t220 = t142 * qJ(3) + t245;
t216 = t137 * t239;
t206 = t236 / 0.2e1;
t205 = -t235 / 0.2e1;
t127 = t142 * qJ(2);
t204 = -pkin(7) * t142 + t127;
t203 = -t232 / 0.2e1;
t202 = t108 * (t138 ^ 2 + t141 ^ 2);
t200 = t100 * t236 + t137 * t88 + t140 * t327;
t194 = -t98 * mrSges(5,1) + t97 * mrSges(5,2);
t189 = mrSges(5,1) * t141 - mrSges(5,2) * t138;
t186 = mrSges(6,1) * t137 + mrSges(6,2) * t140;
t184 = t137 * mrSges(7,1) - t140 * mrSges(7,3);
t177 = -Ifges(6,2) * t137 + t273;
t176 = Ifges(6,2) * t140 + t274;
t173 = -Ifges(5,5) * t138 - Ifges(5,6) * t141;
t170 = Ifges(7,3) * t137 + t271;
t169 = -Ifges(7,3) * t140 + t272;
t168 = pkin(5) * t140 + qJ(6) * t137;
t166 = t137 * t22 - t140 * t21;
t165 = t137 * t25 + t140 * t24;
t74 = -qJD(4) * pkin(4) - t260;
t109 = qJD(1) * t136 - qJD(2);
t163 = qJD(3) * t109 + t136 * t94;
t162 = -pkin(7) * t139 + t220;
t101 = -pkin(4) - t168;
t44 = -qJDD(4) * pkin(4) - t51;
t161 = -t135 + t167;
t160 = t109 * t189;
t159 = t138 * (-Ifges(5,2) * t141 - t276);
t158 = t141 * (-Ifges(5,1) * t138 - t275);
t154 = -t137 * t278 - t140 * t277;
t149 = Ifges(6,6) * t141 - t138 * t177;
t148 = Ifges(7,6) * t141 - t138 * t170;
t122 = qJDD(2) - t262;
t95 = t188 * qJD(1);
t82 = t186 * t141;
t78 = t138 * t248 - t253;
t77 = t137 * t254 + t252;
t76 = t137 * t142 + t138 * t252;
t75 = t138 * t253 - t248;
t66 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t98;
t63 = t161 * t141;
t57 = -t135 * t257 + t261;
t55 = -t261 + (t135 * t137 - pkin(5)) * t138;
t54 = qJ(6) * t138 + t58;
t48 = pkin(5) * t92 + qJ(6) * t91;
t46 = -t108 * t256 + t268;
t41 = -t268 + (-pkin(5) * qJD(1) + t108 * t137) * t141;
t39 = qJ(6) * t243 + t47;
t31 = -Ifges(6,2) * t91 + Ifges(6,6) * t112 + t293;
t23 = pkin(5) * t91 - qJ(6) * t92 + t74;
t16 = -t161 * t239 + (qJD(5) * t168 - qJD(6) * t140 - qJD(2)) * t141;
t14 = -t137 * t217 + t200;
t13 = -pkin(5) * t238 - t15;
t12 = qJ(6) * t238 + (-t135 * t237 + qJD(6)) * t138 + t200;
t10 = mrSges(7,1) * t43 - mrSges(7,3) * t42;
t5 = pkin(5) * t43 - qJ(6) * t42 - qJD(6) * t92 + t44;
t6 = [-t353 * t239 / 0.2e1 + (t14 * t25 + t15 * t24 + t3 * t58 + t354 * t74 + t4 * t57) * m(6) + t159 * t203 + (mrSges(7,2) * t2 - mrSges(6,3) * t4 + t358) * t250 + (t24 * mrSges(6,1) - t21 * mrSges(7,1) - t25 * mrSges(6,2) + t22 * mrSges(7,3) - t352) * t238 + t94 * t188 + t340 * t137 * t205 + (t94 + t233) * mrSges(4,3) + qJD(3) * t95 + t44 * t82 + t57 * t17 + t58 * t19 + t14 * t59 + t15 * t60 + t13 * t61 + t12 * t62 + t63 * t10 + t16 * t49 + t136 * t194 + t54 * t20 + t55 * t18 + (t148 * t305 + t149 * t306 + t173 * t364 + t333 * t300 + t332 * t303 + t160) * qJD(4) + m(3) * (-pkin(1) * t122 + (t335 + t131) * qJ(2)) + (-m(5) * t204 - t347 * (pkin(8) * t251 + t204) + t315 * t76 + t314 * t75 + (-m(3) - m(4)) * t127 + t317 * t142 + (m(3) * pkin(1) + t136 * t348 + t199 * t347 - t316 - t322) * t139) * g(1) + m(7) * (t1 * t54 + t12 * t22 + t13 * t21 + t16 * t23 + t2 * t55 + t5 * t63) + (-Ifges(5,6) * qJDD(4) + t324 / 0.2e1 - Ifges(5,4) * t97 / 0.2e1 - Ifges(5,2) * t98 / 0.2e1 + Ifges(7,6) * t308 + Ifges(6,6) * t309 + t343 * t310 + t360 * t307 + t311) * t138 + (-m(3) * t245 - m(4) * t220 - m(5) * t162 - t347 * (pkin(4) * t254 - pkin(8) * t249 + t162) - t315 * t78 - t314 * t77 + t363 * t249 + t316 * t142 + t317 * t139) * g(2) + (mrSges(4,3) * t136 + Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1) + t158 * t232 / 0.2e1 + (Ifges(5,1) * t350 + Ifges(5,4) * t349 + Ifges(5,5) * qJDD(4) + t170 * t308 + t177 * t309 + t5 * t184 + t28 * t206 + t330 * t307 + t328 * t310 + (-m(6) * t44 + t342) * t135) * t141 + m(4) * (qJ(2) * t106 + qJD(2) * t121 + t163) + (-mrSges(6,1) * t74 - mrSges(7,1) * t23 + mrSges(7,2) * t22 + mrSges(6,3) * t25) * (-t140 * t235 + t216) + (-mrSges(6,2) * t74 - mrSges(7,2) * t21 + mrSges(6,3) * t24 + mrSges(7,3) * t23) * (t137 * t235 + t138 * t234) + t178 * t349 + t183 * t350 + t357 * t354 + (t122 - t262) * mrSges(3,2) + 0.2e1 * t335 * mrSges(3,3) + (t140 * t205 + t216 / 0.2e1) * t31 + (-t1 * mrSges(7,2) - t3 * mrSges(6,3) - t319) * t256 + t321 * mrSges(5,3) + m(5) * (qJD(2) * t202 - t135 * t321 + t163) + t327 * t103 + (-t169 * t305 - t176 * t306 - t300 * t331 - t303 * t329) * t235 + (t106 + t335) * mrSges(4,2) + t66 * t259; t280 * t140 - t279 * t137 + t344 * qJDD(1) + (-m(3) * qJ(2) + t345) * t143 + t155 * qJD(5) + m(6) * (-t3 * t137 - t4 * t140 + (t24 * t137 - t25 * t140) * qJD(5)) + m(7) * (-t1 * t137 + t2 * t140 + (-t21 * t137 - t22 * t140) * qJD(5)) + m(3) * t122 - m(5) * t94 + (t229 * t141 + t336 * t138 - m(6) * (-t141 * t74 - t24 * t257 + t25 * t255) - m(7) * (-t141 * t23 + t21 * t257 + t22 * t255) - m(5) * t202) * qJD(1) - t194 + (-g(1) * t139 + g(2) * t142) * (m(3) - t226) + (-qJD(1) * t121 - t94) * m(4); -t143 * mrSges(4,3) + m(4) * t106 + qJDD(1) * mrSges(4,2) + (-m(6) * t165 - m(7) * t166 - t109 * t348 + t154 - t95) * qJD(1) + (-t10 - t336 * qJD(4) + m(7) * (t21 * t240 + t22 * t234 - t5) + m(6) * (t234 * t25 - t24 * t240 - t44) + m(5) * t51 + t342) * t141 + (t66 + t229 * qJD(4) + t154 * qJD(5) + m(7) * (qJD(4) * t23 + t355) + m(6) * (qJD(4) * t74 + t356) + m(5) * t52 + t351) * t138 + t320 * t226; t353 * t244 / 0.2e1 + t355 * mrSges(7,2) + t356 * mrSges(6,3) + t173 * t203 + t352 * t243 + (-t160 - t21 * (-mrSges(7,1) * t141 - mrSges(7,2) * t255) - t24 * (mrSges(6,1) * t141 + mrSges(6,3) * t255) - t25 * (-mrSges(6,2) * t141 + mrSges(6,3) * t257) - t22 * (mrSges(7,2) * t257 + mrSges(7,3) * t141) + (t149 / 0.2e1 - t148 / 0.2e1) * t91) * qJD(1) + ((-m(7) * t168 + t359) * t141 - t363 * t138 - t189) * t320 + (-pkin(4) * t44 - t24 * t46 - t25 * t47 - t74 * t99 + t313) * m(6) + (t101 * t5 - t21 * t41 - t22 * t39 + t339 * t23 + t313) * m(7) + (-t277 * t236 - t278 * t237 + (-qJD(5) * t165 + t190) * m(6) + (-qJD(5) * t166 + t191) * m(7) + t351) * pkin(8) + t101 * t10 + Ifges(5,5) * t97 + Ifges(5,6) * t98 - t47 * t59 - t46 * t60 - t41 * t61 - t39 * t62 + t51 * mrSges(5,1) - t52 * mrSges(5,2) - t5 * t185 - t44 * t187 - pkin(4) * t11 - t103 * t260 - t357 * t99 + (t188 + (m(6) * pkin(4) - m(7) * t101 - t359) * t138 - t347 * t128 - t322) * g(3) + t28 * t237 / 0.2e1 + (t281 * t328 + t282 * t332) * t92 + t137 * t358 + (t184 * t23 + t186 * t74 - t137 * t31 / 0.2e1 + t281 * t330 + t282 * t333) * t112 + (t159 / 0.2e1 - t158 / 0.2e1) * t143 + Ifges(5,3) * qJDD(4) + t319 * t140 + t329 * t310 + t331 * t307 + (-t177 / 0.2e1 + t170 / 0.2e1) * t263 + t339 * t49 + t340 * t206 + t169 * t308 + t176 * t309; (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t22 - t23 * t48) * m(7) + (-m(7) * t21 + t277 + t294) * t25 + (-m(7) * t22 - t278 - t295) * t24 + t311 + t324 - t74 * (mrSges(6,1) * t92 - mrSges(6,2) * t91) - t23 * (mrSges(7,1) * t92 + mrSges(7,3) * t91) + qJD(6) * t62 - t48 * t49 + qJ(6) * t20 - pkin(5) * t18 + (-t343 * t91 - t361 * t92) * t365 + (-t362 * t91 + t28 - t293 + t84) * t366 + (t82 - (-m(7) * t167 - t184) * t141) * g(3) + (-Ifges(6,2) * t92 + t340 - t85) * t305 + (-t314 * t78 + t315 * t77) * g(1) + (-t314 * t76 + t315 * t75) * g(2) + t22 * t296 + t21 * t297 + t31 * t303 + (Ifges(7,3) * t92 - t292) * t306; -t112 * t62 + t92 * t49 + (-g(1) * t77 - g(2) * t75 - g(3) * t256 - t22 * t112 + t23 * t92 + t2) * m(7) + t18;];
tau  = t6;
