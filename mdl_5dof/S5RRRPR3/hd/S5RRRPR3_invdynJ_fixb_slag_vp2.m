% Calculate vector of inverse dynamics joint torques for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR3_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR3_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:42:18
% EndTime: 2019-12-05 18:42:29
% DurationCPUTime: 4.83s
% Computational Cost: add. (6393->437), mult. (9682->580), div. (0->0), fcn. (6446->16), ass. (0->225)
t257 = sin(qJ(3));
t337 = t257 / 0.2e1;
t261 = cos(qJ(3));
t240 = t261 * qJD(4);
t255 = -qJ(4) - pkin(7);
t286 = qJD(3) * t255;
t173 = t257 * t286 + t240;
t174 = -qJD(4) * t257 + t261 * t286;
t253 = sin(pkin(9));
t254 = cos(pkin(9));
t188 = t253 * t261 + t254 * t257;
t262 = cos(qJ(2));
t316 = qJD(1) * pkin(1);
t293 = t262 * t316;
t348 = -t173 * t253 + t254 * t174 + t188 * t293;
t187 = -t253 * t257 + t254 * t261;
t347 = t254 * t173 + t253 * t174 - t187 * t293;
t210 = -t261 * mrSges(4,1) + mrSges(4,2) * t257;
t371 = -mrSges(3,1) + t210;
t177 = t187 * qJD(3);
t328 = pkin(8) * t177;
t370 = -t328 + t348;
t176 = t188 * qJD(3);
t170 = t176 * pkin(8);
t369 = t170 - t347;
t249 = qJ(3) + pkin(9);
t237 = sin(t249);
t238 = cos(t249);
t368 = t238 * mrSges(5,1) - t237 * mrSges(5,2);
t256 = sin(qJ(5));
t260 = cos(qJ(5));
t248 = qJD(1) + qJD(2);
t160 = t188 * t248;
t329 = pkin(8) * t160;
t258 = sin(qJ(2));
t294 = t258 * t316;
t201 = pkin(7) * t248 + t294;
t283 = qJ(4) * t248 + t201;
t145 = t283 * t261;
t129 = t253 * t145;
t144 = t283 * t257;
t315 = qJD(3) * pkin(3);
t135 = -t144 + t315;
t75 = t254 * t135 - t129;
t50 = qJD(3) * pkin(4) - t329 + t75;
t159 = t187 * t248;
t330 = pkin(8) * t159;
t299 = t254 * t145;
t76 = t253 * t135 + t299;
t52 = t76 + t330;
t22 = -t256 * t52 + t260 * t50;
t23 = t256 * t50 + t260 * t52;
t252 = qJ(1) + qJ(2);
t241 = sin(t252);
t244 = qJDD(3) + qJDD(5);
t247 = qJD(3) + qJD(5);
t282 = t260 * t159 - t160 * t256;
t245 = qJDD(1) + qJDD(2);
t296 = qJD(3) * t257;
t178 = t245 * t261 - t248 * t296;
t295 = qJD(3) * t261;
t179 = t245 * t257 + t248 * t295;
t106 = t178 * t253 + t179 * t254;
t290 = qJD(2) * t316;
t308 = qJDD(1) * pkin(1);
t186 = t258 * t308 + t262 * t290;
t166 = pkin(7) * t245 + t186;
t289 = t201 * t295;
t65 = -t289 + qJDD(3) * pkin(3) - qJ(4) * t179 + (-qJD(4) * t248 - t166) * t257;
t108 = t261 * t166 - t201 * t296;
t70 = qJ(4) * t178 + t240 * t248 + t108;
t31 = -t253 * t70 + t254 * t65;
t18 = qJDD(3) * pkin(4) - pkin(8) * t106 + t31;
t105 = t178 * t254 - t179 * t253;
t32 = t253 * t65 + t254 * t70;
t19 = pkin(8) * t105 + t32;
t3 = qJD(5) * t22 + t18 * t256 + t19 * t260;
t239 = qJ(5) + t249;
t225 = sin(t239);
t304 = t225 * t241;
t226 = cos(t239);
t321 = mrSges(6,2) * t226;
t92 = t159 * t256 + t160 * t260;
t335 = Ifges(6,4) * t92;
t36 = qJD(5) * t282 + t105 * t256 + t106 * t260;
t37 = -qJD(5) * t92 + t105 * t260 - t106 * t256;
t4 = -qJD(5) * t23 + t18 * t260 - t19 * t256;
t86 = Ifges(6,4) * t282;
t41 = Ifges(6,1) * t92 + Ifges(6,5) * t247 + t86;
t331 = pkin(3) * t261;
t230 = pkin(2) + t331;
t156 = -t230 * t248 + qJD(4) - t293;
t99 = -pkin(4) * t159 + t156;
t367 = t4 * mrSges(6,1) + Ifges(6,3) * t244 + Ifges(6,6) * t37 + Ifges(6,5) * t36 - g(2) * (mrSges(6,1) * t304 + t241 * t321) - t3 * mrSges(6,2) - (Ifges(6,5) * t282 - Ifges(6,6) * t92) * t247 / 0.2e1 + (t22 * t282 + t23 * t92) * mrSges(6,3) - (-Ifges(6,2) * t92 + t41 + t86) * t282 / 0.2e1 - t99 * (mrSges(6,1) * t92 + mrSges(6,2) * t282) - (Ifges(6,1) * t282 - t335) * t92 / 0.2e1;
t324 = qJD(3) / 0.2e1;
t227 = pkin(3) * t254 + pkin(4);
t333 = pkin(3) * t253;
t172 = t227 * t256 + t260 * t333;
t80 = t144 * t253 - t299;
t54 = t80 - t330;
t81 = -t254 * t144 - t129;
t55 = t81 - t329;
t363 = -t172 * qJD(5) + t256 * t55 - t260 * t54;
t171 = t227 * t260 - t256 * t333;
t362 = t171 * qJD(5) - t256 * t54 - t260 * t55;
t217 = t226 * mrSges(6,1);
t287 = pkin(4) * t238 + t331;
t361 = m(4) * pkin(2) + m(5) * t230 + m(6) * (pkin(2) + t287) + t217 + t368 - t371;
t360 = -m(4) * pkin(7) - mrSges(6,3) - mrSges(5,3) - mrSges(4,3) + mrSges(3,2) + m(6) * (-pkin(8) + t255) + m(5) * t255;
t40 = Ifges(6,2) * t282 + Ifges(6,6) * t247 + t335;
t357 = t40 / 0.2e1;
t209 = t255 * t257;
t243 = t261 * qJ(4);
t211 = pkin(7) * t261 + t243;
t133 = t254 * t209 - t211 * t253;
t327 = pkin(8) * t188;
t103 = t133 - t327;
t134 = t253 * t209 + t254 * t211;
t182 = t187 * pkin(8);
t104 = t182 + t134;
t45 = t103 * t260 - t104 * t256;
t354 = qJD(5) * t45 + t370 * t256 - t369 * t260;
t46 = t103 * t256 + t104 * t260;
t353 = -qJD(5) * t46 + t369 * t256 + t370 * t260;
t322 = mrSges(4,2) * t261;
t332 = pkin(3) * t257;
t344 = m(5) * pkin(3);
t346 = -t257 * (mrSges(4,1) + t344) + m(6) * (-pkin(4) * t237 - t332) - mrSges(5,1) * t237 - mrSges(5,2) * t238 - t322;
t302 = t248 * t257;
t199 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t302;
t301 = t248 * t261;
t200 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t301;
t202 = -pkin(2) * t248 - t293;
t345 = (mrSges(3,2) * t248 + t199 * t257 - t200 * t261) * t262 - m(4) * (t202 * t258 + (t257 ^ 2 + t261 ^ 2) * t262 * t201);
t341 = t92 / 0.2e1;
t339 = t160 / 0.2e1;
t334 = pkin(1) * t262;
t242 = cos(t252);
t326 = g(3) * t242;
t325 = -qJD(3) / 0.2e1;
t320 = Ifges(4,4) * t257;
t319 = Ifges(4,4) * t261;
t318 = Ifges(5,4) * t160;
t317 = pkin(1) * qJD(2);
t314 = t225 * mrSges(6,2);
t310 = t261 * Ifges(4,2);
t307 = t108 * t261;
t109 = -t166 * t257 - t289;
t306 = t109 * t257;
t229 = pkin(1) * t258 + pkin(7);
t303 = t229 * t261;
t297 = -qJ(4) - t229;
t281 = qJD(3) * t297;
t292 = t262 * t317;
t121 = t257 * t281 + t261 * t292 + t240;
t122 = (-qJD(4) - t292) * t257 + t261 * t281;
t67 = t254 * t121 + t253 * t122;
t183 = t297 * t257;
t184 = t243 + t303;
t112 = t253 * t183 + t254 * t184;
t233 = pkin(3) * t296;
t9 = -t37 * mrSges(6,1) + t36 * mrSges(6,2);
t288 = t295 / 0.2e1;
t142 = pkin(4) * t176 + t233;
t285 = t371 * t248;
t49 = -t105 * mrSges(5,1) + t106 * mrSges(5,2);
t284 = t217 - t314;
t66 = -t121 * t253 + t254 * t122;
t111 = t254 * t183 - t184 * t253;
t185 = -t258 * t290 + t262 * t308;
t280 = -g(2) * t242 - g(3) * t241;
t279 = -mrSges(6,1) * t225 - t321;
t278 = t310 + t320;
t277 = Ifges(4,5) * t261 - Ifges(4,6) * t257;
t84 = t111 - t327;
t85 = t182 + t112;
t38 = -t256 * t85 + t260 * t84;
t39 = t256 * t84 + t260 * t85;
t115 = t187 * t260 - t188 * t256;
t116 = t187 * t256 + t188 * t260;
t276 = t199 * t261 + t200 * t257;
t152 = -pkin(4) * t187 - t230;
t275 = mrSges(2,1) + (m(3) + m(4) + m(5) + m(6)) * pkin(1);
t274 = t202 * (mrSges(4,1) * t257 + t322);
t273 = t257 * (Ifges(4,1) * t261 - t320);
t165 = -pkin(2) * t245 - t185;
t110 = -pkin(3) * t178 + qJDD(4) + t165;
t268 = -qJD(3) * t276 + t261 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t178) - t257 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t179);
t266 = (-t314 + t361) * t242 - t360 * t241;
t265 = -mrSges(6,2) * t304 + t241 * t361 + t242 * t360;
t163 = Ifges(4,6) * qJD(3) + t248 * t278;
t216 = Ifges(4,4) * t301;
t164 = Ifges(4,1) * t302 + Ifges(4,5) * qJD(3) + t216;
t51 = -pkin(4) * t105 + t110;
t57 = qJD(5) * t115 - t176 * t256 + t177 * t260;
t58 = -qJD(5) * t116 - t176 * t260 - t177 * t256;
t87 = Ifges(5,2) * t159 + Ifges(5,6) * qJD(3) + t318;
t149 = Ifges(5,4) * t159;
t88 = Ifges(5,1) * t160 + Ifges(5,5) * qJD(3) + t149;
t264 = (mrSges(6,2) * t51 - mrSges(6,3) * t4 + Ifges(6,1) * t36 + Ifges(6,4) * t37 + Ifges(6,5) * t244) * t116 + (t273 * t324 + (-Ifges(4,2) * t257 + t319) * t288) * t248 + (-t22 * t57 + t23 * t58) * mrSges(6,3) + t178 * t278 / 0.2e1 + (-t176 * t76 - t177 * t75 + t187 * t32 - t188 * t31) * mrSges(5,3) + (-mrSges(6,1) * t51 + mrSges(6,3) * t3 + Ifges(6,4) * t36 + Ifges(6,2) * t37 + Ifges(6,6) * t244) * t115 + (Ifges(4,1) * t179 + Ifges(4,4) * t178) * t337 + t58 * t357 + (Ifges(6,1) * t57 + Ifges(6,4) * t58) * t341 + mrSges(4,3) * t307 + (Ifges(5,1) * t177 - Ifges(5,4) * t176) * t339 + (Ifges(5,5) * t177 - Ifges(5,6) * t176) * t324 + t156 * (mrSges(5,1) * t176 + mrSges(5,2) * t177) + t159 * (Ifges(5,4) * t177 - Ifges(5,2) * t176) / 0.2e1 + t282 * (Ifges(6,4) * t57 + Ifges(6,2) * t58) / 0.2e1 + t261 * (Ifges(4,4) * t179 + Ifges(4,2) * t178) / 0.2e1 + t164 * t288 + (t277 * t324 + t274) * qJD(3) + t57 * t41 / 0.2e1 - t163 * t296 / 0.2e1 + t99 * (-mrSges(6,1) * t58 + mrSges(6,2) * t57) + t179 * (t257 * Ifges(4,1) + t319) / 0.2e1 - t176 * t87 / 0.2e1 + t177 * t88 / 0.2e1 + t185 * mrSges(3,1) - t186 * mrSges(3,2) + t110 * (-mrSges(5,1) * t187 + mrSges(5,2) * t188) + t106 * (Ifges(5,1) * t188 + Ifges(5,4) * t187) + t105 * (Ifges(5,4) * t188 + Ifges(5,2) * t187) + t165 * t210 + Ifges(3,3) * t245 + t247 * (Ifges(6,5) * t57 + Ifges(6,6) * t58) / 0.2e1 + (0.2e1 * Ifges(4,5) * t337 + Ifges(5,5) * t188 + Ifges(4,6) * t261 + Ifges(5,6) * t187) * qJDD(3);
t263 = cos(qJ(1));
t259 = sin(qJ(1));
t234 = t258 * t317;
t231 = -pkin(2) - t334;
t206 = -t230 - t334;
t198 = t234 + t233;
t140 = t152 - t334;
t137 = qJD(3) * mrSges(5,1) - t160 * mrSges(5,3);
t136 = -qJD(3) * mrSges(5,2) + t159 * mrSges(5,3);
t125 = t142 + t234;
t120 = pkin(3) * t302 + pkin(4) * t160;
t113 = -mrSges(4,1) * t178 + mrSges(4,2) * t179;
t95 = -mrSges(5,1) * t159 + mrSges(5,2) * t160;
t94 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t106;
t93 = -qJDD(3) * mrSges(5,2) + mrSges(5,3) * t105;
t79 = mrSges(6,1) * t247 - mrSges(6,3) * t92;
t78 = -mrSges(6,2) * t247 + mrSges(6,3) * t282;
t48 = -t170 + t67;
t47 = t66 - t328;
t44 = -mrSges(6,1) * t282 + mrSges(6,2) * t92;
t30 = -mrSges(6,2) * t244 + mrSges(6,3) * t37;
t29 = mrSges(6,1) * t244 - mrSges(6,3) * t36;
t8 = -qJD(5) * t39 - t256 * t48 + t260 * t47;
t7 = qJD(5) * t38 + t256 * t47 + t260 * t48;
t1 = [m(4) * (t108 * t303 + t165 * t231 - t229 * t306) + (mrSges(2,2) * t263 + t259 * t275 + t265) * g(3) + t268 * t229 + m(6) * (t125 * t99 + t140 * t51 + t22 * t8 + t23 * t7 + t3 * t39 + t38 * t4) + m(5) * (t110 * t206 + t111 * t31 + t112 * t32 + t156 * t198 + t66 * t75 + t67 * t76) + t264 + (-t259 * mrSges(2,2) + t263 * t275 + t266) * g(2) + t38 * t29 + t39 * t30 - mrSges(4,3) * t306 + t7 * t78 + t8 * t79 + t111 * t94 + t112 * t93 + t125 * t44 + t67 * t136 + t66 * t137 + t140 * t9 + t198 * t95 + t206 * t49 + t231 * t113 + Ifges(2,3) * qJDD(1) + ((mrSges(3,1) * t262 - mrSges(3,2) * t258) * t245 + m(3) * (t185 * t262 + t186 * t258) + (t258 * t285 - t345) * qJD(2)) * pkin(1); t354 * t78 + m(4) * (-pkin(2) * t165 + (-t306 + t307) * pkin(7)) + t353 * t79 + (-t109 * mrSges(4,3) + t315 * t95) * t257 + t348 * t137 + t347 * t136 + t264 + t45 * t29 + t46 * t30 + t265 * g(3) + t268 * pkin(7) + t266 * g(2) - pkin(2) * t113 + t133 * t94 + t134 * t93 + t142 * t44 + t152 * t9 - t230 * t49 + ((-t285 - t44 - t95) * t258 + t345) * t316 + (t152 * t51 + t3 * t46 + t4 * t45 + (t142 - t294) * t99 + t354 * t23 + t353 * t22) * m(6) + (-t110 * t230 + t133 * t31 + t134 * t32 + t347 * t76 + t348 * t75 + (t233 - t294) * t156) * m(5); t362 * t78 + (-g(1) * t287 - t120 * t99 + t171 * t4 + t172 * t3 + t22 * t363 + t23 * t362) * m(6) + t363 * t79 + (-m(5) * t331 + t210 - t284 - t368) * g(1) + t367 + t276 * t201 - m(5) * (t75 * t80 + t76 * t81) + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (t253 * t93 + t254 * t94) * pkin(3) + (Ifges(5,5) * t159 - Ifges(5,6) * t160) * t325 - (-Ifges(5,2) * t160 + t149 + t88) * t159 / 0.2e1 + (t159 * t75 + t160 * t76) * mrSges(5,3) - t156 * (mrSges(5,1) * t160 + mrSges(5,2) * t159) - t160 * (Ifges(5,1) * t159 - t318) / 0.2e1 + (t253 * t32 + t254 * t31) * t344 + t87 * t339 + t92 * t357 + (t163 * t337 - t274 + t277 * t325 + (t310 * t337 - t273 / 0.2e1) * t248 + (-m(5) * t156 - t95) * t332 - (t164 + t216) * t261 / 0.2e1) * t248 + (-t279 - t346) * t326 + t346 * g(2) * t241 + t31 * mrSges(5,1) - t32 * mrSges(5,2) + Ifges(5,6) * t105 + Ifges(5,5) * t106 - t108 * mrSges(4,2) + t109 * mrSges(4,1) - t120 * t44 - t81 * t136 - t80 * t137 + t171 * t29 + t172 * t30 + Ifges(4,6) * t178 + Ifges(4,5) * t179; -t159 * t136 + t160 * t137 - t282 * t78 + t92 * t79 + t49 + t9 + (t22 * t92 - t23 * t282 + t280 + t51) * m(6) + (-t159 * t76 + t160 * t75 + t110 + t280) * m(5); -g(1) * t284 - t22 * t78 + t23 * t79 - t279 * t326 + t40 * t341 + t367;];
tau = t1;
