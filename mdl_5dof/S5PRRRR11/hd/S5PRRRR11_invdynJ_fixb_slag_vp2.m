% Calculate vector of inverse dynamics joint torques for
% S5PRRRR11
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% m [6x1]
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
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRR11_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR11_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR11_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR11_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR11_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR11_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR11_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:14
% EndTime: 2024-09-27 21:45:20
% DurationCPUTime: 3.45s
% Computational Cost: add. (5757->472), mult. (16000->651), div. (0->0), fcn. (12287->22), ass. (0->242)
t330 = m(5) * pkin(3);
t350 = -mrSges(4,1) - t330;
t235 = cos(qJ(5));
t232 = sin(qJ(5));
t230 = sin(pkin(5));
t236 = cos(qJ(4));
t237 = cos(qJ(3));
t233 = sin(qJ(4));
t234 = sin(qJ(3));
t295 = t233 * t234;
t250 = t236 * t237 - t295;
t153 = t250 * t230;
t143 = qJD(2) * t153;
t313 = pkin(9) * t143;
t325 = pkin(7) + pkin(8);
t336 = t230 * t325;
t258 = t234 * t336;
t231 = cos(pkin(5));
t297 = t231 * t237;
t209 = pkin(2) * t297;
t284 = qJD(1) * t230;
t286 = qJD(2) * t209 + t237 * t284;
t108 = -qJD(2) * t258 + t286;
t206 = qJD(2) * t231 + qJD(3);
t93 = pkin(3) * t206 + t108;
t264 = t234 * t284;
t298 = t231 * t234;
t208 = pkin(2) * t298;
t347 = t237 * t336 + t208;
t109 = qJD(2) * t347 + t264;
t99 = t236 * t109;
t53 = t233 * t93 + t99;
t44 = t53 + t313;
t304 = t232 * t44;
t195 = qJD(4) + t206;
t293 = t233 * t237;
t251 = t234 * t236 + t293;
t154 = t251 * t230;
t144 = qJD(2) * t154;
t140 = t144 * pkin(9);
t97 = t233 * t109;
t52 = t236 * t93 - t97;
t43 = -t140 + t52;
t38 = pkin(4) * t195 + t43;
t13 = t235 * t38 - t304;
t303 = t235 * t44;
t14 = t232 * t38 + t303;
t183 = qJD(5) + t195;
t217 = t231 * qJDD(2);
t201 = t217 + qJDD(3);
t186 = qJDD(4) + t201;
t182 = qJDD(5) + t186;
t276 = qJD(2) * qJD(3);
t162 = (qJDD(2) * t234 + t237 * t276) * t230;
t299 = t230 * t237;
t285 = pkin(7) * t299 + t208;
t136 = qJD(2) * t285 + t264;
t274 = qJDD(2) * t230;
t268 = pkin(7) * t274;
t348 = pkin(2) * t217 + qJDD(1) * t230;
t79 = -qJD(3) * t136 - t234 * t268 + t237 * t348;
t54 = pkin(3) * t201 - pkin(8) * t162 + t79;
t161 = (qJDD(2) * t237 - t234 * t276) * t230;
t282 = qJD(3) * t234;
t265 = t230 * t282;
t259 = pkin(7) * t265;
t78 = -qJD(2) * t259 + qJD(3) * t286 + t234 * t348 + t237 * t268;
t64 = pkin(8) * t161 + t78;
t12 = -qJD(4) * t53 - t233 * t64 + t236 * t54;
t69 = qJD(4) * t143 + t161 * t233 + t162 * t236;
t7 = pkin(4) * t186 - pkin(9) * t69 + t12;
t279 = qJD(4) * t236;
t280 = qJD(4) * t233;
t11 = -t109 * t280 + t233 * t54 + t236 * t64 + t279 * t93;
t241 = t251 * qJD(4);
t283 = qJD(2) * t230;
t70 = t161 * t236 - t162 * t233 - t241 * t283;
t8 = pkin(9) * t70 + t11;
t2 = qJD(5) * t13 + t232 * t7 + t235 * t8;
t261 = t143 * t235 - t144 * t232;
t23 = qJD(5) * t261 + t232 * t70 + t235 * t69;
t85 = t143 * t232 + t144 * t235;
t24 = -qJD(5) * t85 - t232 * t69 + t235 * t70;
t3 = -qJD(5) * t14 - t232 * t8 + t235 * t7;
t248 = mrSges(6,1) * t3 - t2 * mrSges(6,2) + Ifges(6,5) * t23 + Ifges(6,6) * t24 + Ifges(6,3) * t182;
t318 = Ifges(6,4) * t85;
t81 = Ifges(6,4) * t261;
t41 = Ifges(6,1) * t85 + Ifges(6,5) * t183 + t81;
t216 = pkin(3) * t237 + pkin(2);
t175 = t216 * t230;
t220 = t231 * qJD(1);
t145 = -qJD(2) * t175 + t220;
t92 = -t143 * pkin(4) + t145;
t349 = t248 - (Ifges(6,5) * t261 - Ifges(6,6) * t85) * t183 / 0.2e1 + (t13 * t261 + t14 * t85) * mrSges(6,3) - (-Ifges(6,2) * t85 + t41 + t81) * t261 / 0.2e1 - t92 * (mrSges(6,1) * t85 + mrSges(6,2) * t261) - (Ifges(6,1) * t261 - t318) * t85 / 0.2e1;
t40 = Ifges(6,2) * t261 + Ifges(6,6) * t183 + t318;
t345 = t40 / 0.2e1;
t314 = pkin(4) * t231;
t344 = m(6) * t250 * t314;
t316 = pkin(3) * t236;
t215 = pkin(4) + t316;
t277 = qJD(5) * t235;
t278 = qJD(5) * t232;
t294 = t233 * t235;
t62 = -t108 * t233 - t99;
t48 = t62 - t313;
t63 = t108 * t236 - t97;
t49 = -t140 + t63;
t339 = t232 * t49 - t235 * t48 - t215 * t278 + (-t233 * t277 + (-t232 * t236 - t294) * qJD(4)) * pkin(3);
t296 = t232 * t233;
t338 = -t232 * t48 - t235 * t49 + t215 * t277 + (-t233 * t278 + (t235 * t236 - t296) * qJD(4)) * pkin(3);
t177 = -pkin(2) * t274 + qJDD(1) * t231;
t337 = -m(4) * t177 + mrSges(4,1) * t161 - mrSges(4,2) * t162;
t132 = pkin(3) * t231 + t209 - t258;
t142 = pkin(8) * t299 + t285;
t77 = t132 * t233 + t142 * t236;
t229 = qJ(3) + qJ(4);
t223 = pkin(5) + t229;
t204 = cos(t223) / 0.2e1;
t262 = pkin(5) - t229;
t210 = cos(t262);
t253 = sin(t262);
t272 = sin(t223) / 0.2e1;
t211 = qJ(5) + t223;
t190 = cos(t211) / 0.2e1;
t254 = -qJ(5) + t262;
t205 = cos(t254);
t249 = sin(t254);
t273 = sin(t211) / 0.2e1;
t288 = (t273 + t249 / 0.2e1) * mrSges(6,1) + (t190 - t205 / 0.2e1) * mrSges(6,2);
t335 = -(t272 + t253 / 0.2e1) * mrSges(5,1) - (t204 - t210 / 0.2e1) * mrSges(5,2) - t288;
t173 = t210 / 0.2e1 + t204;
t228 = pkin(10) + qJ(2);
t221 = sin(t228);
t222 = cos(t228);
t224 = sin(t229);
t300 = t222 * t224;
t121 = -t173 * t221 - t300;
t172 = t272 - t253 / 0.2e1;
t225 = cos(t229);
t164 = t205 / 0.2e1 + t190;
t226 = qJ(5) + t229;
t212 = sin(t226);
t107 = -t164 * t221 - t212 * t222;
t163 = t273 - t249 / 0.2e1;
t213 = cos(t226);
t291 = t107 * mrSges(6,1) + (t163 * t221 - t213 * t222) * mrSges(6,2);
t334 = -t121 * mrSges(5,1) - (t172 * t221 - t222 * t225) * mrSges(5,2) - t291;
t301 = t221 * t224;
t120 = -t173 * t222 + t301;
t106 = -t164 * t222 + t212 * t221;
t292 = -t106 * mrSges(6,1) + (-t163 * t222 - t213 * t221) * mrSges(6,2);
t333 = t120 * mrSges(5,1) - (-t172 * t222 - t221 * t225) * mrSges(5,2) - t292;
t332 = mrSges(5,1) * t12 - t11 * mrSges(5,2) + Ifges(5,5) * t69 + Ifges(5,6) * t70 + Ifges(5,3) * t186;
t331 = t79 * mrSges(4,1) - t78 * mrSges(4,2) + Ifges(4,5) * t162 + Ifges(4,6) * t161 + Ifges(4,3) * t201;
t327 = t85 / 0.2e1;
t326 = m(3) + m(2);
t322 = t144 / 0.2e1;
t320 = t237 / 0.2e1;
t317 = pkin(3) * t234;
t315 = pkin(4) * t144;
t311 = t153 * pkin(4);
t308 = mrSges(5,3) * t143;
t307 = Ifges(4,4) * t234;
t306 = t144 * mrSges(5,3);
t305 = t144 * Ifges(5,4);
t281 = qJD(3) * t237;
t267 = t234 * t283;
t266 = t237 * t283;
t76 = t132 * t236 - t142 * t233;
t260 = pkin(3) * t267;
t257 = mrSges(4,3) * t267;
t256 = mrSges(4,3) * t266;
t55 = -pkin(9) * t154 + t314 + t76;
t61 = pkin(9) * t153 + t77;
t27 = -t232 * t61 + t235 * t55;
t28 = t232 * t55 + t235 * t61;
t90 = t153 * t235 - t154 * t232;
t91 = t153 * t232 + t154 * t235;
t214 = pkin(4) * t236 + pkin(3);
t247 = -pkin(4) * t295 + t214 * t237;
t246 = -t221 * t234 + t222 * t297;
t149 = -t221 * t297 - t222 * t234;
t245 = (-pkin(2) * t283 + t220) * (mrSges(4,1) * t234 + mrSges(4,2) * t237);
t244 = t206 * (Ifges(4,5) * t237 - Ifges(4,6) * t234);
t243 = t234 * (Ifges(4,1) * t237 - t307);
t197 = qJD(3) * t209;
t133 = -qJD(3) * t258 + t197;
t134 = t347 * qJD(3);
t36 = t132 * t279 + t133 * t236 - t134 * t233 - t142 * t280;
t119 = -pkin(3) * t161 + t177;
t240 = m(5) * (pkin(3) * t298 - t336) + m(6) * (-t230 * (pkin(9) + t325) + (pkin(4) * t293 + t214 * t234) * t231) + t172 * mrSges(5,1) + t163 * mrSges(6,1) + mrSges(3,2);
t37 = -qJD(4) * t77 - t133 * t233 - t134 * t236;
t239 = m(4) * pkin(2) + m(5) * t216 + m(6) * (pkin(4) * t225 + t216) + t225 * mrSges(5,1) + t213 * mrSges(6,1) + mrSges(3,1);
t139 = Ifges(5,4) * t143;
t74 = t143 * Ifges(5,2) + t195 * Ifges(5,6) + t305;
t75 = t144 * Ifges(5,1) + t195 * Ifges(5,5) + t139;
t238 = t85 * t345 - t195 * (Ifges(5,5) * t143 - Ifges(5,6) * t144) / 0.2e1 + t306 * t53 - t145 * (mrSges(5,1) * t144 + mrSges(5,2) * t143) + t74 * t322 + t52 * t308 + t332 - t144 * (Ifges(5,1) * t143 - t305) / 0.2e1 - (-Ifges(5,2) * t144 + t139 + t75) * t143 / 0.2e1 + t349;
t194 = Ifges(4,4) * t266;
t180 = -pkin(4) * t224 - t317;
t170 = -pkin(7) * t230 * t234 + t209;
t169 = (-mrSges(4,1) * t237 + mrSges(4,2) * t234) * t230;
t166 = pkin(3) * t294 + t215 * t232;
t165 = -pkin(3) * t296 + t215 * t235;
t160 = t285 * qJD(3);
t159 = t197 - t259;
t158 = -mrSges(4,2) * t206 + t256;
t157 = mrSges(4,1) * t206 - t257;
t150 = -t221 * t298 + t222 * t237;
t148 = -t221 * t237 - t222 * t298;
t141 = t247 * t231;
t135 = -pkin(7) * t267 + t286;
t131 = Ifges(4,1) * t267 + Ifges(4,5) * t206 + t194;
t130 = Ifges(4,6) * t206 + (Ifges(4,2) * t237 + t307) * t283;
t129 = mrSges(4,1) * t201 - mrSges(4,3) * t162;
t128 = -mrSges(4,2) * t201 + mrSges(4,3) * t161;
t117 = t260 + t315;
t112 = -t175 - t311;
t111 = mrSges(5,1) * t195 - t306;
t110 = -mrSges(5,2) * t195 + t308;
t95 = (-qJD(3) * t251 - t241) * t230;
t94 = (qJD(3) + qJD(4)) * t153;
t88 = -mrSges(5,1) * t143 + mrSges(5,2) * t144;
t80 = pkin(3) * t265 - pkin(4) * t95;
t72 = mrSges(6,1) * t183 - mrSges(6,3) * t85;
t71 = -mrSges(6,2) * t183 + mrSges(6,3) * t261;
t59 = -mrSges(5,2) * t186 + mrSges(5,3) * t70;
t58 = mrSges(5,1) * t186 - mrSges(5,3) * t69;
t47 = -pkin(4) * t70 + t119;
t46 = -mrSges(6,1) * t261 + mrSges(6,2) * t85;
t35 = -mrSges(5,1) * t70 + mrSges(5,2) * t69;
t34 = -qJD(5) * t91 - t232 * t94 + t235 * t95;
t33 = qJD(5) * t90 + t232 * t95 + t235 * t94;
t30 = -pkin(9) * t94 + t37;
t29 = pkin(9) * t95 + t36;
t18 = -mrSges(6,2) * t182 + mrSges(6,3) * t24;
t17 = mrSges(6,1) * t182 - mrSges(6,3) * t23;
t16 = t235 * t43 - t304;
t15 = -t232 * t43 - t303;
t6 = -mrSges(6,1) * t24 + mrSges(6,2) * t23;
t5 = -qJD(5) * t28 - t232 * t29 + t235 * t30;
t4 = qJD(5) * t27 + t232 * t30 + t235 * t29;
t1 = [t94 * t110 + t95 * t111 + t153 * t58 + t154 * t59 + t90 * t17 + t91 * t18 + t33 * t71 + t34 * t72 + t326 * qJDD(1) + m(6) * (t13 * t34 + t14 * t33 + t2 * t91 + t3 * t90) + m(5) * (t11 * t154 + t12 * t153 + t52 * t95 + t53 * t94) + (-t157 * t282 + t158 * t281 + t234 * t128 + t237 * t129 + m(4) * (-t135 * t282 + t136 * t281 + t234 * t78 + t237 * t79)) * t230 + (-m(4) - m(5) - m(6) - t326) * g(3) + (m(5) * t119 + m(6) * t47 - t337 + t35 + t6) * t231; t27 * t17 + t28 * t18 + t33 * t41 / 0.2e1 + (-mrSges(6,1) * t47 + mrSges(6,3) * t2 + Ifges(6,4) * t23 + Ifges(6,2) * t24 + Ifges(6,6) * t182) * t90 + t4 * t71 + t5 * t72 + t76 * t58 + t77 * t59 + t80 * t46 + t92 * (-mrSges(6,1) * t34 + mrSges(6,2) * t33) + t94 * t75 / 0.2e1 + t95 * t74 / 0.2e1 + t36 * t110 + t37 * t111 + t112 * t6 + t143 * (Ifges(5,4) * t94 + Ifges(5,2) * t95) / 0.2e1 + t145 * (-mrSges(5,1) * t95 + mrSges(5,2) * t94) + t159 * t158 - t160 * t157 + t170 * t129 - t175 * t35 + t177 * t169 + (mrSges(6,2) * t47 - mrSges(6,3) * t3 + Ifges(6,1) * t23 + Ifges(6,4) * t24 + Ifges(6,5) * t182) * t91 + t183 * (Ifges(6,5) * t33 + Ifges(6,6) * t34) / 0.2e1 + t195 * (Ifges(5,5) * t94 + Ifges(5,6) * t95) / 0.2e1 + Ifges(3,3) * qJDD(2) + (-t13 * t33 + t14 * t34) * mrSges(6,3) + m(6) * (t112 * t47 + t13 * t5 + t14 * t4 + t2 * t28 + t27 * t3 + t80 * t92) + t261 * (Ifges(6,4) * t33 + Ifges(6,2) * t34) / 0.2e1 + t285 * t128 + m(4) * (-t135 * t160 + t136 * t159 + t170 * t79 + t285 * t78) + (Ifges(5,1) * t94 + Ifges(5,4) * t95) * t322 + (Ifges(6,1) * t33 + Ifges(6,4) * t34) * t327 + (-t150 * mrSges(4,1) - t149 * mrSges(4,2) - t121 * mrSges(5,2) - t107 * mrSges(6,2) + t221 * t240 - t222 * t239) * g(2) + (t337 * pkin(2) + (t78 * mrSges(4,3) + Ifges(4,4) * t162 + Ifges(4,2) * t161 + Ifges(4,6) * t201) * t237 + (-t79 * mrSges(4,3) + Ifges(4,1) * t162 + Ifges(4,4) * t161 + Ifges(4,5) * t201) * t234 + (-t234 * t130 / 0.2e1 + t131 * t320 + t245 + t244 / 0.2e1 + (t243 / 0.2e1 + (Ifges(4,4) * t237 - Ifges(4,2) * t234) * t320) * t283 + (m(5) * t145 + t88) * t317 + (-t135 * t237 - t136 * t234) * mrSges(4,3)) * qJD(3) + (g(1) * t222 + g(2) * t221) * (-m(4) * pkin(7) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3))) * t230 + (mrSges(5,2) * t119 - mrSges(5,3) * t12 + Ifges(5,1) * t69 + Ifges(5,4) * t70 + Ifges(5,5) * t186) * t154 + t34 * t345 + (-mrSges(5,1) * t119 + mrSges(5,3) * t11 + Ifges(5,4) * t69 + Ifges(5,2) * t70 + Ifges(5,6) * t186) * t153 + (t331 + t332 + t248) * t231 + (-t52 * t94 + t53 * t95) * mrSges(5,3) + (-t148 * mrSges(4,1) + mrSges(4,2) * t246 - t120 * mrSges(5,2) - t106 * mrSges(6,2) + t221 * t239 + t222 * t240) * g(1) + m(5) * (t11 * t77 - t119 * t175 + t12 * t76 + t36 * t53 + t37 * t52); (t157 + t257) * t136 - t63 * t110 - t62 * t111 - t117 * t46 + t165 * t17 + t166 * t18 + (t110 * t279 - t111 * t280 + t233 * t59) * pkin(3) + (t11 * t233 + t12 * t236 + (-t233 * t52 + t236 * t53) * qJD(4)) * t330 + t58 * t316 + (-t158 + t256) * t135 + t338 * t71 + t339 * t72 + (-t117 * t92 + t339 * t13 + t338 * t14 + t165 * t3 + t166 * t2) * m(6) + (-mrSges(4,2) * t148 - m(6) * (t141 * t222 + t180 * t221) + t333 + t350 * t246) * g(2) + (mrSges(4,2) * t150 - m(6) * (-t141 * t221 + t180 * t222) + t350 * t149 + t334) * g(1) + (-m(6) * t230 * t247 - t299 * t330 + t169 + t335) * g(3) + t331 - t245 * t283 + t130 * t267 / 0.2e1 - t230 ^ 2 * qJD(2) ^ 2 * t243 / 0.2e1 + t238 - t88 * t260 - m(5) * (t145 * t260 + t52 * t62 + t53 * t63) - (t244 + (-Ifges(4,2) * t267 + t131 + t194) * t237) * t283 / 0.2e1; -t46 * t315 - t16 * t71 - t15 * t72 - t52 * t110 + t53 * t111 - m(6) * (t13 * t15 + t14 * t16 + t315 * t92) + t238 + (-m(6) * t311 + t335) * g(3) + (-t222 * t344 + t333) * g(2) + (t221 * t344 + t334) * g(1) + (t235 * t17 + t232 * t18 + t71 * t277 - t72 * t278 + (g(1) * t300 + g(2) * t301 + t2 * t232 + t235 * t3 + (-t13 * t232 + t14 * t235) * qJD(5)) * m(6)) * pkin(4); -g(1) * t291 - g(2) * t292 - g(3) * t288 - t13 * t71 + t14 * t72 + t40 * t327 + t349;];
tau = t1;
