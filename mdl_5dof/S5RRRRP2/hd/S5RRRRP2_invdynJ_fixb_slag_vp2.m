% Calculate vector of inverse dynamics joint torques for
% S5RRRRP2
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:48:45
% EndTime: 2022-01-20 11:48:55
% DurationCPUTime: 4.50s
% Computational Cost: add. (4858->401), mult. (7332->497), div. (0->0), fcn. (4451->12), ass. (0->197)
t357 = Ifges(5,4) + Ifges(6,4);
t358 = Ifges(5,1) + Ifges(6,1);
t356 = Ifges(6,5) + Ifges(5,5);
t355 = Ifges(5,2) + Ifges(6,2);
t354 = Ifges(6,6) + Ifges(5,6);
t241 = qJ(3) + qJ(4);
t228 = sin(t241);
t230 = cos(t241);
t359 = -mrSges(6,2) - mrSges(5,2);
t360 = -mrSges(6,1) - mrSges(5,1);
t263 = -t359 * t228 + t360 * t230;
t244 = sin(qJ(3));
t248 = cos(qJ(3));
t196 = -t248 * mrSges(4,1) + t244 * mrSges(4,2);
t350 = -mrSges(3,1) + t196;
t363 = t350 + t263;
t243 = sin(qJ(4));
t247 = cos(qJ(4));
t176 = -t243 * t244 + t247 * t248;
t238 = qJD(1) + qJD(2);
t150 = t176 * t238;
t362 = t357 * t150;
t177 = t243 * t248 + t244 * t247;
t151 = t177 * t238;
t361 = t357 * t151;
t325 = t244 / 0.2e1;
t237 = qJD(3) + qJD(4);
t353 = t355 * t150 + t354 * t237 + t361;
t352 = t358 * t151 + t356 * t237 + t362;
t251 = -pkin(8) - pkin(7);
t197 = t251 * t244;
t232 = t248 * pkin(8);
t198 = pkin(7) * t248 + t232;
t130 = t243 * t197 + t247 * t198;
t277 = qJD(3) * t251;
t182 = t244 * t277;
t183 = t248 * t277;
t249 = cos(qJ(2));
t306 = pkin(1) * qJD(1);
t280 = t249 * t306;
t345 = -qJD(4) * t130 + t177 * t280 - t182 * t243 + t247 * t183;
t282 = qJD(4) * t247;
t283 = qJD(4) * t243;
t343 = -t176 * t280 + t247 * t182 + t243 * t183 + t197 * t282 - t198 * t283;
t291 = t238 * t244;
t185 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t291;
t290 = t238 * t248;
t186 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t290;
t351 = (t238 * mrSges(3,2) + t244 * t185 - t248 * t186) * t249;
t298 = qJ(5) * t150;
t245 = sin(qJ(2));
t281 = t245 * t306;
t187 = pkin(7) * t238 + t281;
t273 = pkin(8) * t238 + t187;
t136 = t273 * t248;
t125 = t247 * t136;
t135 = t273 * t244;
t304 = qJD(3) * pkin(3);
t126 = -t135 + t304;
t76 = t126 * t243 + t125;
t42 = t76 + t298;
t269 = t76 * mrSges(5,3) + t42 * mrSges(6,3);
t242 = qJ(1) + qJ(2);
t229 = sin(t242);
t231 = cos(t242);
t342 = g(1) * t231 + g(2) * t229;
t348 = mrSges(3,2) - mrSges(6,3) - mrSges(5,3) - mrSges(4,3);
t188 = -pkin(2) * t238 - t280;
t268 = mrSges(4,1) * t244 + mrSges(4,2) * t248;
t347 = t188 * t268 + qJD(3) * (Ifges(4,5) * t248 - Ifges(4,6) * t244) / 0.2e1;
t259 = t176 * qJD(4);
t109 = qJD(3) * t176 + t259;
t264 = -qJ(5) * t109 - qJD(5) * t177;
t346 = t264 + t345;
t260 = t177 * qJD(4);
t110 = -qJD(3) * t177 - t260;
t287 = t110 * qJ(5) + t176 * qJD(5);
t344 = t287 + t343;
t142 = t151 * qJ(5);
t123 = t243 * t136;
t75 = t247 * t126 - t123;
t41 = -t142 + t75;
t218 = pkin(1) * t245 + pkin(7);
t314 = -pkin(8) - t218;
t170 = t314 * t244;
t294 = t218 * t248;
t171 = t232 + t294;
t104 = t243 * t170 + t247 * t171;
t341 = t350 * t238;
t339 = t188 * t245 + (t244 ^ 2 + t248 ^ 2) * t187 * t249;
t333 = t150 / 0.2e1;
t331 = t151 / 0.2e1;
t326 = t237 / 0.2e1;
t323 = pkin(1) * t249;
t322 = pkin(4) * t151;
t217 = pkin(4) * t230;
t319 = g(3) * t248;
t313 = mrSges(5,3) * t150;
t312 = mrSges(6,3) * t150;
t311 = Ifges(4,4) * t244;
t310 = Ifges(4,4) * t248;
t307 = Ifges(4,2) * t248;
t305 = pkin(1) * qJD(2);
t278 = qJD(1) * t305;
t299 = pkin(1) * qJDD(1);
t174 = t245 * t299 + t249 * t278;
t235 = qJDD(1) + qJDD(2);
t154 = pkin(7) * t235 + t174;
t284 = qJD(3) * t248;
t99 = -t154 * t244 - t187 * t284;
t302 = t244 * t99;
t285 = qJD(3) * t244;
t98 = t248 * t154 - t187 * t285;
t300 = t248 * t98;
t297 = qJ(5) * t177;
t78 = -t247 * t135 - t123;
t286 = t231 * pkin(2) + t229 * pkin(7);
t279 = t249 * t305;
t223 = pkin(3) * t285;
t220 = pkin(3) * t248 + pkin(2);
t161 = t235 * t248 - t238 * t285;
t162 = t235 * t244 + t238 * t284;
t68 = t161 * t243 + t162 * t247 + t238 * t259;
t69 = t161 * t247 - t162 * t243 - t238 * t260;
t16 = -t69 * mrSges(6,1) + t68 * mrSges(6,2);
t275 = t359 * t230;
t100 = -pkin(4) * t110 + t223;
t272 = qJD(3) * t314;
t77 = t135 * t243 - t125;
t103 = t247 * t170 - t171 * t243;
t181 = t217 + t220;
t236 = qJ(5) - t251;
t271 = t231 * t181 + t229 * t236;
t129 = t247 * t197 - t198 * t243;
t270 = t231 * t220 - t229 * t251;
t173 = -t245 * t278 + t249 * t299;
t267 = t307 + t311;
t265 = t185 * t248 + t186 * t244;
t147 = -pkin(4) * t176 - t220;
t261 = t244 * (Ifges(4,1) * t248 - t311);
t71 = qJDD(3) * pkin(3) - pkin(8) * t162 + t99;
t74 = pkin(8) * t161 + t98;
t7 = t126 * t282 - t136 * t283 + t243 * t71 + t247 * t74;
t131 = t244 * t272 + t248 * t279;
t132 = -t244 * t279 + t248 * t272;
t24 = t247 * t131 + t243 * t132 + t170 * t282 - t171 * t283;
t153 = -pkin(2) * t235 - t173;
t152 = -t220 * t238 - t280;
t101 = -pkin(3) * t161 + t153;
t8 = -qJD(4) * t76 - t243 * t74 + t247 * t71;
t25 = -qJD(4) * t104 - t131 * t243 + t247 * t132;
t257 = -qJD(3) * t265 + t248 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t161) - t244 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t162);
t255 = t348 * t229 + t363 * t231;
t234 = qJDD(3) + qJDD(4);
t3 = pkin(4) * t234 - qJ(5) * t68 - qJD(5) * t151 + t8;
t36 = pkin(4) * t237 + t41;
t4 = qJ(5) * t69 + qJD(5) * t150 + t7;
t95 = -pkin(4) * t150 + qJD(5) + t152;
t254 = t8 * mrSges(5,1) + t3 * mrSges(6,1) - t7 * mrSges(5,2) - t4 * mrSges(6,2) - t152 * (mrSges(5,1) * t151 + mrSges(5,2) * t150) + t36 * t312 - t95 * (mrSges(6,1) * t151 + mrSges(6,2) * t150) + t75 * t313 + t354 * t69 + t356 * t68 - (t358 * t150 - t361) * t151 / 0.2e1 + t353 * t331 - (t356 * t150 - t354 * t151) * t237 / 0.2e1 + (Ifges(6,3) + Ifges(5,3)) * t234 - (-t355 * t151 + t352 + t362) * t150 / 0.2e1;
t253 = (-m(4) * pkin(7) + m(5) * t251 - m(6) * t236 + t348) * t231 + (m(4) * pkin(2) + m(5) * t220 + m(6) * t181 - t363) * t229;
t148 = Ifges(4,6) * qJD(3) + t238 * t267;
t202 = Ifges(4,4) * t290;
t149 = Ifges(4,1) * t291 + Ifges(4,5) * qJD(3) + t202;
t23 = -pkin(4) * t69 + qJDD(5) + t101;
t252 = Ifges(3,3) * t235 + t153 * t196 + (0.2e1 * Ifges(4,5) * t325 + Ifges(4,6) * t248) * qJDD(3) + (t238 * (-Ifges(4,2) * t244 + t310) + t149) * t284 / 0.2e1 + t173 * mrSges(3,1) - t174 * mrSges(3,2) + t248 * (Ifges(4,4) * t162 + Ifges(4,2) * t161) / 0.2e1 + (t261 * t238 / 0.2e1 + t347) * qJD(3) + (Ifges(4,1) * t162 + Ifges(4,4) * t161) * t325 + t161 * t267 / 0.2e1 + t162 * (Ifges(4,1) * t244 + t310) / 0.2e1 - t148 * t285 / 0.2e1 + mrSges(4,3) * t300 + (t353 / 0.2e1 - mrSges(5,1) * t152 - mrSges(6,1) * t95 + t354 * t326 + t357 * t331 + t355 * t333 + t269) * t110 + (t352 / 0.2e1 + mrSges(5,2) * t152 + mrSges(6,2) * t95 + t356 * t326 + t358 * t331 + t357 * t333 - t36 * mrSges(6,3) - t75 * mrSges(5,3)) * t109 + (t101 * mrSges(5,2) + t23 * mrSges(6,2) - t8 * mrSges(5,3) - t3 * mrSges(6,3) + t356 * t234 + t357 * t69 + t358 * t68) * t177 + (-t101 * mrSges(5,1) - t23 * mrSges(6,1) + t7 * mrSges(5,3) + t4 * mrSges(6,3) + t354 * t234 + t355 * t69 + t357 * t68) * t176;
t250 = cos(qJ(1));
t246 = sin(qJ(1));
t233 = t250 * pkin(1);
t224 = t245 * t305;
t221 = -pkin(2) - t323;
t195 = -t220 - t323;
t184 = t224 + t223;
t169 = t176 * qJ(5);
t134 = t147 - t323;
t122 = mrSges(5,1) * t237 - mrSges(5,3) * t151;
t121 = mrSges(6,1) * t237 - mrSges(6,3) * t151;
t120 = -mrSges(5,2) * t237 + t313;
t119 = -mrSges(6,2) * t237 + t312;
t111 = pkin(3) * t291 + t322;
t102 = -mrSges(4,1) * t161 + mrSges(4,2) * t162;
t97 = t169 + t130;
t96 = t129 - t297;
t93 = -mrSges(5,1) * t150 + mrSges(5,2) * t151;
t92 = -mrSges(6,1) * t150 + mrSges(6,2) * t151;
t89 = t100 + t224;
t80 = t169 + t104;
t79 = t103 - t297;
t55 = -mrSges(5,2) * t234 + mrSges(5,3) * t69;
t54 = -mrSges(6,2) * t234 + mrSges(6,3) * t69;
t53 = mrSges(5,1) * t234 - mrSges(5,3) * t68;
t52 = mrSges(6,1) * t234 - mrSges(6,3) * t68;
t45 = -t142 + t78;
t44 = t77 - t298;
t17 = -mrSges(5,1) * t69 + mrSges(5,2) * t68;
t14 = t25 + t264;
t13 = t24 + t287;
t1 = [t221 * t102 + t195 * t17 + t184 * t93 + t257 * t218 + t134 * t16 + t13 * t119 + t24 * t120 + t14 * t121 + t25 * t122 + t103 * t53 + t104 * t55 + t79 * t52 + t80 * t54 + t89 * t92 + m(4) * (t153 * t221 - t218 * t302 + t98 * t294) + t252 + m(5) * (t101 * t195 + t103 * t8 + t104 * t7 + t152 * t184 + t24 * t76 + t25 * t75) + m(6) * (t13 * t42 + t134 * t23 + t14 * t36 + t3 * t79 + t4 * t80 + t89 * t95) + ((mrSges(3,1) * t249 - mrSges(3,2) * t245) * t235 + (-g(2) * t250 + t173 * t249 + t174 * t245) * m(3) + (m(3) + m(4) + m(5) + m(6)) * t246 * g(1) + (m(4) * t339 + t245 * t341 - t351) * qJD(2)) * pkin(1) + (-mrSges(2,1) * t250 + mrSges(2,2) * t246 - m(5) * (t233 + t270) - m(4) * (t233 + t286) - m(6) * (t233 + t271) + t255) * g(2) + (mrSges(2,1) * t246 + mrSges(2,2) * t250 + t253) * g(1) - mrSges(4,3) * t302 + Ifges(2,3) * qJDD(1); t255 * g(2) - t220 * t17 + t147 * t16 + t129 * t53 + t130 * t55 + t96 * t52 + t97 * t54 + t100 * t92 - pkin(2) * t102 + t346 * t121 + (t351 + (-t341 - t92 - t93) * t245) * t306 + t343 * t120 + t344 * t119 + t257 * pkin(7) + t345 * t122 + t253 * g(1) + t252 + (-t99 * mrSges(4,3) + t304 * t93) * t244 + (-t271 * g(2) + t147 * t23 + t3 * t96 + t4 * t97 + (t100 - t281) * t95 + t344 * t42 + t346 * t36) * m(6) + (-t270 * g(2) - t101 * t220 + t129 * t8 + t130 * t7 + t343 * t76 + t345 * t75 + (t223 - t281) * t152) * m(5) + (-t286 * g(2) - pkin(2) * t153 + (t300 - t302) * pkin(7) - t339 * t306) * m(4); Ifges(4,6) * t161 + Ifges(4,5) * t162 - t45 * t119 - t78 * t120 - t44 * t121 - t77 * t122 - t111 * t92 - t98 * mrSges(4,2) + t99 * mrSges(4,1) + (-m(6) * t217 + t196 + t263) * g(3) - m(6) * (t111 * t95 + t36 * t44 + t42 * t45) - m(5) * (t75 * t77 + t76 * t78) + t265 * t187 + (t148 * t325 + (-t261 / 0.2e1 + t307 * t325) * t238 - (t149 + t202) * t248 / 0.2e1 - t347) * t238 + t269 * t151 + t254 + (-m(6) * t319 - t93 * t291 + t247 * t53 + (m(6) * t4 + t54 + t55) * t243 + ((m(6) * t42 + t119 + t120) * t247 + (-m(6) * t36 - t121 - t122) * t243) * qJD(4) + (-t319 - t75 * t283 + t76 * t282 + t243 * t7 + t247 * t8 + (-t152 * t238 + t342) * t244) * m(5)) * pkin(3) + Ifges(4,3) * qJDD(3) + t342 * (-m(6) * (-pkin(3) * t244 - pkin(4) * t228) - t228 * t360 + t268 - t275) + (m(6) * t3 + t52) * (pkin(3) * t247 + pkin(4)); -t41 * t119 - t75 * t120 + t42 * t121 + t76 * t122 + pkin(4) * t52 + t263 * g(3) + t254 + (-pkin(4) * t92 + t269) * t151 + (-g(3) * t217 - t95 * t322 - (-t36 + t41) * t42 + t3 * pkin(4)) * m(6) + t342 * (-t275 + (m(6) * pkin(4) - t360) * t228); -t150 * t119 + t151 * t121 + (-g(1) * t229 + g(2) * t231 - t150 * t42 + t36 * t151 + t23) * m(6) + t16;];
tau = t1;
