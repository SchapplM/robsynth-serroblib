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
% Datum: 2019-12-05 18:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:47:33
% EndTime: 2019-12-05 18:47:43
% DurationCPUTime: 4.21s
% Computational Cost: add. (4858->396), mult. (7332->489), div. (0->0), fcn. (4451->12), ass. (0->195)
t348 = Ifges(5,4) + Ifges(6,4);
t349 = Ifges(5,1) + Ifges(6,1);
t347 = Ifges(6,5) + Ifges(5,5);
t346 = Ifges(5,2) + Ifges(6,2);
t345 = Ifges(6,6) + Ifges(5,6);
t241 = sin(qJ(3));
t245 = cos(qJ(3));
t198 = -t245 * mrSges(4,1) + mrSges(4,2) * t241;
t354 = -mrSges(3,1) + t198;
t238 = qJ(3) + qJ(4);
t226 = sin(t238);
t310 = mrSges(5,2) + mrSges(6,2);
t353 = t310 * t226;
t228 = cos(t238);
t275 = t310 * t228;
t311 = mrSges(5,1) + mrSges(6,1);
t352 = t311 * t228;
t240 = sin(qJ(4));
t244 = cos(qJ(4));
t175 = -t240 * t241 + t244 * t245;
t235 = qJD(1) + qJD(2);
t150 = t175 * t235;
t351 = t348 * t150;
t176 = t240 * t245 + t241 * t244;
t151 = t176 * t235;
t350 = t348 * t151;
t318 = t241 / 0.2e1;
t234 = qJD(3) + qJD(4);
t344 = t346 * t150 + t345 * t234 + t350;
t343 = t349 * t151 + t347 * t234 + t351;
t248 = -pkin(8) - pkin(7);
t199 = t248 * t241;
t230 = t245 * pkin(8);
t200 = pkin(7) * t245 + t230;
t130 = t240 * t199 + t244 * t200;
t276 = qJD(3) * t248;
t181 = t241 * t276;
t182 = t245 * t276;
t246 = cos(qJ(2));
t299 = qJD(1) * pkin(1);
t279 = t246 * t299;
t335 = -qJD(4) * t130 + t176 * t279 - t181 * t240 + t244 * t182;
t281 = qJD(4) * t244;
t282 = qJD(4) * t240;
t333 = -t175 * t279 + t244 * t181 + t240 * t182 + t199 * t281 - t200 * t282;
t215 = pkin(4) * t228;
t218 = pkin(3) * t245 + pkin(2);
t342 = m(4) * pkin(2) + m(5) * t218 + m(6) * (t215 + t218) + t352 - t354;
t294 = qJ(5) * t150;
t242 = sin(qJ(2));
t280 = t242 * t299;
t186 = pkin(7) * t235 + t280;
t273 = pkin(8) * t235 + t186;
t136 = t273 * t245;
t125 = t244 * t136;
t135 = t273 * t241;
t298 = qJD(3) * pkin(3);
t126 = -t135 + t298;
t76 = t126 * t240 + t125;
t42 = t76 + t294;
t269 = t76 * mrSges(5,3) + t42 * mrSges(6,3);
t187 = -pkin(2) * t235 - t279;
t268 = mrSges(4,1) * t241 + mrSges(4,2) * t245;
t341 = t187 * t268 + qJD(3) * (Ifges(4,5) * t245 - Ifges(4,6) * t241) / 0.2e1;
t340 = -m(4) * pkin(7) - mrSges(6,3) - mrSges(5,3) - mrSges(4,3) + mrSges(3,2) + m(6) * (-qJ(5) + t248) + m(5) * t248;
t258 = t175 * qJD(4);
t109 = t175 * qJD(3) + t258;
t264 = -qJ(5) * t109 - qJD(5) * t176;
t336 = t264 + t335;
t259 = t176 * qJD(4);
t110 = -t176 * qJD(3) - t259;
t285 = t110 * qJ(5) + t175 * qJD(5);
t334 = t285 + t333;
t142 = t151 * qJ(5);
t123 = t240 * t136;
t75 = t244 * t126 - t123;
t41 = -t142 + t75;
t216 = pkin(1) * t242 + pkin(7);
t309 = -pkin(8) - t216;
t169 = t309 * t241;
t291 = t216 * t245;
t170 = t230 + t291;
t104 = t240 * t169 + t244 * t170;
t287 = t235 * t241;
t184 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t287;
t286 = t235 * t245;
t185 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t286;
t332 = (mrSges(3,2) * t235 + t184 * t241 - t185 * t245) * t246 - m(4) * (t187 * t242 + (t241 ^ 2 + t245 ^ 2) * t246 * t186);
t326 = t150 / 0.2e1;
t324 = t151 / 0.2e1;
t319 = t234 / 0.2e1;
t316 = pkin(1) * t246;
t315 = g(1) * t245;
t239 = qJ(1) + qJ(2);
t229 = cos(t239);
t314 = g(3) * t229;
t307 = mrSges(5,3) * t150;
t306 = mrSges(6,3) * t150;
t305 = Ifges(4,4) * t241;
t304 = Ifges(4,4) * t245;
t301 = Ifges(4,2) * t245;
t300 = pkin(1) * qJD(2);
t277 = qJD(2) * t299;
t292 = qJDD(1) * pkin(1);
t173 = t242 * t292 + t246 * t277;
t232 = qJDD(1) + qJDD(2);
t154 = pkin(7) * t232 + t173;
t283 = qJD(3) * t245;
t99 = -t154 * t241 - t186 * t283;
t297 = t241 * t99;
t284 = qJD(3) * t241;
t98 = t245 * t154 - t186 * t284;
t295 = t245 * t98;
t293 = qJ(5) * t176;
t227 = sin(t239);
t290 = t226 * t227;
t78 = -t244 * t135 - t123;
t278 = t246 * t300;
t221 = pkin(3) * t284;
t160 = t232 * t245 - t235 * t284;
t161 = t232 * t241 + t235 * t283;
t68 = t160 * t240 + t161 * t244 + t235 * t258;
t69 = t160 * t244 - t161 * t240 - t235 * t259;
t16 = -t69 * mrSges(6,1) + t68 * mrSges(6,2);
t100 = -pkin(4) * t110 + t221;
t272 = qJD(3) * t309;
t271 = t354 * t235;
t77 = t135 * t240 - t125;
t103 = t244 * t169 - t170 * t240;
t129 = t244 * t199 - t200 * t240;
t270 = -t227 * t275 - t290 * t311;
t172 = -t242 * t277 + t246 * t292;
t267 = t301 + t305;
t265 = t184 * t245 + t185 * t241;
t147 = -pkin(4) * t175 - t218;
t263 = mrSges(2,1) + (m(3) + m(4) + m(5) + m(6)) * pkin(1);
t262 = -t352 + t353;
t260 = t241 * (Ifges(4,1) * t245 - t305);
t71 = qJDD(3) * pkin(3) - pkin(8) * t161 + t99;
t74 = pkin(8) * t160 + t98;
t7 = t126 * t281 - t136 * t282 + t240 * t71 + t244 * t74;
t131 = t241 * t272 + t245 * t278;
t132 = -t241 * t278 + t245 * t272;
t24 = t244 * t131 + t240 * t132 + t169 * t281 - t170 * t282;
t153 = -pkin(2) * t232 - t172;
t257 = m(6) * (-pkin(3) * t241 - pkin(4) * t226) - t268;
t152 = -t218 * t235 - t279;
t101 = -pkin(3) * t160 + t153;
t8 = -t76 * qJD(4) - t240 * t74 + t244 * t71;
t25 = -qJD(4) * t104 - t131 * t240 + t244 * t132;
t254 = -t265 * qJD(3) + t245 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t160) - t241 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t161);
t231 = qJDD(3) + qJDD(4);
t3 = pkin(4) * t231 - qJ(5) * t68 - qJD(5) * t151 + t8;
t36 = pkin(4) * t234 + t41;
t4 = qJ(5) * t69 + qJD(5) * t150 + t7;
t95 = -pkin(4) * t150 + qJD(5) + t152;
t252 = t8 * mrSges(5,1) + t3 * mrSges(6,1) - t7 * mrSges(5,2) - t4 * mrSges(6,2) - t152 * (mrSges(5,1) * t151 + mrSges(5,2) * t150) + t36 * t306 - t95 * (mrSges(6,1) * t151 + mrSges(6,2) * t150) + t75 * t307 + t345 * t69 + t347 * t68 - (t349 * t150 - t350) * t151 / 0.2e1 + t344 * t324 - (t150 * t347 - t151 * t345) * t234 / 0.2e1 + (Ifges(6,3) + Ifges(5,3)) * t231 - (-t346 * t151 + t343 + t351) * t150 / 0.2e1;
t251 = -t227 * t340 + (t342 - t353) * t229;
t250 = t227 * t342 + t229 * t340 - t290 * t310;
t148 = Ifges(4,6) * qJD(3) + t267 * t235;
t204 = Ifges(4,4) * t286;
t149 = Ifges(4,1) * t287 + Ifges(4,5) * qJD(3) + t204;
t23 = -pkin(4) * t69 + qJDD(5) + t101;
t249 = t245 * (Ifges(4,4) * t161 + Ifges(4,2) * t160) / 0.2e1 + (Ifges(4,1) * t161 + Ifges(4,4) * t160) * t318 + (0.2e1 * Ifges(4,5) * t318 + Ifges(4,6) * t245) * qJDD(3) + (t260 * t235 / 0.2e1 + t341) * qJD(3) + Ifges(3,3) * t232 + t153 * t198 + t172 * mrSges(3,1) - t173 * mrSges(3,2) + t161 * (Ifges(4,1) * t241 + t304) / 0.2e1 + t160 * t267 / 0.2e1 + (t235 * (-Ifges(4,2) * t241 + t304) + t149) * t283 / 0.2e1 - t148 * t284 / 0.2e1 + mrSges(4,3) * t295 + (t344 / 0.2e1 - t152 * mrSges(5,1) - t95 * mrSges(6,1) + t345 * t319 + t348 * t324 + t346 * t326 + t269) * t110 + (-t36 * mrSges(6,3) - t75 * mrSges(5,3) + t152 * mrSges(5,2) + t95 * mrSges(6,2) + t347 * t319 + t349 * t324 + t348 * t326 + t343 / 0.2e1) * t109 + (mrSges(5,2) * t101 + mrSges(6,2) * t23 - mrSges(5,3) * t8 - mrSges(6,3) * t3 + t231 * t347 + t348 * t69 + t349 * t68) * t176 + (-mrSges(5,1) * t101 - mrSges(6,1) * t23 + mrSges(5,3) * t7 + mrSges(6,3) * t4 + t231 * t345 + t346 * t69 + t348 * t68) * t175;
t247 = cos(qJ(1));
t243 = sin(qJ(1));
t222 = t242 * t300;
t219 = -pkin(2) - t316;
t197 = -t218 - t316;
t183 = t222 + t221;
t168 = t175 * qJ(5);
t134 = t147 - t316;
t122 = mrSges(5,1) * t234 - mrSges(5,3) * t151;
t121 = mrSges(6,1) * t234 - mrSges(6,3) * t151;
t120 = -mrSges(5,2) * t234 + t307;
t119 = -mrSges(6,2) * t234 + t306;
t111 = pkin(3) * t287 + pkin(4) * t151;
t102 = -mrSges(4,1) * t160 + mrSges(4,2) * t161;
t97 = t168 + t130;
t96 = t129 - t293;
t93 = -mrSges(5,1) * t150 + mrSges(5,2) * t151;
t92 = -mrSges(6,1) * t150 + mrSges(6,2) * t151;
t89 = t100 + t222;
t80 = t168 + t104;
t79 = t103 - t293;
t55 = -mrSges(5,2) * t231 + mrSges(5,3) * t69;
t54 = -mrSges(6,2) * t231 + mrSges(6,3) * t69;
t53 = mrSges(5,1) * t231 - mrSges(5,3) * t68;
t52 = mrSges(6,1) * t231 - mrSges(6,3) * t68;
t45 = -t142 + t78;
t44 = t77 - t294;
t17 = -mrSges(5,1) * t69 + mrSges(5,2) * t68;
t14 = t25 + t264;
t13 = t24 + t285;
t1 = [t219 * t102 + t183 * t93 + t197 * t17 + t254 * t216 + t249 + t134 * t16 + t13 * t119 + t24 * t120 + t14 * t121 + t25 * t122 + t103 * t53 + t104 * t55 + t79 * t52 + t80 * t54 + t89 * t92 + m(4) * (t153 * t219 - t216 * t297 + t98 * t291) + (-mrSges(2,2) * t243 + t263 * t247 + t251) * g(2) + (mrSges(2,2) * t247 + t263 * t243 + t250) * g(3) + m(5) * (t101 * t197 + t103 * t8 + t104 * t7 + t152 * t183 + t24 * t76 + t25 * t75) + m(6) * (t13 * t42 + t134 * t23 + t14 * t36 + t3 * t79 + t4 * t80 + t89 * t95) - mrSges(4,3) * t297 + Ifges(2,3) * qJDD(1) + ((mrSges(3,1) * t246 - mrSges(3,2) * t242) * t232 + m(3) * (t172 * t246 + t173 * t242) + (t271 * t242 - t332) * qJD(2)) * pkin(1); -t218 * t17 + t249 + (-t99 * mrSges(4,3) + t93 * t298) * t241 + t147 * t16 + t130 * t55 + t129 * t53 + t96 * t52 + t97 * t54 + t100 * t92 - pkin(2) * t102 + m(4) * (-pkin(2) * t153 + (t295 - t297) * pkin(7)) + t250 * g(3) + t254 * pkin(7) + t251 * g(2) + t334 * t119 + t333 * t120 + t336 * t121 + t335 * t122 + ((-t271 - t92 - t93) * t242 + t332) * t299 + (t147 * t23 + t3 * t96 + t4 * t97 + (t100 - t280) * t95 + t334 * t42 + t336 * t36) * m(6) + (-t101 * t218 + t129 * t8 + t130 * t7 + t333 * t76 + t335 * t75 + (t221 - t280) * t152) * m(5); (-m(6) * t315 - t93 * t287 + t244 * t53 + (m(6) * t4 + t54 + t55) * t240 + ((m(6) * t42 + t119 + t120) * t244 + (-m(6) * t36 - t121 - t122) * t240) * qJD(4) + (-t315 - t75 * t282 + t76 * t281 + t240 * t7 + t244 * t8 + (-g(2) * t227 - t152 * t235 + t314) * t241) * m(5)) * pkin(3) + Ifges(4,6) * t160 + Ifges(4,5) * t161 + (-m(6) * t215 + t198 + t262) * g(1) + (t257 * t227 + t270) * g(2) + (t311 * t226 - t257 + t275) * t314 - m(5) * (t75 * t77 + t76 * t78) - t45 * t119 - t78 * t120 - t44 * t121 - t77 * t122 - t111 * t92 - t98 * mrSges(4,2) + t99 * mrSges(4,1) + t252 + (t148 * t318 + (t301 * t318 - t260 / 0.2e1) * t235 - (t149 + t204) * t245 / 0.2e1 - t341) * t235 - m(6) * (t111 * t95 + t36 * t44 + t42 * t45) + t265 * t186 + t269 * t151 + Ifges(4,3) * qJDD(3) + (m(6) * t3 + t52) * (pkin(3) * t244 + pkin(4)); (-(-t36 + t41) * t42 + (-g(1) * t228 - g(2) * t290 - t151 * t95 + t3) * pkin(4)) * m(6) + t262 * g(1) - t41 * t119 - t75 * t120 + t42 * t121 + t76 * t122 + pkin(4) * t52 + t252 + (-pkin(4) * t92 + t269) * t151 + t270 * g(2) + (t275 + (m(6) * pkin(4) + t311) * t226) * t314; -t150 * t119 + t151 * t121 + (-g(2) * t229 - g(3) * t227 - t150 * t42 + t36 * t151 + t23) * m(6) + t16;];
tau = t1;
