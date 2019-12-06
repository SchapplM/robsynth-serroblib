% Calculate vector of inverse dynamics joint torques for
% S5RRRRR2
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-12-05 18:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_invdynJ_fixb_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:52:56
% EndTime: 2019-12-05 18:53:08
% DurationCPUTime: 3.99s
% Computational Cost: add. (4529->432), mult. (8540->625), div. (0->0), fcn. (5813->14), ass. (0->213)
t209 = sin(qJ(3));
t329 = t209 / 0.2e1;
t208 = sin(qJ(4));
t213 = cos(qJ(4));
t210 = sin(qJ(2));
t306 = pkin(1) * qJD(1);
t266 = t210 * t306;
t228 = qJD(3) * pkin(2) - t209 * t266;
t214 = cos(qJ(3));
t256 = t214 * t266;
t125 = t208 * t228 + t213 * t256;
t215 = cos(qJ(2));
t277 = qJD(1) * t215;
t265 = pkin(1) * t277;
t202 = qJD(1) + qJD(2);
t282 = t202 * t214;
t157 = -pkin(2) * t282 - t265;
t207 = sin(qJ(5));
t212 = cos(qJ(5));
t90 = -t125 * t207 + t157 * t212;
t356 = t90 * mrSges(6,1);
t91 = t125 * t212 + t157 * t207;
t355 = t91 * mrSges(6,2);
t354 = mrSges(6,1) * t212 + mrSges(5,1);
t317 = mrSges(4,2) * t209;
t319 = mrSges(4,1) * t214;
t178 = t317 - t319;
t353 = t178 - mrSges(3,1);
t281 = t208 * t209;
t165 = -t213 * t214 + t281;
t229 = t165 * qJD(4);
t116 = -qJD(3) * t165 - t229;
t280 = t208 * t214;
t166 = t209 * t213 + t280;
t270 = qJD(5) * t212;
t236 = t116 * t207 + t166 * t270;
t276 = qJD(2) * t215;
t261 = qJD(1) * t276;
t161 = (qJDD(1) * t210 + t261) * pkin(1);
t255 = qJD(3) * t266;
t122 = t161 * t214 - t209 * t255;
t123 = -t161 * t209 - t214 * t255;
t352 = t122 * t214 - t123 * t209;
t200 = qJDD(1) + qJDD(2);
t275 = qJD(3) * t209;
t152 = t200 * t214 - t202 * t275;
t326 = pkin(1) * t210;
t264 = qJD(2) * t326;
t325 = pkin(1) * t215;
t160 = qJD(1) * t264 - qJDD(1) * t325;
t113 = -pkin(2) * t152 + t160;
t224 = qJDD(3) * pkin(2) + t123;
t124 = t208 * t256 - t213 * t228;
t273 = qJD(4) * t124;
t52 = t213 * t122 + t208 * t224 - t273;
t12 = qJD(5) * t90 + t113 * t207 + t212 * t52;
t13 = -qJD(5) * t91 + t113 * t212 - t207 * t52;
t351 = t12 * t212 - t13 * t207;
t206 = qJ(1) + qJ(2);
t196 = sin(t206);
t198 = cos(t206);
t254 = g(1) * t198 + g(2) * t196;
t343 = t166 * qJD(4);
t117 = qJD(3) * t166 + t343;
t350 = qJD(3) * (Ifges(4,5) * t214 - Ifges(4,6) * t209);
t349 = mrSges(3,2) - mrSges(5,3) - mrSges(4,3);
t205 = qJ(3) + qJ(4);
t195 = sin(t205);
t316 = mrSges(5,2) * t195;
t348 = t316 + t317 - mrSges(3,1);
t274 = qJD(3) * t214;
t153 = t200 * t209 + t202 * t274;
t81 = t152 * t213 - t153 * t208 - t202 * t343;
t347 = t13 * mrSges(6,1) - t12 * mrSges(6,2);
t150 = t166 * t202;
t201 = qJD(3) + qJD(4);
t114 = -t150 * t207 + t201 * t212;
t199 = qJDD(3) + qJDD(4);
t80 = t152 * t208 + t153 * t213 - t202 * t229;
t44 = qJD(5) * t114 + t199 * t207 + t212 * t80;
t341 = t44 / 0.2e1;
t115 = t150 * t212 + t201 * t207;
t45 = -qJD(5) * t115 + t199 * t212 - t207 * t80;
t340 = t45 / 0.2e1;
t149 = t165 * t202;
t138 = qJD(5) + t149;
t310 = Ifges(6,4) * t115;
t58 = Ifges(6,2) * t114 + Ifges(6,6) * t138 + t310;
t339 = -t58 / 0.2e1;
t79 = qJDD(5) - t81;
t338 = t79 / 0.2e1;
t184 = t195 * mrSges(6,3);
t197 = cos(t205);
t346 = t197 * mrSges(5,1) + t184;
t345 = -t207 * t90 + t212 * t91;
t342 = Ifges(6,1) * t341 + Ifges(6,4) * t340 + Ifges(6,5) * t338;
t337 = -m(5) - m(6);
t336 = -t114 / 0.2e1;
t335 = -t115 / 0.2e1;
t334 = t115 / 0.2e1;
t333 = -t138 / 0.2e1;
t330 = t150 / 0.2e1;
t328 = t212 / 0.2e1;
t324 = pkin(2) * t214;
t320 = -mrSges(5,1) * t199 - mrSges(6,1) * t45 + mrSges(6,2) * t44 + mrSges(5,3) * t80;
t315 = mrSges(6,2) * t207;
t314 = mrSges(5,3) * t149;
t313 = Ifges(4,4) * t209;
t312 = Ifges(4,4) * t214;
t311 = Ifges(5,4) * t150;
t309 = Ifges(6,4) * t207;
t308 = Ifges(6,4) * t212;
t307 = Ifges(4,2) * t214;
t303 = t150 * mrSges(5,3);
t272 = qJD(4) * t125;
t53 = t208 * t122 - t213 * t224 + t272;
t302 = t166 * t53;
t299 = mrSges(5,1) * t201 + mrSges(6,1) * t114 - mrSges(6,2) * t115 - t303;
t133 = t165 * t266;
t295 = t124 * t133;
t234 = pkin(1) * t166;
t134 = t234 * t277;
t294 = t124 * t134;
t293 = t149 * t207;
t292 = t149 * t212;
t291 = t166 * t207;
t290 = t166 * t212;
t289 = t197 * t198;
t288 = t197 * t207;
t287 = t197 * t212;
t286 = t198 * t207;
t285 = t198 * t212;
t284 = t198 * t214;
t283 = t202 * t209;
t267 = t195 * t315;
t279 = (mrSges(6,3) * t197 + t267) * t196;
t278 = mrSges(6,3) * t289 + t198 * t267;
t271 = qJD(5) * t207;
t269 = Ifges(6,5) * t44 + Ifges(6,6) * t45 + Ifges(6,3) * t79;
t268 = pkin(2) * t283;
t111 = Ifges(6,4) * t114;
t59 = Ifges(6,1) * t115 + Ifges(6,5) * t138 + t111;
t263 = t59 * t328;
t259 = -t271 / 0.2e1;
t258 = t353 * t202;
t253 = g(1) * t196 - g(2) * t198;
t211 = sin(qJ(1));
t216 = cos(qJ(1));
t252 = g(1) * t211 - g(2) * t216;
t251 = mrSges(4,1) * t209 + mrSges(4,2) * t214;
t250 = mrSges(6,1) * t207 + mrSges(6,2) * t212;
t249 = Ifges(6,1) * t212 - t309;
t248 = t307 + t313;
t247 = -Ifges(6,2) * t207 + t308;
t245 = Ifges(6,5) * t212 - Ifges(6,6) * t207;
t147 = t210 * t234;
t72 = (t276 * t280 + (t201 * t210 * t214 + t209 * t276) * t213) * pkin(1) - t201 * t281 * t326;
t244 = t124 * t72 + t147 * t53;
t84 = -mrSges(6,2) * t138 + mrSges(6,3) * t114;
t85 = mrSges(6,1) * t138 - mrSges(6,3) * t115;
t243 = t207 * t85 - t212 * t84;
t242 = -t207 * t84 - t212 * t85;
t241 = t207 * t91 + t212 * t90;
t148 = t165 * t326;
t179 = -t324 - t325;
t105 = -t148 * t212 + t179 * t207;
t104 = t148 * t207 + t179 * t212;
t169 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t283;
t170 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t282;
t239 = t169 * t214 + t170 * t209;
t238 = m(4) * t210 * pkin(1) ^ 2 * (t209 ^ 2 + t214 ^ 2 - 0.1e1);
t120 = -mrSges(5,2) * t201 - t314;
t237 = t120 - t243;
t235 = -t116 * t212 + t166 * t271;
t233 = t209 * (Ifges(4,1) * t214 - t313);
t232 = t165 * t215;
t227 = -mrSges(3,2) * t202 - t169 * t209 + t170 * t214;
t226 = mrSges(5,2) * t197 + t195 * t354;
t225 = -mrSges(6,1) * t287 + mrSges(6,2) * t288 + t316 - t346;
t141 = t196 * t212 - t197 * t286;
t142 = t196 * t207 + t197 * t285;
t222 = -mrSges(4,1) * t284 - mrSges(5,1) * t289 - t142 * mrSges(6,1) - t141 * mrSges(6,2) + (-t184 + t348) * t198 + t349 * t196;
t139 = t196 * t288 + t285;
t140 = -t196 * t287 + t286;
t221 = -t140 * mrSges(6,1) - t139 * mrSges(6,2) + t349 * t198 + (t319 + t346 - t348) * t196;
t137 = Ifges(5,4) * t149;
t57 = Ifges(6,5) * t115 + Ifges(6,6) * t114 + Ifges(6,3) * t138;
t6 = Ifges(6,4) * t44 + Ifges(6,2) * t45 + Ifges(6,6) * t79;
t87 = -Ifges(5,2) * t149 + Ifges(5,6) * t201 + t311;
t88 = Ifges(5,1) * t150 + Ifges(5,5) * t201 - t137;
t220 = -t150 * t356 + (t315 - t354) * t53 + t6 * t328 + t87 * t330 + Ifges(5,5) * t80 + Ifges(5,6) * t81 - t52 * mrSges(5,2) + (-qJD(5) * t241 - t292 * t90 - t293 * t91 + t351) * mrSges(6,3) + (Ifges(6,5) * t207 + Ifges(6,6) * t212) * t338 + (Ifges(6,2) * t212 + t309) * t340 + (Ifges(6,1) * t207 + t308) * t341 + t207 * t342 + (t114 * t247 + t115 * t249 + t138 * t245) * qJD(5) / 0.2e1 - (-Ifges(5,1) * t149 - t311 + t57) * t150 / 0.2e1 + (-Ifges(5,2) * t150 - t137 + t88) * t149 / 0.2e1 + t293 * t339 + t150 * t355 + t124 * t314 + Ifges(5,3) * t199 + t58 * t259 + (t124 * t250 + t263) * qJD(5) + (Ifges(6,3) * t150 - t149 * t245) * t333 + (Ifges(6,5) * t150 - t149 * t249) * t335 + (Ifges(6,6) * t150 - t149 * t247) * t336 - t201 * (-Ifges(5,5) * t149 - Ifges(5,6) * t150) / 0.2e1 - t157 * (mrSges(5,1) * t150 - mrSges(5,2) * t149) + t59 * t292 / 0.2e1;
t145 = Ifges(4,6) * qJD(3) + t202 * t248;
t183 = Ifges(4,4) * t282;
t146 = Ifges(4,1) * t283 + Ifges(4,5) * qJD(3) + t183;
t219 = -t117 * t355 + t353 * t160 + (Ifges(5,1) * t116 - Ifges(5,4) * t117) * t330 + (-Ifges(6,1) * t235 - Ifges(6,4) * t236 + Ifges(6,5) * t117) * t334 - t161 * mrSges(3,2) + t157 * (mrSges(5,1) * t117 + mrSges(5,2) * t116) - t149 * (Ifges(5,4) * t116 - Ifges(5,2) * t117) / 0.2e1 + t116 * t88 / 0.2e1 + t117 * t57 / 0.2e1 - t117 * t87 / 0.2e1 + t352 * mrSges(4,3) + t236 * t339 + t214 * (Ifges(4,4) * t153 + Ifges(4,2) * t152) / 0.2e1 + t290 * t342 + (t233 * t202 / 0.2e1 + t350 / 0.2e1) * qJD(3) + (0.2e1 * Ifges(4,5) * t329 + Ifges(4,6) * t214) * qJDD(3) + (t116 * t124 - t125 * t117 + t302) * mrSges(5,3) + (mrSges(5,2) * t113 + Ifges(5,1) * t80 + Ifges(5,4) * t81 + Ifges(5,5) * t199 + t245 * t338 + t247 * t340 + t249 * t341 + t259 * t59) * t166 + (t202 * (-Ifges(4,2) * t209 + t312) + t146) * t274 / 0.2e1 + (t113 * mrSges(5,1) - Ifges(5,2) * t81 - Ifges(5,4) * t80 + Ifges(6,3) * t338 + Ifges(6,6) * t340 + Ifges(6,5) * t341 - t52 * mrSges(5,3) - Ifges(5,6) * t199 + t269 / 0.2e1 + t347) * t165 + (-t12 * t291 - t13 * t290 + t90 * t235 - t91 * t236) * mrSges(6,3) + t117 * t356 + (Ifges(4,1) * t153 + Ifges(4,4) * t152) * t329 + t138 * (-Ifges(6,5) * t235 - Ifges(6,6) * t236 + Ifges(6,3) * t117) / 0.2e1 + t114 * (-Ifges(6,4) * t235 - Ifges(6,2) * t236 + Ifges(6,6) * t117) / 0.2e1 + t124 * (mrSges(6,1) * t236 - mrSges(6,2) * t235) + t250 * t302 + Ifges(3,3) * t200 + t201 * (Ifges(5,5) * t116 - Ifges(5,6) * t117) / 0.2e1 + t152 * t248 / 0.2e1 + t116 * t263 - t145 * t275 / 0.2e1 - t6 * t291 / 0.2e1 + t153 * (Ifges(4,1) * t209 + t312) / 0.2e1;
t168 = pkin(2) * t275 + t264;
t167 = t251 * qJD(3);
t136 = t232 * t306;
t135 = t166 * t266;
t109 = -t136 * t212 + t207 * t266;
t108 = t136 * t207 + t212 * t266;
t103 = -t135 * t212 + t207 * t268;
t102 = t135 * t207 + t212 * t268;
t100 = mrSges(5,1) * t149 + mrSges(5,2) * t150;
t89 = t250 * t149;
t71 = (-qJD(2) * t232 - t117 * t210) * pkin(1);
t70 = -mrSges(5,2) * t199 + mrSges(5,3) * t81;
t29 = -qJD(5) * t105 + t168 * t212 - t207 * t71;
t28 = qJD(5) * t104 + t168 * t207 + t212 * t71;
t27 = -mrSges(5,1) * t81 + mrSges(5,2) * t80;
t15 = -mrSges(6,2) * t79 + mrSges(6,3) * t45;
t14 = mrSges(6,1) * t79 - mrSges(6,3) * t44;
t1 = [m(5) * (t113 * t179 + t125 * t71 - t148 * t52 + t157 * t168 + t244) + m(6) * (t104 * t13 + t105 * t12 + t28 * t91 + t29 * t90 + t244) + t168 * t100 + t179 * t27 - t148 * t70 + t71 * t120 + t104 * t14 + t105 * t15 + t28 * t84 + t29 * t85 - t299 * t72 + t320 * t147 + t219 + (t252 * m(4) + (-m(4) * t160 + t200 * mrSges(3,1) + mrSges(4,1) * t152 - mrSges(4,2) * t153 - qJD(1) * t167) * t215 + (m(4) * t352 - t200 * mrSges(3,2) - t209 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t153) + t214 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t152) - t239 * qJD(3)) * t210 + (-t160 * t215 + t161 * t210 + t252) * m(3) + (t210 * t258 + t215 * t227) * qJD(2)) * pkin(1) + (mrSges(2,1) * t211 + mrSges(2,2) * t216 + t337 * (-pkin(1) * t211 - t196 * t324) + t221) * g(1) + (-mrSges(2,1) * t216 + mrSges(2,2) * t211 + t337 * (pkin(1) * t216 + pkin(2) * t284) + t222) * g(2) + t238 * t261 + Ifges(2,3) * qJDD(1); t136 * t120 - t108 * t85 - t109 * t84 + (-t100 - t258) * t266 + t221 * g(1) + t222 * g(2) + t299 * t134 + t219 + (-qJD(1) ^ 2 * t238 + (-t167 - t227) * t306) * t215 + ((m(5) * t157 + m(6) * t241 + t100 - t242) * t275 + (-t212 * t14 - t207 * t15 - t27 + t243 * qJD(5) + (-t12 * t207 - t13 * t212 - t270 * t91 + t271 * t90 + t253) * m(6) + (-t113 + t253) * m(5)) * t214) * pkin(2) - m(6) * (t108 * t90 + t109 * t91 + t294) - m(5) * (-t125 * t136 + t157 * t266 + t294); -m(5) * (-t125 * t135 - t295) + Ifges(4,6) * t152 + Ifges(4,5) * t153 + t135 * t120 + t124 * t89 - t122 * mrSges(4,2) + t123 * mrSges(4,1) - t102 * t85 - t103 * t84 - g(2) * t279 - g(1) * t278 + (t145 * t329 - t350 / 0.2e1 + t251 * t265 + (t307 * t329 - t233 / 0.2e1) * t202 - (t146 + t183) * t214 / 0.2e1) * t202 - t299 * t133 + (t178 + t225) * g(3) + t125 * t303 - m(6) * (t102 * t90 + t103 * t91 - t295) + t220 + (t337 * t214 * g(3) + (-t202 * t100 + t254 * m(6) + (-t157 * t202 + t254) * m(5)) * t209 + (m(5) * (-t53 + t272) - m(6) * t53 - t320 + (m(6) * t345 + t237) * qJD(4)) * t213 + (-t207 * t14 + t212 * t15 + t70 + t242 * qJD(5) - t299 * qJD(4) + m(5) * (t52 + t273) + m(6) * (-t270 * t90 - t271 * t91 + t273 + t351)) * t208) * pkin(2) + t239 * t266 + Ifges(4,3) * qJDD(3) + t254 * (t226 + t251); (t299 + t303) * t125 + t220 + (t89 - m(6) * (t125 - t345) + t237) * t124 + (t198 * t226 - t278) * g(1) + (t196 * t226 - t279) * g(2) + t225 * g(3); -t124 * (mrSges(6,1) * t115 + mrSges(6,2) * t114) + (Ifges(6,1) * t114 - t310) * t335 + t58 * t334 + (Ifges(6,5) * t114 - Ifges(6,6) * t115) * t333 - t90 * t84 + t91 * t85 - g(1) * (mrSges(6,1) * t141 - mrSges(6,2) * t142) - g(2) * (-mrSges(6,1) * t139 + mrSges(6,2) * t140) + g(3) * t250 * t195 + (t114 * t90 + t115 * t91) * mrSges(6,3) + t269 + (-Ifges(6,2) * t115 + t111 + t59) * t336 + t347;];
tau = t1;
