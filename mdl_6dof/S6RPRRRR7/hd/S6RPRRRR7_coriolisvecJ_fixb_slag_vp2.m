% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR7_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:16:57
% EndTime: 2019-03-09 07:17:11
% DurationCPUTime: 6.26s
% Computational Cost: add. (13466->510), mult. (29342->710), div. (0->0), fcn. (20306->8), ass. (0->244)
t228 = sin(qJ(4));
t229 = sin(qJ(3));
t232 = cos(qJ(4));
t233 = cos(qJ(3));
t193 = -t228 * t233 - t232 * t229;
t227 = sin(qJ(5));
t231 = cos(qJ(5));
t291 = t232 * t233;
t255 = t228 * t229 - t291;
t150 = t193 * t227 - t231 * t255;
t223 = qJD(3) + qJD(4);
t287 = qJD(4) * t228;
t288 = qJD(3) * t229;
t151 = t223 * t291 - t228 * t288 - t229 * t287;
t152 = t223 * t193;
t243 = qJD(5) * t150 + t231 * t151 + t152 * t227;
t234 = -pkin(1) - pkin(7);
t204 = qJD(1) * t234 + qJD(2);
t290 = qJD(1) * t229;
t170 = -pkin(8) * t290 + t204 * t229;
t164 = t232 * t170;
t289 = qJD(1) * t233;
t171 = -pkin(8) * t289 + t233 * t204;
t314 = qJD(3) * pkin(3);
t166 = t171 + t314;
t121 = t166 * t228 + t164;
t185 = t193 * qJD(1);
t335 = pkin(9) * t185;
t106 = t121 + t335;
t302 = t106 * t227;
t163 = t228 * t170;
t120 = t232 * t166 - t163;
t186 = -t228 * t290 + t232 * t289;
t175 = t186 * pkin(9);
t105 = t120 - t175;
t97 = pkin(4) * t223 + t105;
t52 = t231 * t97 - t302;
t292 = t231 * t106;
t53 = t227 * t97 + t292;
t257 = t231 * t193 + t227 * t255;
t74 = qJD(5) * t257 - t151 * t227 + t231 * t152;
t365 = t243 * t53 + t52 * t74;
t141 = t152 * qJD(1);
t333 = t141 * pkin(9);
t142 = t223 * t255 * qJD(1);
t278 = pkin(8) * qJD(1) - t204;
t260 = t278 * t229;
t261 = t278 * t233;
t286 = qJD(4) * t232;
t81 = t166 * t286 - t170 * t287 + (t228 * t260 - t232 * t261) * qJD(3);
t55 = pkin(9) * t142 + t81;
t241 = (t228 * t261 + t232 * t260) * qJD(3);
t82 = -qJD(4) * t121 + t241;
t13 = t227 * t55 - t231 * (t82 - t333) + t53 * qJD(5);
t221 = qJD(5) + t223;
t226 = sin(qJ(6));
t230 = cos(qJ(6));
t258 = t185 * t227 + t231 * t186;
t114 = t221 * t230 - t226 * t258;
t276 = t231 * t185 - t186 * t227;
t68 = qJD(5) * t276 + t141 * t231 + t142 * t227;
t41 = qJD(6) * t114 + t230 * t68;
t115 = t221 * t226 + t230 * t258;
t42 = -qJD(6) * t115 - t226 * t68;
t14 = -mrSges(7,1) * t42 + mrSges(7,2) * t41;
t364 = m(7) * t13 + t14;
t130 = qJD(6) - t276;
t272 = mrSges(7,1) * t226 + mrSges(7,2) * t230;
t50 = -pkin(5) * t221 - t52;
t252 = t50 * t272;
t267 = Ifges(7,5) * t230 - Ifges(7,6) * t226;
t321 = Ifges(7,4) * t230;
t269 = -Ifges(7,2) * t226 + t321;
t322 = Ifges(7,4) * t226;
t271 = Ifges(7,1) * t230 - t322;
t338 = t230 / 0.2e1;
t340 = -t226 / 0.2e1;
t346 = t115 / 0.2e1;
t313 = t115 * Ifges(7,4);
t48 = t114 * Ifges(7,2) + t130 * Ifges(7,6) + t313;
t113 = Ifges(7,4) * t114;
t49 = t115 * Ifges(7,1) + t130 * Ifges(7,5) + t113;
t363 = t114 * t269 / 0.2e1 + t271 * t346 + t130 * t267 / 0.2e1 + t252 + t48 * t340 + t49 * t338;
t323 = Ifges(6,4) * t258;
t129 = Ifges(6,4) * t276;
t307 = mrSges(6,1) * t221 + mrSges(7,1) * t114 - mrSges(7,2) * t115 - mrSges(6,3) * t258;
t125 = t232 * t171 - t163;
t107 = -t175 + t125;
t216 = pkin(3) * t232 + pkin(4);
t124 = -t171 * t228 - t164;
t254 = t124 - t335;
t284 = qJD(5) * t231;
t285 = qJD(5) * t227;
t294 = t228 * t231;
t359 = t107 * t227 - t231 * t254 - t216 * t285 - (t228 * t284 + (t227 * t232 + t294) * qJD(4)) * pkin(3);
t295 = t227 * t228;
t139 = t216 * t284 + (-t228 * t285 + (t231 * t232 - t295) * qJD(4)) * pkin(3);
t62 = t231 * t107 + t227 * t254;
t358 = -t62 + t139;
t357 = (qJ(2) * (m(3) + m(4)));
t329 = pkin(8) - t234;
t198 = t329 * t229;
t199 = t329 * t233;
t154 = -t232 * t198 - t228 * t199;
t51 = pkin(10) * t221 + t53;
t202 = pkin(3) * t290 + qJD(1) * qJ(2);
t155 = -pkin(4) * t185 + t202;
t70 = -pkin(5) * t276 - pkin(10) * t258 + t155;
t20 = -t226 * t51 + t230 * t70;
t21 = t226 * t70 + t230 * t51;
t356 = -t20 * t226 + t21 * t230;
t12 = t52 * qJD(5) + t231 * t55 + (-t166 * t287 - t170 * t286 + t241 - t333) * t227;
t219 = pkin(3) * t289;
t197 = qJD(1) * qJD(2) + qJD(3) * t219;
t116 = -pkin(4) * t142 + t197;
t69 = qJD(5) * t258 + t141 * t227 - t231 * t142;
t19 = pkin(5) * t69 - pkin(10) * t68 + t116;
t2 = qJD(6) * t20 + t12 * t230 + t19 * t226;
t303 = qJD(6) * t21;
t3 = -t12 * t226 + t19 * t230 - t303;
t355 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t41 + Ifges(7,6) * t42;
t96 = pkin(5) * t258 - pkin(10) * t276;
t264 = t20 * t230 + t21 * t226;
t320 = Ifges(6,5) * t221;
t327 = Ifges(6,1) * t258;
t90 = t129 + t320 + t327;
t354 = t155 * mrSges(6,2) + t90 / 0.2e1 + t320 / 0.2e1 + t129 / 0.2e1 - t264 * mrSges(7,3) + t363;
t353 = t41 / 0.2e1;
t352 = t42 / 0.2e1;
t351 = t69 / 0.2e1;
t348 = -t114 / 0.2e1;
t347 = -t115 / 0.2e1;
t345 = -t130 / 0.2e1;
t344 = t152 / 0.2e1;
t343 = -t185 / 0.2e1;
t342 = -t186 / 0.2e1;
t341 = t186 / 0.2e1;
t339 = t229 / 0.2e1;
t337 = -t233 / 0.2e1;
t336 = pkin(4) * t186;
t153 = t198 * t228 - t232 * t199;
t122 = pkin(9) * t255 + t153;
t123 = pkin(9) * t193 + t154;
t259 = t231 * t122 - t123 * t227;
t334 = t13 * t259;
t332 = t2 * t230;
t331 = t226 * t3;
t330 = t52 * mrSges(6,3);
t328 = mrSges(5,3) * t185;
t326 = Ifges(4,4) * t229;
t325 = Ifges(4,4) * t233;
t324 = Ifges(5,4) * t186;
t319 = Ifges(7,5) * t115;
t318 = Ifges(6,2) * t276;
t317 = Ifges(6,6) * t221;
t316 = Ifges(7,6) * t114;
t315 = Ifges(7,3) * t130;
t312 = t13 * t150;
t311 = t258 * t53;
t310 = t186 * mrSges(5,3);
t305 = Ifges(4,5) * qJD(3);
t304 = Ifges(4,6) * qJD(3);
t301 = t276 * t226;
t300 = t276 * t230;
t298 = t193 * t142;
t296 = t255 * t141;
t211 = t229 * pkin(3) + qJ(2);
t180 = pkin(3) * t294 + t227 * t216;
t283 = qJD(6) * t226;
t282 = qJD(6) * t230;
t205 = t233 * t314 + qJD(2);
t144 = -mrSges(5,1) * t185 + mrSges(5,2) * t186;
t279 = -m(5) * t202 - t144;
t167 = -pkin(4) * t193 + t211;
t131 = pkin(4) * t151 + t205;
t274 = mrSges(4,1) * t229 + mrSges(4,2) * t233;
t273 = mrSges(7,1) * t230 - mrSges(7,2) * t226;
t270 = Ifges(7,1) * t226 + t321;
t268 = Ifges(7,2) * t230 + t322;
t266 = Ifges(7,5) * t226 + Ifges(7,6) * t230;
t15 = mrSges(7,1) * t69 - mrSges(7,3) * t41;
t16 = -mrSges(7,2) * t69 + mrSges(7,3) * t42;
t265 = -t226 * t15 + t230 * t16;
t83 = -mrSges(7,2) * t130 + mrSges(7,3) * t114;
t84 = mrSges(7,1) * t130 - mrSges(7,3) * t115;
t262 = -t226 * t83 - t230 * t84;
t86 = t122 * t227 + t123 * t231;
t87 = -pkin(5) * t257 - pkin(10) * t150 + t167;
t32 = t226 * t87 + t230 * t86;
t31 = -t226 * t86 + t230 * t87;
t200 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t290;
t201 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t289;
t256 = t233 * t200 - t229 * t201;
t179 = -pkin(3) * t295 + t216 * t231;
t117 = -mrSges(6,2) * t221 + mrSges(6,3) * t276;
t253 = -t226 * t84 + t230 * t83 + t117;
t251 = qJ(2) * (mrSges(4,1) * t233 - mrSges(4,2) * t229);
t191 = t329 * t288;
t192 = qJD(3) * t199;
t99 = t228 * t191 - t232 * t192 + t198 * t287 - t199 * t286;
t79 = t336 + t96;
t247 = -qJD(6) * t264 - t331;
t246 = t120 * t152 + t121 * t151 - t193 * t81 - t255 * t82;
t245 = t247 * mrSges(7,3);
t100 = -qJD(4) * t154 + t232 * t191 + t192 * t228;
t244 = t247 + t332;
t242 = -pkin(9) * t152 + t100;
t8 = t41 * Ifges(7,4) + t42 * Ifges(7,2) + t69 * Ifges(7,6);
t9 = t41 * Ifges(7,1) + t42 * Ifges(7,4) + t69 * Ifges(7,5);
t240 = -t12 * mrSges(6,2) + mrSges(7,3) * t332 + t270 * t353 + t268 * t352 + t266 * t351 + t226 * t9 / 0.2e1 - Ifges(6,6) * t69 + Ifges(6,5) * t68 + t8 * t338 + (-t273 - mrSges(6,1)) * t13 + t363 * qJD(6);
t47 = t315 + t316 + t319;
t89 = t317 + t318 + t323;
t239 = t155 * mrSges(6,1) + t20 * mrSges(7,1) + t47 / 0.2e1 - t89 / 0.2e1 - t323 / 0.2e1 + t319 / 0.2e1 - t317 / 0.2e1 + t316 / 0.2e1 + t315 / 0.2e1 - t21 * mrSges(7,2);
t238 = m(7) * (-t20 * t282 - t21 * t283 - t331 + t332) - t83 * t283 - t84 * t282 + t265;
t127 = t185 * Ifges(5,2) + t223 * Ifges(5,6) + t324;
t174 = Ifges(5,4) * t185;
t128 = t186 * Ifges(5,1) + t223 * Ifges(5,5) + t174;
t236 = -t81 * mrSges(5,2) + t82 * mrSges(5,1) + (-Ifges(5,2) * t186 + t128 + t174) * t343 + t127 * t341 + (Ifges(5,1) * t185 - t324) * t342 + t120 * t328 - t21 * (-mrSges(7,2) * t258 - mrSges(7,3) * t301) - t20 * (mrSges(7,1) * t258 - mrSges(7,3) * t300) + t258 * t89 / 0.2e1 + t276 * t330 + (Ifges(7,3) * t258 + t267 * t276) * t345 + (Ifges(7,5) * t258 + t271 * t276) * t347 + (Ifges(7,6) * t258 + t269 * t276) * t348 - t155 * (mrSges(6,1) * t258 + mrSges(6,2) * t276) - t221 * (Ifges(6,5) * t276 - Ifges(6,6) * t258) / 0.2e1 - t202 * (mrSges(5,1) * t186 + mrSges(5,2) * t185) - t276 * t252 + t240 + Ifges(5,5) * t141 + Ifges(5,6) * t142 - t223 * (Ifges(5,5) * t185 - Ifges(5,6) * t186) / 0.2e1 - t49 * t300 / 0.2e1 + t48 * t301 / 0.2e1 - (-Ifges(6,2) * t258 + t129 + t90) * t276 / 0.2e1 - (Ifges(6,1) * t276 - t323 + t47) * t258 / 0.2e1;
t195 = t274 * qJD(1);
t184 = t305 + (Ifges(4,1) * t233 - t326) * qJD(1);
t183 = t304 + (-Ifges(4,2) * t229 + t325) * qJD(1);
t173 = pkin(10) + t180;
t172 = -pkin(5) - t179;
t162 = mrSges(5,1) * t223 - t310;
t161 = -mrSges(5,2) * t223 + t328;
t160 = t219 + t336;
t95 = -mrSges(6,1) * t276 + mrSges(6,2) * t258;
t80 = -pkin(9) * t151 + t99;
t72 = t219 + t79;
t65 = Ifges(7,3) * t69;
t57 = t105 * t231 - t302;
t56 = t105 * t227 + t292;
t28 = t226 * t96 + t230 * t52;
t27 = -t226 * t52 + t230 * t96;
t26 = t226 * t72 + t230 * t62;
t25 = -t226 * t62 + t230 * t72;
t24 = t226 * t79 + t230 * t57;
t23 = -t226 * t57 + t230 * t79;
t22 = pkin(5) * t243 - pkin(10) * t74 + t131;
t18 = qJD(5) * t86 + t227 * t80 - t231 * t242;
t17 = qJD(5) * t259 + t227 * t242 + t231 * t80;
t5 = -qJD(6) * t32 - t17 * t226 + t22 * t230;
t4 = qJD(6) * t31 + t17 * t230 + t22 * t226;
t1 = [(t327 / 0.2e1 + t354) * t74 + (-t141 * t153 + t142 * t154 - t246) * mrSges(5,3) + t4 * t83 + t5 * t84 + (t8 * t340 + Ifges(6,1) * t68 - Ifges(6,4) * t69 + t116 * mrSges(6,2) + t9 * t338 + t271 * t353 + t269 * t352 + t267 * t351 + (mrSges(6,3) + t272) * t13 + (-t2 * t226 - t230 * t3) * mrSges(7,3) + (t49 * t340 + t268 * t348 + t270 * t347 + t50 * t273 + t266 * t345 - t230 * t48 / 0.2e1 - t356 * mrSges(7,3)) * qJD(6)) * t150 + (-t259 * t68 - t69 * t86 - t365) * mrSges(6,3) + t128 * t344 + m(5) * (t100 * t120 + t121 * t99 + t153 * t82 + t154 * t81 + t197 * t211 + t202 * t205) + t31 * t15 + t32 * t16 + (t152 * t341 - t296) * Ifges(5,1) + m(7) * (t18 * t50 + t2 * t32 + t20 * t5 + t21 * t4 + t3 * t31 - t334) + m(6) * (t116 * t167 + t12 * t86 + t131 * t155 + t17 * t53 - t18 * t52 - t334) + (-t318 / 0.2e1 + t239) * t243 + t17 * t117 + t131 * t95 - t151 * t127 / 0.2e1 + t99 * t161 + t100 * t162 + t167 * (mrSges(6,1) * t69 + mrSges(6,2) * t68) + qJD(2) * t195 + t202 * (mrSges(5,1) * t151 + mrSges(5,2) * t152) + (((2 * mrSges(3,3)) + t274 + (2 * t357)) * qJD(2) + (0.2e1 * t251 + (-0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(4,2)) * t229 * t233 + (-0.3e1 / 0.2e1 * t233 ^ 2 + 0.3e1 / 0.2e1 * t229 ^ 2) * Ifges(4,4)) * qJD(3)) * qJD(1) + t205 * t144 + t211 * (-mrSges(5,1) * t142 + mrSges(5,2) * t141) + t223 * (Ifges(5,5) * t152 - Ifges(5,6) * t151) / 0.2e1 + ((-t183 / 0.2e1 + t234 * t200 - t304 / 0.2e1) * t233 + (-t184 / 0.2e1 - t234 * t201 - t305 / 0.2e1) * t229) * qJD(3) - t307 * t18 - t259 * t14 - (-t12 * mrSges(6,3) + t65 / 0.2e1 - Ifges(6,4) * t68 + t116 * mrSges(6,1) + (Ifges(6,2) + Ifges(7,3) / 0.2e1) * t69 + t355) * t257 + t197 * (-mrSges(5,1) * t193 - mrSges(5,2) * t255) + (t193 * t141 - t142 * t255 + t151 * t342 + t185 * t344) * Ifges(5,4) + (t151 * t343 + t298) * Ifges(5,2); t151 * t161 + t152 * t162 + t307 * t74 - (t68 * mrSges(6,3) + t14) * t150 + t256 * qJD(3) + (t296 - t298) * mrSges(5,3) + t253 * t243 - (-t69 * mrSges(6,3) + qJD(6) * t262 + t265) * t257 + m(6) * (-t12 * t257 - t312 + t365) + m(5) * t246 + m(7) * (t243 * t356 - t244 * t257 - t50 * t74 - t312) + (-m(6) * t155 - m(7) * t264 - t195 + t262 + t279 - t95 + (-mrSges(3,3) - t357) * qJD(1)) * qJD(1); (-t179 * t68 - t180 * t69 + t311) * mrSges(6,3) + (t184 * t339 + t233 * t183 / 0.2e1 + ((-Ifges(4,1) * t229 - t325) * t337 + (-Ifges(4,2) * t233 - t326) * t339 - t251) * qJD(1) + t279 * t233 * pkin(3) + (-Ifges(4,5) * t229 / 0.2e1 + Ifges(4,6) * t337) * qJD(3)) * qJD(1) - t26 * t83 - t25 * t84 + (-t139 * t84 + (-qJD(6) * t83 - t15) * t173 + (-t3 - t303) * mrSges(7,3)) * t226 + t236 + (t139 * t83 + t173 * t16 + (-t20 * mrSges(7,3) - t173 * t84) * qJD(6)) * t230 + (-qJD(3) * t274 - t256) * t204 + ((t161 * t232 - t162 * t228) * qJD(4) + (-t141 * t232 + t142 * t228) * mrSges(5,3)) * pkin(3) + t121 * t310 - t160 * t95 - t125 * t161 - t124 * t162 + t172 * t14 + t358 * t117 + t307 * t359 + (t12 * t180 - t13 * t179 - t155 * t160 + t358 * t53 + t359 * t52) * m(6) + ((t228 * t81 + t232 * t82 + (-t120 * t228 + t121 * t232) * qJD(4)) * pkin(3) - t120 * t124 - t121 * t125) * m(5) + (t13 * t172 + t139 * t356 + t173 * t244 - t20 * t25 - t21 * t26 - t359 * t50) * m(7); t307 * t56 + (-t186 * t95 + (-t227 * t69 - t231 * t68) * mrSges(6,3) + (-t307 * t227 + t253 * t231 + m(7) * (t227 * t50 + t231 * t356)) * qJD(5) + (0.2e1 * t155 * t342 + t12 * t227 - t13 * t231 + (-t227 * t52 + t231 * t53) * qJD(5)) * m(6)) * pkin(4) + t238 * (pkin(4) * t227 + pkin(10)) - t24 * t83 - t23 * t84 + t245 + t236 + (t162 + t310) * t121 - m(6) * (-t52 * t56 + t53 * t57) - m(7) * (t20 * t23 + t21 * t24 + t50 * t56) + mrSges(6,3) * t311 - t57 * t117 - t120 * t161 + t364 * (-pkin(4) * t231 - pkin(5)); -t28 * t83 - t27 * t84 + t245 + t307 * t53 + (t330 + (Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t258 - t354) * t276 + t238 * pkin(10) + (t53 * mrSges(6,3) - t239) * t258 - m(7) * (t20 * t27 + t21 * t28 + t50 * t53) + t240 - t52 * t117 - t364 * pkin(5); t65 - t50 * (mrSges(7,1) * t115 + mrSges(7,2) * t114) + (Ifges(7,1) * t114 - t313) * t347 + t48 * t346 + (Ifges(7,5) * t114 - Ifges(7,6) * t115) * t345 - t20 * t83 + t21 * t84 + (t114 * t20 + t115 * t21) * mrSges(7,3) + (-Ifges(7,2) * t115 + t113 + t49) * t348 + t355;];
tauc  = t1(:);
