% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRPRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 18:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRPRP8_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP8_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP8_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP8_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:02:16
% EndTime: 2019-05-05 18:02:32
% DurationCPUTime: 6.33s
% Computational Cost: add. (13701->389), mult. (30222->493), div. (0->0), fcn. (20309->8), ass. (0->256)
t213 = sin(pkin(9));
t214 = cos(pkin(9));
t217 = sin(qJ(3));
t220 = cos(qJ(3));
t282 = t213 * t220;
t192 = (t214 * t217 + t282) * qJD(1);
t188 = qJD(5) + t192;
t312 = t188 ^ 2;
t274 = qJD(1) * t220;
t194 = -qJD(1) * t213 * t217 + t214 * t274;
t216 = sin(qJ(5));
t219 = cos(qJ(5));
t172 = -t219 * qJD(3) + t194 * t216;
t313 = t172 ^ 2;
t146 = t313 - t312;
t266 = qJD(1) * qJD(3);
t256 = t220 * t266;
t264 = t217 * qJDD(1);
t198 = -t256 - t264;
t206 = t220 * qJDD(1);
t257 = t217 * t266;
t199 = t206 - t257;
t252 = -t214 * t198 + t199 * t213;
t162 = qJDD(5) + t252;
t174 = qJD(3) * t216 + t194 * t219;
t288 = t174 * t172;
t110 = -t288 - t162;
t298 = t110 * t216;
t82 = -t146 * t219 - t298;
t150 = t188 * t174;
t165 = t198 * t213 + t199 * t214;
t253 = t219 * qJDD(3) - t216 * t165;
t235 = qJD(5) * t174 - t253;
t93 = -t150 + t235;
t375 = t217 * (t213 * t82 - t214 * t93) - t220 * (t213 * t93 + t214 * t82);
t171 = t174 ^ 2;
t325 = -t171 - t312;
t73 = t219 * t325 + t298;
t374 = pkin(3) * t73;
t373 = pkin(4) * t73;
t372 = pkin(8) * t73;
t297 = t110 * t219;
t75 = -t216 * t325 + t297;
t371 = pkin(8) * t75;
t370 = qJ(2) * t73;
t369 = t213 * t75;
t368 = t214 * t75;
t324 = t171 - t313;
t241 = -t216 * qJDD(3) - t219 * t165;
t232 = -qJD(5) * t172 - t241;
t289 = t172 * t188;
t321 = -t289 + t232;
t304 = t216 * t321;
t326 = t150 + t235;
t57 = t219 * t326 + t304;
t365 = t217 * (t213 * t57 + t214 * t324) - t220 * (-t213 * t324 + t214 * t57);
t311 = pkin(7) + pkin(1);
t322 = -t288 + t162;
t105 = t219 * t322;
t319 = -t312 - t313;
t329 = t216 * t319 + t105;
t296 = t322 * t216;
t328 = t219 * t319 - t296;
t341 = t213 * t326 + t214 * t328;
t342 = t213 * t328 - t214 * t326;
t355 = t217 * t341 + t220 * t342;
t362 = qJ(2) * t329 - t311 * t355;
t361 = qJ(4) * t342;
t285 = t194 * t192;
t343 = qJDD(3) - t285;
t360 = t213 * t343;
t359 = t214 * t343;
t358 = -t146 * t216 + t297;
t357 = pkin(3) * t342 + pkin(8) * t328;
t356 = -pkin(3) * t329 + qJ(4) * t341;
t320 = t289 + t232;
t147 = -t171 + t312;
t345 = -t147 * t216 + t105;
t354 = t220 * (t213 * t320 + t214 * t345) - t217 * (t213 * t345 - t214 * t320);
t352 = pkin(4) * t329;
t350 = pkin(8) * t329;
t348 = qJ(6) * t321;
t344 = t219 * t147 + t296;
t273 = qJD(3) * t192;
t141 = t165 - t273;
t323 = t171 + t313;
t340 = pkin(4) * t323;
t338 = t213 * t323;
t334 = t214 * t323;
t272 = qJD(3) * t194;
t139 = t252 + t272;
t327 = -t216 * t326 + t219 * t321;
t218 = sin(qJ(1));
t221 = cos(qJ(1));
t246 = t218 * g(1) - t221 * g(2);
t238 = qJDD(2) - t246;
t223 = qJD(1) ^ 2;
t278 = t223 * qJ(2);
t231 = t238 - t278;
t229 = -t311 * qJDD(1) + t231;
t167 = t220 * g(3) - t217 * t229;
t237 = qJD(3) * pkin(3) - qJ(4) * t274;
t210 = t217 ^ 2;
t284 = t210 * t223;
t130 = -pkin(3) * t284 + t198 * qJ(4) - qJD(3) * t237 - t167;
t226 = t220 * t229;
t279 = t220 * t223;
t224 = t226 - t199 * qJ(4) + qJDD(3) * pkin(3) + (-pkin(3) * t279 - qJ(4) * t266 + g(3)) * t217;
t79 = -0.2e1 * qJD(4) * t192 + t214 * t130 + t213 * t224;
t318 = pkin(5) * t235 - t348;
t132 = pkin(5) * t172 - qJ(6) * t174;
t154 = pkin(4) * t192 - pkin(8) * t194;
t222 = qJD(3) ^ 2;
t64 = -pkin(4) * t222 + qJDD(3) * pkin(8) - t154 * t192 + t79;
t265 = qJD(2) * qJD(1);
t208 = 0.2e1 * t265;
t209 = qJDD(1) * qJ(2);
t247 = t221 * g(1) + t218 * g(2);
t240 = -t209 + t247;
t316 = -t198 * pkin(3) - (qJ(4) * t210 + t311) * t223 + t237 * t274 + qJDD(4) - t240;
t68 = t139 * pkin(4) - t141 * pkin(8) + t208 + t316;
t38 = t216 * t68 + t219 * t64;
t251 = t162 * qJ(6) - t172 * t132 + t38;
t317 = -(t325 + t312) * pkin(5) - qJ(6) * t110 + t251;
t286 = t188 * t219;
t260 = t172 * t286;
t236 = t216 * t235 + t260;
t261 = t214 * t288;
t262 = t213 * t288;
t315 = t220 * (t214 * t236 - t262) - t217 * (t213 * t236 + t261);
t287 = t188 * t216;
t144 = t174 * t287;
t244 = t144 - t260;
t314 = t220 * (t162 * t213 + t214 * t244) - t217 * (-t214 * t162 + t213 * t244);
t190 = t192 ^ 2;
t191 = t194 ^ 2;
t310 = pkin(5) * t219;
t254 = t213 * t130 - t214 * t224;
t239 = -qJDD(3) * pkin(4) - t222 * pkin(8) + t254;
t63 = (0.2e1 * qJD(4) + t154) * t194 + t239;
t307 = t216 * t63;
t305 = t216 * t320;
t303 = t219 * t63;
t268 = qJD(4) * t194;
t78 = t254 + 0.2e1 * t268;
t45 = t213 * t79 - t214 * t78;
t301 = t220 * t45;
t300 = qJ(6) * t219;
t299 = qJDD(1) * pkin(1);
t263 = -0.2e1 * t265;
t131 = t263 - t316;
t295 = t131 * t213;
t294 = t131 * t214;
t158 = qJDD(3) + t285;
t292 = t158 * t213;
t291 = t158 * t214;
t211 = t220 ^ 2;
t283 = t211 * t223;
t259 = t217 * t279;
t281 = t217 * (qJDD(3) + t259);
t280 = t220 * (qJDD(3) - t259);
t275 = t210 + t211;
t271 = qJD(3) * t213;
t270 = qJD(3) * t214;
t267 = qJD(6) * t188;
t258 = -pkin(4) * t214 - pkin(3);
t255 = -qJ(6) * t216 - pkin(4);
t46 = t213 * t78 + t214 * t79;
t37 = t216 * t64 - t219 * t68;
t18 = t216 * t37 + t219 * t38;
t179 = 0.2e1 * t267;
t243 = t179 + t251;
t32 = -pkin(5) * t312 + t243;
t33 = -t162 * pkin(5) - qJ(6) * t312 + t132 * t174 + qJDD(6) + t37;
t249 = -pkin(5) * t33 + qJ(6) * t32;
t248 = -pkin(5) * t320 - qJ(6) * t93;
t245 = t172 * t287 - t219 * t235;
t11 = t18 * t213 - t214 * t63;
t4 = t11 * t220 + (t18 * t214 + t213 * t63) * t217;
t17 = t216 * t38 - t219 * t37;
t166 = t217 * g(3) + t226;
t124 = t220 * t166 - t217 * t167;
t140 = -t252 + t272;
t233 = (-t172 * t216 - t174 * t219) * t188;
t90 = t219 * t232 - t144;
t228 = t220 * (t214 * t90 + t262) - t217 * (t213 * t90 - t261);
t185 = -0.2e1 * t268;
t227 = 0.2e1 * qJD(6) * t174 - t194 * t154 + t185 - t239 - t318;
t225 = pkin(5) * t322 + qJ(6) * t319 - t33;
t201 = t275 * qJDD(1);
t200 = t206 - 0.2e1 * t257;
t197 = 0.2e1 * t256 + t264;
t189 = -t231 + t299;
t182 = -t191 - t222;
t181 = -t191 + t222;
t180 = t190 - t222;
t178 = t311 * t223 + t240 + t263;
t176 = -t281 + t220 * (-t222 - t283);
t175 = t217 * (-t222 - t284) + t280;
t155 = -t222 - t190;
t142 = t165 + t273;
t135 = -t190 - t191;
t127 = -t182 * t213 - t291;
t126 = t182 * t214 - t292;
t113 = t155 * t214 - t360;
t112 = t155 * t213 + t359;
t102 = t140 * t214 + t142 * t213;
t101 = t140 * t213 - t142 * t214;
t100 = (qJD(5) + t188) * t172 + t241;
t95 = (-qJD(5) + t188) * t174 + t253;
t91 = t219 * t320;
t89 = t174 * t286 + t216 * t232;
t84 = t126 * t220 + t127 * t217;
t65 = t112 * t220 + t113 * t217;
t61 = t101 * t220 + t102 * t217;
t60 = t219 * t95 + t305;
t58 = -t219 * t93 + t305;
t56 = t216 * t95 - t91;
t55 = -t216 * t93 - t91;
t53 = -t100 * t213 + t368;
t51 = t100 * t214 + t369;
t49 = -t213 * t321 - t368;
t47 = t214 * t321 - t369;
t44 = t214 * t60 - t338;
t43 = t214 * t58 - t338;
t42 = t213 * t60 + t334;
t41 = t213 * t58 + t334;
t40 = t303 - t372;
t39 = t307 - t350;
t35 = -pkin(4) * t55 - t248;
t34 = (pkin(5) * t188 - 0.2e1 * qJD(6)) * t174 + t63 + t318;
t31 = t38 - t373;
t30 = t37 - t352;
t28 = t217 * t53 + t220 * t51;
t26 = t217 * t49 + t220 * t47;
t25 = (-t326 - t150) * pkin(5) + t227;
t24 = -pkin(5) * t150 + t227 + t348;
t23 = qJ(6) * t323 + t33;
t22 = (t323 - t312) * pkin(5) + t243;
t21 = t217 * t46 + t301;
t20 = t217 * t44 + t220 * t42;
t19 = t217 * t43 + t220 * t41;
t16 = -t225 - t352;
t15 = -t216 * t25 - t300 * t326 - t350;
t14 = -0.2e1 * t267 - t317 + t373;
t13 = -pkin(5) * t304 + t219 * t24 + t372;
t10 = -pkin(8) * t56 - t17;
t9 = t216 * t33 + t219 * t32;
t8 = t216 * t32 - t219 * t33;
t7 = -pkin(8) * t55 - t216 * t22 + t219 * t23;
t6 = t213 * t34 + t214 * t9;
t5 = t213 * t9 - t214 * t34;
t3 = -pkin(8) * t8 + (pkin(5) * t216 - t300) * t34;
t2 = -pkin(4) * t8 - t249;
t1 = t217 * t6 + t220 * t5;
t12 = [0, 0, 0, 0, 0, qJDD(1), t246, t247, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t238 - 0.2e1 * t299, t208 + 0.2e1 * t209 - t247, pkin(1) * t189 + qJ(2) * (-t223 * pkin(1) + t208 - t240), (t199 - t257) * t220, -t197 * t220 - t200 * t217, t280 - t217 * (t222 - t283), (-t198 + t256) * t217, t220 * (-t222 + t284) - t281, 0, qJ(2) * t197 - t311 * t175 - t217 * t178, qJ(2) * t200 - t311 * t176 - t220 * t178, t311 * t201 - t275 * t278 - t124, -qJ(2) * t178 - t311 * t124, t220 * (t165 * t214 - t194 * t271) - t217 * (t165 * t213 + t194 * t270), t220 * (-t139 * t214 - t141 * t213) - t217 * (-t139 * t213 + t141 * t214), t220 * (-t181 * t213 + t359) - t217 * (t181 * t214 + t360), t220 * (t192 * t270 + t213 * t252) - t217 * (t192 * t271 - t214 * t252), t220 * (t180 * t214 - t292) - t217 * (t180 * t213 + t291), (t220 * (-t192 * t214 + t194 * t213) - t217 * (-t192 * t213 - t194 * t214)) * qJD(3), t220 * (-qJ(4) * t112 - t295) - t217 * (-pkin(3) * t139 + qJ(4) * t113 + t294) + qJ(2) * t139 - t311 * t65, t220 * (-qJ(4) * t126 - t294) - t217 * (-pkin(3) * t141 + qJ(4) * t127 - t295) + qJ(2) * t141 - t311 * t84, t220 * (-qJ(4) * t101 - t45) - t217 * (-pkin(3) * t135 + qJ(4) * t102 + t46) + qJ(2) * t135 - t311 * t61, -qJ(4) * t301 - t217 * (pkin(3) * t131 + qJ(4) * t46) - qJ(2) * t131 - t311 * t21, t228, t365, t354, t315, t375, t314, t220 * (-t213 * t30 + t214 * t39 - t361) - t217 * (t213 * t39 + t214 * t30 + t356) + t362, t220 * (-qJ(4) * t51 - t213 * t31 + t214 * t40) - t217 * (qJ(4) * t53 + t213 * t40 + t214 * t31 - t374) + t370 - t311 * t28, t220 * (-qJ(4) * t42 + t10 * t214) - t217 * (qJ(4) * t44 + t213 * t10) + (pkin(4) * t282 - t217 * t258 + qJ(2)) * t56 - t311 * t20, (t220 * (pkin(4) * t213 - pkin(8) * t214) - t217 * (-pkin(8) * t213 + t258) + qJ(2)) * t17 + (-t311 - qJ(4)) * t4, t228, t354, -t365, t314, -t375, t315, t220 * (t15 * t214 - t16 * t213 - t361) - t217 * (t15 * t213 + t16 * t214 + t356) + t362, t220 * (-qJ(4) * t41 - t213 * t35 + t214 * t7) - t217 * (-pkin(3) * t55 + qJ(4) * t43 + t213 * t7 + t214 * t35) + qJ(2) * t55 - t311 * t19, t220 * (-qJ(4) * t47 + t13 * t214 - t14 * t213) - t217 * (qJ(4) * t49 + t13 * t213 + t14 * t214 + t374) - t370 - t311 * t26, t220 * (-qJ(4) * t5 - t2 * t213 + t214 * t3) - t217 * (-pkin(3) * t8 + qJ(4) * t6 + t2 * t214 + t213 * t3) + qJ(2) * t8 - t311 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t223, -t189, 0, 0, 0, 0, 0, 0, t175, t176, -t201, t124, 0, 0, 0, 0, 0, 0, t65, t84, t61, t21, 0, 0, 0, 0, 0, 0, t355, t28, t20, t4, 0, 0, 0, 0, 0, 0, t355, t19, t26, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t259, (-t210 + t211) * t223, t206, -t259, -t264, qJDD(3), t166, t167, 0, 0, t285, t191 - t190, t142, -t285, t140, qJDD(3), pkin(3) * t112 + t185 - t254, pkin(3) * t126 - t79, pkin(3) * t101, pkin(3) * t45, t89, t327, t344, t245, -t358, t233, -pkin(4) * t326 - t303 + t357, pkin(3) * t51 + pkin(4) * t100 + t307 + t371, pkin(3) * t42 + pkin(8) * t60 + t18 + t340, pkin(3) * t11 - pkin(4) * t63 + pkin(8) * t18, t89, t344, -t327, t233, t358, t245, t219 * t25 + t255 * t326 + t357, pkin(3) * t41 + pkin(8) * t58 + t216 * t23 + t219 * t22 + t340, pkin(3) * t47 - t371 + t216 * t24 + (pkin(4) + t310) * t321, pkin(3) * t5 + pkin(8) * t9 + (t255 - t310) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, t141, t135, -t131, 0, 0, 0, 0, 0, 0, t329, t73, t56, t17, 0, 0, 0, 0, 0, 0, t329, t55, -t73, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t288, t324, t320, -t288, -t93, t162, -t37, -t38, 0, 0, t288, t320, -t324, t162, t93, -t288, t225, t248, t179 + t317, t249; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t322, t320, t325, t33;];
tauJ_reg  = t12;
