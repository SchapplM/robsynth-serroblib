% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRPRP9
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
% Datum: 2019-05-05 18:08
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRPRP9_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP9_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP9_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP9_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:06:55
% EndTime: 2019-05-05 18:07:16
% DurationCPUTime: 7.22s
% Computational Cost: add. (15316->373), mult. (32387->473), div. (0->0), fcn. (21459->8), ass. (0->255)
t312 = pkin(7) + pkin(1);
t229 = sin(qJ(3));
t232 = cos(qJ(3));
t226 = sin(pkin(9));
t227 = cos(pkin(9));
t280 = qJD(1) * t232;
t202 = -qJD(3) * t227 + t226 * t280;
t204 = qJD(3) * t226 + t227 * t280;
t228 = sin(qJ(5));
t231 = cos(qJ(5));
t174 = t202 * t231 + t204 * t228;
t219 = t232 * qJDD(1);
t278 = qJD(1) * qJD(3);
t267 = t229 * t278;
t208 = t219 - t267;
t187 = qJDD(3) * t226 + t208 * t227;
t262 = -qJDD(3) * t227 + t208 * t226;
t122 = -qJD(5) * t174 + t187 * t231 - t228 * t262;
t281 = qJD(1) * t229;
t216 = qJD(5) + t281;
t291 = t174 * t216;
t327 = t122 - t291;
t266 = t232 * t278;
t276 = t229 * qJDD(1);
t207 = t266 + t276;
t205 = qJDD(5) + t207;
t176 = -t202 * t228 + t204 * t231;
t290 = t176 * t174;
t128 = -t290 - t205;
t299 = t128 * t228;
t173 = t176 ^ 2;
t313 = t216 ^ 2;
t324 = -t173 - t313;
t72 = -t231 * t324 - t299;
t298 = t128 * t231;
t79 = -t228 * t324 + t298;
t52 = t226 * t72 + t227 * t79;
t41 = t229 * t52 - t232 * t327;
t50 = t226 * t79 - t227 * t72;
t377 = -qJ(2) * t50 + t312 * t41;
t375 = qJ(4) * t50;
t374 = -pkin(3) * t50 + pkin(4) * t72;
t373 = pkin(3) * t327 - qJ(4) * t52;
t160 = t216 * t176;
t263 = -t187 * t228 - t231 * t262;
t244 = qJD(5) * t176 - t263;
t100 = -t160 + t244;
t314 = t174 ^ 2;
t152 = t314 - t313;
t93 = -t152 * t228 + t298;
t97 = -t152 * t231 - t299;
t372 = -t229 * t100 + t232 * (t226 * t93 - t227 * t97);
t370 = pkin(8) * t72;
t369 = pkin(8) * t79;
t323 = t173 - t314;
t325 = t160 + t244;
t61 = -t228 * t325 + t231 * t327;
t303 = t327 * t228;
t63 = t231 * t325 + t303;
t368 = -t229 * t323 + t232 * (t226 * t61 + t227 * t63);
t367 = t226 * t97 + t227 * t93;
t365 = t226 * t63 - t227 * t61;
t321 = -t290 + t205;
t297 = t321 * t228;
t320 = -t313 - t314;
t329 = t231 * t320 - t297;
t123 = t231 * t321;
t330 = t228 * t320 + t123;
t342 = t226 * t329 + t227 * t330;
t343 = -t226 * t330 + t227 * t329;
t356 = t229 * t343 - t232 * t325;
t363 = qJ(2) * t342 - t312 * t356;
t359 = qJ(4) * t342;
t358 = -pkin(3) * t342 - pkin(4) * t330;
t357 = -pkin(3) * t325 + qJ(4) * t343;
t154 = -t173 + t313;
t344 = t154 * t231 + t297;
t345 = -t154 * t228 + t123;
t355 = t226 * t345 + t227 * t344;
t326 = t122 + t291;
t354 = t232 * (-t226 * t344 + t227 * t345) + t229 * t326;
t112 = -t314 - t173;
t353 = pkin(3) * t112;
t352 = pkin(4) * t112;
t350 = pkin(8) * t329;
t349 = pkin(8) * t330;
t348 = t112 * t232;
t180 = t204 * t202;
t322 = -t180 + t207;
t339 = t226 * t322;
t338 = t227 * t322;
t331 = t327 * qJ(6);
t235 = qJD(1) ^ 2;
t328 = t312 * t235;
t133 = pkin(5) * t174 - qJ(6) * t176;
t271 = t202 * t281;
t165 = -t187 - t271;
t277 = qJD(2) * qJD(1);
t221 = 0.2e1 * t277;
t223 = qJDD(1) * qJ(2);
t230 = sin(qJ(1));
t233 = cos(qJ(1));
t258 = g(1) * t233 + g(2) * t230;
t248 = -t223 + t258;
t241 = t221 - t248;
t250 = -t208 + t267;
t251 = t207 + t266;
t144 = pkin(3) * t251 + qJ(4) * t250 + t241 - t328;
t264 = g(1) * t230 - g(2) * t233;
t255 = qJDD(2) - t264;
t239 = -qJ(2) * t235 + t255;
t188 = -qJDD(1) * t312 + t239;
t179 = t232 * g(3) - t188 * t229;
t234 = qJD(3) ^ 2;
t253 = pkin(3) * t229 - qJ(4) * t232;
t242 = t235 * t253;
t147 = -t234 * pkin(3) + qJDD(3) * qJ(4) - t229 * t242 - t179;
t81 = 0.2e1 * qJD(4) * t204 - t144 * t227 + t147 * t226;
t69 = pkin(4) * t322 + pkin(8) * t165 - t81;
t200 = t202 ^ 2;
t247 = pkin(4) * t281 - pkin(8) * t204;
t82 = -0.2e1 * qJD(4) * t202 + t144 * t226 + t147 * t227;
t71 = -t200 * pkin(4) - pkin(8) * t262 - t247 * t281 + t82;
t38 = t228 * t69 + t231 * t71;
t261 = qJ(6) * t205 - t133 * t174 + t38;
t319 = -pkin(5) * (t324 + t313) - qJ(6) * t128 + t261;
t240 = (-t174 * t228 - t176 * t231) * t216;
t289 = t216 * t228;
t151 = t176 * t289;
t288 = t216 * t231;
t275 = t174 * t288;
t256 = t151 - t275;
t318 = t226 * t256 + t227 * t240;
t245 = t228 * t244 + t275;
t257 = t174 * t289 - t231 * t244;
t317 = t226 * t245 + t227 * t257;
t316 = t232 * (-t226 * t240 + t227 * t256) + t229 * t205;
t274 = t229 * t290;
t315 = t232 * (-t226 * t257 + t227 * t245) - t274;
t201 = t204 ^ 2;
t311 = pkin(5) * t244;
t310 = pkin(5) * t231;
t37 = t228 * t71 - t231 * t69;
t18 = t228 * t38 - t231 * t37;
t309 = t18 * t226;
t308 = t18 * t227;
t305 = qJ(6) * t231;
t304 = qJDD(1) * pkin(1);
t302 = t326 * t228;
t178 = g(3) * t229 + t188 * t232;
t146 = qJDD(3) * pkin(3) + qJ(4) * t234 - t232 * t242 - qJDD(4) + t178;
t109 = -pkin(4) * t262 + pkin(8) * t200 - t204 * t247 + t146;
t301 = t109 * t228;
t300 = t109 * t231;
t296 = t146 * t226;
t295 = t146 * t227;
t167 = t180 + t207;
t293 = t167 * t226;
t292 = t167 * t227;
t225 = t232 ^ 2;
t287 = t225 * t235;
t272 = t229 * t235 * t232;
t284 = t229 * (qJDD(3) + t272);
t283 = t232 * (qJDD(3) - t272);
t279 = qJD(6) * t216;
t273 = t229 * t180;
t270 = t204 * t281;
t269 = t226 * t281;
t268 = t227 * t281;
t265 = -qJ(6) * t228 - pkin(4);
t55 = t226 * t81 + t227 * t82;
t19 = t228 * t37 + t231 * t38;
t88 = t122 * t228 + t176 * t288;
t89 = t122 * t231 - t151;
t260 = t232 * (-t226 * t88 + t227 * t89) + t274;
t213 = 0.2e1 * t279;
t254 = t213 + t261;
t26 = -pkin(5) * t313 + t254;
t27 = -pkin(5) * t205 - qJ(6) * t313 + t133 * t176 + qJDD(6) + t37;
t259 = -pkin(5) * t27 + qJ(6) * t26;
t252 = -pkin(5) * t326 - qJ(6) * t100;
t249 = t226 * t82 - t227 * t81;
t143 = t232 * t178 - t229 * t179;
t246 = qJ(2) + t253;
t163 = -t262 + t270;
t238 = pkin(5) * t321 + qJ(6) * t320 - t27;
t237 = -pkin(5) * t160 + 0.2e1 * qJD(6) * t176 + t109;
t236 = t237 + t331;
t224 = t229 ^ 2;
t220 = t224 * t235;
t210 = (t224 + t225) * qJDD(1);
t209 = t219 - 0.2e1 * t267;
t206 = 0.2e1 * t266 + t276;
t192 = -t239 + t304;
t191 = -t201 - t220;
t190 = -t201 + t220;
t189 = t200 - t220;
t185 = t248 - 0.2e1 * t277 + t328;
t183 = -t284 + t232 * (-t234 - t287);
t182 = t229 * (-t220 - t234) + t283;
t177 = -t220 - t200;
t164 = t187 - t271;
t162 = t262 + t270;
t156 = -t200 - t201;
t140 = -t191 * t226 - t292;
t139 = t191 * t227 - t293;
t130 = t177 * t227 - t339;
t129 = t177 * t226 + t338;
t120 = t163 * t227 - t165 * t226;
t110 = t140 * t229 - t164 * t232;
t108 = t130 * t229 - t162 * t232;
t102 = (-qJD(5) + t216) * t176 + t263;
t98 = t231 * t326;
t83 = t120 * t229 - t156 * t232;
t66 = t102 * t231 + t302;
t64 = -t100 * t231 + t302;
t62 = t102 * t228 - t98;
t60 = -t100 * t228 - t98;
t58 = -t300 + t370;
t57 = t226 * t89 + t227 * t88;
t49 = -t301 - t349;
t44 = t146 * t232 + t229 * t55;
t43 = -pkin(4) * t327 - t301 + t369;
t40 = t236 - t311;
t39 = -pkin(4) * t325 + t300 + t350;
t33 = -t226 * t62 + t227 * t66;
t32 = -t226 * t60 + t227 * t64;
t31 = t226 * t66 + t227 * t62;
t30 = t226 * t64 + t227 * t60;
t29 = (-t325 - t244) * pkin(5) + t236;
t28 = t237 - t311 + 0.2e1 * t331;
t25 = t229 * t33 - t348;
t24 = t229 * t32 - t348;
t23 = -qJ(6) * t112 + t27;
t22 = (-t112 - t313) * pkin(5) + t254;
t21 = -t228 * t29 - t305 * t325 - t349;
t20 = -pkin(5) * t303 + t231 * t28 - t370;
t17 = t231 * t29 + t265 * t325 + t350;
t16 = -t369 + t228 * t28 + (pkin(4) + t310) * t327;
t15 = pkin(4) * t109 + pkin(8) * t19;
t14 = -pkin(8) * t62 - t18;
t13 = t228 * t27 + t231 * t26;
t12 = t228 * t26 - t231 * t27;
t11 = pkin(8) * t66 + t19 - t352;
t10 = -pkin(8) * t60 - t22 * t228 + t23 * t231;
t9 = pkin(8) * t64 + t22 * t231 + t228 * t23 - t352;
t8 = t19 * t227 - t309;
t7 = t19 * t226 + t308;
t6 = t109 * t232 + t229 * t8;
t5 = -pkin(8) * t12 + (-pkin(5) * t228 + t305) * t40;
t4 = -t12 * t226 + t13 * t227;
t3 = t12 * t227 + t13 * t226;
t2 = pkin(8) * t13 + (-t265 + t310) * t40;
t1 = t229 * t4 + t232 * t40;
t34 = [0, 0, 0, 0, 0, qJDD(1), t264, t258, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t255 - 0.2e1 * t304, t221 + 0.2e1 * t223 - t258, pkin(1) * t192 + qJ(2) * (-pkin(1) * t235 + t241), -t250 * t232, -t206 * t232 - t209 * t229, t283 - t229 * (t234 - t287), t251 * t229, t232 * (t220 - t234) - t284, 0, qJ(2) * t206 - t182 * t312 - t185 * t229, qJ(2) * t209 - t183 * t312 - t185 * t232, qJ(2) * (-t220 - t287) + t312 * t210 - t143, -qJ(2) * t185 - t143 * t312, t232 * (t187 * t227 - t204 * t269) + t273, t232 * (-t162 * t227 - t164 * t226) - t229 * (-t201 + t200), t232 * (-t190 * t226 + t338) - t229 * t165, t232 * (t202 * t268 + t226 * t262) - t273, t232 * (t189 * t227 - t293) + t229 * t163, (t207 + (-t202 * t227 + t204 * t226) * t280) * t229, t232 * (-qJ(4) * t129 - t296) - t229 * (-pkin(3) * t129 + t81) + qJ(2) * t129 - t312 * t108, t232 * (-qJ(4) * t139 - t295) - t229 * (-pkin(3) * t139 + t82) + qJ(2) * t139 - t312 * t110, -t232 * t249 - t312 * t83 + t246 * (t163 * t226 + t165 * t227), t246 * t249 - t312 * t44, t260, -t368, t354, t315, t372, t316, t232 * (-t226 * t39 + t227 * t49 - t359) - t229 * (t358 + t37) + t363, t232 * (-t226 * t43 + t227 * t58 - t375) - t229 * (t374 + t38) - t377, t232 * (-qJ(4) * t31 - t11 * t226 + t14 * t227) - t229 * (-pkin(3) * t31 - pkin(4) * t62) + qJ(2) * t31 - t312 * t25, t232 * (-pkin(8) * t308 - qJ(4) * t7 - t15 * t226) - t229 * (-pkin(3) * t7 - pkin(4) * t18) + qJ(2) * t7 - t312 * t6, t260, t354, t368, t316, -t372, t315, t232 * (-t17 * t226 + t21 * t227 - t359) - t229 * (-t238 + t358) + t363, t232 * (-qJ(4) * t30 + t10 * t227 - t226 * t9) - t229 * (-pkin(3) * t30 - pkin(4) * t60 - t252) + qJ(2) * t30 - t312 * t24, t232 * (-t16 * t226 + t20 * t227 + t375) - t229 * (-0.2e1 * t279 - t319 - t374) + t377, t232 * (-qJ(4) * t3 - t2 * t226 + t227 * t5) - t229 * (-pkin(3) * t3 - pkin(4) * t12 - t259) + qJ(2) * t3 - t312 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t235, -t192, 0, 0, 0, 0, 0, 0, t182, t183, -t210, t143, 0, 0, 0, 0, 0, 0, t108, t110, t83, t44, 0, 0, 0, 0, 0, 0, t356, t41, t25, t6, 0, 0, 0, 0, 0, 0, t356, t24, -t41, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t272, -t220 + t287, t219, -t272, -t276, qJDD(3), t178, t179, 0, 0, t187 * t226 + t204 * t268, -t162 * t226 + t164 * t227, t190 * t227 + t339, t202 * t269 - t227 * t262, t189 * t226 + t292, (-t202 * t226 - t204 * t227) * t281, -pkin(3) * t162 + qJ(4) * t130 + t295, -pkin(3) * t164 + qJ(4) * t140 - t296, -pkin(3) * t156 + qJ(4) * t120 + t55, pkin(3) * t146 + qJ(4) * t55, t57, -t365, t355, t317, -t367, t318, t226 * t49 + t227 * t39 + t357, t226 * t58 + t227 * t43 - t373, qJ(4) * t33 + t11 * t227 + t14 * t226 - t353, pkin(3) * t109 - pkin(8) * t309 + qJ(4) * t8 + t15 * t227, t57, t355, t365, t318, t367, t317, t17 * t227 + t21 * t226 + t357, qJ(4) * t32 + t10 * t226 + t227 * t9 - t353, t16 * t227 + t20 * t226 + t373, pkin(3) * t40 + qJ(4) * t4 + t2 * t227 + t226 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t162, t164, t156, -t146, 0, 0, 0, 0, 0, 0, t325, t327, t112, -t109, 0, 0, 0, 0, 0, 0, t325, t112, -t327, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t290, t323, t326, -t290, -t100, t205, -t37, -t38, 0, 0, t290, t326, -t323, t205, t100, -t290, t238, t252, t213 + t319, t259; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t321, t326, t324, t27;];
tauJ_reg  = t34;
