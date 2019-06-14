% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RRPRRP11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 18:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RRPRRP11_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP11_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP11_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP11_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:56:42
% EndTime: 2019-05-06 18:56:59
% DurationCPUTime: 5.57s
% Computational Cost: add. (23806->414), mult. (49378->508), div. (0->0), fcn. (31541->8), ass. (0->271)
t241 = sin(qJ(4));
t245 = cos(qJ(4));
t246 = cos(qJ(2));
t301 = qJD(1) * t246;
t200 = qJD(2) * t241 + t245 * t301;
t202 = qJD(2) * t245 - t241 * t301;
t240 = sin(qJ(5));
t244 = cos(qJ(5));
t170 = t244 * t200 + t202 * t240;
t172 = -t200 * t240 + t202 * t244;
t135 = t172 * t170;
t297 = qJD(1) * qJD(2);
t228 = t246 * t297;
t242 = sin(qJ(2));
t230 = t242 * qJDD(1);
t207 = t230 + t228;
t196 = qJDD(4) + t207;
t192 = qJDD(5) + t196;
t361 = t135 - t192;
t366 = pkin(5) * t361;
t365 = 2 * qJD(3);
t347 = -pkin(2) - pkin(8);
t236 = t242 ^ 2;
t249 = qJD(1) ^ 2;
t231 = t236 * t249;
t248 = qJD(2) ^ 2;
t218 = -t231 - t248;
t306 = t242 * t249;
t290 = t246 * t306;
t213 = -qJDD(2) + t290;
t304 = t246 * t213;
t364 = pkin(7) * (-t218 * t242 + t304);
t177 = t202 * t200;
t355 = -t177 + t196;
t363 = t241 * t355;
t362 = t245 * t355;
t324 = t361 * t240;
t323 = t361 * t244;
t168 = t170 ^ 2;
t302 = qJD(1) * t242;
t224 = qJD(4) + t302;
t216 = qJD(5) + t224;
t215 = t216 ^ 2;
t130 = -t215 - t168;
t87 = t130 * t240 - t323;
t88 = t130 * t244 + t324;
t56 = t241 * t88 + t245 * t87;
t57 = -t241 * t87 + t245 * t88;
t360 = pkin(3) * t56 - qJ(3) * t57;
t169 = t172 ^ 2;
t145 = -t169 - t215;
t124 = t135 + t192;
t326 = t124 * t240;
t104 = t145 * t244 - t326;
t325 = t124 * t244;
t105 = -t145 * t240 - t325;
t70 = t104 * t245 + t105 * t241;
t71 = -t104 * t241 + t105 * t245;
t359 = pkin(3) * t70 - qJ(3) * t71;
t272 = t207 + t228;
t358 = t272 * qJ(3);
t227 = t242 * t297;
t294 = t246 * qJDD(1);
t208 = -t227 + t294;
t281 = t241 * qJDD(2) + t245 * t208;
t163 = -qJD(4) * t202 - t281;
t267 = t245 * qJDD(2) - t241 * t208;
t164 = -qJD(4) * t200 + t267;
t115 = -t170 * qJD(5) + t163 * t240 + t164 * t244;
t152 = t216 * t170;
t357 = t115 - t152;
t185 = t224 * t200;
t142 = t164 + t185;
t243 = sin(qJ(1));
t247 = cos(qJ(1));
t275 = g(1) * t247 + g(2) * t243;
t328 = qJDD(1) * pkin(7);
t190 = -pkin(1) * t249 - t275 + t328;
t308 = t242 * qJ(3);
t342 = pkin(2) * t246;
t273 = -t308 - t342;
t203 = t273 * qJD(1);
t233 = t242 * g(3);
t256 = -t248 * pkin(2) + (qJD(1) * t203 + t190) * t246 - t233;
t254 = qJD(2) * t365 + t256;
t214 = pkin(3) * t302 - qJD(2) * pkin(8);
t223 = pkin(2) * t227;
t287 = qJD(3) * t302;
t226 = -0.2e1 * t287;
t237 = t246 ^ 2;
t285 = g(1) * t243 - t247 * g(2);
t263 = -qJDD(1) * pkin(1) - t285;
t118 = -t214 * t302 + t223 + t226 + (-pkin(3) * t237 - pkin(7)) * t249 + t347 * t208 - t358 + t263;
t318 = t190 * t242;
t259 = -qJDD(2) * pkin(2) - t248 * qJ(3) + t203 * t302 + qJDD(3) + t318;
t129 = pkin(3) * t207 - qJDD(2) * pkin(8) + (-pkin(3) * t297 - pkin(8) * t306 + g(3)) * t246 + t259;
t82 = t118 * t241 - t245 * t129;
t74 = pkin(4) * t355 - pkin(9) * t142 - t82;
t182 = pkin(4) * t224 - pkin(9) * t202;
t194 = t200 ^ 2;
t83 = t118 * t245 + t129 * t241;
t78 = -pkin(4) * t194 + pkin(9) * t163 - t182 * t224 + t83;
t38 = t240 * t78 - t244 * t74;
t354 = qJ(6) * t152 + 0.2e1 * qJD(6) * t172 + t366 + t38;
t284 = -t244 * t163 + t240 * t164;
t114 = -qJD(5) * t172 - t284;
t148 = pkin(5) * t216 - qJ(6) * t172;
t39 = t240 * t74 + t244 * t78;
t27 = -t168 * pkin(5) + t114 * qJ(6) - 0.2e1 * qJD(6) * t170 - t216 * t148 + t39;
t310 = t237 * t249;
t353 = t304 - (-t248 + t310) * t242;
t295 = qJDD(2) * qJ(3);
t351 = t208 * pkin(3) - pkin(8) * t310 + t295;
t350 = -t163 * pkin(4) - t194 * pkin(9) + t182 * t202;
t205 = -0.2e1 * t227 + t294;
t219 = t248 + t310;
t212 = qJDD(2) + t290;
t317 = t212 * t242;
t349 = pkin(7) * (t219 * t246 + t317) - pkin(1) * t205;
t348 = -t114 * pkin(5) - t168 * qJ(6) + t148 * t172 + qJDD(6);
t195 = t202 ^ 2;
t221 = t224 ^ 2;
t100 = t115 + t152;
t257 = (-qJD(5) + t216) * t172 - t284;
t63 = -t100 * t244 + t240 * t257;
t346 = pkin(9) * t63;
t345 = pkin(9) * t87;
t26 = -qJ(6) * t115 - t354;
t330 = t244 * t26;
t10 = t240 * t27 + t330;
t25 = pkin(5) * t26;
t344 = pkin(4) * t10 + t25;
t341 = pkin(9) * t104;
t340 = g(3) * t246;
t117 = -t168 - t169;
t65 = t100 * t240 + t244 * t257;
t33 = t241 * t65 + t245 * t63;
t35 = -t241 * t63 + t245 * t65;
t339 = pkin(7) * (t117 * t246 + t242 * t33) - pkin(1) * t35;
t96 = (qJD(5) + t216) * t172 + t284;
t338 = pkin(7) * (t242 * t56 + t246 * t96) - pkin(1) * t57;
t337 = pkin(7) * (t242 * t70 + t246 * t357) - pkin(1) * t71;
t17 = t240 * t39 - t244 * t38;
t334 = t17 * t241;
t333 = t17 * t245;
t332 = t240 * t26;
t128 = (t365 + t214) * qJD(2) + t256 + t351;
t84 = t128 + t350;
t331 = t240 * t84;
t329 = t244 * t84;
t255 = (-qJD(4) + t224) * t202 - t281;
t110 = -t142 * t245 + t241 * t255;
t327 = t110 * t242;
t156 = t177 + t196;
t321 = t156 * t241;
t320 = t156 * t245;
t316 = t216 * t240;
t315 = t216 * t244;
t312 = t224 * t241;
t311 = t224 * t245;
t309 = t241 * t128;
t307 = t242 * t205;
t305 = t245 * t128;
t210 = t231 + t310;
t303 = (t236 + t237) * t328 + pkin(1) * t210;
t298 = qJD(4) + t224;
t103 = pkin(4) * t104;
t293 = t103 - t39;
t292 = t242 * t135;
t291 = t242 * t177;
t289 = -pkin(4) * t96 + pkin(9) * t88;
t288 = -pkin(4) * t117 + pkin(9) * t65;
t286 = -pkin(4) * t357 + pkin(9) * t105;
t18 = t240 * t38 + t244 * t39;
t180 = t318 + t340;
t181 = t190 * t246 - t233;
t283 = t180 * t242 + t246 * t181;
t280 = qJ(3) * t96 + t347 * t56;
t279 = qJ(3) * t357 + t347 * t70;
t278 = qJ(3) * t117 + t33 * t347;
t61 = pkin(4) * t63;
t277 = pkin(3) * t33 - qJ(3) * t35 + t61;
t86 = pkin(4) * t87;
t276 = -t38 + t86;
t50 = t241 * t83 - t245 * t82;
t271 = t241 * t82 + t245 * t83;
t269 = (-t231 + t248) * t246 + t317;
t266 = pkin(3) * t96 + t347 * t57;
t265 = pkin(3) * t357 + t347 * t71;
t264 = pkin(3) * t117 + t347 * t35;
t261 = pkin(5) * t145 - t27;
t260 = t246 * t347 - pkin(1) - t308;
t258 = t103 + t261;
t189 = pkin(7) * t249 - t263;
t147 = t259 + t340;
t253 = t26 - t366;
t252 = pkin(2) * t208 + t189 - t223;
t251 = t253 + t86;
t146 = t254 + t295;
t250 = t252 + 0.2e1 * t287;
t49 = t84 + t348;
t211 = t231 - t310;
t206 = t230 + 0.2e1 * t228;
t184 = -t195 + t221;
t183 = t194 - t221;
t179 = t272 * t242;
t178 = (t208 - t227) * t246;
t175 = t195 - t194;
t174 = -t195 - t221;
t173 = t206 * t246 + t307;
t166 = -t221 - t194;
t154 = -t194 - t195;
t150 = -t169 + t215;
t149 = t168 - t215;
t141 = t164 - t185;
t140 = -t298 * t200 + t267;
t137 = t298 * t202 + t281;
t133 = t169 - t168;
t131 = t174 * t245 - t321;
t126 = t166 * t241 + t362;
t120 = (-t170 * t244 + t172 * t240) * t216;
t119 = (-t170 * t240 - t172 * t244) * t216;
t109 = t149 * t244 - t326;
t108 = -t150 * t240 - t323;
t107 = t149 * t240 + t325;
t106 = t150 * t244 - t324;
t95 = pkin(5) * t100;
t92 = t115 * t244 - t172 * t316;
t91 = t115 * t240 + t172 * t315;
t90 = -t114 * t240 + t170 * t315;
t89 = t114 * t244 + t170 * t316;
t81 = -t119 * t241 + t120 * t245;
t80 = t242 * t192 + t246 * (-t119 * t245 - t120 * t241);
t79 = -pkin(5) * t357 - qJ(6) * t124;
t76 = -t107 * t241 + t109 * t245;
t75 = -t106 * t241 + t108 * t245;
t66 = t329 - t341;
t64 = -t240 * t357 - t244 * t96;
t62 = -t240 * t96 + t244 * t357;
t59 = -t241 * t91 + t245 * t92;
t58 = -t241 * t89 + t245 * t90;
t52 = t331 - t345;
t48 = t292 + t246 * (-t241 * t92 - t245 * t91);
t47 = -t292 + t246 * (-t241 * t90 - t245 * t89);
t46 = t242 * t100 + t246 * (-t106 * t245 - t108 * t241);
t45 = t242 * t257 + t246 * (-t107 * t245 - t109 * t241);
t43 = -qJ(6) * t145 + t49;
t42 = t286 + t331;
t41 = t289 - t329;
t36 = -pkin(5) * t96 + qJ(6) * t130 - qJD(2) * t214 - t254 - t348 - t350 - t351;
t34 = -t241 * t62 + t245 * t64;
t29 = t242 * t133 + t246 * (-t241 * t64 - t245 * t62);
t24 = -t240 * t79 + t244 * t43 - t341;
t23 = qJ(6) * t323 - t240 * t36 - t345;
t22 = (t100 + t115) * qJ(6) + t354;
t21 = t240 * t43 + t244 * t79 + t286;
t20 = qJ(6) * t324 + t244 * t36 + t289;
t19 = -pkin(5) * t117 + qJ(6) * t257 + t27;
t16 = pkin(4) * t17;
t15 = -pkin(5) * t49 + qJ(6) * t27;
t14 = -pkin(4) * t84 + pkin(9) * t18;
t13 = -t17 - t346;
t12 = t18 + t288;
t11 = t244 * t27 - t332;
t7 = t18 * t241 + t333;
t6 = -t19 * t240 + t22 * t244 - t346;
t5 = t19 * t244 + t22 * t240 + t288;
t3 = t10 * t245 + t11 * t241;
t2 = -pkin(9) * t10 - qJ(6) * t330 - t15 * t240;
t1 = -pkin(4) * t49 + pkin(9) * t11 - qJ(6) * t332 + t15 * t244;
t4 = [0, 0, 0, 0, 0, qJDD(1), t285, t275, 0, 0, t179, t173, t269, t178, -t353, 0, t246 * t189 - t349, -pkin(1) * t206 - t242 * t189 + t364, t283 + t303, pkin(1) * t189 + pkin(7) * t283, 0, -t269, t353, t179, t173, t178, t242 * (qJ(3) * t210 + t259) + (pkin(2) * t210 + t146 + t233) * t246 + t303, t246 * (-pkin(2) * t205 + t226 - t252) + (-t246 * t272 - t307) * qJ(3) + t349, t242 * t250 - t364 + (pkin(1) + t342) * t206 + (t206 + t272) * t308, pkin(7) * (t146 * t246 + t147 * t242) + (pkin(1) - t273) * (t250 + t358), t291 + t246 * (-t164 * t241 - t202 * t311), t242 * t175 + t246 * (t137 * t241 - t141 * t245), t242 * t142 + t246 * (-t184 * t245 - t363), -t291 + t246 * (-t163 * t245 - t200 * t312), t242 * t255 + t246 * (-t183 * t241 - t320), t242 * t196 + t246 * (t200 * t241 + t202 * t245) * t224, t242 * (pkin(3) * t126 - t82) + t246 * (pkin(3) * t137 + t305) + pkin(7) * (t126 * t242 + t137 * t246) + t260 * (t166 * t245 - t363), t242 * (pkin(3) * t131 - t83) + t246 * (pkin(3) * t140 - t309) + pkin(7) * (t131 * t242 + t140 * t246) + t260 * (-t174 * t241 - t320), pkin(3) * t327 + t246 * (pkin(3) * t154 - t271) + pkin(7) * (t154 * t246 + t327) + t260 * (t142 * t241 + t245 * t255), t260 * t271 + (pkin(3) + pkin(7)) * (t128 * t246 + t242 * t50), t48, t29, t46, t47, t45, t80, t242 * (t276 + t360) + t246 * (-t241 * t52 - t245 * t41 + t266) + t338, t242 * (t293 + t359) + t246 * (-t241 * t66 - t245 * t42 + t265) + t337, t242 * t277 + t246 * (-t12 * t245 - t13 * t241 + t264) + t339, t242 * (pkin(3) * t7 + t16) + t246 * (pkin(3) * t84 + pkin(9) * t334 - t14 * t245) + pkin(7) * (t242 * t7 + t246 * t84) + t260 * (t18 * t245 - t334), t48, t29, t46, t47, t45, t80, t242 * (t251 + t360) + t246 * (-t20 * t245 - t23 * t241 + t266) + t338, t242 * (t258 + t359) + t246 * (-t21 * t245 - t24 * t241 + t265) + t337, t242 * (t277 - t95) + t246 * (-t241 * t6 - t245 * t5 + t264) + t339, t242 * (pkin(3) * t3 + t344) + t246 * (pkin(3) * t49 - t1 * t245 - t2 * t241) + pkin(7) * (t242 * t3 + t246 * t49) + t260 * (-t10 * t241 + t11 * t245); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t290, t211, t230, t290, t294, qJDD(2), -t180, -t181, 0, 0, qJDD(2), -t230, -t294, -t290, t211, t290, (-pkin(2) * t242 + qJ(3) * t246) * qJDD(1), -pkin(2) * t212 + qJ(3) * t219 + t147, -pkin(2) * t218 + (qJDD(2) - t213) * qJ(3) + t254, -pkin(2) * t147 + qJ(3) * t146, t164 * t245 - t202 * t312, -t137 * t245 - t141 * t241, -t184 * t241 + t362, -t163 * t241 + t200 * t311, t183 * t245 - t321, (-t200 * t245 + t202 * t241) * t224, qJ(3) * t137 + t347 * t126 + t309, qJ(3) * t140 + t347 * t131 + t305, qJ(3) * t154 + t347 * t110 - t50, qJ(3) * t128 + t347 * t50, t59, t34, t75, t58, t76, t81, -t241 * t41 + t245 * t52 + t280, -t241 * t42 + t245 * t66 + t279, -t12 * t241 + t13 * t245 + t278, -pkin(9) * t333 + qJ(3) * t84 - t14 * t241 + t347 * t7, t59, t34, t75, t58, t76, t81, -t20 * t241 + t23 * t245 + t280, -t21 * t241 + t24 * t245 + t279, -t241 * t5 + t245 * t6 + t278, qJ(3) * t49 - t1 * t241 + t2 * t245 + t347 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t230, t212, t218, t147, 0, 0, 0, 0, 0, 0, t126, t131, t110, t50, 0, 0, 0, 0, 0, 0, t56, t70, t33, t7, 0, 0, 0, 0, 0, 0, t56, t70, t33, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t177, t175, t142, -t177, t255, t196, -t82, -t83, 0, 0, t135, t133, t100, -t135, t257, t192, t276, t293, t61, t16, t135, t133, t100, -t135, t257, t192, t251, t258, -t95 + t61, t344; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, t133, t100, -t135, t257, t192, -t38, -t39, 0, 0, t135, t133, t100, -t135, t257, t192, t253, t261, -t95, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, t357, t117, t49;];
tauJ_reg  = t4;
