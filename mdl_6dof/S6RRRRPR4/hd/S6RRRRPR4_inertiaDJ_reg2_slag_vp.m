% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRRPR4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR4_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_inertiaDJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:09:20
% EndTime: 2019-03-09 22:09:36
% DurationCPUTime: 6.27s
% Computational Cost: add. (14094->441), mult. (31184->744), div. (0->0), fcn. (30477->10), ass. (0->242)
t218 = sin(qJ(4));
t221 = cos(qJ(4));
t341 = sin(pkin(11));
t342 = cos(pkin(11));
t270 = t341 * t218 - t342 * t221;
t287 = qJD(4) * t341;
t288 = qJD(4) * t342;
t264 = -t218 * t287 + t221 * t288;
t184 = t342 * t218 + t341 * t221;
t222 = cos(qJ(3));
t223 = cos(qJ(2));
t318 = qJD(3) * t222;
t320 = qJD(2) * t223;
t219 = sin(qJ(3));
t220 = sin(qJ(2));
t334 = t219 * t220;
t363 = qJD(2) + qJD(3);
t141 = -t222 * t320 - t223 * t318 + t334 * t363;
t333 = t219 * t223;
t190 = t222 * t220 + t333;
t215 = qJD(4) * t221;
t303 = t190 * t215;
t263 = -t141 * t218 + t303;
t189 = -t222 * t223 + t334;
t212 = -pkin(2) * t223 - pkin(1);
t257 = t189 * pkin(3) - t190 * pkin(9) + t212;
t360 = -pkin(8) - pkin(7);
t199 = t360 * t223;
t319 = qJD(3) * t219;
t183 = t199 * t319;
t306 = qJD(2) * t360;
t192 = t220 * t306;
t198 = t360 * t220;
t281 = t223 * t306;
t362 = -t219 * t281 - (qJD(3) * t198 + t192) * t222 - t183;
t365 = -qJD(4) * t257 + t362;
t216 = t218 ^ 2;
t217 = t221 ^ 2;
t323 = t216 - t217;
t364 = t323 * qJD(4);
t357 = pkin(2) * t219;
t307 = pkin(9) + t357;
t324 = t218 * t288 + t221 * t287;
t305 = t342 * pkin(4);
t279 = t305 + pkin(5);
t304 = t341 * pkin(4);
t358 = sin(qJ(6));
t359 = cos(qJ(6));
t175 = t359 * t279 - t358 * t304;
t316 = qJD(4) * t218;
t332 = t221 * t141;
t262 = t190 * t316 + t332;
t142 = t363 * t190;
t361 = 0.2e1 * t142;
t356 = pkin(9) * t142;
t355 = t184 * pkin(10);
t354 = t222 * pkin(2);
t353 = -qJ(5) - pkin(9);
t256 = t359 * t270;
t135 = t358 * t184 + t256;
t251 = t184 * t141 - t264 * t190;
t152 = t198 * t219 - t222 * t199;
t102 = t152 * qJD(3) + t219 * t192 - t222 * t281;
t60 = t263 * pkin(4) + t102;
t34 = -pkin(5) * t251 + t60;
t151 = -t222 * t198 - t219 * t199;
t337 = t190 * t218;
t118 = pkin(4) * t337 + t151;
t326 = t184 * t190;
t80 = t326 * pkin(5) + t118;
t136 = t359 * t184 - t358 * t270;
t84 = t136 * qJD(6) + t358 * t264 + t359 * t324;
t352 = t34 * t135 + t80 * t84;
t297 = qJD(6) * t358;
t83 = qJD(6) * t256 + t184 * t297 - t359 * t264 + t358 * t324;
t351 = t34 * t136 - t80 * t83;
t336 = t190 * t221;
t91 = -t152 * t218 + t221 * t257;
t68 = pkin(4) * t189 - qJ(5) * t336 + t91;
t92 = t221 * t152 + t218 * t257;
t79 = -qJ(5) * t337 + t92;
t43 = t341 * t68 + t342 * t79;
t213 = pkin(4) * t316;
t154 = t324 * pkin(5) + t213;
t214 = pkin(2) * t319;
t143 = t214 + t154;
t211 = -t221 * pkin(4) - pkin(3);
t157 = t270 * pkin(5) + t211;
t153 = t157 - t354;
t350 = t143 * t135 + t153 * t84;
t349 = t143 * t136 - t153 * t83;
t348 = pkin(2) * qJD(3);
t347 = t118 * t324 + t60 * t270;
t346 = t118 * t264 + t60 * t184;
t345 = t154 * t135 + t157 * t84;
t344 = t154 * t136 - t157 * t83;
t343 = t102 * t218 + t151 * t215;
t339 = t151 * t102;
t338 = t190 * t141;
t335 = t218 * t142;
t331 = t221 * t142;
t193 = t214 + t213;
t197 = t211 - t354;
t330 = t193 * t270 + t197 * t324;
t329 = t193 * t184 + t197 * t264;
t328 = t211 * t324 + t270 * t213;
t327 = t184 * t213 + t211 * t264;
t282 = qJ(5) + t307;
t266 = t282 * t221;
t267 = t218 * t282;
t129 = t342 * t266 - t341 * t267;
t146 = t270 * t353;
t210 = -pkin(3) - t354;
t325 = t210 * t215 + t218 * t214;
t322 = t216 + t217;
t321 = qJD(2) * t220;
t315 = t221 * qJD(5);
t314 = -0.2e1 * pkin(1) * qJD(2);
t114 = t189 * t361;
t310 = pkin(2) * t321;
t252 = t142 * pkin(3) + t141 * pkin(9) + t310;
t313 = t218 * t252 - t365 * t221;
t312 = pkin(3) * t316;
t311 = pkin(3) * t215;
t309 = pkin(2) * t318;
t308 = t218 * t332;
t301 = t218 * t215;
t300 = t220 * t320;
t125 = t270 * t190;
t42 = -t341 * t79 + t342 * t68;
t248 = t189 * pkin(5) + t125 * pkin(10) + t42;
t242 = t359 * t248;
t33 = -t326 * pkin(10) + t43;
t18 = -t358 * t33 + t242;
t241 = t358 * t248;
t19 = t359 * t33 + t241;
t41 = -t152 * t215 + t365 * t218 + t221 * t252;
t227 = qJ(5) * t262 - t190 * t315 + t41;
t226 = t142 * pkin(4) + t227;
t237 = qJ(5) * t303 - (qJ(5) * t141 - qJD(4) * t152 - qJD(5) * t190) * t218 - t313;
t235 = t342 * t237;
t12 = t341 * t226 - t235;
t224 = pkin(10) * t251 + t12;
t234 = t341 * t237;
t11 = t342 * t226 + t234;
t59 = -t270 * t141 + t324 * t190;
t225 = t142 * pkin(5) + t59 * pkin(10) + t11;
t3 = -qJD(6) * t242 - t359 * t224 - t358 * t225 + t33 * t297;
t298 = qJD(6) * t359;
t4 = -qJD(6) * t241 - t358 * t224 + t359 * t225 - t33 * t298;
t299 = t3 * t135 - t4 * t136 + t18 * t83 - t19 * t84;
t296 = t324 * pkin(10);
t295 = -t11 * t184 - t12 * t270 - t42 * t264 - t43 * t324;
t294 = qJD(4) * t353;
t182 = t270 * pkin(10);
t110 = -t182 + t129;
t258 = t264 * pkin(10);
t265 = qJD(4) * t282;
t275 = qJD(5) + t309;
t244 = -t218 * t265 + t275 * t221;
t245 = -t275 * t218 - t221 * t265;
t89 = -t341 * t244 + t342 * t245;
t228 = -t258 + t89;
t90 = t342 * t244 + t341 * t245;
t229 = -t296 + t90;
t128 = -t341 * t266 - t342 * t267;
t236 = t128 - t355;
t233 = t359 * t236;
t26 = -qJD(6) * t233 + t110 * t297 - t358 * t228 - t359 * t229;
t232 = t358 * t236;
t27 = -qJD(6) * t232 - t110 * t298 + t359 * t228 - t358 * t229;
t63 = -t358 * t110 + t233;
t64 = t359 * t110 + t232;
t290 = t26 * t135 - t27 * t136 + t63 * t83 - t64 * t84;
t122 = -t182 + t146;
t254 = t218 * t294 + t315;
t255 = -t218 * qJD(5) + t221 * t294;
t120 = -t341 * t254 + t342 * t255;
t230 = -t258 + t120;
t121 = t342 * t254 + t341 * t255;
t231 = -t296 + t121;
t145 = t184 * t353;
t243 = t145 - t355;
t239 = t359 * t243;
t36 = -qJD(6) * t239 + t122 * t297 - t358 * t230 - t359 * t231;
t238 = t358 * t243;
t37 = -qJD(6) * t238 - t122 * t298 + t359 * t230 - t358 * t231;
t71 = -t358 * t122 + t239;
t72 = t359 * t122 + t238;
t289 = t36 * t135 - t37 * t136 + t71 * t83 - t72 * t84;
t286 = t322 * t222;
t285 = -t128 * t264 - t129 * t324 - t89 * t184 - t90 * t270;
t283 = -t120 * t184 - t121 * t270 - t145 * t264 - t146 * t324;
t185 = t190 ^ 2;
t280 = t185 * t301;
t278 = t218 * t307;
t277 = t221 * t307;
t274 = qJD(4) * t307;
t273 = t359 * t326;
t272 = t92 * t218 + t91 * t221;
t271 = t218 * t91 - t221 * t92;
t268 = t210 * t316 - t221 * t214;
t261 = t189 * t316 - t331;
t253 = t307 * t189 - t210 * t190;
t78 = -t359 * t125 - t358 * t326;
t249 = -t219 * t142 + (-t189 * t222 + t190 * t219) * qJD(3);
t40 = t152 * t316 - t313;
t13 = -t272 * qJD(4) - t41 * t218 - t40 * t221;
t176 = t358 * t279 + t359 * t304;
t240 = pkin(2) * t249 - t210 * t141 - t356;
t202 = -0.2e1 * t301;
t201 = 0.2e1 * t301;
t188 = -0.2e1 * t364;
t178 = t286 * t348;
t162 = t176 * qJD(6);
t161 = t175 * qJD(6);
t139 = t151 * t316;
t138 = 0.2e1 * t184 * t264;
t137 = 0.2e1 * t270 * t324;
t119 = (-t342 * t264 - t341 * t324) * pkin(4);
t112 = t189 * t215 + t335;
t95 = -0.2e1 * t184 * t324 - 0.2e1 * t264 * t270;
t88 = t190 * t364 + t308;
t86 = t184 * t142 + t264 * t189;
t85 = -t270 * t142 - t324 * t189;
t77 = -t358 * t125 + t273;
t76 = t323 * t141 - 0.4e1 * t190 * t301;
t62 = -0.2e1 * t136 * t83;
t61 = 0.2e1 * t135 * t84;
t53 = -t135 * t142 - t189 * t84;
t52 = t136 * t142 - t189 * t83;
t47 = -t125 * t264 - t59 * t184;
t46 = -t251 * t270 + t326 * t324;
t35 = 0.2e1 * t135 * t83 - 0.2e1 * t136 * t84;
t25 = -t135 * t161 + t136 * t162 + t175 * t83 - t176 * t84;
t24 = t78 * qJD(6) - t359 * t251 - t358 * t59;
t23 = qJD(6) * t273 - t125 * t297 - t358 * t251 + t359 * t59;
t20 = t125 * t324 + t184 * t251 - t264 * t326 + t59 * t270;
t17 = t135 * t24 + t77 * t84;
t16 = -t136 * t23 - t78 * t83;
t5 = t135 * t23 - t136 * t24 + t77 * t83 - t78 * t84;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t300, 0.2e1 * (-t220 ^ 2 + t223 ^ 2) * qJD(2), 0, -0.2e1 * t300, 0, 0, t220 * t314, t223 * t314, 0, 0, -0.2e1 * t338, 0.2e1 * t141 * t189 - 0.2e1 * t142 * t190, 0, t114, 0, 0, 0.2e1 * t142 * t212 + 0.2e1 * t189 * t310, -0.2e1 * t141 * t212 + 0.2e1 * t190 * t310, 0.2e1 * t102 * t190 - 0.2e1 * t151 * t141 - 0.2e1 * t152 * t142 + 0.2e1 * t189 * t362, 0.2e1 * t152 * (t222 * t192 + t198 * t318 + t183) + 0.2e1 * t339 + 0.2e1 * (t212 * t220 * pkin(2) + t152 * t360 * t333) * qJD(2), -0.2e1 * t217 * t338 - 0.2e1 * t280, 0.2e1 * t185 * t364 + 0.4e1 * t190 * t308, -0.2e1 * t189 * t262 + 0.2e1 * t190 * t331, -0.2e1 * t216 * t338 + 0.2e1 * t280, -0.2e1 * t189 * t263 - 0.2e1 * t190 * t335, t114, 0.2e1 * t102 * t337 + 0.2e1 * t142 * t91 + 0.2e1 * t151 * t263 + 0.2e1 * t189 * t41, 0.2e1 * t102 * t336 - 0.2e1 * t142 * t92 - 0.2e1 * t151 * t262 + 0.2e1 * t189 * t40, 0.2e1 * t272 * t141 + 0.2e1 * (qJD(4) * t271 + t218 * t40 - t221 * t41) * t190, -0.2e1 * t40 * t92 + 0.2e1 * t41 * t91 + 0.2e1 * t339, 0.2e1 * t125 * t59, -0.2e1 * t125 * t251 + 0.2e1 * t326 * t59, -0.2e1 * t125 * t142 - 0.2e1 * t189 * t59, -0.2e1 * t326 * t251, -0.2e1 * t142 * t326 + 0.2e1 * t189 * t251, t114, 0.2e1 * t11 * t189 - 0.2e1 * t118 * t251 + 0.2e1 * t42 * t142 + 0.2e1 * t326 * t60, -0.2e1 * t118 * t59 - 0.2e1 * t12 * t189 - 0.2e1 * t125 * t60 - 0.2e1 * t142 * t43, 0.2e1 * t11 * t125 - 0.2e1 * t12 * t326 + 0.2e1 * t251 * t43 + 0.2e1 * t42 * t59, 0.2e1 * t11 * t42 + 0.2e1 * t118 * t60 + 0.2e1 * t12 * t43, -0.2e1 * t78 * t23, 0.2e1 * t23 * t77 - 0.2e1 * t24 * t78, 0.2e1 * t142 * t78 - 0.2e1 * t189 * t23, 0.2e1 * t77 * t24, -0.2e1 * t142 * t77 - 0.2e1 * t189 * t24, t114, 0.2e1 * t142 * t18 + 0.2e1 * t189 * t4 + 0.2e1 * t24 * t80 + 0.2e1 * t34 * t77, -0.2e1 * t142 * t19 + 0.2e1 * t189 * t3 - 0.2e1 * t23 * t80 + 0.2e1 * t34 * t78, 0.2e1 * t18 * t23 - 0.2e1 * t19 * t24 + 0.2e1 * t3 * t77 - 0.2e1 * t4 * t78, 0.2e1 * t18 * t4 - 0.2e1 * t19 * t3 + 0.2e1 * t34 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t320, 0, -t321, 0, -pkin(7) * t320, pkin(7) * t321, 0, 0, 0, 0, -t141, 0, -t142, 0, -t102, t362 (t141 * t222 + t249) * pkin(2), -t102 * t354 + t151 * t214 + t152 * t309 - t357 * t362, -t88, t76, t112, t88, -t261, 0, t139 + (-qJD(4) * t253 - t102) * t221 + t240 * t218, t221 * t240 + t253 * t316 + t343, t13, -t40 * t277 - t41 * t278 + t102 * t210 + (-t277 * t91 - t278 * t92) * qJD(4) + (t151 * t219 - t222 * t271) * t348, t47, t20, t86, t46, t85, 0, t128 * t142 + t89 * t189 + t193 * t326 - t197 * t251 + t347, -t125 * t193 - t129 * t142 - t189 * t90 - t197 * t59 + t346, t89 * t125 + t128 * t59 + t129 * t251 - t326 * t90 + t295, t11 * t128 + t118 * t193 + t12 * t129 + t197 * t60 + t42 * t89 + t43 * t90, t16, t5, t52, t17, t53, 0, t142 * t63 + t143 * t77 + t153 * t24 + t189 * t27 + t352, -t142 * t64 + t143 * t78 - t153 * t23 + t189 * t26 + t351, t23 * t63 - t24 * t64 + t26 * t77 - t27 * t78 + t299, t143 * t80 + t153 * t34 + t18 * t27 - t19 * t26 - t3 * t64 + t4 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t214, -0.2e1 * t309, 0, 0, t201, t188, 0, t202, 0, 0, 0.2e1 * t268, 0.2e1 * t325, 0.2e1 * t178, 0.2e1 * (t307 * t354 * t322 + t210 * t357) * qJD(3), t138, t95, 0, t137, 0, 0, 0.2e1 * t330, 0.2e1 * t329, 0.2e1 * t285, 0.2e1 * t128 * t89 + 0.2e1 * t129 * t90 + 0.2e1 * t193 * t197, t62, t35, 0, t61, 0, 0, 0.2e1 * t350, 0.2e1 * t349, 0.2e1 * t290, 0.2e1 * t143 * t153 - 0.2e1 * t26 * t64 + 0.2e1 * t27 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141, 0, -t142, 0, -t102, t362, 0, 0, -t88, t76, t112, t88, -t261, 0, t139 + (pkin(3) * t141 - t356) * t218 + (-t102 + (-pkin(3) * t190 - pkin(9) * t189) * qJD(4)) * t221, pkin(3) * t262 + pkin(9) * t261 + t343, t13, -pkin(3) * t102 + pkin(9) * t13, t47, t20, t86, t46, t85, 0, t120 * t189 + t145 * t142 - t211 * t251 + t213 * t326 + t347, -t121 * t189 - t125 * t213 - t142 * t146 - t211 * t59 + t346, t120 * t125 - t121 * t326 + t145 * t59 + t146 * t251 + t295, t11 * t145 + t118 * t213 + t12 * t146 + t120 * t42 + t121 * t43 + t211 * t60, t16, t5, t52, t17, t53, 0, t142 * t71 + t154 * t77 + t157 * t24 + t189 * t37 + t352, -t142 * t72 + t154 * t78 - t157 * t23 + t189 * t36 + t351, t23 * t71 - t24 * t72 + t36 * t77 - t37 * t78 + t299, t154 * t80 + t157 * t34 + t18 * t37 - t19 * t36 - t3 * t72 + t4 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t214, -t309, 0, 0, t201, t188, 0, t202, 0, 0, t268 - t312, -t311 + t325, t178 (-pkin(3) * t219 + pkin(9) * t286) * t348, t138, t95, 0, t137, 0, 0, t328 + t330, t327 + t329, t283 + t285, t120 * t128 + t121 * t129 + t145 * t89 + t146 * t90 + t193 * t211 + t197 * t213, t62, t35, 0, t61, 0, 0, t345 + t350, t344 + t349, t289 + t290, t143 * t157 + t153 * t154 - t26 * t72 + t27 * t71 - t36 * t64 + t37 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t201, t188, 0, t202, 0, 0, -0.2e1 * t312, -0.2e1 * t311, 0, 0, t138, t95, 0, t137, 0, 0, 0.2e1 * t328, 0.2e1 * t327, 0.2e1 * t283, 0.2e1 * t120 * t145 + 0.2e1 * t121 * t146 + 0.2e1 * t211 * t213, t62, t35, 0, t61, 0, 0, 0.2e1 * t345, 0.2e1 * t344, 0.2e1 * t289, 0.2e1 * t154 * t157 - 0.2e1 * t36 * t72 + 0.2e1 * t37 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t262, 0, -t263, t142, t41, t40, 0, 0, 0, 0, -t59, 0, t251, t142, t342 * t227 + t305 * t361 + t234, -0.2e1 * t142 * t304 - t341 * t227 + t235 (t341 * t251 + t342 * t59) * pkin(4) (t342 * t11 + t341 * t12) * pkin(4), 0, 0, -t23, 0, -t24, t142, t175 * t142 - t162 * t189 + t4, -t176 * t142 - t161 * t189 + t3, -t161 * t77 + t162 * t78 + t175 * t23 - t176 * t24, t161 * t19 - t162 * t18 + t175 * t4 - t176 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t215, 0, -t316, 0, -t218 * t309 - t221 * t274, t218 * t274 - t221 * t309, 0, 0, 0, 0, t264, 0, -t324, 0, t89, -t90, t119 (t341 * t90 + t342 * t89) * pkin(4), 0, 0, -t83, 0, -t84, 0, t27, t26, t25, t161 * t64 - t162 * t63 + t175 * t27 - t176 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t215, 0, -t316, 0, -pkin(9) * t215, pkin(9) * t316, 0, 0, 0, 0, t264, 0, -t324, 0, t120, -t121, t119 (t342 * t120 + t341 * t121) * pkin(4), 0, 0, -t83, 0, -t84, 0, t37, t36, t25, t161 * t72 - t162 * t71 + t175 * t37 - t176 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t162, -0.2e1 * t161, 0, 0.2e1 * t161 * t176 - 0.2e1 * t162 * t175; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t251, -t59, 0, t60, 0, 0, 0, 0, 0, 0, t24, -t23, 0, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t324, t264, 0, t193, 0, 0, 0, 0, 0, 0, t84, -t83, 0, t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t324, t264, 0, t213, 0, 0, 0, 0, 0, 0, t84, -t83, 0, t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, -t24, t142, t4, t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, 0, -t84, 0, t27, t26, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, 0, -t84, 0, t37, t36, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t162, -t161, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t1;
