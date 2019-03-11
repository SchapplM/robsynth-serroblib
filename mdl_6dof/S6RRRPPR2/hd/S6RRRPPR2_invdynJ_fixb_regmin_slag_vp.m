% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRPPR2
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% tau_reg [6x30]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRPPR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:26:46
% EndTime: 2019-03-09 15:26:55
% DurationCPUTime: 3.92s
% Computational Cost: add. (7290->446), mult. (17433->551), div. (0->0), fcn. (13038->14), ass. (0->251)
t200 = qJ(2) + qJ(3);
t190 = sin(t200);
t191 = cos(t200);
t206 = sin(qJ(1));
t210 = cos(qJ(1));
t248 = g(1) * t210 + g(2) * t206;
t344 = -g(3) * t191 + t190 * t248;
t208 = cos(qJ(3));
t209 = cos(qJ(2));
t283 = qJD(1) * t209;
t270 = t208 * t283;
t204 = sin(qJ(3));
t205 = sin(qJ(2));
t284 = qJD(1) * t205;
t271 = t204 * t284;
t128 = -t270 + t271;
t130 = -t204 * t283 - t208 * t284;
t201 = sin(pkin(10));
t202 = cos(pkin(10));
t238 = -t128 * t201 - t202 * t130;
t334 = qJD(6) + t238;
t203 = sin(qJ(6));
t339 = t203 * t334;
t207 = cos(qJ(6));
t196 = qJD(2) + qJD(3);
t276 = t209 * qJDD(1);
t278 = qJD(1) * qJD(2);
t267 = t209 * t278;
t277 = t205 * qJDD(1);
t338 = t267 + t277;
t81 = qJD(3) * t270 - t196 * t271 + t204 * t276 + t338 * t208;
t141 = t204 * t209 + t205 * t208;
t106 = t196 * t141;
t239 = t204 * t277 - t208 * t276;
t82 = qJD(1) * t106 + t239;
t56 = -t201 * t82 + t202 * t81;
t52 = qJDD(6) + t56;
t49 = t207 * t52;
t343 = -t334 * t339 + t49;
t97 = t128 * t202 - t130 * t201;
t328 = pkin(5) * t97;
t342 = t196 * t97;
t83 = t196 * t203 - t207 * t97;
t341 = t334 * t83;
t172 = t201 * t204 * pkin(2);
t281 = qJD(3) * t208;
t330 = pkin(7) + pkin(8);
t155 = t330 * t209;
t148 = qJD(1) * t155;
t135 = t208 * t148;
t154 = t330 * t205;
t146 = qJD(1) * t154;
t255 = t146 * t204 - t135;
t299 = qJ(4) * t128;
t86 = t255 + t299;
t124 = t130 * qJ(4);
t131 = t204 * t148;
t288 = -t208 * t146 - t131;
t87 = t124 + t288;
t300 = -qJD(3) * t172 - t201 * t86 + (pkin(2) * t281 - t87) * t202;
t337 = t238 ^ 2;
t329 = pkin(5) * t238;
t192 = t209 * pkin(2);
t314 = pkin(1) + t192;
t303 = qJD(5) + t300;
t296 = t202 * t204;
t301 = -t201 * t87 + t202 * t86 + (t201 * t208 + t296) * qJD(3) * pkin(2);
t313 = qJD(2) * pkin(2);
t138 = -t146 + t313;
t237 = -t138 * t204 - t135;
t80 = -t237 - t299;
t74 = t201 * t80;
t256 = t208 * t138 - t131;
t79 = t124 + t256;
t54 = t202 * t79 - t74;
t290 = qJD(5) - t54;
t336 = -qJD(6) + t334;
t287 = -t204 * t154 + t208 * t155;
t189 = pkin(10) + t200;
t177 = sin(t189);
t178 = cos(t189);
t240 = t178 * pkin(4) + t177 * qJ(5);
t335 = g(1) * t206 - g(2) * t210;
t165 = g(3) * t178;
t194 = qJDD(2) + qJDD(3);
t107 = qJDD(2) * pkin(2) - t330 * t338;
t268 = t205 * t278;
t109 = t330 * (-t268 + t276);
t224 = qJD(3) * t237 + t208 * t107 - t204 * t109;
t26 = t194 * pkin(3) - t81 * qJ(4) + t130 * qJD(4) + t224;
t282 = qJD(3) * t204;
t332 = (qJD(3) * t138 + t109) * t208 + t204 * t107 - t148 * t282;
t31 = -t82 * qJ(4) - t128 * qJD(4) + t332;
t10 = -t201 * t31 + t202 * t26;
t243 = qJDD(5) - t10;
t153 = t314 * qJD(1);
t108 = t128 * pkin(3) + qJD(4) - t153;
t222 = -qJ(5) * t238 + t108;
t58 = t97 * pkin(4) + t222;
t217 = -t177 * t248 + t238 * t58 + t165 + t243;
t183 = pkin(2) * t208 + pkin(3);
t122 = t183 * t202 - t172;
t118 = -pkin(4) - t122;
t111 = -pkin(9) + t118;
t333 = (t328 + t301) * t334 + t111 * t52;
t331 = pkin(4) + pkin(9);
t327 = pkin(3) * t130;
t326 = pkin(3) * t190;
t325 = pkin(4) * t177;
t320 = t194 * pkin(4);
t140 = t204 * t205 - t208 * t209;
t102 = t202 * t140 + t141 * t201;
t103 = -t140 * t201 + t141 * t202;
t251 = t140 * pkin(3) - t314;
t235 = -t103 * qJ(5) + t251;
t42 = t331 * t102 + t235;
t319 = t42 * t52;
t310 = t202 * t80;
t53 = t201 * t79 + t310;
t318 = t53 * t238;
t85 = t196 * t207 + t203 * t97;
t317 = t85 * t97;
t316 = t334 * t97;
t315 = t97 * t83;
t11 = t201 * t26 + t202 * t31;
t73 = pkin(3) * t196 + t79;
t48 = t201 * t73 + t310;
t46 = -qJ(5) * t196 - t48;
t34 = -t46 - t328;
t312 = t102 * t34;
t309 = t203 * t52;
t308 = t203 * t85;
t306 = t207 * t334;
t279 = qJD(6) * t207;
t280 = qJD(6) * t203;
t55 = t201 * t81 + t202 * t82;
t29 = t207 * t194 - t196 * t280 + t203 * t55 + t97 * t279;
t305 = t29 * t207;
t304 = t329 + t303;
t298 = qJ(5) * t178;
t297 = t130 * t128;
t295 = t203 * t206;
t294 = t203 * t210;
t293 = t206 * t207;
t292 = t207 * t210;
t291 = t329 + t290;
t181 = pkin(3) * t191;
t286 = t181 + t192;
t198 = t205 ^ 2;
t285 = -t209 ^ 2 + t198;
t187 = t205 * t313;
t275 = t181 + t240;
t179 = -pkin(3) * t202 - pkin(4);
t47 = t202 * t73 - t74;
t250 = qJD(5) - t47;
t32 = -t331 * t196 + t250 + t329;
t37 = t331 * t97 + t222;
t16 = t203 * t32 + t207 * t37;
t253 = t194 * qJ(5) + t196 * qJD(5) + t11;
t6 = -pkin(5) * t55 + t253;
t273 = -t16 * t97 + t6 * t207;
t272 = qJD(2) * t330;
t269 = t301 * t238;
t266 = t106 * pkin(3) + t187;
t150 = -pkin(2) * t205 - t326;
t265 = t150 - t325;
t125 = pkin(2) * t268 - qJDD(1) * t314;
t225 = t82 * pkin(3) + qJDD(4) + t125;
t216 = -t56 * qJ(5) - qJD(5) * t238 + t225;
t7 = t331 * t55 + t216;
t263 = qJD(6) * t32 + t7;
t5 = t56 * pkin(5) - t331 * t194 + t243;
t262 = -qJD(6) * t37 + t5;
t147 = t205 * t272;
t149 = t209 * t272;
t232 = -t208 * t147 - t204 * t149 - t154 * t281 - t155 * t282;
t59 = -qJ(4) * t106 - qJD(4) * t140 + t232;
t105 = t196 * t140;
t223 = -t287 * qJD(3) + t204 * t147 - t208 * t149;
t60 = t105 * qJ(4) - t141 * qJD(4) + t223;
t21 = t201 * t59 - t202 * t60;
t254 = -t208 * t154 - t155 * t204;
t94 = -qJ(4) * t141 + t254;
t95 = -qJ(4) * t140 + t287;
t66 = t201 * t95 - t202 * t94;
t261 = t334 * t34;
t260 = t334 ^ 2;
t259 = t203 * t194 - t207 * t55;
t249 = -t325 - t326;
t45 = -pkin(4) * t196 + t250;
t246 = -t238 * t46 + t45 * t97;
t245 = t238 * t48 - t47 * t97;
t244 = -t97 ^ 2 - t337;
t174 = -pkin(9) + t179;
t242 = t174 * t52 - (t53 - t328) * t334;
t15 = -t203 * t37 + t207 * t32;
t241 = t15 * t97 + t6 * t203 + (t207 * t238 + t279) * t34;
t22 = t201 * t60 + t202 * t59;
t67 = t201 * t94 + t202 * t95;
t123 = pkin(2) * t296 + t183 * t201;
t236 = -0.2e1 * pkin(1) * t278 - pkin(7) * qJDD(2);
t64 = pkin(4) * t238 + qJ(5) * t97 - t327;
t43 = pkin(5) * t103 + t66;
t70 = -t105 * t201 + t202 * t106;
t231 = t6 * t102 + t34 * t70 - t43 * t52;
t230 = -t306 * t334 - t309;
t186 = pkin(2) * t284;
t63 = t186 + t64;
t71 = -t105 * t202 - t106 * t201;
t229 = -t71 * qJ(5) - t103 * qJD(5) + t266;
t228 = -g(3) * t177 - t178 * t248;
t211 = qJD(2) ^ 2;
t227 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t211 + t335;
t212 = qJD(1) ^ 2;
t226 = pkin(1) * t212 - pkin(7) * qJDD(1) + t248;
t220 = t21 * t238 - t22 * t97 - t55 * t67 + t56 * t66 - t248;
t92 = t238 * pkin(9);
t219 = (-qJD(6) * t111 + t63 + t92) * t334 + t228;
t218 = (-qJD(6) * t174 + t64 + t92) * t334 + t228;
t215 = g(3) * t190 - t153 * t128 + t248 * t191 - t332;
t214 = -t58 * t97 + t228 + t253;
t12 = t55 * pkin(4) + t216;
t213 = -t153 * t130 + t224 + t344;
t195 = -qJ(4) - t330;
t176 = pkin(3) * t201 + qJ(5);
t152 = t210 * t298;
t151 = t206 * t298;
t145 = pkin(1) + t286;
t137 = t210 * t145;
t117 = qJ(5) + t123;
t116 = -t177 * t295 + t292;
t115 = t177 * t293 + t294;
t114 = t177 * t294 + t293;
t113 = t177 * t292 - t295;
t88 = -t128 ^ 2 + t130 ^ 2;
t69 = -t239 + (-qJD(1) * t141 - t130) * t196;
t68 = t128 * t196 + t81;
t65 = pkin(4) * t102 + t235;
t44 = -pkin(5) * t102 + t67;
t30 = qJD(6) * t85 + t259;
t23 = pkin(4) * t70 + t229;
t20 = t331 * t70 + t229;
t19 = -pkin(5) * t70 + t22;
t18 = pkin(5) * t71 + t21;
t17 = -t339 * t85 + t305;
t14 = t230 - t315;
t13 = t317 + t343;
t9 = t243 - t320;
t2 = t207 * t5;
t1 = (-t334 * t85 - t30) * t207 + (-t29 + t341) * t203;
t3 = [qJDD(1), t335, t248, qJDD(1) * t198 + 0.2e1 * t205 * t267, 0.2e1 * t205 * t276 - 0.2e1 * t285 * t278, qJDD(2) * t205 + t209 * t211, qJDD(2) * t209 - t205 * t211, 0, t205 * t236 + t209 * t227, -t205 * t227 + t209 * t236, t105 * t130 + t141 * t81, t105 * t128 + t106 * t130 - t140 * t81 - t141 * t82, -t105 * t196 + t141 * t194, -t106 * t196 - t140 * t194, 0, -t153 * t106 + t125 * t140 + t128 * t187 + t191 * t335 + t194 * t254 + t196 * t223 - t314 * t82, t153 * t105 + t125 * t141 - t130 * t187 - t190 * t335 - t194 * t287 - t196 * t232 - t314 * t81, -t10 * t103 - t102 * t11 - t47 * t71 - t48 * t70 + t220, t11 * t67 + t48 * t22 - t10 * t66 - t47 * t21 + t225 * t251 + t108 * t266 - g(1) * (-t145 * t206 - t210 * t195) - g(2) * (-t206 * t195 + t137) -t102 * t253 + t103 * t9 + t45 * t71 + t46 * t70 + t220, -t12 * t102 - t178 * t335 + t66 * t194 + t21 * t196 - t23 * t97 - t65 * t55 - t58 * t70, -t12 * t103 + t177 * t335 + t67 * t194 + t22 * t196 - t23 * t238 - t65 * t56 - t58 * t71, -g(2) * t137 + t12 * t65 + t45 * t21 - t46 * t22 + t58 * t23 + t9 * t66 + t253 * t67 + (g(1) * t195 - g(2) * t240) * t210 + (-g(1) * (-t145 - t240) + g(2) * t195) * t206, t70 * t308 + (t29 * t203 + t279 * t85) * t102 (-t203 * t83 + t207 * t85) * t70 + (-t203 * t30 + t305 + (-t207 * t83 - t308) * qJD(6)) * t102, t70 * t339 + t29 * t103 + t85 * t71 + (t279 * t334 + t309) * t102, t70 * t306 - t30 * t103 - t83 * t71 + (-t280 * t334 + t49) * t102, t103 * t52 + t334 * t71, -g(1) * t116 - g(2) * t114 + t2 * t103 + t15 * t71 + t19 * t83 + t44 * t30 + (-t7 * t103 - t20 * t334 - t319) * t203 + (t18 * t334 - t231) * t207 + ((-t203 * t43 - t207 * t42) * t334 - t16 * t103 + t203 * t312) * qJD(6), g(1) * t115 - g(2) * t113 - t16 * t71 + t19 * t85 + t44 * t29 + (-(qJD(6) * t43 + t20) * t334 - t319 - t263 * t103 + qJD(6) * t312) * t207 + (-(-qJD(6) * t42 + t18) * t334 - t262 * t103 + t231) * t203; 0, 0, 0, -t205 * t212 * t209, t285 * t212, t277, t276, qJDD(2), -g(3) * t209 + t205 * t226, g(3) * t205 + t209 * t226, -t297, t88, t68, t69, t194, -t255 * t196 + (-t128 * t284 + t194 * t208 - t196 * t282) * pkin(2) + t213, t288 * t196 + (t130 * t284 - t194 * t204 - t196 * t281) * pkin(2) + t215, -t122 * t56 - t123 * t55 - t300 * t97 + t245 + t269, t11 * t123 + t10 * t122 - t108 * (t186 - t327) - g(3) * t286 + t300 * t48 - t301 * t47 - t248 * t150, -t117 * t55 + t118 * t56 - t303 * t97 + t246 + t269, t63 * t97 + t301 * t196 + (-pkin(4) + t118) * t194 + t217, t117 * t194 + t196 * t303 + t238 * t63 + t214, t253 * t117 + t9 * t118 - t58 * t63 - g(1) * (t210 * t265 + t152) - g(2) * (t206 * t265 + t151) - g(3) * (t192 + t275) - t303 * t46 + t301 * t45, t17, t1, t13, t14, t316, t117 * t30 + t219 * t203 + t333 * t207 + t304 * t83 + t241, t117 * t29 + t304 * t85 + t219 * t207 + (-t261 - t333) * t203 + t273; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t297, t88, t68, t69, t194, -t196 * t237 + t213, t196 * t256 + t215, -t318 + t54 * t97 + (-t201 * t55 - t202 * t56) * pkin(3) + t245, t47 * t53 - t48 * t54 + (t10 * t202 + t108 * t130 + t11 * t201 + t344) * pkin(3), -t176 * t55 + t179 * t56 - t290 * t97 + t246 - t318, -t53 * t196 + t64 * t97 + (-pkin(4) + t179) * t194 + t217, t176 * t194 + t196 * t290 + t238 * t64 + t214, t253 * t176 + t9 * t179 - t58 * t64 - t45 * t53 - g(1) * (t210 * t249 + t152) - g(2) * (t206 * t249 + t151) - g(3) * t275 - t290 * t46, t17, t1, t13, t14, t316, t176 * t30 + t203 * t218 + t207 * t242 + t291 * t83 + t241, t176 * t29 + t291 * t85 + (-t261 - t242) * t203 + t218 * t207 + t273; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t244, t238 * t47 + t48 * t97 + t225 - t335, t244, -t196 * t238 - t55, -t56 + t342, -t238 * t45 - t46 * t97 + t12 - t335, 0, 0, 0, 0, 0, t230 + t315, t317 - t343; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56 + t342, -t238 * t97 + t194, -t196 ^ 2 - t337, t46 * t196 + t217 - t320, 0, 0, 0, 0, 0, -t196 * t83 - t203 * t260 + t49, -t196 * t85 - t207 * t260 - t309; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85 * t83, -t83 ^ 2 + t85 ^ 2, t29 + t341, t336 * t85 - t259, t52, -g(1) * t113 - g(2) * t115 + t336 * t16 + t165 * t207 - t203 * t7 - t34 * t85 + t2, g(1) * t114 - g(2) * t116 + t15 * t334 + t34 * t83 - t263 * t207 + (-t262 - t165) * t203;];
tau_reg  = t3;
