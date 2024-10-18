% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
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
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRR11_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR11_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR11_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR11_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:57
% EndTime: 2024-09-27 21:46:03
% DurationCPUTime: 3.23s
% Computational Cost: add. (6116->439), mult. (17361->603), div. (0->0), fcn. (13541->22), ass. (0->252)
t215 = sin(pkin(5));
t218 = sin(qJ(4));
t221 = cos(qJ(4));
t222 = cos(qJ(3));
t292 = qJDD(2) * t222;
t274 = t215 * t292;
t219 = sin(qJ(3));
t293 = qJDD(2) * t219;
t312 = t218 * t222;
t247 = t219 * t221 + t312;
t343 = t247 * qJD(4);
t348 = qJD(3) * t247 + t343;
t59 = (qJD(2) * t348 + t218 * t293) * t215 - t221 * t274;
t295 = qJD(2) * qJD(3);
t276 = t222 * t295;
t352 = t276 + t293;
t311 = t221 * t222;
t282 = t215 * t311;
t261 = qJD(2) * t282;
t303 = qJD(2) * t219;
t280 = t215 * t303;
t262 = t218 * t280;
t126 = -t261 + t262;
t137 = t247 * t215;
t128 = qJD(2) * t137;
t217 = sin(qJ(5));
t220 = cos(qJ(5));
t248 = t126 * t217 - t128 * t220;
t70 = t126 * t220 + t128 * t217;
t350 = t70 * t248;
t216 = cos(pkin(5));
t349 = t352 * t216;
t291 = qJD(3) + qJD(4);
t19 = t248 ^ 2 - t70 ^ 2;
t210 = pkin(10) + qJ(2);
t202 = sin(t210);
t203 = cos(t210);
t316 = t216 * t222;
t133 = -t202 * t316 - t203 * t219;
t241 = -t202 * t219 + t203 * t316;
t347 = -g(1) * t133 - g(2) * t241;
t304 = qJD(2) * t216;
t189 = qJD(3) + t304;
t179 = qJD(4) + t189;
t168 = qJD(5) + t179;
t297 = qJD(5) * t220;
t298 = qJD(5) * t217;
t314 = t218 * t219;
t283 = t215 * t314;
t243 = t291 * t283;
t275 = t215 * t293;
t263 = -t218 * t274 - t221 * t275 - t261 * t291;
t58 = qJD(2) * t243 + t263;
t21 = t126 * t297 + t128 * t298 + t217 * t59 + t220 * t58;
t8 = t168 * t70 - t21;
t121 = t128 * pkin(9);
t342 = pkin(7) + pkin(8);
t285 = t342 * t219;
t265 = t215 * t285;
t318 = t215 * t222;
t184 = qJD(1) * t318;
t192 = pkin(2) * t316;
t308 = qJD(2) * t192 + t184;
t94 = -qJD(2) * t265 + t308;
t82 = pkin(3) * t189 + t94;
t317 = t216 * t219;
t191 = pkin(2) * t317;
t305 = qJD(1) * t219;
t281 = t215 * t305;
t286 = t215 * t342;
t95 = t281 + (t222 * t286 + t191) * qJD(2);
t85 = t218 * t95;
t45 = t221 * t82 - t85;
t35 = -t121 + t45;
t33 = pkin(4) * t179 + t35;
t337 = pkin(9) * t126;
t87 = t221 * t95;
t46 = t218 * t82 + t87;
t36 = t46 - t337;
t185 = qJDD(2) * t216 + qJDD(3);
t277 = t219 * t295;
t272 = t216 * t292;
t294 = qJDD(1) * t215;
t309 = pkin(2) * t272 + t222 * t294;
t235 = -pkin(2) * t216 * t277 + t309;
t284 = t342 * t222;
t47 = pkin(3) * t185 + (-qJDD(2) * t285 + (-qJD(2) * t284 - t305) * qJD(3)) * t215 + t235;
t238 = -t277 + t292;
t250 = pkin(2) * t349 + pkin(7) * t274 + qJD(3) * t184 + t219 * t294;
t53 = (-pkin(7) * t277 + pkin(8) * t238) * t215 + t250;
t11 = -qJD(4) * t46 - t218 * t53 + t221 * t47;
t171 = qJDD(4) + t185;
t6 = pkin(4) * t171 + pkin(9) * t58 + t11;
t299 = qJD(4) * t221;
t300 = qJD(4) * t218;
t267 = -t218 * t47 - t221 * t53 - t299 * t82 + t300 * t95;
t7 = -pkin(9) * t59 - t267;
t1 = t220 * (qJD(5) * t33 + t7) + t217 * t6 - t36 * t298;
t214 = qJ(3) + qJ(4);
t270 = pkin(5) - t214;
t257 = -qJ(5) + t270;
t245 = cos(t257);
t204 = pkin(5) + t214;
t194 = qJ(5) + t204;
t289 = cos(t194) / 0.2e1;
t199 = pkin(3) * t222 + pkin(2);
t152 = t199 * t215;
t201 = t216 * qJD(1);
t129 = -qJD(2) * t152 + t201;
t81 = pkin(4) * t126 + t129;
t244 = sin(t257);
t288 = sin(t194) / 0.2e1;
t141 = t288 - t244 / 0.2e1;
t207 = qJ(5) + t214;
t196 = cos(t207);
t91 = -t141 * t203 - t196 * t202;
t93 = t141 * t202 - t196 * t203;
t227 = t70 * t81 - g(1) * t93 - g(2) * t91 - g(3) * (t289 - t245 / 0.2e1) - t1;
t322 = t220 * t36;
t14 = t217 * t33 + t322;
t2 = -qJD(5) * t14 - t217 * t7 + t220 * t6;
t195 = sin(t207);
t234 = t245 / 0.2e1 + t289;
t90 = t195 * t202 - t203 * t234;
t92 = -t195 * t203 - t202 * t234;
t226 = t248 * t81 - g(1) * t92 + g(2) * t90 - g(3) * (t288 + t244 / 0.2e1) + t2;
t233 = qJD(5) * t248 + t217 * t58 - t220 * t59;
t9 = -t168 * t248 + t233;
t301 = qJD(3) * t222;
t346 = t215 * (t185 * t219 + t189 * t301);
t111 = pkin(3) * t216 + t192 - t265;
t307 = pkin(7) * t318 + t191;
t124 = pkin(8) * t318 + t307;
t63 = t111 * t218 + t124 * t221;
t116 = qJD(2) * t307 + t281;
t345 = qJD(3) * t116;
t139 = t307 * qJD(3);
t209 = t215 ^ 2;
t341 = pkin(2) * t209;
t340 = pkin(2) * t215;
t339 = pkin(4) * t216;
t338 = pkin(7) * t215;
t335 = g(3) * t215;
t334 = t219 * pkin(3);
t136 = -t282 + t283;
t229 = t348 * t215;
t83 = -t282 * t291 + t243;
t29 = t136 * t297 + t137 * t298 + t217 * t229 + t220 * t83;
t80 = -t136 * t217 + t137 * t220;
t330 = t233 * t80 + t29 * t70;
t161 = qJDD(5) + t171;
t30 = qJD(5) * t80 - t217 * t83 + t220 * t229;
t79 = t136 * t220 + t137 * t217;
t329 = -t161 * t79 - t168 * t30;
t328 = t126 * t83 - t137 * t59;
t52 = t221 * t94 - t85;
t198 = pkin(3) * t221 + pkin(4);
t315 = t217 * t218;
t51 = -t218 * t94 - t87;
t40 = t51 + t337;
t41 = -t121 + t52;
t327 = -t217 * t40 - t220 * t41 + t198 * t297 + (-t218 * t298 + (t220 * t221 - t315) * qJD(4)) * pkin(3);
t313 = t218 * t220;
t326 = t217 * t41 - t220 * t40 - t198 * t298 + (-t218 * t297 + (-t217 * t221 - t313) * qJD(4)) * pkin(3);
t325 = t216 * t233;
t324 = t216 * t59;
t323 = t217 * t36;
t321 = -t136 * t171 - t179 * t229;
t320 = t128 * t126;
t319 = t209 * qJD(2) ^ 2;
t310 = t349 * t215;
t211 = t219 ^ 2;
t212 = t222 ^ 2;
t306 = t211 - t212;
t302 = qJD(3) * t219;
t296 = qJD(3) - t189;
t188 = cos(t204) / 0.2e1;
t193 = cos(t270);
t290 = t193 / 0.2e1 + t188;
t287 = sin(t204) / 0.2e1;
t279 = t215 * t302;
t62 = t111 * t221 - t124 * t218;
t266 = pkin(3) * t279;
t264 = t219 * t222 * t319;
t260 = t215 * t277;
t258 = t219 * t276;
t256 = sin(t270);
t255 = g(1) * t203 + g(2) * t202;
t254 = g(1) * t202 - g(2) * t203;
t13 = t220 * t33 - t323;
t253 = -t13 * t70 - t14 * t248;
t252 = -t21 * t79 - t248 * t30;
t251 = -t161 * t80 + t168 * t29;
t48 = -pkin(9) * t137 + t339 + t62;
t50 = -pkin(9) * t136 + t63;
t24 = -t217 * t50 + t220 * t48;
t25 = t217 * t48 + t220 * t50;
t249 = -t137 * t171 + t179 * t83;
t148 = t287 - t256 / 0.2e1;
t206 = cos(t214);
t101 = -t148 * t203 - t202 * t206;
t103 = t148 * t202 - t203 * t206;
t246 = t311 - t314;
t197 = pkin(4) * t221 + pkin(3);
t242 = -pkin(4) * t314 + t197 * t222;
t181 = qJD(3) * t192;
t112 = -qJD(3) * t265 + t181;
t113 = (-t215 * t284 - t191) * qJD(3);
t31 = t111 * t299 + t112 * t221 + t113 * t218 - t124 * t300;
t200 = t216 * qJDD(1);
t99 = pkin(3) * t260 - qJDD(2) * t152 + t200;
t32 = -qJD(4) * t63 - t112 * t218 + t113 * t221;
t232 = -g(1) * t103 - g(2) * t101 - g(3) * (t188 - t193 / 0.2e1) + t126 * t129 + t267;
t205 = sin(t214);
t100 = t202 * t205 - t203 * t290;
t102 = -t202 * t290 - t203 * t205;
t225 = -g(1) * t102 + g(2) * t100 - g(3) * (t287 + t256 / 0.2e1) - t128 * t129 + t11;
t224 = t128 * t229 - t136 * t58;
t213 = qJDD(1) - g(3);
t159 = -pkin(4) * t205 - t334;
t156 = -qJD(2) * t340 + t201;
t155 = pkin(4) * t206 + t199;
t153 = -qJDD(2) * t340 + t200;
t151 = t185 * t318;
t150 = pkin(3) * t317 - t286;
t145 = -t219 * t338 + t192;
t144 = pkin(3) * t313 + t198 * t217;
t143 = -pkin(3) * t315 + t198 * t220;
t138 = -pkin(7) * t279 + t181;
t134 = -t202 * t317 + t203 * t222;
t132 = -t202 * t222 - t203 * t317;
t122 = t242 * t216;
t115 = -pkin(7) * t280 + t308;
t109 = -t215 * (pkin(9) + t342) + (pkin(4) * t312 + t197 * t219) * t216;
t97 = pkin(3) * t280 + pkin(4) * t128;
t96 = pkin(4) * t136 - t152;
t67 = (pkin(4) * t343 + (pkin(4) * t247 + t334) * qJD(3)) * t215;
t66 = -pkin(7) * t275 + t309 - t345;
t65 = -pkin(7) * t260 + t250;
t60 = -t126 ^ 2 + t128 ^ 2;
t54 = t58 * t216;
t39 = t128 * t179 - t59;
t38 = t126 * t179 - t262 * t291 - t263;
t37 = pkin(4) * t59 + t99;
t27 = pkin(9) * t83 + t32;
t26 = -pkin(9) * t229 + t31;
t20 = t21 * t216;
t16 = t220 * t35 - t323;
t15 = -t217 * t35 - t322;
t4 = -qJD(5) * t25 - t217 * t26 + t220 * t27;
t3 = qJD(5) * t24 + t217 * t27 + t220 * t26;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t213, 0, 0, 0, 0, 0, 0, 0, 0, 0, t213, 0, 0, 0, 0, 0, 0, t151 + (-t272 + (-t189 + t304) * t302) * t215, t310 - t346, 0, t153 * t216 - g(3) + (t219 * t65 + t222 * t66 + (-t115 * t219 + t116 * t222) * qJD(3)) * t215, 0, 0, 0, 0, 0, 0, t321 + t324, t249 - t54, t224 + t328, -t11 * t136 - t137 * t267 + t99 * t216 - t229 * t45 - t46 * t83 - g(3), 0, 0, 0, 0, 0, 0, -t325 + t329, -t20 + t251, t252 + t330, t1 * t80 - t13 * t30 - t14 * t29 - t2 * t79 + t216 * t37 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t254, t255, 0, 0, (qJDD(2) * t211 + 0.2e1 * t258) * t209, 0.2e1 * (t219 * t292 - t295 * t306) * t209, t310 + t346, (qJDD(2) * t212 - 0.2e1 * t258) * t209, t151 + (t272 + (-t189 - t304) * t302) * t215, t185 * t216, -g(1) * t132 - g(2) * t134 - t139 * t189 + t145 * t185 + t216 * t66 + (-t153 * t222 + t156 * t302) * t215 + t238 * t341, g(1) * t241 - g(2) * t133 - t138 * t189 - t307 * t185 - t216 * t65 + (t153 * t219 + t156 * t301) * t215 - t352 * t341, ((-qJD(3) * t115 + qJDD(2) * t307 + t65 + (-qJD(3) * t145 + t138) * qJD(2)) * t222 + (-qJDD(2) * t145 - t345 - t66) * t219 - t255) * t215, t65 * t307 + t116 * t138 + t66 * t145 - t115 * t139 - t153 * t340 - g(1) * (-pkin(2) * t202 + t203 * t338) - g(2) * (pkin(2) * t203 + t202 * t338), -t128 * t83 - t137 * t58, -t224 + t328, -t249 - t54, t126 * t229 + t59 * t136, t321 - t324, t171 * t216, t32 * t179 + t62 * t171 + t11 * t216 - t152 * t59 + t99 * t136 - g(1) * t101 + g(2) * t103 + (t129 * t343 + (t126 * t334 + t129 * t247) * qJD(3)) * t215, -g(1) * t100 - g(2) * t102 + t128 * t266 - t129 * t83 + t137 * t99 + t152 * t58 - t171 * t63 - t179 * t31 + t216 * t267, t267 * t136 - t11 * t137 - t31 * t126 - t32 * t128 + t45 * t83 + t62 * t58 - t63 * t59 + (t46 * (-t218 * t301 - t219 * t299 - t221 * t302 - t222 * t300) - t255) * t215, -t267 * t63 + t46 * t31 + t11 * t62 + t45 * t32 - t99 * t152 + t129 * t266 - g(1) * (-t150 * t203 - t199 * t202) - g(2) * (-t150 * t202 + t199 * t203), -t21 * t80 + t248 * t29, -t252 + t330, -t20 - t251, -t233 * t79 + t30 * t70, t325 + t329, t161 * t216, -g(1) * t91 + g(2) * t93 + t24 * t161 + t4 * t168 + t2 * t216 - t233 * t96 + t81 * t30 + t37 * t79 + t67 * t70, -g(1) * t90 - g(2) * t92 - t1 * t216 - t161 * t25 - t168 * t3 - t21 * t96 - t248 * t67 - t29 * t81 + t37 * t80, -t1 * t79 + t13 * t29 - t14 * t30 - t2 * t80 + t21 * t24 - t215 * t255 + t233 * t25 + t248 * t4 - t3 * t70, t1 * t25 + t14 * t3 + t2 * t24 + t13 * t4 + t37 * t96 + t81 * t67 - g(1) * (-t109 * t203 - t155 * t202) - g(2) * (-t109 * t202 + t155 * t203); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t264, t306 * t319, (qJD(2) * t222 * t296 + t293) * t215, t264, (-t296 * t303 + t292) * t215, t185, t116 * t189 + ((-pkin(7) * t295 - g(3)) * t222 + (-pkin(7) * qJDD(2) - qJD(1) * qJD(3) - qJD(2) * t156) * t219) * t215 + t235 + t347, g(1) * t134 - g(2) * t132 + t115 * t189 + (g(3) * t219 + (pkin(7) * t302 - t156 * t222) * qJD(2)) * t215 - t250, 0, 0, t320, t60, t38, -t320, t39, t171, -t179 * t51 + (-t126 * t280 + t171 * t221 - t179 * t300) * pkin(3) + t225, t179 * t52 + (-t128 * t280 - t171 * t218 - t179 * t299) * pkin(3) + t232, (t46 + t51) * t128 + (-t45 + t52) * t126 + (-t218 * t59 + t221 * t58 + (-t126 * t221 + t128 * t218) * qJD(4)) * pkin(3), -t45 * t51 - t46 * t52 + (-t267 * t218 + t11 * t221 + (-g(3) * t222 - t129 * t303) * t215 + (-t218 * t45 + t221 * t46) * qJD(4) + t347) * pkin(3), -t350, t19, t8, t350, t9, t161, t143 * t161 + t168 * t326 - t70 * t97 + t226, -t144 * t161 - t168 * t327 + t248 * t97 + t227, t143 * t21 + t144 * t233 + t248 * t326 - t327 * t70 + t253, t1 * t144 + t2 * t143 - t81 * t97 - g(1) * (-t122 * t202 + t159 * t203) - g(2) * (t122 * t203 + t159 * t202) - t242 * t335 + t327 * t14 + t326 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t320, t60, t38, -t320, t39, t171, t179 * t46 + t225, t179 * t45 + t232, 0, 0, -t350, t19, t8, t350, t9, t161, -t15 * t168 + (-t128 * t70 + t161 * t220 - t168 * t298) * pkin(4) + t226, t16 * t168 + (t128 * t248 - t161 * t217 - t168 * t297) * pkin(4) + t227, -t15 * t248 + t16 * t70 + (t21 * t220 + t217 * t233 + (-t217 * t248 - t220 * t70) * qJD(5)) * pkin(4) + t253, -t13 * t15 - t14 * t16 + t254 * t246 * t339 + (t1 * t217 - t81 * t128 + t2 * t220 - t246 * t335 + t255 * t205 + (-t13 * t217 + t14 * t220) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t350, t19, t8, t350, t9, t161, t14 * t168 + t226, t13 * t168 + t227, 0, 0;];
tau_reg = t5;
