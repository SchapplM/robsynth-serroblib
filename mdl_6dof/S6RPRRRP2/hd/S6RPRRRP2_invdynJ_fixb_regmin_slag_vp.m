% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% tau_reg [6x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRRP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:01:24
% EndTime: 2019-03-09 06:01:31
% DurationCPUTime: 4.01s
% Computational Cost: add. (5287->430), mult. (11171->570), div. (0->0), fcn. (7611->14), ass. (0->234)
t199 = sin(qJ(4));
t203 = cos(qJ(3));
t285 = qJD(1) * t203;
t263 = t199 * t285;
t329 = pkin(8) + pkin(9);
t264 = qJD(4) * t329;
t196 = sin(pkin(10));
t174 = pkin(1) * t196 + pkin(7);
t153 = t174 * qJD(1);
t200 = sin(qJ(3));
t103 = qJD(2) * t203 - t200 * t153;
t242 = pkin(3) * t200 - pkin(8) * t203;
t139 = t242 * qJD(1);
t202 = cos(qJ(4));
t308 = t202 * t103 + t199 * t139;
t355 = pkin(9) * t263 - t199 * t264 - t308;
t119 = t202 * t139;
t293 = t202 * t203;
t234 = pkin(4) * t200 - pkin(9) * t293;
t354 = qJD(1) * t234 - t199 * t103 + t202 * t264 + t119;
t280 = qJD(3) * t202;
t286 = qJD(1) * t200;
t130 = -t199 * t286 + t280;
t282 = qJD(3) * t199;
t131 = t202 * t286 + t282;
t198 = sin(qJ(5));
t325 = cos(qJ(5));
t226 = t198 * t130 + t325 * t131;
t330 = t226 ^ 2;
t77 = -t325 * t130 + t131 * t198;
t74 = t77 ^ 2;
t353 = -t74 + t330;
t298 = t198 * t199;
t225 = t325 * t202 - t298;
t335 = qJD(4) + qJD(5);
t257 = t325 * qJD(5);
t336 = t325 * qJD(4) + t257;
t311 = -t336 * t202 + t225 * t285 + t335 * t298;
t133 = t198 * t202 + t325 * t199;
t83 = t335 * t133;
t310 = -t133 * t285 + t83;
t352 = t226 * t77;
t351 = t77 * qJ(6);
t148 = t174 * qJDD(1);
t350 = qJD(2) * qJD(3) + t148;
t172 = -qJD(4) + t285;
t192 = qJ(1) + pkin(10);
t183 = sin(t192);
t184 = cos(t192);
t241 = g(1) * t184 + g(2) * t183;
t229 = t241 * t200;
t216 = g(3) * t203 - t229;
t279 = qJD(3) * t203;
t265 = t153 * t279 + t200 * t350;
t307 = qJDD(3) * pkin(3);
t233 = t265 - t307;
t270 = t203 * qJDD(2);
t59 = t233 - t270;
t349 = -qJD(4) * pkin(8) * t172 + t216 + t59;
t261 = t199 * t279;
t276 = qJD(4) * t202;
t348 = t200 * t276 + t261;
t347 = t103 * qJD(3);
t274 = qJD(1) * qJD(3);
t256 = t203 * t274;
t271 = t200 * qJDD(1);
t346 = qJD(3) * qJD(4) + t256 + t271;
t163 = -qJD(5) + t172;
t277 = qJD(4) * t200;
t255 = qJD(1) * t277;
t244 = t199 * t346 + t202 * t255;
t222 = t202 * qJDD(3) - t244;
t275 = qJD(5) * t198;
t71 = (qJDD(3) - t255) * t199 + t346 * t202;
t21 = -t130 * t257 + t131 * t275 - t198 * t222 - t325 * t71;
t345 = -t163 * t77 - t21;
t58 = qJDD(3) * pkin(8) + t200 * qJDD(2) + t203 * t148 + t347;
t197 = cos(pkin(10));
t175 = -pkin(1) * t197 - pkin(2);
t122 = -pkin(3) * t203 - pkin(8) * t200 + t175;
t142 = t242 * qJD(3);
t73 = qJD(1) * t142 + qJDD(1) * t122;
t95 = t122 * qJD(1);
t269 = t199 * t73 + t202 * t58 + t95 * t276;
t278 = qJD(4) * t199;
t104 = t200 * qJD(2) + t203 * t153;
t94 = qJD(3) * pkin(8) + t104;
t228 = -t94 * t278 + t269;
t13 = pkin(9) * t222 + t228;
t50 = -t199 * t94 + t202 * t95;
t43 = -pkin(9) * t131 + t50;
t36 = -pkin(4) * t172 + t43;
t51 = t199 * t95 + t202 * t94;
t44 = pkin(9) * t130 + t51;
t186 = t203 * qJDD(1);
t126 = t200 * t274 + qJDD(4) - t186;
t68 = t202 * t73;
t8 = t126 * pkin(4) - t71 * pkin(9) - qJD(4) * t51 - t199 * t58 + t68;
t254 = -t325 * t13 - t198 * t8 - t36 * t257 + t44 * t275;
t195 = qJ(4) + qJ(5);
t189 = cos(t195);
t300 = t189 * t200;
t93 = -qJD(3) * pkin(3) - t103;
t72 = -pkin(4) * t130 + t93;
t188 = sin(t195);
t299 = t189 * t203;
t89 = -t183 * t299 + t184 * t188;
t91 = t183 * t188 + t184 * t299;
t344 = g(1) * t91 - g(2) * t89 + g(3) * t300 + t72 * t77 + t254;
t37 = pkin(5) * t77 + qJD(6) + t72;
t343 = t226 * t37;
t341 = -t104 + (-t263 + t278) * pkin(4);
t340 = qJ(6) * t226;
t292 = qJDD(2) - g(3);
t339 = t203 * t292;
t156 = t329 * t199;
t157 = t329 * t202;
t289 = -t198 * t156 + t325 * t157;
t338 = t289 * qJD(5) + t198 * t355 + t354 * t325;
t337 = -t156 * t257 - t157 * t275 - t354 * t198 + t325 * t355;
t322 = g(3) * t200;
t301 = t188 * t203;
t88 = t183 * t301 + t184 * t189;
t90 = t183 * t189 - t184 * t301;
t334 = -g(1) * t90 + g(2) * t88 + t188 * t322;
t42 = t325 * t44;
t17 = t198 * t36 + t42;
t214 = -qJD(5) * t17 - t198 * t13 + t325 * t8;
t333 = -t72 * t226 + t214 + t334;
t22 = qJD(5) * t226 + t198 * t71 - t325 * t222;
t332 = -t163 * t226 - t22;
t260 = t199 * t277;
t220 = t202 * t279 - t260;
t294 = t202 * t126;
t331 = -t172 * t220 + t200 * t294;
t40 = t198 * t44;
t16 = t325 * t36 - t40;
t11 = t16 - t340;
t9 = -pkin(5) * t163 + t11;
t326 = t11 - t9;
t190 = t202 * pkin(4);
t320 = pkin(3) + t190;
t143 = pkin(4) * t199 + pkin(5) * t188;
t319 = pkin(7) + t143;
t318 = -qJ(6) * t310 + qJD(6) * t225 + t337;
t317 = -pkin(5) * t286 + qJ(6) * t311 - t133 * qJD(6) - t338;
t106 = t225 * t200;
t243 = t325 * t279;
t45 = t198 * t261 + t83 * t200 - t202 * t243;
t316 = -t106 * t22 + t45 * t77;
t105 = t133 * t200;
t121 = qJDD(5) + t126;
t297 = t199 * t200;
t46 = t199 * t243 - t198 * t260 - t275 * t297 + (t198 * t279 + t336 * t200) * t202;
t315 = -t105 * t121 + t46 * t163;
t314 = t325 * t43 - t40;
t108 = t202 * t122;
t295 = t200 * t202;
t302 = t174 * t199;
t60 = -pkin(9) * t295 + t108 + (-pkin(4) - t302) * t203;
t137 = t174 * t293;
t288 = t199 * t122 + t137;
t70 = -pkin(9) * t297 + t288;
t312 = t198 * t60 + t325 * t70;
t309 = t71 * t199;
t306 = t130 * t172;
t305 = t131 * t172;
t304 = t143 * t203;
t303 = t172 * t202;
t296 = t199 * t203;
t291 = t122 * t276 + t199 * t142;
t281 = qJD(3) * t200;
t290 = t202 * t142 + t281 * t302;
t111 = pkin(4) * t297 + t200 * t174;
t144 = pkin(5) * t189 + t190;
t193 = t200 ^ 2;
t287 = -t203 ^ 2 + t193;
t154 = qJD(1) * t175;
t283 = qJD(3) * t130;
t84 = pkin(4) * t348 + t174 * t279;
t262 = t172 * t282;
t253 = -t198 * t43 - t42;
t251 = -t198 * t70 + t325 * t60;
t249 = t203 * t21 + t226 * t281;
t248 = -qJD(4) * t95 - t58;
t247 = t131 * t281 - t71 * t203;
t246 = t172 * t174 + t94;
t245 = -t325 * t156 - t157 * t198;
t240 = g(1) * t183 - g(2) * t184;
t201 = sin(qJ(1));
t204 = cos(qJ(1));
t239 = g(1) * t201 - g(2) * t204;
t238 = -t105 * t21 + t226 * t46;
t237 = -t106 * t121 - t163 * t45;
t138 = pkin(3) + t144;
t191 = -qJ(6) - t329;
t236 = t138 * t203 - t191 * t200;
t232 = t239 * pkin(1);
t231 = pkin(2) + t236;
t230 = t203 * t22 - t77 * t281;
t30 = t234 * qJD(3) + (-t137 + (pkin(9) * t200 - t122) * t199) * qJD(4) + t290;
t33 = (-t200 * t280 - t203 * t278) * t174 - t348 * pkin(9) + t291;
t227 = t198 * t30 + t60 * t257 - t70 * t275 + t325 * t33;
t224 = -t199 * t126 + t172 * t276;
t221 = -qJD(1) * t154 + t241;
t219 = -pkin(8) * t126 - t172 * t93;
t217 = 0.2e1 * t154 * qJD(3) - t174 * qJDD(3);
t205 = qJD(3) ^ 2;
t215 = -0.2e1 * qJDD(1) * t175 - t174 * t205 + t240;
t213 = -t312 * qJD(5) - t198 * t33 + t325 * t30;
t210 = -pkin(4) * t222 + t233;
t208 = t22 * pkin(5) + qJDD(6) + t210;
t206 = qJD(1) ^ 2;
t181 = t325 * pkin(4) + pkin(5);
t147 = qJDD(3) * t203 - t200 * t205;
t146 = qJDD(3) * t200 + t203 * t205;
t101 = t183 * t199 + t184 * t293;
t100 = t183 * t202 - t184 * t296;
t99 = -t183 * t293 + t184 * t199;
t98 = t183 * t296 + t184 * t202;
t63 = qJ(6) * t225 + t289;
t62 = -qJ(6) * t133 + t245;
t32 = t210 - t270;
t24 = -qJ(6) * t105 + t312;
t23 = -pkin(5) * t203 - qJ(6) * t106 + t251;
t15 = t314 - t340;
t14 = t253 + t351;
t12 = t17 - t351;
t5 = t208 - t270;
t4 = -qJ(6) * t46 - qJD(6) * t105 + t227;
t3 = pkin(5) * t281 + t45 * qJ(6) - t106 * qJD(6) + t213;
t2 = -qJ(6) * t22 - qJD(6) * t77 - t254;
t1 = t121 * pkin(5) + t21 * qJ(6) - qJD(6) * t226 + t214;
t6 = [qJDD(1), t239, g(1) * t204 + g(2) * t201 (t196 ^ 2 + t197 ^ 2) * pkin(1) ^ 2 * qJDD(1) + t232, qJDD(1) * t193 + 0.2e1 * t200 * t256, 0.2e1 * t200 * t186 - 0.2e1 * t287 * t274, t146, t147, 0, t200 * t217 + t203 * t215, -t200 * t215 + t203 * t217, t131 * t220 + t295 * t71 (t130 * t202 - t131 * t199) * t279 + (t202 * t222 - t309 + (-t130 * t199 - t131 * t202) * qJD(4)) * t200, t247 + t331 (-t222 + t262) * t203 + (t224 + t283) * t200, -t126 * t203 - t172 * t281 -(-t122 * t278 + t290) * t172 + t108 * t126 - g(1) * t99 - g(2) * t101 + (-t174 * t283 - t68 + t246 * t276 + (qJD(3) * t93 - t126 * t174 - t248) * t199) * t203 + (t50 * qJD(3) - t174 * t222 + t59 * t199 + t276 * t93) * t200, t291 * t172 - t288 * t126 - g(1) * t98 - g(2) * t100 + (-t246 * t278 + (t131 * t174 + t202 * t93) * qJD(3) + t269) * t203 + (-t93 * t278 + t174 * t71 + t59 * t202 + (-t174 * t303 - t51) * qJD(3)) * t200, -t106 * t21 - t226 * t45, -t238 + t316, -t237 + t249, t230 + t315, -t121 * t203 - t163 * t281, -g(1) * t89 - g(2) * t91 + t32 * t105 + t111 * t22 + t121 * t251 + t16 * t281 - t163 * t213 - t203 * t214 + t72 * t46 + t84 * t77, -g(1) * t88 - g(2) * t90 + t32 * t106 - t111 * t21 - t121 * t312 + t163 * t227 - t17 * t281 - t203 * t254 + t226 * t84 - t72 * t45, -t1 * t106 - t105 * t2 - t12 * t46 + t200 * t240 + t21 * t23 - t22 * t24 - t226 * t3 - t4 * t77 + t45 * t9, t2 * t24 + t12 * t4 + t1 * t23 + t9 * t3 + t5 * (pkin(5) * t105 + t111) + t37 * (pkin(5) * t46 + t84) + t232 + (-g(1) * t319 - g(2) * t231) * t184 + (g(1) * t231 - g(2) * t319) * t183; 0, 0, 0, t292, 0, 0, 0, 0, 0, t147, -t146, 0, 0, 0, 0, 0 (t222 + t262) * t203 + (t224 - t283) * t200, t247 - t331, 0, 0, 0, 0, 0, -t230 + t315, t237 + t249, t238 + t316, -t1 * t105 + t106 * t2 - t12 * t45 - t203 * t5 + t281 * t37 - t46 * t9 - g(3); 0, 0, 0, 0, -t203 * t206 * t200, t287 * t206, t271, t186, qJDD(3), qJD(3) * t104 + t200 * t221 - t265 + t339, t347 + (qJD(3) * t153 - t292) * t200 + (t221 - t350) * t203, -t131 * t303 + t309 (t71 - t306) * t202 + (t222 + t305) * t199 (-t131 * t200 + t172 * t293) * qJD(1) - t224, t172 * t278 + t294 + (-t130 * t200 - t172 * t296) * qJD(1), t172 * t286, -pkin(3) * t244 + t119 * t172 - t50 * t286 + t104 * t130 + (-t103 * t172 + t219) * t199 + (t307 - t349) * t202, -pkin(3) * t71 - t104 * t131 - t308 * t172 + t199 * t349 + t219 * t202 + t51 * t286, -t21 * t133 - t226 * t311, -t133 * t22 - t21 * t225 - t226 * t310 + t311 * t77, t133 * t121 + t311 * t163 - t226 * t286, t121 * t225 + t310 * t163 + t77 * t286, t163 * t286, -g(3) * t299 + t121 * t245 - t16 * t286 + t338 * t163 - t22 * t320 - t225 * t32 + t241 * t300 + t310 * t72 + t341 * t77, -t289 * t121 + t32 * t133 + t337 * t163 + t17 * t286 + t216 * t188 + t21 * t320 + t341 * t226 - t311 * t72, -t1 * t133 - t12 * t310 + t2 * t225 - t203 * t241 + t21 * t62 - t22 * t63 - t226 * t317 + t311 * t9 - t318 * t77 - t322, t2 * t63 + t1 * t62 + t5 * (-pkin(5) * t225 - t320) - g(3) * t236 + t317 * t9 + (pkin(5) * t310 + t341) * t37 + t318 * t12 + t241 * (t138 * t200 + t191 * t203); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131 * t130, -t130 ^ 2 + t131 ^ 2, t71 + t306, t222 - t305, t126, -t94 * t276 - g(1) * t100 + g(2) * t98 - t131 * t93 - t172 * t51 + t68 + (t248 + t322) * t199, g(1) * t101 - g(2) * t99 + g(3) * t295 - t130 * t93 - t172 * t50 - t228, t352, t353, t345, t332, t121, t253 * t163 + (t325 * t121 - t131 * t77 + t163 * t275) * pkin(4) + t333, -t314 * t163 + (-t198 * t121 - t131 * t226 + t163 * t257) * pkin(4) + t344, t12 * t226 + t14 * t226 + t15 * t77 + t181 * t21 - t9 * t77 + (-t198 * t22 + (t198 * t226 - t325 * t77) * qJD(5)) * pkin(4), t1 * t181 - t12 * t15 - t9 * t14 - pkin(5) * t343 - g(1) * (t144 * t183 - t184 * t304) - g(2) * (-t144 * t184 - t183 * t304) + t143 * t322 + (-t37 * t131 + t2 * t198 + (t325 * t12 - t198 * t9) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t352, t353, t345, t332, t121, -t17 * t163 + t333, -t16 * t163 + t344, pkin(5) * t21 + t326 * t77, -t326 * t12 + (t1 + t334 - t343) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74 - t330, t12 * t77 + t226 * t9 + t208 - t229 - t339;];
tau_reg  = t6;
