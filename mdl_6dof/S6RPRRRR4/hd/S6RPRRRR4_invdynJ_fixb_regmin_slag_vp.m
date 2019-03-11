% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRRR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% tau_reg [6x35]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRRR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:06:07
% EndTime: 2019-03-09 07:06:15
% DurationCPUTime: 5.54s
% Computational Cost: add. (12133->420), mult. (30675->532), div. (0->0), fcn. (25752->18), ass. (0->243)
t210 = cos(qJ(6));
t283 = qJD(6) * t210;
t206 = sin(qJ(5));
t211 = cos(qJ(5));
t204 = cos(pkin(11));
t213 = cos(qJ(3));
t296 = t213 * t204;
t203 = sin(pkin(11));
t208 = sin(qJ(3));
t303 = t203 * t208;
t152 = -t296 + t303;
t142 = t152 * qJD(1);
t153 = t203 * t213 + t204 * t208;
t143 = t153 * qJD(1);
t207 = sin(qJ(4));
t212 = cos(qJ(4));
t240 = -t142 * t207 + t212 * t143;
t241 = -t142 * t212 - t207 * t143;
t74 = t206 * t240 - t211 * t241;
t373 = t210 * t74;
t378 = t283 + t373;
t198 = qJDD(3) + qJDD(4);
t193 = qJDD(5) + t198;
t202 = qJD(3) + qJD(4);
t196 = qJD(5) + t202;
t205 = sin(qJ(6));
t285 = qJD(5) * t211;
t286 = qJD(5) * t206;
t273 = qJD(1) * t303;
t280 = t204 * qJDD(1);
t281 = t203 * qJDD(1);
t289 = qJD(3) * t213;
t275 = t204 * qJD(1) * t289 + t208 * t280 + t213 * t281;
t121 = -qJD(3) * t273 + t275;
t145 = t153 * qJD(3);
t245 = t208 * t281 - t213 * t280;
t122 = qJD(1) * t145 + t245;
t287 = qJD(4) * t212;
t288 = qJD(4) * t207;
t62 = t212 * t121 - t207 * t122 - t142 * t287 - t143 * t288;
t63 = qJD(4) * t240 + t207 * t121 + t212 * t122;
t23 = -t206 * t63 + t211 * t62 - t240 * t286 + t241 * t285;
t284 = qJD(6) * t205;
t332 = t206 * t241 + t211 * t240;
t17 = t205 * t193 + t196 * t283 + t210 * t23 - t284 * t332;
t68 = t196 * t205 + t210 * t332;
t18 = qJD(6) * t68 - t210 * t193 + t205 * t23;
t66 = -t210 * t196 + t205 * t332;
t371 = t17 * t210 - t205 * t18 - t378 * t66;
t364 = qJD(6) + t74;
t375 = t205 * t364;
t377 = -t375 * t68 + t371;
t24 = qJD(5) * t332 + t206 * t62 + t211 * t63;
t22 = qJDD(6) + t24;
t20 = t210 * t22;
t324 = t66 * t332;
t376 = -t364 * t375 + t20 + t324;
t317 = t196 * t74;
t360 = t23 + t317;
t201 = pkin(11) + qJ(3);
t197 = qJ(4) + t201;
t188 = qJ(5) + t197;
t181 = sin(t188);
t209 = sin(qJ(1));
t214 = cos(qJ(1));
t248 = g(1) * t214 + g(2) * t209;
t374 = t248 * t181;
t15 = t17 * t205;
t370 = t378 * t68 + t15;
t19 = t205 * t22;
t323 = t68 * t332;
t70 = t364 * t283;
t369 = t364 * t373 + t19 - t323 + t70;
t322 = pkin(7) + qJ(2);
t167 = t322 * t203;
t154 = qJD(1) * t167;
t168 = t322 * t204;
t155 = qJD(1) * t168;
t239 = t154 * t208 - t155 * t213;
t99 = -pkin(8) * t142 - t239;
t96 = t212 * t99;
t347 = -t213 * t154 - t155 * t208;
t98 = -pkin(8) * t143 + t347;
t97 = qJD(3) * pkin(3) + t98;
t244 = -t207 * t97 - t96;
t354 = pkin(9) * t241;
t51 = -t244 + t354;
t315 = t206 * t51;
t94 = t207 * t99;
t264 = t212 * t97 - t94;
t355 = pkin(9) * t240;
t50 = t264 - t355;
t48 = pkin(4) * t202 + t50;
t29 = t211 * t48 - t315;
t27 = -pkin(5) * t196 - t29;
t368 = t27 * t74;
t365 = t332 * t74;
t318 = t196 * t332;
t336 = -t24 + t318;
t361 = t332 ^ 2 - t74 ^ 2;
t356 = pkin(5) * t332;
t372 = pkin(10) * t74 + t356;
t172 = g(3) * t181;
t182 = cos(t188);
t282 = qJD(1) * qJD(2);
t333 = t322 * qJDD(1) + t282;
t133 = t333 * t203;
t134 = t333 * t204;
t252 = -t213 * t133 - t208 * t134;
t57 = qJDD(3) * pkin(3) - pkin(8) * t121 + qJD(3) * t239 + t252;
t242 = -t208 * t133 + t213 * t134;
t61 = -pkin(8) * t122 + qJD(3) * t347 + t242;
t227 = qJD(4) * t244 - t207 * t61 + t212 * t57;
t12 = pkin(4) * t198 - pkin(9) * t62 + t227;
t334 = (qJD(4) * t97 + t61) * t212 + t207 * t57 - t99 * t288;
t13 = -pkin(9) * t63 + t334;
t335 = (qJD(5) * t48 + t13) * t211 + t206 * t12 - t51 * t286;
t187 = -pkin(2) * t204 - pkin(1);
t162 = t187 * qJD(1) + qJD(2);
t125 = t142 * pkin(3) + t162;
t84 = -pkin(4) * t241 + t125;
t359 = t248 * t182 + t74 * t84 + t172 - t335;
t311 = t211 * t51;
t30 = t206 * t48 + t311;
t339 = qJD(5) * t30 - t211 * t12 + t206 * t13;
t3 = -pkin(5) * t193 + t339;
t327 = g(3) * t182;
t367 = t3 + t327;
t350 = t364 * t332;
t306 = t240 * t202;
t362 = -t63 + t306;
t28 = pkin(10) * t196 + t30;
t39 = pkin(5) * t74 - pkin(10) * t332 + t84;
t8 = -t205 * t28 + t210 * t39;
t343 = t210 * t374 + t27 * t284 - t8 * t332;
t9 = t205 * t39 + t210 * t28;
t342 = t367 * t205 + t27 * t283 + t9 * t332;
t341 = -t332 * t84 - t327 - t339 + t374;
t357 = pkin(4) * t240;
t189 = pkin(4) * t206 + pkin(10);
t353 = (qJD(6) * t189 + t357 + t372) * t364;
t191 = pkin(3) * t212 + pkin(4);
t299 = t207 * t211;
t291 = pkin(3) * t299 + t206 * t191;
t139 = pkin(10) + t291;
t90 = pkin(3) * t143 + t357;
t352 = (qJD(6) * t139 + t372 + t90) * t364;
t351 = (t364 * pkin(10) + t356) * t364;
t305 = t241 * t202;
t349 = t62 - t305;
t348 = t240 * t241;
t251 = -t213 * t167 - t168 * t208;
t109 = -pkin(8) * t153 + t251;
t292 = -t208 * t167 + t213 * t168;
t110 = -pkin(8) * t152 + t292;
t293 = t207 * t109 + t212 * t110;
t346 = qJ(2) * qJDD(1);
t345 = t240 ^ 2 - t241 ^ 2;
t185 = sin(t197);
t186 = cos(t197);
t338 = -g(3) * t186 - t125 * t240 + t248 * t185 + t227;
t337 = g(3) * t185 - t125 * t241 + t248 * t186 - t334;
t271 = pkin(10) * t193 + qJD(6) * t39 + t335;
t124 = -t152 * t207 + t153 * t212;
t254 = t212 * t109 - t110 * t207;
t54 = -pkin(9) * t124 + t254;
t123 = t212 * t152 + t153 * t207;
t55 = -pkin(9) * t123 + t293;
t38 = t206 * t54 + t211 * t55;
t228 = -t167 * t289 + qJD(2) * t296 + (-qJD(2) * t203 - qJD(3) * t168) * t208;
t88 = -pkin(8) * t145 + t228;
t144 = t152 * qJD(3);
t217 = -t153 * qJD(2) - qJD(3) * t292;
t89 = pkin(8) * t144 + t217;
t234 = t109 * t287 - t110 * t288 + t207 * t89 + t212 * t88;
t81 = qJD(4) * t124 - t207 * t144 + t212 * t145;
t35 = -pkin(9) * t81 + t234;
t226 = -qJD(4) * t293 - t207 * t88 + t212 * t89;
t80 = -qJD(4) * t123 - t212 * t144 - t207 * t145;
t36 = -pkin(9) * t80 + t226;
t37 = t206 * t55 - t211 * t54;
t4 = -qJD(5) * t37 + t206 * t36 + t211 * t35;
t82 = t211 * t123 + t124 * t206;
t40 = -qJD(5) * t82 - t206 * t81 + t211 * t80;
t83 = -t123 * t206 + t124 * t211;
t127 = pkin(3) * t152 + t187;
t92 = pkin(4) * t123 + t127;
t44 = pkin(5) * t82 - pkin(10) * t83 + t92;
t331 = -t38 * t22 - (qJD(6) * t44 + t4) * t364 + t27 * t40 - t271 * t82 + t3 * t83;
t330 = pkin(3) * t145;
t326 = t27 * t83;
t325 = t44 * t22;
t319 = t212 * t98 - t94;
t316 = t205 * t68;
t300 = t206 * t207;
t263 = -t207 * t98 - t96;
t52 = t263 - t354;
t53 = t319 - t355;
t310 = t206 * t52 + t211 * t53 - t191 * t285 - (-t207 * t286 + (t211 * t212 - t300) * qJD(4)) * pkin(3);
t309 = -t206 * t53 + t211 * t52 + t191 * t286 + (t207 * t285 + (t206 * t212 + t299) * qJD(4)) * pkin(3);
t307 = qJDD(1) * pkin(1);
t302 = t205 * t209;
t301 = t205 * t214;
t298 = t209 * t210;
t297 = t210 * t214;
t290 = t203 ^ 2 + t204 ^ 2;
t276 = t83 * t284;
t156 = t187 * qJDD(1) + qJDD(2);
t100 = t122 * pkin(3) + t156;
t47 = pkin(4) * t63 + t100;
t7 = pkin(5) * t24 - pkin(10) * t23 + t47;
t270 = qJD(6) * t28 - t7;
t257 = t290 * qJD(1) ^ 2;
t250 = 0.2e1 * t290;
t31 = t206 * t50 + t311;
t249 = pkin(4) * t286 - t31;
t69 = pkin(4) * t81 + t330;
t247 = g(1) * t209 - g(2) * t214;
t246 = t22 * t83 + t364 * t40;
t238 = t20 - (t205 * t74 + t284) * t364;
t237 = -t271 + t172;
t236 = -pkin(3) * t300 + t191 * t211;
t233 = -pkin(10) * t22 + t29 * t364 + t368;
t231 = -t247 - t307;
t192 = qJDD(2) - t307;
t230 = -t192 - t231;
t229 = -t139 * t22 + t310 * t364 + t368;
t32 = t211 * t50 - t315;
t224 = -t189 * t22 + t368 + (-pkin(4) * t285 + t32) * t364;
t222 = t250 * t282 - t248;
t195 = cos(t201);
t194 = sin(t201);
t190 = -pkin(4) * t211 - pkin(5);
t138 = -pkin(5) - t236;
t132 = t182 * t297 + t302;
t131 = -t182 * t301 + t298;
t130 = -t182 * t298 + t301;
t129 = t182 * t302 + t297;
t41 = qJD(5) * t83 + t206 * t80 + t211 * t81;
t10 = pkin(5) * t41 - pkin(10) * t40 + t69;
t6 = t210 * t7;
t5 = qJD(5) * t38 + t206 * t35 - t211 * t36;
t1 = [qJDD(1), t247, t248, t230 * t204, -t230 * t203, t250 * t346 + t222 (-t192 + t247) * pkin(1) + (t290 * t346 + t222) * qJ(2), t121 * t153 - t143 * t144, -t121 * t152 - t122 * t153 + t142 * t144 - t143 * t145, -qJD(3) * t144 + qJDD(3) * t153, -qJD(3) * t145 - qJDD(3) * t152, 0, qJD(3) * t217 + qJDD(3) * t251 + t187 * t122 + t162 * t145 + t156 * t152 + t195 * t247, -t228 * qJD(3) - t292 * qJDD(3) + t187 * t121 - t162 * t144 + t156 * t153 - t247 * t194, t124 * t62 + t240 * t80, -t123 * t62 - t124 * t63 - t240 * t81 + t241 * t80, t124 * t198 + t202 * t80, -t123 * t198 - t202 * t81, 0, t100 * t123 + t125 * t81 + t127 * t63 + t186 * t247 + t198 * t254 + t202 * t226 - t241 * t330, t100 * t124 + t125 * t80 + t127 * t62 - t247 * t185 - t293 * t198 - t234 * t202 + t240 * t330, t23 * t83 + t332 * t40, -t23 * t82 - t24 * t83 - t332 * t41 - t40 * t74, t193 * t83 + t196 * t40, -t193 * t82 - t196 * t41, 0, t182 * t247 - t193 * t37 - t196 * t5 + t24 * t92 + t41 * t84 + t47 * t82 + t69 * t74, -t181 * t247 - t193 * t38 - t196 * t4 + t23 * t92 + t332 * t69 + t40 * t84 + t47 * t83, -t68 * t276 + (t17 * t83 + t40 * t68) * t210 (-t210 * t66 - t316) * t40 + (-t15 - t18 * t210 + (t205 * t66 - t210 * t68) * qJD(6)) * t83, t17 * t82 + t210 * t246 - t276 * t364 + t41 * t68, -t18 * t82 - t205 * t246 - t41 * t66 - t70 * t83, t22 * t82 + t364 * t41, -g(1) * t130 - g(2) * t132 + t37 * t18 + t8 * t41 + t5 * t66 + t6 * t82 + (t10 * t364 + t325 + (-t28 * t82 - t364 * t38 + t326) * qJD(6)) * t210 + t331 * t205, -g(1) * t129 - g(2) * t131 + t37 * t17 - t9 * t41 + t5 * t68 + (-(-qJD(6) * t38 + t10) * t364 - t325 + t270 * t82 - qJD(6) * t326) * t205 + t331 * t210; 0, 0, 0, -t280, t281, -t257, -qJ(2) * t257 + qJDD(2) + t231, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t143 + t245 (-t142 - t273) * qJD(3) + t275, 0, 0, 0, 0, 0, t63 + t306, t62 + t305, 0, 0, 0, 0, 0, t24 + t318, t23 - t317, 0, 0, 0, 0, 0, t238 - t324, -t210 * t364 ^ 2 - t19 - t323; 0, 0, 0, 0, 0, 0, 0, t143 * t142, -t142 ^ 2 + t143 ^ 2 (t142 - t273) * qJD(3) + t275, -t245, qJDD(3), -g(3) * t195 - t162 * t143 + t194 * t248 + t252, g(3) * t194 + t162 * t142 + t248 * t195 - t242, -t348, t345, t349, t362, t198, -t263 * t202 + (t143 * t241 + t198 * t212 - t202 * t288) * pkin(3) + t338, t319 * t202 + (-t143 * t240 - t198 * t207 - t202 * t287) * pkin(3) + t337, t365, t361, t360, t336, t193, t236 * t193 - t309 * t196 - t90 * t74 + t341, -t291 * t193 + t310 * t196 - t332 * t90 + t359, t370, -t316 * t364 + t371, t369, t238 + t324, -t350, t138 * t18 + t309 * t66 + (-t367 - t352) * t210 + t229 * t205 + t343, t138 * t17 + t309 * t68 + t229 * t210 + (-t374 + t352) * t205 + t342; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t348, t345, t349, t362, t198, -t202 * t244 + t338, t264 * t202 + t337, t365, t361, t360, t336, t193, t196 * t31 + (t193 * t211 - t196 * t286 - t240 * t74) * pkin(4) + t341, t196 * t32 + (-t193 * t206 - t196 * t285 - t240 * t332) * pkin(4) + t359, t370, t377, t369, t376, -t350, t190 * t18 + t249 * t66 + (-t367 - t353) * t210 + t224 * t205 + t343, t190 * t17 + t249 * t68 + t224 * t210 + (-t374 + t353) * t205 + t342; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t365, t361, t360, t336, t193, t196 * t30 + t341, t196 * t29 + t359, t370, t377, t369, t376, -t350, -pkin(5) * t18 - t30 * t66 + t233 * t205 + (-t367 - t351) * t210 + t343, -pkin(5) * t17 - t30 * t68 + t233 * t210 + (-t374 + t351) * t205 + t342; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68 * t66, -t66 ^ 2 + t68 ^ 2, t364 * t66 + t17, t364 * t68 - t18, t22, -g(1) * t131 + g(2) * t129 + t205 * t237 - t27 * t68 - t28 * t283 + t364 * t9 + t6, g(1) * t132 - g(2) * t130 + t205 * t270 + t210 * t237 + t27 * t66 + t364 * t8;];
tau_reg  = t1;
