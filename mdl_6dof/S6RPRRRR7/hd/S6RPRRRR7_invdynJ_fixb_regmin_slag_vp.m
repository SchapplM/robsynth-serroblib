% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% tau_reg [6x34]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRRR7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR7_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:17:59
% EndTime: 2019-03-09 07:18:10
% DurationCPUTime: 4.24s
% Computational Cost: add. (7013->384), mult. (14258->502), div. (0->0), fcn. (10589->14), ass. (0->229)
t203 = cos(qJ(6));
t279 = qJD(6) * t203;
t200 = sin(qJ(4));
t205 = cos(qJ(4));
t206 = cos(qJ(3));
t287 = qJD(1) * t206;
t201 = sin(qJ(3));
t288 = qJD(1) * t201;
t127 = t200 * t288 - t205 * t287;
t128 = -t200 * t287 - t205 * t288;
t199 = sin(qJ(5));
t204 = cos(qJ(5));
t86 = t127 * t199 + t204 * t128;
t346 = t203 * t86;
t355 = t279 - t346;
t190 = qJDD(3) + qJDD(4);
t182 = qJDD(5) + t190;
t191 = qJD(3) + qJD(4);
t183 = qJD(5) + t191;
t198 = sin(qJ(6));
t237 = t204 * t127 - t128 * t199;
t280 = qJD(6) * t198;
t281 = qJD(5) * t204;
t282 = qJD(5) * t199;
t275 = t206 * qJDD(1);
t276 = t201 * qJDD(1);
t135 = t200 * t206 + t201 * t205;
t97 = t191 * t135;
t69 = -qJD(1) * t97 - t200 * t276 + t205 * t275;
t284 = qJD(4) * t200;
t269 = t201 * t284;
t278 = qJD(1) * qJD(3);
t340 = t201 * t278 - t275;
t231 = -qJD(1) * t269 - t200 * t340;
t245 = t191 * t206;
t70 = (qJD(1) * t245 + t276) * t205 + t231;
t30 = t127 * t282 + t128 * t281 - t199 * t70 + t204 * t69;
t18 = t198 * t182 + t183 * t279 + t203 * t30 + t237 * t280;
t73 = t183 * t198 - t203 * t237;
t19 = qJD(6) * t73 - t203 * t182 + t198 * t30;
t345 = qJD(6) - t86;
t354 = t198 * t345;
t71 = -t203 * t183 - t198 * t237;
t1 = t18 * t203 - t198 * t19 - t354 * t73 - t355 * t71;
t221 = qJD(5) * t237 - t199 * t69 - t204 * t70;
t29 = qJDD(6) - t221;
t27 = t203 * t29;
t7 = -t237 * t71 - t345 * t354 + t27;
t207 = cos(qJ(1));
t189 = g(2) * t207;
t202 = sin(qJ(1));
t321 = g(1) * t202;
t342 = t189 - t321;
t136 = -t200 * t201 + t205 * t206;
t333 = -t135 * t199 + t204 * t136;
t236 = t204 * t135 + t136 * t199;
t286 = qJD(3) * t201;
t98 = -t200 * t286 + t205 * t245 - t269;
t46 = qJD(5) * t236 + t199 * t98 + t204 * t97;
t353 = t18 * t333 - t46 * t73;
t352 = t182 * t333 - t183 * t46;
t16 = t18 * t198;
t9 = t355 * t73 + t16;
t26 = t198 * t29;
t76 = t345 * t279;
t8 = t237 * t73 - t345 * t346 + t26 + t76;
t350 = pkin(5) * t237;
t208 = -pkin(1) - pkin(7);
t157 = t208 * qJD(1) + qJD(2);
t116 = -pkin(8) * t288 + t157 * t201;
t104 = t205 * t116;
t117 = -pkin(8) * t287 + t206 * t157;
t106 = qJD(3) * pkin(3) + t117;
t238 = -t106 * t200 - t104;
t322 = pkin(9) * t128;
t61 = -t238 + t322;
t310 = t199 * t61;
t121 = t127 * pkin(9);
t103 = t200 * t116;
t249 = t205 * t106 - t103;
t60 = t121 + t249;
t58 = pkin(4) * t191 + t60;
t36 = t204 * t58 - t310;
t34 = -pkin(5) * t183 - t36;
t320 = t34 * t86;
t197 = qJ(3) + qJ(4);
t187 = qJ(5) + t197;
t173 = cos(t187);
t272 = t173 * t321;
t155 = t208 * qJDD(1) + qJDD(2);
t139 = t206 * t155;
t79 = qJDD(3) * pkin(3) + pkin(8) * t340 - t157 * t286 + t139;
t285 = qJD(3) * t206;
t267 = t206 * t278;
t341 = t267 + t276;
t82 = -pkin(8) * t341 + t155 * t201 + t157 * t285;
t222 = qJD(4) * t238 - t200 * t82 + t205 * t79;
t23 = pkin(4) * t190 - pkin(9) * t69 + t222;
t328 = (qJD(4) * t106 + t82) * t205 - t116 * t284 + t200 * t79;
t24 = -pkin(9) * t70 + t328;
t307 = t204 * t61;
t37 = t199 * t58 + t307;
t330 = qJD(5) * t37 + t199 * t24 - t204 * t23;
t4 = -pkin(5) * t182 + t330;
t348 = t4 + t272;
t316 = t237 * t86;
t317 = t345 * t237;
t172 = sin(t187);
t344 = g(3) * t172 + t173 * t189;
t343 = qJDD(2) + t342;
t192 = qJDD(1) * qJ(2);
t193 = qJD(1) * qJD(2);
t242 = g(1) * t207 + g(2) * t202;
t224 = -t242 + 0.2e1 * t193;
t339 = 0.2e1 * t192 + t224;
t25 = t237 ^ 2 - t86 ^ 2;
t21 = -t183 * t86 + t30;
t166 = g(3) * t173;
t329 = (qJD(5) * t58 + t24) * t204 + t199 * t23 - t61 * t282;
t146 = pkin(3) * t288 + qJD(1) * qJ(2);
t99 = -pkin(4) * t128 + t146;
t214 = -t172 * t342 - t99 * t86 + t166 - t329;
t35 = pkin(10) * t183 + t37;
t42 = -pkin(5) * t86 + pkin(10) * t237 + t99;
t12 = -t198 * t35 + t203 * t42;
t235 = t12 * t237 + t344 * t203 + t34 * t280;
t13 = t198 * t42 + t203 * t35;
t240 = -t13 * t237 + t348 * t198 + t34 * t279;
t212 = t237 * t99 - t272 - t330 + t344;
t22 = -t183 * t237 + t221;
t174 = pkin(4) * t199 + pkin(10);
t323 = pkin(4) * t127;
t48 = -pkin(10) * t86 - t323 - t350;
t336 = (qJD(6) * t174 + t48) * t345;
t176 = pkin(3) * t205 + pkin(4);
t297 = t200 * t204;
t291 = pkin(3) * t297 + t199 * t176;
t119 = pkin(10) + t291;
t180 = pkin(3) * t287;
t335 = (qJD(6) * t119 + t180 + t48) * t345;
t334 = (t345 * pkin(10) - t350) * t345;
t315 = pkin(8) - t208;
t141 = t315 * t201;
t142 = t315 * t206;
t292 = -t205 * t141 - t200 * t142;
t305 = pkin(1) * qJDD(1);
t331 = t305 - t343;
t218 = t333 * qJD(5) - t199 * t97 + t204 * t98;
t327 = -t182 * t236 - t183 * t218;
t133 = t315 * t286;
t134 = qJD(3) * t142;
t283 = qJD(4) * t205;
t226 = t200 * t133 - t205 * t134 + t141 * t284 - t142 * t283;
t49 = -pkin(9) * t98 + t226;
t220 = -t292 * qJD(4) + t205 * t133 + t134 * t200;
t50 = pkin(9) * t97 + t220;
t247 = t141 * t200 - t205 * t142;
t77 = -pkin(9) * t136 + t247;
t78 = -pkin(9) * t135 + t292;
t51 = t199 * t78 - t204 * t77;
t10 = -qJD(5) * t51 + t199 * t50 + t204 * t49;
t265 = pkin(10) * t182 + qJD(6) * t42 + t329;
t52 = t199 * t77 + t204 * t78;
t169 = t201 * pkin(3) + qJ(2);
t108 = pkin(4) * t135 + t169;
t53 = pkin(5) * t236 - pkin(10) * t333 + t108;
t326 = -t236 * t265 - (qJD(6) * t53 + t10) * t345 - t52 * t29 - t34 * t46 + t4 * t333;
t319 = t34 * t333;
t318 = t53 * t29;
t298 = t199 * t200;
t248 = -t117 * t200 - t104;
t62 = t248 - t322;
t293 = t205 * t117 - t103;
t63 = t121 + t293;
t312 = t199 * t62 + t204 * t63 - t176 * t281 - (-t200 * t282 + (t204 * t205 - t298) * qJD(4)) * pkin(3);
t311 = -t199 * t63 + t204 * t62 + t176 * t282 + (t200 * t281 + (t199 * t205 + t297) * qJD(4)) * pkin(3);
t306 = t136 * t190 - t97 * t191;
t210 = qJD(1) ^ 2;
t304 = qJ(2) * t210;
t302 = t127 * t128;
t300 = t198 * t202;
t299 = t198 * t207;
t296 = t202 * t203;
t295 = t203 * t207;
t196 = t206 ^ 2;
t290 = t201 ^ 2 - t196;
t209 = qJD(3) ^ 2;
t289 = -t209 - t210;
t158 = pkin(3) * t285 + qJD(2);
t277 = qJDD(3) * t201;
t270 = t333 * t280;
t107 = pkin(3) * t341 + t192 + t193;
t55 = pkin(4) * t70 + t107;
t6 = -pkin(5) * t221 - pkin(10) * t30 + t55;
t264 = qJD(6) * t35 - t6;
t246 = qJD(6) * t236 + qJD(1);
t38 = t199 * t60 + t307;
t243 = pkin(4) * t282 - t38;
t81 = pkin(4) * t98 + t158;
t241 = t29 * t333 - t345 * t46;
t239 = -t135 * t190 - t191 * t98;
t233 = -t265 + t166;
t232 = -pkin(3) * t298 + t176 * t204;
t229 = 0.2e1 * qJ(2) * t278 + qJDD(3) * t208;
t228 = -pkin(10) * t29 + t345 * t36 - t320;
t225 = -t304 + t342;
t223 = -t119 * t29 + t312 * t345 - t320;
t39 = t204 * t60 - t310;
t219 = -t174 * t29 - t320 + (-pkin(4) * t281 + t39) * t345;
t217 = -t208 * t209 + t339;
t185 = sin(t197);
t186 = cos(t197);
t213 = g(3) * t186 - t146 * t128 - t185 * t342 - t328;
t211 = g(3) * t185 + t146 * t127 + t186 * t342 + t222;
t184 = qJDD(3) * t206;
t175 = -pkin(4) * t204 - pkin(5);
t118 = -pkin(5) - t232;
t115 = t172 * t295 - t300;
t114 = t172 * t299 + t296;
t113 = t172 * t296 + t299;
t112 = -t172 * t300 + t295;
t101 = t180 - t323;
t74 = t127 ^ 2 - t128 ^ 2;
t57 = -t127 * t191 + (-t191 * t287 - t276) * t205 - t231;
t56 = -t128 * t191 + t69;
t14 = pkin(5) * t218 + pkin(10) * t46 + t81;
t11 = qJD(5) * t52 + t199 * t49 - t204 * t50;
t5 = t203 * t6;
t2 = [qJDD(1), -t342, t242, -0.2e1 * t305 + t343, t339, t331 * pkin(1) + (t224 + t192) * qJ(2), qJDD(1) * t196 - 0.2e1 * t201 * t267, -0.2e1 * t201 * t275 + 0.2e1 * t290 * t278, -t201 * t209 + t184, -t206 * t209 - t277, 0, t201 * t217 + t206 * t229, -t201 * t229 + t206 * t217, t127 * t97 + t136 * t69, t127 * t98 - t128 * t97 - t135 * t69 - t136 * t70, t306, t239, 0, t107 * t135 - t158 * t128 + t146 * t98 + t169 * t70 - t242 * t185 + t247 * t190 + t220 * t191, t107 * t136 - t158 * t127 - t146 * t97 + t169 * t69 - t186 * t242 - t190 * t292 - t191 * t226, t237 * t46 + t30 * t333, t218 * t237 + t221 * t333 - t236 * t30 - t46 * t86, t352, t327, 0, -t108 * t221 - t11 * t183 - t172 * t242 - t182 * t51 + t218 * t99 + t236 * t55 - t81 * t86, -t10 * t183 + t108 * t30 - t173 * t242 - t182 * t52 - t237 * t81 + t333 * t55 - t46 * t99, t203 * t353 - t73 * t270 -(-t198 * t73 - t203 * t71) * t46 + (-t16 - t19 * t203 + (t198 * t71 - t203 * t73) * qJD(6)) * t333, t18 * t236 + t203 * t241 + t218 * t73 - t270 * t345, -t19 * t236 - t198 * t241 - t218 * t71 - t333 * t76, t218 * t345 + t236 * t29, -g(1) * t115 - g(2) * t113 + t11 * t71 + t12 * t218 + t51 * t19 + t5 * t236 + (t14 * t345 + t318 + (-t236 * t35 - t345 * t52 + t319) * qJD(6)) * t203 + t326 * t198, g(1) * t114 - g(2) * t112 + t11 * t73 - t13 * t218 + t51 * t18 + (-(-qJD(6) * t52 + t14) * t345 - t318 + t264 * t236 - qJD(6) * t319) * t198 + t326 * t203; 0, 0, 0, qJDD(1), -t210, -t304 - t331, 0, 0, 0, 0, 0, t289 * t201 + t184, t289 * t206 - t277, 0, 0, 0, 0, 0, qJD(1) * t128 + t306, qJD(1) * t127 + t239, 0, 0, 0, 0, 0, qJD(1) * t86 + t352, qJD(1) * t237 + t327, 0, 0, 0, 0, 0, -t236 * t26 - t19 * t333 + t46 * t71 + (-t198 * t218 - t203 * t246) * t345, -t236 * t27 + (t198 * t246 - t203 * t218) * t345 - t353; 0, 0, 0, 0, 0, 0, t206 * t210 * t201, -t290 * t210, t275, -t276, qJDD(3), g(3) * t201 + t206 * t225 + t139, g(3) * t206 + (-t155 - t225) * t201, t302, t74, t56, t57, t190, -t248 * t191 + (t128 * t287 + t190 * t205 - t191 * t284) * pkin(3) + t211, t293 * t191 + (t127 * t287 - t190 * t200 - t191 * t283) * pkin(3) + t213, t316, t25, t21, t22, t182, t101 * t86 + t182 * t232 - t183 * t311 + t212, t101 * t237 - t182 * t291 + t183 * t312 + t214, t9, t1, t8, t7, t317, t118 * t19 + t311 * t71 + (-t348 - t335) * t203 + t223 * t198 + t235, t118 * t18 + t311 * t73 + t223 * t203 + (-t344 + t335) * t198 + t240; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t302, t74, t56, t57, t190, -t191 * t238 + t211, t191 * t249 + t213, t316, t25, t21, t22, t182, t183 * t38 + (-t127 * t86 + t182 * t204 - t183 * t282) * pkin(4) + t212, t183 * t39 + (-t127 * t237 - t182 * t199 - t183 * t281) * pkin(4) + t214, t9, t1, t8, t7, t317, t175 * t19 + t243 * t71 + (-t348 - t336) * t203 + t219 * t198 + t235, t175 * t18 + t243 * t73 + t219 * t203 + (-t344 + t336) * t198 + t240; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t316, t25, t21, t22, t182, t183 * t37 + t212, t183 * t36 + t214, t9, t1, t8, t7, t317, -pkin(5) * t19 - t37 * t71 + t228 * t198 + (-t348 - t334) * t203 + t235, -pkin(5) * t18 - t37 * t73 + t228 * t203 + (-t344 + t334) * t198 + t240; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73 * t71, -t71 ^ 2 + t73 ^ 2, t345 * t71 + t18, t345 * t73 - t19, t29, -g(1) * t112 - g(2) * t114 + t13 * t345 + t198 * t233 - t279 * t35 - t34 * t73 + t5, g(1) * t113 - g(2) * t115 + t12 * t345 + t198 * t264 + t203 * t233 + t34 * t71;];
tau_reg  = t2;
