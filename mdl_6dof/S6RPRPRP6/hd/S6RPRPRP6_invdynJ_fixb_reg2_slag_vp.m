% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPRPRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPRP6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP6_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:19:52
% EndTime: 2019-03-09 03:20:00
% DurationCPUTime: 5.19s
% Computational Cost: add. (7174->505), mult. (16650->583), div. (0->0), fcn. (12393->10), ass. (0->264)
t187 = cos(pkin(9));
t327 = cos(qJ(3));
t259 = t327 * t187;
t241 = qJD(1) * t259;
t186 = sin(pkin(9));
t191 = sin(qJ(3));
t284 = t186 * t191;
t258 = qJD(1) * t284;
t133 = -t241 + t258;
t190 = sin(qJ(5));
t193 = cos(qJ(5));
t105 = qJD(3) * t190 - t193 * t133;
t142 = t327 * t186 + t191 * t187;
t341 = t142 * qJD(1);
t348 = qJD(5) + t341;
t353 = t105 * t348;
t138 = t142 * qJD(3);
t254 = qJDD(1) * t327;
t264 = t186 * qJDD(1);
t232 = -t187 * t254 + t191 * t264;
t199 = qJD(1) * t138 + t232;
t269 = qJD(5) * t193;
t270 = qJD(5) * t190;
t59 = qJD(3) * t270 - t193 * qJDD(3) - t133 * t269 - t190 * t199;
t358 = -t59 + t353;
t357 = t59 + t353;
t107 = qJD(3) * t193 + t133 * t190;
t352 = t107 * t348;
t60 = qJD(5) * t107 + qJDD(3) * t190 - t193 * t199;
t42 = -t60 + t352;
t172 = t187 * pkin(2) + pkin(1);
t146 = -qJD(1) * t172 + qJD(2);
t209 = -qJ(4) * t341 + t146;
t329 = pkin(3) + pkin(8);
t65 = t329 * t133 + t209;
t316 = pkin(7) + qJ(2);
t149 = t316 * t186;
t143 = qJD(1) * t149;
t150 = t316 * t187;
t144 = qJD(1) * t150;
t100 = t327 * t143 + t144 * t191;
t340 = qJD(4) + t100;
t277 = pkin(4) * t341 + t340;
t71 = -t329 * qJD(3) + t277;
t32 = -t190 * t65 + t193 * t71;
t356 = t32 * t348;
t33 = t190 * t71 + t193 * t65;
t355 = t33 * t348;
t20 = -qJ(6) * t105 + t33;
t354 = t348 * t20;
t182 = pkin(9) + qJ(3);
t177 = cos(t182);
t192 = sin(qJ(1));
t194 = cos(qJ(1));
t237 = g(1) * t194 + g(2) * t192;
t176 = sin(t182);
t319 = g(3) * t176;
t210 = -t237 * t177 - t319;
t170 = g(3) * t177;
t289 = t176 * t194;
t290 = t176 * t192;
t260 = -g(1) * t289 - g(2) * t290 + t170;
t19 = -qJ(6) * t107 + t32;
t16 = pkin(5) * t348 + t19;
t304 = qJDD(1) * pkin(1);
t175 = qJDD(2) - t304;
t342 = t341 * qJD(3);
t200 = t232 + t342;
t263 = t187 * qJDD(1);
t271 = qJD(3) * t191;
t257 = t186 * t271;
t261 = qJD(3) * t241 + t186 * t254 + t191 * t263;
t99 = qJD(1) * t257 - t261;
t46 = -pkin(2) * t263 + pkin(3) * t200 + t99 * qJ(4) - qJD(4) * t341 + t175;
t29 = pkin(8) * t199 + t46;
t265 = qJD(1) * qJD(2);
t332 = qJDD(1) * t316 + t265;
t114 = t332 * t186;
t115 = t332 * t187;
t256 = qJD(3) * t327;
t244 = t327 * t114 + t191 * t115 - t143 * t271 + t144 * t256;
t229 = qJDD(4) + t244;
t36 = -pkin(4) * t99 - t329 * qJDD(3) + t229;
t3 = t190 * t36 + t193 * t29 + t71 * t269 - t65 * t270;
t313 = qJ(6) * t60;
t2 = -qJD(6) * t105 + t3 - t313;
t336 = t16 * t348 - t2;
t314 = qJ(6) * t59;
t97 = -qJDD(5) + t99;
t328 = pkin(5) * t97;
t4 = -qJD(5) * t33 - t190 * t29 + t193 * t36;
t1 = -qJD(6) * t107 + t314 - t328 + t4;
t337 = t1 + t354;
t351 = -t190 * t336 + t193 * t337 + t260;
t250 = t190 * t348;
t218 = -t193 * t97 - t250 * t348;
t345 = g(1) * t192 - g(2) * t194;
t222 = -t175 + t345;
t245 = t191 * t114 - t327 * t115 + t143 * t256 + t144 * t271;
t350 = -qJD(3) * t100 - t210 + t245;
t349 = -t232 - 0.2e1 * t342;
t249 = t193 * t348;
t275 = t177 * pkin(3) + t176 * qJ(4);
t131 = t341 ^ 2;
t330 = t133 ^ 2;
t344 = -t330 - t131;
t343 = -t330 + t131;
t339 = qJ(2) * qJDD(1);
t279 = t193 * t194;
t282 = t190 * t192;
t125 = t176 * t279 - t282;
t280 = t192 * t193;
t281 = t190 * t194;
t127 = t176 * t280 + t281;
t338 = -g(1) * t125 - g(2) * t127 + t193 * t170;
t335 = t4 + t355;
t185 = qJD(3) * qJ(4);
t101 = -t191 * t143 + t327 * t144;
t82 = -t133 * pkin(4) + t101;
t74 = t185 + t82;
t334 = t329 * t97 + t348 * t74;
t103 = -t191 * t149 + t327 * t150;
t85 = (qJD(2) * t186 + qJD(3) * t150) * t191 - qJD(2) * t259 + t149 * t256;
t333 = qJD(3) * t85 - qJDD(3) * t103 - t176 * t345;
t331 = t107 ^ 2;
t326 = pkin(5) * t190;
t325 = pkin(8) * t177;
t317 = pkin(4) + t316;
t315 = -t19 + t16;
t306 = qJ(4) * t133;
t77 = t329 * t341 + t306;
t44 = t190 * t82 + t193 * t77;
t141 = -t259 + t284;
t223 = -qJ(4) * t142 - t172;
t80 = t329 * t141 + t223;
t102 = t327 * t149 + t150 * t191;
t89 = pkin(4) * t142 + t102;
t48 = t190 * t89 + t193 * t80;
t312 = t190 * t60;
t311 = t193 * t59;
t268 = qJD(6) * t193;
t278 = qJ(6) + t329;
t305 = qJ(6) * t341;
t76 = t193 * t82;
t309 = t270 * t278 - t268 + pkin(5) * t133 - t76 - (-t77 - t305) * t190;
t148 = t278 * t193;
t308 = -qJD(5) * t148 - qJD(6) * t190 - t193 * t305 - t44;
t174 = pkin(5) * t193 + pkin(4);
t307 = pkin(5) * t269 + t174 * t341 + t340;
t303 = qJDD(3) * pkin(3);
t301 = t105 * t133;
t300 = t107 * t105;
t298 = t107 * t133;
t297 = t348 * t133;
t296 = t133 * t341;
t295 = t133 * t138;
t294 = t138 * t190;
t293 = t138 * t193;
t292 = t141 * t190;
t291 = t141 * t193;
t188 = -qJ(6) - pkin(8);
t288 = t177 * t188;
t287 = t177 * t190;
t286 = t177 * t192;
t285 = t177 * t194;
t283 = t316 * t194;
t274 = t174 + t316;
t180 = t186 ^ 2;
t181 = t187 ^ 2;
t273 = t180 + t181;
t272 = qJD(3) * t101;
t262 = pkin(5) * t287;
t253 = t105 * pkin(5) + qJD(6);
t252 = -qJ(6) * t141 - t80;
t251 = t273 * qJD(1) ^ 2;
t243 = 0.2e1 * t273;
t154 = t194 * t172;
t242 = g(2) * (pkin(3) * t285 + qJ(4) * t289 + t154);
t240 = -g(1) * t286 + g(2) * t285;
t239 = g(1) * t127 - g(2) * t125;
t126 = t176 * t281 + t280;
t128 = -t176 * t282 + t279;
t238 = -g(1) * t128 - g(2) * t126;
t233 = t3 - t356;
t231 = -t190 * t32 + t193 * t33;
t137 = -t187 * t256 + t257;
t230 = -t137 * t341 - t142 * t99;
t228 = qJ(4) * t137 - qJD(4) * t142;
t225 = qJD(3) * t137 - qJDD(3) * t142;
t224 = qJD(3) * t138 + qJDD(3) * t141;
t221 = t176 * t326 - t288;
t220 = -t172 - t275;
t183 = qJDD(3) * qJ(4);
t184 = qJD(3) * qJD(4);
t49 = -t183 - t184 + t245;
t61 = t329 * t138 + t228;
t86 = qJD(2) * t142 + qJD(3) * t103;
t67 = -t137 * pkin(4) + t86;
t13 = t190 * t67 + t193 * t61 + t89 * t269 - t270 * t80;
t217 = t141 * t269 + t294;
t216 = t141 * t270 - t293;
t145 = -qJDD(1) * t172 + qJDD(2);
t215 = -t244 - t260;
t214 = t137 * t133 - t138 * t341 + t99 * t141;
t213 = t190 * t97 - t249 * t348;
t211 = t222 + t304;
t207 = -t86 * qJD(3) - t102 * qJDD(3) - t240;
t39 = -pkin(4) * t199 - t49;
t15 = t60 * pkin(5) + qJDD(6) + t39;
t206 = t15 + t210;
t84 = pkin(3) * t133 + t209;
t205 = t341 * t84 + qJDD(4) - t215;
t204 = g(1) * t126 - g(2) * t128 - g(3) * t287 - t3;
t203 = t243 * t265 - t237;
t66 = -pkin(4) * t138 - t85;
t202 = -t102 * t99 + t85 * t133 + t341 * t86 - t237;
t201 = qJD(5) * t329 * t348 + t210 + t39;
t198 = t4 + t338;
t173 = qJ(4) + t326;
t153 = qJ(4) * t285;
t151 = qJ(4) * t286;
t147 = t278 * t190;
t119 = qJD(3) * t133;
t104 = t105 ^ 2;
t98 = pkin(3) * t141 + t223;
t96 = pkin(3) * t341 + t306;
t95 = -t185 - t101;
t94 = -qJD(3) * pkin(3) + t340;
t90 = -t141 * pkin(4) + t103;
t88 = t193 * t89;
t79 = pkin(3) * t138 + t228;
t78 = -t119 + t99;
t70 = -t141 * t174 + t103;
t64 = t193 * t67;
t56 = t193 * t60;
t55 = -t104 + t331;
t54 = t253 + t74;
t51 = -t137 * t348 - t142 * t97;
t50 = t229 - t303;
t47 = -t190 * t80 + t88;
t45 = pkin(5) * t216 + t66;
t43 = -t190 * t77 + t76;
t40 = qJ(6) * t291 + t48;
t38 = -qJD(3) * t107 + t213;
t37 = -qJD(3) * t105 + t218;
t30 = pkin(5) * t142 + t190 * t252 + t88;
t28 = t213 - t301;
t27 = t213 + t301;
t26 = t298 - t218;
t25 = t218 + t298;
t23 = t193 * t353 + t312;
t22 = -t107 * t250 - t311;
t18 = t105 * t216 - t291 * t60;
t17 = t107 * t217 - t292 * t59;
t14 = -qJD(5) * t48 - t190 * t61 + t64;
t12 = -t107 * t137 - t142 * t59 + t217 * t348 - t292 * t97;
t11 = t105 * t137 - t142 * t60 - t216 * t348 - t291 * t97;
t10 = t190 * t42 - t193 * t358;
t9 = -t107 * t249 + t190 * t357 - t56;
t8 = t190 * t358 + t193 * t352 - t56;
t7 = -qJ(6) * t216 + t141 * t268 + t13;
t6 = -pkin(5) * t137 + t64 + t252 * t269 + (-qJ(6) * t138 - qJD(5) * t89 - qJD(6) * t141 - t61) * t190;
t5 = (-t105 * t190 + t107 * t193) * t138 + (-t312 - t311 + (-t105 * t193 - t107 * t190) * qJD(5)) * t141;
t21 = [0, 0, 0, 0, 0, qJDD(1), t345, t237, 0, 0, t180 * qJDD(1), 0.2e1 * t186 * t263, 0, t181 * qJDD(1), 0, 0, t211 * t187, -t211 * t186, t243 * t339 + t203, pkin(1) * t222 + (t273 * t339 + t203) * qJ(2), t230, -t142 * t200 + t214, -t225, t141 * t200 + t295, -t224, 0, t146 * t138 + t145 * t141 - t172 * t199 + t207, -t137 * t146 + t142 * t145 + t172 * t99 + t333, -t100 * t137 - t101 * t138 - t103 * t200 + t141 * t245 + t142 * t244 + t202, -t245 * t103 - t101 * t85 + t244 * t102 + t100 * t86 - t145 * t172 - g(1) * (-t172 * t192 + t283) - g(2) * (t192 * t316 + t154) 0, t225, t224, t230, -t142 * t199 + t214, t141 * t199 + t295, -t103 * t199 - t94 * t137 + t95 * t138 + t49 * t141 + t50 * t142 + t202, -t79 * t133 - t84 * t138 - t46 * t141 - t200 * t98 - t207, t137 * t84 - t142 * t46 - t341 * t79 + t98 * t99 - t333, t46 * t98 + t84 * t79 - t49 * t103 + t95 * t85 + t50 * t102 + t94 * t86 - g(1) * t283 - t242 + (-g(1) * t220 - g(2) * t316) * t192, t17, t5, t12, t18, t11, t51, -t74 * t293 + t105 * t66 + t348 * t14 - t137 * t32 + t142 * t4 - t47 * t97 + t60 * t90 + (-t193 * t39 + t270 * t74) * t141 + t238, t74 * t294 + t107 * t66 - t348 * t13 + t137 * t33 - t142 * t3 + t48 * t97 - t59 * t90 + (t190 * t39 + t269 * t74) * t141 + t239, -t105 * t13 - t107 * t14 + t47 * t59 - t48 * t60 + t231 * t138 + (-t190 * t4 + t193 * t3 + (-t190 * t33 - t193 * t32) * qJD(5)) * t141 - t240, t3 * t48 + t33 * t13 + t4 * t47 + t32 * t14 + t39 * t90 + t74 * t66 - t242 + (-g(1) * t317 - g(2) * t325) * t194 + (-g(1) * (t220 - t325) - g(2) * t317) * t192, t17, t5, t12, t18, t11, t51, -t54 * t293 + t1 * t142 + t105 * t45 + t348 * t6 - t137 * t16 - t30 * t97 + t60 * t70 + (-t15 * t193 + t270 * t54) * t141 + t238, t54 * t294 + t107 * t45 - t348 * t7 + t137 * t20 - t142 * t2 + t40 * t97 - t59 * t70 + (t15 * t190 + t269 * t54) * t141 + t239, -t105 * t7 - t107 * t6 + t30 * t59 - t40 * t60 + (-t16 * t190 + t193 * t20) * t138 + (-t1 * t190 + t193 * t2 + (-t16 * t193 - t190 * t20) * qJD(5)) * t141 - t240, t2 * t40 + t20 * t7 + t1 * t30 + t16 * t6 + t15 * t70 + t54 * t45 - t242 + (-g(1) * t274 - g(2) * t221) * t194 + (-g(1) * (t220 - t221) - g(2) * t274) * t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t263, t264, -t251, -qJ(2) * t251 - t222, 0, 0, 0, 0, 0, 0, -t349 (-t133 - t258) * qJD(3) + t261, t344, -t100 * t341 + t101 * t133 + t145 - t345, 0, 0, 0, 0, 0, 0, t344, t349, t119 + t99, -t95 * t133 - t341 * t94 - t345 + t46, 0, 0, 0, 0, 0, 0, t27, t26, t8, t133 * t74 - t190 * t335 + t233 * t193 - t345, 0, 0, 0, 0, 0, 0, t27, t26, t8, t133 * t54 - t190 * t337 - t193 * t336 - t345; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t296, t343 (t133 - t258) * qJD(3) + t261, -t296, -t232, qJDD(3), -t146 * t341 + t215 + t272, t133 * t146 + t350, 0, 0, qJDD(3), t78, t232, t296, t343, -t296, pkin(3) * t99 - qJ(4) * t199 - (t101 + t95) * t341 + (t94 - t340) * t133, t133 * t96 + t205 - t272 - 0.2e1 * t303, -t133 * t84 + t341 * t96 + 0.2e1 * t183 + 0.2e1 * t184 - t350, -t49 * qJ(4) - t50 * pkin(3) - t84 * t96 - t94 * t101 - g(1) * (-pkin(3) * t289 + t153) - g(2) * (-pkin(3) * t290 + t151) - g(3) * t275 - t340 * t95, t22, t9, t25, t23, t28, t297, qJ(4) * t60 + t277 * t105 + t133 * t32 + t201 * t190 + t193 * t334 - t348 * t43, -qJ(4) * t59 + t277 * t107 - t133 * t33 - t190 * t334 + t201 * t193 + t348 * t44, t105 * t44 + t107 * t43 + (-t341 * t33 - t329 * t59 - t4 + (t105 * t329 - t33) * qJD(5)) * t193 + (t341 * t32 + t329 * t60 - t3 + (-t107 * t329 + t32) * qJD(5)) * t190 - t260, t39 * qJ(4) - t33 * t44 - t32 * t43 - g(1) * t153 - g(2) * t151 - g(3) * (t275 + t325) + t277 * t74 + (-qJD(5) * t231 + t176 * t237 - t3 * t190 - t4 * t193) * t329, t22, t9, t25, t23, t28, t297, t105 * t307 + t133 * t16 + t148 * t97 + t173 * t60 + t190 * t206 + t249 * t54 + t309 * t348, t107 * t307 - t133 * t20 - t147 * t97 - t173 * t59 + t193 * t206 - t250 * t54 - t308 * t348, -t308 * t105 - t309 * t107 + t147 * t60 - t148 * t59 - t351, -t2 * t147 - t1 * t148 + t15 * t173 - g(1) * (t194 * t262 + t153) - g(2) * (t192 * t262 + t151) - g(3) * (t275 - t288) + t307 * t54 + t308 * t20 + (-g(3) * t326 + t237 * (pkin(3) - t188)) * t176 + t309 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, qJDD(3) - t296, -qJD(3) ^ 2 - t131, qJD(3) * t95 + t205 - t303, 0, 0, 0, 0, 0, 0, t37, t38, t10, -qJD(3) * t74 + t233 * t190 + t193 * t335 + t260, 0, 0, 0, 0, 0, 0, t37, t38, t10, -qJD(3) * t54 + t351; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t300, t55, t358, -t300, t42, -t97, -t107 * t74 + t198 + t355, t105 * t74 + t204 + t356, 0, 0, t300, t55, t358, -t300, t42, -t97, -0.2e1 * t328 + t314 + t354 + (-t253 - t54) * t107 + t198, -pkin(5) * t331 + t313 + t348 * t19 + (qJD(6) + t54) * t105 + t204, pkin(5) * t59 - t105 * t315, t315 * t20 + (-t54 * t107 + t1 + t338) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60 + t352, -t357, -t104 - t331, -g(1) * t285 - g(2) * t286 + t20 * t105 + t16 * t107 + t15 - t319;];
tau_reg  = t21;
