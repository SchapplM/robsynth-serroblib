% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPPPR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:19:48
% EndTime: 2019-03-09 08:19:58
% DurationCPUTime: 4.88s
% Computational Cost: add. (3096->494), mult. (6612->602), div. (0->0), fcn. (4187->8), ass. (0->254)
t188 = sin(pkin(9));
t195 = cos(qJ(2));
t294 = qJD(1) * t195;
t271 = t188 * t294;
t189 = cos(pkin(9));
t292 = qJD(2) * t189;
t109 = -t271 + t292;
t284 = t195 * qJDD(1);
t161 = pkin(7) * t284;
t183 = qJDD(2) * qJ(3);
t184 = qJD(2) * qJD(3);
t192 = sin(qJ(2));
t286 = qJD(1) * qJD(2);
t269 = t192 * t286;
t81 = pkin(7) * t269 - t161 - t183 - t184;
t229 = pkin(3) * t284 + qJDD(4) - t81;
t204 = pkin(3) * t269 - t229;
t260 = t189 * qJDD(2) - t188 * t284;
t74 = t188 * t269 + t260;
t323 = t74 * qJ(5);
t200 = t109 * qJD(5) + t204 + t323;
t73 = qJDD(2) * t188 + (-t269 + t284) * t189;
t336 = pkin(4) * t73;
t11 = -t200 + t336;
t193 = sin(qJ(1));
t196 = cos(qJ(1));
t252 = g(1) * t196 + g(2) * t193;
t233 = t252 * t195;
t331 = g(3) * t192;
t347 = t233 + t331;
t348 = t11 - t347;
t338 = pkin(3) + pkin(7);
t268 = t195 * t286;
t285 = t192 * qJDD(1);
t220 = t268 + t285;
t273 = t189 * t294;
t293 = qJD(2) * t188;
t107 = t273 + t293;
t191 = sin(qJ(6));
t194 = cos(qJ(6));
t54 = -t107 * t194 + t109 * t191;
t295 = qJD(1) * t192;
t163 = pkin(7) * t295;
t287 = pkin(3) * t295 + qJD(3) + t163;
t106 = t109 ^ 2;
t346 = -t107 ^ 2 - t106;
t345 = qJD(5) * t189 - t287;
t304 = t194 * t188;
t310 = t189 * t195;
t344 = t191 * t310 - t195 * t304;
t147 = -qJD(6) + t295;
t343 = -qJD(6) - t147;
t56 = t107 * t191 + t109 * t194;
t13 = qJD(6) * t56 + t191 * t74 - t194 * t73;
t111 = -qJDD(6) + t220;
t238 = -t189 * t191 + t304;
t237 = t188 * t191 + t189 * t194;
t222 = t192 * t237;
t95 = t237 * qJD(6);
t327 = -qJD(1) * t222 + t95;
t342 = t111 * t238 - t147 * t327;
t274 = t189 * t295;
t275 = t188 * t295;
t328 = qJD(6) * t238 + t191 * t274 - t194 * t275;
t341 = t111 * t237 + t147 * t328;
t318 = qJ(5) * t188;
t337 = pkin(4) + pkin(5);
t232 = t189 * t337 + t318;
t306 = t193 * t195;
t140 = qJ(3) * t306;
t303 = t195 * t196;
t143 = qJ(3) * t303;
t172 = t192 * qJ(3);
t178 = t195 * pkin(2);
t298 = t178 + t172;
t278 = qJ(4) * t195 + t298;
t330 = pkin(2) + qJ(4);
t340 = t192 * t252 * t330 - g(1) * t143 - g(2) * t140 - g(3) * t278;
t267 = t188 * t285;
t186 = t192 ^ 2;
t198 = qJD(1) ^ 2;
t313 = t186 * t198;
t339 = qJD(2) * (t109 + t271) + t189 * t313 + t267;
t335 = g(1) * t193;
t332 = g(2) * t196;
t182 = g(3) * t195;
t329 = pkin(8) - t330;
t146 = pkin(2) * t269;
t319 = qJ(3) * t195;
t241 = qJ(4) * t192 - t319;
t288 = t192 * qJD(3);
t201 = qJD(2) * t241 - t195 * qJD(4) - t288;
t264 = -pkin(1) - t172;
t217 = -t195 * t330 + t264;
t23 = qJD(1) * t201 + qJDD(1) * t217 + t146;
t145 = pkin(7) * t268;
t160 = pkin(7) * t285;
t265 = qJDD(3) + t145 + t160;
t47 = pkin(3) * t220 - qJD(2) * qJD(4) - qJDD(2) * t330 + t265;
t10 = t188 * t47 + t189 * t23;
t79 = t217 * qJD(1);
t83 = -qJD(2) * t330 + t287;
t29 = t188 * t83 + t189 * t79;
t326 = t147 * t54;
t65 = t188 * t73;
t325 = t189 * t47;
t324 = t56 * t147;
t322 = -t232 * t295 + t345;
t247 = pkin(4) * t189 + t318;
t321 = t247 * t295 - t345;
t290 = qJD(2) * t195;
t119 = t338 * t290;
t291 = qJD(2) * t192;
t167 = pkin(2) * t291;
t64 = t167 + t201;
t32 = t119 * t188 + t189 * t64;
t165 = pkin(7) * t294;
t118 = pkin(3) * t294 + t165;
t168 = pkin(2) * t295;
t89 = qJD(1) * t241 + t168;
t46 = t118 * t188 + t189 * t89;
t104 = -pkin(1) - t278;
t130 = t338 * t192;
t58 = t104 * t189 + t130 * t188;
t320 = pkin(7) * qJDD(2);
t317 = qJDD(2) * pkin(2);
t312 = t188 * t195;
t311 = t189 * qJ(5);
t309 = t192 * t193;
t308 = t192 * t196;
t307 = t192 * t198;
t185 = qJD(2) * qJ(3);
t93 = qJD(4) + t185 + t118;
t219 = t109 * qJ(5) - t93;
t33 = pkin(4) * t107 - t219;
t302 = -qJD(4) + t33;
t301 = qJD(4) - t93;
t266 = t189 * t285;
t300 = (-t189 * t268 - t266) * t330;
t131 = t338 * t195;
t179 = t196 * pkin(7);
t297 = pkin(3) * t196 + t179;
t187 = t195 ^ 2;
t296 = t186 - t187;
t289 = qJD(4) * t107;
t170 = t192 * qJD(5);
t25 = qJ(5) * t295 + t29;
t35 = qJ(5) * t294 + t46;
t50 = qJ(5) * t192 + t58;
t283 = t337 * t192;
t281 = t195 * t307;
t279 = t25 * t295;
t277 = -g(1) * t308 - g(2) * t309 + t182;
t19 = t188 * t23;
t207 = -pkin(4) * t220 + qJDD(5) + t19;
t8 = t207 - t325;
t4 = -pkin(5) * t220 - t74 * pkin(8) + t8;
t7 = qJ(5) * t220 + qJD(1) * t170 + t10;
t5 = pkin(8) * t73 + t7;
t276 = -t191 * t5 + t194 * t4;
t272 = t330 * t290;
t270 = qJD(5) * t312;
t28 = -t188 * t79 + t189 * t83;
t45 = t189 * t118 - t188 * t89;
t31 = t189 * t119 - t188 * t64;
t57 = -t104 * t188 + t189 * t130;
t262 = -qJD(2) * pkin(2) + qJD(3);
t261 = t330 * t65 - t277;
t21 = qJ(5) * t290 + t170 + t32;
t259 = pkin(1) * t196 + pkin(2) * t303 + pkin(7) * t193 + qJ(3) * t308;
t258 = -t160 - t277;
t257 = t330 * t267;
t97 = t188 * t193 - t189 * t308;
t99 = t188 * t196 + t189 * t309;
t256 = -g(1) * t99 - g(2) * t97;
t100 = -t188 * t309 + t189 * t196;
t98 = t188 * t308 + t189 * t193;
t255 = -g(1) * t100 - g(2) * t98;
t254 = qJD(5) - t28;
t197 = qJD(2) ^ 2;
t253 = pkin(7) * t197 + t332;
t251 = t7 * t188 - t8 * t189;
t250 = t191 * t4 + t194 * t5;
t9 = -t19 + t325;
t248 = t10 * t188 + t9 * t189;
t246 = pkin(4) * t188 - t311;
t14 = -t109 * pkin(8) - qJD(1) * t283 + t254;
t15 = pkin(8) * t107 + t25;
t2 = t14 * t194 - t15 * t191;
t3 = t14 * t191 + t15 * t194;
t30 = pkin(8) * t312 - t283 - t57;
t36 = pkin(8) * t310 + t50;
t245 = -t191 * t36 + t194 * t30;
t244 = t191 * t30 + t194 * t36;
t243 = t100 * t194 + t191 * t99;
t242 = t100 * t191 - t194 * t99;
t240 = qJD(4) * t109 + t330 * t74;
t124 = t163 + t262;
t129 = -t165 - t185;
t239 = t124 * t195 + t129 * t192;
t236 = t147 ^ 2;
t235 = t264 - t178;
t234 = pkin(3) * t193 + qJ(4) * t303 + t259;
t103 = t235 * qJD(1);
t231 = t103 * t295 + qJDD(3) - t258;
t230 = -0.2e1 * pkin(1) * t286 - t320;
t121 = t329 * t188;
t218 = -pkin(8) * t188 * t192 - t195 * t337;
t228 = -qJD(1) * t218 + qJD(4) * t189 - qJD(6) * t121 + t45;
t122 = t329 * t189;
t227 = pkin(8) * t274 - qJD(4) * t188 - qJD(6) * t122 - t35;
t226 = qJD(6) * t54 - t191 * t73 - t194 * t74;
t24 = -pkin(4) * t295 + t254;
t225 = (t188 * t24 + t189 * t25) * t192;
t224 = (-t188 * t28 + t189 * t29) * t192;
t223 = -qJ(3) * t290 - t288;
t221 = t238 * t192;
t216 = 0.2e1 * qJDD(1) * pkin(1) - t253;
t126 = -pkin(1) - t298;
t214 = t320 + (-qJD(1) * t126 - t103) * qJD(2);
t213 = -t107 * t274 + t109 * t275 - t189 * t74 - t65;
t210 = t217 * t335;
t208 = -t204 - t347;
t51 = qJD(1) * t223 + qJDD(1) * t235 + t146;
t91 = t167 + t223;
t205 = qJD(1) * t91 + qJDD(1) * t126 + t253 + t51;
t203 = t109 * t295 + t73;
t88 = t265 - t317;
t202 = qJD(2) * t239 + t88 * t192 - t81 * t195;
t199 = (-pkin(3) * t286 - g(3)) * t192 - t233 + t229;
t150 = g(1) * t306;
t123 = qJ(3) + t246;
t117 = t338 * t291;
t115 = -qJ(3) * t294 + t168;
t96 = -t188 * t337 - qJ(3) + t311;
t84 = t237 * t195;
t71 = t195 * t247 + t131;
t61 = -t195 * t232 - t131;
t52 = -pkin(4) * t192 - t57;
t49 = t270 + (-t247 - t338) * t291;
t44 = t191 * t97 + t194 * t98;
t43 = -t191 * t98 + t194 * t97;
t42 = t266 - t188 * t313 + (-t107 + t273) * qJD(2);
t41 = (-t107 + t293) * t295 + t260;
t39 = qJD(2) * t221 + t195 * t95;
t38 = qJD(2) * t222 + qJD(6) * t344;
t37 = -pkin(4) * t294 - t45;
t34 = -t270 + (t232 + t338) * t291;
t26 = -pkin(4) * t290 - t31;
t18 = -t107 * t337 + t219;
t17 = -pkin(8) * t189 * t291 + t21;
t16 = qJD(2) * t218 - t31;
t6 = -t337 * t73 + t200;
t1 = [qJDD(1), -t332 + t335, t252, qJDD(1) * t186 + 0.2e1 * t192 * t268, 0.2e1 * t192 * t284 - 0.2e1 * t286 * t296, qJDD(2) * t192 + t195 * t197, qJDD(2) * t195 - t192 * t197, 0, t192 * t230 + t195 * t216 + t150, t230 * t195 + (-t216 - t335) * t192 (t186 + t187) * qJDD(1) * pkin(7) + t202 - t252, t192 * t214 + t195 * t205 - t150, t214 * t195 + (-t205 + t335) * t192, pkin(7) * t202 - g(1) * t179 - g(2) * t259 + t103 * t91 + t51 * t126 - t235 * t335, -t117 * t107 + t131 * t73 + (-t204 * t189 + (qJD(1) * t57 + t28) * qJD(2)) * t195 + (qJD(1) * t31 + qJDD(1) * t57 - t292 * t93 + t9) * t192 + t255, -t117 * t109 + t131 * t74 + (t204 * t188 + (-qJD(1) * t58 - t29) * qJD(2)) * t195 + (-qJD(1) * t32 - qJDD(1) * t58 + t293 * t93 - t10) * t192 - t256, -t32 * t107 - t31 * t109 - t57 * t74 - t58 * t73 + t150 + qJD(2) * t224 + (-t10 * t189 + t188 * t9 - t332) * t195, -g(1) * t297 - g(2) * t234 + t10 * t58 - t93 * t117 - t131 * t204 + t28 * t31 + t29 * t32 + t9 * t57 - t210, t49 * t107 + t71 * t73 + (t11 * t189 + (-qJD(1) * t52 - t24) * qJD(2)) * t195 + (-qJD(1) * t26 - qJDD(1) * t52 - t292 * t33 - t8) * t192 + t255, -t21 * t107 + t26 * t109 - t50 * t73 + t52 * t74 + t150 + qJD(2) * t225 + (-t188 * t8 - t189 * t7 - t332) * t195, -t49 * t109 - t71 * t74 + (t11 * t188 + (qJD(1) * t50 + t25) * qJD(2)) * t195 + (qJD(1) * t21 + qJDD(1) * t50 - t293 * t33 + t7) * t192 + t256, t7 * t50 + t25 * t21 + t11 * t71 + t33 * t49 + t8 * t52 + t24 * t26 - g(1) * (t100 * pkin(4) + t99 * qJ(5) + t297) - g(2) * (pkin(4) * t98 + qJ(5) * t97 + t234) - t210, -t226 * t344 + t39 * t56, -t13 * t344 - t226 * t84 - t38 * t56 - t39 * t54, -t111 * t344 - t147 * t39 + t192 * t226 - t290 * t56, -t111 * t84 + t13 * t192 + t147 * t38 + t290 * t54, t111 * t192 + t147 * t290 -(t194 * t16 - t191 * t17) * t147 - t245 * t111 - t276 * t192 - t2 * t290 + t34 * t54 + t61 * t13 - t6 * t84 + t18 * t38 - g(1) * t243 - g(2) * t44 + (t147 * t244 + t192 * t3) * qJD(6) (t191 * t16 + t17 * t194) * t147 + t244 * t111 + t250 * t192 + t3 * t290 + t34 * t56 - t61 * t226 + t6 * t344 + t18 * t39 + g(1) * t242 - g(2) * t43 + (t147 * t245 + t192 * t2) * qJD(6); 0, 0, 0, -t281, t296 * t198, t285, t284, qJDD(2), pkin(1) * t307 + t258, t331 - t161 + (pkin(1) * t198 + t252) * t195 (-pkin(2) * t192 + t319) * qJDD(1) + ((-t129 - t185) * t192 + (-t124 + t262) * t195) * qJD(1), -t115 * t294 + t231 - 0.2e1 * t317, t161 + 0.2e1 * t183 + 0.2e1 * t184 + (qJD(1) * t115 - g(3)) * t192 + (qJD(1) * t103 - t252) * t195, -t81 * qJ(3) - t129 * qJD(3) - t88 * pkin(2) - t103 * t115 - g(1) * (-pkin(2) * t308 + t143) - g(2) * (-pkin(2) * t309 + t140) - g(3) * t298 - t239 * qJD(1) * pkin(7), qJ(3) * t73 + t287 * t107 + t208 * t188 + (-t195 * t28 + (-t189 * t301 - t45) * t192) * qJD(1) + t300, t257 + qJ(3) * t74 + t287 * t109 + t208 * t189 + (t192 * t46 + t195 * t29 + (t192 * t301 + t272) * t188) * qJD(1), t107 * t46 + t109 * t45 + (t28 * t295 - t10 + t289) * t188 + (-t29 * t295 + t240 - t9) * t189 + t261, -t204 * qJ(3) - t29 * t46 - t28 * t45 + t287 * t93 - t248 * t330 + (-t188 * t29 - t189 * t28) * qJD(4) + t340, t123 * t73 + t321 * t107 + t348 * t188 + (t195 * t24 + (t189 * t302 + t37) * t192) * qJD(1) + t300, t107 * t35 - t109 * t37 + (-t24 * t295 + t289 - t7) * t188 + (t240 + t8 - t279) * t189 + t261, -t257 - t123 * t74 - t321 * t109 - t348 * t189 + (-t192 * t35 - t195 * t25 + (t192 * t302 - t272) * t188) * qJD(1), t11 * t123 - t25 * t35 - t24 * t37 + t321 * t33 - t251 * t330 + (-t188 * t25 + t189 * t24) * qJD(4) - t347 * t246 + t340, -t226 * t237 + t328 * t56, -t13 * t237 - t226 * t238 - t327 * t56 - t328 * t54, t294 * t56 - t341, -t294 * t54 - t342, -t147 * t294 -(-t121 * t191 - t122 * t194) * t111 + t96 * t13 - t6 * t238 + t322 * t54 - g(3) * t221 + t327 * t18 + (t191 * t227 - t194 * t228) * t147 + (t2 * qJD(1) - t238 * t252) * t195 (t121 * t194 - t122 * t191) * t111 - t96 * t226 + t322 * t56 + t328 * t18 + (t191 * t228 + t194 * t227) * t147 - t3 * t294 + (t6 + t347) * t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t285, qJDD(2) + t281, -t197 - t313, qJD(2) * t129 + t145 + t231 - t317, t42, -t339, t213, qJD(1) * t224 - t93 * qJD(2) + t248 + t277, t42, t213, t339, qJD(1) * t225 - t33 * qJD(2) + t251 + t277, 0, 0, 0, 0, 0, qJD(2) * t54 + t341, qJD(2) * t56 + t342; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t203, t41, t346, t29 * t107 + t28 * t109 + t199, t203, t346, -t41, t336 - t323 + t107 * t25 + (-qJD(5) - t24) * t109 + t199, 0, 0, 0, 0, 0, -t13 + t324, t226 - t326; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107 * t109 - t220 (t107 + t293) * t295 + t260, -t106 - t313, -t279 - g(1) * t97 + g(2) * t99 + t33 * t109 + (-t47 - t182) * t189 + t207, 0, 0, 0, 0, 0, -t109 * t54 - t194 * t111 - t191 * t236, -t109 * t56 + t191 * t111 - t194 * t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56 * t54, -t54 ^ 2 + t56 ^ 2, -t226 - t326, -t13 - t324, -t111, -g(1) * t43 - g(2) * t242 - g(3) * t84 - t18 * t56 + t3 * t343 + t276, g(1) * t44 - g(2) * t243 + g(3) * t344 + t18 * t54 + t2 * t343 - t250;];
tau_reg  = t1;
