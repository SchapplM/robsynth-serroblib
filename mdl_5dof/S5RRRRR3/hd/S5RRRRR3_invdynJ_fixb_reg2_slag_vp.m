% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRRR3
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:56:35
% EndTime: 2019-12-05 18:56:47
% DurationCPUTime: 5.65s
% Computational Cost: add. (7925->511), mult. (18252->682), div. (0->0), fcn. (14252->14), ass. (0->254)
t192 = sin(qJ(3));
t197 = cos(qJ(3));
t198 = cos(qJ(2));
t284 = qJD(1) * t198;
t193 = sin(qJ(2));
t285 = qJD(1) * t193;
t126 = t192 * t285 - t197 * t284;
t191 = sin(qJ(4));
t195 = cos(qJ(5));
t296 = t191 * t195;
t190 = sin(qJ(5));
t196 = cos(qJ(4));
t298 = t190 * t196;
t135 = t296 + t298;
t268 = qJD(4) + qJD(5);
t356 = t268 * t135;
t358 = t135 * t126 + t356;
t274 = qJD(5) * t195;
t277 = qJD(4) * t196;
t299 = t190 * t191;
t345 = t195 * t196 - t299;
t357 = t345 * t126 + t195 * t277 + t196 * t274 - t268 * t299;
t136 = t192 * t198 + t193 * t197;
t128 = t136 * qJD(1);
t185 = qJD(2) + qJD(3);
t100 = t128 * t196 + t185 * t191;
t98 = t128 * t191 - t196 * t185;
t233 = t195 * t100 - t190 * t98;
t62 = t100 * t190 + t195 * t98;
t331 = t62 * t233;
t341 = pkin(1) * t192;
t178 = qJD(2) * t341;
t138 = pkin(5) * t185 + t178;
t266 = pkin(1) * t284;
t81 = pkin(2) * t126 - pkin(5) * t128 - t266;
t66 = -t138 * t191 + t196 * t81;
t355 = t66 * qJD(4);
t189 = qJ(2) + qJ(3);
t180 = sin(t189);
t199 = cos(qJ(1));
t302 = t180 * t199;
t194 = sin(qJ(1));
t303 = t180 * t194;
t354 = g(1) * t302 + g(2) * t303;
t353 = t233 ^ 2 - t62 ^ 2;
t182 = cos(t189);
t335 = g(3) * t182;
t352 = -t335 + t354;
t123 = qJD(4) + t126;
t117 = qJD(5) + t123;
t184 = qJDD(2) + qJDD(3);
t294 = t192 * t193;
t134 = -t197 * t198 + t294;
t224 = t134 * qJD(3);
t203 = qJDD(1) * t136 + (-qJD(2) * t134 - t224) * qJD(1);
t278 = qJD(4) * t191;
t260 = t128 * t277 + t185 * t278 + t191 * t203;
t231 = t196 * t184 - t260;
t276 = qJD(5) * t190;
t47 = t128 * t278 - t191 * t184 - t185 * t277 - t196 * t203;
t15 = t100 * t276 - t190 * t231 + t195 * t47 + t98 * t274;
t351 = t117 * t62 - t15;
t188 = qJ(4) + qJ(5);
t179 = sin(t188);
t181 = cos(t188);
t301 = t182 * t194;
t107 = t179 * t199 - t181 * t301;
t300 = t182 * t199;
t109 = t179 * t194 + t181 * t300;
t171 = g(3) * t180;
t271 = qJDD(2) * t192;
t280 = qJD(3) * t197;
t116 = pkin(5) * t184 + (qJD(2) * t280 + t271) * pkin(1);
t332 = t193 * pkin(1);
t183 = t198 * pkin(1);
t347 = -pkin(5) * t136 - t183;
t269 = t198 * qJDD(1);
t270 = t193 * qJDD(1);
t236 = t192 * t270 - t197 * t269;
t95 = t185 * t136;
t315 = qJD(1) * t95;
t73 = t236 + t315;
t36 = t73 * pkin(2) + t347 * qJDD(1) + (pkin(5) * t224 + (pkin(5) * t134 + t332) * qJD(2)) * qJD(1);
t35 = t196 * t36;
t67 = t138 * t196 + t191 * t81;
t19 = -qJD(4) * t67 - t116 * t191 + t35;
t71 = qJDD(4) + t73;
t12 = pkin(3) * t71 + t19;
t18 = t116 * t196 + t191 * t36 + t355;
t50 = pkin(3) * t123 + t66;
t4 = (qJD(5) * t50 + t18) * t195 + t12 * t190 - t67 * t276;
t282 = qJD(2) * t197;
t265 = pkin(1) * t282;
t139 = -pkin(2) * t185 - t265;
t82 = pkin(3) * t98 + t139;
t350 = g(1) * t109 - g(2) * t107 + t181 * t171 + t62 * t82 - t4;
t174 = pkin(5) + t341;
t87 = pkin(2) * t128 + pkin(5) * t126;
t80 = pkin(1) * t285 + t87;
t348 = (qJD(4) * t174 + t80) * t123;
t308 = t126 * t191;
t346 = (t278 + t308) * pkin(3);
t291 = t196 * t199;
t297 = t191 * t194;
t118 = t182 * t297 + t291;
t293 = t194 * t196;
t295 = t191 * t199;
t120 = -t182 * t295 + t293;
t344 = -g(1) * t120 + g(2) * t118;
t106 = t179 * t301 + t181 * t199;
t108 = -t179 * t300 + t181 * t194;
t319 = t195 * t67;
t29 = t190 * t50 + t319;
t5 = -qJD(5) * t29 + t195 * t12 - t18 * t190;
t343 = -g(1) * t108 + g(2) * t106 + t179 * t171 - t233 * t82 + t5;
t16 = qJD(5) * t233 - t190 * t47 - t195 * t231;
t342 = t117 * t233 - t16;
t340 = pkin(1) * t197;
t337 = g(1) * t194;
t334 = g(3) * t191;
t169 = t180 * pkin(5);
t170 = t182 * pkin(2);
t333 = t184 * pkin(2);
t263 = pkin(1) * t280;
t122 = t128 * pkin(3);
t60 = t196 * t80 + t122;
t330 = -t174 * t356 - t190 * t60 + t263 * t345 - t80 * t296;
t207 = t268 * t345;
t329 = -t135 * t263 - t174 * t207 - t195 * t60 + t80 * t299;
t75 = -t191 * t265 + t196 * t87;
t58 = t122 + t75;
t76 = t191 * t87 + t196 * t265;
t326 = -pkin(5) * t356 - t190 * t58 - t195 * t76;
t325 = -t207 * pkin(5) + t190 * t76 - t195 * t58;
t324 = t100 * t98;
t323 = t123 * t98;
t17 = t18 * t196;
t322 = t190 * t67;
t321 = t191 * t71;
t238 = t185 * t294;
t94 = t238 + (-t280 - t282) * t198;
t320 = t191 * t94;
t318 = t47 * t191;
t317 = t66 * t196;
t316 = t98 * t191;
t314 = t100 * t123;
t313 = t100 * t196;
t312 = t117 * t128;
t311 = t123 * t128;
t310 = t126 * t139;
t309 = t126 * t185;
t307 = t128 * t126;
t305 = t136 * t196;
t175 = pkin(3) * t196 + pkin(2);
t142 = t182 * t175;
t290 = t142 + t169;
t289 = pkin(5) * t302 + t199 * t183;
t288 = t170 + t169;
t287 = -qJD(3) * t178 + qJDD(2) * t340;
t186 = t193 ^ 2;
t187 = t198 ^ 2;
t286 = t186 - t187;
t283 = qJD(2) * t193;
t281 = qJD(3) * t192;
t279 = qJD(4) * t123;
t275 = qJD(5) * t191;
t273 = qJD(1) * qJD(2);
t272 = qJDD(1) * t187;
t264 = pkin(1) * t281;
t261 = pkin(5) * t279;
t201 = qJD(1) ^ 2;
t259 = t193 * t201 * t198;
t258 = g(1) * t300 + g(2) * t301 + t171;
t257 = t136 * t278;
t256 = t136 * t277;
t255 = t193 * t273;
t115 = -t287 - t333;
t254 = -t115 - t335;
t249 = t123 * t196;
t248 = t198 * t185;
t247 = -t126 * t266 + t258;
t246 = t198 * t255;
t245 = t264 + t346;
t244 = -g(1) * t303 + g(2) * t302;
t243 = -t178 + t346;
t57 = pkin(1) * t283 + pkin(2) * t95 + pkin(5) * t94;
t90 = pkin(2) * t134 + t347;
t242 = pkin(3) * t95 + t196 * t57 + (-t275 - t278) * t90;
t241 = -t183 - t169;
t147 = g(1) * t199 + g(2) * t194;
t240 = -g(2) * t199 + t337;
t239 = -pkin(5) * t71 + t310;
t237 = -t175 * t180 - t332;
t235 = t67 * t191 + t317;
t234 = -t191 * t66 + t196 * t67;
t232 = t115 * t191 + t128 * t67 + t139 * t277 + t182 * t334;
t230 = -t66 * t128 + t139 * t278 + t354 * t196;
t229 = -t126 * t317 - t67 * t308 + t17 - t258;
t228 = -qJD(4) * t81 - t116 + t171;
t227 = t147 * t180;
t226 = t147 * t193;
t225 = -t196 * t94 - t257;
t222 = t231 * t196;
t221 = t128 * t266 + t287 + t352;
t220 = qJD(4) * t136 * t139 + t123 * t57 + t71 * t90;
t72 = pkin(3) * t134 + t196 * t90;
t219 = qJD(5) * t72 + t191 * t57 + t277 * t90;
t218 = t115 * t136 - t139 * t94 - t279 * t90;
t217 = -g(3) * t198 + t226;
t216 = -t123 * t263 - t174 * t71 + t310;
t215 = -qJD(4) * t235 - t19 * t191;
t28 = t195 * t50 - t322;
t213 = -t5 * t135 - t357 * t28 - t358 * t29 + t345 * t4 - t258;
t41 = -pkin(3) * t231 + t115;
t212 = -t128 * t28 + t352 * t181 - t345 * t41 + t358 * t82;
t211 = t215 + t17;
t148 = pkin(5) * t301;
t150 = pkin(5) * t300;
t210 = -g(1) * (-pkin(2) * t302 + t150) - g(2) * (-pkin(2) * t303 + t148);
t206 = t128 * t29 + t41 * t135 + t357 * t82 + (-t227 + t335) * t179;
t202 = pkin(1) ^ 2;
t200 = qJD(2) ^ 2;
t176 = -pkin(2) - t340;
t146 = -t175 - t340;
t130 = t345 * pkin(5);
t129 = t135 * pkin(5);
t121 = t182 * t291 + t297;
t119 = -t182 * t293 + t295;
t103 = t345 * t174;
t102 = t135 * t174;
t85 = t345 * t136;
t84 = t135 * t136;
t74 = -t126 ^ 2 + t128 ^ 2;
t70 = qJDD(5) + t71;
t53 = t203 + t309;
t43 = t190 * t72 + t90 * t296;
t42 = t195 * t72 - t90 * t299;
t32 = t195 * t66 - t322;
t31 = -t190 * t66 - t319;
t27 = -t94 * t298 + (t268 * t305 - t320) * t195 + (-t136 * t275 - t257) * t190;
t26 = t136 * t356 + t345 * t94;
t25 = -t100 * t128 + t123 * t249 + t321;
t24 = -t123 ^ 2 * t191 + t128 * t98 + t196 * t71;
t23 = t123 * t316 + t222;
t22 = t100 * t249 - t318;
t14 = -t117 * t358 + t128 * t62 + t345 * t70;
t13 = t117 * t357 - t128 * t233 + t135 * t70;
t10 = (-t47 - t323) * t196 + (t231 - t314) * t191;
t9 = -t190 * t219 + t195 * t242;
t8 = t190 * t242 + t195 * t219;
t7 = -t16 * t345 + t358 * t62;
t6 = -t135 * t15 + t233 * t357;
t1 = -t135 * t16 - t15 * t345 - t233 * t358 - t357 * t62;
t2 = [0, 0, 0, 0, 0, qJDD(1), t240, t147, 0, 0, qJDD(1) * t186 + 0.2e1 * t246, 0.2e1 * t193 * t269 - 0.2e1 * t273 * t286, qJDD(2) * t193 + t198 * t200, -0.2e1 * t246 + t272, qJDD(2) * t198 - t193 * t200, 0, t240 * t198, -t240 * t193, -t147, 0, -t128 * t94 + t203 * t136, t94 * t126 - t128 * t95 - t203 * t134 - t136 * t73, t136 * t184 - t185 * t94, t126 * t95 + t134 * t73, -t134 * t184 - t185 * t95, 0, t240 * t182 + ((qJD(1) * t134 + t126) * t283 + (-qJDD(1) * t134 - t315 - t73) * t198) * pkin(1), (t128 * t283 - 0.2e1 * t136 * t269 + (t136 * t283 + (-t197 * t248 + t238 + t94) * t198) * qJD(1)) * pkin(1) + t244, ((-t134 * t192 - t136 * t197) * qJDD(2) + (-t192 * t95 + t197 * t94 + (-t134 * t197 + t136 * t192) * qJD(3)) * qJD(2)) * pkin(1) - t147, t202 * t272 + (pkin(1) * t240 - 0.2e1 * t202 * t255) * t198, t100 * t225 - t305 * t47, (t100 * t191 + t196 * t98) * t94 + (t222 + t318 + (-t313 + t316) * qJD(4)) * t136, t100 * t95 + t123 * t225 - t134 * t47 + t305 * t71, -t94 * t316 + (-t191 * t231 + t277 * t98) * t136, -t136 * t321 + t231 * t134 - t98 * t95 + (-t256 + t320) * t123, t123 * t95 + t134 * t71, -g(1) * t119 - g(2) * t121 + t134 * t19 + t191 * t218 + t196 * t220 + t66 * t95, -g(1) * t118 - g(2) * t120 - t134 * t18 - t191 * t220 + t196 * t218 - t67 * t95, (-t100 * t57 - t136 * t19 + t47 * t90 + t66 * t94 + (-t136 * t67 - t90 * t98) * qJD(4)) * t196 + (-t57 * t98 + t90 * t231 - t18 * t136 + t67 * t94 + (t90 * t100 + t66 * t136) * qJD(4)) * t191 - t244, -g(2) * (pkin(2) * t300 + t289) + t235 * t57 - (t241 - t170) * t337 + (qJD(4) * t234 + t18 * t191 + t19 * t196) * t90, -t15 * t85 - t233 * t26, t15 * t84 - t16 * t85 - t233 * t27 + t26 * t62, -t117 * t26 - t134 * t15 + t233 * t95 + t70 * t85, t16 * t84 + t27 * t62, -t117 * t27 - t134 * t16 - t62 * t95 - t70 * t84, t117 * t95 + t134 * t70, -g(1) * t107 - g(2) * t109 + t117 * t9 + t134 * t5 + t27 * t82 + t28 * t95 + t41 * t84 + t42 * t70 + (-t62 * t320 + (t16 * t191 + t277 * t62) * t136) * pkin(3), -g(1) * t106 - g(2) * t108 - t117 * t8 - t134 * t4 - t26 * t82 - t29 * t95 + t41 * t85 - t43 * t70 + (-t233 * t320 + (-t15 * t191 + t233 * t277) * t136) * pkin(3), t15 * t42 - t16 * t43 - t233 * t9 + t26 * t28 - t27 * t29 - t4 * t84 - t5 * t85 - t62 * t8 - t244, t4 * t43 + t29 * t8 + t5 * t42 + t28 * t9 - g(2) * (t175 * t300 + t289) - (t241 - t142) * t337 + (t82 * t256 + (t136 * t41 - t82 * t94 - t147) * t191) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t259, t286 * t201, t270, t259, t269, qJDD(2), t217, g(3) * t193 + t147 * t198, 0, 0, t307, t74, t53, -t307, -t236, t184, (-t126 * t285 + t184 * t197 - t185 * t281) * pkin(1) + t221, (-t128 * t285 + (-qJDD(2) - t184) * t192 + (-qJD(2) - t185) * t280) * pkin(1) + t247, (((-qJD(1) * t248 - t270) * t197 - t309) * t197 + (-t73 + t197 * (qJD(3) * t285 + t255 - t269) + t128 * t185) * t192) * pkin(1), (t259 + (t192 ^ 2 + t197 ^ 2) * qJDD(2)) * t202 + t217 * pkin(1), t22, t10, t25, t23, t24, -t311, t98 * t264 + t176 * t260 + t216 * t191 + (-t176 * t184 + t254 - t348) * t196 + t230, t100 * t264 - t176 * t47 + t216 * t196 + (-t227 + t348) * t191 + t232, (-t98 * t263 - t355 + t80 * t100 + (qJD(4) * t100 + t231) * t174) * t196 + (t100 * t263 - t174 * t47 + t80 * t98 - t19 + (t174 * t98 - t67) * qJD(4)) * t191 + t229, t115 * t176 - g(3) * (t183 + t288) - t235 * t80 + (t226 + (t139 * t192 + t197 * t234) * qJD(3)) * pkin(1) + t211 * t174 + t210, t6, t1, t13, t7, t14, -t312, -t102 * t70 + t117 * t329 + t146 * t16 + t245 * t62 + t212, -t103 * t70 - t117 * t330 - t146 * t15 + t233 * t245 + t206, -t102 * t15 - t103 * t16 - t233 * t329 - t330 * t62 + t213, t4 * t103 - t5 * t102 + t41 * t146 - g(1) * (t199 * t237 + t150) - g(2) * (t194 * t237 + t148) - g(3) * (t183 + t290) + t245 * t82 + t330 * t29 + t329 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t307, t74, t53, -t307, -t236, t184, t178 * t185 + t221, (-t271 + (-qJD(3) + t185) * t282) * pkin(1) + t247, 0, 0, t22, t10, t25, t23, t24, -t311, -pkin(2) * t260 - t75 * t123 - t98 * t178 + t239 * t191 + (t254 - t261 + t333) * t196 + t230, -t100 * t178 + pkin(2) * t47 + t123 * t76 + t239 * t196 + (-t227 + t261) * t191 + t232, t75 * t100 + t76 * t98 + (t222 - t318 + (t313 + t316) * qJD(4)) * pkin(5) + t215 + t229, -t115 * pkin(2) + pkin(5) * t211 - g(3) * t288 - t139 * t178 - t66 * t75 - t67 * t76 + t210, t6, t1, t13, t7, t14, -t312, t117 * t325 - t129 * t70 - t16 * t175 + t243 * t62 + t212, -t117 * t326 - t130 * t70 + t15 * t175 + t233 * t243 + t206, -t129 * t15 - t130 * t16 - t233 * t325 - t326 * t62 + t213, t4 * t130 - t5 * t129 - t41 * t175 - g(1) * (-t175 * t302 + t150) - g(2) * (-t175 * t303 + t148) - g(3) * t290 + t243 * t82 + t326 * t29 + t325 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t324, t100 ^ 2 - t98 ^ 2, -t47 + t323, -t324, t231 + t314, t71, -t100 * t139 + t123 * t67 - t138 * t277 + t191 * t228 + t344 + t35, g(1) * t121 - g(2) * t119 + t123 * t66 + t139 * t98 + (qJD(4) * t138 - t36) * t191 + t228 * t196, 0, 0, t331, t353, t351, -t331, t342, t70, -t117 * t31 + (-t100 * t62 - t117 * t276 + t195 * t70) * pkin(3) + t343, t117 * t32 + (-t100 * t233 - t117 * t274 - t190 * t70) * pkin(3) + t350, -t28 * t62 + t29 * t233 + t31 * t233 + t32 * t62 + (t15 * t195 - t16 * t190 + (t190 * t233 - t195 * t62) * qJD(5)) * pkin(3), -t28 * t31 - t29 * t32 + (t4 * t190 + t5 * t195 - t82 * t100 + t180 * t334 + (-t190 * t28 + t195 * t29) * qJD(5) + t344) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t331, t353, t351, -t331, t342, t70, t117 * t29 + t343, t117 * t28 + t350, 0, 0;];
tau_reg = t2;
