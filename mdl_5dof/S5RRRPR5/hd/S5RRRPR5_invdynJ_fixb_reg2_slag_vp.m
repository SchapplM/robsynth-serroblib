% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRPR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPR5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR5_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:14:23
% EndTime: 2019-12-31 21:14:33
% DurationCPUTime: 5.56s
% Computational Cost: add. (9649->503), mult. (23285->641), div. (0->0), fcn. (17107->14), ass. (0->263)
t210 = qJ(2) + qJ(3);
t201 = pkin(9) + t210;
t190 = sin(t201);
t215 = sin(qJ(1));
t218 = cos(qJ(1));
t251 = g(1) * t218 + g(2) * t215;
t356 = t251 * t190;
t212 = sin(qJ(5));
t211 = sin(pkin(9));
t301 = cos(pkin(9));
t207 = qJD(2) + qJD(3);
t213 = sin(qJ(3));
t214 = sin(qJ(2));
t293 = t213 * t214;
t247 = t207 * t293;
t217 = cos(qJ(2));
t341 = cos(qJ(3));
t270 = t341 * t217;
t255 = qJD(1) * t270;
t262 = qJDD(1) * t341;
t279 = t217 * qJDD(1);
t256 = -t207 * t255 - t213 * t279 - t214 * t262;
t84 = qJD(1) * t247 + t256;
t149 = t213 * t217 + t341 * t214;
t109 = t207 * t149;
t280 = t214 * qJDD(1);
t246 = t213 * t280 - t217 * t262;
t85 = t109 * qJD(1) + t246;
t53 = -t211 * t84 + t301 * t85;
t50 = qJDD(5) + t53;
t314 = t212 * t50;
t285 = qJD(1) * t214;
t269 = t213 * t285;
t136 = -t255 + t269;
t138 = t149 * qJD(1);
t95 = t301 * t136 + t211 * t138;
t348 = qJD(5) + t95;
t216 = cos(qJ(5));
t352 = t216 * t348;
t355 = -t348 * t352 - t314;
t132 = t138 * qJ(4);
t342 = pkin(7) + pkin(6);
t168 = t342 * t217;
t155 = qJD(1) * t168;
t139 = t213 * t155;
t167 = t342 * t214;
t153 = qJD(1) * t167;
t317 = qJD(2) * pkin(2);
t144 = -t153 + t317;
t99 = t341 * t144 - t139;
t82 = -t132 + t99;
t76 = t207 * pkin(3) + t82;
t271 = t341 * t155;
t100 = t213 * t144 + t271;
t298 = t136 * qJ(4);
t83 = t100 - t298;
t77 = t301 * t83;
t47 = t211 * t76 + t77;
t45 = t207 * pkin(8) + t47;
t331 = t217 * pkin(2);
t198 = pkin(1) + t331;
t166 = t198 * qJD(1);
t112 = t136 * pkin(3) + qJD(4) - t166;
t236 = -t211 * t136 + t301 * t138;
t56 = t95 * pkin(4) - pkin(8) * t236 + t112;
t20 = -t212 * t45 + t216 * t56;
t354 = t20 * t348;
t21 = t212 * t56 + t216 * t45;
t353 = t21 * t348;
t202 = sin(t210);
t203 = cos(t210);
t351 = -g(3) * t203 + t251 * t202;
t350 = t95 * t236;
t282 = qJD(5) * t216;
t315 = t211 * t83;
t46 = t301 * t76 - t315;
t44 = -t207 * pkin(4) - t46;
t205 = qJDD(2) + qJDD(3);
t281 = qJD(1) * qJD(2);
t266 = t217 * t281;
t111 = qJDD(2) * pkin(2) + t342 * (-t266 - t280);
t267 = t214 * t281;
t113 = t342 * (-t267 + t279);
t62 = -t100 * qJD(3) + t341 * t111 - t213 * t113;
t34 = t205 * pkin(3) + t84 * qJ(4) - t138 * qJD(4) + t62;
t268 = qJD(3) * t341;
t284 = qJD(3) * t213;
t257 = -t213 * t111 - t341 * t113 - t144 * t268 + t155 * t284;
t38 = -qJ(4) * t85 - qJD(4) * t136 - t257;
t10 = -t211 * t38 + t301 * t34;
t8 = -t205 * pkin(4) - t10;
t346 = t8 * t212 + t44 * t282;
t283 = qJD(5) * t212;
t345 = t216 * t50 - t283 * t348;
t343 = g(1) * t215 - g(2) * t218;
t344 = t343 * t190;
t115 = -t213 * t167 + t341 * t168;
t191 = cos(t201);
t252 = t191 * pkin(4) + t190 * pkin(8);
t340 = pkin(3) * t202;
t339 = pkin(4) * t190;
t178 = g(3) * t190;
t335 = g(3) * t191;
t333 = g(3) * t217;
t332 = t138 * pkin(3);
t54 = -t211 * t85 - t301 * t84;
t133 = pkin(2) * t267 - t198 * qJDD(1);
t68 = t85 * pkin(3) + qJDD(4) + t133;
t17 = t53 * pkin(4) - t54 * pkin(8) + t68;
t11 = t211 * t34 + t301 * t38;
t9 = t205 * pkin(8) + t11;
t2 = qJD(5) * t20 + t212 * t17 + t216 * t9;
t1 = t2 * t216;
t16 = t216 * t17;
t3 = -qJD(5) * t21 - t212 * t9 + t16;
t330 = t3 * t212;
t329 = t44 * t95;
t86 = -t216 * t207 + t212 * t236;
t328 = t86 * t95;
t88 = t212 * t207 + t216 * t236;
t327 = t88 * t86;
t326 = t88 * t236;
t325 = t348 * t236;
t324 = t236 * t86;
t323 = t236 ^ 2;
t321 = t95 ^ 2;
t260 = -t216 * t205 + t212 * t54;
t37 = t88 * qJD(5) + t260;
t319 = -t212 * t37 - t86 * t282;
t318 = pkin(2) * qJD(3);
t316 = t20 * t216;
t108 = -qJD(2) * t270 - t217 * t268 + t247;
t72 = -t301 * t108 - t211 * t109;
t313 = t212 * t72;
t312 = t212 * t86;
t311 = t212 * t88;
t310 = t212 * t95;
t309 = t216 * t88;
t36 = -t212 * t205 - t207 * t282 - t216 * t54 + t236 * t283;
t307 = t36 * t212;
t306 = t37 * t216;
t305 = t236 * t207;
t304 = t95 * t207;
t106 = t213 * t153 - t271;
t231 = t106 + t298;
t261 = t301 * t213;
t107 = -t341 * t153 - t139;
t89 = -t132 + t107;
t303 = -t211 * t89 + t301 * t231 + (t341 * t211 + t261) * t318;
t294 = t211 * t213;
t129 = (t301 * t341 - t294) * t318;
t59 = t211 * t231 + t301 * t89;
t302 = t129 - t59;
t300 = pkin(6) * qJDD(1);
t299 = qJD(5) * t348;
t297 = t138 * t136;
t296 = t191 * t215;
t295 = t191 * t218;
t292 = t215 * t212;
t291 = t215 * t216;
t290 = t218 * t212;
t289 = t218 * t216;
t197 = t341 * pkin(2) + pkin(3);
t131 = pkin(2) * t261 + t211 * t197;
t208 = t214 ^ 2;
t209 = t217 ^ 2;
t287 = t208 - t209;
t286 = t208 + t209;
t200 = t214 * t317;
t126 = pkin(8) + t131;
t278 = t126 * t299;
t192 = t211 * pkin(3) + pkin(8);
t277 = t192 * t299;
t221 = qJD(1) ^ 2;
t275 = t214 * t221 * t217;
t42 = t44 * t283;
t274 = g(1) * t295 + g(2) * t296 + t178;
t273 = -t8 - t335;
t272 = qJD(2) * t342;
t195 = pkin(3) * t203;
t265 = t195 + t331;
t98 = t109 * pkin(3) + t200;
t157 = -t214 * pkin(2) - t340;
t264 = t157 - t339;
t259 = t212 * t348;
t254 = t214 * t266;
t253 = -t339 - t340;
t249 = t236 * t47 - t46 * t95;
t248 = -t192 * t50 + t329;
t245 = -t21 * t212 - t316;
t244 = t20 * t212 - t21 * t216;
t148 = -t270 + t293;
t102 = t301 * t148 + t211 * t149;
t103 = -t211 * t148 + t301 * t149;
t117 = t148 * pkin(3) - t198;
t64 = t102 * pkin(4) - t103 * pkin(8) + t117;
t114 = -t341 * t167 - t213 * t168;
t230 = -t149 * qJ(4) + t114;
t93 = -qJ(4) * t148 + t115;
t66 = t211 * t230 + t301 * t93;
t28 = -t212 * t66 + t216 * t64;
t29 = t212 * t64 + t216 * t66;
t243 = t21 * t236 + t212 * t335 + t346;
t242 = -t20 * t236 + t216 * t356 + t42;
t241 = -t310 * t348 + t345;
t240 = -t21 * t310 - t316 * t95 + t1 - t274;
t239 = -qJD(5) * t56 + t178 - t9;
t237 = -0.2e1 * pkin(1) * t281 - pkin(6) * qJDD(2);
t63 = pkin(4) * t236 + pkin(8) * t95 + t332;
t154 = t214 * t272;
t156 = t217 * t272;
t74 = -t341 * t154 - t213 * t156 - t167 * t268 - t168 * t284;
t233 = -t126 * t50 - t129 * t348 + t329;
t232 = t112 * t95 - t11 + t274;
t130 = -pkin(2) * t294 + t301 * t197;
t220 = qJD(2) ^ 2;
t229 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t220 + t343;
t228 = pkin(1) * t221 + t251 - t300;
t227 = g(3) * t202 - t166 * t136 + t251 * t203 + t257;
t226 = -t112 * t236 + t10 - t335 + t356;
t225 = t245 * qJD(5) + t1 - t330;
t75 = -qJD(3) * t115 + t213 * t154 - t341 * t156;
t224 = t166 * t138 + t351 + t62;
t223 = t108 * qJ(4) - t149 * qJD(4) + t75;
t206 = -qJ(4) - t342;
t199 = pkin(2) * t285;
t193 = -t301 * pkin(3) - pkin(4);
t159 = pkin(8) * t295;
t158 = pkin(8) * t296;
t152 = pkin(1) + t265;
t143 = t218 * t152;
t125 = -pkin(4) - t130;
t124 = t191 * t289 + t292;
t123 = -t191 * t290 + t291;
t122 = -t191 * t291 + t290;
t121 = t191 * t292 + t289;
t116 = t199 + t332;
t90 = -t136 ^ 2 + t138 ^ 2;
t71 = -t211 * t108 + t301 * t109;
t69 = -t256 + (t136 - t269) * t207;
t65 = t211 * t93 - t301 * t230;
t60 = t199 + t63;
t57 = -t109 * qJ(4) - t148 * qJD(4) + t74;
t52 = t301 * t82 - t315;
t51 = t211 * t82 + t77;
t41 = t54 + t304;
t40 = -t53 + t305;
t39 = -t321 + t323;
t32 = t71 * pkin(4) - t72 * pkin(8) + t98;
t27 = t212 * t60 + t216 * t59;
t26 = -t212 * t59 + t216 * t60;
t25 = t212 * t63 + t216 * t52;
t24 = -t212 * t52 + t216 * t63;
t23 = t211 * t223 + t301 * t57;
t22 = t211 * t57 - t301 * t223;
t15 = t86 * t259 - t306;
t14 = t352 * t88 - t307;
t13 = -t326 - t355;
t12 = t241 + t324;
t6 = -t29 * qJD(5) - t212 * t23 + t216 * t32;
t5 = t28 * qJD(5) + t212 * t32 + t216 * t23;
t4 = (-t36 - t328) * t216 - t348 * t311 + t319;
t7 = [0, 0, 0, 0, 0, qJDD(1), t343, t251, 0, 0, t208 * qJDD(1) + 0.2e1 * t254, 0.2e1 * t214 * t279 - 0.2e1 * t287 * t281, qJDD(2) * t214 + t220 * t217, t209 * qJDD(1) - 0.2e1 * t254, qJDD(2) * t217 - t220 * t214, 0, t237 * t214 + t229 * t217, -t229 * t214 + t237 * t217, 0.2e1 * t286 * t300 - t251, -g(1) * (-t215 * pkin(1) + t218 * pkin(6)) - g(2) * (t218 * pkin(1) + t215 * pkin(6)) + (t286 * pkin(6) ^ 2 + pkin(1) ^ 2) * qJDD(1), -t108 * t138 - t149 * t84, t108 * t136 - t109 * t138 + t148 * t84 - t149 * t85, -t108 * t207 + t149 * t205, t109 * t136 + t148 * t85, -t109 * t207 - t148 * t205, 0, -t166 * t109 + t114 * t205 + t133 * t148 + t136 * t200 - t198 * t85 + t203 * t343 + t75 * t207, t166 * t108 - t115 * t205 + t133 * t149 + t138 * t200 + t198 * t84 - t202 * t343 - t74 * t207, -t100 * t109 + t99 * t108 + t114 * t84 - t115 * t85 - t74 * t136 - t75 * t138 + t148 * t257 - t62 * t149 - t251, -t257 * t115 + t100 * t74 + t62 * t114 + t99 * t75 - t133 * t198 - t166 * t200 - g(1) * (-t215 * t198 + t218 * t342) - g(2) * (t218 * t198 + t215 * t342), t103 * t54 + t236 * t72, -t102 * t54 - t103 * t53 - t236 * t71 - t72 * t95, t103 * t205 + t72 * t207, t102 * t53 + t71 * t95, -t102 * t205 - t71 * t207, 0, t68 * t102 + t112 * t71 + t117 * t53 + t191 * t343 - t65 * t205 - t22 * t207 + t98 * t95, t68 * t103 + t112 * t72 + t117 * t54 - t66 * t205 - t23 * t207 + t236 * t98 - t344, -t10 * t103 - t11 * t102 + t22 * t236 - t23 * t95 - t46 * t72 - t47 * t71 - t66 * t53 + t65 * t54 - t251, t11 * t66 + t47 * t23 - t10 * t65 - t46 * t22 + t68 * t117 + t112 * t98 - g(1) * (-t215 * t152 - t218 * t206) - g(2) * (-t215 * t206 + t143), t72 * t309 + (-t216 * t36 - t283 * t88) * t103, (-t216 * t86 - t311) * t72 + (t307 - t306 + (-t309 + t312) * qJD(5)) * t103, -t36 * t102 + t103 * t345 + t352 * t72 + t88 * t71, -t319 * t103 + t72 * t312, -t348 * t313 - t37 * t102 - t86 * t71 + (-t282 * t348 - t314) * t103, t102 * t50 + t348 * t71, -g(1) * t122 - g(2) * t124 + t3 * t102 + t103 * t346 + t20 * t71 + t22 * t86 + t28 * t50 + t44 * t313 + t348 * t6 + t65 * t37, t44 * t216 * t72 - g(1) * t121 - g(2) * t123 - t2 * t102 - t21 * t71 + t22 * t88 - t29 * t50 - t65 * t36 - t5 * t348 + (t8 * t216 - t42) * t103, t28 * t36 - t29 * t37 - t5 * t86 - t6 * t88 + t245 * t72 + t344 + (qJD(5) * t244 - t2 * t212 - t3 * t216) * t103, -g(2) * t143 + t2 * t29 + t20 * t6 + t21 * t5 + t44 * t22 + t3 * t28 + t8 * t65 + (g(1) * t206 - g(2) * t252) * t218 + (-g(1) * (-t152 - t252) + g(2) * t206) * t215; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t275, t287 * t221, t280, t275, t279, qJDD(2), t228 * t214 - t333, g(3) * t214 + t228 * t217, 0, 0, t297, t90, t69, -t297, -t246, t205, -t106 * t207 + (-t136 * t285 + t341 * t205 - t207 * t284) * pkin(2) + t224, t107 * t207 + (-t138 * t285 - t205 * t213 - t207 * t268) * pkin(2) + t227, (t100 + t106) * t138 + (t107 - t99) * t136 + (t341 * t84 - t213 * t85 + (-t341 * t136 + t138 * t213) * qJD(3)) * pkin(2), -t100 * t107 - t99 * t106 + (t341 * t62 - t333 - t213 * t257 + (t341 * t100 - t213 * t99) * qJD(3) + (qJD(1) * t166 + t251) * t214) * pkin(2), t350, t39, t41, -t350, t40, t205, -t116 * t95 + t130 * t205 - t303 * t207 + t226, -t116 * t236 - t131 * t205 - t302 * t207 + t232, -t130 * t54 - t131 * t53 + t236 * t303 - t302 * t95 + t249, -g(3) * t265 + t10 * t130 + t11 * t131 - t112 * t116 - t251 * t157 + t302 * t47 - t303 * t46, t14, t4, t13, t15, t12, -t325, t125 * t37 - t26 * t348 + t303 * t86 + (t273 - t278) * t216 + t233 * t212 + t242, -t125 * t36 + t27 * t348 + t303 * t88 + t233 * t216 + (-t356 + t278) * t212 + t243, t26 * t88 + t27 * t86 + (-t126 * t37 - t129 * t86 + (t126 * t88 - t20) * qJD(5)) * t216 + (-t126 * t36 + t129 * t88 - t3 + (t126 * t86 - t21) * qJD(5)) * t212 + t240, t8 * t125 - t21 * t27 - t20 * t26 - g(1) * (t218 * t264 + t159) - g(2) * (t215 * t264 + t158) - g(3) * (t265 + t252) + t303 * t44 - t244 * t129 + t225 * t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t297, t90, t69, -t297, -t246, t205, t100 * t207 + t224, t99 * t207 + t227, 0, 0, t350, t39, t41, -t350, t40, t205, t51 * t207 + (-t138 * t95 + t301 * t205) * pkin(3) + t226, t52 * t207 + (-t138 * t236 - t205 * t211) * pkin(3) + t232, -t51 * t236 + t52 * t95 + (-t211 * t53 - t301 * t54) * pkin(3) + t249, t46 * t51 - t47 * t52 + (t301 * t10 + t11 * t211 - t112 * t138 + t351) * pkin(3), t14, t4, t13, t15, t12, -t325, t193 * t37 - t24 * t348 - t51 * t86 + t248 * t212 + (t273 - t277) * t216 + t242, -t193 * t36 + t25 * t348 - t51 * t88 + t248 * t216 + (-t356 + t277) * t212 + t243, -t330 + t24 * t88 + t25 * t86 + (-t306 - t307) * t192 + ((t309 + t312) * t192 + t245) * qJD(5) + t240, t8 * t193 - t21 * t25 - t20 * t24 - t44 * t51 - g(1) * (t218 * t253 + t159) - g(2) * (t215 * t253 + t158) - g(3) * (t195 + t252) + t225 * t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53 + t305, t54 - t304, -t321 - t323, t236 * t46 + t47 * t95 - t343 + t68, 0, 0, 0, 0, 0, 0, t241 - t324, -t326 + t355, (t36 - t328) * t216 + t88 * t259 + t319, -t44 * t236 + (t3 + t353) * t216 + (t2 - t354) * t212 - t343; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t327, -t86 ^ 2 + t88 ^ 2, t348 * t86 - t36, -t327, -t260 + (-qJD(5) + t348) * t88, t50, -g(1) * t123 + g(2) * t121 + t212 * t239 - t282 * t45 - t44 * t88 + t16 + t353, g(1) * t124 - g(2) * t122 + t354 + t44 * t86 + (qJD(5) * t45 - t17) * t212 + t239 * t216, 0, 0;];
tau_reg = t7;
