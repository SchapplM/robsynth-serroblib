% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRPP7
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPP7_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP7_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP7_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP7_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:05:59
% EndTime: 2019-12-31 21:06:10
% DurationCPUTime: 5.38s
% Computational Cost: add. (3387->511), mult. (7613->587), div. (0->0), fcn. (4817->6), ass. (0->237)
t168 = cos(qJ(3));
t166 = sin(qJ(2));
t266 = qJD(1) * t166;
t242 = t168 * t266;
t165 = sin(qJ(3));
t264 = qJD(2) * t165;
t102 = t242 + t264;
t169 = cos(qJ(2));
t265 = qJD(1) * t169;
t340 = -t265 + qJD(3);
t289 = t102 * t340;
t252 = t166 * qJDD(1);
t260 = qJD(3) * t166;
t339 = qJD(1) * t260 - qJDD(2);
t43 = ((qJD(3) + t265) * qJD(2) + t252) * t165 + t339 * t168;
t347 = t43 - t289;
t346 = t43 + t289;
t243 = t165 * t266;
t256 = t168 * qJD(2);
t100 = t243 - t256;
t291 = t100 * t340;
t253 = qJD(1) * qJD(2);
t239 = t169 * t253;
t341 = -t239 - t252;
t42 = -qJD(3) * t256 + t339 * t165 + t341 * t168;
t343 = t42 + t291;
t4 = t346 * t165 + t343 * t168;
t153 = t169 * qJDD(1);
t328 = -t166 * t253 + t153;
t98 = qJDD(3) - t328;
t344 = t98 * qJ(4) + qJD(4) * t340;
t170 = cos(qJ(1));
t283 = t166 * t170;
t167 = sin(qJ(1));
t285 = t166 * t167;
t331 = g(1) * t283 + g(2) * t285;
t262 = qJD(2) * t169;
t296 = t43 * t168;
t297 = t42 * t165;
t338 = ((t100 * t165 - t102 * t168) * qJD(3) - t296 + t297) * t166 - (t100 * t168 + t102 * t165) * t262;
t337 = 2 * pkin(1);
t288 = t165 * qJ(4);
t318 = pkin(3) + pkin(4);
t198 = -t318 * t168 - t288;
t94 = pkin(2) - t198;
t336 = t43 * qJ(5) + t100 * qJD(5);
t335 = 0.2e1 * t344;
t150 = pkin(6) * t265;
t117 = qJD(2) * pkin(7) + t150;
t155 = t166 * pkin(7);
t158 = t169 * pkin(2);
t246 = -pkin(1) - t158;
t207 = t246 - t155;
t88 = t207 * qJD(1);
t48 = -t165 * t117 + t168 * t88;
t275 = qJD(4) - t48;
t319 = t102 ^ 2;
t332 = -t340 ^ 2 - t319;
t263 = qJD(2) * t166;
t330 = qJ(4) * t263 - t169 * qJD(4);
t329 = t165 * qJD(4) + t150;
t222 = g(1) * t170 + g(2) * t167;
t320 = t100 ^ 2;
t327 = t319 - t320;
t89 = t98 * pkin(3);
t326 = t89 - qJDD(4);
t317 = pkin(7) * t98;
t116 = -qJD(2) * pkin(2) + pkin(6) * t266;
t200 = t102 * qJ(4) - t116;
t40 = pkin(3) * t100 - t200;
t325 = -t340 * t40 + t317;
t259 = qJD(3) * t168;
t261 = qJD(3) * t165;
t223 = pkin(2) * t166 - pkin(7) * t169;
t106 = t223 * qJD(2);
t54 = qJD(1) * t106 + t207 * qJDD(1);
t78 = t328 * pkin(6) + qJDD(2) * pkin(7);
t235 = t117 * t259 + t165 * t78 - t168 * t54 + t88 * t261;
t287 = t165 * t166;
t281 = t167 * t169;
t81 = t165 * t281 + t168 * t170;
t277 = t170 * t165;
t83 = -t167 * t168 + t169 * t277;
t187 = g(1) * t83 + g(2) * t81 + g(3) * t287 - t235;
t184 = t187 + t326;
t24 = -t318 * t100 + qJD(5) + t200;
t298 = t42 * qJ(5);
t324 = (qJD(5) + t24) * t102 + t184 - t298;
t197 = -t165 * t98 - t259 * t340;
t279 = t168 * t169;
t20 = qJD(1) * (-t102 * t166 - t279 * t340) - t197;
t196 = t168 * t98 - t261 * t340;
t286 = t165 * t169;
t323 = qJD(1) * (-t100 * t166 - t286 * t340) - t196;
t8 = (qJD(2) * t102 + t196) * t166 - (-t256 * t340 - t42) * t169;
t176 = (-qJD(2) * t100 + t197) * t166 + (-t264 * t340 + t43) * t169;
t316 = t98 * pkin(4);
t315 = pkin(3) * t165;
t314 = pkin(6) * t165;
t313 = pkin(7) * t102;
t312 = g(1) * t167;
t309 = g(3) * t169;
t308 = pkin(7) - qJ(5);
t251 = t318 * t165;
t293 = qJ(4) * t168;
t199 = -t251 + t293;
t307 = t340 * t199 + t329;
t144 = qJ(4) * t266;
t255 = t168 * qJD(5);
t284 = t166 * t168;
t105 = t223 * qJD(1);
t85 = t165 * t105;
t306 = -t308 * t261 - t255 - t144 - t85 - (-pkin(6) * t284 + qJ(5) * t286) * qJD(1);
t216 = -t293 + t315;
t305 = t340 * t216 - t329;
t111 = t308 * t168;
t245 = -pkin(3) - t314;
t234 = -pkin(4) + t245;
t280 = t168 * t105;
t304 = qJD(3) * t111 - t165 * qJD(5) + t280 - (-qJ(5) * t279 + t234 * t166) * qJD(1);
t270 = t158 + t155;
t108 = -pkin(1) - t270;
t303 = t165 * t106 + t108 * t259;
t302 = qJ(4) * t43;
t120 = t340 * qJ(4);
t49 = t117 * t168 + t165 * t88;
t36 = t120 + t49;
t301 = t340 * t36;
t300 = t340 * t49;
t30 = qJ(5) * t100 + t49;
t23 = t120 + t30;
t299 = t23 * t340;
t131 = pkin(6) * t279;
t295 = qJD(3) * t131 + t108 * t261;
t64 = t165 * t108 + t131;
t294 = pkin(6) * qJDD(1);
t292 = t100 * qJ(4);
t290 = t102 * t100;
t173 = qJD(1) ^ 2;
t282 = t166 * t173;
t278 = t169 * t170;
t29 = t102 * qJ(5) + t48;
t276 = qJD(4) - t29;
t273 = t331 * t165;
t272 = t331 * t168;
t271 = g(1) * t285 - g(2) * t283;
t269 = t170 * pkin(1) + t167 * pkin(6);
t162 = t166 ^ 2;
t163 = t169 ^ 2;
t268 = t162 - t163;
t267 = t162 + t163;
t258 = qJD(4) * t168;
t250 = t169 * t282;
t249 = t24 * t261;
t248 = t24 * t259;
t148 = pkin(6) * t252;
t79 = -qJDD(2) * pkin(2) + pkin(6) * t239 + t148;
t247 = g(1) * t278 + g(2) * t281 + g(3) * t166;
t7 = t43 * pkin(3) + t42 * qJ(4) - t102 * qJD(4) + t79;
t5 = -pkin(4) * t43 + qJDD(5) - t7;
t244 = t5 - t309;
t241 = t340 * t266;
t82 = t167 * t279 - t277;
t237 = -t81 * pkin(3) + qJ(4) * t82;
t84 = t165 * t167 + t168 * t278;
t236 = -t83 * pkin(3) + qJ(4) * t84;
t129 = pkin(6) * t286;
t63 = t108 * t168 - t129;
t233 = -pkin(7) * t296 - t247;
t231 = pkin(3) * t279 + qJ(4) * t286 + t270;
t230 = pkin(2) * t278 + pkin(7) * t283 + t269;
t229 = t166 * t239;
t228 = g(1) * t81 - g(2) * t83;
t227 = g(1) * t82 - g(2) * t84;
t159 = t170 * pkin(6);
t226 = -t82 * pkin(3) - t81 * qJ(4) + t159;
t57 = -qJ(4) * t169 + t64;
t225 = -t168 * t106 + t295;
t224 = t245 * t166;
t221 = -g(2) * t170 + t312;
t218 = (qJD(3) * t100 - t42) * pkin(7);
t217 = pkin(3) * t168 + t288;
t215 = pkin(6) * t100 + t116 * t165;
t214 = pkin(6) * t102 + t116 * t168;
t213 = qJD(3) * t116 - t317;
t33 = -pkin(3) * t340 + t275;
t212 = -t165 * t36 + t168 * t33;
t210 = -t165 * t49 - t168 * t48;
t60 = -pkin(6) * t242 + t85;
t172 = qJD(2) ^ 2;
t204 = qJDD(2) * t169 - t166 * t172;
t9 = t235 - t326;
t203 = pkin(2) + t217;
t202 = -pkin(7) * qJD(3) * t340 - t309;
t55 = -t169 * t98 + t263 * t340;
t10 = -t117 * t261 + t165 * t54 + t168 * t78 + t88 * t259;
t193 = t202 - t7;
t192 = t202 - t79;
t190 = t84 * pkin(3) + qJ(4) * t83 + t230;
t189 = t207 * t312;
t6 = t10 + t344;
t188 = -t98 + t290;
t17 = t165 * t291 - t296;
t27 = (-t166 * t256 - t169 * t261) * pkin(6) + t303;
t183 = t42 - t291;
t13 = t43 * t287 + (t165 * t262 + t166 * t259) * t100;
t181 = t102 * t40 - t184;
t180 = g(1) * t84 + g(2) * t82 + g(3) * t284 - t10;
t178 = t340 * t48 + t180;
t157 = t169 * pkin(3);
t134 = pkin(7) * t278;
t130 = pkin(7) * t281;
t124 = qJ(4) * t284;
t110 = t308 * t165;
t69 = -t124 + (pkin(6) + t315) * t166;
t59 = pkin(6) * t243 + t280;
t58 = t157 - t63;
t56 = t124 + (-pkin(6) - t251) * t166;
t53 = pkin(3) * t102 + t292;
t52 = qJD(1) * t224 - t280;
t51 = t144 + t60;
t47 = qJ(5) * t287 + t57;
t44 = t169 * pkin(4) + t129 + t157 + (-qJ(5) * t166 - t108) * t168;
t31 = -t318 * t102 - t292;
t28 = t263 * t314 - t225;
t26 = (t217 * qJD(3) - t258) * t166 + (pkin(6) + t216) * t262;
t25 = qJD(2) * t224 + t225;
t21 = t27 + t330;
t19 = (qJD(3) * t198 + t258) * t166 + (-pkin(6) + t199) * t262;
t18 = -t318 * t340 + t276;
t16 = t168 * t289 - t297;
t15 = (-pkin(6) * qJD(2) + qJ(5) * qJD(3)) * t284 + (qJD(5) * t166 + (-pkin(6) * qJD(3) + qJ(5) * qJD(2)) * t169) * t165 + t303 + t330;
t14 = (-qJ(5) * t262 - t106) * t168 + (qJ(5) * t261 + t234 * qJD(2) - t255) * t166 + t295;
t12 = -t42 * t284 + (-t165 * t260 + t169 * t256) * t102;
t3 = t6 + t336;
t1 = -qJD(5) * t102 + t298 - t316 + t9;
t2 = [0, 0, 0, 0, 0, qJDD(1), t221, t222, 0, 0, qJDD(1) * t162 + 0.2e1 * t229, 0.2e1 * t166 * t153 - 0.2e1 * t268 * t253, qJDD(2) * t166 + t169 * t172, qJDD(1) * t163 - 0.2e1 * t229, t204, 0, (-0.2e1 * pkin(1) * t253 - pkin(6) * qJDD(2)) * t166 + (-pkin(6) * t172 + qJDD(1) * t337 + t221) * t169, -t204 * pkin(6) + t341 * t337 - t271, 0.2e1 * t267 * t294 - t222, -g(1) * (-t167 * pkin(1) + t159) - g(2) * t269 + (pkin(6) ^ 2 * t267 + (pkin(1) ^ 2)) * qJDD(1), t12, t338, t8, t13, t176, t55, t28 * t340 + t63 * t98 + (qJD(2) * t215 + t235) * t169 + (pkin(6) * t43 + qJD(2) * t48 + t116 * t259 + t79 * t165) * t166 + t227, -t27 * t340 - t64 * t98 + (qJD(2) * t214 + t10) * t169 + (-pkin(6) * t42 - qJD(2) * t49 - t116 * t261 + t79 * t168) * t166 - t228, -t27 * t100 - t28 * t102 + t63 * t42 - t64 * t43 + t210 * t262 + (-t10 * t165 + t235 * t168 + (t165 * t48 - t168 * t49) * qJD(3)) * t166 + t271, t10 * t64 + t49 * t27 - t235 * t63 + t48 * t28 - g(1) * t159 - g(2) * t230 - t189 + (t116 * t262 + t166 * t79) * pkin(6), t12, t8, -t338, t55, -t176, t13, t26 * t100 - t25 * t340 + t69 * t43 - t58 * t98 + (t264 * t40 + t9) * t169 + (-qJD(2) * t33 + t7 * t165 + t259 * t40) * t166 + t227, -t21 * t100 + t25 * t102 - t58 * t42 - t57 * t43 + t212 * t262 + (-t165 * t6 + t168 * t9 + (-t165 * t33 - t168 * t36) * qJD(3)) * t166 + t271, -t26 * t102 + t21 * t340 + t69 * t42 + t57 * t98 + (-t256 * t40 - t6) * t169 + (qJD(2) * t36 - t7 * t168 + t261 * t40) * t166 + t228, -g(1) * t226 - g(2) * t190 + t36 * t21 + t33 * t25 + t40 * t26 + t6 * t57 + t9 * t58 + t7 * t69 - t189, t12, -t338, -t8, t13, t176, t55, -t19 * t100 - t14 * t340 - t56 * t43 - t44 * t98 + (-t24 * t264 + t1) * t169 + (-qJD(2) * t18 - t5 * t165 - t248) * t166 + t227, t19 * t102 + t15 * t340 - t56 * t42 + t47 * t98 + (t24 * t256 - t3) * t169 + (qJD(2) * t23 + t5 * t168 - t249) * t166 + t228, t15 * t100 - t14 * t102 + t44 * t42 + t47 * t43 + (t165 * t23 - t168 * t18) * t262 + (-t1 * t168 + t165 * t3 + (t165 * t18 + t168 * t23) * qJD(3)) * t166 - t271, t3 * t47 + t23 * t15 + t1 * t44 + t18 * t14 + t5 * t56 + t24 * t19 - g(1) * (-t82 * pkin(4) + t226) - g(2) * (pkin(4) * t84 - qJ(5) * t283 + t190) - (-t166 * t308 + t246) * t312; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t250, t268 * t173, t252, t250, t153, qJDD(2), pkin(1) * t282 - t148 - t309 + t331, (pkin(1) * t173 - t294) * t169 + t247, 0, 0, t16, -t4, t20, t17, -t323, -t241, -pkin(2) * t43 - t59 * t340 + t213 * t165 + t192 * t168 + (-t166 * t48 - t169 * t215) * qJD(1) + t272, pkin(2) * t42 + t60 * t340 + t213 * t168 - t192 * t165 + (t166 * t49 - t169 * t214) * qJD(1) - t273, t60 * t100 + t59 * t102 + (t48 * t265 + t10 + (-t48 + t313) * qJD(3)) * t168 + (t235 + t218 - t300) * t165 + t233, -t79 * pkin(2) - t49 * t60 - t48 * t59 - t116 * t150 - g(1) * (-pkin(2) * t283 + t134) - g(2) * (-pkin(2) * t285 + t130) - g(3) * t270 + (qJD(3) * t210 + t10 * t168 + t165 * t235) * pkin(7), t16, t20, t4, -t241, t323, t17, t305 * t100 - t325 * t165 + t193 * t168 - t203 * t43 + t33 * t266 + t340 * t52 + t272, t100 * t51 - t102 * t52 + (-t33 * t265 + t6 + (t33 + t313) * qJD(3)) * t168 + (t218 + t9 - t301) * t165 + t233, -t305 * t102 + t193 * t165 + t325 * t168 - t203 * t42 - t36 * t266 - t340 * t51 + t273, -t36 * t51 - t33 * t52 - g(1) * t134 - g(2) * t130 - g(3) * t231 + t305 * t40 + (qJD(3) * t212 + t9 * t165 + t6 * t168) * pkin(7) + (t222 * t166 - t7) * t203, t16, t4, -t20, t17, -t323, -t241, -t249 - t110 * t98 - t94 * t43 + t244 * t168 - t304 * t340 - t307 * t100 + (t166 * t18 + t24 * t286) * qJD(1) + t272, t248 + t111 * t98 - t94 * t42 + t244 * t165 + t306 * t340 + t307 * t102 + (-t166 * t23 - t24 * t279) * qJD(1) + t273, t110 * t42 + t111 * t43 - t304 * t102 + t306 * t100 + (-t18 * t340 - t3) * t168 + (-t1 + t299) * t165 + t247, t3 * t111 + t1 * t110 + t5 * t94 - g(1) * (-qJ(5) * t278 + t134) - g(2) * (-qJ(5) * t281 + t130) - g(3) * (pkin(4) * t279 + t231) + t307 * t24 + t306 * t23 + t304 * t18 + (g(3) * qJ(5) + t222 * t94) * t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t290, t327, -t183, -t290, -t347, t98, -t102 * t116 + t187 + t300, t100 * t116 + t178, 0, 0, t290, -t183, -t327, t98, t347, -t290, -t100 * t53 - t181 + t300 + t89, pkin(3) * t42 - t302 + (t36 - t49) * t102 + (t33 - t275) * t100, -t100 * t40 + t102 * t53 - t178 + t335, t6 * qJ(4) - t9 * pkin(3) - t40 * t53 - t33 * t49 - g(1) * t236 - g(2) * t237 - g(3) * (-pkin(3) * t287 + t124) + t275 * t36, t290, -t327, t183, -t290, -t347, t98, t100 * t31 + t340 * t30 + (pkin(4) + t318) * t98 + t324, t100 * t24 - t102 * t31 - t29 * t340 - t180 + t335 + t336, t302 - t318 * t42 + (-t23 + t30) * t102 + (-t18 + t276) * t100, t3 * qJ(4) - t1 * t318 - t18 * t30 - t24 * t31 - g(1) * (-pkin(4) * t83 + t236) - g(2) * (-pkin(4) * t81 + t237) - g(3) * (-t166 * t251 + t124) + t276 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t188, -t183, t332, t181 - t301, 0, 0, 0, 0, 0, 0, t188, t332, t183, -t299 - t316 - t324; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t346, -t343, -t319 - t320, -t100 * t23 + t102 * t18 + t244 + t331;];
tau_reg = t2;
