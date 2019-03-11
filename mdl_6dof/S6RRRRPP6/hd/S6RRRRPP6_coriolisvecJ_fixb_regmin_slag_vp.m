% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRRPP6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:15:09
% EndTime: 2019-03-09 21:15:21
% DurationCPUTime: 5.13s
% Computational Cost: add. (7658->498), mult. (18375->613), div. (0->0), fcn. (12390->6), ass. (0->235)
t199 = sin(qJ(2));
t201 = cos(qJ(2));
t230 = pkin(2) * t199 - pkin(8) * t201;
t140 = t230 * qJD(1);
t198 = sin(qJ(3));
t124 = t198 * t140;
t325 = -pkin(9) - pkin(8);
t258 = qJD(3) * t325;
t200 = cos(qJ(3));
t294 = t199 * t200;
t347 = t198 * t258 - t124 - (-pkin(9) * t198 * t201 - pkin(7) * t294) * qJD(1);
t293 = t200 * t201;
t225 = pkin(3) * t199 - pkin(9) * t293;
t276 = qJD(1) * t199;
t251 = t198 * t276;
t279 = pkin(7) * t251 + t200 * t140;
t346 = -qJD(1) * t225 + t200 * t258 - t279;
t270 = qJD(2) * t200;
t136 = -t251 + t270;
t272 = qJD(2) * t198;
t137 = t200 * t276 + t272;
t197 = sin(qJ(4));
t322 = cos(qJ(4));
t90 = -t322 * t136 + t137 * t197;
t326 = t90 ^ 2;
t220 = t197 * t136 + t322 * t137;
t88 = t220 ^ 2;
t345 = t88 - t326;
t275 = qJD(1) * t201;
t256 = t198 * t275;
t257 = t322 * t200;
t297 = t197 * t198;
t327 = qJD(3) + qJD(4);
t248 = t322 * qJD(4);
t329 = t322 * qJD(3) + t248;
t285 = -t197 * t256 - t329 * t200 + t257 * t275 + t327 * t297;
t307 = qJD(2) * pkin(2);
t153 = pkin(7) * t276 - t307;
t110 = -t136 * pkin(3) + t153;
t213 = -qJ(5) * t220 + t110;
t316 = pkin(4) + qJ(6);
t21 = t316 * t90 + t213;
t344 = t21 * t90;
t36 = pkin(4) * t90 + t213;
t343 = t36 * t90;
t318 = t36 * t220;
t263 = qJD(1) * qJD(2);
t183 = t199 * t263;
t177 = pkin(4) * t183;
t265 = qJD(4) * t197;
t216 = t225 * qJD(2);
t143 = t230 * qJD(2);
t130 = qJD(1) * t143;
t237 = pkin(7) * t183;
t282 = t200 * t130 + t198 * t237;
t189 = pkin(7) * t275;
t154 = qJD(2) * pkin(8) + t189;
t147 = -pkin(2) * t201 - pkin(8) * t199 - pkin(1);
t129 = t147 * qJD(1);
t296 = t198 * t129;
t99 = t200 * t154 + t296;
t76 = pkin(9) * t136 + t99;
t35 = qJD(1) * t216 - t76 * qJD(3) + t282;
t268 = qJD(3) * t198;
t266 = qJD(3) * t200;
t283 = t129 * t266 + t198 * t130;
t208 = -t154 * t268 - t200 * t237 + t283;
t267 = qJD(3) * t199;
t246 = qJD(1) * t267;
t247 = t201 * t263;
t262 = qJD(2) * qJD(3);
t259 = t200 * t246 + (t247 + t262) * t198;
t45 = -pkin(9) * t259 + t208;
t179 = -qJD(3) + t275;
t98 = t200 * t129 - t154 * t198;
t75 = -pkin(9) * t137 + t98;
t69 = -pkin(3) * t179 + t75;
t243 = t197 * t45 + t76 * t248 + t69 * t265 - t322 * t35;
t5 = -t177 + t243;
t342 = t318 + t5;
t317 = t220 * t90;
t304 = t90 * qJ(5);
t155 = t325 * t198;
t156 = t325 * t200;
t341 = t155 * t248 + t156 * t265 + t346 * t197 + t347 * t322;
t139 = t197 * t200 + t322 * t198;
t103 = t327 * t139;
t284 = -t139 * t275 + t103;
t232 = -t189 + (-t256 + t268) * pkin(3);
t269 = qJD(2) * t201;
t255 = t198 * t269;
t340 = t199 * t266 + t255;
t305 = t197 * t76;
t26 = -t322 * t69 + t305;
t336 = pkin(5) * t220;
t226 = t26 + t336;
t287 = qJD(5) + t226;
t328 = pkin(5) * t90 - qJD(6);
t161 = -qJD(4) + t179;
t254 = t198 * t267;
t215 = t200 * t269 - t254;
t207 = qJD(1) * t215 + t200 * t262;
t40 = -t136 * t248 + t137 * t265 + t197 * t259 - t322 * t207;
t19 = -t161 * t90 - t40;
t244 = -t197 * t35 - t69 * t248 + t76 * t265 - t322 * t45;
t339 = t110 * t90 + t244;
t338 = -0.2e1 * t263;
t337 = pkin(4) * t220;
t311 = qJ(5) * t276 - t341;
t31 = t322 * t75 - t305;
t301 = -pkin(3) * t248 - qJD(5) + t31;
t73 = t322 * t76;
t30 = t197 * t75 + t73;
t235 = pkin(3) * t265 - t30;
t109 = t197 * t155 - t322 * t156;
t335 = t109 * qJD(4) + t347 * t197 - t346 * t322;
t157 = t161 ^ 2;
t334 = -t88 - t157;
t333 = t110 * t220;
t332 = t220 * t316;
t146 = qJD(5) * t161;
t172 = qJ(5) * t183;
t331 = t172 - t146;
t330 = qJ(5) * t285 - t139 * qJD(5) + t232;
t194 = t199 ^ 2;
t264 = t194 * qJD(1);
t1 = -t40 * pkin(5) - qJ(6) * t183 + qJD(6) * t161 + t5;
t211 = -t21 * t220 - t1;
t41 = qJD(4) * t220 + t197 * t207 + t322 * t259;
t204 = -t161 * t220 - t41;
t323 = t41 * pkin(5);
t321 = pkin(7) * t198;
t320 = pkin(8) * t179;
t138 = -t257 + t297;
t315 = t138 * qJD(6) + t284 * t316 + t330;
t250 = t199 * t316;
t314 = -t285 * pkin(5) + qJD(1) * t250 + t335;
t313 = t284 * pkin(5) + t311;
t312 = t284 * pkin(4) + t330;
t310 = -pkin(4) * t276 - t335;
t27 = t197 * t69 + t73;
t308 = qJ(5) * t41;
t182 = pkin(3) * t197 + qJ(5);
t306 = t182 * t41;
t181 = pkin(7) * t293;
t278 = t198 * t147 + t181;
t295 = t198 * t199;
t104 = -pkin(9) * t295 + t278;
t135 = t200 * t147;
t97 = -pkin(9) * t294 + t135 + (-pkin(3) - t321) * t201;
t303 = t322 * t104 + t197 * t97;
t302 = t235 + t328;
t300 = -t301 + t336;
t299 = t137 * t179;
t298 = t153 * t198;
t292 = t201 * t179;
t203 = qJD(1) ^ 2;
t291 = t201 * t203;
t202 = qJD(2) ^ 2;
t290 = t202 * t199;
t289 = t202 * t201;
t288 = -qJD(5) - t26;
t286 = -t27 + t328;
t281 = t198 * t143 + t147 * t266;
t271 = qJD(2) * t199;
t280 = t200 * t143 + t271 * t321;
t144 = pkin(3) * t295 + t199 * pkin(7);
t277 = -t201 ^ 2 + t194;
t108 = -t322 * t155 - t197 * t156;
t274 = qJD(2) * t108;
t273 = qJD(2) * t109;
t260 = t197 * t295;
t111 = t340 * pkin(3) + pkin(7) * t269;
t187 = -t200 * pkin(3) - pkin(2);
t252 = t179 * t266;
t245 = -t183 + t317;
t242 = pkin(1) * t338;
t241 = t179 + t275;
t240 = -t136 + t270;
t239 = -t137 + t272;
t238 = qJD(3) + t275;
t186 = -t322 * pkin(3) - pkin(4);
t236 = t322 * t269;
t24 = qJ(5) * t161 - t27;
t233 = -t146 - t244;
t52 = qJ(5) * t201 - t303;
t231 = -t197 * t104 + t322 * t97;
t119 = t199 * t257 - t260;
t228 = -t119 * qJ(5) + t144;
t227 = t137 * pkin(3) + t304;
t4 = -t172 - t233;
t53 = t201 * pkin(4) - t231;
t224 = -t139 * qJ(5) + t187;
t223 = -t244 - t323;
t87 = pkin(3) * t259 + pkin(7) * t247;
t222 = -t161 * t27 - t243;
t56 = t216 + (-t181 + (pkin(9) * t199 - t147) * t198) * qJD(3) + t280;
t59 = -t340 * pkin(9) + (-t199 * t270 - t201 * t268) * pkin(7) + t281;
t221 = t104 * t248 + t197 * t59 + t97 * t265 - t322 * t56;
t218 = -t104 * t265 + t197 * t56 + t97 * t248 + t322 * t59;
t217 = t238 * t272;
t214 = t223 - t344;
t60 = t103 * t199 + t197 * t255 - t200 * t236;
t212 = t60 * qJ(5) - t119 * qJD(5) + t111;
t210 = t182 * t183 - t4;
t209 = t40 * qJ(5) - qJD(5) * t220 + t87;
t10 = -qJ(5) * t271 + qJD(5) * t201 - t218;
t176 = -qJ(6) + t186;
t168 = 0.2e1 * t172;
t118 = t139 * t199;
t86 = pkin(4) * t138 + t224;
t79 = -t138 * pkin(5) + t109;
t78 = t139 * pkin(5) + t108;
t74 = t316 * t138 + t224;
t70 = pkin(4) * t118 + t228;
t61 = t198 * t236 - t197 * t254 - qJD(4) * t260 + (t197 * t269 + t329 * t199) * t200;
t55 = t304 + t337;
t54 = t316 * t118 + t228;
t47 = t227 + t337;
t42 = -pkin(5) * t118 - t52;
t34 = t119 * pkin(5) + t201 * qJ(6) + t53;
t29 = t304 + t332;
t23 = pkin(4) * t161 - t288;
t22 = t227 + t332;
t14 = -t24 - t328;
t13 = t316 * t161 + t287;
t12 = pkin(4) * t61 + t212;
t11 = -pkin(4) * t271 + t221;
t9 = t118 * qJD(6) + t316 * t61 + t212;
t8 = t41 * pkin(4) + t209;
t7 = -pkin(5) * t61 - t10;
t6 = -t60 * pkin(5) - qJD(2) * t250 + t201 * qJD(6) + t221;
t3 = t90 * qJD(6) + t316 * t41 + t209;
t2 = t223 + t331;
t15 = [0, 0, 0, 0.2e1 * t201 * t183, t277 * t338, t289, -t290, 0, -pkin(7) * t289 + t199 * t242, pkin(7) * t290 + t201 * t242, t137 * t215 + t207 * t294 (t200 * t136 - t137 * t198) * t269 + ((-t136 + t251) * t268 + (-t137 * qJD(3) - t217 - t259) * t200) * t199, t241 * t254 + (t137 * t199 + (t264 + (-t179 - t238) * t201) * t200) * qJD(2), t199 * t252 + t259 * t201 + (t136 * t199 + (-t264 + t292) * t198) * qJD(2), -t241 * t271 -(-t147 * t268 + t280) * t179 + (pkin(7) * t259 + t153 * t266 + (qJD(1) * t135 + t98) * qJD(2)) * t199 + ((-pkin(7) * t136 + t298) * qJD(2) + (t296 + (pkin(7) * t179 + t154) * t200) * qJD(3) - t282) * t201, t281 * t179 + t283 * t201 + (-t153 * t199 - t154 * t201 + (-t292 - t264) * pkin(7)) * t268 + ((pkin(7) * t137 + t153 * t200) * t201 + (-t278 * qJD(1) - t99 + (-t179 + t238) * pkin(7) * t200) * t199) * qJD(2), -t119 * t40 - t220 * t60, t118 * t40 - t119 * t41 - t220 * t61 + t60 * t90, t60 * t161 + t40 * t201 + (qJD(1) * t119 + t220) * t271, t61 * t161 + t41 * t201 + (-qJD(1) * t118 - t90) * t271 (-t161 - t275) * t271, t221 * t161 + t243 * t201 + t111 * t90 + t144 * t41 + t87 * t118 + t110 * t61 + (qJD(1) * t231 - t26) * t271, t218 * t161 - t244 * t201 + t111 * t220 - t144 * t40 + t87 * t119 - t110 * t60 + (-qJD(1) * t303 - t27) * t271, t10 * t90 + t11 * t220 + t118 * t4 + t119 * t5 - t23 * t60 + t24 * t61 - t40 * t53 + t41 * t52, -t11 * t161 - t8 * t118 - t12 * t90 - t5 * t201 - t36 * t61 - t70 * t41 + (qJD(1) * t53 + t23) * t271, t10 * t161 - t8 * t119 - t12 * t220 + t4 * t201 + t36 * t60 + t70 * t40 + (-qJD(1) * t52 - t24) * t271, t10 * t24 + t11 * t23 + t12 * t36 + t4 * t52 + t5 * t53 + t70 * t8, t1 * t119 - t118 * t2 - t13 * t60 - t14 * t61 + t220 * t6 - t34 * t40 - t41 * t42 - t7 * t90, -t3 * t119 - t7 * t161 - t2 * t201 + t21 * t60 + t54 * t40 - t9 * t220 + (qJD(1) * t42 + t14) * t271, t1 * t201 + t3 * t118 + t6 * t161 + t21 * t61 + t54 * t41 + t9 * t90 + (-qJD(1) * t34 - t13) * t271, t1 * t34 + t13 * t6 + t14 * t7 + t2 * t42 + t21 * t9 + t3 * t54; 0, 0, 0, -t199 * t291, t277 * t203, 0, 0, 0, t203 * pkin(1) * t199, pkin(1) * t291, -t198 ^ 2 * t246 + (t217 - t299) * t200 (-t259 + t299) * t198 + ((t136 + t270) * qJD(3) + (t201 * t240 - t254) * qJD(1)) * t200, -t252 + (t199 * t239 + t200 * t292) * qJD(1), t179 * t268 + (-t198 * t292 + t199 * t240) * qJD(1), t179 * t276, -pkin(2) * t259 + t279 * t179 + (t200 * t320 + t298) * qJD(3) + ((-pkin(8) * t272 - t98) * t199 + (-pkin(7) * t240 - t298) * t201) * qJD(1), -t124 * t179 + (-t198 * t320 + (t153 - t307) * t200) * qJD(3) + ((-t153 - t307) * t293 + (pkin(2) * t268 - pkin(8) * t270 + t99) * t199 + (t179 * t294 + t201 * t239) * pkin(7)) * qJD(1), -t40 * t139 - t220 * t285, t40 * t138 - t139 * t41 - t220 * t284 + t285 * t90, t285 * t161 + (qJD(2) * t139 - t220) * t276, t284 * t161 + (-qJD(2) * t138 + t90) * t276, t161 * t276, t87 * t138 + t187 * t41 + t232 * t90 + t335 * t161 + t284 * t110 + (t26 - t274) * t276, t87 * t139 - t187 * t40 + t232 * t220 + t341 * t161 - t285 * t110 + (t27 - t273) * t276, -t108 * t40 - t109 * t41 + t4 * t138 + t5 * t139 - t220 * t310 - t285 * t23 + t284 * t24 + t311 * t90, -t8 * t138 - t86 * t41 - t312 * t90 - t284 * t36 + t310 * t161 + (-t23 + t274) * t276, -t8 * t139 + t86 * t40 - t312 * t220 + t285 * t36 + t311 * t161 + (t24 + t273) * t276, t5 * t108 - t4 * t109 - t310 * t23 + t311 * t24 + t312 * t36 + t8 * t86, t1 * t139 - t285 * t13 - t2 * t138 - t284 * t14 + t220 * t314 + t313 * t90 - t78 * t40 - t79 * t41, -t3 * t139 + t74 * t40 - t315 * t220 + t285 * t21 + t313 * t161 + (qJD(2) * t79 - t14) * t276, t3 * t138 + t74 * t41 + t315 * t90 + t284 * t21 + t314 * t161 + (-qJD(2) * t78 + t13) * t276, t1 * t78 + t314 * t13 - t313 * t14 + t2 * t79 + t315 * t21 + t3 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t137 * t136, -t136 ^ 2 + t137 ^ 2, t136 * t179 + t207, -t259 - t299, t183, -t153 * t137 + t282 + (-qJD(3) - t179) * t99, -t136 * t153 - t179 * t98 - t208, t317, t345, t19, t204, t183, -t333 - t30 * t161 + (-t137 * t90 + t161 * t265 + t322 * t183) * pkin(3) - t243, -t31 * t161 + (-t137 * t220 + t161 * t248 - t183 * t197) * pkin(3) + t339, -t186 * t40 - t306 + (t235 - t24) * t220 + (t23 + t301) * t90, -t161 * t235 + t183 * t186 + t47 * t90 + t342, t161 * t301 + t220 * t47 + t210 - t343, -t4 * t182 + t5 * t186 + t23 * t235 + t24 * t301 - t36 * t47, -t176 * t40 - t306 + (t14 + t302) * t220 + (t13 - t300) * t90, -t161 * t300 + t22 * t220 + t210 - t323 - t344, t161 * t302 - t176 * t183 - t22 * t90 + t211, t1 * t176 + t13 * t302 + t14 * t300 + t2 * t182 - t21 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t317, t345, t19, t204, t183, t222 - t333, t161 * t26 + t339, pkin(4) * t40 - t308 + (-t24 - t27) * t220 + (t23 + t288) * t90, t55 * t90 - 0.2e1 * t177 - t222 + t318, t161 * t288 + t220 * t55 + t168 + t233 - t343, -t5 * pkin(4) - t4 * qJ(5) - t23 * t27 + t24 * t288 - t36 * t55, -t308 + t316 * t40 + (t14 + t286) * t220 + (t13 - t287) * t90, -t161 * t226 + t220 * t29 - 0.2e1 * t146 + t168 + t214, t161 * t286 + t183 * t316 - t29 * t90 + t211, t2 * qJ(5) - t1 * t316 + t13 * t286 + t14 * t287 - t21 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t245, t334, -t161 * t24 + t342, t19, t334, t245, t14 * t161 - t211; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t204, t183 + t317, -t157 - t326, -t13 * t161 + t214 + t331;];
tauc_reg  = t15;
