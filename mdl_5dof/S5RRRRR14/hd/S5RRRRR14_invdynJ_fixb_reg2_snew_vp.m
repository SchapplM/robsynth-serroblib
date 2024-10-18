% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRRR14
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
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRRR14_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR14_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR14_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR14_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:43:39
% EndTime: 2024-09-27 18:43:44
% DurationCPUTime: 2.47s
% Computational Cost: add. (45129->362), mult. (69795->542), div. (0->0), fcn. (51307->12), ass. (0->254)
t270 = sin(qJ(5));
t271 = sin(qJ(4));
t276 = cos(qJ(4));
t265 = qJD(1) + qJD(2);
t268 = sin(pkin(5));
t277 = cos(qJ(3));
t325 = t268 * t277;
t307 = t265 * t325;
t272 = sin(qJ(3));
t326 = t268 * t272;
t308 = t265 * t326;
t227 = t271 * t308 - t276 * t307;
t327 = t265 * t268;
t229 = (t277 * t271 + t272 * t276) * t327;
t275 = cos(qJ(5));
t201 = t275 * t227 + t270 * t229;
t203 = -t270 * t227 + t275 * t229;
t156 = t203 * t201;
t263 = qJDD(1) + qJDD(2);
t269 = cos(pkin(5));
t257 = t269 * t263 + qJDD(3);
t250 = qJDD(4) + t257;
t243 = qJDD(5) + t250;
t344 = -t156 + t243;
t348 = t270 * t344;
t207 = t229 * t227;
t342 = -t207 + t250;
t347 = t271 * t342;
t346 = t275 * t344;
t345 = t276 * t342;
t274 = sin(qJ(1));
t279 = cos(qJ(1));
t302 = t274 * g(1) - t279 * g(2);
t247 = qJDD(1) * pkin(1) + t302;
t286 = t279 * g(1) + t274 * g(2);
t249 = -qJD(1) ^ 2 * pkin(1) - t286;
t273 = sin(qJ(2));
t278 = cos(qJ(2));
t221 = t273 * t247 + t278 * t249;
t262 = t265 ^ 2;
t339 = pkin(8) * t268;
t212 = -t262 * pkin(2) + t263 * t339 + t221;
t220 = t278 * t247 - t273 * t249;
t329 = t262 * t268;
t211 = t263 * pkin(2) + pkin(8) * t329 + t220;
t336 = t211 * t269;
t299 = -t272 * t212 + t277 * t336;
t169 = g(3) * t325 - t299;
t170 = -g(3) * t326 + t277 * t212 + t272 * t336;
t129 = t272 * t169 + t277 * t170;
t303 = qJD(3) * t327;
t235 = t263 * t326 + t277 * t303;
t236 = t263 * t325 - t272 * t303;
t298 = t271 * t235 - t276 * t236;
t181 = -t229 * qJD(4) - t298;
t182 = -t227 * qJD(4) + t276 * t235 + t271 * t236;
t127 = -t201 * qJD(5) + t270 * t181 + t275 * t182;
t258 = t269 * t265 + qJD(3);
t253 = qJD(4) + t258;
t248 = qJD(5) + t253;
t186 = t248 * t201;
t343 = -t186 + t127;
t218 = t253 * t227;
t341 = -t218 + t182;
t161 = t218 + t182;
t199 = t201 ^ 2;
t200 = t203 ^ 2;
t225 = t227 ^ 2;
t226 = t229 ^ 2;
t244 = t248 ^ 2;
t251 = t253 ^ 2;
t256 = t258 ^ 2;
t340 = pkin(3) * t269;
t264 = t268 ^ 2;
t328 = t264 * t272;
t311 = t262 * t328;
t331 = t258 * t265;
t139 = t257 * pkin(3) - t235 * pkin(9) + (pkin(3) * t311 + (pkin(9) * t331 - g(3)) * t268) * t277 + t299;
t283 = t258 * pkin(3) - pkin(9) * t308;
t267 = t277 ^ 2;
t330 = t262 * t264;
t309 = t267 * t330;
t140 = -pkin(3) * t309 + t236 * pkin(9) - t258 * t283 + t170;
t94 = -t276 * t139 + t271 * t140;
t76 = t342 * pkin(4) - t161 * pkin(10) - t94;
t287 = t253 * pkin(4) - t229 * pkin(10);
t95 = t271 * t139 + t276 * t140;
t78 = -t225 * pkin(4) + t181 * pkin(10) - t253 * t287 + t95;
t51 = t270 * t78 - t275 * t76;
t52 = t270 * t76 + t275 * t78;
t23 = t270 * t52 - t275 * t51;
t338 = t271 * t23;
t337 = t276 * t23;
t335 = t248 * t270;
t334 = t248 * t275;
t333 = t253 * t271;
t332 = t253 * t276;
t205 = t269 * g(3) + t268 * t211;
t150 = t236 * pkin(3) + pkin(9) * t309 - t283 * t308 + t205;
t110 = t181 * pkin(4) + t225 * pkin(10) - t229 * t287 + t150;
t324 = t270 * t110;
t147 = t156 + t243;
t323 = t270 * t147;
t322 = t271 * t150;
t193 = t207 + t250;
t321 = t271 * t193;
t245 = t277 * t311;
t233 = t245 + t257;
t319 = t272 * t233;
t318 = t275 * t110;
t317 = t275 * t147;
t316 = t276 * t150;
t315 = t276 * t193;
t234 = -t245 + t257;
t313 = t277 * t234;
t112 = t268 * t205 + (-t169 * t277 + t170 * t272) * t269;
t312 = pkin(2) * t112 + t129 * t339;
t266 = t272 ^ 2;
t310 = t266 * t330;
t306 = t269 * t156;
t305 = t269 * t207;
t240 = t258 * t307;
t304 = t240 + t235;
t24 = t270 * t51 + t275 * t52;
t62 = t271 * t94 + t276 * t95;
t11 = t271 * t24 + t337;
t22 = pkin(4) * t23;
t10 = pkin(3) * t11 + t22;
t12 = t276 * t24 - t338;
t18 = pkin(4) * t110 + pkin(10) * t24;
t3 = t268 * t110 + (t11 * t277 + t12 * t272) * t269;
t6 = -t272 * t11 + t277 * t12;
t301 = (-pkin(9) * t11 - pkin(10) * t337 - t271 * t18) * t326 + (pkin(3) * t110 + pkin(9) * t12 - pkin(10) * t338 + t276 * t18) * t325 + t269 * t10 + pkin(2) * t3 + t6 * t339;
t300 = -t275 * t181 + t270 * t182;
t213 = -t240 + t235;
t239 = t258 * t308;
t214 = t236 + t239;
t152 = (-t277 * t213 + t272 * t214) * t269 - (-t266 - t267) * t264 * t329;
t176 = t272 * t213 + t277 * t214;
t297 = pkin(2) * t152 + t129 * t268 + t176 * t339;
t222 = -t310 - t256;
t168 = -t268 * t304 + (t277 * t222 - t272 * t234) * t269;
t198 = -t272 * t222 - t313;
t296 = pkin(2) * t168 - t269 * t170 + t198 * t339 - t205 * t326;
t215 = t236 - t239;
t237 = -t256 - t309;
t174 = t268 * t215 + (t233 * t277 + t237 * t272) * t269;
t209 = t277 * t237 - t319;
t295 = pkin(2) * t174 - t269 * t169 + t205 * t325 + t209 * t339;
t138 = -t199 - t200;
t108 = t186 + t127;
t281 = (-qJD(5) + t248) * t203 - t300;
t68 = t270 * t108 + t275 * t281;
t13 = -pkin(4) * t138 + pkin(10) * t68 + t24;
t66 = -t275 * t108 + t270 * t281;
t14 = -pkin(10) * t66 - t23;
t44 = t271 * t68 + t276 * t66;
t45 = -t271 * t66 + t276 * t68;
t16 = -t268 * t138 + (t272 * t45 + t277 * t44) * t269;
t20 = -t272 * t44 + t277 * t45;
t64 = pkin(4) * t66;
t30 = pkin(3) * t44 + t64;
t294 = (-pkin(9) * t44 - t271 * t13 + t276 * t14) * t326 + (-pkin(3) * t138 + pkin(9) * t45 + t276 * t13 + t271 * t14) * t325 + t269 * t30 + pkin(2) * t16 + t20 * t339;
t104 = (qJD(5) + t248) * t203 + t300;
t145 = -t244 - t199;
t116 = t270 * t145 + t346;
t291 = pkin(4) * t116 - t51;
t117 = t275 * t145 - t348;
t73 = t276 * t116 + t271 * t117;
t28 = pkin(3) * t73 + t291;
t74 = -t271 * t116 + t276 * t117;
t39 = -t268 * t104 + (t272 * t74 + t277 * t73) * t269;
t49 = -t272 * t73 + t277 * t74;
t57 = -pkin(4) * t104 + pkin(10) * t117 + t318;
t71 = -pkin(10) * t116 - t324;
t293 = (-pkin(9) * t73 - t271 * t57 + t276 * t71) * t326 + (-pkin(3) * t104 + pkin(9) * t74 + t271 * t71 + t276 * t57) * t325 + t269 * t28 + pkin(2) * t39 + t49 * t339;
t177 = -t200 - t244;
t130 = t275 * t177 - t323;
t284 = pkin(4) * t130 - t52;
t131 = -t270 * t177 - t317;
t85 = t276 * t130 + t271 * t131;
t33 = pkin(3) * t85 + t284;
t86 = -t271 * t130 + t276 * t131;
t43 = -t268 * t343 + (t272 * t86 + t277 * t85) * t269;
t56 = -t272 * t85 + t277 * t86;
t59 = -pkin(4) * t343 + pkin(10) * t131 - t324;
t77 = -pkin(10) * t130 - t318;
t292 = (-pkin(9) * t85 - t271 * t59 + t276 * t77) * t326 + (-pkin(3) * t343 + pkin(9) * t86 + t271 * t77 + t276 * t59) * t325 + t269 * t33 + pkin(2) * t43 + t56 * t339;
t280 = (-qJD(4) + t253) * t229 - t298;
t123 = -t276 * t161 + t271 * t280;
t124 = t271 * t161 + t276 * t280;
t183 = -t225 - t226;
t61 = t271 * t95 - t276 * t94;
t70 = -t268 * t183 + (t123 * t277 + t124 * t272) * t269;
t82 = -t272 * t123 + t277 * t124;
t290 = (-pkin(9) * t123 - t61) * t326 + (-pkin(3) * t183 + pkin(9) * t124 + t62) * t325 + t123 * t340 + pkin(2) * t70 + t82 * t339;
t190 = -t251 - t225;
t143 = t271 * t190 + t345;
t144 = t276 * t190 - t347;
t114 = -t272 * t143 + t277 * t144;
t157 = (qJD(4) + t253) * t229 + t298;
t80 = pkin(3) * t143 - t94;
t88 = -t268 * t157 + (t143 * t277 + t144 * t272) * t269;
t289 = pkin(2) * t88 + t269 * t80 + t114 * t339 + (-pkin(9) * t143 - t322) * t326 + (-pkin(3) * t157 + pkin(9) * t144 + t316) * t325;
t210 = -t226 - t251;
t153 = t276 * t210 - t321;
t154 = -t271 * t210 - t315;
t121 = -t272 * t153 + t277 * t154;
t84 = pkin(3) * t153 - t95;
t91 = -t268 * t341 + (t153 * t277 + t154 * t272) * t269;
t288 = pkin(2) * t91 + t269 * t84 + (-pkin(9) * t153 - t316) * t326 + t121 * t339 + (-pkin(3) * t341 + pkin(9) * t154 - t322) * t325;
t285 = -t262 * t269 + t331;
t35 = t268 * t150 + (t272 * t62 + t277 * t61) * t269;
t37 = -t272 * t61 + t277 * t62;
t282 = pkin(2) * t35 + t37 * t339 + (pkin(3) * t150 + pkin(9) * t62) * t325 + (-pkin(9) * t326 + t340) * t61;
t242 = t269 * t257;
t238 = (t266 - t267) * t330;
t217 = -t226 + t251;
t216 = t225 - t251;
t206 = t226 - t225;
t189 = (t264 * t277 * t285 + t235 * t268) * t272;
t188 = (t236 * t268 - t285 * t328) * t277;
t185 = -t200 + t244;
t184 = t199 - t244;
t173 = t269 * t214 + (t272 * (-t256 + t309) + t313) * t268;
t172 = t269 * t213 + (t319 + t277 * (t256 - t310)) * t268;
t155 = t200 - t199;
t151 = t269 * t238 + (t272 * t215 + t277 * t304) * t268;
t142 = (-t201 * t275 + t203 * t270) * t248;
t141 = (-t201 * t270 - t203 * t275) * t248;
t136 = t269 * t250 + (t272 * (-t227 * t276 + t229 * t271) + t277 * (-t227 * t271 - t229 * t276)) * t268 * t253;
t135 = t275 * t184 - t323;
t134 = -t270 * t185 + t346;
t133 = t270 * t184 + t317;
t132 = t275 * t185 + t348;
t126 = -t203 * qJD(5) - t300;
t103 = t275 * t127 - t203 * t335;
t102 = t270 * t127 + t203 * t334;
t101 = -t270 * t126 + t201 * t334;
t100 = t275 * t126 + t201 * t335;
t99 = t305 + (t272 * (t276 * t182 - t229 * t333) + t277 * (t271 * t182 + t229 * t332)) * t268;
t98 = -t305 + (t272 * (-t271 * t181 + t227 * t332) + t277 * (t276 * t181 + t227 * t333)) * t268;
t97 = t269 * t280 + (t272 * (t276 * t216 - t321) + t277 * (t271 * t216 + t315)) * t268;
t96 = t269 * t161 + (t272 * (-t271 * t217 + t345) + t277 * (t276 * t217 + t347)) * t268;
t72 = t269 * t206 + (t272 * (-t276 * t157 - t271 * t341) + t277 * (-t271 * t157 + t276 * t341)) * t268;
t67 = -t275 * t104 - t270 * t343;
t65 = -t270 * t104 + t275 * t343;
t63 = t269 * t243 + (t272 * (-t271 * t141 + t276 * t142) + t277 * (t276 * t141 + t271 * t142)) * t268;
t47 = t269 * t281 + (t272 * (-t271 * t133 + t276 * t135) + t277 * (t276 * t133 + t271 * t135)) * t268;
t46 = t269 * t108 + (t272 * (-t271 * t132 + t276 * t134) + t277 * (t276 * t132 + t271 * t134)) * t268;
t41 = t306 + (t272 * (-t271 * t102 + t276 * t103) + t277 * (t276 * t102 + t271 * t103)) * t268;
t40 = -t306 + (t272 * (-t271 * t100 + t276 * t101) + t277 * (t276 * t100 + t271 * t101)) * t268;
t17 = t269 * t155 + (t272 * (-t271 * t65 + t276 * t67) + t277 * (t271 * t67 + t276 * t65)) * t268;
t1 = [0, 0, 0, 0, 0, qJDD(1), t302, t286, 0, 0, 0, 0, 0, 0, 0, t263, pkin(1) * (-t273 * t262 + t278 * t263) + t220, pkin(1) * (-t278 * t262 - t273 * t263) - t221, 0, pkin(1) * (t278 * t220 + t273 * t221), t189, t151, t172, t188, t173, t242, pkin(1) * (t278 * t174 + t273 * t209) + t295, pkin(1) * (t278 * t168 + t273 * t198) + t296, pkin(1) * (t278 * t152 + t273 * t176) + t297, pkin(1) * (t278 * t112 + t273 * t129) + t312, t99, t72, t96, t98, t97, t136, pkin(1) * (t273 * t114 + t278 * t88) + t289, pkin(1) * (t273 * t121 + t278 * t91) + t288, pkin(1) * (t273 * t82 + t278 * t70) + t290, pkin(1) * (t273 * t37 + t278 * t35) + t282, t41, t17, t46, t40, t47, t63, pkin(1) * (t273 * t49 + t278 * t39) + t293, pkin(1) * (t273 * t56 + t278 * t43) + t292, pkin(1) * (t278 * t16 + t273 * t20) + t294, pkin(1) * (t273 * t6 + t278 * t3) + t301; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t263, t220, -t221, 0, 0, t189, t151, t172, t188, t173, t242, t295, t296, t297, t312, t99, t72, t96, t98, t97, t136, t289, t288, t290, t282, t41, t17, t46, t40, t47, t63, t293, t292, t294, t301; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t245, t238, t213, t245, t214, t257, -t169, -t170, 0, 0, t207, t206, t161, -t207, t280, t250, t80, t84, pkin(3) * t123, pkin(3) * t61, t156, t155, t108, -t156, t281, t243, t28, t33, t30, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t207, t206, t161, -t207, t280, t250, -t94, -t95, 0, 0, t156, t155, t108, -t156, t281, t243, t291, t284, t64, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t156, t155, t108, -t156, t281, t243, -t51, -t52, 0, 0;];
tauJ_reg = t1;
