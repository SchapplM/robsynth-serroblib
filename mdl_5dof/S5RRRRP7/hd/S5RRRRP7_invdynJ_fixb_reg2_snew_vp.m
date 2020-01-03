% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRRP7
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRRP7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP7_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:57:45
% EndTime: 2019-12-31 21:57:59
% DurationCPUTime: 5.63s
% Computational Cost: add. (13648->384), mult. (27917->485), div. (0->0), fcn. (19544->8), ass. (0->246)
t232 = sin(qJ(3));
t233 = sin(qJ(2));
t236 = cos(qJ(3));
t237 = cos(qJ(2));
t207 = (t237 * t232 + t233 * t236) * qJD(1);
t228 = qJD(2) + qJD(3);
t231 = sin(qJ(4));
t235 = cos(qJ(4));
t192 = t235 * t207 + t231 * t228;
t290 = qJD(1) * t233;
t205 = -t236 * t237 * qJD(1) + t232 * t290;
t202 = qJD(4) + t205;
t174 = t202 * t192;
t223 = t233 * qJDD(1);
t287 = qJD(1) * qJD(2);
t277 = t237 * t287;
t212 = t223 + t277;
t224 = t237 * qJDD(1);
t278 = t233 * t287;
t213 = t224 - t278;
t262 = t236 * t212 + t232 * t213;
t165 = -t205 * qJD(3) + t262;
t286 = qJDD(2) + qJDD(3);
t273 = -t231 * t165 + t235 * t286;
t256 = t192 * qJD(4) - t273;
t108 = -t174 + t256;
t329 = t202 ^ 2;
t190 = t231 * t207 - t235 * t228;
t330 = t190 ^ 2;
t170 = t330 - t329;
t271 = t232 * t212 - t236 * t213;
t164 = -t207 * qJD(3) - t271;
t162 = qJDD(4) - t164;
t318 = t192 * t190;
t123 = -t318 - t162;
t310 = t231 * t123;
t91 = -t235 * t170 - t310;
t390 = t233 * (t232 * t108 + t236 * t91) + t237 * (-t236 * t108 + t232 * t91);
t188 = t192 ^ 2;
t341 = -t188 - t329;
t81 = t235 * t341 + t310;
t389 = pkin(1) * t81;
t388 = pkin(2) * t81;
t387 = pkin(3) * t81;
t386 = pkin(8) * t81;
t298 = t235 * t123;
t83 = -t231 * t341 + t298;
t385 = pkin(8) * t83;
t384 = t232 * t83;
t383 = t236 * t83;
t340 = t188 - t330;
t252 = -t235 * t165 - t231 * t286;
t247 = -t190 * qJD(4) - t252;
t319 = t190 * t202;
t346 = t319 - t247;
t311 = t231 * t346;
t343 = t174 + t256;
t58 = t235 * t343 - t311;
t380 = t233 * (-t232 * t340 + t236 * t58) + t237 * (t232 * t58 + t236 * t340);
t337 = -t318 + t162;
t309 = t231 * t337;
t334 = -t329 - t330;
t344 = t235 * t334 - t309;
t362 = t232 * t344 - t236 * t343;
t377 = pkin(2) * t362;
t376 = pkin(7) * t362;
t375 = t346 * qJ(5);
t87 = -t231 * t170 + t298;
t297 = t235 * t337;
t345 = t231 * t334 + t297;
t361 = t232 * t343 + t236 * t344;
t374 = -pkin(2) * t345 + pkin(7) * t361;
t373 = pkin(6) * (-t233 * t362 + t237 * t361) - pkin(1) * t345;
t335 = t319 + t247;
t171 = -t188 + t329;
t363 = -t231 * t171 + t297;
t372 = t233 * (t232 * t335 + t236 * t363) + t237 * (t232 * t363 - t236 * t335);
t369 = pkin(3) * t345;
t368 = pkin(8) * t344;
t367 = pkin(8) * t345;
t364 = t235 * t171 + t309;
t299 = t235 * t346;
t56 = -t231 * t343 - t299;
t339 = t188 + t330;
t360 = pkin(3) * t339;
t359 = -qJ(5) * t231 - pkin(3);
t356 = t232 * t339;
t184 = t207 * t205;
t342 = -t184 + t286;
t354 = t232 * t342;
t350 = t236 * t339;
t348 = t236 * t342;
t199 = t228 * t205;
t145 = t165 - t199;
t230 = t237 ^ 2;
t239 = qJD(1) ^ 2;
t234 = sin(qJ(1));
t327 = cos(qJ(1));
t276 = t234 * g(1) - t327 * g(2);
t259 = qJDD(1) * pkin(1) + t276;
t260 = qJD(2) * pkin(2) - pkin(7) * t290;
t167 = t213 * pkin(2) + (pkin(7) * t230 + pkin(6)) * t239 - t260 * t290 + t259;
t66 = -t145 * pkin(8) + (t228 * t207 - t164) * pkin(3) - t167;
t261 = t327 * g(1) + t234 * g(2);
t320 = qJDD(1) * pkin(6);
t246 = -t239 * pkin(1) - t261 + t320;
t195 = -t233 * g(3) + t237 * t246;
t226 = t230 * t239;
t155 = -pkin(2) * t226 + t213 * pkin(7) - qJD(2) * t260 + t195;
t295 = t236 * t155;
t244 = t233 * t246;
t302 = t233 * t239;
t325 = t212 * pkin(7);
t338 = qJDD(2) * pkin(2) - t244 + (pkin(2) * t302 + pkin(7) * t287 - g(3)) * t237 - t325;
t118 = t232 * t338 + t295;
t181 = t205 * pkin(3) - t207 * pkin(8);
t328 = t228 ^ 2;
t78 = -t328 * pkin(3) + t286 * pkin(8) - t205 * t181 + t118;
t36 = t231 * t78 - t235 * t66;
t37 = t231 * t66 + t235 * t78;
t18 = t231 * t36 + t235 * t37;
t152 = t190 * pkin(4) - t192 * qJ(5);
t270 = t162 * qJ(5) - t190 * t152 + t37;
t333 = -(t341 + t329) * pkin(4) - qJ(5) * t123 + t270;
t316 = t202 * t235;
t282 = t190 * t316;
t257 = t231 * t256 + t282;
t281 = t236 * t318;
t283 = t232 * t318;
t332 = t233 * (t236 * t257 - t283) + t237 * (t232 * t257 + t281);
t317 = t202 * t231;
t168 = t192 * t317;
t267 = t168 - t282;
t331 = t233 * (t232 * t162 + t236 * t267) + t237 * (-t236 * t162 + t232 * t267);
t203 = t205 ^ 2;
t204 = t207 ^ 2;
t326 = pkin(3) * t232;
t117 = t232 * t155 - t236 * t338;
t77 = -t286 * pkin(3) - t328 * pkin(8) + t207 * t181 + t117;
t324 = -pkin(3) * t77 + pkin(8) * t18;
t73 = t231 * t77;
t62 = -t236 * t117 + t232 * t118;
t323 = t233 * t62;
t74 = t235 * t77;
t321 = qJ(5) * t235;
t315 = t228 * t232;
t314 = t228 * t236;
t312 = t231 * t335;
t305 = t232 * t167;
t179 = t184 + t286;
t304 = t232 * t179;
t218 = t237 * t302;
t215 = qJDD(2) + t218;
t303 = t233 * t215;
t300 = t235 * t335;
t294 = t236 * t167;
t293 = t236 * t179;
t292 = t237 * (qJDD(2) - t218);
t289 = qJD(5) * t202;
t288 = qJD(3) + t228;
t115 = (qJD(4) + t202) * t190 + t252;
t285 = pkin(3) * t115 + t385 + t73;
t284 = -pkin(3) * t343 + t368 - t74;
t280 = -pkin(3) * t236 - pkin(2);
t198 = 0.2e1 * t289;
t266 = t198 + t270;
t20 = (t339 - t329) * pkin(4) + t266;
t31 = -t162 * pkin(4) - qJ(5) * t329 + t192 * t152 + qJDD(5) + t36;
t22 = qJ(5) * t339 + t31;
t59 = -t235 * t108 + t312;
t275 = pkin(8) * t59 + t235 * t20 + t231 * t22 + t360;
t110 = (-qJD(4) + t202) * t192 + t273;
t61 = t235 * t110 + t312;
t274 = pkin(8) * t61 + t18 + t360;
t63 = t232 * t117 + t236 * t118;
t194 = t237 * g(3) + t244;
t272 = t233 * t194 + t237 * t195;
t29 = -pkin(4) * t329 + t266;
t269 = -pkin(4) * t31 + qJ(5) * t29;
t268 = t190 * t317 - t235 * t256;
t265 = -pkin(4) * t335 - qJ(5) * t108;
t263 = t231 * t37 - t235 * t36;
t245 = t256 * pkin(4) + t375 + t77;
t243 = 0.2e1 * qJD(5) * t192 - t245;
t25 = -pkin(4) * t174 + t243 - t375;
t258 = -pkin(3) * t346 - pkin(4) * t299 + t231 * t25 - t385;
t26 = (-t343 - t174) * pkin(4) + t243;
t255 = t235 * t26 + t359 * t343 + t368;
t253 = (-t190 * t231 - t192 * t235) * t202;
t251 = (-qJD(3) + t228) * t207 - t271;
t33 = (pkin(4) * t202 - 0.2e1 * qJD(5)) * t192 + t245;
t8 = t231 * t31 + t235 * t29;
t250 = pkin(8) * t8 + (-pkin(4) * t235 + t359) * t33;
t102 = t235 * t247 - t168;
t242 = t233 * (t236 * t102 + t283) + t237 * (t232 * t102 - t281);
t241 = pkin(4) * t337 + qJ(5) * t334 - t31;
t238 = qJD(2) ^ 2;
t229 = t233 ^ 2;
t225 = t229 * t239;
t214 = t224 - 0.2e1 * t278;
t211 = t223 + 0.2e1 * t277;
t208 = t239 * pkin(6) + t259;
t197 = -t204 + t328;
t196 = t203 - t328;
t193 = -t204 - t328;
t183 = t204 - t203;
t177 = -t328 - t203;
t166 = -t203 - t204;
t148 = -t232 * t193 - t293;
t147 = t236 * t193 - t304;
t146 = t165 + t199;
t144 = -t288 * t205 + t262;
t141 = t288 * t207 + t271;
t138 = t236 * t177 - t354;
t137 = t232 * t177 + t348;
t101 = t192 * t316 + t231 * t247;
t94 = t232 * t146 + t236 * t251;
t93 = -t236 * t146 + t232 * t251;
t57 = t231 * t110 - t300;
t55 = -t231 * t108 - t300;
t50 = -t232 * t115 + t383;
t48 = t236 * t115 + t384;
t46 = t232 * t346 - t383;
t44 = -t236 * t346 - t384;
t43 = t74 - t386;
t42 = t236 * t61 - t356;
t41 = t236 * t59 - t356;
t40 = t232 * t61 + t350;
t39 = t232 * t59 + t350;
t38 = t73 - t367;
t30 = -pkin(3) * t55 - t265;
t28 = t37 - t387;
t27 = t36 - t369;
t15 = -t241 - t369;
t14 = -0.2e1 * t289 - t333 + t387;
t13 = -t231 * t26 - t321 * t343 - t367;
t12 = pkin(4) * t311 + t235 * t25 + t386;
t10 = t232 * t18 - t236 * t77;
t9 = -pkin(8) * t57 - t263;
t7 = t231 * t29 - t235 * t31;
t5 = -pkin(8) * t55 - t231 * t20 + t235 * t22;
t4 = t232 * t33 + t236 * t8;
t3 = t232 * t8 - t236 * t33;
t2 = -pkin(8) * t7 + (pkin(4) * t231 - t321) * t33;
t1 = -pkin(3) * t7 - t269;
t6 = [0, 0, 0, 0, 0, qJDD(1), t276, t261, 0, 0, (t212 + t277) * t233, t237 * t211 + t233 * t214, t303 + t237 * (-t225 + t238), (t213 - t278) * t237, t233 * (t226 - t238) + t292, 0, t237 * t208 + pkin(1) * t214 + pkin(6) * (t237 * (-t226 - t238) - t303), -t233 * t208 - pkin(1) * t211 + pkin(6) * (-t292 - t233 * (-t225 - t238)), pkin(1) * (t225 + t226) + (t229 + t230) * t320 + t272, pkin(1) * t208 + pkin(6) * t272, t233 * (t236 * t165 - t207 * t315) + t237 * (t232 * t165 + t207 * t314), t233 * (-t236 * t141 - t232 * t145) + t237 * (-t232 * t141 + t236 * t145), t233 * (-t232 * t197 + t348) + t237 * (t236 * t197 + t354), t233 * (-t232 * t164 + t205 * t314) + t237 * (t236 * t164 + t205 * t315), t233 * (t236 * t196 - t304) + t237 * (t232 * t196 + t293), (t233 * (-t205 * t236 + t207 * t232) + t237 * (-t205 * t232 - t207 * t236)) * t228, t233 * (-pkin(7) * t137 - t305) + t237 * (-pkin(2) * t141 + pkin(7) * t138 + t294) - pkin(1) * t141 + pkin(6) * (-t233 * t137 + t237 * t138), t233 * (-pkin(7) * t147 - t294) + t237 * (-pkin(2) * t144 + pkin(7) * t148 - t305) - pkin(1) * t144 + pkin(6) * (-t233 * t147 + t237 * t148), t233 * (-pkin(7) * t93 - t62) + t237 * (-pkin(2) * t166 + pkin(7) * t94 + t63) - pkin(1) * t166 + pkin(6) * (-t233 * t93 + t237 * t94), -pkin(7) * t323 + t237 * (pkin(2) * t167 + pkin(7) * t63) + pkin(1) * t167 + pkin(6) * (t237 * t63 - t323), t242, -t380, t372, t332, -t390, t331, t233 * (-t232 * t27 + t236 * t38 - t376) + t237 * (t232 * t38 + t236 * t27 + t374) + t373, t233 * (-pkin(7) * t48 - t232 * t28 + t236 * t43) + t237 * (pkin(7) * t50 + t232 * t43 + t236 * t28 - t388) - t389 + pkin(6) * (-t233 * t48 + t237 * t50), t233 * (-pkin(7) * t40 + t236 * t9) + t237 * (pkin(7) * t42 + t232 * t9) + pkin(6) * (-t233 * t40 + t237 * t42) + (t233 * t326 + t237 * t280 - pkin(1)) * t57, (t233 * (-pkin(8) * t236 + t326) + t237 * (-pkin(8) * t232 + t280) - pkin(1)) * t263 + (pkin(6) + pkin(7)) * (-t233 * t10 + t237 * (t236 * t18 + t232 * t77)), t242, t372, t380, t331, t390, t332, t233 * (t236 * t13 - t232 * t15 - t376) + t237 * (t232 * t13 + t236 * t15 + t374) + t373, t233 * (-pkin(7) * t39 - t232 * t30 + t236 * t5) + t237 * (-pkin(2) * t55 + pkin(7) * t41 + t232 * t5 + t236 * t30) - pkin(1) * t55 + pkin(6) * (-t233 * t39 + t237 * t41), t233 * (-pkin(7) * t44 + t236 * t12 - t232 * t14) + t237 * (pkin(7) * t46 + t232 * t12 + t236 * t14 + t388) + t389 + pkin(6) * (-t233 * t44 + t237 * t46), t233 * (-pkin(7) * t3 - t232 * t1 + t236 * t2) + t237 * (-pkin(2) * t7 + pkin(7) * t4 + t236 * t1 + t232 * t2) - pkin(1) * t7 + pkin(6) * (-t233 * t3 + t237 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t218, t225 - t226, t223, t218, t224, qJDD(2), -t194, -t195, 0, 0, t184, t183, t146, -t184, t251, t286, pkin(2) * t137 - t117, -t295 - t232 * (pkin(7) * t277 - t194 - t325) + (-t232 * t215 + t147) * pkin(2), pkin(2) * t93, pkin(2) * t62, t101, t56, t364, t268, -t87, t253, t284 + t377, pkin(2) * t48 + t285, pkin(2) * t40 + t274, pkin(2) * t10 + t324, t101, t364, -t56, t253, t87, t268, t255 + t377, pkin(2) * t39 + t275, pkin(2) * t44 + t258, pkin(2) * t3 + t250; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t184, t183, t146, -t184, t251, t286, -t117, -t118, 0, 0, t101, t56, t364, t268, -t87, t253, t284, t285, t274, t324, t101, t364, -t56, t253, t87, t268, t255, t275, t258, t250; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t318, t340, t335, -t318, -t108, t162, -t36, -t37, 0, 0, t318, t335, -t340, t162, t108, -t318, t241, t265, t198 + t333, t269; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t337, t335, t341, t31;];
tauJ_reg = t6;
