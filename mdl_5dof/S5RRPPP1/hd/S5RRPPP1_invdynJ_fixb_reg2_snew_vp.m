% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPPP1
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
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPPP1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPP1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPP1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:24:31
% EndTime: 2019-12-31 19:24:46
% DurationCPUTime: 6.14s
% Computational Cost: add. (11917->382), mult. (30797->494), div. (0->0), fcn. (22361->8), ass. (0->253)
t198 = sin(pkin(5));
t200 = cos(pkin(5));
t199 = cos(pkin(8));
t202 = cos(qJ(2));
t273 = qJD(1) * t202;
t179 = -t200 * qJD(2) + t198 * t273;
t178 = t179 ^ 2;
t197 = sin(pkin(8));
t201 = sin(qJ(2));
t272 = qJD(2) * t198;
t277 = t200 * t202;
t161 = t197 * t272 - (-t197 * t277 - t201 * t199) * qJD(1);
t325 = t161 ^ 2;
t333 = -t325 - t178;
t267 = qJD(1) * qJD(2);
t254 = t201 * t267;
t265 = t202 * qJDD(1);
t232 = t254 - t265;
t214 = -t200 * qJDD(2) - t198 * t232;
t335 = t200 * t273 + t272;
t159 = t197 * t201 * qJD(1) - t335 * t199;
t293 = t161 * t159;
t210 = t214 - t293;
t343 = t197 * t210;
t240 = t199 * t333 + t343;
t192 = t201 * qJDD(1);
t253 = t202 * t267;
t182 = t192 + t253;
t266 = t198 * qJDD(2);
t350 = t200 * t232 - t266;
t138 = t199 * t182 - t197 * t350;
t295 = t159 * t179;
t339 = t138 + t295;
t48 = t198 * t240 + t200 * t339;
t413 = pkin(1) * t48;
t412 = pkin(2) * t48;
t209 = t214 + t293;
t281 = t199 * t209;
t326 = t159 ^ 2;
t337 = -t178 - t326;
t238 = -t337 * t197 + t281;
t137 = t197 * t182 + t350 * t199;
t147 = t161 * t179;
t80 = t137 - t147;
t383 = -t198 * t238 + t200 * t80;
t411 = pkin(1) * t383;
t407 = pkin(2) * t383;
t301 = qJ(3) * t198;
t410 = t301 * t48;
t409 = t301 * t383;
t278 = t199 * t200;
t285 = t197 * t200;
t354 = t198 * t339;
t384 = -t210 * t285 - t278 * t333 + t354;
t342 = t199 * t210;
t369 = t197 * t333 - t342;
t391 = t202 * t369;
t408 = pkin(7) * (-t201 * t384 + t391) + t413;
t110 = t326 - t325;
t386 = -t197 * t80 + t199 * t339;
t399 = t202 * (t198 * t110 + t200 * t386) - t201 * (t197 * t339 + t199 * t80);
t406 = pkin(2) * t384;
t363 = t198 * t80;
t385 = -t209 * t278 + t285 * t337 - t363;
t405 = pkin(2) * t385;
t140 = -t325 + t178;
t404 = t201 * (t197 * t140 + t281);
t300 = qJ(3) * t201;
t403 = t300 * t384;
t402 = -t200 * t110 + t198 * t386;
t288 = t197 * t209;
t375 = t199 * t337 + t288;
t390 = t202 * t375;
t398 = pkin(7) * (-t201 * t385 + t390) - t411;
t332 = -t325 - t326;
t334 = -t295 + t138;
t340 = t137 + t147;
t366 = t200 * t332 + (-t197 * t340 - t199 * t334) * t198;
t397 = pkin(2) * t366;
t380 = -t198 * t332 - t334 * t278 - t340 * t285;
t396 = pkin(2) * t380;
t395 = qJ(3) * t369;
t394 = qJ(3) * t375;
t139 = t326 - t178;
t393 = t201 * (t199 * t139 + t343);
t389 = t300 * t380;
t388 = t301 * t366;
t361 = t197 * t334 - t199 * t340;
t382 = pkin(7) * (-t201 * t380 + t202 * t361);
t381 = -pkin(1) * t366 + t382;
t364 = t80 * pkin(3);
t379 = qJ(3) * t361;
t358 = qJ(4) * t339;
t283 = t198 * t199;
t286 = t197 * t198;
t302 = t200 * t334;
t374 = t140 * t283 - t209 * t286 + t302;
t303 = t200 * t340;
t373 = -t139 * t286 + t210 * t283 + t303;
t368 = (-t333 - t178) * pkin(3) - qJ(4) * t210;
t309 = t198 * t340;
t367 = -t393 + t202 * (-t139 * t285 + t210 * t278 - t309);
t308 = t198 * t334;
t365 = -t404 + t202 * (t140 * t278 - t209 * t285 - t308);
t360 = qJ(4) * t332;
t359 = qJ(4) * t337;
t263 = t179 * t285;
t128 = t159 * t263;
t351 = t198 * t214 + t128;
t349 = (-t178 - t332) * pkin(3);
t292 = t161 * t197;
t134 = t179 * t292;
t341 = t201 * (t199 * t138 + t134);
t270 = qJD(4) * t161;
t150 = 0.2e1 * t270;
t269 = qJD(5) * t159;
t338 = -0.2e1 * t269 + t150;
t195 = t201 ^ 2;
t196 = t202 ^ 2;
t336 = -t196 - t195;
t331 = t138 * pkin(4) + t209 * qJ(5);
t299 = qJ(4) * t197;
t149 = -0.2e1 * t270;
t204 = qJD(1) ^ 2;
t323 = sin(qJ(1));
t324 = cos(qJ(1));
t230 = t323 * g(1) - t324 * g(2);
t215 = t204 * pkin(7) + t230;
t318 = t202 * pkin(2);
t256 = -pkin(1) - t318;
t206 = -t182 * t301 + t256 * qJDD(1) + (0.2e1 * t201 * qJD(2) * pkin(2) + (t336 * qJD(1) * t200 - t202 * t272) * qJ(3)) * qJD(1) - t215;
t231 = t324 * g(1) + t323 * g(2);
t297 = qJDD(1) * pkin(7);
t176 = -t204 * pkin(1) - t231 + t297;
t235 = -t198 * t300 - t318;
t177 = t235 * qJD(1);
t317 = t202 * g(3);
t95 = qJDD(2) * pkin(2) - t317 + (-qJD(1) * t177 - t176) * t201 + (qJD(2) * t335 - t182 * t200) * qJ(3);
t57 = t198 * t95 - t200 * t206 - qJDD(3);
t217 = t358 + t57 - t364;
t31 = t149 - t217;
t321 = pkin(3) * t199;
t109 = t159 * pkin(3) - t161 * qJ(4);
t205 = t198 * t206;
t166 = -t201 * g(3) + t202 * t176;
t203 = qJD(2) ^ 2;
t96 = -t203 * pkin(2) + (t200 * t265 + t266) * qJ(3) + t177 * t273 + t166;
t33 = -0.2e1 * qJD(3) * t159 + t197 * t205 + t199 * t96 + t95 * t285;
t233 = -qJ(4) * t214 - 0.2e1 * qJD(4) * t179 - t159 * t109 + t33;
t319 = t178 * pkin(3);
t29 = t233 - t319;
t250 = t197 * t96 - t199 * t205 - t95 * t278;
t223 = pkin(3) * t214 - t178 * qJ(4) + qJDD(4) + t250;
t30 = (0.2e1 * qJD(3) + t109) * t161 + t223;
t9 = t197 * t30 + t199 * t29;
t328 = qJ(3) * t9 + (-t299 - t321) * t31;
t261 = t179 * t278;
t129 = t161 * t261;
t294 = t159 * t199;
t259 = t179 * t294;
t327 = t201 * (-t134 + t259) + t202 * (t129 + t351);
t322 = pkin(3) * t197;
t316 = pkin(3) + qJ(5);
t315 = t197 * t57;
t307 = t199 * t57;
t298 = qJ(4) * t199;
t291 = t161 * t198;
t290 = t161 * t200;
t189 = t202 * t204 * t201;
t276 = t201 * (qJDD(2) + t189);
t275 = t202 * (qJDD(2) - t189);
t271 = qJD(3) * t161;
t268 = qJD(5) * t179;
t264 = t179 * t286;
t262 = t179 * t283;
t260 = t159 * t291;
t258 = t159 * t290;
t126 = t159 * t264;
t127 = t161 * t262;
t257 = -t200 * t214 + t126 + t127;
t165 = t201 * t176 + t317;
t251 = t201 * t165 + t202 * t166;
t249 = -t315 - t395;
t248 = t307 + t394;
t208 = t30 + t331;
t19 = (-pkin(4) * t159 + 0.2e1 * qJD(5)) * t179 + t208;
t132 = t161 * pkin(4) + t179 * qJ(5);
t216 = -t179 * t132 + qJDD(5) + t233;
t213 = -t137 * pkin(4) + t216;
t23 = -qJ(5) * t326 + t213 - t319;
t247 = -t19 * t199 + t197 * t23;
t246 = t197 * t29 - t199 * t30;
t32 = t250 + 0.2e1 * t271;
t245 = t197 * t33 - t199 * t32;
t20 = t197 * t32 + t199 * t33;
t241 = -t140 * t199 + t288;
t239 = t139 * t197 - t342;
t118 = t138 * t286;
t237 = t118 - t127 + t258;
t120 = t137 * t283;
t236 = -t120 - t126 - t258;
t183 = -0.2e1 * t254 + t265;
t211 = pkin(4) * t326 + t161 * t132 + t217;
t207 = t137 * qJ(5) - t211;
t24 = t149 + t207 + 0.2e1 * t269;
t4 = pkin(4) * t23 - t316 * t24;
t7 = t197 * t19 + t199 * t23;
t8 = pkin(4) * t19 - qJ(4) * t24;
t229 = qJ(3) * t7 + t197 * t8 + t199 * t4;
t11 = (-t326 - t332) * qJ(5) + (-t137 - t340) * pkin(4) + t349 + t216;
t17 = -t360 + 0.2e1 * t268 + (t334 - t295) * pkin(4) + t208;
t228 = t11 * t199 + t17 * t197 + t379;
t14 = -t364 + pkin(4) * t337 + (-t137 - t80) * qJ(5) + t211 + t338;
t59 = pkin(4) * t209 - qJ(4) * t80;
t227 = t14 * t199 + t197 * t59 + t394;
t18 = pkin(4) * t333 - t207 + t338 + t358;
t41 = -pkin(4) * t210 + t316 * t339;
t226 = t18 * t197 + t199 * t41 + t395;
t25 = t349 + t233;
t26 = t30 - t360;
t225 = t197 * t26 + t199 * t25 + t379;
t224 = t20 + t379;
t66 = t201 * (t197 * t137 - t259);
t27 = t150 + t217 + t358;
t222 = t197 * t27 + t321 * t339 + t395;
t28 = t31 + t364;
t220 = t199 * t28 + t299 * t80 - t394;
t124 = t137 * t278;
t219 = t66 + t202 * (-t124 - t128 + t260);
t122 = t138 * t285;
t218 = t341 + t202 * (t122 - t129 - t260);
t194 = t196 * t204;
t193 = t195 * t204;
t181 = t192 + 0.2e1 * t253;
t175 = qJDD(1) * pkin(1) + t215;
t77 = qJ(4) * t340;
t58 = -pkin(3) * t334 - t77;
t51 = t200 * t240 - t354;
t47 = t200 * t238 + t363;
t40 = -t316 * t334 - t77;
t22 = pkin(3) * t209 + t30 - t359;
t21 = t233 + t368;
t16 = t198 * t57 + t200 * t245;
t15 = t198 * t245 - t200 * t57;
t13 = (-t333 - t326) * qJ(5) + t213 + t368;
t12 = pkin(4) * t295 - t161 * t109 - t209 * t316 - t223 - 0.2e1 * t268 - 0.2e1 * t271 - t331 + t359;
t10 = -pkin(3) * t30 + qJ(4) * t29;
t6 = -t198 * t31 + t200 * t246;
t5 = t198 * t246 + t200 * t31;
t3 = -t198 * t24 + t200 * t247;
t2 = t198 * t247 + t200 * t24;
t1 = qJ(4) * t23 - t316 * t19;
t34 = [0, 0, 0, 0, 0, qJDD(1), t230, t231, 0, 0, (t182 + t253) * t201, t202 * t181 + t201 * t183, t276 + t202 * (-t193 + t203), t183 * t202, t201 * (t194 - t203) + t275, 0, t202 * t175 + pkin(1) * t183 + pkin(7) * (t202 * (-t194 - t203) - t276), -t201 * t175 - pkin(1) * t181 + pkin(7) * (-t275 - t201 * (-t193 - t203)), pkin(1) * (t193 + t194) - t336 * t297 + t251, pkin(1) * t175 + pkin(7) * t251, t218, t399, t365, t66 + t202 * (-t124 + (-t263 + t291) * t159), t393 + t202 * (t200 * t239 + t309), t327, t201 * (-t315 + (-t198 * t383 - t200 * t385) * qJ(3)) + t202 * (t198 * t32 + t200 * t248 - t407) + t398, t201 * (-t307 + (-t198 * t48 - t200 * t51) * qJ(3)) + t202 * (t198 * t33 + t200 * t249 - t412) - t413 + pkin(7) * (-t201 * t51 - t391), -t201 * t245 + t382 + (-pkin(1) + t235) * t366 + (t202 * t224 - t389) * t200, pkin(7) * (-t201 * t16 + t202 * t20) + t256 * t15 + (t201 * (-t15 * t198 - t16 * t200) + t20 * t277) * qJ(3), t202 * t351 + (t201 * (-t292 + t294) + t161 * t199 * t277) * t179, t404 + t202 * (t200 * t241 + t308), t367, t341 + t202 * (t122 + (-t159 * t198 - t261) * t161), t399, t219, t201 * (-t197 * t25 + t199 * t26 - t388) + t202 * (-t198 * t58 - t397) + (t202 * t225 - t389) * t200 + t381, t201 * (-t197 * t28 + t298 * t80 + t409) + t202 * (-t198 * t22 + t407) + t411 + pkin(7) * (-t201 * t47 - t390) + (t202 * t220 - t47 * t300) * t200, t201 * (t199 * t27 - t322 * t339 + t410) + t202 * (-t198 * t21 + t412) + (t202 * t222 - t403) * t200 + t408, t201 * (-t5 * t301 + (-t298 + t322) * t31) + t202 * (-pkin(2) * t5 - t198 * t10) - pkin(1) * t5 + pkin(7) * (-t201 * t6 + t202 * t9) + (t202 * t328 - t6 * t300) * t200, t327, t367, t365, t219, -t399, t218, t201 * (-t197 * t11 + t199 * t17 - t388) + t202 * (-t198 * t40 - t397) + (t202 * t228 - t389) * t200 + t381, t201 * (t199 * t18 - t197 * t41 + t410) + t202 * (-t198 * t13 + t412) + (t202 * t226 - t403) * t200 + t408, t201 * (-t197 * t14 + t199 * t59 - t409) + t202 * (-t198 * t12 - t407) + (t202 * t227 - t300 * t385) * t200 + t398, t201 * (-t197 * t4 + t199 * t8 - t2 * t301) + t202 * (-pkin(2) * t2 - t198 * t1) - pkin(1) * t2 + pkin(7) * (-t201 * t3 + t202 * t7) + (t202 * t229 - t3 * t300) * t200; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t189, t193 - t194, t192, t189, t265, qJDD(2), -t165, -t166, 0, 0, t237, t402, t374, -t120 + (-t264 - t290) * t159, t198 * t239 - t303, t257, t198 * t248 - t200 * t32 + t405, pkin(2) * t51 + t198 * t249 - t200 * t33, t198 * t224 + t396, pkin(2) * t16 + t20 * t301, t257, t198 * t241 - t302, t373, t118 + (t159 * t200 - t262) * t161, t402, t236, t198 * t225 + t200 * t58 + t396, pkin(2) * t47 + t198 * t220 + t200 * t22, t198 * t222 + t200 * t21 + t406, pkin(2) * t6 + t200 * t10 + t198 * t328, t257, t373, t374, t236, -t402, t237, t198 * t228 + t200 * t40 + t396, t200 * t13 + t198 * t226 + t406, t200 * t12 + t198 * t227 + t405, pkin(2) * t3 + t200 * t1 + t198 * t229; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t339, t332, -t57, 0, 0, 0, 0, 0, 0, t332, -t80, -t339, t31, 0, 0, 0, 0, 0, 0, t332, -t339, t80, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t334, -t209, t333, t30, 0, 0, 0, 0, 0, 0, t334, t333, t209, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t340, -t210, t337, t23;];
tauJ_reg = t34;
