% Calculate minimal parameter regressor of coriolis matrix for
% S5RRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x26]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRRRP1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:46:03
% EndTime: 2019-12-05 18:46:12
% DurationCPUTime: 3.96s
% Computational Cost: add. (5733->262), mult. (10988->352), div. (0->0), fcn. (12165->6), ass. (0->217)
t202 = sin(qJ(3));
t203 = sin(qJ(2));
t205 = cos(qJ(2));
t337 = cos(qJ(3));
t179 = t202 * t203 - t205 * t337;
t348 = pkin(6) + pkin(7);
t187 = t348 * t205;
t183 = t337 * t187;
t186 = t348 * t203;
t301 = t202 * t186;
t353 = -t183 + t301;
t111 = pkin(8) * t179 + t353;
t204 = cos(qJ(4));
t372 = t204 * t111;
t100 = t372 / 0.2e1;
t201 = sin(qJ(4));
t116 = t186 * t337 + t187 * t202;
t181 = -t202 * t205 - t203 * t337;
t355 = pkin(8) * t181 - t116;
t361 = t201 * t355;
t390 = t100 - t361 / 0.2e1;
t396 = 0.2e1 * t390;
t238 = -t179 * t204 + t181 * t201;
t378 = t372 - t361;
t387 = -qJ(5) * t238 + t378;
t347 = t387 * pkin(4);
t268 = t347 / 0.2e1;
t150 = t201 * pkin(3);
t360 = t204 * t355;
t373 = t201 * t111;
t380 = -t360 - t373;
t87 = t179 * t201 + t181 * t204;
t389 = qJ(5) * t87 - t380;
t395 = t389 * t150;
t247 = -t372 / 0.2e1;
t388 = t247 + t100;
t394 = qJD(1) * t388;
t393 = qJD(3) * t388;
t392 = qJD(3) + qJD(4);
t367 = -t360 / 0.2e1;
t375 = t367 - t373 / 0.2e1;
t391 = 0.2e1 * t375;
t386 = qJD(2) * t391;
t385 = qJD(2) * t396;
t269 = qJD(2) + qJD(3);
t248 = t360 / 0.2e1;
t22 = t367 + t248;
t384 = t22 * qJD(1);
t383 = t22 * qJD(4);
t267 = t337 * pkin(2);
t235 = t267 + pkin(3);
t216 = t201 * t235;
t300 = t202 * t204;
t166 = pkin(2) * t300 + t216;
t343 = t166 / 0.2e1;
t381 = t389 * t343;
t134 = t238 ^ 2;
t351 = t87 ^ 2;
t369 = t134 - t351;
t376 = t369 * qJD(1);
t274 = t87 * qJD(4);
t53 = t269 * t87 + t274;
t368 = t87 * pkin(4);
t335 = pkin(4) * t238;
t286 = qJD(1) * t87;
t359 = t238 * t286;
t275 = qJD(4) * t238;
t54 = t238 * t269 + t275;
t188 = t204 * t235;
t302 = t201 * t202;
t165 = pkin(2) * t302 - t188;
t158 = pkin(4) - t165;
t354 = t158 + t165;
t352 = t269 * t116;
t350 = -t389 / 0.2e1;
t344 = t158 / 0.2e1;
t258 = t337 * t201;
t172 = (t258 + t300) * pkin(2);
t342 = -t172 / 0.2e1;
t257 = t337 * t204;
t173 = (t257 - t302) * pkin(2);
t341 = t173 / 0.2e1;
t234 = -t183 / 0.2e1;
t329 = t204 * pkin(3);
t195 = pkin(4) + t329;
t340 = -t195 / 0.2e1;
t339 = t195 / 0.2e1;
t338 = t204 / 0.2e1;
t336 = pkin(2) * t202;
t332 = t172 * pkin(4);
t330 = t181 * pkin(3);
t200 = t203 * pkin(2);
t324 = pkin(3) * qJD(3);
t323 = pkin(3) * qJD(4);
t233 = -t330 - t368;
t196 = -pkin(2) * t205 - pkin(1);
t154 = pkin(3) * t179 + t196;
t75 = t154 - t335;
t10 = t75 * (t200 + t233);
t318 = t10 * qJD(1);
t11 = t233 * t75;
t317 = t11 * qJD(1);
t12 = t368 * t75;
t316 = t12 * qJD(1);
t13 = -t238 * t387 + t389 * t87;
t315 = t13 * qJD(1);
t212 = -t238 * t341 - t342 * t87;
t261 = -t150 / 0.2e1;
t14 = (t340 + t344) * t238 - (t261 + t343) * t87 + t212;
t313 = t14 * qJD(1);
t312 = t158 * t238;
t311 = t166 * t87;
t246 = t165 / 0.2e1 + t344;
t231 = t246 * t238;
t262 = t335 / 0.2e1;
t17 = t262 - t231;
t310 = t17 * qJD(1);
t309 = t195 * t238;
t213 = -t238 * t343 - t344 * t87;
t287 = -t368 / 0.2e1 - t330 / 0.2e1;
t242 = t200 / 0.2e1 + t287;
t30 = t213 + t242;
t297 = t30 * qJD(1);
t249 = t201 * t238 / 0.2e1;
t250 = t87 * t340;
t45 = -t250 + t368 / 0.2e1 + (t249 + t181 / 0.2e1) * pkin(3);
t295 = t45 * qJD(1);
t260 = -t329 / 0.2e1;
t215 = (t260 + t339) * t238;
t63 = t262 - t215;
t294 = t63 * qJD(1);
t155 = t200 - t330;
t76 = t154 * t87;
t66 = -t155 * t238 - t76;
t293 = t66 * qJD(1);
t77 = t154 * t238;
t67 = -t155 * t87 + t77;
t292 = t67 * qJD(1);
t68 = -t238 * t330 + t76;
t291 = t68 * qJD(1);
t69 = -t330 * t87 - t77;
t290 = t69 * qJD(1);
t72 = t134 + t351;
t288 = t72 * qJD(1);
t285 = qJD(1) * t154;
t284 = qJD(1) * t196;
t283 = qJD(1) * t205;
t282 = qJD(3) * t196;
t110 = t179 ^ 2 - t181 ^ 2;
t281 = t110 * qJD(1);
t121 = t179 * t200 - t181 * t196;
t278 = t121 * qJD(1);
t122 = -t179 * t196 - t181 * t200;
t277 = t122 * qJD(1);
t147 = t234 + t183 / 0.2e1;
t273 = t147 * qJD(1);
t157 = t166 * qJD(4);
t191 = -t203 ^ 2 + t205 ^ 2;
t272 = t191 * qJD(1);
t271 = t203 * qJD(2);
t270 = t205 * qJD(2);
t266 = pkin(1) * t203 * qJD(1);
t265 = pkin(1) * t283;
t264 = pkin(4) * t286;
t263 = t201 * t323;
t259 = t389 / 0.2e1 + t350;
t255 = t238 * t285;
t254 = t87 * t285;
t253 = t179 * t284;
t252 = t181 * t284;
t251 = t203 * t283;
t245 = t337 * qJD(2);
t244 = t337 * qJD(3);
t243 = pkin(3) * t392;
t142 = t269 * t181;
t208 = t342 * t389 + t381 + (-t341 + t344) * t387;
t210 = t387 * t339 + t395 / 0.2e1;
t1 = -t208 + t210;
t74 = -t158 * t172 + t166 * t173;
t227 = -qJD(1) * t1 + qJD(2) * t74;
t3 = -t166 * t259 - t246 * t387 + t268;
t71 = t354 * t166;
t226 = qJD(1) * t3 + qJD(2) * t71;
t222 = -t150 * qJD(2) - t394;
t151 = t188 / 0.2e1 + (-t267 / 0.2e1 + pkin(3) / 0.2e1) * t204;
t221 = -t151 * qJD(2) + t384;
t220 = -qJD(2) * t165 - t384;
t21 = t247 + t361 / 0.2e1 + t390;
t219 = qJD(1) * t21 + qJD(2) * t166;
t218 = qJD(2) * t172 - t394;
t26 = t248 + t373 / 0.2e1 + t375;
t217 = -qJD(1) * t26 + qJD(2) * t173;
t159 = (-t195 + t329) * t150;
t207 = (t201 * t259 - t338 * t387) * pkin(3) - t387 * t340;
t6 = -t347 / 0.2e1 + t207;
t206 = (t166 * t338 - t201 * t246) * pkin(3) + t166 * t340;
t65 = t332 / 0.2e1 + t206;
t209 = -qJD(1) * t6 - qJD(2) * t65 - qJD(3) * t159;
t168 = t173 * qJD(3);
t167 = t172 * qJD(3);
t156 = t165 * qJD(4);
t145 = t181 * t179 * qJD(1);
t144 = t260 - t188 / 0.2e1 + (t302 - t257 / 0.2e1) * pkin(2);
t143 = t261 - t216 / 0.2e1 + (-t300 - t258 / 0.2e1) * pkin(2);
t141 = t269 * t179;
t130 = -t335 / 0.2e1;
t117 = 0.2e1 * t234 + t301;
t64 = -t332 / 0.2e1 + t206;
t62 = t130 - t215;
t46 = pkin(3) * t249 - t250 + t287;
t31 = -t213 + t242;
t16 = t130 - t231;
t15 = t311 / 0.2e1 - t312 / 0.2e1 - t87 * t261 - t309 / 0.2e1 - t212;
t5 = t268 + t207;
t4 = t166 * t350 + t381 + t268 + t354 * t387 / 0.2e1;
t2 = t208 + t210;
t7 = [0, 0, 0, t203 * t270, t191 * qJD(2), 0, 0, 0, -pkin(1) * t271, -pkin(1) * t270, t179 * t142, t269 * t110, 0, 0, 0, qJD(2) * t121 - t181 * t282, qJD(2) * t122 - t179 * t282, -t54 * t87, (qJD(4) + t269) * t369, 0, 0, 0, qJD(2) * t66 - qJD(3) * t68 - t154 * t274, qJD(2) * t67 - qJD(3) * t69 + t154 * t275, t72 * qJD(5), qJD(2) * t10 + qJD(3) * t11 - qJD(4) * t12 + qJD(5) * t13; 0, 0, 0, t251, t272, t270, -t271, 0, -pkin(6) * t270 - t266, pkin(6) * t271 - t265, t145, t281, -t141, t142, 0, qJD(2) * t353 + qJD(3) * t117 + t278, t277 + t352, -t359, t376, t54, t53, 0, qJD(2) * t378 + t392 * t396 + t293, qJD(2) * t380 + t391 * t392 + t292, (t311 - t312) * qJD(2) + t15 * qJD(3) + t16 * qJD(4), t318 + (t158 * t387 + t166 * t389) * qJD(2) + t2 * qJD(3) + t4 * qJD(4) + t31 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t145, t281, -t141, t142, 0, qJD(2) * t117 + qJD(3) * t353 - t252, -t253 + t352, -t359, t376, t54, t53, 0, qJD(3) * t378 + qJD(4) * t396 - t291 + t385, qJD(3) * t380 + qJD(4) * t391 - t290 + t386, t15 * qJD(2) + (t150 * t87 - t309) * qJD(3) + t62 * qJD(4), t317 + t2 * qJD(2) + (t195 * t387 + t395) * qJD(3) + t5 * qJD(4) + t46 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t359, t376, t54, t53, 0, qJD(3) * t396 + qJD(4) * t378 - t254 + t385, qJD(3) * t391 + qJD(4) * t380 + t255 + t386, -pkin(4) * t275 + qJD(2) * t16 + qJD(3) * t62, qJD(2) * t4 + qJD(3) * t5 + qJD(4) * t347 - t316; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t288, qJD(2) * t31 + qJD(3) * t46 + t315; 0, 0, 0, -t251, -t272, 0, 0, 0, t266, t265, -t145, -t281, 0, 0, 0, qJD(3) * t147 - t278, -t277, t359, -t376, 0, 0, 0, -qJD(4) * t21 - t293 + t393, qJD(3) * t26 - t292 + t383, -qJD(3) * t14 + qJD(4) * t17, -qJD(3) * t1 - qJD(4) * t3 - qJD(5) * t30 - t318; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t336, -pkin(2) * t244, 0, 0, 0, 0, 0, -t167 - t157, -t168 + t156, 0, qJD(3) * t74 - qJD(4) * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t269 * t336 + t273, (-t245 - t244) * pkin(2), 0, 0, 0, 0, 0, qJD(4) * t143 - t167 - t218, qJD(4) * t144 - t168 - t217, -t313, (t150 * t173 - t172 * t195) * qJD(3) + t64 * qJD(4) + t227; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t143 - t157 - t219, qJD(3) * t144 + t156 - t220, t310, -pkin(4) * t157 + qJD(3) * t64 - t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t297; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t145, -t281, 0, 0, 0, -qJD(2) * t147 + t252, t253, t359, -t376, 0, 0, 0, t291 + (-qJD(2) - qJD(4)) * t388, -qJD(2) * t26 + t290 + t383, qJD(2) * t14 + qJD(4) * t63, qJD(2) * t1 + qJD(4) * t6 + qJD(5) * t45 - t317; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t336 - t273, pkin(2) * t245, 0, 0, 0, 0, 0, -qJD(4) * t150 + t218, -qJD(4) * t151 + t217, t313, qJD(4) * t65 - t227; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t263, -t204 * t323, 0, t159 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t201 * t243 + t222, -t204 * t243 + t221, t294, -pkin(4) * t263 - t209; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t295; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t359, -t376, 0, 0, 0, qJD(2) * t21 + t254 + t393, -t22 * t269 - t255, -qJD(2) * t17 - qJD(3) * t63, qJD(2) * t3 - qJD(3) * t6 + qJD(5) * t368 + t316; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t150 + t219, qJD(3) * t151 + t220, -t310, -qJD(3) * t65 + t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t201 * t324 - t222, t204 * t324 - t221, -t294, t209; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t264; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t288, -pkin(4) * t274 + qJD(2) * t30 - qJD(3) * t45 - t315; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t297; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t295; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t264; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t7;