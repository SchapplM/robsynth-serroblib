% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x27]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRRR13_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR13_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:15:30
% EndTime: 2019-12-31 19:15:39
% DurationCPUTime: 3.86s
% Computational Cost: add. (2514->342), mult. (5727->485), div. (0->0), fcn. (5607->6), ass. (0->291)
t306 = qJD(4) + qJD(5);
t232 = sin(qJ(3));
t227 = t232 ^ 2;
t401 = t227 / 0.2e1;
t231 = sin(qJ(4));
t387 = pkin(7) + pkin(8);
t185 = t387 * t231;
t234 = cos(qJ(4));
t187 = t387 * t234;
t230 = sin(qJ(5));
t233 = cos(qJ(5));
t261 = -t233 * t185 - t230 * t187;
t363 = t261 * t232;
t400 = t363 / 0.2e1;
t126 = -t230 * t185 + t233 * t187;
t399 = t306 * t126;
t398 = t306 * t261;
t397 = 0.2e1 * t231;
t341 = t233 * t231;
t349 = t230 * t234;
t171 = t341 + t349;
t235 = cos(qJ(3));
t336 = t235 * t171;
t396 = t306 * t336;
t221 = t233 * t234;
t350 = t230 * t231;
t167 = -t221 + t350;
t236 = -pkin(1) - pkin(6);
t379 = pkin(4) * t231;
t278 = -t236 + t379;
t265 = t278 * t235;
t258 = -t265 / 0.2e1;
t378 = pkin(4) * t234;
t223 = -pkin(3) - t378;
t355 = t223 * t336;
t395 = t167 * t258 - t355 / 0.2e1;
t257 = t265 / 0.2e1;
t339 = t234 * t235;
t189 = t233 * t339;
t347 = t231 * t235;
t143 = t230 * t347 - t189;
t356 = t223 * t143;
t394 = t171 * t257 - t356 / 0.2e1;
t376 = t235 * pkin(3);
t377 = t232 * pkin(7);
t186 = t376 + t377;
t173 = t231 * t186;
t335 = t235 * t236;
t191 = t234 * t335;
t348 = t231 * t232;
t118 = pkin(8) * t348 + t173 + t191;
t342 = t233 * t118;
t281 = -t342 / 0.2e1;
t174 = t234 * t186;
t346 = t231 * t236;
t279 = pkin(4) - t346;
t102 = t234 * t232 * pkin(8) + t279 * t235 + t174;
t354 = t230 * t102;
t251 = -t354 / 0.2e1 + t281;
t239 = t355 / 0.2e1 + t251;
t10 = t167 * t257 + t239 + t400;
t319 = qJD(3) * t223;
t282 = t350 / 0.2e1;
t383 = -t189 / 0.2e1;
t83 = t383 + (t282 - t167 / 0.2e1) * t235;
t327 = t83 * qJD(2);
t393 = qJD(1) * t10 + t167 * t319 + t327;
t249 = t349 / 0.2e1 + t341 / 0.2e1;
t385 = t171 / 0.2e1;
t85 = (t385 - t249) * t235;
t326 = t85 * qJD(2);
t351 = t230 * t118;
t283 = -t351 / 0.2e1;
t345 = t233 * t102;
t250 = t283 + t345 / 0.2e1;
t240 = t356 / 0.2e1 + t250;
t362 = t126 * t232;
t284 = t362 / 0.2e1;
t9 = t171 * t258 + t240 + t284;
t392 = qJD(1) * t9 - t171 * t319 + t326;
t303 = t379 / 0.2e1;
t237 = -t303 * t336 + t284;
t302 = -t378 / 0.2e1;
t4 = (t233 * pkin(4) / 0.2e1 + t167 * t302 - t278 * t171 / 0.2e1) * t235 + t237 + t240;
t81 = t167 * t379 + t171 * t223;
t391 = t4 * qJD(1) - t81 * qJD(3) + t326;
t238 = t143 * t303 + t400;
t386 = t167 / 0.2e1;
t3 = (-t230 * pkin(4) / 0.2e1 + t171 * t302 + t278 * t386) * t235 + t238 + t239;
t82 = -t167 * t223 + t171 * t379;
t390 = t3 * qJD(1) - t82 * qJD(3) + t327;
t226 = t231 ^ 2;
t228 = t234 ^ 2;
t201 = t228 - t226;
t276 = t339 * t397;
t241 = qJD(1) * t276 - qJD(3) * t201;
t269 = pkin(3) * t232 - pkin(7) * t235;
t176 = qJ(2) + t269;
t164 = t234 * t176;
t305 = pkin(8) * t339;
t99 = t279 * t232 + t164 - t305;
t80 = t233 * t99;
t389 = -t80 / 0.2e1;
t388 = -t99 / 0.2e1;
t188 = t230 * t348;
t384 = t188 / 0.2e1;
t381 = -t232 / 0.2e1;
t380 = t235 / 0.2e1;
t33 = t171 * t143 + t167 * t336;
t375 = t306 * t33;
t45 = t143 * t386 - t336 * t385;
t374 = t306 * t45;
t373 = pkin(4) * qJD(4);
t372 = pkin(4) * qJD(5);
t144 = t171 * t232;
t166 = t278 * t232;
t338 = t234 * t236;
t299 = t232 * t338;
t131 = t231 * t176 + t299;
t115 = -pkin(8) * t347 + t131;
t352 = t230 * t115;
t37 = -t80 + t352;
t1 = (t345 - t351) * t232 - t166 * t336 + (-t278 * t144 - t37) * t235;
t371 = t1 * qJD(1);
t147 = t232 * t221 - t188;
t343 = t233 * t115;
t38 = t230 * t99 + t343;
t2 = (t342 + t354) * t232 + t38 * t235 - t166 * t143 + t147 * t265;
t370 = t2 * qJD(1);
t130 = t232 * t346 - t164;
t114 = -t130 - t305;
t353 = t230 * t114;
t41 = -t343 - t353;
t21 = t41 * t232 + (-t278 * t143 + t336 * t378) * t235;
t369 = qJD(1) * t21;
t108 = t336 * t265;
t337 = t235 * t143;
t344 = t233 * t114;
t42 = t344 - t352;
t22 = t232 * t42 + t337 * t378 + t108;
t368 = qJD(1) * t22;
t23 = -t232 * t37 + t108;
t367 = qJD(1) * t23;
t24 = -t143 * t265 - t38 * t232;
t366 = qJD(1) * t24;
t229 = t235 ^ 2;
t91 = -t130 * t232 - t229 * t346;
t365 = qJD(1) * t91;
t92 = -t131 * t232 - t229 * t338;
t364 = qJD(1) * t92;
t301 = pkin(4) * t381;
t272 = t301 + t114 / 0.2e1;
t13 = (t388 + t272) * t230;
t361 = t13 * qJD(1);
t360 = t144 * t232;
t359 = t336 * t235;
t358 = t147 * t232;
t15 = t272 * t233 + t389;
t357 = t15 * qJD(1);
t222 = t229 * t234;
t340 = t234 * t227;
t36 = -t143 * t144 + t147 * t336;
t334 = t36 * qJD(1);
t298 = t231 * t335;
t46 = -t130 * t235 + (t174 + t298) * t232;
t333 = t46 * qJD(1);
t47 = t131 * t235 + (-t191 + t173) * t232;
t332 = t47 * qJD(1);
t252 = -t358 / 0.2e1 + t337 / 0.2e1;
t266 = t221 / 0.2e1 - t350 / 0.2e1;
t54 = -t252 + t266;
t331 = t54 * qJD(1);
t253 = t360 / 0.2e1 + t359 / 0.2e1;
t55 = -t249 - t253;
t330 = t55 * qJD(1);
t63 = -t359 + t360;
t329 = t63 * qJD(1);
t64 = -t337 - t358;
t328 = t64 * qJD(1);
t84 = (t385 + t249) * t232;
t71 = t84 * qJD(1);
t280 = -t221 / 0.2e1;
t86 = t384 + (t280 + t386) * t232;
t73 = t86 * qJD(1);
t200 = t227 - t229;
t325 = qJD(1) * qJ(2);
t324 = qJD(1) * t143;
t170 = t200 * t231;
t323 = qJD(1) * t170;
t172 = -t222 + t340;
t322 = qJD(1) * t172;
t321 = qJD(2) * t232;
t320 = qJD(3) * t171;
t318 = qJD(3) * t234;
t317 = qJD(4) * t231;
t316 = qJD(4) * t234;
t315 = qJD(5) * t223;
t304 = 0.1e1 / 0.2e1 + t401;
t148 = (-t229 / 0.2e1 - t304) * t231;
t314 = t148 * qJD(1);
t149 = t222 / 0.2e1 + t304 * t234;
t313 = t149 * qJD(1);
t312 = t200 * qJD(1);
t311 = t232 * qJD(1);
t310 = t232 * qJD(3);
t309 = t235 * qJD(1);
t308 = t235 * qJD(3);
t307 = t235 * qJD(4);
t300 = pkin(4) * t380;
t297 = qJ(2) * t311;
t296 = qJ(2) * t309;
t295 = t231 * t318;
t294 = t231 * t308;
t293 = t234 * t308;
t292 = t232 * t317;
t291 = t231 * t307;
t290 = t232 * t316;
t289 = t234 * t307;
t288 = t167 * t311;
t287 = t171 * t311;
t286 = t231 * t316;
t215 = t232 * t308;
t285 = t232 * t309;
t277 = pkin(4) * t306;
t275 = t306 * t171;
t274 = t234 * t300;
t273 = -qJD(4) - t311;
t271 = t311 + qJD(4) / 0.2e1;
t270 = t300 + t102 / 0.2e1;
t268 = qJD(3) * t276;
t262 = qJD(5) - t273;
t58 = -t143 ^ 2 + t336 ^ 2;
t7 = qJD(1) * t58 + qJD(3) * t33;
t65 = t167 ^ 2 - t171 ^ 2;
t17 = qJD(1) * t33 + qJD(3) * t65;
t260 = t273 * t235;
t259 = t377 / 0.2e1 + t376 / 0.2e1;
t244 = t259 * t231;
t116 = t173 / 0.2e1 + t244;
t256 = pkin(3) * t318 - qJD(1) * t116;
t245 = t259 * t234;
t117 = -t174 / 0.2e1 - t245;
t255 = pkin(3) * qJD(3) * t231 - qJD(1) * t117;
t26 = -qJD(3) * t45 - t324 * t336;
t35 = -qJD(1) * t45 + t167 * t320;
t247 = t234 * t260;
t157 = (t226 / 0.2e1 - t228 / 0.2e1) * t235;
t246 = -qJD(1) * t157 + t295;
t243 = qJD(1) * t231 * t222 + qJD(3) * t157;
t169 = t201 * t229;
t242 = qJD(1) * t169 + t268;
t220 = -t309 / 0.2e1;
t219 = t309 / 0.2e1;
t218 = t308 / 0.2e1;
t214 = t234 * t311;
t213 = t231 * t310;
t212 = t231 * t311;
t165 = t271 * t235;
t156 = t157 * qJD(4);
t151 = -t340 / 0.2e1 - t222 / 0.2e1 + t234 / 0.2e1;
t150 = (-0.1e1 / 0.2e1 + t401 + t229 / 0.2e1) * t231;
t142 = (qJD(5) / 0.2e1 + t271) * t235;
t90 = -t144 / 0.2e1 + t249 * t232;
t89 = t167 * t381 + t232 * t280 + t384;
t88 = -t336 / 0.2e1 - t249 * t235;
t87 = t167 * t380 + t235 * t282 + t383;
t77 = -t298 + t174 / 0.2e1 - t245;
t76 = -t191 - t173 / 0.2e1 + t244;
t57 = t252 + t266;
t56 = -t249 + t253;
t53 = t57 * qJD(2);
t52 = t56 * qJD(2);
t51 = t55 * qJD(2);
t50 = t54 * qJD(2);
t49 = qJD(3) * t86 + t311 * t336;
t48 = qJD(3) * t84 - t143 * t311;
t40 = -t275 - t71;
t39 = -t306 * t167 - t73;
t30 = qJD(3) * t83 + t330;
t29 = qJD(3) * t85 + t331;
t28 = t89 * qJD(3) - t262 * t336;
t27 = t90 * qJD(3) + t143 * t262;
t20 = t88 * qJD(3) - t306 * t147 - t331;
t19 = t87 * qJD(3) + t306 * t144 - t330;
t16 = t233 * t301 + t352 + t389 - t344 / 0.2e1;
t14 = -t343 - t353 / 0.2e1 + (t301 + t388) * t230;
t12 = -t362 / 0.2e1 + t250 + t394;
t11 = -t363 / 0.2e1 + t251 + t395;
t6 = t171 * t274 - t270 * t230 - t238 + t281 + t395;
t5 = t167 * t274 + t270 * t233 - t237 + t283 + t394;
t8 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), -t215, t200 * qJD(3), 0, 0, 0, qJ(2) * t308 + t321, -qJ(2) * t310 + qJD(2) * t235, -t215 * t228 - t229 * t286, -qJD(4) * t169 + t232 * t268, -qJD(3) * t172 - t232 * t291, qJD(3) * t170 - t232 * t289, t215, qJD(3) * t46 + qJD(4) * t92 + t234 * t321, -qJD(3) * t47 - qJD(4) * t91 - t231 * t321, (qJD(3) * t147 + t396) * t143, t36 * qJD(3) + t306 * t58, t64 * qJD(3) - t232 * t396, t143 * t232 * t306 + t63 * qJD(3), t215, qJD(3) * t1 + qJD(4) * t21 + qJD(5) * t24 - t167 * t321, -qJD(3) * t2 - qJD(4) * t22 - qJD(5) * t23 - t171 * t321; 0, 0, 0, 0, qJD(1), t325, 0, 0, 0, 0, 0, t311, t309, 0, 0, 0, 0, 0, qJD(4) * t151 + t214, qJD(4) * t150 - t212, 0, 0, 0, 0, 0, t306 * t57 - t288, t306 * t56 - t287; 0, 0, 0, 0, 0, 0, -t285, t312, -t310, -t308, 0, -t236 * t310 + t296, -t236 * t308 - t297, -t156 + (-t228 * t309 - t295) * t232, t241 * t232 - 0.2e1 * t235 * t286, t294 - t322, t293 + t323, t165, t333 + (t231 * t269 - t299) * qJD(3) + t77 * qJD(4), -t332 + (-pkin(7) * t339 + (pkin(3) * t234 + t346) * t232) * qJD(3) + t76 * qJD(4), (-t320 + t324) * t147 + t374, t334 + (t144 * t171 + t147 * t167) * qJD(3) + t375, t171 * t308 + t306 * t89 + t328, -t167 * t308 + t306 * t90 + t329, t142, t371 + (-t144 * t223 - t166 * t167 + t235 * t261) * qJD(3) + t5 * qJD(4) + t12 * qJD(5), -t370 + (-t126 * t235 - t147 * t223 - t166 * t171) * qJD(3) + t6 * qJD(4) + t11 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t243, -t242, t231 * t260, t247, t218, qJD(2) * t151 + qJD(3) * t77 - qJD(4) * t131 + t364, qJD(2) * t150 + qJD(3) * t76 + qJD(4) * t130 - t365, -t26, t7, t28, t27, t218, qJD(3) * t5 + qJD(4) * t41 + qJD(5) * t14 + t369 + t53, qJD(3) * t6 - qJD(4) * t42 + qJD(5) * t16 - t368 + t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t7, t28, t27, t218, qJD(3) * t12 + qJD(4) * t14 - qJD(5) * t38 + t366 + t53, qJD(3) * t11 + qJD(4) * t16 + qJD(5) * t37 - t367 + t52; 0, 0, 0, 0, -qJD(1), -t325, 0, 0, 0, 0, 0, -t311, -t309, 0, 0, 0, 0, 0, -qJD(4) * t149 - t214, -qJD(4) * t148 + t212, 0, 0, 0, 0, 0, -t306 * t54 + t288, -t306 * t55 + t287; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t310, -t308, 0, 0, 0, 0, 0, -t234 * t310 - t291, t213 - t289, 0, 0, 0, 0, 0, t167 * t310 + t306 * t88, t171 * t310 + t306 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t290 - t294 - t313, t292 - t293 - t314, 0, 0, 0, 0, 0, t20, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t19; 0, 0, 0, 0, 0, 0, t285, -t312, 0, 0, 0, -t296, t297, t228 * t285 - t156, t247 * t397, t290 + t322, -t292 - t323, -t165, qJD(4) * t117 - t333, qJD(4) * t116 + t332, -t147 * t324 + t374, -t334 + t375, -t306 * t86 - t328, -t306 * t84 - t329, -t142, -qJD(4) * t4 - qJD(5) * t9 - t371, -qJD(4) * t3 - qJD(5) * t10 + t370; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t306 * t85, -t306 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t286, t201 * qJD(4), 0, 0, 0, -pkin(3) * t317, -pkin(3) * t316, -t167 * t275, t306 * t65, 0, 0, 0, qJD(4) * t81 + t171 * t315, qJD(4) * t82 - t167 * t315; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t246, -t241, t214 + t316, -t212 - t317, t220, -pkin(7) * t316 - t255, pkin(7) * t317 - t256, -t35, t17, t39, t40, t220, -t391 - t399, -t390 - t398; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, t17, t39, t40, t220, -t392 - t399, -t393 - t398; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t243, t242, (t231 * t309 - t318) * t232, t234 * t285 + t213, t218, qJD(2) * t149 - qJD(3) * t117 - t364, qJD(2) * t148 - qJD(3) * t116 + t365, t26, -t7, t49, t48, t218, qJD(3) * t4 + qJD(5) * t13 - t369 + t50, qJD(3) * t3 + qJD(5) * t15 + t368 + t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t313, t314, 0, 0, 0, 0, 0, t29, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t246, t241, -t214, t212, t219, t255, t256, t35, -t17, t73, t71, t219, t391, t390; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t230 * t372, -t233 * t372; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t230 * t277 + t361, -t233 * t277 + t357; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t7, t49, t48, t218, qJD(3) * t9 - qJD(4) * t13 - t366 + t50, qJD(3) * t10 - qJD(4) * t15 + t367 + t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t17, t73, t71, t219, t392, t393; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t230 * t373 - t361, t233 * t373 - t357; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t8;
