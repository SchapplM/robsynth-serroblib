% Calculate inertial parameters regressor of coriolis matrix for
% S6PRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6PRPRPR7_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:53:54
% EndTime: 2019-03-08 19:54:02
% DurationCPUTime: 5.21s
% Computational Cost: add. (3122->356), mult. (6962->506), div. (0->0), fcn. (6845->8), ass. (0->278)
t266 = sin(pkin(6));
t270 = sin(qJ(2));
t394 = t266 * t270;
t273 = cos(qJ(2));
t393 = t266 * t273;
t269 = sin(qJ(4));
t263 = t269 ^ 2;
t272 = cos(qJ(4));
t265 = t272 ^ 2;
t329 = t263 / 0.2e1 + t265 / 0.2e1;
t334 = t394 / 0.2e1;
t136 = t329 * t394 + t334;
t267 = cos(pkin(6));
t178 = t267 * t269 + t272 * t393;
t347 = t269 * t393;
t392 = t267 * t272;
t179 = -t347 + t392;
t277 = (-t178 * t272 + t179 * t269 + t393) * t394;
t428 = t277 * qJD(1);
t438 = t136 * qJD(3) + t428;
t135 = (-0.1e1 / 0.2e1 + t329) * t394;
t437 = -t135 * qJD(3) - t428;
t433 = t265 + t263;
t271 = cos(qJ(6));
t264 = t271 ^ 2;
t268 = sin(qJ(6));
t262 = t268 ^ 2;
t423 = -t262 / 0.2e1;
t330 = t423 - t264 / 0.2e1;
t436 = 0.1e1 / 0.2e1 + t330;
t250 = t272 * qJ(5);
t212 = -t269 * pkin(4) + t250;
t374 = qJD(5) * t269;
t301 = t212 * qJD(4) + t374;
t275 = -pkin(2) - pkin(8);
t244 = t272 * t275;
t412 = t272 * pkin(5);
t328 = -t244 + t412;
t435 = t268 * t328;
t312 = t271 * t328;
t217 = t271 * t394;
t378 = t178 * t268 + t217;
t415 = t272 / 0.2e1;
t434 = t378 * t415;
t259 = qJD(2) * qJ(3);
t365 = t135 * qJD(1);
t432 = t365 - t259;
t203 = qJ(3) - t212;
t356 = t203 * qJD(2);
t431 = t365 - t356;
t182 = t269 * pkin(9) + t203;
t115 = t268 * t182 - t312;
t116 = t271 * t182 + t435;
t385 = t271 * t272;
t388 = t268 * t272;
t34 = t115 * t388 + t116 * t385;
t384 = t273 * t268;
t154 = (-t270 * t385 - t384) * t266;
t383 = t273 * t271;
t155 = (-t270 * t388 + t383) * t266;
t389 = t268 * t270;
t349 = t266 * t389;
t400 = t178 * t271;
t126 = t349 - t400;
t404 = t126 * t272;
t20 = (t434 + t154 / 0.2e1) * t271 + (t404 / 0.2e1 + t155 / 0.2e1) * t268;
t399 = t20 * qJD(1);
t430 = -t34 * qJD(2) - t399;
t204 = (-pkin(5) + t275) * t269;
t186 = t204 * t271;
t251 = qJ(5) * t269;
t253 = pkin(4) * t272;
t213 = t251 + t253;
t187 = t272 * pkin(9) + t213;
t391 = t268 * t187;
t117 = t186 - t391;
t386 = t271 * t187;
t397 = t204 * t268;
t118 = t386 + t397;
t9 = (-t117 * t268 + t118 * t271) * t269 + t34;
t429 = -t9 * qJD(2) - t399;
t427 = t277 * qJD(2);
t224 = t262 - t264;
t387 = t269 * t271;
t324 = 0.2e1 * t268 * t387;
t287 = qJD(2) * t324 - t224 * qJD(4);
t426 = -t126 / 0.2e1;
t425 = t178 / 0.2e1;
t424 = t179 / 0.2e1;
t421 = t264 / 0.2e1;
t420 = t268 / 0.2e1;
t419 = t269 / 0.2e1;
t418 = -t271 / 0.2e1;
t417 = t271 / 0.2e1;
t416 = -t272 / 0.2e1;
t313 = t328 * t269;
t406 = t118 * t272;
t409 = t116 * t269;
t6 = t313 / 0.2e1 + t204 * t416 + (t115 * t419 + t117 * t415) * t271 + (t406 / 0.2e1 - t409 / 0.2e1) * t268;
t414 = t6 * qJD(4);
t327 = t378 * t268;
t405 = t126 * t271;
t8 = (-t327 / 0.2e1 + t405 / 0.2e1 + t425) * t269 - t436 * t272 * t179;
t413 = t8 * qJD(4);
t410 = t116 * t268;
t408 = t117 * t271;
t407 = t118 * t268;
t403 = t155 * t272;
t350 = t269 * t394;
t17 = -t126 * t154 + t378 * t155 + t179 * t350;
t402 = t17 * qJD(1);
t401 = t178 * qJ(5);
t398 = t20 * qJD(2);
t326 = t378 * t271;
t284 = t326 / 0.2e1 + t126 * t420;
t298 = t154 * t417 + t155 * t420;
t21 = t272 * t284 - t298;
t396 = t21 * qJD(2);
t395 = t263 * t271;
t390 = t268 * t269;
t379 = t433 * t275 * t394;
t377 = t262 + t264;
t225 = t263 - t265;
t376 = qJD(2) * t266;
t375 = qJD(3) * t272;
t373 = qJD(5) * t271;
t372 = qJD(5) * t272;
t371 = qJD(6) * t268;
t370 = qJD(6) * t271;
t369 = qJD(6) * t272;
t274 = -pkin(4) - pkin(9);
t368 = qJD(6) * t274;
t133 = t203 * t272 + t213 * t269;
t367 = t133 * qJD(2);
t134 = -t203 * t269 + t213 * t272;
t366 = t134 * qJD(2);
t197 = t377 * t269;
t165 = t272 * t197;
t364 = t165 * qJD(2);
t317 = -0.1e1 / 0.2e1 - t329;
t166 = t317 * t268;
t363 = t166 * qJD(2);
t169 = t317 * t271;
t362 = t169 * qJD(2);
t361 = t179 * qJD(4);
t183 = (t421 + t423) * t269;
t360 = t183 * qJD(6);
t195 = t225 * t268;
t359 = t195 * qJD(2);
t358 = t197 * qJD(2);
t198 = t225 * t271;
t357 = t198 * qJD(2);
t355 = t225 * qJD(2);
t354 = t265 * qJD(2);
t353 = t268 * qJD(4);
t247 = t269 * qJD(2);
t246 = t269 * qJD(4);
t352 = t269 * qJD(6);
t351 = t271 * qJD(4);
t249 = t272 * qJD(2);
t248 = t272 * qJD(4);
t346 = qJ(3) * t247;
t345 = qJ(3) * t249;
t344 = t268 * t369;
t343 = t271 * t369;
t342 = t203 * t249;
t222 = t270 * t376;
t341 = t273 * t376;
t340 = t269 * t351;
t339 = t268 * t370;
t338 = t268 * t351;
t238 = t269 * t248;
t237 = t269 * t249;
t337 = t275 * t248;
t336 = t251 / 0.2e1;
t335 = t179 * t418;
t333 = -t390 / 0.2e1;
t332 = t389 / 0.2e1;
t331 = -t385 / 0.2e1;
t325 = t377 * t274;
t323 = -qJD(3) + t372;
t322 = t263 * t339;
t321 = t271 * t237;
t320 = t263 * t334;
t319 = t266 * t332;
t318 = t271 * t334;
t316 = (t275 / 0.2e1 - pkin(5) / 0.2e1) * t269;
t314 = qJD(4) * t324;
t16 = (-t178 + t327 - t405) * t179;
t311 = t16 * qJD(1) - t8 * qJD(3);
t310 = t274 * t269 + t250;
t14 = t154 * t331 + t320 - t326 / 0.2e1 + (-t403 / 0.2e1 + t426) * t268;
t45 = t115 * t268 + t116 * t271;
t309 = -t14 * qJD(1) + t45 * qJD(2);
t23 = (t117 - t186) * t272 + (t115 + t312) * t269;
t51 = (t319 + t426 - t400 / 0.2e1) * t269;
t308 = -t51 * qJD(1) + t23 * qJD(2);
t24 = -t204 * t388 + t268 * t313 + t406 - t409;
t52 = (t318 - t217 / 0.2e1) * t269;
t307 = -t52 * qJD(1) - t24 * qJD(2);
t280 = t179 * t333 + t434;
t282 = (-t384 / 0.2e1 + t270 * t331) * t266;
t35 = t282 + t280;
t69 = -t116 * t272 + t204 * t390;
t306 = -t35 * qJD(1) + t69 * qJD(2);
t281 = (-t383 / 0.2e1 + t272 * t332) * t266;
t292 = -t404 / 0.2e1 + t269 * t335;
t36 = t281 + t292;
t68 = -t115 * t272 - t204 * t387;
t305 = -t36 * qJD(1) - t68 * qJD(2);
t304 = t407 + t408;
t303 = -t354 - t369;
t302 = t265 * qJD(5) - t375;
t300 = t336 + t253 / 0.2e1;
t299 = -t408 / 0.2e1 - t407 / 0.2e1;
t297 = t274 * t416 + t336;
t27 = (-t213 / 0.2e1 + t300) * t394;
t296 = t27 * qJD(1) - t213 * t356;
t295 = (-qJD(6) - t249) * t390;
t283 = t187 / 0.2e1 + t297;
t106 = t283 * t268;
t294 = -qJ(5) * t351 - t106 * qJD(2);
t107 = t283 * t271;
t293 = qJ(5) * t353 - t107 * qJD(2);
t142 = -t183 * qJD(2) + t338;
t291 = t392 / 0.2e1 - t347 / 0.2e1;
t290 = -t178 * qJD(4) + t222 * t269;
t289 = t222 * t272 - t361;
t137 = t268 * qJD(2) * t395 + t183 * qJD(4);
t196 = t224 * t263;
t288 = t196 * qJD(2) + t314;
t276 = t334 * t251 + t274 * t298;
t278 = t126 * t117 / 0.2e1 + t204 * t425 - t378 * t118 / 0.2e1;
t1 = (-t410 / 0.2e1 + t115 * t417 + t412 / 0.2e1 - t244 / 0.2e1) * t179 + t276 + t278;
t10 = -t115 * t117 + t116 * t118 - t204 * t328;
t286 = -t1 * qJD(1) + t10 * qJD(2) - t6 * qJD(3);
t143 = -t272 * t269 + t165;
t285 = -t8 * qJD(1) - t6 * qJD(2) - t143 * qJD(3);
t168 = t436 * t269;
t258 = qJD(4) * qJ(5);
t26 = t316 + t299;
t67 = t179 * t330 + t291;
t279 = t67 * qJD(1) + t26 * qJD(2) + t168 * qJD(3) + t258;
t261 = qJ(3) * qJD(3);
t260 = qJ(5) * qJD(5);
t242 = -t246 / 0.2e1;
t241 = t275 * t246;
t239 = t271 * t352;
t236 = t271 * t248;
t235 = t271 * t249;
t234 = t268 * t248;
t233 = t268 * t246;
t232 = t268 * t249;
t223 = qJ(3) * t393;
t214 = t225 * qJD(4);
t194 = -t235 - t370;
t193 = -t232 - t371;
t188 = t237 + t352 / 0.2e1;
t171 = t265 * t417 + t395 / 0.2e1 + t418;
t170 = -t268 / 0.2e1 + t433 * t420;
t167 = t262 * t419 + (0.1e1 / 0.2e1 + t421) * t269;
t161 = t433 * t222;
t153 = (-t246 * t270 + t249 * t273) * t266;
t152 = (t247 * t273 + t248 * t270) * t266;
t123 = t136 * qJD(2);
t120 = t135 * qJD(2);
t73 = -t397 - t386 / 0.2e1 + t297 * t271;
t72 = t186 - t391 / 0.2e1 + t297 * t268;
t66 = t179 * t421 + t262 * t424 + t291;
t54 = t178 * t333 + t269 * t318 + t378 * t419;
t53 = t126 * t419 + t269 * t319 + t387 * t425;
t38 = t282 - t280;
t37 = t281 - t292;
t28 = t213 * t334 + t300 * t394;
t25 = t316 - t299;
t15 = -t272 * t298 + t284 + t320;
t2 = t410 * t424 + t115 * t335 - t179 * t328 / 0.2e1 + t276 - t278;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t427, 0, 0, 0, 0, 0, 0, 0, 0, 0, t427, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 * qJD(2) + t16 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t222, -t341, 0, 0, 0, 0, 0, 0, 0, 0, 0, t222, t341 (-pkin(2) * t394 + t223) * qJD(2) + qJD(3) * t394, 0, 0, 0, 0, 0, 0, t152, t153, -t161 (t223 + t379) * qJD(2) + t438, 0, 0, 0, 0, 0, 0, -t161, -t152, -t153 (t203 * t393 + t379) * qJD(2) + t28 * qJD(4) - t372 * t394 + t438, 0, 0, 0, 0, 0, 0 (t154 * t272 - t217 * t263) * qJD(2) + t53 * qJD(4) + t38 * qJD(6) (t263 * t349 - t403) * qJD(2) + t54 * qJD(4) + t37 * qJD(6), t21 * qJD(4) + (-t154 * t268 + t155 * t271) * t247, t402 + (-t154 * t115 + t155 * t116 + t204 * t350) * qJD(2) + t15 * qJD(3) + t2 * qJD(4) - t21 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t222, 0, 0, 0, 0, 0, 0, 0, 0, 0, t123, 0, 0, 0, 0, 0, 0, 0, 0, 0, t123, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 * qJD(2) - t413; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t289, -t290, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t289, t290, t28 * qJD(2) + (-t179 * pkin(4) - t401) * qJD(4) + t179 * qJD(5), 0, 0, 0, 0, 0, 0, t53 * qJD(2) - t178 * t353 + t179 * t370, t54 * qJD(2) - t178 * t351 - t179 * t371, -t361 * t377 + t396, t2 * qJD(2) + (t179 * t325 - t401) * qJD(4) + t66 * qJD(5) + t311; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t289, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66 * qJD(4) - t396; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38 * qJD(2) - qJD(6) * t378 + t179 * t351, t37 * qJD(2) + t126 * qJD(6) - t179 * t353, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t437, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27 * qJD(4) + t437, 0, 0, 0, 0, 0, 0, -t51 * qJD(4) - t35 * qJD(6), -t52 * qJD(4) - t36 * qJD(6), t20 * qJD(4), -t14 * qJD(3) - t1 * qJD(4) - t20 * qJD(5) - t402; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t261, -t238, t214, 0, t238, 0, 0, qJ(3) * t248 + qJD(3) * t269, -qJ(3) * t246 + t375, 0, t261, 0, 0, 0, -t238, t214, t238, 0, -t133 * qJD(4) + t269 * t323, -t134 * qJD(4) + t302 (qJD(4) * t213 - t323) * t203, t238 * t262 + t322, -t196 * qJD(6) + t272 * t314, -t195 * qJD(4) + t269 * t343, t238 * t264 - t322, -t198 * qJD(4) - t269 * t344, -t238, t23 * qJD(4) + t69 * qJD(6) + t268 * t302, -t24 * qJD(4) - t68 * qJD(6) + t271 * t302, t197 * qJD(3) + t9 * qJD(4) - t165 * qJD(5), t45 * qJD(3) + t10 * qJD(4) - t34 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t259, 0, 0, 0, 0, 0, 0, t247, t249, 0, -t432, 0, 0, 0, 0, 0, 0, 0, -t247, -t249, -t431, 0, 0, 0, 0, 0, 0, t170 * qJD(6) - t232, t171 * qJD(6) - t235, t358, t309 - t414; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t237, t355, -t246, t237, -t248, 0, -t241 + t345, -t337 - t346, 0, 0, 0, t246, t248, -t237, t355, t237, -t301, t241 - t367, t337 - t366, t275 * t301 - t296, t360 + (t247 * t262 + t338) * t272, -0.2e1 * t269 * t339 + t272 * t287, -t340 - t359, -t360 + (t247 * t264 - t338) * t272, t233 - t357, -t188 (-t271 * t310 - t435) * qJD(4) - t269 * t373 + t72 * qJD(6) + t308 (t268 * t310 - t312) * qJD(4) + t268 * t374 + t73 * qJD(6) + t307, -qJD(4) * t304 - t429 (-qJ(5) * t328 + t274 * t304) * qJD(4) + t25 * qJD(5) + t286; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t246, t237, t354, t241 - t342, 0, 0, 0, 0, 0, 0, t268 * t354 - t340, t271 * t354 + t233, -t364, t25 * qJD(4) + t430; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137, -t288, t239 + t321, -t137, t295, t242, t170 * qJD(3) + t72 * qJD(4) - t116 * qJD(6) + t306, t171 * qJD(3) + t73 * qJD(4) + t115 * qJD(6) + t305, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 * qJD(2) - t413; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), -t259, 0, 0, 0, 0, 0, 0, -t247, -t249, 0, t432, 0, 0, 0, 0, 0, 0, 0, t247, t249, t431, 0, 0, 0, 0, 0, 0, -t166 * qJD(6) + t232, -t169 * qJD(6) + t235, -t358, -t309 - t414; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t246, -t248, 0, 0, 0, 0, 0, 0, 0, 0, 0, t246, t248, t301, 0, 0, 0, 0, 0, 0, t239 + t234, -t268 * t352 + t236, -t197 * qJD(4) (t269 * t325 + t250) * qJD(4) + t167 * qJD(5) + t285; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t246, 0, 0, 0, 0, 0, 0, 0, 0, 0, t167 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t340 + t344 - t363, -t233 + t343 - t362, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27 * qJD(2), 0, 0, 0, 0, 0, 0, t51 * qJD(2), t52 * qJD(2), -t398, t1 * qJD(2) + t67 * qJD(5) - t311; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t237, -t355, 0, -t237, 0, 0, -t345, t346, 0, 0, 0, 0, 0, t237, -t355, -t237, 0, t367, t366, t296, -t237 * t262 + t360, 0.2e1 * t271 * t295, -t344 + t359, -t237 * t264 - t360, -t343 + t357, t188, t106 * qJD(6) - t308, t107 * qJD(6) - t307, t429, t26 * qJD(5) - t286; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t168 * qJD(5) - t285; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), t260, -t339, t224 * qJD(6), 0, t339, 0, 0, qJ(5) * t370 + qJD(5) * t268, -qJ(5) * t371 + t373, 0, t260; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t258, 0, 0, 0, 0, 0, 0, t353, t351, 0, t279; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t142, -t287, t193, t142, t194, t247 / 0.2e1, -t268 * t368 - t294, -t271 * t368 - t293, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67 * qJD(4) + t398; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t237, -t354, t342, 0, 0, 0, 0, 0, 0, t303 * t268, t303 * t271, t364, -t26 * qJD(4) - t430; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t168 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4), -t258, 0, 0, 0, 0, 0, 0, -t353, -t351, 0, -t279; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, t194, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35 * qJD(2), t36 * qJD(2), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t137, t288, t234 - t321, t137, t237 * t268 + t236, t242, t166 * qJD(3) - t106 * qJD(4) + t268 * t372 - t306, t169 * qJD(3) - t107 * qJD(4) + t271 * t372 - t305, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t363, t362, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t142, t287, t232, -t142, t235, -t247 / 0.2e1, t294, t293, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t232, t235, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t3;
