% Calculate minimal parameter regressor of coriolis matrix for
% S6RPRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x29]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPRPRR4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:46:31
% EndTime: 2019-03-09 03:46:42
% DurationCPUTime: 4.65s
% Computational Cost: add. (3280->358), mult. (6784->490), div. (0->0), fcn. (6687->8), ass. (0->312)
t319 = qJD(5) + qJD(6);
t244 = sin(qJ(5));
t429 = 0.2e1 * t244;
t246 = cos(qJ(6));
t247 = cos(qJ(5));
t358 = t246 * t247;
t243 = sin(qJ(6));
t367 = t243 * t244;
t183 = -t358 + t367;
t174 = t183 * qJD(6);
t245 = sin(qJ(3));
t371 = t183 * t245;
t296 = -t371 / 0.2e1;
t207 = t245 * t358;
t295 = -t367 / 0.2e1;
t348 = t245 * t295 + t207 / 0.2e1;
t281 = t296 + t348;
t416 = t281 * qJD(1);
t428 = qJD(5) * t183 + t174 - t416;
t249 = -pkin(3) - pkin(8);
t426 = -pkin(9) + t249;
t190 = t426 * t244;
t191 = t426 * t247;
t133 = t190 * t246 + t191 * t243;
t427 = t319 * t133;
t359 = t246 * t244;
t366 = t243 * t247;
t184 = t359 + t366;
t425 = t184 * t319;
t248 = cos(qJ(3));
t424 = t248 * t249;
t423 = t281 * t319;
t363 = t244 * t245;
t294 = t363 / 0.2e1;
t114 = t296 + t243 * t294 - t207 / 0.2e1;
t422 = t319 * t114;
t226 = sin(pkin(10)) * pkin(1) + pkin(7);
t401 = pkin(4) + t226;
t404 = pkin(5) * t247;
t287 = t401 + t404;
t272 = t287 * t248;
t261 = -t272 / 0.2e1;
t357 = t247 * t248;
t208 = t246 * t357;
t362 = t244 * t248;
t155 = t243 * t362 - t208;
t405 = pkin(5) * t244;
t228 = qJ(4) + t405;
t369 = t228 * t155;
t421 = t184 * t261 + t369 / 0.2e1;
t157 = t248 * t184;
t368 = t228 * t157;
t420 = t183 * t261 - t368 / 0.2e1;
t177 = t401 * t248;
t159 = t177 * t244;
t356 = t248 * qJ(4);
t204 = -pkin(3) * t245 + t356;
t179 = pkin(8) * t245 - t204;
t167 = t247 * t179;
t350 = t167 + t159;
t360 = t245 * t247;
t81 = pkin(9) * t360 + t350;
t387 = t243 * t81;
t315 = -t387 / 0.2e1;
t160 = t177 * t247;
t402 = t248 * pkin(5);
t75 = t402 + t160 + (-pkin(9) * t245 - t179) * t244;
t384 = t246 * t75;
t269 = t315 + t384 / 0.2e1;
t253 = t368 / 0.2e1 + t269;
t260 = t272 / 0.2e1;
t372 = t133 * t245;
t297 = t372 / 0.2e1;
t11 = t183 * t260 + t253 + t297;
t329 = t114 * qJD(2);
t339 = qJD(3) * t228;
t419 = qJD(1) * t11 + t183 * t339 - t329;
t104 = -t183 * t228 + t184 * t404;
t318 = t404 / 0.2e1;
t251 = t155 * t318 + t297;
t409 = t184 / 0.2e1;
t410 = t183 / 0.2e1;
t4 = (t246 * pkin(5) / 0.2e1 + t405 * t409 + t287 * t410) * t248 + t251 + t253;
t418 = qJD(1) * t4 - qJD(3) * t104 - t329;
t417 = qJD(3) * t281;
t238 = t244 ^ 2;
t240 = t247 ^ 2;
t213 = t238 - t240;
t322 = t248 * qJD(1);
t311 = t247 * t322;
t257 = qJD(3) * t213 + t311 * t429;
t321 = t248 * qJD(4);
t361 = t245 * qJ(4);
t415 = qJD(3) * (t361 - t424) - t321;
t277 = t190 * t243 - t191 * t246;
t268 = -t359 / 0.2e1 - t366 / 0.2e1;
t255 = t409 + t268;
t110 = t255 * t245;
t330 = t110 * qJD(2);
t414 = t277 * t319 - t330;
t403 = t245 * pkin(5);
t227 = -cos(pkin(10)) * pkin(1) - pkin(2);
t271 = t227 - t361;
t150 = t271 + t424;
t176 = t401 * t245;
t158 = t247 * t176;
t82 = t150 * t244 - t158;
t69 = pkin(9) * t362 - t82;
t66 = t69 + t403;
t57 = t246 * t66;
t413 = -t57 / 0.2e1;
t412 = -t66 / 0.2e1;
t411 = -t183 / 0.2e1;
t408 = t208 / 0.2e1;
t407 = -t245 / 0.2e1;
t406 = pkin(3) * t248;
t37 = -t155 * t184 - t157 * t183;
t400 = t319 * t37;
t51 = t155 * t411 + t157 * t409;
t399 = t319 * t51;
t398 = t110 * qJD(4);
t109 = t246 * t294 + t243 * t360 / 0.2e1 + t245 * t409;
t397 = t109 * qJD(4);
t396 = t319 * t110;
t395 = t319 * t109;
t394 = pkin(5) * qJD(5);
t393 = pkin(5) * qJD(6);
t144 = t287 * t245;
t154 = t243 * t363 - t207;
t365 = t244 * t176;
t83 = t150 * t247 + t365;
t70 = -pkin(9) * t357 + t83;
t389 = t243 * t70;
t27 = -t57 + t389;
t1 = (t384 - t387) * t245 + t144 * t155 + (t154 * t287 - t27) * t248;
t392 = t1 * qJD(1);
t156 = t184 * t245;
t385 = t246 * t70;
t28 = t243 * t66 + t385;
t383 = t246 * t81;
t388 = t243 * t75;
t2 = (t383 + t388) * t245 + t28 * t248 - t144 * t157 - t156 * t272;
t391 = t2 * qJD(1);
t390 = t243 * t69;
t386 = t246 * t69;
t317 = -t403 / 0.2e1;
t283 = t317 + t69 / 0.2e1;
t7 = (t412 + t283) * t243;
t382 = t7 * qJD(1);
t9 = t246 * t283 + t413;
t381 = t9 * qJD(1);
t31 = -t385 - t390;
t16 = t31 * t245 + (t155 * t405 - t157 * t287) * t248;
t380 = qJD(1) * t16;
t32 = t386 - t389;
t86 = t155 * t272;
t17 = -pkin(5) * t157 * t362 + t245 * t32 - t86;
t379 = qJD(1) * t17;
t20 = -t245 * t27 - t86;
t378 = qJD(1) * t20;
t21 = -t157 * t272 - t28 * t245;
t377 = qJD(1) * t21;
t49 = t177 * t357 - t245 * t82;
t376 = qJD(1) * t49;
t50 = -t177 * t362 - t245 * t83;
t375 = qJD(1) * t50;
t80 = t156 * t245 - t157 * t248;
t374 = qJD(1) * t80;
t373 = t277 * t245;
t370 = t184 * t228;
t364 = t244 * t179;
t25 = (-t82 - t158) * t248 - t364 * t245;
t355 = t25 * qJD(1);
t26 = t350 * t245 - t177 * t363 + (t83 - t365) * t248;
t354 = t26 * qJD(1);
t43 = t154 * t157 + t155 * t156;
t353 = t43 * qJD(1);
t79 = -t154 * t245 + t155 * t248;
t352 = t79 * qJD(1);
t351 = t319 * t157;
t239 = t245 ^ 2;
t241 = t248 ^ 2;
t214 = t239 - t241;
t347 = qJD(1) * t157;
t186 = t214 * t244;
t346 = qJD(1) * t186;
t188 = t214 * t247;
t345 = qJD(1) * t188;
t344 = qJD(1) * t245;
t343 = qJD(1) * t247;
t342 = qJD(3) * qJ(4);
t341 = qJD(3) * t183;
t340 = qJD(3) * t184;
t338 = qJD(4) * t244;
t337 = qJD(4) * t245;
t336 = qJD(4) * t247;
t335 = qJD(5) * t244;
t334 = qJD(5) * t245;
t333 = qJD(5) * t247;
t332 = qJD(5) * t249;
t90 = t109 * qJD(1);
t331 = t109 * qJD(3);
t91 = t110 * qJD(3);
t173 = t271 - t406;
t128 = t173 * t248 - t204 * t245;
t328 = t128 * qJD(1);
t129 = -t173 * t245 - t204 * t248;
t327 = t129 * qJD(1);
t326 = t214 * qJD(1);
t325 = t239 * qJD(1);
t324 = t244 * qJD(3);
t234 = t245 * qJD(3);
t323 = t247 * qJD(3);
t235 = t248 * qJD(3);
t320 = t248 * qJD(5);
t316 = t402 / 0.2e1;
t314 = -t383 / 0.2e1;
t313 = t156 * t344;
t312 = t157 * t344;
t310 = t226 * t234;
t309 = t244 * t323;
t308 = t244 * t235;
t307 = t184 * t235;
t306 = t244 * t334;
t223 = t244 * t320;
t305 = t245 * t333;
t304 = t247 * t320;
t303 = t154 * t344;
t302 = t173 * t204 * qJD(1);
t301 = t173 * t344;
t300 = t227 * t344;
t299 = t227 * t322;
t218 = t245 * t235;
t221 = t247 * t235;
t298 = t244 * t333;
t222 = t245 * t322;
t293 = -t157 / 0.2e1;
t292 = pkin(5) * t319;
t291 = t357 * t429;
t290 = t319 * t155;
t289 = t319 * t245;
t288 = qJD(5) + t344;
t286 = t244 * t222;
t285 = t244 * t221;
t282 = t316 + t75 / 0.2e1;
t105 = -t183 * t404 - t370;
t250 = t157 * t318 + t277 * t407;
t270 = -t388 / 0.2e1 + t314;
t254 = -t369 / 0.2e1 + t270;
t3 = (-t243 * pkin(5) / 0.2e1 + t405 * t411 + t287 * t409) * t248 + t250 + t254;
t278 = -qJD(1) * t3 + qJD(3) * t105;
t56 = t155 ^ 2 - t157 ^ 2;
t18 = qJD(1) * t56 + qJD(3) * t37;
t84 = -t183 ^ 2 + t184 ^ 2;
t22 = qJD(1) * t37 + qJD(3) * t84;
t108 = t408 + (t295 + t410) * t248;
t276 = qJD(1) * t108 + t340;
t111 = t255 * t248;
t275 = qJD(1) * t111 - t341;
t274 = -t325 - t334;
t273 = qJD(3) * t204 + t337;
t267 = t249 * t407 - t356 / 0.2e1;
t12 = -t373 / 0.2e1 + t184 * t260 + t254;
t265 = qJD(1) * t12 + t184 * t339;
t30 = -qJD(3) * t51 + t155 * t347;
t38 = qJD(1) * t51 + t183 * t340;
t118 = (t179 / 0.2e1 + t267) * t244;
t264 = -qJ(4) * t323 - qJD(1) * t118;
t256 = t267 * t247;
t119 = t167 / 0.2e1 + t256;
t263 = qJ(4) * t324 - qJD(1) * t119;
t178 = (t240 / 0.2e1 - t238 / 0.2e1) * t248;
t262 = qJD(1) * t178 + t309;
t259 = t241 * t244 * t343 - qJD(3) * t178;
t187 = t213 * t241;
t258 = -qJD(1) * t187 + 0.2e1 * t285;
t252 = (-t361 - t406) * qJD(3) + t321;
t232 = -t322 / 0.2e1;
t231 = t322 / 0.2e1;
t230 = t235 / 0.2e1;
t220 = t245 * t323;
t219 = t245 * t343;
t217 = t244 * t344;
t189 = t226 * t235;
t182 = -t219 - t333;
t181 = -t217 - t335;
t180 = t222 + t320 / 0.2e1;
t168 = t178 * qJD(5);
t166 = t183 * t235;
t151 = t222 + (qJD(5) / 0.2e1 + qJD(6) / 0.2e1) * t248;
t113 = t371 / 0.2e1 + t348;
t112 = t248 * t268 + t293;
t107 = t408 + (t295 + t411) * t248;
t102 = t114 * qJD(4);
t100 = t281 * qJD(4);
t59 = -t159 - t167 / 0.2e1 + t256;
t58 = t160 - t364 / 0.2e1 + t267 * t244;
t53 = -t312 + t417;
t52 = -t155 * t344 + t331;
t45 = -t90 - t425;
t42 = -t290 - t331;
t41 = t114 * qJD(3);
t34 = qJD(3) * t113 + t312 + t351;
t33 = -t91 + (qJD(6) + t288) * t155;
t24 = t351 + t417;
t14 = -t372 / 0.2e1 + t269 + t420;
t13 = t373 / 0.2e1 + t270 + t421;
t10 = t246 * t317 + t389 + t413 - t386 / 0.2e1;
t8 = -t385 - t390 / 0.2e1 + (t317 + t412) * t243;
t6 = t183 * t244 * t316 - t243 * t282 - t250 + t314 + t421;
t5 = t246 * t282 + t293 * t405 - t251 + t315 + t420;
t15 = [0, 0, 0, 0, t218, -t214 * qJD(3), 0, 0, 0, t227 * t234, t227 * t235, 0, qJD(3) * t129 - t245 * t321, -qJD(3) * t128 + qJD(4) * t239, -t273 * t173, -t218 * t238 + t241 * t298, -qJD(5) * t187 - 0.2e1 * t245 * t285, qJD(3) * t186 - t245 * t304, qJD(3) * t188 + t223 * t245, t218, qJD(3) * t25 + qJD(5) * t50 + t239 * t338, -qJD(3) * t26 - qJD(5) * t49 + t239 * t336 (-qJD(3) * t156 - t290) * t157, t43 * qJD(3) + t319 * t56, t80 * qJD(3) + t155 * t289, t79 * qJD(3) + t157 * t289, t218, qJD(3) * t1 + qJD(5) * t16 + qJD(6) * t21 + t156 * t337, -qJD(3) * t2 - qJD(5) * t17 - qJD(6) * t20 - t154 * t337; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t222, -t326, t235, -t234, 0, -t189 + t300, t299 + t310, t252, t189 + t327, -t310 - t328, t226 * t252 - t302, -t168 + (-t238 * t322 + t309) * t245, qJD(5) * t291 - t245 * t257, t221 + t346, -t308 + t345, t180, t58 * qJD(5) - t176 * t324 - t247 * t415 + t355, t59 * qJD(5) - t176 * t323 + t244 * t415 - t354 (-t341 - t347) * t156 + t399, t353 + (t154 * t183 - t156 * t184) * qJD(3) + t400, -t166 + t374 - t396, t113 * t319 - t307 + t352, t151, t392 + (-t144 * t184 + t154 * t228 - t248 * t277) * qJD(3) + t107 * qJD(4) + t5 * qJD(5) + t14 * qJD(6), -t391 + (-t133 * t248 + t144 * t183 + t156 * t228) * qJD(3) + t112 * qJD(4) + t6 * qJD(5) + t13 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t235, -t222, t325, t189 - t301, 0, 0, 0, 0, 0, t244 * t325 + t221, t247 * t325 - t308, 0, 0, 0, 0, 0, qJD(3) * t107 + t313 - t396, t112 * qJD(3) - t303 - t422; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t259, t258, -t288 * t357, t223 + t286, t230, qJD(3) * t58 - qJD(5) * t83 + t375, qJD(3) * t59 + qJD(5) * t82 - t376, -t30, t18, t33, t34, t230, qJD(3) * t5 + qJD(5) * t31 + qJD(6) * t8 + t380 - t398, qJD(3) * t6 - qJD(5) * t32 + qJD(6) * t10 - t102 - t379; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t18, t33, t34, t230, qJD(3) * t14 + qJD(5) * t8 - qJD(6) * t28 + t377 - t398, qJD(3) * t13 + qJD(5) * t10 + qJD(6) * t27 - t102 - t378; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t234, -t235, 0, t234, t235, t273, 0, 0, 0, 0, 0, t305 + t308, t221 - t306, 0, 0, 0, 0, 0, t307 + t423, -t166 - t395; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t234, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t223 + t220, -t234 * t244 + t304, 0, 0, 0, 0, 0, t24, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t42; 0, 0, 0, 0, -t222, t326, 0, 0, 0, -t300, -t299, 0, -t327, t328, t302, t222 * t238 - t168, t288 * t291, -t306 - t346, -t305 - t345, -t180, qJD(5) * t118 - t355, qJD(5) * t119 + t354, t156 * t347 + t399, -t353 + t400, -t374 - t395, -t352 - t423, -t151, qJD(4) * t108 - qJD(5) * t4 - qJD(6) * t11 - t392, qJD(4) * t111 - qJD(5) * t3 - qJD(6) * t12 + t391; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t422, -t396; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), qJ(4) * qJD(4), -t298, t213 * qJD(5), 0, 0, 0, qJ(4) * t333 + t338, -qJ(4) * t335 + t336, t183 * t425, t319 * t84, 0, 0, 0, qJD(4) * t184 + qJD(5) * t104 - t174 * t228, -qJD(4) * t183 + qJD(5) * t105 - qJD(6) * t370; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t342, 0, 0, 0, 0, 0, t324, t323, 0, 0, 0, 0, 0, t276, t275; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t262, t257, t181, t182, t232, -t244 * t332 - t264, -t247 * t332 - t263, t38, t22, t45, t428, t232, -t418 - t427, t278 + t414; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t22, t45, t428, t232, -t419 - t427, -t265 + t414; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t222, -t325, t301, 0, 0, 0, 0, 0, t274 * t244, t274 * t247, 0, 0, 0, 0, 0, -qJD(3) * t108 - t313 - t395, -t111 * qJD(3) + t303 - t423; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t342, 0, 0, 0, 0, 0, -t324, -t323, 0, 0, 0, 0, 0, -t276, -t275; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t181, t182, 0, 0, 0, 0, 0, t45, t428; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t428; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t259, -t258 (t311 + t324) * t245, t220 - t286, t230, -qJD(3) * t118 + t244 * t337 - t375, -qJD(3) * t119 + t245 * t336 + t376, t30, -t18, t52, t53, t230, qJD(3) * t4 + qJD(6) * t7 - t380 + t397, qJD(3) * t3 + qJD(6) * t9 + t100 + t379; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t262, -t257, t217, t219, t231, t264, t263, -t38, -t22, t90, t416, t231, t418, -t278 + t330; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t217, t219, 0, 0, 0, 0, 0, t90, t416; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t243 * t393, -t246 * t393; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t243 * t292 + t382, -t246 * t292 + t381; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t18, t52, t53, t230, qJD(3) * t11 - qJD(5) * t7 - t377 + t397, qJD(3) * t12 - qJD(5) * t9 + t100 + t378; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t22, t90, t416, t231, t419, t265 + t330; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, t416; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t243 * t394 - t382, t246 * t394 - t381; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t15;
