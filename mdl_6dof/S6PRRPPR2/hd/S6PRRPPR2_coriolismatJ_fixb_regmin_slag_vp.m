% Calculate minimal parameter regressor of coriolis matrix for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x26]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6PRRPPR2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:21:24
% EndTime: 2021-01-16 02:21:38
% DurationCPUTime: 4.98s
% Computational Cost: add. (3415->345), mult. (8177->509), div. (0->0), fcn. (9365->10), ass. (0->282)
t235 = sin(qJ(3));
t238 = cos(qJ(3));
t233 = sin(pkin(6));
t236 = sin(qJ(2));
t377 = t233 * t236;
t395 = cos(pkin(6));
t185 = t235 * t377 - t238 * t395;
t186 = t235 * t395 + t238 * t377;
t232 = sin(pkin(11));
t394 = cos(pkin(11));
t125 = -t185 * t394 - t186 * t232;
t412 = t125 / 0.2e1;
t442 = -t125 / 0.2e1 + t412;
t223 = t394 * t238;
t379 = t232 * t235;
t199 = -t223 + t379;
t409 = -t199 / 0.2e1;
t172 = t394 * t186;
t381 = t232 * t185;
t429 = t172 - t381;
t310 = t429 * t409;
t301 = t394 * t235;
t378 = t232 * t238;
t201 = t301 + t378;
t391 = t125 * t201;
t246 = t310 - t391 / 0.2e1 + t377 / 0.2e1;
t239 = cos(qJ(2));
t376 = t233 * t239;
t159 = t201 * t376;
t210 = t376 * t379;
t292 = t394 * t376;
t269 = t238 * t292;
t160 = t269 - t210;
t30 = t233 ^ 2 * t236 * t239 + t125 * t159 - t160 * t429;
t437 = t30 * qJD(1);
t441 = t246 * qJD(4) - t437;
t431 = t429 / 0.2e1;
t302 = t431 - t429 / 0.2e1;
t308 = -t377 / 0.2e1;
t262 = t308 + t310;
t245 = t391 / 0.2e1 - t262;
t440 = -qJD(4) * t245 + t437;
t196 = t199 ^ 2;
t414 = t201 ^ 2;
t427 = -t414 - t196;
t439 = t427 * qJD(2);
t438 = t427 * qJD(4);
t436 = t30 * qJD(2);
t237 = cos(qJ(6));
t405 = t237 / 0.2e1;
t309 = t429 * t405;
t289 = t302 * t201;
t244 = t199 * t442 + t289;
t421 = t244 * qJD(3);
t435 = (t159 * t201 - t160 * t199) * qJD(2) + t421;
t234 = sin(qJ(6));
t135 = t234 * t201;
t349 = t135 * qJD(2);
t355 = qJD(6) * t234;
t434 = t349 + t355;
t402 = qJ(4) + pkin(8);
t215 = t402 * t238;
t209 = t394 * t215;
t214 = t402 * t235;
t380 = t232 * t214;
t428 = t209 - t380;
t433 = -t199 * pkin(5) + t428;
t432 = t199 * t412 - t289;
t430 = 0.2e1 * t433;
t155 = -t214 * t394 - t215 * t232;
t273 = -t155 * t201 - t199 * t428;
t423 = qJD(4) * t273;
t422 = t244 * qJD(2);
t420 = t245 * qJD(2);
t419 = t246 * qJD(2);
t230 = t234 ^ 2;
t231 = t237 ^ 2;
t219 = t230 - t231;
t138 = t237 * t199;
t300 = 0.2e1 * t234 * t138;
t255 = qJD(2) * t300 - qJD(3) * t219;
t416 = qJD(1) * t244;
t415 = -qJD(1) * t245 + qJD(2) * t273;
t413 = pkin(4) + pkin(9);
t291 = t172 / 0.2e1;
t408 = -t201 / 0.2e1;
t290 = t209 / 0.2e1;
t225 = -pkin(3) * t394 - pkin(4);
t407 = -t225 / 0.2e1;
t406 = -t234 / 0.2e1;
t403 = t201 * pkin(4);
t228 = t235 * pkin(3);
t401 = qJD(3) * pkin(3);
t385 = t199 * qJ(5);
t299 = t228 + t385;
t90 = t201 * t413 + t299;
t400 = t234 * t90;
t399 = t237 * t90;
t87 = -t125 * t237 + t234 * t376;
t398 = t87 * t201;
t88 = t125 * t234 + t237 * t376;
t397 = t88 * t201;
t190 = t201 * qJD(4);
t253 = t378 / 0.2e1 + t301 / 0.2e1;
t102 = (t408 + t253) * t376;
t353 = t102 * qJD(1);
t396 = t190 + t353;
t119 = -pkin(5) * t201 + t155;
t393 = t119 * t234;
t392 = t119 * t237;
t387 = t160 * t234;
t386 = t160 * t237;
t384 = t199 * t234;
t222 = pkin(3) * t232 + qJ(5);
t383 = t222 * t199;
t382 = t222 * t201;
t375 = t234 * t159;
t374 = t237 * t159;
t226 = -t238 * pkin(3) - pkin(2);
t265 = -t201 * qJ(5) + t226;
t134 = pkin(4) * t199 + t265;
t140 = t299 + t403;
t39 = -t134 * t201 - t140 * t199;
t369 = t39 * qJD(2);
t227 = t228 / 0.2e1;
t47 = t227 + (pkin(4) / 0.2e1 + t407) * t201 + (qJ(5) / 0.2e1 + t222 / 0.2e1) * t199;
t368 = t47 * qJD(2);
t64 = t427 * t234;
t367 = t64 * qJD(2);
t327 = t414 - t196;
t95 = t327 * t234;
t366 = t95 * qJD(2);
t96 = t327 * t237;
t365 = t96 * qJD(2);
t97 = t427 * t237;
t364 = t97 * qJD(2);
t266 = -t292 / 0.2e1;
t362 = t210 / 0.2e1 + t238 * t266;
t361 = -t210 / 0.2e1 + t269 / 0.2e1;
t307 = -t376 / 0.2e1;
t360 = t235 * t266 + t307 * t378;
t359 = qJD(2) * t234;
t358 = qJD(2) * t236;
t357 = qJD(2) * t238;
t356 = qJD(5) * t234;
t354 = qJD(6) * t237;
t324 = -t394 / 0.2e1;
t247 = (t201 * t324 + t232 * t409) * pkin(3);
t116 = -t228 / 0.2e1 + t247;
t352 = t116 * qJD(2);
t351 = t429 * qJD(3);
t350 = t125 * qJD(3);
t348 = t135 * qJD(6);
t347 = t138 * qJD(2);
t346 = t138 * qJD(3);
t343 = t428 * qJD(3);
t342 = t155 * qJD(3);
t193 = t379 / 0.2e1 - t223 / 0.2e1;
t341 = t193 * qJD(2);
t340 = t414 * qJD(2);
t339 = t414 * qJD(5);
t338 = t199 * qJD(2);
t188 = t199 * qJD(3);
t337 = t199 * qJD(4);
t336 = t199 * qJD(5);
t335 = t201 * qJD(2);
t334 = t201 * qJD(3);
t333 = t201 * qJD(5);
t220 = -t235 ^ 2 + t238 ^ 2;
t332 = t220 * qJD(2);
t331 = t234 * qJD(3);
t330 = t235 * qJD(3);
t329 = t237 * qJD(3);
t328 = t238 * qJD(3);
t326 = pkin(2) * t235 * qJD(2);
t325 = pkin(2) * t357;
t323 = t201 * t355;
t322 = t201 * t354;
t321 = t134 * t335;
t320 = t199 * t335;
t319 = t199 * t334;
t318 = t233 * t358;
t317 = qJD(2) * t376;
t316 = t234 * t338;
t315 = t199 * t331;
t314 = t414 * t359;
t313 = t235 * t357;
t312 = t234 * t354;
t181 = t237 * t335;
t311 = t234 * t329;
t306 = t376 / 0.2e1;
t305 = -t375 / 0.2e1;
t304 = t374 / 0.2e1;
t298 = qJD(6) + t335;
t297 = t199 * t318;
t296 = t201 * t318;
t295 = t199 * t307;
t294 = t199 * t306;
t293 = t201 * t306;
t287 = qJD(3) * t300;
t89 = t199 * t413 + t265;
t38 = t237 * t89 - t393;
t4 = -t399 * t201 + (t38 + t393) * t199;
t241 = t234 * t432 + t88 * t409;
t8 = -t386 / 0.2e1 + t241;
t286 = t8 * qJD(1) + t4 * qJD(2);
t242 = -t237 * t432 + t87 * t409;
t10 = -t387 / 0.2e1 + t242;
t37 = t234 * t89 + t392;
t3 = -t400 * t201 + (t37 - t392) * t199;
t285 = t10 * qJD(1) + t3 * qJD(2);
t19 = t134 * t140;
t243 = t155 * t302 + t428 * t442;
t240 = t140 * t307 - t243;
t264 = t159 * t407 - t160 * t222 / 0.2e1;
t2 = t240 + t264;
t284 = qJD(1) * t2 + qJD(2) * t19;
t32 = t226 * t228;
t254 = t160 * t232 / 0.2e1 + t159 * t324;
t5 = (t235 * t306 + t254) * pkin(3) + t243;
t283 = -t5 * qJD(1) + t32 * qJD(2);
t20 = t138 * t433 + t201 * t37;
t25 = t305 + t398 / 0.2e1 + t262 * t237;
t280 = -qJD(1) * t25 + qJD(2) * t20;
t21 = -t201 * t38 + t384 * t433;
t24 = t304 - t397 / 0.2e1 + t262 * t234;
t279 = -qJD(1) * t24 + qJD(2) * t21;
t274 = -t155 * t159 + t160 * t428;
t221 = -pkin(9) + t225;
t272 = t199 * t221 + t382;
t104 = t294 + t361;
t40 = t134 * t199 - t140 * t201;
t271 = -qJD(1) * t104 + qJD(2) * t40;
t270 = t298 * t234;
t103 = t293 + t360;
t127 = t199 * t228 + t201 * t226;
t268 = -qJD(1) * t103 + qJD(2) * t127;
t105 = t295 + t362;
t128 = -t199 * t226 + t201 * t228;
t267 = -qJD(1) * t105 + qJD(2) * t128;
t263 = t221 * t408 + t383 / 0.2e1;
t261 = t199 * t270;
t137 = (t231 / 0.2e1 - t230 / 0.2e1) * t199;
t260 = -qJD(2) * t137 + t311;
t259 = qJD(6) * t193 + t320;
t257 = t196 * t237 * t359 + qJD(3) * t137;
t141 = t219 * t196;
t256 = qJD(2) * t141 + t287;
t252 = t90 / 0.2e1 + t263;
t122 = t291 - t172 / 0.2e1;
t152 = t290 - t209 / 0.2e1;
t251 = qJD(1) * t122 + qJD(2) * t152 + qJD(3) * t222;
t15 = t252 * t234;
t42 = t302 * t237;
t249 = qJD(1) * t42 - qJD(2) * t15 - t222 * t329;
t17 = t252 * t237;
t41 = t302 * t234;
t248 = -qJD(1) * t41 - qJD(2) * t17 + t222 * t331;
t182 = t193 * qJD(3);
t163 = -t181 - t354;
t132 = t137 * qJD(6);
t121 = 0.2e1 * t290 - t380;
t115 = t227 + t247;
t109 = t253 * t376 + t293;
t108 = t201 * t307 + t360;
t107 = t294 + t362;
t106 = t295 + t361;
t93 = t109 * qJD(2);
t91 = t102 * qJD(2);
t63 = 0.2e1 * t291 - t381;
t48 = -t383 / 0.2e1 + t225 * t201 / 0.2e1 + t227 + t385 / 0.2e1 + t403 / 0.2e1;
t44 = 0.2e1 * t309;
t43 = 0.2e1 * t429 * t406;
t27 = t397 / 0.2e1 + t384 * t431 + t234 * t308 + t304;
t26 = -t398 / 0.2e1 + t199 * t309 + t237 * t308 + t305;
t18 = -t399 / 0.2e1 + t263 * t237 + t430 * t406;
t16 = -t400 / 0.2e1 + t263 * t234 + t430 * t405;
t9 = t387 / 0.2e1 + t242;
t7 = t386 / 0.2e1 + t241;
t6 = pkin(3) * t254 + t228 * t307 - t243;
t1 = t240 - t264;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t436, 0, 0, 0, -t436, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t318, -t317, 0, 0, 0, 0, 0, (-t236 * t357 - t239 * t330) * t233, (t235 * t358 - t239 * t328) * t233, qJD(3) * t108 + t297, qJD(3) * t107 + t296, t435, (t226 * t377 + t274) * qJD(2) + t6 * qJD(3) + t441, t435, qJD(3) * t109 - t297, qJD(3) * t106 - t296, (t134 * t377 + t274) * qJD(2) + t1 * qJD(3) + t109 * qJD(5) + t441, 0, 0, 0, 0, 0, ((-t234 * t377 + t374) * t201 - t160 * t138) * qJD(2) + t9 * qJD(3) + t27 * qJD(6), (-(t237 * t377 + t375) * t201 + t160 * t384) * qJD(2) + t7 * qJD(3) + t26 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t186 - t235 * t317, qJD(3) * t185 - t238 * t317, qJD(2) * t108 - t351, qJD(2) * t107 - t350, t422, t6 * qJD(2) + (t125 * t232 - t394 * t429) * t401, t422, t93 + t351, qJD(2) * t106 + t350, t1 * qJD(2) + (t125 * t222 + t225 * t429) * qJD(3) + t63 * qJD(5), 0, 0, 0, 0, 0, qJD(2) * t9 + qJD(6) * t44 + t125 * t331, qJD(2) * t7 + qJD(6) * t43 + t125 * t329; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t419, 0, 0, 0, t419, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t63 + t93, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t27 + qJD(3) * t44 + qJD(6) * t88, qJD(2) * t26 + qJD(3) * t43 - qJD(6) * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t103 * qJD(3), -t105 * qJD(3), t421, -qJD(3) * t5 + t440, t421, -t102 * qJD(3), -t104 * qJD(3), qJD(3) * t2 - qJD(5) * t102 + t440, 0, 0, 0, 0, 0, qJD(3) * t10 - qJD(6) * t24, qJD(3) * t8 - qJD(6) * t25; 0, 0, 0, 0, t235 * t328, t220 * qJD(3), 0, 0, 0, -pkin(2) * t330, -pkin(2) * t328, t127 * qJD(3), t128 * qJD(3), -t438, qJD(3) * t32 + t423, -t438, qJD(3) * t39 + t199 * t333, qJD(3) * t40 + t339, qJD(3) * t19 - t134 * t333 + t423, t196 * t312 + t230 * t319, -qJD(6) * t141 + t201 * t287, qJD(3) * t95 + t199 * t322, qJD(3) * t96 - t199 * t323, -t319, qJD(3) * t3 - qJD(4) * t97 + qJD(6) * t21 + t356 * t414, qJD(3) * t4 + qJD(4) * t64 + qJD(6) * t20 + t237 * t339; 0, 0, 0, 0, t313, t332, t328, -t330, 0, -pkin(8) * t328 - t326, pkin(8) * t330 - t325, t268 - t343, t267 - t342, (t199 * t394 - t201 * t232) * t401 + t416, (t155 * t232 - t394 * t428) * t401 + t115 * qJD(4) + t283, (-t199 * t225 - t382) * qJD(3) - t336 + t416, t343 - t353 + t369, t271 + t342, (t155 * t222 + t225 * t428) * qJD(3) + t48 * qJD(4) + t121 * qJD(5) + t284, t132 + (t230 * t338 + t311) * t201, -0.2e1 * t199 * t312 + t201 * t255, -t199 * t329 + t366, t315 + t365, -t259, (-t237 * t272 + t393) * qJD(3) - t138 * qJD(5) + t16 * qJD(6) + t285, t119 * t329 + t18 * qJD(6) + (qJD(3) * t272 + t336) * t234 + t286; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t439, qJD(3) * t115 + t415, -t439, 0, 0, qJD(3) * t48 + t415, 0, 0, 0, 0, 0, -t364, t367; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t188, t320, t340, qJD(3) * t121 - t321 - t353, 0, 0, 0, 0, 0, t314 - t346, t237 * t340 + t315; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t257, -t256, t298 * t138, -t261, -t182, qJD(3) * t16 - qJD(6) * t38 + t279, qJD(3) * t18 + qJD(6) * t37 + t280; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103 * qJD(2), t105 * qJD(2), -t422, t5 * qJD(2), -t422, t91, t104 * qJD(2), -qJD(2) * t2 + qJD(5) * t122, 0, 0, 0, 0, 0, -qJD(2) * t10 - qJD(6) * t42, -qJD(2) * t8 + qJD(6) * t41; 0, 0, 0, 0, -t313, -t332, 0, 0, 0, t326, t325, -t268 - t190, -t267 + t337, -t416, qJD(4) * t116 - t283, -t416, -t369 + t396, -t271 - t337, -qJD(4) * t47 + qJD(5) * t152 - t284, -t230 * t320 + t132, -0.2e1 * t237 * t261, -t323 - t366, -t322 - t365, t259, qJD(6) * t15 - t234 * t337 - t285, -qJD(4) * t138 + qJD(6) * t17 - t286; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), t222 * qJD(5), -t312, t219 * qJD(6), 0, 0, 0, t222 * t354 + t356, qJD(5) * t237 - t222 * t355; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t335, t338, 0, t352, 0, t335, -t338, -t368, 0, 0, 0, 0, 0, -t316, -t347; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t251, 0, 0, 0, 0, 0, t331, t329; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t260, -t255, -t270, t163, t341, -t221 * t355 - t249, -t221 * t354 - t248; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t420, 0, 0, 0, t420, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t334, -t188, t439, -qJD(3) * t116 - t415, t439, -t334, t188, qJD(3) * t47 - t333 - t415, 0, 0, 0, 0, 0, t315 - t322 + t364, t346 + t348 - t367; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t335, -t338, 0, -t352, 0, -t335, t338, t368, 0, 0, 0, 0, 0, t316, t347; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t335, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t163, t434; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t122 + t91, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t320, -t340, -qJD(3) * t152 + t321 + t396, 0, 0, 0, 0, 0, -t314 - t348, (-qJD(6) * t201 - t340) * t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t251, 0, 0, 0, 0, 0, -t331, -t329; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t335, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t434, -t298 * t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t24 + qJD(3) * t42, qJD(2) * t25 - qJD(3) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t257, t256, (-t237 * t338 + t331) * t201, (t316 + t329) * t201, -t182, -qJD(3) * t15 + qJD(5) * t135 + t190 * t237 - t279, -qJD(3) * t17 - qJD(4) * t135 + t237 * t333 - t280; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t260, t255, t234 * t335, t181, -t341, t249, t248; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t181, -t349; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t349, t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t11;
