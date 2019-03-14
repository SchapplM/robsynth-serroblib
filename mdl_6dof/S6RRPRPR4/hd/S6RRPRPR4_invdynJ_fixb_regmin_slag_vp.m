% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% tau_reg [6x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRPR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:26:22
% EndTime: 2019-03-09 10:26:39
% DurationCPUTime: 7.87s
% Computational Cost: add. (13192->569), mult. (37576->789), div. (0->0), fcn. (31290->16), ass. (0->311)
t264 = sin(pkin(11));
t265 = sin(pkin(6));
t272 = sin(qJ(2));
t370 = qJD(1) * t272;
t347 = t265 * t370;
t267 = cos(pkin(11));
t276 = cos(qJ(2));
t381 = t276 * t267;
t350 = t265 * t381;
t202 = qJD(1) * t350 - t264 * t347;
t196 = qJD(4) - t202;
t271 = sin(qJ(4));
t275 = cos(qJ(4));
t252 = pkin(2) * t264 + pkin(9);
t379 = qJ(5) + t252;
t328 = qJD(4) * t379;
t268 = cos(pkin(6));
t387 = t268 * t272;
t245 = pkin(1) * t387;
t390 = t265 * t276;
t412 = pkin(8) + qJ(3);
t188 = (t390 * t412 + t245) * qJD(1);
t180 = t264 * t188;
t421 = pkin(1) * t268;
t246 = t276 * t421;
t240 = qJD(1) * t246;
t344 = t412 * t272;
t321 = t265 * t344;
t187 = -qJD(1) * t321 + t240;
t121 = t187 * t267 - t180;
t306 = t264 * t276 + t267 * t272;
t296 = qJD(1) * t306;
t205 = t265 * t296;
t141 = pkin(2) * t347 + pkin(3) * t205 - pkin(9) * t202;
t376 = t275 * t121 + t271 * t141;
t396 = t202 * t271;
t442 = qJ(5) * t396 + t275 * qJD(5) - t271 * t328 - t376;
t125 = t275 * t141;
t441 = -pkin(4) * t205 - t125 + (qJ(5) * t202 - t328) * t275 + (-qJD(5) + t121) * t271;
t367 = qJD(4) * t271;
t440 = t367 - t396;
t371 = qJD(1) * t268;
t243 = qJD(2) + t371;
t170 = t205 * t271 - t275 * t243;
t263 = sin(pkin(12));
t266 = cos(pkin(12));
t307 = -t205 * t275 - t243 * t271;
t330 = -t266 * t170 + t263 * t307;
t429 = qJD(6) - t330;
t270 = sin(qJ(6));
t274 = cos(qJ(6));
t308 = -t170 * t263 - t266 * t307;
t82 = -t274 * t196 + t270 * t308;
t439 = t429 * t82;
t220 = t263 * t271 - t266 * t275;
t133 = t220 * t202;
t214 = t220 * qJD(4);
t438 = t133 - t214;
t222 = t263 * t275 + t266 * t271;
t374 = t196 * t222;
t256 = pkin(2) * t276 + pkin(1);
t313 = t256 * qJDD(1);
t369 = qJD(2) * t272;
t346 = t265 * t369;
t320 = qJD(1) * t346;
t358 = pkin(2) * t320 + qJDD(3);
t437 = t265 * t313 - t358;
t336 = t274 * t429;
t362 = qJD(1) * qJD(2);
t343 = t276 * t362;
t155 = -t264 * t320 + (qJDD(1) * t306 + t267 * t343) * t265;
t361 = qJDD(1) * t268;
t242 = qJDD(2) + t361;
t366 = qJD(4) * t275;
t88 = t275 * t155 - t205 * t367 + t271 * t242 + t243 * t366;
t89 = -qJD(4) * t307 + t271 * t155 - t275 * t242;
t51 = -t263 * t88 - t266 * t89;
t50 = qJDD(6) - t51;
t408 = t270 * t50;
t436 = -t336 * t429 - t408;
t209 = t306 * t268;
t221 = t264 * t272 - t381;
t273 = sin(qJ(1));
t277 = cos(qJ(1));
t162 = t209 * t277 - t221 * t273;
t260 = qJ(4) + pkin(12);
t257 = sin(t260);
t258 = cos(t260);
t389 = t265 * t277;
t136 = -t162 * t258 + t257 * t389;
t297 = t268 * t221;
t161 = -t273 * t306 - t277 * t297;
t435 = t136 * t270 - t161 * t274;
t434 = t136 * t274 + t161 * t270;
t405 = t263 * t442 - t441 * t266;
t404 = t441 * t263 + t266 * t442;
t186 = t268 * pkin(2) + t246 - t321;
t373 = pkin(8) * t390 + t245;
t197 = qJ(3) * t390 + t373;
t131 = t264 * t186 + t267 * t197;
t119 = pkin(9) * t268 + t131;
t392 = t265 * t272;
t207 = t264 * t392 - t350;
t208 = t306 * t265;
t318 = t256 * t265;
t143 = t207 * pkin(3) - t208 * pkin(9) - t318;
t377 = t275 * t119 + t271 * t143;
t183 = t208 * t271 - t268 * t275;
t184 = t208 * t275 + t268 * t271;
t112 = -t183 * t263 + t184 * t266;
t200 = t207 * t274;
t431 = -t112 * t270 + t200;
t388 = t267 * t188;
t120 = t187 * t264 + t388;
t430 = pkin(4) * t440 - t120;
t163 = t273 * t209 + t221 * t277;
t391 = t265 * t273;
t144 = t163 * t271 + t275 * t391;
t301 = t162 * t271 + t275 * t389;
t428 = -g(1) * t144 + g(2) * t301 + g(3) * t183;
t251 = pkin(4) * t263 + pkin(10);
t359 = qJDD(1) * t276;
t341 = t265 * t359;
t226 = t267 * t341;
t295 = t306 * qJD(2);
t360 = qJDD(1) * t272;
t342 = t264 * t360;
t153 = qJDD(4) - t226 + (qJD(1) * t295 + t342) * t265;
t353 = pkin(1) * t359;
t239 = t268 * t353;
t356 = qJD(2) * t421;
t323 = qJD(1) * t356;
t338 = qJD(2) * t412;
t368 = qJD(3) * t272;
t105 = -t272 * t323 + t242 * pkin(2) + t239 + (-qJDD(1) * t344 + (-t276 * t338 - t368) * qJD(1)) * t265;
t287 = qJD(3) * t276 - t272 * t338;
t348 = pkin(8) * t341 + qJDD(1) * t245 + t276 * t323;
t114 = (qJ(3) * t359 + qJD(1) * t287) * t265 + t348;
t72 = t264 * t105 + t267 * t114;
t66 = pkin(9) * t242 + t72;
t154 = t226 + (-qJD(2) * t296 - t342) * t265;
t80 = -t154 * pkin(3) - t155 * pkin(9) - t437;
t339 = -t271 * t66 + t275 * t80;
t173 = t243 * pkin(2) + t187;
t107 = t264 * t173 + t388;
t102 = pkin(9) * t243 + t107;
t212 = -qJD(1) * t318 + qJD(3);
t117 = -t202 * pkin(3) - t205 * pkin(9) + t212;
t74 = t102 * t275 + t117 * t271;
t282 = -qJD(4) * t74 + t339;
t16 = t153 * pkin(4) - t88 * qJ(5) + qJD(5) * t307 + t282;
t300 = t102 * t367 - t117 * t366 - t271 * t80 - t275 * t66;
t18 = -qJ(5) * t89 - qJD(5) * t170 - t300;
t5 = t16 * t266 - t18 * t263;
t3 = -pkin(5) * t153 - t5;
t427 = (-pkin(4) * t307 + pkin(5) * t308 - pkin(10) * t330 + qJD(6) * t251) * t429 + g(1) * (t163 * t257 + t258 * t391) + g(2) * (-t162 * t257 - t258 * t389) + g(3) * (-t208 * t257 + t258 * t268) + t3;
t365 = qJD(6) * t270;
t93 = -t133 * t274 + t205 * t270;
t290 = t214 * t274 + t222 * t365 + t93;
t394 = t222 * t274;
t426 = -t290 * t429 + t50 * t394;
t52 = -t263 * t89 + t266 * t88;
t337 = -t274 * t153 + t270 * t52;
t84 = t196 * t270 + t274 * t308;
t30 = qJD(6) * t84 + t337;
t425 = -t220 * t30 - t374 * t82;
t71 = t105 * t267 - t264 * t114;
t65 = -pkin(3) * t242 - t71;
t42 = pkin(4) * t89 + qJDD(5) + t65;
t14 = -pkin(5) * t51 - pkin(10) * t52 + t42;
t73 = -t102 * t271 + t275 * t117;
t57 = qJ(5) * t307 + t73;
t47 = pkin(4) * t196 + t57;
t58 = -qJ(5) * t170 + t74;
t55 = t266 * t58;
t22 = t263 * t47 + t55;
t20 = pkin(10) * t196 + t22;
t106 = t173 * t267 - t180;
t101 = -pkin(3) * t243 - t106;
t81 = pkin(4) * t170 + qJD(5) + t101;
t41 = -pkin(5) * t330 - pkin(10) * t308 + t81;
t312 = t20 * t270 - t274 * t41;
t6 = t263 * t16 + t266 * t18;
t4 = pkin(10) * t153 + t6;
t1 = -t312 * qJD(6) + t270 * t14 + t274 * t4;
t255 = pkin(4) * t275 + pkin(3);
t420 = pkin(2) * t267;
t305 = -t255 - t420;
t159 = pkin(5) * t220 - pkin(10) * t222 + t305;
t219 = t379 * t275;
t335 = t379 * t271;
t169 = t266 * t219 - t263 * t335;
t164 = t273 * t297 - t277 * t306;
t292 = g(1) * t164 + g(2) * t161 - g(3) * t207;
t424 = t258 * t292 + (-pkin(5) * t374 + pkin(10) * t438 + qJD(6) * t169 - t430) * t429 - t159 * t50;
t259 = t265 ^ 2;
t423 = 0.2e1 * t259;
t422 = pkin(1) * t259;
t382 = t273 * t276;
t384 = t272 * t277;
t217 = -t268 * t382 - t384;
t418 = g(1) * t217;
t417 = g(1) * t273;
t415 = g(3) * t276;
t414 = t82 * t308;
t413 = t84 * t308;
t204 = t221 * t265 * qJD(2);
t129 = -qJD(4) * t183 - t204 * t275;
t203 = t265 * t295;
t241 = t276 * t356;
t174 = t265 * t287 + t241;
t345 = t412 * t265;
t175 = -t265 * t368 + (-t276 * t345 - t245) * qJD(2);
t104 = t174 * t267 + t175 * t264;
t142 = pkin(2) * t346 + pkin(3) * t203 + pkin(9) * t204;
t333 = -t271 * t104 + t275 * t142;
t28 = t203 * pkin(4) - t129 * qJ(5) - qJD(4) * t377 - t184 * qJD(5) + t333;
t128 = qJD(4) * t184 - t204 * t271;
t298 = t275 * t104 - t119 * t367 + t271 * t142 + t143 * t366;
t36 = -qJ(5) * t128 - qJD(5) * t183 + t298;
t10 = t263 * t28 + t266 * t36;
t332 = -t119 * t271 + t275 * t143;
t61 = pkin(4) * t207 - qJ(5) * t184 + t332;
t69 = -qJ(5) * t183 + t377;
t38 = t263 * t61 + t266 * t69;
t409 = t263 * t58;
t364 = qJD(6) * t274;
t29 = t270 * t153 + t196 * t364 + t274 * t52 - t308 * t365;
t407 = t29 * t270;
t406 = pkin(5) * t205 + t405;
t400 = t170 * t196;
t399 = t170 * t205;
t398 = t307 * t196;
t397 = t307 * t205;
t395 = t207 * t270;
t278 = qJD(1) ^ 2;
t393 = t259 * t278;
t386 = t271 * t153;
t385 = t272 * t273;
t383 = t273 * t256;
t380 = t276 * t277;
t261 = t272 ^ 2;
t372 = -t276 ^ 2 + t261;
t363 = qJD(2) - t243;
t357 = t272 * t422;
t352 = pkin(8) * t360;
t351 = t276 * t393;
t349 = t268 * t380;
t340 = g(2) * t389 - g(3) * t268;
t329 = t162 * t275 - t271 * t389;
t103 = t174 * t264 - t267 * t175;
t130 = t186 * t267 - t264 * t197;
t327 = t196 * t275;
t326 = t243 + t371;
t325 = qJD(1) * t363;
t324 = t242 + t361;
t322 = t29 * t220 + t374 * t84;
t211 = pkin(2) * t387 - t345;
t319 = pkin(4) * t265 * t271 - t211;
t317 = g(1) * t277 + g(2) * t273;
t12 = t20 * t274 + t270 * t41;
t9 = -t263 * t36 + t266 * t28;
t21 = t266 * t47 - t409;
t37 = -t263 * t69 + t266 * t61;
t33 = pkin(10) * t207 + t38;
t111 = t266 * t183 + t184 * t263;
t118 = -pkin(3) * t268 - t130;
t284 = pkin(4) * t183 + t118;
t45 = pkin(5) * t111 - pkin(10) * t112 + t284;
t311 = t270 * t45 + t274 * t33;
t310 = -t270 * t33 + t274 * t45;
t91 = t112 * t274 + t395;
t304 = t274 * t50 + (t270 * t330 - t365) * t429;
t303 = t275 * t153 - t196 * t440;
t302 = pkin(4) * t128 + t103;
t294 = t101 * t196 - t252 * t153;
t293 = g(1) * t163 - g(2) * t162 - g(3) * t208;
t92 = -t133 * t270 - t274 * t205;
t291 = -t214 * t270 + t222 * t364 - t92;
t285 = t373 * t243;
t19 = -pkin(5) * t196 - t21;
t26 = t266 * t57 - t409;
t283 = -t251 * t50 + (t19 + t26) * t429;
t2 = -qJD(6) * t12 + t274 * t14 - t270 * t4;
t281 = qJD(4) * t196 * t252 + t292 + t65;
t280 = -t222 * t408 - t291 * t429;
t279 = -t169 * t50 + t3 * t222 + (pkin(10) * t205 - qJD(6) * t159 - t404) * t429 + t293;
t269 = -qJ(5) - pkin(9);
t254 = -pkin(3) - t420;
t253 = -pkin(4) * t266 - pkin(5);
t235 = t277 * t256;
t233 = pkin(2) * t349;
t218 = -t268 * t385 + t380;
t216 = -t268 * t384 - t382;
t215 = -t349 + t385;
t177 = t208 * t258 + t257 * t268;
t168 = t219 * t263 + t266 * t335;
t145 = -t163 * t275 + t271 * t391;
t138 = -t163 * t258 + t257 * t391;
t86 = t138 * t274 - t164 * t270;
t85 = -t138 * t270 - t164 * t274;
t79 = -t128 * t263 + t129 * t266;
t78 = t266 * t128 + t129 * t263;
t44 = qJD(6) * t91 - t203 * t274 + t270 * t79;
t43 = qJD(6) * t431 + t203 * t270 + t274 * t79;
t32 = -pkin(5) * t207 - t37;
t25 = t263 * t57 + t55;
t24 = pkin(5) * t78 - pkin(10) * t79 + t302;
t8 = pkin(10) * t203 + t10;
t7 = -pkin(5) * t203 - t9;
t11 = [qJDD(1), -g(2) * t277 + t417, t317 (qJDD(1) * t261 + 0.2e1 * t272 * t343) * t259 (t272 * t359 - t362 * t372) * t423 (qJD(2) * t276 * t326 + t272 * t324) * t265 (t276 * t324 - t326 * t369) * t265, t242 * t268, t353 * t423 + (-pkin(8) * t392 + t246) * t242 + (-t265 * t352 + t239) * t268 - g(1) * t216 - g(2) * t218 + (-t285 + (-t268 * t373 - 0.2e1 * t357) * qJD(1)) * qJD(2) -(-pkin(8) * t346 + t241) * t243 - t373 * t242 - (-pkin(8) * t320 + t348) * t268 - g(1) * t215 - g(2) * t217 + 0.2e1 * (-t343 - t360) * t422, t103 * t205 + t104 * t202 + t106 * t204 - t107 * t203 - t130 * t155 + t131 * t154 - t72 * t207 - t71 * t208 - t265 * t317, t72 * t131 + t107 * t104 + t71 * t130 - t106 * t103 - g(1) * (-t211 * t277 - t383) - g(2) * (-t211 * t273 + t235) + (pkin(2) * t212 * t369 + t256 * t437) * t265, -t129 * t307 + t184 * t88, t128 * t307 - t129 * t170 - t183 * t88 - t184 * t89, t129 * t196 + t153 * t184 - t203 * t307 + t207 * t88, -t128 * t196 - t153 * t183 - t170 * t203 - t207 * t89, t153 * t207 + t196 * t203, t333 * t196 + t332 * t153 + t339 * t207 + t73 * t203 + t103 * t170 + t118 * t89 + t65 * t183 + t101 * t128 + g(1) * t329 - g(2) * t145 + (-t196 * t377 - t207 * t74) * qJD(4), -g(1) * t301 - g(2) * t144 + t101 * t129 - t103 * t307 + t118 * t88 - t153 * t377 + t65 * t184 - t196 * t298 - t74 * t203 + t207 * t300, -g(1) * t161 + g(2) * t164 + t10 * t330 - t111 * t6 - t112 * t5 - t21 * t79 - t22 * t78 - t308 * t9 - t37 * t52 + t38 * t51, t6 * t38 + t22 * t10 + t5 * t37 + t21 * t9 + t42 * t284 + t81 * t302 - g(1) * (-t161 * t269 - t162 * t255 + t277 * t319 - t383) - g(2) * (-t163 * t255 + t164 * t269 + t273 * t319 + t235) t29 * t91 + t43 * t84, t29 * t431 - t30 * t91 - t43 * t82 - t44 * t84, t111 * t29 + t429 * t43 + t50 * t91 + t78 * t84, -t111 * t30 - t429 * t44 + t431 * t50 - t78 * t82, t111 * t50 + t429 * t78 (-qJD(6) * t311 + t274 * t24 - t270 * t8) * t429 + t310 * t50 + t2 * t111 - t312 * t78 + t7 * t82 + t32 * t30 - t3 * t431 + t19 * t44 - g(1) * t434 - g(2) * t86 -(qJD(6) * t310 + t270 * t24 + t274 * t8) * t429 - t311 * t50 - t1 * t111 - t12 * t78 + t7 * t84 + t32 * t29 + t3 * t91 + t19 * t43 + g(1) * t435 - g(2) * t85; 0, 0, 0, -t272 * t351, t372 * t393 (t276 * t325 + t360) * t265 (-t363 * t370 + t359) * t265, t242, t278 * t357 - t418 + g(2) * t215 + t239 + (-t352 - t415) * t265 + (-qJD(2) * t373 + t285) * qJD(1), pkin(1) * t351 + g(1) * t218 - g(2) * t216 + t240 * t243 + (pkin(8) * t325 + g(3)) * t392 - t348 (t107 - t120) * t205 + (t106 - t121) * t202 + (t154 * t264 - t155 * t267) * pkin(2), -g(2) * t233 + t106 * t120 - t107 * t121 + (t72 * t264 + t71 * t267 - t418 + g(2) * t385 + (-t212 * t370 - t415) * t265) * pkin(2), t88 * t271 - t307 * t327 (t88 - t400) * t275 + (-t89 + t398) * t271, t196 * t327 + t386 + t397, t303 + t399, -t196 * t205, -t120 * t170 - t125 * t196 - t73 * t205 + t254 * t89 + (t121 * t196 + t294) * t271 - t281 * t275, t120 * t307 + t196 * t376 + t74 * t205 + t254 * t88 + t271 * t281 + t275 * t294, t168 * t52 + t169 * t51 - t21 * t438 - t22 * t374 - t220 * t6 - t222 * t5 + t308 * t405 + t330 * t404 + t293, t6 * t169 - t5 * t168 + t42 * t305 - g(1) * (pkin(2) * t217 + t163 * t269 + t164 * t255) - g(2) * (-pkin(2) * t385 + t161 * t255 - t162 * t269 + t233) - g(3) * (pkin(2) * t390 - t207 * t255 - t208 * t269) + t430 * t81 + t404 * t22 - t405 * t21, t29 * t394 - t290 * t84, t93 * t82 + t84 * t92 - (-t270 * t84 - t274 * t82) * t214 + (-t407 - t274 * t30 + (t270 * t82 - t274 * t84) * qJD(6)) * t222, t322 + t426, t280 + t425, t50 * t220 + t374 * t429, t168 * t30 + t291 * t19 + t2 * t220 + t279 * t270 - t274 * t424 - t312 * t374 + t406 * t82, -t1 * t220 - t374 * t12 + t168 * t29 - t290 * t19 + t270 * t424 + t279 * t274 + t406 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t202 ^ 2 - t205 ^ 2, t106 * t205 - t107 * t202 + (-t313 - t417) * t265 + t340 + t358, 0, 0, 0, 0, 0, t303 - t399, -t196 ^ 2 * t275 - t386 + t397, t220 * t52 + t222 * t51 + t308 * t374 + t330 * t438, -g(1) * t391 - t81 * t205 - t21 * t374 + t22 * t438 - t5 * t220 + t6 * t222 + t340, 0, 0, 0, 0, 0, t280 - t425, t322 - t426; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t307 * t170, -t170 ^ 2 + t307 ^ 2, t88 + t400, -t89 - t398, t153, t101 * t307 + t74 * t196 + t282 + t428, g(1) * t145 + g(2) * t329 + g(3) * t184 + t101 * t170 + t73 * t196 + t300 (t263 * t51 - t266 * t52) * pkin(4) + (t21 - t26) * t330 + (t22 - t25) * t308, t21 * t25 - t22 * t26 + (t6 * t263 + t5 * t266 + t307 * t81 + t428) * pkin(4), t336 * t84 + t407 (t29 - t439) * t274 + (-t429 * t84 - t30) * t270, -t413 - t436, t304 + t414, -t429 * t308, -t25 * t82 + t253 * t30 + t283 * t270 - t274 * t427 + t308 * t312, t12 * t308 - t25 * t84 + t253 * t29 + t270 * t427 + t283 * t274; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t308 ^ 2 - t330 ^ 2, t21 * t308 - t22 * t330 + t292 + t42, 0, 0, 0, 0, 0, t304 - t414, -t413 + t436; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84 * t82, -t82 ^ 2 + t84 ^ 2, t29 + t439, -t337 + (-qJD(6) + t429) * t84, t50, t12 * t429 - t19 * t84 - g(1) * t85 - g(2) * t435 - g(3) * (-t177 * t270 + t200) + t2, -t312 * t429 + t19 * t82 + g(1) * t86 - g(2) * t434 - g(3) * (-t177 * t274 - t395) - t1;];
tau_reg  = t11;