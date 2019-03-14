% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRRPR5
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% tau_reg [6x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRRPR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:28:22
% EndTime: 2019-03-08 23:28:46
% DurationCPUTime: 10.71s
% Computational Cost: add. (9519->590), mult. (24355->867), div. (0->0), fcn. (20926->18), ass. (0->303)
t263 = sin(pkin(7));
t272 = sin(qJ(3));
t276 = cos(qJ(3));
t305 = t263 * (pkin(3) * t272 - pkin(10) * t276);
t264 = sin(pkin(6));
t273 = sin(qJ(2));
t401 = t264 * t273;
t355 = qJD(1) * t401;
t459 = qJD(3) * t305 - t263 * t355;
t267 = cos(pkin(7));
t277 = cos(qJ(2));
t390 = t276 * t277;
t395 = t272 * t273;
t302 = -t267 * t395 + t390;
t178 = t302 * t264;
t404 = t263 * t272;
t248 = pkin(9) * t404;
t396 = t267 * t276;
t442 = pkin(2) * t396 - t248;
t458 = -qJD(1) * t178 + t442 * qJD(3);
t271 = sin(qJ(4));
t275 = cos(qJ(4));
t397 = t267 * t272;
t403 = t263 * t276;
t384 = pkin(2) * t397 + pkin(9) * t403;
t187 = pkin(10) * t267 + t384;
t317 = -pkin(3) * t276 - pkin(10) * t272;
t188 = (-pkin(2) + t317) * t263;
t443 = t275 * t187 + t271 * t188;
t457 = qJD(4) * t443 + t458 * t271 - t275 * t459;
t374 = qJD(4) * t275;
t375 = qJD(4) * t271;
t456 = -t187 * t375 + t188 * t374 + t271 * t459 + t458 * t275;
t378 = qJD(2) * t276;
t353 = t263 * t378;
t455 = qJD(4) - t353;
t209 = -t275 * t267 + t271 * t404;
t376 = qJD(3) * t276;
t350 = t263 * t376;
t324 = t275 * t350;
t148 = -qJD(4) * t209 + t324;
t359 = t275 * t404;
t210 = t267 * t271 + t359;
t377 = qJD(3) * t272;
t351 = t263 * t377;
t454 = pkin(4) * t351 - qJ(5) * t148 - qJD(5) * t210 - t457;
t348 = t271 * t376;
t149 = qJD(4) * t210 + t263 * t348;
t453 = qJ(5) * t149 + qJD(5) * t209 - t456;
t325 = t271 * t353;
t269 = -qJ(5) - pkin(10);
t342 = qJD(4) * t269;
t380 = qJD(2) * t263;
t212 = pkin(9) * t380 + t355;
t268 = cos(pkin(6));
t382 = qJD(1) * t268;
t356 = t263 * t382;
t381 = qJD(1) * t277;
t224 = qJD(2) * pkin(2) + t264 * t381;
t409 = t224 * t267;
t120 = -t272 * t212 + t276 * (t356 + t409);
t194 = qJD(2) * t305;
t388 = t275 * t120 + t271 * t194;
t452 = qJ(5) * t325 + qJD(5) * t275 + t271 * t342 - t388;
t181 = t275 * t194;
t391 = t275 * t276;
t451 = t275 * t342 - t181 - (pkin(4) * t272 - qJ(5) * t391) * t380 + (-qJD(5) + t120) * t271;
t379 = qJD(2) * t267;
t247 = qJD(3) + t379;
t354 = t272 * t380;
t326 = t271 * t354;
t174 = -t275 * t247 + t326;
t176 = t247 * t271 + t275 * t354;
t262 = sin(pkin(13));
t265 = cos(pkin(13));
t336 = -t265 * t174 - t176 * t262;
t438 = qJD(6) - t336;
t270 = sin(qJ(6));
t274 = cos(qJ(6));
t309 = -t174 * t262 + t265 * t176;
t89 = t270 * t309 - t274 * t455;
t450 = t438 * t89;
t213 = t262 * t271 - t265 * t275;
t151 = t213 * t353;
t208 = t213 * qJD(4);
t449 = t151 - t208;
t448 = pkin(9) * qJDD(2) * t263 + (qJD(2) * t381 + qJDD(1) * t273) * t264 + qJD(3) * t356;
t393 = t273 * t276;
t394 = t272 * t277;
t304 = t267 * t393 + t394;
t177 = t304 * t264;
t349 = t267 * t377;
t386 = pkin(2) * t349 + pkin(9) * t350 - qJD(1) * t177;
t214 = t262 * t275 + t265 * t271;
t385 = t455 * t214;
t333 = t438 * t274;
t367 = qJDD(2) * t267;
t246 = qJDD(3) + t367;
t102 = qJD(2) * t324 - qJD(4) * t326 + qJDD(2) * t359 + t271 * t246 + t247 * t374;
t366 = qJDD(2) * t272;
t103 = -t275 * t246 + t247 * t375 + t263 * (qJD(2) * (t272 * t374 + t348) + t271 * t366);
t62 = -t102 * t262 - t265 * t103;
t61 = qJDD(6) - t62;
t447 = -t270 * t61 - t438 * t333;
t427 = t453 * t262 + t265 * t454;
t426 = t262 * t454 - t453 * t265;
t418 = t262 * t452 - t265 * t451;
t417 = t262 * t451 + t265 * t452;
t439 = t267 * t390 - t395;
t144 = -t264 * t439 - t268 * t403;
t141 = t144 * t274;
t303 = t267 * t394 + t393;
t145 = t264 * t303 + t268 * t404;
t400 = t264 * t277;
t202 = -t263 * t400 + t267 * t268;
t110 = -t145 * t271 + t202 * t275;
t111 = t145 * t275 + t202 * t271;
t66 = t110 * t262 + t111 * t265;
t445 = -t270 * t66 + t141;
t444 = pkin(4) * t149 + t386;
t416 = sin(pkin(12));
t340 = t416 * t273;
t266 = cos(pkin(12));
t398 = t266 * t277;
t203 = t268 * t398 - t340;
t339 = t416 * t277;
t399 = t266 * t273;
t204 = t268 * t399 + t339;
t402 = t264 * t266;
t361 = t263 * t402;
t107 = t203 * t397 + t204 * t276 - t272 * t361;
t205 = -t268 * t339 - t399;
t206 = -t268 * t340 + t398;
t341 = t264 * t416;
t319 = t263 * t341;
t109 = t206 * t276 + (t205 * t267 + t319) * t272;
t146 = -t203 * t263 - t267 * t402;
t147 = -t205 * t263 + t267 * t341;
t441 = -g(1) * (-t109 * t271 + t147 * t275) - g(2) * (-t107 * t271 + t146 * t275) - g(3) * t110;
t316 = g(1) * t206 + g(2) * t204;
t287 = -g(3) * t401 - t316;
t121 = t276 * t212 + t224 * t397 + t272 * t356;
t440 = -t121 + (-t325 + t375) * pkin(4);
t437 = pkin(4) * t103 + qJDD(5);
t253 = pkin(4) * t262 + pkin(11);
t259 = qJ(4) + pkin(13);
t256 = sin(t259);
t257 = cos(t259);
t365 = qJDD(2) * t276;
t244 = t263 * t365;
t369 = qJD(2) * qJD(3);
t346 = t272 * t369;
t190 = t263 * t346 + qJDD(4) - t244;
t245 = qJDD(1) * t400;
t352 = qJD(2) * t401;
t323 = qJD(1) * t352;
t184 = qJDD(2) * pkin(2) + t245 - t323;
t368 = qJDD(1) * t268;
t344 = t263 * t368;
t285 = t184 * t397 - t212 * t377 + t272 * t344 + t276 * t448 + t376 * t409;
t59 = pkin(10) * t246 + t285;
t105 = pkin(10) * t247 + t121;
t243 = t267 * t382;
t137 = t243 + (qJD(2) * t317 - t224) * t263;
t70 = t105 * t275 + t137 * t271;
t242 = t267 * t368;
t298 = t346 - t365;
t345 = t276 * t369;
t299 = t345 + t366;
t88 = t242 + (pkin(3) * t298 - pkin(10) * t299 - t184) * t263;
t84 = t275 * t88;
t282 = -qJD(4) * t70 - t271 * t59 + t84;
t12 = pkin(4) * t190 - qJ(5) * t102 - qJD(5) * t176 + t282;
t300 = t105 * t375 - t137 * t374 - t271 * t88 - t275 * t59;
t14 = -qJ(5) * t103 - qJD(5) * t174 - t300;
t5 = t12 * t265 - t14 * t262;
t3 = -pkin(5) * t190 - t5;
t436 = t438 * (pkin(4) * t176 + pkin(5) * t309 - pkin(11) * t336 + qJD(6) * t253) + g(1) * (-t109 * t256 + t147 * t257) + g(2) * (-t107 * t256 + t146 * t257) + g(3) * (-t145 * t256 + t202 * t257) + t3;
t306 = -t184 * t396 + t212 * t376 + t224 * t349 + t272 * t448 - t276 * t344;
t433 = pkin(3) * t246;
t60 = t306 - t433;
t31 = t60 + t437;
t63 = t102 * t265 - t103 * t262;
t10 = -pkin(5) * t62 - pkin(11) * t63 + t31;
t69 = -t105 * t271 + t275 * t137;
t56 = -qJ(5) * t176 + t69;
t45 = pkin(4) * t455 + t56;
t57 = -qJ(5) * t174 + t70;
t51 = t265 * t57;
t22 = t262 * t45 + t51;
t20 = pkin(11) * t455 + t22;
t104 = -pkin(3) * t247 - t120;
t81 = pkin(4) * t174 + qJD(5) + t104;
t32 = -pkin(5) * t336 - pkin(11) * t309 + t81;
t311 = t20 * t270 - t274 * t32;
t6 = t262 * t12 + t265 * t14;
t4 = pkin(11) * t190 + t6;
t1 = -t311 * qJD(6) + t270 * t10 + t274 * t4;
t428 = -pkin(5) * t351 - t427;
t335 = -t187 * t271 + t275 * t188;
t87 = -pkin(4) * t403 - qJ(5) * t210 + t335;
t95 = -qJ(5) * t209 + t443;
t42 = t262 * t87 + t265 * t95;
t425 = t309 * t89;
t91 = t270 * t455 + t274 * t309;
t424 = t309 * t91;
t423 = t262 * t57;
t372 = qJD(6) * t274;
t373 = qJD(6) * t270;
t27 = t270 * t190 + t274 * t63 - t309 * t373 + t372 * t455;
t422 = t27 * t270;
t419 = pkin(5) * t354 + t418;
t415 = t144 * t270;
t414 = t174 * t455;
t413 = t176 * t455;
t412 = t202 * t263;
t411 = t214 * t270;
t410 = t214 * t274;
t408 = t256 * t263;
t407 = t257 * t270;
t406 = t257 * t274;
t258 = t263 ^ 2;
t278 = qJD(2) ^ 2;
t405 = t258 * t278;
t392 = t273 * t278;
t389 = qJDD(1) - g(3);
t260 = t272 ^ 2;
t383 = -t276 ^ 2 + t260;
t371 = qJD(3) - t247;
t360 = t270 * t403;
t255 = pkin(4) * t275 + pkin(3);
t347 = t271 * t269;
t338 = -t274 * t190 + t270 * t63;
t332 = t247 + t379;
t331 = t246 + t367;
t330 = t258 * t264 * t392;
t327 = t263 * t352;
t138 = t265 * t209 + t210 * t262;
t139 = -t209 * t262 + t210 * t265;
t186 = t248 + (-pkin(2) * t276 - pkin(3)) * t267;
t286 = pkin(4) * t209 + t186;
t64 = pkin(5) * t138 - pkin(11) * t139 + t286;
t320 = -pkin(11) * t351 - qJD(6) * t64 - t426;
t38 = -pkin(11) * t403 + t42;
t93 = t148 * t262 + t265 * t149;
t94 = t148 * t265 - t149 * t262;
t318 = -pkin(5) * t93 + pkin(11) * t94 + qJD(6) * t38 - t444;
t140 = pkin(5) * t213 - pkin(11) * t214 - t255;
t315 = pkin(11) * t354 - qJD(6) * t140 - t417;
t237 = t269 * t275;
t157 = -t265 * t237 + t262 * t347;
t314 = -pkin(5) * t385 + pkin(11) * t449 + qJD(6) * t157 - t440;
t8 = t20 * t274 + t270 * t32;
t21 = t265 * t45 - t423;
t41 = -t262 * t95 + t265 * t87;
t310 = t274 * t66 + t415;
t307 = t274 * t61 + (t270 * t336 - t373) * t438;
t112 = t139 * t270 + t274 * t403;
t106 = -t203 * t396 + t204 * t272 + t276 * t361;
t108 = -t205 * t396 + t206 * t272 - t276 * t319;
t296 = g(1) * t108 + g(2) * t106 + g(3) * t144;
t295 = -g(1) * t109 - g(2) * t107 - g(3) * t145;
t125 = t203 * t272 + t204 * t396;
t127 = t205 * t272 + t206 * t396;
t294 = g(1) * t127 + g(2) * t125 + g(3) * t177;
t126 = t203 * t276 - t204 * t397;
t128 = t205 * t276 - t206 * t397;
t293 = g(1) * t128 + g(2) * t126 + g(3) * t178;
t133 = -t151 * t270 - t274 * t354;
t291 = -t208 * t270 + t214 * t372 - t133;
t135 = -t151 * t274 + t270 * t354;
t290 = -t208 * t274 - t214 * t373 - t135;
t19 = -pkin(5) * t455 - t21;
t26 = t265 * t56 - t423;
t284 = -t253 * t61 + (t19 + t26) * t438;
t283 = -pkin(10) * t190 + t104 * t455;
t2 = -qJD(6) * t8 + t274 * t10 - t270 * t4;
t281 = -pkin(10) * qJD(4) * t455 + t296 - t60;
t280 = t296 - t306;
t254 = -pkin(4) * t265 - pkin(5);
t164 = -t224 * t263 + t243;
t156 = -t237 * t262 - t265 * t347;
t142 = -t184 * t263 + t242;
t131 = t178 * t257 + t401 * t408;
t113 = t139 * t274 - t360;
t100 = t268 * t350 + (t302 * qJD(2) + qJD(3) * t439) * t264;
t99 = t268 * t351 + (qJD(2) * t304 + qJD(3) * t303) * t264;
t98 = t145 * t257 + t202 * t256;
t86 = t128 * t257 + t206 * t408;
t85 = t126 * t257 + t204 * t408;
t74 = t109 * t257 + t147 * t256;
t72 = t107 * t257 + t146 * t256;
t65 = -t265 * t110 + t111 * t262;
t50 = qJD(4) * t110 + t100 * t275 + t271 * t327;
t49 = -qJD(4) * t111 - t100 * t271 + t275 * t327;
t48 = -qJD(6) * t360 + t139 * t372 + t270 * t94 - t274 * t351;
t47 = -qJD(6) * t112 + t270 * t351 + t274 * t94;
t37 = pkin(5) * t403 - t41;
t28 = qJD(6) * t91 + t338;
t25 = t262 * t56 + t51;
t24 = t262 * t49 + t265 * t50;
t23 = t262 * t50 - t265 * t49;
t7 = [t389, 0 (qJDD(2) * t277 - t392) * t264 (-qJDD(2) * t273 - t277 * t278) * t264, 0, 0, 0, 0, 0, -t144 * t246 - t247 * t99 - t276 * t330 + t298 * t412, -t100 * t247 - t145 * t246 + t272 * t330 + t299 * t412, 0, 0, 0, 0, 0, t103 * t144 + t110 * t190 + t174 * t99 + t455 * t49, t102 * t144 - t111 * t190 + t176 * t99 - t455 * t50, t23 * t309 + t24 * t336 + t62 * t66 + t63 * t65, t144 * t31 - t21 * t23 + t22 * t24 - t5 * t65 + t6 * t66 + t81 * t99 - g(3), 0, 0, 0, 0, 0 (-qJD(6) * t310 - t24 * t270 + t274 * t99) * t438 + t445 * t61 + t23 * t89 + t65 * t28 -(qJD(6) * t445 + t24 * t274 + t270 * t99) * t438 - t310 * t61 + t23 * t91 + t65 * t27; 0, qJDD(2), -g(1) * t205 - g(2) * t203 - g(3) * t400 + t245, -t389 * t401 + t316 (qJDD(2) * t260 + 0.2e1 * t272 * t345) * t258, 0.2e1 * (t272 * t365 - t369 * t383) * t258 (t272 * t331 + t332 * t376) * t263 (t276 * t331 - t332 * t377) * t263, t246 * t267, t442 * t246 - t306 * t267 + (-t142 * t276 + t164 * t377) * t263 - t386 * t247 + (-pkin(2) * t298 + t276 * t323) * t258 - t293, -t384 * t246 - t285 * t267 + (t142 * t272 + t164 * t376) * t263 - t458 * t247 + (-pkin(2) * t299 - t272 * t323) * t258 + t294, t102 * t210 + t148 * t176, -t102 * t209 - t103 * t210 - t148 * t174 - t149 * t176, t148 * t455 + t190 * t210 + (-t102 * t276 + t176 * t377) * t263, -t149 * t455 - t190 * t209 + (t103 * t276 - t174 * t377) * t263 (-t190 * t276 + t377 * t455) * t263, t335 * t190 + t186 * t103 + t60 * t209 + t104 * t149 - t293 * t275 + (-(-t105 * t374 + t84) * t276 + t69 * t377 + (-(-qJD(4) * t137 - t59) * t276 + t287) * t271) * t263 - t457 * t455 + t386 * t174, t186 * t102 + t60 * t210 + t104 * t148 - t443 * t190 + t293 * t271 + (t275 * t287 - t276 * t300 - t377 * t70) * t263 - t456 * t455 + t386 * t176, -t138 * t6 - t139 * t5 - t21 * t94 - t22 * t93 - t309 * t427 + t336 * t426 - t41 * t63 + t42 * t62 - t294, t6 * t42 + t5 * t41 + t31 * t286 - g(1) * (pkin(2) * t205 - t127 * t269 + t128 * t255) - g(2) * (pkin(2) * t203 - t125 * t269 + t126 * t255) - g(3) * (pkin(2) * t400 - t177 * t269 + t178 * t255) + t444 * t81 + t426 * t22 + t427 * t21 + t287 * t263 * (pkin(4) * t271 + pkin(9)) t113 * t27 + t47 * t91, -t112 * t27 - t113 * t28 - t47 * t89 - t48 * t91, t113 * t61 + t138 * t27 + t438 * t47 + t91 * t93, -t112 * t61 - t138 * t28 - t438 * t48 - t89 * t93, t138 * t61 + t438 * t93 (-t270 * t38 + t274 * t64) * t61 + t2 * t138 - t311 * t93 + t37 * t28 + t3 * t112 + t19 * t48 - g(1) * (t127 * t270 + t274 * t86) - g(2) * (t125 * t270 + t274 * t85) - g(3) * (t131 * t274 + t177 * t270) + t428 * t89 + (t270 * t320 - t274 * t318) * t438 -(t270 * t64 + t274 * t38) * t61 - t1 * t138 - t8 * t93 + t37 * t27 + t3 * t113 + t19 * t47 - g(1) * (t127 * t274 - t270 * t86) - g(2) * (t125 * t274 - t270 * t85) - g(3) * (-t131 * t270 + t177 * t274) + t428 * t91 + (t270 * t318 + t274 * t320) * t438; 0, 0, 0, 0, -t272 * t276 * t405, t383 * t405 (t371 * t378 + t366) * t263, -t354 * t371 + t244, t246, t121 * t247 - t164 * t354 + t280, t120 * t247 - t164 * t353 - t285 - t295, t102 * t271 + t275 * t413 (t102 - t414) * t275 + (-t103 - t413) * t271, t455 * t374 + t271 * t190 + (-t176 * t272 - t391 * t455) * t380, -t455 * t375 + t275 * t190 + (t271 * t276 * t455 + t174 * t272) * t380, -t455 * t354, -t69 * t354 - pkin(3) * t103 - t121 * t174 - t181 * t455 + (t120 * t455 + t283) * t271 + t281 * t275, -pkin(3) * t102 - t121 * t176 - t271 * t281 + t275 * t283 + t354 * t70 + t388 * t455, t156 * t63 + t157 * t62 - t21 * t449 - t213 * t6 - t214 * t5 - t385 * t22 + t418 * t309 + t417 * t336 + t295, t6 * t157 - t5 * t156 - t31 * t255 - g(1) * (-t108 * t255 - t109 * t269) - g(2) * (-t106 * t255 - t107 * t269) - g(3) * (-t144 * t255 - t145 * t269) + t440 * t81 + t417 * t22 - t418 * t21, t27 * t410 + t290 * t91, t133 * t91 + t135 * t89 - (-t270 * t91 - t274 * t89) * t208 + (-t422 - t274 * t28 + (t270 * t89 - t274 * t91) * qJD(6)) * t214, t213 * t27 + t290 * t438 + t385 * t91 + t410 * t61, -t213 * t28 - t291 * t438 - t385 * t89 - t411 * t61, t213 * t61 + t385 * t438 (t140 * t274 - t157 * t270) * t61 + t2 * t213 + t156 * t28 + t3 * t411 - g(1) * (-t108 * t406 + t109 * t270) - g(2) * (-t106 * t406 + t107 * t270) - g(3) * (-t144 * t406 + t145 * t270) + t419 * t89 - t385 * t311 + (t270 * t315 - t274 * t314) * t438 + t291 * t19 -(t140 * t270 + t157 * t274) * t61 - t1 * t213 + t156 * t27 + t3 * t410 - g(1) * (t108 * t407 + t109 * t274) - g(2) * (t106 * t407 + t107 * t274) - g(3) * (t144 * t407 + t145 * t274) + t419 * t91 - t385 * t8 + (t270 * t314 + t274 * t315) * t438 + t290 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176 * t174, -t174 ^ 2 + t176 ^ 2, t102 + t414, -t103 + t413, t190, -t104 * t176 + t455 * t70 + t282 + t441, t104 * t174 + t69 * t455 - g(1) * (-t109 * t275 - t147 * t271) - g(2) * (-t107 * t275 - t146 * t271) + g(3) * t111 + t300 (t262 * t62 - t265 * t63) * pkin(4) + (-t26 + t21) * t336 + (t22 - t25) * t309, t21 * t25 - t22 * t26 + (-t81 * t176 + t6 * t262 + t5 * t265 + t441) * pkin(4), t333 * t91 + t422 (t27 - t450) * t274 + (-t438 * t91 - t28) * t270, -t424 - t447, t307 + t425, -t438 * t309, -t25 * t89 + t254 * t28 + t284 * t270 - t274 * t436 + t309 * t311, -t25 * t91 + t254 * t27 + t270 * t436 + t284 * t274 + t8 * t309; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t309 ^ 2 - t336 ^ 2, t21 * t309 - t22 * t336 - t280 - t433 + t437, 0, 0, 0, 0, 0, t307 - t425, -t424 + t447; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91 * t89, -t89 ^ 2 + t91 ^ 2, t27 + t450, -t338 + (-qJD(6) + t438) * t91, t61, t8 * t438 - t19 * t91 - g(1) * (t108 * t274 - t270 * t74) - g(2) * (t106 * t274 - t270 * t72) - g(3) * (-t270 * t98 + t141) + t2, -t311 * t438 + t19 * t89 - g(1) * (-t108 * t270 - t274 * t74) - g(2) * (-t106 * t270 - t274 * t72) - g(3) * (-t274 * t98 - t415) - t1;];
tau_reg  = t7;