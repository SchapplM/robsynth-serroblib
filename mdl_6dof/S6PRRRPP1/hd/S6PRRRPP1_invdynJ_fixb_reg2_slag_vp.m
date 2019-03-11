% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6PRRRPP1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRRPP1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP1_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPP1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_invdynJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:47:01
% EndTime: 2019-03-08 22:47:17
% DurationCPUTime: 9.47s
% Computational Cost: add. (8657->675), mult. (20246->879), div. (0->0), fcn. (15601->14), ass. (0->319)
t263 = sin(qJ(3));
t266 = cos(qJ(3));
t323 = pkin(3) * t263 - pkin(9) * t266;
t203 = t323 * qJD(3);
t262 = sin(qJ(4));
t264 = sin(qJ(2));
t265 = cos(qJ(4));
t372 = qJD(3) * t263;
t259 = sin(pkin(6));
t379 = qJD(1) * t259;
t267 = cos(qJ(2));
t392 = t266 * t267;
t446 = pkin(8) * t262;
t469 = t265 * t203 + t372 * t446 - (-t262 * t392 + t264 * t265) * t379;
t324 = pkin(3) * t266 + pkin(9) * t263;
t214 = -pkin(2) - t324;
t368 = qJD(4) * t265;
t396 = t262 * t264;
t468 = t262 * t203 + t214 * t368 - (t265 * t392 + t396) * t379;
t393 = t265 * t266;
t239 = pkin(8) * t393;
t310 = pkin(4) * t263 - qJ(5) * t393;
t367 = qJD(5) * t265;
t467 = -t263 * t367 + t310 * qJD(3) + (-t239 + (qJ(5) * t263 - t214) * t262) * qJD(4) + t469;
t394 = t263 * t265;
t466 = -(-pkin(8) * qJD(3) - qJ(5) * qJD(4)) * t394 - (-qJD(5) * t263 + (-pkin(8) * qJD(4) - qJ(5) * qJD(3)) * t266) * t262 - t468;
t465 = pkin(4) * t262 + pkin(8);
t366 = t265 * qJD(3);
t376 = qJD(2) * t263;
t199 = t262 * t376 - t366;
t373 = qJD(3) * t262;
t201 = t265 * t376 + t373;
t258 = sin(pkin(11));
t427 = cos(pkin(11));
t110 = t427 * t199 + t201 * t258;
t374 = qJD(2) * t266;
t236 = -qJD(4) + t374;
t416 = t110 * t236;
t364 = qJD(2) * qJD(3);
t343 = t266 * t364;
t362 = t263 * qJDD(2);
t369 = qJD(4) * t263;
t459 = qJD(2) * t369 - qJDD(3);
t103 = -qJD(4) * t366 + (-t343 - t362) * t265 + t459 * t262;
t104 = t262 * (qJD(3) * (qJD(4) + t374) + t362) + t459 * t265;
t45 = -t103 * t427 - t258 * t104;
t32 = t45 - t416;
t298 = -t258 * t199 + t201 * t427;
t464 = t110 * t298;
t335 = t427 * t262;
t195 = t258 * t265 + t335;
t180 = t195 * qJD(4);
t386 = t195 * t374 - t180;
t334 = t427 * t265;
t370 = qJD(4) * t262;
t181 = qJD(4) * t334 - t258 * t370;
t353 = t262 * t374;
t385 = -t258 * t353 + t334 * t374 - t181;
t365 = qJD(1) * qJD(2);
t345 = t267 * t365;
t260 = cos(pkin(6));
t378 = qJD(1) * t260;
t424 = qJDD(2) * pkin(8);
t463 = t424 + (qJDD(1) * t264 + t345) * t259 + qJD(3) * t378;
t371 = qJD(3) * t266;
t352 = t262 * t371;
t462 = t263 * t368 + t352;
t245 = pkin(4) * t265 + pkin(3);
t261 = -qJ(5) - pkin(9);
t398 = t261 * t263;
t461 = -t245 * t266 + t398;
t205 = qJD(2) * pkin(8) + t264 * t379;
t377 = qJD(1) * t263;
t149 = t205 * t266 + t260 * t377;
t138 = qJD(3) * pkin(9) + t149;
t356 = t267 * t379;
t151 = qJD(2) * t214 - t356;
t363 = qJDD(1) * t260;
t357 = -t263 * t363 - t266 * t463;
t71 = -t205 * t372 - t357;
t68 = qJDD(3) * pkin(9) + t71;
t346 = t264 * t365;
t399 = t259 * t267;
t319 = -qJDD(1) * t399 + t259 * t346;
t94 = qJD(2) * t203 + qJDD(2) * t214 + t319;
t301 = t138 * t370 - t151 * t368 - t262 * t94 - t265 * t68;
t80 = -t138 * t262 + t265 * t151;
t460 = t80 * t236 - t301;
t420 = t298 ^ 2;
t442 = t466 * t258 + t427 * t467;
t441 = t258 * t467 - t466 * t427;
t340 = qJD(4) * t261;
t173 = t262 * t340 + t367;
t292 = -qJD(5) * t262 + t265 * t340;
t192 = t263 * t205;
t148 = t266 * t378 - t192;
t202 = t323 * qJD(2);
t95 = -t148 * t262 + t265 * t202;
t78 = qJD(2) * t310 + t95;
t96 = t265 * t148 + t262 * t202;
t85 = -qJ(5) * t353 + t96;
t439 = (t292 - t78) * t427 + (-t173 + t85) * t258;
t82 = t138 * t265 + t151 * t262;
t17 = -qJD(4) * t82 - t262 * t68 + t265 * t94;
t458 = t82 * t236 - t17;
t327 = -t149 + (-t353 + t370) * pkin(4);
t311 = pkin(4) * t462 + pkin(8) * t371 - t263 * t356;
t162 = t262 * t214 + t239;
t457 = pkin(4) * t104 + qJDD(5);
t251 = t266 * qJDD(2);
t193 = t263 * t364 + qJDD(4) - t251;
t401 = t258 * t262;
t194 = -t334 + t401;
t456 = -t110 * t376 + t193 * t194 + t236 * t386;
t216 = t261 * t265;
t135 = -t216 * t427 + t261 * t401;
t255 = qJ(4) + pkin(11);
t247 = sin(t255);
t426 = sin(pkin(10));
t332 = t426 * t267;
t428 = cos(pkin(10));
t337 = t428 * t264;
t177 = t260 * t337 + t332;
t339 = t259 * t428;
t122 = t177 * t263 + t266 * t339;
t333 = t426 * t264;
t336 = t428 * t267;
t179 = -t260 * t333 + t336;
t338 = t259 * t426;
t124 = t179 * t263 - t266 * t338;
t400 = t259 * t264;
t183 = -t260 * t266 + t263 * t400;
t294 = g(1) * t124 + g(2) * t122 + g(3) * t183;
t455 = -t135 * t193 - t247 * t294;
t425 = qJDD(2) * pkin(2);
t159 = t319 - t425;
t268 = qJD(3) ^ 2;
t176 = -t260 * t336 + t333;
t178 = t260 * t332 + t337;
t322 = g(1) * t178 + g(2) * t176;
t454 = -pkin(8) * t268 + t259 * (-g(3) * t267 + t346) - t159 + t322 + t425;
t137 = -qJD(3) * pkin(3) - t148;
t101 = pkin(4) * t199 + qJD(5) + t137;
t38 = pkin(5) * t110 - qJ(6) * t298 + t101;
t452 = -t298 * t38 - qJDD(6);
t44 = -t103 * t258 + t427 * t104;
t451 = t110 * t385 - t194 * t45 - t195 * t44 + t298 * t386;
t11 = -qJ(5) * t104 - qJD(5) * t199 - t301;
t8 = pkin(4) * t193 + qJ(5) * t103 - qJD(5) * t201 + t17;
t3 = -t258 * t11 + t427 * t8;
t4 = t427 * t11 + t258 * t8;
t447 = pkin(5) * t193;
t444 = qJ(6) * t372 - qJD(6) * t266 + t441;
t443 = -pkin(5) * t372 - t442;
t60 = -qJ(5) * t201 + t80;
t47 = -pkin(4) * t236 + t60;
t61 = -qJ(5) * t199 + t82;
t55 = t427 * t61;
t23 = t258 * t47 + t55;
t440 = -pkin(5) * t386 + qJ(6) * t385 - qJD(6) * t195 + t327;
t37 = t258 * t78 + t427 * t85;
t438 = pkin(5) * t376 - t439;
t437 = qJD(2) * pkin(2);
t436 = t258 * t61;
t28 = t258 * t60 + t55;
t435 = t28 * t298;
t100 = t173 * t427 + t258 * t292;
t33 = qJ(6) * t376 + t37;
t432 = t100 - t33;
t431 = t100 - t37;
t430 = (-t263 * t366 - t266 * t370) * pkin(8) + t468;
t429 = -qJD(4) * t162 + t469;
t423 = qJDD(3) * pkin(3);
t422 = t103 * t262;
t421 = t104 * t265;
t419 = t298 * t236;
t417 = t110 ^ 2;
t123 = t177 * t266 - t263 * t339;
t414 = t123 * t262;
t125 = t179 * t266 + t263 * t338;
t413 = t125 * t262;
t411 = t176 * t265;
t410 = t178 * t265;
t409 = t199 * t236;
t408 = t199 * t262;
t407 = t201 * t199;
t406 = t201 * t236;
t405 = t201 * t265;
t403 = t247 * t266;
t248 = cos(t255);
t402 = t248 * t266;
t397 = t262 * t263;
t395 = t262 * t266;
t29 = t427 * t60 - t436;
t391 = qJD(6) - t29;
t390 = qJDD(1) - g(3);
t389 = -t122 * t245 - t123 * t261;
t388 = -t124 * t245 - t125 * t261;
t197 = t265 * t214;
t119 = -qJ(5) * t394 + t197 + (-pkin(4) - t446) * t266;
t131 = -qJ(5) * t397 + t162;
t63 = t258 * t119 + t427 * t131;
t184 = t260 * t263 + t266 * t400;
t387 = -t183 * t245 - t184 * t261;
t382 = pkin(2) * t399 + pkin(8) * t400;
t204 = pkin(4) * t397 + t263 * pkin(8);
t256 = t263 ^ 2;
t257 = t266 ^ 2;
t381 = t256 - t257;
t380 = t256 + t257;
t375 = qJD(2) * t264;
t360 = t265 * t399;
t359 = t247 * t399;
t269 = qJD(2) ^ 2;
t358 = t263 * t269 * t266;
t355 = t259 * t375;
t354 = qJD(2) * t399;
t351 = t266 * t366;
t349 = t236 * t376;
t347 = g(3) * t382;
t342 = t267 * t364;
t330 = t205 * t371 + t263 * t463 - t266 * t363;
t328 = t263 * t343;
t326 = t427 * t371;
t321 = g(1) * t179 + g(2) * t177;
t320 = -pkin(5) * t248 - qJ(6) * t247;
t169 = t195 * t263;
t97 = -t181 * t263 - t258 * t351 - t262 * t326;
t318 = -t110 * t97 + t169 * t44;
t316 = -t262 * t82 - t265 * t80;
t315 = -t417 - t420;
t314 = -t417 + t420;
t313 = qJDD(2) * t267 - t264 * t269;
t170 = -t258 * t397 + t263 * t334;
t98 = t180 * t263 + t258 * t352 - t265 * t326;
t312 = pkin(5) * t97 - qJ(6) * t98 + qJD(6) * t170 - t311;
t129 = -t184 * t262 - t360;
t307 = -t184 * t265 + t262 * t399;
t306 = t44 - t419;
t305 = -t44 - t419;
t22 = t427 * t47 - t436;
t303 = t193 * t262 - t236 * t368;
t302 = t193 * t265 + t236 * t370;
t126 = -qJD(3) * t183 + t266 * t354;
t53 = qJD(4) * t307 - t126 * t262 + t265 * t355;
t54 = qJD(4) * t129 + t126 * t265 + t262 * t355;
t25 = t258 * t54 - t427 * t53;
t27 = t258 * t53 + t427 * t54;
t64 = -t129 * t427 - t258 * t307;
t65 = t258 * t129 - t307 * t427;
t300 = -t27 * t110 + t25 * t298 - t65 * t44 + t45 * t64;
t299 = -t110 * t386 + t194 * t44;
t62 = t119 * t427 - t258 * t131;
t140 = -t248 * t400 + t266 * t359;
t90 = -t176 * t403 - t177 * t248;
t92 = -t178 * t403 - t179 * t248;
t297 = -g(1) * t92 - g(2) * t90 - g(3) * t140;
t141 = (t247 * t264 + t248 * t392) * t259;
t91 = -t176 * t402 + t177 * t247;
t93 = -t178 * t402 + t179 * t247;
t296 = -g(1) * t93 - g(2) * t91 - g(3) * t141;
t127 = qJD(3) * t184 + t263 * t354;
t295 = t127 * t110 + t183 * t44 - t193 * t64 + t236 * t25;
t69 = t330 - t423;
t293 = g(1) * t125 + g(2) * t123 + g(3) * t184;
t291 = t294 - t69;
t290 = g(3) * t399 - t322;
t289 = -g(3) * t400 - t321;
t114 = t184 * t247 + t248 * t399;
t73 = t123 * t247 - t176 * t248;
t75 = t125 * t247 - t178 * t248;
t288 = g(1) * t75 + g(2) * t73 + g(3) * t114 + t3;
t115 = t184 * t248 - t359;
t74 = t123 * t248 + t176 * t247;
t76 = t125 * t248 + t178 * t247;
t287 = g(1) * t76 + g(2) * t74 + g(3) * t115 - t4;
t286 = -pkin(9) * t193 - t137 * t236;
t285 = -t398 * t399 + t382 + (pkin(4) * t396 + t245 * t392) * t259;
t171 = t176 * pkin(2);
t284 = t461 * t176 + t177 * t465 - t171;
t172 = t178 * pkin(2);
t283 = t461 * t178 + t179 * t465 - t172;
t282 = -t45 - t416;
t281 = t290 * t263;
t280 = t110 * t98 - t169 * t45 - t170 * t44 + t298 * t97;
t279 = t127 * t298 + t183 * t45 - t193 * t65 + t236 * t27;
t40 = t69 + t457;
t277 = -t28 * t236 + t288;
t276 = t294 - t330;
t134 = -t216 * t258 - t261 * t335;
t275 = -t100 * t110 + t134 * t45 - t135 * t44 - t293;
t274 = t110 * t372 + t169 * t193 + t236 * t97 - t266 * t44;
t273 = pkin(9) * qJD(4) * t236 + t291;
t272 = -t134 * t193 + t248 * t294;
t206 = -t356 - t437;
t271 = -pkin(8) * qJDD(3) + (t206 + t356 - t437) * qJD(3);
t270 = t330 * t263 + t71 * t266 + (-t148 * t266 - t149 * t263) * qJD(3) - t321;
t5 = pkin(5) * t44 - qJ(6) * t45 - qJD(6) * t298 + t40;
t243 = -pkin(4) * t427 - pkin(5);
t240 = pkin(4) * t258 + qJ(6);
t189 = t193 * qJ(6);
t166 = pkin(4) * t410;
t164 = pkin(4) * t411;
t161 = -pkin(8) * t395 + t197;
t128 = -t193 * t266 - t236 * t372;
t107 = pkin(5) * t194 - qJ(6) * t195 - t245;
t84 = pkin(5) * t169 - qJ(6) * t170 + t204;
t58 = t266 * pkin(5) - t62;
t57 = -qJ(6) * t266 + t63;
t42 = pkin(4) * t201 + pkin(5) * t298 + qJ(6) * t110;
t35 = t193 * t195 + t236 * t385 - t298 * t376;
t20 = -qJ(6) * t236 + t23;
t18 = t236 * pkin(5) + qJD(6) - t22;
t15 = t170 * t45 - t298 * t98;
t13 = t195 * t45 - t298 * t385;
t12 = t170 * t193 + t236 * t98 - t266 * t45 + t298 * t372;
t2 = qJDD(6) - t447 - t3;
t1 = -qJD(6) * t236 + t189 + t4;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t390, 0, 0, 0, 0, 0, 0, t313 * t259 (-qJDD(2) * t264 - t267 * t269) * t259, 0, -g(3) + (t260 ^ 2 + (t264 ^ 2 + t267 ^ 2) * t259 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, -qJD(3) * t127 - qJDD(3) * t183 + (-t263 * t342 + t266 * t313) * t259, -qJD(3) * t126 - qJDD(3) * t184 + (-t263 * t313 - t266 * t342) * t259 (t183 * t263 + t184 * t266) * qJDD(2) + (t126 * t266 + t127 * t263 + (t183 * t266 - t184 * t263) * qJD(3)) * qJD(2), t126 * t149 - t127 * t148 + t183 * t330 + t184 * t71 - g(3) + (-t159 * t267 + t206 * t375) * t259, 0, 0, 0, 0, 0, 0, t104 * t183 + t127 * t199 + t129 * t193 - t236 * t53, -t103 * t183 + t127 * t201 + t193 * t307 + t236 * t54, t103 * t129 + t104 * t307 - t199 * t54 - t201 * t53, t127 * t137 + t129 * t17 + t183 * t69 + t301 * t307 + t53 * t80 + t54 * t82 - g(3), 0, 0, 0, 0, 0, 0, t295, t279, t300, t101 * t127 + t183 * t40 - t22 * t25 + t23 * t27 - t3 * t64 + t4 * t65 - g(3), 0, 0, 0, 0, 0, 0, t295, t300, -t279, t1 * t65 + t127 * t38 + t18 * t25 + t183 * t5 + t2 * t64 + t20 * t27 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t390 * t399 + t322, -t390 * t400 + t321, 0, 0, qJDD(2) * t256 + 0.2e1 * t328, 0.2e1 * t251 * t263 - 0.2e1 * t364 * t381, qJDD(3) * t263 + t266 * t268, qJDD(2) * t257 - 0.2e1 * t328, qJDD(3) * t266 - t263 * t268, 0, t271 * t263 + t266 * t454, -t263 * t454 + t271 * t266, t380 * t424 + (-g(3) * t264 - t345 * t380) * t259 + t270, -t159 * pkin(2) + g(1) * t172 + g(2) * t171 - t347 + (-t206 * t264 + (t148 * t263 - t149 * t266) * t267) * t379 + t270 * pkin(8), -t103 * t394 + (-t262 * t369 + t351) * t201 (-t199 * t265 - t201 * t262) * t371 + (t422 - t421 + (-t405 + t408) * qJD(4)) * t263 (-t236 * t366 + t103) * t266 + (qJD(3) * t201 + t302) * t263, t104 * t397 + t199 * t462 (t236 * t373 + t104) * t266 + (-qJD(3) * t199 - t303) * t263, t128, t161 * t193 - t429 * t236 + t289 * t262 + (-t17 + (pkin(8) * t199 + t137 * t262) * qJD(3) - t290 * t265) * t266 + (pkin(8) * t104 + qJD(3) * t80 + t137 * t368 - t199 * t356 + t262 * t69) * t263, -t162 * t193 + t430 * t236 + t289 * t265 + (-t301 + (pkin(8) * t201 + t137 * t265) * qJD(3) + t290 * t262) * t266 + (-pkin(8) * t103 - qJD(3) * t82 - t137 * t370 - t201 * t356 + t265 * t69) * t263, t103 * t161 - t104 * t162 - t429 * t201 - t430 * t199 + t316 * t371 + (t301 * t262 - t17 * t265 + (t262 * t80 - t265 * t82) * qJD(4) - t290) * t263, -t301 * t162 + t17 * t161 - g(1) * (-t178 * t324 - t172) - g(2) * (-t176 * t324 - t171) - t347 + t430 * t82 + t429 * t80 + (-g(3) * t324 - t137 * t377) * t399 + (t137 * t371 + t263 * t69 - t321) * pkin(8), t15, t280, t12, t318, -t274, t128, -t101 * t97 + t110 * t311 + t169 * t40 + t193 * t62 + t204 * t44 + t22 * t372 - t236 * t442 - t266 * t3 + t296, -t101 * t98 + t170 * t40 - t193 * t63 + t204 * t45 - t23 * t372 + t236 * t441 + t266 * t4 + t298 * t311 - t297, -t110 * t441 - t169 * t4 - t170 * t3 + t22 * t98 + t23 * t97 - t298 * t442 - t44 * t63 - t45 * t62 - t281, -g(1) * t283 - g(2) * t284 - g(3) * t285 + t101 * t311 + t40 * t204 + t22 * t442 + t23 * t441 + t3 * t62 + t4 * t63, t15, t12, -t280, t128, t274, t318, -t110 * t312 + t169 * t5 - t18 * t372 - t193 * t58 + t2 * t266 + t236 * t443 - t38 * t97 + t44 * t84 + t296, -t1 * t169 - t110 * t444 + t170 * t2 - t18 * t98 + t20 * t97 + t298 * t443 - t44 * t57 + t45 * t58 - t281, -t1 * t266 - t170 * t5 + t193 * t57 + t20 * t372 - t236 * t444 + t298 * t312 + t38 * t98 - t45 * t84 + t297, t1 * t57 + t5 * t84 + t2 * t58 - g(1) * (pkin(5) * t93 + qJ(6) * t92 + t283) - g(2) * (pkin(5) * t91 + qJ(6) * t90 + t284) - g(3) * (pkin(5) * t141 + qJ(6) * t140 + t285) - t312 * t38 + t444 * t20 + t443 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t358, t381 * t269, t362, t358, t251, qJDD(3), qJD(3) * t149 - t206 * t376 + t276, -t206 * t374 + (t148 + t192) * qJD(3) + t293 + t357, 0, 0, -t236 * t405 - t422 (-t103 + t409) * t265 + (-t104 + t406) * t262 (-t201 * t263 + t236 * t393) * qJD(2) + t303, -t236 * t408 - t421 (t199 * t263 - t236 * t395) * qJD(2) + t302, t349, -pkin(3) * t104 - t149 * t199 + t236 * t95 + t262 * t286 + t265 * t273 - t376 * t80, pkin(3) * t103 - t149 * t201 - t236 * t96 - t262 * t273 + t265 * t286 + t376 * t82, t199 * t96 + t201 * t95 + ((qJD(4) * t201 - t104) * pkin(9) + t460) * t265 + ((qJD(4) * t199 - t103) * pkin(9) + t458) * t262 - t293, -t137 * t149 - t80 * t95 - t82 * t96 + t291 * pkin(3) + (qJD(4) * t316 - t17 * t262 - t265 * t301 - t293) * pkin(9), t13, t451, t35, t299, -t456, t349, -t101 * t386 + t110 * t327 + t194 * t40 - t22 * t376 - t236 * t439 - t245 * t44 + t272, -t101 * t385 + t195 * t40 + t23 * t376 + t236 * t431 - t245 * t45 + t298 * t327 + t455, t110 * t37 - t194 * t4 - t195 * t3 + t22 * t385 + t23 * t386 - t298 * t439 + t275, -g(1) * t388 - g(2) * t389 - g(3) * t387 + t101 * t327 - t3 * t134 + t4 * t135 + t22 * t439 + t23 * t431 - t40 * t245, t13, t35, -t451, t349, t456, t299, t107 * t44 + t110 * t440 + t18 * t376 + t194 * t5 + t236 * t438 - t38 * t386 + t272, -t1 * t194 + t110 * t33 - t18 * t385 + t195 * t2 + t20 * t386 + t298 * t438 + t275, -t107 * t45 - t195 * t5 - t20 * t376 - t236 * t432 - t298 * t440 + t38 * t385 - t455, t1 * t135 + t5 * t107 + t2 * t134 - g(1) * (t124 * t320 + t388) - g(2) * (t122 * t320 + t389) - g(3) * (t183 * t320 + t387) + t440 * t38 + t432 * t20 + t438 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t407, -t199 ^ 2 + t201 ^ 2, -t103 - t409, -t407, -t104 - t406, t193, -t137 * t201 - g(1) * (t410 - t413) - g(2) * (t411 - t414) - g(3) * t129 - t458, t137 * t199 - g(1) * (-t125 * t265 - t178 * t262) - g(2) * (-t123 * t265 - t176 * t262) - g(3) * t307 - t460, 0, 0, t464, t314, t32, -t464, t305, t193, -t101 * t298 + (-t110 * t201 + t193 * t427) * pkin(4) + t277, t101 * t110 - t236 * t29 + (-t193 * t258 - t201 * t298) * pkin(4) + t287, t23 * t298 - t435 + (-t258 * t44 - t427 * t45) * pkin(4) + (t29 - t22) * t110, -g(1) * t166 - g(2) * t164 + t22 * t28 - t23 * t29 + (g(3) * t360 - t101 * t201 + t4 * t258 + t262 * t293 + t3 * t427) * pkin(4), t464, t32, -t314, t193, -t305, -t464, -t110 * t42 + (pkin(5) - t243) * t193 + t277 + t452, t20 * t298 - t240 * t44 + t243 * t45 - t435 + (t18 - t391) * t110, -t110 * t38 + t298 * t42 + t193 * t240 + t189 + (-0.2e1 * qJD(6) + t29) * t236 - t287, t1 * t240 + t2 * t243 - t38 * t42 - t18 * t28 - g(1) * (-pkin(4) * t413 - pkin(5) * t75 + qJ(6) * t76 + t166) - g(2) * (-pkin(4) * t414 - pkin(5) * t73 + qJ(6) * t74 + t164) - g(3) * (pkin(4) * t129 - pkin(5) * t114 + qJ(6) * t115) + t391 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t306, -t282, t315, t110 * t23 + t22 * t298 - t276 - t423 + t457, 0, 0, 0, 0, 0, 0, t306, t315, t282, t110 * t20 - t18 * t298 - t294 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t193 + t464, t32, -t236 ^ 2 - t420, t20 * t236 - t288 - t447 - t452;];
tau_reg  = t6;
