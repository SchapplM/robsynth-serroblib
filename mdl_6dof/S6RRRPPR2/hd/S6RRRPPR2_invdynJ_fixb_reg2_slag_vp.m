% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRRPPR2
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRPPR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:26:50
% EndTime: 2019-03-09 15:27:03
% DurationCPUTime: 7.22s
% Computational Cost: add. (14727->627), mult. (35188->746), div. (0->0), fcn. (26429->14), ass. (0->329)
t246 = qJD(2) + qJD(3);
t254 = sin(qJ(3));
t255 = sin(qJ(2));
t372 = t254 * t255;
t309 = t246 * t372;
t258 = cos(qJ(2));
t423 = cos(qJ(3));
t344 = t423 * t258;
t319 = qJD(1) * t344;
t336 = qJDD(1) * t423;
t351 = t258 * qJDD(1);
t322 = -t246 * t319 - t254 * t351 - t255 * t336;
t107 = qJD(1) * t309 + t322;
t182 = t254 * t258 + t255 * t423;
t138 = t246 * t182;
t352 = t255 * qJDD(1);
t308 = t254 * t352 - t258 * t336;
t108 = qJD(1) * t138 + t308;
t251 = sin(pkin(10));
t252 = cos(pkin(10));
t297 = t107 * t252 + t108 * t251;
t361 = qJD(1) * t255;
t343 = t254 * t361;
t167 = -t319 + t343;
t169 = t182 * qJD(1);
t123 = t252 * t167 + t169 * t251;
t388 = t123 * t246;
t55 = t297 - t388;
t419 = pkin(5) * t123;
t253 = sin(qJ(6));
t257 = cos(qJ(6));
t424 = pkin(8) + pkin(7);
t203 = t424 * t258;
t189 = qJD(1) * t203;
t174 = t423 * t189;
t202 = t424 * t255;
t187 = qJD(1) * t202;
t410 = qJD(2) * pkin(2);
t177 = -t187 + t410;
t129 = t254 * t177 + t174;
t384 = t167 * qJ(4);
t106 = t129 - t384;
t100 = t251 * t106;
t170 = t254 * t189;
t128 = t423 * t177 - t170;
t163 = t169 * qJ(4);
t105 = t128 - t163;
t99 = pkin(3) * t246 + t105;
t68 = t252 * t99 - t100;
t313 = qJD(5) - t68;
t294 = -t167 * t251 + t252 * t169;
t420 = pkin(5) * t294;
t425 = pkin(4) + pkin(9);
t47 = -t246 * t425 + t313 + t420;
t242 = t258 * pkin(2);
t234 = t242 + pkin(1);
t201 = t234 * qJD(1);
t141 = t167 * pkin(3) + qJD(4) - t201;
t270 = -qJ(5) * t294 + t141;
t58 = t123 * t425 + t270;
t21 = -t253 * t58 + t257 * t47;
t427 = qJD(6) + t294;
t441 = t21 * t427;
t135 = t187 * t254 - t174;
t112 = t135 + t384;
t136 = -t423 * t187 - t170;
t113 = -t163 + t136;
t222 = t251 * t254 * pkin(2);
t340 = qJD(3) * t423;
t398 = -qJD(3) * t222 - t112 * t251 + (pkin(2) * t340 - t113) * t252;
t329 = t253 * t427;
t73 = qJDD(6) - t297;
t70 = t257 * t73;
t440 = t427 * t329 - t70;
t109 = -t257 * t123 + t246 * t253;
t327 = t427 * t109;
t244 = qJDD(2) + qJDD(3);
t331 = t107 * t251 - t252 * t108;
t356 = qJD(6) * t257;
t357 = qJD(6) * t253;
t44 = -t123 * t356 - t257 * t244 + t246 * t357 + t253 * t331;
t436 = -t44 + t327;
t111 = t123 * t253 + t246 * t257;
t432 = t427 * t257;
t439 = t111 * t432;
t438 = t123 * t294;
t256 = sin(qJ(1));
t259 = cos(qJ(1));
t315 = g(1) * t259 + g(2) * t256;
t353 = qJD(1) * qJD(2);
t339 = t255 * t353;
t164 = pkin(2) * t339 - qJDD(1) * t234;
t91 = t108 * pkin(3) + qJDD(4) + t164;
t265 = qJ(5) * t297 - qJD(5) * t294 + t91;
t11 = -t331 * t425 + t265;
t22 = t253 * t47 + t257 * t58;
t338 = t258 * t353;
t140 = qJDD(2) * pkin(2) + t424 * (-t338 - t352);
t142 = t424 * (-t339 + t351);
t86 = -qJD(3) * t129 + t423 * t140 - t254 * t142;
t41 = t244 * pkin(3) + t107 * qJ(4) - t169 * qJD(4) + t86;
t360 = qJD(3) * t254;
t323 = -t254 * t140 - t423 * t142 - t177 * t340 + t189 * t360;
t46 = -qJ(4) * t108 - qJD(4) * t167 - t323;
t411 = t251 * t46 - t252 * t41;
t349 = qJDD(5) + t411;
t9 = -pkin(5) * t297 - t244 * t425 + t349;
t2 = -qJD(6) * t22 - t253 * t11 + t257 * t9;
t437 = t22 * t427 + t2;
t359 = qJD(6) * t111;
t45 = t253 * t244 + t257 * t331 + t359;
t435 = t111 * t427 - t45;
t393 = t294 ^ 2;
t16 = t251 * t41 + t252 * t46;
t429 = t244 * qJ(5) + t246 * qJD(5);
t325 = t16 + t429;
t10 = pkin(5) * t331 + t325;
t376 = t252 * t106;
t69 = t251 * t99 + t376;
t67 = -qJ(5) * t246 - t69;
t49 = -t67 - t419;
t434 = t10 * t253 + t49 * t356;
t400 = qJD(5) + t398;
t375 = t252 * t254;
t159 = (t251 * t423 + t375) * qJD(3) * pkin(2);
t82 = -t252 * t112 + t113 * t251;
t399 = t159 - t82;
t431 = t399 * t246;
t75 = t105 * t252 - t100;
t367 = qJD(5) - t75;
t144 = -t254 * t202 + t423 * t203;
t250 = qJ(2) + qJ(3);
t239 = pkin(10) + t250;
t227 = sin(t239);
t228 = cos(t239);
t430 = -t228 * pkin(4) - t227 * qJ(5);
t428 = g(1) * t256 - g(2) * t259;
t215 = g(3) * t228;
t381 = t227 * t259;
t382 = t227 * t256;
t347 = -g(1) * t381 - g(2) * t382 + t215;
t311 = -t347 - t411;
t79 = t123 * pkin(4) + t270;
t279 = t294 * t79 + qJDD(5) - t311;
t181 = -t344 + t372;
t345 = qJD(2) * t424;
t188 = t255 * t345;
t190 = t258 * t345;
t97 = -t423 * t188 - t254 * t190 - t202 * t340 - t203 * t360;
t80 = -qJ(4) * t138 - qJD(4) * t181 + t97;
t137 = -qJD(2) * t344 - t258 * t340 + t309;
t98 = -qJD(3) * t144 + t254 * t188 - t423 * t190;
t81 = t137 * qJ(4) - t182 * qJD(4) + t98;
t35 = t251 * t81 + t252 * t80;
t143 = -t423 * t202 - t203 * t254;
t120 = -qJ(4) * t182 + t143;
t121 = -qJ(4) * t181 + t144;
t90 = t120 * t251 + t121 * t252;
t426 = -t227 * t428 - t90 * t244 - t35 * t246;
t422 = pkin(3) * t169;
t240 = sin(t250);
t421 = pkin(3) * t240;
t214 = g(3) * t227;
t241 = cos(t250);
t415 = g(3) * t241;
t414 = g(3) * t258;
t212 = t228 * pkin(9);
t413 = t244 * pkin(4);
t245 = -qJ(4) - t424;
t412 = pkin(5) - t245;
t8 = t10 * t257;
t409 = t21 * t253;
t408 = t253 * t73;
t94 = -t137 * t251 + t252 * t138;
t407 = t253 * t94;
t406 = t257 * t44;
t405 = t257 * t94;
t404 = t45 * t253;
t74 = t105 * t251 + t376;
t403 = t74 * t294;
t402 = t74 * t246;
t401 = t420 + t400;
t397 = pkin(7) * qJDD(1);
t396 = t111 * t109;
t395 = t111 * t123;
t394 = t427 * t123;
t392 = t294 * t246;
t390 = t123 ^ 2;
t389 = t123 * t109;
t131 = t252 * t181 + t182 * t251;
t386 = t131 * t253;
t385 = t131 * t257;
t383 = t169 * t167;
t380 = t228 * t256;
t379 = t228 * t259;
t378 = t240 * t256;
t377 = t240 * t259;
t374 = t253 * t256;
t373 = t253 * t259;
t371 = t256 * t257;
t370 = t257 * t259;
t369 = t259 * t245;
t368 = t420 + t367;
t191 = -pkin(2) * t255 - t421;
t192 = qJ(5) * t380;
t366 = t256 * t191 + t192;
t194 = qJ(5) * t379;
t365 = t259 * t191 + t194;
t231 = pkin(3) * t241;
t364 = t231 + t242;
t248 = t255 ^ 2;
t249 = t258 ^ 2;
t363 = t248 - t249;
t362 = t248 + t249;
t358 = qJD(6) * t427;
t237 = t255 * t410;
t262 = qJD(1) ^ 2;
t350 = t255 * t262 * t258;
t348 = t231 - t430;
t229 = -pkin(3) * t252 - pkin(4);
t342 = -t22 * t294 - t2;
t341 = -t123 * t22 + t8;
t127 = t138 * pkin(3) + t237;
t34 = t251 * t80 - t252 * t81;
t335 = t399 * t294;
t333 = t427 * t49;
t332 = -t45 + t359;
t89 = -t252 * t120 + t121 * t251;
t233 = pkin(2) * t423 + pkin(3);
t161 = t233 * t252 - t222;
t330 = t427 ^ 2;
t324 = -t21 * t357 + t347;
t321 = t242 + t348;
t146 = t181 * pkin(3) - t234;
t186 = pkin(1) + t364;
t176 = t259 * t186;
t320 = g(2) * (pkin(4) * t379 + qJ(5) * t381 + t176);
t318 = t255 * t338;
t317 = g(1) * t380 - g(2) * t379;
t157 = -pkin(4) - t161;
t316 = -pkin(4) * t227 - t421;
t312 = t294 * t409 - t324;
t310 = g(1) * t379 + g(2) * t380 - t16 + t214;
t66 = -pkin(4) * t246 + t313;
t307 = t123 * t66 - t294 * t67;
t306 = -t123 * t68 + t294 * t69;
t305 = t123 * t94 - t131 * t331;
t132 = -t181 * t251 + t182 * t252;
t95 = -t137 * t252 - t138 * t251;
t304 = -t132 * t297 + t294 * t95;
t303 = t21 * t257 + t22 * t253;
t302 = t22 * t257 - t409;
t287 = -t132 * qJ(5) + t146;
t63 = t131 * t425 + t287;
t64 = pkin(5) * t132 + t89;
t33 = t253 * t64 + t257 * t63;
t32 = -t253 * t63 + t257 * t64;
t300 = t257 * t294 * t49 + t123 * t21 + t434;
t299 = t131 * t244 + t246 * t94;
t298 = t132 * t244 + t246 * t95;
t296 = -t390 - t393;
t295 = -t390 + t393;
t162 = pkin(2) * t375 + t233 * t251;
t293 = -t186 + t430;
t291 = -0.2e1 * pkin(1) * t353 - pkin(7) * qJDD(2);
t290 = t131 * t356 + t407;
t289 = t131 * t357 - t405;
t288 = t331 - t392;
t53 = t331 + t392;
t87 = pkin(4) * t294 + qJ(5) * t123 + t422;
t284 = -t141 * t294 + t311;
t283 = t123 * t141 + t310;
t282 = -t427 * t432 - t408;
t281 = t297 + t388;
t236 = pkin(2) * t361;
t84 = t236 + t87;
t280 = -t95 * qJ(5) - t132 * qJD(5) + t127;
t278 = -t228 * t315 - t214;
t277 = -t244 * t89 - t246 * t34 + t317;
t276 = -t123 * t79 - t310 + t429;
t261 = qJD(2) ^ 2;
t275 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t261 + t428;
t274 = pkin(1) * t262 + t315 - t397;
t273 = g(3) * t240 - t201 * t167 + t241 * t315 + t323;
t272 = -t123 * t95 + t131 * t297 + t132 * t331 - t294 * t94;
t271 = t315 * t425 * t227;
t1 = qJD(6) * t21 + t257 * t11 + t253 * t9;
t269 = qJD(6) * t302 + t1 * t253 + t2 * t257;
t149 = -pkin(9) + t157;
t268 = -t149 * t358 + t278;
t224 = -pkin(9) + t229;
t267 = -t224 * t358 + t278;
t266 = -t123 * t35 + t294 * t34 - t297 * t89 + t331 * t90 - t315;
t17 = -pkin(4) * t331 + t265;
t264 = g(1) * t377 + g(2) * t378 + t201 * t169 - t415 + t86;
t226 = pkin(3) * t251 + qJ(5);
t156 = qJ(5) + t162;
t155 = -t227 * t374 + t370;
t154 = t227 * t371 + t373;
t153 = t227 * t373 + t371;
t152 = t227 * t370 - t374;
t145 = t236 + t422;
t118 = t294 * pkin(9);
t114 = -t167 ^ 2 + t169 ^ 2;
t92 = -t322 + (t167 - t343) * t246;
t88 = pkin(4) * t131 + t287;
t65 = -pkin(5) * t131 + t90;
t62 = t118 + t87;
t61 = t118 + t84;
t59 = t82 - t419;
t56 = t74 - t419;
t40 = t257 * t45;
t36 = pkin(4) * t94 + t280;
t31 = t425 * t94 + t280;
t30 = -pkin(5) * t94 + t35;
t29 = pkin(5) * t95 + t34;
t28 = t253 * t59 + t257 * t61;
t27 = -t253 * t61 + t257 * t59;
t26 = t253 * t56 + t257 * t62;
t25 = -t253 * t62 + t257 * t56;
t24 = t109 * t432 + t404;
t23 = -t111 * t329 - t406;
t19 = t282 - t389;
t18 = t395 - t440;
t14 = t349 - t413;
t5 = -t40 - t439 + (t44 + t327) * t253;
t4 = -qJD(6) * t33 - t253 * t31 + t257 * t29;
t3 = qJD(6) * t32 + t253 * t29 + t257 * t31;
t6 = [0, 0, 0, 0, 0, qJDD(1), t428, t315, 0, 0, qJDD(1) * t248 + 0.2e1 * t318, 0.2e1 * t255 * t351 - 0.2e1 * t353 * t363, qJDD(2) * t255 + t258 * t261, qJDD(1) * t249 - 0.2e1 * t318, qJDD(2) * t258 - t255 * t261, 0, t255 * t291 + t258 * t275, -t255 * t275 + t258 * t291, 0.2e1 * t362 * t397 - t315, -g(1) * (-pkin(1) * t256 + pkin(7) * t259) - g(2) * (pkin(1) * t259 + pkin(7) * t256) + (pkin(7) ^ 2 * t362 + pkin(1) ^ 2) * qJDD(1), -t107 * t182 - t137 * t169, t107 * t181 - t108 * t182 + t137 * t167 - t138 * t169, -t137 * t246 + t182 * t244, t108 * t181 + t138 * t167, -t138 * t246 - t181 * t244, 0, -t234 * t108 - t201 * t138 + t143 * t244 + t164 * t181 + t167 * t237 + t241 * t428 + t98 * t246, t234 * t107 + t201 * t137 - t144 * t244 + t164 * t182 + t169 * t237 - t240 * t428 - t97 * t246, t107 * t143 - t108 * t144 + t128 * t137 - t129 * t138 - t167 * t97 - t169 * t98 + t181 * t323 - t182 * t86 - t315, -t323 * t144 + t129 * t97 + t86 * t143 + t128 * t98 - t164 * t234 - t201 * t237 - g(1) * (-t234 * t256 + t259 * t424) - g(2) * (t234 * t259 + t256 * t424) t304, t272, t298, t305, -t299, 0, t123 * t127 + t131 * t91 + t141 * t94 - t146 * t331 + t277, t127 * t294 + t91 * t132 + t141 * t95 - t146 * t297 + t426, -t131 * t16 + t132 * t411 - t68 * t95 - t69 * t94 + t266, t16 * t90 + t69 * t35 + t411 * t89 - t68 * t34 + t91 * t146 + t141 * t127 - g(1) * (-t256 * t186 - t369) - g(2) * (-t245 * t256 + t176) 0, -t298, t299, t304, t272, t305, -t131 * t325 + t132 * t14 + t66 * t95 + t67 * t94 + t266, -t123 * t36 - t131 * t17 + t331 * t88 - t79 * t94 - t277, -t17 * t132 - t294 * t36 + t297 * t88 - t79 * t95 - t426, t17 * t88 + t79 * t36 + t325 * t90 - t67 * t35 + t14 * t89 + t66 * t34 + g(1) * t369 - t320 + (-g(1) * t293 + g(2) * t245) * t256, t111 * t290 - t386 * t44 (-t109 * t253 + t111 * t257) * t94 + (-t404 - t406 + (-t109 * t257 - t111 * t253) * qJD(6)) * t131, t111 * t95 - t44 * t132 + t290 * t427 + t386 * t73, t109 * t289 - t385 * t45, -t109 * t95 - t45 * t132 - t289 * t427 + t385 * t73, t132 * t73 + t427 * t95, -t49 * t405 - g(1) * t155 - g(2) * t153 + t30 * t109 + t4 * t427 + t2 * t132 + t21 * t95 + t32 * t73 + t65 * t45 + (t357 * t49 - t8) * t131, g(1) * t154 - g(2) * t152 - t1 * t132 + t30 * t111 + t131 * t434 - t22 * t95 - t3 * t427 - t33 * t73 + t49 * t407 - t65 * t44, -t3 * t109 - t4 * t111 + t32 * t44 - t33 * t45 + t302 * t94 + (-qJD(6) * t303 + t1 * t257 - t2 * t253) * t131 + t317, t1 * t33 + t22 * t3 + t2 * t32 + t21 * t4 + t10 * t65 + t49 * t30 - t320 + (-g(1) * t412 - g(2) * t212) * t259 + (-g(1) * (t293 - t212) - g(2) * t412) * t256; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t350, t363 * t262, t352, t350, t351, qJDD(2), t255 * t274 - t414, g(3) * t255 + t258 * t274, 0, 0, t383, t114, t92, -t383, -t308, t244, -t135 * t246 + (-t167 * t361 + t244 * t423 - t246 * t360) * pkin(2) + t264, t136 * t246 + (-t169 * t361 - t244 * t254 - t246 * t340) * pkin(2) + t273 (t129 + t135) * t169 + (-t128 + t136) * t167 + (t423 * t107 - t108 * t254 + (-t167 * t423 + t169 * t254) * qJD(3)) * pkin(2), -t128 * t135 - t129 * t136 + (t423 * t86 - t414 - t254 * t323 + (-t128 * t254 + t129 * t423) * qJD(3) + (qJD(1) * t201 + t315) * t255) * pkin(2), t438, t295, -t55, -t438, t53, t244, -t145 * t123 + t161 * t244 + t284 - t431, -t145 * t294 - t162 * t244 - t246 * t398 + t283, -t123 * t398 + t161 * t297 + t162 * t331 + t306 + t335, -g(3) * t364 - t141 * t145 + t16 * t162 - t161 * t411 - t191 * t315 + t398 * t69 - t399 * t68, t244, t55, -t53, t438, t295, -t438, -t123 * t400 + t156 * t331 - t157 * t297 + t307 + t335, t84 * t123 + t431 + (-pkin(4) + t157) * t244 + t279, t156 * t244 + t246 * t400 + t294 * t84 + t276, t325 * t156 + t14 * t157 - t79 * t84 - g(1) * (-pkin(4) * t381 + t365) - g(2) * (-pkin(4) * t382 + t366) - g(3) * t321 - t400 * t67 + t399 * t66, t23, t5, t18, t24, t19, t394, t149 * t70 + t156 * t45 + (t159 * t257 - t27) * t427 + t401 * t109 + t268 * t253 + t300, t28 * t427 - t156 * t44 + t401 * t111 + t268 * t257 + (-t149 * t73 - t159 * t427 - t333) * t253 + t341, t28 * t109 + t27 * t111 + (-t109 * t159 + t149 * t332 - t1) * t253 + (-t111 * t159 + t149 * t44 + (-t109 * t149 - t22) * qJD(6) + t342) * t257 + t312, t10 * t156 - t22 * t28 - t21 * t27 - g(1) * t365 - g(2) * t366 - g(3) * (t212 + t321) + t401 * t49 + t271 + t303 * t159 + t269 * t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t383, t114, t92, -t383, -t308, t244, t129 * t246 + t264, t128 * t246 + t273, 0, 0, t438, t295, -t55, -t438, t53, t244, t402 + (-t123 * t169 + t244 * t252) * pkin(3) + t284, t75 * t246 + (-t169 * t294 - t244 * t251) * pkin(3) + t283, t75 * t123 - t403 + (t251 * t331 + t252 * t297) * pkin(3) + t306, t68 * t74 - t69 * t75 + (-t141 * t169 + t16 * t251 + t240 * t315 - t252 * t411 - t415) * pkin(3), t244, t55, -t53, t438, t295, -t438, -t123 * t367 + t226 * t331 - t229 * t297 + t307 - t403, t87 * t123 - t402 + (-pkin(4) + t229) * t244 + t279, t226 * t244 + t246 * t367 + t294 * t87 + t276, t325 * t226 + t14 * t229 - t79 * t87 - t66 * t74 - g(1) * (t259 * t316 + t194) - g(2) * (t256 * t316 + t192) - g(3) * t348 - t367 * t67, t23, t5, t18, t24, t19, t394, t109 * t368 + t224 * t70 + t226 * t45 - t25 * t427 + t253 * t267 + t300, t26 * t427 - t226 * t44 + t368 * t111 + (-t224 * t73 - t333) * t253 + t267 * t257 + t341, t26 * t109 + t25 * t111 + (t224 * t332 - t1) * t253 + (t224 * t44 + (-t109 * t224 - t22) * qJD(6) + t342) * t257 + t312, t10 * t226 - t22 * t26 - t21 * t25 - g(1) * (-pkin(3) * t377 + t194) - g(2) * (-pkin(3) * t378 + t192) - g(3) * (t212 + t348) + t368 * t49 + t271 + t269 * t224; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t288, -t281, t296, t123 * t69 + t294 * t68 - t428 + t91, 0, 0, 0, 0, 0, 0, t296, t288, t281, -t123 * t67 - t294 * t66 + t17 - t428, 0, 0, 0, 0, 0, 0, t282 + t389, t395 + t440, t253 * t436 - t40 + t439, t49 * t123 + (t1 - t441) * t257 - t437 * t253 - t428; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, t244 - t438, -t246 ^ 2 - t393, t246 * t67 + t279 - t413, 0, 0, 0, 0, 0, 0, -t246 * t109 - t253 * t330 + t70, -t246 * t111 - t257 * t330 - t408, t253 * t435 - t257 * t436, -t49 * t246 + (-t21 * t294 + t1) * t253 + t437 * t257 + t324; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t396, -t109 ^ 2 + t111 ^ 2, t436, -t396, t435, t73, -g(1) * t152 - g(2) * t154 - t49 * t111 + t215 * t257 + t437, g(1) * t153 - g(2) * t155 + t49 * t109 + t441 + (-qJD(6) * t47 - t11) * t257 + (qJD(6) * t58 - t215 - t9) * t253, 0, 0;];
tau_reg  = t6;