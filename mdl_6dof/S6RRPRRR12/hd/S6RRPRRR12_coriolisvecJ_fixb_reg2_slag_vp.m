% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RRPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRRR12_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR12_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:41:43
% EndTime: 2019-03-09 14:42:12
% DurationCPUTime: 12.61s
% Computational Cost: add. (20807->676), mult. (51477->919), div. (0->0), fcn. (38837->10), ass. (0->296)
t281 = sin(qJ(2));
t276 = sin(pkin(6));
t362 = qJD(1) * t276;
t338 = t281 * t362;
t255 = pkin(2) * t338;
t284 = cos(qJ(2));
t308 = pkin(9) * t281 - qJ(3) * t284;
t178 = t308 * t362 + t255;
t337 = t284 * t362;
t277 = cos(pkin(6));
t361 = qJD(1) * t277;
t344 = pkin(1) * t361;
t205 = pkin(8) * t337 + t281 * t344;
t182 = pkin(3) * t337 + t205;
t280 = sin(qJ(4));
t283 = cos(qJ(4));
t112 = t283 * t178 + t280 * t182;
t285 = -pkin(2) - pkin(9);
t407 = pkin(10) - t285;
t231 = t407 * t283;
t316 = t283 * t338;
t449 = -pkin(10) * t316 - qJD(4) * t231 - t112;
t111 = -t178 * t280 + t283 * t182;
t357 = qJD(4) * t280;
t368 = t280 * t281;
t448 = t407 * t357 - (pkin(4) * t284 - pkin(10) * t368) * t362 - t111;
t282 = cos(qJ(6));
t352 = qJD(6) * t282;
t259 = qJD(2) + t361;
t315 = t280 * t337;
t194 = t259 * t283 - t315;
t279 = sin(qJ(5));
t295 = t280 * t259 + t283 * t337;
t410 = cos(qJ(5));
t119 = t194 * t279 + t295 * t410;
t437 = t119 * t282;
t447 = t352 + t437;
t219 = t279 * t283 + t280 * t410;
t347 = qJD(4) + qJD(5);
t367 = (t338 + t347) * t219;
t332 = t410 * qJD(5);
t355 = qJD(5) * t279;
t369 = t279 * t280;
t366 = -t279 * t357 - t280 * t355 + (qJD(4) * t410 + t332) * t283 + t316 * t410 - t338 * t369;
t253 = t284 * t344;
t422 = qJD(3) - t253;
t278 = sin(qJ(6));
t244 = qJD(4) + t338;
t302 = qJD(5) + t244;
t292 = t282 * t302;
t301 = t194 * t410 - t279 * t295;
t103 = t278 * t301 - t292;
t105 = t278 * t302 + t282 * t301;
t348 = qJD(1) * qJD(2);
t330 = t284 * t348;
t246 = t276 * t330;
t353 = qJD(6) * t278;
t360 = qJD(2) * t281;
t336 = t276 * t360;
t314 = qJD(1) * t336;
t140 = qJD(4) * t295 - t280 * t314;
t356 = qJD(4) * t283;
t141 = -qJD(4) * t315 + t259 * t356 - t283 * t314;
t71 = t140 * t410 + t279 * t141 + t194 * t355 + t295 * t332;
t40 = -qJD(6) * t292 - t278 * t246 + t282 * t71 + t301 * t353;
t354 = qJD(6) * t105;
t41 = -t282 * t246 - t278 * t71 + t354;
t436 = qJD(6) + t119;
t444 = t436 * t278;
t446 = -t103 * t447 - t105 * t444 - t278 * t41 - t40 * t282;
t321 = t279 * t140 - t410 * t141;
t72 = qJD(5) * t301 - t321;
t70 = t282 * t72;
t445 = t103 * t301 - t436 * t444 + t70;
t230 = t407 * t280;
t300 = t279 * t230 - t231 * t410;
t392 = qJD(5) * t300 + t448 * t279 + t410 * t449;
t412 = pkin(3) + pkin(8);
t364 = pkin(4) * t356 - (-pkin(4) * t283 - t412) * t338 + t422;
t37 = t40 * t278;
t442 = t105 * t447 - t37;
t68 = t278 * t72;
t389 = t352 * t436 + t68;
t441 = -t105 * t301 + t436 * t437 + t389;
t125 = t259 * t285 + t412 * t338 + t422;
t329 = -qJ(3) * t281 - pkin(1);
t176 = (t284 * t285 + t329) * t276;
t150 = qJD(1) * t176;
t91 = t125 * t280 + t150 * t283;
t81 = -pkin(10) * t295 + t91;
t393 = t279 * t81;
t90 = t283 * t125 - t150 * t280;
t80 = -pkin(10) * t194 + t90;
t76 = pkin(4) * t244 + t80;
t32 = t410 * t76 - t393;
t30 = -pkin(5) * t302 - t32;
t440 = t119 * t30;
t439 = -t366 * pkin(5) - pkin(11) * t367 - t364;
t438 = -pkin(11) * t337 + t392;
t382 = t119 * t301;
t435 = -t119 ^ 2 + t301 ^ 2;
t84 = pkin(5) * t301 + pkin(11) * t119;
t434 = t119 * t302 - t71;
t247 = t259 * qJ(3);
t143 = t247 + t182;
t113 = pkin(4) * t295 + t143;
t241 = pkin(2) * t314;
t358 = qJD(3) * t281;
t288 = (qJD(2) * t308 - t358) * t276;
t133 = qJD(1) * t288 + t241;
t409 = pkin(1) * t281;
t264 = t277 * t409;
t370 = t276 * t284;
t183 = (t370 * t412 + t264) * qJD(2);
t157 = qJD(1) * t183;
t62 = -qJD(4) * t91 - t133 * t280 + t283 * t157;
t46 = pkin(4) * t246 + pkin(10) * t140 + t62;
t61 = t125 * t356 + t283 * t133 - t150 * t357 + t280 * t157;
t48 = -pkin(10) * t141 + t61;
t327 = -t279 * t46 - t76 * t332 + t81 * t355 - t410 * t48;
t433 = t113 * t119 + t327;
t273 = t276 ^ 2;
t431 = -0.2e1 * t273 * t348;
t343 = t410 * t81;
t33 = t279 * t76 + t343;
t31 = pkin(11) * t302 + t33;
t63 = pkin(5) * t119 - pkin(11) * t301 + t113;
t16 = t278 * t63 + t282 * t31;
t365 = pkin(8) * t314 - qJD(2) * t253;
t163 = -t259 * qJD(3) + t365;
t134 = -pkin(3) * t314 - t163;
t95 = pkin(4) * t141 + t134;
t22 = pkin(5) * t72 + pkin(11) * t71 + t95;
t7 = pkin(11) * t246 - t327;
t3 = -qJD(6) * t16 + t282 * t22 - t278 * t7;
t429 = -t16 * t436 - t3;
t171 = -t230 * t410 - t279 * t231;
t390 = qJD(5) * t171 + t279 * t449 - t448 * t410;
t428 = -t244 * t90 + t61;
t427 = t244 * t91 + t62;
t426 = t436 * t301;
t381 = t295 * t244;
t425 = t140 - t381;
t379 = t194 * t244;
t424 = -t141 + t379;
t371 = t276 * t281;
t261 = pkin(8) * t371;
t408 = pkin(1) * t284;
t339 = -pkin(2) - t408;
t156 = pkin(3) * t371 + t261 + (-pkin(9) + t339) * t277;
t107 = t280 * t156 + t283 * t176;
t306 = t278 * t31 - t282 * t63;
t423 = t16 * t282 + t278 * t306;
t328 = t279 * t48 - t410 * t46;
t10 = -qJD(5) * t33 - t328;
t8 = -pkin(5) * t246 - t10;
t420 = -t8 * t282 + t30 * t353 + t301 * t306;
t419 = t16 * t301 + t8 * t278 + t30 * t352;
t417 = -t113 * t301 - t328;
t397 = t219 * t72;
t416 = t119 * t366 + t397;
t415 = t244 * t301 + t321;
t220 = t283 * t410 - t369;
t414 = -t220 * t71 - t301 * t367;
t106 = t283 * t156 - t176 * t280;
t211 = t277 * t283 - t280 * t370;
t88 = pkin(4) * t371 - pkin(10) * t211 + t106;
t210 = t277 * t280 + t283 * t370;
t93 = -pkin(10) * t210 + t107;
t404 = t279 * t88 + t410 * t93;
t168 = qJD(4) * t210 - t280 * t336;
t359 = qJD(2) * t284;
t335 = t276 * t359;
t250 = pkin(2) * t336;
t146 = t250 + t288;
t67 = -qJD(4) * t107 - t146 * t280 + t283 * t183;
t54 = pkin(4) * t335 + pkin(10) * t168 + t67;
t169 = qJD(4) * t211 - t283 * t336;
t66 = t283 * t146 + t156 * t356 - t176 * t357 + t280 * t183;
t60 = -pkin(10) * t169 + t66;
t14 = -qJD(5) * t404 - t279 * t60 + t410 * t54;
t2 = -qJD(6) * t306 + t22 * t278 + t282 * t7;
t413 = t10 * t220 - t219 * t327 - t32 * t367 + t33 * t366;
t1 = t2 * t282;
t265 = t280 * pkin(4) + qJ(3);
t147 = pkin(5) * t219 - pkin(11) * t220 + t265;
t102 = t147 * t278 + t171 * t282;
t406 = qJD(6) * t102 + t278 * t438 + t282 * t439;
t101 = t147 * t282 - t171 * t278;
t405 = -qJD(6) * t101 + t278 * t439 - t282 * t438;
t402 = t436 * t306;
t136 = t210 * t410 + t211 * t279;
t400 = t136 * t72;
t39 = t41 * t282;
t391 = pkin(5) * t337 + t390;
t388 = t103 * t278;
t387 = t105 * t103;
t386 = t105 * t282;
t380 = t194 * t295;
t214 = pkin(8) * t370 + t264;
t207 = t214 * qJD(2);
t197 = qJD(1) * t207;
t378 = t197 * t281;
t377 = t205 * t259;
t376 = t207 * t281;
t375 = t220 * t278;
t374 = t220 * t282;
t373 = t244 * t285;
t372 = t273 * qJD(1) ^ 2;
t274 = t281 ^ 2;
t363 = -t284 ^ 2 + t274;
t317 = t412 * t371;
t351 = qJD(1) * t317 + t422;
t204 = pkin(8) * t338 - t253;
t350 = qJD(3) + t204;
t345 = t277 * t408;
t342 = t244 * t281 * t283;
t341 = t284 * t372;
t340 = t282 * t371;
t199 = -t277 * qJ(3) - t214;
t334 = t244 * t356;
t325 = t259 * t282 + t278 * t366;
t324 = -t259 * t278 + t282 * t366;
t323 = t278 * t367 - t282 * t337;
t322 = t278 * t337 + t282 * t367;
t318 = pkin(4) * t332;
t175 = pkin(3) * t370 - t199;
t34 = t279 * t80 + t343;
t313 = pkin(4) * t355 - t34;
t311 = pkin(1) * t431;
t267 = pkin(4) * t279 + pkin(11);
t307 = -t267 * t72 + t440;
t50 = pkin(11) * t371 + t404;
t126 = pkin(4) * t210 + t175;
t137 = -t279 * t210 + t211 * t410;
t74 = pkin(5) * t136 - pkin(11) * t137 + t126;
t24 = t278 * t74 + t282 * t50;
t23 = -t278 * t50 + t282 * t74;
t305 = t197 * t277 + t207 * t259;
t254 = qJD(2) * t345;
t206 = -pkin(8) * t336 + t254;
t189 = -t259 * t337 + t246;
t51 = -t279 * t93 + t410 * t88;
t116 = t137 * t282 + t278 * t371;
t13 = t279 * t54 + t88 * t332 - t355 * t93 + t410 * t60;
t299 = t143 * t281 + t285 * t359;
t296 = qJD(2) * t276 * (t259 + t361);
t200 = (-pkin(2) * t284 + t329) * t276;
t294 = t220 * t352 - t323;
t293 = -t220 * t353 - t322;
t270 = t277 * qJD(3);
t155 = -qJD(2) * t317 + t254 + t270;
t291 = (-qJ(3) * t359 - t358) * t276;
t290 = t302 * t370;
t110 = pkin(4) * t169 + t155;
t289 = t220 * t246 - t302 * t367;
t287 = t1 - t3 * t278 + (-t16 * t278 + t282 * t306) * qJD(6);
t268 = -pkin(4) * t410 - pkin(5);
t234 = t281 * t341;
t228 = t283 * t246;
t225 = t273 * t281 * t330;
t222 = -0.2e1 * t225;
t221 = 0.2e1 * t225;
t213 = -t261 + t345;
t209 = t363 * t372;
t203 = -qJ(3) * t337 + t255;
t201 = t277 * t339 + t261;
t195 = t363 * t431;
t191 = -t206 - t270;
t190 = qJD(1) * t200;
t188 = (qJD(2) - t259) * t338;
t185 = t250 + t291;
t177 = -t247 - t205;
t174 = t284 * t296;
t173 = t281 * t296;
t172 = -pkin(2) * t259 + t350;
t161 = qJD(1) * t291 + t241;
t154 = t190 * t338;
t115 = t137 * t278 - t340;
t83 = qJD(5) * t137 - t279 * t168 + t169 * t410;
t82 = t168 * t410 + t279 * t169 + t210 * t332 + t211 * t355;
t73 = pkin(4) * t194 + t84;
t56 = qJD(6) * t116 - t278 * t82 - t282 * t335;
t55 = -qJD(6) * t340 + t137 * t353 - t278 * t335 + t282 * t82;
t49 = -pkin(5) * t371 - t51;
t35 = t410 * t80 - t393;
t25 = pkin(5) * t83 + pkin(11) * t82 + t110;
t21 = t278 * t84 + t282 * t32;
t20 = -t278 * t32 + t282 * t84;
t18 = t278 * t73 + t282 * t35;
t17 = -t278 * t35 + t282 * t73;
t12 = -pkin(5) * t335 - t14;
t11 = pkin(11) * t335 + t13;
t5 = -qJD(6) * t24 - t11 * t278 + t25 * t282;
t4 = qJD(6) * t23 + t11 * t282 + t25 * t278;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t221, t195, t174, t222, -t173, 0, t281 * t311 - t305, -t206 * t259 + t277 * t365 + t284 * t311 (-t365 * t284 + t378 + (t204 * t284 - t205 * t281) * qJD(2) + (t206 * t284 + t376 + (-t213 * t284 - t214 * t281) * qJD(2)) * qJD(1)) * t276, -t197 * t213 + t204 * t207 + t205 * t206 - t214 * t365, 0, -t174, t173, t221, t195, t222 (-t163 * t284 + t378 + (t172 * t284 + t177 * t281) * qJD(2) + (-t191 * t284 + t376 + (t199 * t281 + t201 * t284) * qJD(2)) * qJD(1)) * t276 (-t190 * t360 + t161 * t284 + (t185 * t284 - t200 * t360) * qJD(1)) * t276 + t305, -t163 * t277 - t191 * t259 + (-t190 * t359 - t161 * t281 + (-t185 * t281 - t200 * t359) * qJD(1)) * t276, t161 * t200 + t163 * t199 + t172 * t207 + t177 * t191 + t185 * t190 + t197 * t201, -t140 * t211 - t168 * t194, t140 * t210 - t141 * t211 + t168 * t295 - t169 * t194, -t168 * t244 + (-t140 * t281 + (qJD(1) * t211 + t194) * t359) * t276, t141 * t210 + t169 * t295, -t169 * t244 + (-t141 * t281 + (-qJD(1) * t210 - t295) * t359) * t276, t244 * t335 + t225, t134 * t210 + t141 * t175 + t143 * t169 + t155 * t295 + t244 * t67 + (t281 * t62 + (qJD(1) * t106 + t90) * t359) * t276, t134 * t211 - t140 * t175 - t143 * t168 + t155 * t194 - t244 * t66 + (-t281 * t61 + (-qJD(1) * t107 - t91) * t359) * t276, t106 * t140 - t107 * t141 + t168 * t90 - t169 * t91 - t194 * t67 - t210 * t61 - t211 * t62 - t295 * t66, t106 * t62 + t107 * t61 + t134 * t175 + t143 * t155 + t66 * t91 + t67 * t90, -t137 * t71 - t301 * t82, t119 * t82 + t136 * t71 - t137 * t72 - t301 * t83, -t82 * t347 + (t301 * t359 - t71 * t281 + (t137 * t359 - t281 * t82) * qJD(1)) * t276, t119 * t83 + t400, -t83 * t347 + (-t119 * t359 - t72 * t281 + (-t136 * t359 - t281 * t83) * qJD(1)) * t276, qJD(2) * t290 + t225, t14 * t347 + t110 * t119 + t126 * t72 + t95 * t136 + t113 * t83 + (t32 * t359 + t10 * t281 + (t14 * t281 + t359 * t51) * qJD(1)) * t276, -t13 * t347 + t110 * t301 - t126 * t71 + t95 * t137 - t113 * t82 + (-t33 * t359 + t327 * t281 + (-t13 * t281 - t359 * t404) * qJD(1)) * t276, -t10 * t137 - t119 * t13 + t136 * t327 - t14 * t301 + t32 * t82 - t33 * t83 - t404 * t72 + t51 * t71, t10 * t51 + t110 * t113 + t126 * t95 + t13 * t33 + t14 * t32 - t327 * t404, -t105 * t55 - t116 * t40, t103 * t55 - t105 * t56 + t115 * t40 - t116 * t41, t105 * t83 + t116 * t72 - t136 * t40 - t436 * t55, t103 * t56 + t115 * t41, -t103 * t83 - t115 * t72 - t136 * t41 - t436 * t56, t436 * t83 + t400, t103 * t12 + t115 * t8 + t136 * t3 + t23 * t72 + t30 * t56 - t306 * t83 + t41 * t49 + t436 * t5, t105 * t12 + t116 * t8 - t136 * t2 - t16 * t83 - t24 * t72 - t30 * t55 - t4 * t436 - t40 * t49, -t103 * t4 - t105 * t5 - t115 * t2 - t116 * t3 - t16 * t56 + t23 * t40 - t24 * t41 - t306 * t55, t12 * t30 + t16 * t4 + t2 * t24 + t23 * t3 - t306 * t5 + t49 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t234, t209, t189, t234, -t188, 0, -pkin(8) * t246 + t377 + (-t277 * t348 + t372) * t409, pkin(1) * t341 - t204 * t259 + t365, 0, 0, 0, -t189, t188, -t234, t209, t234 ((-qJ(3) * qJD(2) - t177 - t205) * t281 + (-pkin(2) * qJD(2) - t172 + t350) * t284) * t362, -t377 + t154 + (-t203 * t370 + t207) * qJD(1), t350 * t259 + (t190 * t284 + t203 * t281) * t362 - t163, -pkin(2) * t197 - qJ(3) * t163 - t172 * t205 - t177 * t350 - t190 * t203, -t140 * t283 - t280 * t379 (-t141 - t379) * t283 + (t140 + t381) * t280, -t244 * t357 + t228 + (-t194 * t284 - t244 * t368) * t362, t141 * t280 + t283 * t381, -t334 + (-t342 + (-qJD(2) * t280 + t295) * t284) * t362, -t244 * t337, qJ(3) * t141 - t111 * t244 + t134 * t280 + t351 * t295 + (t143 * t283 - t280 * t373) * qJD(4) + (t283 * t299 - t284 * t90) * t362, -qJ(3) * t140 + t112 * t244 + t134 * t283 + t351 * t194 + (-t143 * t280 - t283 * t373) * qJD(4) + (-t280 * t299 + t284 * t91) * t362, t111 * t194 + t112 * t295 + (-t91 * t338 + t140 * t285 - t62 + (-t285 * t295 - t91) * qJD(4)) * t283 + (t90 * t338 - t141 * t285 - t61 + (t194 * t285 + t90) * qJD(4)) * t280, qJ(3) * t134 - t111 * t90 - t112 * t91 + t351 * t143 + (t280 * t61 + t283 * t62 + (-t280 * t90 + t283 * t91) * qJD(4)) * t285, t414, t119 * t367 + t219 * t71 - t220 * t72 - t301 * t366, -t301 * t337 + t289, t416 ((-qJD(2) * t219 + t119) * t284 - t366 * t281) * t362 - t366 * t347, -qJD(1) * t290, t265 * t72 + t95 * t219 + t364 * t119 + t366 * t113 + ((qJD(2) * t300 - t32) * t284 - t390 * t281) * t362 - t390 * t347, -t265 * t71 + t95 * t220 + t364 * t301 - t367 * t113 + ((-qJD(2) * t171 + t33) * t284 - t392 * t281) * t362 - t392 * t347, -t119 * t392 - t171 * t72 + t300 * t71 + t301 * t390 - t413, t10 * t300 + t113 * t364 - t171 * t327 + t265 * t95 - t32 * t390 + t33 * t392, t105 * t293 - t374 * t40, t323 * t105 + t322 * t103 + (t37 - t39 + (-t386 + t388) * qJD(6)) * t220, t105 * t366 - t219 * t40 + t293 * t436 + t374 * t72, t103 * t294 + t375 * t41, -t103 * t366 - t219 * t41 - t294 * t436 - t375 * t72, t366 * t436 + t397, t101 * t72 + t103 * t391 + t219 * t3 + t294 * t30 - t300 * t41 - t306 * t366 + t375 * t8 - t406 * t436, -t102 * t72 + t105 * t391 - t16 * t366 - t2 * t219 + t293 * t30 + t300 * t40 + t374 * t8 + t405 * t436, t101 * t40 - t102 * t41 + t323 * t16 - t322 * t306 + t406 * t105 + t405 * t103 + (-qJD(6) * t423 - t2 * t278 - t282 * t3) * t220, t101 * t3 + t102 * t2 - t16 * t405 + t30 * t391 - t300 * t8 + t306 * t406; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, t234, -t259 ^ 2 - t274 * t372, t177 * t259 + t154 + t197, 0, 0, 0, 0, 0, 0, -t244 ^ 2 * t280 - t259 * t295 + t228, -t334 - t194 * t259 + (-t280 * t359 - t342) * t362, t280 * t424 + t283 * t425, -t143 * t259 + t280 * t428 + t283 * t427, 0, 0, 0, 0, 0, 0, -t259 * t119 + t289, -t219 * t246 - t259 * t301 - t302 * t366, -t414 - t416, -t113 * t259 + t413, 0, 0, 0, 0, 0, 0, -t219 * t68 - t220 * t41 + t367 * t103 + (-t219 * t352 - t325) * t436, -t219 * t70 + t220 * t40 + t367 * t105 + (t219 * t353 - t324) * t436, t325 * t105 - t324 * t103 + (-t37 - t39 + (t386 + t388) * qJD(6)) * t219, t16 * t324 + t219 * t287 - t220 * t8 + t30 * t367 + t306 * t325; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t380, t194 ^ 2 - t295 ^ 2, -t425, -t380, t424, t246, -t143 * t194 + t427, t143 * t295 - t428, 0, 0, t382, t435, t434, -t382, t415, t246, t34 * t244 + (t34 - t33) * qJD(5) + (-t194 * t119 + t246 * t410 - t302 * t355) * pkin(4) + t417, t35 * t302 + (-t194 * t301 - t246 * t279 - t302 * t332) * pkin(4) + t433, t33 * t301 + t35 * t119 - t32 * t119 - t34 * t301 + (t410 * t71 - t279 * t72 + (-t119 * t410 + t279 * t301) * qJD(5)) * pkin(4), t32 * t34 - t33 * t35 + (t410 * t10 - t113 * t194 - t279 * t327 + (-t279 * t32 + t33 * t410) * qJD(5)) * pkin(4), t442, t446, t441, t388 * t436 - t39, t445, -t426, t268 * t41 + t307 * t278 + t313 * t103 + (-t267 * t352 - t278 * t318 - t17) * t436 + t420, -t268 * t40 + t307 * t282 + t313 * t105 + (t267 * t353 - t282 * t318 + t18) * t436 + t419, t18 * t103 + t17 * t105 + t1 + (-t103 * t318 + t119 * t306 - t267 * t41 + (t105 * t267 + t306) * qJD(6)) * t282 + (t105 * t318 - t119 * t16 - t267 * t40 - t3 + (t103 * t267 - t16) * qJD(6)) * t278, t306 * t17 - t16 * t18 + t8 * t268 - t30 * t34 + (t279 * t30 + t410 * t423) * qJD(5) * pkin(4) + t287 * t267; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t382, t435, t434, -t382, t415, t246, t244 * t33 + t417, t302 * t32 + t433, 0, 0, t442, t446, t441, t103 * t444 - t39, t445, -t426, -pkin(5) * t41 - pkin(11) * t389 - t103 * t33 - t20 * t436 + t278 * t440 + t420, t30 * t437 + pkin(5) * t40 - t105 * t33 + t436 * t21 + (t353 * t436 - t70) * pkin(11) + t419, t103 * t21 + t105 * t20 + t1 + (t402 + (-t41 + t354) * pkin(11)) * t282 + ((qJD(6) * t103 - t40) * pkin(11) + t429) * t278, -pkin(5) * t8 + pkin(11) * t287 - t16 * t21 + t20 * t306 - t30 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t387, -t103 ^ 2 + t105 ^ 2, t103 * t436 - t40, -t387, t105 * t436 - t41, t72, -t105 * t30 - t429, t103 * t30 - t2 - t402, 0, 0;];
tauc_reg  = t6;
