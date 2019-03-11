% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% tau_reg [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRRP8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP8_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:25:20
% EndTime: 2019-03-09 12:25:39
% DurationCPUTime: 8.34s
% Computational Cost: add. (12451->600), mult. (28052->751), div. (0->0), fcn. (20951->14), ass. (0->286)
t257 = sin(pkin(10));
t260 = sin(qJ(4));
t258 = cos(pkin(10));
t263 = cos(qJ(4));
t375 = t258 * t263;
t200 = t257 * t260 - t375;
t264 = cos(qJ(2));
t298 = t264 * t200;
t438 = -qJD(1) * t298 + t200 * qJD(4);
t261 = sin(qJ(2));
t357 = qJD(1) * t261;
t338 = t257 * t357;
t348 = t258 * qJD(2);
t192 = t338 - t348;
t337 = t258 * t357;
t355 = qJD(2) * t257;
t194 = t337 + t355;
t121 = t192 * t260 - t194 * t263;
t259 = sin(qJ(5));
t315 = -t192 * t263 - t260 * t194;
t404 = cos(qJ(5));
t76 = t259 * t121 + t315 * t404;
t440 = t76 ^ 2;
t356 = qJD(1) * t264;
t231 = -qJD(4) + t356;
t224 = -qJD(5) + t231;
t439 = t224 * t76;
t201 = t257 * t263 + t258 * t260;
t421 = t201 * t264;
t155 = qJD(1) * t421;
t291 = t201 * qJD(4);
t415 = t155 - t291;
t317 = pkin(2) * t261 - qJ(3) * t264;
t206 = t317 * qJD(1);
t140 = pkin(7) * t338 + t258 * t206;
t374 = t258 * t264;
t313 = pkin(3) * t261 - pkin(8) * t374;
t110 = qJD(1) * t313 + t140;
t183 = t257 * t206;
t376 = t258 * t261;
t377 = t257 * t264;
t300 = -pkin(7) * t376 - pkin(8) * t377;
t127 = qJD(1) * t300 + t183;
t392 = pkin(8) + qJ(3);
t215 = t392 * t257;
t216 = t392 * t258;
t361 = -t260 * t215 + t263 * t216;
t437 = t201 * qJD(3) + qJD(4) * t361 + t263 * t110 - t127 * t260;
t350 = qJD(4) * t263;
t436 = qJD(3) * t375 - t260 * t110 - t263 * t127 - t215 * t350;
t246 = t264 * qJDD(1);
t347 = qJD(1) * qJD(2);
t294 = t261 * t347 - t246;
t199 = qJDD(4) + t294;
t186 = qJDD(5) + t199;
t254 = pkin(10) + qJ(4);
t247 = qJ(5) + t254;
t233 = sin(t247);
t252 = g(3) * t264;
t262 = sin(qJ(1));
t265 = cos(qJ(1));
t320 = g(1) * t265 + g(2) * t262;
t306 = t320 * t261;
t281 = t306 - t252;
t372 = t260 * t216;
t325 = -t263 * t215 - t372;
t401 = pkin(9) * t201;
t104 = t325 - t401;
t105 = -pkin(9) * t200 + t361;
t67 = t259 * t104 + t105 * t404;
t435 = t67 * t186 + t233 * t281;
t426 = -t121 * t404 + t259 * t315;
t405 = t426 ^ 2;
t434 = pkin(4) * t357 - pkin(9) * t438 + t437;
t352 = qJD(3) * t257;
t433 = pkin(9) * t155 - t260 * t352 + (-t372 - t401) * qJD(4) + t436;
t432 = t224 * t426;
t422 = t200 * t261;
t431 = pkin(9) * t422;
t328 = -qJD(2) * pkin(2) + qJD(3);
t209 = pkin(7) * t357 + t328;
t139 = pkin(3) * t192 + t209;
t85 = -pkin(4) * t315 + t139;
t26 = -pkin(5) * t76 - qJ(6) * t426 + t85;
t430 = t26 * t76;
t429 = t76 * t85;
t335 = qJD(5) * t404;
t349 = qJD(5) * t259;
t387 = t200 * t335 + t201 * t349 - t259 * t415 + t404 * t438;
t126 = -t259 * t200 + t201 * t404;
t386 = qJD(5) * t126 - t259 * t438 - t404 * t415;
t393 = t426 * t76;
t428 = t121 * t231;
t427 = t405 - t440;
t242 = t258 * qJDD(2);
t333 = t264 * t347;
t345 = t261 * qJDD(1);
t295 = t333 + t345;
t419 = t257 * t295;
t147 = -t242 + t419;
t346 = t257 * qJDD(2);
t273 = t258 * t295 + t346;
t351 = qJD(4) * t260;
t312 = t147 * t263 - t192 * t351 + t194 * t350 + t260 * t273;
t69 = -t260 * t147 - t192 * t350 - t194 * t351 + t263 * t273;
t20 = -t121 * t349 + t259 * t312 - t315 * t335 - t404 * t69;
t13 = -t20 + t439;
t37 = pkin(5) * t426 - qJ(6) * t76;
t302 = t104 * t404 - t259 * t105;
t424 = qJD(5) * t302 - t259 * t434 + t404 * t433;
t423 = qJD(5) * t67 + t259 * t433 + t404 * t434;
t420 = t231 * t315;
t240 = pkin(7) * t356;
t343 = pkin(3) * t356;
t187 = t257 * t343 + t240;
t416 = -pkin(4) * t415 - t187;
t318 = pkin(2) * t264 + qJ(3) * t261;
t211 = -pkin(1) - t318;
t191 = t258 * t211;
t128 = -pkin(8) * t376 + t191 + (-pkin(7) * t257 - pkin(3)) * t264;
t153 = pkin(7) * t374 + t257 * t211;
t378 = t257 * t261;
t135 = -pkin(8) * t378 + t153;
t364 = t260 * t128 + t263 * t135;
t179 = t186 * qJ(6);
t214 = t224 * qJD(6);
t414 = t179 - t214;
t413 = pkin(7) * t333 + qJDD(3);
t181 = t186 * pkin(5);
t412 = t181 - qJDD(6);
t244 = sin(t254);
t245 = cos(t254);
t369 = t262 * t264;
t161 = t244 * t369 + t245 * t265;
t368 = t264 * t265;
t163 = -t244 * t368 + t245 * t262;
t394 = g(3) * t261;
t411 = -g(1) * t163 + g(2) * t161 + t244 * t394;
t234 = cos(t247);
t148 = t233 * t369 + t234 * t265;
t367 = t265 * t233;
t370 = t262 * t234;
t150 = t264 * t367 - t370;
t177 = qJD(2) * t317 - t261 * qJD(3);
t120 = qJD(1) * t177 + qJDD(1) * t211;
t160 = -pkin(7) * t294 + qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t81 = t258 * t120 - t257 * t160;
t56 = pkin(3) * t294 - pkin(8) * t273 + t81;
t82 = t257 * t120 + t258 * t160;
t70 = -pkin(8) * t147 + t82;
t185 = t211 * qJD(1);
t217 = qJD(2) * qJ(3) + t240;
t129 = t258 * t185 - t217 * t257;
t89 = -pkin(8) * t194 + t129 - t343;
t130 = t257 * t185 + t258 * t217;
t92 = -pkin(8) * t192 + t130;
t304 = t260 * t56 + t263 * t70 + t89 * t350 - t351 * t92;
t11 = -pkin(9) * t312 + t304;
t50 = -t260 * t92 + t263 * t89;
t41 = pkin(9) * t121 + t50;
t38 = -pkin(4) * t231 + t41;
t51 = t260 * t89 + t263 * t92;
t42 = pkin(9) * t315 + t51;
t330 = -t260 * t70 + t263 * t56;
t276 = -qJD(4) * t51 + t330;
t8 = t199 * pkin(4) - t69 * pkin(9) + t276;
t336 = t259 * t11 + t42 * t335 + t38 * t349 - t404 * t8;
t380 = t233 * t261;
t284 = g(1) * t150 + g(2) * t148 + g(3) * t380 - t336;
t271 = t26 * t426 - t284 - t412;
t410 = -t426 * t85 + t284;
t21 = qJD(5) * t426 + t259 * t69 + t404 * t312;
t409 = -t21 - t432;
t106 = -qJD(2) * t298 - t261 * t291;
t354 = qJD(2) * t261;
t342 = pkin(7) * t354;
t133 = t258 * t177 + t257 * t342;
t100 = qJD(2) * t313 + t133;
t158 = t257 * t177;
t112 = qJD(2) * t300 + t158;
t329 = t263 * t100 - t260 * t112;
t28 = pkin(4) * t354 - t106 * pkin(9) - qJD(4) * t364 + t329;
t373 = t260 * t135;
t326 = t263 * t128 - t373;
t59 = -pkin(4) * t264 + t326 + t431;
t170 = t201 * t261;
t68 = -pkin(9) * t170 + t364;
t307 = t259 * t59 + t404 * t68;
t286 = qJD(2) * t421;
t341 = t260 * t100 + t263 * t112 + t128 * t350;
t31 = -pkin(9) * t286 + (-t373 + t431) * qJD(4) + t341;
t407 = -qJD(5) * t307 - t259 * t31 + t28 * t404;
t406 = -0.2e1 * pkin(1);
t403 = pkin(3) * t257;
t402 = pkin(7) * t192;
t399 = g(1) * t262;
t395 = g(2) * t265;
t391 = -qJ(6) * t357 + t424;
t390 = pkin(5) * t357 + t423;
t389 = pkin(5) * t386 + qJ(6) * t387 - t126 * qJD(6) + t416;
t340 = t404 * t42;
t17 = t259 * t38 + t340;
t385 = t17 * t224;
t384 = t259 * t42;
t19 = t404 * t41 - t384;
t382 = pkin(4) * t335 + qJD(6) - t19;
t381 = qJDD(2) * pkin(2);
t379 = t234 * t261;
t371 = t261 * t265;
t16 = t38 * t404 - t384;
t366 = qJD(6) - t16;
t238 = pkin(7) * t345;
t360 = -t238 - t252;
t353 = qJD(2) * t264;
t188 = (pkin(7) + t403) * t353;
t207 = pkin(3) * t378 + t261 * pkin(7);
t359 = t265 * pkin(1) + t262 * pkin(7);
t255 = t261 ^ 2;
t358 = -t264 ^ 2 + t255;
t339 = pkin(3) * t258 + pkin(2);
t332 = t404 * t11 + t259 * t8 + t38 * t335 - t42 * t349;
t18 = t259 * t41 + t340;
t324 = pkin(4) * t349 - t18;
t232 = t261 * t399;
t323 = -g(2) * t371 + t232;
t132 = pkin(4) * t170 + t207;
t322 = -g(1) * t148 + g(2) * t150;
t149 = t234 * t369 - t367;
t151 = t233 * t262 + t234 * t368;
t321 = g(1) * t149 - g(2) * t151;
t319 = -t395 + t399;
t157 = pkin(4) * t200 - t339;
t205 = pkin(4) * t245 + t339;
t311 = pkin(5) * t234 + qJ(6) * t233 + t205;
t172 = t238 - t381 + t413;
t309 = -t259 * t68 + t404 * t59;
t305 = -pkin(7) * qJDD(2) + t347 * t406;
t303 = t259 * t28 + t404 * t31 + t59 * t335 - t349 * t68;
t99 = -t259 * t170 - t404 * t422;
t293 = t258 * t345 + t346;
t267 = qJD(1) ^ 2;
t290 = pkin(1) * t267 + t320;
t266 = qJD(2) ^ 2;
t289 = pkin(7) * t266 + qJDD(1) * t406 + t395;
t287 = g(2) * t261 * t370 + t302 * t186 + (g(1) * t371 - t252) * t234;
t285 = qJD(4) * t422;
t101 = t147 * pkin(3) + t172;
t283 = t20 + t439;
t282 = g(1) * t151 + g(2) * t149 + g(3) * t379 - t332;
t279 = -t264 * t320 - t394;
t278 = -t306 - t381;
t277 = -t172 + t281;
t274 = -t16 * t224 + t282;
t270 = t21 - t432;
t46 = pkin(4) * t312 + t101;
t269 = t285 - t286;
t268 = -g(1) * (-t150 * pkin(5) + t151 * qJ(6)) - g(2) * (-t148 * pkin(5) + t149 * qJ(6)) - g(3) * (-pkin(5) * t380 + qJ(6) * t379);
t86 = -pkin(4) * t269 + t188;
t3 = t21 * pkin(5) + t20 * qJ(6) - qJD(6) * t426 + t46;
t253 = -pkin(9) - t392;
t250 = t265 * pkin(7);
t237 = -pkin(4) * t404 - pkin(5);
t235 = pkin(4) * t259 + qJ(6);
t208 = pkin(4) * t244 + t403;
t164 = t244 * t262 + t245 * t368;
t162 = t244 * t265 - t245 * t369;
t152 = -pkin(7) * t377 + t191;
t141 = -pkin(7) * t337 + t183;
t134 = -t258 * t342 + t158;
t125 = t200 * t404 + t201 * t259;
t98 = t170 * t404 - t259 * t422;
t65 = pkin(5) * t125 - qJ(6) * t126 + t157;
t47 = pkin(5) * t98 - qJ(6) * t99 + t132;
t45 = qJD(5) * t99 + t259 * t106 - t269 * t404;
t44 = -t106 * t404 + t170 * t335 - t259 * t269 - t349 * t422;
t34 = -pkin(4) * t121 + t37;
t33 = t264 * pkin(5) - t309;
t32 = -qJ(6) * t264 + t307;
t15 = -t224 * qJ(6) + t17;
t14 = t224 * pkin(5) + t366;
t12 = t45 * pkin(5) + t44 * qJ(6) - t99 * qJD(6) + t86;
t5 = -pkin(5) * t354 - t407;
t4 = qJ(6) * t354 - qJD(6) * t264 + t303;
t2 = t336 - t412;
t1 = t332 + t414;
t6 = [qJDD(1), t319, t320, qJDD(1) * t255 + 0.2e1 * t261 * t333, 0.2e1 * t246 * t261 - 0.2e1 * t347 * t358, qJDD(2) * t261 + t264 * t266, qJDD(2) * t264 - t261 * t266, 0, t305 * t261 + (-t289 + t399) * t264, t261 * t289 + t264 * t305 - t232, -t320 * t257 + (pkin(7) * t147 + t172 * t257 + (qJD(1) * t152 + t129) * qJD(2)) * t261 + (-t133 * qJD(1) - t152 * qJDD(1) - t81 + t319 * t258 + (t209 * t257 + t402) * qJD(2)) * t264, -t320 * t258 + (t172 * t258 + (-qJD(1) * t153 - t130) * qJD(2) + t293 * pkin(7)) * t261 + (t134 * qJD(1) + t153 * qJDD(1) + t82 - t319 * t257 + (t209 * t258 + (t194 + t337) * pkin(7)) * qJD(2)) * t264, -t133 * t194 - t134 * t192 - t153 * t147 + (-qJDD(2) * t152 - t130 * t353 - t261 * t82) * t257 + (-t129 * t353 - t152 * t295 - t81 * t261) * t258 + t323, t82 * t153 + t130 * t134 + t81 * t152 + t129 * t133 - g(1) * t250 - g(2) * (t265 * t318 + t359) - t211 * t399 + (t172 * t261 + t209 * t353) * pkin(7), -t106 * t121 - t422 * t69, t106 * t315 - t121 * t269 - t69 * t170 + t312 * t422, -t106 * t231 - t121 * t354 - t199 * t422 - t264 * t69, -t170 * t199 + t312 * t264 - t231 * t285 + (t231 * t421 + t261 * t315) * qJD(2), -t199 * t264 - t231 * t354, -t329 * t231 + t326 * t199 - t330 * t264 - t188 * t315 + t207 * t312 + t101 * t170 - g(1) * t162 - g(2) * t164 + (t139 * t421 + t50 * t261) * qJD(2) + (-t139 * t422 + t231 * t364 + t264 * t51) * qJD(4) (-t135 * t351 + t341) * t231 - t364 * t199 + t304 * t264 - t51 * t354 - t188 * t121 + t207 * t69 - t101 * t422 + t139 * t106 - g(1) * t161 - g(2) * t163, -t20 * t99 - t426 * t44, t20 * t98 - t21 * t99 - t426 * t45 - t44 * t76, t186 * t99 + t20 * t264 + t224 * t44 + t354 * t426, -t186 * t98 + t21 * t264 + t224 * t45 + t354 * t76, -t186 * t264 - t224 * t354, t132 * t21 + t16 * t354 + t309 * t186 - t224 * t407 + t336 * t264 + t85 * t45 + t46 * t98 - t76 * t86 + t321, -t132 * t20 - t17 * t354 - t186 * t307 + t224 * t303 + t264 * t332 + t426 * t86 - t85 * t44 + t46 * t99 + t322, -t12 * t76 - t14 * t354 - t186 * t33 + t2 * t264 + t21 * t47 + t224 * t5 + t26 * t45 + t3 * t98 + t321, -t1 * t98 - t14 * t44 - t15 * t45 + t2 * t99 - t20 * t33 - t21 * t32 + t4 * t76 + t426 * t5 + t323, -t1 * t264 - t12 * t426 + t15 * t354 + t186 * t32 + t20 * t47 - t224 * t4 + t26 * t44 - t3 * t99 - t322, t1 * t32 + t15 * t4 + t3 * t47 + t26 * t12 + t2 * t33 + t14 * t5 - g(1) * (-t149 * pkin(5) - t148 * qJ(6) + t265 * t208 + t250) - g(2) * (t151 * pkin(5) + t150 * qJ(6) + t205 * t368 - t253 * t371 + t359) + (-g(1) * (-t205 * t264 + t261 * t253 - pkin(1)) - g(2) * t208) * t262; 0, 0, 0, -t261 * t267 * t264, t358 * t267, t345, t246, qJDD(2), t261 * t290 + t360, t394 + (-pkin(7) * qJDD(1) + t290) * t264, t257 * qJ(3) * t246 - pkin(2) * t147 + t277 * t258 + ((-qJ(3) * t355 - t129) * t261 + (-t402 + t140 + (qJD(3) - t209) * t257) * t264) * qJD(1), -t317 * t258 * qJDD(1) + (t172 + t278 + t252) * t257 + ((-qJ(3) * t348 + t130) * t261 + (-pkin(7) * t194 - t141 + (-t209 + t328) * t258) * t264) * qJD(1), t140 * t194 + t141 * t192 + (qJ(3) * t346 + qJD(3) * t194 + t130 * t356 - t81) * t257 + (t129 * t356 - qJD(3) * t192 + t82 + (-t147 + t419) * qJ(3)) * t258 + t279, -t209 * t240 - t129 * t140 - t130 * t141 + (-t129 * t257 + t130 * t258) * qJD(3) + t277 * pkin(2) + (-t81 * t257 + t82 * t258 + t279) * qJ(3), t121 * t438 + t69 * t201, -t121 * t415 - t69 * t200 - t201 * t312 - t315 * t438, t121 * t357 + t201 * t199 + t231 * t438, -t200 * t199 - t231 * t415 - t315 * t357, t231 * t357, t101 * t200 - t415 * t139 + t187 * t315 + t325 * t199 + t231 * t437 + t281 * t245 - t339 * t312 - t50 * t357, -t361 * t199 - t339 * t69 + t101 * t201 + t51 * t357 + t187 * t121 + ((-qJD(4) * t216 - t352) * t260 + t436) * t231 - t438 * t139 - t281 * t244, -t20 * t126 - t387 * t426, t20 * t125 - t126 * t21 - t386 * t426 - t387 * t76, t126 * t186 + t224 * t387 - t357 * t426, -t125 * t186 + t224 * t386 - t357 * t76, t224 * t357, t46 * t125 + t157 * t21 - t16 * t357 + t224 * t423 + t386 * t85 - t416 * t76 + t287, t46 * t126 - t157 * t20 + t17 * t357 + t224 * t424 - t387 * t85 + t416 * t426 - t435, t3 * t125 + t14 * t357 + t65 * t21 + t224 * t390 + t26 * t386 - t389 * t76 + t287, -t1 * t125 + t2 * t126 - t14 * t387 - t15 * t386 + t20 * t302 - t67 * t21 + t390 * t426 + t391 * t76 + t279, -t3 * t126 - t15 * t357 + t65 * t20 - t224 * t391 + t26 * t387 - t389 * t426 + t435, t1 * t67 - t2 * t302 + t3 * t65 + t389 * t26 + t391 * t15 + t390 * t14 + (-g(3) * t311 + t253 * t320) * t264 + (g(3) * t253 + t311 * t320) * t261; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t257 * t345 - t242 + (-t194 + t355) * t356 (t192 + t348) * t356 + t293, -t192 ^ 2 - t194 ^ 2, t129 * t194 + t130 * t192 + t278 - t360 + t413, 0, 0, 0, 0, 0, t312 + t428, t69 - t420, 0, 0, 0, 0, 0, t270, -t283, t270, -t405 - t440, t283, -t14 * t426 - t15 * t76 - t281 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121 * t315, t121 ^ 2 - t315 ^ 2, t69 + t420, -t312 + t428, t199, t121 * t139 - t51 * t231 + t276 + t411, g(1) * t164 - g(2) * t162 - t139 * t315 - t231 * t50 + t245 * t394 - t304, -t393, t427, t13, t409, t186, -t18 * t224 + (-t121 * t76 + t186 * t404 + t224 * t349) * pkin(4) + t410, -t19 * t224 - t429 + (t121 * t426 - t259 * t186 + t224 * t335) * pkin(4) + t282, -t237 * t186 + t224 * t324 + t34 * t76 - t271, -t237 * t20 - t235 * t21 + (t15 + t324) * t426 - (t14 - t382) * t76, t235 * t186 - t224 * t382 + t34 * t426 - t282 + t414 + t430, t1 * t235 + t2 * t237 - t26 * t34 - t14 * t18 + t382 * t15 + (t14 * t349 + t411) * pkin(4) + t268; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t393, t427, t13, t409, t186, -t385 + t410, t274 - t429, t37 * t76 + t181 - t271 - t385, pkin(5) * t20 - qJ(6) * t21 + (t15 - t17) * t426 - (t14 - t366) * t76, t37 * t426 + 0.2e1 * t179 - 0.2e1 * t214 - t274 + t430, -t2 * pkin(5) + t1 * qJ(6) - t14 * t17 + t15 * t366 - t26 * t37 + t268; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t186 - t393, t13, -t224 ^ 2 - t405, t15 * t224 + t271;];
tau_reg  = t6;
