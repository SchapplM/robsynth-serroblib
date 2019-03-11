% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
% 
% Output:
% tau_reg [6x30]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRPPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:23:09
% EndTime: 2019-03-09 15:23:25
% DurationCPUTime: 7.05s
% Computational Cost: add. (11790->507), mult. (28440->671), div. (0->0), fcn. (21729->18), ass. (0->290)
t259 = qJDD(2) + qJDD(3);
t267 = sin(pkin(10));
t379 = cos(pkin(10));
t262 = qJD(2) + qJD(3);
t270 = sin(qJ(3));
t274 = cos(qJ(3));
t271 = sin(qJ(2));
t356 = qJD(1) * t271;
t338 = t270 * t356;
t275 = cos(qJ(2));
t346 = t275 * qJDD(1);
t349 = qJD(1) * qJD(2);
t335 = t275 * t349;
t347 = t271 * qJDD(1);
t299 = -t335 - t347;
t348 = qJD(1) * qJD(3);
t406 = -t275 * t348 + t299;
t130 = -t262 * t338 + t270 * t346 - t274 * t406;
t355 = qJD(1) * t275;
t182 = -t270 * t355 - t274 * t356;
t400 = pkin(7) + pkin(8);
t155 = qJDD(2) * pkin(2) + t400 * t299;
t336 = t271 * t349;
t298 = -t336 + t346;
t159 = t400 * t298;
t219 = t400 * t275;
t209 = qJD(1) * t219;
t191 = t274 * t209;
t218 = t400 * t271;
t207 = qJD(1) * t218;
t387 = qJD(2) * pkin(2);
t194 = -t207 + t387;
t307 = -t194 * t270 - t191;
t289 = qJD(3) * t307 + t274 * t155 - t270 * t159;
t61 = pkin(3) * t259 - qJ(4) * t130 + qJD(4) * t182 + t289;
t354 = qJD(3) * t270;
t177 = t209 * t354;
t181 = t274 * t355 - t338;
t324 = -qJD(3) * t194 - t159;
t65 = t181 * qJD(4) - t177 + (qJ(4) * t406 + t155) * t270 + ((-t271 * t348 + t298) * qJ(4) - t324) * t274;
t26 = -t267 * t65 + t379 * t61;
t25 = -t259 * pkin(4) + qJDD(5) - t26;
t265 = qJ(2) + qJ(3);
t253 = pkin(10) + t265;
t241 = sin(t253);
t272 = sin(qJ(1));
t276 = cos(qJ(1));
t319 = g(1) * t276 + g(2) * t272;
t415 = t319 * t241;
t417 = t25 - t415;
t242 = cos(t253);
t416 = -g(3) * t242 - t417;
t147 = t379 * t181 + t182 * t267;
t142 = qJD(6) - t147;
t200 = t270 * t271 - t274 * t275;
t257 = t275 * pkin(2);
t388 = pkin(1) + t257;
t323 = pkin(3) * t200 - t388;
t254 = sin(t265);
t255 = cos(t265);
t414 = -g(3) * t255 + t254 * t319;
t217 = qJD(1) * t388;
t413 = qJD(6) - t142;
t412 = qJDD(1) * t388;
t266 = sin(pkin(11));
t268 = cos(pkin(11));
t300 = t267 * t181 - t379 * t182;
t131 = -t268 * t262 + t266 * t300;
t133 = t262 * t266 + t268 * t300;
t269 = sin(qJ(6));
t273 = cos(qJ(6));
t91 = t273 * t131 + t133 * t269;
t411 = t142 * t91;
t308 = t131 * t269 - t133 * t273;
t410 = t142 * t308;
t201 = t270 * t275 + t271 * t274;
t154 = t262 * t201;
t197 = t266 * t273 + t268 * t269;
t180 = t197 * qJD(6);
t408 = t197 * t147 - t180;
t196 = t266 * t269 - t273 * t268;
t407 = t142 * t196;
t326 = t207 * t270 - t191;
t378 = qJ(4) * t181;
t135 = t326 - t378;
t175 = t182 * qJ(4);
t187 = t270 * t209;
t360 = -t274 * t207 - t187;
t136 = t175 + t360;
t328 = t379 * t270;
t381 = t379 * t135 - t136 * t267 + (t267 * t274 + t328) * qJD(3) * pkin(2);
t399 = pkin(3) * t182;
t103 = pkin(4) * t300 - qJ(5) * t147 - t399;
t129 = -t307 + t378;
t122 = t267 * t129;
t327 = t274 * t194 - t187;
t128 = t175 + t327;
t85 = t379 * t128 - t122;
t50 = t266 * t103 + t268 * t85;
t404 = -qJD(5) * t268 + t50;
t235 = t267 * t270 * pkin(2);
t353 = qJD(3) * t274;
t172 = t379 * pkin(2) * t353 - qJD(3) * t235;
t165 = qJD(5) + t172;
t249 = pkin(2) * t356;
t100 = t103 + t249;
t99 = t267 * t135 + t379 * t136;
t52 = t266 * t100 + t268 * t99;
t403 = -t165 * t268 + t52;
t359 = -t270 * t218 + t274 * t219;
t314 = t242 * pkin(4) + t241 * qJ(5);
t402 = g(1) * t272 - g(2) * t276;
t282 = qJD(1) * t154;
t280 = -t200 * qJDD(1) - t282;
t86 = t130 * t267 - t379 * t280;
t83 = qJDD(6) + t86;
t401 = -t142 * t407 + t197 * t83;
t87 = t379 * t130 + t280 * t267;
t76 = -t268 * t259 + t266 * t87;
t77 = t259 * t266 + t268 * t87;
t24 = -qJD(6) * t308 + t269 * t77 + t273 * t76;
t397 = pkin(3) * t254;
t396 = pkin(4) * t241;
t261 = pkin(11) + qJ(6);
t251 = sin(t261);
t392 = g(3) * t251;
t252 = cos(t261);
t391 = g(3) * t252;
t389 = t268 * pkin(5);
t256 = t268 * pkin(9);
t27 = t267 * t61 + t379 * t65;
t21 = qJ(5) * t259 + qJD(5) * t262 + t27;
t238 = pkin(2) * t336;
t279 = pkin(3) * t282 + t323 * qJDD(1) + qJDD(4) + t238;
t30 = t86 * pkin(4) - t87 * qJ(5) - qJD(5) * t300 + t279;
t7 = t268 * t21 + t266 * t30;
t5 = t7 * t268;
t339 = qJD(2) * t400;
t208 = t271 * t339;
t210 = t275 * t339;
t297 = -t274 * t208 - t270 * t210 - t218 * t353 - t219 * t354;
t96 = -qJ(4) * t154 - qJD(4) * t200 + t297;
t153 = t262 * t200;
t288 = -t359 * qJD(3) + t208 * t270 - t274 * t210;
t97 = qJ(4) * t153 - qJD(4) * t201 + t288;
t48 = t267 * t97 + t379 * t96;
t116 = -t153 * t267 + t379 * t154;
t117 = -t379 * t153 - t267 * t154;
t151 = -t267 * t200 + t379 * t201;
t250 = t271 * t387;
t333 = pkin(3) * t154 + t250;
t55 = pkin(4) * t116 - qJ(5) * t117 - qJD(5) * t151 + t333;
t19 = t266 * t55 + t268 * t48;
t120 = pkin(3) * t262 + t128;
t329 = t379 * t129;
t79 = t267 * t120 + t329;
t75 = qJ(5) * t262 + t79;
t156 = -pkin(3) * t181 + qJD(4) - t217;
t89 = -pkin(4) * t147 - qJ(5) * t300 + t156;
t42 = t266 * t89 + t268 * t75;
t386 = t300 * t91;
t385 = t300 * t308;
t78 = t379 * t120 - t122;
t74 = -t262 * pkin(4) + qJD(5) - t78;
t384 = t147 * t74;
t374 = t147 * t266;
t139 = pkin(5) * t374;
t382 = -t139 + t381;
t380 = t172 - t99;
t377 = qJ(5) * t242;
t376 = t117 * t266;
t375 = t142 * t300;
t373 = t147 * t268;
t372 = t151 * t266;
t371 = t151 * t268;
t247 = pkin(2) * t274 + pkin(3);
t174 = pkin(2) * t328 + t267 * t247;
t168 = qJ(5) + t174;
t369 = t168 * t268;
t368 = t182 * t181;
t239 = pkin(3) * t267 + qJ(5);
t367 = t239 * t268;
t366 = t242 * t266;
t365 = t251 * t272;
t364 = t251 * t276;
t363 = t252 * t272;
t362 = t252 * t276;
t150 = t379 * t200 + t201 * t267;
t109 = pkin(4) * t150 - qJ(5) * t151 + t323;
t325 = -t274 * t218 - t219 * t270;
t143 = -qJ(4) * t201 + t325;
t144 = -qJ(4) * t200 + t359;
t111 = t267 * t143 + t379 * t144;
t59 = t266 * t109 + t268 * t111;
t245 = pkin(3) * t255;
t358 = t245 + t257;
t263 = t271 ^ 2;
t357 = -t275 ^ 2 + t263;
t351 = qJD(6) * t269;
t350 = qJD(6) * t273;
t345 = pkin(9) * t374;
t342 = t245 + t314;
t6 = -t21 * t266 + t268 * t30;
t3 = pkin(5) * t86 - pkin(9) * t77 + t6;
t4 = -pkin(9) * t76 + t7;
t340 = -t269 * t4 + t273 * t3;
t211 = -pkin(2) * t271 - t397;
t332 = t211 - t396;
t18 = -t266 * t48 + t268 * t55;
t41 = -t266 * t75 + t268 * t89;
t51 = t268 * t100 - t266 * t99;
t47 = t267 * t96 - t379 * t97;
t49 = t268 * t103 - t266 * t85;
t58 = t268 * t109 - t111 * t266;
t84 = t128 * t267 + t329;
t110 = -t379 * t143 + t144 * t267;
t322 = pkin(5) * t300 - pkin(9) * t373;
t243 = -t379 * pkin(3) - pkin(4);
t321 = t142 * t408 - t196 * t83;
t320 = -t396 - t397;
t317 = -t266 * t7 - t268 * t6;
t316 = -t6 * t266 + t5;
t315 = t269 * t3 + t273 * t4;
t313 = t78 * t147 + t300 * t79;
t312 = t266 * t41 - t268 * t42;
t31 = -pkin(5) * t147 - pkin(9) * t133 + t41;
t35 = -pkin(9) * t131 + t42;
t311 = t269 * t35 - t273 * t31;
t10 = t269 * t31 + t273 * t35;
t40 = pkin(5) * t150 - pkin(9) * t371 + t58;
t45 = -pkin(9) * t372 + t59;
t310 = -t269 * t45 + t273 * t40;
t309 = t269 * t40 + t273 * t45;
t173 = t379 * t247 - t235;
t169 = -pkin(4) - t173;
t305 = -0.2e1 * pkin(1) * t349 - pkin(7) * qJDD(2);
t157 = (-pkin(9) - t168) * t266;
t304 = -qJD(6) * t157 - t345 + t403;
t158 = t256 + t369;
t303 = qJD(6) * t158 + t165 * t266 + t322 + t51;
t183 = (-pkin(9) - t239) * t266;
t302 = -qJD(6) * t183 - t345 + t404;
t184 = t256 + t367;
t301 = qJD(5) * t266 + qJD(6) * t184 + t322 + t49;
t23 = -t131 * t350 - t133 * t351 - t269 * t76 + t273 * t77;
t295 = t147 * t165 - t168 * t86 - t384;
t294 = qJD(5) * t147 - t239 * t86 - t384;
t277 = qJD(2) ^ 2;
t293 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t277 + t402;
t278 = qJD(1) ^ 2;
t292 = pkin(1) * t278 - pkin(7) * qJDD(1) + t319;
t291 = t117 * t74 + t151 * t25 - t319;
t290 = t268 * t416 - t300 * t41;
t287 = -g(3) * t241 - t242 * t319 + t41 * t373 + t42 * t374 + t5;
t15 = t76 * pkin(5) + t25;
t62 = t131 * pkin(5) + t74;
t286 = t15 * t196 - t242 * t391 + t252 * t415 + t300 * t311 - t408 * t62;
t285 = g(3) * t366 + t417 * t266 + t300 * t42;
t284 = g(3) * t254 - t270 * t155 + t217 * t181 + t255 * t319 + t274 * t324 + t177;
t283 = t10 * t300 + t15 * t197 + t242 * t392 - t251 * t415 - t407 * t62;
t281 = -t217 * t182 + t289 + t414;
t260 = -qJ(4) - t400;
t214 = t276 * t377;
t213 = t272 * t377;
t212 = t243 - t389;
t206 = pkin(1) + t358;
t193 = t276 * t206;
t176 = t238 - t412;
t164 = t169 - t389;
t163 = t242 * t362 + t365;
t162 = -t242 * t364 + t363;
t161 = -t242 * t363 + t364;
t160 = t242 * t365 + t362;
t137 = -t181 ^ 2 + t182 ^ 2;
t115 = -t182 * t262 + t280;
t114 = -t181 * t262 + t130;
t113 = t196 * t151;
t112 = t197 * t151;
t68 = pkin(5) * t372 + t110;
t66 = t139 + t84;
t44 = t117 * t197 + t350 * t371 - t351 * t372;
t43 = -t117 * t196 - t151 * t180;
t34 = pkin(5) * t376 + t47;
t17 = t321 + t386;
t16 = t385 + t401;
t12 = -pkin(9) * t376 + t19;
t11 = pkin(5) * t116 - t117 * t256 + t18;
t8 = t197 * t23 + t308 * t407;
t1 = -t196 * t23 - t197 * t24 - t308 * t408 + t407 * t91;
t2 = [qJDD(1), t402, t319, qJDD(1) * t263 + 0.2e1 * t271 * t335, 0.2e1 * t271 * t346 - 0.2e1 * t357 * t349, qJDD(2) * t271 + t275 * t277, qJDD(2) * t275 - t271 * t277, 0, t271 * t305 + t275 * t293, -t271 * t293 + t275 * t305, t130 * t201 + t153 * t182, -t130 * t200 - t153 * t181 + t182 * t154 + t280 * t201, -t153 * t262 + t201 * t259, -t154 * t262 - t200 * t259, 0, -t181 * t250 + t255 * t402 + t259 * t325 + t262 * t288 + (t176 - t412) * t200 - 0.2e1 * t217 * t154, -t130 * t388 + t217 * t153 + t176 * t201 - t182 * t250 - t254 * t402 - t259 * t359 - t262 * t297, t110 * t87 - t111 * t86 - t116 * t79 - t117 * t78 + t147 * t48 - t150 * t27 - t151 * t26 + t300 * t47 - t319, t27 * t111 + t79 * t48 - t26 * t110 - t78 * t47 + t279 * t323 + t156 * t333 - g(1) * (-t206 * t272 - t260 * t276) - g(2) * (-t260 * t272 + t193) t242 * t268 * t402 + t110 * t76 + t41 * t116 + t47 * t131 - t147 * t18 + t6 * t150 + t266 * t291 + t58 * t86, t110 * t77 - t42 * t116 + t47 * t133 + t147 * t19 - t7 * t150 + t268 * t291 - t366 * t402 - t59 * t86, -t131 * t19 - t133 * t18 - t58 * t77 - t59 * t76 + t402 * t241 + t317 * t151 + (-t266 * t42 - t268 * t41) * t117, -g(2) * t193 + t25 * t110 + t41 * t18 + t42 * t19 + t74 * t47 + t6 * t58 + t7 * t59 + (g(1) * t260 - g(2) * t314) * t276 + (-g(1) * (-t206 - t314) + g(2) * t260) * t272, -t113 * t23 - t308 * t43, -t112 * t23 + t113 * t24 + t308 * t44 - t43 * t91, -t113 * t83 - t116 * t308 + t142 * t43 + t150 * t23, -t112 * t83 - t116 * t91 - t142 * t44 - t150 * t24, t116 * t142 + t150 * t83 (t273 * t11 - t269 * t12) * t142 + t310 * t83 + t340 * t150 - t311 * t116 + t34 * t91 + t68 * t24 + t15 * t112 + t62 * t44 - g(1) * t161 - g(2) * t163 + (-t10 * t150 - t142 * t309) * qJD(6) -(t11 * t269 + t12 * t273) * t142 - t309 * t83 - t315 * t150 - t10 * t116 - t34 * t308 + t68 * t23 - t15 * t113 + t62 * t43 - g(1) * t160 - g(2) * t162 + (-t142 * t310 + t150 * t311) * qJD(6); 0, 0, 0, -t271 * t278 * t275, t357 * t278, t347, t346, qJDD(2), -g(3) * t275 + t271 * t292, g(3) * t271 + t275 * t292, t368, t137, t114, t115, t259, -t326 * t262 + (t181 * t356 + t259 * t274 - t262 * t354) * pkin(2) + t281, t360 * t262 + (t182 * t356 - t259 * t270 - t262 * t353) * pkin(2) + t284, t147 * t380 - t173 * t87 - t174 * t86 + t300 * t381 + t313, t27 * t174 + t26 * t173 - t156 * (t249 - t399) - g(3) * t358 + t380 * t79 - t381 * t78 - t319 * t211, t131 * t381 + t147 * t51 + t169 * t76 + t266 * t295 + t290, t133 * t381 - t147 * t52 + t169 * t77 + t268 * t295 + t285, -t76 * t369 + t133 * t51 + t403 * t131 + (t133 * t165 + t168 * t77 - t6) * t266 + t287, t25 * t169 - t42 * t52 - t41 * t51 - g(1) * (t276 * t332 + t214) - g(2) * (t272 * t332 + t213) - g(3) * (t257 + t342) + t381 * t74 + t316 * t168 - t312 * t165, t8, t1, t16, t17, -t375 (t157 * t273 - t158 * t269) * t83 + t164 * t24 + t382 * t91 + (t269 * t304 - t273 * t303) * t142 + t286 -(t157 * t269 + t158 * t273) * t83 + t164 * t23 - t382 * t308 + (t269 * t303 + t273 * t304) * t142 + t283; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t368, t137, t114, t115, t259, -t262 * t307 + t281, t262 * t327 + t284, -t85 * t147 - t84 * t300 + (-t267 * t86 - t379 * t87) * pkin(3) + t313, t78 * t84 - t79 * t85 + (t156 * t182 + t26 * t379 + t267 * t27 + t414) * pkin(3), -t131 * t84 + t147 * t49 + t243 * t76 + t266 * t294 + t290, -t133 * t84 - t147 * t50 + t243 * t77 + t268 * t294 + t285, -t76 * t367 + t133 * t49 + t404 * t131 + (qJD(5) * t133 + t239 * t77 - t6) * t266 + t287, t25 * t243 - t42 * t50 - t41 * t49 - t74 * t84 - g(1) * (t276 * t320 + t214) - g(2) * (t272 * t320 + t213) - g(3) * t342 + t316 * t239 - t312 * qJD(5), t8, t1, t16, t17, -t375 (t183 * t273 - t184 * t269) * t83 + t212 * t24 - t66 * t91 + (t269 * t302 - t273 * t301) * t142 + t286 -(t183 * t269 + t184 * t273) * t83 + t212 * t23 + t66 * t308 + (t269 * t301 + t273 * t302) * t142 + t283; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t147 ^ 2 - t300 ^ 2, -t79 * t147 + t300 * t78 + t279 - t402, -t131 * t300 - t147 * t374 + t268 * t86, -t133 * t300 - t147 * t373 - t266 * t86, -t266 * t76 - t268 * t77 + (t131 * t268 - t133 * t266) * t147, t147 * t312 - t300 * t74 - t317 - t402, 0, 0, 0, 0, 0, t321 - t386, t385 - t401; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133 * t147 + t76, t131 * t147 + t77, -t131 ^ 2 - t133 ^ 2, t131 * t42 + t133 * t41 - t416, 0, 0, 0, 0, 0, t24 - t410, t23 - t411; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t308 * t91, t308 ^ 2 - t91 ^ 2, t23 + t411, -t24 - t410, t83, -g(1) * t162 + g(2) * t160 - t413 * t10 + t241 * t392 + t62 * t308 + t340, g(1) * t163 - g(2) * t161 + t241 * t391 + t413 * t311 + t62 * t91 - t315;];
tau_reg  = t2;
