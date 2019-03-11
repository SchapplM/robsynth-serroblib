% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% tau_reg [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRPR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:10:20
% EndTime: 2019-03-09 05:10:36
% DurationCPUTime: 7.10s
% Computational Cost: add. (10661->476), mult. (26127->619), div. (0->0), fcn. (21246->18), ass. (0->262)
t239 = cos(pkin(10));
t245 = cos(qJ(3));
t323 = t245 * t239;
t237 = sin(pkin(10));
t242 = sin(qJ(3));
t324 = t237 * t242;
t176 = -t323 + t324;
t161 = t176 * qJD(1);
t178 = t237 * t245 + t239 * t242;
t162 = t178 * qJD(1);
t241 = sin(qJ(4));
t359 = cos(qJ(4));
t136 = t359 * t161 + t162 * t241;
t134 = qJD(6) + t136;
t236 = sin(pkin(11));
t238 = cos(pkin(11));
t240 = sin(qJ(6));
t244 = cos(qJ(6));
t177 = t236 * t244 + t238 * t240;
t175 = t236 * t240 - t244 * t238;
t400 = t134 * t175;
t313 = qJD(1) * qJD(3);
t302 = t242 * t313;
t301 = t245 * t313;
t311 = t239 * qJDD(1);
t312 = t237 * qJDD(1);
t307 = t239 * t301 + t242 * t311 + t245 * t312;
t143 = -t237 * t302 + t307;
t308 = -t237 * t301 - t239 * t302 - t242 * t312;
t263 = t245 * t311 + t308;
t267 = -t241 * t161 + t359 * t162;
t77 = qJD(4) * t267 + t241 * t143 - t359 * t263;
t73 = qJDD(6) + t77;
t402 = -t400 * t134 + t177 * t73;
t164 = t177 * qJD(6);
t401 = t177 * t136 + t164;
t386 = t136 * t236;
t399 = pkin(5) * t386;
t398 = pkin(9) * t386;
t288 = -t401 * t134 - t175 * t73;
t235 = qJD(3) + qJD(4);
t124 = -t238 * t235 + t236 * t267;
t126 = t235 * t236 + t238 * t267;
t86 = t244 * t124 + t126 * t240;
t349 = t267 * t86;
t397 = t288 + t349;
t396 = qJD(6) - t134;
t351 = pkin(7) + qJ(2);
t196 = t351 * t237;
t181 = qJD(1) * t196;
t198 = t351 * t239;
t182 = qJD(1) * t198;
t369 = -t245 * t181 - t182 * t242;
t128 = -pkin(8) * t162 + t369;
t123 = qJD(3) * pkin(3) + t128;
t275 = t181 * t242 - t182 * t245;
t129 = -pkin(8) * t161 - t275;
t303 = qJD(4) * t359;
t317 = qJD(4) * t241;
t314 = qJD(1) * qJD(2);
t360 = t351 * qJDD(1) + t314;
t154 = t360 * t237;
t155 = t360 * t239;
t295 = -t245 * t154 - t242 * t155;
t71 = qJDD(3) * pkin(3) - pkin(8) * t143 + qJD(3) * t275 + t295;
t276 = -t242 * t154 + t245 * t155;
t75 = t263 * pkin(8) + t369 * qJD(3) + t276;
t292 = -t123 * t317 - t129 * t303 - t241 * t75 + t359 * t71;
t229 = qJDD(3) + qJDD(4);
t366 = -pkin(4) * t229 + qJDD(5);
t24 = -t292 + t366;
t76 = t359 * t143 - t161 * t303 - t162 * t317 + t241 * t263;
t64 = -t238 * t229 + t236 * t76;
t13 = pkin(5) * t64 + t24;
t234 = pkin(10) + qJ(3);
t226 = qJ(4) + t234;
t216 = cos(t226);
t233 = pkin(11) + qJ(6);
t224 = cos(t233);
t122 = t359 * t129;
t84 = t241 * t123 + t122;
t78 = qJ(5) * t235 + t84;
t306 = pkin(2) * t239 + pkin(1);
t186 = -qJD(1) * t306 + qJD(2);
t146 = pkin(3) * t161 + t186;
t85 = pkin(4) * t136 - qJ(5) * t267 + t146;
t38 = -t236 * t78 + t238 * t85;
t27 = pkin(5) * t136 - pkin(9) * t126 + t38;
t39 = t236 * t85 + t238 * t78;
t31 = -pkin(9) * t124 + t39;
t280 = t240 * t31 - t244 * t27;
t353 = g(3) * t224;
t215 = sin(t226);
t246 = cos(qJ(1));
t332 = t215 * t246;
t243 = sin(qJ(1));
t333 = t215 * t243;
t382 = g(1) * t332 + g(2) * t333;
t121 = t241 * t129;
t83 = t359 * t123 - t121;
t74 = -t235 * pkin(4) + qJD(5) - t83;
t55 = t124 * pkin(5) + t74;
t395 = t13 * t175 - t216 * t353 + t382 * t224 + t280 * t267 + t401 * t55;
t222 = sin(t233);
t286 = g(1) * t246 + g(2) * t243;
t272 = t215 * t286;
t354 = g(3) * t222;
t8 = t240 * t27 + t244 * t31;
t394 = t13 * t177 + t216 * t354 - t222 * t272 + t8 * t267 - t400 * t55;
t277 = t124 * t240 - t126 * t244;
t348 = t267 * t277;
t393 = t348 + t402;
t315 = qJD(6) * t244;
t316 = qJD(6) * t240;
t65 = t229 * t236 + t238 * t76;
t19 = -t124 * t315 - t126 * t316 - t240 * t64 + t244 * t65;
t392 = t19 * t177 + t400 * t277;
t20 = -qJD(6) * t277 + t240 * t65 + t244 * t64;
t391 = -t19 * t175 - t177 * t20 + t401 * t277 + t400 * t86;
t390 = t134 * t86;
t389 = t136 * t38;
t339 = t136 * t235;
t388 = t76 + t339;
t387 = t134 * t277;
t383 = t267 * t136;
t225 = cos(t234);
t320 = t216 * pkin(4) + t215 * qJ(5);
t381 = pkin(3) * t225 + t320;
t343 = qJDD(1) * pkin(1);
t368 = g(1) * t243 - g(2) * t246;
t274 = -qJDD(2) + t343 + t368;
t341 = t267 * t235;
t380 = -t77 + t341;
t378 = -t136 ^ 2 + t267 ^ 2;
t377 = -g(3) * t216 + t382;
t252 = t123 * t303 - t129 * t317 + t241 * t71 + t359 * t75;
t22 = t229 * qJ(5) + t235 * qJD(5) + t252;
t130 = qJDD(2) - t308 * pkin(3) + (-pkin(1) + (-pkin(3) * t245 - pkin(2)) * t239) * qJDD(1);
t28 = t77 * pkin(4) - t76 * qJ(5) - qJD(5) * t267 + t130;
t5 = -t22 * t236 + t238 * t28;
t376 = -t136 * t39 - t5;
t106 = pkin(4) * t267 + qJ(5) * t136;
t329 = t216 * t246;
t330 = t216 * t243;
t309 = -g(1) * t329 - g(2) * t330 - g(3) * t215;
t375 = t146 * t136 - t252 - t309;
t372 = pkin(5) * t267;
t89 = t241 * t128 + t122;
t289 = pkin(3) * t317 - t89;
t371 = t134 * t267;
t370 = t368 * t215;
t321 = -t242 * t196 + t245 * t198;
t367 = qJ(2) * qJDD(1);
t331 = t216 * t236;
t363 = t267 * t39 + g(3) * t331 + (t24 - t272) * t236;
t362 = -t267 * t38 + (-t24 + t377) * t238;
t259 = t292 + t377;
t361 = -t146 * t267 + t259;
t166 = t178 * qJD(3);
t358 = pkin(3) * t166;
t352 = t238 * pkin(5);
t227 = t238 * pkin(9);
t6 = t238 * t22 + t236 * t28;
t4 = t6 * t238;
t173 = t245 * t196;
t256 = -qJD(3) * t173 + qJD(2) * t323 + (-qJD(2) * t237 - qJD(3) * t198) * t242;
t111 = -pkin(8) * t166 + t256;
t165 = t176 * qJD(3);
t250 = -t178 * qJD(2) - t321 * qJD(3);
t112 = pkin(8) * t165 + t250;
t294 = -t198 * t242 - t173;
t131 = -pkin(8) * t178 + t294;
t132 = -pkin(8) * t176 + t321;
t268 = t359 * t131 - t241 * t132;
t42 = t268 * qJD(4) + t359 * t111 + t241 * t112;
t266 = -t359 * t176 - t241 * t178;
t109 = t266 * qJD(4) - t359 * t165 - t241 * t166;
t145 = -t241 * t176 + t359 * t178;
t110 = t145 * qJD(4) - t241 * t165 + t359 * t166;
t48 = pkin(4) * t110 - qJ(5) * t109 - qJD(5) * t145 + t358;
t15 = t236 * t48 + t238 * t42;
t90 = t359 * t128 - t121;
t93 = pkin(3) * t162 + t106;
t45 = t236 * t93 + t238 * t90;
t103 = t241 * t131 + t359 * t132;
t152 = pkin(3) * t176 - t306;
t97 = -pkin(4) * t266 - qJ(5) * t145 + t152;
t52 = t238 * t103 + t236 * t97;
t347 = t236 * t77;
t50 = t236 * t106 + t238 * t83;
t344 = qJ(5) * t238;
t342 = t109 * t236;
t337 = t136 * t238;
t336 = t145 * t236;
t335 = t145 * t238;
t328 = t222 * t243;
t327 = t222 * t246;
t326 = t224 * t243;
t325 = t224 * t246;
t322 = -qJD(5) + t74;
t319 = t237 ^ 2 + t239 ^ 2;
t318 = qJD(3) * t162;
t2 = pkin(5) * t77 - pkin(9) * t65 + t5;
t3 = -pkin(9) * t64 + t6;
t305 = t244 * t2 - t240 * t3;
t304 = qJD(1) * t324;
t14 = -t236 * t42 + t238 * t48;
t44 = -t236 * t90 + t238 * t93;
t51 = -t103 * t236 + t238 * t97;
t49 = t238 * t106 - t236 * t83;
t297 = t319 * qJD(1) ^ 2;
t293 = t4 + t309;
t291 = 0.2e1 * t319;
t220 = -t359 * pkin(3) - pkin(4);
t290 = t289 + t399;
t223 = sin(t234);
t287 = -pkin(3) * t223 - pkin(4) * t215;
t284 = t2 * t240 + t244 * t3;
t283 = -t236 * t6 - t238 * t5;
t282 = -t5 * t236 + t4;
t281 = t38 * t236 - t39 * t238;
t33 = -pkin(5) * t266 - pkin(9) * t335 + t51;
t41 = -pkin(9) * t336 + t52;
t279 = -t240 * t41 + t244 * t33;
t278 = t240 * t33 + t244 * t41;
t273 = t306 + t381;
t271 = t368 * t216;
t217 = pkin(3) * t241 + qJ(5);
t169 = (-pkin(9) - t217) * t236;
t206 = pkin(3) * t303 + qJD(5);
t270 = -qJD(6) * t169 - t238 * t206 + t398 + t45;
t170 = t217 * t238 + t227;
t269 = pkin(9) * t337 + qJD(6) * t170 + t236 * t206 + t372 + t44;
t195 = (-pkin(9) - qJ(5)) * t236;
t265 = -qJD(5) * t238 - qJD(6) * t195 + t398 + t50;
t197 = t227 + t344;
t264 = qJD(5) * t236 + qJD(6) * t197 + t136 * t227 + t372 + t49;
t260 = -t217 * t77 + (-t206 + t74) * t136;
t258 = t274 + t343;
t257 = t109 * t74 + t145 * t24 - t286;
t253 = t291 * t314 - t286;
t43 = t103 * qJD(4) + t241 * t111 - t359 * t112;
t230 = -pkin(8) - t351;
t218 = -pkin(4) - t352;
t190 = t220 - t352;
t188 = qJ(5) * t329;
t187 = qJ(5) * t330;
t185 = -qJDD(1) * t306 + qJDD(2);
t150 = t216 * t325 + t328;
t149 = -t216 * t327 + t326;
t148 = -t216 * t326 + t327;
t147 = t216 * t328 + t325;
t108 = t175 * t145;
t107 = t177 * t145;
t63 = pkin(5) * t336 - t268;
t56 = t84 - t399;
t36 = t109 * t177 + t315 * t335 - t316 * t336;
t35 = -t109 * t175 - t145 * t164;
t29 = pkin(5) * t342 + t43;
t10 = -pkin(9) * t342 + t15;
t9 = pkin(5) * t110 - t109 * t227 + t14;
t1 = [qJDD(1), t368, t286, t258 * t239, -t258 * t237, t291 * t367 + t253, pkin(1) * t274 + (t319 * t367 + t253) * qJ(2), t143 * t178 - t162 * t165, -t143 * t176 + t165 * t161 - t162 * t166 + t178 * t263, -qJD(3) * t165 + qJDD(3) * t178, -qJD(3) * t166 - qJDD(3) * t176, 0, qJD(3) * t250 + qJDD(3) * t294 + t186 * t166 + t185 * t176 + t225 * t368 + t263 * t306, -t256 * qJD(3) - t321 * qJDD(3) - t143 * t306 - t186 * t165 + t185 * t178 - t223 * t368, t109 * t267 + t145 * t76, -t109 * t136 - t110 * t267 - t145 * t77 + t266 * t76, t109 * t235 + t145 * t229, -t110 * t235 + t229 * t266, 0, t110 * t146 - t130 * t266 + t136 * t358 + t152 * t77 + t229 * t268 - t235 * t43 + t271, -t103 * t229 + t109 * t146 + t130 * t145 + t152 * t76 - t235 * t42 + t267 * t358 - t370, t38 * t110 + t43 * t124 + t14 * t136 + t236 * t257 + t238 * t271 - t266 * t5 - t268 * t64 + t51 * t77, -t39 * t110 + t43 * t126 - t15 * t136 + t238 * t257 + t266 * t6 - t268 * t65 - t331 * t368 - t52 * t77, -t124 * t15 - t126 * t14 - t51 * t65 - t52 * t64 + t370 + t283 * t145 + (-t236 * t39 - t238 * t38) * t109, -t24 * t268 + t38 * t14 + t39 * t15 + t74 * t43 + t5 * t51 + t6 * t52 + (g(1) * t230 - g(2) * t273) * t246 + (g(1) * t273 + g(2) * t230) * t243, -t108 * t19 - t277 * t35, -t107 * t19 + t108 * t20 + t277 * t36 - t35 * t86, -t108 * t73 - t110 * t277 + t134 * t35 - t19 * t266, -t107 * t73 - t110 * t86 - t134 * t36 + t20 * t266, t110 * t134 - t266 * t73 (-t10 * t240 + t244 * t9) * t134 + t279 * t73 - t305 * t266 - t280 * t110 + t29 * t86 + t63 * t20 + t13 * t107 + t55 * t36 - g(1) * t148 - g(2) * t150 + (-t134 * t278 + t266 * t8) * qJD(6) -(t10 * t244 + t240 * t9) * t134 - t278 * t73 + t284 * t266 - t8 * t110 - t29 * t277 + t63 * t19 - t13 * t108 + t55 * t35 - g(1) * t147 - g(2) * t149 + (-t134 * t279 - t266 * t280) * qJD(6); 0, 0, 0, -t311, t312, -t297, -qJ(2) * t297 - t274, 0, 0, 0, 0, 0, -t263 + t318 (-t161 - t304) * qJD(3) + t307, 0, 0, 0, 0, 0, t77 + t341, t76 - t339, -t124 * t267 - t136 * t386 + t238 * t77, -t126 * t267 - t136 * t337 - t347, -t236 * t64 - t238 * t65 - (t124 * t238 - t126 * t236) * t136, -t136 * t281 - t267 * t74 - t283 - t368, 0, 0, 0, 0, 0, t288 - t349, t348 - t402; 0, 0, 0, 0, 0, 0, 0, t162 * t161, -t161 ^ 2 + t162 ^ 2 (t161 - t304) * qJD(3) + t307, t263 + t318, qJDD(3), -g(3) * t225 - t186 * t162 + t223 * t286 + t295, g(3) * t223 + t186 * t161 + t286 * t225 - t276, t383, t378, t388, t380, t229, t89 * t235 + (-t136 * t162 + t359 * t229 - t235 * t317) * pkin(3) + t361, t90 * t235 + (-t162 * t267 - t229 * t241 - t235 * t303) * pkin(3) + t375, t124 * t289 - t136 * t44 + t220 * t64 + t236 * t260 + t362, t126 * t289 + t136 * t45 + t220 * t65 + t238 * t260 + t363, t124 * t45 + t126 * t44 + (-t124 * t206 - t217 * t64 - t389) * t238 + (t126 * t206 + t217 * t65 + t376) * t236 + t293, t24 * t220 - t39 * t45 - t38 * t44 - g(1) * (t246 * t287 + t188) - g(2) * (t243 * t287 + t187) - g(3) * t381 + t289 * t74 + t282 * t217 - t281 * t206, t392, t391, t393, t397, -t371 (t169 * t244 - t170 * t240) * t73 + t190 * t20 + t290 * t86 + (t240 * t270 - t244 * t269) * t134 + t395 -(t169 * t240 + t170 * t244) * t73 + t190 * t19 - t290 * t277 + (t240 * t269 + t244 * t270) * t134 + t394; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t383, t378, t388, t380, t229, t235 * t84 + t361, t83 * t235 + t375, -qJ(5) * t347 - pkin(4) * t64 - t124 * t84 + (t236 * t322 - t49) * t136 + t362, -t77 * t344 - pkin(4) * t65 - t126 * t84 + (t238 * t322 + t50) * t136 + t363, t124 * t50 + t126 * t49 + (-qJ(5) * t64 - qJD(5) * t124 - t389) * t238 + (qJ(5) * t65 + qJD(5) * t126 + t376) * t236 + t293, -t24 * pkin(4) - t39 * t50 - t38 * t49 - t74 * t84 - g(1) * (-pkin(4) * t332 + t188) - g(2) * (-pkin(4) * t333 + t187) - g(3) * t320 - t281 * qJD(5) + t282 * qJ(5), t392, t391, t393, t397, -t371 (t195 * t244 - t197 * t240) * t73 + t218 * t20 - t56 * t86 + (t240 * t265 - t244 * t264) * t134 + t395 -(t195 * t240 + t197 * t244) * t73 + t218 * t19 + t56 * t277 + (t240 * t264 + t244 * t265) * t134 + t394; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126 * t136 + t64, -t124 * t136 + t65, -t124 ^ 2 - t126 ^ 2, t124 * t39 + t126 * t38 - t259 + t366, 0, 0, 0, 0, 0, t20 - t387, t19 - t390; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t277 * t86, t277 ^ 2 - t86 ^ 2, t19 + t390, -t20 - t387, t73, -g(1) * t149 + g(2) * t147 + t215 * t354 + t55 * t277 - t396 * t8 + t305, g(1) * t150 - g(2) * t148 + t215 * t353 + t280 * t396 + t55 * t86 - t284;];
tau_reg  = t1;
