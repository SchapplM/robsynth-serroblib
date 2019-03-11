% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRPPRP5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPPRP5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP5_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:44:05
% EndTime: 2019-03-09 08:44:19
% DurationCPUTime: 7.74s
% Computational Cost: add. (7680->604), mult. (16146->711), div. (0->0), fcn. (10519->10), ass. (0->312)
t240 = cos(qJ(2));
t349 = qJD(1) * qJD(2);
t328 = t240 * t349;
t238 = sin(qJ(2));
t348 = t238 * qJDD(1);
t273 = t328 + t348;
t148 = qJDD(5) + t273;
t361 = qJD(1) * t238;
t198 = qJD(5) + t361;
t233 = sin(pkin(9));
t234 = cos(pkin(9));
t237 = sin(qJ(5));
t412 = cos(qJ(5));
t280 = -t233 * t412 - t237 * t234;
t331 = qJD(5) * t412;
t356 = qJD(5) * t237;
t135 = -t233 * t356 + t234 * t331;
t339 = t412 * t234;
t315 = t238 * t339;
t338 = t233 * t361;
t368 = -qJD(1) * t315 + t237 * t338 - t135;
t274 = t148 * t280 + t198 * t368;
t360 = qJD(1) * t240;
t332 = t234 * t360;
t354 = t233 * qJD(2);
t145 = -t332 - t354;
t333 = t233 * t360;
t359 = qJD(2) * t234;
t146 = -t333 + t359;
t85 = -t412 * t145 + t146 * t237;
t434 = t85 * t360 + t274;
t282 = t237 * t145 + t146 * t412;
t433 = -qJD(2) * t282 + t274;
t150 = -t237 * t233 + t339;
t134 = t280 * qJD(5);
t265 = t238 * t280;
t426 = qJD(1) * t265 + t134;
t308 = t150 * t148 + t198 * t426;
t432 = -qJD(2) * t85 + t308;
t227 = pkin(9) + qJ(5);
t217 = cos(t227);
t239 = sin(qJ(1));
t241 = cos(qJ(1));
t302 = g(1) * t241 + g(2) * t239;
t284 = t302 * t240;
t405 = g(3) * t238;
t425 = t284 + t405;
t402 = pkin(2) + qJ(4);
t401 = -pkin(8) - t402;
t158 = t401 * t233;
t159 = t401 * t234;
t95 = t158 * t412 + t237 * t159;
t431 = t95 * t148 + t217 * t425;
t347 = t240 * qJDD(1);
t183 = t233 * t347;
t318 = t234 * qJDD(2) - t183;
t329 = t238 * t349;
t109 = t233 * t329 + t318;
t185 = t234 * t347;
t214 = t233 * qJDD(2);
t366 = t185 + t214;
t272 = t234 * t329 - t366;
t40 = qJD(5) * t282 + t237 * t109 - t412 * t272;
t277 = -t40 * t280 - t368 * t85;
t430 = t85 ^ 2;
t414 = pkin(3) + pkin(7);
t429 = t198 * t85;
t209 = pkin(7) * t361;
t424 = qJD(3) + t209;
t415 = t282 ^ 2;
t379 = t233 * t238;
t290 = pkin(4) * t240 - pkin(8) * t379;
t213 = pkin(2) * t361;
t387 = qJ(3) * t240;
t294 = qJ(4) * t238 - t387;
t127 = qJD(1) * t294 + t213;
t210 = pkin(7) * t360;
t156 = pkin(3) * t360 + t210;
t78 = -t233 * t127 + t234 * t156;
t59 = qJD(1) * t290 + t78;
t337 = t234 * t361;
t79 = t234 * t127 + t233 * t156;
t69 = pkin(8) * t337 + t79;
t397 = -qJD(4) * t150 - qJD(5) * t95 + t237 * t69 - t412 * t59;
t390 = t282 * t198;
t423 = t40 - t390;
t218 = t233 * pkin(4);
t203 = qJ(3) + t218;
t235 = -pkin(8) - qJ(4);
t346 = t240 * t218;
t376 = t238 * t239;
t422 = t235 * t376 + t239 * t346;
t375 = t238 * t241;
t421 = t235 * t375 + t241 * t346;
t205 = pkin(4) * t234 + pkin(3);
t352 = t205 * t361 + t424;
t281 = -t237 * t158 + t159 * t412;
t39 = -t412 * t109 - t145 * t331 + t146 * t356 - t237 * t272;
t62 = qJD(4) * t280 + qJD(5) * t281;
t420 = t281 * t39 - t95 * t40 - t62 * t85;
t419 = t40 * pkin(5) + t39 * qJ(6) - t282 * qJD(6);
t11 = -t150 * t39 + t282 * t426;
t219 = t238 * qJ(3);
t223 = t240 * pkin(2);
t365 = t223 + t219;
t386 = qJ(4) * t240;
t299 = t365 + t386;
t142 = -pkin(1) - t299;
t175 = t414 * t238;
t152 = t234 * t175;
t71 = t238 * pkin(4) + t152 + (pkin(8) * t240 - t142) * t233;
t378 = t234 * t240;
t91 = t234 * t142 + t233 * t175;
t77 = -pkin(8) * t378 + t91;
t395 = t237 * t71 + t412 * t77;
t357 = qJD(2) * t240;
t157 = t414 * t357;
t358 = qJD(2) * t238;
t212 = pkin(2) * t358;
t353 = t238 * qJD(3);
t249 = qJD(2) * t294 - t240 * qJD(4) - t353;
t97 = t212 + t249;
t67 = t234 * t157 - t233 * t97;
t50 = qJD(2) * t290 + t67;
t68 = t233 * t157 + t234 * t97;
t55 = pkin(8) * t234 * t358 + t68;
t9 = -qJD(5) * t395 - t237 * t55 + t412 * t50;
t197 = pkin(2) * t329;
t325 = -pkin(1) - t219;
t270 = -t240 * t402 + t325;
t56 = qJD(1) * t249 + qJDD(1) * t270 + t197;
t196 = pkin(7) * t328;
t206 = pkin(7) * t348;
t326 = qJDD(3) + t196 + t206;
t80 = pkin(3) * t273 - qJD(2) * qJD(4) - qJDD(2) * t402 + t326;
t33 = -t233 * t56 + t234 * t80;
t19 = pkin(4) * t273 - t109 * pkin(8) + t33;
t34 = t233 * t80 + t234 * t56;
t24 = pkin(8) * t272 + t34;
t113 = t270 * qJD(1);
t351 = pkin(3) * t361 + t424;
t122 = -qJD(2) * t402 + t351;
t64 = -t113 * t233 + t234 * t122;
t45 = pkin(4) * t361 - pkin(8) * t146 + t64;
t65 = t234 * t113 + t233 * t122;
t49 = pkin(8) * t145 + t65;
t3 = t237 * t19 + t412 * t24 + t45 * t331 - t356 * t49;
t381 = t148 * qJ(6);
t1 = qJD(6) * t198 + t3 + t381;
t17 = -t237 * t49 + t412 * t45;
t369 = qJD(6) - t17;
t13 = -t198 * pkin(5) + t369;
t18 = t237 * t45 + t412 * t49;
t14 = t198 * qJ(6) + t18;
t324 = -t412 * t19 + t237 * t24 + t49 * t331 + t45 * t356;
t410 = pkin(5) * t148;
t2 = qJDD(6) + t324 - t410;
t226 = g(3) * t240;
t341 = -g(1) * t375 - g(2) * t376 + t226;
t418 = -t1 * t280 - t13 * t426 - t14 * t368 - t2 * t150 + t341;
t417 = -t150 * t324 + t17 * t426 - t18 * t368 - t280 * t3 + t341;
t416 = -t150 * t40 - t280 * t39 + t282 * t368 - t426 * t85;
t231 = t238 ^ 2;
t411 = pkin(2) * t238;
t409 = g(1) * t239;
t406 = g(2) * t241;
t404 = t282 * t85;
t400 = -pkin(5) * t368 - qJ(6) * t426 - qJD(6) * t150 + t352;
t30 = t237 * t59 + t412 * t69;
t27 = qJ(6) * t360 + t30;
t399 = t62 - t27;
t398 = t62 - t30;
t396 = pkin(5) * t360 - t397;
t393 = t18 * t198;
t392 = t233 * t64;
t388 = pkin(7) * qJDD(2);
t383 = qJDD(2) * pkin(2);
t382 = t109 * t234;
t243 = qJD(1) ^ 2;
t380 = t231 * t243;
t374 = t238 * t243;
t373 = t239 * t240;
t372 = t240 * t235;
t371 = t240 * t241;
t176 = t414 * t240;
t364 = t241 * pkin(1) + t239 * pkin(7);
t232 = t240 ^ 2;
t363 = t231 - t232;
t362 = t231 + t232;
t230 = qJD(2) * qJ(3);
t132 = qJD(4) + t230 + t156;
t355 = t132 * qJD(2);
t350 = qJD(4) - t132;
t191 = pkin(4) * t379;
t345 = -t415 + t430;
t224 = t241 * pkin(7);
t343 = t241 * t205 + t239 * t372 + t224;
t207 = pkin(7) * t347;
t228 = qJDD(2) * qJ(3);
t229 = qJD(2) * qJD(3);
t342 = t207 + t228 + t229;
t133 = pkin(4) * t378 + t176;
t340 = t414 * qJD(2);
t336 = t402 * t357;
t334 = t198 * t360;
t330 = t366 * pkin(4);
t327 = t233 * t348;
t184 = t234 * t348;
t323 = qJD(1) * t91 + t65;
t322 = -qJD(2) * pkin(2) + qJD(3);
t320 = t145 + t354;
t319 = -t146 + t359;
t317 = pkin(2) * t371 + qJ(3) * t375 + t364;
t316 = -t206 - t341;
t313 = qJD(1) * t340;
t312 = t238 * t328;
t192 = qJ(3) * t373;
t311 = -pkin(2) * t376 + t192;
t195 = qJ(3) * t371;
t310 = -pkin(2) * t375 + t195;
t201 = g(1) * t373;
t309 = -g(2) * t371 + t201;
t306 = t362 * qJDD(1) * pkin(7);
t242 = qJD(2) ^ 2;
t305 = pkin(7) * t242 + t406;
t216 = sin(t227);
t118 = t216 * t239 - t217 * t375;
t120 = t216 * t241 + t217 * t376;
t304 = -g(1) * t120 - g(2) * t118;
t119 = t216 * t375 + t217 * t239;
t121 = -t216 * t376 + t217 * t241;
t303 = -g(1) * t121 - g(2) * t119;
t301 = -t406 + t409;
t300 = -t415 - t430;
t123 = t150 * t240;
t74 = -t237 * t238 * t354 + qJD(2) * t315 - t134 * t240;
t298 = t123 * t40 - t74 * t85;
t297 = pkin(5) * t216 - qJ(6) * t217;
t295 = t34 * t233 + t33 * t234;
t161 = t209 + t322;
t172 = -t210 - t230;
t293 = t161 * t240 + t172 * t238;
t292 = pkin(3) * t347 + qJDD(4) + t342;
t291 = (-pkin(7) - t205) * qJD(2);
t289 = t325 - t223;
t288 = -t233 * t380 + t234 * t328 + t184;
t287 = t191 + t365 - t372;
t37 = -t237 * t77 + t412 * t71;
t283 = -0.2e1 * pkin(1) * t349 - t388;
t8 = t237 * t50 + t71 * t331 - t356 * t77 + t412 * t55;
t140 = t289 * qJD(1);
t279 = t140 * t361 + qJDD(3) - t316;
t278 = -t109 * t240 + t146 * t358;
t116 = t238 * t291;
t276 = -qJ(3) * t357 - t353;
t92 = -pkin(4) * t145 + t132;
t275 = qJD(1) * t291;
t83 = -t238 * t313 + t292;
t271 = -t240 * t83 + t302;
t269 = 0.2e1 * qJDD(1) * pkin(1) - t305;
t268 = -t234 * t380 - t327;
t162 = -pkin(1) - t365;
t264 = t388 + (-qJD(1) * t162 - t140) * qJD(2);
t262 = g(1) * t118 - g(2) * t120 + t217 * t226 - t324;
t261 = t39 + t429;
t124 = t280 * t240;
t73 = -qJD(2) * t265 - t135 * t240;
t259 = -t123 * t39 + t124 * t40 - t282 * t74 + t73 * t85;
t258 = t241 * t191 + t239 * t205 - t235 * t371 + t317;
t257 = t83 - t425;
t256 = -t11 - t277;
t255 = (-t203 * t238 - pkin(1) - t223) * t409;
t254 = t123 * t148 - t198 * t74 + t238 * t40 + t357 * t85;
t129 = t212 + t276;
t82 = qJD(1) * t276 + qJDD(1) * t289 + t197;
t253 = qJD(1) * t129 + qJDD(1) * t162 + t305 + t82;
t252 = -t284 + t292;
t114 = pkin(7) * t329 - t342;
t126 = t326 - t383;
t251 = qJD(2) * t293 - t114 * t240 + t126 * t238;
t31 = pkin(5) * t85 - qJ(6) * t282 + t92;
t250 = t282 * t31 + qJDD(6) - t262;
t248 = t148 * t281 - t216 * t425;
t247 = -g(1) * t119 + g(2) * t121 + t216 * t226 + t3;
t246 = t40 + t390;
t53 = t238 * t275 + t292 + t330;
t245 = (-g(3) + t275) * t238 + t330 + t252;
t190 = t240 * t374;
t168 = t363 * t243;
t164 = qJDD(2) * t240 - t238 * t242;
t163 = qJDD(2) * t238 + t240 * t242;
t155 = t238 * t340;
t153 = -qJ(3) * t360 + t213;
t144 = qJDD(1) * t232 - 0.2e1 * t312;
t143 = qJDD(1) * t231 + 0.2e1 * t312;
t108 = 0.2e1 * t238 * t347 - 0.2e1 * t349 * t363;
t96 = t148 * t238 + t198 * t357;
t90 = -t142 * t233 + t152;
t81 = -pkin(5) * t280 - qJ(6) * t150 + t203;
t54 = pkin(5) * t123 - qJ(6) * t124 + t133;
t44 = pkin(5) * t282 + qJ(6) * t85;
t36 = -t238 * pkin(5) - t37;
t35 = qJ(6) * t238 + t395;
t26 = -t282 * t360 + t308;
t23 = -t74 * pkin(5) - t73 * qJ(6) - t124 * qJD(6) + t116;
t20 = -t39 + t429;
t12 = -t124 * t39 + t282 * t73;
t10 = t124 * t148 + t198 * t73 - t238 * t39 + t282 * t357;
t7 = -pkin(5) * t357 - t9;
t6 = qJ(6) * t357 + qJD(6) * t238 + t8;
t5 = t53 + t419;
t4 = [0, 0, 0, 0, 0, qJDD(1), t301, t302, 0, 0, t143, t108, t163, t144, t164, 0, t238 * t283 + t240 * t269 + t201, t283 * t240 + (-t269 - t409) * t238, -t302 + 0.2e1 * t306, -g(1) * (-t239 * pkin(1) + t224) - g(2) * t364 + (pkin(7) ^ 2 * t362 + pkin(1) ^ 2) * qJDD(1), 0, -t163, -t164, t143, t108, t144, t306 + t251 - t302, t238 * t264 + t240 * t253 - t201, t264 * t240 + (-t253 + t409) * t238, pkin(7) * t251 - g(1) * t224 - g(2) * t317 + t140 * t129 + t82 * t162 - t289 * t409, t278 * t233, t278 * t234 + (t145 * t358 - t240 * t272) * t233 (t109 - t183) * t238 + (qJD(1) * t233 * t363 + t146 * t240) * qJD(2) (t366 * t240 + (t145 - t332) * t358) * t234 (-0.2e1 * t185 - t214) * t238 + (t240 * t145 + (t363 + t231) * t234 * qJD(1)) * qJD(2), t143, t155 * t145 + t176 * t366 + (qJD(1) * t90 + t64) * t357 - t271 * t234 + (-t234 * t355 + t90 * qJDD(1) + t33 + t301 * t233 + (-t176 * t359 + t67) * qJD(1)) * t238, t176 * t109 - t155 * t146 - t323 * t357 + t271 * t233 + (-t68 * qJD(1) - t91 * qJDD(1) + t132 * t354 + t234 * t301 - t34) * t238, t68 * t145 - t91 * t366 - t67 * t146 - t90 * t109 + t201 + (t234 * t323 - t392) * t358 + (t233 * t33 - t234 * t34 - t406) * t240, t34 * t91 + t65 * t68 + t33 * t90 + t64 * t67 + t83 * t176 - t132 * t155 - g(1) * (t241 * pkin(3) + t224) - g(2) * (qJ(4) * t371 + t317) + (-g(1) * (t289 - t386) - g(2) * pkin(3)) * t239, t12, -t259, t10, t298, -t254, t96, t116 * t85 + t123 * t53 + t133 * t40 + t148 * t37 + t17 * t357 + t198 * t9 - t238 * t324 - t74 * t92 + t303, t116 * t282 + t124 * t53 - t133 * t39 - t148 * t395 - t18 * t357 - t198 * t8 - t238 * t3 + t73 * t92 - t304, -t123 * t3 + t124 * t324 - t17 * t73 + t18 * t74 - t282 * t9 + t37 * t39 - t395 * t40 - t8 * t85 + t309, -g(1) * t343 - g(2) * t258 + t92 * t116 + t53 * t133 + t17 * t9 + t18 * t8 + t3 * t395 - t324 * t37 - t255, t12, t10, t259, t96, t254, t298, t123 * t5 - t13 * t357 - t148 * t36 - t198 * t7 - t2 * t238 + t23 * t85 - t31 * t74 + t40 * t54 + t303, -t1 * t123 + t124 * t2 + t13 * t73 + t14 * t74 + t282 * t7 - t35 * t40 - t36 * t39 - t6 * t85 + t309, t1 * t238 - t124 * t5 + t14 * t357 + t148 * t35 + t198 * t6 - t23 * t282 - t31 * t73 + t39 * t54 + t304, t1 * t35 + t14 * t6 + t5 * t54 + t31 * t23 + t2 * t36 + t13 * t7 - g(1) * (t121 * pkin(5) + t120 * qJ(6) + t343) - g(2) * (pkin(5) * t119 + qJ(6) * t118 + t258) - t255; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t190, t168, t348, t190, t347, qJDD(2), pkin(1) * t374 + t316, t405 - t207 + (pkin(1) * t243 + t302) * t240, 0, 0, qJDD(2), -t348, -t347, -t190, t168, t190 (t387 - t411) * qJDD(1) + ((-t172 - t230) * t238 + (-t161 + t322) * t240) * qJD(1), -t153 * t360 + t279 - 0.2e1 * t383, t207 + 0.2e1 * t228 + 0.2e1 * t229 + (qJD(1) * t153 - g(3)) * t238 + (qJD(1) * t140 - t302) * t240, -pkin(7) * qJD(1) * t293 - t126 * pkin(2) - g(1) * t310 - g(2) * t311 - g(3) * t365 - t114 * qJ(3) - t172 * qJD(3) - t140 * t153, -t146 * t338 + t382, -t366 * t234 - t109 * t233 + (-t145 * t233 + t234 * t319) * t361, -t146 * t360 + t288, -t145 * t337 - t233 * t272, -t320 * t360 + t268, -t190, -t402 * t184 + qJ(3) * t366 - t351 * t145 + t257 * t233 + (-t78 * t238 - t64 * t240 + (-t336 + (-t350 - t230) * t238) * t234) * qJD(1), t402 * t327 + qJ(3) * t109 + t351 * t146 + t257 * t234 + (t238 * t79 + t65 * t240 + (t238 * t350 + t336) * t233) * qJD(1), -t79 * t145 + t78 * t146 + (qJD(4) * t146 + t109 * t402 - t361 * t65 - t33) * t234 + (-qJD(4) * t145 - t272 * t402 + t361 * t64 - t34) * t233 - t341, t83 * qJ(3) - t65 * t79 - t64 * t78 - g(1) * t195 - g(2) * t192 - g(3) * t299 + t351 * t132 + (-t233 * t65 - t234 * t64) * qJD(4) + (t238 * t302 - t295) * t402, t11, t416, t26, t277, t434, -t334, -t17 * t360 + t198 * t397 + t203 * t40 - t280 * t53 + t352 * t85 - t368 * t92 + t248, t53 * t150 + t18 * t360 - t198 * t398 - t203 * t39 + t282 * t352 + t426 * t92 - t431, -t282 * t397 + t30 * t85 - t417 + t420, t3 * t95 - t324 * t281 + t53 * t203 - g(1) * (t310 + t421) - g(2) * (t311 + t422) - g(3) * t287 + t352 * t92 + t398 * t18 + t397 * t17, t11, t26, -t416, -t334, -t434, t277, t13 * t360 - t198 * t396 - t280 * t5 - t31 * t368 + t81 * t40 + t400 * t85 + t248, t27 * t85 + t282 * t396 - t418 + t420, -t14 * t360 - t5 * t150 + t399 * t198 - t400 * t282 - t31 * t426 + t81 * t39 + t431, t1 * t95 + t5 * t81 - t2 * t281 - g(1) * (t195 + t421) - g(2) * (t192 + t422) - g(3) * (t238 * t297 + t287) + t400 * t31 + t399 * t14 + t396 * t13 + t302 * (-t240 * t297 + t411); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t348, qJDD(2) + t190, -t242 - t380, qJD(2) * t172 + t196 + t279 - t383, 0, 0, 0, 0, 0, 0, qJD(2) * t145 + t288 (-t146 - t333) * qJD(2) + t268, -t366 * t233 - t382 + (t145 * t234 + (t146 + t359) * t233) * t361, -t355 + (t234 * t65 - t392) * t361 + t295 + t341, 0, 0, 0, 0, 0, 0, t432, t433, t256, -t92 * qJD(2) + t417, 0, 0, 0, 0, 0, 0, t432, t256, -t433, -t31 * qJD(2) + t418; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t319 * t361 + t366, t320 * t361 + t318, -t145 ^ 2 - t146 ^ 2, -t65 * t145 + t64 * t146 + (-g(3) - t313) * t238 + t252, 0, 0, 0, 0, 0, 0, t246, -t261, t300, t17 * t282 + t18 * t85 + t245, 0, 0, 0, 0, 0, 0, t246, t300, t261, -t13 * t282 + t14 * t85 + t245 + t419; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t404, -t345, t20, -t404, -t423, t148, -t282 * t92 + t262 + t393, t17 * t198 + t85 * t92 - t247, 0, 0, t404, t20, t345, t148, t423, -t404, -t44 * t85 - t250 + t393 + 0.2e1 * t410, pkin(5) * t39 - qJ(6) * t40 + (t14 - t18) * t282 + (t13 - t369) * t85, 0.2e1 * t381 - t31 * t85 + t44 * t282 + (0.2e1 * qJD(6) - t17) * t198 + t247, t1 * qJ(6) - t2 * pkin(5) - t31 * t44 - t13 * t18 - g(1) * (-pkin(5) * t118 + qJ(6) * t119) - g(2) * (pkin(5) * t120 - qJ(6) * t121) - (-pkin(5) * t217 - qJ(6) * t216) * t226 + t369 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148 + t404, t20, -t198 ^ 2 - t415, -t14 * t198 + t250 - t410;];
tau_reg  = t4;
