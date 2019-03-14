% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPRRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRPP5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP5_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:44:49
% EndTime: 2019-03-09 04:45:03
% DurationCPUTime: 7.09s
% Computational Cost: add. (8118->590), mult. (19153->657), div. (0->0), fcn. (14378->10), ass. (0->278)
t219 = sin(pkin(9));
t223 = sin(qJ(3));
t220 = cos(pkin(9));
t392 = cos(qJ(3));
t309 = t392 * t220;
t260 = -t223 * t219 + t309;
t433 = t260 * qJD(1);
t146 = qJD(4) - t433;
t435 = t433 * qJD(3);
t168 = t219 * t392 + t223 * t220;
t159 = t168 * qJD(1);
t225 = cos(qJ(4));
t222 = sin(qJ(4));
t327 = qJD(3) * t222;
t121 = t159 * t225 + t327;
t358 = t121 * t222;
t119 = -t225 * qJD(3) + t159 * t222;
t362 = t119 * t225;
t275 = t358 + t362;
t325 = qJD(4) * t222;
t247 = t168 * qJDD(1);
t233 = t247 + t435;
t431 = qJD(3) * qJD(4) + t233;
t255 = t222 * qJDD(3) - t159 * t325 + t431 * t225;
t301 = t121 * t325 - t225 * t255;
t324 = qJD(4) * t225;
t284 = -t225 * qJDD(3) + t159 * t324;
t70 = t222 * t247 + (qJD(4) + t433) * t327 + t284;
t373 = t119 * t324 + t222 * t70;
t7 = -t275 * t433 + t301 + t373;
t296 = t146 * t222;
t434 = t373 - t121 * t296 + (-t119 * t433 + t255) * t225;
t162 = t168 * qJD(3);
t319 = t219 * qJDD(1);
t279 = -qJDD(1) * t309 + t223 * t319;
t104 = qJD(1) * t162 + t279;
t98 = qJDD(4) + t104;
t432 = t98 * qJ(5) + t146 * qJD(5);
t217 = pkin(9) + qJ(3);
t209 = sin(t217);
t226 = cos(qJ(1));
t347 = t209 * t226;
t224 = sin(qJ(1));
t349 = t209 * t224;
t416 = g(1) * t347 + g(2) * t349;
t366 = qJDD(1) * pkin(1);
t208 = qJDD(2) - t366;
t415 = g(1) * t224 - g(2) * t226;
t269 = -t208 + t415;
t294 = t146 * t225;
t430 = t146 * t324 + t222 * t98 - t294 * t433;
t395 = pkin(4) + pkin(5);
t315 = t395 * t222;
t367 = qJ(5) * t225;
t267 = t315 - t367;
t363 = t119 * t159;
t356 = t433 * t222;
t407 = t225 * t98 + (-t325 + t356) * t146;
t242 = t407 - t363;
t429 = t407 + t363;
t304 = qJD(3) * t392;
t326 = qJD(3) * t223;
t161 = t219 * t326 - t220 * t304;
t374 = t70 * t225;
t375 = t255 * t222;
t427 = ((t119 * t222 - t121 * t225) * qJD(4) - t374 - t375) * t168 + t275 * t161;
t343 = t222 * qJ(5);
t265 = -t225 * t395 - t343;
t165 = pkin(3) - t265;
t135 = t146 * qJ(5);
t385 = pkin(7) + qJ(2);
t178 = t385 * t219;
t169 = qJD(1) * t178;
t179 = t385 * t220;
t170 = qJD(1) * t179;
t106 = -t223 * t169 + t170 * t392;
t100 = qJD(3) * pkin(8) + t106;
t204 = t220 * pkin(2) + pkin(1);
t176 = -qJD(1) * t204 + qJD(2);
t78 = -pkin(3) * t433 - t159 * pkin(8) + t176;
t41 = t100 * t225 + t222 * t78;
t35 = t135 + t41;
t380 = t146 * t35;
t321 = qJD(1) * qJD(2);
t399 = qJDD(1) * t385 + t321;
t132 = t399 * t219;
t133 = t399 * t220;
t313 = t223 * t132 - t392 * t133 + t169 * t304;
t53 = -t170 * t326 - t313;
t49 = qJDD(3) * pkin(8) + t53;
t318 = t220 * qJDD(1);
t52 = -pkin(2) * t318 + t104 * pkin(3) - pkin(8) * t233 + t208;
t302 = t100 * t324 + t222 * t49 - t225 * t52 + t78 * t325;
t94 = t98 * pkin(4);
t409 = t94 - qJDD(5);
t6 = t302 - t409;
t426 = -t6 + t380;
t378 = t146 * t41;
t425 = -t302 + t378;
t424 = t70 * qJ(6) + t119 * qJD(6);
t423 = 0.2e1 * t432;
t40 = -t222 * t100 + t225 * t78;
t334 = qJD(5) - t40;
t420 = -t222 * qJD(5) - t106;
t397 = t121 ^ 2;
t419 = -t146 ^ 2 - t397;
t418 = t162 * qJ(5) - qJD(5) * t260;
t417 = -t392 * t178 - t223 * t179;
t210 = cos(t217);
t329 = t210 * pkin(3) + t209 * pkin(8);
t286 = g(1) * t226 + g(2) * t224;
t398 = t119 ^ 2;
t413 = t397 - t398;
t114 = -t223 * t178 + t179 * t392;
t411 = t114 * qJD(3);
t410 = t159 * qJD(3);
t408 = qJ(2) * qJDD(1);
t359 = t121 * t159;
t25 = -t359 + t430;
t248 = t359 + t430;
t8 = -t100 * t325 + t222 * t52 + t225 * t49 + t78 * t324;
t5 = t8 + t432;
t2 = t5 + t424;
t30 = t121 * qJ(6) + t40;
t335 = qJD(5) - t30;
t21 = -t146 * t395 + t335;
t406 = t146 * t21 + t2;
t393 = pkin(8) * t98;
t155 = t223 * t170;
t105 = -t392 * t169 - t155;
t289 = qJD(3) * pkin(3) + t105;
t251 = t121 * qJ(5) + t289;
t42 = pkin(4) * t119 - t251;
t405 = t146 * t42 - t393;
t338 = t225 * t226;
t342 = t222 * t224;
t147 = t210 * t342 + t338;
t336 = t226 * t222;
t339 = t224 * t225;
t149 = t210 * t336 - t339;
t350 = t209 * t222;
t243 = g(1) * t149 + g(2) * t147 + g(3) * t350 - t302;
t241 = t243 + t409;
t29 = -t119 * t395 + qJD(6) + t251;
t382 = qJ(6) * t255;
t404 = (qJD(6) + t29) * t121 + t241 + t382;
t308 = t168 * t325;
t353 = t161 * t225;
t258 = t308 + t353;
t351 = t168 * t225;
t403 = -t121 * t162 + t146 * t258 + t255 * t260 - t351 * t98;
t307 = t168 * t324;
t354 = t161 * t222;
t259 = t307 - t354;
t352 = t168 * t222;
t402 = t162 * t119 + t146 * t259 - t260 * t70 + t352 * t98;
t396 = t159 ^ 2;
t394 = pkin(5) * t98;
t391 = pkin(8) * t121;
t388 = g(2) * t385;
t386 = g(3) * t210;
t384 = pkin(8) - qJ(6);
t101 = pkin(3) * t159 - pkin(8) * t433;
t58 = t222 * t101 + t225 * t105;
t383 = qJ(5) * t70;
t31 = qJ(6) * t119 + t41;
t26 = t135 + t31;
t381 = t146 * t26;
t379 = t146 * t40;
t103 = -pkin(3) * t260 - pkin(8) * t168 - t204;
t68 = t222 * t103 + t225 * t114;
t371 = -t146 * t267 - t420;
t37 = t159 * qJ(5) + t58;
t370 = -qJ(6) * t356 - t225 * qJD(6) - t325 * t384 - t37;
t280 = pkin(4) * t222 - t367;
t369 = t146 * t280 + t420;
t181 = t384 * t225;
t96 = t222 * t105;
t368 = qJD(4) * t181 - t222 * qJD(6) - t96 - (-qJ(6) * t433 - t101) * t225 + t395 * t159;
t365 = t119 * qJ(5);
t364 = t119 * t146;
t361 = t121 * t119;
t360 = t121 * t146;
t357 = t146 * t159;
t355 = t159 * t433;
t348 = t209 * t225;
t346 = t210 * t224;
t345 = t210 * t225;
t344 = t210 * t226;
t337 = t226 * t385;
t332 = t416 * t222;
t331 = t416 * t225;
t330 = g(1) * t349 - g(2) * t347;
t215 = t219 ^ 2;
t216 = t220 ^ 2;
t328 = t215 + t216;
t323 = qJD(5) * t225;
t102 = pkin(3) * t162 + pkin(8) * t161;
t79 = t260 * qJD(2) + t417 * qJD(3);
t317 = t222 * t102 + t103 * t324 + t225 * t79;
t316 = t103 * t325 + t114 * t324 + t222 * t79;
t47 = -qJ(5) * t260 + t68;
t184 = t226 * t204;
t312 = pkin(3) * t344 + pkin(8) * t347 + t184;
t311 = g(1) * t344 + g(2) * t346 + g(3) * t209;
t292 = t392 * t132 + t223 * t133 - t169 * t326 + t170 * t304;
t50 = -qJDD(3) * pkin(3) + t292;
t10 = t70 * pkin(4) - qJ(5) * t255 - t121 * qJD(5) + t50;
t4 = -pkin(5) * t70 + qJDD(6) - t10;
t310 = t4 - t386;
t57 = t225 * t101 - t96;
t299 = t328 * qJD(1) ^ 2;
t148 = t210 * t339 - t336;
t298 = -t147 * pkin(4) + qJ(5) * t148;
t150 = t210 * t338 + t342;
t297 = -t149 * pkin(4) + qJ(5) * t150;
t108 = t222 * t114;
t67 = t225 * t103 - t108;
t293 = -pkin(8) * t374 - t311;
t291 = pkin(4) * t345 + t210 * t343 + t329;
t290 = 0.2e1 * t328;
t288 = g(1) * t147 - g(2) * t149;
t287 = g(1) * t148 - g(2) * t150;
t283 = (qJD(4) * t119 + t255) * pkin(8);
t1 = -qJD(6) * t121 - t382 - t394 + t6;
t282 = -t1 + t381;
t281 = pkin(4) * t225 + t343;
t34 = -pkin(4) * t146 + t334;
t278 = t222 * t35 - t225 * t34;
t277 = t222 * t41 + t225 * t40;
t51 = t146 * t162 - t260 * t98;
t273 = qJ(6) * t161 - qJD(6) * t168;
t20 = t102 * t225 - t316;
t271 = pkin(3) + t281;
t270 = -t204 - t329;
t268 = -pkin(8) * qJD(4) * t146 - t386;
t262 = -t148 * pkin(4) - t147 * qJ(5) + t337;
t261 = -t146 * t289 - t393;
t19 = -t114 * t325 + t317;
t257 = -t10 + t268;
t256 = -t268 + t50;
t254 = t150 * pkin(4) + qJ(5) * t149 + t312;
t175 = -qJDD(1) * t204 + qJDD(2);
t23 = t119 * t296 - t374;
t246 = t269 + t366;
t33 = t255 + t364;
t240 = (-g(1) * t270 - t388) * t224;
t239 = t290 * t321 - t286;
t18 = t119 * t259 + t352 * t70;
t237 = t121 * t42 - t241;
t236 = g(1) * t150 + g(2) * t148 + g(3) * t348 - t8;
t235 = t236 + t379;
t80 = t168 * qJD(2) + t411;
t232 = -qJDD(4) - t279 + t361 - t410;
t231 = t431 * t222 + t284;
t230 = -t231 + t360;
t188 = pkin(8) * t344;
t185 = pkin(8) * t346;
t182 = qJ(5) * t348;
t180 = t384 * t222;
t154 = t433 ^ 2;
t72 = pkin(4) * t121 + t365;
t71 = t168 * t280 - t417;
t56 = -t121 * t395 - t365;
t55 = -t168 * t267 + t417;
t48 = pkin(4) * t260 - t67;
t39 = -t159 * pkin(4) - t57;
t36 = qJ(6) * t352 + t47;
t28 = t108 + (-qJ(6) * t168 - t103) * t225 + t395 * t260;
t24 = -t280 * t161 + (qJD(4) * t281 - t323) * t168 + t80;
t22 = t121 * t294 + t375;
t17 = -t121 * t258 + t255 * t351;
t16 = -pkin(4) * t162 - t20;
t15 = -t411 + t267 * t161 + (qJD(4) * t265 - qJD(2) + t323) * t168;
t14 = t19 + t418;
t12 = qJ(6) * t307 + (-qJD(4) * t114 - t273) * t222 + t317 + t418;
t11 = qJ(6) * t308 - t395 * t162 + (-t102 + t273) * t225 + t316;
t3 = [0, 0, 0, 0, 0, qJDD(1), t415, t286, 0, 0, t215 * qJDD(1), 0.2e1 * t219 * t318, 0, t216 * qJDD(1), 0, 0, t246 * t220, -t246 * t219, t290 * t408 + t239, t269 * pkin(1) + (t328 * t408 + t239) * qJ(2), -t159 * t161 + t168 * t233, -t168 * t104 - t159 * t162 - t161 * t433 + t233 * t260, -qJD(3) * t161 + qJDD(3) * t168, -t104 * t260 - t162 * t433, -qJD(3) * t162 + qJDD(3) * t260, 0, -t80 * qJD(3) + qJDD(3) * t417 - t204 * t104 + t176 * t162 - t175 * t260 + t210 * t415, -t79 * qJD(3) - t114 * qJDD(3) - t176 * t161 + t175 * t168 - t204 * t233 - t330, -t114 * t104 + t105 * t161 - t106 * t162 + t80 * t159 + t292 * t168 - t233 * t417 + t53 * t260 + t433 * t79 - t286, t53 * t114 + t106 * t79 - t292 * t417 - t105 * t80 - t175 * t204 - g(1) * (-t224 * t204 + t337) - g(2) * (t224 * t385 + t184) t17, t427, -t403, t18, -t402, t51, t289 * t354 - t417 * t70 + t119 * t80 + t146 * t20 + t162 * t40 + t260 * t302 + t67 * t98 + (t222 * t50 - t289 * t324) * t168 + t287, t289 * t353 - t417 * t255 + t121 * t80 - t146 * t19 - t162 * t41 + t260 * t8 - t68 * t98 + (t225 * t50 + t289 * t325) * t168 - t288, -t119 * t19 - t121 * t20 - t67 * t255 - t68 * t70 + t277 * t161 + (-t222 * t8 + t225 * t302 + (t222 * t40 - t225 * t41) * qJD(4)) * t168 + t330, -g(1) * t337 - g(2) * t312 + t41 * t19 + t40 * t20 - t289 * t80 - t302 * t67 - t417 * t50 + t8 * t68 + t240, t17, -t403, -t427, t51, t402, t18, -t42 * t354 + t119 * t24 - t146 * t16 - t162 * t34 + t260 * t6 - t48 * t98 + t70 * t71 + (t10 * t222 + t324 * t42) * t168 + t287, -t119 * t14 + t121 * t16 - t47 * t70 + t48 * t255 + t278 * t161 + (-t222 * t5 + t225 * t6 + (-t222 * t34 - t225 * t35) * qJD(4)) * t168 + t330, t42 * t353 - t121 * t24 + t14 * t146 + t162 * t35 - t260 * t5 + t47 * t98 - t255 * t71 + (-t10 * t225 + t325 * t42) * t168 + t288, -g(1) * t262 - g(2) * t254 + t10 * t71 + t35 * t14 + t34 * t16 + t42 * t24 + t5 * t47 + t6 * t48 + t240, t17, -t427, t403, t18, -t402, t51, t29 * t354 + t1 * t260 - t11 * t146 - t119 * t15 - t162 * t21 - t28 * t98 - t55 * t70 + (-t222 * t4 - t29 * t324) * t168 + t287, -t29 * t353 + t12 * t146 + t121 * t15 + t162 * t26 - t260 * t2 + t36 * t98 + t55 * t255 + (t225 * t4 - t29 * t325) * t168 + t288, -t11 * t121 + t119 * t12 - t28 * t255 + t36 * t70 + (t21 * t225 - t222 * t26) * t161 + (-t1 * t225 + t2 * t222 + (t21 * t222 + t225 * t26) * qJD(4)) * t168 - t330, t2 * t36 + t26 * t12 + t1 * t28 + t21 * t11 + t4 * t55 + t29 * t15 - g(1) * (-t148 * pkin(5) + t262) - g(2) * (pkin(5) * t150 - qJ(6) * t347 + t254) + (-g(1) * (qJ(6) * t209 + t270) - t388) * t224; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t318, t319, -t299, -qJ(2) * t299 - t269, 0, 0, 0, 0, 0, 0, t279 + 0.2e1 * t410, t247 + 0.2e1 * t435, -t154 - t396, t105 * t159 - t106 * t433 + t175 - t415, 0, 0, 0, 0, 0, 0, t242, -t248, -t434, t159 * t289 + t425 * t225 + (t8 - t379) * t222 - t415, 0, 0, 0, 0, 0, 0, t242 -(t358 - t362) * t433 + t301 - t373, t248, -t159 * t42 + t426 * t225 + (t146 * t34 + t5) * t222 - t415, 0, 0, 0, 0, 0, 0, t242, t248, t434, t159 * t29 + t406 * t222 + t282 * t225 - t415; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t355, -t154 + t396, t247, t355, -t279, qJDD(3), qJD(3) * t106 - t159 * t176 - t292 - t386 + t416, -t176 * t433 + (t105 + t155) * qJD(3) + t311 + t313, 0, 0, t22, -t7, t25, t23, t429, -t357, -pkin(3) * t70 - t106 * t119 - t146 * t57 - t159 * t40 + t222 * t261 - t225 * t256 + t331, -pkin(3) * t255 - t106 * t121 + t146 * t58 + t159 * t41 + t222 * t256 + t225 * t261 - t332, t119 * t58 + t121 * t57 + (t433 * t40 + t8 + (-t40 + t391) * qJD(4)) * t225 + (t283 - t425) * t222 + t293, -t50 * pkin(3) - t41 * t58 - t40 * t57 + t289 * t106 - g(1) * (-pkin(3) * t347 + t188) - g(2) * (-pkin(3) * t349 + t185) - g(3) * t329 + (-qJD(4) * t277 + t222 * t302 + t8 * t225) * pkin(8), t22, t25, t7, -t357, -t429, t23, t369 * t119 + t146 * t39 + t159 * t34 + t405 * t222 + t257 * t225 - t271 * t70 + t331, t119 * t37 - t121 * t39 + (-t433 * t34 + t5 + (t34 + t391) * qJD(4)) * t225 + (t283 - t426) * t222 + t293, -t369 * t121 - t146 * t37 - t159 * t35 + t257 * t222 - t405 * t225 + t255 * t271 + t332, -t35 * t37 - t34 * t39 - g(1) * t188 - g(2) * t185 - g(3) * t291 + t369 * t42 + (-qJD(4) * t278 + t6 * t222 + t5 * t225) * pkin(8) + (t286 * t209 - t10) * t271, t22, t7, -t25, t23, t429, -t357, -t119 * t371 - t146 * t368 + t159 * t21 - t165 * t70 - t180 * t98 + t225 * t310 - t29 * t296 + t331, t121 * t371 + t146 * t370 - t159 * t26 + t165 * t255 + t181 * t98 + t222 * t310 + t29 * t294 + t332, t370 * t119 - t368 * t121 - t180 * t255 + t181 * t70 + t282 * t222 - t406 * t225 + t311, t2 * t181 + t1 * t180 + t4 * t165 - g(1) * (-qJ(6) * t344 + t188) - g(2) * (-qJ(6) * t346 + t185) - g(3) * (pkin(5) * t345 + t291) + t371 * t29 + t370 * t26 + t368 * t21 + (g(3) * qJ(6) + t286 * t165) * t209; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t361, t413, t33, -t361, t230, t98, t121 * t289 + t243 + t378, -t119 * t289 + t235, 0, 0, t361, t33, -t413, t98, -t230, -t361, -t119 * t72 - t237 + t378 + t94, -pkin(4) * t255 - t383 + (t35 - t41) * t121 + (t34 - t334) * t119, -t119 * t42 + t121 * t72 - t235 + t423, t5 * qJ(5) - t6 * pkin(4) - t42 * t72 - t34 * t41 - g(1) * t297 - g(2) * t298 - g(3) * (-pkin(4) * t350 + t182) + t334 * t35, t361, -t413, -t33, -t361, t230, t98, t119 * t56 + t146 * t31 + (pkin(5) + t395) * t98 + t404, t119 * t29 - t121 * t56 - t146 * t30 - t236 + t423 + t424, t383 + t395 * t255 + (-t26 + t31) * t121 + (-t21 + t335) * t119, t2 * qJ(5) - t1 * t395 - t21 * t31 - t29 * t56 - g(1) * (-pkin(5) * t149 + t297) - g(2) * (-pkin(5) * t147 + t298) - g(3) * (-t209 * t315 + t182) + t335 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t232, t33, t419, t237 - t380, 0, 0, 0, 0, 0, 0, t232, t419, -t33, -t381 - t394 - t404; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t231 - t360, t255 - t364, -t397 - t398, -t119 * t26 + t121 * t21 + t310 + t416;];
tau_reg  = t3;