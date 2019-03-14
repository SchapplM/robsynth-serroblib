% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% 
% Output:
% tauc_reg [6x35]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRRPR13_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR13_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:04:04
% EndTime: 2019-03-10 00:04:33
% DurationCPUTime: 11.64s
% Computational Cost: add. (11098->638), mult. (28682->863), div. (0->0), fcn. (22342->10), ass. (0->292)
t257 = sin(pkin(6));
t265 = cos(qJ(2));
t370 = qJD(1) * t265;
t244 = t257 * t370;
t307 = t244 - qJD(3);
t261 = sin(qJ(2));
t409 = cos(pkin(6));
t337 = t409 * qJD(1);
t321 = pkin(1) * t337;
t241 = t261 * t321;
t204 = pkin(8) * t244 + t241;
t260 = sin(qJ(3));
t264 = cos(qJ(3));
t458 = t204 + t307 * (pkin(3) * t260 - pkin(10) * t264);
t295 = t337 + qJD(2);
t371 = qJD(1) * t261;
t348 = t257 * t371;
t182 = t260 * t348 - t264 * t295;
t176 = qJD(4) + t182;
t259 = sin(qJ(4));
t372 = qJD(1) * t257;
t263 = cos(qJ(4));
t390 = t263 * t265;
t169 = (t259 * t261 + t264 * t390) * t372;
t365 = qJD(3) * t264;
t310 = t263 * t365 - t169;
t364 = qJD(4) * t259;
t342 = t260 * t364;
t457 = t310 - t342;
t323 = t260 * t244;
t367 = qJD(3) * t260;
t456 = t323 - t367;
t201 = -pkin(8) * t348 + t265 * t321;
t289 = (pkin(2) * t261 - pkin(9) * t265) * t257;
t202 = qJD(1) * t289;
t377 = t264 * t201 + t260 * t202;
t123 = pkin(10) * t348 + t377;
t233 = -pkin(3) * t264 - pkin(10) * t260 - pkin(2);
t391 = t263 * t264;
t250 = pkin(9) * t391;
t439 = -qJD(4) * t250 + t259 * t123 - t233 * t364 - t263 * t458;
t184 = t260 * t295 + t264 * t348;
t132 = t184 * t259 + t263 * t307;
t134 = t263 * t184 - t259 * t307;
t158 = pkin(9) * t295 + t204;
t198 = (-pkin(2) * t265 - pkin(9) * t261 - pkin(1)) * t257;
t174 = qJD(1) * t198;
t108 = -t260 * t158 + t264 * t174;
t94 = pkin(3) * t307 - t108;
t276 = t134 * qJ(5) - t94;
t430 = pkin(4) + pkin(5);
t29 = -t132 * t430 + t276;
t258 = sin(qJ(6));
t262 = cos(qJ(6));
t77 = t132 * t258 + t134 * t262;
t455 = t29 * t77;
t322 = t264 * t244;
t168 = t259 * t322 - t263 * t348;
t453 = -t259 * t365 + t168;
t356 = qJD(1) * qJD(2);
t340 = t257 * t356;
t319 = t265 * t340;
t296 = t264 * t319;
t268 = -qJD(3) * t182 + t296;
t452 = -qJD(4) * t307 + t268;
t363 = qJD(4) * t263;
t451 = t263 * t123 - t233 * t363 + t259 * t458;
t450 = qJD(6) - qJD(4);
t300 = -t262 * t132 + t134 * t258;
t449 = -t300 ^ 2 + t77 ^ 2;
t447 = t29 * t300;
t446 = t77 * t300;
t203 = qJD(2) * t289;
t193 = qJD(1) * t203;
t349 = pkin(1) * t409;
t397 = t257 * t261;
t440 = -pkin(8) * t397 + t265 * t349;
t205 = t440 * qJD(2);
t194 = qJD(1) * t205;
t281 = t158 * t367 - t174 * t365 - t260 * t193 - t264 * t194;
t320 = t261 * t340;
t54 = pkin(10) * t320 - t281;
t297 = t260 * t319;
t138 = t184 * qJD(3) + t297;
t303 = pkin(8) * t319;
t67 = t138 * pkin(3) - pkin(10) * t268 + qJD(2) * t241 + t303;
t157 = -pkin(2) * t295 - t201;
t91 = t182 * pkin(3) - t184 * pkin(10) + t157;
t109 = t264 * t158 + t260 * t174;
t95 = -pkin(10) * t307 + t109;
t338 = t259 * t54 - t263 * t67 + t95 * t363 + t91 * t364;
t428 = pkin(4) * t138;
t14 = t338 - t428;
t163 = t176 * qJ(5);
t41 = t259 * t91 + t263 * t95;
t35 = t163 + t41;
t445 = -t176 * t35 + t14;
t444 = -qJ(5) * t456 - t451;
t443 = t260 * t307;
t442 = t261 * t265;
t40 = -t259 * t95 + t263 * t91;
t388 = qJD(5) - t40;
t441 = qJD(5) * t259 + t109;
t438 = -t258 * t263 + t259 * t262;
t359 = qJD(6) * t262;
t360 = qJD(6) * t258;
t71 = t184 * t364 - t259 * t320 - t263 * t452;
t72 = t184 * t363 + t259 * t452 - t263 * t320;
t22 = t132 * t359 - t134 * t360 + t258 * t72 - t262 * t71;
t357 = qJD(6) - t176;
t436 = t300 * t357 + t22;
t39 = t132 * pkin(4) - t276;
t426 = pkin(10) * t138;
t435 = t176 * t39 - t426;
t406 = qJ(5) * t263;
t434 = t259 * t430 - t406;
t407 = qJ(5) * t259;
t433 = -t263 * t430 - t407;
t432 = t134 ^ 2;
t431 = t176 ^ 2;
t267 = qJD(1) ^ 2;
t429 = pkin(10) - pkin(11);
t427 = pkin(9) * t264;
t425 = pkin(11) * t182;
t424 = pkin(11) * t260;
t378 = -t260 * t201 + t264 * t202;
t283 = pkin(3) * t348 + t378;
t274 = qJ(5) * t169 + t283;
t279 = -pkin(9) - t434;
t361 = qJD(5) * t263;
t423 = t430 * t168 - t274 + (qJD(4) * t433 + t361) * t260 + t279 * t365;
t422 = t134 * t39;
t326 = -t158 * t365 - t174 * t367 + t264 * t193 - t260 * t194;
t272 = pkin(3) * t320 + t326;
t269 = -qJ(5) * t71 + qJD(5) * t134 + t272;
t19 = pkin(4) * t72 - t269;
t420 = t19 * t259;
t419 = t19 * t263;
t31 = pkin(11) * t132 + t41;
t28 = t163 + t31;
t418 = t258 * t28;
t417 = t259 * t71;
t416 = t272 * t259;
t415 = t272 * t263;
t366 = qJD(3) * t263;
t414 = -qJD(5) * t264 + (-t260 * t366 - t264 * t364) * pkin(9) + t444;
t351 = -pkin(9) * t259 - pkin(4);
t413 = pkin(4) * t323 + t351 * t367 - t439;
t308 = pkin(4) * t259 - t406;
t292 = pkin(9) + t308;
t309 = pkin(4) * t263 + t407;
t412 = -pkin(4) * t168 + (qJD(4) * t309 - t361) * t260 + t292 * t365 + t274;
t411 = -t176 * t434 + t441;
t410 = t176 * t308 - t441;
t408 = qJ(5) * t132;
t405 = t132 * t176;
t404 = t134 * t132;
t403 = t134 * t176;
t211 = t260 * t397 - t264 * t409;
t402 = t138 * t211;
t401 = t138 * t259;
t400 = t138 * t263;
t399 = t138 * t264;
t254 = t257 ^ 2;
t398 = t254 * t267;
t396 = t257 * t265;
t393 = t260 * t262;
t392 = t260 * t263;
t389 = -pkin(11) * t134 + t388;
t126 = pkin(3) * t184 + pkin(10) * t182;
t387 = t263 * t108 + t259 * t126;
t220 = t258 * t259 + t262 * t263;
t386 = t260 * t438 * t450 - t168 * t258 - t169 * t262 + t220 * t365;
t208 = t220 * t260;
t385 = qJD(6) * t208 + t258 * t457 + t453 * t262 - t363 * t393;
t196 = -pkin(2) * t409 - t440;
t212 = t260 * t409 + t264 * t397;
t118 = t211 * pkin(3) - t212 * pkin(10) + t196;
t278 = pkin(8) * t396 + t261 * t349;
t197 = pkin(9) * t409 + t278;
t379 = t264 * t197 + t260 * t198;
t120 = -pkin(10) * t396 + t379;
t384 = t259 * t118 + t263 * t120;
t382 = (t182 - t450) * t220;
t381 = -t182 * t438 + t258 * t363 + t259 * t359 - t262 * t364 - t263 * t360;
t380 = -t260 * t197 + t264 * t198;
t374 = t259 * t233 + t250;
t373 = t261 ^ 2 - t265 ^ 2;
t369 = qJD(2) * t261;
t368 = qJD(3) * t259;
t358 = t157 * qJD(3);
t25 = -t176 * t430 + t389;
t5 = pkin(11) * t71 - t138 * t430 + t338;
t136 = t138 * qJ(5);
t162 = t176 * qJD(5);
t287 = -t259 * t67 - t263 * t54 - t91 * t363 + t364 * t95;
t12 = t136 + t162 - t287;
t7 = pkin(11) * t72 + t12;
t355 = t25 * t359 + t258 * t5 + t262 * t7;
t354 = pkin(10) * t364;
t238 = t429 * t263;
t353 = t259 * t396;
t45 = t184 * qJ(5) + t387;
t48 = t211 * qJ(5) + t384;
t119 = pkin(3) * t396 - t380;
t350 = -t258 * t7 + t262 * t5;
t347 = t257 * t369;
t346 = qJD(2) * t396;
t343 = t176 * t364;
t341 = t254 * t356;
t339 = -t258 * t71 - t262 * t72;
t336 = t118 * t263 - t259 * t120;
t249 = t259 * t427;
t334 = t233 * t263 - t249;
t333 = t357 ^ 2;
t332 = t265 * t307;
t331 = t176 * t263;
t330 = t307 * t257;
t329 = qJD(3) * t307;
t327 = 0.2e1 * t341;
t325 = -t197 * t365 - t198 * t367 + t264 * t203 - t260 * t205;
t253 = t264 * pkin(4);
t137 = pkin(5) * t264 + t249 + t253 + (-t233 - t424) * t263;
t318 = -pkin(11) * t168 + qJD(6) * t137 + (-pkin(9) * qJD(3) + pkin(11) * qJD(4)) * t392 + (-qJD(5) + (-pkin(9) * qJD(4) + pkin(11) * qJD(3)) * t259) * t264 + t444;
t160 = -qJ(5) * t264 + t374;
t139 = t259 * t424 + t160;
t317 = qJD(6) * t139 - t430 * t323 - (-pkin(5) + t351) * t367 + t439 + (qJD(3) * t391 - t169 - t342) * pkin(11);
t316 = t257 * t267 * t409;
t315 = -0.2e1 * pkin(1) * t341;
t237 = t429 * t259;
t313 = qJD(6) * t237 + t259 * t425 - t429 * t364 - t45;
t99 = t259 * t108;
t312 = t99 + (-t126 + t425) * t263 - t430 * t184 + t450 * t238;
t9 = t258 * t25 + t262 * t28;
t150 = t212 * t263 - t353;
t33 = -pkin(11) * t150 - t211 * t430 - t336;
t149 = t212 * t259 + t257 * t390;
t36 = pkin(11) * t149 + t48;
t306 = -t258 * t36 + t262 * t33;
t305 = t258 * t33 + t262 * t36;
t34 = -pkin(4) * t176 + t388;
t304 = -t259 * t35 + t263 * t34;
t302 = qJ(5) * t262 - t258 * t430;
t301 = -qJ(5) * t258 - t262 * t430;
t299 = t149 * t262 - t150 * t258;
t97 = t149 * t258 + t150 * t262;
t294 = 0.2e1 * t337 + qJD(2);
t280 = -t197 * t367 + t198 * t365 + t260 * t203 + t264 * t205;
t63 = pkin(10) * t347 + t280;
t147 = qJD(3) * t212 + t260 * t346;
t148 = -qJD(3) * t211 + t264 * t346;
t206 = t278 * qJD(2);
t80 = t147 * pkin(3) - t148 * pkin(10) + t206;
t293 = -t118 * t364 - t120 * t363 - t259 * t63 + t263 * t80;
t291 = qJ(5) * t150 - t119;
t1 = -t28 * t360 + t355;
t290 = t176 * t41 - t338;
t284 = t176 * t94 - t426;
t282 = t118 * t363 - t120 * t364 + t259 * t80 + t263 * t63;
t277 = pkin(3) * t347 + t325;
t18 = t147 * qJ(5) + t211 * qJD(5) + t282;
t275 = t176 * t40 + t287;
t273 = pkin(1) * (-qJD(2) * t337 + t398);
t2 = -qJD(6) * t9 + t350;
t23 = qJD(6) * t77 + t339;
t88 = -qJD(4) * t149 + t148 * t263 + t259 * t347;
t270 = qJ(5) * t88 + qJD(5) * t150 + t277;
t229 = -pkin(3) - t309;
t216 = pkin(3) - t433;
t207 = t258 * t392 - t259 * t393;
t199 = t292 * t260;
t195 = qJD(1) * t206;
t161 = t253 - t334;
t159 = t279 * t260;
t87 = -qJD(4) * t353 + t148 * t259 + t212 * t363 - t263 * t347;
t79 = pkin(4) * t134 + t408;
t56 = pkin(4) * t149 - t291;
t51 = -t134 * t430 - t408;
t49 = -pkin(4) * t211 - t336;
t47 = -pkin(4) * t184 - t126 * t263 + t99;
t43 = -t149 * t430 + t291;
t38 = -t71 + t405;
t27 = qJD(6) * t299 + t258 * t87 + t262 * t88;
t26 = qJD(6) * t97 + t258 * t88 - t87 * t262;
t21 = pkin(4) * t87 - t270;
t20 = -pkin(4) * t147 - t293;
t17 = -t430 * t87 + t270;
t13 = pkin(11) * t87 + t18;
t11 = -pkin(11) * t88 - t147 * t430 - t293;
t10 = -t430 * t72 + t269;
t8 = t25 * t262 - t418;
t3 = [0, 0, 0, t327 * t442, -t373 * t327, t294 * t346, -t294 * t347, 0, -t195 * t409 - t206 * t295 + t261 * t315, -t194 * t409 - t205 * t295 + t265 * t315, t184 * t148 + t212 * t268, -t212 * t138 - t184 * t147 - t148 * t182 - t211 * t268, -t148 * t307 + t184 * t347 + t212 * t320 - t268 * t396, t147 * t307 + (t138 * t265 + (-qJD(1) * t211 - t182) * t369) * t257 (-t254 * t370 - t330) * t369, -t325 * t307 + t206 * t182 + t196 * t138 + t195 * t211 + t157 * t147 + (-t326 * t265 + (qJD(1) * t380 + t108) * t369) * t257, -t109 * t347 + t157 * t148 + t206 * t184 + t195 * t212 + t196 * t268 + t280 * t307 - t281 * t396 - t320 * t379, t134 * t88 - t150 * t71, -t132 * t88 - t134 * t87 + t149 * t71 - t150 * t72, t134 * t147 + t138 * t150 + t176 * t88 - t211 * t71, -t132 * t147 - t138 * t149 - t176 * t87 - t211 * t72, t147 * t176 + t402, t119 * t72 - t132 * t277 + t138 * t336 + t40 * t147 - t149 * t272 + t176 * t293 - t211 * t338 + t94 * t87, -t119 * t71 - t134 * t277 - t138 * t384 - t41 * t147 - t150 * t272 - t176 * t282 + t211 * t287 + t94 * t88, t132 * t21 - t138 * t49 - t14 * t211 - t147 * t34 + t149 * t19 - t176 * t20 + t39 * t87 + t56 * t72, -t12 * t149 - t132 * t18 + t134 * t20 + t14 * t150 + t34 * t88 - t35 * t87 - t48 * t72 - t49 * t71, t12 * t211 - t134 * t21 + t138 * t48 + t147 * t35 - t150 * t19 + t176 * t18 - t39 * t88 + t56 * t71, t12 * t48 + t14 * t49 + t18 * t35 + t19 * t56 + t20 * t34 + t21 * t39, t22 * t97 + t27 * t77, t22 * t299 - t23 * t97 - t26 * t77 - t27 * t300, -t138 * t97 - t147 * t77 - t211 * t22 + t27 * t357, -t138 * t299 + t147 * t300 + t211 * t23 - t26 * t357, -t147 * t357 + t402 (-qJD(6) * t305 + t11 * t262 - t13 * t258) * t357 - t306 * t138 - t2 * t211 - t8 * t147 + t17 * t300 + t43 * t23 - t10 * t299 + t29 * t26 -(qJD(6) * t306 + t11 * t258 + t13 * t262) * t357 + t305 * t138 + t1 * t211 + t9 * t147 + t17 * t77 + t43 * t22 + t10 * t97 + t29 * t27; 0, 0, 0, -t398 * t442, t373 * t398, -t265 * t316, t261 * t316, 0, t204 * t295 + t261 * t273 - t303, pkin(8) * t320 + t201 * t295 + t265 * t273, -qJD(3) * t260 ^ 2 * t348 + ((qJD(3) * t295 + t319) * t260 - t307 * t184) * t264, -t260 * t138 + t264 * t268 + t456 * t184 + (t322 - t365) * t182, -t264 * t329 + (t264 * t332 + (t260 * qJD(2) - t184) * t261) * t372, t260 * t329 + (-t260 * t332 + (t264 * qJD(2) + t182) * t261) * t372, t330 * t371, -pkin(2) * t138 + t260 * t358 + t378 * t307 - t204 * t182 + (pkin(9) * t329 - t195) * t264 + (-t108 * t261 + (-pkin(9) * t369 - t157 * t265) * t260) * t372, -pkin(2) * t268 + t109 * t348 - t157 * t322 - t204 * t184 + t195 * t260 + t264 * t358 - t320 * t427 + (-pkin(9) * t367 - t377) * t307, t134 * t457 - t71 * t392, t132 * t169 + t134 * t168 + (-t132 * t263 - t134 * t259) * t365 + (t417 - t263 * t72 + (t132 * t259 - t134 * t263) * qJD(4)) * t260, t264 * t71 + t310 * t176 + (-t134 * t307 - t343 + t400) * t260, t264 * t72 + t453 * t176 + (t132 * t307 - t176 * t363 - t401) * t260, -t176 * t443 - t399, t334 * t138 + t283 * t132 - t94 * t168 + t439 * t176 + (t338 + (pkin(9) * t132 + t259 * t94) * qJD(3)) * t264 + (t94 * t363 - t416 - t307 * t40 + (t176 * t368 + t72) * pkin(9)) * t260, -t374 * t138 + t283 * t134 - t94 * t169 + t451 * t176 + (t94 * t366 - t287 + (qJD(3) * t134 + t343) * pkin(9)) * t264 + (-t94 * t364 - t415 + t307 * t41 + (t176 * t366 - t71) * pkin(9)) * t260, -t138 * t161 - t168 * t39 + t199 * t72 + (t368 * t39 + t14) * t264 - t413 * t176 + t412 * t132 + (t307 * t34 + t363 * t39 + t420) * t260, -t160 * t72 - t161 * t71 + t168 * t35 - t169 * t34 + t413 * t134 - t414 * t132 + t304 * t365 + (-t12 * t259 + t14 * t263 + (-t259 * t34 - t263 * t35) * qJD(4)) * t260, t138 * t160 + t169 * t39 + t199 * t71 + (-t366 * t39 - t12) * t264 + t414 * t176 - t412 * t134 + (-t307 * t35 + t364 * t39 - t419) * t260, t12 * t160 + t14 * t161 + t19 * t199 + t34 * t413 + t35 * t414 + t39 * t412, t208 * t22 + t386 * t77, -t207 * t22 - t208 * t23 - t300 * t386 - t385 * t77, -t138 * t208 + t22 * t264 + t357 * t386 + t443 * t77, t138 * t207 - t23 * t264 - t300 * t443 - t357 * t385, t357 * t443 - t399 -(t137 * t262 - t139 * t258) * t138 + t2 * t264 + t159 * t23 + t10 * t207 + t423 * t300 + t385 * t29 + t8 * t443 - (t258 * t318 + t262 * t317) * t357 (t137 * t258 + t139 * t262) * t138 - t1 * t264 + t159 * t22 + t10 * t208 + t423 * t77 + t386 * t29 - t9 * t443 - (-t258 * t317 + t262 * t318) * t357; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t184 * t182, -t182 ^ 2 + t184 ^ 2, -t182 * t244 + t296, -t184 * t244 - t297, t320, -t109 * t307 - t157 * t184 + t326, -t108 * t307 + t157 * t182 + t281, t134 * t331 - t417 (-t71 - t405) * t263 + (-t72 - t403) * t259, -t134 * t184 + t176 * t331 + t401, t132 * t184 - t259 * t431 + t400, -t176 * t184, -pkin(3) * t72 - t109 * t132 - t40 * t184 + t415 + (t99 + (-pkin(10) * qJD(4) - t126) * t263) * t176 + t284 * t259, pkin(3) * t71 - t109 * t134 + t41 * t184 - t416 + (t354 + t387) * t176 + t284 * t263, t184 * t34 - t419 + t229 * t72 + (-pkin(10) * t363 + t47) * t176 + t410 * t132 + t435 * t259, t132 * t45 - t134 * t47 + (t12 + t176 * t34 + (qJD(4) * t134 - t72) * pkin(10)) * t263 + ((qJD(4) * t132 - t71) * pkin(10) + t445) * t259, -t184 * t35 - t420 + t229 * t71 + (-t45 - t354) * t176 - t410 * t134 - t435 * t263, t19 * t229 - t34 * t47 - t35 * t45 + t410 * t39 + (qJD(4) * t304 + t12 * t263 + t14 * t259) * pkin(10), t22 * t438 + t382 * t77, -t22 * t220 - t23 * t438 - t300 * t382 - t381 * t77, -t138 * t438 + t184 * t77 + t357 * t382, t138 * t220 - t184 * t300 - t357 * t381, t357 * t184 -(t237 * t262 - t238 * t258) * t138 + t216 * t23 + t10 * t220 + t8 * t184 + t411 * t300 + t381 * t29 - (t258 * t313 + t262 * t312) * t357 (t237 * t258 + t238 * t262) * t138 + t216 * t22 + t10 * t438 - t9 * t184 + t411 * t77 + t382 * t29 - (-t258 * t312 + t262 * t313) * t357; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t404, -t132 ^ 2 + t432, t38, -t72 + t403, t138, -t134 * t94 + t290, t132 * t94 + t275, -t132 * t79 + t290 - t422 + 0.2e1 * t428, pkin(4) * t71 - qJ(5) * t72 + (t35 - t41) * t134 + (t34 - t388) * t132, -t132 * t39 + t134 * t79 + 0.2e1 * t136 + 0.2e1 * t162 - t275, -pkin(4) * t14 + qJ(5) * t12 - t34 * t41 + t35 * t388 - t39 * t79, -t446, -t449, -t436, -t357 * t77 + t23, t138, -t301 * t138 - t51 * t300 + t455 - (t258 * t389 + t262 * t31) * t357 + (-t302 * t357 + t9) * qJD(6) - t350, t302 * t138 - t51 * t77 - t447 - (-t258 * t31 + t262 * t389) * t357 + (-t301 * t357 - t418) * qJD(6) + t355; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t138 + t404, t38, -t431 - t432, t422 + t445, 0, 0, 0, 0, 0, -t134 * t300 - t138 * t262 - t258 * t333, -t134 * t77 + t138 * t258 - t262 * t333; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t446, t449, t436, -t339 + (-qJD(6) + t357) * t77, -t138, t357 * t9 + t2 - t455, t357 * t8 - t1 + t447;];
tauc_reg  = t3;