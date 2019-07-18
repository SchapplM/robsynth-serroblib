% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RRRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 05:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RRRPPR7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR7_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 05:57:03
% EndTime: 2019-05-07 05:57:33
% DurationCPUTime: 12.44s
% Computational Cost: add. (51769->544), mult. (106264->716), div. (0->0), fcn. (74166->10), ass. (0->328)
t300 = cos(qJ(3));
t296 = sin(qJ(3));
t297 = sin(qJ(2));
t353 = qJD(1) * qJD(2);
t282 = t297 * t353;
t301 = cos(qJ(2));
t352 = t301 * qJDD(1);
t268 = -t282 + t352;
t261 = -qJDD(3) + t268;
t359 = qJD(1) * t297;
t262 = -qJD(2) * t300 + t296 * t359;
t264 = qJD(2) * t296 + t300 * t359;
t372 = t264 * t262;
t312 = t261 - t372;
t378 = t312 * t296;
t260 = t264 ^ 2;
t280 = qJD(1) * t301 - qJD(3);
t407 = t280 ^ 2;
t422 = -t260 - t407;
t174 = t300 * t422 + t378;
t377 = t312 * t300;
t176 = -t296 * t422 + t377;
t283 = t297 * qJDD(1);
t346 = t301 * t353;
t267 = t283 + t346;
t319 = t296 * qJDD(2) + t300 * t267;
t307 = (qJD(3) - t280) * t262 - t319;
t458 = pkin(7) * (t176 * t301 - t297 * t307) - pkin(1) * t174;
t456 = pkin(2) * t174;
t455 = pkin(8) * t174;
t454 = pkin(8) * t176;
t408 = t262 ^ 2;
t234 = t408 - t260;
t323 = -qJD(3) * t264 + qJDD(2) * t300 - t267 * t296;
t371 = t264 * t280;
t315 = -t323 - t371;
t383 = t307 * t296;
t452 = t301 * t234 + t297 * (-t300 * t315 + t383);
t242 = t408 - t407;
t314 = t323 - t371;
t450 = t297 * (t242 * t300 + t378) - t301 * t314;
t210 = t261 + t372;
t375 = t210 * t300;
t420 = -t407 - t408;
t431 = t296 * t420 - t375;
t449 = pkin(2) * t431;
t376 = t210 * t296;
t430 = t300 * t420 + t376;
t448 = pkin(8) * t430;
t447 = pkin(8) * t431;
t243 = -t260 + t407;
t445 = t243 * t300 - t376;
t292 = sin(pkin(10));
t293 = cos(pkin(10));
t228 = -t262 * t293 + t264 * t292;
t226 = t228 ^ 2;
t173 = -t407 - t226;
t230 = t262 * t292 + t264 * t293;
t184 = t230 * t228;
t424 = -t184 + t261;
t437 = t293 * t424;
t121 = t173 * t292 + t437;
t438 = t292 * t424;
t122 = t173 * t293 - t438;
t303 = qJD(1) ^ 2;
t298 = sin(qJ(1));
t302 = cos(qJ(1));
t344 = t298 * g(1) - g(2) * t302;
t250 = qJDD(1) * pkin(1) + t303 * pkin(7) + t344;
t324 = -t268 + t282;
t325 = t267 + t346;
t189 = pkin(2) * t324 - pkin(8) * t325 - t250;
t331 = g(1) * t302 + g(2) * t298;
t390 = qJDD(1) * pkin(7);
t251 = -pkin(1) * t303 - t331 + t390;
t334 = -pkin(2) * t301 - pkin(8) * t297;
t340 = t303 * t334 + t251;
t403 = t297 * g(3);
t406 = qJD(2) ^ 2;
t202 = -pkin(2) * t406 + qJDD(2) * pkin(8) + t301 * t340 - t403;
t150 = -t189 * t300 + t296 * t202;
t231 = pkin(3) * t262 - qJ(4) * t264;
t118 = pkin(3) * t261 - qJ(4) * t407 + t264 * t231 + qJDD(4) + t150;
t197 = (qJD(3) + t280) * t262 - t319;
t93 = pkin(4) * t210 + qJ(5) * t197 + t118;
t358 = qJD(4) * t280;
t272 = -0.2e1 * t358;
t151 = t189 * t296 + t202 * t300;
t339 = qJ(4) * t261 + t231 * t262 - t151;
t322 = t272 - t339;
t116 = -pkin(3) * t407 + t322;
t326 = pkin(4) * t280 - qJ(5) * t264;
t96 = -pkin(4) * t408 - qJ(5) * t323 - t280 * t326 + t116;
t336 = -0.2e1 * qJD(5) * t230 - t292 * t96 + t293 * t93;
t405 = pkin(3) + pkin(4);
t444 = qJ(4) * t122 - t121 * t405 - t336;
t443 = t242 * t296 - t377;
t441 = t296 * t315 + t300 * t307;
t440 = t297 * (-t243 * t296 - t375) + t301 * t197;
t439 = pkin(7) * (t297 * t315 + t301 * t430) - pkin(1) * t431;
t295 = sin(qJ(6));
t299 = cos(qJ(6));
t180 = t228 * t299 + t230 * t295;
t182 = -t228 * t295 + t230 * t299;
t128 = t182 * t180;
t254 = qJDD(6) + t261;
t425 = -t128 + t254;
t436 = t295 * t425;
t434 = t299 * t425;
t402 = t301 * g(3);
t201 = -qJDD(2) * pkin(2) - pkin(8) * t406 + t297 * t340 + t402;
t309 = -qJD(3) * t262 + t319;
t432 = t201 - t323 * pkin(3) + (-t262 * t280 - t309) * qJ(4);
t409 = -t197 * t296 + t300 * t314;
t421 = t260 + t408;
t429 = pkin(2) * t421 + pkin(8) * t409;
t428 = pkin(7) * (-t297 * t421 + t301 * t409);
t170 = -t292 * t323 + t293 * t309;
t374 = t228 * t280;
t143 = -t170 - t374;
t169 = -t292 * t309 - t293 * t323;
t206 = t280 * t230;
t423 = t206 + t169;
t418 = -pkin(3) * (t422 + t407) - qJ(4) * t312 - t339;
t43 = pkin(5) * t424 + pkin(9) * t143 + t336;
t333 = pkin(5) * t280 - pkin(9) * t230;
t357 = qJD(5) * t228;
t222 = -0.2e1 * t357;
t401 = t292 * t93 + t293 * t96;
t59 = t222 + t401;
t44 = -t226 * pkin(5) + t169 * pkin(9) - t280 * t333 + t59;
t25 = t295 * t44 - t299 * t43;
t26 = t295 * t43 + t299 * t44;
t13 = -t25 * t299 + t26 * t295;
t14 = t25 * t295 + t26 * t299;
t398 = t13 * t293;
t7 = t14 * t292 + t398;
t399 = t13 * t292;
t8 = t14 * t293 - t399;
t417 = -pkin(5) * t13 + qJ(4) * t8 - t405 * t7;
t179 = t182 ^ 2;
t275 = qJD(6) + t280;
t273 = t275 ^ 2;
t153 = -t179 - t273;
t124 = t128 + t254;
t389 = t124 * t295;
t97 = t153 * t299 - t389;
t388 = t124 * t299;
t98 = -t153 * t295 - t388;
t61 = t292 * t98 + t293 * t97;
t62 = -t292 * t97 + t293 * t98;
t416 = -pkin(5) * t97 + qJ(4) * t62 - t405 * t61 + t26;
t178 = t180 ^ 2;
t126 = -t273 - t178;
t81 = t126 * t295 + t434;
t82 = t126 * t299 - t436;
t51 = t292 * t82 + t293 * t81;
t52 = -t292 * t81 + t293 * t82;
t415 = -pkin(5) * t81 + qJ(4) * t52 - t405 * t51 + t25;
t342 = -t169 * t299 + t295 * t170;
t311 = (-qJD(6) + t275) * t182 - t342;
t347 = qJD(6) * t180 - t169 * t295 - t170 * t299;
t370 = t275 * t180;
t87 = t347 - t370;
t54 = t295 * t311 + t299 * t87;
t56 = -t295 * t87 + t299 * t311;
t32 = t292 * t56 + t293 * t54;
t34 = -t292 * t54 + t293 * t56;
t414 = -pkin(5) * t54 + qJ(4) * t34 - t32 * t405;
t35 = t292 * t59 + t293 * t336;
t36 = -t292 * t336 + t293 * t59;
t413 = qJ(4) * t36 - t35 * t405;
t104 = t143 * t293 + t292 * t423;
t106 = -t143 * t292 + t293 * t423;
t412 = qJ(4) * t106 - t104 * t405;
t227 = t230 ^ 2;
t200 = -t227 - t407;
t167 = t184 + t261;
t387 = t167 * t292;
t129 = t200 * t293 - t387;
t386 = t167 * t293;
t130 = -t200 * t292 - t386;
t410 = qJ(4) * t130 - t129 * t405 + t401;
t146 = t197 * t300 + t296 * t314;
t404 = pkin(3) * t300;
t117 = (-pkin(3) * t280 - 0.2e1 * qJD(4)) * t264 + t432;
t95 = -pkin(4) * t323 + qJ(5) * t408 - t264 * t326 - qJDD(5) + t117;
t397 = t292 * t95;
t396 = t293 * t95;
t70 = t169 * pkin(5) + t226 * pkin(9) - t230 * t333 + t95;
t395 = t295 * t70;
t394 = t299 * t70;
t391 = qJ(4) * t300;
t380 = t201 * t296;
t379 = t201 * t300;
t369 = t275 * t295;
t368 = t275 * t299;
t367 = t280 * t292;
t366 = t280 * t293;
t365 = t280 * t296;
t364 = t280 * t300;
t279 = t301 * t303 * t297;
t363 = t297 * (qJDD(2) + t279);
t361 = t301 * (-t279 + qJDD(2));
t351 = t262 * t364;
t350 = t301 * t128;
t349 = t301 * t184;
t348 = t301 * t372;
t345 = -qJ(4) * t296 - pkin(2);
t113 = t150 * t296 + t151 * t300;
t239 = t297 * t251 + t402;
t240 = t251 * t301 - t403;
t341 = t297 * t239 + t240 * t301;
t241 = t264 * t365;
t330 = t297 * (t300 * t309 + t241) - t348;
t329 = -t262 * t365 + t300 * t323;
t142 = -t170 + t374;
t328 = -pkin(3) * t118 + qJ(4) * t116;
t327 = pkin(3) * t197 + qJ(4) * t314;
t321 = t150 * t300 - t151 * t296;
t320 = -pkin(1) + t334;
t318 = -t347 - t370;
t316 = (t262 * t296 + t264 * t300) * t280;
t247 = t301 * t261;
t313 = t297 * (-t241 + t351) + t247;
t306 = 0.2e1 * qJD(4) * t264 - t432;
t305 = -pkin(3) * t210 + qJ(4) * t420 - t118;
t304 = t297 * (-t296 * t323 - t351) + t348;
t290 = t301 ^ 2;
t289 = t297 ^ 2;
t287 = t290 * t303;
t285 = t289 * t303;
t269 = -0.2e1 * t282 + t352;
t266 = t283 + 0.2e1 * t346;
t205 = -t227 + t407;
t204 = t226 - t407;
t186 = -t264 * t364 + t296 * t309;
t183 = t227 - t226;
t157 = -t179 + t273;
t156 = t178 - t273;
t155 = (-t228 * t293 + t230 * t292) * t280;
t154 = (t228 * t292 + t230 * t293) * t280;
t152 = -t226 - t227;
t139 = -t206 + t169;
t138 = t170 * t293 - t230 * t367;
t137 = -t170 * t292 - t230 * t366;
t136 = -t169 * t292 + t228 * t366;
t135 = -t169 * t293 - t228 * t367;
t134 = t204 * t293 - t387;
t133 = -t205 * t292 + t437;
t132 = -t204 * t292 - t386;
t131 = -t205 * t293 - t438;
t127 = t179 - t178;
t120 = (-t180 * t299 + t182 * t295) * t275;
t119 = (-t180 * t295 - t182 * t299) * t275;
t115 = qJ(4) * t421 + t118;
t114 = -t178 - t179;
t111 = (t421 - t407) * pkin(3) + t322;
t110 = (-t315 + t371) * pkin(3) + t306;
t109 = pkin(3) * t371 - qJ(4) * t307 + t306;
t107 = -qJD(6) * t182 - t342;
t105 = t139 * t293 + t142 * t292;
t103 = -t139 * t292 + t142 * t293;
t102 = t156 * t299 - t389;
t101 = -t157 * t295 + t434;
t100 = t156 * t295 + t388;
t99 = t157 * t299 + t436;
t90 = t129 * t296 + t130 * t300;
t89 = -t129 * t300 + t130 * t296;
t83 = (qJD(6) + t275) * t182 + t342;
t80 = -t182 * t369 - t299 * t347;
t79 = t182 * t368 - t295 * t347;
t78 = -t107 * t295 + t180 * t368;
t77 = t107 * t299 + t180 * t369;
t76 = t121 * t296 + t122 * t300;
t75 = -t121 * t300 + t122 * t296;
t74 = -t119 * t292 + t120 * t293;
t73 = -t119 * t293 - t120 * t292;
t72 = t116 * t300 + t118 * t296;
t71 = t116 * t296 - t118 * t300;
t69 = t104 * t296 + t106 * t300;
t68 = -t104 * t300 + t106 * t296;
t67 = -t100 * t292 + t102 * t293;
t66 = t101 * t293 - t292 * t99;
t65 = -t100 * t293 - t102 * t292;
t64 = -t101 * t292 - t293 * t99;
t63 = -qJ(4) * t142 - qJ(5) * t129 - t396;
t60 = -qJ(4) * t139 - qJ(5) * t121 - t397;
t57 = -qJ(5) * t130 - t142 * t405 + t397;
t55 = -t295 * t318 - t299 * t83;
t53 = -t295 * t83 + t299 * t318;
t50 = -t292 * t79 + t293 * t80;
t49 = -t292 * t77 + t293 * t78;
t48 = -t292 * t80 - t293 * t79;
t47 = -t292 * t78 - t293 * t77;
t46 = -qJ(5) * t122 - t139 * t405 - t396;
t45 = -pkin(9) * t97 - t394;
t41 = -pkin(9) * t81 - t395;
t40 = -pkin(5) * t318 + pkin(9) * t98 - t395;
t39 = t296 * t61 + t300 * t62;
t38 = t296 * t62 - t300 * t61;
t37 = -pkin(5) * t83 + pkin(9) * t82 + t394;
t33 = -t292 * t53 + t293 * t55;
t31 = -t292 * t55 - t293 * t53;
t30 = t296 * t51 + t300 * t52;
t29 = t296 * t52 - t300 * t51;
t28 = -qJ(4) * t95 - qJ(5) * t35;
t27 = qJ(4) * t152 - qJ(5) * t104 - t35;
t23 = -qJ(5) * t106 + t152 * t405 - t36;
t22 = -qJ(5) * t36 - t405 * t95;
t21 = t296 * t35 + t300 * t36;
t20 = t296 * t36 - t300 * t35;
t19 = t296 * t32 + t300 * t34;
t18 = t296 * t34 - t300 * t32;
t17 = qJ(4) * t318 - qJ(5) * t61 - t292 * t40 + t293 * t45;
t16 = qJ(4) * t83 - qJ(5) * t51 - t292 * t37 + t293 * t41;
t15 = -qJ(5) * t62 - t292 * t45 - t293 * t40 + t318 * t405;
t12 = -qJ(5) * t52 - t292 * t41 - t293 * t37 + t405 * t83;
t11 = pkin(5) * t70 + pkin(9) * t14;
t10 = -pkin(9) * t54 - t13;
t9 = -pkin(5) * t114 + pkin(9) * t56 + t14;
t6 = qJ(4) * t114 - qJ(5) * t32 + t10 * t293 - t292 * t9;
t5 = -qJ(5) * t34 - t292 * t10 + t114 * t405 - t293 * t9;
t4 = t296 * t7 + t300 * t8;
t3 = t296 * t8 - t300 * t7;
t2 = -pkin(9) * t398 - qJ(4) * t70 - qJ(5) * t7 - t11 * t292;
t1 = pkin(9) * t399 - qJ(5) * t8 - t293 * t11 - t405 * t70;
t24 = [0, 0, 0, 0, 0, qJDD(1), t344, t331, 0, 0, t325 * t297, t266 * t301 + t269 * t297, t363 + t301 * (-t285 + t406), -t324 * t301, t297 * (t287 - t406) + t361, 0, t301 * t250 + pkin(1) * t269 + pkin(7) * (t301 * (-t287 - t406) - t363), -t297 * t250 - pkin(1) * t266 + pkin(7) * (-t361 - t297 * (-t285 - t406)), pkin(1) * (t285 + t287) + (t289 + t290) * t390 + t341, pkin(1) * t250 + pkin(7) * t341, t330, t452, t440, t304, t450, t313, t297 * (t380 - t447) + t301 * (t150 - t449) + t439, t297 * (t379 - t455) + t301 * (t151 - t456) + t458, t320 * t146 + t297 * t321 + t428, pkin(7) * (t113 * t301 + t201 * t297) - t320 * t321, t330, t440, -t452, t313, -t450, t304, t297 * (-t110 * t296 - t315 * t391 - t447) + t301 * (-t305 - t449) + t439, t297 * (-pkin(8) * t146 - t111 * t296 + t115 * t300) + t301 * (-pkin(2) * t146 - t327) - pkin(1) * t146 + t428, t297 * (pkin(3) * t383 + t109 * t300 + t455) + t301 * (0.2e1 * t358 - t418 + t456) - t458, t297 * (-pkin(8) * t71 + (pkin(3) * t296 - t391) * t117) + t301 * (-pkin(2) * t71 - t328) - pkin(1) * t71 + pkin(7) * (t117 * t297 + t301 * t72), t297 * (-t137 * t296 + t138 * t300) + t349, t297 * (-t103 * t296 + t105 * t300) + t301 * t183, t297 * (-t131 * t296 + t133 * t300) - t301 * t143, t297 * (-t135 * t296 + t136 * t300) - t349, t297 * (-t132 * t296 + t134 * t300) + t301 * t423, t297 * (-t154 * t296 + t155 * t300) + t247, t297 * (-pkin(8) * t75 - t296 * t46 + t300 * t60) + t301 * (-pkin(2) * t75 - t444) - pkin(1) * t75 + pkin(7) * (t139 * t297 + t301 * t76), t297 * (-pkin(8) * t89 - t296 * t57 + t300 * t63) + t301 * (-pkin(2) * t89 + 0.2e1 * t357 - t410) - pkin(1) * t89 + pkin(7) * (t142 * t297 + t301 * t90), t297 * (-pkin(8) * t68 - t23 * t296 + t27 * t300) + t301 * (-pkin(2) * t68 - t412) - pkin(1) * t68 + pkin(7) * (-t152 * t297 + t301 * t69), t297 * (-pkin(8) * t20 - t22 * t296 + t28 * t300) + t301 * (-pkin(2) * t20 - t413) - pkin(1) * t20 + pkin(7) * (t21 * t301 + t297 * t95), t297 * (-t296 * t48 + t300 * t50) + t350, t297 * (-t296 * t31 + t300 * t33) + t301 * t127, t297 * (-t296 * t64 + t300 * t66) - t301 * t87, t297 * (-t296 * t47 + t300 * t49) - t350, t297 * (-t296 * t65 + t300 * t67) + t301 * t311, t297 * (-t296 * t73 + t300 * t74) + t301 * t254, t297 * (-pkin(8) * t29 - t12 * t296 + t16 * t300) + t301 * (-pkin(2) * t29 - t415) - pkin(1) * t29 + pkin(7) * (-t297 * t83 + t30 * t301), t297 * (-pkin(8) * t38 - t15 * t296 + t17 * t300) + t301 * (-pkin(2) * t38 - t416) - pkin(1) * t38 + pkin(7) * (-t297 * t318 + t301 * t39), t297 * (-pkin(8) * t18 - t296 * t5 + t300 * t6) + t301 * (-pkin(2) * t18 - t414) - pkin(1) * t18 + pkin(7) * (-t114 * t297 + t19 * t301), t297 * (-pkin(8) * t3 - t1 * t296 + t2 * t300) + t301 * (-pkin(2) * t3 - t417) - pkin(1) * t3 + pkin(7) * (t297 * t70 + t301 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t279, t285 - t287, t283, t279, t352, qJDD(2), -t239, -t240, 0, 0, t186, -t441, t445, t329, t443, t316, -pkin(2) * t315 - t379 + t448, pkin(2) * t307 + t380 + t454, t113 + t429, -pkin(2) * t201 + pkin(8) * t113, t186, t445, t441, t316, -t443, t329, t300 * t110 + t315 * t345 + t448, t111 * t300 + t115 * t296 + t429, -t454 + t296 * t109 + (-pkin(2) - t404) * t307, pkin(8) * t72 + (t345 - t404) * t117, t137 * t300 + t138 * t296, t103 * t300 + t105 * t296, t131 * t300 + t133 * t296, t135 * t300 + t136 * t296, t132 * t300 + t134 * t296, t154 * t300 + t155 * t296, -pkin(2) * t139 + pkin(8) * t76 + t296 * t60 + t300 * t46, -pkin(2) * t142 + pkin(8) * t90 + t296 * t63 + t300 * t57, pkin(2) * t152 + pkin(8) * t69 + t23 * t300 + t27 * t296, -pkin(2) * t95 + pkin(8) * t21 + t22 * t300 + t28 * t296, t296 * t50 + t300 * t48, t296 * t33 + t300 * t31, t296 * t66 + t300 * t64, t296 * t49 + t300 * t47, t296 * t67 + t300 * t65, t296 * t74 + t300 * t73, pkin(2) * t83 + pkin(8) * t30 + t12 * t300 + t16 * t296, pkin(2) * t318 + pkin(8) * t39 + t15 * t300 + t17 * t296, pkin(2) * t114 + pkin(8) * t19 + t296 * t6 + t300 * t5, -pkin(2) * t70 + pkin(8) * t4 + t1 * t300 + t2 * t296; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t372, -t234, -t197, -t372, t314, -t261, -t150, -t151, 0, 0, t372, -t197, t234, -t261, -t314, -t372, t305, t327, t272 + t418, t328, -t184, -t183, t143, t184, -t423, -t261, t444, t222 + t410, t412, t413, -t128, -t127, t87, t128, -t311, -t254, t415, t416, t414, t417; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t210, -t197, t422, t118, 0, 0, 0, 0, 0, 0, t121, t129, t104, t35, 0, 0, 0, 0, 0, 0, t51, t61, t32, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t139, -t142, t152, -t95, 0, 0, 0, 0, 0, 0, t83, t318, t114, -t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, t127, -t87, -t128, t311, t254, -t25, -t26, 0, 0;];
tauJ_reg  = t24;