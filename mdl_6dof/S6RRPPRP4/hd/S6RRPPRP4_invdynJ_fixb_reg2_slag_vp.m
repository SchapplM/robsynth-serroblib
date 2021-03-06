% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRPPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPPRP4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP4_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:39:49
% EndTime: 2019-03-09 08:40:05
% DurationCPUTime: 8.45s
% Computational Cost: add. (7164->649), mult. (16267->783), div. (0->0), fcn. (11288->8), ass. (0->311)
t252 = cos(pkin(9));
t234 = t252 * qJDD(2);
t251 = sin(pkin(9));
t257 = cos(qJ(2));
t374 = qJD(1) * qJD(2);
t356 = t257 * t374;
t255 = sin(qJ(2));
t373 = t255 * qJDD(1);
t294 = t356 + t373;
t131 = t251 * t294 - t234;
t132 = t251 * qJDD(2) + t252 * t294;
t254 = sin(qJ(5));
t433 = cos(qJ(5));
t384 = qJD(1) * t255;
t365 = t251 * t384;
t376 = t252 * qJD(2);
t174 = -t365 + t376;
t363 = t252 * t384;
t382 = qJD(2) * t251;
t175 = t363 + t382;
t98 = -t174 * t254 + t175 * t433;
t35 = qJD(5) * t98 - t131 * t433 + t254 * t132;
t383 = qJD(1) * t257;
t217 = qJD(5) + t383;
t415 = t98 * t217;
t459 = t35 - t415;
t355 = t251 * t373;
t248 = t255 ^ 2;
t249 = t257 ^ 2;
t386 = t248 - t249;
t444 = qJD(1) * t386;
t460 = qJD(2) * (-t174 * t255 + t251 * t444) - (t131 + t355) * t257;
t437 = t98 ^ 2;
t358 = qJD(5) * t433;
t377 = qJD(5) * t254;
t159 = t251 * t358 - t252 * t377;
t340 = t433 * t383;
t362 = t252 * t383;
t393 = -t251 * t340 + t254 * t362 - t159;
t300 = -t251 * t254 - t252 * t433;
t147 = t300 * t257;
t158 = t300 * qJD(5);
t392 = -qJD(1) * t147 - t158;
t256 = sin(qJ(1));
t258 = cos(qJ(1));
t329 = g(1) * t258 + g(2) * t256;
t305 = t329 * t255;
t229 = pkin(7) * t373;
t149 = -qJDD(2) * pkin(2) + pkin(7) * t356 + qJDD(3) + t229;
t247 = g(3) * t257;
t351 = -t149 - t247;
t271 = -t305 - t351;
t281 = t340 + t358;
t352 = -qJ(4) * t251 - pkin(2);
t458 = pkin(3) * t252 - t352;
t448 = (t174 * t252 - t175 * t251) * t257;
t456 = qJD(1) * t448 + t251 * t131 - t132 * t252;
t412 = t132 * t251;
t413 = t131 * t252;
t455 = qJD(2) * t448 - t255 * (t412 + t413);
t454 = -2 * pkin(1);
t302 = t174 * t433 + t175 * t254;
t453 = t302 ^ 2;
t388 = pkin(2) * t257 + qJ(3) * t255;
t452 = -pkin(1) - t388;
t425 = -pkin(8) + qJ(3);
t191 = t425 * t251;
t192 = t425 * t252;
t112 = t191 * t254 + t192 * t433;
t367 = t433 * t251;
t379 = qJD(3) * t252;
t369 = -pkin(7) * t251 - pkin(3);
t403 = t252 * t257;
t280 = -pkin(8) * t403 + (-pkin(4) + t369) * t255;
t324 = pkin(2) * t255 - qJ(3) * t257;
t183 = t324 * qJD(1);
t405 = t252 * t183;
t67 = qJD(1) * t280 - t405;
t164 = t251 * t183;
t227 = qJ(4) * t384;
t404 = t252 * t255;
t407 = t251 * t257;
t298 = -pkin(7) * t404 + pkin(8) * t407;
t81 = qJD(1) * t298 + t164 + t227;
t421 = qJD(3) * t367 - qJD(5) * t112 - t433 * t67 + (-t379 + t81) * t254;
t451 = t217 * t302;
t402 = t254 * t252;
t180 = t367 - t402;
t287 = t255 * t180;
t450 = t255 * t300;
t364 = t251 * t383;
t378 = qJD(4) * t251;
t232 = pkin(7) * t383;
t391 = qJ(4) * t362 - t232;
t435 = -pkin(3) - pkin(4);
t346 = -t364 * t435 + t378 - t391;
t172 = t175 ^ 2;
t446 = -t174 ^ 2 - t172;
t381 = qJD(2) * t255;
t445 = qJ(4) * t381 - t257 * qJD(4);
t237 = t257 * qJDD(1);
t357 = t255 * t374;
t293 = -t357 + t237;
t166 = t175 * qJD(4);
t443 = pkin(3) * t131 - qJ(4) * t132 - t166;
t34 = -t131 * t254 - t132 * t433 + t174 * t358 + t175 * t377;
t442 = t35 * pkin(5) + t34 * qJ(6) - t98 * qJD(6);
t178 = -qJDD(5) - t293;
t441 = t175 * t98 - t178 * t254 + t217 * t281;
t440 = -t178 * t300 + t217 * t393 - t302 * t384;
t353 = t252 * t237;
t260 = qJD(1) ^ 2;
t409 = t249 * t260;
t439 = -t251 * t409 + (-t174 + t376) * t384 - t353;
t134 = pkin(7) * t403 + t251 * t452;
t114 = -qJ(4) * t257 + t134;
t408 = t251 * t255;
t100 = pkin(8) * t408 + t114;
t209 = pkin(7) * t407;
t242 = t257 * pkin(3);
t90 = t257 * pkin(4) + t209 + t242 + (-pkin(8) * t255 - t452) * t252;
t419 = t100 * t433 + t254 * t90;
t156 = qJD(2) * t324 - qJD(3) * t255;
t406 = t252 * t156;
t60 = qJD(2) * t280 - t406;
t137 = t251 * t156;
t61 = qJD(2) * t298 + t137 + t445;
t11 = -qJD(5) * t419 - t254 * t61 + t433 * t60;
t438 = -t180 * t35 - t300 * t34 + t302 * t392 + t393 * t98;
t436 = t217 ^ 2;
t432 = pkin(1) * t260;
t431 = pkin(5) * t178;
t430 = pkin(7) * t175;
t429 = g(1) * t256;
t426 = t98 * t302;
t424 = -pkin(5) * t393 + qJ(6) * t392 - qJD(6) * t180 + t346;
t33 = t254 * t67 + t433 * t81;
t29 = -qJ(6) * t384 + t33;
t301 = t191 * t433 - t192 * t254;
t64 = -qJD(3) * t300 + qJD(5) * t301;
t423 = t64 - t29;
t422 = t64 - t33;
t420 = -pkin(5) * t384 - t421;
t168 = t452 * qJD(1);
t195 = qJD(2) * qJ(3) + t232;
t101 = t168 * t252 - t195 * t251;
t79 = pkin(3) * t383 + qJD(4) - t101;
t52 = pkin(4) * t383 - pkin(8) * t175 + t79;
t102 = t168 * t251 + t195 * t252;
t89 = -qJ(4) * t383 + t102;
t56 = -pkin(8) * t174 + t89;
t18 = -t254 * t56 + t433 * t52;
t418 = t18 * t217;
t19 = t254 * t52 + t433 * t56;
t417 = t19 * t217;
t138 = pkin(7) * t293 + qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t92 = qJD(1) * t156 + qJDD(1) * t452;
t51 = t138 * t252 + t251 * t92;
t414 = pkin(7) * qJDD(1);
t410 = t178 * qJ(6);
t401 = t255 * t256;
t400 = t255 * t258;
t399 = t256 * t252;
t398 = t256 * t257;
t397 = t257 * t258;
t396 = t258 * t251;
t394 = qJD(6) - t18;
t361 = t257 * t376;
t390 = -qJ(4) * t361 - qJD(4) * t404;
t389 = g(1) * t401 - g(2) * t400;
t387 = pkin(1) * t258 + pkin(7) * t256;
t385 = t248 + t249;
t380 = qJD(2) * t257;
t372 = pkin(7) * t381;
t371 = -t437 + t453;
t370 = g(1) * t397 + g(2) * t398 + g(3) * t255;
t368 = pkin(3) * t251 + pkin(7);
t366 = qJ(3) * t381;
t150 = t174 * t383;
t186 = -qJD(2) * pkin(2) + pkin(7) * t384 + qJD(3);
t37 = t149 + t443;
t359 = -t37 - t247;
t354 = t251 * t237;
t50 = -t138 * t251 + t252 * t92;
t307 = pkin(3) * t237 + qJDD(4) - t50;
t41 = -pkin(3) * t357 + t307;
t24 = pkin(4) * t293 - t132 * pkin(8) + t41;
t36 = qJ(4) * t357 + (-qJ(4) * qJDD(1) - qJD(1) * qJD(4)) * t257 + t51;
t27 = t131 * pkin(8) + t36;
t4 = t24 * t433 - t254 * t27 - t358 * t56 - t377 * t52;
t133 = t252 * t452 - t209;
t165 = pkin(4) * t252 + t458;
t345 = g(1) * t252 * t400 + g(2) * t255 * t399 + qJ(3) * t354 + qJD(3) * t364;
t344 = pkin(3) * t403 + qJ(4) * t407 + t388;
t343 = pkin(2) * t397 + qJ(3) * t400 + t387;
t342 = t257 * t367;
t339 = t255 * t356;
t338 = qJ(3) * t353;
t337 = t251 * t435 - pkin(7);
t210 = qJ(3) * t398;
t336 = -pkin(8) * t398 + t210;
t215 = qJ(3) * t397;
t335 = -pkin(8) * t397 + t215;
t160 = t251 * t398 + t252 * t258;
t161 = t252 * t398 - t396;
t82 = -t160 * t433 + t161 * t254;
t162 = t257 * t396 - t399;
t163 = t251 * t256 + t252 * t397;
t86 = -t162 * t433 + t163 * t254;
t334 = g(1) * t82 - g(2) * t86;
t83 = t160 * t254 + t161 * t433;
t87 = t162 * t254 + t163 * t433;
t333 = g(1) * t83 - g(2) * t87;
t332 = t255 * t369;
t331 = -g(1) * t160 + g(2) * t162;
t330 = g(1) * t161 - g(2) * t163;
t328 = -g(2) * t258 + t429;
t327 = t453 + t437;
t244 = t258 * pkin(7);
t326 = -pkin(3) * t161 - t160 * qJ(4) + t244;
t71 = -qJD(2) * t342 - t158 * t255 + t254 * t361;
t325 = -t287 * t35 + t302 * t71;
t323 = pkin(5) * t287 - qJ(6) * t450;
t322 = -pkin(7) * t174 + t186 * t251;
t318 = -qJ(3) * t413 + t174 * t379 - t370;
t317 = pkin(4) * t403 + t344;
t74 = -pkin(3) * t174 - qJ(4) * t175 + t186;
t312 = qJ(3) * t132 + qJD(3) * t175;
t311 = qJD(1) * (-t175 + t382);
t259 = qJD(2) ^ 2;
t310 = qJDD(2) * t257 - t255 * t259;
t120 = -pkin(7) * t363 + t164;
t108 = -t252 * t372 + t137;
t43 = -t100 * t254 + t433 * t90;
t3 = t24 * t254 + t27 * t433 + t358 * t52 - t377 * t56;
t10 = -t100 * t377 + t254 * t60 + t358 * t90 + t433 * t61;
t299 = -t300 * t35 - t302 * t393;
t297 = -pkin(4) * t161 + pkin(8) * t401 + t326;
t291 = -t112 * t35 + t301 * t34 - t302 * t64 + t370;
t290 = t150 * t251 - t413;
t205 = qJ(4) * t404;
t113 = t255 * t337 + t205;
t55 = pkin(4) * t174 - t74;
t288 = pkin(3) * t163 + qJ(4) * t162 + t343;
t285 = (t131 * t255 - t174 * t380) * t251;
t284 = t452 * t429;
t283 = -g(1) * t162 - g(2) * t160 - g(3) * t408;
t122 = t131 * pkin(4);
t282 = t122 + t37;
t279 = t34 + t451;
t128 = t256 * t450;
t130 = t258 * t450;
t278 = -g(1) * t130 - g(2) * t128 + g(3) * t147 - t178 * t301;
t127 = t256 * t287;
t129 = t258 * t287;
t146 = t257 * t402 - t342;
t277 = -g(1) * t129 - g(2) * t127 - g(3) * t146 - t112 * t178;
t72 = -qJD(2) * t147 + t159 * t255;
t276 = -t287 * t34 - t302 * t72 + t35 * t450 - t71 * t98;
t78 = t337 * t380 - t390;
t275 = g(1) * t86 + g(2) * t82 - g(3) * t287 + t4;
t274 = -t178 * t287 - t217 * t71 - t257 * t35 + t302 * t381;
t273 = t257 * t311 - t234 + t355;
t272 = pkin(4) * t163 - pkin(8) * t400 + t288;
t269 = -g(1) * t87 - g(2) * t83 + g(3) * t450 + t3;
t268 = -t254 * t459 - t281 * t302 + t34 * t433;
t267 = -t175 * t302 - t178 * t433 - t254 * t436;
t17 = pkin(5) * t302 - qJ(6) * t98 + t55;
t266 = t17 * t98 + qJDD(6) - t275;
t265 = -t35 - t415;
t264 = t271 + t443;
t263 = t122 + t264;
t262 = (g(3) * pkin(8) + t329 * (-t252 * t435 - t352)) * t255;
t208 = t255 * t260 * t257;
t184 = t217 * t384;
t170 = qJDD(1) * t249 - 0.2e1 * t339;
t139 = t255 * t368 - t205;
t126 = pkin(3) * t364 - t391;
t119 = pkin(7) * t365 + t405;
t118 = -t133 + t242;
t109 = -t178 * t257 - t217 * t381;
t107 = t251 * t372 + t406;
t106 = qJD(1) * t332 - t405;
t105 = t120 + t227;
t104 = t368 * t380 + t390;
t93 = qJD(2) * t332 - t406;
t77 = t108 + t445;
t75 = -t150 + t132;
t68 = t252 * t409 + t255 * t311 - t354;
t66 = -pkin(5) * t300 - qJ(6) * t180 + t165;
t62 = -t175 * t362 + t412;
t59 = (t132 * t255 + t175 * t380) * t252;
t47 = t113 - t323;
t46 = pkin(5) * t98 + qJ(6) * t302;
t45 = (-t252 * t373 - t132) * t257 + (t175 * t255 + t252 * t444) * qJD(2);
t39 = -pkin(5) * t257 - t43;
t38 = qJ(6) * t257 + t419;
t22 = -t178 * t180 - t217 * t392 + t384 * t98;
t16 = qJ(6) * t217 + t19;
t15 = -t34 + t451;
t14 = -pkin(5) * t217 + t394;
t13 = t71 * pkin(5) - t72 * qJ(6) + qJD(6) * t450 + t78;
t12 = t34 * t450 + t72 * t98;
t9 = pkin(5) * t381 - t11;
t8 = -qJ(6) * t381 + qJD(6) * t257 + t10;
t7 = -t34 * t180 - t392 * t98;
t6 = t178 * t450 + t217 * t72 - t257 * t34 - t381 * t98;
t5 = -t282 + t442;
t2 = qJDD(6) - t4 + t431;
t1 = qJD(6) * t217 + t3 - t410;
t20 = [0, 0, 0, 0, 0, qJDD(1), t328, t329, 0, 0, qJDD(1) * t248 + 0.2e1 * t339, -0.2e1 * qJD(2) * t444 + 0.2e1 * t237 * t255, qJDD(2) * t255 + t257 * t259, t170, t310, 0 (-pkin(7) * qJDD(2) + t374 * t454) * t255 + (0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t259 + t328) * t257, -pkin(7) * t310 + t294 * t454 - t389, 0.2e1 * t385 * t414 - t329, -g(1) * (-pkin(1) * t256 + t244) - g(2) * t387 + (pkin(7) ^ 2 * t385 + (pkin(1) ^ 2)) * qJDD(1), t59, t455, t45, t285, -t460, t170 (pkin(7) * t131 + t149 * t251 + (qJD(1) * t133 + t101) * qJD(2)) * t255 + (-qJD(1) * t107 + qJD(2) * t322 - qJDD(1) * t133 - t50) * t257 + t330 (pkin(7) * t132 + t149 * t252 + (-qJD(1) * t134 - t102) * qJD(2)) * t255 + (qJD(1) * t108 + qJDD(1) * t134 + t51 + (t186 * t252 + t430) * qJD(2)) * t257 + t331, -t107 * t175 + t108 * t174 - t134 * t131 - t133 * t132 + (-t251 * t51 - t252 * t50) * t255 + (-t101 * t252 - t102 * t251) * t380 + t389, t51 * t134 + t102 * t108 + t50 * t133 + t101 * t107 - g(1) * t244 - g(2) * t343 - t284 + (t149 * t255 + t186 * t380) * pkin(7), t59, t45, -t455, t170, t460, t285, -t104 * t174 + t139 * t131 + (t37 * t251 + (-qJD(1) * t118 - t79) * qJD(2)) * t255 + (qJD(1) * t93 + qJDD(1) * t118 + t382 * t74 + t41) * t257 + t330, -t114 * t131 + t118 * t132 + t77 * t174 + t93 * t175 + (-t251 * t36 + t252 * t41) * t255 + (-t251 * t89 + t252 * t79) * t380 + t389, -t104 * t175 - t139 * t132 + (-t37 * t252 + (qJD(1) * t114 + t89) * qJD(2)) * t255 + (-qJD(1) * t77 - qJDD(1) * t114 - t376 * t74 - t36) * t257 - t331, -g(1) * t326 - g(2) * t288 + t74 * t104 + t36 * t114 + t41 * t118 + t37 * t139 + t89 * t77 + t79 * t93 - t284, t12, t276, t6, t325, t274, t109, t11 * t217 + t113 * t35 - t178 * t43 - t18 * t381 + t257 * t4 + t282 * t287 + t302 * t78 + t55 * t71 + t333, -t10 * t217 - t113 * t34 + t178 * t419 + t19 * t381 - t257 * t3 + t282 * t450 + t55 * t72 + t78 * t98 - t334, -t10 * t302 - t11 * t98 - t18 * t72 - t19 * t71 + t287 * t3 + t34 * t43 - t35 * t419 + t4 * t450 - t389, -g(1) * t297 - g(2) * t272 + t19 * t10 + t18 * t11 - t113 * t282 + t3 * t419 + t4 * t43 + t55 * t78 - t284, t12, t6, -t276, t109, -t274, t325, t13 * t302 + t14 * t381 + t17 * t71 + t178 * t39 - t2 * t257 - t217 * t9 - t287 * t5 + t35 * t47 + t333, t1 * t287 + t14 * t72 - t16 * t71 - t2 * t450 - t302 * t8 - t34 * t39 - t35 * t38 + t9 * t98 - t389, t1 * t257 - t13 * t98 - t16 * t381 - t17 * t72 - t178 * t38 + t217 * t8 + t34 * t47 + t450 * t5 + t334, t1 * t38 + t16 * t8 + t5 * t47 + t17 * t13 + t2 * t39 + t14 * t9 - g(1) * (-pkin(5) * t83 - qJ(6) * t82 + t297) - g(2) * (pkin(5) * t87 + qJ(6) * t86 + t272) - t284; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t208, t386 * t260, t373, t208, t237, qJDD(2), -t229 - t247 + (t329 + t432) * t255 (-t414 + t432) * t257 + t370, 0, 0, t62, -t456, t68, t290, t439, t208, -pkin(2) * t131 + t351 * t252 + ((-qJ(3) * t382 - t101) * t255 + (t119 - t322) * t257) * qJD(1) + t345, t338 - pkin(2) * t132 + t271 * t251 + ((-qJ(3) * t376 + t102) * t255 + (-t430 - t120 + (qJD(3) - t186) * t252) * t257) * qJD(1), t119 * t175 - t120 * t174 + (t101 * t383 + t51) * t252 + (t102 * t383 + t312 - t50) * t251 + t318, -t149 * pkin(2) - t102 * t120 - t101 * t119 - t186 * t232 - g(1) * (-pkin(2) * t400 + t215) - g(2) * (-pkin(2) * t401 + t210) - g(3) * t388 + (-t101 * t251 + t102 * t252) * qJD(3) + (-t50 * t251 + t51 * t252) * qJ(3), t62, t68, t456, t208, -t439, t290, -t458 * t131 + t359 * t252 + (t126 + t378) * t174 + (-t106 * t257 + t255 * t79 + (-t257 * t74 - t366) * t251) * qJD(1) + t345, -t105 * t174 - t106 * t175 + (-t383 * t79 + t36) * t252 + (t383 * t89 + t312 + t41) * t251 + t318, -t338 + t126 * t175 + t458 * t132 + (t166 + t359 + t305) * t251 + (t105 * t257 - t255 * t89 + (t366 + (-qJD(3) + t74) * t257) * t252) * qJD(1), -t89 * t105 - t74 * t126 - t79 * t106 - g(1) * t215 - g(2) * t210 - g(3) * t344 + (qJ(3) * t36 + qJD(3) * t89) * t252 + (qJ(3) * t41 + qJD(3) * t79 - qJD(4) * t74) * t251 + (-t37 + t305) * t458, t7, t438, t22, t299, t440, t184, t165 * t35 + t18 * t384 + t217 * t421 + t282 * t300 + t302 * t346 - t393 * t55 + t278, -t165 * t34 - t180 * t282 - t19 * t384 - t217 * t422 + t346 * t98 - t392 * t55 - t277, t18 * t392 - t4 * t180 + t19 * t393 + t3 * t300 + t302 * t33 - t421 * t98 + t291, -g(1) * t335 - g(2) * t336 - g(3) * t317 + t3 * t112 - t165 * t282 + t18 * t421 + t19 * t422 + t301 * t4 + t346 * t55 + t262, t7, t22, -t438, t184, -t440, t299, -t14 * t384 - t17 * t393 - t217 * t420 - t300 * t5 + t302 * t424 + t66 * t35 + t278, t1 * t300 - t14 * t392 + t16 * t393 + t2 * t180 + t29 * t302 + t420 * t98 + t291, t16 * t384 + t17 * t392 - t5 * t180 + t217 * t423 + t66 * t34 - t424 * t98 + t277, t1 * t112 + t5 * t66 - t2 * t301 - g(1) * (t130 * pkin(5) + t129 * qJ(6) + t335) - g(2) * (t128 * pkin(5) + t127 * qJ(6) + t336) - g(3) * (-pkin(5) * t147 + qJ(6) * t146 + t317) + t424 * t17 + t423 * t16 + t420 * t14 + t262; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t273, t75, t446, t101 * t175 - t102 * t174 + t271, 0, 0, 0, 0, 0, 0, t273, t446, -t75, -t89 * t174 - t79 * t175 + t264, 0, 0, 0, 0, 0, 0, t265, t279, t327, -t18 * t98 - t19 * t302 + t263, 0, 0, 0, 0, 0, 0, t265, t327, -t279, t14 * t98 - t16 * t302 + t263 - t442; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t174 * t175 + t293, t150 + t132, -t172 - t409, t74 * t175 + (-pkin(3) * t381 + t257 * t89) * qJD(1) + t283 + t307, 0, 0, 0, 0, 0, 0, t267, -t441, t268, t4 * t433 - t55 * t175 + t281 * t19 + (t3 - t418) * t254 + t283, 0, 0, 0, 0, 0, 0, t267, t268, t441, -t2 * t433 - t17 * t175 + t281 * t16 + (t14 * t217 + t1) * t254 + t283; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t426, -t371, t15, -t426, -t459, -t178, -t55 * t98 + t275 + t417, t302 * t55 - t269 + t418, 0, 0, t426, t15, t371, -t178, t459, -t426, -t302 * t46 - t266 + t417 - 0.2e1 * t431, pkin(5) * t34 - t35 * qJ(6) + (t16 - t19) * t98 + (t14 - t394) * t302, -0.2e1 * t410 - t17 * t302 + t46 * t98 + (0.2e1 * qJD(6) - t18) * t217 + t269, t1 * qJ(6) - t2 * pkin(5) - t17 * t46 - t14 * t19 - g(1) * (-pkin(5) * t86 + qJ(6) * t87) - g(2) * (-pkin(5) * t82 + qJ(6) * t83) - g(3) * t323 + t394 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t178 + t426, t15, -t436 - t437, -t16 * t217 + t266 + t431;];
tau_reg  = t20;
