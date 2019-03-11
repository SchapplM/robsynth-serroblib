% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RRRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRRPP7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:26:32
% EndTime: 2019-03-09 21:27:05
% DurationCPUTime: 14.72s
% Computational Cost: add. (22269->746), mult. (58628->980), div. (0->0), fcn. (45787->10), ass. (0->308)
t307 = sin(qJ(2));
t440 = cos(pkin(6));
t383 = pkin(1) * t440;
t290 = t307 * t383;
t306 = sin(qJ(3));
t309 = cos(qJ(3));
t348 = pkin(3) * t306 - pkin(10) * t309;
t304 = sin(pkin(6));
t310 = cos(qJ(2));
t416 = t304 * t310;
t480 = (t290 + (pkin(8) + t348) * t416) * qJD(1) - t348 * qJD(3);
t362 = t310 * t383;
t397 = qJD(1) * t304;
t381 = t307 * t397;
t241 = -pkin(8) * t381 + qJD(1) * t362;
t330 = (pkin(2) * t307 - pkin(9) * t310) * t304;
t242 = qJD(1) * t330;
t171 = t309 * t241 + t306 * t242;
t160 = pkin(10) * t381 + t171;
t305 = sin(qJ(4));
t308 = cos(qJ(4));
t394 = qJD(3) * t306;
t385 = pkin(9) * t394;
t479 = t480 * t308 + (-t160 - t385) * t305;
t277 = -pkin(3) * t309 - pkin(10) * t306 - pkin(2);
t391 = qJD(4) * t308;
t478 = t308 * t160 - t277 * t391 + t480 * t305;
t396 = qJD(1) * t310;
t372 = t304 * t396;
t359 = t309 * t372;
t209 = t305 * t359 - t308 * t381;
t393 = qJD(3) * t309;
t346 = -t305 * t393 + t209;
t410 = t309 * t310;
t210 = (t305 * t307 + t308 * t410) * t397;
t477 = -t308 * t393 + t210;
t476 = t306 * t391 - t346;
t411 = t308 * t309;
t292 = pkin(9) * t411;
t360 = t306 * t372;
t389 = t308 * qJD(5);
t475 = pkin(4) * t360 - qJ(5) * t210 + t306 * t389 - (pkin(4) * t306 - qJ(5) * t411) * qJD(3) - (-t292 + (qJ(5) * t306 - t277) * t305) * qJD(4) + t479;
t413 = t306 * t308;
t474 = qJ(5) * t209 + (-pkin(9) * qJD(3) - qJ(5) * qJD(4)) * t413 + (-qJD(5) * t306 + (-pkin(9) * qJD(4) - qJ(5) * qJD(3)) * t309) * t305 - t478;
t366 = t440 * qJD(1);
t332 = t366 + qJD(2);
t223 = t306 * t332 + t309 * t381;
t278 = -qJD(3) + t372;
t180 = t223 * t305 + t308 * t278;
t182 = t223 * t308 - t278 * t305;
t303 = sin(pkin(11));
t439 = cos(pkin(11));
t114 = t439 * t180 + t182 * t303;
t361 = t306 * t381;
t221 = -t309 * t332 + t361;
t216 = qJD(4) + t221;
t432 = t114 * t216;
t395 = qJD(2) * t307;
t380 = t304 * t395;
t355 = qJD(1) * t380;
t392 = qJD(4) * t305;
t388 = qJD(1) * qJD(2);
t370 = t310 * t388;
t354 = t304 * t370;
t462 = qJD(3) * t361 - t309 * (qJD(3) * t332 + t354);
t109 = t223 * t392 + t278 * t391 - t305 * t355 + t308 * t462;
t110 = t223 * t391 - t278 * t392 - t305 * t462 - t308 * t355;
t61 = -t109 * t439 - t303 * t110;
t44 = t61 + t432;
t322 = -t303 * t180 + t182 * t439;
t473 = t114 * t322;
t368 = t439 * t305;
t259 = t303 * t308 + t368;
t248 = t259 * qJD(4);
t472 = t259 * t221 + t248;
t367 = t439 * t308;
t249 = qJD(4) * t367 - t303 * t392;
t418 = t303 * t305;
t258 = -t367 + t418;
t404 = t258 * t221 - t249;
t352 = t439 * t393;
t403 = t209 * t439 - t249 * t306 + t303 * t477 - t305 * t352;
t402 = t210 * t439 + t248 * t306 - t303 * t346 - t308 * t352;
t471 = t360 - t394;
t390 = t223 * qJD(3);
t184 = t306 * t354 + t390;
t239 = t259 * t306;
t326 = t278 * t306;
t60 = -t109 * t303 + t439 * t110;
t470 = t114 * t326 - t184 * t239 + t216 * t403 + t309 * t60;
t436 = t322 ^ 2;
t300 = t304 ^ 2;
t469 = -0.2e1 * t300 * t388;
t454 = t474 * t303 + t475 * t439;
t453 = t475 * t303 - t474 * t439;
t399 = pkin(8) * t416 + t290;
t237 = t440 * pkin(9) + t399;
t204 = qJD(2) * pkin(9) + qJD(1) * t237;
t238 = (-pkin(2) * t310 - pkin(9) * t307 - pkin(1)) * t304;
t215 = qJD(1) * t238;
t146 = -t306 * t204 + t215 * t309;
t136 = pkin(3) * t278 - t146;
t95 = pkin(4) * t180 + qJD(5) + t136;
t45 = pkin(5) * t114 - qJ(6) * t322 + t95;
t468 = t322 * t45;
t457 = -qJ(5) - pkin(10);
t369 = qJD(4) * t457;
t247 = t305 * t369 + t389;
t321 = -t305 * qJD(5) + t308 * t369;
t165 = pkin(3) * t223 + pkin(10) * t221;
t93 = -t146 * t305 + t308 * t165;
t71 = qJ(5) * t221 * t308 + pkin(4) * t223 + t93;
t421 = t221 * t305;
t94 = t308 * t146 + t305 * t165;
t82 = qJ(5) * t421 + t94;
t443 = (-t321 + t71) * t439 + (t247 - t82) * t303;
t339 = pkin(8) * t354;
t347 = qJD(2) * t366;
t108 = pkin(1) * t307 * t347 + t184 * pkin(3) + pkin(10) * t462 + t339;
t203 = -pkin(2) * t332 - t241;
t132 = t221 * pkin(3) - t223 * pkin(10) + t203;
t147 = t204 * t309 + t215 * t306;
t137 = -pkin(10) * t278 + t147;
t243 = qJD(2) * t330;
t232 = qJD(1) * t243;
t417 = t304 * t307;
t254 = -pkin(8) * t417 + t362;
t245 = t254 * qJD(2);
t233 = qJD(1) * t245;
t96 = -t204 * t394 + t215 * t393 + t306 * t232 + t309 * t233;
t89 = pkin(10) * t355 + t96;
t32 = t305 * t108 + t132 * t391 - t137 * t392 + t308 * t89;
t77 = t308 * t132 - t137 * t305;
t467 = -t216 * t77 + t32;
t78 = t132 * t305 + t137 * t308;
t33 = -qJD(4) * t78 + t308 * t108 - t305 * t89;
t466 = -t78 * t216 - t33;
t353 = -t147 + (t392 + t421) * pkin(4);
t170 = -t306 * t241 + t242 * t309;
t159 = -pkin(3) * t381 - t170;
t406 = pkin(4) * t476 + pkin(9) * t393 - t159;
t236 = -pkin(2) * t440 - t254;
t250 = t306 * t417 - t309 * t440;
t251 = t306 * t440 + t309 * t417;
t155 = t250 * pkin(3) - t251 * pkin(10) + t236;
t167 = t309 * t237 + t306 * t238;
t157 = -pkin(10) * t416 + t167;
t92 = t305 * t155 + t308 * t157;
t231 = t305 * t277 + t292;
t465 = t359 - t393;
t464 = t223 * t114 - t184 * t258 - t216 * t472;
t461 = t114 * t404 - t61 * t258 - t259 * t60 - t322 * t472;
t414 = t305 * t306;
t240 = -t303 * t414 + t306 * t367;
t460 = t114 * t402 - t61 * t239 - t240 * t60 + t322 * t403;
t459 = t216 ^ 2;
t311 = qJD(1) ^ 2;
t458 = pkin(9) * t309;
t299 = t306 * pkin(9);
t15 = pkin(4) * t184 + qJ(5) * t109 - qJD(5) * t182 + t33;
t18 = -qJ(5) * t110 - qJD(5) * t180 + t32;
t3 = t439 * t15 - t303 * t18;
t4 = t303 * t15 + t439 * t18;
t379 = qJD(2) * t416;
t193 = -qJD(3) * t250 + t309 * t379;
t194 = t251 * t305 + t308 * t416;
t131 = -qJD(4) * t194 + t193 * t308 + t305 * t380;
t192 = qJD(3) * t251 + t306 * t379;
t384 = t305 * t416;
t195 = t251 * t308 - t384;
t111 = -t237 * t394 + t238 * t393 + t306 * t243 + t309 * t245;
t101 = pkin(10) * t380 + t111;
t246 = t399 * qJD(2);
t123 = t192 * pkin(3) - t193 * pkin(10) + t246;
t43 = -qJD(4) * t92 - t305 * t101 + t308 * t123;
t23 = pkin(4) * t192 - qJ(5) * t131 - qJD(5) * t195 + t43;
t130 = -qJD(4) * t384 + t193 * t305 + t251 * t391 - t308 * t380;
t42 = t308 * t101 + t305 * t123 + t155 * t391 - t157 * t392;
t29 = -qJ(5) * t130 - qJD(5) * t194 + t42;
t8 = t303 * t23 + t439 * t29;
t456 = t471 * qJ(6) + qJD(6) * t309 + t453;
t455 = -t471 * pkin(5) - t454;
t63 = -qJ(5) * t182 + t77;
t55 = pkin(4) * t216 + t63;
t64 = -qJ(5) * t180 + t78;
t58 = t439 * t64;
t28 = t303 * t55 + t58;
t452 = t403 * pkin(5) - t402 * qJ(6) + qJD(6) * t240 - t406;
t91 = t308 * t155 - t157 * t305;
t69 = pkin(4) * t250 - qJ(5) * t195 + t91;
t79 = -qJ(5) * t194 + t92;
t39 = t303 * t69 + t439 * t79;
t41 = t303 * t71 + t439 * t82;
t30 = t303 * t63 + t58;
t450 = t30 * t322;
t449 = t303 * t64;
t363 = t204 * t393 + t215 * t394 - t309 * t232 + t306 * t233;
t90 = -pkin(3) * t355 + t363;
t448 = t305 * t90;
t447 = t308 * t90;
t445 = t472 * pkin(5) + t404 * qJ(6) - qJD(6) * t259 + t353;
t444 = t223 * pkin(5) + t443;
t177 = t247 * t439 + t303 * t321;
t36 = qJ(6) * t223 + t41;
t442 = t177 - t36;
t441 = t177 - t41;
t438 = t109 * t305;
t437 = t110 * t308;
t435 = t322 * t216;
t433 = t114 ^ 2;
t430 = t180 * t216;
t429 = t180 * t305;
t428 = t182 * t180;
t427 = t182 * t216;
t279 = t457 * t308;
t201 = -t279 * t303 - t368 * t457;
t426 = t184 * t201;
t202 = -t279 * t439 + t418 * t457;
t425 = t184 * t202;
t424 = t184 * t250;
t423 = t184 * t309;
t422 = t216 * t223;
t420 = t223 * t221;
t419 = t300 * t311;
t415 = t305 * t184;
t412 = t308 * t184;
t31 = t439 * t63 - t449;
t409 = qJD(6) - t31;
t408 = qJD(4) * t231 + t479;
t407 = -(-t308 * t394 - t309 * t392) * pkin(9) + t478;
t261 = t308 * t277;
t188 = -qJ(5) * t413 + t261 + (-pkin(9) * t305 - pkin(4)) * t309;
t198 = -qJ(5) * t414 + t231;
t139 = t303 * t188 + t439 * t198;
t265 = pkin(4) * t414 + t299;
t398 = t307 ^ 2 - t310 ^ 2;
t387 = t184 * qJ(6) + t4;
t297 = -pkin(4) * t308 - pkin(3);
t377 = t278 * t394;
t374 = t309 * t390;
t166 = -t306 * t237 + t238 * t309;
t365 = t216 * t308;
t364 = t307 * t310 * t419;
t2 = -pkin(5) * t184 - t3;
t358 = t223 * t372;
t357 = t30 * t216 + t3;
t356 = -t177 * t114 + t201 * t61 - t202 * t60;
t351 = t304 * t311 * t440;
t349 = pkin(1) * t469;
t156 = pkin(3) * t416 - t166;
t140 = t194 * t439 + t195 * t303;
t73 = t130 * t439 + t131 * t303;
t344 = t114 * t73 + t140 * t60;
t341 = -t305 * t78 - t308 * t77;
t340 = t306 * t363 + t96 * t309;
t337 = -t433 - t436;
t336 = -t433 + t436;
t334 = t300 * t307 * t370;
t331 = 0.2e1 * t366 + qJD(2);
t112 = -t237 * t393 - t238 * t394 + t243 * t309 - t306 * t245;
t329 = t60 + t435;
t328 = -t60 + t435;
t7 = t23 * t439 - t303 * t29;
t27 = t439 * t55 - t449;
t38 = -t303 * t79 + t439 * t69;
t325 = -pkin(10) * t184 + t136 * t216;
t138 = t188 * t439 - t303 * t198;
t324 = -t114 * t403 + t60 * t239;
t323 = t114 * t472 + t60 * t258;
t118 = pkin(4) * t194 + t156;
t319 = pkin(1) * (-t347 + t419);
t318 = -t61 + t432;
t141 = -t303 * t194 + t195 * t439;
t75 = -t303 * t130 + t131 * t439;
t317 = t114 * t75 + t140 * t61 + t141 * t60 + t322 * t73;
t316 = t114 * t192 + t140 * t184 + t216 * t73 + t250 * t60;
t102 = -pkin(3) * t380 - t112;
t56 = pkin(4) * t110 + t90;
t65 = pkin(4) * t130 + t102;
t9 = pkin(5) * t60 - qJ(6) * t61 - qJD(6) * t322 + t56;
t296 = -pkin(4) * t439 - pkin(5);
t293 = pkin(4) * t303 + qJ(6);
t244 = t399 * qJD(1);
t234 = qJD(1) * t246;
t230 = -t305 * t458 + t261;
t185 = pkin(5) * t258 - qJ(6) * t259 + t297;
t158 = pkin(5) * t239 - qJ(6) * t240 + t265;
t134 = t309 * pkin(5) - t138;
t133 = -qJ(6) * t309 + t139;
t100 = -t216 * t326 - t423;
t98 = t192 * t216 + t424;
t54 = pkin(4) * t182 + pkin(5) * t322 + qJ(6) * t114;
t52 = pkin(5) * t140 - qJ(6) * t141 + t118;
t48 = t259 * t184 - t216 * t404 - t223 * t322;
t35 = -t250 * pkin(5) - t38;
t34 = qJ(6) * t250 + t39;
t25 = qJ(6) * t216 + t28;
t24 = -t216 * pkin(5) + qJD(6) - t27;
t21 = t61 * t259 - t322 * t404;
t20 = t61 * t240 - t322 * t402;
t19 = t141 * t61 + t322 * t75;
t12 = pkin(5) * t73 - qJ(6) * t75 - qJD(6) * t141 + t65;
t11 = t184 * t240 - t216 * t402 - t309 * t61 - t322 * t326;
t10 = t141 * t184 + t192 * t322 + t216 * t75 + t250 * t61;
t6 = -t192 * pkin(5) - t7;
t5 = qJ(6) * t192 + qJD(6) * t250 + t8;
t1 = qJD(6) * t216 + t387;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t334, t398 * t469, t331 * t379, -0.2e1 * t334, -t331 * t380, 0, -t234 * t440 - t246 * t332 + t307 * t349, -t233 * t440 - t245 * t332 + t310 * t349 (t233 * t310 + t234 * t307 + (-t241 * t310 - t244 * t307) * qJD(2) + (t245 * t310 + t246 * t307 + (-t254 * t310 - t307 * t399) * qJD(2)) * qJD(1)) * t304, t233 * t399 - t234 * t254 - t241 * t246 + t244 * t245, t223 * t193 - t251 * t462, -t251 * t184 - t223 * t192 - t193 * t221 + t250 * t462, -t193 * t278 + t223 * t380 + t251 * t355 + t416 * t462, t192 * t221 + t424, t192 * t278 + (t184 * t310 + (-qJD(1) * t250 - t221) * t395) * t304 (-t278 * t304 - t300 * t396) * t395, -t112 * t278 + t184 * t236 + t192 * t203 + t221 * t246 + t234 * t250 + (t310 * t363 + (qJD(1) * t166 + t146) * t395) * t304, t111 * t278 - t147 * t380 - t167 * t355 + t203 * t193 + t246 * t223 + t234 * t251 - t236 * t462 + t416 * t96, -t111 * t221 - t112 * t223 - t146 * t193 - t147 * t192 + t166 * t462 - t167 * t184 - t96 * t250 + t251 * t363, t111 * t147 + t112 * t146 - t166 * t363 + t167 * t96 + t203 * t246 + t234 * t236, -t109 * t195 + t131 * t182, t109 * t194 - t110 * t195 - t130 * t182 - t131 * t180, -t109 * t250 + t131 * t216 + t182 * t192 + t184 * t195, t110 * t194 + t130 * t180, -t110 * t250 - t130 * t216 - t180 * t192 - t184 * t194, t98, t102 * t180 + t110 * t156 + t130 * t136 + t184 * t91 + t192 * t77 + t194 * t90 + t216 * t43 + t250 * t33, t102 * t182 - t109 * t156 + t131 * t136 - t184 * t92 - t192 * t78 + t195 * t90 - t216 * t42 - t250 * t32, t109 * t91 - t110 * t92 - t130 * t78 - t131 * t77 - t180 * t42 - t182 * t43 - t194 * t32 - t195 * t33, t102 * t136 + t156 * t90 + t32 * t92 + t33 * t91 + t42 * t78 + t43 * t77, t19, -t317, t10, t344, -t316, t98, t114 * t65 + t118 * t60 + t140 * t56 + t184 * t38 + t192 * t27 + t216 * t7 + t250 * t3 + t73 * t95, t118 * t61 + t141 * t56 - t184 * t39 - t192 * t28 - t216 * t8 - t250 * t4 + t322 * t65 + t75 * t95, -t114 * t8 - t140 * t4 - t141 * t3 - t27 * t75 - t28 * t73 - t322 * t7 - t38 * t61 - t39 * t60, t118 * t56 + t27 * t7 + t28 * t8 + t3 * t38 + t39 * t4 + t65 * t95, t19, t10, t317, t98, t316, t344, t114 * t12 + t140 * t9 - t184 * t35 - t192 * t24 - t2 * t250 - t216 * t6 + t45 * t73 + t52 * t60, -t1 * t140 - t114 * t5 + t141 * t2 + t24 * t75 - t25 * t73 + t322 * t6 - t34 * t60 + t35 * t61, t1 * t250 - t12 * t322 - t141 * t9 + t184 * t34 + t192 * t25 + t216 * t5 - t45 * t75 - t52 * t61, t1 * t34 + t12 * t45 + t2 * t35 + t24 * t6 + t25 * t5 + t52 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t364, t398 * t419, -t310 * t351, t364, t307 * t351, 0, t244 * t332 + t307 * t319 - t339, pkin(8) * t355 + t241 * t332 + t310 * t319, 0, 0, -t306 * t462 - t309 * t358 + t374, -t309 * t462 + t465 * t221 + (-t184 + t358 - t390) * t306, -t278 * t393 + (t278 * t410 + (t306 * qJD(2) - t223) * t307) * t397, -t221 * t326 - t423, t377 + (-t310 * t326 + (qJD(2) * t309 + t221) * t307) * t397, t278 * t381, -pkin(2) * t184 + t170 * t278 - t221 * t244 - t234 * t309 + (t203 * t306 + t278 * t458) * qJD(3) + (-t146 * t307 + (-pkin(9) * t395 - t203 * t310) * t306) * t397, pkin(2) * t462 - pkin(9) * t377 + t147 * t381 - t171 * t278 - t203 * t465 - t244 * t223 + t234 * t306 - t355 * t458, t170 * t223 - t462 * t299 + t340 + (t171 + t385) * t221 + t471 * t147 + t465 * t146 + (t374 - t423) * pkin(9), -pkin(2) * t234 - t146 * t170 - t147 * t171 - t203 * t244 + ((-t146 * t309 - t147 * t306) * qJD(3) + t340) * pkin(9), -t109 * t413 + (-t306 * t392 - t477) * t182, t180 * t210 + t182 * t209 + (-t180 * t308 - t182 * t305) * t393 + (t438 - t437 + (-t182 * t308 + t429) * qJD(4)) * t306, t109 * t309 - t477 * t216 + (-t182 * t278 - t216 * t392 + t412) * t306, t110 * t414 + t180 * t476, t110 * t309 + t346 * t216 + (t180 * t278 - t216 * t391 - t415) * t306, t100, -t136 * t209 - t159 * t180 + t184 * t230 - t408 * t216 + (-t33 + (pkin(9) * t180 + t136 * t305) * qJD(3)) * t309 + (pkin(9) * t110 + t136 * t391 - t278 * t77 + t448) * t306, -t136 * t210 - t159 * t182 - t184 * t231 + t407 * t216 + (t32 + (pkin(9) * t182 + t136 * t308) * qJD(3)) * t309 + (-pkin(9) * t109 - t136 * t392 + t278 * t78 + t447) * t306, t109 * t230 - t110 * t231 + t209 * t78 + t210 * t77 + t408 * t182 + t407 * t180 + t341 * t393 + (-t305 * t32 - t308 * t33 + (t305 * t77 - t308 * t78) * qJD(4)) * t306, -t136 * t159 + t230 * t33 + t231 * t32 - t407 * t78 - t408 * t77 + (t136 * t393 + t306 * t90) * pkin(9), t20, t460, t11, t324, t470, t100, t114 * t406 + t138 * t184 - t216 * t454 + t239 * t56 + t265 * t60 - t27 * t326 - t3 * t309 - t403 * t95, -t139 * t184 + t216 * t453 + t240 * t56 + t265 * t61 + t28 * t326 + t309 * t4 + t322 * t406 - t402 * t95, t114 * t453 - t138 * t61 - t139 * t60 - t239 * t4 - t240 * t3 + t27 * t402 + t28 * t403 + t322 * t454, t138 * t3 + t139 * t4 + t265 * t56 - t27 * t454 - t28 * t453 + t406 * t95, t20, t11, -t460, t100, -t470, t324, -t114 * t452 - t134 * t184 + t158 * t60 + t2 * t309 + t216 * t455 + t239 * t9 + t24 * t326 - t403 * t45, -t1 * t239 + t114 * t456 - t133 * t60 + t134 * t61 + t2 * t240 - t24 * t402 + t25 * t403 - t322 * t455, -t1 * t309 + t133 * t184 - t158 * t61 - t216 * t456 - t240 * t9 - t25 * t326 + t322 * t452 + t402 * t45, t1 * t133 + t134 * t2 + t158 * t9 - t24 * t455 - t25 * t456 - t45 * t452; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t420, -t221 ^ 2 + t223 ^ 2, -t221 * t278 - t462, -t420, -t223 * t278 - t184, t355, -t147 * t278 - t203 * t223 - t363, -t146 * t278 + t203 * t221 - t96, 0, 0, t182 * t365 - t438 (-t109 - t430) * t308 + (-t110 - t427) * t305, -t182 * t223 + t216 * t365 + t415, t216 * t429 - t437, t180 * t223 - t305 * t459 + t412, -t422, -pkin(3) * t110 - t147 * t180 - t223 * t77 - t447 + (-pkin(10) * t391 - t93) * t216 + t325 * t305, pkin(3) * t109 - t147 * t182 + t223 * t78 + t448 + (pkin(10) * t392 + t94) * t216 + t325 * t308, t180 * t94 + t182 * t93 + ((qJD(4) * t182 - t110) * pkin(10) + t467) * t308 + ((qJD(4) * t180 - t109) * pkin(10) + t466) * t305, -pkin(3) * t90 - t136 * t147 - t77 * t93 - t78 * t94 + (qJD(4) * t341 - t305 * t33 + t308 * t32) * pkin(10), t21, t461, t48, t323, t464, -t422, t114 * t353 - t216 * t443 - t223 * t27 + t258 * t56 + t297 * t60 + t472 * t95 - t426, -t216 * t441 + t223 * t28 + t259 * t56 + t297 * t61 + t322 * t353 - t404 * t95 - t425, t114 * t41 - t258 * t4 - t259 * t3 + t27 * t404 - t28 * t472 + t322 * t443 + t356, -t201 * t3 + t202 * t4 - t27 * t443 + t28 * t441 + t297 * t56 + t353 * t95, t21, t48, -t461, -t422, -t464, t323, t114 * t445 + t185 * t60 - t216 * t444 + t223 * t24 + t258 * t9 + t45 * t472 - t426, -t1 * t258 + t114 * t36 + t2 * t259 - t24 * t404 - t25 * t472 + t322 * t444 + t356, -t185 * t61 + t216 * t442 - t223 * t25 - t259 * t9 - t322 * t445 + t404 * t45 + t425, t1 * t202 + t185 * t9 + t2 * t201 + t24 * t444 + t25 * t442 + t445 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t428, -t180 ^ 2 + t182 ^ 2, -t109 + t430, -t428, -t110 + t427, t184, -t136 * t182 - t466, t136 * t180 - t467, 0, 0, t473, t336, t44, -t473, t328, t184, -t95 * t322 + (-t114 * t182 + t184 * t439) * pkin(4) + t357, t95 * t114 + t31 * t216 + (-t182 * t322 - t184 * t303) * pkin(4) - t4, t28 * t322 - t450 + (-t303 * t60 - t439 * t61) * pkin(4) + (t31 - t27) * t114, t27 * t30 - t28 * t31 + (-t182 * t95 + t3 * t439 + t303 * t4) * pkin(4), t473, t44, -t336, t184, -t328, -t473, -t468 - t114 * t54 + (pkin(5) - t296) * t184 + t357, t25 * t322 - t293 * t60 + t296 * t61 - t450 + (t24 - t409) * t114, -t114 * t45 + t322 * t54 + t184 * t293 + (0.2e1 * qJD(6) - t31) * t216 + t387, t1 * t293 + t2 * t296 - t24 * t30 + t25 * t409 - t45 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t329, -t318, t337, t114 * t28 + t27 * t322 + t56, 0, 0, 0, 0, 0, 0, t329, t337, t318, t114 * t25 - t24 * t322 + t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t184 + t473, t44, -t459 - t436, -t216 * t25 + t2 + t468;];
tauc_reg  = t13;
