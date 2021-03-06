% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% tau_reg [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRPR13_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR13_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR13_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR13_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:29:38
% EndTime: 2019-03-09 11:30:01
% DurationCPUTime: 10.24s
% Computational Cost: add. (9764->691), mult. (23639->936), div. (0->0), fcn. (18280->14), ass. (0->316)
t293 = cos(qJ(2));
t284 = sin(pkin(6));
t289 = sin(qJ(2));
t402 = qJD(1) * t289;
t377 = t284 * t402;
t286 = cos(pkin(6));
t403 = qJD(1) * t286;
t383 = pkin(1) * t403;
t408 = -pkin(8) * t377 + t293 * t383;
t390 = qJD(3) - t408;
t262 = qJD(2) + t403;
t288 = sin(qJ(4));
t292 = cos(qJ(4));
t401 = qJD(1) * t293;
t376 = t284 * t401;
t181 = t262 * t288 + t292 * t376;
t176 = qJD(6) + t181;
t353 = t288 * t376;
t183 = t262 * t292 - t353;
t238 = qJD(4) + t377;
t283 = sin(pkin(11));
t285 = cos(pkin(11));
t119 = t183 * t283 - t285 * t238;
t121 = t183 * t285 + t238 * t283;
t287 = sin(qJ(6));
t291 = cos(qJ(6));
t59 = t291 * t119 + t121 * t287;
t471 = t176 * t59;
t341 = pkin(4) * t292 + qJ(5) * t288;
t470 = (-pkin(3) - t341) * t377 - qJD(4) * t341 + qJD(5) * t292 - t390;
t389 = qJD(1) * qJD(2);
t369 = t289 * t389;
t352 = t284 * t369;
t386 = qJDD(1) * t293;
t467 = -t284 * t386 + t352;
t465 = -t283 * t287 + t291 * t285;
t464 = t465 * qJD(6);
t335 = t119 * t287 - t121 * t291;
t469 = t176 * t335;
t170 = t283 * t288 * t377 - t285 * t376;
t396 = qJD(4) * t288;
t468 = -t283 * t396 - t170;
t279 = t284 ^ 2;
t385 = 0.2e1 * t279;
t295 = -pkin(2) - pkin(9);
t394 = qJD(4) * t295;
t371 = t292 * t394;
t259 = pkin(2) * t377;
t438 = qJ(3) * t293;
t339 = pkin(9) * t289 - t438;
t404 = qJD(1) * t284;
t163 = t339 * t404 + t259;
t255 = pkin(8) * t376;
t196 = t289 * t383 + t255;
t165 = pkin(3) * t376 + t196;
t409 = t292 * t163 + t288 * t165;
t80 = qJ(5) * t376 + t409;
t442 = (t371 - t80) * t285 - t470 * t283;
t466 = -t283 * t80 + t285 * t470;
t294 = cos(qJ(1));
t413 = t293 * t294;
t290 = sin(qJ(1));
t418 = t289 * t290;
t207 = -t286 * t413 + t418;
t416 = t290 * t293;
t417 = t289 * t294;
t209 = t286 * t416 + t417;
t424 = t284 * t293;
t306 = g(1) * t209 + g(2) * t207 - g(3) * t424;
t208 = t286 * t417 + t416;
t210 = -t286 * t418 + t413;
t348 = -g(1) * t210 - g(2) * t208;
t422 = t285 * t287;
t216 = t283 * t291 + t422;
t313 = t216 * qJD(6);
t391 = pkin(3) * t377 + t390;
t366 = -qJ(3) * t289 - pkin(1);
t331 = -pkin(2) * t293 + t366;
t172 = t331 * t404;
t463 = t172 * t377 + qJDD(3);
t368 = t293 * t389;
t387 = qJDD(1) * t289;
t315 = t368 + t387;
t304 = t315 * t284;
t187 = qJDD(4) + t304;
t462 = -pkin(4) * t187 + qJDD(5);
t280 = pkin(11) + qJ(6);
t275 = sin(t280);
t276 = cos(t280);
t423 = t284 * t294;
t327 = -t207 * t288 + t292 * t423;
t461 = t208 * t276 + t275 * t327;
t460 = -t208 * t275 + t276 * t327;
t388 = qJDD(1) * t286;
t260 = qJDD(2) + t388;
t89 = -qJD(4) * t353 + t260 * t288 + (qJD(4) * t262 - t467) * t292;
t88 = -qJD(4) * t181 + t292 * t260 + t467 * t288;
t66 = -t285 * t187 + t283 * t88;
t67 = t187 * t283 + t285 * t88;
t19 = -qJD(6) * t335 + t287 * t67 + t291 * t66;
t459 = pkin(3) + pkin(8);
t458 = pkin(1) * t289;
t457 = pkin(2) * t260;
t455 = pkin(9) * t293;
t454 = pkin(10) * t285;
t453 = pkin(10) * t292;
t450 = pkin(10) + qJ(5);
t110 = t262 * t295 + t391;
t312 = t293 * t295 + t366;
t136 = t312 * t404;
t395 = qJD(4) * t292;
t426 = t284 * t289;
t264 = pkin(8) * t426;
t382 = pkin(1) * qJD(2) * t286;
t358 = qJD(1) * t382;
t381 = pkin(1) * t388;
t355 = qJD(2) * t255 + qJDD(1) * t264 + t289 * t358 - t293 * t381;
t334 = qJDD(3) + t355;
t72 = pkin(3) * t304 + t260 * t295 + t334;
t233 = pkin(2) * t352;
t398 = qJD(3) * t289;
t300 = qJD(2) * t339 - t398;
t79 = t233 + (qJD(1) * t300 + qJDD(1) * t312) * t284;
t317 = t110 * t395 - t136 * t396 + t288 * t72 + t292 * t79;
t14 = qJ(5) * t187 + qJD(5) * t238 + t317;
t239 = t260 * qJ(3);
t241 = t262 * qJD(3);
t356 = t467 * pkin(8) - t289 * t381 - t293 * t358;
t90 = -t239 - t241 + t356;
t75 = -t467 * pkin(3) - t90;
t24 = pkin(4) * t89 - qJ(5) * t88 - qJD(5) * t183 + t75;
t7 = t285 * t14 + t283 * t24;
t400 = qJD(2) * t289;
t374 = t284 * t400;
t254 = pkin(2) * t374;
t132 = t284 * t300 + t254;
t378 = -pkin(1) * t293 - pkin(2);
t140 = pkin(3) * t426 + t264 + (-pkin(9) + t378) * t286;
t407 = pkin(2) * t424 + qJ(3) * t426;
t161 = (-pkin(1) - t455) * t284 - t407;
t270 = t286 * t458;
t166 = (t424 * t459 + t270) * qJD(2);
t314 = t292 * t132 + t140 * t395 - t161 * t396 + t288 * t166;
t399 = qJD(2) * t293;
t37 = (qJ(5) * t399 + qJD(5) * t289) * t284 + t314;
t258 = t293 * t382;
t274 = t286 * qJD(3);
t139 = -t374 * t459 + t258 + t274;
t205 = t286 * t288 + t292 * t424;
t150 = -qJD(4) * t205 + t288 * t374;
t379 = t288 * t424;
t151 = -qJD(4) * t379 + t286 * t395 - t292 * t374;
t206 = t286 * t292 - t379;
t49 = pkin(4) * t151 - qJ(5) * t150 - qJD(5) * t206 + t139;
t16 = t283 * t49 + t285 * t37;
t65 = t288 * t110 + t292 * t136;
t53 = qJ(5) * t238 + t65;
t242 = t262 * qJ(3);
t129 = t242 + t165;
t68 = pkin(4) * t181 - qJ(5) * t183 + t129;
t29 = t283 * t68 + t285 * t53;
t410 = t288 * t140 + t292 * t161;
t76 = qJ(5) * t426 + t410;
t406 = pkin(8) * t424 + t270;
t184 = -t286 * qJ(3) - t406;
t160 = pkin(3) * t424 - t184;
t86 = pkin(4) * t205 - qJ(5) * t206 + t160;
t40 = t283 * t86 + t285 * t76;
t87 = qJDD(6) + t89;
t449 = t465 * t87;
t448 = t216 * t87;
t447 = t283 * t89;
t446 = t285 * t89;
t420 = t288 * t289;
t171 = (t283 * t293 + t285 * t420) * t404;
t445 = t170 * t287 - t171 * t291 - t292 * t313 - t396 * t465;
t444 = -t171 * t287 + t468 * t291 + t292 * t464 - t396 * t422;
t443 = -t283 * t371 - t466;
t364 = pkin(5) * t283 - t295;
t144 = t288 * t163;
t81 = -pkin(4) * t376 - t165 * t292 + t144;
t441 = -pkin(5) * t170 - t364 * t396 - t81;
t440 = t181 * t465 + t464;
t439 = t216 * t181 + t313;
t108 = pkin(4) * t183 + qJ(5) * t181;
t64 = t110 * t292 - t288 * t136;
t42 = t283 * t108 + t285 * t64;
t436 = t181 * t238;
t435 = t181 * t283;
t434 = t183 * t238;
t321 = t238 * t292;
t431 = t260 * t286;
t430 = t275 * t288;
t429 = t276 * t288;
t428 = t279 * qJD(1) ^ 2;
t425 = t284 * t290;
t421 = t288 * t187;
t419 = t288 * t295;
t412 = t295 * t187;
t52 = -pkin(4) * t238 + qJD(5) - t64;
t411 = -qJD(5) + t52;
t340 = pkin(4) * t288 - qJ(5) * t292;
t221 = qJ(3) + t340;
t174 = t283 * t221 + t285 * t419;
t281 = t289 ^ 2;
t405 = -t293 ^ 2 + t281;
t397 = qJD(4) * t176;
t392 = qJD(2) - t262;
t384 = g(3) * t426;
t380 = t293 * t428;
t375 = t292 * t402;
t373 = t284 * t399;
t370 = g(3) * t407;
t365 = -t283 * t295 + pkin(5);
t6 = -t14 * t283 + t285 * t24;
t15 = -t283 * t37 + t285 * t49;
t28 = -t283 * t53 + t285 * t68;
t39 = -t283 * t76 + t285 * t86;
t41 = t285 * t108 - t283 * t64;
t362 = t140 * t292 - t288 * t161;
t361 = -t110 * t396 - t136 * t395 - t288 * t79 + t292 * t72;
t360 = t262 + t403;
t359 = t260 + t388;
t357 = t289 * t380;
t354 = t284 * t375;
t152 = -t209 * t292 + t288 * t425;
t326 = t207 * t292 + t288 * t423;
t351 = g(1) * t326 + g(2) * t152;
t349 = g(1) * t207 - g(2) * t209;
t347 = g(1) * t294 + g(2) * t290;
t346 = -t283 * t6 + t285 * t7;
t4 = pkin(5) * t89 - pkin(10) * t67 + t6;
t5 = -pkin(10) * t66 + t7;
t345 = t287 * t4 + t291 * t5;
t141 = -t283 * t453 + t174;
t344 = -pkin(5) * t354 - pkin(10) * t171 + qJD(6) * t141 - (t288 * t454 + t292 * t365) * qJD(4) + t466;
t212 = t285 * t221;
t130 = -t285 * t453 + t288 * t365 + t212;
t343 = t468 * pkin(10) - qJD(6) * t130 - t442;
t20 = pkin(5) * t181 - pkin(10) * t121 + t28;
t21 = -pkin(10) * t119 + t29;
t8 = t20 * t291 - t21 * t287;
t9 = t20 * t287 + t21 * t291;
t149 = t206 * t285 + t283 * t426;
t25 = pkin(5) * t205 - pkin(10) * t149 + t39;
t148 = t206 * t283 - t285 * t426;
t30 = -pkin(10) * t148 + t40;
t338 = t25 * t291 - t287 * t30;
t337 = t25 * t287 + t291 * t30;
t336 = -t28 * t283 + t285 * t29;
t82 = t291 * t148 + t149 * t287;
t83 = -t148 * t287 + t149 * t291;
t185 = -pkin(1) * t284 - t407;
t333 = qJD(2) * (-qJD(1) * t185 - t172);
t332 = -pkin(8) * t374 + t258;
t329 = t294 * pkin(1) + t210 * pkin(2) + pkin(8) * t425 + qJ(3) * t209;
t328 = -t288 * t132 - t140 * t396 - t161 * t395 + t166 * t292;
t230 = t450 * t283;
t325 = pkin(10) * t435 - qJD(5) * t285 + qJD(6) * t230 + t42;
t231 = t450 * t285;
t324 = pkin(5) * t183 + qJD(5) * t283 + qJD(6) * t231 + t181 * t454 + t41;
t323 = t238 * t288;
t320 = -qJ(3) * t399 - t398;
t167 = t284 * t320 + t254;
t97 = t233 + (qJD(1) * t320 + qJDD(1) * t331) * t284;
t319 = qJD(1) * t167 + qJDD(1) * t185 + t97;
t18 = -t59 * qJD(6) - t287 * t66 + t291 * t67;
t78 = -pkin(4) * t426 - t362;
t311 = g(1) * t152 - g(2) * t326 + g(3) * t205;
t153 = t209 * t288 + t292 * t425;
t310 = -g(1) * t153 + g(2) * t327 - g(3) * t206;
t309 = -pkin(1) * t290 - t208 * pkin(2) + pkin(8) * t423 - qJ(3) * t207;
t197 = t406 * qJD(2);
t308 = -g(1) * t208 + g(2) * t210 + t197 * t262;
t17 = -t361 + t462;
t307 = -t17 + t311;
t305 = -t348 + t384;
t303 = t348 - t356;
t302 = t306 - t355;
t301 = (t392 * t401 + t387) * t284;
t38 = -pkin(4) * t373 - t328;
t2 = -qJD(6) * t9 - t287 * t5 + t291 * t4;
t299 = t311 + t361;
t298 = t196 * t262 + t302;
t297 = -t238 * t394 - t305 + t75;
t273 = -pkin(5) * t285 - pkin(4);
t213 = t364 * t292;
t200 = t209 * pkin(2);
t198 = t207 * pkin(2);
t194 = -qJ(3) * t376 + t259;
t192 = t465 * t292;
t191 = t216 * t292;
t186 = t286 * t378 + t264;
t178 = t292 * t187;
t175 = -t274 - t332;
t173 = -t283 * t419 + t212;
t169 = t262 * t283 - t285 * t354;
t168 = t262 * t285 + t283 * t354;
t162 = -t242 - t196;
t158 = -pkin(2) * t262 + t390;
t112 = t150 * t285 + t283 * t373;
t111 = t150 * t283 - t285 * t373;
t105 = t334 - t457;
t96 = t153 * t276 + t210 * t275;
t95 = -t153 * t275 + t210 * t276;
t50 = pkin(5) * t148 + t78;
t48 = -pkin(5) * t435 + t65;
t43 = pkin(5) * t119 + t52;
t33 = qJD(6) * t83 + t291 * t111 + t112 * t287;
t32 = -qJD(6) * t82 - t111 * t287 + t112 * t291;
t26 = pkin(5) * t111 + t38;
t12 = -pkin(10) * t111 + t16;
t11 = pkin(5) * t151 - pkin(10) * t112 + t15;
t10 = pkin(5) * t66 + t17;
t1 = t8 * qJD(6) + t345;
t3 = [qJDD(1), g(1) * t290 - g(2) * t294, t347 (qJDD(1) * t281 + 0.2e1 * t289 * t368) * t279 (t289 * t386 - t389 * t405) * t385 (t289 * t359 + t360 * t399) * t284 (t293 * t359 - t360 * t400) * t284, t431, -t264 * t260 - t355 * t286 + (t293 * t431 + (-t369 + t386) * t385) * pkin(1) - t308, -pkin(1) * t315 * t385 - t260 * t406 - t262 * t332 + t286 * t356 - t349 ((qJD(2) * t158 - qJDD(1) * t184 - t90 + (qJD(2) * t186 - t175) * qJD(1)) * t293 + (qJD(2) * t162 + qJDD(1) * t186 + t105 + (qJD(2) * t184 + t197) * qJD(1)) * t289 - t347) * t284, t105 * t286 + t186 * t260 + (t289 * t333 + t293 * t319) * t284 + t308, -t175 * t262 - t184 * t260 - t286 * t90 + (-t289 * t319 + t293 * t333) * t284 + t349, -g(1) * t309 - g(2) * t329 + t105 * t186 + t158 * t197 + t162 * t175 + t172 * t167 + t90 * t184 + t97 * t185, t150 * t183 + t206 * t88, -t150 * t181 - t151 * t183 - t205 * t88 - t206 * t89, t150 * t238 + t187 * t206 + (t183 * t399 + t289 * t88) * t284, -t151 * t238 - t187 * t205 + (-t181 * t399 - t289 * t89) * t284 (t187 * t289 + t238 * t399) * t284, t328 * t238 + t362 * t187 + t139 * t181 + t160 * t89 + t75 * t205 + t129 * t151 - g(1) * t327 - g(2) * t153 + (t289 * t361 + t399 * t64) * t284, -t314 * t238 - t410 * t187 + t139 * t183 + t160 * t88 + t75 * t206 + t129 * t150 + (-t289 * t317 - t399 * t65) * t284 + t351, t15 * t181 + t39 * t89 + t6 * t205 + t28 * t151 + t38 * t119 + t78 * t66 + t17 * t148 + t52 * t111 - g(1) * (-t208 * t283 + t285 * t327) - g(2) * (t153 * t285 + t210 * t283) -t16 * t181 - t40 * t89 - t7 * t205 - t29 * t151 + t38 * t121 + t78 * t67 + t17 * t149 + t52 * t112 - g(1) * (-t208 * t285 - t283 * t327) - g(2) * (-t153 * t283 + t210 * t285) -t111 * t29 - t112 * t28 - t119 * t16 - t121 * t15 - t148 * t7 - t149 * t6 - t39 * t67 - t40 * t66 - t351, t7 * t40 + t29 * t16 + t6 * t39 + t28 * t15 + t17 * t78 + t52 * t38 - g(1) * (pkin(3) * t423 + pkin(4) * t327 - pkin(9) * t208 + qJ(5) * t326 + t309) - g(2) * (pkin(3) * t425 + pkin(4) * t153 + pkin(9) * t210 + qJ(5) * t152 + t329) t18 * t83 - t32 * t335, -t18 * t82 - t19 * t83 - t32 * t59 + t33 * t335, -t151 * t335 + t176 * t32 + t18 * t205 + t83 * t87, -t151 * t59 - t176 * t33 - t19 * t205 - t82 * t87, t151 * t176 + t205 * t87 (-qJD(6) * t337 + t11 * t291 - t12 * t287) * t176 + t338 * t87 + t2 * t205 + t8 * t151 + t26 * t59 + t50 * t19 + t10 * t82 + t43 * t33 - g(1) * t460 - g(2) * t96 -(qJD(6) * t338 + t11 * t287 + t12 * t291) * t176 - t337 * t87 - t1 * t205 - t9 * t151 - t26 * t335 + t50 * t18 + t10 * t83 + t43 * t32 + g(1) * t461 - g(2) * t95; 0, 0, 0, -t357, t405 * t428, t301 (-t392 * t402 + t386) * t284, t260, t428 * t458 + t298, pkin(1) * t380 + t262 * t408 - t303 + t384 ((-pkin(2) * t289 + t438) * qJDD(1) + ((-qJ(3) * qJD(2) - t162 - t196) * t289 + (-pkin(2) * qJD(2) - t158 + t390) * t293) * qJD(1)) * t284, -t194 * t376 - t298 - 0.2e1 * t457 + t463, 0.2e1 * t239 + t241 + t390 * t262 + (-g(3) * t289 + (t172 * t293 + t194 * t289) * qJD(1)) * t284 + t303, -t90 * qJ(3) - t105 * pkin(2) - t172 * t194 - t158 * t196 - g(1) * (qJ(3) * t210 - t200) - g(2) * (qJ(3) * t208 - t198) - t370 - t390 * t162, -t183 * t323 + t292 * t88 (-t89 - t434) * t292 + (-t88 + t436) * t288, -t238 * t396 + t178 + (-t183 * t293 - t238 * t420) * t404, -t238 * t395 - t421 + (t181 * t293 - t289 * t321) * t404, -t238 * t376, qJ(3) * t89 + t391 * t181 + t297 * t288 + t412 * t292 - t64 * t376 + (t144 + (t129 - t165) * t292) * t238, qJ(3) * t88 + t409 * t238 + t65 * t376 + t391 * t183 + (-t129 * t238 - t412) * t288 + t297 * t292, -t81 * t119 - t52 * t170 + t173 * t89 + t443 * t181 + t306 * t283 + (t6 + (t119 * t295 - t283 * t52) * qJD(4) - t305 * t285) * t288 + (t17 * t283 + t238 * t28 - t295 * t66) * t292, -t81 * t121 - t52 * t171 - t174 * t89 - t442 * t181 + t306 * t285 + (-t7 + (t121 * t295 - t285 * t52) * qJD(4) + t305 * t283) * t288 + (t17 * t285 - t238 * t29 - t295 * t67) * t292, t170 * t29 + t171 * t28 - t173 * t67 - t174 * t66 - t443 * t121 - t442 * t119 + (t28 * t285 + t283 * t29) * t396 + (-t283 * t7 - t285 * t6 + t305) * t292, t7 * t174 + t6 * t173 - t52 * t81 - g(1) * (-pkin(9) * t209 - t200) - g(2) * (-pkin(9) * t207 - t198) - t370 + (-t17 * t292 + t396 * t52) * t295 + t442 * t29 - g(3) * (t289 * t340 + t455) * t284 + t443 * t28 + t348 * t221, t18 * t192 - t335 * t445, -t18 * t191 - t19 * t192 + t335 * t444 - t445 * t59, t176 * t445 + t18 * t288 + t192 * t87 - t321 * t335, -t176 * t444 - t19 * t288 - t191 * t87 - t321 * t59, t176 * t321 + t288 * t87 (t130 * t291 - t141 * t287) * t87 + t2 * t288 + t8 * t395 + t213 * t19 + t10 * t191 - g(1) * (-t209 * t275 + t210 * t429) - g(2) * (-t207 * t275 + t208 * t429) + t441 * t59 + t444 * t43 + (t287 * t343 - t291 * t344) * t176 + (t8 * t375 - g(3) * (t275 * t293 + t276 * t420)) * t284 -(t130 * t287 + t141 * t291) * t87 - t1 * t288 - t9 * t395 + t213 * t18 + t10 * t192 - g(1) * (-t209 * t276 - t210 * t430) - g(2) * (-t207 * t276 - t208 * t430) - t441 * t335 + t445 * t43 + (t287 * t344 + t291 * t343) * t176 + (-t9 * t375 - g(3) * (-t275 * t420 + t276 * t293)) * t284; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t301, t260 + t357, -t262 ^ 2 - t281 * t428, t162 * t262 - t302 - t457 + t463, 0, 0, 0, 0, 0, -t181 * t262 - t238 * t323 + t178, -t183 * t262 - t238 * t321 - t421, -t292 * t66 + (-t283 * t395 - t168) * t181 + (t119 * t238 - t447) * t288, -t292 * t67 + (-t285 * t395 + t169) * t181 + (t121 * t238 - t446) * t288, t119 * t169 + t121 * t168 + (t283 * t67 - t285 * t66) * t288 + (-t119 * t285 + t121 * t283) * t395, -t168 * t28 - t169 * t29 + (qJD(4) * t336 - t17) * t292 + (t238 * t52 + t346) * t288 - t306, 0, 0, 0, 0, 0 -(t168 * t291 - t169 * t287) * t176 + (-t216 * t397 - t19) * t292 + (-t176 * t464 + t238 * t59 - t448) * t288 (t168 * t287 + t169 * t291) * t176 + (-t397 * t465 - t18) * t292 + (t176 * t313 - t238 * t335 - t449) * t288; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t183 * t181, -t181 ^ 2 + t183 ^ 2, t88 + t436, t434 - t89, t187, -t129 * t183 + t238 * t65 + t299, t129 * t181 + t238 * t64 - t310 - t317, -qJ(5) * t447 - pkin(4) * t66 - t119 * t65 - t183 * t28 + (t283 * t411 - t41) * t181 + t307 * t285, -qJ(5) * t446 - pkin(4) * t67 - t121 * t65 + t183 * t29 + (t285 * t411 + t42) * t181 - t307 * t283, t119 * t42 + t121 * t41 + (-qJ(5) * t66 - qJD(5) * t119 - t181 * t28 + t7) * t285 + (qJ(5) * t67 + qJD(5) * t121 - t181 * t29 - t6) * t283 + t310, -t28 * t41 - t29 * t42 - t52 * t65 + t336 * qJD(5) + t307 * pkin(4) + (t310 + t346) * qJ(5), t18 * t216 - t335 * t440, t18 * t465 - t19 * t216 + t335 * t439 - t440 * t59, t176 * t440 + t183 * t335 + t448, -t176 * t439 + t183 * t59 + t449, -t176 * t183 (-t230 * t291 - t231 * t287) * t87 + t273 * t19 - t10 * t465 - t8 * t183 - t48 * t59 + t439 * t43 + (t287 * t325 - t291 * t324) * t176 + t311 * t276 -(-t230 * t287 + t231 * t291) * t87 + t273 * t18 + t10 * t216 + t9 * t183 + t48 * t335 + t440 * t43 + (t287 * t324 + t291 * t325) * t176 - t311 * t275; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121 * t181 + t66, -t119 * t181 + t67, -t119 ^ 2 - t121 ^ 2, t119 * t29 + t121 * t28 - t299 + t462, 0, 0, 0, 0, 0, t19 - t469, t18 - t471; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t335 * t59, t335 ^ 2 - t59 ^ 2, t18 + t471, -t19 - t469, t87, t9 * t176 + t43 * t335 - g(1) * t95 - g(2) * t461 - g(3) * (-t206 * t275 + t276 * t426) + t2, t43 * t59 + g(1) * t96 - g(2) * t460 - g(3) * (-t206 * t276 - t275 * t426) - t345 + (t176 - qJD(6)) * t8;];
tau_reg  = t3;
