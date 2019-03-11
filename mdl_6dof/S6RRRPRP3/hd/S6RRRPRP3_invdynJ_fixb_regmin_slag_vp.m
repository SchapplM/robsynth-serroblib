% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% tau_reg [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRPRP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:41:44
% EndTime: 2019-03-09 16:41:57
% DurationCPUTime: 7.43s
% Computational Cost: add. (12923->571), mult. (28664->705), div. (0->0), fcn. (21237->14), ass. (0->293)
t274 = sin(qJ(2));
t428 = cos(qJ(3));
t351 = qJD(1) * t428;
t273 = sin(qJ(3));
t276 = cos(qJ(2));
t377 = t273 * t276;
t194 = -qJD(1) * t377 - t274 * t351;
t265 = qJD(2) + qJD(3);
t269 = sin(pkin(10));
t270 = cos(pkin(10));
t167 = t194 * t270 - t265 * t269;
t272 = sin(qJ(5));
t427 = cos(qJ(5));
t181 = t269 * t194;
t440 = t265 * t270 + t181;
t311 = t427 * t440;
t112 = t272 * t167 + t311;
t456 = t112 ^ 2;
t369 = qJD(1) * t274;
t192 = t273 * t369 - t276 * t351;
t186 = qJD(5) + t192;
t455 = t112 * t186;
t211 = t274 * t428 + t377;
t161 = t265 * t211;
t354 = t428 * t276;
t365 = t274 * qJDD(1);
t121 = qJD(1) * t161 - qJDD(1) * t354 + t273 * t365;
t118 = qJDD(5) + t121;
t325 = t272 * t440;
t452 = -t167 * t427 + t325;
t429 = t452 ^ 2;
t454 = t186 * t452;
t268 = qJ(2) + qJ(3);
t261 = cos(t268);
t254 = g(3) * t261;
t260 = sin(t268);
t277 = cos(qJ(1));
t386 = t260 * t277;
t275 = sin(qJ(1));
t387 = t260 * t275;
t447 = g(1) * t386 + g(2) * t387;
t431 = t254 - t447;
t278 = -pkin(8) - pkin(7);
t366 = qJD(1) * qJD(2);
t347 = t276 * t366;
t163 = qJDD(2) * pkin(2) - t278 * (-t347 - t365);
t348 = t274 * t366;
t364 = t276 * qJDD(1);
t170 = t278 * (-t348 + t364);
t228 = t278 * t274;
t213 = qJD(1) * t228;
t415 = qJD(2) * pkin(2);
t201 = t213 + t415;
t229 = t278 * t276;
t215 = qJD(1) * t229;
t350 = qJD(3) * t428;
t368 = qJD(3) * t273;
t342 = -t428 * t163 - t273 * t170 + t201 * t368 - t215 * t350;
t263 = qJDD(2) + qJDD(3);
t421 = t263 * pkin(3);
t74 = qJDD(4) + t342 - t421;
t453 = t74 + t431;
t209 = t269 * t427 + t272 * t270;
t191 = t209 * qJD(5);
t449 = t209 * t192 + t191;
t307 = -t272 * t269 + t427 * t270;
t349 = qJD(5) * t427;
t367 = qJD(5) * t272;
t434 = -t269 * t367 + t270 * t349;
t448 = t307 * t192 + t434;
t247 = pkin(2) * t348;
t310 = -t273 * t274 + t354;
t160 = t265 * t310;
t286 = t160 * qJD(1);
t419 = t276 * pkin(2);
t257 = pkin(1) + t419;
t437 = -qJ(4) * t211 - t257;
t57 = t121 * pkin(3) - qJ(4) * t286 + t194 * qJD(4) + qJDD(1) * t437 + t247;
t292 = t273 * t163 - t170 * t428 + t201 * t350 + t215 * t368;
t70 = t263 * qJ(4) + t265 * qJD(4) + t292;
t27 = -t269 * t70 + t270 * t57;
t28 = t269 * t57 + t270 * t70;
t329 = -t27 * t269 + t28 * t270;
t446 = -t74 - t254;
t241 = pkin(2) * t350 + qJD(4);
t250 = pkin(2) * t273 + qJ(4);
t199 = (-pkin(9) - t250) * t269;
t262 = t270 * pkin(9);
t379 = t270 * t250;
t200 = t262 + t379;
t309 = t199 * t427 - t272 * t200;
t394 = t192 * t270;
t340 = -t194 * pkin(4) + pkin(9) * t394;
t149 = -pkin(3) * t194 + qJ(4) * t192;
t135 = pkin(2) * t369 + t149;
t197 = t273 * t215;
t158 = t213 * t428 + t197;
t92 = t270 * t135 - t158 * t269;
t65 = t340 + t92;
t395 = t192 * t269;
t363 = pkin(9) * t395;
t93 = t269 * t135 + t270 * t158;
t79 = t363 + t93;
t445 = -qJD(5) * t309 - t241 * t307 + t272 * t65 + t427 * t79;
t146 = t272 * t199 + t200 * t427;
t444 = -qJD(5) * t146 - t209 * t241 + t272 * t79 - t427 * t65;
t179 = pkin(4) * t395;
t443 = pkin(5) * t449 - qJ(6) * t448 - qJD(6) * t209 + t179;
t271 = -pkin(9) - qJ(4);
t225 = t271 * t269;
t380 = t270 * qJ(4);
t226 = t262 + t380;
t308 = t225 * t427 - t272 * t226;
t153 = t201 * t428 + t197;
t94 = t270 * t149 - t153 * t269;
t69 = t340 + t94;
t95 = t269 * t149 + t270 * t153;
t81 = t363 + t95;
t442 = -qJD(4) * t307 - qJD(5) * t308 + t272 * t69 + t427 * t81;
t166 = t272 * t225 + t226 * t427;
t441 = -qJD(4) * t209 - qJD(5) * t166 + t272 * t81 - t427 * t69;
t438 = t261 * pkin(3) + t260 * qJ(4);
t284 = t211 * qJDD(1) + t286;
t283 = t269 * t263 + t270 * t284;
t116 = t269 * t284;
t344 = t263 * t270 - t116;
t436 = t269 * t283 + t270 * t344;
t435 = -t167 * t269 + t270 * t440;
t433 = t428 * t228 + t273 * t229;
t198 = t428 * t215;
t157 = t273 * t213 - t198;
t432 = -pkin(2) * t368 + t157;
t333 = g(1) * t277 + g(2) * t275;
t152 = -pkin(3) * t310 + t437;
t172 = t273 * t228 - t229 * t428;
t101 = t270 * t152 - t172 * t269;
t391 = t211 * t270;
t80 = -pkin(4) * t310 - pkin(9) * t391 + t101;
t102 = t269 * t152 + t270 * t172;
t392 = t211 * t269;
t90 = -pkin(9) * t392 + t102;
t316 = t272 * t80 + t427 * t90;
t399 = t160 * t270;
t356 = qJD(2) * t278;
t214 = t274 * t356;
t216 = t276 * t356;
t113 = qJD(3) * t433 + t428 * t214 + t273 * t216;
t361 = t274 * t415;
t91 = pkin(3) * t161 - qJ(4) * t160 - qJD(4) * t211 + t361;
t49 = -t113 * t269 + t270 * t91;
t36 = pkin(4) * t161 - pkin(9) * t399 + t49;
t400 = t160 * t269;
t50 = t270 * t113 + t269 * t91;
t45 = -pkin(9) * t400 + t50;
t430 = -qJD(5) * t316 - t272 * t45 + t36 * t427;
t426 = pkin(2) * t274;
t425 = pkin(5) * t118;
t253 = g(3) * t260;
t422 = t194 * pkin(5);
t420 = t270 * pkin(4);
t183 = t194 * qJ(6);
t416 = t183 - t445;
t227 = t257 * qJD(1);
t131 = pkin(3) * t192 + qJ(4) * t194 - t227;
t154 = t273 * t201 - t198;
t137 = qJ(4) * t265 + t154;
t83 = t270 * t131 - t137 * t269;
t58 = pkin(4) * t192 + pkin(9) * t167 + t83;
t84 = t269 * t131 + t270 * t137;
t63 = pkin(9) * t440 + t84;
t25 = t272 * t58 + t427 * t63;
t414 = t186 * t25;
t412 = -t422 - t444;
t411 = -t432 + t443;
t410 = -t154 + t443;
t409 = t183 - t442;
t408 = -t422 - t441;
t406 = qJ(6) * t118;
t405 = t452 * t112;
t404 = t121 * t269;
t403 = t309 * t118;
t402 = t146 * t118;
t401 = t154 * t265;
t398 = t308 * t118;
t397 = t166 * t118;
t396 = t186 * t194;
t393 = t194 * t192;
t264 = pkin(10) + qJ(5);
t258 = sin(t264);
t390 = t258 * t261;
t259 = cos(t264);
t389 = t259 * t261;
t388 = t260 * t271;
t251 = pkin(3) + t420;
t220 = t261 * t251;
t385 = t261 * t271;
t384 = t261 * t275;
t383 = t261 * t277;
t376 = t275 * t259;
t375 = t277 * t258;
t24 = -t272 * t63 + t427 * t58;
t374 = qJD(6) - t24;
t266 = t274 ^ 2;
t370 = -t276 ^ 2 + t266;
t362 = t428 * pkin(2);
t359 = pkin(5) * t389 + qJ(6) * t390 + t220;
t358 = g(3) * t390 - t258 * t447;
t357 = g(1) * t383 + g(2) * t384 + t253;
t346 = pkin(4) * t269 - t278;
t14 = t121 * pkin(4) - pkin(9) * t283 + t27;
t18 = pkin(9) * t344 + t28;
t345 = -t427 * t14 + t272 * t18 + t63 * t349 + t58 * t367;
t256 = -t362 - pkin(3);
t339 = -g(1) * t387 + g(2) * t386;
t338 = t179 - t432;
t337 = t83 * t194 + t270 * t447;
t336 = -pkin(3) * t260 - t426;
t173 = t258 * t384 + t259 * t277;
t175 = t261 * t375 - t376;
t335 = -g(1) * t173 + g(2) * t175;
t174 = t261 * t376 - t375;
t176 = t258 * t275 + t259 * t383;
t334 = g(1) * t174 - g(2) * t176;
t332 = g(1) * t275 - g(2) * t277;
t328 = -t269 * t83 + t270 * t84;
t324 = t220 - t388;
t323 = -g(3) * t389 + t259 * t447;
t133 = pkin(4) * t392 - t433;
t321 = pkin(5) * t259 + qJ(6) * t258 + t251;
t318 = -t272 * t90 + t427 * t80;
t315 = t332 * t261;
t314 = -0.2e1 * pkin(1) * t366 - pkin(7) * qJDD(2);
t313 = t272 * t14 + t427 * t18 + t58 * t349 - t367 * t63;
t312 = t272 * t36 + t80 * t349 - t367 * t90 + t427 * t45;
t114 = t273 * t214 - t216 * t428 + t228 * t368 - t229 * t350;
t41 = -qJD(5) * t311 - t167 * t367 - t272 * t344 - t427 * t283;
t82 = pkin(4) * t400 + t114;
t148 = -pkin(5) * t307 - t209 * qJ(6) - t251;
t304 = -t83 * t394 - t84 * t395 + t329 - t357;
t134 = -t265 * pkin(3) + qJD(4) - t153;
t302 = g(1) * t175 + g(2) * t173 + t258 * t253 - t345;
t301 = t41 - t455;
t279 = qJD(2) ^ 2;
t300 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t279 + t332;
t280 = qJD(1) ^ 2;
t299 = pkin(1) * t280 - pkin(7) * qJDD(1) + t333;
t42 = qJD(5) * t325 - t167 * t349 + t272 * t283 - t427 * t344;
t48 = -pkin(4) * t344 + t74;
t10 = t42 * pkin(5) + t41 * qJ(6) - qJD(6) * t452 + t48;
t22 = t186 * qJ(6) + t25;
t98 = -pkin(4) * t440 + t134;
t43 = -pkin(5) * t112 - qJ(6) * t452 + t98;
t298 = -t10 * t209 + t194 * t22 - t43 * t448 - t358;
t297 = -t25 * t194 + t48 * t209 + t448 * t98 + t358;
t296 = -t227 * t194 - t342 - t431;
t2 = qJD(6) * t186 + t313 + t406;
t21 = -t186 * pkin(5) + t374;
t4 = qJDD(6) + t345 - t425;
t295 = t2 * t307 + t4 * t209 + t21 * t448 - t22 * t449 - t357;
t294 = -t10 * t307 - t194 * t21 + t43 * t449 + t323;
t293 = t24 * t194 - t48 * t307 + t449 * t98 + t323;
t291 = t43 * t452 + qJDD(6) - t302;
t289 = -g(1) * t176 - g(2) * t174 - t253 * t259 + t313;
t288 = t134 * t394 - t84 * t194 + t269 * t453;
t287 = -t227 * t192 - t292 + t357;
t281 = t42 + t454;
t233 = t277 * t257;
t231 = qJ(4) * t383;
t230 = qJ(4) * t384;
t224 = t256 - t420;
t187 = -qJDD(1) * t257 + t247;
t141 = t307 * t211;
t140 = t209 * t211;
t136 = -t362 + t148;
t122 = -t192 ^ 2 + t194 ^ 2;
t117 = -t179 + t154;
t97 = -t194 * t265 - t121;
t96 = t192 * t265 + t284;
t68 = t160 * t209 + t211 * t434;
t67 = -t160 * t307 + t211 * t191;
t61 = t140 * pkin(5) - t141 * qJ(6) + t133;
t59 = pkin(5) * t452 - qJ(6) * t112;
t40 = pkin(5) * t310 - t318;
t39 = -qJ(6) * t310 + t316;
t32 = t112 * t194 + t118 * t307 - t186 * t449;
t31 = t118 * t209 + t186 * t448 + t194 * t452;
t23 = -t41 - t455;
t15 = t68 * pkin(5) + t67 * qJ(6) - t141 * qJD(6) + t82;
t11 = -t209 * t41 + t448 * t452;
t9 = -t161 * pkin(5) - t430;
t6 = qJ(6) * t161 - qJD(6) * t310 + t312;
t5 = t112 * t448 - t209 * t42 - t307 * t41 - t449 * t452;
t1 = [qJDD(1), t332, t333, qJDD(1) * t266 + 0.2e1 * t274 * t347, 0.2e1 * t274 * t364 - 0.2e1 * t366 * t370, qJDD(2) * t274 + t276 * t279, qJDD(2) * t276 - t274 * t279, 0, t274 * t314 + t276 * t300, -t274 * t300 + t276 * t314, -t194 * t160 + t211 * t284, -t211 * t121 - t160 * t192 + t194 * t161 + t284 * t310, t160 * t265 + t211 * t263, -t161 * t265 + t263 * t310, 0, -t114 * t265 - t121 * t257 - t161 * t227 - t187 * t310 + t192 * t361 + t263 * t433 + t315, -t113 * t265 - t227 * t160 - t172 * t263 + t187 * t211 - t194 * t361 - t257 * t284 + t339, t49 * t192 + t101 * t121 - t27 * t310 + t83 * t161 - t114 * t440 + t433 * t344 + t270 * t315 + (t134 * t160 + t74 * t211 - t333) * t269, -t50 * t192 - t102 * t121 + t28 * t310 - t84 * t161 - t114 * t167 - t433 * t283 + t74 * t391 + t134 * t399 - g(1) * (t269 * t384 + t270 * t277) - g(2) * (-t269 * t383 + t270 * t275) -t101 * t283 + t102 * t344 + t167 * t49 - t27 * t391 - t28 * t392 - t399 * t83 - t400 * t84 + t440 * t50 - t339, -g(2) * t233 + t27 * t101 + t28 * t102 + t134 * t114 - t74 * t433 + t83 * t49 + t84 * t50 + (g(1) * t278 - g(2) * t438) * t277 + (-g(1) * (-t257 - t438) + g(2) * t278) * t275, -t141 * t41 - t452 * t67, -t112 * t67 + t140 * t41 - t141 * t42 - t452 * t68, t118 * t141 + t161 * t452 - t186 * t67 + t310 * t41, t112 * t161 - t118 * t140 - t186 * t68 + t310 * t42, -t118 * t310 + t161 * t186, -t112 * t82 + t318 * t118 + t133 * t42 + t48 * t140 + t24 * t161 + t186 * t430 + t310 * t345 + t98 * t68 + t334, -t118 * t316 - t133 * t41 + t48 * t141 - t25 * t161 - t186 * t312 + t310 * t313 + t452 * t82 - t98 * t67 + t335, t10 * t140 - t112 * t15 - t118 * t40 - t161 * t21 - t186 * t9 + t310 * t4 + t42 * t61 + t43 * t68 + t334, t112 * t6 - t140 * t2 + t141 * t4 - t21 * t67 - t22 * t68 - t39 * t42 - t40 * t41 + t452 * t9 - t339, -t10 * t141 + t118 * t39 - t15 * t452 + t161 * t22 + t186 * t6 - t2 * t310 + t41 * t61 + t43 * t67 - t335, t2 * t39 + t22 * t6 + t10 * t61 + t43 * t15 + t4 * t40 + t21 * t9 - g(1) * (-pkin(5) * t174 - qJ(6) * t173) - g(2) * (pkin(5) * t176 + qJ(6) * t175 + t233) + (-g(1) * t346 - g(2) * t324) * t277 + (-g(1) * (-t257 - t324) - g(2) * t346) * t275; 0, 0, 0, -t274 * t280 * t276, t370 * t280, t365, t364, qJDD(2), -g(3) * t276 + t274 * t299, g(3) * t274 + t276 * t299, -t393, t122, t96, t97, t263, t157 * t265 + (-t192 * t369 + t263 * t428 - t265 * t368) * pkin(2) + t296, t158 * t265 + (t194 * t369 - t263 * t273 - t265 * t350) * pkin(2) + t287, -t250 * t404 - t256 * t344 + (-t92 + (t134 - t241) * t269) * t192 + t337 + t432 * t440 + t446 * t270, -t121 * t379 + t167 * t432 + t93 * t192 - t241 * t394 + t256 * t283 + t288, -t167 * t92 + t241 * t435 + t250 * t436 - t440 * t93 + t304, t74 * t256 - t84 * t93 - t83 * t92 - g(1) * (t277 * t336 + t231) - g(2) * (t275 * t336 + t230) - g(3) * (t438 + t419) + t329 * t250 + t328 * t241 - t432 * t134, t11, t5, t31, t32, t396, -t112 * t338 + t186 * t444 + t224 * t42 + t293 + t403, t186 * t445 - t224 * t41 + t338 * t452 + t297 - t402, -t112 * t411 + t136 * t42 - t186 * t412 + t294 + t403, t112 * t416 - t146 * t42 + t309 * t41 + t412 * t452 + t295, t136 * t41 + t186 * t416 - t411 * t452 + t298 + t402, t2 * t146 + t10 * t136 - t4 * t309 - g(3) * (t359 - t388 + t419) + t411 * t43 + t416 * t22 + t412 * t21 + t333 * (t260 * t321 + t385 + t426); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t393, t122, t96, t97, t263, t296 + t401, t153 * t265 + t287, -qJ(4) * t404 - pkin(3) * t116 + t154 * t181 + (t401 + t421 + t446) * t270 + (-t94 + (-qJD(4) + t134) * t269) * t192 + t337, -pkin(3) * t283 - qJD(4) * t394 - t121 * t380 + t154 * t167 + t95 * t192 + t288, qJ(4) * t436 + qJD(4) * t435 - t167 * t94 - t440 * t95 + t304, -t74 * pkin(3) - t84 * t95 - t83 * t94 - t134 * t154 - g(1) * (-pkin(3) * t386 + t231) - g(2) * (-pkin(3) * t387 + t230) - g(3) * t438 + t328 * qJD(4) + t329 * qJ(4), t11, t5, t31, t32, t396, t112 * t117 + t186 * t441 - t251 * t42 + t293 + t398, -t117 * t452 + t186 * t442 + t251 * t41 + t297 - t397, -t112 * t410 + t148 * t42 - t186 * t408 + t294 + t398, t112 * t409 - t166 * t42 + t308 * t41 + t408 * t452 + t295, t148 * t41 + t186 * t409 - t410 * t452 + t298 + t397, t2 * t166 + t10 * t148 - t4 * t308 - g(3) * t359 + t410 * t43 + t333 * t385 + t409 * t22 + t408 * t21 + (g(3) * t271 + t321 * t333) * t260; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t167 * t192 - t344, t192 * t440 + t283, -t167 ^ 2 - t440 ^ 2, -t167 * t83 - t440 * t84 + t453, 0, 0, 0, 0, 0, t281, -t301, t281, -t429 - t456, t301, -t22 * t112 - t21 * t452 + t10 + t431; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t405, t429 - t456, t23, -t42 + t454, t118, -t452 * t98 + t302 + t414, -t112 * t98 + t186 * t24 - t289, t112 * t59 - t291 + t414 + 0.2e1 * t425, pkin(5) * t41 - qJ(6) * t42 + (t22 - t25) * t452 - (t21 - t374) * t112, 0.2e1 * t406 + t112 * t43 + t452 * t59 + (0.2e1 * qJD(6) - t24) * t186 + t289, t2 * qJ(6) - t4 * pkin(5) - t43 * t59 - t21 * t25 - g(1) * (-pkin(5) * t175 + qJ(6) * t176) - g(2) * (-pkin(5) * t173 + qJ(6) * t174) - (-pkin(5) * t258 + qJ(6) * t259) * t253 + t374 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t405 - t118, t23, -t186 ^ 2 - t429, -t186 * t22 + t291 - t425;];
tau_reg  = t1;
