% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% tau_reg [6x38]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRRRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:29:51
% EndTime: 2019-03-10 03:30:10
% DurationCPUTime: 8.07s
% Computational Cost: add. (17214->491), mult. (41845->633), div. (0->0), fcn. (33222->18), ass. (0->283)
t261 = qJD(2) + qJD(3);
t268 = sin(qJ(3));
t274 = cos(qJ(2));
t421 = cos(qJ(3));
t354 = qJD(1) * t421;
t269 = sin(qJ(2));
t372 = qJD(1) * t269;
t463 = -t268 * t372 + t274 * t354;
t481 = t463 * t261;
t271 = cos(qJ(6));
t266 = sin(qJ(5));
t272 = cos(qJ(5));
t343 = qJDD(1) * t421;
t362 = t274 * qJDD(1);
t126 = t268 * t362 + t269 * t343 + t481;
t384 = t268 * t274;
t197 = t421 * t269 + t384;
t159 = t261 * t197;
t363 = t269 * qJDD(1);
t321 = t268 * t363 - t274 * t343;
t127 = qJD(1) * t159 + t321;
t267 = sin(qJ(4));
t273 = cos(qJ(4));
t185 = -qJD(1) * t384 - t269 * t354;
t314 = t273 * t185 - t267 * t463;
t292 = qJD(4) * t314 - t267 * t126 - t273 * t127;
t315 = t267 * t185 + t273 * t463;
t443 = t266 * t315 - t272 * t314;
t369 = qJD(4) * t273;
t370 = qJD(4) * t267;
t72 = t273 * t126 - t267 * t127 + t185 * t370 + t369 * t463;
t30 = qJD(5) * t443 + t266 * t72 - t272 * t292;
t28 = qJDD(6) + t30;
t26 = t271 * t28;
t101 = t266 * t314 + t272 * t315;
t448 = qJD(6) - t101;
t265 = sin(qJ(6));
t449 = t265 * t448;
t256 = qJD(4) + t261;
t242 = qJD(5) + t256;
t86 = -t271 * t242 + t265 * t443;
t480 = t86 * t443 - t448 * t449 + t26;
t260 = qJDD(2) + qJDD(3);
t255 = qJDD(4) + t260;
t239 = qJDD(5) + t255;
t367 = qJD(5) * t272;
t368 = qJD(5) * t266;
t29 = t266 * t292 + t272 * t72 + t314 * t368 + t315 * t367;
t365 = qJD(6) * t271;
t366 = qJD(6) * t265;
t21 = t265 * t239 + t242 * t365 + t271 * t29 - t366 * t443;
t19 = t21 * t265;
t476 = t101 * t271;
t88 = t242 * t265 + t271 * t443;
t479 = t19 + (t365 - t476) * t88;
t411 = t265 * t28 + t365 * t448;
t478 = -t443 * t88 - t448 * t476 + t411;
t178 = t185 * pkin(9);
t422 = pkin(7) + pkin(8);
t220 = t422 * t274;
t203 = qJD(1) * t220;
t186 = t268 * t203;
t219 = t422 * t269;
t201 = qJD(1) * t219;
t408 = qJD(2) * pkin(2);
t192 = -t201 + t408;
t329 = t421 * t192 - t186;
t124 = t178 + t329;
t112 = pkin(3) * t261 + t124;
t190 = t421 * t203;
t307 = -t268 * t192 - t190;
t415 = t463 * pkin(9);
t125 = -t307 + t415;
t116 = t273 * t125;
t318 = -t112 * t267 - t116;
t419 = pkin(10) * t315;
t65 = -t318 + t419;
t405 = t266 * t65;
t143 = t314 * pkin(10);
t114 = t267 * t125;
t334 = t273 * t112 - t114;
t64 = t334 + t143;
t59 = pkin(4) * t256 + t64;
t38 = t272 * t59 - t405;
t36 = -pkin(5) * t242 - t38;
t477 = t36 * t101;
t264 = qJ(2) + qJ(3);
t259 = qJ(4) + t264;
t252 = qJ(5) + t259;
t237 = sin(t252);
t270 = sin(qJ(1));
t275 = cos(qJ(1));
t323 = g(1) * t275 + g(2) * t270;
t475 = t323 * t237;
t474 = t443 * t101;
t472 = -t101 ^ 2 + t443 ^ 2;
t440 = pkin(5) * t443;
t471 = -pkin(11) * t101 + t440;
t470 = -t101 * t242 + t29;
t251 = -pkin(2) * t274 - pkin(1);
t218 = t251 * qJD(1);
t162 = -pkin(3) * t463 + t218;
t108 = -pkin(4) * t315 + t162;
t229 = g(3) * t237;
t238 = cos(t252);
t364 = qJD(1) * qJD(2);
t351 = t274 * t364;
t160 = qJDD(2) * pkin(2) + t422 * (-t351 - t363);
t352 = t269 * t364;
t161 = t422 * (-t352 + t362);
t289 = qJD(3) * t307 + t421 * t160 - t268 * t161;
t62 = t260 * pkin(3) - t126 * pkin(9) + t289;
t353 = qJD(3) * t421;
t371 = qJD(3) * t268;
t287 = t268 * t160 + t421 * t161 + t192 * t353 - t203 * t371;
t69 = -t127 * pkin(9) + t287;
t295 = qJD(4) * t318 - t267 * t69 + t273 * t62;
t14 = t255 * pkin(4) - t72 * pkin(10) + t295;
t424 = (qJD(4) * t112 + t69) * t273 - t125 * t370 + t267 * t62;
t15 = pkin(10) * t292 + t424;
t425 = (qJD(5) * t59 + t15) * t272 + t266 * t14 - t65 * t368;
t469 = -t108 * t101 + t323 * t238 + t229 - t425;
t320 = t265 * t88 + t271 * t86;
t358 = t88 * t366;
t20 = t21 * t271;
t22 = qJD(6) * t88 - t271 * t239 + t265 * t29;
t360 = -t265 * t22 - t86 * t365 + t20;
t1 = t101 * t320 - t358 + t360;
t466 = t448 * t443;
t328 = t201 * t268 - t190;
t129 = t328 - t415;
t376 = -t421 * t201 - t186;
t130 = t178 + t376;
t250 = t421 * pkin(2) + pkin(3);
t387 = t267 * t268;
t465 = -t250 * t369 - (-t268 * t370 + (t421 * t273 - t387) * qJD(3)) * pkin(2) + t267 * t129 + t273 * t130;
t385 = t268 * t273;
t464 = -t250 * t370 + (-t268 * t369 + (-t421 * t267 - t385) * qJD(3)) * pkin(2) - t273 * t129 + t130 * t267;
t402 = t272 * t65;
t39 = t266 * t59 + t402;
t37 = pkin(11) * t242 + t39;
t50 = -pkin(5) * t101 - pkin(11) * t443 + t108;
t16 = -t265 * t37 + t271 * t50;
t462 = -t16 * t443 + t271 * t475 + t36 * t366;
t17 = t265 * t50 + t271 * t37;
t428 = qJD(5) * t39 - t272 * t14 + t266 * t15;
t4 = -t239 * pkin(5) + t428;
t416 = g(3) * t238;
t454 = t4 + t416;
t461 = t17 * t443 + t454 * t265 + t36 * t365;
t460 = t242 * t443 - t30;
t458 = -t108 * t443 - t416 - t428 + t475;
t456 = pkin(4) * t314;
t451 = -t143 - t465;
t450 = -t419 - t464;
t63 = t314 ^ 2 - t315 ^ 2;
t245 = sin(t259);
t246 = cos(t259);
t280 = -g(3) * t246 + t162 * t314 + t323 * t245 + t295;
t57 = -t256 * t314 + t292;
t247 = pkin(4) * t266 + pkin(11);
t439 = (qJD(6) * t247 - t456 + t471) * t448;
t249 = pkin(3) * t273 + pkin(4);
t386 = t267 * t272;
t374 = pkin(3) * t386 + t266 * t249;
t177 = pkin(11) + t374;
t420 = pkin(3) * t185;
t111 = -t420 - t456;
t52 = t111 + t471;
t438 = (qJD(6) * t177 + t52) * t448;
t325 = -pkin(2) * t387 + t273 * t250;
t176 = pkin(4) + t325;
t180 = pkin(2) * t385 + t250 * t267;
t377 = t266 * t176 + t272 * t180;
t137 = pkin(11) + t377;
t253 = pkin(2) * t372;
t437 = (qJD(6) * t137 + t253 + t52) * t448;
t436 = (t448 * pkin(11) + t440) * t448;
t394 = t314 * t315;
t327 = -t421 * t219 - t220 * t268;
t141 = -pkin(9) * t197 + t327;
t306 = -t268 * t269 + t421 * t274;
t375 = -t268 * t219 + t421 * t220;
t142 = pkin(9) * t306 + t375;
t378 = t267 * t141 + t273 * t142;
t56 = -t256 * t315 + t72;
t283 = g(3) * t245 - t162 * t315 + t323 * t246 - t424;
t156 = t197 * t267 - t273 * t306;
t157 = t197 * t273 + t267 * t306;
t103 = t272 * t156 + t157 * t266;
t104 = -t156 * t266 + t157 * t272;
t349 = t239 * pkin(11) + qJD(6) * t50 + t425;
t158 = t261 * t306;
t82 = -qJD(4) * t156 + t273 * t158 - t267 * t159;
t83 = qJD(4) * t157 + t267 * t158 + t273 * t159;
t46 = -qJD(5) * t103 - t266 * t83 + t272 * t82;
t330 = t273 * t141 - t142 * t267;
t79 = -pkin(10) * t157 + t330;
t80 = -pkin(10) * t156 + t378;
t49 = t266 * t79 + t272 * t80;
t356 = qJD(2) * t422;
t202 = t269 * t356;
t204 = t274 * t356;
t301 = -t421 * t202 - t268 * t204 - t219 * t353 - t220 * t371;
t92 = -pkin(9) * t159 + t301;
t288 = -t375 * qJD(3) + t268 * t202 - t421 * t204;
t93 = -t158 * pkin(9) + t288;
t305 = t141 * t369 - t142 * t370 + t267 * t93 + t273 * t92;
t31 = -pkin(10) * t83 + t305;
t294 = -t378 * qJD(4) - t267 * t92 + t273 * t93;
t32 = -t82 * pkin(10) + t294;
t48 = t266 * t80 - t272 * t79;
t5 = -qJD(5) * t48 + t266 * t32 + t272 * t31;
t166 = -pkin(3) * t306 + t251;
t117 = pkin(4) * t156 + t166;
t55 = pkin(5) * t103 - pkin(11) * t104 + t117;
t423 = -t103 * t349 + t4 * t104 - t49 * t28 - (qJD(6) * t55 + t5) * t448 + t36 * t46;
t414 = t55 * t28;
t316 = t176 * t272 - t180 * t266;
t410 = -qJD(5) * t316 + t450 * t266 - t451 * t272;
t409 = t377 * qJD(5) + t451 * t266 + t450 * t272;
t407 = t104 * t36;
t404 = t271 * t88;
t388 = t266 * t267;
t333 = -t124 * t267 - t116;
t70 = t333 - t419;
t380 = t273 * t124 - t114;
t71 = t143 + t380;
t400 = t266 * t70 + t272 * t71 - t249 * t367 - (-t267 * t368 + (t272 * t273 - t388) * qJD(4)) * pkin(3);
t399 = -t266 * t71 + t272 * t70 + t249 * t368 + (t267 * t367 + (t266 * t273 + t386) * qJD(4)) * pkin(3);
t393 = t185 * t463;
t258 = cos(t264);
t392 = t258 * t270;
t391 = t258 * t275;
t390 = t265 * t270;
t389 = t265 * t275;
t383 = t270 * t271;
t382 = t271 * t275;
t262 = t269 ^ 2;
t373 = -t274 ^ 2 + t262;
t254 = t269 * t408;
t148 = pkin(3) * t159 + t254;
t179 = pkin(2) * t352 + qJDD(1) * t251;
t107 = t127 * pkin(3) + t179;
t53 = -pkin(4) * t292 + t107;
t8 = t30 * pkin(5) - t29 * pkin(11) + t53;
t348 = qJD(6) * t37 - t8;
t40 = t266 * t64 + t402;
t324 = pkin(4) * t368 - t40;
t322 = g(1) * t270 - g(2) * t275;
t76 = pkin(4) * t83 + t148;
t312 = -t349 + t229;
t311 = -pkin(3) * t388 + t249 * t272;
t309 = -0.2e1 * pkin(1) * t364 - pkin(7) * qJDD(2);
t302 = -pkin(11) * t28 + t38 * t448 - t477;
t299 = -t137 * t28 + t410 * t448 - t477;
t298 = -t177 * t28 + t400 * t448 - t477;
t276 = qJD(2) ^ 2;
t297 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t276 + t322;
t277 = qJD(1) ^ 2;
t296 = pkin(1) * t277 - pkin(7) * qJDD(1) + t323;
t41 = t272 * t64 - t405;
t291 = -t477 - t247 * t28 + (-pkin(4) * t367 + t41) * t448;
t257 = sin(t264);
t279 = g(1) * t391 + g(2) * t392 + g(3) * t257 - t218 * t463 - t287;
t278 = -g(3) * t258 + t218 * t185 + t323 * t257 + t289;
t248 = -pkin(4) * t272 - pkin(5);
t175 = -pkin(5) - t311;
t171 = t238 * t382 + t390;
t170 = -t238 * t389 + t383;
t169 = -t238 * t383 + t389;
t168 = t238 * t390 + t382;
t163 = t253 - t420;
t136 = -pkin(5) - t316;
t128 = t185 ^ 2 - t463 ^ 2;
t109 = t111 + t253;
t106 = -t321 + (-qJD(1) * t197 - t185) * t261;
t105 = t126 - t481;
t47 = qJD(5) * t104 + t266 * t82 + t272 * t83;
t12 = pkin(5) * t47 - pkin(11) * t46 + t76;
t7 = t271 * t8;
t6 = qJD(5) * t49 + t266 * t31 - t272 * t32;
t2 = [qJDD(1), t322, t323, qJDD(1) * t262 + 0.2e1 * t269 * t351, 0.2e1 * t269 * t362 - 0.2e1 * t364 * t373, qJDD(2) * t269 + t274 * t276, qJDD(2) * t274 - t269 * t276, 0, t269 * t309 + t274 * t297, -t269 * t297 + t274 * t309, t126 * t197 - t158 * t185, t126 * t306 - t127 * t197 + t158 * t463 + t159 * t185, t158 * t261 + t197 * t260, -t159 * t261 + t260 * t306, 0, g(1) * t392 - g(2) * t391 + t251 * t127 + t218 * t159 - t179 * t306 - t254 * t463 + t260 * t327 + t261 * t288, t251 * t126 + t218 * t158 + t179 * t197 - t185 * t254 - t257 * t322 - t260 * t375 - t261 * t301, t157 * t72 - t314 * t82, -t156 * t72 + t157 * t292 + t314 * t83 + t315 * t82, t157 * t255 + t256 * t82, -t156 * t255 - t256 * t83, 0, t107 * t156 - t148 * t315 + t162 * t83 - t166 * t292 + t246 * t322 + t255 * t330 + t256 * t294, t107 * t157 - t148 * t314 + t162 * t82 + t166 * t72 - t245 * t322 - t255 * t378 - t256 * t305, t104 * t29 + t443 * t46, t101 * t46 - t103 * t29 - t104 * t30 - t443 * t47, t104 * t239 + t242 * t46, -t103 * t239 - t242 * t47, 0, -t101 * t76 + t53 * t103 + t108 * t47 + t117 * t30 + t238 * t322 - t48 * t239 - t6 * t242, t53 * t104 + t108 * t46 + t117 * t29 - t237 * t322 - t49 * t239 - t5 * t242 + t443 * t76, t46 * t404 + (-t358 + t20) * t104, -t320 * t46 + (-t19 - t22 * t271 + (t265 * t86 - t404) * qJD(6)) * t104, t271 * t46 * t448 + t21 * t103 + t88 * t47 + (-t366 * t448 + t26) * t104, -t22 * t103 - t104 * t411 - t449 * t46 - t86 * t47, t103 * t28 + t448 * t47, -g(1) * t169 - g(2) * t171 + t7 * t103 + t16 * t47 + t48 * t22 + t6 * t86 + (t12 * t448 + t414 + (-t103 * t37 - t448 * t49 + t407) * qJD(6)) * t271 + t423 * t265, -g(1) * t168 - g(2) * t170 - t17 * t47 + t48 * t21 + t6 * t88 + (-(-qJD(6) * t49 + t12) * t448 - t414 + t348 * t103 - qJD(6) * t407) * t265 + t423 * t271; 0, 0, 0, -t269 * t277 * t274, t373 * t277, t363, t362, qJDD(2), -g(3) * t274 + t269 * t296, g(3) * t269 + t274 * t296, t393, t128, t105, t106, t260, -t328 * t261 + (t421 * t260 - t261 * t371 + t372 * t463) * pkin(2) + t278, t376 * t261 + (t185 * t372 - t268 * t260 - t261 * t353) * pkin(2) + t279, t394, t63, t56, t57, t255, t163 * t315 + t325 * t255 + t464 * t256 + t280, t163 * t314 - t180 * t255 + t465 * t256 + t283, -t474, t472, t470, t460, t239, t101 * t109 + t239 * t316 - t242 * t409 + t458, -t109 * t443 - t239 * t377 + t242 * t410 + t469, t479, t1, t478, t480, -t466, t136 * t22 + t409 * t86 + (-t454 - t437) * t271 + t299 * t265 + t462, t136 * t21 + t409 * t88 + t299 * t271 + (-t475 + t437) * t265 + t461; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t393, t128, t105, t106, t260, -t261 * t307 + t278, t261 * t329 + t279, t394, t63, t56, t57, t255, -t333 * t256 + (-t185 * t315 + t255 * t273 - t256 * t370) * pkin(3) + t280, t380 * t256 + (-t185 * t314 - t255 * t267 - t256 * t369) * pkin(3) + t283, -t474, t472, t470, t460, t239, t101 * t111 + t239 * t311 - t242 * t399 + t458, -t111 * t443 - t239 * t374 + t242 * t400 + t469, t479, t1, t478, t480, -t466, t175 * t22 + t399 * t86 + (-t454 - t438) * t271 + t298 * t265 + t462, t175 * t21 + t399 * t88 + t298 * t271 + (-t475 + t438) * t265 + t461; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t394, t63, t56, t57, t255, -t256 * t318 + t280, t256 * t334 + t283, -t474, t472, t470, t460, t239, t242 * t40 + (-t101 * t314 + t239 * t272 - t242 * t368) * pkin(4) + t458, t41 * t242 + (-t239 * t266 - t242 * t367 + t314 * t443) * pkin(4) + t469, t479, t1, t478, t480, -t466, t248 * t22 + t324 * t86 + (-t454 - t439) * t271 + t291 * t265 + t462, t248 * t21 + t324 * t88 + t291 * t271 + (-t475 + t439) * t265 + t461; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t474, t472, t470, t460, t239, t39 * t242 + t458, t38 * t242 + t469, t479, -t449 * t88 + t476 * t86 + t360, t478, t480, -t466, -pkin(5) * t22 - t39 * t86 + t302 * t265 + (-t454 - t436) * t271 + t462, -pkin(5) * t21 - t39 * t88 + t302 * t271 + (-t475 + t436) * t265 + t461; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88 * t86, -t86 ^ 2 + t88 ^ 2, t448 * t86 + t21, t448 * t88 - t22, t28, -g(1) * t170 + g(2) * t168 + t17 * t448 + t265 * t312 - t36 * t88 - t365 * t37 + t7, g(1) * t171 - g(2) * t169 + t16 * t448 + t265 * t348 + t271 * t312 + t36 * t86;];
tau_reg  = t2;
