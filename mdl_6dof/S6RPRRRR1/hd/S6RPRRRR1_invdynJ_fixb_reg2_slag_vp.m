% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRRR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:55:24
% EndTime: 2019-03-09 06:55:40
% DurationCPUTime: 9.25s
% Computational Cost: add. (17801->624), mult. (38868->776), div. (0->0), fcn. (28064->18), ass. (0->319)
t257 = cos(qJ(6));
t357 = qJD(6) * t257;
t258 = cos(qJ(3));
t422 = cos(qJ(4));
t339 = t422 * t258;
t315 = qJD(1) * t339;
t255 = sin(qJ(3));
t420 = sin(qJ(4));
t337 = t420 * t255;
t176 = -qJD(1) * t337 + t315;
t184 = t255 * t422 + t258 * t420;
t177 = t184 * qJD(1);
t254 = sin(qJ(5));
t421 = cos(qJ(5));
t125 = -t421 * t176 + t177 * t254;
t442 = t125 * t257;
t445 = t357 + t442;
t253 = sin(qJ(6));
t291 = t254 * t176 + t177 * t421;
t245 = qJD(3) + qJD(4);
t327 = qJD(5) + t245;
t300 = t257 * t327;
t113 = t253 * t291 - t300;
t115 = t253 * t327 + t257 * t291;
t244 = qJDD(3) + qJDD(4);
t320 = qJDD(5) + t244;
t358 = qJD(6) * t253;
t328 = qJDD(1) * t420;
t329 = qJDD(1) * t422;
t334 = t420 * qJD(4);
t363 = qJD(1) * t255;
t111 = (qJD(3) * t420 + t334) * t363 - t255 * t329 - t258 * t328 - t245 * t315;
t440 = t245 * t184;
t112 = qJD(1) * t440 + t255 * t328 - t258 * t329;
t335 = qJD(5) * t421;
t360 = qJD(5) * t254;
t68 = t111 * t421 + t254 * t112 - t176 * t335 + t177 * t360;
t42 = -qJD(6) * t300 - t253 * t320 + t257 * t68 + t291 * t358;
t359 = qJD(6) * t115;
t43 = -t253 * t68 - t257 * t320 + t359;
t437 = qJD(6) + t125;
t444 = t437 * t253;
t6 = -t113 * t445 - t115 * t444 - t253 * t43 - t42 * t257;
t324 = t254 * t111 - t421 * t112;
t69 = qJD(5) * t291 - t324;
t67 = qJDD(6) + t69;
t14 = t113 * t291 + t257 * t67 - t437 * t444;
t39 = t42 * t253;
t16 = t115 * t445 - t39;
t15 = -t115 * t291 + t253 * t67 + t437 * t445;
t238 = t258 * qJD(2);
t251 = sin(pkin(11));
t221 = pkin(1) * t251 + pkin(7);
t203 = t221 * qJD(1);
t326 = pkin(8) * qJD(1) + t203;
t147 = -t255 * t326 + t238;
t145 = qJD(3) * pkin(3) + t147;
t362 = qJD(2) * t255;
t148 = t258 * t326 + t362;
t340 = t422 * t148;
t100 = t145 * t420 + t340;
t414 = t176 * pkin(9);
t91 = t100 + t414;
t402 = t254 * t91;
t172 = t177 * pkin(9);
t143 = t420 * t148;
t99 = t422 * t145 - t143;
t90 = -t172 + t99;
t83 = t245 * pkin(4) + t90;
t55 = t421 * t83 - t402;
t53 = -pkin(5) * t327 - t55;
t407 = t125 * t53;
t441 = t125 * t291;
t250 = qJ(3) + qJ(4);
t241 = qJ(5) + t250;
t226 = sin(t241);
t246 = qJ(1) + pkin(11);
t235 = cos(t246);
t378 = t226 * t235;
t234 = sin(t246);
t379 = t226 * t234;
t439 = g(1) * t378 + g(2) * t379;
t200 = t221 * qJDD(1);
t438 = qJD(2) * qJD(3) + t200;
t60 = -t125 ^ 2 + t291 ^ 2;
t84 = pkin(5) * t291 + pkin(10) * t125;
t46 = t125 * t327 - t68;
t237 = t258 * qJDD(2);
t102 = qJDD(3) * pkin(3) + t237 + (-pkin(8) * qJDD(1) - t200) * t255 - t148 * qJD(3);
t344 = t255 * qJDD(2) + t258 * t438;
t361 = qJD(3) * t255;
t122 = -t203 * t361 + t344;
t355 = qJD(1) * qJD(3);
t333 = t255 * t355;
t353 = t258 * qJDD(1);
t109 = (-t333 + t353) * pkin(8) + t122;
t50 = -qJD(4) * t100 + t422 * t102 - t109 * t420;
t35 = t244 * pkin(4) + t111 * pkin(9) + t50;
t336 = qJD(4) * t422;
t319 = -t420 * t102 - t422 * t109 - t145 * t336 + t148 * t334;
t37 = -pkin(9) * t112 - t319;
t10 = t254 * t35 + t83 * t335 - t91 * t360 + t421 * t37;
t252 = cos(pkin(11));
t222 = -t252 * pkin(1) - pkin(2);
t242 = t258 * pkin(3);
t430 = t222 - t242;
t178 = t430 * qJD(1);
t139 = -pkin(4) * t176 + t178;
t218 = g(3) * t226;
t227 = cos(t241);
t376 = t227 * t235;
t377 = t227 * t234;
t343 = g(1) * t376 + g(2) * t377 + t218;
t282 = t139 * t125 - t10 + t343;
t331 = t254 * t37 - t421 * t35;
t346 = t421 * t91;
t56 = t254 * t83 + t346;
t11 = -qJD(5) * t56 - t331;
t9 = -pkin(5) * t320 - t11;
t435 = t9 * t253 + t53 * t357;
t433 = t245 * t99;
t230 = pkin(3) * t422 + pkin(4);
t338 = t420 * t254;
t132 = t230 * t335 + (-qJD(5) * t338 + (t421 * t422 - t338) * qJD(4)) * pkin(3);
t103 = -t147 * t420 - t340;
t280 = t103 - t414;
t104 = t422 * t147 - t143;
t93 = -t172 + t104;
t66 = t254 * t280 + t421 * t93;
t397 = t132 - t66;
t316 = t421 * t420;
t396 = -t254 * t93 + t280 * t421 + t230 * t360 + (qJD(5) * t316 + (t254 * t422 + t316) * qJD(4)) * pkin(3);
t395 = pkin(1) * qJDD(1);
t391 = t437 * t291;
t310 = g(1) * t234 - g(2) * t235;
t431 = t310 * t226;
t366 = t227 * pkin(5) + t226 * pkin(10);
t54 = pkin(10) * t327 + t56;
t73 = pkin(5) * t125 - pkin(10) * t291 + t139;
t22 = -t253 * t54 + t257 * t73;
t23 = t253 * t73 + t257 * t54;
t429 = -t22 * t253 + t23 * t257;
t51 = t53 * t358;
t427 = -t291 * t22 + t257 * t439 + t51;
t417 = g(3) * t227;
t426 = t291 * t23 + t253 * t417 + t435;
t287 = t337 - t339;
t283 = t254 * t287;
t138 = t184 * t421 - t283;
t268 = t245 * t287;
t279 = t421 * t287;
t78 = qJD(5) * t279 + t184 * t360 + t254 * t440 + t268 * t421;
t401 = t257 * t78;
t292 = t138 * t358 + t401;
t383 = t138 * t257;
t424 = -t292 * t437 + t67 * t383;
t278 = -t139 * t291 - t331 - t417 + t439;
t47 = t245 * t291 + t324;
t423 = t138 * t320 - t327 * t78;
t260 = -pkin(8) - pkin(7);
t256 = sin(qJ(1));
t419 = pkin(1) * t256;
t418 = pkin(4) * t177;
t240 = cos(t250);
t416 = g(3) * t240;
t415 = g(3) * t258;
t154 = pkin(3) * t333 + qJDD(1) * t430;
t92 = pkin(4) * t112 + t154;
t13 = pkin(5) * t69 + pkin(10) * t68 + t92;
t8 = pkin(10) * t320 + t10;
t2 = qJD(6) * t22 + t13 * t253 + t257 * t8;
t1 = t2 * t257;
t413 = pkin(8) + t221;
t412 = t113 * t401 - t43 * t383;
t137 = t184 * t254 + t279;
t79 = -qJD(5) * t283 + t184 * t335 - t254 * t268 + t421 * t440;
t411 = t115 * t79 - t42 * t137;
t410 = t78 * t125 - t138 * t69;
t409 = t437 * t22;
t408 = t437 * t23;
t405 = t23 * t253;
t403 = t253 * t78;
t399 = -t184 * t112 - t268 * t176;
t394 = t113 * t253;
t393 = t115 * t113;
t392 = t115 * t257;
t384 = t138 * t253;
t382 = t177 * t176;
t381 = t203 * t255;
t380 = t203 * t258;
t239 = sin(t250);
t375 = t234 * t239;
t374 = t234 * t240;
t373 = t234 * t253;
t372 = t234 * t257;
t371 = t235 * t239;
t370 = t235 * t240;
t369 = t235 * t253;
t368 = t235 * t257;
t181 = t413 * t255;
t182 = t413 * t258;
t131 = -t420 * t181 + t422 * t182;
t225 = pkin(4) * t240;
t365 = t225 + t242;
t191 = pkin(2) + t365;
t259 = cos(qJ(1));
t243 = t259 * pkin(1);
t367 = t235 * t191 + t243;
t174 = pkin(3) * t316 + t254 * t230;
t248 = t255 ^ 2;
t249 = t258 ^ 2;
t364 = t248 - t249;
t204 = qJD(1) * t222;
t201 = qJDD(1) * t222;
t350 = pkin(10) * qJD(6) * t437;
t233 = pkin(3) * t361;
t349 = t115 * t403;
t262 = qJD(1) ^ 2;
t345 = t255 * t262 * t258;
t342 = t225 + t366;
t341 = -t9 - t417;
t192 = -pkin(3) * t255 - pkin(4) * t239;
t332 = -pkin(5) * t226 + t192;
t330 = qJD(3) * t413;
t321 = t1 - t343;
t318 = pkin(4) * t335;
t314 = t258 * t333;
t57 = t254 * t90 + t346;
t313 = pkin(4) * t360 - t57;
t312 = -pkin(10) * t67 + t407;
t311 = g(1) * t235 + g(2) * t234;
t309 = g(1) * t256 - g(2) * t259;
t130 = -t422 * t181 - t182 * t420;
t247 = -pkin(9) + t260;
t308 = -t235 * t247 - t419;
t307 = -t113 * t79 - t137 * t43;
t306 = -t125 * t55 + t291 * t56;
t167 = pkin(10) + t174;
t305 = -t167 * t67 + t407;
t228 = pkin(4) * t254 + pkin(10);
t304 = -t228 * t67 + t407;
t303 = t137 * t68 - t291 * t79;
t302 = t22 * t257 + t405;
t116 = -t184 * pkin(9) + t130;
t117 = -pkin(9) * t287 + t131;
t76 = t254 * t116 + t117 * t421;
t146 = pkin(4) * t287 + t430;
t81 = t137 * pkin(5) - t138 * pkin(10) + t146;
t45 = t253 * t81 + t257 * t76;
t44 = -t253 * t76 + t257 * t81;
t162 = t362 + t380;
t298 = -t125 * t405 - t22 * t442 + t321;
t297 = -qJD(6) * t73 + t218 - t8;
t295 = t311 * t226;
t294 = t311 * t239;
t293 = t138 * t357 - t403;
t80 = t418 + t84;
t289 = -qJD(1) * t204 + t311;
t170 = t255 * t330;
t171 = t258 * t330;
t288 = t170 * t420 - t171 * t422;
t86 = -t422 * t170 - t420 * t171 - t181 * t336 - t182 * t334;
t173 = -pkin(3) * t338 + t230 * t421;
t284 = 0.2e1 * qJD(3) * t204 - qJDD(3) * t221;
t281 = g(1) * t370 + g(2) * t374 + g(3) * t239 - t178 * t176 + t319;
t277 = t257 * t341 + t427;
t261 = qJD(3) ^ 2;
t276 = -t221 * t261 - 0.2e1 * t201 + t310;
t275 = -t293 * t437 - t384 * t67;
t12 = t257 * t13;
t3 = -qJD(6) * t23 - t253 * t8 + t12;
t274 = -qJD(6) * t302 - t3 * t253 + t1;
t185 = pkin(10) * t377;
t186 = pkin(10) * t376;
t273 = -g(1) * (-pkin(5) * t378 + t186) - g(2) * (-pkin(5) * t379 + t185);
t123 = -qJD(3) * t162 - t200 * t255 + t237;
t161 = t238 - t381;
t272 = t122 * t258 - t123 * t255 + (-t161 * t258 - t162 * t255) * qJD(3);
t271 = -t253 * t295 + t426;
t267 = g(1) * t371 + g(2) * t375 - t178 * t177 - t416 + t50;
t129 = pkin(4) * t440 + t233;
t266 = t184 * t244 - t245 * t268;
t265 = -t111 * t287 + t177 * t440;
t264 = pkin(9) * t268 + t181 * t334 - t182 * t336 + t288;
t232 = pkin(3) * t363;
t231 = t242 + pkin(2);
t229 = -pkin(4) * t421 - pkin(5);
t199 = qJDD(3) * t258 - t255 * t261;
t198 = qJDD(3) * t255 + t258 * t261;
t166 = -pkin(5) - t173;
t153 = t227 * t368 + t373;
t152 = -t227 * t369 + t372;
t151 = -t227 * t372 + t369;
t150 = t227 * t373 + t368;
t149 = t232 + t418;
t118 = -t176 ^ 2 + t177 ^ 2;
t110 = -t244 * t287 - t245 * t440;
t95 = t177 * t245 - t112;
t94 = -t176 * t245 - t111;
t87 = -qJD(4) * t131 + t288;
t77 = t232 + t80;
t75 = -t116 * t421 + t117 * t254;
t74 = -pkin(9) * t440 + t86;
t70 = -t137 * t320 - t327 * t79;
t58 = t421 * t90 - t402;
t41 = t43 * t257;
t32 = t79 * pkin(5) + t78 * pkin(10) + t129;
t29 = t253 * t84 + t257 * t55;
t28 = -t253 * t55 + t257 * t84;
t27 = t253 * t77 + t257 * t66;
t26 = -t253 * t66 + t257 * t77;
t25 = t253 * t80 + t257 * t58;
t24 = -t253 * t58 + t257 * t80;
t19 = qJD(5) * t76 + t254 * t74 - t264 * t421;
t18 = t116 * t335 - t117 * t360 + t254 * t264 + t421 * t74;
t17 = t394 * t437 - t41;
t5 = -qJD(6) * t45 - t18 * t253 + t257 * t32;
t4 = qJD(6) * t44 + t18 * t257 + t253 * t32;
t7 = [0, 0, 0, 0, 0, qJDD(1), t309, g(1) * t259 + g(2) * t256, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t252 * t395 + t310, -0.2e1 * t251 * t395 + t311, 0 (t309 + (t251 ^ 2 + t252 ^ 2) * t395) * pkin(1), qJDD(1) * t248 + 0.2e1 * t314, 0.2e1 * t255 * t353 - 0.2e1 * t355 * t364, t198, qJDD(1) * t249 - 0.2e1 * t314, t199, 0, t255 * t284 + t258 * t276, -t255 * t276 + t258 * t284 (t248 + t249) * t200 + t272 - t311, t201 * t222 - g(1) * (-pkin(2) * t234 + pkin(7) * t235 - t419) - g(2) * (pkin(2) * t235 + pkin(7) * t234 + t243) + t272 * t221, -t111 * t184 - t177 * t268, -t265 + t399, t266, t112 * t287 - t176 * t440, t110, 0, g(1) * t374 - g(2) * t370 + t112 * t430 + t130 * t244 + t154 * t287 - t176 * t233 + t178 * t440 + t87 * t245, -g(1) * t375 + g(2) * t371 - t111 * t430 - t131 * t244 + t154 * t184 + t177 * t233 - t178 * t268 - t86 * t245, -t100 * t440 + t130 * t111 - t131 * t112 + t86 * t176 - t87 * t177 - t50 * t184 - t311 + (t319 + t433) * t287, -t319 * t131 + t100 * t86 + t50 * t130 + t99 * t87 + t154 * t430 + t178 * t233 - g(1) * (-t231 * t234 - t235 * t260 - t419) - g(2) * (t231 * t235 - t234 * t260 + t243) -t138 * t68 - t291 * t78, t303 + t410, t423, t125 * t79 + t137 * t69, t70, 0, t129 * t125 + t92 * t137 + t139 * t79 + t146 * t69 - t19 * t327 + t227 * t310 - t320 * t75, t129 * t291 + t92 * t138 - t139 * t78 - t146 * t68 - t18 * t327 - t320 * t76 - t431, -t10 * t137 - t11 * t138 - t125 * t18 + t19 * t291 + t55 * t78 - t56 * t79 - t68 * t75 - t69 * t76 - t311, t10 * t76 + t56 * t18 - t11 * t75 - t55 * t19 + t92 * t146 + t139 * t129 - g(1) * (-t191 * t234 + t308) - g(2) * (-t234 * t247 + t367) -t115 * t292 - t383 * t42, t349 + (t39 + (-t392 + t394) * qJD(6)) * t138 + t412, t411 + t424, t113 * t293 + t384 * t43, t275 + t307, t137 * t67 + t437 * t79, -g(1) * t151 - g(2) * t153 + t113 * t19 + t137 * t3 + t138 * t435 + t22 * t79 - t53 * t403 + t43 * t75 + t437 * t5 + t44 * t67, -t53 * t401 - g(1) * t150 - g(2) * t152 + t115 * t19 - t437 * t4 - t137 * t2 - t23 * t79 - t42 * t75 - t45 * t67 + (t9 * t257 - t51) * t138, -t113 * t4 - t115 * t5 + t42 * t44 - t43 * t45 + t302 * t78 + t431 + (-qJD(6) * t429 - t2 * t253 - t257 * t3) * t138, t2 * t45 + t23 * t4 + t3 * t44 + t22 * t5 + t9 * t75 + t53 * t19 - g(1) * t308 - g(2) * (pkin(5) * t376 + pkin(10) * t378 + t367) + (-g(1) * (-t191 - t366) + g(2) * t247) * t234; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, t199, -t198, 0, t122 * t255 + t123 * t258 - g(3) + (-t161 * t255 + t162 * t258) * qJD(3), 0, 0, 0, 0, 0, 0, t110, -t266, t265 + t399, -t100 * t268 - t184 * t319 - t287 * t50 - t440 * t99 - g(3), 0, 0, 0, 0, 0, 0, t70, -t423, -t303 + t410, t10 * t138 - t11 * t137 - t55 * t79 - t56 * t78 - g(3), 0, 0, 0, 0, 0, 0, t275 - t307, t411 - t424, -t349 + (-t39 + (t392 + t394) * qJD(6)) * t138 + t412, t137 * t9 + t138 * t274 - t429 * t78 + t53 * t79 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t345, t364 * t262, t255 * qJDD(1), t345, t353, qJDD(3), -t415 + t237 + (t162 - t380) * qJD(3) + (t289 - t438) * t255, g(3) * t255 + (t161 + t381) * qJD(3) + t289 * t258 - t344, 0, 0, -t382, t118, t94, t382, t95, t244, -t103 * t245 + (t176 * t363 + t244 * t422 - t245 * t334) * pkin(3) + t267, t104 * t245 + (-t177 * t363 - t244 * t420 - t245 * t336) * pkin(3) + t281 (t100 + t103) * t177 + (-t104 + t99) * t176 + (t422 * t111 - t112 * t420 + (t176 * t422 + t177 * t420) * qJD(4)) * pkin(3), -t100 * t104 - t99 * t103 + (t50 * t422 - t319 * t420 - t415 + (t100 * t422 - t420 * t99) * qJD(4) + (-qJD(1) * t178 + t311) * t255) * pkin(3), t441, t60, t46, -t441, t47, t320, t173 * t320 - t149 * t125 + (-t56 - t396) * qJD(5) + t278 - t396 * t245, -t149 * t291 - t174 * t320 - t327 * t397 + t282, -t125 * t397 + t173 * t68 - t174 * t69 + t291 * t396 + t306, -g(3) * t365 + t10 * t174 + t11 * t173 - t139 * t149 - t192 * t311 - t396 * t55 + t397 * t56, t16, t6, t15, t17, t14, -t391, t166 * t43 + t305 * t253 + t396 * t113 + (-t132 * t253 - t167 * t357 - t26) * t437 + t277, -t166 * t42 + t305 * t257 + t396 * t115 + (-t257 * t132 + t167 * t358 + t27) * t437 + t271, t113 * t27 + t115 * t26 + (-t113 * t132 - t167 * t43 + (t115 * t167 - t22) * qJD(6)) * t257 + (t115 * t132 - t167 * t42 - t3 + (t113 * t167 - t23) * qJD(6)) * t253 + t298, t9 * t166 - t23 * t27 - t22 * t26 - g(1) * (t235 * t332 + t186) - g(2) * (t234 * t332 + t185) - g(3) * (t242 + t342) + t396 * t53 + t429 * t132 + t274 * t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t382, t118, t94, t382, t95, t244, t100 * t245 + t267, t281 + t433, 0, 0, t441, t60, t46, -t441, t47, t320, t57 * t245 + (t57 - t56) * qJD(5) + (-t177 * t125 + t320 * t421 - t327 * t360) * pkin(4) + t278, t58 * t327 + (-t177 * t291 - t254 * t320 - t327 * t335) * pkin(4) + t282, t58 * t125 - t57 * t291 + (t421 * t68 - t254 * t69 + (-t125 * t421 + t254 * t291) * qJD(5)) * pkin(4) + t306, t55 * t57 - t56 * t58 + (t421 * t11 - t416 + t10 * t254 - t139 * t177 + t294 + (-t254 * t55 + t421 * t56) * qJD(5)) * pkin(4), t16, t6, t15, t17, t14, -t391, t229 * t43 + t304 * t253 + t313 * t113 + (-t228 * t357 - t253 * t318 - t24) * t437 + t277, -t229 * t42 + t304 * t257 + t313 * t115 + (t228 * t358 - t257 * t318 + t25) * t437 + t271, t25 * t113 + t24 * t115 + (-t113 * t318 - t228 * t43 + (t115 * t228 - t22) * qJD(6)) * t257 + (t115 * t318 - t228 * t42 - t3 + (t113 * t228 - t23) * qJD(6)) * t253 + t298, t9 * t229 - t23 * t25 - t22 * t24 - t53 * t57 - g(3) * t342 + (t294 + (t254 * t53 + t421 * t429) * qJD(5)) * pkin(4) + t274 * t228 + t273; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t441, t60, t46, -t441, t47, t320, t245 * t56 + t278, t327 * t55 + t282, 0, 0, t16, t6, t15, t113 * t444 - t41, t14, -t391, -pkin(5) * t43 - t113 * t56 - t437 * t28 + t312 * t253 + (t341 - t350) * t257 + t427, pkin(5) * t42 - t115 * t56 + t437 * t29 + t312 * t257 + (-t295 + t350) * t253 + t426, t113 * t29 + t115 * t28 + (-t409 + (-t43 + t359) * pkin(10)) * t257 + (-t3 - t408 + (qJD(6) * t113 - t42) * pkin(10)) * t253 + t321, -t9 * pkin(5) + pkin(10) * t274 - g(3) * t366 - t22 * t28 - t23 * t29 - t53 * t56 + t273; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t393, -t113 ^ 2 + t115 ^ 2, t113 * t437 - t42, -t393, t115 * t437 - t43, t67, -g(1) * t152 + g(2) * t150 - t115 * t53 + t253 * t297 - t357 * t54 + t12 + t408, g(1) * t153 - g(2) * t151 + t113 * t53 + t409 + (qJD(6) * t54 - t13) * t253 + t297 * t257, 0, 0;];
tau_reg  = t7;
