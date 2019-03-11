% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% tau_reg [6x33]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRRR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:19:14
% EndTime: 2019-03-09 13:19:32
% DurationCPUTime: 8.78s
% Computational Cost: add. (11449->502), mult. (27376->663), div. (0->0), fcn. (22252->16), ass. (0->288)
t271 = sin(pkin(11));
t272 = cos(pkin(11));
t277 = sin(qJ(2));
t281 = cos(qJ(2));
t218 = -t271 * t277 + t272 * t281;
t206 = t218 * qJD(1);
t219 = t271 * t281 + t272 * t277;
t208 = t219 * qJD(1);
t276 = sin(qJ(4));
t410 = cos(qJ(4));
t155 = t410 * t206 - t208 * t276;
t412 = qJD(5) + qJD(6);
t460 = -t155 + t412;
t350 = -qJD(5) + t155;
t142 = qJD(6) - t350;
t275 = sin(qJ(5));
t279 = cos(qJ(6));
t274 = sin(qJ(6));
t280 = cos(qJ(5));
t363 = t274 * t280;
t226 = t275 * t279 + t363;
t225 = t274 * t275 - t279 * t280;
t455 = t460 * t225;
t349 = qJD(1) * qJD(2);
t332 = t281 * t349;
t333 = t277 * t349;
t164 = qJDD(1) * t219 - t271 * t333 + t272 * t332;
t207 = t219 * qJD(2);
t295 = qJD(1) * t207;
t285 = qJDD(1) * t218 - t295;
t303 = -t276 * t206 - t208 * t410;
t288 = qJD(4) * t303 - t276 * t164 + t410 * t285;
t73 = qJDD(5) - t288;
t70 = qJDD(6) + t73;
t459 = -t142 * t455 + t226 * t70;
t458 = t460 * t226;
t353 = qJD(5) * t280;
t431 = t155 * t280;
t457 = t353 - t431;
t354 = qJD(5) * t275;
t432 = t155 * t275;
t456 = t354 - t432;
t262 = qJ(2) + pkin(11) + qJ(4);
t254 = sin(t262);
t278 = sin(qJ(1));
t282 = cos(qJ(1));
t317 = g(1) * t282 + g(2) * t278;
t454 = t317 * t254;
t397 = qJ(3) + pkin(7);
t240 = t397 * t277;
t227 = qJD(1) * t240;
t241 = t397 * t281;
t228 = qJD(1) * t241;
t364 = t272 * t228;
t168 = t227 * t271 - t364;
t407 = pkin(8) * t206;
t139 = t168 - t407;
t211 = t271 * t228;
t169 = -t272 * t227 - t211;
t406 = pkin(8) * t208;
t140 = t169 - t406;
t256 = pkin(2) * t272 + pkin(3);
t409 = pkin(2) * t271;
t311 = t410 * t256 - t276 * t409;
t435 = t311 * qJD(4) - t276 * t139 - t140 * t410;
t109 = -pkin(4) * t303 - pkin(9) * t155;
t178 = t277 * qJD(1) * pkin(2) + pkin(3) * t208;
t94 = t109 + t178;
t453 = -t275 * t435 - t280 * t94;
t452 = pkin(10) * t432;
t267 = qJD(2) + qJD(4);
t376 = t303 * t267;
t451 = t288 - t376;
t319 = -t142 * t458 - t225 * t70;
t136 = -t280 * t267 - t275 * t303;
t138 = t267 * t275 - t280 * t303;
t88 = t279 * t136 + t138 * t274;
t392 = t303 * t88;
t450 = t319 - t392;
t449 = t456 * pkin(5);
t448 = -pkin(5) * t303 - pkin(10) * t431;
t266 = qJDD(2) + qJDD(4);
t393 = qJD(2) * pkin(2);
t215 = -t227 + t393;
t162 = t272 * t215 - t211;
t126 = qJD(2) * pkin(3) + t162 - t406;
t163 = t271 * t215 + t364;
t133 = t163 + t407;
t334 = qJD(4) * t410;
t355 = qJD(4) * t276;
t329 = qJD(2) * t397;
t204 = -t277 * qJD(3) - t281 * t329;
t161 = qJDD(2) * pkin(2) + qJD(1) * t204 - qJDD(1) * t240;
t203 = t281 * qJD(3) - t277 * t329;
t167 = qJD(1) * t203 + qJDD(1) * t241;
t112 = t272 * t161 - t167 * t271;
t77 = qJDD(2) * pkin(3) - pkin(8) * t164 + t112;
t113 = t271 * t161 + t272 * t167;
t82 = pkin(8) * t285 + t113;
t322 = -t126 * t355 - t133 * t334 - t276 * t82 + t410 * t77;
t19 = -pkin(4) * t266 - t322;
t74 = t410 * t164 + t206 * t334 - t208 * t355 + t276 * t285;
t325 = -t280 * t266 + t275 * t74;
t50 = qJD(5) * t138 + t325;
t10 = pkin(5) * t50 + t19;
t255 = cos(t262);
t270 = qJ(5) + qJ(6);
t264 = cos(t270);
t401 = g(3) * t264;
t78 = t126 * t410 - t276 * t133;
t67 = -t267 * pkin(4) - t78;
t54 = t136 * pkin(5) + t67;
t79 = t276 * t126 + t410 * t133;
t68 = pkin(9) * t267 + t79;
t258 = t281 * pkin(2) + pkin(1);
t234 = -qJD(1) * t258 + qJD(3);
t172 = -t206 * pkin(3) + t234;
t84 = -pkin(4) * t155 + pkin(9) * t303 + t172;
t35 = -t275 * t68 + t280 * t84;
t30 = -pkin(10) * t138 + t35;
t21 = -pkin(5) * t350 + t30;
t36 = t275 * t84 + t280 * t68;
t31 = -pkin(10) * t136 + t36;
t6 = t21 * t279 - t274 * t31;
t447 = t10 * t225 - t255 * t401 + t264 * t454 + t6 * t303 + t458 * t54;
t263 = sin(t270);
t402 = g(3) * t263;
t391 = t279 * t31;
t7 = t21 * t274 + t391;
t446 = t10 * t226 + t255 * t402 - t263 * t454 - t7 * t303 - t455 * t54;
t49 = t275 * t266 + t267 * t353 + t280 * t74 + t303 * t354;
t445 = -t136 * t457 - t275 * t50 + t49 * t280;
t47 = t49 * t275;
t444 = t138 * t457 + t47;
t351 = qJD(6) * t279;
t352 = qJD(6) * t274;
t14 = -t136 * t351 - t138 * t352 - t274 * t50 + t279 * t49;
t309 = t136 * t274 - t279 * t138;
t443 = t14 * t226 + t309 * t455;
t390 = t309 * t303;
t442 = -t390 + t459;
t378 = t138 * t303;
t65 = t275 * t73;
t441 = -t350 * t457 + t378 + t65;
t290 = qJD(6) * t309 - t274 * t49 - t279 * t50;
t440 = -t14 * t225 + t226 * t290 + t309 * t458 + t455 * t88;
t438 = t155 * t67;
t437 = t309 * t88;
t403 = g(3) * t255;
t436 = t19 + t403;
t373 = t155 * t267;
t434 = t74 - t373;
t379 = t136 * t303;
t433 = t142 * t303;
t428 = t303 * t155;
t427 = t350 * t303;
t425 = t309 ^ 2 - t88 ^ 2;
t424 = -t155 ^ 2 + t303 ^ 2;
t423 = t142 * t88 + t14;
t366 = t264 * t278;
t367 = t263 * t282;
t181 = -t255 * t366 + t367;
t365 = t264 * t282;
t368 = t263 * t278;
t183 = t255 * t365 + t368;
t27 = t31 * t352;
t420 = g(1) * t183 - g(2) * t181 + t254 * t401 + t54 * t88 + t27;
t180 = t255 * t368 + t365;
t182 = -t255 * t367 + t366;
t289 = t126 * t334 - t133 * t355 + t276 * t77 + t410 * t82;
t18 = t266 * pkin(9) + t289;
t186 = -pkin(3) * t218 - t258;
t346 = pkin(2) * t333 + qJDD(3);
t130 = pkin(3) * t295 + qJDD(1) * t186 + t346;
t28 = -pkin(4) * t288 - t74 * pkin(9) + t130;
t26 = t280 * t28;
t2 = t73 * pkin(5) - t49 * pkin(10) - qJD(5) * t36 - t275 * t18 + t26;
t305 = t280 * t18 + t275 * t28 + t84 * t353 - t354 * t68;
t3 = -pkin(10) * t50 + t305;
t340 = t279 * t2 - t274 * t3;
t419 = -g(1) * t182 + g(2) * t180 - qJD(6) * t7 + t254 * t402 + t54 * t309 + t340;
t63 = t67 * t354;
t418 = t280 * t454 + t35 * t303 + t63;
t417 = t275 * t436 - t36 * t303 + t67 * t353;
t416 = -t142 * t309 + t290;
t247 = g(3) * t254;
t415 = -t172 * t155 + t255 * t317 + t247 - t289;
t414 = t172 * t303 + t322 - t403 + t454;
t357 = t276 * t256 + t410 * t409;
t383 = t357 * qJD(4) + t139 * t410 - t276 * t140;
t166 = t276 * t218 + t219 * t410;
t110 = t226 * t166;
t413 = t275 * t94 - t280 * t435;
t411 = -pkin(9) - pkin(10);
t400 = g(3) * t281;
t399 = t280 * pkin(5);
t202 = pkin(9) + t357;
t398 = -pkin(10) - t202;
t173 = -t272 * t240 - t241 * t271;
t145 = -pkin(8) * t219 + t173;
t174 = -t271 * t240 + t272 * t241;
t146 = pkin(8) * t218 + t174;
t105 = t276 * t145 + t146 * t410;
t98 = t280 * t105;
t302 = t218 * t410 - t276 * t219;
t99 = -pkin(4) * t302 - pkin(9) * t166 + t186;
t394 = t275 * t99 + t98;
t385 = t383 + t449;
t382 = t275 * t109 + t280 * t78;
t210 = t218 * qJD(2);
t114 = qJD(4) * t302 - t276 * t207 + t210 * t410;
t381 = t114 * t275;
t380 = t114 * t280;
t377 = t138 * t275;
t370 = t166 * t275;
t369 = t166 * t280;
t362 = t275 * t278;
t361 = t275 * t282;
t360 = t278 * t280;
t358 = t280 * t282;
t144 = t272 * t203 + t271 * t204;
t268 = t277 ^ 2;
t356 = -t281 ^ 2 + t268;
t348 = t277 * qJDD(1);
t347 = t281 * qJDD(1);
t261 = t277 * t393;
t343 = qJD(5) * pkin(9) * t350;
t339 = qJD(5) * t411;
t337 = t166 * t354;
t336 = t166 * t353;
t179 = pkin(3) * t207 + t261;
t331 = qJD(6) * t21 + t3;
t328 = qJD(5) * t398;
t326 = -qJD(5) * t84 - t18;
t143 = -t203 * t271 + t272 * t204;
t323 = t350 * t275;
t321 = t449 - t79;
t320 = -t353 * t68 + t26;
t318 = -pkin(9) * t73 - t438;
t316 = g(1) * t278 - g(2) * t282;
t175 = t398 * t275;
t315 = -qJD(6) * t175 - t275 * t328 + t413 - t452;
t265 = t280 * pkin(10);
t176 = t202 * t280 + t265;
t314 = qJD(6) * t176 - t280 * t328 + t448 - t453;
t243 = t411 * t275;
t313 = -qJD(6) * t243 - t275 * t339 + t382 - t452;
t107 = t280 * t109;
t244 = pkin(9) * t280 + t265;
t312 = qJD(6) * t244 - t275 * t78 - t280 * t339 + t107 + t448;
t310 = -t202 * t73 - t438;
t201 = -pkin(4) - t311;
t66 = t280 * t73;
t308 = t350 * t456 + t66;
t306 = -0.2e1 * pkin(1) * t349 - pkin(7) * qJDD(2);
t304 = t145 * t410 - t276 * t146;
t301 = t336 + t381;
t300 = -t337 + t380;
t119 = -pkin(8) * t210 + t143;
t120 = -pkin(8) * t207 + t144;
t39 = qJD(4) * t304 + t276 * t119 + t120 * t410;
t115 = qJD(4) * t166 + t207 * t410 + t276 * t210;
t53 = pkin(4) * t115 - pkin(9) * t114 + t179;
t299 = -t105 * t354 + t275 * t53 + t280 * t39 + t99 * t353;
t296 = -qJDD(1) * t258 + t346;
t283 = qJD(2) ^ 2;
t293 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t283 + t316;
t284 = qJD(1) ^ 2;
t292 = pkin(1) * t284 - pkin(7) * qJDD(1) + t317;
t40 = qJD(4) * t105 - t119 * t410 + t276 * t120;
t259 = -pkin(4) - t399;
t200 = t255 * t358 + t362;
t199 = -t255 * t361 + t360;
t198 = -t255 * t360 + t361;
t197 = t255 * t362 + t358;
t184 = t201 - t399;
t111 = t225 * t166;
t97 = t280 * t99;
t62 = pkin(5) * t370 - t304;
t52 = t280 * t53;
t37 = -pkin(10) * t370 + t394;
t33 = -pkin(5) * t302 - pkin(10) * t369 - t105 * t275 + t97;
t23 = t114 * t363 - t274 * t337 - t352 * t370 + (t369 * t412 + t381) * t279;
t22 = -t110 * t412 - t225 * t114;
t20 = pkin(5) * t301 + t40;
t5 = -pkin(10) * t301 + t299;
t4 = -pkin(10) * t380 + t115 * pkin(5) - t275 * t39 + t52 + (-t98 + (pkin(10) * t166 - t99) * t275) * qJD(5);
t1 = [qJDD(1), t316, t317, qJDD(1) * t268 + 0.2e1 * t277 * t332, 0.2e1 * t277 * t347 - 0.2e1 * t349 * t356, qJDD(2) * t277 + t281 * t283, qJDD(2) * t281 - t277 * t283, 0, t277 * t306 + t281 * t293, -t277 * t293 + t281 * t306, -t112 * t219 + t113 * t218 - t143 * t208 + t144 * t206 - t162 * t210 - t163 * t207 - t173 * t164 + t174 * t285 - t317, t113 * t174 + t163 * t144 + t112 * t173 + t162 * t143 - t296 * t258 + t234 * t261 - g(1) * (-t258 * t278 + t282 * t397) - g(2) * (t258 * t282 + t278 * t397) -t114 * t303 + t166 * t74, t114 * t155 + t115 * t303 + t166 * t288 + t302 * t74, t114 * t267 + t166 * t266, -t115 * t267 + t266 * t302, 0, t172 * t115 - t130 * t302 - t155 * t179 - t186 * t288 + t255 * t316 + t266 * t304 - t40 * t267, -t105 * t266 + t172 * t114 + t130 * t166 - t179 * t303 + t186 * t74 - t254 * t316 - t39 * t267, t138 * t300 + t369 * t49 (-t136 * t280 - t377) * t114 + (-t47 - t280 * t50 + (t136 * t275 - t138 * t280) * qJD(5)) * t166, t138 * t115 - t300 * t350 - t302 * t49 + t369 * t73, -t136 * t115 + t301 * t350 + t302 * t50 - t370 * t73, -t115 * t350 - t302 * t73 -(-t105 * t353 + t52) * t350 + t97 * t73 - t320 * t302 + t35 * t115 + t40 * t136 - t304 * t50 + t67 * t336 - g(1) * t198 - g(2) * t200 + (-(-qJD(5) * t99 - t39) * t350 - t105 * t73 - t326 * t302 + t19 * t166 + t67 * t114) * t275, t299 * t350 - t394 * t73 + t305 * t302 - t36 * t115 + t40 * t138 - t304 * t49 + t67 * t380 - g(1) * t197 - g(2) * t199 + (t19 * t280 - t63) * t166, -t111 * t14 - t22 * t309, -t110 * t14 - t111 * t290 - t22 * t88 + t23 * t309, -t111 * t70 - t115 * t309 - t14 * t302 + t142 * t22, -t110 * t70 - t115 * t88 - t142 * t23 - t290 * t302, t115 * t142 - t302 * t70 (-t274 * t5 + t279 * t4) * t142 + (-t274 * t37 + t279 * t33) * t70 - t340 * t302 + t6 * t115 + t20 * t88 - t62 * t290 + t10 * t110 + t54 * t23 - g(1) * t181 - g(2) * t183 + ((-t274 * t33 - t279 * t37) * t142 + t7 * t302) * qJD(6), -g(1) * t180 - g(2) * t182 - t10 * t111 - t7 * t115 + t62 * t14 - t27 * t302 - t20 * t309 + t54 * t22 + (-(-qJD(6) * t37 + t4) * t142 - t33 * t70 + t2 * t302) * t274 + (-(qJD(6) * t33 + t5) * t142 - t37 * t70 + t331 * t302) * t279; 0, 0, 0, -t277 * t284 * t281, t356 * t284, t348, t347, qJDD(2), t277 * t292 - t400, g(3) * t277 + t281 * t292 (t163 + t168) * t208 + (-t169 + t162) * t206 + (-t272 * t164 + ((-t332 - t348) * t271 + (-t333 + t347) * t272) * t271) * pkin(2), -t162 * t168 - t163 * t169 + (-t400 + t112 * t272 + t113 * t271 + (-qJD(1) * t234 + t317) * t277) * pkin(2), t428, t424, t434, t451, t266, t155 * t178 + t266 * t311 - t267 * t383 + t414, t178 * t303 - t266 * t357 - t267 * t435 + t415, t444, t350 * t377 + t445, t441, t308 - t379, -t427, t201 * t50 - t436 * t280 + t310 * t275 + t383 * t136 - (-t202 * t353 + t453) * t350 + t418, t201 * t49 + t310 * t280 - t275 * t454 + t383 * t138 - (t202 * t354 + t413) * t350 + t417, t443, t440, t442, t450, t433 (t175 * t279 - t176 * t274) * t70 - t184 * t290 + t385 * t88 + (t274 * t315 - t279 * t314) * t142 + t447 -(t175 * t274 + t176 * t279) * t70 + t184 * t14 - t385 * t309 + (t274 * t314 + t279 * t315) * t142 + t446; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t206 ^ 2 - t208 ^ 2, t162 * t208 - t163 * t206 + t296 - t316, 0, 0, 0, 0, 0, -t288 - t376, t74 + t373, 0, 0, 0, 0, 0, t308 + t379, -t280 * t350 ^ 2 + t378 - t65, 0, 0, 0, 0, 0, t319 + t392, -t390 - t459; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t428, t424, t434, t451, t266, t267 * t79 + t414, t78 * t267 + t415, t444, t138 * t323 + t445, t441, -t323 * t350 - t379 + t66, -t427, -pkin(4) * t50 + t107 * t350 - t79 * t136 + (-t350 * t78 + t318) * t275 + (-t436 + t343) * t280 + t418, -pkin(4) * t49 - t382 * t350 - t79 * t138 + t318 * t280 + (-t454 - t343) * t275 + t417, t443, t440, t442, t450, t433 (t243 * t279 - t244 * t274) * t70 - t259 * t290 + t321 * t88 + (t274 * t313 - t279 * t312) * t142 + t447 -(t243 * t274 + t244 * t279) * t70 + t259 * t14 - t321 * t309 + (t274 * t312 + t279 * t313) * t142 + t446; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t138 * t136, -t136 ^ 2 + t138 ^ 2, -t136 * t350 + t49, -t325 + (-qJD(5) - t350) * t138, t73, -g(1) * t199 + g(2) * t197 - t67 * t138 - t36 * t350 + (t326 + t247) * t275 + t320, g(1) * t200 - g(2) * t198 + t136 * t67 + t247 * t280 - t35 * t350 - t305, -t437, t425, t423, t416, t70 -(-t274 * t30 - t391) * t142 + (-t138 * t88 - t142 * t352 + t279 * t70) * pkin(5) + t419 (-t142 * t31 - t2) * t274 + (t142 * t30 - t331) * t279 + (t138 * t309 - t142 * t351 - t274 * t70) * pkin(5) + t420; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t437, t425, t423, t416, t70, t7 * t142 + t419, t6 * t142 - t274 * t2 - t279 * t331 + t420;];
tau_reg  = t1;
