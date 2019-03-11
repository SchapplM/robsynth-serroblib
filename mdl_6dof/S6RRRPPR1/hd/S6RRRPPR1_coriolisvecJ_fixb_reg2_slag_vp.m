% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRPPR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:23:06
% EndTime: 2019-03-09 15:23:24
% DurationCPUTime: 7.53s
% Computational Cost: add. (19314->512), mult. (50247->680), div. (0->0), fcn. (38334->10), ass. (0->287)
t276 = qJD(2) + qJD(3);
t279 = sin(pkin(11));
t281 = cos(pkin(11));
t283 = sin(qJ(3));
t284 = sin(qJ(2));
t343 = qJD(1) * t284;
t333 = t283 * t343;
t286 = cos(qJ(3));
t287 = cos(qJ(2));
t355 = t286 * t287;
t224 = -qJD(1) * t355 + t333;
t242 = t283 * t287 + t284 * t286;
t226 = qJD(1) * t242;
t280 = sin(pkin(10));
t385 = cos(pkin(10));
t304 = -t280 * t224 + t226 * t385;
t165 = t276 * t279 + t281 * t304;
t282 = sin(qJ(6));
t285 = cos(qJ(6));
t173 = t279 * t304;
t403 = t276 * t281 - t173;
t309 = t285 * t403;
t113 = -t165 * t282 + t309;
t418 = t113 ^ 2;
t181 = t385 * t224 + t226 * t280;
t176 = qJD(6) + t181;
t417 = t113 * t176;
t112 = t285 * t165 + t282 * t403;
t416 = t112 ^ 2;
t356 = t285 * t281;
t401 = -t279 * t282 + t356;
t414 = t181 * t304;
t240 = t279 * t285 + t281 * t282;
t223 = t240 * qJD(6);
t413 = t240 * t181 + t223;
t402 = t401 * qJD(6);
t412 = t401 * t181 + t402;
t357 = t283 * t284;
t318 = t276 * t357;
t340 = qJD(1) * qJD(2);
t332 = t287 * t340;
t341 = qJD(3) * t286;
t336 = t287 * t341;
t345 = -qJD(1) * t336 - t286 * t332;
t187 = qJD(1) * t318 + t345;
t406 = t226 * t276;
t139 = -t187 * t385 - t280 * t406;
t379 = t139 * t279;
t411 = t165 * t181 + t379;
t358 = t281 * t139;
t410 = t181 * t403 + t358;
t393 = qJD(2) * pkin(2);
t272 = t284 * t393;
t409 = 0.2e1 * t272;
t408 = -0.2e1 * t340;
t398 = -pkin(8) - pkin(7);
t253 = t398 * t284;
t246 = qJD(1) * t253;
t235 = t246 + t393;
t337 = qJD(2) * t398;
t325 = qJD(1) * t337;
t236 = t284 * t325;
t237 = t287 * t325;
t254 = t398 * t287;
t248 = qJD(1) * t254;
t342 = qJD(3) * t283;
t141 = t235 * t341 + t286 * t236 + t283 * t237 + t248 * t342;
t87 = -qJ(4) * t406 - t224 * qJD(4) + t141;
t233 = t286 * t248;
t189 = t235 * t283 - t233;
t326 = -t283 * t236 + t286 * t237;
t142 = -qJD(3) * t189 + t326;
t88 = qJ(4) * t187 - qJD(4) * t226 + t142;
t50 = t280 * t88 + t385 * t87;
t47 = qJD(5) * t276 + t50;
t138 = -t187 * t280 + t385 * t406;
t271 = pkin(2) * t343;
t170 = pkin(3) * t406 + qJD(2) * t271;
t58 = t138 * pkin(4) - t139 * qJ(5) - qJD(5) * t304 + t170;
t23 = -t279 * t47 + t281 * t58;
t24 = t279 * t58 + t281 * t47;
t315 = -t23 * t279 + t24 * t281;
t407 = t139 * t240;
t405 = t304 * t403;
t241 = -t355 + t357;
t191 = t241 * t385 + t242 * t280;
t101 = t138 * t191;
t195 = -qJD(2) * t355 + t318 - t336;
t300 = t242 * qJD(3);
t196 = qJD(2) * t242 + t300;
t143 = -t195 * t280 + t196 * t385;
t404 = t143 * t181 + t101;
t229 = t283 * t248;
t188 = t286 * t235 + t229;
t219 = t226 * qJ(4);
t162 = t188 - t219;
t194 = t286 * t246 + t229;
t167 = -t219 + t194;
t193 = -t246 * t283 + t233;
t384 = qJ(4) * t224;
t305 = t193 + t384;
t329 = t385 * t283;
t346 = -t167 * t280 + t385 * t305 + (t280 * t286 + t329) * qJD(3) * pkin(2);
t201 = t283 * t253 - t286 * t254;
t42 = t282 * (qJD(6) * t165 + t379) - qJD(6) * t309 - t139 * t356;
t400 = t112 * t413 + t401 * t42;
t399 = t240 * t138 + t176 * t412;
t16 = -pkin(9) * t379 + t24;
t270 = -pkin(2) * t287 - pkin(1);
t252 = qJD(1) * t270;
t197 = t224 * pkin(3) + qJD(4) + t252;
t107 = pkin(4) * t181 - qJ(5) * t304 + t197;
t151 = pkin(3) * t276 + t162;
t163 = t189 - t384;
t330 = t385 * t163;
t99 = t280 * t151 + t330;
t95 = qJ(5) * t276 + t99;
t59 = t281 * t107 - t279 * t95;
t36 = pkin(5) * t181 - pkin(9) * t165 + t59;
t60 = t279 * t107 + t281 * t95;
t41 = pkin(9) * t403 + t60;
t313 = t282 * t41 - t285 * t36;
t9 = pkin(5) * t138 - pkin(9) * t358 + t23;
t3 = -qJD(6) * t313 + t285 * t16 + t282 * t9;
t397 = pkin(3) * t226;
t396 = t281 * pkin(5);
t273 = t281 * pkin(9);
t269 = pkin(2) * t286 + pkin(3);
t218 = pkin(2) * t329 + t280 * t269;
t209 = qJ(5) + t218;
t198 = (-pkin(9) - t209) * t279;
t199 = t209 * t281 + t273;
t146 = t198 * t282 + t199 * t285;
t262 = t280 * t283 * pkin(2);
t216 = t385 * pkin(2) * t341 - qJD(3) * t262;
t206 = qJD(5) + t216;
t366 = t181 * t281;
t323 = pkin(5) * t304 + pkin(9) * t366;
t118 = t167 * t385 + t280 * t305;
t122 = pkin(4) * t304 + qJ(5) * t181 + t397;
t119 = t122 + t271;
t69 = -t118 * t279 + t281 * t119;
t38 = t323 + t69;
t367 = t181 * t279;
t339 = pkin(9) * t367;
t70 = t281 * t118 + t279 * t119;
t52 = t339 + t70;
t395 = qJD(6) * t146 + t206 * t240 - t282 * t52 + t285 * t38;
t145 = t198 * t285 - t199 * t282;
t394 = -qJD(6) * t145 - t206 * t401 + t282 * t38 + t285 * t52;
t247 = t284 * t337;
t249 = t287 * t337;
t148 = t286 * t247 + t283 * t249 + t253 * t341 + t254 * t342;
t116 = -qJ(4) * t196 - qJD(4) * t241 + t148;
t149 = -qJD(3) * t201 - t283 * t247 + t286 * t249;
t291 = qJ(4) * t195 - qJD(4) * t242 + t149;
t66 = t116 * t385 + t280 * t291;
t144 = -t195 * t385 - t280 * t196;
t186 = pkin(3) * t196 + t272;
t192 = -t280 * t241 + t242 * t385;
t73 = pkin(4) * t143 - qJ(5) * t144 - qJD(5) * t192 + t186;
t29 = t279 * t73 + t281 * t66;
t200 = t286 * t253 + t254 * t283;
t177 = -qJ(4) * t242 + t200;
t178 = -qJ(4) * t241 + t201;
t130 = -t385 * t177 + t178 * t280;
t49 = t280 * t87 - t385 * t88;
t392 = t130 * t49;
t157 = t280 * t163;
t98 = t151 * t385 - t157;
t92 = -t276 * pkin(4) + qJD(5) - t98;
t391 = t181 * t92;
t265 = pkin(3) * t280 + qJ(5);
t227 = (-pkin(9) - t265) * t279;
t228 = t265 * t281 + t273;
t184 = t227 * t285 - t228 * t282;
t105 = t162 * t385 - t157;
t67 = -t105 * t279 + t281 * t122;
t37 = t323 + t67;
t68 = t281 * t105 + t279 * t122;
t51 = t339 + t68;
t388 = qJD(5) * t401 + qJD(6) * t184 - t282 * t37 - t285 * t51;
t185 = t227 * t282 + t228 * t285;
t387 = -qJD(5) * t240 - qJD(6) * t185 + t282 * t51 - t285 * t37;
t172 = pkin(5) * t367;
t386 = t172 + t346;
t383 = t113 * t304;
t382 = t112 * t113;
t381 = t112 * t304;
t380 = t139 * t192;
t377 = t144 * t279;
t376 = t144 * t281;
t375 = t165 * t304;
t374 = t165 * t279;
t373 = t176 * t304;
t372 = t304 ^ 2;
t371 = t304 * t276;
t369 = t181 ^ 2;
t368 = t181 * t276;
t364 = t192 * t279;
t363 = t226 * t224;
t362 = t226 * t252;
t289 = qJD(1) ^ 2;
t354 = t287 * t289;
t288 = qJD(2) ^ 2;
t353 = t288 * t284;
t352 = t288 * t287;
t351 = t118 - t216;
t204 = pkin(3) * t241 + t270;
t129 = pkin(4) * t191 - qJ(5) * t192 + t204;
t131 = t280 * t177 + t178 * t385;
t77 = t279 * t129 + t281 * t131;
t347 = t281 * t138 - t181 * t367;
t344 = t284 ^ 2 - t287 ^ 2;
t338 = t284 * t354;
t331 = t49 * t279 + t304 * t60;
t28 = -t279 * t66 + t281 * t73;
t328 = pkin(1) * t408;
t65 = t116 * t280 - t385 * t291;
t76 = t281 * t129 - t131 * t279;
t104 = t162 * t280 + t330;
t324 = t284 * t332;
t43 = qJD(6) * t112 + t407;
t322 = t113 * t412 - t240 * t43;
t35 = pkin(5) * t379 + t49;
t321 = t181 * t197 - t50;
t267 = -pkin(3) * t385 - pkin(4);
t320 = t138 * t401 - t413 * t176;
t317 = -t181 * t98 + t304 * t99;
t316 = -t49 * t281 - t304 * t59;
t314 = t279 * t59 - t281 * t60;
t13 = t282 * t36 + t285 * t41;
t53 = pkin(5) * t191 - t192 * t273 + t76;
t63 = -pkin(9) * t364 + t77;
t26 = -t282 * t63 + t285 * t53;
t27 = t282 * t53 + t285 * t63;
t217 = t269 * t385 - t262;
t312 = t130 * t139 + t192 * t49;
t311 = -t59 * t366 - t60 * t367 + t315;
t310 = t281 * t403;
t210 = -pkin(4) - t217;
t307 = -t197 * t304 - t49;
t306 = t138 * t279 + t181 * t366;
t80 = -pkin(5) * t403 + t92;
t303 = t313 * t304 - t35 * t401 + t413 * t80;
t302 = t13 * t304 + t35 * t240 + t412 * t80;
t301 = t224 * t252 - t141;
t299 = t144 * t92 + t312;
t298 = t138 * t192 + t139 * t191 + t144 * t181;
t4 = -qJD(6) * t13 - t282 * t16 + t285 * t9;
t297 = -t413 * t13 - t4 * t240 + t3 * t401 + t412 * t313;
t296 = t310 - t374;
t295 = t310 + t374;
t293 = -t138 * t209 + t139 * t210 - t181 * t206 + t391;
t292 = -qJD(5) * t181 - t138 * t265 + t139 * t267 + t391;
t275 = t281 ^ 2;
t274 = t279 ^ 2;
t250 = t267 - t396;
t203 = t210 - t396;
t202 = t271 + t397;
t168 = -t224 ^ 2 + t226 ^ 2;
t154 = -t345 + (t224 - t333) * t276;
t137 = t401 * t192;
t136 = t240 * t192;
t94 = t139 + t368;
t93 = -t138 + t371;
t89 = pkin(5) * t364 + t130;
t82 = t104 - t172;
t81 = -t369 + t372;
t79 = t410 * t279;
t78 = t411 * t281;
t62 = t144 * t240 + t192 * t402;
t61 = -t144 * t401 + t192 * t223;
t55 = t347 - t405;
t54 = t306 - t375;
t40 = pkin(5) * t377 + t65;
t34 = t296 * t181 + (-t274 + t275) * t139;
t31 = t320 - t383;
t30 = -t381 + t399;
t25 = -pkin(9) * t377 + t29;
t17 = pkin(5) * t143 - pkin(9) * t376 + t28;
t15 = -t113 * t413 - t401 * t43;
t14 = t112 * t412 - t42 * t240;
t7 = t322 - t400;
t6 = -qJD(6) * t27 + t285 * t17 - t282 * t25;
t5 = qJD(6) * t26 + t282 * t17 + t285 * t25;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t324, t344 * t408, t352, -0.2e1 * t324, -t353, 0, -pkin(7) * t352 + t284 * t328, pkin(7) * t353 + t287 * t328, 0, 0, -t187 * t242 - t195 * t226, t187 * t241 + t195 * t224 - t226 * t196 - t242 * t406, -t195 * t276, t224 * t196 + t241 * t406, -t196 * t276, 0, t224 * t272 + t149 * t276 + t252 * t196 + (t270 * t300 + (t284 * pkin(2) * t241 + t242 * t270) * qJD(2)) * qJD(1), -t148 * t276 - t187 * t270 - t195 * t252 + t226 * t409, -t141 * t241 - t142 * t242 - t148 * t224 - t149 * t226 + t200 * t187 + t188 * t195 - t189 * t196 - t201 * t406, t141 * t201 + t142 * t200 + t148 * t189 + t149 * t188 + t252 * t409, t144 * t304 + t380, -t143 * t304 - t298, t144 * t276, t404, -t143 * t276, 0, t138 * t204 + t143 * t197 + t170 * t191 + t181 * t186 - t276 * t65, t139 * t204 + t144 * t197 + t170 * t192 + t186 * t304 - t276 * t66, -t131 * t138 - t143 * t99 - t144 * t98 - t181 * t66 - t191 * t50 + t304 * t65 + t312, t131 * t50 + t170 * t204 + t186 * t197 - t65 * t98 + t66 * t99 + t392, t165 * t376 + t275 * t380, t144 * t296 - 0.2e1 * t358 * t364, t165 * t143 + t281 * t298, t274 * t380 - t377 * t403, t143 * t403 - t279 * t298, t404, t76 * t138 + t59 * t143 + t28 * t181 + t23 * t191 + t279 * t299 - t403 * t65, -t138 * t77 - t143 * t60 + t165 * t65 - t181 * t29 - t191 * t24 + t281 * t299, -t28 * t165 - t29 * t173 + (-t139 * t76 - t144 * t59 - t192 * t23 + t276 * t29) * t281 + (-t139 * t77 - t144 * t60 - t192 * t24) * t279, t23 * t76 + t24 * t77 + t28 * t59 + t29 * t60 + t65 * t92 + t392, -t112 * t61 - t137 * t42, -t112 * t62 - t113 * t61 + t136 * t42 - t137 * t43, t112 * t143 + t137 * t138 - t176 * t61 - t191 * t42, -t113 * t62 + t136 * t43, t113 * t143 - t136 * t138 - t176 * t62 - t191 * t43, t143 * t176 + t101, -t113 * t40 + t136 * t35 + t138 * t26 - t143 * t313 + t176 * t6 + t191 * t4 + t43 * t89 + t62 * t80, t112 * t40 - t13 * t143 + t137 * t35 - t138 * t27 - t176 * t5 - t191 * t3 - t42 * t89 - t61 * t80, -t112 * t6 + t113 * t5 - t13 * t62 - t136 * t3 - t137 * t4 + t26 * t42 - t27 * t43 - t313 * t61, t13 * t5 + t26 * t4 + t27 * t3 - t313 * t6 + t35 * t89 + t40 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t338, t344 * t289, 0, t338, 0, 0, t289 * pkin(1) * t284, pkin(1) * t354, 0, 0, t363, t168, t154, -t363, 0, 0, -t224 * t271 - t193 * t276 - t362 + (t233 + (-pkin(2) * t276 - t235) * t283) * qJD(3) + t326, t194 * t276 + (-t226 * t343 - t276 * t341) * pkin(2) + t301 (t189 + t193) * t226 + (-t188 + t194) * t224 + (t286 * t187 + (-t286 * t224 + t283 * t226) * qJD(3) - t283 * t406) * pkin(2), -t188 * t193 - t189 * t194 + (-t252 * t343 + t141 * t283 + t142 * t286 + (-t188 * t283 + t189 * t286) * qJD(3)) * pkin(2), t414, t81, t94, -t414, t93, 0, -t181 * t202 - t276 * t346 + t307, -t202 * t304 + t276 * t351 + t321, -t138 * t218 - t139 * t217 + t181 * t351 + t304 * t346 + t317, -t197 * t202 - t217 * t49 + t218 * t50 - t346 * t98 - t351 * t99, t78, t34, t54, -t79, t55, -t414, -t69 * t181 + t279 * t293 - t346 * t403 + t316, t165 * t346 + t181 * t70 + t281 * t293 + t331, t206 * t310 - t70 * t403 + (t206 * t279 + t69) * t165 + t311, -t206 * t314 + t209 * t315 + t210 * t49 + t346 * t92 - t59 * t69 - t60 * t70, t14, t7, t30, t15, t31, -t373, -t113 * t386 + t138 * t145 - t176 * t395 + t203 * t43 + t303, t112 * t386 - t138 * t146 + t176 * t394 - t203 * t42 + t302, t112 * t395 - t113 * t394 + t145 * t42 - t146 * t43 + t297, -t13 * t394 + t145 * t4 + t146 * t3 + t203 * t35 + t313 * t395 + t386 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t363, t168, t154, -t363, 0, 0, t189 * t276 + t142 - t362, t188 * t276 + t301, 0, 0, t414, t81, t94, -t414, t93, 0, t104 * t276 - t181 * t397 + t307, t105 * t276 - t304 * t397 + t321, -t104 * t304 + t105 * t181 + (-t138 * t280 - t139 * t385) * pkin(3) + t317, t98 * t104 - t99 * t105 + (-t197 * t226 + t280 * t50 - t385 * t49) * pkin(3), t78, t34, t54, -t79, t55, -t414, t104 * t403 - t67 * t181 + t279 * t292 + t316, -t104 * t165 + t181 * t68 + t281 * t292 + t331, qJD(5) * t295 + t67 * t165 - t403 * t68 + t311, -qJD(5) * t314 - t104 * t92 + t265 * t315 + t267 * t49 - t59 * t67 - t60 * t68, t14, t7, t30, t15, t31, -t373, t113 * t82 + t138 * t184 + t176 * t387 + t250 * t43 + t303, -t112 * t82 - t138 * t185 - t176 * t388 - t250 * t42 + t302, -t112 * t387 + t113 * t388 + t184 * t42 - t185 * t43 + t297, t13 * t388 + t184 * t4 + t185 * t3 + t250 * t35 - t313 * t387 - t80 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t138 + t371, t139 - t368, -t369 - t372, t181 * t99 + t304 * t98 + t170, 0, 0, 0, 0, 0, 0, t347 + t405, -t306 - t375, t295 * t181 + (-t274 - t275) * t139, -t181 * t314 + t23 * t281 + t24 * t279 - t304 * t92, 0, 0, 0, 0, 0, 0, t320 + t383, -t381 - t399, t322 + t400, t13 * t412 + t240 * t3 - t304 * t80 + t313 * t413 + t4 * t401; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t411, t410, -t165 ^ 2 - t403 ^ 2, t165 * t59 - t403 * t60 + t49, 0, 0, 0, 0, 0, 0, t112 * t176 + t43, -t42 + t417, -t416 - t418, -t112 * t313 - t113 * t13 + t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t382, t416 - t418, -t42 - t417, t382, -t407 + (-qJD(6) + t176) * t112, t138, -t112 * t80 + t13 * t176 + t4, -t113 * t80 - t176 * t313 - t3, 0, 0;];
tauc_reg  = t1;
