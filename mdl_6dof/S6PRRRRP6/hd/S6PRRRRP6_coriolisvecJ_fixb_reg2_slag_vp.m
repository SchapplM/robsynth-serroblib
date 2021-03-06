% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6PRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRRRP6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:34:26
% EndTime: 2019-03-09 00:34:53
% DurationCPUTime: 11.60s
% Computational Cost: add. (12704->623), mult. (34516->844), div. (0->0), fcn. (28050->12), ass. (0->288)
t227 = cos(qJ(3));
t369 = cos(pkin(7));
t311 = t227 * t369;
t218 = sin(pkin(7));
t223 = sin(qJ(3));
t356 = t218 * t223;
t183 = pkin(2) * t311 - pkin(9) * t356;
t175 = qJD(3) * t183;
t224 = sin(qJ(2));
t228 = cos(qJ(2));
t312 = t223 * t369;
t246 = -t224 * t312 + t227 * t228;
t219 = sin(pkin(6));
t345 = qJD(1) * t219;
t348 = t246 * t345 - t175;
t283 = pkin(3) * t223 - pkin(10) * t227;
t257 = t283 * qJD(3);
t322 = t224 * t345;
t427 = (t257 - t322) * t218;
t226 = cos(qJ(4));
t394 = pkin(10) * t226;
t355 = t218 * t227;
t184 = pkin(2) * t312 + pkin(9) * t355;
t170 = pkin(10) * t369 + t184;
t284 = -pkin(3) * t227 - pkin(10) * t223;
t171 = (-pkin(2) + t284) * t218;
t222 = sin(qJ(4));
t335 = qJD(4) * t226;
t337 = qJD(4) * t222;
t412 = -t170 * t337 + t171 * t335 + t427 * t222 - t348 * t226;
t176 = qJD(3) * t184;
t244 = t223 * t228 + t224 * t311;
t347 = -t244 * t345 + t176;
t426 = pkin(10) * t337;
t339 = qJD(3) * t223;
t319 = t218 * t339;
t425 = -pkin(11) * t319 - t412;
t180 = t222 * t356 - t226 * t369;
t338 = qJD(3) * t227;
t318 = t218 * t338;
t140 = qJD(4) * t180 - t226 * t318;
t181 = t222 * t369 + t226 * t356;
t141 = qJD(4) * t181 + t222 * t318;
t424 = t141 * pkin(4) + t140 * pkin(11) + t347;
t220 = cos(pkin(6));
t387 = qJD(2) * pkin(2);
t197 = t228 * t345 + t387;
t309 = t369 * t197;
t343 = qJD(2) * t218;
t187 = pkin(9) * t343 + t322;
t350 = t227 * t187;
t247 = t223 * t309 + t350;
t282 = pkin(4) * t222 - pkin(11) * t226;
t340 = qJD(2) * t227;
t344 = qJD(1) * t223;
t423 = -qJD(5) * t394 + t282 * qJD(4) - (t220 * t344 + t282 * t340) * t218 - t247;
t411 = -t170 * t335 - t171 * t337 + t348 * t222 + t427 * t226;
t208 = t218 * t340;
t422 = t222 * t208 - t337;
t307 = t369 * qJD(2);
t262 = t307 + qJD(3);
t341 = qJD(2) * t223;
t321 = t218 * t341;
t295 = t222 * t321;
t159 = -t226 * t262 + t295;
t155 = qJD(5) + t159;
t221 = sin(qJ(5));
t225 = cos(qJ(5));
t276 = t208 - qJD(4);
t357 = t218 * t220;
t242 = qJD(1) * t357 + t309;
t119 = t223 * t242 + t350;
t104 = pkin(10) * t262 + t119;
t313 = t220 * t369;
t207 = qJD(1) * t313;
t126 = t207 + (qJD(2) * t284 - t197) * t218;
t63 = t226 * t104 + t222 * t126;
t56 = -pkin(11) * t276 + t63;
t177 = t223 * t187;
t118 = t227 * t242 - t177;
t103 = -pkin(3) * t262 - t118;
t161 = t222 * t262 + t226 * t321;
t65 = t159 * pkin(4) - t161 * pkin(11) + t103;
t25 = t221 * t65 + t225 * t56;
t136 = (t257 + t322) * t343;
t241 = t246 * qJD(2);
t290 = t220 * t318;
t84 = (t227 * t309 - t177) * qJD(3) + (t219 * t241 + t290) * qJD(1);
t251 = t104 * t337 - t126 * t335 - t222 * t136 - t226 * t84;
t329 = qJD(2) * qJD(3);
t314 = t218 * t329;
t289 = t223 * t314;
t22 = pkin(11) * t289 - t251;
t333 = qJD(5) * t225;
t334 = qJD(5) * t221;
t288 = t227 * t314;
t263 = t222 * t288;
t330 = t161 * qJD(4);
t132 = t263 + t330;
t249 = t307 * t322;
t285 = qJD(3) * t309;
t291 = t220 * t319;
t342 = qJD(2) * t219;
t320 = t228 * t342;
t403 = qJD(4) * t295 - t226 * (qJD(4) * t262 + t288);
t45 = t132 * pkin(4) + pkin(11) * t403 + qJD(1) * t291 + t187 * t338 + t223 * t285 + t227 * t249 + t320 * t344;
t308 = t221 * t22 - t225 * t45 + t56 * t333 + t65 * t334;
t259 = t155 * t25 - t308;
t421 = qJD(5) * t276 + t403;
t420 = -t222 * t403 + t226 * t330;
t128 = t161 * t221 + t225 * t276;
t130 = t225 * t161 - t221 * t276;
t293 = t226 * t208;
t152 = t221 * t293 - t225 * t321;
t352 = t225 * t227;
t153 = (t221 * t223 + t226 * t352) * t343;
t73 = t161 * t333 - t221 * t421 - t225 * t289;
t379 = t225 * t73;
t72 = t161 * t334 - t221 * t289 + t225 * t421;
t381 = t221 * t72;
t419 = (qJD(5) * (t128 * t221 - t130 * t225) - t379 + t381) * t222 - (t128 * t225 + t130 * t221) * t335 + t128 * t153 + t130 * t152;
t281 = t221 * t335 - t152;
t364 = t132 * t221;
t418 = -t281 * t155 + (t128 * t276 - t155 * t333 - t364) * t222 + t226 * t73;
t305 = t128 * t155;
t366 = t130 * t155;
t417 = (t73 + t366) * t221 + (t72 + t305) * t225;
t169 = -pkin(3) * t369 - t183;
t107 = t180 * pkin(4) - t181 * pkin(11) + t169;
t121 = t226 * t170 + t222 * t171;
t109 = -pkin(11) * t355 + t121;
t389 = t107 * t333 - t109 * t334 + t221 * t424 - t225 * t425;
t18 = qJ(6) * t155 + t25;
t395 = pkin(5) * t132;
t2 = t308 - t395;
t416 = -t155 * t18 + t2;
t204 = -pkin(4) * t226 - pkin(11) * t222 - pkin(3);
t336 = qJD(4) * t225;
t173 = t283 * t343;
t92 = t226 * t118 + t222 * t173;
t82 = pkin(11) * t321 + t92;
t374 = -t222 * t336 * pkin(10) + t204 * t333 + t221 * t423 - t225 * t82;
t372 = -pkin(4) * t319 - t411;
t415 = t92 + t426;
t410 = t223 * t227;
t409 = t221 * t107 + t225 * t109;
t310 = t228 * t369;
t408 = -t223 * t224 + t227 * t310;
t407 = t293 - t335;
t62 = -t222 * t104 + t226 * t126;
t55 = pkin(4) * t276 - t62;
t29 = t128 * pkin(5) - t130 * qJ(6) + t55;
t392 = pkin(11) * t132;
t406 = t155 * t29 - t392;
t363 = t132 * t225;
t398 = t155 ^ 2;
t402 = -t128 * t161 + t221 * t398 - t363;
t373 = -t204 * t334 + (t82 + t426) * t221 + t423 * t225;
t388 = -qJD(5) * t409 + t221 * t425 + t225 * t424;
t399 = t130 ^ 2;
t229 = qJD(2) ^ 2;
t397 = qJ(6) * t141 + qJD(6) * t180 + t389;
t396 = pkin(5) * t141 + t388;
t393 = pkin(11) * t130;
t27 = -t104 * t335 - t126 * t337 + t226 * t136 - t222 * t84;
t23 = -pkin(4) * t289 - t27;
t5 = pkin(5) * t73 + qJ(6) * t72 - qJD(6) * t130 + t23;
t391 = t221 * t5;
t390 = t225 * t5;
t386 = t130 * t29;
t138 = -t408 * t219 - t220 * t355;
t240 = t244 * qJD(2);
t85 = t247 * qJD(3) + (t219 * t240 + t291) * qJD(1);
t385 = t138 * t85;
t382 = t221 * t23;
t380 = t225 * t23;
t377 = -qJ(6) * t422 - qJD(6) * t226 + t374;
t376 = pkin(5) * t422 - t373;
t277 = pkin(5) * t221 - qJ(6) * t225;
t260 = pkin(10) + t277;
t278 = pkin(5) * t225 + qJ(6) * t221;
t91 = -t222 * t118 + t173 * t226;
t81 = -pkin(4) * t321 - t91;
t375 = -pkin(5) * t152 + qJ(6) * t153 + (qJD(5) * t278 - qJD(6) * t225) * t222 + t260 * t335 - t81;
t323 = t221 * t355;
t143 = t181 * t225 - t323;
t142 = t181 * t221 + t218 * t352;
t86 = qJD(5) * t142 + t225 * t140 - t221 * t319;
t87 = -qJD(5) * t323 - t140 * t221 + t181 * t333 - t225 * t319;
t371 = -pkin(5) * t87 - qJ(6) * t86 + qJD(6) * t143 - t372;
t370 = -qJD(6) * t221 + t155 * t277 - t63;
t117 = pkin(4) * t161 + pkin(11) * t159;
t38 = t221 * t117 + t225 * t62;
t368 = qJ(6) * t132;
t367 = t130 * t128;
t365 = t132 * t180;
t154 = -t197 * t218 + t207;
t362 = t154 * t218;
t361 = t155 * t161;
t360 = t161 * t159;
t292 = t161 * t208;
t359 = t204 * t225;
t215 = t218 ^ 2;
t358 = t215 * t229;
t354 = t219 * t229;
t351 = t226 * t132;
t24 = -t221 * t56 + t225 * t65;
t349 = qJD(6) - t24;
t167 = t221 * t204 + t225 * t394;
t346 = t223 ^ 2 - t227 ^ 2;
t331 = t103 * qJD(4);
t328 = t215 * t387;
t326 = pkin(11) * t334;
t325 = pkin(11) * t333;
t324 = t224 * t354;
t316 = t128 ^ 2 - t399;
t315 = t215 * t329;
t120 = -t222 * t170 + t171 * t226;
t306 = t227 * t276;
t304 = t155 * t225;
t303 = t276 * t218;
t302 = qJD(4) * t276;
t300 = t215 * t324;
t299 = t358 * t410;
t296 = t218 * t224 * t342;
t286 = t218 * t229 * t369;
t108 = pkin(4) * t355 - t120;
t280 = t225 * t335 - t153;
t279 = (qJD(5) * t128 - t72) * pkin(11);
t275 = t128 * t87 + t142 * t73;
t17 = -pkin(5) * t155 + t349;
t274 = t17 * t225 - t18 * t221;
t273 = -t221 * t25 - t225 * t24;
t270 = -t27 * t222 - t226 * t251;
t37 = t117 * t225 - t221 * t62;
t245 = t223 * t310 + t224 * t227;
t139 = t219 * t245 + t220 * t356;
t179 = -t219 * t228 * t218 + t313;
t106 = t139 * t226 + t179 * t222;
t71 = t106 * t225 + t138 * t221;
t70 = t106 * t221 - t138 * t225;
t53 = t107 * t225 - t109 * t221;
t105 = t139 * t222 - t179 * t226;
t264 = t315 * t410;
t261 = 0.2e1 * t307 + qJD(3);
t258 = -t328 + t362;
t3 = t225 * t22 + t221 * t45 + t65 * t333 - t334 * t56;
t254 = t276 * t222;
t253 = t155 * t55 - t392;
t99 = t290 + (t408 * qJD(3) + t241) * t219;
t50 = -qJD(4) * t105 + t222 * t296 + t226 * t99;
t98 = t291 + (qJD(3) * t245 + t240) * t219;
t13 = qJD(5) * t71 + t221 * t50 - t98 * t225;
t14 = -qJD(5) * t70 + t221 * t98 + t225 * t50;
t252 = -t14 * t128 + t13 * t130 - t70 * t72 - t71 * t73;
t49 = qJD(4) * t106 + t222 * t99 - t226 * t296;
t250 = t105 * t73 + t49 * t128 - t13 * t155 - t132 * t70;
t243 = t221 * t305 - t379;
t239 = t105 * t72 - t130 * t49 + t132 * t71 + t14 * t155;
t238 = t128 * t86 - t130 * t87 + t142 * t72 - t143 * t73;
t237 = qJD(3) * t187 + t249;
t236 = t128 * t141 + t132 * t142 + t155 * t87 + t180 * t73;
t234 = t221 * t222 * t73 + (t222 * t333 + t281) * t128;
t231 = -t285 - t154 * t343 + (-qJD(3) * t357 - t320) * qJD(1);
t230 = t73 - t366;
t199 = -pkin(4) - t278;
t172 = t260 * t222;
t166 = -t221 * t394 + t359;
t146 = -t359 + (pkin(10) * t221 + pkin(5)) * t226;
t145 = -qJ(6) * t226 + t167;
t78 = pkin(5) * t130 + qJ(6) * t128;
t69 = pkin(11) * t379;
t66 = -t155 * t254 - t351;
t61 = t141 * t155 + t365;
t57 = pkin(5) * t142 - qJ(6) * t143 + t108;
t48 = -pkin(5) * t180 - t53;
t47 = qJ(6) * t180 + t409;
t44 = -t72 + t305;
t33 = -t130 * t161 + t155 * t304 + t364;
t32 = -pkin(5) * t161 - t37;
t31 = qJ(6) * t161 + t38;
t28 = t130 * t304 - t381;
t19 = -t130 * t86 - t143 * t72;
t16 = -t222 * t225 * t72 + (-t222 * t334 + t280) * t130;
t8 = t226 * t72 + t280 * t155 + (-t130 * t276 - t155 * t334 + t363) * t222;
t7 = t130 * t141 + t132 * t143 - t155 * t86 - t180 * t72;
t1 = qJD(6) * t155 + t3 + t368;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t324, -t228 * t354, 0, 0, 0, 0, 0, 0, 0, 0, t179 * t289 - t227 * t300 - t262 * t98, t179 * t288 + t223 * t300 - t262 * t99 (t223 * t98 + t227 * t99 + (t138 * t227 - t139 * t223) * qJD(3)) * t343, -t118 * t98 + t119 * t99 + t385 + t139 * t84 + (qJD(1) * t179 + t154) * t296, 0, 0, 0, 0, 0, 0, -t105 * t289 + t138 * t132 + t98 * t159 + t276 * t49, -t106 * t289 - t138 * t403 + t98 * t161 + t276 * t50, -t105 * t403 - t106 * t132 - t50 * t159 + t49 * t161, t103 * t98 - t105 * t27 - t106 * t251 - t49 * t62 + t50 * t63 + t385, 0, 0, 0, 0, 0, 0, t250, -t239, t252, t105 * t23 - t13 * t24 + t14 * t25 + t3 * t71 + t308 * t70 + t49 * t55, 0, 0, 0, 0, 0, 0, t250, t252, t239, t1 * t71 + t105 * t5 + t13 * t17 + t14 * t18 + t2 * t70 + t29 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t264, -0.2e1 * t346 * t315, t261 * t318, -0.2e1 * t264, -t261 * t319, 0 (t223 * t258 - t347) * qJD(3) + (-t347 * qJD(2) - t85) * t369 (t227 * t258 + t348) * qJD(3) + (t348 * qJD(2) - t84) * t369 (t223 * t85 + t227 * t84 + (-t118 * t227 - t119 * t223) * qJD(3) + ((-t348 - t175) * t227 + (t347 - t176) * t223) * qJD(2)) * t218, -t183 * t85 + t184 * t84 - t348 * t119 - t347 * t118 + (-t328 - t362) * t322, -t161 * t140 - t181 * t403, -t181 * t132 + t140 * t159 - t161 * t141 + t180 * t403, t140 * t276 + t161 * t319 + t181 * t289 + t355 * t403, t141 * t159 + t365, t141 * t276 + (t132 * t227 + (-qJD(2) * t180 - t159) * t339) * t218 (-t215 * t340 - t303) * t339, t169 * t132 + t85 * t180 + t103 * t141 + t347 * t159 + (-t27 * t227 + (qJD(2) * t120 + t62) * t339) * t218 - t411 * t276, -t103 * t140 - t121 * t289 + t347 * t161 - t169 * t403 + t85 * t181 - t251 * t355 + t412 * t276 - t319 * t63, t120 * t403 - t121 * t132 + t62 * t140 - t63 * t141 - t412 * t159 - t411 * t161 + t180 * t251 - t27 * t181, t347 * t103 + t120 * t27 - t121 * t251 + t169 * t85 + t411 * t62 + t412 * t63, t19, t238, t7, t275, -t236, t61, t108 * t73 + t128 * t372 + t132 * t53 + t141 * t24 + t142 * t23 + t155 * t388 - t180 * t308 + t55 * t87, -t108 * t72 + t130 * t372 - t132 * t409 - t141 * t25 + t143 * t23 - t155 * t389 - t180 * t3 - t55 * t86, -t128 * t389 - t130 * t388 - t142 * t3 + t143 * t308 + t24 * t86 - t25 * t87 - t409 * t73 + t53 * t72, t108 * t23 + t24 * t388 + t25 * t389 + t3 * t409 - t308 * t53 + t372 * t55, t19, t7, -t238, t61, t236, t275, -t128 * t371 - t132 * t48 - t141 * t17 + t142 * t5 + t155 * t396 - t180 * t2 + t29 * t87 + t57 * t73, -t1 * t142 - t128 * t397 - t130 * t396 + t143 * t2 - t17 * t86 - t18 * t87 - t47 * t73 - t48 * t72, t1 * t180 + t130 * t371 + t132 * t47 + t141 * t18 - t143 * t5 + t155 * t397 + t29 * t86 + t57 * t72, t1 * t47 - t17 * t396 + t18 * t397 + t2 * t48 - t29 * t371 + t5 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t299, t346 * t358, -t227 * t286, t299, t223 * t286, 0, t119 * t262 + t223 * t231 - t227 * t237, t118 * t262 + t223 * t237 + t227 * t231, 0, 0, -t226 * t292 + t420, -t226 * t403 + t407 * t159 + (-t132 + t292 - t330) * t222, -t226 * t302 + (t226 * t306 + (t222 * qJD(3) - t161) * t223) * t343, -t159 * t254 - t351, t222 * t302 + (-t222 * t306 + (qJD(3) * t226 + t159) * t223) * t343, t303 * t341, -pkin(3) * t132 + t222 * t331 + t91 * t276 - t119 * t159 + (pkin(10) * t302 - t85) * t226 + (-t223 * t62 + (-pkin(10) * t339 - t103 * t227) * t222) * t343, pkin(3) * t403 - t103 * t293 - t119 * t161 + t85 * t222 + t226 * t331 - t415 * t276 - t289 * t394 + t321 * t63, t91 * t161 + t270 + t422 * t63 + t407 * t62 + t415 * t159 + (-t351 + t420) * pkin(10), -pkin(3) * t85 - t103 * t119 - t62 * t91 - t63 * t92 + ((-t222 * t63 - t226 * t62) * qJD(4) + t270) * pkin(10), t16, t419, t8, t234, t418, t66, -t128 * t81 + t132 * t166 - t152 * t55 + t373 * t155 + (t308 + (pkin(10) * t128 + t221 * t55) * qJD(4)) * t226 + (pkin(10) * t73 - t24 * t276 + t333 * t55 + t382) * t222, -t130 * t81 - t132 * t167 - t153 * t55 - t374 * t155 + (t3 + (pkin(10) * t130 + t225 * t55) * qJD(4)) * t226 + (-pkin(10) * t72 + t25 * t276 - t334 * t55 + t380) * t222, t152 * t25 + t153 * t24 + t166 * t72 - t167 * t73 - t373 * t130 - t374 * t128 + t273 * t335 + (-t221 * t3 + t225 * t308 + (t221 * t24 - t225 * t25) * qJD(5)) * t222, -t166 * t308 + t167 * t3 - t55 * t81 + t374 * t25 + t373 * t24 + (t222 * t23 + t335 * t55) * pkin(10), t16, t8, -t419, t66, -t418, t234, -t132 * t146 - t152 * t29 + t172 * t73 + (qJD(4) * t221 * t29 + t2) * t226 - t376 * t155 + t375 * t128 + (t17 * t276 + t29 * t333 + t391) * t222, -t145 * t73 - t146 * t72 + t152 * t18 - t153 * t17 + t376 * t130 - t377 * t128 + t274 * t335 + (-t1 * t221 + t2 * t225 + (-t17 * t221 - t18 * t225) * qJD(5)) * t222, t132 * t145 + t153 * t29 + t172 * t72 + (-t29 * t336 - t1) * t226 + t377 * t155 - t375 * t130 + (-t18 * t276 + t29 * t334 - t390) * t222, t1 * t145 + t146 * t2 + t17 * t376 + t172 * t5 + t18 * t377 + t29 * t375; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t360, -t159 ^ 2 + t161 ^ 2, -t159 * t276 - t403, -t360, -t263 - t292, t289, -t103 * t161 - t276 * t63 + t27, t103 * t159 - t276 * t62 + t251, 0, 0, t28, -t417, t33, t243, -t402, -t361, -pkin(4) * t73 - t128 * t63 - t161 * t24 - t380 + (-t37 - t325) * t155 + t253 * t221, pkin(4) * t72 - t130 * t63 + t161 * t25 + t382 + (t38 + t326) * t155 + t253 * t225, t128 * t38 + t130 * t37 - t69 + (-t159 * t24 + t3 + (-t24 + t393) * qJD(5)) * t225 + (t279 - t259) * t221, -pkin(4) * t23 - t24 * t37 - t25 * t38 - t55 * t63 + (qJD(5) * t273 + t221 * t308 + t225 * t3) * pkin(11), t28, t33, t417, -t361, t402, t243, t161 * t17 + t199 * t73 - t390 + (t32 - t325) * t155 + t370 * t128 + t406 * t221, t128 * t31 - t130 * t32 - t69 + (t159 * t17 + t1 + (t17 + t393) * qJD(5)) * t225 + (t279 + t416) * t221, -t161 * t18 + t199 * t72 - t391 + (-t31 - t326) * t155 - t370 * t130 - t406 * t225, -t17 * t32 - t18 * t31 + t199 * t5 + t370 * t29 + (qJD(5) * t274 + t1 * t225 + t2 * t221) * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t367, -t316, t44, -t367, -t230, t132, -t130 * t55 + t259, t128 * t55 + t155 * t24 - t3, 0, 0, t367, t44, t316, t132, t230, -t367, -t128 * t78 + t259 - t386 + 0.2e1 * t395, pkin(5) * t72 - qJ(6) * t73 + (t18 - t25) * t130 + (t17 - t349) * t128, 0.2e1 * t368 - t128 * t29 + t130 * t78 + (0.2e1 * qJD(6) - t24) * t155 + t3, -pkin(5) * t2 + qJ(6) * t1 - t17 * t25 + t18 * t349 - t29 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t367 - t132, t44, -t398 - t399, t386 + t416;];
tauc_reg  = t4;
