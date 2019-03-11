% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6PRRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRPRP5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP5_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:49:23
% EndTime: 2019-03-08 21:49:33
% DurationCPUTime: 5.72s
% Computational Cost: add. (4856->558), mult. (10775->681), div. (0->0), fcn. (7793->10), ass. (0->296)
t205 = sin(qJ(3));
t330 = qJD(2) * qJD(3);
t304 = t205 * t330;
t208 = cos(qJ(3));
t324 = t208 * qJDD(2);
t421 = -t304 + t324;
t204 = sin(qJ(5));
t207 = cos(qJ(5));
t343 = qJD(2) * t208;
t309 = t204 * t343;
t337 = qJD(3) * t207;
t144 = -t309 + t337;
t345 = qJD(2) * t205;
t185 = qJD(5) + t345;
t259 = t144 * t185;
t339 = qJD(3) * t204;
t142 = t207 * t343 + t339;
t260 = t142 * t185;
t56 = qJD(5) * t142 - t207 * qJDD(3) + t421 * t204;
t57 = -qJD(5) * t309 + t204 * qJDD(3) + (qJD(3) * qJD(5) + t421) * t207;
t427 = (t57 + t259) * t207 - (t56 + t260) * t204;
t395 = pkin(4) + pkin(8);
t334 = qJD(5) * t207;
t55 = t57 * t204;
t426 = -t142 * t334 - t55;
t206 = sin(qJ(2));
t202 = sin(pkin(6));
t348 = qJD(1) * t202;
t316 = t206 * t348;
t152 = qJD(2) * pkin(8) + t316;
t203 = cos(pkin(6));
t347 = qJD(1) * t203;
t178 = t208 * t347;
t93 = t205 * t152 - t178;
t423 = -qJD(4) - t93;
t209 = cos(qJ(2));
t331 = qJD(1) * qJD(2);
t306 = t209 * t331;
t340 = qJD(3) * t203;
t422 = qJDD(2) * pkin(8) + (qJDD(1) * t206 + t306) * t202 + qJD(1) * t340;
t359 = t202 * t206;
t133 = t203 * t205 + t208 * t359;
t212 = qJD(2) ^ 2;
t256 = qJDD(2) * t209 - t206 * t212;
t302 = t209 * t330;
t358 = t202 * t209;
t314 = qJD(2) * t358;
t322 = t205 * t359;
t76 = -qJD(3) * t322 + (t314 + t340) * t208;
t420 = t76 * qJD(3) + t133 * qJDD(3) + (t205 * t256 + t208 * t302) * t202;
t132 = -t203 * t208 + t322;
t77 = qJD(3) * t133 + t205 * t314;
t419 = -t77 * qJD(3) - t132 * qJDD(3) + (-t205 * t302 + t208 * t256) * t202;
t378 = t207 * t56;
t16 = -t204 * t259 - t378;
t338 = qJD(3) * t205;
t418 = ((t142 * t207 + t144 * t204) * qJD(5) + t55 + t378) * t208 - (t142 * t204 - t144 * t207) * t338;
t354 = t205 * t209;
t101 = (t204 * t354 + t206 * t207) * t202;
t210 = -pkin(3) - pkin(9);
t356 = t205 * qJ(4);
t301 = -pkin(2) - t356;
t135 = t208 * t210 + t301;
t336 = qJD(3) * t208;
t149 = t395 * t336;
t168 = t395 * t205;
t335 = qJD(5) * t204;
t191 = pkin(3) * t338;
t368 = qJ(4) * t208;
t273 = pkin(9) * t205 - t368;
t332 = t205 * qJD(4);
t230 = qJD(3) * t273 - t332;
t97 = t191 + t230;
t384 = -qJD(1) * t101 - t135 * t335 + t204 * t149 + t168 * t334 + t207 * t97;
t303 = t208 * t330;
t325 = t205 * qJDD(2);
t240 = t303 + t325;
t140 = qJDD(5) + t240;
t389 = pkin(5) * t140;
t329 = qJDD(1) * t203;
t288 = t152 * t336 + t422 * t205 - t208 * t329;
t266 = qJDD(4) + t288;
t15 = pkin(4) * t240 + qJDD(3) * t210 + t266;
t307 = t206 * t331;
t272 = -qJDD(1) * t358 + t202 * t307;
t253 = pkin(3) * t304 + t272;
t27 = qJD(2) * t230 + qJDD(2) * t135 + t253;
t408 = qJD(4) - t178;
t353 = (pkin(4) * qJD(2) + t152) * t205 + t408;
t58 = qJD(3) * t210 + t353;
t310 = t209 * t348;
t83 = qJD(2) * t135 - t310;
t4 = t207 * t15 - t204 * t27 - t83 * t334 - t58 * t335;
t2 = qJDD(6) - t4 - t389;
t23 = t204 * t58 + t207 * t83;
t13 = qJ(6) * t185 + t23;
t381 = t13 * t185;
t417 = -t2 + t381;
t379 = t185 * t23;
t416 = t4 + t379;
t94 = t152 * t208 + t205 * t347;
t86 = -qJD(3) * qJ(4) - t94;
t85 = -qJD(3) * pkin(3) - t423;
t414 = t57 - t259;
t412 = t207 * t135 + t204 * t168;
t124 = t204 * t140;
t411 = -t185 * t334 - t124;
t125 = t207 * t140;
t342 = qJD(3) * t142;
t396 = t185 ^ 2;
t410 = -t204 * t396 + t125 - t342;
t409 = -t185 * t335 + t125;
t155 = -pkin(3) * t208 + t301;
t346 = qJD(2) * t155;
t96 = -t310 + t346;
t407 = t96 * t345 + qJDD(4);
t289 = t152 * t338 - t205 * t329 - t422 * t208;
t406 = (-t205 * t94 + t208 * t93) * qJD(3) + t288 * t205 - t289 * t208;
t197 = qJDD(3) * qJ(4);
t198 = qJD(3) * qJD(4);
t28 = -t197 - t198 + t289;
t365 = qJDD(3) * pkin(3);
t29 = t266 - t365;
t405 = (t205 * t86 + t208 * t85) * qJD(3) + t29 * t205 - t28 * t208;
t370 = sin(pkin(10));
t293 = t370 * t209;
t371 = cos(pkin(10));
t296 = t371 * t206;
t127 = t203 * t296 + t293;
t298 = t202 * t371;
t73 = t127 * t208 - t205 * t298;
t294 = t370 * t206;
t295 = t371 * t209;
t129 = -t203 * t294 + t295;
t297 = t202 * t370;
t75 = t129 * t208 + t205 * t297;
t242 = g(1) * t75 + g(2) * t73 + g(3) * t133;
t362 = t140 * t210;
t190 = pkin(4) * t343;
t65 = t190 - t86;
t404 = t185 * t65 + t362;
t366 = qJDD(2) * pkin(2);
t103 = t272 - t366;
t126 = -t203 * t295 + t294;
t128 = t203 * t293 + t296;
t279 = g(1) * t128 + g(2) * t126;
t211 = qJD(3) ^ 2;
t388 = pkin(8) * t211;
t403 = t202 * (-g(3) * t209 + t307) - t103 + t279 + t366 - t388;
t248 = -qJ(4) * t336 - t332;
t123 = t191 + t248;
t236 = -g(3) * t358 + t279;
t328 = qJDD(2) * t155;
t40 = qJD(2) * t248 + t253 + t328;
t400 = qJD(2) * (-t123 + t316) + t236 - t328 - t388 - t40;
t321 = t207 * t358;
t287 = t205 * t321;
t383 = -qJD(1) * t287 - qJD(5) * t412 + t207 * t149 + (t316 - t97) * t204;
t398 = t205 * (t185 * t337 - t57) + t208 * (-t342 - t409);
t397 = t144 ^ 2;
t72 = t127 * t205 + t208 * t298;
t392 = t72 * pkin(9);
t74 = t129 * t205 - t208 * t297;
t391 = t74 * pkin(9);
t390 = qJ(6) * t336 + qJD(6) * t205 + t384;
t387 = pkin(9) * t132;
t385 = -pkin(5) * t336 - t383;
t382 = qJD(2) * pkin(2);
t22 = -t204 * t83 + t207 * t58;
t380 = t185 * t22;
t377 = t210 * t56;
t192 = pkin(3) * t345;
t115 = qJD(2) * t273 + t192;
t82 = t190 + t94;
t38 = t207 * t115 + t204 * t82;
t275 = pkin(5) * t207 + qJ(6) * t204;
t255 = -pkin(4) - t275;
t372 = qJD(5) * t275 - t207 * qJD(6) - (qJD(2) * t255 - t152) * t205 + t408;
t369 = pkin(8) * qJDD(3);
t367 = qJ(6) * t140;
t364 = t126 * t208;
t363 = t128 * t208;
t361 = t144 * t142;
t360 = t144 * t210;
t357 = t204 * t205;
t355 = t205 * t207;
t352 = qJD(6) - t22;
t351 = qJDD(1) - g(3);
t350 = pkin(2) * t358 + pkin(8) * t359;
t169 = t395 * t208;
t200 = t205 ^ 2;
t201 = t208 ^ 2;
t349 = t200 - t201;
t344 = qJD(2) * t206;
t341 = qJD(3) * t144;
t333 = qJD(5) * t208;
t327 = qJDD(2) * t200;
t326 = qJDD(2) * t201;
t313 = t207 * t345;
t323 = t142 * t313 - t426;
t320 = t208 * t358;
t117 = t126 * pkin(2);
t319 = -pkin(3) * t364 - t126 * t356 - t117;
t119 = t128 * pkin(2);
t318 = -pkin(3) * t363 - t128 * t356 - t119;
t317 = t185 * t313 - t411;
t315 = t202 * t344;
t311 = t185 * t343;
t308 = t142 ^ 2 - t397;
t67 = t72 * pkin(3);
t300 = qJ(4) * t73 - t67;
t68 = t74 * pkin(3);
t299 = qJ(4) * t75 - t68;
t122 = t132 * pkin(3);
t290 = qJ(4) * t133 - t122;
t286 = t202 * qJ(4) * t354 + pkin(3) * t320 + t350;
t285 = t142 * t310;
t284 = t144 * t310;
t283 = t208 * t310;
t282 = t205 * t303;
t278 = g(1) * t129 + g(2) * t127;
t274 = -pkin(5) * t204 + qJ(6) * t207;
t12 = -pkin(5) * t185 + t352;
t271 = t12 * t204 + t13 * t207;
t270 = -t204 * t22 + t207 * t23;
t37 = -t115 * t204 + t207 * t82;
t69 = -t135 * t204 + t168 * t207;
t154 = qJ(4) - t274;
t251 = pkin(4) * t359 + pkin(9) * t320 + t286;
t78 = t132 * t207 + t204 * t358;
t3 = t204 * t15 + t207 * t27 + t58 * t334 - t335 * t83;
t20 = -qJD(5) * t321 - t77 * t207 + (qJD(5) * t132 + t315) * t204;
t21 = qJD(5) * t78 + t77 * t204 + t207 * t315;
t79 = -t132 * t204 + t321;
t249 = -t21 * t142 + t144 * t20 + t56 * t78 + t79 * t57;
t100 = t204 * t359 - t287;
t47 = t126 * t355 + t127 * t204;
t49 = t128 * t355 + t129 * t204;
t245 = -g(1) * t49 - g(2) * t47 - g(3) * t100;
t48 = -t126 * t357 + t127 * t207;
t50 = -t128 * t357 + t129 * t207;
t244 = -g(1) * t50 - g(2) * t48 - g(3) * t101;
t243 = g(1) * t74 + g(2) * t72 + g(3) * t132;
t241 = t133 * t57 + t140 * t78 + t76 * t142 - t185 * t20;
t239 = t142 * t343 - t317;
t238 = -pkin(9) * t364 + t395 * t127 + t319;
t237 = -pkin(9) * t363 + t395 * t129 + t318;
t30 = pkin(5) * t142 - qJ(6) * t144 + t65;
t235 = t185 * t30 + t362;
t234 = t426 * t210 + t243;
t231 = t133 * t56 - t140 * t79 - t144 * t76 + t185 * t21;
t41 = t126 * t204 - t72 * t207;
t43 = t128 * t204 - t74 * t207;
t229 = g(1) * t43 + g(2) * t41 - g(3) * t78 + t4;
t228 = -qJD(5) * t185 * t210 - t242;
t226 = t243 - t288;
t17 = t421 * pkin(4) - t28;
t5 = pkin(5) * t57 + qJ(6) * t56 - qJD(6) * t144 + t17;
t225 = -t228 - t5;
t224 = t17 + t228;
t153 = -t310 - t382;
t223 = -t369 + (t153 + t310 - t382) * qJD(3);
t222 = t369 + (-t310 - t96 - t346) * qJD(3);
t221 = -t16 - t323;
t42 = t126 * t207 + t204 * t72;
t44 = t128 * t207 + t204 * t74;
t220 = -g(1) * t44 - g(2) * t42 + g(3) * t79 + t3;
t219 = t57 * t207 * t208 + (-t204 * t333 - t205 * t337) * t142;
t218 = qJD(3) * t93 - t242 - t289;
t217 = qJD(3) * t94 + t226;
t216 = t144 * t30 + qJDD(6) - t229;
t214 = (-g(3) * t206 + (-t200 - t201) * t306) * t202 - t278 + (t327 + t326) * pkin(8);
t213 = (t132 * t205 + t133 * t208) * qJDD(2) + (t205 * t77 + t208 * t76 + (t132 * t208 - t133 * t205) * qJD(3)) * qJD(2);
t180 = t205 * t212 * t208;
t160 = t349 * t212;
t157 = qJDD(3) * t208 - t205 * t211;
t156 = qJDD(3) * t205 + t208 * t211;
t148 = t395 * t338;
t147 = -qJ(4) * t343 + t192;
t137 = -0.2e1 * t282 + t326;
t136 = 0.2e1 * t282 + t327;
t95 = 0.2e1 * t205 * t324 - 0.2e1 * t330 * t349;
t92 = t208 * t275 + t169;
t80 = t140 * t205 + t185 * t336;
t71 = pkin(5) * t144 + qJ(6) * t142;
t60 = -pkin(5) * t205 - t69;
t59 = qJ(6) * t205 + t412;
t45 = (qJD(5) * t274 + qJD(6) * t204) * t208 + (-pkin(8) + t255) * t338;
t36 = t260 - t56;
t33 = -pkin(5) * t343 - t37;
t32 = qJ(6) * t343 + t38;
t31 = (-t144 * t208 - t185 * t357) * qJD(2) + t409;
t9 = t56 * t204 * t208 + (t204 * t338 - t207 * t333) * t144;
t6 = (t185 * t339 - t56) * t205 + (t341 + t411) * t208;
t1 = qJD(6) * t185 + t3 + t367;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t351, 0, 0, 0, 0, 0, 0, t256 * t202 (-qJDD(2) * t206 - t209 * t212) * t202, 0, -g(3) + (t203 ^ 2 + (t206 ^ 2 + t209 ^ 2) * t202 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t419, -t420, t213, t288 * t132 - t289 * t133 + t94 * t76 + t93 * t77 - g(3) + (-t103 * t209 + t153 * t344) * t202, 0, 0, 0, 0, 0, 0, t213, -t419, t420, t29 * t132 - t28 * t133 - t86 * t76 + t85 * t77 - g(3) + (-t209 * t40 + t344 * t96) * t202, 0, 0, 0, 0, 0, 0, t241, -t231, t249, t133 * t17 - t20 * t22 + t21 * t23 - t3 * t79 + t4 * t78 + t65 * t76 - g(3), 0, 0, 0, 0, 0, 0, t241, t249, t231, -t1 * t79 + t12 * t20 + t13 * t21 + t133 * t5 - t2 * t78 + t30 * t76 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t351 * t358 + t279, -t351 * t359 + t278, 0, 0, t136, t95, t156, t137, t157, 0, t223 * t205 + t403 * t208, -t403 * t205 + t223 * t208, t214 + t406, -t103 * pkin(2) + g(1) * t119 + g(2) * t117 - g(3) * t350 + (-t153 * t206 + (-t205 * t93 - t208 * t94) * t209) * t348 + (-t278 + t406) * pkin(8), 0, -t156, -t157, t136, t95, t137, t214 + t405, t222 * t205 - t400 * t208, t400 * t205 + t222 * t208, t40 * t155 + t96 * t123 - g(1) * t318 - g(2) * t319 - g(3) * t286 + (-t206 * t96 + (-t205 * t85 + t208 * t86) * t209) * t348 + (-t278 + t405) * pkin(8), t9, t418, t6, t219, t398, t80, t140 * t69 - t142 * t148 + t169 * t57 + (-t337 * t65 + t4) * t205 + t383 * t185 + (qJD(3) * t22 + t17 * t207 - t335 * t65 - t285) * t208 + t244, -t140 * t412 - t144 * t148 - t169 * t56 + (t339 * t65 - t3) * t205 - t384 * t185 + (-qJD(3) * t23 - t17 * t204 - t334 * t65 - t284) * t208 - t245, t56 * t69 - t57 * t412 - t383 * t144 - t384 * t142 + t270 * t338 + (t204 * t4 - t207 * t3 + (t204 * t23 + t207 * t22) * qJD(5) + t236) * t208, t3 * t412 + t4 * t69 + t17 * t169 - g(1) * t237 - g(2) * t238 - g(3) * t251 + (-t148 - t283) * t65 + t384 * t23 + t383 * t22, t9, t6, -t418, t80, -t398, t219, -t140 * t60 + t142 * t45 + t57 * t92 + (-t30 * t337 - t2) * t205 - t385 * t185 + (-qJD(3) * t12 + t207 * t5 - t30 * t335 - t285) * t208 + t244, -t56 * t60 - t57 * t59 + t385 * t144 - t390 * t142 + t271 * t338 + (-t1 * t207 - t2 * t204 + (-t12 * t207 + t13 * t204) * qJD(5) + t236) * t208, t140 * t59 - t144 * t45 + t56 * t92 + (-t30 * t339 + t1) * t205 + t390 * t185 + (qJD(3) * t13 + t204 * t5 + t30 * t334 + t284) * t208 + t245, t1 * t59 + t5 * t92 + t2 * t60 - g(1) * (pkin(5) * t50 + qJ(6) * t49 + t237) - g(2) * (pkin(5) * t48 + qJ(6) * t47 + t238) - g(3) * (pkin(5) * t101 + qJ(6) * t100 + t251) + (t45 - t283) * t30 + t390 * t13 + t385 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t180, t160, t325, t180, t324, qJDD(3), -t153 * t345 + t217, -t153 * t343 - t218, 0, 0, qJDD(3), -t325, -t324, -t180, t160, t180 (-pkin(3) * t205 + t368) * qJDD(2), -t147 * t343 - t217 - 0.2e1 * t365 + t407, 0.2e1 * t197 + 0.2e1 * t198 + (t147 * t205 + t208 * t96) * qJD(2) + t218, -t29 * pkin(3) - g(1) * t299 - g(2) * t300 - g(3) * t290 - t28 * qJ(4) - t96 * t147 + t423 * t86 - t85 * t94, t16, -t427, t31, t323, t239, -t311, qJ(4) * t57 + t353 * t142 - t185 * t37 + t224 * t204 + t404 * t207 - t22 * t343, -qJ(4) * t56 + t353 * t144 + t185 * t38 - t404 * t204 + t224 * t207 + t23 * t343, t142 * t38 + t144 * t37 + (t377 - t416) * t207 + (t22 * t345 - t3 + (t22 + t360) * qJD(5)) * t204 + t234, t17 * qJ(4) - t23 * t38 - t22 * t37 - g(1) * (t299 - t391) - g(2) * (t300 - t392) - g(3) * (t290 - t387) + t353 * t65 + (qJD(5) * t270 + t3 * t204 + t4 * t207) * t210, t16, t31, t427, -t311, -t239, t207 * t260 + t55, t12 * t343 + t142 * t372 + t154 * t57 + t185 * t33 - t204 * t225 + t207 * t235, t142 * t32 - t144 * t33 + (t377 - t417) * t207 + (-t12 * t345 - t1 + (-t12 + t360) * qJD(5)) * t204 + t234, -t13 * t343 - t144 * t372 + t154 * t56 - t185 * t32 + t204 * t235 + t207 * t225, -t13 * t32 - t12 * t33 - g(1) * (-t68 - t391) - g(2) * (-t67 - t392) - g(3) * (-t122 - t387) + t372 * t30 + (qJD(5) * t271 + t1 * t204 - t2 * t207) * t210 + (t5 - t242) * t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t325, qJDD(3) + t180, -t200 * t212 - t211, qJD(3) * t86 - t226 - t365 + t407, 0, 0, 0, 0, 0, 0, t410, -t207 * t396 - t124 - t341, t221, -qJD(3) * t65 + t416 * t207 + (t3 - t380) * t204 - t243, 0, 0, 0, 0, 0, 0, t410, t221, t317 + t341, -qJD(3) * t30 + t417 * t207 + (t12 * t185 + t1) * t204 - t243; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t361, -t308, t36, -t361, -t414, t140, -t144 * t65 + t229 + t379, t142 * t65 - t220 + t380, 0, 0, t361, t36, t308, t140, t414, -t361, -t142 * t71 - t216 + t379 + 0.2e1 * t389, pkin(5) * t56 - t57 * qJ(6) + (t13 - t23) * t144 + (t12 - t352) * t142, 0.2e1 * t367 - t142 * t30 + t144 * t71 + (0.2e1 * qJD(6) - t22) * t185 + t220, t1 * qJ(6) - t2 * pkin(5) - t30 * t71 - t12 * t23 - g(1) * (-pkin(5) * t43 + qJ(6) * t44) - g(2) * (-pkin(5) * t41 + qJ(6) * t42) - g(3) * (pkin(5) * t78 - qJ(6) * t79) + t352 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140 + t361, t36, -t396 - t397, t216 - t381 - t389;];
tau_reg  = t7;
