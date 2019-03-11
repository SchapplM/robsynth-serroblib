% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RRRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 04:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRRRR7_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR7_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_inertiaDJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:41:07
% EndTime: 2019-03-10 04:41:54
% DurationCPUTime: 21.52s
% Computational Cost: add. (27638->769), mult. (70847->1333), div. (0->0), fcn. (70418->12), ass. (0->300)
t404 = sin(qJ(5));
t323 = t404 * qJD(5);
t324 = t404 * qJD(4);
t425 = t324 + t323;
t406 = cos(qJ(5));
t325 = t406 * qJD(5);
t424 = qJD(4) * t406 + t325;
t177 = sin(qJ(4));
t180 = cos(qJ(4));
t175 = sin(pkin(6));
t178 = sin(qJ(3));
t394 = cos(pkin(6));
t322 = t394 * t178;
t179 = sin(qJ(2));
t181 = cos(qJ(3));
t388 = t179 * t181;
t264 = t175 * t388 + t322;
t182 = cos(qJ(2));
t391 = t175 * t182;
t246 = -t177 * t391 + t180 * t264;
t149 = t180 * t391;
t384 = t177 * t264 + t149;
t274 = t246 * t404 + t384 * t406;
t392 = t175 * t179;
t123 = t178 * t392 - t181 * t394;
t350 = pkin(1) * t394;
t306 = t179 * t350;
t255 = pkin(9) * t394 + t306;
t247 = t181 * t255;
t353 = -pkin(9) * t179 - pkin(1);
t206 = t247 + (t178 * t353 + (-pkin(2) * t178 + pkin(8) * t181 - pkin(10)) * t182) * t175;
t305 = t182 * t350;
t256 = -pkin(2) * t394 - t305;
t365 = pkin(8) * t392;
t116 = t256 + t365;
t399 = t123 * pkin(3);
t208 = -pkin(10) * t264 + t116 + t399;
t192 = t123 * pkin(4) - pkin(11) * t246 - t177 * t206 + t180 * t208;
t218 = -pkin(10) * t322 + t256 + t399;
t279 = -pkin(2) * t182 + t353;
t385 = t181 * t182;
t244 = pkin(8) * t385 + t178 * t279;
t230 = -t182 * pkin(10) + t244;
t277 = pkin(8) * t179 - pkin(10) * t388;
t45 = t180 * t247 + t177 * t218 + (t177 * t277 + t180 * t230) * t175;
t194 = -pkin(11) * t384 + t45;
t19 = t192 * t406 - t194 * t404;
t402 = pkin(3) * t181;
t407 = -pkin(11) - pkin(10);
t271 = t178 * t407 - pkin(2) - t402;
t257 = t271 * t180;
t398 = t177 * pkin(9);
t351 = -pkin(4) - t398;
t223 = t181 * t351 + t257;
t387 = t180 * t181;
t363 = pkin(9) * t387;
t235 = t177 * t271 + t363;
t63 = t223 * t406 - t235 * t404;
t66 = t246 * t406 - t384 * t404;
t258 = -t305 + t365;
t248 = t258 * qJD(2);
t270 = t279 * t175;
t423 = -qJD(3) * t270 + t248;
t308 = t407 * t406;
t347 = t404 * t180;
t101 = t177 * t308 + t347 * t407;
t372 = qJD(4) * t180;
t375 = qJD(3) * t181;
t422 = t177 * t375 + t178 * t372;
t130 = t177 * t406 + t347;
t307 = t407 * t404;
t131 = t177 * t307;
t348 = t406 * t180;
t102 = -t348 * t407 + t131;
t421 = t102 * qJD(5);
t328 = qJD(3) * t404;
t297 = t181 * t328;
t329 = qJD(3) * t406;
t298 = t181 * t329;
t304 = t178 * t348;
t310 = t177 * t298 + t180 * t297 + (qJD(4) + qJD(5)) * t304;
t419 = t177 * t425;
t420 = -t178 * t419 + t310;
t309 = t177 * t424 + t180 * t425;
t20 = t192 * t404 + t194 * t406;
t64 = t223 * t404 + t235 * t406;
t252 = qJD(4) * t264;
t378 = qJD(2) * t179;
t333 = t175 * t378;
t345 = qJD(3) * t392;
t136 = t178 * t345;
t319 = t394 * qJD(3);
t377 = qJD(2) * t182;
t332 = t175 * t377;
t410 = -t181 * (t319 + t332) + t136;
t212 = -qJD(4) * t149 - t180 * t410 + (-t252 + t333) * t177;
t211 = t212 * t177;
t362 = pkin(8) * t391;
t238 = t255 + t362;
t229 = qJD(3) * t238;
t295 = pkin(2) * t179 - pkin(9) * t182;
t379 = qJD(2) * t175;
t268 = t295 * t379;
t52 = t423 * t181 + (t229 - t268) * t178;
t418 = -pkin(10) * t333 - qJD(4) * t208 + t52;
t249 = t264 * qJD(3);
t100 = t178 * t332 + t249;
t320 = qJD(2) * t394;
t296 = t179 * t320;
t417 = -pkin(1) * t296 - t100 * pkin(3) - pkin(8) * t332 - pkin(10) * t410 + qJD(4) * t206;
t232 = qJD(4) * t246;
t416 = -t177 * t232 + t180 * t212;
t415 = t180 * t232 + t211;
t316 = t384 * qJD(4);
t373 = qJD(4) * t177;
t336 = t182 * t373;
t395 = t177 * t410 - t180 * t252;
t409 = -t175 * (t180 * t378 + t336) - t395;
t414 = t177 * t316 - t180 * t409;
t338 = t178 * t373;
t341 = t180 * t375;
t412 = t338 - t341;
t411 = -t100 * t180 + t123 * t373;
t168 = qJD(3) * t178;
t73 = -t181 * t100 + t123 * t168;
t171 = t177 ^ 2;
t173 = t180 ^ 2;
t381 = t171 - t173;
t318 = qJD(4) * t381;
t157 = t177 * t324;
t408 = t177 * t323 - t180 * t424 + t157;
t174 = t181 ^ 2;
t405 = cos(qJ(6));
t403 = pkin(3) * t178;
t401 = pkin(9) * t175;
t400 = t100 * pkin(5);
t397 = t178 * pkin(4);
t169 = t178 * pkin(9);
t396 = t180 * pkin(3);
t390 = t177 * t178;
t389 = t178 * t180;
t383 = t130 * t178;
t133 = pkin(4) * t390 + t169;
t172 = t178 ^ 2;
t380 = t172 - t174;
t376 = qJD(3) * t180;
t374 = qJD(3) * t182;
t371 = qJD(4) * t181;
t176 = sin(qJ(6));
t370 = qJD(6) * t176;
t369 = t116 * qJD(3);
t368 = -0.2e1 * pkin(2) * qJD(3);
t367 = -0.2e1 * pkin(3) * qJD(4);
t78 = 0.2e1 * t123 * t100;
t366 = t181 * t398;
t364 = pkin(9) * t389;
t361 = t406 * pkin(4);
t360 = t404 * pkin(4);
t359 = pkin(5) * t168;
t167 = pkin(9) * t375;
t358 = pkin(4) * t373;
t357 = pkin(10) * t372;
t356 = pkin(5) * t370;
t355 = t407 * t181;
t354 = -t178 * t423 + t181 * t229;
t107 = pkin(4) * t422 + t167;
t166 = -pkin(4) * t180 - pkin(3);
t165 = t404 * t177;
t170 = t175 ^ 2;
t346 = t170 * t377;
t342 = t178 * t374;
t340 = t123 * t375;
t337 = t177 * t371;
t334 = t180 * t371;
t331 = t177 * t372;
t330 = t178 * t375;
t327 = qJD(6) * t405;
t321 = -0.2e1 * t356;
t317 = t384 * qJD(3);
t315 = t380 * qJD(3);
t314 = 0.2e1 * t330;
t313 = pkin(4) * t325;
t312 = pkin(4) * t323;
t311 = pkin(5) * t327;
t303 = t180 * t330;
t302 = t172 * t331;
t301 = t179 * t346;
t300 = t361 + pkin(5);
t299 = t405 * t404;
t294 = -pkin(10) * t178 - t402;
t293 = -pkin(10) * t181 + t403;
t291 = t405 * t383;
t289 = t181 * t317;
t287 = t165 - t348;
t21 = t177 * t417 + t180 * t418;
t22 = t177 * t418 - t180 * t417;
t286 = -t22 * t177 - t21 * t180;
t53 = t181 * t268 - t354;
t285 = -t53 * t178 - t52 * t181;
t278 = pkin(2) - t294;
t112 = -t177 * t278 + t363;
t245 = t180 * t278 + t366;
t282 = -t112 * t177 + t180 * t245;
t281 = qJD(4) * t308;
t280 = qJD(4) * t307;
t17 = -pkin(12) * t274 + t20;
t196 = qJD(5) * t274 - t212 * t406 + t404 * t409;
t184 = t100 * pkin(4) - pkin(11) * t212 + t22;
t185 = -pkin(11) * t409 - t21;
t7 = -qJD(5) * t20 + t406 * t184 - t404 * t185;
t183 = pkin(12) * t196 + t400 + t7;
t188 = t123 * pkin(5) - t66 * pkin(12) + t19;
t186 = t405 * t188;
t214 = qJD(5) * t66 + t404 * t212 + t406 * t409;
t6 = -qJD(5) * t19 - t404 * t184 - t406 * t185;
t190 = -pkin(12) * t214 - t6;
t1 = -qJD(6) * t186 + t17 * t370 - t176 * t183 - t190 * t405;
t48 = (-t179 * pkin(3) - t181 * t295) * t379 + t354;
t76 = -t178 * t238 + t181 * t270;
t72 = pkin(3) * t391 - t76;
t276 = t48 * t177 + t372 * t72;
t275 = t177 * t293;
t243 = t177 * t297 + t178 * t309 - t180 * t298;
t202 = -t235 * qJD(4) + (t180 * t355 + (-t351 + t396) * t178) * qJD(3);
t203 = (t257 - t366) * qJD(4) + (-t364 + (t355 + t403) * t177) * qJD(3);
t32 = -qJD(5) * t64 + t406 * t202 - t404 * t203;
t189 = pkin(12) * t243 + t32 + t359;
t119 = -t165 * t178 + t304;
t200 = -t181 * pkin(5) - t119 * pkin(12) + t63;
t198 = t405 * t200;
t31 = -qJD(5) * t63 - t404 * t202 - t406 * t203;
t201 = -pkin(12) * t420 - t31;
t56 = -pkin(12) * t383 + t64;
t10 = -qJD(6) * t198 - t176 * t189 - t201 * t405 + t370 * t56;
t273 = t176 * t287;
t272 = t100 * t177 + t123 * t372;
t269 = t405 * t300;
t267 = t405 * t287;
t266 = t178 * t378 - t181 * t374;
t263 = t176 * t274;
t81 = t119 * t405 - t176 * t383;
t128 = t306 + t362;
t253 = t405 * t274;
t70 = -qJD(5) * t101 - t177 * t281 - t180 * t280;
t74 = t245 * qJD(4) + (-t275 + t364) * qJD(3);
t75 = -t112 * qJD(4) + (pkin(9) * t390 + t180 * t293) * qJD(3);
t234 = qJD(4) * t282 - t177 * t75 - t180 * t74;
t233 = qJD(3) * t246;
t85 = -qJD(6) * t269 - t405 * t313 + (qJD(6) * t404 + t323) * t176 * pkin(4);
t227 = t181 * t233;
t224 = -t130 * pkin(12) + t101;
t58 = pkin(4) * t384 + t72;
t219 = t176 * t224;
t217 = -pkin(12) * t309 - t70;
t216 = t405 * t224;
t213 = t177 * t409 + t180 * t316;
t209 = (qJD(5) + qJD(6)) * (-t176 * t406 - t299) * pkin(4);
t199 = t176 * t200;
t197 = pkin(12) * t408 - t177 * t280 + t180 * t281 - t421;
t11 = -qJD(6) * t199 - t176 * t201 + t189 * t405 - t327 * t56;
t187 = t176 * t188;
t2 = -qJD(6) * t187 - t17 * t327 - t176 * t190 + t183 * t405;
t156 = -0.2e1 * t330;
t134 = -0.2e1 * t301;
t122 = pkin(4) * t299 + t176 * t300;
t121 = -t176 * t360 + t269;
t118 = qJD(2) * t128;
t115 = pkin(5) * t287 + t166;
t103 = -t177 * t341 + t178 * t318;
t98 = pkin(5) * t383 + t133;
t95 = t130 * t405 - t273;
t94 = t130 * t176 + t267;
t86 = t209 - t356;
t83 = pkin(5) * t309 + t358;
t82 = -pkin(12) * t287 + t102;
t80 = t119 * t176 + t291;
t77 = t175 * t244 + t247;
t71 = -t421 + (t180 * t308 - t131) * qJD(4);
t57 = pkin(5) * t420 + t107;
t55 = t405 * t82 + t219;
t54 = -t176 * t82 + t216;
t47 = -qJD(6) * t273 + t130 * t327 - t176 * t408 + t309 * t405;
t46 = qJD(6) * t267 + t130 * t370 + t176 * t309 + t405 * t408;
t44 = -t177 * t247 + t180 * t218 + (-t177 * t230 + t180 * t277) * t175;
t43 = t405 * t66 - t263;
t42 = t176 * t66 + t253;
t40 = pkin(5) * t274 + t58;
t37 = qJD(6) * t81 - t176 * t243 + t405 * t420;
t36 = qJD(6) * t291 + t119 * t370 + t176 * t420 + t243 * t405;
t35 = -t395 * pkin(4) + (-pkin(4) * t336 + (pkin(9) * t385 + (-pkin(2) * t181 + t166) * t179) * qJD(2)) * t175 + t354;
t34 = t405 * t56 + t199;
t33 = -t176 * t56 + t198;
t26 = -qJD(6) * t219 - t176 * t217 + t197 * t405 - t327 * t82;
t25 = -qJD(6) * t216 - t176 * t197 - t217 * t405 + t370 * t82;
t16 = -pkin(3) * t333 + pkin(4) * t409 + pkin(5) * t214 - t53;
t13 = -qJD(6) * t263 - t176 * t196 + t214 * t405 + t327 * t66;
t12 = qJD(6) * t253 + t176 * t214 + t196 * t405 + t370 * t66;
t9 = t17 * t405 + t187;
t8 = -t17 * t176 + t186;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t301, 0.2e1 * (-t179 ^ 2 + t182 ^ 2) * t170 * qJD(2), 0.2e1 * t320 * t391, t134, -0.2e1 * t175 * t296, 0, -0.2e1 * pkin(1) * t170 * t378 - 0.2e1 * t118 * t394, -0.2e1 * pkin(1) * t346 + 0.2e1 * t248 * t394, 0.2e1 * t118 * t392 - 0.2e1 * t128 * t333 - 0.2e1 * t248 * t391 + 0.2e1 * t258 * t332, 0.2e1 * t118 * t258 - 0.2e1 * t128 * t248, -0.2e1 * t264 * t410, -0.2e1 * t100 * t264 + 0.2e1 * t123 * t410, 0.2e1 * t264 * t333 + 0.2e1 * t391 * t410, t78, 0.2e1 * (t100 * t182 - t123 * t378) * t175, t134, 0.2e1 * t116 * t100 + 0.2e1 * t118 * t123 + 0.2e1 * (-t182 * t53 + t378 * t76) * t175, -0.2e1 * t116 * t410 + 0.2e1 * t118 * t264 - 0.2e1 * t333 * t77 - 0.2e1 * t391 * t52, -0.2e1 * t77 * t100 + 0.2e1 * t52 * t123 - 0.2e1 * t264 * t53 + 0.2e1 * t410 * t76, 0.2e1 * t116 * t118 - 0.2e1 * t52 * t77 + 0.2e1 * t53 * t76, 0.2e1 * t246 * t212, -0.2e1 * t212 * t384 - 0.2e1 * t246 * t409, 0.2e1 * t100 * t246 + 0.2e1 * t123 * t212, 0.2e1 * t384 * t409, -0.2e1 * t100 * t384 - 0.2e1 * t123 * t409, t78, 0.2e1 * t44 * t100 + 0.2e1 * t22 * t123 + 0.2e1 * t384 * t48 + 0.2e1 * t409 * t72, -0.2e1 * t100 * t45 + 0.2e1 * t123 * t21 + 0.2e1 * t212 * t72 + 0.2e1 * t246 * t48, 0.2e1 * t21 * t384 - 0.2e1 * t212 * t44 - 0.2e1 * t22 * t246 - 0.2e1 * t409 * t45, -0.2e1 * t21 * t45 + 0.2e1 * t22 * t44 + 0.2e1 * t48 * t72, -0.2e1 * t66 * t196, 0.2e1 * t196 * t274 - 0.2e1 * t214 * t66, 0.2e1 * t100 * t66 - 0.2e1 * t123 * t196, 0.2e1 * t274 * t214, -0.2e1 * t100 * t274 - 0.2e1 * t123 * t214, t78, 0.2e1 * t100 * t19 + 0.2e1 * t123 * t7 + 0.2e1 * t214 * t58 + 0.2e1 * t274 * t35, -0.2e1 * t100 * t20 + 0.2e1 * t123 * t6 - 0.2e1 * t196 * t58 + 0.2e1 * t35 * t66, 0.2e1 * t19 * t196 - 0.2e1 * t20 * t214 + 0.2e1 * t274 * t6 - 0.2e1 * t66 * t7, 0.2e1 * t19 * t7 - 0.2e1 * t20 * t6 + 0.2e1 * t35 * t58, -0.2e1 * t43 * t12, 0.2e1 * t12 * t42 - 0.2e1 * t13 * t43, 0.2e1 * t100 * t43 - 0.2e1 * t12 * t123, 0.2e1 * t42 * t13, -0.2e1 * t100 * t42 - 0.2e1 * t123 * t13, t78, 0.2e1 * t100 * t8 + 0.2e1 * t123 * t2 + 0.2e1 * t13 * t40 + 0.2e1 * t16 * t42, 0.2e1 * t1 * t123 - 0.2e1 * t100 * t9 - 0.2e1 * t12 * t40 + 0.2e1 * t16 * t43, 0.2e1 * t1 * t42 + 0.2e1 * t12 * t8 - 0.2e1 * t13 * t9 - 0.2e1 * t2 * t43, -0.2e1 * t1 * t9 + 0.2e1 * t16 * t40 + 0.2e1 * t2 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t332, 0, -t333, 0, -t118, t248, 0, 0, t174 * t345 + (-t136 + (0.2e1 * t319 + t332) * t181) * t178, -t181 * t410 - t340 + (-t100 - t249) * t178, t266 * t175, t73 (t181 * t378 + t342) * t175, 0, -pkin(2) * t100 - t118 * t181 + t178 * t369 - t266 * t401, pkin(2) * t410 + t118 * t178 - t342 * t401 + (-pkin(9) * t333 + t369) * t181, pkin(9) * t73 + t167 * t264 - t168 * t77 - t169 * t410 - t375 * t76 + t285, -pkin(2) * t118 + ((-t77 * t178 - t76 * t181) * qJD(3) + t285) * pkin(9), t178 * t416 + t180 * t227, -t177 * t227 - t180 * t289 + (t414 - t415) * t178, t100 * t389 - t123 * t338 + t178 * t233 + t180 * t340 - t181 * t212, t177 * t289 + t178 * t213 (-qJD(3) * t123 * t177 + t409) * t181 + (-t317 - t272) * t178, t73, -t245 * t100 + t75 * t123 + (-t22 + (pkin(9) * t384 + t72 * t177) * qJD(3)) * t181 + (pkin(9) * t409 + t44 * qJD(3) + t276) * t178, pkin(9) * t227 - t112 * t100 + t74 * t123 - t168 * t45 + t212 * t169 - t21 * t181 + t389 * t48 - t412 * t72, -t112 * t409 + t21 * t390 + t212 * t245 - t22 * t389 - t246 * t75 + t384 * t74 + t412 * t44 - t422 * t45, -t245 * t22 - t112 * t21 + t44 * t75 - t45 * t74 + (t48 * t178 + t375 * t72) * pkin(9), -t119 * t196 - t243 * t66, -t119 * t214 + t196 * t383 + t243 * t274 - t420 * t66, t100 * t119 - t123 * t243 + t168 * t66 + t181 * t196, t214 * t383 + t274 * t420, -t310 * t123 - t383 * t100 + t214 * t181 + (-qJD(3) * t274 + t123 * t419) * t178, t73, t63 * t100 + t107 * t274 + t32 * t123 + t133 * t214 + t168 * t19 - t7 * t181 + t35 * t383 + t420 * t58, -t64 * t100 + t107 * t66 + t35 * t119 + t31 * t123 - t133 * t196 - t168 * t20 - t6 * t181 - t243 * t58, -t7 * t119 + t19 * t243 + t196 * t63 - t20 * t420 - t214 * t64 + t274 * t31 - t32 * t66 + t383 * t6, t107 * t58 + t133 * t35 + t19 * t32 - t20 * t31 - t6 * t64 + t63 * t7, -t12 * t81 - t36 * t43, t12 * t80 - t13 * t81 + t36 * t42 - t37 * t43, t100 * t81 + t12 * t181 - t123 * t36 + t168 * t43, t13 * t80 + t37 * t42, -t100 * t80 - t123 * t37 + t13 * t181 - t168 * t42, t73, t100 * t33 + t11 * t123 + t13 * t98 + t16 * t80 + t168 * t8 - t181 * t2 + t37 * t40 + t42 * t57, -t1 * t181 + t10 * t123 - t100 * t34 - t12 * t98 + t16 * t81 - t168 * t9 - t36 * t40 + t43 * t57, t1 * t80 + t10 * t42 - t11 * t43 + t12 * t33 - t13 * t34 - t2 * t81 + t36 * t8 - t37 * t9, -t1 * t34 - t10 * t9 + t11 * t8 + t16 * t98 + t2 * t33 + t40 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t314, -0.2e1 * t315, 0, t156, 0, 0, t178 * t368, t181 * t368, 0, 0, 0.2e1 * t173 * t330 - 0.2e1 * t302, 0.2e1 * t172 * t318 - 0.4e1 * t177 * t303, 0.2e1 * t178 * t337 + 0.2e1 * t376 * t380, 0.2e1 * t171 * t330 + 0.2e1 * t302, -0.2e1 * t177 * t315 + 0.2e1 * t178 * t334, t156, -0.2e1 * t245 * t168 - 0.2e1 * t181 * t75 + 0.2e1 * (t172 * t372 + t177 * t314) * pkin(9), -0.2e1 * t112 * t168 - 0.2e1 * t181 * t74 + 0.2e1 * (-t172 * t373 + 0.2e1 * t303) * pkin(9), 0.2e1 * t282 * t375 + 0.2e1 * (t177 * t74 - t180 * t75 + (-t112 * t180 - t177 * t245) * qJD(4)) * t178, 0.2e1 * pkin(9) ^ 2 * t330 - 0.2e1 * t112 * t74 - 0.2e1 * t245 * t75, -0.2e1 * t119 * t243, -0.2e1 * t119 * t420 + 0.2e1 * t243 * t383, 0.2e1 * t119 * t168 + 0.2e1 * t181 * t243, 0.2e1 * t383 * t420, 0.2e1 * t310 * t181 + 0.2e1 * (-qJD(3) * t383 - t181 * t419) * t178, t156, -0.2e1 * t32 * t181 + 0.2e1 * t107 * t383 + 0.2e1 * t133 * t310 + 0.2e1 * (t63 * qJD(3) - t133 * t419) * t178, 0.2e1 * t107 * t119 - 0.2e1 * t133 * t243 - 0.2e1 * t168 * t64 - 0.2e1 * t181 * t31, -0.2e1 * t119 * t32 + 0.2e1 * t243 * t63 + 0.2e1 * t31 * t383 - 0.2e1 * t420 * t64, 0.2e1 * t107 * t133 - 0.2e1 * t31 * t64 + 0.2e1 * t32 * t63, -0.2e1 * t81 * t36, 0.2e1 * t36 * t80 - 0.2e1 * t37 * t81, 0.2e1 * t168 * t81 + 0.2e1 * t181 * t36, 0.2e1 * t80 * t37, -0.2e1 * t168 * t80 + 0.2e1 * t181 * t37, t156, -0.2e1 * t11 * t181 + 0.2e1 * t168 * t33 + 0.2e1 * t37 * t98 + 0.2e1 * t57 * t80, -0.2e1 * t10 * t181 - 0.2e1 * t168 * t34 - 0.2e1 * t36 * t98 + 0.2e1 * t57 * t81, 0.2e1 * t10 * t80 - 0.2e1 * t11 * t81 + 0.2e1 * t33 * t36 - 0.2e1 * t34 * t37, -0.2e1 * t10 * t34 + 0.2e1 * t11 * t33 + 0.2e1 * t57 * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t410, 0, -t100, t333, t53, t52, 0, 0, t415, -t213 + t416, t272, t414, -t411, 0, -pkin(3) * t409 - pkin(10) * t272 - t48 * t180 + t373 * t72, -pkin(3) * t212 + pkin(10) * t411 + t276, t246 * t357 - t372 * t44 - t373 * t45 + t286 + (t414 + t211) * pkin(10), -pkin(3) * t48 + ((-t45 * t177 - t44 * t180) * qJD(4) + t286) * pkin(10), -t130 * t196 - t408 * t66, -t130 * t214 + t196 * t287 + t274 * t408 - t309 * t66, t100 * t130 - t123 * t408, t214 * t287 + t274 * t309, -t100 * t287 - t123 * t309, 0, t101 * t100 + t71 * t123 + t166 * t214 + t274 * t358 + t287 * t35 + t309 * t58, -t102 * t100 + t70 * t123 + t35 * t130 - t166 * t196 + t358 * t66 - t408 * t58, t101 * t196 - t102 * t214 - t7 * t130 + t19 * t408 - t20 * t309 + t274 * t70 + t287 * t6 - t71 * t66, t101 * t7 - t102 * t6 + t166 * t35 + t19 * t71 - t20 * t70 + t358 * t58, -t12 * t95 - t43 * t46, t12 * t94 - t13 * t95 + t42 * t46 - t43 * t47, t100 * t95 - t123 * t46, t13 * t94 + t42 * t47, -t100 * t94 - t123 * t47, 0, t100 * t54 + t115 * t13 + t123 * t26 + t16 * t94 + t40 * t47 + t42 * t83, -t100 * t55 - t115 * t12 + t123 * t25 + t16 * t95 - t40 * t46 + t43 * t83, t1 * t94 + t12 * t54 - t13 * t55 - t2 * t95 + t25 * t42 - t26 * t43 + t46 * t8 - t47 * t9, -t1 * t55 + t115 * t16 + t2 * t54 - t25 * t9 + t26 * t8 + t40 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t375, 0, -t168, 0, -t167, pkin(9) * t168, 0, 0, -t103, -0.4e1 * t178 * t331 - t375 * t381, t168 * t177 - t334, t103, t178 * t376 + t337, 0 (pkin(10) * t387 + (-t396 + t398) * t178) * qJD(4) + (t177 * t294 - t363) * qJD(3) (t275 + t364) * qJD(4) + (t180 * t294 + t366) * qJD(3), t234, -pkin(3) * t167 + pkin(10) * t234, -t119 * t408 - t130 * t243, -t119 * t309 - t130 * t420 + t243 * t287 + t383 * t408, t130 * t168 + t181 * t408, t287 * t420 + t309 * t383, -t168 * t287 + t181 * t309, 0, t101 * t168 + t107 * t287 + t133 * t309 + t166 * t420 - t71 * t181 + t358 * t383, -t102 * t168 + t107 * t130 + t119 * t358 - t133 * t408 - t166 * t243 - t70 * t181, t101 * t243 - t102 * t420 - t71 * t119 - t32 * t130 + t287 * t31 - t309 * t64 + t383 * t70 + t408 * t63, t101 * t32 - t102 * t31 + t107 * t166 + t133 * t358 + t63 * t71 - t64 * t70, -t36 * t95 - t46 * t81, t36 * t94 - t37 * t95 + t46 * t80 - t47 * t81, t168 * t95 + t181 * t46, t37 * t94 + t47 * t80, -t168 * t94 + t181 * t47, 0, t115 * t37 + t168 * t54 - t181 * t26 + t47 * t98 + t57 * t94 + t80 * t83, -t115 * t36 - t168 * t55 - t181 * t25 - t46 * t98 + t57 * t95 + t81 * t83, t10 * t94 - t11 * t95 + t25 * t80 - t26 * t81 + t33 * t46 - t34 * t47 + t36 * t54 - t37 * t55, -t10 * t55 + t11 * t54 + t115 * t57 - t25 * t34 + t26 * t33 + t83 * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t331, -0.2e1 * t318, 0, -0.2e1 * t331, 0, 0, t177 * t367, t180 * t367, 0, 0, -0.2e1 * t130 * t408, -0.2e1 * t130 * t309 + 0.2e1 * t287 * t408, 0, 0.2e1 * t287 * t309, 0, 0, 0.2e1 * t166 * t309 + 0.2e1 * t287 * t358, 0.2e1 * t130 * t358 - 0.2e1 * t166 * t408, 0.2e1 * t101 * t408 - 0.2e1 * t102 * t309 - 0.2e1 * t71 * t130 + 0.2e1 * t287 * t70, 0.2e1 * t101 * t71 - 0.2e1 * t102 * t70 + 0.2e1 * t166 * t358, -0.2e1 * t95 * t46, 0.2e1 * t46 * t94 - 0.2e1 * t47 * t95, 0, 0.2e1 * t94 * t47, 0, 0, 0.2e1 * t115 * t47 + 0.2e1 * t83 * t94, -0.2e1 * t115 * t46 + 0.2e1 * t83 * t95, 0.2e1 * t25 * t94 - 0.2e1 * t26 * t95 + 0.2e1 * t46 * t54 - 0.2e1 * t47 * t55, 0.2e1 * t115 * t83 - 0.2e1 * t25 * t55 + 0.2e1 * t26 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t212, 0, -t409, t100, t22, t21, 0, 0, 0, 0, -t196, 0, -t214, t100, t100 * t361 - t123 * t312 + t7, -t100 * t360 - t123 * t313 + t6, t196 * t361 - t214 * t360 - t274 * t313 + t312 * t66 (-t404 * t6 + t406 * t7 + (-t19 * t404 + t20 * t406) * qJD(5)) * pkin(4), 0, 0, -t12, 0, -t13, t100, t121 * t100 + t86 * t123 + t2, -t100 * t122 + t123 * t85 + t1, t12 * t121 - t122 * t13 + t42 * t85 - t43 * t86, -t1 * t122 + t121 * t2 + t8 * t86 - t85 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t412, 0, -t422, t168, t75, t74, 0, 0, 0, 0, -t243, 0, -t420, t168, t181 * t312 + t329 * t397 + t32, t181 * t313 - t328 * t397 + t31 (t404 * (t157 * t178 - t310) + t406 * t243 + (-t406 * t383 + (t390 * t404 + t119) * t404) * qJD(5)) * pkin(4) (-t404 * t31 + t406 * t32 + (-t404 * t63 + t406 * t64) * qJD(5)) * pkin(4), 0, 0, -t36, 0, -t37, t168, t121 * t168 - t86 * t181 + t11, -t122 * t168 - t181 * t85 + t10, t121 * t36 - t122 * t37 + t80 * t85 - t81 * t86, -t10 * t122 + t11 * t121 + t33 * t86 - t34 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t372, 0, -t373, 0, -t357, pkin(10) * t373, 0, 0, 0, 0, -t408, 0, -t309, 0, t71, t70, t130 * t312 - t287 * t313 - t309 * t360 + t361 * t408 (-t404 * t70 + t406 * t71 + (-t101 * t404 + t102 * t406) * qJD(5)) * pkin(4), 0, 0, -t46, 0, -t47, 0, t26, t25, t121 * t46 - t122 * t47 + t85 * t94 - t86 * t95, t121 * t26 - t122 * t25 + t54 * t86 - t55 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t312, -0.2e1 * t313, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t86, 0.2e1 * t85, 0, 0.2e1 * t121 * t86 - 0.2e1 * t122 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t196, 0, -t214, t100, t7, t6, 0, 0, 0, 0, -t12, 0, -t13, t100, -t123 * t356 + t400 * t405 + t2 (-t100 * t176 - t123 * t327) * pkin(5) + t1 (t405 * t12 - t13 * t176 + (t176 * t43 - t405 * t42) * qJD(6)) * pkin(5) (t405 * t2 - t1 * t176 + (-t176 * t8 + t405 * t9) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t243, 0, -t420, t168, t32, t31, 0, 0, 0, 0, -t36, 0, -t37, t168, t181 * t356 + t359 * t405 + t11 (-t168 * t176 + t181 * t327) * pkin(5) + t10 (t405 * t36 - t176 * t37 + (t176 * t81 - t405 * t80) * qJD(6)) * pkin(5) (t405 * t11 - t10 * t176 + (-t176 * t33 + t34 * t405) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t408, 0, -t309, 0, t71, t70, 0, 0, 0, 0, -t46, 0, -t47, 0, t26, t25 (t405 * t46 - t176 * t47 + (t176 * t95 - t405 * t94) * qJD(6)) * pkin(5) (t405 * t26 - t176 * t25 + (-t176 * t54 + t405 * t55) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t312, -t313, 0, 0, 0, 0, 0, 0, 0, 0, t209 + t321, -t311 + t85, 0 (t405 * t86 - t176 * t85 + (-t121 * t176 + t122 * t405) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t321, -0.2e1 * t311, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, -t13, t100, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, 0, -t37, t168, t11, t10, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, 0, -t47, 0, t26, t25, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t85, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t356, -t311, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
