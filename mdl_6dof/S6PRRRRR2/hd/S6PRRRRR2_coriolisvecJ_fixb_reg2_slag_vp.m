% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6PRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRRRR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:45:23
% EndTime: 2019-03-09 00:45:44
% DurationCPUTime: 9.14s
% Computational Cost: add. (14787->566), mult. (35767->768), div. (0->0), fcn. (27555->12), ass. (0->290)
t254 = sin(qJ(2));
t248 = sin(pkin(6));
t342 = qJD(1) * t248;
t316 = t254 * t342;
t253 = sin(qJ(3));
t336 = qJD(3) * t253;
t270 = pkin(3) * t336 - t316;
t245 = qJD(3) + qJD(4);
t257 = cos(qJ(3));
t400 = cos(qJ(4));
t308 = t400 * qJD(3);
t293 = t257 * t308;
t307 = t400 * qJD(4);
t252 = sin(qJ(4));
t355 = t252 * t253;
t166 = t245 * t355 - t257 * t307 - t293;
t354 = t252 * t257;
t210 = t253 * t400 + t354;
t167 = t245 * t210;
t425 = -pkin(4) * t167 - pkin(10) * t166 - t270;
t404 = -pkin(9) - pkin(8);
t319 = qJD(3) * t404;
t212 = t253 * t319;
t317 = t400 * t257;
t275 = t317 - t355;
t258 = cos(qJ(2));
t311 = t258 * t342;
t221 = t404 * t253;
t223 = t404 * t257;
t406 = t400 * t221 + t252 * t223;
t349 = qJD(4) * t406 + t400 * t212 - t275 * t311 + t319 * t354;
t251 = sin(qJ(5));
t297 = pkin(3) * t307;
t217 = qJD(2) * pkin(8) + t316;
t303 = pkin(9) * qJD(2) + t217;
t249 = cos(pkin(6));
t341 = qJD(1) * t253;
t310 = t249 * t341;
t174 = t257 * t303 + t310;
t160 = t252 * t174;
t361 = t249 * t257;
t230 = qJD(1) * t361;
t284 = t303 * t253;
t173 = t230 - t284;
t107 = t173 * t400 - t160;
t339 = qJD(2) * t253;
t313 = t252 * t339;
t201 = -qJD(2) * t317 + t313;
t203 = t210 * qJD(2);
t156 = pkin(4) * t203 + pkin(10) * t201;
t326 = pkin(3) * t339;
t137 = t156 + t326;
t256 = cos(qJ(5));
t68 = -t107 * t251 + t256 * t137;
t424 = -t251 * t297 - t68;
t243 = -pkin(3) * t257 - pkin(2);
t157 = -pkin(4) * t275 - pkin(10) * t210 + t243;
t183 = t252 * t221 - t223 * t400;
t333 = qJD(5) * t256;
t334 = qJD(5) * t251;
t383 = -t157 * t333 + t183 * t334 + t251 * t425 - t349 * t256;
t69 = t256 * t107 + t251 * t137;
t423 = -t256 * t297 + t69;
t422 = -t349 * t251 - t256 * t425;
t172 = t256 * t183;
t379 = t166 * t256;
t421 = pkin(11) * t379 + t167 * pkin(5) + (-t172 + (pkin(11) * t210 - t157) * t251) * qJD(5) + t422;
t357 = t251 * t166;
t274 = t210 * t333 - t357;
t420 = pkin(11) * t274 + t383;
t240 = pkin(3) * t252 + pkin(10);
t398 = -pkin(11) - t240;
t304 = qJD(5) * t398;
t371 = t201 * t251;
t327 = pkin(11) * t371;
t419 = t251 * t304 - t327 - t423;
t370 = t201 * t256;
t289 = t203 * pkin(5) + pkin(11) * t370;
t418 = -t256 * t304 + t289 - t424;
t403 = -pkin(11) - pkin(10);
t318 = qJD(5) * t403;
t162 = qJD(3) * pkin(3) + t173;
t102 = t162 * t400 - t160;
t74 = t256 * t102 + t251 * t156;
t417 = t251 * t318 - t327 - t74;
t73 = -t102 * t251 + t256 * t156;
t416 = t256 * t318 - t289 - t73;
t177 = t203 * t251 - t256 * t245;
t179 = t203 * t256 + t245 * t251;
t250 = sin(qJ(6));
t255 = cos(qJ(6));
t110 = t255 * t177 + t179 * t250;
t280 = t177 * t250 - t255 * t179;
t380 = t110 * t280;
t359 = t250 * t256;
t209 = t251 * t255 + t359;
t329 = qJD(5) + qJD(6);
t165 = t329 * t209;
t415 = t209 * t201 + t165;
t353 = t255 * t256;
t360 = t250 * t251;
t207 = -t353 + t360;
t331 = qJD(6) * t255;
t346 = t207 * t201 - t255 * t333 - t256 * t331 + t329 * t360;
t414 = -t110 ^ 2 + t280 ^ 2;
t195 = qJD(2) * t243 - t311;
t124 = t201 * pkin(4) - t203 * pkin(10) + t195;
t340 = qJD(2) * t248;
t306 = qJD(1) * t340;
t292 = t258 * t306;
t345 = qJD(3) * t230 + t257 * t292;
t128 = -qJD(3) * t284 + t345;
t129 = -qJD(3) * t174 - t253 * t292;
t335 = qJD(4) * t252;
t47 = t128 * t400 + t252 * t129 + t162 * t307 - t174 * t335;
t330 = qJD(2) * qJD(3);
t305 = t253 * t330;
t337 = qJD(2) * t257;
t152 = -(t308 + t307) * t337 + qJD(4) * t313 + t252 * t305;
t153 = t167 * qJD(2);
t198 = pkin(3) * t305 + t254 * t306;
t83 = pkin(4) * t153 + pkin(10) * t152 + t198;
t161 = t400 * t174;
t103 = t252 * t162 + t161;
t94 = t245 * pkin(10) + t103;
t16 = t124 * t333 + t251 * t83 + t256 * t47 - t334 * t94;
t320 = -t251 * t152 + t203 * t333 + t245 * t334;
t12 = -pkin(11) * t320 + t16;
t66 = t124 * t251 + t256 * t94;
t51 = -pkin(11) * t177 + t66;
t392 = t255 * t51;
t197 = qJD(5) + t201;
t65 = t256 * t124 - t251 * t94;
t50 = -pkin(11) * t179 + t65;
t42 = pkin(5) * t197 + t50;
t21 = t250 * t42 + t392;
t100 = t256 * t152 + t203 * t334 - t245 * t333;
t17 = -qJD(5) * t66 - t251 * t47 + t256 * t83;
t9 = t153 * pkin(5) + t100 * pkin(11) + t17;
t4 = -qJD(6) * t21 - t250 * t12 + t255 * t9;
t93 = -t245 * pkin(4) - t102;
t77 = t177 * pkin(5) + t93;
t413 = t77 * t280 + t4;
t194 = qJD(6) + t197;
t332 = qJD(6) * t250;
t36 = t255 * t100 + t177 * t331 + t179 * t332 + t250 * t320;
t412 = t110 * t194 - t36;
t3 = (qJD(6) * t42 + t12) * t255 + t250 * t9 - t51 * t332;
t411 = t77 * t110 - t3;
t264 = qJD(6) * t280 + t250 * t100 - t255 * t320;
t410 = -t194 * t280 + t264;
t299 = t252 * t128 - t400 * t129;
t48 = qJD(4) * t103 + t299;
t409 = t48 * t251 + t93 * t333;
t408 = t48 * t256 - t93 * t334;
t106 = t173 * t252 + t161;
t288 = pkin(3) * t335 - t106;
t348 = qJD(4) * t183 - t210 * t311 + t252 * t212 - t293 * t404;
t105 = t251 * t157 + t172;
t407 = (t334 + t371) * pkin(5);
t405 = -t251 * t65 + t256 * t66;
t104 = t256 * t157 - t183 * t251;
t367 = t210 * t256;
t78 = -pkin(5) * t275 - pkin(11) * t367 + t104;
t368 = t210 * t251;
t88 = -pkin(11) * t368 + t105;
t40 = -t250 * t88 + t255 * t78;
t402 = qJD(6) * t40 + t250 * t421 - t420 * t255;
t41 = t250 * t78 + t255 * t88;
t401 = -qJD(6) * t41 + t420 * t250 + t255 * t421;
t399 = t256 * pkin(5);
t205 = t398 * t251;
t244 = t256 * pkin(11);
t206 = t240 * t256 + t244;
t155 = t205 * t250 + t206 * t255;
t397 = qJD(6) * t155 + t250 * t419 + t255 * t418;
t154 = t205 * t255 - t206 * t250;
t396 = -qJD(6) * t154 + t250 * t418 - t255 * t419;
t395 = qJD(2) * pkin(2);
t15 = t16 * t256;
t394 = t250 * t51;
t364 = t248 * t254;
t199 = -t253 * t364 + t361;
t200 = t249 * t253 + t257 * t364;
t276 = t199 * t400 - t252 * t200;
t390 = t48 * t276;
t389 = t48 * t406;
t220 = t403 * t251;
t222 = pkin(10) * t256 + t244;
t180 = t220 * t255 - t222 * t250;
t387 = qJD(6) * t180 + t250 * t416 + t255 * t417;
t182 = t220 * t250 + t222 * t255;
t386 = -qJD(6) * t182 - t250 * t417 + t255 * t416;
t385 = pkin(5) * t274 + t348;
t384 = t407 + t288;
t382 = -qJD(5) * t105 + t422;
t381 = t100 * t251;
t125 = t153 * t275;
t378 = t177 * t197;
t377 = t177 * t251;
t376 = t179 * t177;
t375 = t179 * t197;
t374 = t179 * t256;
t373 = t194 * t203;
t372 = t197 * t203;
t369 = t203 * t201;
t366 = t217 * t253;
t365 = t217 * t257;
t363 = t248 * t258;
t260 = qJD(2) ^ 2;
t362 = t248 * t260;
t358 = t251 * t153;
t352 = t256 * t153;
t259 = qJD(3) ^ 2;
t351 = t259 * t253;
t350 = t259 * t257;
t246 = t253 ^ 2;
t247 = t257 ^ 2;
t343 = t246 - t247;
t338 = qJD(2) * t254;
t328 = -t65 * t370 - t66 * t371 + t15;
t322 = t254 * t362;
t321 = t253 * t260 * t257;
t315 = t248 * t338;
t314 = t258 * t340;
t312 = t210 * t334;
t298 = t197 * t256;
t296 = t253 * t314;
t295 = t257 * t314;
t241 = -pkin(3) * t400 - pkin(4);
t294 = t66 * t203 + t409;
t291 = t257 * t305;
t290 = -t103 + t407;
t287 = t320 * t256;
t218 = -t311 - t395;
t286 = -t218 - t311;
t285 = t251 * t66 + t256 * t65;
t281 = -t153 * t240 + t201 * t93;
t139 = t252 * t199 + t200 * t400;
t122 = -t139 * t251 - t256 * t363;
t278 = -t139 * t256 + t251 * t363;
t70 = t122 * t255 + t250 * t278;
t71 = t122 * t250 - t255 * t278;
t279 = -t65 * t203 - t408;
t277 = -t195 * t203 - t299;
t187 = t310 + t365;
t273 = -t312 - t379;
t20 = t255 * t42 - t394;
t35 = pkin(5) * t320 + t48;
t272 = -t20 * t203 + t35 * t207 + t415 * t77;
t271 = t21 * t203 + t35 * t209 - t346 * t77;
t268 = t20 * t346 - t3 * t207 - t4 * t209 - t21 * t415;
t266 = qJD(3) * (-t286 - t395);
t265 = -qJD(5) * t285 - t17 * t251;
t263 = t265 + t15;
t147 = -t217 * t336 + t345;
t148 = -qJD(3) * t365 + (-qJD(3) * t249 - t314) * t341;
t186 = t230 - t366;
t262 = t147 * t257 - t148 * t253 + (-t186 * t257 - t187 * t253) * qJD(3);
t261 = t195 * t201 - t47;
t242 = -pkin(4) - t399;
t219 = t241 - t399;
t171 = -qJD(3) * t200 - t296;
t170 = qJD(3) * t199 + t295;
t146 = t207 * t210;
t145 = t209 * t210;
t134 = pkin(5) * t368 - t406;
t126 = -t201 ^ 2 + t203 ^ 2;
t119 = t201 * t245 - t152;
t76 = qJD(4) * t139 + t252 * t170 - t171 * t400;
t75 = qJD(4) * t276 + t170 * t400 + t252 * t171;
t61 = -t179 * t203 + t197 * t298 + t358;
t60 = -t197 ^ 2 * t251 + t177 * t203 + t352;
t58 = t197 * t377 - t287;
t57 = t179 * t298 - t381;
t53 = -t166 * t359 - t250 * t312 - t332 * t368 + (t329 * t367 - t357) * t255;
t52 = t165 * t210 + t166 * t353 - t250 * t357;
t44 = qJD(5) * t278 - t251 * t75 + t256 * t315;
t43 = qJD(5) * t122 + t251 * t315 + t256 * t75;
t34 = t110 * t203 - t207 * t153 - t194 * t415;
t33 = t209 * t153 - t194 * t346 + t203 * t280;
t24 = (-t100 - t378) * t256 + (-t320 - t375) * t251;
t23 = t255 * t50 - t394;
t22 = -t250 * t50 - t392;
t14 = t110 * t415 - t207 * t264;
t13 = -t36 * t209 + t280 * t346;
t11 = -qJD(6) * t71 - t250 * t43 + t255 * t44;
t10 = qJD(6) * t70 + t250 * t44 + t255 * t43;
t5 = t110 * t346 + t36 * t207 + t209 * t264 + t280 * t415;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t322, -t258 * t362, 0, 0, 0, 0, 0, 0, 0, 0, -t257 * t322 + (t171 - t296) * qJD(3), t253 * t322 + (-t170 - t295) * qJD(3) (t170 * t257 - t171 * t253 + (-t199 * t257 - t200 * t253) * qJD(3)) * qJD(2), t147 * t200 + t148 * t199 + t187 * t170 + t186 * t171 + (t218 - t311) * t315, 0, 0, 0, 0, 0, 0, -t76 * t245 + (-t153 * t258 + t201 * t338) * t248, -t75 * t245 + (t152 * t258 + t203 * t338) * t248, -t139 * t153 + t152 * t276 - t201 * t75 + t203 * t76, -t102 * t76 + t103 * t75 - t390 + t47 * t139 + (t195 * t338 - t198 * t258) * t248, 0, 0, 0, 0, 0, 0, t122 * t153 + t76 * t177 + t44 * t197 - t276 * t320, t100 * t276 + t153 * t278 + t179 * t76 - t197 * t43, t122 * t100 - t43 * t177 - t44 * t179 + t278 * t320, t122 * t17 - t16 * t278 + t43 * t66 + t44 * t65 + t76 * t93 - t390, 0, 0, 0, 0, 0, 0, t11 * t194 + t110 * t76 + t153 * t70 + t264 * t276, -t10 * t194 - t153 * t71 + t276 * t36 - t280 * t76, -t10 * t110 + t11 * t280 + t264 * t71 + t36 * t70, t10 * t21 + t11 * t20 - t276 * t35 + t3 * t71 + t4 * t70 + t76 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t291, -0.2e1 * t343 * t330, t350, -0.2e1 * t291, -t351, 0, -pkin(8) * t350 + t253 * t266, pkin(8) * t351 + t257 * t266 (-t246 - t247) * t292 + t262 ((t186 * t253 - t187 * t257) * t258 + (-t218 - t395) * t254) * t342 + t262 * pkin(8), -t152 * t210 - t166 * t203, -t152 * t275 - t153 * t210 + t166 * t201 - t167 * t203, -t166 * t245, t167 * t201 - t125, -t167 * t245, 0, t243 * t153 + t195 * t167 - t198 * t275 + t201 * t270 - t245 * t348, -t243 * t152 - t195 * t166 + t198 * t210 + t203 * t270 - t245 * t349, t102 * t166 - t103 * t167 + t152 * t406 - t183 * t153 - t201 * t349 + t203 * t348 + t48 * t210 + t275 * t47, -t102 * t348 + t103 * t349 + t47 * t183 + t195 * t270 + t198 * t243 - t389, -t100 * t367 + t179 * t273 (t177 * t256 + t179 * t251) * t166 + (-t287 + t381 + (-t374 + t377) * qJD(5)) * t210, t100 * t275 + t179 * t167 + t197 * t273 + t210 * t352, t177 * t274 + t320 * t368, -t177 * t167 - t197 * t274 - t210 * t358 + t275 * t320, t167 * t197 - t125, t104 * t153 + t65 * t167 - t17 * t275 + t348 * t177 + t382 * t197 + t210 * t409 - t320 * t406 - t93 * t357, t100 * t406 - t105 * t153 + t16 * t275 - t66 * t167 + t348 * t179 + t383 * t197 + t210 * t408 - t93 * t379, -t105 * t320 + t104 * t100 - t382 * t179 + t383 * t177 + t285 * t166 + (-qJD(5) * t405 - t16 * t251 - t17 * t256) * t210, t17 * t104 + t16 * t105 + t348 * t93 + t382 * t65 - t383 * t66 - t389, t146 * t36 + t280 * t52, t110 * t52 + t145 * t36 - t146 * t264 + t280 * t53, -t146 * t153 - t167 * t280 - t194 * t52 + t275 * t36, t110 * t53 - t145 * t264, -t110 * t167 - t145 * t153 - t194 * t53 - t264 * t275, t167 * t194 - t125, t110 * t385 - t134 * t264 + t35 * t145 + t40 * t153 + t20 * t167 + t194 * t401 - t275 * t4 + t77 * t53, -t134 * t36 - t35 * t146 - t41 * t153 - t21 * t167 - t194 * t402 + t275 * t3 - t280 * t385 - t77 * t52, -t110 * t402 - t3 * t145 + t4 * t146 + t20 * t52 - t21 * t53 + t264 * t41 + t280 * t401 + t40 * t36, t35 * t134 + t20 * t401 + t21 * t402 + t3 * t41 + t385 * t77 + t4 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t321, t343 * t260, 0, t321, 0, 0, t286 * t339, -t218 * t337 + (t186 + t366) * qJD(3) - t345, 0, 0, t369, t126, t119, -t369, 0, 0, -t201 * t326 + t106 * t245 + (-t161 + (-pkin(3) * t245 - t162) * t252) * qJD(4) + t277, t107 * t245 + (-t203 * t339 - t245 * t307) * pkin(3) + t261 (t103 - t106) * t203 + (-t102 + t107) * t201 + (t400 * t152 - t153 * t252 + (-t201 * t400 + t203 * t252) * qJD(4)) * pkin(3), t102 * t106 - t103 * t107 + (-t195 * t339 - t400 * t48 + t252 * t47 + (-t102 * t252 + t103 * t400) * qJD(4)) * pkin(3), t57, t24, t61, t58, t60, -t372, t241 * t320 + t281 * t251 + t288 * t177 + (-t240 * t333 + t424) * t197 + t279, -t241 * t100 + t281 * t256 + t288 * t179 + (t240 * t334 + t423) * t197 + t294, t69 * t177 + t68 * t179 + (-t177 * t297 - t240 * t320 + (t179 * t240 - t65) * qJD(5)) * t256 + (t179 * t297 - t240 * t100 - t17 + (t177 * t240 - t66) * qJD(5)) * t251 + t328, -t93 * t106 + t48 * t241 - t65 * t68 - t66 * t69 + (t252 * t93 + t400 * t405) * qJD(4) * pkin(3) + t263 * t240, t13, t5, t33, t14, t34, -t373, t110 * t384 + t154 * t153 - t194 * t397 - t219 * t264 + t272, -t155 * t153 + t194 * t396 - t219 * t36 - t280 * t384 + t271, t110 * t396 + t154 * t36 + t155 * t264 - t280 * t397 + t268, t4 * t154 + t3 * t155 - t20 * t397 - t21 * t396 + t35 * t219 + t384 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t369, t126, t119, -t369, 0, 0, t277 + (-qJD(4) + t245) * t103, t102 * t245 + t261, 0, 0, t57, t24, t61, t58, t60, -t372, -pkin(4) * t320 - t73 * t197 - t103 * t177 + t93 * t371 + (-t197 * t333 - t358) * pkin(10) + t279, t93 * t370 + pkin(4) * t100 - t103 * t179 + t74 * t197 + (t197 * t334 - t352) * pkin(10) + t294, t74 * t177 + t73 * t179 + (-t287 - t381 + (t374 + t377) * qJD(5)) * pkin(10) + t265 + t328, -t48 * pkin(4) + pkin(10) * t263 - t93 * t103 - t65 * t73 - t66 * t74, t13, t5, t33, t14, t34, -t373, t110 * t290 + t180 * t153 + t194 * t386 - t242 * t264 + t272, -t182 * t153 - t194 * t387 - t242 * t36 - t280 * t290 + t271, -t110 * t387 + t180 * t36 + t182 * t264 + t280 * t386 + t268, t4 * t180 + t3 * t182 + t20 * t386 + t21 * t387 + t35 * t242 + t290 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t376, -t177 ^ 2 + t179 ^ 2, -t100 + t378, -t376, -t320 + t375, t153, -t93 * t179 + t66 * t197 + t17, t177 * t93 + t197 * t65 - t16, 0, 0, -t380, t414, t412, t380, t410, t153, -t22 * t194 + (-t110 * t179 + t153 * t255 - t194 * t332) * pkin(5) + t413, t23 * t194 + (-t153 * t250 + t179 * t280 - t194 * t331) * pkin(5) + t411, -t21 * t280 + t23 * t110 - t20 * t110 - t22 * t280 + (t250 * t264 + t255 * t36 + (-t110 * t255 - t250 * t280) * qJD(6)) * pkin(5), -t20 * t22 - t21 * t23 + (-t179 * t77 + t250 * t3 + t255 * t4 + (-t20 * t250 + t21 * t255) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t380, t414, t412, t380, t410, t153, t21 * t194 + t413, t20 * t194 + t411, 0, 0;];
tauc_reg  = t1;