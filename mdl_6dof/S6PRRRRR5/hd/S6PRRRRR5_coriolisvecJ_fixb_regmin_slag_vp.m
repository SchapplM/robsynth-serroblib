% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRRRR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:10:02
% EndTime: 2019-03-09 01:10:25
% DurationCPUTime: 9.63s
% Computational Cost: add. (8500->537), mult. (22572->788), div. (0->0), fcn. (18884->14), ass. (0->267)
t226 = sin(pkin(7));
t233 = sin(qJ(3));
t354 = t226 * t233;
t217 = pkin(9) * t354;
t228 = cos(pkin(7));
t238 = cos(qJ(3));
t239 = cos(qJ(2));
t342 = t238 * t239;
t234 = sin(qJ(2));
t347 = t233 * t234;
t257 = -t228 * t347 + t342;
t227 = sin(pkin(6));
t331 = qJD(1) * t227;
t350 = t228 * t238;
t338 = t257 * t331 - (pkin(2) * t350 - t217) * qJD(3);
t271 = pkin(3) * t233 - pkin(10) * t238;
t256 = t271 * qJD(3);
t309 = t234 * t331;
t406 = (t256 - t309) * t226;
t326 = qJD(2) * t238;
t307 = t226 * t326;
t405 = qJD(4) - t307;
t353 = t226 * t238;
t314 = pkin(9) * t353;
t171 = t314 + (pkin(2) * t233 + pkin(10)) * t228;
t272 = -pkin(3) * t238 - pkin(10) * t233;
t172 = (-pkin(2) + t272) * t226;
t232 = sin(qJ(4));
t237 = cos(qJ(4));
t320 = qJD(4) * t237;
t322 = qJD(4) * t232;
t404 = -t171 * t322 + t172 * t320 + t406 * t232 - t237 * t338;
t345 = t234 * t238;
t346 = t233 * t239;
t259 = t228 * t345 + t346;
t351 = t228 * t233;
t336 = -t259 * t331 + (pkin(2) * t351 + t314) * qJD(3);
t324 = qJD(3) * t233;
t305 = t226 * t324;
t403 = -pkin(11) * t305 - t404;
t183 = -t237 * t228 + t232 * t354;
t325 = qJD(3) * t226;
t304 = t238 * t325;
t143 = -qJD(4) * t183 + t237 * t304;
t184 = t228 * t232 + t237 * t354;
t144 = qJD(4) * t184 + t232 * t304;
t402 = -pkin(4) * t144 + pkin(11) * t143 - t336;
t229 = cos(pkin(6));
t330 = qJD(1) * t229;
t310 = t226 * t330;
t329 = qJD(2) * t226;
t189 = pkin(9) * t329 + t309;
t376 = qJD(2) * pkin(2);
t203 = t239 * t331 + t376;
t384 = t238 * t189 + t203 * t351;
t118 = t233 * t310 + t384;
t401 = -t118 + t405 * (pkin(4) * t232 - pkin(11) * t237);
t231 = sin(qJ(5));
t236 = cos(qJ(5));
t308 = t233 * t329;
t155 = t231 * t237 * t307 - t236 * t308;
t395 = -t231 * t320 + t155;
t343 = t237 * t238;
t156 = (t231 * t233 + t236 * t343) * t329;
t267 = t236 * t320 - t156;
t319 = qJD(5) * t231;
t400 = -t232 * t319 + t267;
t327 = qJD(2) * t228;
t287 = qJD(3) + t327;
t165 = t232 * t287 + t237 * t308;
t129 = t165 * t231 - t236 * t405;
t103 = pkin(10) * t287 + t118;
t216 = t228 * t330;
t127 = t216 + (qJD(2) * t272 - t203) * t226;
t48 = -t232 * t103 + t127 * t237;
t42 = -pkin(4) * t405 - t48;
t36 = pkin(5) * t129 + t42;
t131 = t165 * t236 + t231 * t405;
t230 = sin(qJ(6));
t235 = cos(qJ(6));
t62 = t235 * t129 + t131 * t230;
t399 = t36 * t62;
t264 = t129 * t230 - t235 * t131;
t398 = t264 * t62;
t397 = t171 * t320 + t172 * t322 - t338 * t232 - t237 * t406;
t318 = qJD(5) * t236;
t396 = t232 * t318 - t395;
t382 = qJD(5) + qJD(6);
t393 = -t232 * t308 + t237 * t287;
t394 = -t393 + t382;
t392 = t264 ^ 2 - t62 ^ 2;
t315 = qJD(2) * qJD(3);
t297 = t226 * t315;
t278 = t238 * t297;
t133 = qJD(4) * t165 + t232 * t278;
t132 = qJD(4) * t393 + t237 * t278;
t279 = t233 * t297;
t59 = t236 * t132 - t165 * t319 + t231 * t279 + t318 * t405;
t136 = (t256 + t309) * t329;
t180 = t233 * t189;
t245 = t257 * qJD(2);
t280 = t229 * t304;
t76 = (t203 * t350 - t180) * qJD(3) + (t227 * t245 + t280) * qJD(1);
t249 = -t103 * t322 + t127 * t320 + t232 * t136 + t237 * t76;
t22 = pkin(11) * t279 + t249;
t49 = t237 * t103 + t232 * t127;
t43 = pkin(11) * t405 + t49;
t357 = t203 * t228;
t117 = t238 * (t310 + t357) - t180;
t102 = -pkin(3) * t287 - t117;
t51 = -pkin(4) * t393 - t165 * pkin(11) + t102;
t25 = t231 * t51 + t236 * t43;
t246 = t259 * qJD(2);
t281 = t229 * t305;
t77 = t384 * qJD(3) + (t227 * t246 + t281) * qJD(1);
t37 = pkin(4) * t133 - pkin(11) * t132 + t77;
t7 = -qJD(5) * t25 - t22 * t231 + t236 * t37;
t4 = pkin(5) * t133 - pkin(12) * t59 + t7;
t6 = t236 * t22 + t231 * t37 + t51 * t318 - t319 * t43;
t60 = qJD(5) * t131 + t132 * t231 - t236 * t279;
t5 = -pkin(12) * t60 + t6;
t158 = qJD(5) - t393;
t24 = -t231 * t43 + t236 * t51;
t19 = -pkin(12) * t131 + t24;
t15 = pkin(5) * t158 + t19;
t20 = -pkin(12) * t129 + t25;
t375 = t20 * t235;
t9 = t15 * t230 + t375;
t2 = -qJD(6) * t9 - t230 * t5 + t235 * t4;
t391 = t36 * t264 + t2;
t149 = qJD(6) + t158;
t316 = qJD(6) * t235;
t317 = qJD(6) * t230;
t16 = -t129 * t316 - t131 * t317 - t230 * t60 + t235 * t59;
t390 = t149 * t62 + t16;
t242 = qJD(6) * t264 - t230 * t59 - t235 * t60;
t389 = -t149 * t264 + t242;
t170 = t217 + (-pkin(2) * t238 - pkin(3)) * t228;
t108 = pkin(4) * t183 - pkin(11) * t184 + t170;
t337 = t237 * t171 + t232 * t172;
t110 = -pkin(11) * t353 + t337;
t388 = -t108 * t318 + t110 * t319 + t402 * t231 + t236 * t403;
t368 = -pkin(4) * t305 + t397;
t377 = t231 * t108 + t236 * t110;
t387 = -qJD(5) * t377 + t231 * t403 - t402 * t236;
t379 = pkin(10) * t231;
t386 = t236 * t401 + t322 * t379;
t194 = t230 * t236 + t231 * t235;
t177 = t194 * t232;
t385 = t233 * t238;
t209 = -pkin(4) * t237 - pkin(11) * t232 - pkin(3);
t173 = t271 * t329;
t341 = t237 * t117 + t232 * t173;
t73 = pkin(11) * t308 + t341;
t383 = -t209 * t318 - t231 * t401 + t236 * t73;
t18 = t20 * t317;
t296 = qJD(6) * t15 + t5;
t1 = t230 * t4 + t235 * t296 - t18;
t240 = qJD(2) ^ 2;
t381 = pkin(11) + pkin(12);
t380 = pkin(5) * t232;
t289 = t103 * t320 + t127 * t322 - t237 * t136 + t232 * t76;
t23 = -pkin(4) * t279 + t289;
t374 = t23 * t231;
t373 = t23 * t236;
t372 = t231 * t59;
t293 = -t232 * t117 + t173 * t237;
t72 = -pkin(4) * t308 - t293;
t371 = pkin(5) * t396 + pkin(10) * t320 - t72;
t313 = t231 * t353;
t79 = -qJD(5) * t313 + t143 * t231 + t184 * t318 - t236 * t305;
t370 = pkin(5) * t79 + t368;
t116 = pkin(4) * t165 - pkin(11) * t393;
t369 = t231 * t116 + t236 * t48;
t193 = t230 * t231 - t235 * t236;
t367 = t155 * t230 - t156 * t235 - t177 * t382 - t193 * t320;
t348 = t232 * t236;
t349 = t231 * t232;
t366 = -t317 * t349 + (t348 * t382 - t395) * t235 + t400 * t230;
t365 = t129 * t158;
t364 = t131 * t158;
t363 = t133 * t231;
t362 = t133 * t236;
t361 = t133 * t237;
t360 = t393 * t405;
t359 = t393 * t231;
t358 = t165 * t405;
t252 = t405 * t232;
t356 = t405 * t237;
t223 = t226 ^ 2;
t355 = t223 * t240;
t352 = t227 * t240;
t344 = t236 * t237;
t340 = t394 * t193;
t339 = t394 * t194;
t219 = pkin(10) * t344;
t333 = t231 * t209 + t219;
t332 = t233 ^ 2 - t238 ^ 2;
t328 = qJD(2) * t227;
t323 = qJD(3) * t237;
t321 = qJD(4) * t236;
t312 = t234 * t352;
t311 = qJD(5) * t381;
t306 = t234 * t328;
t303 = t158 * t319;
t301 = t226 * t228 * t240;
t294 = t236 * t108 - t110 * t231;
t291 = -t232 * t171 + t172 * t237;
t290 = t158 * t236;
t288 = 0.2e1 * t223 * t315;
t286 = qJD(3) + 0.2e1 * t327;
t285 = t223 * t312;
t283 = t226 * t306;
t277 = -t49 + (t319 - t359) * pkin(5);
t145 = t184 * t231 + t236 * t353;
t38 = -pkin(12) * t145 + t377;
t78 = -qJD(5) * t145 + t143 * t236 + t231 * t305;
t276 = -pkin(5) * t144 + pkin(12) * t78 + qJD(6) * t38 - t387;
t146 = t184 * t236 - t313;
t33 = pkin(5) * t183 - pkin(12) * t146 + t294;
t275 = pkin(12) * t79 - qJD(6) * t33 + t388;
t191 = t236 * t209;
t138 = -pkin(12) * t348 + t191 + (-pkin(5) - t379) * t237;
t274 = -qJD(6) * t138 - (-t232 * t321 - t237 * t319) * pkin(10) + t383 + t396 * pkin(12);
t147 = -pkin(12) * t349 + t333;
t273 = -pkin(12) * t156 + qJD(6) * t147 - t231 * t73 + t307 * t380 - (-pkin(12) * t344 + t380) * qJD(4) - (-t219 + (pkin(12) * t232 - t209) * t231) * qJD(5) - t386;
t213 = t381 * t231;
t269 = -pkin(12) * t359 + qJD(6) * t213 + t231 * t311 + t369;
t113 = t236 * t116;
t214 = t381 * t236;
t268 = pkin(5) * t165 + qJD(6) * t214 - t231 * t48 + t113 + (-pkin(12) * t393 + t311) * t236;
t109 = pkin(4) * t353 - t291;
t258 = t228 * t346 + t345;
t140 = t227 * t258 + t229 * t354;
t182 = -t226 * t227 * t239 + t228 * t229;
t107 = t140 * t237 + t182 * t232;
t260 = t228 * t342 - t347;
t139 = -t227 * t260 - t229 * t353;
t57 = -t107 * t231 + t139 * t236;
t58 = t107 * t236 + t139 * t231;
t266 = -t230 * t58 + t235 * t57;
t265 = t230 * t57 + t235 * t58;
t106 = t140 * t232 - t182 * t237;
t86 = t235 * t145 + t146 * t230;
t87 = -t145 * t230 + t146 * t235;
t157 = -t203 * t226 + t216;
t261 = t157 * t226 - t223 * t376;
t254 = -t158 * t318 - t363;
t250 = -pkin(11) * t133 + t158 * t42;
t244 = qJD(1) * t228 * t306 + qJD(3) * t189;
t241 = -t157 * t329 - qJD(3) * t357 + (-t229 * t325 - t239 * t328) * qJD(1);
t222 = -pkin(5) * t236 - pkin(4);
t199 = (pkin(5) * t231 + pkin(10)) * t232;
t178 = t193 * t232;
t111 = t133 * t183;
t96 = t280 + (qJD(3) * t260 + t245) * t227;
t95 = t281 + (qJD(3) * t258 + t246) * t227;
t66 = pkin(5) * t145 + t109;
t40 = -qJD(4) * t106 + t232 * t283 + t237 * t96;
t39 = qJD(4) * t107 + t232 * t96 - t237 * t283;
t27 = qJD(6) * t87 + t230 * t78 + t235 * t79;
t26 = -qJD(6) * t86 - t230 * t79 + t235 * t78;
t14 = qJD(5) * t57 + t231 * t95 + t236 * t40;
t13 = -qJD(5) * t58 - t231 * t40 + t236 * t95;
t12 = pkin(5) * t60 + t23;
t8 = t15 * t235 - t20 * t230;
t3 = [0, 0, -t312, -t239 * t352, 0, 0, 0, 0, 0, t182 * t279 - t238 * t285 - t287 * t95, t182 * t278 + t233 * t285 - t287 * t96, 0, 0, 0, 0, 0, -t106 * t279 + t133 * t139 - t39 * t405 - t393 * t95, -t107 * t279 + t132 * t139 + t165 * t95 - t40 * t405, 0, 0, 0, 0, 0, t106 * t60 + t129 * t39 + t13 * t158 + t57 * t133, t106 * t59 + t131 * t39 - t58 * t133 - t14 * t158, 0, 0, 0, 0, 0 (-qJD(6) * t265 + t13 * t235 - t14 * t230) * t149 + t266 * t133 + t39 * t62 - t106 * t242 -(qJD(6) * t266 + t13 * t230 + t14 * t235) * t149 - t265 * t133 - t39 * t264 + t106 * t16; 0, 0, 0, 0, t288 * t385, -t332 * t288, t286 * t304, -t286 * t305, 0 (-qJD(2) * t336 - t77) * t228 + (t233 * t261 - t336) * qJD(3) (qJD(2) * t338 - t76) * t228 + (t238 * t261 + t338) * qJD(3), t132 * t184 + t143 * t165, -t132 * t183 - t133 * t184 + t143 * t393 - t144 * t165, t143 * t405 + (-t132 * t238 + (qJD(2) * t184 + t165) * t324) * t226, -t144 * t405 + (t133 * t238 + (-qJD(2) * t183 + t393) * t324) * t226 (-t223 * t326 + t226 * t405) * t324, t102 * t144 + t170 * t133 + t77 * t183 - t397 * t405 - t336 * t393 + (t289 * t238 + (qJD(2) * t291 + t48) * t324) * t226, t102 * t143 + t170 * t132 + t77 * t184 - t404 * t405 + t336 * t165 + (t249 * t238 + (-qJD(2) * t337 - t49) * t324) * t226, t131 * t78 + t146 * t59, -t129 * t78 - t131 * t79 - t145 * t59 - t146 * t60, t131 * t144 + t133 * t146 + t158 * t78 + t183 * t59, -t129 * t144 - t133 * t145 - t158 * t79 - t183 * t60, t144 * t158 + t111, t109 * t60 + t368 * t129 + t294 * t133 + t24 * t144 + t23 * t145 + t158 * t387 + t7 * t183 + t42 * t79, t109 * t59 + t368 * t131 - t377 * t133 - t25 * t144 + t23 * t146 + t158 * t388 - t6 * t183 + t42 * t78, t16 * t87 - t26 * t264, -t16 * t86 + t242 * t87 - t26 * t62 + t264 * t27, t133 * t87 - t144 * t264 + t149 * t26 + t16 * t183, -t133 * t86 - t144 * t62 - t149 * t27 + t183 * t242, t144 * t149 + t111 (-t230 * t38 + t235 * t33) * t133 + t2 * t183 + t8 * t144 - t66 * t242 + t12 * t86 + t36 * t27 + t370 * t62 + (t230 * t275 - t235 * t276) * t149 -(t230 * t33 + t235 * t38) * t133 - t1 * t183 - t9 * t144 + t66 * t16 + t12 * t87 + t36 * t26 - t370 * t264 + (t230 * t276 + t235 * t275) * t149; 0, 0, 0, 0, -t355 * t385, t332 * t355, -t238 * t301, t233 * t301, 0, t118 * t287 + t233 * t241 - t238 * t244, t117 * t287 + t233 * t244 + t238 * t241, t132 * t232 + t165 * t356 (t132 + t360) * t237 + (-t133 - t358) * t232, t405 * t320 + (-t405 * t343 + (t232 * qJD(3) - t165) * t233) * t329, -t405 * t322 + (t238 * t252 + (-t393 + t323) * t233) * t329, -t405 * t308, -pkin(3) * t133 - t77 * t237 - t293 * t405 + t118 * t393 + (-pkin(10) * t356 + t102 * t232) * qJD(4) + (-t233 * t48 + (-pkin(10) * t324 - t102 * t238) * t232) * t329, -pkin(3) * t132 + t77 * t232 - t118 * t165 + t341 * t405 + (pkin(10) * t252 + t102 * t237) * qJD(4) + (-t102 * t343 + (-pkin(10) * t323 + t49) * t233) * t329, t131 * t400 + t59 * t348, t129 * t156 + t131 * t155 + (-t129 * t236 - t131 * t231) * t320 + (-t372 - t236 * t60 + (t129 * t231 - t131 * t236) * qJD(5)) * t232, -t237 * t59 + t267 * t158 + (t131 * t405 - t303 + t362) * t232, t237 * t60 + t395 * t158 + (-t129 * t405 + t254) * t232, t158 * t252 - t361, -t72 * t129 + t191 * t133 - t42 * t155 + ((-qJD(5) * t209 + t73) * t231 + t386) * t158 + (t42 * t231 * qJD(4) - t7 + (qJD(4) * t129 + t254) * pkin(10)) * t237 + (pkin(10) * t60 + t24 * t405 + t318 * t42 + t374) * t232, -t333 * t133 - t72 * t131 - t42 * t156 + t383 * t158 + (t42 * t321 + t6 + (qJD(4) * t131 + t303) * pkin(10)) * t237 + (-t42 * t319 + t373 - t405 * t25 + (t158 * t321 + t59) * pkin(10)) * t232, -t16 * t178 - t264 * t367, -t16 * t177 - t178 * t242 + t264 * t366 - t367 * t62, -t133 * t178 + t149 * t367 - t16 * t237 - t252 * t264, -t133 * t177 - t149 * t366 - t237 * t242 - t252 * t62, t149 * t252 - t361 (t138 * t235 - t147 * t230) * t133 - t2 * t237 - t199 * t242 + t12 * t177 + t371 * t62 + t366 * t36 + t8 * t252 + (t230 * t274 - t235 * t273) * t149 -(t138 * t230 + t147 * t235) * t133 + t1 * t237 + t199 * t16 - t12 * t178 - t371 * t264 + t367 * t36 - t9 * t252 + (t230 * t273 + t235 * t274) * t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t165 * t393, t165 ^ 2 - t393 ^ 2, t132 - t360, -t133 + t358, t279, -t102 * t165 + t405 * t49 - t289, -t102 * t393 + t405 * t48 - t249, t131 * t290 + t372 (t59 - t365) * t236 + (-t60 - t364) * t231, -t131 * t165 + t158 * t290 + t363, -t158 ^ 2 * t231 + t129 * t165 + t362, -t158 * t165, -pkin(4) * t60 - t49 * t129 - t24 * t165 - t373 + (-pkin(11) * t318 - t113) * t158 + (t48 * t158 + t250) * t231, -pkin(4) * t59 - t49 * t131 + t25 * t165 + t374 + (pkin(11) * t319 + t369) * t158 + t250 * t236, t16 * t194 + t264 * t340, -t16 * t193 + t194 * t242 + t264 * t339 + t340 * t62, t133 * t194 - t149 * t340 + t165 * t264, -t133 * t193 - t149 * t339 + t165 * t62, -t149 * t165 (-t213 * t235 - t214 * t230) * t133 - t222 * t242 + t12 * t193 - t8 * t165 + t277 * t62 + t339 * t36 + (t230 * t269 - t235 * t268) * t149 -(-t213 * t230 + t214 * t235) * t133 + t222 * t16 + t12 * t194 + t9 * t165 - t277 * t264 - t340 * t36 + (t230 * t268 + t235 * t269) * t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131 * t129, -t129 ^ 2 + t131 ^ 2, t59 + t365, t364 - t60, t133, -t131 * t42 + t158 * t25 + t7, t129 * t42 + t158 * t24 - t6, -t398, t392, t390, t389, t133 -(-t19 * t230 - t375) * t149 + (-t131 * t62 + t235 * t133 - t149 * t317) * pkin(5) + t391, t399 + t18 + (-t149 * t20 - t4) * t230 + (t149 * t19 - t296) * t235 + (t131 * t264 - t230 * t133 - t149 * t316) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t398, t392, t390, t389, t133, t149 * t9 + t391, t149 * t8 - t1 + t399;];
tauc_reg  = t3;
