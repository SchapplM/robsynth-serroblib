% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RRRRRR6
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
% Datum: 2019-03-10 04:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRRRR6_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR6_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_inertiaDJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:21:36
% EndTime: 2019-03-10 04:22:13
% DurationCPUTime: 16.32s
% Computational Cost: add. (26479->690), mult. (66654->1181), div. (0->0), fcn. (67574->12), ass. (0->308)
t223 = sin(qJ(3));
t431 = cos(qJ(4));
t368 = t431 * t223;
t433 = -pkin(10) - pkin(9);
t187 = t433 * t368;
t226 = cos(qJ(3));
t429 = sin(qJ(4));
t364 = t429 * t226;
t457 = t433 * t364 + t187;
t221 = sin(pkin(6));
t224 = sin(qJ(2));
t395 = qJD(2) * t224;
t203 = t221 * t395;
t406 = t221 * t224;
t414 = cos(pkin(6));
t292 = t223 * t414 + t226 * t406;
t342 = t414 * t226;
t293 = t223 * t406 - t342;
t134 = t429 * t292 + t431 * t293;
t227 = cos(qJ(2));
t371 = pkin(1) * t414;
t287 = pkin(8) * t406 - t227 * t371;
t164 = -pkin(2) * t414 + t287;
t139 = pkin(3) * t293 + t164;
t257 = t431 * t292 - t429 * t293;
t241 = t134 * pkin(4) - pkin(11) * t257 + t139;
t334 = pkin(3) * t203;
t393 = qJD(3) * t221;
t361 = t224 * t393;
t192 = t223 * t361;
t338 = t414 * qJD(3);
t394 = qJD(2) * t227;
t354 = t221 * t394;
t453 = t338 + t354;
t438 = -t226 * t453 + t192;
t281 = t287 * qJD(2);
t372 = -pkin(2) * t227 - pkin(1);
t305 = -pkin(9) * t224 + t372;
t446 = -t305 * t393 + t281;
t325 = t224 * t371;
t285 = pkin(9) * t414 + t325;
t405 = t221 * t227;
t383 = pkin(8) * t405;
t447 = qJD(3) * (t285 + t383) - qJD(2) * t221 * (pkin(2) * t224 - pkin(9) * t227);
t91 = t223 * t446 - t226 * t447;
t237 = pkin(10) * t438 + t334 + t91;
t373 = t223 * t453 + t226 * t361;
t90 = t223 * t447 + t226 * t446;
t239 = -pkin(10) * t373 - t90;
t265 = pkin(10) * t414 + t285;
t298 = t224 * t433 + t372;
t252 = -t265 * t223 + ((-pkin(8) * t223 - pkin(3)) * t227 + t298 * t226) * t221;
t427 = pkin(8) * t227;
t384 = t226 * t427;
t254 = t265 * t226 + (t223 * t298 + t384) * t221;
t448 = -t431 * t252 + t429 * t254;
t33 = qJD(4) * t448 - t429 * t237 - t431 * t239;
t456 = -pkin(11) * t203 - qJD(5) * t241 + t33;
t72 = t429 * t252 + t431 * t254;
t236 = -pkin(11) * t405 + t72;
t246 = qJD(4) * t134 + t373 * t429 + t438 * t431;
t340 = qJD(2) * t414;
t318 = t224 * t340;
t321 = t373 * pkin(3);
t85 = qJD(4) * t257 + t431 * t373 - t429 * t438;
t455 = -pkin(1) * t318 - t85 * pkin(4) - pkin(8) * t354 - pkin(11) * t246 + qJD(5) * t236 - t321;
t367 = t431 * t226;
t320 = qJD(3) * t367;
t345 = t429 * qJD(4);
t350 = qJD(4) * t431;
t151 = -t320 - t226 * t350 + (qJD(3) * t429 + t345) * t223;
t186 = t368 + t364;
t440 = qJD(3) + qJD(4);
t152 = t440 * t186;
t311 = t433 * t367;
t328 = t433 * t429;
t449 = t223 * t328;
t270 = -t311 + t449;
t392 = qJD(3) * t223;
t379 = pkin(3) * t392;
t454 = -pkin(4) * t152 - pkin(11) * t151 + qJD(5) * t270 - t379;
t184 = t223 * t429 - t367;
t215 = -pkin(3) * t226 - pkin(2);
t291 = pkin(4) * t184 - pkin(11) * t186 + t215;
t436 = -qJD(3) * (t226 * t328 + t187) - t457 * qJD(4);
t452 = -qJD(5) * t291 + t436;
t222 = sin(qJ(5));
t225 = cos(qJ(5));
t199 = t225 * t405;
t256 = qJD(5) * t257;
t245 = qJD(5) * t199 + t225 * t246 + (-t203 + t256) * t222;
t119 = -t222 * t405 + t225 * t257;
t387 = t119 * qJD(5);
t43 = -t222 * t245 + t225 * t387;
t219 = t222 ^ 2;
t220 = t225 ^ 2;
t396 = t219 - t220;
t441 = t396 * qJD(5);
t10 = t222 * t455 + t225 * t456;
t11 = t222 * t456 - t225 * t455;
t451 = -t10 * t225 - t11 * t222;
t402 = t222 * t257 + t199;
t336 = qJD(5) * t402;
t388 = qJD(5) * t222;
t294 = t225 * t395 + t227 * t388;
t415 = t222 * t246 - t225 * t256;
t437 = -t221 * t294 - t415;
t44 = t222 * t336 - t225 * t437;
t380 = t429 * pkin(3);
t212 = t380 + pkin(11);
t333 = pkin(3) * t350;
t312 = t225 * t333;
t445 = t212 * t388 - t312;
t283 = qJD(3) * t293;
t444 = t223 * t283 - t373 * t226;
t430 = cos(qJ(6));
t346 = t430 * qJD(6);
t443 = t430 * qJD(5) + t346;
t357 = t134 * t388;
t442 = -t225 * t85 + t357;
t216 = qJD(5) * t225;
t413 = t151 * t222;
t302 = t186 * t216 - t413;
t412 = t151 * t225;
t301 = t186 * t388 + t412;
t428 = sin(qJ(6));
t363 = t428 * t222;
t365 = t430 * t225;
t183 = t363 - t365;
t439 = qJD(5) + qJD(6);
t435 = t44 + t43;
t434 = t226 ^ 2;
t432 = pkin(12) + pkin(11);
t426 = pkin(9) * t221;
t425 = t225 * pkin(5);
t424 = -pkin(12) - t212;
t366 = t430 * t222;
t185 = t225 * t428 + t366;
t150 = t439 * t185;
t337 = qJD(4) * t72 - t431 * t237 + t429 * t239;
t24 = -t415 * pkin(5) + (-pkin(4) * t395 - pkin(5) * t294) * t221 + t337;
t69 = pkin(4) * t405 + t448;
t53 = pkin(5) * t402 + t69;
t423 = t53 * t150 + t24 * t183;
t349 = qJD(6) * t428;
t149 = (qJD(5) * t428 + t349) * t222 - t443 * t225;
t422 = -t53 * t149 + t24 * t185;
t32 = -pkin(4) * t203 + t337;
t421 = t69 * t216 + t32 * t222;
t40 = t222 * t241 + t225 * t236;
t408 = t186 * t222;
t131 = pkin(5) * t408 - t457;
t114 = -qJD(4) * t311 - t433 * t320 + t440 * t449;
t81 = pkin(5) * t302 + t114;
t420 = t131 * t150 + t81 * t183;
t419 = -t131 * t149 + t81 * t185;
t418 = pkin(3) * qJD(4);
t411 = t457 * t114;
t177 = t325 + t383;
t170 = qJD(2) * t177;
t410 = t170 * t223;
t409 = t186 * t151;
t407 = t186 * t225;
t404 = t225 * t212;
t403 = t114 * t222 - t216 * t457;
t107 = t222 * t291 + t225 * t270;
t332 = pkin(3) * t345;
t375 = pkin(5) * t388;
t188 = t332 + t375;
t382 = t431 * pkin(3);
t213 = -t382 - pkin(4);
t196 = t213 - t425;
t401 = t196 * t150 + t188 * t183;
t400 = -t196 * t149 + t188 * t185;
t214 = -pkin(4) - t425;
t399 = t214 * t150 + t183 * t375;
t398 = -t214 * t149 + t185 * t375;
t397 = t213 * t216 + t222 * t332;
t391 = qJD(3) * t226;
t390 = qJD(3) * t227;
t386 = t164 * qJD(3);
t385 = -0.2e1 * pkin(2) * qJD(3);
t70 = 0.2e1 * t134 * t85;
t126 = 0.2e1 * t184 * t152;
t381 = t430 * pkin(5);
t378 = pkin(4) * t388;
t377 = pkin(4) * t216;
t376 = pkin(9) * t391;
t374 = t222 * t412;
t370 = t222 * t431;
t369 = t225 * t431;
t218 = t221 ^ 2;
t362 = t218 * t394;
t360 = t223 * t390;
t353 = t222 * t216;
t352 = t223 * t391;
t39 = -t222 * t236 + t225 * t241;
t232 = t134 * pkin(5) - t119 * pkin(12) + t39;
t231 = t430 * t232;
t36 = -pkin(12) * t402 + t40;
t18 = -t36 * t428 + t231;
t230 = t428 * t232;
t19 = t36 * t430 + t230;
t228 = t85 * pkin(5) + pkin(12) * t245 + t11;
t233 = -pkin(12) * t437 - t10;
t3 = -qJD(6) * t231 - t428 * t228 - t430 * t233 + t349 * t36;
t4 = -qJD(6) * t230 + t430 * t228 - t233 * t428 - t346 * t36;
t351 = t18 * t149 - t19 * t150 + t3 * t183 - t4 * t185;
t348 = t424 * t222;
t344 = -t32 * t225 + t69 * t388;
t343 = qJD(5) * t424;
t52 = t222 * t452 - t225 * t454;
t235 = t152 * pkin(5) + pkin(12) * t301 + t52;
t106 = -t222 * t270 + t225 * t291;
t253 = t184 * pkin(5) - pkin(12) * t407 + t106;
t250 = t430 * t253;
t51 = t222 * t454 + t225 * t452;
t255 = -pkin(12) * t302 - t51;
t96 = -pkin(12) * t408 + t107;
t14 = -qJD(6) * t250 - t428 * t235 - t430 * t255 + t349 * t96;
t249 = t428 * t253;
t15 = -qJD(6) * t249 + t430 * t235 - t255 * t428 - t346 * t96;
t57 = -t428 * t96 + t250;
t58 = t430 * t96 + t249;
t341 = t14 * t183 + t57 * t149 - t15 * t185 - t58 * t150;
t339 = t151 * t402;
t217 = t225 * pkin(12);
t178 = t217 + t404;
t304 = t430 * t348;
t144 = -t178 * t428 + t304;
t303 = t428 * t348;
t145 = t178 * t430 + t303;
t273 = t222 * t343 + t312;
t313 = t222 * t333;
t274 = t225 * t343 - t313;
t88 = -qJD(6) * t304 + t178 * t349 - t273 * t430 - t274 * t428;
t89 = -qJD(6) * t303 - t178 * t346 - t273 * t428 + t274 * t430;
t335 = t144 * t149 - t145 * t150 + t88 * t183 - t89 * t185;
t331 = pkin(5) * t346;
t330 = pkin(5) * t349;
t197 = t225 * pkin(11) + t217;
t327 = t432 * t430;
t310 = t222 * t327;
t326 = t432 * t428;
t112 = t197 * t349 + t326 * t216 + t310 * t439;
t308 = t222 * t326;
t158 = t197 * t430 - t308;
t113 = -t158 * qJD(6) + (-t225 * t327 + t308) * qJD(5);
t156 = -t197 * t428 - t310;
t329 = t112 * t183 - t113 * t185 + t156 * t149 - t158 * t150;
t324 = t186 * t363;
t181 = t186 ^ 2;
t323 = t181 * t353;
t322 = t224 * t362;
t317 = t430 * t402;
t316 = t428 * t402;
t314 = -t91 * t223 - t90 * t226;
t56 = t134 * t152 + t184 * t85;
t307 = t106 * t225 + t107 * t222;
t306 = t184 * t212 - t186 * t213;
t63 = t134 * t216 + t222 * t85;
t300 = -t152 * t225 + t184 * t388;
t299 = t213 * t388 - t225 * t332;
t296 = (t219 + t220) * t431;
t295 = t223 * t395 - t226 * t390;
t284 = (-t184 * t431 + t186 * t429) * qJD(4);
t275 = -t216 * t39 - t388 * t40 + t451;
t266 = (-t222 * t40 - t225 * t39) * qJD(5) + t451;
t25 = -qJD(5) * t307 - t222 * t52 - t225 * t51;
t261 = pkin(3) * t284 - t151 * t213 - t152 * t212;
t260 = t222 * t437 + t225 * t336;
t242 = t245 * t225;
t202 = -0.2e1 * t353;
t201 = 0.2e1 * t353;
t189 = -0.2e1 * t322;
t182 = -0.2e1 * t441;
t172 = t296 * t418;
t147 = t457 * t388;
t137 = t186 * t365 - t324;
t136 = t185 * t186;
t127 = -0.2e1 * t185 * t149;
t125 = 0.2e1 * t183 * t150;
t124 = t226 * t285 + (t223 * t305 + t384) * t221;
t123 = -t223 * t285 + (-t223 * t427 + t226 * t305) * t221;
t122 = t321 + t170;
t121 = t152 * t222 + t184 * t216;
t103 = t186 * t441 + t374;
t95 = -t150 * t184 - t152 * t183;
t94 = -t149 * t184 + t152 * t185;
t93 = t151 * t396 - 0.4e1 * t186 * t353;
t92 = 0.2e1 * t149 * t183 - 0.2e1 * t150 * t185;
t86 = (t430 * t149 - t428 * t150 + (-t183 * t430 + t185 * t428) * qJD(6)) * pkin(5);
t80 = t119 * t430 - t316;
t79 = t119 * t428 + t317;
t68 = -t151 * t366 + (-t428 * t151 + t186 * t443) * t225 - t439 * t324;
t67 = t150 * t186 - t151 * t183;
t55 = -t134 * t150 - t183 * t85;
t54 = -t134 * t149 + t185 * t85;
t46 = t136 * t150 + t183 * t68;
t45 = -t137 * t149 - t185 * t67;
t29 = t136 * t149 - t137 * t150 + t183 * t67 - t185 * t68;
t28 = -t222 * t387 - t242 - t260;
t27 = -qJD(6) * t316 + t119 * t346 - t245 * t428 + t430 * t437;
t26 = qJD(6) * t317 + t119 * t349 + t245 * t430 + t428 * t437;
t21 = t150 * t79 + t183 * t27;
t20 = -t149 * t80 - t185 * t26;
t5 = t149 * t79 - t150 * t80 + t183 * t26 - t185 * t27;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t322, 0.2e1 * (-t224 ^ 2 + t227 ^ 2) * t218 * qJD(2), 0.2e1 * t340 * t405, t189, -0.2e1 * t221 * t318, 0, -0.2e1 * pkin(1) * t218 * t395 - 0.2e1 * t170 * t414, -0.2e1 * pkin(1) * t362 + 0.2e1 * t281 * t414, 0.2e1 * t170 * t406 - 0.2e1 * t177 * t203 - 0.2e1 * t281 * t405 + 0.2e1 * t287 * t354, 0.2e1 * t170 * t287 - 0.2e1 * t177 * t281, -0.2e1 * t292 * t438, -0.2e1 * t292 * t373 + 0.2e1 * t293 * t438, 0.2e1 * t203 * t292 + 0.2e1 * t405 * t438, 0.2e1 * t293 * t373, -0.2e1 * t203 * t293 + 0.2e1 * t373 * t405, t189, -0.2e1 * t170 * t342 + 0.2e1 * t164 * t373 + 0.2e1 * (-t91 * t227 + (qJD(2) * t123 + t410) * t224) * t221, -0.2e1 * t124 * t203 - 0.2e1 * t164 * t438 + 0.2e1 * t170 * t292 - 0.2e1 * t405 * t90, 0.2e1 * t123 * t438 - 0.2e1 * t124 * t373 - 0.2e1 * t292 * t91 + 0.2e1 * t293 * t90, 0.2e1 * t123 * t91 - 0.2e1 * t124 * t90 + 0.2e1 * t164 * t170, -0.2e1 * t257 * t246, 0.2e1 * t134 * t246 - 0.2e1 * t257 * t85, 0.2e1 * t203 * t257 + 0.2e1 * t246 * t405, t70, 0.2e1 * (-t134 * t395 + t227 * t85) * t221, t189, 0.2e1 * t122 * t134 + 0.2e1 * t139 * t85 + 0.2e1 * (t227 * t337 - t395 * t448) * t221, 0.2e1 * t122 * t257 - 0.2e1 * t139 * t246 - 0.2e1 * t203 * t72 - 0.2e1 * t33 * t405, 0.2e1 * t33 * t134 - 0.2e1 * t246 * t448 + 0.2e1 * t257 * t337 - 0.2e1 * t72 * t85, 0.2e1 * t122 * t139 - 0.2e1 * t33 * t72 + 0.2e1 * t337 * t448, -0.2e1 * t119 * t245, -0.2e1 * t119 * t437 + 0.2e1 * t245 * t402, 0.2e1 * t119 * t85 - 0.2e1 * t134 * t245, 0.2e1 * t402 * t437, -0.2e1 * t134 * t437 - 0.2e1 * t402 * t85, t70, 0.2e1 * t11 * t134 + 0.2e1 * t32 * t402 + 0.2e1 * t39 * t85 + 0.2e1 * t437 * t69, 0.2e1 * t10 * t134 + 0.2e1 * t32 * t119 - 0.2e1 * t245 * t69 - 0.2e1 * t40 * t85, 0.2e1 * t10 * t402 - 0.2e1 * t11 * t119 + 0.2e1 * t245 * t39 - 0.2e1 * t40 * t437, -0.2e1 * t10 * t40 + 0.2e1 * t11 * t39 + 0.2e1 * t32 * t69, -0.2e1 * t80 * t26, 0.2e1 * t26 * t79 - 0.2e1 * t27 * t80, -0.2e1 * t134 * t26 + 0.2e1 * t80 * t85, 0.2e1 * t79 * t27, -0.2e1 * t134 * t27 - 0.2e1 * t79 * t85, t70, 0.2e1 * t134 * t4 + 0.2e1 * t18 * t85 + 0.2e1 * t24 * t79 + 0.2e1 * t27 * t53, 0.2e1 * t134 * t3 - 0.2e1 * t19 * t85 + 0.2e1 * t24 * t80 - 0.2e1 * t26 * t53, 0.2e1 * t18 * t26 - 0.2e1 * t19 * t27 + 0.2e1 * t3 * t79 - 0.2e1 * t4 * t80, 0.2e1 * t18 * t4 - 0.2e1 * t19 * t3 + 0.2e1 * t24 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t354, 0, -t203, 0, -t170, t281, 0, 0, t434 * t361 + (-t192 + (0.2e1 * t338 + t354) * t226) * t223, -t223 * t373 - t292 * t392 + (-t438 - t283) * t226, t295 * t221, t444 (t226 * t395 + t360) * t221, 0, -pkin(2) * t373 - t170 * t226 + t223 * t386 - t295 * t426, pkin(2) * t438 - t360 * t426 + t410 + (-pkin(9) * t203 + t386) * t226, -t123 * t391 - t124 * t392 + t292 * t376 + t314 + (-t223 * t438 + t444) * pkin(9), -pkin(2) * t170 + ((-t123 * t226 - t124 * t223) * qJD(3) + t314) * pkin(9), -t151 * t257 - t186 * t246, t151 * t134 - t152 * t257 + t184 * t246 - t186 * t85 (t151 * t227 + t186 * t395) * t221, t56 (t152 * t227 - t184 * t395) * t221, 0, t134 * t379 + t122 * t184 + t139 * t152 + t215 * t85 + (t114 * t227 + t395 * t457) * t221, t122 * t186 - t139 * t151 - t203 * t270 - t215 * t246 + t257 * t379 - t405 * t436, t114 * t257 + t134 * t436 - t151 * t448 - t72 * t152 + t33 * t184 + t186 * t337 + t246 * t457 - t270 * t85, t114 * t448 + t122 * t215 + t139 * t379 - t270 * t33 - t337 * t457 - t436 * t72, -t119 * t301 - t186 * t242, t119 * t413 + t225 * t339 + (t44 - t43) * t186, t119 * t152 - t134 * t412 - t184 * t245 - t186 * t357 + t407 * t85, t186 * t260 - t222 * t339, -t134 * t302 - t152 * t402 - t184 * t437 - t408 * t85, t56, t106 * t85 + t11 * t184 + t114 * t402 + t52 * t134 + t39 * t152 + t186 * t421 - t413 * t69 - t437 * t457, t10 * t184 - t107 * t85 + t114 * t119 + t51 * t134 - t40 * t152 + t245 * t457 - t301 * t69 + t32 * t407, t10 * t408 + t106 * t245 - t107 * t437 - t11 * t407 - t52 * t119 + t301 * t39 - t302 * t40 + t402 * t51, -t10 * t107 + t106 * t11 + t114 * t69 - t32 * t457 + t39 * t52 - t40 * t51, -t137 * t26 - t67 * t80, t136 * t26 - t137 * t27 + t67 * t79 - t68 * t80, -t134 * t67 + t137 * t85 + t152 * t80 - t184 * t26, t136 * t27 + t68 * t79, -t134 * t68 - t136 * t85 - t152 * t79 - t184 * t27, t56, t131 * t27 + t134 * t15 + t136 * t24 + t152 * t18 + t184 * t4 + t53 * t68 + t57 * t85 + t79 * t81, -t131 * t26 + t134 * t14 + t137 * t24 - t152 * t19 + t184 * t3 - t53 * t67 - t58 * t85 + t80 * t81, t136 * t3 - t137 * t4 + t14 * t79 - t15 * t80 + t18 * t67 - t19 * t68 + t26 * t57 - t27 * t58, t131 * t24 - t14 * t19 + t15 * t18 - t3 * t58 + t4 * t57 + t53 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t352, 0.2e1 * (-t223 ^ 2 + t434) * qJD(3), 0, -0.2e1 * t352, 0, 0, t223 * t385, t226 * t385, 0, 0, -0.2e1 * t409, 0.2e1 * t151 * t184 - 0.2e1 * t152 * t186, 0, t126, 0, 0, 0.2e1 * t152 * t215 + 0.2e1 * t184 * t379, -0.2e1 * t151 * t215 + 0.2e1 * t186 * t379, 0.2e1 * t114 * t186 + 0.2e1 * t151 * t457 - 0.2e1 * t152 * t270 + 0.2e1 * t184 * t436, 0.2e1 * t215 * t379 - 0.2e1 * t270 * t436 - 0.2e1 * t411, -0.2e1 * t220 * t409 - 0.2e1 * t323, 0.2e1 * t181 * t441 + 0.4e1 * t186 * t374, 0.2e1 * t152 * t407 - 0.2e1 * t184 * t301, -0.2e1 * t219 * t409 + 0.2e1 * t323, -0.2e1 * t152 * t408 - 0.2e1 * t184 * t302, t126, 0.2e1 * t106 * t152 + 0.2e1 * t114 * t408 + 0.2e1 * t184 * t52 - 0.2e1 * t302 * t457, -0.2e1 * t107 * t152 + 0.2e1 * t114 * t407 + 0.2e1 * t184 * t51 + 0.2e1 * t301 * t457, 0.2e1 * t307 * t151 + 0.2e1 * (t222 * t51 - t225 * t52 + (t106 * t222 - t107 * t225) * qJD(5)) * t186, 0.2e1 * t106 * t52 - 0.2e1 * t107 * t51 - 0.2e1 * t411, -0.2e1 * t137 * t67, 0.2e1 * t136 * t67 - 0.2e1 * t137 * t68, 0.2e1 * t137 * t152 - 0.2e1 * t184 * t67, 0.2e1 * t136 * t68, -0.2e1 * t136 * t152 - 0.2e1 * t184 * t68, t126, 0.2e1 * t131 * t68 + 0.2e1 * t136 * t81 + 0.2e1 * t15 * t184 + 0.2e1 * t152 * t57, -0.2e1 * t131 * t67 + 0.2e1 * t137 * t81 + 0.2e1 * t14 * t184 - 0.2e1 * t152 * t58, 0.2e1 * t136 * t14 - 0.2e1 * t137 * t15 + 0.2e1 * t57 * t67 - 0.2e1 * t58 * t68, 0.2e1 * t131 * t81 - 0.2e1 * t14 * t58 + 0.2e1 * t15 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t438, 0, -t373, t203, t91, t90, 0, 0, 0, 0, -t246, 0, -t85, t203 (t227 * t345 + t395 * t431) * t221 * pkin(3) - t337, t333 * t405 - t334 * t429 + t33, -t134 * t333 + t246 * t382 + t257 * t332 - t380 * t85 (-t429 * t33 - t431 * t337 + (t429 * t448 + t431 * t72) * qJD(4)) * pkin(3), t43, t28, t63, t44, -t442, 0, t213 * t437 - t63 * t212 + (-t134 * t370 + t402 * t429) * t418 + t344, t119 * t332 + t134 * t445 - t213 * t245 - t404 * t85 + t421, t119 * t313 + t212 * t435 - t312 * t402 + t275, t32 * t213 + (t369 * t40 - t370 * t39 + t429 * t69) * t418 + t266 * t212, t20, t5, t54, t21, t55, 0, t134 * t89 + t144 * t85 + t188 * t79 + t196 * t27 + t423, t134 * t88 - t145 * t85 + t188 * t80 - t196 * t26 + t422, t144 * t26 - t145 * t27 + t79 * t88 - t80 * t89 + t351, t144 * t4 - t145 * t3 + t18 * t89 + t188 * t53 - t19 * t88 + t196 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t391, 0, -t392, 0, -t376, pkin(9) * t392, 0, 0, 0, 0, -t151, 0, -t152, 0, -t114, t436 (t151 * t431 - t152 * t429 + t284) * pkin(3), -t114 * t382 + t270 * t333 - t332 * t457 - t380 * t436, -t103, t93, t121, t103, -t300, 0, -t147 + (-qJD(5) * t306 - t114) * t225 + t261 * t222, t225 * t261 + t306 * t388 + t403, t25, t114 * t213 + (-t106 * t370 + t107 * t369 - t429 * t457) * t418 + t25 * t212, t45, t29, t94, t46, t95, 0, t136 * t188 + t144 * t152 + t184 * t89 + t196 * t68 + t420, t137 * t188 - t145 * t152 + t184 * t88 - t196 * t67 + t419, t136 * t88 - t137 * t89 + t144 * t67 - t145 * t68 + t341, t131 * t188 - t14 * t145 + t144 * t15 + t196 * t81 + t57 * t89 - t58 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t332, -0.2e1 * t333, 0, 0, t201, t182, 0, t202, 0, 0, 0.2e1 * t299, 0.2e1 * t397, 0.2e1 * t172, 0.2e1 * (t212 * t296 + t213 * t429) * t418, t127, t92, 0, t125, 0, 0, 0.2e1 * t401, 0.2e1 * t400, 0.2e1 * t335, 0.2e1 * t144 * t89 - 0.2e1 * t145 * t88 + 0.2e1 * t188 * t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t246, 0, -t85, t203, -t337, t33, 0, 0, t43, t28, t63, t44, -t442, 0, -pkin(4) * t437 - pkin(11) * t63 + t344, pkin(4) * t245 + pkin(11) * t442 + t421, pkin(11) * t435 + t275, -pkin(4) * t32 + pkin(11) * t266, t20, t5, t54, t21, t55, 0, t113 * t134 + t156 * t85 + t214 * t27 + t375 * t79 + t423, t112 * t134 - t158 * t85 - t214 * t26 + t375 * t80 + t422, t112 * t79 - t113 * t80 + t156 * t26 - t158 * t27 + t351, -t112 * t19 + t113 * t18 + t156 * t4 - t158 * t3 + t214 * t24 + t375 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t151, 0, -t152, 0, -t114, t436, 0, 0, -t103, t93, t121, t103, -t300, 0, -t147 + (pkin(4) * t151 - pkin(11) * t152) * t222 + (-t114 + (-pkin(4) * t186 - pkin(11) * t184) * qJD(5)) * t225, pkin(4) * t301 + pkin(11) * t300 + t403, t25, -pkin(4) * t114 + pkin(11) * t25, t45, t29, t94, t46, t95, 0, t113 * t184 + t136 * t375 + t152 * t156 + t214 * t68 + t420, t112 * t184 + t137 * t375 - t152 * t158 - t214 * t67 + t419, t112 * t136 - t113 * t137 + t156 * t67 - t158 * t68 + t341, -t112 * t58 + t113 * t57 + t131 * t375 - t14 * t158 + t15 * t156 + t214 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t332, -t333, 0, 0, t201, t182, 0, t202, 0, 0, t299 - t378, -t377 + t397, t172 (-pkin(4) * t429 + pkin(11) * t296) * t418, t127, t92, 0, t125, 0, 0, t399 + t401, t398 + t400, t329 + t335, -t112 * t145 + t113 * t144 + t156 * t89 - t158 * t88 + t188 * t214 + t196 * t375; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t201, t182, 0, t202, 0, 0, -0.2e1 * t378, -0.2e1 * t377, 0, 0, t127, t92, 0, t125, 0, 0, 0.2e1 * t399, 0.2e1 * t398, 0.2e1 * t329, -0.2e1 * t112 * t158 + 0.2e1 * t113 * t156 + 0.2e1 * t214 * t375; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t245, 0, -t437, t85, t11, t10, 0, 0, 0, 0, -t26, 0, -t27, t85, -t134 * t330 + t381 * t85 + t4 (-t134 * t346 - t428 * t85) * pkin(5) + t3 (t430 * t26 - t428 * t27 + (t428 * t80 - t430 * t79) * qJD(6)) * pkin(5) (t430 * t4 - t428 * t3 + (-t18 * t428 + t19 * t430) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t301, 0, -t302, t152, t52, t51, 0, 0, 0, 0, -t67, 0, -t68, t152, t152 * t381 - t184 * t330 + t15 (-t152 * t428 - t184 * t346) * pkin(5) + t14 (t430 * t67 - t428 * t68 + (-t136 * t430 + t137 * t428) * qJD(6)) * pkin(5) (t430 * t15 - t428 * t14 + (-t428 * t57 + t430 * t58) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t216, 0, -t388, 0, -t212 * t216 - t313, t445, 0, 0, 0, 0, -t149, 0, -t150, 0, t89, t88, t86 (t430 * t89 - t428 * t88 + (-t144 * t428 + t145 * t430) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t216, 0, -t388, 0, -pkin(11) * t216, pkin(11) * t388, 0, 0, 0, 0, -t149, 0, -t150, 0, t113, t112, t86 (t430 * t113 - t428 * t112 + (-t156 * t428 + t158 * t430) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t330, -0.2e1 * t331, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, 0, -t27, t85, t4, t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, 0, -t68, t152, t15, t14, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t149, 0, -t150, 0, t89, t88, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t149, 0, -t150, 0, t113, t112, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t330, -t331, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t1;
