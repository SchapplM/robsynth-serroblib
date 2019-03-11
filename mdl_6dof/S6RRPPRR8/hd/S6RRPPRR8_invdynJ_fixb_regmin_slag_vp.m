% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPPRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% 
% Output:
% tau_reg [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPPRR8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR8_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:26:12
% EndTime: 2019-03-09 09:26:28
% DurationCPUTime: 7.06s
% Computational Cost: add. (5388->585), mult. (12286->760), div. (0->0), fcn. (9133->12), ass. (0->275)
t254 = cos(qJ(2));
t348 = qJD(1) * t254;
t214 = qJD(5) + t348;
t200 = qJD(6) + t214;
t248 = sin(qJ(6));
t252 = cos(qJ(6));
t245 = sin(pkin(10));
t250 = sin(qJ(2));
t349 = qJD(1) * t250;
t326 = t245 * t349;
t246 = cos(pkin(10));
t335 = t246 * qJD(2);
t172 = t326 - t335;
t324 = t246 * t349;
t347 = qJD(2) * t245;
t174 = t324 + t347;
t249 = sin(qJ(5));
t253 = cos(qJ(5));
t288 = t172 * t249 + t174 * t253;
t289 = -t253 * t172 + t174 * t249;
t226 = t246 * qJDD(2);
t334 = qJD(1) * qJD(2);
t317 = t254 * t334;
t333 = t250 * qJDD(1);
t271 = t317 + t333;
t132 = t245 * t271 - t226;
t133 = t245 * qJDD(2) + t246 * t271;
t338 = qJD(5) * t253;
t339 = qJD(5) * t249;
t33 = t249 * t132 + t253 * t133 + t172 * t338 - t174 * t339;
t336 = qJD(6) * t252;
t337 = qJD(6) * t248;
t34 = qJD(5) * t288 - t253 * t132 + t249 * t133;
t274 = t248 * t34 - t252 * t33 + t288 * t337 + t289 * t336;
t40 = t248 * t288 + t252 * t289;
t400 = t200 * t40;
t409 = -t274 + t400;
t292 = t248 * t289 - t252 * t288;
t408 = t292 * t40;
t407 = t200 * t292;
t180 = t245 * t249 + t246 * t253;
t272 = t180 * t254;
t357 = -qJD(1) * t272 - t180 * qJD(5);
t323 = t246 * t348;
t325 = t245 * t348;
t356 = t245 * t338 - t246 * t339 - t249 * t323 + t253 * t325;
t240 = g(3) * t254;
t251 = sin(qJ(1));
t255 = cos(qJ(1));
t297 = g(1) * t255 + g(2) * t251;
t396 = t297 * t250;
t406 = t240 - t396;
t221 = pkin(7) * t333;
t149 = -qJDD(2) * pkin(2) + pkin(7) * t317 + qJDD(3) + t221;
t316 = -t149 - t240;
t259 = -t396 - t316;
t405 = t292 ^ 2 - t40 ^ 2;
t352 = t254 * pkin(2) + t250 * qJ(3);
t401 = -pkin(1) - t352;
t165 = t401 * qJD(1);
t224 = pkin(7) * t348;
t194 = qJD(2) * qJ(3) + t224;
t107 = t165 * t246 - t245 * t194;
t86 = pkin(3) * t348 + qJD(4) - t107;
t52 = pkin(4) * t348 - pkin(8) * t174 + t86;
t108 = t245 * t165 + t246 * t194;
t92 = -qJ(4) * t348 + t108;
t55 = pkin(8) * t172 + t92;
t18 = t249 * t52 + t253 * t55;
t14 = -pkin(9) * t289 + t18;
t12 = t14 * t337;
t243 = qJ(5) + qJ(6);
t230 = sin(t243);
t231 = cos(t243);
t286 = t230 * t245 + t231 * t246;
t186 = -qJD(2) * pkin(2) + pkin(7) * t349 + qJD(3);
t80 = t172 * pkin(3) - t174 * qJ(4) + t186;
t53 = -pkin(4) * t172 - t80;
t32 = pkin(5) * t289 + t53;
t380 = g(3) * t250;
t360 = t251 * t254;
t157 = t245 * t360 + t246 * t255;
t358 = t255 * t245;
t158 = t246 * t360 - t358;
t71 = -t157 * t230 - t158 * t231;
t361 = t251 * t246;
t159 = t254 * t358 - t361;
t359 = t254 * t255;
t160 = t245 * t251 + t246 * t359;
t73 = t159 * t230 + t160 * t231;
t404 = g(1) * t73 - g(2) * t71 + t286 * t380 + t32 * t40 + t12;
t287 = t230 * t246 - t231 * t245;
t229 = t254 * qJDD(1);
t318 = t250 * t334;
t270 = -t318 + t229;
t179 = -qJDD(5) - t270;
t142 = pkin(7) * t270 + qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t296 = pkin(2) * t250 - qJ(3) * t254;
t152 = qJD(2) * t296 - t250 * qJD(3);
t97 = qJD(1) * t152 + qJDD(1) * t401;
t47 = -t245 * t142 + t246 * t97;
t279 = pkin(3) * t229 + qJDD(4) - t47;
t38 = -pkin(3) * t318 + t279;
t22 = pkin(4) * t270 - pkin(8) * t133 + t38;
t48 = t246 * t142 + t245 * t97;
t35 = qJ(4) * t318 + (-qJ(4) * qJDD(1) - qJD(1) * qJD(4)) * t254 + t48;
t25 = pkin(8) * t132 + t35;
t314 = t253 * t22 - t249 * t25;
t263 = -qJD(5) * t18 + t314;
t2 = -pkin(5) * t179 - pkin(9) * t33 + t263;
t276 = -t249 * t22 - t253 * t25 - t52 * t338 + t339 * t55;
t3 = -pkin(9) * t34 - t276;
t329 = t252 * t2 - t248 * t3;
t17 = -t249 * t55 + t253 * t52;
t13 = -pkin(9) * t288 + t17;
t11 = pkin(5) * t214 + t13;
t374 = t14 * t252;
t5 = t11 * t248 + t374;
t70 = t157 * t231 - t158 * t230;
t72 = t159 * t231 - t160 * t230;
t403 = -g(1) * t72 - g(2) * t70 - qJD(6) * t5 + t287 * t380 + t32 * t292 + t329;
t262 = qJD(6) * t292 - t248 * t33 - t252 * t34;
t402 = t262 - t407;
t342 = qJD(3) * t253;
t331 = -pkin(7) * t245 - pkin(3);
t364 = t246 * t254;
t264 = -pkin(8) * t364 + (-pkin(4) + t331) * t250;
t184 = t296 * qJD(1);
t366 = t246 * t184;
t67 = qJD(1) * t264 - t366;
t399 = t245 * t342 - t253 * t67;
t398 = t214 * t288;
t397 = t214 * t289;
t341 = qJD(4) * t245;
t354 = qJ(4) * t323 - t224;
t386 = -pkin(3) - pkin(4);
t308 = -t325 * t386 + t341 - t354;
t170 = t174 ^ 2;
t395 = -t172 ^ 2 - t170;
t379 = -pkin(8) + qJ(3);
t190 = t379 * t245;
t191 = t379 * t246;
t355 = t249 * t190 + t253 * t191;
t346 = qJD(2) * t250;
t394 = qJ(4) * t346 - qJD(4) * t254;
t343 = qJD(3) * t249;
t161 = t245 * t184;
t219 = qJ(4) * t349;
t365 = t246 * t250;
t368 = t245 * t254;
t273 = -pkin(7) * t365 + pkin(8) * t368;
t88 = qJD(1) * t273 + t161 + t219;
t393 = -t190 * t338 + t191 * t339 - t245 * t343 - t246 * t342 + t249 * t67 + t253 * t88;
t392 = t246 * pkin(3) + t245 * qJ(4) + pkin(2);
t163 = t174 * qJD(4);
t391 = t132 * pkin(3) - t133 * qJ(4) - t163;
t388 = t200 ^ 2;
t387 = -0.2e1 * pkin(1);
t385 = pkin(7) * t174;
t384 = g(1) * t251;
t381 = g(2) * t255;
t236 = t254 * pkin(3);
t369 = t245 * t253;
t181 = -t246 * t249 + t369;
t104 = t252 * t180 + t181 * t248;
t378 = -qJD(6) * t104 - t248 * t356 + t252 * t357;
t105 = -t180 * t248 + t181 * t252;
t377 = qJD(6) * t105 + t248 * t357 + t252 * t356;
t137 = pkin(7) * t364 + t245 * t401;
t371 = qJ(4) * t254;
t118 = t137 - t371;
t370 = t245 * t250;
t106 = pkin(8) * t370 + t118;
t208 = pkin(7) * t368;
t95 = pkin(4) * t254 + t208 + t236 + (-pkin(8) * t250 - t401) * t246;
t375 = t253 * t106 + t249 * t95;
t373 = t246 * t48;
t372 = pkin(5) * t356 + t308;
t367 = t246 * t152;
t363 = t250 * t251;
t362 = t250 * t255;
t322 = t254 * t335;
t353 = -qJ(4) * t322 - qJD(4) * t365;
t241 = t250 ^ 2;
t242 = t254 ^ 2;
t351 = t241 - t242;
t350 = qJD(1) * t246;
t345 = qJD(2) * t254;
t344 = qJD(3) * t246;
t332 = pkin(7) * t346;
t330 = pkin(3) * t245 + pkin(7);
t328 = qJ(3) * t346;
t37 = t149 + t391;
t321 = -t37 - t240;
t320 = g(3) * t352;
t319 = qJ(3) * t229;
t315 = qJD(6) * t11 + t3;
t61 = qJD(2) * t264 - t367;
t141 = t245 * t152;
t62 = qJD(2) * t273 + t141 + t394;
t312 = -t249 * t62 + t253 * t61;
t311 = -t106 * t249 + t253 * t95;
t309 = t253 * t190 - t191 * t249;
t136 = t246 * t401 - t208;
t162 = t246 * pkin(4) + t392;
t307 = g(1) * t246 * t362 + g(2) * t250 * t361 + qJD(3) * t325 + t245 * t319;
t306 = t255 * pkin(1) + pkin(2) * t359 + t251 * pkin(7) + qJ(3) * t362;
t305 = t246 * t319;
t304 = t245 * t386 - pkin(7);
t303 = -t246 * qJ(3) * t132 - t172 * t344 - t380;
t76 = -pkin(9) * t180 + t355;
t302 = -pkin(5) * t349 + pkin(9) * t357 + t355 * qJD(5) + qJD(6) * t76 + t246 * t343 - t249 * t88 - t399;
t75 = -pkin(9) * t181 + t309;
t301 = pkin(9) * t356 - qJD(6) * t75 + t393;
t300 = t331 * t250;
t299 = -g(1) * t157 + g(2) * t159;
t298 = g(1) * t158 - g(2) * t160;
t295 = pkin(7) * t172 + t186 * t245;
t146 = t249 * t365 - t250 * t369;
t147 = t180 * t250;
t68 = t252 * t146 + t147 * t248;
t69 = -t146 * t248 + t147 * t252;
t291 = t157 * t253 - t158 * t249;
t290 = t157 * t249 + t158 * t253;
t285 = t248 * t253 + t249 * t252;
t284 = t248 * t249 - t252 * t253;
t283 = qJ(3) * t133 + qJD(3) * t174;
t282 = t214 ^ 2;
t124 = -pkin(7) * t324 + t161;
t115 = -t246 * t332 + t141;
t277 = -pkin(7) * qJDD(2) + t334 * t387;
t275 = -t106 * t339 + t249 * t61 + t253 * t62 + t95 * t338;
t257 = qJD(1) ^ 2;
t268 = pkin(1) * t257 + t297;
t256 = qJD(2) ^ 2;
t267 = pkin(7) * t256 + qJDD(1) * t387 + t381;
t206 = qJ(4) * t365;
t117 = t250 * t304 + t206;
t265 = t401 * t384;
t85 = t304 * t345 - t353;
t26 = -pkin(4) * t132 - t37;
t260 = t245 * t333 - t226 + (-t174 + t347) * t348;
t238 = t255 * pkin(7);
t216 = g(1) * t363;
t212 = qJ(3) * t359;
t209 = qJ(3) * t360;
t167 = -qJDD(6) + t179;
t150 = t172 * t348;
t143 = t250 * t330 - t206;
t130 = pkin(3) * t325 - t354;
t123 = pkin(7) * t326 + t366;
t122 = -t136 + t236;
t114 = t245 * t332 + t367;
t113 = pkin(5) * t180 + t162;
t112 = qJD(1) * t300 - t366;
t111 = t124 + t219;
t110 = t330 * t345 + t353;
t98 = qJD(2) * t300 - t367;
t90 = t159 * t249 + t160 * t253;
t89 = t159 * t253 - t160 * t249;
t84 = t115 + t394;
t81 = t150 + t133;
t79 = qJD(5) * t181 * t250 + qJD(2) * t272;
t78 = qJD(5) * t147 + t249 * t322 - t345 * t369;
t65 = pkin(5) * t146 + t117;
t36 = pkin(5) * t78 + t85;
t28 = -pkin(9) * t146 + t375;
t27 = pkin(5) * t254 - pkin(9) * t147 + t311;
t16 = qJD(6) * t69 + t248 * t79 + t252 * t78;
t15 = -qJD(6) * t68 - t248 * t78 + t252 * t79;
t10 = pkin(5) * t34 + t26;
t9 = -pkin(9) * t78 + t275;
t8 = -pkin(5) * t346 - pkin(9) * t79 - qJD(5) * t375 + t312;
t4 = t11 * t252 - t14 * t248;
t1 = [qJDD(1), -t381 + t384, t297, qJDD(1) * t241 + 0.2e1 * t250 * t317, 0.2e1 * t229 * t250 - 0.2e1 * t334 * t351, qJDD(2) * t250 + t254 * t256, qJDD(2) * t254 - t250 * t256, 0, t277 * t250 + (-t267 + t384) * t254, t250 * t267 + t254 * t277 - t216 (pkin(7) * t132 + t149 * t245 + (qJD(1) * t136 + t107) * qJD(2)) * t250 + (-qJD(1) * t114 + qJD(2) * t295 - qJDD(1) * t136 - t47) * t254 + t298 (pkin(7) * t133 + t149 * t246 + (-qJD(1) * t137 - t108) * qJD(2)) * t250 + (qJD(1) * t115 + qJDD(1) * t137 + t48 + (t186 * t246 + t385) * qJD(2)) * t254 + t299, -t114 * t174 - t115 * t172 - t132 * t137 - t133 * t136 + t216 + (-t107 * t246 - t108 * t245) * t345 + (-t245 * t48 - t246 * t47 - t381) * t250, t48 * t137 + t108 * t115 + t47 * t136 + t107 * t114 - g(1) * t238 - g(2) * t306 - t265 + (t149 * t250 + t186 * t345) * pkin(7), t110 * t172 + t132 * t143 + (t245 * t37 + (-qJD(1) * t122 - t86) * qJD(2)) * t250 + (qJD(1) * t98 + qJDD(1) * t122 + t347 * t80 + t38) * t254 + t298, -t118 * t132 + t122 * t133 - t172 * t84 + t174 * t98 + t216 + (-t245 * t92 + t246 * t86) * t345 + (-t245 * t35 + t246 * t38 - t381) * t250, -t110 * t174 - t133 * t143 + (-t246 * t37 + (qJD(1) * t118 + t92) * qJD(2)) * t250 + (-qJD(1) * t84 - qJDD(1) * t118 - t335 * t80 - t35) * t254 - t299, t35 * t118 + t92 * t84 + t37 * t143 + t80 * t110 + t38 * t122 + t86 * t98 - g(1) * (-pkin(3) * t158 - qJ(4) * t157 + t238) - g(2) * (pkin(3) * t160 + qJ(4) * t159 + t306) - t265, t147 * t33 + t288 * t79, -t146 * t33 - t147 * t34 - t288 * t78 - t289 * t79, -t147 * t179 + t214 * t79 + t254 * t33 - t288 * t346, t146 * t179 - t214 * t78 - t254 * t34 + t289 * t346, -t179 * t254 - t214 * t346, t312 * t214 - t311 * t179 + t314 * t254 - t17 * t346 + t85 * t289 + t117 * t34 + t26 * t146 + t53 * t78 + g(1) * t290 - g(2) * t90 + (-t18 * t254 - t214 * t375) * qJD(5), g(1) * t291 - g(2) * t89 + t117 * t33 + t26 * t147 + t179 * t375 + t18 * t346 - t214 * t275 + t254 * t276 + t288 * t85 + t53 * t79, -t15 * t292 - t274 * t69, -t15 * t40 + t16 * t292 + t262 * t69 + t274 * t68, t15 * t200 - t167 * t69 - t254 * t274 + t292 * t346, -t16 * t200 + t167 * t68 + t254 * t262 + t346 * t40, -t167 * t254 - t200 * t346 (-t248 * t9 + t252 * t8) * t200 - (-t248 * t28 + t252 * t27) * t167 + t329 * t254 - t4 * t346 + t36 * t40 - t65 * t262 + t10 * t68 + t32 * t16 - g(1) * t71 - g(2) * t73 + ((-t248 * t27 - t252 * t28) * t200 - t5 * t254) * qJD(6), t12 * t254 + t5 * t346 - t36 * t292 - t65 * t274 + t10 * t69 + t32 * t15 + g(1) * t70 - g(2) * t72 + (-(-qJD(6) * t28 + t8) * t200 + t27 * t167 - t2 * t254) * t248 + (-(qJD(6) * t27 + t9) * t200 + t28 * t167 - t315 * t254) * t252; 0, 0, 0, -t250 * t257 * t254, t351 * t257, t333, t229, qJDD(2), t250 * t268 - t221 - t240, t380 + (-pkin(7) * qJDD(1) + t268) * t254, -pkin(2) * t132 + t316 * t246 + ((-qJ(3) * t347 - t107) * t250 + (t123 - t295) * t254) * qJD(1) + t307, t305 - pkin(2) * t133 + t259 * t245 + ((-qJ(3) * t335 + t108) * t250 + (-t385 - t124 + (qJD(3) - t186) * t246) * t254) * qJD(1), t123 * t174 + t124 * t172 + t373 + (t107 * t350 - t297) * t254 + (t108 * t348 + t283 - t47) * t245 + t303, -t149 * pkin(2) - t108 * t124 - t107 * t123 - t186 * t224 - g(1) * (-pkin(2) * t362 + t212) - g(2) * (-pkin(2) * t363 + t209) - t320 + (-t107 * t245 + t108 * t246) * qJD(3) + (-t245 * t47 + t373) * qJ(3), -t132 * t392 + t321 * t246 + (-t130 - t341) * t172 + (-t112 * t254 + t250 * t86 + (-t254 * t80 - t328) * t245) * qJD(1) + t307, t111 * t172 - t112 * t174 + t246 * t35 + (-t350 * t86 - t297) * t254 + (t348 * t92 + t283 + t38) * t245 + t303, -t305 + t130 * t174 + t133 * t392 + (t163 + t321 + t396) * t245 + (t111 * t254 - t250 * t92 + (t328 + (-qJD(3) + t80) * t254) * t246) * qJD(1), -t92 * t111 - t80 * t130 - t86 * t112 - g(1) * t212 - g(2) * t209 - t320 + (-g(3) * t236 + qJ(3) * t35 + qJD(3) * t92) * t246 + (-g(3) * t371 + qJ(3) * t38 + qJD(3) * t86 - qJD(4) * t80) * t245 + (-t37 + t396) * t392, t33 * t181 + t288 * t357, -t180 * t33 - t181 * t34 - t288 * t356 - t289 * t357, -t179 * t181 + t214 * t357 + t288 * t349, t180 * t179 - t214 * t356 - t289 * t349, t214 * t349, -t309 * t179 + t162 * t34 + t356 * t53 - g(3) * t272 + (-t191 * t338 + (-qJD(5) * t190 - t344 + t88) * t249 + t399) * t214 + t308 * t289 + t17 * t349 + (t26 + t396) * t180, t355 * t179 + t162 * t33 + t357 * t53 + t393 * t214 + t308 * t288 - t18 * t349 + (t26 - t406) * t181, -t105 * t274 - t292 * t378, t104 * t274 + t105 * t262 + t292 * t377 - t378 * t40, -t105 * t167 + t200 * t378 - t292 * t349, t104 * t167 - t200 * t377 - t349 * t40, t200 * t349 -(-t248 * t76 + t252 * t75) * t167 - t113 * t262 + t10 * t104 + t372 * t40 + t377 * t32 + (t248 * t301 - t252 * t302) * t200 + t4 * t349 - t406 * t286 (t248 * t75 + t252 * t76) * t167 - t113 * t274 + t10 * t105 - t372 * t292 + t378 * t32 + (t248 * t302 + t252 * t301) * t200 - t5 * t349 + t406 * t287; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t260, t81, t395, t107 * t174 + t108 * t172 + t259, t260, t395, -t81, t172 * t92 - t174 * t86 + t259 + t391, 0, 0, 0, 0, 0, -t34 - t398, -t33 + t397, 0, 0, 0, 0, 0, t262 + t407, t274 + t400; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t172 * t174 + t270, -t150 + t133, -t242 * t257 - t170, -g(3) * t370 - g(1) * t159 - g(2) * t157 + t174 * t80 + (-pkin(3) * t346 + t254 * t92) * qJD(1) + t279, 0, 0, 0, 0, 0, -t174 * t289 - t179 * t253 - t249 * t282, -t174 * t288 + t179 * t249 - t253 * t282, 0, 0, 0, 0, 0, t284 * t167 - t174 * t40 - t285 * t388, t285 * t167 + t174 * t292 + t284 * t388; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t288 * t289, t288 ^ 2 - t289 ^ 2, t33 + t397, -t34 + t398, -t179, -g(1) * t89 - g(2) * t291 + g(3) * t146 + t18 * t214 - t288 * t53 + t263, g(1) * t90 + g(2) * t290 + g(3) * t147 + t17 * t214 + t289 * t53 + t276, -t408, t405, t409, t402, -t167 -(-t13 * t248 - t374) * t200 + (-t252 * t167 - t200 * t337 - t288 * t40) * pkin(5) + t403 (-t14 * t200 - t2) * t248 + (t13 * t200 - t315) * t252 + (t248 * t167 - t200 * t336 + t288 * t292) * pkin(5) + t404; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t408, t405, t409, t402, -t167, t200 * t5 + t403, -t248 * t2 + t200 * t4 - t252 * t315 + t404;];
tau_reg  = t1;
