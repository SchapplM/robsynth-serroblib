% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6PRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRPRR5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:20:47
% EndTime: 2019-03-08 22:21:07
% DurationCPUTime: 8.36s
% Computational Cost: add. (11157->527), mult. (27980->758), div. (0->0), fcn. (21757->12), ass. (0->244)
t250 = sin(qJ(3));
t253 = cos(qJ(3));
t271 = pkin(3) * t250 - qJ(4) * t253;
t181 = qJD(3) * t271 - t250 * qJD(4);
t244 = sin(pkin(12));
t305 = qJD(3) * t250;
t298 = pkin(8) * t305;
t224 = t244 * t298;
t246 = cos(pkin(12));
t251 = sin(qJ(2));
t245 = sin(pkin(6));
t311 = qJD(1) * t245;
t254 = cos(qJ(2));
t322 = t253 * t254;
t319 = t246 * t181 - (-t244 * t322 + t246 * t251) * t311 + t224;
t378 = t244 * t181 - (t244 * t251 + t246 * t322) * t311;
t324 = t246 * t253;
t268 = pkin(4) * t250 - pkin(9) * t324;
t264 = t268 * qJD(3);
t377 = t264 + t319;
t325 = t246 * t250;
t330 = t244 * t253;
t376 = (-pkin(8) * t325 - pkin(9) * t330) * qJD(3) + t378;
t300 = t246 * qJD(3);
t308 = qJD(2) * t250;
t196 = t244 * t308 - t300;
t292 = t246 * t308;
t306 = qJD(3) * t244;
t198 = t292 + t306;
t249 = sin(qJ(5));
t252 = cos(qJ(5));
t128 = t196 * t249 - t198 * t252;
t129 = t252 * t196 + t198 * t249;
t248 = sin(qJ(6));
t351 = cos(qJ(6));
t63 = t351 * t128 + t248 * t129;
t375 = t63 ^ 2;
t67 = t248 * t128 - t351 * t129;
t374 = t67 ^ 2;
t307 = qJD(2) * t253;
t232 = -qJD(5) + t307;
t226 = -qJD(6) + t232;
t373 = t226 * t67;
t372 = t63 * t226;
t326 = t246 * t249;
t204 = t244 * t252 + t326;
t185 = t204 * qJD(5);
t265 = t204 * t253;
t317 = qJD(2) * t265 - t185;
t323 = t252 * t246;
t203 = t244 * t249 - t323;
t302 = qJD(5) * t252;
t303 = qJD(5) * t249;
t357 = -t244 * t303 + t246 * t302;
t363 = t203 * t307 + t357;
t214 = -pkin(3) * t253 - qJ(4) * t250 - pkin(2);
t195 = t246 * t214;
t137 = -pkin(9) * t325 + t195 + (-pkin(8) * t244 - pkin(4)) * t253;
t165 = pkin(8) * t324 + t244 * t214;
t331 = t244 * t250;
t149 = -pkin(9) * t331 + t165;
t346 = t137 * t302 - t149 * t303 + t249 * t377 + t376 * t252;
t73 = t249 * t137 + t252 * t149;
t345 = -t73 * qJD(5) - t376 * t249 + t252 * t377;
t349 = pkin(9) + qJ(4);
t217 = t349 * t244;
t218 = t349 * t246;
t288 = t251 * t311;
t212 = qJD(2) * pkin(8) + t288;
t193 = t250 * t212;
t247 = cos(pkin(6));
t310 = qJD(1) * t247;
t161 = t253 * t310 - t193;
t205 = t271 * qJD(2);
t110 = -t244 * t161 + t246 * t205;
t86 = qJD(2) * t268 + t110;
t111 = t246 * t161 + t244 * t205;
t295 = t244 * t307;
t93 = -pkin(9) * t295 + t111;
t344 = -qJD(4) * t323 + t217 * t302 + t252 * t93 + (qJD(4) * t244 + qJD(5) * t218 + t86) * t249;
t151 = -t249 * t217 + t252 * t218;
t339 = -t204 * qJD(4) - t151 * qJD(5) + t249 * t93 - t252 * t86;
t371 = t128 ^ 2;
t370 = t129 ^ 2;
t304 = qJD(3) * t253;
t291 = t244 * t304;
t115 = -t252 * t253 * t300 + t185 * t250 + t249 * t291;
t369 = pkin(5) * t305 + t115 * pkin(10) + t345;
t116 = qJD(3) * t265 + t357 * t250;
t368 = -pkin(10) * t116 + t346;
t367 = pkin(5) * t308 + pkin(10) * t363 - t339;
t366 = -pkin(10) * t317 + t344;
t350 = t67 * t63;
t365 = t128 * t232;
t364 = t129 * t232;
t287 = t254 * t311;
t275 = qJD(2) * t287;
t362 = qJD(3) * t310 + t275;
t361 = -t374 + t375;
t285 = qJD(6) * t351;
t301 = qJD(6) * t248;
t299 = qJD(2) * qJD(3);
t283 = t253 * t299;
t273 = t252 * t283;
t274 = t244 * t283;
t83 = t249 * (qJD(5) * t198 + t274) + t196 * t302 - t246 * t273;
t84 = -t196 * t303 + t198 * t302 + t244 * t273 + t283 * t326;
t21 = -t128 * t301 + t129 * t285 + t248 * t84 + t351 * t83;
t360 = -t21 + t373;
t162 = t212 * t253 + t250 * t310;
t156 = qJD(3) * qJ(4) + t162;
t163 = qJD(2) * t214 - t287;
t87 = -t156 * t244 + t246 * t163;
t61 = -pkin(4) * t307 - pkin(9) * t198 + t87;
t88 = t246 * t156 + t244 * t163;
t71 = -pkin(9) * t196 + t88;
t34 = -t249 * t71 + t252 * t61;
t26 = pkin(10) * t128 + t34;
t24 = -pkin(5) * t232 + t26;
t35 = t249 * t61 + t252 * t71;
t27 = -pkin(10) * t129 + t35;
t314 = t362 * t253;
t118 = (qJD(4) - t193) * qJD(3) + t314;
t141 = (t181 + t288) * qJD(2);
t58 = -t244 * t118 + t246 * t141;
t49 = qJD(2) * t264 + t58;
t59 = t246 * t118 + t244 * t141;
t55 = -pkin(9) * t274 + t59;
t15 = -qJD(5) * t35 - t249 * t55 + t252 * t49;
t236 = t250 * t299;
t8 = pkin(5) * t236 + t83 * pkin(10) + t15;
t14 = t249 * t49 + t252 * t55 + t61 * t302 - t303 * t71;
t9 = -pkin(10) * t84 + t14;
t260 = -t24 * t285 - t248 * t8 + t27 * t301 - t351 * t9;
t281 = -qJD(3) * pkin(3) + qJD(4);
t152 = -t161 + t281;
t112 = pkin(4) * t196 + t152;
t60 = pkin(5) * t129 + t112;
t359 = -t60 * t67 + t260;
t242 = t250 ^ 2;
t243 = t253 ^ 2;
t356 = qJD(2) * (t242 - 0.2e1 * t243);
t297 = t351 * t27;
t11 = t248 * t24 + t297;
t2 = -qJD(6) * t11 - t248 * t9 + t351 * t8;
t355 = t60 * t63 + t2;
t22 = -qJD(6) * t63 - t248 * t83 + t351 * t84;
t354 = -t22 + t372;
t176 = t203 * t250;
t72 = t252 * t137 - t149 * t249;
t50 = -pkin(5) * t253 + pkin(10) * t176 + t72;
t175 = t204 * t250;
t54 = -pkin(10) * t175 + t73;
t28 = -t248 * t54 + t351 * t50;
t353 = t28 * qJD(6) + t369 * t248 + t368 * t351;
t29 = t248 * t50 + t351 * t54;
t352 = -t29 * qJD(6) - t368 * t248 + t369 * t351;
t150 = -t252 * t217 - t218 * t249;
t113 = -pkin(10) * t204 + t150;
t114 = -pkin(10) * t203 + t151;
t53 = t248 * t113 + t351 * t114;
t348 = t53 * qJD(6) - t366 * t248 + t367 * t351;
t52 = t351 * t113 - t248 * t114;
t347 = -t52 * qJD(6) + t367 * t248 + t366 * t351;
t343 = t203 * t285 + t204 * t301 - t317 * t248 - t363 * t351;
t134 = -t248 * t203 + t351 * t204;
t342 = t134 * qJD(6) + t363 * t248 - t317 * t351;
t341 = qJD(2) * pkin(2);
t340 = t248 * t27;
t123 = t212 * t304 + t362 * t250;
t329 = t245 * t251;
t187 = -t247 * t253 + t250 * t329;
t338 = t123 * t187;
t337 = t123 * t244;
t336 = t123 * t246;
t335 = t123 * t250;
t334 = t128 * t129;
t333 = t196 * t246;
t256 = qJD(2) ^ 2;
t332 = t243 * t256;
t328 = t245 * t254;
t327 = t245 * t256;
t255 = qJD(3) ^ 2;
t321 = t255 * t250;
t320 = t255 * t253;
t279 = t246 * t298;
t318 = -t279 + t378;
t191 = pkin(4) * t291 + pkin(8) * t304;
t206 = pkin(4) * t331 + t250 * pkin(8);
t312 = t242 - t243;
t309 = qJD(2) * t245;
t296 = t251 * t327;
t235 = -pkin(4) * t246 - pkin(3);
t294 = t251 * t309;
t293 = t254 * t309;
t140 = pkin(4) * t295 + t162;
t286 = -t317 * pkin(5) - t140;
t280 = t196 + t300;
t96 = pkin(4) * t274 + t123;
t278 = t250 * t287;
t277 = t250 * t293;
t276 = t253 * t293;
t272 = t253 * t236;
t188 = t247 * t250 + t253 * t329;
t145 = -t188 * t244 - t246 * t328;
t146 = t188 * t246 - t244 * t328;
t74 = t145 * t252 - t146 * t249;
t75 = t145 * t249 + t146 * t252;
t270 = qJD(2) * t280;
t269 = qJD(2) * (-t198 + t306);
t37 = -t248 * t75 + t351 * t74;
t38 = t248 * t74 + t351 * t75;
t46 = pkin(5) * t84 + t96;
t106 = -t248 * t175 - t351 * t176;
t266 = t253 * t269;
t213 = -t287 - t341;
t262 = qJD(3) * (t213 + t287 - t341);
t261 = -qJ(4) * t305 + (-t152 + t281) * t253;
t122 = -t212 * t305 + t314;
t257 = t122 * t253 + t335 + (-t161 * t253 - t162 * t250) * qJD(3);
t241 = t246 ^ 2;
t240 = t244 ^ 2;
t229 = t250 * t256 * t253;
t220 = -0.2e1 * t272;
t169 = pkin(5) * t203 + t235;
t164 = -pkin(8) * t330 + t195;
t148 = qJD(3) * t188 + t277;
t147 = -qJD(3) * t187 + t276;
t139 = pkin(5) * t175 + t206;
t133 = t351 * t203 + t204 * t248;
t109 = t147 * t246 + t244 * t294;
t108 = -t147 * t244 + t246 * t294;
t105 = t351 * t175 - t176 * t248;
t85 = pkin(5) * t116 + t191;
t41 = t106 * qJD(6) - t248 * t115 + t351 * t116;
t40 = t351 * t115 + t248 * t116 + t175 * t285 - t176 * t301;
t31 = -qJD(5) * t75 + t252 * t108 - t249 * t109;
t30 = qJD(5) * t74 + t249 * t108 + t252 * t109;
t13 = t351 * t26 - t340;
t12 = -t248 * t26 - t297;
t10 = t351 * t24 - t340;
t6 = -t38 * qJD(6) - t248 * t30 + t351 * t31;
t5 = t37 * qJD(6) + t248 * t31 + t351 * t30;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t296, -t254 * t327, 0, 0, 0, 0, 0, 0, 0, 0, -t253 * t296 + (-t148 - t277) * qJD(3), t250 * t296 + (-t147 - t276) * qJD(3) (t147 * t253 + t148 * t250 + (t187 * t253 - t188 * t250) * qJD(3)) * qJD(2), t122 * t188 + t338 + t162 * t147 - t161 * t148 + (t213 - t287) * t294, 0, 0, 0, 0, 0, 0, t148 * t196 + (-t108 * t253 + (t145 * t250 + t187 * t330) * qJD(3)) * qJD(2), t148 * t198 + (t109 * t253 + (-t146 * t250 + t187 * t324) * qJD(3)) * qJD(2), -t108 * t198 - t109 * t196 + (-t145 * t246 - t146 * t244) * t283, t108 * t87 + t109 * t88 + t145 * t58 + t146 * t59 + t148 * t152 + t338, 0, 0, 0, 0, 0, 0, t129 * t148 + t187 * t84 - t232 * t31 + t236 * t74, -t128 * t148 - t187 * t83 + t232 * t30 - t236 * t75, t128 * t31 - t129 * t30 + t74 * t83 - t75 * t84, t112 * t148 + t14 * t75 + t15 * t74 + t187 * t96 + t30 * t35 + t31 * t34, 0, 0, 0, 0, 0, 0, -t148 * t67 + t187 * t22 - t226 * t6 + t236 * t37, -t148 * t63 - t187 * t21 + t226 * t5 - t236 * t38, t21 * t37 - t22 * t38 + t5 * t67 + t6 * t63, t10 * t6 + t11 * t5 + t148 * t60 + t187 * t46 + t2 * t37 - t260 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t272, -0.2e1 * t312 * t299, t320, t220, -t321, 0, -pkin(8) * t320 + t250 * t262, pkin(8) * t321 + t253 * t262 (-t242 - t243) * t275 + t257 ((t161 * t250 - t162 * t253) * t254 + (-t213 - t341) * t251) * t311 + t257 * pkin(8) (t198 * t246 + t241 * t308) * t304 (-t333 + (-t198 - 0.2e1 * t292) * t244) * t304 (t198 * t250 + t246 * t356) * qJD(3) (t196 * t244 + t240 * t308) * t304 (-t196 * t250 - t244 * t356) * qJD(3), t220 (-t196 * t287 + t337 + (qJD(2) * t164 + t87) * qJD(3)) * t250 + (-t58 + (pkin(8) * t196 + t152 * t244) * qJD(3) + (t224 - t319) * qJD(2)) * t253 (-t198 * t287 + t336 + (-qJD(2) * t165 - t88) * qJD(3)) * t250 + (t59 + (pkin(8) * t198 + t152 * t246) * qJD(3) + (t279 + t318) * qJD(2)) * t253 (-t244 * t59 - t246 * t58) * t250 - t319 * t198 - t318 * t196 + (-t244 * t88 - t246 * t87 + (-t164 * t246 - t165 * t244) * qJD(2)) * t304, -t152 * t278 + t58 * t164 + t59 * t165 + t318 * t88 + t319 * t87 + (t152 * t304 + t335) * pkin(8), t115 * t128 + t176 * t83, t115 * t129 + t116 * t128 + t175 * t83 + t176 * t84, t115 * t232 + t83 * t253 + (-qJD(2) * t176 - t128) * t305, t116 * t129 + t175 * t84, t116 * t232 + t84 * t253 + (-qJD(2) * t175 - t129) * t305 (-t232 - t307) * t305, t112 * t116 + t191 * t129 - t15 * t253 + t96 * t175 + t206 * t84 - t345 * t232 + (-t129 * t287 + (qJD(2) * t72 + t34) * qJD(3)) * t250, -t112 * t115 - t191 * t128 + t14 * t253 - t96 * t176 - t206 * t83 + t346 * t232 + (t128 * t287 + (-qJD(2) * t73 - t35) * qJD(3)) * t250, t34 * t115 - t35 * t116 + t128 * t345 - t129 * t346 - t14 * t175 + t15 * t176 + t72 * t83 - t73 * t84, t14 * t73 + t15 * t72 + t96 * t206 + t346 * t35 + t345 * t34 + (t191 - t278) * t112, -t106 * t21 + t40 * t63, t105 * t21 - t106 * t22 - t40 * t67 + t41 * t63, t21 * t253 + t40 * t226 + (qJD(2) * t106 - t63) * t305, t105 * t22 - t41 * t67, t22 * t253 + t41 * t226 + (-qJD(2) * t105 + t67) * t305 (-t226 - t307) * t305, t46 * t105 + t139 * t22 - t2 * t253 + t60 * t41 - t85 * t67 - t352 * t226 + (t67 * t287 + (qJD(2) * t28 + t10) * qJD(3)) * t250, -t260 * t253 + t46 * t106 - t139 * t21 - t60 * t40 - t85 * t63 + t353 * t226 + (t63 * t287 + (-qJD(2) * t29 - t11) * qJD(3)) * t250, t10 * t40 + t105 * t260 - t2 * t106 - t11 * t41 + t28 * t21 - t29 * t22 + t352 * t63 + t353 * t67, -t260 * t29 + t46 * t139 + t2 * t28 + (t85 - t278) * t60 + t353 * t11 + t352 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t229, t312 * t256, 0, t229, 0, 0, qJD(3) * t162 - t213 * t308 - t123, -t213 * t307 + (t161 + t193) * qJD(3) - t314, 0, 0, t246 * t266 (t333 + t198 * t244 + (-t240 + t241) * qJD(3)) * t307, t246 * t332 + t250 * t269, -t280 * t295, -t244 * t332 + t250 * t270, t229, -t336 - t162 * t196 + (t110 * t253 + t244 * t261 - t250 * t87) * qJD(2), t337 - t162 * t198 + (-t111 * t253 + t246 * t261 + t250 * t88) * qJD(2), t110 * t198 + t111 * t196 + (-qJD(4) * t196 + t307 * t87 + t59) * t246 + (qJD(4) * t198 + t307 * t88 - t58) * t244, -t123 * pkin(3) - t87 * t110 - t88 * t111 - t152 * t162 + (-t244 * t87 + t246 * t88) * qJD(4) + (-t58 * t244 + t59 * t246) * qJ(4), -t128 * t363 - t83 * t204, -t128 * t317 - t129 * t363 + t83 * t203 - t204 * t84, -t363 * t232 + (qJD(3) * t204 + t128) * t308, -t129 * t317 + t84 * t203, -t317 * t232 + (-qJD(3) * t203 + t129) * t308, t232 * t308, -t140 * t129 + t96 * t203 + t235 * t84 - t339 * t232 - t317 * t112 + (qJD(3) * t150 - t34) * t308, t140 * t128 + t96 * t204 - t235 * t83 - t344 * t232 + t363 * t112 + (-qJD(3) * t151 + t35) * t308, t128 * t339 + t129 * t344 - t14 * t203 - t15 * t204 + t150 * t83 - t151 * t84 + t317 * t35 - t34 * t363, -t112 * t140 + t14 * t151 + t15 * t150 + t96 * t235 + t339 * t34 - t344 * t35, -t21 * t134 + t343 * t63, t21 * t133 - t134 * t22 + t342 * t63 - t343 * t67, t343 * t226 + (qJD(3) * t134 + t63) * t308, t22 * t133 - t342 * t67, t342 * t226 + (-qJD(3) * t133 - t67) * t308, t226 * t308, t46 * t133 + t169 * t22 - t286 * t67 + t342 * t60 + t348 * t226 + (qJD(3) * t52 - t10) * t308, t46 * t134 - t169 * t21 - t286 * t63 - t343 * t60 - t347 * t226 + (-qJD(3) * t53 + t11) * t308, t10 * t343 - t11 * t342 + t133 * t260 - t2 * t134 + t52 * t21 - t53 * t22 - t347 * t67 - t348 * t63, -t10 * t348 - t11 * t347 + t46 * t169 + t2 * t52 - t260 * t53 + t286 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t266, t253 * t270, -t196 ^ 2 - t198 ^ 2, t196 * t88 + t198 * t87 + t123, 0, 0, 0, 0, 0, 0, t84 + t365, -t83 + t364, -t370 - t371, -t128 * t34 + t129 * t35 + t96, 0, 0, 0, 0, 0, 0, t22 + t372, -t21 - t373, -t374 - t375, -t10 * t63 - t11 * t67 + t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t334, -t370 + t371, -t83 - t364, t334, -t84 + t365, t236, t112 * t128 - t35 * t232 + t15, t112 * t129 - t232 * t34 - t14, 0, 0, t350, t361, t360, -t350, t354, t236, t12 * t226 + (-t128 * t67 + t226 * t301 + t236 * t351) * pkin(5) + t355, -t13 * t226 + (-t128 * t63 + t226 * t285 - t236 * t248) * pkin(5) + t359, t10 * t67 - t11 * t63 - t12 * t63 - t13 * t67 + (t351 * t21 - t22 * t248 + (-t248 * t63 + t351 * t67) * qJD(6)) * pkin(5), -t10 * t12 - t11 * t13 + (t351 * t2 - t260 * t248 + t128 * t60 + (-t10 * t248 + t11 * t351) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t350, t361, t360, -t350, t354, t236, -t11 * t226 + t355, -t10 * t226 + t359, 0, 0;];
tauc_reg  = t1;
