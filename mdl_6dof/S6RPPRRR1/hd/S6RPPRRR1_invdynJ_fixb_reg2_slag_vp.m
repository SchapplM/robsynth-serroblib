% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRRR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:18:47
% EndTime: 2019-03-09 02:18:57
% DurationCPUTime: 6.66s
% Computational Cost: add. (12231->510), mult. (27255->628), div. (0->0), fcn. (20887->18), ass. (0->269)
t210 = cos(qJ(6));
t207 = sin(qJ(6));
t202 = sin(pkin(11));
t343 = sin(qJ(4));
t273 = t343 * t202;
t256 = qJD(1) * t273;
t204 = cos(pkin(11));
t264 = qJDD(1) * t343;
t345 = cos(qJ(4));
t265 = qJDD(1) * t345;
t270 = qJD(4) * t345;
t296 = qJD(1) * t204;
t276 = t202 * t265 + t204 * t264 + t270 * t296;
t102 = qJD(4) * t256 - t276;
t239 = t202 * t264 - t204 * t265;
t230 = -t345 * t202 - t343 * t204;
t350 = t230 * qJD(1);
t361 = qJD(4) * t350;
t103 = t239 - t361;
t231 = -t345 * t204 + t273;
t144 = t231 * qJD(1);
t208 = sin(qJ(5));
t344 = cos(qJ(5));
t268 = qJD(5) * t344;
t295 = qJD(5) * t208;
t233 = -t344 * t102 - t208 * t103 - t144 * t268 + t295 * t350;
t235 = t208 * t144 + t344 * t350;
t289 = qJD(4) + qJD(5);
t258 = t210 * t289;
t288 = qJDD(4) + qJDD(5);
t294 = qJD(6) * t207;
t29 = -qJD(6) * t258 - t207 * t288 - t210 * t233 - t235 * t294;
t27 = t210 * t29;
t293 = qJD(6) * t210;
t87 = t207 * t289 - t210 * t235;
t312 = qJD(6) * t87;
t30 = t207 * t233 - t210 * t288 + t312;
t85 = -t207 * t235 - t258;
t329 = -t207 * t30 - t85 * t293;
t94 = -t344 * t144 + t208 * t350;
t364 = qJD(6) - t94;
t374 = t207 * t364;
t353 = t87 * t374;
t366 = t210 * t94;
t376 = t366 * t85 - t27 + t329 - t353;
t260 = t208 * t102 - t344 * t103;
t360 = qJD(5) * t235;
t56 = -t260 - t360;
t205 = cos(pkin(10));
t183 = -t205 * pkin(1) - pkin(2);
t195 = t204 * pkin(3);
t352 = t183 - t195;
t139 = qJDD(1) * t352 + qJDD(3);
t340 = t103 * pkin(4);
t82 = t139 + t340;
t12 = t56 * pkin(5) - pkin(9) * t233 + t82;
t11 = t210 * t12;
t203 = sin(pkin(10));
t177 = pkin(1) * t203 + qJ(3);
t164 = t177 * qJD(1);
t189 = t204 * qJD(2);
t117 = t189 + (-pkin(7) * qJD(1) - t164) * t202;
t128 = t202 * qJD(2) + t204 * t164;
t118 = pkin(7) * t296 + t128;
t232 = -t343 * t117 - t345 * t118;
t71 = -t144 * pkin(8) - t232;
t281 = t344 * t71;
t274 = t343 * t118;
t75 = t345 * t117 - t274;
t70 = pkin(8) * t350 + t75;
t68 = qJD(4) * pkin(4) + t70;
t42 = t208 * t68 + t281;
t37 = pkin(9) * t289 + t42;
t141 = qJD(1) * t352 + qJD(3);
t101 = t144 * pkin(4) + t141;
t57 = -pkin(5) * t94 + pkin(9) * t235 + t101;
t14 = t207 * t57 + t210 * t37;
t151 = qJD(1) * qJD(3) + qJDD(1) * t177;
t187 = t204 * qJDD(2);
t113 = t187 + (-pkin(7) * qJDD(1) - t151) * t202;
t120 = t202 * qJDD(2) + t204 * t151;
t290 = t204 * qJDD(1);
t114 = pkin(7) * t290 + t120;
t250 = t345 * t113 - t343 * t114;
t48 = qJD(4) * t232 + t250;
t39 = qJDD(4) * pkin(4) + t102 * pkin(8) + t48;
t269 = qJD(4) * t343;
t278 = t343 * t113 + t345 * t114 + t117 * t270;
t47 = -t118 * t269 + t278;
t40 = -t103 * pkin(8) + t47;
t9 = t208 * t39 + t68 * t268 - t71 * t295 + t344 * t40;
t7 = pkin(9) * t288 + t9;
t3 = -qJD(6) * t14 - t207 * t7 + t11;
t338 = t14 * t364;
t375 = t3 + t338;
t13 = -t207 * t37 + t210 * t57;
t339 = t13 * t364;
t291 = qJDD(1) * t183;
t161 = qJDD(3) + t291;
t201 = qJ(1) + pkin(10);
t191 = sin(t201);
t193 = cos(t201);
t267 = -g(1) * t191 + g(2) * t193;
t373 = -t161 - t267;
t354 = t289 * t94;
t372 = t233 - t354;
t332 = t94 ^ 2;
t333 = t235 ^ 2;
t371 = -t332 + t333;
t26 = t29 * t207;
t369 = -t26 + (t293 - t366) * t87;
t53 = qJDD(6) + t56;
t50 = t207 * t53;
t326 = t293 * t364 + t50;
t334 = t87 * t235;
t368 = -t364 * t366 + t326 + t334;
t319 = t208 * t71;
t41 = t344 * t68 - t319;
t36 = -pkin(5) * t289 - t41;
t367 = t36 * t94;
t336 = t85 * t235;
t365 = t364 * t235;
t331 = t94 * t235;
t200 = pkin(11) + qJ(4);
t194 = qJ(5) + t200;
t180 = sin(t194);
t308 = t180 * t193;
t309 = t180 * t191;
t363 = g(1) * t308 + g(2) * t309;
t301 = t235 * qJD(4);
t362 = -t301 + t260;
t173 = g(3) * t180;
t181 = cos(t194);
t306 = t181 * t193;
t307 = t181 * t191;
t277 = -g(1) * t306 - g(2) * t307 - t173;
t359 = -t101 * t94 - t277 - t9;
t32 = t36 * t294;
t358 = t13 * t235 + t210 * t363 + t32;
t341 = g(3) * t181;
t266 = t208 * t40 - t344 * t39;
t10 = -qJD(5) * t42 - t266;
t8 = -pkin(5) * t288 - t10;
t355 = t8 * t207 + t36 * t293;
t357 = -t14 * t235 + t207 * t341 + t355;
t356 = t101 * t235 - t266 - t341 + t363;
t65 = -pkin(5) * t235 - pkin(9) * t94;
t51 = t210 * t53;
t238 = t294 * t364 - t51;
t313 = pkin(1) * qJDD(1);
t298 = t181 * pkin(5) + t180 * pkin(9);
t351 = -t13 * t207 + t14 * t210;
t227 = t208 * t231;
t107 = -t230 * t344 - t227;
t221 = t344 * t231;
t224 = qJD(4) * t230;
t72 = -t208 * t224 + t221 * t289 - t230 * t295;
t318 = t210 * t72;
t349 = -t107 * t238 - t318 * t364;
t348 = t107 * t288 - t289 * t72;
t347 = t350 ^ 2;
t346 = qJD(4) ^ 2;
t2 = qJD(6) * t13 + t207 * t12 + t210 * t7;
t1 = t2 * t210;
t209 = sin(qJ(1));
t337 = t209 * pkin(1);
t335 = t87 * t85;
t206 = -pkin(7) - qJ(3);
t330 = pkin(7) + t177;
t182 = t195 + pkin(2);
t28 = t30 * t210;
t328 = -t107 * t28 + t85 * t318;
t106 = -t208 * t230 + t221;
t225 = qJD(4) * t231;
t73 = -qJD(5) * t227 - t208 * t225 - t344 * t224 - t230 * t268;
t327 = -t29 * t106 + t87 * t73;
t325 = -t107 * t56 - t72 * t94;
t321 = t207 * t72;
t320 = t207 * t85;
t317 = t210 * t87;
t314 = t103 * t230 + t144 * t225;
t311 = qJD(6) * t364;
t310 = t350 * t144;
t305 = t191 * t207;
t304 = t191 * t210;
t303 = t193 * t207;
t302 = t193 * t210;
t147 = t330 * t202;
t148 = t330 * t204;
t97 = -t343 * t147 + t345 * t148;
t192 = cos(t200);
t179 = pkin(4) * t192;
t154 = t179 + t182;
t211 = cos(qJ(1));
t196 = t211 * pkin(1);
t299 = t193 * t154 + t196;
t198 = t202 ^ 2;
t199 = t204 ^ 2;
t297 = t198 + t199;
t285 = pkin(9) * t311;
t283 = t87 * t321;
t184 = pkin(4) * t208 + pkin(9);
t280 = t184 * t311;
t275 = -t8 - t341;
t272 = qJD(3) * t345;
t271 = qJD(3) * t343;
t259 = t1 + t277;
t257 = pkin(4) * t268;
t43 = t208 * t70 + t281;
t255 = pkin(4) * t295 - t43;
t254 = -pkin(9) * t53 - t367;
t253 = g(1) * t193 + g(2) * t191;
t251 = g(1) * t209 - g(2) * t211;
t96 = -t345 * t147 - t343 * t148;
t249 = -t106 * t30 - t73 * t85;
t248 = t106 * t233 - t235 * t73;
t197 = -pkin(8) + t206;
t247 = -t193 * t197 - t337;
t246 = t13 * t210 + t14 * t207;
t83 = pkin(8) * t230 + t96;
t84 = -pkin(8) * t231 + t97;
t60 = t208 * t83 + t344 * t84;
t116 = pkin(4) * t231 + t352;
t62 = t106 * pkin(5) - t107 * pkin(9) + t116;
t23 = t207 * t62 + t210 * t60;
t22 = -t207 * t60 + t210 * t62;
t119 = -t151 * t202 + t187;
t244 = -t119 * t202 + t120 * t204;
t243 = (-t164 * t202 + t189) * t202 - t128 * t204;
t241 = t374 * t94 - t238;
t240 = -qJD(6) * t57 + t173 - t7;
t237 = t253 * t180;
t190 = sin(t200);
t236 = t253 * t190;
t228 = -t291 + t373;
t226 = -g(3) * t192 + t236;
t222 = pkin(4) * t224;
t220 = -t184 * t53 - t257 * t364 - t367;
t219 = -t326 * t107 + t374 * t72;
t218 = -qJD(6) * t246 - t3 * t207 + t1;
t77 = -t147 * t270 - t148 * t269 - t202 * t271 + t204 * t272;
t217 = -g(1) * (-pkin(5) * t308 + pkin(9) * t306) - g(2) * (-pkin(5) * t309 + pkin(9) * t307);
t216 = t139 + t267;
t215 = t102 * t231 - t224 * t350;
t214 = pkin(8) * t225 + t147 * t269 - t148 * t270 - t202 * t272 - t204 * t271;
t185 = -t344 * pkin(4) - pkin(5);
t142 = t144 ^ 2;
t124 = t181 * t302 + t305;
t123 = -t181 * t303 + t304;
t122 = -t181 * t304 + t303;
t121 = t181 * t305 + t302;
t105 = -qJDD(4) * t231 + t346 * t230;
t104 = -qJDD(4) * t230 - t346 * t231;
t78 = t230 * qJD(3) - t97 * qJD(4);
t74 = pkin(8) * t224 + t77;
t61 = -pkin(4) * t350 + t65;
t59 = t208 * t84 - t344 * t83;
t54 = -t106 * t288 - t289 * t73;
t44 = t344 * t70 - t319;
t31 = t73 * pkin(5) + t72 * pkin(9) - t222;
t20 = t207 * t65 + t210 * t41;
t19 = -t207 * t41 + t210 * t65;
t18 = t207 * t61 + t210 * t44;
t17 = -t207 * t44 + t210 * t61;
t16 = t60 * qJD(5) + t208 * t74 - t344 * t214;
t15 = t208 * t214 + t83 * t268 - t84 * t295 + t344 * t74;
t5 = -qJD(6) * t23 - t207 * t15 + t210 * t31;
t4 = qJD(6) * t22 + t210 * t15 + t207 * t31;
t6 = [0, 0, 0, 0, 0, qJDD(1), t251, g(1) * t211 + g(2) * t209, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t205 * t313 - t267, -0.2e1 * t203 * t313 + t253, 0 (t251 + (t203 ^ 2 + t205 ^ 2) * t313) * pkin(1), t198 * qJDD(1), 0.2e1 * t202 * t290, 0, t199 * qJDD(1), 0, 0, t228 * t204, -t228 * t202, t151 * t297 + t244 - t253, t161 * t183 - g(1) * (-pkin(2) * t191 + qJ(3) * t193 - t337) - g(2) * (pkin(2) * t193 + qJ(3) * t191 + t196) + t244 * t177 - t243 * qJD(3), t102 * t230 + t225 * t350, t215 + t314, t104, t103 * t231 - t144 * t224, t105, 0, t96 * qJDD(4) + t352 * t103 + t139 * t231 - t267 * t192 + (-t141 * t230 + t78) * qJD(4), -t77 * qJD(4) - t97 * qJDD(4) - t102 * t352 - t139 * t230 - t141 * t225 + t190 * t267, t96 * t102 - t97 * t103 - t77 * t144 + t78 * t350 + t48 * t230 - t224 * t232 - t253 + (t75 * qJD(4) - t47) * t231, t47 * t97 - t232 * t77 + t48 * t96 + t75 * t78 + t139 * t352 - g(1) * (-t182 * t191 - t193 * t206 - t337) - g(2) * (t182 * t193 - t191 * t206 + t196) t107 * t233 + t235 * t72, -t248 + t325, t348, t106 * t56 - t73 * t94, t54, 0, g(1) * t307 - g(2) * t306 + t101 * t73 + t82 * t106 + t116 * t56 - t16 * t289 + t222 * t94 - t288 * t59, -g(1) * t309 + g(2) * t308 - t101 * t72 + t82 * t107 + t116 * t233 - t15 * t289 + t222 * t235 - t288 * t60, -t10 * t107 - t106 * t9 + t15 * t94 - t16 * t235 + t233 * t59 + t41 * t72 - t42 * t73 - t56 * t60 - t253, t9 * t60 + t42 * t15 - t10 * t59 - t41 * t16 + t82 * t116 - t101 * t222 - g(1) * (-t154 * t191 + t247) - g(2) * (-t191 * t197 + t299) -t72 * t317 + (-t294 * t87 - t27) * t107, t283 + (t26 + (-t317 + t320) * qJD(6)) * t107 + t328, t327 + t349, -t329 * t107 - t72 * t320, t219 + t249, t106 * t53 + t364 * t73, -g(1) * t122 - g(2) * t124 + t3 * t106 + t107 * t355 + t13 * t73 + t16 * t85 + t22 * t53 + t59 * t30 - t36 * t321 + t364 * t5, -t36 * t318 - g(1) * t121 - g(2) * t123 - t2 * t106 - t14 * t73 + t16 * t87 - t23 * t53 - t59 * t29 - t4 * t364 + (t8 * t210 - t32) * t107, t22 * t29 - t23 * t30 - t4 * t85 - t5 * t87 + t246 * t72 - t267 * t180 + (-qJD(6) * t351 - t2 * t207 - t3 * t210) * t107, t2 * t23 + t14 * t4 + t3 * t22 + t13 * t5 + t8 * t59 + t36 * t16 - g(1) * t247 - g(2) * (pkin(5) * t306 + pkin(9) * t308 + t299) + (-g(1) * (-t154 - t298) + g(2) * t197) * t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t119 * t204 + t120 * t202 - g(3), 0, 0, 0, 0, 0, 0, t105, -t104, -t215 + t314, t224 * t75 + t225 * t232 - t230 * t47 - t231 * t48 - g(3), 0, 0, 0, 0, 0, 0, t54, -t348, t248 + t325, -t10 * t106 + t107 * t9 - t41 * t73 - t42 * t72 - g(3), 0, 0, 0, 0, 0, 0, t219 - t249, t327 - t349, -t283 + (-t26 + (t317 + t320) * qJD(6)) * t107 + t328, t8 * t106 + t107 * t218 - t351 * t72 + t36 * t73 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t290, t202 * qJDD(1), -t297 * qJD(1) ^ 2, qJD(1) * t243 - t373, 0, 0, 0, 0, 0, 0, t239 - 0.2e1 * t361 (-t256 - t144) * qJD(4) + t276, -t142 - t347, -t144 * t232 - t350 * t75 + t216, 0, 0, 0, 0, 0, 0, -t260 - t301 - 0.2e1 * t360, t233 + t354, -t332 - t333, -t235 * t41 - t42 * t94 + t216 + t340, 0, 0, 0, 0, 0, 0, t241 + t336, -t210 * t364 ^ 2 + t334 - t50 (t85 * t94 + t29) * t210 + t353 + t329, t36 * t235 + t375 * t210 + (t2 - t339) * t207 + t267; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t310, -t142 + t347 (-t256 + t144) * qJD(4) + t276, t310, -t239, qJDD(4), t141 * t350 + t226 + t250, g(3) * t190 + t141 * t144 + t253 * t192 + (t274 + t75) * qJD(4) - t278, 0, 0, t331, t371, t372, -t331, t362, t288, t43 * qJD(4) + (t43 - t42) * qJD(5) + (t344 * t288 - t289 * t295 - t350 * t94) * pkin(4) + t356, t44 * t289 + (-t208 * t288 - t235 * t350 - t268 * t289) * pkin(4) + t359, t41 * t94 - t42 * t235 + t43 * t235 - t44 * t94 + (-t344 * t233 - t208 * t56 + (-t208 * t235 + t344 * t94) * qJD(5)) * pkin(4), t41 * t43 - t42 * t44 + (t344 * t10 + t101 * t350 + t208 * t9 + (-t208 * t41 + t344 * t42) * qJD(5) + t226) * pkin(4), t369, t376, t368, t374 * t85 - t28, t241 - t336, t365, -t17 * t364 + t185 * t30 + t255 * t85 + (t275 - t280) * t210 + t220 * t207 + t358, t18 * t364 - t185 * t29 + t255 * t87 + t220 * t210 + (-t237 + t280) * t207 + t357, t17 * t87 + t18 * t85 + (-t85 * t257 + t13 * t94 - t184 * t30 + (t184 * t87 - t13) * qJD(6)) * t210 + (t87 * t257 + t14 * t94 - t184 * t29 - t3 + (t184 * t85 - t14) * qJD(6)) * t207 + t259, t8 * t185 - t14 * t18 - t13 * t17 - t36 * t43 - g(3) * (t179 + t298) + (t236 + (t208 * t36 + t344 * t351) * qJD(5)) * pkin(4) + t218 * t184 + t217; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t331, t371, t372, -t331, t362, t288, t42 * qJD(4) + t356, t289 * t41 + t359, 0, 0, t369, t376, t368, t320 * t364 - t28, -t364 * t374 - t336 + t51, t365, -pkin(5) * t30 - t19 * t364 - t42 * t85 + t254 * t207 + (t275 - t285) * t210 + t358, pkin(5) * t29 + t20 * t364 - t42 * t87 + t254 * t210 + (-t237 + t285) * t207 + t357, t19 * t87 + t20 * t85 + (-t339 + (-t30 + t312) * pkin(9)) * t210 + ((qJD(6) * t85 - t29) * pkin(9) - t375) * t207 + t259, -t8 * pkin(5) + pkin(9) * t218 - g(3) * t298 - t13 * t19 - t14 * t20 - t36 * t42 + t217; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t335, -t85 ^ 2 + t87 ^ 2, t364 * t85 - t29, -t335, t364 * t87 - t30, t53, -g(1) * t123 + g(2) * t121 + t207 * t240 - t293 * t37 - t36 * t87 + t11 + t338, g(1) * t124 - g(2) * t122 + t339 + t36 * t85 + (qJD(6) * t37 - t12) * t207 + t240 * t210, 0, 0;];
tau_reg  = t6;
