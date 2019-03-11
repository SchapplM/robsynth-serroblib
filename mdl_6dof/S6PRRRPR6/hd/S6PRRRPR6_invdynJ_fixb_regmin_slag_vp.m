% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRRPR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:36:50
% EndTime: 2019-03-08 23:37:03
% DurationCPUTime: 5.97s
% Computational Cost: add. (4325->547), mult. (9891->743), div. (0->0), fcn. (7828->12), ass. (0->263)
t194 = sin(pkin(11));
t199 = sin(qJ(2));
t203 = cos(qJ(2));
t338 = cos(pkin(11));
t339 = cos(pkin(6));
t252 = t339 * t338;
t122 = t194 * t203 + t199 * t252;
t272 = t194 * t339;
t124 = -t199 * t272 + t203 * t338;
t198 = sin(qJ(3));
t202 = cos(qJ(3));
t195 = sin(pkin(6));
t329 = t195 * t199;
t127 = t198 * t329 - t202 * t339;
t271 = t195 * t338;
t328 = t195 * t202;
t374 = g(3) * t127 - g(2) * (-t122 * t198 - t202 * t271) - g(1) * (-t124 * t198 + t194 * t328);
t201 = cos(qJ(4));
t197 = sin(qJ(4));
t307 = qJD(3) * t197;
t310 = qJD(2) * t198;
t142 = t201 * t310 + t307;
t295 = qJD(1) * qJD(2);
t109 = qJDD(2) * pkin(8) + (qJDD(1) * t199 + t203 * t295) * t195;
t312 = qJD(1) * t195;
t150 = qJD(2) * pkin(8) + t199 * t312;
t268 = qJD(1) * t339;
t174 = t198 * t268;
t266 = qJDD(1) * t339;
t305 = qJD(3) * t202;
t264 = -qJD(3) * t174 - t198 * t109 - t150 * t305 + t202 * t266;
t226 = qJDD(3) * pkin(3) + t264;
t294 = qJD(2) * qJD(3);
t276 = t202 * t294;
t293 = t198 * qJDD(2);
t297 = t201 * qJD(3);
t303 = qJD(4) * t198;
t387 = qJD(2) * t303 - qJDD(3);
t65 = -qJD(4) * t297 + (-t276 - t293) * t201 + t387 * t197;
t210 = -qJ(5) * t65 + qJD(5) * t142 + t226;
t365 = pkin(4) + pkin(5);
t309 = qJD(2) * t202;
t66 = (qJD(3) * (qJD(4) + t309) + t293) * t197 + t387 * t201;
t7 = -t365 * t66 + t210;
t391 = t7 + t374;
t319 = t202 * t203;
t107 = (t197 * t199 + t201 * t319) * t195;
t259 = pkin(3) * t198 - pkin(9) * t202;
t149 = t259 * qJD(3);
t237 = pkin(3) * t202 + pkin(9) * t198 + pkin(2);
t302 = qJD(4) * t201;
t390 = -qJD(1) * t107 + t197 * t149 - t237 * t302;
t389 = qJD(4) - t309;
t388 = qJD(6) - qJD(4);
t140 = t197 * t310 - t297;
t196 = sin(qJ(6));
t200 = cos(qJ(6));
t239 = -t200 * t140 + t142 * t196;
t72 = t140 * t196 + t142 * t200;
t386 = -t239 ^ 2 + t72 ^ 2;
t379 = -t198 * t150 + t202 * t268;
t260 = qJD(3) * pkin(3) + t379;
t218 = qJ(5) * t142 + t260;
t23 = -t140 * t365 + t218;
t128 = t198 * t339 + t199 * t328;
t327 = t195 * t203;
t86 = t128 * t197 + t201 * t327;
t87 = t128 * t201 - t197 * t327;
t243 = t196 * t87 - t200 * t86;
t189 = t202 * qJDD(2);
t138 = t198 * t294 + qJDD(4) - t189;
t286 = t203 * t312;
t104 = -qJD(2) * t237 - t286;
t304 = qJD(4) * t197;
t254 = t198 * t266;
t32 = qJDD(3) * pkin(9) + qJD(3) * t379 + t202 * t109 + t254;
t277 = t199 * t295;
t249 = -qJDD(1) * t327 + t195 * t277;
t54 = qJD(2) * t149 - qJDD(2) * t237 + t249;
t103 = t202 * t150 + t174;
t93 = qJD(3) * pkin(9) + t103;
t270 = t104 * t304 + t197 * t32 - t201 * t54 + t93 * t302;
t253 = qJDD(5) + t270;
t3 = pkin(10) * t65 - t138 * t365 + t253;
t131 = t138 * qJ(5);
t167 = t389 * qJD(5);
t229 = t104 * t302 + t197 * t54 + t201 * t32 - t304 * t93;
t6 = t131 + t167 + t229;
t5 = pkin(10) * t66 + t6;
t287 = -t196 * t5 + t200 * t3;
t121 = t194 * t199 - t203 * t252;
t79 = t122 * t202 - t198 * t271;
t43 = -t121 * t201 + t197 * t79;
t44 = t121 * t197 + t201 * t79;
t123 = t199 * t338 + t203 * t272;
t81 = t194 * t195 * t198 + t124 * t202;
t45 = -t123 * t201 + t197 * t81;
t46 = t123 * t197 + t201 * t81;
t385 = -g(3) * t243 + t23 * t72 - g(1) * (t196 * t46 - t200 * t45) - g(2) * (t196 * t44 - t200 * t43) - t287;
t383 = t72 * t239;
t39 = t201 * t104 - t197 * t93;
t317 = qJD(5) - t39;
t320 = t201 * t202;
t182 = pkin(8) * t320;
t377 = t197 * t319 - t199 * t201;
t382 = qJD(4) * t182 - t149 * t201 - t237 * t304 - t377 * t312;
t306 = qJD(3) * t198;
t381 = -qJ(5) * t306 - t390;
t380 = qJD(5) * t197 + t103;
t278 = t197 * t303;
t378 = t202 * t297 - t278;
t296 = -qJD(6) + t389;
t376 = -qJD(6) - t296;
t298 = qJD(6) * t200;
t299 = qJD(6) * t196;
t12 = t140 * t298 - t142 * t299 + t196 * t66 - t200 * t65;
t373 = t239 * t296 - t12;
t359 = pkin(9) * t138;
t41 = pkin(4) * t140 - t218;
t372 = -t389 * t41 + t359;
t205 = qJD(3) ^ 2;
t356 = g(2) * t121;
t357 = g(1) * t123;
t258 = t356 + t357;
t371 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t205 + t195 * (-g(3) * t203 + t277) - t249 + t258;
t336 = qJ(5) * t197;
t370 = -t201 * t365 - t336;
t26 = t196 * t86 + t200 * t87;
t367 = -g(1) * (t196 * t45 + t200 * t46) - g(2) * (t196 * t43 + t200 * t44) - t23 * t239 - g(3) * t26;
t366 = t142 ^ 2;
t364 = pkin(9) - pkin(10);
t360 = pkin(4) * t138;
t358 = pkin(10) * t198;
t40 = t197 * t104 + t201 * t93;
t354 = pkin(9) * qJD(4);
t353 = qJD(2) * pkin(2);
t169 = t389 * qJ(5);
t27 = t169 + t40;
t352 = t389 * t27;
t351 = t389 * t40;
t22 = pkin(10) * t140 + t40;
t16 = t169 + t22;
t350 = t196 * t16;
t349 = t197 * t65;
t288 = -pkin(8) * t197 - pkin(4);
t291 = pkin(10) * t320;
t348 = pkin(10) * t278 + (-t291 + (-pkin(5) + t288) * t198) * qJD(3) + t382;
t322 = t198 * t201;
t347 = -(-pkin(8) * qJD(3) + pkin(10) * qJD(4)) * t322 - (-qJD(5) + (-pkin(8) * qJD(4) + pkin(10) * qJD(3)) * t197) * t202 + t381;
t335 = qJ(5) * t201;
t231 = -t197 * t365 + t335;
t346 = t389 * t231 + t380;
t250 = pkin(4) * t197 - t335;
t345 = t389 * t250 - t380;
t344 = -qJD(5) * t202 + (-t198 * t297 - t202 * t304) * pkin(8) - t381;
t343 = t288 * t306 + t382;
t144 = t196 * t197 + t200 * t201;
t220 = t144 * t202;
t342 = -qJD(2) * t220 - t144 * t388;
t283 = t197 * t309;
t326 = t196 * t201;
t341 = t196 * t302 + t197 * t298 - t201 * t299 - t309 * t326 + (t283 - t304) * t200;
t146 = t259 * qJD(2);
t340 = t197 * t146 + t201 * t379;
t337 = qJ(5) * t140;
t333 = t140 * t389;
t332 = t142 * t140;
t331 = t142 * t389;
t330 = t142 * t201;
t325 = t197 * t200;
t324 = t197 * t202;
t323 = t198 * t200;
t318 = pkin(10) * t142 - t317;
t316 = qJDD(1) - g(3);
t314 = -t197 * t237 + t182;
t192 = t198 ^ 2;
t313 = -t202 ^ 2 + t192;
t311 = qJD(2) * t195;
t308 = qJD(3) * t142;
t300 = qJD(5) * t201;
t15 = -t365 * t389 - t318;
t292 = t15 * t298 + t196 * t3 + t200 * t5;
t50 = qJ(5) * t310 + t340;
t162 = t364 * t201;
t285 = t199 * t311;
t284 = t203 * t311;
t282 = t389 * t307;
t281 = t389 * t297;
t279 = t389 * t304;
t275 = t203 * t294;
t273 = -t196 * t65 - t200 * t66;
t90 = t197 * t379;
t269 = -t146 * t201 + t90;
t181 = pkin(8) * t324;
t267 = -t201 * t237 - t181;
t265 = t296 ^ 2;
t262 = t140 * t286;
t261 = t142 * t286;
t97 = -qJ(5) * t202 + t314;
t161 = t364 * t197;
t257 = pkin(10) * t283 - qJD(6) * t161 + t364 * t304 + t50;
t256 = (-t198 * t365 - t291) * qJD(2) + t269 + t388 * t162;
t251 = pkin(4) * t201 + t336;
t10 = t196 * t15 + t200 * t16;
t191 = t202 * pkin(4);
t67 = pkin(5) * t202 + t181 + t191 + (t237 - t358) * t201;
t73 = t197 * t358 + t97;
t244 = t196 * t67 + t200 * t73;
t24 = -pkin(4) * t389 + t317;
t242 = -t197 * t27 + t201 * t24;
t241 = qJ(5) * t200 - t196 * t365;
t240 = qJ(5) * t196 + t200 * t365;
t238 = -t325 + t326;
t206 = qJD(2) ^ 2;
t236 = qJDD(2) * t203 - t199 * t206;
t235 = pkin(3) + t251;
t234 = pkin(8) + t250;
t232 = t16 * t299 - t292;
t228 = t138 * t197 + t302 * t389;
t227 = t138 * t201 - t279;
t106 = t377 * t195;
t56 = -t121 * t324 - t122 * t201;
t58 = -t123 * t324 - t124 * t201;
t225 = g(1) * t58 + g(2) * t56 + g(3) * t106;
t57 = -t121 * t320 + t122 * t197;
t59 = -t123 * t320 + t124 * t197;
t224 = -g(1) * t59 - g(2) * t57 - g(3) * t107;
t222 = g(1) * t81 + g(2) * t79 + g(3) * t128;
t84 = -qJD(3) * t127 + t202 * t284;
t18 = qJD(4) * t87 + t197 * t84 - t201 * t285;
t85 = qJD(3) * t128 + t198 * t284;
t221 = t127 * t66 - t138 * t86 + t85 * t140 - t18 * t389;
t118 = t144 * t198;
t219 = -pkin(8) + t231;
t216 = -t260 * t389 - t359;
t19 = -qJD(4) * t86 + t197 * t285 + t201 * t84;
t215 = t127 * t65 + t138 * t87 - t142 * t85 + t19 * t389;
t214 = g(1) * t45 + g(2) * t43 + g(3) * t86 - t270;
t213 = -t354 * t389 + t374;
t13 = qJD(6) * t72 + t273;
t11 = pkin(4) * t66 - t210;
t211 = -t11 + t213;
t151 = -t286 - t353;
t209 = -pkin(8) * qJDD(3) + (t151 + t286 - t353) * qJD(3);
t208 = t142 * t41 + qJDD(5) - t214;
t207 = g(1) * t46 + g(2) * t44 + g(3) * t87 + t389 * t39 - t229;
t137 = pkin(3) - t370;
t132 = -qJDD(6) + t138;
t117 = t196 * t322 - t197 * t323;
t116 = t234 * t198;
t98 = t191 - t267;
t96 = t219 * t198;
t77 = pkin(4) * t142 + t337;
t55 = -t142 * t365 - t337;
t53 = (qJD(4) * t251 - t300) * t198 + t234 * t305;
t52 = -pkin(4) * t310 + t269;
t37 = -t65 + t333;
t35 = qJD(6) * t118 + t378 * t196 - t302 * t323 - t305 * t325;
t34 = -t198 * t238 * t388 + qJD(3) * t220;
t30 = (t370 * qJD(4) + t300) * t198 + t219 * t305;
t9 = t15 * t200 - t350;
t8 = t253 - t360;
t1 = [t316, 0, t236 * t195 (-qJDD(2) * t199 - t203 * t206) * t195, 0, 0, 0, 0, 0, -qJD(3) * t85 - qJDD(3) * t127 + (-t198 * t275 + t202 * t236) * t195, -qJD(3) * t84 - qJDD(3) * t128 + (-t198 * t236 - t202 * t275) * t195, 0, 0, 0, 0, 0, t221, -t215, t221, -t140 * t19 + t142 * t18 - t65 * t86 - t66 * t87, t215, t11 * t127 + t18 * t24 + t19 * t27 + t41 * t85 + t6 * t87 + t8 * t86 - g(3), 0, 0, 0, 0, 0 -(-qJD(6) * t26 + t18 * t200 - t19 * t196) * t296 + t243 * t132 - t85 * t239 - t127 * t13 (-qJD(6) * t243 + t18 * t196 + t19 * t200) * t296 + t26 * t132 - t85 * t72 - t127 * t12; 0, qJDD(2), t316 * t327 + t258, g(1) * t124 + g(2) * t122 - t316 * t329, qJDD(2) * t192 + 0.2e1 * t198 * t276, 0.2e1 * t189 * t198 - 0.2e1 * t294 * t313, qJDD(3) * t198 + t202 * t205, qJDD(3) * t202 - t198 * t205, 0, t209 * t198 + t371 * t202, -t371 * t198 + t209 * t202, t378 * t142 - t65 * t322 (-t140 * t201 - t142 * t197) * t305 + (t349 - t201 * t66 + (t140 * t197 - t330) * qJD(4)) * t198 (t65 + t281) * t202 + (t227 + t308) * t198 (t66 - t282) * t202 + (-qJD(3) * t140 - t228) * t198, -t138 * t202 + t306 * t389, t267 * t138 - t382 * t389 + ((pkin(8) * t140 - t197 * t260) * qJD(3) + t270) * t202 + (-t262 - t260 * t302 + t39 * qJD(3) - t226 * t197 + (t66 + t282) * pkin(8)) * t198 + t224, -t314 * t138 - t390 * t389 + (-t260 * t297 + (t279 + t308) * pkin(8) + t229) * t202 + (-t261 + t260 * t304 - t40 * qJD(3) - t226 * t201 + (-t65 + t281) * pkin(8)) * t198 + t225, t116 * t66 - t138 * t98 + t140 * t53 + (t307 * t41 + t8) * t202 - t343 * t389 + (-qJD(3) * t24 + t11 * t197 + t302 * t41 - t262) * t198 + t224, -t65 * t98 - t66 * t97 + t343 * t142 - t344 * t140 + t242 * t305 + (-g(3) * t327 - t197 * t6 + t201 * t8 + (-t197 * t24 - t201 * t27) * qJD(4) + t258) * t198, t116 * t65 + t138 * t97 - t142 * t53 + (-t297 * t41 - t6) * t202 + t344 * t389 + (qJD(3) * t27 - t11 * t201 + t304 * t41 + t261) * t198 - t225, t6 * t97 + t11 * t116 + t41 * t53 + t8 * t98 - g(1) * (pkin(4) * t59 + pkin(8) * t124 + qJ(5) * t58) - g(2) * (pkin(4) * t57 + pkin(8) * t122 + qJ(5) * t56) - g(3) * (pkin(4) * t107 + qJ(5) * t106) + t344 * t27 + t343 * t24 + t237 * t357 + t237 * t356 + (-g(3) * pkin(8) * t199 + (-t41 * t198 * qJD(1) - g(3) * t237) * t203) * t195, t118 * t12 + t34 * t72, -t117 * t12 - t118 * t13 - t239 * t34 - t35 * t72, -t118 * t132 + t12 * t202 - t296 * t34 - t306 * t72, t117 * t132 - t13 * t202 + t239 * t306 + t296 * t35, -t132 * t202 + t296 * t306 -(-t196 * t73 + t200 * t67) * t132 + t287 * t202 + t30 * t239 + t96 * t13 + t7 * t117 + t23 * t35 - g(1) * (t196 * t58 + t200 * t59) - g(2) * (t196 * t56 + t200 * t57) - g(3) * (t106 * t196 + t107 * t200) + (-qJD(3) * t9 + t239 * t286) * t198 - (t196 * t347 + t200 * t348) * t296 + (-t10 * t202 + t244 * t296) * qJD(6), t244 * t132 + t232 * t202 + t30 * t72 + t96 * t12 + t7 * t118 + t23 * t34 - g(1) * (-t196 * t59 + t200 * t58) - g(2) * (-t196 * t57 + t200 * t56) - g(3) * (t106 * t200 - t107 * t196) + (qJD(3) * t10 + t286 * t72) * t198 - ((-qJD(6) * t67 + t347) * t200 + (qJD(6) * t73 - t348) * t196) * t296; 0, 0, 0, 0, -t198 * t206 * t202, t313 * t206, t293, t189, qJDD(3), qJD(3) * t103 - t151 * t310 + t264 + t374, -t254 + (-qJD(2) * t151 - t109) * t202 + t222, t330 * t389 - t349 (-t65 - t333) * t201 + (-t331 - t66) * t197 (-t142 * t198 - t320 * t389) * qJD(2) + t228 (t140 * t198 + t324 * t389) * qJD(2) + t227, -t389 * t310, -t39 * t310 - pkin(3) * t66 - t103 * t140 + t90 * t389 + t216 * t197 + (t226 - (t146 + t354) * t389 + t374) * t201, pkin(3) * t65 + t340 * t389 + t40 * t310 - t103 * t142 + t216 * t201 + (-t213 - t226) * t197, t345 * t140 - t372 * t197 + t211 * t201 - t235 * t66 + t24 * t310 + t389 * t52, t140 * t50 - t142 * t52 + (t6 + t389 * t24 + (qJD(4) * t142 - t66) * pkin(9)) * t201 + (t8 - t352 + (qJD(4) * t140 - t65) * pkin(9)) * t197 - t222, -t345 * t142 + t211 * t197 + t372 * t201 - t235 * t65 - t27 * t310 - t389 * t50, -t24 * t52 - t27 * t50 + t345 * t41 + (qJD(4) * t242 + t8 * t197 + t6 * t201 - t222) * pkin(9) + (-t11 + t374) * t235, -t12 * t238 + t342 * t72, -t12 * t144 + t13 * t238 - t239 * t342 - t341 * t72, t132 * t238 - t296 * t342 + t310 * t72, t132 * t144 - t239 * t310 + t296 * t341, -t296 * t310 -(t161 * t200 - t162 * t196) * t132 + t137 * t13 + t9 * t310 + t346 * t239 + t341 * t23 - (t196 * t257 - t200 * t256) * t296 + t391 * t144 (t161 * t196 + t162 * t200) * t132 + t137 * t12 - t10 * t310 + t346 * t72 + t342 * t23 - (t196 * t256 + t200 * t257) * t296 - t391 * t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t332, -t140 ^ 2 + t366, t37, -t66 + t331, t138, t142 * t260 + t214 + t351, -t140 * t260 + t207, -t140 * t77 - t208 + t351 + 0.2e1 * t360, pkin(4) * t65 - qJ(5) * t66 + (t27 - t40) * t142 + (t24 - t317) * t140, -t140 * t41 + t142 * t77 + 0.2e1 * t131 + 0.2e1 * t167 - t207, t6 * qJ(5) - t8 * pkin(4) - t41 * t77 - t24 * t40 - g(1) * (-pkin(4) * t45 + qJ(5) * t46) - g(2) * (-pkin(4) * t43 + qJ(5) * t44) - g(3) * (-pkin(4) * t86 + qJ(5) * t87) + t317 * t27, -t383, -t386, t373, t296 * t72 + t13, t132, t240 * t132 - t55 * t239 - (t196 * t318 - t200 * t22) * t296 + (t241 * t296 + t10) * qJD(6) + t385, t241 * t132 - t55 * t72 - (t196 * t22 + t200 * t318) * t296 + (-t240 * t296 - t350) * qJD(6) + t292 + t367; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t138 + t332, t37, -t389 ^ 2 - t366, t208 - t352 - t360, 0, 0, 0, 0, 0, -t200 * t132 - t142 * t239 - t196 * t265, t196 * t132 - t142 * t72 - t200 * t265; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t383, t386, -t373, t376 * t72 - t273, -t132, t376 * t10 - t385, -t296 * t9 + t232 - t367;];
tau_reg  = t1;
