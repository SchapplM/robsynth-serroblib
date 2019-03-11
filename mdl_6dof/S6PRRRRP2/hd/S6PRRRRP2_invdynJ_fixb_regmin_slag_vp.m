% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRRRP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:04:23
% EndTime: 2019-03-09 00:04:39
% DurationCPUTime: 6.69s
% Computational Cost: add. (7630->520), mult. (17094->682), div. (0->0), fcn. (13446->14), ass. (0->288)
t207 = sin(qJ(3));
t206 = sin(qJ(4));
t210 = cos(qJ(3));
t324 = t206 * t210;
t378 = cos(qJ(4));
t167 = t207 * t378 + t324;
t304 = qJD(3) + qJD(4);
t122 = t304 * t167;
t291 = t378 * t210;
t306 = t207 * qJDD(2);
t82 = qJD(2) * t122 - qJDD(2) * t291 + t206 * t306;
t78 = qJDD(5) + t82;
t284 = qJD(2) * t378;
t315 = qJD(2) * t207;
t162 = t206 * t315 - t210 * t284;
t205 = sin(qJ(5));
t311 = qJD(5) * t205;
t397 = t162 * t205 + t311;
t247 = -t206 * t207 + t291;
t211 = cos(qJ(2));
t203 = sin(pkin(6));
t318 = qJD(1) * t203;
t289 = t211 * t318;
t129 = t247 * t289;
t212 = -pkin(9) - pkin(8);
t292 = qJD(3) * t212;
t168 = t207 * t292;
t169 = t210 * t292;
t176 = t212 * t207;
t177 = t212 * t210;
t248 = t176 * t378 + t206 * t177;
t68 = qJD(4) * t248 + t168 * t378 + t206 * t169;
t396 = -t68 + t129;
t209 = cos(qJ(5));
t263 = t209 * pkin(5) + t205 * qJ(6);
t173 = -pkin(4) - t263;
t313 = qJD(4) * t206;
t208 = sin(qJ(2));
t290 = t208 * t318;
t273 = -qJD(2) * t212 + t290;
t204 = cos(pkin(6));
t317 = qJD(1) * t204;
t127 = t207 * t317 + t210 * t273;
t119 = t378 * t127;
t126 = -t273 * t207 + t210 * t317;
t64 = t206 * t126 + t119;
t267 = pkin(3) * t313 - t64;
t310 = qJD(5) * t209;
t336 = t162 * t209;
t395 = -qJD(6) * t205 + (-t310 - t336) * qJ(6) + t397 * pkin(5);
t349 = cos(pkin(11));
t275 = t349 * t208;
t202 = sin(pkin(11));
t329 = t202 * t211;
t157 = t204 * t275 + t329;
t201 = qJ(3) + qJ(4);
t197 = sin(t201);
t198 = cos(t201);
t276 = t203 * t349;
t110 = -t157 * t197 - t198 * t276;
t274 = t349 * t211;
t330 = t202 * t208;
t159 = -t204 * t330 + t274;
t331 = t202 * t203;
t112 = -t159 * t197 + t198 * t331;
t328 = t203 * t208;
t145 = -t197 * t328 + t198 * t204;
t393 = g(1) * t112 + g(2) * t110 + g(3) * t145;
t155 = qJD(5) + t162;
t364 = qJD(3) * pkin(3);
t120 = t126 + t364;
t283 = qJD(4) * t378;
t307 = t204 * qJDD(1);
t182 = t210 * t307;
t309 = qJD(1) * qJD(2);
t144 = qJDD(2) * pkin(8) + (qJDD(1) * t208 + t211 * t309) * t203;
t272 = pkin(9) * qJDD(2) + t144;
t52 = qJDD(3) * pkin(3) - qJD(3) * t127 - t272 * t207 + t182;
t57 = qJD(3) * t126 + t207 * t307 + t272 * t210;
t226 = t120 * t283 - t127 * t313 + t206 * t52 + t378 * t57;
t303 = qJDD(3) + qJDD(4);
t12 = pkin(10) * t303 + t226;
t282 = t208 * t309;
t326 = t203 * t211;
t261 = -qJDD(1) * t326 + t203 * t282;
t348 = qJDD(2) * pkin(2);
t143 = t261 - t348;
t308 = qJD(2) * qJD(3);
t281 = t207 * t308;
t190 = pkin(3) * t281;
t121 = t304 * t247;
t216 = t121 * qJD(2) + t167 * qJDD(2);
t305 = t210 * qJDD(2);
t27 = -pkin(3) * t305 + t82 * pkin(4) - pkin(10) * t216 + t143 + t190;
t63 = t206 * t120 + t119;
t55 = pkin(10) * t304 + t63;
t195 = t210 * pkin(3) + pkin(2);
t153 = -qJD(2) * t195 - t289;
t164 = -qJD(2) * t324 - t207 * t284;
t81 = pkin(4) * t162 + pkin(10) * t164 + t153;
t250 = t209 * t12 + t205 * t27 + t81 * t310 - t311 * t55;
t366 = qJ(6) * t78;
t2 = qJD(6) * t155 + t250 + t366;
t279 = t205 * t12 - t209 * t27 + t55 * t310 + t81 * t311;
t382 = pkin(5) * t78;
t4 = qJDD(6) + t279 - t382;
t392 = t2 * t209 + t4 * t205;
t268 = pkin(3) * t283;
t356 = t209 * t78;
t391 = pkin(10) * (t155 * t311 - t356);
t352 = t395 + t267;
t108 = -pkin(4) * t247 - pkin(10) * t167 - t195;
t135 = t206 * t176 - t177 * t378;
t298 = t207 * t364;
t60 = pkin(4) * t122 - pkin(10) * t121 + t298;
t390 = -t108 * t310 + t135 * t311 + t396 * t209 + (t290 - t60) * t205;
t193 = pkin(3) * t206 + pkin(10);
t230 = -t193 * t311 + t209 * t268;
t152 = t164 * qJ(6);
t118 = t206 * t127;
t65 = t126 * t378 - t118;
t107 = -pkin(4) * t164 + pkin(10) * t162;
t91 = pkin(3) * t315 + t107;
t368 = t205 * t91 + t209 * t65;
t34 = -t152 + t368;
t389 = -t34 + t230;
t231 = t167 * t326;
t350 = -qJD(1) * t231 + qJD(4) * t135 + t206 * t168 - t169 * t378;
t156 = -t204 * t274 + t330;
t158 = t204 * t329 + t275;
t266 = g(1) * t158 + g(2) * t156;
t235 = -g(3) * t326 + t266;
t388 = t235 * t197;
t320 = t205 * t108 + t209 * t135;
t180 = t205 * t326;
t160 = t204 * t210 - t207 * t328;
t327 = t203 * t210;
t161 = t204 * t207 + t208 * t327;
t93 = t206 * t160 + t161 * t378;
t80 = t209 * t93 - t180;
t387 = t290 - t298;
t111 = t157 * t198 - t197 * t276;
t113 = t159 * t198 + t197 * t331;
t146 = t197 * t204 + t198 * t328;
t240 = g(1) * t113 + g(2) * t111 + g(3) * t146;
t386 = -t240 + t392;
t213 = qJD(3) ^ 2;
t385 = -pkin(8) * t213 + t203 * (-g(3) * t211 + t282) - t143 + t266 + t348;
t239 = t209 * t164 - t205 * t304;
t44 = -qJD(5) * t239 + t205 * t216 - t209 * t303;
t384 = t239 ^ 2;
t383 = t155 ^ 2;
t381 = pkin(10) * t78;
t380 = qJ(6) * t122 - qJD(6) * t247 - t390;
t89 = t129 * t205 - t209 * t290;
t379 = -pkin(5) * t122 + qJD(5) * t320 + t205 * t68 - t209 * t60 - t89;
t377 = pkin(5) * t164;
t270 = t120 * t313 + t127 * t283 + t206 * t57 - t378 * t52;
t13 = -pkin(4) * t303 + t270;
t269 = t209 * t304;
t43 = -qJD(5) * t269 - t164 * t311 - t205 * t303 - t209 * t216;
t5 = t44 * pkin(5) + t43 * qJ(6) + qJD(6) * t239 + t13;
t370 = t5 * t205;
t62 = t120 * t378 - t118;
t369 = t205 * t107 + t209 * t62;
t367 = pkin(10) * qJD(5);
t365 = qJD(2) * pkin(2);
t32 = t205 * t81 + t209 * t55;
t363 = t155 * t32;
t130 = -t164 * t205 - t269;
t54 = -pkin(4) * t304 - t62;
t36 = t130 * pkin(5) + qJ(6) * t239 + t54;
t362 = t162 * t36;
t361 = t193 * t78;
t31 = -t205 * t55 + t209 * t81;
t322 = qJD(6) - t31;
t22 = -pkin(5) * t155 + t322;
t360 = t205 * t22;
t359 = t205 * t43;
t358 = t205 * t78;
t357 = t209 * t44;
t354 = t54 * t162;
t262 = pkin(5) * t205 - qJ(6) * t209;
t353 = t262 * t121 + (qJD(5) * t263 - qJD(6) * t209) * t167 + t350;
t351 = -t63 + t395;
t347 = t121 * t205;
t346 = t121 * t209;
t345 = t130 * t155;
t344 = t130 * t205;
t343 = t130 * t209;
t342 = t239 * t130;
t341 = t239 * t155;
t340 = t239 * t205;
t339 = t239 * t209;
t338 = t155 * t164;
t335 = t164 * t162;
t334 = t167 * t209;
t333 = t198 * t205;
t332 = t198 * t209;
t323 = t209 * t211;
t321 = qJDD(1) - g(3);
t199 = t207 ^ 2;
t319 = -t210 ^ 2 + t199;
t316 = qJD(2) * t203;
t314 = qJD(2) * t208;
t312 = qJD(5) * t193;
t301 = t173 * t110;
t300 = t173 * t112;
t299 = t378 * pkin(3);
t296 = t203 * t323;
t21 = t22 * t310;
t30 = t36 * t311;
t295 = t36 * t310;
t47 = t54 * t311;
t294 = t393 * t205;
t293 = t173 * t145;
t287 = t203 * t314;
t286 = t211 * t316;
t280 = t211 * t308;
t278 = -t22 * t164 + t30;
t277 = t31 * t164 + t47;
t271 = t155 * t209;
t265 = g(1) * t159 + g(2) * t157;
t259 = t354 - t361;
t23 = qJ(6) * t155 + t32;
t258 = -t205 * t23 + t209 * t22;
t214 = qJD(2) ^ 2;
t255 = qJDD(2) * t211 - t208 * t214;
t254 = pkin(4) * t198 + pkin(10) * t197 + t195;
t253 = -g(1) * t202 + g(2) * t349;
t79 = t205 * t93 + t296;
t249 = t160 * t378 - t206 * t161;
t246 = -t167 * t311 + t346;
t123 = qJD(3) * t160 + t210 * t286;
t124 = -qJD(3) * t161 - t207 * t286;
t40 = qJD(4) * t249 + t123 * t378 + t206 * t124;
t16 = qJD(5) * t80 + t205 * t40 - t209 * t287;
t41 = qJD(4) * t93 + t206 * t123 - t124 * t378;
t244 = t41 * t130 - t155 * t16 - t249 * t44 - t78 * t79;
t139 = t180 * t198 - t209 * t328;
t84 = -t156 * t333 - t157 * t209;
t86 = -t158 * t333 - t159 * t209;
t243 = g(1) * t86 + g(2) * t84 + g(3) * t139;
t140 = (t198 * t323 + t205 * t208) * t203;
t85 = -t156 * t332 + t157 * t205;
t87 = -t158 * t332 + t159 * t205;
t242 = -g(1) * t87 - g(2) * t85 - g(3) * t140;
t238 = t13 * t205 - t32 * t164 + t54 * t310 + t294;
t237 = -t393 - t5;
t236 = -t13 - t393;
t171 = -t289 - t365;
t233 = -qJD(2) * t171 - t144 + t265;
t232 = t23 * t164 - t36 * t336 - t294 - t370;
t15 = -qJD(5) * t79 + t205 * t287 + t209 * t40;
t229 = -t15 * t155 - t239 * t41 + t249 * t43 - t78 * t80;
t116 = t146 * t205 + t296;
t70 = t111 * t205 - t156 * t209;
t72 = t113 * t205 - t158 * t209;
t228 = g(1) * t72 + g(2) * t70 + g(3) * t116 - t279;
t225 = -pkin(8) * qJDD(3) + (t171 + t289 - t365) * qJD(3);
t224 = -t359 - t357 + (-t339 + t344) * qJD(5);
t117 = t146 * t209 - t180;
t71 = t111 * t209 + t156 * t205;
t73 = t113 * t209 + t158 * t205;
t223 = -g(1) * t73 - g(2) * t71 - g(3) * t117 + t250;
t222 = t153 * t164 - t270 - t393;
t221 = -t239 * t36 + qJDD(6) - t228;
t220 = t22 * t336 - t23 * t397 + t21 + t386;
t218 = t153 * t162 - t226 + t240;
t194 = -t299 - pkin(4);
t165 = -t299 + t173;
t109 = -qJDD(2) * t195 + t190 + t261;
t83 = -t162 ^ 2 + t164 ^ 2;
t74 = -pkin(5) * t239 + qJ(6) * t130;
t67 = t167 * t262 - t248;
t59 = -t164 * t304 - t82;
t58 = t162 * t304 + t216;
t46 = pkin(5) * t247 - t108 * t209 + t135 * t205;
t45 = -qJ(6) * t247 + t320;
t38 = -t107 * t209 + t205 * t62 + t377;
t37 = -t152 + t369;
t35 = t205 * t65 - t209 * t91 + t377;
t28 = -t43 + t345;
t20 = t155 * t271 - t164 * t239 + t358;
t19 = -t130 * t164 - t205 * t383 + t356;
t17 = -t239 * t271 - t359;
t6 = (-t43 - t345) * t209 + (-t44 + t341) * t205;
t1 = [t321, 0, t255 * t203 (-qJDD(2) * t208 - t211 * t214) * t203, 0, 0, 0, 0, 0, qJD(3) * t124 + qJDD(3) * t160 + (-t207 * t280 + t210 * t255) * t203, -qJD(3) * t123 - qJDD(3) * t161 + (-t207 * t255 - t210 * t280) * t203, 0, 0, 0, 0, 0, -t41 * t304 + t249 * t303 + (t162 * t314 - t211 * t82) * t203, -t40 * t304 - t93 * t303 - qJDD(2) * t231 + (-t121 * t211 - t208 * t164) * t316, 0, 0, 0, 0, 0, t244, t229, t244, -t130 * t15 - t16 * t239 - t43 * t79 - t44 * t80, -t229, t15 * t23 + t16 * t22 + t2 * t80 - t249 * t5 + t36 * t41 + t4 * t79 - g(3); 0, qJDD(2), t321 * t326 + t266, -t321 * t328 + t265, qJDD(2) * t199 + 0.2e1 * t210 * t281, 0.2e1 * t207 * t305 - 0.2e1 * t308 * t319, qJDD(3) * t207 + t210 * t213, qJDD(3) * t210 - t207 * t213, 0, t225 * t207 + t210 * t385, -t207 * t385 + t225 * t210, -t164 * t121 + t167 * t216, -t121 * t162 + t164 * t122 - t167 * t82 + t216 * t247, t121 * t304 + t167 * t303, -t122 * t304 + t247 * t303, 0, -t109 * t247 + t153 * t122 - t162 * t387 - t195 * t82 + t235 * t198 + t248 * t303 - t304 * t350, t109 * t167 + t153 * t121 - t135 * t303 + t387 * t164 - t195 * t216 + t304 * t396 - t388, -t239 * t246 - t334 * t43 (t340 - t343) * t121 + (t359 - t357 + (t339 + t344) * qJD(5)) * t167, -t122 * t239 + t155 * t246 + t247 * t43 + t334 * t78, -t167 * t358 - t122 * t130 + t247 * t44 + (-t167 * t310 - t347) * t155, t122 * t155 - t247 * t78, t279 * t247 + t31 * t122 - t248 * t44 + t89 * t155 + t350 * t130 + ((-qJD(5) * t135 + t60) * t155 + t108 * t78 + t54 * qJD(5) * t167) * t209 + ((-qJD(5) * t108 - t68) * t155 - t135 * t78 + t13 * t167 + t54 * t121) * t205 + t242, -t320 * t78 + t250 * t247 - t32 * t122 + t248 * t43 + t54 * t346 + (t13 * t209 - t47) * t167 + t390 * t155 - t350 * t239 + t243, t36 * t347 - t122 * t22 + t247 * t4 + t44 * t67 - t46 * t78 + (t295 + t370) * t167 - t379 * t155 + t353 * t130 + t242, -t43 * t46 - t44 * t45 - t379 * t239 - t380 * t130 + t258 * t121 + t388 + (-t2 * t205 + t209 * t4 + (-t209 * t23 - t360) * qJD(5)) * t167, -t36 * t346 + t122 * t23 - t247 * t2 + t43 * t67 + t45 * t78 + (-t5 * t209 + t30) * t167 + t380 * t155 + t353 * t239 - t243, t2 * t45 + t5 * t67 + t4 * t46 - g(1) * (pkin(5) * t87 + qJ(6) * t86 - t159 * t212) - g(2) * (pkin(5) * t85 + qJ(6) * t84 - t157 * t212) + t353 * t36 + t380 * t23 + t379 * t22 + t266 * t254 + (-pkin(5) * t140 - qJ(6) * t139 - (-t208 * t212 + t211 * t254) * t203) * g(3); 0, 0, 0, 0, -t207 * t214 * t210, t319 * t214, t306, t305, qJDD(3), -g(3) * t160 + t207 * t233 + t253 * t327 + t182, g(3) * t161 + (-t203 * t253 - t307) * t207 + t233 * t210, -t335, t83, t58, t59, t303, t64 * t304 + (-t162 * t315 + t303 * t378 - t304 * t313) * pkin(3) + t222, t65 * t304 + (t164 * t315 - t206 * t303 - t283 * t304) * pkin(3) + t218, t17, t6, t20, t19, t338, t194 * t44 + t267 * t130 + ((-t268 + t65) * t155 + t259) * t205 + ((-t91 - t312) * t155 + t236) * t209 + t277, -t194 * t43 + t259 * t209 - t267 * t239 + (-t230 + t368) * t155 + t238, t35 * t155 + t165 * t44 + t352 * t130 + (-t155 * t268 - t361 + t362) * t205 + (-t155 * t312 + t237) * t209 + t278, t34 * t130 + t35 * t239 + t224 * t193 + t220 + (-t340 - t343) * t268, t165 * t43 + (-qJD(5) * t36 + t361) * t209 + t352 * t239 + t389 * t155 + t232, t5 * t165 + t268 * t360 - t22 * t35 - g(1) * (pkin(10) * t113 + (-t159 * t207 + t202 * t327) * pkin(3) - t300) - g(2) * (t111 * pkin(10) + (-t157 * t207 - t210 * t276) * pkin(3) - t301) - g(3) * (pkin(3) * t160 + pkin(10) * t146 - t293) + t352 * t36 + t389 * t23 + (t21 + t392) * t193; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t335, t83, t58, t59, t303, t304 * t63 + t222, t304 * t62 + t218, t17, t6, t20, t19, t338, -pkin(4) * t44 - t63 * t130 + (t155 * t62 + t354 - t381) * t205 + ((-t107 - t367) * t155 + t236) * t209 + t277, pkin(4) * t43 + t155 * t369 + t239 * t63 + t336 * t54 + t238 + t391, t155 * t38 + t173 * t44 + (t362 - t381) * t205 + t351 * t130 + (-t155 * t367 + t237) * t209 + t278, pkin(10) * t224 + t130 * t37 + t239 * t38 + t220, -t155 * t37 + t173 * t43 + t239 * t351 + t232 - t295 - t391, t5 * t173 - t23 * t37 - t22 * t38 + g(1) * t300 + g(2) * t301 + g(3) * t293 + t351 * t36 + (qJD(5) * t258 + t386) * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t342, -t130 ^ 2 + t384, t28, -t44 - t341, t78, t239 * t54 + t228 + t363, t130 * t54 + t155 * t31 - t223, -t130 * t74 - t221 + t363 + 0.2e1 * t382, pkin(5) * t43 - qJ(6) * t44 - (t23 - t32) * t239 + (t22 - t322) * t130, 0.2e1 * t366 - t130 * t36 - t239 * t74 + (0.2e1 * qJD(6) - t31) * t155 + t223, t2 * qJ(6) - t4 * pkin(5) - t36 * t74 - t22 * t32 - g(1) * (-pkin(5) * t72 + qJ(6) * t73) - g(2) * (-pkin(5) * t70 + qJ(6) * t71) - g(3) * (-pkin(5) * t116 + qJ(6) * t117) + t322 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t342 - t78, t28, -t383 - t384, -t155 * t23 + t221 - t382;];
tau_reg  = t1;
