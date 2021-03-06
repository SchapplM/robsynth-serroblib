% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRR10_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR10_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR10_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR10_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_invdynJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:26:13
% EndTime: 2019-12-05 17:26:37
% DurationCPUTime: 10.87s
% Computational Cost: add. (8452->618), mult. (22238->896), div. (0->0), fcn. (18955->14), ass. (0->283)
t214 = sin(pkin(6));
t220 = sin(qJ(3));
t224 = cos(qJ(3));
t257 = t214 * (pkin(3) * t220 - pkin(9) * t224);
t215 = sin(pkin(5));
t221 = sin(qJ(2));
t351 = t215 * t221;
t309 = qJD(1) * t351;
t402 = qJD(3) * t257 - t214 * t309;
t216 = cos(pkin(6));
t225 = cos(qJ(2));
t340 = t224 * t225;
t346 = t220 * t221;
t254 = -t216 * t346 + t340;
t140 = t254 * t215;
t354 = t214 * t220;
t205 = pkin(8) * t354;
t348 = t216 * t224;
t166 = pkin(2) * t348 - t205;
t156 = qJD(3) * t166;
t338 = -qJD(1) * t140 + t156;
t333 = qJD(1) * t225;
t182 = qJD(2) * pkin(2) + t215 * t333;
t217 = cos(pkin(5));
t334 = qJD(1) * t217;
t250 = t182 * t216 + t214 * t334;
t223 = cos(qJ(4));
t401 = pkin(9) * t223;
t349 = t216 * t220;
t352 = t214 * t224;
t167 = pkin(2) * t349 + pkin(8) * t352;
t149 = pkin(9) * t216 + t167;
t270 = -pkin(3) * t224 - pkin(9) * t220;
t150 = (-pkin(2) + t270) * t214;
t219 = sin(qJ(4));
t326 = qJD(4) * t223;
t327 = qJD(4) * t219;
t375 = -t149 * t327 + t150 * t326 + t402 * t219 + t338 * t223;
t344 = t221 * t224;
t345 = t220 * t225;
t256 = t216 * t344 + t345;
t139 = t256 * t215;
t157 = qJD(3) * t167;
t337 = -qJD(1) * t139 + t157;
t319 = qJDD(2) * t216;
t204 = qJDD(3) + t319;
t350 = t215 * t225;
t201 = qJDD(1) * t350;
t307 = qJD(2) * t351;
t275 = qJD(1) * t307;
t144 = qJDD(2) * pkin(2) + t201 - t275;
t332 = qJD(2) * t214;
t169 = pkin(8) * t332 + t309;
t320 = qJDD(1) * t217;
t299 = t214 * t320;
t328 = qJD(3) * t224;
t378 = pkin(8) * t214;
t394 = -qJDD(2) * t378 - (qJD(2) * t333 + qJDD(1) * t221) * t215 - t250 * qJD(3);
t27 = t144 * t348 - t169 * t328 + t220 * t394 + t224 * t299;
t23 = -pkin(3) * t204 - t27;
t331 = qJD(2) * t216;
t284 = qJD(3) + t331;
t258 = qJD(4) * t284;
t318 = qJDD(2) * t220;
t298 = t214 * t318;
t304 = t214 * t328;
t312 = t219 * t354;
t396 = qJD(4) * t312 - t223 * t304;
t66 = -t219 * t204 + (-t258 - t298) * t223 + t396 * qJD(2);
t303 = t219 * t328;
t67 = -t223 * t204 + t214 * (qJD(2) * (t220 * t326 + t303) + t219 * t318) + t219 * t258;
t12 = pkin(4) * t67 + pkin(10) * t66 + t23;
t218 = sin(qJ(5));
t222 = cos(qJ(5));
t330 = qJD(2) * t224;
t202 = t214 * t330;
t264 = t202 - qJD(4);
t341 = t224 * t169;
t86 = t220 * t250 + t341;
t71 = pkin(9) * t284 + t86;
t199 = t216 * t334;
t98 = t199 + (qJD(2) * t270 - t182) * t214;
t33 = t219 * t98 + t223 * t71;
t31 = -pkin(10) * t264 + t33;
t308 = t220 * t332;
t136 = t219 * t308 - t223 * t284;
t138 = t219 * t284 + t223 * t308;
t85 = -t220 * t169 + t224 * t250;
t70 = -pkin(3) * t284 - t85;
t34 = t136 * pkin(4) - t138 * pkin(10) + t70;
t262 = t218 * t31 - t222 * t34;
t321 = qJD(2) * qJD(3);
t301 = t220 * t321;
t317 = qJDD(2) * t224;
t200 = t214 * t317;
t316 = qJDD(4) - t200;
t240 = t214 * t301 + t316;
t329 = qJD(3) * t220;
t233 = -t144 * t349 + t169 * t329 - t220 * t299 + t224 * t394;
t22 = pkin(9) * t204 - t233;
t198 = t216 * t320;
t241 = t301 - t317;
t300 = t224 * t321;
t242 = t300 + t318;
t56 = t198 + (pkin(3) * t241 - pkin(9) * t242 - t144) * t214;
t253 = -t219 * t56 - t223 * t22 - t98 * t326 + t327 * t71;
t3 = pkin(10) * t240 - t253;
t1 = -t262 * qJD(5) + t218 * t12 + t222 * t3;
t127 = qJD(5) + t136;
t400 = t262 * t127 + t1;
t305 = t214 * t329;
t399 = pkin(10) * t305 + t375;
t111 = -t216 * t326 + t396;
t353 = t214 * t223;
t165 = t216 * t219 + t220 * t353;
t112 = qJD(4) * t165 + t214 * t303;
t398 = t112 * pkin(4) + t111 * pkin(10) + t337;
t385 = t216 * t340 - t346;
t107 = -t215 * t385 - t217 * t352;
t361 = sin(pkin(11));
t290 = t361 * t225;
t362 = cos(pkin(11));
t293 = t362 * t221;
t162 = t217 * t293 + t290;
t291 = t361 * t221;
t292 = t362 * t225;
t231 = -t217 * t292 + t291;
t228 = t231 * t224;
t295 = t215 * t362;
t272 = t214 * t295;
t72 = t162 * t220 + t216 * t228 + t224 * t272;
t163 = -t217 * t291 + t292;
t232 = t217 * t290 + t293;
t230 = t232 * t224;
t294 = t215 * t361;
t271 = t214 * t294;
t74 = t163 * t220 + t216 * t230 - t224 * t271;
t245 = g(1) * t74 + g(2) * t72 + g(3) * t107;
t285 = qJD(4) * t264;
t397 = pkin(9) * t285 - t23 + t245;
t32 = -t219 * t71 + t223 * t98;
t395 = t264 * t32 - t253;
t154 = qJD(2) * t257;
t55 = t219 * t154 + t223 * t85;
t393 = pkin(9) * t327 + pkin(10) * t308 + t55;
t269 = pkin(4) * t219 - pkin(10) * t223;
t392 = -qJD(5) * t401 + t269 * qJD(4) - t182 * t349 - t341 - (t220 * t334 + t269 * t330) * t214;
t14 = t218 * t34 + t222 * t31;
t2 = -qJD(5) * t14 + t222 * t12 - t218 * t3;
t390 = -t14 * t127 - t2;
t386 = t223 * t149 + t219 * t150;
t374 = -qJD(4) * t386 - t338 * t219 + t402 * t223;
t388 = t136 * t264;
t387 = t138 * t264;
t101 = t222 * t138 - t218 * t264;
t325 = qJD(5) * t101;
t25 = -t218 * t66 - t222 * t240 + t325;
t384 = -pkin(4) * t223 - pkin(10) * t219;
t383 = -pkin(9) * t240 - t264 * t70;
t296 = t219 * t22 - t223 * t56;
t6 = -qJD(4) * t33 - t296;
t226 = qJD(2) ^ 2;
t148 = t205 + (-pkin(2) * t224 - pkin(3)) * t216;
t164 = -t223 * t216 + t312;
t78 = t164 * pkin(4) - t165 * pkin(10) + t148;
t80 = -pkin(10) * t352 + t386;
t28 = -t218 * t80 + t222 * t78;
t381 = qJD(5) * t28 + t218 * t398 + t222 * t399;
t29 = t218 * t78 + t222 * t80;
t380 = -qJD(5) * t29 - t218 * t399 + t222 * t398;
t376 = -pkin(4) * t305 - t374;
t195 = -pkin(3) + t384;
t323 = qJD(5) * t222;
t373 = t195 * t323 + t218 * t392 - t222 * t393;
t324 = qJD(5) * t218;
t372 = -t195 * t324 + t218 * t393 + t222 * t392;
t288 = t222 * t264;
t99 = t138 * t218 + t288;
t371 = t101 * t99;
t370 = t127 * t99;
t63 = qJDD(5) + t67;
t367 = t218 * t63;
t366 = t218 * t99;
t365 = t222 * t63;
t24 = qJD(5) * t288 + t138 * t324 - t218 * t240 + t222 * t66;
t364 = t24 * t218;
t363 = t25 * t222;
t360 = t101 * t127;
t359 = t138 * t136;
t161 = -t214 * t350 + t216 * t217;
t358 = t161 * t214;
t210 = t214 ^ 2;
t356 = t210 * t226;
t355 = t214 * t219;
t347 = t218 * t223;
t343 = t221 * t226;
t342 = t222 * t224;
t339 = qJDD(1) - g(3);
t314 = t214 * t351;
t336 = pkin(2) * t350 + pkin(8) * t314;
t212 = t220 ^ 2;
t213 = t224 ^ 2;
t335 = t212 - t213;
t313 = t218 * t352;
t302 = t214 * t216 * t226;
t287 = t224 * t264;
t286 = t127 * t222;
t283 = qJD(3) + 0.2e1 * t331;
t282 = t204 + t319;
t281 = t210 * t215 * t343;
t280 = t220 * t224 * t356;
t278 = t214 * t307;
t273 = t220 * t300;
t268 = g(1) * t163 + g(2) * t162;
t122 = t202 * t347 - t222 * t308;
t267 = -t218 * t326 + t122;
t123 = (t218 * t220 + t223 * t342) * t332;
t266 = t222 * t326 - t123;
t263 = -t14 * t218 + t222 * t262;
t255 = t216 * t345 + t344;
t108 = t215 * t255 + t217 * t354;
t77 = t108 * t223 + t161 * t219;
t37 = t107 * t222 - t218 * t77;
t38 = t107 * t218 + t222 * t77;
t54 = t154 * t223 - t219 * t85;
t76 = t108 * t219 - t161 * t223;
t87 = -t149 * t219 + t150 * t223;
t259 = t140 * pkin(3) + pkin(9) * t139 + t336;
t113 = t165 * t218 + t214 * t342;
t30 = pkin(4) * t264 - t32;
t251 = -pkin(10) * t63 + t127 * t30;
t249 = t219 * t264;
t109 = t214 * t231 - t216 * t295;
t110 = t214 * t232 + t216 * t294;
t229 = t231 * t220;
t73 = t162 * t224 - t216 * t229 - t220 * t272;
t75 = t163 * t224 + (-t216 * t232 + t271) * t220;
t248 = -g(1) * (t110 * t223 - t219 * t75) - g(2) * (t109 * t223 - t219 * t73) + g(3) * t76;
t40 = t109 * t219 + t223 * t73;
t42 = t110 * t219 + t223 * t75;
t247 = -g(1) * t42 - g(2) * t40 - g(3) * t77;
t102 = t140 * t219 - t223 * t314;
t93 = -t162 * t349 - t228;
t57 = -t162 * t353 + t219 * t93;
t95 = -t163 * t349 - t230;
t59 = -t163 * t353 + t219 * t95;
t246 = -g(1) * t59 - g(2) * t57 - g(3) * t102;
t244 = g(1) * t75 + g(2) * t73 + g(3) * t108;
t92 = t162 * t348 - t229;
t94 = t163 * t348 - t220 * t232;
t243 = g(1) * t94 + g(2) * t92 + g(3) * t139;
t4 = -pkin(4) * t240 - t6;
t239 = t248 - t4;
t158 = t231 * pkin(2);
t238 = t93 * pkin(3) + pkin(9) * t92 + t162 * t378 - t158;
t159 = t232 * pkin(2);
t237 = t95 * pkin(3) + pkin(9) * t94 + t163 * t378 - t159;
t227 = pkin(10) * qJD(5) * t127 - t239;
t146 = t195 * t218 + t222 * t401;
t145 = -pkin(9) * t347 + t195 * t222;
t125 = -t214 * t182 + t199;
t114 = t165 * t222 - t313;
t106 = t107 * pkin(3);
t104 = -t214 * t144 + t198;
t103 = t140 * t223 + t219 * t314;
t84 = pkin(4) * t138 + pkin(10) * t136;
t79 = pkin(4) * t352 - t87;
t69 = t74 * pkin(3);
t68 = t72 * pkin(3);
t65 = t217 * t304 + (t254 * qJD(2) + qJD(3) * t385) * t215;
t64 = t217 * t305 + (qJD(2) * t256 + qJD(3) * t255) * t215;
t60 = t163 * t355 + t223 * t95;
t58 = t162 * t355 + t223 * t93;
t49 = -qJD(5) * t313 - t111 * t218 + t165 * t323 - t222 * t305;
t48 = qJD(5) * t113 + t222 * t111 - t218 * t305;
t45 = -pkin(4) * t308 - t54;
t20 = -qJD(4) * t76 + t219 * t278 + t65 * t223;
t19 = qJD(4) * t77 + t65 * t219 - t223 * t278;
t16 = t218 * t84 + t222 * t32;
t15 = -t218 * t32 + t222 * t84;
t10 = qJD(5) * t37 + t20 * t222 + t64 * t218;
t9 = -qJD(5) * t38 - t20 * t218 + t64 * t222;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t339, 0, 0, 0, 0, 0, 0, (qJDD(2) * t225 - t343) * t215, (-qJDD(2) * t221 - t225 * t226) * t215, 0, -g(3) + (t217 ^ 2 + (t221 ^ 2 + t225 ^ 2) * t215 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, -t107 * t204 - t224 * t281 + t241 * t358 - t284 * t64, -t108 * t204 + t220 * t281 + t242 * t358 - t284 * t65, ((t107 * t220 + t108 * t224) * qJDD(2) + (t220 * t64 + t224 * t65 + (t107 * t224 - t108 * t220) * qJD(3)) * qJD(2)) * t214, t104 * t161 - t107 * t27 - t108 * t233 + t125 * t278 - t64 * t85 + t65 * t86 - g(3), 0, 0, 0, 0, 0, 0, t107 * t67 + t64 * t136 + t19 * t264 - t240 * t76, -t107 * t66 + t64 * t138 + t20 * t264 - t240 * t77, -t136 * t20 + t138 * t19 - t66 * t76 - t67 * t77, t107 * t23 - t19 * t32 + t20 * t33 - t253 * t77 - t6 * t76 + t64 * t70 - g(3), 0, 0, 0, 0, 0, 0, t127 * t9 + t19 * t99 + t25 * t76 + t37 * t63, -t10 * t127 + t101 * t19 - t24 * t76 - t38 * t63, -t10 * t99 - t101 * t9 + t24 * t37 - t25 * t38, t1 * t38 + t10 * t14 + t19 * t30 + t2 * t37 - t262 * t9 + t4 * t76 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), g(1) * t232 + g(2) * t231 - g(3) * t350 + t201, -t339 * t351 + t268, 0, 0, (qJDD(2) * t212 + 0.2e1 * t273) * t210, 0.2e1 * (t220 * t317 - t321 * t335) * t210, (t220 * t282 + t283 * t328) * t214, (qJDD(2) * t213 - 0.2e1 * t273) * t210, (t224 * t282 - t283 * t329) * t214, t204 * t216, t166 * t204 + t27 * t216 - g(1) * t95 - g(2) * t93 - g(3) * t140 + (-t104 * t224 + t125 * t329) * t214 + (-pkin(2) * t241 + t224 * t275) * t210 - t337 * t284, -t167 * t204 + t233 * t216 + (t104 * t220 + t125 * t328) * t214 + (-pkin(2) * t242 - t220 * t275) * t210 + t243 - t338 * t284, (-g(3) * t351 + (-qJD(3) * t85 + qJDD(2) * t167 - t233) * t224 + (-qJD(3) * t86 - qJDD(2) * t166 - t27) * t220 + ((t338 - t156) * t224 + (t337 - t157) * t220) * qJD(2) - t268) * t214, -t233 * t167 + t27 * t166 + g(1) * t159 + g(2) * t158 - g(3) * t336 + t338 * t86 - t337 * t85 + (-t104 * pkin(2) - pkin(8) * t268 - t125 * t309) * t214, -t111 * t138 - t165 * t66, t111 * t136 - t112 * t138 + t164 * t66 - t165 * t67, t111 * t264 + t165 * t316 + (t66 * t224 + (qJD(2) * t165 + t138) * t329) * t214, t112 * t136 + t164 * t67, t112 * t264 - t164 * t316 + (t67 * t224 + (-qJD(2) * t164 - t136) * t329) * t214, (-t316 * t224 + (-t202 - t264) * t329) * t214, t87 * t316 + t148 * t67 + t23 * t164 + t70 * t112 - g(1) * t60 - g(2) * t58 - g(3) * t103 + (-t6 * t224 + (qJD(2) * t87 + t32) * t329) * t214 + t337 * t136 - t374 * t264, -t386 * t316 - t148 * t66 + t23 * t165 - t70 * t111 + (-t253 * t224 + (-qJD(2) * t386 - t33) * t329) * t214 + t337 * t138 - t246 + t375 * t264, t32 * t111 - t33 * t112 - t136 * t375 - t138 * t374 + t164 * t253 - t6 * t165 - t386 * t67 + t87 * t66 - t243, -g(1) * t237 - g(2) * t238 - g(3) * t259 + t23 * t148 - t253 * t386 + t32 * t374 + t33 * t375 + t337 * t70 + t6 * t87, -t101 * t48 - t114 * t24, -t101 * t49 + t113 * t24 - t114 * t25 + t48 * t99, t101 * t112 + t114 * t63 - t127 * t48 - t164 * t24, t113 * t25 + t49 * t99, -t112 * t99 - t113 * t63 - t127 * t49 - t164 * t25, t112 * t127 + t164 * t63, t28 * t63 + t2 * t164 - t262 * t112 + t79 * t25 + t4 * t113 + t30 * t49 - g(1) * (t218 * t94 + t222 * t60) - g(2) * (t218 * t92 + t222 * t58) - g(3) * (t103 * t222 + t139 * t218) + t376 * t99 + t380 * t127, -t29 * t63 - t1 * t164 - t14 * t112 - t79 * t24 + t4 * t114 - t30 * t48 - g(1) * (-t218 * t60 + t222 * t94) - g(2) * (-t218 * t58 + t222 * t92) - g(3) * (-t103 * t218 + t139 * t222) - t381 * t127 + t376 * t101, -t1 * t113 - t101 * t380 - t114 * t2 - t14 * t49 + t24 * t28 - t25 * t29 - t262 * t48 - t381 * t99 + t246, t1 * t29 + t2 * t28 + t4 * t79 - g(1) * (pkin(4) * t60 + pkin(10) * t59 + t237) - g(2) * (pkin(4) * t58 + pkin(10) * t57 + t238) - g(3) * (pkin(4) * t103 + pkin(10) * t102 + t259) + t376 * t30 + t381 * t14 - t380 * t262; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t280, t335 * t356, -t224 * t302 + t298, t280, t220 * t302 + t200, t204, -t125 * t308 + t284 * t86 + t245 + t27, -t125 * t202 + t284 * t85 + t233 + t244, 0, 0, -t66 * t219 - t223 * t387, (-t66 + t388) * t223 + (-t67 + t387) * t219, -t223 * t285 + t219 * t316 + (t223 * t287 + (qJD(3) * t219 - t138) * t220) * t332, -t136 * t249 - t67 * t223, t219 * t285 + t223 * t316 + (-t219 * t287 + (qJD(3) * t223 + t136) * t220) * t332, t264 * t308, -pkin(3) * t67 - t86 * t136 + t383 * t219 + t223 * t397 + t54 * t264 - t32 * t308, pkin(3) * t66 - t86 * t138 - t219 * t397 + t383 * t223 - t55 * t264 + t33 * t308, t55 * t136 + t54 * t138 + ((qJD(4) * t138 - t67) * pkin(9) + t395) * t223 + (-t6 + t264 * t33 + (qJD(4) * t136 - t66) * pkin(9)) * t219 - t244, -t23 * pkin(3) + g(1) * t69 + g(2) * t68 + g(3) * t106 - t32 * t54 - t33 * t55 - t70 * t86 + (-t6 * t219 - t253 * t223 + (-t219 * t33 - t223 * t32) * qJD(4) - t244) * pkin(9), -t24 * t222 * t219 + (-t219 * t324 + t266) * t101, t101 * t122 + t123 * t99 + (-t101 * t218 - t222 * t99) * t326 + (t364 - t363 + (-t101 * t222 + t366) * qJD(5)) * t219, t24 * t223 + t266 * t127 + (-t101 * t264 - t127 * t324 + t365) * t219, t25 * t218 * t219 + (t219 * t323 - t267) * t99, t25 * t223 + t267 * t127 + (-t127 * t323 + t264 * t99 - t367) * t219, -t127 * t249 - t63 * t223, -t30 * t122 + t145 * t63 - t45 * t99 + t372 * t127 - t244 * t218 + (-t2 + (pkin(9) * t99 + t218 * t30) * qJD(4) + t245 * t222) * t223 + (pkin(9) * t25 + t4 * t218 + t262 * t264 + t30 * t323) * t219, -t45 * t101 - t30 * t123 - t146 * t63 - t373 * t127 - t244 * t222 + (t1 + (pkin(9) * t101 + t222 * t30) * qJD(4) - t245 * t218) * t223 + (-pkin(9) * t24 + t14 * t264 + t4 * t222 - t30 * t324) * t219, t122 * t14 - t123 * t262 + t145 * t24 - t146 * t25 - t373 * t99 - t372 * t101 + t263 * t326 + (-t1 * t218 - t2 * t222 + (-t14 * t222 - t218 * t262) * qJD(5) + t245) * t219, t1 * t146 + t2 * t145 - t30 * t45 - g(1) * (t384 * t74 - t69) - g(2) * (t384 * t72 - t68) - g(3) * (t107 * t384 - t106) + t373 * t14 - t372 * t262 + (t219 * t4 + t30 * t326 - t244) * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t359, -t136 ^ 2 + t138 ^ 2, -t66 - t388, -t359, -t387 - t67, t240, -t70 * t138 - t33 * t202 + t248 - t296, t70 * t136 - t247 - t395, 0, 0, t101 * t286 - t364, (-t24 - t370) * t222 + (-t25 - t360) * t218, -t101 * t138 + t127 * t286 + t367, t127 * t366 - t363, -t127 ^ 2 * t218 + t99 * t138 + t365, -t127 * t138, -pkin(4) * t25 - t15 * t127 + t138 * t262 + t218 * t251 - t222 * t227 - t33 * t99, pkin(4) * t24 - t33 * t101 + t16 * t127 + t14 * t138 + t218 * t227 + t222 * t251, t15 * t101 + t16 * t99 + ((-t25 + t325) * pkin(10) + t400) * t222 + ((qJD(5) * t99 - t24) * pkin(10) + t390) * t218 + t247, t262 * t15 - t14 * t16 - t30 * t33 + t239 * pkin(4) + (qJD(5) * t263 + t1 * t222 - t2 * t218 + t247) * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t371, t101 ^ 2 - t99 ^ 2, -t24 + t370, -t371, t360 - t25, t63, -t30 * t101 - g(1) * (-t218 * t42 + t222 * t74) - g(2) * (-t218 * t40 + t222 * t72) - g(3) * t37 - t390, t30 * t99 - g(1) * (-t218 * t74 - t222 * t42) - g(2) * (-t218 * t72 - t222 * t40) + g(3) * t38 - t400, 0, 0;];
tau_reg = t5;
