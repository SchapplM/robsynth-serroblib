% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% tau_reg [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRRR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:02:33
% EndTime: 2019-03-09 07:02:49
% DurationCPUTime: 6.34s
% Computational Cost: add. (7165->496), mult. (15329->678), div. (0->0), fcn. (11177->18), ass. (0->255)
t241 = sin(qJ(4));
t246 = cos(qJ(4));
t320 = t246 * qJD(3);
t242 = sin(qJ(3));
t334 = qJD(1) * t242;
t176 = t241 * t334 - t320;
t330 = qJD(3) * t241;
t178 = t246 * t334 + t330;
t240 = sin(qJ(5));
t245 = cos(qJ(5));
t113 = t245 * t176 + t178 * t240;
t239 = sin(qJ(6));
t244 = cos(qJ(6));
t276 = t176 * t240 - t245 * t178;
t277 = t113 * t239 + t244 * t276;
t52 = t244 * t113 - t239 * t276;
t395 = t277 * t52;
t388 = t277 ^ 2 - t52 ^ 2;
t247 = cos(qJ(3));
t333 = qJD(1) * t247;
t212 = -qJD(4) + t333;
t207 = -qJD(5) + t212;
t198 = -qJD(6) + t207;
t319 = qJD(1) * qJD(3);
t302 = t247 * t319;
t317 = t242 * qJDD(1);
t326 = qJD(4) * t242;
t389 = -qJD(1) * t326 + qJDD(3);
t102 = qJD(4) * t320 + (t302 + t317) * t246 + t389 * t241;
t103 = ((qJD(4) + t333) * qJD(3) + t317) * t241 - t389 * t246;
t255 = t276 * qJD(5) - t102 * t240 - t245 * t103;
t323 = qJD(5) * t245;
t324 = qJD(5) * t240;
t32 = t245 * t102 - t240 * t103 - t176 * t323 - t178 * t324;
t321 = qJD(6) * t244;
t322 = qJD(6) * t239;
t8 = -t113 * t321 + t239 * t255 + t244 * t32 + t276 * t322;
t386 = -t198 * t52 + t8;
t236 = qJ(4) + qJ(5);
t232 = qJ(6) + t236;
t219 = sin(t232);
t220 = cos(t232);
t233 = qJ(1) + pkin(11);
t224 = cos(t233);
t223 = sin(t233);
t354 = t223 * t247;
t124 = t219 * t224 - t220 * t354;
t353 = t224 * t247;
t126 = t219 * t223 + t220 * t353;
t237 = sin(pkin(11));
t214 = pkin(1) * t237 + pkin(7);
t199 = t214 * qJD(1);
t148 = t242 * qJD(2) + t247 * t199;
t136 = qJD(3) * pkin(8) + t148;
t238 = cos(pkin(11));
t215 = -pkin(1) * t238 - pkin(2);
t168 = -pkin(3) * t247 - pkin(8) * t242 + t215;
t138 = t168 * qJD(1);
t73 = -t136 * t241 + t246 * t138;
t64 = -pkin(9) * t178 + t73;
t57 = -pkin(4) * t212 + t64;
t74 = t136 * t246 + t138 * t241;
t65 = -pkin(9) * t176 + t74;
t63 = t245 * t65;
t26 = t240 * t57 + t63;
t397 = pkin(10) * t113;
t19 = t26 - t397;
t16 = t19 * t322;
t370 = g(3) * t242;
t147 = qJD(2) * t247 - t242 * t199;
t135 = -qJD(3) * pkin(3) - t147;
t108 = pkin(4) * t176 + t135;
t58 = pkin(5) * t113 + t108;
t385 = g(1) * t126 - g(2) * t124 + t220 * t370 + t52 * t58 + t16;
t123 = t219 * t354 + t220 * t224;
t125 = -t219 * t353 + t220 * t223;
t227 = t247 * qJDD(1);
t172 = t242 * t319 + qJDD(4) - t227;
t167 = qJDD(5) + t172;
t192 = t214 * qJDD(1);
t390 = qJD(3) * t147;
t87 = qJDD(3) * pkin(8) + qJDD(2) * t242 + t192 * t247 + t390;
t283 = pkin(3) * t242 - pkin(8) * t247;
t188 = t283 * qJD(3);
t109 = qJD(1) * t188 + qJDD(1) * t168;
t97 = t246 * t109;
t13 = pkin(4) * t172 - pkin(9) * t102 - t74 * qJD(4) - t241 * t87 + t97;
t325 = qJD(4) * t246;
t316 = t241 * t109 + t138 * t325 + t246 * t87;
t327 = qJD(4) * t241;
t268 = -t136 * t327 + t316;
t20 = -pkin(9) * t103 + t268;
t299 = t245 * t13 - t240 * t20;
t257 = -t26 * qJD(5) + t299;
t2 = pkin(5) * t167 - pkin(10) * t32 + t257;
t294 = -t240 * t13 - t245 * t20 - t57 * t323 + t65 * t324;
t3 = pkin(10) * t255 - t294;
t311 = t244 * t2 - t239 * t3;
t400 = -g(1) * t125 + g(2) * t123 + t219 * t370 + t58 * t277 + t311;
t256 = t277 * qJD(6) - t239 * t32 + t244 * t255;
t380 = t198 * t277 + t256;
t309 = t241 * t333;
t372 = pkin(8) + pkin(9);
t310 = qJD(4) * t372;
t185 = t283 * qJD(1);
t340 = t246 * t147 + t241 * t185;
t394 = pkin(9) * t309 - t241 * t310 - t340;
t165 = t246 * t185;
t344 = t246 * t247;
t275 = pkin(4) * t242 - pkin(9) * t344;
t399 = qJD(1) * t275 - t147 * t241 + t246 * t310 + t165;
t282 = g(1) * t224 + g(2) * t223;
t260 = -g(3) * t247 + t282 * t242;
t328 = qJD(3) * t247;
t392 = -qJD(2) * qJD(3) - t192;
t312 = -t199 * t328 + t392 * t242;
t88 = -qJDD(3) * pkin(3) - qJDD(2) * t247 - t312;
t398 = qJD(4) * pkin(8) * t212 + t260 - t88;
t396 = pkin(10) * t276;
t393 = t276 * t113;
t179 = t240 * t241 - t245 * t246;
t267 = t179 * t247;
t375 = qJD(4) + qJD(5);
t342 = qJD(1) * t267 - t375 * t179;
t180 = t240 * t246 + t241 * t245;
t341 = (-t333 + t375) * t180;
t307 = t241 * t328;
t391 = t242 * t325 + t307;
t387 = -t113 ^ 2 + t276 ^ 2;
t384 = -t113 * t207 + t32;
t61 = t240 * t65;
t25 = t245 * t57 - t61;
t18 = t25 + t396;
t14 = -pkin(5) * t207 + t18;
t361 = t244 * t19;
t5 = t239 * t14 + t361;
t383 = -t5 * qJD(6) + t400;
t230 = sin(t236);
t231 = cos(t236);
t351 = t231 * t247;
t131 = -t223 * t351 + t224 * t230;
t133 = t223 * t230 + t224 * t351;
t382 = g(1) * t133 - g(2) * t131 + t108 * t113 + t231 * t370 + t294;
t352 = t230 * t247;
t130 = t223 * t352 + t224 * t231;
t132 = t223 * t231 - t224 * t352;
t381 = -g(1) * t132 + g(2) * t130 + t108 * t276 + t230 * t370 + t257;
t379 = t207 * t276 + t255;
t378 = t399 * t245;
t149 = t180 * t242;
t284 = -t148 + (-t309 + t327) * pkin(4);
t202 = t372 * t241;
t203 = t372 * t246;
t337 = -t240 * t202 + t245 * t203;
t376 = -t202 * t323 - t203 * t324 - t399 * t240 + t394 * t245;
t264 = -t241 * t326 + t247 * t320;
t346 = t242 * t246;
t374 = t172 * t346 - t212 * t264;
t371 = pkin(4) * t240;
t156 = qJDD(6) + t167;
t67 = -qJD(3) * t267 - t375 * t149;
t348 = t241 * t242;
t68 = -t324 * t348 + (t375 * t346 + t307) * t245 + t264 * t240;
t150 = t179 * t242;
t86 = -t149 * t239 - t150 * t244;
t22 = t86 * qJD(6) + t239 * t67 + t244 * t68;
t85 = t244 * t149 - t150 * t239;
t368 = -t85 * t156 + t22 * t198;
t117 = t244 * t179 + t180 * t239;
t367 = -t117 * qJD(6) - t341 * t239 + t342 * t244;
t118 = -t179 * t239 + t180 * t244;
t366 = t118 * qJD(6) + t342 * t239 + t341 * t244;
t365 = t245 * t64 - t61;
t184 = t214 * t344;
t336 = t241 * t168 + t184;
t101 = -pkin(9) * t348 + t336;
t152 = t246 * t168;
t355 = t214 * t241;
t89 = -pkin(9) * t346 + t152 + (-pkin(4) - t355) * t247;
t363 = t245 * t101 + t240 * t89;
t362 = t341 * pkin(5) + t284;
t360 = -t149 * t167 + t68 * t207;
t359 = t102 * t241;
t358 = t176 * t212;
t357 = t178 * t212;
t356 = t212 * t246;
t350 = t239 * t156;
t349 = t240 * t244;
t347 = t241 * t247;
t345 = t244 * t156;
t343 = qJDD(2) - g(3);
t339 = t168 * t325 + t241 * t188;
t329 = qJD(3) * t242;
t338 = t246 * t188 + t329 * t355;
t155 = pkin(4) * t348 + t242 * t214;
t234 = t242 ^ 2;
t335 = -t247 ^ 2 + t234;
t200 = qJD(1) * t215;
t331 = qJD(3) * t176;
t121 = t391 * pkin(4) + t214 * t328;
t222 = -pkin(4) * t246 - pkin(3);
t308 = t212 * t330;
t303 = -t247 * t8 - t277 * t329;
t300 = qJD(6) * t14 + t3;
t43 = t275 * qJD(3) + (-t184 + (pkin(9) * t242 - t168) * t241) * qJD(4) + t338;
t46 = (-t242 * t320 - t247 * t327) * t214 - t391 * pkin(9) + t339;
t297 = -t240 * t46 + t245 * t43;
t296 = -t240 * t64 - t63;
t295 = -t247 * t32 - t276 * t329;
t293 = -t101 * t240 + t245 * t89;
t291 = -qJD(4) * t138 - t87;
t290 = -t102 * t247 + t178 * t329;
t289 = t212 * t214 + t136;
t288 = -t245 * t202 - t203 * t240;
t95 = -pkin(10) * t179 + t337;
t286 = pkin(5) * t334 + t342 * pkin(10) + t337 * qJD(5) + qJD(6) * t95 + t394 * t240 + t378;
t94 = -pkin(10) * t180 + t288;
t285 = -t341 * pkin(10) + qJD(6) * t94 + t376;
t243 = sin(qJ(1));
t248 = cos(qJ(1));
t281 = g(1) * t243 - g(2) * t248;
t21 = -t85 * qJD(6) - t239 * t68 + t244 * t67;
t280 = -t156 * t86 + t198 * t21;
t278 = t150 * t167 + t207 * t67;
t273 = -t247 * t256 - t52 * t329;
t272 = -t113 * t329 - t247 * t255;
t271 = -t172 * t241 + t212 * t325;
t270 = -t101 * t324 + t240 * t43 + t245 * t46 + t89 * t323;
t265 = -qJD(1) * t200 + t282;
t263 = -pkin(8) * t172 - t212 * t135;
t261 = 0.2e1 * t200 * qJD(3) - qJDD(3) * t214;
t45 = pkin(4) * t103 + t88;
t249 = qJD(3) ^ 2;
t254 = g(1) * t223 - g(2) * t224 - 0.2e1 * qJDD(1) * t215 - t214 * t249;
t250 = qJD(1) ^ 2;
t221 = pkin(4) * t245 + pkin(5);
t191 = qJDD(3) * t247 - t242 * t249;
t190 = qJDD(3) * t242 + t247 * t249;
t145 = t223 * t241 + t224 * t344;
t144 = t223 * t246 - t224 * t347;
t143 = -t223 * t344 + t224 * t241;
t142 = t223 * t347 + t224 * t246;
t141 = pkin(5) * t179 + t222;
t104 = pkin(5) * t149 + t155;
t75 = pkin(4) * t178 - pkin(5) * t276;
t40 = pkin(5) * t68 + t121;
t35 = -pkin(10) * t149 + t363;
t34 = -pkin(5) * t247 + pkin(10) * t150 + t293;
t24 = t365 + t396;
t23 = t296 + t397;
t10 = -pkin(5) * t255 + t45;
t7 = -pkin(10) * t68 + t270;
t6 = pkin(5) * t329 - pkin(10) * t67 - qJD(5) * t363 + t297;
t4 = t244 * t14 - t19 * t239;
t1 = [qJDD(1), t281, g(1) * t248 + g(2) * t243 (t281 + (t237 ^ 2 + t238 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), qJDD(1) * t234 + 0.2e1 * t242 * t302, 0.2e1 * t242 * t227 - 0.2e1 * t335 * t319, t190, t191, 0, t242 * t261 + t247 * t254, -t242 * t254 + t247 * t261, t102 * t346 + t178 * t264 (-t176 * t246 - t178 * t241) * t328 + (-t359 - t103 * t246 + (t176 * t241 - t178 * t246) * qJD(4)) * t242, t290 + t374 (t103 + t308) * t247 + (t271 - t331) * t242, -t172 * t247 - t212 * t329 -(-t168 * t327 + t338) * t212 + t152 * t172 - g(1) * t143 - g(2) * t145 + (t214 * t331 - t97 + t289 * t325 + (qJD(3) * t135 - t172 * t214 - t291) * t241) * t247 + (qJD(3) * t73 + t103 * t214 + t135 * t325 + t88 * t241) * t242, t339 * t212 - t336 * t172 - g(1) * t142 - g(2) * t144 + (-t289 * t327 + (t135 * t246 + t178 * t214) * qJD(3) + t316) * t247 + (-t135 * t327 + t214 * t102 + t88 * t246 + (-t214 * t356 - t74) * qJD(3)) * t242, -t150 * t32 - t276 * t67, -t113 * t67 - t149 * t32 - t150 * t255 + t276 * t68, -t278 + t295, t272 + t360, -t167 * t247 - t207 * t329, -t297 * t207 + t293 * t167 - t299 * t247 + t25 * t329 + t121 * t113 - t155 * t255 + t45 * t149 + t108 * t68 - g(1) * t131 - g(2) * t133 + (t207 * t363 + t247 * t26) * qJD(5), -g(1) * t130 - g(2) * t132 + t108 * t67 - t121 * t276 - t45 * t150 + t155 * t32 - t363 * t167 + t270 * t207 - t294 * t247 - t26 * t329, -t21 * t277 + t8 * t86, -t21 * t52 + t22 * t277 + t256 * t86 - t8 * t85, -t280 + t303, t273 + t368, -t156 * t247 - t198 * t329 -(-t239 * t7 + t244 * t6) * t198 + (-t239 * t35 + t244 * t34) * t156 - t311 * t247 + t4 * t329 + t40 * t52 - t104 * t256 + t10 * t85 + t58 * t22 - g(1) * t124 - g(2) * t126 + (-(-t239 * t34 - t244 * t35) * t198 + t5 * t247) * qJD(6), -t5 * t329 - g(1) * t123 - g(2) * t125 + t10 * t86 + t104 * t8 - t16 * t247 + t58 * t21 - t40 * t277 + ((-qJD(6) * t35 + t6) * t198 - t34 * t156 + t2 * t247) * t239 + ((qJD(6) * t34 + t7) * t198 - t35 * t156 + t300 * t247) * t244; 0, 0, 0, t343, 0, 0, 0, 0, 0, t191, -t190, 0, 0, 0, 0, 0 (-t103 + t308) * t247 + (t271 + t331) * t242, t290 - t374, 0, 0, 0, 0, 0, -t272 + t360, t278 + t295, 0, 0, 0, 0, 0, -t273 + t368, t280 + t303; 0, 0, 0, 0, -t242 * t250 * t247, t335 * t250, t317, t227, qJDD(3), qJD(3) * t148 + t265 * t242 + t343 * t247 + t312, t390 + (qJD(3) * t199 - t343) * t242 + (t265 + t392) * t247, -t178 * t356 + t359 (t102 + t358) * t246 + (-t103 + t357) * t241 (-t178 * t242 + t212 * t344) * qJD(1) - t271, t212 * t327 + t172 * t246 + (t176 * t242 - t212 * t347) * qJD(1), t212 * t334, -t73 * t334 - pkin(3) * t103 - t148 * t176 + t165 * t212 + (-t147 * t212 + t263) * t241 + t398 * t246, -pkin(3) * t102 - t148 * t178 - t340 * t212 - t398 * t241 + t263 * t246 + t74 * t334, t180 * t32 - t276 * t342, -t342 * t113 - t179 * t32 + t180 * t255 + t276 * t341, t167 * t180 - t342 * t207 + t276 * t334, t113 * t334 - t167 * t179 + t341 * t207, t207 * t334, t288 * t167 - t222 * t255 + t45 * t179 - t25 * t334 + (t203 * t323 + (-qJD(5) * t202 + t394) * t240 + t378) * t207 + t284 * t113 + t341 * t108 + t260 * t231, t342 * t108 - t337 * t167 + t45 * t180 + t376 * t207 + t222 * t32 - t230 * t260 + t26 * t334 - t276 * t284, t118 * t8 - t277 * t367, -t117 * t8 + t118 * t256 + t277 * t366 - t367 * t52, t118 * t156 - t367 * t198 + t277 * t334, -t117 * t156 + t366 * t198 + t52 * t334, t198 * t334 (-t239 * t95 + t244 * t94) * t156 - t141 * t256 + t10 * t117 - t4 * t334 + t366 * t58 + t362 * t52 + (t239 * t285 + t244 * t286) * t198 + t260 * t220 -(t239 * t94 + t244 * t95) * t156 + t141 * t8 + t10 * t118 + t5 * t334 + t367 * t58 - t362 * t277 + (-t239 * t286 + t244 * t285) * t198 - t260 * t219; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t178 * t176, -t176 ^ 2 + t178 ^ 2, t102 - t358, -t103 - t357, t172, -t136 * t325 - g(1) * t144 + g(2) * t142 - t135 * t178 - t212 * t74 + t97 + (t291 + t370) * t241, g(1) * t145 - g(2) * t143 + g(3) * t346 + t135 * t176 - t212 * t73 - t268, -t393, t387, t384, t379, t167, t296 * t207 + (-t113 * t178 + t167 * t245 + t207 * t324) * pkin(4) + t381, -t365 * t207 + (-t240 * t167 + t178 * t276 + t207 * t323) * pkin(4) + t382, -t395, t388, t386, t380, t156, t221 * t345 + (t23 * t244 - t239 * t24) * t198 - t75 * t52 + (-t240 * t350 - (-t239 * t245 - t349) * t198 * qJD(5)) * pkin(4) + (-(-pkin(4) * t349 - t221 * t239) * t198 - t5) * qJD(6) + t400, t75 * t277 + (-t221 * t156 - t2 + (-t23 + (-qJD(5) - qJD(6)) * t371) * t198) * t239 + (-t156 * t371 + (pkin(4) * t323 + qJD(6) * t221 - t24) * t198 - t300) * t244 + t385; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t393, t387, t384, t379, t167, -t207 * t26 + t381, -t207 * t25 + t382, -t395, t388, t386, t380, t156 (-t18 * t239 - t361) * t198 + (t198 * t322 + t276 * t52 + t345) * pkin(5) + t383 (t19 * t198 - t2) * t239 + (-t18 * t198 - t300) * t244 + (t198 * t321 - t276 * t277 - t350) * pkin(5) + t385; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t395, t388, t386, t380, t156, -t198 * t5 + t383, -t198 * t4 - t239 * t2 - t244 * t300 + t385;];
tau_reg  = t1;
