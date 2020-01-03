% Calculate minimal parameter regressor of coriolis matrix for
% S5RRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x28]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPRR15_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR15_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:43:12
% EndTime: 2019-12-31 20:43:21
% DurationCPUTime: 3.86s
% Computational Cost: add. (2696->332), mult. (5763->479), div. (0->0), fcn. (5627->6), ass. (0->286)
t281 = qJD(4) + qJD(5);
t216 = cos(qJ(5));
t217 = cos(qJ(4));
t315 = t216 * t217;
t213 = sin(qJ(5));
t214 = sin(qJ(4));
t322 = t213 * t214;
t155 = -t315 + t322;
t147 = t155 * qJD(5);
t215 = sin(qJ(2));
t262 = -t322 / 0.2e1;
t182 = t215 * t315;
t362 = t182 / 0.2e1;
t365 = -t155 / 0.2e1;
t226 = t362 + (t262 + t365) * t215;
t371 = t226 * qJD(1);
t380 = t155 * qJD(4) + t147 - t371;
t219 = -pkin(2) - pkin(7);
t378 = -pkin(8) + t219;
t164 = t378 * t214;
t165 = t378 * t217;
t118 = t216 * t164 + t213 * t165;
t379 = t281 * t118;
t364 = t155 / 0.2e1;
t372 = t262 + t364;
t225 = t372 * t215 + t362;
t377 = t225 * t281;
t376 = t226 * t281;
t316 = t216 * t214;
t321 = t213 * t217;
t156 = t316 + t321;
t375 = t281 * t156;
t218 = cos(qJ(2));
t358 = pkin(4) * t217;
t366 = pkin(3) + pkin(6);
t259 = t358 + t366;
t246 = t259 * t218;
t233 = -t246 / 0.2e1;
t314 = t217 * t218;
t183 = t216 * t314;
t317 = t214 * t218;
t131 = t213 * t317 - t183;
t359 = pkin(4) * t214;
t199 = qJ(3) + t359;
t324 = t199 * t131;
t374 = t156 * t233 + t324 / 0.2e1;
t132 = t156 * t218;
t323 = t199 * t132;
t373 = t155 * t233 - t323 / 0.2e1;
t336 = qJ(3) * t215;
t248 = -t218 * t219 + t336;
t209 = t214 ^ 2;
t211 = t217 ^ 2;
t188 = t209 - t211;
t284 = t218 * qJD(1);
t270 = t217 * t284;
t228 = qJD(2) * t188 + 0.2e1 * t214 * t270;
t283 = t218 * qJD(3);
t370 = qJD(2) * t248 - t283;
t247 = t213 * t164 - t216 * t165;
t369 = t281 * t247;
t356 = t215 * pkin(4);
t149 = -pkin(1) - t248;
t178 = t366 * t215;
t161 = t217 * t178;
t111 = t149 * t214 - t161;
t83 = pkin(8) * t317 - t111;
t67 = t83 + t356;
t56 = t216 * t67;
t368 = -t56 / 0.2e1;
t367 = -t67 / 0.2e1;
t363 = t156 / 0.2e1;
t361 = t183 / 0.2e1;
t360 = -t215 / 0.2e1;
t357 = pkin(8) * t215;
t355 = t218 * pkin(4);
t31 = -t131 * t156 - t132 * t155;
t354 = t281 * t31;
t44 = t131 * t365 + t132 * t363;
t353 = t281 * t44;
t240 = t316 / 0.2e1 + t321 / 0.2e1;
t96 = t156 * t360 + t240 * t215;
t352 = t281 * t96;
t92 = (t363 + t240) * t215;
t351 = t281 * t92;
t350 = pkin(4) * qJD(4);
t349 = pkin(4) * qJD(5);
t276 = -t356 / 0.2e1;
t254 = t276 + t83 / 0.2e1;
t7 = (t367 + t254) * t213;
t348 = qJD(1) * t7;
t9 = t254 * t216 + t368;
t347 = qJD(1) * t9;
t318 = t214 * t215;
t130 = t213 * t318 - t182;
t141 = t259 * t215;
t319 = t214 * t178;
t112 = t149 * t217 + t319;
t84 = -pkin(8) * t314 + t112;
t342 = t213 * t84;
t34 = -t56 + t342;
t313 = t218 * qJ(3);
t177 = t215 * pkin(2) - t313;
t151 = pkin(7) * t215 + t177;
t179 = t366 * t218;
t163 = t179 * t217;
t70 = t355 + t163 + (-t151 - t357) * t214;
t340 = t216 * t70;
t139 = t217 * t151;
t162 = t179 * t214;
t309 = t139 + t162;
t85 = t217 * t357 + t309;
t341 = t213 * t85;
t1 = (t340 - t341) * t215 + t141 * t131 + (t259 * t130 - t34) * t218;
t346 = t1 * qJD(1);
t133 = t156 * t215;
t337 = t216 * t85;
t344 = t213 * t70;
t338 = t216 * t84;
t35 = t213 * t67 + t338;
t2 = (t337 + t344) * t215 + t35 * t218 - t141 * t132 - t133 * t246;
t345 = t2 * qJD(1);
t343 = t213 * t83;
t339 = t216 * t83;
t36 = -t338 - t343;
t17 = t36 * t215 + (t131 * t359 - t259 * t132) * t218;
t335 = qJD(1) * t17;
t101 = t131 * t246;
t37 = t339 - t342;
t18 = -pkin(4) * t132 * t317 + t215 * t37 - t101;
t334 = qJD(1) * t18;
t21 = -t215 * t34 - t101;
t333 = qJD(1) * t21;
t22 = -t132 * t246 - t35 * t215;
t332 = qJD(1) * t22;
t48 = -t111 * t215 + t179 * t314;
t331 = qJD(1) * t48;
t49 = -t112 * t215 - t179 * t317;
t330 = qJD(1) * t49;
t59 = -t130 * t215 + t131 * t218;
t329 = qJD(1) * t59;
t60 = -t132 * t218 + t133 * t215;
t328 = qJD(1) * t60;
t327 = t247 * t215;
t326 = t118 * t215;
t325 = t156 * t199;
t320 = t214 * t151;
t27 = (-t111 - t161) * t218 - t320 * t215;
t312 = t27 * qJD(1);
t28 = t309 * t215 - t179 * t318 + (t112 - t319) * t218;
t311 = t28 * qJD(1);
t38 = t130 * t132 + t131 * t133;
t310 = t38 * qJD(1);
t76 = t92 * qJD(1);
t210 = t215 ^ 2;
t212 = t218 ^ 2;
t189 = t212 - t210;
t307 = qJD(1) * t132;
t158 = t189 * t214;
t306 = qJD(1) * t158;
t160 = t189 * t217;
t305 = qJD(1) * t160;
t304 = qJD(1) * t215;
t303 = qJD(1) * t217;
t302 = qJD(2) * qJ(3);
t301 = qJD(2) * t155;
t300 = qJD(2) * t156;
t299 = qJD(2) * t199;
t298 = qJD(3) * t214;
t297 = qJD(3) * t215;
t296 = qJD(3) * t217;
t295 = qJD(4) * t214;
t294 = qJD(4) * t215;
t293 = qJD(4) * t217;
t292 = qJD(4) * t219;
t252 = -pkin(2) * t218 - t336;
t166 = -pkin(1) + t252;
t119 = t166 * t218 + t177 * t215;
t291 = t119 * qJD(1);
t120 = -t166 * t215 + t177 * t218;
t290 = t120 * qJD(1);
t289 = t189 * qJD(1);
t288 = t210 * qJD(1);
t287 = t214 * qJD(2);
t286 = t215 * qJD(2);
t285 = t217 * qJD(2);
t206 = t218 * qJD(2);
t282 = t218 * qJD(4);
t280 = pkin(1) * t304;
t279 = pkin(1) * t284;
t278 = pkin(6) * t286;
t277 = t358 / 0.2e1;
t275 = t355 / 0.2e1;
t274 = -t341 / 0.2e1;
t273 = -t337 / 0.2e1;
t272 = t130 * t304;
t271 = t133 * t304;
t269 = t214 * t285;
t268 = t214 * t206;
t267 = t215 * t282;
t266 = t166 * t177 * qJD(1);
t265 = t166 * t304;
t194 = t215 * t206;
t193 = t215 * t284;
t196 = t217 * t206;
t264 = t214 * t293;
t263 = t326 / 0.2e1;
t261 = -t317 / 0.2e1;
t260 = pkin(4) * t281;
t258 = t281 * t215;
t257 = qJD(4) + t304;
t256 = t214 * t196;
t253 = t275 + t70 / 0.2e1;
t220 = t132 * t277 + t247 * t360;
t242 = -t344 / 0.2e1 + t273;
t224 = -t324 / 0.2e1 + t242;
t3 = (-t213 * pkin(4) / 0.2e1 + t359 * t365 + t259 * t363) * t218 + t220 + t224;
t87 = -t155 * t358 - t325;
t251 = t3 * qJD(1) - t87 * qJD(2);
t221 = t131 * t277 + t263;
t241 = t274 + t340 / 0.2e1;
t223 = t323 / 0.2e1 + t241;
t4 = (t216 * pkin(4) / 0.2e1 + t359 * t363 + t259 * t364) * t218 + t221 + t223;
t86 = -t155 * t199 + t156 * t358;
t250 = t4 * qJD(1) - t86 * qJD(2);
t249 = qJD(5) + t257;
t47 = t131 ^ 2 - t132 ^ 2;
t15 = qJD(1) * t47 + qJD(2) * t31;
t71 = -t155 ^ 2 + t156 ^ 2;
t19 = qJD(1) * t31 + qJD(2) * t71;
t93 = (t363 - t240) * t218;
t245 = qJD(1) * t93 - t301;
t91 = t372 * t218 + t361;
t244 = qJD(1) * t91 + t300;
t243 = -t288 - t294;
t239 = t219 * t360 - t313 / 0.2e1;
t227 = t239 * t217;
t99 = t139 / 0.2e1 + t227;
t238 = qJ(3) * t287 - t99 * qJD(1);
t98 = (t151 / 0.2e1 + t239) * t214;
t237 = -qJ(3) * t285 - t98 * qJD(1);
t232 = t246 / 0.2e1;
t11 = t155 * t232 + t223 + t263;
t236 = qJD(1) * t11 + t155 * t299;
t12 = -t327 / 0.2e1 + t156 * t232 + t224;
t235 = qJD(1) * t12 + t156 * t299;
t24 = -qJD(2) * t44 + t131 * t307;
t32 = qJD(1) * t44 + t155 * t300;
t234 = t257 * t317;
t150 = (t211 / 0.2e1 - t209 / 0.2e1) * t218;
t231 = qJD(1) * t150 + t269;
t230 = t212 * t214 * t303 - qJD(2) * t150;
t159 = t188 * t212;
t229 = -qJD(1) * t159 + 0.2e1 * t256;
t222 = t252 * qJD(2) + t283;
t205 = pkin(6) * t206;
t203 = -t284 / 0.2e1;
t202 = t284 / 0.2e1;
t201 = t206 / 0.2e1;
t195 = t215 * t303;
t192 = t214 * t304;
t154 = -t195 - t293;
t153 = -t192 - t295;
t152 = t193 + t282 / 0.2e1;
t140 = t150 * qJD(4);
t127 = t193 + (qJD(4) / 0.2e1 + qJD(5) / 0.2e1) * t218;
t94 = -t132 / 0.2e1 - t240 * t218;
t90 = t213 * t261 + t218 * t365 + t361;
t82 = t225 * qJD(3);
t79 = t96 * qJD(3);
t75 = t92 * qJD(3);
t73 = t226 * qJD(3);
t51 = -t162 - t139 / 0.2e1 + t227;
t50 = t163 - t320 / 0.2e1 + t239 * t214;
t46 = qJD(2) * t226 - t132 * t304;
t45 = qJD(2) * t92 - t131 * t304;
t40 = -t76 - t375;
t26 = qJD(2) * t225 + t249 * t132;
t25 = qJD(2) * t96 + t249 * t131;
t14 = -t326 / 0.2e1 + t241 + t373;
t13 = t327 / 0.2e1 + t242 + t374;
t10 = t216 * t276 + t342 + t368 - t339 / 0.2e1;
t8 = -t338 - t343 / 0.2e1 + (t276 + t367) * t213;
t6 = t214 * t155 * t275 - t253 * t213 - t220 + t273 + t374;
t5 = pkin(4) * t156 * t261 + t253 * t216 - t221 + t274 + t373;
t16 = [0, 0, 0, t194, t189 * qJD(2), 0, 0, 0, -pkin(1) * t286, -pkin(1) * t206, 0, qJD(2) * t120 - t215 * t283, -qJD(2) * t119 + qJD(3) * t210, (qJD(2) * t177 - t297) * t166, -t194 * t209 + t212 * t264, -qJD(4) * t159 - 0.2e1 * t215 * t256, -qJD(2) * t158 - t217 * t267, -qJD(2) * t160 + t214 * t267, t194, qJD(2) * t27 + qJD(4) * t49 + t210 * t298, -qJD(2) * t28 - qJD(4) * t48 + t210 * t296, (-qJD(2) * t133 - t131 * t281) * t132, qJD(2) * t38 + t281 * t47, qJD(2) * t60 + t131 * t258, qJD(2) * t59 + t132 * t258, t194, qJD(2) * t1 + qJD(4) * t17 + qJD(5) * t22 + t133 * t297, -qJD(2) * t2 - qJD(4) * t18 - qJD(5) * t21 - t130 * t297; 0, 0, 0, t193, t289, t206, -t286, 0, -t205 - t280, t278 - t279, t222, t205 + t290, -t278 - t291, pkin(6) * t222 + t266, -t140 + (-t209 * t284 + t269) * t215, -t228 * t215 + 0.2e1 * t218 * t264, t196 - t306, -t268 - t305, t152, t50 * qJD(4) - t178 * t287 - t370 * t217 + t312, t51 * qJD(4) - t178 * t285 + t370 * t214 - t311, (-t301 - t307) * t133 + t353, t310 + (t130 * t155 - t133 * t156) * qJD(2) + t354, -t155 * t206 + t328 + t352, -t156 * t206 + t329 + t377, t127, t346 + (t130 * t199 - t141 * t156 - t218 * t247) * qJD(2) + t90 * qJD(3) + t5 * qJD(4) + t14 * qJD(5), -t345 + (-t118 * t218 + t133 * t199 + t141 * t155) * qJD(2) + t94 * qJD(3) + t6 * qJD(4) + t13 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t206, -t193, t288, t205 - t265, 0, 0, 0, 0, 0, t214 * t288 + t196, t217 * t288 - t268, 0, 0, 0, 0, 0, qJD(2) * t90 + t271 + t352, qJD(2) * t94 - t272 + t377; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t230, t229, -t257 * t314, t234, t201, qJD(2) * t50 - qJD(4) * t112 + t330, qJD(2) * t51 + qJD(4) * t111 - t331, -t24, t15, t25, t26, t201, qJD(2) * t5 + qJD(4) * t36 + qJD(5) * t8 + t335 + t79, qJD(2) * t6 - qJD(4) * t37 + qJD(5) * t10 - t334 + t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, t15, t25, t26, t201, qJD(2) * t14 + qJD(4) * t8 - qJD(5) * t35 + t332 + t79, qJD(2) * t13 + qJD(4) * t10 + qJD(5) * t34 - t333 + t82; 0, 0, 0, -t193, -t289, 0, 0, 0, t280, t279, 0, -t290, t291, -t266, t193 * t209 - t140, 0.2e1 * t217 * t234, -t214 * t294 + t306, -t215 * t293 + t305, -t152, qJD(4) * t98 - t312, qJD(4) * t99 + t311, t133 * t307 + t353, -t310 + t354, -t328 - t351, -t329 - t376, -t127, qJD(3) * t91 - qJD(4) * t4 - qJD(5) * t11 - t346, qJD(3) * t93 - qJD(4) * t3 - qJD(5) * t12 + t345; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), qJ(3) * qJD(3), -t264, t188 * qJD(4), 0, 0, 0, qJ(3) * t293 + t298, -qJ(3) * t295 + t296, t155 * t375, t281 * t71, 0, 0, 0, qJD(3) * t156 + qJD(4) * t86 - t147 * t199, -qJD(3) * t155 + qJD(4) * t87 - qJD(5) * t325; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t302, 0, 0, 0, 0, 0, t287, t285, 0, 0, 0, 0, 0, t244, t245; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t231, t228, t153, t154, t203, -t214 * t292 - t237, -t217 * t292 - t238, t32, t19, t40, t380, t203, -t250 - t379, -t251 + t369; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t19, t40, t380, t203, -t236 - t379, -t235 + t369; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, -t288, t265, 0, 0, 0, 0, 0, t243 * t214, t243 * t217, 0, 0, 0, 0, 0, -qJD(2) * t91 - t271 - t351, -qJD(2) * t93 + t272 - t376; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), -t302, 0, 0, 0, 0, 0, -t287, -t285, 0, 0, 0, 0, 0, -t244, -t245; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, t154, 0, 0, 0, 0, 0, t40, t380; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t380; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t230, -t229, (t270 + t287) * t215, (-t214 * t284 + t285) * t215, t201, -qJD(2) * t98 + t214 * t297 - t330, -qJD(2) * t99 + t215 * t296 + t331, t24, -t15, t45, t46, t201, qJD(2) * t4 + qJD(5) * t7 - t335 + t75, qJD(2) * t3 + qJD(5) * t9 + t334 + t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t231, -t228, t192, t195, t202, t237, t238, -t32, -t19, t76, t371, t202, t250, t251; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t192, t195, 0, 0, 0, 0, 0, t76, t371; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t213 * t349, -t216 * t349; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t213 * t260 + t348, -t216 * t260 + t347; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t15, t45, t46, t201, qJD(2) * t11 - qJD(4) * t7 - t332 + t75, qJD(2) * t12 - qJD(4) * t9 + t333 + t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t19, t76, t371, t202, t236, t235; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, t371; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t213 * t350 - t348, t216 * t350 - t347; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t16;
