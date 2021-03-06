% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRRR9_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR9_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:25:19
% EndTime: 2019-03-09 07:25:41
% DurationCPUTime: 8.82s
% Computational Cost: add. (14458->565), mult. (31031->790), div. (0->0), fcn. (20914->8), ass. (0->246)
t252 = sin(qJ(3));
t255 = cos(qJ(3));
t273 = pkin(3) * t255 + pkin(8) * t252;
t214 = t273 * qJD(1);
t256 = -pkin(1) - pkin(7);
t234 = qJD(1) * t256 + qJD(2);
t254 = cos(qJ(4));
t251 = sin(qJ(4));
t328 = t251 * t255;
t147 = t254 * t214 - t234 * t328;
t361 = -pkin(9) - pkin(8);
t292 = qJD(4) * t361;
t327 = t252 * t254;
t298 = pkin(9) * t327;
t388 = (pkin(4) * t255 + t298) * qJD(1) + t147 - t254 * t292;
t324 = t254 * t255;
t148 = t251 * t214 + t234 * t324;
t315 = qJD(1) * t252;
t290 = t251 * t315;
t387 = pkin(9) * t290 - t251 * t292 + t148;
t302 = t254 * qJD(3);
t314 = qJD(1) * t255;
t205 = t251 * t314 - t302;
t289 = t254 * t314;
t313 = qJD(3) * t251;
t207 = t289 + t313;
t250 = sin(qJ(5));
t253 = cos(qJ(5));
t137 = t205 * t253 + t207 * t250;
t249 = sin(qJ(6));
t270 = -t205 * t250 + t253 * t207;
t360 = cos(qJ(6));
t373 = -t249 * t137 + t270 * t360;
t72 = t360 * t137 + t249 * t270;
t359 = t72 * t373;
t300 = qJD(4) + qJD(5);
t304 = qJD(5) * t253;
t307 = qJD(4) * t254;
t325 = t253 * t254;
t329 = t250 * t251;
t319 = t250 * t290 - t253 * t307 - t254 * t304 + t300 * t329 - t315 * t325;
t210 = t250 * t254 + t251 * t253;
t146 = t300 * t210;
t188 = t210 * qJD(1);
t318 = t252 * t188 + t146;
t379 = t373 ^ 2 - t72 ^ 2;
t238 = qJD(4) + t315;
t232 = qJD(5) + t238;
t370 = pkin(10) * t270;
t218 = pkin(3) * t252 - pkin(8) * t255 + qJ(2);
t190 = t218 * qJD(1);
t217 = t252 * t234;
t196 = qJD(3) * pkin(8) + t217;
t122 = t254 * t190 - t196 * t251;
t100 = -pkin(9) * t207 + t122;
t91 = pkin(4) * t238 + t100;
t123 = t190 * t251 + t196 * t254;
t101 = -pkin(9) * t205 + t123;
t95 = t250 * t101;
t39 = t253 * t91 - t95;
t34 = t39 - t370;
t32 = pkin(5) * t232 + t34;
t383 = pkin(10) * t137;
t97 = t253 * t101;
t40 = t250 * t91 + t97;
t35 = t40 - t383;
t353 = t249 * t35;
t10 = t32 * t360 - t353;
t295 = t360 * t35;
t11 = t249 * t32 + t295;
t386 = -t10 * t72 + t11 * t373;
t306 = qJD(4) * t255;
t284 = t251 * t306;
t286 = t252 * t302;
t264 = t284 + t286;
t156 = qJD(1) * t264 - qJD(4) * t302;
t312 = qJD(3) * t252;
t287 = t251 * t312;
t309 = qJD(4) * t207;
t157 = -qJD(1) * t287 + t309;
t305 = qJD(5) * t250;
t267 = -t250 * t156 + t157 * t253 - t205 * t305 + t207 * t304;
t282 = qJD(6) * t360;
t303 = qJD(6) * t249;
t56 = t253 * t156 + t250 * t157 + t205 * t304 + t207 * t305;
t20 = t137 * t282 + t249 * t267 + t270 * t303 + t360 * t56;
t226 = qJD(6) + t232;
t377 = t226 * t72 - t20;
t301 = qJD(1) * qJD(3);
t240 = t255 * t301;
t203 = qJD(3) * t273 + qJD(2);
t172 = t203 * qJD(1);
t262 = -qJD(4) * t123 + t254 * t172;
t311 = qJD(3) * t255;
t47 = pkin(9) * t156 + (pkin(4) * qJD(1) - t234 * t251) * t311 + t262;
t288 = t234 * t311;
t308 = qJD(4) * t251;
t65 = t251 * t172 + t190 * t307 - t196 * t308 + t254 * t288;
t49 = -pkin(9) * t157 + t65;
t9 = -qJD(5) * t40 - t250 * t49 + t253 * t47;
t6 = pkin(5) * t240 + pkin(10) * t56 + t9;
t279 = t101 * t305 - t250 * t47 - t253 * t49 - t91 * t304;
t7 = -pkin(10) * t267 - t279;
t261 = -t249 * t6 - t282 * t32 + t35 * t303 - t360 * t7;
t334 = t234 * t255;
t354 = qJD(3) * pkin(3);
t197 = -t334 - t354;
t150 = pkin(4) * t205 + t197;
t82 = pkin(5) * t137 + t150;
t376 = t82 * t72 + t261;
t227 = t361 * t251;
t228 = t361 * t254;
t356 = -t227 * t304 - t228 * t305 + t388 * t250 + t387 * t253;
t155 = t250 * t227 - t253 * t228;
t355 = -t155 * qJD(5) + t387 * t250 - t388 * t253;
t2 = -qJD(6) * t11 - t249 * t7 + t360 * t6;
t365 = -t373 * t82 + t2;
t21 = qJD(6) * t373 - t249 * t56 + t360 * t267;
t363 = t226 * t373 - t21;
t382 = -pkin(5) * t314 + t319 * pkin(10) + t355;
t381 = t318 * pkin(10) + t356;
t173 = t210 * t252;
t209 = -t325 + t329;
t176 = t209 * t255;
t346 = qJD(3) * t176 + t300 * t173 + t188;
t175 = t209 * t252;
t345 = t209 * qJD(1) + t300 * t175 - t210 * t311;
t342 = t137 * t270;
t283 = t254 * t306;
t380 = t283 - t287;
t378 = -t137 ^ 2 + t270 ^ 2;
t375 = t137 * t232 - t56;
t374 = t137 * t150 + t279;
t299 = 0.2e1 * qJD(1);
t369 = -t122 * t238 + t65;
t66 = -t251 * t288 + t262;
t368 = -t123 * t238 - t66;
t202 = t254 * t218;
t281 = -t251 * t256 + pkin(4);
t132 = -pkin(9) * t324 + t252 * t281 + t202;
t326 = t252 * t256;
t231 = t254 * t326;
t161 = t251 * t218 + t231;
t149 = -pkin(9) * t328 + t161;
t79 = t250 * t132 + t253 * t149;
t274 = -t217 + (t290 + t308) * pkin(4);
t364 = -t270 * t150 + t9;
t362 = t232 * t270 - t267;
t154 = t253 * t227 + t228 * t250;
t114 = -pkin(10) * t210 + t154;
t115 = -pkin(10) * t209 + t155;
t62 = t114 * t360 - t249 * t115;
t358 = qJD(6) * t62 + t382 * t249 - t381 * t360;
t63 = t249 * t114 + t115 * t360;
t357 = -qJD(6) * t63 + t381 * t249 + t382 * t360;
t46 = t253 * t100 - t95;
t144 = -t249 * t209 + t210 * t360;
t352 = -qJD(6) * t144 + t319 * t249 - t318 * t360;
t351 = t209 * t282 + t210 * t303 + t318 * t249 + t319 * t360;
t110 = -t249 * t173 - t175 * t360;
t350 = qJD(6) * t110 - t346 * t249 - t345 * t360;
t108 = -t173 * t360 + t249 * t175;
t349 = -qJD(6) * t108 - t345 * t249 + t346 * t360;
t243 = pkin(4) * t253 + pkin(5);
t330 = t249 * t250;
t45 = -t100 * t250 - t97;
t36 = t45 + t383;
t37 = t46 - t370;
t348 = -t249 * t36 - t360 * t37 + t243 * t282 + (-t250 * t303 + (t253 * t360 - t330) * qJD(5)) * pkin(4);
t291 = t360 * t250;
t347 = t249 * t37 - t360 * t36 - t243 * t303 + (-t250 * t282 + (-t249 * t253 - t291) * qJD(5)) * pkin(4);
t341 = t156 * t251;
t340 = t156 * t252;
t339 = t157 * t252;
t338 = t157 * t254;
t337 = t205 * t238;
t336 = t207 * t205;
t335 = t207 * t238;
t333 = t238 * t251;
t332 = t238 * t252;
t331 = t238 * t254;
t257 = qJD(3) ^ 2;
t323 = t257 * t252;
t322 = t257 * t255;
t258 = qJD(1) ^ 2;
t321 = t258 * qJ(2);
t320 = t318 * pkin(5) + t274;
t248 = t255 ^ 2;
t317 = t252 ^ 2 - t248;
t316 = -t257 - t258;
t310 = qJD(3) * t256;
t296 = qJD(2) * t299;
t294 = t251 * t326;
t293 = t255 * t258 * t252;
t244 = -pkin(4) * t254 - pkin(3);
t285 = t255 * t310;
t116 = t157 * pkin(4) + t234 * t312;
t78 = t253 * t132 - t149 * t250;
t278 = -t197 + t334;
t204 = pkin(4) * t328 - t255 * t256;
t277 = t205 + t302;
t276 = -t207 + t313;
t275 = qJD(4) * t252 + qJD(1);
t272 = t122 * t254 + t123 * t251;
t271 = t122 * t251 - t123 * t254;
t269 = t278 * qJD(3);
t268 = qJD(1) * t248 - t332;
t52 = pkin(5) * t252 + pkin(10) * t176 + t78;
t174 = t210 * t255;
t59 = -pkin(10) * t174 + t79;
t26 = -t249 * t59 + t360 * t52;
t27 = t249 * t52 + t360 * t59;
t266 = -pkin(8) * t311 + t197 * t252;
t111 = -t249 * t174 - t176 * t360;
t158 = t380 * pkin(4) + t252 * t310;
t184 = t254 * t203;
t69 = t184 + (-t231 + (pkin(9) * t255 - t218) * t251) * qJD(4) + (t255 * t281 + t298) * qJD(3);
t98 = -qJD(4) * t294 + t251 * t203 + t218 * t307 + t254 * t285;
t80 = -pkin(9) * t380 + t98;
t24 = t132 * t304 - t149 * t305 + t250 * t69 + t253 * t80;
t25 = -t79 * qJD(5) - t250 * t80 + t253 * t69;
t259 = -qJD(4) * t272 - t251 * t66 + t254 * t65;
t245 = qJ(2) * t296;
t230 = t252 * t240;
t182 = pkin(4) * t291 + t249 * t243;
t181 = -pkin(4) * t330 + t243 * t360;
t169 = pkin(5) * t209 + t244;
t160 = t202 - t294;
t143 = t209 * t360 + t210 * t249;
t142 = pkin(5) * t174 + t204;
t109 = t174 * t360 - t176 * t249;
t102 = pkin(4) * t207 + pkin(5) * t270;
t99 = -t161 * qJD(4) - t251 * t285 + t184;
t90 = -t305 * t328 + (t300 * t324 - t287) * t253 - t264 * t250;
t88 = t146 * t255 - t250 * t287 + t253 * t286;
t64 = pkin(5) * t90 + t158;
t38 = pkin(5) * t267 + t116;
t31 = qJD(6) * t111 - t249 * t88 + t360 * t90;
t29 = t174 * t282 - t176 * t303 + t249 * t90 + t360 * t88;
t17 = -pkin(10) * t90 + t24;
t16 = pkin(5) * t311 + pkin(10) * t88 + t25;
t13 = t34 * t360 - t353;
t12 = -t249 * t34 - t295;
t4 = -qJD(6) * t27 + t16 * t360 - t249 * t17;
t3 = qJD(6) * t26 + t249 * t16 + t17 * t360;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t296, t245, -0.2e1 * t230, 0.2e1 * t317 * t301, -t323, 0.2e1 * t230, -t322, 0, -t256 * t323 + (qJ(2) * t311 + qJD(2) * t252) * t299, -t256 * t322 + (-qJ(2) * t312 + qJD(2) * t255) * t299, 0, t245, -t156 * t324 - t207 * t264 (t205 * t254 + t207 * t251) * t312 + (t341 - t338 + (t205 * t251 - t207 * t254) * qJD(4)) * t255, -t238 * t284 - t340 + (t207 * t255 + t254 * t268) * qJD(3), t157 * t328 + t205 * t380, -t238 * t283 - t339 + (-t205 * t255 - t251 * t268) * qJD(3), t238 * t311 + t230, t238 * t99 + t252 * t66 + (-t157 * t256 + t197 * t307) * t255 + ((qJD(1) * t160 + t122) * t255 + (t205 * t256 + t251 * t278) * t252) * qJD(3), -t238 * t98 - t252 * t65 + (t156 * t256 - t197 * t308) * t255 + ((-qJD(1) * t161 - t123) * t255 + (t207 * t256 + t254 * t278) * t252) * qJD(3), t156 * t160 - t157 * t161 - t205 * t98 - t207 * t99 + t272 * t312 + (qJD(4) * t271 - t251 * t65 - t254 * t66) * t255, t122 * t99 + t123 * t98 + t160 * t66 + t161 * t65 - t269 * t326, t176 * t56 - t270 * t88, t88 * t137 + t56 * t174 + t176 * t267 - t270 * t90, -t232 * t88 - t252 * t56 + (-qJD(1) * t176 + t270) * t311, t137 * t90 + t174 * t267, -t90 * t232 - t267 * t252 + (-qJD(1) * t174 - t137) * t311, t232 * t311 + t230, t25 * t232 + t9 * t252 + t158 * t137 + t204 * t267 + t116 * t174 + t150 * t90 + (qJD(1) * t78 + t39) * t311, -t116 * t176 + t270 * t158 - t150 * t88 - t204 * t56 - t232 * t24 + t252 * t279 + (-qJD(1) * t79 - t40) * t311, -t24 * t137 + t174 * t279 + t9 * t176 - t25 * t270 - t267 * t79 + t39 * t88 - t40 * t90 + t78 * t56, t116 * t204 + t150 * t158 + t24 * t40 + t25 * t39 - t279 * t79 + t78 * t9, -t111 * t20 - t29 * t373, t109 * t20 - t111 * t21 + t29 * t72 - t31 * t373, -t20 * t252 - t226 * t29 + (qJD(1) * t111 + t373) * t311, t109 * t21 + t31 * t72, -t21 * t252 - t226 * t31 + (-qJD(1) * t109 - t72) * t311, t226 * t311 + t230, t109 * t38 + t142 * t21 + t2 * t252 + t226 * t4 + t31 * t82 + t64 * t72 + (qJD(1) * t26 + t10) * t311, t261 * t252 + t111 * t38 - t142 * t20 - t226 * t3 - t29 * t82 + t64 * t373 + (-qJD(1) * t27 - t11) * t311, t10 * t29 + t109 * t261 - t11 * t31 - t111 * t2 + t20 * t26 - t21 * t27 - t3 * t72 - t373 * t4, t10 * t4 + t11 * t3 + t142 * t38 + t2 * t26 - t261 * t27 + t64 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t258, -t321, 0, 0, 0, 0, 0, 0, t316 * t252, t316 * t255, 0, -t321, 0, 0, 0, 0, 0, 0, -t157 * t255 - t275 * t331 + (t205 * t252 + (-t238 - t315) * t328) * qJD(3), t156 * t255 + t275 * t333 + (-t238 * t324 + (t207 - t289) * t252) * qJD(3) (-t205 * t311 + t207 * t275 - t339) * t254 + (t205 * t275 + t207 * t311 - t340) * t251, -t271 * t311 - t272 * qJD(1) + (-t269 + t259) * t252, 0, 0, 0, 0, 0, 0, -t255 * t267 + t345 * t232 + (t137 * t252 - t173 * t314) * qJD(3), t255 * t56 + t346 * t232 + (t175 * t314 + t252 * t270) * qJD(3), t137 * t346 - t173 * t56 + t175 * t267 - t270 * t345, -t116 * t255 + t150 * t312 - t173 * t9 + t175 * t279 + t345 * t39 - t346 * t40, 0, 0, 0, 0, 0, 0, -t21 * t255 - t350 * t226 + (t108 * t314 + t252 * t72) * qJD(3), t20 * t255 + t349 * t226 + (-t110 * t314 + t252 * t373) * qJD(3), t108 * t20 - t110 * t21 + t349 * t72 + t350 * t373, -t10 * t350 + t108 * t2 - t11 * t349 - t110 * t261 - t255 * t38 + t312 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t293, -t317 * t258, 0, -t293, 0, 0, -t255 * t321, t252 * t321, 0, 0, t207 * t331 - t341 (-t156 - t337) * t254 + (-t157 - t335) * t251, t238 * t307 + (t238 * t327 + t255 * t276) * qJD(1), t205 * t333 - t338, -t238 * t308 + (-t251 * t332 + t255 * t277) * qJD(1), -t238 * t314, -pkin(3) * t157 - t147 * t238 - t277 * t217 + (-pkin(8) * t331 + t197 * t251) * qJD(4) + (-t122 * t255 + t251 * t266) * qJD(1), pkin(3) * t156 + t148 * t238 + t276 * t217 + (pkin(8) * t333 + t197 * t254) * qJD(4) + (t123 * t255 + t254 * t266) * qJD(1), t147 * t207 + t148 * t205 + ((-t157 + t309) * pkin(8) + t369) * t254 + ((qJD(4) * t205 - t156) * pkin(8) + t368) * t251, -t122 * t147 - t123 * t148 + (-t197 - t354) * t217 + t259 * pkin(8), -t210 * t56 - t270 * t319, t137 * t319 + t56 * t209 - t210 * t267 - t270 * t318, -t319 * t232 + (qJD(3) * t210 - t270) * t314, t137 * t318 + t209 * t267, -t318 * t232 + (-qJD(3) * t209 + t137) * t314, -t232 * t314, t244 * t267 + t116 * t209 + t355 * t232 + t318 * t150 + t274 * t137 + (qJD(3) * t154 - t39) * t314, t116 * t210 - t244 * t56 + t356 * t232 - t319 * t150 + t274 * t270 + (-qJD(3) * t155 + t40) * t314, t137 * t356 + t154 * t56 - t155 * t267 + t209 * t279 - t9 * t210 - t270 * t355 - t318 * t40 + t319 * t39, t116 * t244 + t150 * t274 + t154 * t9 - t155 * t279 + t355 * t39 - t356 * t40, -t144 * t20 - t351 * t373, t143 * t20 - t144 * t21 + t351 * t72 + t352 * t373, -t351 * t226 + (qJD(3) * t144 - t373) * t314, t143 * t21 - t352 * t72, t352 * t226 + (-qJD(3) * t143 + t72) * t314, -t226 * t314, t143 * t38 + t169 * t21 - t352 * t82 + t320 * t72 + t357 * t226 + (qJD(3) * t62 - t10) * t314, t144 * t38 - t169 * t20 - t351 * t82 + t320 * t373 - t358 * t226 + (-qJD(3) * t63 + t11) * t314, t10 * t351 + t11 * t352 + t143 * t261 - t144 * t2 + t20 * t62 - t21 * t63 - t357 * t373 - t358 * t72, t10 * t357 + t11 * t358 + t169 * t38 + t2 * t62 - t261 * t63 + t320 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t336, -t205 ^ 2 + t207 ^ 2, -t156 + t337, -t336, t335 - t157, t240, -t197 * t207 - t368, t197 * t205 - t369, 0, 0, t342, t378, t375, -t342, t362, t240, -t232 * t45 + (-t137 * t207 - t232 * t305 + t240 * t253) * pkin(4) + t364, t232 * t46 + (-t207 * t270 - t232 * t304 - t240 * t250) * pkin(4) + t374, t40 * t270 + t46 * t137 - t39 * t137 + t45 * t270 + (-t250 * t267 + t253 * t56 + (-t137 * t253 + t250 * t270) * qJD(5)) * pkin(4), -t39 * t45 - t40 * t46 + (-t150 * t207 - t250 * t279 + t253 * t9 + (-t250 * t39 + t253 * t40) * qJD(5)) * pkin(4), t359, t379, t377, -t359, t363, t240, -t102 * t72 + t181 * t240 + t226 * t347 + t365, -t102 * t373 - t182 * t240 - t226 * t348 + t376, t181 * t20 - t182 * t21 - t347 * t373 - t348 * t72 + t386, t10 * t347 - t102 * t82 + t11 * t348 + t181 * t2 - t182 * t261; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t342, t378, t375, -t342, t362, t240, t232 * t40 + t364, t232 * t39 + t374, 0, 0, t359, t379, t377, -t359, t363, t240, -t12 * t226 + (-t226 * t303 + t240 * t360 - t270 * t72) * pkin(5) + t365, t13 * t226 + (-t226 * t282 - t240 * t249 - t270 * t373) * pkin(5) + t376, t12 * t373 + t13 * t72 + (t360 * t20 - t21 * t249 + (t249 * t373 - t360 * t72) * qJD(6)) * pkin(5) + t386, -t10 * t12 - t11 * t13 + (t360 * t2 - t261 * t249 - t270 * t82 + (-t10 * t249 + t11 * t360) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t359, t379, t377, -t359, t363, t240, t11 * t226 + t365, t10 * t226 + t376, 0, 0;];
tauc_reg  = t1;
