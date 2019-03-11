% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% 
% Output:
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRRPP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:02:35
% EndTime: 2019-03-09 21:02:52
% DurationCPUTime: 6.38s
% Computational Cost: add. (11931->485), mult. (29210->646), div. (0->0), fcn. (20438->8), ass. (0->251)
t242 = sin(qJ(3));
t243 = sin(qJ(2));
t309 = qJD(1) * t243;
t294 = t242 * t309;
t245 = cos(qJ(3));
t300 = t245 * qJD(2);
t185 = t294 - t300;
t307 = qJD(2) * t242;
t187 = t245 * t309 + t307;
t241 = sin(qJ(4));
t244 = cos(qJ(4));
t131 = t244 * t185 + t187 * t241;
t204 = -qJD(2) * pkin(2) + pkin(7) * t309;
t149 = pkin(3) * t185 + t204;
t100 = pkin(4) * t131 + qJD(5) + t149;
t239 = sin(pkin(10));
t240 = cos(pkin(10));
t265 = t185 * t241 - t244 * t187;
t356 = -t131 * t239 - t240 * t265;
t85 = t240 * t131 - t239 * t265;
t38 = pkin(5) * t85 - qJ(6) * t356 + t100;
t374 = t38 * t85;
t188 = t241 * t242 - t244 * t245;
t352 = qJD(3) + qJD(4);
t140 = t352 * t188;
t246 = cos(qJ(2));
t256 = t188 * t246;
t153 = qJD(1) * t256;
t373 = t140 - t153;
t189 = t241 * t245 + t242 * t244;
t141 = t352 * t189;
t308 = qJD(1) * t246;
t152 = t189 * t308;
t317 = t141 - t152;
t372 = pkin(5) * t356 + qJ(6) * t85;
t221 = -qJD(3) + t308;
t210 = -qJD(4) + t221;
t363 = qJ(5) * t265;
t199 = -pkin(2) * t246 - pkin(8) * t243 - pkin(1);
t178 = t199 * qJD(1);
t232 = pkin(7) * t308;
t205 = qJD(2) * pkin(8) + t232;
t137 = t245 * t178 - t205 * t242;
t110 = -pkin(9) * t187 + t137;
t103 = -pkin(3) * t221 + t110;
t330 = t242 * t178;
t138 = t205 * t245 + t330;
t111 = -pkin(9) * t185 + t138;
t107 = t241 * t111;
t62 = t244 * t103 - t107;
t51 = t62 + t363;
t48 = -pkin(4) * t210 + t51;
t364 = qJ(5) * t131;
t109 = t244 * t111;
t63 = t103 * t241 + t109;
t52 = t63 - t364;
t49 = t240 * t52;
t20 = t239 * t48 + t49;
t18 = -qJ(6) * t210 + t20;
t370 = t18 * t356;
t369 = t20 * t356;
t304 = qJD(3) * t242;
t271 = -t232 + (-t242 * t308 + t304) * pkin(3);
t368 = pkin(4) * t317 + t271;
t299 = qJD(1) * qJD(2);
t227 = t243 * t299;
t286 = t246 * t299;
t292 = t243 * t304;
t298 = qJD(2) * qJD(3);
t146 = -qJD(1) * t292 + (t286 + t298) * t245;
t270 = pkin(2) * t243 - pkin(8) * t246;
t196 = t270 * qJD(2);
t179 = qJD(1) * t196;
t274 = pkin(7) * t227;
t316 = -t245 * t179 - t242 * t274;
t253 = -qJD(3) * t138 - t316;
t70 = pkin(3) * t227 - t146 * pkin(9) + t253;
t305 = qJD(2) * t246;
t288 = t242 * t305;
t303 = qJD(3) * t245;
t254 = t243 * t303 + t288;
t147 = qJD(1) * t254 + t242 * t298;
t260 = t178 * t303 + t242 * t179 - t205 * t304;
t250 = -t245 * t274 + t260;
t81 = -pkin(9) * t147 + t250;
t284 = -t241 * t81 + t244 * t70;
t252 = -qJD(4) * t63 + t284;
t301 = qJD(4) * t244;
t302 = qJD(4) * t241;
t74 = t244 * t146 - t241 * t147 - t185 * t301 - t187 * t302;
t12 = pkin(4) * t227 - t74 * qJ(5) + qJD(5) * t265 + t252;
t251 = qJD(4) * t265 - t241 * t146 - t244 * t147;
t277 = -t103 * t301 + t111 * t302 - t241 * t70 - t244 * t81;
t15 = qJ(5) * t251 - qJD(5) * t131 - t277;
t3 = t240 * t12 - t239 * t15;
t2 = -pkin(5) * t227 - t3;
t351 = t356 * t38 + t2;
t367 = pkin(4) * t265;
t365 = t210 * t85;
t342 = t373 * t239 - t317 * t240;
t341 = t317 * t239 + t373 * t240;
t362 = t265 * t131;
t361 = -t131 ^ 2 + t265 ^ 2;
t360 = -t131 * t210 + t74;
t359 = t131 * t149 + t277;
t358 = t149 * t265 + t252;
t357 = t210 * t265 + t251;
t355 = t356 ^ 2;
t354 = -0.2e1 * t299;
t350 = pkin(8) + pkin(9);
t295 = qJD(3) * t350;
t195 = t245 * t295;
t175 = t244 * t195;
t194 = t242 * t295;
t206 = t350 * t242;
t207 = t350 * t245;
t313 = -t241 * t206 + t244 * t207;
t249 = t140 * qJ(5) - t313 * qJD(4) - t189 * qJD(5) + t241 * t194 - t175;
t326 = t245 * t246;
t261 = pkin(3) * t243 - pkin(9) * t326;
t193 = t270 * qJD(1);
t312 = pkin(7) * t294 + t245 * t193;
t118 = qJD(1) * t261 + t312;
t116 = t244 * t118;
t172 = t242 * t193;
t327 = t243 * t245;
t328 = t242 * t246;
t135 = t172 + (-pkin(7) * t327 - pkin(9) * t328) * qJD(1);
t64 = pkin(4) * t309 + qJ(5) * t153 - t135 * t241 + t116;
t255 = -t244 * t194 - t241 * t195 - t206 * t301 - t207 * t302;
t65 = -qJ(5) * t141 - qJD(5) * t188 + t255;
t320 = t241 * t118 + t244 * t135;
t67 = -qJ(5) * t152 + t320;
t345 = (-t249 + t64) * t240 + (t65 - t67) * t239;
t281 = -t110 * t241 - t109;
t259 = t281 + t364;
t297 = pkin(3) * t239 * t241;
t321 = t244 * t110 - t107;
t55 = t321 + t363;
t338 = -qJD(4) * t297 - t239 * t259 + (t301 * pkin(3) - t55) * t240;
t163 = t189 * t243;
t184 = t245 * t199;
t349 = pkin(7) * t242;
t136 = -pkin(9) * t327 + t184 + (-pkin(3) - t349) * t246;
t223 = pkin(7) * t326;
t311 = t242 * t199 + t223;
t329 = t242 * t243;
t142 = -pkin(9) * t329 + t311;
t319 = t241 * t136 + t244 * t142;
t353 = t246 * t300 - t292;
t25 = t239 * t51 + t49;
t348 = t25 * t356;
t4 = t239 * t12 + t240 * t15;
t164 = t188 * t243;
t306 = qJD(2) * t243;
t314 = t245 * t196 + t306 * t349;
t92 = t261 * qJD(2) + (-t223 + (pkin(9) * t243 - t199) * t242) * qJD(3) + t314;
t291 = t246 * t304;
t315 = t242 * t196 + t199 * t303;
t97 = -t254 * pkin(9) + (-t243 * t300 - t291) * pkin(7) + t315;
t283 = -t241 * t97 + t244 * t92;
t98 = -qJD(2) * t256 - t352 * t163;
t22 = pkin(4) * t306 - t98 * qJ(5) - qJD(4) * t319 + t164 * qJD(5) + t283;
t257 = t136 * t301 - t142 * t302 + t241 * t92 + t244 * t97;
t99 = -t302 * t329 + (t352 * t327 + t288) * t244 + t353 * t241;
t27 = -qJ(5) * t99 - qJD(5) * t163 + t257;
t8 = t239 * t22 + t240 * t27;
t36 = t239 * t64 + t240 * t67;
t31 = qJ(6) * t309 + t36;
t34 = t239 * t249 + t240 * t65;
t347 = t31 - t34;
t346 = pkin(5) * t309 + t345;
t129 = -t188 * t239 + t189 * t240;
t344 = -t342 * pkin(5) + t341 * qJ(6) - qJD(6) * t129 + t368;
t280 = t244 * t136 - t142 * t241;
t76 = -pkin(4) * t246 + qJ(5) * t164 + t280;
t82 = -qJ(5) * t163 + t319;
t46 = t239 * t76 + t240 * t82;
t343 = t239 * t52;
t340 = qJD(6) + t338;
t331 = t240 * t241;
t339 = -t239 * t55 + t240 * t259 + (t239 * t244 + t331) * qJD(4) * pkin(3);
t337 = t146 * t242;
t336 = t185 * t221;
t335 = t187 * t221;
t334 = t204 * t242;
t333 = t204 * t245;
t332 = t221 * t245;
t248 = qJD(1) ^ 2;
t325 = t246 * t248;
t247 = qJD(2) ^ 2;
t324 = t247 * t243;
t323 = t247 * t246;
t26 = t240 * t51 - t343;
t322 = qJD(6) - t26;
t229 = pkin(3) * t244 + pkin(4);
t167 = pkin(3) * t331 + t239 * t229;
t197 = pkin(3) * t329 + t243 * pkin(7);
t237 = t243 ^ 2;
t310 = -t246 ^ 2 + t237;
t296 = qJ(6) * t227 + t4;
t150 = pkin(3) * t254 + pkin(7) * t305;
t230 = -pkin(3) * t245 - pkin(2);
t289 = t221 * t303;
t287 = t339 * t356;
t127 = pkin(3) * t147 + pkin(7) * t286;
t41 = t239 * t74 - t240 * t251;
t282 = pkin(1) * t354;
t278 = -t244 * t206 - t207 * t241;
t198 = t210 * qJD(6);
t1 = -t198 + t296;
t276 = t185 + t300;
t275 = -t187 + t307;
t42 = t239 * t251 + t240 * t74;
t117 = -qJ(5) * t188 + t313;
t258 = -qJ(5) * t189 + t278;
t79 = t239 * t117 - t240 * t258;
t80 = t240 * t117 + t239 * t258;
t273 = -t34 * t85 - t80 * t41 + t79 * t42;
t272 = pkin(4) * t163 + t197;
t269 = pkin(3) * t187 - t367;
t268 = -t85 ^ 2 - t355;
t7 = t22 * t240 - t239 * t27;
t19 = t240 * t48 - t343;
t45 = -t239 * t82 + t240 * t76;
t264 = pkin(4) * t99 + t150;
t263 = qJD(1) * t237 - t221 * t246;
t262 = pkin(4) * t188 + t230;
t58 = -pkin(4) * t251 + t127;
t166 = t229 * t240 - t297;
t9 = pkin(5) * t41 - qJ(6) * t42 - qJD(6) * t356 + t58;
t225 = -pkin(4) * t240 - pkin(5);
t224 = pkin(4) * t239 + qJ(6);
t157 = -pkin(5) - t166;
t156 = qJ(6) + t167;
t128 = t240 * t188 + t189 * t239;
t113 = -t163 * t239 - t164 * t240;
t112 = t240 * t163 - t164 * t239;
t78 = pkin(5) * t128 - qJ(6) * t129 + t262;
t59 = pkin(5) * t112 - qJ(6) * t113 + t272;
t57 = -t239 * t99 + t240 * t98;
t56 = t239 * t98 + t240 * t99;
t44 = -t367 + t372;
t43 = pkin(5) * t246 - t45;
t40 = -qJ(6) * t246 + t46;
t39 = t269 + t372;
t17 = pkin(5) * t210 + qJD(6) - t19;
t16 = pkin(5) * t56 - qJ(6) * t57 - qJD(6) * t113 + t264;
t6 = -pkin(5) * t306 - t7;
t5 = qJ(6) * t306 - qJD(6) * t246 + t8;
t10 = [0, 0, 0, 0.2e1 * t246 * t227, t310 * t354, t323, -t324, 0, -pkin(7) * t323 + t243 * t282, pkin(7) * t324 + t246 * t282, t146 * t327 + t353 * t187 (-t185 * t245 - t187 * t242) * t305 + (-t337 - t147 * t245 + (t185 * t242 - t187 * t245) * qJD(3)) * t243, t221 * t292 - t146 * t246 + (t187 * t243 + t245 * t263) * qJD(2), t243 * t289 + t147 * t246 + (-t185 * t243 - t242 * t263) * qJD(2) (-t221 - t308) * t306 -(-t199 * t304 + t314) * t221 + (t204 * t303 + pkin(7) * t147 + (qJD(1) * t184 + t137) * qJD(2)) * t243 + ((pkin(7) * t185 + t334) * qJD(2) + (t330 + (pkin(7) * t221 + t205) * t245) * qJD(3) + t316) * t246 (-pkin(7) * t291 + t315) * t221 + t260 * t246 + (pkin(7) * t146 - t204 * t304) * t243 + ((pkin(7) * t187 + t333) * t246 + (-pkin(7) * t332 - qJD(1) * t311 - t138) * t243) * qJD(2), -t164 * t74 - t265 * t98, -t131 * t98 - t163 * t74 - t164 * t251 + t265 * t99, -t98 * t210 - t74 * t246 + (-qJD(1) * t164 - t265) * t306, t99 * t210 - t251 * t246 + (-qJD(1) * t163 - t131) * t306 (-t210 - t308) * t306, -t283 * t210 - t284 * t246 + t150 * t131 - t197 * t251 + t127 * t163 + t149 * t99 + (t210 * t319 + t246 * t63) * qJD(4) + (qJD(1) * t280 + t62) * t306, t257 * t210 - t277 * t246 - t150 * t265 + t197 * t74 - t127 * t164 + t149 * t98 + (-qJD(1) * t319 - t63) * t306, -t112 * t4 - t113 * t3 - t19 * t57 - t20 * t56 - t356 * t7 - t41 * t46 - t42 * t45 - t8 * t85, t100 * t264 + t19 * t7 + t20 * t8 + t272 * t58 + t3 * t45 + t4 * t46, t9 * t112 + t16 * t85 + t2 * t246 + t6 * t210 + t38 * t56 + t59 * t41 + (-qJD(1) * t43 - t17) * t306, -t1 * t112 + t113 * t2 + t17 * t57 - t18 * t56 + t356 * t6 - t40 * t41 + t42 * t43 - t5 * t85, -t1 * t246 - t9 * t113 - t16 * t356 - t5 * t210 - t38 * t57 - t59 * t42 + (qJD(1) * t40 + t18) * t306, t1 * t40 + t16 * t38 + t17 * t6 + t18 * t5 + t2 * t43 + t59 * t9; 0, 0, 0, -t243 * t325, t310 * t248, 0, 0, 0, t248 * pkin(1) * t243, pkin(1) * t325, -t187 * t332 + t337 (t146 + t336) * t245 + (-t147 + t335) * t242, -t289 + (t221 * t326 + t243 * t275) * qJD(1), t221 * t304 + (-t221 * t328 + t243 * t276) * qJD(1), t221 * t309, -pkin(2) * t147 + t312 * t221 + (pkin(8) * t332 + t334) * qJD(3) + ((-pkin(8) * t307 - t137) * t243 + (-pkin(7) * t276 - t334) * t246) * qJD(1), -pkin(2) * t146 - t172 * t221 + (-pkin(8) * t221 * t242 + t333) * qJD(3) + (-t204 * t326 + (-pkin(8) * t300 + t138) * t243 + (t221 * t327 + t246 * t275) * pkin(7)) * qJD(1), t74 * t189 + t265 * t373, t131 * t373 - t74 * t188 + t189 * t251 + t265 * t317, t373 * t210 + (qJD(2) * t189 + t265) * t309, t317 * t210 + (-qJD(2) * t188 + t131) * t309, t210 * t309, t127 * t188 - t230 * t251 + (t207 * t301 + t116 + t175 + (-qJD(4) * t206 - t135 - t194) * t241) * t210 + t317 * t149 + t271 * t131 + (qJD(2) * t278 - t62) * t309, t127 * t189 + t230 * t74 + (t255 - t320) * t210 - t373 * t149 - t271 * t265 + (-qJD(2) * t313 + t63) * t309, -t4 * t128 - t3 * t129 + t19 * t341 + t20 * t342 + t345 * t356 + t36 * t85 + t273, t4 * t80 - t3 * t79 + t58 * t262 + (t34 - t36) * t20 - t345 * t19 + t368 * t100, t9 * t128 + t78 * t41 + t344 * t85 - t342 * t38 + t346 * t210 + (-qJD(2) * t79 + t17) * t309, -t1 * t128 + t2 * t129 - t17 * t341 + t18 * t342 + t31 * t85 + t346 * t356 + t273, -t9 * t129 - t78 * t42 - t344 * t356 + t341 * t38 + t347 * t210 + (qJD(2) * t80 - t18) * t309, t1 * t80 + t17 * t346 - t18 * t347 + t2 * t79 + t344 * t38 + t9 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t187 * t185, -t185 ^ 2 + t187 ^ 2, t146 - t336, -t147 - t335, t227, -t138 * t221 - t204 * t187 + t253, -t137 * t221 + t185 * t204 - t250, -t362, t361, t360, t357, t227, t281 * t210 + (-t131 * t187 + t210 * t302 + t227 * t244) * pkin(3) + t358, -t321 * t210 + (t187 * t265 + t210 * t301 - t227 * t241) * pkin(3) + t359, -t166 * t42 - t167 * t41 + t287 + t369 + (-t19 - t338) * t85, -t100 * t269 + t3 * t166 + t4 * t167 - t19 * t339 + t20 * t338, -t157 * t227 + t210 * t339 - t39 * t85 - t351, -t156 * t41 + t157 * t42 + t287 + t370 + (t17 - t340) * t85, t156 * t227 - t210 * t340 + t356 * t39 + t1 - t374, t1 * t156 + t2 * t157 + t17 * t339 + t18 * t340 - t38 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t362, t361, t360, t357, t227, -t63 * t210 + t358, -t210 * t62 + t359, t369 - t348 + (-t239 * t41 - t240 * t42) * pkin(4) + (-t19 + t26) * t85, t19 * t25 - t20 * t26 + (t100 * t265 + t239 * t4 + t240 * t3) * pkin(4), -t210 * t25 - t225 * t227 - t44 * t85 - t351, t370 - t224 * t41 + t225 * t42 - t348 + (t17 - t322) * t85, t210 * t26 + t224 * t227 + t356 * t44 - 0.2e1 * t198 + t296 - t374, t1 * t224 - t17 * t25 + t18 * t322 + t2 * t225 - t38 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t268, t19 * t356 + t20 * t85 + t58, -t210 * t356 + t41, t268, -t42 - t365, -t17 * t356 + t18 * t85 + t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t356 * t85 - t227, t42 - t365, -t210 ^ 2 - t355, t18 * t210 + t351;];
tauc_reg  = t10;
