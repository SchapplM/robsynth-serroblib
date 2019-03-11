% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRPRP5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:51:33
% EndTime: 2019-03-09 16:51:50
% DurationCPUTime: 6.60s
% Computational Cost: add. (10511->472), mult. (25868->627), div. (0->0), fcn. (18551->8), ass. (0->244)
t224 = sin(qJ(3));
t225 = sin(qJ(2));
t307 = qJD(1) * t225;
t289 = t224 * t307;
t226 = cos(qJ(3));
t297 = t226 * qJD(2);
t185 = t289 - t297;
t305 = qJD(2) * t224;
t187 = t226 * t307 + t305;
t222 = sin(pkin(10));
t332 = cos(pkin(10));
t134 = -t222 * t185 + t187 * t332;
t223 = sin(qJ(5));
t294 = qJD(2) * qJD(3);
t300 = qJD(3) * t226;
t227 = cos(qJ(2));
t303 = qJD(2) * t227;
t364 = t224 * t303 + t225 * t300;
t146 = qJD(1) * t364 + t224 * t294;
t301 = qJD(3) * t225;
t287 = t224 * t301;
t255 = t227 * t297 - t287;
t237 = qJD(1) * t255 + t226 * t294;
t235 = -t222 * t146 + t237 * t332;
t236 = -t146 * t332 - t222 * t237;
t259 = t185 * t332 + t222 * t187;
t353 = cos(qJ(5));
t242 = t353 * t259;
t299 = qJD(5) * t223;
t33 = qJD(5) * t242 + t134 * t299 - t223 * t236 - t353 * t235;
t306 = qJD(1) * t227;
t209 = -qJD(3) + t306;
t199 = -qJD(5) + t209;
t83 = t134 * t223 + t242;
t339 = t199 * t83;
t16 = -t33 - t339;
t323 = t226 * t227;
t268 = pkin(3) * t225 - qJ(4) * t323;
t345 = -qJ(4) - pkin(8);
t282 = qJD(3) * t345;
t271 = pkin(2) * t225 - pkin(8) * t227;
t188 = t271 * qJD(1);
t313 = pkin(7) * t289 + t226 * t188;
t381 = -qJD(1) * t268 - t224 * qJD(4) + t226 * t282 - t313;
t172 = t224 * t188;
t296 = t226 * qJD(4);
t324 = t225 * t226;
t380 = t172 + (-qJ(4) * t224 * t227 - pkin(7) * t324) * qJD(1) - t224 * t282 - t296;
t347 = t83 ^ 2;
t363 = t134 * t353 - t223 * t259;
t375 = t363 ^ 2;
t379 = -t347 + t375;
t340 = qJD(2) * pkin(2);
t197 = pkin(7) * t307 - t340;
t144 = t185 * pkin(3) + qJD(4) + t197;
t92 = pkin(4) * t259 + t144;
t27 = t83 * pkin(5) - qJ(6) * t363 + t92;
t378 = t27 * t83;
t377 = t83 * t92;
t346 = t363 * t83;
t179 = t222 * t226 + t224 * t332;
t249 = t227 * t179;
t151 = qJD(1) * t249;
t247 = qJD(3) * t179;
t369 = t151 - t247;
t258 = t222 * t224 - t226 * t332;
t250 = t227 * t258;
t152 = qJD(1) * t250;
t356 = qJD(3) * t258;
t376 = t152 - t356;
t41 = pkin(5) * t363 + qJ(6) * t83;
t361 = t380 * t222 + t381 * t332;
t360 = t381 * t222 - t380 * t332;
t373 = pkin(4) * t307 + t376 * pkin(9) - t361;
t372 = t369 * pkin(9) + t360;
t349 = t27 * t363;
t371 = t363 * t92;
t338 = t363 * t199;
t34 = qJD(5) * t363 + t223 * t235 - t353 * t236;
t370 = -t34 - t338;
t368 = pkin(9) * t134;
t367 = pkin(9) * t259;
t241 = t353 * t258;
t334 = -t179 * t299 + (-qJD(3) - qJD(5)) * t241 + t152 * t353 + t369 * t223;
t252 = t223 * t258;
t284 = qJD(5) * t353;
t290 = t353 * t179;
t333 = qJD(3) * t290 - qJD(5) * t252 - t151 * t353 + t179 * t284 + t376 * t223;
t215 = pkin(7) * t306;
t302 = qJD(3) * t224;
t365 = -t215 + (-t224 * t306 + t302) * pkin(3);
t295 = qJD(1) * qJD(2);
t362 = -0.2e1 * t295;
t212 = pkin(3) * t332 + pkin(4);
t352 = pkin(3) * t222;
t293 = t223 * t352;
t192 = -pkin(2) * t227 - pkin(8) * t225 - pkin(1);
t176 = t192 * qJD(1);
t198 = qJD(2) * pkin(8) + t215;
t139 = t226 * t176 - t198 * t224;
t108 = -qJ(4) * t187 + t139;
t326 = t224 * t176;
t140 = t226 * t198 + t326;
t109 = -qJ(4) * t185 + t140;
t281 = t332 * t109;
t63 = -t108 * t222 - t281;
t49 = t63 + t367;
t104 = t222 * t109;
t64 = t332 * t108 - t104;
t50 = t64 - t368;
t359 = -qJD(5) * t293 + t212 * t284 - t223 * t49 - t353 * t50;
t194 = t345 * t224;
t195 = t345 * t226;
t142 = t332 * t194 + t195 * t222;
t119 = -pkin(9) * t179 + t142;
t143 = t222 * t194 - t332 * t195;
t120 = -pkin(9) * t258 + t143;
t261 = t119 * t353 - t223 * t120;
t358 = qJD(5) * t261 - t373 * t223 + t372 * t353;
t72 = t223 * t119 + t120 * t353;
t357 = qJD(5) * t72 + t372 * t223 + t373 * t353;
t318 = -t369 * pkin(4) + t365;
t312 = t223 * t212 + t353 * t352;
t220 = t225 ^ 2;
t298 = t220 * qJD(1);
t163 = t258 * t225;
t181 = t226 * t192;
t351 = pkin(7) * t224;
t136 = -qJ(4) * t324 + t181 + (-pkin(3) - t351) * t227;
t211 = pkin(7) * t323;
t311 = t224 * t192 + t211;
t325 = t224 * t225;
t141 = -qJ(4) * t325 + t311;
t89 = t332 * t136 - t141 * t222;
t67 = -pkin(4) * t227 + pkin(9) * t163 + t89;
t251 = t225 * t179;
t90 = t222 * t136 + t332 * t141;
t73 = -pkin(9) * t251 + t90;
t263 = t223 * t67 + t353 * t73;
t240 = t225 * t247;
t253 = t268 * qJD(2);
t189 = t271 * qJD(2);
t304 = qJD(2) * t225;
t314 = t226 * t189 + t304 * t351;
t78 = -t225 * t296 + t253 + (-t211 + (qJ(4) * t225 - t192) * t224) * qJD(3) + t314;
t315 = t224 * t189 + t192 * t300;
t88 = (-pkin(7) * qJD(2) - qJ(4) * qJD(3)) * t324 + (-qJD(4) * t225 + (-pkin(7) * qJD(3) - qJ(4) * qJD(2)) * t227) * t224 + t315;
t39 = -t222 * t88 + t332 * t78;
t29 = pkin(9) * t240 + (t225 * pkin(4) + pkin(9) * t250) * qJD(2) + t39;
t233 = -qJD(2) * t249 + t225 * t356;
t40 = t222 * t78 + t332 * t88;
t32 = pkin(9) * t233 + t40;
t355 = -qJD(5) * t263 - t223 * t32 + t29 * t353;
t350 = pkin(8) * t209;
t344 = qJ(6) * t307 - t358;
t343 = pkin(5) * t307 + t357;
t130 = t290 - t252;
t342 = t333 * pkin(5) - t334 * qJ(6) - t130 * qJD(6) + t318;
t177 = qJD(1) * t189;
t213 = t225 * t295;
t273 = pkin(7) * t213;
t316 = t226 * t177 + t224 * t273;
t55 = qJD(1) * t253 - t109 * qJD(3) - t187 * qJD(4) + t316;
t317 = t176 * t300 + t224 * t177;
t238 = -t198 * t302 - t226 * t273 + t317;
t61 = -qJ(4) * t146 - qJD(4) * t185 + t238;
t23 = t222 * t55 + t332 * t61;
t100 = -pkin(3) * t209 + t108;
t56 = t332 * t100 - t104;
t42 = -pkin(4) * t209 - t368 + t56;
t57 = t222 * t100 + t281;
t45 = t57 - t367;
t9 = -t223 * t45 + t353 * t42;
t337 = qJD(6) - t9;
t336 = qJD(6) + t359;
t335 = t312 * qJD(5) - t223 * t50 + t353 * t49;
t330 = qJD(2) * t261;
t329 = qJD(2) * t72;
t328 = t187 * t209;
t327 = t197 * t224;
t322 = t227 * t209;
t229 = qJD(1) ^ 2;
t321 = t227 * t229;
t228 = qJD(2) ^ 2;
t320 = t228 * t225;
t319 = t228 * t227;
t309 = pkin(3) * t325 + t225 * pkin(7);
t308 = -t227 ^ 2 + t220;
t292 = t364 * pkin(3) + pkin(7) * t303;
t291 = -t226 * pkin(3) - pkin(2);
t285 = t209 * t300;
t131 = t146 * pkin(3) + qJD(2) * t215;
t22 = -t222 * t61 + t332 * t55;
t13 = pkin(4) * t213 - pkin(9) * t235 + t22;
t17 = pkin(9) * t236 + t23;
t280 = -t223 * t13 - t353 * t17 - t42 * t284 + t45 * t299;
t279 = -t353 * t13 + t223 * t17 + t45 * t284 + t42 * t299;
t278 = pkin(1) * t362;
t277 = t209 + t306;
t276 = t185 + t297;
t275 = -t187 + t305;
t274 = qJD(3) + t306;
t272 = pkin(5) * t213;
t101 = pkin(3) * t187 + pkin(4) * t134;
t191 = t199 * qJD(6);
t206 = qJ(6) * t213;
t1 = t206 - t191 - t280;
t269 = -t199 * t9 + t280;
t10 = t223 * t42 + t353 * t45;
t267 = -t10 * t199 - t279;
t265 = -t223 * t73 + t353 * t67;
t262 = t223 * t29 + t67 * t284 - t299 * t73 + t353 * t32;
t260 = t274 * t305;
t257 = t212 * t353 - t293;
t256 = t199 * t335 - t279;
t2 = -t272 + t279;
t245 = t33 - t339;
t244 = t223 * t251;
t239 = t225 * t290;
t153 = pkin(4) * t258 + t291;
t137 = pkin(4) * t251 + t309;
t234 = -qJD(2) * t250 - t240;
t91 = -pkin(4) * t233 + t292;
t69 = -pkin(4) * t236 + t131;
t5 = t34 * pkin(5) + t33 * qJ(6) - qJD(6) * t363 + t69;
t230 = t34 - t338;
t164 = -pkin(5) - t257;
t162 = qJ(6) + t312;
t129 = t179 * t223 + t241;
t114 = -t163 * t353 - t244;
t113 = -t163 * t223 + t239;
t68 = t129 * pkin(5) - t130 * qJ(6) + t153;
t51 = t113 * pkin(5) - t114 * qJ(6) + t137;
t48 = -qJD(5) * t244 - t163 * t284 + t223 * t234 - t233 * t353;
t47 = qJD(5) * t239 - t163 * t299 - t223 * t233 - t234 * t353;
t36 = t227 * pkin(5) - t265;
t35 = -qJ(6) * t227 + t263;
t31 = t101 + t41;
t8 = -t199 * qJ(6) + t10;
t7 = t199 * pkin(5) + t337;
t6 = t48 * pkin(5) + t47 * qJ(6) - t114 * qJD(6) + t91;
t4 = -pkin(5) * t304 - t355;
t3 = qJ(6) * t304 - qJD(6) * t227 + t262;
t11 = [0, 0, 0, 0.2e1 * t227 * t213, t308 * t362, t319, -t320, 0, -pkin(7) * t319 + t225 * t278, pkin(7) * t320 + t227 * t278, t187 * t255 + t237 * t324 (-t226 * t185 - t187 * t224) * t303 + ((t185 + t289) * t302 + (-t187 * qJD(3) - t146 - t260) * t226) * t225, t277 * t287 + (t187 * t225 + (t298 + (-t209 - t274) * t227) * t226) * qJD(2), t225 * t285 + t146 * t227 + (-t185 * t225 + (-t298 + t322) * t224) * qJD(2), -t277 * t304 -(-t192 * t302 + t314) * t209 + (t197 * t300 + pkin(7) * t146 + (qJD(1) * t181 + t139) * qJD(2)) * t225 + ((pkin(7) * t185 + t327) * qJD(2) + (t326 + (pkin(7) * t209 + t198) * t226) * qJD(3) - t316) * t227, t315 * t209 + t317 * t227 + (-t197 * t225 - t198 * t227 + (-t322 - t298) * pkin(7)) * t302 + ((pkin(7) * t187 + t197 * t226) * t227 + (-t311 * qJD(1) - t140 + (-t209 + t274) * pkin(7) * t226) * t225) * qJD(2), -t40 * t259 + t90 * t236 - t23 * t251 + t57 * t233 - t39 * t134 - t89 * t235 + t22 * t163 + t56 * (t179 * t301 + t258 * t303) t131 * t309 + t144 * t292 + t22 * t89 + t23 * t90 + t56 * t39 + t57 * t40, -t114 * t33 - t363 * t47, t113 * t33 - t114 * t34 - t363 * t48 + t47 * t83, t47 * t199 + t33 * t227 + (qJD(1) * t114 + t363) * t304, t48 * t199 + t34 * t227 + (-qJD(1) * t113 - t83) * t304 (-t199 - t306) * t304, t69 * t113 + t137 * t34 - t355 * t199 + t265 * t213 + t279 * t227 + t9 * t304 + t92 * t48 + t91 * t83, t262 * t199 - t280 * t227 + t91 * t363 - t137 * t33 + t69 * t114 - t92 * t47 + (-qJD(1) * t263 - t10) * t304, t5 * t113 + t4 * t199 + t2 * t227 + t27 * t48 + t51 * t34 + t6 * t83 + (-qJD(1) * t36 - t7) * t304, -t1 * t113 + t114 * t2 - t3 * t83 - t33 * t36 - t34 * t35 + t363 * t4 - t47 * t7 - t48 * t8, -t1 * t227 - t5 * t114 - t3 * t199 + t27 * t47 + t51 * t33 - t6 * t363 + (qJD(1) * t35 + t8) * t304, t1 * t35 + t2 * t36 + t27 * t6 + t3 * t8 + t4 * t7 + t5 * t51; 0, 0, 0, -t225 * t321, t308 * t229, 0, 0, 0, t229 * pkin(1) * t225, pkin(1) * t321, -t224 ^ 2 * qJD(1) * t301 + (t260 - t328) * t226 (-t146 + t328) * t224 + ((-t185 + t297) * qJD(3) + (t227 * t276 - t287) * qJD(1)) * t226, -t285 + (t225 * t275 + t226 * t322) * qJD(1), t209 * t302 + (-t224 * t322 + t225 * t276) * qJD(1), t209 * t307, -pkin(2) * t146 + t313 * t209 + (t226 * t350 + t327) * qJD(3) + ((-pkin(8) * t305 - t139) * t225 + (-pkin(7) * t276 - t327) * t227) * qJD(1), -t172 * t209 + (-t224 * t350 + (t197 - t340) * t226) * qJD(3) + ((-t197 - t340) * t323 + (pkin(2) * t302 - pkin(8) * t297 + t140) * t225 + (t209 * t324 + t227 * t275) * pkin(7)) * qJD(1), -t361 * t134 - t142 * t235 + t143 * t236 - t22 * t179 - t23 * t258 - t360 * t259 + t369 * t57 - t376 * t56, t131 * t291 + t22 * t142 + t23 * t143 + t365 * t144 + t360 * t57 + t361 * t56, -t130 * t33 + t334 * t363, t129 * t33 - t130 * t34 - t333 * t363 - t334 * t83, -t334 * t199 + (qJD(2) * t130 - t363) * t307, t333 * t199 + (-qJD(2) * t129 + t83) * t307, t199 * t307, t69 * t129 + t153 * t34 + t333 * t92 + t318 * t83 + t357 * t199 + (-t9 + t330) * t307, t69 * t130 - t153 * t33 + t334 * t92 + t318 * t363 + t358 * t199 + (t10 - t329) * t307, t5 * t129 + t68 * t34 + t342 * t83 + t333 * t27 + t343 * t199 + (t7 + t330) * t307, -t1 * t129 + t130 * t2 + t261 * t33 - t333 * t8 + t334 * t7 - t34 * t72 + t343 * t363 + t344 * t83, -t5 * t130 + t68 * t33 - t342 * t363 - t334 * t27 + t344 * t199 + (-t8 + t329) * t307, t1 * t72 - t2 * t261 + t27 * t342 + t343 * t7 - t344 * t8 + t5 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t187 * t185, -t185 ^ 2 + t187 ^ 2, -t185 * t209 + t237, -t146 - t328, t213, -t197 * t187 + t316 + (-qJD(3) - t209) * t140, -t139 * t209 + t185 * t197 - t238 (-t222 ^ 2 - t332 ^ 2) * pkin(3) * t237 + (-t56 + t64) * t259 + (t57 + t63) * t134, -t56 * t63 - t57 * t64 + (-t144 * t187 + t22 * t332 + t222 * t23) * pkin(3), t346, t379, t16, t370, t213, -t101 * t83 + t213 * t257 + t256 - t371, -t101 * t363 + t359 * t199 - t312 * t213 + t280 + t377, -t349 - t31 * t83 + (pkin(5) - t164) * t213 + t256, -t162 * t34 - t164 * t33 + (t335 + t8) * t363 + (-t336 + t7) * t83, t162 * t213 - t199 * t336 + t31 * t363 + t1 - t378, t1 * t162 + t164 * t2 - t27 * t31 + t335 * t7 + t336 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134 ^ 2 - t259 ^ 2, t134 * t56 + t259 * t57 + t131, 0, 0, 0, 0, 0, t230, -t245, t230, -t347 - t375, t245, -t363 * t7 + t8 * t83 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t346, t379, t16, t370, t213, t267 - t371, t269 + t377, -t41 * t83 + t267 + 0.2e1 * t272 - t349, pkin(5) * t33 - qJ(6) * t34 + (-t10 + t8) * t363 + (t7 - t337) * t83, t363 * t41 - 0.2e1 * t191 + 0.2e1 * t206 - t269 - t378, -pkin(5) * t2 + qJ(6) * t1 - t10 * t7 - t27 * t41 + t337 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t213 + t346, t16, -t199 ^ 2 - t375, t199 * t8 + t2 + t349;];
tauc_reg  = t11;
