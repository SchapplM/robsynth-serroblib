% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRPRP9_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP9_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:27:15
% EndTime: 2019-03-09 17:27:30
% DurationCPUTime: 5.65s
% Computational Cost: add. (5790->514), mult. (13275->648), div. (0->0), fcn. (8518->6), ass. (0->242)
t199 = cos(qJ(2));
t196 = sin(qJ(2));
t335 = pkin(8) * t196;
t231 = -pkin(2) * t199 - t335;
t144 = -pkin(1) + t231;
t117 = t144 * qJD(1);
t284 = qJD(1) * t199;
t182 = pkin(7) * t284;
t151 = qJD(2) * pkin(8) + t182;
t195 = sin(qJ(3));
t198 = cos(qJ(3));
t78 = t198 * t117 - t195 * t151;
t296 = qJD(4) - t78;
t277 = qJD(3) * t198;
t364 = t198 * t284 - t277;
t261 = t195 * t284;
t279 = qJD(3) * t195;
t363 = t261 - t279;
t230 = pkin(2) * t196 - pkin(8) * t199;
t139 = t230 * qJD(1);
t115 = t195 * t139;
t285 = qJD(1) * t196;
t293 = qJ(4) * t285 + t115;
t302 = t196 * t198;
t303 = t195 * t199;
t336 = pkin(8) - pkin(9);
t362 = t336 * t279 + (-pkin(7) * t302 + pkin(9) * t303) * qJD(1) + t293;
t153 = t336 * t198;
t263 = -pkin(7) * t195 - pkin(3);
t301 = t198 * t199;
t208 = -pkin(9) * t301 + (-pkin(4) + t263) * t196;
t309 = t139 * t198;
t361 = -qJD(1) * t208 + qJD(3) * t153 + t309;
t260 = t198 * t285;
t283 = qJD(2) * t195;
t131 = t260 + t283;
t358 = -pkin(9) * t131 + t296;
t262 = t195 * t285;
t273 = t198 * qJD(2);
t129 = t262 - t273;
t194 = sin(qJ(5));
t197 = cos(qJ(5));
t222 = -t197 * t129 + t131 * t194;
t75 = t129 * t194 + t131 * t197;
t340 = t75 ^ 2;
t360 = -t222 ^ 2 + t340;
t193 = qJD(2) * pkin(2);
t150 = pkin(7) * t285 - t193;
t68 = t129 * pkin(3) - t131 * qJ(4) + t150;
t53 = -pkin(4) * t129 - t68;
t12 = pkin(5) * t222 - qJ(6) * t75 + t53;
t334 = t12 * t75;
t359 = t53 * t75;
t275 = qJD(5) * t197;
t276 = qJD(5) * t194;
t327 = -t364 * t194 + t195 * t275 + t363 * t197 - t198 * t276;
t133 = t194 * t195 + t197 * t198;
t213 = t133 * t199;
t81 = t133 * qJD(5) - t194 * t279 - t197 * t277;
t314 = -qJD(1) * t213 - t81;
t357 = t364 * qJ(4) - t195 * qJD(4) - t182;
t37 = pkin(5) * t75 + qJ(6) * t222;
t171 = -qJD(3) + t284;
t160 = qJD(5) + t171;
t253 = t196 * t277;
t281 = qJD(2) * t199;
t258 = t195 * t281;
t209 = t253 + t258;
t271 = qJD(2) * qJD(3);
t249 = t195 * t271;
t204 = qJD(1) * t209 + t249;
t272 = qJD(1) * qJD(2);
t250 = t199 * t272;
t278 = qJD(3) * t196;
t255 = t195 * t278;
t96 = qJD(1) * t255 + (-t250 - t271) * t198;
t27 = qJD(5) * t75 - t194 * t96 - t197 * t204;
t349 = t160 * t75 - t27;
t175 = t196 * t272;
t338 = t160 ^ 2;
t356 = t131 * t75 - t175 * t194 + t338 * t197;
t354 = -0.2e1 * t272;
t170 = pkin(5) * t175;
t257 = t196 * t273;
t142 = t230 * qJD(2);
t119 = qJD(1) * t142;
t294 = t117 * t277 + t195 * t119;
t325 = pkin(9) * qJD(2);
t154 = t171 * qJD(4);
t162 = qJ(4) * t175;
t344 = t162 - t154;
t18 = (-t151 + t325) * t279 + (-pkin(7) * t257 + pkin(9) * t209) * qJD(1) + t294 + t344;
t236 = pkin(7) * t175;
t235 = -t117 * t279 + t198 * t119 - t151 * t277 + t195 * t236;
t337 = pkin(3) + pkin(4);
t19 = pkin(9) * t96 - t337 * t175 - t235;
t42 = t337 * t171 + t358;
t156 = t171 * qJ(4);
t79 = t195 * t117 + t198 * t151;
t57 = pkin(9) * t129 + t79;
t51 = -t156 + t57;
t246 = -t194 * t18 + t197 * t19 - t51 * t275 - t42 * t276;
t2 = t170 - t246;
t11 = t194 * t42 + t197 * t51;
t8 = qJ(6) * t160 + t11;
t353 = t160 * t8 - t2;
t265 = t337 * t195;
t326 = -qJD(3) * t265 + t337 * t261 - t357;
t332 = t75 * t222;
t352 = qJ(4) * t276 + t194 * t57 - t197 * t358 + t275 * t337;
t152 = t336 * t195;
t93 = t152 * t194 + t153 * t197;
t351 = -qJD(5) * t93 + t194 * t362 + t361 * t197;
t26 = -t129 * t275 + t131 * t276 - t194 * t204 + t197 * t96;
t350 = -t160 * t222 + t26;
t237 = pkin(3) * t175;
t38 = -t235 - t237;
t65 = -t156 + t79;
t348 = t171 * t65 + t38;
t221 = t152 * t197 - t153 * t194;
t347 = -qJD(5) * t221 - t361 * t194 + t197 * t362;
t239 = -t131 + t283;
t345 = t199 * t239;
t287 = t197 * qJ(4) - t194 * t337;
t247 = -t197 * t18 - t194 * t19 - t42 * t275 + t51 * t276;
t343 = -t12 * t222 - t247;
t342 = -t53 * t222 - t247;
t173 = pkin(7) * t303;
t190 = t199 * pkin(3);
t69 = pkin(4) * t199 + t173 + t190 + (-pkin(9) * t196 - t144) * t198;
t305 = t195 * t196;
t174 = pkin(7) * t301;
t289 = t195 * t144 + t174;
t94 = -qJ(4) * t199 + t289;
t77 = pkin(9) * t305 + t94;
t224 = t194 * t69 + t197 * t77;
t229 = -qJD(3) * t174 + t142 * t198 - t144 * t279;
t34 = pkin(9) * t255 + qJD(2) * t208 - t229;
t282 = qJD(2) * t196;
t292 = t195 * t142 + t144 * t277;
t264 = qJ(4) * t282 + t292;
t36 = (-pkin(7) * qJD(2) + pkin(9) * qJD(3)) * t302 + (-qJD(4) + (-pkin(7) * qJD(3) + t325) * t195) * t199 + t264;
t341 = -qJD(5) * t224 - t194 * t36 + t197 * t34;
t339 = t131 ^ 2;
t304 = t195 * t197;
t220 = t194 * t198 - t304;
t331 = t327 * pkin(5) - t314 * qJ(6) + qJD(6) * t220 + t326;
t330 = -qJ(6) * t285 + t347;
t329 = pkin(5) * t285 + t351;
t324 = t131 * t68;
t268 = -pkin(3) * t204 - t96 * qJ(4) + t131 * qJD(4);
t35 = pkin(7) * t250 - t268;
t320 = t35 * t195;
t319 = t35 * t198;
t318 = t68 * t199;
t317 = -qJD(6) - t352;
t316 = t287 * qJD(5) + t358 * t194 + t197 * t57;
t315 = -t363 * pkin(3) + t357;
t313 = qJD(2) * t221;
t312 = qJD(2) * t93;
t311 = t129 * t171;
t310 = t131 * t129;
t308 = t150 * t198;
t307 = t171 * t195;
t306 = t171 * t198;
t186 = t195 * qJ(4);
t300 = t199 * t171;
t202 = qJD(1) ^ 2;
t299 = t199 * t202;
t201 = qJD(2) ^ 2;
t298 = t201 * t196;
t297 = t201 * t199;
t10 = -t194 * t51 + t197 * t42;
t295 = qJD(6) - t10;
t85 = t131 * pkin(3) + t129 * qJ(4);
t256 = t199 * t273;
t290 = qJ(4) * t256 + qJD(4) * t302;
t191 = t196 ^ 2;
t286 = -t199 ^ 2 + t191;
t280 = qJD(3) * t129;
t274 = t191 * qJD(1);
t270 = pkin(8) * t307;
t269 = pkin(8) * t306;
t143 = -t198 * pkin(3) - pkin(2) - t186;
t267 = pkin(8) * t273;
t266 = pkin(7) * t281;
t254 = t199 * t279;
t245 = qJD(1) * t94 + t65;
t244 = pkin(1) * t354;
t243 = t144 * t198 - t173;
t124 = t198 * pkin(4) - t143;
t241 = t171 + t284;
t240 = t129 + t273;
t238 = qJD(3) + t284;
t234 = qJ(6) * t175;
t233 = -pkin(7) - t265;
t59 = -pkin(4) * t131 - t85;
t232 = t263 * t196;
t228 = pkin(7) * t131 + t308;
t225 = -t194 * t77 + t197 * t69;
t223 = -qJ(4) * t194 - t197 * t337;
t219 = t198 * t238;
t217 = -t151 * t279 + t294;
t216 = t11 * t160 + t246;
t215 = t194 * t34 + t197 * t36 + t69 * t275 - t276 * t77;
t214 = -t171 * t79 + t235;
t168 = qJ(4) * t302;
t91 = t196 * t233 + t168;
t210 = -t316 * t160 - t246;
t207 = -t198 * t236 + t217;
t206 = -t131 * t222 - t175 * t197 - t194 * t338;
t203 = -t171 * t78 - t207;
t43 = (-t337 * t198 - t186) * t278 + t233 * t281 + t290;
t21 = -pkin(4) * t249 + (-pkin(4) * t209 - t266) * qJD(1) + t268;
t138 = pkin(5) - t223;
t137 = -qJ(6) + t287;
t111 = t133 * t196;
t110 = t194 * t302 - t196 * t304;
t104 = -t168 + (pkin(3) * t195 + pkin(7)) * t196;
t95 = t190 - t243;
t84 = qJD(1) * t232 - t309;
t83 = -pkin(7) * t260 + t293;
t63 = pkin(3) * t171 + t296;
t61 = pkin(5) * t133 + qJ(6) * t220 + t124;
t58 = -t96 - t311;
t54 = pkin(3) * t209 + qJ(4) * t255 + t266 - t290;
t52 = qJD(2) * t232 - t229;
t50 = -qJD(4) * t199 + (-t254 - t257) * pkin(7) + t264;
t45 = qJD(2) * t213 + (qJD(3) - qJD(5)) * t196 * t220;
t44 = t194 * t256 + t81 * t196 - t197 * t258;
t41 = pkin(5) * t110 - qJ(6) * t111 + t91;
t30 = t207 + t344;
t29 = -pkin(5) * t199 - t225;
t28 = qJ(6) * t199 + t224;
t13 = -t37 + t59;
t7 = -pkin(5) * t160 + t295;
t6 = pkin(5) * t44 - qJ(6) * t45 - qJD(6) * t111 + t43;
t5 = pkin(5) * t282 - t341;
t4 = -qJ(6) * t282 + qJD(6) * t199 + t215;
t3 = t27 * pkin(5) + t26 * qJ(6) - t75 * qJD(6) + t21;
t1 = qJD(6) * t160 - t234 - t247;
t9 = [0, 0, 0, 0.2e1 * t199 * t175, t286 * t354, t297, -t298, 0, -pkin(7) * t297 + t196 * t244, pkin(7) * t298 + t199 * t244, -t96 * t302 + (-t255 + t256) * t131 (-t129 * t281 + (-t131 - t260) * t278) * t198 + ((t96 + t280) * t196 + (-t131 * t199 - t196 * t219) * qJD(2)) * t195, t171 * t255 + t199 * t96 + (t131 * t196 + (t274 - t300) * t198) * qJD(2), t241 * t253 + (-t129 * t196 + (-t274 + (t171 + t238) * t199) * t195) * qJD(2), -t241 * t282, -t229 * t171 + (-pkin(7) * t307 + qJD(1) * t243 + t78) * t282 + t228 * t278 + ((t150 * t195 + (t129 + 0.2e1 * t262) * pkin(7)) * qJD(2) - t235) * t199 (-pkin(7) * t254 + t292) * t171 + t217 * t199 + (-pkin(7) * t96 - t150 * t279) * t196 + (t228 * t199 + (-pkin(7) * t306 - qJD(1) * t289 - t79) * t196) * qJD(2), t54 * t129 + t52 * t171 + t38 * t199 + (t320 + (qJD(1) * t104 + t68) * t277) * t196 + ((-qJD(1) * t95 - t63) * t196 + (t104 * t238 + t318) * t195) * qJD(2), -t50 * t129 + t52 * t131 - t95 * t96 + (t63 * t281 + (-qJD(3) * t245 + t38) * t196) * t198 + ((-t63 * qJD(3) - t30) * t196 + (-t65 * t199 - t238 * t94) * qJD(2)) * t195, t104 * t96 - t131 * t54 - t171 * t50 + (-t273 * t68 - t30) * t199 + (qJD(2) * t245 + t279 * t68 - t319) * t196, t104 * t35 + t30 * t94 + t38 * t95 + t50 * t65 + t52 * t63 + t54 * t68, -t111 * t26 + t45 * t75, t110 * t26 - t111 * t27 - t222 * t45 - t44 * t75, t160 * t45 - t199 * t26 + (-qJD(1) * t111 - t75) * t282, -t160 * t44 - t199 * t27 + (qJD(1) * t110 + t222) * t282 (-t160 - t284) * t282, t341 * t160 + t246 * t199 + t43 * t222 + t91 * t27 + t21 * t110 + t53 * t44 + (-qJD(1) * t225 - t10) * t282, -t215 * t160 + t247 * t199 + t43 * t75 - t91 * t26 + t21 * t111 + t53 * t45 + (t224 * qJD(1) + t11) * t282, t110 * t3 + t12 * t44 - t160 * t5 - t199 * t2 + t27 * t41 + t6 * t222 + (qJD(1) * t29 + t7) * t282, -t1 * t110 + t111 * t2 - t222 * t4 - t26 * t29 - t27 * t28 - t44 * t8 + t45 * t7 + t5 * t75, t1 * t199 - t111 * t3 - t12 * t45 + t160 * t4 + t26 * t41 - t6 * t75 + (-qJD(1) * t28 - t8) * t282, t1 * t28 + t12 * t6 + t2 * t29 + t3 * t41 + t4 * t8 + t5 * t7; 0, 0, 0, -t196 * t299, t286 * t202, 0, 0, 0, t202 * pkin(1) * t196, pkin(1) * t299, -t131 * t306 - t96 * t195 (-t96 + t311) * t198 + ((-t131 - t283) * qJD(3) + (-t253 - t345) * qJD(1)) * t195, -t171 * t277 + (t196 * t239 + t198 * t300) * qJD(1), t171 * t279 + (-t195 * t300 + t196 * t240) * qJD(1), t171 * t285, t139 * t306 + (t269 + (t150 - t193) * t195) * qJD(3) + ((-pkin(2) * t277 - t78) * t196 + (qJD(2) * t231 - t150 * t199) * t195 + (t171 * t305 - t199 * t240) * pkin(7)) * qJD(1), pkin(2) * t96 - t115 * t171 + (-t270 + t308) * qJD(3) + (-t150 * t301 + (t79 - t267) * t196 + (t171 * t302 + t345) * pkin(7)) * qJD(1), -t84 * t171 - t319 + t315 * t129 + (t269 + (qJD(2) * t143 + t68) * t195) * qJD(3) + ((t143 * t277 + t63) * t196 + (-t318 + (t143 * t199 - t335) * qJD(2)) * t195) * qJD(1), t83 * t129 - t84 * t131 + (-t63 * t284 + t30 + (t63 + (t131 - t260) * pkin(8)) * qJD(3)) * t198 + ((-qJD(2) * t219 + t280 - t96) * pkin(8) + t348) * t195, t143 * t96 + t171 * t83 - t320 - t315 * t131 + (-t198 * t68 + t270) * qJD(3) + (t68 * t301 + (-t65 + t267) * t196) * qJD(1), t143 * t35 - t63 * t84 - t65 * t83 + t315 * t68 + (t38 * t195 + t30 * t198 + (-t65 * t195 + t63 * t198) * qJD(3)) * pkin(8), t220 * t26 + t314 * t75, t133 * t26 + t220 * t27 - t222 * t314 - t327 * t75, t314 * t160 + (qJD(2) * t220 + t75) * t285, -t327 * t160 + (qJD(2) * t133 - t222) * t285, t160 * t285, t124 * t27 + t21 * t133 + t326 * t222 + t327 * t53 + t351 * t160 + (t10 - t313) * t285, -t124 * t26 - t21 * t220 + t326 * t75 + t314 * t53 + t347 * t160 + (-t11 + t312) * t285, t133 * t3 + t27 * t61 + t331 * t222 + t329 * t160 + t327 * t12 + (-t7 - t313) * t285, -t1 * t133 - t2 * t220 + t221 * t26 + t222 * t330 - t27 * t93 + t314 * t7 - t327 * t8 - t329 * t75, t220 * t3 + t26 * t61 - t331 * t75 - t330 * t160 - t314 * t12 + (t8 - t312) * t285, t1 * t93 + t331 * t12 - t2 * t221 + t3 * t61 - t329 * t7 - t330 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t310, -t129 ^ 2 + t339, t58, -t131 * t171 - t204, t175, -t131 * t150 + t214, t129 * t150 + t203, -t129 * t85 + t214 + 0.2e1 * t237 - t324, -t204 * qJ(4) + pkin(3) * t96 + (t65 - t79) * t131 + (t63 - t296) * t129, -t129 * t68 + t131 * t85 - 0.2e1 * t154 + 0.2e1 * t162 - t203, -pkin(3) * t38 + qJ(4) * t30 + t296 * t65 - t63 * t79 - t68 * t85, -t332, -t360, t350, -t349, t175, -t175 * t223 - t222 * t59 + t210 + t359, t352 * t160 + t287 * t175 - t59 * t75 + t342, -t13 * t222 + t138 * t175 + t170 + t210 + t334, -t137 * t27 - t138 * t26 + (-t317 - t7) * t222 + (t316 - t8) * t75, t13 * t75 + (qJ(6) - t137) * t175 + (-qJD(6) + t317) * t160 - t343, t1 * t137 - t12 * t13 + t138 * t2 + t316 * t7 + t317 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t175 + t310, t58, -t171 ^ 2 - t339, t324 + t348, 0, 0, 0, 0, 0, t206, -t356, t206, t349 * t194 + t350 * t197, t356, -t12 * t131 + t353 * t197 + (t160 * t7 + t1) * t194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t332, t360, -t350, t349, -t175, t216 - t359, t10 * t160 - t342, -t222 * t37 - 0.2e1 * t170 + t216 - t334, pkin(5) * t26 - qJ(6) * t27 + (-t11 + t8) * t75 + (t7 - t295) * t222, -0.2e1 * t234 + t37 * t75 + (0.2e1 * qJD(6) - t10) * t160 + t343, -pkin(5) * t2 + qJ(6) * t1 - t11 * t7 - t12 * t37 + t295 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t175 + t332, -t350, -t338 - t340, t334 - t353;];
tauc_reg  = t9;
