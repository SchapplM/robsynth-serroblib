% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% tauc_reg [6x35]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRRR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:09:50
% EndTime: 2019-03-09 07:10:05
% DurationCPUTime: 6.12s
% Computational Cost: add. (9139->396), mult. (23928->535), div. (0->0), fcn. (19856->10), ass. (0->227)
t203 = cos(pkin(11));
t210 = cos(qJ(3));
t278 = t203 * t210;
t202 = sin(pkin(11));
t207 = sin(qJ(3));
t279 = t202 * t207;
t229 = -t278 + t279;
t165 = t229 * qJD(1);
t173 = t202 * t210 + t203 * t207;
t166 = t173 * qJD(1);
t206 = sin(qJ(4));
t313 = cos(qJ(4));
t143 = t313 * t165 + t166 * t206;
t322 = qJD(5) + qJD(6);
t368 = t143 + t322;
t265 = qJD(1) * qJD(3);
t252 = t210 * t265;
t189 = t203 * t252;
t253 = t207 * t265;
t162 = -t202 * t253 + t189;
t224 = -t206 * t165 + t313 * t166;
t276 = t202 * t252 + t203 * t253;
t100 = qJD(4) * t224 + t206 * t162 + t313 * t276;
t339 = -qJD(5) - t143;
t138 = qJD(6) - t339;
t205 = sin(qJ(5));
t208 = cos(qJ(6));
t204 = sin(qJ(6));
t209 = cos(qJ(5));
t277 = t204 * t209;
t177 = t205 * t208 + t277;
t176 = t204 * t205 - t208 * t209;
t366 = t368 * t176;
t367 = t177 * t100 - t138 * t366;
t365 = t368 * t177;
t271 = qJD(5) * t209;
t346 = t143 * t209;
t364 = t271 + t346;
t347 = t143 * t205;
t363 = pkin(10) * t347;
t254 = qJD(4) * t313;
t311 = pkin(7) + qJ(2);
t182 = t311 * t202;
t174 = qJD(1) * t182;
t183 = t311 * t203;
t175 = qJD(1) * t183;
t231 = t174 * t207 - t175 * t210;
t130 = -pkin(8) * t165 - t231;
t124 = t206 * t130;
t325 = -t210 * t174 - t175 * t207;
t129 = -pkin(8) * t166 + t325;
t74 = t313 * t129 - t124;
t349 = pkin(3) * t254 - t74;
t239 = -t176 * t100 - t138 * t365;
t201 = qJD(3) + qJD(4);
t131 = -t209 * t201 + t205 * t224;
t133 = t201 * t205 + t209 * t224;
t76 = t208 * t131 + t133 * t204;
t305 = t224 * t76;
t362 = t239 + t305;
t272 = qJD(5) * t205;
t361 = (t272 + t347) * pkin(5);
t360 = pkin(5) * t224 + pkin(10) * t346;
t104 = pkin(4) * t224 + pkin(9) * t143;
t82 = pkin(3) * t166 + t104;
t359 = -t349 * t205 - t209 * t82;
t111 = -t276 * pkin(8) - qJD(2) * t165 + t325 * qJD(3);
t220 = t173 * qJD(2);
t219 = qJD(1) * t220;
t112 = -pkin(8) * t162 + qJD(3) * t231 - t219;
t126 = qJD(3) * pkin(3) + t129;
t273 = qJD(4) * t206;
t28 = t206 * t111 - t313 * t112 + t126 * t273 + t130 * t254;
t99 = t313 * t162 - t165 * t254 - t166 * t273 - t206 * t276;
t303 = t205 * t99;
t56 = qJD(5) * t133 + t303;
t13 = pkin(5) * t56 + t28;
t70 = t313 * t126 - t124;
t64 = -t201 * pkin(4) - t70;
t46 = t131 * pkin(5) + t64;
t125 = t313 * t130;
t71 = t206 * t126 + t125;
t65 = pkin(9) * t201 + t71;
t193 = -pkin(2) * t203 - pkin(1);
t181 = t193 * qJD(1) + qJD(2);
t151 = pkin(3) * t165 + t181;
t72 = pkin(4) * t143 - pkin(9) * t224 + t151;
t32 = -t205 * t65 + t209 * t72;
t22 = -pkin(10) * t133 + t32;
t17 = -pkin(5) * t339 + t22;
t33 = t205 * t72 + t209 * t65;
t23 = -pkin(10) * t131 + t33;
t302 = t208 * t23;
t7 = t204 * t17 + t302;
t358 = t13 * t177 + t7 * t224 - t366 * t46;
t6 = t208 * t17 - t204 * t23;
t357 = t13 * t176 - t6 * t224 + t365 * t46;
t55 = t201 * t271 + t209 * t99 - t224 * t272;
t356 = -t131 * t364 - t205 * t56 + t55 * t209;
t53 = t55 * t205;
t355 = t133 * t364 + t53;
t232 = t131 * t204 - t208 * t133;
t304 = t224 * t232;
t354 = t304 + t367;
t269 = qJD(6) * t208;
t270 = qJD(6) * t204;
t14 = -t131 * t269 - t133 * t270 - t204 * t56 + t208 * t55;
t353 = t14 * t177 + t232 * t366;
t291 = t133 * t224;
t96 = t205 * t100;
t297 = -t271 * t339 + t96;
t352 = -t339 * t346 - t291 + t297;
t215 = qJD(6) * t232 - t204 * t55 - t208 * t56;
t351 = -t14 * t176 + t177 * t215 + t365 * t232 + t366 * t76;
t350 = t232 * t76;
t286 = t143 * t201;
t348 = t99 + t286;
t343 = t224 * t143;
t289 = t224 * t201;
t342 = -t100 + t289;
t338 = -t143 ^ 2 + t224 ^ 2;
t337 = t232 ^ 2 - t76 ^ 2;
t21 = t23 * t270;
t336 = t46 * t76 + t21;
t335 = t138 * t76 + t14;
t213 = -t313 * t111 - t206 * t112 - t126 * t254 + t130 * t273;
t42 = t276 * pkin(3) + t100 * pkin(4) - t99 * pkin(9);
t39 = t209 * t42;
t216 = -qJD(5) * t33 + t205 * t213 + t39;
t2 = pkin(5) * t100 - pkin(10) * t55 + t216;
t226 = t205 * t42 - t209 * t213 + t72 * t271 - t65 * t272;
t3 = -pkin(10) * t56 + t226;
t261 = t208 * t2 - t204 * t3;
t334 = -qJD(6) * t7 + t46 * t232 + t261;
t333 = -t138 * t232 + t215;
t332 = t151 * t143 + t213;
t98 = t209 * t100;
t328 = -t272 * t339 - t98;
t249 = -t28 * t209 + t64 * t272;
t73 = t206 * t129 + t125;
t241 = pkin(3) * t273 - t73;
t292 = t131 * t224;
t327 = t138 * t224;
t326 = t339 * t224;
t148 = t313 * t173 - t206 * t229;
t105 = t177 * t148;
t280 = t182 * t210;
t136 = -pkin(8) * t173 - t183 * t207 - t280;
t230 = t182 * t207 - t183 * t210;
t137 = -pkin(8) * t229 - t230;
t324 = t313 * t136 - t206 * t137;
t323 = t205 * t82 - t349 * t209;
t319 = -t224 * t151 - t28;
t318 = -t32 * t224 + t249;
t317 = t28 * t205 + t33 * t224 + t64 * t271;
t316 = -pkin(9) - pkin(10);
t312 = t209 * pkin(5);
t195 = pkin(3) * t206 + pkin(9);
t310 = -pkin(10) - t195;
t90 = t206 * t136 + t313 * t137;
t86 = t209 * t90;
t157 = pkin(3) * t229 + t193;
t223 = -t206 * t173 - t229 * t313;
t91 = -pkin(4) * t223 - pkin(9) * t148 + t157;
t306 = t205 * t91 + t86;
t296 = t241 + t361;
t295 = t205 * t104 + t209 * t70;
t167 = t229 * qJD(3);
t168 = t173 * qJD(3);
t109 = t223 * qJD(4) - t313 * t167 - t206 * t168;
t294 = t109 * t205;
t293 = t109 * t209;
t290 = t133 * t205;
t283 = t148 * t205;
t282 = t148 * t209;
t275 = t202 ^ 2 + t203 ^ 2;
t274 = qJD(3) * t166;
t266 = qJD(1) * qJD(2);
t260 = qJD(5) * t316;
t258 = qJD(1) * t279;
t257 = t148 * t272;
t256 = t148 * t271;
t251 = qJD(6) * t17 + t3;
t248 = qJD(5) * t310;
t246 = t209 * t104 - t205 * t70;
t245 = t275 * qJD(1) ^ 2;
t243 = t339 * t205;
t196 = -t313 * pkin(3) - pkin(4);
t240 = t361 - t71;
t170 = t310 * t205;
t238 = -qJD(6) * t170 - t205 * t248 + t323 + t363;
t198 = t209 * pkin(10);
t171 = t195 * t209 + t198;
t237 = qJD(6) * t171 - t209 * t248 - t359 + t360;
t185 = t316 * t205;
t236 = -qJD(6) * t185 - t205 * t260 + t295 + t363;
t186 = pkin(9) * t209 + t198;
t235 = qJD(6) * t186 - t209 * t260 + t246 + t360;
t233 = -t195 * t100 + t143 * t64;
t228 = 0.2e1 * t275 * t266;
t227 = t339 * t347 - t328;
t218 = -qJD(3) * t280 + qJD(2) * t278 + (-qJD(2) * t202 - qJD(3) * t183) * t207;
t115 = -pkin(8) * t168 + t218;
t212 = qJD(3) * t230 - t220;
t116 = pkin(8) * t167 + t212;
t36 = t324 * qJD(4) + t313 * t115 + t206 * t116;
t110 = t148 * qJD(4) - t206 * t167 + t313 * t168;
t45 = pkin(3) * t168 + pkin(4) * t110 - pkin(9) * t109;
t225 = t205 * t45 + t209 * t36 + t91 * t271 - t90 * t272;
t222 = t256 + t294;
t221 = -t257 + t293;
t37 = t90 * qJD(4) + t206 * t115 - t313 * t116;
t197 = -pkin(4) - t312;
t184 = t196 - t312;
t106 = t176 * t148;
t88 = t209 * t91;
t63 = t100 * t223;
t60 = pkin(5) * t283 - t324;
t44 = t209 * t45;
t35 = -pkin(10) * t283 + t306;
t29 = -pkin(5) * t223 - pkin(10) * t282 - t205 * t90 + t88;
t19 = t109 * t277 - t204 * t257 - t270 * t283 + (t322 * t282 + t294) * t208;
t18 = -t322 * t105 - t176 * t109;
t16 = pkin(5) * t222 + t37;
t5 = -pkin(10) * t222 + t225;
t4 = -pkin(10) * t293 + pkin(5) * t110 - t205 * t36 + t44 + (-t86 + (pkin(10) * t148 - t91) * t205) * qJD(5);
t1 = [0, 0, 0, 0, 0, t228, qJ(2) * t228, t162 * t173 - t166 * t167, -t162 * t229 + t167 * t165 - t166 * t168 - t173 * t276, -t167 * qJD(3), -t168 * qJD(3), 0, t212 * qJD(3) + t181 * t168 + t193 * t276, -qJD(3) * t218 + t193 * t162 - t181 * t167, t109 * t224 + t148 * t99, -t100 * t148 - t109 * t143 - t110 * t224 + t223 * t99, t109 * t201, -t110 * t201, 0, t157 * t100 + t151 * t110 - t37 * t201 + (t168 * t143 - t223 * t276) * pkin(3), t151 * t109 + t157 * t99 - t36 * t201 + (t148 * t276 + t168 * t224) * pkin(3), t133 * t221 + t55 * t282 (-t131 * t209 - t290) * t109 + (-t53 - t209 * t56 + (t131 * t205 - t133 * t209) * qJD(5)) * t148, t110 * t133 + t148 * t98 - t221 * t339 - t223 * t55, -t110 * t131 - t148 * t96 + t222 * t339 + t223 * t56, -t110 * t339 - t63 -(-t271 * t90 + t44) * t339 + t88 * t100 - (-t271 * t65 + t39) * t223 + t32 * t110 + t37 * t131 - t324 * t56 + t64 * t256 + (-(-qJD(5) * t91 - t36) * t339 - t90 * t100 - (-qJD(5) * t72 + t213) * t223 + t28 * t148 + t64 * t109) * t205, -t306 * t100 - t33 * t110 + t37 * t133 - t249 * t148 + t223 * t226 + t225 * t339 + t64 * t293 - t324 * t55, -t106 * t14 - t18 * t232, -t105 * t14 - t106 * t215 - t18 * t76 + t19 * t232, -t100 * t106 - t110 * t232 + t138 * t18 - t14 * t223, -t100 * t105 - t110 * t76 - t138 * t19 - t215 * t223, t110 * t138 - t63 (-t204 * t5 + t208 * t4) * t138 + (-t204 * t35 + t208 * t29) * t100 - t261 * t223 + t6 * t110 + t16 * t76 - t60 * t215 + t13 * t105 + t46 * t19 + ((-t204 * t29 - t208 * t35) * t138 + t7 * t223) * qJD(6), -t13 * t106 - t7 * t110 + t60 * t14 - t21 * t223 - t16 * t232 + t46 * t18 + (-(-qJD(6) * t35 + t4) * t138 - t29 * t100 + t2 * t223) * t204 + (-(qJD(6) * t29 + t5) * t138 - t35 * t100 + t251 * t223) * t208; 0, 0, 0, 0, 0, -t245, -qJ(2) * t245, 0, 0, 0, 0, 0, t274 + t276, t189 + (-t165 - t258) * qJD(3), 0, 0, 0, 0, 0, t100 + t289, t99 - t286, 0, 0, 0, 0, 0, t227 - t292, -t209 * t339 ^ 2 - t291 - t96, 0, 0, 0, 0, 0, t239 - t305, t304 - t367; 0, 0, 0, 0, 0, 0, 0, t166 * t165, -t165 ^ 2 + t166 ^ 2, t189 + (t165 - t258) * qJD(3), t274 - t276, 0, -t181 * t166 - t219, t181 * t165 + t229 * t266, t343, t338, t348, t342, 0, t201 * t73 + (-t143 * t166 - t201 * t273) * pkin(3) + t319, t74 * t201 + (-t166 * t224 - t201 * t254) * pkin(3) + t332, t355, t290 * t339 + t356, t352, t227 + t292, t326, t196 * t56 + t233 * t205 + t241 * t131 - (-t195 * t271 + t359) * t339 + t318, t196 * t55 + t233 * t209 + t241 * t133 - (t195 * t272 + t323) * t339 + t317, t353, t351, t354, t362, -t327 (t170 * t208 - t171 * t204) * t100 - t184 * t215 + t296 * t76 + (t204 * t238 - t208 * t237) * t138 + t357 -(t170 * t204 + t171 * t208) * t100 + t184 * t14 - t296 * t232 + (t204 * t237 + t208 * t238) * t138 + t358; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t343, t338, t348, t342, 0, t201 * t71 + t319, t70 * t201 + t332, t355, t133 * t243 + t356, t352, -t243 * t339 + t292 + t98, t326, -pkin(4) * t56 - t297 * pkin(9) - t71 * t131 + t246 * t339 + t347 * t64 + t318, -pkin(4) * t55 + t328 * pkin(9) - t71 * t133 - t295 * t339 + t346 * t64 + t317, t353, t351, t354, t362, -t327 (t185 * t208 - t186 * t204) * t100 - t197 * t215 + t240 * t76 + (t204 * t236 - t208 * t235) * t138 + t357 -(t185 * t204 + t186 * t208) * t100 + t197 * t14 - t240 * t232 + (t204 * t235 + t208 * t236) * t138 + t358; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133 * t131, -t131 ^ 2 + t133 ^ 2, -t131 * t339 + t55, -t303 + (-qJD(5) - t339) * t133, t100, -t133 * t64 - t33 * t339 + t216, t131 * t64 - t32 * t339 - t226, -t350, t337, t335, t333, t100 -(-t204 * t22 - t302) * t138 + (t208 * t100 - t133 * t76 - t138 * t270) * pkin(5) + t334 (-t23 * t138 - t2) * t204 + (t22 * t138 - t251) * t208 + (-t204 * t100 + t133 * t232 - t138 * t269) * pkin(5) + t336; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t350, t337, t335, t333, t100, t138 * t7 + t334, t138 * t6 - t204 * t2 - t208 * t251 + t336;];
tauc_reg  = t1;
