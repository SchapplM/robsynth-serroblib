% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPRRRR7
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
% Datum: 2019-03-09 07:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRRR7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:17:55
% EndTime: 2019-03-09 07:18:11
% DurationCPUTime: 6.43s
% Computational Cost: add. (14858->483), mult. (31368->628), div. (0->0), fcn. (21883->8), ass. (0->248)
t195 = sin(qJ(5));
t189 = qJD(3) + qJD(4);
t196 = sin(qJ(4));
t197 = sin(qJ(3));
t199 = cos(qJ(3));
t329 = cos(qJ(4));
t337 = t196 * t199 + t329 * t197;
t352 = t189 * t337;
t108 = t352 * qJD(1);
t200 = -pkin(1) - pkin(7);
t168 = t200 * qJD(1) + qJD(2);
t249 = pkin(8) * qJD(1) - t168;
t275 = qJD(3) * t197;
t131 = t249 * t275;
t274 = qJD(3) * t199;
t132 = t249 * t274;
t244 = t329 * t131 + t196 * t132;
t276 = qJD(1) * t199;
t138 = -pkin(8) * t276 + t199 * t168;
t130 = qJD(3) * pkin(3) + t138;
t277 = qJD(1) * t197;
t137 = -pkin(8) * t277 + t168 * t197;
t260 = t329 * t137;
t92 = t196 * t130 + t260;
t67 = -qJD(4) * t92 + t244;
t209 = t108 * pkin(9) + t67;
t328 = cos(qJ(5));
t256 = qJD(5) * t328;
t272 = qJD(5) * t195;
t254 = t329 * qJD(4);
t215 = t329 * qJD(3) + t254;
t258 = t196 * t277;
t109 = -t189 * t258 + t215 * t276;
t273 = qJD(4) * t196;
t238 = -t130 * t254 - t196 * t131 + t329 * t132 + t137 * t273;
t47 = -pkin(9) * t109 - t238;
t257 = qJD(1) * t329;
t150 = t199 * t257 - t258;
t142 = t150 * pkin(9);
t128 = t196 * t137;
t91 = t329 * t130 - t128;
t78 = -t142 + t91;
t72 = pkin(4) * t189 + t78;
t149 = -t196 * t276 - t197 * t257;
t325 = t149 * pkin(9);
t79 = t92 + t325;
t205 = -t195 * t209 - t72 * t256 + t79 * t272 - t328 * t47;
t157 = -t196 * t197 + t329 * t199;
t339 = t195 * t337;
t114 = -t328 * t157 + t339;
t307 = t195 * t47;
t46 = t328 * t209;
t251 = -t46 + t307;
t264 = t328 * t79;
t45 = t195 * t72 + t264;
t9 = qJD(5) * t45 + t251;
t326 = t114 * t9;
t234 = t328 * t337;
t351 = t195 * t157 + t234;
t281 = t215 * t199;
t340 = t189 * t197;
t349 = -t196 * t340 + t281;
t356 = -qJD(5) * t339 + t157 * t256 - t195 * t352 + t328 * t349;
t364 = -t205 * t351 + t356 * t45 + t326;
t210 = -t328 * t108 - t195 * t109 + t149 * t256 - t150 * t272;
t363 = t114 * t210;
t250 = qJD(5) + t189;
t362 = t356 * t250;
t102 = -t328 * t149 + t150 * t195;
t220 = t195 * t149 + t328 * t150;
t245 = t195 * t108 - t328 * t109;
t56 = qJD(5) * t220 - t245;
t315 = t351 * t56;
t360 = t102 * t356 + t315;
t358 = t349 * t189;
t198 = cos(qJ(6));
t270 = qJD(6) * t198;
t194 = sin(qJ(6));
t226 = t198 * t250;
t271 = qJD(6) * t194;
t35 = -qJD(6) * t226 - t198 * t210 + t220 * t271;
t32 = t35 * t194;
t354 = t102 * t198;
t87 = t194 * t250 + t198 * t220;
t12 = -t32 + (t270 + t354) * t87;
t357 = t195 * t349 + t328 * t352;
t350 = qJD(6) + t102;
t52 = t194 * t56;
t319 = t270 * t350 + t52;
t11 = -t220 * t87 + t350 * t354 + t319;
t355 = t67 * t157 - t238 * t337 + t349 * t92 - t352 * t91;
t353 = t102 * t220;
t39 = -t102 ^ 2 + t220 ^ 2;
t71 = pkin(5) * t220 + pkin(10) * t102;
t162 = pkin(3) * t277 + qJD(1) * qJ(2);
t120 = -pkin(4) * t149 + t162;
t204 = t120 * t102 + t205;
t37 = t102 * t250 + t210;
t85 = t194 * t220 - t226;
t303 = t198 * t85;
t309 = t194 * t87;
t228 = t303 + t309;
t343 = -t35 * t198 - t87 * t271;
t297 = qJD(6) * t87;
t36 = t194 * t210 + t297;
t344 = t194 * t36 + t85 * t270;
t6 = -t102 * t228 + t343 - t344;
t43 = pkin(10) * t250 + t45;
t57 = pkin(5) * t102 - pkin(10) * t220 + t120;
t21 = t194 * t57 + t198 * t43;
t190 = qJD(1) * qJD(2);
t267 = qJD(1) * qJD(3);
t253 = t199 * t267;
t159 = pkin(3) * t253 + t190;
t88 = pkin(4) * t109 + t159;
t19 = pkin(5) * t56 - pkin(10) * t210 + t88;
t3 = -qJD(6) * t21 + t198 * t19 + t194 * t205;
t345 = -t21 * t350 - t3;
t300 = t350 * t220;
t182 = t329 * pkin(3) + pkin(4);
t286 = t195 * t196;
t106 = t182 * t256 + (-t196 * t272 + (t328 * t329 - t286) * qJD(4)) * pkin(3);
t96 = -t196 * t138 - t260;
t214 = t96 - t325;
t97 = t329 * t138 - t128;
t80 = -t142 + t97;
t51 = t195 * t214 + t328 * t80;
t299 = t106 - t51;
t259 = t328 * t196;
t298 = -t195 * t80 + t328 * t214 + t182 * t272 + (t196 * t256 + (t329 * t195 + t259) * qJD(4)) * pkin(3);
t341 = t120 * t220;
t338 = t337 * t109;
t320 = pkin(8) - t200;
t160 = t320 * t197;
t161 = t320 * t199;
t119 = -t329 * t160 - t196 * t161;
t229 = t194 * t43 - t198 * t57;
t336 = t194 * t229 + t198 * t21;
t306 = t195 * t79;
t44 = t328 * t72 - t306;
t42 = -pkin(5) * t250 - t44;
t40 = t42 * t271;
t252 = t220 * t229 + t40;
t54 = t198 * t56;
t335 = t220 * t85 + t54;
t330 = t9 * t194 + t42 * t270;
t237 = t21 * t220 + t330;
t38 = t220 * t189 + t245;
t2 = -qJD(6) * t229 + t19 * t194 - t198 * t205;
t332 = 0.2e1 * t190;
t118 = t160 * t196 - t329 * t161;
t94 = -pkin(9) * t157 + t118;
t95 = -pkin(9) * t337 + t119;
t68 = t195 * t95 - t328 * t94;
t331 = t68 * t9;
t327 = pkin(4) * t150;
t1 = t2 * t198;
t324 = t229 * t350;
t322 = t87 * t85;
t321 = t9 * t198;
t317 = pkin(4) * qJD(5);
t316 = t102 * t42;
t313 = t194 * t21;
t312 = t194 * t42;
t310 = t194 * t85;
t308 = t194 * t350;
t60 = qJD(5) * t234 + t157 * t272 + t357;
t304 = t198 * t60;
t302 = t198 * t87;
t34 = t36 * t198;
t296 = qJD(6) * t350;
t290 = t150 * t149;
t289 = t157 * t108;
t288 = t162 * t150;
t201 = qJD(3) ^ 2;
t285 = t201 * t197;
t284 = t201 * t199;
t202 = qJD(1) ^ 2;
t283 = t202 * qJ(2);
t177 = t197 * pkin(3) + qJ(2);
t146 = pkin(3) * t259 + t195 * t182;
t279 = t197 ^ 2 - t199 ^ 2;
t278 = -t201 - t202;
t269 = t162 * qJD(1);
t169 = pkin(3) * t274 + qJD(2);
t266 = -t102 * t313 + t229 * t354 + t1;
t265 = 0.2e1 * qJD(1);
t185 = pkin(3) * t276;
t262 = t199 * t202 * t197;
t155 = t320 * t275;
t156 = qJD(3) * t161;
t261 = t196 * t155 - t329 * t156 - t161 * t254;
t248 = t46 - t341;
t243 = t199 * t189;
t240 = qJD(6) * t351 + qJD(1);
t239 = pkin(4) * t256;
t236 = t197 * t253;
t48 = t195 * t78 + t264;
t235 = pkin(4) * t272 - t48;
t232 = -t102 * t44 + t220 * t45;
t231 = -t198 * t229 + t313;
t69 = t195 * t94 + t328 * t95;
t133 = pkin(4) * t337 + t177;
t70 = pkin(5) * t351 + pkin(10) * t114 + t133;
t30 = t194 * t70 + t198 * t69;
t29 = -t194 * t69 + t198 * t70;
t227 = t302 + t310;
t224 = t271 * t350 - t54;
t222 = t329 * t155 + t196 * t156;
t219 = -t162 * t149 + t238;
t64 = t327 + t71;
t145 = -pkin(3) * t286 + t328 * t182;
t140 = pkin(10) + t146;
t218 = -t106 * t350 - t140 * t56 + t316;
t213 = t197 * t215;
t180 = t195 * pkin(4) + pkin(10);
t212 = -t180 * t56 - t239 * t350 + t316;
t99 = pkin(4) * t349 + t169;
t211 = -t231 * qJD(6) - t3 * t194 + t1;
t203 = pkin(9) * t352 + t160 * t254 + t161 * t273 + t222;
t187 = qJ(2) * t332;
t181 = -t328 * pkin(4) - pkin(5);
t139 = -pkin(5) - t145;
t126 = t185 + t327;
t110 = t352 * t189;
t89 = -t149 ^ 2 + t150 ^ 2;
t82 = t150 * t189 - t109;
t81 = -t149 * t189 - t108;
t76 = -t119 * qJD(4) + t222;
t75 = t160 * t273 + t261;
t65 = -t281 * pkin(9) + (pkin(9) * t340 + qJD(4) * t160) * t196 + t261;
t61 = qJD(5) * t351 + t357;
t58 = t185 + t64;
t49 = t328 * t78 - t306;
t28 = t194 * t71 + t198 * t44;
t27 = -t194 * t44 + t198 * t71;
t26 = t194 * t58 + t198 * t51;
t25 = -t194 * t51 + t198 * t58;
t24 = t194 * t64 + t198 * t49;
t23 = -t194 * t49 + t198 * t64;
t22 = pkin(5) * t356 + t60 * pkin(10) + t99;
t17 = t69 * qJD(5) + t195 * t65 - t328 * t203;
t16 = t195 * t203 + t94 * t256 - t95 * t272 + t328 * t65;
t13 = t310 * t350 - t34;
t10 = -t308 * t350 + t335;
t5 = -qJD(6) * t30 - t16 * t194 + t198 * t22;
t4 = qJD(6) * t29 + t16 * t198 + t194 * t22;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t332, t187, -0.2e1 * t236, 0.2e1 * t279 * t267, -t285, 0.2e1 * t236, -t284, 0, -t200 * t285 + (qJ(2) * t274 + qJD(2) * t197) * t265, -t200 * t284 + (-qJ(2) * t275 + qJD(2) * t199) * t265, 0, t187, -t150 * t352 - t289, -t157 * t109 + t108 * t337 - t150 * t281 - t149 * t213 + (-t149 * t243 + t150 * t340) * t196, -t110, -t149 * t349 + t338, -t358, 0, t177 * t109 - t169 * t149 + t159 * t337 + t162 * t349 + t76 * t189, -t177 * t108 + t169 * t150 + t159 * t157 - t162 * t352 - t75 * t189, t118 * t108 - t119 * t109 + t75 * t149 - t76 * t150 - t355, t118 * t67 - t119 * t238 + t159 * t177 + t162 * t169 + t75 * t92 + t76 * t91, -t220 * t60 - t363, t102 * t60 + t114 * t56 - t210 * t351 - t220 * t356, -t60 * t250, t360, -t362, 0, t99 * t102 + t120 * t356 + t133 * t56 - t17 * t250 + t351 * t88, -t114 * t88 - t120 * t60 + t133 * t210 - t16 * t250 + t220 * t99, -t102 * t16 + t17 * t220 + t210 * t68 + t44 * t60 - t56 * t69 - t364, t120 * t99 + t133 * t88 + t16 * t45 - t17 * t44 - t205 * t69 + t331, -t114 * t343 - t60 * t302, t228 * t60 - (t32 - t34 + (-t302 + t310) * qJD(6)) * t114, t114 * t224 - t304 * t350 - t35 * t351 + t356 * t87, -t114 * t344 - t60 * t310, t114 * t319 + t308 * t60 - t351 * t36 - t356 * t85, t350 * t356 + t315, -t114 * t330 + t17 * t85 - t229 * t356 + t29 * t56 + t3 * t351 - t60 * t312 + t350 * t5 + t36 * t68, -t42 * t304 - t351 * t2 + t17 * t87 - t21 * t356 - t30 * t56 - t35 * t68 - t4 * t350 - (-t40 + t321) * t114, t29 * t35 - t30 * t36 - t4 * t85 - t5 * t87 + t231 * t60 - (-qJD(6) * t336 - t194 * t2 - t198 * t3) * t114, t17 * t42 + t2 * t30 + t21 * t4 - t229 * t5 + t29 * t3 + t331; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t202, -t283, 0, 0, 0, 0, 0, 0, t278 * t197, t278 * t199, 0, -t283, 0, 0, 0, 0, 0, 0, qJD(1) * t149 - t110, -qJD(1) * t150 - t358, t281 * t149 - t338 + t289 + t150 * t213 + (-t149 * t340 + t150 * t243) * t196, -t269 + t355, 0, 0, 0, 0, 0, 0, -qJD(1) * t102 - t250 * t61, -qJD(1) * t220 - t362, t220 * t61 - t360 + t363, -qJD(1) * t120 - t44 * t61 + t364, 0, 0, 0, 0, 0, 0, -t351 * t52 + t114 * t36 + t61 * t85 + (-t194 * t356 - t198 * t240) * t350, -t351 * t54 - t114 * t35 + t61 * t87 + (t194 * t240 - t198 * t356) * t350 (-t303 + t309) * t356 + t227 * qJD(1) + (qJD(6) * t227 - t32 - t34) * t351, -qJD(1) * t231 + t211 * t351 + t336 * t356 + t42 * t61 + t326; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t262, -t279 * t202, 0, -t262, 0, 0, -t199 * t283, t197 * t283, 0, 0, -t290, t89, t81, t290, t82, 0, t149 * t185 - t288 - t96 * t189 + (-t260 + (-pkin(3) * t189 - t130) * t196) * qJD(4) + t244, t97 * t189 + (-t150 * t276 - t189 * t254) * pkin(3) + t219 (t92 + t96) * t150 + (t91 - t97) * t149 + (t329 * t108 - t109 * t196 + (t329 * t149 + t150 * t196) * qJD(4)) * pkin(3), -t91 * t96 - t92 * t97 + (-t199 * t269 + t329 * t67 - t196 * t238 + (-t196 * t91 + t329 * t92) * qJD(4)) * pkin(3), t353, t39, t37, -t353, t38, 0, -t307 - t126 * t102 + (-t45 - t298) * qJD(5) + t248 - t298 * t189, -t126 * t220 - t299 * t250 + t204, -t102 * t299 - t145 * t210 - t146 * t56 + t220 * t298 + t232, -t120 * t126 - t145 * t9 - t146 * t205 - t298 * t44 + t299 * t45, t12, t6, t11, t13, t10, -t300, t139 * t36 - t25 * t350 + t298 * t85 + (-t140 * t296 - t9) * t198 + t218 * t194 + t252, -t139 * t35 + (t140 * t271 + t26) * t350 + t298 * t87 + t218 * t198 + t237, t25 * t87 + t26 * t85 + (-t106 * t85 - t140 * t36 + (t140 * t87 + t229) * qJD(6)) * t198 + (t106 * t87 - t140 * t35 - t3 + (t140 * t85 - t21) * qJD(6)) * t194 + t266, t106 * t336 + t139 * t9 + t140 * t211 - t21 * t26 + t229 * t25 + t298 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t290, t89, t81, t290, t82, 0, t92 * t189 - t288 + t67, t189 * t91 + t219, 0, 0, t353, t39, t37, -t353, t38, 0, -t79 * t256 - t102 * t327 + t48 * t250 + (-qJD(5) * t72 - t250 * t317 - t47) * t195 + t248, -t220 * t327 + t204 + (-t239 + t49) * t250, t49 * t102 - t48 * t220 + (-t328 * t210 - t195 * t56 + (-t328 * t102 + t195 * t220) * qJD(5)) * pkin(4) + t232, t44 * t48 - t45 * t49 + (-t328 * t9 - t120 * t150 - t195 * t205 + (-t195 * t44 + t328 * t45) * qJD(5)) * pkin(4), t12, t6, t11, t13, t10, -t300, t181 * t36 - t23 * t350 + t235 * t85 + (-t180 * t296 - t9) * t198 + t212 * t194 + t252, -t181 * t35 + (t180 * t271 + t24) * t350 + t235 * t87 + t212 * t198 + t237, t23 * t87 + t24 * t85 + (-t85 * t239 - t180 * t36 + (t180 * t87 + t229) * qJD(6)) * t198 + (t87 * t239 - t180 * t35 - t3 + (t180 * t85 - t21) * qJD(6)) * t194 + t266, t9 * t181 + t229 * t23 - t21 * t24 - t42 * t48 + (t195 * t42 + t336 * t328) * t317 + t211 * t180; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t353, t39, t37, -t353, t38, 0, t189 * t45 - t251 - t341, t250 * t44 + t204, 0, 0, t12, t6, t11, t308 * t85 - t34, -t194 * t350 ^ 2 + t335, -t300, -pkin(5) * t36 - pkin(10) * t319 + t102 * t312 - t27 * t350 - t45 * t85 + t252 - t321, pkin(5) * t35 + pkin(10) * t224 + t28 * t350 + t354 * t42 - t45 * t87 + t237, t27 * t87 + t28 * t85 + t1 + (t324 + (-t36 + t297) * pkin(10)) * t198 + ((qJD(6) * t85 - t35) * pkin(10) + t345) * t194, -pkin(5) * t9 + pkin(10) * t211 - t21 * t28 + t229 * t27 - t42 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t322, -t85 ^ 2 + t87 ^ 2, t350 * t85 - t35, -t322, t350 * t87 - t36, t56, -t42 * t87 - t345, t42 * t85 - t2 - t324, 0, 0;];
tauc_reg  = t7;
