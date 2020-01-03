% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRRP7
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRP7_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP7_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP7_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP7_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:57:50
% EndTime: 2019-12-31 21:58:01
% DurationCPUTime: 5.32s
% Computational Cost: add. (7494->530), mult. (17013->637), div. (0->0), fcn. (11844->10), ass. (0->267)
t195 = cos(qJ(4));
t298 = qJD(4) * t195;
t196 = cos(qJ(2));
t356 = cos(qJ(3));
t280 = t356 * t196;
t169 = qJD(1) * t280;
t192 = sin(qJ(3));
t193 = sin(qJ(2));
t301 = qJD(1) * t193;
t279 = t192 * t301;
t129 = -t169 + t279;
t323 = t129 * t195;
t380 = t298 + t323;
t190 = qJ(2) + qJ(3);
t186 = sin(t190);
t197 = cos(qJ(1));
t318 = t186 * t197;
t194 = sin(qJ(1));
t319 = t186 * t194;
t379 = g(1) * t318 + g(2) * t319;
t191 = sin(qJ(4));
t255 = t195 * pkin(4) + t191 * qJ(5);
t276 = t356 * qJD(2);
t264 = t196 * t276;
t293 = qJD(2) + qJD(3);
t268 = t192 * t293;
t275 = t356 * qJD(3);
t102 = t193 * t268 - t196 * t275 - t264;
t312 = t192 * t196;
t138 = t193 * t356 + t312;
t131 = t138 * qJD(1);
t267 = t195 * t293;
t109 = t131 * t191 - t267;
t111 = t131 * t195 + t191 * t293;
t326 = t111 * t195;
t329 = t109 * t191;
t272 = qJDD(1) * t356;
t297 = qJD(1) * qJD(2);
t274 = t193 * t297;
t305 = qJD(3) * t279 + t192 * t274;
t205 = t193 * t272 + (qJD(1) * (t276 + t275) + qJDD(1) * t192) * t196 - t305;
t292 = qJDD(2) + qJDD(3);
t52 = qJD(4) * t111 + t191 * t205 - t195 * t292;
t337 = t52 * t195;
t299 = qJD(4) * t191;
t51 = -qJD(4) * t267 + t131 * t299 - t191 * t292 - t195 * t205;
t338 = t51 * t191;
t378 = (qJD(4) * (-t326 + t329) - t337 + t338) * t138 + (t109 * t195 + t111 * t191) * t102;
t125 = qJD(4) + t129;
t327 = t111 * t125;
t331 = t109 * t125;
t5 = (t52 + t327) * t191 + (t51 + t331) * t195;
t359 = pkin(7) + pkin(6);
t103 = t293 * t138;
t295 = t193 * qJDD(1);
t253 = t192 * t295 - t196 * t272;
t76 = qJD(1) * t103 + t253;
t75 = qJDD(4) + t76;
t346 = qJ(5) * t75;
t174 = pkin(2) * t274;
t184 = t196 * pkin(2) + pkin(1);
t365 = -pkin(8) * t138 - t184;
t28 = t174 - (t169 * t293 - t305) * pkin(8) + t76 * pkin(3) + t365 * qJDD(1);
t273 = t196 * t297;
t106 = qJDD(2) * pkin(2) + t359 * (-t273 - t295);
t294 = t196 * qJDD(1);
t108 = t359 * (-t274 + t294);
t149 = t359 * t193;
t140 = qJD(1) * t149;
t345 = qJD(2) * pkin(2);
t135 = -t140 + t345;
t150 = t359 * t196;
t142 = qJD(1) * t150;
t300 = qJD(3) * t192;
t43 = t106 * t192 + t108 * t356 + t135 * t275 - t142 * t300;
t39 = pkin(8) * t292 + t43;
t148 = t184 * qJD(1);
t84 = t129 * pkin(3) - t131 * pkin(8) - t148;
t133 = t356 * t142;
t97 = t135 * t192 + t133;
t87 = pkin(8) * t293 + t97;
t8 = t191 * t28 + t195 * t39 + t298 * t84 - t299 * t87;
t2 = qJD(5) * t125 + t346 + t8;
t271 = t191 * t39 - t195 * t28 + t298 * t87 + t299 * t84;
t358 = pkin(4) * t75;
t4 = qJDD(5) + t271 - t358;
t377 = t4 * t191 + t2 * t195;
t339 = t195 * t75;
t376 = (t125 * t299 - t339) * pkin(8);
t355 = pkin(2) * t192;
t182 = pkin(8) + t355;
t265 = pkin(2) * t275;
t247 = t195 * t265;
t375 = -t109 * t247 - t182 * t337;
t266 = -t106 * t356 + t108 * t192 + t135 * t300 + t142 * t275;
t40 = -pkin(3) * t292 + t266;
t132 = t192 * t142;
t96 = t135 * t356 - t132;
t86 = -pkin(3) * t293 - t96;
t374 = t191 * t40 + t298 * t86;
t114 = -t149 * t192 + t150 * t356;
t137 = t192 * t193 - t280;
t95 = pkin(3) * t137 + t365;
t373 = t114 * t195 + t191 * t95;
t371 = t52 - t327;
t369 = t131 * t293;
t100 = -t140 * t192 + t133;
t261 = pkin(2) * t300 - t100;
t324 = t129 * t191;
t367 = -qJD(5) * t191 - t380 * qJ(5) + (t299 + t324) * pkin(4);
t366 = -t149 * t356 - t192 * t150;
t187 = cos(t190);
t304 = pkin(3) * t187 + pkin(8) * t186;
t258 = g(1) * t197 + g(2) * t194;
t333 = t102 * t191;
t234 = t138 * t298 - t333;
t321 = t138 * t191;
t364 = t103 * t109 + t125 * t234 + t137 * t52 + t321 * t75;
t290 = t193 * t345;
t66 = pkin(3) * t103 + pkin(8) * t102 + t290;
t283 = qJD(2) * t359;
t141 = t193 * t283;
t71 = qJD(3) * t366 - t141 * t356 - t283 * t312;
t16 = -qJD(4) * t373 - t191 * t71 + t195 * t66;
t361 = t111 ^ 2;
t360 = t125 ^ 2;
t357 = pkin(8) * t75;
t354 = pkin(4) * t131;
t178 = g(3) * t186;
t350 = g(3) * t187;
t349 = g(3) * t191;
t348 = g(3) * t196;
t7 = t8 * t195;
t92 = pkin(3) * t131 + pkin(8) * t129;
t60 = t191 * t92 + t195 * t96;
t101 = -t140 * t356 - t132;
t85 = pkin(2) * t301 + t92;
t58 = t101 * t195 + t191 * t85;
t10 = pkin(4) * t52 + qJ(5) * t51 - qJD(5) * t111 + t40;
t344 = t10 * t191;
t54 = t191 * t84 + t195 * t87;
t343 = t125 * t54;
t47 = pkin(4) * t109 - qJ(5) * t111 + t86;
t342 = t129 * t47;
t341 = t129 * t86;
t340 = t182 * t75;
t336 = t367 + t261;
t335 = -t97 + t367;
t334 = pkin(6) * qJDD(1);
t332 = t102 * t195;
t330 = t109 * t182;
t328 = t111 * t109;
t325 = t125 * t131;
t322 = t131 * t129;
t320 = t138 * t195;
t317 = t187 * t194;
t316 = t187 * t197;
t314 = t191 * t194;
t311 = t194 * t195;
t310 = t195 * t197;
t309 = t197 * t191;
t308 = t197 * t359;
t53 = -t191 * t87 + t195 * t84;
t307 = qJD(5) - t53;
t306 = t379 * t195;
t188 = t193 ^ 2;
t189 = t196 ^ 2;
t303 = t188 - t189;
t302 = t188 + t189;
t291 = t356 * pkin(2);
t200 = qJD(1) ^ 2;
t288 = t193 * t200 * t196;
t41 = t47 * t299;
t287 = t47 * t298;
t80 = t86 * t299;
t152 = t197 * t184;
t286 = pkin(3) * t316 + pkin(8) * t318 + t152;
t285 = t187 * t349 - t191 * t379;
t284 = g(1) * t316 + g(2) * t317 + t178;
t282 = t191 * t356;
t281 = t195 * t356;
t278 = t182 * t298;
t277 = t109 ^ 2 - t361;
t270 = t125 * t191;
t269 = t125 * t195;
t263 = t193 * t273;
t262 = g(1) * t319 - g(2) * t318;
t120 = t187 * t314 + t310;
t122 = t187 * t309 - t311;
t260 = g(1) * t120 - g(2) * t122;
t121 = t187 * t311 - t309;
t123 = t187 * t310 + t314;
t259 = g(1) * t121 - g(2) * t123;
t257 = g(1) * t194 - g(2) * t197;
t254 = pkin(4) * t191 - qJ(5) * t195;
t252 = -t340 + t341;
t32 = -pkin(4) * t125 + t307;
t33 = qJ(5) * t125 + t54;
t251 = t191 * t33 - t195 * t32;
t249 = t191 * t54 + t195 * t53;
t59 = -t191 * t96 + t195 * t92;
t57 = -t101 * t191 + t195 * t85;
t246 = -t131 * t109 - t339;
t67 = -t114 * t191 + t195 * t95;
t242 = t326 + t329;
t241 = t32 * t131 + t306 + t41;
t240 = -t53 * t131 + t306 + t80;
t239 = -t323 * t53 - t324 * t54 - t284 + t7;
t238 = pkin(3) + t255;
t237 = -pkin(8) * qJD(4) * t125 - t350;
t236 = -0.2e1 * pkin(1) * t297 - pkin(6) * qJDD(2);
t233 = -t138 * t299 - t332;
t15 = -t114 * t299 + t191 * t66 + t195 * t71 + t298 * t95;
t229 = t54 * t131 + t285 + t374;
t228 = t109 * t270 - t337;
t226 = t111 * t265 - t182 * t51;
t224 = -t33 * t131 - t323 * t47 - t285 - t344;
t223 = t193 * t258 - t348;
t222 = g(1) * t122 + g(2) * t120 + t186 * t349 - t271;
t221 = -t191 * t265 - t278;
t220 = t182 * t299 - t247;
t199 = qJD(2) ^ 2;
t219 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t199 + t257;
t218 = pkin(1) * t200 + t258 - t334;
t217 = -t148 * t129 + t284 - t43;
t216 = t32 * t380 - t33 * t324 - t284 + t377;
t215 = t148 * t131 - t266 - t350 + t379;
t214 = (-g(1) * (-t184 - t304) - g(2) * t359) * t194;
t213 = -qJD(4) * t251 + t377;
t212 = -qJD(4) * t249 + t191 * t271 + t7;
t211 = t109 * t234 + t321 * t52;
t72 = qJD(3) * t114 - t192 * t141 + t264 * t359;
t210 = t111 * t47 + qJDD(5) - t222;
t209 = -g(1) * t123 - g(2) * t121 - t178 * t195 + t8;
t155 = pkin(8) * t317;
t158 = pkin(8) * t316;
t207 = -g(1) * (-pkin(3) * t318 + t158) - g(2) * (-pkin(3) * t319 + t155) - g(3) * t304;
t204 = -t253 - t369;
t202 = -g(1) * t158 - g(2) * t155 - g(3) * (t187 * t255 + t304) + t258 * t186 * t238;
t183 = -t291 - pkin(3);
t136 = -t291 - t238;
t126 = -qJDD(1) * t184 + t174;
t119 = t131 * qJ(5);
t77 = -t129 ^ 2 + t131 ^ 2;
t73 = pkin(4) * t111 + qJ(5) * t109;
t70 = t138 * t254 - t366;
t65 = t204 + t369;
t64 = t129 * t293 + t205;
t56 = -pkin(4) * t137 - t67;
t55 = qJ(5) * t137 + t373;
t50 = -t59 - t354;
t49 = t119 + t60;
t48 = pkin(8) * t337;
t46 = -t57 - t354;
t45 = t119 + t58;
t30 = t103 * t125 + t137 * t75;
t24 = -t51 + t331;
t23 = -t111 * t131 + t125 * t269 + t191 * t75;
t22 = -t191 * t360 - t246;
t21 = t125 * t270 + t246;
t18 = t111 * t269 - t338;
t17 = -t254 * t102 + (qJD(4) * t255 - qJD(5) * t195) * t138 + t72;
t14 = t111 * t233 - t320 * t51;
t13 = -t103 * pkin(4) - t16;
t12 = qJ(5) * t103 + qJD(5) * t137 + t15;
t11 = t111 * t103 + t125 * t233 - t51 * t137 + t320 * t75;
t1 = [0, 0, 0, 0, 0, qJDD(1), t257, t258, 0, 0, qJDD(1) * t188 + 0.2e1 * t263, 0.2e1 * t193 * t294 - 0.2e1 * t297 * t303, qJDD(2) * t193 + t196 * t199, qJDD(1) * t189 - 0.2e1 * t263, qJDD(2) * t196 - t193 * t199, 0, t193 * t236 + t196 * t219, -t193 * t219 + t196 * t236, 0.2e1 * t302 * t334 - t258, -g(1) * (-pkin(1) * t194 + pkin(6) * t197) - g(2) * (pkin(1) * t197 + pkin(6) * t194) + (pkin(6) ^ 2 * t302 + pkin(1) ^ 2) * qJDD(1), -t131 * t102 + t138 * t205, t102 * t129 - t131 * t103 - t137 * t205 - t138 * t76, -t102 * t293 + t138 * t292, t103 * t129 + t137 * t76, -t103 * t293 - t137 * t292, 0, -t148 * t103 + t126 * t137 + t129 * t290 - t184 * t76 + t187 * t257 + t292 * t366 - t293 * t72, t148 * t102 - t114 * t292 + t126 * t138 + t131 * t290 - t184 * t205 - t293 * t71 - t262, t96 * t102 - t97 * t103 - t114 * t76 - t71 * t129 + t72 * t131 - t43 * t137 + t138 * t266 - t205 * t366 - t258, t43 * t114 + t97 * t71 - t266 * t366 - t96 * t72 - t126 * t184 - t148 * t290 - g(1) * (-t184 * t194 + t308) - g(2) * (t194 * t359 + t152), t14, t378, t11, t211, -t364, t30, t103 * t53 + t109 * t72 + t125 * t16 - t137 * t271 + t138 * t374 - t333 * t86 - t366 * t52 + t67 * t75 + t259, -t86 * t332 - t103 * t54 + t111 * t72 + t366 * t51 - t125 * t15 - t137 * t8 - t373 * t75 + (t40 * t195 - t80) * t138 - t260, -t109 * t15 - t111 * t16 + t51 * t67 - t52 * t373 + t249 * t102 + (-t191 * t8 + t195 * t271 + (t191 * t53 - t195 * t54) * qJD(4)) * t138 + t262, -g(1) * t308 - g(2) * t286 + t54 * t15 + t53 * t16 - t271 * t67 - t366 * t40 + t373 * t8 + t86 * t72 + t214, t14, t11, -t378, t30, t364, t211, -t47 * t333 - t103 * t32 + t109 * t17 - t125 * t13 - t137 * t4 + t52 * t70 - t56 * t75 + (t287 + t344) * t138 + t259, -t109 * t12 + t111 * t13 - t51 * t56 - t52 * t55 + t251 * t102 + (-t191 * t2 + t195 * t4 + (-t191 * t32 - t195 * t33) * qJD(4)) * t138 + t262, t47 * t332 + t103 * t33 - t111 * t17 + t12 * t125 + t137 * t2 + t51 * t70 + t55 * t75 + (-t10 * t195 + t41) * t138 + t260, t2 * t55 + t33 * t12 + t10 * t70 + t47 * t17 + t4 * t56 + t32 * t13 - g(1) * (-pkin(4) * t121 - qJ(5) * t120 + t308) - g(2) * (pkin(4) * t123 + qJ(5) * t122 + t286) + t214; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t288, t303 * t200, t295, t288, t294, qJDD(2), t193 * t218 - t348, g(3) * t193 + t196 * t218, 0, 0, t322, t77, t64, -t322, t65, t292, t100 * t293 + (-qJD(3) * t268 - t129 * t301 + t292 * t356) * pkin(2) + t215, t101 * t293 + (-t131 * t301 - t192 * t292 - t275 * t293) * pkin(2) + t217, -t205 * t291 - t355 * t76 + (t97 + t261) * t131 + (t101 - t265 - t96) * t129, t96 * t100 - t97 * t101 + (-t356 * t266 - t348 + t192 * t43 + (-t192 * t96 + t356 * t97) * qJD(3) + (qJD(1) * t148 + t258) * t193) * pkin(2), t18, -t5, t23, t228, t22, -t325, t183 * t52 + (-t40 - t350) * t195 + t252 * t191 + t261 * t109 + (-t57 + t221) * t125 + t240, -t183 * t51 + t252 * t195 + t261 * t111 + (t58 + t220) * t125 + t229, t58 * t109 + t57 * t111 + (t111 * t182 - t53) * t298 + (t271 + (-t54 + t330) * qJD(4) + t226) * t191 + t239 + t375, t40 * t183 - t54 * t58 - t53 * t57 - t86 * t100 + ((t192 * t86 + t281 * t54 - t282 * t53) * qJD(3) + t223) * pkin(2) + t212 * t182 + t207, t18, t23, t5, -t325, t21, t228, t136 * t52 + (-t10 - t350) * t195 + (-t340 + t342) * t191 + t336 * t109 + (t46 + t221) * t125 + t241, t45 * t109 + (-t46 + t278) * t111 + ((-t33 + t330) * qJD(4) + t226) * t191 + t216 + t375, t136 * t51 + (-qJD(4) * t47 + t340) * t195 - t336 * t111 + (-t45 - t220) * t125 + t224, t10 * t136 - t33 * t45 - t32 * t46 + t336 * t47 + ((t281 * t33 + t282 * t32) * qJD(3) + t223) * pkin(2) + t213 * t182 + t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t322, t77, t64, -t322, t65, t292, t293 * t97 + t215, t293 * t96 + t217, 0, 0, t18, -t5, t23, t228, t22, -t325, -pkin(3) * t52 - t109 * t97 - t125 * t59 + (t341 - t357) * t191 + (t237 - t40) * t195 + t240, pkin(3) * t51 - t111 * t97 + t125 * t60 + t323 * t86 + t229 + t376, t109 * t60 + t111 * t59 - t48 + (-pkin(8) * t51 + t271) * t191 + (pkin(8) * t242 - t249) * qJD(4) + t239, -t40 * pkin(3) + pkin(8) * t212 - t53 * t59 - t54 * t60 - t86 * t97 + t207, t18, t23, t5, -t325, t21, t228, t125 * t50 - t238 * t52 + (t342 - t357) * t191 + t335 * t109 + (-t10 + t237) * t195 + t241, -t33 * t299 + t109 * t49 - t111 * t50 - t48 + (qJD(4) * t242 - t338) * pkin(8) + t216, -t111 * t335 - t125 * t49 - t238 * t51 + t224 - t287 - t376, pkin(8) * t213 - t10 * t238 - t32 * t50 - t33 * t49 + t335 * t47 + t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t328, -t277, t24, -t328, -t371, t75, -t111 * t86 + t222 + t343, t109 * t86 + t125 * t53 - t209, 0, 0, t328, t24, t277, t75, t371, -t328, -t109 * t73 - t210 + t343 + 0.2e1 * t358, pkin(4) * t51 - qJ(5) * t52 + (t33 - t54) * t111 + (t32 - t307) * t109, 0.2e1 * t346 - t109 * t47 + t111 * t73 + (0.2e1 * qJD(5) - t53) * t125 + t209, t2 * qJ(5) - t4 * pkin(4) - t47 * t73 - t32 * t54 - g(1) * (-pkin(4) * t122 + qJ(5) * t123) - g(2) * (-pkin(4) * t120 + qJ(5) * t121) + t307 * t33 + t254 * t178; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) + t204 + t328, t24, -t360 - t361, -t125 * t33 + t210 - t358;];
tau_reg = t1;
