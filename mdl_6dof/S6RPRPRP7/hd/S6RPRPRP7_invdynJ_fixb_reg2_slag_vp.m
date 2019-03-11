% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPRPRP7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPRP7_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP7_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP7_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP7_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:22:53
% EndTime: 2019-03-09 03:23:02
% DurationCPUTime: 4.88s
% Computational Cost: add. (7125->526), mult. (14247->618), div. (0->0), fcn. (9720->10), ass. (0->276)
t194 = sin(qJ(5));
t197 = cos(qJ(5));
t200 = -pkin(1) - pkin(7);
t150 = t200 * qJD(1) + qJD(2);
t198 = cos(qJ(3));
t282 = qJD(1) * t198;
t114 = -qJ(4) * t282 + t198 * t150;
t110 = qJD(3) * pkin(3) + t114;
t191 = sin(pkin(9));
t195 = sin(qJ(3));
t246 = -qJ(4) * qJD(1) + t150;
t113 = t246 * t195;
t316 = cos(pkin(9));
t255 = t316 * t113;
t71 = t191 * t110 + t255;
t62 = qJD(3) * pkin(8) + t71;
t214 = -t191 * t198 - t195 * t316;
t121 = t214 * qJD(1);
t254 = t316 * t198;
t283 = qJD(1) * t195;
t261 = t191 * t283;
t124 = qJD(1) * t254 - t261;
t137 = pkin(3) * t283 + qJD(1) * qJ(2) + qJD(4);
t72 = -pkin(4) * t121 - pkin(8) * t124 + t137;
t38 = t194 * t72 + t197 * t62;
t276 = t197 * qJD(3);
t98 = t124 * t194 - t276;
t29 = -qJ(6) * t98 + t38;
t355 = qJD(5) - t121;
t361 = t355 * t29;
t37 = -t194 * t62 + t197 * t72;
t360 = t355 * t37;
t359 = t355 * t38;
t199 = cos(qJ(1));
t182 = g(2) * t199;
t196 = sin(qJ(1));
t352 = g(1) * t196 - t182;
t274 = qJD(1) * qJD(3);
t259 = t198 * t274;
t271 = t195 * qJDD(1);
t358 = t259 + t271;
t202 = qJD(1) ^ 2;
t357 = -qJ(2) * t202 - t352;
t247 = t197 * t355;
t260 = t195 * t274;
t253 = qJD(3) * t316;
t233 = qJD(1) * t253;
t249 = qJDD(1) * t316;
t270 = t198 * qJDD(1);
t264 = t191 * t270 + t195 * t249 + t198 * t233;
t92 = t191 * t260 - t264;
t90 = -qJDD(5) + t92;
t322 = t194 * t90;
t356 = t247 * t355 - t322;
t186 = qJDD(1) * qJ(2);
t149 = t200 * qJDD(1) + qJDD(2);
t189 = t195 ^ 2;
t190 = t198 ^ 2;
t286 = t189 + t190;
t251 = t286 * t149;
t180 = t195 * pkin(3);
t193 = -qJ(4) - pkin(7);
t353 = t199 * t180 + t196 * t193;
t100 = qJD(3) * t194 + t124 * t197;
t351 = qJD(5) * t100;
t184 = qJ(3) + pkin(9);
t172 = sin(t184);
t296 = t197 * t199;
t299 = t194 * t196;
t116 = -t172 * t299 + t296;
t297 = t196 * t197;
t298 = t194 * t199;
t118 = t172 * t298 + t297;
t173 = cos(t184);
t335 = g(3) * t173;
t350 = -g(1) * t116 - g(2) * t118 + t194 * t335;
t237 = g(1) * t199 + g(2) * t196;
t187 = qJD(1) * qJD(2);
t266 = 0.2e1 * t187;
t349 = 0.2e1 * t186 + t266 - t237;
t209 = t358 * t191 + t195 * t233 - t198 * t249;
t279 = qJD(5) * t194;
t54 = -qJD(5) * t276 - t194 * qJDD(3) + t124 * t279 + t197 * t209;
t329 = qJ(6) * t54;
t342 = pkin(5) * t90;
t133 = t198 * t149;
t273 = qJD(1) * qJD(4);
t281 = qJD(3) * t195;
t67 = -t198 * t273 - t150 * t281 + qJDD(3) * pkin(3) + t133 + (t260 - t270) * qJ(4);
t280 = qJD(3) * t198;
t77 = t246 * t280 + (-qJ(4) * qJDD(1) + t149 - t273) * t195;
t40 = t191 * t67 + t316 * t77;
t36 = qJDD(3) * pkin(8) + t40;
t106 = t358 * pkin(3) + qJDD(4) + t186 + t187;
t43 = -t92 * pkin(4) + pkin(8) * t209 + t106;
t6 = -qJD(5) * t38 - t194 * t36 + t197 * t43;
t1 = -qJD(6) * t100 + t329 - t342 + t6;
t348 = t1 + t361;
t126 = t214 * qJD(3);
t132 = -t191 * t195 + t254;
t347 = t124 * t126 - t132 * t209;
t336 = g(3) * t172;
t346 = -t173 * t352 + t336;
t55 = -t197 * qJDD(3) - t194 * t209 + t351;
t345 = t100 ^ 2;
t344 = t124 ^ 2;
t341 = t98 * pkin(5);
t340 = pkin(3) * t198;
t339 = pkin(5) * t194;
t334 = g(3) * t195;
t333 = g(3) * t197;
t332 = t197 * pkin(5);
t28 = -qJ(6) * t100 + t37;
t20 = pkin(5) * t355 + t28;
t331 = -t28 + t20;
t278 = qJD(5) * t197;
t330 = -t194 * t55 - t98 * t278;
t39 = -t191 * t77 + t316 * t67;
t107 = t191 * t113;
t81 = t114 * t316 - t107;
t83 = pkin(3) * t282 + pkin(4) * t124 - pkin(8) * t121;
t45 = t194 * t83 + t197 * t81;
t164 = qJ(2) + t180;
t91 = -pkin(4) * t214 - pkin(8) * t132 + t164;
t295 = qJ(4) - t200;
t138 = t295 * t195;
t252 = t295 * t198;
t95 = -t138 * t316 - t191 * t252;
t93 = t197 * t95;
t53 = t194 * t91 + t93;
t328 = qJ(6) * t55;
t327 = t100 * t98;
t326 = t355 * t98;
t325 = t121 * t98;
t324 = t124 * t98;
t323 = t194 * t54;
t321 = t194 * t98;
t320 = t197 * t55;
t85 = t197 * t90;
t319 = t197 * t98;
t162 = pkin(3) * t191 + pkin(8);
t294 = qJ(6) + t162;
t250 = qJD(5) * t294;
t308 = t121 * t194;
t318 = qJ(6) * t308 + qJD(6) * t197 - t194 * t250 - t45;
t44 = -t194 * t81 + t197 * t83;
t317 = -pkin(5) * t124 - qJD(6) * t194 - t44 + (qJ(6) * t121 - t250) * t197;
t315 = pkin(1) * qJDD(1);
t313 = t100 * t355;
t312 = t100 * t124;
t311 = t100 * t194;
t310 = t100 * t197;
t309 = t355 * t124;
t307 = t124 * t121;
t305 = t126 * t194;
t304 = t126 * t197;
t303 = t132 * t194;
t302 = t132 * t197;
t301 = t173 * t196;
t300 = t173 * t199;
t293 = t126 * qJD(3) + t132 * qJDD(3);
t267 = g(2) * t300;
t292 = t172 * t333 + t197 * t267;
t291 = g(1) * t300 + g(2) * t301;
t290 = (t266 + t186) * qJ(2);
t289 = t199 * pkin(1) + t196 * qJ(2);
t287 = t189 - t190;
t201 = qJD(3) ^ 2;
t285 = -t201 - t202;
t284 = qJD(1) * t137;
t277 = t121 * qJD(3);
t154 = pkin(3) * t280 + qJD(2);
t272 = qJDD(3) * t195;
t111 = -qJD(3) * t252 - qJD(4) * t195;
t210 = -qJD(4) * t198 + t281 * t295;
t76 = t111 * t316 + t191 * t210;
t123 = t191 * t281 - t198 * t253;
t82 = -pkin(4) * t123 - pkin(8) * t126 + t154;
t269 = t194 * t82 + t197 * t76 + t91 * t278;
t268 = g(1) * t301;
t265 = t198 * t202 * t195;
t263 = t196 * t180 + t289;
t262 = t132 * t278;
t179 = t199 * qJ(2);
t258 = -pkin(1) * t196 + t179;
t257 = qJD(6) + t341;
t256 = -t194 * t76 + t197 * t82;
t52 = -t194 * t95 + t197 * t91;
t5 = t194 * t43 + t197 * t36 + t72 * t278 - t62 * t279;
t75 = t111 * t191 - t316 * t210;
t80 = t114 * t191 + t255;
t94 = -t138 * t191 + t316 * t252;
t248 = t194 * t355;
t245 = t286 * qJDD(1);
t244 = qJDD(2) - t315;
t243 = -qJD(5) * t214 + qJD(1);
t242 = t195 * t259;
t35 = -qJDD(3) * pkin(4) - t39;
t163 = -pkin(3) * t316 - pkin(4);
t241 = pkin(4) * t173 + pkin(8) * t172;
t240 = pkin(4) * t172 - pkin(8) * t173;
t239 = g(1) * t118 - g(2) * t116;
t117 = t172 * t297 + t298;
t119 = t172 * t296 - t299;
t238 = -g(1) * t119 - g(2) * t117;
t235 = -qJD(5) * t162 * t355 - t35;
t2 = -qJD(6) * t98 - t328 + t5;
t232 = -t20 * t355 + t2;
t231 = t194 * t29 + t197 * t20;
t230 = t194 * t20 - t197 * t29;
t229 = t194 * t38 + t197 * t37;
t228 = t194 * t37 - t197 * t38;
t70 = t110 * t316 - t107;
t227 = t310 + t321;
t226 = t121 * t123 + t214 * t92;
t167 = pkin(4) + t332;
t192 = -qJ(6) - pkin(8);
t225 = t167 * t173 - t172 * t192;
t224 = t167 * t172 + t173 * t192;
t223 = -qJ(6) * t126 - qJD(6) * t132;
t222 = t258 + t353;
t221 = -t85 + (-t279 + t308) * t355;
t220 = qJD(3) * t123 + qJDD(3) * t214;
t219 = -t267 - t336;
t218 = -t193 * t199 + t263;
t217 = t262 + t305;
t216 = -t132 * t279 + t304;
t61 = -qJD(3) * pkin(4) - t70;
t215 = t162 * t90 + t355 * t61;
t213 = t352 * t172 + t335;
t212 = 0.2e1 * qJ(2) * t274 + qJDD(3) * t200;
t15 = pkin(5) * t55 + qJDD(6) + t35;
t208 = g(1) * t117 - g(2) * t119 + t173 * t333 - t5;
t206 = -qJD(5) * t229 - t194 * t6 + t197 * t5;
t205 = t71 * t123 - t70 * t126 - t39 * t132 + t214 * t40 + t352;
t204 = -t200 * t201 + t349;
t203 = t6 + t350;
t176 = qJDD(3) * t198;
t159 = t196 * t340;
t139 = t163 - t332;
t135 = t194 * t268;
t129 = t294 * t197;
t128 = t294 * t194;
t120 = t121 ^ 2;
t97 = t98 ^ 2;
t68 = pkin(5) * t303 + t94;
t56 = pkin(5) * t308 + t80;
t50 = -t97 + t345;
t49 = t257 + t61;
t48 = pkin(5) * t217 + t75;
t47 = -t123 * t355 + t214 * t90;
t46 = -qJ(6) * t303 + t53;
t34 = -pkin(5) * t214 - qJ(6) * t302 + t52;
t32 = t313 - t55;
t31 = -t54 + t326;
t26 = -t312 - t356;
t25 = -t312 + t356;
t24 = t221 + t324;
t23 = t221 - t324;
t22 = t248 * t98 - t320;
t21 = t100 * t247 - t323;
t19 = -qJD(5) * t53 + t256;
t18 = -t279 * t95 + t269;
t17 = -t132 * t330 + t305 * t98;
t16 = t100 * t216 - t302 * t54;
t14 = -qJ(6) * t262 + (-qJD(5) * t95 + t223) * t194 + t269;
t13 = -pkin(5) * t123 + t223 * t197 + (-t93 + (qJ(6) * t132 - t91) * t194) * qJD(5) + t256;
t12 = t123 * t98 + t214 * t55 - t217 * t355 + t303 * t90;
t11 = -t100 * t123 + t214 * t54 + t216 * t355 - t302 * t90;
t10 = -t214 * t322 - t126 * t98 - t132 * t55 + (t123 * t194 - t197 * t243) * t355;
t9 = -t214 * t85 - t100 * t126 + t132 * t54 + (t123 * t197 + t194 * t243) * t355;
t8 = (-t54 + t325) * t197 - t355 * t311 + t330;
t7 = (t54 + t325) * t197 + t100 * t248 + t330;
t4 = -(t311 + t319) * t126 + (t323 - t320 + (-t310 + t321) * qJD(5)) * t132;
t3 = (-t311 + t319) * t123 + t227 * qJD(1) - (qJD(5) * t227 - t320 - t323) * t214;
t27 = [0, 0, 0, 0, 0, qJDD(1), t352, t237, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - t352 - 0.2e1 * t315, t349, -pkin(1) * t244 - g(1) * t258 - g(2) * t289 + t290, qJDD(1) * t190 - 0.2e1 * t242, -0.2e1 * t195 * t270 + 0.2e1 * t274 * t287, -t195 * t201 + t176, qJDD(1) * t189 + 0.2e1 * t242, -t198 * t201 - t272, 0, t195 * t204 + t198 * t212, -t195 * t212 + t198 * t204, -t200 * t245 - t251 + t352, -g(1) * (t196 * t200 + t179) - g(2) * (pkin(7) * t199 + t289) + t200 * t251 + t290, t347, t126 * t121 + t124 * t123 + t132 * t92 - t209 * t214, t293, t226, t220, 0, -qJD(3) * t75 - qJDD(3) * t94 - t106 * t214 - t121 * t154 - t123 * t137 - t164 * t92 - t172 * t237, -t76 * qJD(3) - t95 * qJDD(3) + t106 * t132 + t154 * t124 + t137 * t126 - t164 * t209 - t291, t76 * t121 + t75 * t124 - t209 * t94 + t95 * t92 + t205, -g(1) * t222 - g(2) * t218 + t106 * t164 + t137 * t154 - t39 * t94 + t40 * t95 - t70 * t75 + t71 * t76, t16, t4, t11, t17, t12, t47, t61 * t305 + t355 * t19 - t123 * t37 - t214 * t6 - t52 * t90 + t55 * t94 + t75 * t98 + (t194 * t35 + t278 * t61) * t132 + t238, t61 * t304 + t100 * t75 - t355 * t18 + t123 * t38 + t214 * t5 + t53 * t90 - t54 * t94 + (t197 * t35 - t279 * t61) * t132 + t239, -t100 * t19 - t18 * t98 + t52 * t54 - t53 * t55 - t229 * t126 + (qJD(5) * t228 - t194 * t5 - t197 * t6) * t132 + t291, t5 * t53 + t38 * t18 + t6 * t52 + t37 * t19 + t35 * t94 + t61 * t75 - g(1) * (t199 * t240 + t222) - g(2) * (t196 * t240 + t218) t16, t4, t11, t17, t12, t47, t49 * t305 - t1 * t214 + t355 * t13 - t123 * t20 - t34 * t90 + t48 * t98 + t55 * t68 + (t15 * t194 + t278 * t49) * t132 + t238, t49 * t304 + t100 * t48 - t355 * t14 + t123 * t29 + t214 * t2 + t46 * t90 - t54 * t68 + (t15 * t197 - t279 * t49) * t132 + t239, -t100 * t13 - t14 * t98 + t34 * t54 - t46 * t55 - t231 * t126 + (qJD(5) * t230 - t1 * t197 - t194 * t2) * t132 + t291, t2 * t46 + t29 * t14 + t1 * t34 + t20 * t13 + t15 * t68 + t49 * t48 - g(1) * (t179 + t353) - g(2) * t263 + (-g(1) * t224 - g(2) * (-t193 + t339)) * t199 + (-g(1) * (-pkin(1) - t339) - g(2) * t224) * t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t202, t357 + t244, 0, 0, 0, 0, 0, 0, t195 * t285 + t176, t198 * t285 - t272, -t245, t251 + t357, 0, 0, 0, 0, 0, 0, qJD(1) * t121 + t293, -qJD(1) * t124 + t220, -t226 - t347, -t205 - t284, 0, 0, 0, 0, 0, 0, t10, t9, t3, -qJD(1) * t229 + t123 * t228 - t126 * t61 - t132 * t35 - t206 * t214 - t352, 0, 0, 0, 0, 0, 0, t10, t9, t3, -t126 * t49 - t132 * t15 + t230 * t123 - t231 * qJD(1) - (-qJD(5) * t231 - t1 * t194 + t197 * t2) * t214 - t352; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t265, -t287 * t202, t270, -t265, -t271, qJDD(3), t198 * t357 + t133 + t334, g(3) * t198 + (-t149 - t357) * t195, 0, 0, -t307, -t120 + t344, -t209 - t277, t307 (t124 + t261) * qJD(3) - t264, qJDD(3), t80 * qJD(3) - t137 * t124 + (qJDD(3) * t316 + t121 * t282) * pkin(3) + t39 + t346, qJD(3) * t81 - t121 * t137 + (-qJDD(3) * t191 - t124 * t282) * pkin(3) + t213 - t40 (t71 - t80) * t124 - (t81 - t70) * t121 + (t191 * t92 + t209 * t316) * pkin(3), t70 * t80 - t71 * t81 + (t316 * t39 + t334 + t191 * t40 + (-t352 - t284) * t198) * pkin(3), t21, t8, t25, t22, t24, -t309, -t355 * t44 - t124 * t37 + t163 * t55 - t80 * t98 + (t235 - t268) * t197 + t215 * t194 + t292, -t100 * t80 + t355 * t45 + t124 * t38 - t163 * t54 + t135 + t215 * t197 + (t219 - t235) * t194, t100 * t44 + t45 * t98 + (t121 * t37 - t162 * t55 + t5 + (t100 * t162 - t37) * qJD(5)) * t197 + (t121 * t38 - t162 * t54 - t6 + (t162 * t98 - t38) * qJD(5)) * t194 - t213, t35 * t163 - t38 * t45 - t37 * t44 - t61 * t80 - g(1) * (t196 * t241 + t159) - g(3) * (-t180 - t240) - (-t241 - t340) * t182 + t206 * t162, t21, t8, t25, t22, t24, -t309, -t124 * t20 + t128 * t90 + t139 * t55 - t56 * t98 + (-t15 - t268) * t197 + t317 * t355 + (-t121 * t49 + (t49 + t341) * qJD(5)) * t194 + t292, -t100 * t56 + t124 * t29 + t129 * t90 - t139 * t54 + t135 + t49 * t247 - t318 * t355 + (pkin(5) * t351 + t15 + t219) * t194, -t317 * t100 - t128 * t54 - t129 * t55 - t194 * t348 + t232 * t197 - t318 * t98 - t213, t2 * t129 - t1 * t128 + t15 * t139 - g(1) * (t196 * t225 + t159) - g(3) * (-t180 - t224) + (pkin(5) * t279 - t56) * t49 + t318 * t29 + t317 * t20 - (-t225 - t340) * t182; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t124 - t261) * qJD(3) + t264, -t209 + t277, -t120 - t344, -t121 * t71 + t124 * t70 + t106 - t237, 0, 0, 0, 0, 0, 0, t23, t26, t7, -t124 * t61 + (t6 + t359) * t197 + (t5 - t360) * t194 - t237, 0, 0, 0, 0, 0, 0, t23, t26, t7, -t124 * t49 + t232 * t194 + t197 * t348 - t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t327, t50, t31, -t327, t32, -t90, -t100 * t61 + t203 + t359, t61 * t98 + t208 + t360, 0, 0, t327, t50, t31, -t327, t32, -t90, -0.2e1 * t342 + t329 + t361 + (-t257 - t49) * t100 + t203, -pkin(5) * t345 + t328 + t355 * t28 + (qJD(6) + t49) * t98 + t208, pkin(5) * t54 - t331 * t98, t331 * t29 + (-t49 * t100 + t1 + t350) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55 + t313, -t54 - t326, -t97 - t345, t100 * t20 + t29 * t98 + t15 - t346;];
tau_reg  = t27;
