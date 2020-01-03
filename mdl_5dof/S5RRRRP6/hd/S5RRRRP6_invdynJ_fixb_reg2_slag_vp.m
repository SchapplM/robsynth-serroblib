% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRRP6
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
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRP6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:54:40
% EndTime: 2019-12-31 21:54:50
% DurationCPUTime: 5.01s
% Computational Cost: add. (7409->514), mult. (16799->631), div. (0->0), fcn. (11725->10), ass. (0->266)
t195 = qJ(2) + qJ(3);
t190 = sin(t195);
t203 = cos(qJ(1));
t310 = t190 * t203;
t200 = sin(qJ(1));
t311 = t190 * t200;
t367 = g(1) * t310 + g(2) * t311;
t198 = sin(qJ(3));
t199 = sin(qJ(2));
t202 = cos(qJ(2));
t355 = cos(qJ(3));
t263 = qJDD(1) * t355;
t269 = t355 * qJD(3);
t270 = t355 * qJD(2);
t292 = qJD(1) * qJD(2);
t268 = t199 * t292;
t296 = qJD(1) * t199;
t274 = t198 * t296;
t299 = qJD(3) * t274 + t198 * t268;
t209 = t199 * t263 + (qJD(1) * (t270 + t269) + qJDD(1) * t198) * t202 - t299;
t357 = pkin(7) + pkin(6);
t267 = t202 * t292;
t290 = t199 * qJDD(1);
t110 = qJDD(2) * pkin(2) + t357 * (-t267 - t290);
t289 = t202 * qJDD(1);
t113 = t357 * (-t268 + t289);
t157 = t357 * t199;
t146 = qJD(1) * t157;
t341 = qJD(2) * pkin(2);
t140 = -t146 + t341;
t158 = t357 * t202;
t148 = qJD(1) * t158;
t295 = qJD(3) * t198;
t255 = -t355 * t110 + t198 * t113 + t140 * t295 + t148 * t269;
t287 = qJDD(2) + qJDD(3);
t42 = -pkin(3) * t287 + t255;
t197 = sin(qJ(4));
t201 = cos(qJ(4));
t304 = t198 * t202;
t144 = t355 * t199 + t304;
t136 = t144 * qJD(1);
t288 = qJD(2) + qJD(3);
t116 = t201 * t136 + t197 * t288;
t363 = qJD(4) * t116;
t49 = t209 * t197 - t201 * t287 + t363;
t18 = t49 * pkin(4) + qJDD(5) + t42;
t293 = qJD(4) * t201;
t257 = t201 * t288;
t114 = t136 * t197 - t257;
t262 = t114 * pkin(4) + qJD(5);
t137 = t198 * t148;
t98 = t355 * t140 - t137;
t87 = -pkin(3) * t288 - t98;
t65 = t262 + t87;
t374 = t18 * t197 + t65 * t293;
t373 = t42 * t197 + t87 * t293;
t119 = -t198 * t157 + t355 * t158;
t111 = t201 * t119;
t275 = t355 * t202;
t143 = t198 * t199 - t275;
t343 = t202 * pkin(2);
t188 = pkin(1) + t343;
t366 = -pkin(8) * t144 - t188;
t97 = pkin(3) * t143 + t366;
t67 = t197 * t97 + t111;
t174 = qJD(1) * t275;
t134 = -t174 + t274;
t130 = qJD(4) + t134;
t323 = t114 * t130;
t294 = qJD(4) * t197;
t48 = -qJD(4) * t257 + t136 * t294 - t197 * t287 - t201 * t209;
t372 = -t48 - t323;
t320 = t116 * t130;
t371 = t49 + t320;
t256 = pkin(2) * t269;
t103 = -t355 * t146 - t137;
t95 = pkin(3) * t136 + pkin(8) * t134;
t86 = pkin(2) * t296 + t95;
t55 = t201 * t103 + t197 * t86;
t370 = -t201 * t256 + t55;
t138 = t355 * t148;
t102 = -t146 * t198 + t138;
t250 = pkin(2) * t295 - t102;
t369 = -t355 * t157 - t198 * t158;
t344 = t201 * pkin(4);
t186 = pkin(3) + t344;
t191 = cos(t195);
t196 = -qJ(5) - pkin(8);
t368 = t191 * t186 - t190 * t196;
t249 = t191 * pkin(3) + t190 * pkin(8);
t317 = t134 * t197;
t365 = -qJ(5) * t317 + t201 * qJD(5);
t156 = t188 * qJD(1);
t84 = t134 * pkin(3) - t136 * pkin(8) - t156;
t99 = t198 * t140 + t138;
t88 = pkin(8) * t288 + t99;
t50 = -t197 * t88 + t201 * t84;
t51 = t197 * t84 + t201 * t88;
t364 = -t197 * t50 + t201 * t51;
t245 = g(1) * t203 + g(2) * t200;
t362 = t51 * t136 + t373;
t302 = t201 * t203;
t307 = t197 * t200;
t123 = t191 * t307 + t302;
t303 = t200 * t201;
t306 = t197 * t203;
t125 = -t191 * t306 + t303;
t346 = g(3) * t197;
t361 = -g(1) * t125 + g(2) * t123 + t190 * t346;
t316 = t134 * t201;
t32 = -qJ(5) * t114 + t51;
t360 = t136 * t32 + t65 * t316 + t374;
t358 = t116 ^ 2;
t107 = t288 * t144;
t242 = t198 * t290 - t202 * t263;
t74 = qJD(1) * t107 + t242;
t73 = qJDD(4) + t74;
t356 = t73 * pkin(4);
t354 = pkin(2) * t198;
t159 = t203 * t188;
t349 = g(2) * t159;
t183 = g(3) * t190;
t347 = g(3) * t191;
t345 = g(3) * t202;
t179 = pkin(2) * t268;
t30 = t179 - (t288 * t174 - t299) * pkin(8) + t74 * pkin(3) + t366 * qJDD(1);
t43 = t198 * t110 + t355 * t113 + t140 * t269 - t148 * t295;
t41 = pkin(8) * t287 + t43;
t7 = t197 * t30 + t201 * t41 + t84 * t293 - t88 * t294;
t6 = t7 * t201;
t31 = -qJ(5) * t116 + t50;
t25 = pkin(4) * t130 + t31;
t342 = -t31 + t25;
t57 = t197 * t95 + t201 * t98;
t340 = t134 * t87;
t336 = t48 * qJ(5);
t335 = t48 * t197;
t334 = t49 * qJ(5);
t333 = t49 * t201;
t185 = pkin(8) + t354;
t301 = -qJ(5) - t185;
t260 = qJD(4) * t301;
t331 = t197 * t260 + t365 - t370;
t192 = t201 * qJ(5);
t243 = t136 * pkin(4) + t134 * t192;
t54 = -t103 * t197 + t201 * t86;
t330 = (-t256 - qJD(5)) * t197 + t201 * t260 - t243 - t54;
t264 = qJD(4) * t196;
t329 = t197 * t264 + t365 - t57;
t56 = -t197 * t98 + t201 * t95;
t328 = -t197 * qJD(5) + t201 * t264 - t243 - t56;
t120 = pkin(4) * t317;
t282 = pkin(4) * t294;
t327 = t282 + t120 + t250;
t326 = pkin(6) * qJDD(1);
t254 = t202 * t270;
t258 = t198 * t288;
t106 = t199 * t258 - t202 * t269 - t254;
t325 = t106 * t197;
t324 = t106 * t201;
t322 = t114 * t197;
t321 = t116 * t114;
t319 = t116 * t201;
t318 = t130 * t136;
t315 = t136 * t134;
t314 = t144 * t197;
t313 = t144 * t201;
t309 = t191 * t200;
t308 = t191 * t203;
t300 = t367 * t201;
t193 = t199 ^ 2;
t194 = t202 ^ 2;
t298 = t193 - t194;
t297 = t193 + t194;
t284 = t199 * t341;
t64 = pkin(3) * t107 + pkin(8) * t106 + t284;
t276 = qJD(2) * t357;
t147 = t199 * t276;
t69 = t369 * qJD(3) - t355 * t147 - t276 * t304;
t286 = t197 * t64 + t201 * t69 + t97 * t293;
t285 = t355 * pkin(2);
t281 = pkin(8) * qJD(4) * t130;
t206 = qJD(1) ^ 2;
t278 = t199 * t206 * t202;
t58 = t65 * t294;
t79 = t87 * t294;
t277 = g(1) * t308 + g(2) * t309 + t183;
t273 = t144 * t293;
t272 = -t18 - t347;
t271 = -t42 - t347;
t266 = pkin(4) * t197 + t357;
t265 = -t197 * t69 + t201 * t64;
t66 = -t119 * t197 + t201 * t97;
t259 = t130 * t201;
t187 = -t285 - pkin(3);
t253 = t199 * t267;
t72 = -t120 + t99;
t252 = -t72 + t282;
t251 = g(1) * t311 - g(2) * t310;
t248 = -pkin(8) * t73 + t340;
t247 = -g(1) * t123 - g(2) * t125;
t124 = -t191 * t303 + t306;
t126 = t191 * t302 + t307;
t246 = -g(1) * t124 - g(2) * t126;
t244 = g(1) * t200 - g(2) * t203;
t241 = -t185 * t73 + t340;
t240 = t197 * t32 + t201 * t25;
t239 = t197 * t51 + t201 * t50;
t236 = t186 * t190 + t191 * t196;
t235 = qJ(5) * t106 - qJD(5) * t144;
t234 = -t50 * t136 + t300 + t79;
t233 = -t50 * t316 - t51 * t317 - t277 + t6;
t232 = t245 * t190;
t231 = -0.2e1 * pkin(1) * t292 - pkin(6) * qJDD(2);
t230 = t273 - t325;
t229 = -t144 * t294 - t324;
t166 = t191 * t346;
t223 = -t197 * t232 + t166;
t222 = g(1) * t126 - g(2) * t124 + t201 * t183 - t7;
t8 = -qJD(4) * t51 - t197 * t41 + t201 * t30;
t221 = -qJD(4) * t239 - t8 * t197;
t205 = qJD(2) ^ 2;
t220 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t205 + t244;
t219 = pkin(1) * t206 + t245 - t326;
t218 = -t156 * t134 + t277 - t43;
t217 = t156 * t136 - t255 - t347 + t367;
t216 = t221 + t6;
t70 = qJD(3) * t119 - t198 * t147 + t357 * t254;
t215 = -t136 * t25 + t201 * t272 + t65 * t317 + t300 + t58;
t1 = -t116 * qJD(5) + t336 + t356 + t8;
t3 = -qJD(5) * t114 - t334 + t7;
t213 = -qJD(4) * t240 - t1 * t197 + t3 * t201 - t25 * t316 - t32 * t317 - t277;
t212 = t8 + t361;
t211 = -g(1) * (-pkin(3) * t310 + pkin(8) * t308) - g(2) * (-pkin(3) * t311 + pkin(8) * t309) - g(3) * t249;
t155 = pkin(8) * t201 + t192;
t154 = t196 * t197;
t153 = t187 - t344;
t142 = t185 * t201 + t192;
t141 = t301 * t197;
t131 = -qJDD(1) * t188 + t179;
t112 = t114 ^ 2;
t85 = pkin(4) * t314 - t369;
t76 = -t134 ^ 2 + t136 ^ 2;
t62 = t134 * t288 + t209;
t53 = -t112 + t358;
t52 = -qJ(5) * t314 + t67;
t46 = pkin(4) * t143 - t144 * t192 + t66;
t36 = pkin(4) * t230 + t70;
t34 = t107 * t130 + t143 * t73;
t27 = t320 - t49;
t26 = -t48 + t323;
t22 = -t116 * t136 + t130 * t259 + t197 * t73;
t21 = -t130 ^ 2 * t197 + t114 * t136 + t201 * t73;
t20 = t130 * t322 - t333;
t19 = t116 * t259 - t335;
t16 = -t67 * qJD(4) + t265;
t15 = -t119 * t294 + t286;
t14 = t114 * t230 + t49 * t314;
t13 = t116 * t229 - t48 * t313;
t12 = -qJ(5) * t273 + (-qJD(4) * t119 + t235) * t197 + t286;
t11 = t107 * pkin(4) + t235 * t201 + (-t111 + (qJ(5) * t144 - t97) * t197) * qJD(4) + t265;
t10 = -t114 * t107 - t130 * t230 - t49 * t143 - t73 * t314;
t9 = t116 * t107 + t130 * t229 - t48 * t143 + t73 * t313;
t5 = -t371 * t197 + t372 * t201;
t4 = (t114 * t201 + t116 * t197) * t106 + (t335 - t333 + (-t319 + t322) * qJD(4)) * t144;
t2 = [0, 0, 0, 0, 0, qJDD(1), t244, t245, 0, 0, qJDD(1) * t193 + 0.2e1 * t253, 0.2e1 * t199 * t289 - 0.2e1 * t292 * t298, qJDD(2) * t199 + t202 * t205, qJDD(1) * t194 - 0.2e1 * t253, qJDD(2) * t202 - t199 * t205, 0, t199 * t231 + t202 * t220, -t199 * t220 + t202 * t231, 0.2e1 * t297 * t326 - t245, -g(1) * (-pkin(1) * t200 + pkin(6) * t203) - g(2) * (pkin(1) * t203 + pkin(6) * t200) + (pkin(6) ^ 2 * t297 + pkin(1) ^ 2) * qJDD(1), -t136 * t106 + t209 * t144, t106 * t134 - t136 * t107 - t209 * t143 - t144 * t74, -t106 * t288 + t144 * t287, t107 * t134 + t143 * t74, -t107 * t288 - t143 * t287, 0, -t156 * t107 + t131 * t143 + t134 * t284 - t188 * t74 + t191 * t244 + t287 * t369 - t288 * t70, t156 * t106 - t119 * t287 + t131 * t144 + t136 * t284 - t188 * t209 - t288 * t69 - t251, t98 * t106 - t99 * t107 - t119 * t74 - t69 * t134 + t70 * t136 - t43 * t143 + t144 * t255 - t209 * t369 - t245, t43 * t119 + t99 * t69 - t255 * t369 - t98 * t70 - t131 * t188 - t156 * t284 - g(1) * (-t200 * t188 + t203 * t357) - g(2) * (t200 * t357 + t159), t13, t4, t9, t14, t10, t34, t50 * t107 + t70 * t114 + t16 * t130 + t8 * t143 + t373 * t144 - t87 * t325 - t369 * t49 + t66 * t73 + t246, -t87 * t324 - t51 * t107 + t70 * t116 + t369 * t48 - t15 * t130 - t7 * t143 - t67 * t73 + (t42 * t201 - t79) * t144 + t247, -t15 * t114 - t16 * t116 + t66 * t48 - t67 * t49 + t239 * t106 + (-t364 * qJD(4) - t197 * t7 - t201 * t8) * t144 + t251, -t349 - t42 * t369 + t51 * t15 + t50 * t16 + t8 * t66 + t7 * t67 + t87 * t70 + (-g(1) * t357 - g(2) * t249) * t203 + (-g(1) * (-t188 - t249) - g(2) * t357) * t200, t13, t4, t9, t14, t10, t34, t1 * t143 + t25 * t107 + t11 * t130 + t36 * t114 + t374 * t144 - t65 * t325 + t46 * t73 + t85 * t49 + t246, -t65 * t324 - t32 * t107 + t36 * t116 - t12 * t130 - t3 * t143 - t85 * t48 - t52 * t73 + (t18 * t201 - t58) * t144 + t247, -t11 * t116 - t12 * t114 + t46 * t48 - t52 * t49 + t240 * t106 + (-t1 * t201 - t197 * t3 + (t197 * t25 - t201 * t32) * qJD(4)) * t144 + t251, -t349 + t1 * t46 + t25 * t11 + t32 * t12 + t18 * t85 + t3 * t52 + t65 * t36 + (-g(1) * t266 - g(2) * t368) * t203 + (-g(1) * (-t188 - t368) - g(2) * t266) * t200; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t278, t298 * t206, t290, t278, t289, qJDD(2), t199 * t219 - t345, g(3) * t199 + t202 * t219, 0, 0, t315, t76, t62, -t315, -t242, t287, t102 * t288 + (-qJD(3) * t258 - t134 * t296 + t355 * t287) * pkin(2) + t217, t103 * t288 + (-t136 * t296 - t198 * t287 - t269 * t288) * pkin(2) + t218, -t209 * t285 - t74 * t354 + (t99 + t250) * t136 + (t103 - t256 - t98) * t134, t98 * t102 - t99 * t103 + (-t355 * t255 - t345 + t198 * t43 + (-t198 * t98 + t355 * t99) * qJD(3) + (qJD(1) * t156 + t245) * t199) * pkin(2), t19, t5, t22, t20, t21, -t318, t187 * t49 + t271 * t201 + t241 * t197 + t250 * t114 + (-t185 * t293 - t197 * t256 - t54) * t130 + t234, -t187 * t48 + t241 * t201 + t250 * t116 + (t185 * t294 + t370) * t130 + t223 + t362, t55 * t114 + t54 * t116 + (-t114 * t256 - t185 * t49 + (t116 * t185 - t50) * qJD(4)) * t201 + (t116 * t256 - t185 * t48 - t8 + (t114 * t185 - t51) * qJD(4)) * t197 + t233, t42 * t187 - t51 * t55 - t50 * t54 - t87 * t102 + (-t345 + t245 * t199 + (t198 * t87 + t364 * t355) * qJD(3)) * pkin(2) + t216 * t185 + t211, t19, t5, t22, t20, t21, -t318, t114 * t327 + t130 * t330 + t141 * t73 + t153 * t49 + t215, t116 * t327 - t130 * t331 - t142 * t73 - t153 * t48 + t223 + t360, -t114 * t331 - t116 * t330 + t141 * t48 - t142 * t49 + t213, t3 * t142 + t1 * t141 + t18 * t153 - g(3) * (t368 + t343) + t327 * t65 + t331 * t32 + t330 * t25 + t245 * (pkin(2) * t199 + t236); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t315, t76, t62, -t315, -t242, t287, t288 * t99 + t217, t288 * t98 + t218, 0, 0, t19, t5, t22, t20, t21, -t318, -pkin(3) * t49 - t99 * t114 - t56 * t130 + t248 * t197 + (t271 - t281) * t201 + t234, pkin(3) * t48 - t116 * t99 + t130 * t57 + t166 + t248 * t201 + (-t232 + t281) * t197 + t362, t114 * t57 + t116 * t56 + (-t335 - t333 + (t319 + t322) * qJD(4)) * pkin(8) + t221 + t233, -t42 * pkin(3) + pkin(8) * t216 - t50 * t56 - t51 * t57 - t87 * t99 + t211, t19, t5, t22, t20, t21, -t318, t114 * t252 + t130 * t328 + t154 * t73 - t186 * t49 + t215, -t116 * t72 - t155 * t73 + t186 * t48 + t166 - t329 * t130 + (pkin(4) * t363 - t232) * t197 + t360, -t114 * t329 - t116 * t328 + t154 * t48 - t155 * t49 + t213, -g(3) * t368 + t1 * t154 + t3 * t155 - t18 * t186 + t245 * t236 + t25 * t328 + t252 * t65 + t32 * t329; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t321, t53, t26, -t321, t27, t73, -t87 * t116 + t51 * t130 + t212, t114 * t87 + t130 * t50 + t222, 0, 0, t321, t53, t26, -t321, t27, t73, 0.2e1 * t356 + t336 + t130 * t32 + (-t262 - t65) * t116 + t212, -t358 * pkin(4) + t334 + t31 * t130 + (qJD(5) + t65) * t114 + t222, t48 * pkin(4) - t114 * t342, t342 * t32 + (-t65 * t116 + t1 + t361) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t371, t372, -t112 - t358, t32 * t114 + t25 * t116 - t272 - t367;];
tau_reg = t2;
