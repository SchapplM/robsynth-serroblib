% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRPR9
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRPR9_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR9_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR9_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR9_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:24:42
% EndTime: 2019-12-31 21:24:56
% DurationCPUTime: 6.37s
% Computational Cost: add. (35460->447), mult. (74477->631), div. (0->0), fcn. (53071->10), ass. (0->282)
t250 = sin(qJ(3));
t254 = cos(qJ(3));
t251 = sin(qJ(2));
t292 = qJD(1) * t251;
t220 = -t254 * qJD(2) + t250 * t292;
t221 = qJD(2) * t250 + t254 * t292;
t246 = sin(pkin(9));
t247 = cos(pkin(9));
t199 = t247 * t220 + t221 * t246;
t197 = t199 ^ 2;
t255 = cos(qJ(2));
t234 = qJD(1) * t255 - qJD(3);
t328 = t234 ^ 2;
t158 = -t328 - t197;
t286 = qJD(1) * qJD(2);
t236 = t251 * t286;
t285 = t255 * qJDD(1);
t225 = -t236 + t285;
t219 = -qJDD(3) + t225;
t201 = -t220 * t246 + t221 * t247;
t304 = t201 * t199;
t260 = -t219 - t304;
t335 = t247 * t260;
t111 = t158 * t246 + t335;
t257 = qJD(1) ^ 2;
t252 = sin(qJ(1));
t256 = cos(qJ(1));
t280 = t252 * g(1) - t256 * g(2);
t213 = qJDD(1) * pkin(1) + t257 * pkin(6) + t280;
t266 = -t225 + t236;
t237 = t251 * qJDD(1);
t281 = t255 * t286;
t224 = t237 + t281;
t267 = t224 + t281;
t169 = t266 * pkin(2) - t267 * pkin(7) - t213;
t269 = g(1) * t256 + g(2) * t252;
t315 = qJDD(1) * pkin(6);
t214 = -pkin(1) * t257 - t269 + t315;
t270 = -pkin(2) * t255 - pkin(7) * t251;
t274 = t257 * t270 + t214;
t323 = t251 * g(3);
t327 = qJD(2) ^ 2;
t179 = -t327 * pkin(2) + qJDD(2) * pkin(7) + t274 * t255 - t323;
t139 = -t254 * t169 + t250 * t179;
t263 = -t250 * qJDD(2) - t254 * t224;
t194 = -qJD(3) * t220 - t263;
t210 = t220 * t234;
t174 = t194 - t210;
t302 = t221 * t220;
t259 = -t219 - t302;
t106 = t259 * pkin(3) - t174 * qJ(4) - t139;
t140 = t250 * t169 + t254 * t179;
t205 = -pkin(3) * t234 - qJ(4) * t221;
t275 = -t254 * qJDD(2) + t250 * t224;
t262 = -qJD(3) * t221 - t275;
t329 = t220 ^ 2;
t108 = -t329 * pkin(3) + t262 * qJ(4) + t234 * t205 + t140;
t268 = -0.2e1 * qJD(4) * t201 + t247 * t106 - t246 * t108;
t339 = pkin(3) * t111 + t268;
t249 = sin(qJ(5));
t253 = cos(qJ(5));
t163 = t253 * t199 + t201 * t249;
t165 = -t199 * t249 + t201 * t253;
t118 = t165 * t163;
t215 = -qJDD(5) + t219;
t331 = -t118 - t215;
t338 = t249 * t331;
t337 = t253 * t331;
t336 = t246 * t260;
t334 = t250 * t259;
t333 = t254 * t259;
t157 = t247 * t194 + t246 * t262;
t277 = t194 * t246 - t247 * t262;
t100 = -qJD(5) * t163 + t157 * t253 - t249 * t277;
t230 = -qJD(5) + t234;
t147 = t163 * t230;
t332 = t100 + t147;
t184 = t199 * t234;
t135 = -t184 + t157;
t330 = t184 + t157;
t278 = t249 * t157 + t253 * t277;
t80 = (qJD(5) + t230) * t165 + t278;
t170 = (qJD(3) + t234) * t221 + t275;
t161 = t163 ^ 2;
t162 = t165 ^ 2;
t198 = t201 ^ 2;
t218 = t221 ^ 2;
t229 = t230 ^ 2;
t291 = qJD(4) * t199;
t192 = -0.2e1 * t291;
t293 = t246 * t106 + t247 * t108;
t65 = t192 + t293;
t38 = t246 * t65 + t247 * t268;
t326 = pkin(3) * t38;
t303 = t201 * t234;
t131 = t277 + t303;
t96 = -t131 * t246 - t135 * t247;
t325 = pkin(3) * t96;
t322 = t255 * g(3);
t52 = t260 * pkin(4) - t135 * pkin(8) + t268;
t180 = -pkin(4) * t234 - pkin(8) * t201;
t53 = -pkin(4) * t197 - pkin(8) * t277 + t180 * t234 + t65;
t29 = t249 * t53 - t253 * t52;
t30 = t249 * t52 + t253 * t53;
t15 = t249 * t30 - t253 * t29;
t321 = t15 * t246;
t320 = t15 * t247;
t178 = -qJDD(2) * pkin(2) - t327 * pkin(7) + t274 * t251 + t322;
t121 = -t262 * pkin(3) - t329 * qJ(4) + t221 * t205 + qJDD(4) + t178;
t78 = pkin(4) * t277 - t197 * pkin(8) + t201 * t180 + t121;
t319 = t249 * t78;
t318 = t250 * t38;
t317 = t253 * t78;
t316 = t254 * t38;
t113 = -t118 + t215;
t314 = t113 * t249;
t313 = t113 * t253;
t312 = t121 * t246;
t311 = t121 * t247;
t152 = t219 - t304;
t310 = t152 * t246;
t309 = t152 * t247;
t308 = t178 * t250;
t307 = t178 * t254;
t186 = t219 - t302;
t306 = t186 * t250;
t305 = t186 * t254;
t301 = t230 * t249;
t300 = t230 * t253;
t299 = t234 * t246;
t298 = t234 * t247;
t297 = t234 * t250;
t296 = t234 * t254;
t233 = t255 * t257 * t251;
t295 = t251 * (qJDD(2) + t233);
t294 = t255 * (-t233 + qJDD(2));
t289 = qJD(3) - t234;
t284 = t255 * t118;
t283 = t255 * t304;
t282 = t255 * t302;
t39 = -t246 * t268 + t247 * t65;
t16 = t249 * t29 + t253 * t30;
t102 = t139 * t250 + t254 * t140;
t206 = t214 * t251 + t322;
t207 = t214 * t255 - t323;
t276 = t251 * t206 + t255 * t207;
t7 = t16 * t246 + t320;
t273 = pkin(3) * t7 + pkin(4) * t15;
t83 = t100 - t147;
t48 = -t249 * t80 - t253 * t83;
t50 = t249 * t83 - t253 * t80;
t25 = t246 * t50 + t247 * t48;
t272 = pkin(3) * t25 + pkin(4) * t48;
t177 = -t198 - t328;
t119 = t177 * t247 + t310;
t271 = pkin(3) * t119 - t293;
t265 = t139 * t254 - t140 * t250;
t264 = -pkin(1) + t270;
t116 = -t229 - t161;
t76 = t116 * t249 + t337;
t77 = t116 * t253 - t338;
t44 = t246 * t77 + t247 * t76;
t261 = pkin(3) * t44 + pkin(4) * t76 - t29;
t142 = -t162 - t229;
t89 = t142 * t253 + t314;
t90 = -t142 * t249 + t313;
t55 = t246 * t90 + t247 * t89;
t258 = pkin(3) * t55 + pkin(4) * t89 - t30;
t244 = t255 ^ 2;
t243 = t251 ^ 2;
t241 = t244 * t257;
t239 = t243 * t257;
t226 = -0.2e1 * t236 + t285;
t223 = t237 + 0.2e1 * t281;
t211 = t255 * t219;
t209 = -t218 + t328;
t208 = -t328 + t329;
t203 = t218 - t329;
t202 = -t218 - t328;
t196 = -t328 - t329;
t185 = t218 + t329;
t182 = -t198 + t328;
t181 = t197 - t328;
t175 = t289 * t220 + t263;
t173 = t194 + t210;
t171 = -t289 * t221 - t275;
t166 = t198 - t197;
t160 = -t202 * t250 + t305;
t159 = t202 * t254 + t306;
t151 = t196 * t254 - t334;
t150 = t196 * t250 + t333;
t146 = -t162 + t229;
t145 = t161 - t229;
t144 = (t199 * t247 - t201 * t246) * t234;
t143 = (t199 * t246 + t201 * t247) * t234;
t141 = -t197 - t198;
t137 = -t170 * t254 + t174 * t250;
t130 = t277 - t303;
t129 = t157 * t247 + t201 * t299;
t128 = t157 * t246 - t201 * t298;
t127 = -t199 * t298 + t246 * t277;
t126 = -t199 * t299 - t247 * t277;
t125 = t181 * t247 + t310;
t124 = -t182 * t246 + t335;
t123 = t181 * t246 - t309;
t122 = t182 * t247 + t336;
t120 = -t177 * t246 + t309;
t117 = t162 - t161;
t112 = t158 * t247 - t336;
t110 = (t163 * t253 - t165 * t249) * t230;
t109 = (t163 * t249 + t165 * t253) * t230;
t103 = -t161 - t162;
t99 = -qJD(5) * t165 - t278;
t98 = -t131 * t247 + t135 * t246;
t97 = -t130 * t247 - t246 * t330;
t95 = -t130 * t246 + t247 * t330;
t94 = t145 * t253 + t314;
t93 = -t146 * t249 + t337;
t92 = t145 * t249 - t313;
t91 = t146 * t253 + t338;
t88 = -qJ(4) * t119 + t311;
t87 = -t119 * t250 + t120 * t254;
t86 = t119 * t254 + t120 * t250;
t85 = -qJ(4) * t111 + t312;
t79 = (qJD(5) - t230) * t165 + t278;
t75 = t100 * t253 + t165 * t301;
t74 = t100 * t249 - t165 * t300;
t73 = -t163 * t300 - t249 * t99;
t72 = -t163 * t301 + t253 * t99;
t71 = -t111 * t250 + t112 * t254;
t70 = t111 * t254 + t112 * t250;
t69 = -t109 * t246 + t110 * t247;
t68 = t109 * t247 + t110 * t246;
t67 = -pkin(3) * t330 + qJ(4) * t120 + t312;
t66 = -pkin(3) * t130 + qJ(4) * t112 - t311;
t62 = -t250 * t96 + t254 * t98;
t61 = t250 * t98 + t254 * t96;
t60 = -t246 * t92 + t247 * t94;
t59 = -t246 * t91 + t247 * t93;
t58 = t246 * t94 + t247 * t92;
t57 = t246 * t93 + t247 * t91;
t56 = -t246 * t89 + t247 * t90;
t54 = -pkin(8) * t89 + t317;
t49 = -t249 * t332 - t253 * t79;
t47 = -t249 * t79 + t253 * t332;
t46 = -pkin(8) * t76 + t319;
t45 = -t246 * t76 + t247 * t77;
t43 = -t246 * t74 + t247 * t75;
t42 = -t246 * t72 + t247 * t73;
t41 = t246 * t75 + t247 * t74;
t40 = t246 * t73 + t247 * t72;
t37 = -pkin(4) * t332 + pkin(8) * t90 + t319;
t36 = -pkin(4) * t79 + pkin(8) * t77 - t317;
t35 = -pkin(3) * t121 + qJ(4) * t39;
t34 = -t250 * t55 + t254 * t56;
t33 = t250 * t56 + t254 * t55;
t32 = -qJ(4) * t96 - t38;
t31 = -pkin(3) * t141 + qJ(4) * t98 + t39;
t27 = -t246 * t48 + t247 * t50;
t26 = -t246 * t47 + t247 * t49;
t24 = t246 * t49 + t247 * t47;
t23 = -t250 * t44 + t254 * t45;
t22 = t250 * t45 + t254 * t44;
t21 = t254 * t39 - t318;
t20 = t250 * t39 + t316;
t19 = -qJ(4) * t55 - t246 * t37 + t247 * t54;
t18 = -qJ(4) * t44 - t246 * t36 + t247 * t46;
t17 = -pkin(3) * t332 + qJ(4) * t56 + t246 * t54 + t247 * t37;
t14 = -t25 * t250 + t254 * t27;
t13 = t25 * t254 + t250 * t27;
t12 = -pkin(3) * t79 + qJ(4) * t45 + t246 * t46 + t247 * t36;
t11 = -pkin(4) * t78 + pkin(8) * t16;
t10 = -pkin(8) * t48 - t15;
t9 = -pkin(4) * t103 + pkin(8) * t50 + t16;
t8 = t16 * t247 - t321;
t6 = -qJ(4) * t25 + t10 * t247 - t246 * t9;
t5 = -pkin(3) * t103 + qJ(4) * t27 + t10 * t246 + t247 * t9;
t4 = -t250 * t7 + t254 * t8;
t3 = t250 * t8 + t254 * t7;
t2 = -pkin(8) * t320 - qJ(4) * t7 - t11 * t246;
t1 = -pkin(3) * t78 - pkin(8) * t321 + qJ(4) * t8 + t11 * t247;
t28 = [0, 0, 0, 0, 0, qJDD(1), t280, t269, 0, 0, t267 * t251, t223 * t255 + t226 * t251, t295 + t255 * (-t239 + t327), -t266 * t255, t251 * (t241 - t327) + t294, 0, t255 * t213 + pkin(1) * t226 + pkin(6) * (t255 * (-t241 - t327) - t295), -t251 * t213 - pkin(1) * t223 + pkin(6) * (-t294 - t251 * (-t239 - t327)), pkin(1) * (t239 + t241) + (t243 + t244) * t315 + t276, pkin(1) * t213 + pkin(6) * t276, t251 * (t194 * t254 + t221 * t297) - t282, t251 * (t171 * t254 - t173 * t250) - t255 * t203, t251 * (-t209 * t250 + t333) - t255 * t174, t251 * (-t220 * t296 - t250 * t262) + t282, t251 * (t208 * t254 + t306) + t255 * t170, t211 + t251 * (t220 * t254 - t221 * t250) * t234, t251 * (-pkin(7) * t150 + t308) + t255 * (-pkin(2) * t150 + t139) - pkin(1) * t150 + pkin(6) * (t151 * t255 - t171 * t251), t251 * (-pkin(7) * t159 + t307) + t255 * (-pkin(2) * t159 + t140) - pkin(1) * t159 + pkin(6) * (t160 * t255 - t175 * t251), t251 * t265 + pkin(6) * (t137 * t255 - t185 * t251) + t264 * (-t170 * t250 - t174 * t254), pkin(6) * (t102 * t255 + t178 * t251) - t264 * t265, t251 * (-t128 * t250 + t129 * t254) - t283, t251 * (-t250 * t95 + t254 * t97) - t255 * t166, t251 * (-t122 * t250 + t124 * t254) - t255 * t135, t251 * (-t126 * t250 + t127 * t254) + t283, t251 * (-t123 * t250 + t125 * t254) + t255 * t131, t251 * (-t143 * t250 + t144 * t254) + t211, t251 * (-pkin(7) * t70 - t250 * t66 + t254 * t85) + t255 * (-pkin(2) * t70 - t339) - pkin(1) * t70 + pkin(6) * (t130 * t251 + t255 * t71), t251 * (-pkin(7) * t86 - t250 * t67 + t254 * t88) + t255 * (-pkin(2) * t86 + t192 - t271) - pkin(1) * t86 + pkin(6) * (t251 * t330 + t255 * t87), t251 * (-pkin(7) * t61 - t250 * t31 + t254 * t32) + t255 * (-pkin(2) * t61 - t325) - pkin(1) * t61 + pkin(6) * (t141 * t251 + t255 * t62), t251 * (-pkin(7) * t20 - qJ(4) * t316 - t250 * t35) + t255 * (-pkin(2) * t20 - t326) - pkin(1) * t20 + pkin(6) * (t121 * t251 + t21 * t255), t251 * (-t250 * t41 + t254 * t43) - t284, t251 * (-t24 * t250 + t254 * t26) - t255 * t117, t251 * (-t250 * t57 + t254 * t59) - t255 * t83, t251 * (-t250 * t40 + t254 * t42) + t284, t251 * (-t250 * t58 + t254 * t60) + t255 * t80, t251 * (-t250 * t68 + t254 * t69) + t255 * t215, t251 * (-pkin(7) * t22 - t12 * t250 + t18 * t254) + t255 * (-pkin(2) * t22 - t261) - pkin(1) * t22 + pkin(6) * (t23 * t255 + t251 * t79), t251 * (-pkin(7) * t33 - t17 * t250 + t19 * t254) + t255 * (-pkin(2) * t33 - t258) - pkin(1) * t33 + pkin(6) * (t251 * t332 + t255 * t34), t251 * (-pkin(7) * t13 - t250 * t5 + t254 * t6) + t255 * (-pkin(2) * t13 - t272) - pkin(1) * t13 + pkin(6) * (t103 * t251 + t14 * t255), t251 * (-pkin(7) * t3 - t1 * t250 + t2 * t254) + t255 * (-pkin(2) * t3 - t273) - pkin(1) * t3 + pkin(6) * (t251 * t78 + t255 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t233, t239 - t241, t237, t233, t285, qJDD(2), -t206, -t207, 0, 0, t194 * t250 - t221 * t296, t171 * t250 + t173 * t254, t209 * t254 + t334, -t220 * t297 + t254 * t262, t208 * t250 - t305, (t220 * t250 + t221 * t254) * t234, pkin(2) * t171 + pkin(7) * t151 - t307, pkin(2) * t175 + pkin(7) * t160 + t308, pkin(2) * t185 + pkin(7) * t137 + t102, -pkin(2) * t178 + pkin(7) * t102, t128 * t254 + t129 * t250, t250 * t97 + t254 * t95, t122 * t254 + t124 * t250, t126 * t254 + t127 * t250, t123 * t254 + t125 * t250, t143 * t254 + t144 * t250, -pkin(2) * t130 + pkin(7) * t71 + t250 * t85 + t254 * t66, -pkin(2) * t330 + pkin(7) * t87 + t250 * t88 + t254 * t67, -pkin(2) * t141 + pkin(7) * t62 + t250 * t32 + t254 * t31, -pkin(2) * t121 + pkin(7) * t21 - qJ(4) * t318 + t254 * t35, t250 * t43 + t254 * t41, t24 * t254 + t250 * t26, t250 * t59 + t254 * t57, t250 * t42 + t254 * t40, t250 * t60 + t254 * t58, t250 * t69 + t254 * t68, -pkin(2) * t79 + pkin(7) * t23 + t12 * t254 + t18 * t250, -pkin(2) * t332 + pkin(7) * t34 + t17 * t254 + t19 * t250, -pkin(2) * t103 + pkin(7) * t14 + t250 * t6 + t254 * t5, -pkin(2) * t78 + pkin(7) * t4 + t1 * t254 + t2 * t250; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t302, t203, t174, -t302, -t170, -t219, -t139, -t140, 0, 0, t304, t166, t135, -t304, -t131, -t219, t339, t271 + 0.2e1 * t291, t325, t326, t118, t117, t83, -t118, -t80, -t215, t261, t258, t272, t273; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, t330, t141, t121, 0, 0, 0, 0, 0, 0, t79, t332, t103, t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, t117, t83, -t118, -t80, -t215, -t29, -t30, 0, 0;];
tauJ_reg = t28;
