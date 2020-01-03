% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPPR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPPR6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:33:08
% EndTime: 2019-12-31 19:33:22
% DurationCPUTime: 4.82s
% Computational Cost: add. (22626->435), mult. (53927->630), div. (0->0), fcn. (38036->10), ass. (0->253)
t229 = sin(pkin(9));
t230 = sin(pkin(8));
t232 = cos(pkin(8));
t237 = cos(qJ(2));
t267 = qJD(1) * t237;
t234 = sin(qJ(2));
t268 = qJD(1) * t234;
t204 = t230 * t267 + t232 * t268;
t231 = cos(pkin(9));
t186 = -qJD(2) * t231 + t204 * t229;
t188 = qJD(2) * t229 + t204 * t231;
t154 = t188 * t186;
t220 = t234 * qJDD(1);
t261 = qJD(1) * qJD(2);
t254 = t237 * t261;
t209 = t220 + t254;
t221 = t237 * qJDD(1);
t255 = t234 * t261;
t210 = t221 - t255;
t178 = t209 * t230 - t210 * t232;
t302 = -t154 + t178;
t310 = t229 * t302;
t202 = t230 * t268 - t232 * t267;
t177 = t204 * t202;
t300 = qJDD(2) - t177;
t309 = t230 * t300;
t308 = t231 * t302;
t307 = t232 * t300;
t233 = sin(qJ(5));
t236 = cos(qJ(5));
t148 = t186 * t236 + t188 * t233;
t150 = -t186 * t233 + t188 * t236;
t117 = t150 * t148;
t175 = qJDD(5) + t178;
t304 = -t117 + t175;
t306 = t233 * t304;
t305 = t236 * t304;
t179 = t209 * t232 + t210 * t230;
t168 = qJDD(2) * t229 + t179 * t231;
t249 = -qJDD(2) * t231 + t179 * t229;
t109 = -qJD(5) * t148 + t168 * t236 - t233 * t249;
t198 = qJD(5) + t202;
t138 = t198 * t148;
t303 = -t138 + t109;
t166 = t202 * t186;
t129 = -t166 - t168;
t301 = -t166 + t168;
t266 = qJD(2) * t204;
t157 = t178 + t266;
t196 = qJD(2) * t202;
t159 = t179 - t196;
t239 = qJD(1) ^ 2;
t235 = sin(qJ(1));
t296 = cos(qJ(1));
t246 = g(1) * t296 + g(2) * t235;
t287 = qJDD(1) * pkin(6);
t242 = -pkin(1) * t239 - t246 + t287;
t190 = -t234 * g(3) + t237 * t242;
t227 = t237 ^ 2;
t223 = t227 * t239;
t244 = qJD(2) * pkin(2) - qJ(3) * t268;
t152 = -pkin(2) * t223 + t210 * qJ(3) - qJD(2) * t244 + t190;
t241 = t234 * t242;
t272 = t234 * t239;
t240 = -t241 - t209 * qJ(3) + qJDD(2) * pkin(2) + (pkin(2) * t272 + qJ(3) * t261 - g(3)) * t237;
t112 = -0.2e1 * qJD(3) * t202 + t152 * t232 + t230 * t240;
t250 = t233 * t168 + t236 * t249;
t84 = (qJD(5) - t198) * t150 + t250;
t252 = t235 * g(1) - g(2) * t296;
t245 = qJDD(1) * pkin(1) + t252;
t155 = t210 * pkin(2) + (qJ(3) * t227 + pkin(6)) * t239 - t244 * t268 - qJDD(3) + t245;
t146 = t148 ^ 2;
t147 = t150 ^ 2;
t299 = t186 ^ 2;
t185 = t188 ^ 2;
t197 = t198 ^ 2;
t298 = t202 ^ 2;
t201 = t204 ^ 2;
t297 = 0.2e1 * qJD(3);
t295 = pkin(3) * t230;
t102 = pkin(3) * t157 - qJ(4) * t159 - t155;
t169 = pkin(3) * t202 - qJ(4) * t204;
t238 = qJD(2) ^ 2;
t91 = -pkin(3) * t238 + qJDD(2) * qJ(4) - t169 * t202 + t112;
t57 = 0.2e1 * qJD(4) * t188 - t102 * t231 + t229 * t91;
t42 = pkin(4) * t302 + pkin(7) * t129 - t57;
t161 = pkin(4) * t202 - pkin(7) * t188;
t58 = -0.2e1 * qJD(4) * t186 + t102 * t229 + t231 * t91;
t45 = -pkin(4) * t299 - pkin(7) * t249 - t161 * t202 + t58;
t21 = t233 * t45 - t236 * t42;
t22 = t233 * t42 + t236 * t45;
t11 = -t21 * t236 + t22 * t233;
t294 = t229 * t11;
t251 = t230 * t152 - t232 * t240;
t90 = qJDD(4) - t238 * qJ(4) - qJDD(2) * pkin(3) + (t297 + t169) * t204 + t251;
t293 = t229 * t90;
t292 = t231 * t11;
t291 = t231 * t90;
t61 = pkin(4) * t249 - pkin(7) * t299 + t188 * t161 + t90;
t290 = t233 * t61;
t111 = t204 * t297 + t251;
t65 = -t111 * t232 + t112 * t230;
t289 = t234 * t65;
t288 = t236 * t61;
t286 = t198 * t233;
t285 = t198 * t236;
t284 = t202 * t188;
t283 = t202 * t229;
t282 = t202 * t231;
t131 = t154 + t178;
t281 = t229 * t131;
t280 = t230 * t155;
t172 = qJDD(2) + t177;
t279 = t230 * t172;
t278 = t230 * t178;
t277 = t231 * t131;
t276 = t232 * t155;
t275 = t232 * t172;
t106 = t117 + t175;
t274 = t233 * t106;
t217 = t237 * t272;
t273 = t234 * (qJDD(2) + t217);
t271 = t236 * t106;
t270 = t237 * (qJDD(2) - t217);
t265 = qJD(2) * t230;
t264 = qJD(2) * t232;
t260 = t230 * t117;
t259 = t230 * t154;
t258 = t232 * t117;
t257 = t232 * t154;
t256 = -pkin(3) * t232 - pkin(2);
t12 = t21 * t233 + t22 * t236;
t34 = t229 * t57 + t231 * t58;
t66 = t111 * t230 + t112 * t232;
t189 = t237 * g(3) + t241;
t248 = t234 * t189 + t190 * t237;
t33 = t229 * t58 - t231 * t57;
t125 = t249 - t284;
t158 = -t178 + t266;
t226 = t234 ^ 2;
t222 = t226 * t239;
t211 = t221 - 0.2e1 * t255;
t208 = t220 + 0.2e1 * t254;
t206 = t239 * pkin(6) + t245;
t193 = -t201 - t238;
t192 = -t201 + t238;
t191 = t298 - t238;
t174 = t232 * t178;
t170 = -t298 - t238;
t164 = -t185 + t298;
t163 = -t298 + t299;
t160 = t179 + t196;
t156 = -t298 - t201;
t153 = -t185 + t299;
t143 = -t185 - t298;
t142 = -t193 * t230 - t275;
t141 = t193 * t232 - t279;
t140 = -t298 - t299;
t137 = -t147 + t197;
t136 = t146 - t197;
t135 = -t185 - t299;
t134 = t170 * t232 - t309;
t133 = t170 * t230 + t307;
t124 = t249 + t284;
t123 = t168 * t231 - t188 * t283;
t122 = t186 * t282 + t229 * t249;
t121 = (-t186 * t231 + t188 * t229) * t202;
t120 = -t147 - t197;
t119 = t158 * t232 + t160 * t230;
t118 = t158 * t230 - t160 * t232;
t116 = t147 - t146;
t115 = t163 * t231 - t281;
t114 = -t164 * t229 + t308;
t113 = -t197 - t146;
t108 = -qJD(5) * t150 - t250;
t104 = -t143 * t229 - t277;
t103 = t143 * t231 - t281;
t99 = (-t148 * t236 + t150 * t233) * t198;
t98 = (-t148 * t233 - t150 * t236) * t198;
t97 = t140 * t231 - t310;
t96 = t140 * t229 + t308;
t95 = -t125 * t231 - t129 * t229;
t94 = -t124 * t231 - t229 * t301;
t93 = -t125 * t229 + t129 * t231;
t92 = -t146 - t147;
t87 = t138 + t109;
t83 = (qJD(5) + t198) * t150 + t250;
t82 = t109 * t236 - t150 * t286;
t81 = t109 * t233 + t150 * t285;
t80 = -t108 * t233 + t148 * t285;
t79 = t108 * t236 + t148 * t286;
t78 = t136 * t236 - t274;
t77 = -t137 * t233 + t305;
t76 = t136 * t233 + t271;
t75 = t137 * t236 + t306;
t74 = t104 * t232 + t230 * t301;
t73 = t104 * t230 - t232 * t301;
t72 = t124 * t230 + t232 * t97;
t71 = -t124 * t232 + t230 * t97;
t70 = -t120 * t233 - t271;
t69 = t120 * t236 - t274;
t68 = t135 * t230 + t232 * t95;
t67 = -t135 * t232 + t230 * t95;
t64 = t113 * t236 - t306;
t63 = t113 * t233 + t305;
t62 = -t229 * t98 + t231 * t99;
t60 = -qJ(4) * t103 + t291;
t59 = -qJ(4) * t96 + t293;
t55 = t233 * t87 - t236 * t84;
t54 = -t233 * t303 - t236 * t83;
t53 = -t233 * t84 - t236 * t87;
t52 = -t233 * t83 + t236 * t303;
t51 = -t229 * t81 + t231 * t82;
t50 = -t229 * t79 + t231 * t80;
t49 = -t229 * t76 + t231 * t78;
t48 = -t229 * t75 + t231 * t77;
t47 = -t229 * t69 + t231 * t70;
t46 = t229 * t70 + t231 * t69;
t44 = -pkin(3) * t103 + t58;
t43 = -pkin(3) * t96 + t57;
t40 = -t229 * t63 + t231 * t64;
t39 = t229 * t64 + t231 * t63;
t38 = -pkin(7) * t69 + t288;
t37 = -pkin(7) * t63 + t290;
t36 = t230 * t303 + t232 * t47;
t35 = t230 * t47 - t232 * t303;
t32 = -pkin(4) * t303 + pkin(7) * t70 + t290;
t31 = t230 * t83 + t232 * t40;
t30 = t230 * t40 - t232 * t83;
t29 = -pkin(4) * t83 + pkin(7) * t64 - t288;
t28 = -t229 * t53 + t231 * t55;
t27 = -t229 * t52 + t231 * t54;
t26 = t229 * t55 + t231 * t53;
t25 = -qJ(4) * t93 - t33;
t23 = t230 * t34 - t232 * t90;
t19 = t230 * t92 + t232 * t28;
t18 = t230 * t28 - t232 * t92;
t17 = -pkin(3) * t26 - pkin(4) * t53;
t16 = -pkin(3) * t46 - pkin(4) * t69 + t22;
t15 = -qJ(4) * t46 - t229 * t32 + t231 * t38;
t14 = -pkin(3) * t39 - pkin(4) * t63 + t21;
t13 = -qJ(4) * t39 - t229 * t29 + t231 * t37;
t10 = -pkin(4) * t61 + pkin(7) * t12;
t9 = -pkin(7) * t53 - t11;
t8 = -pkin(4) * t92 + pkin(7) * t55 + t12;
t7 = t12 * t231 - t294;
t6 = t12 * t229 + t292;
t5 = t230 * t61 + t232 * t7;
t4 = t230 * t7 - t232 * t61;
t3 = -qJ(4) * t26 - t229 * t8 + t231 * t9;
t2 = -pkin(3) * t6 - pkin(4) * t11;
t1 = -pkin(7) * t292 - qJ(4) * t6 - t10 * t229;
t20 = [0, 0, 0, 0, 0, qJDD(1), t252, t246, 0, 0, (t209 + t254) * t234, t208 * t237 + t211 * t234, t273 + t237 * (-t222 + t238), (t210 - t255) * t237, t234 * (t223 - t238) + t270, 0, t237 * t206 + pkin(1) * t211 + pkin(6) * (t237 * (-t223 - t238) - t273), -t234 * t206 - pkin(1) * t208 + pkin(6) * (-t270 - t234 * (-t222 - t238)), pkin(1) * (t222 + t223) + (t226 + t227) * t287 + t248, pkin(1) * t206 + pkin(6) * t248, t234 * (t179 * t232 - t204 * t265) + t237 * (t179 * t230 + t204 * t264), t234 * (-t157 * t232 - t159 * t230) + t237 * (-t157 * t230 + t159 * t232), t234 * (-t192 * t230 + t307) + t237 * (t192 * t232 + t309), t234 * (t202 * t264 + t278) + t237 * (t202 * t265 - t174), t234 * (t191 * t232 - t279) + t237 * (t191 * t230 + t275), (t234 * (-t202 * t232 + t204 * t230) + t237 * (-t202 * t230 - t204 * t232)) * qJD(2), t234 * (-qJ(3) * t133 - t280) + t237 * (-pkin(2) * t157 + qJ(3) * t134 + t276) - pkin(1) * t157 + pkin(6) * (-t133 * t234 + t134 * t237), t234 * (-qJ(3) * t141 - t276) + t237 * (-pkin(2) * t159 + qJ(3) * t142 - t280) - pkin(1) * t159 + pkin(6) * (-t141 * t234 + t142 * t237), t234 * (-qJ(3) * t118 - t65) + t237 * (-pkin(2) * t156 + qJ(3) * t119 + t66) - pkin(1) * t156 + pkin(6) * (-t118 * t234 + t119 * t237), -qJ(3) * t289 + t237 * (pkin(2) * t155 + qJ(3) * t66) + pkin(1) * t155 + pkin(6) * (t237 * t66 - t289), t234 * (t123 * t232 + t259) + t237 * (t123 * t230 - t257), t234 * (-t153 * t230 + t232 * t94) + t237 * (t153 * t232 + t230 * t94), t234 * (t114 * t232 - t129 * t230) + t237 * (t114 * t230 + t129 * t232), t234 * (t122 * t232 - t259) + t237 * (t122 * t230 + t257), t234 * (t115 * t232 - t125 * t230) + t237 * (t115 * t230 + t125 * t232), t234 * (t121 * t232 + t278) + t237 * (t121 * t230 - t174), t234 * (-qJ(3) * t71 - t230 * t43 + t232 * t59) + t237 * (-pkin(2) * t96 + qJ(3) * t72 + t230 * t59 + t232 * t43) - pkin(1) * t96 + pkin(6) * (-t234 * t71 + t237 * t72), t234 * (-qJ(3) * t73 - t230 * t44 + t232 * t60) + t237 * (-pkin(2) * t103 + qJ(3) * t74 + t230 * t60 + t232 * t44) - pkin(1) * t103 + pkin(6) * (-t234 * t73 + t237 * t74), t234 * (-qJ(3) * t67 + t232 * t25) + t237 * (qJ(3) * t68 + t230 * t25) + pkin(6) * (-t234 * t67 + t237 * t68) + (t234 * t295 + t237 * t256 - pkin(1)) * t93, (t234 * (-qJ(4) * t232 + t295) + t237 * (-qJ(4) * t230 + t256) - pkin(1)) * t33 + (pkin(6) + qJ(3)) * (-t234 * t23 + t237 * (t230 * t90 + t232 * t34)), t234 * (t232 * t51 + t260) + t237 * (t230 * t51 - t258), t234 * (t116 * t230 + t232 * t27) + t237 * (-t116 * t232 + t230 * t27), t234 * (t230 * t87 + t232 * t48) + t237 * (t230 * t48 - t232 * t87), t234 * (t232 * t50 - t260) + t237 * (t230 * t50 + t258), t234 * (-t230 * t84 + t232 * t49) + t237 * (t230 * t49 + t232 * t84), t234 * (t175 * t230 + t232 * t62) + t237 * (-t175 * t232 + t230 * t62), t234 * (-qJ(3) * t30 + t13 * t232 - t14 * t230) + t237 * (-pkin(2) * t39 + qJ(3) * t31 + t13 * t230 + t14 * t232) - pkin(1) * t39 + pkin(6) * (-t234 * t30 + t237 * t31), t234 * (-qJ(3) * t35 + t15 * t232 - t16 * t230) + t237 * (-pkin(2) * t46 + qJ(3) * t36 + t15 * t230 + t16 * t232) - pkin(1) * t46 + pkin(6) * (-t234 * t35 + t237 * t36), t234 * (-qJ(3) * t18 - t17 * t230 + t232 * t3) + t237 * (-pkin(2) * t26 + qJ(3) * t19 + t17 * t232 + t230 * t3) - pkin(1) * t26 + pkin(6) * (-t18 * t234 + t19 * t237), t234 * (-qJ(3) * t4 + t1 * t232 - t2 * t230) + t237 * (-pkin(2) * t6 + qJ(3) * t5 + t1 * t230 + t2 * t232) - pkin(1) * t6 + pkin(6) * (-t234 * t4 + t237 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t217, t222 - t223, t220, t217, t221, qJDD(2), -t189, -t190, 0, 0, t177, t201 - t298, t160, -t177, t158, qJDD(2), pkin(2) * t133 - t111, pkin(2) * t141 - t112, pkin(2) * t118, pkin(2) * t65, t168 * t229 + t188 * t282, -t124 * t229 + t231 * t301, t164 * t231 + t310, t186 * t283 - t231 * t249, t163 * t229 + t277, (-t186 * t229 - t188 * t231) * t202, pkin(2) * t71 - pkin(3) * t124 + qJ(4) * t97 - t291, pkin(2) * t73 - pkin(3) * t301 + qJ(4) * t104 + t293, pkin(2) * t67 - pkin(3) * t135 + qJ(4) * t95 + t34, pkin(2) * t23 - pkin(3) * t90 + qJ(4) * t34, t229 * t82 + t231 * t81, t229 * t54 + t231 * t52, t229 * t77 + t231 * t75, t229 * t80 + t231 * t79, t229 * t78 + t231 * t76, t229 * t99 + t231 * t98, pkin(2) * t30 - pkin(3) * t83 + qJ(4) * t40 + t229 * t37 + t231 * t29, pkin(2) * t35 - pkin(3) * t303 + qJ(4) * t47 + t229 * t38 + t231 * t32, pkin(2) * t18 - pkin(3) * t92 + qJ(4) * t28 + t229 * t9 + t231 * t8, pkin(2) * t4 - pkin(3) * t61 - pkin(7) * t294 + qJ(4) * t7 + t10 * t231; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t157, t159, t156, -t155, 0, 0, 0, 0, 0, 0, t96, t103, t93, t33, 0, 0, 0, 0, 0, 0, t39, t46, t26, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, t301, t135, t90, 0, 0, 0, 0, 0, 0, t83, t303, t92, t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117, t116, t87, -t117, -t84, t175, -t21, -t22, 0, 0;];
tauJ_reg = t20;
