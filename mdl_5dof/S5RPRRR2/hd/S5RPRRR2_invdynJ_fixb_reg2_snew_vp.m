% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRRR2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRRR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:12:19
% EndTime: 2019-12-05 18:12:35
% DurationCPUTime: 5.74s
% Computational Cost: add. (28329->416), mult. (71003->609), div. (0->0), fcn. (55277->10), ass. (0->263)
t237 = qJD(1) ^ 2;
t232 = sin(qJ(1));
t296 = cos(qJ(1));
t244 = t296 * g(1) + t232 * g(2);
t302 = -t237 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) - t244;
t229 = sin(qJ(5));
t227 = sin(pkin(9));
t228 = cos(pkin(9));
t231 = sin(qJ(3));
t235 = cos(qJ(3));
t207 = (t227 * t231 - t228 * t235) * qJD(1);
t245 = t227 * t235 + t228 * t231;
t209 = t245 * qJD(1);
t230 = sin(qJ(4));
t234 = cos(qJ(4));
t180 = t234 * t207 + t209 * t230;
t182 = -t207 * t230 + t209 * t234;
t233 = cos(qJ(5));
t154 = t233 * t180 + t182 * t229;
t156 = -t180 * t229 + t182 * t233;
t120 = t156 * t154;
t223 = qJDD(3) + qJDD(4);
t217 = qJDD(5) + t223;
t301 = -t120 + t217;
t308 = t229 * t301;
t162 = t182 * t180;
t300 = -t162 + t223;
t307 = t230 * t300;
t195 = t209 * t207;
t298 = qJDD(3) - t195;
t306 = t231 * t298;
t305 = t233 * t301;
t304 = t234 * t300;
t303 = t235 * t298;
t262 = t228 * qJDD(1);
t263 = t227 * qJDD(1);
t169 = -t231 * t263 + t235 * t262;
t267 = t209 * qJD(3);
t191 = t169 - t267;
t206 = t245 * qJDD(1);
t268 = qJD(3) * t207;
t193 = t206 - t268;
t246 = t191 * t230 + t193 * t234;
t141 = -qJD(4) * t180 + t246;
t226 = qJD(3) + qJD(4);
t175 = t226 * t180;
t130 = t141 + t175;
t294 = t228 * pkin(2);
t295 = g(3) * t228;
t299 = -t295 + (-pkin(6) * qJDD(1) + t237 * t294 - t302) * t227;
t258 = g(1) * t232 - t296 * g(2);
t249 = -qJDD(2) + t258;
t286 = qJDD(1) * pkin(1);
t287 = qJ(2) * t237;
t203 = t249 + t286 + t287;
t224 = t227 ^ 2;
t225 = t228 ^ 2;
t272 = t225 * t237;
t297 = qJ(2) * t272 + t224 * t287 - t203 - t286;
t152 = t154 ^ 2;
t153 = t156 ^ 2;
t178 = t180 ^ 2;
t179 = t182 ^ 2;
t204 = t207 ^ 2;
t205 = t209 ^ 2;
t220 = qJD(5) + t226;
t216 = t220 ^ 2;
t222 = t226 ^ 2;
t251 = -g(3) * t227 + t302 * t228;
t176 = -pkin(2) * t272 + pkin(6) * t262 + t251;
t144 = t176 * t231 - t235 * t299;
t111 = (-t193 - t268) * pkin(7) + t298 * pkin(3) - t144;
t145 = t235 * t176 + t299 * t231;
t248 = qJD(3) * pkin(3) - pkin(7) * t209;
t119 = -t204 * pkin(3) + t191 * pkin(7) - qJD(3) * t248 + t145;
t81 = -t234 * t111 + t119 * t230;
t57 = t300 * pkin(4) - t130 * pkin(8) - t81;
t254 = -t234 * t191 + t193 * t230;
t140 = -qJD(4) * t182 - t254;
t250 = pkin(4) * t226 - pkin(8) * t182;
t82 = t230 * t111 + t234 * t119;
t58 = -t178 * pkin(4) + t140 * pkin(8) - t226 * t250 + t82;
t30 = t229 * t58 - t233 * t57;
t31 = t229 * t57 + t233 * t58;
t18 = t229 * t31 - t233 * t30;
t293 = t18 * t230;
t292 = t18 * t234;
t269 = t224 + t225;
t186 = (pkin(1) + t294) * qJDD(1) + (t269 * pkin(6) + qJ(2)) * t237 + t249;
t137 = pkin(3) * t191 + pkin(7) * t204 - t209 * t248 + t186;
t79 = pkin(4) * t140 + pkin(8) * t178 - t182 * t250 + t137;
t291 = t229 * t79;
t49 = t230 * t82 - t234 * t81;
t290 = t231 * t49;
t289 = t233 * t79;
t288 = t235 * t49;
t106 = -t144 * t235 + t145 * t231;
t285 = t106 * t227;
t114 = t120 + t217;
t284 = t114 * t229;
t283 = t114 * t233;
t282 = t137 * t230;
t281 = t137 * t234;
t159 = t162 + t223;
t280 = t159 * t230;
t279 = t159 * t234;
t278 = t186 * t231;
t277 = t186 * t235;
t188 = qJDD(3) + t195;
t276 = t188 * t231;
t275 = t188 * t235;
t274 = t220 * t229;
t273 = t220 * t233;
t271 = t226 * t230;
t270 = t226 * t234;
t266 = qJD(4) + t226;
t265 = qJD(5) + t220;
t19 = t229 * t30 + t233 * t31;
t50 = t230 * t81 + t234 * t82;
t255 = -t233 * t140 + t141 * t229;
t107 = t144 * t231 + t235 * t145;
t253 = t227 * (t302 * t227 + t295) + t228 * t251;
t112 = -t216 - t152;
t84 = t112 * t229 + t305;
t252 = pkin(4) * t84 - t30;
t247 = t140 * t229 + t141 * t233;
t143 = -t153 - t216;
t98 = t143 * t233 - t284;
t243 = pkin(4) * t98 - t31;
t242 = (-qJD(5) + t220) * t156 - t255;
t241 = (-qJD(4) + t226) * t182 - t254;
t89 = -qJD(5) * t154 + t247;
t236 = qJD(3) ^ 2;
t219 = t225 * qJDD(1);
t218 = t224 * qJDD(1);
t211 = t269 * t237;
t199 = -t205 - t236;
t198 = -t205 + t236;
t197 = t204 - t236;
t192 = t206 - 0.2e1 * t268;
t190 = -t169 + 0.2e1 * t267;
t185 = -t236 - t204;
t174 = -t179 + t222;
t173 = t178 - t222;
t172 = -t179 - t222;
t170 = -t204 - t205;
t168 = -t199 * t231 - t275;
t167 = t199 * t235 - t276;
t166 = t169 * t235 + t206 * t231;
t165 = t169 * t231 - t206 * t235;
t164 = t185 * t235 - t306;
t163 = t185 * t231 + t303;
t161 = t179 - t178;
t157 = -t222 - t178;
t150 = t220 * t154;
t149 = -t153 + t216;
t148 = t152 - t216;
t147 = (-t180 * t234 + t182 * t230) * t226;
t146 = (-t180 * t230 - t182 * t234) * t226;
t139 = -t178 - t179;
t136 = t173 * t234 - t280;
t135 = -t174 * t230 + t304;
t134 = t173 * t230 + t279;
t133 = t174 * t234 + t307;
t132 = -t172 * t230 - t279;
t131 = t172 * t234 - t280;
t129 = t141 - t175;
t128 = -t266 * t180 + t246;
t125 = t266 * t182 + t254;
t124 = t141 * t234 - t182 * t271;
t123 = t141 * t230 + t182 * t270;
t122 = -t140 * t230 + t180 * t270;
t121 = t140 * t234 + t180 * t271;
t118 = t153 - t152;
t117 = t157 * t234 - t307;
t116 = t157 * t230 + t304;
t110 = (-t154 * t233 + t156 * t229) * t220;
t109 = (-t154 * t229 - t156 * t233) * t220;
t105 = -t152 - t153;
t104 = -pkin(7) * t131 - t281;
t103 = t148 * t233 - t284;
t102 = -t149 * t229 + t305;
t101 = t148 * t229 + t283;
t100 = t149 * t233 + t308;
t99 = -t143 * t229 - t283;
t96 = -t131 * t231 + t132 * t235;
t95 = t131 * t235 + t132 * t231;
t94 = -pkin(7) * t116 - t282;
t93 = t130 * t230 + t234 * t241;
t92 = -t125 * t234 - t129 * t230;
t91 = -t130 * t234 + t230 * t241;
t90 = -t125 * t230 + t129 * t234;
t88 = -qJD(5) * t156 - t255;
t87 = -t116 * t231 + t117 * t235;
t86 = t116 * t235 + t117 * t231;
t85 = t112 * t233 - t308;
t78 = -t109 * t230 + t110 * t234;
t77 = t109 * t234 + t110 * t230;
t76 = -t265 * t154 + t247;
t75 = t150 + t89;
t74 = -t150 + t89;
t71 = t265 * t156 + t255;
t70 = -t156 * t274 + t233 * t89;
t69 = t156 * t273 + t229 * t89;
t68 = t154 * t273 - t229 * t88;
t67 = t154 * t274 + t233 * t88;
t66 = -pkin(3) * t128 + pkin(7) * t132 - t282;
t65 = -pkin(3) * t125 + pkin(7) * t117 + t281;
t64 = -t101 * t230 + t103 * t234;
t63 = -t100 * t230 + t102 * t234;
t62 = t101 * t234 + t103 * t230;
t61 = t100 * t234 + t102 * t230;
t60 = -t230 * t98 + t234 * t99;
t59 = t230 * t99 + t234 * t98;
t55 = -pkin(8) * t98 - t289;
t54 = -t231 * t91 + t235 * t93;
t53 = t231 * t93 + t235 * t91;
t52 = -t230 * t84 + t234 * t85;
t51 = t230 * t85 + t234 * t84;
t48 = -pkin(8) * t84 - t291;
t47 = t229 * t75 + t233 * t242;
t46 = -t229 * t74 - t233 * t71;
t45 = t229 * t242 - t233 * t75;
t44 = -t229 * t71 + t233 * t74;
t43 = pkin(4) * t45;
t42 = -t230 * t69 + t234 * t70;
t41 = -t230 * t67 + t234 * t68;
t40 = t230 * t70 + t234 * t69;
t39 = t230 * t68 + t234 * t67;
t38 = pkin(3) * t137 + pkin(7) * t50;
t37 = -pkin(7) * t91 - t49;
t36 = -pkin(4) * t76 + pkin(8) * t99 - t291;
t35 = -t231 * t59 + t235 * t60;
t34 = t231 * t60 + t235 * t59;
t33 = -pkin(4) * t71 + pkin(8) * t85 + t289;
t32 = -pkin(3) * t139 + pkin(7) * t93 + t50;
t28 = -t231 * t51 + t235 * t52;
t27 = t231 * t52 + t235 * t51;
t26 = t235 * t50 - t290;
t25 = t231 * t50 + t288;
t24 = -t230 * t45 + t234 * t47;
t23 = -t230 * t44 + t234 * t46;
t22 = t230 * t47 + t234 * t45;
t21 = t230 * t46 + t234 * t44;
t20 = -pkin(7) * t59 - t230 * t36 + t234 * t55;
t17 = pkin(4) * t18;
t16 = -pkin(7) * t51 - t230 * t33 + t234 * t48;
t15 = -pkin(3) * t76 + pkin(7) * t60 + t230 * t55 + t234 * t36;
t14 = pkin(4) * t79 + pkin(8) * t19;
t13 = -pkin(3) * t71 + pkin(7) * t52 + t230 * t48 + t234 * t33;
t12 = -pkin(8) * t45 - t18;
t11 = -t22 * t231 + t235 * t24;
t10 = t22 * t235 + t231 * t24;
t9 = -pkin(4) * t105 + pkin(8) * t47 + t19;
t8 = t19 * t234 - t293;
t7 = t19 * t230 + t292;
t6 = -pkin(7) * t22 + t12 * t234 - t230 * t9;
t5 = -pkin(3) * t105 + pkin(7) * t24 + t12 * t230 + t234 * t9;
t4 = -t231 * t7 + t235 * t8;
t3 = t231 * t8 + t235 * t7;
t2 = -pkin(7) * t7 - pkin(8) * t292 - t14 * t230;
t1 = pkin(3) * t79 + pkin(7) * t8 - pkin(8) * t293 + t14 * t234;
t29 = [0, 0, 0, 0, 0, qJDD(1), t258, t244, 0, 0, t218, 0.2e1 * t227 * t262, 0, t219, 0, 0, -t297 * t228, t297 * t227, pkin(1) * t211 + qJ(2) * (t219 + t218) + t253, pkin(1) * t203 + qJ(2) * t253, t227 * (t193 * t235 - t231 * t267) + t228 * (t193 * t231 + t235 * t267), t227 * (-t190 * t235 - t192 * t231) + t228 * (-t190 * t231 + t192 * t235), t227 * (-t198 * t231 + t303) + t228 * (t198 * t235 + t306), t227 * (-t191 * t231 + t235 * t268) + t228 * (t191 * t235 + t231 * t268), t227 * (t197 * t235 - t276) + t228 * (t197 * t231 + t275), (t227 * (-t207 * t235 + t209 * t231) + t228 * (-t207 * t231 - t209 * t235)) * qJD(3), t227 * (-pkin(6) * t163 - t278) + t228 * (-pkin(2) * t190 + pkin(6) * t164 + t277) - pkin(1) * t190 + qJ(2) * (-t163 * t227 + t164 * t228), t227 * (-pkin(6) * t167 - t277) + t228 * (-pkin(2) * t192 + pkin(6) * t168 - t278) - pkin(1) * t192 + qJ(2) * (-t167 * t227 + t168 * t228), t227 * (-pkin(6) * t165 - t106) + t228 * (-pkin(2) * t170 + pkin(6) * t166 + t107) - pkin(1) * t170 + qJ(2) * (-t165 * t227 + t166 * t228), -pkin(6) * t285 + t228 * (pkin(2) * t186 + pkin(6) * t107) + pkin(1) * t186 + qJ(2) * (t107 * t228 - t285), t227 * (-t123 * t231 + t124 * t235) + t228 * (t123 * t235 + t124 * t231), t227 * (-t231 * t90 + t235 * t92) + t228 * (t231 * t92 + t235 * t90), t227 * (-t133 * t231 + t135 * t235) + t228 * (t133 * t235 + t135 * t231), t227 * (-t121 * t231 + t122 * t235) + t228 * (t121 * t235 + t122 * t231), t227 * (-t134 * t231 + t136 * t235) + t228 * (t134 * t235 + t136 * t231), t227 * (-t146 * t231 + t147 * t235) + t228 * (t146 * t235 + t147 * t231), t227 * (-pkin(6) * t86 - t231 * t65 + t235 * t94) + t228 * (-pkin(2) * t125 + pkin(6) * t87 + t231 * t94 + t235 * t65) - pkin(1) * t125 + qJ(2) * (-t227 * t86 + t228 * t87), t227 * (-pkin(6) * t95 + t104 * t235 - t231 * t66) + t228 * (-pkin(2) * t128 + pkin(6) * t96 + t104 * t231 + t235 * t66) - pkin(1) * t128 + qJ(2) * (-t227 * t95 + t228 * t96), t227 * (-pkin(6) * t53 - t231 * t32 + t235 * t37) + t228 * (-pkin(2) * t139 + pkin(6) * t54 + t231 * t37 + t235 * t32) - pkin(1) * t139 + qJ(2) * (-t227 * t53 + t228 * t54), t227 * (-pkin(6) * t25 - pkin(7) * t288 - t231 * t38) + t228 * (pkin(2) * t137 + pkin(6) * t26 - pkin(7) * t290 + t235 * t38) + pkin(1) * t137 + qJ(2) * (-t227 * t25 + t228 * t26), t227 * (-t231 * t40 + t235 * t42) + t228 * (t231 * t42 + t235 * t40), t227 * (-t21 * t231 + t23 * t235) + t228 * (t21 * t235 + t23 * t231), t227 * (-t231 * t61 + t235 * t63) + t228 * (t231 * t63 + t235 * t61), t227 * (-t231 * t39 + t235 * t41) + t228 * (t231 * t41 + t235 * t39), t227 * (-t231 * t62 + t235 * t64) + t228 * (t231 * t64 + t235 * t62), t227 * (-t231 * t77 + t235 * t78) + t228 * (t231 * t78 + t235 * t77), t227 * (-pkin(6) * t27 - t13 * t231 + t16 * t235) + t228 * (-pkin(2) * t71 + pkin(6) * t28 + t13 * t235 + t16 * t231) - pkin(1) * t71 + qJ(2) * (-t227 * t27 + t228 * t28), t227 * (-pkin(6) * t34 - t15 * t231 + t20 * t235) + t228 * (-pkin(2) * t76 + pkin(6) * t35 + t15 * t235 + t20 * t231) - pkin(1) * t76 + qJ(2) * (-t227 * t34 + t228 * t35), t227 * (-pkin(6) * t10 - t231 * t5 + t235 * t6) + t228 * (-pkin(2) * t105 + pkin(6) * t11 + t231 * t6 + t235 * t5) - pkin(1) * t105 + qJ(2) * (-t10 * t227 + t11 * t228), t227 * (-pkin(6) * t3 - t1 * t231 + t2 * t235) + t228 * (pkin(2) * t79 + pkin(6) * t4 + t1 * t235 + t2 * t231) + pkin(1) * t79 + qJ(2) * (-t227 * t3 + t228 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t262, t263, -t211, -t203, 0, 0, 0, 0, 0, 0, t190, t192, t170, -t186, 0, 0, 0, 0, 0, 0, t125, t128, t139, -t137, 0, 0, 0, 0, 0, 0, t71, t76, t105, -t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t195, t205 - t204, t206, -t195, t169, qJDD(3), -t144, -t145, 0, 0, t162, t161, t130, -t162, t241, t223, pkin(3) * t116 - t81, pkin(3) * t131 - t82, pkin(3) * t91, pkin(3) * t49, t120, t118, t75, -t120, t242, t217, pkin(3) * t51 + t252, pkin(3) * t59 + t243, pkin(3) * t22 + t43, pkin(3) * t7 + t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t162, t161, t130, -t162, t241, t223, -t81, -t82, 0, 0, t120, t118, t75, -t120, t242, t217, t252, t243, t43, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, t118, t75, -t120, t242, t217, -t30, -t31, 0, 0;];
tauJ_reg = t29;
