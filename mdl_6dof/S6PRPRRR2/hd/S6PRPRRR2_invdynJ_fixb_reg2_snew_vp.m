% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRPRRR2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 00:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRPRRR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 00:24:58
% EndTime: 2019-05-05 00:25:07
% DurationCPUTime: 3.55s
% Computational Cost: add. (17028->367), mult. (32113->540), div. (0->0), fcn. (23741->14), ass. (0->230)
t201 = sin(qJ(6));
t202 = sin(qJ(5));
t206 = cos(qJ(5));
t203 = sin(qJ(4));
t240 = qJD(2) * t203;
t160 = -t206 * qJD(4) + t202 * t240;
t162 = t202 * qJD(4) + t206 * t240;
t205 = cos(qJ(6));
t139 = t205 * t160 + t201 * t162;
t141 = -t201 * t160 + t205 * t162;
t109 = t141 * t139;
t236 = qJD(2) * qJD(4);
t186 = t203 * t236;
t207 = cos(qJ(4));
t234 = t207 * qJDD(2);
t166 = -t186 + t234;
t159 = -qJDD(5) + t166;
t156 = -qJDD(6) + t159;
t265 = -t109 - t156;
t272 = t201 * t265;
t271 = t205 * t265;
t196 = sin(pkin(6));
t241 = -g(3) + qJDD(1);
t198 = cos(pkin(6));
t253 = sin(pkin(11));
t254 = cos(pkin(11));
t215 = t253 * g(1) - t254 * g(2);
t268 = t198 * t215;
t270 = t196 * t241 + t268;
t209 = qJD(2) ^ 2;
t183 = t207 * qJD(2) - qJD(5);
t177 = -qJD(6) + t183;
t125 = t139 * t177;
t231 = t207 * t236;
t235 = t203 * qJDD(2);
t165 = t231 + t235;
t220 = -t202 * qJDD(4) - t206 * t165;
t135 = -t160 * qJD(5) - t220;
t226 = -t206 * qJDD(4) + t202 * t165;
t218 = t162 * qJD(5) + t226;
t92 = -t139 * qJD(6) + t205 * t135 - t201 * t218;
t269 = t125 + t92;
t252 = t162 * t160;
t214 = -t159 - t252;
t267 = t202 * t214;
t266 = t206 * t214;
t153 = t160 * t183;
t114 = t135 - t153;
t227 = t201 * t135 + t205 * t218;
t65 = (qJD(6) + t177) * t141 + t227;
t110 = (qJD(5) + t183) * t162 + t226;
t137 = t139 ^ 2;
t138 = t141 ^ 2;
t264 = t160 ^ 2;
t158 = t162 ^ 2;
t176 = t177 ^ 2;
t181 = t183 ^ 2;
t263 = qJD(4) ^ 2;
t221 = t198 * t241 + qJDD(3);
t211 = -t196 * t215 + t221;
t210 = t203 * t211;
t225 = -t207 * pkin(4) - t203 * pkin(9);
t172 = -t254 * g(1) - t253 * g(2);
t204 = sin(qJ(2));
t208 = cos(qJ(2));
t132 = t208 * t172 + t270 * t204;
t124 = -t209 * pkin(2) + t132;
t195 = sin(pkin(12));
t197 = cos(pkin(12));
t131 = -t204 * t172 + t270 * t208;
t216 = qJDD(2) * pkin(2) + t131;
t97 = t197 * t124 + t195 * t216;
t95 = -t209 * pkin(3) + qJDD(2) * pkin(8) + t97;
t229 = t209 * t225 + t95;
t72 = -t263 * pkin(4) + qJDD(4) * pkin(9) + t229 * t207 + t210;
t223 = -t166 + t186;
t224 = t165 + t231;
t228 = t195 * t124 - t197 * t216;
t94 = -qJDD(2) * pkin(3) - t209 * pkin(8) + t228;
t79 = t223 * pkin(4) - t224 * pkin(9) + t94;
t45 = t202 * t72 - t206 * t79;
t34 = t214 * pkin(5) - t114 * pkin(10) - t45;
t150 = -t183 * pkin(5) - t162 * pkin(10);
t46 = t202 * t79 + t206 * t72;
t36 = -t264 * pkin(5) - pkin(10) * t218 + t183 * t150 + t46;
t15 = t201 * t36 - t205 * t34;
t16 = t201 * t34 + t205 * t36;
t9 = -t205 * t15 + t201 * t16;
t262 = pkin(5) * t9;
t68 = -t125 + t92;
t40 = -t201 * t65 - t205 * t68;
t261 = pkin(5) * t40;
t260 = t202 * t9;
t259 = t206 * t9;
t145 = t207 * t211;
t71 = -qJDD(4) * pkin(4) - t263 * pkin(9) + t229 * t203 - t145;
t52 = pkin(5) * t218 - t264 * pkin(10) + t162 * t150 + t71;
t258 = t201 * t52;
t257 = t202 * t71;
t256 = t205 * t52;
t255 = t206 * t71;
t251 = t177 * t201;
t250 = t177 * t205;
t249 = t183 * t202;
t100 = -t109 + t156;
t248 = t201 * t100;
t128 = t159 - t252;
t247 = t202 * t128;
t182 = t203 * t209 * t207;
t173 = qJDD(4) + t182;
t246 = t203 * t173;
t245 = t205 * t100;
t244 = t206 * t128;
t243 = t206 * t183;
t174 = qJDD(4) - t182;
t242 = t207 * t174;
t239 = qJD(5) - t183;
t233 = t207 * t109;
t232 = t207 * t252;
t10 = t201 * t15 + t205 * t16;
t25 = t202 * t45 + t206 * t46;
t87 = t203 * t95 - t145;
t88 = t207 * t95 + t210;
t51 = t203 * t87 + t207 * t88;
t24 = t202 * t46 - t206 * t45;
t222 = -pkin(3) + t225;
t105 = -t176 - t137;
t58 = t201 * t105 + t271;
t219 = pkin(5) * t58 - t15;
t116 = -t138 - t176;
t77 = t205 * t116 + t248;
t217 = pkin(5) * t77 - t16;
t192 = t207 ^ 2;
t191 = t203 ^ 2;
t190 = t192 * t209;
t188 = t191 * t209;
t180 = -t190 - t263;
t179 = -t188 - t263;
t171 = t188 + t190;
t170 = (t191 + t192) * qJDD(2);
t169 = -t195 * qJDD(2) - t197 * t209;
t168 = t197 * qJDD(2) - t195 * t209;
t167 = -0.2e1 * t186 + t234;
t164 = 0.2e1 * t231 + t235;
t152 = -t158 + t181;
t151 = -t181 + t264;
t149 = -t203 * t179 - t242;
t148 = t207 * t180 - t246;
t147 = -t203 * t174 + t207 * t179;
t146 = t207 * t173 + t203 * t180;
t144 = t158 - t264;
t143 = -t158 - t181;
t142 = t195 * t170 + t197 * t171;
t136 = -t181 - t264;
t127 = t158 + t264;
t123 = -t138 + t176;
t122 = t137 - t176;
t121 = t195 * t149 - t197 * t164;
t120 = t195 * t148 + t197 * t167;
t115 = t239 * t160 + t220;
t113 = t135 + t153;
t111 = -t239 * t162 - t226;
t108 = t138 - t137;
t107 = -t202 * t143 + t244;
t106 = t206 * t143 + t247;
t104 = t206 * t136 - t267;
t103 = t202 * t136 + t266;
t99 = (t139 * t205 - t141 * t201) * t177;
t98 = (t139 * t201 + t141 * t205) * t177;
t93 = -t137 - t138;
t91 = -t141 * qJD(6) - t227;
t90 = -t110 * t206 + t202 * t114;
t89 = -t110 * t202 - t206 * t114;
t86 = t205 * t122 + t248;
t85 = -t201 * t123 + t271;
t84 = t201 * t122 - t245;
t83 = t205 * t123 + t272;
t81 = t207 * t107 - t203 * t115;
t80 = t203 * t107 + t207 * t115;
t78 = -t201 * t116 + t245;
t74 = t207 * t104 - t203 * t111;
t73 = t203 * t104 + t207 * t111;
t64 = (qJD(6) - t177) * t141 + t227;
t63 = t141 * t251 + t205 * t92;
t62 = -t141 * t250 + t201 * t92;
t61 = -t139 * t250 - t201 * t91;
t60 = -t139 * t251 + t205 * t91;
t59 = t205 * t105 - t272;
t57 = -t203 * t127 + t207 * t90;
t56 = t207 * t127 + t203 * t90;
t55 = t195 * t97 - t197 * t228;
t54 = -t197 * t106 + t195 * t81;
t53 = -t197 * t103 + t195 * t74;
t50 = t203 * t88 - t207 * t87;
t49 = -t202 * t77 + t206 * t78;
t48 = t202 * t78 + t206 * t77;
t47 = t195 * t57 - t197 * t89;
t43 = -t201 * t269 - t205 * t64;
t42 = t201 * t68 - t205 * t65;
t41 = -t201 * t64 + t205 * t269;
t39 = -t202 * t58 + t206 * t59;
t38 = t202 * t59 + t206 * t58;
t37 = t195 * t51 - t197 * t94;
t35 = -pkin(10) * t77 + t256;
t32 = -pkin(10) * t58 + t258;
t31 = t203 * t269 + t207 * t49;
t30 = t203 * t49 - t207 * t269;
t29 = t203 * t64 + t207 * t39;
t28 = t203 * t39 - t207 * t64;
t27 = -pkin(5) * t269 + pkin(10) * t78 + t258;
t26 = -pkin(5) * t64 + pkin(10) * t59 - t256;
t23 = -t202 * t40 + t206 * t42;
t22 = t202 * t42 + t206 * t40;
t21 = t195 * t31 - t197 * t48;
t20 = t203 * t93 + t207 * t23;
t19 = t203 * t23 - t207 * t93;
t18 = t203 * t71 + t207 * t25;
t17 = t203 * t25 - t207 * t71;
t13 = t195 * t29 - t197 * t38;
t12 = t195 * t18 - t197 * t24;
t11 = t195 * t20 - t197 * t22;
t8 = -pkin(5) * t52 + pkin(10) * t10;
t7 = -pkin(10) * t40 - t9;
t6 = -pkin(5) * t93 + pkin(10) * t42 + t10;
t5 = t206 * t10 - t260;
t4 = t202 * t10 + t259;
t3 = t203 * t52 + t207 * t5;
t2 = t203 * t5 - t207 * t52;
t1 = t195 * t3 - t197 * t4;
t14 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t241, 0, 0, 0, 0, 0, 0, (qJDD(2) * t208 - t204 * t209) * t196, (-qJDD(2) * t204 - t208 * t209) * t196, 0, t198 ^ 2 * t241 + (t208 * t131 + t204 * t132 - t268) * t196, 0, 0, 0, 0, 0, 0, (t168 * t208 + t169 * t204) * t196, (-t168 * t204 + t169 * t208) * t196, 0, t198 * t221 + (t204 * (t195 * t228 + t197 * t97) + t208 * t55 - t268) * t196, 0, 0, 0, 0, 0, 0, t198 * t146 + (t204 * (t197 * t148 - t195 * t167) + t208 * t120) * t196, t198 * t147 + (t204 * (t197 * t149 + t195 * t164) + t208 * t121) * t196, (t204 * (t197 * t170 - t195 * t171) + t208 * t142) * t196, t198 * t50 + (t204 * (t195 * t94 + t197 * t51) + t208 * t37) * t196, 0, 0, 0, 0, 0, 0, t198 * t73 + (t204 * (t195 * t103 + t197 * t74) + t208 * t53) * t196, t198 * t80 + (t204 * (t195 * t106 + t197 * t81) + t208 * t54) * t196, t198 * t56 + (t204 * (t195 * t89 + t197 * t57) + t208 * t47) * t196, t198 * t17 + (t204 * (t197 * t18 + t195 * t24) + t208 * t12) * t196, 0, 0, 0, 0, 0, 0, t198 * t28 + (t204 * (t195 * t38 + t197 * t29) + t208 * t13) * t196, t198 * t30 + (t204 * (t195 * t48 + t197 * t31) + t208 * t21) * t196, t198 * t19 + (t204 * (t195 * t22 + t197 * t20) + t208 * t11) * t196, t198 * t2 + (t204 * (t195 * t4 + t197 * t3) + t208 * t1) * t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t131, -t132, 0, 0, 0, 0, 0, 0, 0, qJDD(2), pkin(2) * t168 - t228, pkin(2) * t169 - t97, 0, pkin(2) * t55, t224 * t203, t207 * t164 + t203 * t167, t246 + t207 * (-t188 + t263), -t223 * t207, t203 * (t190 - t263) + t242, 0, pkin(2) * t120 + pkin(3) * t167 + pkin(8) * t148 - t207 * t94, pkin(2) * t121 - pkin(3) * t164 + pkin(8) * t149 + t203 * t94, pkin(2) * t142 + pkin(3) * t171 + pkin(8) * t170 + t51, pkin(2) * t37 - pkin(3) * t94 + pkin(8) * t51, t203 * (t206 * t135 + t162 * t249) - t232, t203 * (t206 * t111 - t202 * t113) - t207 * t144, t203 * (-t202 * t152 + t266) - t207 * t114, t203 * (-t160 * t243 + t202 * t218) + t232, t203 * (t206 * t151 + t247) + t207 * t110, t207 * t159 + t203 * (t160 * t206 - t162 * t202) * t183, t203 * (-pkin(9) * t103 + t257) + t207 * (-pkin(4) * t103 + t45) - pkin(3) * t103 + pkin(8) * t74 + pkin(2) * t53, t203 * (-pkin(9) * t106 + t255) + t207 * (-pkin(4) * t106 + t46) - pkin(3) * t106 + pkin(8) * t81 + pkin(2) * t54, pkin(2) * t47 + pkin(8) * t57 - t203 * t24 + t222 * t89, pkin(2) * t12 + pkin(8) * t18 + t222 * t24, t203 * (-t202 * t62 + t206 * t63) - t233, t203 * (-t202 * t41 + t206 * t43) - t207 * t108, t203 * (-t202 * t83 + t206 * t85) - t207 * t68, t203 * (-t202 * t60 + t206 * t61) + t233, t203 * (-t202 * t84 + t206 * t86) + t207 * t65, t203 * (-t202 * t98 + t206 * t99) + t207 * t156, t203 * (-pkin(9) * t38 - t202 * t26 + t206 * t32) + t207 * (-pkin(4) * t38 - t219) - pkin(3) * t38 + pkin(8) * t29 + pkin(2) * t13, t203 * (-pkin(9) * t48 - t202 * t27 + t206 * t35) + t207 * (-pkin(4) * t48 - t217) - pkin(3) * t48 + pkin(8) * t31 + pkin(2) * t21, t203 * (-pkin(9) * t22 - t202 * t6 + t206 * t7) + t207 * (-pkin(4) * t22 - t261) - pkin(3) * t22 + pkin(8) * t20 + pkin(2) * t11, t203 * (-pkin(9) * t4 - pkin(10) * t259 - t202 * t8) + t207 * (-pkin(4) * t4 - t262) - pkin(3) * t4 + pkin(8) * t3 + pkin(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t211, 0, 0, 0, 0, 0, 0, t146, t147, 0, t50, 0, 0, 0, 0, 0, 0, t73, t80, t56, t17, 0, 0, 0, 0, 0, 0, t28, t30, t19, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t182, t188 - t190, t235, t182, t234, qJDD(4), -t87, -t88, 0, 0, t202 * t135 - t162 * t243, t202 * t111 + t206 * t113, t206 * t152 + t267, -t160 * t249 - t206 * t218, t202 * t151 - t244, (t160 * t202 + t162 * t206) * t183, pkin(4) * t111 + pkin(9) * t104 - t255, pkin(4) * t115 + pkin(9) * t107 + t257, pkin(4) * t127 + pkin(9) * t90 + t25, -pkin(4) * t71 + pkin(9) * t25, t202 * t63 + t206 * t62, t202 * t43 + t206 * t41, t202 * t85 + t206 * t83, t202 * t61 + t206 * t60, t202 * t86 + t206 * t84, t202 * t99 + t206 * t98, -pkin(4) * t64 + pkin(9) * t39 + t202 * t32 + t206 * t26, -pkin(4) * t269 + pkin(9) * t49 + t202 * t35 + t206 * t27, -pkin(4) * t93 + pkin(9) * t23 + t202 * t7 + t206 * t6, -pkin(4) * t52 + pkin(9) * t5 - pkin(10) * t260 + t206 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t252, t144, t114, -t252, -t110, -t159, -t45, -t46, 0, 0, t109, t108, t68, -t109, -t65, -t156, t219, t217, t261, t262; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, t108, t68, -t109, -t65, -t156, -t15, -t16, 0, 0;];
tauJ_reg  = t14;
