% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRRP7
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP7_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP7_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP7_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP7_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:45:42
% EndTime: 2019-12-31 18:45:48
% DurationCPUTime: 3.04s
% Computational Cost: add. (2946->388), mult. (6107->484), div. (0->0), fcn. (3774->10), ass. (0->197)
t130 = sin(qJ(4));
t133 = cos(qJ(4));
t134 = cos(qJ(3));
t208 = t134 * qJD(1);
t103 = -qJD(4) + t208;
t211 = t130 * qJD(3);
t131 = sin(qJ(3));
t217 = qJD(1) * t131;
t80 = t133 * t217 + t211;
t229 = t80 * t103;
t209 = t133 * qJD(3);
t78 = t130 * t217 - t209;
t230 = t78 * t103;
t197 = t134 * t209;
t205 = t131 * qJDD(1);
t213 = qJD(4) * t131;
t263 = qJD(1) * t213 - qJDD(3);
t36 = -qJD(1) * t197 - qJD(4) * t209 + t263 * t130 - t133 * t205;
t37 = ((qJD(4) + t208) * qJD(3) + t205) * t130 + t263 * t133;
t268 = (t36 - t230) * t133 + (t37 - t229) * t130;
t128 = sin(pkin(8));
t109 = t128 * pkin(1) + pkin(6);
t90 = t109 * qJDD(1);
t267 = -qJD(2) * qJD(3) - t90;
t264 = t37 + t229;
t125 = qJ(1) + pkin(8);
t114 = sin(t125);
t115 = cos(t125);
t179 = g(1) * t115 + g(2) * t114;
t228 = pkin(1) * qJDD(1);
t262 = t131 * t179;
t219 = t134 * pkin(3) + t131 * pkin(7);
t92 = t109 * qJD(1);
t77 = t131 * t92;
t63 = t134 * qJD(2) - t77;
t51 = -qJD(3) * pkin(3) - t63;
t18 = t78 * pkin(4) - t80 * qJ(5) + t51;
t120 = t134 * qJDD(1);
t207 = qJD(1) * qJD(3);
t191 = t131 * t207;
t74 = qJDD(4) - t120 + t191;
t255 = pkin(7) * t74;
t261 = t103 * t18 + t255;
t212 = qJD(4) * t133;
t194 = t131 * t212;
t198 = t134 * t211;
t224 = t130 * t131;
t260 = -t103 * (t194 + t198) + t74 * t224;
t214 = qJD(4) * t130;
t196 = t103 * t214;
t223 = t130 * t134;
t259 = qJD(1) * (t103 * t223 - t131 * t78) - t133 * t74 - t196;
t221 = t133 * t134;
t129 = cos(pkin(8));
t245 = t129 * pkin(1);
t110 = -pkin(2) - t245;
t71 = t110 - t219;
t239 = t109 * t221 + t130 * t71;
t252 = pkin(7) * t134;
t180 = pkin(3) * t131 - t252;
t84 = t180 * qJD(3);
t257 = qJD(4) * t239 - t133 * t84;
t256 = t80 ^ 2;
t254 = pkin(7) * t80;
t253 = t74 * pkin(4);
t251 = g(1) * t114;
t248 = g(2) * t115;
t247 = g(3) * t131;
t246 = g(3) * t134;
t244 = t80 * t78;
t222 = t131 * t133;
t54 = t78 * t197;
t243 = -t37 * t222 - t54;
t83 = t180 * qJD(1);
t35 = t130 * t83 + t133 * t63;
t61 = t74 * t222;
t242 = -t103 * t197 + t61;
t241 = t130 * t84 + t71 * t212;
t174 = pkin(4) * t130 - qJ(5) * t133;
t210 = t131 * qJD(2);
t240 = t174 * qJD(4) - t130 * qJD(5) - t210 - (t174 * qJD(1) + t92) * t134;
t238 = t179 * t222;
t237 = t133 * t71;
t235 = t134 * t36;
t64 = t134 * t92 + t210;
t52 = qJD(3) * pkin(7) + t64;
t167 = -pkin(2) - t219;
t152 = t167 - t245;
t53 = t152 * qJD(1);
t20 = t130 * t53 + t133 * t52;
t15 = -t103 * qJ(5) + t20;
t234 = t15 * t103;
t233 = t20 * t103;
t232 = t37 * t133;
t231 = t74 * qJ(5);
t227 = t114 * t131;
t226 = t115 * t131;
t225 = t115 * t134;
t19 = -t130 * t52 + t133 * t53;
t220 = qJD(5) - t19;
t126 = t131 ^ 2;
t127 = t134 ^ 2;
t218 = t126 - t127;
t93 = qJD(1) * t110;
t216 = qJD(3) * t131;
t215 = qJD(3) * t134;
t91 = qJDD(1) * t110;
t56 = t80 * t194;
t204 = t80 * t198 - t36 * t224 + t56;
t202 = -t131 * qJDD(2) + t267 * t134;
t201 = t78 ^ 2 - t256;
t137 = qJD(1) ^ 2;
t200 = t131 * t137 * t134;
t135 = cos(qJ(1));
t199 = t135 * pkin(1) + t115 * pkin(2) + t114 * pkin(6);
t195 = t130 * t213;
t193 = t103 * t217;
t132 = sin(qJ(1));
t190 = -t132 * pkin(1) + t115 * pkin(6);
t189 = t109 * t130 + pkin(4);
t188 = t80 * t216 + t235;
t187 = -t134 * t37 + t78 * t216;
t30 = -t92 * t216 - t202;
t27 = qJDD(3) * pkin(7) + t30;
t38 = qJD(1) * t84 + t152 * qJDD(1);
t186 = t130 * t27 - t133 * t38 + t52 * t212 + t53 * t214;
t185 = qJD(4) * t78 - t36;
t31 = t134 * qJDD(2) + t267 * t131 - t92 * t215;
t183 = t134 * t191;
t57 = t114 * t223 + t115 * t133;
t59 = -t114 * t133 + t115 * t223;
t182 = -g(1) * t57 + g(2) * t59;
t58 = t114 * t221 - t115 * t130;
t60 = t114 * t130 + t115 * t221;
t181 = g(1) * t58 - g(2) * t60;
t178 = g(1) * t132 - g(2) * t135;
t177 = t185 * pkin(7);
t176 = pkin(3) * t225 + pkin(7) * t226 + t199;
t175 = t133 * pkin(4) + t130 * qJ(5);
t14 = t103 * pkin(4) + t220;
t173 = -t130 * t15 + t133 * t14;
t172 = t130 * t14 + t133 * t15;
t171 = -t130 * t20 - t133 * t19;
t170 = t130 * t19 - t133 * t20;
t34 = -t130 * t63 + t133 * t83;
t166 = pkin(3) + t175;
t165 = pkin(7) * qJD(4) * t103 - t246;
t163 = t109 + t174;
t160 = -t103 * t212 + t130 * t74;
t3 = t130 * t38 + t133 * t27 + t53 * t212 - t52 * t214;
t158 = t78 * t195 + t243;
t28 = -qJDD(3) * pkin(3) - t31;
t5 = t37 * pkin(4) + t36 * qJ(5) - t80 * qJD(5) + t28;
t157 = t165 - t5;
t156 = t165 - t28;
t155 = -qJD(1) * t93 + t179;
t153 = t167 * t251;
t151 = -t103 * t51 - t255;
t136 = qJD(3) ^ 2;
t150 = t109 * t136 + t248 + 0.2e1 * t91;
t149 = 0.2e1 * t93 * qJD(3) - qJDD(3) * t109;
t148 = g(1) * t59 + g(2) * t57 + g(3) * t224 - t186;
t147 = -t130 * t230 - t232;
t146 = -pkin(7) * t232 - t179 * t134 - t247;
t1 = -t103 * qJD(5) + t231 + t3;
t2 = qJDD(5) + t186 - t253;
t145 = t173 * qJD(4) + t1 * t133 + t2 * t130;
t144 = t171 * qJD(4) + t130 * t186 + t3 * t133;
t143 = -t31 * t131 + t30 * t134 + (-t131 * t64 - t134 * t63) * qJD(3);
t142 = t18 * t80 + qJDD(5) - t148;
t141 = t78 * t194 + (t131 * t37 + t78 * t215) * t130;
t140 = -g(1) * t60 - g(2) * t58 - g(3) * t222 + t3;
t139 = t187 - t260;
t98 = g(1) * t227;
t96 = pkin(7) * t225;
t94 = t114 * t252;
t89 = qJDD(3) * t134 - t136 * t131;
t88 = qJDD(3) * t131 + t136 * t134;
t49 = t163 * t131;
t45 = -t103 * t216 - t74 * t134;
t44 = t80 * pkin(4) + t78 * qJ(5);
t42 = -t109 * t223 + t237;
t41 = t189 * t134 - t237;
t40 = -t134 * qJ(5) + t239;
t24 = -pkin(4) * t217 - t34;
t23 = qJ(5) * t217 + t35;
t17 = (t175 * qJD(4) - qJD(5) * t133) * t131 + t163 * t215;
t16 = -t36 - t230;
t13 = (t103 * t221 - t131 * t80) * qJD(1) + t160;
t12 = t131 * t109 * t211 - t257;
t11 = (-t131 * t209 - t134 * t214) * t109 + t241;
t10 = -t189 * t216 + t257;
t9 = (-t109 * t214 - qJD(5)) * t134 + (-t109 * t133 + qJ(5)) * t216 + t241;
t8 = -t36 * t130 - t133 * t229;
t7 = t80 * t197 + (-t36 * t133 - t80 * t214) * t131;
t6 = t103 * t195 + t188 + t242;
t4 = [0, 0, 0, 0, 0, qJDD(1), t178, g(1) * t135 + g(2) * t132, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t129 * t228 - t248 + t251, -0.2e1 * t128 * t228 + t179, 0, (t178 + (t128 ^ 2 + t129 ^ 2) * t228) * pkin(1), t126 * qJDD(1) + 0.2e1 * t183, 0.2e1 * t131 * t120 - 0.2e1 * t218 * t207, t88, t127 * qJDD(1) - 0.2e1 * t183, t89, 0, t149 * t131 + (-t150 + t251) * t134, t131 * t150 + t134 * t149 - t98, (t126 + t127) * t90 + t143 - t179, t91 * t110 - g(1) * (-t114 * pkin(2) + t190) - g(2) * t199 + t143 * t109, t7, t158 - t204, t6, t141, (t103 * t211 + t37) * t134 + (-qJD(3) * t78 - t160) * t131, t45, -t12 * t103 + t42 * t74 + (t186 + (t109 * t78 + t130 * t51) * qJD(3)) * t134 + (qJD(3) * t19 + t109 * t37 + t28 * t130 + t212 * t51) * t131 + t181, t11 * t103 - t239 * t74 + (t3 + (t109 * t80 + t133 * t51) * qJD(3)) * t134 + (-qJD(3) * t20 - t109 * t36 + t28 * t133 - t214 * t51) * t131 + t182, -t11 * t78 - t12 * t80 + t42 * t36 - t239 * t37 + t98 + t171 * t215 + (qJD(4) * t170 - t130 * t3 + t133 * t186 - t248) * t131, t3 * t239 + t20 * t11 - t186 * t42 + t19 * t12 - g(1) * t190 - g(2) * t176 - t153 + (t28 * t131 + t215 * t51) * t109, t7, t6, t54 + (-t214 * t78 + t232) * t131 + t204, t45, t187 + t260, t141, t10 * t103 + t17 * t78 + t49 * t37 - t41 * t74 + (t18 * t211 + t2) * t134 + (-qJD(3) * t14 + t5 * t130 + t18 * t212) * t131 + t181, t10 * t80 - t41 * t36 - t40 * t37 - t9 * t78 + t98 + t173 * t215 + (-qJD(4) * t172 - t1 * t130 + t133 * t2 - t248) * t131, -t9 * t103 - t17 * t80 + t49 * t36 + t40 * t74 + (-t18 * t209 - t1) * t134 + (qJD(3) * t15 - t5 * t133 + t18 * t214) * t131 - t182, t1 * t40 + t15 * t9 + t5 * t49 + t18 * t17 + t2 * t41 + t14 * t10 - g(1) * (-t58 * pkin(4) - t57 * qJ(5) + t190) - g(2) * (t60 * pkin(4) + t59 * qJ(5) + t176) - t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, t89, -t88, 0, t30 * t131 + t31 * t134 - g(3) + (-t131 * t63 + t134 * t64) * qJD(3), 0, 0, 0, 0, 0, 0, t139, -t61 + (-t195 + t197) * t103 + t188, t56 + (t131 * t185 + t215 * t80) * t130 + t243, -g(3) + (-qJD(3) * t170 - t28) * t134 + (qJD(3) * t51 + t144) * t131, 0, 0, 0, 0, 0, 0, t139, t158 + t204, -t235 + (-qJD(3) * t80 + t196) * t131 + t242, -g(3) + (qJD(3) * t172 - t5) * t134 + (qJD(3) * t18 + t145) * t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t200, t218 * t137, t205, t200, t120, qJDD(3), t64 * qJD(3) + t131 * t155 - t246 + t31, t247 + (t63 + t77) * qJD(3) + t155 * t134 + t202, 0, 0, t8, -t268, t13, t147, -t259, t193, -pkin(3) * t37 + t34 * t103 + t130 * t151 + t133 * t156 - t19 * t217 - t64 * t78 + t238, t20 * t217 + pkin(3) * t36 - t35 * t103 - t64 * t80 + t151 * t133 + (-t156 - t262) * t130, t34 * t80 + t35 * t78 + (t19 * t208 + t3 + (-t19 + t254) * qJD(4)) * t133 + (t177 + t186 + t233) * t130 + t146, -t28 * pkin(3) - t20 * t35 - t19 * t34 - t51 * t64 - g(1) * (-pkin(3) * t226 + t96) - g(2) * (-pkin(3) * t227 + t94) - g(3) * t219 + t144 * pkin(7), t8, t13, t268, t193, t259, t147, -t24 * t103 - t261 * t130 + t157 * t133 + t14 * t217 - t166 * t37 + t240 * t78 + t238, t23 * t78 - t24 * t80 + (-t14 * t208 + t1 + (t14 + t254) * qJD(4)) * t133 + (t177 + t2 + t234) * t130 + t146, -t15 * t217 + t23 * t103 - t166 * t36 - t240 * t80 + t261 * t133 + (t157 + t262) * t130, -t15 * t23 - t14 * t24 - g(1) * t96 - g(2) * t94 - g(3) * (t134 * t175 + t219) + t240 * t18 + t145 * pkin(7) + (-t5 + t262) * t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t244, -t201, t16, -t244, -t264, t74, -t51 * t80 + t148 - t233, -t19 * t103 + t51 * t78 - t140, 0, 0, t244, t16, t201, t74, t264, -t244, -t44 * t78 - t142 - t233 + 0.2e1 * t253, pkin(4) * t36 - t37 * qJ(5) + (t15 - t20) * t80 + (t14 - t220) * t78, 0.2e1 * t231 - t18 * t78 + t44 * t80 + (-0.2e1 * qJD(5) + t19) * t103 + t140, t1 * qJ(5) - t2 * pkin(4) - t18 * t44 - t14 * t20 - g(1) * (-t59 * pkin(4) + t60 * qJ(5)) - g(2) * (-t57 * pkin(4) + t58 * qJ(5)) + t220 * t15 + t174 * t247; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74 + t244, t16, -t103 ^ 2 - t256, t142 + t234 - t253;];
tau_reg = t4;
