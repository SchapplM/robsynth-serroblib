% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
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
% tau_reg [5x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:33:08
% EndTime: 2019-12-31 19:33:17
% DurationCPUTime: 2.35s
% Computational Cost: add. (3175->335), mult. (7561->464), div. (0->0), fcn. (5617->14), ass. (0->171)
t149 = qJ(2) + pkin(8);
t144 = sin(t149);
t146 = cos(t149);
t158 = sin(qJ(1));
t161 = cos(qJ(1));
t183 = g(1) * t161 + g(2) * t158;
t168 = -g(3) * t146 + t144 * t183;
t153 = sin(pkin(8));
t209 = cos(pkin(8));
t157 = sin(qJ(2));
t160 = cos(qJ(2));
t218 = qJ(3) + pkin(6);
t187 = qJD(2) * t218;
t104 = -t157 * qJD(3) - t160 * t187;
t129 = t218 * t157;
t74 = qJDD(2) * pkin(2) + qJD(1) * t104 - qJDD(1) * t129;
t103 = t160 * qJD(3) - t157 * t187;
t130 = t218 * t160;
t83 = qJD(1) * t103 + qJDD(1) * t130;
t33 = -t153 * t83 + t209 * t74;
t32 = -qJDD(2) * pkin(3) + qJDD(4) - t33;
t232 = t32 - t168;
t185 = t209 * t160;
t134 = qJD(1) * t185;
t198 = qJD(1) * t157;
t106 = t153 * t198 - t134;
t102 = qJD(5) + t106;
t156 = sin(qJ(5));
t159 = cos(qJ(5));
t122 = t153 * t160 + t157 * t209;
t109 = t122 * qJD(1);
t152 = sin(pkin(9));
t154 = cos(pkin(9));
t93 = qJD(2) * t152 + t109 * t154;
t94 = qJD(2) * t154 - t109 * t152;
t41 = t156 * t93 - t159 * t94;
t234 = t41 * t102;
t123 = t152 * t159 + t154 * t156;
t113 = t123 * qJD(5);
t210 = t123 * t106 + t113;
t175 = -t156 * t94 - t159 * t93;
t233 = t102 * t175;
t230 = g(1) * t158 - g(2) * t161;
t231 = t146 * t230;
t229 = -qJD(5) + t102;
t121 = t152 * t156 - t154 * t159;
t211 = t102 * t121;
t108 = t122 * qJD(2);
t194 = t157 * qJDD(1);
t80 = qJD(1) * t108 - qJDD(1) * t185 + t153 * t194;
t75 = qJDD(5) + t80;
t228 = t102 * t211 - t123 * t75;
t195 = qJD(1) * qJD(2);
t190 = t157 * t195;
t81 = qJD(2) * t134 + qJDD(1) * t122 - t153 * t190;
t59 = -qJDD(2) * t154 + t152 * t81;
t60 = qJDD(2) * t152 + t154 * t81;
t9 = -qJD(5) * t175 + t156 * t60 + t159 * t59;
t105 = t106 ^ 2;
t224 = g(3) * t144;
t222 = g(3) * t160;
t221 = t154 * pkin(7);
t220 = t160 * pkin(2);
t136 = pkin(2) * t153 + qJ(4);
t219 = pkin(7) + t136;
t140 = pkin(1) + t220;
t169 = pkin(2) * t190 - qJDD(1) * t140 + qJDD(3);
t18 = t80 * pkin(3) - t81 * qJ(4) - t109 * qJD(4) + t169;
t34 = t153 * t74 + t209 * t83;
t29 = qJDD(2) * qJ(4) + qJD(2) * qJD(4) + t34;
t7 = t152 * t18 + t154 * t29;
t171 = -t153 * t157 + t185;
t111 = t171 * qJD(2);
t217 = qJD(2) * pkin(2);
t192 = t157 * t217;
t44 = pkin(3) * t108 - qJ(4) * t111 - qJD(4) * t122 + t192;
t64 = t103 * t209 + t104 * t153;
t20 = t152 * t44 + t154 * t64;
t128 = -qJD(1) * t140 + qJD(3);
t52 = t106 * pkin(3) - t109 * qJ(4) + t128;
t125 = qJD(1) * t129;
t119 = -t125 + t217;
t126 = qJD(1) * t130;
t186 = t209 * t126;
t79 = t119 * t153 + t186;
t71 = qJD(2) * qJ(4) + t79;
t24 = t152 * t52 + t154 * t71;
t62 = pkin(2) * t198 + pkin(3) * t109 + qJ(4) * t106;
t114 = t153 * t126;
t85 = -t125 * t209 - t114;
t31 = t152 * t62 + t154 * t85;
t77 = -pkin(3) * t171 - qJ(4) * t122 - t140;
t90 = -t129 * t153 + t130 * t209;
t36 = t152 * t77 + t154 * t90;
t216 = t109 * t41;
t214 = t152 * t80;
t213 = t154 * t80;
t212 = t175 * t109;
t208 = t106 * t152;
t207 = t111 * t152;
t206 = t122 * t152;
t205 = t122 * t154;
t148 = pkin(9) + qJ(5);
t143 = sin(t148);
t204 = t158 * t143;
t145 = cos(t148);
t203 = t158 * t145;
t202 = t161 * t143;
t201 = t161 * t145;
t78 = t119 * t209 - t114;
t65 = -qJD(2) * pkin(3) + qJD(4) - t78;
t200 = -qJD(4) + t65;
t150 = t157 ^ 2;
t199 = -t160 ^ 2 + t150;
t197 = qJD(5) * t156;
t196 = qJD(5) * t159;
t193 = t160 * qJDD(1);
t6 = -t152 * t29 + t154 * t18;
t2 = pkin(4) * t80 - pkin(7) * t60 + t6;
t5 = -pkin(7) * t59 + t7;
t191 = -t156 * t5 + t159 * t2;
t19 = -t152 * t64 + t154 * t44;
t23 = -t152 * t71 + t154 * t52;
t30 = -t152 * t85 + t154 * t62;
t35 = -t152 * t90 + t154 * t77;
t63 = t103 * t153 - t104 * t209;
t84 = -t125 * t153 + t186;
t89 = t129 * t209 + t130 * t153;
t184 = -t102 * t210 - t121 * t75;
t139 = -pkin(2) * t209 - pkin(3);
t181 = -t7 * t152 - t6 * t154;
t180 = t156 * t2 + t159 * t5;
t179 = pkin(3) * t146 + qJ(4) * t144;
t11 = pkin(4) * t106 - pkin(7) * t93 + t23;
t15 = pkin(7) * t94 + t24;
t3 = t11 * t159 - t15 * t156;
t4 = t11 * t156 + t15 * t159;
t178 = -t152 * t23 + t154 * t24;
t21 = -pkin(4) * t171 - pkin(7) * t205 + t35;
t25 = -pkin(7) * t206 + t36;
t177 = -t156 * t25 + t159 * t21;
t176 = t156 * t21 + t159 * t25;
t174 = -0.2e1 * pkin(1) * t195 - pkin(6) * qJDD(2);
t8 = -t156 * t59 + t159 * t60 + t196 * t94 - t197 * t93;
t117 = t219 * t152;
t173 = pkin(7) * t208 - qJD(4) * t154 + qJD(5) * t117 + t31;
t118 = t219 * t154;
t172 = pkin(4) * t109 + qJD(4) * t152 + qJD(5) * t118 + t106 * t221 + t30;
t162 = qJD(2) ^ 2;
t166 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t162 + t230;
t163 = qJD(1) ^ 2;
t165 = pkin(1) * t163 - pkin(6) * qJDD(1) + t183;
t164 = t111 * t65 + t122 * t32 - t183;
t132 = t161 * t140;
t127 = -pkin(4) * t154 + t139;
t98 = t146 * t201 + t204;
t97 = -t146 * t202 + t203;
t96 = -t146 * t203 + t202;
t95 = t146 * t204 + t201;
t67 = t121 * t122;
t66 = t123 * t122;
t58 = pkin(4) * t206 + t89;
t45 = -pkin(4) * t208 + t84;
t38 = pkin(4) * t207 + t63;
t37 = -pkin(4) * t94 + t65;
t28 = t111 * t123 + t196 * t205 - t197 * t206;
t27 = -t111 * t121 - t113 * t122;
t13 = pkin(4) * t59 + t32;
t12 = -pkin(7) * t207 + t20;
t10 = pkin(4) * t108 - t111 * t221 + t19;
t1 = [qJDD(1), t230, t183, qJDD(1) * t150 + 0.2e1 * t160 * t190, 0.2e1 * t157 * t193 - 0.2e1 * t195 * t199, qJDD(2) * t157 + t160 * t162, qJDD(2) * t160 - t157 * t162, 0, t157 * t174 + t160 * t166, -t157 * t166 + t160 * t174, -t106 * t64 - t108 * t79 + t109 * t63 - t111 * t78 - t122 * t33 + t171 * t34 - t80 * t90 + t81 * t89 - t183, t34 * t90 + t79 * t64 - t33 * t89 - t78 * t63 - t169 * t140 + t128 * t192 - g(1) * (-t140 * t158 + t161 * t218) - g(2) * (t158 * t218 + t132), t19 * t106 + t23 * t108 + t152 * t164 + t154 * t231 - t171 * t6 + t35 * t80 + t89 * t59 - t63 * t94, -t20 * t106 - t24 * t108 - t152 * t231 + t154 * t164 + t171 * t7 - t36 * t80 + t89 * t60 + t63 * t93, -t19 * t93 + t20 * t94 - t35 * t60 - t36 * t59 + t230 * t144 + t181 * t122 + (-t152 * t24 - t154 * t23) * t111, -g(2) * t132 + t23 * t19 + t24 * t20 + t32 * t89 + t6 * t35 + t7 * t36 + t65 * t63 + (-g(1) * t218 - g(2) * t179) * t161 + (-g(1) * (-t140 - t179) - g(2) * t218) * t158, -t175 * t27 - t67 * t8, t175 * t28 - t27 * t41 - t66 * t8 + t67 * t9, t102 * t27 - t108 * t175 - t171 * t8 - t67 * t75, -t102 * t28 - t108 * t41 + t171 * t9 - t66 * t75, t102 * t108 - t171 * t75, (t159 * t10 - t156 * t12) * t102 + t177 * t75 - t191 * t171 + t3 * t108 + t38 * t41 + t58 * t9 + t13 * t66 + t37 * t28 - g(1) * t96 - g(2) * t98 + (-t102 * t176 + t171 * t4) * qJD(5), -(t156 * t10 + t159 * t12) * t102 - t176 * t75 + t180 * t171 - t4 * t108 - t38 * t175 + t58 * t8 - t13 * t67 + t37 * t27 - g(1) * t95 - g(2) * t97 + (-t102 * t177 + t171 * t3) * qJD(5); 0, 0, 0, -t157 * t163 * t160, t199 * t163, t194, t193, qJDD(2), t157 * t165 - t222, g(3) * t157 + t160 * t165, (t79 - t84) * t109 + (-t78 + t85) * t106 + (-t153 * t80 - t209 * t81) * pkin(2), t78 * t84 - t79 * t85 + (t209 * t33 - t222 + t153 * t34 + (-qJD(1) * t128 + t183) * t157) * pkin(2), -t136 * t214 - t23 * t109 + t139 * t59 + t84 * t94 + (t152 * t200 - t30) * t106 - t232 * t154, -t136 * t213 + t24 * t109 + t139 * t60 - t84 * t93 + (t154 * t200 + t31) * t106 + t232 * t152, -t224 + t30 * t93 - t31 * t94 - t183 * t146 + (qJD(4) * t94 - t106 * t23 - t136 * t59 + t7) * t154 + (qJD(4) * t93 - t106 * t24 + t136 * t60 - t6) * t152, t32 * t139 - t24 * t31 - t23 * t30 - t65 * t84 - g(3) * (t179 + t220) + (-t6 * t152 + t7 * t154) * t136 + t178 * qJD(4) + t183 * (pkin(2) * t157 + pkin(3) * t144 - qJ(4) * t146), t8 * t123 + t175 * t211, -t8 * t121 - t123 * t9 + t175 * t210 + t211 * t41, t212 - t228, t184 + t216, -t102 * t109, (-t117 * t159 - t118 * t156) * t75 + t127 * t9 + t13 * t121 - t3 * t109 - t45 * t41 + t210 * t37 + (t156 * t173 - t159 * t172) * t102 + t168 * t145, -(-t117 * t156 + t118 * t159) * t75 + t127 * t8 + t13 * t123 + t4 * t109 + t45 * t175 - t211 * t37 + (t156 * t172 + t159 * t173) * t102 - t168 * t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109 ^ 2 - t105, t79 * t106 + t78 * t109 + t169 - t230, -t105 * t152 + t109 * t94 + t213, -t105 * t154 - t109 * t93 - t214, -t152 * t59 - t154 * t60 + (t152 * t93 + t154 * t94) * t106, t106 * t178 - t65 * t109 - t181 - t230, 0, 0, 0, 0, 0, t184 - t216, t212 + t228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106 * t93 + t59, t106 * t94 + t60, -t93 ^ 2 - t94 ^ 2, t23 * t93 - t24 * t94 + t232, 0, 0, 0, 0, 0, t9 - t233, t8 - t234; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t175 * t41, t175 ^ 2 - t41 ^ 2, t8 + t234, -t9 - t233, t75, -g(1) * t97 + g(2) * t95 + t143 * t224 + t175 * t37 + t229 * t4 + t191, g(1) * t98 - g(2) * t96 + t145 * t224 + t229 * t3 + t37 * t41 - t180;];
tau_reg = t1;
