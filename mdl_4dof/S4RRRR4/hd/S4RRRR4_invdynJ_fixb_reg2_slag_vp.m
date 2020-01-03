% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRRR4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR4_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:26:18
% EndTime: 2019-12-31 17:26:22
% DurationCPUTime: 2.37s
% Computational Cost: add. (3444->340), mult. (8054->468), div. (0->0), fcn. (5463->10), ass. (0->180)
t117 = qJ(2) + qJ(3);
t113 = sin(t117);
t124 = cos(qJ(1));
t196 = t113 * t124;
t121 = sin(qJ(1));
t197 = t113 * t121;
t225 = g(1) * t196 + g(2) * t197;
t178 = qJDD(2) + qJDD(3);
t119 = sin(qJ(3));
t219 = cos(qJ(3));
t125 = -pkin(6) - pkin(5);
t123 = cos(qJ(2));
t182 = qJD(1) * qJD(2);
t164 = t123 * t182;
t120 = sin(qJ(2));
t181 = t120 * qJDD(1);
t57 = qJDD(2) * pkin(2) - t125 * (-t164 - t181);
t165 = t120 * t182;
t180 = t123 * qJDD(1);
t58 = t125 * (-t165 + t180);
t163 = -t119 * t58 - t219 * t57;
t88 = t125 * t123;
t83 = qJD(1) * t88;
t174 = t219 * t83;
t209 = qJD(2) * pkin(2);
t87 = t125 * t120;
t81 = qJD(1) * t87;
t77 = t81 + t209;
t49 = t119 * t77 - t174;
t19 = -t49 * qJD(3) - t163;
t17 = -t178 * pkin(3) - t19;
t114 = cos(t117);
t218 = g(3) * t114;
t224 = t17 + t218;
t118 = sin(qJ(4));
t122 = cos(qJ(4));
t168 = t219 * t123;
t186 = qJD(1) * t120;
t72 = -qJD(1) * t168 + t119 * t186;
t191 = t119 * t123;
t79 = t219 * t120 + t191;
t74 = t79 * qJD(1);
t112 = pkin(2) * t123 + pkin(1);
t86 = t112 * qJD(1);
t42 = pkin(3) * t72 - pkin(7) * t74 - t86;
t179 = qJD(2) + qJD(3);
t45 = t179 * pkin(7) + t49;
t22 = -t118 * t45 + t122 * t42;
t23 = t118 * t42 + t122 * t45;
t147 = t118 * t23 + t122 * t22;
t151 = g(1) * t121 - g(2) * t124;
t223 = t151 * t113;
t153 = t114 * pkin(3) + t113 * pkin(7);
t222 = -t118 * t22 + t122 * t23;
t221 = t79 * qJDD(1);
t139 = -t119 * t120 + t168;
t220 = t179 * qJD(1);
t37 = -t139 * t220 - t221;
t61 = t118 * t179 + t122 * t74;
t21 = t61 * qJD(4) - t118 * t37 - t122 * t178;
t108 = g(3) * t113;
t217 = g(3) * t123;
t148 = -qJDD(1) * t168 + t119 * t181;
t55 = t179 * t79;
t38 = t55 * qJD(1) + t148;
t70 = pkin(2) * t165 - t112 * qJDD(1);
t12 = pkin(3) * t38 + pkin(7) * t37 + t70;
t166 = t219 * qJD(3);
t185 = qJD(3) * t119;
t18 = t119 * t57 + t77 * t166 + t83 * t185 - t219 * t58;
t16 = t178 * pkin(7) + t18;
t3 = qJD(4) * t22 + t118 * t12 + t122 * t16;
t2 = t3 * t122;
t204 = t119 * t83;
t48 = t219 * t77 + t204;
t44 = -t179 * pkin(3) - t48;
t216 = t44 * t72;
t159 = t122 * t179;
t59 = t118 * t74 - t159;
t69 = qJD(4) + t72;
t215 = t59 * t69;
t214 = t61 * t59;
t213 = t61 * t69;
t212 = t69 * t74;
t211 = t74 * t72;
t52 = t119 * t81 - t174;
t210 = t49 - t52;
t184 = qJD(4) * t118;
t20 = -qJD(4) * t159 - t118 * t178 + t122 * t37 + t74 * t184;
t208 = t118 * t20;
t205 = t118 * t59;
t203 = t122 * t21;
t200 = t122 * t61;
t199 = pkin(5) * qJDD(1);
t198 = qJD(4) * t69;
t195 = t114 * t121;
t194 = t114 * t124;
t193 = t118 * t121;
t192 = t118 * t124;
t190 = t121 * t122;
t189 = t122 * t124;
t115 = t120 ^ 2;
t116 = t123 ^ 2;
t188 = t115 - t116;
t187 = t115 + t116;
t183 = qJD(4) * t122;
t177 = pkin(7) * t198;
t176 = g(1) * t194 + g(2) * t195 + t108;
t175 = t120 * t209;
t110 = t119 * pkin(2) + pkin(7);
t173 = t110 * t198;
t172 = t79 * t184;
t171 = t79 * t183;
t127 = qJD(1) ^ 2;
t170 = t120 * t127 * t123;
t40 = t44 * t184;
t41 = t44 * t183;
t169 = qJD(2) * t125;
t161 = t122 * t69;
t160 = t119 * t179;
t158 = pkin(2) * t166;
t157 = qJD(2) * t168;
t156 = t120 * t164;
t155 = pkin(2) * t185 - t52;
t46 = pkin(3) * t74 + pkin(7) * t72;
t36 = qJDD(4) + t38;
t154 = -pkin(7) * t36 + t216;
t152 = g(1) * t124 + g(2) * t121;
t54 = t120 * t160 - t123 * t166 - t157;
t150 = t17 * t79 - t44 * t54;
t149 = t36 * t79 - t54 * t69;
t47 = -pkin(3) * t139 - pkin(7) * t79 - t112;
t63 = t119 * t87 - t219 * t88;
t31 = -t118 * t63 + t122 * t47;
t32 = t118 * t47 + t122 * t63;
t146 = t224 * t118 + t23 * t74 + t41;
t145 = t225 * t122 - t22 * t74 + t40;
t144 = -t147 * t72 - t176 + t2;
t143 = -qJD(4) * t42 + t108 - t16;
t142 = t119 * t88 + t219 * t87;
t141 = t152 * t113;
t140 = -0.2e1 * pkin(1) * t182 - pkin(5) * qJDD(2);
t136 = -t86 * t72 + t176 - t18;
t11 = t122 * t12;
t4 = -qJD(4) * t23 - t118 * t16 + t11;
t135 = -t147 * qJD(4) - t4 * t118;
t126 = qJD(2) ^ 2;
t134 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t126 + t151;
t133 = pkin(1) * t127 + t152 - t199;
t132 = t86 * t74 - t163 - t218 + t225;
t131 = -t110 * t36 - t69 * t158 + t216;
t130 = t135 + t2;
t129 = -g(1) * (-pkin(3) * t196 + pkin(7) * t194) - g(2) * (-pkin(3) * t197 + pkin(7) * t195) - g(3) * t153;
t111 = -t219 * pkin(2) - pkin(3);
t89 = t124 * t112;
t82 = t120 * t169;
t68 = t114 * t189 + t193;
t67 = -t114 * t192 + t190;
t66 = -t114 * t190 + t192;
t65 = t114 * t193 + t189;
t53 = t219 * t81 + t204;
t43 = pkin(2) * t186 + t46;
t39 = -t72 ^ 2 + t74 ^ 2;
t34 = t63 * qJD(3) + t119 * t82 - t125 * t157;
t33 = t142 * qJD(3) + t169 * t191 + t219 * t82;
t30 = pkin(3) * t55 + pkin(7) * t54 + t175;
t29 = t74 * t179 - t79 * t220 - t148;
t28 = t221 + (qJD(1) * t139 + t72) * t179;
t27 = t118 * t46 + t122 * t48;
t26 = -t118 * t48 + t122 * t46;
t25 = t118 * t43 + t122 * t53;
t24 = -t118 * t53 + t122 * t43;
t10 = t118 * t36 + t69 * t161 - t61 * t74;
t9 = -t69 ^ 2 * t118 + t122 * t36 + t59 * t74;
t8 = t69 * t205 - t203;
t7 = t61 * t161 - t208;
t6 = -t32 * qJD(4) - t118 * t33 + t122 * t30;
t5 = t31 * qJD(4) + t118 * t30 + t122 * t33;
t1 = (-t20 - t215) * t122 + (-t21 - t213) * t118;
t13 = [0, 0, 0, 0, 0, qJDD(1), t151, t152, 0, 0, qJDD(1) * t115 + 0.2e1 * t156, 0.2e1 * t120 * t180 - 0.2e1 * t188 * t182, qJDD(2) * t120 + t123 * t126, qJDD(1) * t116 - 0.2e1 * t156, qJDD(2) * t123 - t120 * t126, 0, t140 * t120 + t134 * t123, -t134 * t120 + t140 * t123, 0.2e1 * t187 * t199 - t152, -g(1) * (-pkin(1) * t121 + pkin(5) * t124) - g(2) * (pkin(1) * t124 + pkin(5) * t121) + (t187 * pkin(5) ^ 2 + pkin(1) ^ 2) * qJDD(1), -t37 * t79 - t54 * t74, -t139 * t37 - t38 * t79 + t54 * t72 - t55 * t74, t178 * t79 - t179 * t54, -t139 * t38 + t55 * t72, t139 * t178 - t179 * t55, 0, -t112 * t38 + t114 * t151 - t139 * t70 + t142 * t178 + t175 * t72 - t179 * t34 - t86 * t55, t112 * t37 + t175 * t74 - t178 * t63 - t179 * t33 + t86 * t54 + t70 * t79 - t223, t139 * t18 + t142 * t37 - t19 * t79 - t33 * t72 + t34 * t74 - t38 * t63 + t48 * t54 - t49 * t55 - t152, t18 * t63 + t49 * t33 + t19 * t142 - t48 * t34 - t70 * t112 - t86 * t175 - g(1) * (-t112 * t121 - t124 * t125) - g(2) * (-t121 * t125 + t89), -t61 * t172 + (-t20 * t79 - t54 * t61) * t122, (t118 * t61 + t122 * t59) * t54 + (t208 - t203 + (-t200 + t205) * qJD(4)) * t79, t122 * t149 + t139 * t20 - t172 * t69 + t55 * t61, t59 * t171 + (t21 * t79 - t54 * t59) * t118, -t118 * t149 + t139 * t21 - t171 * t69 - t55 * t59, -t139 * t36 + t55 * t69, -g(1) * t66 - g(2) * t68 + t118 * t150 - t139 * t4 - t142 * t21 + t22 * t55 + t31 * t36 + t34 * t59 + t41 * t79 + t6 * t69, -g(1) * t65 - g(2) * t67 + t122 * t150 + t139 * t3 + t142 * t20 - t23 * t55 - t32 * t36 + t34 * t61 - t40 * t79 - t5 * t69, t20 * t31 - t21 * t32 - t5 * t59 - t6 * t61 + t147 * t54 + t223 + (-t222 * qJD(4) - t118 * t3 - t122 * t4) * t79, -g(2) * t89 - t17 * t142 + t22 * t6 + t23 * t5 + t3 * t32 + t4 * t31 + t44 * t34 + (g(1) * t125 - g(2) * t153) * t124 + (-g(1) * (-t112 - t153) + g(2) * t125) * t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t170, t188 * t127, t181, t170, t180, qJDD(2), t133 * t120 - t217, g(3) * t120 + t123 * t133, 0, 0, t211, t39, t28, -t211, t29, t178, t52 * qJD(2) - t210 * qJD(3) + (-qJD(3) * t160 + t219 * t178 - t72 * t186) * pkin(2) + t132, t53 * t179 + (-t119 * t178 - t166 * t179 - t74 * t186) * pkin(2) + t136, t210 * t74 + (-t48 + t53) * t72 + (t219 * t37 - t119 * t38 + (t119 * t74 - t219 * t72) * qJD(3)) * pkin(2), t48 * t52 - t49 * t53 + (t219 * t19 - t217 + t119 * t18 + (-t119 * t48 + t219 * t49) * qJD(3) + (qJD(1) * t86 + t152) * t120) * pkin(2), t7, t1, t10, t8, t9, -t212, t111 * t21 - t24 * t69 + t155 * t59 + (-t224 - t173) * t122 + t131 * t118 + t145, -t111 * t20 + t25 * t69 + t155 * t61 + t131 * t122 + (-t141 + t173) * t118 + t146, t24 * t61 + t25 * t59 + (-t59 * t158 - t110 * t21 + (t110 * t61 - t22) * qJD(4)) * t122 + (t61 * t158 - t110 * t20 - t4 + (t110 * t59 - t23) * qJD(4)) * t118 + t144, t17 * t111 - t23 * t25 - t22 * t24 - t44 * t52 + (-t217 + t152 * t120 + (t119 * t44 + t222 * t219) * qJD(3)) * pkin(2) + t130 * t110 + t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t211, t39, t28, -t211, t29, t178, t49 * qJD(2) + t132, t179 * t48 + t136, 0, 0, t7, t1, t10, t8, t9, -t212, -pkin(3) * t21 - t26 * t69 - t49 * t59 + t154 * t118 + (-t224 - t177) * t122 + t145, pkin(3) * t20 + t27 * t69 - t49 * t61 + t154 * t122 + (-t141 + t177) * t118 + t146, t26 * t61 + t27 * t59 + (-t208 - t203 + (t200 + t205) * qJD(4)) * pkin(7) + t135 + t144, -t17 * pkin(3) + pkin(7) * t130 - t22 * t26 - t23 * t27 - t44 * t49 + t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t214, -t59 ^ 2 + t61 ^ 2, -t20 + t215, -t214, -t21 + t213, t36, -g(1) * t67 + g(2) * t65 + t118 * t143 - t183 * t45 + t23 * t69 - t44 * t61 + t11, g(1) * t68 - g(2) * t66 + t22 * t69 + t44 * t59 + (qJD(4) * t45 - t12) * t118 + t143 * t122, 0, 0;];
tau_reg = t13;
