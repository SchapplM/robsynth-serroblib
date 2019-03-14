% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRPRP6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP6_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:01:14
% EndTime: 2019-03-09 17:01:21
% DurationCPUTime: 2.77s
% Computational Cost: add. (5788->323), mult. (15144->596), div. (0->0), fcn. (14833->10), ass. (0->161)
t130 = sin(qJ(5));
t133 = cos(qJ(5));
t127 = sin(pkin(11));
t134 = cos(qJ(3));
t131 = sin(qJ(3));
t194 = cos(pkin(11));
t160 = t194 * t131;
t103 = t127 * t134 + t160;
t159 = t194 * t134;
t181 = qJD(3) * t131;
t96 = qJD(3) * t159 - t127 * t181;
t148 = -qJ(6) * t96 - qJD(6) * t103;
t123 = pkin(3) * t181;
t95 = t103 * qJD(3);
t146 = pkin(4) * t95 - pkin(10) * t96 + t123;
t202 = -qJ(4) - pkin(9);
t161 = qJD(3) * t202;
t138 = -qJD(4) * t131 + t134 * t161;
t94 = qJD(4) * t134 + t131 * t161;
t63 = t127 * t138 + t194 * t94;
t163 = -t130 * t63 + t133 * t146;
t189 = t127 * t131;
t102 = -t159 + t189;
t122 = -pkin(3) * t134 - pkin(2);
t71 = pkin(4) * t102 - pkin(10) * t103 + t122;
t111 = t202 * t134;
t76 = -t194 * t111 + t202 * t189;
t73 = t133 * t76;
t15 = pkin(5) * t95 + t148 * t133 + (-t73 + (qJ(6) * t103 - t71) * t130) * qJD(5) + t163;
t177 = qJD(5) * t133;
t167 = t103 * t177;
t175 = t130 * t146 + t133 * t63 + t71 * t177;
t16 = -qJ(6) * t167 + (-qJD(5) * t76 + t148) * t130 + t175;
t192 = t103 * t133;
t47 = -t130 * t76 + t133 * t71;
t40 = pkin(5) * t102 - qJ(6) * t192 + t47;
t193 = t103 * t130;
t48 = t130 * t71 + t73;
t43 = -qJ(6) * t193 + t48;
t210 = -t130 * t16 - t133 * t15 + (t130 * t40 - t133 * t43) * qJD(5);
t129 = cos(pkin(6));
t128 = sin(pkin(6));
t135 = cos(qJ(2));
t182 = qJD(2) * t135;
t165 = t128 * t182;
t132 = sin(qJ(2));
t188 = t128 * t132;
t171 = t134 * t188;
t169 = qJD(3) * t171 + t129 * t181 + t131 * t165;
t97 = -t129 * t134 + t131 * t188;
t74 = -t97 * qJD(3) + t134 * t165;
t137 = -t127 * t169 + t194 * t74;
t183 = qJD(2) * t132;
t166 = t128 * t183;
t187 = t128 * t135;
t98 = t129 * t131 + t171;
t66 = -t127 * t97 + t194 * t98;
t55 = t130 * t66 + t133 * t187;
t33 = t55 * qJD(5) - t130 * t166 - t133 * t137;
t174 = pkin(8) * t187;
t206 = pkin(1) * t132;
t88 = t174 + (pkin(9) + t206) * t129;
t89 = (-pkin(2) * t135 - pkin(9) * t132 - pkin(1)) * t128;
t200 = t131 * t89 + t134 * t88;
t91 = (pkin(2) * t132 - pkin(9) * t135) * t128 * qJD(2);
t92 = -t129 * pkin(1) * t182 + pkin(8) * t166;
t42 = -t200 * qJD(3) + t131 * t92 + t134 * t91;
t26 = pkin(3) * t166 - qJ(4) * t74 - qJD(4) * t98 + t42;
t180 = qJD(3) * t134;
t41 = -t131 * t91 + t134 * t92 - t89 * t180 + t88 * t181;
t31 = -t169 * qJ(4) - t97 * qJD(4) - t41;
t12 = t127 * t26 + t194 * t31;
t10 = pkin(10) * t166 + t12;
t52 = t127 * t74 + t194 * t169;
t93 = (t129 * t206 + t174) * qJD(2);
t57 = t169 * pkin(3) + t93;
t136 = t52 * pkin(4) - t137 * pkin(10) + t57;
t162 = -t131 * t88 + t134 * t89;
t46 = -pkin(3) * t187 - t98 * qJ(4) + t162;
t53 = -qJ(4) * t97 + t200;
t30 = t127 * t46 + t194 * t53;
t25 = -pkin(10) * t187 + t30;
t65 = t127 * t98 + t194 * t97;
t87 = pkin(8) * t188 + (-pkin(1) * t135 - pkin(2)) * t129;
t70 = t97 * pkin(3) + t87;
t38 = t65 * pkin(4) - t66 * pkin(10) + t70;
t14 = t130 * t38 + t133 * t25;
t4 = -t14 * qJD(5) - t10 * t130 + t133 * t136;
t170 = t130 * t187;
t56 = t133 * t66 - t170;
t1 = pkin(5) * t52 + qJ(6) * t33 - qJD(6) * t56 + t4;
t178 = qJD(5) * t130;
t3 = -t133 * t10 - t130 * t136 - t38 * t177 + t25 * t178;
t34 = -qJD(5) * t170 + t130 * t137 - t133 * t166 + t66 * t177;
t2 = -qJ(6) * t34 - qJD(6) * t55 - t3;
t13 = -t130 * t25 + t133 * t38;
t6 = pkin(5) * t65 - qJ(6) * t56 + t13;
t7 = -qJ(6) * t55 + t14;
t209 = -t1 * t133 - t130 * t2 + (t130 * t6 - t133 * t7) * qJD(5);
t208 = 0.2e1 * t128;
t207 = 0.2e1 * qJD(5);
t205 = pkin(9) * t128;
t11 = -t127 * t31 + t194 * t26;
t9 = -pkin(4) * t166 - t11;
t204 = t130 * t9;
t203 = t133 * t9;
t201 = -t130 * t34 - t55 * t177;
t199 = t130 * t33;
t62 = t127 * t94 - t194 * t138;
t198 = t130 * t62;
t197 = t130 * t96;
t196 = t133 * t62;
t195 = t133 * t96;
t120 = pkin(3) * t127 + pkin(10);
t191 = t120 * t130;
t190 = t120 * t133;
t186 = t130 * t133;
t185 = qJ(6) + t120;
t125 = t130 ^ 2;
t126 = t133 ^ 2;
t184 = t125 - t126;
t179 = qJD(3) * t135;
t176 = -0.2e1 * pkin(2) * qJD(3);
t121 = -t194 * pkin(3) - pkin(4);
t173 = t121 * t207;
t172 = pkin(5) * t178;
t124 = t128 ^ 2;
t168 = t124 * t182;
t164 = t130 * t177;
t158 = -0.4e1 * t103 * t186;
t75 = -t111 * t127 - t202 * t160;
t157 = qJD(5) * t185;
t156 = t184 * qJD(5);
t155 = t132 * t168;
t153 = -t130 * t7 - t133 * t6;
t29 = -t127 * t53 + t194 * t46;
t151 = -t120 * t95 + t121 * t96;
t147 = t102 * t120 - t103 * t121;
t24 = pkin(4) * t187 - t29;
t145 = -t133 * t33 - t56 * t178;
t144 = t130 * t52 + t65 * t177;
t39 = t133 * t52 - t65 * t178;
t143 = t102 * t177 + t130 * t95;
t142 = t167 + t197;
t141 = t103 * t178 - t195;
t140 = t131 * t179 + t134 * t183;
t139 = t131 * t183 - t134 * t179;
t105 = -t133 * pkin(5) + t121;
t101 = t103 ^ 2;
t100 = t185 * t133;
t99 = t185 * t130;
t82 = -qJD(6) * t130 - t133 * t157;
t81 = qJD(6) * t133 - t130 * t157;
t69 = -t102 * t178 + t133 * t95;
t61 = pkin(5) * t193 + t75;
t44 = t142 * pkin(5) + t62;
t21 = -t48 * qJD(5) + t163;
t20 = t76 * t178 - t175;
t17 = t55 * pkin(5) + t24;
t5 = t34 * pkin(5) + t9;
t8 = [0, 0, 0, 0.2e1 * t155, 0.2e1 * (-t132 ^ 2 + t135 ^ 2) * t124 * qJD(2), 0.2e1 * t129 * t165, -0.2e1 * t129 * t166, 0, -0.2e1 * pkin(1) * t124 * t183 - 0.2e1 * t129 * t93, -0.2e1 * pkin(1) * t168 + 0.2e1 * t92 * t129, 0.2e1 * t98 * t74, -0.2e1 * t98 * t169 - 0.2e1 * t74 * t97 (-t135 * t74 + t98 * t183) * t208 (t169 * t135 - t97 * t183) * t208, -0.2e1 * t155, 0.2e1 * t93 * t97 + 0.2e1 * t87 * t169 + 0.2e1 * (-t42 * t135 + t162 * t183) * t128, 0.2e1 * t87 * t74 + 0.2e1 * t93 * t98 + 0.2e1 * (-t41 * t135 - t200 * t183) * t128, -0.2e1 * t11 * t66 - 0.2e1 * t12 * t65 - 0.2e1 * t137 * t29 - 0.2e1 * t30 * t52, 0.2e1 * t11 * t29 + 0.2e1 * t12 * t30 + 0.2e1 * t57 * t70, -0.2e1 * t56 * t33, 0.2e1 * t33 * t55 - 0.2e1 * t34 * t56, -0.2e1 * t33 * t65 + 0.2e1 * t52 * t56, -0.2e1 * t34 * t65 - 0.2e1 * t52 * t55, 0.2e1 * t65 * t52, 0.2e1 * t13 * t52 + 0.2e1 * t24 * t34 + 0.2e1 * t4 * t65 + 0.2e1 * t55 * t9, -0.2e1 * t14 * t52 - 0.2e1 * t24 * t33 + 0.2e1 * t3 * t65 + 0.2e1 * t56 * t9, -0.2e1 * t1 * t56 - 0.2e1 * t2 * t55 + 0.2e1 * t33 * t6 - 0.2e1 * t34 * t7, 0.2e1 * t1 * t6 + 0.2e1 * t17 * t5 + 0.2e1 * t2 * t7; 0, 0, 0, 0, 0, t165, -t166, 0, -t93, t92, t131 * t74 + t98 * t180, -t131 * t169 + t74 * t134 + (-t98 * t131 - t134 * t97) * qJD(3), t139 * t128, t140 * t128, 0, -pkin(2) * t169 - t93 * t134 - t139 * t205 + t87 * t181, -pkin(2) * t74 + t93 * t131 - t140 * t205 + t87 * t180, -t12 * t102 - t11 * t103 + t137 * t75 - t29 * t96 - t30 * t95 - t76 * t52 + t62 * t66 - t63 * t65, -t11 * t75 + t12 * t76 + t122 * t57 + t70 * t123 - t29 * t62 + t30 * t63, t103 * t145 + t56 * t195 (-t130 * t56 - t133 * t55) * t96 + (t199 - t133 * t34 + (t130 * t55 - t133 * t56) * qJD(5)) * t103, -t102 * t33 + t103 * t39 + t65 * t195 + t56 * t95, -t102 * t34 - t103 * t144 - t65 * t197 - t55 * t95, t102 * t52 + t65 * t95, t24 * t197 + t102 * t4 + t13 * t95 + t21 * t65 + t34 * t75 + t47 * t52 + t55 * t62 + (t24 * t177 + t204) * t103, t24 * t195 + t102 * t3 - t14 * t95 + t20 * t65 - t33 * t75 - t48 * t52 + t56 * t62 + (-t24 * t178 + t203) * t103, t209 * t103 - t15 * t56 + t153 * t96 - t16 * t55 + t33 * t40 - t34 * t43, t1 * t40 + t15 * t6 + t16 * t7 + t17 * t44 + t2 * t43 + t5 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t131 * t180, 0.2e1 * (-t131 ^ 2 + t134 ^ 2) * qJD(3), 0, 0, 0, t131 * t176, t134 * t176, -0.2e1 * t102 * t63 + 0.2e1 * t103 * t62 + 0.2e1 * t75 * t96 - 0.2e1 * t76 * t95, 0.2e1 * t122 * t123 + 0.2e1 * t62 * t75 + 0.2e1 * t63 * t76, 0.2e1 * t103 * t126 * t96 - 0.2e1 * t101 * t164, t184 * t101 * t207 + t96 * t158, -0.2e1 * t102 * t141 + 0.2e1 * t95 * t192, -0.2e1 * t102 * t142 - 0.2e1 * t95 * t193, 0.2e1 * t102 * t95, 0.2e1 * t75 * t197 + 0.2e1 * t102 * t21 + 0.2e1 * t47 * t95 + 0.2e1 * (t75 * t177 + t198) * t103, 0.2e1 * t75 * t195 + 0.2e1 * t102 * t20 - 0.2e1 * t48 * t95 + 0.2e1 * (-t75 * t178 + t196) * t103, 0.2e1 * (-t130 * t43 - t133 * t40) * t96 + 0.2e1 * t210 * t103, 0.2e1 * t15 * t40 + 0.2e1 * t16 * t43 + 0.2e1 * t44 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, -t169, t166, t42, t41 (-t127 * t52 - t194 * t137) * pkin(3) (t194 * t11 + t12 * t127) * pkin(3), t56 * t177 - t199, t145 + t201, t144, t39, 0, -t52 * t191 + t121 * t34 - t203 + (t130 * t24 - t65 * t190) * qJD(5), -t52 * t190 - t121 * t33 + t204 + (t133 * t24 + t65 * t191) * qJD(5), qJD(5) * t153 - t1 * t130 - t100 * t34 + t133 * t2 - t33 * t99 - t55 * t81 - t56 * t82, -t1 * t99 + t100 * t2 + t105 * t5 + t17 * t172 + t6 * t82 + t7 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t180, -t181, 0, -pkin(9) * t180, pkin(9) * t181 (-t127 * t95 - t194 * t96) * pkin(3) (t127 * t63 - t194 * t62) * pkin(3), -t103 * t156 + t96 * t186, qJD(5) * t158 - t184 * t96, t143, t69, 0, -t196 + t151 * t130 + (t130 * t75 - t133 * t147) * qJD(5), t198 + t151 * t133 + (t130 * t147 + t133 * t75) * qJD(5) (-t103 * t82 + t96 * t99 + t16 + (-t100 * t103 - t40) * qJD(5)) * t133 + (-t100 * t96 - t103 * t81 - t15 + (-t103 * t99 - t43) * qJD(5)) * t130, t100 * t16 + t105 * t44 - t15 * t99 + t172 * t61 + t40 * t82 + t43 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t164, -0.2e1 * t156, 0, 0, 0, t130 * t173, t133 * t173, -0.2e1 * t130 * t82 + 0.2e1 * t133 * t81 + 0.2e1 * (-t100 * t130 + t133 * t99) * qJD(5), 0.2e1 * t100 * t81 + 0.2e1 * t105 * t172 - 0.2e1 * t82 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, 0, 0, 0, 0, t39, -t144, -t145 + t201, -t209; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t123, 0, 0, 0, 0, 0, t69, -t143 (-t125 - t126) * t96, -t210; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130 * t81 + t133 * t82 + (t100 * t133 + t130 * t99) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, -t34, t52, t4, t3, pkin(5) * t33, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141, -t142, t95, t21, t20, t141 * pkin(5), t15 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t177, -t178, 0, -t120 * t177, t120 * t178, -pkin(5) * t177, t82 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t178, -t177, 0, -t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t8;