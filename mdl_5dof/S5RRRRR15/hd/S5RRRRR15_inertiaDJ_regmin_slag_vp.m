% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x31]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRR15_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR15_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:26:25
% EndTime: 2024-09-27 22:26:28
% DurationCPUTime: 2.41s
% Computational Cost: add. (5935->225), mult. (17972->420), div. (0->0), fcn. (17708->12), ass. (0->155)
t104 = sin(qJ(3));
t105 = sin(qJ(2));
t108 = cos(qJ(3));
t109 = cos(qJ(2));
t99 = sin(pkin(5));
t70 = (-t104 * t105 + t108 * t109) * t99;
t100 = cos(pkin(6));
t103 = sin(qJ(4));
t107 = cos(qJ(4));
t101 = cos(pkin(5));
t154 = t101 * t105 * pkin(1);
t170 = t109 * t99;
t189 = pkin(8) + pkin(9);
t65 = t170 * t189 + t154;
t171 = t108 * t65;
t151 = t189 * t99;
t128 = t105 * t151;
t61 = -t128 + (pkin(1) * t109 + pkin(2)) * t101;
t122 = -t104 * t61 - t171;
t40 = pkin(10) * t70 - t122;
t172 = t107 * t40;
t174 = t104 * t65;
t186 = t101 * pkin(3);
t121 = t104 * t109 + t105 * t108;
t71 = t121 * t99;
t39 = -pkin(10) * t71 + t108 * t61 - t174 + t186;
t123 = -t103 * t39 - t172;
t115 = t121 * qJD(3);
t111 = (-qJD(2) * t121 - t115) * t99;
t160 = qJD(3) * t108;
t162 = qJD(2) * t109;
t140 = t101 * t162;
t87 = pkin(1) * t140;
t62 = -qJD(2) * t128 + t87;
t63 = (-t109 * t151 - t154) * qJD(2);
t155 = -t104 * t63 - t108 * t62 - t160 * t61;
t161 = qJD(3) * t104;
t34 = t161 * t65 + t155;
t24 = pkin(10) * t111 - t34;
t131 = -t104 * t62 + t108 * t63;
t35 = qJD(3) * t122 + t131;
t52 = (qJD(3) + qJD(2)) * t70;
t25 = -t52 * pkin(10) + t35;
t133 = -t103 * t24 + t107 * t25;
t13 = qJD(4) * t123 + t133;
t187 = pkin(11) * t100;
t48 = t103 * t71 - t107 * t70;
t27 = -qJD(4) * t48 + t103 * t111 + t107 * t52;
t49 = t103 * t70 + t107 * t71;
t28 = qJD(4) * t49 + t103 * t52 - t107 * t111;
t185 = t105 * pkin(2);
t42 = (pkin(3) * t115 + (pkin(3) * t121 + t185) * qJD(2)) * t99;
t98 = sin(pkin(6));
t95 = t98 * pkin(11);
t193 = t98 * (t28 * pkin(4) - t27 * t95 + t42) + t100 * (-t187 * t27 + t13);
t159 = qJD(4) * t103;
t132 = -t103 * t25 + t159 * t40;
t12 = -(qJD(4) * t39 + t24) * t107 + t132;
t191 = -0.2e1 * t101;
t190 = 0.2e1 * t101;
t102 = sin(qJ(5));
t157 = qJD(5) * t102;
t142 = t98 * t157;
t22 = pkin(4) * t101 - t103 * t40 + t107 * t39 - t187 * t49;
t78 = (-pkin(2) * t109 - pkin(1)) * t99;
t56 = -t70 * pkin(3) + t78;
t33 = t48 * pkin(4) - t49 * t95 + t56;
t18 = t100 * t33 - t22 * t98;
t106 = cos(qJ(5));
t179 = t100 * t28;
t114 = pkin(11) * t179 + t12;
t139 = t100 * t157;
t156 = qJD(5) * t106;
t124 = -t100 * t48 + t101 * t98;
t21 = pkin(11) * t124 - t123;
t3 = t102 * t114 + t106 * t193 - t22 * t139 - t33 * t142 - t21 * t156;
t188 = t100 * t3 + t142 * t18;
t184 = t108 * pkin(2);
t6 = pkin(4) * t179 + t100 * t42 - t13 * t98;
t183 = t6 * t106;
t182 = t98 * t28;
t96 = t98 ^ 2;
t173 = t106 * t96;
t168 = t100 * t102;
t165 = t104 * t107;
t94 = pkin(3) + t184;
t68 = pkin(2) * t165 + t103 * t94 + t95;
t153 = t103 * t104 * pkin(2);
t76 = t107 * t94 + pkin(4) - t153;
t119 = t106 * t68 + t168 * t76;
t167 = t100 * t106;
t149 = pkin(2) * t160;
t150 = pkin(2) * t161;
t180 = qJD(4) * t153 + t103 * t150;
t54 = (-qJD(4) * t94 - t149) * t107 + t180;
t158 = qJD(4) * t107;
t113 = (-t104 * t158 + (-t103 * t108 - t165) * qJD(3)) * pkin(2);
t55 = -t159 * t94 + t113;
t32 = -qJD(5) * t119 + t102 * t54 + t167 * t55;
t181 = t100 * t32 + t173 * t55;
t178 = t102 * t49;
t177 = t102 * t55;
t176 = t102 * t96;
t175 = t102 * t98;
t169 = t98 * t106;
t163 = qJD(2) * t105;
t152 = t98 * t179;
t148 = pkin(3) * t159;
t147 = pkin(3) * t158;
t146 = t28 * t167;
t97 = t99 ^ 2;
t145 = t97 * t162;
t144 = t96 * t157;
t143 = t96 * t156;
t91 = t98 * t156;
t141 = t99 * t163;
t138 = t100 * t156;
t137 = -t39 - t186;
t136 = qJD(5) * (-pkin(4) - t76);
t93 = pkin(3) * t107 + pkin(4);
t135 = qJD(5) * (-pkin(4) - t93);
t134 = (-pkin(3) - t94) * qJD(4);
t130 = qJD(5) * (-t76 - t93);
t127 = t106 * t148;
t2 = -t102 * t193 + t106 * t114 - t22 * t138 + t21 * t157 - t33 * t91;
t126 = t2 * t100 + t175 * t6 + t18 * t91;
t125 = t100 * t22 + t33 * t98;
t83 = pkin(3) * t103 + t95;
t118 = t106 * t83 + t168 * t93;
t117 = pkin(4) * t168 + pkin(11) * t169;
t116 = t78 * t121;
t37 = t102 * t124 + t106 * t49;
t82 = 0.2e1 * t102 * t143;
t81 = t148 * t176;
t80 = 0.2e1 * t98 * t138;
t79 = -0.2e1 * t98 * t139;
t75 = (-pkin(8) * t170 - t154) * qJD(2);
t74 = pkin(8) * t141 - t87;
t73 = t117 * qJD(5);
t72 = -pkin(4) * t138 + pkin(11) * t142;
t69 = 0.2e1 * (-t102 ^ 2 + t106 ^ 2) * t96 * qJD(5);
t67 = t73 * t100;
t46 = -t118 * qJD(5) + (-t102 * t107 - t103 * t167) * qJD(4) * pkin(3);
t45 = -t93 * t138 - t106 * t147 + (qJD(5) * t83 + t100 * t148) * t102;
t43 = t46 * t100;
t41 = -t100 * t101 - t48 * t98;
t36 = -t101 * t169 + t167 * t48 + t178;
t31 = t106 * t54 - t138 * t76 + t157 * t68 - t168 * t55;
t15 = qJD(5) * t37 + t102 * t27 + t146;
t14 = -t28 * t168 + t106 * t27 + (t106 * t124 - t178) * qJD(5);
t11 = (t102 * t14 + t156 * t37) * t98;
t10 = t100 * t14 + t176 * t28 - t41 * t91;
t9 = -t100 * t15 + t142 * t41 + t173 * t28;
t4 = (-t102 * t15 + t106 * t14 + (-t102 * t37 - t106 * t36) * qJD(5)) * t98;
t1 = [0, 0, 0, 0.2e1 * t105 * t145, 0.2e1 * (-t105 ^ 2 + t109 ^ 2) * t97 * qJD(2), 0.2e1 * t99 * t140, t141 * t191, 0, -0.2e1 * pkin(1) * t163 * t97 + 0.2e1 * t101 * t75, -0.2e1 * pkin(1) * t145 + 0.2e1 * t101 * t74, 0.2e1 * t71 * t52, 0.2e1 * t111 * t71 + 0.2e1 * t52 * t70, t52 * t190, t111 * t190, 0, 0.2e1 * t35 * t101 + 0.2e1 * (qJD(3) * t116 + (-t185 * t70 + t116) * qJD(2)) * t99, 0.2e1 * pkin(2) * t141 * t71 + 0.2e1 * t101 * t34 + 0.2e1 * t52 * t78, 0.2e1 * t49 * t27, -0.2e1 * t27 * t48 - 0.2e1 * t28 * t49, t27 * t190, t28 * t191, 0, 0.2e1 * t101 * t13 + 0.2e1 * t28 * t56 + 0.2e1 * t42 * t48, 0.2e1 * t101 * t12 + 0.2e1 * t27 * t56 + 0.2e1 * t42 * t49, 0.2e1 * t37 * t14, -0.2e1 * t14 * t36 - 0.2e1 * t15 * t37, -0.2e1 * t14 * t41 + 0.2e1 * t182 * t37, 0.2e1 * t15 * t41 - 0.2e1 * t182 * t36, -0.2e1 * t41 * t182, -0.2e1 * t3 * t41 + 0.2e1 * (-t102 * t21 + t106 * t125) * t182 + 0.2e1 * t6 * t36 + 0.2e1 * t18 * t15, -0.2e1 * t2 * t41 - 0.2e1 * (t102 * t125 + t106 * t21) * t182 + 0.2e1 * t6 * t37 + 0.2e1 * t18 * t14; 0, 0, 0, 0, 0, t99 * t162, -t141, 0, t75, t74, 0, 0, t52, t111, 0, (-t171 + (-t101 * pkin(2) - t61) * t104) * qJD(3) + t131, (-t101 * t184 + t174) * qJD(3) + t155, 0, 0, t27, -t28, 0, t55 * t101 + t13, t54 * t101 + t12, t11, t4, t10, t9, t152, -t32 * t41 + ((-t102 * t68 + t167 * t76) * t28 - t55 * t36 - t76 * t15 - t183) * t98 + t188, -t31 * t41 + (-t119 * t28 - t76 * t14 - t55 * t37) * t98 + t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t150, -0.2e1 * t149, 0, 0, 0, 0, 0, 0.2e1 * t55, 0.2e1 * t54, t82, t69, t80, t79, 0, -0.2e1 * t144 * t76 + 0.2e1 * t181, 0.2e1 * t31 * t100 + 0.2e1 * (-t156 * t76 - t177) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t111, 0, t35, t34, 0, 0, t27, -t28, 0, (t103 * t137 - t172) * qJD(4) + t133, (qJD(4) * t137 - t24) * t107 + t132, t11, t4, t10, t9, t152, -t46 * t41 + ((-t102 * t83 + t167 * t93) * t28 + t36 * t148 - t93 * t15 - t183) * t98 + t188, -t45 * t41 + (-t118 * t28 - t93 * t14 + t148 * t37) * t98 + t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t150, -t149, 0, 0, 0, 0, 0, t103 * t134 + t113, (t134 - t149) * t107 + t180, t82, t69, t80, t79, 0, t43 + (t102 * t130 - t127) * t96 + t181, t81 + (t31 + t45) * t100 + (t106 * t130 - t177) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t148, -0.2e1 * t147, t82, t69, t80, t79, 0, 0.2e1 * t43 + 0.2e1 * (-t157 * t93 - t127) * t96, 0.2e1 * t45 * t100 - 0.2e1 * t143 * t93 + 0.2e1 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t28, 0, t13, t12, t11, t4, t10, t9, t152, t73 * t41 + (-pkin(11) * t28 * t175 - t183 + (-t15 + t146) * pkin(4)) * t98 + t188, -t98 * pkin(4) * t14 - t117 * t182 - t72 * t41 + t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t54, t82, t69, t80, t79, 0, t136 * t176 + t181 - t67, (t31 + t72) * t100 + (t106 * t136 - t177) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148, -t147, t82, t69, t80, t79, 0, t43 - t67 + (t102 * t135 - t127) * t96, t81 + (t45 + t72) * t100 + t135 * t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, t69, t80, t79, 0, -0.2e1 * pkin(4) * t144 - 0.2e1 * t67, -0.2e1 * pkin(4) * t143 + 0.2e1 * t72 * t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t15, t182, t3, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, -t142, 0, t32, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, -t142, 0, t46, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, -t142, 0, -t73, t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
