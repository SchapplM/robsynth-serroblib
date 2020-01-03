% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRR8_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:06:06
% EndTime: 2019-12-31 19:06:10
% DurationCPUTime: 1.42s
% Computational Cost: add. (2793->196), mult. (4455->281), div. (0->0), fcn. (2472->6), ass. (0->146)
t82 = qJD(4) + qJD(5);
t139 = qJD(1) - qJD(3);
t182 = t139 ^ 2;
t88 = sin(qJ(3));
t93 = qJD(4) ^ 2;
t183 = (t93 + t182) * t88;
t142 = t88 * qJD(2);
t91 = cos(qJ(3));
t144 = qJD(3) * t91;
t92 = -pkin(1) - pkin(2);
t71 = t92 * qJD(1) + qJD(2);
t159 = t71 * t88;
t36 = (qJ(2) * t144 + t142) * qJD(1) + qJD(3) * t159;
t141 = qJ(2) * qJD(1);
t52 = t91 * t141 + t159;
t179 = -t52 * t139 - t36;
t86 = sin(qJ(5));
t87 = sin(qJ(4));
t89 = cos(qJ(5));
t90 = cos(qJ(4));
t59 = t86 * t90 + t87 * t89;
t181 = t82 * t59;
t26 = t181 * t139;
t170 = pkin(3) * t139;
t75 = t88 * t141;
t51 = t91 * t71 - t75;
t40 = -t51 + t170;
t180 = t139 * t40;
t174 = t139 * t91;
t156 = t89 * t90;
t177 = t82 * t156;
t176 = qJD(4) * t139;
t140 = qJD(1) * qJD(2);
t35 = -qJD(3) * t75 + t91 * t140 + t71 * t144;
t175 = t139 * t51 + t35;
t157 = t86 * t87;
t58 = -t156 + t157;
t54 = t58 * t88;
t147 = t91 * qJ(2) + t88 * t92;
t41 = -pkin(7) * t139 + t52;
t128 = -pkin(8) * t139 + t41;
t108 = qJD(4) * t128;
t12 = -t87 * t108 + t90 * t35;
t30 = t128 * t87;
t28 = qJD(4) * pkin(4) - t30;
t173 = (qJD(5) * t28 + t12) * t89;
t84 = t87 ^ 2;
t85 = t90 ^ 2;
t145 = t84 + t85;
t124 = t145 * t35;
t50 = t147 * qJD(3) + t142;
t172 = t139 * t50 + t36;
t120 = -t88 * qJ(2) + t91 * t92;
t171 = -pkin(8) - pkin(7);
t169 = pkin(4) * t90;
t61 = -pkin(7) + t147;
t168 = pkin(8) - t61;
t31 = t128 * t90;
t167 = t31 * t86;
t166 = t31 * t89;
t77 = -pkin(3) - t169;
t32 = -t139 * t77 - t51;
t48 = t59 * t139;
t165 = t32 * t48;
t113 = t82 * t157;
t33 = t113 - t177;
t164 = t33 * t82;
t163 = t181 * t82;
t137 = t139 * t157;
t46 = t139 * t156 - t137;
t162 = t48 * t46;
t158 = t139 * t87;
t154 = t93 * t87;
t153 = t93 * t90;
t152 = t58 * t174 - t181 * t88;
t151 = t59 * t174 + t82 * t54;
t69 = t171 * t87;
t70 = t171 * t90;
t37 = t69 * t89 + t70 * t86;
t126 = qJD(4) * t171;
t62 = t87 * t126;
t63 = t90 * t126;
t150 = t37 * qJD(5) + t58 * t51 + t89 * t62 + t86 * t63;
t38 = t69 * t86 - t70 * t89;
t149 = -t38 * qJD(5) + t59 * t51 - t86 * t62 + t89 * t63;
t148 = t177 * t139;
t146 = t84 - t85;
t143 = qJD(4) * t87;
t138 = pkin(4) * t158;
t136 = t87 * t182 * t90;
t135 = pkin(4) * t143;
t133 = t139 * t143;
t130 = 0.2e1 * t140;
t129 = -pkin(4) * t82 - t28;
t13 = -t90 * t108 - t87 * t35;
t127 = -t86 * t12 + t89 * t13;
t125 = t145 * t51;
t49 = t91 * qJD(2) + t120 * qJD(3);
t123 = t145 * t49;
t122 = -qJD(5) * t167 + t86 * t13;
t121 = -t35 + t180;
t119 = qJD(4) * t168;
t117 = t139 * t88;
t115 = t90 * t133;
t60 = pkin(3) - t120;
t114 = -t52 + t135;
t25 = -t113 * t139 + t148;
t112 = t25 * t59 - t33 * t48;
t111 = t181 * t46 - t26 * t58;
t6 = t28 * t86 + t166;
t29 = -pkin(4) * t133 + t36;
t110 = t29 * t59 - t32 * t33;
t109 = t181 * t32 + t29 * t58;
t44 = t168 * t87;
t45 = t168 * t90;
t21 = t44 * t89 + t45 * t86;
t22 = t44 * t86 - t45 * t89;
t107 = pkin(7) * t93 - t179;
t106 = t32 * t46 - t122;
t105 = t61 * t93 - t172;
t104 = qJD(4) * (t40 + t51 + t170);
t103 = qJD(4) * (-t139 * t60 - t40 - t49);
t102 = 0.2e1 * t91 * t176;
t98 = t139 * t145;
t1 = t122 + t173;
t2 = -t6 * qJD(5) + t127;
t5 = t28 * t89 - t167;
t97 = -t1 * t58 - t181 * t6 - t2 * t59 + t33 * t5;
t96 = t181 * t48 + t25 * t58 + t26 * t59 + t33 * t46;
t94 = qJD(1) ^ 2;
t65 = 0.2e1 * t115;
t64 = -0.2e1 * t115;
t56 = t60 + t169;
t55 = t146 * t176;
t53 = t59 * t88;
t39 = t50 - t135;
t24 = t90 * t119 - t87 * t49;
t23 = t87 * t119 + t90 * t49;
t14 = -t46 ^ 2 + t48 ^ 2;
t11 = -t48 * t82 + t26;
t10 = -t148 + (t137 + t46) * t82;
t8 = -t30 * t89 - t167;
t7 = t30 * t86 - t166;
t4 = -t22 * qJD(5) - t86 * t23 + t89 * t24;
t3 = t21 * qJD(5) + t89 * t23 + t86 * t24;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, qJ(2) * t130, 0, 0, 0, 0, 0, 0, t172, t139 * t49 + t35, 0, -t36 * t120 + t35 * t147 + t52 * t49 - t51 * t50, t65, -0.2e1 * t55, -t153, t64, t154, 0, t87 * t103 - t105 * t90, t90 * t103 + t105 * t87, -t123 * t139 - t124, t41 * t123 + t61 * t124 + t36 * t60 + t40 * t50, t112, -t96, t164, -t111, t163, 0, -t26 * t56 + t39 * t46 + t4 * t82 - t109, -t25 * t56 - t3 * t82 - t39 * t48 - t110, t21 * t25 + t22 * t26 - t3 * t46 + t4 * t48 - t97, t1 * t22 + t2 * t21 + t29 * t56 + t3 * t6 + t32 * t39 + t4 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, -t94 * qJ(2), 0, 0, 0, 0, 0, 0, -t88 * t182, -t91 * t182, 0, t175 * t88 + t179 * t91, 0, 0, 0, 0, 0, 0, t87 * t102 - t90 * t183, t90 * t102 + t87 * t183, t98 * t174, (t124 - t180) * t88 + (-t98 * t41 - t36) * t91, 0, 0, 0, 0, 0, 0, -t46 * t117 + t151 * t82 + t91 * t26, t48 * t117 - t152 * t82 + t91 * t25, t151 * t48 - t152 * t46 - t53 * t25 - t54 * t26, -t1 * t54 - t32 * t117 + t151 * t5 + t152 * t6 - t2 * t53 - t29 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, -t175, 0, 0, t64, 0.2e1 * t55, t153, t65, -t154, 0, t87 * t104 - t107 * t90, t90 * t104 + t107 * t87, t125 * t139 + t124, -t36 * pkin(3) + pkin(7) * t124 - t41 * t125 - t40 * t52, -t112, t96, -t164, t111, -t163, 0, t114 * t46 + t149 * t82 - t77 * t26 + t109, -t114 * t48 - t150 * t82 - t77 * t25 + t110, t149 * t48 - t150 * t46 + t25 * t37 + t26 * t38 + t97, t1 * t38 + t114 * t32 + t149 * t5 + t150 * t6 + t2 * t37 + t29 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t136, t146 * t182, 0, t136, 0, 0, t121 * t87, t121 * t90, 0, 0, -t162, t14, t10, t162, t11, 0, t46 * t138 + t165 - t7 * t82 + (t129 * t86 - t166) * qJD(5) + t127, -t48 * t138 + t8 * t82 + (qJD(5) * t129 - t12) * t89 + t106, -(t6 + t7) * t48 + (-t5 + t8) * t46 + (t25 * t89 + t26 * t86 + (-t46 * t89 - t48 * t86) * qJD(5)) * pkin(4), -t5 * t7 - t6 * t8 + (t32 * t158 + t1 * t86 + t2 * t89 + (-t5 * t86 + t6 * t89) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t162, t14, t10, t162, t11, 0, t6 * t82 + t165 + t2, t5 * t82 + t106 - t173, 0, 0;];
tauc_reg = t9;
