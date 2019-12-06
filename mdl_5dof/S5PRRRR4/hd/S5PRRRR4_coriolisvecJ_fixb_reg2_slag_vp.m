% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRR4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:58
% EndTime: 2019-12-05 17:08:02
% DurationCPUTime: 0.97s
% Computational Cost: add. (1847->183), mult. (3314->254), div. (0->0), fcn. (2037->6), ass. (0->140)
t99 = cos(qJ(4));
t143 = qJD(4) * t99;
t98 = cos(qJ(5));
t158 = t98 * t99;
t175 = -qJD(5) * t158 - t98 * t143;
t95 = sin(qJ(5));
t96 = sin(qJ(4));
t64 = t95 * t99 + t98 * t96;
t92 = qJD(2) + qJD(3);
t56 = t64 * t92;
t147 = pkin(2) * qJD(2);
t97 = sin(qJ(3));
t134 = t97 * t147;
t70 = t92 * pkin(7) + t134;
t123 = pkin(8) * t92 + t70;
t114 = t123 * t96;
t100 = cos(qJ(3));
t146 = pkin(2) * qJD(3);
t126 = qJD(2) * t146;
t115 = t100 * t126;
t87 = t99 * qJD(1);
t150 = qJD(4) * t87 + t99 * t115;
t28 = -qJD(4) * t114 + t150;
t43 = t87 - t114;
t38 = qJD(4) * pkin(4) + t43;
t174 = (qJD(5) * t38 + t28) * t98;
t93 = t96 ^ 2;
t94 = t99 ^ 2;
t173 = (t93 + t94) * t92;
t144 = qJD(4) * t96;
t35 = -t70 * t144 + t150;
t108 = t96 * t115;
t142 = t96 * qJD(1);
t52 = t99 * t70 + t142;
t36 = -t52 * qJD(4) - t108;
t161 = t96 * t70;
t51 = t87 - t161;
t172 = t35 * t99 - t36 * t96 + (-t51 * t99 - t52 * t96) * qJD(4);
t109 = t51 * t96 - t52 * t99;
t171 = t109 * t100;
t91 = qJD(4) + qJD(5);
t44 = t123 * t99 + t142;
t170 = -pkin(8) - pkin(7);
t169 = t99 * pkin(4);
t84 = t97 * pkin(2) + pkin(7);
t168 = -pkin(8) - t84;
t162 = t95 * t96;
t111 = t91 * t162;
t39 = t111 + t175;
t33 = t39 * t91;
t137 = t92 * t162;
t54 = -t92 * t158 + t137;
t167 = t56 * t54;
t129 = t100 * t147;
t86 = -pkin(3) - t169;
t57 = t86 * t92 - t129;
t166 = t57 * t56;
t71 = -t92 * pkin(3) - t129;
t165 = t71 * t92;
t164 = t92 * t96;
t163 = t95 * t44;
t160 = t97 * t99;
t159 = t98 * t44;
t78 = t170 * t96;
t89 = t99 * pkin(8);
t79 = t99 * pkin(7) + t89;
t45 = t98 * t78 - t95 * t79;
t63 = -t158 + t162;
t122 = qJD(4) * t170;
t65 = t96 * t122;
t66 = t99 * t122;
t157 = t45 * qJD(5) + t63 * t129 + t98 * t65 + t95 * t66;
t46 = t95 * t78 + t98 * t79;
t156 = -t46 * qJD(5) + t64 * t129 - t95 * t65 + t98 * t66;
t40 = t91 * t64;
t27 = t40 * t92;
t155 = -t64 * t27 + t39 * t54;
t116 = t97 * t126;
t132 = t92 * t144;
t58 = pkin(4) * t132 + t116;
t154 = t57 * t40 + t58 * t63;
t153 = -t57 * t39 + t58 * t64;
t152 = t96 * t116 + t71 * t143;
t151 = t175 * t92;
t149 = t93 - t94;
t101 = qJD(4) ^ 2;
t145 = t101 * t96;
t88 = t101 * t99;
t141 = -qJD(2) - t92;
t140 = -qJD(3) + t92;
t139 = qJD(4) * t100;
t138 = pkin(4) * t164;
t90 = t92 ^ 2;
t136 = t96 * t90 * t99;
t135 = pkin(4) * t144;
t133 = t100 * t146;
t131 = t92 * t143;
t14 = t98 * t38 - t163;
t15 = t95 * t38 + t159;
t29 = -t44 * qJD(4) - t108;
t120 = -qJD(5) * t163 + t95 * t29;
t4 = t120 + t174;
t121 = -t95 * t28 + t98 * t29;
t5 = -t15 * qJD(5) + t121;
t127 = t14 * t39 - t15 * t40 - t4 * t63 - t5 * t64;
t125 = t96 * t139;
t85 = -t100 * pkin(2) - pkin(3);
t124 = -pkin(4) * t91 - t38;
t119 = qJD(4) * t168;
t117 = t96 * t131;
t113 = t140 * t147;
t112 = t141 * t146;
t26 = t92 * t111 + t151;
t110 = -t63 * t26 + t56 * t40;
t61 = t168 * t96;
t62 = t99 * t84 + t89;
t31 = t98 * t61 - t95 * t62;
t32 = t95 * t61 + t98 * t62;
t107 = t57 * t54 - t120;
t106 = -t134 + t135;
t105 = t99 * t139 - t97 * t164;
t74 = t85 - t169;
t69 = -0.2e1 * t117;
t68 = 0.2e1 * t117;
t67 = t97 * t146 + t135;
t59 = t71 * t144;
t53 = -0.2e1 * t149 * t92 * qJD(4);
t48 = t99 * t119 - t96 * t133;
t47 = t96 * t119 + t99 * t133;
t34 = t40 * t91;
t18 = -t54 ^ 2 + t56 ^ 2;
t17 = t98 * t43 - t163;
t16 = -t95 * t43 - t159;
t12 = -t151 + (-t137 + t54) * t91;
t11 = -t32 * qJD(5) - t95 * t47 + t98 * t48;
t10 = t31 * qJD(5) + t98 * t47 + t95 * t48;
t7 = t27 * t63 + t54 * t40;
t6 = -t26 * t64 - t56 * t39;
t1 = -t110 + t155;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t145, -t88, 0, -t109 * qJD(4) + t35 * t96 + t36 * t99, 0, 0, 0, 0, 0, 0, -t34, t33, t110 + t155, -t14 * t40 - t15 * t39 + t4 * t64 - t5 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97 * t112, t100 * t112, 0, 0, t68, t53, t88, t69, -t145, 0, t85 * t132 - t84 * t88 + t59 + (t141 * t160 - t125) * t146, -t105 * t146 + t85 * t131 + t84 * t145 + t152, t133 * t173 + t172, t172 * t84 + ((qJD(2) * t85 + t71) * t97 - t171) * t146, t6, t1, -t33, t7, -t34, 0, t11 * t91 + t74 * t27 + t67 * t54 + t154, -t10 * t91 - t74 * t26 + t67 * t56 + t153, -t10 * t54 - t11 * t56 + t31 * t26 - t32 * t27 + t127, t15 * t10 + t14 * t11 + t5 * t31 + t4 * t32 + t57 * t67 + t58 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97 * t113, t100 * t113, 0, 0, t68, t53, t88, t69, -t145, 0, -pkin(3) * t132 - pkin(7) * t88 + t59 + (t140 * t160 + t125) * t147, -pkin(3) * t131 + pkin(7) * t145 + t105 * t147 + t152, -t129 * t173 + t172, t172 * pkin(7) + ((-pkin(3) * qJD(3) - t71) * t97 + t171) * t147, t6, t1, -t33, t7, -t34, 0, t106 * t54 + t156 * t91 + t86 * t27 + t154, t106 * t56 - t157 * t91 - t86 * t26 + t153, -t156 * t56 - t157 * t54 + t45 * t26 - t46 * t27 + t127, t106 * t57 + t156 * t14 + t157 * t15 + t4 * t46 + t5 * t45 + t58 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t136, t149 * t90, 0, t136, 0, 0, (-t115 - t165) * t96, -t99 * t165 + (t51 + t161) * qJD(4) - t150, 0, 0, t167, t18, t12, -t167, 0, 0, -t54 * t138 - t16 * t91 - t166 + (t124 * t95 - t159) * qJD(5) + t121, -t56 * t138 + t17 * t91 + (t124 * qJD(5) - t28) * t98 + t107, (t15 + t16) * t56 + (-t14 + t17) * t54 + (t26 * t98 - t27 * t95 + (-t54 * t98 + t56 * t95) * qJD(5)) * pkin(4), -t14 * t16 - t15 * t17 + (-t57 * t164 + t4 * t95 + t5 * t98 + (-t14 * t95 + t15 * t98) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t167, t18, t12, -t167, 0, 0, t15 * t91 - t166 + t5, t14 * t91 + t107 - t174, 0, 0;];
tauc_reg = t2;
