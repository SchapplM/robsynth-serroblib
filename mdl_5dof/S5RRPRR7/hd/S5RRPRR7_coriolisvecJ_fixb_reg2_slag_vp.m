% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRR7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:40
% EndTime: 2019-12-31 20:15:44
% DurationCPUTime: 1.01s
% Computational Cost: add. (2180->196), mult. (3201->245), div. (0->0), fcn. (1824->6), ass. (0->137)
t141 = pkin(1) * qJD(1);
t97 = sin(qJ(2));
t126 = t97 * t141;
t92 = qJD(1) + qJD(2);
t54 = (qJD(2) - t92) * t126;
t95 = sin(qJ(5));
t132 = qJD(5) * t95;
t96 = sin(qJ(4));
t134 = qJD(4) * t96;
t165 = -t96 * t132 - t95 * t134;
t98 = cos(qJ(5));
t99 = cos(qJ(4));
t66 = t95 * t99 + t96 * t98;
t49 = t66 * t92;
t140 = pkin(1) * qJD(2);
t120 = qJD(1) * t140;
t112 = t97 * t120;
t101 = -pkin(2) - pkin(7);
t100 = cos(qJ(2));
t123 = t100 * t141;
t109 = qJD(3) - t123;
t53 = t101 * t92 + t109;
t117 = pkin(8) * t92 - t53;
t133 = qJD(4) * t99;
t27 = t112 * t96 - t117 * t133;
t155 = t92 * t99;
t39 = -pkin(8) * t155 + t99 * t53;
t35 = qJD(4) * pkin(4) + t39;
t164 = (qJD(5) * t35 + t27) * t98;
t89 = t96 * pkin(4);
t86 = qJ(3) + t89;
t127 = t97 * t140;
t93 = t96 ^ 2;
t94 = t99 ^ 2;
t142 = t93 + t94;
t163 = t127 * t142;
t162 = qJD(1) + t92;
t91 = qJD(4) + qJD(5);
t156 = t92 * t96;
t38 = -pkin(8) * t156 + t53 * t96;
t154 = t95 * t38;
t12 = t35 * t98 - t154;
t152 = t98 * t38;
t13 = t35 * t95 + t152;
t36 = t91 * t66;
t151 = t98 * t99;
t110 = t91 * t151;
t37 = t110 + t165;
t76 = t99 * t112;
t26 = t117 * t134 + t76;
t115 = -t38 * t132 + t95 * t26;
t4 = t115 + t164;
t116 = t26 * t98 - t95 * t27;
t5 = -qJD(5) * t13 + t116;
t153 = t95 * t96;
t67 = t151 - t153;
t121 = t12 * t36 - t13 * t37 - t4 * t66 - t5 * t67;
t90 = t92 ^ 2;
t119 = -pkin(1) * t100 - pkin(2);
t84 = -pkin(7) + t119;
t160 = -pkin(8) + t84;
t34 = t37 * t91;
t129 = t92 * t151;
t51 = -t153 * t92 + t129;
t159 = t51 * t49;
t52 = t86 * t92 + t126;
t158 = t52 * t51;
t136 = t92 * qJ(3);
t68 = t126 + t136;
t157 = t68 * t92;
t150 = -pkin(8) + t101;
t73 = t150 * t96;
t74 = t150 * t99;
t40 = -t73 * t95 + t74 * t98;
t87 = pkin(8) * t134;
t61 = -t101 * t134 + t87;
t62 = qJD(4) * t74;
t149 = qJD(5) * t40 - t66 * t126 + t95 * t61 + t98 * t62;
t41 = t73 * t98 + t74 * t95;
t148 = -qJD(5) * t41 - t67 * t126 + t98 * t61 - t95 * t62;
t88 = pkin(4) * t133;
t79 = qJD(3) + t88;
t83 = t100 * t120;
t44 = t79 * t92 + t83;
t147 = t37 * t52 + t44 * t66;
t146 = -t36 * t52 + t44 * t67;
t64 = qJD(3) * t92 + t83;
t145 = t68 * t133 + t64 * t96;
t144 = t165 * t92;
t143 = t93 - t94;
t139 = t100 * t68;
t102 = qJD(4) ^ 2;
t138 = t102 * t96;
t137 = t102 * t99;
t135 = -t102 - t90;
t131 = t101 * t102;
t130 = pkin(4) * t155;
t128 = t99 * t90 * t96;
t125 = t100 * t140;
t60 = t160 * t99;
t118 = -pkin(4) * t91 - t35;
t85 = pkin(1) * t97 + qJ(3);
t113 = t133 * t156;
t111 = t79 - t123;
t78 = qJD(3) + t125;
t24 = t36 * t92;
t6 = -t24 * t67 - t36 * t51;
t25 = t110 * t92 + t144;
t7 = t25 * t66 + t37 * t49;
t59 = t160 * t96;
t29 = t59 * t98 + t60 * t95;
t28 = -t59 * t95 + t60 * t98;
t108 = t64 * t85 + t68 * t78;
t107 = -t102 * t84 + t78 * t92;
t106 = t64 * qJ(3) + t68 * qJD(3);
t105 = t52 * t49 - t115;
t104 = t85 * t92 + t127;
t103 = t112 - t157;
t75 = t85 + t89;
t70 = -0.2e1 * t113;
t69 = 0.2e1 * t113;
t65 = -pkin(2) * t92 + t109;
t63 = t78 + t88;
t57 = t64 * t99;
t55 = t162 * t127;
t48 = 0.2e1 * t143 * t92 * qJD(4);
t43 = qJD(4) * t60 + t127 * t96;
t42 = t127 * t99 - t134 * t84 + t87;
t33 = t36 * t91;
t20 = -t49 ^ 2 + t51 ^ 2;
t17 = t39 * t98 - t154;
t16 = -t39 * t95 - t152;
t15 = -t144 + (-t129 + t51) * t91;
t11 = -qJD(5) * t29 + t98 * t42 - t95 * t43;
t10 = qJD(5) * t28 + t95 * t42 + t98 * t43;
t1 = t24 * t66 - t25 * t67 + t36 * t49 - t37 * t51;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, -t125 * t92 - t83, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t83 + (qJD(3) + t78) * t92, (qJD(1) * t119 + t65) * t127 + t108, t70, t48, -t138, t69, -t137, 0, t104 * t133 + t107 * t96 + t145, t57 + t107 * t99 + (-t104 - t68) * t134, -t162 * t163, t108 + (qJD(1) * t84 + t53) * t163, t6, t1, -t33, t7, -t34, 0, t11 * t91 + t25 * t75 + t49 * t63 + t147, -t10 * t91 - t24 * t75 + t51 * t63 + t146, -t10 * t49 - t11 * t51 + t24 * t28 - t25 * t29 + t121, t10 * t13 + t11 * t12 + t28 * t5 + t29 * t4 + t44 * t75 + t52 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t123 * t92 - t83, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t83 + (0.2e1 * qJD(3) - t123) * t92, (-t139 + (-pkin(2) * qJD(2) - t65) * t97) * t141 + t106, t70, t48, -t138, t69, -t137, 0, -t126 * t133 - t96 * t131 + (qJ(3) * t133 + t109 * t96) * t92 + t145, t57 + (t109 * t92 - t131) * t99 + (t126 - t68 - t136) * t134, -t142 * t54, (-t139 + (qJD(2) * t101 - t53) * t97 * t142) * t141 + t106, t6, t1, -t33, t7, -t34, 0, t111 * t49 + t148 * t91 + t86 * t25 + t147, t111 * t51 - t149 * t91 - t86 * t24 + t146, -t148 * t51 - t149 * t49 + t40 * t24 - t41 * t25 + t121, t111 * t52 + t12 * t148 + t13 * t149 + t4 * t41 + t5 * t40 + t44 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t90, t103, 0, 0, 0, 0, 0, 0, t135 * t96, t135 * t99, 0, t112 * t142 - t157, 0, 0, 0, 0, 0, 0, -t49 * t92 - t33, -t51 * t92 - t34, -t6 - t7, -t52 * t92 - t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, -t143 * t90, 0, -t128, 0, 0, -t155 * t68 + t76, -t103 * t96, 0, 0, t159, t20, 0, -t159, t15, 0, -t49 * t130 - t16 * t91 - t158 + (t118 * t95 - t152) * qJD(5) + t116, -t51 * t130 + t17 * t91 + (qJD(5) * t118 - t27) * t98 + t105, (t13 + t16) * t51 + (-t12 + t17) * t49 + (t24 * t98 - t25 * t95 + (-t49 * t98 + t51 * t95) * qJD(5)) * pkin(4), -t12 * t16 - t13 * t17 + (-t52 * t155 + t4 * t95 + t5 * t98 + (-t12 * t95 + t13 * t98) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159, t20, 0, -t159, t15, 0, t13 * t91 - t158 + t5, t12 * t91 + t105 - t164, 0, 0;];
tauc_reg = t2;
