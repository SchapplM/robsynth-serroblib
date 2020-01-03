% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRP1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:59:12
% EndTime: 2020-01-03 11:59:15
% DurationCPUTime: 0.75s
% Computational Cost: add. (1170->166), mult. (2398->221), div. (0->0), fcn. (1328->6), ass. (0->117)
t119 = pkin(1) * qJD(1);
t78 = cos(qJ(2));
t106 = t78 * t119;
t76 = sin(qJ(2));
t108 = t76 * t119;
t73 = sin(pkin(8));
t58 = t73 * t108;
t74 = cos(pkin(8));
t44 = t74 * t106 - t58;
t38 = t44 * qJD(2);
t139 = -qJD(3) * qJD(4) - t38;
t75 = sin(qJ(4));
t112 = t75 * qJD(3);
t77 = cos(qJ(4));
t70 = qJD(1) + qJD(2);
t52 = t70 * pkin(2) + t106;
t28 = t74 * t108 + t73 * t52;
t23 = t70 * pkin(7) + t28;
t98 = qJ(5) * t70 + t23;
t86 = t98 * t77;
t11 = t86 + t112;
t114 = qJD(4) * t75;
t109 = t23 * t114 + t139 * t77;
t66 = t77 * qJD(3);
t15 = -t75 * t23 + t66;
t16 = t77 * t23 + t112;
t6 = -t16 * qJD(4) - t75 * t38;
t138 = -t109 * t77 + (-t15 * t77 - t16 * t75) * qJD(4) - t6 * t75;
t71 = t75 ^ 2;
t137 = pkin(4) * t71;
t136 = t74 * pkin(2);
t135 = t77 * pkin(4);
t10 = -t98 * t75 + t66;
t117 = qJD(4) * pkin(4);
t9 = t10 + t117;
t134 = -t10 + t9;
t127 = t74 * t76;
t87 = t73 * t78 + t127;
t84 = pkin(1) * t87;
t43 = qJD(2) * t84;
t133 = t43 * t70;
t118 = pkin(1) * qJD(2);
t128 = t73 * t76;
t45 = (t74 * t78 - t128) * t118;
t132 = t45 * t70;
t69 = t70 ^ 2;
t131 = t69 * t77;
t130 = t70 * t75;
t129 = t70 * t77;
t79 = qJD(4) ^ 2;
t126 = t79 * t75;
t113 = qJD(4) * t77;
t103 = -pkin(3) - t135;
t27 = t74 * t52 - t58;
t14 = t103 * t70 + qJD(5) - t27;
t105 = t70 * t114;
t37 = qJD(1) * t43;
t21 = pkin(4) * t105 + t37;
t125 = t14 * t113 + t21 * t75;
t22 = -t70 * pkin(3) - t27;
t124 = t22 * t113 + t37 * t75;
t42 = qJD(1) * t84;
t123 = t44 * t114 + t42 * t129;
t64 = t78 * pkin(1) + pkin(2);
t122 = pkin(1) * t127 + t73 * t64;
t72 = t77 ^ 2;
t121 = t71 - t72;
t120 = t71 + t72;
t41 = pkin(7) + t122;
t116 = -qJ(5) - t41;
t61 = t73 * pkin(2) + pkin(7);
t115 = -qJ(5) - t61;
t65 = t77 * qJD(5);
t111 = -qJD(5) - t14;
t107 = pkin(4) * t114;
t104 = t70 * t113;
t101 = qJ(5) * t114;
t2 = (-t101 + t65) * t70 - t109;
t3 = (-qJD(5) * t70 - t38) * t75 - t11 * qJD(4);
t102 = t2 * t77 - t3 * t75;
t100 = t120 * t44;
t99 = -pkin(1) * t128 + t74 * t64;
t97 = qJD(4) * t116;
t96 = qJD(4) * t115;
t94 = t75 * t104;
t40 = -pkin(3) - t99;
t93 = (-qJD(2) + t70) * t119;
t92 = (-qJD(1) - t70) * t118;
t91 = t11 * t77 - t75 * t9;
t90 = -t11 * t75 - t77 * t9;
t89 = t15 * t75 - t16 * t77;
t88 = t41 * t79 + t133;
t85 = qJD(4) * (t40 * t70 - t45);
t83 = t87 * qJD(1) * t118;
t68 = t77 * qJ(5);
t67 = t79 * t77;
t62 = -pkin(3) - t136;
t55 = t75 * t131;
t53 = t103 - t136;
t51 = -0.2e1 * t94;
t50 = 0.2e1 * t94;
t49 = t121 * t69;
t48 = t77 * t61 + t68;
t47 = t115 * t75;
t39 = -0.2e1 * t121 * t70 * qJD(4);
t36 = -t75 * qJD(5) + t77 * t96;
t35 = t75 * t96 + t65;
t34 = t40 - t135;
t33 = t44 * t113;
t31 = t43 + t107;
t25 = t77 * t41 + t68;
t24 = t116 * t75;
t17 = t22 * t114;
t12 = t14 * t114;
t8 = (-qJD(5) - t45) * t75 + t77 * t97;
t7 = t77 * t45 + t75 * t97 + t65;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76 * t92, t78 * t92, 0, 0, 0, 0, 0, 0, 0, 0, -t83 - t133, -t38 - t132, 0, t38 * t122 - t27 * t43 + t28 * t45 - t37 * t99, t50, t39, t67, t51, -t126, 0, t17 + t75 * t85 + (-t37 - t88) * t77, t88 * t75 + t77 * t85 + t124, t120 * t132 + t138, t138 * t41 + t22 * t43 + t37 * t40 - t89 * t45, t50, t39, t67, t51, -t126, 0, t12 + (-t31 * t70 - t21) * t77 + (t34 * t130 + t8) * qJD(4), t31 * t130 + (t34 * t129 - t7) * qJD(4) + t125, (t7 * t77 - t75 * t8) * t70 + ((-t24 * t77 - t25 * t75) * t70 + t90) * qJD(4) + t102, t11 * t7 + t14 * t31 + t2 * t25 + t21 * t34 + t3 * t24 + t9 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76 * t93, t78 * t93, 0, 0, 0, 0, 0, 0, 0, 0, t42 * t70 - t83, t44 * t70 - t38, 0, t27 * t42 - t28 * t44 + (-t37 * t74 + t38 * t73) * pkin(2), t50, t39, t67, t51, -t126, 0, t62 * t105 + t17 + (-t61 * t79 - t37) * t77 + t123, t61 * t126 + t33 + (t62 * t113 - t42 * t75) * t70 + t124, -t70 * t100 + t138, t138 * t61 - t22 * t42 + t37 * t62 + t89 * t44, t50, t39, t67, t51, -t126, 0, -t21 * t77 + t12 + (t36 + (t53 - t135) * t130) * qJD(4) + t123, -t42 * t130 + t33 + (-t35 + (t53 * t77 + t137) * t70) * qJD(4) + t125, t90 * qJD(4) + (t35 * t77 - t36 * t75 - t100 + (-t47 * t77 - t48 * t75) * qJD(4)) * t70 + t102, t2 * t48 + t21 * t53 + t3 * t47 + (t44 * t75 + t36) * t9 + (-t42 + t107) * t14 + (-t44 * t77 + t35) * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126, -t67, 0, -t89 * qJD(4) - t109 * t75 + t6 * t77, 0, 0, 0, 0, 0, 0, -t126, -t67, 0, t91 * qJD(4) + t2 * t75 + t3 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, t49, 0, t55, 0, 0, (-t22 * t70 - t38) * t75, t15 * qJD(4) - t22 * t129 + t109, 0, 0, -t55, t49, 0, t55, 0, 0, (t11 - t86) * qJD(4) + (pkin(4) * t131 + t111 * t70 + t139) * t75, -t69 * t137 + t10 * qJD(4) + (t111 * t77 + t101) * t70 + t109, (-t117 + t134) * t129, t134 * t11 + (-t14 * t130 + t3) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t105, 0.2e1 * t104, -t120 * t69, -t91 * t70 + t21;];
tauc_reg = t1;
