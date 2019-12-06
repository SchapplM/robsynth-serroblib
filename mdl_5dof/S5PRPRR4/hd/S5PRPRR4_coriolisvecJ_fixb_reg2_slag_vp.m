% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRR4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:51:35
% EndTime: 2019-12-05 15:51:41
% DurationCPUTime: 1.44s
% Computational Cost: add. (1857->224), mult. (4891->343), div. (0->0), fcn. (3757->10), ass. (0->135)
t67 = sin(pkin(5));
t66 = sin(pkin(10));
t68 = cos(pkin(10));
t72 = sin(qJ(2));
t75 = cos(qJ(2));
t87 = t66 * t75 + t68 * t72;
t45 = t87 * t67;
t39 = qJD(1) * t45;
t71 = sin(qJ(4));
t74 = cos(qJ(4));
t94 = pkin(4) * t71 - pkin(8) * t74;
t53 = t94 * qJD(4);
t25 = (t53 + t39) * qJD(2);
t44 = (t66 * t72 - t68 * t75) * t67;
t41 = qJD(2) * t44;
t36 = qJD(1) * t41;
t125 = qJD(1) * t67;
t112 = t72 * t125;
t105 = t75 * t125;
t54 = qJD(2) * pkin(2) + t105;
t34 = t68 * t112 + t66 * t54;
t32 = qJD(2) * pkin(7) + t34;
t69 = cos(pkin(5));
t57 = t69 * qJD(1) + qJD(3);
t88 = t32 * t71 - t57 * t74;
t7 = -t88 * qJD(4) - t74 * t36;
t70 = sin(qJ(5));
t73 = cos(qJ(5));
t22 = t32 * t74 + t57 * t71;
t20 = qJD(4) * pkin(8) + t22;
t55 = t66 * t112;
t33 = t68 * t54 - t55;
t84 = -t74 * pkin(4) - t71 * pkin(8) - pkin(3);
t24 = t84 * qJD(2) - t33;
t90 = t20 * t70 - t24 * t73;
t1 = -t90 * qJD(5) + t70 * t25 + t73 * t7;
t116 = t74 * qJD(2);
t58 = -qJD(5) + t116;
t156 = -t90 * t58 + t1;
t155 = -t39 + t53;
t6 = t20 * t73 + t24 * t70;
t2 = -qJD(5) * t6 + t73 * t25 - t70 * t7;
t154 = t6 * t58 - t2;
t120 = qJD(5) * t70;
t108 = t71 * t120;
t117 = t73 * qJD(4);
t64 = t71 ^ 2;
t85 = qJD(2) * t64 - t58 * t74;
t153 = -t58 * t108 - t85 * t117;
t114 = qJD(4) * qJD(5);
t119 = qJD(5) * t73;
t107 = t71 * t119;
t118 = t70 * qJD(4);
t80 = t74 * t118 + t107;
t38 = t80 * qJD(2) + t70 * t114;
t29 = t45 * t71 - t69 * t74;
t8 = t22 * qJD(4) - t71 * t36;
t152 = t29 * t8;
t149 = t68 * pkin(2);
t148 = t8 * t70;
t147 = t8 * t71;
t146 = t8 * t73;
t19 = -qJD(4) * pkin(4) + t88;
t145 = t19 * t70;
t144 = t19 * t73;
t40 = qJD(2) * t45;
t35 = qJD(1) * t40;
t143 = t35 * t44;
t124 = qJD(2) * t71;
t104 = t70 * t124;
t49 = t104 - t117;
t142 = t49 * t58;
t141 = t49 * t71;
t51 = t73 * t124 + t118;
t140 = t51 * t49;
t139 = t51 * t58;
t138 = t58 * t70;
t137 = t58 * t73;
t77 = qJD(2) ^ 2;
t136 = t67 * t77;
t135 = t70 * t74;
t134 = t71 * t73;
t133 = t73 * t74;
t132 = t74 * t38;
t76 = qJD(4) ^ 2;
t131 = t76 * t71;
t130 = t76 * t74;
t123 = qJD(4) * t71;
t60 = t66 * pkin(2) + pkin(7);
t110 = t60 * t123;
t48 = t84 - t149;
t26 = -t60 * t135 + t73 * t48;
t42 = t68 * t105 - t55;
t129 = t26 * qJD(5) - t73 * t110 - t42 * t133 + t155 * t70;
t27 = t60 * t133 + t70 * t48;
t128 = -t27 * qJD(5) + t70 * t110 + t42 * t135 + t155 * t73;
t109 = t74 * t117;
t127 = -t49 * t109 - t38 * t134;
t65 = t74 ^ 2;
t126 = t64 - t65;
t122 = qJD(4) * t74;
t121 = qJD(5) * t49;
t115 = qJD(2) * qJD(4);
t113 = t71 * t77 * t74;
t111 = t51 * t122;
t106 = t58 * t119;
t102 = t71 * t115;
t37 = -qJD(2) * t109 + qJD(5) * t104 - t73 * t114;
t101 = t51 * t123 + t37 * t74;
t31 = -qJD(2) * pkin(3) - t33;
t100 = -qJD(2) * t31 + t36;
t99 = -t37 + t121;
t97 = t71 * t106;
t96 = t51 * t107;
t95 = t74 * t102;
t93 = -t6 * t70 + t73 * t90;
t92 = -t6 * t73 - t70 * t90;
t89 = -t22 * t74 - t71 * t88;
t30 = t45 * t74 + t69 * t71;
t16 = t30 * t73 + t44 * t70;
t15 = -t30 * t70 + t44 * t73;
t83 = t39 * qJD(2) - t60 * t76 - t35;
t61 = -pkin(3) - t149;
t82 = qJD(4) * (qJD(2) * t61 + t31 + t42);
t81 = t85 * t70;
t79 = t93 * qJD(5) + t1 * t73 - t2 * t70;
t78 = t7 * t74 + t147 + (-t22 * t71 + t74 * t88) * qJD(4);
t52 = t94 * qJD(2);
t14 = -t29 * qJD(4) - t41 * t74;
t13 = t30 * qJD(4) - t41 * t71;
t10 = t52 * t70 - t73 * t88;
t9 = t52 * t73 + t70 * t88;
t4 = t15 * qJD(5) + t14 * t73 + t40 * t70;
t3 = -t16 * qJD(5) - t14 * t70 + t40 * t73;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72 * t136, -t75 * t136, 0, 0, 0, 0, 0, 0, 0, 0, -t40 * qJD(2), t41 * qJD(2), 0, -t33 * t40 - t34 * t41 - t36 * t45 + t143, 0, 0, 0, 0, 0, 0, -t13 * qJD(4) + (t44 * t123 - t40 * t74) * qJD(2), -t14 * qJD(4) + (t44 * t122 + t40 * t71) * qJD(2), (t13 * t71 + t14 * t74 + (t29 * t74 - t30 * t71) * qJD(4)) * qJD(2), t13 * t88 + t14 * t22 + t30 * t7 + t31 * t40 + t143 + t152, 0, 0, 0, 0, 0, 0, t102 * t15 + t13 * t49 + t29 * t38 - t3 * t58, -t102 * t16 + t13 * t51 - t29 * t37 + t4 * t58, t15 * t37 - t16 * t38 - t3 * t51 - t4 * t49, t1 * t16 + t13 * t19 + t15 * t2 - t3 * t90 + t4 * t6 + t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t87 * t125 + t39) * qJD(2), (qJD(1) * t44 + t42) * qJD(2), 0, t33 * t39 - t34 * t42 + (-t35 * t68 - t36 * t66) * pkin(2), 0.2e1 * t95, -0.2e1 * t126 * t115, t130, -0.2e1 * t95, -t131, 0, t71 * t82 + t83 * t74, -t83 * t71 + t74 * t82, (-t64 - t65) * t42 * qJD(2) + t78, -t31 * t39 + t35 * t61 + t89 * t42 + t78 * t60, -t37 * t134 + (-t108 + t109) * t51, -t96 + (-t111 + (t37 + t121) * t71) * t70 + t127, t101 - t153, t38 * t70 * t71 + t49 * t80, t97 + t132 + (-t81 - t141) * qJD(4), (-t58 - t116) * t123, -t128 * t58 + (-t2 + (t49 * t60 + t145) * qJD(4)) * t74 + (t19 * t119 + t38 * t60 - t42 * t49 + t148 + (qJD(2) * t26 - t90) * qJD(4)) * t71, t129 * t58 + (t1 + (t51 * t60 + t144) * qJD(4)) * t74 + (-t19 * t120 - t37 * t60 - t42 * t51 + t146 + (-qJD(2) * t27 - t6) * qJD(4)) * t71, t26 * t37 - t27 * t38 - t128 * t51 - t129 * t49 + t93 * t122 + (qJD(5) * t92 - t1 * t70 - t2 * t73) * t71, t60 * t147 + t1 * t27 + t2 * t26 + t129 * t6 - t128 * t90 + (t122 * t60 - t42 * t71) * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131, -t130, 0, -t89 * qJD(4) + t7 * t71 - t8 * t74, 0, 0, 0, 0, 0, 0, t97 - t132 + (-t81 + t141) * qJD(4), t101 + t153, t96 + (t71 * t99 + t111) * t70 + t127, (-qJD(4) * t92 - t8) * t74 + (qJD(4) * t19 + t79) * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113, t126 * t77, 0, t113, 0, 0, t100 * t71, t100 * t74, 0, 0, -t51 * t137 - t37 * t70, (-t37 + t142) * t73 + (-t38 + t139) * t70, -t106 + (t58 * t133 + (-t51 + t118) * t71) * qJD(2), -t49 * t138 - t38 * t73, t58 * t120 + (-t58 * t135 + (t49 + t117) * t71) * qJD(2), t58 * t124, -pkin(4) * t38 - t22 * t49 + t9 * t58 - t146 + (pkin(8) * t137 + t145) * qJD(5) + (t90 * t71 + (-pkin(8) * t123 - t19 * t74) * t70) * qJD(2), pkin(4) * t37 - t10 * t58 - t22 * t51 + t148 + (-pkin(8) * t138 + t144) * qJD(5) + (-t19 * t133 + (-pkin(8) * t117 + t6) * t71) * qJD(2), t10 * t49 + t9 * t51 + ((qJD(5) * t51 - t38) * pkin(8) + t156) * t73 + (pkin(8) * t99 + t154) * t70, -t8 * pkin(4) + pkin(8) * t79 - t6 * t10 - t19 * t22 + t9 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140, -t49 ^ 2 + t51 ^ 2, -t37 - t142, -t140, -t139 - t38, t102, -t19 * t51 - t154, t19 * t49 - t156, 0, 0;];
tauc_reg = t5;
