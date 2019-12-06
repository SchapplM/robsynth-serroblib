% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRR3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:47:32
% EndTime: 2019-12-05 15:47:40
% DurationCPUTime: 1.09s
% Computational Cost: add. (1582->178), mult. (3771->257), div. (0->0), fcn. (2771->8), ass. (0->124)
t81 = sin(qJ(5));
t82 = sin(qJ(4));
t84 = cos(qJ(5));
t85 = cos(qJ(4));
t61 = t81 * t85 + t84 * t82;
t76 = qJD(4) + qJD(5);
t144 = t76 * t61;
t32 = t144 * qJD(2);
t79 = sin(pkin(9));
t80 = cos(pkin(9));
t83 = sin(qJ(2));
t86 = cos(qJ(2));
t59 = t79 * t86 + t80 * t83;
t130 = t84 * t85;
t134 = t81 * t82;
t60 = -t130 + t134;
t26 = t60 * t59;
t58 = t79 * t83 - t80 * t86;
t51 = t58 * qJD(2);
t46 = qJD(1) * t51;
t75 = t85 * qJD(3);
t125 = qJD(4) * t75 - t85 * t46;
t121 = qJD(1) * t83;
t114 = t86 * qJD(1);
t66 = qJD(2) * pkin(2) + t114;
t44 = t80 * t121 + t79 * t66;
t41 = qJD(2) * pkin(6) + t44;
t101 = pkin(7) * qJD(2) + t41;
t94 = t101 * t82;
t12 = -qJD(4) * t94 + t125;
t27 = t75 - t94;
t24 = qJD(4) * pkin(4) + t27;
t143 = (qJD(5) * t24 + t12) * t84;
t50 = t59 * qJD(1);
t115 = t82 * qJD(3);
t28 = t101 * t85 + t115;
t142 = t80 * pkin(2);
t71 = t79 * pkin(2) + pkin(6);
t141 = pkin(7) + t71;
t56 = t141 * t82;
t57 = t141 * t85;
t29 = -t84 * t56 - t81 * t57;
t102 = qJD(4) * t141;
t47 = t82 * t102;
t48 = t85 * t102;
t67 = t79 * t121;
t52 = t80 * t114 - t67;
t140 = -t29 * qJD(5) + t84 * t47 + t81 * t48 - t60 * t52;
t117 = qJD(4) * t85;
t97 = t76 * t134;
t36 = -qJD(5) * t130 - t84 * t117 + t97;
t139 = t36 * t76;
t107 = -t85 * pkin(4) - pkin(3);
t43 = t80 * t66 - t67;
t38 = t107 * qJD(2) - t43;
t55 = t61 * qJD(2);
t138 = t38 * t55;
t49 = t59 * qJD(2);
t45 = qJD(1) * t49;
t137 = t45 * t58;
t119 = qJD(2) * t85;
t108 = t84 * t119;
t120 = qJD(2) * t82;
t109 = t81 * t120;
t53 = -t108 + t109;
t136 = t55 * t53;
t135 = t81 * t28;
t133 = t82 * t41;
t132 = t82 * t46;
t131 = t84 * t28;
t87 = qJD(4) ^ 2;
t129 = t87 * t82;
t128 = t87 * t85;
t30 = -t81 * t56 + t84 * t57;
t127 = -t30 * qJD(5) + t81 * t47 - t84 * t48 + t61 * t52;
t126 = -t61 * t32 + t36 * t53;
t113 = qJD(2) * qJD(4);
t105 = t85 * t113;
t124 = -qJD(5) * t108 - t84 * t105;
t77 = t82 ^ 2;
t78 = t85 ^ 2;
t123 = t77 - t78;
t122 = t77 + t78;
t118 = qJD(4) * t82;
t116 = t51 * qJD(2);
t88 = qJD(2) ^ 2;
t112 = t82 * t88 * t85;
t111 = pkin(4) * t120;
t110 = pkin(4) * t118;
t106 = -pkin(4) * t76 - t24;
t13 = -t28 * qJD(4) + t132;
t104 = -t81 * t12 + t84 * t13;
t103 = -qJD(5) * t135 + t81 * t13;
t99 = t82 * t105;
t98 = -t50 + t110;
t6 = t81 * t24 + t131;
t31 = qJD(2) * t97 + t124;
t96 = t144 * t55 - t60 * t31;
t34 = t75 - t133;
t35 = t85 * t41 + t115;
t95 = t34 * t82 - t35 * t85;
t93 = t38 * t53 - t103;
t92 = t50 * qJD(2) - t71 * t87 - t45;
t40 = -qJD(2) * pkin(3) - t43;
t72 = -pkin(3) - t142;
t91 = qJD(4) * (qJD(2) * t72 + t40 + t52);
t2 = -t6 * qJD(5) + t104;
t14 = -t41 * t118 + t125;
t15 = -t35 * qJD(4) + t132;
t90 = t14 * t85 - t15 * t82 + (-t34 * t85 - t35 * t82) * qJD(4);
t63 = t107 - t142;
t39 = (t50 + t110) * qJD(2);
t33 = t144 * t76;
t25 = t61 * t59;
t19 = -t53 ^ 2 + t55 ^ 2;
t17 = t55 * t76 - t32;
t16 = -t124 + (-t109 + t53) * t76;
t8 = t84 * t27 - t135;
t7 = -t81 * t27 - t131;
t5 = t84 * t24 - t135;
t4 = t76 * t26 + t61 * t51;
t3 = -t144 * t59 + t60 * t51;
t1 = t103 + t143;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88 * t83, -t88 * t86, 0, 0, 0, 0, 0, 0, 0, 0, -t49 * qJD(2), t116, 0, -t43 * t49 - t44 * t51 - t46 * t59 + t137, 0, 0, 0, 0, 0, 0, t51 * t118 - t59 * t128 + (t58 * t118 - t49 * t85) * qJD(2), t51 * t117 + t59 * t129 + (t58 * t117 + t49 * t82) * qJD(2), -t122 * t116, t40 * t49 + t95 * t51 + t90 * t59 + t137, 0, 0, 0, 0, 0, 0, t58 * t32 + t4 * t76 + t49 * t53, -t3 * t76 - t58 * t31 + t49 * t55, -t25 * t31 + t26 * t32 - t3 * t53 - t4 * t55, -t1 * t26 - t2 * t25 + t6 * t3 + t38 * t49 + t39 * t58 + t5 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t58 * qJD(1) + t52) * qJD(2), 0, t43 * t50 - t44 * t52 + (-t45 * t80 - t46 * t79) * pkin(2), 0.2e1 * t99, -0.2e1 * t123 * t113, t128, -0.2e1 * t99, -t129, 0, t82 * t91 + t92 * t85, -t92 * t82 + t85 * t91, -t122 * t52 * qJD(2) + t90, -t40 * t50 + t45 * t72 + t95 * t52 + t90 * t71, -t31 * t61 - t55 * t36, -t96 + t126, -t139, t144 * t53 + t32 * t60, -t33, 0, t127 * t76 + t144 * t38 + t63 * t32 + t39 * t60 + t98 * t53, t140 * t76 - t63 * t31 - t38 * t36 + t39 * t61 + t98 * t55, -t1 * t60 - t127 * t55 + t140 * t53 - t144 * t6 - t2 * t61 + t29 * t31 - t30 * t32 + t5 * t36, t1 * t30 + t127 * t5 - t140 * t6 + t2 * t29 + t98 * t38 + t39 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t129, -t128, 0, -t95 * qJD(4) + t14 * t82 + t15 * t85, 0, 0, 0, 0, 0, 0, -t33, t139, t96 + t126, t1 * t61 - t144 * t5 - t2 * t60 - t6 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112, t123 * t88, 0, t112, 0, 0, (-qJD(2) * t40 + t46) * t82, -t40 * t119 + (t34 + t133) * qJD(4) - t125, 0, 0, t136, t19, t16, -t136, t17, 0, -t53 * t111 - t138 - t7 * t76 + (t106 * t81 - t131) * qJD(5) + t104, -t55 * t111 + t8 * t76 + (t106 * qJD(5) - t12) * t84 + t93, (t6 + t7) * t55 + (-t5 + t8) * t53 + (t31 * t84 - t32 * t81 + (-t53 * t84 + t55 * t81) * qJD(5)) * pkin(4), -t5 * t7 - t6 * t8 + (-t38 * t120 + t1 * t81 + t2 * t84 + (-t5 * t81 + t6 * t84) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, t19, t16, -t136, t17, 0, t6 * t76 - t138 + t2, t5 * t76 - t143 + t93, 0, 0;];
tauc_reg = t9;
