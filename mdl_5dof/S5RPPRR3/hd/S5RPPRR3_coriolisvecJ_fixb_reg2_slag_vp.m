% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:15
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRR3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:14:36
% EndTime: 2022-01-23 09:14:39
% DurationCPUTime: 1.07s
% Computational Cost: add. (2384->179), mult. (5806->248), div. (0->0), fcn. (4334->8), ass. (0->107)
t102 = sin(qJ(5));
t142 = cos(qJ(5));
t104 = cos(qJ(4));
t100 = cos(pkin(9));
t125 = qJD(1) * t100;
t119 = t104 * t125;
t103 = sin(qJ(4));
t98 = sin(pkin(9));
t129 = t103 * t98;
t121 = qJD(1) * t129;
t74 = -t119 + t121;
t84 = t103 * t100 + t104 * t98;
t76 = t84 * qJD(1);
t111 = t102 * t74 - t142 * t76;
t88 = qJD(4) * t119;
t64 = qJD(4) * t121 - t88;
t79 = t84 * qJD(4);
t65 = qJD(1) * t79;
t107 = t111 * qJD(5) + t102 * t64 - t142 * t65;
t97 = qJD(4) + qJD(5);
t140 = t111 * t97;
t150 = t107 - t140;
t118 = qJD(5) * t142;
t123 = qJD(5) * t102;
t110 = -t102 * t65 - t74 * t118 - t76 * t123 - t142 * t64;
t38 = -t102 * t76 - t142 * t74;
t138 = t38 * t97;
t149 = t110 - t138;
t137 = t111 ^ 2;
t139 = t38 ^ 2;
t148 = t137 - t139;
t136 = t38 * t111;
t85 = -cos(pkin(8)) * pkin(1) - t100 * pkin(3) - pkin(2);
t71 = t85 * qJD(1) + qJD(3);
t42 = t74 * pkin(4) + t71;
t147 = t42 * t111;
t128 = qJD(3) * t98;
t120 = qJD(1) * t128;
t124 = qJD(4) * t104;
t92 = sin(pkin(8)) * pkin(1) + qJ(3);
t86 = t92 * qJD(1);
t94 = t100 * qJD(2);
t50 = t94 + (-pkin(6) * qJD(1) - t86) * t98;
t127 = t104 * t100;
t89 = qJD(3) * t127;
t133 = qJD(1) * t89 + t50 * t124;
t58 = t98 * qJD(2) + t100 * t86;
t51 = pkin(6) * t125 + t58;
t20 = (-qJD(4) * t51 - t120) * t103 + t133;
t14 = -t65 * pkin(7) + t20;
t109 = t84 * qJD(3);
t108 = qJD(1) * t109;
t30 = t103 * t50 + t104 * t51;
t21 = -t30 * qJD(4) - t108;
t15 = t64 * pkin(7) + t21;
t130 = t103 * t51;
t29 = t104 * t50 - t130;
t23 = -t76 * pkin(7) + t29;
t22 = qJD(4) * pkin(4) + t23;
t24 = -t74 * pkin(7) + t30;
t106 = -t102 * t15 - t22 * t118 + t24 * t123 - t142 * t14;
t146 = -t42 * t38 + t106;
t145 = t76 ^ 2;
t144 = pkin(6) + t92;
t78 = qJD(4) * t129 - t100 * t124;
t83 = t127 - t129;
t25 = t102 * t79 - t83 * t118 + t84 * t123 + t142 * t78;
t44 = t102 * t83 + t142 * t84;
t143 = t107 * t44 - t25 * t38;
t141 = t25 * t97;
t135 = t76 * t74;
t134 = -t84 * t65 + t78 * t74;
t80 = t144 * t98;
t81 = t144 * t100;
t41 = -t103 * t80 + t104 * t81;
t132 = t100 ^ 2 + t98 ^ 2;
t131 = t102 * t24;
t126 = t78 * qJD(4);
t122 = t142 * t24;
t117 = -t102 * t14 + t142 * t15;
t40 = -t103 * t81 - t104 * t80;
t115 = qJD(1) * t132;
t26 = t44 * qJD(5) - t102 * t78 + t142 * t79;
t43 = t102 * t84 - t142 * t83;
t114 = t110 * t43 - t111 * t26;
t113 = -t83 * t64 - t76 * t79;
t112 = -t100 * t58 + (-t98 * t86 + t94) * t98;
t33 = -t84 * pkin(7) + t40;
t34 = t83 * pkin(7) + t41;
t11 = -t102 * t34 + t142 * t33;
t6 = t102 * t22 + t122;
t12 = t102 * t33 + t142 * t34;
t31 = -t80 * t124 + t89 + (-qJD(4) * t81 - t128) * t103;
t2 = -t6 * qJD(5) + t117;
t32 = -t41 * qJD(4) - t109;
t72 = t74 ^ 2;
t69 = t79 * qJD(4);
t49 = -t83 * pkin(4) + t85;
t28 = t78 * pkin(7) + t32;
t27 = -t79 * pkin(7) + t31;
t19 = t26 * t97;
t8 = t142 * t23 - t131;
t7 = -t102 * t23 - t122;
t5 = t142 * t22 - t131;
t4 = -t12 * qJD(5) - t102 * t27 + t142 * t28;
t3 = t11 * qJD(5) + t102 * t28 + t142 * t27;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t115, (t92 * t115 - t112) * qJD(3), -t64 * t84 - t76 * t78, t113 + t134, -t126, -t65 * t83 + t74 * t79, -t69, 0, t32 * qJD(4) + t85 * t65 + t71 * t79, -t31 * qJD(4) - t85 * t64 - t71 * t78, t20 * t83 - t21 * t84 + t29 * t78 - t30 * t79 - t31 * t74 - t32 * t76 + t40 * t64 - t41 * t65, t20 * t41 + t21 * t40 + t29 * t32 + t30 * t31, t110 * t44 + t111 * t25, -t114 + t143, -t141, -t107 * t43 - t26 * t38, -t19, 0, -t49 * t107 + t42 * t26 + t4 * t97 + (-t38 * t79 + t43 * t65) * pkin(4), t49 * t110 - t42 * t25 - t3 * t97 + (-t111 * t79 + t44 * t65) * pkin(4), t106 * t43 + t107 * t12 - t11 * t110 + t111 * t4 - t2 * t44 + t5 * t25 - t6 * t26 + t3 * t38, -t106 * t12 + t2 * t11 + t6 * t3 + t5 * t4 + (t42 * t79 + t49 * t65) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, t126, -t113 + t134, t20 * t84 + t21 * t83 - t29 * t79 - t30 * t78, 0, 0, 0, 0, 0, 0, -t19, t141, t114 + t143, -t106 * t44 - t2 * t43 - t6 * t25 - t5 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t132 * qJD(1) ^ 2, t112 * qJD(1), 0, 0, 0, 0, 0, 0, 0.2e1 * t76 * qJD(4), t88 + (-t74 - t121) * qJD(4), -t72 - t145, t29 * t76 + t30 * t74, 0, 0, 0, 0, 0, 0, -t107 - t140, t110 + t138, -t137 - t139, t65 * pkin(4) - t111 * t5 - t6 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, -t72 + t145, t88 + (t74 - t121) * qJD(4), -t135, 0, 0, -t71 * t76 - t108, t103 * t120 + t71 * t74 + (t29 + t130) * qJD(4) - t133, 0, 0, t136, t148, t149, -t136, t150, 0, t76 * pkin(4) * t38 + t147 - t7 * t97 + (-t122 + (-pkin(4) * t97 - t22) * t102) * qJD(5) + t117, t8 * t97 + (t111 * t76 - t97 * t118) * pkin(4) + t146, -t6 * t111 - t8 * t38 + t5 * t38 - t7 * t111 + (-t142 * t110 + t102 * t107 + (-t102 * t111 + t142 * t38) * qJD(5)) * pkin(4), -t5 * t7 - t6 * t8 + (t142 * t2 - t106 * t102 - t42 * t76 + (-t102 * t5 + t142 * t6) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, t148, t149, -t136, t150, 0, t6 * t97 + t147 + t2, t5 * t97 + t146, 0, 0;];
tauc_reg = t1;
