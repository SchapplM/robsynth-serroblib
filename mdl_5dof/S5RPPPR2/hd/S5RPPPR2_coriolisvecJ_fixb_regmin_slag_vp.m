% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPPR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 08:59:46
% EndTime: 2022-01-23 08:59:49
% DurationCPUTime: 0.95s
% Computational Cost: add. (775->162), mult. (2304->286), div. (0->0), fcn. (1800->8), ass. (0->110)
t74 = sin(pkin(8));
t120 = qJD(1) * t74;
t80 = cos(qJ(5));
t107 = t80 * t120;
t75 = sin(pkin(7));
t119 = qJD(1) * t75;
t76 = cos(pkin(9));
t102 = t76 * t119;
t78 = cos(pkin(7));
t118 = qJD(1) * t78;
t73 = sin(pkin(9));
t104 = t73 * t118;
t77 = cos(pkin(8));
t42 = t77 * t102 - t104;
t79 = sin(qJ(5));
t136 = t75 * t107 - t79 * t42;
t126 = t78 * t76;
t127 = t75 * t77;
t49 = t73 * t127 + t126;
t39 = t49 * qJD(1);
t35 = qJD(5) + t39;
t135 = qJD(5) - t35;
t116 = qJD(3) * t75;
t117 = qJD(2) * t78;
t51 = t77 * t116 + t74 * t117;
t45 = t51 * qJD(1);
t134 = t45 * t77;
t114 = qJ(2) * qJD(1);
t60 = t75 * t114 + qJD(3);
t133 = t60 * t75;
t71 = t75 ^ 2;
t81 = qJD(1) ^ 2;
t132 = t71 * t81;
t131 = t73 * t74;
t130 = t74 * t75;
t129 = t74 * t79;
t128 = t74 * t80;
t101 = t78 * t114;
t55 = -pkin(2) * t78 - t75 * qJ(3) - pkin(1);
t48 = t55 * qJD(1) + qJD(2);
t25 = t77 * t101 + t74 * t48;
t15 = -qJ(4) * t118 + t25;
t91 = pkin(3) * t74 - qJ(4) * t77;
t30 = t91 * t119 + t60;
t5 = t76 * t15 + t73 * t30;
t121 = qJ(2) * t78;
t123 = t77 * t121 + t74 * t55;
t29 = -t78 * qJ(4) + t123;
t36 = (qJ(2) + t91) * t75;
t124 = t76 * t29 + t73 * t36;
t72 = t78 ^ 2;
t122 = t71 + t72;
t115 = qJD(5) * t35;
t113 = qJD(1) * qJD(2);
t70 = t74 ^ 2;
t112 = t70 * t132;
t111 = t75 * t78 * t81;
t110 = 0.2e1 * qJD(2) * t71;
t109 = t74 * t119;
t108 = t79 * t120;
t106 = t77 * t119;
t105 = t73 * t115;
t103 = t74 * t116;
t63 = t77 * t117;
t100 = t122 * t81;
t57 = qJD(1) * t63;
t84 = -t78 * qJD(4) - t103;
t31 = t84 * qJD(1) + t57;
t87 = (-qJD(4) * t77 + qJD(2)) * t75;
t83 = qJD(1) * t87;
t10 = t76 * t31 + t73 * t83;
t99 = -t79 * t10 + t80 * t45;
t24 = -t74 * t101 + t77 * t48;
t98 = -t74 * t121 + t77 * t55;
t97 = qJ(2) * t113;
t96 = t35 ^ 2;
t94 = t78 * pkin(3) - t98;
t16 = t136 * qJD(5);
t14 = pkin(3) * t118 + qJD(4) - t24;
t1 = t39 * pkin(4) - t42 * pkin(6) + t14;
t3 = pkin(6) * t109 + t5;
t93 = t80 * t1 - t79 * t3;
t92 = -t79 * t1 - t80 * t3;
t90 = t80 * t10 + t79 * t45;
t4 = -t73 * t15 + t76 * t30;
t89 = -t73 * t29 + t76 * t36;
t46 = -qJD(1) * t103 + t57;
t88 = t46 * t74 - t134;
t50 = t76 * t127 - t78 * t73;
t86 = t75 * t128 - t79 * t50;
t27 = t75 * t129 + t50 * t80;
t85 = t76 * t129 + t77 * t80;
t20 = t75 * t108 + t80 * t42;
t82 = (-t76 * t128 + t77 * t79) * t35;
t66 = t71 * t97;
t52 = t63 - t103;
t44 = (t77 * t126 + t73 * t75) * qJD(1);
t41 = t77 * t104 - t102;
t34 = t63 + t84;
t22 = t27 * qJD(5);
t21 = t86 * qJD(5);
t17 = t20 * qJD(5);
t12 = t76 * t34 + t73 * t87;
t11 = t73 * t34 - t76 * t87;
t9 = t73 * t31 - t76 * t83;
t8 = pkin(6) * t130 + t124;
t7 = -pkin(4) * t130 - t89;
t6 = t49 * pkin(4) - t50 * pkin(6) + t94;
t2 = -pkin(4) * t109 - t4;
t13 = [0, 0, 0, 0, 0, 0.2e1 * t122 * t113, 0.2e1 * t72 * t97 + 0.2e1 * t66, t45 * t78 + (t74 * t110 + t51 * t78) * qJD(1), t46 * t78 + (t77 * t110 + t52 * t78) * qJD(1), ((t51 * t77 - t52 * t74) * qJD(1) - t88) * t75, qJD(2) * t133 + t46 * t123 - t24 * t51 + t25 * t52 - t45 * t98 + t66, t51 * t39 + t45 * t49 + (-qJD(1) * t11 - t9) * t130, t51 * t42 + t45 * t50 + (-qJD(1) * t12 - t10) * t130, -t10 * t49 + t11 * t42 - t12 * t39 + t9 * t50, t10 * t124 - t4 * t11 + t5 * t12 + t14 * t51 + t45 * t94 - t9 * t89, t16 * t27 + t20 * t21, t136 * t21 + t16 * t86 - t27 * t17 - t20 * t22, t16 * t49 + t21 * t35, -t17 * t49 - t22 * t35, 0, (-t79 * t12 + t80 * t51) * t35 + t99 * t49 - t11 * t136 + t7 * t17 - t9 * t86 + t2 * t22 + ((-t6 * t79 - t8 * t80) * t35 + t92 * t49) * qJD(5), -(t80 * t12 + t79 * t51) * t35 - t90 * t49 + t11 * t20 + t7 * t16 + t9 * t27 + t2 * t21 + (-(t6 * t80 - t79 * t8) * t35 - t93 * t49) * qJD(5); 0, 0, 0, 0, 0, -t100, -qJ(2) * t100, -t74 * t100, -t77 * t100, 0, (-t133 + (t24 * t74 - t25 * t77) * t78) * qJD(1) + t88, (-t39 * t78 + t41 * t75) * t120, (-t42 * t78 + t44 * t75) * t120, t44 * t39 - t41 * t42, t4 * t41 - t5 * t44 - t134 + (t10 * t76 - t14 * t118 + t73 * t9) * t74, 0, 0, 0, 0, 0, t17 * t131 - (t78 * t107 - t79 * t44) * t35 + t41 * t136 + qJD(5) * t82, t16 * t131 + (t78 * t108 + t80 * t44) * t35 - t41 * t20 + t85 * t115; 0, 0, 0, 0, 0, 0, 0, -t77 * t111, t74 * t111, (-t77 ^ 2 - t70) * t132, (t24 * t77 + t25 * t74 + qJD(2)) * t119, -t39 * t106 - t73 * t112, -t42 * t106 - t76 * t112, (-t39 * t76 + t42 * t73) * t109, t10 * t73 - t9 * t76 + (-t14 * t77 + (-t4 * t73 + t5 * t76) * t74) * t119, 0, 0, 0, 0, 0, -t80 * t105 - t76 * t17 + (-t131 * t136 - t85 * t35) * t119, t79 * t105 - t76 * t16 + (t20 * t131 + t82) * t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42 * t109, -t39 * t109, -t39 ^ 2 - t42 ^ 2, t5 * t39 + t4 * t42 + t45, 0, 0, 0, 0, 0, t136 * t42 - t79 * t96, -t42 * t20 - t80 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20 * t136, -t136 ^ 2 + t20 ^ 2, -t136 * t35 + t16, -t135 * t20, 0, t135 * t92 - t2 * t20 + t99, -t135 * t93 - t136 * t2 - t90;];
tauc_reg = t13;
