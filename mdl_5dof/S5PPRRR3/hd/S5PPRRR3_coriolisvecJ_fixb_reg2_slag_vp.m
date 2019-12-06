% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPRRR3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:59
% EndTime: 2019-12-05 15:17:05
% DurationCPUTime: 1.15s
% Computational Cost: add. (1384->190), mult. (3592->288), div. (0->0), fcn. (2711->8), ass. (0->127)
t76 = sin(qJ(5));
t77 = sin(qJ(4));
t79 = cos(qJ(5));
t80 = cos(qJ(4));
t56 = t76 * t80 + t79 * t77;
t71 = qJD(4) + qJD(5);
t147 = t71 * t56;
t28 = t147 * qJD(3);
t74 = sin(pkin(9));
t122 = qJD(1) * t74;
t78 = sin(qJ(3));
t81 = cos(qJ(3));
t143 = t81 * qJD(2) - t78 * t122;
t47 = t143 * qJD(3);
t113 = qJD(3) * qJD(4);
t146 = -0.2e1 * t113;
t134 = t79 * t80;
t136 = t76 * t77;
t55 = -t134 + t136;
t43 = t55 * t78;
t75 = cos(pkin(9));
t121 = qJD(1) * t75;
t105 = t80 * t121;
t54 = t78 * qJD(2) + t81 * t122;
t45 = qJD(3) * pkin(6) + t54;
t98 = pkin(7) * qJD(3) + t45;
t29 = -t98 * t77 - t105;
t41 = t80 * t47;
t14 = t29 * qJD(4) + t41;
t123 = qJD(4) * pkin(4);
t26 = t29 + t123;
t145 = (qJD(5) * t26 + t14) * t79;
t82 = qJD(4) ^ 2;
t83 = qJD(3) ^ 2;
t144 = (t82 + t83) * t78;
t142 = -pkin(7) - pkin(6);
t70 = -t80 * pkin(4) - pkin(3);
t40 = t70 * qJD(3) - t143;
t52 = t56 * qJD(3);
t141 = t40 * t52;
t46 = t54 * qJD(3);
t140 = t46 * t78;
t117 = qJD(3) * t80;
t106 = t79 * t117;
t119 = qJD(3) * t77;
t107 = t76 * t119;
t50 = -t106 + t107;
t139 = t52 * t50;
t138 = t74 * t81;
t64 = t77 * t121;
t30 = t98 * t80 - t64;
t137 = t76 * t30;
t135 = t79 * t30;
t133 = t80 * t45;
t132 = t83 * t78;
t131 = t83 * t81;
t59 = t142 * t77;
t60 = t142 * t80;
t37 = t79 * t59 + t76 * t60;
t103 = qJD(4) * t142;
t57 = t77 * t103;
t58 = t80 * t103;
t130 = t37 * qJD(5) + t143 * t55 + t79 * t57 + t76 * t58;
t38 = t76 * t59 - t79 * t60;
t129 = -t38 * qJD(5) + t143 * t56 - t76 * t57 + t79 * t58;
t102 = t80 * t113;
t128 = -qJD(5) * t106 - t79 * t102;
t72 = t77 ^ 2;
t73 = t80 ^ 2;
t127 = t72 - t73;
t126 = t72 + t73;
t124 = qJD(3) * pkin(3);
t44 = -t143 - t124;
t120 = qJD(3) * t44;
t118 = qJD(3) * t78;
t116 = qJD(3) * t81;
t115 = qJD(4) * t80;
t112 = t77 * t83 * t80;
t111 = t74 * t131;
t110 = pkin(4) * t119;
t109 = t77 * t123;
t108 = t74 * t118;
t104 = -pkin(4) * t71 - t26;
t61 = qJD(4) * t64;
t99 = -t77 * t47 + t61;
t15 = -t98 * t115 + t99;
t101 = -t76 * t14 + t79 * t15;
t100 = -qJD(5) * t137 + t76 * t15;
t96 = t77 * t108;
t95 = t80 * t108;
t94 = t81 * t146;
t93 = t77 * t102;
t92 = -t54 + t109;
t91 = t71 * t136;
t6 = t76 * t26 + t135;
t35 = -t77 * t45 - t105;
t36 = -t64 + t133;
t90 = t35 * t77 - t36 * t80;
t48 = -t77 * t138 - t75 * t80;
t49 = t80 * t138 - t75 * t77;
t21 = t79 * t48 - t76 * t49;
t22 = t76 * t48 + t79 * t49;
t89 = t40 * t50 - t100;
t88 = pkin(6) * t82;
t87 = qJD(4) * (t44 + t143 - t124);
t2 = -t6 * qJD(5) + t101;
t18 = t35 * qJD(4) + t41;
t19 = -t45 * t115 + t99;
t85 = t18 * t80 - t19 * t77 + (-t35 * t80 - t36 * t77) * qJD(4);
t42 = t56 * t78;
t39 = (t54 + t109) * qJD(3);
t34 = t48 * qJD(4) - t95;
t33 = -t49 * qJD(4) + t96;
t31 = -qJD(5) * t134 - t79 * t115 + t91;
t27 = qJD(3) * t91 + t128;
t20 = -t50 ^ 2 + t52 ^ 2;
t17 = t52 * t71 - t28;
t16 = -t128 + (-t107 + t50) * t71;
t10 = -t56 * t116 + t71 * t43;
t9 = -t55 * t116 - t147 * t78;
t8 = t79 * t29 - t137;
t7 = -t76 * t29 - t135;
t5 = t79 * t26 - t137;
t4 = -t22 * qJD(5) + t79 * t33 - t76 * t34;
t3 = t21 * qJD(5) + t76 * t33 + t79 * t34;
t1 = t100 + t145;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111, t74 * t132, 0, (t140 + t47 * t81 + (-t143 * t81 - t54 * t78) * qJD(3)) * t74, 0, 0, 0, 0, 0, 0, -t80 * t111 + (t33 + t96) * qJD(4), t77 * t111 + (-t34 + t95) * qJD(4), (-t33 * t77 + t34 * t80 + (-t48 * t80 - t49 * t77) * qJD(4)) * qJD(3), t18 * t49 + t19 * t48 + t35 * t33 + t36 * t34 + (t44 * t116 + t140) * t74, 0, 0, 0, 0, 0, 0, t4 * t71 + (t50 * t116 + t28 * t78) * t74, -t3 * t71 + (t52 * t116 - t27 * t78) * t74, t21 * t27 - t22 * t28 - t3 * t50 - t4 * t52, t1 * t22 + t2 * t21 + t6 * t3 + t5 * t4 + (t40 * t116 + t39 * t78) * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t132, -t131, 0, -t46 * t81 + t47 * t78 + (-t143 * t78 + t54 * t81) * qJD(3), 0, 0, 0, 0, 0, 0, -t80 * t144 + t77 * t94, t77 * t144 + t80 * t94, t126 * t131, (-t90 * qJD(3) - t46) * t81 + (t85 + t120) * t78, 0, 0, 0, 0, 0, 0, t10 * t71 + t50 * t118 - t81 * t28, t52 * t118 + t81 * t27 - t9 * t71, -t10 * t52 - t42 * t27 + t43 * t28 - t9 * t50, -t1 * t43 + t5 * t10 + t118 * t40 - t2 * t42 - t39 * t81 + t6 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t93, t127 * t146, t82 * t80, -0.2e1 * t93, -t82 * t77, 0, t77 * t87 - t88 * t80, t88 * t77 + t80 * t87, -t126 * t47 + t85, -t46 * pkin(3) + t85 * pkin(6) + t143 * t90 - t44 * t54, -t27 * t56 - t52 * t31, -t147 * t52 + t27 * t55 - t56 * t28 + t31 * t50, -t31 * t71, t147 * t50 + t28 * t55, -t147 * t71, 0, t129 * t71 + t147 * t40 + t70 * t28 + t39 * t55 + t92 * t50, -t130 * t71 - t70 * t27 - t40 * t31 + t39 * t56 + t92 * t52, -t1 * t55 - t129 * t52 - t130 * t50 - t147 * t6 - t2 * t56 + t37 * t27 - t38 * t28 + t5 * t31, t1 * t38 + t129 * t5 + t130 * t6 + t2 * t37 + t39 * t70 + t40 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112, t127 * t83, 0, t112, 0, 0, t61 + (-t47 - t120) * t77 + (t36 - t133) * qJD(4), -t44 * t117 - t41, 0, 0, t139, t20, t16, -t139, t17, 0, -t50 * t110 - t141 - t7 * t71 + (t104 * t76 - t135) * qJD(5) + t101, -t52 * t110 + t8 * t71 + (t104 * qJD(5) - t14) * t79 + t89, (t6 + t7) * t52 + (-t5 + t8) * t50 + (t27 * t79 - t28 * t76 + (-t50 * t79 + t52 * t76) * qJD(5)) * pkin(4), -t5 * t7 - t6 * t8 + (-t40 * t119 + t1 * t76 + t2 * t79 + (-t5 * t76 + t6 * t79) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, t20, t16, -t139, t17, 0, t6 * t71 - t141 + t2, t5 * t71 - t145 + t89, 0, 0;];
tauc_reg = t11;
