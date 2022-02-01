% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPPR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:51:43
% EndTime: 2022-01-20 09:51:46
% DurationCPUTime: 0.75s
% Computational Cost: add. (1541->138), mult. (3008->191), div. (0->0), fcn. (1974->8), ass. (0->107)
t91 = sin(pkin(9));
t93 = cos(pkin(9));
t121 = t91 ^ 2 + t93 ^ 2;
t120 = pkin(1) * qJD(1);
t98 = cos(qJ(2));
t116 = t98 * t120;
t96 = sin(qJ(2));
t117 = t96 * t120;
t92 = sin(pkin(8));
t76 = t92 * t117;
t94 = cos(pkin(8));
t60 = t94 * t116 - t76;
t49 = t60 * qJD(2);
t90 = qJD(1) + qJD(2);
t40 = qJD(4) * t90 + t49;
t142 = t40 * t121;
t139 = t60 - qJD(4);
t97 = cos(qJ(5));
t129 = t97 * t93;
t95 = sin(qJ(5));
t130 = t95 * t91;
t67 = -t129 + t130;
t141 = t40 * t67;
t68 = t91 * t97 + t93 * t95;
t52 = t68 * t90;
t140 = t121 * t90;
t138 = t52 ^ 2;
t137 = t93 * pkin(4);
t86 = t93 * pkin(7);
t136 = t94 * pkin(2);
t118 = t90 * t130;
t50 = -t129 * t90 + t118;
t135 = t52 * t50;
t131 = t94 * t96;
t102 = t92 * t98 + t131;
t101 = pkin(1) * t102;
t58 = qJD(1) * t101;
t134 = t58 * t90;
t59 = qJD(2) * t101;
t133 = t59 * t90;
t132 = t90 * t91;
t81 = pkin(2) * t92 + qJ(4);
t64 = (-pkin(7) - t81) * t91;
t65 = t81 * t93 + t86;
t27 = t64 * t97 - t65 * t95;
t128 = qJD(5) * t27 + t139 * t67;
t28 = t64 * t95 + t65 * t97;
t127 = -qJD(5) * t28 + t139 * t68;
t69 = pkin(2) * t90 + t116;
t41 = t94 * t69 - t76;
t107 = qJD(4) - t41;
t112 = -pkin(3) - t137;
t26 = t112 * t90 + t107;
t48 = qJD(1) * t59;
t63 = t68 * qJD(5);
t126 = t26 * t63 + t48 * t67;
t45 = t90 * t63;
t114 = qJD(5) * t129;
t115 = qJD(5) * t130;
t62 = -t114 + t115;
t125 = -t45 * t68 + t50 * t62;
t124 = -t26 * t62 + t48 * t68;
t42 = t117 * t94 + t69 * t92;
t37 = qJ(4) * t90 + t42;
t25 = t91 * qJD(3) + t37 * t93;
t83 = pkin(1) * t98 + pkin(2);
t122 = pkin(1) * t131 + t92 * t83;
t119 = pkin(1) * qJD(2);
t79 = t92 * t96 * pkin(1);
t85 = t93 * qJD(3);
t22 = t85 + (-pkin(7) * t90 - t37) * t91;
t23 = t86 * t90 + t25;
t10 = t22 * t97 - t23 * t95;
t11 = t22 * t95 + t23 * t97;
t3 = t10 * qJD(5) - t141;
t100 = t68 * t40;
t4 = -t11 * qJD(5) - t100;
t113 = t10 * t62 - t11 * t63 - t3 * t67 - t4 * t68;
t110 = t83 * t94 - t79;
t108 = -pkin(3) - t110;
t106 = (-qJD(2) + t90) * t120;
t105 = (-qJD(1) - t90) * t119;
t104 = (-t37 * t91 + t85) * t91 - t25 * t93;
t55 = qJ(4) + t122;
t38 = (-pkin(7) - t55) * t91;
t39 = t55 * t93 + t86;
t14 = t38 * t97 - t39 * t95;
t15 = t38 * t95 + t39 * t97;
t70 = t90 * t114;
t44 = t115 * t90 - t70;
t103 = -t44 * t67 + t52 * t63;
t61 = t94 * t98 * t119 - qJD(2) * t79;
t99 = t102 * qJD(1) * t119;
t71 = t112 - t136;
t57 = t63 * qJD(5);
t56 = t62 * qJD(5);
t54 = qJD(4) + t61;
t47 = t50 ^ 2;
t46 = t108 - t137;
t43 = t48 * t91;
t34 = -pkin(3) * t90 + t107;
t13 = t45 * t67 + t50 * t63;
t12 = -t44 * t68 - t52 * t62;
t9 = -qJD(5) * t15 - t54 * t68;
t8 = qJD(5) * t14 - t54 * t67;
t5 = -t103 + t125;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96 * t105, t98 * t105, 0, 0, 0, 0, 0, 0, 0, 0, -t99 - t133, -t61 * t90 - t49, 0, -t110 * t48 + t122 * t49 - t41 * t59 + t42 * t61, 0, 0, 0, 0, 0, 0, (-t48 - t133) * t93, t132 * t59 + t43, t140 * t54 + t142, -t104 * t54 + t108 * t48 + t142 * t55 + t34 * t59, t12, t5, -t56, t13, -t57, 0, qJD(5) * t9 + t45 * t46 + t50 * t59 + t126, -qJD(5) * t8 - t44 * t46 + t52 * t59 + t124, t14 * t44 - t15 * t45 - t50 * t8 - t52 * t9 + t113, t10 * t9 + t11 * t8 + t14 * t4 + t15 * t3 + t26 * t59 + t46 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96 * t106, t98 * t106, 0, 0, 0, 0, 0, 0, 0, 0, -t99 + t134, t60 * t90 - t49, 0, t41 * t58 - t42 * t60 + (-t48 * t94 + t49 * t92) * pkin(2), 0, 0, 0, 0, 0, 0, (-t48 + t134) * t93, -t132 * t58 + t43, -t139 * t140 + t142, t48 * (-pkin(3) - t136) - t34 * t58 + t81 * t142 + t139 * t104, t12, t5, -t56, t13, -t57, 0, qJD(5) * t127 + t71 * t45 - t58 * t50 + t126, -qJD(5) * t128 - t71 * t44 - t58 * t52 + t124, -t127 * t52 - t128 * t50 + t27 * t44 - t28 * t45 + t113, t10 * t127 + t11 * t128 - t26 * t58 + t4 * t27 + t3 * t28 + t48 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, t56, t103 + t125, -t10 * t63 - t11 * t62 + t3 * t68 - t4 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121 * t90 ^ 2, t104 * t90 + t48, 0, 0, 0, 0, 0, 0, 0.2e1 * t52 * qJD(5), t70 + (-t50 - t118) * qJD(5), -t47 - t138, t10 * t52 + t11 * t50 + t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, -t47 + t138, t70 + (t50 - t118) * qJD(5), -t135, 0, 0, -t26 * t52 - t100, t26 * t50 + t141, 0, 0;];
tauc_reg = t1;
