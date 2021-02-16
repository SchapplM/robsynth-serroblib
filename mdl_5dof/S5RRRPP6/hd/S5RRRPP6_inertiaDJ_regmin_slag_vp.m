% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:38
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPP6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:37:06
% EndTime: 2021-01-15 22:37:12
% DurationCPUTime: 1.18s
% Computational Cost: add. (1552->190), mult. (3907->350), div. (0->0), fcn. (3090->6), ass. (0->104)
t81 = cos(qJ(2));
t114 = t81 * qJD(2);
t78 = sin(qJ(3));
t104 = t78 * t114;
t80 = cos(qJ(3));
t118 = qJD(3) * t80;
t79 = sin(qJ(2));
t138 = t79 * t118 + t104;
t137 = -0.4e1 * t79;
t91 = -t81 * pkin(2) - t79 * pkin(7);
t56 = -pkin(1) + t91;
t129 = t80 * t81;
t64 = pkin(6) * t129;
t125 = t78 * t56 + t64;
t119 = qJD(3) * t78;
t77 = sin(pkin(8));
t120 = cos(pkin(8));
t98 = t120 * t80;
t41 = qJD(3) * t98 - t77 * t119;
t74 = t79 ^ 2;
t95 = (-t81 ^ 2 + t74) * qJD(2);
t75 = t80 ^ 2;
t124 = t78 ^ 2 - t75;
t96 = t124 * qJD(3);
t136 = 2 * qJD(5);
t90 = pkin(2) * t79 - pkin(7) * t81;
t50 = t90 * qJD(2);
t127 = -t56 * t118 - t78 * t50;
t130 = t79 * t80;
t12 = (-pkin(6) * qJD(2) - qJ(4) * qJD(3)) * t130 + (-qJD(4) * t79 + (-pkin(6) * qJD(3) - qJ(4) * qJD(2)) * t81) * t78 - t127;
t115 = t80 * qJD(4);
t121 = qJ(4) * t80;
t122 = qJ(4) * t79;
t116 = t79 * qJD(2);
t105 = t78 * t116;
t126 = pkin(6) * t105 + t80 * t50;
t135 = pkin(3) * t79;
t8 = -t79 * t115 + (-t81 * t121 + t135) * qJD(2) + (-t64 + (-t56 + t122) * t78) * qJD(3) + t126;
t4 = t120 * t12 + t77 * t8;
t134 = pkin(6) * t78;
t133 = t77 * pkin(3);
t132 = t77 * t78;
t131 = t78 * t79;
t128 = qJ(4) + pkin(7);
t48 = t80 * t56;
t30 = -t79 * t121 + t48 + (-pkin(3) - t134) * t81;
t32 = -t78 * t122 + t125;
t16 = t120 * t32 + t77 * t30;
t51 = pkin(3) * t131 + t79 * pkin(6);
t117 = qJD(3) * t81;
t113 = t81 * qJD(5);
t112 = qJ(5) * t116 + t4;
t111 = -0.2e1 * pkin(1) * qJD(2);
t110 = -0.2e1 * pkin(2) * qJD(3);
t70 = pkin(6) * t114;
t35 = t138 * pkin(3) + t70;
t71 = pkin(3) * t119;
t109 = t78 * t117;
t107 = t80 * t117;
t103 = t78 * t118;
t102 = t79 * t114;
t101 = t80 * t114;
t69 = -t80 * pkin(3) - pkin(2);
t3 = -t77 * t12 + t120 * t8;
t97 = qJD(3) * t128;
t39 = -t78 * t97 + t115;
t83 = -t78 * qJD(4) - t80 * t97;
t26 = -t120 * t83 + t77 * t39;
t27 = t120 * t39 + t77 * t83;
t57 = t128 * t80;
t99 = t120 * t78;
t33 = t128 * t99 + t77 * t57;
t34 = t120 * t57 - t128 * t132;
t100 = t33 * t26 + t34 * t27;
t94 = 0.2e1 * t102;
t93 = qJD(2) * t120;
t92 = t78 * t101;
t89 = t81 * t93;
t87 = -t33 * t116 + t26 * t81;
t86 = t34 * t116 - t27 * t81;
t15 = t120 * t30 - t77 * t32;
t46 = t77 * t80 + t99;
t24 = -t77 * t101 - t41 * t79 - t78 * t89;
t40 = t46 * qJD(3);
t25 = t77 * t104 + t79 * t40 - t80 * t89;
t37 = t46 * t79;
t38 = -t77 * t131 + t79 * t98;
t85 = t34 * t24 - t33 * t25 + t26 * t38 - t27 * t37;
t84 = t80 * t116 + t109;
t45 = -t98 + t132;
t82 = 0.2e1 * t26 * t46 - 0.2e1 * t27 * t45 + 0.2e1 * t33 * t41 - 0.2e1 * t34 * t40;
t67 = -t120 * pkin(3) - pkin(4);
t65 = qJ(5) + t133;
t28 = t45 * pkin(4) - t46 * qJ(5) + t69;
t22 = -t125 * qJD(3) + t126;
t21 = t84 * pkin(6) + t127;
t20 = t37 * pkin(4) - t38 * qJ(5) + t51;
t18 = t40 * pkin(4) - t41 * qJ(5) - t46 * qJD(5) + t71;
t11 = t81 * pkin(4) - t15;
t10 = -t81 * qJ(5) + t16;
t5 = -t24 * pkin(4) + t25 * qJ(5) - t38 * qJD(5) + t35;
t2 = -pkin(4) * t116 - t3;
t1 = t112 - t113;
t6 = [0, 0, 0, t94, -0.2e1 * t95, 0, 0, 0, t79 * t111, t81 * t111, 0.2e1 * t75 * t102 - 0.2e1 * t74 * t103, t92 * t137 + 0.2e1 * t74 * t96, 0.2e1 * t79 * t109 + 0.2e1 * t80 * t95, 0.2e1 * t79 * t107 - 0.2e1 * t78 * t95, -0.2e1 * t102, 0.2e1 * t48 * t116 - 0.2e1 * t22 * t81 + 0.2e1 * (t78 * t102 + t74 * t118) * pkin(6), -0.2e1 * t21 * t81 - 0.2e1 * t125 * t116 + 0.2e1 * (-t74 * t119 + t80 * t94) * pkin(6), 0.2e1 * t15 * t116 - 0.2e1 * t51 * t24 - 0.2e1 * t3 * t81 + 0.2e1 * t35 * t37, -0.2e1 * t16 * t116 - 0.2e1 * t51 * t25 + 0.2e1 * t35 * t38 + 0.2e1 * t4 * t81, 0.2e1 * t15 * t25 + 0.2e1 * t16 * t24 - 0.2e1 * t3 * t38 - 0.2e1 * t4 * t37, 0.2e1 * t15 * t3 + 0.2e1 * t16 * t4 + 0.2e1 * t51 * t35, -0.2e1 * t11 * t116 + 0.2e1 * t2 * t81 - 0.2e1 * t20 * t24 + 0.2e1 * t5 * t37, -0.2e1 * t1 * t37 + 0.2e1 * t10 * t24 - 0.2e1 * t11 * t25 + 0.2e1 * t2 * t38, -0.2e1 * t1 * t81 + 0.2e1 * t10 * t116 + 0.2e1 * t20 * t25 - 0.2e1 * t5 * t38, 0.2e1 * t10 * t1 + 0.2e1 * t11 * t2 + 0.2e1 * t20 * t5; 0, 0, 0, 0, 0, t114, -t116, 0, -t70, pkin(6) * t116, -t79 * t96 + t92, t103 * t137 - t124 * t114, t105 - t107, t84, 0, (pkin(7) * t129 + (-pkin(2) * t80 + t134) * t79) * qJD(3) + (t91 * t78 - t64) * qJD(2), (pkin(6) * t130 + t90 * t78) * qJD(3) + (t81 * t134 + t91 * t80) * qJD(2), -t69 * t24 + t35 * t45 + t37 * t71 + t51 * t40 + t87, -t69 * t25 + t35 * t46 + t38 * t71 + t51 * t41 - t86, -t15 * t41 - t16 * t40 - t3 * t46 - t4 * t45 + t85, -t15 * t26 + t16 * t27 - t3 * t33 + t4 * t34 + t35 * t69 + t51 * t71, t18 * t37 + t20 * t40 - t28 * t24 + t5 * t45 + t87, -t1 * t45 - t10 * t40 + t11 * t41 + t2 * t46 + t85, -t18 * t38 - t20 * t41 + t28 * t25 - t5 * t46 + t86, t1 * t34 + t10 * t27 + t11 * t26 + t20 * t18 + t2 * t33 + t5 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t103, -0.2e1 * t96, 0, 0, 0, t78 * t110, t80 * t110, 0.2e1 * t69 * t40 + 0.2e1 * t45 * t71, 0.2e1 * t69 * t41 + 0.2e1 * t46 * t71, t82, 0.2e1 * t69 * t71 + 0.2e1 * t100, 0.2e1 * t18 * t45 + 0.2e1 * t28 * t40, t82, -0.2e1 * t18 * t46 - 0.2e1 * t28 * t41, 0.2e1 * t28 * t18 + 0.2e1 * t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79 * t119 + t101, -t138, t116, t22, t21, t93 * t135 + t3, -t116 * t133 - t4, (t120 * t25 + t24 * t77) * pkin(3), (t120 * t3 + t4 * t77) * pkin(3), (pkin(4) - t67) * t116 + t3, -qJD(5) * t37 + t65 * t24 - t67 * t25, t65 * t116 + t112 - 0.2e1 * t113, t10 * qJD(5) + t1 * t65 + t2 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, -t119, 0, -pkin(7) * t118, pkin(7) * t119, -t26, -t27, (-t120 * t41 - t40 * t77) * pkin(3), (-t120 * t26 + t27 * t77) * pkin(3), -t26, -qJD(5) * t45 - t65 * t40 + t67 * t41, t27, t34 * qJD(5) + t26 * t67 + t27 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, t65 * t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t25, 0, t35, -t24, 0, t25, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t41, 0, t71, t40, 0, -t41, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116, -t25, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
