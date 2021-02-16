% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% tauc_reg [5x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:27
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:27:19
% EndTime: 2021-01-15 12:27:24
% DurationCPUTime: 0.76s
% Computational Cost: add. (1285->158), mult. (2737->215), div. (0->0), fcn. (1695->4), ass. (0->111)
t78 = sin(qJ(3));
t118 = qJD(1) * t78;
t77 = sin(qJ(4));
t107 = t77 * t118;
t72 = qJD(3) + qJD(4);
t125 = t72 * t107;
t79 = cos(qJ(4));
t80 = cos(qJ(3));
t129 = t79 * t80;
t97 = t72 * t129;
t20 = qJD(1) * t97 - t125;
t114 = qJD(4) * t77;
t117 = qJD(1) * t80;
t81 = -pkin(1) - pkin(6);
t60 = t81 * qJD(1) + qJD(2);
t37 = -pkin(7) * t117 + t80 * t60;
t32 = qJD(3) * pkin(3) + t37;
t116 = qJD(3) * t78;
t99 = pkin(7) * qJD(1) - t60;
t33 = t99 * t116;
t115 = qJD(3) * t80;
t34 = t99 * t115;
t36 = -pkin(7) * t118 + t78 * t60;
t92 = -(qJD(4) * t32 - t34) * t79 + t36 * t114 - t77 * t33;
t140 = -t20 * qJ(5) - t92;
t29 = t77 * t36;
t104 = t79 * t32 - t29;
t106 = t79 * t117;
t45 = t106 - t107;
t38 = t45 * qJ(5);
t6 = t104 - t38;
t50 = t77 * t80 + t79 * t78;
t43 = t50 * qJD(1);
t139 = t45 ^ 2;
t73 = qJD(1) * qJD(2);
t138 = 0.2e1 * t73;
t5 = t72 * pkin(4) + t6;
t137 = t5 - t6;
t136 = pkin(3) * t72;
t135 = pkin(7) - t81;
t23 = -t78 * t114 - t77 * t116 + t97;
t134 = t23 * t72;
t133 = t43 * t45;
t131 = t45 * t72;
t56 = pkin(3) * t118 + qJD(1) * qJ(2);
t130 = t56 * t45;
t30 = t79 * t36;
t82 = qJD(3) ^ 2;
t128 = t82 * t78;
t127 = t82 * t80;
t126 = t79 * t37 - t29;
t110 = qJD(1) * qJD(3);
t105 = t80 * t110;
t53 = pkin(3) * t105 + t73;
t124 = t78 ^ 2 - t80 ^ 2;
t83 = qJD(1) ^ 2;
t123 = -t82 - t83;
t22 = t72 * t50;
t19 = t22 * qJD(1);
t122 = t19 * qJ(5);
t120 = t43 * qJ(5);
t119 = t83 * qJ(2);
t66 = t78 * pkin(3) + qJ(2);
t113 = qJD(4) * t79;
t100 = -t43 * pkin(4) - qJD(5);
t24 = -t100 + t56;
t112 = qJD(5) + t24;
t61 = pkin(3) * t115 + qJD(2);
t111 = qJ(2) * qJD(3);
t109 = 0.2e1 * qJD(1);
t108 = pkin(3) * t117;
t55 = t135 * t80;
t103 = t79 * t33 + t77 * t34;
t101 = -t77 * t37 - t30;
t14 = t20 * pkin(4) + t53;
t51 = -t77 * t78 + t129;
t96 = -t51 * t19 - t45 * t22;
t95 = -t77 * t32 - t30;
t54 = t135 * t78;
t94 = t79 * t54 + t77 * t55;
t48 = t135 * t116;
t49 = qJD(3) * t55;
t93 = -t55 * t113 + t54 * t114 + t77 * t48 - t79 * t49;
t1 = -t43 * qJD(5) + t140;
t90 = t95 * qJD(4) + t103;
t87 = t90 + t122;
t2 = -t45 * qJD(5) + t87;
t7 = -t95 - t120;
t91 = t1 * t50 + t2 * t51 - t5 * t22 + t7 * t23;
t89 = t94 * qJD(4) + t79 * t48 + t77 * t49;
t88 = t56 * t43 + t92;
t86 = (-t30 + (-t32 - t136) * t77) * qJD(4) + t103;
t85 = t112 * t43 - t140;
t67 = t79 * pkin(3) + pkin(4);
t57 = t113 * t136;
t42 = t43 ^ 2;
t35 = t50 * pkin(4) + t66;
t27 = t45 * pkin(4) + t108;
t21 = t72 * t22;
t18 = t23 * pkin(4) + t61;
t17 = -t50 * qJ(5) - t94;
t16 = -t51 * qJ(5) + t77 * t54 - t79 * t55;
t15 = -t42 + t139;
t13 = -qJD(1) * t43 - t21;
t12 = -qJD(1) * t45 - t134;
t11 = -t72 * t106 + t125 + t131;
t9 = -t38 + t126;
t8 = t101 + t120;
t4 = t22 * qJ(5) - t51 * qJD(5) + t89;
t3 = -t23 * qJ(5) - t50 * qJD(5) + t93;
t10 = [0, 0, 0, 0, t138, qJ(2) * t138, -0.2e1 * t78 * t105, 0.2e1 * t124 * t110, -t128, -t127, 0, -t81 * t128 + (qJD(2) * t78 + t80 * t111) * t109, -t81 * t127 + (qJD(2) * t80 - t78 * t111) * t109, t96, t50 * t19 - t20 * t51 + t43 * t22 - t23 * t45, -t21, -t134, 0, t66 * t20 + t56 * t23 + t61 * t43 + t53 * t50 + t72 * t89, -t66 * t19 - t56 * t22 + t61 * t45 + t53 * t51 - t72 * t93, t14 * t50 + t18 * t43 + t35 * t20 + t24 * t23 + t4 * t72, t14 * t51 + t18 * t45 - t35 * t19 - t24 * t22 - t3 * t72, t16 * t19 - t17 * t20 - t3 * t43 - t4 * t45 - t91, t1 * t17 + t14 * t35 + t2 * t16 + t24 * t18 + t7 * t3 + t5 * t4; 0, 0, 0, 0, -t83, -t119, 0, 0, 0, 0, 0, t123 * t78, t123 * t80, 0, 0, 0, 0, 0, t13, t12, t13, t12, -t50 * t20 - t23 * t43 - t96, -t24 * qJD(1) + t91; 0, 0, 0, 0, 0, 0, t80 * t83 * t78, -t124 * t83, 0, 0, 0, -t80 * t119, t78 * t119, t133, t15, 0, t11, 0, -t101 * t72 - t43 * t108 - t130 + t86, -t45 * t108 + t126 * t72 - t57 + t88, -t112 * t45 - t27 * t43 - t8 * t72 + t122 + t86, -t27 * t45 + t9 * t72 - t57 + t85, t67 * t19 + (t7 + t8) * t45 + (-t5 + t9) * t43 + (-t20 * t77 + (-t43 * t79 + t45 * t77) * qJD(4)) * pkin(3), t2 * t67 - t24 * t27 - t5 * t8 - t7 * t9 + (t1 * t77 + (-t5 * t77 + t7 * t79) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133, t15, 0, t11, 0, -t72 * t95 - t130 + t90, t104 * t72 + t88, t7 * t72 + (t100 - t24) * t45 + t87, -t139 * pkin(4) + t6 * t72 + t85, t19 * pkin(4) - t137 * t43, t137 * t7 + (-t24 * t45 + t2) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20 + t131, -0.2e1 * t43 * t72, -t42 - t139, t7 * t43 + t5 * t45 + t14;];
tauc_reg = t10;
