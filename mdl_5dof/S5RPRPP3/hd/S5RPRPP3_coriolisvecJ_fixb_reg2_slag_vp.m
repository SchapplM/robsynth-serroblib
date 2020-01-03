% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPP3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:54
% EndTime: 2019-12-31 18:12:56
% DurationCPUTime: 0.76s
% Computational Cost: add. (1198->186), mult. (3299->204), div. (0->0), fcn. (2325->4), ass. (0->107)
t132 = cos(qJ(3));
t85 = sin(pkin(7));
t86 = cos(pkin(7));
t88 = sin(qJ(3));
t65 = t132 * t85 + t88 * t86;
t61 = t65 * qJD(3);
t42 = qJD(1) * t61;
t113 = t132 * t86;
t107 = qJD(1) * t113;
t128 = t88 * t85;
t114 = qJD(1) * t128;
t56 = -t107 + t114;
t139 = t42 * qJ(5) + t56 * qJD(5);
t138 = t65 * qJD(1);
t134 = t138 ^ 2;
t135 = t56 ^ 2;
t137 = -t135 - t134;
t136 = -t135 + t134;
t109 = qJD(3) * t132;
t110 = qJD(2) * t132;
t126 = pkin(6) + qJ(2);
t69 = t126 * t85;
t70 = t126 * t86;
t22 = (qJD(2) * t85 + qJD(3) * t70) * t88 + t69 * t109 - t86 * t110;
t133 = t56 * pkin(4);
t106 = qJD(1) * t110;
t115 = qJD(1) * qJD(2);
t111 = t88 * t115;
t121 = qJD(3) * t88;
t66 = qJD(1) * t69;
t67 = qJD(1) * t70;
t16 = t85 * t106 + t67 * t109 + t86 * t111 - t66 * t121;
t37 = t132 * t69 + t88 * t70;
t131 = t16 * t37;
t80 = -t86 * pkin(2) - pkin(1);
t68 = t80 * qJD(1) + qJD(2);
t93 = -qJ(4) * t138 + t68;
t21 = t56 * pkin(3) + t93;
t130 = t21 * t138;
t129 = t56 * t138;
t127 = pkin(3) + qJ(5);
t124 = t85 ^ 2 + t86 ^ 2;
t123 = qJ(4) * t42;
t122 = t56 * qJ(4);
t60 = -t86 * t109 + t85 * t121;
t44 = qJD(3) * t60;
t46 = qJD(3) * t61;
t120 = t22 * qJD(3);
t38 = t132 * t70 - t88 * t69;
t23 = t65 * qJD(2) + t38 * qJD(3);
t119 = t23 * qJD(3);
t35 = t132 * t66 + t88 * t67;
t118 = -qJD(4) - t35;
t19 = -pkin(4) * t138 - t35;
t117 = qJD(4) - t19;
t36 = t132 * t67 - t88 * t66;
t20 = t36 - t133;
t116 = -qJD(5) - t20;
t112 = t124 * qJD(1) ^ 2;
t108 = -t86 * t106 + t66 * t109 + t85 * t111 + t67 * t121;
t75 = qJD(3) * t107;
t41 = qJD(3) * t114 - t75;
t105 = t42 * pkin(3) + t41 * qJ(4);
t5 = -t138 * t60 - t41 * t65;
t64 = -t113 + t128;
t6 = t42 * t64 + t56 * t61;
t104 = t60 * qJ(4) - t65 * qJD(4);
t103 = 0.2e1 * t124 * t115;
t101 = -t65 * qJ(4) + t80;
t100 = t41 * pkin(4) - t16;
t99 = -t42 * pkin(4) - t108;
t98 = -t35 * qJD(3) + t108;
t97 = t36 * qJD(3) - t16;
t10 = -qJD(4) * t138 + t105;
t7 = t127 * t56 + t93;
t95 = t138 * t7 - t100;
t94 = -t7 * t56 + t99;
t32 = -qJD(3) * qJ(4) - t36;
t92 = -t138 * t61 + t41 * t64 - t42 * t65 + t60 * t56;
t30 = 0.2e1 * t138 * qJD(3);
t91 = t138 * t23 + t16 * t65 + t22 * t56 - t37 * t41 - t38 * t42;
t89 = qJD(3) ^ 2;
t84 = qJD(3) * qJD(4);
t81 = 0.2e1 * t84;
t45 = qJD(3) * t56;
t40 = -t134 - t89;
t34 = t64 * pkin(3) + t101;
t33 = pkin(3) * t138 + t122;
t29 = t75 + (t56 - t114) * qJD(3);
t28 = t41 - t45;
t27 = -t75 + (t56 + t114) * qJD(3);
t26 = -qJD(3) * pkin(3) - t118;
t25 = -t64 * pkin(4) + t38;
t24 = t65 * pkin(4) + t37;
t18 = t127 * t64 + t101;
t17 = t61 * pkin(3) + t104;
t14 = t127 * t138 + t122;
t13 = qJD(5) - t32 - t133;
t12 = -t84 + t108;
t11 = -t127 * qJD(3) + t117;
t9 = -t60 * pkin(4) + t23;
t8 = -t61 * pkin(4) - t22;
t4 = -qJD(3) * qJD(5) - t100;
t3 = t84 + t99;
t2 = t64 * qJD(5) + t127 * t61 + t104;
t1 = t10 + t139;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, qJ(2) * t103, t5, t92, -t44, t6, -t46, 0, t80 * t42 + t68 * t61 - t119, -t80 * t41 - t68 * t60 + t120, t108 * t64 - t35 * t60 - t36 * t61 + t91, -t108 * t38 - t36 * t22 + t35 * t23 + t131, 0, t44, t46, t5, t92, t6, t12 * t64 - t26 * t60 + t32 * t61 + t91, -t10 * t64 - t17 * t56 - t21 * t61 - t34 * t42 + t119, -t10 * t65 - t138 * t17 + t21 * t60 + t34 * t41 - t120, t10 * t34 - t12 * t38 + t21 * t17 + t32 * t22 + t26 * t23 + t131, 0, t46, -t44, t6, -t92, t5, -t11 * t60 - t13 * t61 + t138 * t9 - t24 * t41 - t25 * t42 - t3 * t64 + t4 * t65 - t8 * t56, t8 * qJD(3) - t1 * t65 - t138 * t2 + t18 * t41 + t7 * t60, -t9 * qJD(3) + t1 * t64 + t18 * t42 + t2 * t56 + t7 * t61, t1 * t18 + t11 * t9 + t13 * t8 + t7 * t2 + t4 * t24 + t3 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112, -qJ(2) * t112, 0, 0, 0, 0, 0, 0, t30, -t27, t137, -t138 * t35 + t36 * t56, 0, 0, 0, 0, 0, 0, t137, -t30, t41 + t45, -t32 * t56 + (-qJD(4) - t26) * t138 + t105, 0, 0, 0, 0, 0, 0, t137, t27, t30, t13 * t56 + (-qJD(4) - t11) * t138 + t105 + t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, t136, t29, -t129, 0, 0, -t138 * t68 + t97, t68 * t56 + t98, 0, 0, 0, t28, 0, t129, t136, -t129, pkin(3) * t41 - t123 + (-t32 - t36) * t138 + (t26 + t118) * t56, t33 * t56 + t130 - t97, t138 * t33 - t21 * t56 + t81 - t98, -t16 * pkin(3) - t12 * qJ(4) + t118 * t32 - t21 * t33 - t26 * t36, 0, 0, t29, -t129, -t136, t129, -t123 + t127 * t41 + (t13 + t116) * t138 + (t11 - t117) * t56, -t19 * qJD(3) + t138 * t14 + t81 + t94, -t14 * t56 + (0.2e1 * qJD(5) + t20) * qJD(3) - t95, t3 * qJ(4) + t116 * t11 + t117 * t13 - t127 * t4 - t7 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t129, t40, t32 * qJD(3) + t130 + t16, 0, 0, 0, 0, 0, 0, t29, t40, t129, (-qJD(5) - t13) * qJD(3) + t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, -t135 - t89, t11 * qJD(3) + t84 + t94;];
tauc_reg = t15;
