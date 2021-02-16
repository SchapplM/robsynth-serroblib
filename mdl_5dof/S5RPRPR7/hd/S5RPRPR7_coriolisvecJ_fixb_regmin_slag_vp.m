% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:07
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:06:27
% EndTime: 2021-01-15 12:06:32
% DurationCPUTime: 1.00s
% Computational Cost: add. (1278->181), mult. (3154->261), div. (0->0), fcn. (2170->8), ass. (0->106)
t63 = sin(pkin(8)) * pkin(1) + pkin(6);
t119 = qJ(4) + t63;
t75 = sin(qJ(3));
t77 = cos(qJ(3));
t92 = t119 * qJD(1);
t37 = t75 * qJD(2) + t92 * t77;
t71 = sin(pkin(9));
t127 = t71 * t37;
t36 = t77 * qJD(2) - t92 * t75;
t107 = qJD(1) * qJD(4);
t140 = -t37 * qJD(3) - t75 * t107;
t115 = qJD(1) * t75;
t118 = cos(pkin(9));
t99 = t118 * t77;
t59 = qJD(1) * t99;
t46 = t71 * t115 - t59;
t44 = qJD(5) + t46;
t139 = -qJD(5) + t44;
t104 = pkin(3) * t115;
t100 = t118 * t75;
t54 = t71 * t77 + t100;
t117 = qJD(1) * t54;
t27 = t36 * qJD(3) + t77 * t107;
t3 = -t118 * t140 + t71 * t27;
t62 = t71 * pkin(3) + pkin(7);
t138 = (pkin(4) * t117 + t46 * pkin(7) + qJD(5) * t62 + t104) * t44 + t3;
t48 = t54 * qJD(3);
t41 = qJD(1) * t48;
t76 = cos(qJ(5));
t39 = t76 * t41;
t74 = sin(qJ(5));
t113 = qJD(5) * t74;
t126 = t71 * t75;
t85 = t99 - t126;
t51 = t85 * qJD(3);
t86 = t54 * t113 - t76 * t51;
t137 = t54 * t39 - t86 * t44;
t35 = t74 * qJD(3) + t117 * t76;
t108 = qJD(1) * qJD(3);
t103 = t75 * t108;
t58 = t71 * t103;
t42 = qJD(3) * t59 - t58;
t17 = t35 * qJD(5) + t74 * t42;
t93 = qJD(3) * t119;
t38 = t77 * qJD(4) - t75 * t93;
t83 = -t75 * qJD(4) - t77 * t93;
t13 = t118 * t38 + t71 * t83;
t65 = -cos(pkin(8)) * pkin(1) - pkin(2);
t55 = -t77 * pkin(3) + t65;
t116 = qJD(1) * t55;
t45 = qJD(4) + t116;
t14 = t46 * pkin(4) - pkin(7) * t117 + t45;
t22 = -pkin(4) * t85 - t54 * pkin(7) + t55;
t101 = t118 * t27;
t4 = t140 * t71 + t101;
t120 = qJD(3) * pkin(3);
t31 = t36 + t120;
t8 = t118 * t31 - t127;
t6 = -qJD(3) * pkin(4) - t8;
t52 = t119 * t77;
t25 = t118 * t52 - t119 * t126;
t91 = -t25 * t41 + t3 * t54;
t136 = -(qJD(5) * t22 + t13) * t44 + t6 * t51 + (qJD(5) * t14 + t4) * t85 + t91;
t135 = pkin(3) * t75;
t110 = t76 * qJD(3);
t16 = qJD(5) * t110 - t113 * t117 + t76 * t42;
t134 = -t16 * t85 + t35 * t48;
t133 = t16 * t74;
t132 = t22 * t41;
t33 = t117 * t74 - t110;
t131 = t33 * t44;
t130 = t35 * t44;
t129 = t35 * t117;
t128 = t117 * t33;
t125 = t74 * t41;
t78 = qJD(3) ^ 2;
t123 = t78 * t75;
t122 = t78 * t77;
t29 = t118 * t37;
t9 = t71 * t31 + t29;
t121 = t75 ^ 2 - t77 ^ 2;
t57 = qJD(1) * t65;
t114 = qJD(5) * t54;
t112 = t57 * qJD(1);
t105 = t75 * t120;
t97 = t44 * t76;
t96 = 0.2e1 * t117;
t7 = qJD(3) * pkin(7) + t9;
t1 = t76 * t14 - t74 * t7;
t2 = t74 * t14 + t76 * t7;
t90 = t17 * t85 - t48 * t33;
t88 = 0.2e1 * qJD(3) * t57;
t87 = t39 + (-t46 * t74 - t113) * t44;
t11 = t118 * t36 - t127;
t82 = -t62 * t41 + (t11 + t6) * t44;
t80 = (-t76 * t114 - t74 * t51) * t44 - t54 * t125;
t79 = qJD(1) ^ 2;
t64 = -t118 * pkin(3) - pkin(4);
t60 = pkin(3) * t103;
t24 = t119 * t100 + t71 * t52;
t20 = t48 * pkin(4) - t51 * pkin(7) + t105;
t18 = t41 * pkin(4) - t42 * pkin(7) + t60;
t15 = t76 * t18;
t12 = -t118 * t83 + t71 * t38;
t10 = t71 * t36 + t29;
t5 = [0, 0, 0, 0, 0.2e1 * t77 * t103, -0.2e1 * t121 * t108, t122, -t123, 0, -t63 * t122 + t75 * t88, t63 * t123 + t77 * t88, t55 * t41 + t45 * t48 + (-t12 + (-qJD(1) * t85 + t46) * t135) * qJD(3), t55 * t42 + t45 * t51 + (t96 * t135 - t13) * qJD(3), t117 * t12 - t13 * t46 + t24 * t42 + t4 * t85 - t9 * t48 - t8 * t51 + t91, -t8 * t12 + t9 * t13 + t3 * t24 + t4 * t25 + (t45 + t116) * t105, t16 * t76 * t54 - t86 * t35, (-t33 * t76 - t35 * t74) * t51 + (-t133 - t17 * t76 + (t33 * t74 - t35 * t76) * qJD(5)) * t54, t134 + t137, t80 + t90, -t41 * t85 + t44 * t48, t1 * t48 + t12 * t33 - t15 * t85 + t24 * t17 + (t20 * t44 + t132 + (-t25 * t44 + t6 * t54 + t7 * t85) * qJD(5)) * t76 + t136 * t74, t12 * t35 + t24 * t16 - t2 * t48 + (-(-qJD(5) * t25 + t20) * t44 - t132 + (-qJD(5) * t7 + t18) * t85 - t6 * t114) * t74 + t136 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t123, -t122, -t48 * qJD(3), -t51 * qJD(3), t117 * t48 - t54 * t41 - t42 * t85 - t51 * t46, -t3 * t85 + t4 * t54 - t8 * t48 + t9 * t51, 0, 0, 0, 0, 0, t80 - t90, t134 - t137; 0, 0, 0, 0, -t75 * t79 * t77, t121 * t79, 0, 0, 0, -t75 * t112, -t77 * t112, t10 * qJD(3) - t46 * t104 - t117 * t45 - t3, -t101 + t45 * t46 + (-pkin(3) * t117 + qJD(4) * t71) * t115 + (t11 + t127) * qJD(3), (-t10 + t9) * t117 + (t11 - t8) * t46 + (-t118 * t42 - t41 * t71) * pkin(3), t8 * t10 - t9 * t11 + (-t45 * t115 - t118 * t3 + t4 * t71) * pkin(3), t35 * t97 + t133, (t16 - t131) * t76 + (-t17 - t130) * t74, t44 * t97 + t125 - t129, t87 + t128, -t44 * t117, -t1 * t117 - t10 * t33 - t138 * t76 + t64 * t17 + t82 * t74, -t10 * t35 + t117 * t2 + t138 * t74 + t64 * t16 + t82 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96 * qJD(3), -t58 + (t59 - t46) * qJD(3), -t117 ^ 2 - t46 ^ 2, t117 * t8 + t9 * t46 + t60, 0, 0, 0, 0, 0, t87 - t128, -t44 ^ 2 * t76 - t125 - t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35 * t33, -t33 ^ 2 + t35 ^ 2, t16 + t131, t130 - t17, t41, t139 * t2 - t6 * t35 - t74 * t4 + t15, t139 * t1 - t74 * t18 + t6 * t33 - t76 * t4;];
tauc_reg = t5;
