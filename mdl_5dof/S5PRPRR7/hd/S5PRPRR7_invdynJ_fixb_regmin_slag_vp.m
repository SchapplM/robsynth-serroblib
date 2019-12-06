% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% tau_reg [5x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRR7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR7_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:00:48
% EndTime: 2019-12-05 16:00:52
% DurationCPUTime: 0.93s
% Computational Cost: add. (638->183), mult. (1237->241), div. (0->0), fcn. (878->10), ass. (0->119)
t69 = cos(qJ(2));
t56 = g(3) * t69;
t66 = sin(qJ(2));
t62 = sin(pkin(8));
t63 = cos(pkin(8));
t98 = g(1) * t63 + g(2) * t62;
t148 = t98 * t66 - t56;
t64 = sin(qJ(5));
t65 = sin(qJ(4));
t67 = cos(qJ(5));
t68 = cos(qJ(4));
t34 = t64 * t68 + t67 * t65;
t58 = qJD(4) + qJD(5);
t83 = t34 * t58;
t152 = qJD(2) * t83;
t57 = qJDD(4) + qJDD(5);
t90 = t64 * t65 - t67 * t68;
t151 = t90 * t57;
t117 = qJDD(1) - g(3);
t85 = t98 * t69;
t150 = -t117 * t66 + t85;
t116 = qJD(2) * qJ(3);
t119 = t66 * qJD(1);
t41 = t116 + t119;
t120 = t41 * qJD(2);
t149 = -t120 - t148;
t106 = qJD(2) * qJD(4);
t110 = t68 * qJDD(2);
t147 = t65 * t106 - t110;
t125 = qJD(4) * t68;
t144 = qJD(5) * t68 + t125;
t70 = -pkin(2) - pkin(6);
t143 = (t116 + t41 - t119) * qJD(4) + qJDD(4) * t70;
t140 = g(3) * t66;
t139 = pkin(7) - t70;
t138 = -t58 * t83 - t151;
t128 = qJD(2) * t65;
t104 = t64 * t128;
t127 = qJD(2) * t68;
t27 = -t67 * t127 + t104;
t28 = t34 * qJD(2);
t137 = t27 * t28;
t136 = t34 * t57;
t135 = t62 * t66;
t134 = t63 * t66;
t118 = t69 * qJD(1);
t100 = qJD(3) - t118;
t33 = t70 * qJD(2) + t100;
t16 = -pkin(7) * t128 + t65 * t33;
t133 = t67 * t16;
t60 = t68 ^ 2;
t132 = t65 ^ 2 - t60;
t71 = qJD(4) ^ 2;
t72 = qJD(2) ^ 2;
t131 = t71 + t72;
t130 = qJD(2) * t27;
t129 = qJD(2) * t28;
t126 = qJD(4) * t65;
t124 = qJD(5) * t64;
t123 = qJD(5) * t67;
t121 = qJDD(2) * pkin(2);
t115 = qJDD(2) * t66;
t114 = qJDD(4) * t65;
t112 = t65 * qJDD(2);
t111 = t66 * qJDD(1);
t109 = t69 * qJDD(1);
t108 = qJD(1) * qJD(2);
t107 = qJD(2) * qJD(3);
t105 = qJDD(2) * qJ(3);
t37 = t139 * t68;
t102 = t68 * t106;
t50 = t65 * pkin(4) + qJ(3);
t101 = t58 * t68;
t99 = qJDD(3) - t109;
t17 = -pkin(7) * t127 + t68 * t33;
t97 = g(1) * t62 - g(2) * t63;
t45 = pkin(4) * t125 + qJD(3);
t78 = -t65 * t124 - t64 * t126;
t12 = t67 * t101 + t78;
t95 = -t12 * t58 - t136;
t14 = qJD(4) * pkin(4) + t17;
t94 = -t64 * t14 - t133;
t36 = t139 * t65;
t93 = -t67 * t36 - t64 * t37;
t92 = t64 * t36 - t67 * t37;
t91 = (-qJD(2) * pkin(2) + t100) * t66 + t41 * t69;
t51 = t66 * t108;
t88 = t51 + t99;
t87 = -qJD(5) * t104 - t147 * t64;
t82 = t90 * t58;
t80 = t131 * t69 + t115;
t79 = -qJDD(4) * t69 + 0.2e1 * t66 * t106;
t77 = -t85 - t140;
t76 = t99 - t148;
t30 = t50 * qJD(2) + t119;
t22 = t70 * qJDD(2) + t88;
t15 = t68 * t22;
t4 = qJDD(4) * pkin(4) + t147 * pkin(7) - t33 * t126 + t15;
t61 = qJ(4) + qJ(5);
t54 = sin(t61);
t55 = cos(t61);
t75 = t16 * t124 + (-t16 * t58 - t4) * t64 - g(1) * (-t54 * t134 - t62 * t55) - g(2) * (-t54 * t135 + t63 * t55) + t30 * t28 - t54 * t56;
t5 = t33 * t125 + t65 * t22 + (-t102 - t112) * pkin(7);
t74 = -g(1) * (t55 * t134 - t62 * t54) - g(2) * (t55 * t135 + t63 * t54) + t94 * qJD(5) + t30 * t27 + t67 * t4 - t64 * t5 + t55 * t56;
t23 = t105 + t111 + (qJD(3) + t118) * qJD(2);
t73 = -t140 + t105 + t107 - t70 * t71 + t23 + (-t98 - t108) * t69;
t6 = t67 * t110 - t64 * t112 - t152;
t53 = qJDD(4) * t68;
t39 = -t69 * qJDD(2) + t72 * t66;
t38 = t72 * t69 + t115;
t32 = qJD(4) * t37;
t31 = t139 * t126;
t25 = t88 - t121;
t10 = t111 + t50 * qJDD(2) + (t45 + t118) * qJD(2);
t8 = t27 ^ 2 - t28 ^ 2;
t7 = (qJD(2) * t101 + t112) * t67 + t87;
t2 = -t27 * t58 + (-t58 * t127 - t112) * t67 - t87;
t1 = t28 * t58 + t6;
t3 = [t117, 0, -t39, -t38, t39, t38, t91 * qJD(2) + t23 * t66 - t25 * t69 - g(3), 0, 0, 0, 0, 0, t80 * t65 + t79 * t68, -t79 * t65 + t80 * t68, 0, 0, 0, 0, 0, (-qJD(2) * t82 + t7) * t66 + ((t65 * t123 + t67 * t126 + t144 * t64) * t58 + t151 + t129) * t69, (t6 - t152) * t66 + (-(-t144 * t67 - t78) * t58 + t136 - t130) * t69; 0, qJDD(2), t109 + t148, t150, t76 - 0.2e1 * t121, 0.2e1 * t105 + 0.2e1 * t107 - t150, t23 * qJ(3) + t41 * qJD(3) - t25 * pkin(2) - g(3) * (t69 * pkin(2) + t66 * qJ(3)) - t91 * qJD(1) + t98 * (pkin(2) * t66 - qJ(3) * t69), t60 * qJDD(2) - 0.2e1 * t65 * t102, 0.2e1 * t132 * t106 - 0.2e1 * t65 * t110, -t71 * t65 + t53, -t71 * t68 - t114, 0, t143 * t68 + t73 * t65, -t143 * t65 + t73 * t68, t27 * t83 - t6 * t90, t27 * t12 + t28 * t83 - t6 * t34 + t7 * t90, t138, t95, 0, (-t93 * qJD(5) + t67 * t31 + t64 * t32) * t58 + t92 * t57 + t45 * t28 + t50 * t7 + t10 * t34 + t30 * t12 + t77 * t54 + (-t69 * t28 + t66 * t82) * qJD(1), -(t92 * qJD(5) + t64 * t31 - t67 * t32) * t58 - t93 * t57 - t45 * t27 + t50 * t6 - t10 * t90 - t30 * t83 + t77 * t55 + (t69 * t27 + t66 * t83) * qJD(1); 0, 0, 0, 0, qJDD(2), -t72, -t120 + t51 + t76 - t121, 0, 0, 0, 0, 0, -t131 * t65 + t53, -t131 * t68 - t114, 0, 0, 0, 0, 0, -t129 + t138, t95 + t130; 0, 0, 0, 0, 0, 0, 0, t68 * t72 * t65, -t132 * t72, t110, -t112, qJDD(4), t149 * t68 + t97 * t65 + t15, t97 * t68 + (-t22 - t149) * t65, -t137, t8, t1, t2, t57, -(-t64 * t17 - t133) * t58 + (-t58 * t124 - t28 * t127 + t67 * t57) * pkin(4) + t74, (-qJD(5) * t14 + t17 * t58 - t5) * t67 + (-t58 * t123 + t27 * t127 - t64 * t57) * pkin(4) + t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t137, t8, t1, t2, t57, -t94 * t58 + t74, (-t5 + (-qJD(5) + t58) * t14) * t67 + t75;];
tau_reg = t3;
