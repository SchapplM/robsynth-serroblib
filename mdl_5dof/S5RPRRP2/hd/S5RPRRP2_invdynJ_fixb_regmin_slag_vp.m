% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tau_reg [5x16]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:01:47
% EndTime: 2019-12-05 18:01:49
% DurationCPUTime: 0.64s
% Computational Cost: add. (996->154), mult. (1691->193), div. (0->0), fcn. (947->12), ass. (0->100)
t133 = qJ(5) + pkin(7);
t68 = sin(pkin(8));
t128 = pkin(1) * t68;
t105 = qJD(1) * t128;
t69 = cos(pkin(8));
t55 = pkin(1) * t69 + pkin(2);
t132 = -qJD(3) * t105 + t55 * qJDD(1);
t65 = qJ(1) + pkin(8);
t58 = qJ(3) + t65;
t53 = sin(t58);
t54 = cos(t58);
t130 = g(2) * t54 + g(3) * t53;
t39 = t55 * qJD(1);
t131 = qJD(3) * t39 + qJDD(1) * t128;
t116 = g(2) * t53 - g(3) * t54;
t72 = sin(qJ(3));
t75 = cos(qJ(3));
t115 = t75 * t128 + t72 * t55;
t99 = -t72 * t128 + t55 * t75;
t129 = -t131 * t75 - t132 * t72;
t63 = qJDD(1) + qJDD(3);
t127 = pkin(3) * t63;
t64 = qJD(1) + qJD(3);
t126 = pkin(3) * t64;
t74 = cos(qJ(4));
t125 = pkin(4) * t74;
t124 = g(1) * t74;
t21 = -t72 * t105 + t75 * t39;
t17 = -t21 - t126;
t71 = sin(qJ(4));
t96 = -t131 * t72 + t132 * t75;
t9 = -t96 - t127;
t123 = t17 * qJD(4) * t74 + t9 * t71;
t22 = t75 * t105 + t39 * t72;
t122 = t22 * t64;
t24 = t115 * qJD(3);
t121 = t24 * t64;
t119 = t71 * t63;
t118 = t74 * t63;
t98 = t133 * t64 + t22;
t11 = t74 * qJD(2) - t98 * t71;
t112 = qJD(4) * pkin(4);
t10 = t11 + t112;
t117 = t10 - t11;
t66 = t71 ^ 2;
t67 = t74 ^ 2;
t114 = -t66 - t67;
t113 = t66 - t67;
t27 = pkin(7) + t115;
t111 = -qJ(5) - t27;
t109 = qJD(4) * t71;
t108 = qJDD(2) - g(1);
t8 = t63 * pkin(7) - t129;
t80 = qJ(5) * t63 + qJD(2) * qJD(4) + qJD(5) * t64 + t8;
t88 = t98 * qJD(4);
t3 = (qJDD(2) - t88) * t71 + t80 * t74;
t107 = t3 * t74 + t116;
t106 = t17 * t109 + t130 * t74;
t104 = pkin(4) * t109;
t102 = t64 * t109;
t56 = pkin(3) + t125;
t100 = -t133 * t53 - t54 * t56;
t97 = qJD(4) * t133;
t95 = qJD(4) * t111;
t26 = -pkin(3) - t99;
t73 = sin(qJ(1));
t76 = cos(qJ(1));
t91 = g(2) * t76 + g(3) * t73;
t12 = t71 * qJD(2) + t98 * t74;
t90 = t10 * t71 - t12 * t74;
t89 = t133 * t54 - t53 * t56;
t87 = t96 + t130;
t77 = qJD(4) ^ 2;
t85 = -pkin(7) * t77 + t122 + t127;
t84 = -t26 * t63 - t27 * t77 - t121;
t83 = -t17 * t64 - t116 - t8;
t82 = -pkin(7) * qJDD(4) + (t21 - t126) * qJD(4);
t23 = t99 * qJD(3);
t81 = -qJDD(4) * t27 + (t26 * t64 - t23) * qJD(4);
t4 = pkin(4) * t102 - t56 * t63 + qJDD(5) - t96;
t79 = -t116 + t129;
t62 = t64 ^ 2;
t61 = t74 * qJ(5);
t59 = t74 * qJD(5);
t57 = t74 * qJDD(2);
t42 = pkin(7) * t74 + t61;
t41 = t133 * t71;
t33 = qJDD(4) * t74 - t71 * t77;
t32 = qJDD(4) * t71 + t74 * t77;
t29 = -t71 * qJD(5) - t74 * t97;
t28 = -t71 * t97 + t59;
t25 = 0.2e1 * t74 * t102 + t63 * t66;
t20 = t27 * t74 + t61;
t19 = t111 * t71;
t16 = -0.2e1 * t113 * t64 * qJD(4) + 0.2e1 * t71 * t118;
t13 = -t56 * t64 + qJD(5) - t21;
t6 = (-qJD(5) - t23) * t71 + t74 * t95;
t5 = t74 * t23 + t71 * t95 + t59;
t2 = qJDD(4) * pkin(4) - t80 * t71 - t74 * t88 + t57;
t1 = [qJDD(1), t91, -g(2) * t73 + g(3) * t76, (t91 + (t68 ^ 2 + t69 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t63, t99 * t63 - t121 + t87, -t115 * t63 - t23 * t64 + t79, t25, t16, t32, t33, 0, t81 * t71 + (t84 - t9) * t74 + t106, t81 * t74 + (-t84 - t130) * t71 + t123, (t20 * t63 + t5 * t64 + (-t19 * t64 - t10) * qJD(4)) * t74 + (-t19 * t63 - t6 * t64 - t2 + (-t20 * t64 - t12) * qJD(4)) * t71 + t107, t3 * t20 + t12 * t5 + t2 * t19 + t10 * t6 + t4 * (t26 - t125) + t13 * (t24 + t104) - g(2) * (-pkin(2) * cos(t65) - t76 * pkin(1) + t100) - g(3) * (-pkin(2) * sin(t65) - t73 * pkin(1) + t89); 0, 0, 0, t108, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t32, 0, -t90 * qJD(4) + t2 * t74 + t3 * t71 - g(1); 0, 0, 0, 0, t63, t87 + t122, t21 * t64 + t79, t25, t16, t32, t33, 0, t82 * t71 + (t85 - t9) * t74 + t106, t82 * t74 + (-t85 - t130) * t71 + t123, (-qJD(4) * t10 + t42 * t63) * t74 + (-qJD(4) * t12 + t41 * t63 - t2) * t71 + (t28 * t74 - t29 * t71 + t114 * t21 + (t41 * t74 - t42 * t71) * qJD(4)) * t64 + t107, t3 * t42 - t2 * t41 - t4 * t56 - g(2) * t100 - g(3) * t89 + (-t22 + t104) * t13 + (-t21 * t74 + t28) * t12 + (t21 * t71 + t29) * t10; 0, 0, 0, 0, 0, 0, 0, -t71 * t62 * t74, t113 * t62, t119, t118, qJDD(4), t83 * t71 - t124 + t57, -t108 * t71 + t83 * t74, -pkin(4) * t119 + (-t112 + t117) * t74 * t64, t117 * t12 + (-t124 + t2 + (-t13 * t64 - t116) * t71) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114 * t62, t90 * t64 - t130 + t4;];
tau_reg = t1;
