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
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:45:32
% EndTime: 2020-01-03 11:45:35
% DurationCPUTime: 0.82s
% Computational Cost: add. (996->151), mult. (1691->192), div. (0->0), fcn. (947->12), ass. (0->100)
t65 = qJ(1) + pkin(8);
t58 = qJ(3) + t65;
t53 = sin(t58);
t54 = cos(t58);
t93 = g(2) * t54 + g(3) * t53;
t68 = sin(pkin(8));
t130 = pkin(1) * t68;
t69 = cos(pkin(8));
t55 = pkin(1) * t69 + pkin(2);
t41 = t55 * qJD(1);
t132 = qJD(3) * t41 + qJDD(1) * t130;
t106 = qJD(1) * t130;
t133 = -qJD(3) * t106 + t55 * qJDD(1);
t72 = sin(qJ(3));
t75 = cos(qJ(3));
t96 = -t132 * t72 + t133 * t75;
t83 = t93 - t96;
t63 = qJDD(1) + qJDD(3);
t129 = pkin(3) * t63;
t88 = t129 - t83;
t70 = -qJ(5) - pkin(7);
t116 = -g(2) * t53 + g(3) * t54;
t115 = t75 * t130 + t72 * t55;
t100 = -t72 * t130 + t55 * t75;
t131 = -t132 * t75 - t133 * t72;
t64 = qJD(1) + qJD(3);
t128 = pkin(3) * t64;
t74 = cos(qJ(4));
t127 = pkin(4) * t74;
t126 = g(1) * t74;
t22 = t75 * t106 + t41 * t72;
t123 = t22 * t64;
t24 = t115 * qJD(3);
t122 = t24 * t64;
t71 = sin(qJ(4));
t120 = t71 * t63;
t119 = t74 * t63;
t98 = -t70 * t64 + t22;
t11 = t74 * qJD(2) - t98 * t71;
t112 = qJD(4) * pkin(4);
t10 = t11 + t112;
t118 = t10 - t11;
t56 = pkin(3) + t127;
t117 = t53 * t56 + t54 * t70;
t66 = t71 ^ 2;
t67 = t74 ^ 2;
t114 = -t66 - t67;
t113 = t66 - t67;
t27 = pkin(7) + t115;
t111 = -qJ(5) - t27;
t109 = qJD(4) * t71;
t108 = qJDD(2) - g(1);
t8 = t63 * pkin(7) - t131;
t80 = qJ(5) * t63 + qJD(2) * qJD(4) + qJD(5) * t64 + t8;
t89 = t98 * qJD(4);
t3 = (qJDD(2) - t89) * t71 + t80 * t74;
t107 = t3 * t74 + t116;
t105 = pkin(4) * t109;
t103 = t64 * t109;
t102 = -t53 * t70 + t54 * t56;
t21 = -t72 * t106 + t75 * t41;
t17 = -t21 - t128;
t99 = t17 * qJD(4) * t74 - t88 * t71;
t97 = qJD(4) * t70;
t95 = qJD(4) * t111;
t26 = -pkin(3) - t100;
t73 = sin(qJ(1));
t76 = cos(qJ(1));
t91 = -g(2) * t76 - g(3) * t73;
t12 = t71 * qJD(2) + t98 * t74;
t90 = t10 * t71 - t12 * t74;
t77 = qJD(4) ^ 2;
t86 = pkin(7) * t77 - t123 - t129;
t85 = t26 * t63 + t27 * t77 + t122;
t84 = -t17 * t64 - t116 - t8;
t82 = -pkin(7) * qJDD(4) + (t21 - t128) * qJD(4);
t23 = t100 * qJD(3);
t81 = -qJDD(4) * t27 + (t26 * t64 - t23) * qJD(4);
t4 = pkin(4) * t103 - t56 * t63 + qJDD(5) - t96;
t79 = -t116 + t131;
t62 = t64 ^ 2;
t61 = t74 * qJ(5);
t59 = t74 * qJD(5);
t57 = t74 * qJDD(2);
t44 = pkin(7) * t74 + t61;
t43 = t70 * t71;
t35 = qJDD(4) * t74 - t71 * t77;
t34 = qJDD(4) * t71 + t74 * t77;
t29 = -t71 * qJD(5) + t74 * t97;
t28 = t71 * t97 + t59;
t25 = 0.2e1 * t74 * t103 + t63 * t66;
t20 = t27 * t74 + t61;
t19 = t111 * t71;
t16 = -0.2e1 * t113 * t64 * qJD(4) + 0.2e1 * t71 * t119;
t14 = t17 * t109;
t13 = -t56 * t64 + qJD(5) - t21;
t6 = (-qJD(5) - t23) * t71 + t74 * t95;
t5 = t74 * t23 + t71 * t95 + t59;
t2 = qJDD(4) * pkin(4) - t80 * t71 - t74 * t89 + t57;
t1 = [qJDD(1), t91, g(2) * t73 - g(3) * t76, (t91 + (t68 ^ 2 + t69 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t63, t100 * t63 - t122 - t83, -t115 * t63 - t23 * t64 + t79, t25, t16, t34, t35, 0, t14 + t81 * t71 + (-t85 + t88) * t74, t85 * t71 + t81 * t74 + t99, (t20 * t63 + t5 * t64 + (-t19 * t64 - t10) * qJD(4)) * t74 + (-t19 * t63 - t6 * t64 - t2 + (-t20 * t64 - t12) * qJD(4)) * t71 + t107, t3 * t20 + t12 * t5 + t2 * t19 + t10 * t6 + t4 * (t26 - t127) + t13 * (t24 + t105) - g(2) * (pkin(2) * cos(t65) + t76 * pkin(1) + t102) - g(3) * (pkin(2) * sin(t65) + t73 * pkin(1) + t117); 0, 0, 0, t108, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t34, 0, -t90 * qJD(4) + t2 * t74 + t3 * t71 - g(1); 0, 0, 0, 0, t63, -t83 + t123, t21 * t64 + t79, t25, t16, t34, t35, 0, t14 + t82 * t71 + (-t86 + t88) * t74, t86 * t71 + t82 * t74 + t99, (-qJD(4) * t10 + t44 * t63) * t74 + (-qJD(4) * t12 - t43 * t63 - t2) * t71 + (t28 * t74 - t29 * t71 + t114 * t21 + (-t43 * t74 - t44 * t71) * qJD(4)) * t64 + t107, t3 * t44 + t2 * t43 - t4 * t56 - g(2) * t102 - g(3) * t117 + (-t22 + t105) * t13 + (-t21 * t74 + t28) * t12 + (t21 * t71 + t29) * t10; 0, 0, 0, 0, 0, 0, 0, -t71 * t62 * t74, t113 * t62, t120, t119, qJDD(4), t84 * t71 - t126 + t57, -t108 * t71 + t84 * t74, -pkin(4) * t120 + (-t112 + t118) * t74 * t64, t118 * t12 + (-t126 + t2 + (-t13 * t64 - t116) * t71) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114 * t62, t90 * t64 + t4 + t93;];
tau_reg = t1;
