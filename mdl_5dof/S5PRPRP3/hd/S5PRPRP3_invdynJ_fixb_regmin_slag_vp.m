% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRPRP3
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% tau_reg [5x16]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:14
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:14:01
% EndTime: 2021-01-15 15:14:04
% DurationCPUTime: 0.75s
% Computational Cost: add. (789->183), mult. (1606->234), div. (0->0), fcn. (1131->10), ass. (0->112)
t105 = qJD(1) * qJD(2);
t74 = sin(qJ(2));
t76 = cos(qJ(2));
t139 = t74 * qJDD(1) + t76 * t105;
t75 = cos(qJ(4));
t108 = qJ(5) * qJD(2);
t110 = t76 * qJD(1);
t49 = qJD(2) * pkin(2) + t110;
t116 = qJD(1) * t74;
t70 = cos(pkin(8));
t53 = t70 * t116;
t68 = sin(pkin(8));
t23 = t68 * t49 + t53;
t16 = qJD(2) * pkin(6) + t23;
t97 = t16 + t108;
t92 = t97 * t75;
t69 = sin(pkin(7));
t71 = cos(pkin(7));
t103 = -g(1) * t69 + g(2) * t71;
t94 = g(1) * t71 + g(2) * t69;
t65 = qJ(2) + pkin(8);
t58 = sin(t65);
t132 = g(3) * t58;
t59 = cos(t65);
t138 = -t94 * t59 - t132;
t118 = qJD(4) * pkin(4);
t63 = t75 * qJD(3);
t73 = sin(qJ(4));
t9 = -t97 * t73 + t63;
t6 = t9 + t118;
t137 = -t9 + t6;
t66 = t73 ^ 2;
t136 = pkin(4) * t66;
t131 = g(3) * t59;
t130 = g(3) * t73;
t129 = t70 * pkin(2);
t128 = t75 * pkin(4);
t127 = t69 * t73;
t126 = t69 * t75;
t125 = t71 * t73;
t124 = t71 * t75;
t52 = t68 * t116;
t32 = t70 * t110 - t52;
t123 = t75 * t32;
t78 = qJD(2) ^ 2;
t122 = t75 * t78;
t121 = qJD(4) * t123 + t59 * t130;
t67 = t75 ^ 2;
t120 = t66 - t67;
t119 = t66 + t67;
t55 = t68 * pkin(2) + pkin(6);
t117 = qJ(5) + t55;
t36 = t68 * t74 - t70 * t76;
t115 = qJD(2) * t36;
t114 = qJD(2) * t75;
t113 = qJDD(4) * pkin(4);
t22 = t70 * t49 - t52;
t15 = -qJD(2) * pkin(3) - t22;
t112 = t15 * qJD(2);
t111 = t73 * qJD(4);
t109 = qJDD(1) - g(3);
t60 = t73 * qJDD(2);
t106 = t75 * qJDD(2);
t104 = qJD(2) * qJD(4);
t62 = t76 * qJDD(1);
t33 = qJDD(2) * pkin(2) - t74 * t105 + t62;
t12 = t139 * t70 + t68 * t33;
t57 = pkin(3) + t128;
t101 = t73 * t104;
t30 = t68 * t110 + t53;
t100 = t32 * t111 + t30 * t114 + (g(1) * t124 + g(2) * t126) * t58;
t41 = -t57 - t129;
t11 = -t139 * t68 + t70 * t33;
t5 = pkin(4) * t101 - t57 * qJDD(2) + qJDD(5) - t11;
t99 = qJDD(2) * t41 + t5;
t98 = qJD(4) * t117;
t8 = qJDD(2) * pkin(6) + t12;
t96 = -qJD(4) * qJD(3) - t8;
t95 = 0.2e1 * t75 * t104;
t10 = t73 * qJD(3) + t92;
t93 = t10 * t75 - t6 * t73;
t37 = t68 * t76 + t70 * t74;
t91 = t94 * t58;
t56 = -pkin(3) - t129;
t77 = qJD(4) ^ 2;
t90 = t55 * t77 - t11 + (-pkin(3) + t56) * qJDD(2);
t61 = t75 * qJDD(3);
t89 = -g(1) * (-t59 * t125 + t126) - g(2) * (-t59 * t127 - t124) + t58 * t130 + t61;
t88 = -qJ(5) * qJDD(2) + t96;
t29 = t37 * qJD(2);
t86 = qJD(2) * t29 + qJDD(2) * t36 + t37 * t77;
t85 = 0.2e1 * t115 * qJD(4) - qJDD(4) * t37;
t84 = -qJDD(4) * t55 + (qJD(2) * t56 + t15) * qJD(4);
t83 = -g(3) * t76 + t74 * t94;
t14 = t16 * t111;
t82 = -g(1) * (-t59 * t124 - t127) - g(2) * (-t59 * t126 + t125) - t73 * qJDD(3) + t14 + t75 * t132;
t81 = -qJD(2) * t30 - t91;
t80 = qJD(2) * qJD(5) - t88;
t13 = -t57 * qJD(2) + qJD(5) - t22;
t79 = (-qJD(5) - t13) * qJD(2) + t88;
t72 = -qJ(5) - pkin(6);
t43 = qJDD(4) * t75 - t77 * t73;
t42 = qJDD(4) * t73 + t77 * t75;
t35 = t117 * t75;
t34 = t117 * t73;
t18 = -t73 * qJD(5) - t75 * t98;
t17 = t75 * qJD(5) - t73 * t98;
t4 = -t14 + (-qJ(5) * t104 + qJDD(3)) * t73 + t80 * t75;
t3 = -qJD(4) * t92 - t73 * t80 + t113 + t61;
t2 = t73 * t85 - t75 * t86;
t1 = t73 * t86 + t75 * t85;
t7 = [t109, 0, t76 * qJDD(2) - t78 * t74, -qJDD(2) * t74 - t78 * t76, -t11 * t36 - t115 * t23 + t12 * t37 - t22 * t29 - g(3), 0, 0, 0, 0, 0, t2, t1, t2, t1, (-qJD(2) * t115 + qJDD(2) * t37) * t119, t13 * t29 + t5 * t36 - g(3) - t93 * t115 + (-t3 * t73 + t4 * t75 + (-t10 * t73 - t6 * t75) * qJD(4)) * t37; 0, qJDD(2), t62 + t83, -t109 * t74 + t94 * t76, t22 * t30 - t23 * t32 + (t11 * t70 + t12 * t68 + t83) * pkin(2), t66 * qJDD(2) + t73 * t95, -0.2e1 * t120 * t104 + 0.2e1 * t73 * t106, t42, t43, 0, t84 * t73 + (-t90 - t131) * t75 + t100, t84 * t75 + (t81 + t90) * t73 + t121, -t34 * qJDD(4) + (-t99 - t131) * t75 + (t18 + (t13 + (t41 - t128) * qJD(2)) * t73) * qJD(4) + t100, -t35 * qJDD(4) + (t13 * t75 - t17 + (t41 * t75 + t136) * qJD(2)) * qJD(4) + (t81 + t99) * t73 + t121, (-qJD(4) * t6 + qJDD(2) * t35 + t4) * t75 + (-t10 * qJD(4) + qJDD(2) * t34 - t3) * t73 + (t17 * t75 - t18 * t73 - t119 * t32 + (t34 * t75 - t35 * t73) * qJD(4)) * qJD(2) + t138, t4 * t35 - t3 * t34 + t5 * t41 - g(3) * (t76 * pkin(2) + t59 * t57 - t58 * t72) + (t73 * t32 + t18) * t6 + (pkin(4) * t111 - t30) * t13 + (t17 - t123) * t10 + t94 * (pkin(2) * t74 + t57 * t58 + t59 * t72); 0, 0, 0, 0, qJDD(3) + t103, 0, 0, 0, 0, 0, t43, -t42, t43, -t42, 0, qJD(4) * t93 + t3 * t75 + t4 * t73 + t103; 0, 0, 0, 0, 0, -t73 * t122, t120 * t78, t60, t106, qJDD(4), (-t8 - t112) * t73 + t89, (-t73 * t16 + t63) * qJD(4) + (t96 - t112) * t75 + t82, 0.2e1 * t113 + (t10 - t92) * qJD(4) + (pkin(4) * t122 + t79) * t73 + t89, -t78 * t136 + (t73 * t108 + t9) * qJD(4) + t79 * t75 + t82, -pkin(4) * t60 + (-t118 + t137) * t114, t137 * t10 + (t3 + t103 * t75 + (-t13 * qJD(2) - t138) * t73) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t101 - t106, t60 + t95, -t119 * t78, -qJD(2) * t93 + t131 + t5 - t91;];
tau_reg = t7;
