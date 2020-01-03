% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% tau_reg [5x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:25:51
% EndTime: 2020-01-03 11:25:57
% DurationCPUTime: 0.87s
% Computational Cost: add. (857->189), mult. (1699->263), div. (0->0), fcn. (1071->10), ass. (0->110)
t64 = cos(pkin(7));
t48 = -pkin(1) * t64 - pkin(2);
t61 = sin(pkin(8));
t63 = cos(pkin(8));
t30 = -pkin(3) * t63 - pkin(6) * t61 + t48;
t20 = t30 * qJD(1) + qJD(3);
t62 = sin(pkin(7));
t43 = pkin(1) * t62 + qJ(3);
t35 = t43 * qJD(1);
t27 = qJD(2) * t61 + t35 * t63;
t66 = sin(qJ(4));
t68 = cos(qJ(4));
t111 = qJD(1) * t61;
t89 = qJ(5) * t111;
t7 = (t20 - t89) * t66 + t27 * t68;
t137 = qJD(4) * t7;
t58 = qJ(1) + pkin(7);
t52 = sin(t58);
t53 = cos(t58);
t134 = g(2) * t53 + g(3) * t52;
t101 = qJDD(1) * t48;
t34 = qJDD(3) + t101;
t136 = t134 + t34;
t96 = qJD(1) * qJD(3);
t31 = qJDD(1) * t43 + t96;
t120 = t63 * t66;
t21 = -t52 * t120 - t53 * t68;
t23 = t53 * t120 - t52 * t68;
t133 = -g(2) * t21 - g(3) * t23;
t110 = qJD(1) * t63;
t41 = -qJD(4) + t110;
t87 = t68 * t20 - t27 * t66;
t6 = -t68 * t89 + t87;
t3 = -pkin(4) * t41 + t6;
t132 = -t6 + t3;
t131 = pkin(4) * t66;
t126 = t31 * t61;
t50 = t63 * qJDD(2);
t18 = -t50 + t126;
t128 = t18 * t61;
t97 = t63 * qJDD(1);
t40 = -qJDD(4) + t97;
t125 = t40 * t63;
t124 = (pkin(4) * t68 + pkin(3)) * t63;
t123 = t52 * t66;
t122 = t53 * t66;
t56 = t61 ^ 2;
t70 = qJD(1) ^ 2;
t121 = t56 * t70;
t119 = t63 * t68;
t105 = qJD(4) * t68;
t107 = qJD(3) * t68;
t118 = t30 * t105 + t63 * t107;
t32 = t43 * t119;
t117 = t66 * t30 + t32;
t115 = t63 ^ 2 + t56;
t59 = t66 ^ 2;
t60 = t68 ^ 2;
t114 = -t59 - t60;
t113 = t59 - t60;
t112 = qJ(5) * t61;
t109 = qJD(1) * t66;
t108 = qJD(3) * t63;
t106 = qJD(4) * t66;
t104 = qJD(5) * t61;
t103 = -qJD(4) - t41;
t100 = qJDD(1) * t66;
t99 = qJDD(1) * t68;
t98 = t61 * qJDD(1);
t95 = qJD(1) * qJD(4);
t94 = qJD(1) * qJD(5);
t69 = cos(qJ(1));
t93 = t69 * pkin(1) + t53 * pkin(2) + t52 * qJ(3);
t92 = t68 * t112;
t91 = t41 * t106;
t90 = t63 * t106;
t88 = t68 * t95;
t86 = t40 - t97;
t85 = t40 + t97;
t84 = qJD(1) * t103;
t83 = t61 * pkin(4) * t88 + t98 * t131 + qJDD(5) - t50;
t67 = sin(qJ(1));
t82 = t67 * pkin(1) + t52 * pkin(2) - qJ(3) * t53;
t80 = -g(2) * t52 + g(3) * t53;
t79 = -g(2) * t69 - g(3) * t67;
t78 = t3 * t68 + t66 * t7;
t77 = t3 * t66 - t68 * t7;
t19 = qJDD(2) * t61 + t31 * t63;
t76 = t19 * t63 + t128;
t51 = t63 * qJD(2);
t26 = t35 * t61 - t51;
t75 = t26 * t61 + t27 * t63;
t17 = t30 * qJDD(1) + qJDD(3);
t74 = -t20 * t105 + t27 * t106 - t66 * t17 - t68 * t19;
t73 = t101 + t136;
t72 = -t41 ^ 2 - t121;
t33 = t90 * t111;
t29 = t68 * t30;
t24 = t53 * t119 + t123;
t22 = t52 * t119 - t122;
t14 = t68 * t17;
t11 = qJD(5) - t51 + (pkin(4) * t109 + t35) * t61;
t10 = -t66 * t112 + t117;
t9 = t83 + t126;
t8 = -t92 + t29 + (-t43 * t66 - pkin(4)) * t63;
t5 = -t66 * t108 - t68 * t104 + (-t32 + (-t30 + t112) * t66) * qJD(4);
t4 = -t66 * t104 + (-t43 * t120 - t92) * qJD(4) + t118;
t2 = (-t66 * t94 + (-t88 - t100) * qJ(5)) * t61 - t74;
t1 = -pkin(4) * t40 - t66 * t19 + t14 + (-qJ(5) * qJDD(1) - t94) * t68 * t61 - t137;
t12 = [qJDD(1), t79, g(2) * t67 - g(3) * t69, (t79 + (t62 ^ 2 + t64 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), -t73 * t63, t73 * t61, t31 * t115 + t76 + t80, -g(2) * t93 - g(3) * t82 + t75 * qJD(3) + t34 * t48 + t76 * t43, (qJDD(1) * t60 - 0.2e1 * t66 * t88) * t56, 0.2e1 * (t113 * t95 - t66 * t99) * t56, t33 + (-t85 * t68 + t91) * t61, (t85 * t66 + (t41 + t110) * t105) * t61, t125, -g(2) * t24 - g(3) * t22 - t14 * t63 - t29 * t40 + ((qJD(1) * t56 + t41 * t63) * t43 + t75) * t105 + (-(-qJD(4) * t30 - t108) * t41 - (-qJD(4) * t20 - t19) * t63 + t56 * t96 + t128 + (qJDD(1) * t56 + t125) * t43) * t66, (-t43 * t90 + t118) * t41 + t117 * t40 - t74 * t63 + g(2) * t23 - g(3) * t21 + (-t26 * t106 + t18 * t68) * t61 + (t43 * t99 + (-t43 * t106 + t107) * qJD(1)) * t56, ((-t137 - qJDD(1) * t8 - t1 + (-qJD(4) * t10 - t5) * qJD(1)) * t68 + (qJD(4) * t3 - qJDD(1) * t10 - t2 + (qJD(4) * t8 - t4) * qJD(1)) * t66 - t134) * t61, t2 * t10 + t7 * t4 + t1 * t8 + t3 * t5 - g(2) * (pkin(4) * t123 + t53 * t124 + t93) - g(3) * (-pkin(4) * t122 + t52 * t124 + t82) + (t9 * (t43 + t131) + t11 * (pkin(4) * t105 + qJD(3)) + t134 * (-qJ(5) - pkin(6))) * t61; 0, 0, 0, qJDD(2) - g(1), 0, 0, 0, -t18 * t63 + t19 * t61 - g(1), 0, 0, 0, 0, 0, (t86 * t66 + (t41 - t110) * t105) * t61, t33 + (t86 * t68 - t91) * t61, 0, -t63 * t9 - g(1) + (-t78 * qJD(4) - t1 * t66 + t2 * t68) * t61; 0, 0, 0, 0, -t97, t98, -t115 * t70, -t75 * qJD(1) + t136, 0, 0, 0, 0, 0, -t40 * t68 + t72 * t66, t40 * t66 + t72 * t68, t114 * t98, t1 * t68 + t2 * t66 - t77 * qJD(4) + (-t11 * t61 + t77 * t63) * qJD(1) + t134; 0, 0, 0, 0, 0, 0, 0, 0, t68 * t66 * t121, -t113 * t121, (t66 * t84 + t99) * t61, (t68 * t84 - t100) * t61, -t40, t14 + (t103 * t27 - t26 * t111) * t68 + (g(1) * t61 + t103 * t20 - t19) * t66 + t133, -t87 * t41 + g(2) * t22 - g(3) * t24 + (g(1) * t68 + t26 * t109) * t61 + t74, (-pkin(4) * t99 + (pkin(4) * qJD(4) - t132) * t109) * t61, t132 * t7 + (t1 + (-t11 * t68 * qJD(1) + g(1) * t66) * t61 + t133) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114 * t121, g(1) * t63 + (t78 * qJD(1) + t31 + t80) * t61 + t83;];
tau_reg = t12;
