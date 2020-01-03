% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPPRP3
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRP3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:51:06
% EndTime: 2019-12-31 17:51:07
% DurationCPUTime: 0.66s
% Computational Cost: add. (813->177), mult. (1305->194), div. (0->0), fcn. (658->8), ass. (0->110)
t69 = sin(qJ(4));
t104 = t69 * qJDD(1);
t102 = qJD(1) * qJD(4);
t71 = cos(qJ(4));
t95 = t71 * t102;
t133 = t95 + t104;
t67 = cos(pkin(7));
t47 = -t67 * pkin(1) - pkin(2);
t39 = -pkin(6) + t47;
t116 = qJ(5) - t39;
t19 = t116 * t69;
t132 = qJD(4) * t19;
t60 = qJ(1) + pkin(7);
t53 = cos(t60);
t45 = g(2) * t53;
t52 = sin(t60);
t46 = g(1) * t52;
t121 = t45 - t46;
t125 = t69 * pkin(4);
t66 = sin(pkin(7));
t43 = t66 * pkin(1) + qJ(3);
t28 = t43 + t125;
t114 = qJD(1) * t28;
t18 = qJD(5) + t114;
t131 = -t18 * qJD(1) + t121;
t115 = pkin(1) * qJDD(1);
t48 = t66 * t115;
t120 = qJDD(1) * qJ(3) + t48;
t62 = qJD(3) * qJD(1);
t27 = t62 + t120;
t8 = t133 * pkin(4) + qJDD(5) + t27;
t90 = g(1) * t53 + g(2) * t52;
t130 = t8 - t90;
t129 = (qJD(5) + t18) * qJD(1);
t128 = t47 * qJDD(1);
t26 = t39 * qJD(1) + qJD(3);
t12 = -t69 * qJD(2) + t71 * t26;
t123 = t69 * t26;
t13 = t71 * qJD(2) + t123;
t100 = qJD(2) * qJD(4);
t112 = qJD(4) * t71;
t23 = t39 * qJDD(1) + qJDD(3);
t3 = t71 * qJDD(2) + t26 * t112 + (-t100 + t23) * t69;
t16 = t71 * t23;
t92 = -t69 * qJDD(2) + t16;
t4 = -t13 * qJD(4) + t92;
t77 = -(t12 * t69 - t13 * t71) * qJD(4) + t3 * t69 + t4 * t71;
t33 = qJD(1) * t43;
t127 = 0.2e1 * qJD(4) * t33 + qJDD(4) * t39;
t117 = qJD(4) * pkin(4);
t107 = qJ(5) * qJD(1);
t98 = t71 * t107;
t6 = t12 - t98;
t5 = t6 + t117;
t126 = t5 - t6;
t58 = g(3) * t69;
t124 = t52 * t69;
t74 = qJD(1) ^ 2;
t122 = t74 * t69;
t63 = t69 ^ 2;
t64 = t71 ^ 2;
t119 = -t63 - t64;
t118 = t63 - t64;
t20 = t116 * t71;
t113 = qJD(4) * t20;
t111 = qJDD(4) * pkin(4);
t109 = t33 * qJD(1);
t105 = qJDD(4) * t69;
t55 = t71 * qJDD(1);
t103 = qJ(5) * qJDD(1);
t101 = qJD(1) * qJD(5);
t72 = cos(qJ(1));
t99 = t72 * pkin(1) + t53 * pkin(2) + t52 * qJ(3);
t70 = sin(qJ(1));
t97 = -t70 * pkin(1) + t53 * qJ(3);
t96 = t69 * t102;
t94 = t18 + t114;
t29 = t119 * qJDD(1);
t73 = qJD(4) ^ 2;
t31 = qJDD(4) * t71 - t73 * t69;
t91 = t69 * t95;
t89 = g(1) * t70 - g(2) * t72;
t7 = -t69 * t107 + t13;
t87 = t5 * t71 + t69 * t7;
t86 = -t109 + t121;
t84 = t33 * qJD(3) + t27 * t43;
t83 = g(1) * t124 + g(3) * t71 - t3;
t82 = t71 * t45 + t58 + t92;
t81 = -t100 - t103;
t80 = qJDD(3) + t128;
t79 = t43 * qJDD(1) - t90;
t36 = pkin(4) * t112 + qJD(3);
t78 = qJD(1) * t36 + qJDD(1) * t28 + t130;
t76 = -t39 * t73 + t27 + t62 + t79;
t68 = -qJ(5) - pkin(6);
t65 = qJDD(2) - g(3);
t38 = t71 * t122;
t37 = qJ(5) * t96;
t32 = t118 * t74;
t30 = -t73 * t71 - t105;
t25 = t64 * qJDD(1) - 0.2e1 * t91;
t24 = t63 * qJDD(1) + 0.2e1 * t91;
t22 = t31 - t122;
t21 = -t105 + (-t73 - t74) * t71;
t11 = 0.2e1 * t118 * t102 - 0.2e1 * t69 * t55;
t10 = -t69 * qJD(5) - t113;
t9 = -t71 * qJD(5) + t132;
t2 = -t133 * qJ(5) - t69 * t101 + t3;
t1 = t111 + t16 + t37 + (-qJD(4) * t26 - qJDD(2)) * t69 + (t81 - t101) * t71;
t14 = [0, 0, 0, 0, 0, qJDD(1), t89, g(1) * t72 + g(2) * t70, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t67 * t115 - t121, -0.2e1 * t48 + t90, 0, (t89 + (t66 ^ 2 + t67 ^ 2) * t115) * pkin(1), qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(3) + t121 + 0.2e1 * t128, 0.2e1 * t62 + t79 + t120, t80 * t47 - g(1) * (-t52 * pkin(2) + t97) - g(2) * t99 + t84, t25, t11, t31, t24, t30, 0, t127 * t71 + t76 * t69, -t127 * t69 + t76 * t71, t39 * t29 - t121 - t77, -g(1) * ((-pkin(2) - pkin(6)) * t52 + t97) - g(2) * (t53 * pkin(6) + t99) + t77 * t39 + t84, t25, t11, t31, t24, t30, 0, -t20 * qJDD(4) + (t94 * t71 + t9) * qJD(4) + t78 * t69, t19 * qJDD(4) + (-t94 * t69 - t10) * qJD(4) + t78 * t71, (-t7 * qJD(4) + qJDD(1) * t20 - t1 + (-t9 + t132) * qJD(1)) * t71 + (qJD(4) * t5 + qJDD(1) * t19 - t2 + (-t10 - t113) * qJD(1)) * t69 - t121, -t2 * t19 + t7 * t10 - t1 * t20 + t5 * t9 + t8 * t28 + t18 * t36 - g(1) * (t53 * t125 + (-pkin(2) + t68) * t52 + t97) - g(2) * (pkin(4) * t124 - t53 * t68 + t99); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, 0, 0, 0, 0, 0, t30, -t31, 0, t3 * t71 - t4 * t69 - g(3) + (-t12 * t71 - t13 * t69) * qJD(4), 0, 0, 0, 0, 0, 0, t30, -t31, 0, -t87 * qJD(4) - t1 * t69 + t2 * t71 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t74, t80 + t86, 0, 0, 0, 0, 0, 0, t22, t21, t29, t77 + t86, 0, 0, 0, 0, 0, 0, t22, t21, t29, t1 * t71 + t2 * t69 + (-t5 * t69 + t7 * t71) * qJD(4) + t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t32, t55, -t38, -t104, qJDD(4), (-t109 - t46) * t71 + t82, t12 * qJD(4) + (t109 - t45) * t69 + t83, 0, 0, t38, -t32, t55, -t38, -t104, qJDD(4), 0.2e1 * t111 + t37 + (t7 - t123) * qJD(4) + (-pkin(4) * t122 - t129 - t46 + t81) * t71 + t82, -t64 * t74 * pkin(4) + (t6 + t98) * qJD(4) + (t103 - t45 + t129) * t69 + t83, -pkin(4) * t55 + (t117 - t126) * t69 * qJD(1), t126 * t7 + (t131 * t71 + t1 + t58) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t95 + t104, t55 - 0.2e1 * t96, t119 * t74, t87 * qJD(1) + t130;];
tau_reg = t14;
