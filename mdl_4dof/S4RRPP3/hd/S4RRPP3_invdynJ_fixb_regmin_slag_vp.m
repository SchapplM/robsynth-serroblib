% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% tau_reg [4x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:36
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:36:05
% EndTime: 2021-01-15 10:36:08
% DurationCPUTime: 0.68s
% Computational Cost: add. (924->194), mult. (2150->246), div. (0->0), fcn. (1373->8), ass. (0->98)
t81 = sin(qJ(1));
t83 = cos(qJ(1));
t128 = g(1) * t81 - g(2) * t83;
t98 = g(1) * t83 + g(2) * t81;
t77 = sin(pkin(6));
t78 = cos(pkin(6));
t80 = sin(qJ(2));
t82 = cos(qJ(2));
t49 = t77 * t82 + t78 * t80;
t42 = t49 * qJD(1);
t38 = t42 ^ 2;
t116 = t78 * t82;
t105 = qJD(1) * t116;
t113 = qJD(1) * t80;
t39 = t113 * t77 - t105;
t129 = -t39 ^ 2 - t38;
t79 = -qJ(3) - pkin(5);
t103 = t79 * t80;
t100 = qJD(2) * t79;
t93 = -t80 * qJD(3) + t100 * t82;
t20 = qJDD(2) * pkin(2) + qJD(1) * t93 + qJDD(1) * t103;
t36 = t82 * qJD(3) + t100 * t80;
t56 = t79 * t82;
t26 = qJD(1) * t36 - qJDD(1) * t56;
t5 = t77 * t20 + t78 * t26;
t73 = qJ(2) + pkin(6);
t70 = sin(t73);
t71 = cos(t73);
t127 = g(3) * t70 + t98 * t71 - t5;
t32 = t103 * t77 - t78 * t56;
t126 = -t32 * qJDD(2) - t128 * t70;
t125 = pkin(2) * t80;
t121 = g(3) * t71;
t120 = g(3) * t82;
t119 = t82 * pkin(2);
t69 = pkin(1) + t119;
t54 = -qJD(1) * t69 + qJD(3);
t10 = t39 * pkin(3) - t42 * qJ(4) + t54;
t118 = t10 * t42;
t52 = qJD(1) * t56;
t117 = t77 * t52;
t45 = t78 * t52;
t115 = t81 * t79;
t4 = t78 * t20 - t77 * t26;
t51 = qJD(1) * t103;
t47 = qJD(2) * pkin(2) + t51;
t23 = t77 * t47 - t45;
t75 = t80 ^ 2;
t114 = -t82 ^ 2 + t75;
t112 = qJD(2) * t80;
t28 = t78 * t51 + t117;
t111 = qJD(4) - t28;
t109 = t80 * qJDD(1);
t108 = t82 * qJDD(1);
t107 = qJD(1) * qJD(2);
t106 = pkin(2) * t112;
t102 = t80 * t107;
t101 = t82 * t107;
t96 = -t78 * t108 + t109 * t77;
t22 = t78 * t47 + t117;
t3 = -qJDD(2) * pkin(3) + qJDD(4) - t4;
t95 = -0.2e1 * pkin(1) * t107 - pkin(5) * qJDD(2);
t41 = t49 * qJD(2);
t31 = -t103 * t78 - t77 * t56;
t94 = -t31 * qJDD(2) + t128 * t71;
t35 = pkin(2) * t102 - qJDD(1) * t69 + qJDD(3);
t92 = qJDD(1) * t49 - t77 * t102;
t27 = t77 * t51 - t45;
t91 = t27 * qJD(2) + t98 * t70 - t121 + t4;
t84 = qJD(2) ^ 2;
t90 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t84 + t128;
t85 = qJD(1) ^ 2;
t89 = pkin(1) * t85 - pkin(5) * qJDD(1) + t98;
t12 = t77 * t36 - t78 * t93;
t13 = t78 * t36 + t77 * t93;
t24 = qJD(1) * t41 + t96;
t25 = t101 * t78 + t92;
t88 = t12 * t42 - t13 * t39 - t32 * t24 + t31 * t25 - t98;
t87 = t24 * pkin(3) - t25 * qJ(4) + t35;
t86 = 0.2e1 * t42 * qJD(2) + t96;
t74 = qJDD(2) * qJ(4);
t68 = t79 * t83;
t67 = -t78 * pkin(2) - pkin(3);
t65 = t77 * pkin(2) + qJ(4);
t55 = -t77 * pkin(3) + qJ(4) * t78;
t53 = pkin(3) * t78 + qJ(4) * t77 + pkin(2);
t48 = t77 * t80 - t116;
t44 = qJD(2) * t116 - t112 * t77;
t29 = t53 * t82 + t55 * t80 + pkin(1);
t21 = t48 * pkin(3) - t49 * qJ(4) - t69;
t17 = qJD(2) * qJ(4) + t23;
t14 = -qJD(2) * pkin(3) + qJD(4) - t22;
t11 = pkin(2) * t113 + t42 * pkin(3) + t39 * qJ(4);
t9 = (-t39 + t105) * qJD(2) + t92;
t8 = t41 * pkin(3) - t44 * qJ(4) - t49 * qJD(4) + t106;
t2 = qJD(2) * qJD(4) + t5 + t74;
t1 = -t42 * qJD(4) + t87;
t6 = [qJDD(1), t128, t98, t75 * qJDD(1) + 0.2e1 * t101 * t80, -0.2e1 * t107 * t114 + 0.2e1 * t108 * t80, qJDD(2) * t80 + t84 * t82, qJDD(2) * t82 - t84 * t80, 0, t80 * t95 + t82 * t90, -t80 * t90 + t82 * t95, -t69 * t24 + t35 * t48 + t54 * t41 + (t39 * t125 - t12) * qJD(2) + t94, -t69 * t25 + t35 * t49 + t54 * t44 + (t42 * t125 - t13) * qJD(2) + t126, -t22 * t44 - t23 * t41 - t4 * t49 - t5 * t48 + t88, t5 * t32 + t23 * t13 - t4 * t31 - t22 * t12 - t35 * t69 + t54 * t106 - g(1) * (-t81 * t69 - t68) - g(2) * (t83 * t69 - t115), -t12 * qJD(2) + t1 * t48 + t10 * t41 + t21 * t24 + t8 * t39 + t94, t14 * t44 - t17 * t41 - t2 * t48 + t3 * t49 + t88, t13 * qJD(2) - t1 * t49 - t10 * t44 - t21 * t25 - t8 * t42 - t126, t2 * t32 + t17 * t13 + t1 * t21 + t10 * t8 + t3 * t31 + t14 * t12 - g(1) * (-t29 * t81 - t68) - g(2) * (t29 * t83 - t115); 0, 0, 0, -t80 * t85 * t82, t114 * t85, t109, t108, qJDD(2), t80 * t89 - t120, g(3) * t80 + t82 * t89, -t54 * t42 + (qJDD(2) * t78 - t39 * t113) * pkin(2) + t91, t28 * qJD(2) + t54 * t39 + (-qJDD(2) * t77 - t113 * t42) * pkin(2) + t127, (t23 - t27) * t42 + (-t22 + t28) * t39 + (-t24 * t77 - t25 * t78) * pkin(2), t22 * t27 - t23 * t28 + (-t120 + t4 * t78 + t5 * t77 + (-qJD(1) * t54 + t98) * t80) * pkin(2), -t118 - t11 * t39 - qJDD(4) + (pkin(3) - t67) * qJDD(2) + t91, -t65 * t24 + t67 * t25 + (t17 - t27) * t42 + (t14 - t111) * t39, t65 * qJDD(2) - t10 * t39 + t11 * t42 + t74 + (0.2e1 * qJD(4) - t28) * qJD(2) - t127, t2 * t65 + t3 * t67 - t10 * t11 - t14 * t27 - g(3) * (t71 * pkin(3) + t70 * qJ(4) + t119) - t98 * (-t53 * t80 + t55 * t82) + t111 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t9, t129, t22 * t42 + t23 * t39 - t128 + t35, t86, t129, -t9, t17 * t39 + (-qJD(4) - t14) * t42 + t87 - t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42 * t39 - qJDD(2), (t39 + t105) * qJD(2) + t92, -t38 - t84, -t17 * qJD(2) - t98 * t49 + t118 + t121 + t3;];
tau_reg = t6;
