% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRPRP4
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
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:36:08
% EndTime: 2019-12-05 15:36:11
% DurationCPUTime: 0.66s
% Computational Cost: add. (802->165), mult. (1591->215), div. (0->0), fcn. (1135->10), ass. (0->109)
t107 = qJD(1) * qJD(2);
t73 = sin(qJ(2));
t75 = cos(qJ(2));
t138 = qJDD(1) * t73 + t75 * t107;
t63 = t75 * qJDD(1);
t35 = qJDD(2) * pkin(2) - t73 * t107 + t63;
t68 = sin(pkin(8));
t70 = cos(pkin(8));
t12 = t138 * t70 + t68 * t35;
t140 = qJDD(2) * pkin(6) + qJD(3) * qJD(4) + t12;
t118 = qJD(1) * t75;
t47 = qJD(2) * pkin(2) + t118;
t119 = qJD(1) * t73;
t51 = t70 * t119;
t23 = t68 * t47 + t51;
t19 = qJD(2) * pkin(6) + t23;
t72 = sin(qJ(4));
t129 = t19 * t72;
t74 = cos(qJ(4));
t15 = qJD(3) * t74 - t129;
t139 = qJD(5) - t15;
t65 = qJ(2) + pkin(8);
t58 = sin(t65);
t69 = sin(pkin(7));
t71 = cos(pkin(7));
t92 = g(1) * t71 + g(2) * t69;
t137 = t92 * t58;
t10 = -qJD(4) * pkin(4) + t139;
t128 = t19 * t74;
t16 = qJD(3) * t72 + t128;
t13 = qJD(4) * qJ(5) + t16;
t66 = t72 ^ 2;
t67 = t74 ^ 2;
t120 = t66 + t67;
t136 = qJD(2) * t120;
t113 = qJDD(4) * pkin(4);
t135 = qJDD(5) - t113;
t134 = pkin(2) * t70;
t131 = g(3) * t58;
t59 = cos(t65);
t130 = g(3) * t59;
t127 = t69 * t72;
t126 = t69 * t74;
t125 = t71 * t72;
t124 = t71 * t74;
t123 = t72 * t74;
t90 = pkin(4) * t72 - qJ(5) * t74;
t29 = t90 * qJD(4) - qJD(5) * t72;
t31 = t68 * t118 + t51;
t122 = t29 - t31;
t121 = t66 - t67;
t117 = qJD(2) * t29;
t116 = qJD(2) * t31;
t36 = t68 * t73 - t70 * t75;
t115 = qJD(2) * t36;
t114 = qJD(2) * t72;
t112 = qJDD(1) - g(3);
t54 = pkin(2) * t68 + pkin(6);
t109 = qJDD(4) * t54;
t108 = t74 * qJDD(2);
t106 = qJD(2) * qJD(4);
t104 = qJDD(4) * qJ(5);
t103 = t72 * qJDD(3) + t140 * t74;
t77 = qJD(2) ^ 2;
t102 = t77 * t123;
t101 = -g(1) * t69 + g(2) * t71;
t99 = t15 + t129;
t50 = t68 * t119;
t22 = t47 * t70 - t50;
t98 = qJD(4) * t128 - t74 * qJDD(3) + t140 * t72;
t33 = t70 * t118 - t50;
t97 = t72 * t33 * qJD(4) + t74 * t116 + (g(1) * t124 + g(2) * t126) * t58;
t88 = pkin(4) * t74 + qJ(5) * t72 + pkin(3);
t14 = -t88 * qJD(2) - t22;
t34 = -t88 - t134;
t96 = qJD(2) * t34 + t14;
t18 = -qJD(2) * pkin(3) - t22;
t55 = -pkin(3) - t134;
t95 = qJD(2) * t55 + t18;
t94 = qJDD(2) * t120;
t11 = -t138 * t68 + t35 * t70;
t76 = qJD(4) ^ 2;
t91 = t54 * t76 + t130;
t89 = t10 * t72 + t13 * t74;
t37 = t68 * t75 + t70 * t73;
t26 = t59 * t126 - t125;
t28 = t59 * t124 + t127;
t87 = g(1) * t28 + g(2) * t26 - t103;
t30 = t37 * qJD(2);
t86 = qJD(2) * t30 + qJDD(2) * t36 + t37 * t76;
t5 = -t88 * qJDD(2) - t11 + t117;
t85 = -qJDD(2) * t34 - t5 - t91;
t84 = -t11 + t91 + (-pkin(3) + t55) * qJDD(2);
t83 = 0.2e1 * t115 * qJD(4) - qJDD(4) * t37;
t82 = -g(3) * t75 + t92 * t73;
t25 = t59 * t127 + t124;
t27 = t59 * t125 - t126;
t81 = g(1) * t27 + g(2) * t25 + t72 * t131 - t98;
t80 = qJD(4) * t16 + t81;
t3 = t104 + (qJD(5) - t129) * qJD(4) + t103;
t4 = t98 + t135;
t79 = t3 * t74 + t4 * t72 + (t10 * t74 - t13 * t72) * qJD(4);
t61 = t72 * qJDD(2);
t43 = qJDD(4) * t74 - t72 * t76;
t42 = qJDD(4) * t72 + t74 * t76;
t38 = t90 * qJD(2);
t2 = t83 * t72 - t86 * t74;
t1 = t86 * t72 + t83 * t74;
t6 = [t112, 0, qJDD(2) * t75 - t73 * t77, -qJDD(2) * t73 - t75 * t77, -t11 * t36 - t115 * t23 + t12 * t37 - t22 * t30 - g(3), 0, 0, 0, 0, 0, t2, t1, t2, -t115 * t136 + t37 * t94, -t1, -t115 * t89 + t14 * t30 + t36 * t5 + t79 * t37 - g(3); 0, qJDD(2), t63 + t82, -t112 * t73 + t92 * t75, t22 * t31 - t23 * t33 + (t11 * t70 + t12 * t68 + t82) * pkin(2), qJDD(2) * t66 + 0.2e1 * t106 * t123, -0.2e1 * t121 * t106 + 0.2e1 * t72 * t108, t42, t43, 0, (t95 * qJD(4) - t109) * t72 - t84 * t74 + t97, (-t109 + (t33 + t95) * qJD(4)) * t74 + (-t116 + t84 - t137) * t72, (t96 * qJD(4) - t109) * t72 + (t85 - t117) * t74 + t97, -t33 * t136 + t54 * t94 - t92 * t59 - t131 + t79, (t109 + (-t33 - t96) * qJD(4)) * t74 + (-t122 * qJD(2) + t137 + t85) * t72, t5 * t34 - g(3) * (pkin(2) * t75 + pkin(6) * t58) - t88 * t130 - t89 * t33 + t122 * t14 + t79 * t54 + t92 * (pkin(2) * t73 - pkin(6) * t59 + t88 * t58); 0, 0, 0, 0, qJDD(3) + t101, 0, 0, 0, 0, 0, t43, -t42, t43, 0, t42, t89 * qJD(4) + t3 * t72 - t4 * t74 + t101; 0, 0, 0, 0, 0, -t102, t121 * t77, t61, t108, qJDD(4), -t18 * t114 + t80, (-qJD(2) * t18 + t131) * t74 + t99 * qJD(4) + t87, 0.2e1 * t113 - qJDD(5) + (-t14 * t72 + t38 * t74) * qJD(2) + t80, -t90 * qJDD(2), -t74 * t131 + 0.2e1 * t104 + (t14 * t74 + t38 * t72) * qJD(2) + (0.2e1 * qJD(5) - t99) * qJD(4) - t87, t3 * qJ(5) - t4 * pkin(4) - t14 * t38 - t10 * t16 - g(1) * (-pkin(4) * t27 + qJ(5) * t28) - g(2) * (-pkin(4) * t25 + qJ(5) * t26) + t90 * t131 + t139 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) - t102, t61, -t66 * t77 - t76, -qJD(4) * t13 + t14 * t114 + t135 - t81;];
tau_reg = t6;
