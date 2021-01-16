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
% tau_reg [5x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:36
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
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
% StartTime: 2021-01-15 12:36:27
% EndTime: 2021-01-15 12:36:31
% DurationCPUTime: 0.83s
% Computational Cost: add. (1303->203), mult. (2181->238), div. (0->0), fcn. (1229->12), ass. (0->132)
t81 = sin(pkin(8));
t159 = pkin(1) * t81;
t126 = qJD(1) * t159;
t82 = cos(pkin(8));
t67 = t82 * pkin(1) + pkin(2);
t164 = -qJD(3) * t126 + t67 * qJDD(1);
t163 = qJDD(2) - g(1);
t52 = t67 * qJD(1);
t162 = qJD(3) * t52 + qJDD(1) * t159;
t78 = qJ(1) + pkin(8);
t70 = qJ(3) + t78;
t66 = cos(t70);
t60 = g(3) * t66;
t65 = sin(t70);
t61 = g(2) * t65;
t140 = t60 - t61;
t85 = sin(qJ(3));
t88 = cos(qJ(3));
t139 = t88 * t159 + t85 * t67;
t77 = qJD(1) + qJD(3);
t136 = qJ(5) * t77;
t31 = t88 * t126 + t85 * t52;
t23 = t77 * pkin(7) + t31;
t117 = t23 + t136;
t87 = cos(qJ(4));
t105 = t117 * t87;
t84 = sin(qJ(4));
t130 = t84 * qJD(4);
t123 = t77 * t130;
t154 = t87 * pkin(4);
t68 = pkin(3) + t154;
t76 = qJDD(1) + qJDD(3);
t150 = t68 * t76;
t161 = -pkin(4) * t123 + t150;
t120 = -t85 * t159 + t88 * t67;
t160 = -t162 * t88 - t164 * t85;
t79 = t84 ^ 2;
t158 = pkin(4) * t79;
t157 = g(1) * t87;
t156 = g(2) * t66;
t155 = t76 * pkin(3);
t153 = t31 * t77;
t33 = t139 * qJD(3);
t152 = t33 * t77;
t151 = t65 * t84;
t75 = t77 ^ 2;
t149 = t75 * t87;
t148 = t77 * t84;
t147 = t77 * t87;
t64 = t84 * t76;
t146 = t87 * t76;
t83 = -qJ(5) - pkin(7);
t72 = t87 * qJD(2);
t12 = -t117 * t84 + t72;
t135 = qJD(4) * pkin(4);
t11 = t12 + t135;
t144 = t11 - t12;
t30 = -t85 * t126 + t88 * t52;
t143 = t30 * t130 + t31 * t147;
t142 = t65 * t68 + t66 * t83;
t141 = g(3) * t151 + t84 * t156;
t80 = t87 ^ 2;
t138 = -t79 - t80;
t137 = t79 - t80;
t36 = pkin(7) + t139;
t134 = -qJ(5) - t36;
t132 = qJD(4) * t77;
t131 = qJDD(4) * pkin(4);
t129 = t87 * qJD(4);
t128 = qJD(4) * qJD(2);
t19 = t23 * t130;
t9 = t76 * pkin(7) - t160;
t98 = -qJ(5) * t76 - t128 - t9;
t94 = qJD(5) * t77 - t98;
t3 = -t19 + (-qJ(5) * t132 + qJDD(2)) * t84 + t94 * t87;
t127 = t3 * t87 + t140;
t125 = pkin(4) * t130;
t122 = -t65 * t83 + t66 * t68;
t16 = -t68 * t77 + qJD(5) - t30;
t114 = -t162 * t85 + t164 * t88;
t5 = qJDD(5) - t114 - t161;
t119 = t16 * t129 + t5 * t84 + t141;
t10 = -t114 - t155;
t22 = -t77 * pkin(3) - t30;
t118 = t10 * t84 + t22 * t129 + t141;
t116 = qJD(4) * t83;
t115 = 0.2e1 * t77 * t129;
t113 = qJD(4) * t134;
t35 = -pkin(3) - t120;
t69 = t87 * qJDD(2);
t111 = g(2) * t151 - t157 + t69;
t90 = qJD(4) ^ 2;
t110 = pkin(7) * t90 - t155;
t109 = g(3) * t65 + t156;
t86 = sin(qJ(1));
t89 = cos(qJ(1));
t108 = -g(2) * t89 - g(3) * t86;
t13 = t84 * qJD(2) + t105;
t107 = t11 * t84 - t13 * t87;
t28 = t33 + t125;
t29 = t35 - t154;
t106 = t28 * t77 + t29 * t76;
t104 = -t109 - t5;
t103 = -t10 - t109;
t102 = -t22 * t77 - t60 - t9;
t101 = -t163 * t84 + t87 * t61 + t19;
t100 = -pkin(3) * t132 - pkin(7) * qJDD(4);
t97 = t35 * t76 + t36 * t90 + t152;
t96 = t109 - t114;
t32 = t120 * qJD(3);
t95 = -qJDD(4) * t36 + (t35 * t77 - t32) * qJD(4);
t93 = -t60 + (-qJD(5) - t16) * t77 + t98;
t92 = -t140 + t160;
t73 = t87 * qJ(5);
t71 = t87 * qJD(5);
t55 = t87 * pkin(7) + t73;
t54 = t83 * t84;
t44 = qJDD(4) * t87 - t90 * t84;
t43 = qJDD(4) * t84 + t90 * t87;
t38 = -t84 * qJD(5) + t87 * t116;
t37 = t84 * t116 + t71;
t34 = t84 * t115 + t79 * t76;
t27 = t87 * t36 + t73;
t26 = t134 * t84;
t25 = t30 * t129;
t20 = -0.2e1 * t137 * t132 + 0.2e1 * t84 * t146;
t17 = t22 * t130;
t14 = t16 * t130;
t7 = (-qJD(5) - t32) * t84 + t87 * t113;
t6 = t84 * t113 + t87 * t32 + t71;
t2 = -qJD(4) * t105 - t94 * t84 + t131 + t69;
t1 = [qJDD(1), t108, g(2) * t86 - g(3) * t89, (t108 + (t81 ^ 2 + t82 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t76, t120 * t76 - t152 - t96, -t139 * t76 - t32 * t77 + t92, t34, t20, t43, t44, 0, t17 + t95 * t84 + (t103 - t97) * t87, t97 * t84 + t95 * t87 + t118, t26 * qJDD(4) + t14 + (t29 * t148 + t7) * qJD(4) + (t104 - t106) * t87, -t27 * qJDD(4) + t106 * t84 + (t29 * t147 - t6) * qJD(4) + t119, (t27 * t76 + t6 * t77 + (-t26 * t77 - t11) * qJD(4)) * t87 + (-t26 * t76 - t7 * t77 - t2 + (-t27 * t77 - t13) * qJD(4)) * t84 + t127, t3 * t27 + t13 * t6 + t2 * t26 + t11 * t7 + t5 * t29 + t16 * t28 - g(2) * (pkin(2) * cos(t78) + t89 * pkin(1) + t122) - g(3) * (pkin(2) * sin(t78) + t86 * pkin(1) + t142); 0, 0, 0, t163, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t43, t44, -t43, 0, -t107 * qJD(4) + t2 * t87 + t3 * t84 - g(1); 0, 0, 0, 0, t76, -t96 + t153, t30 * t77 + t92, t34, t20, t43, t44, 0, t17 + t100 * t84 + (t103 - t110) * t87 + t143, t25 + t100 * t87 + (t110 - t153) * t84 + t118, t54 * qJDD(4) + t14 + (-t68 * t148 + t38) * qJD(4) + (t104 + t161) * t87 + t143, -t55 * qJDD(4) + t25 + (-t150 - t153) * t84 + (-t37 + (-t68 * t87 + t158) * t77) * qJD(4) + t119, (-qJD(4) * t11 + t55 * t76) * t87 + (-t13 * qJD(4) - t54 * t76 - t2) * t84 + (t37 * t87 - t38 * t84 + t138 * t30 + (-t54 * t87 - t55 * t84) * qJD(4)) * t77 + t127, t3 * t55 + t2 * t54 - t5 * t68 - g(2) * t122 - g(3) * t142 + (-t31 + t125) * t16 + (-t87 * t30 + t37) * t13 + (t84 * t30 + t38) * t11; 0, 0, 0, 0, 0, 0, 0, -t84 * t149, t137 * t75, t64, t146, qJDD(4), t102 * t84 + t111, (-t84 * t23 + t72) * qJD(4) + (t102 - t128) * t87 + t101, 0.2e1 * t131 + (t13 - t105) * qJD(4) + (pkin(4) * t149 + t93) * t84 + t111, -t75 * t158 + (t84 * t136 + t12) * qJD(4) + t93 * t87 + t101, -pkin(4) * t64 + (-t135 + t144) * t147, t144 * t13 + (-t157 + t2 + (-t16 * t77 - t140) * t84) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t123 - t146, t64 + t115, t138 * t75, t107 * t77 - t104;];
tau_reg = t1;
