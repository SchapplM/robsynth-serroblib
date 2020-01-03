% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPRR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:28:23
% EndTime: 2020-01-03 11:28:30
% DurationCPUTime: 0.95s
% Computational Cost: add. (1216->182), mult. (2601->247), div. (0->0), fcn. (2041->16), ass. (0->116)
t101 = sin(qJ(5));
t104 = cos(qJ(5));
t105 = cos(qJ(4));
t99 = cos(pkin(9));
t142 = t105 * t99;
t102 = sin(qJ(4));
t97 = sin(pkin(9));
t144 = t102 * t97;
t64 = -t142 + t144;
t56 = t64 * qJD(1);
t65 = t102 * t99 + t105 * t97;
t57 = t65 * qJD(1);
t19 = t101 * t57 + t104 * t56;
t95 = qJD(4) + qJD(5);
t148 = t19 * t95;
t138 = qJD(5) * t104;
t139 = qJD(5) * t101;
t131 = qJD(1) * t144;
t133 = qJDD(1) * t105;
t134 = qJDD(1) * t102;
t140 = qJD(4) * t105;
t141 = qJD(1) * t99;
t132 = t97 * t133 + t99 * t134 + t140 * t141;
t27 = -qJD(4) * t131 + t132;
t124 = -t99 * t133 + t97 * t134;
t59 = t65 * qJD(4);
t28 = qJD(1) * t59 + t124;
t5 = -t101 * t28 + t104 * t27 - t56 * t138 - t57 * t139;
t162 = t5 + t148;
t117 = -t101 * t56 + t104 * t57;
t161 = t117 * t19;
t149 = t117 * t95;
t6 = t117 * qJD(5) + t101 * t27 + t104 * t28;
t160 = -t6 + t149;
t96 = qJ(1) + pkin(8);
t87 = sin(t96);
t89 = cos(t96);
t126 = g(2) * t87 - g(3) * t89;
t159 = t117 ^ 2 - t19 ^ 2;
t98 = sin(pkin(8));
t77 = t98 * pkin(1) + qJ(3);
t70 = t77 * qJD(1);
t85 = t99 * qJD(2);
t38 = t85 + (-pkin(6) * qJD(1) - t70) * t97;
t45 = t97 * qJD(2) + t99 * t70;
t39 = pkin(6) * t141 + t45;
t115 = -t102 * t38 - t105 * t39;
t10 = -t56 * pkin(7) - t115;
t62 = qJD(1) * qJD(3) + t77 * qJDD(1);
t83 = t99 * qJDD(2);
t34 = t83 + (-pkin(6) * qJDD(1) - t62) * t97;
t136 = t99 * qJDD(1);
t41 = t97 * qJDD(2) + t99 * t62;
t35 = pkin(6) * t136 + t41;
t129 = -t102 * t35 + t105 * t34;
t2 = qJDD(4) * pkin(4) - t27 * pkin(7) + t115 * qJD(4) + t129;
t100 = cos(pkin(8));
t81 = -t100 * pkin(1) - pkin(2);
t69 = -t99 * pkin(3) + t81;
t54 = t69 * qJD(1) + qJD(3);
t26 = t56 * pkin(4) + t54;
t94 = pkin(9) + qJ(4);
t90 = qJ(5) + t94;
t79 = sin(t90);
t80 = cos(t90);
t158 = t26 * t19 + t10 * t139 + g(1) * t79 + (-t10 * t95 - t2) * t101 + t126 * t80;
t127 = g(2) * t89 + g(3) * t87;
t137 = qJDD(1) * t81;
t68 = qJDD(3) + t137;
t156 = t127 + t68;
t155 = -t102 * t39 + t105 * t38;
t143 = t104 * t10;
t9 = -t57 * pkin(7) + t155;
t8 = qJD(4) * pkin(4) + t9;
t120 = -t101 * t8 - t143;
t116 = t102 * t34 + t105 * t35;
t3 = -t28 * pkin(7) + t155 * qJD(4) + t116;
t154 = -g(1) * t80 + t120 * qJD(5) - t101 * t3 + t104 * t2 - t26 * t117 + t126 * t79;
t153 = pkin(4) * t59;
t150 = pkin(6) + t77;
t60 = t150 * t97;
t61 = t150 * t99;
t147 = -t102 * t60 + t105 * t61;
t146 = t97 ^ 2 + t99 ^ 2;
t128 = -t102 * t61 - t105 * t60;
t103 = sin(qJ(1));
t106 = cos(qJ(1));
t125 = -g(2) * t106 - g(3) * t103;
t31 = t101 * t65 + t104 * t64;
t58 = t64 * qJD(4);
t11 = -t31 * qJD(5) - t101 * t59 - t104 * t58;
t32 = -t101 * t64 + t104 * t65;
t91 = qJDD(4) + qJDD(5);
t123 = t11 * t95 + t32 * t91;
t40 = -t97 * t62 + t83;
t122 = -t40 * t97 + t41 * t99;
t121 = (-t97 * t70 + t85) * t97 - t45 * t99;
t16 = -t65 * pkin(7) + t128;
t17 = -t64 * pkin(7) + t147;
t119 = -t101 * t17 + t104 * t16;
t118 = t101 * t16 + t104 * t17;
t114 = t137 + t156;
t52 = t69 * qJDD(1) + qJDD(3);
t113 = -t60 * t140 + qJD(3) * t142 + (-qJD(3) * t97 - qJD(4) * t61) * t102;
t110 = -t65 * qJD(3) - t147 * qJD(4);
t88 = cos(t94);
t86 = sin(t94);
t37 = t64 * pkin(4) + t69;
t30 = -t59 * qJD(4) - t64 * qJDD(4);
t29 = -t58 * qJD(4) + t65 * qJDD(4);
t15 = t28 * pkin(4) + t52;
t14 = t58 * pkin(7) + t110;
t13 = -t59 * pkin(7) + t113;
t12 = t32 * qJD(5) - t101 * t58 + t104 * t59;
t4 = -t12 * t95 - t31 * t91;
t1 = [qJDD(1), t125, g(2) * t103 - g(3) * t106, (t125 + (t100 ^ 2 + t98 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), -t114 * t99, t114 * t97, t62 * t146 + t122 - t126, t68 * t81 - g(2) * (t106 * pkin(1) + t89 * pkin(2) + t87 * qJ(3)) - g(3) * (t103 * pkin(1) + t87 * pkin(2) - t89 * qJ(3)) + t122 * t77 - t121 * qJD(3), t27 * t65 - t57 * t58, -t27 * t64 - t65 * t28 + t58 * t56 - t57 * t59, t29, t30, 0, t110 * qJD(4) + t128 * qJDD(4) - t127 * t88 + t69 * t28 + t52 * t64 + t54 * t59, -t113 * qJD(4) - t147 * qJDD(4) + t127 * t86 + t69 * t27 + t52 * t65 - t54 * t58, t11 * t117 + t5 * t32, -t11 * t19 - t117 * t12 - t5 * t31 - t32 * t6, t123, t4, 0, t19 * t153 + t37 * t6 + t15 * t31 + t26 * t12 + (-t118 * qJD(5) - t101 * t13 + t104 * t14) * t95 + t119 * t91 - t127 * t80, t117 * t153 + t37 * t5 + t15 * t32 + t26 * t11 - (t119 * qJD(5) + t101 * t14 + t104 * t13) * t95 - t118 * t91 + t127 * t79; 0, 0, 0, qJDD(2) - g(1), 0, 0, 0, t40 * t99 + t41 * t97 - g(1), 0, 0, 0, 0, 0, t30, -t29, 0, 0, 0, 0, 0, t4, -t123; 0, 0, 0, 0, -t136, t97 * qJDD(1), -t146 * qJD(1) ^ 2, t121 * qJD(1) + t156, 0, 0, 0, 0, 0, 0.2e1 * t57 * qJD(4) + t124, (-t56 - t131) * qJD(4) + t132, 0, 0, 0, 0, 0, t6 + t149, t5 - t148; 0, 0, 0, 0, 0, 0, 0, 0, t57 * t56, -t56 ^ 2 + t57 ^ 2, (t56 - t131) * qJD(4) + t132, -t124, qJDD(4), -g(1) * t88 + t126 * t86 - t54 * t57 + t129, g(1) * t86 + t126 * t88 + t54 * t56 - t116, t161, t159, t162, t160, t91, -(-t101 * t9 - t143) * t95 + (t104 * t91 - t95 * t139 - t57 * t19) * pkin(4) + t154, (-qJD(5) * t8 + t9 * t95 - t3) * t104 + (-t101 * t91 - t117 * t57 - t95 * t138) * pkin(4) + t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t161, t159, t162, t160, t91, -t120 * t95 + t154, (-t3 + (-qJD(5) + t95) * t8) * t104 + t158;];
tau_reg = t1;
