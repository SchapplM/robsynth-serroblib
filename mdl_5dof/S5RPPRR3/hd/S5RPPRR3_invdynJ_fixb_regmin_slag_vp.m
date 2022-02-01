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
% tau_reg [5x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:15
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:14:36
% EndTime: 2022-01-23 09:14:38
% DurationCPUTime: 0.92s
% Computational Cost: add. (1208->182), mult. (2589->244), div. (0->0), fcn. (2032->16), ass. (0->114)
t103 = sin(qJ(5));
t106 = cos(qJ(5));
t101 = cos(pkin(9));
t107 = cos(qJ(4));
t143 = t107 * t101;
t104 = sin(qJ(4));
t99 = sin(pkin(9));
t145 = t104 * t99;
t116 = t143 - t145;
t56 = t116 * qJD(1);
t65 = t104 * t101 + t107 * t99;
t57 = t65 * qJD(1);
t19 = t103 * t57 - t106 * t56;
t97 = qJD(4) + qJD(5);
t149 = t19 * t97;
t139 = qJD(5) * t106;
t140 = qJD(5) * t103;
t133 = qJD(1) * t145;
t135 = t101 * qJDD(1);
t137 = qJDD(1) * t99;
t141 = qJD(4) * t107;
t142 = qJD(1) * t101;
t134 = t104 * t135 + t107 * t137 + t141 * t142;
t27 = -qJD(4) * t133 + t134;
t126 = t104 * t137 - t107 * t135;
t59 = t65 * qJD(4);
t28 = qJD(1) * t59 + t126;
t5 = -t103 * t28 + t106 * t27 + t56 * t139 - t57 * t140;
t163 = t5 + t149;
t119 = t103 * t56 + t106 * t57;
t162 = t119 * t19;
t150 = t119 * t97;
t6 = t119 * qJD(5) + t103 * t27 + t106 * t28;
t161 = -t6 + t150;
t98 = qJ(1) + pkin(8);
t89 = sin(t98);
t91 = cos(t98);
t129 = g(1) * t91 + g(2) * t89;
t160 = t119 ^ 2 - t19 ^ 2;
t100 = sin(pkin(8));
t79 = t100 * pkin(1) + qJ(3);
t72 = t79 * qJD(1);
t87 = t101 * qJD(2);
t38 = t87 + (-pkin(6) * qJD(1) - t72) * t99;
t45 = t99 * qJD(2) + t101 * t72;
t39 = pkin(6) * t142 + t45;
t117 = -t104 * t38 - t107 * t39;
t10 = t56 * pkin(7) - t117;
t62 = qJD(1) * qJD(3) + t79 * qJDD(1);
t85 = t101 * qJDD(2);
t34 = t85 + (-pkin(6) * qJDD(1) - t62) * t99;
t41 = t99 * qJDD(2) + t101 * t62;
t35 = pkin(6) * t135 + t41;
t131 = -t104 * t35 + t107 * t34;
t2 = qJDD(4) * pkin(4) - t27 * pkin(7) + t117 * qJD(4) + t131;
t102 = cos(pkin(8));
t83 = -t102 * pkin(1) - pkin(2);
t71 = -t101 * pkin(3) + t83;
t54 = t71 * qJD(1) + qJD(3);
t26 = -t56 * pkin(4) + t54;
t96 = pkin(9) + qJ(4);
t92 = qJ(5) + t96;
t81 = sin(t92);
t82 = cos(t92);
t159 = t26 * t19 + t10 * t140 + g(3) * t81 + (-t10 * t97 - t2) * t103 + t129 * t82;
t128 = g(1) * t89 - g(2) * t91;
t138 = qJDD(1) * t83;
t70 = qJDD(3) + t138;
t157 = t128 - t70;
t156 = -t104 * t39 + t107 * t38;
t144 = t106 * t10;
t9 = -t57 * pkin(7) + t156;
t8 = qJD(4) * pkin(4) + t9;
t124 = -t103 * t8 - t144;
t118 = t104 * t34 + t107 * t35;
t3 = -t28 * pkin(7) + t156 * qJD(4) + t118;
t155 = -g(3) * t82 + t124 * qJD(5) - t103 * t3 + t106 * t2 - t26 * t119 + t129 * t81;
t154 = pkin(4) * t59;
t151 = pkin(6) + t79;
t60 = t151 * t99;
t61 = t151 * t101;
t148 = -t104 * t60 + t107 * t61;
t147 = t101 ^ 2 + t99 ^ 2;
t130 = -t104 * t61 - t107 * t60;
t105 = sin(qJ(1));
t108 = cos(qJ(1));
t127 = g(1) * t105 - g(2) * t108;
t31 = t103 * t65 - t106 * t116;
t58 = t116 * qJD(4);
t11 = -t31 * qJD(5) - t103 * t59 + t106 * t58;
t32 = t103 * t116 + t106 * t65;
t93 = qJDD(4) + qJDD(5);
t125 = t11 * t97 + t32 * t93;
t40 = -t99 * t62 + t85;
t123 = t41 * t101 - t40 * t99;
t122 = t45 * t101 - (-t99 * t72 + t87) * t99;
t16 = -t65 * pkin(7) + t130;
t17 = pkin(7) * t116 + t148;
t121 = -t103 * t17 + t106 * t16;
t120 = t103 * t16 + t106 * t17;
t52 = t71 * qJDD(1) + qJDD(3);
t114 = -t60 * t141 + qJD(3) * t143 + (-qJD(3) * t99 - qJD(4) * t61) * t104;
t111 = -t65 * qJD(3) - t148 * qJD(4);
t90 = cos(t96);
t88 = sin(t96);
t37 = -pkin(4) * t116 + t71;
t30 = -t59 * qJD(4) + qJDD(4) * t116;
t29 = t58 * qJD(4) + t65 * qJDD(4);
t15 = t28 * pkin(4) + t52;
t14 = -t58 * pkin(7) + t111;
t13 = -t59 * pkin(7) + t114;
t12 = t32 * qJD(5) + t103 * t58 + t106 * t59;
t4 = -t12 * t97 - t31 * t93;
t1 = [qJDD(1), t127, g(1) * t108 + g(2) * t105, (t127 + (t100 ^ 2 + t102 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), (-t138 + t157) * t101, t62 * t147 + t123 - t129, t70 * t83 - g(1) * (-t105 * pkin(1) - t89 * pkin(2) + t91 * qJ(3)) - g(2) * (t108 * pkin(1) + t91 * pkin(2) + t89 * qJ(3)) + t123 * t79 + t122 * qJD(3), t27 * t65 + t57 * t58, t116 * t27 - t65 * t28 + t58 * t56 - t57 * t59, t29, t30, 0, t111 * qJD(4) + t130 * qJDD(4) - t116 * t52 + t128 * t90 + t71 * t28 + t54 * t59, -t114 * qJD(4) - t148 * qJDD(4) - t128 * t88 + t71 * t27 + t52 * t65 + t54 * t58, t11 * t119 + t5 * t32, -t11 * t19 - t119 * t12 - t5 * t31 - t32 * t6, t125, t4, 0, t19 * t154 + t37 * t6 + t15 * t31 + t26 * t12 + (-t120 * qJD(5) - t103 * t13 + t106 * t14) * t97 + t121 * t93 + t128 * t82, t119 * t154 + t37 * t5 + t15 * t32 + t26 * t11 - (t121 * qJD(5) + t103 * t14 + t106 * t13) * t97 - t120 * t93 - t128 * t81; 0, 0, 0, qJDD(2) - g(3), 0, 0, t40 * t101 + t41 * t99 - g(3), 0, 0, 0, 0, 0, t30, -t29, 0, 0, 0, 0, 0, t4, -t125; 0, 0, 0, 0, -t135, -t147 * qJD(1) ^ 2, -t122 * qJD(1) - t157, 0, 0, 0, 0, 0, 0.2e1 * t57 * qJD(4) + t126, (t56 - t133) * qJD(4) + t134, 0, 0, 0, 0, 0, t6 + t150, t5 - t149; 0, 0, 0, 0, 0, 0, 0, -t57 * t56, -t56 ^ 2 + t57 ^ 2, (-t56 - t133) * qJD(4) + t134, -t126, qJDD(4), -g(3) * t90 + t129 * t88 - t54 * t57 + t131, g(3) * t88 + t129 * t90 - t54 * t56 - t118, t162, t160, t163, t161, t93, -(-t103 * t9 - t144) * t97 + (t106 * t93 - t97 * t140 - t57 * t19) * pkin(4) + t155, (-qJD(5) * t8 + t9 * t97 - t3) * t106 + (-t103 * t93 - t119 * t57 - t97 * t139) * pkin(4) + t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t162, t160, t163, t161, t93, -t124 * t97 + t155, (-t3 + (-qJD(5) + t97) * t8) * t106 + t159;];
tau_reg = t1;
