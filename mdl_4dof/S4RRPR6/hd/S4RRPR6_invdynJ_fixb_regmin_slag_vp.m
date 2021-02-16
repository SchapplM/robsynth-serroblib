% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRPR6
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [4x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:46
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:46:22
% EndTime: 2021-01-15 10:46:27
% DurationCPUTime: 1.00s
% Computational Cost: add. (1122->198), mult. (2702->281), div. (0->0), fcn. (1963->12), ass. (0->117)
t108 = sin(qJ(4));
t111 = cos(qJ(4));
t105 = sin(pkin(7));
t106 = cos(pkin(7));
t112 = cos(qJ(2));
t144 = t106 * t112;
t136 = qJD(1) * t144;
t109 = sin(qJ(2));
t142 = qJD(1) * t109;
t62 = t105 * t142 - t136;
t74 = t105 * t112 + t106 * t109;
t65 = t74 * qJD(1);
t124 = t108 * t62 - t111 * t65;
t138 = t112 * qJDD(1);
t139 = t109 * qJDD(1);
t128 = t105 * t139 - t106 * t138;
t64 = t74 * qJD(2);
t36 = qJD(1) * t64 + t128;
t140 = qJD(1) * qJD(2);
t135 = t109 * t140;
t119 = t74 * qJDD(1) - t105 * t135;
t134 = t112 * t140;
t37 = t106 * t134 + t119;
t116 = t124 * qJD(4) - t108 * t37 - t111 * t36;
t101 = qJD(2) + qJD(4);
t146 = t124 * t101;
t162 = t116 - t146;
t53 = t111 * t62;
t27 = -t108 * t65 - t53;
t145 = t27 * t101;
t141 = qJD(4) * t108;
t2 = -qJD(4) * t53 - t108 * t36 + t111 * t37 - t65 * t141;
t161 = t2 - t145;
t160 = t124 * t27;
t110 = sin(qJ(1));
t113 = cos(qJ(1));
t130 = g(1) * t113 + g(2) * t110;
t159 = t124 ^ 2 - t27 ^ 2;
t155 = t62 * pkin(6);
t107 = -qJ(3) - pkin(5);
t86 = t107 * t112;
t79 = qJD(1) * t86;
t147 = t106 * t79;
t148 = qJD(2) * pkin(2);
t85 = t107 * t109;
t78 = qJD(1) * t85;
t72 = t78 + t148;
t35 = t105 * t72 - t147;
t15 = t35 - t155;
t96 = t112 * pkin(2) + pkin(1);
t80 = -t96 * qJD(1) + qJD(3);
t43 = t62 * pkin(3) + t80;
t102 = qJ(2) + pkin(7);
t99 = qJ(4) + t102;
t93 = sin(t99);
t94 = cos(t99);
t158 = g(3) * t93 + t130 * t94 + t15 * t141 - t43 * t27;
t132 = qJD(2) * t107;
t60 = -t109 * qJD(3) + t112 * t132;
t33 = qJDD(2) * pkin(2) + t60 * qJD(1) + qJDD(1) * t85;
t59 = t112 * qJD(3) + t109 * t132;
t40 = t59 * qJD(1) - qJDD(1) * t86;
t6 = -t105 * t40 + t106 * t33;
t4 = qJDD(2) * pkin(3) - t37 * pkin(6) + t6;
t7 = t105 * t33 + t106 * t40;
t5 = -t36 * pkin(6) + t7;
t157 = -g(3) * t94 - t108 * t5 + t111 * t4 + t43 * t124 + t130 * t93;
t156 = qJD(4) - t101;
t154 = t65 * pkin(6);
t153 = pkin(2) * t105;
t152 = pkin(2) * t109;
t149 = g(3) * t112;
t19 = t105 * t60 + t106 * t59;
t68 = t105 * t79;
t42 = t106 * t78 + t68;
t45 = t105 * t85 - t106 * t86;
t103 = t109 ^ 2;
t143 = -t112 ^ 2 + t103;
t137 = t109 * t148;
t18 = -t105 * t59 + t106 * t60;
t34 = t106 * t72 + t68;
t41 = -t105 * t78 + t147;
t44 = t105 * t86 + t106 * t85;
t129 = g(1) * t110 - g(2) * t113;
t13 = qJD(2) * pkin(3) - t154 + t34;
t127 = -t108 * t13 - t111 * t15;
t20 = -t74 * pkin(6) + t44;
t73 = t105 * t109 - t144;
t21 = -t73 * pkin(6) + t45;
t126 = -t108 * t21 + t111 * t20;
t125 = t108 * t20 + t111 * t21;
t38 = t108 * t74 + t111 * t73;
t39 = -t108 * t73 + t111 * t74;
t95 = t106 * pkin(2) + pkin(3);
t123 = t108 * t95 + t111 * t153;
t122 = -t108 * t153 + t111 * t95;
t121 = -0.2e1 * pkin(1) * t140 - pkin(5) * qJDD(2);
t58 = pkin(2) * t135 - t96 * qJDD(1) + qJDD(3);
t114 = qJD(2) ^ 2;
t118 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t114 + t129;
t115 = qJD(1) ^ 2;
t117 = pkin(1) * t115 - pkin(5) * qJDD(1) + t130;
t100 = qJDD(2) + qJDD(4);
t98 = cos(t102);
t97 = sin(t102);
t67 = t73 * qJD(2);
t49 = t73 * pkin(3) - t96;
t47 = t64 * pkin(3) + t137;
t46 = pkin(2) * t142 + t65 * pkin(3);
t17 = t42 - t154;
t16 = t41 + t155;
t14 = t36 * pkin(3) + t58;
t11 = -t64 * pkin(6) + t19;
t10 = t67 * pkin(6) + t18;
t9 = t39 * qJD(4) - t108 * t67 + t111 * t64;
t8 = -t38 * qJD(4) - t108 * t64 - t111 * t67;
t1 = [qJDD(1), t129, t130, t103 * qJDD(1) + 0.2e1 * t109 * t134, 0.2e1 * t109 * t138 - 0.2e1 * t143 * t140, qJDD(2) * t109 + t114 * t112, qJDD(2) * t112 - t114 * t109, 0, t121 * t109 + t118 * t112, -t118 * t109 + t121 * t112, t44 * qJDD(2) - t96 * t36 + t58 * t73 + t80 * t64 + t129 * t98 + (t62 * t152 + t18) * qJD(2), -t45 * qJDD(2) - t96 * t37 + t58 * t74 - t80 * t67 - t129 * t97 + (t65 * t152 - t19) * qJD(2), -t18 * t65 - t19 * t62 + t34 * t67 - t35 * t64 - t45 * t36 - t44 * t37 - t6 * t74 - t7 * t73 - t130, t7 * t45 + t35 * t19 + t6 * t44 + t34 * t18 - t58 * t96 + t80 * t137 - g(1) * (-t107 * t113 - t110 * t96) - g(2) * (-t110 * t107 + t113 * t96), -t124 * t8 + t2 * t39, t116 * t39 + t124 * t9 - t2 * t38 + t27 * t8, t39 * t100 + t8 * t101, -t38 * t100 - t9 * t101, 0, -t47 * t27 - t49 * t116 + t14 * t38 + t43 * t9 + (-t125 * qJD(4) + t111 * t10 - t108 * t11) * t101 + t126 * t100 + t129 * t94, -t47 * t124 + t49 * t2 + t14 * t39 + t43 * t8 - (t126 * qJD(4) + t108 * t10 + t111 * t11) * t101 - t125 * t100 - t129 * t93; 0, 0, 0, -t109 * t115 * t112, t143 * t115, t139, t138, qJDD(2), t109 * t117 - t149, g(3) * t109 + t112 * t117, -g(3) * t98 - t41 * qJD(2) - t80 * t65 + t130 * t97 + (qJDD(2) * t106 - t62 * t142) * pkin(2) + t6, g(3) * t97 + t42 * qJD(2) + t80 * t62 + t130 * t98 + (-qJDD(2) * t105 - t65 * t142) * pkin(2) - t7, (t35 + t41) * t65 + (-t34 + t42) * t62 + (-t105 * t36 - t106 * t37) * pkin(2), -t34 * t41 - t35 * t42 + (-t149 + t105 * t7 + t106 * t6 + (-qJD(1) * t80 + t130) * t109) * pkin(2), t160, t159, t161, t162, t100, t122 * t100 + t46 * t27 - (-t108 * t17 + t111 * t16) * t101 + (-t123 * t101 + t127) * qJD(4) + t157, -t123 * t100 - t111 * t5 - t108 * t4 + t46 * t124 + (t108 * t16 + t111 * t17) * t101 + (-t122 * t101 - t111 * t13) * qJD(4) + t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t65 * qJD(2) + t128, (-t62 + t136) * qJD(2) + t119, -t62 ^ 2 - t65 ^ 2, t34 * t65 + t35 * t62 - t129 + t58, 0, 0, 0, 0, 0, -t116 - t146, t2 + t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t160, t159, t161, t162, t100, t156 * t127 + t157, (-t15 * t101 - t4) * t108 + (-t156 * t13 - t5) * t111 + t158;];
tau_reg = t1;
