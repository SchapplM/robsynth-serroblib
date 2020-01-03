% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRP6
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
% tau_reg [5x20]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:43:18
% EndTime: 2019-12-31 18:43:22
% DurationCPUTime: 1.43s
% Computational Cost: add. (1516->256), mult. (3161->348), div. (0->0), fcn. (1960->10), ass. (0->147)
t89 = sin(pkin(8));
t74 = t89 * pkin(1) + pkin(6);
t65 = t74 * qJDD(1);
t195 = -qJD(2) * qJD(3) - t65;
t67 = t74 * qJD(1);
t93 = sin(qJ(3));
t96 = cos(qJ(3));
t37 = t96 * qJD(2) - t93 * t67;
t194 = t37 * qJD(3);
t86 = qJ(1) + pkin(8);
t79 = sin(t86);
t80 = cos(t86);
t123 = g(1) * t80 + g(2) * t79;
t111 = t123 * t93;
t155 = qJD(3) * t96;
t138 = -t67 * t155 + t195 * t93;
t116 = -qJDD(3) * pkin(3) - t138;
t143 = t96 * qJDD(2);
t16 = t116 - t143;
t147 = t96 * qJD(1);
t73 = -qJD(4) + t147;
t193 = -qJD(4) * pkin(7) * t73 + g(3) * t96 - t111 + t16;
t152 = qJD(4) * t93;
t192 = qJD(1) * t152 - qJDD(3);
t144 = t93 * qJDD(1);
t92 = sin(qJ(4));
t95 = cos(qJ(4));
t23 = ((qJD(4) + t147) * qJD(3) + t144) * t92 + t192 * t95;
t145 = qJDD(2) - g(3);
t191 = t145 * t96;
t173 = t92 * t96;
t32 = t79 * t173 + t80 * t95;
t34 = -t80 * t173 + t79 * t95;
t190 = -g(1) * t34 + g(2) * t32;
t149 = t95 * qJD(3);
t135 = t96 * t149;
t107 = -t92 * t152 + t135;
t142 = qJD(1) * qJD(3);
t131 = t93 * t142;
t83 = t96 * qJDD(1);
t52 = qJDD(4) - t83 + t131;
t171 = t95 * t52;
t189 = -t107 * t73 + t93 * t171;
t150 = t92 * qJD(3);
t159 = qJD(1) * t93;
t58 = t95 * t159 + t150;
t187 = t58 ^ 2;
t179 = t73 * pkin(4);
t38 = t93 * qJD(2) + t96 * t67;
t29 = qJD(3) * pkin(7) + t38;
t90 = cos(pkin(8));
t75 = -t90 * pkin(1) - pkin(2);
t48 = -t96 * pkin(3) - t93 * pkin(7) + t75;
t30 = t48 * qJD(1);
t10 = -t92 * t29 + t95 * t30;
t7 = -t58 * qJ(5) + t10;
t6 = t7 - t179;
t186 = -t7 + t6;
t185 = pkin(4) * t92;
t180 = g(3) * t93;
t160 = qJ(5) * t95;
t114 = pkin(4) * t93 - t96 * t160;
t169 = qJ(5) + pkin(7);
t128 = qJD(4) * t169;
t124 = pkin(3) * t93 - pkin(7) * t96;
t61 = t124 * qJD(1);
t46 = t95 * t61;
t178 = -t114 * qJD(1) - t95 * t128 - t46 + (-qJD(5) + t37) * t92;
t56 = t92 * t159 - t149;
t177 = t56 * t73;
t176 = t58 * t73;
t175 = t58 * t93;
t174 = t73 * t95;
t172 = t93 * t95;
t170 = t95 * t96;
t168 = -t56 * t135 - t23 * t172;
t167 = t95 * t37 + t92 * t61;
t151 = qJD(4) * t95;
t62 = t124 * qJD(3);
t166 = t48 * t151 + t92 * t62;
t148 = t95 * qJD(5);
t165 = t148 - t167 + (qJ(5) * t147 - t128) * t92;
t164 = t93 * t74 * t150 + t95 * t62;
t60 = t74 * t170;
t163 = t92 * t48 + t60;
t87 = t93 ^ 2;
t162 = -t96 ^ 2 + t87;
t161 = qJ(5) * t93;
t68 = qJD(1) * t75;
t158 = qJD(3) * t56;
t157 = qJD(3) * t74;
t156 = qJD(3) * t93;
t154 = qJD(4) * t56;
t153 = qJD(4) * t92;
t15 = qJDD(3) * pkin(7) + t93 * qJDD(2) + t96 * t65 + t194;
t24 = qJD(1) * t62 + t48 * qJDD(1);
t139 = t95 * t15 + t30 * t151 + t92 * t24;
t137 = t73 * t150;
t136 = t58 * t155;
t134 = pkin(6) + t185;
t133 = t74 + t185;
t22 = -qJD(1) * t135 - qJD(4) * t149 - t95 * t144 + t192 * t92;
t130 = t58 * t156 + t22 * t96;
t129 = t73 * t74 + t29;
t127 = -qJD(4) * t30 - t15;
t126 = t151 * t175;
t122 = g(1) * t79 - g(2) * t80;
t94 = sin(qJ(1));
t97 = cos(qJ(1));
t121 = g(1) * t94 - g(2) * t97;
t11 = t95 * t29 + t92 * t30;
t8 = -t56 * qJ(5) + t11;
t120 = t6 * t95 + t8 * t92;
t119 = t6 * t92 - t8 * t95;
t78 = t95 * pkin(4) + pkin(3);
t118 = t169 * t93 + t96 * t78;
t115 = pkin(2) + t118;
t112 = t121 * pkin(1);
t110 = t73 * t151 - t92 * t52;
t28 = -qJD(3) * pkin(3) - t37;
t109 = -t29 * t153 + t139;
t106 = -qJD(1) * t68 + t123;
t105 = -pkin(7) * t52 - t73 * t28;
t104 = 0.2e1 * t68 * qJD(3) - qJDD(3) * t74;
t103 = t23 * pkin(4) + qJDD(5) + t116;
t98 = qJD(3) ^ 2;
t102 = -0.2e1 * qJDD(1) * t75 - t74 * t98 + t122;
t19 = t95 * t24;
t1 = t52 * pkin(4) + t22 * qJ(5) - t11 * qJD(4) - t58 * qJD(5) - t92 * t15 + t19;
t2 = -t23 * qJ(5) - t56 * qJD(5) + t109;
t101 = -t120 * qJD(4) - t1 * t92 + t2 * t95;
t99 = qJD(1) ^ 2;
t70 = t169 * t95;
t69 = t169 * t92;
t64 = qJDD(3) * t96 - t98 * t93;
t63 = qJDD(3) * t93 + t98 * t96;
t51 = t56 ^ 2;
t40 = t95 * t48;
t35 = t80 * t170 + t79 * t92;
t33 = -t79 * t170 + t80 * t92;
t21 = t56 * pkin(4) + qJD(5) + t28;
t20 = -t92 * t161 + t163;
t17 = -t93 * t160 + t40 + (-t74 * t92 - pkin(4)) * t96;
t5 = t103 - t143;
t4 = (-qJ(5) * qJD(4) - t157) * t172 + (-qJD(5) * t93 + (-qJ(5) * qJD(3) - qJD(4) * t74) * t96) * t92 + t166;
t3 = -t93 * t148 + t114 * qJD(3) + (-t60 + (-t48 + t161) * t92) * qJD(4) + t164;
t9 = [qJDD(1), t121, g(1) * t97 + g(2) * t94, (t89 ^ 2 + t90 ^ 2) * pkin(1) ^ 2 * qJDD(1) + t112, t87 * qJDD(1) + 0.2e1 * t96 * t131, -0.2e1 * t162 * t142 + 0.2e1 * t93 * t83, t63, t64, 0, t102 * t96 + t104 * t93, -t102 * t93 + t104 * t96, t107 * t58 - t22 * t172, -t126 + (-t136 + (t22 + t154) * t93) * t92 + t168, t130 + t189, (t23 + t137) * t96 + (t110 - t158) * t93, -t73 * t156 - t52 * t96, -(-t48 * t153 + t164) * t73 + t40 * t52 - g(1) * t33 - g(2) * t35 + (t56 * t157 - t19 + t129 * t151 + (qJD(3) * t28 - t52 * t74 - t127) * t92) * t96 + (t10 * qJD(3) + t28 * t151 + t16 * t92 + t74 * t23) * t93, t166 * t73 - t163 * t52 - g(1) * t32 - g(2) * t34 + (-t129 * t153 + (t28 * t95 + t58 * t74) * qJD(3) + t139) * t96 + (-t28 * t153 + t16 * t95 - t74 * t22 + (-t74 * t174 - t11) * qJD(3)) * t93, t17 * t22 - t20 * t23 - t3 * t58 - t4 * t56 - t120 * t155 + (qJD(4) * t119 - t1 * t95 - t2 * t92 + t122) * t93, t1 * t17 + t2 * t20 + t6 * t3 + t8 * t4 + t21 * t133 * t155 + (t21 * pkin(4) * t151 + t133 * t5) * t93 + t112 + (-g(1) * t134 - g(2) * t115) * t80 + (g(1) * t115 - g(2) * t134) * t79; 0, 0, 0, t145, 0, 0, 0, 0, 0, t64, -t63, 0, 0, 0, 0, 0, (-t23 + t137) * t96 + (t110 + t158) * t93, t130 - t189, t126 + (t136 + (-t22 + t154) * t93) * t92 + t168, -g(3) + (-t119 * qJD(3) - t5) * t96 + (qJD(3) * t21 + t101) * t93; 0, 0, 0, 0, -t93 * t99 * t96, t162 * t99, t144, t83, qJDD(3), t38 * qJD(3) + t106 * t93 + t138 + t191, t194 + (qJD(3) * t67 - t145) * t93 + (t106 + t195) * t96, -t58 * t174 - t22 * t92, (-t22 + t177) * t95 + (-t23 + t176) * t92, (t73 * t170 - t175) * qJD(1) - t110, t73 * t153 + t171 + (-t73 * t173 + t56 * t93) * qJD(1), t73 * t159, -t10 * t159 - pkin(3) * t23 - t38 * t56 + t46 * t73 + (-t37 * t73 + t105) * t92 - t193 * t95, pkin(3) * t22 + t105 * t95 + t11 * t159 - t167 * t73 + t193 * t92 - t38 * t58, -t180 - t69 * t22 - t70 * t23 - t178 * t58 - t165 * t56 + (qJD(1) * t120 - t123) * t96 + t101, t2 * t70 - t1 * t69 - t5 * t78 - g(3) * t118 + t165 * t8 + t178 * t6 + (-t92 * t179 - t38) * t21 + t123 * (-t169 * t96 + t78 * t93); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58 * t56, -t51 + t187, -t22 - t177, -t176 - t23, t52, -t29 * t151 - t11 * t73 - t28 * t58 + t19 + (t127 + t180) * t92 + t190, g(1) * t35 - g(2) * t33 + g(3) * t172 - t10 * t73 + t28 * t56 - t109, pkin(4) * t22 - t186 * t56, t186 * t8 + (t92 * t180 - t21 * t58 + t1 + t190) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51 - t187, t8 * t56 + t6 * t58 + t103 - t111 - t191;];
tau_reg = t9;
