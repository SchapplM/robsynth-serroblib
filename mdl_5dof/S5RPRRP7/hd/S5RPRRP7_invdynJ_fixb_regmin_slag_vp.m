% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRP7
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
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP7_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:45:39
% EndTime: 2019-12-31 18:45:44
% DurationCPUTime: 1.94s
% Computational Cost: add. (2078->310), mult. (4190->396), div. (0->0), fcn. (2602->10), ass. (0->155)
t96 = sin(pkin(8));
t82 = pkin(1) * t96 + pkin(6);
t71 = t82 * qJDD(1);
t207 = -qJD(2) * qJD(3) - t71;
t102 = cos(qJ(3));
t73 = t82 * qJD(1);
t99 = sin(qJ(3));
t44 = qJD(2) * t102 - t99 * t73;
t206 = qJD(3) * t44;
t158 = qJD(1) * t102;
t205 = qJD(4) - t158;
t162 = qJD(4) * t99;
t204 = qJD(1) * t162 - qJDD(3);
t101 = cos(qJ(4));
t153 = t99 * qJDD(1);
t98 = sin(qJ(4));
t23 = t204 * t101 + ((qJD(4) + t158) * qJD(3) + t153) * t98;
t161 = t101 * t102;
t124 = t102 * pkin(3) + pkin(7) * t99 + pkin(2);
t97 = cos(pkin(8));
t199 = pkin(1) * t97;
t52 = -t124 - t199;
t180 = t82 * t161 + t98 * t52;
t93 = qJ(1) + pkin(8);
t87 = cos(t93);
t196 = g(1) * t87;
t86 = sin(t93);
t132 = g(2) * t86 + t196;
t203 = t132 * t99;
t34 = -qJD(3) * pkin(3) - t44;
t154 = t101 * qJD(3);
t168 = qJD(1) * t99;
t59 = t98 * t168 - t154;
t166 = qJD(3) * t98;
t61 = t101 * t168 + t166;
t10 = pkin(4) * t59 - qJ(5) * t61 + t34;
t152 = qJD(1) * qJD(3);
t90 = t102 * qJDD(1);
t55 = t99 * t152 + qJDD(4) - t90;
t197 = pkin(7) * t55;
t202 = -t10 * t205 + t197;
t200 = t61 ^ 2;
t198 = pkin(4) * t55;
t194 = g(3) * t99;
t45 = t99 * qJD(2) + t102 * t73;
t35 = qJD(3) * pkin(7) + t45;
t36 = t52 * qJD(1);
t12 = t101 * t35 + t36 * t98;
t7 = qJ(5) * t205 + t12;
t193 = t7 * t205;
t192 = g(3) * t102;
t191 = t12 * t205;
t190 = t59 * t205;
t189 = t61 * t59;
t188 = t61 * t205;
t187 = t61 * t99;
t186 = t82 * t98;
t142 = t102 * t154;
t172 = t101 * t99;
t185 = -t59 * t142 - t23 * t172;
t130 = pkin(3) * t99 - pkin(7) * t102;
t64 = t130 * qJD(1);
t184 = t101 * t44 + t98 * t64;
t183 = t142 * t205 + t55 * t172;
t155 = qJD(4) * t101;
t65 = t130 * qJD(3);
t182 = t52 * t155 + t98 * t65;
t126 = pkin(4) * t98 - qJ(5) * t101;
t181 = -qJD(5) * t98 + t205 * t126 - t45;
t169 = t86 * t101;
t179 = g(2) * t99 * t169 + t172 * t196;
t94 = t99 ^ 2;
t178 = -t102 ^ 2 + t94;
t177 = pkin(7) * qJD(4);
t176 = qJ(5) * t55;
t175 = t101 * t55;
t174 = t101 * t65;
t173 = t101 * t82;
t141 = t102 * t152;
t22 = -qJD(4) * t154 + t204 * t98 + (-t141 - t153) * t101;
t171 = t102 * t22;
t170 = t102 * t98;
t83 = -pkin(2) - t199;
t74 = qJD(1) * t83;
t167 = qJD(3) * t59;
t165 = qJD(3) * t99;
t164 = qJD(4) * t59;
t163 = qJD(4) * t98;
t11 = t101 * t36 - t35 * t98;
t160 = qJD(5) - t11;
t159 = qJDD(2) - g(3);
t156 = qJD(3) * t102;
t18 = qJDD(3) * pkin(7) + qJDD(2) * t99 + t102 * t71 + t206;
t24 = qJD(1) * t65 + t52 * qJDD(1);
t150 = t101 * t18 + t36 * t155 + t98 * t24;
t149 = -t73 * t156 + t207 * t99;
t148 = t205 * t166;
t147 = t205 * t163;
t146 = t98 * t162;
t145 = t61 * t156;
t144 = pkin(4) + t186;
t140 = t61 * t165 + t171;
t139 = -t101 * t24 + t35 * t155 + t36 * t163 + t98 * t18;
t138 = -t22 + t164;
t137 = t205 * t146;
t136 = t155 * t187;
t38 = t101 * t87 + t86 * t170;
t40 = t87 * t170 - t169;
t134 = -g(1) * t38 + g(2) * t40;
t39 = t86 * t161 - t87 * t98;
t41 = t87 * t161 + t86 * t98;
t133 = g(1) * t39 - g(2) * t41;
t131 = g(1) * t86 - g(2) * t87;
t6 = -pkin(4) * t205 + t160;
t129 = t101 * t7 + t6 * t98;
t128 = t101 * t6 - t7 * t98;
t100 = sin(qJ(1));
t103 = cos(qJ(1));
t127 = g(1) * t100 - g(2) * t103;
t125 = pkin(4) * t101 + qJ(5) * t98;
t123 = pkin(3) + t125;
t122 = t177 * t205 + t192;
t120 = t126 + t82;
t119 = -t155 * t205 - t98 * t55;
t118 = -t35 * t163 + t150;
t19 = -qJDD(3) * pkin(3) - qJDD(2) * t102 - t149;
t3 = pkin(4) * t23 + qJ(5) * t22 - qJD(5) * t61 + t19;
t117 = -t122 - t3;
t116 = -qJD(1) * t74 + t132;
t115 = t205 * t34 - t197;
t114 = 0.2e1 * t74 * qJD(3) - qJDD(3) * t82;
t113 = -t132 * t102 - t194;
t112 = g(1) * t40 + g(2) * t38 + t98 * t194 - t139;
t104 = qJD(3) ^ 2;
t111 = -0.2e1 * qJDD(1) * t83 - t104 * t82 + t131;
t1 = qJD(5) * t205 + t118 + t176;
t2 = qJDD(5) + t139 - t198;
t110 = t128 * qJD(4) + t1 * t101 + t2 * t98;
t109 = t10 * t61 + qJDD(5) - t112;
t108 = -g(1) * t41 - g(2) * t39 - g(3) * t172 + t118;
t107 = t59 * t165 + t119 * t99 + (-t23 - t148) * t102;
t105 = qJD(1) ^ 2;
t70 = qJDD(3) * t102 - t104 * t99;
t69 = qJDD(3) * t99 + t102 * t104;
t32 = t120 * t99;
t28 = pkin(4) * t61 + qJ(5) * t59;
t27 = -t101 * t52 + t144 * t102;
t26 = -qJ(5) * t102 + t180;
t16 = -pkin(4) * t168 - t101 * t64 + t44 * t98;
t15 = qJ(5) * t168 + t184;
t9 = (t125 * qJD(4) - qJD(5) * t101) * t99 + t120 * t156;
t8 = -t22 + t190;
t5 = t180 * qJD(4) - t144 * t165 - t174;
t4 = (-t82 * t163 - qJD(5)) * t102 + (qJ(5) - t173) * t165 + t182;
t13 = [qJDD(1), t127, g(1) * t103 + g(2) * t100, (t127 + (t96 ^ 2 + t97 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), qJDD(1) * t94 + 0.2e1 * t99 * t141, -0.2e1 * t178 * t152 + 0.2e1 * t99 * t90, t69, t70, 0, t111 * t102 + t114 * t99, t114 * t102 - t111 * t99, -t61 * t146 + (-t99 * t22 + t145) * t101, -t136 + (-t145 + (t22 + t164) * t99) * t98 + t185, -t137 + t140 + t183, (t23 - t148) * t102 + (t119 - t167) * t99, -t102 * t55 + t165 * t205, (-t52 * t163 + t174) * t205 + t52 * t175 + (t34 * t166 + (t119 + t167) * t82 + t139) * t102 + (t34 * t155 + t19 * t98 + t82 * t23 + (t186 * t205 + t11) * qJD(3)) * t99 + t133, -t182 * t205 - t180 * t55 + ((t205 * t82 - t35) * t163 + (t34 * t101 + t61 * t82) * qJD(3) + t150) * t102 + (-t34 * t163 + t19 * t101 - t82 * t22 + (t173 * t205 - t12) * qJD(3)) * t99 + t134, t23 * t32 - t27 * t55 - t5 * t205 + t59 * t9 + (t10 * t166 + t2) * t102 + (-qJD(3) * t6 + t10 * t155 + t3 * t98) * t99 + t133, -t22 * t27 - t23 * t26 - t4 * t59 + t5 * t61 + t128 * t156 + (-t129 * qJD(4) - t1 * t98 + t101 * t2 + t131) * t99, t22 * t32 + t26 * t55 + t4 * t205 - t61 * t9 + (-t10 * t154 - t1) * t102 + (qJD(3) * t7 + t10 * t163 - t101 * t3) * t99 - t134, t1 * t26 + t7 * t4 + t3 * t32 + t10 * t9 + t2 * t27 + t6 * t5 - g(1) * (-pkin(1) * t100 - pkin(4) * t39 - qJ(5) * t38) - g(2) * (pkin(1) * t103 + pkin(4) * t41 + qJ(5) * t40) + (-g(1) * pkin(6) - g(2) * t124) * t87 + (-g(2) * pkin(6) + g(1) * t124) * t86; 0, 0, 0, t159, 0, 0, 0, 0, 0, t70, -t69, 0, 0, 0, 0, 0, t107, t137 + (-t156 * t205 - t55 * t99) * t101 + t140, t107, t136 + (t138 * t99 + t145) * t98 + t185, -t171 + (-qJD(3) * t61 - t147) * t99 + t183, -g(3) + (qJD(3) * t129 - t3) * t102 + (qJD(3) * t10 + t110) * t99; 0, 0, 0, 0, -t99 * t105 * t102, t178 * t105, t153, t90, qJDD(3), qJD(3) * t45 + t159 * t102 + t116 * t99 + t149, t206 + (qJD(3) * t73 - t159) * t99 + (t116 + t207) * t102, t101 * t188 - t22 * t98, (-t188 - t23) * t98 + (-t22 - t190) * t101, (-t161 * t205 - t187) * qJD(1) - t119, -t147 + t175 + (t170 * t205 + t59 * t99) * qJD(1), -t205 * t168, -t11 * t168 - pkin(3) * t23 - t45 * t59 + (t205 * t44 + t115) * t98 + (-t192 - t19 - (t64 + t177) * t205) * t101 + t179, pkin(3) * t22 + t184 * t205 + t12 * t168 - t45 * t61 + t115 * t101 + (t122 + t19 - t203) * t98, t117 * t101 - t123 * t23 + t16 * t205 + t6 * t168 + t181 * t59 - t202 * t98 + t179, t15 * t59 - t16 * t61 + (t138 * pkin(7) - t193 + t2) * t98 + (t1 + t205 * t6 + (qJD(4) * t61 - t23) * pkin(7)) * t101 + t113, -t7 * t168 - t15 * t205 - t22 * t123 - t181 * t61 + t202 * t101 + (t117 + t203) * t98, -t7 * t15 - t6 * t16 + t181 * t10 + (t110 + t113) * pkin(7) + (-t3 - t192 + t203) * t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, -t59 ^ 2 + t200, t8, -t23 + t188, t55, -t34 * t61 + t112 + t191, t11 * t205 + t34 * t59 - t108, -t28 * t59 - t109 + t191 + 0.2e1 * t198, pkin(4) * t22 - qJ(5) * t23 + (-t12 + t7) * t61 + (t6 - t160) * t59, 0.2e1 * t176 - t10 * t59 + t28 * t61 - (-0.2e1 * qJD(5) + t11) * t205 + t108, t1 * qJ(5) - t2 * pkin(4) - t10 * t28 - t6 * t12 - g(1) * (-pkin(4) * t40 + qJ(5) * t41) - g(2) * (-pkin(4) * t38 + qJ(5) * t39) + t126 * t194 + t160 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55 + t189, t8, -t205 ^ 2 - t200, t109 - t193 - t198;];
tau_reg = t13;
