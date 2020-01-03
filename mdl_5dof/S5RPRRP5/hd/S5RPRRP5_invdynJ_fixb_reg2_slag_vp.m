% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRRP5
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
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP5_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:41:02
% EndTime: 2019-12-31 18:41:04
% DurationCPUTime: 1.18s
% Computational Cost: add. (1765->227), mult. (3124->248), div. (0->0), fcn. (1737->12), ass. (0->143)
t107 = sin(qJ(4));
t102 = t107 ^ 2;
t110 = cos(qJ(4));
t103 = t110 ^ 2;
t207 = t102 + t103;
t133 = t110 * pkin(4) + t107 * qJ(5);
t99 = qJDD(1) + qJDD(3);
t180 = t103 * t99;
t181 = t102 * t99;
t206 = t180 + t181;
t108 = sin(qJ(3));
t111 = cos(qJ(3));
t105 = sin(pkin(8));
t169 = pkin(1) * qJDD(1);
t150 = t105 * t169;
t106 = cos(pkin(8));
t87 = t106 * pkin(1) + pkin(2);
t67 = t87 * qJD(1);
t202 = qJD(3) * t67 + t150;
t192 = pkin(1) * t105;
t151 = qJD(1) * t192;
t204 = qJD(3) * t151 - t87 * qJDD(1);
t144 = -t204 * t108 + t202 * t111;
t205 = t99 * pkin(7) + qJD(2) * qJD(4) + t144;
t100 = qJD(1) + qJD(3);
t38 = t108 * t67 + t111 * t151;
t29 = t100 * pkin(7) + t38;
t179 = t107 * t29;
t20 = t110 * qJD(2) - t179;
t203 = qJD(5) - t20;
t153 = t107 * qJDD(2) + t205 * t110;
t154 = qJDD(4) * qJ(5);
t3 = t154 + (qJD(5) - t179) * qJD(4) + t153;
t158 = qJD(4) * t110;
t148 = -t110 * qJDD(2) + t205 * t107 + t29 * t158;
t166 = qJDD(4) * pkin(4);
t199 = qJDD(5) - t166;
t4 = t148 + t199;
t201 = t4 * t107 + t3 * t110;
t21 = t107 * qJD(2) + t110 * t29;
t159 = qJD(4) * t107;
t6 = -t29 * t159 + t153;
t116 = t6 * t110 + (-t107 * t21 - t110 * t20) * qJD(4) + t148 * t107;
t101 = qJ(1) + pkin(8);
t95 = qJ(3) + t101;
t86 = cos(t95);
t200 = t133 * t86;
t85 = sin(t95);
t184 = g(1) * t86 + g(2) * t85;
t17 = -qJD(4) * pkin(4) + t203;
t18 = qJD(4) * qJ(5) + t21;
t45 = -t108 * t192 + t111 * t87;
t37 = -t108 * t151 + t111 * t67;
t173 = t37 * t100;
t198 = t206 * pkin(7) - t207 * t173;
t39 = t45 * qJD(3);
t171 = t39 * t100;
t46 = t108 * t87 + t111 * t192;
t44 = pkin(7) + t46;
t197 = t207 * t171 + t206 * t44;
t160 = qJD(4) * t100;
t162 = -t102 + t103;
t163 = t107 * t110;
t196 = 0.2e1 * t162 * t160 + 0.2e1 * t99 * t163;
t79 = g(1) * t85;
t195 = g(2) * t86;
t194 = t85 * pkin(3);
t193 = t99 * pkin(3);
t113 = qJD(4) ^ 2;
t191 = pkin(7) * t113;
t190 = t100 * pkin(3);
t55 = -pkin(3) - t133;
t188 = t55 * t99;
t156 = t107 * qJD(5);
t47 = pkin(4) * t159 - qJ(5) * t158 - t156;
t187 = t47 - t38;
t177 = t107 * t86;
t178 = t107 * t85;
t186 = g(1) * t178 - g(2) * t177;
t185 = t86 * pkin(3) + t85 * pkin(7);
t112 = cos(qJ(1));
t92 = cos(t101);
t183 = t112 * pkin(1) + pkin(2) * t92;
t182 = t100 * t55;
t176 = t110 * t99;
t174 = t113 * t44;
t172 = t38 * t100;
t40 = t46 * qJD(3);
t170 = t40 * t100;
t168 = pkin(7) * qJDD(4);
t165 = t100 * t107;
t157 = qJDD(4) * t44;
t62 = t110 * t79;
t152 = t110 * t172 + t37 * t159 + t62;
t143 = t202 * t108 + t204 * t111;
t13 = t143 - t193;
t149 = -t13 - t195;
t30 = -t45 + t55;
t147 = t100 * t30 - t39;
t146 = t20 + t179;
t28 = -t37 - t190;
t145 = t13 * t107 + t28 * t158 - t186;
t142 = t183 + t185;
t139 = t158 * t165;
t109 = sin(qJ(1));
t91 = sin(t101);
t136 = -t109 * pkin(1) - pkin(2) * t91;
t135 = t191 - t193;
t134 = g(1) * t109 - g(2) * t112;
t132 = pkin(4) * t107 - qJ(5) * t110;
t131 = t17 * t107 + t18 * t110;
t130 = t20 * t107 - t21 * t110;
t129 = -t144 + t184;
t76 = t86 * pkin(7);
t128 = t136 + t76;
t8 = t188 + (t132 * qJD(4) - t156) * t100 + t143;
t127 = -t188 - t8 - t191;
t126 = -t143 + t79 - t195;
t43 = -pkin(3) - t45;
t125 = t43 * t99 + t170 + t174;
t124 = g(1) * t177 + g(2) * t178 - g(3) * t110 - t148;
t123 = t17 * t158 - t18 * t159 - t184 + t201;
t122 = t55 * t79;
t19 = t40 + t47;
t121 = -t100 * t19 - t30 * t99 - t174 - t8;
t120 = -t157 + (t100 * t43 - t39) * qJD(4);
t118 = t21 * qJD(4) + t124;
t117 = (-t107 * t18 + t110 * t17) * qJD(4) + t201;
t115 = -t184 + t116;
t104 = qJDD(2) - g(3);
t98 = t100 ^ 2;
t83 = t107 * t99;
t68 = t98 * t163;
t57 = qJDD(4) * t110 - t113 * t107;
t56 = qJDD(4) * t107 + t113 * t110;
t49 = t162 * t98;
t48 = t132 * t100;
t42 = -0.2e1 * t139 + t180;
t41 = 0.2e1 * t139 + t181;
t22 = t28 * t159;
t16 = -t37 + t182;
t14 = t16 * t159;
t1 = [0, 0, 0, 0, 0, qJDD(1), t134, g(1) * t112 + g(2) * t109, 0, 0, 0, 0, 0, 0, 0, qJDD(1), g(1) * t91 - g(2) * t92 + 0.2e1 * t106 * t169, g(1) * t92 + g(2) * t91 - 0.2e1 * t150, 0, (t134 + (t105 ^ 2 + t106 ^ 2) * t169) * pkin(1), 0, 0, 0, 0, 0, t99, t45 * t99 + t126 - t170, -t46 * t99 + t129 - t171, 0, -g(1) * t136 - g(2) * t183 - t143 * t45 + t144 * t46 - t37 * t40 + t38 * t39, t41, t196, t56, t42, t57, 0, t22 + t62 + t120 * t107 + (-t125 + t149) * t110, t107 * t125 + t110 * t120 + t145, t115 + t197, t13 * t43 + t28 * t40 - g(1) * (t128 - t194) - g(2) * t142 - t130 * t39 + t116 * t44, t41, t56, -t196, 0, -t57, t42, t14 + t62 + (qJD(4) * t147 - t157) * t107 + (t121 - t195) * t110, t123 + t197, (t157 + (-t147 - t16) * qJD(4)) * t110 + t121 * t107 + t186, t8 * t30 + t16 * t19 - g(1) * t128 - g(2) * (t142 + t200) - t122 + t131 * t39 + t117 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, 0, 0, 0, 0, 0, 0, t57, -t56, 0, -qJD(4) * t130 + t6 * t107 - t110 * t148 - g(3), 0, 0, 0, 0, 0, 0, t57, 0, t56, qJD(4) * t131 + t3 * t107 - t4 * t110 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, t126 + t172, t129 + t173, 0, 0, t41, t196, t56, t42, t57, 0, t22 + (-pkin(3) * t160 - t168) * t107 + (-t135 + t149) * t110 + t152, (-t168 + (t37 - t190) * qJD(4)) * t110 + (t135 - t172) * t107 + t145, t115 + t198, -t13 * pkin(3) - t28 * t38 - g(1) * (t76 - t194) - g(2) * t185 + t130 * t37 + t116 * pkin(7), t41, t56, -t196, 0, -t57, t42, t14 + (t160 * t55 - t168) * t107 + (-t100 * t47 + t127 - t195) * t110 + t152, t123 + t198, (t168 + (-t16 - t37 - t182) * qJD(4)) * t110 + (-t100 * t187 + t127) * t107 + t186, t8 * t55 - g(1) * t76 - g(2) * (t185 + t200) - t122 - t131 * t37 + t187 * t16 + t117 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, -t49, t83, t68, t176, qJDD(4), -t165 * t28 + t118, g(3) * t107 + t146 * qJD(4) + (-t100 * t28 + t184) * t110 - t153, 0, 0, -t68, t83, t49, qJDD(4), -t176, t68, 0.2e1 * t166 - qJDD(5) + (-t107 * t16 + t110 * t48) * t100 + t118, -t132 * t99, 0.2e1 * t154 + (t100 * t48 - g(3)) * t107 + (t100 * t16 - t184) * t110 + (0.2e1 * qJD(5) - t146) * qJD(4) + t153, -t4 * pkin(4) - g(3) * t133 + t3 * qJ(5) + t184 * t132 - t16 * t48 - t17 * t21 + t203 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) - t68, t83, -t102 * t98 - t113, -t18 * qJD(4) + t16 * t165 - t124 + t199;];
tau_reg = t1;
