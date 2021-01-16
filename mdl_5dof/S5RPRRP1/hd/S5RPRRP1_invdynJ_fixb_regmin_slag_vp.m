% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRP1
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% tau_reg [5x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:27
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:27:19
% EndTime: 2021-01-15 12:27:26
% DurationCPUTime: 1.38s
% Computational Cost: add. (1748->227), mult. (3323->271), div. (0->0), fcn. (2111->8), ass. (0->142)
t101 = sin(qJ(1));
t104 = cos(qJ(1));
t183 = g(1) * t101 - g(2) * t104;
t98 = qJ(3) + qJ(4);
t86 = sin(t98);
t87 = cos(t98);
t190 = g(3) * t86 - t183 * t87;
t189 = qJDD(2) - t183;
t103 = cos(qJ(3));
t149 = qJD(1) * qJD(3);
t136 = t103 * t149;
t100 = sin(qJ(3));
t147 = t100 * qJDD(1);
t188 = t136 + t147;
t146 = t103 * qJDD(1);
t187 = t100 * t149 - t146;
t131 = g(1) * t104 + g(2) * t101;
t94 = (qJD(1) * qJD(2));
t119 = -t131 + (2 * t94);
t93 = qJDD(1) * qJ(2);
t186 = 0.2e1 * t93 + t119;
t158 = qJDD(1) * pkin(1);
t185 = t158 - t189;
t102 = cos(qJ(4));
t99 = sin(qJ(4));
t55 = t102 * t100 + t99 * t103;
t45 = t55 * qJD(1);
t92 = qJD(3) + qJD(4);
t172 = t45 * t92;
t105 = -pkin(1) - pkin(6);
t169 = pkin(7) - t105;
t59 = t169 * t100;
t60 = t169 * t103;
t165 = -t102 * t59 - t99 * t60;
t150 = qJD(4) * t102;
t175 = pkin(3) * t92;
t91 = qJDD(3) + qJDD(4);
t184 = -t99 * pkin(3) * t91 - t150 * t175;
t154 = qJD(1) * t100;
t68 = t105 * qJD(1) + qJD(2);
t37 = -pkin(7) * t154 + t100 * t68;
t31 = t99 * t37;
t153 = qJD(1) * t103;
t38 = -pkin(7) * t153 + t103 * t68;
t34 = qJD(3) * pkin(3) + t38;
t140 = t102 * t34 - t31;
t143 = t99 * t154;
t47 = t102 * t153 - t143;
t39 = t47 * qJ(5);
t11 = t140 - t39;
t129 = -t102 * t146 + t99 * t147;
t25 = t92 * t55;
t15 = t25 * qJD(1) + t129;
t162 = t15 * qJ(5);
t84 = t91 * pkin(4);
t182 = t162 + t84;
t152 = qJD(3) * t100;
t67 = t105 * qJDD(1) + qJDD(2);
t57 = t103 * t67;
t21 = qJDD(3) * pkin(3) + t187 * pkin(7) - t68 * t152 + t57;
t151 = qJD(3) * t103;
t23 = -t188 * pkin(7) + t100 * t67 + t68 * t151;
t181 = t102 * t21 - t99 * t23;
t159 = qJD(4) * t99;
t179 = (qJD(4) * t34 + t23) * t102 - t37 * t159 + t99 * t21;
t178 = t47 ^ 2;
t174 = g(3) * t100;
t173 = t45 * t47;
t171 = t47 * t92;
t10 = t92 * pkin(4) + t11;
t168 = t10 - t11;
t167 = t102 * t38 - t31;
t56 = -t99 * t100 + t102 * t103;
t166 = -t92 * t25 + t91 * t56;
t61 = pkin(3) * t154 + qJD(1) * qJ(2);
t97 = t103 ^ 2;
t163 = t100 ^ 2 - t97;
t32 = t102 * t37;
t125 = -qJD(4) * t143 - t187 * t99;
t133 = t92 * t103;
t16 = (qJD(1) * t133 + t147) * t102 + t125;
t161 = t16 * qJ(5);
t160 = t45 * qJ(5);
t107 = qJD(1) ^ 2;
t157 = t107 * qJ(2);
t135 = -t45 * pkin(4) - qJD(5);
t27 = -t135 + t61;
t156 = qJD(5) + t27;
t106 = qJD(3) ^ 2;
t155 = -t106 - t107;
t70 = pkin(3) * t151 + qJD(2);
t148 = qJDD(3) * t100;
t144 = pkin(3) * t153;
t76 = t100 * pkin(3) + qJ(2);
t142 = pkin(4) * t86 + t76;
t139 = -t99 * t38 - t32;
t138 = -t102 * t60 + t99 * t59;
t35 = t188 * pkin(3) + t93 + t94;
t128 = -t56 * t15 - t47 * t25;
t26 = -t100 * t159 + t102 * t133 - t99 * t152;
t127 = -t26 * t92 - t55 * t91;
t126 = -t99 * t34 - t32;
t124 = t131 * t86;
t123 = t131 * t87;
t53 = t169 * t152;
t54 = qJD(3) * t60;
t122 = -t102 * t54 - t60 * t150 + t59 * t159 + t99 * t53;
t121 = 0.2e1 * qJ(2) * t149 + qJDD(3) * t105;
t5 = t16 * pkin(4) + qJDD(5) + t35;
t120 = -t183 - t157;
t117 = t126 * qJD(4) + t181;
t116 = -t165 * qJD(4) + t102 * t53 + t99 * t54;
t1 = -t47 * qJD(5) + t117 + t182;
t12 = -t126 - t160;
t2 = -t45 * qJD(5) - t161 + t179;
t115 = t1 * t56 - t10 * t25 + t12 * t26 + t2 * t55 - t183;
t114 = -t105 * t106 + t186;
t113 = g(3) * t87 + t183 * t86 - t179;
t112 = t117 + t190;
t111 = t61 * t45 + t113;
t110 = -t61 * t47 + t112;
t109 = t156 * t45 + t113 + t161;
t108 = -t129 - t172;
t88 = qJ(5) + t169;
t85 = qJDD(3) * t103;
t79 = t102 * pkin(3) + pkin(4);
t44 = t45 ^ 2;
t36 = t55 * pkin(4) + t76;
t29 = t47 * pkin(4) + t144;
t22 = t26 * pkin(4) + t70;
t20 = -t55 * qJ(5) + t165;
t19 = -t56 * qJ(5) + t138;
t17 = -t44 + t178;
t14 = -t39 + t167;
t13 = t139 + t160;
t9 = -qJD(1) * t45 + t166;
t8 = -qJD(1) * t47 + t127;
t7 = t171 + (-t92 * t153 - t147) * t102 - t125;
t6 = t108 + t172;
t4 = t25 * qJ(5) - t56 * qJD(5) + t116;
t3 = -t26 * qJ(5) - t55 * qJD(5) + t122;
t18 = [qJDD(1), t183, t131, -0.2e1 * t158 + t189, t186, t185 * pkin(1) + (t119 + t93) * qJ(2), t97 * qJDD(1) - 0.2e1 * t100 * t136, -0.2e1 * t100 * t146 + 0.2e1 * t163 * t149, -t106 * t100 + t85, -t106 * t103 - t148, 0, t114 * t100 + t121 * t103, -t121 * t100 + t114 * t103, t128, t55 * t15 - t16 * t56 + t45 * t25 - t26 * t47, t166, t127, 0, t116 * t92 + t138 * t91 + t76 * t16 + t61 * t26 + t35 * t55 + t70 * t45 - t124, -t122 * t92 - t76 * t15 - t165 * t91 - t61 * t25 + t35 * t56 + t70 * t47 - t123, t36 * t16 + t19 * t91 + t22 * t45 + t27 * t26 + t4 * t92 + t5 * t55 - t124, -t36 * t15 - t20 * t91 + t22 * t47 - t27 * t25 - t3 * t92 + t5 * t56 - t123, t19 * t15 - t20 * t16 - t3 * t45 - t4 * t47 - t115, t2 * t20 + t12 * t3 + t1 * t19 + t10 * t4 + t5 * t36 + t27 * t22 - g(1) * (-t88 * t101 + t104 * t142) - g(2) * (t101 * t142 + t88 * t104); 0, 0, 0, qJDD(1), -t107, -t157 - t185, 0, 0, 0, 0, 0, t155 * t100 + t85, t155 * t103 - t148, 0, 0, 0, 0, 0, t9, t8, t9, t8, -t55 * t16 - t26 * t45 - t128, -t27 * qJD(1) + t115; 0, 0, 0, 0, 0, 0, t103 * t107 * t100, -t163 * t107, t146, -t147, qJDD(3), t120 * t103 + t174 + t57, g(3) * t103 + (-t120 - t67) * t100, t173, t17, t6, t7, t91, -t139 * t92 + (t102 * t91 - t153 * t45 - t159 * t92) * pkin(3) + t110, -t144 * t47 + t167 * t92 + t111 + t184, -t13 * t92 - t29 * t45 + t79 * t91 - t156 * t47 + (-t32 + (-t34 - t175) * t99) * qJD(4) + t190 + t181 + t182, t14 * t92 - t29 * t47 + t109 + t184, t79 * t15 + (t12 + t13) * t47 + (-t10 + t14) * t45 + (-t16 * t99 + (-t102 * t45 + t47 * t99) * qJD(4)) * pkin(3), t1 * t79 - t10 * t13 - t12 * t14 - t27 * t29 + t190 * pkin(4) + (t174 + t2 * t99 - t183 * t103 + (-t10 * t99 + t102 * t12) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, t17, t6, t7, t91, -t126 * t92 + t110, t140 * t92 + t111, t162 + t12 * t92 + 0.2e1 * t84 + (t135 - t27) * t47 + t112, -t178 * pkin(4) + t11 * t92 + t109, t15 * pkin(4) - t168 * t45, t168 * t12 + (-t27 * t47 + t1 + t190) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 + t171, t108 - t172, -t44 - t178, t10 * t47 + t12 * t45 - t131 + t5;];
tau_reg = t18;
