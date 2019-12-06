% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
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
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:01:51
% EndTime: 2019-12-05 18:01:54
% DurationCPUTime: 1.35s
% Computational Cost: add. (1782->257), mult. (3137->294), div. (0->0), fcn. (1761->12), ass. (0->150)
t106 = sin(qJ(3));
t109 = cos(qJ(3));
t102 = sin(pkin(8));
t159 = pkin(1) * qJDD(1);
t148 = t102 * t159;
t103 = cos(pkin(8));
t83 = pkin(1) * t103 + pkin(2);
t62 = t83 * qJD(1);
t187 = qJD(3) * t62 + t148;
t181 = pkin(1) * t102;
t149 = qJD(1) * t181;
t188 = qJD(3) * t149 - t83 * qJDD(1);
t141 = -t106 * t188 + t109 * t187;
t96 = qJDD(1) + qJDD(3);
t183 = t96 * pkin(7);
t13 = t141 + t183;
t190 = qJD(2) * qJD(4) + t13;
t98 = qJ(1) + pkin(8);
t90 = qJ(3) + t98;
t81 = sin(t90);
t73 = g(3) * t81;
t82 = cos(t90);
t134 = g(2) * t82 + t73;
t140 = t106 * t187 + t109 * t188;
t184 = t96 * pkin(3);
t14 = t140 - t184;
t189 = t134 - t14;
t105 = sin(qJ(4));
t108 = cos(qJ(4));
t37 = t106 * t62 + t109 * t149;
t97 = qJD(1) + qJD(3);
t29 = pkin(7) * t97 + t37;
t92 = t108 * qJD(2);
t21 = -t105 * t29 + t92;
t152 = t105 * qJD(2);
t167 = t108 * t29;
t22 = t152 + t167;
t154 = qJD(4) * t105;
t5 = t105 * qJDD(2) + t108 * t190 - t29 * t154;
t89 = t108 * qJDD(2);
t6 = -t22 * qJD(4) - t105 * t13 + t89;
t114 = -t6 * t105 + t5 * t108 + (-t105 * t22 - t108 * t21) * qJD(4);
t74 = g(3) * t82;
t75 = g(2) * t81;
t171 = t75 - t74;
t100 = t108 ^ 2;
t99 = t105 ^ 2;
t160 = t100 + t99;
t44 = -t106 * t181 + t109 * t83;
t144 = qJ(5) * t97 + t29;
t125 = t144 * t108;
t147 = t97 * t154;
t186 = pkin(4) * t147 + qJDD(5);
t185 = pkin(4) * t99;
t153 = qJD(4) * t108;
t36 = -t106 * t149 + t109 * t62;
t179 = t108 * pkin(4);
t85 = pkin(3) + t179;
t20 = -t85 * t97 + qJD(5) - t36;
t175 = t85 * t96;
t8 = t140 - t175 + t186;
t182 = t8 * t105 + t153 * t20;
t180 = g(1) * t108;
t178 = t37 * t97;
t38 = t44 * qJD(3);
t177 = t38 * t97;
t45 = t106 * t83 + t109 * t181;
t39 = t45 * qJD(3);
t176 = t39 * t97;
t104 = -qJ(5) - pkin(7);
t28 = -pkin(3) * t97 - t36;
t174 = t14 * t105 + t153 * t28;
t16 = -t144 * t105 + t92;
t169 = qJD(4) * pkin(4);
t15 = t16 + t169;
t173 = t15 - t16;
t166 = t108 * t82;
t172 = g(2) * t166 + t108 * t73;
t170 = qJ(5) * t96;
t79 = t105 * t96;
t168 = t105 * t97;
t80 = t108 * t96;
t165 = t108 * t97;
t163 = t15 * t105;
t43 = pkin(7) + t45;
t162 = -qJ(5) - t43;
t161 = t100 - t99;
t157 = qJD(4) * t97;
t156 = qJDD(4) * pkin(4);
t155 = t105 * t108;
t91 = t108 * qJD(5);
t3 = t97 * t91 + (-t147 + t80) * qJ(5) + t5;
t150 = t108 * t3 + t171;
t146 = t104 * t81 - t82 * t85;
t145 = t160 * t36;
t143 = qJD(4) * t104;
t142 = t154 * t36 + t165 * t37 + t172;
t139 = qJD(4) * t162;
t136 = t108 * t147;
t42 = -pkin(3) - t44;
t135 = t105 * t74 - t180 + t89;
t107 = sin(qJ(1));
t86 = sin(t98);
t133 = -pkin(1) * t107 - pkin(2) * t86;
t110 = cos(qJ(1));
t87 = cos(t98);
t132 = -pkin(1) * t110 - pkin(2) * t87;
t111 = qJD(4) ^ 2;
t131 = -pkin(7) * t111 + t184;
t130 = -t28 * t97 - t75;
t129 = g(2) * t110 + g(3) * t107;
t34 = pkin(4) * t154 + t39;
t35 = t42 - t179;
t128 = -t34 * t97 - t35 * t96;
t126 = -t104 * t82 - t81 * t85;
t124 = t21 * t105 - t108 * t22;
t123 = g(1) * t105 + g(3) * t166 - t5;
t122 = -t140 + t134;
t121 = -t141 - t171;
t120 = -pkin(3) * t157 - pkin(7) * qJDD(4);
t119 = -t134 - t178;
t118 = -t111 * t43 - t42 * t96 - t176;
t117 = -qJDD(4) * t43 + (t42 * t97 - t38) * qJD(4);
t116 = -t75 - t170 + (-qJD(5) - t20) * t97;
t113 = t171 + t114;
t101 = qJDD(2) - g(1);
t95 = t97 ^ 2;
t93 = t108 * qJ(5);
t72 = t82 * pkin(7);
t66 = t95 * t155;
t65 = pkin(7) * t108 + t93;
t64 = t104 * t105;
t54 = qJDD(4) * t108 - t111 * t105;
t53 = qJDD(4) * t105 + t108 * t111;
t48 = t161 * t95;
t47 = -t105 * qJD(5) + t108 * t143;
t46 = t105 * t143 + t91;
t41 = t100 * t96 - 0.2e1 * t136;
t40 = t96 * t99 + 0.2e1 * t136;
t33 = t108 * t43 + t93;
t32 = t162 * t105;
t31 = t36 * t153;
t26 = 0.2e1 * t155 * t96 + 0.2e1 * t157 * t161;
t23 = t28 * t154;
t18 = t20 * t154;
t17 = t152 + t125;
t10 = (-qJD(5) - t38) * t105 + t108 * t139;
t9 = t105 * t139 + t108 * t38 + t91;
t2 = t156 + t89 - qJD(4) * t125 + (-qJD(5) * t97 - t170 - t190) * t105;
t1 = [0, 0, 0, 0, 0, qJDD(1), t129, -g(2) * t107 + g(3) * t110, 0, 0, 0, 0, 0, 0, 0, qJDD(1), g(2) * t87 + g(3) * t86 + 0.2e1 * t103 * t159, -g(2) * t86 + g(3) * t87 - 0.2e1 * t148, 0, (t129 + (t102 ^ 2 + t103 ^ 2) * t159) * pkin(1), 0, 0, 0, 0, 0, t96, t44 * t96 + t122 - t176, -t45 * t96 + t121 - t177, 0, -g(2) * t132 - g(3) * t133 - t140 * t44 + t141 * t45 - t36 * t39 + t37 * t38, t40, t26, t53, t41, t54, 0, t23 + t117 * t105 + (t118 - t14) * t108 + t172, t117 * t108 + (-t118 - t134) * t105 + t174, t113 + t160 * (t43 * t96 + t177), t14 * t42 + t28 * t39 - g(2) * (-pkin(3) * t82 - pkin(7) * t81 + t132) - g(3) * (-pkin(3) * t81 + t133 + t72) - t124 * t38 + t114 * t43, t40, t26, t53, t41, t54, 0, t32 * qJDD(4) + t18 + (t168 * t35 + t10) * qJD(4) + (t128 - t8) * t108 + t172, -t33 * qJDD(4) + (t165 * t35 - t9) * qJD(4) + (-t128 - t134) * t105 + t182, (t33 * t96 + t9 * t97 + (-t32 * t97 - t15) * qJD(4)) * t108 + (-t10 * t97 - t32 * t96 - t2 + (-t33 * t97 - t17) * qJD(4)) * t105 + t150, t3 * t33 + t17 * t9 + t2 * t32 + t15 * t10 + t8 * t35 + t20 * t34 - g(2) * (t132 + t146) - g(3) * (t126 + t133); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, 0, 0, 0, 0, 0, 0, t54, -t53, 0, -qJD(4) * t124 + t5 * t105 + t6 * t108 - g(1), 0, 0, 0, 0, 0, 0, t54, -t53, 0, t3 * t105 + t2 * t108 - g(1) + (t108 * t17 - t163) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, t122 + t178, t36 * t97 + t121, 0, 0, t40, t26, t53, t41, t54, 0, t23 + t120 * t105 + (t131 - t14) * t108 + t142, t31 + t120 * t108 + (t119 - t131) * t105 + t174, -t145 * t97 + t160 * t183 + t113, -g(3) * t72 - t28 * t37 + t124 * t36 + t189 * pkin(3) + (t114 + t75) * pkin(7), t40, t26, t53, t41, t54, 0, t64 * qJDD(4) + t18 + (-t8 + t175) * t108 + (t47 + (-t85 - t179) * t168) * qJD(4) + t142, -t65 * qJDD(4) + t31 + (-t46 + (-t108 * t85 + t185) * t97) * qJD(4) + (t119 - t175) * t105 + t182, (-qJD(4) * t15 + t65 * t96) * t108 + (-t17 * qJD(4) - t64 * t96 - t2) * t105 + (-t105 * t47 + t108 * t46 - t145 + (-t105 * t65 - t108 * t64) * qJD(4)) * t97 + t150, t3 * t65 + t2 * t64 + t15 * t47 - t8 * t85 - t20 * t37 - g(2) * t146 - g(3) * t126 + (-t108 * t36 + t46) * t17 + (t15 * t36 + t169 * t20) * t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, -t48, t79, t66, t80, qJDD(4), (t22 - t167) * qJD(4) + (t130 - t190) * t105 + t135, t21 * qJD(4) + t108 * t130 + t123, 0, 0, -t66, -t48, t79, t66, t80, qJDD(4), 0.2e1 * t156 + (t17 - t125) * qJD(4) + (t179 * t95 + t116 - t190) * t105 + t135, -t95 * t185 + (qJ(5) * t168 + t16) * qJD(4) + t116 * t108 + t123, -pkin(4) * t79 + (-t169 + t173) * t165, t173 * t17 + (-t180 + t2 + (-t20 * t97 - t171) * t105) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80 + 0.2e1 * t147, 0.2e1 * t153 * t97 + t79, -t160 * t95, t97 * t163 + (-pkin(4) * t96 - t17 * t97) * t108 + t186 - t189;];
tau_reg = t1;
