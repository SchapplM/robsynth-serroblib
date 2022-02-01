% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:20:55
% EndTime: 2022-01-23 09:20:58
% DurationCPUTime: 1.36s
% Computational Cost: add. (2398->242), mult. (4079->328), div. (0->0), fcn. (2468->14), ass. (0->158)
t106 = sin(pkin(8));
t193 = pkin(1) * t106;
t150 = qJD(3) * t193;
t108 = cos(pkin(8));
t89 = t108 * pkin(1) + pkin(2);
t197 = qJD(1) * t150 - t89 * qJDD(1);
t105 = sin(pkin(9));
t101 = qJ(1) + pkin(8);
t94 = qJ(3) + t101;
t88 = cos(t94);
t181 = t105 * t88;
t107 = cos(pkin(9));
t192 = t107 * pkin(4);
t196 = pkin(7) * t181 + t88 * t192;
t168 = pkin(1) * qJDD(1);
t100 = qJD(1) + qJD(3);
t98 = t105 ^ 2;
t99 = t107 ^ 2;
t183 = t98 + t99;
t195 = t100 * t183;
t110 = sin(qJ(3));
t113 = cos(qJ(3));
t151 = qJD(1) * t193;
t69 = t89 * qJD(1);
t41 = -t110 * t151 + t113 * t69;
t132 = qJD(4) - t41;
t87 = sin(t94);
t80 = g(1) * t87;
t79 = g(2) * t88;
t97 = qJDD(1) + qJDD(3);
t194 = t97 * pkin(3);
t96 = t100 ^ 2;
t191 = t98 * t96;
t109 = sin(qJ(5));
t112 = cos(qJ(5));
t158 = qJD(4) * t107;
t162 = t107 * t112;
t163 = t107 * t109;
t59 = -t105 * pkin(7) - pkin(3) - t192;
t43 = -qJ(4) * t163 + t112 * t59;
t166 = qJD(5) * t43;
t176 = t110 * t69;
t42 = t113 * t151 + t176;
t190 = -t109 * t42 + t112 * t158 - t41 * t162 + t166;
t44 = qJ(4) * t162 + t109 * t59;
t165 = qJD(5) * t44;
t189 = -t109 * t158 - t112 * t42 + t41 * t163 - t165;
t188 = g(2) * t181 - t105 * t80;
t187 = t88 * pkin(3) + t87 * qJ(4);
t186 = -g(1) * t88 - g(2) * t87;
t185 = t79 - t80;
t51 = t110 * t89 + t113 * t193;
t114 = cos(qJ(1));
t93 = cos(t101);
t184 = t114 * pkin(1) + pkin(2) * t93;
t25 = t59 * t100 + t132;
t33 = t100 * qJ(4) + t42;
t29 = t105 * qJD(2) + t107 * t33;
t8 = t109 * t25 + t112 * t29;
t182 = qJD(5) * t8;
t180 = t105 * t97;
t179 = t107 * t87;
t178 = t107 * t97;
t177 = t109 * t97;
t175 = t112 * t97;
t148 = t106 * t168;
t159 = qJD(3) * t113;
t140 = -t197 * t110 + t113 * t148 + t69 * t159;
t169 = t97 * qJ(4);
t18 = t100 * qJD(4) + t140 + t169;
t14 = -t107 * qJDD(2) + t105 * t18;
t174 = t14 * t105;
t28 = -t107 * qJD(2) + t105 * t33;
t173 = t28 * t105;
t172 = t42 * t100;
t47 = t51 * qJD(3);
t171 = t47 * t100;
t68 = -qJDD(5) + t178;
t170 = t68 * t107;
t50 = -t110 * t193 + t113 * t89;
t34 = -t50 + t59;
t48 = qJ(4) + t51;
t20 = t112 * t34 - t48 * t163;
t167 = qJD(5) * t20;
t164 = t107 * t100;
t161 = t109 * t112;
t102 = t109 ^ 2;
t103 = t112 ^ 2;
t160 = t102 - t103;
t157 = qJD(5) * t109;
t156 = qJD(5) * t112;
t155 = pkin(7) * t80;
t144 = t105 * t157;
t7 = -t109 * t29 + t112 * t25;
t154 = t7 * t144 - t188;
t153 = pkin(4) * t179;
t139 = qJD(3) * t176 + t110 * t148 + t197 * t113;
t127 = qJDD(4) + t139;
t24 = t127 - t194;
t152 = t24 * t105 + t188;
t71 = -qJD(5) + t164;
t149 = t71 * t157;
t147 = -t87 * pkin(3) + t88 * qJ(4);
t146 = -t24 - t79;
t145 = t183 * t97;
t143 = -t14 * t107 - g(3);
t142 = t68 - t178;
t141 = t68 + t178;
t138 = t184 + t187;
t137 = t100 * (-qJD(5) - t71);
t136 = t161 * t191;
t134 = t100 * t109 * t156;
t111 = sin(qJ(1));
t92 = sin(t101);
t133 = -t111 * pkin(1) - pkin(2) * t92;
t131 = -t172 - t194;
t130 = g(1) * t111 - g(2) * t114;
t129 = t109 * t7 - t112 * t8;
t46 = -t110 * t150 + t89 * t159;
t45 = qJD(4) + t46;
t128 = t14 * t48 + t28 * t45;
t49 = -pkin(3) - t50;
t126 = t49 * t97 + t171;
t15 = t105 * qJDD(2) + t107 * t18;
t125 = t15 * t107 + t174 + t186;
t124 = -t139 - t185;
t123 = -t140 - t186;
t122 = g(3) * t105 - qJD(5) * t25 - t15;
t21 = t109 * t34 + t48 * t162;
t121 = -qJD(5) * t29 - t100 * t173;
t120 = t14 * qJ(4) + t132 * t28;
t119 = t133 + t147;
t118 = -t71 ^ 2 - t191;
t13 = t59 * t97 + t127;
t11 = t112 * t13;
t3 = -t109 * t15 + t11 - t182;
t38 = t88 * t109 - t87 * t162;
t40 = t87 * t109 + t88 * t162;
t117 = -g(1) * t38 - g(2) * t40 - t3 * t107 + t109 * t174 + t156 * t173;
t2 = t7 * qJD(5) + t109 * t13 + t112 * t15;
t37 = t88 * t112 + t87 * t163;
t39 = t87 * t112 - t88 * t163;
t116 = -g(1) * t37 - g(2) * t39 + t2 * t107 + t112 * t174 - t28 * t144;
t104 = qJDD(2) - g(3);
t83 = t99 * t97;
t82 = t98 * t97;
t64 = g(1) * t179;
t57 = 0.2e1 * t105 * t178;
t54 = t144 * t164;
t36 = (t103 * t97 - 0.2e1 * t134) * t98;
t35 = (t102 * t97 + 0.2e1 * t134) * t98;
t32 = -t100 * pkin(3) + t132;
t27 = 0.2e1 * (t160 * t100 * qJD(5) - t97 * t161) * t98;
t17 = (t141 * t109 + (t71 + t164) * t156) * t105;
t16 = t54 + (-t141 * t112 + t149) * t105;
t5 = -t21 * qJD(5) + t112 * t47 - t45 * t163;
t4 = t109 * t47 + t45 * t162 + t167;
t1 = [0, 0, 0, 0, 0, qJDD(1), t130, g(1) * t114 + g(2) * t111, 0, 0, 0, 0, 0, 0, 0, qJDD(1), g(1) * t92 - g(2) * t93 + 0.2e1 * t108 * t168, g(1) * t93 + g(2) * t92 - 0.2e1 * t148, 0, (t130 + (t106 ^ 2 + t108 ^ 2) * t168) * pkin(1), 0, 0, 0, 0, 0, t97, t50 * t97 + t124 - t171, -t46 * t100 - t51 * t97 + t123, 0, -g(1) * t133 - g(2) * t184 - t139 * t50 + t140 * t51 - t41 * t47 + t42 * t46, t82, t57, 0, t83, 0, 0, t64 + (-t126 + t146) * t107, t126 * t105 + t152, t48 * t145 + t45 * t195 + t125, t24 * t49 + t32 * t47 - g(1) * t119 - g(2) * t138 + (t15 * t48 + t29 * t45) * t107 + t128 * t105, t36, t27, t16, t35, t17, t170, -t20 * t68 - t5 * t71 + (t48 * t177 + (t109 * t45 + t156 * t48) * t100) * t98 + t117, t21 * t68 + t4 * t71 + (t48 * t175 + (t112 * t45 - t157 * t48) * t100) * t98 + t116, ((-t21 * t97 - t2 + (-t4 + t167) * t100) * t109 + (-t100 * t5 - t20 * t97 - t3 + (-t100 * t21 - t8) * qJD(5)) * t112) * t105 + t154, t2 * t21 + t8 * t4 + t3 * t20 + t7 * t5 - g(1) * (t119 - t153) - g(2) * (t138 + t196) + (t128 + t155) * t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 * t105 + t143, 0, 0, 0, 0, 0, 0, (t142 * t109 + (t71 - t164) * t156) * t105, t54 + (t112 * t142 - t149) * t105, 0, (-t109 * t3 + t112 * t2 + (-t109 * t8 - t112 * t7) * qJD(5)) * t105 + t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, t124 + t172, t41 * t100 + t123, 0, 0, t82, t57, 0, t83, 0, 0, t64 + (-t131 + t146) * t107, t131 * t105 + t152, qJ(4) * t145 + t132 * t195 + t125, -t24 * pkin(3) - t32 * t42 - g(1) * t147 - g(2) * t187 + (t15 * qJ(4) + t132 * t29) * t107 + t120 * t105, t36, t27, t16, t35, t17, t170, -t43 * t68 - t189 * t71 + (t109 * t169 + (qJ(4) * t156 + t109 * t132) * t100) * t98 + t117, t44 * t68 + t190 * t71 + (t112 * t169 + (-qJ(4) * t157 + t112 * t132) * t100) * t98 + t116, ((-t44 * t97 - t2) * t109 + (-t43 * t97 - t182 - t3) * t112 + ((-t165 - t189) * t112 + (t166 - t190) * t109) * t100) * t105 + t154, t2 * t44 + t3 * t43 - g(1) * (t147 - t153) - g(2) * (t187 + t196) + t190 * t8 + t189 * t7 + (t120 + t155) * t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t178, t180, -t183 * t96, (-t29 * t107 - t173) * t100 + t24 + t185, 0, 0, 0, 0, 0, 0, t109 * t118 - t112 * t68, t109 * t68 + t112 * t118, (-t102 - t103) * t180, t2 * t109 + t3 * t112 - t129 * qJD(5) + (t107 * t129 - t173) * t100 + t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, -t160 * t191, (t109 * t137 + t175) * t105, -t136, (t112 * t137 - t177) * t105, -t68, -g(1) * t39 + g(2) * t37 + t109 * t122 + t112 * t121 - t8 * t71 + t11, g(1) * t40 - g(2) * t38 - t7 * t71 + t122 * t112 + (-t121 - t13) * t109, 0, 0;];
tau_reg = t1;
