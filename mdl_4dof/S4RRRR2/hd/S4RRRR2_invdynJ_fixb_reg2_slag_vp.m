% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RRRR2
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
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRRR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:21
% EndTime: 2019-12-31 17:23:24
% DurationCPUTime: 1.36s
% Computational Cost: add. (2056->236), mult. (3278->318), div. (0->0), fcn. (2004->12), ass. (0->162)
t121 = sin(qJ(3));
t116 = t121 ^ 2;
t125 = cos(qJ(3));
t117 = t125 ^ 2;
t184 = t116 + t117;
t113 = qJDD(1) + qJDD(2);
t122 = sin(qJ(2));
t176 = qJDD(1) * t122;
t126 = cos(qJ(2));
t181 = qJD(2) * t126;
t56 = t113 * pkin(6) + (qJD(1) * t181 + t176) * pkin(1);
t159 = t184 * t56;
t119 = qJ(1) + qJ(2);
t108 = cos(t119);
t212 = g(2) * t108;
t200 = pkin(1) * qJD(1);
t172 = t122 * t200;
t208 = t126 * pkin(1);
t196 = qJD(2) * t172 - qJDD(1) * t208;
t210 = t113 * pkin(2);
t55 = t196 - t210;
t219 = t55 + t212;
t115 = qJD(1) + qJD(2);
t214 = t184 * t126;
t218 = t115 * t214;
t124 = cos(qJ(4));
t177 = qJD(4) * t124;
t179 = qJD(3) * t125;
t217 = -t124 * t179 - t125 * t177;
t68 = t115 * pkin(6) + t172;
t183 = qJD(1) * t126;
t171 = pkin(1) * t183;
t209 = t115 * pkin(2);
t69 = -t171 - t209;
t216 = t122 * t69 + t214 * t68;
t106 = sin(t119);
t215 = g(1) * t108 + g(2) * t106;
t120 = sin(qJ(4));
t188 = t124 * t121;
t64 = t120 * t125 + t188;
t53 = t64 * t115;
t114 = qJD(3) + qJD(4);
t178 = qJD(4) * t120;
t166 = t115 * t179;
t189 = t121 * t113;
t18 = -t68 * t179 + qJDD(3) * pkin(3) - t121 * t56 + (-t166 - t189) * pkin(7);
t180 = qJD(3) * t121;
t186 = t125 * t113;
t139 = t115 * t180 - t186;
t23 = -t139 * pkin(7) + t125 * t56 - t68 * t180;
t163 = pkin(7) * t115 + t68;
t44 = t163 * t121;
t39 = qJD(3) * pkin(3) - t44;
t45 = t163 * t125;
t4 = (qJD(4) * t39 + t23) * t124 + t120 * t18 - t45 * t178;
t128 = -pkin(7) - pkin(6);
t98 = g(1) * t106;
t123 = sin(qJ(1));
t213 = g(1) * t123;
t211 = g(3) * t125;
t190 = t120 * t121;
t170 = t115 * t190;
t187 = t124 * t125;
t51 = -t115 * t187 + t170;
t207 = t53 * t51;
t101 = t122 * pkin(1) + pkin(6);
t206 = -pkin(7) - t101;
t86 = t128 * t121;
t109 = t125 * pkin(7);
t87 = t125 * pkin(6) + t109;
t40 = -t120 * t87 + t124 * t86;
t63 = -t187 + t190;
t168 = qJD(3) * t128;
t65 = t121 * t168;
t66 = t125 * t168;
t205 = t40 * qJD(4) + t120 * t66 + t124 * t65 + t63 * t171;
t41 = t120 * t86 + t124 * t87;
t204 = -t41 * qJD(4) - t120 * t65 + t124 * t66 + t64 * t171;
t203 = t125 * t98 + t69 * t180;
t202 = t108 * pkin(2) + t106 * pkin(6);
t199 = t120 * t45;
t197 = t124 * t45;
t118 = qJ(3) + qJ(4);
t105 = sin(t118);
t195 = t105 * t106;
t194 = t105 * t108;
t107 = cos(t118);
t193 = t106 * t107;
t192 = t107 * t108;
t191 = t115 * t121;
t185 = t116 - t117;
t182 = qJD(2) * t122;
t175 = t219 * t121 + t69 * t179;
t174 = pkin(1) * t181;
t173 = pkin(3) * t180;
t111 = t115 ^ 2;
t169 = t121 * t111 * t125;
t102 = t125 * pkin(3) + pkin(2);
t167 = t115 * t182;
t160 = qJD(3) * t206;
t158 = -t215 + t159;
t157 = -t113 * t188 + t217 * t115 - t120 * t186;
t155 = t108 * t102 - t106 * t128;
t153 = t184 * t113;
t151 = t115 * t172;
t150 = t121 * t166;
t149 = g(1) * (-t106 * pkin(2) + t108 * pkin(6));
t148 = t120 * t189 - t124 * t186;
t127 = cos(qJ(1));
t146 = -g(2) * t127 + t213;
t145 = t114 * t190;
t21 = t124 * t39 - t199;
t22 = t120 * t39 + t197;
t35 = t145 + t217;
t36 = t114 * t64;
t5 = -t22 * qJD(4) - t120 * t23 + t124 * t18;
t144 = t21 * t35 - t22 * t36 - t4 * t63 - t5 * t64 - t215;
t61 = t206 * t121;
t62 = t125 * t101 + t109;
t33 = -t120 * t62 + t124 * t61;
t34 = t120 * t61 + t124 * t62;
t143 = -t196 + t98 - t212;
t142 = -t106 * t102 - t108 * t128;
t32 = t139 * pkin(3) + t55;
t54 = -t102 * t115 - t171;
t141 = -g(1) * t195 + g(2) * t194 + t32 * t64 - t54 * t35;
t140 = g(1) * t193 - g(2) * t192 + t32 * t63 + t54 * t36;
t138 = -t172 + t173;
t137 = -t115 * t69 + t215 - t56;
t129 = qJD(3) ^ 2;
t136 = pkin(6) * t129 - t151 - t210;
t103 = -pkin(2) - t208;
t135 = pkin(1) * t167 + t101 * t129 + t103 * t113;
t134 = -pkin(6) * qJDD(3) + (t171 - t209) * qJD(3);
t133 = -qJDD(3) * t101 + (t103 * t115 - t174) * qJD(3);
t132 = g(1) * t192 + g(2) * t193 + g(3) * t105 + t54 * t51 - t4;
t131 = g(1) * t194 + g(2) * t195 - g(3) * t107 - t54 * t53 + t5;
t112 = qJDD(3) + qJDD(4);
t110 = t127 * pkin(1);
t81 = -t102 - t208;
t80 = qJDD(3) * t125 - t129 * t121;
t79 = qJDD(3) * t121 + t129 * t125;
t67 = pkin(1) * t182 + t173;
t58 = t117 * t113 - 0.2e1 * t150;
t57 = t116 * t113 + 0.2e1 * t150;
t43 = -t121 * t174 + t125 * t160;
t42 = t121 * t160 + t125 * t174;
t37 = -0.2e1 * t185 * t115 * qJD(3) + 0.2e1 * t121 * t186;
t27 = -t63 * t112 - t36 * t114;
t26 = t64 * t112 - t35 * t114;
t25 = -t124 * t44 - t199;
t24 = t120 * t44 - t197;
t17 = -t51 ^ 2 + t53 ^ 2;
t15 = t36 * t115 + t148;
t14 = t115 * t145 + t157;
t10 = -t157 + (-t170 + t51) * t114;
t9 = -t34 * qJD(4) - t120 * t42 + t124 * t43;
t8 = t33 * qJD(4) + t120 * t43 + t124 * t42;
t7 = t15 * t63 + t51 * t36;
t6 = -t14 * t64 - t53 * t35;
t1 = t14 * t63 - t64 * t15 + t35 * t51 - t53 * t36;
t2 = [0, 0, 0, 0, 0, qJDD(1), t146, g(1) * t127 + g(2) * t123, 0, 0, 0, 0, 0, 0, 0, t113, (t113 * t126 - t167) * pkin(1) + t143, ((-qJDD(1) - t113) * t122 + (-qJD(1) - t115) * t181) * pkin(1) + t215, 0, (t146 + (t122 ^ 2 + t126 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t57, t37, t79, t58, t80, 0, t133 * t121 + (-t135 - t219) * t125 + t203, t133 * t125 + (t135 - t98) * t121 + t175, pkin(1) * qJD(2) * t218 + t101 * t153 + t158, t55 * t103 - t149 - g(2) * (t110 + t202) + t101 * t159 + (t216 * qJD(2) + t213) * pkin(1), t6, t1, t26, t7, t27, 0, t33 * t112 + t9 * t114 + t81 * t15 + t67 * t51 + t140, -t34 * t112 - t8 * t114 - t81 * t14 + t67 * t53 + t141, t33 * t14 - t34 * t15 - t8 * t51 - t9 * t53 + t144, t4 * t34 + t22 * t8 + t5 * t33 + t21 * t9 + t32 * t81 + t54 * t67 - g(1) * (-t123 * pkin(1) + t142) - g(2) * (t110 + t155); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113, t143 + t151, (-t176 + (-qJD(2) + t115) * t183) * pkin(1) + t215, 0, 0, t57, t37, t79, t58, t80, 0, t134 * t121 + (-t136 - t219) * t125 + t203, t134 * t125 + (t136 - t98) * t121 + t175, pkin(6) * t153 - t200 * t218 + t158, -t55 * pkin(2) + pkin(6) * t159 - g(2) * t202 - t216 * t200 - t149, t6, t1, t26, t7, t27, 0, -t102 * t15 + t40 * t112 + t204 * t114 + t138 * t51 + t140, t102 * t14 - t41 * t112 - t205 * t114 + t138 * t53 + t141, t40 * t14 - t41 * t15 - t204 * t53 - t205 * t51 + t144, -g(1) * t142 - g(2) * t155 - t32 * t102 + t138 * t54 + t204 * t21 + t205 * t22 + t4 * t41 + t5 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t169, t185 * t111, t189, t169, t186, qJDD(3), t137 * t121 - t211, g(3) * t121 + t137 * t125, 0, 0, t207, t17, t10, -t207, -t148, t112, -t24 * t114 + (t112 * t124 - t114 * t178 - t191 * t51) * pkin(3) + t131, t25 * t114 + (-t112 * t120 - t114 * t177 - t191 * t53) * pkin(3) + t132, (t22 + t24) * t53 + (-t21 + t25) * t51 + (-t120 * t15 + t124 * t14 + (t120 * t53 - t124 * t51) * qJD(4)) * pkin(3), -t21 * t24 - t22 * t25 + (-t211 + t120 * t4 + t124 * t5 + (-t120 * t21 + t124 * t22) * qJD(4) + (-t115 * t54 + t215) * t121) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t207, t17, t10, -t207, -t148, t112, t22 * t114 + t131, t21 * t114 + t132, 0, 0;];
tau_reg = t2;
