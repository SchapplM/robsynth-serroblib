% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRR12
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% tau_reg [5x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRR12_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR12_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR12_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR12_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:13:11
% EndTime: 2019-12-31 19:13:17
% DurationCPUTime: 1.65s
% Computational Cost: add. (2112->263), mult. (4151->352), div. (0->0), fcn. (2890->10), ass. (0->160)
t117 = sin(qJ(4));
t118 = sin(qJ(3));
t121 = cos(qJ(4));
t122 = cos(qJ(3));
t64 = t117 * t122 + t121 * t118;
t59 = t64 * qJD(1);
t214 = qJD(5) + t59;
t109 = qJD(3) + qJD(4);
t116 = sin(qJ(5));
t120 = cos(qJ(5));
t177 = qJD(1) * t118;
t163 = t117 * t177;
t176 = qJD(1) * t122;
t58 = -t121 * t176 + t163;
t44 = -t120 * t109 - t116 * t58;
t217 = t214 * t44;
t123 = cos(qJ(1));
t107 = g(2) * t123;
t119 = sin(qJ(1));
t199 = g(1) * t119;
t206 = -t107 + t199;
t115 = qJ(3) + qJ(4);
t105 = cos(t115);
t165 = t105 * t199;
t108 = qJDD(3) + qJDD(4);
t124 = -pkin(1) - pkin(6);
t79 = t124 * qJD(1) + qJD(2);
t50 = -pkin(7) * t177 + t118 * t79;
t187 = t121 * t50;
t51 = -pkin(7) * t176 + t122 * t79;
t48 = qJD(3) * pkin(3) + t51;
t29 = t117 * t48 + t187;
t175 = qJD(3) * t118;
t166 = t122 * qJDD(1);
t169 = qJD(1) * qJD(3);
t211 = t118 * t169 - t166;
t77 = t124 * qJDD(1) + qJDD(2);
t65 = t122 * t77;
t34 = qJDD(3) * pkin(3) + t211 * pkin(7) - t79 * t175 + t65;
t174 = qJD(3) * t122;
t161 = t122 * t169;
t167 = t118 * qJDD(1);
t212 = t161 + t167;
t36 = -t212 * pkin(7) + t118 * t77 + t79 * t174;
t204 = t29 * qJD(4) + t117 * t36 - t121 * t34;
t7 = -t108 * pkin(4) + t204;
t216 = t7 + t165;
t104 = sin(t115);
t215 = g(3) * t104 + t105 * t107;
t213 = qJDD(2) - t206;
t110 = qJDD(1) * qJ(2);
t111 = qJD(1) * qJD(2);
t145 = g(1) * t123 + g(2) * t119;
t131 = -t145 + 0.2e1 * t111;
t210 = 0.2e1 * t110 + t131;
t37 = -t58 * pkin(4) + t59 * pkin(8);
t97 = t117 * pkin(3) + pkin(8);
t209 = (pkin(3) * t176 + qJD(5) * t97 + t37) * t214;
t208 = (pkin(8) * qJD(5) + t37) * t214;
t207 = t59 * t109;
t186 = pkin(1) * qJDD(1);
t205 = t186 - t213;
t173 = qJD(4) * t117;
t203 = (qJD(4) * t48 + t36) * t121 + t117 * t34 - t50 * t173;
t71 = pkin(3) * t177 + qJD(1) * qJ(2);
t31 = t59 * pkin(4) + t58 * pkin(8) + t71;
t158 = t108 * pkin(8) + qJD(5) * t31 + t203;
t193 = pkin(7) - t124;
t69 = t193 * t118;
t70 = t193 * t122;
t42 = -t117 * t69 + t121 * t70;
t61 = t193 * t175;
t62 = qJD(3) * t70;
t17 = -t42 * qJD(4) + t117 * t61 - t121 * t62;
t137 = -qJD(4) * t163 - t211 * t117;
t148 = t109 * t122;
t24 = (qJD(1) * t148 + t167) * t121 + t137;
t22 = qJDD(5) + t24;
t190 = t117 * t50;
t28 = t121 * t48 - t190;
t25 = -t109 * pkin(4) - t28;
t63 = t117 * t118 - t121 * t122;
t93 = t118 * pkin(3) + qJ(2);
t38 = t64 * pkin(4) + t63 * pkin(8) + t93;
t172 = qJD(4) * t121;
t40 = -t117 * t174 - t118 * t172 - t121 * t175 - t122 * t173;
t43 = -t117 * t70 - t121 * t69;
t202 = -(qJD(5) * t38 + t17) * t214 - t158 * t64 - t43 * t22 + t25 * t40 - t7 * t63;
t95 = g(3) * t105;
t198 = t25 * t59;
t197 = t25 * t63;
t196 = t38 * t22;
t195 = t214 * t58;
t194 = t58 * t59;
t192 = -t63 * t108 + t40 * t109;
t191 = t116 * t22;
t170 = qJD(5) * t120;
t171 = qJD(5) * t116;
t23 = -t117 * t167 + t121 * t166 - t207;
t12 = t116 * t108 + t109 * t170 + t120 * t23 + t58 * t171;
t189 = t12 * t116;
t188 = t120 * t22;
t185 = t119 * t116;
t184 = t119 * t120;
t183 = t123 * t116;
t182 = t123 * t120;
t126 = qJD(1) ^ 2;
t181 = t126 * qJ(2);
t114 = t122 ^ 2;
t179 = t118 ^ 2 - t114;
t125 = qJD(3) ^ 2;
t178 = -t125 - t126;
t80 = pkin(3) * t174 + qJD(2);
t168 = qJDD(3) * t118;
t164 = t63 * t171;
t26 = t109 * pkin(8) + t29;
t49 = t212 * pkin(3) + t110 + t111;
t9 = t24 * pkin(4) - t23 * pkin(8) + t49;
t159 = qJD(5) * t26 - t9;
t157 = -t120 * t108 + t116 * t23;
t151 = t120 * t214;
t149 = qJD(5) * t64 + qJD(1);
t32 = t117 * t51 + t187;
t146 = pkin(3) * t173 - t32;
t46 = t116 * t109 - t120 * t58;
t144 = -t63 * t12 + t40 * t46;
t143 = -t214 * t40 + t22 * t63;
t11 = t116 * t31 + t120 * t26;
t142 = -t11 * t58 + t216 * t116 + t25 * t170;
t41 = -t117 * t175 - t118 * t173 + t121 * t148;
t141 = -t64 * t108 - t41 * t109;
t10 = -t116 * t26 + t120 * t31;
t140 = t10 * t58 + t215 * t120 + t25 * t171;
t138 = -t158 + t95;
t135 = 0.2e1 * qJ(2) * t169 + qJDD(3) * t124;
t134 = -pkin(8) * t22 + t214 * t28 + t198;
t132 = -t181 - t206;
t33 = t121 * t51 - t190;
t130 = -t97 * t22 + t198 + (-pkin(3) * t172 + t33) * t214;
t129 = -t124 * t125 + t210;
t128 = t206 * t104 + t71 * t59 - t203 + t95;
t127 = t71 * t58 - t165 - t204 + t215;
t103 = qJDD(3) * t122;
t98 = -t121 * pkin(3) - pkin(4);
t55 = t104 * t182 - t185;
t54 = t104 * t183 + t184;
t53 = t104 * t184 + t183;
t52 = -t104 * t185 + t182;
t27 = t58 ^ 2 - t59 ^ 2;
t18 = t43 * qJD(4) - t117 * t62 - t121 * t61;
t16 = -t58 * t109 + (-t109 * t176 - t167) * t121 - t137;
t15 = t23 + t207;
t14 = t41 * pkin(4) - t40 * pkin(8) + t80;
t13 = t46 * qJD(5) + t157;
t8 = t120 * t9;
t4 = t151 * t214 + t46 * t58 + t191;
t3 = -t116 * t214 ^ 2 - t44 * t58 + t188;
t2 = t46 * t151 + t189;
t1 = (t12 - t217) * t120 + (-t214 * t46 - t13) * t116;
t5 = [qJDD(1), t206, t145, -0.2e1 * t186 + t213, t210, t205 * pkin(1) + (t131 + t110) * qJ(2), t114 * qJDD(1) - 0.2e1 * t118 * t161, -0.2e1 * t118 * t166 + 0.2e1 * t179 * t169, -t125 * t118 + t103, -t125 * t122 - t168, 0, t129 * t118 + t135 * t122, -t118 * t135 + t122 * t129, -t23 * t63 - t58 * t40, -t23 * t64 + t63 * t24 - t40 * t59 + t58 * t41, t192, t141, 0, -t104 * t145 - t42 * t108 - t18 * t109 + t93 * t24 + t71 * t41 + t49 * t64 + t80 * t59, -t105 * t145 - t43 * t108 - t17 * t109 + t93 * t23 + t71 * t40 - t49 * t63 - t80 * t58, t120 * t144 + t164 * t46, (-t116 * t46 - t120 * t44) * t40 + (t189 + t120 * t13 + (-t116 * t44 + t120 * t46) * qJD(5)) * t63, t12 * t64 - t120 * t143 + t164 * t214 + t46 * t41, t170 * t214 * t63 + t116 * t143 - t13 * t64 - t44 * t41, t214 * t41 + t22 * t64, -g(1) * t55 - g(2) * t53 + t10 * t41 + t42 * t13 + t18 * t44 + t8 * t64 + (t14 * t214 + t196 + (-t214 * t43 - t26 * t64 - t197) * qJD(5)) * t120 + t202 * t116, g(1) * t54 - g(2) * t52 - t11 * t41 + t42 * t12 + t18 * t46 + (-(-qJD(5) * t43 + t14) * t214 - t196 + t159 * t64 + qJD(5) * t197) * t116 + t202 * t120; 0, 0, 0, qJDD(1), -t126, -t181 - t205, 0, 0, 0, 0, 0, t178 * t118 + t103, t178 * t122 - t168, 0, 0, 0, 0, 0, -qJD(1) * t59 + t192, qJD(1) * t58 + t141, 0, 0, 0, 0, 0, -t64 * t191 + t63 * t13 - t40 * t44 + (-t116 * t41 - t120 * t149) * t214, -t64 * t188 + (t116 * t149 - t120 * t41) * t214 - t144; 0, 0, 0, 0, 0, 0, t122 * t126 * t118, -t179 * t126, t166, -t167, qJDD(3), g(3) * t118 + t122 * t132 + t65, g(3) * t122 + (-t132 - t77) * t118, -t194, t27, t15, t16, t108, t32 * t109 + (t108 * t121 - t109 * t173 - t59 * t176) * pkin(3) + t127, t33 * t109 + (-t108 * t117 - t109 * t172 + t58 * t176) * pkin(3) + t128, t2, t1, t4, t3, t195, t98 * t13 + t146 * t44 + (-t216 - t209) * t120 + t130 * t116 + t140, t98 * t12 + t146 * t46 + t130 * t120 + (-t215 + t209) * t116 + t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t194, t27, t15, t16, t108, t29 * t109 + t127, t28 * t109 + t128, t2, t1, t4, t3, t195, -pkin(4) * t13 - t29 * t44 + t134 * t116 + (-t216 - t208) * t120 + t140, -pkin(4) * t12 - t29 * t46 + t134 * t120 + (-t215 + t208) * t116 + t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46 * t44, -t44 ^ 2 + t46 ^ 2, t12 + t217, -t157 + (-qJD(5) + t214) * t46, t22, -g(1) * t52 - g(2) * t54 + t11 * t214 + t116 * t138 - t26 * t170 - t25 * t46 + t8, g(1) * t53 - g(2) * t55 + t10 * t214 + t116 * t159 + t120 * t138 + t25 * t44;];
tau_reg = t5;
