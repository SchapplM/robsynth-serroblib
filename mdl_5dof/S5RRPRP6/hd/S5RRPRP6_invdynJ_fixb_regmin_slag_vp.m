% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRP6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [5x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRP6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:58:29
% EndTime: 2019-12-31 19:58:35
% DurationCPUTime: 1.99s
% Computational Cost: add. (2650->309), mult. (6222->410), div. (0->0), fcn. (4381->10), ass. (0->164)
t123 = sin(qJ(4));
t119 = sin(pkin(8));
t124 = sin(qJ(2));
t173 = qJD(1) * t124;
t120 = cos(pkin(8));
t127 = cos(qJ(2));
t180 = t120 * t127;
t80 = qJD(1) * t180 - t119 * t173;
t73 = qJD(4) - t80;
t153 = t123 * t73;
t126 = cos(qJ(4));
t92 = t119 * t127 + t120 * t124;
t82 = t92 * qJD(1);
t67 = qJD(2) * t123 + t126 * t82;
t214 = t67 * t153;
t194 = t127 * pkin(2);
t109 = pkin(1) + t194;
t213 = -pkin(7) * t92 - t109;
t125 = sin(qJ(1));
t128 = cos(qJ(1));
t212 = g(1) * t125 - g(2) * t128;
t149 = g(1) * t128 + g(2) * t125;
t116 = qJ(2) + pkin(8);
t111 = cos(t116);
t176 = t126 * t128;
t179 = t123 * t125;
t74 = t111 * t179 + t176;
t177 = t125 * t126;
t178 = t123 * t128;
t76 = -t111 * t178 + t177;
t211 = -g(1) * t76 + g(2) * t74;
t105 = pkin(2) * t119 + pkin(7);
t110 = sin(t116);
t135 = g(3) * t111 - t110 * t149;
t122 = -qJ(3) - pkin(6);
t156 = qJD(2) * t122;
t79 = -qJD(3) * t124 + t127 * t156;
t97 = t122 * t124;
t50 = qJDD(2) * pkin(2) + qJD(1) * t79 + qJDD(1) * t97;
t78 = qJD(3) * t127 + t124 * t156;
t98 = t122 * t127;
t57 = t78 * qJD(1) - qJDD(1) * t98;
t21 = -t119 * t57 + t120 * t50;
t19 = -qJDD(2) * pkin(3) - t21;
t210 = qJD(4) * t105 * t73 + t135 + t19;
t169 = qJD(1) * qJD(2);
t159 = t127 * t169;
t160 = t124 * t169;
t99 = t119 * t160;
t141 = t120 * t159 - t99;
t166 = pkin(2) * t160 + qJDD(3);
t167 = t127 * qJDD(1);
t168 = t124 * qJDD(1);
t145 = -t119 * t168 + t120 * t167;
t56 = -qJD(2) * t82 + t145;
t16 = -t56 * pkin(3) - t141 * pkin(7) + t213 * qJDD(1) + t166;
t13 = t126 * t16;
t96 = -t109 * qJD(1) + qJD(3);
t29 = -pkin(3) * t80 - pkin(7) * t82 + t96;
t95 = qJD(1) * t98;
t185 = t120 * t95;
t186 = qJD(2) * pkin(2);
t94 = qJD(1) * t97;
t88 = t94 + t186;
t55 = t119 * t88 - t185;
t45 = qJD(2) * pkin(7) + t55;
t15 = t123 * t29 + t126 * t45;
t22 = t119 * t50 + t120 * t57;
t20 = qJDD(2) * pkin(7) + t22;
t132 = qJDD(1) * t92 + t141;
t170 = t126 * qJD(2);
t172 = qJD(4) * t123;
t24 = -qJD(4) * t170 - t123 * qJDD(2) - t126 * t132 + t82 * t172;
t81 = t92 * qJD(2);
t51 = qJD(1) * t81 + qJDD(4) - t145;
t1 = pkin(4) * t51 + qJ(5) * t24 - qJD(4) * t15 - qJD(5) * t67 - t123 * t20 + t13;
t65 = t123 * t82 - t170;
t9 = -qJ(5) * t65 + t15;
t209 = t73 * t9 + t1;
t208 = t67 ^ 2;
t14 = -t123 * t45 + t126 * t29;
t8 = -qJ(5) * t67 + t14;
t6 = pkin(4) * t73 + t8;
t207 = -t8 + t6;
t175 = qJ(5) + t105;
t152 = qJD(4) * t175;
t181 = qJ(5) * t126;
t37 = pkin(2) * t173 + pkin(3) * t82 - pkin(7) * t80;
t33 = t126 * t37;
t85 = t119 * t95;
t59 = t120 * t94 + t85;
t203 = -pkin(4) * t82 - t126 * t152 + t80 * t181 - t33 + (-qJD(5) + t59) * t123;
t202 = pkin(2) * t120;
t201 = pkin(4) * t123;
t197 = g(3) * t110;
t195 = g(3) * t127;
t193 = t65 * t80;
t192 = t65 * t82;
t191 = t67 * t82;
t171 = qJD(4) * t126;
t131 = -t126 * qJDD(2) + t123 * t132;
t25 = qJD(4) * t67 + t131;
t190 = -t123 * t25 - t65 * t171;
t189 = t123 * t37 + t126 * t59;
t91 = t119 * t124 - t180;
t53 = pkin(3) * t91 + t213;
t63 = t119 * t97 - t120 * t98;
t60 = t126 * t63;
t188 = t123 * t53 + t60;
t182 = qJ(5) * t123;
t187 = qJD(5) * t126 - t123 * t152 + t80 * t182 - t189;
t184 = t123 * t24;
t183 = t123 * t51;
t117 = t124 ^ 2;
t174 = -t127 ^ 2 + t117;
t36 = t119 * t79 + t120 * t78;
t164 = t124 * t186;
t84 = t91 * qJD(2);
t38 = pkin(3) * t81 + pkin(7) * t84 + t164;
t165 = t123 * t38 + t126 * t36 + t53 * t171;
t162 = t92 * t172;
t161 = t92 * t171;
t108 = pkin(4) * t126 + pkin(3);
t158 = -t122 + t201;
t35 = t119 * t78 - t120 * t79;
t58 = t119 * t94 - t185;
t62 = -t119 * t98 - t120 * t97;
t54 = t120 * t88 + t85;
t155 = -qJD(4) * t29 - t20;
t154 = t126 * t73;
t151 = -t45 * t171 + t13;
t139 = t123 * t16 + t126 * t20 + t29 * t171 - t45 * t172;
t2 = -qJ(5) * t25 - qJD(5) * t65 + t139;
t150 = -t73 * t6 + t2;
t44 = -qJD(2) * pkin(3) - t54;
t147 = t19 * t92 - t44 * t84;
t146 = t51 * t92 - t73 * t84;
t144 = qJ(5) * t84 - qJD(5) * t92;
t121 = -qJ(5) - pkin(7);
t143 = t108 * t111 - t110 * t121;
t142 = t126 * t51 + t80 * t153 - t73 * t172;
t140 = -0.2e1 * pkin(1) * t169 - pkin(6) * qJDD(2);
t138 = -t105 * t51 + t73 * t44;
t136 = -t109 * qJDD(1) + t166;
t5 = pkin(4) * t25 + qJDD(5) + t19;
t129 = qJD(2) ^ 2;
t134 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t129 + t212;
t130 = qJD(1) ^ 2;
t133 = pkin(1) * t130 - pkin(6) * qJDD(1) + t149;
t106 = -pkin(3) - t202;
t100 = t128 * t109;
t90 = t175 * t126;
t89 = t175 * t123;
t77 = t111 * t176 + t179;
t75 = -t111 * t177 + t178;
t64 = t65 ^ 2;
t43 = t126 * t53;
t34 = t126 * t38;
t26 = pkin(4) * t65 + qJD(5) + t44;
t17 = -t92 * t182 + t188;
t10 = pkin(4) * t91 - t123 * t63 - t92 * t181 + t43;
t4 = -qJ(5) * t161 + (-qJD(4) * t63 + t144) * t123 + t165;
t3 = pkin(4) * t81 - t123 * t36 + t34 + t144 * t126 + (-t60 + (qJ(5) * t92 - t53) * t123) * qJD(4);
t7 = [qJDD(1), t212, t149, qJDD(1) * t117 + 0.2e1 * t124 * t159, 0.2e1 * t124 * t167 - 0.2e1 * t174 * t169, qJDD(2) * t124 + t127 * t129, qJDD(2) * t127 - t124 * t129, 0, t124 * t140 + t127 * t134, -t124 * t134 + t127 * t140, t62 * t132 - t21 * t92 - t22 * t91 + t35 * t82 + t36 * t80 + t54 * t84 - t55 * t81 + t63 * t56 - t149, t22 * t63 + t55 * t36 - t21 * t62 - t54 * t35 - t136 * t109 + t96 * t164 - g(1) * (-t109 * t125 - t122 * t128) - g(2) * (-t122 * t125 + t100), -t67 * t162 + (-t24 * t92 - t67 * t84) * t126, -(-t123 * t67 - t126 * t65) * t84 + (t184 - t126 * t25 + (t123 * t65 - t126 * t67) * qJD(4)) * t92, t126 * t146 - t162 * t73 - t24 * t91 + t67 * t81, -t123 * t146 - t161 * t73 - t25 * t91 - t65 * t81, t51 * t91 + t73 * t81, (-t171 * t63 + t34) * t73 + t43 * t51 + t151 * t91 + t14 * t81 + t35 * t65 + t62 * t25 + t44 * t161 - g(1) * t75 - g(2) * t77 + ((-qJD(4) * t53 - t36) * t73 - t63 * t51 + t155 * t91 + t147) * t123, -(-t172 * t63 + t165) * t73 - t188 * t51 - t139 * t91 - t15 * t81 + t35 * t67 - t62 * t24 - t44 * t162 - g(1) * t74 - g(2) * t76 + t147 * t126, t10 * t24 - t17 * t25 - t3 * t67 - t4 * t65 - (-t123 * t9 - t126 * t6) * t84 + t212 * t110 + (-t1 * t126 - t123 * t2 + (t123 * t6 - t126 * t9) * qJD(4)) * t92, t2 * t17 + t9 * t4 + t1 * t10 + t6 * t3 + t5 * (t201 * t92 + t62) + t26 * ((-t123 * t84 + t161) * pkin(4) + t35) - g(2) * t100 + (-g(1) * t158 - g(2) * t143) * t128 + (-g(1) * (-t109 - t143) - g(2) * t158) * t125; 0, 0, 0, -t124 * t130 * t127, t174 * t130, t168, t167, qJDD(2), t124 * t133 - t195, g(3) * t124 + t127 * t133, (t55 - t58) * t82 + (-t59 + t54) * t80 + (t119 * t56 + (-t119 * t167 + t99 + (-t159 - t168) * t120) * t120) * pkin(2), t54 * t58 - t55 * t59 + (-t195 + t119 * t22 + t120 * t21 + (-qJD(1) * t96 + t149) * t124) * pkin(2), t154 * t67 - t184, (-t24 + t193) * t126 - t214 + t190, t154 * t73 + t183 - t191, t142 + t192, -t73 * t82, t106 * t25 - t14 * t82 - t33 * t73 - t58 * t65 + (t59 * t73 + t138) * t123 - t210 * t126, -t106 * t24 + t210 * t123 + t138 * t126 + t15 * t82 + t189 * t73 - t58 * t67, -t149 * t111 - t209 * t123 + t150 * t126 - t187 * t65 - t203 * t67 - t24 * t89 - t25 * t90 - t197, t2 * t90 - t1 * t89 + t5 * (-t108 - t202) - g(3) * (t143 + t194) + t187 * t9 + t203 * t6 + (pkin(4) * t153 - t58) * t26 + t149 * (pkin(2) * t124 + t108 * t110 + t111 * t121); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80 ^ 2 - t82 ^ 2, t54 * t82 - t55 * t80 + t136 - t212, 0, 0, 0, 0, 0, t142 - t192, -t126 * t73 ^ 2 - t183 - t191, (t24 + t193) * t126 + t214 + t190, t150 * t123 + t209 * t126 - t26 * t82 - t212; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67 * t65, -t64 + t208, t65 * t73 - t24, -t131 + (-qJD(4) + t73) * t67, t51, t15 * t73 - t44 * t67 + (t155 + t197) * t123 + t151 + t211, g(1) * t77 - g(2) * t75 + t126 * t197 + t14 * t73 + t44 * t65 - t139, pkin(4) * t24 - t207 * t65, t207 * t9 + (t123 * t197 - t26 * t67 + t1 + t211) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64 - t208, t6 * t67 + t65 * t9 + t135 + t5;];
tau_reg = t7;
