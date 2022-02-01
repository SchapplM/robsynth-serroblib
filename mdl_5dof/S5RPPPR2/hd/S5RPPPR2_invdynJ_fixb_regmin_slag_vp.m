% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPPR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 08:59:46
% EndTime: 2022-01-23 08:59:51
% DurationCPUTime: 2.03s
% Computational Cost: add. (1409->288), mult. (3576->444), div. (0->0), fcn. (2868->10), ass. (0->163)
t120 = sin(pkin(9));
t123 = cos(pkin(9));
t125 = cos(pkin(7));
t193 = t125 * t123;
t122 = sin(pkin(7));
t124 = cos(pkin(8));
t196 = t122 * t124;
t70 = t120 * t196 + t193;
t61 = t70 * qJD(1);
t55 = qJD(5) + t61;
t127 = sin(qJ(1));
t129 = cos(qJ(1));
t156 = g(1) * t127 - g(2) * t129;
t128 = cos(qJ(5));
t121 = sin(pkin(8));
t126 = sin(qJ(5));
t198 = t121 * t126;
t72 = t120 * t122 + t124 * t193;
t138 = t125 * t198 + t128 * t72;
t197 = t121 * t128;
t76 = t123 * t197 - t124 * t126;
t213 = t127 * t138 - t129 * t76;
t202 = qJDD(1) * pkin(1);
t212 = t202 + t156;
t41 = t125 * t197 - t126 * t72;
t75 = t123 * t198 + t124 * t128;
t211 = -t127 * t75 + t41 * t129;
t210 = -t55 + qJD(5);
t178 = qJ(2) * qJDD(1);
t176 = qJD(1) * qJD(4);
t172 = t125 * qJDD(1);
t158 = t124 * t172;
t177 = qJD(1) * qJD(2);
t163 = t125 * t177;
t180 = qJD(3) * t122;
t157 = qJ(3) * t122 + pkin(1);
t85 = pkin(2) * t125 + t157;
t49 = -qJD(1) * t180 - qJDD(1) * t85 + qJDD(2);
t23 = qJ(2) * t158 + t121 * t49 + t124 * t163;
t18 = (-qJ(4) * qJDD(1) - t176) * t125 + t23;
t144 = pkin(3) * t121 - qJ(4) * t124;
t173 = t122 * qJDD(1);
t81 = qJ(2) * t173 + t122 * t177 + qJDD(3);
t27 = (qJDD(1) * t144 - t124 * t176) * t122 + t81;
t4 = t120 * t27 + t123 * t18;
t182 = qJD(1) * t125;
t169 = qJ(2) * t182;
t69 = -qJD(1) * t85 + qJD(2);
t38 = t121 * t69 + t124 * t169;
t30 = -qJ(4) * t182 + t38;
t184 = qJD(1) * t122;
t97 = qJ(2) * t184 + qJD(3);
t47 = t144 * t184 + t97;
t11 = t120 * t47 + t123 * t30;
t204 = qJ(2) * t125;
t52 = -t121 * t85 + t124 * t204;
t44 = -qJ(4) * t125 + t52;
t140 = qJ(2) + t144;
t57 = t140 * t122;
t20 = t120 * t57 + t123 * t44;
t59 = t70 * qJDD(1);
t54 = qJDD(5) + t59;
t208 = t126 * t54;
t206 = t128 * t54;
t203 = qJD(5) * t55;
t117 = t122 ^ 2;
t130 = qJD(1) ^ 2;
t201 = t117 * t130;
t200 = t120 * t121;
t199 = t121 * t122;
t195 = t122 * t127;
t194 = t122 * t129;
t192 = t125 * t130;
t191 = t127 * t121;
t190 = t127 * t124;
t189 = t129 * t121;
t188 = t129 * t124;
t116 = t121 ^ 2;
t187 = -t124 ^ 2 - t116;
t186 = t125 ^ 2 + t117;
t185 = qJD(1) * t121;
t183 = qJD(1) * t124;
t181 = qJD(2) * t125;
t179 = qJD(5) * t126;
t175 = qJDD(1) * t121;
t174 = qJDD(1) * t124;
t171 = t116 * t201;
t161 = t121 * t173;
t2 = pkin(6) * t161 + t4;
t159 = t121 * t172;
t22 = -qJ(2) * t159 - t121 * t163 + t124 * t49;
t21 = pkin(3) * t172 + qJDD(4) - t22;
t160 = t123 * t174;
t98 = t120 * t172;
t60 = t122 * t160 - t98;
t6 = pkin(4) * t59 - pkin(6) * t60 + t21;
t170 = -t126 * t2 + t128 * t6;
t168 = t121 * t184;
t167 = t126 * t185;
t166 = t128 * t185;
t165 = t120 * t182;
t164 = t123 * t184;
t162 = t116 * t173;
t37 = -t121 * t169 + t124 * t69;
t155 = t126 * t60 - t128 * t161;
t154 = t55 ^ 2;
t51 = -t121 * t204 - t124 * t85;
t153 = t186 * t130;
t152 = 0.2e1 * t186;
t151 = t122 * t166;
t77 = t125 * t191 + t188;
t79 = t125 * t189 - t190;
t150 = -g(1) * t77 + g(2) * t79;
t46 = pkin(3) * t125 - t51;
t149 = g(1) * t129 + g(2) * t127;
t148 = t126 * t6 + t128 * t2;
t29 = pkin(3) * t182 + qJD(4) - t37;
t64 = t124 * t164 - t165;
t7 = pkin(4) * t61 - pkin(6) * t64 + t29;
t9 = pkin(6) * t168 + t11;
t147 = t126 * t9 - t128 * t7;
t146 = -t126 * t7 - t128 * t9;
t74 = -t121 * t180 + t124 * t181;
t3 = -t120 * t18 + t123 * t27;
t10 = -t120 * t30 + t123 * t47;
t19 = -t120 * t44 + t123 * t57;
t71 = -t120 * t125 + t123 * t196;
t14 = pkin(4) * t70 - pkin(6) * t71 + t46;
t16 = pkin(6) * t199 + t20;
t143 = -t126 * t16 + t128 * t14;
t142 = t126 * t14 + t128 * t16;
t141 = (-qJD(4) * t124 + qJD(2)) * t122;
t139 = t122 * t197 - t126 * t71;
t40 = t122 * t198 + t128 * t71;
t73 = t121 * t181 + t124 * t180;
t137 = qJD(1) * t73 - qJDD(1) * t51 - t22;
t136 = qJD(1) * t74 + qJDD(1) * t52 + t23;
t12 = qJD(5) * t151 + t126 * t161 + t128 * t60 - t179 * t64;
t135 = t152 * t177;
t33 = t122 * t167 + t128 * t64;
t134 = t76 * t55;
t110 = qJDD(2) - t202;
t133 = -t110 + t212;
t132 = t122 * t81 + (t177 + t178) * t117;
t114 = g(3) * t125;
t113 = t129 * qJ(2);
t111 = t127 * qJ(2);
t80 = t125 * t188 + t191;
t78 = -t125 * t190 + t189;
t66 = t72 * qJD(1);
t63 = t124 * t165 - t164;
t56 = (pkin(3) * t124 + qJ(4) * t121 + pkin(2)) * t125 + t157;
t53 = -qJD(4) * t125 + t74;
t35 = t40 * qJD(5);
t34 = t139 * qJD(5);
t31 = t126 * t64 - t151;
t26 = t120 * t141 + t123 * t53;
t25 = t120 * t53 - t123 * t141;
t15 = -pkin(4) * t199 - t19;
t13 = qJD(5) * t33 + t155;
t8 = -pkin(4) * t168 - t10;
t1 = -pkin(4) * t161 - t3;
t5 = [qJDD(1), t156, t149, t133 * t125, -t133 * t122, t152 * t178 + t135 - t149, -t110 * pkin(1) - g(1) * (-pkin(1) * t127 + t113) - g(2) * (pkin(1) * t129 + t111) + (t178 * t186 + t135) * qJ(2), -g(1) * t78 - g(2) * t80 + t121 * t132 + t125 * t137, t124 * t132 + t125 * t136 + t150, (-t121 * t136 + t124 * t137 + t156) * t122, t23 * t52 + t38 * t74 + t22 * t51 - t37 * t73 - g(1) * (-t127 * t85 + t113) - g(2) * (t129 * t85 + t111) + (qJ(2) * t81 + qJD(2) * t97) * t122, t73 * t61 + t46 * t59 + t21 * t70 - g(1) * (-t120 * t195 + t123 * t78) - g(2) * (t120 * t194 + t123 * t80) + (-qJD(1) * t25 + qJDD(1) * t19 + t3) * t199, t73 * t64 + t46 * t60 + t21 * t71 - g(1) * (-t120 * t78 - t123 * t195) - g(2) * (-t120 * t80 + t123 * t194) + (-qJD(1) * t26 - qJDD(1) * t20 - t4) * t199, -t19 * t60 - t20 * t59 + t25 * t64 - t26 * t61 - t3 * t71 - t4 * t70 - t150, t4 * t20 + t11 * t26 + t3 * t19 - t10 * t25 + t21 * t46 + t29 * t73 - g(1) * (-t127 * t56 + t129 * t140) - g(2) * (t127 * t140 + t129 * t56), t12 * t40 + t33 * t34, t12 * t139 - t13 * t40 - t31 * t34 - t33 * t35, t12 * t70 + t34 * t55 + t40 * t54, -t13 * t70 + t139 * t54 - t35 * t55, t54 * t70, (-t126 * t26 + t128 * t73) * t55 + t143 * t54 + t170 * t70 + t25 * t31 + t15 * t13 - t1 * t139 + t8 * t35 + g(1) * t213 - g(2) * ((t123 * t191 + t129 * t72) * t128 + t79 * t126) + (-t142 * t55 + t146 * t70) * qJD(5), -(t126 * t73 + t128 * t26) * t55 - t142 * t54 - t148 * t70 + t25 * t33 + t15 * t12 + t1 * t40 + t8 * t34 - g(1) * (-t127 * t41 - t129 * t75) - g(2) * t211 + (-t143 * t55 + t147 * t70) * qJD(5); 0, 0, 0, -t172, t173, -t153, -qJ(2) * t153 + qJDD(2) - t212, -t121 * t153 - t158, -t124 * t153 + t159, t187 * t173, t23 * t121 + t22 * t124 + (-t122 * t97 + (t121 * t37 - t124 * t38) * t125) * qJD(1) - t156, -t120 * t162 - t124 * t59 + (t122 * t63 - t125 * t61) * t185, -t123 * t162 - t124 * t60 + (t122 * t66 - t125 * t64) * t185, t66 * t61 - t63 * t64 + (t120 * t60 - t123 * t59) * t121, t10 * t63 - t11 * t66 - t21 * t124 + (-t120 * t3 + t123 * t4 - t182 * t29) * t121 - t156, 0, 0, 0, 0, 0, -t75 * t54 + t13 * t200 - (t125 * t166 - t126 * t66) * t55 - t63 * t31 - qJD(5) * t134, -t76 * t54 + t12 * t200 + (t125 * t167 + t128 * t66) * t55 - t63 * t33 + t75 * t203; 0, 0, 0, 0, 0, 0, 0, (-t124 * t192 + t175) * t122, (t121 * t192 + t174) * t122, t187 * t201, t114 + ((t121 * t38 + t124 * t37) * qJD(1) - t149) * t122 + t81, -t120 * t171 + (t123 * t175 - t183 * t61) * t122, -t123 * t171 + (-t120 * t175 - t183 * t64) * t122, -t120 * t59 - t123 * t60 + (t120 * t64 - t123 * t61) * t168, t4 * t120 + t3 * t123 + t114 + ((-t124 * t29 + (-t10 * t120 + t11 * t123) * t121) * qJD(1) - t149) * t122, 0, 0, 0, 0, 0, -t123 * t13 + (-t128 * t203 - t208) * t120 + (t200 * t31 - t55 * t75) * t184, -t123 * t12 + (t179 * t55 - t206) * t120 + (t200 * t33 - t134) * t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t168 * t64 + t59, -t98 + (-t185 * t61 + t160) * t122, -t61 ^ 2 - t64 ^ 2, -g(1) * t79 - g(2) * t77 - g(3) * t199 + t10 * t64 + t11 * t61 + t21, 0, 0, 0, 0, 0, -t126 * t154 - t64 * t31 + t206, -t128 * t154 - t64 * t33 - t208; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33 * t31, -t31 ^ 2 + t33 ^ 2, t31 * t55 + t12, -t210 * t33 - t155, t54, -t8 * t33 - g(1) * t211 - g(2) * (-(-t123 * t189 + t127 * t72) * t126 + t77 * t128) - g(3) * t139 + t170 + t210 * t146, t8 * t31 - g(1) * (-t127 * t76 - t129 * t138) + g(2) * t213 + g(3) * t40 - t148 + t210 * t147;];
tau_reg = t5;
