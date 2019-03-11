% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
% 
% Output:
% tau_reg [6x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPPRR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:34:02
% EndTime: 2019-03-09 01:34:06
% DurationCPUTime: 1.80s
% Computational Cost: add. (2111->292), mult. (4058->376), div. (0->0), fcn. (2974->12), ass. (0->163)
t124 = sin(qJ(6));
t119 = sin(pkin(10));
t125 = sin(qJ(5));
t192 = t119 * t125;
t127 = cos(qJ(5));
t121 = cos(pkin(10));
t184 = qJD(1) * t121;
t80 = t127 * t184;
t62 = qJD(1) * t192 - t80;
t189 = qJD(6) - t62;
t126 = cos(qJ(6));
t218 = t126 * t189;
t174 = t121 * qJDD(1);
t175 = t119 * qJDD(1);
t152 = t125 * t175 - t127 * t174;
t68 = t119 * t127 + t121 * t125;
t65 = t68 * qJD(5);
t36 = qJD(1) * t65 + t152;
t33 = -qJDD(6) + t36;
t220 = t124 * t33 - t189 * t218;
t180 = t126 * qJD(5);
t63 = t68 * qJD(1);
t45 = -t124 * t63 - t180;
t219 = t189 * t45;
t117 = pkin(10) + qJ(5);
t101 = sin(t117);
t102 = cos(t117);
t120 = sin(pkin(9));
t122 = cos(pkin(9));
t208 = sin(qJ(1));
t209 = cos(qJ(1));
t67 = -t208 * t120 - t209 * t122;
t212 = g(1) * t67;
t69 = t209 * t120 - t208 * t122;
t159 = g(2) * t69 + t212;
t134 = g(3) * t102 - t159 * t101;
t100 = t121 * qJD(3);
t183 = qJD(1) * t122;
t128 = -pkin(1) - pkin(2);
t78 = t128 * qJD(1) + qJD(2);
t59 = qJ(2) * t183 + t120 * t78;
t49 = -qJD(1) * qJ(4) + t59;
t28 = t100 + (pkin(7) * qJD(1) - t49) * t119;
t39 = t119 * qJD(3) + t121 * t49;
t29 = -pkin(7) * t184 + t39;
t10 = t125 * t28 + t127 * t29;
t176 = qJDD(1) * t122;
t77 = t128 * qJDD(1) + qJDD(2);
t201 = qJ(2) * t176 + t120 * t77;
t178 = qJD(1) * qJD(2);
t86 = t122 * t178;
t43 = t86 + t201;
t37 = -qJDD(1) * qJ(4) - qJD(1) * qJD(4) + t43;
t98 = t121 * qJDD(3);
t18 = t98 + (pkin(7) * qJDD(1) - t37) * t119;
t22 = t119 * qJDD(3) + t121 * t37;
t19 = -pkin(7) * t174 + t22;
t148 = t125 * t19 - t127 * t18;
t2 = -qJDD(5) * pkin(5) + t10 * qJD(5) + t148;
t217 = t134 - (-pkin(5) * t63 + t189 * pkin(8)) * t189 - t2;
t216 = qJDD(1) * pkin(3) + qJDD(4);
t170 = qJD(5) * t192;
t190 = t127 * t121;
t64 = -qJD(5) * t190 + t170;
t153 = t189 * t64 + t33 * t68;
t193 = qJD(6) * t68;
t214 = -t153 * t124 + t193 * t218;
t74 = t122 * qJ(2) + t120 * t128;
t70 = -qJ(4) + t74;
t210 = pkin(7) - t70;
t54 = t210 * t119;
t55 = t210 * t121;
t16 = -t125 * t55 - t127 * t54;
t66 = -t190 + t192;
t82 = qJD(2) * t122 - qJD(4);
t11 = -t16 * qJD(5) - t66 * t82;
t146 = t125 * t29 - t127 * t28;
t147 = t125 * t18 + t127 * t19;
t185 = qJD(1) * t120;
t58 = -qJ(2) * t185 + t122 * t78;
t48 = qJD(1) * pkin(3) + qJD(4) - t58;
t41 = pkin(4) * t184 + t48;
t15 = -pkin(5) * t62 + pkin(8) * t63 + t41;
t168 = qJDD(5) * pkin(8) - t146 * qJD(5) + qJD(6) * t15 + t147;
t17 = t125 * t54 - t127 * t55;
t73 = -t120 * qJ(2) + t122 * t128;
t71 = pkin(3) - t73;
t61 = t121 * pkin(4) + t71;
t20 = -pkin(5) * t66 + pkin(8) * t68 + t61;
t7 = -qJD(5) * pkin(5) + t146;
t213 = -(qJD(6) * t20 + t11) * t189 + t168 * t66 + t17 * t33 - t2 * t68 + t7 * t64 - t212;
t211 = g(2) * t67;
t207 = g(3) * t101;
t205 = t45 * t63;
t47 = qJD(5) * t124 - t126 * t63;
t204 = t47 * t63;
t203 = t120 * t65 - t66 * t183;
t57 = t66 * t120;
t202 = -qJD(5) * t57 - t122 * t63;
t8 = qJD(5) * pkin(8) + t10;
t200 = qJD(6) * t8;
t199 = t102 * t67;
t198 = t102 * t69;
t181 = qJD(6) * t124;
t135 = qJD(1) * t170 - t68 * qJDD(1);
t35 = -qJD(5) * t80 + t135;
t13 = qJD(6) * t180 + t124 * qJDD(5) + t126 * t35 + t63 * t181;
t195 = t13 * t124;
t194 = pkin(1) * qJDD(1);
t129 = qJD(1) ^ 2;
t191 = t122 * t129;
t188 = t209 * pkin(1) + t208 * qJ(2);
t187 = g(1) * t208 - g(2) * t209;
t186 = t119 ^ 2 + t121 ^ 2;
t182 = qJD(2) * t120;
t179 = qJ(2) * qJDD(1);
t177 = qJDD(1) * t120;
t173 = 0.2e1 * t178;
t172 = t68 * t181;
t171 = t209 * pkin(2) + t188;
t84 = t120 * t178;
t167 = -qJ(2) * t177 + t122 * t77;
t165 = -t126 * qJDD(5) + t124 * t35;
t163 = qJDD(1) * t186;
t161 = qJDD(2) - t194;
t42 = t167 - t84;
t160 = g(1) * t69 - t211;
t158 = -t208 * pkin(1) + t209 * qJ(2);
t157 = -qJD(6) * t122 - t203;
t156 = -t200 + t211;
t155 = -t66 * t13 - t65 * t47;
t14 = t47 * qJD(6) + t165;
t154 = t66 * t14 + t65 * t45;
t21 = -t119 * t37 + t98;
t151 = -t21 * t119 + t22 * t121;
t38 = -t119 * t49 + t100;
t150 = t119 * t38 - t121 * t39;
t149 = t120 * t58 - t122 * t59;
t145 = -qJD(6) * t57 + t185;
t144 = qJD(5) * t64 - qJDD(5) * t68;
t143 = qJD(5) * t65 + qJDD(5) * t66;
t142 = -t126 * t33 + (t124 * t62 - t181) * t189;
t141 = t120 * t129 + t176;
t140 = -g(1) * t198 - t20 * t33;
t40 = -t42 + t216;
t139 = g(1) * t209 + g(2) * t208;
t31 = pkin(4) * t174 + t40;
t138 = -t160 - t167;
t137 = -t208 * pkin(2) + t158;
t136 = pkin(8) * t33 + (-t146 + t7) * t189;
t132 = -g(2) * t198 - t168 - t207;
t131 = t153 * t126 + t172 * t189;
t130 = qJDD(1) * t71 - t160 + t40 + t84;
t56 = t68 * t120;
t26 = t124 * t69 - t126 * t199;
t25 = t124 * t199 + t126 * t69;
t23 = -pkin(5) * t65 - pkin(8) * t64 + t182;
t12 = t17 * qJD(5) + t68 * t82;
t6 = -pkin(5) * t36 - pkin(8) * t35 + t31;
t5 = t126 * t6;
t4 = t124 * t15 + t126 * t8;
t3 = -t124 * t8 + t126 * t15;
t1 = [qJDD(1), t187, t139, -qJDD(2) + t187 + 0.2e1 * t194, -t139 + t173 + 0.2e1 * t179, -t161 * pkin(1) - g(1) * t158 - g(2) * t188 + (t173 + t179) * qJ(2), -qJDD(1) * t73 + t138 + 0.2e1 * t84, qJDD(1) * t74 + t159 + t201 + 0.2e1 * t86, -g(1) * t137 - g(2) * t171 - t149 * qJD(2) + t42 * t73 + t43 * t74, t130 * t121, -t130 * t119, -t186 * t82 * qJD(1) - t70 * t163 - t151 - t159, t40 * t71 + t48 * t182 - g(1) * (t69 * pkin(3) + t67 * qJ(4) + t137) - g(2) * (-pkin(3) * t67 + qJ(4) * t69 + t171) + (t22 * t70 + t39 * t82) * t121 + (-t21 * t70 - t38 * t82) * t119, -t35 * t68 - t63 * t64, t35 * t66 - t36 * t68 + t62 * t64 - t63 * t65, t144, t143, 0, -t12 * qJD(5) - t16 * qJDD(5) - t160 * t102 - t62 * t182 - t31 * t66 - t61 * t36 - t41 * t65, -t11 * qJD(5) - t17 * qJDD(5) + t160 * t101 - t63 * t182 - t31 * t68 + t61 * t35 + t41 * t64, t47 * t172 + (-t13 * t68 + t47 * t64) * t126 (-t124 * t47 - t126 * t45) * t64 + (t195 + t126 * t14 + (-t124 * t45 + t126 * t47) * qJD(6)) * t68, t131 + t155, t154 + t214, -t189 * t65 + t33 * t66, -g(2) * t26 + t12 * t45 + t16 * t14 - t3 * t65 - t5 * t66 + (t23 * t189 + (-t17 * t189 + t66 * t8 - t68 * t7) * qJD(6) + t140) * t126 + t213 * t124, -g(2) * t25 + t12 * t47 + t16 * t13 + t4 * t65 + (-(-qJD(6) * t17 + t23) * t189 + (t6 - t200) * t66 + t7 * t193 - t140) * t124 + t213 * t126; 0, 0, 0, -qJDD(1), -t129, -qJ(2) * t129 + t161 - t187, -t141, t177 - t191, t149 * qJD(1) + t43 * t120 + t42 * t122 - t187, -t141 * t121, t141 * t119, -t120 * t163 + t186 * t191, -t40 * t122 + t151 * t120 + (-t120 * t48 + t150 * t122) * qJD(1) - t187, 0, 0, 0, 0, 0, -t202 * qJD(5) - t56 * qJDD(5) + t122 * t36 + t62 * t185, t203 * qJD(5) + t57 * qJDD(5) - t122 * t35 + t63 * t185, 0, 0, 0, 0, 0 -(-t122 * t126 + t124 * t57) * t33 + t56 * t14 - (t157 * t124 + t145 * t126) * t189 + t202 * t45 (-t122 * t124 - t126 * t57) * t33 + t56 * t13 - (-t145 * t124 + t157 * t126) * t189 + t202 * t47; 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3) + g(3), 0, 0, 0, t119 * t22 + t121 * t21 + g(3), 0, 0, 0, 0, 0, -t143, t144, 0, 0, 0, 0, 0, t154 - t214, t131 - t155; 0, 0, 0, 0, 0, 0, 0, 0, 0, t174, -t175, -t186 * t129, -t150 * qJD(1) + t138 + t216 + t84, 0, 0, 0, 0, 0, -0.2e1 * t63 * qJD(5) - t152 (t62 - t80) * qJD(5) + t135, 0, 0, 0, 0, 0, t142 + t205, t204 + t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63 * t62, -t62 ^ 2 + t63 ^ 2 (-t62 - t80) * qJD(5) + t135, t152, qJDD(5), t41 * t63 + t134 - t148, -t159 * t102 - t41 * t62 - t147 - t207, t218 * t47 + t195 (t13 - t219) * t126 + (-t189 * t47 - t14) * t124, t204 - t220, t142 - t205, t189 * t63, -pkin(5) * t14 - t10 * t45 + t136 * t124 + t217 * t126 + t3 * t63, -pkin(5) * t13 - t10 * t47 - t217 * t124 + t136 * t126 - t4 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47 * t45, -t45 ^ 2 + t47 ^ 2, t13 + t219, -t165 + (-qJD(6) + t189) * t47, -t33, -g(1) * t25 + t132 * t124 + t156 * t126 + t189 * t4 - t7 * t47 + t5, g(1) * t26 + t3 * t189 + t7 * t45 + (-t156 - t6) * t124 + t132 * t126;];
tau_reg  = t1;
