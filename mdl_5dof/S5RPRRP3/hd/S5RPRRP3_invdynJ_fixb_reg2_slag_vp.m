% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRRP3
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
% Datum: 2020-01-03 11:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:47:50
% EndTime: 2020-01-03 11:47:53
% DurationCPUTime: 1.79s
% Computational Cost: add. (2619->303), mult. (5621->352), div. (0->0), fcn. (3733->12), ass. (0->172)
t138 = sin(qJ(4));
t139 = sin(qJ(3));
t141 = cos(qJ(3));
t219 = cos(qJ(4));
t84 = t138 * t141 + t219 * t139;
t77 = t84 * qJD(1);
t202 = t77 * qJ(5);
t136 = sin(pkin(8));
t110 = pkin(1) * t136 + pkin(6);
t96 = t110 * qJD(1);
t173 = pkin(7) * qJD(1) + t96;
t192 = qJD(2) * t139;
t58 = t173 * t141 + t192;
t52 = t138 * t58;
t123 = t141 * qJD(2);
t57 = -t173 * t139 + t123;
t55 = qJD(3) * pkin(3) + t57;
t24 = t219 * t55 - t52;
t15 = t24 - t202;
t94 = t110 * qJDD(1);
t226 = qJD(2) * qJD(3) + t94;
t132 = qJ(1) + pkin(8);
t119 = sin(t132);
t120 = cos(t132);
t224 = g(2) * t119 - g(3) * t120;
t131 = qJD(3) + qJD(4);
t129 = qJDD(3) + qJDD(4);
t176 = t219 * qJD(4);
t167 = pkin(3) * t176;
t218 = pkin(3) * t138;
t225 = -t129 * t218 - t131 * t167;
t200 = pkin(1) * qJDD(1);
t137 = cos(pkin(8));
t111 = -pkin(1) * t137 - pkin(2);
t127 = t141 * pkin(3);
t223 = t111 - t127;
t118 = t129 * pkin(4);
t197 = t138 * t139;
t160 = t131 * t197;
t178 = t219 * t141;
t166 = qJD(1) * t178;
t169 = qJDD(1) * t219;
t186 = t141 * qJDD(1);
t168 = -t131 * t166 - t138 * t186 - t139 * t169;
t35 = qJD(1) * t160 + t168;
t204 = t35 * qJ(5);
t222 = t118 + t204;
t135 = qJ(3) + qJ(4);
t124 = sin(t135);
t198 = t120 * t124;
t199 = t119 * t124;
t125 = cos(t135);
t217 = g(1) * t125;
t221 = g(2) * t199 - g(3) * t198 - t217;
t220 = t77 ^ 2;
t143 = -pkin(7) - pkin(6);
t216 = g(1) * t141;
t193 = qJD(1) * t139;
t177 = t138 * t193;
t75 = -t166 + t177;
t214 = t77 * t75;
t213 = pkin(7) + t110;
t11 = pkin(4) * t131 + t15;
t212 = t11 - t15;
t187 = t139 * qJDD(1);
t159 = t138 * t187 - t141 * t169;
t50 = t131 * t84;
t36 = qJD(1) * t50 + t159;
t211 = -t75 * t167 - t36 * t218;
t49 = -qJD(3) * t178 - t141 * t176 + t160;
t210 = -t84 * t36 + t49 * t75;
t29 = t219 * t57 - t52;
t80 = t213 * t139;
t81 = t213 * t141;
t45 = -t138 * t80 + t219 * t81;
t209 = g(2) * t198 + g(3) * t199;
t208 = qJ(5) * t36;
t207 = t131 * t75;
t206 = t139 * t96;
t205 = t141 * t96;
t203 = t75 * qJ(5);
t174 = pkin(4) * t75 + qJD(5);
t79 = t223 * qJD(1);
t46 = t174 + t79;
t196 = qJD(5) + t46;
t195 = pkin(4) * t125 + t127;
t133 = t139 ^ 2;
t134 = t141 ^ 2;
t194 = t133 - t134;
t97 = qJD(1) * t111;
t191 = qJD(3) * t139;
t190 = qJD(4) * t138;
t189 = qJD(1) * qJD(3);
t95 = qJDD(1) * t111;
t184 = t219 * pkin(3);
t183 = pkin(3) * t193;
t182 = pkin(3) * t191;
t181 = pkin(3) * t190;
t180 = t139 * qJDD(2) + t226 * t141;
t54 = t219 * t58;
t145 = qJD(1) ^ 2;
t179 = t139 * t145 * t141;
t175 = t139 * t189;
t122 = t141 * qJDD(2);
t27 = qJDD(3) * pkin(3) + t122 + (-pkin(7) * qJDD(1) - t94) * t139 - t58 * qJD(3);
t41 = -t96 * t191 + t180;
t31 = (-t175 + t186) * pkin(7) + t41;
t171 = -t138 * t31 + t219 * t27;
t28 = -t138 * t57 - t54;
t44 = -t138 * t81 - t219 * t80;
t170 = qJD(3) * t213;
t5 = t138 * t27 + t55 * t176 - t58 * t190 + t219 * t31;
t165 = t141 * t175;
t164 = g(2) * t120 + g(3) * t119;
t140 = sin(qJ(1));
t142 = cos(qJ(1));
t162 = -g(2) * t142 - g(3) * t140;
t83 = -t178 + t197;
t161 = -t35 * t83 + t50 * t77;
t32 = t129 * t84 - t131 * t49;
t67 = t192 + t205;
t25 = t138 * t55 + t54;
t157 = t164 * t125;
t72 = t139 * t170;
t73 = t141 * t170;
t12 = -t138 * t73 - t80 * t176 - t81 * t190 - t219 * t72;
t156 = -qJD(1) * t97 + t224;
t155 = g(1) * t124 + t224 * t125 - t5;
t60 = pkin(3) * t175 + qJDD(1) * t223;
t154 = 0.2e1 * t97 * qJD(3) - qJDD(3) * t110;
t153 = -t131 * t177 - t168;
t144 = qJD(3) ^ 2;
t152 = t110 * t144 + t164 + 0.2e1 * t95;
t151 = t75 * t79 + t155;
t6 = -qJD(4) * t25 + t171;
t13 = -t45 * qJD(4) + t138 * t72 - t219 * t73;
t42 = -qJD(3) * t67 - t139 * t94 + t122;
t66 = t123 - t206;
t150 = -t42 * t139 + t41 * t141 + (-t139 * t67 - t141 * t66) * qJD(3);
t14 = pkin(4) * t36 + qJDD(5) + t60;
t149 = t196 * t75 + t155 + t208;
t148 = t6 + t221;
t147 = -t79 * t77 + t148;
t130 = -qJ(5) + t143;
t128 = t142 * pkin(1);
t126 = t140 * pkin(1);
t117 = t127 + pkin(2);
t116 = t184 + pkin(4);
t93 = qJDD(3) * t141 - t139 * t144;
t92 = qJDD(3) * t139 + t141 * t144;
t85 = pkin(2) + t195;
t74 = t75 ^ 2;
t59 = pkin(4) * t77 + t183;
t56 = pkin(4) * t83 + t223;
t43 = pkin(4) * t50 + t182;
t39 = -t74 + t220;
t38 = -qJ(5) * t83 + t45;
t37 = -qJ(5) * t84 + t44;
t34 = -t129 * t83 - t131 * t50;
t19 = t153 + t207;
t18 = -t202 + t29;
t17 = t28 + t203;
t16 = t25 - t203;
t10 = t36 * t83 + t50 * t75;
t9 = -t35 * t84 - t49 * t77;
t8 = t49 * qJ(5) - t84 * qJD(5) + t13;
t7 = -qJ(5) * t50 - qJD(5) * t83 + t12;
t4 = t161 + t210;
t3 = -t161 + t210;
t2 = -qJD(5) * t75 - t208 + t5;
t1 = -t77 * qJD(5) + t222 + t6;
t20 = [0, 0, 0, 0, 0, qJDD(1), t162, g(2) * t140 - g(3) * t142, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t137 * t200 - t164, -0.2e1 * t136 * t200 + t224, 0, (t162 + (t136 ^ 2 + t137 ^ 2) * t200) * pkin(1), qJDD(1) * t133 + 0.2e1 * t165, 0.2e1 * t139 * t186 - 0.2e1 * t194 * t189, t92, qJDD(1) * t134 - 0.2e1 * t165, t93, 0, t139 * t154 - t141 * t152, t139 * t152 + t141 * t154, (t133 + t134) * t94 + t150 - t224, t95 * t111 - g(2) * (pkin(2) * t120 + pkin(6) * t119 + t128) - g(3) * (pkin(2) * t119 - pkin(6) * t120 + t126) + t150 * t110, t9, t3, t32, t10, t34, 0, t129 * t44 + t13 * t131 + t182 * t75 + t223 * t36 + t50 * t79 + t60 * t83 - t157, -t12 * t131 - t129 * t45 + t182 * t77 - t223 * t35 - t49 * t79 + t60 * t84 + t209, -t12 * t75 - t13 * t77 + t24 * t49 - t25 * t50 + t35 * t44 - t36 * t45 - t5 * t83 - t6 * t84 - t224, t5 * t45 + t25 * t12 + t6 * t44 + t24 * t13 + t60 * t223 + t79 * t182 - g(2) * (t117 * t120 - t119 * t143 + t128) - g(3) * (t117 * t119 + t120 * t143 + t126), t9, t3, t32, t10, t34, 0, t129 * t37 + t131 * t8 + t14 * t83 + t36 * t56 + t43 * t75 + t46 * t50 - t157, -t129 * t38 - t131 * t7 + t14 * t84 - t35 * t56 + t43 * t77 - t46 * t49 + t209, -t1 * t84 + t11 * t49 - t16 * t50 - t2 * t83 + t35 * t37 - t36 * t38 - t7 * t75 - t77 * t8 - t224, t2 * t38 + t16 * t7 + t1 * t37 + t11 * t8 + t14 * t56 + t46 * t43 - g(2) * (-t119 * t130 + t120 * t85 + t128) - g(3) * (t119 * t85 + t120 * t130 + t126); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(1), 0, 0, 0, 0, 0, 0, t93, -t92, 0, t139 * t41 + t141 * t42 - g(1) + (-t139 * t66 + t141 * t67) * qJD(3), 0, 0, 0, 0, 0, 0, t34, -t32, t4, -t24 * t50 - t25 * t49 + t5 * t84 - t6 * t83 - g(1), 0, 0, 0, 0, 0, 0, t34, -t32, t4, -t1 * t83 - t11 * t50 - t16 * t49 + t2 * t84 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t179, t194 * t145, t187, t179, t186, qJDD(3), -t216 + t122 + (t67 - t205) * qJD(3) + (t156 - t226) * t139, g(1) * t139 + (t66 + t206) * qJD(3) + t156 * t141 - t180, 0, 0, t214, t39, t19, -t214, -t159, t129, -t28 * t131 + (t219 * t129 - t131 * t190 - t75 * t193) * pkin(3) + t147, t131 * t29 - t183 * t77 + t151 + t225, t35 * t184 + (-t24 + t29) * t75 + (t25 + t28 + t181) * t77 + t211, -t24 * t28 - t25 * t29 + (t219 * t6 - t216 + t138 * t5 + (-t138 * t24 + t219 * t25) * qJD(4) + (-qJD(1) * t79 + t224) * t139) * pkin(3), t214, t39, t19, -t214, -t159, t129, t116 * t129 - t17 * t131 - t59 * t75 - t196 * t77 + (-t54 + (-pkin(3) * t131 - t55) * t138) * qJD(4) + t171 + t221 + t222, t131 * t18 - t59 * t77 + t149 + t225, t116 * t35 + (-t11 + t18) * t75 + (t16 + t17 + t181) * t77 + t211, t1 * t116 - t16 * t18 - t11 * t17 - t46 * t59 - g(1) * t195 - t224 * (-pkin(3) * t139 - pkin(4) * t124) + (t2 * t138 + (-t11 * t138 + t219 * t16) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t214, t39, t19, -t214, -t159, t129, t25 * t131 + t147, t131 * t24 + t151, 0, 0, t214, t39, t19, -t214, -t159, t129, t204 + t16 * t131 + 0.2e1 * t118 + (-t174 - t46) * t77 + t148, -pkin(4) * t220 + t131 * t15 + t149, pkin(4) * t35 - t212 * t75, t212 * t16 + (t124 * t224 - t46 * t77 + t1 - t217) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77 * t131 + t36, t153 - t207, -t74 - t220, t11 * t77 + t16 * t75 + t14 + t164;];
tau_reg = t20;
