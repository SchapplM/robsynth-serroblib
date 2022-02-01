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
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:30:19
% EndTime: 2022-01-23 09:30:24
% DurationCPUTime: 1.91s
% Computational Cost: add. (2619->305), mult. (5621->356), div. (0->0), fcn. (3733->12), ass. (0->174)
t138 = sin(qJ(4));
t139 = sin(qJ(3));
t141 = cos(qJ(3));
t219 = cos(qJ(4));
t84 = t138 * t141 + t139 * t219;
t77 = t84 * qJD(1);
t203 = t77 * qJ(5);
t136 = sin(pkin(8));
t111 = pkin(1) * t136 + pkin(6);
t98 = t111 * qJD(1);
t173 = pkin(7) * qJD(1) + t98;
t189 = t139 * qJD(2);
t58 = t141 * t173 + t189;
t52 = t138 * t58;
t124 = t141 * qJD(2);
t57 = -t139 * t173 + t124;
t55 = qJD(3) * pkin(3) + t57;
t24 = t219 * t55 - t52;
t15 = t24 - t203;
t96 = t111 * qJDD(1);
t225 = qJD(2) * qJD(3) + t96;
t131 = qJD(3) + qJD(4);
t201 = pkin(1) * qJDD(1);
t129 = qJDD(3) + qJDD(4);
t176 = qJD(4) * t219;
t169 = pkin(3) * t176;
t217 = pkin(3) * t138;
t224 = -t129 * t217 - t131 * t169;
t137 = cos(pkin(8));
t112 = -pkin(1) * t137 - pkin(2);
t127 = t141 * pkin(3);
t223 = t112 - t127;
t119 = t129 * pkin(4);
t196 = t138 * t139;
t159 = t131 * t196;
t178 = t219 * t141;
t165 = qJD(1) * t178;
t170 = qJDD(1) * t219;
t185 = t141 * qJDD(1);
t168 = -t131 * t165 - t138 * t185 - t139 * t170;
t35 = qJD(1) * t159 + t168;
t205 = t35 * qJ(5);
t222 = t119 + t205;
t132 = qJ(1) + pkin(8);
t121 = cos(t132);
t135 = qJ(3) + qJ(4);
t125 = sin(t135);
t198 = t121 * t125;
t120 = sin(t132);
t200 = t120 * t125;
t126 = cos(t135);
t216 = g(3) * t126;
t221 = g(1) * t198 + g(2) * t200 - t216;
t220 = t77 ^ 2;
t143 = -pkin(7) - pkin(6);
t140 = sin(qJ(1));
t218 = pkin(1) * t140;
t215 = g(3) * t141;
t192 = qJD(1) * t139;
t177 = t138 * t192;
t75 = -t165 + t177;
t214 = t77 * t75;
t213 = pkin(7) + t111;
t11 = pkin(4) * t131 + t15;
t212 = t11 - t15;
t186 = t139 * qJDD(1);
t158 = t138 * t186 - t141 * t170;
t50 = t131 * t84;
t36 = qJD(1) * t50 + t158;
t211 = -t169 * t75 - t217 * t36;
t49 = -qJD(3) * t178 - t141 * t176 + t159;
t210 = -t36 * t84 + t49 * t75;
t29 = t219 * t57 - t52;
t80 = t213 * t139;
t81 = t213 * t141;
t45 = -t138 * t80 + t219 * t81;
t209 = qJ(5) * t36;
t208 = t131 * t75;
t207 = t139 * t98;
t206 = t141 * t98;
t204 = t75 * qJ(5);
t199 = t120 * t126;
t197 = t121 * t126;
t174 = pkin(4) * t75 + qJD(5);
t79 = t223 * qJD(1);
t46 = t174 + t79;
t195 = qJD(5) + t46;
t194 = pkin(4) * t126 + t127;
t133 = t139 ^ 2;
t134 = t141 ^ 2;
t193 = t133 - t134;
t99 = qJD(1) * t112;
t191 = qJD(3) * t139;
t190 = qJD(4) * t138;
t188 = qJD(1) * qJD(3);
t97 = qJDD(1) * t112;
t184 = t219 * pkin(3);
t183 = pkin(3) * t192;
t182 = pkin(3) * t191;
t181 = pkin(3) * t190;
t180 = t139 * qJDD(2) + t141 * t225;
t54 = t219 * t58;
t145 = qJD(1) ^ 2;
t179 = t139 * t145 * t141;
t175 = t139 * t188;
t123 = t141 * qJDD(2);
t27 = qJDD(3) * pkin(3) + t123 + (-pkin(7) * qJDD(1) - t96) * t139 - t58 * qJD(3);
t41 = -t191 * t98 + t180;
t31 = (-t175 + t185) * pkin(7) + t41;
t172 = -t138 * t31 + t219 * t27;
t28 = -t138 * t57 - t54;
t44 = -t138 * t81 - t219 * t80;
t171 = qJD(3) * t213;
t5 = t138 * t27 + t176 * t55 - t190 * t58 + t219 * t31;
t167 = -g(1) * t200 + g(2) * t198;
t166 = g(1) * t199 - g(2) * t197;
t164 = t141 * t175;
t163 = g(1) * t121 + g(2) * t120;
t162 = g(1) * t120 - g(2) * t121;
t142 = cos(qJ(1));
t161 = g(1) * t140 - g(2) * t142;
t83 = -t178 + t196;
t160 = -t35 * t83 + t50 * t77;
t32 = t129 * t84 - t131 * t49;
t67 = t189 + t206;
t25 = t138 * t55 + t54;
t156 = g(1) * t197 + g(2) * t199 + g(3) * t125 - t5;
t72 = t139 * t171;
t73 = t141 * t171;
t12 = -t138 * t73 - t176 * t80 - t190 * t81 - t219 * t72;
t155 = -qJD(1) * t99 + t163;
t154 = t75 * t79 + t156;
t60 = pkin(3) * t175 + qJDD(1) * t223;
t153 = 0.2e1 * qJD(3) * t99 - qJDD(3) * t111;
t152 = -t131 * t177 - t168;
t144 = qJD(3) ^ 2;
t151 = -t111 * t144 + t162 - 0.2e1 * t97;
t6 = -qJD(4) * t25 + t172;
t13 = -qJD(4) * t45 + t138 * t72 - t219 * t73;
t42 = -qJD(3) * t67 - t139 * t96 + t123;
t66 = t124 - t207;
t150 = -t42 * t139 + t41 * t141 + (-t139 * t67 - t141 * t66) * qJD(3);
t14 = pkin(4) * t36 + qJDD(5) + t60;
t149 = t195 * t75 + t156 + t209;
t148 = t6 + t221;
t147 = -t79 * t77 + t148;
t130 = -qJ(5) + t143;
t128 = t142 * pkin(1);
t118 = t127 + pkin(2);
t117 = t184 + pkin(4);
t95 = qJDD(3) * t141 - t139 * t144;
t94 = qJDD(3) * t139 + t141 * t144;
t85 = pkin(2) + t194;
t74 = t75 ^ 2;
t59 = pkin(4) * t77 + t183;
t56 = pkin(4) * t83 + t223;
t43 = pkin(4) * t50 + t182;
t39 = -t74 + t220;
t38 = -qJ(5) * t83 + t45;
t37 = -qJ(5) * t84 + t44;
t34 = -t129 * t83 - t131 * t50;
t19 = t152 + t208;
t18 = -t203 + t29;
t17 = t28 + t204;
t16 = t25 - t204;
t10 = t36 * t83 + t50 * t75;
t9 = -t35 * t84 - t49 * t77;
t8 = t49 * qJ(5) - t84 * qJD(5) + t13;
t7 = -qJ(5) * t50 - qJD(5) * t83 + t12;
t4 = t160 + t210;
t3 = -t160 + t210;
t2 = -qJD(5) * t75 - t209 + t5;
t1 = -t77 * qJD(5) + t222 + t6;
t20 = [0, 0, 0, 0, 0, qJDD(1), t161, g(1) * t142 + g(2) * t140, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t137 * t201 + t162, -0.2e1 * t136 * t201 + t163, 0, (t161 + (t136 ^ 2 + t137 ^ 2) * t201) * pkin(1), qJDD(1) * t133 + 0.2e1 * t164, 0.2e1 * t139 * t185 - 0.2e1 * t188 * t193, t94, qJDD(1) * t134 - 0.2e1 * t164, t95, 0, t139 * t153 + t141 * t151, -t139 * t151 + t141 * t153, (t133 + t134) * t96 + t150 - t163, t97 * t112 - g(1) * (-pkin(2) * t120 + pkin(6) * t121 - t218) - g(2) * (pkin(2) * t121 + pkin(6) * t120 + t128) + t150 * t111, t9, t3, t32, t10, t34, 0, t129 * t44 + t13 * t131 + t182 * t75 + t223 * t36 + t50 * t79 + t60 * t83 + t166, -t12 * t131 - t129 * t45 + t182 * t77 - t223 * t35 - t49 * t79 + t60 * t84 + t167, -t12 * t75 - t13 * t77 + t24 * t49 - t25 * t50 + t35 * t44 - t36 * t45 - t5 * t83 - t6 * t84 - t163, t5 * t45 + t25 * t12 + t6 * t44 + t24 * t13 + t60 * t223 + t79 * t182 - g(1) * (-t118 * t120 - t121 * t143 - t218) - g(2) * (t118 * t121 - t120 * t143 + t128), t9, t3, t32, t10, t34, 0, t129 * t37 + t131 * t8 + t14 * t83 + t36 * t56 + t43 * t75 + t46 * t50 + t166, -t129 * t38 - t131 * t7 + t14 * t84 - t35 * t56 + t43 * t77 - t46 * t49 + t167, -t1 * t84 + t11 * t49 - t16 * t50 - t2 * t83 + t35 * t37 - t36 * t38 - t7 * t75 - t77 * t8 - t163, t2 * t38 + t16 * t7 + t1 * t37 + t11 * t8 + t14 * t56 + t46 * t43 - g(1) * (-t120 * t85 - t121 * t130 - t218) - g(2) * (-t120 * t130 + t121 * t85 + t128); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, t95, -t94, 0, t139 * t41 + t141 * t42 - g(3) + (-t139 * t66 + t141 * t67) * qJD(3), 0, 0, 0, 0, 0, 0, t34, -t32, t4, -t24 * t50 - t25 * t49 + t5 * t84 - t6 * t83 - g(3), 0, 0, 0, 0, 0, 0, t34, -t32, t4, -t1 * t83 - t11 * t50 - t16 * t49 + t2 * t84 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t179, t193 * t145, t186, t179, t185, qJDD(3), -t215 + t123 + (t67 - t206) * qJD(3) + (t155 - t225) * t139, g(3) * t139 + (t66 + t207) * qJD(3) + t155 * t141 - t180, 0, 0, t214, t39, t19, -t214, -t158, t129, -t28 * t131 + (t129 * t219 - t131 * t190 - t192 * t75) * pkin(3) + t147, t131 * t29 - t183 * t77 + t154 + t224, t35 * t184 + (-t24 + t29) * t75 + (t25 + t28 + t181) * t77 + t211, -t24 * t28 - t25 * t29 + (t219 * t6 - t215 + t138 * t5 + (-t138 * t24 + t219 * t25) * qJD(4) + (-qJD(1) * t79 + t163) * t139) * pkin(3), t214, t39, t19, -t214, -t158, t129, t117 * t129 - t17 * t131 - t59 * t75 - t195 * t77 + (-t54 + (-pkin(3) * t131 - t55) * t138) * qJD(4) + t172 + t221 + t222, t131 * t18 - t59 * t77 + t149 + t224, t117 * t35 + (-t11 + t18) * t75 + (t16 + t17 + t181) * t77 + t211, t1 * t117 - t16 * t18 - t11 * t17 - t46 * t59 - g(3) * t194 - t163 * (-pkin(3) * t139 - pkin(4) * t125) + (t2 * t138 + (-t11 * t138 + t16 * t219) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t214, t39, t19, -t214, -t158, t129, t25 * t131 + t147, t131 * t24 + t154, 0, 0, t214, t39, t19, -t214, -t158, t129, t205 + t16 * t131 + 0.2e1 * t119 + (-t174 - t46) * t77 + t148, -pkin(4) * t220 + t131 * t15 + t149, t35 * pkin(4) - t212 * t75, t212 * t16 + (t125 * t163 - t46 * t77 + t1 - t216) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77 * t131 + t36, t152 - t208, -t74 - t220, t11 * t77 + t16 * t75 + t14 - t162;];
tau_reg = t20;
