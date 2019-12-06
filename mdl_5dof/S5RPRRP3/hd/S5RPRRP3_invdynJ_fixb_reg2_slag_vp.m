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
% Datum: 2019-12-05 18:04
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:04:02
% EndTime: 2019-12-05 18:04:06
% DurationCPUTime: 1.97s
% Computational Cost: add. (2619->302), mult. (5621->351), div. (0->0), fcn. (3733->12), ass. (0->171)
t133 = qJ(3) + qJ(4);
t124 = sin(t133);
t125 = cos(t133);
t130 = qJ(1) + pkin(8);
t119 = sin(t130);
t120 = cos(t130);
t161 = g(2) * t119 - g(3) * t120;
t221 = -g(1) * t125 - t124 * t161;
t136 = sin(qJ(4));
t137 = sin(qJ(3));
t139 = cos(qJ(3));
t219 = cos(qJ(4));
t84 = t136 * t139 + t219 * t137;
t77 = t84 * qJD(1);
t200 = t77 * qJ(5);
t134 = sin(pkin(8));
t110 = pkin(1) * t134 + pkin(6);
t96 = t110 * qJD(1);
t171 = pkin(7) * qJD(1) + t96;
t190 = qJD(2) * t137;
t58 = t171 * t139 + t190;
t52 = t136 * t58;
t123 = t139 * qJD(2);
t57 = -t171 * t137 + t123;
t55 = qJD(3) * pkin(3) + t57;
t24 = t219 * t55 - t52;
t15 = t24 - t200;
t94 = t110 * qJDD(1);
t225 = qJD(2) * qJD(3) + t94;
t129 = qJD(3) + qJD(4);
t127 = qJDD(3) + qJDD(4);
t174 = t219 * qJD(4);
t165 = pkin(3) * t174;
t216 = pkin(3) * t136;
t224 = -t127 * t216 - t129 * t165;
t198 = pkin(1) * qJDD(1);
t135 = cos(pkin(8));
t111 = -pkin(1) * t135 - pkin(2);
t126 = t139 * pkin(3);
t223 = t111 - t126;
t118 = t127 * pkin(4);
t195 = t136 * t137;
t158 = t129 * t195;
t176 = t219 * t139;
t164 = qJD(1) * t176;
t167 = qJDD(1) * t219;
t184 = t139 * qJDD(1);
t166 = -t129 * t164 - t136 * t184 - t137 * t167;
t35 = qJD(1) * t158 + t166;
t202 = t35 * qJ(5);
t222 = t118 + t202;
t220 = t77 ^ 2;
t141 = -pkin(7) - pkin(6);
t138 = sin(qJ(1));
t218 = pkin(1) * t138;
t140 = cos(qJ(1));
t217 = pkin(1) * t140;
t214 = g(1) * t139;
t191 = qJD(1) * t137;
t175 = t136 * t191;
t75 = -t164 + t175;
t212 = t77 * t75;
t211 = pkin(7) + t110;
t11 = pkin(4) * t129 + t15;
t210 = t11 - t15;
t185 = t137 * qJDD(1);
t157 = t136 * t185 - t139 * t167;
t50 = t129 * t84;
t36 = qJD(1) * t50 + t157;
t209 = -t75 * t165 - t36 * t216;
t49 = -qJD(3) * t176 - t139 * t174 + t158;
t208 = -t84 * t36 + t49 * t75;
t29 = t219 * t57 - t52;
t80 = t211 * t137;
t81 = t211 * t139;
t45 = -t136 * t80 + t219 * t81;
t196 = t120 * t125;
t197 = t119 * t125;
t207 = g(2) * t196 + g(3) * t197;
t206 = qJ(5) * t36;
t205 = t129 * t75;
t204 = t137 * t96;
t203 = t139 * t96;
t201 = t75 * qJ(5);
t172 = pkin(4) * t75 + qJD(5);
t79 = t223 * qJD(1);
t46 = t172 + t79;
t194 = qJD(5) + t46;
t193 = pkin(4) * t125 + t126;
t131 = t137 ^ 2;
t132 = t139 ^ 2;
t192 = t131 - t132;
t97 = qJD(1) * t111;
t189 = qJD(3) * t137;
t188 = qJD(4) * t136;
t187 = qJD(1) * qJD(3);
t95 = qJDD(1) * t111;
t182 = t219 * pkin(3);
t181 = pkin(3) * t191;
t180 = pkin(3) * t189;
t179 = pkin(3) * t188;
t178 = t137 * qJDD(2) + t225 * t139;
t54 = t219 * t58;
t143 = qJD(1) ^ 2;
t177 = t137 * t143 * t139;
t173 = t137 * t187;
t122 = t139 * qJDD(2);
t27 = qJDD(3) * pkin(3) + t122 + (-pkin(7) * qJDD(1) - t94) * t137 - t58 * qJD(3);
t41 = -t96 * t189 + t178;
t31 = (-t173 + t184) * pkin(7) + t41;
t169 = -t136 * t31 + t219 * t27;
t28 = -t136 * t57 - t54;
t44 = -t136 * t81 - t219 * t80;
t168 = qJD(3) * t211;
t5 = t136 * t27 + t55 * t174 - t58 * t188 + t219 * t31;
t163 = t139 * t173;
t162 = g(2) * t120 + g(3) * t119;
t160 = g(2) * t140 + g(3) * t138;
t83 = -t176 + t195;
t159 = -t35 * t83 + t50 * t77;
t32 = t127 * t84 - t129 * t49;
t67 = t190 + t203;
t25 = t136 * t55 + t54;
t155 = t162 * t124;
t72 = t137 * t168;
t73 = t139 * t168;
t12 = -t136 * t73 - t80 * t174 - t81 * t188 - t219 * t72;
t154 = -qJD(1) * t97 - t161;
t153 = g(1) * t124 - g(2) * t197 + g(3) * t196 - t5;
t60 = pkin(3) * t173 + qJDD(1) * t223;
t152 = 0.2e1 * t97 * qJD(3) - qJDD(3) * t110;
t151 = -t129 * t175 - t166;
t142 = qJD(3) ^ 2;
t150 = -t110 * t142 + t162 - 0.2e1 * t95;
t149 = t75 * t79 + t153;
t6 = -t25 * qJD(4) + t169;
t13 = -t45 * qJD(4) + t136 * t72 - t219 * t73;
t42 = -qJD(3) * t67 - t137 * t94 + t122;
t66 = t123 - t204;
t148 = -t137 * t42 + t139 * t41 + (-t137 * t67 - t139 * t66) * qJD(3);
t14 = pkin(4) * t36 + qJDD(5) + t60;
t147 = t194 * t75 + t153 + t206;
t146 = t6 + t221;
t145 = -t79 * t77 + t146;
t128 = -qJ(5) + t141;
t117 = t126 + pkin(2);
t116 = t182 + pkin(4);
t93 = qJDD(3) * t139 - t137 * t142;
t92 = qJDD(3) * t137 + t139 * t142;
t85 = pkin(2) + t193;
t74 = t75 ^ 2;
t59 = pkin(4) * t77 + t181;
t56 = pkin(4) * t83 + t223;
t43 = pkin(4) * t50 + t180;
t39 = -t74 + t220;
t38 = -qJ(5) * t83 + t45;
t37 = -qJ(5) * t84 + t44;
t34 = -t127 * t83 - t129 * t50;
t19 = t151 + t205;
t18 = -t200 + t29;
t17 = t28 + t201;
t16 = t25 - t201;
t10 = t36 * t83 + t50 * t75;
t9 = -t35 * t84 - t49 * t77;
t8 = t49 * qJ(5) - t84 * qJD(5) + t13;
t7 = -qJ(5) * t50 - qJD(5) * t83 + t12;
t4 = t159 + t208;
t3 = -t159 + t208;
t2 = -qJD(5) * t75 - t206 + t5;
t1 = -t77 * qJD(5) + t222 + t6;
t20 = [0, 0, 0, 0, 0, qJDD(1), t160, -g(2) * t138 + g(3) * t140, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t135 * t198 + t162, -0.2e1 * t134 * t198 - t161, 0, (t160 + (t134 ^ 2 + t135 ^ 2) * t198) * pkin(1), qJDD(1) * t131 + 0.2e1 * t163, 0.2e1 * t137 * t184 - 0.2e1 * t192 * t187, t92, qJDD(1) * t132 - 0.2e1 * t163, t93, 0, t152 * t137 + t150 * t139, -t150 * t137 + t152 * t139, (t131 + t132) * t94 + t148 + t161, t95 * t111 - g(2) * (-pkin(2) * t120 - t119 * pkin(6) - t217) - g(3) * (-pkin(2) * t119 + pkin(6) * t120 - t218) + t148 * t110, t9, t3, t32, t10, t34, 0, t127 * t44 + t129 * t13 + t180 * t75 + t223 * t36 + t50 * t79 + t60 * t83 + t207, -t12 * t129 - t127 * t45 + t180 * t77 - t223 * t35 - t49 * t79 + t60 * t84 - t155, -t12 * t75 - t13 * t77 + t24 * t49 - t25 * t50 + t35 * t44 - t36 * t45 - t5 * t83 - t6 * t84 + t161, t5 * t45 + t25 * t12 + t6 * t44 + t24 * t13 + t60 * t223 + t79 * t180 - g(2) * (-t117 * t120 + t119 * t141 - t217) - g(3) * (-t117 * t119 - t120 * t141 - t218), t9, t3, t32, t10, t34, 0, t127 * t37 + t129 * t8 + t14 * t83 + t36 * t56 + t43 * t75 + t46 * t50 + t207, -t127 * t38 - t129 * t7 + t14 * t84 - t35 * t56 + t43 * t77 - t46 * t49 - t155, -t1 * t84 + t11 * t49 - t16 * t50 - t2 * t83 + t35 * t37 - t36 * t38 - t7 * t75 - t77 * t8 + t161, t2 * t38 + t16 * t7 + t1 * t37 + t11 * t8 + t14 * t56 + t46 * t43 - g(2) * (t119 * t128 - t120 * t85 - t217) - g(3) * (-t119 * t85 - t120 * t128 - t218); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(1), 0, 0, 0, 0, 0, 0, t93, -t92, 0, t137 * t41 + t139 * t42 - g(1) + (-t137 * t66 + t139 * t67) * qJD(3), 0, 0, 0, 0, 0, 0, t34, -t32, t4, -t24 * t50 - t25 * t49 + t5 * t84 - t6 * t83 - g(1), 0, 0, 0, 0, 0, 0, t34, -t32, t4, -t1 * t83 - t11 * t50 - t16 * t49 + t2 * t84 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t177, t192 * t143, t185, t177, t184, qJDD(3), -t214 + t122 + (t67 - t203) * qJD(3) + (t154 - t225) * t137, g(1) * t137 + (t66 + t204) * qJD(3) + t154 * t139 - t178, 0, 0, t212, t39, t19, -t212, -t157, t127, -t28 * t129 + (t219 * t127 - t129 * t188 - t75 * t191) * pkin(3) + t145, t129 * t29 - t181 * t77 + t149 + t224, t35 * t182 + (-t24 + t29) * t75 + (t25 + t28 + t179) * t77 + t209, -t24 * t28 - t25 * t29 + (t219 * t6 - t214 + t136 * t5 + (-t136 * t24 + t219 * t25) * qJD(4) + (-qJD(1) * t79 - t161) * t137) * pkin(3), t212, t39, t19, -t212, -t157, t127, t116 * t127 - t17 * t129 - t59 * t75 - t194 * t77 + (-t54 + (-pkin(3) * t129 - t55) * t136) * qJD(4) + t169 + t221 + t222, t129 * t18 - t59 * t77 + t147 + t224, t116 * t35 + (-t11 + t18) * t75 + (t16 + t17 + t179) * t77 + t209, t1 * t116 - t16 * t18 - t11 * t17 - t46 * t59 - g(1) * t193 + t161 * (-pkin(3) * t137 - pkin(4) * t124) + (t2 * t136 + (-t11 * t136 + t219 * t16) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t212, t39, t19, -t212, -t157, t127, t25 * t129 + t145, t129 * t24 + t149, 0, 0, t212, t39, t19, -t212, -t157, t127, t202 + t16 * t129 + 0.2e1 * t118 + (-t172 - t46) * t77 + t146, -pkin(4) * t220 + t129 * t15 + t147, pkin(4) * t35 - t210 * t75, t210 * t16 + (-t46 * t77 + t1 + t221) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77 * t129 + t36, t151 - t205, -t74 - t220, t11 * t77 + t16 * t75 + t14 - t162;];
tau_reg = t20;
