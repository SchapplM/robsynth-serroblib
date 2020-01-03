% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRRP6
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
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:43:21
% EndTime: 2019-12-31 18:43:26
% DurationCPUTime: 2.49s
% Computational Cost: add. (2905->372), mult. (6124->469), div. (0->0), fcn. (3809->10), ass. (0->195)
t132 = sin(pkin(8));
t114 = t132 * pkin(1) + pkin(6);
t98 = t114 * qJDD(1);
t254 = -qJD(2) * qJD(3) - t98;
t129 = qJ(1) + pkin(8);
t120 = sin(t129);
t121 = cos(t129);
t167 = g(1) * t121 + g(2) * t120;
t136 = sin(qJ(3));
t203 = qJD(4) * t136;
t253 = qJD(1) * t203 - qJDD(3);
t135 = sin(qJ(4));
t138 = cos(qJ(4));
t194 = t136 * qJDD(1);
t139 = cos(qJ(3));
t197 = t139 * qJD(1);
t47 = (qJD(3) * (qJD(4) + t197) + t194) * t135 + t253 * t138;
t107 = -qJD(4) + t197;
t199 = t138 * qJD(3);
t207 = qJD(1) * t136;
t87 = t135 * t207 - t199;
t221 = t87 * t107;
t186 = t139 * t199;
t46 = -qJD(1) * t186 - qJD(4) * t199 + t253 * t135 - t138 * t194;
t252 = -t46 + t221;
t201 = t135 * qJD(3);
t89 = t138 * t207 + t201;
t220 = t89 * t107;
t251 = t47 - t220;
t219 = pkin(1) * qJDD(1);
t156 = t167 * t136;
t119 = t138 * pkin(4) + pkin(3);
t233 = qJ(5) + pkin(7);
t161 = t139 * t119 + t136 * t233;
t214 = t135 * t136;
t213 = t135 * t139;
t62 = t120 * t213 + t121 * t138;
t64 = t120 * t138 - t121 * t213;
t250 = -g(1) * t64 + g(2) * t62 + g(3) * t214;
t211 = t136 * t138;
t126 = t139 * qJDD(1);
t196 = qJD(1) * qJD(3);
t181 = t136 * t196;
t83 = qJDD(4) - t126 + t181;
t249 = -t107 * (-t135 * t203 + t186) + t83 * t211;
t247 = t89 ^ 2;
t244 = t83 * pkin(4);
t243 = t87 * pkin(4);
t242 = pkin(4) * t135;
t241 = g(1) * t120;
t238 = g(2) * t121;
t237 = g(3) * t136;
t236 = g(3) * t139;
t133 = cos(pkin(8));
t235 = t133 * pkin(1);
t234 = t89 * t87;
t100 = t114 * qJD(1);
t200 = t136 * qJD(2);
t68 = t139 * t100 + t200;
t59 = qJD(3) * pkin(7) + t68;
t169 = t139 * pkin(3) + t136 * pkin(7);
t159 = -pkin(2) - t169;
t79 = t159 - t235;
t60 = t79 * qJD(1);
t29 = -t135 * t59 + t138 * t60;
t20 = -t89 * qJ(5) + t29;
t19 = -t107 * pkin(4) + t20;
t232 = -t20 + t19;
t231 = -t87 * t186 - t47 * t211;
t86 = t136 * t100;
t67 = t139 * qJD(2) - t86;
t168 = pkin(3) * t136 - pkin(7) * t139;
t92 = t168 * qJD(1);
t44 = t135 * t92 + t138 * t67;
t202 = qJD(4) * t138;
t93 = t168 * qJD(3);
t230 = t135 * t93 + t79 * t202;
t177 = qJD(4) * t233;
t198 = t138 * qJD(5);
t229 = t198 - t44 + (qJ(5) * t197 - t177) * t135;
t210 = t138 * t139;
t158 = pkin(4) * t136 - qJ(5) * t210;
t43 = -t135 * t67 + t138 * t92;
t228 = -t158 * qJD(1) - t135 * qJD(5) - t138 * t177 - t43;
t227 = t136 * t114 * t201 + t138 * t93;
t91 = t114 * t210;
t50 = t135 * t79 + t91;
t226 = t167 * t211;
t225 = t29 * t107;
t30 = t135 * t60 + t138 * t59;
t224 = t30 * t107;
t223 = t46 * qJ(5);
t222 = t47 * qJ(5);
t218 = qJD(3) * t87;
t217 = qJD(4) * t87;
t216 = qJD(4) * t89;
t215 = t121 * t135;
t130 = t136 ^ 2;
t131 = t139 ^ 2;
t208 = t130 - t131;
t115 = -pkin(2) - t235;
t101 = qJD(1) * t115;
t206 = qJD(3) * t136;
t205 = qJD(3) * t139;
t204 = qJD(4) * t135;
t99 = qJDD(1) * t115;
t192 = -t136 * qJDD(2) + t254 * t139;
t191 = t89 * t205;
t142 = qJD(1) ^ 2;
t190 = t136 * t142 * t139;
t179 = -qJD(5) - t243;
t58 = -qJD(3) * pkin(3) - t67;
t45 = -t179 + t58;
t189 = t45 * t202;
t140 = cos(qJ(1));
t188 = t140 * pkin(1) + t121 * pkin(2) + t120 * pkin(6);
t187 = t107 * t201;
t185 = t136 * t202;
t184 = t107 * t207;
t39 = t139 * qJDD(2) - t100 * t205 + t254 * t136;
t35 = -qJDD(3) * pkin(3) - t39;
t16 = t47 * pkin(4) + qJDD(5) + t35;
t183 = -t16 - t236;
t137 = sin(qJ(1));
t180 = -t137 * pkin(1) + t121 * pkin(6);
t178 = t46 * t139 + t89 * t206;
t38 = -t100 * t206 - t192;
t34 = qJDD(3) * pkin(7) + t38;
t48 = qJD(1) * t93 + qJDD(1) * t79;
t6 = t135 * t48 + t138 * t34 + t60 * t202 - t59 * t204;
t176 = -t46 + t217;
t174 = t89 * t185;
t173 = t139 * t181;
t172 = -pkin(7) * qJD(4) * t107 + t35;
t171 = -g(1) * t62 - g(2) * t64;
t63 = -t120 * t210 + t215;
t65 = t120 * t135 + t121 * t210;
t170 = -g(1) * t63 - g(2) * t65;
t166 = g(1) * t137 - g(2) * t140;
t21 = -t87 * qJ(5) + t30;
t165 = t135 * t21 + t138 * t19;
t164 = t135 * t19 - t138 * t21;
t163 = -t135 * t30 - t138 * t29;
t162 = t135 * t29 - t138 * t30;
t155 = t107 * t202 - t135 * t83;
t154 = -qJD(1) * t101 + t167;
t152 = -pkin(7) * t83 - t107 * t58;
t141 = qJD(3) ^ 2;
t151 = t114 * t141 + t238 + 0.2e1 * t99;
t150 = g(1) * t65 - g(2) * t63 + g(3) * t211 - t6;
t149 = 0.2e1 * t101 * qJD(3) - qJDD(3) * t114;
t148 = -t139 * t167 - t237;
t7 = -qJD(4) * t30 - t135 * t34 + t138 * t48;
t1 = -t89 * qJD(5) + t223 + t244 + t7;
t2 = -t87 * qJD(5) - t222 + t6;
t147 = -qJD(4) * t165 - t1 * t135 + t2 * t138;
t146 = qJD(4) * t163 - t7 * t135 + t6 * t138;
t145 = -t39 * t136 + t38 * t139 + (-t136 * t68 - t139 * t67) * qJD(3);
t144 = t7 + t250;
t110 = g(3) * t213;
t104 = t136 * t241;
t103 = t233 * t138;
t102 = t233 * t135;
t97 = qJDD(3) * t139 - t141 * t136;
t96 = qJDD(3) * t136 + t141 * t139;
t82 = t87 ^ 2;
t74 = (t114 + t242) * t136;
t70 = t138 * t79;
t53 = t200 + (qJD(1) * t242 + t100) * t139;
t52 = t114 * t205 + (t139 * t201 + t185) * pkin(4);
t51 = -t107 * t206 - t83 * t139;
t49 = -t114 * t213 + t70;
t42 = -qJ(5) * t214 + t50;
t37 = -t82 + t247;
t36 = -qJ(5) * t211 + t70 + (-t114 * t135 - pkin(4)) * t139;
t27 = -t220 - t47;
t26 = -t46 - t221;
t25 = (t107 * t210 - t136 * t89) * qJD(1) - t155;
t24 = t107 * t204 + t138 * t83 + (-t107 * t213 + t136 * t87) * qJD(1);
t23 = -t50 * qJD(4) + t227;
t22 = (-t136 * t199 - t139 * t204) * t114 + t230;
t18 = -t135 * t221 - t47 * t138;
t17 = -t46 * t135 - t138 * t220;
t15 = t87 * t185 + (t136 * t47 + t87 * t205) * t135;
t14 = t89 * t186 + (-t46 * t138 - t89 * t204) * t136;
t13 = (-qJ(5) * qJD(4) - qJD(3) * t114) * t211 + (-qJD(5) * t136 + (-qJ(5) * qJD(3) - qJD(4) * t114) * t139) * t135 + t230;
t12 = -t136 * t198 + t158 * qJD(3) + (-t91 + (qJ(5) * t136 - t79) * t135) * qJD(4) + t227;
t11 = (t47 + t187) * t139 + (t155 - t218) * t136;
t10 = (-t47 + t187) * t139 + (t155 + t218) * t136;
t9 = t178 + t249;
t8 = t178 - t249;
t5 = -t251 * t135 + t252 * t138;
t4 = -t174 + (-t191 + (t46 + t217) * t136) * t135 + t231;
t3 = t174 + (t136 * t176 + t191) * t135 + t231;
t28 = [0, 0, 0, 0, 0, qJDD(1), t166, g(1) * t140 + g(2) * t137, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t133 * t219 - t238 + t241, -0.2e1 * t132 * t219 + t167, 0, (t166 + (t132 ^ 2 + t133 ^ 2) * t219) * pkin(1), t130 * qJDD(1) + 0.2e1 * t173, 0.2e1 * t136 * t126 - 0.2e1 * t208 * t196, t96, t131 * qJDD(1) - 0.2e1 * t173, t97, 0, t149 * t136 + (-t151 + t241) * t139, t136 * t151 + t139 * t149 - t104, (t130 + t131) * t98 + t145 - t167, t99 * t115 - g(1) * (-t120 * pkin(2) + t180) - g(2) * t188 + t145 * t114, t14, t4, t9, t15, t11, t51, -t23 * t107 + t49 * t83 + (-t7 + (t114 * t87 + t135 * t58) * qJD(3)) * t139 + (qJD(3) * t29 + t114 * t47 + t35 * t135 + t58 * t202) * t136 + t170, t22 * t107 - t50 * t83 + (t6 + (t114 * t89 + t138 * t58) * qJD(3)) * t139 + (-qJD(3) * t30 - t114 * t46 + t35 * t138 - t58 * t204) * t136 + t171, -t22 * t87 - t23 * t89 + t49 * t46 - t50 * t47 + t104 + t163 * t205 + (qJD(4) * t162 - t135 * t6 - t138 * t7 - t238) * t136, t6 * t50 + t30 * t22 + t7 * t49 + t29 * t23 - g(1) * t180 - g(2) * (t121 * t169 + t188) - t159 * t241 + (t35 * t136 + t58 * t205) * t114, t14, t4, t9, t15, t11, t51, -t12 * t107 + t36 * t83 + t74 * t47 + t52 * t87 + (t45 * t201 - t1) * t139 + (qJD(3) * t19 + t16 * t135 + t189) * t136 + t170, t13 * t107 - t42 * t83 - t74 * t46 + t52 * t89 + (t45 * t199 + t2) * t139 + (-qJD(3) * t21 + t16 * t138 - t45 * t204) * t136 + t171, -t12 * t89 - t13 * t87 + t36 * t46 - t42 * t47 + t104 - t165 * t205 + (qJD(4) * t164 - t1 * t138 - t135 * t2 - t238) * t136, t2 * t42 + t21 * t13 + t1 * t36 + t19 * t12 + t16 * t74 + t45 * t52 - g(1) * (pkin(4) * t215 + t180) - g(2) * (t161 * t121 + t188) + (-g(1) * (-pkin(2) - t161) - g(2) * t242) * t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, t97, -t96, 0, t38 * t136 + t39 * t139 - g(3) + (-t136 * t67 + t139 * t68) * qJD(3), 0, 0, 0, 0, 0, 0, t10, t8, t3, -g(3) + (-qJD(3) * t162 - t35) * t139 + (qJD(3) * t58 + t146) * t136, 0, 0, 0, 0, 0, 0, t10, t8, t3, -g(3) + (-qJD(3) * t164 - t16) * t139 + (qJD(3) * t45 + t147) * t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t190, t208 * t142, t194, t190, t126, qJDD(3), t68 * qJD(3) + t136 * t154 - t236 + t39, t237 + (t67 + t86) * qJD(3) + t154 * t139 + t192, 0, 0, t17, t5, t25, t18, t24, t184, -t29 * t207 - pkin(3) * t47 + t43 * t107 - t68 * t87 + (-t172 - t236) * t138 + t152 * t135 + t226, t30 * t207 + pkin(3) * t46 - t44 * t107 - t68 * t89 + t110 + t152 * t138 + (-t156 + t172) * t135, t43 * t89 + t44 * t87 + (t6 + t225 + (-t47 + t216) * pkin(7)) * t138 + (pkin(7) * t176 + t224 - t7) * t135 + t148, -t29 * t43 - t30 * t44 - t58 * t68 + (-t236 - t35 + t156) * pkin(3) + (t146 + t148) * pkin(7), t17, t5, t25, t18, t24, t184, -t19 * t207 - t102 * t83 - t119 * t47 - t53 * t87 + t183 * t138 - t228 * t107 + (-t45 * t197 + (t45 + t243) * qJD(4)) * t135 + t226, t189 - t103 * t83 + t119 * t46 - t53 * t89 + t110 + t229 * t107 + (t136 * t21 - t45 * t210) * qJD(1) + (pkin(4) * t216 - t156 + t16) * t135, -t237 - t102 * t46 - t103 * t47 - t228 * t89 - t229 * t87 + (qJD(1) * t165 - t167) * t139 + t147, t2 * t103 - t1 * t102 - t16 * t119 - g(3) * t161 + (pkin(4) * t204 - t53) * t45 + t229 * t21 + t228 * t19 + t167 * (t119 * t136 - t139 * t233); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t234, t37, t26, -t234, t27, t83, -t58 * t89 + t144 - t224, t58 * t87 + t150 - t225, 0, 0, t234, t37, t26, -t234, t27, t83, 0.2e1 * t244 + t223 - t21 * t107 + (t179 - t45) * t89 + t144, -t247 * pkin(4) + t222 - t20 * t107 + (qJD(5) + t45) * t87 + t150, t46 * pkin(4) - t232 * t87, t232 * t21 + (-t45 * t89 + t1 + t250) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t251, t252, -t82 - t247, t19 * t89 + t21 * t87 - t156 - t183;];
tau_reg = t28;
