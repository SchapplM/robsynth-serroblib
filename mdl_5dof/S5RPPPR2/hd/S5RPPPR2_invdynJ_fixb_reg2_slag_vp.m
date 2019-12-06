% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
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
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPPR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:31:40
% EndTime: 2019-12-05 17:31:47
% DurationCPUTime: 2.94s
% Computational Cost: add. (2911->337), mult. (7297->507), div. (0->0), fcn. (5685->10), ass. (0->184)
t138 = sin(pkin(9));
t141 = cos(pkin(9));
t143 = cos(pkin(7));
t194 = qJD(1) * qJD(4);
t139 = sin(pkin(8));
t142 = cos(pkin(8));
t189 = t143 * qJDD(1);
t175 = t142 * t189;
t195 = qJD(1) * qJD(2);
t180 = t143 * t195;
t140 = sin(pkin(7));
t235 = -t143 * pkin(2) - pkin(1);
t159 = t140 * qJ(3) - t235;
t198 = qJD(3) * t140;
t66 = -qJD(1) * t198 - qJDD(1) * t159 + qJDD(2);
t33 = qJ(2) * t175 + t139 * t66 + t142 * t180;
t26 = (-qJ(4) * qJDD(1) - t194) * t143 + t33;
t163 = pkin(3) * t139 - qJ(4) * t142;
t190 = t140 * qJDD(1);
t99 = qJ(2) * t190 + t140 * t195 + qJDD(3);
t37 = (qJDD(1) * t163 - t142 * t194) * t140 + t99;
t12 = t138 * t37 + t141 * t26;
t178 = t139 * t190;
t10 = pkin(6) * t178 + t12;
t176 = t139 * t189;
t32 = -qJ(2) * t176 - t139 * t180 + t142 * t66;
t29 = pkin(3) * t189 + qJDD(4) - t32;
t213 = t143 * t141;
t216 = t140 * t142;
t89 = t138 * t216 + t213;
t75 = qJDD(1) * t89;
t115 = t138 * t189;
t192 = qJDD(1) * t142;
t177 = t141 * t192;
t76 = t140 * t177 - t115;
t14 = t75 * pkin(4) - t76 * pkin(6) + t29;
t144 = sin(qJ(5));
t146 = cos(qJ(5));
t200 = qJD(1) * t143;
t186 = qJ(2) * t200;
t88 = -qJD(1) * t159 + qJD(2);
t49 = -t139 * t186 + t142 * t88;
t39 = pkin(3) * t200 + qJD(4) - t49;
t77 = t89 * qJD(1);
t202 = qJD(1) * t140;
t181 = t141 * t202;
t182 = t138 * t200;
t80 = t142 * t181 - t182;
t15 = t77 * pkin(4) - t80 * pkin(6) + t39;
t185 = t139 * t202;
t50 = t139 * t88 + t142 * t186;
t40 = -qJ(4) * t200 + t50;
t114 = qJ(2) * t202 + qJD(3);
t62 = t163 * t202 + t114;
t19 = t138 * t62 + t141 * t40;
t17 = pkin(6) * t185 + t19;
t6 = t144 * t15 + t146 * t17;
t2 = -qJD(5) * t6 - t144 * t10 + t146 * t14;
t238 = qJD(5) + t77;
t241 = t238 * t6 + t2;
t161 = t144 * t17 - t146 * t15;
t1 = -t161 * qJD(5) + t146 * t10 + t144 * t14;
t240 = t161 * t238 + t1;
t203 = qJD(1) * t139;
t183 = t146 * t203;
t168 = t140 * t183;
t197 = qJD(5) * t144;
t20 = -qJD(5) * t168 - t144 * t178 - t146 * t76 + t197 * t80;
t41 = t144 * t80 - t168;
t239 = t238 * t41 - t20;
t147 = cos(qJ(1));
t214 = t140 * t147;
t208 = t147 * t142;
t145 = sin(qJ(1));
t211 = t145 * t139;
t98 = t143 * t208 + t211;
t56 = t138 * t214 + t98 * t141;
t209 = t147 * t139;
t210 = t145 * t142;
t97 = t143 * t209 - t210;
t237 = t56 * t144 - t97 * t146;
t236 = t97 * t144 + t56 * t146;
t196 = qJ(2) * qJDD(1);
t234 = g(2) * t145;
t233 = g(2) * t147;
t232 = g(3) * t147;
t184 = t144 * t203;
t44 = t140 * t184 + t146 * t80;
t230 = t44 * t41;
t82 = (t138 * t140 + t142 * t213) * qJD(1);
t217 = t139 * t146;
t96 = t141 * t217 - t144 * t142;
t229 = t96 * qJD(5) + t143 * t183 - t144 * t82;
t218 = t139 * t144;
t158 = t141 * t218 + t146 * t142;
t228 = t158 * qJD(5) + t143 * t184 + t146 * t82;
t223 = qJ(2) * t143;
t69 = -t139 * t159 + t142 * t223;
t59 = -t143 * qJ(4) + t69;
t73 = (qJ(2) + t163) * t140;
t28 = t138 * t73 + t141 * t59;
t71 = qJDD(5) + t75;
t227 = t144 * t71;
t226 = t146 * t71;
t222 = qJDD(1) * pkin(1);
t135 = t140 ^ 2;
t148 = qJD(1) ^ 2;
t221 = t135 * t148;
t220 = t138 * t139;
t219 = t139 * t140;
t215 = t140 * t145;
t212 = t143 * t148;
t206 = g(1) * t143 + g(2) * t215;
t134 = t139 ^ 2;
t136 = t142 ^ 2;
t205 = -t134 - t136;
t137 = t143 ^ 2;
t204 = t135 + t137;
t201 = qJD(1) * t142;
t199 = qJD(2) * t143;
t193 = qJDD(1) * t139;
t191 = t135 * qJDD(1);
t188 = t134 * t221;
t187 = t140 * t217;
t179 = t134 * t190;
t174 = -t98 * pkin(3) - t97 * qJ(4);
t173 = t238 ^ 2;
t172 = t144 * t76 - t146 * t178;
t171 = t204 * t148;
t68 = -t139 * t223 - t142 * t159;
t170 = 0.2e1 * t140 * t189;
t169 = 0.2e1 * t204;
t95 = -t143 * t210 + t209;
t54 = t95 * t138 + t141 * t215;
t57 = -t98 * t138 + t141 * t214;
t167 = g(2) * t57 + g(3) * t54;
t94 = t143 * t211 + t208;
t166 = -g(2) * t97 - g(3) * t94;
t132 = t147 * qJ(2);
t165 = t95 * pkin(3) - t94 * qJ(4) + t132;
t164 = g(3) * t145 + t233;
t61 = t143 * pkin(3) - t68;
t92 = -t139 * t198 + t142 * t199;
t11 = -t138 * t26 + t141 * t37;
t18 = -t138 * t40 + t141 * t62;
t27 = -t138 * t59 + t141 * t73;
t90 = -t143 * t138 + t141 * t216;
t22 = t89 * pkin(4) - t90 * pkin(6) + t61;
t24 = pkin(6) * t219 + t28;
t7 = -t144 * t24 + t146 * t22;
t8 = t144 * t22 + t146 * t24;
t160 = (-qJD(4) * t142 + qJD(2)) * t140;
t52 = t140 * t218 + t90 * t146;
t91 = t139 * t199 + t142 * t198;
t157 = qJD(1) * t91 - qJDD(1) * t68 - t32;
t156 = qJD(1) * t92 + qJDD(1) * t69 + t33;
t155 = -t164 - t222;
t154 = -g(1) * t219 + g(2) * t94 - g(3) * t97;
t129 = qJDD(2) - t222;
t153 = -t129 - t155;
t152 = t169 * t195 + t234;
t151 = t99 * t140 + (t195 + t196) * t135;
t150 = (g(2) * qJ(2) + g(3) * t159) * t145 + t159 * t233;
t130 = t137 * qJDD(1);
t116 = t134 * t191;
t79 = t142 * t182 - t181;
t70 = -t143 * qJD(4) + t92;
t64 = t96 * t202;
t63 = t158 * t202;
t55 = -t138 * t215 + t95 * t141;
t51 = t90 * t144 - t187;
t47 = t52 * qJD(5);
t46 = -qJD(5) * t187 + t197 * t90;
t36 = t138 * t160 + t141 * t70;
t35 = t138 * t70 - t141 * t160;
t31 = -t94 * t144 + t55 * t146;
t30 = -t55 * t144 - t94 * t146;
t23 = -pkin(4) * t219 - t27;
t21 = qJD(5) * t44 + t172;
t16 = -pkin(4) * t185 - t18;
t9 = -pkin(4) * t178 - t11;
t4 = -qJD(5) * t8 - t144 * t36 + t146 * t91;
t3 = qJD(5) * t7 + t144 * t91 + t146 * t36;
t5 = [0, 0, 0, 0, 0, qJDD(1), t164, t232 - t234, 0, 0, t191, t170, 0, t130, 0, 0, t153 * t143, -t153 * t140, t169 * t196 + t152 - t232, -g(3) * t132 + (-t129 + t164) * pkin(1) + (t204 * t196 + t152) * qJ(2), t136 * t191, -0.2e1 * t142 * t139 * t191, -0.2e1 * t140 * t175, t116, t139 * t170, t130, g(2) * t98 - g(3) * t95 + t139 * t151 + t143 * t157, t142 * t151 + t143 * t156 + t166, (-t139 * t156 + t142 * t157 + t164) * t140, t33 * t69 + t50 * t92 + t32 * t68 - t49 * t91 - g(2) * (-t145 * qJ(2) + t235 * t147) - g(3) * (t235 * t145 + t132) + (t99 * qJ(2) + qJ(3) * t164 + t114 * qJD(2)) * t140, t76 * t90, -t90 * t75 - t76 * t89, (qJDD(1) * t90 + t76) * t219, t75 * t89, -0.2e1 * t75 * t219, t116, g(2) * t56 - g(3) * t55 + t29 * t89 + t61 * t75 + t91 * t77 + (-qJD(1) * t35 + qJDD(1) * t27 + t11) * t219, t29 * t90 + t61 * t76 + t91 * t80 + (-qJD(1) * t36 - qJDD(1) * t28 - t12) * t219 + t167, -t11 * t90 - t12 * t89 - t27 * t76 - t28 * t75 + t35 * t80 - t36 * t77 - t166, -g(2) * t174 - g(3) * t165 + t11 * t27 + t12 * t28 - t18 * t35 + t19 * t36 + t29 * t61 + t39 * t91 + t150, -t20 * t52 - t44 * t46, t20 * t51 - t52 * t21 + t46 * t41 - t44 * t47, -t20 * t89 - t238 * t46 + t52 * t71, t21 * t51 + t41 * t47, -t21 * t89 - t238 * t47 - t51 * t71, t71 * t89, g(2) * t236 - g(3) * t31 + t16 * t47 + t2 * t89 + t23 * t21 + t238 * t4 + t35 * t41 + t9 * t51 + t7 * t71, -g(2) * t237 - g(3) * t30 - t1 * t89 - t16 * t46 - t23 * t20 - t238 * t3 + t35 * t44 + t9 * t52 - t8 * t71, -t1 * t51 - t161 * t46 - t2 * t52 + t7 * t20 - t8 * t21 - t3 * t41 - t4 * t44 - t6 * t47 - t167, t1 * t8 + t6 * t3 + t2 * t7 - t161 * t4 + t9 * t23 + t16 * t35 - g(2) * (-pkin(4) * t56 + t57 * pkin(6) + t174) - g(3) * (t55 * pkin(4) + t54 * pkin(6) + t165) + t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t189, t190, -t171, -qJ(2) * t171 + qJDD(2) + t155, 0, 0, 0, 0, 0, 0, -t139 * t171 - t175, -t142 * t171 + t176, t205 * t190, t33 * t139 + t32 * t142 + (-t114 * t140 + (t139 * t49 - t142 * t50) * t143) * qJD(1) - t164, 0, 0, 0, 0, 0, 0, -t138 * t179 - t142 * t75 + (t140 * t79 - t143 * t77) * t203, -t141 * t179 - t142 * t76 + (t140 * t82 - t143 * t80) * t203, t82 * t77 - t79 * t80 + (t138 * t76 - t141 * t75) * t139, -t29 * t142 + t18 * t79 - t19 * t82 + (-t11 * t138 + t12 * t141 - t200 * t39) * t139 - t164, 0, 0, 0, 0, 0, 0, -t158 * t71 + t21 * t220 - t229 * t238 - t79 * t41, -t20 * t220 + t228 * t238 - t79 * t44 - t96 * t71, -t158 * t20 - t96 * t21 + t228 * t41 + t229 * t44, t1 * t96 - t158 * t2 - t16 * t79 + t161 * t229 + t9 * t220 - t228 * t6 - t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t142 * t212 + t193) * t140, (t139 * t212 + t192) * t140, t205 * t221, (-t232 + (t139 * t50 + t142 * t49) * qJD(1)) * t140 + t99 + t206, 0, 0, 0, 0, 0, 0, -t138 * t188 + (t141 * t193 - t201 * t77) * t140, -t141 * t188 + (-t138 * t193 - t201 * t80) * t140, -t138 * t75 - t141 * t76 + (t138 * t80 - t141 * t77) * t185, t11 * t141 + t12 * t138 + (-t232 + (-t142 * t39 + (-t138 * t18 + t141 * t19) * t139) * qJD(1)) * t140 + t206, 0, 0, 0, 0, 0, 0, -t141 * t21 - t63 * t238 + (-qJD(5) * t146 * t238 + t185 * t41 - t227) * t138, t141 * t20 - t64 * t238 + (t185 * t44 + t197 * t238 - t226) * t138, -t64 * t41 + t63 * t44 + (-t144 * t20 - t146 * t21 + (t144 * t41 + t146 * t44) * qJD(5)) * t138, -g(3) * t214 - t9 * t141 + t161 * t63 + t6 * t64 + (t16 * t185 + t1 * t146 - t144 * t2 + (-t144 * t6 + t146 * t161) * qJD(5)) * t138 + t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t185 * t80 + t75, -t115 + (-t203 * t77 + t177) * t140, -t77 ^ 2 - t80 ^ 2, t18 * t80 + t19 * t77 + t154 + t29, 0, 0, 0, 0, 0, 0, -t144 * t173 - t80 * t41 + t226, -t146 * t173 - t80 * t44 - t227, -t239 * t146 + (t238 * t44 - t21) * t144, t240 * t144 + t241 * t146 - t16 * t80 + t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t230, -t41 ^ 2 + t44 ^ 2, t239, -t230, -t172 + (-qJD(5) + t238) * t44, t71, g(1) * t51 - g(2) * t30 + g(3) * t237 - t16 * t44 + t241, g(1) * t52 + g(2) * t31 + g(3) * t236 + t16 * t41 - t240, 0, 0;];
tau_reg = t5;
