% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRPRR6
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRR6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR6_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:58:00
% EndTime: 2019-12-05 15:58:05
% DurationCPUTime: 2.76s
% Computational Cost: add. (3669->366), mult. (8816->513), div. (0->0), fcn. (7227->14), ass. (0->192)
t135 = cos(pkin(10));
t249 = cos(qJ(4));
t193 = t249 * t135;
t132 = sin(pkin(10));
t139 = sin(qJ(4));
t218 = t139 * t132;
t156 = t193 - t218;
t134 = sin(pkin(5));
t142 = cos(qJ(2));
t221 = t134 * t142;
t147 = t156 * t221;
t240 = pkin(7) + qJ(3);
t104 = t240 * t132;
t105 = t240 * t135;
t157 = -t249 * t104 - t139 * t105;
t238 = -qJD(1) * t147 + t156 * qJD(3) + qJD(4) * t157;
t140 = sin(qJ(2));
t213 = qJD(1) * t134;
t191 = t140 * t213;
t188 = qJD(4) * t249;
t210 = qJD(4) * t139;
t94 = t132 * t210 - t135 * t188;
t99 = t249 * t132 + t139 * t135;
t95 = t99 * qJD(4);
t261 = t95 * pkin(4) + t94 * pkin(8) - t191;
t138 = sin(qJ(5));
t141 = cos(qJ(5));
t189 = qJD(2) * t218;
t117 = qJD(2) * t193;
t184 = qJDD(2) * t249;
t202 = t135 * qJDD(2);
t194 = qJD(4) * t117 + t132 * t184 + t139 * t202;
t52 = qJD(4) * t189 - t194;
t203 = t132 * qJDD(2);
t170 = -t135 * t184 + t139 * t203;
t53 = qJD(2) * t95 + t170;
t121 = t135 * pkin(3) + pkin(2);
t206 = qJDD(1) * t134;
t186 = t142 * t206;
t211 = qJD(2) * t140;
t187 = qJD(1) * t211;
t256 = t134 * t187 + qJDD(3);
t158 = -t186 + t256;
t61 = -t121 * qJDD(2) + t158;
t16 = t53 * pkin(4) + t52 * pkin(8) + t61;
t102 = qJD(2) * qJ(3) + t191;
t136 = cos(pkin(5));
t212 = qJD(1) * t136;
t113 = t135 * t212;
t234 = pkin(7) * qJD(2);
t58 = t113 + (-t102 - t234) * t132;
t72 = t135 * t102 + t132 * t212;
t59 = t135 * t234 + t72;
t27 = t139 * t58 + t249 * t59;
t23 = qJD(4) * pkin(8) + t27;
t190 = t142 * t213;
t171 = qJD(3) - t190;
t81 = -t121 * qJD(2) + t171;
t90 = -t117 + t189;
t92 = t99 * qJD(2);
t28 = t90 * pkin(4) - t92 * pkin(8) + t81;
t167 = t138 * t23 - t141 * t28;
t205 = qJDD(1) * t136;
t111 = t135 * t205;
t204 = qJDD(2) * qJ(3);
t70 = t140 * t206 + t204 + (qJD(3) + t190) * qJD(2);
t36 = t111 + (-pkin(7) * qJDD(2) - t70) * t132;
t45 = t132 * t205 + t135 * t70;
t37 = pkin(7) * t202 + t45;
t201 = -t139 * t36 - t58 * t188 - t249 * t37;
t5 = -t59 * t210 - t201;
t3 = qJDD(4) * pkin(8) + t5;
t1 = -t167 * qJD(5) + t138 * t16 + t141 * t3;
t82 = qJD(5) + t90;
t173 = t167 * t82 + t1;
t12 = t138 * t28 + t141 * t23;
t2 = -qJD(5) * t12 - t138 * t3 + t141 * t16;
t260 = t12 * t82 + t2;
t179 = t138 * t82;
t69 = t138 * qJD(4) + t141 * t92;
t259 = t69 * t179;
t131 = pkin(10) + qJ(4);
t123 = sin(t131);
t229 = cos(pkin(9));
t181 = t229 * t142;
t133 = sin(pkin(9));
t224 = t133 * t140;
t86 = -t136 * t181 + t224;
t182 = t229 * t140;
t223 = t133 * t142;
t88 = t136 * t223 + t182;
t177 = g(1) * t88 + g(2) * t86;
t257 = g(3) * t221 - t177;
t258 = t257 * t123;
t143 = qJD(2) ^ 2;
t153 = (qJDD(2) * t142 - t140 * t143) * t134;
t228 = qJD(5) * t69;
t21 = -t141 * qJDD(4) - t138 * t52 + t228;
t185 = t139 * t37 - t249 * t36;
t6 = -t27 * qJD(4) - t185;
t124 = cos(t131);
t232 = t139 * t59;
t26 = t249 * t58 - t232;
t22 = -qJD(4) * pkin(4) - t26;
t255 = t22 * qJD(5) * t99 + t124 * t177;
t227 = qJDD(2) * pkin(2);
t75 = t158 - t227;
t166 = t177 - t75;
t254 = (-g(3) * t142 + t187) * t134 + t166 + t227;
t168 = (-t132 * t102 + t113) * t132 - t72 * t135;
t253 = t168 * t142 - (-qJD(2) * pkin(2) + t171) * t140;
t252 = t92 ^ 2;
t50 = -pkin(4) * t156 - t99 * pkin(8) - t121;
t64 = -t139 * t104 + t249 * t105;
t18 = -t138 * t64 + t141 * t50;
t251 = t18 * qJD(5) + t261 * t138 + t238 * t141;
t19 = t138 * t50 + t141 * t64;
t250 = -t19 * qJD(5) - t238 * t138 + t261 * t141;
t248 = g(3) * t134;
t207 = t141 * qJD(4);
t67 = t138 * t92 - t207;
t245 = t67 * t90;
t244 = t69 * t67;
t243 = t69 * t92;
t242 = t92 * t67;
t241 = t92 * t90;
t208 = qJD(5) * t141;
t239 = -t138 * t21 - t67 * t208;
t148 = t99 * t221;
t237 = -qJD(1) * t148 + qJD(3) * t99 + qJD(4) * t64;
t87 = t136 * t182 + t223;
t236 = -t86 * t121 + t240 * t87;
t89 = -t136 * t224 + t181;
t235 = -t88 * t121 + t240 * t89;
t47 = qJDD(5) + t53;
t233 = t138 * t47;
t209 = qJD(5) * t138;
t20 = -qJD(5) * t207 - t138 * qJDD(4) + t141 * t52 + t92 * t209;
t231 = t20 * t138;
t230 = t21 * t141;
t225 = t133 * t134;
t222 = t134 * t140;
t220 = t240 * t140;
t219 = t138 * t142;
t217 = t141 * t142;
t216 = t142 * t143;
t215 = qJDD(1) - g(3);
t128 = t132 ^ 2;
t130 = t135 ^ 2;
t214 = t128 + t130;
t200 = g(3) * t222;
t198 = t99 * t209;
t197 = t99 * t208;
t196 = t134 * t219;
t195 = t134 * t217;
t192 = t134 * t211;
t183 = t134 * t229;
t178 = t141 * t82;
t176 = g(1) * t89 + g(2) * t87;
t175 = pkin(4) * t124 + pkin(8) * t123;
t174 = t47 * t99 - t82 * t94;
t169 = t12 * t138 - t141 * t167;
t164 = t141 * t47 - t90 * t179 - t82 * t209;
t84 = -t132 * t222 + t136 * t135;
t85 = t136 * t132 + t135 * t222;
t162 = -t139 * t85 + t249 * t84;
t41 = t139 * t84 + t249 * t85;
t29 = -t138 * t41 - t195;
t161 = -t141 * t41 + t196;
t160 = -pkin(8) * t47 + t82 * t22;
t155 = -g(1) * (-t89 * t123 + t124 * t225) - g(2) * (-t87 * t123 - t124 * t183) - g(3) * (-t123 * t222 + t136 * t124);
t55 = -t123 * t183 + t87 * t124;
t57 = t123 * t225 + t89 * t124;
t77 = t136 * t123 + t124 * t222;
t154 = -g(1) * t57 - g(2) * t55 - g(3) * t77;
t4 = -qJDD(4) * pkin(4) - t6;
t152 = t155 - t4;
t149 = -t22 * t94 + t4 * t99 - t176;
t146 = -t257 + t186;
t44 = -t132 * t70 + t111;
t145 = -t44 * t132 + t45 * t135 - t176;
t144 = pkin(8) * qJD(5) * t82 - t152;
t101 = t121 * t221;
t83 = t90 ^ 2;
t48 = t92 * pkin(4) + t90 * pkin(8);
t25 = qJD(2) * t148 + t41 * qJD(4);
t24 = qJD(2) * t147 + t162 * qJD(4);
t14 = t138 * t48 + t141 * t26;
t13 = -t138 * t26 + t141 * t48;
t10 = t161 * qJD(5) - t138 * t24 + t141 * t192;
t9 = t29 * qJD(5) + t138 * t192 + t141 * t24;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t215, 0, 0, 0, 0, 0, 0, t153, (-qJDD(2) * t140 - t216) * t134, 0, -g(3) + (t136 ^ 2 + (t140 ^ 2 + t142 ^ 2) * t134 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t135 * t153, -t132 * t153, t214 * t134 * t216 + (-t132 * t84 + t135 * t85) * qJDD(2), t44 * t84 + t45 * t85 - g(3) + (-qJD(2) * t253 - t142 * t75) * t134, 0, 0, 0, 0, 0, 0, -t25 * qJD(4) + t162 * qJDD(4) + (-t142 * t53 + t90 * t211) * t134, -t24 * qJD(4) - t41 * qJDD(4) + (t142 * t52 + t92 * t211) * t134, t162 * t52 - t24 * t90 + t25 * t92 - t41 * t53, t27 * t24 - t26 * t25 + t6 * t162 + t5 * t41 - g(3) + (-t142 * t61 + t81 * t211) * t134, 0, 0, 0, 0, 0, 0, t10 * t82 - t162 * t21 + t25 * t67 + t29 * t47, t161 * t47 + t162 * t20 + t25 * t69 - t9 * t82, -t10 * t69 + t161 * t21 + t29 * t20 - t9 * t67, -t1 * t161 - t10 * t167 + t12 * t9 - t162 * t4 + t2 * t29 + t22 * t25 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t146, -t215 * t222 + t176, 0, 0, t128 * qJDD(2), 0.2e1 * t132 * t202, 0, t130 * qJDD(2), 0, 0, t254 * t135, -t254 * t132, -t200 + t145 + (t171 * qJD(2) + t204) * t214, -t168 * qJD(3) + t166 * pkin(2) + t145 * qJ(3) + (-g(3) * (pkin(2) * t142 + qJ(3) * t140) + t253 * qJD(1)) * t134, -t52 * t99 - t92 * t94, -t156 * t52 - t99 * t53 + t94 * t90 - t92 * t95, -t94 * qJD(4) + t99 * qJDD(4), -t156 * t53 + t90 * t95, -t95 * qJD(4) + qJDD(4) * t156, 0, -t237 * qJD(4) + qJDD(4) * t157 - t121 * t53 - t124 * t257 - t156 * t61 - t90 * t191 + t81 * t95, -t238 * qJD(4) - t64 * qJDD(4) + t121 * t52 - t92 * t191 + t61 * t99 - t81 * t94 + t258, t156 * t5 + t157 * t52 + t237 * t92 - t238 * t90 + t26 * t94 - t27 * t95 - t64 * t53 - t6 * t99 - t176 - t200, t5 * t64 + t6 * t157 - t61 * t121 - t81 * t191 - g(1) * t235 - g(2) * t236 - g(3) * (t134 * t220 + t101) + t238 * t27 - t237 * t26, -t69 * t198 + (-t20 * t99 - t69 * t94) * t141, (t138 * t69 + t141 * t67) * t94 + (t231 - t230 + (t138 * t67 - t141 * t69) * qJD(5)) * t99, t141 * t174 + t156 * t20 - t198 * t82 + t69 * t95, t67 * t197 + (t21 * t99 - t67 * t94) * t138, -t138 * t174 + t156 * t21 - t197 * t82 - t67 * t95, -t156 * t47 + t82 * t95, -t167 * t95 + t18 * t47 - t2 * t156 - t157 * t21 + t250 * t82 + t237 * t67 + t255 * t141 + t149 * t138 - (t124 * t217 + t138 * t140) * t248, t1 * t156 - t12 * t95 - t19 * t47 + t157 * t20 - t251 * t82 + t237 * t69 - t255 * t138 + t149 * t141 - (-t124 * t219 + t140 * t141) * t248, t18 * t20 - t19 * t21 + t169 * t94 - t250 * t69 - t251 * t67 - t258 + (-t1 * t138 - t2 * t141 + (-t12 * t141 - t138 * t167) * qJD(5)) * t99, t1 * t19 + t2 * t18 - t4 * t157 - g(1) * (-t175 * t88 + t235) - g(2) * (-t175 * t86 + t236) - g(3) * t101 + t237 * t22 - (t142 * t175 + t220) * t248 + t251 * t12 - t250 * t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t202, t203, -t214 * t143, t168 * qJD(2) - t146 - t227 + t256, 0, 0, 0, 0, 0, 0, 0.2e1 * t92 * qJD(4) + t170, (-t90 - t189) * qJD(4) + t194, -t83 - t252, t26 * t92 + t27 * t90 + t257 + t61, 0, 0, 0, 0, 0, 0, t164 - t242, -t141 * t82 ^ 2 - t233 - t243, (t20 - t245) * t141 + t259 + t239, t173 * t138 + t260 * t141 - t22 * t92 + t257; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t241, -t83 + t252, (t90 - t189) * qJD(4) + t194, -t241, -t170, qJDD(4), -t81 * t92 + t155 - t185, t81 * t90 + (t26 + t232) * qJD(4) - t154 + t201, 0, 0, t178 * t69 - t231, (-t20 - t245) * t141 - t259 + t239, t178 * t82 + t233 - t243, t179 * t67 - t230, t164 + t242, -t82 * t92, -pkin(4) * t21 - t13 * t82 + t138 * t160 - t141 * t144 + t167 * t92 - t27 * t67, pkin(4) * t20 + t12 * t92 + t138 * t144 + t14 * t82 + t141 * t160 - t27 * t69, t13 * t69 + t14 * t67 + ((-t21 + t228) * pkin(8) + t173) * t141 + ((qJD(5) * t67 - t20) * pkin(8) - t260) * t138 + t154, t167 * t13 - t12 * t14 - t22 * t27 + t152 * pkin(4) + (-qJD(5) * t169 + t1 * t141 - t2 * t138 + t154) * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t244, -t67 ^ 2 + t69 ^ 2, t67 * t82 - t20, -t244, t69 * t82 - t21, t47, -t22 * t69 - g(1) * (-t57 * t138 + t88 * t141) - g(2) * (-t55 * t138 + t86 * t141) - g(3) * (-t77 * t138 - t195) + t260, t22 * t67 - g(1) * (-t88 * t138 - t57 * t141) - g(2) * (-t86 * t138 - t55 * t141) - g(3) * (-t77 * t141 + t196) - t173, 0, 0;];
tau_reg = t7;
