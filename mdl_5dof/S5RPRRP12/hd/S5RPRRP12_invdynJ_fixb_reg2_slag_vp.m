% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRRP12
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP12_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP12_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP12_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:57:26
% EndTime: 2019-12-31 18:57:30
% DurationCPUTime: 2.77s
% Computational Cost: add. (2691->381), mult. (5258->471), div. (0->0), fcn. (3100->6), ass. (0->209)
t117 = sin(qJ(4));
t120 = cos(qJ(4));
t192 = t120 * qJD(3);
t121 = cos(qJ(3));
t201 = qJD(1) * t121;
t75 = t117 * t201 - t192;
t118 = sin(qJ(3));
t202 = qJD(1) * t118;
t93 = qJD(4) + t202;
t241 = t75 * t93;
t194 = qJD(4) * t121;
t171 = t117 * t194;
t173 = t118 * t192;
t188 = t121 * qJDD(1);
t33 = (t171 + t173) * qJD(1) - qJD(4) * t192 - t117 * qJDD(3) - t120 * t188;
t257 = -t33 - t241;
t114 = t118 ^ 2;
t115 = t121 ^ 2;
t204 = t114 + t115;
t123 = -pkin(1) - pkin(6);
t89 = t123 * qJDD(1) + qJDD(2);
t164 = t204 * t89;
t200 = qJD(3) * t117;
t77 = t120 * t201 + t200;
t239 = t77 * t93;
t199 = qJD(3) * t118;
t174 = t117 * t199;
t222 = qJD(4) * t77;
t34 = -qJD(1) * t174 - t120 * qJDD(3) + t117 * t188 + t222;
t256 = t34 + t239;
t122 = cos(qJ(1));
t111 = g(2) * t122;
t119 = sin(qJ(1));
t112 = g(1) * t119;
t206 = t112 - t111;
t255 = -pkin(4) * t117 + t123;
t221 = t117 * t121;
t212 = t120 * t122;
t216 = t119 * t117;
t59 = -t118 * t216 + t212;
t215 = t119 * t120;
t218 = t118 * t122;
t61 = t117 * t218 + t215;
t254 = -g(1) * t59 - g(2) * t61 + g(3) * t221;
t253 = t77 ^ 2;
t252 = 2 * qJ(2);
t191 = qJD(1) * qJD(3);
t168 = t121 * t191;
t189 = t118 * qJDD(1);
t72 = qJDD(4) + t168 + t189;
t251 = pkin(4) * t72;
t250 = pkin(4) * t75;
t246 = g(3) * t118;
t245 = g(3) * t121;
t152 = pkin(3) * t118 - pkin(7) * t121;
t82 = qJ(2) + t152;
t58 = t82 * qJD(1);
t91 = t123 * qJD(1) + qJD(2);
t81 = t118 * t91;
t65 = qJD(3) * pkin(7) + t81;
t32 = t117 * t58 + t120 * t65;
t25 = -qJ(5) * t75 + t32;
t244 = t25 * t93;
t31 = -t117 * t65 + t120 * t58;
t243 = t31 * t93;
t242 = t32 * t93;
t240 = t77 * t75;
t116 = -qJ(5) - pkin(7);
t24 = -qJ(5) * t77 + t31;
t19 = pkin(4) * t93 + t24;
t238 = -t24 + t19;
t165 = qJD(4) * t116;
t175 = t117 * t202;
t193 = qJD(5) * t120;
t213 = t120 * t121;
t153 = pkin(3) * t121 + pkin(7) * t118;
t80 = t153 * qJD(1);
t43 = t117 * t80 + t91 * t213;
t237 = -qJ(5) * t175 + t117 * t165 + t193 - t43;
t219 = t118 * t120;
t42 = t120 * t80 - t91 * t221;
t236 = -qJD(5) * t117 + t120 * t165 - (pkin(4) * t121 + qJ(5) * t219) * qJD(1) - t42;
t211 = t121 * t122;
t185 = g(2) * t211;
t235 = g(3) * t219 + t120 * t185;
t217 = t118 * t123;
t49 = t117 * t82 + t120 * t217;
t214 = t119 * t121;
t234 = g(1) * t211 + g(2) * t214;
t233 = qJ(5) * t33;
t232 = qJ(5) * t34;
t231 = t117 * t33;
t230 = t117 * t93;
t229 = t120 * t34;
t228 = t120 * t77;
t227 = t121 * t89;
t226 = pkin(1) * qJDD(1);
t225 = qJD(1) * t93;
t224 = qJD(3) * t75;
t223 = qJD(3) * t77;
t220 = t117 * t122;
t124 = qJD(3) ^ 2;
t210 = t123 * t124;
t125 = qJD(1) ^ 2;
t209 = t125 * qJ(2);
t184 = 0.2e1 * qJD(1) * qJD(2);
t208 = (qJDD(1) * qJ(2) + t184) * qJ(2);
t207 = t122 * pkin(1) + t119 * qJ(2);
t205 = t114 - t115;
t203 = -t124 - t125;
t198 = qJD(3) * t121;
t197 = qJD(3) * t123;
t196 = qJD(4) * t117;
t195 = qJD(4) * t120;
t190 = qJDD(3) * t118;
t172 = t121 * t197;
t73 = t153 * qJD(3) + qJD(2);
t187 = t117 * t73 + t120 * t172 + t82 * t195;
t186 = g(1) * t214;
t182 = t93 * t200;
t181 = t93 * t192;
t180 = t117 * t217;
t179 = t121 * t125 * t118;
t167 = -qJD(5) - t250;
t66 = -qJD(3) * pkin(3) - t121 * t91;
t40 = -t167 + t66;
t178 = t40 * t195;
t177 = t93 * t201;
t176 = t122 * pkin(6) + t207;
t170 = t120 * t194;
t169 = -g(2) * t218 + t245;
t166 = -t117 * t123 + pkin(4);
t38 = t73 * qJD(1) + t82 * qJDD(1);
t46 = qJDD(3) * pkin(7) + t118 * t89 + t91 * t198;
t8 = t117 * t38 + t120 * t46 + t58 * t195 - t65 * t196;
t163 = -qJDD(3) * pkin(3) + t91 * t199;
t162 = qJD(4) * t75 - t33;
t161 = -t34 + t222;
t160 = t204 * qJDD(1);
t159 = qJDD(2) - t226;
t45 = t163 - t227;
t158 = -pkin(7) * qJD(4) * t93 - t45;
t157 = g(2) * t176;
t156 = t118 * t168;
t155 = g(1) * t61 - g(2) * t59;
t60 = t118 * t215 + t220;
t62 = t118 * t212 - t216;
t154 = -g(1) * t62 - g(2) * t60;
t151 = g(1) * t122 + g(2) * t119;
t149 = -t206 - t209;
t148 = t117 * t25 + t120 * t19;
t147 = t117 * t19 - t120 * t25;
t146 = t117 * t32 + t120 * t31;
t145 = t117 * t31 - t120 * t32;
t101 = pkin(4) * t120 + pkin(3);
t143 = t101 * t118 + t116 * t121;
t142 = -t185 - t246;
t141 = qJDD(1) * t252 + t184;
t140 = t209 - t89 + t112;
t138 = t117 * t72 + t93 * t195;
t137 = t120 * t72 - t93 * t196;
t136 = -t118 * t112 - t169;
t135 = qJDD(3) * t123 + t191 * t252;
t134 = t34 * pkin(4) + qJDD(5) + t163;
t132 = -pkin(7) * t72 + t93 * t66;
t131 = g(1) * t60 - g(2) * t62 + g(3) * t213 - t8;
t130 = -t49 * qJD(4) + t120 * t73;
t129 = t141 - t151;
t9 = -qJD(4) * t32 - t117 * t46 + t120 * t38;
t128 = -t146 * qJD(4) - t117 * t9 + t120 * t8;
t126 = t9 + t254;
t108 = t122 * qJ(2);
t105 = qJDD(3) * t121;
t87 = t117 * t186;
t85 = t116 * t120;
t84 = t116 * t117;
t74 = t255 * t121;
t71 = t75 ^ 2;
t70 = t120 * t82;
t50 = -pkin(4) * t175 + t81;
t48 = t70 - t180;
t47 = t118 * t197 + (t170 - t174) * pkin(4);
t44 = t118 * t72 + t93 * t198;
t41 = -qJ(5) * t221 + t49;
t37 = -qJ(5) * t213 + t166 * t118 + t70;
t29 = -t71 + t253;
t27 = -t117 * t172 + t130;
t26 = -qJD(4) * t180 + t187;
t23 = t239 - t34;
t22 = -t33 + t241;
t21 = (-t121 * t77 + t93 * t219) * qJD(1) + t138;
t20 = (-t118 * t230 + t121 * t75) * qJD(1) + t137;
t18 = t134 - t227;
t17 = -qJ(5) * t170 + (-qJD(5) * t121 + (qJ(5) * qJD(3) - qJD(4) * t123) * t118) * t117 + t187;
t16 = t75 * t230 - t229;
t15 = t93 * t228 - t231;
t14 = qJ(5) * t173 + (qJ(5) * t196 + t166 * qJD(3) - t193) * t121 + t130;
t13 = t75 * t170 + (t121 * t34 - t75 * t199) * t117;
t12 = -t77 * t171 + (-t121 * t33 - t77 * t199) * t120;
t11 = (-t34 + t182) * t118 + (-t138 - t224) * t121;
t10 = (-t33 - t181) * t118 + (t137 + t223) * t121;
t7 = -t120 * t225 + (-t34 - t182) * t121 + (-t138 + t224) * t118;
t6 = t117 * t225 + (t33 - t181) * t121 + (-t137 + t223) * t118;
t5 = -t256 * t117 + t257 * t120;
t4 = -qJD(5) * t75 - t232 + t8;
t3 = (t117 * t77 + t120 * t75) * t199 + (t231 - t229 + (t117 * t75 - t228) * qJD(4)) * t121;
t2 = -qJD(5) * t77 + t233 + t251 + t9;
t1 = (qJD(1) * t77 + t161 * t118 - t75 * t198) * t120 + (qJD(1) * t75 + t162 * t118 + t77 * t198) * t117;
t28 = [0, 0, 0, 0, 0, qJDD(1), t206, t151, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - t206 - 0.2e1 * t226, t129, -t159 * pkin(1) - g(1) * (-t119 * pkin(1) + t108) - g(2) * t207 + t208, qJDD(1) * t115 - 0.2e1 * t156, -0.2e1 * t118 * t188 + 0.2e1 * t205 * t191, -t124 * t118 + t105, qJDD(1) * t114 + 0.2e1 * t156, -t124 * t121 - t190, 0, t135 * t121 + (t129 - t210) * t118, -t135 * t118 + (t141 - t210) * t121 - t234, -t123 * t160 - t164 + t206, -g(1) * (t123 * t119 + t108) - t157 + t123 * t164 + t208, t12, t3, t10, t13, t11, t44, t27 * t93 + t48 * t72 + (t9 + (-t117 * t66 + t123 * t75) * qJD(3)) * t118 + (qJD(3) * t31 + t45 * t117 - t123 * t34 + t66 * t195) * t121 + t154, -t26 * t93 - t49 * t72 + (-t8 + (-t120 * t66 + t123 * t77) * qJD(3)) * t118 + (-qJD(3) * t32 + t120 * t45 + t123 * t33 - t66 * t196) * t121 + t155, -t26 * t75 - t27 * t77 + t33 * t48 - t34 * t49 + t146 * t199 + (t145 * qJD(4) - t117 * t8 - t120 * t9) * t121 + t234, t8 * t49 + t32 * t26 + t9 * t48 + t31 * t27 - g(1) * (pkin(3) * t218 - pkin(7) * t211 + t108) - t157 + (-t121 * t45 + t66 * t199) * t123 + (-g(1) * t123 - g(2) * t152) * t119, t12, t3, t10, t13, t11, t44, t14 * t93 - t34 * t74 + t37 * t72 + t47 * t75 + (-t40 * t200 + t2) * t118 + (qJD(3) * t19 + t18 * t117 + t178) * t121 + t154, -t17 * t93 + t33 * t74 - t41 * t72 + t47 * t77 + (-t192 * t40 - t4) * t118 + (-qJD(3) * t25 + t120 * t18 - t196 * t40) * t121 + t155, -t14 * t77 - t17 * t75 + t33 * t37 - t34 * t41 + t148 * t199 + (qJD(4) * t147 - t117 * t4 - t120 * t2) * t121 + t234, t4 * t41 + t25 * t17 + t2 * t37 + t19 * t14 - t18 * t74 + t40 * t47 - g(1) * (t101 * t218 + t116 * t211 + t108) - g(2) * (pkin(4) * t220 + t176) + (-g(1) * t255 - g(2) * t143) * t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t125, t149 + t159, 0, 0, 0, 0, 0, 0, t203 * t118 + t105, t203 * t121 - t190, -t160, t164 + t149, 0, 0, 0, 0, 0, 0, t7, t6, t1, -t146 * qJD(1) + (-qJD(3) * t145 - t45) * t121 + (qJD(3) * t66 + t128) * t118 - t206, 0, 0, 0, 0, 0, 0, t7, t6, t1, -t148 * qJD(1) + (-qJD(3) * t147 - t18) * t121 + (qJD(3) * t40 - qJD(4) * t148 - t117 * t2 + t120 * t4) * t118 - t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, -t205 * t125, t188, -t179, -t189, qJDD(3), t246 + (-t140 + t111) * t121, t118 * t140 + t169, 0, 0, t15, t5, t21, t16, t20, -t177, -t31 * t201 - t75 * t81 - pkin(3) * t34 - t42 * t93 + (t158 - t186) * t120 + t132 * t117 + t235, t32 * t201 - t77 * t81 + pkin(3) * t33 + t43 * t93 + t87 + t132 * t120 + (t142 - t158) * t117, t42 * t77 + t43 * t75 + (pkin(7) * t161 - t243 + t8) * t120 + (pkin(7) * t162 - t242 - t9) * t117 + t136, -t66 * t81 - t31 * t42 - t32 * t43 + (-t121 * t206 + t246 - t45) * pkin(3) + (-t118 * t206 + t128 - t245) * pkin(7), t15, t5, t21, t16, t20, -t177, -t19 * t201 - t101 * t34 - t50 * t75 + t72 * t84 + t236 * t93 + (-t18 - t186) * t120 + (t40 * t202 + (t40 + t250) * qJD(4)) * t117 + t235, t178 + t101 * t33 - t50 * t77 + t85 * t72 + t87 - t237 * t93 + (t121 * t25 + t219 * t40) * qJD(1) + (pkin(4) * t222 + t142 + t18) * t117, t33 * t84 + t34 * t85 - t236 * t77 - t237 * t75 + (-t19 * t93 + t4) * t120 + (-t2 - t244) * t117 + t136, -t4 * t85 + t2 * t84 - t18 * t101 + g(3) * t143 + (pkin(4) * t196 - t50) * t40 + t237 * t25 + t236 * t19 - t206 * (t101 * t121 - t116 * t118); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t240, t29, t22, -t240, t23, t72, -t66 * t77 + t126 + t242, t66 * t75 + t131 + t243, 0, 0, t240, t29, t22, -t240, t23, t72, 0.2e1 * t251 + t233 + t244 + (t167 - t40) * t77 + t126, -pkin(4) * t253 + t232 + t24 * t93 + (qJD(5) + t40) * t75 + t131, pkin(4) * t33 - t238 * t75, t238 * t25 + (-t40 * t77 + t2 + t254) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t256, t257, -t71 - t253, -t246 + t19 * t77 + t25 * t75 + (t206 - t89) * t121 + t134;];
tau_reg = t28;
