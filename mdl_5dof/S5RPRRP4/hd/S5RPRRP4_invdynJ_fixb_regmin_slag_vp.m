% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRP4
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
% tau_reg [5x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:32:52
% EndTime: 2022-01-23 09:32:59
% DurationCPUTime: 2.90s
% Computational Cost: add. (2804->328), mult. (6757->426), div. (0->0), fcn. (4829->10), ass. (0->196)
t136 = sin(pkin(8));
t138 = sin(qJ(4));
t142 = cos(qJ(3));
t197 = qJDD(1) * t142;
t139 = sin(qJ(3));
t198 = qJDD(1) * t139;
t141 = cos(qJ(4));
t255 = t136 * t141;
t195 = qJD(3) + qJD(4);
t87 = t138 * t142 + t139 * t141;
t258 = t195 * t87;
t22 = (qJD(1) * t258 + t138 * t198) * t136 - t197 * t255;
t131 = t136 ^ 2;
t250 = 0.2e1 * t131;
t140 = sin(qJ(1));
t143 = cos(qJ(1));
t183 = g(1) * t140 - g(2) * t143;
t137 = cos(pkin(8));
t196 = t137 * qJDD(1);
t111 = -qJDD(3) + t196;
t100 = -qJDD(4) + t111;
t204 = qJD(4) * t141;
t209 = qJD(1) * t137;
t113 = -qJD(3) + t209;
t106 = -qJD(4) + t113;
t246 = pkin(3) * t106;
t257 = t138 * pkin(3) * t100 + t204 * t246;
t222 = t136 * t142;
t194 = pkin(7) * t222;
t229 = qJ(2) * t139;
t156 = -t137 * t229 - t194;
t93 = -pkin(2) * t137 - pkin(6) * t136 - pkin(1);
t73 = t93 * qJD(1) + qJD(2);
t65 = t142 * t73;
t42 = t156 * qJD(1) + t65;
t28 = -pkin(3) * t113 + t42;
t228 = qJ(2) * t142;
t112 = t137 * t228;
t210 = qJD(1) * t136;
t191 = t139 * t210;
t43 = -pkin(7) * t191 + qJD(1) * t112 + t139 * t73;
t35 = t138 * t43;
t181 = t141 * t28 - t35;
t169 = t138 * t191;
t208 = qJD(1) * t142;
t190 = t136 * t208;
t60 = t141 * t190 - t169;
t54 = t60 * qJ(5);
t8 = t181 - t54;
t236 = qJ(5) * t22;
t94 = t100 * pkin(4);
t256 = t236 - t94;
t201 = qJD(1) * qJD(2);
t202 = qJ(2) * qJDD(1);
t160 = t201 + t202;
t254 = t160 * t137;
t226 = qJDD(1) * pkin(1);
t253 = t226 + t183;
t135 = qJ(3) + qJ(4);
t123 = sin(t135);
t243 = g(3) * t136;
t124 = cos(t135);
t221 = t137 * t140;
t66 = t123 * t221 + t124 * t143;
t220 = t137 * t143;
t68 = -t123 * t220 + t124 * t140;
t252 = -g(1) * t68 + g(2) * t66 + t123 * t243;
t251 = t60 ^ 2;
t132 = t137 ^ 2;
t5 = -pkin(4) * t106 + t8;
t249 = -t8 + t5;
t245 = pkin(7) * t136;
t157 = qJD(1) * t87;
t57 = t136 * t157;
t242 = t60 * t57;
t241 = t141 * t42 - t35;
t83 = t142 * t93;
t46 = -t194 + t83 + (-pkin(3) - t229) * t137;
t223 = t136 * t139;
t230 = t139 * t93 + t112;
t53 = -pkin(7) * t223 + t230;
t240 = t138 * t46 + t141 * t53;
t219 = t138 * t139;
t86 = t141 * t142 - t219;
t239 = (t195 - t209) * t86;
t238 = t137 * t157 - t258;
t237 = t195 * t169;
t173 = t195 * t142;
t185 = t138 * t197;
t23 = (t185 + (qJD(1) * t173 + t198) * t141) * t136 - t237;
t235 = qJ(5) * t23;
t234 = qJ(5) * t57;
t233 = t106 * t57;
t232 = t106 * t60;
t37 = t141 * t43;
t206 = qJD(3) * t142;
t207 = qJD(2) * t137;
t231 = t142 * t207 + t93 * t206;
t227 = qJD(3) * t73;
t225 = (-qJ(5) - pkin(7) - pkin(6)) * t136;
t144 = qJD(1) ^ 2;
t224 = t131 * t144;
t218 = t139 * t140;
t217 = t139 * t143;
t216 = t140 * t142;
t215 = t142 * t143;
t184 = -pkin(4) * t57 - qJD(5);
t74 = pkin(3) * t191 + qJ(2) * t210;
t40 = -t184 + t74;
t214 = qJD(5) + t40;
t92 = t142 * pkin(3) + pkin(4) * t124;
t109 = t136 * pkin(3) * t206;
t80 = t136 * qJD(2) + t109;
t114 = pkin(3) * t223;
t85 = t136 * qJ(2) + t114;
t213 = t143 * pkin(1) + t140 * qJ(2);
t212 = t131 + t132;
t134 = t142 ^ 2;
t211 = t139 ^ 2 - t134;
t205 = qJD(4) * t138;
t203 = qJD(3) + t113;
t200 = qJD(1) * qJD(3);
t199 = qJDD(1) * t136;
t192 = qJ(2) * qJD(3) * t137;
t189 = qJ(2) * t196;
t188 = t142 * t201;
t187 = t142 * t200;
t168 = qJD(1) * t192;
t72 = t93 * qJDD(1) + qJDD(2);
t64 = t142 * t72;
t15 = -pkin(3) * t111 + t64 + (-pkin(7) * t199 - t168) * t142 + (-t189 - t227 + (qJD(3) * t245 - t207) * qJD(1)) * t139;
t175 = t137 * t188 + t139 * t72 + t142 * t189 + t73 * t206;
t153 = -t139 * t168 + t175;
t21 = (-t187 - t198) * t245 + t153;
t182 = -t138 * t21 + t141 * t15;
t38 = t156 * qJD(3) + t231;
t39 = -t139 * t207 + (-t112 + (-t93 + t245) * t139) * qJD(3);
t180 = -t138 * t38 + t141 * t39;
t179 = -t138 * t42 - t37;
t178 = -t138 * t53 + t141 * t46;
t177 = -t138 * t15 - t141 * t21 - t28 * t204 + t43 * t205;
t176 = t212 * t144;
t174 = qJD(1) * t203;
t172 = pkin(3) * t190;
t49 = qJ(2) * t199 + qJD(1) * t109 + qJDD(1) * t114 + t136 * t201;
t171 = t111 + t196;
t170 = 0.2e1 * t212;
t167 = -g(1) * t66 - g(2) * t68;
t67 = t123 * t143 - t124 * t221;
t69 = t123 * t140 + t124 * t220;
t166 = -g(1) * t67 - g(2) * t69;
t165 = g(1) * t143 + g(2) * t140;
t163 = -t138 * t28 - t37;
t161 = qJD(3) * (t113 + t209);
t159 = t138 * t39 + t141 * t38 + t46 * t204 - t53 * t205;
t158 = t170 * t201;
t12 = pkin(4) * t23 + qJDD(5) + t49;
t155 = -t113 ^ 2 - t224;
t154 = g(1) * t69 - g(2) * t67 + t124 * t243 + t177;
t152 = t163 * qJD(4) + t182;
t151 = t74 * t57 + t154;
t149 = t214 * t57 + t154 + t235;
t148 = t152 + t252;
t147 = -t74 * t60 + t148;
t126 = t143 * qJ(2);
t120 = qJDD(2) - t226;
t119 = pkin(3) * t141 + pkin(4);
t91 = pkin(3) * t139 + pkin(4) * t123;
t88 = pkin(2) + t92;
t78 = t137 * t215 + t218;
t77 = -t137 * t217 + t216;
t76 = -t137 * t216 + t217;
t75 = t137 * t218 + t215;
t71 = t86 * t136;
t70 = t87 * t136;
t55 = t57 ^ 2;
t50 = pkin(4) * t60 + t172;
t47 = pkin(4) * t70 + t85;
t34 = -t195 * t136 * t219 + t173 * t255;
t33 = t258 * t136;
t25 = pkin(4) * t34 + t80;
t24 = -t55 + t251;
t20 = -qJ(5) * t70 + t240;
t18 = -pkin(4) * t137 - qJ(5) * t71 + t178;
t17 = -t232 + (-t185 + (-t195 * t208 - t198) * t141) * t136 + t237;
t16 = -t22 - t233;
t11 = -t54 + t241;
t10 = t179 + t234;
t9 = -t163 - t234;
t7 = -t100 * t86 - t238 * t106 - t57 * t210;
t6 = t100 * t87 + t239 * t106 - t60 * t210;
t4 = qJ(5) * t33 - qJD(4) * t240 - qJD(5) * t71 + t180;
t3 = -qJ(5) * t34 - qJD(5) * t70 + t159;
t2 = -qJD(5) * t57 - t177 - t235;
t1 = -qJD(5) * t60 + t152 + t256;
t13 = [qJDD(1), t183, t165, (-t120 + t253) * t137, t170 * t202 + t158 - t165, -t120 * pkin(1) - g(1) * (-pkin(1) * t140 + t126) - g(2) * t213 + (t212 * t202 + t158) * qJ(2), (qJDD(1) * t134 - 0.2e1 * t139 * t187) * t131, (-t139 * t197 + t211 * t200) * t250, (t139 * t161 - t142 * t171) * t136, (t139 * t171 + t142 * t161) * t136, t111 * t137, -g(1) * t76 - g(2) * t78 - t83 * t111 - t64 * t137 + (t113 * t137 + (t250 + t132) * qJD(1)) * qJ(2) * t206 + (qJD(3) * t93 * t113 + t160 * t250 + (qJ(2) * t111 + qJD(2) * t113 + t227 + t254) * t137) * t139, (-t139 * t192 + t231) * t113 + t230 * t111 + t153 * t137 - g(1) * t75 - g(2) * t77 + (t188 + (-t139 * t200 + t197) * qJ(2)) * t250, -t22 * t71 - t33 * t60, t22 * t70 - t23 * t71 + t33 * t57 - t34 * t60, -t100 * t71 + t106 * t33 + t137 * t22, t100 * t70 + t106 * t34 + t137 * t23, t100 * t137, -t180 * t106 - t178 * t100 - t182 * t137 + t80 * t57 + t85 * t23 + t49 * t70 + t74 * t34 + (t106 * t240 - t137 * t163) * qJD(4) + t166, t240 * t100 + t159 * t106 - t177 * t137 - t85 * t22 - t74 * t33 + t49 * t71 + t80 * t60 + t167, -t1 * t137 - t100 * t18 - t106 * t4 + t12 * t70 + t23 * t47 + t25 * t57 + t34 * t40 + t166, t100 * t20 + t106 * t3 + t12 * t71 + t137 * t2 - t22 * t47 + t25 * t60 - t33 * t40 + t167, -t1 * t71 + t183 * t136 + t18 * t22 - t2 * t70 - t20 * t23 - t3 * t57 + t33 * t5 - t34 * t9 - t4 * t60, t2 * t20 + t9 * t3 + t1 * t18 + t5 * t4 + t12 * t47 + t40 * t25 - g(1) * (t143 * t91 + t126) - g(2) * (-t143 * t225 + t88 * t220 + t213) + (-g(1) * (-t137 * t88 - pkin(1) + t225) - g(2) * t91) * t140; 0, 0, 0, -t196, -t176, -qJ(2) * t176 + qJDD(2) - t253, 0, 0, 0, 0, 0, -t111 * t142 + t139 * t155, t111 * t139 + t142 * t155, 0, 0, 0, 0, 0, t7, t6, t7, t6, t22 * t86 - t23 * t87 - t238 * t60 - t239 * t57, t1 * t86 + t2 * t87 - t40 * t210 + t238 * t5 + t239 * t9 - t183; 0, 0, 0, 0, 0, 0, t142 * t139 * t224, -t211 * t224, (-t139 * t174 + t197) * t136, (-t142 * t174 - t198) * t136, -t111, -g(1) * t77 + g(2) * t75 + t64 + (-t137 * t174 - t224) * t228 + (-t203 * t73 + t243 - t254) * t139, g(3) * t222 + g(1) * t78 - g(2) * t76 - t65 * t113 + (t203 * t209 + t224) * t229 - t175, t242, t24, t16, t17, -t100, t179 * t106 + (-t100 * t141 + t106 * t205 - t190 * t57) * pkin(3) + t147, -t241 * t106 - t60 * t172 + t151 + t257, t10 * t106 - t100 * t119 - t50 * t57 - t214 * t60 + (-t37 + (-t28 + t246) * t138) * qJD(4) + t182 + t252 + t256, -t106 * t11 - t50 * t60 + t149 + t257, t119 * t22 + (t10 + t9) * t60 + (t11 - t5) * t57 + (-t138 * t23 + (t138 * t60 - t141 * t57) * qJD(4)) * pkin(3), t1 * t119 - t9 * t11 - t5 * t10 - t40 * t50 - g(1) * (t140 * t92 - t91 * t220) - g(2) * (-t143 * t92 - t91 * t221) + t91 * t243 + (t2 * t138 + (-t138 * t5 + t141 * t9) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t242, t24, t16, t17, -t100, t106 * t163 + t147, -t106 * t181 + t151, t236 - t106 * t9 - 0.2e1 * t94 + (t184 - t40) * t60 + t148, -pkin(4) * t251 - t106 * t8 + t149, pkin(4) * t22 - t249 * t57, t249 * t9 + (-t40 * t60 + t1 + t252) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 - t232, -t22 + t233, -t55 - t251, g(3) * t137 - t136 * t165 + t5 * t60 + t57 * t9 + t12;];
tau_reg = t13;
