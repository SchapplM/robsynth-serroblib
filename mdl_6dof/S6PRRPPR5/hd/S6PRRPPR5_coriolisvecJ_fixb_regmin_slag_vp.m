% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRPPR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:21:05
% EndTime: 2019-03-08 21:21:12
% DurationCPUTime: 2.58s
% Computational Cost: add. (2083->329), mult. (5253->476), div. (0->0), fcn. (3715->10), ass. (0->192)
t144 = sin(qJ(3));
t223 = qJD(2) * t144;
t127 = qJD(6) + t223;
t264 = qJD(6) - t127;
t254 = pkin(4) + pkin(8);
t138 = sin(pkin(11));
t140 = cos(pkin(11));
t143 = sin(qJ(6));
t146 = cos(qJ(6));
t98 = t138 * t146 + t140 * t143;
t157 = t98 * t144;
t83 = t98 * qJD(6);
t248 = -qJD(2) * t157 - t83;
t147 = cos(qJ(3));
t221 = qJD(2) * t147;
t200 = t140 * t221;
t220 = qJD(3) * t138;
t94 = t200 + t220;
t198 = t138 * t221;
t219 = qJD(3) * t140;
t96 = -t198 + t219;
t170 = t143 * t94 - t146 * t96;
t263 = t127 * t170;
t145 = sin(qJ(2));
t139 = sin(pkin(6));
t227 = qJD(1) * t139;
t205 = t145 * t227;
t112 = qJD(2) * pkin(8) + t205;
t141 = cos(pkin(6));
t226 = qJD(1) * t141;
t67 = t144 * t112 - t147 * t226;
t262 = qJD(4) + t67;
t148 = cos(qJ(2));
t224 = qJD(2) * t139;
t195 = qJD(1) * t224;
t183 = t148 * t195;
t218 = qJD(3) * t141;
t261 = qJD(1) * t218 + t183;
t216 = qJD(3) * t147;
t104 = t254 * t216;
t234 = t144 * t148;
t217 = qJD(3) * t144;
t130 = pkin(3) * t217;
t169 = -qJ(4) * t147 + qJ(5) * t144;
t215 = qJD(4) * t144;
t151 = qJD(3) * t169 - qJD(5) * t147 - t215;
t60 = t130 + t151;
t250 = t140 * t104 - t138 * t60 - (-t138 * t145 + t140 * t234) * t227;
t249 = t138 * t104 + t140 * t60 - (t138 * t234 + t140 * t145) * t227;
t135 = qJD(3) * qJ(4);
t68 = t147 * t112 + t144 * t226;
t63 = -t135 - t68;
t62 = -qJD(3) * pkin(3) + t262;
t230 = pkin(4) * t223 + t262;
t259 = qJD(5) + t135;
t142 = -pkin(3) - qJ(5);
t129 = pkin(4) * t221;
t59 = t129 + t68;
t45 = t59 + t259;
t258 = t144 * (-t45 + t259) - t142 * t216;
t149 = qJD(3) ^ 2;
t158 = -qJ(4) * t216 - t215;
t212 = qJD(2) * qJD(3);
t193 = t144 * t212;
t229 = pkin(3) * t193 + t145 * t195;
t52 = qJD(2) * t158 + t229;
t80 = t130 + t158;
t257 = qJD(2) * (-t80 + t205) - pkin(8) * t149 - t52;
t202 = t148 * t224;
t150 = qJD(2) ^ 2;
t236 = t139 * t150;
t209 = t145 * t236;
t238 = t139 * t145;
t210 = t144 * t238;
t55 = -qJD(3) * t210 + (t202 + t218) * t147;
t256 = qJD(3) * (t147 * t202 + t55) - t144 * t209;
t185 = t144 * t202;
t88 = t141 * t144 + t147 * t238;
t56 = qJD(3) * t88 + t185;
t255 = qJD(3) * (t56 + t185) + t147 * t209;
t181 = t146 * t193;
t182 = t143 * t193;
t16 = -qJD(6) * t170 + t138 * t182 - t140 * t181;
t253 = -pkin(9) + t142;
t239 = t138 * t144;
t164 = pkin(5) * t147 - pkin(9) * t239;
t155 = t164 * qJD(3);
t252 = t155 + t250;
t189 = pkin(9) * t140 * t217;
t251 = -t189 - t249;
t39 = t112 * t216 + t261 * t144;
t25 = (-qJD(5) + t129) * qJD(3) + t39;
t31 = qJD(2) * t151 + t229;
t7 = t138 * t25 + t140 * t31;
t41 = t142 * qJD(3) + t230;
t204 = t148 * t227;
t191 = -t144 * qJ(4) - pkin(2);
t91 = t142 * t147 + t191;
t57 = qJD(2) * t91 - t204;
t11 = t138 * t41 + t140 * t57;
t131 = pkin(3) * t223;
t78 = qJD(2) * t169 + t131;
t20 = t138 * t59 + t140 * t78;
t201 = t140 * t223;
t213 = qJD(6) * t146;
t214 = qJD(6) * t143;
t240 = t138 * t143;
t247 = -t138 * t214 + t140 * t213 + t146 * t201 - t223 * t240;
t246 = qJD(2) * pkin(2);
t42 = t143 * t96 + t146 * t94;
t245 = t127 * t42;
t134 = qJD(3) * qJD(4);
t211 = t112 * t217 - t261 * t147;
t188 = t134 - t211;
t24 = -pkin(4) * t193 + t188;
t244 = t138 * t24;
t243 = t140 * t24;
t242 = t144 * t88;
t117 = t254 * t144;
t50 = t138 * t117 + t140 * t91;
t136 = t144 ^ 2;
t241 = t136 * t150;
t237 = t139 * t148;
t235 = t140 * t147;
t233 = t149 * t144;
t232 = t149 * t147;
t207 = -pkin(5) * t140 - pkin(4);
t161 = t207 * t223;
t231 = -t161 + t262;
t118 = t254 * t147;
t137 = t147 ^ 2;
t228 = t136 - t137;
t114 = -pkin(3) * t147 + t191;
t225 = qJD(2) * t114;
t222 = qJD(2) * t145;
t208 = t144 * t150 * t147;
t6 = -t138 * t31 + t140 * t25;
t4 = qJD(2) * t155 + t6;
t5 = qJD(2) * t189 + t7;
t206 = -t143 * t5 + t146 * t4;
t203 = t139 * t222;
t196 = t247 * t127;
t192 = t147 * t212;
t10 = -t138 * t57 + t140 * t41;
t19 = -t138 * t78 + t140 * t59;
t168 = -t140 * t146 + t240;
t184 = t248 * t127 - t168 * t192;
t179 = t138 * t7 + t140 * t6;
t178 = t143 * t4 + t146 * t5;
t8 = pkin(5) * t223 - pkin(9) * t96 + t10;
t9 = -pkin(9) * t94 + t11;
t177 = t143 * t9 - t146 * t8;
t2 = t143 * t8 + t146 * t9;
t175 = -t10 * t138 + t11 * t140;
t101 = t140 * t117;
t32 = t144 * pkin(5) + t101 + (pkin(9) * t147 - t91) * t138;
t38 = -pkin(9) * t235 + t50;
t174 = -t143 * t38 + t146 * t32;
t173 = t143 * t32 + t146 * t38;
t87 = -t141 * t147 + t210;
t53 = t138 * t237 + t140 * t87;
t54 = t138 * t87 - t140 * t237;
t172 = -t143 * t54 + t146 * t53;
t171 = t143 * t53 + t146 * t54;
t166 = qJD(3) * t68 - t39;
t165 = qJD(3) * t67 - t211;
t105 = t253 * t138;
t160 = t164 * qJD(2) + qJD(5) * t140 + qJD(6) * t105 + t19;
t106 = t253 * t140;
t159 = pkin(9) * t201 + qJD(5) * t138 - qJD(6) * t106 + t20;
t15 = t138 * t181 + t140 * t182 - t94 * t213 - t96 * t214;
t75 = t168 * t147;
t113 = -t204 - t246;
t154 = qJD(3) * (t113 + t204 - t246);
t69 = -t204 + t225;
t153 = qJD(3) * (-t204 - t69 - t225);
t152 = t144 * t39 + t147 * t188 + (t144 * t63 + t147 * t62) * qJD(3);
t128 = pkin(5) * t138 + qJ(4);
t103 = t254 * t217;
t102 = -qJ(4) * t221 + t131;
t82 = pkin(5) * t235 + t118;
t76 = t98 * t147;
t73 = (-pkin(8) + t207) * t217;
t61 = t69 * t223;
t49 = -t138 * t91 + t101;
t34 = t147 * t83 - t168 * t217;
t33 = qJD(3) * t157 + qJD(6) * t75;
t30 = t138 * t56 + t140 * t203;
t29 = -t138 * t203 + t140 * t56;
t23 = pkin(5) * t94 + t45;
t17 = qJD(3) * t161 + t188;
t1 = [0, 0, -t209, -t148 * t236, 0, 0, 0, 0, 0, -t255, -t256 (t144 * t56 + t147 * t55 + (t147 * t87 - t242) * qJD(3)) * qJD(2), t255, t256, t188 * t88 + t39 * t87 - t55 * t63 + t56 * t62 + (-t148 * t52 + t69 * t222) * t139, t55 * t94 + (t144 * t29 + (-t140 * t242 + t147 * t53) * qJD(3)) * qJD(2), t55 * t96 + (-t144 * t30 + (-t147 * t54 + t88 * t239) * qJD(3)) * qJD(2), -t29 * t96 - t30 * t94 + (-t138 * t53 + t140 * t54) * t193, t10 * t29 + t11 * t30 + t24 * t88 + t45 * t55 + t53 * t6 + t54 * t7, 0, 0, 0, 0, 0 (-qJD(6) * t171 - t143 * t30 + t146 * t29) * t127 + t172 * t192 + t55 * t42 + t88 * t16 -(qJD(6) * t172 + t143 * t29 + t146 * t30) * t127 - t171 * t192 - t55 * t170 + t88 * t15; 0, 0, 0, 0, 0.2e1 * t144 * t192, -0.2e1 * t228 * t212, t232, -t233, 0, -pkin(8) * t232 + t144 * t154, pkin(8) * t233 + t147 * t154 (-t136 - t137) * t183 + t152, t144 * t153 - t257 * t147, t257 * t144 + t147 * t153, t114 * t52 + t69 * t80 + (-t145 * t69 + (-t144 * t62 + t147 * t63) * t148) * t227 + t152 * pkin(8), -t103 * t94 + (-t94 * t204 + t243 + (qJD(2) * t49 + t10) * qJD(3)) * t147 + (-t45 * t219 + t6 + (-t118 * t219 + t250) * qJD(2)) * t144, -t103 * t96 + (-t96 * t204 - t244 + (-qJD(2) * t50 - t11) * qJD(3)) * t147 + (t45 * t220 - t7 + (t118 * t220 - t249) * qJD(2)) * t144, -t250 * t96 - t249 * t94 + (t138 * t6 - t140 * t7) * t147 + ((-t138 * t49 + t140 * t50) * qJD(2) + t175) * t217, t118 * t24 + t49 * t6 + t50 * t7 + (-t147 * t204 - t103) * t45 + t249 * t11 + t250 * t10, -t15 * t76 - t170 * t33, t15 * t75 + t16 * t76 - t170 * t34 - t33 * t42, t127 * t33 + t144 * t15 + (-qJD(2) * t76 - t170) * t216, t127 * t34 - t144 * t16 + (qJD(2) * t75 - t42) * t216 (t127 + t223) * t216, t206 * t144 + t73 * t42 + t82 * t16 - t17 * t75 - t23 * t34 + (t251 * t143 + t252 * t146) * t127 + (-t127 * t173 - t144 * t2) * qJD(6) + (-t42 * t204 + (qJD(2) * t174 - t177) * qJD(3)) * t147, -t178 * t144 - t73 * t170 + t82 * t15 - t17 * t76 + t23 * t33 + (-t252 * t143 + t251 * t146) * t127 + (-t127 * t174 + t144 * t177) * qJD(6) + (t170 * t204 + (-qJD(2) * t173 - t2) * qJD(3)) * t147; 0, 0, 0, 0, -t208, t228 * t150, 0, 0, 0, -t113 * t223 + t166, -t113 * t221 - t165, 0, -t102 * t221 - t166 + t61, 0.2e1 * t134 + (t102 * t144 + t147 * t69) * qJD(2) + t165, -pkin(3) * t39 + qJ(4) * t188 - t102 * t69 - t262 * t63 - t62 * t68, t244 + t230 * t94 + (-t10 * t147 - t258 * t140 - t144 * t19) * qJD(2), t243 + t230 * t96 + (t11 * t147 + t258 * t138 + t144 * t20) * qJD(2), t19 * t96 + t20 * t94 + (qJD(5) * t96 - t11 * t223 - t6) * t140 + (qJD(5) * t94 + t10 * t223 - t7) * t138, qJ(4) * t24 - t10 * t19 - t11 * t20 + t230 * t45 + t179 * t142 + (-t10 * t140 - t11 * t138) * qJD(5), -t15 * t168 - t170 * t248, -t15 * t98 + t16 * t168 + t170 * t247 - t248 * t42, t170 * t221 + t184, -t196 + (-qJD(3) * t98 + t42) * t221, -t127 * t221, t128 * t16 + t17 * t98 + t231 * t42 + t247 * t23 + (t143 * t159 - t146 * t160) * t127 + ((-t105 * t143 + t106 * t146) * qJD(3) + t177) * t221, t128 * t15 - t17 * t168 - t231 * t170 + t248 * t23 + (t143 * t160 + t146 * t159) * t127 + (-(t105 * t146 + t106 * t143) * qJD(3) + t2) * t221; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t208, -t149 - t241, qJD(3) * t63 + t39 + t61, -t138 * t241 + (-t94 + t200) * qJD(3), -t140 * t241 + (-t96 - t198) * qJD(3) (t138 * t96 - t140 * t94) * t223, -qJD(3) * t45 + t175 * t223 + t179, 0, 0, 0, 0, 0, -qJD(3) * t42 + t184, -t196 + (-t221 * t98 + t170) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t96 - t219) * t223 (-t94 + t220) * t223, -t94 ^ 2 - t96 ^ 2, t10 * t96 + t11 * t94 + t24, 0, 0, 0, 0, 0, t16 - t263, t15 - t245; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t170 * t42, t170 ^ 2 - t42 ^ 2, t15 + t245, -t16 - t263, t192, t23 * t170 - t2 * t264 + t206, t177 * t264 + t23 * t42 - t178;];
tauc_reg  = t1;
