% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRR6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:05:54
% EndTime: 2020-01-03 12:06:01
% DurationCPUTime: 2.04s
% Computational Cost: add. (4317->258), mult. (7305->381), div. (0->0), fcn. (4658->8), ass. (0->191)
t146 = sin(qJ(4));
t147 = sin(qJ(2));
t149 = cos(qJ(4));
t144 = cos(pkin(9));
t205 = qJD(3) * t144;
t150 = cos(qJ(2));
t213 = t144 * t150;
t232 = pkin(1) * qJD(1);
t143 = sin(pkin(9));
t120 = -t144 * pkin(3) - t143 * pkin(7) - pkin(2);
t214 = t144 * t149;
t83 = qJ(3) * t214 + t146 * t120;
t247 = qJD(4) * t83;
t237 = -t146 * t205 - t247 - (-t146 * t213 + t147 * t149) * t232;
t203 = qJD(4) * t149;
t251 = t120 * t203 + t149 * t205 - (t146 * t147 + t149 * t213) * t232;
t215 = t144 * t146;
t190 = qJ(3) * t215;
t217 = t143 * t149;
t199 = pkin(8) * t217;
t250 = (-t190 - t199) * qJD(4) + t251;
t204 = qJD(4) * t146;
t183 = t143 * t204;
t128 = pkin(8) * t183;
t249 = t128 + t237;
t145 = sin(qJ(5));
t212 = t145 * t146;
t191 = t143 * t212;
t248 = -qJD(5) * t191 - t145 * t183;
t196 = t150 * t232;
t165 = qJD(3) - t196;
t148 = cos(qJ(5));
t115 = t145 * t149 + t148 * t146;
t92 = t115 * t143;
t140 = qJD(1) + qJD(2);
t118 = t140 * qJ(3) + t147 * t232;
t189 = t118 * t215;
t192 = t140 * t217;
t153 = -pkin(8) * t192 - t189;
t231 = pkin(1) * qJD(2);
t194 = qJD(1) * t231;
t111 = t140 * qJD(3) + t150 * t194;
t174 = t147 * t194;
t72 = t120 * t140 + t165;
t198 = t111 * t214 + t146 * t174 + t72 * t203;
t22 = t153 * qJD(4) + t198;
t216 = t144 * t140;
t127 = -qJD(4) + t216;
t63 = t149 * t72;
t40 = t153 + t63;
t32 = -t127 * pkin(4) + t40;
t246 = (qJD(5) * t32 + t22) * t148;
t202 = qJD(4) + qJD(5);
t245 = t147 * pkin(1);
t244 = t150 * pkin(1);
t73 = (pkin(4) * t140 * t146 + t118) * t143;
t211 = t148 * t149;
t187 = t143 * t211;
t173 = t140 * t187;
t78 = -t140 * t191 + t173;
t243 = t73 * t78;
t154 = t140 * t115;
t75 = t143 * t154;
t242 = t78 * t75;
t108 = t149 * t120;
t59 = -t199 + t108 + (-qJ(3) * t146 - pkin(4)) * t144;
t218 = t143 * t146;
t200 = pkin(8) * t218;
t71 = -t200 + t83;
t33 = -t145 * t71 + t148 * t59;
t241 = t33 * qJD(5) + t249 * t145 + t250 * t148;
t34 = t145 * t59 + t148 * t71;
t240 = -t34 * qJD(5) - t250 * t145 + t249 * t148;
t138 = t143 ^ 2;
t220 = t138 * t149;
t182 = t144 * t204;
t29 = -t118 * t182 + t198;
t239 = t111 * t220 + t29 * t144;
t238 = -qJ(3) * t182 + t251;
t114 = t211 - t212;
t236 = (t202 - t216) * t114;
t69 = t202 * t115;
t235 = t144 * t154 - t69;
t132 = t150 * t231 + qJD(3);
t136 = qJ(3) + t245;
t225 = t118 * t138;
t98 = t138 * t111;
t234 = t132 * t225 + t136 * t98;
t139 = t144 ^ 2;
t233 = t139 * t111 + t98;
t186 = t118 * t214;
t48 = t146 * t72 + t186;
t41 = -t140 * t200 + t48;
t230 = t145 * t41;
t229 = t148 * t41;
t228 = t29 * t146;
t206 = qJD(3) * t118;
t227 = qJ(3) * t98 + t138 * t206;
t102 = t120 - t244;
t67 = t146 * t102 + t136 * t214;
t226 = qJD(4) * t48;
t224 = t132 * t140;
t137 = t140 ^ 2;
t223 = t138 * t137;
t222 = t138 * t140;
t221 = t138 * t146;
t219 = t140 * t143;
t209 = t248 * t140;
t208 = t138 + t139;
t207 = t146 ^ 2 - t149 ^ 2;
t164 = -t111 * t215 + t149 * t174;
t23 = (-t186 + (pkin(8) * t219 - t72) * t146) * qJD(4) + t164;
t179 = -qJD(5) * t230 + t145 * t23;
t4 = t179 + t246;
t54 = t202 * t92;
t184 = t140 * t203;
t70 = (pkin(4) * t184 + t111) * t143;
t93 = t114 * t143;
t201 = t4 * t144 - t73 * t54 + t70 * t93;
t197 = t147 * t231;
t195 = t102 * t203 + t132 * t214 + t146 * t197;
t193 = t118 * t220;
t188 = t136 * t215;
t15 = t148 * t32 - t230;
t16 = t145 * t32 + t229;
t180 = -t145 * t22 + t148 * t23;
t5 = -t16 * qJD(5) + t180;
t156 = t202 * t187;
t55 = t156 + t248;
t185 = t15 * t54 - t16 * t55 - t4 * t92 - t5 * t93;
t119 = -qJD(5) + t127;
t181 = pkin(4) * t119 - t32;
t177 = t208 * t150;
t176 = t219 * t245;
t175 = pkin(4) * t192;
t172 = t137 * t146 * t220;
t171 = -t5 * t144 + t73 * t55 + t70 * t92;
t30 = t164 - t226;
t168 = qJD(4) * t193 + t111 * t221 - t30 * t144;
t167 = (-qJD(2) + t140) * t232;
t166 = (-qJD(1) - t140) * t231;
t163 = (-qJD(4) - t127) * t219;
t96 = t149 * t102;
t53 = -t199 + t96 + (-t136 * t146 - pkin(4)) * t144;
t58 = -t200 + t67;
t26 = -t145 * t58 + t148 * t53;
t27 = t145 * t53 + t148 * t58;
t47 = t63 - t189;
t162 = t146 * t47 - t149 * t48;
t161 = t184 * t221;
t129 = t143 * pkin(4) * t203;
t160 = t165 * t143 + t129;
t159 = t147 * t166;
t158 = t147 * t167;
t157 = t73 * t75 - t179;
t155 = qJD(4) * t143 * (t127 + t216);
t152 = -t228 + (-t30 - t226) * t149;
t151 = -t127 ^ 2 - t223;
t44 = -t67 * qJD(4) - t132 * t215 + t149 * t197;
t135 = pkin(4) * t218;
t124 = t143 * t174;
t113 = t143 * qJ(3) + t135;
t112 = -t140 * pkin(2) + t165;
t104 = -0.2e1 * t161;
t103 = 0.2e1 * t161;
t97 = t143 * t136 + t135;
t84 = t143 * t132 + t129;
t82 = t108 - t190;
t74 = 0.2e1 * t207 * qJD(4) * t222;
t66 = t96 - t188;
t65 = t149 * t155;
t64 = t146 * t155;
t46 = t140 * t156 + t209;
t45 = t140 * t54;
t43 = -t136 * t182 + t195;
t42 = t47 * t183;
t39 = t128 + t44;
t38 = (-t188 - t199) * qJD(4) + t195;
t31 = -t75 ^ 2 + t78 ^ 2;
t25 = -t78 * t119 - t202 * t173 - t209;
t24 = -t75 * t119 - t69 * t219;
t20 = t55 * t119 + t46 * t144;
t19 = t54 * t119 + t45 * t144;
t18 = t148 * t40 - t230;
t17 = -t145 * t40 - t229;
t14 = t46 * t92 + t75 * t55;
t13 = -t45 * t93 - t78 * t54;
t8 = -t27 * qJD(5) - t145 * t38 + t148 * t39;
t7 = t26 * qJD(5) + t145 * t39 + t148 * t38;
t6 = t45 * t92 - t93 * t46 + t54 * t75 - t78 * t55;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159, t150 * t166, 0, 0, 0, 0, 0, 0, 0, 0, t144 * t159, qJD(2) * t176 + t124, t208 * t224 + t233, (t111 * t136 + t118 * t132) * t139 + (t112 + (-pkin(2) - t244) * qJD(1)) * t197 + t234, t104, t74, t64, t103, t65, 0, -t44 * t127 + (t132 * t146 + t136 * t203) * t222 + t168, t43 * t127 + (t149 * t224 + (-t136 * t140 - t118) * t204) * t138 + t239, t42 + ((-t146 * t43 - t149 * t44 + (t146 * t66 - t149 * t67) * qJD(4)) * t140 + t152) * t143, t29 * t67 + t30 * t66 + t48 * t43 + t47 * t44 + t234, t13, t6, t19, t14, t20, 0, -t8 * t119 + t97 * t46 + t84 * t75 + t171, t7 * t119 - t97 * t45 + t84 * t78 + t201, t26 * t45 - t27 * t46 - t7 * t75 - t8 * t78 + t185, t15 * t8 + t16 * t7 + t5 * t26 + t4 * t27 + t70 * t97 + t73 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158, t150 * t167, 0, 0, 0, 0, 0, 0, 0, 0, t144 * t158, -qJD(1) * t176 + t124, (t208 * qJD(3) - t177 * t232) * t140 + t233, (qJ(3) * t111 + t206) * t139 + ((-pkin(2) * qJD(2) - t112) * t147 - t118 * t177) * t232 + t227, t104, t74, t64, t103, t65, 0, -t237 * t127 + (qJ(3) * t203 + t165 * t146) * t222 + t168, t238 * t127 + (-t118 * t204 + (-qJ(3) * t204 + t165 * t149) * t140) * t138 + t239, t42 + (((-t237 - t247) * t149 + (qJD(4) * t82 - t238) * t146) * t140 + t152) * t143, -t196 * t225 + t237 * t47 + t238 * t48 + t29 * t83 + t30 * t82 + t227, t13, t6, t19, t14, t20, 0, t113 * t46 - t240 * t119 + t160 * t75 + t171, -t113 * t45 + t241 * t119 + t160 * t78 + t201, -t240 * t78 - t241 * t75 + t33 * t45 - t34 * t46 + t185, t70 * t113 + t240 * t15 + t241 * t16 + t160 * t73 + t5 * t33 + t4 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t208 * t137, -t208 * t140 * t118 + t174, 0, 0, 0, 0, 0, 0, t151 * t146, t151 * t149, 0, t228 + t30 * t149 - t162 * qJD(4) + (t162 * t144 - t225) * t140, 0, 0, 0, 0, 0, 0, -t235 * t119 - t75 * t219, t236 * t119 - t78 * t219, t114 * t45 - t115 * t46 - t235 * t78 - t236 * t75, t5 * t114 + t4 * t115 + t235 * t15 + t236 * t16 - t73 * t219; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t172, -t207 * t223, t146 * t163, -t172, t149 * t163, 0, -t48 * t127 - t140 * t193 + t30, -t47 * t127 + (qJD(4) * t144 + t222) * t146 * t118 - t198, 0, 0, t242, t31, t24, -t242, t25, 0, -t75 * t175 + t17 * t119 - t243 + (t181 * t145 - t229) * qJD(5) + t180, -t78 * t175 - t18 * t119 + (t181 * qJD(5) - t22) * t148 + t157, (t16 + t17) * t78 + (-t15 + t18) * t75 + (-t145 * t46 + t148 * t45 + (t145 * t78 - t148 * t75) * qJD(5)) * pkin(4), -t15 * t17 - t16 * t18 + (-t73 * t192 + t145 * t4 + t148 * t5 + (-t145 * t15 + t148 * t16) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t242, t31, t24, -t242, t25, 0, -t16 * t119 - t243 + t5, -t15 * t119 + t157 - t246, 0, 0;];
tauc_reg = t1;
