% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% tauc_reg [5x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPR11_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR11_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:35:14
% EndTime: 2019-12-31 21:35:23
% DurationCPUTime: 3.04s
% Computational Cost: add. (2180->362), mult. (5356->504), div. (0->0), fcn. (3419->6), ass. (0->190)
t134 = cos(qJ(2));
t204 = qJD(1) * t134;
t257 = qJD(3) - t204;
t247 = qJD(3) - qJD(5);
t129 = sin(qJ(5));
t132 = cos(qJ(5));
t130 = sin(qJ(3));
t133 = cos(qJ(3));
t194 = t133 * qJD(2);
t131 = sin(qJ(2));
t205 = qJD(1) * t131;
t83 = t130 * t205 - t194;
t182 = t133 * t205;
t203 = qJD(2) * t130;
t85 = t182 + t203;
t156 = t129 * t85 - t132 * t83;
t40 = t129 * t83 + t132 * t85;
t256 = -t156 ^ 2 + t40 ^ 2;
t106 = t257 * qJD(4);
t192 = qJD(1) * qJD(2);
t119 = t131 * t192;
t112 = qJ(4) * t119;
t124 = pkin(6) * t204;
t103 = qJD(2) * pkin(7) + t124;
t199 = qJD(3) * t133;
t200 = qJD(3) * t130;
t96 = -pkin(2) * t134 - pkin(7) * t131 - pkin(1);
t76 = t96 * qJD(1);
t161 = pkin(2) * t131 - pkin(7) * t134;
t94 = t161 * qJD(2);
t77 = qJD(1) * t94;
t150 = -t103 * t200 + t130 * t77 + t76 * t199;
t166 = pkin(6) * t119;
t139 = t133 * t166 - t150;
t11 = t106 + t112 - t139;
t177 = t131 * t199;
t201 = qJD(2) * t134;
t142 = t130 * t201 + t177;
t191 = qJD(2) * qJD(3);
t57 = t142 * qJD(1) + t130 * t191;
t6 = pkin(8) * t57 + t11;
t171 = t103 * t199 - t130 * t166 - t133 * t77 + t76 * t200;
t243 = pkin(3) + pkin(4);
t176 = t134 * t192;
t179 = t131 * t200;
t56 = qJD(1) * t179 + (-t176 - t191) * t133;
t7 = pkin(8) * t56 - t243 * t119 + t171;
t183 = -t129 * t6 + t132 * t7;
t102 = -qJD(2) * pkin(2) + pkin(6) * t205;
t148 = qJ(4) * t85 - t102;
t24 = -t243 * t83 + t148;
t255 = t24 * t40 - t183;
t253 = -0.2e1 * t192;
t252 = t24 * t156;
t251 = t40 * t156;
t167 = pkin(3) * t119;
t15 = -t167 + t171;
t108 = t257 * qJ(4);
t43 = t133 * t103 + t130 * t76;
t32 = t108 + t43;
t250 = -t257 * t32 + t15;
t42 = -t130 * t103 + t133 * t76;
t207 = qJD(4) - t42;
t249 = qJD(4) * t130 + t124;
t215 = t130 * t132;
t152 = t129 * t133 - t215;
t193 = -qJD(5) + t257;
t248 = -qJD(5) - t193;
t195 = qJD(5) * t132;
t196 = qJD(5) * t129;
t9 = t129 * t57 - t132 * t56 + t83 * t195 - t85 * t196;
t246 = t156 * t193 - t9;
t222 = qJ(4) * t130;
t245 = -t243 * t133 - t222;
t244 = t85 ^ 2;
t242 = pkin(7) - pkin(8);
t241 = pkin(6) * t130;
t34 = pkin(3) * t83 - t148;
t240 = t34 * t85;
t239 = t85 * t83;
t87 = t129 * t130 + t132 * t133;
t144 = t87 * t134;
t238 = -qJD(1) * t144 + t247 * t87;
t237 = t129 * t199 + t130 * t195 - t132 * t200 - t133 * t196 - t152 * t204;
t221 = qJ(4) * t133;
t159 = pkin(3) * t130 - t221;
t236 = -t257 * t159 + t249;
t145 = -t243 * t130 + t221;
t235 = t257 * t145 + t249;
t234 = t130 * t94 + t96 * t199;
t233 = qJ(4) * t83;
t231 = t257 * t83;
t230 = t257 * t85;
t27 = pkin(8) * t83 + t43;
t22 = t108 + t27;
t229 = t129 * t22;
t140 = -pkin(6) * t176 - qJ(4) * t56 + qJD(4) * t85;
t13 = pkin(3) * t57 - t140;
t228 = t13 * t130;
t227 = t13 * t133;
t226 = t130 * t56;
t91 = t161 * qJD(1);
t225 = t133 * t91;
t212 = t133 * t134;
t118 = pkin(6) * t212;
t224 = t130 * t96 + t118;
t74 = t130 * t91;
t223 = qJ(4) * t205 + t74;
t220 = t102 * t130;
t219 = t102 * t133;
t218 = t257 * t133;
t216 = t130 * t131;
t214 = t130 * t134;
t213 = t131 * t133;
t137 = qJD(1) ^ 2;
t211 = t134 * t137;
t136 = qJD(2) ^ 2;
t210 = t136 * t131;
t209 = t136 * t134;
t208 = pkin(8) * t85 - t207;
t127 = t131 ^ 2;
t206 = -t134 ^ 2 + t127;
t202 = qJD(2) * t131;
t197 = qJD(4) * t133;
t17 = -t243 * t257 - t208;
t190 = t129 * t7 + t132 * t6 + t17 * t195;
t189 = pkin(7) * t257 * t130;
t188 = pkin(7) * t218;
t187 = qJ(4) * t202 + t234;
t186 = pkin(7) * t202;
t185 = pkin(7) * t194;
t105 = t242 * t133;
t184 = -pkin(3) - t241;
t180 = t134 * t194;
t178 = t134 * t200;
t174 = -t129 * t56 - t132 * t57;
t117 = pkin(6) * t214;
t173 = t133 * t96 - t117;
t172 = pkin(1) * t253;
t170 = -t85 + t203;
t169 = t83 + t194;
t168 = t193 ^ 2;
t165 = -qJD(3) * t118 + t133 * t94 - t96 * t200;
t54 = -qJ(4) * t134 + t224;
t104 = t242 * t130;
t164 = -qJD(5) * t104 + (-pkin(6) * t213 + pkin(8) * t214) * qJD(1) + t223 + t242 * t200;
t141 = -pkin(8) * t212 + (-pkin(4) + t184) * t131;
t163 = t141 * qJD(1) - t247 * t105 - t225;
t162 = t184 * t131;
t160 = pkin(3) * t133 + t222;
t2 = t129 * t17 + t132 * t22;
t126 = t134 * pkin(3);
t35 = pkin(4) * t134 + t117 + t126 + (-pkin(8) * t131 - t96) * t133;
t41 = pkin(8) * t216 + t54;
t158 = -t129 * t41 + t132 * t35;
t157 = t129 * t35 + t132 * t41;
t31 = -pkin(3) * t257 + t207;
t155 = -t130 * t32 + t133 * t31;
t154 = qJ(4) * t132 - t129 * t243;
t153 = qJ(4) * t129 + t132 * t243;
t151 = qJD(1) * t127 + t134 * t257;
t149 = pkin(6) + t159;
t147 = -t22 * t196 + t190;
t146 = t257 * t43 - t171;
t70 = t87 * t131;
t143 = -pkin(6) + t145;
t10 = t40 * qJD(5) + t174;
t138 = t257 * t42 + t139;
t95 = -pkin(2) - t160;
t80 = pkin(2) - t245;
t69 = t129 * t213 - t131 * t215;
t64 = t149 * t131;
t55 = t126 - t173;
t53 = t143 * t131;
t49 = pkin(3) * t85 + t233;
t48 = qJD(1) * t162 - t225;
t47 = -pkin(6) * t182 + t223;
t29 = -t243 * t85 - t233;
t28 = -t56 + t231;
t25 = (t160 * qJD(3) - t197) * t131 + t149 * t201;
t23 = qJD(2) * t162 - t165;
t21 = -qJD(4) * t134 + (-t131 * t194 - t178) * pkin(6) + t187;
t20 = t131 * t152 * t247 + qJD(2) * t144;
t19 = qJD(5) * t70 - t142 * t132 + (-t179 + t180) * t129;
t18 = (qJD(3) * t245 + t197) * t131 + t143 * t201;
t14 = (-pkin(6) * qJD(2) + pkin(8) * qJD(3)) * t213 + (-qJD(4) + (-pkin(6) * qJD(3) + pkin(8) * qJD(2)) * t130) * t134 + t187;
t12 = pkin(8) * t179 + t141 * qJD(2) - t165;
t8 = -t243 * t57 + t140;
t1 = t132 * t17 - t229;
t3 = [0, 0, 0, 0.2e1 * t134 * t119, t206 * t253, t209, -t210, 0, -pkin(6) * t209 + t131 * t172, pkin(6) * t210 + t134 * t172, t85 * t180 + (-t133 * t56 - t85 * t200) * t131, (-t130 * t85 - t133 * t83) * t201 + (t226 - t133 * t57 + (t130 * t83 - t133 * t85) * qJD(3)) * t131, -t257 * t179 + t134 * t56 + (t131 * t85 + t151 * t133) * qJD(2), -t257 * t177 + t134 * t57 + (-t151 * t130 - t131 * t83) * qJD(2), (t257 - t204) * t202, t165 * t257 + t171 * t134 + (pkin(6) * t57 + t102 * t199) * t131 + ((pkin(6) * t83 + t220) * t134 + (t173 * qJD(1) + t42 + (t257 + t204) * t241) * t131) * qJD(2), -(-pkin(6) * t178 + t234) * t257 + t150 * t134 + (-pkin(6) * t56 - t102 * t200) * t131 + ((pkin(6) * t85 + t219) * t134 + (pkin(6) * t218 - t224 * qJD(1) - t43) * t131) * qJD(2), -t257 * t23 + t25 * t83 + t57 * t64 + (t34 * t203 + t15) * t134 + (t34 * t199 + t228 + (-qJD(1) * t55 - t31) * qJD(2)) * t131, -t21 * t83 + t23 * t85 - t54 * t57 - t55 * t56 + t155 * t201 + (-t11 * t130 + t133 * t15 + (-t130 * t31 - t133 * t32) * qJD(3)) * t131, t257 * t21 - t25 * t85 + t56 * t64 + (-t194 * t34 - t11) * t134 + (t34 * t200 - t227 + (qJD(1) * t54 + t32) * qJD(2)) * t131, t11 * t54 + t13 * t64 + t15 * t55 + t21 * t32 + t23 * t31 + t25 * t34, t20 * t40 + t70 * t9, -t10 * t70 - t156 * t20 - t19 * t40 - t69 * t9, -t193 * t20 + t134 * t9 + (-qJD(1) * t70 - t40) * t202, -t10 * t134 + t193 * t19 + (qJD(1) * t69 + t156) * t202, (t193 - t204) * t202, -(t12 * t132 - t129 * t14) * t193 + t183 * t134 + t18 * t156 + t53 * t10 + t8 * t69 + t24 * t19 + (-t134 * t2 + t157 * t193) * qJD(5) + (-qJD(1) * t158 - t1) * t202, (qJD(5) * t158 + t12 * t129 + t132 * t14) * t193 - t147 * t134 + t18 * t40 + t53 * t9 + t8 * t70 + t24 * t20 + (qJD(1) * t157 + t2) * t202; 0, 0, 0, -t131 * t211, t206 * t137, 0, 0, 0, t137 * pkin(1) * t131, pkin(1) * t211, t85 * t218 - t226, (-t56 - t231) * t133 + (-t230 - t57) * t130, t257 * t199 + (t170 * t131 - t212 * t257) * qJD(1), -t257 * t200 + (t131 * t169 + t214 * t257) * qJD(1), -t257 * t205, -t91 * t218 - pkin(2) * t57 + (-t188 + t220) * qJD(3) + (-t42 * t131 + (-t102 * t134 - t186) * t130 + (-t134 * t169 - t216 * t257) * pkin(6)) * qJD(1), pkin(2) * t56 + t74 * t257 + (t189 + t219) * qJD(3) + (-t102 * t212 + (t43 - t185) * t131 + (t134 * t170 - t213 * t257) * pkin(6)) * qJD(1), t257 * t48 - t227 + t57 * t95 - t236 * t83 + (t130 * t34 - t188) * qJD(3) + (t131 * t31 + (-t134 * t34 - t186) * t130) * qJD(1), t47 * t83 - t48 * t85 + (t11 + t257 * t31 + (qJD(3) * t85 - t57) * pkin(7)) * t133 + ((qJD(3) * t83 - t56) * pkin(7) + t250) * t130, -t257 * t47 - t228 + t56 * t95 + t236 * t85 + (-t133 * t34 - t189) * qJD(3) + (t34 * t212 + (-t32 + t185) * t131) * qJD(1), t13 * t95 - t31 * t48 - t32 * t47 - t236 * t34 + (qJD(3) * t155 + t11 * t133 + t130 * t15) * pkin(7), -t152 * t9 + t238 * t40, t10 * t152 - t156 * t238 - t237 * t40 - t87 * t9, -t238 * t193 + (qJD(2) * t152 + t40) * t205, t237 * t193 + (qJD(2) * t87 - t156) * t205, -t193 * t205, t80 * t10 + t8 * t87 + t235 * t156 + t237 * t24 - (t129 * t164 - t132 * t163) * t193 + (-(t104 * t132 - t105 * t129) * qJD(2) + t1) * t205, -t8 * t152 + t80 * t9 + t235 * t40 + t238 * t24 - (t129 * t163 + t132 * t164) * t193 + ((t104 * t129 + t105 * t132) * qJD(2) - t2) * t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t239, -t83 ^ 2 + t244, t28, -t57 + t230, t119, -t102 * t85 + t146, t102 * t83 + t138, -t49 * t83 + t146 + 0.2e1 * t167 - t240, pkin(3) * t56 - qJ(4) * t57 + (t32 - t43) * t85 + (t31 - t207) * t83, -t34 * t83 + t49 * t85 + 0.2e1 * t106 + 0.2e1 * t112 - t138, -pkin(3) * t15 + qJ(4) * t11 + t207 * t32 - t31 * t43 - t34 * t49, -t251, -t256, t246, t193 * t40 + t10, t119, t153 * t119 - t29 * t156 - (t129 * t208 - t132 * t27) * t193 + (t154 * t193 + t2) * qJD(5) + t255, t154 * t119 - t29 * t40 - t252 - (t129 * t27 + t132 * t208) * t193 + (-t153 * t193 - t229) * qJD(5) + t190; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119 + t239, t28, -t257 ^ 2 - t244, t240 + t250, 0, 0, 0, 0, 0, -t119 * t132 - t129 * t168 - t156 * t85, t119 * t129 - t132 * t168 - t40 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t251, t256, -t246, t248 * t40 - t174, -t119, t2 * t248 - t255, -t1 * t193 - t147 + t252;];
tauc_reg = t3;
