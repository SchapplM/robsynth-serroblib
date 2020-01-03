% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauc_reg [5x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPR7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:17:16
% EndTime: 2019-12-31 21:17:25
% DurationCPUTime: 2.78s
% Computational Cost: add. (3750->293), mult. (9433->428), div. (0->0), fcn. (6836->8), ass. (0->175)
t158 = qJD(2) + qJD(3);
t164 = sin(qJ(3));
t167 = cos(qJ(2));
t238 = cos(qJ(3));
t201 = qJD(1) * t238;
t165 = sin(qJ(2));
t213 = qJD(1) * t165;
t249 = -t164 * t213 + t167 * t201;
t87 = t249 * t158;
t114 = qJD(5) - t249;
t248 = qJD(5) - t114;
t162 = cos(pkin(9));
t166 = cos(qJ(5));
t218 = t166 * t162;
t161 = sin(pkin(9));
t163 = sin(qJ(5));
t221 = t163 * t161;
t129 = -t218 + t221;
t230 = t114 * t129;
t219 = t166 * t161;
t130 = t163 * t162 + t219;
t247 = t114 * t130;
t220 = t164 * t167;
t121 = -qJD(1) * t220 - t165 * t201;
t102 = t162 * t121 - t161 * t158;
t104 = t161 * t121 + t162 * t158;
t183 = t166 * t102 - t163 * t104;
t246 = t114 * t183;
t245 = t166 * t104;
t209 = qJD(1) * qJD(2);
t244 = -0.2e1 * t209;
t199 = t165 * t209;
t132 = t238 * t165 + t220;
t99 = t158 * t132;
t88 = t99 * qJD(1);
t32 = pkin(2) * t199 + t88 * pkin(3) - t87 * qJ(4) + t121 * qJD(4);
t239 = -pkin(7) - pkin(6);
t142 = t239 * t165;
t134 = qJD(1) * t142;
t235 = qJD(2) * pkin(2);
t126 = t134 + t235;
t204 = qJD(2) * t239;
t196 = qJD(1) * t204;
t127 = t165 * t196;
t128 = t167 * t196;
t143 = t239 * t167;
t136 = qJD(1) * t143;
t200 = t238 * qJD(3);
t212 = qJD(3) * t164;
t171 = t126 * t200 + t238 * t127 + t164 * t128 + t136 * t212;
t48 = t158 * qJD(4) + t171;
t12 = -t161 * t48 + t162 * t32;
t13 = t161 * t32 + t162 * t48;
t191 = -t12 * t161 + t13 * t162;
t146 = pkin(2) * t200 + qJD(4);
t89 = -t121 * pkin(3) - qJ(4) * t249;
t78 = pkin(2) * t213 + t89;
t122 = t164 * t136;
t96 = t238 * t134 + t122;
t46 = -t161 * t96 + t162 * t78;
t243 = t146 * t161 + t46;
t47 = t161 * t78 + t162 * t96;
t242 = -t146 * t162 + t47;
t123 = t238 * t136;
t95 = t164 * t134 - t123;
t194 = pkin(2) * t212 - t95;
t241 = t238 * t142 + t164 * t143;
t240 = qJD(1) * t132;
t18 = -t183 * qJD(5) + t130 * t87;
t237 = t162 * pkin(4);
t157 = t162 * pkin(8);
t207 = t165 * t235;
t180 = -t164 * t165 + t238 * t167;
t98 = t158 * t180;
t44 = t99 * pkin(3) - t98 * qJ(4) - t132 * qJD(4) + t207;
t135 = t165 * t204;
t137 = t167 * t204;
t60 = t241 * qJD(3) + t238 * t135 + t164 * t137;
t16 = t161 * t44 + t162 * t60;
t156 = -t167 * pkin(2) - pkin(1);
t141 = t156 * qJD(1);
t73 = -pkin(3) * t249 + t121 * qJ(4) + t141;
t93 = t164 * t126 - t123;
t79 = t158 * qJ(4) + t93;
t38 = t161 * t73 + t162 * t79;
t210 = qJD(5) * t166;
t236 = t104 * t210 + t87 * t218;
t92 = t238 * t126 + t122;
t50 = t161 * t89 + t162 * t92;
t233 = t161 * t87;
t232 = t161 * t98;
t231 = t162 * t87;
t106 = t164 * t142 - t238 * t143;
t91 = -pkin(3) * t180 - t132 * qJ(4) + t156;
t55 = t162 * t106 + t161 * t91;
t228 = t114 * t121;
t227 = t249 * t161;
t226 = t249 * t162;
t225 = t121 * t249;
t224 = t132 * t161;
t169 = qJD(1) ^ 2;
t217 = t167 * t169;
t168 = qJD(2) ^ 2;
t216 = t168 * t165;
t215 = t168 * t167;
t214 = t165 ^ 2 - t167 ^ 2;
t211 = qJD(5) * t132;
t208 = pkin(8) * t227;
t4 = t88 * pkin(4) - pkin(8) * t231 + t12;
t7 = -pkin(8) * t233 + t13;
t205 = -t163 * t7 + t166 * t4;
t52 = t126 * t212 + t164 * t127 - t238 * t128 - t136 * t200;
t198 = -t38 * t121 + t52 * t161;
t15 = -t161 * t60 + t162 * t44;
t37 = -t161 * t79 + t162 * t73;
t49 = -t161 * t92 + t162 * t89;
t54 = -t161 * t106 + t162 * t91;
t197 = pkin(1) * t244;
t155 = -t238 * pkin(2) - pkin(3);
t108 = pkin(4) * t227;
t195 = -t108 + t194;
t193 = -t121 * pkin(4) - pkin(8) * t226;
t192 = t163 * t4 + t166 * t7;
t190 = t37 * t121 - t52 * t162;
t189 = -t161 * t37 + t162 * t38;
t19 = -pkin(4) * t249 + pkin(8) * t102 + t37;
t22 = pkin(8) * t104 + t38;
t188 = t163 * t22 - t166 * t19;
t6 = t163 * t19 + t166 * t22;
t34 = -pkin(4) * t180 - t132 * t157 + t54;
t43 = -pkin(8) * t224 + t55;
t187 = -t163 * t43 + t166 * t34;
t186 = t163 * t34 + t166 * t43;
t185 = t37 * t226 + t38 * t227 + t191;
t184 = qJD(5) * t102 - t233;
t152 = t164 * pkin(2) + qJ(4);
t124 = (-pkin(8) - t152) * t161;
t182 = -qJD(5) * t124 - t208 + t242;
t125 = t162 * t152 + t157;
t181 = qJD(5) * t125 + t193 + t243;
t139 = (-pkin(8) - qJ(4)) * t161;
t179 = -qJD(4) * t162 - qJD(5) * t139 - t208 + t50;
t140 = t162 * qJ(4) + t157;
t178 = qJD(4) * t161 + qJD(5) * t140 + t193 + t49;
t27 = pkin(4) * t233 + t52;
t77 = -t158 * pkin(3) + qJD(4) - t92;
t53 = -pkin(4) * t104 + t77;
t177 = -t188 * t121 + t27 * t129 + t247 * t53;
t176 = -t6 * t121 + t27 * t130 - t230 * t53;
t175 = t141 * t121 - t52;
t174 = t132 * t52 - t241 * t87 + t77 * t98;
t173 = -pkin(3) * t87 - qJ(4) * t88 - (-qJD(4) + t77) * t249;
t172 = -t152 * t88 + t155 * t87 - (-t146 + t77) * t249;
t17 = t184 * t163 + t236;
t170 = -t141 * t249 - t171;
t61 = t106 * qJD(3) + t164 * t135 - t238 * t137;
t153 = -pkin(3) - t237;
t138 = t155 - t237;
t83 = t129 * t132;
t82 = t130 * t132;
t75 = pkin(4) * t224 - t241;
t66 = t121 ^ 2 - t249 ^ 2;
t64 = t108 + t93;
t63 = (-t121 - t240) * t158;
t57 = -t102 * t163 - t245;
t36 = pkin(4) * t232 + t61;
t25 = t98 * t219 - t211 * t221 + (t132 * t210 + t163 * t98) * t162;
t24 = -t129 * t98 - t130 * t211;
t14 = -pkin(8) * t232 + t16;
t10 = t99 * pkin(4) - t98 * t157 + t15;
t9 = -t114 * t247 - t57 * t121 - t129 * t88;
t8 = -t230 * t114 - t121 * t183 + t130 * t88;
t2 = t17 * t130 + t183 * t230;
t1 = -t17 * t129 - t130 * t18 + t183 * t247 + t230 * t57;
t3 = [0, 0, 0, 0.2e1 * t167 * t199, t214 * t244, t215, -t216, 0, -pkin(6) * t215 + t165 * t197, pkin(6) * t216 + t167 * t197, -t121 * t98 + t87 * t132, t121 * t99 - t132 * t88 + t180 * t87 + t249 * t98, t98 * t158, -t99 * t158, 0, t141 * t99 + t156 * t88 - t61 * t158 + (-qJD(1) * t180 - t249) * t207, t141 * t98 + t156 * t87 - t60 * t158 + (-t121 + t240) * t207, -t104 * t61 - t12 * t180 - t15 * t249 + t161 * t174 + t37 * t99 + t54 * t88, -t102 * t61 + t13 * t180 + t16 * t249 + t162 * t174 - t38 * t99 - t55 * t88, t16 * t104 + t15 * t102 + (-t12 * t132 - t37 * t98 - t54 * t87) * t162 + (-t13 * t132 - t38 * t98 - t55 * t87) * t161, t12 * t54 + t13 * t55 + t37 * t15 + t38 * t16 - t241 * t52 + t77 * t61, -t17 * t83 - t183 * t24, -t17 * t82 + t83 * t18 + t183 * t25 - t24 * t57, t24 * t114 - t17 * t180 - t183 * t99 - t83 * t88, -t25 * t114 + t18 * t180 - t57 * t99 - t82 * t88, t114 * t99 - t180 * t88, (t166 * t10 - t163 * t14) * t114 + t187 * t88 - t205 * t180 - t188 * t99 + t36 * t57 + t75 * t18 + t27 * t82 + t53 * t25 + (-t114 * t186 + t180 * t6) * qJD(5), -(t163 * t10 + t166 * t14) * t114 - t186 * t88 + t192 * t180 - t6 * t99 - t36 * t183 + t75 * t17 - t27 * t83 + t53 * t24 + (-t114 * t187 - t180 * t188) * qJD(5); 0, 0, 0, -t165 * t217, t214 * t169, 0, 0, 0, t169 * pkin(1) * t165, pkin(1) * t217, t225, t66, 0, t63, 0, t95 * t158 + (-t158 * t212 + t213 * t249) * pkin(2) + t175, t96 * t158 + (t121 * t213 - t158 * t200) * pkin(2) + t170, -t104 * t194 + t161 * t172 + t249 * t46 + t190, -t102 * t194 + t162 * t172 - t249 * t47 + t198, -t102 * t243 - t104 * t242 + t185, t146 * t189 + t152 * t191 + t52 * t155 + t194 * t77 - t37 * t46 - t38 * t47, t2, t1, t8, t9, t228, (t166 * t124 - t163 * t125) * t88 + t138 * t18 + t195 * t57 + (t163 * t182 - t166 * t181) * t114 + t177, -(t163 * t124 + t166 * t125) * t88 + t138 * t17 - t195 * t183 + (t163 * t181 + t166 * t182) * t114 + t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t225, t66, 0, t63, 0, t93 * t158 + t175, t92 * t158 + t170, t104 * t93 + t161 * t173 + t249 * t49 + t190, t102 * t93 + t162 * t173 - t249 * t50 + t198, -t50 * t104 - t49 * t102 + (-t102 * t161 + t104 * t162) * qJD(4) + t185, -t52 * pkin(3) + qJ(4) * t191 + qJD(4) * t189 - t37 * t49 - t38 * t50 - t77 * t93, t2, t1, t8, t9, t228, (t166 * t139 - t163 * t140) * t88 + t153 * t18 - t64 * t57 + (t163 * t179 - t166 * t178) * t114 + t177, -(t163 * t139 + t166 * t140) * t88 + t153 * t17 + t64 * t183 + (t163 * t178 + t166 * t179) * t114 + t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102 * t249 + t233, -t104 * t249 + t231, -t102 ^ 2 - t104 ^ 2, -t37 * t102 - t38 * t104 + t52, 0, 0, 0, 0, 0, t18 - t246, t114 * t245 + (t102 * t114 + t184) * t163 + t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t183 * t57, t183 ^ 2 - t57 ^ 2, t57 * t114 + t17, -t18 - t246, t88, t53 * t183 - t248 * t6 + t205, t248 * t188 + t53 * t57 - t192;];
tauc_reg = t3;
