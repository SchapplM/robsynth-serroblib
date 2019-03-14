% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRRRR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR4_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_inertiaDJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:00:01
% EndTime: 2019-03-09 01:00:10
% DurationCPUTime: 3.21s
% Computational Cost: add. (4536->341), mult. (12925->621), div. (0->0), fcn. (13337->14), ass. (0->182)
t144 = sin(qJ(5));
t141 = cos(pkin(7));
t145 = sin(qJ(4));
t150 = cos(qJ(4));
t139 = sin(pkin(7));
t146 = sin(qJ(3));
t214 = t139 * t146;
t187 = t150 * t214;
t107 = t141 * t145 + t187;
t151 = cos(qJ(3));
t100 = (-pkin(3) * t151 - pkin(10) * t146 - pkin(2)) * t139;
t213 = t139 * t151;
t194 = pkin(9) * t213;
t228 = pkin(2) * t146;
t99 = t194 + (pkin(10) + t228) * t141;
t179 = t150 * t100 - t145 * t99;
t195 = pkin(4) * t213;
t54 = -pkin(11) * t107 + t179 - t195;
t149 = cos(qJ(5));
t106 = -t141 * t150 + t145 * t214;
t225 = t145 * t100 + t150 * t99;
t58 = -pkin(11) * t106 + t225;
t55 = t149 * t58;
t168 = t144 * t54 + t55;
t204 = qJD(3) * t146;
t128 = t139 * t204;
t176 = pkin(4) * t128;
t102 = (pkin(3) * t146 - pkin(10) * t151) * t139 * qJD(3);
t203 = qJD(3) * t151;
t103 = -t141 * pkin(2) * t203 + pkin(9) * t128;
t48 = -t225 * qJD(4) + t150 * t102 + t145 * t103;
t184 = t139 * t203;
t87 = -qJD(4) * t106 + t150 * t184;
t38 = -pkin(11) * t87 + t176 + t48;
t201 = qJD(4) * t150;
t202 = qJD(4) * t145;
t47 = -t100 * t201 - t145 * t102 + t150 * t103 + t99 * t202;
t88 = qJD(4) * t187 + t141 * t202 + t145 * t184;
t41 = -pkin(11) * t88 - t47;
t180 = t144 * t41 - t149 * t38;
t12 = -t168 * qJD(5) - t180;
t10 = -pkin(5) * t128 - t12;
t148 = cos(qJ(6));
t143 = sin(qJ(6));
t197 = qJD(6) * t143;
t169 = -t144 * t58 + t149 * t54;
t30 = pkin(5) * t213 - t169;
t181 = -t10 * t148 + t30 * t197;
t138 = t148 ^ 2;
t206 = t143 ^ 2 - t138;
t177 = t206 * qJD(6);
t232 = qJD(4) + qJD(5);
t231 = 0.2e1 * t139;
t230 = pkin(10) + pkin(11);
t135 = qJD(6) * t148;
t229 = t10 * t143 + t30 * t135;
t227 = pkin(10) * t139;
t186 = qJD(4) * t230;
t174 = t150 * t186;
t175 = t145 * t186;
t123 = t230 * t145;
t124 = t230 * t150;
t90 = -t123 * t144 + t124 * t149;
t64 = t90 * qJD(5) - t144 * t175 + t149 * t174;
t89 = t123 * t149 + t124 * t144;
t226 = t89 * t135 + t64 * t143;
t224 = pkin(4) * qJD(5);
t115 = t144 * t150 + t145 * t149;
t114 = t144 * t145 - t149 * t150;
t85 = t232 * t114;
t222 = t115 * t85;
t74 = t149 * t106 + t107 * t144;
t45 = -t74 * qJD(5) - t144 * t88 + t149 * t87;
t75 = -t106 * t144 + t107 * t149;
t69 = t143 * t75 + t148 * t213;
t22 = -qJD(6) * t69 + t143 * t128 + t148 * t45;
t221 = t143 * t22;
t220 = t143 * t85;
t86 = t232 * t115;
t219 = t143 * t86;
t218 = t148 * t85;
t217 = t148 * t86;
t216 = t149 * t74;
t215 = t22 * t148;
t212 = t143 * t148;
t147 = sin(qJ(2));
t211 = t146 * t147;
t152 = cos(qJ(2));
t210 = t146 * t152;
t209 = t147 * t151;
t208 = t151 * t152;
t133 = -pkin(4) * t149 - pkin(5);
t199 = qJD(5) * t144;
t190 = pkin(4) * t199;
t207 = t133 * t135 + t143 * t190;
t140 = sin(pkin(6));
t205 = qJD(2) * t140;
t200 = qJD(4) * t151;
t198 = qJD(5) * t149;
t196 = -0.2e1 * pkin(3) * qJD(4);
t193 = pkin(4) * t202;
t192 = pkin(5) * t197;
t191 = pkin(5) * t135;
t189 = pkin(4) * t198;
t188 = t143 * t213;
t83 = t89 * t197;
t134 = -pkin(4) * t150 - pkin(3);
t136 = t139 ^ 2;
t185 = t136 * t203;
t183 = t147 * t205;
t182 = t143 * t135;
t178 = -0.4e1 * t115 * t212;
t173 = t136 * t183;
t172 = t139 * t183;
t171 = t146 * t185;
t31 = -pkin(12) * t213 + t168;
t98 = pkin(9) * t214 + (-pkin(2) * t151 - pkin(3)) * t141;
t77 = t106 * pkin(4) + t98;
t42 = pkin(5) * t74 - pkin(12) * t75 + t77;
t17 = t143 * t42 + t148 * t31;
t142 = cos(pkin(6));
t105 = -t139 * t140 * t152 + t142 * t141;
t160 = t141 * t210 + t209;
t82 = t140 * t160 + t142 * t214;
t67 = t105 * t150 - t145 * t82;
t68 = t105 * t145 + t150 * t82;
t44 = t144 * t67 + t149 * t68;
t161 = t141 * t208 - t211;
t81 = -t140 * t161 - t142 * t213;
t25 = t143 * t81 + t148 * t44;
t70 = t148 * t75 - t188;
t170 = -t143 * t70 - t148 * t69;
t79 = pkin(5) * t114 - pkin(12) * t115 + t134;
t57 = t143 * t79 + t148 * t90;
t132 = pkin(4) * t144 + pkin(12);
t167 = t114 * t132 - t115 * t133;
t166 = t133 * t197 - t148 * t190;
t66 = t142 * t184 + (t161 * qJD(3) + (-t141 * t211 + t208) * qJD(2)) * t140;
t153 = -t68 * qJD(4) - t66 * t145 + t150 * t172;
t37 = t67 * qJD(4) + t145 * t172 + t66 * t150;
t15 = t44 * qJD(5) + t144 * t37 - t149 * t153;
t43 = t144 * t68 - t149 * t67;
t5 = t43 * t135 + t143 * t15;
t6 = -t148 * t15 + t43 * t197;
t46 = t75 * qJD(5) + t144 * t87 + t149 * t88;
t27 = t74 * t135 + t143 * t46;
t165 = -t148 * t46 + t74 * t197;
t164 = t114 * t197 - t217;
t163 = -t115 * t135 + t220;
t162 = t115 * t197 + t218;
t11 = -t144 * t38 - t149 * t41 - t54 * t198 + t58 * t199;
t159 = t145 * t200 + t150 * t204;
t158 = t145 * t204 - t150 * t200;
t157 = pkin(5) * t86 + pkin(12) * t85 + t193;
t104 = (t141 * t228 + t194) * qJD(3);
t156 = pkin(12) * t128 - t11;
t73 = t88 * pkin(4) + t104;
t155 = -t132 * t86 - t133 * t85 + (-t114 * t149 + t115 * t144) * t224;
t154 = t46 * pkin(5) - t45 * pkin(12) + t73;
t127 = 0.2e1 * t182;
t116 = -0.2e1 * t171;
t113 = -0.2e1 * t177;
t112 = t115 ^ 2;
t72 = t114 * t135 + t219;
t65 = t142 * t128 + (t160 * qJD(3) + (t141 * t209 + t210) * qJD(2)) * t140;
t63 = t123 * t198 + t124 * t199 + t144 * t174 + t149 * t175;
t56 = -t143 * t90 + t148 * t79;
t53 = -t115 * t177 - t85 * t212;
t49 = qJD(6) * t178 + t206 * t85;
t24 = -t143 * t44 + t148 * t81;
t23 = -qJD(6) * t188 - t148 * t128 + t75 * t135 + t143 * t45;
t21 = -t57 * qJD(6) + t143 * t63 + t148 * t157;
t20 = -t79 * t135 - t143 * t157 + t148 * t63 + t90 * t197;
t19 = t70 * t135 + t221;
t16 = -t143 * t31 + t148 * t42;
t14 = -t144 * t153 - t149 * t37 - t67 * t198 + t68 * t199;
t7 = t170 * qJD(6) - t143 * t23 + t215;
t4 = -t25 * qJD(6) + t143 * t14 + t65 * t148;
t3 = -t81 * t135 + t148 * t14 - t65 * t143 + t44 * t197;
t2 = -t17 * qJD(6) - t143 * t156 + t148 * t154;
t1 = -t42 * t135 - t143 * t154 - t148 * t156 + t31 * t197;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t183, -t152 * t205, 0, 0, 0, 0, 0, t105 * t128 - t141 * t65 - t151 * t173, t105 * t184 - t141 * t66 + t146 * t173, 0, 0, 0, 0, 0, t65 * t106 + t128 * t67 - t153 * t213 + t81 * t88, t65 * t107 + t81 * t87 + (t151 * t37 - t204 * t68) * t139, 0, 0, 0, 0, 0, t81 * t46 + t65 * t74 + (t15 * t151 - t204 * t43) * t139, t81 * t45 + t65 * t75 + (-t14 * t151 - t204 * t44) * t139, 0, 0, 0, 0, 0, t15 * t69 + t23 * t43 + t24 * t46 + t4 * t74, t15 * t70 + t22 * t43 - t25 * t46 + t3 * t74; 0, 0, 0, 0, 0.2e1 * t171, 0.2e1 * (-t146 ^ 2 + t151 ^ 2) * t136 * qJD(3), 0.2e1 * t141 * t184, -0.2e1 * t141 * t128, 0, -0.2e1 * pkin(2) * t136 * t204 - 0.2e1 * t104 * t141, -0.2e1 * pkin(2) * t185 + 0.2e1 * t103 * t141, 0.2e1 * t107 * t87, -0.2e1 * t106 * t87 - 0.2e1 * t107 * t88 (t107 * t204 - t151 * t87) * t231 (-t106 * t204 + t151 * t88) * t231, t116, 0.2e1 * t104 * t106 + 0.2e1 * t98 * t88 + 0.2e1 * (-t48 * t151 + t179 * t204) * t139, 0.2e1 * t104 * t107 + 0.2e1 * t98 * t87 + 0.2e1 * (-t47 * t151 - t225 * t204) * t139, 0.2e1 * t75 * t45, -0.2e1 * t45 * t74 - 0.2e1 * t46 * t75 (-t151 * t45 + t204 * t75) * t231 (t151 * t46 - t204 * t74) * t231, t116, 0.2e1 * t77 * t46 + 0.2e1 * t73 * t74 + 0.2e1 * (-t12 * t151 + t169 * t204) * t139, 0.2e1 * t77 * t45 + 0.2e1 * t73 * t75 + 0.2e1 * (-t11 * t151 - t168 * t204) * t139, 0.2e1 * t70 * t22, -0.2e1 * t22 * t69 - 0.2e1 * t23 * t70, 0.2e1 * t22 * t74 + 0.2e1 * t46 * t70, -0.2e1 * t23 * t74 - 0.2e1 * t46 * t69, 0.2e1 * t74 * t46, 0.2e1 * t10 * t69 + 0.2e1 * t16 * t46 + 0.2e1 * t2 * t74 + 0.2e1 * t23 * t30, 0.2e1 * t1 * t74 + 0.2e1 * t10 * t70 - 0.2e1 * t17 * t46 + 0.2e1 * t22 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, -t66, 0, 0, 0, 0, 0, -t150 * t65 + t202 * t81, t145 * t65 + t201 * t81, 0, 0, 0, 0, 0, t114 * t65 + t81 * t86, t115 * t65 - t81 * t85, 0, 0, 0, 0, 0, t4 * t114 + t115 * t5 - t43 * t220 + t24 * t86, t3 * t114 - t115 * t6 - t43 * t218 - t25 * t86; 0, 0, 0, 0, 0, 0, t184, -t128, 0, -t104, t103, t107 * t201 + t145 * t87, -t145 * t88 + t87 * t150 + (-t106 * t150 - t107 * t145) * qJD(4), t158 * t139, t159 * t139, 0, -pkin(3) * t88 - t104 * t150 - t158 * t227 + t202 * t98, -pkin(3) * t87 + t104 * t145 - t159 * t227 + t201 * t98, t115 * t45 - t75 * t85, -t114 * t45 - t115 * t46 + t74 * t85 - t75 * t86 (t115 * t204 + t151 * t85) * t139 (-t114 * t204 + t151 * t86) * t139, 0, t74 * t193 + t114 * t73 + t134 * t46 + t77 * t86 + (t151 * t64 - t204 * t89) * t139, t75 * t193 + t115 * t73 + t134 * t45 - t77 * t85 + (-t151 * t63 - t204 * t90) * t139, -t70 * t218 + (-t197 * t70 + t215) * t115, -t170 * t85 + (-t221 - t148 * t23 + (t143 * t69 - t148 * t70) * qJD(6)) * t115, t22 * t114 - t115 * t165 - t74 * t218 + t70 * t86, -t23 * t114 - t115 * t27 + t74 * t220 - t69 * t86, t114 * t46 + t74 * t86, t114 * t2 + t229 * t115 + t16 * t86 + t21 * t74 - t30 * t220 + t23 * t89 + t46 * t56 + t64 * t69, t1 * t114 - t181 * t115 - t17 * t86 + t20 * t74 - t30 * t218 + t22 * t89 - t46 * t57 + t64 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t145 * t201, 0.2e1 * (-t145 ^ 2 + t150 ^ 2) * qJD(4), 0, 0, 0, t145 * t196, t150 * t196, -0.2e1 * t222, 0.2e1 * t114 * t85 - 0.2e1 * t115 * t86, 0, 0, 0, 0.2e1 * t114 * t193 + 0.2e1 * t134 * t86, 0.2e1 * t115 * t193 - 0.2e1 * t134 * t85, -0.2e1 * t112 * t182 - 0.2e1 * t138 * t222, 0.2e1 * t112 * t177 - t178 * t85, -0.2e1 * t114 * t162 + 0.2e1 * t115 * t217, 0.2e1 * t114 * t163 - 0.2e1 * t115 * t219, 0.2e1 * t114 * t86, 0.2e1 * t21 * t114 + 0.2e1 * t226 * t115 - 0.2e1 * t89 * t220 + 0.2e1 * t56 * t86, -0.2e1 * t89 * t218 + 0.2e1 * t20 * t114 - 0.2e1 * t57 * t86 + 0.2e1 * (t148 * t64 - t83) * t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, -t37, 0, 0, 0, 0, 0, -t15, t14, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, -t88, t128, t48, t47, 0, 0, t45, -t46, t128, t149 * t176 + (-t55 + (-t54 + t195) * t144) * qJD(5) - t180 (-t144 * t204 + t151 * t198) * t139 * pkin(4) + t11, t19, t7, t27, -t165, 0, t133 * t23 - t27 * t132 + (-t143 * t216 + t144 * t69) * t224 + t181, t133 * t22 + t165 * t132 + (t144 * t70 - t148 * t216) * t224 + t229; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t201, -t202, 0, -pkin(10) * t201, pkin(10) * t202, 0, 0, -t85, -t86, 0, -t64, t63, t53, t49, t72, -t164, 0, t83 + (-qJD(6) * t167 - t64) * t148 + t155 * t143, t148 * t155 + t167 * t197 + t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t190, -0.2e1 * t189, t127, t113, 0, 0, 0, 0.2e1 * t166, 0.2e1 * t207; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, t14, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t46, t128, t12, t11, t19, t7, t27, -t165, 0, -pkin(5) * t23 - pkin(12) * t27 + t181, -pkin(5) * t22 + pkin(12) * t165 + t229; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, -t86, 0, -t64, t63, t53, t49, t72, -t164, 0, t83 + (pkin(5) * t85 - pkin(12) * t86) * t143 + (-t64 + (-pkin(5) * t115 - pkin(12) * t114) * qJD(6)) * t148, pkin(5) * t162 + pkin(12) * t164 + t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t190, -t189, t127, t113, 0, 0, 0, t166 - t192, -t191 + t207; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127, t113, 0, 0, 0, -0.2e1 * t192, -0.2e1 * t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t23, t46, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t162, t163, t86, t21, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, -t197, 0, -t132 * t135 - t143 * t189, t132 * t197 - t148 * t189; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, -t197, 0, -pkin(12) * t135, pkin(12) * t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t8;