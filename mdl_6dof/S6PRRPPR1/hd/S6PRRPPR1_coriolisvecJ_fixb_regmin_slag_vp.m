% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:08
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRPPR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:06:04
% EndTime: 2021-01-16 02:06:18
% DurationCPUTime: 3.32s
% Computational Cost: add. (3708->347), mult. (9969->518), div. (0->0), fcn. (7853->12), ass. (0->185)
t159 = sin(qJ(2));
t154 = sin(pkin(6));
t212 = qJD(1) * t154;
t198 = t159 * t212;
t158 = sin(qJ(3));
t206 = t158 * qJD(3);
t253 = pkin(3) * t206 - t198;
t161 = cos(qJ(3));
t228 = cos(pkin(11));
t192 = t228 * t161;
t141 = qJD(2) * t192;
t153 = sin(pkin(11));
t210 = qJD(2) * t158;
t118 = t153 * t210 - t141;
t114 = qJD(6) + t118;
t243 = qJ(4) + pkin(8);
t193 = qJD(3) * t243;
t115 = t161 * qJD(4) - t158 * t193;
t169 = -t158 * qJD(4) - t161 * t193;
t171 = -t153 * t158 + t192;
t162 = cos(qJ(2));
t197 = t162 * t212;
t232 = -t115 * t228 - t153 * t169 + t171 * t197;
t131 = t153 * t161 + t158 * t228;
t120 = t131 * qJD(3);
t123 = t171 * qJD(3);
t252 = pkin(4) * t120 - qJ(5) * t123 - qJD(5) * t131 + t253;
t251 = qJD(6) - t114;
t152 = sin(pkin(12));
t155 = cos(pkin(12));
t157 = sin(qJ(6));
t160 = cos(qJ(6));
t132 = t152 * t160 + t155 * t157;
t229 = t114 * t132;
t121 = t131 * qJD(2);
t105 = qJD(3) * t152 + t121 * t155;
t106 = qJD(3) * t155 - t121 * t152;
t177 = -t105 * t160 - t106 * t157;
t250 = t114 * t177;
t249 = t160 * t106;
t240 = t152 * t232 + t155 * t252;
t239 = t152 * t252 - t155 * t232;
t237 = t115 * t153 - t131 * t197 - t169 * t228;
t111 = qJD(2) * t120;
t216 = t160 * t155;
t217 = t157 * t152;
t130 = -t216 + t217;
t230 = t114 * t130;
t248 = -t132 * t111 + t114 * t230;
t205 = qJD(2) * qJD(3);
t194 = t158 * t205;
t139 = t153 * t194;
t112 = qJD(3) * t141 - t139;
t22 = -qJD(6) * t177 + t112 * t132;
t117 = t118 ^ 2;
t247 = pkin(3) * t158;
t246 = t155 * pkin(9);
t134 = qJD(2) * pkin(8) + t198;
t156 = cos(pkin(6));
t211 = qJD(1) * t156;
t142 = t161 * t211;
t184 = qJD(4) + t197;
t65 = (-t134 * t158 + t142) * qJD(3) + (-qJ(4) * t206 + t161 * t184) * qJD(2);
t196 = t158 * t211;
t66 = (-t134 * t161 - t196) * qJD(3) + (-qJD(3) * t161 * qJ(4) - t158 * t184) * qJD(2);
t23 = t153 * t65 - t228 * t66;
t221 = t154 * t159;
t126 = t156 * t158 + t161 * t221;
t175 = t156 * t161 - t158 * t221;
t80 = t126 * t153 - t175 * t228;
t245 = t23 * t80;
t144 = pkin(3) * t153 + qJ(5);
t244 = pkin(9) + t144;
t242 = pkin(5) * t120 - t123 * t246 + t240;
t226 = t123 * t152;
t241 = pkin(9) * t226 - t239;
t24 = t153 * t66 + t228 * t65;
t20 = qJD(3) * qJD(5) + t24;
t116 = pkin(3) * t194 + qJD(2) * t198;
t37 = pkin(4) * t111 - qJ(5) * t112 - qJD(5) * t121 + t116;
t7 = t152 * t37 + t155 * t20;
t191 = qJ(4) * qJD(2) + t134;
t97 = t161 * t191 + t196;
t195 = t228 * t97;
t96 = -t158 * t191 + t142;
t92 = qJD(3) * pkin(3) + t96;
t45 = t153 * t92 + t195;
t40 = qJD(3) * qJ(5) + t45;
t148 = -pkin(3) * t161 - pkin(2);
t113 = qJD(2) * t148 + qJD(4) - t197;
t62 = t118 * pkin(4) - t121 * qJ(5) + t113;
t14 = t152 * t62 + t155 * t40;
t88 = t153 * t97;
t49 = t228 * t96 - t88;
t203 = pkin(3) * t210;
t75 = pkin(4) * t121 + qJ(5) * t118 + t203;
t19 = t152 * t75 + t155 * t49;
t238 = pkin(5) * t226 + t237;
t137 = t243 * t158;
t138 = t243 * t161;
t102 = -t137 * t153 + t138 * t228;
t84 = -pkin(4) * t171 - qJ(5) * t131 + t148;
t42 = t102 * t155 + t152 * t84;
t236 = qJD(2) * pkin(2);
t53 = t105 * t157 - t249;
t235 = t121 * t53;
t101 = t137 * t228 + t138 * t153;
t234 = t23 * t101;
t233 = t177 * t121;
t207 = qJD(6) * t160;
t231 = t106 * t207 + t112 * t216;
t227 = t118 * t152;
t225 = t131 * t152;
t224 = t131 * t155;
t222 = t152 * t112;
t220 = t154 * t162;
t164 = qJD(2) ^ 2;
t219 = t154 * t164;
t218 = t155 * t112;
t163 = qJD(3) ^ 2;
t215 = t163 * t158;
t214 = t163 * t161;
t213 = t158 ^ 2 - t161 ^ 2;
t209 = qJD(2) * t159;
t208 = qJD(6) * t131;
t202 = t159 * t219;
t6 = -t152 * t20 + t155 * t37;
t4 = pkin(5) * t111 - pkin(9) * t218 + t6;
t5 = -pkin(9) * t222 + t7;
t201 = -t157 * t5 + t160 * t4;
t200 = t154 * t209;
t199 = qJD(2) * t220;
t13 = -t152 * t40 + t155 * t62;
t18 = -t152 * t49 + t155 * t75;
t47 = t153 * t96 + t195;
t41 = -t102 * t152 + t155 * t84;
t190 = t158 * t199;
t189 = t161 * t199;
t188 = -t130 * t111 - t114 * t229;
t147 = -pkin(3) * t228 - pkin(4);
t187 = t157 * t4 + t160 * t5;
t8 = pkin(5) * t118 - pkin(9) * t105 + t13;
t9 = pkin(9) * t106 + t14;
t186 = t157 * t9 - t160 * t8;
t2 = t157 * t8 + t160 * t9;
t44 = t228 * t92 - t88;
t183 = -t13 * t152 + t14 * t155;
t27 = -pkin(5) * t171 - pkin(9) * t224 + t41;
t29 = -pkin(9) * t225 + t42;
t182 = -t157 * t29 + t160 * t27;
t181 = t157 * t27 + t160 * t29;
t81 = t126 * t228 + t153 * t175;
t63 = -t152 * t81 - t155 * t220;
t64 = -t152 * t220 + t155 * t81;
t180 = -t157 * t64 + t160 * t63;
t179 = t157 * t63 + t160 * t64;
t178 = t101 * t112 + t23 * t131;
t176 = -qJD(6) * t105 - t222;
t128 = t244 * t155;
t174 = pkin(5) * t121 + qJD(5) * t152 + qJD(6) * t128 + t118 * t246 + t18;
t127 = t244 * t152;
t173 = pkin(9) * t227 - qJD(5) * t155 + qJD(6) * t127 + t19;
t172 = t236 * qJD(2);
t39 = -qJD(3) * pkin(4) + qJD(5) - t44;
t170 = t123 * t39 + t178;
t168 = t126 * qJD(3);
t167 = -0.2e1 * qJD(3) * t236;
t21 = t157 * t176 + t231;
t166 = -t111 * t144 + t112 * t147 + (-qJD(5) + t39) * t118;
t165 = -t168 - t190;
t136 = -pkin(5) * t155 + t147;
t95 = qJD(3) * t175 + t189;
t79 = t130 * t131;
t78 = t132 * t131;
t72 = pkin(5) * t225 + t101;
t48 = t153 * t165 + t228 * t95;
t46 = t153 * t95 - t165 * t228;
t34 = t152 * t200 + t155 * t48;
t33 = -t152 * t48 + t155 * t200;
t32 = -pkin(5) * t227 + t47;
t31 = t123 * t132 + t207 * t224 - t208 * t217;
t30 = -t123 * t130 - t132 * t208;
t28 = -pkin(5) * t106 + t39;
t15 = pkin(5) * t222 + t23;
t1 = [0, 0, -t202, -t162 * t219, 0, 0, 0, 0, 0, -t161 * t202 + (-t168 - 0.2e1 * t190) * qJD(3), t158 * t202 + (-t95 - t189) * qJD(3), -t46 * qJD(3) + (-t111 * t162 + t118 * t209) * t154, -t48 * qJD(3) + (-t112 * t162 + t121 * t209) * t154, -t111 * t81 + t112 * t80 - t118 * t48 + t121 * t46, t245 + t24 * t81 - t44 * t46 + t45 * t48 + (t113 * t209 - t116 * t162) * t154, -t106 * t46 + t63 * t111 + t118 * t33 + t222 * t80, t105 * t46 - t64 * t111 - t118 * t34 + t218 * t80, t34 * t106 - t33 * t105 + (-t152 * t64 - t155 * t63) * t112, t13 * t33 + t14 * t34 + t39 * t46 + t6 * t63 + t64 * t7 + t245, 0, 0, 0, 0, 0, (-qJD(6) * t179 - t157 * t34 + t160 * t33) * t114 + t180 * t111 + t46 * t53 + t80 * t22, -(qJD(6) * t180 + t157 * t33 + t160 * t34) * t114 - t179 * t111 - t46 * t177 + t80 * t21; 0, 0, 0, 0, 0.2e1 * t161 * t194, -0.2e1 * t213 * t205, t214, -t215, 0, -pkin(8) * t214 + t158 * t167, pkin(8) * t215 + t161 * t167, -t118 * t198 + t148 * t111 + t113 * t120 - t116 * t171 + (t118 * t247 - t237) * qJD(3), -t121 * t198 + t148 * t112 + t113 * t123 + t116 * t131 + (t121 * t247 + t232) * qJD(3), -t102 * t111 + t118 * t232 - t45 * t120 + t121 * t237 - t44 * t123 + t171 * t24 + t178, t24 * t102 + t113 * t253 + t116 * t148 - t232 * t45 - t237 * t44 + t234, -t106 * t237 + t41 * t111 + t118 * t240 + t13 * t120 + t152 * t170 - t171 * t6, t105 * t237 - t42 * t111 - t118 * t239 - t14 * t120 + t155 * t170 + t171 * t7, -t240 * t105 + t239 * t106 + (-t112 * t41 - t123 * t13 - t131 * t6) * t155 + (-t112 * t42 - t123 * t14 - t131 * t7) * t152, t13 * t240 + t14 * t239 + t237 * t39 + t6 * t41 + t7 * t42 + t234, -t177 * t30 - t21 * t79, t177 * t31 - t21 * t78 + t22 * t79 - t30 * t53, -t111 * t79 + t114 * t30 - t120 * t177 - t171 * t21, -t111 * t78 - t114 * t31 - t120 * t53 + t171 * t22, -t111 * t171 + t114 * t120, t182 * t111 - t201 * t171 - t186 * t120 + t72 * t22 + t15 * t78 + t28 * t31 + t238 * t53 + (t157 * t241 + t160 * t242) * t114 + (-t114 * t181 + t171 * t2) * qJD(6), -t181 * t111 + t187 * t171 - t2 * t120 + t72 * t21 - t15 * t79 + t28 * t30 - t238 * t177 + (-t157 * t242 + t160 * t241) * t114 + (-t114 * t182 - t171 * t186) * qJD(6); 0, 0, 0, 0, -t158 * t164 * t161, t213 * t164, 0, 0, 0, t158 * t172, t161 * t172, qJD(3) * t47 - t113 * t121 - t118 * t203 - t23, qJD(3) * t49 + t113 * t118 - t121 * t203 - t24, (t45 - t47) * t121 + (-t44 + t49) * t118 + (-t111 * t153 - t112 * t228) * pkin(3), t44 * t47 - t45 * t49 + (-t113 * t210 + t153 * t24 - t228 * t23) * pkin(3), t106 * t47 - t18 * t118 - t13 * t121 + t152 * t166 - t23 * t155, -t47 * t105 + t19 * t118 + t14 * t121 + t23 * t152 + t155 * t166, -t19 * t106 + t18 * t105 + (qJD(5) * t106 - t118 * t13 + t7) * t155 + (qJD(5) * t105 - t118 * t14 - t6) * t152, -t13 * t18 - t14 * t19 + t23 * t147 - t39 * t47 + (-t6 * t152 + t7 * t155) * t144 + t183 * qJD(5), t21 * t132 + t177 * t230, -t21 * t130 - t132 * t22 + t177 * t229 + t230 * t53, t233 - t248, t188 + t235, -t114 * t121, (-t127 * t160 - t128 * t157) * t111 + t136 * t22 + t15 * t130 + t186 * t121 - t32 * t53 + t229 * t28 + (t157 * t173 - t160 * t174) * t114, -(-t127 * t157 + t128 * t160) * t111 + t136 * t21 + t15 * t132 + t2 * t121 + t32 * t177 - t230 * t28 + (t157 * t174 + t160 * t173) * t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t121 * qJD(3), -t139 + (t141 - t118) * qJD(3), -t121 ^ 2 - t117, t118 * t45 + t121 * t44 + t116, t106 * t121 + t111 * t155 - t117 * t152, -t105 * t121 - t111 * t152 - t117 * t155, (t105 * t152 + t106 * t155) * t118 + (-t152 ^ 2 - t155 ^ 2) * t112, t118 * t183 - t39 * t121 + t7 * t152 + t6 * t155, 0, 0, 0, 0, 0, t188 - t235, t233 + t248; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105 * t118 + t222, t106 * t118 + t218, -t105 ^ 2 - t106 ^ 2, t105 * t13 - t106 * t14 + t23, 0, 0, 0, 0, 0, t22 - t250, t114 * t249 + (-t105 * t114 + t176) * t157 + t231; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t177 * t53, t177 ^ 2 - t53 ^ 2, t53 * t114 + t21, -t22 - t250, t111, t28 * t177 - t2 * t251 + t201, t186 * t251 + t28 * t53 - t187;];
tauc_reg = t1;
