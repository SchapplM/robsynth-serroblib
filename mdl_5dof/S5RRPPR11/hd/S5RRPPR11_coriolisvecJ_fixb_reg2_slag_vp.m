% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPPR11_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR11_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:47:59
% EndTime: 2019-12-31 19:48:06
% DurationCPUTime: 2.06s
% Computational Cost: add. (2908->296), mult. (6799->427), div. (0->0), fcn. (4206->6), ass. (0->183)
t140 = sin(qJ(5));
t142 = cos(qJ(5));
t138 = cos(pkin(8));
t143 = cos(qJ(2));
t193 = qJD(1) * t143;
t178 = t138 * t193;
t137 = sin(pkin(8));
t186 = t137 * qJD(2);
t87 = -t178 - t186;
t179 = t137 * t193;
t192 = qJD(2) * t138;
t88 = -t179 + t192;
t39 = t140 * t88 - t142 * t87;
t231 = t39 ^ 2;
t42 = t140 * t87 + t142 * t88;
t230 = t42 ^ 2;
t220 = pkin(3) + pkin(6);
t141 = sin(qJ(2));
t90 = t137 * t142 + t138 * t140;
t150 = t90 * t141;
t78 = t90 * qJD(5);
t212 = qJD(1) * t150 + t78;
t194 = qJD(1) * t141;
t121 = qJD(5) + t194;
t229 = t121 * t39;
t185 = qJD(1) * qJD(2);
t228 = -0.2e1 * t185;
t123 = pkin(6) * t194;
t227 = qJD(3) + t123;
t134 = qJD(2) * qJ(3);
t226 = qJD(4) + t134;
t225 = -qJD(5) + t121;
t156 = t137 * t140 - t138 * t142;
t176 = t141 * t185;
t165 = t142 * t176;
t166 = t140 * t176;
t187 = qJD(5) * t142;
t188 = qJD(5) * t140;
t19 = -t137 * t165 - t138 * t166 - t87 * t187 + t188 * t88;
t224 = t156 * t19 - t212 * t42;
t149 = t137 * t166 - t138 * t165;
t20 = qJD(5) * t42 + t149;
t181 = t138 * t194;
t182 = t137 * t194;
t213 = t137 * t188 - t138 * t187 + t140 * t182 - t142 * t181;
t223 = t90 * t20 - t213 * t39;
t139 = -pkin(2) - qJ(4);
t190 = qJD(2) * t143;
t124 = pkin(6) * t193;
t125 = pkin(3) * t193;
t97 = t124 + t125;
t76 = t97 + t226;
t222 = t141 * (-t76 + t226) - t139 * t190;
t203 = t137 * t141;
t153 = pkin(4) * t143 - pkin(7) * t203;
t148 = t153 * qJD(2);
t120 = pkin(2) * t176;
t157 = -qJ(3) * t143 + qJ(4) * t141;
t189 = qJD(3) * t141;
t147 = qJD(2) * t157 - qJD(4) * t143 - t189;
t37 = qJD(1) * t147 + t120;
t175 = t143 * t185;
t119 = pkin(6) * t175;
t71 = t119 + (-qJD(4) + t125) * qJD(2);
t17 = -t137 * t37 + t138 * t71;
t12 = qJD(1) * t148 + t17;
t191 = qJD(2) * t141;
t168 = pkin(7) * t138 * t191;
t18 = t137 * t71 + t138 * t37;
t13 = qJD(1) * t168 + t18;
t160 = t12 * t140 + t13 * t142;
t174 = -qJ(3) * t141 - pkin(1);
t86 = t139 * t143 + t174;
t60 = t86 * qJD(1);
t196 = pkin(3) * t194 + t227;
t65 = t139 * qJD(2) + t196;
t25 = -t137 * t60 + t138 * t65;
t14 = pkin(4) * t194 - pkin(7) * t88 + t25;
t26 = t137 * t65 + t138 * t60;
t15 = pkin(7) * t87 + t26;
t5 = t14 * t142 - t140 * t15;
t1 = qJD(5) * t5 + t160;
t6 = t14 * t140 + t142 * t15;
t2 = -qJD(5) * t6 + t142 * t12 - t13 * t140;
t221 = t1 * t90 - t156 * t2 - t212 * t5 - t213 * t6;
t127 = pkin(2) * t194;
t72 = qJD(1) * t157 + t127;
t35 = -t137 * t72 + t138 * t97;
t22 = qJD(1) * t153 + t35;
t36 = t137 * t97 + t138 * t72;
t29 = pkin(7) * t181 + t36;
t214 = -pkin(7) + t139;
t100 = t214 * t138;
t99 = t214 * t137;
t47 = t100 * t142 - t140 * t99;
t219 = -qJD(4) * t90 + qJD(5) * t47 - t140 * t22 - t142 * t29;
t48 = t100 * t140 + t142 * t99;
t218 = qJD(4) * t156 - qJD(5) * t48 + t140 * t29 - t142 * t22;
t217 = t42 * t39;
t126 = pkin(2) * t191;
t50 = t126 + t147;
t98 = t220 * t190;
t28 = t137 * t98 + t138 * t50;
t112 = t220 * t141;
t44 = t137 * t112 + t138 * t86;
t211 = qJD(2) * pkin(2);
t133 = qJD(2) * qJD(3);
t96 = t220 * t191;
t70 = -qJD(1) * t96 + t133;
t210 = t137 * t70;
t209 = t137 * t88;
t208 = t138 * t87;
t207 = t138 * t88;
t206 = t143 * t25;
t205 = t143 * t26;
t135 = t141 ^ 2;
t145 = qJD(1) ^ 2;
t204 = t135 * t145;
t202 = t138 * t141;
t201 = t138 * t143;
t200 = t143 * t145;
t144 = qJD(2) ^ 2;
t199 = t144 * t141;
t198 = t144 * t143;
t183 = -pkin(4) * t138 - pkin(3);
t197 = -t183 * t194 + t227;
t113 = t220 * t143;
t136 = t143 ^ 2;
t195 = t135 - t136;
t106 = -pkin(2) * t143 + t174;
t85 = qJD(1) * t106;
t184 = t138 * t204;
t177 = t213 * t121;
t27 = -t137 * t50 + t138 * t98;
t173 = pkin(1) * t228;
t172 = qJD(3) - t211;
t171 = t87 + t186;
t170 = -t88 + t192;
t169 = qJD(1) * (0.2e1 * t135 - t136);
t167 = -t212 * t121 - t156 * t175;
t164 = t141 * t175;
t162 = -t137 * t204 + t138 * t175;
t159 = t137 * t18 + t138 * t17;
t158 = -t137 * t25 + t138 * t26;
t93 = t138 * t112;
t30 = pkin(4) * t141 + t93 + (pkin(7) * t143 - t86) * t137;
t33 = -pkin(7) * t201 + t44;
t9 = -t140 * t33 + t142 * t30;
t10 = t140 * t30 + t142 * t33;
t155 = qJD(1) * t171;
t154 = -0.2e1 * qJD(2) * t85;
t151 = -qJ(3) * t190 - t189;
t61 = qJD(1) * t151 + t120;
t74 = t126 + t151;
t152 = pkin(6) * t144 + qJD(1) * t74 + t61;
t63 = (-pkin(6) + t183) * t191;
t67 = t156 * t143;
t49 = qJD(1) * t63 + t133;
t104 = pkin(6) * t176 - t133;
t105 = t123 + t172;
t109 = -t124 - t134;
t146 = -t104 * t143 + (t105 * t143 + (t109 + t124) * t141) * qJD(2);
t132 = t138 ^ 2;
t131 = t137 ^ 2;
t122 = pkin(4) * t137 + qJ(3);
t118 = t141 * t200;
t111 = -0.2e1 * t164;
t110 = 0.2e1 * t164;
t108 = t195 * t145;
t94 = -qJ(3) * t193 + t127;
t84 = t195 * t228;
t77 = pkin(4) * t201 + t113;
t68 = t90 * t143;
t64 = t85 * t194;
t45 = -pkin(4) * t87 + t76;
t43 = -t137 * t86 + t93;
t32 = t143 * t78 - t156 * t191;
t31 = qJD(2) * t150 + qJD(5) * t67;
t21 = t168 + t28;
t16 = t148 + t27;
t4 = -qJD(5) * t10 - t140 * t21 + t142 * t16;
t3 = qJD(5) * t9 + t140 * t16 + t142 * t21;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, t84, t198, t111, -t199, 0, -pkin(6) * t198 + t141 * t173, pkin(6) * t199 + t143 * t173, 0, 0, 0, -t198, t199, t110, t84, t111, t146, t141 * t154 + t143 * t152, -t141 * t152 + t143 * t154, pkin(6) * t146 + t106 * t61 + t74 * t85, (-t131 * t193 + t209) * t191, (t207 + (t87 - 0.2e1 * t178) * t137) * t191, (t137 * t169 + t143 * t88) * qJD(2), (-t132 * t193 + t208) * t191, (t138 * t169 + t143 * t87) * qJD(2), t110, t70 * t201 + t87 * t96 + (qJD(1) * t27 + t17) * t141 + (-t76 * t202 + t206 + (-t113 * t202 + t143 * t43) * qJD(1)) * qJD(2), -t143 * t210 - t88 * t96 + (-qJD(1) * t28 - t18) * t141 + (t76 * t203 - t205 + (t113 * t203 - t143 * t44) * qJD(1)) * qJD(2), -t27 * t88 + t28 * t87 + (t137 * t17 - t138 * t18) * t143 + ((-t137 * t43 + t138 * t44) * qJD(1) + t158) * t191, t113 * t70 + t17 * t43 + t18 * t44 + t25 * t27 + t26 * t28 - t76 * t96, t19 * t68 + t31 * t42, -t19 * t67 + t20 * t68 - t31 * t39 + t32 * t42, t121 * t31 - t141 * t19 + (-qJD(1) * t68 + t42) * t190, -t20 * t67 - t32 * t39, t121 * t32 - t141 * t20 + (qJD(1) * t67 - t39) * t190, (t121 + t194) * t190, t121 * t4 + t141 * t2 + t20 * t77 - t32 * t45 + t39 * t63 - t49 * t67 + (qJD(1) * t9 + t5) * t190, -t1 * t141 - t121 * t3 - t19 * t77 + t31 * t45 + t42 * t63 - t49 * t68 + (-qJD(1) * t10 - t6) * t190, t1 * t67 - t10 * t20 + t19 * t9 + t2 * t68 - t3 * t39 - t31 * t5 + t32 * t6 - t4 * t42, t1 * t10 + t2 * t9 + t3 * t6 + t4 * t5 + t45 * t63 + t49 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t118, t108, 0, t118, 0, 0, t145 * pkin(1) * t141, pkin(1) * t200, 0, 0, 0, 0, 0, -t118, t108, t118, ((-t109 - t134) * t141 + (-t105 + t172) * t143) * qJD(1), -t193 * t94 + t64, 0.2e1 * t133 + (t141 * t94 + t143 * t85) * qJD(1), -qJ(3) * t104 - qJD(3) * t109 - t85 * t94 + (-t109 * t141 + (-t105 - t211) * t143) * qJD(1) * pkin(6), t170 * t182, (-t137 * t87 - t207 + (-t131 + t132) * qJD(2)) * t194, -t193 * t88 + t162, -t155 * t202, -t143 * t155 - t184, -t118, t210 - t196 * t87 + (-t222 * t138 - t141 * t35 - t206) * qJD(1), t138 * t70 + t196 * t88 + (t222 * t137 + t141 * t36 + t205) * qJD(1), t35 * t88 - t36 * t87 + (qJD(4) * t88 - t194 * t26 - t17) * t138 + (-qJD(4) * t87 + t194 * t25 - t18) * t137, qJ(3) * t70 - t25 * t35 - t26 * t36 + t196 * t76 + t159 * t139 + (-t137 * t26 - t138 * t25) * qJD(4), t224, t156 * t20 + t19 * t90 + t212 * t39 + t213 * t42, -t193 * t42 + t167, t223, t177 + (-qJD(2) * t90 + t39) * t193, -t121 * t193, t122 * t20 + t49 * t90 - t213 * t45 + t197 * t39 + t218 * t121 + (qJD(2) * t47 - t5) * t193, -t122 * t19 - t49 * t156 - t212 * t45 + t197 * t42 - t219 * t121 + (-qJD(2) * t48 + t6) * t193, t19 * t47 - t20 * t48 - t218 * t42 - t219 * t39 - t221, t1 * t48 + t122 * t49 + t197 * t45 + t2 * t47 + t218 * t5 + t219 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, -t144 - t204, qJD(2) * t109 + t119 + t64, 0, 0, 0, 0, 0, 0, qJD(2) * t87 + t162, -t184 + (-t88 - t179) * qJD(2), (t208 + t209) * t194, -qJD(2) * t76 + t158 * t194 + t159, 0, 0, 0, 0, 0, 0, -qJD(2) * t39 + t167, t177 + (-t193 * t90 - t42) * qJD(2), -t223 - t224, -qJD(2) * t45 + t221; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t170 * t194, t171 * t194, -t87 ^ 2 - t88 ^ 2, t25 * t88 - t26 * t87 + t70, 0, 0, 0, 0, 0, 0, t121 * t42 + t20, -t19 - t229, -t230 - t231, t39 * t6 + t42 * t5 + t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t217, t230 - t231, -t19 + t229, -t217, t225 * t42 - t149, t175, t121 * t6 - t42 * t45 + t2, t225 * t5 + t39 * t45 - t160, 0, 0;];
tauc_reg = t7;
