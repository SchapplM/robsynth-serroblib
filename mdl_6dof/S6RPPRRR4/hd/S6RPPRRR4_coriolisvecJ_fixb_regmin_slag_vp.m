% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRRR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:26:52
% EndTime: 2019-03-09 02:26:57
% DurationCPUTime: 2.32s
% Computational Cost: add. (2396->326), mult. (4984->476), div. (0->0), fcn. (3278->8), ass. (0->177)
t123 = sin(qJ(6));
t126 = cos(qJ(6));
t124 = sin(qJ(5));
t127 = cos(qJ(5));
t190 = qJD(4) * t127;
t125 = sin(qJ(4));
t195 = qJD(1) * t125;
t84 = t124 * t195 + t190;
t184 = t124 * qJD(4);
t85 = t127 * t195 - t184;
t144 = t123 * t84 - t126 * t85;
t38 = t123 * t85 + t126 * t84;
t243 = t144 * t38;
t242 = t144 ^ 2 - t38 ^ 2;
t128 = cos(qJ(4));
t129 = -pkin(1) - pkin(2);
t104 = t129 * qJD(1) + qJD(2);
t120 = sin(pkin(10));
t121 = cos(pkin(10));
t196 = qJ(2) * qJD(1);
t72 = t120 * t104 + t121 * t196;
t61 = -qJD(1) * pkin(7) + t72;
t238 = qJD(3) * t128 - t125 * t61;
t32 = -qJD(4) * pkin(4) - t238;
t23 = -pkin(5) * t84 + t32;
t117 = t125 * qJD(3);
t41 = t128 * t61 + t117;
t33 = qJD(4) * pkin(8) + t41;
t148 = pkin(4) * t128 + pkin(8) * t125;
t71 = t104 * t121 - t120 * t196;
t60 = qJD(1) * pkin(3) - t71;
t34 = qJD(1) * t148 + t60;
t15 = t124 * t34 + t127 * t33;
t13 = pkin(9) * t84 + t15;
t186 = qJD(6) * t123;
t9 = t13 * t186;
t241 = -t23 * t38 + t9;
t185 = qJD(6) * t126;
t182 = qJD(1) * qJD(4);
t163 = t128 * t182;
t188 = qJD(5) * t124;
t168 = t125 * t188;
t181 = qJD(4) * qJD(5);
t55 = qJD(1) * t168 + (-t163 + t181) * t127;
t170 = t128 * t184;
t187 = qJD(5) * t127;
t136 = t125 * t187 + t170;
t56 = qJD(1) * t136 - t124 * t181;
t10 = t123 * t56 + t126 * t55 + t84 * t185 + t186 * t85;
t194 = qJD(1) * t128;
t105 = qJD(5) + t194;
t103 = qJD(6) + t105;
t240 = -t103 * t38 + t10;
t130 = qJD(4) ^ 2;
t131 = qJD(1) ^ 2;
t237 = t120 * (t130 + t131);
t118 = t125 ^ 2;
t236 = (qJD(1) * t118 - t105 * t128) * t124;
t235 = qJD(5) + qJD(6);
t183 = qJD(1) * qJD(2);
t162 = t121 * t183;
t25 = qJD(4) * t238 + t128 * t162;
t147 = -pkin(4) * t125 + pkin(8) * t128;
t77 = qJD(2) * t120 + qJD(4) * t147;
t62 = t77 * qJD(1);
t54 = t127 * t62;
t132 = -qJD(5) * t15 - t124 * t25 + t54;
t161 = t125 * t182;
t4 = -pkin(5) * t161 - pkin(9) * t55 + t132;
t180 = -t124 * t62 - t127 * t25 - t34 * t187;
t139 = -t188 * t33 - t180;
t5 = pkin(9) * t56 + t139;
t175 = -t123 * t5 + t126 * t4;
t215 = t126 * t13;
t14 = -t124 * t33 + t127 * t34;
t12 = pkin(9) * t85 + t14;
t8 = pkin(5) * t105 + t12;
t2 = t123 * t8 + t215;
t234 = -qJD(6) * t2 - t23 * t144 + t175;
t11 = qJD(6) * t144 + t123 * t55 - t126 * t56;
t233 = t103 * t144 - t11;
t232 = pkin(8) + pkin(9);
t90 = t147 * qJD(1);
t231 = t124 * t90 + t127 * t238;
t203 = t126 * t127;
t88 = t123 * t124 - t203;
t137 = t88 * t128;
t230 = -qJD(1) * t137 - t235 * t88;
t89 = t123 * t127 + t124 * t126;
t43 = t235 * t89;
t229 = t89 * t194 + t43;
t169 = t125 * t190;
t202 = t127 * t128;
t205 = t124 * t128;
t80 = -t120 * t205 - t121 * t127;
t228 = -qJD(5) * t80 + t120 * t169 + (t120 * t124 + t121 * t202) * qJD(1);
t171 = t125 * t184;
t81 = t120 * t202 - t121 * t124;
t227 = -qJD(5) * t81 + t120 * t171 - (t120 * t127 - t121 * t205) * qJD(1);
t199 = t121 * qJ(2) + t120 * t129;
t87 = -pkin(7) + t199;
t226 = t127 * t77 + t87 * t171;
t157 = -t120 * qJ(2) + t121 * t129;
t86 = pkin(3) - t157;
t63 = t148 + t86;
t73 = t87 * t202;
t225 = t124 * t63 + t73;
t224 = t10 * t128;
t18 = qJD(4) * t137 + t125 * t43;
t223 = t103 * t18;
t189 = qJD(4) * t128;
t165 = t127 * t189;
t206 = t124 * t125;
t176 = t123 * t206;
t204 = t125 * t127;
t19 = -qJD(6) * t176 + (t204 * t235 + t170) * t126 + (t165 - t168) * t123;
t222 = t103 * t19;
t221 = t105 * t84;
t220 = t105 * t85;
t219 = t11 * t128;
t218 = t124 * t55;
t217 = t125 * t84;
t216 = t125 * t85;
t214 = t127 * t32;
t213 = t128 * t55;
t212 = t128 * t56;
t26 = qJD(4) * t117 + t125 * t162 + t61 * t189;
t211 = t26 * t124;
t210 = t26 * t127;
t74 = t89 * t125;
t209 = qJD(1) * t74;
t75 = t125 * t203 - t176;
t208 = qJD(1) * t75;
t207 = t105 * t127;
t201 = t130 * t125;
t200 = t130 * t128;
t198 = -t128 ^ 2 + t118;
t193 = qJD(2) * t121;
t191 = qJD(4) * t125;
t173 = t128 * t193;
t179 = t124 * t77 + t127 * t173 + t63 * t187;
t178 = 0.2e1 * t183;
t177 = t105 * t202;
t174 = qJD(5) * t232;
t172 = t120 * t189;
t167 = t105 * t187;
t166 = t124 * t194;
t164 = qJD(6) * t8 + t5;
t160 = t105 * t87 + t33;
t158 = -t124 * t238 + t127 * t90;
t156 = t120 * t178;
t155 = 0.2e1 * t161;
t154 = t125 * t167;
t153 = -t41 + (t166 + t188) * pkin(5);
t141 = -pkin(5) * t125 + pkin(9) * t202;
t98 = t232 * t127;
t152 = qJD(1) * t141 + qJD(6) * t98 + t127 * t174 + t158;
t97 = t232 * t124;
t151 = pkin(9) * t166 + qJD(6) * t97 + t124 * t174 + t231;
t150 = -qJD(6) * t80 + t228;
t149 = qJD(6) * t81 - t227;
t146 = t71 * t120 - t72 * t121;
t58 = t127 * t63;
t21 = pkin(9) * t204 + t58 + (-t124 * t87 + pkin(5)) * t128;
t22 = pkin(9) * t206 + t225;
t145 = t123 * t21 + t126 * t22;
t140 = t127 * t118 * t182 + t105 * t168;
t138 = -t130 * t87 + t156;
t135 = qJD(4) * (-qJD(1) * t86 - t193 - t60);
t113 = -pkin(5) * t127 - pkin(4);
t59 = (-pkin(5) * t124 + t87) * t125;
t27 = -pkin(5) * t136 + t125 * t193 + t189 * t87;
t16 = -pkin(5) * t56 + t26;
t7 = (-t128 * t188 - t169) * t87 + t136 * pkin(9) + t179;
t6 = -t124 * t173 + t141 * qJD(4) + (-t73 + (-pkin(9) * t125 - t63) * t124) * qJD(5) + t226;
t1 = -t123 * t13 + t126 * t8;
t3 = [0, 0, 0, 0, t178, qJ(2) * t178, t156, 0.2e1 * t162 ((-t120 * t157 + t121 * t199) * qJD(1) - t146) * qJD(2), t128 * t155, -0.2e1 * t198 * t182, -t200, t201, 0, t125 * t135 + t128 * t138, -t125 * t138 + t128 * t135, t85 * t165 + (-t127 * t55 - t188 * t85) * t125 (-t124 * t85 - t127 * t84) * t189 + (t218 - t127 * t56 + (t124 * t84 - t127 * t85) * qJD(5)) * t125, t213 + (-t177 + t216) * qJD(4) + t140, t154 + t212 + (-t217 - t236) * qJD(4) (-t105 - t194) * t191 (-t188 * t63 + t226) * t105 + (-qJD(4) * t87 * t84 + t54 - t160 * t187 + (-t32 * qJD(4) - qJD(5) * t34 - t105 * t193 - t25) * t124) * t128 + (-t84 * t193 - t32 * t187 - t211 - t87 * t56 + (-(-t87 * t205 + t58) * qJD(1) - t14) * qJD(4)) * t125, -t179 * t105 + (t160 * t188 + (-t85 * t87 - t214) * qJD(4) + t180) * t128 + (-t85 * t193 + t32 * t188 - t210 + t87 * t55 + (t225 * qJD(1) + t87 * t207 + t15) * qJD(4)) * t125, -t10 * t75 + t144 * t18, t10 * t74 + t11 * t75 + t144 * t19 + t18 * t38, t224 + t223 + (-t144 + t208) * t191, t222 - t219 + (-t38 - t209) * t191 (-t103 - t194) * t191 (-t123 * t7 + t126 * t6) * t103 + t175 * t128 - t27 * t38 + t59 * t11 - t16 * t74 - t23 * t19 + (-t103 * t145 - t128 * t2) * qJD(6) + (-(-t123 * t22 + t126 * t21) * qJD(1) - t1) * t191, t59 * t10 + t9 * t128 - t16 * t75 + t23 * t18 + t27 * t144 + (-(-qJD(6) * t22 + t6) * t103 - t4 * t128) * t123 + (-(qJD(6) * t21 + t7) * t103 - t164 * t128) * t126 + (qJD(1) * t145 + t2) * t191; 0, 0, 0, 0, -t131, -t131 * qJ(2), -t120 * t131, -t121 * t131, t146 * qJD(1), 0, 0, 0, 0, 0, t121 * t155 - t128 * t237, 0.2e1 * t121 * t163 + t125 * t237, 0, 0, 0, 0, 0, -t84 * t172 + t227 * t105 + (-t120 * t56 + (-qJD(4) * t80 + t121 * t84) * qJD(1)) * t125, -t85 * t172 + t228 * t105 + (t120 * t55 + (qJD(4) * t81 + t121 * t85) * qJD(1)) * t125, 0, 0, 0, 0, 0, -t38 * t172 + (t123 * t150 - t126 * t149) * t103 + (t120 * t11 + (-(-t123 * t81 + t126 * t80) * qJD(4) + t121 * t38) * qJD(1)) * t125, t144 * t172 + (t123 * t149 + t126 * t150) * t103 + (t120 * t10 + ((t123 * t80 + t126 * t81) * qJD(4) - t121 * t144) * qJD(1)) * t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t201, -t200, 0, 0, 0, 0, 0, -t154 + t212 + (-t217 + t236) * qJD(4), -t213 + (-t177 - t216) * qJD(4) + t140, 0, 0, 0, 0, 0, -t222 - t219 + (-t38 + t209) * t191, -t224 + t223 + (t144 + t208) * t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t125 * t131 * t128, t198 * t131, 0, 0, 0, qJD(4) * t41 + t195 * t60 - t26 (t60 - t193) * t194, -t85 * t207 + t218 (t55 + t221) * t127 + (t56 + t220) * t124, t167 + (t177 + (-t85 - t184) * t125) * qJD(1), -t105 * t188 + (-t105 * t205 + (t84 - t190) * t125) * qJD(1), t105 * t195, pkin(4) * t56 - t210 - t158 * t105 + t41 * t84 + (-pkin(8) * t207 + t124 * t32) * qJD(5) + (t14 * t125 + (pkin(8) * t191 + t128 * t32) * t124) * qJD(1), -pkin(4) * t55 + t211 + t231 * t105 + t41 * t85 + (pkin(8) * t105 * t124 + t214) * qJD(5) + (t32 * t202 + (pkin(8) * t190 - t15) * t125) * qJD(1), t10 * t89 + t144 * t230, -t10 * t88 - t11 * t89 - t144 * t229 + t230 * t38, t230 * t103 + (-qJD(4) * t89 + t144) * t195, -t229 * t103 + (qJD(4) * t88 + t38) * t195, t103 * t195, t113 * t11 + t16 * t88 - t153 * t38 + t229 * t23 + (t123 * t151 - t126 * t152) * t103 + (-(-t123 * t98 - t126 * t97) * qJD(4) + t1) * t195, t113 * t10 + t16 * t89 + t153 * t144 + t230 * t23 + (t123 * t152 + t126 * t151) * t103 + ((-t123 * t97 + t126 * t98) * qJD(4) - t2) * t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85 * t84, -t84 ^ 2 + t85 ^ 2, t55 - t221, t56 - t220, -t161, t105 * t15 + t32 * t85 + t132, t105 * t14 - t32 * t84 - t139, -t243, t242, t240, t233, -t161 -(-t12 * t123 - t215) * t103 + (-t103 * t186 - t126 * t161 - t38 * t85) * pkin(5) + t234 (-t103 * t13 - t4) * t123 + (t103 * t12 - t164) * t126 + (-t103 * t185 + t123 * t161 + t144 * t85) * pkin(5) + t241; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t243, t242, t240, t233, -t161, t103 * t2 + t234, t1 * t103 - t123 * t4 - t126 * t164 + t241;];
tauc_reg  = t3;
