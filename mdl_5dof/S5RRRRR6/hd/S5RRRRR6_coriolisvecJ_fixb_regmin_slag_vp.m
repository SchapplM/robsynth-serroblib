% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% tauc_reg [5x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 12:09
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:08:22
% EndTime: 2022-01-20 12:08:29
% DurationCPUTime: 1.74s
% Computational Cost: add. (2668->224), mult. (4556->304), div. (0->0), fcn. (3270->8), ass. (0->175)
t139 = sin(qJ(5));
t195 = qJD(5) * t139;
t140 = sin(qJ(4));
t144 = cos(qJ(4));
t145 = cos(qJ(3));
t136 = qJD(1) + qJD(2);
t142 = sin(qJ(2));
t213 = pkin(1) * qJD(1);
t189 = t142 * t213;
t235 = pkin(7) + pkin(8);
t177 = t235 * t136 + t189;
t82 = t177 * t145;
t75 = t144 * t82;
t141 = sin(qJ(3));
t81 = t177 * t141;
t76 = qJD(3) * pkin(3) - t81;
t163 = -t140 * t76 - t75;
t203 = t144 * t145;
t185 = t136 * t203;
t208 = t140 * t141;
t186 = t136 * t208;
t90 = -t185 + t186;
t228 = t90 * pkin(9);
t28 = -t163 - t228;
t143 = cos(qJ(5));
t83 = t143 * t90;
t106 = t140 * t145 + t144 * t141;
t92 = t106 * t136;
t45 = -t139 * t92 - t83;
t128 = -t145 * pkin(3) - pkin(2);
t146 = cos(qJ(2));
t188 = t146 * t213;
t93 = t128 * t136 - t188;
t50 = t90 * pkin(4) + t93;
t237 = t28 * t195 - t50 * t45;
t192 = -qJD(4) - qJD(5);
t131 = qJD(3) - t192;
t162 = qJD(3) * t177;
t212 = pkin(1) * qJD(2);
t187 = qJD(1) * t212;
t169 = t146 * t187;
t55 = -t141 * t162 + t145 * t169;
t56 = -t141 * t169 - t145 * t162;
t175 = -t140 * t55 + t144 * t56;
t151 = t163 * qJD(4) + t175;
t135 = qJD(3) + qJD(4);
t199 = qJD(3) * t145;
t181 = t136 * t199;
t47 = qJD(4) * t185 - t135 * t186 + t144 * t181;
t9 = -t47 * pkin(9) + t151;
t236 = (-t28 * t131 - t9) * t139 + t237;
t164 = t139 * t90 - t143 * t92;
t223 = t164 * t45;
t6 = t164 ^ 2 - t45 ^ 2;
t64 = t135 * t106;
t48 = t64 * t136;
t10 = -qJD(5) * t83 - t139 * t48 + t143 * t47 - t92 * t195;
t4 = -t45 * t131 + t10;
t197 = qJD(4) * t140;
t174 = t140 * t56 - t82 * t197;
t232 = (qJD(4) * t76 + t55) * t144;
t8 = -t48 * pkin(9) + t174 + t232;
t159 = -t139 * t8 + t143 * t9 + t164 * t50;
t150 = t164 * qJD(5) - t139 * t47 - t143 * t48;
t5 = -t131 * t164 + t150;
t105 = -t203 + t208;
t183 = qJD(3) * t235;
t107 = t141 * t183;
t108 = t145 * t183;
t120 = t235 * t141;
t133 = t145 * pkin(8);
t121 = t145 * pkin(7) + t133;
t196 = qJD(4) * t144;
t234 = -t105 * t188 + t144 * t107 + t140 * t108 + t120 * t196 + t121 * t197;
t73 = t140 * t82;
t173 = t144 * t76 - t73;
t87 = t92 * pkin(9);
t27 = t173 - t87;
t160 = t140 * t120 - t144 * t121;
t233 = t160 * qJD(4) + t106 * t188 + t140 * t107 - t144 * t108;
t231 = qJD(5) - t131;
t63 = t135 * t105;
t229 = t63 * pkin(9);
t227 = t92 * pkin(4);
t226 = pkin(3) * t131;
t225 = t106 * pkin(9);
t224 = t146 * pkin(1);
t221 = t92 * t90;
t220 = t93 * t92;
t125 = t142 * pkin(1) + pkin(7);
t219 = -pkin(8) - t125;
t60 = t143 * t105 + t139 * t106;
t19 = -t60 * qJD(5) - t139 * t64 - t143 * t63;
t123 = t142 * t187;
t200 = qJD(3) * t141;
t182 = t136 * t200;
t95 = pkin(3) * t182 + t123;
t33 = t48 * pkin(4) + t95;
t61 = -t139 * t105 + t143 * t106;
t218 = t50 * t19 + t33 * t61;
t20 = t61 * qJD(5) - t139 * t63 + t143 * t64;
t217 = t50 * t20 + t33 * t60;
t216 = t95 * t106 - t93 * t63;
t215 = t95 * t105 + t93 * t64;
t214 = -t144 * t81 - t73;
t211 = t143 * t28;
t112 = -t136 * pkin(2) - t188;
t210 = t112 * t199 + t141 * t123;
t209 = t136 * t141;
t207 = t140 * t143;
t205 = t142 * t145;
t147 = qJD(3) ^ 2;
t202 = t147 * t141;
t132 = t147 * t145;
t201 = t141 ^ 2 - t145 ^ 2;
t198 = qJD(3) * t146;
t194 = -qJD(1) - t136;
t193 = -qJD(2) + t136;
t191 = pkin(3) * t209;
t190 = t146 * t212;
t130 = t142 * t212;
t129 = pkin(3) * t200;
t180 = t141 * t198;
t179 = -pkin(3) * t135 - t76;
t54 = t64 * pkin(4) + t129;
t25 = t135 * pkin(4) + t27;
t178 = -pkin(4) * t131 - t25;
t172 = t140 * t81 - t75;
t171 = qJD(3) * t219;
t168 = t54 - t189;
t62 = t64 * pkin(9);
t167 = -qJD(5) * (-t144 * t120 - t140 * t121 - t225) + t62 + t234;
t101 = t105 * pkin(9);
t166 = qJD(5) * (-t101 - t160) - t229 - t233;
t165 = -t139 * t25 - t211;
t102 = t219 * t141;
t103 = t145 * t125 + t133;
t161 = -t140 * t102 - t144 * t103;
t88 = t105 * pkin(4) + t128;
t158 = t93 * t90 - t174;
t77 = t141 * t171 + t145 * t190;
t78 = -t141 * t190 + t145 * t171;
t156 = t102 * t196 - t103 * t197 + t140 * t78 + t144 * t77;
t154 = t189 - t129;
t153 = -t112 * t136 - t169;
t152 = -t142 * t209 + t145 * t198;
t149 = t161 * qJD(4) - t140 * t77 + t144 * t78;
t134 = t136 ^ 2;
t127 = -pkin(2) - t224;
t126 = t144 * pkin(3) + pkin(4);
t115 = t128 - t224;
t110 = 0.2e1 * t141 * t181;
t109 = t130 + t129;
t96 = t112 * t200;
t89 = -0.2e1 * t201 * t136 * qJD(3);
t80 = t88 - t224;
t67 = t191 + t227;
t58 = t64 * t135;
t57 = t63 * t135;
t49 = t130 + t54;
t38 = -t101 - t161;
t37 = t144 * t102 - t140 * t103 - t225;
t34 = -t90 ^ 2 + t92 ^ 2;
t31 = t90 * t135 + t47;
t30 = -t87 + t214;
t29 = t172 + t228;
t18 = t20 * t131;
t17 = t19 * t131;
t16 = t47 * t106 - t92 * t63;
t15 = t149 + t229;
t14 = t156 - t62;
t3 = -t47 * t105 - t106 * t48 + t63 * t90 - t92 * t64;
t2 = t10 * t61 - t164 * t19;
t1 = -t10 * t60 + t150 * t61 + t164 * t20 + t19 * t45;
t7 = [0, 0, 0, 0, -t136 * t130 - t123, t194 * t190, t110, t89, t132, -t202, 0, t127 * t182 - t125 * t132 + t96 + (t194 * t205 - t180) * t212, t125 * t202 + t127 * t181 - t152 * t212 + t210, t16, t3, -t57, -t58, 0, t109 * t90 + t115 * t48 + t149 * t135 + t215, t109 * t92 + t115 * t47 - t156 * t135 + t216, t2, t1, t17, -t18, 0, -t49 * t45 - t80 * t150 + (-t139 * t14 + t143 * t15 + (-t139 * t37 - t143 * t38) * qJD(5)) * t131 + t217, -t49 * t164 + t80 * t10 - (t139 * t15 + t143 * t14 + (-t139 * t38 + t143 * t37) * qJD(5)) * t131 + t218; 0, 0, 0, 0, t136 * t189 - t123, t193 * t188, t110, t89, t132, -t202, 0, -pkin(2) * t182 - pkin(7) * t132 + t96 + (t193 * t205 + t180) * t213, -pkin(2) * t181 + pkin(7) * t202 + t152 * t213 + t210, t16, t3, -t57, -t58, 0, t128 * t48 + t233 * t135 - t154 * t90 + t215, t128 * t47 + t234 * t135 - t154 * t92 + t216, t2, t1, t17, -t18, 0, -t88 * t150 - t168 * t45 + (t167 * t139 - t166 * t143) * t131 + t217, t88 * t10 - t168 * t164 + (t139 * t166 + t143 * t167) * t131 + t218; 0, 0, 0, 0, 0, 0, -t141 * t134 * t145, t201 * t134, 0, 0, 0, t153 * t141, t153 * t145, t221, t34, t31, 0, 0, -t90 * t191 - t220 - t172 * t135 + (t179 * t140 - t75) * qJD(4) + t175, -t92 * t191 + t214 * t135 + (t179 * qJD(4) - t55) * t144 + t158, t223, t6, t4, t5, 0, t67 * t45 - (-t139 * t30 + t143 * t29) * t131 + (-t139 * t144 - t207) * qJD(4) * t226 + ((-pkin(3) * t207 - t126 * t139) * t131 + t165) * qJD(5) + t159, t67 * t164 + (-t140 * t192 * t226 + t29 * t131 - t9) * t139 + (-qJD(5) * t25 - t8 + (-pkin(3) * t196 - qJD(5) * t126 + t30) * t131) * t143 + t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t221, t34, t31, 0, 0, -t163 * t135 + t151 - t220, t173 * t135 + t158 - t232, t223, t6, t4, t5, 0, t45 * t227 - (-t139 * t27 - t211) * t131 + (t178 * t139 - t211) * qJD(5) + t159, t164 * t227 + (qJD(5) * t178 + t27 * t131 - t8) * t143 + t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t223, t6, t4, t5, 0, t231 * t165 + t159, (-t231 * t25 - t8) * t143 + t236;];
tauc_reg = t7;
