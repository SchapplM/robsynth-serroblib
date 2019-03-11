% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
% 
% Output:
% tauc_reg [6x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRRR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:23:48
% EndTime: 2019-03-09 02:23:52
% DurationCPUTime: 1.86s
% Computational Cost: add. (1963->281), mult. (4258->407), div. (0->0), fcn. (2817->8), ass. (0->156)
t114 = sin(qJ(6));
t117 = cos(qJ(6));
t115 = sin(qJ(5));
t118 = cos(qJ(5));
t165 = t118 * qJD(4);
t119 = cos(qJ(4));
t174 = qJD(1) * t119;
t81 = t115 * t174 - t165;
t155 = t118 * t174;
t173 = qJD(4) * t115;
t83 = t155 + t173;
t137 = t114 * t81 - t117 * t83;
t28 = t114 * t83 + t117 * t81;
t220 = t137 * t28;
t211 = qJD(5) + qJD(6);
t84 = t114 * t115 - t117 * t118;
t219 = t211 * t84;
t85 = t114 * t118 + t115 * t117;
t122 = t211 * t85;
t116 = sin(qJ(4));
t172 = qJD(4) * t116;
t154 = t115 * t172;
t168 = qJD(5) * t119;
t126 = -t118 * t168 + t154;
t218 = t137 ^ 2 - t28 ^ 2;
t166 = qJD(6) * t117;
t167 = qJD(6) * t114;
t152 = t115 * t168;
t212 = -t116 * t165 - t152;
t48 = t212 * qJD(1) + qJD(5) * t165;
t49 = -qJD(1) * t154 + t83 * qJD(5);
t7 = -t114 * t49 + t117 * t48 - t81 * t166 - t83 * t167;
t175 = qJD(1) * t116;
t100 = qJD(5) + t175;
t98 = qJD(6) + t100;
t217 = t28 * t98 + t7;
t101 = sin(pkin(10)) * pkin(1) + qJ(3);
t75 = pkin(4) * t116 - pkin(8) * t119 + t101;
t55 = t75 * qJD(1);
t197 = t115 * t55;
t109 = t119 * qJD(2);
t99 = -cos(pkin(10)) * pkin(1) - pkin(2) - pkin(7);
t78 = t99 * qJD(1) + qJD(3);
t52 = t116 * t78 + t109;
t44 = qJD(4) * pkin(8) + t52;
t20 = t118 * t44 + t197;
t15 = -pkin(9) * t81 + t20;
t12 = t15 * t167;
t213 = -t116 * qJD(2) + t119 * t78;
t43 = -qJD(4) * pkin(4) - t213;
t24 = pkin(5) * t81 + t43;
t216 = t24 * t28 + t12;
t123 = t137 * qJD(6) - t114 * t48 - t117 * t49;
t215 = -t137 * t98 + t123;
t164 = qJD(1) * qJD(4);
t102 = t119 * t164;
t45 = t213 * qJD(4);
t139 = pkin(4) * t119 + pkin(8) * t116;
t79 = t139 * qJD(4) + qJD(3);
t60 = t79 * qJD(1);
t146 = -t115 * t45 + t118 * t60;
t124 = -t20 * qJD(5) + t146;
t4 = pkin(5) * t102 - pkin(9) * t48 + t124;
t169 = qJD(5) * t118;
t163 = -t115 * t60 - t118 * t45 - t55 * t169;
t170 = qJD(5) * t115;
t131 = -t44 * t170 - t163;
t5 = -pkin(9) * t49 + t131;
t159 = -t114 * t5 + t117 * t4;
t19 = -t115 * t44 + t118 * t55;
t14 = -pkin(9) * t83 + t19;
t11 = pkin(5) * t100 + t14;
t195 = t117 * t15;
t2 = t11 * t114 + t195;
t214 = -t2 * qJD(6) + t24 * t137 + t159;
t111 = t119 ^ 2;
t176 = qJD(1) * t111;
t210 = -(-t100 * t116 + t176) * t165 + t100 * t152;
t209 = 0.2e1 * qJD(3);
t208 = pkin(8) + pkin(9);
t171 = qJD(4) * t119;
t207 = t7 * t116 - t137 * t171;
t181 = t118 * t119;
t183 = t115 * t119;
t17 = -t167 * t183 + (t211 * t181 - t154) * t117 + t212 * t114;
t63 = t85 * t119;
t206 = -t63 * t102 - t17 * t98;
t130 = t84 * t116;
t205 = -qJD(1) * t130 - t219;
t129 = qJD(1) * t85;
t204 = t116 * t129 + t122;
t203 = t48 * t116 + t83 * t171;
t86 = t139 * qJD(1);
t202 = t115 * t86 + t118 * t213;
t182 = t116 * t118;
t80 = t99 * t182;
t201 = t115 * t75 + t80;
t200 = t100 * t81;
t199 = t100 * t83;
t198 = t115 * t43;
t196 = t116 * t49;
t194 = t118 * t43;
t193 = t119 * t48;
t191 = t119 * t81;
t46 = qJD(4) * t109 + t78 * t172;
t190 = t46 * t115;
t189 = t46 * t118;
t188 = t48 * t115;
t187 = qJD(4) * t98;
t186 = t100 * t115;
t185 = t100 * t118;
t184 = t115 * t116;
t120 = qJD(4) ^ 2;
t180 = t120 * t116;
t179 = t120 * t119;
t178 = t116 ^ 2 - t111;
t121 = qJD(1) ^ 2;
t177 = -t120 - t121;
t89 = qJD(1) * t101;
t162 = t119 * t99 * t165 + t115 * t79 + t75 * t169;
t161 = pkin(9) * t182;
t160 = t99 * t184;
t158 = qJD(5) * t208;
t157 = t115 * t176;
t156 = t115 * t175;
t150 = -t115 * t99 + pkin(5);
t149 = qJD(6) * t11 + t5;
t148 = t100 * t99 + t44;
t145 = -t115 * t213 + t118 * t86;
t144 = qJD(5) * t116 + qJD(1);
t142 = -t52 + (t156 + t170) * pkin(5);
t95 = t208 * t118;
t141 = qJD(6) * t95 + (pkin(5) * t119 + t161) * qJD(1) + t145 + t118 * t158;
t94 = t208 * t115;
t140 = pkin(9) * t156 + qJD(6) * t94 + t115 * t158 + t202;
t62 = t118 * t75;
t23 = -pkin(9) * t181 + t150 * t116 + t62;
t25 = -pkin(9) * t183 + t201;
t138 = t114 * t23 + t117 * t25;
t134 = t126 * t100;
t133 = -pkin(8) * t171 + t116 * t43;
t132 = t116 * t123 - t28 * t171;
t128 = t84 * qJD(1);
t16 = qJD(4) * t130 - t122 * t119;
t64 = t84 * t119;
t127 = t64 * t102 - t16 * t98;
t106 = -pkin(5) * t118 - pkin(4);
t97 = t116 * t102;
t70 = t118 * t79;
t65 = (pkin(5) * t115 - t99) * t119;
t35 = -t126 * pkin(5) + t99 * t172;
t21 = pkin(5) * t49 + t46;
t10 = t126 * pkin(9) - qJD(5) * t160 + t162;
t9 = t70 + (-t80 + (pkin(9) * t119 - t75) * t115) * qJD(5) + (t150 * t119 + t161) * qJD(4);
t1 = t11 * t117 - t114 * t15;
t3 = [0, 0, 0, 0, 0, qJD(1) * t209, t89 * t209, -0.2e1 * t97, 0.2e1 * t178 * t164, -t180, -t179, 0, t89 * t171 - t99 * t180 + (t101 * t171 + t116 * t209) * qJD(1), -t89 * t172 - t99 * t179 + (-t101 * t172 + t119 * t209) * qJD(1), -t83 * t152 + (-t83 * t172 + t193) * t118 (t115 * t83 + t118 * t81) * t172 + (-t188 - t118 * t49 + (t115 * t81 - t118 * t83) * qJD(5)) * t119, t203 - t210, -t196 + (-t157 - t191) * qJD(4) + t134, t100 * t171 + t97 (-t75 * t170 + t70) * t100 + ((t81 * t99 - t198) * qJD(4) + (-t148 * t118 - t197) * qJD(5) + t146) * t116 + (t43 * t169 + t190 - t99 * t49 + (-t99 * t186 + (t62 - t160) * qJD(1) + t19) * qJD(4)) * t119, -t162 * t100 + (t148 * t170 + (t83 * t99 - t194) * qJD(4) + t163) * t116 + (-t43 * t170 + t189 - t99 * t48 + (-t201 * qJD(1) - t20) * qJD(4)) * t119, -t137 * t16 - t64 * t7, -t123 * t64 + t137 * t17 - t16 * t28 - t63 * t7, -t127 + t207, t132 + t206, t171 * t98 + t97 (-t10 * t114 + t117 * t9) * t98 + t159 * t116 + t35 * t28 - t65 * t123 + t21 * t63 + t24 * t17 + (-t116 * t2 - t138 * t98) * qJD(6) + ((-t114 * t25 + t117 * t23) * qJD(1) + t1) * t171, t12 * t116 + t24 * t16 - t21 * t64 - t35 * t137 + t65 * t7 + (-(-qJD(6) * t25 + t9) * t98 - t4 * t116) * t114 + (-(qJD(6) * t23 + t10) * t98 - t149 * t116) * t117 + (-qJD(1) * t138 - t2) * t171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t179, t180, 0, 0, 0, 0, 0, t196 + (-t157 + t191) * qJD(4) + t134, t203 + t210, 0, 0, 0, 0, 0, -t132 + t206, t127 + t207; 0, 0, 0, 0, 0, -t121, -t89 * qJD(1), 0, 0, 0, 0, 0, t177 * t116, t177 * t119, 0, 0, 0, 0, 0, -t119 * t49 - t144 * t185 + (t116 * t81 + (-t100 - t175) * t183) * qJD(4), -t193 + t144 * t186 + (-t100 * t181 + (t83 - t155) * t116) * qJD(4), 0, 0, 0, 0, 0, t98 * t128 + (-t85 * t187 + t123) * t119 + ((-t174 * t85 + t28) * qJD(4) + t98 * t219) * t116, t98 * t129 + (t84 * t187 - t7) * t119 + (t122 * t98 + (t119 * t128 - t137) * qJD(4)) * t116; 0, 0, 0, 0, 0, 0, 0, t119 * t121 * t116, -t178 * t121, 0, 0, 0, qJD(4) * t52 - t89 * t174 - t46, t89 * t175, t83 * t185 + t188 (t48 - t200) * t118 + (-t49 - t199) * t115, t100 * t169 + (t100 * t182 + (-t83 + t173) * t119) * qJD(1), -t100 * t170 + (-t100 * t184 + (t81 + t165) * t119) * qJD(1), -t100 * t174, -pkin(4) * t49 - t189 - t145 * t100 - t52 * t81 + (-pkin(8) * t185 + t198) * qJD(5) + (t133 * t115 - t19 * t119) * qJD(1), -pkin(4) * t48 + t190 + t202 * t100 - t52 * t83 + (pkin(8) * t186 + t194) * qJD(5) + (t133 * t118 + t20 * t119) * qJD(1), -t137 * t205 + t7 * t85, t123 * t85 + t137 * t204 - t205 * t28 - t7 * t84, t205 * t98 + (qJD(4) * t85 + t137) * t174, -t204 * t98 + (-qJD(4) * t84 + t28) * t174, -t98 * t174, -t106 * t123 + t21 * t84 + (t114 * t140 - t117 * t141) * t98 + t142 * t28 + t204 * t24 + ((-t114 * t95 - t117 * t94) * qJD(4) - t1) * t174, t106 * t7 + t21 * t85 + (t114 * t141 + t117 * t140) * t98 - t142 * t137 + t205 * t24 + (-(-t114 * t94 + t117 * t95) * qJD(4) + t2) * t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83 * t81, -t81 ^ 2 + t83 ^ 2, t48 + t200, t199 - t49, t102, t100 * t20 - t43 * t83 + t124, t100 * t19 + t43 * t81 - t131, -t220, t218, t217, t215, t102 -(-t114 * t14 - t195) * t98 + (t102 * t117 - t167 * t98 - t28 * t83) * pkin(5) + t214 (-t15 * t98 - t4) * t114 + (t14 * t98 - t149) * t117 + (-t102 * t114 + t137 * t83 - t166 * t98) * pkin(5) + t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t220, t218, t217, t215, t102, t2 * t98 + t214, t1 * t98 - t114 * t4 - t117 * t149 + t216;];
tauc_reg  = t3;
