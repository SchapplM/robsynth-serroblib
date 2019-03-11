% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRRR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:29:55
% EndTime: 2019-03-08 20:30:03
% DurationCPUTime: 2.50s
% Computational Cost: add. (2411->316), mult. (6124->464), div. (0->0), fcn. (4788->12), ass. (0->181)
t138 = sin(qJ(4));
t142 = cos(qJ(4));
t169 = pkin(4) * t138 - pkin(9) * t142;
t109 = t169 * qJD(4);
t134 = cos(pkin(12));
t139 = sin(qJ(2));
t133 = sin(pkin(6));
t204 = qJD(1) * t133;
t187 = t139 * t204;
t113 = t134 * t187;
t132 = sin(pkin(12));
t143 = cos(qJ(2));
t186 = t143 * t204;
t75 = t132 * t186 + t113;
t251 = t109 - t75;
t137 = sin(qJ(5));
t141 = cos(qJ(5));
t194 = t141 * qJD(4);
t203 = qJD(2) * t138;
t100 = t137 * t203 - t194;
t135 = cos(pkin(6));
t121 = qJD(1) * t135 + qJD(3);
t111 = qJD(2) * pkin(2) + t186;
t69 = t132 * t111 + t113;
t62 = qJD(2) * pkin(8) + t69;
t242 = t121 * t142 - t138 * t62;
t36 = -qJD(4) * pkin(4) - t242;
t25 = pkin(5) * t100 + t36;
t201 = qJD(4) * t137;
t102 = t141 * t203 + t201;
t136 = sin(qJ(6));
t140 = cos(qJ(6));
t51 = t140 * t100 + t102 * t136;
t250 = t25 * t51;
t161 = t100 * t136 - t140 * t102;
t249 = t161 * t51;
t198 = qJD(5) * t137;
t185 = t138 * t198;
t241 = t142 * t194 - t185;
t248 = t161 ^ 2 - t51 ^ 2;
t202 = qJD(2) * t142;
t122 = -qJD(5) + t202;
t120 = -qJD(6) + t122;
t195 = qJD(6) * t140;
t196 = qJD(6) * t136;
t192 = qJD(4) * qJD(5);
t72 = t241 * qJD(2) + t141 * t192;
t199 = qJD(4) * t142;
t183 = t137 * t199;
t197 = qJD(5) * t141;
t184 = t138 * t197;
t149 = t183 + t184;
t73 = qJD(2) * t149 + t137 * t192;
t15 = -t100 * t195 - t102 * t196 - t136 * t73 + t140 * t72;
t247 = -t120 * t51 + t15;
t193 = qJD(2) * qJD(4);
t125 = t138 * t193;
t158 = -pkin(4) * t142 - pkin(9) * t138 - pkin(3);
t112 = t132 * t187;
t68 = t111 * t134 - t112;
t45 = t158 * qJD(2) - t68;
t225 = t137 * t45;
t40 = t138 * t121 + t142 * t62;
t37 = qJD(4) * pkin(9) + t40;
t14 = t141 * t37 + t225;
t81 = (t132 * t139 - t134 * t143) * t133;
t77 = qJD(2) * t81;
t71 = qJD(1) * t77;
t21 = qJD(4) * t242 - t142 * t71;
t82 = (t132 * t143 + t134 * t139) * t133;
t46 = (qJD(1) * t82 + t109) * qJD(2);
t177 = t137 * t21 - t141 * t46;
t147 = -t14 * qJD(5) - t177;
t2 = pkin(5) * t125 - pkin(10) * t72 + t147;
t191 = t137 * t46 + t141 * t21 + t45 * t197;
t153 = -t37 * t198 + t191;
t3 = -pkin(10) * t73 + t153;
t190 = -t136 * t3 + t140 * t2;
t11 = -pkin(10) * t100 + t14;
t227 = t11 * t140;
t13 = -t137 * t37 + t141 * t45;
t10 = -pkin(10) * t102 + t13;
t8 = -pkin(5) * t122 + t10;
t5 = t136 * t8 + t227;
t246 = -t5 * qJD(6) + t25 * t161 + t190;
t146 = t161 * qJD(6) - t136 * t72 - t140 * t73;
t245 = t120 * t161 + t146;
t200 = qJD(4) * t138;
t210 = t137 * t142;
t123 = pkin(2) * t132 + pkin(8);
t213 = t123 * t137;
t78 = t134 * t186 - t112;
t244 = -t251 * t141 - t200 * t213 - t78 * t210;
t208 = t141 * t142;
t236 = pkin(2) * t134;
t98 = t158 - t236;
t243 = t251 * t137 + t98 * t197 - t78 * t208;
t104 = t136 * t141 + t137 * t140;
t84 = t104 * t138;
t240 = qJD(5) + qJD(6);
t180 = qJD(6) * t8 + t3;
t9 = t11 * t196;
t239 = t136 * t2 + t140 * t180 - t9;
t130 = t138 ^ 2;
t159 = qJD(2) * t130 - t122 * t142;
t238 = -t122 * t185 - t159 * t194;
t237 = pkin(9) + pkin(10);
t105 = t123 * t208;
t157 = pkin(5) * t138 - pkin(10) * t208;
t235 = -t157 * qJD(4) - (-t105 + (pkin(10) * t138 - t98) * t137) * qJD(5) + t244;
t234 = (-t138 * t194 - t142 * t198) * t123 - t149 * pkin(10) + t243;
t209 = t138 * t141;
t211 = t137 * t138;
t27 = -t196 * t211 + (t240 * t209 + t183) * t140 + t241 * t136;
t233 = t27 * t120 - t84 * t125;
t106 = t169 * qJD(2);
t232 = t137 * t106 + t141 * t242;
t103 = t136 * t137 - t140 * t141;
t152 = t103 * t142;
t231 = qJD(2) * t152 - t240 * t103;
t230 = (-t202 + t240) * t104;
t226 = t137 * t36;
t224 = t137 * t72;
t223 = t141 * t36;
t222 = t142 * t73;
t22 = t121 * t200 - t138 * t71 + t62 * t199;
t221 = t22 * t137;
t220 = t22 * t141;
t219 = t137 * t98 + t105;
t218 = t100 * t122;
t217 = t100 * t138;
t216 = t102 * t122;
t214 = t122 * t141;
t145 = qJD(2) ^ 2;
t212 = t133 * t145;
t144 = qJD(4) ^ 2;
t207 = t144 * t138;
t206 = t144 * t142;
t205 = -t142 ^ 2 + t130;
t188 = qJD(5) * t237;
t182 = t137 * t202;
t176 = t141 * t106 - t137 * t242;
t175 = -t142 * t15 - t161 * t200;
t174 = t102 * t200 - t142 * t72;
t173 = t122 * t123 + t37;
t171 = t122 * t184;
t170 = -t40 + (-t182 + t198) * pkin(5);
t117 = t237 * t137;
t168 = pkin(10) * t182 - qJD(6) * t117 - t137 * t188 - t232;
t118 = t237 * t141;
t167 = t157 * qJD(2) + qJD(6) * t118 + t141 * t188 + t176;
t60 = t135 * t138 + t142 * t82;
t30 = -t137 * t60 + t141 * t81;
t31 = t137 * t81 + t141 * t60;
t166 = -t136 * t31 + t140 * t30;
t165 = t136 * t30 + t140 * t31;
t87 = t141 * t98;
t44 = -pkin(10) * t209 + t87 + (-pkin(5) - t213) * t142;
t48 = -pkin(10) * t211 + t219;
t164 = t136 * t44 + t140 * t48;
t162 = t135 * t142 - t138 * t82;
t156 = -t142 * t146 - t51 * t200;
t76 = qJD(2) * t82;
t70 = qJD(1) * t76;
t155 = qJD(2) * t75 - t123 * t144 - t70;
t61 = -qJD(2) * pkin(3) - t68;
t154 = qJD(4) * (qJD(2) * (-pkin(3) - t236) + t61 + t78);
t151 = t159 * t137;
t26 = -qJD(4) * t152 - t240 * t84;
t85 = t103 * t138;
t150 = t120 * t26 + t85 * t125;
t128 = -pkin(5) * t141 - pkin(4);
t90 = (pkin(5) * t137 + t123) * t138;
t64 = pkin(5) * t149 + t123 * t199;
t29 = t162 * qJD(4) - t142 * t77;
t28 = t60 * qJD(4) - t138 * t77;
t12 = pkin(5) * t73 + t22;
t7 = t30 * qJD(5) + t137 * t76 + t141 * t29;
t6 = -t31 * qJD(5) - t137 * t29 + t141 * t76;
t4 = -t11 * t136 + t140 * t8;
t1 = [0, 0, -t139 * t212, -t143 * t212, -t68 * t76 - t69 * t77 + t70 * t81 - t71 * t82, 0, 0, 0, 0, 0, -qJD(4) * t28 + (-t142 * t76 + t81 * t200) * qJD(2), -qJD(4) * t29 + (t138 * t76 + t81 * t199) * qJD(2), 0, 0, 0, 0, 0, t100 * t28 - t122 * t6 + t125 * t30 - t162 * t73, t102 * t28 + t122 * t7 - t125 * t31 - t162 * t72, 0, 0, 0, 0, 0 -(-qJD(6) * t165 - t136 * t7 + t140 * t6) * t120 + t166 * t125 + t28 * t51 + t162 * t146 (qJD(6) * t166 + t136 * t6 + t140 * t7) * t120 - t165 * t125 - t28 * t161 - t162 * t15; 0, 0, 0, 0, t68 * t75 - t69 * t78 + (-t132 * t71 - t134 * t70) * pkin(2), 0.2e1 * t142 * t125, -0.2e1 * t205 * t193, t206, -t207, 0, t138 * t154 + t142 * t155, -t138 * t155 + t142 * t154, t241 * t102 + t72 * t209 (-t100 * t141 - t102 * t137) * t199 + (-t224 - t141 * t73 + (t100 * t137 - t102 * t141) * qJD(5)) * t138, t174 - t238, t171 + t222 + (-t151 - t217) * qJD(4) (-t122 - t202) * t200 (t198 * t98 + t244) * t122 + ((t100 * t123 + t226) * qJD(4) + (t141 * t173 + t225) * qJD(5) + t177) * t142 + (t36 * t197 - t78 * t100 + t123 * t73 + t221 + ((-t123 * t210 + t87) * qJD(2) + t13) * qJD(4)) * t138, t243 * t122 + (-t173 * t198 + (t102 * t123 + t223) * qJD(4) + t191) * t142 + (-t36 * t198 - t78 * t102 + t123 * t72 + t220 + (-t219 * qJD(2) - t123 * t214 - t14) * qJD(4)) * t138, -t15 * t85 - t161 * t26, -t146 * t85 - t15 * t84 + t161 * t27 - t26 * t51, -t150 + t175, t156 + t233 (-t120 - t202) * t200, -t190 * t142 + t64 * t51 - t90 * t146 + t12 * t84 + t25 * t27 + (t234 * t136 + t235 * t140) * t120 + (t120 * t164 + t142 * t5) * qJD(6) + (-t78 * t51 + ((-t136 * t48 + t140 * t44) * qJD(2) + t4) * qJD(4)) * t138, t239 * t142 - t64 * t161 + t90 * t15 - t12 * t85 + t25 * t26 + ((qJD(6) * t44 + t234) * t140 + (-qJD(6) * t48 - t235) * t136) * t120 + (t78 * t161 + (-qJD(2) * t164 - t5) * qJD(4)) * t138; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t207, -t206, 0, 0, 0, 0, 0, t171 - t222 + (-t151 + t217) * qJD(4), t174 + t238, 0, 0, 0, 0, 0, -t156 + t233, t150 + t175; 0, 0, 0, 0, 0, -t138 * t145 * t142, t205 * t145, 0, 0, 0, qJD(4) * t40 - t61 * t203 - t22 (-qJD(2) * t61 + t71) * t142, -t102 * t214 + t224 (t72 + t218) * t141 + (-t73 + t216) * t137, -t122 * t197 + (t122 * t208 + (-t102 + t201) * t138) * qJD(2), t122 * t198 + (-t122 * t210 + (t100 + t194) * t138) * qJD(2), t122 * t203, -pkin(4) * t73 - t220 + t176 * t122 - t40 * t100 + (pkin(9) * t214 + t226) * qJD(5) + (-t13 * t138 + (-pkin(9) * t200 - t142 * t36) * t137) * qJD(2), -pkin(4) * t72 + t221 - t232 * t122 - t40 * t102 + (-pkin(9) * t122 * t137 + t223) * qJD(5) + (-t36 * t208 + (-pkin(9) * t194 + t14) * t138) * qJD(2), t104 * t15 - t161 * t231, -t103 * t15 + t104 * t146 + t161 * t230 - t231 * t51, -t231 * t120 + (qJD(4) * t104 + t161) * t203, t230 * t120 + (-qJD(4) * t103 + t51) * t203, t120 * t203, t12 * t103 - t128 * t146 + t170 * t51 + t230 * t25 + (t136 * t168 + t140 * t167) * t120 + ((-t117 * t140 - t118 * t136) * qJD(4) - t4) * t203, t12 * t104 + t128 * t15 - t170 * t161 + t231 * t25 + (-t136 * t167 + t140 * t168) * t120 + (-(-t117 * t136 + t118 * t140) * qJD(4) + t5) * t203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102 * t100, -t100 ^ 2 + t102 ^ 2, t72 - t218, -t216 - t73, t125, -t102 * t36 - t122 * t14 + t147, t100 * t36 - t122 * t13 - t153, -t249, t248, t247, t245, t125 (-t10 * t136 - t227) * t120 + (-t102 * t51 + t120 * t196 + t125 * t140) * pkin(5) + t246, t250 + t9 + (t11 * t120 - t2) * t136 + (-t10 * t120 - t180) * t140 + (t102 * t161 + t120 * t195 - t125 * t136) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t249, t248, t247, t245, t125, -t120 * t5 + t246, -t120 * t4 - t239 + t250;];
tauc_reg  = t1;
