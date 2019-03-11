% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPPPR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:19:44
% EndTime: 2019-03-09 08:19:53
% DurationCPUTime: 2.76s
% Computational Cost: add. (2177->352), mult. (5157->478), div. (0->0), fcn. (3113->6), ass. (0->201)
t150 = sin(qJ(2));
t219 = qJD(1) * t150;
t123 = -qJD(6) + t219;
t261 = qJD(6) + t123;
t252 = pkin(3) + pkin(7);
t131 = pkin(7) * t219;
t210 = pkin(3) * t219 + qJD(3) + t131;
t209 = qJD(1) * qJD(2);
t260 = -0.2e1 * t209;
t146 = sin(pkin(9));
t152 = cos(qJ(2));
t218 = qJD(1) * t152;
t198 = t146 * t218;
t147 = cos(pkin(9));
t216 = qJD(2) * t147;
t94 = -t198 + t216;
t91 = t94 ^ 2;
t201 = t147 * t218;
t217 = qJD(2) * t146;
t92 = t201 + t217;
t259 = -t92 ^ 2 - t91;
t149 = sin(qJ(6));
t151 = cos(qJ(6));
t48 = t149 * t92 + t151 * t94;
t258 = qJD(6) * t48;
t175 = t151 * t146 - t147 * t149;
t73 = t175 * t152;
t174 = t149 * t146 + t151 * t147;
t163 = t174 * t150;
t84 = t174 * qJD(6);
t248 = -qJD(1) * t163 + t84;
t257 = t248 * t123;
t203 = t146 * t219;
t191 = t151 * t203;
t202 = t147 * t219;
t249 = t175 * qJD(6) + t149 * t202 - t191;
t256 = t249 * t123;
t255 = t147 * qJD(5) - t210;
t143 = qJD(2) * qJ(3);
t254 = qJD(4) + t143;
t238 = qJ(5) * t146;
t251 = pkin(4) + pkin(5);
t169 = t251 * t147 + t238;
t144 = t150 ^ 2;
t154 = qJD(1) ^ 2;
t232 = t144 * t154;
t253 = qJD(2) * (t94 + t198) + t147 * t232;
t148 = -pkin(2) - qJ(4);
t250 = pkin(8) + t148;
t197 = t150 * t209;
t122 = pkin(2) * t197;
t176 = -qJ(3) * t152 + qJ(4) * t150;
t213 = qJD(3) * t150;
t157 = qJD(2) * t176 - qJD(4) * t152 - t213;
t45 = qJD(1) * t157 + t122;
t196 = t152 * t209;
t121 = pkin(7) * t196;
t134 = pkin(3) * t218;
t77 = t121 + (-qJD(4) + t134) * qJD(2);
t15 = t146 * t77 + t147 * t45;
t214 = qJD(2) * t152;
t103 = t252 * t214;
t215 = qJD(2) * t150;
t135 = pkin(2) * t215;
t57 = t135 + t157;
t29 = t146 * t103 + t147 * t57;
t195 = -t150 * qJ(3) - pkin(1);
t89 = t148 * t152 + t195;
t67 = t89 * qJD(1);
t71 = t148 * qJD(2) + t210;
t26 = t146 * t71 + t147 * t67;
t133 = pkin(7) * t218;
t102 = t133 + t134;
t136 = pkin(2) * t219;
t78 = qJD(1) * t176 + t136;
t40 = t146 * t102 + t147 * t78;
t115 = t252 * t150;
t51 = t146 * t115 + t147 * t89;
t247 = qJD(2) * pkin(2);
t246 = t123 * t48;
t25 = -t146 * t67 + t147 * t71;
t189 = qJD(5) - t25;
t20 = -pkin(4) * t219 + t189;
t245 = t152 * t20;
t21 = qJ(5) * t219 + t26;
t244 = t152 * t21;
t243 = t152 * t25;
t242 = t152 * t26;
t212 = qJD(6) * t151;
t241 = qJD(2) * t191 + t92 * t212;
t240 = -t169 * t219 + t255;
t182 = pkin(4) * t147 + t238;
t239 = t182 * t219 - t255;
t237 = qJD(4) * t92;
t236 = qJD(4) * t94;
t235 = qJD(5) * t94;
t234 = t123 * t149;
t233 = t123 * t151;
t231 = t146 * t150;
t230 = t146 * t152;
t228 = t147 * t150;
t227 = t147 * t152;
t225 = t152 * t154;
t153 = qJD(2) ^ 2;
t224 = t153 * t150;
t223 = t153 * t152;
t116 = t252 * t152;
t142 = qJD(2) * qJD(3);
t221 = -pkin(7) * t197 + t142;
t220 = -t152 ^ 2 + t144;
t111 = -pkin(2) * t152 + t195;
t88 = qJD(1) * t111;
t138 = t150 * qJD(5);
t32 = qJ(5) * t218 + t40;
t43 = t150 * qJ(5) + t51;
t208 = pkin(4) * t214;
t207 = t251 * t150;
t205 = t150 * t225;
t14 = -t146 * t45 + t147 * t77;
t162 = -pkin(8) * t231 - t251 * t152;
t159 = t162 * qJD(2);
t4 = qJD(1) * t159 - t14;
t190 = t147 * t197;
t7 = qJ(5) * t196 + qJD(1) * t138 + t15;
t5 = -pkin(8) * t190 + t7;
t204 = -t149 * t5 + t151 * t4;
t200 = t148 * t214;
t199 = qJD(5) * t230;
t39 = t102 * t147 - t146 * t78;
t28 = t103 * t147 - t146 * t57;
t50 = t115 * t147 - t146 * t89;
t194 = pkin(1) * t260;
t193 = qJD(3) - t247;
t192 = t147 * qJ(5) - qJ(3);
t18 = qJ(5) * t214 + t138 + t29;
t82 = t102 + t254;
t9 = -pkin(4) * t196 - t14;
t187 = t146 * t7 - t147 * t9;
t186 = t149 * t4 + t151 * t5;
t6 = -pkin(8) * t94 - qJD(1) * t207 + t189;
t8 = pkin(8) * t92 + t21;
t185 = t149 * t8 - t151 * t6;
t2 = t149 * t6 + t151 * t8;
t184 = -t82 + t254;
t108 = pkin(4) * t146 - t192;
t164 = qJ(5) * t94 - t82;
t30 = pkin(4) * t92 - t164;
t183 = -qJD(2) * t108 - qJD(4) + t30;
t181 = t14 * t147 + t146 * t15;
t180 = t146 * t20 + t147 * t21;
t179 = -t146 * t25 + t147 * t26;
t27 = pkin(8) * t230 - t207 - t50;
t33 = pkin(8) * t227 + t43;
t178 = -t149 * t33 + t151 * t27;
t177 = t149 * t27 + t151 * t33;
t173 = -0.2e1 * qJD(2) * t88;
t172 = -t92 * t202 + t94 * t203;
t171 = -pkin(3) - t182;
t165 = -qJ(3) * t214 - t213;
t68 = qJD(1) * t165 + t122;
t80 = t135 + t165;
t170 = pkin(7) * t153 + qJD(1) * t80 + t68;
t168 = (t94 - t216) * t219;
t105 = t250 * t146;
t167 = -qJD(1) * t162 + qJD(4) * t147 - qJD(6) * t105 + t39;
t106 = t250 * t147;
t166 = pkin(8) * t202 - qJD(4) * t146 - qJD(6) * t106 - t32;
t76 = -pkin(3) * t197 + t221;
t161 = pkin(3) + t169;
t160 = qJD(2) * t163;
t16 = (-qJD(6) * t94 - t190) * t149 + t241;
t158 = t171 * t197 + t221;
t110 = t131 + t193;
t114 = -t133 - t143;
t156 = t221 * t152 + (t110 * t152 + (t114 + t133) * t150) * qJD(2);
t155 = -t174 * t197 - t258;
t104 = t147 * t148 * t196;
t101 = t252 * t215;
t99 = -qJ(3) * t218 + t136;
t85 = -t251 * t146 + t192;
t72 = t174 * t152;
t70 = t88 * t219;
t60 = t152 * t182 + t116;
t56 = (-t92 + t217) * t219;
t54 = -t169 * t152 - t116;
t49 = -t146 * t232 + (-t92 + t201) * qJD(2);
t46 = t149 * t94 - t151 * t92;
t44 = -pkin(4) * t150 - t50;
t42 = t199 + (-pkin(7) + t171) * t215;
t36 = t152 * t84 + t175 * t215;
t35 = -qJD(6) * t73 + t160;
t34 = -pkin(4) * t218 - t39;
t31 = -t199 + (pkin(7) + t161) * t215;
t23 = -t28 - t208;
t22 = t158 - t235;
t17 = qJD(1) * t160 + t258;
t13 = -t251 * t92 + t164;
t12 = t161 * t197 - t221 + t235;
t11 = -pkin(8) * t147 * t215 + t18;
t10 = t159 - t28;
t1 = [0, 0, 0, 0.2e1 * t150 * t196, t220 * t260, t223, -t224, 0, -pkin(7) * t223 + t150 * t194, pkin(7) * t224 + t152 * t194, t156, t150 * t173 + t152 * t170, -t150 * t170 + t152 * t173, pkin(7) * t156 + t111 * t68 + t80 * t88, t76 * t227 - t101 * t92 + (qJD(1) * t28 + t14) * t150 + (-t82 * t228 + t243 + (-t116 * t228 + t152 * t50) * qJD(1)) * qJD(2), -t76 * t230 - t101 * t94 + (-qJD(1) * t29 - t15) * t150 + (t82 * t231 - t242 + (t116 * t231 - t152 * t51) * qJD(1)) * qJD(2), -t28 * t94 - t29 * t92 + (t14 * t146 - t147 * t15) * t152 + ((-t146 * t50 + t147 * t51) * qJD(1) + t179) * t215, -t101 * t82 + t116 * t76 + t14 * t50 + t15 * t51 + t25 * t28 + t26 * t29, t22 * t227 + t42 * t92 + (-qJD(1) * t23 - t9) * t150 + (-t30 * t228 - t245 + (-t152 * t44 - t228 * t60) * qJD(1)) * qJD(2), -t18 * t92 + t23 * t94 + (-t146 * t9 - t147 * t7) * t152 + ((t146 * t44 + t147 * t43) * qJD(1) + t180) * t215, t22 * t230 - t42 * t94 + (qJD(1) * t18 + t7) * t150 + (-t30 * t231 + t244 + (t152 * t43 - t231 * t60) * qJD(1)) * qJD(2), t18 * t21 + t20 * t23 + t22 * t60 + t30 * t42 + t43 * t7 + t44 * t9, -t16 * t73 + t36 * t48, t16 * t72 + t17 * t73 - t35 * t48 - t36 * t46, -t36 * t123 - t16 * t150 + (qJD(1) * t73 - t48) * t214, t35 * t123 + t17 * t150 + (-qJD(1) * t72 + t46) * t214 (t123 + t219) * t214 -(t151 * t10 - t149 * t11) * t123 - t204 * t150 + t31 * t46 + t54 * t17 - t12 * t72 + t13 * t35 + (t123 * t177 + t150 * t2) * qJD(6) + (-qJD(1) * t178 + t185) * t214 (t149 * t10 + t151 * t11) * t123 + t186 * t150 + t31 * t48 + t54 * t16 - t12 * t73 + t13 * t36 + (t123 * t178 - t150 * t185) * qJD(6) + (qJD(1) * t177 + t2) * t214; 0, 0, 0, -t205, t220 * t154, 0, 0, 0, t154 * pkin(1) * t150, pkin(1) * t225 ((-t114 - t143) * t150 + (-t110 + t193) * t152) * qJD(1), -t218 * t99 + t70, 0.2e1 * t142 + (t150 * t99 + t152 * t88) * qJD(1), qJ(3) * t221 - qJD(3) * t114 - t88 * t99 + (-t114 * t150 + (-t110 - t247) * t152) * qJD(1) * pkin(7), t146 * t76 + t104 + t210 * t92 + (-t243 + (-t147 * t184 - t39) * t150) * qJD(1), t147 * t76 + t210 * t94 + (t150 * t40 + t242 + (t150 * t184 - t200) * t146) * qJD(1), t39 * t94 + t40 * t92 + (-t219 * t26 - t14 + t236) * t147 + (t219 * t25 - t15 + t237) * t146, qJ(3) * t76 - t25 * t39 - t26 * t40 + t210 * t82 + t181 * t148 + (-t146 * t26 - t147 * t25) * qJD(4), t146 * t22 + t104 + t239 * t92 + (t245 + (t147 * t183 + t34) * t150) * qJD(1), t32 * t92 - t34 * t94 + (-t21 * t219 + t236 + t9) * t147 + (-t20 * t219 + t237 - t7) * t146, -t147 * t22 - t239 * t94 + (-t150 * t32 - t244 + (t150 * t183 + t200) * t146) * qJD(1), t108 * t22 - t20 * t34 - t21 * t32 + t239 * t30 + t187 * t148 + (-t146 * t21 + t147 * t20) * qJD(4), t16 * t174 + t249 * t48, t16 * t175 - t17 * t174 - t248 * t48 - t249 * t46, -t256 + (-qJD(2) * t174 + t48) * t218, t257 + (-qJD(2) * t175 - t46) * t218, -t123 * t218, -t12 * t175 + t85 * t17 + t240 * t46 + t248 * t13 + (t149 * t166 - t151 * t167) * t123 + (-(-t105 * t149 - t106 * t151) * qJD(2) - t185) * t218, t12 * t174 + t85 * t16 + t240 * t48 + t249 * t13 + (t149 * t167 + t151 * t166) * t123 + ((t105 * t151 - t106 * t149) * qJD(2) - t2) * t218; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t205, -t153 - t232, qJD(2) * t114 + t121 + t70, t49, -t253, t172, -qJD(2) * t82 + t179 * t219 + t181, t49, t172, t253, -qJD(2) * t30 + t180 * t219 + t187, 0, 0, 0, 0, 0, t256 + (t174 * t218 + t46) * qJD(2), -t257 + (t175 * t218 + t48) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t168, t56, t259, t25 * t94 + t26 * t92 + t76, t168, t259, -t56, t21 * t92 + (-qJD(5) - t20) * t94 + t158, 0, 0, 0, 0, 0, t155 + t246, t92 * t233 + (t190 + (qJD(6) - t123) * t94) * t149 - t241; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92 * t94 - t196 (t92 + t217) * t219, -t91 - t232, t30 * t94 + (-t150 * t21 - t208) * qJD(1) - t14, 0, 0, 0, 0, 0, qJD(6) * t234 - t94 * t46 + (-t150 * t234 - t151 * t214) * qJD(1), t123 * t212 - t48 * t94 + (t149 * t214 - t150 * t233) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48 * t46, -t46 ^ 2 + t48 ^ 2, -t123 * t46 + t16, t155 - t246, -t196, -t13 * t48 - t261 * t2 + t204, t13 * t46 + t261 * t185 - t186;];
tauc_reg  = t1;
