% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tauc_reg [5x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRP8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:02:08
% EndTime: 2019-12-31 22:02:16
% DurationCPUTime: 2.63s
% Computational Cost: add. (3054->309), mult. (7625->449), div. (0->0), fcn. (5144->6), ass. (0->168)
t149 = sin(qJ(3));
t237 = pkin(7) + pkin(8);
t189 = qJD(3) * t237;
t150 = sin(qJ(2));
t151 = cos(qJ(3));
t211 = t150 * t151;
t152 = cos(qJ(2));
t212 = t149 * t152;
t168 = pkin(2) * t150 - pkin(7) * t152;
t106 = t168 * qJD(1);
t86 = t149 * t106;
t257 = -t149 * t189 - t86 - (-pkin(6) * t211 - pkin(8) * t212) * qJD(1);
t210 = t151 * t152;
t165 = pkin(3) * t150 - pkin(8) * t210;
t204 = qJD(1) * t150;
t184 = t149 * t204;
t221 = pkin(6) * t184 + t151 * t106;
t256 = t165 * qJD(1) + t151 * t189 + t221;
t202 = qJD(2) * t149;
t100 = t151 * t204 + t202;
t148 = sin(qJ(4));
t235 = cos(qJ(4));
t200 = qJD(2) * t151;
t99 = -t184 + t200;
t164 = t235 * t100 + t148 * t99;
t238 = t164 ^ 2;
t55 = t100 * t148 - t235 * t99;
t53 = t55 ^ 2;
t255 = -t53 + t238;
t214 = t148 * t149;
t162 = t235 * t151 - t214;
t203 = qJD(1) * t152;
t241 = qJD(3) + qJD(4);
t182 = t235 * qJD(4);
t242 = t235 * qJD(3) + t182;
t228 = -t242 * t151 + t162 * t203 + t241 * t214;
t102 = t148 * t151 + t235 * t149;
t66 = t241 * t102;
t227 = -t102 * t203 + t66;
t254 = t164 * t55;
t253 = t55 * qJ(5);
t193 = qJD(1) * qJD(2);
t181 = t152 * t193;
t252 = qJD(2) * qJD(3) + t181;
t199 = qJD(2) * t152;
t185 = t149 * t199;
t196 = qJD(3) * t151;
t186 = t150 * t196;
t251 = t185 + t186;
t133 = -qJD(3) + t203;
t123 = -qJD(4) + t133;
t197 = qJD(3) * t150;
t180 = qJD(1) * t197;
t190 = t149 * t252 + t151 * t180;
t195 = qJD(4) * t148;
t71 = -t149 * t180 + t151 * t252;
t21 = t100 * t195 + t148 * t190 - t99 * t182 - t235 * t71;
t250 = -t123 * t55 - t21;
t136 = t150 * t193;
t170 = pkin(6) * t136;
t109 = t168 * qJD(2);
t93 = qJD(1) * t109;
t222 = -t149 * t170 - t151 * t93;
t143 = pkin(6) * t203;
t117 = qJD(2) * pkin(7) + t143;
t111 = -pkin(2) * t152 - pkin(7) * t150 - pkin(1);
t92 = t111 * qJD(1);
t224 = t149 * t92;
t63 = t117 * t151 + t224;
t160 = -qJD(3) * t63 - t222;
t19 = pkin(3) * t136 - pkin(8) * t71 + t160;
t198 = qJD(3) * t149;
t166 = -t117 * t198 + t149 * t93 + t92 * t196;
t159 = -t151 * t170 + t166;
t24 = -pkin(8) * t190 + t159;
t62 = -t117 * t149 + t151 * t92;
t41 = -pkin(8) * t100 + t62;
t36 = -pkin(3) * t133 + t41;
t42 = pkin(8) * t99 + t63;
t175 = -t148 * t19 - t36 * t182 + t42 * t195 - t235 * t24;
t116 = -qJD(2) * pkin(2) + pkin(6) * t204;
t72 = -pkin(3) * t99 + t116;
t249 = t55 * t72 + t175;
t248 = -0.2e1 * t193;
t33 = pkin(4) * t55 + qJD(5) + t72;
t247 = t164 * t33;
t246 = -t143 + (-t149 * t203 + t198) * pkin(3);
t245 = qJ(5) * t164;
t118 = t237 * t149;
t119 = t237 * t151;
t206 = -t148 * t118 + t235 * t119;
t244 = t206 * qJD(4) + t257 * t148 + t256 * t235;
t243 = -t118 * t182 - t119 * t195 - t256 * t148 + t257 * t235;
t40 = t235 * t42;
t12 = t148 * t36 + t40;
t158 = -qJD(4) * t12 - t148 * t24 + t235 * t19;
t240 = -t72 * t164 + t158;
t22 = qJD(4) * t164 + t148 * t71 + t235 * t190;
t239 = -t123 * t164 - t22;
t38 = t148 * t42;
t11 = t235 * t36 - t38;
t6 = t11 - t245;
t5 = -pkin(4) * t123 + t6;
t236 = t5 - t6;
t234 = pkin(6) * t149;
t233 = -qJ(5) * t227 + qJD(5) * t162 + t243;
t232 = -pkin(4) * t204 + qJ(5) * t228 - t102 * qJD(5) - t244;
t231 = t235 * t41 - t38;
t98 = t151 * t111;
t61 = -pkin(8) * t211 + t98 + (-pkin(3) - t234) * t152;
t213 = t149 * t150;
t135 = pkin(6) * t210;
t219 = t149 * t111 + t135;
t67 = -pkin(8) * t213 + t219;
t229 = t148 * t61 + t235 * t67;
t226 = t149 * t109 + t111 * t196;
t225 = t133 * t99;
t223 = t71 * t149;
t201 = qJD(2) * t150;
t220 = t151 * t109 + t201 * t234;
t218 = t100 * t133;
t217 = t116 * t149;
t216 = t116 * t151;
t215 = t133 * t151;
t154 = qJD(1) ^ 2;
t209 = t152 * t154;
t153 = qJD(2) ^ 2;
t208 = t153 * t150;
t207 = t153 * t152;
t110 = pkin(3) * t213 + t150 * pkin(6);
t146 = t150 ^ 2;
t205 = -t152 ^ 2 + t146;
t194 = t116 * qJD(3);
t73 = pkin(3) * t251 + pkin(6) * t199;
t141 = -pkin(3) * t151 - pkin(2);
t188 = t149 * t197;
t187 = t152 * t198;
t179 = -t148 * t41 - t40;
t177 = -t148 * t67 + t235 * t61;
t174 = pkin(1) * t248;
t173 = -t99 + t200;
t172 = -t235 * t118 - t119 * t148;
t171 = -t100 + t202;
t169 = t235 * t199;
t167 = qJD(1) * t146 - t133 * t152;
t52 = pkin(3) * t190 + pkin(6) * t181;
t28 = t165 * qJD(2) + (-t135 + (pkin(8) * t150 - t111) * t149) * qJD(3) + t220;
t30 = -t251 * pkin(8) + (-t150 * t200 - t187) * pkin(6) + t226;
t163 = t148 * t28 + t61 * t182 - t195 * t67 + t235 * t30;
t10 = t22 * pkin(4) + t52;
t157 = -t229 * qJD(4) - t148 * t30 + t235 * t28;
t140 = t235 * pkin(3) + pkin(4);
t81 = t162 * t150;
t80 = t102 * t150;
t46 = qJ(5) * t162 + t206;
t45 = -qJ(5) * t102 + t172;
t32 = t149 * t169 - t148 * t188 - t195 * t213 + (t148 * t199 + t242 * t150) * t151;
t31 = t148 * t185 + t150 * t66 - t151 * t169;
t25 = -qJ(5) * t80 + t229;
t23 = -pkin(4) * t152 - qJ(5) * t81 + t177;
t9 = t231 - t245;
t8 = t179 + t253;
t7 = t12 - t253;
t4 = -qJ(5) * t32 - qJD(5) * t80 + t163;
t3 = pkin(4) * t201 + t31 * qJ(5) - t81 * qJD(5) + t157;
t2 = -qJ(5) * t22 - qJD(5) * t55 - t175;
t1 = pkin(4) * t136 + t21 * qJ(5) - qJD(5) * t164 + t158;
t13 = [0, 0, 0, 0.2e1 * t152 * t136, t205 * t248, t207, -t208, 0, -pkin(6) * t207 + t150 * t174, pkin(6) * t208 + t152 * t174, t71 * t211 + (t151 * t199 - t188) * t100, (-t100 * t149 + t151 * t99) * t199 + (-t151 * t190 - t223 + (-t100 * t151 - t149 * t99) * qJD(3)) * t150, t133 * t188 - t152 * t71 + (t100 * t150 + t151 * t167) * qJD(2), t133 * t186 + t190 * t152 + (-t149 * t167 + t99 * t150) * qJD(2), (-t133 - t203) * t201, -(-t111 * t198 + t220) * t133 + (pkin(6) * t190 + t151 * t194 + (t98 * qJD(1) + t62) * qJD(2)) * t150 + ((-pkin(6) * t99 + t217) * qJD(2) + (t224 + (pkin(6) * t133 + t117) * t151) * qJD(3) + t222) * t152, (-pkin(6) * t187 + t226) * t133 + t166 * t152 + (pkin(6) * t71 - t149 * t194) * t150 + ((pkin(6) * t100 + t216) * t152 + (-pkin(6) * t215 - qJD(1) * t219 - t63) * t150) * qJD(2), -t164 * t31 - t21 * t81, -t164 * t32 + t21 * t80 - t22 * t81 + t31 * t55, t123 * t31 + t152 * t21 + (qJD(1) * t81 + t164) * t201, t123 * t32 + t152 * t22 + (-qJD(1) * t80 - t55) * t201, (-t123 - t203) * t201, t11 * t201 + t110 * t22 - t123 * t157 + t136 * t177 - t152 * t158 + t72 * t32 + t52 * t80 + t73 * t55, t163 * t123 - t175 * t152 + t73 * t164 - t110 * t21 + t52 * t81 - t72 * t31 + (-t229 * qJD(1) - t12) * t201, -t1 * t81 - t164 * t3 - t2 * t80 + t21 * t23 - t22 * t25 + t31 * t5 - t32 * t7 - t4 * t55, t2 * t25 + t7 * t4 + t1 * t23 + t5 * t3 + t10 * (pkin(4) * t80 + t110) + t33 * (pkin(4) * t32 + t73); 0, 0, 0, -t150 * t209, t205 * t154, 0, 0, 0, t154 * pkin(1) * t150, pkin(1) * t209, -t100 * t215 + t223, (t71 - t225) * t151 + (-t190 + t218) * t149, -t133 * t196 + (t133 * t210 + t150 * t171) * qJD(1), t133 * t198 + (-t133 * t212 + t150 * t173) * qJD(1), t133 * t204, -pkin(2) * t190 + t221 * t133 + (pkin(7) * t215 + t217) * qJD(3) + ((-pkin(7) * t202 - t62) * t150 + (-pkin(6) * t173 - t217) * t152) * qJD(1), -pkin(2) * t71 - t86 * t133 + (-t149 * pkin(7) * t133 + t216) * qJD(3) + (-t116 * t210 + (-pkin(7) * t200 + t63) * t150 + (t133 * t211 + t152 * t171) * pkin(6)) * qJD(1), -t102 * t21 - t164 * t228, -t102 * t22 - t162 * t21 - t164 * t227 + t228 * t55, t228 * t123 + (qJD(2) * t102 - t164) * t204, t227 * t123 + (qJD(2) * t162 + t55) * t204, t123 * t204, -t11 * t204 + t244 * t123 + t136 * t172 + t141 * t22 - t162 * t52 + t227 * t72 + t246 * t55, t52 * t102 - t141 * t21 - t228 * t72 + t246 * t164 + t243 * t123 + (-qJD(2) * t206 + t12) * t204, -t1 * t102 + t162 * t2 - t164 * t232 + t21 * t45 - t22 * t46 - t227 * t7 + t228 * t5 - t233 * t55, t2 * t46 + t1 * t45 + t10 * (-pkin(4) * t162 + t141) + t233 * t7 + t232 * t5 + (t227 * pkin(4) + t246) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100 * t99, t100 ^ 2 - t99 ^ 2, t71 + t225, -t190 - t218, t136, -t100 * t116 - t133 * t63 + t160, -t116 * t99 - t133 * t62 - t159, t254, t255, t250, t239, t136, t179 * t123 + (-t100 * t55 + t123 * t195 + t235 * t136) * pkin(3) + t240, -t231 * t123 + (-t100 * t164 + t123 * t182 - t136 * t148) * pkin(3) + t249, t140 * t21 - t5 * t55 + t7 * t164 + t9 * t55 + t8 * t164 + (-t148 * t22 + (t148 * t164 - t235 * t55) * qJD(4)) * pkin(3), -pkin(4) * t247 + t1 * t140 - t5 * t8 - t7 * t9 + (-t33 * t100 + t2 * t148 + (-t148 * t5 + t235 * t7) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t254, t255, t250, t239, t136, -t12 * t123 + t240, -t11 * t123 + t249, pkin(4) * t21 - t236 * t55, t236 * t7 + (t1 - t247) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53 - t238, t164 * t5 + t7 * t55 + t10;];
tauc_reg = t13;
