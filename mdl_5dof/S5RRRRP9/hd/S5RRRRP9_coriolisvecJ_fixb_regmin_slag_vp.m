% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRRP9
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
% tauc_reg [5x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRP9_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP9_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:06:04
% EndTime: 2019-12-31 22:06:14
% DurationCPUTime: 3.00s
% Computational Cost: add. (3958->362), mult. (9685->491), div. (0->0), fcn. (6413->6), ass. (0->178)
t148 = sin(qJ(3));
t245 = -pkin(8) - pkin(7);
t192 = qJD(3) * t245;
t149 = sin(qJ(2));
t150 = cos(qJ(3));
t216 = t149 * t150;
t151 = cos(qJ(2));
t217 = t148 * t151;
t172 = pkin(2) * t149 - pkin(7) * t151;
t100 = t172 * qJD(1);
t86 = t148 * t100;
t262 = t148 * t192 - t86 - (-pkin(6) * t216 - pkin(8) * t217) * qJD(1);
t214 = t150 * t151;
t169 = pkin(3) * t149 - pkin(8) * t214;
t207 = qJD(1) * t149;
t187 = t148 * t207;
t228 = pkin(6) * t187 + t150 * t100;
t261 = -qJD(1) * t169 + t150 * t192 - t228;
t147 = sin(qJ(4));
t244 = cos(qJ(4));
t199 = t150 * qJD(2);
t96 = -t187 + t199;
t200 = t148 * qJD(2);
t97 = t150 * t207 + t200;
t162 = t147 * t96 + t244 * t97;
t55 = t147 * t97 - t244 * t96;
t111 = -qJD(2) * pkin(2) + pkin(6) * t207;
t73 = -t96 * pkin(3) + t111;
t17 = t55 * pkin(4) - qJ(5) * t162 + t73;
t260 = t17 * t55;
t259 = t73 * t55;
t219 = t147 * t148;
t159 = t244 * t150 - t219;
t198 = t151 * qJD(1);
t249 = qJD(3) + qJD(4);
t185 = t244 * qJD(4);
t250 = t244 * qJD(3) + t185;
t236 = -t250 * t150 + t159 * t198 + t249 * t219;
t99 = t147 * t150 + t244 * t148;
t66 = t249 * t99;
t235 = -t99 * t198 + t66;
t241 = t162 * t55;
t197 = qJD(1) * qJD(2);
t184 = t151 * t197;
t258 = qJD(2) * qJD(3) + t184;
t188 = t151 * t200;
t202 = qJD(3) * t150;
t189 = t149 * t202;
t257 = t188 + t189;
t246 = t162 ^ 2;
t256 = -t55 ^ 2 + t246;
t131 = -qJD(3) + t198;
t118 = -qJD(4) + t131;
t203 = qJD(3) * t149;
t183 = qJD(1) * t203;
t209 = t258 * t150;
t158 = t148 * t183 - t209;
t193 = t258 * t148 + t150 * t183;
t201 = qJD(4) * t147;
t20 = t147 * t193 + t244 * t158 - t96 * t185 + t97 * t201;
t7 = -t55 * t118 - t20;
t30 = pkin(4) * t162 + t55 * qJ(5);
t255 = -0.2e1 * t197;
t242 = t17 * t162;
t253 = t73 * t162;
t113 = t245 * t148;
t114 = t245 * t150;
t160 = t244 * t113 + t147 * t114;
t252 = t160 * qJD(4) + t261 * t147 + t262 * t244;
t72 = t147 * t113 - t244 * t114;
t251 = t72 * qJD(4) + t262 * t147 - t261 * t244;
t142 = pkin(6) * t198;
t204 = qJD(3) * t148;
t174 = -t142 + (-t148 * t198 + t204) * pkin(3);
t21 = qJD(4) * t162 - t147 * t158 + t244 * t193;
t248 = -t118 * t162 - t21;
t106 = -t151 * pkin(2) - t149 * pkin(7) - pkin(1);
t95 = t150 * t106;
t61 = -pkin(8) * t216 + t95 + (-pkin(6) * t148 - pkin(3)) * t151;
t218 = t148 * t149;
t133 = pkin(6) * t214;
t225 = t148 * t106 + t133;
t67 = -pkin(8) * t218 + t225;
t163 = t147 * t61 + t244 * t67;
t103 = t172 * qJD(2);
t144 = t149 * pkin(6);
t227 = t150 * t103 + t200 * t144;
t31 = t169 * qJD(2) + (-t133 + (pkin(8) * t149 - t106) * t148) * qJD(3) + t227;
t190 = t151 * t204;
t234 = t148 * t103 + t106 * t202;
t33 = -t257 * pkin(8) + (-t149 * t199 - t190) * pkin(6) + t234;
t247 = -qJD(4) * t163 - t147 * t33 + t244 * t31;
t243 = pkin(6) * t131;
t240 = t235 * pkin(4) + t236 * qJ(5) - t99 * qJD(5) + t174;
t239 = qJ(5) * t207 - t252;
t238 = pkin(4) * t207 + t251;
t233 = qJD(2) * pkin(7);
t112 = t142 + t233;
t215 = t150 * t112;
t90 = t106 * qJD(1);
t231 = t148 * t90;
t63 = t215 + t231;
t46 = t96 * pkin(8) + t63;
t232 = t147 * t46;
t230 = t96 * t131;
t135 = t149 * t197;
t177 = pkin(6) * t135;
t91 = qJD(1) * t103;
t229 = t148 * t177 + t150 * t91;
t62 = -t148 * t112 + t150 * t90;
t45 = -t97 * pkin(8) + t62;
t13 = t244 * t45 - t232;
t226 = pkin(3) * t185 + qJD(5) - t13;
t224 = qJD(2) * t160;
t223 = qJD(2) * t72;
t222 = t111 * t148;
t221 = t111 * t150;
t220 = t131 * t150;
t153 = qJD(1) ^ 2;
t213 = t151 * t153;
t152 = qJD(2) ^ 2;
t212 = t152 * t149;
t211 = t152 * t151;
t41 = -t131 * pkin(3) + t45;
t10 = t244 * t41 - t232;
t210 = qJD(5) - t10;
t104 = pkin(3) * t218 + t144;
t145 = t149 ^ 2;
t208 = -t151 ^ 2 + t145;
t206 = qJD(2) * t149;
t205 = qJD(2) * t151;
t194 = t244 * t46;
t74 = t257 * pkin(3) + pkin(6) * t205;
t140 = -t150 * pkin(3) - pkin(2);
t191 = t148 * t203;
t16 = -t209 * pkin(8) + pkin(3) * t135 + (-t215 + (pkin(8) * t207 - t90) * t148) * qJD(3) + t229;
t170 = -t112 * t204 + t148 * t91 + t90 * t202;
t155 = -t150 * t177 + t170;
t22 = -t193 * pkin(8) + t155;
t182 = -t147 * t16 - t41 * t185 + t46 * t201 - t244 * t22;
t181 = t147 * t22 - t244 * t16 + t46 * t185 + t41 * t201;
t180 = pkin(1) * t255;
t179 = -t97 + t200;
t178 = -t96 + t199;
t176 = pkin(4) * t135;
t175 = t244 * t205;
t12 = t147 * t45 + t194;
t173 = pkin(3) * t201 - t12;
t171 = qJD(1) * t145 - t131 * t151;
t105 = t118 * qJD(5);
t126 = qJ(5) * t135;
t1 = t126 - t105 - t182;
t53 = t193 * pkin(3) + pkin(6) * t184;
t168 = -t10 * t118 + t182;
t11 = t147 * t41 + t194;
t167 = -t11 * t118 - t181;
t165 = -t147 * t67 + t244 * t61;
t161 = t147 * t31 + t61 * t185 - t67 * t201 + t244 * t33;
t2 = -t176 + t181;
t157 = t158 * t149;
t139 = -t244 * pkin(3) - pkin(4);
t134 = t147 * pkin(3) + qJ(5);
t82 = t159 * t149;
t81 = t99 * t149;
t52 = -pkin(4) * t159 - t99 * qJ(5) + t140;
t42 = t81 * pkin(4) - t82 * qJ(5) + t104;
t35 = t148 * t175 - t147 * t191 - t201 * t218 + (t147 * t205 + t250 * t149) * t150;
t34 = t147 * t188 + t149 * t66 - t150 * t175;
t29 = t151 * pkin(4) - t165;
t28 = -t151 * qJ(5) + t163;
t24 = t97 * pkin(3) + t30;
t9 = -t118 * qJ(5) + t11;
t8 = t118 * pkin(4) + t210;
t6 = t35 * pkin(4) + t34 * qJ(5) - t82 * qJD(5) + t74;
t5 = -pkin(4) * t206 - t247;
t4 = qJ(5) * t206 - t151 * qJD(5) + t161;
t3 = t21 * pkin(4) + t20 * qJ(5) - qJD(5) * t162 + t53;
t14 = [0, 0, 0, 0.2e1 * t151 * t135, t208 * t255, t211, -t212, 0, -pkin(6) * t211 + t149 * t180, pkin(6) * t212 + t151 * t180, -t97 * t191 + (t97 * t205 - t157) * t150, (-t97 * t148 + t150 * t96) * t205 + (-t150 * t193 - t209 * t148 + (-t97 * t150 + (-t96 + t187) * t148) * qJD(3)) * t149, t131 * t191 + t158 * t151 + (t97 * t149 + t150 * t171) * qJD(2), t131 * t189 + t193 * t151 + (-t148 * t171 + t96 * t149) * qJD(2), (-t131 - t198) * t206, -(-t106 * t204 + t227) * t131 + (pkin(6) * t193 + t111 * t202 + (t95 * qJD(1) + t62) * qJD(2)) * t149 + ((-pkin(6) * t96 + t222) * qJD(2) + (t231 + (t112 + t243) * t150) * qJD(3) - t229) * t151, -pkin(6) * t157 - t111 * t191 + (-pkin(6) * t190 + t234) * t131 + t170 * t151 + ((pkin(6) * t97 + t221) * t151 + (-pkin(6) * t220 - t225 * qJD(1) - t63) * t149) * qJD(2), -t162 * t34 - t20 * t82, -t162 * t35 + t20 * t81 - t82 * t21 + t34 * t55, t34 * t118 + t20 * t151 + (qJD(1) * t82 + t162) * t206, t35 * t118 + t21 * t151 + (-qJD(1) * t81 - t55) * t206, (-t118 - t198) * t206, t10 * t206 + t104 * t21 - t247 * t118 + t165 * t135 + t181 * t151 + t73 * t35 + t53 * t81 + t74 * t55, t161 * t118 - t182 * t151 + t74 * t162 - t104 * t20 + t53 * t82 - t73 * t34 + (-t163 * qJD(1) - t11) * t206, t5 * t118 + t2 * t151 + t17 * t35 + t42 * t21 + t3 * t81 + t6 * t55 + (-qJD(1) * t29 - t8) * t206, -t1 * t81 + t162 * t5 + t2 * t82 - t29 * t20 - t28 * t21 - t8 * t34 - t9 * t35 - t4 * t55, -t1 * t151 - t4 * t118 + t17 * t34 + t42 * t20 - t3 * t82 - t6 * t162 + (qJD(1) * t28 + t9) * t206, t1 * t28 + t17 * t6 + t2 * t29 + t3 * t42 + t9 * t4 + t8 * t5; 0, 0, 0, -t149 * t213, t208 * t153, 0, 0, 0, t153 * pkin(1) * t149, pkin(1) * t213, -t148 * t158 - t97 * t220, (t209 - t230) * t150 + (-t97 * qJD(3) + (t151 * t97 - t189) * qJD(1) - t193) * t148, -t131 * t202 + (t131 * t214 + t149 * t179) * qJD(1), t131 * t204 + (-t131 * t217 + t149 * t178) * qJD(1), t131 * t207, -pkin(2) * t193 + t228 * t131 + (pkin(7) * t220 + t222) * qJD(3) + ((-pkin(7) * t200 - t62) * t149 + (-pkin(6) * t178 - t222) * t151) * qJD(1), -pkin(2) * t209 - t86 * t131 + (-t148 * pkin(7) * t131 + t221) * qJD(3) + ((pkin(6) * t179 - t221) * t151 + (pkin(2) * t204 + t63 + (-t233 + t243) * t150) * t149) * qJD(1), -t162 * t236 - t20 * t99, -t159 * t20 - t162 * t235 - t99 * t21 + t236 * t55, t236 * t118 + (qJD(2) * t99 - t162) * t207, t235 * t118 + (qJD(2) * t159 + t55) * t207, t118 * t207, t140 * t21 - t53 * t159 + t235 * t73 + t174 * t55 + t251 * t118 + (-t10 + t224) * t207, -t140 * t20 + t53 * t99 - t236 * t73 + t174 * t162 + t252 * t118 + (t11 - t223) * t207, t52 * t21 - t3 * t159 + t240 * t55 + t235 * t17 + t238 * t118 + (t8 + t224) * t207, t1 * t159 + t160 * t20 + t162 * t238 + t2 * t99 - t72 * t21 - t235 * t9 - t236 * t8 + t239 * t55, t52 * t20 - t3 * t99 - t240 * t162 + t236 * t17 + t239 * t118 + (-t9 + t223) * t207, t1 * t72 - t160 * t2 + t17 * t240 + t238 * t8 - t239 * t9 + t3 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97 * t96, -t96 ^ 2 + t97 ^ 2, -t158 + t230, -t97 * t131 - t193, t135, -t111 * t97 + t229 + (-qJD(3) - t131) * t63, -t111 * t96 - t62 * t131 - t155, t241, t256, t7, t248, t135, -t12 * t118 - t253 + (t118 * t201 + t244 * t135 - t55 * t97) * pkin(3) - t181, -t13 * t118 + t259 + (t118 * t185 - t135 * t147 - t162 * t97) * pkin(3) + t182, -t242 - t24 * t55 + t173 * t118 + (pkin(4) - t139) * t135 - t181, -t134 * t21 - t139 * t20 + (t173 + t9) * t162 + (-t226 + t8) * t55, -t118 * t226 + t134 * t135 + t162 * t24 + t1 - t260, t1 * t134 + t2 * t139 - t17 * t24 + t173 * t8 + t226 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t241, t256, t7, t248, t135, t167 - t253, t168 + t259, -t30 * t55 + t167 + 0.2e1 * t176 - t242, pkin(4) * t20 - t21 * qJ(5) + (-t11 + t9) * t162 + (t8 - t210) * t55, t162 * t30 - 0.2e1 * t105 + 0.2e1 * t126 - t168 - t260, -t2 * pkin(4) + t1 * qJ(5) - t8 * t11 - t17 * t30 + t210 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t135 + t241, t7, -t118 ^ 2 - t246, t9 * t118 + t2 + t242;];
tauc_reg = t14;
