% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% tauc_reg [5x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRP11_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP11_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:18:31
% EndTime: 2019-12-31 22:18:47
% DurationCPUTime: 5.14s
% Computational Cost: add. (5744->435), mult. (15320->606), div. (0->0), fcn. (11481->8), ass. (0->192)
t150 = sin(pkin(5));
t157 = cos(qJ(2));
t223 = qJD(1) * t157;
t142 = t150 * t223;
t177 = t142 - qJD(3);
t154 = sin(qJ(2));
t151 = cos(pkin(5));
t225 = qJD(1) * t151;
t211 = pkin(1) * t225;
t141 = t154 * t211;
t119 = pkin(7) * t142 + t141;
t153 = sin(qJ(3));
t156 = cos(qJ(3));
t269 = -t119 - t177 * (pkin(3) * t153 - pkin(9) * t156);
t220 = qJD(3) * t153;
t270 = pkin(8) * t220;
t191 = qJD(2) + t225;
t224 = qJD(1) * t154;
t206 = t150 * t224;
t160 = -t153 * t206 + t156 * t191;
t214 = qJD(1) * qJD(2);
t201 = t150 * t214;
t184 = t157 * t201;
t171 = t156 * t184;
t159 = t160 * qJD(3) + t171;
t268 = -qJD(4) * t177 + t159;
t267 = t153 * t142 - t220;
t96 = qJD(4) - t160;
t152 = sin(qJ(4));
t155 = cos(qJ(4));
t102 = t153 * t191 + t156 * t206;
t116 = -pkin(7) * t206 + t157 * t211;
t82 = -t191 * pkin(2) - t116;
t44 = -pkin(3) * t160 - t102 * pkin(9) + t82;
t83 = t191 * pkin(8) + t119;
t114 = (-pkin(2) * t157 - pkin(8) * t154 - pkin(1)) * t150;
t95 = qJD(1) * t114;
t52 = t153 * t95 + t156 * t83;
t46 = -t177 * pkin(9) + t52;
t14 = t152 * t44 + t155 * t46;
t10 = qJ(5) * t96 + t14;
t166 = (pkin(2) * t154 - pkin(8) * t157) * t150;
t118 = qJD(2) * t166;
t109 = qJD(1) * t118;
t234 = t150 * t154;
t143 = pkin(7) * t234;
t255 = pkin(1) * t157;
t120 = (t151 * t255 - t143) * qJD(2);
t110 = qJD(1) * t120;
t218 = qJD(3) * t156;
t165 = -t153 * t109 - t156 * t110 - t95 * t218 + t83 * t220;
t185 = t154 * t201;
t21 = pkin(9) * t185 - t165;
t216 = qJD(4) * t155;
t217 = qJD(4) * t152;
t173 = pkin(7) * t184;
t172 = t153 * t184;
t71 = t102 * qJD(3) + t172;
t33 = t71 * pkin(3) - t159 * pkin(9) + qJD(2) * t141 + t173;
t200 = t152 * t21 - t155 * t33 + t46 * t216 + t44 * t217;
t258 = pkin(4) * t71;
t2 = t200 - t258;
t266 = -t10 * t96 + t2;
t265 = t154 * t157;
t187 = t156 * t142;
t264 = t187 - t218;
t138 = -pkin(3) * t156 - pkin(9) * t153 - pkin(2);
t117 = qJD(1) * t166;
t230 = t156 * t116 + t153 * t117;
t60 = pkin(9) * t206 + t230;
t263 = -t138 * t216 - t269 * t152 + t155 * t60;
t51 = -t153 * t83 + t156 * t95;
t45 = t177 * pkin(3) - t51;
t67 = t102 * t152 + t155 * t177;
t69 = t155 * t102 - t152 * t177;
t12 = t67 * pkin(4) - t69 * qJ(5) + t45;
t257 = pkin(9) * t71;
t262 = t12 * t96 - t257;
t112 = t143 + (-pkin(2) - t255) * t151;
t123 = -t151 * t156 + t153 * t234;
t124 = t151 * t153 + t156 * t234;
t55 = pkin(3) * t123 - pkin(9) * t124 + t112;
t233 = t150 * t157;
t213 = pkin(7) * t233;
t256 = pkin(1) * t154;
t113 = t213 + (pkin(8) + t256) * t151;
t243 = t156 * t113 + t153 * t114;
t57 = -pkin(9) * t233 + t243;
t174 = t152 * t55 + t155 * t57;
t162 = -t113 * t220 + t114 * t218 + t153 * t118 + t156 * t120;
t222 = qJD(2) * t154;
t205 = t150 * t222;
t29 = pkin(9) * t205 + t162;
t121 = (t151 * t256 + t213) * qJD(2);
t204 = qJD(2) * t233;
t75 = t124 * qJD(3) + t153 * t204;
t76 = -t123 * qJD(3) + t156 * t204;
t38 = pkin(3) * t75 - pkin(9) * t76 + t121;
t261 = -t174 * qJD(4) - t152 * t29 + t155 * t38;
t260 = t69 ^ 2;
t259 = t96 ^ 2;
t158 = qJD(1) ^ 2;
t254 = pkin(8) * t156;
t252 = t12 * t69;
t251 = t67 * t96;
t250 = t69 * t67;
t249 = t69 * t96;
t215 = qJD(4) * t156;
t219 = qJD(3) * t155;
t248 = qJD(5) * t156 - (-t152 * t215 - t153 * t219) * pkin(8) + t263 + t267 * qJ(5);
t247 = -t138 * t217 + (t60 + t270) * t152 + (-t215 * pkin(8) + t269) * t155 - t267 * pkin(4);
t178 = pkin(4) * t152 - qJ(5) * t155;
t169 = pkin(8) + t178;
t179 = pkin(4) * t155 + qJ(5) * t152;
t103 = t153 * t116;
t59 = -pkin(3) * t206 - t117 * t156 + t103;
t89 = t152 * t187 - t155 * t206;
t226 = qJD(1) * t150;
t232 = t155 * t157;
t90 = (t152 * t154 + t156 * t232) * t226;
t246 = pkin(4) * t89 - qJ(5) * t90 + t59 - (t179 * qJD(4) - qJD(5) * t155) * t153 - t169 * t218;
t63 = pkin(3) * t102 - pkin(9) * t160;
t245 = t152 * t63 + t155 * t51;
t242 = pkin(9) * qJD(4);
t241 = qJ(5) * t71;
t34 = t102 * t217 - t152 * t185 - t268 * t155;
t240 = t152 * t34;
t239 = t152 * t71;
t238 = t155 * t71;
t237 = -qJD(5) * t152 + t96 * t178 - t52;
t236 = t138 * t155;
t147 = t150 ^ 2;
t235 = t147 * t158;
t13 = -t152 * t46 + t155 * t44;
t231 = qJD(5) - t13;
t228 = t152 * t138 + t155 * t254;
t227 = t154 ^ 2 - t157 ^ 2;
t221 = qJD(3) * t152;
t212 = t96 * t242;
t210 = t96 * t221;
t209 = t96 * t219;
t208 = t96 * t217;
t207 = t152 * t233;
t203 = t150 * t151 * t158;
t202 = t147 * t214;
t199 = t156 * t109 - t153 * t110 - t83 * t218 - t95 * t220;
t198 = -t153 * t113 + t114 * t156;
t197 = t155 * t96;
t196 = t157 * t177;
t195 = t177 * t150;
t194 = qJD(3) * t177;
t192 = 0.2e1 * t202;
t190 = qJD(2) + 0.2e1 * t225;
t22 = -pkin(3) * t185 - t199;
t35 = t102 * t216 + t268 * t152 - t155 * t185;
t6 = pkin(4) * t35 + qJ(5) * t34 - qJD(5) * t69 + t22;
t189 = -t6 - t212;
t56 = pkin(3) * t233 - t198;
t183 = -0.2e1 * pkin(1) * t202;
t181 = -t155 * t218 + t90;
t9 = -pkin(4) * t96 + t231;
t180 = -t10 * t152 + t155 * t9;
t175 = -t152 * t57 + t155 * t55;
t170 = -t113 * t218 - t114 * t220 + t118 * t156 - t153 * t120;
t168 = t14 * t96 - t200;
t167 = -t96 * t216 - t239;
t77 = t124 * t152 + t150 * t232;
t3 = t152 * t33 + t155 * t21 + t44 * t216 - t46 * t217;
t164 = t152 * t38 + t155 * t29 + t55 * t216 - t57 * t217;
t163 = t96 * t45 - t257;
t161 = pkin(1) * (-t151 * t214 + t235);
t30 = -pkin(3) * t205 - t170;
t134 = -pkin(3) - t179;
t115 = t169 * t153;
t111 = qJD(1) * t121;
t85 = -t236 + (pkin(8) * t152 + pkin(4)) * t156;
t84 = -qJ(5) * t156 + t228;
t78 = t124 * t155 - t207;
t42 = -t77 * qJD(4) + t152 * t205 + t155 * t76;
t41 = -qJD(4) * t207 + t124 * t216 + t152 * t76 - t155 * t205;
t37 = pkin(4) * t69 + qJ(5) * t67;
t23 = pkin(4) * t77 - qJ(5) * t78 + t56;
t18 = -pkin(4) * t123 - t175;
t17 = qJ(5) * t123 + t174;
t16 = -pkin(4) * t102 + t152 * t51 - t155 * t63;
t15 = qJ(5) * t102 + t245;
t11 = -t34 + t251;
t8 = pkin(4) * t41 - qJ(5) * t42 - qJD(5) * t78 + t30;
t7 = -pkin(4) * t75 - t261;
t5 = qJ(5) * t75 + qJD(5) * t123 + t164;
t1 = qJD(5) * t96 + t241 + t3;
t4 = [0, 0, 0, t192 * t265, -t227 * t192, t190 * t204, -t190 * t205, 0, -t111 * t151 - t121 * t191 + t154 * t183, -t110 * t151 - t120 * t191 + t157 * t183, t102 * t76 + t124 * t159, -t102 * t75 - t123 * t159 - t124 * t71 + t160 * t76, t102 * t205 + t124 * t185 - t159 * t233 - t76 * t177, t75 * t177 + (t71 * t157 + (-qJD(1) * t123 + t160) * t222) * t150, (-t147 * t223 - t195) * t222, -t170 * t177 - t121 * t160 + t112 * t71 + t111 * t123 + t82 * t75 + (-t199 * t157 + (t198 * qJD(1) + t51) * t222) * t150, t121 * t102 + t111 * t124 + t112 * t159 + t162 * t177 - t165 * t233 - t243 * t185 - t52 * t205 + t82 * t76, -t34 * t78 + t42 * t69, t34 * t77 - t35 * t78 - t41 * t69 - t42 * t67, -t123 * t34 + t42 * t96 + t69 * t75 + t71 * t78, -t123 * t35 - t41 * t96 - t67 * t75 - t71 * t77, t123 * t71 + t75 * t96, -t123 * t200 + t13 * t75 + t175 * t71 + t22 * t77 + t261 * t96 + t30 * t67 + t56 * t35 + t45 * t41, -t3 * t123 - t14 * t75 - t164 * t96 - t174 * t71 + t22 * t78 + t30 * t69 - t56 * t34 + t45 * t42, t12 * t41 - t123 * t2 - t18 * t71 + t23 * t35 + t6 * t77 + t67 * t8 - t7 * t96 - t75 * t9, -t1 * t77 - t10 * t41 - t17 * t35 - t18 * t34 + t2 * t78 + t42 * t9 - t5 * t67 + t69 * t7, t1 * t123 + t10 * t75 - t12 * t42 + t17 * t71 + t23 * t34 + t5 * t96 - t6 * t78 - t69 * t8, t1 * t17 + t10 * t5 + t12 * t8 + t18 * t2 + t23 * t6 + t7 * t9; 0, 0, 0, -t235 * t265, t227 * t235, -t157 * t203, t154 * t203, 0, t119 * t191 + t154 * t161 - t173, pkin(7) * t185 + t116 * t191 + t157 * t161, -qJD(3) * t153 ^ 2 * t206 + ((qJD(3) * t191 + t184) * t153 - t177 * t102) * t156, t267 * t102 - t153 * t71 + t156 * t159 - t160 * t264, -t156 * t194 + (t156 * t196 + (t153 * qJD(2) - t102) * t154) * t226, t153 * t194 + (-t153 * t196 + (t156 * qJD(2) - t160) * t154) * t226, t195 * t224, -pkin(2) * t71 + t82 * t220 - t103 * t177 + t119 * t160 + (pkin(8) * t194 + t117 * t177 - t111) * t156 + (-t51 * t154 + (-pkin(8) * t222 - t82 * t157) * t153) * t226, -pkin(2) * t159 - t119 * t102 + t111 * t153 - t185 * t254 + t52 * t206 - t264 * t82 + (-t230 - t270) * t177, -t153 * t155 * t34 + (-t153 * t217 - t181) * t69, t67 * t90 + t69 * t89 + (-t152 * t69 - t155 * t67) * t218 + (t240 - t155 * t35 + (t152 * t67 - t155 * t69) * qJD(4)) * t153, -t90 * t96 + (t34 + t209) * t156 + (-t177 * t69 - t208 + t238) * t153, t89 * t96 + (t35 - t210) * t156 + (t177 * t67 + t167) * t153, -t153 * t177 * t96 - t156 * t71, t71 * t236 - t45 * t89 - t59 * t67 + (t269 * t155 + (-qJD(4) * t138 + t60) * t152) * t96 + (t45 * t221 + t200 + (qJD(3) * t67 + t167) * pkin(8)) * t156 + (t45 * t216 + t22 * t152 - t177 * t13 + (t35 + t210) * pkin(8)) * t153, -t228 * t71 - t59 * t69 - t45 * t90 + t263 * t96 + (t45 * t219 + t3 + (qJD(3) * t69 + t208) * pkin(8)) * t156 + (-t45 * t217 + t22 * t155 + t177 * t14 + (-t34 + t209) * pkin(8)) * t153, t115 * t35 + t156 * t2 - t71 * t85 + t247 * t96 - t246 * t67 + (t152 * t218 - t89) * t12 + (t12 * t216 + t152 * t6 + t177 * t9) * t153, t10 * t89 - t34 * t85 - t35 * t84 - t9 * t90 - t247 * t69 + t248 * t67 + t180 * t218 + (-t1 * t152 + t155 * t2 + (-t10 * t155 - t152 * t9) * qJD(4)) * t153, -t1 * t156 + t115 * t34 + t71 * t84 - t248 * t96 + t246 * t69 + t181 * t12 + (-t10 * t177 + t12 * t217 - t155 * t6) * t153, t1 * t84 - t248 * t10 + t115 * t6 - t246 * t12 + t2 * t85 - t247 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102 * t160, t102 ^ 2 - t160 ^ 2, t160 * t142 + t171, -t102 * t142 - t172, t185, -t82 * t102 - t52 * t177 + t199, -t160 * t82 - t51 * t177 + t165, t69 * t197 - t240, (-t34 - t251) * t155 + (-t35 - t249) * t152, -t102 * t69 + t197 * t96 + t239, t102 * t67 - t152 * t259 + t238, -t96 * t102, -pkin(3) * t35 - t13 * t102 - t52 * t67 + (-t22 + (-t63 - t242) * t96) * t155 + (t51 * t96 + t163) * t152, pkin(3) * t34 + t245 * t96 + t14 * t102 - t52 * t69 + (t22 + t212) * t152 + t163 * t155, t102 * t9 + t134 * t35 + t262 * t152 + t189 * t155 + t16 * t96 + t237 * t67, t15 * t67 - t16 * t69 + (t1 + t96 * t9 + (qJD(4) * t69 - t35) * pkin(9)) * t155 + ((qJD(4) * t67 - t34) * pkin(9) + t266) * t152, -t10 * t102 + t134 * t34 - t15 * t96 + t189 * t152 - t262 * t155 - t237 * t69, -t10 * t15 + t134 * t6 - t16 * t9 + t237 * t12 + (qJD(4) * t180 + t1 * t155 + t152 * t2) * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t250, -t67 ^ 2 + t260, t11, -t35 + t249, t71, -t45 * t69 + t168, t13 * t96 + t45 * t67 - t3, -t37 * t67 + t168 - t252 + 0.2e1 * t258, pkin(4) * t34 - qJ(5) * t35 + (t10 - t14) * t69 + (t9 - t231) * t67, 0.2e1 * t241 - t12 * t67 + t37 * t69 + (0.2e1 * qJD(5) - t13) * t96 + t3, -pkin(4) * t2 + qJ(5) * t1 + t10 * t231 - t12 * t37 - t14 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t250 - t71, t11, -t259 - t260, t252 + t266;];
tauc_reg = t4;
