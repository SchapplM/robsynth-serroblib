% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRPPR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:16:26
% EndTime: 2019-03-08 21:16:35
% DurationCPUTime: 3.08s
% Computational Cost: add. (2260->377), mult. (5982->526), div. (0->0), fcn. (4388->10), ass. (0->202)
t150 = cos(qJ(3));
t225 = qJD(2) * t150;
t133 = qJD(6) + t225;
t270 = qJD(6) - t133;
t144 = cos(pkin(11));
t142 = sin(pkin(11));
t224 = qJD(3) * t142;
t147 = sin(qJ(3));
t226 = qJD(2) * t147;
t103 = t144 * t226 + t224;
t235 = t144 * t150;
t148 = sin(qJ(2));
t143 = sin(pkin(6));
t227 = qJD(2) * t143;
t201 = t148 * t227;
t151 = cos(qJ(2));
t200 = t151 * t227;
t185 = t150 * t200;
t145 = cos(pkin(6));
t238 = t143 * t148;
t92 = -t145 * t150 + t147 * t238;
t59 = -qJD(3) * t92 + t185;
t33 = t142 * t201 + t144 * t59;
t93 = t145 * t147 + t150 * t238;
t58 = -t142 * t143 * t151 + t144 * t93;
t186 = t147 * t200;
t60 = qJD(3) * t93 + t186;
t269 = (qJD(3) * (-t147 * t58 + t92 * t235) + t150 * t33) * qJD(2) + t103 * t60;
t214 = t144 * qJD(3);
t101 = t142 * t226 - t214;
t99 = t103 ^ 2;
t268 = -t101 ^ 2 - t99;
t229 = qJD(1) * t143;
t203 = t151 * t229;
t204 = t148 * t229;
t239 = t142 * t150;
t177 = pkin(3) * t147 - qJ(4) * t150;
t87 = t177 * qJD(3) - qJD(4) * t147;
t267 = t203 * t239 + (-t204 + t87) * t144;
t217 = qJD(5) * t150;
t223 = qJD(3) * t147;
t234 = t144 * t151;
t72 = (t142 * t148 + t150 * t234) * t229;
t266 = qJ(5) * t223 - t217 - t72;
t213 = qJD(2) * qJD(3);
t198 = t150 * t213;
t124 = t142 * t198;
t146 = sin(qJ(6));
t149 = cos(qJ(6));
t183 = t144 * t198;
t48 = t101 * t146 + t103 * t149;
t19 = t48 * qJD(6) - t149 * t124 + t146 * t183;
t264 = pkin(4) + pkin(5);
t263 = pkin(9) * t147;
t262 = -pkin(9) + qJ(4);
t207 = -pkin(8) * t142 - pkin(4);
t212 = pkin(9) * t235;
t261 = (-t212 + (-pkin(5) + t207) * t147) * qJD(3) - t267;
t236 = t144 * t147;
t80 = t142 * t87;
t260 = -t80 - (-pkin(8) * t236 + pkin(9) * t239) * qJD(3) - t266;
t184 = qJD(1) * t200;
t112 = qJD(2) * pkin(8) + t204;
t228 = qJD(1) * t145;
t243 = -t147 * t112 + t150 * t228;
t39 = t150 * t184 + (qJD(4) + t243) * qJD(3);
t52 = (t87 + t204) * qJD(2);
t11 = t142 * t52 + t144 * t39;
t210 = pkin(8) * t223;
t191 = t144 * t210;
t56 = t80 - t191;
t259 = t56 + t266;
t258 = t207 * t223 - t267;
t127 = t147 * t228;
t74 = t150 * t112 + t127;
t68 = qJD(3) * qJ(4) + t74;
t242 = qJ(4) * t147;
t178 = -pkin(3) * t150 - t242;
t117 = -pkin(2) + t178;
t75 = qJD(2) * t117 - t203;
t24 = t142 * t75 + t144 * t68;
t192 = t142 * t210;
t257 = t192 + t267;
t256 = t56 - t72;
t107 = t177 * qJD(2);
t35 = t142 * t107 + t144 * t243;
t105 = t142 * t146 + t144 * t149;
t163 = t105 * t150;
t255 = -qJD(2) * t163 - t105 * qJD(6);
t202 = t142 * t225;
t233 = t146 * t150;
t208 = t144 * t233;
t215 = qJD(6) * t149;
t216 = qJD(6) * t146;
t254 = -qJD(2) * t208 + t142 * t215 - t144 * t216 + t149 * t202;
t253 = qJD(2) * pkin(2);
t222 = qJD(3) * t150;
t43 = qJD(3) * t127 + t112 * t222 + t147 * t184;
t190 = -pkin(4) * t124 - t43;
t157 = -qJ(5) * t183 - t190;
t219 = qJD(5) * t103;
t13 = t157 - t219;
t251 = t13 * t142;
t250 = t13 * t144;
t46 = -t149 * t101 + t103 * t146;
t249 = t133 * t46;
t248 = t133 * t48;
t247 = t142 * t43;
t246 = t144 * t43;
t241 = qJ(5) * t144;
t240 = t142 * t149;
t153 = qJD(2) ^ 2;
t237 = t143 * t153;
t152 = qJD(3) ^ 2;
t232 = t152 * t147;
t231 = t152 * t150;
t77 = pkin(8) * t235 + t142 * t117;
t141 = t150 ^ 2;
t230 = t147 ^ 2 - t141;
t221 = qJD(4) * t103;
t220 = qJD(4) * t144;
t218 = qJD(5) * t142;
t197 = t147 * t213;
t211 = qJ(5) * t197 + t11;
t26 = qJ(5) * t226 + t35;
t209 = t148 * t237;
t10 = -t142 * t39 + t144 * t52;
t159 = (-t264 * t147 - t212) * qJD(2);
t4 = qJD(3) * t159 - t10;
t5 = (pkin(9) * t224 - qJD(5)) * t225 + t211;
t206 = -t146 * t5 + t149 * t4;
t205 = qJ(4) * t223;
t199 = qJD(5) * t236;
t196 = qJ(5) * t142 + pkin(3);
t23 = -t142 * t68 + t144 * t75;
t34 = t107 * t144 - t142 * t243;
t195 = qJD(3) * pkin(3) - qJD(4);
t168 = -t264 * t142 + t241;
t194 = -t168 * t225 + t218 + t74;
t176 = pkin(4) * t142 - t241;
t41 = t176 * t225 + t74;
t193 = t41 + t218;
t131 = pkin(8) * t239;
t76 = t117 * t144 - t131;
t189 = t101 * t203;
t188 = t103 * t203;
t187 = t147 * t203;
t182 = t146 * t4 + t149 * t5;
t67 = -qJ(5) * t150 + t77;
t113 = -t203 - t253;
t181 = -t113 - t203;
t180 = -t103 * t225 + t124;
t16 = -qJ(5) * t225 + t24;
t12 = pkin(9) * t101 + t16;
t15 = pkin(4) * t225 + qJD(5) - t23;
t8 = pkin(5) * t225 - pkin(9) * t103 + t15;
t2 = t12 * t149 + t146 * t8;
t179 = t12 * t146 - t149 * t8;
t139 = t150 * pkin(4);
t44 = pkin(5) * t150 + t131 + t139 + (-t117 - t263) * t144;
t50 = t142 * t263 + t67;
t175 = t146 * t50 - t149 * t44;
t174 = t146 * t44 + t149 * t50;
t57 = t142 * t93 + t143 * t234;
t173 = -t146 * t58 + t149 * t57;
t172 = t146 * t57 + t149 * t58;
t106 = -t144 * t146 + t240;
t171 = t195 + t243;
t170 = pkin(8) + t176;
t120 = t262 * t142;
t166 = pkin(9) * t202 - qJD(6) * t120 - t220 + t26;
t121 = t262 * t144;
t165 = qJD(4) * t142 - qJD(6) * t121 - t159 + t34;
t164 = -t101 * t215 + t103 * t216 - t146 * t124 - t149 * t183;
t83 = t105 * t147;
t162 = -pkin(8) + t168;
t160 = qJ(5) * t103 + t171;
t158 = qJD(3) * (-t181 - t253);
t32 = t142 * t59 - t144 * t201;
t156 = t60 * t101 + t92 * t124 + (t150 * t32 - t57 * t223) * qJD(2);
t155 = t103 * t32 - t33 * t101 + (-t142 * t58 + t144 * t57) * t198;
t123 = qJD(4) * t202;
t114 = -pkin(4) * t144 - t196;
t95 = t264 * t144 + t196;
t85 = t101 * t225;
t84 = t101 * t220;
t82 = t146 * t236 - t147 * t240;
t81 = t170 * t147;
t69 = t139 - t76;
t63 = t162 * t147;
t61 = t85 + t183;
t51 = t170 * t222 - t199;
t42 = t162 * t222 + t199;
t38 = qJD(3) * t208 + qJD(6) * t83 - t222 * t240;
t37 = t106 * t147 * qJD(6) + qJD(3) * t163;
t29 = -pkin(4) * t226 - t34;
t21 = pkin(4) * t101 - t160;
t14 = -t264 * t101 + t160;
t9 = t219 + (-pkin(5) * t142 + t241) * t198 + t190;
t7 = -pkin(4) * t197 - t10;
t6 = -qJD(2) * t217 + t211;
t1 = [0, 0, -t209, -t151 * t237, 0, 0, 0, 0, 0, -t150 * t209 + (-t60 - t186) * qJD(3), t147 * t209 + (-t59 - t185) * qJD(3), t156, t269, t155, -t10 * t57 + t11 * t58 - t171 * t60 - t23 * t32 + t24 * t33 + t43 * t92, t156, t155, -t269, t13 * t92 + t15 * t32 + t16 * t33 + t21 * t60 + t57 * t7 + t58 * t6, 0, 0, 0, 0, 0 (-qJD(6) * t172 - t146 * t33 + t149 * t32) * t133 - t173 * t197 - t60 * t46 - t92 * t19 -(qJD(6) * t173 + t146 * t32 + t149 * t33) * t133 + t172 * t197 - t60 * t48 + t92 * t164; 0, 0, 0, 0, 0.2e1 * t150 * t197, -0.2e1 * t230 * t213, t231, -t232, 0, -pkin(8) * t231 + t147 * t158, pkin(8) * t232 + t150 * t158 (-t189 + t247 + (qJD(2) * t76 + t23) * qJD(3)) * t147 + (-t10 + (pkin(8) * t101 - t142 * t171) * qJD(3) + (t192 - t257) * qJD(2)) * t150 (-t188 + t246 + (-qJD(2) * t77 - t24) * qJD(3)) * t147 + (t11 + (pkin(8) * t103 - t144 * t171) * qJD(3) + (t191 + t256) * qJD(2)) * t150 (-t10 * t144 - t11 * t142) * t147 - t257 * t103 - t256 * t101 + (-t142 * t24 - t144 * t23 + (-t142 * t77 - t144 * t76) * qJD(2)) * t222, t171 * t187 + t10 * t76 + t11 * t77 + t256 * t24 + t257 * t23 + (t147 * t43 - t171 * t222) * pkin(8), t101 * t51 + (-t189 + t251 + (-qJD(2) * t69 - t15) * qJD(3)) * t147 + (t21 * t224 + t7 + (t224 * t81 + t258) * qJD(2)) * t150 (-t142 * t6 + t144 * t7) * t147 + t258 * t103 - t259 * t101 + (-t142 * t16 + t144 * t15 + (-t142 * t67 + t144 * t69) * qJD(2)) * t222, -t103 * t51 + (t188 - t250 + (qJD(2) * t67 + t16) * qJD(3)) * t147 + (-t21 * t214 - t6 + (-t214 * t81 - t259) * qJD(2)) * t150, t13 * t81 + t6 * t67 + t69 * t7 + (t51 - t187) * t21 + t259 * t16 + t258 * t15, -t164 * t83 + t37 * t48, t164 * t82 - t19 * t83 - t37 * t46 - t38 * t48, t133 * t37 - t150 * t164 + (-qJD(2) * t83 - t48) * t223, -t133 * t38 - t150 * t19 + (qJD(2) * t82 + t46) * t223 (-t133 - t225) * t223, t206 * t150 + t42 * t46 + t63 * t19 + t9 * t82 + t14 * t38 + (t260 * t146 + t261 * t149) * t133 + (-t133 * t174 - t150 * t2) * qJD(6) + (t46 * t203 + (qJD(2) * t175 + t179) * qJD(3)) * t147, -t182 * t150 + t42 * t48 - t63 * t164 + t9 * t83 + t14 * t37 + (-t261 * t146 + t260 * t149) * t133 + (t133 * t175 + t150 * t179) * qJD(6) + (t48 * t203 + (qJD(2) * t174 + t2) * qJD(3)) * t147; 0, 0, 0, 0, -t147 * t153 * t150, t230 * t153, 0, 0, 0, qJD(3) * t74 - t113 * t226 - t43, t181 * t225, -t101 * t74 - t246 + t123 + (-t147 * t23 + t150 * t34 + (t178 * qJD(3) + t150 * t171) * t142) * qJD(2), -t103 * t74 + t247 + (t147 * t24 - t150 * t35 + (-t205 + (-t195 + t171) * t150) * t144) * qJD(2), t101 * t35 + t103 * t34 - t84 + (t225 * t23 + t11) * t144 + (t225 * t24 - t10 + t221) * t142, -pkin(3) * t43 - t23 * t34 - t24 * t35 + t171 * t74 + (-t142 * t23 + t144 * t24) * qJD(4) + (-t10 * t142 + t11 * t144) * qJ(4), -t250 + t123 - t193 * t101 + (t147 * t15 - t150 * t29 + (-t150 * t21 + (t114 * t150 - t242) * qJD(3)) * t142) * qJD(2), t101 * t26 - t103 * t29 - t84 + (-t15 * t225 + t6) * t144 + (t16 * t225 + t221 + t7) * t142, -t251 + t193 * t103 + (-t147 * t16 + t150 * t26 + (t205 + (-qJD(3) * t114 - qJD(4) + t21) * t150) * t144) * qJD(2), t114 * t13 - t15 * t29 - t16 * t26 - t21 * t41 + (qJ(4) * t6 + qJD(4) * t16) * t144 + (qJ(4) * t7 + qJD(4) * t15 - qJD(5) * t21) * t142, -t106 * t164 + t255 * t48, t105 * t164 - t106 * t19 - t254 * t48 - t255 * t46, t255 * t133 + (-qJD(3) * t106 + t48) * t226, -t254 * t133 + (qJD(3) * t105 - t46) * t226, t133 * t226, t9 * t105 + t95 * t19 + t194 * t46 + t254 * t14 + (t146 * t166 + t149 * t165) * t133 + (-(t120 * t149 - t121 * t146) * qJD(3) - t179) * t226, t9 * t106 - t95 * t164 + t194 * t48 + t255 * t14 + (-t146 * t165 + t149 * t166) * t133 + ((t120 * t146 + t121 * t149) * qJD(3) - t2) * t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t180, t61, t268, t101 * t24 + t103 * t23 + t43, t180, t268, -t61, t101 * t16 + (-qJD(5) - t15) * t103 + t157, 0, 0, 0, 0, 0, -t19 - t248, t164 + t249; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101 * t103 - t197, -t85 + t183, -t141 * t153 - t99, t103 * t21 + (-pkin(4) * t223 + t150 * t16) * qJD(2) - t10, 0, 0, 0, 0, 0, -t133 * t216 - t103 * t46 + (-t133 * t233 - t149 * t223) * qJD(2), -t133 * t215 - t103 * t48 + (-t133 * t149 * t150 + t146 * t223) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48 * t46, -t46 ^ 2 + t48 ^ 2, -t164 + t249, -t19 + t248, -t197, -t14 * t48 - t2 * t270 + t206, t14 * t46 + t179 * t270 - t182;];
tauc_reg  = t1;
