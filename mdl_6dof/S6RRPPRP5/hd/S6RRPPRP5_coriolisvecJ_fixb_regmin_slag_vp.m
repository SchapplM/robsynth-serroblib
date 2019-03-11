% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPPRP5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:43:57
% EndTime: 2019-03-09 08:44:06
% DurationCPUTime: 3.20s
% Computational Cost: add. (4167->361), mult. (9512->479), div. (0->0), fcn. (6040->6), ass. (0->193)
t260 = cos(qJ(5));
t213 = qJD(5) * t260;
t167 = sin(qJ(2));
t219 = t167 * t260;
t279 = qJD(1) * t219 + t213;
t166 = sin(qJ(5));
t226 = t167 * qJD(1);
t229 = qJD(5) * t166;
t278 = -t166 * t226 - t229;
t147 = qJD(5) + t226;
t163 = sin(pkin(9));
t164 = cos(pkin(9));
t246 = -t279 * t163 + t278 * t164;
t277 = t246 * t147;
t168 = cos(qJ(2));
t232 = qJD(1) * t168;
t214 = t164 * t232;
t228 = t163 * qJD(2);
t111 = -t214 - t228;
t215 = t163 * t232;
t227 = t164 * qJD(2);
t112 = -t215 + t227;
t63 = -t260 * t111 + t166 * t112;
t276 = t63 ^ 2;
t181 = -t166 * t111 - t260 * t112;
t262 = t181 ^ 2;
t261 = pkin(3) + pkin(7);
t275 = t63 * t147;
t245 = t278 * t163 + t279 * t164;
t274 = t147 * t181;
t178 = -t166 * t163 + t260 * t164;
t92 = t178 * t168;
t179 = -t260 * t163 - t166 * t164;
t209 = t245 * t147;
t273 = (-t179 * t232 - t181) * qJD(2) + t209;
t150 = pkin(7) * t226;
t272 = qJD(3) + t150;
t223 = qJD(1) * qJD(2);
t271 = -0.2e1 * t223;
t165 = -pkin(2) - qJ(4);
t257 = -pkin(8) + t165;
t123 = t257 * t163;
t124 = t257 * t164;
t180 = -t166 * t123 + t260 * t124;
t241 = t163 * t167;
t190 = pkin(4) * t168 - pkin(8) * t241;
t151 = pkin(7) * t232;
t152 = pkin(3) * t232;
t121 = t151 + t152;
t154 = pkin(2) * t226;
t193 = -qJ(3) * t168 + qJ(4) * t167;
t97 = qJD(1) * t193 + t154;
t58 = t164 * t121 - t163 * t97;
t41 = qJD(1) * t190 + t58;
t240 = t164 * t167;
t222 = pkin(8) * t240;
t205 = qJD(1) * t222;
t59 = t163 * t121 + t164 * t97;
t50 = t205 + t59;
t270 = -qJD(4) * t179 - qJD(5) * t180 + t166 * t41 + t260 * t50;
t73 = t260 * t123 + t166 * t124;
t269 = -qJD(4) * t178 - qJD(5) * t73 + t166 * t50 - t260 * t41;
t220 = -pkin(4) * t164 - pkin(3);
t234 = -t220 * t226 + t272;
t160 = qJD(2) * qJ(3);
t267 = qJD(4) + t160;
t201 = qJD(2) * t219;
t191 = qJD(1) * t201;
t212 = t167 * t223;
t200 = t166 * t212;
t34 = -t111 * t213 + t112 * t229 - t163 * t191 - t164 * t200;
t266 = -t178 * t34 - t181 * t246;
t101 = t121 + t267;
t230 = qJD(2) * t168;
t265 = t167 * (-t101 + t267) - t165 * t230;
t35 = -qJD(5) * t181 + t163 * t200 - t164 * t191;
t210 = -t167 * qJ(3) - pkin(1);
t110 = t165 * t168 + t210;
t136 = t261 * t167;
t117 = t164 * t136;
t52 = t167 * pkin(4) + t117 + (pkin(8) * t168 - t110) * t163;
t239 = t164 * t168;
t69 = t164 * t110 + t163 * t136;
t56 = -pkin(8) * t239 + t69;
t184 = t166 * t52 + t260 * t56;
t175 = t190 * qJD(2);
t122 = t261 * t230;
t231 = qJD(2) * t167;
t153 = pkin(2) * t231;
t225 = t167 * qJD(3);
t173 = t193 * qJD(2) - t168 * qJD(4) - t225;
t75 = t153 + t173;
t48 = t164 * t122 - t163 * t75;
t31 = t175 + t48;
t49 = t163 * t122 + t164 * t75;
t38 = qJD(2) * t222 + t49;
t264 = -qJD(5) * t184 - t166 * t38 + t260 * t31;
t146 = pkin(2) * t212;
t61 = qJD(1) * t173 + t146;
t211 = t168 * t223;
t145 = pkin(7) * t211;
t96 = t145 + (-qJD(4) + t152) * qJD(2);
t32 = -t163 * t61 + t164 * t96;
t19 = qJD(1) * t175 + t32;
t33 = t163 * t96 + t164 * t61;
t25 = qJD(2) * t205 + t33;
t85 = t110 * qJD(1);
t224 = pkin(3) * t226 + t272;
t90 = t165 * qJD(2) + t224;
t45 = -t163 * t85 + t164 * t90;
t26 = pkin(4) * t226 - t112 * pkin(8) + t45;
t46 = t163 * t90 + t164 * t85;
t30 = t111 * pkin(8) + t46;
t183 = -t166 * t19 - t26 * t213 + t229 * t30 - t260 * t25;
t199 = qJ(6) * t211;
t1 = t147 * qJD(6) - t183 + t199;
t204 = pkin(5) * t211;
t208 = t166 * t25 - t260 * t19 + t30 * t213 + t26 * t229;
t2 = -t204 + t208;
t8 = -t166 * t30 + t260 * t26;
t247 = qJD(6) - t8;
t6 = -t147 * pkin(5) + t247;
t9 = t166 * t26 + t260 * t30;
t7 = t147 * qJ(6) + t9;
t263 = -t1 * t179 - t178 * t2 + t245 * t7 - t246 * t6;
t70 = -t111 * pkin(4) + t101;
t13 = t63 * pkin(5) + qJ(6) * t181 + t70;
t259 = t13 * t181;
t258 = t181 * t63;
t256 = qJ(6) * t232 + t270;
t255 = -pkin(5) * t232 + t269;
t254 = -pkin(5) * t245 + qJ(6) * t246 + t178 * qJD(6) - t234;
t252 = qJD(2) * pkin(2);
t250 = t168 * t45;
t249 = t168 * t46;
t120 = t261 * t231;
t159 = qJD(2) * qJD(3);
t95 = -qJD(1) * t120 + t159;
t248 = t95 * t163;
t244 = qJD(2) * t180;
t243 = qJD(2) * t73;
t161 = t167 ^ 2;
t170 = qJD(1) ^ 2;
t242 = t161 * t170;
t237 = t168 * t170;
t169 = qJD(2) ^ 2;
t236 = t169 * t167;
t235 = t169 * t168;
t148 = t163 * pkin(4) + qJ(3);
t137 = t261 * t168;
t233 = -t168 ^ 2 + t161;
t130 = -t168 * pkin(2) + t210;
t109 = qJD(1) * t130;
t221 = t167 * t237;
t102 = pkin(4) * t239 + t137;
t207 = pkin(1) * t271;
t206 = qJD(3) - t252;
t203 = t178 * t211 + t277;
t196 = t33 * t163 + t32 * t164;
t195 = -t163 * t45 + t164 * t46;
t192 = -0.2e1 * qJD(2) * t109;
t177 = -qJ(3) * t230 - t225;
t86 = qJD(1) * t177 + t146;
t99 = t153 + t177;
t189 = pkin(7) * t169 + qJD(1) * t99 + t86;
t188 = t9 * t147 - t208;
t186 = -t166 * t56 + t260 * t52;
t182 = t166 * t31 + t52 * t213 - t229 * t56 + t260 * t38;
t88 = (-pkin(7) + t220) * t231;
t93 = t179 * t168;
t174 = t34 + t275;
t74 = qJD(1) * t88 + t159;
t128 = pkin(7) * t212 - t159;
t129 = t150 + t206;
t135 = -t151 - t160;
t172 = -t128 * t168 + (t129 * t168 + (t135 + t151) * t167) * qJD(2);
t5 = t35 * pkin(5) + t34 * qJ(6) + qJD(6) * t181 + t74;
t171 = t35 - t274;
t118 = -qJ(3) * t232 + t154;
t89 = t109 * t226;
t68 = -t163 * t110 + t117;
t60 = -pkin(5) * t179 - qJ(6) * t178 + t148;
t55 = -t166 * t167 * t228 - qJD(5) * t93 + t164 * t201;
t54 = -qJD(5) * t92 - t179 * t231;
t37 = t92 * pkin(5) - t93 * qJ(6) + t102;
t24 = -pkin(5) * t181 + t63 * qJ(6);
t16 = -t167 * pkin(5) - t186;
t15 = t167 * qJ(6) + t184;
t14 = -t34 + t275;
t10 = -t55 * pkin(5) - t54 * qJ(6) - t93 * qJD(6) + t88;
t4 = -pkin(5) * t230 - t264;
t3 = qJ(6) * t230 + t167 * qJD(6) + t182;
t11 = [0, 0, 0, 0.2e1 * t167 * t211, t233 * t271, t235, -t236, 0, -pkin(7) * t235 + t167 * t207, pkin(7) * t236 + t168 * t207, t172, t167 * t192 + t168 * t189, -t189 * t167 + t168 * t192, pkin(7) * t172 + t109 * t99 + t86 * t130, t95 * t239 + t120 * t111 + (qJD(1) * t48 + t32) * t167 + (-t101 * t240 + t250 + (-t137 * t240 + t168 * t68) * qJD(1)) * qJD(2), -t168 * t248 - t120 * t112 + (-qJD(1) * t49 - t33) * t167 + (t101 * t241 - t249 + (t137 * t241 - t168 * t69) * qJD(1)) * qJD(2), t49 * t111 - t48 * t112 + (t163 * t32 - t164 * t33) * t168 + ((-t163 * t68 + t164 * t69) * qJD(1) + t195) * t231, -t101 * t120 + t95 * t137 + t32 * t68 + t33 * t69 + t45 * t48 + t46 * t49, -t181 * t54 - t34 * t93, -t181 * t55 + t34 * t92 - t93 * t35 - t54 * t63, t54 * t147 - t34 * t167 + (qJD(1) * t93 - t181) * t230, t55 * t147 - t35 * t167 + (-qJD(1) * t92 - t63) * t230 (t147 + t226) * t230, t102 * t35 + t147 * t264 - t208 * t167 + t186 * t211 + t8 * t230 - t70 * t55 + t88 * t63 + t74 * t92, -t182 * t147 + t183 * t167 - t88 * t181 - t102 * t34 + t74 * t93 + t70 * t54 + (-t184 * qJD(1) - t9) * t230, t10 * t63 - t13 * t55 - t4 * t147 - t2 * t167 + t37 * t35 + t5 * t92 + (-qJD(1) * t16 - t6) * t230, -t1 * t92 - t15 * t35 - t16 * t34 - t181 * t4 + t2 * t93 - t3 * t63 + t6 * t54 + t7 * t55, t1 * t167 + t10 * t181 - t13 * t54 + t3 * t147 + t37 * t34 - t5 * t93 + (qJD(1) * t15 + t7) * t230, t1 * t15 + t13 * t10 + t2 * t16 + t7 * t3 + t5 * t37 + t6 * t4; 0, 0, 0, -t221, t233 * t170, 0, 0, 0, t170 * pkin(1) * t167, pkin(1) * t237 ((-t135 - t160) * t167 + (-t129 + t206) * t168) * qJD(1), -t118 * t232 + t89, 0.2e1 * t159 + (t109 * t168 + t118 * t167) * qJD(1), -t128 * qJ(3) - t135 * qJD(3) - t109 * t118 + (-t135 * t167 + (-t129 - t252) * t168) * qJD(1) * pkin(7), t248 - t224 * t111 + (-t164 * t265 - t167 * t58 - t250) * qJD(1), t95 * t164 + t224 * t112 + (t163 * t265 + t167 * t59 + t249) * qJD(1), -t59 * t111 + t58 * t112 + (qJD(4) * t112 - t226 * t46 - t32) * t164 + (-qJD(4) * t111 + t226 * t45 - t33) * t163, t95 * qJ(3) - t45 * t58 - t46 * t59 + t196 * t165 + t224 * t101 + (-t163 * t46 - t164 * t45) * qJD(4), t266, -t178 * t35 - t179 * t34 + t181 * t245 - t246 * t63, t181 * t232 + t203, -t209 + (qJD(2) * t179 + t63) * t232, -t147 * t232, -t74 * t179 + t148 * t35 + t245 * t70 + t234 * t63 + t269 * t147 + (-t8 + t244) * t232, t74 * t178 - t148 * t34 + t246 * t70 - t234 * t181 + t270 * t147 + (t9 - t243) * t232, -t5 * t179 + t60 * t35 - t254 * t63 + t255 * t147 + t245 * t13 + (t6 + t244) * t232, t180 * t34 + t181 * t255 + t256 * t63 - t73 * t35 - t263, -t5 * t178 + t60 * t34 - t254 * t181 - t256 * t147 - t246 * t13 + (-t7 + t243) * t232, t1 * t73 - t254 * t13 - t180 * t2 - t255 * t6 - t256 * t7 + t5 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t221, -t169 - t242, t135 * qJD(2) + t145 + t89, -t163 * t242 + (t111 + t214) * qJD(2), -t164 * t242 + (-t112 - t215) * qJD(2) (t111 * t164 + t112 * t163) * t226, -t101 * qJD(2) + t195 * t226 + t196, 0, 0, 0, 0, 0, -qJD(2) * t63 + t203, -t273, t277 + (t178 * t232 - t63) * qJD(2), t179 * t35 - t245 * t63 - t266, t273, -t13 * qJD(2) + t263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t112 - t227) * t226 (t111 + t228) * t226, -t111 ^ 2 - t112 ^ 2, -t46 * t111 + t45 * t112 + t95, 0, 0, 0, 0, 0, t171, -t174, t171, -t262 - t276, t174, t181 * t6 + t63 * t7 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t258, t262 - t276, t14, -t35 - t274, t211, t181 * t70 + t188, t8 * t147 + t70 * t63 + t183, -t24 * t63 + t188 + 0.2e1 * t204 + t259, pkin(5) * t34 - t35 * qJ(6) - (t7 - t9) * t181 + (t6 - t247) * t63, 0.2e1 * t199 - t13 * t63 - t24 * t181 + (0.2e1 * qJD(6) - t8) * t147 - t183, -t2 * pkin(5) + t1 * qJ(6) - t13 * t24 + t247 * t7 - t6 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t211 - t258, t14, -t147 ^ 2 - t262, -t7 * t147 + t2 - t259;];
tauc_reg  = t11;
