% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% tauc_reg [6x31]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRRP10_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP10_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:32:32
% EndTime: 2019-03-09 06:32:40
% DurationCPUTime: 3.57s
% Computational Cost: add. (5234->407), mult. (11129->545), div. (0->0), fcn. (7152->6), ass. (0->201)
t152 = sin(qJ(4));
t153 = sin(qJ(3));
t224 = qJD(1) * t153;
t199 = t152 * t224;
t267 = -pkin(9) - pkin(8);
t200 = qJD(4) * t267;
t156 = cos(qJ(3));
t181 = pkin(3) * t156 + pkin(8) * t153;
t114 = t181 * qJD(1);
t157 = -pkin(1) - pkin(7);
t134 = t157 * qJD(1) + qJD(2);
t155 = cos(qJ(4));
t233 = t155 * t156;
t251 = t152 * t114 + t134 * t233;
t282 = pkin(9) * t199 - t152 * t200 + t251;
t100 = t155 * t114;
t240 = t152 * t156;
t202 = t134 * t240;
t237 = t153 * t155;
t207 = pkin(9) * t237;
t281 = t155 * t200 + t202 - t100 - (pkin(4) * t156 + t207) * qJD(1);
t220 = qJD(3) * t155;
t223 = qJD(1) * t156;
t109 = -t152 * t223 + t220;
t198 = t155 * t223;
t222 = qJD(3) * t152;
t110 = t198 + t222;
t151 = sin(qJ(5));
t154 = cos(qJ(5));
t175 = t109 * t151 + t154 * t110;
t62 = -t154 * t109 + t110 * t151;
t256 = qJD(3) * pkin(3);
t102 = -t134 * t156 - t256;
t71 = -pkin(4) * t109 + t102;
t20 = pkin(5) * t62 - qJ(6) * t175 + t71;
t280 = t20 * t62;
t279 = t62 * t71;
t210 = qJD(4) + qJD(5);
t213 = qJD(5) * t154;
t216 = qJD(4) * t155;
t235 = t154 * t155;
t241 = t151 * t152;
t258 = t151 * t199 - t154 * t216 - t155 * t213 + t210 * t241 - t224 * t235;
t113 = t151 * t155 + t152 * t154;
t69 = t210 * t113;
t96 = t113 * qJD(1);
t257 = t153 * t96 + t69;
t265 = t175 * t62;
t215 = qJD(4) * t156;
t194 = t155 * t215;
t221 = qJD(3) * t153;
t197 = t152 * t221;
t278 = t194 - t197;
t268 = t175 ^ 2;
t277 = -t62 ^ 2 + t268;
t141 = qJD(4) + t224;
t132 = qJD(5) + t141;
t195 = t152 * t215;
t196 = t153 * t220;
t165 = -t195 - t196;
t211 = qJD(3) * qJD(4);
t162 = qJD(1) * t165 + t155 * t211;
t192 = qJD(1) * t215;
t228 = t152 * t211 + t155 * t192;
t166 = qJD(1) * t197 - t228;
t214 = qJD(5) * t151;
t23 = -t109 * t213 + t110 * t214 - t151 * t166 - t154 * t162;
t7 = t132 * t62 - t23;
t32 = pkin(5) * t175 + qJ(6) * t62;
t208 = 0.2e1 * qJD(1);
t275 = t175 * t71;
t266 = t20 * t175;
t127 = t267 * t152;
t128 = t267 * t155;
t76 = t127 * t151 - t128 * t154;
t274 = -qJD(5) * t76 + t282 * t151 + t281 * t154;
t174 = t127 * t154 + t128 * t151;
t273 = -qJD(5) * t174 - t281 * t151 + t282 * t154;
t117 = t153 * t134;
t218 = qJD(4) * t152;
t183 = -t117 + (t199 + t218) * pkin(4);
t150 = t156 ^ 2;
t226 = t153 ^ 2 - t150;
t24 = qJD(5) * t175 + t151 * t162 - t154 * t166;
t272 = t132 * t175 - t24;
t180 = pkin(3) * t153 - pkin(8) * t156;
t118 = qJ(2) + t180;
t106 = t155 * t118;
t239 = t152 * t157;
t191 = pkin(4) - t239;
t60 = -pkin(9) * t233 + t191 * t153 + t106;
t236 = t153 * t157;
t131 = t155 * t236;
t227 = t152 * t118 + t131;
t70 = -pkin(9) * t240 + t227;
t176 = t151 * t60 + t154 * t70;
t107 = t181 * qJD(3) + qJD(2);
t92 = t155 * t107;
t30 = t92 + (-t131 + (pkin(9) * t156 - t118) * t152) * qJD(4) + (t191 * t156 + t207) * qJD(3);
t219 = qJD(3) * t156;
t193 = t155 * t219;
t204 = t152 * t107 + t118 * t216 + t157 * t193;
t217 = qJD(4) * t153;
t33 = -t278 * pkin(9) - t217 * t239 + t204;
t271 = -qJD(5) * t176 - t151 * t33 + t154 * t30;
t270 = t210 * t153;
t112 = -t235 + t241;
t89 = t112 * t156;
t261 = -qJD(3) * t89 - t113 * t270 - t96;
t88 = t112 * t153;
t269 = (-t153 * t175 - t88 * t223) * qJD(3) + t261 * t132 - t156 * t23;
t264 = qJ(6) * t223 + t273;
t263 = -pkin(5) * t223 + t274;
t262 = t257 * pkin(5) + t258 * qJ(6) - qJD(6) * t113 + t183;
t87 = t113 * t156;
t260 = qJD(3) * t87 + (-qJD(1) - t270) * t112;
t101 = qJD(3) * pkin(8) + t117;
t234 = t155 * t101;
t98 = t118 * qJD(1);
t254 = t152 * t98;
t56 = t234 + t254;
t48 = pkin(9) * t109 + t56;
t255 = t151 * t48;
t253 = t154 * t48;
t55 = -t101 * t152 + t155 * t98;
t47 = -pkin(9) * t110 + t55;
t15 = t154 * t47 - t255;
t250 = pkin(4) * t213 + qJD(6) - t15;
t249 = qJD(3) * t174;
t248 = qJD(3) * t76;
t247 = t102 * t152;
t246 = t102 * t153;
t245 = t109 * t156;
t244 = t110 * t141;
t243 = t141 * t152;
t242 = t141 * t155;
t238 = t153 * t141;
t158 = qJD(3) ^ 2;
t232 = t158 * t153;
t231 = t158 * t156;
t159 = qJD(1) ^ 2;
t230 = t159 * qJ(2);
t41 = pkin(4) * t141 + t47;
t10 = t154 * t41 - t255;
t229 = qJD(6) - t10;
t225 = -t158 - t159;
t212 = qJD(1) * qJD(3);
t85 = t107 * qJD(1);
t209 = t134 * t193 + t152 * t85 + t98 * t216;
t205 = qJD(2) * t208;
t203 = t152 * t238;
t201 = t141 * t240;
t148 = -pkin(4) * t155 - pkin(3);
t144 = t156 * t212;
t185 = -qJD(4) + t224;
t79 = t155 * t85;
t16 = t79 + (-t234 + (pkin(9) * t223 - t98) * t152) * qJD(4) + ((pkin(4) * qJD(1) - t134 * t152) * t156 + t185 * pkin(9) * t155) * qJD(3);
t167 = -t101 * t218 + t209;
t19 = pkin(9) * t166 + t167;
t190 = -t151 * t16 - t154 * t19 - t41 * t213 + t48 * t214;
t189 = t151 * t19 - t154 * t16 + t48 * t213 + t41 * t214;
t188 = t141 * t157 + t101;
t108 = pkin(4) * t240 - t156 * t157;
t187 = -t109 + t220;
t186 = -t110 + t222;
t184 = pkin(5) * t144;
t14 = t151 * t47 + t253;
t182 = pkin(4) * t214 - t14;
t11 = t151 * t41 + t253;
t177 = -t151 * t70 + t154 * t60;
t173 = t187 * t153;
t172 = t186 * t153;
t119 = t132 * qJD(6);
t136 = qJ(6) * t144;
t1 = t136 + t119 - t190;
t171 = t10 * t132 + t190;
t170 = t11 * t132 - t189;
t168 = t151 * t30 + t154 * t33 + t60 * t213 - t70 * t214;
t77 = t278 * pkin(4) + t157 * t221;
t2 = -t184 + t189;
t164 = t180 * qJD(3) + t246;
t52 = -pkin(4) * t166 + t134 * t221;
t86 = t113 * t153;
t163 = t62 * t221 + (-t86 * t212 - t24) * t156 - t260 * t132;
t147 = -pkin(4) * t154 - pkin(5);
t143 = pkin(4) * t151 + qJ(6);
t129 = t153 * t144;
t59 = pkin(5) * t112 - qJ(6) * t113 + t148;
t46 = pkin(5) * t87 + qJ(6) * t89 + t108;
t40 = -t214 * t240 + (t210 * t233 - t197) * t154 + t165 * t151;
t38 = -t151 * t197 + t154 * t196 + t156 * t69;
t29 = -pkin(5) * t153 - t177;
t28 = qJ(6) * t153 + t176;
t26 = pkin(4) * t110 + t32;
t9 = qJ(6) * t132 + t11;
t8 = -pkin(5) * t132 + t229;
t6 = pkin(5) * t40 + qJ(6) * t38 + qJD(6) * t89 + t77;
t5 = -pkin(5) * t219 - t271;
t4 = qJ(6) * t219 + qJD(6) * t153 + t168;
t3 = t24 * pkin(5) + t23 * qJ(6) - qJD(6) * t175 + t52;
t12 = [0, 0, 0, 0, t205, qJ(2) * t205, -0.2e1 * t129, 0.2e1 * t226 * t212, -t232, -t231, 0, -t157 * t232 + (qJ(2) * t219 + qJD(2) * t153) * t208, -t157 * t231 + (-qJ(2) * t221 + qJD(2) * t156) * t208, t110 * t165 + t162 * t233, -t228 * t233 + (-0.2e1 * t152 * t109 - t110 * t155) * t215 + (-t155 * t109 + (t110 + 0.2e1 * t198) * t152) * t221 (-t141 - t224) * t195 + (t110 * t156 + (t150 * qJD(1) + (-t141 - t185) * t153) * t155) * qJD(3), -t141 * t194 - t228 * t153 + (t245 + (t226 * qJD(1) + t238) * t152) * qJD(3), t141 * t219 + t129 (-t118 * t218 + t92) * t141 + (-t157 * t228 + t102 * t216 + (qJD(1) * t106 - t141 * t239 + t55) * qJD(3)) * t156 + (t79 + (-t109 * t157 - t247) * qJD(3) + (-t155 * t188 - t254) * qJD(4)) * t153, -t204 * t141 - t209 * t153 + ((t157 * t223 - t102) * t156 + t188 * t153) * t218 + (t110 * t236 + (-t227 * qJD(1) - t56) * t156 + (-t246 + (t157 * t185 + t117) * t156) * t155) * qJD(3), -t175 * t38 + t23 * t89, -t175 * t40 + t23 * t87 + t24 * t89 + t38 * t62, -t38 * t132 - t23 * t153 + (-qJD(1) * t89 + t175) * t219, -t40 * t132 - t24 * t153 + (-qJD(1) * t87 - t62) * t219, t132 * t219 + t129, t271 * t132 - t189 * t153 + t77 * t62 + t108 * t24 + t52 * t87 + t71 * t40 + (qJD(1) * t177 + t10) * t219, -t168 * t132 + t190 * t153 + t77 * t175 - t108 * t23 - t52 * t89 - t71 * t38 + (-t176 * qJD(1) - t11) * t219, -t5 * t132 - t2 * t153 + t20 * t40 + t46 * t24 + t3 * t87 + t6 * t62 + (-qJD(1) * t29 - t8) * t219, -t1 * t87 + t175 * t5 - t2 * t89 - t23 * t29 - t24 * t28 - t38 * t8 - t4 * t62 - t40 * t9, t1 * t153 + t4 * t132 + t20 * t38 + t46 * t23 + t3 * t89 - t6 * t175 + (qJD(1) * t28 + t9) * t219, t1 * t28 + t2 * t29 + t20 * t6 + t3 * t46 + t4 * t9 + t5 * t8; 0, 0, 0, 0, -t159, -t230, 0, 0, 0, 0, 0, t225 * t153, t225 * t156, 0, 0, 0, 0, 0, -t156 * t228 + (-qJD(1) - t217) * t242 + (-t109 * t153 - t201) * qJD(3), qJD(1) * t243 + (t110 * t153 - t141 * t233) * qJD(3) + (t203 - t245) * qJD(4), 0, 0, 0, 0, 0, t163, -t269, t163, t175 * t260 - t86 * t23 + t88 * t24 - t261 * t62, t269, -t1 * t88 - t3 * t156 + t2 * t86 + t20 * t221 + t260 * t8 + t261 * t9; 0, 0, 0, 0, 0, 0, t156 * t159 * t153, -t226 * t159, 0, 0, 0, -t156 * t230, t153 * t230, -t152 ^ 2 * t192 + (-t185 * t222 + t244) * t155 (qJD(1) * t172 - t110 * qJD(4) - t228) * t152 + ((t109 + t220) * qJD(4) + (-t173 - t195) * qJD(1)) * t155, t141 * t216 + (t141 * t237 + t156 * t186) * qJD(1), -t141 * t218 + (t156 * t187 - t203) * qJD(1), -t141 * t223, -pkin(3) * t228 - t100 * t141 + (-t173 + t201) * t134 + (-pkin(8) * t242 + t247) * qJD(4) + (t152 * t164 - t55 * t156) * qJD(1), t251 * t141 + t134 * t172 + (pkin(8) * t243 + (t102 - t256) * t155) * qJD(4) + ((pkin(3) * t218 + t56) * t156 + t164 * t155) * qJD(1), -t23 * t113 - t175 * t258, t23 * t112 - t113 * t24 - t175 * t257 + t258 * t62, -t258 * t132 + (qJD(3) * t113 - t175) * t223, -t257 * t132 + (-qJD(3) * t112 + t62) * t223, -t132 * t223, t52 * t112 + t148 * t24 + t257 * t71 + t183 * t62 + t274 * t132 + (-t10 + t249) * t223, t52 * t113 - t148 * t23 - t258 * t71 + t183 * t175 + t273 * t132 + (t11 - t248) * t223, t112 * t3 + t24 * t59 + t262 * t62 + t257 * t20 + t263 * t132 + (t8 + t249) * t223, -t1 * t112 + t113 * t2 + t174 * t23 - t175 * t263 - t24 * t76 - t257 * t9 - t258 * t8 + t264 * t62, -t113 * t3 + t23 * t59 - t262 * t175 + t258 * t20 - t264 * t132 + (-t9 + t248) * t223, t1 * t76 - t174 * t2 + t262 * t20 - t263 * t8 - t264 * t9 + t3 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110 * t109, -t109 ^ 2 + t110 ^ 2, -t109 * t141 + t162, t166 + t244, t144, -qJD(3) * t202 - t102 * t110 + t79 + (-qJD(4) + t141) * t56, -t102 * t109 + t141 * t55 - t167, t265, t277, t7, t272, t144, t132 * t14 - t275 + (-t110 * t62 - t132 * t214 + t144 * t154) * pkin(4) - t189, t132 * t15 + t279 + (-t110 * t175 - t132 * t213 - t144 * t151) * pkin(4) + t190, -t266 - t26 * t62 - t182 * t132 + (pkin(5) - t147) * t144 - t189, -t143 * t24 - t147 * t23 + (t182 + t9) * t175 + (-t250 + t8) * t62, t250 * t132 + t143 * t144 + t175 * t26 + t1 - t280, t1 * t143 + t2 * t147 + t182 * t8 - t20 * t26 + t250 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t265, t277, t7, t272, t144, t170 - t275, t171 + t279, -t32 * t62 + t170 + 0.2e1 * t184 - t266, pkin(5) * t23 - t24 * qJ(6) + (-t11 + t9) * t175 + (t8 - t229) * t62, t175 * t32 + 0.2e1 * t119 + 0.2e1 * t136 - t171 - t280, -t2 * pkin(5) + t1 * qJ(6) - t8 * t11 - t20 * t32 + t229 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t144 + t265, t7, -t132 ^ 2 - t268, -t132 * t9 + t2 + t266;];
tauc_reg  = t12;
