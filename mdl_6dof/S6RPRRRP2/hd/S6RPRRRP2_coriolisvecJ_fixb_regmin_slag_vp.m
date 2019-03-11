% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauc_reg [6x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRRP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:01:16
% EndTime: 2019-03-09 06:01:24
% DurationCPUTime: 2.88s
% Computational Cost: add. (4177->335), mult. (9757->461), div. (0->0), fcn. (6536->8), ass. (0->183)
t167 = sin(qJ(4));
t170 = cos(qJ(3));
t221 = t170 * qJD(1);
t213 = t167 * t221;
t265 = pkin(8) + pkin(9);
t214 = qJD(4) * t265;
t168 = sin(qJ(3));
t191 = pkin(3) * t168 - pkin(8) * t170;
t127 = t191 * qJD(1);
t169 = cos(qJ(4));
t152 = sin(pkin(10)) * pkin(1) + pkin(7);
t135 = t152 * qJD(1);
t274 = t170 * qJD(2) - t168 * t135;
t245 = t167 * t127 + t169 * t274;
t287 = pkin(9) * t213 - t167 * t214 - t245;
t236 = t169 * t170;
t186 = pkin(4) * t168 - pkin(9) * t236;
t197 = t169 * t127 - t167 * t274;
t286 = t186 * qJD(1) + t169 * t214 + t197;
t222 = t169 * qJD(3);
t230 = qJD(1) * t168;
t119 = -t167 * t230 + t222;
t223 = t167 * qJD(3);
t120 = t169 * t230 + t223;
t166 = sin(qJ(5));
t263 = cos(qJ(5));
t182 = t166 * t119 + t263 * t120;
t266 = t182 ^ 2;
t66 = -t263 * t119 + t166 * t120;
t64 = t66 ^ 2;
t285 = -t64 + t266;
t240 = t166 * t167;
t181 = t263 * t169 - t240;
t270 = qJD(4) + qJD(5);
t208 = t263 * qJD(5);
t271 = t263 * qJD(4) + t208;
t255 = -t271 * t169 + t181 * t221 + t270 * t240;
t122 = t166 * t169 + t263 * t167;
t72 = t270 * t122;
t254 = -t122 * t221 + t72;
t284 = t182 * t66;
t283 = t66 * qJ(6);
t219 = qJD(1) * qJD(3);
t282 = qJD(3) * qJD(4) + t170 * t219;
t210 = t170 * t223;
t225 = qJD(4) * t169;
t211 = t168 * t225;
t281 = t210 + t211;
t150 = -qJD(4) + t221;
t143 = -qJD(5) + t150;
t226 = qJD(4) * t168;
t206 = qJD(1) * t226;
t215 = t167 * t282 + t169 * t206;
t224 = qJD(5) * t166;
t80 = -t167 * t206 + t169 * t282;
t26 = -t119 * t208 + t120 * t224 + t166 * t215 - t263 * t80;
t280 = -t66 * t143 - t26;
t154 = t168 * t219;
t130 = t191 * qJD(3);
t112 = qJD(1) * t130;
t86 = t274 * qJD(3);
t201 = -t169 * t112 + t167 * t86;
t153 = -cos(pkin(10)) * pkin(1) - pkin(2);
t113 = -t170 * pkin(3) - t168 * pkin(8) + t153;
t88 = t113 * qJD(1);
t252 = t167 * t88;
t161 = t168 * qJD(2);
t93 = t170 * t135 + t161;
t85 = qJD(3) * pkin(8) + t93;
t50 = t169 * t85 + t252;
t177 = -t50 * qJD(4) - t201;
t16 = pkin(4) * t154 - t80 * pkin(9) + t177;
t217 = t167 * t112 + t169 * t86 + t88 * t225;
t227 = qJD(4) * t167;
t184 = -t85 * t227 + t217;
t19 = -t215 * pkin(9) + t184;
t49 = -t167 * t85 + t169 * t88;
t42 = -t120 * pkin(9) + t49;
t35 = -t150 * pkin(4) + t42;
t43 = t119 * pkin(9) + t50;
t199 = -t166 * t16 - t263 * t19 - t35 * t208 + t43 * t224;
t84 = -qJD(3) * pkin(3) - t274;
t63 = -t119 * pkin(4) + t84;
t279 = t63 * t66 + t199;
t36 = t66 * pkin(5) + qJD(6) + t63;
t278 = t182 * t36;
t277 = -t93 + (-t213 + t227) * pkin(4);
t276 = qJ(6) * t182;
t275 = t170 * t215;
t138 = t265 * t167;
t139 = t265 * t169;
t232 = -t166 * t138 + t263 * t139;
t273 = t232 * qJD(5) + t166 * t287 + t286 * t263;
t272 = -t138 * t208 - t139 * t224 - t286 * t166 + t263 * t287;
t41 = t263 * t43;
t11 = t166 * t35 + t41;
t176 = -t11 * qJD(5) + t263 * t16 - t166 * t19;
t269 = -t63 * t182 + t176;
t27 = t182 * qJD(5) + t166 * t80 + t263 * t215;
t268 = -t143 * t182 - t27;
t162 = t168 ^ 2;
t189 = qJD(1) * t162 - t150 * t170;
t212 = t167 * t226;
t267 = -t150 * t212 - t189 * t222;
t39 = t166 * t43;
t10 = t263 * t35 - t39;
t6 = t10 - t276;
t5 = -t143 * pkin(5) + t6;
t264 = t5 - t6;
t262 = -t254 * qJ(6) + t181 * qJD(6) + t272;
t261 = -pkin(5) * t230 + t255 * qJ(6) - t122 * qJD(6) - t273;
t228 = qJD(3) * t170;
t192 = t263 * t228;
t44 = t166 * t210 + t72 * t168 - t169 * t192;
t97 = t181 * t168;
t260 = -t97 * t27 + t44 * t66;
t239 = t167 * t168;
t45 = t167 * t192 - t166 * t212 - t224 * t239 + (t166 * t228 + t271 * t168) * t169;
t96 = t122 * t168;
t259 = t45 * t143 - t96 * t154;
t258 = t263 * t42 - t39;
t237 = t168 * t169;
t99 = t169 * t113;
t57 = -pkin(9) * t237 + t99 + (-t152 * t167 - pkin(4)) * t170;
t126 = t152 * t236;
t246 = t167 * t113 + t126;
t62 = -pkin(9) * t239 + t246;
t256 = t166 * t57 + t263 * t62;
t253 = t167 * t84;
t251 = t169 * t84;
t250 = t80 * t167;
t87 = qJD(3) * t161 + t135 * t228;
t249 = t87 * t167;
t248 = t87 * t169;
t247 = t113 * t225 + t167 * t130;
t244 = t119 * t150;
t243 = t119 * t168;
t242 = t120 * t150;
t241 = t150 * t169;
t238 = t167 * t170;
t137 = t168 * t152;
t171 = qJD(3) ^ 2;
t235 = t171 * t168;
t234 = t171 * t170;
t233 = t169 * t130 + t223 * t137;
t102 = pkin(4) * t239 + t137;
t231 = -t170 ^ 2 + t162;
t136 = qJD(1) * t153;
t229 = qJD(3) * t168;
t73 = pkin(4) * t281 + t152 * t228;
t160 = -t169 * pkin(4) - pkin(3);
t205 = -t166 * t42 - t41;
t203 = -t166 * t62 + t263 * t57;
t200 = t26 * t170 + t182 * t229;
t198 = t120 * t229 - t80 * t170;
t196 = t150 * t152 + t85;
t195 = -t263 * t138 - t166 * t139;
t193 = t150 * t211;
t190 = t182 * t45 - t96 * t26;
t187 = 0.2e1 * qJD(3) * t136;
t185 = t170 * t27 - t66 * t229;
t30 = t186 * qJD(3) + (-t126 + (pkin(9) * t168 - t113) * t167) * qJD(4) + t233;
t32 = (-t168 * t222 - t170 * t227) * t152 - t281 * pkin(9) + t247;
t183 = t166 * t30 + t57 * t208 - t62 * t224 + t263 * t32;
t179 = t189 * t167;
t52 = t215 * pkin(4) + t87;
t178 = -t44 * t143 - t97 * t154;
t12 = t27 * pkin(5) + t52;
t175 = -t256 * qJD(5) - t166 * t32 + t263 * t30;
t172 = qJD(1) ^ 2;
t159 = t263 * pkin(4) + pkin(5);
t60 = qJ(6) * t181 + t232;
t59 = -t122 * qJ(6) + t195;
t22 = -t96 * qJ(6) + t256;
t20 = -t170 * pkin(5) - t97 * qJ(6) + t203;
t9 = t258 - t276;
t8 = t205 + t283;
t7 = t11 - t283;
t4 = -t45 * qJ(6) - t96 * qJD(6) + t183;
t3 = pkin(5) * t229 + t44 * qJ(6) - t97 * qJD(6) + t175;
t2 = -t27 * qJ(6) - t66 * qJD(6) - t199;
t1 = pkin(5) * t154 + t26 * qJ(6) - qJD(6) * t182 + t176;
t13 = [0, 0, 0, 0, 0.2e1 * t170 * t154, -0.2e1 * t231 * t219, t234, -t235, 0, -t152 * t234 + t168 * t187, t152 * t235 + t170 * t187, t80 * t237 + (t170 * t222 - t212) * t120 (t169 * t119 - t120 * t167) * t228 + (-t169 * t215 - t250 + (-t167 * t119 - t120 * t169) * qJD(4)) * t168, t198 - t267, t193 + t275 + (-t179 + t243) * qJD(3) (-t150 - t221) * t229 -(-t113 * t227 + t233) * t150 + ((-t119 * t152 + t253) * qJD(3) + (t196 * t169 + t252) * qJD(4) + t201) * t170 + (t152 * t215 + t249 + t84 * t225 + ((-t152 * t238 + t99) * qJD(1) + t49) * qJD(3)) * t168, t247 * t150 + (-t196 * t227 + (t120 * t152 + t251) * qJD(3) + t217) * t170 + (-t84 * t227 + t152 * t80 + t248 + (-t246 * qJD(1) - t152 * t241 - t50) * qJD(3)) * t168, -t182 * t44 - t26 * t97, -t190 + t260, -t178 + t200, t185 + t259 (-t143 - t221) * t229, t10 * t229 + t102 * t27 - t143 * t175 + t154 * t203 - t170 * t176 + t63 * t45 + t52 * t96 + t73 * t66, t183 * t143 - t199 * t170 + t73 * t182 - t102 * t26 + t52 * t97 - t63 * t44 + (-t256 * qJD(1) - t11) * t229, -t1 * t97 - t182 * t3 - t2 * t96 + t20 * t26 - t22 * t27 - t4 * t66 + t5 * t44 - t7 * t45, t2 * t22 + t7 * t4 + t1 * t20 + t5 * t3 + t12 * (t96 * pkin(5) + t102) + t36 * (t45 * pkin(5) + t73); 0, 0, 0, 0, 0, 0, 0, 0, 0, -t235, -t234, 0, 0, 0, 0, 0, t193 - t275 + (-t179 - t243) * qJD(3), t198 + t267, 0, 0, 0, 0, 0, -t185 + t259, t178 + t200, t190 + t260, -t1 * t96 - t12 * t170 + t2 * t97 + t229 * t36 - t7 * t44 - t5 * t45; 0, 0, 0, 0, -t170 * t172 * t168, t231 * t172, 0, 0, 0, t93 * qJD(3) - t136 * t230 - t87, -t136 * t221, -t120 * t241 + t250 (t80 - t244) * t169 + (-t215 + t242) * t167, -t150 * t225 + (t150 * t236 + (-t120 + t223) * t168) * qJD(1), t150 * t227 + (-t150 * t238 + (-t119 + t222) * t168) * qJD(1), t150 * t230, -pkin(3) * t215 - t248 + t197 * t150 + t93 * t119 + (pkin(8) * t241 + t253) * qJD(4) + (-t49 * t168 + (-pkin(8) * t229 - t170 * t84) * t167) * qJD(1), -pkin(3) * t80 + t249 - t245 * t150 - t93 * t120 + (-t167 * pkin(8) * t150 + t251) * qJD(4) + (-t84 * t236 + (-pkin(8) * t222 + t50) * t168) * qJD(1), -t26 * t122 - t182 * t255, -t122 * t27 - t181 * t26 - t182 * t254 + t255 * t66, t255 * t143 + (qJD(3) * t122 - t182) * t230, t254 * t143 + (qJD(3) * t181 + t66) * t230, t143 * t230, -t10 * t230 + t273 * t143 + t154 * t195 + t160 * t27 - t181 * t52 + t254 * t63 + t277 * t66, t52 * t122 - t160 * t26 + t277 * t182 - t255 * t63 + t272 * t143 + (-qJD(3) * t232 + t11) * t230, -t1 * t122 + t181 * t2 - t182 * t261 - t254 * t7 + t255 * t5 + t59 * t26 - t262 * t66 - t60 * t27, t2 * t60 + t1 * t59 + t12 * (-pkin(5) * t181 + t160) + t262 * t7 + t261 * t5 + (t254 * pkin(5) + t277) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120 * t119, -t119 ^ 2 + t120 ^ 2, t80 + t244, -t215 - t242, t154, -t84 * t120 - t50 * t150 + t177, -t84 * t119 - t49 * t150 - t184, t284, t285, t280, t268, t154, t205 * t143 + (-t120 * t66 + t143 * t224 + t263 * t154) * pkin(4) + t269, -t258 * t143 + (-t120 * t182 + t143 * t208 - t154 * t166) * pkin(4) + t279, t159 * t26 - t5 * t66 + t7 * t182 + t9 * t66 + t8 * t182 + (-t166 * t27 + (t166 * t182 - t263 * t66) * qJD(5)) * pkin(4), -pkin(5) * t278 + t1 * t159 - t5 * t8 - t7 * t9 + (-t36 * t120 + t2 * t166 + (-t166 * t5 + t263 * t7) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t284, t285, t280, t268, t154, -t11 * t143 + t269, -t10 * t143 + t279, pkin(5) * t26 - t264 * t66, t264 * t7 + (t1 - t278) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64 - t266, t182 * t5 + t7 * t66 + t12;];
tauc_reg  = t13;
