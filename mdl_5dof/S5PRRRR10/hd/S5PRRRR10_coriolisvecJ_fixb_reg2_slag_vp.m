% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRR10_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR10_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:26:06
% EndTime: 2019-12-05 17:26:26
% DurationCPUTime: 6.89s
% Computational Cost: add. (6442->446), mult. (18097->670), div. (0->0), fcn. (14689->12), ass. (0->216)
t147 = sin(pkin(6));
t153 = sin(qJ(3));
t247 = t147 * t153;
t140 = pkin(8) * t247;
t149 = cos(pkin(6));
t157 = cos(qJ(3));
t242 = t149 * t157;
t115 = pkin(2) * t242 - t140;
t109 = qJD(3) * t115;
t158 = cos(qJ(2));
t235 = t157 * t158;
t154 = sin(qJ(2));
t240 = t153 * t154;
t168 = -t149 * t240 + t235;
t148 = sin(pkin(5));
t233 = qJD(1) * t148;
t252 = -t168 * t233 + t109;
t184 = pkin(3) * t153 - pkin(9) * t157;
t166 = t184 * qJD(3);
t206 = t154 * t233;
t291 = (t166 - t206) * t147;
t156 = cos(qJ(4));
t290 = pkin(9) * t156;
t243 = t149 * t153;
t246 = t147 * t157;
t116 = pkin(2) * t243 + pkin(8) * t246;
t105 = t149 * pkin(9) + t116;
t185 = -pkin(3) * t157 - pkin(9) * t153;
t106 = (-pkin(2) + t185) * t147;
t152 = sin(qJ(4));
t223 = qJD(4) * t156;
t225 = qJD(4) * t152;
t264 = -t105 * t225 + t106 * t223 + t291 * t152 + t252 * t156;
t110 = qJD(3) * t116;
t238 = t154 * t157;
t239 = t153 * t158;
t170 = t149 * t238 + t239;
t251 = -t170 * t233 + t110;
t151 = sin(qJ(5));
t155 = cos(qJ(5));
t228 = qJD(2) * t157;
t138 = t147 * t228;
t180 = t138 - qJD(4);
t219 = t149 * qJD(2);
t197 = qJD(3) + t219;
t150 = cos(pkin(5));
t232 = qJD(1) * t150;
t261 = qJD(2) * pkin(2);
t128 = t158 * t233 + t261;
t244 = t149 * t128;
t165 = t147 * t232 + t244;
t231 = qJD(2) * t147;
t118 = pkin(8) * t231 + t206;
t236 = t157 * t118;
t63 = t153 * t165 + t236;
t52 = pkin(9) * t197 + t63;
t137 = t149 * t232;
t69 = t137 + (qJD(2) * t185 - t128) * t147;
t26 = t152 * t69 + t156 * t52;
t24 = -pkin(10) * t180 + t26;
t111 = t153 * t118;
t62 = t157 * t165 - t111;
t51 = -pkin(3) * t197 - t62;
t229 = qJD(2) * t153;
t210 = t147 * t229;
t95 = t152 * t210 - t156 * t197;
t97 = t152 * t197 + t156 * t210;
t27 = t95 * pkin(4) - t97 * pkin(10) + t51;
t179 = t151 * t24 - t155 * t27;
t163 = t170 * qJD(2);
t171 = t128 * t243 + t236;
t226 = qJD(3) * t153;
t209 = t147 * t226;
t190 = t150 * t209;
t41 = t171 * qJD(3) + (t148 * t163 + t190) * qJD(1);
t186 = qJD(3) * t138;
t285 = qJD(4) * t95;
t74 = -t156 * t186 + t285;
t174 = t152 * t186;
t250 = qJD(4) * t97;
t75 = t174 + t250;
t18 = t75 * pkin(4) + t74 * pkin(10) + t41;
t162 = t168 * qJD(2);
t227 = qJD(3) * t147;
t208 = t157 * t227;
t189 = t150 * t208;
t40 = (t128 * t242 - t111) * qJD(3) + (t148 * t162 + t189) * qJD(1);
t78 = (t166 + t206) * t231;
t167 = -t152 * t78 - t156 * t40 - t69 * t223 + t225 * t52;
t217 = qJD(2) * qJD(3);
t205 = t153 * t217;
t187 = t147 * t205;
t7 = pkin(10) * t187 - t167;
t1 = -qJD(5) * t179 + t151 * t18 + t155 * t7;
t92 = qJD(5) + t95;
t289 = t179 * t92 + t1;
t288 = pkin(10) * t209 + t264;
t211 = t152 * t247;
t82 = qJD(4) * t211 - t149 * t223 - t156 * t208;
t114 = t152 * t149 + t156 * t247;
t83 = qJD(4) * t114 + t152 * t208;
t287 = t83 * pkin(4) + t82 * pkin(10) + t251;
t198 = qJD(4) * t180;
t286 = -pkin(9) * t198 + t41;
t25 = -t152 * t52 + t156 * t69;
t284 = t180 * t25 - t167;
t183 = pkin(4) * t152 - pkin(10) * t156;
t283 = qJD(5) * t290 - t183 * qJD(4) + (t153 * t232 + t183 * t228) * t147 + t171;
t10 = t151 * t27 + t155 * t24;
t2 = -qJD(5) * t10 - t151 * t7 + t155 * t18;
t281 = -t10 * t92 - t2;
t262 = t156 * t105 + t152 * t106;
t263 = -qJD(4) * t262 - t252 * t152 + t291 * t156;
t279 = t180 * t95;
t278 = t180 * t97;
t73 = -t151 * t180 + t155 * t97;
t249 = qJD(5) * t73;
t33 = -t151 * t74 - t155 * t187 + t249;
t277 = t149 * t235 - t240;
t204 = t152 * t40 - t156 * t78;
t12 = -qJD(4) * t26 - t204;
t159 = qJD(2) ^ 2;
t104 = t140 + (-pkin(2) * t157 - pkin(3)) * t149;
t113 = -t156 * t149 + t211;
t55 = t113 * pkin(4) - t114 * pkin(10) + t104;
t57 = -pkin(10) * t246 + t262;
t21 = -t151 * t57 + t155 * t55;
t275 = qJD(5) * t21 + t287 * t151 + t288 * t155;
t22 = t151 * t55 + t155 * t57;
t274 = -qJD(5) * t22 - t288 * t151 + t287 * t155;
t80 = -t148 * t277 - t150 * t246;
t272 = t41 * t80;
t201 = t155 * t180;
t71 = t151 * t97 + t201;
t271 = t71 * t92;
t270 = t73 * t71;
t269 = t73 * t92;
t268 = t97 * t95;
t135 = -t156 * pkin(4) - t152 * pkin(10) - pkin(3);
t222 = qJD(5) * t151;
t107 = t184 * t231;
t47 = t152 * t107 + t156 * t62;
t38 = pkin(10) * t210 + t47;
t267 = t135 * t222 + (-t225 * pkin(9) - t38) * t151 + t283 * t155;
t221 = qJD(5) * t155;
t224 = qJD(4) * t155;
t266 = t152 * t224 * pkin(9) - t135 * t221 + t283 * t151 + t155 * t38;
t265 = -pkin(4) * t209 - t263;
t91 = -t147 * t128 + t137;
t260 = t147 * t91;
t259 = t151 * t75;
t258 = t151 * t92;
t257 = t155 * t75;
t32 = qJD(5) * t201 - t151 * t187 + t155 * t74 + t222 * t97;
t256 = t32 * t151;
t255 = t33 * t155;
t254 = t75 * t113;
t253 = t75 * t156;
t144 = t147 ^ 2;
t248 = t144 * t159;
t245 = t148 * t159;
t241 = t151 * t156;
t237 = t156 * t157;
t234 = t153 ^ 2 - t157 ^ 2;
t230 = qJD(2) * t148;
t218 = t156 * qJD(3);
t215 = t144 * t261;
t214 = t154 * t245;
t212 = t151 * t246;
t207 = t147 * t149 * t159;
t202 = t155 * t92;
t200 = t157 * t180;
t199 = t180 * t147;
t196 = qJD(3) + 0.2e1 * t219;
t195 = t144 * t214;
t194 = t153 * t157 * t248;
t8 = -pkin(4) * t187 - t12;
t193 = pkin(10) * qJD(5) * t92 + t8;
t191 = t147 * t154 * t230;
t182 = -t10 * t151 + t155 * t179;
t112 = -t148 * t158 * t147 + t150 * t149;
t169 = t149 * t239 + t238;
t81 = t148 * t169 + t150 * t247;
t54 = t112 * t152 + t81 * t156;
t31 = t80 * t151 + t54 * t155;
t30 = -t54 * t151 + t80 * t155;
t46 = t156 * t107 - t152 * t62;
t178 = t112 * t156 - t81 * t152;
t64 = -t152 * t105 + t156 * t106;
t175 = t144 * t157 * t205;
t173 = -t215 + t260;
t23 = pkin(4) * t180 - t25;
t172 = -pkin(10) * t75 + t23 * t92;
t84 = t151 * t114 + t155 * t246;
t164 = t152 * t180;
t161 = qJD(3) * t118 + t206 * t219;
t160 = -t91 * t231 - qJD(3) * t244 + (-t150 * t227 - t158 * t230) * qJD(1);
t102 = t151 * t135 + t155 * t290;
t101 = -pkin(9) * t241 + t155 * t135;
t90 = (t151 * t153 + t155 * t237) * t231;
t89 = t138 * t241 - t155 * t210;
t85 = t155 * t114 - t212;
t61 = t97 * pkin(4) + t95 * pkin(10);
t56 = pkin(4) * t246 - t64;
t50 = t189 + (qJD(3) * t277 + t162) * t148;
t49 = t190 + (qJD(3) * t169 + t163) * t148;
t43 = -qJD(5) * t212 + t114 * t221 - t151 * t82 - t155 * t209;
t42 = qJD(5) * t84 - t151 * t209 + t155 * t82;
t37 = -pkin(4) * t210 - t46;
t20 = qJD(4) * t178 + t152 * t191 + t50 * t156;
t19 = qJD(4) * t54 + t50 * t152 - t156 * t191;
t14 = t151 * t61 + t155 * t25;
t13 = -t151 * t25 + t155 * t61;
t6 = qJD(5) * t30 + t49 * t151 + t20 * t155;
t5 = -qJD(5) * t31 - t20 * t151 + t49 * t155;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t214, -t158 * t245, 0, 0, 0, 0, 0, 0, 0, 0, t112 * t187 - t157 * t195 - t197 * t49, t112 * t186 + t153 * t195 - t197 * t50, (t153 * t49 + t157 * t50 + (-t153 * t81 + t157 * t80) * qJD(3)) * t231, t40 * t81 + t272 - t62 * t49 + t63 * t50 + (qJD(1) * t112 + t91) * t191, 0, 0, 0, 0, 0, 0, t178 * t187 + t180 * t19 + t49 * t95 + t80 * t75, t180 * t20 - t187 * t54 + t49 * t97 - t80 * t74, t178 * t74 + t19 * t97 - t20 * t95 - t54 * t75, t12 * t178 - t167 * t54 - t25 * t19 + t26 * t20 + t51 * t49 + t272, 0, 0, 0, 0, 0, 0, -t178 * t33 + t19 * t71 + t30 * t75 + t5 * t92, t178 * t32 + t19 * t73 - t31 * t75 - t6 * t92, t30 * t32 - t31 * t33 - t5 * t73 - t6 * t71, t1 * t31 + t10 * t6 - t178 * t8 - t179 * t5 + t19 * t23 + t2 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t175, -0.2e1 * t234 * t144 * t217, t196 * t208, -0.2e1 * t175, -t196 * t209, 0, (-qJD(2) * t251 - t41) * t149 + (t153 * t173 - t251) * qJD(3), (-qJD(2) * t252 - t40) * t149 + (t157 * t173 - t252) * qJD(3), (t153 * t41 + t157 * t40 + (-t153 * t63 - t157 * t62) * qJD(3) + ((t252 - t109) * t157 + (t251 - t110) * t153) * qJD(2)) * t147, -t41 * t115 + t40 * t116 + t252 * t63 - t251 * t62 + (-t215 - t260) * t206, -t74 * t114 - t97 * t82, t74 * t113 - t114 * t75 + t82 * t95 - t97 * t83, t82 * t180 + (t74 * t157 + (qJD(2) * t114 + t97) * t226) * t147, t95 * t83 + t254, t83 * t180 + (t75 * t157 + (-qJD(2) * t113 - t95) * t226) * t147, (-t144 * t228 - t199) * t226, t104 * t75 + t41 * t113 + t51 * t83 + t251 * t95 + (-t12 * t157 + (qJD(2) * t64 + t25) * t226) * t147 - t263 * t180, -t104 * t74 + t41 * t114 - t51 * t82 + t251 * t97 + (-t167 * t157 + (-qJD(2) * t262 - t26) * t226) * t147 + t264 * t180, t113 * t167 - t12 * t114 + t25 * t82 - t26 * t83 - t262 * t75 - t263 * t97 - t264 * t95 + t64 * t74, t41 * t104 + t12 * t64 - t167 * t262 + t25 * t263 + t251 * t51 + t26 * t264, -t32 * t85 - t73 * t42, t32 * t84 - t85 * t33 + t42 * t71 - t73 * t43, -t32 * t113 - t42 * t92 + t73 * t83 + t85 * t75, t33 * t84 + t71 * t43, -t33 * t113 - t43 * t92 - t71 * t83 - t84 * t75, t92 * t83 + t254, t2 * t113 - t179 * t83 + t21 * t75 + t23 * t43 + t265 * t71 + t274 * t92 + t56 * t33 + t8 * t84, -t1 * t113 - t10 * t83 - t22 * t75 - t23 * t42 + t265 * t73 - t275 * t92 - t56 * t32 + t8 * t85, -t1 * t84 - t10 * t43 - t179 * t42 - t2 * t85 + t21 * t32 - t22 * t33 - t274 * t73 - t275 * t71, t1 * t22 + t10 * t275 - t179 * t274 + t2 * t21 + t265 * t23 + t8 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t194, t234 * t248, -t157 * t207, t194, t153 * t207, 0, t153 * t160 - t157 * t161 + t197 * t63, t153 * t161 + t157 * t160 + t197 * t62, 0, 0, -t74 * t152 - t156 * t278, (-t74 + t279) * t156 + (-t75 + t278) * t152, -t156 * t198 + (t156 * t200 + (t152 * qJD(3) - t97) * t153) * t231, -t164 * t95 - t253, t152 * t198 + (-t152 * t200 + (t95 + t218) * t153) * t231, t199 * t229, -pkin(3) * t75 + t51 * t225 + t46 * t180 - t63 * t95 - t286 * t156 + (-t153 * t25 + (-pkin(9) * t226 - t157 * t51) * t152) * t231, pkin(3) * t74 + t51 * t223 - t47 * t180 - t63 * t97 + t286 * t152 + (-t51 * t237 + (-pkin(9) * t218 + t26) * t153) * t231, t46 * t97 + t47 * t95 + ((-t75 + t250) * pkin(9) + t284) * t156 + (-t12 + t180 * t26 + (-t74 + t285) * pkin(9)) * t152, -t41 * pkin(3) - t25 * t46 - t26 * t47 - t51 * t63 + (-t167 * t156 - t12 * t152 + (-t152 * t26 - t156 * t25) * qJD(4)) * pkin(9), -t32 * t155 * t152 + (-t152 * t222 + t155 * t223 - t90) * t73, t90 * t71 + t73 * t89 + (-t151 * t73 - t155 * t71) * t223 + (t256 - t255 + (t151 * t71 - t155 * t73) * qJD(5)) * t152, -t90 * t92 + (t224 * t92 + t32) * t156 + (-t180 * t73 - t222 * t92 + t257) * t152, t33 * t151 * t152 + (t151 * t223 + t152 * t221 - t89) * t71, t89 * t92 + (-qJD(4) * t258 + t33) * t156 + (t180 * t71 - t221 * t92 - t259) * t152, -t164 * t92 - t253, t101 * t75 - t23 * t89 - t37 * t71 - t267 * t92 + (-t2 + (pkin(9) * t71 + t151 * t23) * qJD(4)) * t156 + (pkin(9) * t33 + t8 * t151 + t179 * t180 + t221 * t23) * t152, -t102 * t75 - t23 * t90 - t37 * t73 + t266 * t92 + (t1 + (pkin(9) * t73 + t155 * t23) * qJD(4)) * t156 + (-pkin(9) * t32 + t10 * t180 + t8 * t155 - t222 * t23) * t152, t10 * t89 + t101 * t32 - t102 * t33 - t179 * t90 + t267 * t73 + t266 * t71 + t182 * t223 + (-t1 * t151 - t155 * t2 + (-t10 * t155 - t151 * t179) * qJD(5)) * t152, t1 * t102 + t2 * t101 - t23 * t37 + t267 * t179 - t266 * t10 + (t152 * t8 + t223 * t23) * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t268, -t95 ^ 2 + t97 ^ 2, -t74 - t279, -t268, -t97 * t138 - t174, t187, -t26 * t138 - t51 * t97 - t204, t51 * t95 - t284, 0, 0, t202 * t73 - t256, (-t32 - t271) * t155 + (-t33 - t269) * t151, t202 * t92 - t73 * t97 + t259, t258 * t71 - t255, -t151 * t92 ^ 2 + t71 * t97 + t257, -t92 * t97, -pkin(4) * t33 - t13 * t92 + t151 * t172 - t155 * t193 + t179 * t97 - t26 * t71, pkin(4) * t32 + t10 * t97 + t14 * t92 + t151 * t193 + t155 * t172 - t26 * t73, t13 * t73 + t14 * t71 + ((-t33 + t249) * pkin(10) + t289) * t155 + ((qJD(5) * t71 - t32) * pkin(10) + t281) * t151, -t8 * pkin(4) - t10 * t14 + t179 * t13 - t23 * t26 + (qJD(5) * t182 + t1 * t155 - t2 * t151) * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t270, -t71 ^ 2 + t73 ^ 2, -t32 + t271, -t270, t269 - t33, t75, -t23 * t73 - t281, t23 * t71 - t289, 0, 0;];
tauc_reg = t3;