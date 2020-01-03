% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPR8_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR8_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR8_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:20:31
% EndTime: 2019-12-31 21:20:40
% DurationCPUTime: 4.63s
% Computational Cost: add. (5738->460), mult. (12824->565), div. (0->0), fcn. (8831->10), ass. (0->242)
t164 = sin(qJ(3));
t165 = sin(qJ(2));
t168 = cos(qJ(2));
t300 = cos(qJ(3));
t113 = t164 * t168 + t300 * t165;
t104 = t113 * qJD(1);
t307 = qJD(5) + t104;
t242 = t300 * t168;
t221 = qJD(1) * t242;
t253 = qJD(1) * t165;
t241 = t164 * t253;
t102 = -t221 + t241;
t158 = qJD(2) + qJD(3);
t163 = sin(qJ(5));
t167 = cos(qJ(5));
t83 = -t167 * t102 + t158 * t163;
t229 = t307 * t83;
t157 = qJDD(2) + qJDD(3);
t250 = qJD(5) * t167;
t251 = qJD(5) * t163;
t233 = qJDD(1) * t300;
t247 = t165 * qJDD(1);
t211 = t164 * t247 - t168 * t233;
t78 = t158 * t113;
t57 = qJD(1) * t78 + t211;
t34 = -t102 * t250 - t167 * t157 + t158 * t251 - t163 * t57;
t310 = t34 - t229;
t85 = t102 * t163 + t158 * t167;
t305 = t307 * t85;
t274 = qJD(5) * t85;
t35 = t157 * t163 - t167 * t57 + t274;
t309 = -t35 + t305;
t170 = -pkin(7) - pkin(6);
t118 = t170 * t165;
t115 = qJD(1) * t118;
t290 = qJD(2) * pkin(2);
t110 = t115 + t290;
t119 = t170 * t168;
t116 = qJD(1) * t119;
t237 = qJD(3) * t300;
t252 = qJD(3) * t164;
t248 = qJD(1) * qJD(2);
t235 = t168 * t248;
t81 = qJDD(2) * pkin(2) - t170 * (-t235 - t247);
t236 = t165 * t248;
t246 = t168 * qJDD(1);
t82 = t170 * (-t236 + t246);
t227 = t110 * t252 - t116 * t237 - t164 * t82 - t300 * t81;
t212 = qJDD(4) + t227;
t301 = pkin(3) + pkin(8);
t262 = t164 * t165;
t213 = t158 * t262;
t224 = -t158 * t221 - t164 * t246 - t165 * t233;
t56 = qJD(1) * t213 + t224;
t12 = -pkin(4) * t56 - t301 * t157 + t212;
t263 = t164 * t116;
t71 = -t300 * t110 - t263;
t257 = qJD(4) + t71;
t293 = t104 * pkin(4);
t258 = t293 + t257;
t42 = -t301 * t158 + t258;
t299 = pkin(2) * t168;
t150 = pkin(1) + t299;
t117 = t150 * qJD(1);
t190 = -qJ(4) * t104 - t117;
t44 = t301 * t102 + t190;
t20 = t163 * t42 + t167 * t44;
t99 = pkin(2) * t236 - qJDD(1) * t150;
t180 = qJ(4) * t56 - qJD(4) * t104 + t99;
t6 = t301 * t57 + t180;
t2 = -qJD(5) * t20 + t167 * t12 - t163 * t6;
t308 = t20 * t307 + t2;
t228 = -t110 * t237 - t116 * t252 - t164 * t81 + t300 * t82;
t151 = t157 * qJ(4);
t304 = -t158 * qJD(4) - t151;
t24 = t228 + t304;
t13 = -pkin(4) * t57 - t24;
t294 = t102 * pkin(4);
t108 = t300 * t116;
t72 = t164 * t110 - t108;
t68 = -t158 * qJ(4) - t72;
t48 = -t68 - t294;
t306 = t13 * t163 + t48 * t250;
t76 = t300 * t115 + t263;
t276 = pkin(2) * t237 + qJD(4) - t76;
t162 = qJ(2) + qJ(3);
t155 = sin(t162);
t156 = cos(t162);
t256 = t156 * pkin(3) + t155 * qJ(4);
t166 = sin(qJ(1));
t169 = cos(qJ(1));
t216 = g(1) * t169 + g(2) * t166;
t215 = g(1) * t166 - g(2) * t169;
t238 = qJD(2) * t300;
t202 = t170 * t238;
t223 = qJD(2) * t164 * t170;
t45 = -t118 * t237 - t119 * t252 - t165 * t202 - t168 * t223;
t87 = t164 * t118 - t300 * t119;
t303 = -t155 * t215 - t157 * t87 + t158 * t45;
t302 = t104 ^ 2;
t298 = pkin(3) * t157;
t146 = g(3) * t155;
t147 = g(3) * t156;
t295 = g(3) * t168;
t144 = t156 * pkin(8);
t292 = t85 * t83;
t291 = pkin(4) - t170;
t11 = t13 * t167;
t289 = t158 * t72;
t19 = -t163 * t44 + t167 * t42;
t288 = t163 * t19;
t287 = t163 * t35;
t52 = -qJDD(5) + t56;
t286 = t163 * t52;
t285 = t163 * t78;
t284 = t163 * t85;
t283 = t167 * t34;
t282 = t167 * t48;
t49 = t167 * t52;
t281 = t167 * t83;
t280 = t167 * t307;
t279 = t301 * t52;
t278 = t307 * t102;
t277 = t293 + t276;
t275 = pkin(6) * qJDD(1);
t273 = qJD(5) * t307;
t272 = t102 * t104;
t271 = t102 * t158;
t269 = t155 * t166;
t268 = t155 * t169;
t267 = t156 * t166;
t266 = t156 * t169;
t265 = t163 * t166;
t264 = t163 * t169;
t261 = t166 * t167;
t260 = t167 * t169;
t259 = t169 * t170;
t160 = t165 ^ 2;
t161 = t168 ^ 2;
t255 = t160 - t161;
t254 = t160 + t161;
t245 = pkin(2) * t252;
t153 = t165 * t290;
t173 = qJD(1) ^ 2;
t244 = t165 * t173 * t168;
t243 = -g(1) * t268 - g(2) * t269 + t147;
t240 = -t102 ^ 2 + t302;
t239 = -t104 * t20 - t2;
t234 = -t102 * t20 + t11;
t232 = -t35 + t274;
t230 = t307 * t48;
t75 = t115 * t164 - t108;
t86 = -t300 * t118 - t119 * t164;
t226 = t163 * t307;
t225 = -t19 * t251 + t243;
t149 = -t300 * pkin(2) - pkin(3);
t123 = t169 * t150;
t222 = g(2) * (pkin(3) * t266 + qJ(4) * t268 + t123);
t220 = t165 * t235;
t219 = -t75 + t245;
t218 = g(1) * t267 - g(2) * t266;
t217 = -pkin(2) * t165 - pkin(3) * t155;
t214 = t104 * t288 - t225;
t69 = pkin(3) * t104 + qJ(4) * t102;
t112 = -t242 + t262;
t210 = t102 * t78 + t112 * t57;
t77 = t213 + (-t237 - t238) * t168;
t209 = -t104 * t77 - t113 * t56;
t207 = t163 * t20 + t167 * t19;
t206 = t167 * t20 - t288;
t201 = -qJ(4) * t113 - t150;
t58 = t301 * t112 + t201;
t66 = pkin(4) * t113 + t86;
t31 = t163 * t66 + t167 * t58;
t30 = -t163 * t58 + t167 * t66;
t205 = t112 * t157 + t158 * t78;
t204 = t113 * t157 - t158 * t77;
t203 = t102 * t19 + t104 * t282 + t306;
t65 = pkin(2) * t253 + t69;
t200 = -t150 - t256;
t199 = -t226 * t307 - t49;
t198 = -0.2e1 * pkin(1) * t248 - pkin(6) * qJDD(2);
t197 = qJ(4) * t77 - qJD(4) * t113 + t153;
t196 = -t227 - t243;
t195 = -g(1) * t266 - g(2) * t267 - t146 - t228;
t141 = -pkin(8) + t149;
t192 = -t141 * t52 + t245 * t307;
t191 = -t280 * t307 + t286;
t189 = -t156 * t216 - t146;
t46 = qJD(3) * t87 + t165 * t223 - t168 * t202;
t188 = t157 * t86 + t158 * t46 - t218;
t187 = -t117 * t102 - t195;
t186 = t117 * t104 + t196;
t172 = qJD(2) ^ 2;
t185 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t172 + t215;
t184 = pkin(1) * t173 + t216 - t275;
t183 = t102 * t77 - t104 * t78 + t112 * t56 - t113 * t57;
t63 = pkin(3) * t102 + t190;
t182 = t104 * t63 + qJDD(4) - t196;
t181 = -t102 * t63 + t195 - t304;
t1 = qJD(5) * t19 + t12 * t163 + t167 * t6;
t179 = qJD(5) * t206 + t1 * t163 + t2 * t167;
t178 = -t141 * t273 + t189;
t177 = t273 * t301 + t189;
t39 = -t158 * t241 - t224 + t271;
t176 = t102 * t45 + t104 * t46 - t56 * t86 - t57 * t87 - t216;
t120 = qJ(4) * t267;
t122 = qJ(4) * t266;
t175 = -g(1) * t122 - g(2) * t120 - g(3) * (t144 + t256) + t216 * t155 * t301;
t143 = pkin(2) * t164 + qJ(4);
t97 = t104 * pkin(8);
t96 = -t155 * t265 + t260;
t95 = t155 * t261 + t264;
t94 = t155 * t264 + t261;
t93 = t155 * t260 - t265;
t70 = pkin(3) * t112 + t201;
t67 = -t112 * pkin(4) + t87;
t64 = -pkin(3) * t158 + t257;
t61 = t75 - t294;
t55 = t72 - t294;
t53 = t69 + t97;
t47 = t65 + t97;
t40 = -t104 * t158 + t57;
t38 = t56 - t271;
t36 = pkin(3) * t78 + t197;
t33 = -t77 * pkin(4) + t46;
t32 = -pkin(4) * t78 - t45;
t27 = t212 - t298;
t26 = t163 * t55 + t167 * t53;
t25 = -t163 * t53 + t167 * t55;
t23 = t163 * t61 + t167 * t47;
t22 = -t163 * t47 + t167 * t61;
t21 = t301 * t78 + t197;
t17 = pkin(3) * t57 + t180;
t15 = -t102 * t83 + t191;
t14 = t102 * t85 + t199;
t8 = t167 * t229 + t287;
t7 = -t226 * t85 - t283;
t5 = -qJD(5) * t31 - t163 * t21 + t167 * t33;
t4 = qJD(5) * t30 + t163 * t33 + t167 * t21;
t3 = (-t35 - t305) * t167 + (t34 + t229) * t163;
t9 = [0, 0, 0, 0, 0, qJDD(1), t215, t216, 0, 0, qJDD(1) * t160 + 0.2e1 * t220, 0.2e1 * t165 * t246 - 0.2e1 * t255 * t248, qJDD(2) * t165 + t168 * t172, qJDD(1) * t161 - 0.2e1 * t220, qJDD(2) * t168 - t165 * t172, 0, t165 * t198 + t168 * t185, -t165 * t185 + t168 * t198, 0.2e1 * t254 * t275 - t216, -g(1) * (-pkin(1) * t166 + pkin(6) * t169) - g(2) * (pkin(1) * t169 + pkin(6) * t166) + (t254 * pkin(6) ^ 2 + pkin(1) ^ 2) * qJDD(1), t209, t183, t204, t210, -t205, 0, t102 * t153 + t112 * t99 - t117 * t78 - t150 * t57 - t188, t104 * t153 + t113 * t99 + t117 * t77 + t150 * t56 + t303, t112 * t228 + t113 * t227 - t71 * t77 - t72 * t78 + t176, -t228 * t87 - t72 * t45 + t227 * t86 + t71 * t46 - t99 * t150 - t117 * t153 - g(1) * (-t150 * t166 - t259) - g(2) * (-t166 * t170 + t123), 0, -t204, t205, t209, t183, t210, t112 * t24 + t113 * t27 - t64 * t77 + t68 * t78 + t176, -t102 * t36 - t112 * t17 - t57 * t70 - t63 * t78 + t188, -t104 * t36 - t113 * t17 + t56 * t70 + t63 * t77 - t303, t17 * t70 + t63 * t36 - t24 * t87 + t68 * t45 + t27 * t86 + t64 * t46 + g(1) * t259 - t222 + (-g(1) * t200 + g(2) * t170) * t166, t78 * t284 + (-t163 * t34 + t85 * t250) * t112, (-t163 * t83 + t167 * t85) * t78 + (-t287 - t283 + (-t281 - t284) * qJD(5)) * t112, t307 * t285 - t113 * t34 - t77 * t85 + (t250 * t307 - t286) * t112, -t78 * t281 + (-t167 * t35 + t251 * t83) * t112, t78 * t280 - t113 * t35 + t77 * t83 + (-t251 * t307 - t49) * t112, -t113 * t52 - t307 * t77, -t78 * t282 - g(1) * t96 - g(2) * t94 + t113 * t2 - t19 * t77 - t30 * t52 + t32 * t83 + t35 * t67 + t5 * t307 + (t251 * t48 - t11) * t112, g(1) * t95 - g(2) * t93 - t1 * t113 + t306 * t112 + t20 * t77 + t48 * t285 - t307 * t4 + t31 * t52 + t32 * t85 - t34 * t67, t30 * t34 - t31 * t35 - t4 * t83 - t5 * t85 + t206 * t78 + (-qJD(5) * t207 + t1 * t167 - t163 * t2) * t112 + t218, t1 * t31 + t20 * t4 + t2 * t30 + t19 * t5 + t13 * t67 + t48 * t32 - t222 + (-g(1) * t291 - g(2) * t144) * t169 + (-g(1) * (t200 - t144) - g(2) * t291) * t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t244, t255 * t173, t247, t244, t246, qJDD(2), t165 * t184 - t295, g(3) * t165 + t168 * t184, 0, 0, t272, t240, t39, -t272, -t211, t157, t75 * t158 + (-t102 * t253 + t300 * t157 - t158 * t252) * pkin(2) + t186, t76 * t158 + (-t104 * t253 - t157 * t164 - t158 * t237) * pkin(2) + t187, (t72 - t75) * t104 + (t71 + t76) * t102 + (t300 * t56 - t164 * t57 + (-t300 * t102 + t104 * t164) * qJD(3)) * pkin(2), -t71 * t75 - t72 * t76 + (-t300 * t227 - t295 - t164 * t228 + (t164 * t71 + t300 * t72) * qJD(3) + (qJD(1) * t117 + t216) * t165) * pkin(2), t157, t38, t40, t272, t240, -t272, -t143 * t57 - t149 * t56 + (t219 - t68) * t104 + (t64 - t276) * t102, t102 * t65 + t219 * t158 + (-pkin(3) + t149) * t157 + t182, t104 * t65 + t143 * t157 + t276 * t158 + t181, -t24 * t143 + t27 * t149 - t63 * t65 - g(1) * (t169 * t217 + t122) - g(2) * (t166 * t217 + t120) - g(3) * (t256 + t299) - t276 * t68 + t219 * t64, t7, t3, t14, t8, t15, t278, t143 * t35 + t163 * t178 + t167 * t192 - t22 * t307 + t277 * t83 + t203, -t143 * t34 + t23 * t307 + t277 * t85 + t178 * t167 + (-t230 - t192) * t163 + t234, t22 * t85 + t23 * t83 + (t141 * t232 - t245 * t83 - t1) * t163 + (-t85 * t245 + t141 * t34 + (-t141 * t83 - t20) * qJD(5) + t239) * t167 + t214, t13 * t143 - t20 * t23 - t19 * t22 + t277 * t48 + (t165 * t216 + t207 * t252 - t295) * pkin(2) + t179 * t141 + t175; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t272, t240, t39, -t272, -t211, t157, t186 + t289, -t158 * t71 + t187, 0, 0, t157, t38, t40, t272, t240, -t272, pkin(3) * t56 - qJ(4) * t57 + (-t68 - t72) * t104 + (t64 - t257) * t102, t102 * t69 + t182 - t289 - 0.2e1 * t298, t104 * t69 + t257 * t158 + t151 + t181, -t24 * qJ(4) - t27 * pkin(3) - t63 * t69 - t64 * t72 - g(1) * (-pkin(3) * t268 + t122) - g(2) * (-pkin(3) * t269 + t120) - g(3) * t256 - t257 * t68, t7, t3, t14, t8, t15, t278, qJ(4) * t35 + t163 * t177 + t167 * t279 - t25 * t307 + t258 * t83 + t203, -qJ(4) * t34 + t26 * t307 + t258 * t85 + (-t230 - t279) * t163 + t177 * t167 + t234, t25 * t85 + t26 * t83 + (-t232 * t301 - t1) * t163 + (-t301 * t34 + (t301 * t83 - t20) * qJD(5) + t239) * t167 + t214, t13 * qJ(4) - t179 * t301 - t19 * t25 - t20 * t26 + t258 * t48 + t175; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t157 - t272, -t158 ^ 2 - t302, t158 * t68 + t182 - t298, 0, 0, 0, 0, 0, 0, -t158 * t83 + t199, -t158 * t85 + t191, t309 * t163 + t310 * t167, -t158 * t48 + (-t104 * t19 + t1) * t163 + t308 * t167 + t225; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t292, -t83 ^ 2 + t85 ^ 2, -t310, -t292, t309, -t52, -g(1) * t93 - g(2) * t95 + t147 * t167 - t48 * t85 + t308, g(1) * t94 - g(2) * t96 + t19 * t307 + t48 * t83 + (-qJD(5) * t42 - t6) * t167 + (qJD(5) * t44 - t12 - t147) * t163, 0, 0;];
tau_reg = t9;
