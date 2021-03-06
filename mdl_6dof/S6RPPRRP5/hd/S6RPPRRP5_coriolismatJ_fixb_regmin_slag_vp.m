% Calculate minimal parameter regressor of coriolis matrix for
% S6RPPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x25]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPPRRP5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:08:48
% EndTime: 2019-03-09 02:08:55
% DurationCPUTime: 2.74s
% Computational Cost: add. (2320->309), mult. (4093->419), div. (0->0), fcn. (3316->4), ass. (0->236)
t177 = sin(qJ(5));
t313 = 0.2e1 * t177;
t179 = cos(qJ(5));
t180 = cos(qJ(4));
t266 = t179 * t180;
t203 = t266 * t313;
t293 = -qJ(6) - pkin(8);
t122 = t293 * t179;
t178 = sin(qJ(4));
t175 = -pkin(7) + qJ(2);
t272 = t175 * t177;
t205 = pkin(5) - t272;
t221 = qJ(6) * t266;
t176 = pkin(1) + qJ(3);
t199 = t178 * pkin(4) - t180 * pkin(8);
t119 = t199 + t176;
t94 = t179 * t119;
t55 = t178 * t205 - t221 + t94;
t311 = -t55 / 0.2e1;
t268 = t178 * t175;
t70 = t177 * t268 - t94;
t59 = -t70 - t221;
t224 = t59 / 0.2e1 + t311;
t312 = t224 * t122;
t172 = t178 ^ 2;
t174 = t180 ^ 2;
t251 = t172 + t174;
t112 = t251 * t177;
t171 = t177 ^ 2;
t173 = t179 ^ 2;
t134 = t171 + t173;
t294 = t180 * pkin(4);
t296 = t178 * pkin(8);
t123 = t294 + t296;
t117 = t179 * t123;
t267 = t178 * t179;
t57 = qJ(6) * t267 + t180 * t205 + t117;
t310 = t57 / 0.2e1;
t308 = t171 / 0.2e1;
t307 = t172 / 0.2e1;
t306 = -t173 / 0.2e1;
t305 = t174 / 0.2e1;
t304 = -t177 / 0.2e1;
t302 = t178 / 0.2e1;
t301 = t179 / 0.2e1;
t300 = -t180 / 0.2e1;
t299 = t180 / 0.2e1;
t298 = t177 * pkin(5);
t297 = t178 * pkin(5);
t295 = t179 * pkin(5);
t292 = pkin(5) * qJD(5);
t291 = t55 * t177;
t290 = t57 * t179;
t227 = -t297 / 0.2e1;
t6 = (t227 + t224) * t179;
t289 = t6 * qJD(1);
t269 = t177 * t180;
t222 = t175 * t267;
t71 = t177 * t119 + t222;
t60 = -qJ(6) * t269 + t71;
t288 = t60 * t177;
t287 = t60 * t179;
t116 = t177 * t123;
t264 = t180 * t175;
t125 = t179 * t264;
t270 = t177 * t178;
t62 = qJ(6) * t270 + t116 + t125;
t286 = t62 * t177;
t219 = t288 / 0.2e1;
t32 = t55 * t267;
t44 = t60 * t270;
t225 = -t32 / 0.2e1 - t44 / 0.2e1;
t228 = -t174 / 0.2e1 - 0.1e1 / 0.2e1;
t8 = t178 * t219 + (t228 * pkin(5) + t59 * t302) * t179 + t225;
t285 = t8 * qJD(1);
t204 = -t175 + t298;
t103 = t204 * t178;
t104 = t204 * t180;
t9 = -t104 * t103 + t55 * t57 + t60 * t62;
t284 = t9 * qJD(1);
t31 = t55 * t269;
t10 = t59 * t269 - t31;
t283 = t10 * qJD(1);
t282 = t104 * t177;
t11 = -t32 - t44 + (t286 + t290) * t180;
t281 = t11 * qJD(1);
t126 = t177 * t227;
t12 = t177 * t224 + t126;
t280 = t12 * qJD(1);
t121 = t293 * t177;
t279 = t121 * t177;
t278 = t121 * t178;
t277 = t122 * t178;
t276 = t122 * t179;
t230 = pkin(5) * t266;
t14 = t104 * t230 + (-t55 + t59) * t60;
t275 = t14 * qJD(1);
t157 = -pkin(4) - t295;
t265 = t180 * t157;
t208 = -t265 / 0.2e1;
t17 = t208 + (-t277 / 0.2e1 + t310) * t179 + (-t278 / 0.2e1 + t62 / 0.2e1) * t177;
t274 = t17 * qJD(1);
t273 = t174 * t179;
t271 = t177 * t122;
t19 = -t104 * t180 + (t287 - t291) * t178;
t263 = t19 * qJD(1);
t21 = t55 * t179 + t288;
t20 = t21 * t180;
t262 = t20 * qJD(1);
t261 = t21 * qJD(1);
t220 = t177 * t264;
t22 = -t70 * t180 + (t117 + t220) * t178;
t260 = t22 * qJD(1);
t23 = t71 * t180 + (-t125 + t116) * t178;
t259 = t23 * qJD(1);
t45 = -t174 * t272 - t178 * t70;
t258 = t45 * qJD(1);
t46 = -t175 * t273 - t178 * t71;
t257 = t46 * qJD(1);
t68 = t305 + t134 * (t307 + 0.1e1 / 0.2e1);
t256 = t68 * qJD(1);
t201 = t307 - t228;
t76 = t201 * t177;
t255 = t76 * qJD(1);
t78 = t201 * t179;
t254 = t78 * qJD(1);
t79 = t171 * t300 + (-0.1e1 / 0.2e1 + t306) * t180;
t253 = t79 * qJD(1);
t232 = t180 * qJD(4);
t149 = t177 * t232;
t245 = qJD(5) * t179;
t152 = t178 * t245;
t252 = t149 + t152;
t135 = t172 - t174;
t136 = t173 - t171;
t250 = qJD(2) * t178;
t249 = qJD(3) * t178;
t248 = qJD(4) * t177;
t247 = qJD(4) * t179;
t246 = qJD(5) * t177;
t155 = t171 * t178;
t156 = t173 * t178;
t108 = -t155 - t156;
t244 = t108 * qJD(1);
t109 = t134 * t174;
t243 = t109 * qJD(1);
t111 = t134 * t180;
t242 = t111 * qJD(1);
t241 = t112 * qJD(1);
t113 = t135 * t177;
t240 = t113 * qJD(1);
t114 = t251 * t179;
t239 = t114 * qJD(1);
t115 = t135 * t179;
t238 = t115 * qJD(1);
t237 = t134 * qJD(4);
t236 = t135 * qJD(1);
t235 = t176 * qJD(1);
t163 = t178 * qJD(1);
t234 = t178 * qJD(4);
t164 = t180 * qJD(1);
t233 = t180 * qJD(2);
t231 = t180 * qJD(5);
t226 = t295 / 0.2e1;
t218 = t177 * t231;
t217 = t179 * t231;
t216 = t177 * t245;
t215 = t177 * t247;
t214 = t177 * t164;
t213 = t179 * t232;
t212 = t179 * t164;
t211 = t178 * t232;
t210 = t178 * t164;
t209 = -t273 / 0.2e1;
t207 = t155 / 0.2e1 + t156 / 0.2e1;
t206 = t173 / 0.2e1 + t308;
t202 = -qJD(2) + t235;
t200 = -t122 * t300 + t311;
t198 = qJD(4) * t203;
t151 = t178 * t246;
t197 = t151 - t213;
t181 = (t62 * t301 + t57 * t304 + t104 / 0.2e1) * t178 - t31 / 0.2e1;
t4 = -t103 * t300 + t271 / 0.2e1 + (t60 * t299 - t121 / 0.2e1) * t179 + t181;
t72 = (-0.1e1 + t134) * t180 * t178;
t196 = t4 * qJD(1) + t72 * qJD(3);
t1 = t312 + (t310 + t179 * t208 - t282 / 0.2e1) * pkin(5);
t24 = t157 * t298;
t195 = -qJD(1) * t1 + qJD(4) * t24;
t67 = -t276 - t279;
t194 = (-qJD(5) - t163) * t180;
t193 = t296 / 0.2e1 + t294 / 0.2e1;
t187 = t193 * t179;
t64 = -t117 / 0.2e1 - t187;
t192 = pkin(4) * t248 - qJD(1) * t64;
t186 = t193 * t177;
t63 = t116 / 0.2e1 + t186;
t191 = pkin(4) * t247 - qJD(1) * t63;
t190 = (t60 / 0.2e1 + t121 * t300) * t179;
t99 = (t308 + t306) * t180;
t189 = -qJD(1) * t99 + t215;
t188 = t179 * t194;
t185 = qJD(1) * t177 * t273 + qJD(4) * t99;
t110 = t136 * t174;
t184 = qJD(1) * t110 + t198;
t183 = qJD(1) * t203 - qJD(4) * t136;
t16 = -t268 / 0.2e1 + t190 + (t297 / 0.2e1 + t200) * t177;
t77 = t302 - t207;
t182 = qJD(1) * t16 - qJD(3) * t77 + qJD(4) * t67;
t170 = qJ(2) * qJD(1);
t169 = qJ(2) * qJD(2);
t159 = pkin(5) * t246;
t154 = t232 / 0.2e1;
t150 = t179 * t163;
t148 = t177 * t234;
t147 = t177 * t163;
t107 = t150 + t245;
t106 = t147 + t246;
t105 = (t163 + qJD(5) / 0.2e1) * t180;
t93 = (t212 + t248) * pkin(5);
t92 = t99 * qJD(5);
t83 = -t179 * t172 / 0.2e1 + t209 + t301;
t82 = t302 + t207;
t81 = t304 + t112 / 0.2e1;
t80 = t173 * t299 + (t308 - 0.1e1 / 0.2e1) * t180;
t69 = t172 * t206 - t206 + t305;
t43 = -t220 + t117 / 0.2e1 - t187;
t42 = -t125 - t116 / 0.2e1 + t186;
t36 = pkin(5) * t269;
t18 = -t286 / 0.2e1 - t290 / 0.2e1 + t208 + (-t276 / 0.2e1 - t279 / 0.2e1) * t178;
t15 = t268 / 0.2e1 + t126 + t190 + t200 * t177;
t13 = t59 * t304 + t291 / 0.2e1 + t126;
t7 = pkin(5) * t209 + t226 + (t301 * t59 + t219) * t178 + t225;
t5 = t178 * t226 + t179 * t224;
t3 = -t271 / 0.2e1 + t121 * t301 + (t287 / 0.2e1 + t103 / 0.2e1) * t180 + t181;
t2 = t226 * t265 - t312 + (t282 / 0.2e1 + t310) * pkin(5);
t25 = [0, 0, 0, 0, qJD(2), t169, qJD(2), qJD(3), qJD(3) * t176 + t169, -t211, t135 * qJD(4), 0, 0, 0, t176 * t232 + t249, qJD(3) * t180 - t176 * t234, -t173 * t211 - t174 * t216, -qJD(5) * t110 + t178 * t198, -qJD(4) * t115 - t178 * t218, qJD(4) * t113 - t178 * t217, t211, -qJD(2) * t112 + qJD(4) * t22 + qJD(5) * t46 + t179 * t249, -qJD(2) * t114 - qJD(4) * t23 - qJD(5) * t45 - t177 * t249, -qJD(3) * t111 - qJD(4) * t11 - qJD(5) * t10 + qJD(6) * t109, qJD(2) * t19 + qJD(3) * t21 + qJD(4) * t9 + qJD(5) * t14 - qJD(6) * t20; 0, 0, 0, 0, qJD(1), t170, qJD(1), 0, t170, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t241, -t239, 0, qJD(3) * t69 + qJD(4) * t18 + qJD(5) * t13 + qJD(6) * t80 + t263; 0, 0, 0, 0, 0, 0, 0, qJD(1), t235, 0, 0, 0, 0, 0, t163, t164, 0, 0, 0, 0, 0, qJD(5) * t83 + t150, qJD(5) * t81 - t147, -t242, qJD(2) * t69 + qJD(4) * t3 + qJD(5) * t7 + t261; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t210, t236, -t234, -t232, 0, t164 * t176 - t175 * t234, -t163 * t176 - t175 * t232, -t92 + (-t164 * t173 - t215) * t178 (t155 - t156) * qJD(4) + (-qJD(5) + t163) * t203, t149 - t238, t213 + t240, t105, t260 + (t177 * t199 - t222) * qJD(4) + t43 * qJD(5), -t259 + (-pkin(8) * t266 + (pkin(4) * t179 + t272) * t178) * qJD(4) + t42 * qJD(5), -t281 + ((t62 + t278) * t179 + (-t57 - t277) * t177) * qJD(4) + t5 * qJD(5), t284 + t18 * qJD(2) + t3 * qJD(3) + (-t103 * t157 + t121 * t57 - t122 * t62) * qJD(4) + t2 * qJD(5) + t15 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t185, -t184, t177 * t194, t188, t154, qJD(3) * t83 + qJD(4) * t43 - qJD(5) * t71 + t257, qJD(3) * t81 + qJD(4) * t42 + qJD(5) * t70 - t258, pkin(5) * t218 + qJD(4) * t5 - t283, qJD(2) * t13 + qJD(3) * t7 + qJD(4) * t2 - t292 * t60 + t275; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t243, qJD(2) * t80 + qJD(4) * t15 - t262; 0, 0, 0, 0, -qJD(1), -t170, -qJD(1), 0, -qJD(3) - t170, 0, 0, 0, 0, 0, -t232, t234, 0, 0, 0, 0, 0, t197 + t241, t239 + t252, t108 * qJD(4), -qJD(3) * t68 - qJD(4) * t17 - qJD(5) * t12 - qJD(6) * t79 - t263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -qJD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t256; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t164, t163, 0, 0, 0, 0, 0, -t212, t214, t244, -t274; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, t107, 0, t159 - t280; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t253; 0, 0, 0, 0, 0, 0, 0, -qJD(1), -t202, 0, 0, 0, 0, 0, -t163, -t164, 0, 0, 0, 0, 0, -qJD(5) * t78 - t150, qJD(5) * t76 + t147, t242, qJD(2) * t68 + qJD(4) * t4 + qJD(5) * t8 - t261; 0, 0, 0, 0, 0, 0, 0, 0, qJD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t256; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t234, -t232, 0, 0, 0, 0, 0, -t179 * t234 - t218, t148 - t217, t111 * qJD(4) (t178 * t157 + t180 * t67) * qJD(4) - t36 * qJD(5) + t82 * qJD(6) + t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t252 - t254, t197 + t255, 0, -pkin(5) * t152 - qJD(4) * t36 + t285; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, t210, -t236, 0, 0, 0, -t202 * t180, t202 * t178, t173 * t210 - t92, t188 * t313, t152 + t238, -t151 - t240, -t105, qJD(5) * t64 + t179 * t233 - t260, qJD(5) * t63 - t177 * t233 + t259, -qJD(2) * t108 + qJD(5) * t6 + t281, qJD(2) * t17 - qJD(3) * t4 - qJD(5) * t1 + qJD(6) * t16 - t284; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t164, -t163, 0, 0, 0, 0, 0, t212, -t214, -t244, t274; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(6) * t77 - t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t216, t136 * qJD(5), 0, 0, 0, -pkin(4) * t246, -pkin(4) * t245, t134 * qJD(6), qJD(5) * t24 + qJD(6) * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, -t183, t107, -t106, -t164 / 0.2e1, -pkin(8) * t245 - t192, pkin(8) * t246 - t191, -pkin(5) * t245 + t289, t122 * t292 + t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t237, t182; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t185, t184 (t214 - t247) * t178, t179 * t210 + t148, t154, qJD(3) * t78 - qJD(4) * t64 - t177 * t250 - t257, -qJD(3) * t76 - qJD(4) * t63 - t179 * t250 + t258, -qJD(4) * t6 + t283, qJD(2) * t12 - qJD(3) * t8 + qJD(4) * t1 - qJD(6) * t230 - t275; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t147, -t150, 0, t280; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t254, -t255, 0, -t285; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t189, t183, -t150, t147, t164 / 0.2e1, t192, t191, -t289, -qJD(6) * t298 - t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t243, pkin(5) * t217 + qJD(2) * t79 - qJD(4) * t16 + t262; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t253; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t237, t159 - t182; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t25;
