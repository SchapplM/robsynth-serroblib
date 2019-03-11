% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRRP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:55:39
% EndTime: 2019-03-09 11:55:51
% DurationCPUTime: 3.74s
% Computational Cost: add. (9171->423), mult. (22837->547), div. (0->0), fcn. (16992->8), ass. (0->214)
t191 = sin(qJ(4));
t189 = sin(pkin(10));
t181 = pkin(2) * t189 + pkin(8);
t287 = pkin(9) + t181;
t227 = qJD(4) * t287;
t194 = cos(qJ(2));
t270 = cos(pkin(10));
t225 = t270 * t194;
t176 = qJD(1) * t225;
t192 = sin(qJ(2));
t249 = qJD(1) * t192;
t152 = -t189 * t249 + t176;
t264 = t152 * t191;
t286 = -qJ(3) - pkin(7);
t173 = t286 * t194;
t169 = qJD(1) * t173;
t157 = t189 * t169;
t172 = t286 * t192;
t168 = qJD(1) * t172;
t113 = t168 * t270 + t157;
t193 = cos(qJ(4));
t164 = t189 * t194 + t192 * t270;
t154 = t164 * qJD(1);
t94 = pkin(2) * t249 + pkin(3) * t154 - pkin(8) * t152;
t274 = t113 * t193 + t191 * t94;
t312 = -pkin(9) * t264 + t191 * t227 + t274;
t263 = t152 * t193;
t87 = t193 * t94;
t311 = pkin(4) * t154 - pkin(9) * t263 - t113 * t191 + t193 * t227 + t87;
t146 = qJD(4) - t152;
t142 = qJD(5) + t146;
t153 = t164 * qJD(2);
t143 = qJD(1) * t153;
t190 = sin(qJ(5));
t289 = cos(qJ(5));
t239 = t289 * t191;
t167 = t190 * t193 + t239;
t258 = t190 * t191;
t205 = t193 * t289 - t258;
t296 = qJD(4) + qJD(5);
t231 = t289 * qJD(5);
t298 = qJD(4) * t289 + t231;
t273 = t152 * t205 - t193 * t298 + t258 * t296;
t218 = -t142 * t273 + t167 * t143;
t124 = qJD(2) * t193 - t154 * t191;
t125 = qJD(2) * t191 + t154 * t193;
t207 = t124 * t190 + t125 * t289;
t275 = t207 * t154;
t310 = t218 + t275;
t118 = t296 * t167;
t272 = -t152 * t167 + t118;
t219 = -t142 * t272 + t143 * t205;
t64 = -t124 * t289 + t125 * t190;
t279 = t154 * t64;
t309 = t219 - t279;
t281 = qJD(2) * pkin(2);
t160 = t168 + t281;
t109 = t160 * t270 + t157;
t102 = -qJD(2) * pkin(3) - t109;
t62 = -pkin(4) * t124 + t102;
t28 = pkin(5) * t64 - qJ(6) * t207 + t62;
t308 = t28 * t64;
t307 = t62 * t64;
t288 = t207 * t64;
t201 = -t189 * t192 + t225;
t156 = t201 * qJD(2);
t247 = qJD(4) * t193;
t306 = t156 * t191 + t164 * t247;
t290 = t207 ^ 2;
t305 = -t64 ^ 2 + t290;
t245 = qJD(1) * qJD(2);
t230 = t192 * t245;
t174 = t189 * t230;
t199 = qJD(2) * t176 - t174;
t198 = qJD(4) * t124 + t193 * t199;
t248 = qJD(4) * t191;
t241 = qJD(2) * t248 + t154 * t247 + t191 * t199;
t246 = qJD(5) * t190;
t30 = -t124 * t231 + t125 * t246 + t190 * t241 - t198 * t289;
t21 = t142 * t64 - t30;
t46 = pkin(5) * t207 + qJ(6) * t64;
t304 = -0.2e1 * t245;
t161 = t287 * t191;
t162 = t287 * t193;
t206 = -t161 * t289 - t162 * t190;
t302 = -t206 * qJD(5) + t190 * t311 + t289 * t312;
t106 = -t161 * t190 + t162 * t289;
t301 = -t106 * qJD(5) + t190 * t312 - t289 * t311;
t226 = t270 * t169;
t112 = t168 * t189 - t226;
t221 = -t112 + (t248 - t264) * pkin(4);
t257 = t191 * t143;
t128 = t193 * t143;
t300 = t146 * t248 - t128;
t299 = -t125 * t248 + t193 * t198;
t236 = t164 * t248;
t261 = t156 * t193;
t297 = -t236 + t261;
t138 = t143 * pkin(5);
t110 = t160 * t189 - t226;
t103 = qJD(2) * pkin(8) + t110;
t179 = pkin(2) * t230;
t82 = pkin(3) * t143 - pkin(8) * t199 + t179;
t75 = t193 * t82;
t217 = -t103 * t247 + t75;
t228 = qJD(2) * t286;
t150 = t194 * qJD(3) + t192 * t228;
t130 = t150 * qJD(1);
t151 = -t192 * qJD(3) + t194 * t228;
t131 = t151 * qJD(1);
t78 = t130 * t270 + t131 * t189;
t276 = t191 * t78;
t240 = -pkin(2) * t194 - pkin(1);
t216 = t240 * qJD(1);
t170 = qJD(3) + t216;
t83 = -t152 * pkin(3) - t154 * pkin(8) + t170;
t11 = t143 * pkin(4) - pkin(9) * t198 - t248 * t83 + t217 - t276;
t203 = -t103 * t248 + t191 * t82 + t193 * t78 + t247 * t83;
t16 = -pkin(9) * t241 + t203;
t53 = -t103 * t191 + t193 * t83;
t44 = -pkin(9) * t125 + t53;
t36 = pkin(4) * t146 + t44;
t54 = t193 * t103 + t191 * t83;
t45 = pkin(9) * t124 + t54;
t224 = -t11 * t289 + t16 * t190 + t231 * t45 + t246 * t36;
t2 = -t138 + t224;
t200 = t207 * t28 + t2;
t31 = t124 * t246 + t125 * t231 + t190 * t198 + t241 * t289;
t295 = t142 * t207 - t31;
t294 = -t207 * t62 - t224;
t293 = -t146 ^ 2 * t193 - t257;
t292 = -t205 * t30 - t207 * t272;
t108 = -pkin(3) * t201 - pkin(8) * t164 + t240;
t101 = t193 * t108;
t123 = t172 * t189 - t173 * t270;
t259 = t164 * t193;
t50 = -pkin(4) * t201 - pkin(9) * t259 - t123 * t191 + t101;
t116 = t193 * t123;
t251 = t108 * t191 + t116;
t260 = t164 * t191;
t57 = -pkin(9) * t260 + t251;
t209 = t190 * t50 + t289 * t57;
t244 = t192 * t281;
t95 = pkin(3) * t153 - pkin(8) * t156 + t244;
t88 = t193 * t95;
t93 = t150 * t270 + t151 * t189;
t23 = -pkin(9) * t261 + t153 * pkin(4) - t191 * t93 + t88 + (-t116 + (pkin(9) * t164 - t108) * t191) * qJD(4);
t202 = t108 * t247 - t123 * t248 + t191 * t95 + t193 * t93;
t27 = -pkin(9) * t306 + t202;
t291 = -qJD(5) * t209 - t190 * t27 + t23 * t289;
t285 = qJ(6) * t154 + t302;
t284 = -pkin(5) * t154 + t301;
t283 = -pkin(5) * t272 - qJ(6) * t273 + qJD(6) * t167 - t221;
t242 = t289 * t45;
t13 = t190 * t36 + t242;
t280 = t13 * t142;
t277 = t190 * t45;
t20 = t289 * t44 - t277;
t271 = pkin(4) * t231 + qJD(6) - t20;
t269 = t206 * t143;
t268 = t106 * t143;
t267 = t124 * t154;
t266 = t125 * t146;
t265 = t125 * t154;
t256 = t193 * t124;
t196 = qJD(1) ^ 2;
t255 = t194 * t196;
t195 = qJD(2) ^ 2;
t254 = t195 * t192;
t253 = t195 * t194;
t12 = t289 * t36 - t277;
t252 = qJD(6) - t12;
t250 = t192 ^ 2 - t194 ^ 2;
t234 = t102 * t247;
t229 = -t11 * t190 - t16 * t289 - t231 * t36 + t246 * t45;
t223 = pkin(1) * t304;
t77 = t130 * t189 - t131 * t270;
t92 = t150 * t189 - t151 * t270;
t122 = -t172 * t270 - t173 * t189;
t222 = -t167 * t31 + t273 * t64;
t19 = t190 * t44 + t242;
t220 = pkin(4) * t246 - t19;
t183 = -pkin(2) * t270 - pkin(3);
t91 = pkin(4) * t260 + t122;
t215 = -t123 * t143 + t77 * t164;
t132 = t142 * qJD(6);
t134 = t143 * qJ(6);
t1 = t134 + t132 - t229;
t59 = pkin(4) * t306 + t92;
t214 = t146 * t264 - t300;
t213 = t12 * t142 + t229;
t211 = -t190 * t57 + t289 * t50;
t208 = t190 * t23 + t231 * t50 - t246 * t57 + t27 * t289;
t171 = -pkin(4) * t193 + t183;
t51 = pkin(4) * t241 + t77;
t186 = -pkin(4) * t289 - pkin(5);
t182 = pkin(4) * t190 + qJ(6);
t111 = t143 * t201;
t99 = -pkin(5) * t205 - qJ(6) * t167 + t171;
t98 = t205 * t164;
t97 = t167 * t164;
t39 = pkin(5) * t97 - qJ(6) * t98 + t91;
t38 = t156 * t239 - t190 * t236 - t246 * t260 + (t156 * t190 + t164 * t298) * t193;
t37 = t118 * t164 - t156 * t205;
t33 = pkin(4) * t125 + t46;
t25 = pkin(5) * t201 - t211;
t24 = -qJ(6) * t201 + t209;
t8 = qJ(6) * t142 + t13;
t7 = -pkin(5) * t142 + t252;
t6 = pkin(5) * t38 + qJ(6) * t37 - qJD(6) * t98 + t59;
t5 = pkin(5) * t31 + qJ(6) * t30 - qJD(6) * t207 + t51;
t4 = -t153 * pkin(5) - t291;
t3 = qJ(6) * t153 - qJD(6) * t201 + t208;
t9 = [0, 0, 0, 0.2e1 * t194 * t230, t250 * t304, t253, -t254, 0, -pkin(7) * t253 + t192 * t223, pkin(7) * t254 + t194 * t223, -t109 * t156 - t110 * t153 + t122 * t199 + t152 * t93 + t154 * t92 + t201 * t78 + t215, -t109 * t92 + t110 * t93 + t77 * t122 + t78 * t123 + (t170 + t216) * t244, t125 * t261 + t164 * t299, -t124 * t236 - t125 * t306 + t156 * t256 - t198 * t260 - t241 * t259, t125 * t153 + t128 * t164 + t146 * t297 - t198 * t201, t124 * t153 - t146 * t306 - t164 * t257 + t201 * t241, t146 * t153 - t111 (-t123 * t247 + t88) * t146 + t101 * t143 - t217 * t201 + t53 * t153 - t92 * t124 + t122 * t241 + t164 * t234 + ((-qJD(4) * t108 - t93) * t146 - (-qJD(4) * t83 - t78) * t201 + t102 * t156 + t215) * t191, t102 * t297 + t122 * t198 + t92 * t125 - t143 * t251 - t146 * t202 - t54 * t153 + t201 * t203 + t259 * t77, -t207 * t37 - t30 * t98, -t207 * t38 + t30 * t97 - t31 * t98 + t37 * t64, -t142 * t37 + t143 * t98 + t153 * t207 + t201 * t30, -t142 * t38 - t143 * t97 - t153 * t64 + t201 * t31, t142 * t153 - t111, t12 * t153 + t142 * t291 + t143 * t211 + t201 * t224 + t91 * t31 + t62 * t38 + t51 * t97 + t59 * t64, -t13 * t153 - t142 * t208 - t143 * t209 - t201 * t229 + t207 * t59 - t30 * t91 - t37 * t62 + t51 * t98, -t142 * t4 - t143 * t25 - t153 * t7 + t2 * t201 + t28 * t38 + t31 * t39 + t5 * t97 + t6 * t64, -t1 * t97 + t2 * t98 + t207 * t4 - t24 * t31 - t25 * t30 - t3 * t64 - t37 * t7 - t38 * t8, -t1 * t201 + t142 * t3 + t143 * t24 + t153 * t8 - t207 * t6 + t28 * t37 + t30 * t39 - t5 * t98, t1 * t24 + t2 * t25 + t28 * t6 + t3 * t8 + t39 * t5 + t4 * t7; 0, 0, 0, -t192 * t255, t250 * t196, 0, 0, 0, t196 * pkin(1) * t192, pkin(1) * t255 (t110 - t112) * t154 + (-t113 + t109) * t152 + (-t189 * t143 - t199 * t270) * pkin(2), t109 * t112 - t110 * t113 + (-t170 * t249 + t189 * t78 - t270 * t77) * pkin(2), -qJD(4) * t191 ^ 2 * t154 + ((-t174 + (t176 + qJD(4)) * qJD(2)) * t191 + t266) * t193, t124 * t247 + t125 * t264 - t152 * t256 - t191 * t241 + t299, -t265 - t293, t214 - t267, -t146 * t154, -t181 * t257 + t112 * t124 - t53 * t154 + t183 * t241 - t77 * t193 + (-t181 * t247 - t87 + (t102 + t113) * t191) * t146, -t102 * t263 - t112 * t125 + t274 * t146 + t54 * t154 + t181 * t300 + t183 * t198 + t77 * t191 + t234, -t30 * t167 - t207 * t273, t222 + t292, t218 - t275, t219 + t279, -t142 * t154, -t12 * t154 + t142 * t301 + t171 * t31 - t205 * t51 + t221 * t64 + t272 * t62 + t269, t13 * t154 + t142 * t302 + t51 * t167 - t171 * t30 + t207 * t221 - t273 * t62 - t268, t142 * t284 + t154 * t7 - t205 * t5 + t272 * t28 - t283 * t64 + t31 * t99 + t269, t1 * t205 - t106 * t31 + t167 * t2 + t206 * t30 - t207 * t284 - t272 * t8 - t273 * t7 + t285 * t64, -t142 * t285 - t154 * t8 - t167 * t5 + t207 * t283 + t273 * t28 + t30 * t99 + t268, t1 * t106 - t2 * t206 - t28 * t283 - t284 * t7 - t285 * t8 + t5 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152 ^ 2 - t154 ^ 2, t109 * t154 - t110 * t152 + t179, 0, 0, 0, 0, 0, t214 + t267, -t265 + t293, 0, 0, 0, 0, 0, t309, -t310, t309, t222 - t292, t310, t1 * t167 - t154 * t28 - t2 * t205 + t272 * t7 - t273 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t125 * t124, -t124 ^ 2 + t125 ^ 2, -t124 * t146 + t198, -t241 + t266, t143, -t102 * t125 - t276 + t75 + (-qJD(4) + t146) * t54, -t102 * t124 + t146 * t53 - t203, t288, t305, t21, t295, t143, t19 * t142 + (-t125 * t64 - t142 * t246 + t143 * t289) * pkin(4) + t294, t20 * t142 + t307 + (-t125 * t207 - t142 * t231 - t143 * t190) * pkin(4) + t229, -t142 * t220 - t143 * t186 - t33 * t64 - t200, -t182 * t31 - t186 * t30 + (t220 + t8) * t207 + (-t271 + t7) * t64, t142 * t271 + t143 * t182 + t207 * t33 + t1 - t308, t1 * t182 + t186 * t2 + t220 * t7 + t271 * t8 - t28 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t288, t305, t21, t295, t143, t280 + t294, t213 + t307, -t46 * t64 + t138 - t200 + t280, pkin(5) * t30 - t31 * qJ(6) + (-t13 + t8) * t207 + (t7 - t252) * t64, t207 * t46 + 0.2e1 * t132 + 0.2e1 * t134 - t213 - t308, -t2 * pkin(5) + t1 * qJ(6) - t7 * t13 + t252 * t8 - t28 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2) * t154 + t288, t21, -t142 ^ 2 - t290, -t142 * t8 + t200;];
tauc_reg  = t9;
