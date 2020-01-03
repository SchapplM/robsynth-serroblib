% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% tauc_reg [5x31]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRR11_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR11_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:43:42
% EndTime: 2019-12-31 22:44:01
% DurationCPUTime: 6.79s
% Computational Cost: add. (6444->460), mult. (16996->669), div. (0->0), fcn. (13472->10), ass. (0->218)
t202 = sin(pkin(5));
t211 = cos(qJ(2));
t279 = qJD(1) * t211;
t257 = t202 * t279;
t337 = qJD(3) - t257;
t207 = sin(qJ(2));
t203 = cos(pkin(5));
t269 = t203 * qJD(1);
t266 = pkin(1) * t269;
t155 = pkin(7) * t257 + t207 * t266;
t206 = sin(qJ(3));
t210 = cos(qJ(3));
t336 = -t155 + t337 * (pkin(3) * t206 - pkin(9) * t210);
t205 = sin(qJ(4));
t209 = cos(qJ(4));
t280 = qJD(1) * t202;
t263 = t207 * t280;
t125 = t205 * t210 * t257 - t209 * t263;
t274 = qJD(3) * t210;
t331 = -t205 * t274 + t125;
t287 = t210 * t211;
t126 = (t205 * t207 + t209 * t287) * t280;
t230 = t209 * t274 - t126;
t273 = qJD(4) * t205;
t335 = -t206 * t273 + t230;
t243 = qJD(2) + t269;
t120 = t243 * pkin(8) + t155;
t151 = (-pkin(2) * t211 - pkin(8) * t207 - pkin(1)) * t202;
t131 = qJD(1) * t151;
t74 = -t206 * t120 + t210 * t131;
t63 = -pkin(3) * t337 - t74;
t139 = t206 * t243 + t210 * t263;
t98 = t205 * t139 - t209 * t337;
t33 = t98 * pkin(4) + t63;
t100 = t209 * t139 + t205 * t337;
t204 = sin(qJ(5));
t208 = cos(qJ(5));
t46 = t204 * t100 + t208 * t98;
t334 = t33 * t46;
t227 = t208 * t100 - t204 * t98;
t333 = t227 * t46;
t272 = qJD(4) * t209;
t332 = t206 * t272 - t331;
t322 = qJD(4) + qJD(5);
t329 = -t206 * t263 + t210 * t243;
t330 = -t329 + t322;
t328 = t227 ^ 2 - t46 ^ 2;
t270 = qJD(5) * t208;
t271 = qJD(5) * t204;
t268 = qJD(1) * qJD(2);
t252 = t202 * t268;
t238 = t211 * t252;
t101 = t329 * qJD(3) + t210 * t238;
t239 = t207 * t252;
t43 = t209 * t101 - t139 * t273 + t205 * t239 + t272 * t337;
t44 = qJD(4) * t100 + t205 * t101 - t209 * t239;
t12 = -t100 * t271 - t204 * t44 + t208 * t43 - t98 * t270;
t132 = qJD(4) - t329;
t121 = qJD(5) + t132;
t327 = t46 * t121 + t12;
t276 = qJD(3) * t206;
t315 = pkin(8) * t205;
t325 = t336 * t209 + t276 * t315;
t170 = t204 * t209 + t208 * t205;
t158 = t170 * t206;
t324 = t207 * t211;
t183 = -t210 * pkin(3) - t206 * pkin(9) - pkin(2);
t152 = -pkin(7) * t263 + t211 * t266;
t224 = (pkin(2) * t207 - pkin(8) * t211) * t202;
t153 = qJD(1) * t224;
t285 = t210 * t152 + t206 * t153;
t87 = pkin(9) * t263 + t285;
t323 = -t183 * t272 - t336 * t205 + t209 * t87;
t102 = qJD(3) * t139 + t206 * t238;
t119 = -t243 * pkin(2) - t152;
t61 = -pkin(3) * t329 - t139 * pkin(9) + t119;
t75 = t210 * t120 + t206 * t131;
t64 = pkin(9) * t337 + t75;
t25 = t205 * t61 + t209 * t64;
t154 = qJD(2) * t224;
t146 = qJD(1) * t154;
t294 = t202 * t207;
t193 = pkin(7) * t294;
t317 = pkin(1) * t211;
t156 = (t203 * t317 - t193) * qJD(2);
t147 = qJD(1) * t156;
t218 = -t120 * t276 + t131 * t274 + t206 * t146 + t210 * t147;
t31 = pkin(9) * t239 + t218;
t293 = t202 * t211;
t267 = pkin(7) * t293;
t318 = pkin(1) * t207;
t157 = (t203 * t318 + t267) * qJD(2);
t148 = qJD(1) * t157;
t39 = t102 * pkin(3) - t101 * pkin(9) + t148;
t11 = -t25 * qJD(4) - t205 * t31 + t209 * t39;
t4 = t102 * pkin(4) - t43 * pkin(10) + t11;
t10 = t205 * t39 + t209 * t31 + t61 * t272 - t64 * t273;
t5 = -t44 * pkin(10) + t10;
t24 = -t205 * t64 + t209 * t61;
t19 = -t100 * pkin(10) + t24;
t16 = t132 * pkin(4) + t19;
t20 = -t98 * pkin(10) + t25;
t309 = t208 * t20;
t8 = t204 * t16 + t309;
t2 = -t8 * qJD(5) - t204 * t5 + t208 * t4;
t321 = -t33 * t227 + t2;
t13 = t227 * qJD(5) + t204 * t43 + t208 * t44;
t320 = t121 * t227 - t13;
t18 = t20 * t271;
t251 = qJD(5) * t16 + t5;
t1 = t204 * t4 + t251 * t208 - t18;
t212 = qJD(1) ^ 2;
t319 = pkin(9) + pkin(10);
t316 = pkin(4) * t206;
t90 = pkin(3) * t139 - pkin(9) * t329;
t314 = t205 * t90 + t209 * t74;
t169 = t204 * t205 - t208 * t209;
t313 = t204 * t125 - t208 * t126 - t322 * t158 - t169 * t274;
t290 = t206 * t209;
t291 = t205 * t206;
t312 = -t271 * t291 + (t322 * t290 - t331) * t208 + t335 * t204;
t149 = t193 + (-pkin(2) - t317) * t203;
t160 = -t203 * t210 + t206 * t294;
t161 = t203 * t206 + t210 * t294;
t82 = t160 * pkin(3) - t161 * pkin(9) + t149;
t150 = t267 + (pkin(8) + t318) * t203;
t286 = t210 * t150 + t206 * t151;
t84 = -pkin(9) * t293 + t286;
t311 = t205 * t82 + t209 * t84;
t241 = t120 * t274 + t131 * t276 - t210 * t146 + t206 * t147;
t32 = -pkin(3) * t239 + t241;
t308 = t32 * t205;
t307 = t32 * t209;
t306 = t43 * t205;
t305 = t98 * t132;
t304 = t330 * t169;
t303 = t330 * t170;
t246 = -t206 * t152 + t210 * t153;
t86 = -pkin(3) * t263 - t246;
t302 = t332 * pkin(4) + pkin(8) * t274 - t86;
t301 = t100 * t132;
t300 = t102 * t210;
t299 = t329 * t337;
t298 = t329 * t205;
t297 = t139 * t337;
t221 = t337 * t206;
t296 = t337 * t210;
t199 = t202 ^ 2;
t295 = t199 * t212;
t292 = t205 * t102;
t289 = t209 * t102;
t288 = t209 * t210;
t195 = pkin(8) * t288;
t282 = t205 * t183 + t195;
t281 = t207 ^ 2 - t211 ^ 2;
t278 = qJD(2) * t207;
t277 = qJD(2) * t210;
t275 = qJD(3) * t209;
t265 = t205 * t293;
t264 = qJD(4) * t319;
t262 = t202 * t278;
t261 = qJD(2) * t293;
t260 = t132 * t273;
t258 = t202 * t203 * t212;
t253 = t199 * t268;
t249 = -t205 * t84 + t209 * t82;
t247 = -t206 * t150 + t210 * t151;
t245 = t132 * t209;
t244 = 0.2e1 * t253;
t242 = qJD(2) + 0.2e1 * t269;
t237 = -t75 + (t273 - t298) * pkin(4);
t167 = t209 * t183;
t104 = -pkin(10) * t290 + t167 + (-pkin(4) - t315) * t210;
t236 = -qJD(5) * t104 - (-t206 * t275 - t210 * t273) * pkin(8) + t323 + t332 * pkin(10);
t116 = -pkin(10) * t291 + t282;
t235 = -t126 * pkin(10) + qJD(5) * t116 - t205 * t87 + t257 * t316 - (-pkin(10) * t288 + t316) * qJD(3) - (-t195 + (pkin(10) * t206 - t183) * t205) * qJD(4) - t325;
t234 = -0.2e1 * pkin(1) * t253;
t187 = t319 * t205;
t232 = -pkin(10) * t298 + qJD(5) * t187 + t205 * t264 + t314;
t188 = t319 * t209;
t89 = t209 * t90;
t231 = t139 * pkin(4) + qJD(5) * t188 - t205 * t74 + t89 + (-pkin(10) * t329 + t264) * t209;
t83 = pkin(3) * t293 - t247;
t113 = t161 * t209 - t265;
t22 = t160 * pkin(4) - t113 * pkin(10) + t249;
t112 = t161 * t205 + t209 * t293;
t26 = -pkin(10) * t112 + t311;
t229 = -t204 * t26 + t208 * t22;
t228 = t204 * t22 + t208 * t26;
t65 = t208 * t112 + t204 * t113;
t66 = -t204 * t112 + t208 * t113;
t225 = -t150 * t274 - t151 * t276 + t210 * t154 - t206 * t156;
t217 = -t150 * t276 + t151 * t274 + t206 * t154 + t210 * t156;
t35 = pkin(9) * t262 + t217;
t110 = qJD(3) * t161 + t206 * t261;
t111 = -qJD(3) * t160 + t210 * t261;
t54 = t110 * pkin(3) - t111 * pkin(9) + t157;
t223 = t205 * t54 + t209 * t35 + t82 * t272 - t84 * t273;
t222 = -t132 * t272 - t292;
t219 = -pkin(9) * t102 + t132 * t63;
t215 = pkin(1) * (-t203 * t268 + t295);
t36 = -pkin(3) * t262 - t225;
t214 = -t311 * qJD(4) - t205 * t35 + t209 * t54;
t198 = -t209 * pkin(4) - pkin(3);
t175 = (pkin(4) * t205 + pkin(8)) * t206;
t159 = t169 * t206;
t85 = t102 * t160;
t60 = -qJD(4) * t112 + t111 * t209 + t205 * t262;
t59 = -qJD(4) * t265 + t111 * t205 + t161 * t272 - t209 * t262;
t50 = t112 * pkin(4) + t83;
t21 = t59 * pkin(4) + t36;
t17 = t44 * pkin(4) + t32;
t15 = qJD(5) * t66 + t204 * t60 + t208 * t59;
t14 = -qJD(5) * t65 - t204 * t59 + t208 * t60;
t9 = -t59 * pkin(10) + t223;
t7 = t208 * t16 - t204 * t20;
t6 = t110 * pkin(4) - t60 * pkin(10) + t214;
t3 = [0, 0, 0, t244 * t324, -t281 * t244, t242 * t261, -t242 * t262, 0, -t148 * t203 - t157 * t243 + t207 * t234, -t147 * t203 - t156 * t243 + t211 * t234, t101 * t161 + t139 * t111, -t101 * t160 - t161 * t102 - t139 * t110 + t111 * t329, t111 * t337 + (-t101 * t211 + (qJD(1) * t161 + t139) * t278) * t202, -t110 * t337 + (t102 * t211 + (-qJD(1) * t160 + t329) * t278) * t202, (-t199 * t279 + t202 * t337) * t278, t225 * t337 - t157 * t329 + t149 * t102 + t148 * t160 + t119 * t110 + (t241 * t211 + (t247 * qJD(1) + t74) * t278) * t202, -t217 * t337 + t157 * t139 + t149 * t101 + t148 * t161 + t119 * t111 + (t218 * t211 + (-t286 * qJD(1) - t75) * t278) * t202, t100 * t60 + t113 * t43, -t100 * t59 - t112 * t43 - t113 * t44 - t60 * t98, t100 * t110 + t113 * t102 + t60 * t132 + t43 * t160, -t112 * t102 - t98 * t110 - t59 * t132 - t44 * t160, t110 * t132 + t85, t102 * t249 + t11 * t160 + t24 * t110 + t32 * t112 + t132 * t214 + t36 * t98 + t83 * t44 + t63 * t59, -t10 * t160 + t36 * t100 - t311 * t102 - t25 * t110 + t32 * t113 - t223 * t132 + t83 * t43 + t63 * t60, t12 * t66 + t14 * t227, -t12 * t65 - t13 * t66 - t14 * t46 - t15 * t227, t66 * t102 + t110 * t227 + t12 * t160 + t14 * t121, -t65 * t102 - t46 * t110 - t15 * t121 - t13 * t160, t110 * t121 + t85, t21 * t46 + t50 * t13 + t17 * t65 + t33 * t15 + (-qJD(5) * t228 - t204 * t9 + t208 * t6) * t121 + t229 * t102 + t2 * t160 + t7 * t110, t21 * t227 + t50 * t12 + t17 * t66 + t33 * t14 - (qJD(5) * t229 + t204 * t6 + t208 * t9) * t121 - t228 * t102 - t1 * t160 - t8 * t110; 0, 0, 0, -t295 * t324, t281 * t295, -t211 * t258, t207 * t258, 0, -pkin(7) * t238 + t155 * t243 + t207 * t215, pkin(7) * t239 + t152 * t243 + t211 * t215, t101 * t206 + t139 * t296, (t101 + t299) * t210 + (-t102 - t297) * t206, t337 * t274 + (-t337 * t287 + (t206 * qJD(2) - t139) * t207) * t280, -t337 * t276 + (t211 * t221 + (-t329 + t277) * t207) * t280, -t337 * t263, -pkin(2) * t102 - t148 * t210 - t246 * t337 + t155 * t329 + (-pkin(8) * t296 + t119 * t206) * qJD(3) + (-t74 * t207 + (-pkin(8) * t278 - t119 * t211) * t206) * t280, -pkin(2) * t101 + t148 * t206 + t285 * t337 - t155 * t139 + (pkin(8) * t221 + t119 * t210) * qJD(3) + (-t119 * t287 + (-pkin(8) * t277 + t75) * t207) * t280, t335 * t100 + t43 * t290, t100 * t125 + t126 * t98 + (-t100 * t205 - t209 * t98) * t274 + (-t306 - t209 * t44 + (-t100 * t209 + t205 * t98) * qJD(4)) * t206, -t43 * t210 + t230 * t132 + (t100 * t337 - t260 + t289) * t206, t44 * t210 + t331 * t132 + (-t337 * t98 + t222) * t206, t132 * t221 - t300, t167 * t102 - t63 * t125 - t86 * t98 + ((-qJD(4) * t183 + t87) * t205 + t325) * t132 + (t63 * t205 * qJD(3) - t11 + (qJD(3) * t98 + t222) * pkin(8)) * t210 + (pkin(8) * t44 + t24 * t337 + t272 * t63 + t308) * t206, -t282 * t102 - t86 * t100 - t63 * t126 + t323 * t132 + (t63 * t275 + t10 + (qJD(3) * t100 + t260) * pkin(8)) * t210 + (-t63 * t273 + t307 - t337 * t25 + (t132 * t275 + t43) * pkin(8)) * t206, -t12 * t159 + t227 * t313, -t12 * t158 + t159 * t13 - t227 * t312 - t313 * t46, -t159 * t102 - t12 * t210 + t313 * t121 + t221 * t227, -t158 * t102 - t312 * t121 + t13 * t210 - t46 * t221, t121 * t221 - t300, t175 * t13 + t17 * t158 + (t208 * t104 - t204 * t116) * t102 - t2 * t210 + t302 * t46 + t312 * t33 + t7 * t221 + (t204 * t236 - t208 * t235) * t121, t175 * t12 - t17 * t159 - (t204 * t104 + t208 * t116) * t102 + t1 * t210 + t302 * t227 + t313 * t33 - t8 * t221 + (t204 * t235 + t208 * t236) * t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t139 * t329, t139 ^ 2 - t329 ^ 2, t101 - t299, -t102 + t297, t239, -t119 * t139 + t337 * t75 - t241, -t119 * t329 + t337 * t74 - t218, t100 * t245 + t306, (t43 - t305) * t209 + (-t44 - t301) * t205, -t100 * t139 + t132 * t245 + t292, -t132 ^ 2 * t205 + t98 * t139 + t289, -t132 * t139, -pkin(3) * t44 - t24 * t139 - t307 - t75 * t98 + (-pkin(9) * t272 - t89) * t132 + (t74 * t132 + t219) * t205, -pkin(3) * t43 - t75 * t100 + t25 * t139 + t308 + (pkin(9) * t273 + t314) * t132 + t219 * t209, t12 * t170 - t227 * t304, -t12 * t169 - t170 * t13 - t227 * t303 + t304 * t46, t170 * t102 - t304 * t121 - t139 * t227, -t169 * t102 - t303 * t121 + t46 * t139, -t121 * t139, t198 * t13 + t17 * t169 + (-t208 * t187 - t204 * t188) * t102 - t7 * t139 + t237 * t46 + t303 * t33 + (t204 * t232 - t208 * t231) * t121, t198 * t12 + t17 * t170 - (-t204 * t187 + t208 * t188) * t102 + t8 * t139 + t237 * t227 - t304 * t33 + (t204 * t231 + t208 * t232) * t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100 * t98, t100 ^ 2 - t98 ^ 2, t43 + t305, t301 - t44, t102, -t63 * t100 + t25 * t132 + t11, t24 * t132 + t63 * t98 - t10, t333, t328, t327, t320, t102, -(-t204 * t19 - t309) * t121 + (-t100 * t46 + t208 * t102 - t121 * t271) * pkin(4) + t321, t334 + t18 + (-t20 * t121 - t4) * t204 + (t19 * t121 - t251) * t208 + (-t100 * t227 - t204 * t102 - t121 * t270) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t333, t328, t327, t320, t102, t8 * t121 + t321, t7 * t121 - t1 + t334;];
tauc_reg = t3;
