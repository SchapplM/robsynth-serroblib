% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
% 
% Output:
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRPPR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:23:06
% EndTime: 2019-03-09 15:23:15
% DurationCPUTime: 4.10s
% Computational Cost: add. (9562->397), mult. (25090->554), div. (0->0), fcn. (19073->10), ass. (0->240)
t220 = sin(qJ(3));
t221 = sin(qJ(2));
t286 = qJD(1) * t221;
t275 = t220 * t286;
t223 = cos(qJ(3));
t224 = cos(qJ(2));
t285 = qJD(1) * t224;
t276 = t223 * t285;
t167 = -t275 + t276;
t168 = -t220 * t285 - t223 * t286;
t217 = sin(pkin(10));
t309 = cos(pkin(10));
t139 = t309 * t167 + t168 * t217;
t133 = qJD(6) - t139;
t218 = cos(pkin(11));
t222 = cos(qJ(6));
t293 = t222 * t218;
t216 = sin(pkin(11));
t219 = sin(qJ(6));
t294 = t216 * t219;
t181 = -t293 + t294;
t310 = t133 * t181;
t182 = t216 * t222 + t218 * t219;
t332 = t133 * t182;
t213 = qJD(2) + qJD(3);
t237 = t217 * t167 - t309 * t168;
t122 = -t218 * t213 + t216 * t237;
t124 = t213 * t216 + t218 * t237;
t249 = t122 * t219 - t124 * t222;
t331 = t133 * t249;
t330 = t222 * t122;
t280 = qJD(1) * qJD(2);
t329 = -0.2e1 * t280;
t274 = t224 * t280;
t141 = qJD(3) * t276 - t213 * t275 + t223 * t274;
t322 = pkin(7) + pkin(8);
t193 = t322 * t224;
t187 = qJD(1) * t193;
t175 = t223 * t187;
t192 = t322 * t221;
t185 = qJD(1) * t192;
t319 = qJD(2) * pkin(2);
t177 = -t185 + t319;
t248 = -t220 * t177 - t175;
t277 = qJD(2) * t322;
t262 = qJD(1) * t277;
t178 = t221 * t262;
t179 = t224 * t262;
t266 = t220 * t178 - t223 * t179;
t231 = qJD(3) * t248 + t266;
t229 = -qJ(4) * t141 + qJD(4) * t168 + t231;
t184 = t220 * t224 + t221 * t223;
t325 = qJD(1) * t184;
t228 = t213 * t325;
t284 = qJD(3) * t220;
t265 = -t220 * t179 - t187 * t284;
t326 = (qJD(3) * t177 - t178) * t223;
t65 = -qJ(4) * t228 + t167 * qJD(4) + t265 + t326;
t33 = t217 * t229 + t309 * t65;
t30 = qJD(5) * t213 + t33;
t105 = t141 * t217 + t309 * t228;
t106 = t309 * t141 - t217 * t228;
t210 = pkin(2) * t286;
t227 = pkin(3) * t228 + qJD(2) * t210;
t39 = t105 * pkin(4) - t106 * qJ(5) - qJD(5) * t237 + t227;
t10 = -t216 * t30 + t218 * t39;
t11 = t216 * t39 + t218 * t30;
t258 = -t10 * t216 + t11 * t218;
t264 = t185 * t220 - t175;
t308 = qJ(4) * t167;
t126 = t264 - t308;
t162 = t168 * qJ(4);
t171 = t220 * t187;
t288 = -t223 * t185 - t171;
t127 = t162 + t288;
t269 = t309 * t220;
t312 = t309 * t126 - t127 * t217 + (t217 * t223 + t269) * qJD(3) * pkin(2);
t200 = t217 * t220 * pkin(2);
t283 = qJD(3) * t223;
t159 = t309 * pkin(2) * t283 - qJD(3) * t200;
t151 = qJD(5) + t159;
t88 = t217 * t126 + t309 * t127;
t321 = pkin(3) * t168;
t92 = pkin(4) * t237 - qJ(5) * t139 - t321;
t89 = t210 + t92;
t50 = -t216 * t88 + t218 * t89;
t328 = t151 * t216 + t50;
t51 = t216 * t89 + t218 * t88;
t327 = -t151 * t218 + t51;
t267 = t223 * t177 - t171;
t120 = t162 + t267;
t324 = -qJD(6) + t133;
t323 = -t105 * t182 + t310 * t133;
t26 = -qJD(6) * t249 + t106 * t182;
t320 = t218 * pkin(5);
t212 = t218 * pkin(9);
t235 = t184 * qJD(3);
t146 = qJD(2) * t184 + t235;
t183 = t220 * t221 - t223 * t224;
t186 = t221 * t277;
t188 = t224 * t277;
t236 = -t223 * t186 - t220 * t188 - t192 * t283 - t193 * t284;
t85 = -qJ(4) * t146 - qJD(4) * t183 + t236;
t145 = t213 * t183;
t247 = t192 * t220 - t193 * t223;
t230 = qJD(3) * t247 + t186 * t220 - t223 * t188;
t86 = qJ(4) * t145 - qJD(4) * t184 + t230;
t47 = t217 * t86 + t309 * t85;
t107 = -t145 * t217 + t309 * t146;
t108 = -t309 * t145 - t217 * t146;
t144 = -t217 * t183 + t309 * t184;
t211 = t221 * t319;
t272 = pkin(3) * t146 + t211;
t54 = pkin(4) * t107 - qJ(5) * t108 - qJD(5) * t144 + t272;
t14 = t216 * t54 + t218 * t47;
t111 = pkin(3) * t213 + t120;
t121 = -t248 + t308;
t270 = t309 * t121;
t71 = t217 * t111 + t270;
t69 = qJ(5) * t213 + t71;
t209 = -pkin(2) * t224 - pkin(1);
t191 = t209 * qJD(1);
t147 = -pkin(3) * t167 + qJD(4) + t191;
t78 = -pkin(4) * t139 - qJ(5) * t237 + t147;
t41 = t216 * t78 + t218 * t69;
t115 = t217 * t121;
t76 = t309 * t120 - t115;
t49 = t216 * t92 + t218 * t76;
t134 = -qJ(4) * t184 - t192 * t223 - t193 * t220;
t135 = -qJ(4) * t183 - t247;
t101 = t217 * t134 + t309 * t135;
t143 = t309 * t183 + t184 * t217;
t245 = pkin(3) * t183 + t209;
t99 = pkin(4) * t143 - qJ(5) * t144 + t245;
t58 = t218 * t101 + t216 * t99;
t80 = t124 * t219 + t330;
t317 = t237 * t80;
t316 = t237 * t249;
t70 = t309 * t111 - t115;
t68 = -t213 * pkin(4) + qJD(5) - t70;
t315 = t139 * t68;
t100 = -t309 * t134 + t135 * t217;
t32 = t217 * t65 - t309 * t229;
t314 = t32 * t100;
t302 = t139 * t216;
t130 = pkin(5) * t302;
t313 = -t130 + t312;
t311 = t159 - t88;
t306 = t106 * t216;
t305 = t106 * t218;
t304 = t108 * t216;
t303 = t133 * t237;
t301 = t139 * t218;
t300 = t144 * t216;
t299 = t144 * t218;
t296 = t168 * t167;
t295 = t191 * t168;
t226 = qJD(1) ^ 2;
t292 = t224 * t226;
t225 = qJD(2) ^ 2;
t291 = t225 * t221;
t290 = t225 * t224;
t281 = qJD(6) * t222;
t289 = t106 * t293 - t122 * t281;
t208 = pkin(2) * t223 + pkin(3);
t161 = pkin(2) * t269 + t217 * t208;
t287 = t221 ^ 2 - t224 ^ 2;
t282 = qJD(6) * t144;
t279 = pkin(9) * t302;
t3 = pkin(5) * t105 - pkin(9) * t305 + t10;
t7 = -pkin(9) * t306 + t11;
t278 = -t219 * t7 + t222 * t3;
t273 = -pkin(2) * t213 - t177;
t271 = t32 * t216 + t237 * t41;
t13 = -t216 * t47 + t218 * t54;
t40 = -t216 * t69 + t218 * t78;
t48 = -t216 * t76 + t218 * t92;
t46 = t217 * t85 - t309 * t86;
t57 = -t101 * t216 + t218 * t99;
t268 = pkin(1) * t329;
t75 = t120 * t217 + t270;
t261 = pkin(5) * t237 - pkin(9) * t301;
t260 = -t181 * t105 - t332 * t133;
t206 = -t309 * pkin(3) - pkin(4);
t259 = t219 * t3 + t222 * t7;
t257 = t70 * t139 + t237 * t71;
t256 = -t218 * t32 - t237 * t40;
t20 = -pkin(5) * t139 - pkin(9) * t124 + t40;
t24 = -pkin(9) * t122 + t41;
t4 = t222 * t20 - t219 * t24;
t5 = t20 * t219 + t222 * t24;
t255 = t216 * t40 - t218 * t41;
t36 = pkin(5) * t143 - pkin(9) * t299 + t57;
t44 = -pkin(9) * t300 + t58;
t254 = -t219 * t44 + t222 * t36;
t253 = t219 * t36 + t222 * t44;
t252 = t40 * t301 + t41 * t302 + t258;
t160 = t309 * t208 - t200;
t251 = t100 * t106 + t144 * t32;
t250 = t122 * t218 - t124 * t216;
t246 = -qJD(6) * t124 - t306;
t155 = -pkin(4) - t160;
t154 = qJ(5) + t161;
t148 = (-pkin(9) - t154) * t216;
t244 = -qJD(6) * t148 - t279 + t327;
t149 = t154 * t218 + t212;
t243 = qJD(6) * t149 + t261 + t328;
t242 = -t191 * t167 - t265;
t204 = pkin(3) * t217 + qJ(5);
t169 = (-pkin(9) - t204) * t216;
t241 = -qJD(5) * t218 - qJD(6) * t169 - t279 + t49;
t170 = t204 * t218 + t212;
t240 = qJD(5) * t216 + qJD(6) * t170 + t261 + t48;
t19 = pkin(5) * t306 + t32;
t59 = t122 * pkin(5) + t68;
t239 = t19 * t181 - t237 * t4 + t332 * t59;
t238 = t19 * t182 + t237 * t5 - t310 * t59;
t234 = t108 * t68 + t251;
t25 = t219 * t246 + t289;
t233 = -t105 * t154 + t106 * t155 + t139 * t151 - t315;
t232 = qJD(5) * t139 - t105 * t204 + t106 * t206 - t315;
t189 = t206 - t320;
t150 = t155 - t320;
t128 = -t167 ^ 2 + t168 ^ 2;
t113 = -t168 * t213 - t228;
t112 = -t167 * t213 + t141;
t104 = t181 * t144;
t103 = t182 * t144;
t66 = pkin(5) * t300 + t100;
t60 = t130 + t75;
t43 = t108 * t182 + t281 * t299 - t282 * t294;
t42 = -t108 * t181 - t182 * t282;
t23 = pkin(5) * t304 + t46;
t16 = t260 + t317;
t15 = t316 - t323;
t12 = -pkin(9) * t304 + t14;
t8 = pkin(5) * t107 - t108 * t212 + t13;
t6 = t182 * t25 + t249 * t310;
t1 = -t181 * t25 - t182 * t26 + t332 * t249 + t310 * t80;
t2 = [0, 0, 0, 0.2e1 * t221 * t274, t287 * t329, t290, -t291, 0, -pkin(7) * t290 + t221 * t268, pkin(7) * t291 + t224 * t268, t141 * t184 + t145 * t168, -t141 * t183 - t145 * t167 + t168 * t146 - t184 * t228, -t145 * t213, -t146 * t213, 0, -t167 * t211 + t191 * t146 + t230 * t213 + (t209 * t235 + (t221 * pkin(2) * t183 + t184 * t209) * qJD(2)) * qJD(1), t209 * t141 - t191 * t145 - t236 * t213 + (-t168 + t325) * t211, -t101 * t105 - t107 * t71 - t108 * t70 + t139 * t47 - t143 * t33 + t237 * t46 + t251, t33 * t101 + t147 * t272 + t227 * t245 - t70 * t46 + t71 * t47 + t314, t10 * t143 + t105 * t57 + t107 * t40 + t122 * t46 - t13 * t139 + t216 * t234, -t105 * t58 - t107 * t41 - t11 * t143 + t124 * t46 + t139 * t14 + t218 * t234, -t122 * t14 - t124 * t13 + (-t10 * t144 - t106 * t57 - t108 * t40) * t218 + (-t106 * t58 - t108 * t41 - t11 * t144) * t216, t10 * t57 + t11 * t58 + t13 * t40 + t14 * t41 + t46 * t68 + t314, -t104 * t25 - t249 * t42, -t103 * t25 + t104 * t26 + t249 * t43 - t42 * t80, -t104 * t105 - t107 * t249 + t133 * t42 + t143 * t25, -t103 * t105 - t107 * t80 - t133 * t43 - t143 * t26, t105 * t143 + t107 * t133 (-t12 * t219 + t222 * t8) * t133 + t254 * t105 + t278 * t143 + t4 * t107 + t23 * t80 + t66 * t26 + t19 * t103 + t59 * t43 + (-t133 * t253 - t143 * t5) * qJD(6) -(t222 * t12 + t219 * t8) * t133 - t253 * t105 - t259 * t143 - t5 * t107 - t23 * t249 + t66 * t25 - t19 * t104 + t59 * t42 + (-t133 * t254 - t143 * t4) * qJD(6); 0, 0, 0, -t221 * t292, t287 * t226, 0, 0, 0, t226 * pkin(1) * t221, pkin(1) * t292, t296, t128, t112, t113, 0, t167 * t210 + t295 - t264 * t213 + (t273 * t220 - t175) * qJD(3) + t266, t168 * t210 + t288 * t213 + (t273 * qJD(3) + t178) * t223 + t242, -t105 * t161 - t106 * t160 + t139 * t311 + t237 * t312 + t257, t33 * t161 - t32 * t160 - t147 * (t210 - t321) + t311 * t71 - t312 * t70, t312 * t122 + t139 * t50 + t233 * t216 + t256, t312 * t124 - t139 * t51 + t233 * t218 + t271, t122 * t327 + t124 * t328 + t252, -t255 * t151 + t258 * t154 + t155 * t32 + t312 * t68 - t40 * t50 - t41 * t51, t6, t1, t15, t16, -t303 (t148 * t222 - t149 * t219) * t105 + t150 * t26 + t313 * t80 + (t219 * t244 - t222 * t243) * t133 + t239 -(t148 * t219 + t149 * t222) * t105 + t150 * t25 - t313 * t249 + (t219 * t243 + t222 * t244) * t133 + t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t296, t128, t112, t113, 0, -t213 * t248 + t231 + t295, t213 * t267 + t242 - t326, -t76 * t139 - t75 * t237 + (-t105 * t217 - t309 * t106) * pkin(3) + t257, t70 * t75 - t71 * t76 + (t147 * t168 + t217 * t33 - t309 * t32) * pkin(3), -t122 * t75 + t139 * t48 + t216 * t232 + t256, -t124 * t75 - t139 * t49 + t218 * t232 + t271, -qJD(5) * t250 + t122 * t49 + t124 * t48 + t252, -qJD(5) * t255 + t204 * t258 + t206 * t32 - t40 * t48 - t41 * t49 - t68 * t75, t6, t1, t15, t16, -t303 (t169 * t222 - t170 * t219) * t105 + t189 * t26 - t60 * t80 + (t219 * t241 - t222 * t240) * t133 + t239 -(t169 * t219 + t170 * t222) * t105 + t189 * t25 + t60 * t249 + (t219 * t240 + t222 * t241) * t133 + t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t139 ^ 2 - t237 ^ 2, -t71 * t139 + t237 * t70 + t227, t105 * t218 - t122 * t237 - t139 * t302, -t105 * t216 - t124 * t237 - t139 * t301, t250 * t139 + (-t216 ^ 2 - t218 ^ 2) * t106, t10 * t218 + t11 * t216 + t139 * t255 - t237 * t68, 0, 0, 0, 0, 0, t260 - t317, t316 + t323; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t124 * t139 + t306, t122 * t139 + t305, -t122 ^ 2 - t124 ^ 2, t122 * t41 + t124 * t40 + t32, 0, 0, 0, 0, 0, t26 - t331, -t133 * t330 + (-t124 * t133 + t246) * t219 + t289; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t249 * t80, t249 ^ 2 - t80 ^ 2, t133 * t80 + t25, -t26 - t331, t105, t249 * t59 + t324 * t5 + t278, t324 * t4 + t59 * t80 - t259;];
tauc_reg  = t2;
