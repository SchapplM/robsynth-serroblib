% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% 
% Output:
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRPR11_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR11_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:15:03
% EndTime: 2019-03-09 11:15:16
% DurationCPUTime: 4.61s
% Computational Cost: add. (5006->402), mult. (11596->566), div. (0->0), fcn. (7769->8), ass. (0->230)
t218 = cos(qJ(4));
t215 = sin(qJ(4));
t286 = qJD(2) * t215;
t219 = cos(qJ(2));
t287 = qJD(1) * t219;
t161 = t218 * t287 + t286;
t259 = t215 * t287;
t284 = qJD(2) * t218;
t163 = -t259 + t284;
t212 = sin(pkin(10));
t213 = cos(pkin(10));
t101 = t161 * t212 - t163 * t213;
t214 = sin(qJ(6));
t217 = cos(qJ(6));
t240 = -t161 * t213 - t163 * t212;
t278 = qJD(6) * t217;
t279 = qJD(6) * t214;
t216 = sin(qJ(2));
t273 = qJD(1) * qJD(2);
t257 = t216 * t273;
t114 = -qJD(4) * t161 + t215 * t257;
t281 = qJD(4) * t218;
t115 = qJD(2) * t281 - qJD(4) * t259 - t218 * t257;
t67 = -t114 * t212 - t115 * t213;
t68 = t114 * t213 - t115 * t212;
t12 = t101 * t279 + t214 * t67 + t217 * t68 + t240 * t278;
t288 = qJD(1) * t216;
t193 = qJD(4) + t288;
t188 = qJD(6) + t193;
t49 = t101 * t214 + t217 * t240;
t315 = t188 * t49;
t343 = t12 - t315;
t332 = -t217 * t101 + t214 * t240;
t342 = t332 * t49;
t341 = t332 ^ 2 - t49 ^ 2;
t220 = -pkin(2) - pkin(8);
t256 = -t216 * qJ(3) - pkin(1);
t157 = t220 * t219 + t256;
t128 = t157 * qJD(1);
t199 = pkin(7) * t288;
t333 = qJD(3) + t199;
t275 = pkin(3) * t288 + t333;
t132 = t220 * qJD(2) + t275;
t77 = t128 * t218 + t132 * t215;
t63 = -qJ(5) * t161 + t77;
t314 = t213 * t63;
t76 = -t128 * t215 + t218 * t132;
t62 = -qJ(5) * t163 + t76;
t54 = pkin(4) * t193 + t62;
t25 = t212 * t54 + t314;
t330 = pkin(9) * t240;
t16 = t25 + t330;
t15 = t16 * t279;
t198 = t219 * t273;
t192 = pkin(2) * t257;
t243 = pkin(8) * t216 - qJ(3) * t219;
t277 = t216 * qJD(3);
t226 = qJD(2) * t243 - t277;
t109 = qJD(1) * t226 + t192;
t191 = pkin(7) * t198;
t153 = pkin(3) * t198 + t191;
t250 = -t215 * t109 + t218 * t153;
t224 = -qJD(4) * t77 + t250;
t21 = pkin(4) * t198 - qJ(5) * t114 - qJD(5) * t163 + t224;
t270 = -t218 * t109 - t132 * t281 - t215 * t153;
t282 = qJD(4) * t215;
t228 = -t128 * t282 - t270;
t23 = -qJ(5) * t115 - qJD(5) * t161 + t228;
t6 = t213 * t21 - t212 * t23;
t4 = pkin(5) * t198 - pkin(9) * t68 + t6;
t200 = pkin(7) * t287;
t169 = pkin(3) * t287 + t200;
t209 = qJD(2) * qJ(3);
t147 = t209 + t169;
t106 = pkin(4) * t161 + qJD(5) + t147;
t55 = -pkin(5) * t240 + t106;
t340 = -t214 * t4 - t55 * t49 + t15;
t13 = qJD(6) * t332 + t214 * t68 - t217 * t67;
t312 = t332 * t188;
t338 = -t13 + t312;
t204 = pkin(2) * t288;
t140 = qJD(1) * t243 + t204;
t249 = -t215 * t140 + t218 * t169;
t276 = t218 * qJD(5);
t294 = qJ(5) - t220;
t300 = t215 * t216;
t337 = t294 * t282 - t276 - (pkin(4) * t219 - qJ(5) * t300) * qJD(1) - t249;
t173 = t294 * t218;
t266 = t218 * t288;
t291 = t218 * t140 + t215 * t169;
t336 = -qJ(5) * t266 - qJD(4) * t173 - t215 * qJD(5) - t291;
t239 = t212 * t218 + t213 * t215;
t292 = t193 * t239;
t7 = t212 * t21 + t213 * t23;
t5 = pkin(9) * t67 + t7;
t267 = -t214 * t5 + t217 * t4;
t335 = -t55 * t332 + t267;
t322 = pkin(3) + pkin(7);
t334 = pkin(9) * t101;
t260 = t213 * t281;
t301 = t212 * t215;
t293 = -t212 * t282 + t213 * t266 - t288 * t301 + t260;
t331 = -0.2e1 * t273;
t318 = -t336 * t212 + t337 * t213;
t317 = t337 * t212 + t336 * t213;
t327 = t292 * t217;
t181 = t322 * t216;
t290 = t218 * t157 + t215 * t181;
t238 = -t213 * t218 + t301;
t100 = -t214 * t239 - t217 * t238;
t280 = qJD(4) * t219;
t263 = t215 * t280;
t268 = pkin(4) * t218 + pkin(3);
t285 = qJD(2) * t216;
t326 = (-pkin(7) - t268) * t285 - pkin(4) * t263;
t325 = t216 * t284 + t263;
t324 = pkin(4) * t281 + t268 * t288 + t333;
t56 = t212 * t63;
t24 = t213 * t54 - t56;
t323 = -t238 * t6 + t239 * t7 - t292 * t24 + t293 * t25;
t321 = pkin(4) * t212;
t203 = pkin(2) * t285;
t117 = t203 + t226;
t283 = qJD(2) * t219;
t170 = t322 * t283;
t152 = t218 * t170;
t251 = qJ(5) * t219 - t157;
t36 = pkin(4) * t283 + t152 + t251 * t281 + (-qJ(5) * t285 - qJD(4) * t181 + qJD(5) * t219 - t117) * t215;
t227 = t218 * t117 - t157 * t282 + t215 * t170 + t181 * t281;
t41 = qJ(5) * t325 - t219 * t276 + t227;
t11 = t212 * t36 + t213 * t41;
t99 = -t214 * t238 + t217 * t239;
t320 = -qJD(6) * t99 - t214 * t293 - t327;
t319 = qJD(6) * t100 - t214 * t292 + t217 * t293;
t30 = t213 * t62 - t56;
t165 = t218 * t181;
t90 = t216 * pkin(4) + t215 * t251 + t165;
t299 = t218 * t219;
t95 = -qJ(5) * t299 + t290;
t43 = t212 * t90 + t213 * t95;
t316 = qJD(2) * pkin(2);
t14 = pkin(5) * t193 + t24 + t334;
t313 = t217 * t14;
t311 = pkin(5) * t293 + t324;
t310 = t114 * t218;
t168 = t322 * t285;
t208 = qJD(2) * qJD(3);
t137 = -qJD(1) * t168 + t208;
t309 = t137 * t215;
t308 = t137 * t218;
t306 = t161 * t193;
t305 = t163 * t193;
t304 = t163 * t219;
t303 = t193 * t216;
t302 = t193 * t220;
t222 = qJD(1) ^ 2;
t298 = t219 * t222;
t221 = qJD(2) ^ 2;
t297 = t221 * t216;
t296 = t221 * t219;
t295 = t215 * pkin(4) + qJ(3);
t172 = t294 * t215;
t108 = -t213 * t172 - t212 * t173;
t182 = t322 * t219;
t210 = t216 ^ 2;
t211 = t219 ^ 2;
t289 = t210 - t211;
t175 = -pkin(2) * t219 + t256;
t148 = qJD(1) * t175;
t272 = t218 * t303;
t271 = t216 * t298;
t269 = pkin(4) * t299 + t182;
t265 = t215 * t285;
t262 = t193 * t281;
t261 = t218 * t280;
t255 = qJD(6) * t14 + t5;
t10 = -t212 * t41 + t213 * t36;
t29 = -t212 * t62 - t314;
t42 = -t212 * t95 + t213 * t90;
t253 = pkin(1) * t331;
t252 = qJD(3) - t316;
t107 = t172 * t212 - t213 * t173;
t85 = -pkin(9) * t239 + t108;
t247 = pkin(5) * t287 - pkin(9) * t292 + qJD(6) * t85 - t318;
t84 = pkin(9) * t238 + t107;
t246 = pkin(9) * t293 - qJD(6) * t84 - t317;
t244 = qJD(6) * t238 - t293;
t2 = t214 * t14 + t217 * t16;
t136 = t239 * t219;
t32 = pkin(5) * t216 + pkin(9) * t136 + t42;
t135 = t238 * t219;
t33 = pkin(9) * t135 + t43;
t242 = t214 * t32 + t217 * t33;
t241 = t217 * t135 + t136 * t214;
t87 = t135 * t214 - t136 * t217;
t237 = -qJD(1) * t211 + t303;
t236 = -0.2e1 * qJD(2) * t148;
t235 = t193 * t215;
t197 = pkin(4) * t213 + pkin(5);
t234 = t197 * t214 + t217 * t321;
t233 = t197 * t217 - t214 * t321;
t229 = -qJ(3) * t283 - t277;
t125 = qJD(1) * t229 + t192;
t142 = t203 + t229;
t232 = pkin(7) * t221 + qJD(1) * t142 + t125;
t231 = t147 * t216 + t220 * t283;
t88 = t115 * pkin(4) + t137;
t171 = pkin(7) * t257 - t208;
t174 = t199 + t252;
t180 = -t200 - t209;
t223 = -t171 * t219 + (t174 * t219 + (t180 + t200) * t216) * qJD(2);
t186 = t218 * t198;
t184 = t216 * t198;
t166 = -qJ(3) * t287 + t204;
t129 = t148 * t288;
t127 = pkin(5) * t239 + t295;
t98 = -pkin(5) * t135 + t269;
t92 = -t212 * t325 - t213 * t265 + t219 * t260;
t91 = qJD(4) * t136 - t238 * t285;
t72 = pkin(4) * t163 - pkin(5) * t101;
t53 = -t91 * pkin(5) + t326;
t40 = -t67 * pkin(5) + t88;
t27 = qJD(6) * t87 - t214 * t92 - t217 * t91;
t26 = qJD(6) * t241 + t214 * t91 - t217 * t92;
t18 = t30 + t334;
t17 = t29 - t330;
t9 = pkin(9) * t91 + t11;
t8 = pkin(5) * t283 + pkin(9) * t92 + t10;
t1 = -t16 * t214 + t313;
t3 = [0, 0, 0, 0.2e1 * t184, t289 * t331, t296, -t297, 0, -pkin(7) * t296 + t216 * t253, pkin(7) * t297 + t219 * t253, t223, t216 * t236 + t219 * t232, -t216 * t232 + t219 * t236, pkin(7) * t223 + t125 * t175 + t142 * t148, -t114 * t215 * t219 + (-t261 + t265) * t163 (-t161 * t215 + t163 * t218) * t285 + (-t310 + t115 * t215 + (t161 * t218 + t163 * t215) * qJD(4)) * t219, -t193 * t261 + t114 * t216 + (t215 * t237 + t304) * qJD(2), t193 * t263 - t115 * t216 + (-t161 * t219 + t218 * t237) * qJD(2), t193 * t283 + t184 (-t215 * t117 + t152) * t193 - t168 * t161 + t182 * t115 + (-t147 * t284 + t250) * t216 + (-t193 * t290 - t77 * t216) * qJD(4) + (-t147 * t282 + t308 + ((-t157 * t215 + t165) * qJD(1) + t76) * qJD(2)) * t219, -t227 * t193 - t168 * t163 + t182 * t114 + ((qJD(2) * t147 + qJD(4) * t128) * t215 + t270) * t216 + (-t147 * t281 - t309 + (-qJD(1) * t290 - t77) * qJD(2)) * t219, t10 * t101 + t11 * t240 + t135 * t7 + t136 * t6 + t24 * t92 + t25 * t91 - t42 * t68 + t43 * t67, t24 * t10 + t106 * t326 + t25 * t11 + t269 * t88 + t6 * t42 + t7 * t43, t12 * t87 + t26 * t332, t12 * t241 - t13 * t87 + t26 * t49 - t27 * t332, t12 * t216 + t26 * t188 + (qJD(1) * t87 + t332) * t283, -t13 * t216 - t27 * t188 + (qJD(1) * t241 + t49) * t283, t188 * t283 + t184 (-t214 * t9 + t217 * t8) * t188 + t267 * t216 - t53 * t49 + t98 * t13 - t40 * t241 + t55 * t27 + (-t188 * t242 - t2 * t216) * qJD(6) + ((-t214 * t33 + t217 * t32) * qJD(1) + t1) * t283, t98 * t12 + t15 * t216 + t55 * t26 + t40 * t87 + t53 * t332 + (-(-qJD(6) * t33 + t8) * t188 - t4 * t216) * t214 + (-(qJD(6) * t32 + t9) * t188 - t255 * t216) * t217 + (-qJD(1) * t242 - t2) * t283; 0, 0, 0, -t271, t289 * t222, 0, 0, 0, t222 * pkin(1) * t216, pkin(1) * t298 ((-t180 - t209) * t216 + (-t174 + t252) * t219) * qJD(1), -t166 * t287 + t129, 0.2e1 * t208 + (t148 * t219 + t166 * t216) * qJD(1), -qJ(3) * t171 - qJD(3) * t180 - t148 * t166 + (-t180 * t216 + (-t174 - t316) * t219) * qJD(1) * pkin(7), -t163 * t235 + t310 (-t115 - t305) * t218 + (-t114 + t306) * t215, -t193 * t282 + t186 + (-t193 * t300 - t304) * qJD(1), -t262 + (-t272 + (t161 - t286) * t219) * qJD(1), -t193 * t287, qJ(3) * t115 + t309 - t249 * t193 + t275 * t161 + (t147 * t218 - t215 * t302) * qJD(4) + (t218 * t231 - t76 * t219) * qJD(1), qJ(3) * t114 + t308 + t291 * t193 + t275 * t163 + (-t147 * t215 - t218 * t302) * qJD(4) + (-t215 * t231 + t77 * t219) * qJD(1), t101 * t318 - t107 * t68 + t108 * t67 + t240 * t317 - t323, t106 * t324 + t6 * t107 + t7 * t108 + t318 * t24 + t317 * t25 + t88 * t295, t12 * t100 + t320 * t332, -t100 * t13 - t12 * t99 - t319 * t332 + t320 * t49, t320 * t188 + (qJD(2) * t100 - t332) * t287, -t319 * t188 + (-qJD(2) * t99 - t49) * t287, -t188 * t287, t127 * t13 + t40 * t99 + t319 * t55 - t311 * t49 + (t214 * t246 - t217 * t247) * t188 + ((-t214 * t85 + t217 * t84) * qJD(2) - t1) * t287, t40 * t100 + t127 * t12 + t320 * t55 + t311 * t332 + (t214 * t247 + t217 * t246) * t188 + (-(t214 * t84 + t217 * t85) * qJD(2) + t2) * t287; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t271, -t210 * t222 - t221, qJD(2) * t180 + t129 + t191, 0, 0, 0, 0, 0, -qJD(2) * t161 - t193 * t235 + t186, -t262 - qJD(2) * t163 + (-t215 * t283 - t272) * qJD(1), -t101 * t292 + t238 * t68 + t239 * t67 + t240 * t293, -qJD(2) * t106 + t323, 0, 0, 0, 0, 0 (t214 * t244 - t239 * t278 - t327) * t188 + (t100 * t287 + t49) * qJD(2) (t244 * t217 + (qJD(6) * t239 + t292) * t214) * t188 + (-t287 * t99 - t332) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t163 * t161, -t161 ^ 2 + t163 ^ 2, t114 + t306, -t115 + t305, t198, -t147 * t163 + t77 * t193 + t224, t147 * t161 + t193 * t76 - t228 (t212 * t67 - t213 * t68) * pkin(4) + (-t30 + t24) * t240 + (-t25 - t29) * t101, -t24 * t29 - t25 * t30 + (-t106 * t163 + t212 * t7 + t213 * t6) * pkin(4), -t342, t341, t343, t338, t198, t233 * t198 - (t17 * t217 - t18 * t214) * t188 + t72 * t49 + (-t188 * t234 - t2) * qJD(6) + t335, -t234 * t198 - t217 * t5 + (t17 * t214 + t18 * t217) * t188 - t72 * t332 + (-t188 * t233 - t313) * qJD(6) + t340; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101 ^ 2 - t240 ^ 2, -t101 * t24 - t240 * t25 + t88, 0, 0, 0, 0, 0, t13 + t312, t12 + t315; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t342, t341, t343, t338, t198 (-qJD(6) + t188) * t2 + t335, t1 * t188 - t217 * t255 + t340;];
tauc_reg  = t3;
