% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% tauc_reg [6x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRRPR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:20:28
% EndTime: 2019-03-08 23:20:41
% DurationCPUTime: 4.92s
% Computational Cost: add. (4952->406), mult. (12628->593), div. (0->0), fcn. (9658->12), ass. (0->226)
t208 = sin(qJ(4));
t212 = cos(qJ(4));
t268 = t212 * qJD(3);
t209 = sin(qJ(3));
t277 = qJD(2) * t209;
t172 = t208 * t277 - t268;
t275 = qJD(3) * t208;
t174 = t212 * t277 + t275;
t203 = sin(pkin(12));
t205 = cos(pkin(12));
t109 = t172 * t203 - t174 * t205;
t207 = sin(qJ(6));
t211 = cos(qJ(6));
t269 = qJD(6) * t207;
t229 = -t172 * t205 - t174 * t203;
t291 = t211 * t229;
t265 = qJD(3) * qJD(4);
t271 = qJD(4) * t209;
t250 = t208 * t271;
t213 = cos(qJ(3));
t252 = t213 * t268;
t329 = -t250 + t252;
t131 = qJD(2) * t329 + t212 * t265;
t270 = qJD(4) * t212;
t249 = t209 * t270;
t273 = qJD(3) * t213;
t253 = t208 * t273;
t330 = t249 + t253;
t132 = qJD(2) * t330 + t208 * t265;
t76 = -t131 * t203 - t132 * t205;
t77 = t131 * t205 - t132 * t203;
t12 = qJD(6) * t291 + t109 * t269 + t207 * t76 + t211 * t77;
t276 = qJD(2) * t213;
t192 = -qJD(4) + t276;
t186 = -qJD(6) + t192;
t52 = t109 * t207 + t291;
t309 = t186 * t52;
t346 = t12 + t309;
t210 = sin(qJ(2));
t204 = sin(pkin(6));
t280 = qJD(1) * t204;
t248 = t210 * t280;
t178 = qJD(2) * pkin(8) + t248;
t206 = cos(pkin(6));
t279 = qJD(1) * t213;
t137 = -t209 * t178 + t206 * t279;
t127 = -qJD(3) * pkin(3) - t137;
t96 = pkin(4) * t172 + qJD(5) + t127;
t45 = -pkin(5) * t229 + t96;
t345 = t45 * t52;
t328 = -t211 * t109 + t207 * t229;
t344 = t328 * t52;
t236 = pkin(3) * t209 - pkin(9) * t213;
t176 = t236 * qJD(3);
t181 = -pkin(3) * t213 - pkin(9) * t209 - pkin(2);
t214 = cos(qJ(2));
t289 = t213 * t214;
t343 = -(t208 * t210 + t212 * t289) * t280 + t208 * t176 + t181 * t270;
t274 = qJD(3) * t209;
t321 = pkin(8) * t208;
t342 = t212 * t176 + t274 * t321 - (-t208 * t289 + t210 * t212) * t280;
t341 = t328 ^ 2 - t52 ^ 2;
t13 = qJD(6) * t328 + t207 * t77 - t211 * t76;
t310 = t186 * t328;
t339 = -t13 - t310;
t290 = t212 * t213;
t194 = pkin(8) * t290;
t226 = pkin(4) * t209 - qJ(5) * t290;
t267 = t212 * qJD(5);
t338 = -t209 * t267 + t226 * qJD(3) + (-t194 + (qJ(5) * t209 - t181) * t208) * qJD(4) + t342;
t292 = t209 * t212;
t337 = -(-pkin(8) * qJD(3) - qJ(5) * qJD(4)) * t292 - (-qJD(5) * t209 + (-pkin(8) * qJD(4) - qJ(5) * qJD(3)) * t213) * t208 - t343;
t319 = -qJ(5) - pkin(9);
t243 = qJD(4) * t319;
t254 = t208 * t276;
t175 = t236 * qJD(2);
t286 = t212 * t137 + t208 * t175;
t336 = qJ(5) * t254 + t208 * t243 + t267 - t286;
t241 = -t137 * t208 + t212 * t175;
t335 = -qJD(2) * t226 - t208 * qJD(5) + t212 * t243 - t241;
t266 = qJD(2) * qJD(3);
t196 = t209 * t266;
t278 = qJD(2) * t204;
t255 = t214 * t278;
t101 = -t178 * t274 + (qJD(3) * t206 + t255) * t279;
t136 = (t176 + t248) * qJD(2);
t242 = t208 * t101 - t212 * t136;
t296 = t206 * t209;
t191 = qJD(1) * t296;
t138 = t213 * t178 + t191;
t128 = qJD(3) * pkin(9) + t138;
t258 = t214 * t280;
t139 = qJD(2) * t181 - t258;
t295 = t208 * t139;
t82 = t128 * t212 + t295;
t217 = -qJD(4) * t82 - t242;
t19 = pkin(4) * t196 - qJ(5) * t131 - qJD(5) * t174 + t217;
t272 = qJD(4) * t208;
t222 = t212 * t101 - t128 * t272 + t208 * t136 + t139 * t270;
t21 = -qJ(5) * t132 - qJD(5) * t172 + t222;
t6 = t205 * t19 - t203 * t21;
t2 = pkin(5) * t196 - pkin(10) * t77 + t6;
t7 = t203 * t19 + t205 * t21;
t3 = pkin(10) * t76 + t7;
t259 = t211 * t2 - t207 * t3;
t334 = -t45 * t328 + t259;
t333 = pkin(10) * t109;
t166 = t203 * t212 + t205 * t208;
t221 = t213 * t166;
t332 = qJD(2) * t221 - t166 * qJD(4);
t228 = t203 * t208 - t205 * t212;
t331 = t192 * t228;
t327 = pkin(10) * t229;
t317 = t337 * t203 + t338 * t205;
t316 = t338 * t203 - t337 * t205;
t313 = -t336 * t203 + t335 * t205;
t312 = t335 * t203 + t336 * t205;
t325 = -t138 + (-t254 + t272) * pkin(4);
t62 = -qJ(5) * t172 + t82;
t308 = t205 * t62;
t80 = -t128 * t208 + t212 * t139;
t61 = -qJ(5) * t174 + t80;
t43 = -pkin(4) * t192 + t61;
t23 = t203 * t43 + t308;
t14 = t23 + t327;
t246 = -t14 * t269 + t207 * t2;
t56 = t203 * t62;
t22 = t205 * t43 - t56;
t9 = -pkin(5) * t192 + t22 + t333;
t324 = (qJD(6) * t9 + t3) * t211 + t246;
t93 = t166 * t271 + t203 * t253 - t205 * t252;
t323 = -pkin(5) * t274 - pkin(10) * t93 - t317;
t322 = pkin(4) * t203;
t320 = t211 * t9;
t92 = -qJD(3) * t221 + t228 * t271;
t318 = pkin(10) * t92 + t316;
t230 = -t166 * t207 - t211 * t228;
t315 = qJD(6) * t230 + t207 * t332 + t211 * t331;
t108 = t166 * t211 - t207 * t228;
t314 = qJD(6) * t108 + t207 * t331 - t211 * t332;
t29 = t205 * t61 - t56;
t311 = qJD(2) * pkin(2);
t307 = -pkin(5) * t332 + t325;
t240 = t209 * t255;
t102 = qJD(1) * t240 + qJD(3) * t191 + t178 * t273;
t306 = t102 * t208;
t305 = t102 * t212;
t304 = t127 * t208;
t303 = t131 * t208;
t302 = t172 * t192;
t301 = t174 * t192;
t300 = t192 * t212;
t299 = t204 * t210;
t298 = t204 * t214;
t216 = qJD(2) ^ 2;
t297 = t204 * t216;
t294 = t208 * t209;
t293 = t208 * t213;
t215 = qJD(3) ^ 2;
t288 = t215 * t209;
t287 = t215 * t213;
t168 = t212 * t181;
t113 = -qJ(5) * t292 + t168 + (-pkin(4) - t321) * t213;
t283 = t208 * t181 + t194;
t120 = -qJ(5) * t294 + t283;
t64 = t203 * t113 + t205 * t120;
t182 = t319 * t208;
t183 = t319 * t212;
t125 = t203 * t182 - t205 * t183;
t282 = pkin(4) * t294 + t209 * pkin(8);
t201 = t209 ^ 2;
t281 = -t213 ^ 2 + t201;
t263 = t210 * t297;
t261 = pkin(4) * t330 + pkin(8) * t273;
t260 = -pkin(4) * t212 - pkin(3);
t256 = t210 * t278;
t251 = t192 * t272;
t28 = -t203 * t61 - t308;
t63 = t205 * t113 - t120 * t203;
t124 = t205 * t182 + t183 * t203;
t239 = t213 * t255;
t98 = -pkin(10) * t228 + t125;
t238 = pkin(5) * t277 + pkin(10) * t331 + qJD(6) * t98 - t313;
t97 = -pkin(10) * t166 + t124;
t237 = pkin(10) * t332 + qJD(6) * t97 + t312;
t179 = -t258 - t311;
t235 = -t179 - t258;
t5 = t211 * t14 + t207 * t9;
t149 = t228 * t209;
t37 = -pkin(5) * t213 + pkin(10) * t149 + t63;
t148 = t166 * t209;
t38 = -pkin(10) * t148 + t64;
t234 = t207 * t37 + t211 * t38;
t157 = t213 * t299 + t296;
t118 = -t157 * t208 - t212 * t298;
t223 = -t157 * t212 + t208 * t298;
t66 = t118 * t205 + t203 * t223;
t67 = t118 * t203 - t205 * t223;
t233 = -t207 * t67 + t211 * t66;
t232 = t207 * t66 + t211 * t67;
t231 = -t211 * t148 + t149 * t207;
t91 = -t148 * t207 - t149 * t211;
t227 = qJD(2) * t201 - t192 * t213;
t68 = pkin(4) * t132 + t102;
t195 = pkin(4) * t205 + pkin(5);
t225 = t195 * t207 + t211 * t322;
t224 = t195 * t211 - t207 * t322;
t156 = -t206 * t213 + t209 * t299;
t219 = qJD(3) * (-t235 - t311);
t143 = pkin(5) * t228 + t260;
t117 = qJD(3) * t157 + t240;
t116 = -qJD(3) * t156 + t239;
t114 = pkin(5) * t148 + t282;
t86 = pkin(4) * t174 - pkin(5) * t109;
t65 = -pkin(5) * t92 + t261;
t55 = qJD(4) * t118 + t116 * t212 + t208 * t256;
t54 = qJD(4) * t223 - t116 * t208 + t212 * t256;
t34 = -pkin(5) * t76 + t68;
t33 = qJD(6) * t91 - t207 * t93 - t211 * t92;
t32 = qJD(6) * t231 + t207 * t92 - t211 * t93;
t27 = t203 * t54 + t205 * t55;
t25 = -t203 * t55 + t205 * t54;
t16 = t29 + t333;
t15 = t28 - t327;
t4 = -t14 * t207 + t320;
t1 = [0, 0, -t263, -t214 * t297, 0, 0, 0, 0, 0, -t213 * t263 + (-t117 - t240) * qJD(3), t209 * t263 + (-t116 - t239) * qJD(3), 0, 0, 0, 0, 0, t117 * t172 + t118 * t196 + t132 * t156 - t192 * t54, t117 * t174 + t131 * t156 + t192 * t55 + t196 * t223, t109 * t25 + t229 * t27 - t66 * t77 + t67 * t76, t117 * t96 + t156 * t68 + t22 * t25 + t23 * t27 + t6 * t66 + t67 * t7, 0, 0, 0, 0, 0 -(-qJD(6) * t232 - t207 * t27 + t211 * t25) * t186 + t233 * t196 - t117 * t52 + t156 * t13 (qJD(6) * t233 + t207 * t25 + t211 * t27) * t186 - t232 * t196 + t117 * t328 + t156 * t12; 0, 0, 0, 0, 0.2e1 * t213 * t196, -0.2e1 * t281 * t266, t287, -t288, 0, -pkin(8) * t287 + t209 * t219, pkin(8) * t288 + t213 * t219, t131 * t292 + t174 * t329 (-t172 * t212 - t174 * t208) * t273 + (-t303 - t132 * t212 + (t172 * t208 - t174 * t212) * qJD(4)) * t209, t192 * t250 - t131 * t213 + (t174 * t209 + t212 * t227) * qJD(3), t192 * t249 + t132 * t213 + (-t172 * t209 - t208 * t227) * qJD(3) (-t192 - t276) * t274 (t181 * t272 - t342) * t192 + ((pkin(8) * t172 + t304) * qJD(3) + (t295 + (pkin(8) * t192 + t128) * t212) * qJD(4) + t242) * t213 + (-t172 * t258 + t127 * t270 + pkin(8) * t132 + t306 + ((-pkin(8) * t293 + t168) * qJD(2) + t80) * qJD(3)) * t209, t343 * t192 + (t127 * t268 + (qJD(3) * t174 - t251) * pkin(8) + t222) * t213 + (-t174 * t258 - t127 * t272 + pkin(8) * t131 + t305 + (-pkin(8) * t300 - qJD(2) * t283 - t82) * qJD(3)) * t209, t109 * t317 - t148 * t7 + t149 * t6 + t22 * t93 + t229 * t316 + t23 * t92 - t63 * t77 + t64 * t76, t7 * t64 + t6 * t63 + t68 * t282 + (-t209 * t258 + t261) * t96 + t316 * t23 + t317 * t22, t12 * t91 + t32 * t328, t12 * t231 - t13 * t91 + t32 * t52 - t328 * t33, -t12 * t213 - t186 * t32 + (qJD(2) * t91 + t328) * t274, t13 * t213 + t186 * t33 + (qJD(2) * t231 + t52) * t274 (-t186 - t276) * t274, -t259 * t213 - t65 * t52 + t114 * t13 - t34 * t231 + t45 * t33 + (t318 * t207 + t323 * t211) * t186 + (t186 * t234 + t213 * t5) * qJD(6) + (t52 * t258 + ((-t207 * t38 + t211 * t37) * qJD(2) + t4) * qJD(3)) * t209, t324 * t213 + t65 * t328 + t114 * t12 + t34 * t91 + t45 * t32 + ((qJD(6) * t37 + t318) * t211 + (-qJD(6) * t38 - t323) * t207) * t186 + (-t328 * t258 + (-qJD(2) * t234 - t5) * qJD(3)) * t209; 0, 0, 0, 0, -t209 * t216 * t213, t281 * t216, 0, 0, 0, qJD(3) * t138 - t179 * t277 - t102, t235 * t276, -t174 * t300 + t303 (t131 + t302) * t212 + (-t132 + t301) * t208, -t192 * t270 + (t192 * t290 + (-t174 + t275) * t209) * qJD(2), t251 + (-t192 * t293 + (t172 + t268) * t209) * qJD(2), t192 * t277, -pkin(3) * t132 - t305 + t241 * t192 - t138 * t172 + (pkin(9) * t300 + t304) * qJD(4) + (-t80 * t209 + (-pkin(9) * t274 - t127 * t213) * t208) * qJD(2), -pkin(3) * t131 + t306 - t286 * t192 - t138 * t174 + (-pkin(9) * t192 * t208 + t127 * t212) * qJD(4) + (-t127 * t290 + (-pkin(9) * t268 + t82) * t209) * qJD(2), t313 * t109 - t124 * t77 + t125 * t76 - t166 * t6 - t22 * t331 - t228 * t7 + t312 * t229 + t23 * t332, t6 * t124 + t7 * t125 + t313 * t22 + t312 * t23 + t68 * t260 + t325 * t96, t12 * t108 + t315 * t328, -t108 * t13 + t12 * t230 - t314 * t328 + t315 * t52, -t315 * t186 + (qJD(3) * t108 - t328) * t277, t314 * t186 + (qJD(3) * t230 - t52) * t277, t186 * t277, -t34 * t230 + t143 * t13 - t307 * t52 + t314 * t45 + (t207 * t237 + t211 * t238) * t186 + ((-t207 * t98 + t211 * t97) * qJD(3) - t4) * t277, t34 * t108 + t143 * t12 + t307 * t328 + t315 * t45 + (-t207 * t238 + t211 * t237) * t186 + (-(t207 * t97 + t211 * t98) * qJD(3) + t5) * t277; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t174 * t172, -t172 ^ 2 + t174 ^ 2, t131 - t302, -t132 - t301, t196, -t127 * t174 - t192 * t82 + t217, t127 * t172 - t192 * t80 - t222 (t203 * t76 - t205 * t77) * pkin(4) + (-t29 + t22) * t229 + (-t23 - t28) * t109, -t22 * t28 - t23 * t29 + (-t174 * t96 + t203 * t7 + t205 * t6) * pkin(4), -t344, t341, t346, t339, t196, t224 * t196 + (t15 * t211 - t16 * t207) * t186 + t86 * t52 + (t186 * t225 - t5) * qJD(6) + t334, -t225 * t196 - t211 * t3 - (t15 * t207 + t16 * t211) * t186 - t86 * t328 - t345 + (t186 * t224 - t320) * qJD(6) - t246; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109 ^ 2 - t229 ^ 2, -t109 * t22 - t229 * t23 + t68, 0, 0, 0, 0, 0, t13 - t310, t12 - t309; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t344, t341, t346, t339, t196 (-qJD(6) - t186) * t5 + t334, -t186 * t4 - t324 - t345;];
tauc_reg  = t1;
