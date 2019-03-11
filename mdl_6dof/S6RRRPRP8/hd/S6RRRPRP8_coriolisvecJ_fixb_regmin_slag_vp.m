% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRPRP8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:19:40
% EndTime: 2019-03-09 17:19:52
% DurationCPUTime: 4.70s
% Computational Cost: add. (4657->451), mult. (10879->597), div. (0->0), fcn. (7050->6), ass. (0->217)
t197 = sin(qJ(2));
t200 = cos(qJ(2));
t149 = -pkin(2) * t200 - pkin(8) * t197 - pkin(1);
t118 = t149 * qJD(1);
t271 = qJD(1) * t200;
t185 = pkin(7) * t271;
t156 = qJD(2) * pkin(8) + t185;
t196 = sin(qJ(3));
t199 = cos(qJ(3));
t77 = t199 * t118 - t196 * t156;
t281 = qJD(4) - t77;
t266 = qJD(3) * t199;
t342 = -t199 * t271 + t266;
t253 = t196 * t271;
t267 = qJD(3) * t196;
t341 = t253 - t267;
t224 = pkin(2) * t197 - pkin(8) * t200;
t143 = t224 * qJD(1);
t114 = t196 * t143;
t272 = qJD(1) * t197;
t280 = qJ(4) * t272 + t114;
t287 = t197 * t199;
t288 = t196 * t200;
t319 = pkin(8) - pkin(9);
t336 = t319 * t267 + (-pkin(7) * t287 + pkin(9) * t288) * qJD(1) + t280;
t158 = t319 * t199;
t316 = pkin(7) * t196;
t254 = -pkin(3) - t316;
t285 = t199 * t200;
t207 = -pkin(9) * t285 + (-pkin(4) + t254) * t197;
t286 = t199 * t143;
t340 = -qJD(1) * t207 + qJD(3) * t158 + t286;
t263 = t199 * qJD(2);
t131 = t196 * t272 - t263;
t252 = t199 * t272;
t270 = qJD(2) * t196;
t133 = t252 + t270;
t195 = sin(qJ(5));
t198 = cos(qJ(5));
t74 = t131 * t195 + t133 * t198;
t322 = t74 ^ 2;
t220 = -t198 * t131 + t133 * t195;
t69 = t220 ^ 2;
t339 = t69 - t322;
t338 = t133 * pkin(9) - t281;
t135 = t195 * t196 + t198 * t199;
t213 = t135 * t200;
t80 = t135 * qJD(5) - t195 * t267 - t198 * t266;
t309 = qJD(1) * t213 + t80;
t264 = qJD(5) * t198;
t265 = qJD(5) * t195;
t308 = t342 * t195 + t196 * t264 + t341 * t198 - t199 * t265;
t337 = qJ(6) * t74;
t174 = -qJD(3) + t271;
t165 = qJD(5) + t174;
t262 = qJD(1) * qJD(2);
t243 = t200 * t262;
t248 = t197 * t267;
t261 = qJD(2) * qJD(3);
t93 = qJD(1) * t248 + (-t243 - t261) * t199;
t268 = qJD(2) * t200;
t250 = t196 * t268;
t209 = t197 * t266 + t250;
t94 = qJD(1) * t209 + t196 * t261;
t26 = qJD(5) * t74 - t195 * t93 - t198 * t94;
t330 = t165 * t74 - t26;
t320 = pkin(3) + pkin(4);
t38 = t320 * t174 - t338;
t161 = t174 * qJ(4);
t78 = t196 * t118 + t199 * t156;
t53 = pkin(9) * t131 + t78;
t43 = -t161 + t53;
t12 = t195 * t38 + t198 * t43;
t159 = t174 * qJD(4);
t178 = t197 * t262;
t167 = qJ(4) * t178;
t146 = t224 * qJD(2);
t120 = qJD(1) * t146;
t217 = t118 * t266 + t196 * t120 - t156 * t267;
t227 = pkin(7) * t178;
t205 = t199 * t227 - t217;
t28 = -t159 + t167 - t205;
t20 = pkin(9) * t94 + t28;
t229 = t118 * t267 - t199 * t120 + t156 * t266 - t196 * t227;
t21 = t93 * pkin(9) - t320 * t178 + t229;
t241 = -t195 * t20 + t198 * t21;
t206 = -t12 * qJD(5) + t241;
t155 = -qJD(2) * pkin(2) + pkin(7) * t272;
t66 = t131 * pkin(3) - t133 * qJ(4) + t155;
t45 = -pkin(4) * t131 - t66;
t335 = -t45 * t74 + t206;
t334 = -0.2e1 * t262;
t333 = t74 * t220;
t332 = t340 * t198;
t25 = -t131 * t264 + t133 * t265 - t195 * t94 + t198 * t93;
t331 = -t165 * t220 + t25;
t228 = pkin(3) * t178;
t35 = -t228 + t229;
t60 = -t161 + t78;
t329 = t174 * t60 + t35;
t328 = qJ(6) * t220;
t235 = -qJ(4) * t195 - t198 * t320;
t327 = -qJD(5) * t235 + t195 * t53 + t338 * t198;
t147 = qJ(4) * t198 - t195 * t320;
t326 = -qJD(5) * t147 + t338 * t195 - t198 * t53;
t157 = t319 * t196;
t278 = t195 * t157 + t198 * t158;
t325 = -t157 * t264 + t158 * t265 - t340 * t195 + t336 * t198;
t324 = t342 * qJ(4) + t196 * qJD(4) + t185;
t216 = t195 * t21 + t198 * t20 + t38 * t264 - t43 * t265;
t323 = -t45 * t220 + t216;
t173 = qJ(4) * t287;
t256 = t320 * t196;
t226 = -pkin(7) - t256;
t90 = t197 * t226 + t173;
t188 = t196 * qJ(4);
t249 = t200 * t263;
t277 = qJ(4) * t249 + qJD(4) * t287;
t39 = t197 * qJD(3) * (-t320 * t199 - t188) + t226 * t268 + t277;
t321 = t133 ^ 2;
t11 = -t195 * t43 + t198 * t38;
t7 = t11 - t337;
t6 = pkin(5) * t165 + t7;
t318 = t6 - t7;
t317 = pkin(3) * t174;
t315 = t326 + t328;
t314 = -t308 * qJ(6) - qJD(6) * t135 - t325;
t289 = t196 * t198;
t219 = t195 * t199 - t289;
t313 = pkin(5) * t272 + t309 * qJ(6) - qJD(5) * t278 + t219 * qJD(6) + t336 * t195 + t332;
t176 = pkin(7) * t288;
t191 = t200 * pkin(3);
t67 = t200 * pkin(4) + t176 + t191 + (-pkin(9) * t197 - t149) * t199;
t290 = t196 * t197;
t177 = pkin(7) * t285;
t276 = t196 * t149 + t177;
t91 = -qJ(4) * t200 + t276;
t76 = pkin(9) * t290 + t91;
t310 = t195 * t67 + t198 * t76;
t307 = -qJD(3) * t256 + t320 * t253 + t324;
t306 = t133 * t66;
t33 = t94 * pkin(3) + pkin(7) * t243 + t93 * qJ(4) - t133 * qJD(4);
t302 = t33 * t196;
t301 = t33 * t199;
t300 = t93 * t196;
t299 = -t327 - t337;
t298 = -t341 * pkin(3) - t324;
t296 = t131 * t174;
t295 = t133 * t131;
t294 = t133 * t174;
t293 = t155 * t196;
t292 = t155 * t199;
t291 = t174 * t199;
t203 = qJD(1) ^ 2;
t284 = t200 * t203;
t202 = qJD(2) ^ 2;
t283 = t202 * t197;
t282 = t202 * t200;
t279 = t196 * t146 + t149 * t266;
t84 = t133 * pkin(3) + t131 * qJ(4);
t192 = t197 ^ 2;
t273 = -t200 ^ 2 + t192;
t269 = qJD(2) * t197;
t260 = pkin(8) * t174 * t196;
t259 = pkin(8) * t291;
t148 = -t199 * pkin(3) - pkin(2) - t188;
t258 = pkin(8) * t269;
t257 = pkin(8) * t263;
t255 = qJ(4) * t269 + t279;
t247 = t200 * t267;
t246 = t174 * t266;
t223 = -qJD(3) * t177 + t199 * t146 - t149 * t267;
t32 = pkin(9) * t248 + qJD(2) * t207 - t223;
t34 = (-pkin(7) * qJD(2) + pkin(9) * qJD(3)) * t287 + (-qJD(4) + (-pkin(7) * qJD(3) + pkin(9) * qJD(2)) * t196) * t200 + t255;
t240 = -t195 * t34 + t198 * t32;
t238 = -t195 * t76 + t198 * t67;
t236 = pkin(1) * t334;
t234 = t198 * t157 - t158 * t195;
t233 = t149 * t199 - t176;
t232 = t165 ^ 2;
t125 = t199 * pkin(4) - t148;
t231 = t131 + t263;
t230 = -t133 + t270;
t55 = -pkin(4) * t133 - t84;
t225 = t197 * t254;
t58 = t281 + t317;
t221 = -t196 * t60 + t199 * t58;
t218 = qJD(1) * t192 - t174 * t200;
t215 = t195 * t32 + t198 * t34 + t67 * t264 - t76 * t265;
t214 = -t174 * t78 - t229;
t22 = -t94 * pkin(4) - t33;
t5 = pkin(5) * t26 + t22;
t204 = -t174 * t77 + t205;
t142 = -pkin(5) + t235;
t109 = t135 * t197;
t108 = t195 * t287 - t197 * t289;
t101 = -t173 + (pkin(3) * t196 + pkin(7)) * t197;
t92 = t191 - t233;
t83 = qJD(1) * t225 - t286;
t82 = -pkin(7) * t252 + t280;
t62 = -qJ(6) * t135 + t278;
t61 = qJ(6) * t219 + t234;
t54 = -t93 - t296;
t49 = pkin(3) * t209 + pkin(7) * t268 + qJ(4) * t248 - t277;
t44 = qJD(2) * t225 - t223;
t42 = -t200 * qJD(4) + (-t197 * t263 - t247) * pkin(7) + t255;
t41 = qJD(2) * t213 + (qJD(3) - qJD(5)) * t197 * t219;
t40 = t195 * t249 + t197 * t80 - t198 * t250;
t27 = pkin(5) * t220 + qJD(6) + t45;
t24 = -qJ(6) * t108 + t310;
t23 = pkin(5) * t200 - qJ(6) * t109 + t238;
t8 = t12 - t328;
t4 = -qJ(6) * t40 - qJD(6) * t108 + t215;
t3 = -pkin(5) * t269 - t41 * qJ(6) - qJD(5) * t310 - t109 * qJD(6) + t240;
t2 = -qJ(6) * t26 - qJD(6) * t220 + t216;
t1 = -pkin(5) * t178 + t25 * qJ(6) - t74 * qJD(6) + t206;
t9 = [0, 0, 0, 0.2e1 * t200 * t178, t273 * t334, t282, -t283, 0, -pkin(7) * t282 + t197 * t236, pkin(7) * t283 + t200 * t236, -t93 * t287 + (-t248 + t249) * t133 (-t131 * t199 - t133 * t196) * t268 + (t300 - t199 * t94 + (t131 * t196 - t133 * t199) * qJD(3)) * t197, t174 * t248 + t93 * t200 + (t133 * t197 + t199 * t218) * qJD(2), t197 * t246 + t94 * t200 + (-t131 * t197 - t196 * t218) * qJD(2) (-t174 - t271) * t269, -t223 * t174 + t229 * t200 + (pkin(7) * t94 + t155 * t266) * t197 + ((pkin(7) * t131 + t293) * t200 + (t233 * qJD(1) + t77 + (-t174 + t271) * t316) * t197) * qJD(2) (-pkin(7) * t247 + t279) * t174 + t217 * t200 + (-pkin(7) * t93 - t155 * t267) * t197 + ((pkin(7) * t133 + t292) * t200 + (-pkin(7) * t291 - t276 * qJD(1) - t78) * t197) * qJD(2), t101 * t94 + t49 * t131 + t44 * t174 + (t66 * t270 + t35) * t200 + (t66 * t266 + t302 + (-qJD(1) * t92 - t58) * qJD(2)) * t197, -t131 * t42 + t133 * t44 - t91 * t94 - t92 * t93 + t221 * t268 + (-t196 * t28 + t199 * t35 + (-t196 * t58 - t199 * t60) * qJD(3)) * t197, t101 * t93 - t49 * t133 - t42 * t174 + (-t263 * t66 - t28) * t200 + (t66 * t267 - t301 + (qJD(1) * t91 + t60) * qJD(2)) * t197, t101 * t33 + t28 * t91 + t35 * t92 + t42 * t60 + t44 * t58 + t49 * t66, -t109 * t25 + t41 * t74, t108 * t25 - t109 * t26 - t220 * t41 - t40 * t74, t165 * t41 - t200 * t25 + (-qJD(1) * t109 - t74) * t269, -t165 * t40 - t200 * t26 + (qJD(1) * t108 + t220) * t269 (-t165 - t271) * t269, t240 * t165 + t241 * t200 + t39 * t220 + t90 * t26 + t22 * t108 + t45 * t40 + (-t12 * t200 - t165 * t310) * qJD(5) + (-qJD(1) * t238 - t11) * t269, -t215 * t165 - t216 * t200 + t39 * t74 - t90 * t25 + t22 * t109 + t45 * t41 + (t310 * qJD(1) + t12) * t269, -t1 * t109 - t108 * t2 - t220 * t4 + t23 * t25 - t24 * t26 - t3 * t74 - t40 * t8 - t41 * t6, t1 * t23 + t2 * t24 + t6 * t3 + t8 * t4 + (t108 * pkin(5) + t90) * t5 + (t40 * pkin(5) + t39) * t27; 0, 0, 0, -t197 * t284, t273 * t203, 0, 0, 0, t203 * pkin(1) * t197, pkin(1) * t284, -t133 * t291 - t300 (-t93 + t296) * t199 + (-t94 + t294) * t196, -t246 + (t174 * t285 + t197 * t230) * qJD(1), t174 * t267 + (-t174 * t288 + t197 * t231) * qJD(1), t174 * t272, t174 * t286 - pkin(2) * t94 + (t259 + t293) * qJD(3) + (-t77 * t197 + (-t155 * t200 - t258) * t196 + (t174 * t290 - t200 * t231) * pkin(7)) * qJD(1), pkin(2) * t93 - t114 * t174 + (-t260 + t292) * qJD(3) + (-t155 * t285 + (t78 - t257) * t197 + (t174 * t287 + t200 * t230) * pkin(7)) * qJD(1), t148 * t94 - t83 * t174 - t301 + t298 * t131 + (t196 * t66 + t259) * qJD(3) + (t197 * t58 + (-t200 * t66 - t258) * t196) * qJD(1), t131 * t82 - t133 * t83 + (t28 - t174 * t58 + (qJD(3) * t133 - t94) * pkin(8)) * t199 + ((qJD(3) * t131 - t93) * pkin(8) + t329) * t196, t148 * t93 + t174 * t82 - t302 - t298 * t133 + (-t199 * t66 + t260) * qJD(3) + (t66 * t285 + (-t60 + t257) * t197) * qJD(1), t33 * t148 - t58 * t83 - t60 * t82 + t298 * t66 + (qJD(3) * t221 + t35 * t196 + t28 * t199) * pkin(8), t219 * t25 - t309 * t74, t25 * t135 + t219 * t26 + t220 * t309 - t308 * t74, -t309 * t165 + (qJD(2) * t219 + t74) * t272, -t308 * t165 + (qJD(2) * t135 - t220) * t272, t165 * t272, t125 * t26 + t22 * t135 + t307 * t220 + t308 * t45 + (-t158 * t264 + (-qJD(5) * t157 + t336) * t195 + t332) * t165 + (-qJD(2) * t234 + t11) * t272, -t125 * t25 - t22 * t219 + t307 * t74 - t309 * t45 + t325 * t165 + (qJD(2) * t278 - t12) * t272, t1 * t219 - t135 * t2 - t220 * t314 + t25 * t61 - t26 * t62 - t308 * t8 + t309 * t6 - t313 * t74, t2 * t62 + t1 * t61 + t5 * (pkin(5) * t135 + t125) + t314 * t8 + t313 * t6 + (t308 * pkin(5) + (pkin(4) * t174 + t317) * t196 + t324) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t295, -t131 ^ 2 + t321, t54, -t294 - t94, t178, -t133 * t155 + t214, t131 * t155 + t204, -t131 * t84 + t214 + 0.2e1 * t228 - t306, pkin(3) * t93 - t94 * qJ(4) + (t60 - t78) * t133 + (t58 - t281) * t131, -t131 * t66 + t133 * t84 - 0.2e1 * t159 + 0.2e1 * t167 - t204, -t35 * pkin(3) + t28 * qJ(4) + t281 * t60 - t58 * t78 - t66 * t84, -t333, t339, t331, -t330, t178, t165 * t326 - t235 * t178 - t55 * t220 - t335, t147 * t178 + t327 * t165 - t55 * t74 + t323, t142 * t25 - t147 * t26 + (-t299 + t6) * t220 + (-t315 - t8) * t74, t2 * t147 + t1 * t142 - t27 * (-pkin(5) * t74 + t55) + t299 * t8 + t315 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t178 + t295, t54, -t174 ^ 2 - t321, t306 + t329, 0, 0, 0, 0, 0, -t133 * t220 - t178 * t198 - t195 * t232, -t133 * t74 + t178 * t195 - t198 * t232, t330 * t195 + t331 * t198, -t27 * t133 + (t165 * t8 + t1) * t198 + (-t165 * t6 + t2) * t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t333, -t339, -t331, t330, -t178, t12 * t165 + t335, t11 * t165 - t323, pkin(5) * t25 - t220 * t318, t318 * t8 + (-t27 * t74 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69 - t322, t220 * t8 + t6 * t74 + t5;];
tauc_reg  = t9;
