% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRPRP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:41:31
% EndTime: 2019-03-09 16:41:43
% DurationCPUTime: 4.82s
% Computational Cost: add. (10398->425), mult. (25258->554), div. (0->0), fcn. (18512->8), ass. (0->227)
t233 = sin(qJ(2));
t337 = cos(qJ(3));
t281 = qJD(1) * t337;
t232 = sin(qJ(3));
t235 = cos(qJ(2));
t305 = t232 * t235;
t185 = -qJD(1) * t305 - t233 * t281;
t226 = qJD(2) + qJD(3);
t229 = sin(pkin(10));
t230 = cos(pkin(10));
t162 = t185 * t230 - t226 * t229;
t231 = sin(qJ(5));
t234 = cos(qJ(5));
t171 = t229 * t185;
t344 = t226 * t230 + t171;
t260 = t234 * t344;
t107 = t162 * t231 + t260;
t359 = t107 ^ 2;
t295 = qJD(1) * t233;
t356 = -t232 * t295 + t235 * t281;
t178 = qJD(5) - t356;
t358 = t107 * t178;
t145 = t356 * t226;
t353 = -t162 * t234 + t231 * t344;
t339 = t353 ^ 2;
t304 = t234 * t230;
t317 = t145 * t229;
t46 = t231 * (-qJD(5) * t162 + t317) - qJD(5) * t260 - t145 * t304;
t357 = t46 - t358;
t198 = t229 * t234 + t230 * t231;
t182 = t198 * qJD(5);
t354 = -t198 * t356 + t182;
t197 = t229 * t231 - t304;
t292 = qJD(5) * t234;
t293 = qJD(5) * t231;
t343 = -t229 * t293 + t230 * t292;
t298 = -t197 * t356 - t343;
t291 = qJD(1) * qJD(2);
t352 = -0.2e1 * t291;
t220 = pkin(2) * t232 + qJ(4);
t190 = (-pkin(9) - t220) * t229;
t225 = t230 * pkin(9);
t191 = t220 * t230 + t225;
t142 = t190 * t231 + t191 * t234;
t280 = qJD(3) * t337;
t214 = pkin(2) * t280 + qJD(4);
t311 = t356 * t230;
t273 = -t185 * pkin(4) - pkin(9) * t311;
t147 = -pkin(3) * t185 - qJ(4) * t356;
t130 = pkin(2) * t295 + t147;
t338 = -pkin(8) - pkin(7);
t211 = t338 * t235;
t204 = qJD(1) * t211;
t188 = t232 * t204;
t210 = t338 * t233;
t202 = qJD(1) * t210;
t155 = t337 * t202 + t188;
t89 = t230 * t130 - t155 * t229;
t61 = t273 + t89;
t312 = t356 * t229;
t290 = pkin(9) * t312;
t90 = t229 * t130 + t230 * t155;
t75 = -t290 + t90;
t351 = -qJD(5) * t142 - t198 * t214 + t231 * t75 - t234 * t61;
t200 = t337 * t233 + t305;
t158 = t226 * t200;
t146 = t158 * qJD(1);
t279 = t233 * t291;
t74 = pkin(2) * t279 + pkin(3) * t146 - qJ(4) * t145 + qJD(4) * t185;
t329 = qJD(2) * pkin(2);
t192 = t202 + t329;
t286 = qJD(2) * t338;
t274 = qJD(1) * t286;
t193 = t233 * t274;
t194 = t235 * t274;
t294 = qJD(3) * t232;
t241 = t192 * t280 + t337 * t193 + t232 * t194 + t204 * t294;
t91 = t226 * qJD(4) + t241;
t41 = -t229 * t91 + t230 * t74;
t42 = t229 * t74 + t230 * t91;
t269 = -t41 * t229 + t42 * t230;
t259 = t190 * t234 - t191 * t231;
t350 = -qJD(5) * t259 + t197 * t214 + t231 * t61 + t234 * t75;
t349 = t145 * t198;
t169 = pkin(4) * t312;
t348 = t354 * pkin(5) + t298 * qJ(6) - qJD(6) * t198 - t169;
t207 = (-pkin(9) - qJ(4)) * t229;
t208 = qJ(4) * t230 + t225;
t258 = t207 * t234 - t208 * t231;
t151 = t337 * t192 + t188;
t92 = t230 * t147 - t151 * t229;
t64 = t273 + t92;
t93 = t229 * t147 + t230 * t151;
t77 = -t290 + t93;
t347 = qJD(4) * t197 - qJD(5) * t258 + t231 * t64 + t234 * t77;
t161 = t207 * t231 + t208 * t234;
t346 = -qJD(4) * t198 - qJD(5) * t161 + t231 * t77 - t234 * t64;
t189 = t337 * t204;
t154 = t232 * t202 - t189;
t345 = -pkin(2) * t294 + t154;
t342 = t337 * t210 + t232 * t211;
t341 = qJD(1) * t200;
t252 = -t232 * t233 + t337 * t235;
t224 = -pkin(2) * t235 - pkin(1);
t150 = -pkin(3) * t252 - qJ(4) * t200 + t224;
t166 = t232 * t210 - t337 * t211;
t99 = t230 * t150 - t166 * t229;
t76 = -pkin(4) * t252 - t200 * t225 + t99;
t100 = t229 * t150 + t230 * t166;
t309 = t200 * t229;
t86 = -pkin(9) * t309 + t100;
t263 = t231 * t76 + t234 * t86;
t157 = t226 * t252;
t203 = t233 * t286;
t205 = t235 * t286;
t108 = qJD(3) * t342 + t337 * t203 + t232 * t205;
t288 = t233 * t329;
t87 = pkin(3) * t158 - qJ(4) * t157 - qJD(4) * t200 + t288;
t44 = -t108 * t229 + t230 * t87;
t33 = pkin(4) * t158 - t157 * t225 + t44;
t316 = t157 * t229;
t45 = t230 * t108 + t229 * t87;
t43 = -pkin(9) * t316 + t45;
t340 = -qJD(5) * t263 - t231 * t43 + t234 * t33;
t336 = pkin(5) * t146;
t335 = pkin(5) * t185;
t334 = t230 * pkin(4);
t175 = t185 * qJ(6);
t333 = -t175 + t350;
t332 = t335 + t351;
t129 = -t226 * pkin(3) + qJD(4) - t151;
t96 = -pkin(4) * t344 + t129;
t38 = -pkin(5) * t107 - qJ(6) * t353 + t96;
t328 = t353 * t38;
t326 = -t345 + t348;
t152 = t232 * t192 - t189;
t325 = -t152 + t348;
t324 = t175 - t347;
t323 = -t335 - t346;
t322 = qJ(6) * t146;
t321 = t353 * t107;
t319 = t259 * t146;
t318 = t142 * t146;
t315 = t258 * t146;
t314 = t161 * t146;
t313 = t178 * t185;
t310 = t185 * t356;
t306 = t230 * t145;
t237 = qJD(1) ^ 2;
t303 = t235 * t237;
t236 = qJD(2) ^ 2;
t302 = t236 * t233;
t301 = t236 * t235;
t209 = t224 * qJD(1);
t125 = -pkin(3) * t356 + qJ(4) * t185 + t209;
t133 = qJ(4) * t226 + t152;
t79 = t230 * t125 - t133 * t229;
t52 = -pkin(4) * t356 + pkin(9) * t162 + t79;
t80 = t229 * t125 + t230 * t133;
t59 = pkin(9) * t344 + t80;
t20 = -t231 * t59 + t234 * t52;
t300 = qJD(6) - t20;
t296 = t233 ^ 2 - t235 ^ 2;
t289 = t337 * pkin(2);
t221 = -pkin(3) - t334;
t95 = t192 * t294 + t232 * t193 - t337 * t194 - t204 * t280;
t278 = -t185 * t80 + t95 * t229;
t19 = pkin(4) * t146 - pkin(9) * t306 + t41;
t28 = -pkin(9) * t317 + t42;
t277 = -t234 * t19 + t231 * t28 + t59 * t292 + t52 * t293;
t276 = pkin(1) * t352;
t223 = -t289 - pkin(3);
t272 = -t169 - t345;
t270 = t79 * t185 - t95 * t230;
t268 = -t229 * t79 + t230 * t80;
t21 = t231 * t52 + t234 * t59;
t264 = -t231 * t86 + t234 * t76;
t65 = pkin(4) * t317 + t95;
t262 = t79 * t311 + t80 * t312 + t269;
t261 = t230 * t344;
t127 = pkin(4) * t309 - t342;
t4 = t277 - t336;
t255 = t178 * t21 - t277;
t254 = t231 * t19 + t234 * t28 + t52 * t292 - t59 * t293;
t253 = t231 * t33 + t234 * t43 + t76 * t292 - t86 * t293;
t47 = qJD(5) * t353 + t349;
t10 = pkin(5) * t47 + qJ(6) * t46 - qJD(6) * t353 + t65;
t15 = -pkin(5) * t178 + t300;
t251 = t10 * t197 - t15 * t185 + t354 * t38;
t16 = qJ(6) * t178 + t21;
t250 = -t10 * t198 + t16 * t185 + t298 * t38;
t249 = t20 * t185 + t65 * t197 + t354 * t96;
t248 = -t21 * t185 + t65 * t198 - t298 * t96;
t247 = t185 * t209 - t95;
t109 = t232 * t203 - t337 * t205 + t210 * t294 - t211 * t280;
t246 = t129 * t157 - t145 * t342 + t95 * t200;
t78 = pkin(4) * t316 + t109;
t144 = t197 * pkin(5) - t198 * qJ(6) + t221;
t2 = qJD(6) * t178 + t254 + t322;
t245 = -t15 * t298 - t16 * t354 - t2 * t197 + t4 * t198;
t244 = -pkin(3) * t145 - qJ(4) * t146 - (-qJD(4) + t129) * t356;
t243 = t223 * t145 - t220 * t146 - (t129 - t214) * t356;
t240 = -t209 * t356 - t241;
t239 = t226 * t341;
t238 = t178 * t353 + t47;
t206 = t223 - t334;
t137 = t197 * t200;
t136 = t198 * t200;
t131 = -t289 + t144;
t116 = t185 ^ 2 - t356 ^ 2;
t113 = t169 + t152;
t112 = -t185 * t226 - t239;
t63 = t157 * t198 + t200 * t343;
t62 = t157 * t197 + t182 * t200;
t57 = t136 * pkin(5) + t137 * qJ(6) + t127;
t55 = pkin(5) * t353 - qJ(6) * t107;
t37 = pkin(5) * t252 - t264;
t36 = -qJ(6) * t252 + t263;
t31 = t107 * t185 - t146 * t197 - t178 * t354;
t30 = t146 * t198 - t298 * t178 + t185 * t353;
t29 = -t46 - t358;
t12 = -t198 * t46 - t298 * t353;
t11 = t63 * pkin(5) + t62 * qJ(6) + t137 * qJD(6) + t78;
t7 = -pkin(5) * t158 - t340;
t6 = qJ(6) * t158 - qJD(6) * t252 + t253;
t5 = -t107 * t298 + t197 * t46 - t198 * t47 - t353 * t354;
t1 = [0, 0, 0, 0.2e1 * t235 * t279, t296 * t352, t301, -t302, 0, -pkin(7) * t301 + t233 * t276, pkin(7) * t302 + t235 * t276, t145 * t200 - t157 * t185, t145 * t252 - t146 * t200 + t157 * t356 + t158 * t185, t157 * t226, -t158 * t226, 0, -t109 * t226 + t146 * t224 + t158 * t209 + (-qJD(1) * t252 - t356) * t288, -t108 * t226 + t145 * t224 + t157 * t209 + (-t185 + t341) * t288, -t109 * t344 + t99 * t146 + t79 * t158 + t229 * t246 - t252 * t41 - t356 * t44, -t100 * t146 - t109 * t162 - t158 * t80 + t230 * t246 + t252 * t42 + t356 * t45, t44 * t162 + t45 * t171 + (-t99 * t145 - t79 * t157 - t41 * t200 + t45 * t226) * t230 + (-t100 * t145 - t80 * t157 - t42 * t200) * t229, t100 * t42 + t109 * t129 - t342 * t95 + t41 * t99 + t44 * t79 + t45 * t80, t137 * t46 - t353 * t62, -t107 * t62 + t136 * t46 + t137 * t47 - t353 * t63, -t137 * t146 + t158 * t353 - t178 * t62 + t252 * t46, t107 * t158 - t136 * t146 - t178 * t63 + t252 * t47, -t146 * t252 + t158 * t178, -t107 * t78 + t127 * t47 + t65 * t136 + t264 * t146 + t20 * t158 + t178 * t340 + t252 * t277 + t96 * t63, -t127 * t46 - t65 * t137 - t263 * t146 - t21 * t158 - t253 * t178 + t252 * t254 + t353 * t78 - t96 * t62, t10 * t136 - t107 * t11 - t146 * t37 - t15 * t158 - t178 * t7 + t252 * t4 + t38 * t63 + t47 * t57, t107 * t6 - t136 * t2 - t137 * t4 - t15 * t62 - t16 * t63 + t353 * t7 - t36 * t47 - t37 * t46, t10 * t137 - t11 * t353 + t146 * t36 + t158 * t16 + t178 * t6 - t2 * t252 + t38 * t62 + t46 * t57, t10 * t57 + t11 * t38 + t15 * t7 + t16 * t6 + t2 * t36 + t37 * t4; 0, 0, 0, -t233 * t303, t296 * t237, 0, 0, 0, t237 * pkin(1) * t233, pkin(1) * t303, t310, t116, 0, t112, 0, t154 * t226 + (-t226 * t294 + t295 * t356) * pkin(2) + t247, t155 * t226 + (t185 * t295 - t226 * t280) * pkin(2) + t240, t229 * t243 + t344 * t345 + t356 * t89 + t270, t162 * t345 + t230 * t243 - t356 * t90 + t278, t214 * t261 - t90 * t344 - (t229 * t214 + t89) * t162 + t262, -t129 * t345 + t214 * t268 + t220 * t269 + t223 * t95 - t79 * t89 - t80 * t90, t12, t5, t30, t31, t313, -t107 * t272 + t178 * t351 + t206 * t47 + t249 + t319, t178 * t350 - t206 * t46 + t272 * t353 + t248 - t318, -t107 * t326 + t131 * t47 + t332 * t178 + t251 + t319, -t107 * t333 - t142 * t47 + t259 * t46 - t332 * t353 + t245, t131 * t46 - t333 * t178 - t326 * t353 + t250 + t318, t10 * t131 + t142 * t2 - t332 * t15 - t333 * t16 - t259 * t4 + t326 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t310, t116, 0, t112, 0, t152 * t226 + t247, t151 * t226 + t240, t152 * t344 + t229 * t244 + t356 * t92 + t270, t152 * t162 + t230 * t244 - t356 * t93 + t278, -t93 * t344 - t92 * t162 + (-t162 * t229 + t261) * qJD(4) + t262, -pkin(3) * t95 + qJ(4) * t269 + qJD(4) * t268 - t129 * t152 - t79 * t92 - t80 * t93, t12, t5, t30, t31, t313, t107 * t113 + t178 * t346 + t221 * t47 + t249 + t315, -t113 * t353 + t178 * t347 - t221 * t46 + t248 - t314, -t107 * t325 + t144 * t47 - t178 * t323 + t251 + t315, t107 * t324 - t161 * t47 + t258 * t46 + t323 * t353 + t245, t144 * t46 + t178 * t324 - t325 * t353 + t250 + t314, t10 * t144 + t15 * t323 + t16 * t324 + t161 * t2 - t258 * t4 + t325 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t162 * t356 + t317, -t344 * t356 + t306, -t162 ^ 2 - t344 ^ 2, -t162 * t79 - t344 * t80 + t95, 0, 0, 0, 0, 0, t238, -t357, t238, -t339 - t359, t357, -t107 * t16 - t15 * t353 + t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t321, t339 - t359, t29, -t349 + (-qJD(5) + t178) * t353, t146, -t353 * t96 + t255, -t107 * t96 + t178 * t20 - t254, t107 * t55 + t255 - t328 + 0.2e1 * t336, pkin(5) * t46 - qJ(6) * t47 + (t16 - t21) * t353 - (t15 - t300) * t107, 0.2e1 * t322 + t107 * t38 + t353 * t55 + (0.2e1 * qJD(6) - t20) * t178 + t254, -pkin(5) * t4 + qJ(6) * t2 - t15 * t21 + t16 * t300 - t38 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t239 - t321, t29, -t178 ^ 2 - t339, -t16 * t178 + t328 + t4;];
tauc_reg  = t1;
