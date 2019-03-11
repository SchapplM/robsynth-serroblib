% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPRP5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP5_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:16:23
% EndTime: 2019-03-09 03:16:36
% DurationCPUTime: 5.35s
% Computational Cost: add. (8323->474), mult. (19730->583), div. (0->0), fcn. (15443->14), ass. (0->226)
t195 = sin(pkin(10));
t197 = cos(pkin(10));
t196 = sin(pkin(9));
t198 = cos(pkin(9));
t202 = sin(qJ(3));
t326 = cos(qJ(3));
t161 = t326 * t196 + t202 * t198;
t332 = t161 * qJD(1);
t124 = qJD(3) * t195 + t197 * t332;
t201 = sin(qJ(5));
t325 = cos(qJ(5));
t137 = t195 * t332;
t281 = qJD(3) * t197;
t339 = -t137 + t281;
t226 = t325 * t339;
t74 = -t124 * t201 + t226;
t357 = t74 ^ 2;
t268 = t326 * t198;
t229 = -t202 * t196 + t268;
t350 = t229 * qJD(1);
t141 = qJD(5) - t350;
t356 = t141 * t74;
t243 = t201 * t339;
t352 = t325 * t124 + t243;
t327 = t352 ^ 2;
t194 = pkin(9) + qJ(3);
t189 = cos(t194);
t180 = g(3) * t189;
t187 = sin(t194);
t203 = sin(qJ(1));
t204 = cos(qJ(1));
t249 = g(1) * t204 + g(2) * t203;
t217 = t249 * t187 - t180;
t276 = qJD(1) * qJD(2);
t318 = pkin(7) + qJ(2);
t328 = qJDD(1) * t318 + t276;
t135 = t328 * t196;
t136 = t328 * t198;
t168 = t318 * t196;
t162 = qJD(1) * t168;
t170 = t318 * t198;
t163 = qJD(1) * t170;
t265 = qJD(3) * t326;
t280 = qJD(3) * t202;
t225 = -t326 * t135 - t202 * t136 + t162 * t280 - t163 * t265;
t302 = qJDD(3) * pkin(3);
t61 = qJDD(4) - t225 - t302;
t215 = -t61 + t217;
t355 = t141 * t352;
t152 = t161 * qJD(3);
t275 = t196 * qJDD(1);
t246 = -qJDD(1) * t268 + t202 * t275;
t111 = qJD(1) * t152 + t246;
t108 = qJDD(5) + t111;
t160 = t325 * t195 + t201 * t197;
t227 = -t201 * t195 + t325 * t197;
t264 = qJD(5) * t325;
t278 = qJD(5) * t201;
t335 = -t195 * t278 + t197 * t264;
t349 = -t227 * t350 + t335;
t253 = t160 * t108 + t349 * t141;
t310 = t332 * t352;
t354 = t253 + t310;
t317 = pkin(8) + qJ(4);
t167 = t317 * t195;
t169 = t317 * t197;
t121 = -t201 * t167 + t325 * t169;
t193 = pkin(10) + qJ(5);
t186 = sin(t193);
t353 = t121 * t108 + t186 * t217;
t150 = t160 * qJD(5);
t305 = -t160 * t350 + t150;
t252 = t227 * t108 - t141 * t305;
t311 = t332 * t74;
t351 = t252 + t311;
t348 = t124 * t350;
t303 = qJDD(1) * pkin(1);
t185 = qJDD(2) - t303;
t337 = g(1) * t203 - g(2) * t204;
t242 = -t185 + t337;
t345 = qJD(3) * t332;
t292 = t187 * t204;
t293 = t187 * t203;
t344 = g(1) * t292 + g(2) * t293 - t180;
t107 = pkin(3) * t332 - qJ(4) * t350;
t230 = t326 * t162 + t202 * t163;
t63 = t195 * t107 - t197 * t230;
t343 = qJD(4) * t197 - t63;
t228 = -t325 * t167 - t201 * t169;
t298 = t350 * t197;
t62 = t197 * t107 + t195 * t230;
t42 = pkin(4) * t332 - pkin(8) * t298 + t62;
t299 = t350 * t195;
t52 = -pkin(8) * t299 + t63;
t342 = -qJD(4) * t227 - qJD(5) * t228 + t201 * t42 + t325 * t52;
t341 = -qJD(4) * t160 - qJD(5) * t121 + t201 * t52 - t325 * t42;
t338 = 0.2e1 * t350;
t334 = -t326 * t168 - t202 * t170;
t151 = t229 * qJD(3);
t221 = t161 * qJDD(1);
t211 = qJD(1) * t151 + t221;
t274 = t198 * qJDD(1);
t50 = -pkin(2) * t274 + t111 * pkin(3) - qJ(4) * t211 - qJD(4) * t332 + t185;
t231 = -t202 * t135 + t326 * t136;
t58 = qJDD(3) * qJ(4) + (qJD(4) - t230) * qJD(3) + t231;
t22 = -t195 * t58 + t197 * t50;
t23 = t195 * t50 + t197 * t58;
t333 = -t22 * t195 + t23 * t197;
t331 = qJ(2) * qJDD(1);
t210 = t195 * qJDD(3) + t197 * t211;
t103 = t195 * t211;
t258 = qJDD(3) * t197 - t103;
t26 = -qJD(5) * t226 + t124 * t278 - t201 * t258 - t325 * t210;
t330 = -t227 * t26 - t305 * t352;
t320 = g(3) * t187;
t216 = -t249 * t189 - t320;
t294 = t161 * t197;
t184 = t198 * pkin(2) + pkin(1);
t110 = -pkin(3) * t229 - qJ(4) * t161 - t184;
t122 = -t202 * t168 + t326 * t170;
t65 = t197 * t110 - t122 * t195;
t49 = -pkin(4) * t229 - pkin(8) * t294 + t65;
t295 = t161 * t195;
t66 = t195 * t110 + t197 * t122;
t57 = -pkin(8) * t295 + t66;
t237 = t201 * t49 + t325 * t57;
t296 = t151 * t197;
t79 = pkin(3) * t152 - qJ(4) * t151 - qJD(4) * t161;
t91 = t229 * qJD(2) + qJD(3) * t334;
t47 = -t195 * t91 + t197 * t79;
t30 = pkin(4) * t152 - pkin(8) * t296 + t47;
t297 = t151 * t195;
t48 = t195 * t79 + t197 * t91;
t37 = -pkin(8) * t297 + t48;
t329 = -qJD(5) * t237 - t201 * t37 + t325 * t30;
t142 = t350 ^ 2;
t324 = pkin(5) * t108;
t319 = t352 * t74;
t115 = -t202 * t162 + t326 * t163;
t77 = pkin(4) * t299 + t115;
t315 = t305 * pkin(5) - t349 * qJ(6) - qJD(6) * t160 - t77;
t314 = -qJ(6) * t332 - t342;
t313 = pkin(5) * t332 - t341;
t106 = qJD(3) * qJ(4) + t115;
t166 = -qJD(1) * t184 + qJD(2);
t88 = -pkin(3) * t350 - qJ(4) * t332 + t166;
t56 = t197 * t106 + t195 * t88;
t55 = -t106 * t195 + t197 * t88;
t33 = -pkin(4) * t350 - pkin(8) * t124 + t55;
t41 = pkin(8) * t339 + t56;
t13 = t201 * t33 + t325 * t41;
t312 = t13 * t141;
t304 = qJ(6) * t108;
t301 = t111 * t195;
t291 = t189 * t203;
t290 = t189 * t204;
t289 = t197 * t111;
t188 = cos(t193);
t285 = t203 * t188;
t284 = t204 * t186;
t12 = -t201 * t41 + t325 * t33;
t283 = qJD(6) - t12;
t282 = t196 ^ 2 + t198 ^ 2;
t277 = t115 * qJD(3);
t183 = t197 * pkin(4) + pkin(3);
t10 = t111 * pkin(4) - pkin(8) * t210 + t22;
t16 = pkin(8) * t258 + t23;
t263 = t325 * t10 - t201 * t16 - t41 * t264 - t33 * t278;
t262 = pkin(4) * t195 + t318;
t260 = t282 * qJD(1) ^ 2;
t92 = t161 * qJD(2) - t168 * t280 + t170 * t265;
t209 = t201 * t210 - t325 * t258;
t27 = qJD(5) * t352 + t209;
t257 = -t160 * t27 + t349 * t74;
t256 = 0.2e1 * t282;
t254 = -g(1) * t293 + g(2) * t292;
t128 = t186 * t291 + t188 * t204;
t130 = t189 * t284 - t285;
t251 = -g(1) * t128 + g(2) * t130;
t129 = t189 * t285 - t284;
t131 = t186 * t203 + t188 * t290;
t250 = g(1) * t129 - g(2) * t131;
t247 = pkin(3) * t189 + qJ(4) * t187;
t245 = -t55 * t195 + t56 * t197;
t67 = pkin(4) * t297 + t92;
t244 = t183 * t189 + t187 * t317;
t93 = pkin(4) * t295 - t334;
t241 = pkin(5) * t188 + qJ(6) * t186 + t183;
t239 = -t201 * t57 + t325 * t49;
t236 = t201 * t10 + t325 * t16 + t33 * t264 - t278 * t41;
t233 = t337 * t189;
t232 = t201 * t30 + t49 * t264 - t278 * t57 + t325 * t37;
t224 = t228 * t108 + t344 * t188;
t220 = t242 + t303;
t104 = -qJD(3) * pkin(3) + qJD(4) + t230;
t219 = t26 - t356;
t218 = g(1) * t130 + g(2) * t128 + t186 * t320 + t263;
t214 = t256 * t276 - t249;
t68 = -pkin(4) * t339 + t104;
t24 = -pkin(5) * t74 - qJ(6) * t352 + t68;
t213 = t24 * t352 + qJDD(6) - t218;
t212 = -g(1) * t131 - g(2) * t129 - t188 * t320 + t236;
t34 = -pkin(4) * t258 + t61;
t4 = t27 * pkin(5) + t26 * qJ(6) - qJD(6) * t352 + t34;
t208 = -qJD(5) * t243 - t124 * t264 - t209;
t207 = -t208 + t355;
t171 = t204 * t184;
t164 = -qJDD(1) * t184 + qJDD(2);
t109 = -pkin(5) * t227 - qJ(6) * t160 - t183;
t100 = t227 * t161;
t99 = t160 * t161;
t60 = t151 * t160 + t161 * t335;
t59 = t150 * t161 - t151 * t227;
t39 = pkin(5) * t352 - qJ(6) * t74;
t36 = t99 * pkin(5) - t100 * qJ(6) + t93;
t21 = pkin(5) * t229 - t239;
t20 = -qJ(6) * t229 + t237;
t19 = -t26 - t356;
t11 = t141 * qJ(6) + t13;
t9 = -t141 * pkin(5) + t283;
t6 = pkin(5) * t60 + qJ(6) * t59 - qJD(6) * t100 + t67;
t5 = -t152 * pkin(5) - t329;
t3 = qJ(6) * t152 - qJD(6) * t229 + t232;
t2 = qJDD(6) - t263 - t324;
t1 = qJD(6) * t141 + t236 + t304;
t7 = [qJDD(1), t337, t249, t220 * t198, -t220 * t196, t256 * t331 + t214, pkin(1) * t242 + (t282 * t331 + t214) * qJ(2), t151 * t332 + t161 * t211, -t161 * t111 + t151 * t350 - t152 * t332 + t211 * t229, qJD(3) * t151 + qJDD(3) * t161, -qJD(3) * t152 + qJDD(3) * t229, 0, -qJD(3) * t92 + qJDD(3) * t334 - t111 * t184 + t152 * t166 - t164 * t229 + t233, -t91 * qJD(3) - t122 * qJDD(3) + t166 * t151 + t164 * t161 - t184 * t211 + t254, -t47 * t350 + t65 * t111 - t22 * t229 + t55 * t152 - t92 * t339 + t334 * t258 + t197 * t233 + (t104 * t151 + t61 * t161 - t249) * t195, t48 * t350 - t66 * t111 + t23 * t229 - t56 * t152 + t92 * t124 - t334 * t210 + t61 * t294 + t104 * t296 - g(1) * (t195 * t291 + t197 * t204) - g(2) * (-t195 * t290 + t197 * t203) -t47 * t124 - t210 * t65 - t22 * t294 - t23 * t295 + t258 * t66 - t296 * t55 - t297 * t56 + t339 * t48 - t254, -g(2) * t171 + t104 * t92 - t61 * t334 + t22 * t65 + t23 * t66 + t55 * t47 + t56 * t48 + (-g(1) * t318 - g(2) * t247) * t204 + (-g(1) * (-t184 - t247) - g(2) * t318) * t203, -t100 * t26 - t352 * t59, -t100 * t27 + t26 * t99 - t352 * t60 - t59 * t74, t100 * t108 - t141 * t59 + t152 * t352 + t229 * t26, -t108 * t99 - t141 * t60 + t152 * t74 + t229 * t27, -t108 * t229 + t141 * t152, t239 * t108 + t12 * t152 + t141 * t329 - t229 * t263 + t93 * t27 + t34 * t99 + t68 * t60 - t67 * t74 + t250, t34 * t100 - t108 * t237 - t13 * t152 - t141 * t232 + t229 * t236 - t93 * t26 + t352 * t67 - t68 * t59 + t251, -t108 * t21 - t141 * t5 - t152 * t9 + t2 * t229 + t24 * t60 + t27 * t36 + t4 * t99 - t6 * t74 + t250, -t1 * t99 + t100 * t2 - t11 * t60 - t20 * t27 - t21 * t26 + t3 * t74 + t352 * t5 - t59 * t9 - t254, -t1 * t229 - t100 * t4 + t108 * t20 + t11 * t152 + t141 * t3 + t24 * t59 + t26 * t36 - t352 * t6 - t251, t1 * t20 + t11 * t3 + t4 * t36 + t24 * t6 + t2 * t21 + t9 * t5 - g(1) * (-pkin(5) * t129 - qJ(6) * t128) - g(2) * (pkin(5) * t131 + qJ(6) * t130 + t171) + (-g(1) * t262 - g(2) * t244) * t204 + (-g(1) * (-t184 - t244) - g(2) * t262) * t203; 0, 0, 0, -t274, t275, -t260, -qJ(2) * t260 - t242, 0, 0, 0, 0, 0, t246 + 0.2e1 * t345, qJD(3) * t338 + t221, -t195 * t142 + t332 * t339 + t289, -t124 * t332 - t142 * t197 - t301 (-qJDD(1) * t294 + t137 * t350 - t281 * t338) * t197 + (-t103 - t348) * t195, -t104 * t332 + t23 * t195 + t22 * t197 - t245 * t350 - t337, 0, 0, 0, 0, 0, t351, -t354, t351, t257 - t330, t354, t1 * t160 + t11 * t349 - t2 * t227 - t24 * t332 + t305 * t9 - t337; 0, 0, 0, 0, 0, 0, 0, -t332 * t350, t332 ^ 2 - t142, t221, -t246, qJDD(3), -t166 * t332 + t217 + t225 + t277, -t166 * t350 - t216 - t231, -qJ(4) * t301 - pkin(3) * t103 - t115 * t137 - t55 * t332 - (-t62 + (-qJD(4) + t104) * t195) * t350 + (t215 + t277 + t302) * t197, -pkin(3) * t210 - qJ(4) * t289 - t104 * t298 - t115 * t124 + t56 * t332 + t343 * t350 + (t61 - t344) * t195, -g(1) * t290 - g(2) * t291 + t298 * t55 + t299 * t56 - t320 + t343 * t339 + (qJD(4) * t195 + t62) * t124 + (t195 * t210 + t197 * t258) * qJ(4) + t333, -t104 * t115 - t55 * t62 - t56 * t63 + t245 * qJD(4) + t215 * pkin(3) + (t216 + t333) * qJ(4), -t160 * t26 + t349 * t352, t257 + t330, t253 - t310, t252 - t311, -t141 * t332, -t12 * t332 + t341 * t141 - t183 * t27 - t227 * t34 + t305 * t68 + t74 * t77 + t224, t13 * t332 + t342 * t141 + t34 * t160 + t183 * t26 + t349 * t68 - t352 * t77 - t353, t109 * t27 - t141 * t313 - t227 * t4 + t24 * t305 - t315 * t74 + t332 * t9 + t224, t1 * t227 - t11 * t305 - t121 * t27 + t160 * t2 + t228 * t26 + t313 * t352 + t314 * t74 + t349 * t9 + t216, t109 * t26 - t11 * t332 + t141 * t314 - t160 * t4 - t24 * t349 - t315 * t352 + t353, t1 * t121 + t4 * t109 - t2 * t228 + t313 * t9 + t315 * t24 + t314 * t11 + (-g(3) * t241 - t249 * t317) * t189 + (-g(3) * t317 + t241 * t249) * t187; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t258 - t348, -t339 * t350 + t210, -t124 ^ 2 - t339 ^ 2, t124 * t55 - t339 * t56 - t215, 0, 0, 0, 0, 0, t207, -t219, t207, -t327 - t357, t219, -t11 * t74 - t352 * t9 - t217 + t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t319, t327 - t357, t19, t208 + t355, t108, -t352 * t68 + t218 + t312, t12 * t141 - t68 * t74 - t212, t39 * t74 - t213 + t312 + 0.2e1 * t324, pkin(5) * t26 - qJ(6) * t27 + (t11 - t13) * t352 - (t9 - t283) * t74, 0.2e1 * t304 + t24 * t74 + t39 * t352 + (0.2e1 * qJD(6) - t12) * t141 + t212, t1 * qJ(6) - t2 * pkin(5) - t24 * t39 - t9 * t13 - g(1) * (-pkin(5) * t130 + qJ(6) * t131) - g(2) * (-pkin(5) * t128 + qJ(6) * t129) - (-pkin(5) * t186 + qJ(6) * t188) * t320 + t283 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(5) - t246 - t319 - t345, t19, -t141 ^ 2 - t327, -t11 * t141 + t213 - t324;];
tau_reg  = t7;
