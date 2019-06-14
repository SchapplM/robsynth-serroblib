% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRRRRP3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 09:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRRRRP3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP3_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:47:40
% EndTime: 2019-05-05 09:47:55
% DurationCPUTime: 5.27s
% Computational Cost: add. (22493->405), mult. (44259->542), div. (0->0), fcn. (32639->12), ass. (0->263)
t240 = sin(qJ(4));
t244 = cos(qJ(4));
t241 = sin(qJ(3));
t290 = qJD(2) * t241;
t205 = qJD(3) * t240 + t244 * t290;
t239 = sin(qJ(5));
t243 = cos(qJ(5));
t270 = qJD(3) * t244 - t240 * t290;
t186 = t239 * t205 - t243 * t270;
t188 = t243 * t205 + t239 * t270;
t150 = t188 * t186;
t285 = qJD(2) * qJD(3);
t228 = t241 * t285;
t245 = cos(qJ(3));
t283 = t245 * qJDD(2);
t209 = -t228 + t283;
t204 = -qJDD(4) + t209;
t201 = -qJDD(5) + t204;
t336 = t150 + t201;
t342 = pkin(5) * t336;
t275 = t245 * t285;
t284 = t241 * qJDD(2);
t208 = t275 + t284;
t180 = qJD(4) * t270 + t240 * qJDD(3) + t244 * t208;
t259 = t244 * qJDD(3) - t240 * t208;
t254 = -qJD(4) * t205 + t259;
t125 = -t186 * qJD(5) + t243 * t180 + t239 * t254;
t225 = qJD(2) * t245 - qJD(4);
t218 = -qJD(5) + t225;
t166 = t186 * t218;
t106 = -t166 + t125;
t337 = qJ(6) * t106;
t197 = t270 * t225;
t155 = t197 + t180;
t261 = t270 * t205;
t335 = -t204 + t261;
t235 = sin(pkin(6));
t236 = cos(pkin(6));
t307 = sin(pkin(11));
t308 = cos(pkin(11));
t255 = t307 * g(1) - t308 * g(2);
t291 = -g(3) + qJDD(1);
t249 = -t235 * t255 + t236 * t291;
t248 = t241 * t249;
t213 = -t308 * g(1) - t307 * g(2);
t242 = sin(qJ(2));
t246 = cos(qJ(2));
t252 = t236 * t255;
t334 = t235 * t291 + t252;
t173 = t246 * t213 + t242 * t334;
t247 = qJD(2) ^ 2;
t161 = -t247 * pkin(2) + qJDD(2) * pkin(8) + t173;
t267 = -pkin(3) * t245 - pkin(9) * t241;
t271 = t247 * t267 + t161;
t330 = qJD(3) ^ 2;
t130 = -t330 * pkin(3) + qJDD(3) * pkin(9) + t245 * t271 + t248;
t266 = t242 * t213 - t246 * t334;
t160 = -qJDD(2) * pkin(2) - t247 * pkin(8) + t266;
t263 = -t209 + t228;
t264 = t208 + t275;
t135 = pkin(3) * t263 - pkin(9) * t264 + t160;
t86 = t240 * t130 - t244 * t135;
t74 = pkin(4) * t335 - t155 * pkin(10) - t86;
t194 = -pkin(4) * t225 - pkin(10) * t205;
t268 = t270 ^ 2;
t87 = t244 * t130 + t240 * t135;
t80 = -pkin(4) * t268 + pkin(10) * t254 + t225 * t194 + t87;
t41 = t239 * t80 - t243 * t74;
t253 = 0.2e1 * qJD(6) * t188 + t337 + t342 + t41;
t251 = -t253 - t342;
t184 = t186 ^ 2;
t217 = t218 ^ 2;
t143 = -t217 - t184;
t301 = t336 * t243;
t93 = t143 * t239 - t301;
t92 = pkin(4) * t93;
t341 = t251 + t92;
t340 = t240 * t335;
t339 = t244 * t335;
t185 = t188 ^ 2;
t159 = -t185 - t217;
t137 = -t150 + t201;
t304 = t137 * t239;
t112 = t159 * t243 + t304;
t111 = pkin(4) * t112;
t272 = t239 * t180 - t243 * t254;
t124 = -qJD(5) * t188 - t272;
t162 = -pkin(5) * t218 - qJ(6) * t188;
t42 = t239 * t74 + t243 * t80;
t30 = -t184 * pkin(5) + t124 * qJ(6) - 0.2e1 * qJD(6) * t186 + t218 * t162 + t42;
t257 = pkin(5) * t159 - t30;
t338 = t111 + t257;
t193 = t245 * t249;
t129 = -qJDD(3) * pkin(3) - t330 * pkin(9) + t271 * t241 - t193;
t89 = -t254 * pkin(4) - t268 * pkin(10) + t205 * t194 + t129;
t302 = t336 * t239;
t53 = -t124 * pkin(5) - t184 * qJ(6) + t188 * t162 + qJDD(6) + t89;
t332 = t166 + t125;
t156 = t197 - t180;
t103 = (qJD(5) + t218) * t188 + t272;
t151 = (qJD(4) + t225) * t205 - t259;
t203 = t205 ^ 2;
t222 = t225 ^ 2;
t94 = t143 * t243 + t302;
t61 = t240 * t94 + t244 * t93;
t329 = pkin(3) * t61;
t303 = t137 * t243;
t113 = -t159 * t239 + t303;
t77 = t112 * t244 + t113 * t240;
t328 = pkin(3) * t77;
t20 = t239 * t42 - t243 * t41;
t327 = pkin(4) * t20;
t68 = -t103 * t239 - t106 * t243;
t70 = -t103 * t243 + t106 * t239;
t37 = t240 * t70 + t244 * t68;
t326 = pkin(9) * t37;
t325 = pkin(9) * t61;
t324 = pkin(9) * t77;
t323 = pkin(10) * t68;
t322 = pkin(10) * t93;
t321 = pkin(10) * t112;
t127 = -t184 - t185;
t39 = -t240 * t68 + t244 * t70;
t32 = t127 * t241 + t245 * t39;
t320 = -pkin(2) * t37 + pkin(8) * t32;
t102 = (qJD(5) - t218) * t188 + t272;
t62 = -t240 * t93 + t244 * t94;
t45 = t102 * t241 + t245 * t62;
t319 = -pkin(2) * t61 + pkin(8) * t45;
t78 = -t112 * t240 + t113 * t244;
t50 = t241 * t332 + t245 * t78;
t318 = -pkin(2) * t77 + pkin(8) * t50;
t317 = t20 * t240;
t316 = t20 * t244;
t315 = t239 * t253;
t314 = t239 * t89;
t313 = t243 * t253;
t312 = t243 * t89;
t311 = -pkin(3) * t127 + pkin(9) * t39;
t310 = -pkin(3) * t102 + pkin(9) * t62;
t309 = -pkin(3) * t332 + pkin(9) * t78;
t306 = t129 * t240;
t305 = t129 * t244;
t169 = t204 + t261;
t300 = t169 * t240;
t299 = t169 * t244;
t297 = t218 * t239;
t296 = t218 * t243;
t295 = t240 * t205;
t224 = t241 * t247 * t245;
t215 = qJDD(3) + t224;
t294 = t241 * t215;
t293 = t244 * t205;
t214 = -t224 + qJDD(3);
t292 = t245 * t214;
t282 = t111 - t42;
t281 = t245 * t150;
t67 = pkin(4) * t68;
t279 = -pkin(3) * t37 - t67;
t12 = t239 * t30 - t313;
t28 = pkin(5) * t253;
t278 = pkin(4) * t12 - t28;
t277 = -pkin(4) * t102 + pkin(10) * t94;
t276 = -pkin(4) * t127 + pkin(10) * t70;
t274 = -pkin(4) * t332 + pkin(10) * t113;
t21 = t239 * t41 + t243 * t42;
t57 = t240 * t86 + t244 * t87;
t146 = t161 * t241 - t193;
t147 = t245 * t161 + t248;
t108 = t146 * t241 + t245 * t147;
t269 = t41 - t92;
t56 = t240 * t87 - t244 * t86;
t260 = -pkin(2) + t267;
t258 = t245 * t261;
t232 = t245 ^ 2;
t231 = t241 ^ 2;
t230 = t232 * t247;
t229 = t231 * t247;
t221 = -t230 - t330;
t220 = -t229 - t330;
t212 = t229 + t230;
t211 = (t231 + t232) * qJDD(2);
t210 = -0.2e1 * t228 + t283;
t207 = 0.2e1 * t275 + t284;
t196 = -t203 + t222;
t195 = t268 - t222;
t192 = -t220 * t241 - t292;
t191 = t221 * t245 - t294;
t190 = t203 - t268;
t189 = -t203 - t222;
t182 = -t222 - t268;
t168 = t268 + t203;
t164 = -t185 + t217;
t163 = t184 - t217;
t152 = (-qJD(4) + t225) * t205 + t259;
t148 = t185 - t184;
t145 = -t189 * t240 + t299;
t144 = t189 * t244 + t300;
t141 = t182 * t244 - t340;
t140 = t182 * t240 + t339;
t132 = (t186 * t243 - t188 * t239) * t218;
t131 = (t186 * t239 + t188 * t243) * t218;
t120 = -t151 * t244 + t155 * t240;
t119 = -t151 * t240 - t155 * t244;
t118 = t163 * t243 + t304;
t117 = -t164 * t239 - t301;
t116 = t163 * t239 - t303;
t115 = t164 * t243 - t302;
t114 = t145 * t245 - t156 * t241;
t109 = t141 * t245 - t152 * t241;
t99 = pkin(5) * t106;
t98 = t125 * t243 + t188 * t297;
t97 = t125 * t239 - t188 * t296;
t96 = -t124 * t239 - t186 * t296;
t95 = t124 * t243 - t186 * t297;
t90 = t120 * t245 - t168 * t241;
t88 = t131 * t244 + t132 * t240;
t84 = t241 * (-t131 * t240 + t132 * t244) + t245 * t201;
t83 = -pkin(5) * t332 + qJ(6) * t137;
t82 = t116 * t244 + t118 * t240;
t81 = t115 * t244 + t117 * t240;
t71 = -t102 * t243 - t239 * t332;
t69 = -t102 * t239 + t243 * t332;
t65 = t312 - t321;
t64 = t240 * t98 + t244 * t97;
t63 = t240 * t96 + t244 * t95;
t58 = t314 - t322;
t55 = t241 * (-t240 * t97 + t244 * t98) - t281;
t54 = t241 * (-t240 * t95 + t244 * t96) + t281;
t52 = t241 * (-t116 * t240 + t118 * t244) + t245 * t103;
t51 = t241 * (-t115 * t240 + t117 * t244) - t245 * t106;
t48 = -qJ(6) * t159 + t53;
t47 = t274 + t314;
t46 = t129 * t241 + t245 * t57;
t43 = t277 - t312;
t38 = t240 * t71 + t244 * t69;
t34 = -pkin(5) * t102 + qJ(6) * t143 - t53;
t33 = t241 * (-t240 * t69 + t244 * t71) - t245 * t148;
t27 = -t239 * t83 + t243 * t48 - t321;
t26 = qJ(6) * t301 - t239 * t34 - t322;
t25 = t253 + t337;
t24 = t239 * t48 + t243 * t83 + t274;
t23 = qJ(6) * t302 + t243 * t34 + t277;
t22 = -pkin(5) * t127 - qJ(6) * t103 + t30;
t19 = t236 * (t241 * t78 - t245 * t332) + (t242 * t50 - t246 * t77) * t235;
t18 = -pkin(5) * t53 + qJ(6) * t30;
t17 = -pkin(4) * t89 + pkin(10) * t21;
t16 = t236 * (-t102 * t245 + t241 * t62) + (t242 * t45 - t246 * t61) * t235;
t15 = -t20 - t323;
t14 = t21 + t276;
t13 = t243 * t30 + t315;
t11 = t236 * (-t127 * t245 + t241 * t39) + (t242 * t32 - t246 * t37) * t235;
t10 = t21 * t244 - t317;
t9 = t21 * t240 + t316;
t8 = -t22 * t239 + t243 * t25 - t323;
t7 = t10 * t245 + t241 * t89;
t6 = t22 * t243 + t239 * t25 + t276;
t5 = -t12 * t240 + t13 * t244;
t4 = t12 * t244 + t13 * t240;
t3 = -pkin(10) * t12 + qJ(6) * t313 - t18 * t239;
t2 = t241 * t53 + t245 * t5;
t1 = -pkin(4) * t53 + pkin(10) * t13 + qJ(6) * t315 + t18 * t243;
t29 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t291, 0, 0, 0, 0, 0, 0, (qJDD(2) * t246 - t242 * t247) * t235, (-qJDD(2) * t242 - t246 * t247) * t235, 0, t236 ^ 2 * t291 + (t242 * t173 - t246 * t266 - t252) * t235, 0, 0, 0, 0, 0, 0, t236 * (t215 * t245 + t221 * t241) + (t191 * t242 + t210 * t246) * t235, t236 * (-t214 * t241 + t220 * t245) + (t192 * t242 - t207 * t246) * t235, (t211 * t242 + t212 * t246) * t235, t236 * (-t146 * t245 + t147 * t241) + (t108 * t242 - t160 * t246) * t235, 0, 0, 0, 0, 0, 0, t236 * (t141 * t241 + t152 * t245) + (t109 * t242 - t140 * t246) * t235, t236 * (t145 * t241 + t156 * t245) + (t114 * t242 - t144 * t246) * t235, t236 * (t120 * t241 + t168 * t245) + (-t119 * t246 + t242 * t90) * t235, t236 * (-t129 * t245 + t241 * t57) + (t242 * t46 - t246 * t56) * t235, 0, 0, 0, 0, 0, 0, t16, t19, t11, t236 * (t10 * t241 - t245 * t89) + (t242 * t7 - t246 * t9) * t235, 0, 0, 0, 0, 0, 0, t16, t19, t11, t236 * (t241 * t5 - t245 * t53) + (t2 * t242 - t246 * t4) * t235; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t266, -t173, 0, 0, t264 * t241, t207 * t245 + t210 * t241, t294 + t245 * (-t229 + t330), -t263 * t245, t241 * (t230 - t330) + t292, 0, pkin(2) * t210 + pkin(8) * t191 - t160 * t245, -pkin(2) * t207 + pkin(8) * t192 + t160 * t241, pkin(2) * t212 + pkin(8) * t211 + t108, -pkin(2) * t160 + pkin(8) * t108, t241 * (t180 * t244 + t225 * t295) + t258, t241 * (t152 * t244 + t156 * t240) - t245 * t190, t241 * (-t196 * t240 + t339) - t245 * t155, t241 * (t197 * t244 - t240 * t254) - t258, t241 * (t195 * t244 + t300) + t245 * t151, t245 * t204 + t241 * (-t244 * t270 - t295) * t225, t241 * (-pkin(9) * t140 + t306) + t245 * (-pkin(3) * t140 + t86) - pkin(2) * t140 + pkin(8) * t109, t241 * (-pkin(9) * t144 + t305) + t245 * (-pkin(3) * t144 + t87) - pkin(2) * t144 + pkin(8) * t114, pkin(8) * t90 + t119 * t260 - t241 * t56, pkin(8) * t46 + t260 * t56, t55, t33, t51, t54, t52, t84, t241 * (-t240 * t43 + t244 * t58 - t325) + t245 * (t269 - t329) + t319, t241 * (-t240 * t47 + t244 * t65 - t324) + t245 * (-t282 - t328) + t318, t241 * (-t14 * t240 + t15 * t244 - t326) + t245 * t279 + t320, t241 * (-pkin(9) * t9 - pkin(10) * t316 - t17 * t240) + t245 * (-pkin(3) * t9 - t327) - pkin(2) * t9 + pkin(8) * t7, t55, t33, t51, t54, t52, t84, t241 * (-t23 * t240 + t244 * t26 - t325) + t245 * (-t329 - t341) + t319, t241 * (-t24 * t240 + t244 * t27 - t324) + t245 * (-t328 - t338) + t318, t241 * (-t240 * t6 + t244 * t8 - t326) + t245 * (t279 + t99) + t320, t241 * (-pkin(9) * t4 - t1 * t240 + t244 * t3) + t245 * (-pkin(3) * t4 - t278) - pkin(2) * t4 + pkin(8) * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t224, t229 - t230, t284, t224, t283, qJDD(3), -t146, -t147, 0, 0, t180 * t240 - t225 * t293, t152 * t240 - t156 * t244, t196 * t244 + t340, t197 * t240 + t244 * t254, t195 * t240 - t299, (-t240 * t270 + t293) * t225, pkin(3) * t152 + pkin(9) * t141 - t305, pkin(3) * t156 + pkin(9) * t145 + t306, pkin(3) * t168 + pkin(9) * t120 + t57, -pkin(3) * t129 + pkin(9) * t57, t64, t38, t81, t63, t82, t88, t240 * t58 + t244 * t43 + t310, t240 * t65 + t244 * t47 + t309, t14 * t244 + t15 * t240 + t311, -pkin(3) * t89 + pkin(9) * t10 - pkin(10) * t317 + t17 * t244, t64, t38, t81, t63, t82, t88, t23 * t244 + t240 * t26 + t310, t24 * t244 + t240 * t27 + t309, t240 * t8 + t244 * t6 + t311, -pkin(3) * t53 + pkin(9) * t5 + t1 * t244 + t240 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t261, t190, t155, t261, -t151, -t204, -t86, -t87, 0, 0, t150, t148, t106, -t150, -t103, -t201, -t269, t282, t67, t327, t150, t148, t106, -t150, -t103, -t201, t341, t338, -t99 + t67, t278; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150, t148, t106, -t150, -t103, -t201, -t41, -t42, 0, 0, t150, t148, t106, -t150, -t103, -t201, t251, t257, -t99, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, t332, t127, t53;];
tauJ_reg  = t29;
