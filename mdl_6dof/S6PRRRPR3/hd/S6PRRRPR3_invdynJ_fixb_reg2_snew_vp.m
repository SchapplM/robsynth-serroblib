% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 07:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRRRPR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 07:30:57
% EndTime: 2019-05-05 07:31:09
% DurationCPUTime: 7.21s
% Computational Cost: add. (15459->410), mult. (31377->563), div. (0->0), fcn. (22624->12), ass. (0->243)
t244 = cos(qJ(4));
t232 = qJD(3) + qJD(4);
t230 = t232 ^ 2;
t240 = sin(qJ(4));
t245 = cos(qJ(3));
t241 = sin(qJ(3));
t290 = t241 * t244;
t202 = (t245 * t240 + t290) * qJD(2);
t317 = t202 ^ 2;
t321 = -t317 - t230;
t231 = qJDD(3) + qJDD(4);
t285 = qJD(2) * t241;
t200 = -t244 * t245 * qJD(2) + t240 * t285;
t296 = t202 * t200;
t335 = t296 + t231;
t344 = t335 * t240;
t124 = -t244 * t321 + t344;
t343 = t335 * t244;
t126 = t240 * t321 + t343;
t77 = t241 * t124 - t245 * t126;
t362 = pkin(8) * t77;
t242 = sin(qJ(2));
t361 = t242 * t77;
t237 = cos(pkin(6));
t360 = t237 * (t124 * t245 + t126 * t241);
t359 = pkin(3) * t124;
t358 = pkin(9) * t124;
t357 = pkin(9) * t126;
t183 = -t317 + t230;
t161 = t296 - t231;
t350 = t161 * t244;
t351 = t161 * t240;
t355 = t241 * (t183 * t240 + t350) + t245 * (-t183 * t244 + t351);
t318 = t200 ^ 2;
t182 = t318 - t230;
t354 = t241 * (-t182 * t244 + t344) - t245 * (t182 * t240 + t343);
t155 = -t230 - t318;
t103 = t155 * t240 - t350;
t106 = -t155 * t244 - t351;
t66 = t103 * t241 + t245 * t106;
t353 = pkin(8) * t66;
t352 = t242 * t66;
t345 = t237 * (t103 * t245 - t106 * t241);
t342 = pkin(3) * t103;
t341 = pkin(9) * t103;
t340 = pkin(9) * t106;
t315 = -pkin(4) - pkin(10);
t319 = -t318 - t317;
t333 = pkin(2) * t319;
t332 = pkin(3) * t319;
t239 = sin(qJ(6));
t243 = cos(qJ(6));
t171 = -t243 * t200 + t232 * t239;
t173 = t200 * t239 + t232 * t243;
t130 = t173 * t171;
t282 = qJD(2) * qJD(3);
t276 = t245 * t282;
t281 = t241 * qJDD(2);
t207 = t276 + t281;
t227 = t245 * qJDD(2);
t277 = t241 * t282;
t267 = t227 - t277;
t259 = t244 * t207 + t240 * t267;
t142 = -t200 * qJD(4) + t259;
t138 = qJDD(6) + t142;
t324 = -t130 + t138;
t331 = t239 * t324;
t329 = t243 * t324;
t246 = cos(qJ(2));
t327 = t246 * t319;
t188 = t232 * t200;
t326 = t188 - t142;
t236 = sin(pkin(6));
t305 = sin(pkin(11));
t306 = cos(pkin(11));
t262 = t305 * g(1) - t306 * g(2);
t287 = -g(3) + qJDD(1);
t325 = t236 * t287 + t237 * t262;
t123 = t188 + t142;
t162 = pkin(4) * t200 - qJ(5) * t202;
t211 = -t306 * g(1) - t305 * g(2);
t151 = t246 * t211 + t325 * t242;
t248 = qJD(2) ^ 2;
t146 = -t248 * pkin(2) + qJDD(2) * pkin(8) + t151;
t181 = -t236 * t262 + t237 * t287;
t109 = t146 * t241 - t245 * t181;
t218 = t241 * t248 * t245;
t212 = qJDD(3) + t218;
t93 = (-t207 + t276) * pkin(9) + t212 * pkin(3) - t109;
t110 = t245 * t146 + t241 * t181;
t215 = qJD(3) * pkin(3) - pkin(9) * t285;
t234 = t245 ^ 2;
t229 = t234 * t248;
t94 = -pkin(3) * t229 + pkin(9) * t267 - qJD(3) * t215 + t110;
t57 = t240 * t94 - t244 * t93;
t49 = -t231 * pkin(4) - t230 * qJ(5) + t202 * t162 + qJDD(5) + t57;
t29 = pkin(5) * t123 + pkin(10) * t161 + t49;
t271 = t207 * t240 - t244 * t267;
t141 = qJD(4) * t202 + t271;
t180 = pkin(5) * t202 - pkin(10) * t232;
t269 = t211 * t242 - t325 * t246;
t145 = -qJDD(2) * pkin(2) - t248 * pkin(8) + t269;
t107 = -t267 * pkin(3) - pkin(9) * t229 + t215 * t285 + t145;
t253 = t141 * pkin(4) + t326 * qJ(5) + t107;
t274 = pkin(4) * t232 - (2 * qJD(5));
t35 = (-t180 + t274) * t202 - pkin(5) * t318 + pkin(10) * t141 + t253;
t18 = t239 * t35 - t243 * t29;
t19 = t239 * t29 + t243 * t35;
t8 = -t243 * t18 + t19 * t239;
t100 = -t171 * qJD(6) + t239 * t141 + t243 * t231;
t194 = qJD(6) + t202;
t149 = t194 * t171;
t323 = -t149 + t100;
t320 = t317 - t318;
t169 = t171 ^ 2;
t170 = t173 ^ 2;
t192 = t194 ^ 2;
t316 = 2 * qJD(5);
t314 = pkin(4) * t240;
t313 = pkin(4) * t244;
t58 = t240 * t93 + t244 * t94;
t257 = -t230 * pkin(4) + qJ(5) * t231 - t162 * t200 + t58;
t47 = t232 * t316 + t257;
t312 = -pkin(4) * t49 + qJ(5) * t47;
t32 = -t141 * pkin(5) - t318 * pkin(10) + (t316 + t180) * t232 + t257;
t310 = t239 * t32;
t96 = t130 + t138;
t309 = t239 * t96;
t26 = t240 * t58 - t244 * t57;
t308 = t241 * t26;
t30 = t243 * t32;
t307 = t243 * t96;
t298 = t194 * t239;
t297 = t194 * t243;
t295 = t202 * t232;
t294 = t232 * t240;
t293 = t232 * t244;
t292 = t240 * t107;
t291 = t241 * t212;
t289 = t244 * t107;
t213 = qJDD(3) - t218;
t288 = t245 * t213;
t284 = -qJD(4) + t232;
t117 = t284 * t202 - t271;
t120 = t284 * t200 + t259;
t286 = -pkin(4) * t120 + qJ(5) * t117;
t283 = qJD(4) + t232;
t280 = -t170 - t192;
t279 = t240 * t130;
t278 = t244 * t130;
t275 = qJ(5) * t240 + pkin(3);
t27 = t240 * t57 + t244 * t58;
t71 = t109 * t241 + t245 * t110;
t272 = -t243 * t141 + t231 * t239;
t270 = qJ(5) * t32 + t315 * t8;
t63 = t243 * t280 - t309;
t268 = qJ(5) * t323 + t315 * t63 + t30;
t208 = t227 - 0.2e1 * t277;
t9 = t239 * t18 + t19 * t243;
t23 = t240 * t47 - t244 * t49;
t24 = t240 * t49 + t244 * t47;
t10 = -t241 * t23 + t245 * t24;
t111 = -t192 - t169;
t60 = t239 * t111 + t329;
t85 = (qJD(6) + t194) * t173 + t272;
t264 = qJ(5) * t85 + t315 * t60 + t310;
t102 = -t169 - t170;
t261 = (-qJD(6) + t194) * t173 - t272;
t89 = t149 + t100;
t51 = t239 * t261 - t243 * t89;
t263 = qJ(5) * t102 + t315 * t51 - t8;
t258 = pkin(4) * t161 - qJ(5) * t155 + t49;
t256 = t241 * (t244 * t142 - t202 * t294) + t245 * (t240 * t142 + t202 * t293);
t255 = -t283 * t200 + t259;
t254 = t241 * (t141 * t240 + t200 * t293) + t245 * (-t244 * t141 + t200 * t294);
t252 = -pkin(4) * t321 + qJ(5) * t335 + t47;
t251 = (t241 * (-t200 * t244 + t202 * t240) + t245 * (-t200 * t240 - t202 * t244)) * t232;
t249 = -t202 * t316 + t253;
t247 = qJD(3) ^ 2;
t233 = t241 ^ 2;
t228 = t233 * t248;
t217 = -t229 - t247;
t216 = -t228 - t247;
t210 = t228 + t229;
t209 = (t233 + t234) * qJDD(2);
t206 = 0.2e1 * t276 + t281;
t176 = -t216 * t241 - t288;
t175 = t217 * t245 - t291;
t148 = -t170 + t192;
t147 = t169 - t192;
t129 = t170 - t169;
t116 = t141 + t295;
t115 = t141 - t295;
t114 = t283 * t202 + t271;
t99 = -qJD(6) * t173 - t272;
t98 = (-t171 * t243 + t173 * t239) * t194;
t97 = (t171 * t239 + t173 * t243) * t194;
t81 = t100 * t243 - t173 * t298;
t80 = -t100 * t239 - t173 * t297;
t79 = t171 * t297 - t239 * t99;
t78 = -t171 * t298 - t243 * t99;
t75 = -t115 * t244 + t123 * t240;
t74 = t117 * t244 + t120 * t240;
t73 = -t115 * t240 - t123 * t244;
t72 = t117 * t240 - t120 * t244;
t70 = t147 * t243 - t309;
t69 = -t148 * t239 + t329;
t68 = -t147 * t239 - t307;
t67 = -t148 * t243 - t331;
t64 = -t239 * t280 - t307;
t61 = t111 * t243 - t331;
t55 = t202 * t274 + t253;
t54 = -t239 * t323 - t243 * t85;
t53 = t239 * t89 + t243 * t261;
t52 = t239 * t85 - t243 * t323;
t45 = -t241 * t73 + t245 * t75;
t44 = -t241 * t72 + t245 * t74;
t43 = (t116 + t295) * pkin(4) + t249;
t42 = -pkin(4) * t295 + qJ(5) * t255 - t249;
t41 = t240 * t63 + t244 * t323;
t40 = t240 * t323 - t244 * t63;
t39 = t240 * t60 + t244 * t85;
t38 = t240 * t85 - t244 * t60;
t37 = -qJ(5) * t319 + t49;
t36 = -pkin(4) * t319 + t47;
t34 = t102 * t244 + t240 * t51;
t33 = t102 * t240 - t244 * t51;
t25 = pkin(5) * t51 - qJ(5) * t53;
t22 = -t241 * t40 + t245 * t41;
t21 = -t241 * t38 + t245 * t39;
t20 = -t241 * t33 + t245 * t34;
t15 = pkin(5) * t323 + t315 * t64 - t310;
t14 = pkin(5) * t85 + t315 * t61 + t30;
t13 = t245 * t27 - t308;
t12 = pkin(5) * t63 - qJ(5) * t64 - t19;
t11 = pkin(5) * t60 - qJ(5) * t61 - t18;
t6 = t240 * t8 + t244 * t32;
t5 = t240 * t32 - t244 * t8;
t4 = pkin(5) * t102 + t315 * t53 - t9;
t3 = pkin(5) * t8 - qJ(5) * t9;
t2 = pkin(5) * t32 + t315 * t9;
t1 = -t241 * t5 + t245 * t6;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t287, 0, 0, 0, 0, 0, 0, (qJDD(2) * t246 - t242 * t248) * t236, (-qJDD(2) * t242 - t246 * t248) * t236, 0, t181 * t237 + (t151 * t242 - t246 * t269) * t236, 0, 0, 0, 0, 0, 0, t237 * (t212 * t245 + t217 * t241) + (t175 * t242 + t208 * t246) * t236, t237 * (-t213 * t241 + t216 * t245) + (t176 * t242 - t206 * t246) * t236, (t209 * t242 + t210 * t246) * t236, t237 * (-t109 * t245 + t110 * t241) + (-t145 * t246 + t242 * t71) * t236, 0, 0, 0, 0, 0, 0, t345 + (-t114 * t246 - t352) * t236, -t360 + (t246 * t326 + t361) * t236, t237 * (t241 * t75 + t245 * t73) + (t242 * t45 - t327) * t236, t237 * (t241 * t27 + t245 * t26) + (-t107 * t246 + t13 * t242) * t236, 0, 0, 0, 0, 0, 0, t237 * (t241 * t74 + t245 * t72) + (t242 * t44 - t327) * t236, -t345 + (t116 * t246 + t352) * t236, t360 + (t246 * t255 - t361) * t236, t237 * (t23 * t245 + t24 * t241) + (t10 * t242 - t246 * t55) * t236, 0, 0, 0, 0, 0, 0, t237 * (t241 * t39 + t245 * t38) + (t21 * t242 - t246 * t61) * t236, t237 * (t241 * t41 + t245 * t40) + (t22 * t242 - t246 * t64) * t236, t237 * (t241 * t34 + t245 * t33) + (t20 * t242 - t246 * t53) * t236, t237 * (t241 * t6 + t245 * t5) + (t1 * t242 - t246 * t9) * t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t269, -t151, 0, 0, (t207 + t276) * t241, t206 * t245 + t208 * t241, t291 + t245 * (-t228 + t247), t208 * t245, t241 * (t229 - t247) + t288, 0, pkin(2) * t208 + pkin(8) * t175 - t145 * t245, -pkin(2) * t206 + pkin(8) * t176 + t145 * t241, pkin(2) * t210 + pkin(8) * t209 + t71, -pkin(2) * t145 + pkin(8) * t71, t256, t241 * (-t114 * t244 + t240 * t326) + t245 * (-t114 * t240 - t244 * t326), -t355, t254, -t354, t251, t241 * (t292 - t341) + t245 * (-pkin(3) * t114 - t289 - t340) - pkin(2) * t114 - t353, t241 * (t289 + t358) + t245 * (pkin(3) * t326 + t292 - t357) + pkin(2) * t326 + t362, t241 * (-pkin(9) * t73 - t26) + t245 * (pkin(9) * t75 + t27 - t332) - t333 + pkin(8) * t45, -pkin(9) * t308 + t245 * (-pkin(3) * t107 + pkin(9) * t27) - pkin(2) * t107 + pkin(8) * t13, t251, t355, t354, t256, t241 * (-t116 * t244 - t240 * t255) + t245 * (-t116 * t240 + t244 * t255), t254, t241 * (-pkin(9) * t72 - t240 * t36 + t244 * t37) + t245 * (pkin(9) * t74 + t240 * t37 + t244 * t36 - t332) - t333 + pkin(8) * t44, t241 * (-t240 * t43 + t341) + t245 * (t244 * t43 + t340) + t353 + (qJ(5) * t290 + t245 * t275 + pkin(2)) * t116, t241 * (t244 * t42 - t358) + t245 * (t240 * t42 + t357) - t362 + (-t241 * t314 + t245 * (pkin(3) + t313) + pkin(2)) * t255, (t241 * (-qJ(5) * t244 + t314) + t245 * (-t275 - t313) - pkin(2)) * t55 + (pkin(8) + pkin(9)) * t10, t241 * (-t240 * t80 + t278) + t245 * (t244 * t80 + t279), t241 * (t129 * t244 - t240 * t52) + t245 * (t129 * t240 + t244 * t52), t241 * (-t240 * t67 + t244 * t89) + t245 * (t240 * t89 + t244 * t67), t241 * (-t240 * t78 - t278) + t245 * (t244 * t78 - t279), t241 * (-t240 * t68 + t244 * t261) + t245 * (t240 * t261 + t244 * t68), t241 * (t138 * t244 - t240 * t97) + t245 * (t138 * t240 + t244 * t97), t241 * (-pkin(9) * t38 + t244 * t11 - t240 * t14) + t245 * (-pkin(3) * t61 + pkin(9) * t39 + t240 * t11 + t244 * t14) - pkin(2) * t61 + pkin(8) * t21, t241 * (-pkin(9) * t40 + t12 * t244 - t15 * t240) + t245 * (-pkin(3) * t64 + pkin(9) * t41 + t12 * t240 + t15 * t244) - pkin(2) * t64 + pkin(8) * t22, t241 * (-pkin(9) * t33 - t240 * t4 + t244 * t25) + t245 * (-pkin(3) * t53 + pkin(9) * t34 + t240 * t25 + t244 * t4) - pkin(2) * t53 + pkin(8) * t20, t241 * (-pkin(9) * t5 - t2 * t240 + t244 * t3) + t245 * (-pkin(3) * t9 + pkin(9) * t6 + t2 * t244 + t240 * t3) - pkin(2) * t9 + pkin(8) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t218, -t229 + t228, t281, t218, t227, qJDD(3), -t109, -t110, 0, 0, t296, t320, t123, -t296, t117, t231, -t57 + t342, -t58 - t359, pkin(3) * t73, pkin(3) * t26, t231, -t120, t115, t296, t320, -t296, pkin(3) * t72 + t286, t258 - t342, t252 + t359, pkin(3) * t23 + t312, t81, t54, t69, t79, t70, t98, pkin(3) * t38 + t264, pkin(3) * t40 + t268, pkin(3) * t33 + t263, pkin(3) * t5 + t270; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t296, t320, t123, -t296, t117, t231, -t57, -t58, 0, 0, t231, -t120, t115, t296, t320, -t296, t286, t258, t252, t312, t81, t54, t69, t79, t70, t98, t264, t268, t263, t270; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, -t161, t321, t49, 0, 0, 0, 0, 0, 0, t60, t63, t51, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, t129, t89, -t130, t261, t138, -t18, -t19, 0, 0;];
tauJ_reg  = t7;
