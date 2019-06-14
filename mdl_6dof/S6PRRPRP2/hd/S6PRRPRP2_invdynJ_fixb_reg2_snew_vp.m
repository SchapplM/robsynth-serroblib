% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 03:50
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRRPRP2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP2_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:47:37
% EndTime: 2019-05-05 03:47:53
% DurationCPUTime: 7.08s
% Computational Cost: add. (16063->430), mult. (34970->586), div. (0->0), fcn. (25577->12), ass. (0->260)
t223 = sin(pkin(11));
t225 = cos(pkin(11));
t230 = sin(qJ(3));
t233 = cos(qJ(3));
t279 = qJD(2) * t233;
t280 = qJD(2) * t230;
t190 = t223 * t280 - t225 * t279;
t187 = qJD(5) + t190;
t318 = t187 ^ 2;
t192 = t223 * t279 + t225 * t280;
t229 = sin(qJ(5));
t232 = cos(qJ(5));
t174 = -qJD(3) * t232 + t192 * t229;
t319 = t174 ^ 2;
t148 = t319 - t318;
t273 = qJD(2) * qJD(3);
t266 = t233 * t273;
t272 = t230 * qJDD(2);
t197 = t266 + t272;
t215 = t233 * qJDD(2);
t267 = t230 * t273;
t252 = t215 - t267;
t261 = t197 * t223 - t225 * t252;
t166 = qJDD(5) + t261;
t176 = qJD(3) * t229 + t192 * t232;
t299 = t176 * t174;
t113 = -t299 - t166;
t291 = t229 * t113;
t83 = -t148 * t232 - t291;
t153 = t176 * t187;
t169 = t225 * t197 + t223 * t252;
t262 = qJDD(3) * t232 - t229 * t169;
t248 = qJD(5) * t176 - t262;
t93 = -t153 + t248;
t386 = t230 * (t223 * t93 + t225 * t83) + t233 * (t223 * t83 - t225 * t93);
t173 = t176 ^ 2;
t331 = -t173 - t318;
t77 = t232 * t331 + t291;
t385 = pkin(2) * t77;
t384 = pkin(3) * t77;
t383 = pkin(4) * t77;
t382 = pkin(9) * t77;
t287 = t232 * t113;
t79 = -t229 * t331 + t287;
t381 = pkin(9) * t79;
t380 = t223 * t79;
t379 = t225 * t79;
t234 = cos(qJ(2));
t378 = t234 * t77;
t330 = t173 - t319;
t251 = -t229 * qJDD(3) - t232 * t169;
t244 = -qJD(5) * t174 - t251;
t300 = t174 * t187;
t327 = -t300 + t244;
t309 = t229 * t327;
t332 = t153 + t248;
t63 = t232 * t332 + t309;
t375 = t230 * (-t223 * t330 + t225 * t63) + t233 * (t223 * t63 + t225 * t330);
t328 = -t299 + t166;
t286 = t232 * t328;
t325 = -t318 - t319;
t335 = t229 * t325 + t286;
t290 = t229 * t328;
t334 = t232 * t325 - t290;
t351 = t223 * t332 + t225 * t334;
t352 = t223 * t334 - t225 * t332;
t366 = -t230 * t352 + t233 * t351;
t374 = -pkin(2) * t335 + pkin(8) * t366;
t224 = sin(pkin(6));
t226 = cos(pkin(6));
t231 = sin(qJ(2));
t371 = t226 * (t230 * t351 + t233 * t352) + (t231 * t366 - t234 * t335) * t224;
t370 = qJ(4) * t352;
t369 = -t148 * t229 + t287;
t368 = pkin(3) * t352 + pkin(9) * t334;
t367 = -pkin(3) * t335 + qJ(4) * t351;
t326 = t300 + t244;
t149 = -t173 + t318;
t353 = -t149 * t229 + t286;
t364 = t230 * (t223 * t326 + t225 * t353) + t233 * (t223 * t353 - t225 * t326);
t361 = pkin(4) * t335;
t359 = pkin(9) * t335;
t358 = qJ(6) * t327;
t354 = t149 * t232 + t290;
t329 = t173 + t319;
t350 = pkin(4) * t329;
t167 = t192 * t190;
t324 = qJDD(3) - t167;
t349 = t223 * t324;
t347 = t223 * t329;
t344 = t225 * t324;
t342 = t225 * t329;
t304 = sin(pkin(10));
t305 = cos(pkin(10));
t245 = g(1) * t304 - g(2) * t305;
t243 = t226 * t245;
t283 = -g(3) + qJDD(1);
t337 = t224 * t283 + t243;
t274 = t192 * qJD(3);
t140 = t261 + t274;
t201 = -g(1) * t305 - g(2) * t304;
t156 = t234 * t201 + t231 * t337;
t336 = qJDD(2) * pkin(8) + t156;
t333 = -t229 * t332 + t232 * t327;
t202 = qJD(3) * pkin(3) - qJ(4) * t280;
t220 = t233 ^ 2;
t240 = -t224 * t245 + t226 * t283;
t238 = t230 * t240;
t303 = qJ(4) * t230;
t106 = t233 * t336 + t238 - qJD(3) * t202 + t215 * qJ(4) + (-qJD(3) * t303 + (-pkin(2) * t233 - pkin(3) * t220) * qJD(2)) * qJD(2);
t317 = qJD(2) ^ 2;
t237 = -pkin(2) * t317 + t336;
t121 = t230 * t237 - t233 * t240;
t211 = t230 * t317 * t233;
t203 = qJDD(3) + t211;
t236 = -t121 + (-t197 + t266) * qJ(4) + t203 * pkin(3);
t60 = -0.2e1 * qJD(4) * t190 + t106 * t225 + t223 * t236;
t323 = pkin(5) * t248 - t358;
t133 = pkin(5) * t174 - qJ(6) * t176;
t158 = pkin(4) * t190 - pkin(9) * t192;
t235 = qJD(3) ^ 2;
t49 = -pkin(4) * t235 + qJDD(3) * pkin(9) - t158 * t190 + t60;
t254 = t231 * t201 - t234 * t337;
t146 = -qJDD(2) * pkin(2) - pkin(8) * t317 + t254;
t218 = t220 * t317;
t114 = -pkin(3) * t252 - qJ(4) * t218 + t202 * t280 + qJDD(4) + t146;
t275 = t190 * qJD(3);
t259 = -t169 + t275;
t70 = pkin(4) * t140 + pkin(9) * t259 + t114;
t37 = t229 * t70 + t232 * t49;
t260 = qJ(6) * t166 - t133 * t174 + t37;
t322 = -(t331 + t318) * pkin(5) - qJ(6) * t113 + t260;
t297 = t187 * t232;
t269 = t174 * t297;
t249 = t229 * t248 + t269;
t270 = t225 * t299;
t271 = t223 * t299;
t321 = t230 * (t225 * t249 - t271) + t233 * (t223 * t249 + t270);
t298 = t187 * t229;
t145 = t176 * t298;
t255 = t145 - t269;
t320 = t230 * (t166 * t223 + t225 * t255) + t233 * (-t166 * t225 + t223 * t255);
t188 = t190 ^ 2;
t189 = t192 ^ 2;
t316 = pkin(4) * t223;
t315 = pkin(5) * t232;
t263 = t223 * t106 - t225 * t236;
t250 = -qJDD(3) * pkin(4) - pkin(9) * t235 + t263;
t48 = (0.2e1 * qJD(4) + t158) * t192 + t250;
t312 = t229 * t48;
t310 = t229 * t326;
t308 = t232 * t48;
t306 = t232 * t326;
t302 = qJ(6) * t232;
t296 = t223 * t114;
t162 = qJDD(3) + t167;
t295 = t223 * t162;
t293 = t225 * t114;
t292 = t225 * t162;
t288 = t230 * t203;
t204 = qJDD(3) - t211;
t285 = t233 * t204;
t277 = qJD(4) * t192;
t276 = qJD(6) * t187;
t268 = -pkin(4) * t225 - pkin(3);
t265 = -qJ(6) * t229 - pkin(4);
t59 = t263 + 0.2e1 * t277;
t39 = t223 * t59 + t225 * t60;
t36 = t229 * t49 - t232 * t70;
t14 = t229 * t36 + t232 * t37;
t122 = t233 * t237 + t238;
t76 = t121 * t230 + t122 * t233;
t179 = 0.2e1 * t276;
t253 = t179 + t260;
t26 = -pkin(5) * t318 + t253;
t29 = -pkin(5) * t166 - qJ(6) * t318 + t133 * t176 + qJDD(6) + t36;
t258 = -pkin(5) * t29 + qJ(6) * t26;
t257 = -pkin(5) * t326 - qJ(6) * t93;
t256 = t174 * t298 - t232 * t248;
t8 = t14 * t223 - t225 * t48;
t9 = t14 * t225 + t223 * t48;
t2 = -t230 * t8 + t233 * t9;
t198 = t215 - 0.2e1 * t267;
t38 = t223 * t60 - t225 * t59;
t13 = t229 * t37 - t232 * t36;
t141 = -t261 + t274;
t246 = (-t174 * t229 - t176 * t232) * t187;
t91 = t232 * t244 - t145;
t242 = t230 * (t225 * t91 + t271) + t233 * (t223 * t91 - t270);
t185 = -0.2e1 * t277;
t241 = 0.2e1 * qJD(6) * t176 - t192 * t158 + t185 - t250 - t323;
t239 = pkin(5) * t328 + qJ(6) * t325 - t29;
t219 = t230 ^ 2;
t216 = t219 * t317;
t210 = -t218 - t235;
t209 = -t216 - t235;
t200 = t216 + t218;
t199 = (t219 + t220) * qJDD(2);
t196 = 0.2e1 * t266 + t272;
t182 = -t189 - t235;
t181 = -t189 + t235;
t180 = t188 - t235;
t178 = -t209 * t230 - t285;
t177 = t210 * t233 - t288;
t159 = -t235 - t188;
t143 = t169 + t275;
t136 = -t188 - t189;
t130 = -t182 * t223 - t292;
t129 = t182 * t225 - t295;
t117 = t159 * t225 - t349;
t116 = t159 * t223 + t344;
t105 = t141 * t225 + t143 * t223;
t104 = t141 * t223 - t143 * t225;
t100 = (qJD(5) + t187) * t174 + t251;
t95 = (-qJD(5) + t187) * t176 + t262;
t90 = t176 * t297 + t229 * t244;
t85 = -t129 * t230 + t130 * t233;
t71 = -t116 * t230 + t117 * t233;
t67 = -t104 * t230 + t105 * t233;
t65 = -t232 * t93 + t310;
t64 = t232 * t95 + t310;
t62 = -t229 * t93 - t306;
t61 = t229 * t95 - t306;
t56 = -t100 * t223 + t379;
t54 = t100 * t225 + t380;
t52 = -t223 * t327 - t379;
t50 = t225 * t327 - t380;
t46 = t225 * t65 - t347;
t45 = t225 * t64 - t347;
t44 = t223 * t65 + t342;
t43 = t223 * t64 + t342;
t42 = t308 - t382;
t41 = -pkin(4) * t62 - t257;
t40 = t312 - t359;
t34 = (pkin(5) * t187 - 0.2e1 * qJD(6)) * t176 + t48 + t323;
t32 = -t230 * t54 + t233 * t56;
t30 = -t230 * t50 + t233 * t52;
t28 = -t230 * t44 + t233 * t46;
t27 = -t230 * t43 + t233 * t45;
t25 = t37 - t383;
t24 = t36 - t361;
t23 = (-t332 - t153) * pkin(5) + t241;
t22 = -pkin(5) * t153 + t241 + t358;
t21 = qJ(6) * t329 + t29;
t20 = (t329 - t318) * pkin(5) + t253;
t19 = -t239 - t361;
t18 = -t229 * t23 - t302 * t332 - t359;
t17 = -pkin(5) * t309 + t22 * t232 + t382;
t16 = -0.2e1 * t276 - t322 + t383;
t15 = -t230 * t38 + t233 * t39;
t12 = -pkin(9) * t61 - t13;
t11 = t229 * t29 + t232 * t26;
t10 = t229 * t26 - t232 * t29;
t7 = -pkin(9) * t62 - t20 * t229 + t21 * t232;
t6 = t11 * t225 + t223 * t34;
t5 = t11 * t223 - t225 * t34;
t4 = -pkin(9) * t10 + (pkin(5) * t229 - t302) * t34;
t3 = -pkin(4) * t10 - t258;
t1 = -t230 * t5 + t233 * t6;
t31 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t283, 0, 0, 0, 0, 0, 0, (qJDD(2) * t234 - t231 * t317) * t224, (-qJDD(2) * t231 - t234 * t317) * t224, 0, t226 ^ 2 * t283 + (t156 * t231 - t234 * t254 - t243) * t224, 0, 0, 0, 0, 0, 0, t226 * (t203 * t233 + t210 * t230) + (t177 * t231 + t198 * t234) * t224, t226 * (-t204 * t230 + t209 * t233) + (t178 * t231 - t196 * t234) * t224, (t199 * t231 + t200 * t234) * t224, t226 * (-t121 * t233 + t122 * t230) + (-t146 * t234 + t231 * t76) * t224, 0, 0, 0, 0, 0, 0, t226 * (t116 * t233 + t117 * t230) + (-t140 * t234 + t231 * t71) * t224, t226 * (t129 * t233 + t130 * t230) + (t231 * t85 + t234 * t259) * t224, t226 * (t104 * t233 + t105 * t230) + (-t136 * t234 + t231 * t67) * t224, t226 * (t230 * t39 + t233 * t38) + (-t114 * t234 + t15 * t231) * t224, 0, 0, 0, 0, 0, 0, t371, t226 * (t230 * t56 + t233 * t54) + (t231 * t32 - t378) * t224, t226 * (t230 * t45 + t233 * t43) + (t231 * t27 - t234 * t61) * t224, t226 * (t230 * t9 + t233 * t8) + (-t13 * t234 + t2 * t231) * t224, 0, 0, 0, 0, 0, 0, t371, t226 * (t230 * t46 + t233 * t44) + (t231 * t28 - t234 * t62) * t224, t226 * (t230 * t52 + t233 * t50) + (t231 * t30 + t378) * t224, t226 * (t230 * t6 + t233 * t5) + (t1 * t231 - t10 * t234) * t224; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t254, -t156, 0, 0, (t197 + t266) * t230, t196 * t233 + t198 * t230, t288 + t233 * (-t216 + t235), t198 * t233, t230 * (t218 - t235) + t285, 0, pkin(2) * t198 + pkin(8) * t177 - t146 * t233, -pkin(2) * t196 + pkin(8) * t178 + t146 * t230, pkin(2) * t200 + pkin(8) * t199 + t76, -pkin(2) * t146 + pkin(8) * t76, t230 * (t169 * t225 - t223 * t274) + t233 * (t169 * t223 + t225 * t274), t230 * (-t140 * t225 + t223 * t259) + t233 * (-t140 * t223 - t225 * t259), t230 * (-t181 * t223 + t344) + t233 * (t181 * t225 + t349), t230 * (t223 * t261 + t225 * t275) + t233 * (t223 * t275 - t225 * t261), t230 * (t180 * t225 - t295) + t233 * (t180 * t223 + t292), (t230 * (-t190 * t225 + t192 * t223) + t233 * (-t190 * t223 - t192 * t225)) * qJD(3), t230 * (-qJ(4) * t116 + t296) + t233 * (-pkin(3) * t140 + qJ(4) * t117 - t293) - pkin(2) * t140 + pkin(8) * t71, t230 * (-qJ(4) * t129 + t293) + t233 * (pkin(3) * t259 + qJ(4) * t130 + t296) + pkin(2) * t259 + pkin(8) * t85, t230 * (-qJ(4) * t104 - t38) + t233 * (-pkin(3) * t136 + qJ(4) * t105 + t39) - pkin(2) * t136 + pkin(8) * t67, -t38 * t303 + t233 * (-pkin(3) * t114 + qJ(4) * t39) - pkin(2) * t114 + pkin(8) * t15, t242, -t375, t364, t321, -t386, t320, t230 * (-t223 * t24 + t225 * t40 - t370) + t233 * (t223 * t40 + t225 * t24 + t367) + t374, t230 * (-qJ(4) * t54 - t223 * t25 + t225 * t42) + t233 * (qJ(4) * t56 + t223 * t42 + t225 * t25 - t384) - t385 + pkin(8) * t32, t230 * (-qJ(4) * t43 + t12 * t225) + t233 * (qJ(4) * t45 + t223 * t12) + pkin(8) * t27 + (t230 * t316 + t233 * t268 - pkin(2)) * t61, (t230 * (-pkin(9) * t225 + t316) + t233 * (-pkin(9) * t223 + t268) - pkin(2)) * t13 + (pkin(8) + qJ(4)) * t2, t242, t364, t375, t320, t386, t321, t230 * (t18 * t225 - t19 * t223 - t370) + t233 * (t18 * t223 + t19 * t225 + t367) + t374, t230 * (-qJ(4) * t44 - t223 * t41 + t225 * t7) + t233 * (-pkin(3) * t62 + qJ(4) * t46 + t223 * t7 + t225 * t41) - pkin(2) * t62 + pkin(8) * t28, t230 * (-qJ(4) * t50 - t16 * t223 + t17 * t225) + t233 * (qJ(4) * t52 + t16 * t225 + t17 * t223 + t384) + t385 + pkin(8) * t30, t230 * (-qJ(4) * t5 - t223 * t3 + t225 * t4) + t233 * (-pkin(3) * t10 + qJ(4) * t6 + t223 * t4 + t225 * t3) - pkin(2) * t10 + pkin(8) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t211, t216 - t218, t272, t211, t215, qJDD(3), -t121, -t122, 0, 0, t167, t189 - t188, t143, -t167, t141, qJDD(3), pkin(3) * t116 + t185 - t263, pkin(3) * t129 - t60, pkin(3) * t104, pkin(3) * t38, t90, t333, t354, t256, -t369, t246, -pkin(4) * t332 - t308 + t368, pkin(3) * t54 + pkin(4) * t100 + t312 + t381, pkin(3) * t43 + pkin(9) * t64 + t14 + t350, pkin(3) * t8 - pkin(4) * t48 + pkin(9) * t14, t90, t354, -t333, t246, t369, t256, t232 * t23 + t265 * t332 + t368, pkin(3) * t44 + pkin(9) * t65 + t20 * t232 + t21 * t229 + t350, pkin(3) * t50 - t381 + t229 * t22 + (pkin(4) + t315) * t327, pkin(3) * t5 + pkin(9) * t11 + (t265 - t315) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140, -t259, t136, t114, 0, 0, 0, 0, 0, 0, t335, t77, t61, t13, 0, 0, 0, 0, 0, 0, t335, t62, -t77, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t299, t330, t326, -t299, -t93, t166, -t36, -t37, 0, 0, t299, t326, -t330, t166, t93, -t299, t239, t257, t179 + t322, t258; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t328, t326, t331, t29;];
tauJ_reg  = t31;
