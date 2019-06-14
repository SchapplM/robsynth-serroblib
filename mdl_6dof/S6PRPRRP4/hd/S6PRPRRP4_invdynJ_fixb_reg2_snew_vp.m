% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 23:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRPRRP4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP4_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:50:33
% EndTime: 2019-05-04 23:50:48
% DurationCPUTime: 6.54s
% Computational Cost: add. (14092->394), mult. (31565->550), div. (0->0), fcn. (24250->12), ass. (0->245)
t213 = sin(pkin(11));
t215 = cos(pkin(11));
t219 = sin(qJ(4));
t222 = cos(qJ(4));
t267 = qJD(2) * t215;
t272 = t213 * t219;
t188 = qJD(2) * t272 - t222 * t267;
t184 = qJD(5) + t188;
t304 = t184 ^ 2;
t244 = t213 * t222 + t215 * t219;
t190 = t244 * qJD(2);
t218 = sin(qJ(5));
t221 = cos(qJ(5));
t172 = -t221 * qJD(4) + t190 * t218;
t305 = t172 ^ 2;
t144 = t305 - t304;
t261 = t215 * qJDD(2);
t262 = t213 * qJDD(2);
t245 = t219 * t262 - t222 * t261;
t263 = t190 * qJD(4);
t164 = -t245 - t263;
t156 = qJDD(5) - t164;
t174 = qJD(4) * t218 + t190 * t221;
t275 = t174 * t172;
t109 = -t275 - t156;
t286 = t109 * t218;
t83 = -t144 * t221 - t286;
t149 = t174 * t184;
t187 = t244 * qJDD(2);
t264 = t188 * qJD(4);
t166 = t187 - t264;
t253 = t221 * qJDD(4) - t166 * t218;
t241 = qJD(5) * t174 - t253;
t93 = -t149 + t241;
t372 = t213 * (t219 * t93 + t222 * t83) + t215 * (t219 * t83 - t222 * t93);
t171 = t174 ^ 2;
t317 = -t171 - t304;
t77 = t221 * t317 + t286;
t371 = pkin(2) * t77;
t370 = pkin(3) * t77;
t369 = pkin(4) * t77;
t368 = pkin(9) * t77;
t285 = t109 * t221;
t79 = -t218 * t317 + t285;
t367 = pkin(9) * t79;
t366 = t219 * t79;
t365 = t222 * t79;
t223 = cos(qJ(2));
t364 = t223 * t77;
t316 = t171 - t305;
t243 = -qJDD(4) * t218 - t166 * t221;
t237 = -qJD(5) * t172 - t243;
t276 = t172 * t184;
t312 = -t276 + t237;
t295 = t218 * t312;
t318 = t149 + t241;
t60 = t221 * t318 + t295;
t361 = t213 * (-t219 * t316 + t222 * t60) + t215 * (t219 * t60 + t222 * t316);
t313 = -t275 + t156;
t283 = t313 * t221;
t310 = -t304 - t305;
t323 = t218 * t310 + t283;
t284 = t313 * t218;
t322 = t221 * t310 - t284;
t339 = t219 * t318 + t222 * t322;
t340 = t219 * t322 - t222 * t318;
t353 = -t213 * t340 + t215 * t339;
t358 = -pkin(2) * t323 + qJ(3) * t353;
t214 = sin(pkin(6));
t216 = cos(pkin(6));
t220 = sin(qJ(2));
t357 = t216 * (t213 * t339 + t215 * t340) + (t220 * t353 - t223 * t323) * t214;
t356 = pkin(8) * t340;
t355 = -t144 * t218 + t285;
t354 = -pkin(3) * t323 + pkin(8) * t339;
t311 = t276 + t237;
t145 = -t171 + t304;
t341 = -t145 * t218 + t283;
t351 = t213 * (t219 * t311 + t222 * t341) + t215 * (t219 * t341 - t222 * t311);
t348 = pkin(4) * t323;
t347 = pkin(9) * t322;
t346 = pkin(9) * t323;
t342 = t221 * t145 + t284;
t315 = t171 + t305;
t338 = pkin(4) * t315;
t337 = qJ(6) * t312;
t167 = t190 * t188;
t309 = qJDD(4) - t167;
t335 = t219 * t309;
t333 = t219 * t315;
t329 = t222 * t309;
t327 = t222 * t315;
t208 = t213 ^ 2;
t209 = t215 ^ 2;
t303 = qJD(2) ^ 2;
t195 = (t208 + t209) * t303;
t288 = sin(pkin(10));
t289 = cos(pkin(10));
t238 = t288 * g(1) - t289 * g(2);
t236 = t216 * t238;
t269 = -g(3) + qJDD(1);
t324 = t214 * t269 + t236;
t196 = -t289 * g(1) - t288 * g(2);
t153 = t223 * t196 + t324 * t220;
t228 = qJDD(2) * qJ(3) + t153;
t232 = -t214 * t238 + t216 * t269;
t301 = 2 * qJD(3);
t112 = t215 * (-t303 * pkin(2) + t228) + t213 * t232 + t267 * t301;
t266 = t303 * t209;
t102 = -pkin(3) * t266 + pkin(8) * t261 + t112;
t229 = t215 * t232;
t320 = qJ(3) + pkin(8);
t226 = t229 + (-t320 * qJDD(2) + (-(2 * qJD(3)) + (t215 * pkin(3) + pkin(2)) * qJD(2)) * qJD(2) - t153) * t213;
t66 = t222 * t102 + t219 * t226;
t321 = -t218 * t318 + t221 * t312;
t131 = pkin(5) * t172 - qJ(6) * t174;
t157 = pkin(4) * t188 - pkin(9) * t190;
t302 = qJD(4) ^ 2;
t57 = -t302 * pkin(4) + qJDD(4) * pkin(9) - t188 * t157 + t66;
t212 = qJDD(2) * pkin(2);
t247 = t196 * t220 - t324 * t223;
t141 = -t303 * qJ(3) + qJDD(3) - t212 + t247;
t126 = -pkin(3) * t261 + t141 + (-t208 * t303 - t266) * pkin(8);
t69 = (-t166 + t264) * pkin(9) + (-t164 + t263) * pkin(4) + t126;
t38 = t218 * t69 + t221 * t57;
t252 = t156 * qJ(6) - t172 * t131 + t38;
t308 = -(t317 + t304) * pkin(5) - qJ(6) * t109 + t252;
t273 = t184 * t221;
t258 = t172 * t273;
t242 = t218 * t241 + t258;
t259 = t222 * t275;
t260 = t219 * t275;
t307 = t213 * (t222 * t242 - t260) + t215 * (t219 * t242 + t259);
t274 = t184 * t218;
t142 = t174 * t274;
t248 = t142 - t258;
t306 = t213 * (t156 * t219 + t222 * t248) + t215 * (-t222 * t156 + t219 * t248);
t185 = t188 ^ 2;
t186 = t190 ^ 2;
t300 = pkin(5) * t221;
t65 = t102 * t219 - t222 * t226;
t41 = t219 * t66 - t222 * t65;
t299 = t213 * t41;
t56 = -qJDD(4) * pkin(4) - t302 * pkin(9) + t157 * t190 + t65;
t298 = t218 * t56;
t296 = t218 * t311;
t293 = t221 * t56;
t291 = t221 * t311;
t287 = qJ(6) * t221;
t282 = t126 * t219;
t281 = t126 * t222;
t161 = qJDD(4) + t167;
t278 = t161 * t219;
t277 = t161 * t222;
t265 = qJD(6) * t184;
t257 = -pkin(4) * t222 - pkin(3);
t256 = -qJ(6) * t218 - pkin(4);
t37 = t218 * t57 - t221 * t69;
t18 = t218 * t37 + t221 * t38;
t42 = t219 * t65 + t222 * t66;
t254 = -t141 + t212;
t111 = -t229 + ((-pkin(2) * qJD(2) + t301) * qJD(2) + t228) * t213;
t70 = t111 * t213 + t215 * t112;
t176 = 0.2e1 * t265;
t246 = t176 + t252;
t28 = -pkin(5) * t304 + t246;
t29 = -t156 * pkin(5) - qJ(6) * t304 + t131 * t174 + qJDD(6) + t37;
t251 = -pkin(5) * t29 + qJ(6) * t28;
t250 = -pkin(5) * t311 - qJ(6) * t93;
t249 = t172 * t274 - t221 * t241;
t10 = t18 * t219 - t222 * t56;
t11 = t18 * t222 + t219 * t56;
t3 = -t10 * t213 + t11 * t215;
t17 = t218 * t38 - t221 * t37;
t239 = (-t172 * t218 - t174 * t221) * t184;
t235 = t241 * pkin(5) - t337 + t56;
t234 = 0.2e1 * qJD(6) * t174 - t235;
t91 = t221 * t237 - t142;
t233 = t213 * (t222 * t91 + t260) + t215 * (t219 * t91 - t259);
t230 = pkin(5) * t313 + qJ(6) * t310 - t29;
t205 = t209 * qJDD(2);
t204 = t208 * qJDD(2);
t194 = t205 + t204;
t193 = t215 * t195;
t192 = t213 * t195;
t179 = -t186 - t302;
t178 = -t186 + t302;
t177 = t185 - t302;
t165 = t187 - 0.2e1 * t264;
t163 = t245 + 0.2e1 * t263;
t158 = -t302 - t185;
t135 = -t185 - t186;
t128 = -t179 * t219 - t277;
t127 = t179 * t222 - t278;
t117 = t187 * t219 - t222 * t245;
t116 = -t187 * t222 - t219 * t245;
t115 = t158 * t222 - t335;
t114 = t158 * t219 + t329;
t100 = (qJD(5) + t184) * t172 + t243;
t95 = (-qJD(5) + t184) * t174 + t253;
t90 = t174 * t273 + t218 * t237;
t85 = -t127 * t213 + t128 * t215;
t76 = -t116 * t213 + t117 * t215;
t71 = -t114 * t213 + t115 * t215;
t62 = -t221 * t93 + t296;
t61 = t221 * t95 + t296;
t59 = -t218 * t93 - t291;
t58 = t218 * t95 - t291;
t53 = -t100 * t219 + t365;
t51 = t100 * t222 + t366;
t49 = -t219 * t312 - t365;
t47 = t222 * t312 - t366;
t46 = t222 * t62 - t333;
t45 = t222 * t61 - t333;
t44 = t219 * t62 + t327;
t43 = t219 * t61 + t327;
t40 = t293 - t368;
t39 = t298 - t346;
t35 = -pkin(4) * t59 - t250;
t34 = (pkin(5) * t184 - 0.2e1 * qJD(6)) * t174 + t235;
t32 = -t213 * t51 + t215 * t53;
t30 = -t213 * t47 + t215 * t49;
t27 = t38 - t369;
t26 = t37 - t348;
t25 = (-t318 - t149) * pkin(5) + t234;
t24 = -pkin(5) * t149 + t234 + t337;
t23 = -t213 * t44 + t215 * t46;
t22 = -t213 * t43 + t215 * t45;
t21 = qJ(6) * t315 + t29;
t20 = (t315 - t304) * pkin(5) + t246;
t19 = t215 * t42 - t299;
t16 = -t230 - t348;
t15 = -t218 * t25 - t287 * t318 - t346;
t14 = -pkin(5) * t295 + t221 * t24 + t368;
t13 = -0.2e1 * t265 - t308 + t369;
t12 = -pkin(9) * t58 - t17;
t9 = t218 * t29 + t221 * t28;
t8 = t218 * t28 - t221 * t29;
t7 = -pkin(9) * t59 - t20 * t218 + t21 * t221;
t6 = t219 * t34 + t222 * t9;
t5 = t219 * t9 - t222 * t34;
t4 = -pkin(9) * t8 + (pkin(5) * t218 - t287) * t34;
t2 = -pkin(4) * t8 - t251;
t1 = -t213 * t5 + t215 * t6;
t31 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t269, 0, 0, 0, 0, 0, 0, (qJDD(2) * t223 - t303 * t220) * t214, (-qJDD(2) * t220 - t303 * t223) * t214, 0, t216 ^ 2 * t269 + (t220 * t153 - t223 * t247 - t236) * t214, 0, 0, 0, 0, 0, 0, (-t193 * t220 + t223 * t261) * t214, (t192 * t220 - t223 * t262) * t214, (t194 * t220 + t195 * t223) * t214, t216 * (-t111 * t215 + t112 * t213) + (-t141 * t223 + t220 * t70) * t214, 0, 0, 0, 0, 0, 0, t216 * (t114 * t215 + t115 * t213) + (-t163 * t223 + t220 * t71) * t214, t216 * (t127 * t215 + t128 * t213) + (-t165 * t223 + t220 * t85) * t214, t216 * (t116 * t215 + t117 * t213) + (-t135 * t223 + t220 * t76) * t214, t216 * (t213 * t42 + t215 * t41) + (-t126 * t223 + t19 * t220) * t214, 0, 0, 0, 0, 0, 0, t357, t216 * (t213 * t53 + t215 * t51) + (t220 * t32 - t364) * t214, t216 * (t213 * t45 + t215 * t43) + (t22 * t220 - t223 * t58) * t214, t216 * (t10 * t215 + t11 * t213) + (-t17 * t223 + t220 * t3) * t214, 0, 0, 0, 0, 0, 0, t357, t216 * (t213 * t46 + t215 * t44) + (t220 * t23 - t223 * t59) * t214, t216 * (t213 * t49 + t215 * t47) + (t220 * t30 + t364) * t214, t216 * (t213 * t6 + t215 * t5) + (t1 * t220 - t223 * t8) * t214; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t247, -t153, 0, 0, t204, 0.2e1 * t213 * t261, 0, t205, 0, 0, -qJ(3) * t193 + t254 * t215, qJ(3) * t192 - t213 * t254, pkin(2) * t195 + qJ(3) * t194 + t70, -pkin(2) * t141 + qJ(3) * t70, t213 * (t166 * t222 - t219 * t263) + t215 * (t166 * t219 + t222 * t263), t213 * (-t163 * t222 - t165 * t219) + t215 * (-t163 * t219 + t165 * t222), t213 * (-t178 * t219 + t329) + t215 * (t178 * t222 + t335), t213 * (-t164 * t219 + t222 * t264) + t215 * (t164 * t222 + t219 * t264), t213 * (t177 * t222 - t278) + t215 * (t177 * t219 + t277), (t213 * (-t188 * t222 + t190 * t219) + t215 * (-t188 * t219 - t190 * t222)) * qJD(4), t213 * (-pkin(8) * t114 + t282) + t215 * (-pkin(3) * t163 + pkin(8) * t115 - t281) - pkin(2) * t163 + qJ(3) * t71, t213 * (-pkin(8) * t127 + t281) + t215 * (-pkin(3) * t165 + pkin(8) * t128 + t282) - pkin(2) * t165 + qJ(3) * t85, t213 * (-pkin(8) * t116 - t41) + t215 * (-pkin(3) * t135 + pkin(8) * t117 + t42) - pkin(2) * t135 + qJ(3) * t76, -pkin(8) * t299 + t215 * (-pkin(3) * t126 + pkin(8) * t42) - pkin(2) * t126 + qJ(3) * t19, t233, -t361, t351, t307, -t372, t306, t213 * (-t219 * t26 + t222 * t39 - t356) + t215 * (t219 * t39 + t222 * t26 + t354) + t358, t213 * (-pkin(8) * t51 - t219 * t27 + t222 * t40) + t215 * (pkin(8) * t53 + t219 * t40 + t222 * t27 - t370) - t371 + qJ(3) * t32, t213 * (-pkin(8) * t43 + t12 * t222) + t215 * (pkin(8) * t45 + t12 * t219) + qJ(3) * t22 + (pkin(4) * t272 + t215 * t257 - pkin(2)) * t58, (t213 * (pkin(4) * t219 - pkin(9) * t222) + t215 * (-pkin(9) * t219 + t257) - pkin(2)) * t17 + t320 * t3, t233, t351, t361, t306, t372, t307, t213 * (t15 * t222 - t16 * t219 - t356) + t215 * (t15 * t219 + t16 * t222 + t354) + t358, t213 * (-pkin(8) * t44 - t219 * t35 + t222 * t7) + t215 * (-pkin(3) * t59 + pkin(8) * t46 + t219 * t7 + t222 * t35) - pkin(2) * t59 + qJ(3) * t23, t213 * (-pkin(8) * t47 - t13 * t219 + t14 * t222) + t215 * (pkin(8) * t49 + t13 * t222 + t14 * t219 + t370) + t371 + qJ(3) * t30, t213 * (-pkin(8) * t5 - t2 * t219 + t222 * t4) + t215 * (-pkin(3) * t8 + pkin(8) * t6 + t2 * t222 + t219 * t4) - pkin(2) * t8 + qJ(3) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t261, t262, -t195, t141, 0, 0, 0, 0, 0, 0, t163, t165, t135, t126, 0, 0, 0, 0, 0, 0, t323, t77, t58, t17, 0, 0, 0, 0, 0, 0, t323, t59, -t77, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t167, t186 - t185, t187, -t167, -t245, qJDD(4), -t65, -t66, 0, 0, t90, t321, t342, t249, -t355, t239, -pkin(4) * t318 - t293 + t347, pkin(4) * t100 + t298 + t367, pkin(9) * t61 + t18 + t338, -pkin(4) * t56 + pkin(9) * t18, t90, t342, -t321, t239, t355, t249, t221 * t25 + t256 * t318 + t347, pkin(9) * t62 + t20 * t221 + t21 * t218 + t338, -t367 + t218 * t24 + (pkin(4) + t300) * t312, pkin(9) * t9 + (t256 - t300) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t275, t316, t311, -t275, -t93, t156, -t37, -t38, 0, 0, t275, t311, -t316, t156, t93, -t275, t230, t250, t176 + t308, t251; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t313, t311, t317, t29;];
tauJ_reg  = t31;
