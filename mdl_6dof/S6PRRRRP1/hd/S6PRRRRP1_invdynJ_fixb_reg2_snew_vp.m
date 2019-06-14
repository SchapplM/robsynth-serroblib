% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRRRRP1
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
% Datum: 2019-05-05 09:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRRRRP1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP1_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:31:35
% EndTime: 2019-05-05 09:31:45
% DurationCPUTime: 4.41s
% Computational Cost: add. (20487->405), mult. (40567->549), div. (0->0), fcn. (30012->12), ass. (0->265)
t236 = sin(qJ(3));
t239 = cos(qJ(4));
t235 = sin(qJ(4));
t240 = cos(qJ(3));
t289 = t240 * t235;
t201 = (t236 * t239 + t289) * qJD(2);
t234 = sin(qJ(5));
t238 = cos(qJ(5));
t280 = qJD(3) + qJD(4);
t184 = t201 * t234 - t238 * t280;
t186 = t238 * t201 + t234 * t280;
t156 = t186 * t184;
t282 = qJD(2) * qJD(3);
t269 = t240 * t282;
t281 = t236 * qJDD(2);
t205 = t269 + t281;
t225 = t240 * qJDD(2);
t270 = t236 * t282;
t206 = t225 - t270;
t263 = t235 * t205 - t239 * t206;
t159 = -t201 * qJD(4) - t263;
t158 = qJDD(5) - t159;
t331 = -t156 + t158;
t339 = pkin(5) * t331;
t286 = qJD(2) * t236;
t199 = -t239 * t240 * qJD(2) + t235 * t286;
t160 = -t199 * qJD(4) + t205 * t239 + t206 * t235;
t279 = qJDD(3) + qJDD(4);
t135 = -t184 * qJD(5) + t238 * t160 + t234 * t279;
t195 = qJD(5) + t199;
t167 = t195 * t184;
t122 = t167 + t135;
t338 = qJ(6) * t122;
t177 = t201 * t199;
t329 = -t177 + t279;
t337 = t235 * t329;
t336 = t239 * t329;
t301 = t331 * t234;
t300 = t331 * t238;
t193 = t280 * t199;
t335 = t160 - t193;
t231 = sin(pkin(6));
t232 = cos(pkin(6));
t305 = sin(pkin(11));
t306 = cos(pkin(11));
t250 = t305 * g(1) - t306 * g(2);
t249 = t232 * t250;
t288 = -g(3) + qJDD(1);
t334 = t231 * t288 + t249;
t183 = t186 ^ 2;
t194 = t195 ^ 2;
t153 = -t183 - t194;
t182 = t184 ^ 2;
t264 = t160 * t234 - t238 * t279;
t134 = -qJD(5) * t186 - t264;
t162 = pkin(5) * t195 - qJ(6) * t186;
t175 = pkin(4) * t199 - pkin(10) * t201;
t276 = t280 ^ 2;
t210 = -t306 * g(1) - t305 * g(2);
t237 = sin(qJ(2));
t241 = cos(qJ(2));
t170 = t241 * t210 + t237 * t334;
t243 = qJD(2) ^ 2;
t246 = -t243 * pkin(2) + qJDD(2) * pkin(8) + t170;
t245 = t236 * t246;
t247 = -t231 * t250 + t232 * t288;
t143 = t236 * t247 + t240 * t246;
t214 = qJD(3) * pkin(3) - pkin(9) * t286;
t229 = t240 ^ 2;
t227 = t229 * t243;
t125 = -pkin(3) * t227 + t206 * pkin(9) - qJD(3) * t214 + t143;
t291 = t239 * t125;
t292 = t236 * t243;
t304 = qJDD(3) * pkin(3);
t317 = t205 * pkin(9);
t94 = t291 + t235 * (-t245 + t304 - t317) + (pkin(3) * t292 + pkin(9) * t282 + t247) * t289;
t80 = -t276 * pkin(4) + t279 * pkin(10) - t199 * t175 + t94;
t261 = t237 * t210 - t241 * t334;
t163 = -qJDD(2) * pkin(2) - t243 * pkin(8) + t261;
t140 = -t206 * pkin(3) - pkin(9) * t227 + t214 * t286 + t163;
t262 = t201 * t280;
t92 = -t335 * pkin(10) + (-t159 + t262) * pkin(4) + t140;
t51 = t234 * t92 + t238 * t80;
t252 = t134 * qJ(6) - 0.2e1 * qJD(6) * t184 - t162 * t195 + t51;
t333 = -t252 + (t153 + t182) * pkin(5);
t50 = t234 * t80 - t238 * t92;
t22 = t234 * t50 + t238 * t51;
t330 = -t167 + t135;
t119 = (qJD(5) - t195) * t186 + t264;
t197 = t199 ^ 2;
t198 = t201 ^ 2;
t144 = -t194 - t182;
t97 = t144 * t234 + t300;
t328 = pkin(4) * t97;
t284 = qJD(6) * t186;
t179 = -0.2e1 * t284;
t251 = -t338 - t50 + t339;
t31 = t179 + t251;
t327 = pkin(5) * t31;
t137 = -t182 - t183;
t85 = -t119 * t238 + t122 * t234;
t63 = -t137 * t239 + t235 * t85;
t326 = pkin(9) * t63;
t118 = (qJD(5) + t195) * t186 + t264;
t98 = t144 * t238 - t301;
t67 = -t118 * t239 + t235 * t98;
t325 = pkin(9) * t67;
t127 = t156 + t158;
t302 = t127 * t238;
t102 = -t153 * t234 - t302;
t71 = t102 * t235 - t239 * t330;
t324 = pkin(9) * t71;
t83 = -t119 * t234 - t122 * t238;
t323 = pkin(10) * t83;
t322 = pkin(10) * t97;
t303 = t127 * t234;
t101 = t153 * t238 - t303;
t321 = pkin(4) * t101;
t320 = pkin(4) * t235;
t319 = pkin(5) * t122;
t318 = pkin(10) * t101;
t64 = t137 * t235 + t239 * t85;
t33 = -t236 * t63 + t240 * t64;
t316 = -pkin(2) * t83 + pkin(8) * t33;
t68 = t118 * t235 + t239 * t98;
t41 = -t236 * t67 + t240 * t68;
t315 = -pkin(2) * t97 + pkin(8) * t41;
t217 = t240 * t292;
t142 = -t240 * t247 + t245;
t244 = pkin(9) * t269 - t142 - t317;
t93 = t125 * t235 - t239 * (pkin(3) * t217 + t244 + t304);
t79 = -t279 * pkin(4) - t276 * pkin(10) + t175 * t201 + t93;
t314 = -pkin(4) * t79 + pkin(10) * t22;
t313 = t234 * t31;
t75 = t234 * t79;
t56 = t235 * t94 - t239 * t93;
t312 = t236 * t56;
t311 = t238 * t31;
t76 = t238 * t79;
t72 = t102 * t239 + t235 * t330;
t43 = -t236 * t71 + t240 * t72;
t310 = -pkin(2) * t101 + pkin(8) * t43;
t309 = -pkin(4) * t137 + pkin(10) * t85;
t308 = -pkin(4) * t118 + pkin(10) * t98;
t307 = -pkin(4) * t330 + pkin(10) * t102;
t299 = t140 * t235;
t298 = t140 * t239;
t173 = t177 + t279;
t297 = t173 * t235;
t296 = t173 * t239;
t295 = t195 * t234;
t294 = t195 * t238;
t211 = qJDD(3) + t217;
t293 = t236 * t211;
t212 = qJDD(3) - t217;
t290 = t240 * t212;
t278 = t75 + t307;
t277 = -t76 + t308;
t275 = t235 * t156;
t274 = t239 * t156;
t273 = -pkin(4) * t239 - pkin(3);
t272 = -pkin(3) * t97 + pkin(9) * t68;
t271 = -pkin(3) * t101 + pkin(9) * t72;
t57 = t235 * t93 + t239 * t94;
t27 = -qJ(6) * t119 + (-t137 - t182) * pkin(5) + t252;
t180 = 0.2e1 * t284;
t29 = t180 - t251 + t338;
t267 = t234 * t29 + t238 * t27 + t309;
t53 = -t134 * pkin(5) - t182 * qJ(6) + t162 * t186 + qJDD(6) + t79;
t47 = -qJ(6) * t153 + t53;
t88 = -pkin(5) * t330 - qJ(6) * t127;
t266 = t234 * t47 + t238 * t88 + t307;
t265 = t309 + t22;
t108 = t142 * t236 + t240 * t143;
t260 = t235 * t262;
t259 = t235 * t193;
t258 = t239 * t262;
t257 = t239 * t193;
t14 = t22 * t235 - t239 * t79;
t15 = t22 * t239 + t235 * t79;
t3 = -t236 * t14 + t240 * t15;
t21 = t234 * t51 - t238 * t50;
t34 = -pkin(5) * t182 + t252;
t12 = t238 * t34 - t313;
t18 = -pkin(5) * t53 + qJ(6) * t34;
t255 = -pkin(4) * t53 + pkin(10) * t12 - qJ(6) * t313 + t238 * t18;
t254 = t201 * qJD(3) - t263;
t38 = -pkin(5) * t118 + qJ(6) * t144 - t53;
t253 = -qJ(6) * t301 + t238 * t38 + t308;
t248 = t251 + t339;
t242 = qJD(3) ^ 2;
t228 = t236 ^ 2;
t226 = t228 * t243;
t216 = -t227 - t242;
t215 = -t226 - t242;
t209 = t226 + t227;
t208 = (t228 + t229) * qJDD(2);
t207 = t225 - 0.2e1 * t270;
t204 = 0.2e1 * t269 + t281;
t191 = -t198 + t276;
t190 = t197 - t276;
t189 = -t198 - t276;
t188 = -t215 * t236 - t290;
t187 = t216 * t240 - t293;
t176 = t198 - t197;
t171 = -t276 - t197;
t165 = -t183 + t194;
t164 = t182 - t194;
t161 = -t197 - t198;
t154 = t183 - t182;
t152 = -t189 * t235 - t296;
t151 = t189 * t239 - t297;
t150 = t160 + t193;
t145 = (0.2e1 * qJD(4) + qJD(3)) * t201 + t263;
t139 = t171 * t239 - t337;
t138 = t171 * t235 + t336;
t130 = (-t184 * t238 + t186 * t234) * t195;
t129 = (-t184 * t234 - t186 * t238) * t195;
t115 = t135 * t238 - t186 * t295;
t114 = t135 * t234 + t186 * t294;
t113 = -t134 * t234 + t184 * t294;
t112 = t134 * t238 + t184 * t295;
t111 = -t236 * t151 + t240 * t152;
t110 = t150 * t235 + t239 * t254;
t109 = -t150 * t239 + t235 * t254;
t107 = t164 * t238 - t303;
t106 = -t165 * t234 + t300;
t105 = t164 * t234 + t302;
t104 = t165 * t238 + t301;
t103 = -t236 * t138 + t240 * t139;
t86 = -t118 * t238 - t234 * t330;
t84 = -t118 * t234 + t238 * t330;
t74 = -t236 * t109 + t240 * t110;
t73 = t236 * (t130 * t239 + t158 * t235) + t240 * (t130 * t235 - t158 * t239);
t70 = pkin(3) * t71;
t66 = pkin(3) * t67;
t62 = pkin(3) * t63;
t61 = pkin(9) * t64;
t60 = t236 * (t115 * t239 + t275) + t240 * (t115 * t235 - t274);
t59 = t236 * (t113 * t239 - t275) + t240 * (t113 * t235 + t274);
t58 = -pkin(4) * t83 + t319;
t55 = t76 - t318;
t54 = t75 - t322;
t45 = t236 * (t107 * t239 - t119 * t235) + t240 * (t107 * t235 + t119 * t239);
t44 = t236 * (t106 * t239 + t122 * t235) + t240 * (t106 * t235 - t122 * t239);
t39 = t236 * (t154 * t235 + t239 * t86) + t240 * (-t154 * t239 + t235 * t86);
t37 = t51 - t321;
t35 = t50 - t328;
t30 = t240 * t57 - t312;
t25 = -t321 - t333;
t24 = -t234 * t88 + t238 * t47 - t318;
t23 = -qJ(6) * t300 - t234 * t38 - t322;
t19 = t180 - t248 - t328;
t16 = -t21 - t323;
t13 = t232 * (t236 * t72 + t240 * t71) + (-t241 * t101 + t237 * t43) * t231;
t11 = t234 * t34 + t311;
t10 = t232 * (t236 * t68 + t240 * t67) + (t237 * t41 - t241 * t97) * t231;
t8 = t232 * (t236 * t64 + t240 * t63) + (t237 * t33 - t241 * t83) * t231;
t7 = t12 * t239 + t235 * t53;
t6 = t12 * t235 - t239 * t53;
t5 = -t234 * t27 + t238 * t29 - t323;
t4 = -pkin(4) * t11 - t327;
t2 = -pkin(10) * t11 - qJ(6) * t311 - t18 * t234;
t1 = -t236 * t6 + t240 * t7;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t288, 0, 0, 0, 0, 0, 0, (qJDD(2) * t241 - t237 * t243) * t231, (-qJDD(2) * t237 - t241 * t243) * t231, 0, t232 ^ 2 * t288 + (t237 * t170 - t241 * t261 - t249) * t231, 0, 0, 0, 0, 0, 0, t232 * (t211 * t240 + t216 * t236) + (t237 * t187 + t241 * t207) * t231, t232 * (-t212 * t236 + t215 * t240) + (t237 * t188 - t241 * t204) * t231, (t208 * t237 + t209 * t241) * t231, t232 * (-t142 * t240 + t143 * t236) + (t237 * t108 - t241 * t163) * t231, 0, 0, 0, 0, 0, 0, t232 * (t138 * t240 + t139 * t236) + (t237 * t103 - t241 * t145) * t231, t232 * (t151 * t240 + t152 * t236) + (t237 * t111 - t241 * t335) * t231, t232 * (t109 * t240 + t110 * t236) + (-t241 * t161 + t237 * t74) * t231, t232 * (t236 * t57 + t240 * t56) + (-t241 * t140 + t237 * t30) * t231, 0, 0, 0, 0, 0, 0, t10, t13, t8, t232 * (t14 * t240 + t15 * t236) + (-t241 * t21 + t237 * t3) * t231, 0, 0, 0, 0, 0, 0, t10, t13, t8, t232 * (t236 * t7 + t240 * t6) + (t237 * t1 - t241 * t11) * t231; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t261, -t170, 0, 0, (t205 + t269) * t236, t204 * t240 + t207 * t236, t293 + t240 * (-t226 + t242), (t206 - t270) * t240, t236 * (t227 - t242) + t290, 0, pkin(2) * t207 + pkin(8) * t187 - t163 * t240, -pkin(2) * t204 + pkin(8) * t188 + t163 * t236, pkin(2) * t209 + pkin(8) * t208 + t108, -pkin(2) * t163 + pkin(8) * t108, t236 * (t239 * t160 - t260) + t240 * (t235 * t160 + t258), t236 * (-t145 * t239 - t235 * t335) + t240 * (-t145 * t235 + t239 * t335), t236 * (-t191 * t235 + t336) + t240 * (t191 * t239 + t337), t236 * (-t235 * t159 + t257) + t240 * (t239 * t159 + t259), t236 * (t190 * t239 - t297) + t240 * (t190 * t235 + t296), t236 * (-t257 + t260) + t240 * (-t259 - t258), t236 * (-pkin(9) * t138 + t299) + t240 * (-pkin(3) * t145 + pkin(9) * t139 - t298) - pkin(2) * t145 + pkin(8) * t103, t236 * (-pkin(9) * t151 + t298) + t240 * (-pkin(3) * t335 + pkin(9) * t152 + t299) - pkin(2) * t335 + pkin(8) * t111, t236 * (-pkin(9) * t109 - t56) + t240 * (-pkin(3) * t161 + pkin(9) * t110 + t57) - pkin(2) * t161 + pkin(8) * t74, -pkin(9) * t312 + t240 * (-pkin(3) * t140 + pkin(9) * t57) - pkin(2) * t140 + pkin(8) * t30, t60, t39, t44, t59, t45, t73, t236 * (-t235 * t35 + t239 * t54 - t325) + t240 * (t235 * t54 + t239 * t35 + t272) + t315, t236 * (-t235 * t37 + t239 * t55 - t324) + t240 * (t235 * t55 + t239 * t37 + t271) + t310, t236 * (t16 * t239 + t83 * t320 - t326) + t240 * (t16 * t235 + t273 * t83 + t61) + t316, (t236 * (-pkin(10) * t239 + t320) + t240 * (-pkin(10) * t235 + t273) - pkin(2)) * t21 + (pkin(8) + pkin(9)) * t3, t60, t39, t44, t59, t45, t73, t236 * (-t19 * t235 + t23 * t239 - t325) + t240 * (t19 * t239 + t23 * t235 + t272) + t315, t236 * (-t235 * t25 + t239 * t24 - t324) + t240 * (t235 * t24 + t239 * t25 + t271) + t310, t236 * (-t235 * t58 + t239 * t5 - t326) + t240 * (-pkin(3) * t83 + t235 * t5 + t239 * t58 + t61) + t316, t236 * (-pkin(9) * t6 + t2 * t239 - t235 * t4) + t240 * (-pkin(3) * t11 + pkin(9) * t7 + t2 * t235 + t239 * t4) - pkin(2) * t11 + pkin(8) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t217, t226 - t227, t281, t217, t225, qJDD(3), -t142, -t143, 0, 0, t177, t176, t150, -t177, t254, t279, pkin(3) * t138 - t93, -t291 - t235 * t244 + (-t211 * t235 + t151) * pkin(3), pkin(3) * t109, pkin(3) * t56, t114, t84, t104, t112, t105, t129, t66 + t277, t70 + t278, t62 + t265, pkin(3) * t14 + t314, t114, t84, t104, t112, t105, t129, t253 + t66, t70 + t266, t62 + t267, pkin(3) * t6 + t255; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t177, t176, t150, -t177, t254, t279, -t93, -t94, 0, 0, t114, t84, t104, t112, t105, t129, t277, t278, t265, t314, t114, t84, t104, t112, t105, t129, t253, t266, t267, t255; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t156, t154, t122, -t156, -t119, t158, -t50, -t51, 0, 0, t156, t154, t122, -t156, -t119, t158, t179 + t248, t333, -t319, t327; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, t330, t137, t53;];
tauJ_reg  = t9;
