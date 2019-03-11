% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRPP5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRPP5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP5_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:44:45
% EndTime: 2019-03-09 04:44:55
% DurationCPUTime: 4.35s
% Computational Cost: add. (6415->503), mult. (14864->577), div. (0->0), fcn. (11118->10), ass. (0->242)
t188 = sin(pkin(9));
t192 = sin(qJ(3));
t189 = cos(pkin(9));
t334 = cos(qJ(3));
t264 = t334 * t189;
t220 = -t192 * t188 + t264;
t362 = t220 * qJD(1);
t118 = qJD(4) - t362;
t139 = t334 * t188 + t192 * t189;
t133 = t139 * qJD(3);
t271 = t188 * qJDD(1);
t240 = -qJDD(1) * t264 + t192 * t271;
t83 = qJD(1) * t133 + t240;
t77 = qJDD(4) + t83;
t361 = t77 * qJ(5) + t118 * qJD(5);
t327 = pkin(7) + qJ(2);
t149 = t327 * t188;
t140 = qJD(1) * t149;
t150 = t327 * t189;
t141 = qJD(1) * t150;
t85 = -t192 * t140 + t334 * t141;
t360 = qJD(3) * t85;
t186 = pkin(9) + qJ(3);
t178 = sin(t186);
t195 = cos(qJ(1));
t295 = t178 * t195;
t193 = sin(qJ(1));
t297 = t178 * t193;
t348 = g(1) * t295 + g(2) * t297;
t303 = qJDD(1) * pkin(1);
t177 = qJDD(2) - t303;
t346 = g(1) * t193 - g(2) * t195;
t231 = -t177 + t346;
t132 = t220 * qJD(3);
t212 = t139 * qJDD(1);
t200 = qJD(1) * t132 + t212;
t359 = -qJD(3) * qJD(4) - t200;
t345 = t139 * qJD(1);
t358 = qJD(3) * t345;
t194 = cos(qJ(4));
t191 = sin(qJ(4));
t291 = t191 * qJ(5);
t337 = pkin(4) + pkin(5);
t227 = -t337 * t194 - t291;
t136 = pkin(3) - t227;
t276 = qJD(4) * t194;
t245 = -t194 * qJDD(3) + t276 * t345;
t278 = qJD(3) * t191;
t48 = t191 * t212 + (qJD(4) + t362) * t278 + t245;
t95 = -t194 * qJD(3) + t191 * t345;
t357 = t48 * qJ(6) + t95 * qJD(6);
t108 = t118 * qJ(5);
t173 = t189 * pkin(2) + pkin(1);
t147 = -t173 * qJD(1) + qJD(2);
t57 = -pkin(3) * t362 - pkin(8) * t345 + t147;
t79 = qJD(3) * pkin(8) + t85;
t29 = t191 * t57 + t194 * t79;
t23 = t108 + t29;
t320 = t118 * t23;
t277 = qJD(4) * t191;
t273 = qJD(1) * qJD(2);
t339 = t327 * qJDD(1) + t273;
t105 = t339 * t188;
t106 = t339 * t189;
t222 = -t192 * t105 + t334 * t106;
t350 = -t334 * t140 - t192 * t141;
t37 = qJDD(3) * pkin(8) + qJD(3) * t350 + t222;
t270 = t189 * qJDD(1);
t39 = -pkin(2) * t270 + t83 * pkin(3) - pkin(8) * t200 + t177;
t259 = t191 * t37 - t194 * t39 + t79 * t276 + t57 * t277;
t73 = t77 * pkin(4);
t344 = t73 - qJDD(5);
t5 = t259 - t344;
t356 = -t5 + t320;
t355 = 0.2e1 * t361;
t113 = t118 ^ 2;
t97 = t194 * t345 + t278;
t94 = t97 ^ 2;
t354 = -t113 - t94;
t353 = t191 * qJD(5) + t85;
t318 = t118 * t95;
t47 = -t191 * qJDD(3) + t359 * t194 + t345 * t277;
t210 = t47 - t318;
t91 = -t192 * t149 + t334 * t150;
t352 = t91 * qJD(3);
t28 = -t191 * t79 + t194 * t57;
t283 = qJD(5) - t28;
t351 = t133 * qJ(5) - qJD(5) * t220;
t349 = -t334 * t149 - t192 * t150;
t179 = cos(t186);
t347 = -t179 * pkin(3) - t178 * pkin(8);
t302 = t362 * t191;
t65 = t194 * t77;
t235 = t65 + (-t277 + t302) * t118;
t317 = t345 * t95;
t208 = t235 - t317;
t247 = g(1) * t195 + g(2) * t193;
t343 = qJ(2) * qJDD(1);
t18 = t97 * qJ(6) + t28;
t284 = qJD(5) - t18;
t12 = -t337 * t118 + t284;
t225 = t191 * t39 + t194 * t37 + t57 * t276 - t79 * t277;
t4 = t225 + t361;
t2 = t4 + t357;
t342 = t118 * t12 + t2;
t250 = qJD(3) * pkin(3) + t350;
t216 = t97 * qJ(5) + t250;
t30 = pkin(4) * t95 - t216;
t336 = pkin(8) * t77;
t341 = t118 * t30 - t336;
t17 = -t337 * t95 + qJD(6) + t216;
t286 = t194 * t195;
t290 = t191 * t193;
t119 = t179 * t290 + t286;
t285 = t195 * t191;
t287 = t193 * t194;
t121 = t179 * t285 - t287;
t298 = t178 * t191;
t209 = g(1) * t121 + g(2) * t119 + g(3) * t298 - t259;
t207 = t209 + t344;
t313 = t47 * qJ(6);
t340 = (qJD(6) + t17) * t97 + t207 - t313;
t338 = t95 ^ 2;
t335 = t77 * pkin(5);
t331 = g(2) * t327;
t329 = g(3) * t179;
t328 = t97 * t95;
t326 = pkin(8) - qJ(6);
t325 = -t191 * t48 - t95 * t276;
t80 = pkin(3) * t345 - pkin(8) * t362;
t324 = t191 * t80 + t194 * t350;
t82 = -pkin(3) * t220 - pkin(8) * t139 - t173;
t323 = t191 * t82 + t194 * t91;
t322 = pkin(8) * qJD(4);
t321 = qJ(5) * t48;
t319 = t118 * t29;
t19 = qJ(6) * t95 + t29;
t14 = t108 + t19;
t316 = t14 * t118;
t63 = t191 * t77;
t315 = t191 * t97;
t314 = t194 * t95;
t312 = t47 * t191;
t311 = t95 * qJ(5);
t310 = t97 * t118;
t309 = t97 * t345;
t267 = t337 * t191;
t304 = qJ(5) * t194;
t228 = -t267 + t304;
t308 = t118 * t228 + t353;
t25 = qJ(5) * t345 + t324;
t307 = -qJ(6) * t302 - t194 * qJD(6) - t326 * t277 - t25;
t241 = pkin(4) * t191 - t304;
t306 = t118 * t241 - t353;
t152 = t326 * t194;
t75 = t191 * t350;
t305 = qJD(4) * t152 - t191 * qJD(6) - t75 - (-qJ(6) * t362 - t80) * t194 + t337 * t345;
t255 = t118 * t194;
t301 = t132 * t191;
t300 = t132 * t194;
t299 = t139 * t191;
t296 = t178 * t194;
t294 = t179 * t193;
t293 = t179 * t194;
t292 = t179 * t195;
t281 = t348 * t191;
t280 = t348 * t194;
t279 = t188 ^ 2 + t189 ^ 2;
t275 = qJD(5) * t194;
t58 = t220 * qJD(2) + qJD(3) * t349;
t269 = t191 * t58 + t91 * t276 + t82 * t277;
t81 = pkin(3) * t133 - pkin(8) * t132;
t268 = t191 * t81 + t194 * t58 + t82 * t276;
t35 = -qJ(5) * t220 + t323;
t266 = g(1) * t292 + g(2) * t294 + g(3) * t178;
t254 = t334 * t105 + t192 * t106 + t360;
t38 = -qJDD(3) * pkin(3) + t254;
t6 = t48 * pkin(4) + t47 * qJ(5) - t97 * qJD(5) + t38;
t3 = -pkin(5) * t48 + qJDD(6) - t6;
t265 = t3 - t329;
t263 = t139 * t277;
t262 = t139 * t276;
t87 = t191 * t91;
t260 = t194 * t82 - t87;
t258 = t279 * qJD(1) ^ 2;
t120 = t179 * t287 - t285;
t257 = -t119 * pkin(4) + qJ(5) * t120;
t122 = t179 * t286 + t290;
t256 = -t121 * pkin(4) + qJ(5) * t122;
t253 = pkin(4) * t293 + t179 * t291 - t347;
t252 = 0.2e1 * t279;
t251 = -g(1) * t297 + g(2) * t295;
t249 = g(1) * t119 - g(2) * t121;
t248 = g(1) * t120 - g(2) * t122;
t1 = -t97 * qJD(6) + t313 - t335 + t5;
t244 = -t1 + t316;
t22 = -pkin(4) * t118 + t283;
t243 = t118 * t22 + t4;
t242 = pkin(4) * t194 + t291;
t239 = -t191 * t23 + t194 * t22;
t238 = t314 + t315;
t237 = t194 * t81 - t269;
t236 = -qJ(6) * t132 - qJD(6) * t139;
t234 = t118 * t276 - t255 * t362 + t63;
t233 = pkin(3) + t242;
t232 = -t173 + t347;
t230 = t118 * t322 + t329;
t229 = -t194 * t47 - t97 * t277;
t226 = -t120 * pkin(4) - t119 * qJ(5) + t195 * t327;
t224 = -t91 * t277 + t268;
t223 = -t118 * t250 - t336;
t219 = -t230 - t6;
t218 = pkin(3) * t292 + t122 * pkin(4) + pkin(8) * t295 + qJ(5) * t121 + t195 * t173;
t215 = t234 + t309;
t211 = t231 + t303;
t206 = t252 * t273 - t247;
t205 = t30 * t97 - t207;
t203 = g(1) * t122 + g(2) * t120 + g(3) * t296 - t225;
t202 = t118 * t28 + t203;
t59 = t139 * qJD(2) + t352;
t201 = -qJDD(4) - t240 + t328 - t358;
t199 = t359 * t191 - t245;
t159 = pkin(8) * t292;
t156 = pkin(8) * t294;
t153 = qJ(5) * t296;
t151 = t326 * t191;
t146 = -t173 * qJDD(1) + qJDD(2);
t50 = pkin(4) * t97 + t311;
t49 = t139 * t241 - t349;
t41 = -t337 * t97 - t311;
t40 = t139 * t228 + t349;
t36 = pkin(4) * t220 - t260;
t27 = -pkin(4) * t345 - t194 * t80 + t75;
t24 = qJ(6) * t299 + t35;
t16 = t87 + (-qJ(6) * t139 - t82) * t194 + t337 * t220;
t13 = t241 * t132 + (qJD(4) * t242 - t275) * t139 + t59;
t11 = -t133 * pkin(4) - t237;
t10 = -t352 + t228 * t132 + (qJD(4) * t227 - qJD(2) + t275) * t139;
t9 = t224 + t351;
t8 = qJ(6) * t262 + (-qJD(4) * t91 - t236) * t191 + t268 + t351;
t7 = qJ(6) * t263 - t337 * t133 + (t236 - t81) * t194 + t269;
t15 = [qJDD(1), t346, t247, t211 * t189, -t211 * t188, t252 * t343 + t206, t231 * pkin(1) + (t279 * t343 + t206) * qJ(2), t132 * t345 + t139 * t200, t132 * t362 - t133 * t345 - t139 * t83 + t200 * t220, qJD(3) * t132 + qJDD(3) * t139, -qJD(3) * t133 + qJDD(3) * t220, 0, -t59 * qJD(3) + qJDD(3) * t349 + t147 * t133 - t146 * t220 - t173 * t83 + t179 * t346, -t58 * qJD(3) - t91 * qJDD(3) + t147 * t132 + t146 * t139 - t173 * t200 + t251, t139 * t229 + t97 * t300, -t238 * t132 + (t312 - t194 * t48 + (t191 * t95 - t194 * t97) * qJD(4)) * t139, t139 * t65 + t97 * t133 + t47 * t220 + (-t263 + t300) * t118, -t77 * t299 - t95 * t133 + t48 * t220 + (-t262 - t301) * t118, t118 * t133 - t220 * t77, t237 * t118 + t260 * t77 + t259 * t220 + t28 * t133 + t59 * t95 - t349 * t48 - t250 * t301 + (t38 * t191 - t250 * t276) * t139 + t248, -t224 * t118 - t323 * t77 + t225 * t220 - t29 * t133 + t59 * t97 + t349 * t47 - t250 * t300 + (t38 * t194 + t250 * t277) * t139 - t249, t30 * t301 - t11 * t118 + t13 * t95 - t22 * t133 + t5 * t220 - t36 * t77 + t49 * t48 + (t6 * t191 + t276 * t30) * t139 + t248, t11 * t97 - t35 * t48 - t36 * t47 - t9 * t95 + t239 * t132 + (-t4 * t191 + t5 * t194 + (-t191 * t22 - t194 * t23) * qJD(4)) * t139 - t251, -t30 * t300 + t9 * t118 - t13 * t97 + t23 * t133 - t4 * t220 + t35 * t77 + t49 * t47 + (-t6 * t194 + t277 * t30) * t139 + t249, t4 * t35 + t23 * t9 + t6 * t49 + t30 * t13 + t5 * t36 + t22 * t11 - g(1) * t226 - g(2) * t218 + (-g(1) * t232 - t331) * t193, -t17 * t301 + t1 * t220 - t10 * t95 - t7 * t118 - t12 * t133 - t16 * t77 - t40 * t48 + (-t17 * t276 - t3 * t191) * t139 + t248, t17 * t300 + t10 * t97 + t8 * t118 + t14 * t133 - t2 * t220 + t24 * t77 - t40 * t47 + (-t17 * t277 + t3 * t194) * t139 + t249, t16 * t47 + t24 * t48 - t7 * t97 + t8 * t95 + (-t12 * t194 + t14 * t191) * t132 + (-t1 * t194 + t2 * t191 + (t12 * t191 + t14 * t194) * qJD(4)) * t139 + t251, t2 * t24 + t14 * t8 + t1 * t16 + t12 * t7 + t3 * t40 + t17 * t10 - g(1) * (-t120 * pkin(5) + t226) - g(2) * (pkin(5) * t122 - qJ(6) * t295 + t218) + (-g(1) * (qJ(6) * t178 + t232) - t331) * t193; 0, 0, 0, -t270, t271, -t258, -qJ(2) * t258 - t231, 0, 0, 0, 0, 0, t240 + 0.2e1 * t358, 0.2e1 * t362 * qJD(3) + t212, 0, 0, 0, 0, 0, t208, -t113 * t194 - t309 - t63, t208 (t314 - t315) * t362 - t229 + t325, t215, t243 * t191 + t194 * t356 - t30 * t345 - t346, t208, t215, -t210 * t194 + (t48 - t310) * t191, t17 * t345 + t191 * t342 + t244 * t194 - t346; 0, 0, 0, 0, 0, 0, 0, -t345 * t362, t345 ^ 2 - t362 ^ 2, t212, -t240, qJDD(3), -t147 * t345 - t254 - t329 + t348 + t360, -t147 * t362 - t222 + t266, t255 * t97 - t312, t238 * t362 + t229 + t325, t234 - t309, t235 + t317, -t118 * t345, -pkin(3) * t48 + t75 * t118 - t28 * t345 - t85 * t95 + (-t329 - t38 + (-t80 - t322) * t118) * t194 + t223 * t191 + t280, pkin(3) * t47 + t324 * t118 + t29 * t345 - t85 * t97 + t223 * t194 + (t230 + t38) * t191 - t281, t27 * t118 + t191 * t341 + t219 * t194 + t22 * t345 - t233 * t48 + t306 * t95 + t280, t25 * t95 - t27 * t97 + ((qJD(4) * t97 - t48) * pkin(8) + t243) * t194 + ((qJD(4) * t95 - t47) * pkin(8) - t356) * t191 - t266, -t25 * t118 + t219 * t191 - t194 * t341 - t23 * t345 - t233 * t47 - t306 * t97 + t281, -t23 * t25 - t22 * t27 - g(1) * t159 - g(2) * t156 - g(3) * t253 + t306 * t30 + (qJD(4) * t239 + t5 * t191 + t4 * t194) * pkin(8) + (t178 * t247 - t6) * t233, t12 * t345 - t136 * t48 - t151 * t77 + t194 * t265 - t308 * t95 + t280 + (-t17 * t191 - t305) * t118, t118 * t307 - t136 * t47 - t14 * t345 + t152 * t77 + t17 * t255 + t191 * t265 + t308 * t97 + t281, t151 * t47 + t152 * t48 + t244 * t191 - t194 * t342 - t305 * t97 + t307 * t95 + t266, t2 * t152 + t1 * t151 + t3 * t136 - g(1) * (-qJ(6) * t292 + t159) - g(2) * (-qJ(6) * t294 + t156) - g(3) * (pkin(5) * t293 + t253) + t308 * t17 + t307 * t14 + t305 * t12 + (g(3) * qJ(6) + t136 * t247) * t178; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t328, t94 - t338, -t210, t199 + t310, t77, t250 * t97 + t209 + t319, -t250 * t95 + t202, -t50 * t95 - t205 + t319 + t73, pkin(4) * t47 - t321 + (t23 - t29) * t97 + (t22 - t283) * t95, -t30 * t95 + t50 * t97 - t202 + t355, t4 * qJ(5) - t5 * pkin(4) - t30 * t50 - t22 * t29 - g(1) * t256 - g(2) * t257 - g(3) * (-pkin(4) * t298 + t153) + t283 * t23, t19 * t118 + t41 * t95 + (pkin(5) + t337) * t77 + t340, -t118 * t18 + t17 * t95 - t41 * t97 - t203 + t355 + t357, t321 - t337 * t47 + (-t14 + t19) * t97 + (-t12 + t284) * t95, t2 * qJ(5) - t1 * t337 - t12 * t19 - t17 * t41 - g(1) * (-pkin(5) * t121 + t256) - g(2) * (-pkin(5) * t119 + t257) - g(3) * (-t178 * t267 + t153) + t284 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t201, -t210, t354, t205 - t320, t201, t354, t210, -t316 - t335 - t340; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t199 - t310, -t47 - t318, -t94 - t338, t12 * t97 - t14 * t95 + t265 + t348;];
tau_reg  = t15;
