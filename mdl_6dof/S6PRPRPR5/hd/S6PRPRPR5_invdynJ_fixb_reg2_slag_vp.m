% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6PRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPRPR5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_invdynJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:45:10
% EndTime: 2019-03-08 19:45:19
% DurationCPUTime: 5.50s
% Computational Cost: add. (5750->495), mult. (13422->637), div. (0->0), fcn. (10917->14), ass. (0->256)
t170 = sin(pkin(11));
t173 = cos(pkin(11));
t177 = sin(qJ(4));
t331 = cos(qJ(4));
t132 = t331 * t170 + t177 * t173;
t346 = t132 * qJD(2);
t354 = qJD(6) + t346;
t263 = t331 * t173;
t239 = qJD(2) * t263;
t280 = qJD(2) * t170;
t261 = t177 * t280;
t122 = -t239 + t261;
t176 = sin(qJ(6));
t179 = cos(qJ(6));
t89 = qJD(4) * t176 - t179 * t122;
t247 = t354 * t89;
t276 = qJD(6) * t179;
t277 = qJD(6) * t176;
t127 = t132 * qJD(4);
t252 = qJDD(2) * t331;
t268 = t170 * qJDD(2);
t228 = -t173 * t252 + t177 * t268;
t75 = qJD(2) * t127 + t228;
t30 = qJD(4) * t277 - t179 * qJDD(4) - t122 * t276 - t176 * t75;
t356 = t30 - t247;
t172 = sin(pkin(6));
t180 = cos(qJ(2));
t294 = t172 * t180;
t197 = t132 * t294;
t327 = pkin(8) + qJ(3);
t142 = t327 * t170;
t143 = t327 * t173;
t86 = -t177 * t142 + t331 * t143;
t324 = -qJD(1) * t197 + qJD(3) * t132 + qJD(4) * t86;
t257 = qJD(4) * t331;
t278 = qJD(4) * t177;
t259 = t170 * t278;
t126 = -t173 * t257 + t259;
t178 = sin(qJ(2));
t282 = qJD(1) * t172;
t258 = t178 * t282;
t360 = -qJ(5) * t126 + qJD(5) * t132 + t258;
t267 = t173 * qJDD(2);
t264 = qJD(4) * t239 + t170 * t252 + t177 * t267;
t74 = qJD(2) * t259 - t264;
t322 = qJ(5) * t74;
t272 = qJDD(1) * t172;
t255 = t180 * t272;
t279 = qJD(2) * t178;
t256 = qJD(1) * t279;
t345 = t172 * t256 + qJDD(3);
t206 = -t255 + t345;
t162 = pkin(3) * t173 + pkin(2);
t230 = t162 * qJDD(2);
t83 = -t230 + t206;
t185 = -qJD(5) * t346 + t322 + t83;
t335 = pkin(4) + pkin(9);
t14 = t335 * t75 + t185;
t174 = cos(pkin(6));
t271 = qJDD(1) * t174;
t151 = t173 * t271;
t262 = t180 * t282;
t270 = qJDD(2) * qJ(3);
t92 = t178 * t272 + t270 + (qJD(3) + t262) * qJD(2);
t52 = t151 + (-pkin(8) * qJDD(2) - t92) * t170;
t64 = t170 * t271 + t173 * t92;
t53 = pkin(8) * t267 + t64;
t139 = qJD(2) * qJ(3) + t258;
t281 = qJD(1) * t174;
t153 = t173 * t281;
t323 = pkin(8) * qJD(2);
t81 = t153 + (-t139 - t323) * t170;
t94 = t173 * t139 + t170 * t281;
t82 = t173 * t323 + t94;
t13 = -t177 * t53 - t82 * t257 - t81 * t278 + t331 * t52;
t232 = qJDD(5) - t13;
t6 = -pkin(5) * t74 - t335 * qJDD(4) + t232;
t35 = t177 * t82 - t331 * t81;
t287 = -qJD(5) - t35;
t286 = pkin(5) * t346 - t287;
t22 = -t335 * qJD(4) + t286;
t229 = qJD(3) - t262;
t111 = -qJD(2) * t162 + t229;
t187 = -qJ(5) * t346 + t111;
t27 = t335 * t122 + t187;
t9 = t176 * t22 + t179 * t27;
t2 = -qJD(6) * t9 - t176 * t14 + t179 * t6;
t344 = t9 * t354 + t2;
t226 = t176 * t27 - t179 * t22;
t1 = -t226 * qJD(6) + t179 * t14 + t176 * t6;
t233 = t226 * t354 + t1;
t359 = -t126 * pkin(5) + t324;
t351 = t179 * t354;
t91 = qJD(4) * t179 + t122 * t176;
t358 = t91 * t351;
t357 = -t335 * t127 + t360;
t244 = t176 * t354;
t66 = -qJDD(6) + t74;
t60 = t179 * t66;
t211 = -t244 * t354 - t60;
t169 = pkin(11) + qJ(4);
t164 = cos(t169);
t311 = cos(pkin(10));
t249 = t311 * t180;
t171 = sin(pkin(10));
t297 = t171 * t178;
t118 = -t174 * t249 + t297;
t250 = t311 * t178;
t296 = t171 * t180;
t120 = t174 * t296 + t250;
t238 = g(1) * t120 + g(2) * t118;
t349 = g(3) * t294 - t238;
t196 = t349 * t164;
t85 = t331 * t142 + t143 * t177;
t355 = -qJD(4) * t324 - qJDD(4) * t85 - t196;
t353 = 0.2e1 * t346 * qJD(4) + t228;
t182 = qJD(2) ^ 2;
t203 = (qJDD(2) * t180 - t178 * t182) * t172;
t163 = sin(t169);
t309 = qJ(5) * t163;
t350 = (pkin(4) * t164 + t309) * t294;
t291 = t176 * t180;
t148 = t172 * t291;
t295 = t172 * t178;
t117 = t170 * t174 + t173 * t295;
t293 = t173 * t174;
t207 = t170 * t295 - t293;
t198 = t331 * t207;
t58 = t117 * t177 + t198;
t41 = t179 * t58 + t148;
t116 = t346 ^ 2;
t336 = t122 ^ 2;
t348 = -t336 - t116;
t347 = -t336 + t116;
t36 = t177 * t81 + t331 * t82;
t29 = -qJD(4) * qJ(5) - t36;
t329 = t122 * pkin(5);
t23 = -t29 - t329;
t343 = t23 * t354 + t335 * t66;
t307 = qJDD(2) * pkin(2);
t101 = t206 - t307;
t214 = -t101 + t238;
t342 = t172 * (-g(3) * t180 + t256) + t214 + t307;
t223 = t263 * t294;
t32 = qJD(4) * t198 - qJD(2) * t223 + (qJD(4) * t117 + t280 * t294) * t177;
t59 = t331 * t117 - t177 * t207;
t341 = t172 * (-t180 * t74 - t279 * t346) - qJD(4) * t32 + qJDD(4) * t59;
t33 = qJD(2) * t197 + qJD(4) * t59;
t340 = t172 * (t122 * t279 - t180 * t75) - qJD(4) * t33 - qJDD(4) * t58;
t227 = t170 * (-t139 * t170 + t153) - t173 * t94;
t339 = t180 * t227 - (-qJD(2) * pkin(2) + t229) * t178;
t299 = t170 * t177;
t325 = qJD(1) * t223 + (qJD(3) * t170 + qJD(4) * t143) * t177 - t262 * t299 - qJD(3) * t263 + t142 * t257;
t337 = qJD(4) * t325 - qJDD(4) * t86 + t163 * t349;
t334 = pkin(4) * t75;
t131 = -t263 + t299;
t215 = -qJ(5) * t132 - t162;
t45 = t335 * t131 + t215;
t54 = pkin(5) * t132 + t85;
t20 = -t176 * t45 + t179 * t54;
t333 = qJD(6) * t20 + t359 * t176 - t357 * t179;
t21 = t176 * t54 + t179 * t45;
t332 = -qJD(6) * t21 + t357 * t176 + t359 * t179;
t330 = g(3) * t172;
t328 = t91 * t89;
t326 = -pkin(5) * t127 - t325;
t321 = t122 * t89;
t320 = t122 * t91;
t63 = -t170 * t92 + t151;
t319 = t170 * t63;
t240 = qJDD(4) * t176 - t179 * t75;
t31 = qJD(6) * t91 + t240;
t318 = t176 * t31;
t317 = t176 * t66;
t316 = t179 * t30;
t26 = t179 * t31;
t119 = t174 * t250 + t296;
t313 = -t118 * t162 + t119 * t327;
t121 = -t174 * t297 + t249;
t312 = -t120 * t162 + t121 * t327;
t310 = qJ(5) * t122;
t308 = qJD(4) * t36;
t306 = qJDD(4) * pkin(4);
t305 = t118 * t164;
t304 = t120 * t164;
t303 = t122 * t346;
t301 = t163 * t176;
t300 = t163 * t179;
t298 = t171 * t172;
t292 = t176 * t127;
t290 = t179 * t127;
t289 = t179 * t180;
t288 = t180 * t182;
t285 = qJDD(1) - g(3);
t166 = t170 ^ 2;
t168 = t173 ^ 2;
t283 = t166 + t168;
t269 = qJDD(4) * qJ(5);
t266 = g(3) * t295;
t265 = t172 * t289;
t260 = t172 * t279;
t251 = t172 * t311;
t76 = t119 * t163 + t164 * t251;
t77 = t119 * t164 - t163 * t251;
t254 = -t76 * pkin(4) + qJ(5) * t77;
t78 = t121 * t163 - t164 * t298;
t79 = t121 * t164 + t163 * t298;
t253 = -t78 * pkin(4) + qJ(5) * t79;
t248 = -t177 * t52 - t81 * t257 + t82 * t278 - t331 * t53;
t246 = -pkin(4) * t305 - t118 * t309 + t313;
t105 = t163 * t295 - t174 * t164;
t106 = t163 * t174 + t164 * t295;
t245 = -t105 * pkin(4) + qJ(5) * t106;
t241 = -pkin(4) * t304 - t120 * t309 + t312;
t237 = g(1) * t121 + g(2) * t119;
t236 = -pkin(4) * t127 + t360;
t137 = t162 * t294;
t235 = t295 * t327 + t137;
t234 = t176 * t226 + t179 * t9;
t225 = t122 * t127 + t131 * t75;
t224 = -t126 * t346 - t132 * t74;
t217 = qJD(4) * t126 - qJDD(4) * t132;
t216 = qJD(4) * t127 + qJDD(4) * t131;
t212 = -t176 * t58 + t265;
t205 = g(1) * t78 + g(2) * t76 + g(3) * t105;
t204 = -g(1) * t79 - g(2) * t77 - g(3) * t106;
t202 = -t351 * t354 + t317;
t10 = -qJD(4) * qJD(5) + t248 - t269;
t195 = -t349 + t255;
t194 = t173 * t64 - t237 - t319;
t193 = t122 * t32 + t33 * t346 - t58 * t74 - t59 * t75;
t192 = t122 * t126 - t127 * t346 + t131 * t74 - t132 * t75;
t191 = t205 + t13;
t190 = t204 - t248;
t7 = -pkin(5) * t75 - t10;
t189 = qJD(6) * t335 * t354 + t204 + t7;
t188 = -t195 + t345;
t39 = pkin(4) * t122 + t187;
t186 = t346 * t39 + qJDD(5) - t191;
t184 = -t230 + t188;
t183 = t122 * t325 + t324 * t346 - t74 * t85 - t75 * t86 - t237 - t266;
t112 = qJD(4) * t122;
t67 = pkin(4) * t131 + t215;
t65 = pkin(4) * t346 + t310;
t55 = -t131 * pkin(5) + t86;
t43 = -t112 + t74;
t40 = t335 * t346 + t310;
t28 = -qJD(4) * pkin(4) - t287;
t25 = t36 - t329;
t19 = t185 + t334;
t18 = t41 * qJD(6) + t176 * t33 + t179 * t260;
t17 = qJD(6) * t212 - t176 * t260 + t179 * t33;
t16 = t176 * t25 + t179 * t40;
t15 = -t176 * t40 + t179 * t25;
t11 = t232 - t306;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t285, 0, 0, 0, 0, 0, 0, t203 (-qJDD(2) * t178 - t288) * t172, 0, -g(3) + (t174 ^ 2 + (t178 ^ 2 + t180 ^ 2) * t172 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t173 * t203, -t170 * t203, t283 * t172 * t288 + (t117 * t173 + t170 * t207) * qJDD(2), t63 * t293 + t64 * t117 - g(3) + (-t339 * qJD(2) - t101 * t180 - t178 * t319) * t172, 0, 0, 0, 0, 0, 0, t340, -t341, t193, -t248 * t59 - t13 * t58 - t32 * t36 + t33 * t35 - g(3) + (t111 * t279 - t180 * t83) * t172, 0, 0, 0, 0, 0, 0, t193, -t340, t341, -t10 * t59 + t11 * t58 + t28 * t33 + t29 * t32 - g(3) + (-t180 * t19 + t279 * t39) * t172, 0, 0, 0, 0, 0, 0, t17 * t354 + t31 * t59 - t32 * t89 - t41 * t66, -t18 * t354 - t212 * t66 - t30 * t59 - t32 * t91, -t17 * t91 - t18 * t89 + t212 * t31 + t30 * t41, -t1 * t212 - t17 * t226 + t18 * t9 + t2 * t41 - t23 * t32 + t59 * t7 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t195, -t285 * t295 + t237, 0, 0, t166 * qJDD(2), 0.2e1 * t170 * t267, 0, t168 * qJDD(2), 0, 0, t342 * t173, -t342 * t170, -t266 + t194 + (t229 * qJD(2) + t270) * t283, -t227 * qJD(3) + t214 * pkin(2) + t194 * qJ(3) + (-g(3) * (pkin(2) * t180 + qJ(3) * t178) + t339 * qJD(1)) * t172, t224, t192, -t217, t225, -t216, 0, t111 * t127 - t122 * t258 + t131 * t83 - t162 * t75 + t355, -t111 * t126 + t132 * t83 + t162 * t74 - t258 * t346 + t337, -t126 * t35 - t127 * t36 - t13 * t132 + t131 * t248 + t183, -g(1) * t312 - g(2) * t313 - g(3) * t235 - t111 * t258 - t13 * t85 - t83 * t162 - t248 * t86 + t324 * t35 - t325 * t36, 0, t217, t216, t224, t192, t225, t10 * t131 + t11 * t132 - t126 * t28 + t127 * t29 + t183, t122 * t236 - t127 * t39 - t131 * t19 - t67 * t75 - t355, t126 * t39 - t132 * t19 + t236 * t346 + t67 * t74 - t337, t19 * t67 - t10 * t86 + t11 * t85 - g(1) * t241 - g(2) * t246 - g(3) * (t235 + t350) - t236 * t39 + t325 * t29 + t324 * t28, t91 * t292 + (-t176 * t30 + t276 * t91) * t131 (-t176 * t89 + t179 * t91) * t127 + (-t318 - t316 + (-t176 * t91 - t179 * t89) * qJD(6)) * t131, -t131 * t317 - t126 * t91 - t132 * t30 + (t131 * t276 + t292) * t354, -t89 * t290 + (t277 * t89 - t26) * t131, -t131 * t60 + t126 * t89 - t132 * t31 + (-t131 * t277 + t290) * t354, -t126 * t354 - t132 * t66, -t20 * t66 + t2 * t132 + t226 * t126 + t55 * t31 - t23 * t290 - g(1) * (-t120 * t301 + t121 * t179) - g(2) * (-t118 * t301 + t119 * t179) + t326 * t89 - (t163 * t291 + t178 * t179) * t330 + (-t7 * t179 + t23 * t277) * t131 + t332 * t354, t21 * t66 - t1 * t132 + t9 * t126 - t55 * t30 + t23 * t292 - g(1) * (-t120 * t300 - t121 * t176) - g(2) * (-t118 * t300 - t119 * t176) + t326 * t91 - (t163 * t289 - t176 * t178) * t330 + (t7 * t176 + t23 * t276) * t131 - t333 * t354, t20 * t30 - t21 * t31 - t332 * t91 - t333 * t89 + t234 * t127 - t196 + (t1 * t179 - t176 * t2 + (-t176 * t9 + t179 * t226) * qJD(6)) * t131, t1 * t21 + t2 * t20 + t7 * t55 - g(1) * (pkin(5) * t121 - pkin(9) * t304 + t241) - g(2) * (pkin(5) * t119 - pkin(9) * t305 + t246) - g(3) * (t137 + t350) + t333 * t9 - t332 * t226 + t326 * t23 - (pkin(9) * t164 * t180 + (pkin(5) + t327) * t178) * t330; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t267, t268, -t283 * t182, qJD(2) * t227 + t188 - t307, 0, 0, 0, 0, 0, 0, t353 (-t122 - t261) * qJD(4) + t264, t348, t122 * t36 - t346 * t35 + t184, 0, 0, 0, 0, 0, 0, t348, -t353, t112 + t74, t334 + t322 - t122 * t29 + (-qJD(5) - t28) * t346 + t184, 0, 0, 0, 0, 0, 0, t202 + t321, t320 - t211, -t356 * t176 - t26 + t358, t122 * t23 - t344 * t176 + t233 * t179 + t349; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t303, t347 (t122 - t261) * qJD(4) + t264, -t303, -t228, qJDD(4), -t111 * t346 + t191 + t308, -qJD(4) * t35 + t111 * t122 - t190, 0, 0, qJDD(4), t43, t228, t303, t347, -t303, pkin(4) * t74 - qJ(5) * t75 + (-t29 - t36) * t346 + (t28 + t287) * t122, t122 * t65 + t186 - 0.2e1 * t306 - t308, 0.2e1 * t269 - t122 * t39 + t346 * t65 + (0.2e1 * qJD(5) + t35) * qJD(4) + t190, -t11 * pkin(4) - g(1) * t253 - g(2) * t254 - g(3) * t245 - t10 * qJ(5) - t28 * t36 + t287 * t29 - t39 * t65, -t244 * t91 - t316, -t26 - t358 + (t30 + t247) * t176, t211 + t320, t351 * t89 + t318, t202 - t321, t354 * t122, qJ(5) * t31 - t122 * t226 - t15 * t354 + t189 * t176 + t343 * t179 + t286 * t89, -qJ(5) * t30 - t122 * t9 + t16 * t354 - t343 * t176 + t189 * t179 + t286 * t91, t15 * t91 + t16 * t89 + (-t346 * t9 - t335 * t30 - t2 + (t335 * t89 - t9) * qJD(6)) * t179 + (-t346 * t226 + t335 * t31 - t1 + (-t335 * t91 - t226) * qJD(6)) * t176 + t205, t7 * qJ(5) - t9 * t16 + t226 * t15 - g(1) * (-pkin(9) * t78 + t253) - g(2) * (-pkin(9) * t76 + t254) - g(3) * (-pkin(9) * t105 + t245) + t286 * t23 - (qJD(6) * t234 + t1 * t176 + t2 * t179) * t335; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, qJDD(4) - t303, -qJD(4) ^ 2 - t116, qJD(4) * t29 + t186 - t306, 0, 0, 0, 0, 0, 0, -qJD(4) * t89 + t211, -qJD(4) * t91 + t202, t356 * t179 + (t354 * t91 - t31) * t176, -t23 * qJD(4) + t233 * t176 + t344 * t179 - t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t328, -t89 ^ 2 + t91 ^ 2, -t356, -t328, -t240 + (-qJD(6) + t354) * t91, -t66, -t23 * t91 - g(1) * (-t120 * t176 + t179 * t78) - g(2) * (-t118 * t176 + t179 * t76) - g(3) * (t105 * t179 + t148) + t344, t23 * t89 - g(1) * (-t120 * t179 - t176 * t78) - g(2) * (-t118 * t179 - t176 * t76) - g(3) * (-t105 * t176 + t265) - t233, 0, 0;];
tau_reg  = t3;
