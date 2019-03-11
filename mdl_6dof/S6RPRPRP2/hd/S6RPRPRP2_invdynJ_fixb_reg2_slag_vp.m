% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPRPRP2
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
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPRP2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP2_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:06:20
% EndTime: 2019-03-09 03:06:28
% DurationCPUTime: 5.55s
% Computational Cost: add. (7409->545), mult. (15915->636), div. (0->0), fcn. (11118->14), ass. (0->274)
t195 = cos(qJ(3));
t321 = cos(pkin(10));
t263 = t321 * t195;
t158 = qJD(1) * t263;
t187 = sin(pkin(10));
t192 = sin(qJ(3));
t290 = qJD(1) * t192;
t127 = t187 * t290 - t158;
t123 = qJD(5) + t127;
t138 = t187 * t195 + t192 * t321;
t191 = sin(qJ(5));
t287 = qJD(5) * t191;
t269 = t138 * t287;
t224 = -t187 * t192 + t263;
t132 = t224 * qJD(3);
t194 = cos(qJ(5));
t306 = t132 * t194;
t228 = t269 - t306;
t304 = t138 * t194;
t129 = t138 * qJD(3);
t281 = t192 * qJDD(1);
t242 = -qJDD(1) * t263 + t187 * t281;
t95 = qJD(1) * t129 + t242;
t94 = qJDD(5) + t95;
t226 = -t123 * t228 + t304 * t94;
t130 = t138 * qJD(1);
t109 = qJD(3) * t191 + t130 * t194;
t282 = qJD(3) * qJD(5);
t285 = qJD(1) * qJD(3);
t266 = t192 * t285;
t214 = qJDD(1) * t138 - t187 * t266;
t353 = qJD(3) * t158 + t214;
t47 = -t191 * qJDD(3) + t130 * t287 + (-t353 - t282) * t194;
t333 = t109 * t129 + t224 * t47;
t361 = t226 - t333;
t188 = sin(pkin(9));
t165 = pkin(1) * t188 + pkin(7);
t149 = t165 * qJDD(1);
t255 = -qJD(2) * qJD(3) - t149;
t359 = t130 * t129 - t224 * t353;
t320 = pkin(1) * qJDD(1);
t183 = qJ(3) + pkin(10);
t171 = sin(t183);
t184 = qJ(1) + pkin(9);
t172 = sin(t184);
t174 = cos(t184);
t248 = g(1) * t174 + g(2) * t172;
t358 = t171 * t248;
t151 = t165 * qJD(1);
t256 = qJ(4) * qJD(1) + t151;
t357 = t195 * t256;
t173 = cos(t183);
t356 = -pkin(4) * t173 - pkin(8) * t171;
t342 = g(1) * t172;
t265 = g(2) * t174 - t342;
t189 = cos(pkin(9));
t167 = -t189 * pkin(1) - pkin(2);
t181 = t195 * pkin(3);
t355 = t167 - t181;
t286 = qJD(5) * t194;
t268 = t138 * t286;
t307 = t132 * t191;
t229 = t268 + t307;
t305 = t138 * t191;
t209 = -t123 * t229 - t305 * t94;
t107 = -qJD(3) * t194 + t130 * t191;
t201 = -qJDD(3) * t194 + t191 * t353;
t48 = qJD(5) * t109 + t201;
t335 = t129 * t107 - t224 * t48;
t354 = t209 - t335;
t352 = t130 * qJD(3);
t289 = qJD(2) * t192;
t113 = t289 + t357;
t102 = t187 * t113;
t180 = t195 * qJD(2);
t112 = -t192 * t256 + t180;
t105 = qJD(3) * pkin(3) + t112;
t58 = t105 * t321 - t102;
t52 = -qJD(3) * pkin(4) - t58;
t27 = pkin(5) * t107 - qJ(6) * t109 + t52;
t344 = pkin(3) * t187;
t164 = pkin(8) + t344;
t326 = t164 * t94;
t351 = t123 * t27 - t326;
t339 = g(3) * t171;
t350 = t173 * t248 + t339;
t84 = -pkin(4) * t224 - pkin(8) * t138 + t355;
t295 = qJ(4) + t165;
t134 = t295 * t195;
t260 = t295 * t192;
t91 = t134 * t321 - t187 * t260;
t332 = t191 * t84 + t194 * t91;
t259 = qJD(3) * t295;
t114 = qJD(4) * t195 - t192 * t259;
t220 = -qJD(4) * t192 - t195 * t259;
t69 = t114 * t321 + t187 * t220;
t288 = qJD(3) * t192;
t277 = pkin(3) * t288;
t80 = pkin(4) * t129 - pkin(8) * t132 + t277;
t10 = -qJD(5) * t332 - t191 * t69 + t194 * t80;
t349 = t109 ^ 2;
t348 = t123 ^ 2;
t347 = t130 ^ 2;
t346 = pkin(5) * t94;
t193 = sin(qJ(1));
t345 = pkin(1) * t193;
t343 = pkin(3) * t192;
t338 = g(3) * t173;
t337 = g(3) * t195;
t37 = t48 * t304;
t336 = -t107 * t306 - t37;
t178 = t195 * qJDD(2);
t284 = qJD(1) * qJD(4);
t56 = qJDD(3) * pkin(3) + t178 - qJD(3) * t357 + (-qJ(4) * qJDD(1) + t255 - t284) * t192;
t280 = t195 * qJDD(1);
t273 = -t192 * qJDD(2) + t195 * t255;
t87 = -t151 * t288 - t273;
t60 = t195 * t284 + (-t266 + t280) * qJ(4) + t87;
t21 = -t187 * t60 + t321 * t56;
t22 = t187 * t56 + t321 * t60;
t62 = t112 * t321 - t102;
t79 = pkin(3) * t290 + pkin(4) * t130 + pkin(8) * t127;
t30 = t191 * t79 + t194 * t62;
t334 = -t127 * t132 - t138 * t95;
t331 = qJ(6) * t94;
t264 = t321 * t113;
t59 = t105 * t187 + t264;
t53 = qJD(3) * pkin(8) + t59;
t126 = qJD(1) * t355 + qJD(4);
t70 = pkin(4) * t127 - pkin(8) * t130 + t126;
t26 = t191 * t70 + t194 * t53;
t15 = qJ(6) * t123 + t26;
t330 = t123 * t15;
t25 = -t191 * t53 + t194 * t70;
t329 = t123 * t25;
t328 = t123 * t26;
t327 = t164 * t47;
t325 = t191 * t47;
t85 = t191 * t94;
t324 = t194 * t48;
t86 = t194 * t94;
t323 = -t107 * t286 - t191 * t48;
t243 = pkin(5) * t191 - qJ(6) * t194;
t61 = t112 * t187 + t264;
t322 = -qJD(6) * t191 + t123 * t243 - t61;
t319 = t107 * t123;
t318 = t107 * t130;
t317 = t107 * t164;
t316 = t107 * t194;
t315 = t109 * t107;
t314 = t109 * t130;
t313 = t109 * t164;
t312 = t109 * t191;
t311 = t123 * t130;
t257 = t123 * t191;
t310 = t123 * t194;
t309 = t130 * t127;
t303 = t151 * t192;
t302 = t151 * t195;
t301 = t171 * t174;
t300 = t172 * t191;
t299 = t172 * t194;
t298 = t173 * t174;
t297 = t174 * t191;
t296 = t174 * t194;
t294 = qJD(6) - t25;
t293 = (g(1) * t296 + g(2) * t299) * t171;
t170 = t181 + pkin(2);
t196 = cos(qJ(1));
t182 = t196 * pkin(1);
t292 = t170 * t174 + t182;
t185 = t192 ^ 2;
t186 = t195 ^ 2;
t291 = t185 - t186;
t152 = qJD(1) * t167;
t150 = qJDD(1) * t167;
t279 = pkin(3) * t266 + qJDD(4);
t76 = t109 * t307;
t278 = t109 * t268 - t305 * t47 + t76;
t275 = t189 * t320;
t198 = qJD(1) ^ 2;
t274 = t192 * t198 * t195;
t272 = t181 - t356;
t271 = t321 * pkin(3);
t270 = qJD(5) * t123 * t164;
t267 = t107 ^ 2 - t349;
t20 = qJDD(3) * pkin(8) + t22;
t36 = -qJDD(1) * pkin(2) - pkin(3) * t280 + t95 * pkin(4) - pkin(8) * t353 - t275 + t279;
t262 = t191 * t20 - t194 * t36 + t286 * t53 + t287 * t70;
t100 = t109 * t287;
t261 = -t194 * t47 - t100;
t68 = t114 * t187 - t220 * t321;
t90 = t134 * t187 + t260 * t321;
t254 = pkin(4) * t298 + pkin(8) * t301 + t292;
t253 = t195 * t266;
t252 = -g(2) * t301 + t171 * t342;
t19 = -qJDD(3) * pkin(4) - t21;
t251 = -pkin(4) * t171 - t343;
t116 = t173 * t300 + t296;
t118 = t173 * t297 - t299;
t250 = g(1) * t116 - g(2) * t118;
t117 = t173 * t299 - t297;
t119 = t173 * t296 + t300;
t249 = g(1) * t117 - g(2) * t119;
t246 = g(1) * t193 - g(2) * t196;
t190 = -qJ(4) - pkin(7);
t245 = -t174 * t190 - t345;
t244 = pkin(5) * t194 + qJ(6) * t191;
t14 = -pkin(5) * t123 + t294;
t241 = t14 * t194 - t15 * t191;
t240 = t14 * t191 + t15 * t194;
t239 = -t191 * t26 - t194 * t25;
t238 = t191 * t25 - t194 * t26;
t29 = -t191 * t62 + t194 * t79;
t40 = -t191 * t91 + t194 * t84;
t121 = t289 + t302;
t234 = -t123 * t287 - t127 * t257 + t86;
t233 = t123 * t286 + t127 * t310 + t85;
t232 = pkin(4) + t244;
t230 = t270 + t338;
t3 = t191 * t36 + t194 * t20 + t286 * t70 - t287 * t53;
t9 = t191 * t80 + t194 * t69 + t286 * t84 - t287 * t91;
t227 = t107 * t269 + t336;
t225 = t123 * t52 - t326;
t223 = t19 + t230;
t222 = -qJD(1) * t152 + t248;
t221 = t107 * t257 - t324;
t219 = -t234 - t318;
t218 = 0.2e1 * qJD(3) * t152 - qJDD(3) * t165;
t217 = -t338 + t358;
t216 = g(1) * t118 + g(2) * t116 + t191 * t339 - t262;
t215 = -t164 * t324 - t350;
t111 = qJDD(1) * t355 + t279;
t197 = qJD(3) ^ 2;
t213 = -t165 * t197 - 0.2e1 * t150 - t265;
t212 = (-g(1) * (-t170 + t356) + g(2) * t190) * t172;
t1 = qJD(6) * t123 + t3 + t331;
t2 = qJDD(6) + t262 - t346;
t211 = qJD(5) * t241 + t1 * t194 + t2 * t191;
t210 = qJD(5) * t239 + t191 * t262 + t3 * t194;
t208 = t107 * t229 + t305 * t48;
t120 = t180 - t303;
t88 = -qJD(3) * t121 - t149 * t192 + t178;
t207 = -t88 * t192 + t87 * t195 + (-t120 * t195 - t121 * t192) * qJD(3);
t206 = t109 * t27 + qJDD(6) - t216;
t205 = t209 + t335;
t204 = -g(1) * t119 - g(2) * t117 - t194 * t339 + t3;
t200 = t109 * t123 - t130 * t286 - t191 * t282 - t201;
t166 = -t271 - pkin(4);
t148 = qJDD(3) * t195 - t192 * t197;
t147 = qJDD(3) * t192 + t195 * t197;
t142 = pkin(8) * t298;
t140 = t172 * t173 * pkin(8);
t133 = -t271 - t232;
t125 = t127 ^ 2;
t97 = qJD(3) * t132 + qJDD(3) * t138;
t96 = -qJD(3) * t129 + qJDD(3) * t224;
t66 = pkin(5) * t109 + qJ(6) * t107;
t46 = t138 * t243 + t90;
t39 = t123 * t129 - t224 * t94;
t34 = pkin(5) * t224 - t40;
t31 = -qJ(6) * t224 + t332;
t28 = -t47 + t319;
t24 = -pkin(5) * t130 - t29;
t23 = qJ(6) * t130 + t30;
t16 = t233 - t314;
t13 = t109 * t310 - t325;
t12 = t243 * t132 + (qJD(5) * t244 - qJD(6) * t194) * t138 + t68;
t11 = -t109 * t228 - t304 * t47;
t8 = -pkin(5) * t129 - t10;
t7 = qJ(6) * t129 - qJD(6) * t224 + t9;
t6 = t226 + t333;
t5 = pkin(5) * t48 + qJ(6) * t47 - qJD(6) * t109 + t19;
t4 = [0, 0, 0, 0, 0, qJDD(1), t246, g(1) * t196 + g(2) * t193, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t265 + 0.2e1 * t275, -0.2e1 * t188 * t320 + t248, 0 (t246 + (t188 ^ 2 + t189 ^ 2) * t320) * pkin(1), qJDD(1) * t185 + 0.2e1 * t253, 0.2e1 * t192 * t280 - 0.2e1 * t285 * t291, t147, qJDD(1) * t186 - 0.2e1 * t253, t148, 0, t192 * t218 + t195 * t213, -t192 * t213 + t195 * t218 (t185 + t186) * t149 + t207 - t248, t150 * t167 - g(1) * (-pkin(2) * t172 + pkin(7) * t174 - t345) - g(2) * (pkin(2) * t174 + pkin(7) * t172 + t182) + t207 * t165, t130 * t132 + t138 * t353, t334 - t359, t97, t127 * t129 - t224 * t95, t96, 0, -qJDD(3) * t90 - t111 * t224 + t126 * t129 + t355 * t95 - t265 * t173 + (t127 * t343 - t68) * qJD(3), -t69 * qJD(3) - t91 * qJDD(3) + t111 * t138 + t126 * t132 + t130 * t277 + t353 * t355 - t252, -t69 * t127 - t59 * t129 + t68 * t130 - t58 * t132 - t21 * t138 + t22 * t224 + t353 * t90 - t91 * t95 - t248, t22 * t91 + t59 * t69 - t21 * t90 - t58 * t68 + t111 * t355 + t126 * t277 - g(1) * (-t170 * t172 + t245) - g(2) * (-t172 * t190 + t292) t11, t227 - t278, t6, t208, t354, t39, t52 * t307 + t10 * t123 + t107 * t68 + t129 * t25 + t224 * t262 + t40 * t94 + t48 * t90 + (t19 * t191 + t286 * t52) * t138 + t249, t52 * t306 + t109 * t68 - t123 * t9 - t129 * t26 + t224 * t3 - t332 * t94 - t47 * t90 + (t19 * t194 - t287 * t52) * t138 - t250, -t10 * t109 - t107 * t9 + t40 * t47 - t332 * t48 + t239 * t132 + (qJD(5) * t238 - t191 * t3 + t194 * t262) * t138 + t252, -g(1) * t245 - g(2) * t254 + t25 * t10 + t19 * t90 + t26 * t9 - t262 * t40 + t3 * t332 + t52 * t68 + t212, t11, t6, -t107 * t228 + t278 + t37, t39, -t354, t208, t27 * t307 + t107 * t12 - t123 * t8 - t129 * t14 + t224 * t2 - t34 * t94 + t46 * t48 + (t191 * t5 + t27 * t286) * t138 + t249, -t107 * t7 + t109 * t8 - t31 * t48 - t34 * t47 + t241 * t132 + (-qJD(5) * t240 - t1 * t191 + t194 * t2) * t138 + t252, -t27 * t306 - t1 * t224 - t109 * t12 + t123 * t7 + t129 * t15 + t31 * t94 + t46 * t47 + (-t194 * t5 + t27 * t287) * t138 + t250, t1 * t31 + t15 * t7 + t5 * t46 + t27 * t12 + t2 * t34 + t14 * t8 - g(1) * (-pkin(5) * t117 - qJ(6) * t116 + t245) - g(2) * (pkin(5) * t119 + qJ(6) * t118 + t254) + t212; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, t148, -t147, 0, t192 * t87 + t195 * t88 - g(3) + (-t120 * t192 + t121 * t195) * qJD(3), 0, 0, 0, 0, 0, 0, t96, -t97, t334 + t359, -t129 * t58 + t132 * t59 + t138 * t22 + t21 * t224 - g(3), 0, 0, 0, 0, 0, 0, t205, -t361, t76 + (-t325 + (t107 * t191 + t109 * t194) * qJD(5)) * t138 + t336, t129 * t52 - t132 * t238 + t138 * t210 - t19 * t224 - g(3), 0, 0, 0, 0, 0, 0, t205, t227 + t278, t361, t129 * t27 + t132 * t240 + t138 * t211 - t224 * t5 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t274, t291 * t198, t281, t274, t280, qJDD(3), -t337 + t178 + (t121 - t302) * qJD(3) + (t222 + t255) * t192, g(3) * t192 + (t120 + t303) * qJD(3) + t222 * t195 + t273, 0, 0, t309, -t125 + t347 (t158 + t127) * qJD(3) + t214, -t309, -t242, qJDD(3), t61 * qJD(3) - t126 * t130 + (qJDD(3) * t321 - t127 * t290) * pkin(3) + t217 + t21, qJD(3) * t62 + t126 * t127 + (-qJDD(3) * t187 - t130 * t290) * pkin(3) - t22 + t350, -t353 * t271 - t95 * t344 - (-t59 + t61) * t130 + (t62 - t58) * t127, t58 * t61 - t59 * t62 + (t321 * t21 - t337 + t187 * t22 + (-qJD(1) * t126 + t248) * t192) * pkin(3), t13 (-t312 - t316) * t127 + t261 + t323, t16, t221, -t219, -t311, -t107 * t61 - t123 * t29 - t130 * t25 + t166 * t48 + t191 * t225 - t194 * t223 + t293, -t109 * t61 + t123 * t30 + t130 * t26 - t166 * t47 + t225 * t194 + (t223 - t358) * t191, t107 * t30 + t109 * t29 + (-t127 * t25 + t3 + (-t25 + t313) * qJD(5)) * t194 + (-t127 * t26 - t327 + t262 + (-t26 + t317) * qJD(5)) * t191 + t215, t19 * t166 - t26 * t30 - t25 * t29 - t52 * t61 - g(1) * (t174 * t251 + t142) - g(2) * (t172 * t251 + t140) - g(3) * t272 + t210 * t164, t13, t16, t100 + (t109 * t127 + t48) * t191 + (t47 + t319) * t194, -t311, t219, t221, t123 * t24 + t130 * t14 + t133 * t48 + t322 * t107 + (-t230 - t5) * t194 + t351 * t191 + t293, t107 * t23 - t109 * t24 + (t127 * t14 + t1 + (t14 + t313) * qJD(5)) * t194 + (-t127 * t15 - t327 + t2 + (-t15 + t317) * qJD(5)) * t191 + t215, -t123 * t23 - t130 * t15 + t133 * t47 - t322 * t109 - t351 * t194 + (t217 - t5 - t270) * t191, t5 * t133 - t15 * t23 - t14 * t24 - g(1) * (-t174 * t343 + t142) - g(2) * (-t172 * t343 + t140) - g(3) * (t173 * t244 + t272) + t322 * t27 + t211 * t164 + t232 * t358; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t242 + 0.2e1 * t352 (t158 - t127) * qJD(3) + t214, -t125 - t347, t127 * t59 + t130 * t58 + t111 + t265, 0, 0, 0, 0, 0, 0, t234 - t318, -t194 * t348 - t314 - t85 (-t107 * t127 + t47) * t194 + t109 * t257 + t323, -t130 * t52 + (-t262 + t328) * t194 + (t3 - t329) * t191 + t265, 0, 0, 0, 0, 0, 0, -t191 * t348 - t318 + t86 (t312 - t316) * t127 - t261 + t323, t233 + t314, -t130 * t27 + (-t2 + t330) * t194 + (t123 * t14 + t1) * t191 + t265; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t315, -t267, t28, -t315, t200, t94, -t109 * t52 + t216 + t328, t107 * t52 - t204 + t329, 0, 0, t315, t28, t267, t94, -t200, -t315, -t107 * t66 - t206 + t328 + 0.2e1 * t346, pkin(5) * t47 - qJ(6) * t48 + (t15 - t26) * t109 + (t14 - t294) * t107, 0.2e1 * t331 - t107 * t27 + t109 * t66 + (0.2e1 * qJD(6) - t25) * t123 + t204, t1 * qJ(6) - t2 * pkin(5) - t27 * t66 - t14 * t26 - g(1) * (-pkin(5) * t118 + qJ(6) * t119) - g(2) * (-pkin(5) * t116 + qJ(6) * t117) + t243 * t339 + t294 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(5) - t242 + t315 - t352, t28, -t348 - t349, t206 - t330 - t346;];
tau_reg  = t4;
