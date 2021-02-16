% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRPRP1
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
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:40
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRPRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:38:48
% EndTime: 2021-01-16 02:39:08
% DurationCPUTime: 5.96s
% Computational Cost: add. (5630->501), mult. (13250->659), div. (0->0), fcn. (10499->14), ass. (0->262)
t193 = sin(qJ(3));
t194 = sin(qJ(2));
t188 = sin(pkin(6));
t288 = qJD(1) * t188;
t267 = t194 * t288;
t332 = qJD(3) * pkin(3);
t358 = -t193 * t332 + t267;
t186 = sin(pkin(11));
t196 = cos(qJ(3));
t321 = cos(pkin(11));
t251 = t321 * t193;
t149 = t186 * t196 + t251;
t140 = t149 * qJD(3);
t250 = t321 * t196;
t305 = t186 * t193;
t223 = t250 - t305;
t143 = t223 * qJD(3);
t369 = pkin(4) * t140 - pkin(9) * t143 - t358;
t197 = cos(qJ(2));
t266 = t197 * t288;
t109 = t223 * t266;
t191 = qJ(4) + pkin(8);
t256 = qJD(3) * t191;
t133 = qJD(4) * t196 - t193 * t256;
t216 = -qJD(4) * t193 - t196 * t256;
t71 = t321 * t133 + t186 * t216;
t362 = t71 - t109;
t141 = t149 * qJD(2);
t192 = sin(qJ(5));
t195 = cos(qJ(5));
t115 = qJD(3) * t192 + t141 * t195;
t314 = t115 * t192;
t189 = cos(pkin(6));
t301 = t188 * t194;
t144 = t189 * t196 - t193 * t301;
t183 = qJ(3) + pkin(11);
t178 = sin(t183);
t179 = cos(t183);
t127 = t178 * t301 - t189 * t179;
t322 = cos(pkin(10));
t254 = t322 * t194;
t187 = sin(pkin(10));
t302 = t187 * t197;
t136 = t189 * t254 + t302;
t255 = t188 * t322;
t92 = t136 * t178 + t179 * t255;
t253 = t322 * t197;
t303 = t187 * t194;
t134 = t189 * t303 - t253;
t304 = t187 * t188;
t94 = t134 * t178 + t179 * t304;
t220 = -g(1) * t94 + g(2) * t92 + g(3) * t127;
t368 = g(3) * t188;
t279 = t189 * qJDD(1);
t165 = t196 * t279;
t282 = qJD(1) * qJD(2);
t125 = qJDD(2) * pkin(8) + (qJDD(1) * t194 + t197 * t282) * t188;
t287 = qJD(1) * t189;
t207 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + qJD(3) * t287 + t125;
t246 = t191 * qJD(2) + t267;
t234 = t246 * qJD(3);
t36 = qJDD(3) * pkin(3) - t207 * t193 - t196 * t234 + t165;
t37 = (-t234 + t279) * t193 + t207 * t196;
t13 = t186 * t36 + t321 * t37;
t11 = qJDD(3) * pkin(9) + t13;
t106 = -t246 * t193 + t196 * t287;
t102 = t106 + t332;
t107 = t193 * t287 + t246 * t196;
t252 = t321 * t107;
t46 = t186 * t102 + t252;
t41 = qJD(3) * pkin(9) + t46;
t177 = t196 * pkin(3) + pkin(2);
t131 = -t177 * qJD(2) + qJD(4) - t266;
t166 = qJD(2) * t250;
t286 = qJD(2) * t193;
t138 = t186 * t286 - t166;
t57 = pkin(4) * t138 - pkin(9) * t141 + t131;
t23 = t192 * t57 + t195 * t41;
t299 = t188 * t197;
t258 = qJDD(1) * t299;
t261 = t194 * t282;
t239 = t188 * t261 - t258;
t320 = qJDD(2) * pkin(2);
t124 = t239 - t320;
t277 = t196 * qJDD(2);
t281 = qJD(2) * qJD(3);
t260 = t193 * t281;
t357 = pkin(3) * t260 + qJDD(4);
t210 = t149 * qJDD(2) - t186 * t260;
t359 = qJD(3) * t166 + t210;
t278 = t193 * qJDD(2);
t238 = -qJDD(2) * t250 + t186 * t278;
t90 = qJD(2) * t140 + t238;
t26 = -pkin(3) * t277 + t90 * pkin(4) - pkin(9) * t359 + t124 + t357;
t25 = t195 * t26;
t209 = -t23 * qJD(5) - t11 * t192 + t25;
t280 = qJD(3) * qJD(5);
t284 = qJD(5) * t192;
t218 = t192 * qJDD(3) - t141 * t284 + (t359 + t280) * t195;
t335 = qJ(6) * t218;
t86 = qJDD(5) + t90;
t353 = pkin(5) * t86;
t1 = -qJD(6) * t115 + t209 - t335 + t353;
t132 = qJD(5) + t138;
t113 = -t195 * qJD(3) + t141 * t192;
t15 = -qJ(6) * t113 + t23;
t331 = t132 * t15;
t367 = t1 + t331;
t12 = -t186 * t37 + t321 * t36;
t10 = -qJDD(3) * pkin(4) - t12;
t201 = -t195 * qJDD(3) + t192 * t359;
t43 = t115 * qJD(5) + t201;
t5 = t43 * pkin(5) + qJDD(6) + t10;
t366 = t220 - t5;
t365 = t109 * t192 + t369 * t195;
t344 = pkin(3) * t186;
t173 = pkin(9) + t344;
t292 = qJ(6) + t173;
t248 = qJD(5) * t292;
t311 = t138 * t192;
t98 = t186 * t107;
t50 = t321 * t106 - t98;
t274 = pkin(3) * t286;
t72 = pkin(4) * t141 + pkin(9) * t138 + t274;
t336 = t192 * t72 + t195 * t50;
t364 = -qJ(6) * t311 + qJD(6) * t195 - t192 * t248 - t336;
t310 = t138 * t195;
t63 = t195 * t72;
t325 = -pkin(5) * t141 - qJ(6) * t310 - t195 * t248 - t63 + (-qJD(6) + t50) * t192;
t283 = qJD(5) * t195;
t88 = -pkin(4) * t223 - pkin(9) * t149 - t177;
t363 = t369 * t192 + t362 * t195 + t88 * t283;
t323 = t133 * t186 - t149 * t266 - t321 * t216;
t361 = t132 * t314;
t135 = -t189 * t253 + t303;
t137 = t189 * t302 + t254;
t241 = g(1) * t137 + g(2) * t135;
t214 = -g(3) * t299 + t241;
t360 = t214 * t178;
t128 = t178 * t189 + t179 * t301;
t293 = t195 * t197;
t269 = t188 * t293;
t91 = t134 * t179 - t178 * t304;
t93 = t136 * t179 - t178 * t255;
t356 = -g(3) * (-t128 * t192 - t269) - g(2) * (t135 * t195 - t192 * t93) - g(1) * (t137 * t195 + t192 * t91);
t198 = qJD(3) ^ 2;
t355 = -pkin(8) * t198 + t188 * (-g(3) * t197 + t261) - t124 + t241 + t320;
t354 = t115 ^ 2;
t22 = -t192 * t41 + t195 * t57;
t14 = -qJ(6) * t115 + t22;
t8 = pkin(5) * t132 + t14;
t347 = t14 - t8;
t158 = t191 * t196;
t111 = t321 * t158 - t191 * t305;
t103 = t195 * t111;
t236 = -qJ(6) * t143 - qJD(6) * t149;
t346 = pkin(5) * t140 - t192 * t71 + t236 * t195 + (-t103 + (qJ(6) * t149 - t88) * t192) * qJD(5) + t365;
t263 = t149 * t283;
t345 = -qJ(6) * t263 + (-qJD(5) * t111 + t236) * t192 + t363;
t343 = pkin(3) * t193;
t339 = t113 * pkin(5);
t338 = t192 * t5;
t337 = t195 * pkin(5);
t334 = qJ(6) * t43;
t333 = qJD(2) * pkin(2);
t330 = t192 * t218;
t329 = t192 * t86;
t328 = t192 * t88 + t103;
t327 = -t113 * t283 - t192 * t43;
t297 = t192 * t143;
t226 = t263 + t297;
t324 = t226 * pkin(5) + t323;
t319 = t113 * t132;
t318 = t113 * t138;
t317 = t113 * t141;
t316 = t115 * t132;
t315 = t115 * t141;
t312 = t134 * t192;
t309 = t149 * t192;
t308 = t149 * t195;
t307 = t179 * t192;
t306 = t179 * t195;
t300 = t188 * t196;
t296 = t192 * t194;
t295 = t192 * t197;
t294 = t195 * t143;
t291 = qJDD(1) - g(3);
t290 = -t134 * t191 - t137 * t177;
t184 = t193 ^ 2;
t289 = -t196 ^ 2 + t184;
t285 = qJD(2) * t194;
t275 = t220 * t192;
t249 = qJD(6) + t339;
t45 = t321 * t102 - t98;
t40 = -qJD(3) * pkin(4) - t45;
t28 = t249 + t40;
t272 = t28 * t284;
t271 = t188 * t295;
t268 = t321 * pkin(3);
t265 = t188 * t285;
t264 = qJD(2) * t299;
t262 = g(3) * (t177 * t299 + t191 * t301);
t259 = t196 * t281;
t257 = t195 * t11 + t192 * t26 + t57 * t283 - t41 * t284;
t48 = t106 * t186 + t252;
t110 = t158 * t186 + t191 * t251;
t247 = t132 * t195;
t245 = t193 * t264;
t244 = t187 * pkin(3) * t300 + t134 * t343;
t174 = -t268 - pkin(4);
t242 = g(1) * t134 - g(2) * t136;
t2 = -qJD(6) * t113 + t257 - t334;
t240 = -t132 * t8 + t2;
t176 = pkin(4) + t337;
t190 = -qJ(6) - pkin(9);
t237 = t176 * t179 - t178 * t190;
t235 = t144 * pkin(3);
t199 = qJD(2) ^ 2;
t233 = qJDD(2) * t197 - t194 * t199;
t232 = t195 * t86 + (-t284 - t311) * t132;
t229 = -g(1) * t187 + t322 * g(2);
t145 = t189 * t193 + t194 * t300;
t78 = t186 * t144 + t321 * t145;
t58 = -t192 * t78 - t269;
t228 = -t195 * t78 + t271;
t225 = -t149 * t284 + t294;
t224 = t132 * t40 - t173 * t86;
t222 = -g(1) * (-t134 * t195 + t137 * t307) - g(2) * (t135 * t307 + t136 * t195) - (-t179 * t295 + t194 * t195) * t368;
t221 = -g(1) * (-t137 * t306 - t312) - g(2) * (-t135 * t306 + t136 * t192) - (t179 * t293 + t296) * t368;
t219 = -g(1) * t91 + g(2) * t93 + g(3) * t128;
t217 = t233 * t188;
t215 = t145 * qJD(3);
t155 = -t266 - t333;
t213 = -qJD(2) * t155 - t125 - t242;
t212 = (-t136 * t193 - t196 * t255) * pkin(3);
t211 = -g(1) * (-t137 * t192 + t195 * t91) - g(2) * (-t135 * t192 - t195 * t93) - g(3) * (-t128 * t195 + t271) - t257;
t206 = -pkin(8) * qJDD(3) + (t155 + t266 - t333) * qJD(3);
t89 = -t177 * qJDD(2) + t239 + t357;
t205 = -t215 - t245;
t202 = t209 + t356;
t200 = t141 * t283 + t192 * t280 + t201;
t156 = t174 - t337;
t147 = t292 * t195;
t146 = t292 * t192;
t120 = t135 * t177;
t112 = t113 ^ 2;
t105 = t144 * qJD(3) + t196 * t264;
t81 = t195 * t88;
t77 = -t321 * t144 + t145 * t186;
t69 = pkin(5) * t309 + t110;
t49 = t321 * t105 + t186 * t205;
t47 = t105 * t186 - t321 * t205;
t30 = -pkin(5) * t311 + t48;
t29 = -qJ(6) * t309 + t328;
t27 = -pkin(5) * t223 - qJ(6) * t308 - t111 * t192 + t81;
t21 = t228 * qJD(5) - t192 * t49 + t195 * t265;
t20 = t58 * qJD(5) + t192 * t265 + t195 * t49;
t18 = -t132 ^ 2 * t195 - t315 - t329;
t17 = t232 - t317;
t4 = t113 * t47 + t132 * t21 + t43 * t77 + t58 * t86;
t3 = t115 * t47 - t132 * t20 + t218 * t77 + t228 * t86;
t6 = [t291, 0, t217, (-qJDD(2) * t194 - t197 * t199) * t188, 0, 0, 0, 0, 0, t144 * qJDD(3) + t196 * t217 + (-t215 - 0.2e1 * t245) * qJD(3), -qJD(3) * t105 - qJDD(3) * t145 + (-t193 * t233 - t197 * t259) * t188, -qJD(3) * t47 - qJDD(3) * t77 + (t138 * t285 - t197 * t90) * t188, -t49 * qJD(3) - t78 * qJDD(3) + (t141 * t285 - t197 * t359) * t188, -t49 * t138 + t47 * t141 + t359 * t77 - t78 * t90, -t12 * t77 + t13 * t78 - t45 * t47 + t46 * t49 - g(3) + (t131 * t285 - t197 * t89) * t188, 0, 0, 0, 0, 0, t4, t3, t4, t3, -t113 * t20 - t115 * t21 - t218 * t58 + t228 * t43, t1 * t58 + t15 * t20 - t2 * t228 + t21 * t8 + t28 * t47 + t5 * t77 - g(3); 0, qJDD(2), t214 + t258, -t291 * t301 - t242, qJDD(2) * t184 + 0.2e1 * t193 * t259, 0.2e1 * t193 * t277 - 0.2e1 * t289 * t281, qJDD(3) * t193 + t196 * t198, qJDD(3) * t196 - t193 * t198, 0, t206 * t193 + t196 * t355, -t193 * t355 + t206 * t196, -t138 * t267 - qJDD(3) * t110 + t131 * t140 - t223 * t89 - t177 * t90 + t214 * t179 + (t138 * t343 - t323) * qJD(3), -qJD(3) * t362 - t111 * qJDD(3) + t131 * t143 - t141 * t358 + t89 * t149 - t177 * t359 - t360, -g(3) * t301 + t110 * t359 - t111 * t90 - t12 * t149 + t13 * t223 - t138 * t362 - t46 * t140 + t141 * t323 - t45 * t143 + t242, t13 * t111 - t12 * t110 - t89 * t177 - g(1) * t290 - g(2) * (t136 * t191 - t120) - t262 + t362 * t46 - t323 * t45 - t358 * t131, t115 * t225 + t218 * t308, (-t113 * t195 - t314) * t143 + (-t330 - t195 * t43 + (t113 * t192 - t115 * t195) * qJD(5)) * t149, t115 * t140 + t132 * t225 - t218 * t223 + t86 * t308, -t113 * t140 - t132 * t226 + t223 * t43 - t86 * t309, t132 * t140 - t223 * t86, t81 * t86 - (-t283 * t41 + t25) * t223 + t22 * t140 + t110 * t43 + t40 * t263 + (-t111 * t283 + t365) * t132 + t323 * t113 + ((-qJD(5) * t88 - t71) * t132 - t111 * t86 - (-qJD(5) * t57 - t11) * t223 + t10 * t149 + t40 * t143) * t192 + t221, -t328 * t86 + t257 * t223 - t23 * t140 + t110 * t218 + t40 * t294 + (t10 * t195 - t40 * t284) * t149 + (t111 * t284 - t363) * t132 + t323 * t115 + t222, t28 * t297 - t1 * t223 + t140 * t8 + t27 * t86 + t43 * t69 + (t28 * t283 + t338) * t149 + t346 * t132 + t324 * t113 + t221, t28 * t294 - t140 * t15 + t223 * t2 - t29 * t86 + t218 * t69 + (t195 * t5 - t272) * t149 - t345 * t132 + t324 * t115 + t222, -t27 * t218 - t29 * t43 + (-t15 * t192 - t195 * t8) * t143 - t346 * t115 - t345 * t113 + t360 + (-t1 * t195 - t192 * t2 + (-t15 * t195 + t192 * t8) * qJD(5)) * t149, t2 * t29 + t1 * t27 + t5 * t69 - g(1) * (-pkin(5) * t312 - t137 * t237 + t290) - g(2) * (-t120 + (pkin(5) * t192 + t191) * t136 - t237 * t135) - t262 + t346 * t8 + t324 * t28 - (pkin(5) * t296 + t197 * t237) * t368 + t345 * t15; 0, 0, 0, 0, -t193 * t199 * t196, t289 * t199, t278, t277, qJDD(3), -g(3) * t144 + t193 * t213 + t229 * t300 + t165, g(3) * t145 + (-t188 * t229 - t279) * t193 + t213 * t196, t48 * qJD(3) - t131 * t141 + (t321 * qJDD(3) - t138 * t286) * pkin(3) + t220 + t12, qJD(3) * t50 + t131 * t138 + (-qJDD(3) * t186 - t141 * t286) * pkin(3) + t219 - t13, -t359 * t268 - t90 * t344 - (-t46 + t48) * t141 + (-t45 + t50) * t138, -g(1) * t244 - g(2) * t212 - g(3) * t235 + t12 * t268 + t13 * t344 - t131 * t274 + t45 * t48 - t46 * t50, t115 * t247 + t330, (t218 - t318) * t195 - t361 + t327, t132 * t247 - t315 + t329, t232 + t317, -t132 * t141, -t48 * t113 - t63 * t132 - t22 * t141 + t174 * t43 + (t50 * t132 + t224) * t192 + (-qJD(5) * t132 * t173 - t10 + t220) * t195, t10 * t192 - t48 * t115 + t23 * t141 + t174 * t218 + (t173 * t284 + t336) * t132 + t224 * t195 - t275, -t113 * t30 - t141 * t8 - t146 * t86 + t156 * t43 + t325 * t132 + (t138 * t28 + (t28 + t339) * qJD(5)) * t192 + t366 * t195, t28 * t310 - t115 * t30 + t141 * t15 - t147 * t86 + t156 * t218 + t338 - t364 * t132 + (pkin(5) * t314 + t195 * t28) * qJD(5) - t275, -t113 * t364 - t325 * t115 + t146 * t218 - t147 * t43 - t192 * t367 + t240 * t195 - t219, t2 * t147 - t1 * t146 + t5 * t156 + pkin(5) * t272 - t28 * t30 - g(1) * (t176 * t94 + t190 * t91 + t244) - g(2) * (-t92 * t176 - t93 * t190 + t212) - g(3) * (-t127 * t176 - t128 * t190 + t235) + t325 * t8 + t364 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t141 * qJD(3) + t238, (t166 - t138) * qJD(3) + t210, -t138 ^ 2 - t141 ^ 2, t138 * t46 + t141 * t45 - t214 + t89, 0, 0, 0, 0, 0, t17, t18, t17, t18, (-t218 - t318) * t195 + t361 + t327, -t141 * t28 + t240 * t192 + t195 * t367 - t214; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115 * t113, -t112 + t354, t218 + t319, -t200 + t316, t86, -t115 * t40 + t132 * t23 + t202, t113 * t40 + t132 * t22 + t211, 0.2e1 * t353 - t335 + t331 + (-t249 - t28) * t115 + t202, -pkin(5) * t354 + t334 + t132 * t14 + (qJD(6) + t28) * t113 + t211, -pkin(5) * t218 + t347 * t113, -t347 * t15 + (-t28 * t115 + t1 + t356) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t200 + t316, t218 - t319, -t112 - t354, t15 * t113 + t8 * t115 - t366;];
tau_reg = t6;
