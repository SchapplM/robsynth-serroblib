% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:06:00
% EndTime: 2020-01-03 12:06:07
% DurationCPUTime: 3.96s
% Computational Cost: add. (5893->404), mult. (8918->546), div. (0->0), fcn. (5766->14), ass. (0->246)
t185 = qJDD(1) + qJDD(2);
t188 = qJD(1) + qJD(2);
t193 = sin(pkin(9));
t199 = cos(qJ(5));
t200 = cos(qJ(4));
t294 = t193 * t200;
t252 = t185 * t294;
t195 = sin(qJ(5));
t196 = sin(qJ(4));
t288 = t195 * t196;
t287 = t195 * t200;
t131 = t196 * t199 + t287;
t267 = qJD(4) + qJD(5);
t335 = t267 * t131;
t31 = (t185 * t288 + t188 * t335) * t193 - t199 * t252;
t197 = sin(qJ(2));
t194 = cos(pkin(9));
t272 = qJD(3) * t194;
t201 = cos(qJ(2));
t289 = t194 * t201;
t319 = pkin(1) * qJD(1);
t138 = -pkin(3) * t194 - pkin(7) * t193 - pkin(2);
t290 = t194 * t200;
t100 = qJ(3) * t290 + t196 * t138;
t333 = qJD(4) * t100;
t314 = -t196 * t272 - (-t196 * t289 + t197 * t200) * t319 - t333;
t270 = qJD(4) * t200;
t339 = (t196 * t197 + t200 * t289) * t319 - t138 * t270 - t200 * t272;
t271 = qJD(4) * t196;
t247 = t193 * t271;
t155 = pkin(8) * t247;
t337 = t155 + t314;
t291 = t194 * t196;
t258 = qJ(3) * t291;
t264 = pkin(8) * t294;
t336 = -(-t258 - t264) * qJD(4) + t339;
t275 = qJD(2) * t197;
t262 = pkin(1) * t275;
t330 = pkin(1) * t201;
t279 = -qJD(1) * t262 + qJDD(1) * t330;
t245 = qJDD(3) - t279;
t326 = t185 * pkin(2);
t102 = t245 - t326;
t192 = qJ(1) + qJ(2);
t179 = sin(t192);
t181 = cos(t192);
t281 = g(2) * t181 + g(3) * t179;
t220 = -t102 - t281;
t253 = t193 * t288;
t334 = -qJD(5) * t253 - t195 * t247;
t276 = qJD(1) * t201;
t261 = pkin(1) * t276;
t227 = qJD(3) - t261;
t263 = t197 * t319;
t134 = qJ(3) * t188 + t263;
t186 = t193 ^ 2;
t332 = t134 * (qJD(4) * t194 + t186 * t188);
t118 = t138 - t330;
t168 = pkin(1) * t197 + qJ(3);
t73 = t196 * t118 + t168 * t290;
t106 = -t179 * t291 - t181 * t200;
t108 = -t179 * t200 + t181 * t291;
t331 = -g(2) * t106 - g(3) * t108;
t292 = t194 * t185;
t142 = -qJDD(4) + t292;
t65 = t138 * t185 + t245;
t62 = t200 * t65;
t268 = qJDD(1) * t197;
t274 = qJD(2) * t201;
t90 = t185 * qJ(3) + t188 * qJD(3) + (qJD(1) * t274 + t268) * pkin(1);
t228 = -t291 * t90 + t62;
t251 = t134 * t290;
t300 = t188 * t193;
t81 = t138 * t188 + t227;
t18 = -pkin(8) * t252 - t142 * pkin(4) + (-t251 + (pkin(8) * t300 - t81) * t196) * qJD(4) + t228;
t248 = t188 * t270;
t303 = t185 * t196;
t218 = t248 + t303;
t246 = t194 * t271;
t266 = t196 * t65 + t81 * t270 + t90 * t290;
t25 = -t134 * t246 + t266;
t19 = -pkin(8) * t193 * t218 + t25;
t269 = qJD(5) * t195;
t299 = t188 * t194;
t145 = -qJD(4) + t299;
t254 = t188 * t294;
t52 = -t134 * t291 + t200 * t81;
t44 = -pkin(8) * t254 + t52;
t34 = -pkin(4) * t145 + t44;
t295 = t193 * t196;
t265 = pkin(8) * t295;
t53 = t196 * t81 + t251;
t45 = -t188 * t265 + t53;
t4 = t199 * (qJD(5) * t34 + t19) + t195 * t18 - t45 * t269;
t329 = g(1) * t193;
t216 = t188 * t131;
t85 = t193 * t216;
t88 = -t188 * t253 + t199 * t254;
t325 = t88 * t85;
t124 = t200 * t138;
t67 = -t264 + t124 + (-qJ(3) * t196 - pkin(4)) * t194;
t80 = -t265 + t100;
t37 = -t195 * t80 + t199 * t67;
t324 = qJD(5) * t37 + t337 * t195 - t336 * t199;
t38 = t195 * t67 + t199 * t80;
t323 = -qJD(5) * t38 + t336 * t195 + t337 * t199;
t160 = pkin(1) * t274 + qJD(3);
t312 = t134 * t186;
t82 = t186 * t90;
t322 = t160 * t312 + t168 * t82;
t130 = t199 * t200 - t288;
t321 = (t267 - t299) * t130;
t320 = t194 * t216 - t335;
t318 = t195 * t45;
t317 = t199 * t45;
t273 = qJD(3) * t134;
t316 = qJ(3) * t82 + t186 * t273;
t315 = -qJ(3) * t246 - t339;
t313 = qJD(4) * t53;
t135 = -qJDD(5) + t142;
t311 = t135 * t194;
t310 = t142 * t194;
t309 = t160 * t188;
t308 = t179 * t193;
t307 = t179 * t194;
t306 = t179 * t196;
t305 = t181 * t193;
t304 = t181 * t194;
t302 = t185 * t200;
t184 = t188 ^ 2;
t301 = t186 * t184;
t298 = t188 * t196;
t297 = t188 * t200;
t296 = t193 * t185;
t293 = t193 * (-pkin(8) - pkin(7));
t286 = t196 * t200;
t284 = t334 * t188;
t283 = g(2) * t305 + g(3) * t308;
t282 = t181 * pkin(2) + t179 * qJ(3);
t280 = g(2) * t179 - g(3) * t181;
t187 = t194 ^ 2;
t278 = t186 + t187;
t189 = t196 ^ 2;
t190 = t200 ^ 2;
t277 = t189 - t190;
t260 = t52 * t247 - t283;
t259 = t102 * t193 + t283;
t257 = t168 * t291;
t255 = t185 * t287;
t250 = t118 * t270 + t160 * t290 + t196 * t262;
t249 = t188 * t275;
t242 = t201 * t278;
t241 = t278 * t185;
t169 = t179 * pkin(2);
t240 = -qJ(3) * t181 + t169;
t239 = t142 + t292;
t238 = t188 * (-qJD(4) - t145);
t237 = t187 * t90 - t280 + t82;
t236 = t267 * t200;
t235 = t188 * t263;
t234 = t286 * t301;
t233 = pkin(3) * t304 + pkin(7) * t305 + t282;
t232 = t279 - t281;
t231 = t196 * t248;
t198 = sin(qJ(1));
t202 = cos(qJ(1));
t226 = -g(2) * t202 - g(3) * t198;
t21 = t195 * t34 + t317;
t114 = t200 * t118;
t56 = -t264 + t114 + (-t168 * t196 - pkin(4)) * t194;
t63 = -t265 + t73;
t29 = -t195 * t63 + t199 * t56;
t30 = t195 * t56 + t199 * t63;
t225 = t196 * t52 - t200 * t53;
t111 = t131 * t193;
t112 = t130 * t193;
t20 = t199 * t34 - t318;
t5 = -qJD(5) * t21 + t199 * t18 - t195 * t19;
t57 = t335 * t193;
t58 = t193 * t199 * t236 + t334;
t224 = -t4 * t111 - t5 * t112 + t20 * t57 - t21 * t58 - t283;
t223 = qJD(4) * (t145 + t299);
t156 = t193 * pkin(4) * t270;
t221 = t227 * t193 + t156;
t219 = pkin(3) * t307 + pkin(7) * t308 + t240;
t217 = g(2) * t108 - g(3) * t106 + t25 * t194 + t200 * t82;
t215 = -t235 - t326;
t51 = (pkin(4) * t218 + t90) * t193;
t84 = (pkin(4) * t298 + t134) * t193;
t191 = qJ(4) + qJ(5);
t178 = sin(t191);
t180 = cos(t191);
t94 = -t178 * t307 - t180 * t181;
t96 = t178 * t304 - t179 * t180;
t214 = g(2) * t96 - g(3) * t94 + t51 * t112 + t4 * t194 - t84 * t57;
t176 = -pkin(2) - t330;
t213 = pkin(1) * t249 + t176 * t185;
t175 = pkin(4) * t200 + pkin(3);
t212 = pkin(4) * t306 + t175 * t304 - t181 * t293 + t282;
t95 = -t181 * t178 + t180 * t307;
t97 = t178 * t179 + t180 * t304;
t211 = -g(2) * t97 - g(3) * t95 + t51 * t111 - t194 * t5 + t84 * t58;
t210 = -t145 ^ 2 - t301;
t107 = t179 * t290 - t181 * t196;
t109 = t181 * t290 + t306;
t26 = t228 - t313;
t209 = -g(2) * t109 - g(3) * t107 - t26 * t194 + t196 * t82 + t270 * t312;
t208 = -t179 * t293 + t175 * t307 + t169 + (-pkin(4) * t196 - qJ(3)) * t181;
t50 = -t73 * qJD(4) - t160 * t291 + t200 * t262;
t206 = g(2) * t95 - g(3) * t97 + t180 * t329 + t84 * t85 - t4;
t205 = -g(2) * t94 - g(3) * t96 + t178 * t329 - t84 * t88 + t5;
t183 = t202 * pkin(1);
t182 = t198 * pkin(1);
t165 = pkin(4) * t295;
t164 = t187 * t185;
t163 = t186 * t185;
t137 = -qJD(5) + t145;
t136 = 0.2e1 * t193 * t292;
t128 = qJ(3) * t193 + t165;
t127 = -pkin(2) * t188 + t227;
t115 = t168 * t193 + t165;
t101 = t160 * t193 + t156;
t99 = t124 - t258;
t92 = (t185 * t190 - 0.2e1 * t231) * t186;
t91 = (t185 * t189 + 0.2e1 * t231) * t186;
t72 = t114 - t257;
t64 = 0.2e1 * (qJD(4) * t188 * t277 - t185 * t286) * t186;
t49 = -t168 * t246 + t250;
t47 = (t196 * t239 + t200 * t223) * t193;
t46 = (t196 * t223 - t200 * t239) * t193;
t43 = t155 + t50;
t42 = (-t257 - t264) * qJD(4) + t250;
t33 = -t85 ^ 2 + t88 ^ 2;
t32 = (t255 + (t188 * t236 + t303) * t199) * t193 + t284;
t28 = -t88 * t137 + (-t255 + (-t267 * t297 - t303) * t199) * t193 - t284;
t27 = -t85 * t137 - t31;
t23 = t199 * t44 - t318;
t22 = -t195 * t44 - t317;
t14 = t111 * t32 + t58 * t85;
t13 = -t112 * t31 - t57 * t88;
t12 = t111 * t135 + t137 * t58 + t194 * t32;
t11 = -t112 * t135 + t137 * t57 + t194 * t31;
t8 = -qJD(5) * t30 - t195 * t42 + t199 * t43;
t7 = qJD(5) * t29 + t195 * t43 + t199 * t42;
t6 = t111 * t31 - t112 * t32 + t57 * t85 - t58 * t88;
t1 = [0, 0, 0, 0, 0, qJDD(1), t226, g(2) * t198 - g(3) * t202, 0, 0, 0, 0, 0, 0, 0, t185, (t185 * t201 - t249) * pkin(1) + t232, ((-qJDD(1) - t185) * t197 + (-qJD(1) - t188) * t274) * pkin(1) + t280, 0, (t226 + (t197 ^ 2 + t201 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t163, t136, 0, t164, 0, 0, (-t213 + t220) * t194, t193 * t213 + t259, t168 * t241 + t278 * t309 + t237, t102 * t176 + t127 * t262 - g(2) * (t183 + t282) - g(3) * (t182 + t240) + (t134 * t160 + t168 * t90) * t187 + t322, t92, t64, t46, t91, t47, t310, -t72 * t142 - t50 * t145 + (t160 * t298 + t168 * t218) * t186 + t209, t73 * t142 + t49 * t145 + ((t168 * t185 + t309) * t200 + (-t168 * t188 - t134) * t271) * t186 + t217, ((-t185 * t73 - t25 + (qJD(4) * t72 - t49) * t188) * t196 + (-t185 * t72 - t188 * t50 - t26 + (-t188 * t73 - t53) * qJD(4)) * t200) * t193 + t260, t25 * t73 + t53 * t49 + t26 * t72 + t52 * t50 - g(2) * (t183 + t233) - g(3) * (t182 + t219) + t322, t13, t6, t11, t14, t12, t311, t101 * t85 + t115 * t32 - t135 * t29 - t137 * t8 + t211, t101 * t88 - t115 * t31 + t135 * t30 + t137 * t7 + t214, t29 * t31 - t30 * t32 - t7 * t85 - t8 * t88 + t224, t4 * t30 + t21 * t7 + t5 * t29 + t20 * t8 + t51 * t115 + t84 * t101 - g(2) * (t183 + t212) - g(3) * (t182 + t208); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t185, t232 + t235, (-t268 + (-qJD(2) + t188) * t276) * pkin(1) + t280, 0, 0, t163, t136, 0, t164, 0, 0, (-t215 + t220) * t194, t193 * t215 + t259, qJ(3) * t241 + (qJD(3) * t278 - t242 * t319) * t188 + t237, -t102 * pkin(2) - g(2) * t282 - g(3) * t240 + (qJ(3) * t90 + t273) * t187 + (-t127 * t197 - t134 * t242) * t319 + t316, t92, t64, t46, t91, t47, t310, -t142 * t99 - t314 * t145 + (qJ(3) * t303 + (qJ(3) * t270 + t196 * t227) * t188) * t186 + t209, t100 * t142 + t315 * t145 + (qJ(3) * t302 - t134 * t271 + (-qJ(3) * t271 + t200 * t227) * t188) * t186 + t217, ((-t100 * t185 - t25) * t196 + (-t185 * t99 - t26 - t313) * t200 + ((-t314 - t333) * t200 + (qJD(4) * t99 - t315) * t196) * t188) * t193 + t260, -g(2) * t233 - g(3) * t219 + t25 * t100 + t26 * t99 - t261 * t312 + t314 * t52 + t315 * t53 + t316, t13, t6, t11, t14, t12, t311, t128 * t32 - t135 * t37 - t323 * t137 + t221 * t85 + t211, -t128 * t31 + t38 * t135 + t324 * t137 + t221 * t88 + t214, t31 * t37 - t32 * t38 - t323 * t88 - t324 * t85 + t224, -g(2) * t212 - g(3) * t208 + t51 * t128 + t323 * t20 + t324 * t21 + t221 * t84 + t5 * t37 + t4 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t292, t296, -t278 * t184, -t134 * t188 * t278 - t220, 0, 0, 0, 0, 0, 0, -t200 * t142 + t196 * t210, t196 * t142 + t200 * t210, (-t189 - t190) * t296, t25 * t196 + t26 * t200 - t225 * qJD(4) + (t194 * t225 - t312) * t188 + t281, 0, 0, 0, 0, 0, 0, -t130 * t135 - t320 * t137 - t85 * t300, t131 * t135 + t321 * t137 - t88 * t300, t130 * t31 - t131 * t32 - t320 * t88 - t321 * t85, t5 * t130 + t4 * t131 + t320 * t20 + t321 * t21 - t84 * t300 + t281; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t234, -t277 * t301, (t196 * t238 + t302) * t193, -t234, (t200 * t238 - t303) * t193, -t142, -t53 * t145 + t62 - t200 * t332 + (-qJD(4) * t81 - t194 * t90 + t329) * t196 + t331, g(1) * t294 + g(2) * t107 - g(3) * t109 - t52 * t145 + t196 * t332 - t266, 0, 0, t325, t33, t27, -t325, t28, -t135, t137 * t22 + (-t135 * t199 + t137 * t269 - t254 * t85) * pkin(4) + t205, -t137 * t23 + (qJD(5) * t137 * t199 + t135 * t195 - t254 * t88) * pkin(4) + t206, (t21 + t22) * t88 + (-t20 + t23) * t85 + (-t195 * t32 + t199 * t31 + (t195 * t88 - t199 * t85) * qJD(5)) * pkin(4), -t20 * t22 - t21 * t23 + (t4 * t195 + t5 * t199 + (g(1) * t196 - t297 * t84) * t193 + (-t195 * t20 + t199 * t21) * qJD(5) + t331) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t325, t33, t27, -t325, t28, -t135, -t21 * t137 + t205, -t20 * t137 + t206, 0, 0;];
tau_reg = t1;
