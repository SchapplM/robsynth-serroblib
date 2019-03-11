% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
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
% tau_reg [6x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRRRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:59:21
% EndTime: 2019-03-08 23:59:34
% DurationCPUTime: 5.44s
% Computational Cost: add. (5683->448), mult. (12800->614), div. (0->0), fcn. (10011->14), ass. (0->252)
t188 = sin(qJ(3));
t189 = sin(qJ(2));
t183 = sin(pkin(6));
t288 = qJD(1) * t183;
t262 = t189 * t288;
t320 = qJD(3) * pkin(3);
t337 = -t188 * t320 + t262;
t332 = cos(qJ(4));
t255 = qJD(4) * t332;
t191 = cos(qJ(3));
t193 = -pkin(9) - pkin(8);
t244 = -qJD(2) * t193 + t262;
t184 = cos(pkin(6));
t287 = qJD(1) * t184;
t108 = -t244 * t188 + t191 * t287;
t109 = t188 * t287 + t191 * t244;
t187 = sin(qJ(4));
t98 = t187 * t109;
t56 = t108 * t332 - t98;
t348 = -pkin(3) * t255 + t56;
t263 = t332 * t191;
t219 = -t187 * t188 + t263;
t274 = qJD(3) + qJD(4);
t103 = t274 * t219;
t294 = t187 * t191;
t144 = t188 * t332 + t294;
t104 = t274 * t144;
t350 = pkin(4) * t104 - pkin(10) * t103 - t337;
t192 = cos(qJ(2));
t261 = t192 * t288;
t112 = t219 * t261;
t264 = qJD(3) * t193;
t145 = t188 * t264;
t146 = t191 * t264;
t155 = t193 * t188;
t156 = t193 * t191;
t220 = t155 * t332 + t187 * t156;
t58 = qJD(4) * t220 + t145 * t332 + t187 * t146;
t349 = t58 - t112;
t186 = sin(qJ(5));
t190 = cos(qJ(5));
t347 = t112 * t186 + t190 * t350;
t100 = t108 + t320;
t283 = qJD(4) * t187;
t277 = t184 * qJDD(1);
t160 = t191 * t277;
t279 = qJD(1) * qJD(2);
t123 = qJDD(2) * pkin(8) + (qJDD(1) * t189 + t192 * t279) * t183;
t243 = pkin(9) * qJDD(2) + t123;
t42 = qJDD(3) * pkin(3) - qJD(3) * t109 - t243 * t188 + t160;
t48 = qJD(3) * t108 + t188 * t277 + t243 * t191;
t245 = t100 * t283 + t109 * t255 + t187 * t48 - t332 * t42;
t273 = qJDD(3) + qJDD(4);
t11 = -pkin(4) * t273 + t245;
t181 = qJ(3) + qJ(4);
t176 = sin(t181);
t177 = cos(t181);
t298 = t183 * t189;
t124 = t176 * t298 - t184 * t177;
t312 = cos(pkin(11));
t247 = t312 * t189;
t182 = sin(pkin(11));
t299 = t182 * t192;
t133 = t184 * t247 + t299;
t248 = t183 * t312;
t92 = t133 * t176 + t177 * t248;
t246 = t312 * t192;
t300 = t182 * t189;
t135 = -t184 * t300 + t246;
t301 = t182 * t183;
t94 = t135 * t176 - t177 * t301;
t216 = g(1) * t94 + g(2) * t92 + g(3) * t124;
t212 = -t11 + t216;
t280 = qJD(5) * t190;
t174 = t191 * pkin(3) + pkin(2);
t90 = -pkin(4) * t219 - pkin(10) * t144 - t174;
t346 = t186 * t350 + t190 * t349 + t90 * t280;
t256 = qJD(2) * t332;
t285 = qJD(2) * t188;
t138 = t187 * t285 - t191 * t256;
t307 = t138 * t186;
t345 = -qJ(6) * t307 + t190 * qJD(6);
t140 = -qJD(2) * t294 - t188 * t256;
t89 = -pkin(4) * t140 + pkin(10) * t138;
t73 = pkin(3) * t285 + t89;
t344 = t186 * t73 + t190 * t348;
t119 = t187 * t155 - t156 * t332;
t107 = t190 * t119;
t228 = -qJ(6) * t103 - qJD(6) * t144;
t343 = pkin(5) * t104 - t186 * t58 + t228 * t190 + (-t107 + (qJ(6) * t144 - t90) * t186) * qJD(5) + t347;
t257 = t144 * t280;
t342 = -qJ(6) * t257 + (-qJD(5) * t119 + t228) * t186 + t346;
t171 = t187 * pkin(3) + pkin(10);
t292 = -qJ(6) - t171;
t242 = qJD(5) * t292;
t341 = t186 * t242 - t344 + t345;
t178 = t190 * qJ(6);
t235 = -t140 * pkin(5) + t138 * t178;
t70 = t190 * t73;
t340 = t190 * t242 - t235 - t70 + (-qJD(6) + t348) * t186;
t99 = t332 * t109;
t55 = t187 * t108 + t99;
t238 = pkin(3) * t283 - t55;
t296 = t183 * t192;
t208 = t144 * t296;
t313 = -qJD(1) * t208 + qJD(4) * t119 + t187 * t145 - t146 * t332;
t132 = -t184 * t246 + t300;
t134 = t184 * t299 + t247;
t237 = g(1) * t134 + g(2) * t132;
t211 = -g(3) * t296 + t237;
t339 = t211 * t176;
t281 = qJD(5) * t186;
t338 = (t281 + t307) * pkin(5);
t197 = t103 * qJD(2) + t144 * qJDD(2);
t254 = t189 * t279;
t234 = -qJDD(1) * t296 + t183 * t254;
t311 = qJDD(2) * pkin(2);
t122 = t234 - t311;
t194 = qJD(3) ^ 2;
t336 = -pkin(8) * t194 + t183 * (-g(3) * t192 + t254) - t122 + t237 + t311;
t213 = t190 * t140 - t186 * t274;
t36 = -qJD(5) * t213 + t186 * t197 - t190 * t273;
t335 = t213 ^ 2;
t328 = g(3) * t183;
t327 = t190 * pkin(5);
t185 = -qJ(6) - pkin(10);
t131 = qJD(5) + t138;
t54 = t187 * t100 + t99;
t46 = pkin(10) * t274 + t54;
t127 = -qJD(2) * t174 - t261;
t66 = pkin(4) * t138 + pkin(10) * t140 + t127;
t28 = -t186 * t46 + t190 * t66;
t18 = qJ(6) * t213 + t28;
t14 = pkin(5) * t131 + t18;
t326 = -t18 + t14;
t53 = t100 * t332 - t98;
t325 = t186 * t89 + t190 * t53;
t172 = pkin(4) + t327;
t93 = t133 * t177 - t176 * t248;
t323 = -t92 * t172 - t93 * t185;
t95 = t135 * t177 + t176 * t301;
t322 = -t94 * t172 - t95 * t185;
t321 = qJD(2) * pkin(2);
t240 = t190 * t274;
t35 = -qJD(5) * t240 - t140 * t281 - t186 * t273 - t190 * t197;
t319 = t186 * t35;
t276 = t188 * qJDD(2);
t233 = -qJDD(2) * t263 + t187 * t276;
t67 = qJD(2) * t104 + t233;
t63 = qJDD(5) + t67;
t318 = t190 * t63;
t45 = -pkin(4) * t274 - t53;
t317 = t45 * t138;
t316 = t186 * t90 + t107;
t249 = qJD(5) * t185;
t315 = t186 * t249 - t325 + t345;
t80 = t190 * t89;
t314 = t190 * t249 - t235 - t80 + (-qJD(6) + t53) * t186;
t113 = -t140 * t186 - t240;
t310 = t113 * t131;
t309 = t213 * t131;
t308 = t131 * t140;
t306 = t138 * t190;
t305 = t140 * t138;
t304 = t144 * t186;
t303 = t144 * t190;
t302 = t177 * t186;
t297 = t183 * t191;
t295 = t186 * t192;
t293 = t190 * t103;
t291 = qJDD(1) - g(3);
t125 = t176 * t184 + t177 * t298;
t290 = -t124 * t172 - t125 * t185;
t179 = t188 ^ 2;
t289 = -t191 ^ 2 + t179;
t286 = qJD(2) * t183;
t284 = qJD(2) * t189;
t282 = qJD(5) * t131;
t278 = qJD(2) * qJD(3);
t275 = t191 * qJDD(2);
t266 = t183 * t295;
t265 = t190 * t296;
t38 = t45 * t281;
t259 = t183 * t284;
t258 = t192 * t286;
t253 = t188 * t278;
t252 = t192 * t278;
t251 = pkin(5) * t186 - t193;
t250 = t28 * t140 + t38;
t241 = t131 * t190;
t173 = -pkin(3) * t332 - pkin(4);
t236 = g(1) * t135 + g(2) * t133;
t232 = -t171 * t63 + t317;
t29 = t186 * t66 + t190 * t46;
t19 = -qJ(6) * t113 + t29;
t231 = -t14 * t190 - t186 * t19;
t195 = qJD(2) ^ 2;
t227 = qJDD(2) * t192 - t189 * t195;
t226 = -g(1) * t182 + g(2) * t312;
t136 = t184 * t191 - t188 * t298;
t137 = t184 * t188 + t189 * t297;
t77 = t187 * t136 + t137 * t332;
t64 = -t186 * t77 - t265;
t225 = -t190 * t77 + t266;
t204 = t100 * t255 - t109 * t283 + t187 * t42 + t332 * t48;
t10 = pkin(10) * t273 + t204;
t168 = pkin(3) * t253;
t26 = -pkin(3) * t275 + t67 * pkin(4) - pkin(10) * t197 + t122 + t168;
t224 = -t190 * t10 - t186 * t26 - t66 * t280 + t281 * t46;
t222 = t172 * t177 - t176 * t185 + t174;
t221 = t136 * t332 - t187 * t137;
t218 = t103 * t186 + t257;
t217 = -t144 * t281 + t293;
t215 = g(1) * t95 + g(2) * t93 + g(3) * t125;
t214 = -t29 * t140 - t186 * t212 + t45 * t280;
t148 = -t261 - t321;
t209 = -qJD(2) * t148 - t123 + t236;
t207 = t211 * t177;
t25 = t190 * t26;
t206 = -qJD(5) * t29 - t186 * t10 + t25;
t5 = t36 * pkin(5) + qJDD(6) + t11;
t203 = -pkin(8) * qJDD(3) + (t148 + t261 - t321) * qJD(3);
t202 = t127 * t140 + t216 - t245;
t1 = pkin(5) * t63 + qJ(6) * t35 + qJD(6) * t213 + t206;
t3 = -qJ(6) * t36 - qJD(6) * t113 - t224;
t201 = qJD(5) * t231 - t1 * t186 - t14 * t306 - t19 * t307 + t3 * t190 - t215;
t199 = t127 * t138 - t204 + t215;
t198 = -g(1) * (t134 * t190 - t186 * t95) - g(2) * (t132 * t190 - t186 * t93) - g(3) * (-t125 * t186 - t265);
t152 = pkin(10) * t190 + t178;
t151 = t185 * t186;
t142 = t171 * t190 + t178;
t141 = t292 * t186;
t110 = t113 ^ 2;
t106 = -qJD(3) * t137 - t188 * t258;
t105 = qJD(3) * t136 + t191 * t258;
t91 = -qJDD(2) * t174 + t168 + t234;
t84 = t190 * t90;
t68 = -t138 ^ 2 + t140 ^ 2;
t50 = -t233 + (-qJD(2) * t144 - t140) * t274;
t49 = t138 * t274 + t197;
t37 = -qJ(6) * t304 + t316;
t34 = -pkin(5) * t219 - t119 * t186 - t144 * t178 + t84;
t33 = t113 * pkin(5) + qJD(6) + t45;
t32 = qJD(4) * t77 + t187 * t105 - t106 * t332;
t31 = qJD(4) * t221 + t105 * t332 + t187 * t106;
t21 = t131 * t241 - t140 * t213 + t186 * t63;
t20 = -t131 ^ 2 * t186 - t113 * t140 + t318;
t17 = -t213 * t241 - t319;
t16 = qJD(5) * t225 - t186 * t31 + t190 * t259;
t15 = qJD(5) * t64 + t186 * t259 + t190 * t31;
t4 = (-t35 - t310) * t190 + (-t36 + t309) * t186;
t2 = [t291, 0, t227 * t183 (-qJDD(2) * t189 - t192 * t195) * t183, 0, 0, 0, 0, 0, qJD(3) * t106 + qJDD(3) * t136 + (-t188 * t252 + t191 * t227) * t183, -qJD(3) * t105 - qJDD(3) * t137 + (-t188 * t227 - t191 * t252) * t183, 0, 0, 0, 0, 0, -t32 * t274 + t221 * t273 + (t138 * t284 - t192 * t67) * t183, -t31 * t274 - t77 * t273 - qJDD(2) * t208 + (-t103 * t192 - t189 * t140) * t286, 0, 0, 0, 0, 0, t113 * t32 + t131 * t16 - t221 * t36 + t63 * t64, -t131 * t15 - t213 * t32 + t221 * t35 + t225 * t63, -t113 * t15 + t16 * t213 + t225 * t36 + t35 * t64, t1 * t64 + t14 * t16 + t15 * t19 - t221 * t5 - t225 * t3 + t32 * t33 - g(3); 0, qJDD(2), t291 * t296 + t237, -t291 * t298 + t236, qJDD(2) * t179 + 0.2e1 * t191 * t253, 0.2e1 * t188 * t275 - 0.2e1 * t278 * t289, qJDD(3) * t188 + t191 * t194, qJDD(3) * t191 - t188 * t194, 0, t203 * t188 + t191 * t336, -t188 * t336 + t203 * t191, -t140 * t103 + t144 * t197, -t103 * t138 + t140 * t104 - t144 * t67 + t197 * t219, t103 * t274 + t144 * t273, -t104 * t274 + t219 * t273, 0, t127 * t104 - t138 * t337 - t174 * t67 - t219 * t91 + t220 * t273 - t274 * t313 + t207, t127 * t103 - t119 * t273 + t337 * t140 + t91 * t144 - t174 * t197 - t274 * t349 - t339, -t213 * t217 - t303 * t35 (-t113 * t190 + t186 * t213) * t103 + (t319 - t190 * t36 + (t113 * t186 + t190 * t213) * qJD(5)) * t144, -t104 * t213 + t131 * t217 + t219 * t35 + t303 * t63, -t104 * t113 - t131 * t218 + t219 * t36 - t304 * t63, t104 * t131 - t219 * t63, t28 * t104 - t220 * t36 - t25 * t219 + t84 * t63 + t347 * t131 + t313 * t113 + (t207 + (-t119 * t131 + t144 * t45 + t219 * t46) * qJD(5)) * t190 + ((-qJD(5) * t90 - t58) * t131 - t119 * t63 - (-qJD(5) * t66 - t10) * t219 + t11 * t144 + t45 * t103 - g(3) * t298 - t236) * t186, -t316 * t63 - t224 * t219 - t29 * t104 + t220 * t35 + t45 * t293 - g(1) * (t134 * t302 + t135 * t190) - g(2) * (t132 * t302 + t133 * t190) - (-t177 * t295 + t189 * t190) * t328 + (t11 * t190 - t38) * t144 + (t119 * t281 - t346) * t131 - t313 * t213, t34 * t35 - t36 * t37 + t343 * t213 - t342 * t113 + t231 * t103 + t339 + (-t1 * t190 - t186 * t3 + (t14 * t186 - t19 * t190) * qJD(5)) * t144, t3 * t37 + t1 * t34 + t5 * (pkin(5) * t304 - t220) - g(1) * (-t134 * t222 + t135 * t251) - g(2) * (-t132 * t222 + t133 * t251) - (t189 * t251 + t192 * t222) * t328 + (pkin(5) * t218 + t313) * t33 + t342 * t19 + t343 * t14; 0, 0, 0, 0, -t188 * t195 * t191, t289 * t195, t276, t275, qJDD(3), -g(3) * t136 + t188 * t209 + t226 * t297 + t160, g(3) * t137 + (-t183 * t226 - t277) * t188 + t209 * t191, -t305, t68, t49, t50, t273, t55 * t274 + (-t138 * t285 + t273 * t332 - t274 * t283) * pkin(3) + t202, t56 * t274 + (t140 * t285 - t187 * t273 - t255 * t274) * pkin(3) + t199, t17, t4, t21, t20, t308, -t70 * t131 + t173 * t36 + t238 * t113 + (t131 * t348 + t232) * t186 + (-t171 * t282 + t212) * t190 + t250, -t173 * t35 + t232 * t190 - t238 * t213 + (t171 * t281 + t344) * t131 + t214, -t113 * t341 + t141 * t35 - t142 * t36 + t213 * t340 + t201, t3 * t142 + t1 * t141 + t5 * (t173 - t327) - g(1) * ((-t135 * t188 + t182 * t297) * pkin(3) + t322) - g(2) * ((-t133 * t188 - t191 * t248) * pkin(3) + t323) - g(3) * (pkin(3) * t136 + t290) + (t338 + t238) * t33 + t341 * t19 + t340 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t305, t68, t49, t50, t273, t274 * t54 + t202, t274 * t53 + t199, t17, t4, t21, t20, t308, -pkin(4) * t36 - t54 * t113 - t80 * t131 + (-pkin(10) * t63 + t53 * t131 + t317) * t186 + (-pkin(10) * t282 + t212) * t190 + t250, pkin(4) * t35 + t325 * t131 + t54 * t213 + t45 * t306 + (t131 * t281 - t318) * pkin(10) + t214, -t113 * t315 + t151 * t35 - t152 * t36 + t213 * t314 + t201, t3 * t152 + t1 * t151 - t5 * t172 - g(1) * t322 - g(2) * t323 - g(3) * t290 + (-t54 + t338) * t33 + t315 * t19 + t314 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t213 * t113, -t110 + t335, -t35 + t310, -t36 - t309, t63, t29 * t131 + t213 * t45 + t198 + t206, t28 * t131 + t45 * t113 - g(1) * (-t134 * t186 - t190 * t95) - g(2) * (-t132 * t186 - t190 * t93) - g(3) * (-t125 * t190 + t266) + t224, pkin(5) * t35 - t113 * t326, t326 * t19 + (t213 * t33 + t1 + t198) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110 - t335, t19 * t113 - t14 * t213 - t216 + t5;];
tau_reg  = t2;
