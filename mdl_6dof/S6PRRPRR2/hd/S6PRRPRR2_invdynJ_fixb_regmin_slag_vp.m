% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRPRR2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:49
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRPRR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:46:37
% EndTime: 2021-01-16 03:47:00
% DurationCPUTime: 6.67s
% Computational Cost: add. (5321->509), mult. (12629->722), div. (0->0), fcn. (10376->18), ass. (0->260)
t204 = cos(qJ(3));
t318 = cos(pkin(12));
t255 = t318 * t204;
t172 = qJD(2) * t255;
t193 = sin(pkin(12));
t200 = sin(qJ(3));
t291 = qJD(2) * t200;
t148 = t193 * t291 - t172;
t346 = qJD(5) + qJD(6);
t357 = t148 + t346;
t256 = t318 * t200;
t159 = t193 * t204 + t256;
t150 = t159 * qJD(3);
t305 = t193 * t200;
t227 = t255 - t305;
t153 = t227 * qJD(3);
t195 = sin(pkin(6));
t201 = sin(qJ(2));
t292 = qJD(1) * t201;
t273 = t195 * t292;
t332 = qJD(3) * pkin(3);
t278 = t200 * t332;
t356 = -pkin(4) * t150 + pkin(9) * t153 + t273 - t278;
t205 = cos(qJ(2));
t299 = t195 * t205;
t272 = qJD(1) * t299;
t124 = t227 * t272;
t197 = qJ(4) + pkin(8);
t262 = qJD(3) * t197;
t143 = qJD(4) * t204 - t200 * t262;
t223 = -qJD(4) * t200 - t204 * t262;
t81 = t143 * t318 + t193 * t223;
t320 = t81 - t124;
t151 = t159 * qJD(2);
t199 = sin(qJ(5));
t203 = cos(qJ(5));
t285 = t203 * qJD(3);
t127 = t151 * t199 - t285;
t129 = qJD(3) * t199 + t151 * t203;
t198 = sin(qJ(6));
t202 = cos(qJ(6));
t239 = t127 * t198 - t202 * t129;
t62 = t202 * t127 + t129 * t198;
t355 = t239 * t62;
t162 = t198 * t203 + t199 * t202;
t323 = t357 * t162;
t289 = qJD(5) * t199;
t312 = t148 * t199;
t354 = t289 + t312;
t353 = t239 ^ 2 - t62 ^ 2;
t140 = qJD(5) + t148;
t138 = qJD(6) + t140;
t286 = qJD(6) * t202;
t287 = qJD(6) * t198;
t283 = qJD(2) * qJD(3);
t267 = t200 * t283;
t216 = qJDD(2) * t159 - t193 * t267;
t104 = qJD(3) * t172 + t216;
t47 = qJD(5) * t285 + t199 * qJDD(3) + t203 * t104 - t151 * t289;
t48 = qJD(5) * t129 - t203 * qJDD(3) + t104 * t199;
t8 = -t127 * t286 - t129 * t287 - t198 * t48 + t202 * t47;
t352 = t138 * t62 + t8;
t196 = cos(pkin(6));
t319 = cos(pkin(11));
t259 = t319 * t201;
t194 = sin(pkin(11));
t302 = t194 * t205;
t146 = t196 * t259 + t302;
t189 = qJ(3) + pkin(12);
t182 = sin(t189);
t183 = cos(t189);
t260 = t195 * t319;
t106 = t146 * t183 - t182 * t260;
t258 = t319 * t205;
t303 = t194 * t201;
t144 = t196 * t303 - t258;
t304 = t194 * t195;
t108 = -t144 * t183 + t182 * t304;
t301 = t195 * t201;
t136 = t182 * t196 + t183 * t301;
t145 = -t196 * t258 + t303;
t147 = t196 * t302 + t259;
t251 = t197 * qJD(2) + t273;
t293 = qJD(1) * t196;
t119 = -t200 * t251 + t204 * t293;
t113 = t119 + t332;
t120 = t200 * t293 + t204 * t251;
t257 = t318 * t120;
t51 = t193 * t113 + t257;
t46 = qJD(3) * pkin(9) + t51;
t181 = pkin(3) * t204 + pkin(2);
t139 = -qJD(2) * t181 + qJD(4) - t272;
t67 = pkin(4) * t148 - pkin(9) * t151 + t139;
t25 = t199 * t67 + t203 * t46;
t19 = -pkin(10) * t127 + t25;
t17 = t19 * t287;
t192 = qJ(5) + qJ(6);
t187 = sin(t192);
t188 = cos(t192);
t109 = t193 * t120;
t50 = t113 * t318 - t109;
t45 = -qJD(3) * pkin(4) - t50;
t32 = t127 * pkin(5) + t45;
t351 = t32 * t62 - g(1) * (-t108 * t188 - t147 * t187) - g(2) * (-t106 * t188 - t145 * t187) - g(3) * (-t136 * t188 + t187 * t299) + t17;
t282 = t196 * qJDD(1);
t171 = t204 * t282;
t284 = qJD(1) * qJD(2);
t131 = qJDD(2) * pkin(8) + (qJDD(1) * t201 + t205 * t284) * t195;
t213 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + qJD(3) * t293 + t131;
t238 = t251 * qJD(3);
t40 = qJDD(3) * pkin(3) - t200 * t213 - t204 * t238 + t171;
t41 = (-t238 + t282) * t200 + t213 * t204;
t16 = t193 * t40 + t318 * t41;
t14 = qJDD(3) * pkin(9) + t16;
t268 = t201 * t284;
t169 = t195 * t268;
t218 = pkin(3) * t267 - qJDD(2) * t181 + qJDD(4) + t169;
t265 = qJDD(1) * t299;
t102 = t218 - t265;
t281 = t200 * qJDD(2);
t243 = -qJDD(2) * t255 + t193 * t281;
t103 = qJD(2) * t150 + t243;
t30 = pkin(4) * t103 - pkin(9) * t104 + t102;
t29 = t203 * t30;
t215 = -qJD(5) * t25 - t199 * t14 + t29;
t99 = qJDD(5) + t103;
t2 = pkin(5) * t99 - pkin(10) * t47 + t215;
t288 = qJD(5) * t203;
t232 = -t203 * t14 - t199 * t30 - t67 * t288 + t289 * t46;
t3 = -pkin(10) * t48 - t232;
t275 = -t198 * t3 + t202 * t2;
t24 = -t199 * t46 + t203 * t67;
t18 = -pkin(10) * t129 + t24;
t11 = pkin(5) * t140 + t18;
t328 = t19 * t202;
t5 = t11 * t198 + t328;
t350 = t32 * t239 - g(1) * (-t108 * t187 + t147 * t188) - g(2) * (-t106 * t187 + t145 * t188) - g(3) * (-t136 * t187 - t188 * t299) - qJD(6) * t5 + t275;
t212 = qJD(6) * t239 - t198 * t47 - t202 * t48;
t349 = -t138 * t239 + t212;
t348 = t124 * t199 - t356 * t203;
t101 = -pkin(4) * t227 - pkin(9) * t159 - t181;
t166 = t197 * t204;
t126 = t166 * t318 - t197 * t305;
t347 = -t101 * t288 + t126 * t289 + t356 * t199 - t320 * t203;
t321 = t143 * t193 - t159 * t272 - t318 * t223;
t90 = t162 * t159;
t296 = t203 * t153;
t230 = -t159 * t289 + t296;
t161 = t198 * t199 - t202 * t203;
t324 = t357 * t161;
t95 = qJDD(6) + t99;
t345 = t138 * t324 - t162 * t95;
t206 = qJD(3) ^ 2;
t246 = g(1) * t147 + g(2) * t145;
t344 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t206 + t195 * (-g(3) * t205 + t268) - t169 + t246 + t265;
t264 = qJD(6) * t11 + t3;
t343 = t198 * t2 + t202 * t264;
t114 = t203 * t126;
t342 = -pkin(10) * t296 + pkin(5) * t150 - t199 * t81 + (-t114 + (pkin(10) * t159 - t101) * t199) * qJD(5) + t348;
t341 = pkin(3) * t193;
t340 = pkin(3) * t200;
t154 = t196 * t204 - t200 * t301;
t339 = g(3) * t154;
t338 = g(3) * t195;
t177 = pkin(9) + t341;
t336 = pkin(10) + t177;
t298 = t199 * t153;
t231 = t159 * t288 + t298;
t335 = pkin(10) * t231 + t347;
t15 = -t193 * t41 + t318 * t40;
t55 = t119 * t318 - t109;
t279 = pkin(3) * t291;
t82 = pkin(4) * t151 + pkin(9) * t148 + t279;
t334 = t199 * t82 + t203 * t55;
t333 = qJD(2) * pkin(2);
t331 = t151 * t62;
t330 = t151 * t239;
t327 = t199 * t47;
t326 = t199 * t99;
t325 = t199 * t101 + t114;
t322 = pkin(5) * t231 + t321;
t316 = t127 * t140;
t315 = t127 * t151;
t314 = t129 * t140;
t313 = t129 * t151;
t311 = t159 * t199;
t310 = t159 * t203;
t309 = t183 * t187;
t308 = t183 * t188;
t307 = t183 * t199;
t306 = t183 * t205;
t300 = t195 * t204;
t297 = t199 * t205;
t295 = qJDD(1) - g(3);
t190 = t200 ^ 2;
t294 = -t204 ^ 2 + t190;
t290 = qJD(2) * t201;
t280 = t204 * qJDD(2);
t277 = t195 * t297;
t276 = t203 * t299;
t274 = t318 * pkin(3);
t271 = t195 * t290;
t270 = qJD(2) * t299;
t266 = t204 * t283;
t261 = qJD(5) * t336;
t254 = t195 * t295;
t53 = t119 * t193 + t257;
t125 = t166 * t193 + t197 * t256;
t252 = t140 * t203;
t250 = t200 * t270;
t249 = -t323 * t138 - t161 * t95;
t248 = t354 * pkin(5) - t53;
t13 = -qJDD(3) * pkin(4) - t15;
t178 = -t274 - pkin(4);
t247 = g(1) * t144 - g(2) * t146;
t156 = t336 * t199;
t245 = pkin(10) * t312 + qJD(6) * t156 + t199 * t261 + t334;
t157 = t336 * t203;
t75 = t203 * t82;
t244 = pkin(5) * t151 + qJD(6) * t157 - t199 * t55 + t75 + (pkin(10) * t148 + t261) * t203;
t93 = t203 * t101;
t31 = -pkin(5) * t227 - pkin(10) * t310 - t126 * t199 + t93;
t33 = -pkin(10) * t311 + t325;
t242 = t198 * t31 + t202 * t33;
t155 = t196 * t200 + t201 * t300;
t88 = t193 * t154 + t155 * t318;
t234 = -t203 * t88 + t277;
t68 = -t199 * t88 - t276;
t241 = t198 * t234 + t202 * t68;
t240 = t198 * t68 - t202 * t234;
t207 = qJD(2) ^ 2;
t237 = qJDD(2) * t205 - t201 * t207;
t236 = -t354 * t140 + t203 * t99;
t235 = -g(1) * t194 + g(2) * t319;
t228 = t140 * t45 - t177 * t99;
t226 = g(1) * (t144 * t182 + t183 * t304) + g(2) * (-t146 * t182 - t183 * t260) + g(3) * (-t182 * t301 + t183 * t196);
t225 = t237 * t195;
t222 = t155 * qJD(3);
t221 = -g(3) * t301 + t247;
t220 = g(3) * t299 - t246;
t164 = -t272 - t333;
t219 = -qJD(2) * t164 - t131 - t247;
t217 = t220 * t183;
t211 = qJD(5) * t140 * t177 + t13 + t226;
t210 = -pkin(8) * qJDD(3) + (t164 + t272 - t333) * qJD(3);
t209 = -t222 - t250;
t165 = -t203 * pkin(5) + t178;
t118 = qJD(3) * t154 + t204 * t270;
t91 = t161 * t159;
t87 = -t154 * t318 + t155 * t193;
t79 = pkin(5) * t311 + t125;
t54 = t118 * t318 + t193 * t209;
t52 = t118 * t193 - t209 * t318;
t27 = -t287 * t311 + (t346 * t310 + t298) * t202 + t230 * t198;
t26 = -t161 * t153 - t346 * t90;
t22 = qJD(5) * t234 - t199 * t54 + t203 * t271;
t21 = qJD(5) * t68 + t199 * t271 + t203 * t54;
t6 = pkin(5) * t48 + t13;
t4 = t11 * t202 - t19 * t198;
t1 = [t295, 0, t225, (-qJDD(2) * t201 - t205 * t207) * t195, 0, 0, 0, 0, 0, t154 * qJDD(3) + t204 * t225 + (-t222 - 0.2e1 * t250) * qJD(3), -qJD(3) * t118 - qJDD(3) * t155 + (-t200 * t237 - t205 * t266) * t195, -qJD(3) * t52 - qJDD(3) * t87 + (-t103 * t205 + t148 * t290) * t195, -qJD(3) * t54 - qJDD(3) * t88 + (-t104 * t205 + t151 * t290) * t195, -t103 * t88 + t104 * t87 - t148 * t54 + t151 * t52, -t15 * t87 + t16 * t88 - t50 * t52 + t51 * t54 - g(3) + (-t102 * t205 + t139 * t290) * t195, 0, 0, 0, 0, 0, t127 * t52 + t140 * t22 + t48 * t87 + t68 * t99, t129 * t52 - t140 * t21 + t234 * t99 + t47 * t87, 0, 0, 0, 0, 0, (-qJD(6) * t240 - t198 * t21 + t202 * t22) * t138 + t241 * t95 + t52 * t62 - t87 * t212, -(qJD(6) * t241 + t198 * t22 + t202 * t21) * t138 - t240 * t95 - t52 * t239 + t87 * t8; 0, qJDD(2), t295 * t299 + t246, -t201 * t254 - t247, qJDD(2) * t190 + 0.2e1 * t200 * t266, 0.2e1 * t200 * t280 - 0.2e1 * t283 * t294, qJDD(3) * t200 + t204 * t206, qJDD(3) * t204 - t200 * t206, 0, t210 * t200 + t344 * t204, -t344 * t200 + t210 * t204, -t148 * t273 - qJDD(3) * t125 - t102 * t227 - t103 * t181 + t139 * t150 - t217 + (t148 * t340 - t321) * qJD(3), -t151 * t273 - qJDD(3) * t126 + t102 * t159 - t104 * t181 + t139 * t153 + t220 * t182 + (t151 * t340 - t320) * qJD(3), -t103 * t126 + t104 * t125 - t148 * t320 - t15 * t159 - t150 * t51 + t151 * t321 - t153 * t50 + t16 * t227 + t221, t16 * t126 - t15 * t125 - t102 * t181 + t139 * t278 - g(1) * (-t144 * t197 - t147 * t181) - g(2) * (-t145 * t181 + t146 * t197) + t320 * t51 - t321 * t50 + (-t139 * t292 - g(3) * (t181 * t205 + t197 * t201)) * t195, t129 * t230 + t310 * t47, (-t127 * t203 - t129 * t199) * t153 + (-t327 - t203 * t48 + (t127 * t199 - t129 * t203) * qJD(5)) * t159, t129 * t150 + t140 * t230 - t227 * t47 + t310 * t99, -t127 * t150 - t140 * t231 + t227 * t48 - t311 * t99, t140 * t150 - t227 * t99, t125 * t48 + t24 * t150 - t29 * t227 + t93 * t99 + t348 * t140 + t321 * t127 + (-t217 + (-t126 * t140 + t159 * t45 + t227 * t46) * qJD(5)) * t203 + ((-qJD(5) * t101 - t81) * t140 - t126 * t99 - (-qJD(5) * t67 - t14) * t227 + t13 * t159 + t45 * t153 + t221) * t199, -t325 * t99 - t232 * t227 - t25 * t150 + t125 * t47 + t45 * t296 - g(1) * (-t144 * t203 + t147 * t307) - g(2) * (t145 * t307 + t146 * t203) - (-t183 * t297 + t201 * t203) * t338 + (t13 * t203 - t289 * t45) * t159 + t347 * t140 + t321 * t129, -t239 * t26 - t8 * t91, -t212 * t91 + t239 * t27 - t26 * t62 - t8 * t90, t138 * t26 - t150 * t239 - t227 * t8 - t91 * t95, -t138 * t27 - t150 * t62 - t212 * t227 - t90 * t95, t138 * t150 - t227 * t95, (-t198 * t33 + t202 * t31) * t95 - t275 * t227 + t4 * t150 - t79 * t212 + t6 * t90 + t32 * t27 - g(1) * (-t144 * t187 - t147 * t308) - g(2) * (-t145 * t308 + t146 * t187) + t322 * t62 - (t187 * t201 + t188 * t306) * t338 + (t335 * t198 + t342 * t202) * t138 + (-t138 * t242 + t227 * t5) * qJD(6), -t242 * t95 + (-t17 + t343) * t227 - t5 * t150 + t79 * t8 - t6 * t91 + t32 * t26 - g(1) * (-t144 * t188 + t147 * t309) - g(2) * (t145 * t309 + t146 * t188) - t322 * t239 - (-t187 * t306 + t188 * t201) * t338 + ((-qJD(6) * t31 + t335) * t202 + (qJD(6) * t33 - t342) * t198) * t138; 0, 0, 0, 0, -t200 * t207 * t204, t294 * t207, t281, t280, qJDD(3), t200 * t219 + t235 * t300 + t171 - t339, g(3) * t155 + (-t195 * t235 - t282) * t200 + t219 * t204, t53 * qJD(3) - t139 * t151 + (qJDD(3) * t318 - t148 * t291) * pkin(3) - t226 + t15, t55 * qJD(3) + t139 * t148 + g(1) * t108 + g(2) * t106 + g(3) * t136 + (-qJDD(3) * t193 - t151 * t291) * pkin(3) - t16, (t51 - t53) * t151 + (-t50 + t55) * t148 + (-t103 * t193 - t104 * t318) * pkin(3), -t139 * t279 + t15 * t274 + t16 * t341 + t50 * t53 - t51 * t55 + (-g(1) * (t144 * t200 + t194 * t300) - g(2) * (-t146 * t200 - t204 * t260) - t339) * pkin(3), t129 * t252 + t327, (t47 - t316) * t203 + (-t48 - t314) * t199, t140 * t252 - t313 + t326, t236 + t315, -t140 * t151, -t53 * t127 - t75 * t140 - t24 * t151 + t178 * t48 + (t55 * t140 + t228) * t199 - t211 * t203, -t53 * t129 + t140 * t334 + t25 * t151 + t178 * t47 + t199 * t211 + t203 * t228, t162 * t8 + t239 * t324, -t161 * t8 + t162 * t212 + t239 * t323 + t324 * t62, t330 - t345, t249 + t331, -t138 * t151, (-t156 * t202 - t157 * t198) * t95 - t165 * t212 + t6 * t161 - t4 * t151 + t248 * t62 + t323 * t32 + (t198 * t245 - t202 * t244) * t138 - t226 * t188, -(-t156 * t198 + t157 * t202) * t95 + t165 * t8 + t6 * t162 + t5 * t151 - t248 * t239 - t324 * t32 + (t198 * t244 + t202 * t245) * t138 + t226 * t187; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t151 * qJD(3) + t243, (t172 - t148) * qJD(3) + t216, -t148 ^ 2 - t151 ^ 2, t148 * t51 + t151 * t50 - t205 * t254 + t218 - t246, 0, 0, 0, 0, 0, t236 - t315, -t140 ^ 2 * t203 - t313 - t326, 0, 0, 0, 0, 0, t249 - t331, t330 + t345; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129 * t127, -t127 ^ 2 + t129 ^ 2, t47 + t316, t314 - t48, t99, t25 * t140 - t45 * t129 - g(1) * (-t108 * t199 + t147 * t203) - g(2) * (-t106 * t199 + t145 * t203) - g(3) * (-t136 * t199 - t276) + t215, t24 * t140 + t45 * t127 - g(1) * (-t108 * t203 - t147 * t199) - g(2) * (-t106 * t203 - t145 * t199) - g(3) * (-t136 * t203 + t277) + t232, -t355, t353, t352, t349, t95, -(-t18 * t198 - t328) * t138 + (-t129 * t62 - t138 * t287 + t202 * t95) * pkin(5) + t350, (-t138 * t19 - t2) * t198 + (t138 * t18 - t264) * t202 + (t129 * t239 - t138 * t286 - t198 * t95) * pkin(5) + t351; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t355, t353, t352, t349, t95, t138 * t5 + t350, t138 * t4 - t343 + t351;];
tau_reg = t1;
