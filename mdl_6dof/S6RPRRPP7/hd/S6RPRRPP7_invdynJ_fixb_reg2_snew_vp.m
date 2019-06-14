% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRRPP7
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 21:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRRPP7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP7_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:48:31
% EndTime: 2019-05-05 21:48:41
% DurationCPUTime: 3.74s
% Computational Cost: add. (6288->315), mult. (12199->334), div. (0->0), fcn. (7410->6), ass. (0->212)
t180 = sin(qJ(3));
t183 = cos(qJ(3));
t179 = sin(qJ(4));
t182 = cos(qJ(4));
t238 = qJD(1) * t183;
t155 = qJD(3) * t179 + t182 * t238;
t151 = t155 ^ 2;
t153 = -t182 * qJD(3) + t179 * t238;
t282 = t153 ^ 2;
t87 = -t282 - t151;
t325 = t183 * t87;
t172 = t183 * qJDD(1);
t234 = qJD(1) * qJD(3);
t227 = t180 * t234;
t158 = t172 - t227;
t210 = -t179 * qJDD(3) - t182 * t158;
t102 = -qJD(4) * t153 - t210;
t169 = qJD(1) * t180 + qJD(4);
t251 = t153 * t169;
t294 = t102 + t251;
t272 = t179 * t294;
t222 = -t182 * qJDD(3) + t179 * t158;
t71 = (qJD(4) - t169) * t155 + t222;
t42 = t182 * t71 - t272;
t26 = t180 * t42 + t325;
t278 = pkin(7) + pkin(1);
t349 = t278 * t26;
t295 = t102 - t251;
t281 = t169 ^ 2;
t105 = -t281 - t151;
t116 = t155 * t153;
t226 = t183 * t234;
t232 = t180 * qJDD(1);
t157 = -t226 - t232;
t152 = qJDD(4) - t157;
t292 = t116 + t152;
t303 = t182 * t292;
t334 = t105 * t179 + t303;
t342 = t180 * t334;
t339 = t183 * t295 + t342;
t308 = t179 * t292;
t53 = t105 * t182 - t308;
t343 = qJ(2) * t53;
t348 = -t278 * t339 - t343;
t346 = pkin(3) * t53;
t345 = pkin(8) * t53;
t344 = pkin(8) * t334;
t273 = t179 * t295;
t288 = -t282 + t151;
t101 = qJD(4) * t155 + t222;
t249 = t169 * t155;
t296 = t101 + t249;
t331 = -t180 * t288 + t183 * (t182 * t296 + t273);
t329 = pkin(3) * t87;
t341 = -pkin(8) * t42 - t329;
t340 = pkin(3) * t295 + t344;
t287 = t282 - t281;
t297 = t101 - t249;
t338 = t180 * t297 + t183 * (-t182 * t287 + t308);
t126 = t151 - t281;
t293 = -t116 + t152;
t302 = t182 * t293;
t315 = t180 * t294 + t183 * (t126 * t179 + t302);
t289 = -t281 - t282;
t320 = t179 * t289 + t302;
t307 = t179 * t293;
t318 = t182 * t289 - t307;
t332 = t180 * t318 - t183 * t296;
t330 = qJ(2) * t320 - t278 * t332;
t328 = pkin(3) * t320;
t327 = pkin(8) * t320;
t336 = t126 * t182 - t307;
t264 = t182 * t295;
t322 = -t179 * t296 + t264;
t319 = t179 * t287 + t303;
t317 = -pkin(3) * t296 + pkin(8) * t318;
t326 = qJ(5) * t87;
t64 = t182 * t294;
t39 = t179 * t71 + t64;
t313 = qJ(5) * t292;
t312 = qJ(5) * t295;
t311 = qJ(6) * t294;
t299 = pkin(4) * t293 + qJ(5) * t289;
t186 = qJD(1) ^ 2;
t298 = t278 * t186;
t237 = qJD(5) * t169;
t163 = 0.2e1 * t237;
t236 = qJD(6) * t153;
t291 = 0.2e1 * t236 + t163;
t164 = -0.2e1 * t237;
t290 = -0.2e1 * t236 + t164;
t108 = t180 * t116;
t247 = t169 * t182;
t230 = t153 * t247;
t198 = t183 * (t101 * t179 + t230) - t108;
t233 = qJD(2) * qJD(1);
t174 = 0.2e1 * t233;
t176 = qJDD(1) * qJ(2);
t181 = sin(qJ(1));
t184 = cos(qJ(1));
t215 = t184 * g(1) + t181 * g(2);
t208 = -t176 + t215;
t205 = t174 - t208;
t212 = -t158 + t227;
t213 = -t157 + t226;
t62 = t213 * pkin(3) + t212 * pkin(8) + t205 - t298;
t224 = t181 * g(1) - t184 * g(2);
t214 = qJDD(2) - t224;
t240 = t186 * qJ(2);
t202 = t214 - t240;
t131 = -t278 * qJDD(1) + t202;
t110 = t183 * g(3) - t180 * t131;
t185 = qJD(3) ^ 2;
t218 = pkin(3) * t180 - pkin(8) * t183;
t206 = t186 * t218;
t82 = -t185 * pkin(3) + qJDD(3) * pkin(8) - t180 * t206 - t110;
t37 = t179 * t82 - t182 * t62;
t204 = -t152 * pkin(4) - qJ(5) * t281 + qJDD(5) + t37;
t196 = -t152 * pkin(5) + t204 - t311;
t107 = pkin(4) * t153 - qJ(5) * t155;
t225 = -pkin(5) * t153 - t107;
t207 = (-0.2e1 * qJD(6) - t225) * t155;
t14 = t207 + t196;
t120 = -pkin(5) * t169 - qJ(6) * t155;
t38 = t179 * t62 + t182 * t82;
t219 = pkin(4) * t281 - t152 * qJ(5) + t153 * t107 - t38;
t201 = -pkin(5) * t282 + t101 * qJ(6) + t169 * t120 - t219;
t15 = t201 + t291;
t279 = pkin(4) + pkin(5);
t286 = qJ(5) * t15 - t279 * t14;
t285 = qJ(5) * t71 + t279 * t294;
t248 = t169 * t179;
t122 = t155 * t248;
t259 = t183 * (t182 * t102 - t122) + t108;
t284 = -t105 * t279 + t201 + t313;
t280 = 0.2e1 * qJD(5);
t283 = -t101 * pkin(5) - t282 * qJ(6) + qJDD(6) + (t280 + t120) * t155;
t277 = t101 * pkin(4);
t109 = t180 * g(3) + t183 * t131;
t81 = qJDD(3) * pkin(3) + t185 * pkin(8) - t183 * t206 + t109;
t271 = t179 * t81;
t263 = t182 * t81;
t257 = qJ(5) * t179;
t256 = qJ(5) * t182;
t255 = qJDD(1) * pkin(1);
t250 = t155 * t107;
t177 = t180 ^ 2;
t246 = t177 * t186;
t178 = t183 ^ 2;
t245 = t178 * t186;
t243 = t180 * t152;
t229 = t180 * t186 * t183;
t242 = t180 * (qJDD(3) + t229);
t241 = t183 * (qJDD(3) - t229);
t239 = t177 + t178;
t231 = t155 * t280;
t228 = (-t101 - t296) * pkin(4);
t18 = t179 * t37 + t182 * t38;
t220 = -t182 * t101 + t153 * t248;
t63 = t179 * t102 + t155 * t247;
t23 = t163 - t219;
t24 = t204 + t250;
t217 = -pkin(4) * t24 + qJ(5) * t23;
t216 = -pkin(4) * t294 - qJ(5) * t297;
t211 = t179 * t38 - t182 * t37;
t61 = t183 * t109 - t180 * t110;
t209 = qJ(2) + t218;
t203 = (-t153 * t179 - t155 * t182) * t169;
t199 = t183 * (t122 - t230) + t243;
t197 = -pkin(4) * t105 - t219 + t313;
t193 = -pkin(4) * t249 + t81;
t192 = -t24 + t299;
t191 = -t196 + t299;
t190 = t193 + t312;
t189 = t190 + t231;
t188 = t193 - t277 + 0.2e1 * t312;
t187 = t190 + t283;
t160 = t239 * qJDD(1);
t159 = t172 - 0.2e1 * t227;
t156 = 0.2e1 * t226 + t232;
t141 = 0.2e1 * qJD(6) * t155;
t137 = -t202 + t255;
t123 = t208 - 0.2e1 * t233 + t298;
t119 = -t242 + t183 * (-t185 - t245);
t118 = t180 * (-t185 - t246) + t241;
t78 = (qJD(4) + t169) * t153 + t210;
t45 = -qJ(5) * t296 + qJ(6) * t293;
t43 = -t182 * t297 + t272;
t40 = -t179 * t297 - t64;
t34 = t183 * t78 - t342;
t29 = -qJ(6) * t292 + t279 * t295;
t27 = t180 * t43 - t325;
t25 = t189 - t277;
t22 = t24 - t326;
t21 = -pkin(4) * t87 + t23;
t20 = t228 + t189;
t19 = t188 + t231;
t16 = t187 - t277;
t13 = t18 * t180 + t183 * t81;
t12 = -qJ(6) * t105 + t188 + t283;
t11 = t179 * t24 + t182 * t23;
t10 = t179 * t23 - t182 * t24;
t9 = t225 * t155 + t141 - t196 + t311 + t326;
t8 = -pkin(5) * t296 - qJ(6) * t289 + t187 + t228;
t7 = -qJ(6) * t71 + t279 * t87 - t201 + t290;
t6 = qJ(5) * t16 - qJ(6) * t14;
t5 = t11 * t180 + t183 * t25;
t4 = t14 * t179 + t15 * t182;
t3 = -t14 * t182 + t15 * t179;
t2 = -qJ(6) * t15 + t279 * t16;
t1 = t16 * t183 + t180 * t4;
t17 = [0, 0, 0, 0, 0, qJDD(1), t224, t215, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t214 - 0.2e1 * t255, t174 + 0.2e1 * t176 - t215, pkin(1) * t137 + qJ(2) * (-t186 * pkin(1) + t205), -t212 * t183, -t156 * t183 - t159 * t180, t241 - t180 * (t185 - t245), t213 * t180, t183 * (-t185 + t246) - t242, 0, qJ(2) * t156 - t278 * t118 - t180 * t123, qJ(2) * t159 - t278 * t119 - t183 * t123, t278 * t160 - t239 * t240 - t61, -qJ(2) * t123 - t278 * t61, t259, -t331, t315, t198, -t338, t199, t183 * (-t271 - t327) - t180 * (t37 - t328) + t330, t183 * (-t263 - t345) - t180 * (t38 - t346) + t343 - t278 * t34, -t183 * t211 - t209 * t39 + t349, -t278 * t13 + t209 * t211, t259, t315, t331, t199, t338, t198, t183 * (-t179 * t20 - t256 * t296 - t327) - t180 * (-t192 - t328) + t330, t183 * (-pkin(8) * t40 - t179 * t21 + t182 * t22) - t180 * (-pkin(3) * t40 - t216) + qJ(2) * t40 - t278 * t27, t183 * (-pkin(4) * t273 + t182 * t19 + t345) - t180 * (t164 - t197 + t346) + t348, t183 * (-pkin(8) * t10 + (-pkin(4) * t179 + t256) * t25) - t180 * (-pkin(3) * t10 - t217) + qJ(2) * t10 - t278 * t5, t259, t331, -t315, t198, -t338, t243 + t183 * (-t153 * t182 + t155 * t179) * t169, t183 * (-t179 * t8 + t182 * t45 - t327) + (pkin(5) * t293 + t191 - t207 + t328) * t180 + t330, t183 * (t12 * t182 - t179 * t29 + t345) - t180 * (-t284 + t290 + t346) + t348, t183 * (-pkin(8) * t39 - t179 * t7 + t182 * t9) - t180 * (-pkin(3) * t39 - t285) + qJ(2) * t39 - t349, t183 * (-pkin(8) * t3 - t179 * t2 + t182 * t6) - t180 * (-pkin(3) * t3 - t286) + qJ(2) * t3 - t278 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t186, -t137, 0, 0, 0, 0, 0, 0, t118, t119, -t160, t61, 0, 0, 0, 0, 0, 0, t332, t34, -t26, t13, 0, 0, 0, 0, 0, 0, t332, t27, t339, t5, 0, 0, 0, 0, 0, 0, t332, t339, t26, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t229, (-t177 + t178) * t186, t172, -t229, -t232, qJDD(3), t109, t110, 0, 0, t63, t322, -t336, t220, t319, t203, t263 + t317, pkin(3) * t78 - t271 - t344, t18 + t341, pkin(3) * t81 + pkin(8) * t18, t63, -t336, -t322, t203, -t319, t220, t182 * t20 - t257 * t296 + t317, pkin(8) * t43 + t179 * t22 + t182 * t21 - t329, pkin(4) * t264 + t179 * t19 + t340, pkin(8) * t11 + (pkin(4) * t182 + pkin(3) + t257) * t25, t63, -t322, t336, t220, t319, t203, t179 * t45 + t182 * t8 + t317, t12 * t179 + t182 * t29 + t340, t179 * t9 + t182 * t7 - t341, pkin(3) * t16 + pkin(8) * t4 + t179 * t6 + t182 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116, t288, t294, -t116, -t297, t152, -t37, -t38, 0, 0, t116, t294, -t288, t152, t297, -t116, t192, t216, t163 + t197, t217, t116, -t288, -t294, -t116, -t297, t152, -t250 + t141 + (t293 - t116) * pkin(5) + t191, t284 + t291, t285, t286; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t293, t294, t105, t24, 0, 0, 0, 0, 0, 0, -t293, t105, -t294, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t296, t295, t87, t16;];
tauJ_reg  = t17;
