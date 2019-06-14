% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
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
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 17:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRPRP2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP2_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:34:25
% EndTime: 2019-05-05 17:34:39
% DurationCPUTime: 6.65s
% Computational Cost: add. (15384->427), mult. (33575->567), div. (0->0), fcn. (22786->10), ass. (0->260)
t220 = sin(pkin(10));
t222 = cos(pkin(10));
t227 = sin(qJ(3));
t229 = cos(qJ(3));
t279 = qJD(1) * t229;
t280 = qJD(1) * t227;
t188 = t220 * t280 - t222 * t279;
t185 = qJD(5) + t188;
t317 = t185 ^ 2;
t190 = t220 * t279 + t222 * t280;
t226 = sin(qJ(5));
t228 = cos(qJ(5));
t172 = -t228 * qJD(3) + t190 * t226;
t318 = t172 ^ 2;
t147 = t318 - t317;
t271 = qJD(1) * qJD(3);
t262 = t229 * t271;
t270 = t227 * qJDD(1);
t196 = t262 + t270;
t212 = t229 * qJDD(1);
t263 = t227 * t271;
t249 = t212 - t263;
t257 = t196 * t220 - t222 * t249;
t164 = qJDD(5) + t257;
t174 = qJD(3) * t226 + t190 * t228;
t290 = t174 * t172;
t112 = -t290 - t164;
t301 = t112 * t226;
t82 = -t147 * t228 - t301;
t151 = t185 * t174;
t167 = t222 * t196 + t220 * t249;
t259 = t228 * qJDD(3) - t167 * t226;
t244 = qJD(5) * t174 - t259;
t92 = -t151 + t244;
t382 = t227 * (t220 * t92 + t222 * t82) + t229 * (t220 * t82 - t222 * t92);
t171 = t174 ^ 2;
t330 = -t171 - t317;
t76 = t228 * t330 + t301;
t381 = pkin(2) * t76;
t380 = pkin(3) * t76;
t379 = pkin(4) * t76;
t378 = pkin(8) * t76;
t300 = t112 * t228;
t78 = -t226 * t330 + t300;
t377 = pkin(8) * t78;
t376 = t220 * t78;
t375 = t222 * t78;
t223 = cos(pkin(9));
t374 = t223 * t76;
t329 = t171 - t318;
t247 = -qJDD(3) * t226 - t167 * t228;
t239 = -qJD(5) * t172 - t247;
t291 = t172 * t185;
t326 = -t291 + t239;
t307 = t226 * t326;
t331 = t151 + t244;
t60 = t228 * t331 + t307;
t371 = t227 * (-t220 * t329 + t222 * t60) + t229 * (t220 * t60 + t222 * t329);
t221 = sin(pkin(9));
t327 = -t290 + t164;
t298 = t327 * t228;
t324 = -t317 - t318;
t334 = t226 * t324 + t298;
t299 = t327 * t226;
t333 = t228 * t324 - t299;
t348 = t220 * t331 + t222 * t333;
t349 = t220 * t333 - t222 * t331;
t363 = -t227 * t349 + t229 * t348;
t368 = pkin(1) * (t221 * t363 - t223 * t334) + pkin(7) * t363 - pkin(2) * t334;
t367 = qJ(4) * t349;
t366 = -t147 * t226 + t300;
t365 = pkin(3) * t349 + pkin(8) * t333;
t364 = -pkin(3) * t334 + qJ(4) * t348;
t362 = t227 * t348 + t229 * t349;
t325 = t291 + t239;
t148 = -t171 + t317;
t350 = -t148 * t226 + t298;
t361 = t227 * (t220 * t325 + t222 * t350) + t229 * (t220 * t350 - t222 * t325);
t358 = pkin(4) * t334;
t356 = pkin(8) * t334;
t355 = qJ(6) * t326;
t351 = t228 * t148 + t299;
t328 = t171 + t318;
t347 = pkin(4) * t328;
t165 = t190 * t188;
t323 = qJDD(3) - t165;
t346 = t220 * t323;
t344 = t220 * t328;
t341 = t222 * t323;
t339 = t222 * t328;
t277 = qJD(3) * t190;
t138 = t257 + t277;
t332 = -t226 * t331 + t228 * t326;
t231 = qJD(1) ^ 2;
t315 = sin(qJ(1));
t316 = cos(qJ(1));
t240 = t315 * g(1) - t316 * g(2);
t238 = qJDD(1) * pkin(1) + t240;
t241 = t316 * g(1) + t315 * g(2);
t194 = -t231 * pkin(1) - t241;
t287 = t223 * t194;
t234 = -t231 * pkin(2) + qJDD(1) * pkin(7) + t221 * t238 + t287;
t283 = -g(3) + qJDD(2);
t143 = t227 * t283 + t229 * t234;
t200 = qJD(3) * pkin(3) - qJ(4) * t280;
t217 = t229 ^ 2;
t215 = t217 * t231;
t113 = -pkin(3) * t215 + qJ(4) * t249 - qJD(3) * t200 + t143;
t233 = t227 * t234;
t285 = t227 * t231;
t232 = -t233 - t196 * qJ(4) + qJDD(3) * pkin(3) + (pkin(3) * t285 + qJ(4) * t271 + t283) * t229;
t67 = -0.2e1 * qJD(4) * t188 + t222 * t113 + t220 * t232;
t322 = pkin(5) * t244 - t355;
t130 = pkin(5) * t172 - qJ(6) * t174;
t156 = pkin(4) * t188 - pkin(8) * t190;
t230 = qJD(3) ^ 2;
t57 = -pkin(4) * t230 + qJDD(3) * pkin(8) - t156 * t188 + t67;
t191 = t223 * t238;
t258 = -t221 * t194 + t191;
t154 = -qJDD(1) * pkin(2) - t231 * pkin(7) - t258;
t119 = -t249 * pkin(3) - qJ(4) * t215 + t200 * t280 + qJDD(4) + t154;
t278 = qJD(3) * t188;
t255 = -t167 + t278;
t70 = t138 * pkin(4) + t255 * pkin(8) + t119;
t38 = t226 * t70 + t228 * t57;
t256 = t164 * qJ(6) - t172 * t130 + t38;
t321 = -pkin(5) * (t330 + t317) - qJ(6) * t112 + t256;
t288 = t185 * t228;
t267 = t172 * t288;
t245 = t226 * t244 + t267;
t268 = t222 * t290;
t269 = t220 * t290;
t320 = t227 * (t222 * t245 - t269) + t229 * (t220 * t245 + t268);
t289 = t185 * t226;
t145 = t174 * t289;
t251 = t145 - t267;
t319 = t227 * (t164 * t220 + t222 * t251) + t229 * (-t222 * t164 + t220 * t251);
t186 = t188 ^ 2;
t187 = t190 ^ 2;
t314 = pkin(4) * t220;
t313 = pkin(5) * t228;
t260 = t113 * t220 - t222 * t232;
t246 = -qJDD(3) * pkin(4) - t230 * pkin(8) + t260;
t56 = (0.2e1 * qJD(4) + t156) * t190 + t246;
t310 = t226 * t56;
t308 = t226 * t325;
t273 = qJD(4) * t190;
t66 = t260 + 0.2e1 * t273;
t41 = t220 * t67 - t222 * t66;
t306 = t227 * t41;
t305 = t228 * t56;
t303 = t228 * t325;
t302 = qJ(6) * t228;
t297 = t119 * t220;
t296 = t119 * t222;
t160 = qJDD(3) + t165;
t294 = t160 * t220;
t293 = t160 * t222;
t208 = t229 * t285;
t201 = qJDD(3) + t208;
t286 = t227 * t201;
t202 = qJDD(3) - t208;
t284 = t229 * t202;
t276 = qJD(3) * t220;
t275 = qJD(3) * t222;
t272 = qJD(6) * t185;
t266 = -pkin(1) * t223 - pkin(2);
t265 = t221 * pkin(1) + pkin(7);
t264 = -pkin(4) * t222 - pkin(3);
t261 = -qJ(6) * t226 - pkin(4);
t42 = t220 * t66 + t222 * t67;
t37 = t226 * t57 - t228 * t70;
t17 = t226 * t37 + t228 * t38;
t142 = -t229 * t283 + t233;
t102 = t227 * t142 + t229 * t143;
t177 = 0.2e1 * t272;
t250 = t177 + t256;
t28 = -pkin(5) * t317 + t250;
t29 = -t164 * pkin(5) - qJ(6) * t317 + t130 * t174 + qJDD(6) + t37;
t254 = -pkin(5) * t29 + qJ(6) * t28;
t253 = -pkin(5) * t325 - qJ(6) * t92;
t252 = t172 * t289 - t228 * t244;
t197 = t212 - 0.2e1 * t263;
t16 = t226 * t38 - t228 * t37;
t139 = -t257 + t277;
t242 = (-t172 * t226 - t174 * t228) * t185;
t90 = t228 * t239 - t145;
t237 = t227 * (t222 * t90 + t269) + t229 * (t220 * t90 - t268);
t183 = -0.2e1 * t273;
t236 = 0.2e1 * qJD(6) * t174 - t156 * t190 + t183 - t246 - t322;
t235 = pkin(5) * t327 + qJ(6) * t324 - t29;
t216 = t227 ^ 2;
t213 = t216 * t231;
t207 = -t215 - t230;
t206 = -t213 - t230;
t199 = t213 + t215;
t198 = (t216 + t217) * qJDD(1);
t195 = 0.2e1 * t262 + t270;
t180 = -t187 - t230;
t179 = -t187 + t230;
t178 = t186 - t230;
t176 = -t206 * t227 - t284;
t175 = t207 * t229 - t286;
t157 = -t230 - t186;
t141 = t167 + t278;
t134 = -t186 - t187;
t127 = -t180 * t220 - t293;
t126 = t180 * t222 - t294;
t116 = t157 * t222 - t346;
t115 = t157 * t220 + t341;
t101 = t139 * t222 + t141 * t220;
t100 = t139 * t220 - t141 * t222;
t99 = (qJD(5) + t185) * t172 + t247;
t94 = (-qJD(5) + t185) * t174 + t259;
t89 = t174 * t288 + t226 * t239;
t84 = -t227 * t126 + t229 * t127;
t71 = -t227 * t115 + t229 * t116;
t64 = -t227 * t100 + t229 * t101;
t63 = t228 * t94 + t308;
t61 = -t228 * t92 + t308;
t59 = t226 * t94 - t303;
t58 = -t226 * t92 - t303;
t53 = -t220 * t99 + t375;
t51 = t222 * t99 + t376;
t49 = -t220 * t326 - t375;
t47 = t222 * t326 - t376;
t46 = t222 * t63 - t344;
t45 = t222 * t61 - t344;
t44 = t220 * t63 + t339;
t43 = t220 * t61 + t339;
t40 = t305 - t378;
t39 = t310 - t356;
t35 = -pkin(4) * t58 - t253;
t34 = (pkin(5) * t185 - 0.2e1 * qJD(6)) * t174 + t56 + t322;
t32 = -t227 * t51 + t229 * t53;
t30 = -t227 * t47 + t229 * t49;
t27 = t38 - t379;
t26 = t37 - t358;
t25 = (-t331 - t151) * pkin(5) + t236;
t24 = -pkin(5) * t151 + t236 + t355;
t22 = -t227 * t43 + t229 * t45;
t21 = qJ(6) * t328 + t29;
t20 = (t328 - t317) * pkin(5) + t250;
t19 = t229 * t42 - t306;
t18 = -t235 - t358;
t15 = -t226 * t25 - t302 * t331 - t356;
t14 = -pkin(5) * t307 + t228 * t24 + t378;
t13 = -0.2e1 * t272 - t321 + t379;
t12 = -pkin(8) * t59 - t16;
t11 = t17 * t222 + t220 * t56;
t10 = t17 * t220 - t222 * t56;
t9 = t226 * t29 + t228 * t28;
t8 = t226 * t28 - t228 * t29;
t7 = -pkin(8) * t58 - t20 * t226 + t21 * t228;
t6 = t220 * t34 + t222 * t9;
t5 = t220 * t9 - t222 * t34;
t4 = -pkin(8) * t8 + (pkin(5) * t226 - t302) * t34;
t2 = -pkin(4) * t8 - t254;
t1 = -t227 * t5 + t229 * t6;
t3 = [0, 0, 0, 0, 0, qJDD(1), t240, t241, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t223 - t221 * t231) + t258, -t287 - t221 * t240 + (-0.2e1 * qJDD(1) * t221 - t223 * t231) * pkin(1), 0, pkin(1) * (t221 ^ 2 * t238 + t223 * t191), (t196 + t262) * t227, t195 * t229 + t197 * t227, t286 + t229 * (-t213 + t230), t197 * t229, t227 * (t215 - t230) + t284, 0, -t229 * t154 + pkin(2) * t197 + pkin(7) * t175 + pkin(1) * (t175 * t221 + t197 * t223), t227 * t154 - pkin(2) * t195 + pkin(7) * t176 + pkin(1) * (t176 * t221 - t195 * t223), pkin(2) * t199 + pkin(7) * t198 + pkin(1) * (t198 * t221 + t199 * t223) + t102, -pkin(2) * t154 + pkin(7) * t102 + pkin(1) * (t102 * t221 - t154 * t223), t227 * (t167 * t222 - t190 * t276) + t229 * (t167 * t220 + t190 * t275), t227 * (-t138 * t222 + t220 * t255) + t229 * (-t138 * t220 - t222 * t255), t227 * (-t179 * t220 + t341) + t229 * (t179 * t222 + t346), t227 * (t188 * t275 + t220 * t257) + t229 * (t188 * t276 - t222 * t257), t227 * (t178 * t222 - t294) + t229 * (t178 * t220 + t293), (t227 * (-t188 * t222 + t190 * t220) + t229 * (-t188 * t220 - t190 * t222)) * qJD(3), t227 * (-qJ(4) * t115 + t297) + t229 * (-pkin(3) * t138 + qJ(4) * t116 - t296) - pkin(2) * t138 + pkin(7) * t71 + pkin(1) * (-t138 * t223 + t221 * t71), t227 * (-qJ(4) * t126 + t296) + t229 * (pkin(3) * t255 + qJ(4) * t127 + t297) + pkin(2) * t255 + pkin(7) * t84 + pkin(1) * (t221 * t84 + t223 * t255), t227 * (-qJ(4) * t100 - t41) + t229 * (-pkin(3) * t134 + qJ(4) * t101 + t42) - pkin(2) * t134 + pkin(7) * t64 + pkin(1) * (-t134 * t223 + t221 * t64), -qJ(4) * t306 + t229 * (-pkin(3) * t119 + qJ(4) * t42) - pkin(2) * t119 + pkin(7) * t19 + pkin(1) * (-t119 * t223 + t19 * t221), t237, -t371, t361, t320, -t382, t319, t227 * (-t220 * t26 + t222 * t39 - t367) + t229 * (t220 * t39 + t222 * t26 + t364) + t368, t227 * (-qJ(4) * t51 - t220 * t27 + t222 * t40) + t229 * (qJ(4) * t53 + t220 * t40 + t222 * t27 - t380) - t381 + pkin(7) * t32 + pkin(1) * (t221 * t32 - t374), t227 * (-qJ(4) * t44 + t12 * t222) + t229 * (qJ(4) * t46 + t12 * t220) + t265 * (-t227 * t44 + t229 * t46) + (t227 * t314 + t229 * t264 + t266) * t59, (t227 * (-pkin(8) * t222 + t314) + t229 * (-pkin(8) * t220 + t264) + t266) * t16 + (t265 + qJ(4)) * (-t227 * t10 + t229 * t11), t237, t361, t371, t319, t382, t320, t227 * (t15 * t222 - t18 * t220 - t367) + t229 * (t15 * t220 + t18 * t222 + t364) + t368, t227 * (-qJ(4) * t43 - t220 * t35 + t222 * t7) + t229 * (-pkin(3) * t58 + qJ(4) * t45 + t220 * t7 + t222 * t35) - pkin(2) * t58 + pkin(7) * t22 + pkin(1) * (t22 * t221 - t223 * t58), t227 * (-qJ(4) * t47 - t13 * t220 + t14 * t222) + t229 * (qJ(4) * t49 + t13 * t222 + t14 * t220 + t380) + t381 + pkin(7) * t30 + pkin(1) * (t221 * t30 + t374), t227 * (-qJ(4) * t5 - t2 * t220 + t222 * t4) + t229 * (-pkin(3) * t8 + qJ(4) * t6 + t2 * t222 + t220 * t4) - pkin(2) * t8 + pkin(7) * t1 + pkin(1) * (t1 * t221 - t223 * t8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t283, 0, 0, 0, 0, 0, 0, t201 * t229 + t207 * t227, -t202 * t227 + t206 * t229, 0, -t142 * t229 + t143 * t227, 0, 0, 0, 0, 0, 0, t115 * t229 + t116 * t227, t126 * t229 + t127 * t227, t100 * t229 + t101 * t227, t227 * t42 + t229 * t41, 0, 0, 0, 0, 0, 0, t362, t227 * t53 + t229 * t51, t227 * t46 + t229 * t44, t10 * t229 + t11 * t227, 0, 0, 0, 0, 0, 0, t362, t227 * t45 + t229 * t43, t227 * t49 + t229 * t47, t227 * t6 + t229 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t208, t213 - t215, t270, t208, t212, qJDD(3), -t142, -t143, 0, 0, t165, t187 - t186, t141, -t165, t139, qJDD(3), pkin(3) * t115 + t183 - t260, pkin(3) * t126 - t67, pkin(3) * t100, pkin(3) * t41, t89, t332, t351, t252, -t366, t242, -pkin(4) * t331 - t305 + t365, pkin(3) * t51 + pkin(4) * t99 + t310 + t377, pkin(3) * t44 + pkin(8) * t63 + t17 + t347, pkin(3) * t10 - pkin(4) * t56 + pkin(8) * t17, t89, t351, -t332, t242, t366, t252, t228 * t25 + t261 * t331 + t365, pkin(3) * t43 + pkin(8) * t61 + t20 * t228 + t21 * t226 + t347, pkin(3) * t47 - t377 + t226 * t24 + (pkin(4) + t313) * t326, pkin(3) * t5 + pkin(8) * t9 + (t261 - t313) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t138, -t255, t134, t119, 0, 0, 0, 0, 0, 0, t334, t76, t59, t16, 0, 0, 0, 0, 0, 0, t334, t58, -t76, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t290, t329, t325, -t290, -t92, t164, -t37, -t38, 0, 0, t290, t325, -t329, t164, t92, -t290, t235, t253, t177 + t321, t254; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t327, t325, t330, t29;];
tauJ_reg  = t3;
