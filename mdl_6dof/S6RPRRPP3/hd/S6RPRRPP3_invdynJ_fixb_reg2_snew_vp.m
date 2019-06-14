% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRRPP3
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
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 21:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRRPP3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP3_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:27:24
% EndTime: 2019-05-05 21:27:37
% DurationCPUTime: 4.72s
% Computational Cost: add. (7219->347), mult. (13911->395), div. (0->0), fcn. (8687->8), ass. (0->209)
t190 = sin(pkin(9));
t191 = cos(pkin(9));
t195 = sin(qJ(3));
t198 = cos(qJ(3));
t194 = sin(qJ(4));
t176 = qJD(1) * t198 - qJD(4);
t174 = t176 ^ 2;
t197 = cos(qJ(4));
t245 = qJD(1) * t195;
t154 = qJD(3) * t194 + t197 * t245;
t284 = t154 ^ 2;
t291 = -t284 - t174;
t240 = qJD(1) * qJD(3);
t180 = t195 * t240;
t238 = t198 * qJDD(1);
t161 = -t180 + t238;
t151 = -qJDD(4) + t161;
t152 = -t197 * qJD(3) + t194 * t245;
t254 = t154 * t152;
t212 = t151 - t254;
t297 = t197 * t212;
t311 = t194 * t291 - t297;
t232 = t198 * t240;
t239 = t195 * qJDD(1);
t160 = t232 + t239;
t216 = qJDD(3) * t194 + t160 * t197;
t81 = (qJD(4) - t176) * t152 - t216;
t329 = t195 * t81 + t198 * t311;
t299 = t194 * t212;
t56 = t197 * t291 + t299;
t350 = pkin(1) * (t190 * t329 + t191 * t56) + pkin(2) * t56 + pkin(7) * t329;
t211 = t151 + t254;
t263 = t197 * t211;
t285 = t152 ^ 2;
t292 = -t174 - t285;
t318 = t292 * t194 - t263;
t107 = qJD(4) * t154 - t197 * qJDD(3) + t160 * t194;
t253 = t154 * t176;
t295 = t107 - t253;
t269 = t194 * t211;
t319 = t292 * t197 + t269;
t328 = t195 * t295 + t198 * t319;
t347 = pkin(1) * (t190 * t328 - t191 * t318) - pkin(2) * t318 + pkin(7) * t328;
t243 = qJD(5) * t176;
t345 = pkin(3) * t56;
t354 = 0.2e1 * t243 + t345;
t128 = -t284 + t174;
t314 = t195 * (t128 * t194 + t263);
t77 = (-qJD(4) - t176) * t152 + t216;
t353 = -t198 * t77 - t314;
t344 = pkin(8) * t56;
t342 = pkin(3) * t318;
t341 = pkin(8) * t311;
t340 = pkin(8) * t318;
t339 = pkin(8) * t319;
t320 = t197 * t128 - t269;
t335 = -pkin(3) * t295 + t339;
t334 = pkin(3) * t81 - t341;
t127 = t285 - t174;
t294 = t107 + t253;
t333 = t195 * (t127 * t197 + t299) + t198 * t294;
t332 = t195 * t311 - t198 * t81;
t331 = t195 * t319 - t198 * t295;
t267 = t197 * t295;
t115 = t285 - t284;
t321 = t198 * t115;
t330 = t195 * (t194 * t81 - t267) + t321;
t90 = t284 + t285;
t326 = pkin(3) * t90;
t325 = (t295 - t253) * pkin(4);
t324 = qJ(5) * t90;
t323 = t195 * t90;
t322 = t198 * t90;
t133 = t152 * t176;
t209 = qJD(4) * t152 - t216;
t79 = t133 - t209;
t248 = -g(3) + qJDD(2);
t179 = t198 * t248;
t200 = qJD(1) ^ 2;
t196 = sin(qJ(1));
t199 = cos(qJ(1));
t230 = t196 * g(1) - g(2) * t199;
t156 = qJDD(1) * pkin(1) + t230;
t220 = g(1) * t199 + g(2) * t196;
t157 = -pkin(1) * t200 - t220;
t246 = t190 * t156 + t191 * t157;
t102 = -pkin(2) * t200 + qJDD(1) * pkin(7) + t246;
t222 = -pkin(3) * t198 - pkin(8) * t195;
t226 = t200 * t222 + t102;
t283 = qJD(3) ^ 2;
t63 = -qJDD(3) * pkin(3) - t283 * pkin(8) + t195 * t226 - t179;
t206 = t107 * pkin(4) - t79 * qJ(5) + t63;
t204 = 0.2e1 * qJD(5) * t154 - t206;
t274 = t194 * t294;
t80 = -t133 - t209;
t40 = -t197 * t80 - t274;
t275 = t194 * t295;
t317 = t197 * t81 + t275;
t313 = t127 * t194 - t297;
t69 = t197 * t294;
t307 = t194 * t80 - t69;
t310 = pkin(8) * t307 + t326;
t27 = t198 * t307 - t323;
t309 = t195 * t307 + t322;
t306 = qJ(5) * t292;
t300 = qJ(5) * t81;
t293 = t285 * pkin(5) - 0.2e1 * qJD(6) * t152;
t278 = pkin(4) + qJ(6);
t70 = qJ(5) * t294;
t289 = -t278 * t80 - t70;
t87 = qJ(5) * t212;
t288 = -t278 * t291 - t87;
t111 = pkin(4) * t152 - qJ(5) * t154;
t227 = t191 * t156 - t190 * t157;
t101 = -qJDD(1) * pkin(2) - t200 * pkin(7) - t227;
t218 = -t161 + t180;
t219 = t160 + t232;
t55 = pkin(3) * t218 - pkin(8) * t219 + t101;
t228 = t195 * t248;
t64 = -t283 * pkin(3) + qJDD(3) * pkin(8) + t198 * t226 + t228;
t30 = t194 * t64 - t197 * t55;
t24 = t151 * pkin(4) - t174 * qJ(5) + t111 * t154 + qJDD(5) + t30;
t205 = -t209 * pkin(5) + t211 * qJ(6) + t24;
t282 = 0.2e1 * qJD(6);
t14 = (-pkin(5) * t152 + t282) * t176 + t205;
t169 = -0.2e1 * t243;
t125 = pkin(5) * t154 + qJ(6) * t176;
t31 = t194 * t55 + t197 * t64;
t223 = -t174 * pkin(4) - t151 * qJ(5) - t152 * t111 + t31;
t208 = -t107 * pkin(5) - qJ(6) * t285 - t176 * t125 + qJDD(6) + t223;
t17 = t169 + t208;
t287 = qJ(5) * t17 - t278 * t14;
t252 = t176 * t194;
t126 = t154 * t252;
t236 = t198 * t254;
t224 = t195 * (-t197 * t209 + t126) - t236;
t286 = t211 * t278 - t306;
t281 = pkin(4) * t176;
t280 = pkin(4) * t194;
t279 = pkin(4) * t197;
t276 = t194 * t63;
t268 = t197 * t63;
t259 = qJ(5) * t197;
t258 = qJ(6) * t107;
t251 = t176 * t197;
t175 = t198 * t200 * t195;
t166 = qJDD(3) + t175;
t250 = t195 * t166;
t165 = -t175 + qJDD(3);
t249 = t198 * t165;
t237 = t152 * t251;
t234 = pkin(1) * t190 + pkin(7);
t233 = -pkin(4) * t77 - t70;
t231 = qJ(5) * t194 + pkin(3);
t16 = t194 * t30 + t197 * t31;
t83 = t102 * t195 - t179;
t84 = t198 * t102 + t228;
t45 = t195 * t83 + t198 * t84;
t229 = -0.2e1 * qJD(5) - t281;
t225 = t195 * (t107 * t194 - t237) + t236;
t67 = -t197 * t107 - t152 * t252;
t68 = -t154 * t251 - t194 * t209;
t23 = t169 + t223;
t221 = -pkin(4) * t24 + qJ(5) * t23;
t217 = t194 * t31 - t197 * t30;
t215 = (t152 * t194 + t154 * t197) * t176;
t214 = -pkin(1) * t191 - pkin(2) + t222;
t136 = t198 * t151;
t213 = t195 * (-t126 + t237) + t136;
t210 = -pkin(4) * t291 + t223 - t87;
t203 = pkin(4) * t211 + t24 - t306;
t202 = t204 + t293;
t201 = t176 * t282 + t205;
t187 = t198 ^ 2;
t186 = t195 ^ 2;
t184 = t187 * t200;
t182 = t186 * t200;
t173 = -t184 - t283;
t172 = -t182 - t283;
t164 = t182 + t184;
t163 = (t186 + t187) * qJDD(1);
t162 = -0.2e1 * t180 + t238;
t159 = 0.2e1 * t232 + t239;
t124 = -t172 * t195 - t249;
t123 = t173 * t198 - t250;
t46 = pkin(5) * t211 - qJ(5) * t295;
t42 = t194 * t77 - t69;
t39 = -t197 * t77 - t274;
t32 = -pkin(5) * t212 - t278 * t81;
t26 = t198 * t42 - t323;
t25 = t154 * t229 + t206;
t22 = -t204 + t325;
t21 = pkin(4) * t253 + t204 - t300;
t20 = t24 + t324;
t19 = pkin(4) * t90 + t23;
t18 = t258 + (-t125 + t229) * t154 + t206 - t293;
t13 = pkin(5) * t291 - t300 + (t125 + t281) * t154 + t202 - t258;
t11 = t324 + (t80 - t133) * pkin(5) + t201;
t10 = t202 - t325 + (-t107 - t295) * qJ(6) + pkin(5) * t292 + t125 * t154;
t9 = -pkin(5) * t294 + t278 * t90 + t17;
t8 = t194 * t24 + t197 * t23;
t7 = t194 * t23 - t197 * t24;
t6 = pkin(5) * t14 - qJ(5) * t18;
t5 = t14 * t194 + t17 * t197;
t4 = -t14 * t197 + t17 * t194;
t3 = t195 * t25 + t198 * t8;
t2 = pkin(5) * t17 - t278 * t18;
t1 = t18 * t195 + t198 * t5;
t12 = [0, 0, 0, 0, 0, qJDD(1), t230, t220, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t191 - t190 * t200) + t227, pkin(1) * (-qJDD(1) * t190 - t191 * t200) - t246, 0, pkin(1) * (t190 * t246 + t191 * t227), t219 * t195, t159 * t198 + t162 * t195, t250 + t198 * (-t182 + t283), -t218 * t198, t195 * (t184 - t283) + t249, 0, -t198 * t101 + pkin(2) * t162 + pkin(7) * t123 + pkin(1) * (t123 * t190 + t162 * t191), t195 * t101 - pkin(2) * t159 + pkin(7) * t124 + pkin(1) * (t124 * t190 - t159 * t191), pkin(2) * t164 + pkin(7) * t163 + pkin(1) * (t163 * t190 + t164 * t191) + t45, -pkin(2) * t101 + pkin(7) * t45 + pkin(1) * (-t101 * t191 + t190 * t45), t224, t195 * (-t194 * t79 - t267) + t321, -t198 * t80 - t314, t225, t333, t213, t195 * (t276 - t340) + t198 * (t30 - t342) + t347, t195 * (t268 - t344) + t198 * (t31 - t345) - t350, -t195 * t217 + t214 * t40 + t234 * t27, t234 * (t16 * t198 + t195 * t63) + t214 * t217, t136 + t195 * (t152 * t197 - t154 * t194) * t176, -t353, -t333, t224, t330, t225, t195 * (-pkin(8) * t39 - t19 * t194 + t197 * t20) + t198 * (-pkin(3) * t39 - t233) - pkin(2) * t39 + pkin(7) * t26 + pkin(1) * (t190 * t26 - t191 * t39), t195 * (-t194 * t22 + t259 * t295 + t340) + t198 * (-t203 + t342) - t347, t195 * (t197 * t21 + t280 * t81 + t344) + t198 * (-t210 + t354) + t350, t195 * (-pkin(8) * t7 + (-t259 + t280) * t25) + t198 * (-pkin(3) * t7 - t221) - pkin(2) * t7 + pkin(7) * t3 + pkin(1) * (t190 * t3 - t191 * t7), t213, -t333, t353, t225, -t330, t224, t195 * (-pkin(8) * t40 + t11 * t197 - t194 * t9) + t198 * (-pkin(3) * t40 - t289) - pkin(2) * t40 + pkin(7) * t27 + pkin(1) * (t190 * t27 - t191 * t40), t195 * (t13 * t197 - t194 * t32 + t344) + t198 * (-t208 - t288 + t354) + t350, t195 * (-t10 * t194 + t197 * t46 - t340) + t198 * (t14 + t286 - t342) + t347, t195 * (-pkin(8) * t4 - t194 * t2 + t197 * t6) + t198 * (-pkin(3) * t4 - t287) - pkin(2) * t4 + pkin(7) * t1 + pkin(1) * (t1 * t190 - t191 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t248, 0, 0, 0, 0, 0, 0, t166 * t198 + t173 * t195, -t165 * t195 + t172 * t198, 0, t195 * t84 - t198 * t83, 0, 0, 0, 0, 0, 0, t331, -t332, t309, t16 * t195 - t198 * t63, 0, 0, 0, 0, 0, 0, t195 * t42 + t322, -t331, t332, t195 * t8 - t198 * t25, 0, 0, 0, 0, 0, 0, t309, t332, t331, -t18 * t198 + t195 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t175, t182 - t184, t239, t175, t238, qJDD(3), -t83, -t84, 0, 0, t68, t197 * t79 - t275, t320, t67, t313, t215, -t268 + t335, t276 + t334, t16 + t310, -pkin(3) * t63 + pkin(8) * t16, t215, -t320, -t313, t68, -t317, t67, pkin(8) * t42 + t19 * t197 + t194 * t20 + t326, t197 * t22 + t231 * t295 - t339, t341 + t194 * t21 - (pkin(3) + t279) * t81, pkin(8) * t8 + (-t231 - t279) * t25, t215, -t313, t320, t67, t317, t68, t11 * t194 + t197 * t9 + t310, t13 * t194 + t197 * t32 - t334, t10 * t197 + t194 * t46 + t335, -pkin(3) * t18 + pkin(8) * t5 + t194 * t6 + t197 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t254, -t115, t80, -t254, -t294, -t151, -t30, -t31, 0, 0, -t151, -t77, t294, t254, -t115, -t254, t233, t203, t169 + t210, t221, -t151, t294, t77, -t254, t115, t254, t289, t17 + t288, pkin(5) * t133 - t201 - t286, t287; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, -t211, t291, t24, 0, 0, 0, 0, 0, 0, t80, t291, t211, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t294, -t212, t292, t17;];
tauJ_reg  = t12;
