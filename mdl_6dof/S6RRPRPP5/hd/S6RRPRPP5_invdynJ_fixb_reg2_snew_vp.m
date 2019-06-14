% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RRPRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 12:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RRPRPP5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP5_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:44:49
% EndTime: 2019-05-06 12:45:01
% DurationCPUTime: 4.44s
% Computational Cost: add. (7789->360), mult. (16028->364), div. (0->0), fcn. (9215->6), ass. (0->215)
t314 = -pkin(2) - pkin(8);
t198 = sin(qJ(4));
t201 = cos(qJ(4));
t202 = cos(qJ(2));
t260 = qJD(1) * t202;
t161 = qJD(2) * t201 - t198 * t260;
t154 = t161 ^ 2;
t199 = sin(qJ(2));
t261 = qJD(1) * t199;
t183 = qJD(4) + t261;
t316 = t183 ^ 2;
t107 = -t316 - t154;
t159 = qJD(2) * t198 + t201 * t260;
t121 = t161 * t159;
t256 = qJD(1) * qJD(2);
t247 = t202 * t256;
t253 = t199 * qJDD(1);
t165 = t247 + t253;
t155 = qJDD(4) + t165;
t330 = t121 + t155;
t344 = t198 * t330;
t49 = t107 * t201 - t344;
t381 = t314 * t49;
t188 = t202 * qJDD(1);
t248 = t199 * t256;
t233 = -t188 + t248;
t216 = t201 * qJDD(2) + t198 * t233;
t104 = -t159 * qJD(4) + t216;
t276 = t159 * t183;
t78 = t276 - t104;
t387 = qJ(3) * t78 - t381;
t268 = t199 * qJ(3);
t224 = t202 * t314 - pkin(1) - t268;
t327 = t276 + t104;
t291 = t201 * t327;
t241 = qJDD(2) * t198 - t201 * t233;
t68 = (qJD(4) - t183) * t161 + t241;
t30 = t198 * t68 + t291;
t296 = t199 * t30;
t317 = t159 ^ 2;
t87 = -t317 - t154;
t363 = t202 * t87;
t57 = t198 * t327;
t386 = -pkin(7) * (t296 - t363) - t224 * (t201 * t68 - t57);
t339 = t201 * t330;
t375 = t224 * (t107 * t198 + t339);
t382 = t199 * t49;
t385 = t375 + pkin(7) * (t202 * t78 - t382);
t383 = pkin(3) * t49;
t117 = -t317 + t154;
t292 = t201 * t78;
t103 = qJD(4) * t161 + t241;
t272 = t183 * t161;
t332 = t103 + t272;
t369 = t199 * t117 + t202 * (t198 * t332 + t292);
t366 = qJ(3) * t87;
t380 = -t314 * t30 + t366;
t325 = t317 - t316;
t67 = t103 - t272;
t379 = t199 * t67 + t202 * (t198 * t325 + t339);
t132 = t154 - t316;
t331 = -t121 + t155;
t343 = t198 * t331;
t378 = -t199 * t327 + t202 * (-t132 * t201 + t343);
t326 = -t316 - t317;
t338 = t201 * t331;
t359 = t198 * t326 + t338;
t367 = pkin(3) * t359;
t361 = t132 * t198 + t338;
t300 = t198 * t78;
t374 = t201 * t332 - t300;
t370 = -t201 * t325 + t344;
t113 = pkin(4) * t159 - qJ(5) * t161;
t204 = qJD(2) ^ 2;
t307 = g(3) * t202;
t228 = -qJDD(2) * pkin(2) - t204 * qJ(3) + qJDD(3) + t307;
t205 = qJD(1) ^ 2;
t200 = sin(qJ(1));
t203 = cos(qJ(1));
t237 = g(1) * t203 + g(2) * t200;
t283 = qJDD(1) * pkin(7);
t146 = -pkin(1) * t205 - t237 + t283;
t306 = t202 * pkin(2);
t235 = -t268 - t306;
t242 = t205 * t235 + t146;
t265 = t202 * t205;
t206 = -qJDD(2) * pkin(8) + (t165 - t247) * pkin(3) + (-pkin(8) * t265 + t242) * t199 + t228;
t173 = pkin(3) * t261 - qJD(2) * pkin(8);
t249 = qJD(3) * t261;
t186 = -0.2e1 * t249;
t195 = t202 ^ 2;
t245 = t200 * g(1) - t203 * g(2);
t227 = -qJDD(1) * pkin(1) - t245;
t211 = -t227 + (-t233 - t248) * pkin(2);
t234 = t165 + t247;
t333 = t234 * qJ(3);
t38 = t186 + t233 * pkin(8) - t173 * t261 + (-t195 * pkin(3) - pkin(7)) * t205 - t333 - t211;
t28 = t198 * t206 + t201 * t38;
t238 = t155 * qJ(5) + 0.2e1 * qJD(5) * t183 - t159 * t113 + t28;
t18 = -(t316 + t87) * pkin(4) + t238;
t354 = qJ(3) * t332 + t314 * t359;
t351 = t224 * (t201 * t326 - t343) + pkin(7) * (t199 * t359 + t202 * t332);
t368 = pkin(3) * t87;
t365 = qJ(5) * t78;
t364 = qJ(5) * t87;
t334 = -pkin(4) * t107 + qJ(5) * t330;
t350 = 0.2e1 * qJD(3);
t312 = pkin(4) + pkin(5);
t311 = pkin(3) * t332;
t194 = t199 ^ 2;
t189 = t194 * t205;
t176 = -t189 - t204;
t250 = t199 * t265;
t171 = -qJDD(2) + t250;
t266 = t202 * t171;
t349 = pkin(7) * (-t176 * t199 + t266);
t347 = qJ(6) * t327;
t336 = pkin(4) * t331 + qJ(5) * t326;
t313 = pkin(3) + pkin(7);
t114 = t199 * t121;
t271 = t183 * t198;
t252 = t159 * t271;
t219 = -t114 + t202 * (t103 * t201 - t252);
t191 = t199 * g(3);
t221 = -t204 * pkin(2) + t242 * t202 - t191;
t213 = qJD(2) * t350 + t221;
t269 = t195 * t205;
t329 = t266 - (-t204 + t269) * t199;
t124 = -pkin(5) * t183 - qJ(6) * t161;
t324 = t161 * t124 + qJDD(6);
t323 = t103 * pkin(4) + t365;
t322 = t103 * qJ(6) + 0.2e1 * qJD(6) * t159 + t183 * t124;
t254 = qJDD(2) * qJ(3);
t321 = -t233 * pkin(3) - pkin(8) * t269 + t254;
t320 = -pkin(5) * t103 + t324;
t166 = t188 - 0.2e1 * t248;
t177 = t204 + t269;
t170 = qJDD(2) + t250;
t275 = t170 * t199;
t319 = pkin(7) * (t177 * t202 + t275) - pkin(1) * t166;
t27 = t198 * t38 - t201 * t206;
t225 = -t155 * pkin(4) - qJ(5) * t316 + qJDD(5) + t27;
t218 = -t155 * pkin(5) + t225 - t347;
t315 = 0.2e1 * t161;
t215 = qJD(6) * t315 - t218;
t281 = t113 * t161;
t318 = t215 - t281 + (t331 - t121) * pkin(5) + t336;
t270 = t183 * t201;
t128 = t161 * t270;
t286 = t114 + t202 * (-t104 * t198 - t128);
t310 = pkin(3) * t78;
t305 = t205 * pkin(7);
t23 = -pkin(4) * t316 + t238;
t25 = t225 + t281;
t304 = -pkin(4) * t25 + qJ(5) * t23;
t303 = -pkin(4) * t327 - qJ(5) * t67;
t47 = (t350 + t173) * qJD(2) + t221 + t321;
t302 = t198 * t47;
t294 = t201 * t47;
t127 = t161 * t271;
t56 = t201 * t104 - t127;
t285 = qJ(5) * t198;
t284 = qJ(5) * t201;
t267 = t199 * t166;
t264 = -t107 - t317;
t168 = t189 + t269;
t263 = pkin(1) * t168 + (t194 + t195) * t283;
t251 = t159 * t270;
t246 = -pkin(5) * t159 - t113;
t125 = t146 * t199 + t307;
t126 = t146 * t202 - t191;
t243 = t125 * t199 + t202 * t126;
t14 = (-0.2e1 * qJD(6) - t246) * t161 + t218;
t217 = t23 + t322;
t16 = -pkin(5) * t317 + t217;
t240 = qJ(5) * t16 - t312 * t14;
t239 = qJ(5) * t68 + t312 * t327;
t236 = t127 - t251;
t11 = t198 * t28 - t201 * t27;
t231 = t198 * t27 + t201 * t28;
t230 = (-t189 + t204) * t202 + t275;
t226 = t103 * t198 + t251;
t141 = t199 * t155;
t222 = t141 + t202 * (t128 + t252);
t220 = t23 + t334;
t212 = -t25 + t336;
t80 = t242 * t199 + t228;
t79 = t213 + t254;
t209 = t211 + t305;
t208 = t209 + 0.2e1 * t249;
t207 = -qJD(2) * t173 + qJD(5) * t315 - t213 - t321 - t323;
t26 = (pkin(4) * t183 - 0.2e1 * qJD(5)) * t161 + t47 + t323;
t21 = t207 + (-t332 - t272) * pkin(4);
t20 = -pkin(4) * t272 + t207 - t365;
t169 = t189 - t269;
t164 = 0.2e1 * t247 + t253;
t145 = -t227 + t305;
t123 = t234 * t199;
t122 = t166 * t202;
t110 = t202 * t164 + t267;
t74 = (-qJD(4) - t183) * t159 + t216;
t36 = -qJ(5) * t332 + qJ(6) * t331;
t31 = -t198 * t67 - t291;
t29 = -qJ(6) * t330 - t312 * t78;
t19 = t25 - t364;
t17 = qJ(6) * t317 + t26 - t320;
t10 = t264 * qJ(6) + t20 + t320;
t9 = t21 + (-t103 - t332) * pkin(5) + (-t326 - t317) * qJ(6) + t324;
t8 = t246 * t161 + t215 + t347 + t364;
t7 = -qJ(6) * t68 + (t317 + t87) * pkin(5) - t322 - t18;
t5 = t198 * t23 - t201 * t25;
t4 = -qJ(5) * t17 - qJ(6) * t14;
t2 = -t14 * t201 + t16 * t198;
t1 = -qJ(6) * t16 - t312 * t17;
t3 = [0, 0, 0, 0, 0, qJDD(1), t245, t237, 0, 0, t123, t110, t230, t122, -t329, 0, t202 * t145 - t319, -pkin(1) * t164 - t199 * t145 + t349, t243 + t263, pkin(1) * t145 + pkin(7) * t243, 0, -t230, t329, t123, t110, t122, (pkin(2) * t168 + t79) * t202 + (qJ(3) * t168 + t80) * t199 + t263, t202 * (-pkin(2) * t166 + t186 - t209) + (-t202 * t234 - t267) * qJ(3) + t319, t199 * t208 - t349 + (pkin(1) + t306) * t164 + (t164 + t234) * t268, pkin(7) * (t199 * t80 + t202 * t79) + (pkin(1) - t235) * (t208 + t333), t286, t369, -t378, t219, -t379, t222, t199 * (-t27 + t367) + t202 * (t294 + t311) + t351, t199 * (-t28 + t383) + t202 * (pkin(3) * t74 - t302) + pkin(7) * (t202 * t74 + t382) - t375, -pkin(3) * t296 + t202 * (-t231 + t368) + t386, t224 * t231 + t313 * (t11 * t199 + t202 * t47), t286, -t378, -t369, t222, t379, t219, t199 * (t212 + t367) + t202 * (-t201 * t21 + t285 * t332 + t311) + t351, t199 * (pkin(3) * t31 + t303) + t202 * (-t18 * t201 - t19 * t198 + t368) + pkin(7) * (t199 * t31 + t363) + t224 * (-t201 * t67 + t57), t199 * (t220 - t383) + t202 * (pkin(4) * t292 - t198 * t20 + t310) + t385, (t313 * t5 + t304) * t199 + t224 * (t198 * t25 + t201 * t23) + (pkin(4) * t201 + t285 + t313) * t202 * t26, t286, -t369, t378, t219, -t379, t141 + t202 * (t159 * t198 + t161 * t201) * t183, t199 * (t318 + t367) + t202 * (-t198 * t36 - t201 * t9 + t311) + t351, t199 * (-pkin(5) * t107 + t16 + t334 - t383) + t202 * (-t10 * t198 - t201 * t29 + t310) + t385, t199 * (pkin(3) * t30 + t239) + t202 * (-t198 * t8 - t201 * t7 - t368) - t386, t199 * (pkin(3) * t2 + t240) + t202 * (pkin(3) * t17 - t1 * t201 - t198 * t4) + pkin(7) * (t17 * t202 + t199 * t2) + t224 * (t14 * t198 + t16 * t201); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t250, t169, t253, t250, t188, qJDD(2), -t125, -t126, 0, 0, qJDD(2), -t253, -t188, -t250, t169, t250, (-pkin(2) * t199 + qJ(3) * t202) * qJDD(1), -pkin(2) * t170 + qJ(3) * t177 + t80, -pkin(2) * t176 + (qJDD(2) - t171) * qJ(3) + t213, -pkin(2) * t80 + qJ(3) * t79, t56, -t374, t361, t226, -t370, t236, t302 + t354, qJ(3) * t74 + t294 + t381, -t11 + t380, qJ(3) * t47 + t314 * t11, t56, t361, t374, t236, t370, t226, -t198 * t21 - t284 * t332 + t354, -t18 * t198 + t19 * t201 + t314 * t31 + t366, pkin(4) * t300 + t20 * t201 + t387, t314 * t5 + (pkin(4) * t198 + qJ(3) - t284) * t26, t56, t374, -t361, t226, -t370, (-t159 * t201 + t161 * t198) * t183, -t198 * t9 + t201 * t36 + t354, t10 * t201 - t198 * t29 + t387, -t198 * t7 + t201 * t8 - t380, qJ(3) * t17 - t1 * t198 + t314 * t2 + t201 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t253, t170, t176, t80, 0, 0, 0, 0, 0, 0, t359, t49, -t30, t11, 0, 0, 0, 0, 0, 0, t359, t31, -t49, t5, 0, 0, 0, 0, 0, 0, t359, -t49, t30, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, t117, t327, -t121, -t67, t155, -t27, -t28, 0, 0, t121, t327, -t117, t155, t67, -t121, t212, t303, t220, t304, t121, -t117, -t327, -t121, -t67, t155, t318, t264 * pkin(5) + t217 + t334, t239, t240; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t331, t327, t107, t25, 0, 0, 0, 0, 0, 0, -t331, t107, -t327, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t332, -t78, t87, -t17;];
tauJ_reg  = t3;
