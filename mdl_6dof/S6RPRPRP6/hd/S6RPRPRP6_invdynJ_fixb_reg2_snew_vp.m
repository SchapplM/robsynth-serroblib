% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRPRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 17:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRPRP6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP6_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:54:21
% EndTime: 2019-05-05 17:54:37
% DurationCPUTime: 6.76s
% Computational Cost: add. (15038->403), mult. (36660->515), div. (0->0), fcn. (26603->8), ass. (0->243)
t245 = qJD(3) ^ 2;
t238 = cos(pkin(9));
t241 = sin(qJ(3));
t237 = sin(pkin(9));
t244 = cos(qJ(3));
t300 = t237 * t244;
t264 = t238 * t241 + t300;
t219 = t264 * qJD(1);
t336 = t219 ^ 2;
t198 = t336 + t245;
t301 = t237 * t241;
t217 = (-t238 * t244 + t301) * qJD(1);
t302 = t219 * t217;
t355 = qJDD(3) + t302;
t366 = t355 * t241;
t132 = t198 * t244 + t366;
t365 = t355 * t244;
t134 = -t198 * t241 + t365;
t380 = qJ(2) * (t132 * t237 - t134 * t238);
t379 = pkin(7) * t132;
t378 = pkin(7) * t134;
t199 = t336 - t245;
t168 = t302 - qJDD(3);
t373 = t168 * t244;
t374 = t168 * t241;
t376 = t237 * (-t199 * t241 + t373) + t238 * (t199 * t244 + t374);
t337 = t217 ^ 2;
t194 = t337 - t245;
t375 = t237 * (-t194 * t244 + t366) - t238 * (t194 * t241 + t365);
t240 = sin(qJ(5));
t243 = cos(qJ(5));
t190 = t240 * qJD(3) - t243 * t217;
t192 = t243 * qJD(3) + t240 * t217;
t146 = t192 * t190;
t205 = t217 * qJD(3);
t216 = t264 * qJDD(1);
t177 = t216 - t205;
t164 = qJDD(5) + t177;
t351 = t146 - t164;
t368 = t351 * pkin(5);
t166 = -t245 - t337;
t116 = t166 * t241 - t373;
t119 = -t166 * t244 - t374;
t367 = qJ(2) * (t116 * t237 + t119 * t238);
t364 = pkin(7) * t116;
t363 = pkin(7) * t119;
t141 = -t205 + t177;
t347 = qJ(4) * t141;
t246 = qJD(1) ^ 2;
t242 = sin(qJ(1));
t331 = cos(qJ(1));
t260 = t331 * g(1) + t242 * g(2);
t350 = -t246 * pkin(1) + qJDD(1) * qJ(2) + 0.2e1 * qJD(1) * qJD(2) - t260;
t232 = t237 ^ 2;
t233 = t238 ^ 2;
t298 = t232 + t233;
t279 = t242 * g(1) - t331 * g(2);
t269 = -qJDD(2) + t279;
t329 = t238 * pkin(2);
t167 = (t298 * pkin(7) + qJ(2)) * t246 + (pkin(1) + t329) * qJDD(1) + t269;
t294 = t219 * qJD(3);
t356 = pkin(3) * t294 - 0.2e1 * qJD(4) * t219 - t167;
t335 = pkin(3) + pkin(8);
t341 = -t337 - t336;
t353 = pkin(1) * t341;
t352 = pkin(2) * t341;
t314 = t351 * t240;
t313 = t351 * t243;
t188 = t190 ^ 2;
t209 = qJD(5) + t219;
t207 = t209 ^ 2;
t129 = -t207 - t188;
t80 = t129 * t240 - t313;
t81 = t129 * t243 + t314;
t349 = pkin(4) * t80 - qJ(4) * t81;
t189 = t192 ^ 2;
t136 = -t189 - t207;
t113 = t146 + t164;
t316 = t113 * t240;
t87 = t136 * t243 - t316;
t315 = t113 * t243;
t88 = -t136 * t240 - t315;
t348 = pkin(4) * t87 - qJ(4) * t88;
t290 = t238 * qJDD(1);
t291 = t237 * qJDD(1);
t267 = t241 * t291 - t244 * t290;
t175 = t267 + t294;
t128 = -qJD(5) * t190 + qJDD(3) * t243 + t175 * t240;
t156 = t209 * t190;
t346 = t128 - t156;
t140 = t205 + t177;
t193 = t219 * pkin(4) - qJD(3) * pkin(8);
t249 = -t347 + t356;
t57 = -pkin(4) * t337 + t335 * t175 - t219 * t193 + t249;
t328 = t238 * g(3);
t151 = -t328 + (-pkin(7) * qJDD(1) + t246 * t329 - t350) * t237;
t270 = -t237 * g(3) + t350 * t238;
t155 = -t233 * t246 * pkin(2) + pkin(7) * t290 + t270;
t107 = -t244 * t151 + t241 * t155;
t163 = t217 * pkin(3) - t219 * qJ(4);
t83 = -qJDD(3) * pkin(3) - t245 * qJ(4) + t219 * t163 + qJDD(4) + t107;
t60 = t140 * pkin(4) + pkin(8) * t168 + t83;
t37 = t240 * t57 - t243 * t60;
t343 = qJ(6) * t156 + 0.2e1 * qJD(6) * t192 + t368 + t37;
t275 = t240 * qJDD(3) - t243 * t175;
t127 = -qJD(5) * t192 - t275;
t149 = t209 * pkin(5) - t192 * qJ(6);
t38 = t240 * t60 + t243 * t57;
t26 = -t188 * pkin(5) + t127 * qJ(6) - 0.2e1 * qJD(6) * t190 - t209 * t149 + t38;
t342 = t336 - t337;
t299 = t246 * qJ(2);
t317 = qJDD(1) * pkin(1);
t210 = t269 + t299 + t317;
t340 = t298 * t299 - t210 - t317;
t124 = -t188 - t189;
t104 = t128 + t156;
t257 = (-qJD(5) + t209) * t192 - t275;
t66 = -t104 * t243 + t240 * t257;
t47 = t124 * t241 - t244 * t66;
t334 = pkin(7) * t47;
t100 = (qJD(5) + t209) * t192 + t275;
t50 = t100 * t241 - t244 * t80;
t333 = pkin(7) * t50;
t54 = t241 * t346 - t244 * t87;
t332 = pkin(7) * t54;
t330 = pkin(3) * t244;
t48 = t124 * t244 + t241 * t66;
t68 = t104 * t240 + t243 * t257;
t327 = qJ(2) * (-t237 * t47 + t238 * t48) - pkin(1) * t68;
t51 = t100 * t244 + t241 * t80;
t326 = qJ(2) * (-t237 * t50 + t238 * t51) - pkin(1) * t81;
t55 = t241 * t87 + t244 * t346;
t325 = qJ(2) * (-t237 * t54 + t238 * t55) - pkin(1) * t88;
t108 = t241 * t151 + t244 * t155;
t74 = -t107 * t244 + t108 * t241;
t322 = t237 * t74;
t25 = -qJ(6) * t128 - t343;
t321 = t240 * t25;
t292 = qJD(4) * qJD(3);
t229 = 0.2e1 * t292;
t265 = -t245 * pkin(3) + qJDD(3) * qJ(4) - t217 * t163 + t108;
t256 = -t175 * pkin(4) - pkin(8) * t337 + qJD(3) * t193 + t265;
t61 = t229 + t256;
t320 = t240 * t61;
t319 = t243 * t25;
t318 = t243 * t61;
t310 = t167 * t241;
t309 = t167 * t244;
t304 = t209 * t240;
t303 = t209 * t243;
t286 = t241 * t146;
t285 = t244 * t146;
t283 = -pkin(2) * t68 + pkin(7) * t48;
t282 = -pkin(2) * t81 + pkin(7) * t51;
t281 = -pkin(2) * t88 + pkin(7) * t55;
t280 = qJ(4) * t241 + pkin(2);
t40 = pkin(4) * t66 - qJ(4) * t68;
t75 = t107 * t241 + t244 * t108;
t277 = t237 * (t350 * t237 + t328) + t238 * t270;
t273 = qJ(4) * t100 - t335 * t80;
t272 = qJ(4) * t346 - t335 * t87;
t271 = qJ(4) * t124 - t335 * t66;
t18 = t240 * t38 - t243 * t37;
t19 = t240 * t37 + t243 * t38;
t263 = pkin(4) * t100 - t335 * t81;
t262 = pkin(4) * t346 - t335 * t88;
t261 = pkin(4) * t124 - t335 * t68;
t82 = t229 + t265;
t258 = pkin(5) * t136 - t26;
t255 = t25 - t368;
t254 = t237 * (t244 * t177 - t241 * t294) + t238 * (t241 * t177 + t244 * t294);
t253 = t237 * (t175 * t241 + t244 * t205) + t238 * (-t244 * t175 + t241 * t205);
t252 = -t127 * pkin(5) - t188 * qJ(6) + t149 * t192 + qJDD(6) + t256;
t251 = (t237 * (-t217 * t244 + t219 * t241) + t238 * (-t217 * t241 - t219 * t244)) * qJD(3);
t41 = t229 + t252;
t250 = -t175 * pkin(3) - t356;
t228 = t233 * qJDD(1);
t227 = t232 * qJDD(1);
t222 = t298 * t246;
t176 = t216 - 0.2e1 * t205;
t174 = t267 + 0.2e1 * t294;
t153 = -t189 + t207;
t152 = t188 - t207;
t142 = t189 - t188;
t139 = t175 + t294;
t138 = t175 - t294;
t123 = t140 * t241 - t244 * t267;
t122 = -t138 * t244 + t216 * t241;
t121 = -t140 * t244 - t241 * t267;
t120 = -t138 * t241 - t216 * t244;
t110 = (-t190 * t243 + t192 * t240) * t209;
t109 = (t190 * t240 + t192 * t243) * t209;
t99 = pkin(5) * t104;
t96 = t128 * t243 - t192 * t304;
t95 = -t128 * t240 - t192 * t303;
t94 = -t127 * t240 + t190 * t303;
t93 = -t127 * t243 - t190 * t304;
t92 = t152 * t243 - t316;
t91 = -t153 * t240 - t313;
t90 = -t152 * t240 - t315;
t89 = -t153 * t243 + t314;
t79 = t250 + t347;
t73 = -pkin(5) * t346 - qJ(6) * t113;
t72 = -qJ(4) * t341 + t83;
t71 = -pkin(3) * t341 + t82;
t70 = (t139 + t175) * pkin(3) + t249;
t69 = t250 + 0.2e1 * t347;
t67 = -t100 * t243 - t240 * t346;
t65 = t100 * t240 - t243 * t346;
t52 = t237 * (-t109 * t241 + t164 * t244) + t238 * (t109 * t244 + t164 * t241);
t43 = t237 * (-t241 * t95 + t285) + t238 * (t244 * t95 + t286);
t42 = t237 * (-t241 * t93 - t285) + t238 * (t244 * t93 - t286);
t39 = -qJ(6) * t136 + t41;
t36 = t237 * (t104 * t244 - t241 * t89) + t238 * (t104 * t241 + t244 * t89);
t35 = t237 * (-t241 * t90 + t244 * t257) + t238 * (t241 * t257 + t244 * t90);
t32 = t40 - t99;
t31 = t237 * (t142 * t244 - t241 * t65) + t238 * (t142 * t241 + t244 * t65);
t30 = -pkin(5) * t100 + qJ(6) * t129 - t252 - 0.2e1 * t292;
t28 = t262 - t320;
t27 = t263 + t318;
t24 = pkin(5) * t25;
t23 = -t38 + t348;
t22 = -t37 + t349;
t21 = (t104 + t128) * qJ(6) + t343;
t20 = -pkin(5) * t124 + qJ(6) * t257 + t26;
t17 = t258 + t348;
t16 = -t240 * t39 - t243 * t73 + t262;
t15 = -pkin(5) * t41 + qJ(6) * t26;
t14 = -qJ(6) * t314 - t243 * t30 + t263;
t13 = t255 + t349;
t12 = t18 * t241 + t244 * t61;
t11 = -t18 * t244 + t241 * t61;
t10 = t243 * t26 - t321;
t9 = t240 * t26 + t319;
t8 = -t19 + t261;
t7 = t241 * t9 + t244 * t41;
t6 = t241 * t41 - t244 * t9;
t5 = pkin(4) * t18 - qJ(4) * t19;
t4 = pkin(4) * t61 - t335 * t19;
t3 = -t243 * t20 - t240 * t21 + t261;
t2 = pkin(4) * t9 - qJ(4) * t10 + t24;
t1 = pkin(4) * t41 + qJ(6) * t321 - t335 * t10 - t243 * t15;
t29 = [0, 0, 0, 0, 0, qJDD(1), t279, t260, 0, 0, t227, 0.2e1 * t237 * t290, 0, t228, 0, 0, -t340 * t238, t340 * t237, pkin(1) * t222 + qJ(2) * (t228 + t227) + t277, pkin(1) * t210 + qJ(2) * t277, t254, t237 * (-t174 * t244 - t176 * t241) + t238 * (-t174 * t241 + t176 * t244), -t376, t253, -t375, t251, t237 * (-t310 - t364) + t238 * (-pkin(2) * t174 + t309 - t363) - pkin(1) * t174 - t367, t237 * (-t309 + t379) + t238 * (-pkin(2) * t176 - t310 - t378) - pkin(1) * t176 + t380, t237 * (-pkin(7) * t120 - t74) + t238 * (pkin(7) * t122 - t352 + t75) - t353 + qJ(2) * (-t120 * t237 + t122 * t238), -pkin(7) * t322 + t238 * (pkin(2) * t167 + pkin(7) * t75) + pkin(1) * t167 + qJ(2) * (t238 * t75 - t322), t251, t376, t375, t254, t237 * (-t139 * t244 - t141 * t241) + t238 * (-t139 * t241 + t141 * t244), t253, t237 * (-pkin(7) * t121 - t241 * t71 + t244 * t72) + t238 * (pkin(7) * t123 + t241 * t72 + t244 * t71 - t352) - t353 + qJ(2) * (-t121 * t237 + t123 * t238), t237 * (-t241 * t70 + t364) + t238 * (t244 * t70 + t363) + t367 + (qJ(4) * t300 + t238 * t280 + pkin(1)) * t139, t237 * (t244 * t69 - t379) + t238 * (t241 * t69 + t378) - t380 + (-pkin(3) * t301 + t238 * (pkin(2) + t330) + pkin(1)) * t141, (t237 * (-pkin(3) * t241 + qJ(4) * t244) + t238 * (t280 + t330) + pkin(1)) * t79 + (qJ(2) + pkin(7)) * (-t237 * (t241 * t82 - t244 * t83) + t238 * (t241 * t83 + t244 * t82)), t43, t31, t36, t42, t35, t52, t237 * (t22 * t244 - t241 * t27 - t333) + t238 * (t22 * t241 + t244 * t27 + t282) + t326, t237 * (t23 * t244 - t241 * t28 - t332) + t238 * (t23 * t241 + t244 * t28 + t281) + t325, t237 * (-t241 * t8 + t244 * t40 - t334) + t238 * (t241 * t40 + t244 * t8 + t283) + t327, t237 * (-pkin(7) * t11 - t241 * t4 + t244 * t5) + t238 * (-pkin(2) * t19 + pkin(7) * t12 + t241 * t5 + t244 * t4) - pkin(1) * t19 + qJ(2) * (-t11 * t237 + t12 * t238), t43, t31, t36, t42, t35, t52, t237 * (t13 * t244 - t14 * t241 - t333) + t238 * (t13 * t241 + t14 * t244 + t282) + t326, t237 * (-t16 * t241 + t17 * t244 - t332) + t238 * (t16 * t244 + t17 * t241 + t281) + t325, t237 * (-t241 * t3 + t244 * t32 - t334) + t238 * (t241 * t32 + t244 * t3 + t283) + t327, t237 * (-pkin(7) * t6 - t1 * t241 + t2 * t244) + t238 * (-pkin(2) * t10 + pkin(7) * t7 + t1 * t244 + t2 * t241) - pkin(1) * t10 + qJ(2) * (-t237 * t6 + t238 * t7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t290, t291, -t222, -t210, 0, 0, 0, 0, 0, 0, t174, t176, t341, -t167, 0, 0, 0, 0, 0, 0, t341, -t139, -t141, -t79, 0, 0, 0, 0, 0, 0, t81, t88, t68, t19, 0, 0, 0, 0, 0, 0, t81, t88, t68, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t302, t342, t216, -t302, -t267, qJDD(3), -t107, -t108, 0, 0, qJDD(3), -t140, t138, t302, t342, -t302, -pkin(3) * t140 - qJ(4) * t267, pkin(3) * t168 - qJ(4) * t166 + t83, pkin(3) * t198 + qJ(4) * t355 + t82, -pkin(3) * t83 + qJ(4) * t82, t96, t67, t91, t94, t92, t110, t273 + t320, t272 + t318, -t18 + t271, qJ(4) * t61 - t335 * t18, t96, t67, t91, t94, t92, t110, qJ(6) * t313 - t240 * t30 + t273, -t240 * t73 + t243 * t39 + t272, -t20 * t240 + t21 * t243 + t271, qJ(4) * t41 - qJ(6) * t319 - t240 * t15 - t335 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140, -t168, -t198, t83, 0, 0, 0, 0, 0, 0, t80, t87, t66, t18, 0, 0, 0, 0, 0, 0, t80, t87, t66, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t146, t142, t104, -t146, t257, t164, -t37, -t38, 0, 0, t146, t142, t104, -t146, t257, t164, t255, t258, -t99, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, t346, t124, t41;];
tauJ_reg  = t29;
