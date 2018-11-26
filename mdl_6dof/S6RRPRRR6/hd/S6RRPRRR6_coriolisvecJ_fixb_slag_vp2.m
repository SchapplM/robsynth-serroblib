% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:24
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR6_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:24:26
% EndTime: 2018-11-23 17:24:35
% DurationCPUTime: 8.63s
% Computational Cost: add. (12959->571), mult. (30237->748), div. (0->0), fcn. (20686->8), ass. (0->260)
t238 = sin(qJ(6));
t347 = -t238 / 0.2e1;
t242 = cos(qJ(6));
t346 = t242 / 0.2e1;
t402 = Ifges(7,5) * t346 + Ifges(7,6) * t347;
t240 = sin(qJ(4));
t244 = cos(qJ(4));
t245 = cos(qJ(2));
t301 = qJD(1) * t245;
t241 = sin(qJ(2));
t302 = qJD(1) * t241;
t163 = -t240 * t302 - t244 * t301;
t165 = -t240 * t301 + t244 * t302;
t239 = sin(qJ(5));
t243 = cos(qJ(5));
t113 = -t243 * t163 + t239 * t165;
t106 = qJD(6) + t113;
t281 = mrSges(7,1) * t238 + mrSges(7,2) * t242;
t367 = qJD(2) - qJD(4);
t230 = qJD(5) - t367;
t226 = pkin(7) * t302;
t192 = pkin(8) * t302 - t226;
t265 = qJD(3) - t192;
t246 = -pkin(2) - pkin(3);
t286 = t246 * qJD(2);
t258 = t286 + t265;
t227 = pkin(7) * t301;
t193 = -pkin(8) * t301 + t227;
t236 = qJD(2) * qJ(3);
t166 = t193 + t236;
t305 = t244 * t166;
t100 = t240 * t258 + t305;
t338 = t163 * pkin(9);
t88 = t100 + t338;
t314 = t239 * t88;
t341 = pkin(9) * t165;
t142 = t244 * t258;
t99 = -t166 * t240 + t142;
t87 = t99 - t341;
t83 = -pkin(4) * t367 + t87;
t46 = t243 * t83 - t314;
t43 = -pkin(5) * t230 - t46;
t401 = t106 * t402 + t43 * t281;
t326 = Ifges(7,4) * t242;
t384 = Ifges(7,2) * t347 + t326 / 0.2e1;
t327 = Ifges(7,4) * t238;
t280 = Ifges(7,1) * t242 - t327;
t400 = t280 / 0.2e1;
t382 = -qJD(2) / 0.2e1;
t398 = -mrSges(3,1) - mrSges(4,1);
t269 = t163 * t239 + t243 * t165;
t94 = t230 * t238 + t242 * t269;
t345 = Ifges(7,4) * t94;
t93 = t230 * t242 - t238 * t269;
t41 = Ifges(7,2) * t93 + Ifges(7,6) * t106 + t345;
t91 = Ifges(7,4) * t93;
t42 = Ifges(7,1) * t94 + Ifges(7,5) * t106 + t91;
t397 = t42 * t346 + t41 * t347;
t77 = pkin(5) * t269 + pkin(10) * t113;
t333 = -mrSges(6,1) * t230 - mrSges(7,1) * t93 + mrSges(7,2) * t94 + mrSges(6,3) * t269;
t197 = -qJ(3) * t240 + t244 * t246;
t191 = -pkin(4) + t197;
t198 = t244 * qJ(3) + t240 * t246;
t127 = t239 * t191 + t243 * t198;
t150 = t244 * qJD(3) + qJD(4) * t197;
t151 = -t240 * qJD(3) - qJD(4) * t198;
t124 = -t192 * t240 + t244 * t193;
t262 = t124 + t338;
t125 = t244 * t192 + t240 * t193;
t95 = t125 + t341;
t396 = -qJD(5) * t127 + (t151 - t262) * t243 + (-t150 + t95) * t239;
t266 = t240 * t241 + t244 * t245;
t131 = t367 * t266;
t115 = t131 * qJD(1);
t340 = t115 * pkin(9);
t313 = t243 * t88;
t47 = t239 * t83 + t313;
t267 = t240 * t245 - t241 * t244;
t130 = t367 * t267;
t116 = t130 * qJD(1);
t357 = pkin(7) - pkin(8);
t211 = t357 * t245;
t195 = qJD(2) * t211;
t170 = t240 * t195;
t299 = qJD(4) * t240;
t373 = t357 * t241;
t74 = qJD(4) * t142 - t166 * t299 + t244 * (-qJD(1) * t373 + qJD(3)) * qJD(2) + qJD(1) * t170;
t51 = -pkin(9) * t116 + t74;
t137 = t244 * t211 + t240 * t373;
t254 = ((-qJD(4) * t246 - qJD(3)) * t240 + t137 * qJD(1)) * qJD(2);
t261 = t240 * t265;
t75 = (-t261 - t305) * qJD(4) + t254;
t13 = t239 * t51 - t243 * (t75 - t340) + t47 * qJD(5);
t378 = qJD(6) * t93;
t55 = -qJD(5) * t113 + t115 * t243 - t116 * t239;
t35 = t242 * t55 + t378;
t377 = qJD(6) * t94;
t36 = -t238 * t55 - t377;
t14 = -mrSges(7,1) * t36 + mrSges(7,2) * t35;
t395 = m(7) * t13 + t14;
t394 = -t124 + t151;
t393 = -t125 + t150;
t167 = -qJD(1) * pkin(1) - pkin(2) * t301 - qJ(3) * t302;
t203 = t227 + t236;
t205 = mrSges(4,2) * t301 + qJD(2) * mrSges(4,3);
t329 = Ifges(3,4) * t241;
t386 = -Ifges(4,5) * t302 / 0.2e1;
t387 = Ifges(4,3) / 0.2e1;
t392 = -(m(4) * t203 - qJD(2) * mrSges(3,2) + mrSges(3,3) * t301 + t205) * pkin(7) - t203 * mrSges(4,2) - t301 * t387 - t386 - (t245 * Ifges(3,2) + t329) * qJD(1) / 0.2e1 + t167 * mrSges(4,1) + (-Ifges(4,6) + Ifges(3,6)) * t382;
t319 = qJD(2) * pkin(2);
t199 = qJD(3) + t226 - t319;
t325 = Ifges(4,5) * t245;
t385 = -Ifges(3,4) * t301 / 0.2e1;
t388 = -Ifges(3,1) / 0.2e1;
t391 = (m(4) * t199 + (mrSges(4,2) + mrSges(3,3)) * t302 + t398 * qJD(2)) * pkin(7) - t167 * mrSges(4,3) + (t241 * Ifges(4,1) - t325) * qJD(1) / 0.2e1 - t302 * t388 - t385 + t199 * mrSges(4,2) - (Ifges(4,4) + Ifges(3,5)) * t382;
t155 = Ifges(5,4) * t165;
t101 = t163 * Ifges(5,2) - Ifges(5,6) * t367 + t155;
t154 = Ifges(5,4) * t163;
t102 = t165 * Ifges(5,1) - Ifges(5,5) * t367 + t154;
t146 = pkin(3) * t301 - t167;
t298 = qJD(4) * t244;
t12 = t46 * qJD(5) + t243 * t51 + (-qJD(4) * t261 - t166 * t298 + t254 - t340) * t239;
t56 = qJD(5) * t269 + t115 * t239 + t243 * t116;
t285 = t241 * t286;
t218 = qJ(3) * t301;
t231 = t241 * qJD(3);
t304 = qJD(1) * t231 + qJD(2) * t218;
t133 = qJD(1) * t285 + t304;
t80 = pkin(4) * t116 + t133;
t15 = pkin(5) * t56 - pkin(10) * t55 + t80;
t44 = pkin(10) * t230 + t47;
t114 = -pkin(4) * t163 + t146;
t57 = pkin(5) * t113 - pkin(10) * t269 + t114;
t21 = -t238 * t44 + t242 * t57;
t2 = qJD(6) * t21 + t12 * t242 + t15 * t238;
t337 = t2 * t242;
t264 = -t13 * mrSges(6,1) - t12 * mrSges(6,2) + mrSges(7,3) * t337 + Ifges(6,5) * t55 - Ifges(6,6) * t56;
t348 = t165 / 0.2e1;
t349 = -t165 / 0.2e1;
t390 = t75 * mrSges(5,1) - t74 * mrSges(5,2) + Ifges(5,5) * t115 - Ifges(5,6) * t116 + t101 * t348 + t264 + t367 * (Ifges(5,5) * t163 - Ifges(5,6) * t165) / 0.2e1 - t146 * (mrSges(5,1) * t165 + mrSges(5,2) * t163) - (-Ifges(5,1) * t163 + t155) * t349 - (-Ifges(5,2) * t165 + t102 + t154) * t163 / 0.2e1;
t207 = Ifges(7,5) * t238 + Ifges(7,6) * t242;
t208 = Ifges(7,2) * t242 + t327;
t209 = Ifges(7,1) * t238 + t326;
t282 = mrSges(7,1) * t242 - mrSges(7,2) * t238;
t389 = -t13 * t282 + t377 * t400 + t378 * t384 + t56 * t207 / 0.2e1 + t36 * t208 / 0.2e1 + t35 * t209 / 0.2e1 + t401 * qJD(6);
t22 = t238 * t57 + t242 * t44;
t316 = t21 * t242;
t275 = t22 * t238 + t316;
t380 = mrSges(7,3) * t275;
t181 = t239 * t244 + t240 * t243;
t366 = qJD(4) + qJD(5);
t371 = (-qJD(2) + t366) * t181;
t103 = Ifges(6,4) * t113;
t358 = t94 / 0.2e1;
t330 = Ifges(6,1) * t269;
t364 = -t103 / 0.2e1 + t330 / 0.2e1;
t249 = t114 * mrSges(6,2) + Ifges(6,5) * t230 + t280 * t358 + t93 * t384 + t364 + t397 + t401;
t370 = -t249 + t380 + t103 / 0.2e1;
t369 = -t21 * t238 + t22 * t242;
t3 = -qJD(6) * t22 - t12 * t238 + t15 * t242;
t365 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t35 + Ifges(7,6) * t36;
t8 = Ifges(7,4) * t35 + Ifges(7,2) * t36 + Ifges(7,6) * t56;
t9 = t35 * Ifges(7,1) + t36 * Ifges(7,4) + t56 * Ifges(7,5);
t361 = t238 * t9 / 0.2e1 + t8 * t346 + t389 + t397 * qJD(6);
t360 = -t93 / 0.2e1;
t359 = -t94 / 0.2e1;
t355 = pkin(1) * mrSges(3,1);
t354 = pkin(1) * mrSges(3,2);
t352 = -t106 / 0.2e1;
t350 = t163 / 0.2e1;
t342 = pkin(4) * t165;
t136 = -t211 * t240 + t244 * t373;
t104 = pkin(9) * t267 + t136;
t105 = -pkin(9) * t266 + t137;
t271 = t243 * t104 - t105 * t239;
t339 = t13 * t271;
t336 = t3 * t238;
t335 = t46 * mrSges(6,3);
t334 = t47 * mrSges(6,3);
t322 = Ifges(6,2) * t113;
t178 = t239 * t240 - t243 * t244;
t318 = t13 * t178;
t300 = qJD(2) * t245;
t303 = qJ(3) * t300 + t231;
t297 = qJD(6) * t238;
t296 = qJD(6) * t242;
t200 = -t245 * pkin(2) - t241 * qJ(3) - pkin(1);
t294 = t241 * t319;
t293 = 0.3e1 / 0.2e1 * Ifges(3,4) - 0.3e1 / 0.2e1 * Ifges(4,5);
t292 = Ifges(4,4) / 0.2e1 + Ifges(3,5) / 0.2e1;
t291 = Ifges(4,6) / 0.2e1 - Ifges(3,6) / 0.2e1;
t290 = m(4) * pkin(7) + mrSges(4,2);
t175 = t245 * pkin(3) - t200;
t16 = mrSges(7,1) * t56 - mrSges(7,3) * t35;
t17 = -mrSges(7,2) * t56 + mrSges(7,3) * t36;
t276 = -t238 * t16 + t242 * t17;
t121 = -t239 * t266 - t243 * t267;
t132 = pkin(4) * t266 + t175;
t268 = t239 * t267 - t243 * t266;
t65 = -pkin(5) * t268 - pkin(10) * t121 + t132;
t73 = t104 * t239 + t105 * t243;
t29 = -t238 * t73 + t242 * t65;
t30 = t238 * t65 + t242 * t73;
t273 = t100 * t165 + t163 * t99;
t139 = mrSges(5,2) * t367 + t163 * mrSges(5,3);
t140 = -mrSges(5,1) * t367 - t165 * mrSges(5,3);
t270 = t139 * t244 - t140 * t240;
t126 = t191 * t243 - t198 * t239;
t153 = t246 * t302 + t218;
t66 = -mrSges(7,2) * t106 + mrSges(7,3) * t93;
t67 = mrSges(7,1) * t106 - mrSges(7,3) * t94;
t96 = -mrSges(6,2) * t230 - mrSges(6,3) * t113;
t263 = -t238 * t67 + t242 * t66 + t96;
t141 = t285 + t303;
t194 = qJD(2) * t373;
t85 = -t244 * t194 - t211 * t299 + t298 * t373 + t170;
t118 = t153 - t342;
t89 = pkin(4) * t130 + t141;
t257 = -qJD(6) * t275 - t336;
t86 = -qJD(4) * t137 + t194 * t240 + t244 * t195;
t256 = t257 + t337;
t255 = -pkin(9) * t131 + t86;
t252 = t22 * mrSges(7,2) - Ifges(7,3) * t106 - Ifges(7,6) * t93 - Ifges(7,5) * t94 + Ifges(6,6) * t230 - t322 / 0.2e1 + Ifges(6,4) * t269 - t114 * mrSges(6,1) - t21 * mrSges(7,1);
t251 = -t66 * t297 - t67 * t296 + m(7) * (-t21 * t296 - t22 * t297 - t336 + t337) + t276;
t250 = t322 / 0.2e1 - t252;
t248 = t249 + t364;
t196 = (qJD(3) - t226) * qJD(2);
t185 = (-t245 * mrSges(4,1) - mrSges(4,3) * t241) * qJD(1);
t164 = t178 * qJD(2);
t156 = t294 - t303;
t145 = qJD(1) * t294 - t304;
t135 = -t164 * t242 + t238 * t302;
t134 = t164 * t238 + t242 * t302;
t128 = t366 * t178;
t123 = -pkin(10) + t127;
t122 = pkin(5) - t126;
t117 = -mrSges(5,1) * t163 + mrSges(5,2) * t165;
t78 = qJD(5) * t126 + t150 * t243 + t151 * t239;
t76 = mrSges(6,1) * t113 + mrSges(6,2) * t269;
t69 = -pkin(9) * t130 + t85;
t64 = t342 + t77;
t63 = qJD(5) * t121 + t243 * t130 + t131 * t239;
t62 = qJD(5) * t268 - t130 * t239 + t131 * t243;
t61 = t239 * t262 + t243 * t95;
t58 = t118 - t77;
t52 = Ifges(7,3) * t56;
t49 = t243 * t87 - t314;
t48 = t239 * t87 + t313;
t28 = t238 * t77 + t242 * t46;
t27 = -t238 * t46 + t242 * t77;
t26 = t238 * t58 + t242 * t61;
t25 = -t238 * t61 + t242 * t58;
t24 = t238 * t64 + t242 * t49;
t23 = -t238 * t49 + t242 * t64;
t20 = pkin(5) * t63 - pkin(10) * t62 + t89;
t19 = qJD(5) * t73 + t239 * t69 - t243 * t255;
t18 = qJD(5) * t271 + t239 * t255 + t243 * t69;
t5 = -qJD(6) * t30 - t18 * t238 + t20 * t242;
t4 = qJD(6) * t29 + t18 * t242 + t20 * t238;
t1 = [(Ifges(6,1) * t55 + t80 * mrSges(6,2) + t8 * t347 + t35 * t400 + t36 * t384 + t9 * t346 + (mrSges(6,3) + t281) * t13 + (-t2 * t238 - t242 * t3) * mrSges(7,3) + (t42 * t347 - t242 * t41 / 0.2e1 + t43 * t282 + t208 * t360 + t209 * t359 + t207 * t352 - t369 * mrSges(7,3)) * qJD(6) + (-Ifges(6,4) + t402) * t56) * t121 - t367 * (Ifges(5,5) * t131 - Ifges(5,6) * t130) / 0.2e1 + (t116 * t266 - t130 * t350) * Ifges(5,2) + (-t100 * t130 - t115 * t136 - t116 * t137 - t131 * t99 - t266 * t74 + t267 * t75) * mrSges(5,3) + (-t115 * t266 + t116 * t267 - t130 * t348 + t131 * t350) * Ifges(5,4) + t175 * (mrSges(5,1) * t116 + mrSges(5,2) * t115) + m(4) * (t145 * t200 + t156 * t167) + m(5) * (t100 * t85 + t133 * t175 + t136 * t75 + t137 * t74 + t141 * t146 + t86 * t99) + (-t145 * mrSges(4,1) + t290 * t196 + (t292 * qJD(2) + (-0.2e1 * t354 - t200 * mrSges(4,3) + t293 * t245 + (-0.3e1 / 0.2e1 * Ifges(3,2) - 0.3e1 / 0.2e1 * Ifges(4,3) + 0.3e1 / 0.2e1 * Ifges(3,1) + 0.3e1 / 0.2e1 * Ifges(4,1) + t290 * pkin(7)) * t241) * qJD(1) + t391) * qJD(2)) * t245 + (-t145 * mrSges(4,3) + (t291 * qJD(2) + (t200 * mrSges(4,1) - t241 * t293 - 0.2e1 * t355) * qJD(1) + t392) * qJD(2)) * t241 + t250 * t63 + m(7) * (t19 * t43 + t2 * t30 + t21 * t5 + t22 * t4 + t29 * t3 - t339) + m(6) * (t114 * t89 + t12 * t73 + t132 * t80 + t18 * t47 - t19 * t46 - t339) + t29 * t16 + t30 * t17 + (t248 - t380) * t62 + t146 * (mrSges(5,1) * t130 + mrSges(5,2) * t131) + (-t271 * t55 - t46 * t62 - t47 * t63 - t56 * t73) * mrSges(6,3) - t271 * t14 - (-t12 * mrSges(6,3) + t52 / 0.2e1 - Ifges(6,4) * t55 + t80 * mrSges(6,1) + (Ifges(6,2) + Ifges(7,3) / 0.2e1) * t56 + t365) * t268 + (-t115 * t267 + t131 * t348) * Ifges(5,1) + t133 * (mrSges(5,1) * t266 - mrSges(5,2) * t267) + t4 * t66 + t5 * t67 + t89 * t76 + t18 * t96 + t333 * t19 - t130 * t101 / 0.2e1 + t131 * t102 / 0.2e1 + t132 * (mrSges(6,1) * t56 + mrSges(6,2) * t55) + t85 * t139 + t86 * t140 + t141 * t117 + t156 * t185; ((t385 + (t354 + t325 / 0.2e1) * qJD(1) + (-pkin(2) * mrSges(4,2) + (-m(4) * pkin(2) + t398) * pkin(7) + t292) * qJD(2) - t391) * t245 + (t386 + (t355 + t329 / 0.2e1) * qJD(1) + (t387 + t388 - Ifges(4,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t301 + (pkin(7) * mrSges(3,2) - qJ(3) * mrSges(4,2) + t291) * qJD(2) - t392) * t241) * qJD(1) - t390 - t389 + (-t115 * t197 - t116 * t198 - t273) * mrSges(5,3) + (-t114 * t118 + t12 * t127 - t126 * t13 + (-t61 + t78) * t47 + t396 * t46) * m(6) + (t122 * t13 + t123 * t256 - t21 * t25 - t22 * t26 + t369 * t78 - t396 * t43) * m(7) - t333 * t396 + (-t8 / 0.2e1 + t123 * t17) * t242 + ((t21 * mrSges(7,3) - t123 * t67 - t42 / 0.2e1) * t242 + (t41 / 0.2e1 + t22 * mrSges(7,3) - t123 * t66) * t238) * qJD(6) + m(4) * (qJ(3) * t196 + qJD(3) * t203) + (t3 * mrSges(7,3) - t123 * t16 - t9 / 0.2e1) * t238 + t393 * t139 + t394 * t140 + (t100 * t393 - t146 * t153 + t197 * t75 + t198 * t74 + t394 * t99) * m(5) + (-t330 / 0.2e1 + t370) * t113 + (-m(4) * t167 - t185) * (pkin(2) * t302 - t218) + t250 * t269 + (t113 * t46 - t126 * t55 - t127 * t56 - t269 * t47) * mrSges(6,3) - t26 * t66 - t25 * t67 - t61 * t96 - t118 * t76 + t122 * t14 + t263 * t78 - t153 * t117 + t196 * mrSges(4,3) + qJD(3) * t205; -t134 * t67 - t135 * t66 + t164 * t96 + (t55 * mrSges(6,3) + t14) * t178 + t270 * qJD(4) + (-t115 * t244 - t116 * t240) * mrSges(5,3) - t263 * t128 + (-t205 - t270) * qJD(2) + (-t56 * mrSges(6,3) + (-t238 * t66 - t242 * t67) * qJD(6) + t276) * t181 + (t290 * t300 + (-t117 + t185 - t76) * t241) * qJD(1) - m(4) * (qJD(2) * t203 - t167 * t302) + t371 * t333 + (-t128 * t369 - t134 * t21 - t135 * t22 + t256 * t181 + t371 * t43 + t318) * m(7) + (-t114 * t302 + t12 * t181 + t318 + (-t128 + t164) * t47 - t371 * t46) * m(6) + (-t146 * t302 + t240 * t74 + t244 * t75 - t367 * (t100 * t244 - t240 * t99)) * m(5); t390 - m(6) * (-t46 * t48 + t47 * t49) + t273 * mrSges(5,3) + t361 + t251 * (pkin(4) * t239 + pkin(10)) - m(7) * (t21 * t23 + t22 * t24 + t43 * t48) + (-t250 + t334) * t269 + (-t165 * t76 + (-t239 * t56 - t243 * t55) * mrSges(6,3) + (t333 * t239 + t263 * t243 + m(7) * (t239 * t43 + t243 * t369)) * qJD(5) + (0.2e1 * t114 * t349 + t12 * t239 - t13 * t243 + (-t239 * t46 + t243 * t47) * qJD(5)) * m(6)) * pkin(4) + t395 * (-pkin(4) * t243 - pkin(5)) - t24 * t66 - t23 * t67 - t49 * t96 + (-t106 * t316 + (-t106 * t22 - t3) * t238) * mrSges(7,3) - t333 * t48 - (-t248 + t335) * t113 - t99 * t139 + t100 * t140; -m(7) * (t21 * t27 + t22 * t28 + t43 * t47) - t395 * pkin(5) + t361 + t264 + t251 * pkin(10) + (t252 + t334) * t269 - (t335 + (Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t269 + t370) * t113 + t257 * mrSges(7,3) - t28 * t66 - t27 * t67 - t46 * t96 - t333 * t47; t52 - t43 * (mrSges(7,1) * t94 + mrSges(7,2) * t93) + (Ifges(7,1) * t93 - t345) * t359 + t41 * t358 + (Ifges(7,5) * t93 - Ifges(7,6) * t94) * t352 - t21 * t66 + t22 * t67 + (t21 * t93 + t22 * t94) * mrSges(7,3) + (-Ifges(7,2) * t94 + t42 + t91) * t360 + t365;];
tauc  = t1(:);
