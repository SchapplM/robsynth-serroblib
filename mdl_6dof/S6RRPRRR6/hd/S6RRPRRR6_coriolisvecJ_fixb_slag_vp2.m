% Calculate vector of centrifugal and Coriolis load on the joints for
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

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
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR6_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:50:32
% EndTime: 2019-03-09 13:50:50
% DurationCPUTime: 9.45s
% Computational Cost: add. (12959->572), mult. (30237->751), div. (0->0), fcn. (20686->8), ass. (0->263)
t238 = sin(qJ(6));
t349 = -t238 / 0.2e1;
t242 = cos(qJ(6));
t348 = t242 / 0.2e1;
t403 = Ifges(7,5) * t348 + Ifges(7,6) * t349;
t240 = sin(qJ(4));
t244 = cos(qJ(4));
t245 = cos(qJ(2));
t303 = qJD(1) * t245;
t241 = sin(qJ(2));
t304 = qJD(1) * t241;
t163 = -t240 * t304 - t244 * t303;
t165 = -t240 * t303 + t244 * t304;
t239 = sin(qJ(5));
t243 = cos(qJ(5));
t113 = -t243 * t163 + t239 * t165;
t106 = qJD(6) + t113;
t281 = mrSges(7,1) * t238 + mrSges(7,2) * t242;
t367 = qJD(2) - qJD(4);
t230 = qJD(5) - t367;
t226 = pkin(7) * t304;
t192 = pkin(8) * t304 - t226;
t265 = qJD(3) - t192;
t246 = -pkin(2) - pkin(3);
t286 = t246 * qJD(2);
t258 = t286 + t265;
t227 = pkin(7) * t303;
t193 = -pkin(8) * t303 + t227;
t236 = qJD(2) * qJ(3);
t166 = t193 + t236;
t307 = t244 * t166;
t100 = t240 * t258 + t307;
t340 = t163 * pkin(9);
t88 = t100 + t340;
t316 = t239 * t88;
t343 = pkin(9) * t165;
t142 = t244 * t258;
t99 = -t166 * t240 + t142;
t87 = t99 - t343;
t83 = -pkin(4) * t367 + t87;
t46 = t243 * t83 - t316;
t43 = -pkin(5) * t230 - t46;
t402 = t106 * t403 + t43 * t281;
t328 = Ifges(7,4) * t242;
t385 = Ifges(7,2) * t349 + t328 / 0.2e1;
t329 = Ifges(7,4) * t238;
t280 = Ifges(7,1) * t242 - t329;
t401 = t280 / 0.2e1;
t399 = -mrSges(3,1) - mrSges(4,1);
t269 = t163 * t239 + t243 * t165;
t94 = t230 * t238 + t242 * t269;
t347 = Ifges(7,4) * t94;
t93 = t230 * t242 - t238 * t269;
t41 = Ifges(7,2) * t93 + Ifges(7,6) * t106 + t347;
t91 = Ifges(7,4) * t93;
t42 = t94 * Ifges(7,1) + t106 * Ifges(7,5) + t91;
t398 = t42 * t348 + t41 * t349;
t77 = pkin(5) * t269 + pkin(10) * t113;
t397 = -t303 / 0.2e1;
t335 = -mrSges(6,1) * t230 - mrSges(7,1) * t93 + mrSges(7,2) * t94 + mrSges(6,3) * t269;
t197 = -qJ(3) * t240 + t244 * t246;
t191 = -pkin(4) + t197;
t198 = t244 * qJ(3) + t240 * t246;
t127 = t239 * t191 + t243 * t198;
t150 = t244 * qJD(3) + qJD(4) * t197;
t151 = -t240 * qJD(3) - qJD(4) * t198;
t124 = -t192 * t240 + t244 * t193;
t262 = t124 + t340;
t125 = t244 * t192 + t240 * t193;
t95 = t125 + t343;
t396 = qJD(5) * t127 + (-t151 + t262) * t243 + (t150 - t95) * t239;
t266 = t240 * t241 + t244 * t245;
t131 = t367 * t266;
t115 = t131 * qJD(1);
t342 = t115 * pkin(9);
t315 = t243 * t88;
t47 = t239 * t83 + t315;
t267 = t240 * t245 - t241 * t244;
t130 = t367 * t267;
t116 = t130 * qJD(1);
t358 = pkin(7) - pkin(8);
t211 = t358 * t245;
t195 = qJD(2) * t211;
t170 = t240 * t195;
t301 = qJD(4) * t240;
t373 = t358 * t241;
t74 = qJD(4) * t142 - t166 * t301 + t244 * (-qJD(1) * t373 + qJD(3)) * qJD(2) + qJD(1) * t170;
t51 = -pkin(9) * t116 + t74;
t137 = t244 * t211 + t240 * t373;
t254 = ((-qJD(4) * t246 - qJD(3)) * t240 + t137 * qJD(1)) * qJD(2);
t261 = t240 * t265;
t75 = (-t261 - t307) * qJD(4) + t254;
t13 = t239 * t51 - t243 * (t75 - t342) + t47 * qJD(5);
t378 = qJD(6) * t93;
t55 = -qJD(5) * t113 + t115 * t243 - t116 * t239;
t35 = t242 * t55 + t378;
t377 = qJD(6) * t94;
t36 = -t238 * t55 - t377;
t14 = -mrSges(7,1) * t36 + mrSges(7,2) * t35;
t395 = m(7) * t13 + t14;
t394 = -t124 + t151;
t393 = -t125 + t150;
t167 = -qJD(1) * pkin(1) - pkin(2) * t303 - qJ(3) * t304;
t321 = qJD(2) * pkin(2);
t199 = qJD(3) + t226 - t321;
t327 = Ifges(4,5) * t245;
t382 = -qJD(2) / 0.2e1;
t383 = -qJD(1) / 0.2e1;
t386 = Ifges(3,4) * t397;
t388 = -Ifges(3,1) / 0.2e1;
t392 = (m(4) * t199 + (mrSges(4,2) + mrSges(3,3)) * t304 + t399 * qJD(2)) * pkin(7) - t167 * mrSges(4,3) - (Ifges(4,1) * t241 - t327) * t383 - t304 * t388 - t386 + t199 * mrSges(4,2) - (Ifges(4,4) + Ifges(3,5)) * t382;
t203 = t227 + t236;
t205 = mrSges(4,2) * t303 + qJD(2) * mrSges(4,3);
t223 = Ifges(4,5) * t304;
t331 = Ifges(3,4) * t241;
t387 = Ifges(4,6) / 0.2e1;
t391 = -(m(4) * t203 - qJD(2) * mrSges(3,2) + mrSges(3,3) * t303 + t205) * pkin(7) + t167 * mrSges(4,1) + qJD(2) * t387 + Ifges(4,3) * t397 + t223 / 0.2e1 + Ifges(3,6) * t382 + (t245 * Ifges(3,2) + t331) * t383 - t203 * mrSges(4,2);
t155 = Ifges(5,4) * t165;
t101 = t163 * Ifges(5,2) - Ifges(5,6) * t367 + t155;
t154 = Ifges(5,4) * t163;
t102 = t165 * Ifges(5,1) - Ifges(5,5) * t367 + t154;
t146 = pkin(3) * t303 - t167;
t300 = qJD(4) * t244;
t12 = t46 * qJD(5) + t243 * t51 + (-qJD(4) * t261 - t166 * t300 + t254 - t342) * t239;
t56 = qJD(5) * t269 + t115 * t239 + t243 * t116;
t285 = t241 * t286;
t218 = qJ(3) * t303;
t231 = t241 * qJD(3);
t306 = qJD(1) * t231 + qJD(2) * t218;
t133 = qJD(1) * t285 + t306;
t80 = pkin(4) * t116 + t133;
t15 = pkin(5) * t56 - pkin(10) * t55 + t80;
t44 = pkin(10) * t230 + t47;
t114 = -pkin(4) * t163 + t146;
t57 = pkin(5) * t113 - pkin(10) * t269 + t114;
t21 = -t238 * t44 + t242 * t57;
t2 = qJD(6) * t21 + t12 * t242 + t15 * t238;
t339 = t2 * t242;
t264 = -t13 * mrSges(6,1) - t12 * mrSges(6,2) + mrSges(7,3) * t339 + Ifges(6,5) * t55 - Ifges(6,6) * t56;
t350 = t165 / 0.2e1;
t351 = -t165 / 0.2e1;
t390 = t75 * mrSges(5,1) - t74 * mrSges(5,2) + Ifges(5,5) * t115 - Ifges(5,6) * t116 + t101 * t350 + t264 + t367 * (Ifges(5,5) * t163 - Ifges(5,6) * t165) / 0.2e1 - t146 * (mrSges(5,1) * t165 + mrSges(5,2) * t163) - (-Ifges(5,1) * t163 + t155) * t351 - (-Ifges(5,2) * t165 + t102 + t154) * t163 / 0.2e1;
t207 = Ifges(7,5) * t238 + Ifges(7,6) * t242;
t208 = Ifges(7,2) * t242 + t329;
t209 = Ifges(7,1) * t238 + t328;
t282 = mrSges(7,1) * t242 - mrSges(7,2) * t238;
t389 = -t13 * t282 + t377 * t401 + t378 * t385 + t56 * t207 / 0.2e1 + t36 * t208 / 0.2e1 + t35 * t209 / 0.2e1 + t402 * qJD(6);
t22 = t238 * t57 + t242 * t44;
t318 = t21 * t242;
t275 = t22 * t238 + t318;
t380 = mrSges(7,3) * t275;
t181 = t239 * t244 + t240 * t243;
t366 = qJD(4) + qJD(5);
t371 = (-qJD(2) + t366) * t181;
t103 = Ifges(6,4) * t113;
t289 = t103 / 0.2e1;
t332 = Ifges(6,1) * t269;
t370 = t289 - t332 / 0.2e1;
t369 = -t21 * t238 + t22 * t242;
t3 = -qJD(6) * t22 - t12 * t238 + t15 * t242;
t365 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t35 + Ifges(7,6) * t36;
t8 = t35 * Ifges(7,4) + t36 * Ifges(7,2) + t56 * Ifges(7,6);
t9 = t35 * Ifges(7,1) + t36 * Ifges(7,4) + t56 * Ifges(7,5);
t362 = t238 * t9 / 0.2e1 + t8 * t348 + t389 + t398 * qJD(6);
t361 = -t93 / 0.2e1;
t360 = -t94 / 0.2e1;
t359 = t94 / 0.2e1;
t356 = pkin(1) * mrSges(3,1);
t355 = pkin(1) * mrSges(3,2);
t354 = -t106 / 0.2e1;
t352 = t163 / 0.2e1;
t344 = pkin(4) * t165;
t136 = -t211 * t240 + t244 * t373;
t104 = pkin(9) * t267 + t136;
t105 = -pkin(9) * t266 + t137;
t271 = t243 * t104 - t105 * t239;
t341 = t13 * t271;
t338 = t3 * t238;
t337 = t46 * mrSges(6,3);
t336 = t47 * mrSges(6,3);
t324 = Ifges(6,2) * t113;
t178 = t239 * t240 - t243 * t244;
t320 = t13 * t178;
t302 = qJD(2) * t245;
t305 = qJ(3) * t302 + t231;
t299 = qJD(6) * t238;
t298 = qJD(6) * t242;
t200 = -t245 * pkin(2) - t241 * qJ(3) - pkin(1);
t296 = t241 * t321;
t295 = Ifges(4,4) / 0.2e1 + Ifges(3,5) / 0.2e1;
t294 = 0.3e1 / 0.2e1 * Ifges(4,5) - 0.3e1 / 0.2e1 * Ifges(3,4);
t293 = t387 - Ifges(3,6) / 0.2e1;
t292 = m(4) * pkin(7) + mrSges(4,2);
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
t140 = -mrSges(5,1) * t367 - mrSges(5,3) * t165;
t270 = t139 * t244 - t140 * t240;
t126 = t191 * t243 - t198 * t239;
t153 = t246 * t304 + t218;
t66 = -mrSges(7,2) * t106 + mrSges(7,3) * t93;
t67 = mrSges(7,1) * t106 - mrSges(7,3) * t94;
t96 = -mrSges(6,2) * t230 - mrSges(6,3) * t113;
t263 = -t238 * t67 + t242 * t66 + t96;
t141 = t285 + t305;
t194 = qJD(2) * t373;
t85 = -t244 * t194 - t211 * t301 + t300 * t373 + t170;
t118 = t153 - t344;
t89 = pkin(4) * t130 + t141;
t257 = -qJD(6) * t275 - t338;
t86 = -qJD(4) * t137 + t194 * t240 + t244 * t195;
t256 = t257 + t339;
t255 = -pkin(9) * t131 + t86;
t252 = t22 * mrSges(7,2) - Ifges(7,3) * t106 - Ifges(7,6) * t93 - Ifges(7,5) * t94 + Ifges(6,6) * t230 - t324 / 0.2e1 + Ifges(6,4) * t269 - t114 * mrSges(6,1) - t21 * mrSges(7,1);
t251 = -t66 * t299 - t67 * t298 + m(7) * (-t21 * t298 - t22 * t299 - t338 + t339) + t276;
t250 = t324 / 0.2e1 - t252;
t249 = t114 * mrSges(6,2) - t103 / 0.2e1 + Ifges(6,5) * t230 + t332 / 0.2e1 + t93 * t385 + t280 * t359 + t398 + t402;
t248 = -t249 + t370;
t247 = -t249 + t380;
t196 = (qJD(3) - t226) * qJD(2);
t185 = (-t245 * mrSges(4,1) - mrSges(4,3) * t241) * qJD(1);
t164 = t178 * qJD(2);
t156 = t296 - t305;
t145 = qJD(1) * t296 - t306;
t135 = -t164 * t242 + t238 * t304;
t134 = t164 * t238 + t242 * t304;
t128 = t366 * t178;
t123 = -pkin(10) + t127;
t122 = pkin(5) - t126;
t117 = -mrSges(5,1) * t163 + mrSges(5,2) * t165;
t78 = qJD(5) * t126 + t150 * t243 + t151 * t239;
t76 = mrSges(6,1) * t113 + mrSges(6,2) * t269;
t69 = -pkin(9) * t130 + t85;
t64 = t344 + t77;
t63 = qJD(5) * t121 + t243 * t130 + t131 * t239;
t62 = qJD(5) * t268 - t130 * t239 + t131 * t243;
t61 = t239 * t262 + t243 * t95;
t58 = t118 - t77;
t52 = Ifges(7,3) * t56;
t49 = t243 * t87 - t316;
t48 = t239 * t87 + t315;
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
t1 = [(-t145 * mrSges(4,3) + (t293 * qJD(2) + (t200 * mrSges(4,1) + t241 * t294 - 0.2e1 * t356) * qJD(1) + t391) * qJD(2)) * t241 + (-t145 * mrSges(4,1) + t292 * t196 + (t295 * qJD(2) + (-0.2e1 * t355 - t200 * mrSges(4,3) - t294 * t245 + (0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(4,3) - 0.3e1 / 0.2e1 * Ifges(3,2) + t292 * pkin(7)) * t241) * qJD(1) + t392) * qJD(2)) * t245 + (Ifges(6,1) * t55 + t80 * mrSges(6,2) + t8 * t349 + t9 * t348 + t35 * t401 + t36 * t385 + (mrSges(6,3) + t281) * t13 + (-t2 * t238 - t242 * t3) * mrSges(7,3) + (t42 * t349 - t242 * t41 / 0.2e1 + t207 * t354 + t208 * t361 + t209 * t360 + t43 * t282 - t369 * mrSges(7,3)) * qJD(6) + (-Ifges(6,4) + t403) * t56) * t121 - t367 * (Ifges(5,5) * t131 - Ifges(5,6) * t130) / 0.2e1 + (t116 * t266 - t130 * t352) * Ifges(5,2) + (-t115 * t266 + t116 * t267 - t130 * t350 + t131 * t352) * Ifges(5,4) + (-t100 * t130 - t115 * t136 - t116 * t137 - t131 * t99 - t266 * t74 + t267 * t75) * mrSges(5,3) + t175 * (mrSges(5,1) * t116 + mrSges(5,2) * t115) + (-t248 - t380) * t62 + t18 * t96 + t89 * t76 + t146 * (mrSges(5,1) * t130 + mrSges(5,2) * t131) - t271 * t14 + (-t271 * t55 - t46 * t62 - t47 * t63 - t56 * t73) * mrSges(6,3) - (t52 / 0.2e1 - Ifges(6,4) * t55 + t80 * mrSges(6,1) - t12 * mrSges(6,3) + (Ifges(6,2) + Ifges(7,3) / 0.2e1) * t56 + t365) * t268 + t133 * (mrSges(5,1) * t266 - mrSges(5,2) * t267) + (-t115 * t267 + t131 * t350) * Ifges(5,1) + t4 * t66 + t5 * t67 + t29 * t16 + t30 * t17 + t250 * t63 + m(5) * (t100 * t85 + t133 * t175 + t136 * t75 + t137 * t74 + t141 * t146 + t86 * t99) + m(4) * (t145 * t200 + t156 * t167) + m(7) * (t19 * t43 + t2 * t30 + t21 * t5 + t22 * t4 + t29 * t3 - t341) + m(6) * (t114 * t89 + t12 * t73 + t132 * t80 + t18 * t47 - t19 * t46 - t341) - t130 * t101 / 0.2e1 + t131 * t102 / 0.2e1 + t132 * (mrSges(6,1) * t56 + mrSges(6,2) * t55) + t85 * t139 + t86 * t140 + t141 * t117 + t156 * t185 + t335 * t19; t335 * t396 + (t122 * t13 + t123 * t256 - t21 * t25 - t22 * t26 + t369 * t78 + t396 * t43) * m(7) + (-t114 * t118 + t12 * t127 - t126 * t13 + (-t61 + t78) * t47 - t396 * t46) * m(6) + t393 * t139 + (t100 * t393 - t146 * t153 + t197 * t75 + t198 * t74 + t394 * t99) * m(5) + t394 * t140 + (-t115 * t197 - t116 * t198 - t273) * mrSges(5,3) + (t123 * t17 - t8 / 0.2e1) * t242 + ((t386 + (t355 + t327 / 0.2e1) * qJD(1) + (-pkin(2) * mrSges(4,2) + (-m(4) * pkin(2) + t399) * pkin(7) + t295) * qJD(2) - t392) * t245 + (-t223 / 0.2e1 + (t356 + t331 / 0.2e1) * qJD(1) + (Ifges(4,3) / 0.2e1 + t388 - Ifges(4,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t303 + (pkin(7) * mrSges(3,2) - qJ(3) * mrSges(4,2) + t293) * qJD(2) - t391) * t241) * qJD(1) + (t247 + t370) * t113 - t61 * t96 + (-m(4) * t167 - t185) * (pkin(2) * t304 - t218) + (t113 * t46 - t126 * t55 - t127 * t56 - t269 * t47) * mrSges(6,3) + t250 * t269 - t26 * t66 - t25 * t67 + ((t21 * mrSges(7,3) - t123 * t67 - t42 / 0.2e1) * t242 + (t22 * mrSges(7,3) - t123 * t66 + t41 / 0.2e1) * t238) * qJD(6) + (-t9 / 0.2e1 + t3 * mrSges(7,3) - t123 * t16) * t238 + t263 * t78 + m(4) * (qJ(3) * t196 + qJD(3) * t203) - t389 - t118 * t76 + t122 * t14 - t153 * t117 - t390 + t196 * mrSges(4,3) + qJD(3) * t205; -t134 * t67 - t135 * t66 + t164 * t96 + (t55 * mrSges(6,3) + t14) * t178 + t270 * qJD(4) + (-t115 * t244 - t116 * t240) * mrSges(5,3) - t263 * t128 + (-t205 - t270) * qJD(2) + (-t56 * mrSges(6,3) + (-t238 * t66 - t242 * t67) * qJD(6) + t276) * t181 + (t292 * t302 + (-t117 + t185 - t76) * t241) * qJD(1) - m(4) * (qJD(2) * t203 - t167 * t304) + t371 * t335 + (-t128 * t369 - t134 * t21 - t135 * t22 + t181 * t256 + t371 * t43 + t320) * m(7) + (-t114 * t304 + t12 * t181 + t320 + (-t128 + t164) * t47 - t371 * t46) * m(6) + (-t146 * t304 + t240 * t74 + t244 * t75 - t367 * (t100 * t244 - t240 * t99)) * m(5); (-t165 * t76 + (-t239 * t56 - t243 * t55) * mrSges(6,3) + (t335 * t239 + t263 * t243 + m(7) * (t239 * t43 + t243 * t369)) * qJD(5) + (0.2e1 * t114 * t351 + t12 * t239 - t13 * t243 + (-t239 * t46 + t243 * t47) * qJD(5)) * m(6)) * pkin(4) - t49 * t96 + t395 * (-pkin(4) * t243 - pkin(5)) - t24 * t66 - t23 * t67 + (-t250 + t336) * t269 + (-t106 * t318 + (-t106 * t22 - t3) * t238) * mrSges(7,3) - m(6) * (-t46 * t48 + t47 * t49) + t251 * (pkin(4) * t239 + pkin(10)) - m(7) * (t21 * t23 + t22 * t24 + t43 * t48) + t273 * mrSges(5,3) + t362 - t99 * t139 + t100 * t140 + t390 - t335 * t48 - (t248 + t337) * t113; -t46 * t96 - t28 * t66 - t27 * t67 - (t289 + t337 + (Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t269 + t247) * t113 + (t252 + t336) * t269 - m(7) * (t21 * t27 + t22 * t28 + t43 * t47) + t257 * mrSges(7,3) + t251 * pkin(10) + t264 - t395 * pkin(5) + t362 - t335 * t47; t52 - t43 * (mrSges(7,1) * t94 + mrSges(7,2) * t93) + (Ifges(7,1) * t93 - t347) * t360 + t41 * t359 + (Ifges(7,5) * t93 - Ifges(7,6) * t94) * t354 - t21 * t66 + t22 * t67 + (t21 * t93 + t22 * t94) * mrSges(7,3) + (-Ifges(7,2) * t94 + t42 + t91) * t361 + t365;];
tauc  = t1(:);
