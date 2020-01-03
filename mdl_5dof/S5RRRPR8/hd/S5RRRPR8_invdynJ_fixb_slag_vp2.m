% Calculate vector of inverse dynamics joint torques for
% S5RRRPR8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR8_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR8_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR8_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR8_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:19:08
% EndTime: 2019-12-31 21:19:33
% DurationCPUTime: 11.84s
% Computational Cost: add. (5461->556), mult. (11926->705), div. (0->0), fcn. (7875->10), ass. (0->267)
t411 = mrSges(4,1) - mrSges(5,2);
t220 = sin(qJ(5));
t224 = cos(qJ(5));
t405 = mrSges(6,1) * t220 + mrSges(6,2) * t224;
t410 = -Ifges(4,5) + Ifges(5,4);
t409 = -Ifges(4,6) + Ifges(5,5);
t221 = sin(qJ(3));
t352 = cos(qJ(3));
t353 = cos(qJ(2));
t267 = t352 * t353;
t222 = sin(qJ(2));
t303 = qJD(1) * t222;
t151 = -qJD(1) * t267 + t221 * t303;
t218 = qJD(2) + qJD(3);
t129 = t151 * t220 + t218 * t224;
t145 = Ifges(4,4) * t151;
t165 = t221 * t353 + t222 * t352;
t152 = t165 * qJD(1);
t147 = qJD(5) + t152;
t396 = t147 * Ifges(6,3);
t128 = t151 * t224 - t218 * t220;
t397 = t128 * Ifges(6,6);
t408 = t152 * Ifges(4,1) + t218 * Ifges(4,5) + t129 * Ifges(6,5) - t145 + t396 + t397;
t164 = t221 * t222 - t267;
t298 = qJD(1) * qJD(2);
t407 = t353 * qJDD(1) - t222 * t298;
t219 = qJ(2) + qJ(3);
t216 = cos(t219);
t406 = t405 * t216;
t215 = sin(t219);
t404 = -t411 * t216 + (mrSges(4,2) - mrSges(5,3)) * t215;
t258 = mrSges(6,1) * t224 - mrSges(6,2) * t220;
t334 = Ifges(6,4) * t129;
t51 = Ifges(6,2) * t128 + Ifges(6,6) * t147 + t334;
t294 = t353 * pkin(6);
t180 = pkin(7) * t353 + t294;
t168 = t180 * qJD(1);
t154 = t352 * t168;
t226 = -pkin(7) - pkin(6);
t179 = t226 * t222;
t167 = qJD(1) * t179;
t160 = qJD(2) * pkin(2) + t167;
t116 = t221 * t160 + t154;
t106 = -t218 * qJ(4) - t116;
t344 = t151 * pkin(4);
t66 = -t106 - t344;
t403 = t258 * t66 - t224 * t51 / 0.2e1;
t310 = t221 * t168;
t115 = -t352 * t160 + t310;
t343 = t152 * pkin(4);
t248 = t115 + t343;
t402 = qJD(4) + t248;
t295 = t353 * pkin(2);
t207 = t295 + pkin(1);
t178 = t207 * qJD(1);
t233 = -t152 * qJ(4) - t178;
t92 = t151 * pkin(3) + t233;
t401 = t178 * mrSges(4,1) + t92 * mrSges(5,2);
t367 = pkin(3) + pkin(8);
t54 = -t218 * t367 + t402;
t56 = t151 * t367 + t233;
t18 = -t220 * t56 + t224 * t54;
t19 = t220 * t54 + t224 * t56;
t400 = t18 * mrSges(6,1) - t178 * mrSges(4,2) - t19 * mrSges(6,2) + t115 * mrSges(4,3) - t92 * mrSges(5,3);
t217 = qJDD(2) + qJDD(3);
t279 = qJD(2) * t353;
t170 = qJD(1) * t279 + t222 * qJDD(1);
t232 = t165 * qJD(3);
t87 = qJD(1) * t232 + t221 * t170 - t352 * t407;
t42 = qJD(5) * t128 + t217 * t224 + t220 * t87;
t370 = t42 / 0.2e1;
t43 = -qJD(5) * t129 - t217 * t220 + t224 * t87;
t369 = t43 / 0.2e1;
t86 = qJD(3) * t151 - t352 * t170 - t221 * t407;
t84 = qJDD(5) - t86;
t368 = t84 / 0.2e1;
t356 = t217 / 0.2e1;
t399 = t407 / 0.2e1;
t10 = -mrSges(6,1) * t43 + mrSges(6,2) * t42;
t69 = mrSges(5,1) * t87 - mrSges(5,3) * t217;
t398 = t10 - t69;
t340 = mrSges(5,1) * t151;
t134 = -mrSges(5,3) * t218 + t340;
t65 = -mrSges(6,1) * t128 + mrSges(6,2) * t129;
t394 = t65 - t134;
t132 = -mrSges(4,2) * t218 - mrSges(4,3) * t151;
t393 = t132 - t134;
t330 = t152 * mrSges(5,1);
t337 = mrSges(4,3) * t152;
t392 = t218 * t411 - t330 - t337;
t225 = cos(qJ(1));
t313 = t216 * t225;
t315 = t215 * t225;
t391 = pkin(3) * t313 + qJ(4) * t315;
t183 = qJ(4) * t313;
t388 = -m(6) * t183 - mrSges(5,2) * t315 - mrSges(5,3) * t313 - t225 * t406;
t223 = sin(qJ(1));
t314 = t216 * t223;
t181 = qJ(4) * t314;
t316 = t215 * t223;
t387 = -m(6) * t181 - mrSges(5,2) * t316 - mrSges(5,3) * t314 - t223 * t406;
t162 = t407 * pkin(6);
t163 = t170 * pkin(6);
t386 = t353 * t162 + t163 * t222;
t280 = qJD(1) * t353;
t349 = pkin(6) * t222;
t385 = (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t280) * t349 + (qJD(2) * mrSges(3,1) - mrSges(3,3) * t303) * t294;
t14 = mrSges(6,1) * t84 - mrSges(6,3) * t42;
t15 = -mrSges(6,2) * t84 + mrSges(6,3) * t43;
t384 = t224 * t14 + t220 * t15;
t383 = g(1) * t225 + g(2) * t223;
t382 = -t115 - qJD(4);
t247 = mrSges(3,1) * t222 + mrSges(3,2) * t353;
t335 = Ifges(3,4) * t222;
t381 = pkin(1) * t247 - t222 * (Ifges(3,1) * t353 - t335) / 0.2e1;
t380 = 0.2e1 * t356;
t379 = -t216 * mrSges(6,3) - t215 * t405 + t404;
t346 = g(3) * t216;
t378 = -t215 * t383 + t346;
t177 = -mrSges(3,1) * t353 + t222 * mrSges(3,2);
t377 = -m(3) * pkin(1) - mrSges(2,1) + t177 + t404;
t93 = -pkin(3) * t218 - t382;
t376 = m(5) * t93 - t392;
t125 = qJDD(2) * pkin(2) - t170 * pkin(7) - t163;
t127 = pkin(7) * t407 + t162;
t278 = qJD(3) * t352;
t301 = qJD(3) * t221;
t30 = t125 * t352 - t221 * t127 - t160 * t301 - t168 * t278;
t235 = qJDD(4) - t30;
t12 = -t86 * pkin(4) - t217 * t367 + t235;
t323 = qJDD(1) * pkin(1);
t143 = -pkin(2) * t407 - t323;
t229 = t86 * qJ(4) - t152 * qJD(4) + t143;
t9 = t367 * t87 + t229;
t1 = qJD(5) * t18 + t12 * t220 + t224 * t9;
t2 = -qJD(5) * t19 + t12 * t224 - t220 * t9;
t375 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t374 = -m(3) * pkin(6) - m(6) * (pkin(4) - t226) - mrSges(5,1) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t253 = t18 * t220 - t19 * t224;
t345 = t1 * t220;
t230 = -qJD(5) * t253 + t2 * t224 + t345;
t299 = qJD(5) * t224;
t300 = qJD(5) * t220;
t90 = -mrSges(6,2) * t147 + mrSges(6,3) * t128;
t91 = mrSges(6,1) * t147 - mrSges(6,3) * t129;
t373 = m(6) * t230 + t90 * t299 - t91 * t300 + t384;
t371 = Ifges(6,1) * t370 + Ifges(6,4) * t369 + Ifges(6,5) * t368;
t365 = -t128 / 0.2e1;
t364 = -t129 / 0.2e1;
t363 = t129 / 0.2e1;
t362 = -t147 / 0.2e1;
t361 = -t151 / 0.2e1;
t360 = t151 / 0.2e1;
t359 = -t152 / 0.2e1;
t358 = t152 / 0.2e1;
t355 = -t218 / 0.2e1;
t354 = t218 / 0.2e1;
t351 = pkin(2) * t221;
t350 = pkin(2) * t222;
t204 = t216 * pkin(3);
t336 = mrSges(6,3) * t224;
t333 = Ifges(6,4) * t220;
t332 = Ifges(6,4) * t224;
t329 = t152 * Ifges(4,4);
t328 = t152 * Ifges(5,6);
t122 = qJD(2) * t165 + t232;
t322 = t122 * t220;
t321 = t122 * t224;
t320 = t152 * t220;
t318 = t164 * t220;
t317 = t164 * t224;
t198 = t215 * qJ(4);
t312 = t220 * t223;
t311 = t220 * t225;
t308 = t223 * t224;
t307 = t224 * t225;
t304 = t204 + t198;
t302 = qJD(2) * t222;
t297 = Ifges(6,5) * t42 + Ifges(6,6) * t43 + Ifges(6,3) * t84;
t293 = t352 * pkin(2);
t211 = pkin(2) * t302;
t210 = pkin(2) * t303;
t288 = Ifges(3,4) * t353;
t276 = -t300 / 0.2e1;
t70 = -t86 * mrSges(5,1) + t217 * mrSges(5,2);
t272 = -t207 - t198;
t118 = t167 * t221 + t154;
t130 = -t352 * t179 + t180 * t221;
t184 = t225 * t207;
t271 = -t223 * t226 + t184;
t269 = pkin(2) * t278;
t268 = -m(6) * t367 - mrSges(6,3);
t206 = -t293 - pkin(3);
t262 = -pkin(3) * t215 - t350;
t259 = mrSges(4,1) * t215 + mrSges(4,2) * t216;
t256 = Ifges(6,1) * t220 + t332;
t255 = Ifges(6,2) * t224 + t333;
t254 = Ifges(6,5) * t220 + Ifges(6,6) * t224;
t110 = pkin(3) * t152 + qJ(4) * t151;
t252 = -t220 * t91 + t224 * t90;
t104 = pkin(4) * t165 + t130;
t241 = -t165 * qJ(4) - t207;
t85 = t164 * t367 + t241;
t33 = t104 * t224 - t220 * t85;
t34 = t104 * t220 + t224 * t85;
t250 = t295 + t304;
t249 = t268 * t215;
t94 = t110 + t210;
t246 = t353 * Ifges(3,2) + t335;
t245 = Ifges(3,5) * t353 - Ifges(3,6) * t222;
t119 = t167 * t352 - t310;
t131 = t221 * t179 + t180 * t352;
t244 = t164 * t299 + t322;
t243 = t164 * t300 - t321;
t121 = t164 * t218;
t242 = qJ(4) * t121 - qJD(4) * t165 + t211;
t29 = t221 * t125 + t352 * t127 + t160 * t278 - t168 * t301;
t169 = qJD(2) * t179;
t234 = qJD(2) * t180;
t57 = -t352 * t169 - t179 * t278 + t180 * t301 + t221 * t234;
t231 = -m(6) * t350 + t249;
t23 = -qJ(4) * t217 - qJD(4) * t218 - t29;
t58 = qJD(3) * t131 + t221 * t169 + t352 * t234;
t144 = Ifges(5,6) * t151;
t100 = t218 * Ifges(5,4) - t152 * Ifges(5,2) + t144;
t101 = -t151 * Ifges(4,2) + t218 * Ifges(4,6) + t329;
t13 = -pkin(4) * t87 - t23;
t27 = -t217 * pkin(3) + t235;
t126 = Ifges(6,4) * t128;
t52 = t129 * Ifges(6,1) + t147 * Ifges(6,5) + t126;
t7 = t42 * Ifges(6,4) + t43 * Ifges(6,2) + t84 * Ifges(6,6);
t99 = t218 * Ifges(5,5) + t151 * Ifges(5,3) - t328;
t228 = (-t19 * t299 - t345 + (t300 + t320) * t18) * mrSges(6,3) + (Ifges(5,1) + Ifges(4,3)) * t217 + (-Ifges(4,2) * t152 - t145 + t408) * t360 + (t328 + t101) * t358 - t106 * t330 + (t144 + t100) * t361 + (-t329 + t99) * t359 + (-t320 / 0.2e1 + t276) * t52 - t2 * t336 + t409 * t87 + (Ifges(5,3) * t361 - t19 * t336 + t254 * t362 + t255 * t365 + t256 * t364 + t355 * t409 + t401 + t403) * t152 + (-Ifges(4,1) * t359 - Ifges(6,5) * t364 + Ifges(5,2) * t358 - Ifges(6,6) * t365 - Ifges(6,3) * t362 + t355 * t410 + t400) * t151 + t410 * t86 - t23 * mrSges(5,3) + t27 * mrSges(5,2) - t29 * mrSges(4,2) + t30 * mrSges(4,1) + t403 * qJD(5) + t13 * t405 + t116 * t337 + t93 * t340 + (Ifges(6,5) * t224 - Ifges(6,6) * t220) * t368 + (-Ifges(6,2) * t220 + t332) * t369 + (Ifges(6,1) * t224 - t333) * t370 + t224 * t371 - t220 * t7 / 0.2e1 - (t128 * t255 + t129 * t256 + t147 * t254) * qJD(5) / 0.2e1;
t209 = Ifges(3,4) * t280;
t203 = t216 * pkin(8);
t202 = qJ(4) + t351;
t193 = t269 + qJD(4);
t150 = Ifges(3,1) * t303 + Ifges(3,5) * qJD(2) + t209;
t149 = Ifges(3,6) * qJD(2) + qJD(1) * t246;
t146 = t152 * pkin(8);
t142 = -t215 * t312 + t307;
t141 = t215 * t308 + t311;
t140 = t215 * t311 + t308;
t139 = t215 * t307 - t312;
t114 = t164 * pkin(3) + t241;
t113 = -mrSges(5,2) * t151 - mrSges(5,3) * t152;
t112 = mrSges(4,1) * t151 + mrSges(4,2) * t152;
t105 = -t164 * pkin(4) + t131;
t89 = t119 - t343;
t88 = t118 - t344;
t79 = t116 - t344;
t75 = t110 + t146;
t68 = -mrSges(4,2) * t217 - mrSges(4,3) * t87;
t67 = mrSges(4,1) * t217 + mrSges(4,3) * t86;
t64 = t146 + t94;
t46 = pkin(3) * t122 + t242;
t39 = -t121 * pkin(4) + t58;
t38 = -pkin(4) * t122 - t57;
t25 = t220 * t79 + t224 * t75;
t24 = -t220 * t75 + t224 * t79;
t22 = t220 * t88 + t224 * t64;
t21 = -t220 * t64 + t224 * t88;
t20 = t122 * t367 + t242;
t17 = t87 * pkin(3) + t229;
t4 = -qJD(5) * t34 - t20 * t220 + t224 * t39;
t3 = qJD(5) * t33 + t20 * t224 + t220 * t39;
t5 = [-t177 * t323 + t353 * (Ifges(3,4) * t170 + Ifges(3,2) * t407) / 0.2e1 - pkin(1) * (-mrSges(3,1) * t407 + t170 * mrSges(3,2)) + (t170 * t349 + t294 * t407 + t386) * mrSges(3,3) + (t1 * t317 - t18 * t244 - t19 * t243 - t2 * t318) * mrSges(6,3) + t66 * (mrSges(6,1) * t243 + mrSges(6,2) * t244) + m(5) * (t114 * t17 + t46 * t92) + (t27 * mrSges(5,1) + t143 * mrSges(4,2) - t30 * mrSges(4,3) - t17 * mrSges(5,3) + Ifges(4,5) * t356 + Ifges(6,5) * t370 + Ifges(6,6) * t369 + Ifges(6,3) * t368 + (-Ifges(4,4) / 0.2e1 - Ifges(5,6)) * t87 - Ifges(5,2) * t86 - t380 * Ifges(5,4) + t375) * t165 + (-Ifges(4,4) * t87 + Ifges(4,5) * t217 + t297) * t165 / 0.2e1 + (-m(4) * t30 + m(5) * t27 - t67 + t70) * t130 + t51 * t321 / 0.2e1 + t52 * t322 / 0.2e1 + t112 * t211 + t246 * t399 + t147 * (Ifges(6,5) * t244 - Ifges(6,6) * t243) / 0.2e1 + (m(4) * t29 - m(5) * t23 + t68 - t69) * t131 + m(6) * (t1 * t34 + t105 * t13 + t18 * t4 + t19 * t3 + t2 * t33 + t38 * t66) + (-t408 / 0.2e1 - t396 / 0.2e1 - t397 / 0.2e1 - t93 * mrSges(5,1) - Ifges(4,1) * t358 + Ifges(5,2) * t359 + Ifges(5,6) * t360 - Ifges(4,4) * t361 - Ifges(6,5) * t363 + t100 / 0.2e1 + t410 * t354 - t400) * t121 - t149 * t302 / 0.2e1 - t86 * Ifges(4,1) * t165 + t128 * (Ifges(6,4) * t244 - Ifges(6,2) * t243) / 0.2e1 + (-mrSges(3,1) * t349 - mrSges(3,2) * t294 + Ifges(3,5) * t222 + Ifges(3,6) * t353) * qJDD(2) + Ifges(2,3) * qJDD(1) + (mrSges(5,1) * t106 - mrSges(4,3) * t116 - Ifges(4,4) * t358 + Ifges(5,6) * t359 + Ifges(5,3) * t360 - Ifges(4,2) * t361 + t99 / 0.2e1 - t101 / 0.2e1 + t409 * t354 - t401) * t122 + (t143 * mrSges(4,1) + t23 * mrSges(5,1) - t17 * mrSges(5,2) - t29 * mrSges(4,3) - t13 * t258 + t254 * t368 + t255 * t369 + t256 * t370 + t51 * t276 + (Ifges(4,2) + Ifges(5,3)) * t87 + (Ifges(4,4) + Ifges(5,6)) * t86 + t409 * t380) * t164 + (m(4) * t115 + t376) * t58 + (-t142 * mrSges(6,1) + t141 * mrSges(6,2) + ((m(4) + m(5)) * t226 + t374) * t225 + (m(4) * t207 - m(6) * t272 - t216 * t268 - m(5) * (t272 - t204) - t377) * t223) * g(1) + t33 * t14 + t34 * t15 - t381 * t298 + (t245 * qJD(2) / 0.2e1 - t385) * qJD(2) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t386) + (-m(5) * (t271 + t391) - m(4) * t271 - m(6) * (pkin(8) * t313 + t184 + t391) - t140 * mrSges(6,1) - t139 * mrSges(6,2) - mrSges(6,3) * t313 + t377 * t225 + t374 * t223) * g(2) + (-m(4) * t116 + m(5) * t106 - t393) * t57 + t318 * t371 + t38 * t65 + (qJD(1) * (-Ifges(3,2) * t222 + t288) + t150) * t279 / 0.2e1 + (qJD(5) * t52 + t7) * t317 / 0.2e1 + m(4) * (-t143 * t207 - t178 * t211) + t3 * t90 + t4 * t91 + t105 * t10 + t46 * t113 + t114 * (-mrSges(5,2) * t87 + mrSges(5,3) * t86) + (Ifges(3,1) * t170 + Ifges(3,4) * t399) * t222 - t207 * (mrSges(4,1) * t87 - mrSges(4,2) * t86) + (Ifges(6,1) * t244 - Ifges(6,4) * t243) * t363 + t170 * t288 / 0.2e1; Ifges(3,6) * t407 + t228 - t112 * t210 + (t13 * t202 - t18 * t21 - t19 * t22 + (-t89 + t193) * t66) * m(6) + t132 * t269 + t67 * t293 + (m(4) * t350 + t247 + t259) * t383 + (-t202 * t23 + t206 * t27 - g(1) * (t225 * t262 + t183) - g(2) * (t223 * t262 + t181) - t118 * t93 - t92 * t94 + (t119 - t193) * t106) * m(5) + (m(6) * (t18 * t224 + t19 * t220) + t224 * t91 + t220 * t90 + t376) * pkin(2) * t301 + t149 * t303 / 0.2e1 - t245 * t298 / 0.2e1 + t373 * (-pkin(8) + t206) + (-m(5) * t250 - m(4) * t295 + t177 - m(6) * (t203 + t250) + t379) * g(3) + (-t223 * t231 + t387) * g(2) + (-t225 * t231 + t388) * g(1) + t392 * t118 - t393 * t119 + (qJD(1) * t381 + t385) * qJD(1) + t68 * t351 - (-Ifges(3,2) * t303 + t150 + t209) * t280 / 0.2e1 + t394 * t193 + (-t115 * t118 - t116 * t119 + t178 * t210 + (t352 * t30 + t221 * t29 + (t115 * t221 + t116 * t352) * qJD(3)) * pkin(2)) * m(4) - t89 * t65 - t22 * t90 - t21 * t91 - t94 * t113 + t398 * t202 - t162 * mrSges(3,2) - t163 * mrSges(3,1) + Ifges(3,5) * t170 + t206 * t70 + Ifges(3,3) * qJDD(2); -pkin(3) * t70 - t110 * t113 - t24 * t91 - t25 * t90 + t248 * t65 + t228 + t383 * t259 + t392 * t116 + t393 * t115 + t394 * qJD(4) + t398 * qJ(4) + (-t223 * t249 + t387) * g(2) + (-t225 * t249 + t388) * g(1) + (t13 * qJ(4) - t18 * t24 - t19 * t25 + t402 * t66) * m(6) + (-pkin(3) * t27 - qJ(4) * t23 - t110 * t92 - t116 * t93 - g(1) * (-pkin(3) * t315 + t183) - g(2) * (-pkin(3) * t316 + t181) + t382 * t106) * m(5) + (-m(6) * (t203 + t304) - m(5) * t304 + t379) * g(3) - t373 * t367; -t394 * t218 + t252 * qJD(5) + (t113 + t252) * t152 + t70 + (-t152 * t253 - t218 * t66 + t230 + t378) * m(6) + (t106 * t218 + t152 * t92 + t27 + t378) * m(5) + t384; -t66 * (mrSges(6,1) * t129 + mrSges(6,2) * t128) + (Ifges(6,1) * t128 - t334) * t364 + t51 * t363 + (Ifges(6,5) * t128 - Ifges(6,6) * t129) * t362 - t18 * t90 + t19 * t91 - g(1) * (mrSges(6,1) * t139 - mrSges(6,2) * t140) - g(2) * (mrSges(6,1) * t141 + mrSges(6,2) * t142) + t258 * t346 + (t128 * t18 + t129 * t19) * mrSges(6,3) + t297 + (-Ifges(6,2) * t129 + t126 + t52) * t365 + t375;];
tau = t5;
