% Calculate vector of inverse dynamics joint torques for
% S5RRRRP7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRP7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP7_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP7_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP7_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP7_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP7_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:56:08
% EndTime: 2019-12-31 21:56:38
% DurationCPUTime: 15.61s
% Computational Cost: add. (7048->597), mult. (15588->773), div. (0->0), fcn. (10533->10), ass. (0->275)
t485 = -mrSges(5,3) - mrSges(6,2);
t482 = Ifges(5,1) + Ifges(6,1);
t487 = -Ifges(5,4) + Ifges(6,5);
t457 = Ifges(6,4) + Ifges(5,5);
t456 = Ifges(6,2) + Ifges(5,3);
t481 = Ifges(5,6) - Ifges(6,6);
t250 = sin(qJ(4));
t338 = qJD(4) * t250;
t251 = sin(qJ(3));
t252 = sin(qJ(2));
t255 = cos(qJ(3));
t256 = cos(qJ(2));
t197 = t251 * t252 - t255 * t256;
t188 = t197 * qJD(1);
t367 = t188 * t250;
t486 = t338 + t367;
t254 = cos(qJ(4));
t337 = qJD(4) * t254;
t366 = t188 * t254;
t469 = t366 + t337;
t336 = qJD(1) * qJD(2);
t204 = qJDD(1) * t256 - t252 * t336;
t205 = qJDD(1) * t252 + t256 * t336;
t267 = t197 * qJD(3);
t116 = -qJD(1) * t267 + t204 * t251 + t205 * t255;
t248 = qJDD(2) + qJDD(3);
t198 = t251 * t256 + t252 * t255;
t189 = t198 * qJD(1);
t334 = qJD(2) + qJD(3);
t166 = t250 * t189 - t254 * t334;
t339 = qJD(4) * t166;
t64 = t254 * t116 + t250 * t248 - t339;
t423 = t64 / 0.2e1;
t167 = t254 * t189 + t250 * t334;
t65 = qJD(4) * t167 + t250 * t116 - t254 * t248;
t421 = t65 / 0.2e1;
t268 = t198 * qJD(3);
t117 = -qJD(1) * t268 + t204 * t255 - t205 * t251;
t115 = qJDD(4) - t117;
t419 = t115 / 0.2e1;
t460 = mrSges(5,1) + mrSges(6,1);
t459 = mrSges(5,2) - mrSges(6,3);
t184 = qJD(4) + t188;
t386 = mrSges(6,2) * t166;
t119 = mrSges(6,3) * t184 - t386;
t384 = mrSges(5,3) * t166;
t120 = -mrSges(5,2) * t184 - t384;
t484 = -t120 - t119;
t383 = mrSges(5,3) * t167;
t121 = mrSges(5,1) * t184 - t383;
t385 = mrSges(6,2) * t167;
t122 = -mrSges(6,1) * t184 + t385;
t483 = t122 - t121;
t480 = t457 * t115 + t482 * t64 + t487 * t65;
t479 = -t481 * t166 + t457 * t167 + t456 * t184;
t164 = Ifges(5,4) * t166;
t377 = Ifges(6,5) * t166;
t454 = t482 * t167 + t457 * t184 - t164 + t377;
t478 = t334 * Ifges(4,5);
t477 = t334 * Ifges(4,6);
t372 = t189 * mrSges(4,3);
t476 = -mrSges(4,1) * t334 + mrSges(5,1) * t166 + mrSges(5,2) * t167 + t372;
t475 = t486 * pkin(4) - t469 * qJ(5) - qJD(5) * t250;
t249 = qJ(2) + qJ(3);
t246 = sin(t249);
t474 = t485 * t246;
t281 = pkin(4) * t254 + qJ(5) * t250;
t208 = -pkin(3) - t281;
t291 = -t254 * mrSges(6,1) - t250 * mrSges(6,3);
t387 = mrSges(5,1) * t254;
t442 = (-m(6) * t208 - t291 + t387) * t246;
t473 = -t481 * t250 + t457 * t254;
t376 = Ifges(6,5) * t250;
t379 = Ifges(5,4) * t250;
t472 = -t482 * t254 - t376 + t379;
t258 = -pkin(7) - pkin(6);
t214 = t258 * t256;
t201 = qJD(1) * t214;
t191 = t255 * t201;
t213 = t258 * t252;
t200 = qJD(1) * t213;
t153 = t200 * t251 - t191;
t341 = qJD(3) * t251;
t471 = pkin(2) * t341 - t153;
t402 = pkin(2) * t256;
t241 = pkin(1) + t402;
t212 = t241 * qJD(1);
t124 = pkin(3) * t188 - pkin(8) * t189 - t212;
t193 = qJD(2) * pkin(2) + t200;
t151 = t251 * t193 - t191;
t136 = pkin(8) * t334 + t151;
t370 = qJDD(1) * pkin(1);
t181 = -pkin(2) * t204 - t370;
t30 = -pkin(3) * t117 - pkin(8) * t116 + t181;
t196 = t205 * pkin(6);
t161 = qJDD(2) * pkin(2) - pkin(7) * t205 - t196;
t195 = t204 * pkin(6);
t165 = pkin(7) * t204 + t195;
t340 = qJD(3) * t255;
t40 = t251 * t161 + t255 * t165 + t193 * t340 + t201 * t341;
t37 = pkin(8) * t248 + t40;
t7 = t124 * t337 - t136 * t338 + t250 * t30 + t254 * t37;
t58 = t124 * t250 + t136 * t254;
t8 = -qJD(4) * t58 - t250 * t37 + t254 * t30;
t468 = -t250 * t8 + t254 * t7;
t2 = qJ(5) * t115 + qJD(5) * t184 + t7;
t4 = -pkin(4) * t115 + qJDD(5) - t8;
t467 = t2 * t254 + t250 * t4;
t190 = t251 * t201;
t150 = t255 * t193 + t190;
t135 = -pkin(3) * t334 - t150;
t290 = t250 * mrSges(6,1) - t254 * mrSges(6,3);
t292 = mrSges(5,1) * t250 + mrSges(5,2) * t254;
t52 = t166 * pkin(4) - t167 * qJ(5) + t135;
t466 = t135 * t292 + t290 * t52;
t465 = Ifges(6,5) * t423 + Ifges(6,6) * t419 - t64 * Ifges(5,4) / 0.2e1 - t115 * Ifges(5,6) / 0.2e1 + (Ifges(6,3) + Ifges(5,2)) * t421;
t25 = -mrSges(6,2) * t65 + mrSges(6,3) * t115;
t26 = mrSges(5,1) * t115 - mrSges(5,3) * t64;
t27 = -t115 * mrSges(6,1) + t64 * mrSges(6,2);
t28 = -mrSges(5,2) * t115 - mrSges(5,3) * t65;
t57 = t124 * t254 - t136 * t250;
t33 = -pkin(4) * t184 + qJD(5) - t57;
t35 = qJ(5) * t184 + t58;
t464 = t250 * (-t26 + t27) + t254 * (t25 + t28) + m(5) * ((-t250 * t58 - t254 * t57) * qJD(4) + t468) + m(6) * ((-t250 * t35 + t254 * t33) * qJD(4) + t467) + t484 * t338 + t483 * t337;
t374 = t167 * Ifges(5,4);
t85 = -t166 * Ifges(5,2) + t184 * Ifges(5,6) + t374;
t463 = -t85 / 0.2e1;
t462 = m(6) + m(5);
t461 = t204 / 0.2e1;
t406 = t256 / 0.2e1;
t453 = t256 * Ifges(3,2);
t451 = -t151 + t475;
t450 = t471 + t475;
t149 = pkin(3) * t197 - pkin(8) * t198 - t241;
t169 = t213 * t251 - t214 * t255;
t447 = t250 * t149 + t254 * t169;
t446 = t255 * t213 + t214 * t251;
t247 = cos(t249);
t354 = t247 * t254;
t356 = t247 * t250;
t445 = pkin(4) * t354 + qJ(5) * t356;
t444 = t247 * mrSges(4,1) - t246 * mrSges(4,2);
t345 = t247 * pkin(3) + t246 * pkin(8);
t401 = pkin(3) * t246;
t404 = pkin(2) * t252;
t443 = m(6) * t404 - m(5) * (-t401 - t404) + t442;
t257 = cos(qJ(1));
t353 = t247 * t257;
t227 = pkin(8) * t353;
t358 = t246 * t250;
t329 = mrSges(5,2) * t358;
t441 = -m(6) * t227 - t257 * t329 + t353 * t485;
t253 = sin(qJ(1));
t355 = t247 * t253;
t224 = pkin(8) * t355;
t440 = -m(6) * t224 - t253 * t329 + t355 * t485;
t157 = -qJD(2) * t197 - t267;
t271 = t157 * t250 + t198 * t337;
t438 = t115 * t456 + t457 * t64 - t481 * t65;
t343 = qJD(1) * t256;
t344 = qJD(1) * t252;
t398 = pkin(6) * t256;
t399 = pkin(6) * t252;
t437 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t344) * t398 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t343) * t399;
t436 = t195 * t256 + t196 * t252;
t435 = g(1) * t257 + g(2) * t253;
t433 = -m(3) * pkin(6) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t211 = -mrSges(3,1) * t256 + mrSges(3,2) * t252;
t432 = -m(3) * pkin(1) - mrSges(2,1) + t211 - t444;
t431 = -t354 * t460 + t356 * t459 - t444 + t474;
t158 = qJD(2) * t198 + t268;
t342 = qJD(2) * t252;
t327 = pkin(2) * t342;
t91 = pkin(3) * t158 - pkin(8) * t157 + t327;
t320 = qJD(2) * t258;
t202 = t252 * t320;
t203 = t256 * t320;
t98 = qJD(3) * t446 + t202 * t255 + t203 * t251;
t21 = -qJD(4) * t447 - t250 * t98 + t254 * t91;
t428 = m(6) * pkin(4) + t460;
t427 = m(6) * qJ(5) - t459;
t426 = t8 * mrSges(5,1) - t4 * mrSges(6,1) - t7 * mrSges(5,2) + t2 * mrSges(6,3);
t422 = -t65 / 0.2e1;
t417 = -t166 / 0.2e1;
t416 = t166 / 0.2e1;
t415 = -t167 / 0.2e1;
t414 = t167 / 0.2e1;
t413 = -t184 / 0.2e1;
t411 = t188 / 0.2e1;
t409 = t189 / 0.2e1;
t407 = t250 / 0.2e1;
t405 = pkin(2) * t251;
t403 = pkin(2) * t255;
t400 = pkin(4) * t189;
t393 = g(3) * t246;
t382 = Ifges(3,4) * t252;
t381 = Ifges(3,4) * t256;
t380 = Ifges(4,4) * t189;
t378 = Ifges(5,4) * t254;
t375 = Ifges(6,5) * t254;
t373 = t188 * mrSges(4,3);
t368 = t157 * t254;
t357 = t246 * t257;
t352 = t250 * t253;
t351 = t253 * t254;
t350 = t254 * t257;
t349 = t257 * t250;
t147 = pkin(3) * t189 + pkin(8) * t188;
t81 = t250 * t147 + t254 * t150;
t328 = pkin(2) * t344;
t130 = t147 + t328;
t154 = t200 * t255 + t190;
t77 = t250 * t130 + t254 * t154;
t325 = pkin(2) * t340;
t163 = Ifges(6,5) * t167;
t82 = t184 * Ifges(6,6) + t166 * Ifges(6,3) + t163;
t322 = t82 * t407;
t307 = -t338 / 0.2e1;
t306 = t337 / 0.2e1;
t303 = t336 / 0.2e1;
t301 = t257 * t241 - t253 * t258;
t300 = t250 * t325;
t299 = t254 * t325;
t297 = t345 + t402;
t295 = mrSges(3,1) * t252 + mrSges(3,2) * t256;
t293 = mrSges(4,1) * t246 + mrSges(4,2) * t247;
t287 = t382 + t453;
t286 = -Ifges(5,2) * t250 + t378;
t284 = Ifges(3,5) * t256 - Ifges(3,6) * t252;
t282 = Ifges(6,3) * t250 + t375;
t280 = pkin(4) * t250 - qJ(5) * t254;
t76 = t130 * t254 - t154 * t250;
t80 = t147 * t254 - t150 * t250;
t92 = t149 * t254 - t169 * t250;
t41 = t161 * t255 - t251 * t165 - t193 * t341 + t201 * t340;
t272 = pkin(1) * t295;
t270 = t198 * t338 - t368;
t269 = t252 * (Ifges(3,1) * t256 - t382);
t20 = t149 * t337 - t169 * t338 + t250 * t91 + t254 * t98;
t38 = -pkin(3) * t248 - t41;
t99 = qJD(3) * t169 + t202 * t251 - t255 * t203;
t131 = -Ifges(4,2) * t188 + t380 + t477;
t182 = Ifges(4,4) * t188;
t132 = Ifges(4,1) * t189 - t182 + t478;
t9 = pkin(4) * t65 - qJ(5) * t64 - qJD(5) * t167 + t38;
t259 = (t463 + t82 / 0.2e1) * t367 + (-t212 * mrSges(4,2) + t478 / 0.2e1 - t286 * t416 - t282 * t417 + t472 * t415 - t473 * t413 + t466) * t188 + (t322 + t466) * qJD(4) + (-t286 / 0.2e1 + t282 / 0.2e1) * t339 + (-t167 * t472 + t184 * t473) * qJD(4) / 0.2e1 + (t33 * t469 - t35 * t486 + t467) * mrSges(6,2) + (-t469 * t57 - t486 * t58 + t468) * mrSges(5,3) - t38 * t387 - (-Ifges(4,1) * t188 - t380 + t479) * t189 / 0.2e1 + t480 * t407 + (Ifges(5,2) * t422 - Ifges(6,3) * t421 + t419 * t481 - t465) * t254 + (t212 * mrSges(4,1) + t477 / 0.2e1 + Ifges(5,6) * t416 + Ifges(6,6) * t417 - t57 * mrSges(5,1) + t33 * mrSges(6,1) + t58 * mrSges(5,2) - t35 * mrSges(6,3) - Ifges(4,2) * t411 + t457 * t415 + t456 * t413) * t189 + (t38 * mrSges(5,2) + t457 * t419 + t423 * t482) * t250 + t376 * t421 + t379 * t422 + (-t375 + t378) * t423 + (t306 + t366 / 0.2e1) * t454 + t9 * t291 + (-t182 + t132) * t411 + t131 * t409 + t85 * t307 + t151 * t372 - t150 * t373 - t40 * mrSges(4,2) + t41 * mrSges(4,1) + Ifges(4,5) * t116 + Ifges(4,6) * t117 + Ifges(4,3) * t248;
t243 = Ifges(3,4) * t343;
t240 = -pkin(3) - t403;
t194 = t208 - t403;
t187 = Ifges(3,1) * t344 + Ifges(3,5) * qJD(2) + t243;
t186 = Ifges(3,6) * qJD(2) + qJD(1) * t287;
t180 = t247 * t350 + t352;
t179 = t247 * t349 - t351;
t178 = t247 * t351 - t349;
t177 = t247 * t352 + t350;
t176 = t189 * qJ(5);
t170 = -mrSges(4,2) * t334 - t373;
t146 = mrSges(4,1) * t188 + mrSges(4,2) * t189;
t106 = -mrSges(4,2) * t248 + mrSges(4,3) * t117;
t105 = mrSges(4,1) * t248 - mrSges(4,3) * t116;
t103 = mrSges(6,1) * t166 - mrSges(6,3) * t167;
t102 = pkin(4) * t167 + qJ(5) * t166;
t97 = t198 * t280 - t446;
t75 = -pkin(4) * t197 - t92;
t70 = qJ(5) * t197 + t447;
t54 = -t80 - t400;
t53 = t176 + t81;
t47 = -t76 - t400;
t46 = t176 + t77;
t24 = mrSges(5,1) * t65 + mrSges(5,2) * t64;
t23 = mrSges(6,1) * t65 - mrSges(6,3) * t64;
t22 = t280 * t157 + (qJD(4) * t281 - qJD(5) * t254) * t198 + t99;
t11 = -pkin(4) * t158 - t21;
t10 = qJ(5) * t158 + qJD(5) * t197 + t20;
t1 = [(t428 * t178 + t427 * t177 + (m(4) * t241 - t462 * (-t241 - t345) - t432 - t474) * t253 + (t433 + (m(4) + t462) * t258) * t257) * g(1) + (-m(4) * t150 + m(5) * t135 + t476) * t99 + (t158 * t456 - t270 * t457 - t271 * t481) * t184 / 0.2e1 + (-t150 * t157 - t151 * t158) * mrSges(4,3) + (Ifges(3,4) * t205 + Ifges(3,2) * t204) * t406 + (-m(4) * t301 + t485 * t357 - t462 * (pkin(3) * t353 + pkin(8) * t357 + t301) - t428 * t180 - t427 * t179 + t432 * t257 + t433 * t253) * g(2) - (-m(4) * t41 + m(5) * t38 - t105 + t24) * t446 + t287 * t461 + t271 * t463 + (Ifges(3,1) * t205 + Ifges(3,4) * t461 + Ifges(3,5) * qJDD(2) - t303 * t453) * t252 + (-mrSges(3,1) * t399 - mrSges(3,2) * t398 + 0.2e1 * Ifges(3,6) * t406) * qJDD(2) + t479 * t158 / 0.2e1 + t454 * t368 / 0.2e1 + (t204 * t398 + t205 * t399 + t436) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t436) + (t187 * t406 + t284 * qJD(2) / 0.2e1 - t437) * qJD(2) + (t181 * mrSges(4,2) - t41 * mrSges(4,3) + Ifges(4,1) * t116 + Ifges(4,4) * t117 + Ifges(4,5) * t248 + t282 * t421 + t286 * t422 + t290 * t9 + t292 * t38 + t306 * t82 + t419 * t473 - t423 * t472 + (-t2 * mrSges(6,2) - t7 * mrSges(5,3) + t465) * t250 + t454 * t307 + (mrSges(6,2) * t4 - mrSges(5,3) * t8 + t480 / 0.2e1) * t254) * t198 + m(4) * (t151 * t98 + t169 * t40 - t181 * t241 - t212 * t327) + (t181 * mrSges(4,1) - t40 * mrSges(4,3) - Ifges(4,4) * t116 - Ifges(4,2) * t117 - Ifges(4,6) * t248 + Ifges(5,6) * t422 + Ifges(6,6) * t421 + t419 * t456 + t423 * t457 + t426 + t438 / 0.2e1) * t197 + t334 * (Ifges(4,5) * t157 - Ifges(4,6) * t158) / 0.2e1 - t272 * t336 + t57 * (mrSges(5,1) * t158 + mrSges(5,3) * t270) + (t256 * t381 + t269) * t303 + t146 * t327 + m(5) * (t20 * t58 + t21 * t57 + t447 * t7 + t8 * t92) + t447 * t28 + t157 * t322 - t211 * t370 - t186 * t342 / 0.2e1 + (Ifges(4,1) * t157 - Ifges(4,4) * t158) * t409 + (-Ifges(6,5) * t270 + Ifges(6,6) * t158 + Ifges(6,3) * t271) * t416 + (-Ifges(5,4) * t270 - Ifges(5,2) * t271 + Ifges(5,6) * t158) * t417 + t205 * t381 / 0.2e1 + (t457 * t158 - t482 * t270 + t487 * t271) * t414 + t33 * (-mrSges(6,1) * t158 - mrSges(6,2) * t270) + t52 * (mrSges(6,1) * t271 + mrSges(6,3) * t270) + t135 * (mrSges(5,1) * t271 - mrSges(5,2) * t270) + t58 * (-mrSges(5,2) * t158 - mrSges(5,3) * t271) + t35 * (-mrSges(6,2) * t271 + mrSges(6,3) * t158) + t70 * t25 + t75 * t27 + t92 * t26 + t97 * t23 + t22 * t103 + Ifges(2,3) * qJDD(1) + t10 * t119 + t20 * t120 + t21 * t121 + t11 * t122 + t157 * t132 / 0.2e1 - t158 * t131 / 0.2e1 + t169 * t106 + m(6) * (t10 * t35 + t11 * t33 + t2 * t70 + t22 * t52 + t4 * t75 + t9 * t97) + t98 * t170 - t188 * (Ifges(4,4) * t157 - Ifges(4,2) * t158) / 0.2e1 - pkin(1) * (-mrSges(3,1) * t204 + mrSges(3,2) * t205) - t212 * (mrSges(4,1) * t158 + mrSges(4,2) * t157) - t241 * (-mrSges(4,1) * t117 + mrSges(4,2) * t116); (t240 * t38 + (t135 * t251 + (-t250 * t57 + t254 * t58) * t255) * qJD(3) * pkin(2) - t135 * t153 - t57 * t76 - t58 * t77 - g(2) * t224 - g(1) * t227) * m(5) + (t325 - t154) * t170 + (t300 - t47) * t122 + (-t300 - t76) * t121 + (t299 - t77) * t120 + (t299 - t46) * t119 + (m(4) * t404 + t293 + t295) * t435 + t476 * t471 + t464 * (pkin(8) + t405) + t450 * t103 + (t194 * t9 + (t250 * t33 + t254 * t35) * t325 - t33 * t47 - t35 * t46 + t450 * t52) * m(6) + (-m(4) * t402 - m(5) * t297 + t211 - m(6) * (t297 + t445) + t431) * g(3) + (t253 * t443 + t440) * g(2) + (t257 * t443 + t441) * g(1) + (t150 * t153 - t151 * t154 + t212 * t328 + (t251 * t40 + t255 * t41 + (-t150 * t251 + t151 * t255) * qJD(3)) * pkin(2)) * m(4) - (-Ifges(3,2) * t344 + t187 + t243) * t343 / 0.2e1 - t284 * t336 / 0.2e1 - t146 * t328 + t106 * t405 + t105 * t403 + t259 + (t437 + (-t269 / 0.2e1 + t272) * qJD(1)) * qJD(1) + t186 * t344 / 0.2e1 + Ifges(3,3) * qJDD(2) + t194 * t23 - t195 * mrSges(3,2) - t196 * mrSges(3,1) + Ifges(3,6) * t204 + Ifges(3,5) * t205 + t240 * t24; -pkin(3) * t24 - t53 * t119 - t81 * t120 - t80 * t121 - t54 * t122 - t150 * t170 + t208 * t23 + t259 + t435 * t293 - t476 * t151 + t451 * t103 + (t253 * t442 + t440) * g(2) + (t257 * t442 + t441) * g(1) + (t208 * t9 - t33 * t54 - t35 * t53 + t451 * t52) * m(6) + (-g(2) * (-t253 * t401 + t224) - g(1) * (-pkin(3) * t357 + t227) - t135 * t151 - t57 * t80 - t58 * t81 - pkin(3) * t38) * m(5) + (-m(6) * (t345 + t445) - m(5) * t345 + t431) * g(3) + t464 * pkin(8); (-t166 * t482 + t163 - t374 + t82) * t415 + (t179 * t460 + t180 * t459) * g(1) + (t460 * t177 + t459 * t178) * g(2) + (-t166 * t457 - t167 * t481) * t413 + (t292 + t290) * t393 + (t280 * t393 - t102 * t52 - pkin(4) * t4 + qJ(5) * t2 + qJD(5) * t35 - g(2) * (-pkin(4) * t177 + qJ(5) * t178) - g(1) * (-pkin(4) * t179 + qJ(5) * t180)) * m(6) + t426 + t438 + t35 * t385 + t33 * t386 + (-m(6) * t33 + t383 - t483) * t58 + (-m(6) * t35 - t384 + t484) * t57 + (-Ifges(5,2) * t167 - t164 + t454) * t416 + t85 * t414 + (Ifges(6,3) * t167 - t377) * t417 + qJ(5) * t25 - pkin(4) * t27 - t102 * t103 + qJD(5) * t119 - t52 * (mrSges(6,1) * t167 + mrSges(6,3) * t166) - t135 * (mrSges(5,1) * t167 - mrSges(5,2) * t166); t167 * t103 - t184 * t119 + (-g(1) * t179 - g(2) * t177 - g(3) * t358 + t52 * t167 - t35 * t184 + t4) * m(6) + t27;];
tau = t1;
