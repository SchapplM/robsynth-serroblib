% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRPP2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:31:33
% EndTime: 2019-03-09 04:31:47
% DurationCPUTime: 7.13s
% Computational Cost: add. (6640->541), mult. (13899->693), div. (0->0), fcn. (10881->6), ass. (0->269)
t280 = sin(qJ(4));
t277 = t280 ^ 2;
t282 = cos(qJ(4));
t278 = t282 ^ 2;
t395 = t277 + t278;
t272 = Ifges(6,5) * t280;
t327 = t282 * Ifges(6,3) - t272;
t273 = Ifges(7,4) * t280;
t329 = t282 * Ifges(7,2) - t273;
t493 = t329 + t327;
t400 = t282 * qJ(5);
t463 = pkin(4) + pkin(5);
t511 = t463 * t280;
t512 = t511 - t400;
t504 = Ifges(6,4) + Ifges(5,5);
t503 = Ifges(5,6) + Ifges(7,6);
t424 = Ifges(6,5) * t282;
t425 = Ifges(7,4) * t282;
t494 = t424 + t425 + (Ifges(7,2) + Ifges(6,3)) * t280;
t496 = Ifges(7,1) * t282 + t273;
t497 = Ifges(6,1) * t282 + t272;
t510 = t496 + t497;
t405 = t280 * qJ(5);
t509 = t282 * t463 + t405;
t283 = cos(qJ(3));
t399 = t282 * t283;
t381 = mrSges(7,3) * t399;
t281 = sin(qJ(3));
t417 = t281 * mrSges(7,1);
t211 = -t381 - t417;
t508 = -t211 / 0.2e1;
t507 = -t281 / 0.2e1;
t454 = t281 / 0.2e1;
t450 = -t283 / 0.2e1;
t448 = t283 / 0.2e1;
t442 = m(7) * t281;
t378 = -t442 / 0.2e1;
t502 = pkin(8) - qJ(6);
t217 = t502 * t280;
t222 = t502 * t282;
t506 = m(7) * (t217 * t280 + t222 * t282);
t505 = -mrSges(5,1) - mrSges(6,1);
t276 = t283 * pkin(4);
t370 = -cos(pkin(9)) * pkin(1) - pkin(2);
t439 = pkin(8) * t281;
t195 = -pkin(3) * t283 + t370 - t439;
t262 = sin(pkin(9)) * pkin(1) + pkin(7);
t403 = t280 * t283;
t397 = -t282 * t195 + t262 * t403;
t402 = t281 * t282;
t58 = qJ(6) * t402 - t397;
t43 = pkin(5) * t283 + t276 - t58;
t429 = t58 + t43;
t69 = t276 + t397;
t501 = t69 - t397;
t500 = mrSges(7,3) * t395;
t219 = t282 * mrSges(7,1) + t280 * mrSges(7,2);
t337 = t282 * mrSges(6,1) + t280 * mrSges(6,3);
t499 = -t337 - t219;
t271 = t282 * mrSges(7,2);
t422 = t280 * mrSges(7,1);
t498 = t271 - t422;
t274 = Ifges(5,4) * t282;
t495 = -Ifges(5,2) * t280 + t274;
t336 = mrSges(6,1) * t280 - mrSges(6,3) * t282;
t440 = pkin(4) * t280;
t221 = -t400 + t440;
t446 = m(6) * t221;
t492 = t336 + t446;
t404 = t280 * t281;
t255 = mrSges(7,3) * t404;
t414 = t283 * mrSges(7,2);
t209 = t255 - t414;
t365 = t404 / 0.2e1;
t491 = mrSges(7,3) * t365 - t209 / 0.2e1;
t490 = -t395 * t439 / 0.2e1;
t334 = -t280 * Ifges(5,1) - t274;
t488 = Ifges(6,2) + Ifges(5,3) + Ifges(7,3);
t151 = t283 * Ifges(7,5) + t281 * t496;
t153 = -t283 * Ifges(6,4) + t281 * t497;
t426 = Ifges(5,4) * t280;
t335 = Ifges(5,1) * t282 - t426;
t303 = t335 * t281;
t155 = -t283 * Ifges(5,5) + t303;
t487 = t155 + t153 + t151;
t384 = mrSges(6,2) * t404;
t309 = -t283 * mrSges(6,3) - t384;
t386 = mrSges(5,3) * t404;
t310 = t283 * mrSges(5,2) - t386;
t385 = mrSges(5,3) * t402;
t313 = -t283 * mrSges(5,1) - t385;
t486 = t282 * t313 + (t309 + t310) * t280;
t193 = t336 * t281;
t92 = (t221 + t262) * t281;
t485 = m(6) * t92 + t193;
t411 = qJ(5) * t283;
t79 = t280 * t195 + t262 * t399;
t68 = t79 - t411;
t484 = m(6) * t68 + t309;
t466 = Ifges(6,6) / 0.2e1;
t344 = -Ifges(7,6) / 0.2e1 + t466 - Ifges(5,6) / 0.2e1;
t438 = pkin(8) * t283;
t441 = pkin(3) * t281;
t233 = -t438 + t441;
t91 = t280 * t233 - t262 * t402;
t81 = t281 * qJ(5) + t91;
t483 = -t81 * mrSges(6,2) - t91 * mrSges(5,3) + Ifges(6,6) * t454 + t344 * t281 + t448 * t494 + t495 * t450 + t503 * t507;
t253 = qJ(6) * t404;
t341 = t79 - 0.2e1 * t411;
t431 = -mrSges(6,3) - mrSges(7,2);
t470 = -m(7) / 0.2e1;
t472 = -m(6) / 0.2e1;
t481 = -t283 * t431 + t341 * t472 + (t253 + t341) * t470;
t226 = t280 * Ifges(7,1) - t425;
t227 = t280 * Ifges(6,1) - t424;
t480 = -t334 + t495 + t227 + t226;
t258 = Ifges(6,5) * t402;
t145 = -t283 * Ifges(6,6) + Ifges(6,3) * t404 + t258;
t259 = Ifges(7,4) * t402;
t147 = Ifges(7,2) * t404 + t283 * Ifges(7,6) + t259;
t479 = t145 + t147 + t258 + t259 + (-Ifges(6,1) - Ifges(7,1)) * t404;
t478 = t510 - t493;
t477 = 0.2e1 * m(7);
t476 = t281 ^ 2;
t475 = 2 * qJD(3);
t474 = 2 * qJD(4);
t473 = -m(5) / 0.2e1;
t471 = m(6) / 0.2e1;
t469 = m(7) / 0.2e1;
t467 = mrSges(7,2) / 0.2e1;
t465 = t79 / 0.2e1;
t464 = m(6) + m(7);
t194 = t498 * t281;
t461 = t194 / 0.2e1;
t460 = t512 / 0.4e1;
t256 = mrSges(6,2) * t399;
t418 = t281 * mrSges(6,1);
t212 = t256 - t418;
t458 = -t212 / 0.2e1;
t457 = -t280 / 0.2e1;
t453 = -t282 / 0.2e1;
t452 = t282 / 0.2e1;
t326 = pkin(4) * t282 + t405;
t216 = -pkin(3) - t326;
t447 = m(6) * t216;
t445 = t509 * t442;
t444 = m(7) * t512;
t443 = m(7) * t280;
t354 = -t262 * t280 - pkin(4);
t53 = (-qJ(6) * t283 - t233) * t282 + (-pkin(5) + t354) * t281;
t437 = t53 * mrSges(7,1);
t436 = t81 * mrSges(6,3);
t407 = t233 * t282;
t82 = t281 * t354 - t407;
t435 = t82 * mrSges(6,1);
t90 = t262 * t404 + t407;
t434 = t90 * mrSges(5,1);
t433 = t91 * mrSges(5,2);
t432 = mrSges(6,2) - mrSges(7,3);
t51 = t253 + t68;
t59 = t253 + t79;
t430 = -t51 + t59;
t419 = t280 * t51;
t192 = t326 * t281;
t294 = t501 * t282 + (-t68 + t79) * t280;
t305 = t337 * t281;
t296 = t305 / 0.2e1;
t383 = mrSges(6,2) * t402;
t312 = t283 * mrSges(6,1) + t383;
t298 = t282 * t312;
t304 = t281 * t219;
t339 = t282 * mrSges(5,1) - t280 * mrSges(5,2);
t306 = t339 * t281;
t382 = mrSges(7,3) * t402;
t311 = t283 * mrSges(7,1) - t382;
t362 = -t402 / 0.2e1;
t391 = t476 / 0.2e1;
t7 = t298 * t507 + t311 * t362 + t209 * t365 + t306 * t448 + t283 * t296 - t304 * t450 + (-t192 * t283 + t281 * t294) * t472 + ((-t411 + t430) * t280 + (-t283 * t463 + t429) * t282) * t378 - t476 * t500 / 0.2e1 + t486 * t454 + (mrSges(6,2) + mrSges(5,3)) * t395 * t391;
t413 = t7 * qJD(1);
t412 = -t339 - mrSges(4,1);
t76 = (-t262 - t512) * t281;
t13 = (m(7) * t51 + t209 + t484) * t283 + (-m(7) * t76 - t194 + t485) * t402;
t410 = qJD(1) * t13;
t25 = (t280 * t209 + m(7) * (-t282 * t43 + t419) - mrSges(7,1) * t399 + t278 * mrSges(7,3) * t281) * t281;
t409 = qJD(1) * t25;
t408 = t216 * t281;
t406 = t262 * t281;
t401 = t281 * t283;
t396 = t395 * t438;
t394 = qJD(3) * t280;
t393 = qJD(3) * t283;
t392 = qJD(4) * t281;
t390 = m(6) / 0.4e1 + m(7) / 0.4e1;
t388 = m(7) * t402;
t254 = qJ(5) * t399;
t77 = t254 + (-t262 - t511) * t283;
t387 = t77 * t469;
t380 = t445 / 0.2e1;
t379 = -t444 / 0.2e1;
t376 = -mrSges(6,1) / 0.2e1 - mrSges(7,1) / 0.2e1;
t372 = t271 / 0.2e1;
t360 = t402 / 0.2e1;
t356 = -t193 / 0.2e1 + t461;
t210 = mrSges(7,2) * t281 + mrSges(7,3) * t403;
t213 = -mrSges(6,2) * t403 + mrSges(6,3) * t281;
t355 = -t213 / 0.2e1 - t210 / 0.2e1;
t349 = -t217 * t281 - t51;
t345 = Ifges(7,5) / 0.2e1 - Ifges(6,4) / 0.2e1 - Ifges(5,5) / 0.2e1;
t340 = (t337 / 0.2e1 + t219 / 0.2e1) * t282;
t338 = mrSges(5,1) * t280 + mrSges(5,2) * t282;
t225 = t282 * Ifges(5,2) + t426;
t149 = -t283 * Ifges(5,6) + t281 * t495;
t288 = t82 * mrSges(6,2) - t90 * mrSges(5,3) - t53 * mrSges(7,3) + Ifges(7,5) * t507 - t345 * t281 + t504 * t454 + (t335 + t510) * t448;
t307 = t262 * t338;
t60 = qJ(6) * t403 + t81;
t93 = -t254 + (t262 + t440) * t283;
t3 = t93 * t193 + t77 * t194 + t60 * t209 + t51 * t210 + t43 * t211 + t69 * t212 + t68 * t213 + m(5) * (-t397 * t90 + t79 * t91) + m(6) * (t68 * t81 + t69 * t82 + t92 * t93) + m(7) * (t43 * t53 + t51 * t60 + t76 * t77) + (t370 * mrSges(4,1) - mrSges(5,1) * t397 - t79 * mrSges(5,2) - Ifges(4,4) * t281 + t483 * t280 + t288 * t282) * t281 + (Ifges(4,4) * t283 + t370 * mrSges(4,2) + t437 - t436 + t435 - t434 + t433 + (t76 * mrSges(7,2) + t345 * t283 + t397 * mrSges(5,3) + t155 / 0.2e1 + t153 / 0.2e1 + t151 / 0.2e1 - t92 * mrSges(6,3)) * t282 + (-t76 * mrSges(7,1) - t344 * t283 - t79 * mrSges(5,3) - t149 / 0.2e1 + t147 / 0.2e1 + t145 / 0.2e1 + t92 * mrSges(6,1)) * t280 + (m(5) * t262 ^ 2 + Ifges(4,1) - Ifges(4,2) + 0.2e1 * t307 - t488) * t281) * t283;
t292 = -t414 / 0.2e1 + t491;
t308 = m(5) * (-t280 * t90 + t282 * t91);
t323 = t280 * t82 + t282 * t81;
t6 = (t355 * t282 + (t458 + t508) * t280 - t308 / 0.2e1 + (t323 + t92) * t472 + (t280 * t53 + t282 * t60 - t76) * t470 + t406 * t473 + t356) * t281 + (t292 * t282 + (t280 * t397 + t282 * t79) * t473 + (t280 * t69 + t282 * t68 - t93) * t472 + (t280 * t43 + t282 * t51 + t77) * t470 + m(5) * t262 * t448) * t283;
t325 = t3 * qJD(1) - t6 * qJD(2);
t257 = Ifges(6,6) * t402;
t295 = t281 * t509;
t5 = m(7) * (-t295 * t76 + t43 * t59 + t51 * t58) + t51 * t382 + t43 * t255 + t149 * t362 + (-Ifges(7,5) * t280 + Ifges(7,6) * t282) * t401 / 0.2e1 - (-Ifges(5,5) * t280 - Ifges(5,6) * t282) * t401 / 0.2e1 - t69 * t384 - t68 * t383 - t194 * t295 - t76 * t304 + t92 * t305 + t476 * t262 * t339 + (-Ifges(6,4) * t404 + t257) * t450 + t59 * t311 + t58 * t209 + t485 * t192 + (m(6) * t69 + t312 - t313 - t385) * t79 - (t310 + t386 + t484) * t397 - t487 * t404 / 0.2e1 + t479 * t360 + (t282 * t334 + (t225 + t493) * t280) * t391;
t324 = t5 * qJD(1) - t7 * qJD(2);
t322 = (t222 * t281 - t43) * t280;
t29 = 0.4e1 * (m(5) / 0.4e1 + t390) * (-0.1e1 + t395) * t401;
t321 = -t6 * qJD(1) + t29 * qJD(2);
t196 = pkin(3) + t509;
t36 = -t196 * t443 + (t447 + t499) * t280;
t287 = t356 * t280 + (-t280 * t92 + (-t408 - t438) * t282) * t471 + (t196 * t402 - t222 * t283 + t280 * t76) * t469;
t315 = t469 * t53 + t471 * t82;
t9 = t381 - t256 + (t340 - t376) * t281 + t287 - t315;
t320 = qJD(1) * t9 - qJD(3) * t36;
t35 = t380 + ((-pkin(5) * t282 - t326) * t470 + t219) * t281;
t49 = (-t400 / 0.4e1 + t511 / 0.4e1 + t460) * t477 - t498;
t319 = qJD(1) * t35 + qJD(3) * t49;
t316 = m(6) * t465 + t469 * t59;
t23 = t316 + t481;
t234 = qJ(5) * t464 - t431;
t317 = qJD(1) * t23 - qJD(4) * t234;
t314 = -qJ(5) * t432 - t503;
t302 = pkin(4) * mrSges(6,2) - mrSges(7,3) * t463 + Ifges(7,5) - t504;
t19 = 0.2e1 * t390 * t254 + (t372 - t271 / 0.2e1 + t446 / 0.2e1 + t440 * t472 + (-t511 / 0.4e1 + t460) * t477) * t283;
t291 = -t295 / 0.2e1;
t285 = -t196 * t304 / 0.2e1 + t225 * t362 + t219 * t291 + t216 * t296 - t192 * t337 / 0.2e1 + (-t196 * t295 + t222 * t429 - t512 * t76) * t469 - t512 * t461 + (t192 * t216 + t221 * t92) * t471 + t490 * mrSges(5,3) + (t430 * t469 + t491) * t217 + (t334 + t494) * t404 / 0.4e1 + (-t493 / 0.4e1 + t478 / 0.4e1) * t402 - (t504 * t282 + (-Ifges(5,6) + Ifges(6,6)) * t280) * t283 / 0.4e1 + (t298 / 0.2e1 - t486 / 0.2e1 + t294 * t471) * pkin(8) + (t429 * t453 + t222 * t360 + t419 / 0.2e1 + t59 * t457) * mrSges(7,3) + (t501 * t452 + t490 + (-t68 / 0.2e1 + t465) * t280) * mrSges(6,2) + (t303 + t487) * t282 / 0.4e1 - pkin(3) * t306 / 0.2e1 + t479 * t280 / 0.4e1 - t480 * t404 / 0.4e1 + t307 * t454 + t222 * t311 / 0.2e1 + t76 * t498 / 0.2e1 + t221 * t193 / 0.2e1 + t92 * t336 / 0.2e1 - t280 * t149 / 0.4e1 + t283 * (Ifges(7,5) * t282 + Ifges(7,6) * t280) / 0.4e1;
t286 = (-pkin(4) * t82 + qJ(5) * t81) * t471 + (qJ(5) * t60 - t463 * t53) * t469 + pkin(4) * t458 + t463 * t508 - t437 / 0.2e1 + t60 * t467 + t436 / 0.2e1 - t435 / 0.2e1 + t434 / 0.2e1 - t433 / 0.2e1;
t2 = -t285 + t286 + (t210 + t213) * qJ(5) / 0.2e1 + t488 * t454 + (t466 - t503 / 0.2e1) * t403 + (-Ifges(7,5) / 0.2e1 + t504 / 0.2e1) * t399;
t8 = -t221 * t337 - t512 * t219 + t225 * t457 - pkin(3) * t338 + t492 * t216 + (t498 - t444) * t196 + t494 * t453 + t480 * t452 + (t335 + t478) * t280 / 0.2e1;
t297 = -t2 * qJD(1) - t19 * qJD(2) + t8 * qJD(3);
t103 = (-t277 / 0.2e1 - t278 / 0.2e1 + 0.1e1 / 0.2e1) * t442;
t15 = t387 + t322 * t470 + (t349 * t470 - t292) * t282;
t57 = t500 - t506;
t293 = -qJD(1) * t15 + qJD(2) * t103 + qJD(3) * t57;
t199 = t464 * t402;
t198 = t464 * t403;
t197 = (qJD(1) * t402 + t394) * m(7);
t104 = (t395 + 0.1e1) * t378;
t85 = m(7) * t222 + (m(6) * pkin(8) + t432) * t282;
t71 = t469 * t512 + t379;
t42 = m(7) * t291 + t380;
t21 = t255 + t316 - t384 - t481;
t20 = t283 * t379 + (-mrSges(5,1) / 0.2e1 + t376) * t403 + (-pkin(4) * t403 + t254) * t471 + (-t283 * t511 + t254) * t469 + (mrSges(6,3) / 0.2e1 + t467 - mrSges(5,2) / 0.2e1) * t399 + (t338 - t498 + t492) * t450;
t16 = t209 * t453 + (t282 * t349 + t322) * t469 + t311 * t457 + t387 + (t372 - t422 / 0.2e1) * t283;
t11 = -t417 / 0.2e1 - t418 / 0.2e1 + t281 * t340 + t287 + t315;
t4 = -qJD(3) * t6 - qJD(4) * t7;
t1 = (Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(7,3) / 0.2e1) * t281 - t355 * qJ(5) + t285 + t286 + (t280 * t344 - t282 * t345) * t283;
t10 = [qJD(3) * t3 + qJD(4) * t5 - qJD(5) * t13 + qJD(6) * t25, t4, t1 * qJD(4) + t11 * qJD(5) + t16 * qJD(6) + (t93 * t447 / 0.2e1 + (t196 * t77 + t217 * t53 + t222 * t60) * t469) * t475 + t288 * t394 + (Ifges(4,5) + (t226 / 0.2e1 + t227 / 0.2e1 - t334 / 0.2e1 - pkin(3) * mrSges(5,2) + t196 * mrSges(7,2) - t216 * mrSges(6,3)) * t282 + (-t327 / 0.2e1 - t329 / 0.2e1 - t225 / 0.2e1 - pkin(3) * mrSges(5,1) - t196 * mrSges(7,1) + t216 * mrSges(6,1)) * t280 + (-m(5) * pkin(3) + t412) * t262) * t393 + t325 + (mrSges(4,2) * t406 - Ifges(4,6) * t281 + t222 * t210 + t217 * t211 - t93 * t337 + t77 * t219 + (-t60 * mrSges(7,3) - t483) * t282 + ((-t281 * mrSges(5,2) + t213) * t282 + (-mrSges(5,1) * t281 + t212) * t280 + t308 + m(6) * t323) * pkin(8)) * qJD(3), t1 * qJD(3) + (-t59 * mrSges(7,1) + t58 * mrSges(7,2) + t257 + t505 * t79 - (-mrSges(5,2) + mrSges(6,3)) * t397) * qJD(4) + t21 * qJD(5) + t42 * qJD(6) + ((-pkin(4) * t79 - qJ(5) * t397) * t471 + (qJ(5) * t58 - t463 * t59) * t469) * t474 + (t280 * t302 + t282 * t314) * t392 + t324, qJD(3) * t11 + qJD(4) * t21 - t410, qJD(3) * t16 + qJD(4) * t42 + t409; t4, t29 * qJD(3), t20 * qJD(4) + t198 * qJD(5) + t104 * qJD(6) + (t412 + t499) * qJD(3) * t281 + (m(5) * (t396 - t441) / 0.2e1 + (t396 + t408) * t471 + t196 * t378) * t475 + (t506 - mrSges(4,2) + t395 * (mrSges(5,3) + t432)) * t393 + t321, -t413 + t20 * qJD(3) + t199 * qJD(5) + (t192 * t472 - t445 / 0.2e1) * t474 + ((-mrSges(7,1) + t505) * t282 + (mrSges(5,2) + t431) * t280) * t392, qJD(3) * t198 + qJD(4) * t199, t104 * qJD(3); -qJD(4) * t2 + qJD(5) * t9 - qJD(6) * t15 - t325, -qJD(4) * t19 + qJD(6) * t103 - t321, qJD(4) * t8 - qJD(5) * t36 + qJD(6) * t57, t85 * qJD(5) + t71 * qJD(6) + t297 + (m(7) * (-qJ(5) * t217 - t222 * t463) - t222 * mrSges(7,1) - t217 * mrSges(7,2) - t302 * t282 + (Ifges(6,6) + t314) * t280 + (-m(6) * t326 - t337 - t339) * pkin(8)) * qJD(4), qJD(4) * t85 + t320, qJD(4) * t71 + t293; qJD(3) * t2 - qJD(5) * t23 + qJD(6) * t35 - t324, qJD(3) * t19 + t413, qJD(6) * t49 - t297, t234 * qJD(5), -t317, t319; -qJD(3) * t9 + qJD(4) * t23 - qJD(6) * t388 + t410, 0, -qJD(6) * t443 - t320, t317, 0, -t197; qJD(3) * t15 - qJD(4) * t35 + qJD(5) * t388 - t409, -t103 * qJD(3), -qJD(4) * t49 + qJD(5) * t443 - t293, -t319, t197, 0;];
Cq  = t10;
