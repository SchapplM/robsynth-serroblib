% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-03-09 08:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPPRP1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:25:32
% EndTime: 2019-03-09 08:25:47
% DurationCPUTime: 7.76s
% Computational Cost: add. (20123->553), mult. (40303->743), div. (0->0), fcn. (45824->8), ass. (0->279)
t464 = sin(qJ(2));
t408 = t464 * pkin(2);
t525 = m(4) * t408;
t513 = Ifges(7,4) + Ifges(6,5);
t524 = -Ifges(6,6) + Ifges(7,6);
t432 = sin(pkin(9));
t433 = cos(pkin(9));
t466 = cos(qJ(2));
t310 = -t432 * t466 - t433 * t464;
t333 = sin(pkin(10));
t334 = cos(pkin(10));
t463 = sin(qJ(5));
t465 = cos(qJ(5));
t349 = t333 * t465 + t334 * t463;
t507 = t349 * t310;
t431 = qJ(6) * t507;
t311 = t463 * t333 - t465 * t334;
t345 = t311 * t310;
t462 = pkin(5) * t345;
t111 = t462 - t431;
t308 = t432 * t464 - t433 * t466;
t523 = -t310 * mrSges(4,1) - t308 * mrSges(4,2);
t481 = -t507 / 0.2e1;
t522 = mrSges(6,3) + mrSges(7,2);
t328 = -pkin(2) * t466 - pkin(1);
t243 = t308 * pkin(3) + t310 * qJ(4) + t328;
t407 = t464 * pkin(7);
t318 = -qJ(3) * t464 - t407;
t409 = t466 * pkin(7);
t319 = qJ(3) * t466 + t409;
t506 = t432 * t318 + t433 * t319;
t150 = t333 * t243 + t334 * t506;
t422 = t310 * t333;
t109 = pkin(8) * t422 + t150;
t149 = t334 * t243 - t333 * t506;
t421 = t310 * t334;
t82 = t308 * pkin(4) + pkin(8) * t421 + t149;
t57 = -t109 * t463 + t465 * t82;
t49 = -t308 * pkin(5) - t57;
t521 = t49 + t57;
t406 = mrSges(6,3) / 0.2e1 + mrSges(7,2) / 0.2e1;
t520 = t406 * t507;
t472 = -t311 / 0.2e1;
t519 = t472 * t507;
t225 = t349 * t308;
t158 = t310 * mrSges(6,2) + t225 * mrSges(6,3);
t165 = t225 * mrSges(7,2) - mrSges(7,3) * t310;
t518 = t158 + t165;
t509 = t507 * mrSges(6,3);
t159 = -mrSges(6,2) * t308 + t509;
t297 = t308 * mrSges(7,3);
t164 = mrSges(7,2) * t507 + t297;
t417 = t159 + t164;
t397 = t433 * pkin(2);
t327 = -t397 - pkin(3);
t317 = -t334 * pkin(4) + t327;
t364 = t311 * pkin(5) - qJ(6) * t349;
t223 = t317 + t364;
t258 = mrSges(7,1) * t311 - mrSges(7,3) * t349;
t517 = m(7) * t223 + t258;
t516 = m(7) + m(6);
t484 = t345 / 0.2e1;
t514 = mrSges(7,1) + mrSges(6,1);
t512 = Ifges(7,2) + Ifges(6,3);
t332 = t334 ^ 2;
t411 = t333 ^ 2 + t332;
t511 = mrSges(5,3) * t411;
t510 = mrSges(6,3) * t345;
t425 = t308 * qJ(6);
t399 = t465 * t109;
t401 = t463 * t82;
t58 = t399 + t401;
t48 = t58 + t425;
t508 = t345 * t406;
t162 = mrSges(6,1) * t308 - t510;
t163 = -mrSges(7,1) * t308 + mrSges(7,2) * t345;
t416 = -t162 + t163;
t502 = -t311 * t513 + t349 * t524;
t501 = t345 * t524 + t507 * t513;
t491 = m(7) / 0.2e1;
t493 = m(6) / 0.2e1;
t500 = t491 + t493;
t329 = m(7) * qJ(6) + mrSges(7,3);
t499 = -m(7) * pkin(5) - t514;
t498 = -mrSges(6,2) + t329;
t228 = t308 * t311;
t404 = Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1;
t405 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t497 = -t404 * t225 - t405 * t228;
t443 = t228 * mrSges(7,3);
t445 = t228 * mrSges(6,2);
t446 = t225 * mrSges(7,1);
t447 = t225 * mrSges(6,1);
t496 = t443 / 0.2e1 - t445 / 0.2e1 + t446 / 0.2e1 + t447 / 0.2e1;
t495 = 0.2e1 * m(7);
t494 = m(5) / 0.2e1;
t492 = -m(7) / 0.2e1;
t490 = m(4) * pkin(2);
t248 = -t310 * pkin(3) + t308 * qJ(4) + t408;
t273 = -t433 * t318 + t319 * t432;
t153 = t333 * t248 - t273 * t334;
t424 = t308 * t333;
t110 = pkin(8) * t424 + t153;
t152 = t334 * t248 + t273 * t333;
t423 = t308 * t334;
t83 = -t310 * pkin(4) + pkin(8) * t423 + t152;
t59 = -t110 * t463 + t465 * t83;
t53 = t310 * pkin(5) - t59;
t488 = t53 / 0.2e1;
t441 = t310 * mrSges(7,1);
t444 = t228 * mrSges(7,2);
t161 = t441 + t444;
t487 = t161 / 0.2e1;
t486 = -t165 / 0.2e1;
t483 = t225 / 0.2e1;
t482 = -t225 / 0.2e1;
t480 = t228 / 0.2e1;
t395 = t432 * pkin(2);
t373 = t395 + qJ(4);
t356 = pkin(8) + t373;
t296 = t356 * t334;
t350 = t333 * t356;
t239 = t296 * t463 + t350 * t465;
t479 = -t239 / 0.2e1;
t240 = t296 * t465 - t350 * t463;
t478 = -t240 / 0.2e1;
t420 = t311 * qJ(6);
t461 = t349 * pkin(5);
t365 = -t420 - t461;
t477 = -t365 / 0.2e1;
t476 = t258 / 0.2e1;
t473 = -t310 / 0.2e1;
t471 = t311 / 0.2e1;
t470 = t349 / 0.2e1;
t469 = -t349 / 0.2e1;
t468 = t333 / 0.2e1;
t467 = t334 / 0.2e1;
t245 = m(7) * t349;
t458 = -t48 + t58;
t456 = mrSges(4,3) * t310;
t455 = Ifges(5,4) * t333;
t454 = Ifges(5,4) * t334;
t453 = Ifges(6,4) * t345;
t452 = Ifges(6,4) * t349;
t451 = Ifges(7,5) * t507;
t450 = Ifges(7,5) * t311;
t160 = -t310 * mrSges(6,1) - t228 * mrSges(6,3);
t342 = t411 * t373;
t358 = -t225 * t239 + t228 * t240;
t368 = t311 * mrSges(6,1) + mrSges(6,2) * t349;
t370 = -t334 * mrSges(5,1) + t333 * mrSges(5,2);
t336 = (-t308 * t342 - t327 * t310) * t494 + (-t310 * t317 + t358) * t493 + (-t223 * t310 + t358) * t491 + (-t308 * t432 + t310 * t433) * t490 / 0.2e1 - t308 * t511 / 0.2e1 + (t258 + t370 + t368) * t473 + t522 * (-t225 * t470 + t228 * t472);
t249 = t310 * mrSges(5,2) + mrSges(5,3) * t424;
t251 = -t310 * mrSges(5,1) + mrSges(5,3) * t423;
t60 = t465 * t110 + t463 * t83;
t52 = -qJ(6) * t310 + t60;
t338 = (t152 * t334 + t333 * t153) * t494 + (-t311 * t59 + t349 * t60) * t493 + (t311 * t53 + t349 * t52) * t491 + t249 * t468 + t251 * t467 + t525 / 0.2e1;
t9 = t160 * t472 + t161 * t471 + t470 * t518 - t336 + t338 + t523;
t449 = qJD(1) * t9;
t114 = -t443 - t446;
t115 = t445 - t447;
t116 = -mrSges(7,1) * t507 - mrSges(7,3) * t345;
t117 = -mrSges(6,1) * t507 + mrSges(6,2) * t345;
t438 = t333 * Ifges(5,2);
t181 = -Ifges(5,6) * t310 + (t438 - t454) * t308;
t182 = -Ifges(5,5) * t310 + (-Ifges(5,1) * t334 + t455) * t308;
t188 = -pkin(4) * t424 + t506;
t189 = -pkin(4) * t422 + t273;
t436 = t334 * mrSges(5,2);
t439 = t333 * mrSges(5,1);
t369 = -t436 - t439;
t241 = t369 * t308;
t242 = t369 * t310;
t250 = -t308 * mrSges(5,2) + mrSges(5,3) * t422;
t252 = t308 * mrSges(5,1) + mrSges(5,3) * t421;
t95 = Ifges(7,1) * t345 + t308 * Ifges(7,4) - t451;
t221 = Ifges(6,4) * t507;
t97 = Ifges(6,1) * t345 + t308 * Ifges(6,5) + t221;
t402 = t95 / 0.2e1 + t97 / 0.2e1;
t218 = Ifges(7,5) * t345;
t91 = t308 * Ifges(7,6) - Ifges(7,3) * t507 + t218;
t93 = Ifges(6,2) * t507 + t308 * Ifges(6,6) + t453;
t403 = t91 / 0.2e1 - t93 / 0.2e1;
t435 = t334 * Ifges(5,5);
t437 = t333 * Ifges(5,6);
t367 = pkin(5) * t225 + t228 * qJ(6);
t67 = t188 - t367;
t366 = -pkin(5) * t507 - qJ(6) * t345;
t68 = t189 + t366;
t90 = Ifges(7,5) * t228 - Ifges(7,6) * t310 - Ifges(7,3) * t225;
t92 = Ifges(6,4) * t228 + Ifges(6,2) * t225 - Ifges(6,6) * t310;
t94 = Ifges(7,1) * t228 - Ifges(7,4) * t310 - Ifges(7,5) * t225;
t96 = Ifges(6,1) * t228 + Ifges(6,4) * t225 - Ifges(6,5) * t310;
t1 = m(6) * (t188 * t189 + t57 * t59 + t58 * t60) + m(7) * (t48 * t52 + t49 * t53 + t67 * t68) + m(5) * (t149 * t152 + t150 * t153 + t273 * t506) + (t523 + t525) * t328 + (Ifges(3,1) - Ifges(3,2)) * t466 * t464 - (-t92 / 0.2e1 + t90 / 0.2e1) * t507 + t506 * t242 + (t181 * t468 - t334 * t182 / 0.2e1 - mrSges(4,2) * t408 + (t435 / 0.2e1 - t437 / 0.2e1 - Ifges(4,4)) * t310 - t405 * t345 - t404 * t507) * t310 - pkin(1) * (mrSges(3,1) * t464 + mrSges(3,2) * t466) + t402 * t228 - t403 * t225 + t273 * t241 + t150 * t249 + t153 * t250 + t149 * t251 + t152 * t252 + t188 * t117 + t189 * t115 + t58 * t158 + t60 * t159 + t57 * t160 + t49 * t161 + t59 * t162 + t53 * t163 + t52 * t164 + t48 * t165 + t67 * t116 + t68 * t114 + (-t464 ^ 2 + t466 ^ 2) * Ifges(3,4) + (mrSges(4,1) * t408 + (t332 * Ifges(5,1) / 0.2e1 - Ifges(5,3) + Ifges(4,1) - Ifges(4,2) + (-t454 + t438 / 0.2e1) * t333 - t512) * t310 + (Ifges(4,4) - t435 + t437) * t308 - t497) * t308 + (t94 / 0.2e1 + t96 / 0.2e1) * t345;
t448 = t1 * qJD(1);
t442 = t308 * mrSges(4,3);
t440 = t311 * mrSges(7,2);
t112 = mrSges(7,1) * t345 - mrSges(7,3) * t507;
t113 = mrSges(6,1) * t345 + mrSges(6,2) * t507;
t118 = Ifges(7,3) * t345 + t451;
t119 = -Ifges(6,2) * t345 + t221;
t120 = Ifges(7,1) * t507 + t218;
t121 = Ifges(6,1) * t507 - t453;
t4 = t68 * t112 + t189 * t113 + (-t48 * mrSges(7,2) + t120 / 0.2e1 + t121 / 0.2e1 + t403) * t345 - (-t119 / 0.2e1 - t49 * mrSges(7,2) + t118 / 0.2e1 - t402) * t507 + (m(7) * t68 + t116) * t111 + (m(7) * t49 + t416 - t510) * t58 + (m(7) * t48 + t417 - t509) * t57 + t501 * t308 / 0.2e1;
t434 = t4 * qJD(1);
t12 = t417 * t507 + t416 * t345 + m(7) * (t345 * t49 + t48 * t507) + m(6) * (-t345 * t57 + t507 * t58) + (t334 * t252 + t333 * t250 + m(5) * (t149 * t334 + t150 * t333)) * t310;
t430 = qJD(1) * t12;
t21 = t308 * t164 + m(7) * (t308 * t48 - t345 * t68) - t345 * t116;
t429 = qJD(1) * t21;
t360 = t333 * t149 - t150 * t334;
t418 = t334 * t250;
t419 = t333 * t252;
t426 = t273 * t310;
t11 = t417 * t228 - t416 * t225 + (-t418 + t419 + t442) * t308 + (-t116 - t117 - t242 + t456) * t310 + m(7) * (-t225 * t49 + t228 * t48 - t310 * t68) + m(6) * (-t189 * t310 + t225 * t57 + t228 * t58) + m(5) * (t308 * t360 - t426) + m(4) * (-t308 * t506 - t426);
t428 = t11 * qJD(1);
t388 = t308 * t470;
t130 = (t388 + t483) * m(7);
t410 = t130 * qJD(1);
t387 = -t162 / 0.2e1 + t163 / 0.2e1;
t386 = -t164 / 0.2e1 - t159 / 0.2e1;
t383 = m(5) * t411;
t256 = mrSges(7,1) * t349 + t311 * mrSges(7,3);
t257 = mrSges(6,1) * t349 - t311 * mrSges(6,2);
t375 = t162 * t469 + t163 * t470 + t417 * t472;
t341 = (-t223 * t345 + t240 * t308 - t349 * t68) * t492 + t345 * t476 + t116 * t470;
t354 = m(7) * t488 + t441 / 0.2e1;
t18 = t341 + t354 + t444;
t84 = t517 * t349;
t363 = qJD(1) * t18 + qJD(2) * t84;
t26 = (t111 / 0.4e1 + t462 / 0.4e1 - t431 / 0.4e1) * t495 + 0.2e1 * (mrSges(7,3) - mrSges(6,2)) * t481 + 0.2e1 * t514 * t484;
t351 = -t256 - t257;
t66 = (-t365 / 0.4e1 + t461 / 0.4e1 + t420 / 0.4e1) * t495 - t351;
t362 = qJD(1) * t26 + qJD(2) * t66;
t343 = t500 * (t311 * t345 + t349 * t507);
t34 = 0.2e1 * (t383 / 0.4e1 + m(6) / 0.4e1 + m(7) / 0.4e1 + m(5) / 0.4e1) * t310 + t343;
t361 = qJD(1) * t34;
t359 = t239 * t345 + t240 * t507;
t23 = -t297 + (t58 / 0.4e1 - t399 / 0.4e1 - t401 / 0.4e1 - t425 / 0.2e1) * t495;
t357 = qJD(1) * t23 - qJD(5) * t329;
t104 = t495 * t484;
t355 = qJD(1) * t104 + qJD(2) * t245;
t353 = t333 * t373;
t352 = t334 * t373;
t259 = Ifges(7,3) * t349 - t450;
t303 = Ifges(7,5) * t349;
t260 = Ifges(7,3) * t311 + t303;
t306 = Ifges(6,4) * t311;
t261 = -Ifges(6,2) * t349 - t306;
t262 = -Ifges(6,2) * t311 + t452;
t263 = -Ifges(7,1) * t311 + t303;
t264 = Ifges(7,1) * t349 + t450;
t265 = -Ifges(6,1) * t311 - t452;
t266 = Ifges(6,1) * t349 - t306;
t14 = t223 * t256 + t317 * t257 - (-t263 / 0.2e1 + t262 / 0.2e1 - t260 / 0.2e1 - t265 / 0.2e1) * t349 + (-t264 / 0.2e1 - t261 / 0.2e1 + t259 / 0.2e1 - t266 / 0.2e1) * t311 - t517 * t365;
t335 = t386 * t239 + (-t119 / 0.4e1 - t97 / 0.4e1 - t95 / 0.4e1 + t118 / 0.4e1) * t311 - (-t121 / 0.4e1 - t120 / 0.4e1 - t91 / 0.4e1 + t93 / 0.4e1) * t349 + (t345 * t478 - t507 * t479 - (-t58 / 0.2e1 + t48 / 0.2e1) * t349 + (-t49 / 0.2e1 - t57 / 0.2e1) * t311) * mrSges(7,2) - (mrSges(6,3) * t479 - t264 / 0.4e1 - t266 / 0.4e1 + t259 / 0.4e1 - t261 / 0.4e1) * t507 + (mrSges(6,3) * t478 + t263 / 0.4e1 + t265 / 0.4e1 + t260 / 0.4e1 - t262 / 0.4e1) * t345 + (t111 * t223 + t239 * t458 - t365 * t68) * t491 + t111 * t476 + t189 * t257 / 0.2e1 + t223 * t112 / 0.2e1 + t116 * t477 + t317 * t113 / 0.2e1 + t68 * t256 / 0.2e1 + t502 * t308 / 0.4e1 + (t491 * t521 + t387) * t240;
t339 = (-pkin(5) * t53 + qJ(6) * t52) * t492 + pkin(5) * t487 + qJ(6) * t486 - t52 * mrSges(7,3) / 0.2e1 + mrSges(7,1) * t488 - t59 * mrSges(6,1) / 0.2e1 + t60 * mrSges(6,2) / 0.2e1;
t2 = (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t310 + t335 + t339 + t497;
t348 = t2 * qJD(1) + t14 * qJD(2);
t344 = t367 * t491 + t496;
t5 = -(-t492 * t521 + t387 - t508) * t349 + (t458 * t492 - t386 - t520) * t311 + t344;
t347 = t5 * qJD(1);
t32 = m(5) * t342 + t511 + (t311 ^ 2 + t349 ^ 2) * t522 + t516 * (t239 * t349 - t311 * t240);
t337 = t360 * t494 - m(6) * (-t311 * t58 - t349 * t57 + t359) / 0.2e1 + (-t311 * t48 + t349 * t49 + t359) * t492 + t419 / 0.2e1 - t418 / 0.2e1;
t340 = (-t436 / 0.2e1 - t439 / 0.2e1) * t308 + t506 * t494 + t188 * t493 + t67 * t491 - t496;
t7 = -(t387 + t508) * t349 + (-t386 + t520) * t311 + t337 + t340;
t346 = -qJD(1) * t7 + qJD(2) * t32;
t151 = m(7) * t477 + t365 * t491;
t129 = (t388 + t482) * m(7);
t128 = m(7) * t240 - t440;
t105 = m(7) * t484 + t345 * t492;
t33 = t310 * t383 / 0.2e1 + t343 + (m(5) + t516) * t473;
t25 = (mrSges(7,1) / 0.2e1 + mrSges(6,1) / 0.2e1 - t514 / 0.2e1) * t345 + (t481 + t507 / 0.2e1) * mrSges(6,2);
t24 = 0.2e1 * t48 * t491 + t164;
t19 = -t341 + t354;
t10 = (-t160 / 0.2e1 + t487) * t311 - (t486 - t158 / 0.2e1) * t349 + t336 + t338;
t8 = -t337 + t340 + t375 + t522 * (t345 * t470 + t519);
t6 = (t311 * t458 + t349 * t521) * t491 + t344 + t375 + t522 * (t345 * t469 - t519);
t3 = Ifges(6,6) * t483 + Ifges(7,6) * t482 + t473 * t512 + t480 * t513 + t335 - t339;
t13 = [qJD(2) * t1 + qJD(3) * t11 + qJD(4) * t12 + qJD(5) * t4 + qJD(6) * t21, t10 * qJD(3) + t8 * qJD(4) + t3 * qJD(5) + t19 * qJD(6) + t448 + (t53 * t349 * mrSges(7,2) + (m(6) * t60 + m(7) * t52 + t518) * t240 + (t266 + t264) * t480 + (t96 + t94) * t470 + (-m(6) * t59 + m(7) * t53 - t160 + t161) * t239 + (Ifges(5,5) * t333 + Ifges(5,6) * t334 + t311 * t524 + t513 * t349) * t473 + (m(7) * t67 + t114) * t223 + (m(6) * t317 + t368) * t188 + t260 * t482 + t262 * t483 + t182 * t468 + t90 * t471 + t92 * t472 + t395 * t456 + t181 * t467 + t397 * t442 + (-t152 * t333 + t153 * t334) * mrSges(5,3) + m(5) * (-t152 * t353 + t153 * t352) + mrSges(3,2) * t407 + t249 * t352 + (m(5) * t327 - t433 * t490 - mrSges(4,1) + t370) * t506 + (-t432 * t490 + mrSges(4,2)) * t273 + (-t311 * t60 - t349 * t59) * mrSges(6,3) + Ifges(3,5) * t466 - Ifges(3,6) * t464 - (Ifges(5,1) * t333 + t454) * t423 / 0.2e1 + (Ifges(5,2) * t334 + t455) * t424 / 0.2e1 - t52 * t440 - mrSges(3,1) * t409 + t317 * t115 + t327 * t241 - Ifges(4,5) * t308 + Ifges(4,6) * t310 + t67 * t258 - t251 * t353) * qJD(2), t10 * qJD(2) + t33 * qJD(4) + t6 * qJD(5) + t129 * qJD(6) + t428 + 0.2e1 * t500 * qJD(3) * (-t225 * t311 + t228 * t349) qJD(2) * t8 + qJD(3) * t33 + qJD(5) * t25 + qJD(6) * t105 + t430, t434 + t3 * qJD(2) + t6 * qJD(3) + t25 * qJD(4) + (t366 * mrSges(7,2) + t498 * t57 + t499 * t58 + t501) * qJD(5) + t24 * qJD(6), qJD(2) * t19 + qJD(3) * t129 + qJD(4) * t105 + qJD(5) * t24 + t429; -qJD(3) * t9 - qJD(4) * t7 + qJD(5) * t2 - qJD(6) * t18 - t448, qJD(4) * t32 + qJD(5) * t14 - qJD(6) * t84, -t449, qJD(5) * t151 + t346, t151 * qJD(4) + (t364 * mrSges(7,2) - t239 * t498 + t240 * t499 + t502) * qJD(5) + t128 * qJD(6) + t348, qJD(5) * t128 - t363; qJD(2) * t9 + qJD(4) * t34 - qJD(5) * t5 + qJD(6) * t130 - t428, t449, 0, t361, t351 * qJD(5) + (t365 * qJD(5) / 0.2e1 + qJD(6) * t470) * t495 - t347, qJD(5) * t245 + t410; qJD(2) * t7 - qJD(3) * t34 + qJD(5) * t26 - qJD(6) * t104 - t430, qJD(5) * t66 - qJD(6) * t245 - t346, -t361, 0, t362, -t355; -qJD(2) * t2 + qJD(3) * t5 - qJD(4) * t26 - qJD(6) * t23 - t434, -qJD(4) * t66 - t348, t347, -t362, t329 * qJD(6), -t357; qJD(2) * t18 - qJD(3) * t130 + qJD(4) * t104 + qJD(5) * t23 - t429, qJD(4) * t245 + t363, -t410, t355, t357, 0;];
Cq  = t13;
