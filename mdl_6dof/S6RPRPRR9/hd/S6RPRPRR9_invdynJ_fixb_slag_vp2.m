% Calculate vector of inverse dynamics joint torques for
% S6RPRPRR9
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRR9_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR9_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR9_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR9_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_invdynJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR9_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR9_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR9_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:01:18
% EndTime: 2019-03-09 04:02:19
% DurationCPUTime: 36.03s
% Computational Cost: add. (36590->950), mult. (122895->1328), div. (0->0), fcn. (107413->16), ass. (0->424)
t353 = sin(pkin(12));
t355 = sin(pkin(6));
t357 = cos(pkin(12));
t362 = sin(qJ(3));
t358 = cos(pkin(7));
t366 = cos(qJ(3));
t466 = t358 * t366;
t375 = (-t353 * t466 - t357 * t362) * t355;
t286 = qJD(1) * t375;
t454 = qJD(1) * t355;
t467 = t358 * t362;
t468 = t357 * t366;
t287 = (-t353 * t467 + t468) * t454;
t352 = sin(pkin(13));
t356 = cos(pkin(13));
t459 = t366 * t356;
t323 = t352 * t362 - t459;
t354 = sin(pkin(7));
t590 = -t323 * t354 * qJD(3) - t286 * t352 - t287 * t356;
t359 = cos(pkin(6));
t391 = -t353 * t362 + t357 * t466;
t474 = t354 * t366;
t266 = t355 * t391 + t359 * t474;
t253 = t266 * qJD(1);
t390 = t353 * t366 + t357 * t467;
t475 = t354 * t362;
t267 = t355 * t390 + t359 * t475;
t256 = t267 * qJD(1);
t423 = t356 * t253 - t256 * t352;
t565 = t423 - qJD(5);
t526 = pkin(3) * t352;
t348 = pkin(10) + t526;
t361 = sin(qJ(5));
t451 = qJD(5) * t361;
t398 = t253 * t352 + t356 * t256;
t527 = pkin(3) * t256;
t115 = pkin(4) * t398 - pkin(10) * t423 + t527;
t365 = cos(qJ(5));
t473 = t355 * t357;
t341 = qJ(2) * t473;
t529 = pkin(1) * t359;
t443 = qJD(1) * t529;
t307 = qJD(1) * t341 + t353 * t443;
t472 = t355 * t358;
t476 = t354 * t359;
t378 = (t357 * t472 + t476) * pkin(9);
t240 = qJD(1) * t378 + t307;
t339 = t357 * t443;
t479 = t353 * t355;
t521 = pkin(9) * t358;
t528 = pkin(2) * t359;
t374 = t528 + (-qJ(2) - t521) * t479;
t249 = qJD(1) * t374 + t339;
t522 = pkin(9) * t353;
t294 = (-pkin(2) * t357 - t354 * t522 - pkin(1)) * t355;
t280 = qJD(1) * t294 + qJD(2);
t166 = t362 * (t249 * t358 + t280 * t354) + t240 * t366;
t145 = qJ(4) * t253 + t166;
t133 = t352 * t145;
t489 = t240 * t362;
t165 = t249 * t466 + t280 * t474 - t489;
t144 = -qJ(4) * t256 + t165;
t80 = t144 * t356 - t133;
t53 = t361 * t115 + t365 * t80;
t603 = -t348 * t451 - t53;
t395 = t352 * t366 + t356 * t362;
t303 = t395 * t354;
t276 = t303 * t361 - t365 * t358;
t434 = t353 * t454;
t420 = t354 * t434;
t593 = -qJD(5) * t276 - t361 * t420 + t365 * t590;
t591 = qJD(3) * t303 + t356 * t286 - t287 * t352;
t545 = m(6) + m(7);
t602 = pkin(10) * t545 - mrSges(5,2);
t387 = pkin(3) * t266;
t469 = t356 * t145;
t79 = t144 * t352 + t469;
t601 = -t79 - t565 * (pkin(5) * t361 - pkin(11) * t365);
t600 = pkin(11) * t398 - t603;
t321 = (qJ(2) * qJDD(1) + qJD(1) * qJD(2)) * t355;
t441 = qJDD(1) * t529;
t283 = t357 * t321 + t353 * t441;
t236 = qJDD(1) * t378 + t283;
t282 = -t321 * t353 + t357 * t441;
t237 = (-t472 * t522 + t528) * qJDD(1) + t282;
t275 = qJDD(1) * t294 + qJDD(2);
t452 = qJD(3) * t366;
t431 = t358 * t452;
t432 = t354 * t452;
t100 = -qJD(3) * t489 + t366 * t236 + t237 * t467 + t249 * t431 + t275 * t475 + t280 * t432;
t101 = -qJD(3) * t166 - t236 * t362 + t237 * t466 + t275 * t474;
t446 = qJD(1) * qJD(3);
t210 = (qJDD(1) * t362 + t366 * t446) * t476 + (qJDD(1) * t390 + t391 * t446) * t355;
t211 = (qJDD(1) * t366 - t362 * t446) * t476 + (qJDD(1) * t391 - t390 * t446) * t355;
t146 = -t210 * t352 + t211 * t356;
t147 = t210 * t356 + t211 * t352;
t311 = -t354 * t473 + t358 * t359;
t299 = qJDD(1) * t311 + qJDD(3);
t69 = pkin(3) * t299 - qJ(4) * t210 - qJD(4) * t256 + t101;
t73 = qJ(4) * t211 + qJD(4) * t253 + t100;
t31 = -t352 * t73 + t356 * t69;
t32 = t352 * t69 + t356 * t73;
t599 = t101 * mrSges(4,1) + t31 * mrSges(5,1) - t100 * mrSges(4,2) - t32 * mrSges(5,2) + Ifges(4,5) * t210 + Ifges(5,5) * t147 + Ifges(4,6) * t211 + Ifges(5,6) * t146;
t300 = qJD(1) * t311 + qJD(3);
t171 = t300 * t361 + t365 * t398;
t360 = sin(qJ(6));
t364 = cos(qJ(6));
t117 = -t171 * t360 - t364 * t565;
t118 = t171 * t364 - t360 * t565;
t170 = t300 * t365 - t361 * t398;
t169 = qJD(6) - t170;
t55 = t118 * Ifges(7,5) + t117 * Ifges(7,6) + t169 * Ifges(7,3);
t498 = t171 * Ifges(6,4);
t597 = t565 * Ifges(6,6);
t598 = t170 * Ifges(6,2);
t90 = t498 - t597 + t598;
t587 = t55 / 0.2e1 - t90 / 0.2e1;
t505 = mrSges(6,3) * t171;
t120 = -mrSges(6,1) * t565 - t505;
t70 = -mrSges(7,1) * t117 + mrSges(7,2) * t118;
t576 = t120 - t70;
t277 = t303 * t365 + t358 * t361;
t302 = t352 * t475 - t354 * t459;
t228 = -t277 * t360 + t302 * t364;
t595 = qJD(6) * t228 + t360 * t591 + t364 * t593;
t229 = t277 * t364 + t302 * t360;
t594 = -qJD(6) * t229 - t360 * t593 + t364 * t591;
t458 = -mrSges(5,1) * t300 - mrSges(6,1) * t170 + mrSges(6,2) * t171 + mrSges(5,3) * t398;
t592 = qJD(5) * t277 + t361 * t590 + t365 * t420;
t463 = t360 * t365;
t131 = t364 * t398 - t423 * t463;
t450 = qJD(5) * t365;
t589 = t360 * t450 + t131;
t367 = cos(qJ(1));
t465 = t359 * t367;
t363 = sin(qJ(1));
t477 = t353 * t363;
t312 = -t357 * t465 + t477;
t461 = t363 * t357;
t313 = t353 * t465 + t461;
t484 = t313 * t362;
t588 = t312 * t466 + t484;
t30 = pkin(10) * t299 + t32;
t197 = -t237 * t354 + t358 * t275;
t138 = -pkin(3) * t211 + qJDD(4) + t197;
t60 = -pkin(4) * t146 - pkin(10) * t147 + t138;
t126 = pkin(3) * t300 + t144;
t77 = t352 * t126 + t469;
t75 = pkin(10) * t300 + t77;
t212 = -t249 * t354 + t358 * t280;
t167 = -pkin(3) * t253 + qJD(4) + t212;
t99 = -pkin(4) * t423 - pkin(10) * t398 + t167;
t7 = t365 * t30 + t361 * t60 + t99 * t450 - t451 * t75;
t43 = t361 * t99 + t365 * t75;
t8 = -qJD(5) * t43 - t30 * t361 + t365 * t60;
t586 = t8 * mrSges(6,1) - t7 * mrSges(6,2);
t29 = -pkin(4) * t299 - t31;
t96 = qJD(5) * t170 + t147 * t365 + t299 * t361;
t97 = -qJD(5) * t171 - t147 * t361 + t299 * t365;
t12 = -pkin(5) * t97 - pkin(11) * t96 + t29;
t37 = -pkin(11) * t565 + t43;
t76 = t126 * t356 - t133;
t74 = -pkin(4) * t300 - t76;
t49 = -pkin(5) * t170 - pkin(11) * t171 + t74;
t13 = -t360 * t37 + t364 * t49;
t143 = qJDD(5) - t146;
t5 = pkin(11) * t143 + t7;
t1 = qJD(6) * t13 + t12 * t360 + t364 * t5;
t14 = t360 * t49 + t364 * t37;
t2 = -qJD(6) * t14 + t12 * t364 - t360 * t5;
t585 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t470 = t355 * t367;
t485 = t312 * t354;
t272 = t358 * t470 - t485;
t305 = t395 * t358;
t373 = t303 * t470 + t312 * t305 + t313 * t323;
t584 = -t272 * t365 + t361 * t373;
t583 = t272 * t361 + t365 * t373;
t582 = mrSges(6,1) * t74 + t13 * mrSges(7,1) - t14 * mrSges(7,2);
t540 = t143 / 0.2e1;
t546 = t97 / 0.2e1;
t547 = t96 / 0.2e1;
t557 = Ifges(6,1) * t547 + Ifges(6,4) * t546 + Ifges(6,5) * t540;
t40 = qJD(6) * t117 + t143 * t360 + t364 * t96;
t556 = t40 / 0.2e1;
t41 = -qJD(6) * t118 + t143 * t364 - t360 * t96;
t555 = t41 / 0.2e1;
t95 = qJDD(6) - t97;
t548 = t95 / 0.2e1;
t15 = -mrSges(7,1) * t41 + mrSges(7,2) * t40;
t61 = mrSges(6,1) * t143 - mrSges(6,3) * t96;
t511 = t15 - t61;
t581 = t398 * Ifges(5,1);
t580 = t398 * Ifges(5,4);
t579 = t423 * Ifges(5,4);
t578 = t423 * Ifges(5,2);
t417 = -mrSges(7,1) * t364 + mrSges(7,2) * t360;
t384 = m(7) * pkin(5) - t417;
t577 = mrSges(6,1) + t384;
t442 = m(7) * pkin(11) + mrSges(7,3);
t422 = mrSges(6,2) - t442;
t524 = pkin(3) * t356;
t349 = -pkin(4) - t524;
t320 = -pkin(5) * t365 - pkin(11) * t361 + t349;
t270 = t320 * t364 - t348 * t463;
t575 = qJD(6) * t270 + t601 * t360 - t600 * t364;
t460 = t364 * t365;
t271 = t320 * t360 + t348 * t460;
t574 = -qJD(6) * t271 + t600 * t360 + t601 * t364;
t345 = t357 * t529;
t268 = t345 + t374;
t219 = -t268 * t354 + t358 * t294;
t180 = t219 - t387;
t213 = -t356 * t266 + t267 * t352;
t214 = t266 * t352 + t267 * t356;
t104 = pkin(4) * t213 - pkin(10) * t214 + t180;
t317 = t353 * t529 + t341;
t261 = t378 + t317;
t172 = -t261 * t362 + t268 * t466 + t294 * t474;
t148 = pkin(3) * t311 - qJ(4) * t267 + t172;
t247 = t366 * t261;
t173 = t268 * t467 + t294 * t475 + t247;
t156 = qJ(4) * t266 + t173;
t88 = t352 * t148 + t356 * t156;
t83 = pkin(10) * t311 + t88;
t573 = t361 * t104 + t365 * t83;
t447 = qJD(6) * t364;
t572 = t361 * t447 + t589;
t132 = t360 * t398 + t423 * t460;
t448 = qJD(6) * t361;
t571 = t360 * t448 - t364 * t450 + t132;
t301 = t311 * t365;
t570 = -t214 * t361 + t301;
t492 = t423 * t361;
t569 = t451 - t492;
t24 = mrSges(7,1) * t95 - mrSges(7,3) * t40;
t25 = -mrSges(7,2) * t95 + mrSges(7,3) * t41;
t568 = -t360 * t24 + t364 * t25;
t567 = -t361 * t8 + t365 * t7;
t566 = t1 * t364 - t2 * t360;
t418 = mrSges(6,1) * t365 - mrSges(6,2) * t361;
t564 = -t361 * t442 - t365 * t384 - mrSges(5,1) - t418;
t218 = t355 * (t305 * t357 - t323 * t353) + t303 * t359;
t254 = t266 * qJD(3);
t255 = t267 * qJD(3);
t200 = t254 * t352 + t356 * t255;
t203 = t254 * t356 - t255 * t352;
t453 = qJD(2) * t355;
t433 = t353 * t453;
t329 = t354 * t433;
t235 = pkin(3) * t255 + t329;
t110 = pkin(4) * t200 - pkin(10) * t203 + t235;
t157 = t268 * t431 + t294 * t432 + t453 * t468 + (-qJD(3) * t261 - t358 * t433) * t362;
t112 = -qJ(4) * t255 + qJD(4) * t266 + t157;
t158 = qJD(2) * t375 + (-t247 + (-t268 * t358 - t294 * t354) * t362) * qJD(3);
t113 = -qJ(4) * t254 - qJD(4) * t267 + t158;
t67 = t112 * t356 + t113 * t352;
t19 = -qJD(5) * t573 + t110 * t365 - t361 * t67;
t416 = mrSges(7,1) * t360 + mrSges(7,2) * t364;
t563 = t416 + t602;
t539 = -t169 / 0.2e1;
t542 = -t118 / 0.2e1;
t544 = -t117 / 0.2e1;
t561 = Ifges(7,5) * t542 + Ifges(7,6) * t544 + Ifges(7,3) * t539;
t560 = 0.2e1 * t359;
t10 = t40 * Ifges(7,4) + t41 * Ifges(7,2) + t95 * Ifges(7,6);
t559 = t10 / 0.2e1;
t558 = Ifges(7,1) * t556 + Ifges(7,4) * t555 + Ifges(7,5) * t548;
t501 = Ifges(7,4) * t118;
t56 = Ifges(7,2) * t117 + Ifges(7,6) * t169 + t501;
t553 = -t56 / 0.2e1;
t552 = t56 / 0.2e1;
t116 = Ifges(7,4) * t117;
t57 = Ifges(7,1) * t118 + Ifges(7,5) * t169 + t116;
t551 = -t57 / 0.2e1;
t550 = t57 / 0.2e1;
t543 = t117 / 0.2e1;
t541 = t118 / 0.2e1;
t538 = t169 / 0.2e1;
t537 = -t170 / 0.2e1;
t536 = -t171 / 0.2e1;
t535 = t171 / 0.2e1;
t534 = t565 / 0.2e1;
t532 = t256 / 0.2e1;
t531 = -t300 / 0.2e1;
t523 = pkin(3) * t366;
t6 = -pkin(5) * t143 - t8;
t518 = t361 * t6;
t513 = -mrSges(4,3) - mrSges(5,3);
t512 = pkin(9) + qJ(4);
t510 = mrSges(3,1) * t357;
t509 = mrSges(3,2) * t353;
t508 = mrSges(4,3) * t253;
t507 = mrSges(4,3) * t256;
t506 = mrSges(6,3) * t170;
t504 = Ifges(4,4) * t256;
t503 = Ifges(6,4) * t361;
t502 = Ifges(6,4) * t365;
t500 = Ifges(7,4) * t360;
t499 = Ifges(7,4) * t364;
t494 = t170 * t360;
t493 = t170 * t364;
t491 = t423 * t365;
t486 = t311 * t361;
t314 = -t353 * t367 - t359 * t461;
t483 = t314 * t354;
t315 = t357 * t367 - t359 * t477;
t482 = t315 * t362;
t481 = (-mrSges(3,2) * t359 + mrSges(3,3) * t473) * qJD(1) * t357;
t471 = t355 * t363;
t464 = t360 * t361;
t462 = t361 * t364;
t455 = t367 * pkin(1) + qJ(2) * t471;
t449 = qJD(6) * t360;
t445 = qJDD(1) * t355;
t9 = Ifges(7,5) * t40 + Ifges(7,6) * t41 + Ifges(7,3) * t95;
t444 = Ifges(6,5) * t96 + Ifges(6,6) * t97 + Ifges(6,3) * t143;
t439 = t354 * t471;
t438 = t354 * t470;
t429 = t348 * t450;
t427 = t450 / 0.2e1;
t426 = -t448 / 0.2e1;
t425 = -t363 * pkin(1) + qJ(2) * t470;
t84 = -t146 * mrSges(5,1) + t147 * mrSges(5,2);
t424 = m(3) + m(4) + m(5) + t545;
t66 = t112 * t352 - t356 * t113;
t87 = t148 * t356 - t352 * t156;
t421 = -m(4) * t521 - mrSges(3,3);
t415 = Ifges(6,1) * t365 - t503;
t414 = Ifges(7,1) * t364 - t500;
t413 = Ifges(7,1) * t360 + t499;
t412 = -Ifges(6,2) * t361 + t502;
t411 = -Ifges(7,2) * t360 + t499;
t410 = Ifges(7,2) * t364 + t500;
t409 = Ifges(6,5) * t365 - Ifges(6,6) * t361;
t408 = Ifges(7,5) * t364 - Ifges(7,6) * t360;
t407 = Ifges(7,5) * t360 + Ifges(7,6) * t364;
t45 = pkin(11) * t213 + t573;
t182 = t214 * t365 + t486;
t82 = -pkin(4) * t311 - t87;
t54 = -pkin(5) * t570 - pkin(11) * t182 + t82;
t21 = t360 * t54 + t364 * t45;
t20 = -t360 * t45 + t364 * t54;
t42 = -t361 * t75 + t365 * t99;
t404 = mrSges(6,3) + t602;
t308 = pkin(3) * t475 + t358 * t512;
t309 = pkin(3) * t467 - t354 * t512;
t350 = pkin(2) + t523;
t403 = t308 * t471 + t314 * t309 + t315 * t350 + t455;
t50 = t104 * t365 - t361 * t83;
t52 = t115 * t365 - t361 * t80;
t130 = t182 * t364 + t213 * t360;
t129 = -t182 * t360 + t213 * t364;
t396 = -(-qJ(2) * t434 + t339) * t353 + t307 * t357;
t394 = t439 * t523 + (t314 * t466 - t482) * pkin(3);
t196 = t303 * t471 + t305 * t314 - t315 * t323;
t393 = t196 * pkin(4) + t403;
t392 = t314 * t358 + t439;
t36 = pkin(5) * t565 - t42;
t389 = t36 * t416;
t388 = t74 * (mrSges(6,1) * t361 + mrSges(6,2) * t365);
t385 = t308 * t470 + t312 * t309 - t313 * t350 + t425;
t18 = t104 * t450 + t361 * t110 + t365 * t67 - t451 * t83;
t274 = t358 * t471 - t483;
t382 = -g(1) * t274 + g(2) * t272 - g(3) * t311;
t381 = t312 * t467 - t313 * t366 + t362 * t438;
t304 = t323 * t358;
t190 = t302 * t470 + t312 * t304 - t313 * t395;
t371 = (-t366 * t438 - t588) * pkin(3);
t369 = (-t13 * t364 - t14 * t360) * qJD(6) + t566;
t340 = -pkin(1) * t445 + qJDD(2);
t332 = t445 * t509;
t318 = (mrSges(3,1) * t359 - mrSges(3,3) * t479) * qJD(1);
t316 = -qJ(2) * t479 + t345;
t293 = Ifges(4,3) * t299;
t292 = Ifges(5,3) * t299;
t248 = Ifges(4,4) * t253;
t231 = t315 * t366 + t362 * t392;
t230 = t366 * t392 - t482;
t223 = mrSges(4,1) * t300 - t507;
t222 = -mrSges(4,2) * t300 + t508;
t217 = -t302 * t359 + (-t304 * t357 - t353 * t395) * t355;
t209 = -mrSges(4,1) * t253 + mrSges(4,2) * t256;
t195 = -t302 * t471 - t304 * t314 - t315 * t395;
t188 = t218 * t365 + t486;
t179 = t256 * Ifges(4,1) + t300 * Ifges(4,5) + t248;
t178 = t253 * Ifges(4,2) + t300 * Ifges(4,6) + t504;
t177 = -mrSges(4,2) * t299 + mrSges(4,3) * t211;
t176 = mrSges(4,1) * t299 - mrSges(4,3) * t210;
t174 = -mrSges(5,2) * t300 + mrSges(5,3) * t423;
t168 = Ifges(6,4) * t170;
t164 = t196 * t365 + t274 * t361;
t163 = t196 * t361 - t274 * t365;
t149 = -mrSges(4,1) * t211 + mrSges(4,2) * t210;
t136 = -mrSges(5,1) * t423 + mrSges(5,2) * t398;
t128 = qJD(5) * t182 + t203 * t361;
t127 = qJD(5) * t570 + t203 * t365;
t125 = mrSges(5,1) * t299 - mrSges(5,3) * t147;
t124 = -mrSges(5,2) * t299 + mrSges(5,3) * t146;
t122 = t300 * Ifges(5,5) + t579 + t581;
t121 = t300 * Ifges(5,6) + t578 + t580;
t119 = mrSges(6,2) * t565 + t506;
t108 = pkin(5) * t171 - pkin(11) * t170;
t106 = t164 * t364 - t195 * t360;
t105 = -t164 * t360 - t195 * t364;
t91 = Ifges(6,1) * t171 - Ifges(6,5) * t565 + t168;
t89 = t171 * Ifges(6,5) + t170 * Ifges(6,6) - Ifges(6,3) * t565;
t86 = mrSges(7,1) * t169 - mrSges(7,3) * t118;
t85 = -mrSges(7,2) * t169 + mrSges(7,3) * t117;
t64 = -qJD(6) * t130 - t127 * t360 + t200 * t364;
t63 = qJD(6) * t129 + t127 * t364 + t200 * t360;
t62 = -mrSges(6,2) * t143 + mrSges(6,3) * t97;
t48 = -mrSges(6,1) * t97 + mrSges(6,2) * t96;
t46 = -pkin(5) * t398 - t52;
t44 = -pkin(5) * t213 - t50;
t34 = t96 * Ifges(6,4) + t97 * Ifges(6,2) + t143 * Ifges(6,6);
t33 = pkin(5) * t128 - pkin(11) * t127 + t66;
t27 = t108 * t360 + t364 * t42;
t26 = t108 * t364 - t360 * t42;
t17 = -pkin(5) * t200 - t19;
t16 = pkin(11) * t200 + t18;
t4 = -qJD(6) * t21 - t16 * t360 + t33 * t364;
t3 = qJD(6) * t20 + t16 * t364 + t33 * t360;
t11 = [(Ifges(4,5) * t254 + Ifges(5,5) * t203 - Ifges(4,6) * t255 - Ifges(5,6) * t200) * t300 / 0.2e1 + (-t165 * t254 - t166 * t255) * mrSges(4,3) + (Ifges(4,1) * t254 - Ifges(4,4) * t255) * t532 + t253 * (Ifges(4,4) * t254 - Ifges(4,2) * t255) / 0.2e1 + t212 * (mrSges(4,1) * t255 + mrSges(4,2) * t254) + m(6) * (t18 * t43 + t19 * t42 + t29 * t82 + t50 * t8 + t573 * t7 + t66 * t74) + t573 * t62 + (t1 * t129 - t13 * t63 - t130 * t2 + t14 * t64) * mrSges(7,3) + (t127 * t74 - t200 * t43) * mrSges(6,2) + (mrSges(6,2) * t29 - mrSges(6,3) * t8 + 0.2e1 * t557) * t182 + (Ifges(7,6) * t543 - Ifges(6,4) * t535 + Ifges(7,3) * t538 + Ifges(7,5) * t541 + t597 / 0.2e1 - t598 / 0.2e1 - t43 * mrSges(6,3) + t582 + t587) * t128 + (Ifges(7,1) * t130 + Ifges(7,4) * t129) * t556 + (Ifges(7,1) * t63 + Ifges(7,4) * t64) * t541 - (Ifges(7,6) * t555 + Ifges(7,5) * t556 - Ifges(6,2) * t546 - Ifges(6,4) * t547 + Ifges(7,3) * t548 - Ifges(6,6) * t540 - t7 * mrSges(6,3) + t9 / 0.2e1 - t34 / 0.2e1 + t29 * mrSges(6,1) + t585) * t570 - t565 * (Ifges(6,5) * t127 + Ifges(6,3) * t200) / 0.2e1 + (Ifges(6,1) * t127 + Ifges(6,5) * t200) * t535 + (-Ifges(5,4) * t147 - Ifges(5,2) * t146 + t444 / 0.2e1 - Ifges(5,6) * t299 - t32 * mrSges(5,3) + Ifges(6,6) * t546 + Ifges(6,5) * t547 + Ifges(6,3) * t540 + t138 * mrSges(5,1) + t586) * t213 + (Ifges(7,5) * t130 + Ifges(7,6) * t129) * t548 + (Ifges(7,5) * t63 + Ifges(7,6) * t64) * t538 + (t292 / 0.2e1 + t293 / 0.2e1 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t299 + t599) * t311 + (mrSges(4,2) * t197 - mrSges(4,3) * t101 + Ifges(4,1) * t210 + Ifges(4,4) * t211) * t267 + (Ifges(2,3) + (mrSges(3,1) * t316 - mrSges(3,2) * t317 + Ifges(3,3) * t359) * t359 + ((-mrSges(3,3) * t316 + Ifges(3,1) * t479 + Ifges(3,5) * t560) * t353 + (t317 * mrSges(3,3) + Ifges(3,6) * t560 + (mrSges(3,1) * pkin(1) + 0.2e1 * Ifges(3,4) * t353 + Ifges(3,2) * t357) * t355) * t357) * t355) * qJDD(1) + (-m(5) * t385 - t373 * mrSges(5,1) - m(3) * t425 + t313 * mrSges(3,1) - t312 * mrSges(3,2) + t363 * mrSges(2,1) - m(4) * (-t313 * pkin(2) - pkin(9) * t485 + t425) - t381 * mrSges(4,1) - t588 * mrSges(4,2) + (mrSges(2,2) + (-mrSges(4,2) * t474 + t421) * t355) * t367 + t513 * t272 - t577 * t583 - (t404 + t416) * t190 + t422 * t584 + t545 * (-pkin(4) * t373 - t385)) * g(1) + t170 * (Ifges(6,4) * t127 + Ifges(6,6) * t200) / 0.2e1 + (mrSges(3,1) * t282 - mrSges(3,2) * t283) * t359 + (-m(4) * (pkin(2) * t315 - pkin(9) * t483 + t455) - t231 * mrSges(4,1) - t230 * mrSges(4,2) - m(3) * t455 - t315 * mrSges(3,1) - t314 * mrSges(3,2) - mrSges(2,1) * t367 - m(5) * t403 - t196 * mrSges(5,1) - m(7) * (pkin(5) * t164 + t393) - t106 * mrSges(7,1) - t105 * mrSges(7,2) - m(6) * t393 - t164 * mrSges(6,1) + (t355 * t421 + mrSges(2,2)) * t363 + t513 * t274 + t404 * t195 + t422 * t163) * g(2) + (Ifges(7,4) * t130 + Ifges(7,2) * t129) * t555 + (Ifges(7,4) * t63 + Ifges(7,2) * t64) * t543 + (mrSges(5,2) * t138 - mrSges(5,3) * t31 + Ifges(5,1) * t147 + Ifges(5,4) * t146) * t214 + (-pkin(1) * t332 + t340 * (t509 - t510) + (-t282 * t353 + t283 * t357) * mrSges(3,3) + (t481 + (t209 * t354 - t318) * t353) * qJD(2)) * t355 + (-t200 * t77 - t203 * t76) * mrSges(5,3) + (-mrSges(4,1) * t197 + mrSges(4,3) * t100 + Ifges(4,4) * t210 + Ifges(4,2) * t211) * t266 + m(5) * (t138 * t180 + t167 * t235 + t31 * t87 + t32 * t88 - t66 * t76 + t67 * t77) + m(7) * (t1 * t21 + t13 * t4 + t14 * t3 + t17 * t36 + t2 * t20 + t44 * t6) + m(4) * (t100 * t173 + t101 * t172 + t157 * t166 + t158 * t165 + t197 * t219 + t212 * t329) + t63 * t550 + t64 * t552 + t130 * t558 + t129 * t559 + m(3) * (t282 * t316 + t283 * t317 + (-pkin(1) * t340 + qJD(2) * t396) * t355) + t254 * t179 / 0.2e1 - t255 * t178 / 0.2e1 + t235 * t136 + t158 * t223 + t219 * t149 + t157 * t222 + t167 * (mrSges(5,1) * t200 + mrSges(5,2) * t203) + t203 * t122 / 0.2e1 + t42 * (mrSges(6,1) * t200 - mrSges(6,3) * t127) + t200 * t89 / 0.2e1 - t200 * t121 / 0.2e1 + t180 * t84 + t67 * t174 + t172 * t176 + t173 * t177 + t6 * (-mrSges(7,1) * t129 + mrSges(7,2) * t130) + t88 * t124 + t87 * t125 + t127 * t91 / 0.2e1 + t18 * t119 + t19 * t120 + t4 * t86 + t3 * t85 + t82 * t48 + t458 * t66 + t423 * (Ifges(5,4) * t203 - Ifges(5,2) * t200) / 0.2e1 + t398 * (Ifges(5,1) * t203 - Ifges(5,4) * t200) / 0.2e1 + t20 * t24 + t21 * t25 + (Ifges(4,5) * t267 + Ifges(5,5) * t214 + Ifges(4,6) * t266) * t299 + t44 * t15 + t50 * t61 + t36 * (-mrSges(7,1) * t64 + mrSges(7,2) * t63) + t17 * t70; (t366 * t176 + t362 * t177 + (t222 * t366 - t223 * t362) * qJD(3)) * t354 + t594 * t86 + t595 * t85 + t332 + (-t125 + t48) * t302 + (t149 + t84) * t358 + t511 * t276 + t593 * t119 + t303 * t124 - t286 * t223 - t287 * t222 + t277 * t62 + t228 * t24 + t229 * t25 + t590 * t174 - t424 * t359 * g(3) + (-qJDD(1) * t510 + (-t481 + (t318 + (-t136 - t209) * t354) * t353) * qJD(1) + (-g(1) * t363 + g(2) * t367) * t424) * t355 - t576 * t592 + t458 * t591 + (t1 * t229 + t13 * t594 + t14 * t595 + t2 * t228 + t276 * t6 + t36 * t592) * m(7) + (-t276 * t8 + t277 * t7 + t29 * t302 - t42 * t592 + t43 * t593 + t591 * t74) * m(6) + (t138 * t358 - t167 * t420 - t302 * t31 + t303 * t32 + t590 * t77 - t591 * t76) * m(5) + (t197 * t358 + (t100 * t362 + t101 * t366 + (-t165 * t362 + t166 * t366) * qJD(3)) * t354 - t165 * t286 - t166 * t287 - t212 * t420) * m(4) + (-t396 * t454 + t340) * m(3); (-m(5) * t394 - mrSges(4,1) * t230 + mrSges(4,2) * t231 - t545 * (t195 * pkin(4) + t394) - t563 * t196 + t564 * t195) * g(1) + (t429 - t46) * t70 + (-t167 * t527 + t76 * t79 - t77 * t80 + (t31 * t356 + t32 * t352) * pkin(3)) * m(5) + (-t429 - t52) * t120 - (-Ifges(4,2) * t256 + t179 + t248) * t253 / 0.2e1 + t574 * t86 + t575 * t85 + t603 * t119 + (t561 - t587) * t492 + (Ifges(7,1) * t132 + Ifges(7,4) * t131) * t542 + (Ifges(7,4) * t132 + Ifges(7,2) * t131) * t544 + (t170 * t412 + t171 * t415 - t409 * t565) * qJD(5) / 0.2e1 + (-g(1) * t196 + g(2) * t373 - g(3) * t218 - t569 * t43 + (-t450 + t491) * t42 + t567) * mrSges(6,3) + (-m(5) * t387 - mrSges(4,1) * t266 + mrSges(4,2) * t267 - t545 * (t217 * pkin(4) + t387) - t563 * t218 + t564 * t217) * g(3) + (t426 * t56 + t427 * t57) * t364 + (t427 - t491 / 0.2e1) * t91 + (-(-t484 + (-t312 * t358 - t438) * t366) * mrSges(4,1) - t381 * mrSges(4,2) - m(5) * t371 - t545 * (t190 * pkin(4) + t371) + t563 * t373 + t564 * t190) * g(2) + (-t42 * mrSges(6,1) + t43 * mrSges(6,2) - Ifges(5,6) * t531 + Ifges(6,3) * t534 + Ifges(6,5) * t536 + Ifges(6,6) * t537 + t580 / 0.2e1 + t578 / 0.2e1 - t167 * mrSges(5,1) - t89 / 0.2e1 + t121 / 0.2e1 + t77 * mrSges(5,3)) * t398 + (Ifges(5,5) * t531 + t409 * t534 + t415 * t536 + t412 * t537 - t388 - t581 / 0.2e1 - t579 / 0.2e1 - t167 * mrSges(5,2) - t122 / 0.2e1 + t76 * mrSges(5,3)) * t423 + (Ifges(7,5) * t132 + Ifges(7,6) * t131) * t539 - t136 * t527 + t2 * (-mrSges(7,1) * t365 - mrSges(7,3) * t462) - t10 * t464 / 0.2e1 + t1 * (mrSges(7,2) * t365 - mrSges(7,3) * t464) - t256 * (Ifges(4,1) * t253 - t504) / 0.2e1 + t587 * t451 + t589 * t553 + t599 - t458 * t79 + (mrSges(7,1) * t569 + mrSges(7,3) * t571) * t13 + (mrSges(7,1) * t572 - mrSges(7,2) * t571) * t36 + (-mrSges(7,2) * t569 - mrSges(7,3) * t572) * t14 - t365 * t9 / 0.2e1 + t365 * t34 / 0.2e1 + t292 + t293 + (t1 * t271 + t574 * t13 + t575 * t14 + t2 * t270 - t36 * t46) * m(7) + (t29 * t349 - t42 * t52 - t43 * t53 - t74 * t79) * m(6) + (t511 * t361 + (t36 * t450 + t518) * m(7) + ((-t361 * t43 - t365 * t42) * qJD(5) + t567) * m(6) + t365 * t62) * t348 + t349 * t48 + (t507 + t223) * t166 - t29 * t418 + t132 * t551 + (-Ifges(7,6) * t365 + t361 * t411) * t555 + (-Ifges(7,5) * t365 + t361 * t414) * t556 + t361 * t557 + t462 * t558 + (-t410 * t448 + (Ifges(7,6) * t361 + t365 * t411) * qJD(5)) * t543 + (Ifges(6,2) * t365 + t503) * t546 + (Ifges(6,1) * t361 + t502) * t547 + (-Ifges(7,3) * t365 + t361 * t408) * t548 + (Ifges(4,5) * t253 - Ifges(4,6) * t256) * t531 + t178 * t532 + (-t407 * t448 + (Ifges(7,3) * t361 + t365 * t408) * qJD(5)) * t538 + (Ifges(6,5) * t361 + Ifges(6,6) * t365) * t540 + (-t413 * t448 + (Ifges(7,5) * t361 + t365 * t414) * qJD(5)) * t541 + t125 * t524 + t124 * t526 + t416 * t518 + t270 * t24 + t271 * t25 - t212 * (mrSges(4,1) * t256 + mrSges(4,2) * t253) + (t508 - t222) * t165 - t80 * t174 + t360 * t57 * t426 + qJD(5) * t388; -t131 * t86 - t132 * t85 - t423 * t174 - t458 * t398 + (-t423 * t119 + (-t360 * t86 + t364 * t85 + t119) * qJD(5) - t511) * t365 + (t62 + (-t360 * t85 - t364 * t86) * qJD(6) + t565 * t576 + t568) * t361 + t84 + (t382 + (-t6 + (-t13 * t360 + t14 * t364) * qJD(5)) * t365 + (qJD(5) * t36 + t369) * t361 - t13 * t131 - t132 * t14 - t36 * t492) * m(7) + (t361 * t7 + t365 * t8 - t398 * t74 + t382 - t565 * (-t361 * t42 + t365 * t43)) * m(6) + (t398 * t76 - t423 * t77 + t138 + t382) * m(5); (-t498 + t55) * t536 + (t168 + t91) * t537 + (t117 * t411 + t118 * t414 + t169 * t408) * qJD(6) / 0.2e1 + (-pkin(5) * t6 - t13 * t26 - t14 * t27) * m(7) + (-m(7) * t36 + t505 + t576) * t43 + (t163 * t577 + t164 * t422) * g(1) + (t422 * t188 - t577 * (-t218 * t361 + t301)) * g(3) + (t506 - t119) * t42 + ((-t449 + t494) * t14 + (-t447 + t493) * t13 + t566) * mrSges(7,3) + (m(7) * t369 - t447 * t86 - t449 * t85 + t568) * pkin(11) + (-t422 * t583 - t577 * t584) * g(2) + (-Ifges(6,2) * t537 - Ifges(6,6) * t534 + t561 - t582) * t171 + (-t74 * mrSges(6,2) + Ifges(6,1) * t536 + Ifges(6,5) * t534 + t408 * t539 + t411 * t544 + t414 * t542 - t389) * t170 + t586 - pkin(5) * t15 + t6 * t417 + t447 * t550 + t493 * t551 + t494 * t552 + t449 * t553 + t410 * t555 + t413 * t556 + t360 * t558 + t364 * t559 + t407 * t548 + t90 * t535 - t27 * t85 - t26 * t86 + qJD(6) * t389 + t444; -t36 * (mrSges(7,1) * t118 + mrSges(7,2) * t117) - t13 * t85 + t14 * t86 + (Ifges(7,1) * t117 - t501) * t542 + t56 * t541 + (Ifges(7,5) * t117 - Ifges(7,6) * t118) * t539 - g(1) * (mrSges(7,1) * t105 - mrSges(7,2) * t106) - g(2) * ((-t190 * t364 + t360 * t583) * mrSges(7,1) + (t190 * t360 + t364 * t583) * mrSges(7,2)) - g(3) * ((-t188 * t360 - t217 * t364) * mrSges(7,1) + (-t188 * t364 + t217 * t360) * mrSges(7,2)) + (t117 * t13 + t118 * t14) * mrSges(7,3) + t9 + (-Ifges(7,2) * t118 + t116 + t57) * t544 + t585;];
tau  = t11;
