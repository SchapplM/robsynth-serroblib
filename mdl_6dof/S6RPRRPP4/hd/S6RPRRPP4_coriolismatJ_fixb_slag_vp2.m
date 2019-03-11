% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-03-09 04:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRPP4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP4_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP4_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP4_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP4_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:38:53
% EndTime: 2019-03-09 04:39:10
% DurationCPUTime: 9.65s
% Computational Cost: add. (21134->577), mult. (43311->755), div. (0->0), fcn. (48936->8), ass. (0->303)
t340 = sin(qJ(4));
t474 = sin(pkin(10));
t422 = t474 * t340;
t475 = cos(pkin(10));
t519 = cos(qJ(4));
t312 = -t475 * t519 + t422;
t332 = -pkin(4) * t519 - pkin(3);
t424 = t475 * t340;
t356 = t474 * t519 + t424;
t246 = t312 * pkin(5) - qJ(6) * t356 + t332;
t260 = mrSges(7,1) * t312 - mrSges(7,3) * t356;
t604 = m(7) * t246 + t260;
t261 = mrSges(6,1) * t312 + mrSges(6,2) * t356;
t603 = m(6) * t332 + t261;
t533 = -t312 / 0.2e1;
t530 = t356 / 0.2e1;
t600 = Ifges(6,1) + Ifges(7,1);
t599 = Ifges(7,4) + Ifges(6,5);
t598 = Ifges(6,6) - Ifges(7,6);
t338 = sin(pkin(9));
t339 = cos(pkin(9));
t341 = sin(qJ(3));
t520 = cos(qJ(3));
t318 = t338 * t520 + t341 * t339;
t230 = t312 * t318;
t232 = t356 * t318;
t590 = mrSges(6,3) + mrSges(7,2);
t597 = t590 * (t230 * t530 + t232 * t533);
t315 = t338 * t341 - t339 * t520;
t488 = t232 * mrSges(6,3);
t171 = -mrSges(6,2) * t315 - t488;
t306 = t315 * mrSges(7,3);
t489 = t232 * mrSges(7,2);
t176 = t306 - t489;
t583 = t176 + t171;
t596 = t340 ^ 2 + t519 ^ 2;
t528 = -t356 / 0.2e1;
t442 = t318 * t519;
t399 = t442 / 0.2e1;
t333 = t519 * qJ(5);
t335 = t519 * pkin(8);
t457 = t335 + t333;
t511 = -qJ(5) - pkin(8);
t565 = t511 * t422 + t457 * t475;
t595 = t565 / 0.2e1;
t507 = Ifges(5,4) * t340;
t373 = -Ifges(5,2) * t519 - t507;
t594 = -t373 * t340 / 0.2e1;
t492 = t230 * mrSges(6,3);
t174 = mrSges(6,1) * t315 + t492;
t493 = t230 * mrSges(7,2);
t175 = -mrSges(7,1) * t315 - t493;
t402 = t174 * t528 + t175 * t530 + t533 * t583;
t444 = -pkin(2) * t339 - pkin(1);
t515 = pkin(8) * t318;
t247 = pkin(3) * t315 + t444 - t515;
t512 = pkin(7) + qJ(2);
t322 = t512 * t339;
t425 = t512 * t338;
t270 = t322 * t520 - t341 * t425;
t439 = t519 * t270;
t120 = t439 + (-qJ(5) * t318 + t247) * t340;
t113 = t475 * t120;
t153 = t519 * t247 - t340 * t270;
t119 = -t318 * t333 + t153;
t89 = t315 * pkin(4) + t119;
t54 = t474 * t89 + t113;
t42 = qJ(6) * t315 + t54;
t423 = t474 * t120;
t53 = t475 * t89 - t423;
t43 = -t315 * pkin(5) - t53;
t462 = t318 * t340;
t450 = mrSges(5,3) * t462;
t254 = -mrSges(5,2) * t315 - t450;
t440 = t519 * t254;
t403 = mrSges(5,3) * t442;
t256 = t315 * mrSges(5,1) - t403;
t461 = t340 * t256;
t553 = m(7) / 0.2e1;
t555 = m(6) / 0.2e1;
t62 = t119 * t474 + t113;
t63 = t119 * t475 - t423;
t593 = -(-(t53 - t63) * t356 + (-t54 + t62) * t312) * t555 - (-(-t43 - t63) * t356 + (-t42 + t62) * t312) * t553 + t461 / 0.2e1 - t440 / 0.2e1 - t597 - t402;
t231 = t356 * t315;
t233 = t315 * t312;
t587 = t599 * t318 + t600 * t233 + (Ifges(6,4) - Ifges(7,5)) * t231;
t225 = Ifges(6,4) * t232;
t502 = Ifges(7,5) * t232;
t586 = -t230 * t600 + t315 * t599 - t225 + t502;
t170 = -mrSges(6,2) * t318 + t231 * mrSges(6,3);
t177 = t231 * mrSges(7,2) + mrSges(7,3) * t318;
t585 = t170 + t177;
t584 = t174 - t175;
t309 = Ifges(6,4) * t312;
t501 = Ifges(7,5) * t312;
t582 = t356 * t600 - t309 + t501;
t334 = Ifges(5,4) * t519;
t324 = Ifges(5,1) * t340 + t334;
t580 = -Ifges(5,2) * t340 + t334;
t581 = t324 + t580;
t575 = -t596 * mrSges(5,3) / 0.2e1;
t574 = t230 * t598 - t232 * t599;
t573 = -t519 * mrSges(5,1) + t340 * mrSges(5,2);
t249 = m(7) * t356;
t571 = qJD(6) * t249;
t374 = Ifges(5,5) * t519 - Ifges(5,6) * t340;
t570 = -t312 * t599 - t356 * t598 + t374;
t308 = Ifges(7,5) * t356;
t262 = t312 * Ifges(7,3) + t308;
t505 = Ifges(6,4) * t356;
t569 = -t312 * t600 + t262 + t308 - t505;
t568 = -Ifges(6,2) * t356 - t309 + t582;
t567 = Ifges(6,2) * t230 - t225 + t586;
t222 = Ifges(7,5) * t230;
t100 = t315 * Ifges(7,6) + t232 * Ifges(7,3) - t222;
t506 = Ifges(6,4) * t230;
t566 = -t232 * t600 + t100 - t222 + t506;
t437 = t474 * pkin(4);
t325 = t437 + qJ(6);
t321 = m(7) * t325 + mrSges(7,3);
t277 = -t511 * t424 + t457 * t474;
t563 = m(6) * t54 + m(7) * t42 + t583;
t562 = -m(6) * t53 + m(7) * t43 - t584;
t438 = t475 * pkin(4);
t331 = -t438 - pkin(5);
t552 = m(6) * pkin(4);
t561 = m(7) * t331 - t475 * t552 - mrSges(6,1) - mrSges(7,1);
t560 = t474 * t552 - mrSges(6,2) + t321;
t559 = t318 ^ 2;
t557 = m(5) / 0.2e1;
t556 = -m(6) / 0.2e1;
t554 = -m(7) / 0.2e1;
t551 = -mrSges(6,2) / 0.2e1;
t549 = mrSges(7,3) / 0.2e1;
t516 = pkin(8) * t315;
t518 = pkin(3) * t318;
t266 = t516 + t518;
t269 = t322 * t341 + t520 * t425;
t164 = t519 * t266 + t340 * t269;
t107 = t318 * pkin(4) + t315 * t333 + t164;
t165 = t340 * t266 - t519 * t269;
t463 = t315 * t340;
t130 = qJ(5) * t463 + t165;
t60 = t107 * t475 - t130 * t474;
t48 = -t318 * pkin(5) - t60;
t548 = -t48 / 0.2e1;
t479 = t318 * mrSges(7,1);
t486 = t233 * mrSges(7,2);
t173 = -t479 + t486;
t547 = t173 / 0.2e1;
t546 = -t230 / 0.2e1;
t544 = t231 / 0.2e1;
t543 = -t231 / 0.2e1;
t542 = -t232 / 0.2e1;
t540 = t232 / 0.2e1;
t539 = t233 / 0.2e1;
t517 = pkin(4) * t340;
t250 = pkin(5) * t356 + qJ(6) * t312 + t517;
t538 = t250 / 0.2e1;
t537 = t260 / 0.2e1;
t531 = t312 / 0.2e1;
t527 = -t315 / 0.2e1;
t526 = t315 / 0.2e1;
t524 = -t318 / 0.2e1;
t523 = t318 / 0.2e1;
t521 = t340 / 0.2e1;
t514 = pkin(8) * t340;
t508 = Ifges(4,4) * t318;
t504 = Ifges(7,4) * t233;
t503 = Ifges(6,5) * t233;
t499 = Ifges(7,2) * t318;
t498 = Ifges(6,6) * t231;
t497 = Ifges(7,6) * t231;
t496 = Ifges(5,3) * t318;
t495 = Ifges(6,3) * t318;
t172 = mrSges(6,1) * t318 - t233 * mrSges(6,3);
t383 = -t231 * t277 + t233 * t565;
t344 = (-t516 * t596 - t518) * t557 + (t318 * t332 + t383) * t555 + (t246 * t318 + t383) * t553 - t573 * t524 + (t260 + t261) * t523 + t575 * t315 + t590 * (-t231 * t530 + t233 * t533);
t253 = -mrSges(5,2) * t318 + mrSges(5,3) * t463;
t443 = t315 * t519;
t255 = t318 * mrSges(5,1) + mrSges(5,3) * t443;
t452 = t519 / 0.2e1;
t61 = t474 * t107 + t475 * t130;
t47 = qJ(6) * t318 + t61;
t348 = (t164 * t519 + t340 * t165) * t557 + (-t312 * t60 + t356 * t61) * t555 + (t312 * t48 + t356 * t47) * t553 + t253 * t521 + t255 * t452;
t419 = t318 * mrSges(4,1) - t315 * mrSges(4,2);
t8 = t172 * t533 + t173 * t531 + t530 * t585 - t344 + t348 + t419;
t494 = qJD(1) * t8;
t491 = t231 * mrSges(6,1);
t490 = t231 * mrSges(7,1);
t487 = t233 * mrSges(6,2);
t485 = t233 * mrSges(7,3);
t101 = Ifges(6,4) * t233 + Ifges(6,2) * t231 + t318 * Ifges(6,6);
t102 = -t232 * Ifges(6,2) + t315 * Ifges(6,6) - t506;
t126 = -t485 - t490;
t127 = t487 - t491;
t128 = mrSges(7,1) * t232 + mrSges(7,3) * t230;
t129 = mrSges(6,1) * t232 - mrSges(6,2) * t230;
t154 = t340 * t247 + t439;
t191 = Ifges(5,6) * t318 - t315 * t580;
t192 = Ifges(5,6) * t315 + t318 * t580;
t377 = Ifges(5,1) * t519 - t507;
t193 = Ifges(5,5) * t318 - t315 * t377;
t364 = t377 * t318;
t194 = Ifges(5,5) * t315 + t364;
t202 = pkin(4) * t462 + t269;
t203 = -pkin(4) * t463 + t270;
t370 = t340 * mrSges(5,1) + mrSges(5,2) * t519;
t244 = t370 * t315;
t245 = t370 * t318;
t400 = -t443 / 0.2e1;
t427 = -t462 / 0.2e1;
t428 = t463 / 0.2e1;
t74 = t232 * pkin(5) + t230 * qJ(6) + t202;
t75 = -t231 * pkin(5) - t233 * qJ(6) + t203;
t99 = Ifges(7,5) * t233 + t318 * Ifges(7,6) - Ifges(7,3) * t231;
t3 = t102 * t544 + t586 * t539 + t587 * t546 + t191 * t427 + t192 * t428 + t193 * t399 + t194 * t400 + (t374 * t318 - t508 - t598 * t232 - t599 * t230 + (-Ifges(4,1) + Ifges(7,2) + Ifges(5,3) + Ifges(6,3)) * t315) * t523 + t99 * t540 + t101 * t542 + t100 * t543 + (-Ifges(4,2) * t315 + t508) * t524 - t269 * t244 + t270 * t245 + t154 * t253 + t165 * t254 + t153 * t255 + t164 * t256 + t444 * t419 + t202 * t127 + t203 * t129 + t54 * t170 + t61 * t171 + t53 * t172 + t43 * t173 + t60 * t174 + t48 * t175 + t47 * t176 + t42 * t177 + t74 * t126 + t75 * t128 + (-t315 * t374 + t495 + t496 - t497 + t498 + t499 + t503 + t504) * t526 + (-0.2e1 * Ifges(4,4) * t315 + (Ifges(4,1) - Ifges(4,2)) * t318) * t527 + m(5) * (t153 * t164 + t154 * t165 + t269 * t270) + m(6) * (t202 * t203 + t53 * t60 + t54 * t61) + m(7) * (t42 * t47 + t43 * t48 + t74 * t75);
t484 = t3 * qJD(1);
t483 = t312 * mrSges(7,2);
t482 = t312 * mrSges(6,3);
t481 = t356 * mrSges(7,2);
t480 = t356 * mrSges(6,3);
t365 = t318 * t573;
t372 = -Ifges(5,5) * t340 - Ifges(5,6) * t519;
t390 = -Ifges(7,3) * t230 - t502;
t406 = pkin(4) * t442;
t418 = -t230 * mrSges(6,1) - t232 * mrSges(6,2);
t421 = -t230 * mrSges(7,1) + t232 * mrSges(7,3);
t91 = -t230 * pkin(5) + t232 * qJ(6) + t406;
t4 = t315 * t372 * t523 + t194 * t427 + t390 * t540 + t42 * t493 + t53 * t488 + t54 * t492 - t43 * t489 + t230 * t102 / 0.2e1 - t192 * t442 / 0.2e1 - t269 * t365 + t91 * t128 + t129 * t406 + (m(7) * t91 + t421) * t74 + (m(6) * t406 + t418) * t202 + t563 * t63 + t562 * t62 + (-t324 * t452 + t594) * t559 + t574 * t526 + (-t256 - t403) * t154 + (t450 + t254) * t153 + t566 * t546 + t567 * t542;
t477 = t4 * qJD(1);
t7 = t559 * mrSges(4,3) + t563 * t233 - t562 * t231 + (m(5) * (t153 * t340 - t154 * t519) - m(4) * t270 - t440 + t461 + mrSges(4,3) * t315) * t315 + (m(3) * qJ(2) + mrSges(3,3)) * (t338 ^ 2 + t339 ^ 2) + ((m(5) + m(4)) * t269 + m(6) * t202 + m(7) * t74 + t128 + t129 + t245) * t318;
t476 = t7 * qJD(1);
t12 = -t583 * t232 + t584 * t230 + m(6) * (t230 * t53 - t232 * t54) + m(7) * (-t230 * t43 - t232 * t42);
t473 = qJD(1) * t12;
t20 = t315 * t176 + m(7) * (t230 * t74 + t315 * t42) + t230 * t128;
t472 = qJD(1) * t20;
t430 = t356 * t527;
t135 = (-t430 + t544) * m(7);
t455 = t135 * qJD(1);
t454 = t552 / 0.2e1;
t453 = t62 * t553;
t449 = t517 / 0.2e1;
t429 = t312 * t527;
t420 = mrSges(7,1) * t356 + t312 * mrSges(7,3);
t417 = mrSges(6,1) * t356 - t312 * mrSges(6,2);
t411 = -t230 * t277 - t232 * t565;
t405 = mrSges(6,3) * t438;
t404 = mrSges(6,3) * t437;
t397 = 0.2e1 * (m(6) / 0.4e1 + m(7) / 0.4e1) * t318;
t389 = Ifges(7,3) * t356 - t501;
t388 = t277 * t62 + t565 * t63;
t387 = pkin(4) * t399;
t350 = (-t230 * t331 - t232 * t325) * t553 + (t230 * t475 - t232 * t474) * t454;
t357 = m(6) * t387 + t91 * t553;
t25 = -t350 + t357 + t418 + t421;
t354 = (-t312 * t474 - t356 * t475) * t552;
t368 = m(7) * (-t312 * t325 + t331 * t356);
t349 = t368 / 0.2e1 + t354 / 0.2e1;
t366 = t417 + t420;
t367 = m(6) * t449 + m(7) * t538;
t57 = -t349 + t366 + t367;
t386 = qJD(1) * t25 + qJD(3) * t57;
t352 = (t554 + t556) * (-t230 * t312 - t232 * t356);
t35 = t397 + t352;
t385 = qJD(1) * t35;
t108 = t604 * t356;
t351 = (t246 * t230 + t315 * t565 - t356 * t74) * t553 + t230 * t537 + t128 * t528;
t378 = m(7) * t548 + t479 / 0.2e1;
t19 = (t429 - t233 / 0.2e1) * mrSges(7,2) + t351 + t378;
t382 = qJD(1) * t19 - qJD(3) * t108;
t118 = m(7) * t230;
t381 = qJD(1) * t118 - qJD(3) * t249;
t380 = t203 * t556 + t554 * t75;
t263 = -t312 * Ifges(6,2) + t505;
t342 = t570 * t315 / 0.4e1 + t575 * t515 + t129 * t449 + (t566 / 0.4e1 - t102 / 0.4e1) * t356 - t565 * t174 / 0.2e1 + t91 * t537 + t128 * t538 + t373 * t399 - t254 * t514 / 0.2e1 + t53 * t482 / 0.2e1 - t43 * t483 / 0.2e1 - t54 * t480 / 0.2e1 - t42 * t481 / 0.2e1 - t340 * t192 / 0.4e1 - t256 * t335 / 0.2e1 + t202 * t417 / 0.2e1 + t332 * t418 / 0.2e1 + t74 * t420 / 0.2e1 + t246 * t421 / 0.2e1 + t261 * t387 + (-t567 / 0.4e1 + t390 / 0.4e1) * t312 + (-t568 / 0.4e1 + t389 / 0.4e1) * t232 + (-t569 / 0.4e1 + t263 / 0.4e1) * t230 + (t364 + t194) * t519 / 0.4e1 + pkin(3) * t365 / 0.2e1 - t583 * t277 / 0.2e1 + (t246 * t91 + t250 * t74 - t277 * t42 + t43 * t565 + t388) * t553 + (-t565 * t53 - t277 * t54 + (t202 * t340 + t332 * t442) * pkin(4) + t388) * t555 + (-t581 / 0.4e1 - t324 / 0.4e1) * t462 + t590 * (t230 * t595 + t277 * t542 + t62 * t530 + t63 * t533) + t175 * t595 + t269 * t370 / 0.2e1;
t343 = (t325 * t47 + t331 * t48) * t553 + t504 / 0.2e1 + t503 / 0.2e1 + t499 / 0.2e1 + t498 / 0.2e1 - t497 / 0.2e1 + t496 / 0.2e1 + t495 / 0.2e1 + t164 * mrSges(5,1) / 0.2e1 - t165 * mrSges(5,2) / 0.2e1 + t325 * t177 / 0.2e1 + t331 * t547 + t47 * t549 + mrSges(7,1) * t548 + t60 * mrSges(6,1) / 0.2e1 + t61 * t551 + (t474 * t61 + t475 * t60) * t454 + Ifges(5,5) * t400 + Ifges(5,6) * t428 + t172 * t438 / 0.2e1 + t170 * t437 / 0.2e1;
t1 = t343 - t342;
t13 = -pkin(3) * t370 + t246 * t420 + t250 * t604 + t263 * t528 + t332 * t417 + t377 * t521 + t389 * t531 + t581 * t452 + t603 * t517 + t569 * t530 + t568 * t533 - t594;
t363 = -t1 * qJD(1) + t13 * qJD(3);
t353 = -t491 / 0.2e1 - t490 / 0.2e1 + t487 / 0.2e1 - t485 / 0.2e1;
t347 = (-t231 * t331 + t233 * t325) * t553 + (t231 * t475 + t233 * t474) * t454 + mrSges(5,1) * t428 + mrSges(5,2) * t443 / 0.2e1 - t353;
t5 = t347 + t593;
t362 = t5 * qJD(1);
t346 = (-t312 * t54 - t356 * t53 + t411) * t555 + (-t312 * t42 + t356 * t43 + t411) * t553 + t402 - t597;
t11 = (t549 + t551) * t233 - (-mrSges(7,1) / 0.2e1 - mrSges(6,1) / 0.2e1) * t231 + t346 + t380;
t34 = (t312 ^ 2 + t356 ^ 2) * t590 + (m(6) + m(7)) * (t277 * t356 - t565 * t312);
t360 = qJD(1) * t11 + qJD(3) * t34;
t355 = -t306 + ((qJ(6) + t325) * t315 + t54) * t554;
t23 = t453 + t355;
t359 = qJD(1) * t23 - qJD(4) * t321;
t134 = (-t430 + t543) * m(7);
t109 = 0.2e1 * m(7) * t595 - t483;
t73 = t349 + t367;
t36 = t397 - t352;
t32 = t350 + t357;
t24 = -t355 + t453 - t489;
t18 = mrSges(7,2) * t429 + t486 / 0.2e1 + t351 - t378;
t10 = t346 + t353 - t380;
t9 = -(-t170 / 0.2e1 - t177 / 0.2e1) * t356 + (-t172 / 0.2e1 + t547) * t312 + t344 + t348;
t6 = t347 - t593;
t2 = t343 + t342;
t14 = [qJD(2) * t7 + qJD(3) * t3 + qJD(4) * t4 + qJD(5) * t12 + qJD(6) * t20, t9 * qJD(3) + t6 * qJD(4) + t36 * qJD(5) + t134 * qJD(6) + t476 + 0.2e1 * (t555 + t553) * qJD(2) * (-t231 * t312 + t233 * t356) t484 + t9 * qJD(2) + (t263 * t544 + t582 * t539 + t587 * t530 + t324 * t400 - t373 * t428 + (m(6) * t61 + m(7) * t47 + t585) * t565 + (m(5) * pkin(8) + mrSges(5,3)) * (-t164 * t340 + t519 * t165) + (-t312 * t598 + t356 * t599) * t523 + t101 * t533 + t262 * t543 + t99 * t531 + t372 * t524 + t193 * t521 - t60 * t480 + t191 * t452 + t253 * t335 - t255 * t514 - t47 * t483 + t48 * t481 - t61 * t482 + t332 * t127 - Ifges(4,6) * t318 - Ifges(4,5) * t315 + t269 * mrSges(4,2) + t75 * t260 + pkin(3) * t244 + t603 * t203 + (-m(6) * t60 + m(7) * t48 - t172 + t173) * t277 + (m(7) * t75 + t126) * t246 + (-m(5) * pkin(3) - mrSges(4,1) + t573) * t270) * qJD(3) + t2 * qJD(4) + t10 * qJD(5) + t18 * qJD(6), t477 + t6 * qJD(2) + t2 * qJD(3) + (-t154 * mrSges(5,1) - t153 * mrSges(5,2) - Ifges(5,5) * t462 - Ifges(5,6) * t442 + t230 * t404 + t232 * t405 + t325 * t493 - t331 * t489 + t560 * t63 + t561 * t62 + t574) * qJD(4) + t32 * qJD(5) + t24 * qJD(6), qJD(2) * t36 + qJD(3) * t10 + qJD(4) * t32 + t473, qJD(2) * t134 + qJD(3) * t18 + qJD(4) * t24 + t472; qJD(3) * t8 - qJD(4) * t5 - qJD(5) * t35 + qJD(6) * t135 - t476, 0, t494 (t354 + t368 - t366 - t370) * qJD(4) + t571 - t362, -t385, qJD(4) * t249 + t455; -qJD(2) * t8 - qJD(4) * t1 + qJD(5) * t11 + qJD(6) * t19 - t484, -t494, qJD(4) * t13 + qJD(5) * t34 - qJD(6) * t108 (pkin(8) * t573 - t277 * t560 + t312 * t405 - t325 * t481 - t331 * t483 - t356 * t404 + t561 * t565 + t570) * qJD(4) + t73 * qJD(5) + t109 * qJD(6) + t363, qJD(4) * t73 + t360, qJD(4) * t109 + t382; qJD(2) * t5 + qJD(3) * t1 - qJD(5) * t25 - qJD(6) * t23 - t477, t362, -qJD(5) * t57 - t363, t321 * qJD(6), -t386, -t359; qJD(2) * t35 - qJD(3) * t11 + qJD(4) * t25 + qJD(6) * t118 - t473, t385, qJD(4) * t57 - t360 - t571, t386, 0, t381; -qJD(2) * t135 - qJD(3) * t19 + qJD(4) * t23 - qJD(5) * t118 - t472, -t455, qJD(5) * t249 - t382, t359, -t381, 0;];
Cq  = t14;
