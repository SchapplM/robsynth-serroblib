% Calculate vector of inverse dynamics joint torques for
% S6PRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRPR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_invdynJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:23:04
% EndTime: 2019-03-08 23:24:10
% DurationCPUTime: 39.20s
% Computational Cost: add. (17933->952), mult. (44919->1347), div. (0->0), fcn. (37905->18), ass. (0->421)
t354 = sin(pkin(6));
t362 = sin(qJ(2));
t488 = t354 * t362;
t455 = qJD(1) * t488;
t353 = sin(pkin(7));
t476 = qJD(2) * t353;
t312 = pkin(9) * t476 + t455;
t361 = sin(qJ(3));
t365 = cos(qJ(3));
t535 = cos(qJ(2));
t459 = t354 * t535;
t426 = qJD(1) * t459;
t323 = qJD(2) * pkin(2) + t426;
t356 = cos(pkin(7));
t357 = cos(pkin(6));
t477 = qJD(1) * t357;
t456 = t353 * t477;
t634 = t323 * t356 + t456;
t191 = -t361 * t312 + t365 * t634;
t395 = (pkin(3) * t361 - pkin(10) * t365) * t353;
t289 = qJD(2) * t395;
t360 = sin(qJ(4));
t364 = cos(qJ(4));
t140 = -t191 * t360 + t364 * t289;
t358 = -qJ(5) - pkin(10);
t440 = qJD(4) * t358;
t482 = t364 * t365;
t636 = -(pkin(4) * t361 - qJ(5) * t482) * t476 - t140 - qJD(5) * t360 + t364 * t440;
t141 = t364 * t191 + t360 * t289;
t453 = t365 * t476;
t427 = t360 * t453;
t635 = -qJ(5) * t427 - qJD(5) * t364 - t360 * t440 + t141;
t457 = t535 * t365;
t484 = t361 * t362;
t381 = -t356 * t484 + t457;
t265 = t381 * t354;
t250 = qJD(1) * t265;
t429 = t353 * t455;
t210 = -t250 * t360 + t364 * t429;
t212 = t250 * t364 + t360 * t429;
t352 = sin(pkin(13));
t355 = cos(pkin(13));
t487 = t356 * t361;
t489 = t353 * t365;
t311 = pkin(2) * t487 + pkin(9) * t489;
t277 = pkin(10) * t356 + t311;
t410 = -pkin(3) * t365 - pkin(10) * t361;
t278 = (-pkin(2) + t410) * t353;
t194 = t364 * t277 + t360 * t278;
t290 = qJD(3) * t395;
t491 = t353 * t361;
t341 = pkin(9) * t491;
t486 = t356 * t365;
t310 = pkin(2) * t486 - t341;
t291 = t310 * qJD(3);
t124 = -qJD(4) * t194 + t364 * t290 - t291 * t360;
t308 = t356 * t364 - t360 * t491;
t474 = qJD(3) * t365;
t449 = t353 * t474;
t234 = qJD(4) * t308 + t364 * t449;
t490 = t353 * t364;
t309 = t356 * t360 + t361 * t490;
t475 = qJD(3) * t361;
t450 = t353 * t475;
t68 = pkin(4) * t450 - qJ(5) * t234 - qJD(5) * t309 + t124;
t472 = qJD(4) * t364;
t473 = qJD(4) * t360;
t123 = -t277 * t473 + t278 * t472 + t360 * t290 + t364 * t291;
t235 = -qJD(4) * t309 - t360 * t449;
t75 = qJ(5) * t235 + qJD(5) * t308 + t123;
t608 = (-t212 + t75) * t355 + (-t210 + t68) * t352;
t604 = t352 * t636 - t635 * t355;
t192 = t365 * t312 + t323 * t487 + t361 * t456;
t600 = -t192 + (-t427 + t473) * pkin(4);
t328 = qJD(4) - t453;
t538 = -t328 / 0.2e1;
t340 = qJD(2) * t356 + qJD(3);
t454 = t361 * t476;
t262 = t340 * t364 - t360 * t454;
t263 = t340 * t360 + t364 * t454;
t396 = t262 * t352 + t355 * t263;
t543 = -t396 / 0.2e1;
t435 = t355 * t262 - t263 * t352;
t545 = -t435 / 0.2e1;
t633 = -Ifges(6,4) * t543 - Ifges(6,2) * t545 - Ifges(6,6) * t538;
t181 = qJD(6) - t435;
t546 = t181 / 0.2e1;
t359 = sin(qJ(6));
t363 = cos(qJ(6));
t147 = t328 * t359 + t363 * t396;
t552 = t147 / 0.2e1;
t146 = t328 * t363 - t359 * t396;
t554 = t146 / 0.2e1;
t632 = Ifges(7,5) * t552 + Ifges(7,6) * t554 + Ifges(7,3) * t546;
t614 = Ifges(5,3) + Ifges(6,3);
t407 = -mrSges(7,1) * t363 + mrSges(7,2) * t359;
t631 = m(7) * pkin(5) + mrSges(6,1) - t407;
t587 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t630 = -pkin(11) * t450 - t608;
t629 = -pkin(11) * t454 + t604;
t149 = t234 * t352 - t355 * t235;
t150 = t234 * t355 + t235 * t352;
t292 = t311 * qJD(3);
t189 = -pkin(4) * t235 + t292;
t458 = t535 * t361;
t483 = t362 * t365;
t380 = t356 * t483 + t458;
t264 = t380 * t354;
t249 = qJD(1) * t264;
t628 = pkin(5) * t149 - pkin(11) * t150 + t189 - t249;
t315 = t352 * t364 + t355 * t360;
t237 = t315 * t453;
t314 = t352 * t360 - t355 * t364;
t238 = t314 * t453;
t306 = t315 * qJD(4);
t307 = t314 * qJD(4);
t627 = t600 + (-t238 + t307) * pkin(11) + (-t237 + t306) * pkin(5);
t171 = -pkin(3) * t340 - t191;
t129 = -pkin(4) * t262 + qJD(5) + t171;
t172 = pkin(10) * t340 + t192;
t336 = t356 * t477;
t213 = t336 + (qJD(2) * t410 - t323) * t353;
t115 = -t172 * t360 + t364 * t213;
t89 = -qJ(5) * t263 + t115;
t77 = pkin(4) * t328 + t89;
t116 = t172 * t364 + t213 * t360;
t90 = qJ(5) * t262 + t116;
t85 = t355 * t90;
t39 = t352 * t77 + t85;
t35 = pkin(11) * t328 + t39;
t56 = -pkin(5) * t435 - pkin(11) * t396 + t129;
t14 = -t35 * t359 + t363 * t56;
t15 = t35 * t363 + t359 * t56;
t626 = mrSges(6,1) * t129 + mrSges(7,1) * t14 - mrSges(7,2) * t15 - mrSges(6,3) * t39 + t632 - t633;
t605 = t635 * t352 + t355 * t636;
t379 = t356 * t458 + t483;
t229 = t354 * t379 + t357 * t491;
t303 = -t353 * t459 + t357 * t356;
t177 = -t229 * t360 + t303 * t364;
t511 = cos(pkin(12));
t412 = t511 * t535;
t510 = sin(pkin(12));
t436 = t510 * t362;
t305 = -t357 * t436 + t412;
t411 = t510 * t535;
t437 = t511 * t362;
t372 = t357 * t411 + t437;
t438 = t354 * t510;
t414 = t353 * t438;
t176 = t305 * t365 + (-t356 * t372 + t414) * t361;
t231 = t353 * t372 + t356 * t438;
t625 = -t176 * t360 + t231 * t364;
t304 = t357 * t437 + t411;
t371 = -t357 * t412 + t436;
t369 = t371 * t361;
t439 = t354 * t511;
t415 = t353 * t439;
t174 = t304 * t365 - t356 * t369 - t361 * t415;
t230 = t353 * t371 - t356 * t439;
t624 = -t174 * t360 + t230 * t364;
t537 = t328 / 0.2e1;
t542 = t396 / 0.2e1;
t544 = t435 / 0.2e1;
t623 = -Ifges(6,4) * t542 - Ifges(6,2) * t544 - Ifges(6,6) * t537 + t626 + t632;
t547 = -t181 / 0.2e1;
t553 = -t147 / 0.2e1;
t555 = -t146 / 0.2e1;
t622 = Ifges(7,5) * t553 + Ifges(7,6) * t555 + Ifges(7,3) * t547 - t626 + t633;
t333 = qJD(2) * t426;
t298 = qJDD(1) * t488 + t333;
t621 = pkin(9) * qJDD(2) * t353 + qJD(3) * t634 + t298;
t406 = mrSges(7,1) * t359 + mrSges(7,2) * t363;
t620 = -m(5) * pkin(10) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3) - t406;
t469 = qJD(2) * qJD(3);
t295 = (-qJDD(2) * t365 + t361 * t469) * t353;
t282 = qJDD(4) + t295;
t539 = t282 / 0.2e1;
t296 = (qJDD(2) * t361 + t365 * t469) * t353;
t339 = qJDD(2) * t356 + qJDD(3);
t169 = qJD(4) * t262 + t296 * t364 + t339 * t360;
t170 = -qJD(4) * t263 - t296 * t360 + t339 * t364;
t103 = t169 * t355 + t170 * t352;
t560 = t103 / 0.2e1;
t102 = -t169 * t352 + t170 * t355;
t561 = t102 / 0.2e1;
t572 = Ifges(6,1) * t560 + Ifges(6,4) * t561 + Ifges(6,5) * t539;
t49 = qJD(6) * t146 + t103 * t363 + t282 * t359;
t571 = t49 / 0.2e1;
t50 = -qJD(6) * t147 - t103 * t359 + t282 * t363;
t570 = t50 / 0.2e1;
t618 = -m(7) - m(6);
t101 = qJDD(6) - t102;
t562 = t101 / 0.2e1;
t549 = t169 / 0.2e1;
t548 = t170 / 0.2e1;
t214 = -t355 * t308 + t309 * t352;
t215 = t308 * t352 + t309 * t355;
t276 = t341 + (-pkin(2) * t365 - pkin(3)) * t356;
t218 = -pkin(4) * t308 + t276;
t104 = pkin(5) * t214 - pkin(11) * t215 + t218;
t193 = -t277 * t360 + t364 * t278;
t142 = -pkin(4) * t489 - qJ(5) * t309 + t193;
t153 = qJ(5) * t308 + t194;
t71 = t352 * t142 + t355 * t153;
t67 = -pkin(11) * t489 + t71;
t37 = t104 * t359 + t363 * t67;
t617 = -qJD(6) * t37 + t359 * t630 + t363 * t628;
t36 = t104 * t363 - t359 * t67;
t616 = qJD(6) * t36 + t359 * t628 - t363 * t630;
t16 = -mrSges(7,1) * t50 + mrSges(7,2) * t49;
t73 = mrSges(6,1) * t282 - mrSges(6,3) * t103;
t613 = t16 - t73;
t348 = pkin(4) * t364 + pkin(3);
t219 = pkin(5) * t314 - pkin(11) * t315 - t348;
t329 = t358 * t364;
t445 = t358 * t360;
t245 = -t355 * t329 + t352 * t445;
t144 = t219 * t363 - t245 * t359;
t612 = qJD(6) * t144 + t359 * t627 + t363 * t629;
t145 = t219 * t359 + t245 * t363;
t611 = -qJD(6) * t145 - t359 * t629 + t363 * t627;
t452 = qJD(2) * t488;
t425 = qJD(1) * t452;
t297 = qJDD(1) * t459 - t425;
t273 = qJDD(2) * pkin(2) + t297;
t468 = qJDD(1) * t357;
t444 = t353 * t468;
t111 = t273 * t487 - t312 * t475 + t361 * t444 + t365 * t621;
t610 = t111 * mrSges(4,2);
t112 = t365 * (t273 * t356 + t444) - t312 * t474 - t621 * t361;
t609 = t112 * mrSges(4,1);
t152 = mrSges(6,1) * t328 - mrSges(6,3) * t396;
t76 = -mrSges(7,1) * t146 + mrSges(7,2) * t147;
t607 = t152 - t76;
t606 = pkin(5) * t454 - t605;
t603 = t123 - t212;
t602 = t124 - t210;
t601 = t263 * Ifges(5,5) + Ifges(6,5) * t396 + t262 * Ifges(5,6) + Ifges(6,6) * t435 + t328 * t614;
t434 = mrSges(4,3) * t454;
t599 = -mrSges(4,1) * t340 - mrSges(5,1) * t262 + mrSges(5,2) * t263 + t434;
t209 = t238 * t359 + t363 * t454;
t470 = qJD(6) * t363;
t391 = -t359 * t307 + t315 * t470;
t598 = t209 + t391;
t211 = -t238 * t363 + t359 * t454;
t471 = qJD(6) * t359;
t390 = t363 * t307 + t315 * t471;
t597 = t211 + t390;
t178 = t229 * t364 + t303 * t360;
t109 = t177 * t352 + t178 * t355;
t595 = t356 * t457 - t484;
t228 = -t354 * t595 - t357 * t489;
t225 = t228 * t363;
t63 = -t109 * t359 + t225;
t596 = -t250 + t291;
t594 = Ifges(5,5) * t169 + Ifges(6,5) * t103 + Ifges(5,6) * t170 + Ifges(6,6) * t102 + t282 * t614;
t226 = -t273 * t353 + t356 * t468;
t139 = pkin(3) * t295 - pkin(10) * t296 + t226;
t91 = pkin(10) * t339 + t111;
t28 = t360 * t139 - t172 * t473 + t213 * t472 + t364 * t91;
t29 = -qJD(4) * t116 + t364 * t139 - t360 * t91;
t593 = t28 * t364 - t29 * t360;
t92 = -pkin(3) * t339 - t112;
t57 = -pkin(4) * t170 + qJDD(5) + t92;
t17 = -pkin(5) * t102 - pkin(11) * t103 + t57;
t19 = pkin(4) * t282 - qJ(5) * t169 - qJD(5) * t263 + t29;
t23 = qJ(5) * t170 + qJD(5) * t262 + t28;
t6 = t352 * t19 + t355 * t23;
t4 = pkin(11) * t282 + t6;
t1 = qJD(6) * t14 + t17 * t359 + t363 * t4;
t2 = -qJD(6) * t15 + t17 * t363 - t359 * t4;
t592 = t1 * t363 - t2 * t359;
t253 = -t323 * t353 + t336;
t591 = (t253 * (mrSges(4,1) * t361 + mrSges(4,2) * t365) + t340 * (Ifges(4,5) * t365 - Ifges(4,6) * t361) / 0.2e1) * t353;
t512 = t352 * t90;
t38 = t355 * t77 - t512;
t589 = mrSges(6,2) * t129 - t38 * mrSges(6,3);
t351 = qJ(4) + pkin(13);
t349 = sin(t351);
t350 = cos(t351);
t409 = -mrSges(5,1) * t364 + mrSges(5,2) * t360;
t586 = m(5) * pkin(3) - t349 * t587 + t350 * t631 + mrSges(4,1) - t409;
t585 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t584 = -m(5) * t171 - t599;
t114 = -mrSges(6,1) * t435 + mrSges(6,2) * t396;
t580 = m(4) * t191 - m(6) * t129 - t114 + t584;
t9 = Ifges(7,5) * t49 + Ifges(7,6) * t50 + Ifges(7,3) * t101;
t577 = mrSges(6,1) * t57 - mrSges(6,3) * t6 + Ifges(7,5) * t571 + Ifges(7,6) * t570 + Ifges(7,3) * t562 + t9 / 0.2e1 + t585 + (-t539 - t282 / 0.2e1) * Ifges(6,6) + (-t561 - t102 / 0.2e1) * Ifges(6,2) + (-t560 - t103 / 0.2e1) * Ifges(6,4);
t576 = t353 ^ 2;
t366 = qJD(2) ^ 2;
t574 = Ifges(7,1) * t571 + Ifges(7,4) * t570 + Ifges(7,5) * t562;
t520 = Ifges(7,4) * t147;
t52 = t146 * Ifges(7,2) + t181 * Ifges(7,6) + t520;
t566 = t52 / 0.2e1;
t143 = Ifges(7,4) * t146;
t53 = t147 * Ifges(7,1) + t181 * Ifges(7,5) + t143;
t565 = -t53 / 0.2e1;
t564 = Ifges(5,4) * t549 + Ifges(5,2) * t548 + Ifges(5,6) * t539;
t563 = Ifges(5,1) * t549 + Ifges(5,4) * t548 + Ifges(5,5) * t539;
t107 = Ifges(6,1) * t396 + Ifges(6,4) * t435 + t328 * Ifges(6,5);
t557 = -t107 / 0.2e1;
t556 = t107 / 0.2e1;
t515 = t263 * Ifges(5,4);
t159 = t262 * Ifges(5,2) + t328 * Ifges(5,6) + t515;
t551 = t159 / 0.2e1;
t254 = Ifges(5,4) * t262;
t160 = t263 * Ifges(5,1) + t328 * Ifges(5,5) + t254;
t550 = t160 / 0.2e1;
t540 = t263 / 0.2e1;
t536 = t363 / 0.2e1;
t533 = pkin(2) * t353;
t532 = pkin(4) * t263;
t531 = pkin(4) * t352;
t530 = pkin(4) * t355;
t524 = Ifges(4,4) * t361;
t523 = Ifges(4,4) * t365;
t522 = Ifges(5,4) * t360;
t521 = Ifges(5,4) * t364;
t519 = Ifges(7,4) * t359;
t518 = Ifges(7,4) * t363;
t517 = t115 * mrSges(5,3);
t516 = t116 * mrSges(5,3);
t506 = t435 * t359;
t505 = t435 * t363;
t504 = t228 * t359;
t499 = t304 * t353;
t498 = t305 * t353;
t497 = t315 * t359;
t496 = t315 * t363;
t494 = t349 * t353;
t493 = t350 * t353;
t492 = t353 * t360;
t463 = t353 * t488;
t478 = pkin(2) * t459 + pkin(9) * t463;
t462 = t53 * t536;
t460 = Ifges(4,5) * t296 - Ifges(4,6) * t295 + Ifges(4,3) * t339;
t447 = t491 / 0.2e1;
t441 = -t471 / 0.2e1;
t46 = -t102 * mrSges(6,1) + t103 * mrSges(6,2);
t433 = mrSges(4,3) * t453;
t432 = t360 * t463;
t428 = t353 * t452;
t420 = t624 * pkin(4);
t419 = t625 * pkin(4);
t418 = t177 * pkin(4);
t293 = t371 * pkin(2);
t417 = pkin(9) * t499 - t293;
t294 = t372 * pkin(2);
t416 = pkin(9) * t498 - t294;
t413 = (pkin(4) * t360 + pkin(9)) * t353;
t405 = Ifges(5,1) * t364 - t522;
t404 = Ifges(7,1) * t363 - t519;
t403 = -Ifges(5,2) * t360 + t521;
t402 = -Ifges(7,2) * t359 + t518;
t401 = Ifges(5,5) * t364 - Ifges(5,6) * t360;
t400 = Ifges(7,5) * t363 - Ifges(7,6) * t359;
t5 = t19 * t355 - t23 * t352;
t32 = -t352 * t75 + t355 * t68;
t93 = -mrSges(7,2) * t181 + mrSges(7,3) * t146;
t94 = mrSges(7,1) * t181 - mrSges(7,3) * t147;
t398 = -t359 * t94 + t363 * t93;
t64 = t109 * t363 + t504;
t70 = t142 * t355 - t153 * t352;
t179 = -t215 * t359 - t363 * t489;
t394 = -t215 * t363 + t359 * t489;
t34 = -pkin(5) * t328 - t38;
t393 = t34 * t406;
t388 = t353 * (-mrSges(4,1) * t365 + mrSges(4,2) * t361);
t387 = (t365 * Ifges(4,2) + t524) * t353;
t368 = t371 * t365;
t173 = t304 * t361 + t356 * t368 + t365 * t415;
t370 = t372 * t365;
t175 = t305 * t361 + t356 * t370 - t365 * t414;
t382 = -g(1) * t175 - g(2) * t173 - g(3) * t228;
t376 = t361 * t576 * (Ifges(4,1) * t365 - t524);
t347 = -pkin(5) - t530;
t334 = Ifges(4,4) * t453;
t288 = qJD(2) * t388;
t287 = -mrSges(4,2) * t340 + t433;
t244 = -t329 * t352 - t355 * t445;
t240 = Ifges(4,1) * t454 + t340 * Ifges(4,5) + t334;
t239 = t340 * Ifges(4,6) + qJD(2) * t387;
t233 = mrSges(4,1) * t339 - mrSges(4,3) * t296;
t232 = -mrSges(4,2) * t339 - mrSges(4,3) * t295;
t221 = mrSges(5,1) * t328 - mrSges(5,3) * t263;
t220 = -mrSges(5,2) * t328 + mrSges(5,3) * t262;
t216 = mrSges(4,1) * t295 + mrSges(4,2) * t296;
t203 = -t305 * t487 - t370;
t202 = t305 * t486 - t361 * t372;
t201 = -t304 * t487 - t368;
t200 = t304 * t486 - t369;
t162 = t357 * t449 + (t381 * qJD(2) + qJD(3) * t595) * t354;
t161 = t357 * t450 + (qJD(2) * t380 + qJD(3) * t379) * t354;
t157 = t229 * t350 + t303 * t349;
t151 = -mrSges(6,2) * t328 + mrSges(6,3) * t435;
t131 = -mrSges(5,2) * t282 + mrSges(5,3) * t170;
t130 = mrSges(5,1) * t282 - mrSges(5,3) * t169;
t125 = -t355 * t210 + t212 * t352;
t120 = t176 * t350 + t231 * t349;
t118 = t174 * t350 + t230 * t349;
t110 = -mrSges(5,1) * t170 + mrSges(5,2) * t169;
t87 = pkin(5) * t396 - pkin(11) * t435 + t532;
t84 = qJD(4) * t177 + t162 * t364 + t360 * t428;
t83 = -qJD(4) * t178 - t162 * t360 + t364 * t428;
t82 = qJD(6) * t394 - t150 * t359 + t363 * t450;
t81 = qJD(6) * t179 + t150 * t363 + t359 * t450;
t72 = -mrSges(6,2) * t282 + mrSges(6,3) * t102;
t66 = pkin(5) * t489 - t70;
t45 = t355 * t89 - t512;
t44 = t352 * t89 + t85;
t41 = t352 * t83 + t355 * t84;
t30 = -pkin(5) * t450 - t32;
t25 = -mrSges(7,2) * t101 + mrSges(7,3) * t50;
t24 = mrSges(7,1) * t101 - mrSges(7,3) * t49;
t21 = t359 * t87 + t363 * t45;
t20 = -t359 * t45 + t363 * t87;
t13 = -qJD(6) * t64 + t161 * t363 - t359 * t41;
t12 = qJD(6) * t63 + t161 * t359 + t363 * t41;
t10 = t49 * Ifges(7,4) + t50 * Ifges(7,2) + t101 * Ifges(7,6);
t3 = -pkin(5) * t282 - t5;
t7 = [(-m(6) * t5 + m(7) * t3 + t613) * (-t355 * t177 + t178 * t352) - t580 * t161 + (-m(6) * t38 + m(7) * t34 - t607) * (t352 * t84 - t355 * t83) + t288 * t428 + m(3) * (t357 ^ 2 * qJDD(1) + (t297 * t535 + t298 * t362) * t354) + m(7) * (t1 * t64 + t12 * t15 + t13 * t14 + t2 * t63) + m(6) * (t109 * t6 + t39 * t41) + m(5) * (t115 * t83 + t116 * t84 + t177 * t29 + t178 * t28) + m(4) * (t111 * t229 + t162 * t192 + t226 * t303 + t253 * t428) + t109 * t72 + t12 * t93 + t13 * t94 + t64 * t25 + t63 * t24 + (-m(4) * t112 + m(5) * t92 + m(6) * t57 + t110 - t233 + t46) * t228 + t41 * t151 + (-m(2) - m(3) - m(4) - m(5) + t618) * g(3) + t177 * t130 + t178 * t131 + t84 * t220 + t83 * t221 + t229 * t232 + (qJDD(2) * t459 - t366 * t488) * mrSges(3,1) + (-qJDD(2) * t488 - t366 * t459) * mrSges(3,2) + t162 * t287 + m(2) * qJDD(1) + t303 * t216; (-Ifges(7,4) * t394 + Ifges(7,2) * t179) * t570 + (-Ifges(7,5) * t394 + Ifges(7,6) * t179) * t562 + (-Ifges(7,1) * t394 + Ifges(7,4) * t179) * t571 + (t1 * t179 - t14 * t81 + t15 * t82 + t2 * t394) * mrSges(7,3) + t3 * (-mrSges(7,1) * t179 - mrSges(7,2) * t394) - t394 * t574 + (t447 * t601 + t591) * qJD(3) + t623 * t149 + (-m(4) * t417 - t201 * mrSges(4,1) - mrSges(4,3) * t499 - m(5) * (pkin(3) * t201 + t417) - (t201 * t364 + t304 * t492) * mrSges(5,1) - (-t201 * t360 + t304 * t490) * mrSges(5,2) + mrSges(3,1) * t371 + t304 * mrSges(3,2) + t618 * (-t200 * t358 + t201 * t348 + t304 * t413 - t293) + t587 * (t201 * t349 - t304 * t493) - t631 * (t201 * t350 + t304 * t494) + t620 * t200) * g(2) + (-m(4) * t416 - t203 * mrSges(4,1) - mrSges(4,3) * t498 - m(5) * (pkin(3) * t203 + t416) - (t203 * t364 + t305 * t492) * mrSges(5,1) - (-t203 * t360 + t305 * t490) * mrSges(5,2) + mrSges(3,1) * t372 + t305 * mrSges(3,2) + t618 * (-t202 * t358 + t203 * t348 + t305 * t413 - t294) + t587 * (t203 * t349 - t305 * t493) - t631 * (t203 * t350 + t305 * t494) + t620 * t202) * g(1) + (-m(5) * (pkin(3) * t265 + t478) - (t265 * t364 + t432) * mrSges(5,1) - (-t265 * t360 + t364 * t463) * mrSges(5,2) - m(4) * t478 - t265 * mrSges(4,1) - mrSges(4,3) * t463 - (mrSges(3,1) * t535 - mrSges(3,2) * t362) * t354 + t618 * (pkin(4) * t432 - t264 * t358 + t265 * t348 + t478) + t587 * (t265 * t349 - t350 * t463) - t631 * (t265 * t350 + t349 * t463) + t620 * t264) * g(3) + t616 * t93 + (t1 * t37 + t2 * t36 + t3 * t66 + (-t125 + t30) * t34 + t616 * t15 + t617 * t14) * m(7) + t617 * t94 + (Ifges(5,5) * t234 + Ifges(6,5) * t150 + Ifges(5,6) * t235 + t450 * t614) * t537 + (Ifges(5,5) * t309 + Ifges(6,5) * t215 + Ifges(5,6) * t308 - t489 * t614) * t539 + t356 * t609 + t599 * t292 + t602 * t221 + (t115 * t602 + t116 * t603 + t171 * t292 + t193 * t29 + t194 * t28 + t276 * t92) * m(5) + t603 * t220 + t580 * t249 + (Ifges(5,4) * t309 + Ifges(5,2) * t308 - Ifges(5,6) * t489) * t548 + (Ifges(5,1) * t309 + Ifges(5,4) * t308 - Ifges(5,5) * t489) * t549 + t234 * t550 + t235 * t551 + t150 * t556 + t607 * t125 + (t129 * t189 + t218 * t57 + t5 * t70 + t6 * t71 + t608 * t39 + (t125 + t32) * t38) * m(6) + t608 * t151 + t296 * (Ifges(4,5) * t356 + (t361 * Ifges(4,1) + t523) * t353) / 0.2e1 + (Ifges(7,4) * t81 + Ifges(7,2) * t82) * t554 + (t129 * t150 + t215 * t57 - t39 * t450 + t489 * t6) * mrSges(6,2) + (Ifges(7,5) * t81 + Ifges(7,6) * t82) * t546 + (Ifges(4,4) * t296 - Ifges(4,2) * t295 + Ifges(4,6) * t339) * t489 / 0.2e1 + (Ifges(4,1) * t296 - Ifges(4,4) * t295 + Ifges(4,5) * t339) * t447 + t5 * (-mrSges(6,1) * t489 - mrSges(6,3) * t215) + t29 * (-mrSges(5,1) * t489 - mrSges(5,3) * t309) + t28 * (mrSges(5,2) * t489 + mrSges(5,3) * t308) + (t425 + t297) * mrSges(3,1) - t216 * t533 + t226 * t388 + t356 * t460 / 0.2e1 + (Ifges(5,1) * t234 + Ifges(5,4) * t235 + Ifges(5,5) * t450) * t540 + t309 * t563 + t308 * t564 + t82 * t566 + t215 * t572 + (t333 - t298) * mrSges(3,2) - t288 * t429 - t356 * t610 - t239 * t450 / 0.2e1 + t116 * (-mrSges(5,2) * t450 + mrSges(5,3) * t235) - t594 * t489 / 0.2e1 + t596 * t287 + (t111 * t311 + t112 * t310 - t191 * t292 + t192 * t596 - t226 * t533 - t253 * t429) * m(4) + (Ifges(6,1) * t215 - Ifges(6,5) * t489) * t560 + (Ifges(6,1) * t150 + Ifges(6,5) * t450) * t542 + t81 * t53 / 0.2e1 + t34 * (-mrSges(7,1) * t82 + mrSges(7,2) * t81) + t30 * t76 + t71 * t72 + t70 * t73 + t66 * t16 + Ifges(3,3) * qJDD(2) + t36 * t24 + t37 * t25 + (Ifges(7,1) * t81 + Ifges(7,4) * t82) * t552 + t339 * (Ifges(4,3) * t356 + (Ifges(4,5) * t361 + Ifges(4,6) * t365) * t353) / 0.2e1 - t295 * (Ifges(4,6) * t356 + t387) / 0.2e1 + (Ifges(6,4) * t215 - Ifges(6,6) * t489) * t561 + (Ifges(6,4) * t150 + Ifges(6,6) * t450) * t544 + t32 * t152 + t179 * t10 / 0.2e1 + t189 * t114 + t193 * t130 + t194 * t131 + t218 * t46 + t262 * (Ifges(5,4) * t234 + Ifges(5,2) * t235 + Ifges(5,6) * t450) / 0.2e1 + t38 * (mrSges(6,1) * t450 - mrSges(6,3) * t150) + t115 * (mrSges(5,1) * t450 - mrSges(5,3) * t234) + t577 * t214 + t171 * (-mrSges(5,1) * t235 + mrSges(5,2) * t234) + t276 * t110 + (t111 * t489 - t112 * t491 - t191 * t449 - t192 * t450) * mrSges(4,3) + (t576 * qJD(2) * (-Ifges(4,2) * t361 + t523) + t353 * t240) * t474 / 0.2e1 + t92 * (-mrSges(5,1) * t308 + mrSges(5,2) * t309) + t310 * t233 + t311 * t232 + t376 * t469 / 0.2e1; -(t601 * t361 + (-Ifges(4,2) * t454 + t364 * t160 + t240 + t334) * t365 + t328 * (Ifges(5,3) * t361 + t365 * t401) + t263 * (Ifges(5,5) * t361 + t365 * t405) + t262 * (Ifges(5,6) * t361 + t365 * t403)) * t476 / 0.2e1 + (t57 * mrSges(6,2) - mrSges(6,3) * t5 + t3 * t406 + t400 * t562 + t402 * t570 + t404 * t571 + t441 * t53 + 0.2e1 * t572) * t315 + t622 * t237 + t623 * t306 + (t618 * (-t228 * t348 - t229 * t358) + t620 * t229 + t586 * t228) * g(3) + (t618 * (-t175 * t348 - t176 * t358) + t620 * t176 + t586 * t175) * g(1) + (t618 * (-t173 * t348 - t174 * t358) + t620 * t174 + t586 * t173) * g(2) + t613 * t244 + t611 * t94 + (t1 * t145 + t14 * t611 + t144 * t2 + t15 * t612 + t244 * t3 + t34 * t606) * m(7) + t612 * t93 + (t434 + t584) * t192 + (-pkin(3) * t92 - t115 * t140 - t116 * t141) * m(5) - (Ifges(6,1) * t542 + Ifges(6,4) * t544 + Ifges(6,5) * t537 + t462 + t556 + t589) * t307 + (-t1 * t497 + t14 * t597 - t15 * t598 - t2 * t496) * mrSges(7,3) + t600 * t114 + (Ifges(5,2) * t364 + t522) * t548 + (Ifges(5,1) * t360 + t521) * t549 + t427 * t551 - t238 * t557 + t604 * t151 + (t129 * t600 - t244 * t5 + t245 * t6 - t348 * t57 + t38 * t605 + t39 * t604) * m(6) + t605 * t152 + t606 * t76 + (-t159 / 0.2e1 - t516) * t473 - t10 * t497 / 0.2e1 + (-Ifges(7,4) * t390 - Ifges(7,2) * t391) * t554 + (Ifges(7,4) * t211 + Ifges(7,2) * t209) * t555 + (-Ifges(7,5) * t390 - Ifges(7,6) * t391) * t546 + (Ifges(7,5) * t211 + Ifges(7,6) * t209) * t547 + t460 - t366 * t376 / 0.2e1 + (t550 - t517) * t472 + (-Ifges(7,1) * t390 - Ifges(7,4) * t391) * t552 + (Ifges(7,1) * t211 + Ifges(7,4) * t209) * t553 + t609 + (-t115 * (mrSges(5,1) * t361 - mrSges(5,3) * t482) - t116 * (-mrSges(5,3) * t360 * t365 - mrSges(5,2) * t361)) * t476 + t144 * t24 + t145 * t25 + t328 * t171 * (mrSges(5,1) * t360 + mrSges(5,2) * t364) + (t129 * t238 + t39 * t454) * mrSges(6,2) + (-Ifges(6,1) * t238 + Ifges(6,5) * t454) * t543 + (-Ifges(6,5) * t238 + Ifges(6,3) * t454) * t538 + (-Ifges(6,4) * t238 + Ifges(6,6) * t454) * t545 - t38 * (mrSges(6,1) * t454 + mrSges(6,3) * t238) + (Ifges(5,5) * t360 + Ifges(5,6) * t364) * t539 + (t239 * t447 - t591) * qJD(2) + t360 * t563 + t364 * t564 + t211 * t565 + t496 * t574 + (t433 - t287) * t191 - t610 + (m(5) * ((-t115 * t364 - t116 * t360) * qJD(4) + t593) + t364 * t131 - t360 * t130 - t221 * t472 - t220 * t473) * pkin(10) + t593 * mrSges(5,3) - t598 * t52 / 0.2e1 + (mrSges(7,1) * t598 - mrSges(7,2) * t597) * t34 - pkin(3) * t110 + t92 * t409 - t348 * t46 - t141 * t220 - t140 * t221 + t577 * t314 + t245 * t72 + (t262 * t403 + t263 * t405 + t328 * t401) * qJD(4) / 0.2e1; t622 * t396 + (-m(6) * t420 + t118 * mrSges(6,2) - t624 * mrSges(5,1) - (-t174 * t364 - t230 * t360) * mrSges(5,2) - m(7) * (pkin(11) * t118 + t420) - t631 * (-t174 * t349 + t230 * t350)) * g(2) + (-m(6) * t419 + t120 * mrSges(6,2) - t625 * mrSges(5,1) - (-t176 * t364 - t231 * t360) * mrSges(5,2) - m(7) * (pkin(11) * t120 + t419) - t631 * (-t176 * t349 + t231 * t350)) * g(1) + (-m(6) * t418 + t157 * mrSges(6,2) - m(7) * (pkin(11) * t157 + t418) - mrSges(5,1) * t177 + mrSges(5,2) * t178 - t631 * (-t229 * t349 + t303 * t350)) * g(3) + (Ifges(6,1) * t543 + Ifges(6,4) * t545 + Ifges(6,5) * t538 + t400 * t547 + t402 * t555 + t404 * t553 - t393 + t557 - t589) * t435 + t594 + (Ifges(7,5) * t359 + Ifges(7,6) * t363) * t562 + t607 * t44 - t45 * t151 + t263 * t516 + t262 * t517 + (t393 + t462) * qJD(6) - t263 * (Ifges(5,1) * t262 - t515) / 0.2e1 + t73 * t530 + t72 * t531 + t10 * t536 + (Ifges(5,5) * t262 - Ifges(5,6) * t263) * t538 + t159 * t540 + (-g(1) * t120 - g(2) * t118 - g(3) * t157 + (-t471 + t506) * t15 + (-t470 + t505) * t14 + t592) * mrSges(7,3) + t505 * t565 + t506 * t566 + (Ifges(7,2) * t363 + t519) * t570 + (Ifges(7,1) * t359 + t518) * t571 + t359 * t574 - t114 * t532 + (-t94 * t470 - t93 * t471 + m(7) * ((-t14 * t363 - t15 * t359) * qJD(6) + t592) + t363 * t25 - t359 * t24) * (pkin(11) + t531) + t3 * t407 - t21 * t93 - t20 * t94 - t28 * mrSges(5,2) + t29 * mrSges(5,1) - t6 * mrSges(6,2) + t5 * mrSges(6,1) + t347 * t16 + (-t14 * t20 - t15 * t21 + t3 * t347 - t34 * t44) * m(7) + t52 * t441 - t115 * t220 + t116 * t221 - t171 * (mrSges(5,1) * t263 + mrSges(5,2) * t262) + (t146 * t402 + t147 * t404 + t181 * t400) * qJD(6) / 0.2e1 - (-Ifges(5,2) * t263 + t160 + t254) * t262 / 0.2e1 + (-t129 * t532 + t38 * t44 - t39 * t45 + (t352 * t6 + t355 * t5) * pkin(4)) * m(6); t363 * t24 + t359 * t25 + t607 * t396 + t398 * qJD(6) + (-t151 - t398) * t435 + t46 + (t1 * t359 - t396 * t34 + t2 * t363 + t382 + t181 * (-t14 * t359 + t15 * t363)) * m(7) + (t38 * t396 - t39 * t435 + t382 + t57) * m(6); -t34 * (mrSges(7,1) * t147 + mrSges(7,2) * t146) + (Ifges(7,1) * t146 - t520) * t553 + t52 * t552 + (Ifges(7,5) * t146 - Ifges(7,6) * t147) * t547 - t14 * t93 + t15 * t94 - g(1) * ((-t120 * t359 + t175 * t363) * mrSges(7,1) + (-t120 * t363 - t175 * t359) * mrSges(7,2)) - g(2) * ((-t118 * t359 + t173 * t363) * mrSges(7,1) + (-t118 * t363 - t173 * t359) * mrSges(7,2)) - g(3) * ((-t157 * t359 + t225) * mrSges(7,1) + (-t157 * t363 - t504) * mrSges(7,2)) + (t14 * t146 + t147 * t15) * mrSges(7,3) + t9 + (-Ifges(7,2) * t147 + t143 + t53) * t555 + t585;];
tau  = t7;
