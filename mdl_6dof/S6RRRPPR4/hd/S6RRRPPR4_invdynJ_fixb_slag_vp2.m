% Calculate vector of inverse dynamics joint torques for
% S6RRRPPR4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPPR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:31:52
% EndTime: 2019-03-09 15:33:01
% DurationCPUTime: 44.59s
% Computational Cost: add. (13160->890), mult. (28979->1152), div. (0->0), fcn. (20331->12), ass. (0->401)
t339 = sin(qJ(2));
t343 = cos(qJ(2));
t405 = pkin(2) * t339 - pkin(8) * t343;
t275 = t405 * qJD(1);
t342 = cos(qJ(3));
t338 = sin(qJ(3));
t455 = qJD(1) * t339;
t432 = t338 * t455;
t204 = pkin(7) * t432 + t342 * t275;
t461 = t342 * t343;
t376 = pkin(3) * t339 - qJ(4) * t461;
t336 = -qJ(4) - pkin(8);
t415 = qJD(3) * t336;
t644 = -qJD(1) * t376 - qJD(4) * t338 + t342 * t415 - t204;
t249 = t338 * t275;
t447 = qJD(4) * t342;
t467 = t339 * t342;
t469 = t338 * t343;
t643 = t249 + (-pkin(7) * t467 - qJ(4) * t469) * qJD(1) - t338 * t415 - t447;
t642 = mrSges(5,1) + mrSges(6,1);
t600 = mrSges(6,3) - mrSges(5,2);
t335 = sin(pkin(10));
t478 = cos(pkin(10));
t577 = t335 * t644 - t643 * t478;
t413 = t478 * t338;
t267 = t335 * t342 + t413;
t454 = qJD(1) * t343;
t212 = t267 * t454;
t241 = t267 * qJD(3);
t628 = t212 - t241;
t412 = t478 * t342;
t367 = -t335 * t338 + t412;
t357 = t343 * t367;
t213 = qJD(1) * t357;
t242 = t367 * qJD(3);
t641 = t213 - t242;
t308 = qJD(3) - t454;
t508 = t308 / 0.2e1;
t452 = qJD(2) * t342;
t273 = -t432 + t452;
t430 = t342 * t455;
t274 = qJD(2) * t338 + t430;
t368 = t335 * t273 + t274 * t478;
t520 = t368 / 0.2e1;
t187 = -t478 * t273 + t274 * t335;
t523 = t187 / 0.2e1;
t524 = -t187 / 0.2e1;
t597 = Ifges(6,4) + Ifges(5,5);
t599 = Ifges(5,1) + Ifges(6,1);
t322 = pkin(7) * t455;
t292 = -qJD(2) * pkin(2) + t322;
t364 = pkin(3) * t273 - qJD(4) - t292;
t406 = pkin(2) * t343 + pkin(8) * t339;
t284 = -pkin(1) - t406;
t258 = t284 * qJD(1);
t323 = pkin(7) * t454;
t293 = qJD(2) * pkin(8) + t323;
t194 = t342 * t258 - t293 * t338;
t152 = -qJ(4) * t274 + t194;
t136 = pkin(3) * t308 + t152;
t195 = t258 * t338 + t293 * t342;
t153 = qJ(4) * t273 + t195;
t472 = t335 * t153;
t66 = t136 * t478 - t472;
t371 = qJD(5) - t66;
t63 = -t308 * pkin(4) + t371;
t350 = qJ(5) * t368 + t364;
t76 = pkin(4) * t187 - t350;
t627 = -t364 * mrSges(5,2) + mrSges(6,2) * t63 - mrSges(5,3) * t66 - t76 * mrSges(6,3);
t640 = Ifges(5,4) * t524 + Ifges(6,5) * t523 + t597 * t508 + t599 * t520 + t627;
t491 = Ifges(3,4) * t339;
t506 = t343 / 0.2e1;
t633 = Ifges(3,2) * t506 + t491 / 0.2e1;
t509 = -t308 / 0.2e1;
t521 = -t368 / 0.2e1;
t598 = Ifges(5,4) - Ifges(6,5);
t596 = Ifges(6,2) + Ifges(5,3);
t624 = Ifges(4,3) + t596;
t587 = -qJ(5) * t455 + t577;
t638 = t597 * t521;
t578 = t643 * t335 + t478 * t644;
t414 = t478 * t153;
t67 = t335 * t136 + t414;
t64 = t308 * qJ(5) + t67;
t637 = Ifges(5,4) * t520 + Ifges(6,5) * t521 + Ifges(5,6) * t508 + Ifges(6,6) * t509 + (Ifges(5,2) + Ifges(6,3)) * t524 + mrSges(5,1) * t364 - mrSges(6,1) * t76 + mrSges(6,2) * t64 + mrSges(5,3) * t67;
t70 = t478 * t152 - t472;
t566 = -t70 + qJD(5);
t299 = qJD(6) - t308;
t510 = t299 / 0.2e1;
t337 = sin(qJ(6));
t341 = cos(qJ(6));
t611 = t187 * t337 + t341 * t368;
t528 = t611 / 0.2e1;
t119 = t187 * t341 - t337 * t368;
t530 = t119 / 0.2e1;
t636 = Ifges(7,5) * t528 + Ifges(7,6) * t530 + Ifges(7,3) * t510;
t445 = qJD(1) * qJD(2);
t278 = qJDD(1) * t339 + t343 * t445;
t180 = qJD(3) * t273 + qJDD(2) * t338 + t278 * t342;
t181 = -qJD(3) * t274 + qJDD(2) * t342 - t278 * t338;
t106 = t180 * t335 - t181 * t478;
t534 = -t106 / 0.2e1;
t108 = t180 * t478 + t335 * t181;
t532 = t108 / 0.2e1;
t277 = t343 * qJDD(1) - t339 * t445;
t265 = qJDD(3) - t277;
t515 = t265 / 0.2e1;
t635 = -t273 / 0.2e1;
t634 = -t274 / 0.2e1;
t263 = t277 * pkin(7);
t632 = -mrSges(5,3) - mrSges(6,2);
t595 = Ifges(5,6) - Ifges(6,6);
t537 = pkin(4) + pkin(5);
t438 = t537 * t339;
t631 = pkin(9) * t641 + qJD(1) * t438 - t578;
t630 = -pkin(9) * t628 + t587;
t602 = pkin(9) * t368;
t629 = -t602 + t566;
t431 = t338 * t454;
t450 = qJD(3) * t338;
t572 = -t323 + (-t431 + t450) * pkin(3);
t625 = Ifges(5,2) * t524 - Ifges(6,3) * t523 + t598 * t520 + t637;
t333 = qJ(3) + pkin(10);
t325 = sin(t333);
t326 = cos(t333);
t344 = cos(qJ(1));
t340 = sin(qJ(1));
t463 = t340 * t343;
t226 = t325 * t463 + t326 * t344;
t227 = -t325 * t344 + t326 * t463;
t381 = t226 * t337 + t227 * t341;
t570 = t226 * t341 - t227 * t337;
t623 = mrSges(7,1) * t570 - t381 * mrSges(7,2);
t379 = t325 * t337 + t326 * t341;
t380 = t325 * t341 - t326 * t337;
t400 = -mrSges(4,1) * t342 + mrSges(4,2) * t338;
t622 = m(4) * pkin(2) + t379 * mrSges(7,1) + t380 * mrSges(7,2) + t600 * t325 + t326 * t642 - t400;
t581 = -t187 * t598 + t308 * t597 + t368 * t599;
t621 = t581 / 0.2e1;
t620 = pkin(9) * t187;
t619 = qJ(5) * t641 - qJD(5) * t267 + t572;
t448 = qJD(3) * t342;
t451 = qJD(2) * t343;
t360 = t338 * t451 + t339 * t448;
t460 = t343 * t344;
t247 = -t338 * t460 + t340 * t342;
t618 = g(1) * t344 + g(2) * t340;
t591 = mrSges(3,2) * t343;
t401 = mrSges(3,1) * t339 + t591;
t614 = pkin(1) * t401 - t339 * (Ifges(3,1) * t343 - t491) / 0.2e1;
t613 = -Ifges(5,2) * t523 + Ifges(6,3) * t524 - t509 * t595 - t521 * t598 + t637;
t612 = Ifges(5,4) * t523 + Ifges(6,5) * t524 + t509 * t597 + t521 * t599 - t627;
t46 = -t308 * t537 + t371 - t602;
t52 = t64 + t620;
t12 = -t337 * t52 + t341 * t46;
t476 = qJDD(1) * pkin(1);
t196 = -pkin(2) * t277 - pkin(8) * t278 - t476;
t236 = qJDD(2) * pkin(8) + t263;
t83 = -qJD(3) * t195 + t342 * t196 - t236 * t338;
t50 = pkin(3) * t265 - qJ(4) * t180 - qJD(4) * t274 + t83;
t82 = t338 * t196 + t342 * t236 + t258 * t448 - t293 * t450;
t55 = qJ(4) * t181 + qJD(4) * t273 + t82;
t15 = -t335 * t55 + t478 * t50;
t373 = qJDD(5) - t15;
t8 = -t108 * pkin(9) - t265 * t537 + t373;
t16 = t335 * t50 + t478 * t55;
t11 = t265 * qJ(5) + t308 * qJD(5) + t16;
t9 = pkin(9) * t106 + t11;
t1 = qJD(6) * t12 + t337 * t8 + t341 * t9;
t13 = t337 * t46 + t341 * t52;
t2 = -qJD(6) * t13 - t337 * t9 + t341 * t8;
t610 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t14 = -t265 * pkin(4) + t373;
t607 = -t83 * mrSges(4,1) - t15 * mrSges(5,1) + t14 * mrSges(6,1) + t82 * mrSges(4,2) + t16 * mrSges(5,2) - t11 * mrSges(6,3);
t606 = t12 * mrSges(7,1) + t63 * mrSges(6,1) + t67 * mrSges(5,2) + Ifges(4,5) * t634 + Ifges(4,6) * t635 + t638 + t595 * t523 + t624 * t509 + Ifges(3,6) * qJD(2) / 0.2e1 + qJD(1) * t633 - t13 * mrSges(7,2) - t64 * mrSges(6,3) - t66 * mrSges(5,1) + t636;
t264 = t278 * pkin(7);
t237 = -qJDD(2) * pkin(2) + t264;
t135 = -t181 * pkin(3) + qJDD(4) + t237;
t21 = t106 * pkin(4) - t108 * qJ(5) - qJD(5) * t368 + t135;
t533 = t106 / 0.2e1;
t605 = mrSges(5,2) * t135 + mrSges(6,2) * t14 - mrSges(5,3) * t15 - mrSges(6,3) * t21 + Ifges(6,5) * t533 + 0.2e1 * t597 * t515 + 0.2e1 * t599 * t532 + (Ifges(5,4) + t598) * t534;
t24 = qJD(6) * t119 + t106 * t337 + t108 * t341;
t546 = t24 / 0.2e1;
t25 = -qJD(6) * t611 + t106 * t341 - t108 * t337;
t545 = t25 / 0.2e1;
t604 = -m(4) - m(3);
t603 = -m(7) - m(6);
t526 = t180 / 0.2e1;
t525 = t181 / 0.2e1;
t257 = qJDD(6) - t265;
t516 = t257 / 0.2e1;
t601 = -mrSges(3,3) + mrSges(2,2);
t287 = t336 * t338;
t289 = t336 * t342;
t201 = -t478 * t287 - t289 * t335;
t165 = -pkin(9) * t267 + t201;
t202 = t335 * t287 - t478 * t289;
t166 = -pkin(9) * t367 + t202;
t78 = t165 * t341 - t166 * t337;
t594 = qJD(6) * t78 + t337 * t631 + t341 * t630;
t79 = t165 * t337 + t166 * t341;
t593 = -qJD(6) * t79 - t337 * t630 + t341 * t631;
t589 = t537 * t628 - t619;
t588 = pkin(4) * t455 - t578;
t433 = t478 * pkin(3);
t316 = -t433 - pkin(4);
t307 = -pkin(5) + t316;
t498 = pkin(3) * t335;
t314 = qJ(5) + t498;
t220 = t307 * t341 - t314 * t337;
t69 = t152 * t335 + t414;
t58 = t69 + t620;
t586 = qJD(6) * t220 - t337 * t58 + t341 * t629;
t221 = t307 * t337 + t314 * t341;
t585 = -qJD(6) * t221 - t337 * t629 - t341 * t58;
t580 = -pkin(4) * t628 + t619;
t155 = mrSges(5,1) * t308 - mrSges(5,3) * t368;
t156 = -mrSges(6,1) * t308 + mrSges(6,2) * t368;
t579 = t155 - t156;
t262 = Ifges(4,4) * t273;
t171 = t274 * Ifges(4,1) + t308 * Ifges(4,5) + t262;
t321 = Ifges(3,4) * t454;
t576 = Ifges(3,1) * t455 + Ifges(3,5) * qJD(2) + t342 * t171 + t321;
t575 = -qJD(2) * mrSges(3,1) - mrSges(4,1) * t273 + mrSges(4,2) * t274 + mrSges(3,3) * t455;
t330 = t343 * pkin(4);
t477 = qJ(5) * t343;
t574 = t325 * t477 + t326 * t330;
t312 = pkin(7) * t461;
t219 = t338 * t284 + t312;
t449 = qJD(3) * t339;
t573 = t338 * t449 - t342 * t451;
t571 = t263 * t343 + t264 * t339;
t569 = -t338 * t83 + t342 * t82;
t568 = t632 * t339;
t565 = m(5) - t603;
t402 = t343 * mrSges(3,1) - mrSges(3,2) * t339;
t564 = t339 * mrSges(4,3) + mrSges(2,1) + t402;
t562 = Ifges(4,5) * t180 + Ifges(4,6) * t181 - t106 * t595 + t108 * t597 + t265 * t624;
t60 = -t187 * t537 + t350;
t561 = -mrSges(7,1) * t60 + mrSges(7,3) * t13;
t560 = mrSges(7,2) * t60 - mrSges(7,3) * t12;
t559 = m(7) * pkin(5) + t642;
t551 = mrSges(5,1) * t135 + mrSges(6,1) * t21 - mrSges(6,2) * t11 - mrSges(5,3) * t16 + 0.2e1 * Ifges(6,3) * t533 - t108 * Ifges(5,4) / 0.2e1 - t265 * Ifges(5,6) / 0.2e1 + Ifges(6,6) * t515 + (-t598 + Ifges(6,5)) * t532 + (-t534 + t533) * Ifges(5,2);
t549 = Ifges(7,4) * t546 + Ifges(7,2) * t545 + Ifges(7,6) * t516;
t548 = Ifges(7,1) * t546 + Ifges(7,4) * t545 + Ifges(7,5) * t516;
t547 = m(5) * pkin(3);
t487 = Ifges(7,4) * t611;
t41 = Ifges(7,2) * t119 + Ifges(7,6) * t299 + t487;
t542 = -t41 / 0.2e1;
t541 = t41 / 0.2e1;
t111 = Ifges(7,4) * t119;
t42 = Ifges(7,1) * t611 + Ifges(7,5) * t299 + t111;
t540 = -t42 / 0.2e1;
t539 = t42 / 0.2e1;
t538 = Ifges(4,1) * t526 + Ifges(4,4) * t525 + Ifges(4,5) * t515;
t531 = -t119 / 0.2e1;
t529 = -t611 / 0.2e1;
t512 = t274 / 0.2e1;
t511 = -t299 / 0.2e1;
t499 = pkin(3) * t274;
t497 = pkin(5) * t326;
t494 = g(3) * t339;
t328 = t339 * pkin(7);
t490 = Ifges(3,4) * t343;
t489 = Ifges(4,4) * t338;
t488 = Ifges(4,4) * t342;
t486 = t194 * mrSges(4,3);
t485 = t195 * mrSges(4,3);
t484 = t274 * Ifges(4,4);
t276 = t405 * qJD(2);
t453 = qJD(2) * t339;
t440 = pkin(7) * t453;
t458 = t342 * t276 + t338 * t440;
t110 = -t339 * t447 + t376 * qJD(2) + (-t312 + (qJ(4) * t339 - t284) * t338) * qJD(3) + t458;
t459 = t338 * t276 + t284 * t448;
t125 = (-pkin(7) * qJD(2) - qJ(4) * qJD(3)) * t467 + (-qJD(4) * t339 + (-pkin(7) * qJD(3) - qJ(4) * qJD(2)) * t343) * t338 + t459;
t49 = t335 * t110 + t478 * t125;
t471 = t336 * t339;
t470 = t338 * t339;
t468 = t338 * t344;
t466 = t339 * t344;
t465 = t340 * t338;
t319 = pkin(3) * t342 + pkin(2);
t295 = t343 * t319;
t269 = t342 * t284;
t191 = -qJ(4) * t467 + t269 + (-pkin(7) * t338 - pkin(3)) * t343;
t200 = -qJ(4) * t470 + t219;
t129 = t335 * t191 + t478 * t200;
t309 = pkin(3) * t470;
t457 = t339 * t326 * qJ(5) - t309;
t279 = t328 + t309;
t456 = t344 * pkin(1) + t340 * pkin(7);
t444 = Ifges(7,5) * t24 + Ifges(7,6) * t25 + Ifges(7,3) * t257;
t324 = pkin(7) * t451;
t439 = m(4) * pkin(8) + mrSges(4,3);
t437 = t336 * t466;
t331 = t344 * pkin(7);
t434 = pkin(3) * t468 + t340 * t471 + t331;
t209 = pkin(3) * t360 + t324;
t170 = t273 * Ifges(4,2) + t308 * Ifges(4,6) + t484;
t425 = -t338 * t170 / 0.2e1;
t45 = t106 * mrSges(5,1) + t108 * mrSges(5,2);
t75 = -t265 * mrSges(6,1) + t108 * mrSges(6,2);
t44 = t106 * mrSges(6,1) - t108 * mrSges(6,3);
t409 = t295 - t471;
t408 = pkin(3) * t465 + t319 * t460 + t456;
t407 = m(7) * (-pkin(9) - t336) - mrSges(7,3);
t7 = -t25 * mrSges(7,1) + t24 * mrSges(7,2);
t122 = t129 - t477;
t232 = t367 * t339;
t403 = qJ(5) * t232 - t279;
t399 = mrSges(4,1) * t338 + mrSges(4,2) * t342;
t228 = t325 * t460 - t340 * t326;
t229 = t340 * t325 + t326 * t460;
t150 = t228 * t341 - t229 * t337;
t151 = t228 * t337 + t229 * t341;
t396 = mrSges(7,1) * t150 - mrSges(7,2) * t151;
t395 = (mrSges(7,1) * t380 - mrSges(7,2) * t379) * t339;
t394 = Ifges(4,1) * t342 - t489;
t393 = Ifges(4,1) * t338 + t488;
t391 = -Ifges(4,2) * t338 + t488;
t390 = Ifges(4,2) * t342 + t489;
t389 = Ifges(3,5) * t343 - Ifges(3,6) * t339;
t388 = Ifges(4,5) * t342 - Ifges(4,6) * t338;
t387 = Ifges(4,5) * t338 + Ifges(4,6) * t342;
t386 = -qJ(5) * t187 - t499;
t128 = t191 * t478 - t335 * t200;
t124 = t330 - t128;
t71 = t343 * pkin(5) - t232 * pkin(9) + t124;
t231 = t267 * t339;
t77 = pkin(9) * t231 + t122;
t32 = -t337 * t77 + t341 * t71;
t33 = t337 * t71 + t341 * t77;
t84 = -mrSges(7,2) * t299 + mrSges(7,3) * t119;
t85 = mrSges(7,1) * t299 - mrSges(7,3) * t611;
t384 = -t337 * t85 + t341 * t84;
t383 = -t110 * t478 + t335 * t125;
t158 = t231 * t341 - t232 * t337;
t159 = t231 * t337 + t232 * t341;
t182 = -t267 * t337 - t341 * t367;
t183 = t267 * t341 - t337 * t367;
t378 = t247 * pkin(3);
t377 = qJ(5) * t267 + t319;
t39 = qJ(5) * t453 - qJD(5) * t343 + t49;
t375 = -pkin(4) * t326 - qJ(5) * t325 - t319;
t374 = t407 * t339;
t245 = t338 * t463 + t342 * t344;
t370 = t292 * t399;
t363 = t444 + t610;
t362 = t245 * pkin(3);
t359 = t229 * pkin(4) + t228 * qJ(5) + t408;
t358 = -g(1) * t228 - g(2) * t226 - t325 * t494;
t161 = qJD(2) * t357 - t267 * t449;
t356 = qJ(5) * t161 + qJD(5) * t232 - t209;
t353 = Ifges(4,5) * t339 + t343 * t394;
t352 = Ifges(4,6) * t339 + t343 * t391;
t351 = Ifges(4,3) * t339 + t343 * t388;
t286 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t454;
t259 = t399 * t339;
t248 = t342 * t460 + t465;
t246 = -t340 * t461 + t468;
t218 = -pkin(7) * t469 + t269;
t207 = mrSges(4,1) * t308 - mrSges(4,3) * t274;
t206 = -mrSges(4,2) * t308 + mrSges(4,3) * t273;
t205 = -pkin(7) * t430 + t249;
t177 = -pkin(4) * t367 - t377;
t160 = t335 * t573 - t412 * t449 - t413 * t451;
t157 = -mrSges(6,2) * t187 + mrSges(6,3) * t308;
t154 = -mrSges(5,2) * t308 - mrSges(5,3) * t187;
t149 = -qJD(3) * t219 + t458;
t148 = (-t339 * t452 - t343 * t450) * pkin(7) + t459;
t142 = t212 * t337 + t213 * t341;
t141 = t212 * t341 - t213 * t337;
t140 = t367 * t537 + t377;
t139 = pkin(4) * t231 - t403;
t138 = -mrSges(4,2) * t265 + mrSges(4,3) * t181;
t137 = mrSges(4,1) * t265 - mrSges(4,3) * t180;
t127 = mrSges(5,1) * t187 + mrSges(5,2) * t368;
t126 = mrSges(6,1) * t187 - mrSges(6,3) * t368;
t123 = -t231 * t537 + t403;
t114 = -mrSges(4,1) * t181 + mrSges(4,2) * t180;
t86 = pkin(4) * t368 - t386;
t80 = t180 * Ifges(4,4) + t181 * Ifges(4,2) + t265 * Ifges(4,6);
t74 = mrSges(5,1) * t265 - mrSges(5,3) * t108;
t73 = -mrSges(5,2) * t265 - mrSges(5,3) * t106;
t72 = -mrSges(6,2) * t106 + mrSges(6,3) * t265;
t62 = -t368 * t537 + t386;
t61 = -pkin(4) * t160 - t356;
t57 = -qJD(6) * t159 - t160 * t341 - t161 * t337;
t56 = qJD(6) * t158 - t160 * t337 + t161 * t341;
t51 = -mrSges(7,1) * t119 + mrSges(7,2) * t611;
t43 = -pkin(4) * t453 + t383;
t38 = t160 * t537 + t356;
t31 = -pkin(9) * t160 + t39;
t30 = -t161 * pkin(9) - qJD(2) * t438 + t383;
t20 = -mrSges(7,2) * t257 + mrSges(7,3) * t25;
t19 = mrSges(7,1) * t257 - mrSges(7,3) * t24;
t10 = pkin(5) * t106 + t21;
t4 = -qJD(6) * t33 + t30 * t341 - t31 * t337;
t3 = qJD(6) * t32 + t30 * t337 + t31 * t341;
t5 = [(Ifges(7,4) * t159 + Ifges(7,2) * t158) * t545 + (Ifges(7,4) * t56 + Ifges(7,2) * t57) * t530 - t614 * t445 - (t342 * t170 + t338 * t171) * t449 / 0.2e1 + (qJD(2) * t351 + t160 * t595 - t387 * t449) * t508 + (Ifges(7,1) * t159 + Ifges(7,4) * t158) * t546 + (Ifges(7,1) * t56 + Ifges(7,4) * t57) * t528 + t625 * t160 + t605 * t232 + (t194 * mrSges(4,1) - t195 * mrSges(4,2) + Ifges(5,6) * t524 + Ifges(6,6) * t523 + t596 * t508 + t597 * t520 - t606 - t636) * t453 + t444 * t506 + (qJD(2) * t353 - t393 * t449) * t512 + t339 * t394 * t526 + t467 * t538 + t56 * t539 + t57 * t541 + t159 * t548 + t158 * t549 + t148 * t206 + t149 * t207 + t209 * t127 - t286 * t440 + (t621 + t640) * t161 + (-t595 * t231 + t339 * t388) * t515 + t575 * t324 + t576 * t451 / 0.2e1 + t114 * t328 + (-m(5) * t434 - t246 * mrSges(4,1) + t381 * mrSges(7,1) - t245 * mrSges(4,2) + t570 * mrSges(7,2) + t603 * (-t227 * pkin(4) - qJ(5) * t226 + t434) + t601 * t344 + t604 * t331 + t559 * t227 + t600 * t226 + (m(3) * pkin(1) - m(4) * t284 - m(7) * pkin(9) * t339 - t565 * (-pkin(1) - t295) + t564 - t568) * t340) * g(1) + (Ifges(7,5) * t56 + Ifges(7,6) * t57) * t510 + (Ifges(7,5) * t159 + Ifges(7,6) * t158) * t516 + m(4) * (t148 * t195 + t149 * t194 + t218 * t83 + t219 * t82 + (t237 * t339 + t292 * t451) * pkin(7)) + (pkin(1) * mrSges(3,1) + 0.2e1 * t633) * t277 + t273 * (qJD(2) * t352 - t390 * t449) / 0.2e1 + (-g(1) * t339 * t340 + t1 * t158 - t12 * t56 + t13 * t57 - t159 * t2) * mrSges(7,3) + t551 * t231 + t425 * t451 - t80 * t470 / 0.2e1 + (-t248 * mrSges(4,1) - t247 * mrSges(4,2) - m(5) * (t408 - t437) - m(6) * (t359 - t437) - m(7) * t359 - t151 * mrSges(7,1) - t150 * mrSges(7,2) + t632 * t466 + t604 * t456 + t601 * t340 - t559 * t229 - t600 * t228 + (-m(4) * t406 - t374 - t564) * t344) * g(2) + t39 * t157 - t10 * (-mrSges(7,1) * t158 + mrSges(7,2) * t159) + t49 * t154 + t43 * t156 + t139 * t44 + t61 * t126 + t128 * t74 + t129 * t73 + t122 * t72 + t123 * t7 + t124 * t75 + Ifges(2,3) * qJDD(1) + t3 * t84 + t4 * t85 + t60 * (-mrSges(7,1) * t57 + mrSges(7,2) * t56) + t38 * t51 + t32 * t19 + t33 * t20 + t218 * t137 + t219 * t138 + (-Ifges(6,6) * t533 - Ifges(5,6) * t534 + (-Ifges(3,2) * t339 + t490) * t445 / 0.2e1 + Ifges(7,3) * t516 + Ifges(7,6) * t545 + Ifges(7,5) * t546 - Ifges(4,6) * t525 - Ifges(4,5) * t526 - t624 * t515 - t597 * t532 + mrSges(3,3) * t263 - t562 / 0.2e1 + t607 + t610) * t343 + (Ifges(3,4) * t506 + t339 * Ifges(3,1) + t490 / 0.2e1 + mrSges(3,3) * t328 - pkin(1) * mrSges(3,2)) * t278 + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t571) + t571 * mrSges(3,3) + t237 * t259 + (-mrSges(3,1) * t328 + Ifges(3,5) * t339 + 0.2e1 * Ifges(3,6) * t506 - pkin(7) * t591) * qJDD(2) + t279 * t45 + t292 * (mrSges(4,1) * t360 - mrSges(4,2) * t573) + (t194 * t573 - t195 * t360 - t467 * t83 - t470 * t82) * mrSges(4,3) + qJD(2) ^ 2 * t389 / 0.2e1 + m(7) * (t1 * t33 - t10 * t123 + t12 * t4 + t13 * t3 + t2 * t32 + t38 * t60) + m(6) * (t11 * t122 + t124 * t14 + t139 * t21 + t39 * t64 + t43 * t63 + t61 * t76) - t383 * t155 + m(5) * (t128 * t15 + t129 * t16 + t135 * t279 - t209 * t364 - t383 * t66 + t49 * t67) + t402 * t476 + t339 * t391 * t525; (t75 - t74) * t201 + (t370 + t425) * qJD(3) + (Ifges(3,2) * t454 / 0.2e1 + t606 + Ifges(5,6) * t523 + Ifges(6,6) * t524 - Ifges(7,5) * t529 - Ifges(7,6) * t531 - Ifges(7,3) * t511 + t596 * t509 + t638) * t455 + t612 * t213 + t613 * t212 + t614 * qJD(1) ^ 2 + (Ifges(7,4) * t142 + Ifges(7,2) * t141) * t531 + (Ifges(7,5) * t142 + Ifges(7,6) * t141) * t511 + (t72 + t73) * t202 + (-t595 * t508 - t625) * t241 + (Ifges(7,1) * t528 + Ifges(7,4) * t530 + Ifges(7,5) * t510 + t539 + t560) * (qJD(6) * t182 + t241 * t337 + t242 * t341) + (Ifges(7,4) * t528 + Ifges(7,2) * t530 + Ifges(7,6) * t510 + t541 + t561) * (-qJD(6) * t183 + t241 * t341 - t242 * t337) + t605 * t267 + (-t195 * (-mrSges(4,2) * t339 - mrSges(4,3) * t469) - t194 * (mrSges(4,1) * t339 - mrSges(4,3) * t461)) * qJD(1) + (t1 * t182 + t12 * t142 - t13 * t141 - t183 * t2) * mrSges(7,3) + (-pkin(2) * t237 - t194 * t204 - t195 * t205 - t292 * t323) * m(4) + t587 * t157 + t588 * t156 + (t11 * t202 + t14 * t201 + t177 * t21 + t580 * t76 + t587 * t64 + t588 * t63) * m(6) + t589 * t51 + t387 * t515 + (Ifges(7,5) * t183 + Ifges(7,6) * t182) * t516 + (-t486 + t171 / 0.2e1) * t448 - t450 * t485 + t338 * t538 + t142 * t540 + t141 * t542 + (Ifges(7,4) * t183 + Ifges(7,2) * t182) * t545 + (Ifges(7,1) * t183 + Ifges(7,4) * t182) * t546 + t183 * t548 + t182 * t549 - t205 * t206 - t204 * t207 + t286 * t322 + t640 * t242 + t170 * t431 / 0.2e1 - t575 * t323 + t577 * t154 + t578 * t155 + t580 * t126 + (-t339 * t439 - m(7) * (t295 + t574) - t374 - m(6) * (t409 + t574) - m(5) * t409 - t402 + (-m(7) * t497 - t622) * t343 + t568) * g(3) + (t273 * t391 + t274 * t394 + t308 * t388) * qJD(3) / 0.2e1 - (t273 * t352 + t274 * t353 + t308 * t351) * qJD(1) / 0.2e1 + t581 * (t242 / 0.2e1 - t213 / 0.2e1) - t370 * t454 - t389 * t445 / 0.2e1 - t10 * (-mrSges(7,1) * t182 + mrSges(7,2) * t183) + t177 * t44 + t390 * t525 + t393 * t526 - (t321 + t576) * t454 / 0.2e1 + (Ifges(7,1) * t142 + Ifges(7,4) * t141) * t529 + t618 * (t401 + (-t407 - t439 + (m(5) + m(6)) * t336 + t632) * t343 + (-m(7) * (t375 - t497) + m(5) * t319 - m(6) * t375 + t622) * t339) + t140 * t7 - t60 * (-mrSges(7,1) * t141 + mrSges(7,2) * t142) - pkin(2) * t114 + t78 * t19 + t79 * t20 + Ifges(3,3) * qJDD(2) + t342 * t80 / 0.2e1 + (-t206 * t450 - t207 * t448 + m(4) * ((-t194 * t342 - t195 * t338) * qJD(3) + t569) + t342 * t138 - t338 * t137) * pkin(8) + t569 * mrSges(4,3) + t572 * t127 + t593 * t85 + t594 * t84 + (t1 * t79 - t10 * t140 + t593 * t12 + t594 * t13 + t2 * t78 + t589 * t60) * m(7) - t263 * mrSges(3,2) - t264 * mrSges(3,1) + Ifges(3,6) * t277 + Ifges(3,5) * t278 - t319 * t45 + t237 * t400 + (-t135 * t319 - t15 * t201 + t16 * t202 - t364 * t572 + t577 * t67 + t578 * t66) * m(5) - (-t515 * t595 + t551) * t367; (m(5) * t362 + mrSges(4,1) * t245 - mrSges(4,2) * t246 + t603 * (-t226 * pkin(4) + t227 * qJ(5) - t362) - t600 * t227 + t559 * t226 + t623) * g(2) + t613 * t368 + (t621 - t612) * t187 + t585 * t85 + (t1 * t221 + t12 * t585 + t13 * t586 + t2 * t220 - t60 * t62) * m(7) + t586 * t84 + (-m(7) * (-t325 * t438 + t457) + t395 - m(6) * t457 - (t326 * mrSges(6,3) + (-m(6) * pkin(4) - mrSges(6,1)) * t325) * t339 + t259) * g(3) + (Ifges(4,5) * t273 - Ifges(4,6) * t274) * t509 + t170 * t512 - (-mrSges(5,1) * t325 - mrSges(5,2) * t326 - t338 * t547) * t494 - t127 * t499 + (t15 * t478 + t16 * t335) * t547 - t194 * t206 + t195 * t207 + (m(5) * t66 + t579) * t69 + t274 * t485 + t273 * t486 + t73 * t498 - (Ifges(7,1) * t529 + Ifges(7,4) * t531 + Ifges(7,5) * t511 + t540 - t560) * t119 + (Ifges(4,1) * t273 - t484) * t634 + t562 + (-Ifges(4,2) * t274 + t171 + t262) * t635 + t74 * t433 - t70 * t154 - t86 * t126 - t62 * t51 + t220 * t19 + t221 * t20 + (t11 * t314 + t14 * t316 + t566 * t64 - t63 * t69 - t76 * t86) * m(6) + t566 * t157 + (Ifges(7,4) * t529 + Ifges(7,2) * t531 + Ifges(7,6) * t511 + t542 - t561) * t611 + (-m(5) * t378 - mrSges(4,1) * t247 + mrSges(4,2) * t248 + t396 + t603 * (-t228 * pkin(4) + t229 * qJ(5) + t378) - t600 * t229 + t559 * t228) * g(1) - t363 - t292 * (mrSges(4,1) * t274 + mrSges(4,2) * t273) + t314 * t72 + t316 * t75 - t607 - m(5) * (-t364 * t499 + t67 * t70); -t611 * t85 + t119 * t84 - (-t154 - t157) * t187 + t579 * t368 - t7 + t44 + t45 + (t119 * t13 - t12 * t611 + t10) * m(7) + (t187 * t64 - t368 * t63 + t21) * m(6) + (t187 * t67 + t368 * t66 + t135) * m(5) + (t343 * g(3) - t339 * t618) * t565; t341 * t19 + t337 * t20 + (t126 - t51) * t368 + t384 * qJD(6) + (-t157 - t384) * t308 + t75 + (t1 * t337 - t368 * t60 + t2 * t341 + t358 + t299 * (-t12 * t337 + t13 * t341)) * m(7) + (-t308 * t64 + t368 * t76 + t14 + t358) * m(6); -t60 * (mrSges(7,1) * t611 + mrSges(7,2) * t119) + (Ifges(7,1) * t119 - t487) * t529 + t41 * t528 + (Ifges(7,5) * t119 - Ifges(7,6) * t611) * t511 - t12 * t84 + t13 * t85 - g(1) * t396 - g(2) * t623 - g(3) * t395 + (t119 * t12 + t13 * t611) * mrSges(7,3) + t363 + (-Ifges(7,2) * t611 + t111 + t42) * t531;];
tau  = t5;
