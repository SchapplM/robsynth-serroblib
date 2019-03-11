% Calculate vector of inverse dynamics joint torques for
% S6RRRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
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
% Datum: 2019-03-09 16:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPPR8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR8_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR8_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR8_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR8_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR8_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:03:37
% EndTime: 2019-03-09 16:04:51
% DurationCPUTime: 47.93s
% Computational Cost: add. (11026->959), mult. (26423->1211), div. (0->0), fcn. (20076->10), ass. (0->431)
t534 = pkin(9) - qJ(5);
t591 = -mrSges(6,3) + mrSges(5,2) + mrSges(4,3);
t629 = -m(7) - m(6);
t636 = t534 * t629 + mrSges(3,2) - t591;
t321 = sin(qJ(3));
t324 = cos(qJ(3));
t502 = cos(pkin(6));
t425 = t502 * qJD(1);
t364 = t425 + qJD(2);
t318 = sin(pkin(6));
t322 = sin(qJ(2));
t483 = qJD(1) * t322;
t449 = t318 * t483;
t207 = t321 * t449 - t324 * t364;
t325 = cos(qJ(2));
t484 = qJD(1) * t318;
t448 = t325 * t484;
t275 = -qJD(3) + t448;
t320 = sin(qJ(6));
t323 = cos(qJ(6));
t144 = -t207 * t320 + t275 * t323;
t447 = t324 * t483;
t208 = t318 * t447 + t321 * t364;
t202 = qJD(6) + t208;
t145 = t207 * t323 + t275 * t320;
t521 = Ifges(7,4) * t145;
t55 = Ifges(7,2) * t144 + Ifges(7,6) * t202 + t521;
t657 = t55 / 0.2e1;
t470 = qJD(1) * qJD(2);
t439 = t325 * t470;
t352 = qJDD(1) * t322 + t439;
t339 = t352 * t318;
t348 = qJD(3) * t364;
t420 = t502 * qJDD(1);
t361 = t420 + qJDD(2);
t496 = t318 * t322;
t462 = t321 * t496;
t412 = qJD(3) * t462;
t116 = qJD(1) * t412 - t321 * t361 + (-t339 - t348) * t324;
t564 = t116 / 0.2e1;
t117 = t321 * t348 - t324 * t361 + (qJD(3) * t447 + t321 * t352) * t318;
t563 = -t117 / 0.2e1;
t482 = qJD(2) * t322;
t446 = t318 * t482;
t469 = qJDD(1) * t318;
t241 = -qJD(1) * t446 + t325 * t469;
t230 = qJDD(3) - t241;
t549 = -t230 / 0.2e1;
t656 = -t323 / 0.2e1;
t410 = pkin(1) * t425;
t234 = -pkin(8) * t449 + t325 * t410;
t360 = (pkin(2) * t322 - pkin(9) * t325) * t318;
t235 = qJD(1) * t360;
t138 = -t321 * t234 + t235 * t324;
t476 = qJD(3) * t324;
t655 = -qJD(5) * t321 + t476 * t534 + t138;
t413 = t324 * t448;
t638 = -t321 * qJD(4) + (t413 - t476) * qJ(4);
t565 = -t116 / 0.2e1;
t562 = t117 / 0.2e1;
t548 = t230 / 0.2e1;
t624 = Ifges(4,1) + Ifges(5,1);
t623 = -Ifges(4,4) + Ifges(5,5);
t622 = Ifges(5,4) + Ifges(4,5);
t654 = Ifges(6,4) - Ifges(5,5);
t620 = Ifges(6,5) - Ifges(5,6);
t619 = Ifges(5,2) + Ifges(4,3);
t618 = Ifges(4,6) - Ifges(5,6);
t617 = -Ifges(5,3) - Ifges(6,1);
t453 = pkin(1) * t502;
t305 = t322 * t453;
t567 = -pkin(4) - pkin(10);
t471 = pkin(3) - t567;
t354 = pkin(5) * t324 - t321 * t471;
t495 = t318 * t325;
t653 = qJD(3) * t354 - (-t305 + (-pkin(8) + t354) * t495) * qJD(1) - t638;
t421 = t471 * t322;
t490 = t324 * t325;
t463 = qJ(5) * t490;
t652 = -(-t421 - t463) * t484 + t655;
t139 = t324 * t234 + t321 * t235;
t121 = qJ(4) * t449 + t139;
t477 = qJD(3) * t321;
t468 = pkin(9) * t477;
t650 = -t468 - t121;
t487 = pkin(8) * t495 + t305;
t226 = t502 * pkin(9) + t487;
t185 = qJD(2) * pkin(9) + qJD(1) * t226;
t191 = (-pkin(2) * t325 - pkin(9) * t322 - pkin(1)) * t484;
t489 = t321 * t185 - t324 * t191;
t515 = t489 * mrSges(4,3);
t194 = t208 * qJ(5);
t440 = qJD(4) + t489;
t411 = -t194 + t440;
t568 = pkin(3) + pkin(4);
t57 = t275 * t568 + t411;
t184 = -t364 * pkin(2) - t234;
t83 = t207 * pkin(3) - t208 * qJ(4) + t184;
t359 = qJD(5) - t83;
t35 = pkin(5) * t208 + t207 * t567 + t359;
t47 = t275 * t471 + t411;
t13 = -t320 * t47 + t323 * t35;
t14 = t320 * t35 + t323 * t47;
t582 = t13 * mrSges(7,1) - t14 * mrSges(7,2);
t198 = Ifges(4,4) * t207;
t518 = Ifges(5,5) * t207;
t590 = t145 * Ifges(7,5) + t144 * Ifges(7,6) + t202 * Ifges(7,3) + t208 * t624 - t275 * t622 - t198 + t518;
t88 = pkin(3) * t275 + t440;
t197 = Ifges(6,4) * t207;
t98 = -t208 * Ifges(6,2) + t275 * Ifges(6,6) + t197;
t649 = t515 + t582 - t98 / 0.2e1 + t88 * mrSges(5,2) + t590 / 0.2e1 - t57 * mrSges(6,3);
t141 = Ifges(7,4) * t144;
t56 = t145 * Ifges(7,1) + t202 * Ifges(7,5) + t141;
t648 = t320 * t657 + t56 * t656;
t363 = qJD(2) * t410;
t399 = pkin(1) * t420;
t146 = pkin(8) * t241 + t322 * t399 + t325 * t363;
t129 = pkin(9) * t361 + t146;
t137 = -t241 * pkin(2) + (-qJDD(1) * pkin(1) - pkin(9) * t352) * t318;
t31 = t324 * t129 + t321 * t137 - t185 * t477 + t191 * t476;
t21 = t230 * qJ(4) - t275 * qJD(4) + t31;
t12 = -qJ(5) * t117 - qJD(5) * t207 - t21;
t647 = t654 * t564 - Ifges(4,4) * t565 - t617 * t562 - (Ifges(4,6) + t618) * t548 - t12 * mrSges(6,3) + (t617 - (2 * Ifges(4,2))) * t563 + (t620 + Ifges(6,5)) * t549;
t32 = -t321 * t129 + t137 * t324 - t185 * t476 - t191 * t477;
t349 = qJDD(4) - t32;
t331 = qJ(5) * t116 - qJD(5) * t208 + t349;
t11 = -t230 * t568 + t331;
t112 = qJDD(6) - t116;
t566 = t112 / 0.2e1;
t46 = -qJD(6) * t145 - t117 * t320 - t230 * t323;
t570 = t46 / 0.2e1;
t45 = qJD(6) * t144 + t117 * t323 - t230 * t320;
t571 = t45 / 0.2e1;
t297 = pkin(8) * t496;
t147 = -t318 * pkin(8) * t439 - qJDD(1) * t297 - t322 * t363 + t325 * t399;
t130 = -t361 * pkin(2) - t147;
t23 = t117 * pkin(3) + t116 * qJ(4) - t208 * qJD(4) + t130;
t335 = qJDD(5) - t23;
t5 = -pkin(5) * t116 + t117 * t567 + t335;
t9 = -t230 * t471 + t331;
t1 = qJD(6) * t13 + t320 * t5 + t323 * t9;
t2 = -qJD(6) * t14 - t320 * t9 + t323 * t5;
t585 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t6 = Ifges(7,5) * t45 + Ifges(7,6) * t46 + Ifges(7,3) * t112;
t646 = -t11 * mrSges(6,3) + t585 + Ifges(6,4) * t563 + Ifges(7,5) * t571 + Ifges(7,6) * t570 + Ifges(7,3) * t566 + 0.2e1 * t548 * t622 + 0.2e1 * t565 * t624 + (t548 - t549) * Ifges(6,6) + (-t564 + t565) * Ifges(6,2) + t6 / 0.2e1 + t623 * t562;
t513 = t208 * Ifges(4,4);
t100 = -t207 * Ifges(4,2) - t275 * Ifges(4,6) + t513;
t114 = t324 * t185 + t321 * t191;
t514 = t114 * mrSges(4,3);
t196 = Ifges(5,5) * t208;
t512 = t208 * Ifges(6,4);
t612 = -t207 * t617 + t275 * t620 + t196 - t512;
t264 = t275 * qJ(4);
t70 = qJ(5) * t207 + t114;
t64 = t264 - t70;
t90 = -t264 + t114;
t645 = -t100 / 0.2e1 - t514 + t612 / 0.2e1 - t64 * mrSges(6,3) - t90 * mrSges(5,2);
t644 = t241 / 0.2e1;
t643 = t339 / 0.2e1;
t642 = t361 / 0.2e1;
t640 = qJ(5) * t477 - qJD(5) * t324 + t650;
t533 = mrSges(5,2) * t207;
t153 = -mrSges(5,3) * t275 - t533;
t531 = mrSges(6,3) * t207;
t148 = mrSges(6,1) * t275 - t531;
t68 = -mrSges(7,1) * t144 + mrSges(7,2) * t145;
t503 = t148 - t68;
t639 = t153 - t503;
t388 = t323 * mrSges(7,1) - t320 * mrSges(7,2);
t390 = mrSges(6,1) * t321 - mrSges(6,2) * t324;
t391 = mrSges(5,1) * t324 + mrSges(5,3) * t321;
t393 = mrSges(4,1) * t324 - mrSges(4,2) * t321;
t466 = m(7) * pkin(10) + mrSges(7,3);
t637 = -t321 * (m(7) * pkin(5) + t388) - t324 * t466 - t390 - t391 - t393;
t62 = -pkin(4) * t207 + t359;
t635 = t184 * mrSges(4,1) + t83 * mrSges(5,1) + t62 * mrSges(6,2);
t634 = t62 * mrSges(6,1) + t184 * mrSges(4,2) - t83 * mrSges(5,3);
t26 = -pkin(3) * t230 + t349;
t631 = -t32 * mrSges(4,1) + t26 * mrSges(5,1) + t12 * mrSges(6,1) + t31 * mrSges(4,2) - t11 * mrSges(6,2) - t21 * mrSges(5,3);
t319 = qJ(4) + pkin(5);
t628 = m(7) * t319;
t626 = -mrSges(5,1) + mrSges(6,2);
t311 = t321 * qJ(4);
t269 = -t324 * pkin(3) - pkin(2) - t311;
t261 = t324 * pkin(4) - t269;
t201 = pkin(5) * t321 + pkin(10) * t324 + t261;
t278 = t534 * t321;
t142 = t201 * t323 - t278 * t320;
t615 = qJD(6) * t142 + t320 * t653 + t323 * t652;
t143 = t201 * t320 + t278 * t323;
t614 = -qJD(6) * t143 - t320 * t652 + t323 * t653;
t613 = -t207 * t618 + t208 * t622 - t275 * t619;
t527 = Ifges(3,4) * t322;
t334 = Ifges(3,6) * t502 + (t325 * Ifges(3,2) + t527) * t318;
t611 = t207 * Ifges(6,5) + Ifges(3,6) * qJD(2) - t208 * Ifges(6,6) + t275 * Ifges(6,3) + qJD(1) * t334;
t493 = t321 * t325;
t610 = -(pkin(5) * t322 + qJ(5) * t493) * t484 + t640;
t609 = qJ(5) * t321 * t448 - t640;
t608 = -(-t322 * t568 - t463) * t484 + t655;
t189 = (-t320 * t493 - t322 * t323) * t484;
t472 = qJD(6) * t324;
t607 = t320 * t477 - t323 * t472 + t189;
t491 = t323 * t325;
t190 = (-t320 * t322 + t321 * t491) * t484;
t606 = -t320 * t472 - t323 * t477 + t190;
t465 = t568 * t321;
t605 = -qJD(3) * t465 - (-t305 + (-pkin(8) - t465) * t495) * qJD(1) - t638;
t541 = cos(qJ(1));
t396 = t502 * t541;
t540 = sin(qJ(1));
t254 = t322 * t540 - t325 * t396;
t498 = t254 * t324;
t604 = -pkin(3) * t498 - t254 * t311;
t395 = t502 * t540;
t256 = t322 * t541 + t325 * t395;
t497 = t256 * t324;
t603 = -pkin(3) * t497 - t256 * t311;
t253 = t321 * t502 + t324 * t496;
t481 = qJD(2) * t325;
t445 = t318 * t481;
t168 = qJD(3) * t253 + t321 * t445;
t169 = -t412 + (qJD(3) * t502 + t445) * t324;
t239 = t487 * qJD(2);
t63 = t168 * pkin(3) - t169 * qJ(4) - t253 * qJD(4) + t239;
t252 = -t324 * t502 + t462;
t244 = t252 * pkin(3);
t422 = t253 * qJ(4) - t244;
t602 = pkin(3) * t477 - (t305 + (pkin(3) * t321 + pkin(8)) * t495) * qJD(1) + t638;
t601 = -t321 * t618 + t324 * t622;
t516 = Ifges(5,5) * t324;
t522 = Ifges(6,4) * t324;
t600 = -t321 * t617 + t516 - t522;
t517 = Ifges(5,5) * t321;
t525 = Ifges(4,4) * t321;
t599 = t324 * t624 + t517 - t525;
t598 = -t116 * t622 - t117 * t618 + t230 * t619;
t24 = mrSges(7,1) * t112 - mrSges(7,3) * t45;
t25 = -mrSges(7,2) * t112 + mrSges(7,3) * t46;
t597 = -t320 * t24 + t323 * t25;
t596 = t31 * t324 - t32 * t321;
t595 = t21 * t324 + t26 * t321;
t576 = t318 ^ 2;
t593 = (t322 * (Ifges(3,1) * t325 - t527) / 0.2e1 - pkin(1) * (mrSges(3,1) * t322 + mrSges(3,2) * t325)) * t576;
t592 = mrSges(4,2) - mrSges(6,1) - mrSges(5,3);
t589 = -t62 * (mrSges(6,1) * t324 + mrSges(6,2) * t321) - t83 * (mrSges(5,1) * t321 - mrSges(5,3) * t324) - t184 * (mrSges(4,1) * t321 + mrSges(4,2) * t324);
t337 = t388 + t628;
t587 = -t337 + t592;
t419 = mrSges(3,3) * t449;
t584 = -m(4) * t184 + mrSges(3,1) * t364 - mrSges(4,1) * t207 - mrSges(4,2) * t208 - t419;
t387 = t320 * mrSges(7,1) + t323 * mrSges(7,2);
t583 = t387 + t636;
t581 = mrSges(3,1) - t637;
t580 = t592 - t628;
t579 = mrSges(4,1) + t466 - t626;
t555 = -t202 / 0.2e1;
t559 = -t145 / 0.2e1;
t561 = -t144 / 0.2e1;
t578 = -Ifges(7,5) * t559 - Ifges(7,6) * t561 - Ifges(7,3) * t555 + t582;
t373 = Ifges(7,5) * t323 - Ifges(7,6) * t320;
t519 = Ifges(7,4) * t323;
t378 = -Ifges(7,2) * t320 + t519;
t520 = Ifges(7,4) * t320;
t383 = Ifges(7,1) * t323 - t520;
t58 = -pkin(5) * t275 - t64;
t577 = t373 * t555 + t378 * t561 + t383 * t559 - t387 * t58 + t648;
t7 = Ifges(7,4) * t45 + Ifges(7,2) * t46 + Ifges(7,6) * t112;
t574 = t7 / 0.2e1;
t8 = Ifges(7,1) * t45 + Ifges(7,4) * t46 + Ifges(7,5) * t112;
t573 = -t8 / 0.2e1;
t560 = t144 / 0.2e1;
t558 = t145 / 0.2e1;
t554 = t202 / 0.2e1;
t553 = -t207 / 0.2e1;
t552 = t207 / 0.2e1;
t551 = -t208 / 0.2e1;
t550 = t208 / 0.2e1;
t545 = -t275 / 0.2e1;
t544 = t275 / 0.2e1;
t539 = pkin(1) * t318;
t538 = pkin(9) * t256;
t536 = pkin(9) * t324;
t243 = t252 * pkin(4);
t535 = t254 * pkin(9);
t532 = mrSges(5,2) * t208;
t530 = mrSges(6,3) * t208;
t529 = mrSges(7,3) * t320;
t528 = mrSges(7,3) * t323;
t526 = Ifges(3,4) * t325;
t524 = Ifges(4,4) * t324;
t523 = Ifges(6,4) * t321;
t501 = qJ(5) * t252;
t500 = t254 * t320;
t499 = t254 * t323;
t494 = t320 * t324;
t492 = t323 * t324;
t132 = t208 * pkin(3) + t207 * qJ(4);
t488 = pkin(2) * t495 + pkin(9) * t496;
t227 = -t488 - t539;
t136 = t324 * t226 + t321 * t227;
t259 = t325 * t453 - t297;
t451 = t318 * t540;
t485 = t541 * pkin(1) + pkin(8) * t451;
t475 = qJD(4) * t325;
t474 = qJD(6) * t320;
t473 = qJD(6) * t323;
t467 = pkin(9) * t476;
t464 = qJ(4) * t495;
t461 = t318 * t490;
t455 = Ifges(3,5) * t339 + Ifges(3,6) * t241 + Ifges(3,3) * t361;
t257 = -t322 * t395 + t325 * t541;
t454 = t257 * pkin(2) + t485;
t225 = -t502 * pkin(2) - t259;
t452 = t318 * t541;
t444 = t318 * t475;
t443 = t496 / 0.2e1;
t438 = -t484 / 0.2e1;
t437 = t484 / 0.2e1;
t245 = t254 * pkin(2);
t255 = t322 * t396 + t325 * t540;
t428 = pkin(9) * t255 - t245;
t247 = t256 * pkin(2);
t427 = pkin(9) * t257 - t247;
t51 = -t116 * mrSges(6,1) + t117 * mrSges(6,2);
t426 = t253 * mrSges(6,1) + t252 * mrSges(6,2);
t172 = t255 * t321 + t324 * t452;
t159 = t172 * pkin(3);
t173 = t255 * t324 - t321 * t452;
t424 = qJ(4) * t173 - t159;
t176 = t257 * t321 - t324 * t451;
t163 = t176 * pkin(3);
t177 = t257 * t324 + t321 * t451;
t423 = qJ(4) * t177 - t163;
t135 = -t321 * t226 + t227 * t324;
t418 = mrSges(3,3) * t448;
t417 = t177 * pkin(3) + t454;
t414 = pkin(3) * t461 + t321 * t464 + t488;
t402 = t325 * t437;
t397 = -pkin(1) * t540 + pkin(8) * t452;
t120 = pkin(3) * t495 - t135;
t394 = mrSges(4,1) * t252 + mrSges(4,2) * t253;
t392 = t252 * mrSges(5,1) - t253 * mrSges(5,3);
t170 = -t252 * t320 + t318 * t491;
t171 = t252 * t323 + t320 * t495;
t389 = mrSges(7,1) * t170 - mrSges(7,2) * t171;
t382 = Ifges(7,1) * t320 + t519;
t381 = -Ifges(4,2) * t321 + t524;
t379 = -Ifges(6,2) * t324 + t523;
t377 = Ifges(7,2) * t323 + t520;
t374 = Ifges(6,5) * t321 - Ifges(6,6) * t324;
t372 = Ifges(7,5) * t320 + Ifges(7,6) * t323;
t371 = t13 * t323 + t14 * t320;
t118 = t225 - t422;
t53 = pkin(5) * t253 + t252 * t567 - t118;
t79 = pkin(4) * t495 - qJ(5) * t253 + t120;
t67 = pkin(10) * t495 + t79;
t19 = -t320 * t67 + t323 * t53;
t20 = t320 * t53 + t323 * t67;
t91 = -mrSges(7,2) * t202 + mrSges(7,3) * t144;
t92 = mrSges(7,1) * t202 - mrSges(7,3) * t145;
t369 = -t320 * t92 + t323 * t91;
t368 = -t320 * t91 - t323 * t92;
t362 = -t255 * pkin(2) + t397;
t119 = -t464 + t136;
t236 = qJD(2) * t360;
t238 = t259 * qJD(2);
t66 = -t226 * t476 - t227 * t477 + t236 * t324 - t321 * t238;
t355 = -pkin(3) * t173 + t362;
t353 = qJ(4) * t176 + t417;
t65 = -t226 * t477 + t227 * t476 + t321 * t236 + t324 * t238;
t350 = -g(1) * t176 - g(2) * t172 - g(3) * t252;
t346 = Ifges(6,5) * t117 + Ifges(6,6) * t116 - Ifges(6,3) * t230;
t342 = qJ(4) * t446 + t65;
t341 = g(1) * t256 + g(2) * t254 - g(3) * t495;
t336 = -qJ(4) * t172 + t355;
t333 = -qJ(5) * t169 - qJD(5) * t253 - t66;
t332 = -qJD(6) * t371 + t1 * t323 - t2 * t320;
t330 = t318 * t364 * (Ifges(3,5) * t325 - Ifges(3,6) * t322);
t329 = qJ(5) * t168 + qJD(5) * t252 + t342;
t292 = Ifges(3,4) * t448;
t279 = -t324 * qJ(5) + t536;
t258 = (-mrSges(3,1) * t325 + mrSges(3,2) * t322) * t318;
t237 = t487 * qJD(1);
t232 = -mrSges(3,2) * t364 + t418;
t181 = Ifges(3,1) * t449 + Ifges(3,5) * t364 + t292;
t164 = t177 * pkin(4);
t162 = t176 * pkin(4);
t160 = t173 * pkin(4);
t158 = t172 * pkin(4);
t152 = -mrSges(6,2) * t275 - t530;
t151 = mrSges(5,1) * t275 + t532;
t150 = -mrSges(4,1) * t275 - mrSges(4,3) * t208;
t149 = mrSges(4,2) * t275 - mrSges(4,3) * t207;
t133 = mrSges(5,1) * t207 - mrSges(5,3) * t208;
t131 = mrSges(6,1) * t208 + mrSges(6,2) * t207;
t125 = t176 * t323 - t256 * t320;
t124 = -t176 * t320 - t256 * t323;
t123 = -pkin(3) * t449 - t138;
t106 = t116 * mrSges(5,2);
t105 = t116 * mrSges(6,3);
t89 = -pkin(4) * t208 - t132;
t86 = -t119 - t501;
t84 = -t118 - t243;
t81 = qJD(6) * t170 + t168 * t323 - t320 * t446;
t80 = -qJD(6) * t171 - t168 * t320 - t323 * t446;
t78 = -t319 * t495 + t136 + t501;
t77 = -mrSges(5,2) * t117 + mrSges(5,3) * t230;
t76 = -mrSges(4,2) * t230 - mrSges(4,3) * t117;
t75 = -mrSges(6,1) * t230 - mrSges(6,3) * t117;
t74 = -t230 * mrSges(5,1) - t106;
t73 = mrSges(4,1) * t230 + mrSges(4,3) * t116;
t72 = t230 * mrSges(6,2) + t105;
t69 = -t194 + t489;
t59 = -pkin(3) * t446 - t66;
t52 = t342 - t444;
t50 = mrSges(4,1) * t117 - mrSges(4,2) * t116;
t49 = mrSges(5,1) * t117 + mrSges(5,3) * t116;
t48 = -pkin(5) * t207 + t208 * t567 - t132;
t42 = -t168 * pkin(4) - t63;
t34 = -t329 + t444;
t33 = -t446 * t568 + t333;
t28 = (pkin(5) * t482 - t475) * t318 + t329;
t27 = -qJD(2) * t318 * t421 + t333;
t22 = t169 * pkin(5) + t168 * t567 - t63;
t18 = t320 * t48 + t323 * t70;
t17 = -t320 * t70 + t323 * t48;
t16 = -pkin(4) * t117 + t335;
t15 = -mrSges(7,1) * t46 + mrSges(7,2) * t45;
t10 = pkin(5) * t230 - t12;
t4 = -qJD(6) * t20 + t22 * t323 - t27 * t320;
t3 = qJD(6) * t19 + t22 * t320 + t27 * t323;
t29 = [(t26 * mrSges(5,2) - t32 * mrSges(4,3) + Ifges(4,4) * t563 - t562 * t654 + t646) * t253 + (Ifges(4,4) * t553 + Ifges(7,5) * t558 - Ifges(6,2) * t551 - Ifges(6,6) * t544 + Ifges(7,6) * t560 + Ifges(7,3) * t554 + t622 * t545 + t624 * t550 - t552 * t654 + t634 + t649) * t169 + (t318 * t181 + t576 * qJD(1) * (-Ifges(3,2) * t322 + t526)) * t481 / 0.2e1 + (-t146 * t502 - t339 * t539 - t361 * t487) * mrSges(3,2) + (Ifges(7,5) * t171 + Ifges(7,6) * t170) * t566 + (Ifges(7,5) * t81 + Ifges(7,6) * t80) * t554 + (t147 * t502 + t241 * t539 + t259 * t361) * mrSges(3,1) + (-m(4) * (t362 - t535) - m(5) * (t336 - t535) - m(7) * (-t160 + t355) - t500 * mrSges(7,1) - t499 * mrSges(7,2) - m(6) * (-t160 + t336) - m(3) * t397 + t255 * mrSges(3,1) - mrSges(3,3) * t452 + mrSges(2,1) * t540 + mrSges(2,2) * t541 - (-t388 + t580) * t172 - t636 * t254 + t579 * t173) * g(1) + (-m(3) * t234 - t584) * t239 + (-t21 * mrSges(5,2) - t31 * mrSges(4,3) + Ifges(6,4) * t564 + t623 * t565 + t647) * t252 + m(7) * (t1 * t20 + t10 * t78 + t13 * t4 + t14 * t3 + t19 * t2 + t28 * t58) + m(6) * (t11 * t79 + t12 * t86 + t16 * t84 + t33 * t57 + t34 * t64 + t42 * t62) + m(5) * (t118 * t23 + t119 * t21 + t120 * t26 + t52 * t90 + t59 * t88 + t63 * t83) + (Ifges(7,4) * t171 + Ifges(7,2) * t170) * t570 + (Ifges(7,4) * t81 + Ifges(7,2) * t80) * t560 + (Ifges(3,1) * t339 + Ifges(3,4) * t241 + Ifges(3,5) * t361) * t443 - pkin(1) * t258 * t469 + t16 * t426 + t80 * t657 + (Ifges(7,1) * t171 + Ifges(7,4) * t170) * t571 + (Ifges(7,1) * t81 + Ifges(7,4) * t80) * t558 + (t613 * t443 + t330 / 0.2e1) * qJD(2) + (-m(6) * (t164 + t353) - m(3) * t485 - t257 * mrSges(3,1) - mrSges(3,3) * t451 - m(5) * (t353 + t538) - m(4) * (t454 + t538) - m(7) * (t164 + t417) - t125 * mrSges(7,1) - t124 * mrSges(7,2) - mrSges(2,1) * t541 + mrSges(2,2) * t540 + t580 * t176 + t636 * t256 - t579 * t177) * g(2) - t10 * t389 + (Ifges(3,4) * t643 + Ifges(3,6) * t642 - Ifges(4,6) * t563 + Ifges(6,6) * t564 + t620 * t562 - t619 * t548 - t622 * t565 + t146 * mrSges(3,3) + Ifges(6,3) * t549 + t346 / 0.2e1 - t598 / 0.2e1 + t631 + Ifges(3,2) * t644) * t495 + (t1 * t170 - t13 * t81 + t14 * t80 - t171 * t2) * mrSges(7,3) + (Ifges(3,5) * t502 + (t322 * Ifges(3,1) + t526) * t318) * t643 + t334 * t644 + (Ifges(3,3) * t502 + (Ifges(3,5) * t322 + Ifges(3,6) * t325) * t318) * t642 + t170 * t574 + (-t114 * mrSges(4,2) - t64 * mrSges(6,1) - t88 * mrSges(5,1) + t90 * mrSges(5,3) + t57 * mrSges(6,2) - t489 * mrSges(4,1) - t611 / 0.2e1 - t620 * t552 + t619 * t545 + t622 * t550 - t237 * mrSges(3,3) - Ifges(6,3) * t544 - Ifges(6,6) * t551 + Ifges(4,6) * t553) * t446 + m(4) * (t114 * t65 + t130 * t225 + t135 * t32 + t136 * t31 - t489 * t66) + t502 * t455 / 0.2e1 + Ifges(2,3) * qJDD(1) + t238 * t232 + t225 * t50 + t171 * t8 / 0.2e1 + t52 * t153 + t34 * t148 + t65 * t149 + t66 * t150 + t59 * t151 + t33 * t152 + t63 * t133 + t135 * t73 + t136 * t76 + t42 * t131 + t118 * t49 + t119 * t77 + t120 * t74 + t3 * t91 + t4 * t92 + t58 * (-mrSges(7,1) * t80 + mrSges(7,2) * t81) + t81 * t56 / 0.2e1 + t84 * t51 + t86 * t75 + t78 * t15 + t79 * t72 + t28 * t68 + t19 * t24 + t20 * t25 + t593 * t470 + (Ifges(6,4) * t551 + Ifges(6,5) * t544 - Ifges(4,2) * t553 - t618 * t545 + t623 * t550 - t617 * t552 + t635 + t645) * t168 + (-t147 * t496 - t234 * t445 + t241 * t487 - t259 * t339) * mrSges(3,3) + m(3) * (pkin(1) ^ 2 * qJDD(1) * t576 + t146 * t487 + t147 * t259 + t237 * t238) + t23 * t392 + t130 * t394; (t77 + t76) * t536 + (-t10 * t387 - t373 * t566 - t378 * t570 - t383 * t571 + t98 * t402 - t647) * t324 - t522 * t564 + ((-Ifges(3,2) * t449 + t324 * t590 + t181 + t292) * t325 + t613 * t322 + (t322 * t622 + t325 * t599) * t208 + (-t322 * t620 + t325 * t600) * t207 + t275 * (-Ifges(6,3) * t322 + t325 * t374)) * t438 + (t320 * t56 + t323 * t55) * t472 / 0.2e1 + (-t330 / 0.2e1 - t593 * qJD(1)) * qJD(1) + (-t468 - t139) * t149 + (t419 + t584) * t237 + (-t467 - t138) * t150 + (Ifges(7,4) * t190 + Ifges(7,2) * t189) * t561 + (Ifges(7,1) * t190 + Ifges(7,4) * t189) * t559 + (-pkin(2) * t130 + pkin(9) * t596 - t114 * t139 + t489 * t138) * m(4) + (pkin(9) * t595 - t121 * t90 - t123 * t88 + t23 * t269 + t602 * t83) * m(5) + (-t232 + t418) * t234 + t595 * mrSges(5,2) - t578 * t413 + t650 * t153 + t455 + ((t374 / 0.2e1 - t601 / 0.2e1) * t275 + (-t379 / 0.2e1 + t599 / 0.2e1) * t208 + (-t381 / 0.2e1 + t600 / 0.2e1) * t207 + (Ifges(7,3) * t324 + t321 * t373) * t554 + (Ifges(7,5) * t324 + t321 * t383) * t558 + (Ifges(7,6) * t324 + t321 * t378) * t560 - t589 + ((-t114 * t321 + t324 * t489) * m(4) + (-t321 * t90 + t324 * t88) * m(5)) * pkin(9)) * qJD(3) + (t15 - t75) * t279 + (-m(4) * t488 - m(5) * t414 + t258 + t629 * (pkin(4) * t461 + t414) + (t637 * t325 + (-qJ(5) * t629 + t387 - t591) * t322) * t318) * g(3) + (-t516 + t524) * t565 + (t611 * t322 + (t322 * t619 + t325 * t601) * t275 + t208 * (-Ifges(6,6) * t322 + t325 * t379) + t207 * (Ifges(4,6) * t322 + t325 * t381)) * t437 + t492 * t573 + t494 * t574 + (t100 * t402 + (-t73 + t74) * pkin(9) + t612 * t325 * t438 + t646) * t321 + (t372 * t554 + t377 * t560 + t382 * t558) * t472 + t525 * t563 + (Ifges(7,5) * t190 + Ifges(7,6) * t189) * t555 + (t467 - t123) * t151 + (-t57 * (mrSges(6,2) * t322 - mrSges(6,3) * t490) + t489 * (mrSges(4,1) * t322 - mrSges(4,3) * t490) - t88 * (-mrSges(5,1) * t322 + mrSges(5,2) * t490) - t64 * (-mrSges(6,1) * t322 - mrSges(6,3) * t493) - t114 * (-mrSges(4,2) * t322 - mrSges(4,3) * t493) - t90 * (-mrSges(5,2) * t493 + mrSges(5,3) * t322)) * t484 + (t645 - t648) * t477 + t649 * t476 + (t517 - t523) * t562 + t278 * t72 + t261 * t51 + t269 * t49 - t189 * t55 / 0.2e1 - t190 * t56 / 0.2e1 - t146 * mrSges(3,2) + t147 * mrSges(3,1) + t142 * t24 + t143 * t25 - pkin(2) * t50 + (-m(5) * (t427 + t603) - m(4) * t427 + t629 * (-pkin(4) * t497 - t247 + t603) + t583 * t257 + t581 * t256) * g(1) + (-m(5) * (t428 + t604) - m(4) * t428 + t629 * (-pkin(4) * t498 - t245 + t604) + t583 * t255 + t581 * t254) * g(2) + t589 * t448 + t596 * mrSges(4,3) + t602 * t133 + t605 * t131 + (mrSges(7,1) * t607 - mrSges(7,2) * t606) * t58 + (t1 * t494 + t13 * t606 - t14 * t607 + t2 * t492) * mrSges(7,3) + t608 * t152 + t609 * t148 + (t11 * t278 - t12 * t279 + t16 * t261 + t57 * t608 + t605 * t62 + t609 * t64) * m(6) + t610 * t68 + t614 * t92 + t615 * t91 + (t1 * t143 + t10 * t279 + t13 * t614 + t14 * t615 + t142 * t2 + t58 * t610) * m(7) + t16 * t390 - t23 * t391 - t130 * t393; (t77 - t75) * qJ(4) + (t150 - t151) * t114 + (-t518 + t197 + t98) * t553 + (t10 * t319 - t13 * t17 - t14 * t18 + t58 * t69) * m(7) + t577 * qJD(6) - t1 * t528 + t7 * t656 - t346 + t10 * t388 + t598 + (-m(6) * (-t158 + t424) - m(5) * t424 - m(7) * (-t158 - t159) + t587 * t173 + t579 * t172) * g(2) + (-m(6) * (-t162 + t423) - m(5) * t423 - m(7) * (-t162 - t163) + t587 * t177 + t579 * t176) * g(1) + (-m(5) * t422 + t392 - m(6) * (-t243 + t422) - t426 - m(7) * (-pkin(10) * t252 - t243 - t244) + t252 * mrSges(7,3) - t337 * t253 + t394) * g(3) + t320 * t573 - t112 * t372 / 0.2e1 + (t13 * t473 + t14 * t474) * mrSges(7,3) + (m(5) * t90 - m(6) * t64 + m(7) * t58 + t639) * qJD(4) - t631 - (-t153 - t149) * t489 + (-pkin(3) * t26 + qJ(4) * t21 - t114 * t88 - t132 * t83 + t489 * t90) * m(5) - t568 * t72 + (-t12 * qJ(4) - t11 * t568 - t57 * t70 - t62 * t89 - t64 * t69) * m(6) - t57 * t531 + t319 * t15 + (t512 + t100) * t550 - t70 * t152 - t132 * t133 - t89 * t131 - t18 * t91 - t17 * t92 - pkin(3) * t74 - t46 * t377 / 0.2e1 - t45 * t382 / 0.2e1 + (-t198 + t590) * t552 - (m(7) * t332 - t473 * t92 - t474 * t91 + t597) * t471 - t503 * t69 + (t196 - t513 + t612) * t551 + (Ifges(6,2) * t550 + Ifges(6,6) * t545 - t544 * t622 - t551 * t624 + t515 + t578 + t634) * t207 + (Ifges(6,5) * t545 - Ifges(4,2) * t552 + t13 * t528 + t14 * t529 - t544 * t618 - t553 * t617 + t514 + t577 - t635) * t208 + t2 * t529 + t64 * t530 + t90 * t532 + t88 * t533; t105 - t106 + t626 * t230 + t368 * qJD(6) + t639 * t275 + (-t131 + t133 + t368) * t208 + (-t208 * t371 + t275 * t58 + t332 + t350) * m(7) + (-t208 * t62 - t275 * t64 + t11 + t350) * m(6) + (t208 * t83 + t275 * t90 + t26 + t350) * m(5) + t597; t323 * t24 + t320 * t25 + t503 * t207 + t369 * qJD(6) + (t152 + t369) * t208 + t51 + (t1 * t320 + t2 * t323 - t207 * t58 + t341 - t202 * (t13 * t320 - t14 * t323)) * m(7) + (t207 * t64 + t208 * t57 + t16 + t341) * m(6); -t58 * (mrSges(7,1) * t145 + mrSges(7,2) * t144) + (Ifges(7,1) * t144 - t521) * t559 + t55 * t558 + (Ifges(7,5) * t144 - Ifges(7,6) * t145) * t555 - t13 * t91 + t14 * t92 - g(1) * (mrSges(7,1) * t124 - mrSges(7,2) * t125) - g(2) * ((-t172 * t320 - t499) * mrSges(7,1) + (-t172 * t323 + t500) * mrSges(7,2)) - g(3) * t389 + (t13 * t144 + t14 * t145) * mrSges(7,3) + t6 + (-Ifges(7,2) * t145 + t141 + t56) * t561 + t585;];
tau  = t29;
