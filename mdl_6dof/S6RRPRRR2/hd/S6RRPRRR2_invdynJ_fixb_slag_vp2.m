% Calculate vector of inverse dynamics joint torques for
% S6RRPRRR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:16:42
% EndTime: 2019-03-09 13:17:30
% DurationCPUTime: 27.57s
% Computational Cost: add. (26494->861), mult. (61421->1131), div. (0->0), fcn. (47812->18), ass. (0->418)
t377 = sin(qJ(5));
t460 = qJD(5) * t377;
t373 = sin(pkin(11));
t374 = cos(pkin(11));
t379 = sin(qJ(2));
t384 = cos(qJ(2));
t306 = -t373 * t379 + t374 * t384;
t288 = t306 * qJD(1);
t466 = qJD(1) * t384;
t467 = qJD(1) * t379;
t289 = -t373 * t466 - t374 * t467;
t378 = sin(qJ(4));
t383 = cos(qJ(4));
t600 = t383 * t288 + t378 * t289;
t638 = t600 * t377;
t656 = t460 - t638;
t372 = qJ(5) + qJ(6);
t364 = sin(t372);
t365 = cos(t372);
t382 = cos(qJ(5));
t632 = -mrSges(6,1) * t382 + mrSges(6,2) * t377;
t655 = -mrSges(7,1) * t365 + mrSges(7,2) * t364 + t632;
t649 = -mrSges(6,3) - mrSges(7,3);
t545 = pkin(2) * t374;
t353 = pkin(3) + t545;
t546 = pkin(2) * t373;
t283 = t353 * t383 - t378 * t546;
t267 = t283 * qJD(4);
t375 = -qJ(3) - pkin(7);
t335 = t375 * t379;
t318 = qJD(1) * t335;
t337 = t375 * t384;
t319 = qJD(1) * t337;
t478 = t374 * t319;
t234 = -t318 * t373 + t478;
t538 = pkin(8) * t288;
t202 = t234 - t538;
t292 = t373 * t319;
t235 = t374 * t318 + t292;
t537 = pkin(8) * t289;
t203 = t235 + t537;
t150 = t202 * t378 + t203 * t383;
t409 = t288 * t378 - t383 * t289;
t171 = pkin(4) * t409 - pkin(9) * t600;
t358 = pkin(2) * t467;
t252 = -pkin(3) * t289 + t358;
t151 = t171 + t252;
t81 = t382 * t150 + t377 * t151;
t654 = t267 * t382 - t81;
t80 = -t150 * t377 + t382 * t151;
t653 = -t267 * t377 - t80;
t370 = qJD(2) + qJD(4);
t200 = t370 * t382 - t377 * t409;
t213 = qJD(5) - t600;
t201 = t370 * t377 + t382 * t409;
t510 = t201 * Ifges(6,4);
t110 = t200 * Ifges(6,2) + t213 * Ifges(6,6) + t510;
t651 = -t110 / 0.2e1;
t650 = pkin(10) * t638;
t284 = t378 * t353 + t383 * t546;
t281 = pkin(9) + t284;
t528 = -pkin(10) - t281;
t433 = qJD(5) * t528;
t648 = t377 * t433 + t650 + t654;
t637 = t600 * t382;
t642 = pkin(5) * t409 - pkin(10) * t637;
t647 = t382 * t433 - t642 + t653;
t386 = -pkin(10) - pkin(9);
t446 = qJD(5) * t386;
t299 = qJD(2) * pkin(2) + t318;
t227 = t374 * t299 + t292;
t191 = qJD(2) * pkin(3) + t227 + t537;
t228 = t373 * t299 - t478;
t195 = t228 + t538;
t123 = t191 * t383 - t378 * t195;
t86 = t382 * t123 + t377 * t171;
t646 = t377 * t446 + t650 - t86;
t85 = -t123 * t377 + t382 * t171;
t645 = t382 * t446 - t642 - t85;
t371 = qJ(2) + pkin(11);
t363 = qJ(4) + t371;
t351 = sin(t363);
t352 = cos(t363);
t540 = pkin(5) * t382;
t354 = pkin(4) + t540;
t644 = m(7) * t352 * t386 + (m(7) * t354 - t655) * t351;
t643 = t656 * pkin(5);
t116 = -pkin(4) * t370 - t123;
t102 = -pkin(5) * t200 + t116;
t458 = qJD(1) * qJD(2);
t438 = t379 * t458;
t457 = qJDD(1) * t384;
t322 = -t438 + t457;
t323 = qJDD(1) * t379 + t384 * t458;
t236 = t322 * t374 - t323 * t373;
t237 = t322 * t373 + t323 * t374;
t136 = qJD(4) * t600 + t236 * t378 + t237 * t383;
t461 = qJD(4) * t383;
t462 = qJD(4) * t378;
t137 = t236 * t383 - t378 * t237 - t288 * t462 + t289 * t461;
t381 = cos(qJ(6));
t376 = sin(qJ(6));
t124 = t191 * t378 + t195 * t383;
t117 = pkin(9) * t370 + t124;
t367 = t384 * pkin(2);
t355 = t367 + pkin(1);
t325 = -qJD(1) * t355 + qJD(3);
t240 = -pkin(3) * t288 + t325;
t140 = -pkin(4) * t600 - pkin(9) * t409 + t240;
t76 = t117 * t382 + t140 * t377;
t65 = pkin(10) * t200 + t76;
t505 = t376 * t65;
t75 = -t117 * t377 + t382 * t140;
t64 = -pkin(10) * t201 + t75;
t55 = pkin(5) * t213 + t64;
t20 = t381 * t55 - t505;
t503 = t381 * t65;
t21 = t376 * t55 + t503;
t406 = t376 * t377 - t381 * t382;
t593 = qJD(5) + qJD(6);
t238 = t593 * t406;
t314 = t376 * t382 + t377 * t381;
t239 = t593 * t314;
t459 = qJD(5) * t382;
t368 = qJDD(2) + qJDD(4);
t305 = t323 * pkin(7);
t464 = qJD(3) * t379;
t225 = qJDD(2) * pkin(2) - qJ(3) * t323 - qJD(1) * t464 - t305;
t356 = pkin(7) * t457;
t465 = qJD(2) * t379;
t449 = pkin(7) * t465;
t463 = qJD(3) * t384;
t232 = qJ(3) * t322 + t356 + (-t449 + t463) * qJD(1);
t174 = t374 * t225 - t232 * t373;
t130 = qJDD(2) * pkin(3) - pkin(8) * t237 + t174;
t175 = t373 * t225 + t374 * t232;
t139 = pkin(8) * t236 + t175;
t52 = t378 * t130 + t383 * t139 + t191 * t461 - t195 * t462;
t48 = pkin(9) * t368 + t52;
t501 = qJDD(1) * pkin(1);
t275 = -pkin(2) * t322 + qJDD(3) - t501;
t199 = -pkin(3) * t236 + t275;
t68 = -pkin(4) * t137 - pkin(9) * t136 + t199;
t15 = -t117 * t460 + t140 * t459 + t377 * t68 + t382 * t48;
t97 = -qJD(5) * t201 - t136 * t377 + t368 * t382;
t10 = pkin(10) * t97 + t15;
t135 = qJDD(5) - t137;
t16 = -qJD(5) * t76 - t377 * t48 + t382 * t68;
t96 = qJD(5) * t200 + t136 * t382 + t368 * t377;
t7 = pkin(5) * t135 - pkin(10) * t96 + t16;
t3 = qJD(6) * t20 + t10 * t381 + t376 * t7;
t53 = t130 * t383 - t378 * t139 - t191 * t462 - t195 * t461;
t49 = -pkin(4) * t368 - t53;
t33 = -pkin(5) * t97 + t49;
t380 = sin(qJ(1));
t385 = cos(qJ(1));
t418 = mrSges(6,1) * t377 + mrSges(6,2) * t382;
t399 = t116 * t418;
t4 = -qJD(6) * t21 - t10 * t376 + t381 * t7;
t40 = t96 * Ifges(6,4) + t97 * Ifges(6,2) + t135 * Ifges(6,6);
t412 = Ifges(6,5) * t382 - Ifges(6,6) * t377;
t517 = Ifges(6,4) * t382;
t414 = -Ifges(6,2) * t377 + t517;
t518 = Ifges(6,4) * t377;
t416 = Ifges(6,1) * t382 - t518;
t436 = -t460 / 0.2e1;
t198 = Ifges(6,4) * t200;
t111 = t201 * Ifges(6,1) + t213 * Ifges(6,5) + t198;
t473 = t382 * t111;
t440 = t473 / 0.2e1;
t484 = t352 * t385;
t485 = t352 * t380;
t207 = qJD(6) + t213;
t556 = t207 / 0.2e1;
t557 = -t207 / 0.2e1;
t148 = t200 * t376 + t201 * t381;
t562 = t148 / 0.2e1;
t563 = -t148 / 0.2e1;
t428 = t381 * t200 - t201 * t376;
t564 = t428 / 0.2e1;
t565 = -t428 / 0.2e1;
t566 = t135 / 0.2e1;
t127 = qJDD(6) + t135;
t567 = t127 / 0.2e1;
t568 = t97 / 0.2e1;
t569 = t96 / 0.2e1;
t141 = Ifges(7,4) * t428;
t72 = Ifges(7,1) * t148 + Ifges(7,5) * t207 + t141;
t570 = t72 / 0.2e1;
t571 = -t72 / 0.2e1;
t516 = Ifges(7,4) * t148;
t71 = Ifges(7,2) * t428 + Ifges(7,6) * t207 + t516;
t572 = t71 / 0.2e1;
t573 = -t71 / 0.2e1;
t45 = -qJD(6) * t148 - t376 * t96 + t381 * t97;
t574 = t45 / 0.2e1;
t44 = qJD(6) * t428 + t376 * t97 + t381 * t96;
t575 = t44 / 0.2e1;
t576 = Ifges(6,1) * t569 + Ifges(6,4) * t568 + Ifges(6,5) * t566;
t577 = Ifges(7,1) * t575 + Ifges(7,4) * t574 + Ifges(7,5) * t567;
t578 = Ifges(7,4) * t575 + Ifges(7,2) * t574 + Ifges(7,6) * t567;
t627 = mrSges(7,2) * t102 - t20 * mrSges(7,3);
t628 = -mrSges(7,1) * t102 + t21 * mrSges(7,3);
t633 = t15 * t382 - t16 * t377;
t639 = t406 * t600;
t640 = t314 * t600;
t641 = -(Ifges(7,1) * t562 + Ifges(7,4) * t564 + Ifges(7,5) * t556 + t570 + t627) * t238 - (Ifges(7,4) * t562 + Ifges(7,2) * t564 + Ifges(7,6) * t556 + t572 + t628) * t239 + (t200 * t414 + t201 * t416 + t213 * t412) * qJD(5) / 0.2e1 + t49 * t632 + (t633 - t656 * t76 + (-t459 + t637) * t75) * mrSges(6,3) + t33 * (mrSges(7,1) * t406 + mrSges(7,2) * t314) + (Ifges(7,5) * t314 - Ifges(7,6) * t406) * t567 + (Ifges(7,4) * t314 - Ifges(7,2) * t406) * t574 + (Ifges(7,1) * t314 - Ifges(7,4) * t406) * t575 - t406 * t578 - t640 * t573 + (-Ifges(7,4) * t639 - Ifges(7,2) * t640) * t565 - t102 * (mrSges(7,1) * t640 - mrSges(7,2) * t639) + (-Ifges(7,5) * t639 - Ifges(7,6) * t640) * t557 + (-Ifges(7,1) * t639 - Ifges(7,4) * t640) * t563 - t639 * t571 + (-t20 * t639 + t21 * t640 - t3 * t406 - t314 * t4) * mrSges(7,3) + (t399 + t440) * qJD(5) + t382 * t40 / 0.2e1 + Ifges(5,3) * t368 + (Ifges(6,5) * t377 + Ifges(6,6) * t382) * t566 - t52 * mrSges(5,2) + t53 * mrSges(5,1) + (t385 * t644 + t484 * t649) * g(1) + (t380 * t644 + t485 * t649) * g(2) + t110 * t436 + (Ifges(6,2) * t382 + t518) * t568 + (Ifges(6,1) * t377 + t517) * t569 + t377 * t576 + t314 * t577 + Ifges(5,5) * t136 + Ifges(5,6) * t137;
t579 = m(7) * pkin(5);
t601 = t284 * qJD(4) + t383 * t202 - t203 * t378;
t336 = -t384 * mrSges(3,1) + t379 * mrSges(3,2);
t361 = sin(t371);
t362 = cos(t371);
t636 = -t362 * mrSges(4,1) + t361 * mrSges(4,2) + t336;
t631 = g(1) * t385 + g(2) * t380;
t630 = t240 * mrSges(5,2) - t123 * mrSges(5,3) + t377 * t651;
t629 = -t352 * mrSges(5,1) + (mrSges(5,2) + t649) * t351;
t597 = t352 * pkin(4) + t351 * pkin(9);
t599 = -t351 * t386 + t352 * t354;
t626 = -m(6) * t597 - m(7) * t599;
t625 = t322 / 0.2e1;
t624 = t323 / 0.2e1;
t552 = -t409 / 0.2e1;
t623 = t377 * t579;
t530 = qJD(2) / 0.2e1;
t243 = t528 * t377;
t366 = t382 * pkin(10);
t489 = t281 * t382;
t244 = t366 + t489;
t184 = t243 * t376 + t244 * t381;
t620 = -qJD(6) * t184 - t376 * t648 + t381 * t647;
t183 = t243 * t381 - t244 * t376;
t619 = qJD(6) * t183 + t376 * t647 + t381 * t648;
t212 = Ifges(5,4) * t600;
t307 = t373 * t384 + t374 * t379;
t618 = Ifges(4,5) * t307;
t617 = Ifges(4,6) * t306;
t616 = t200 * Ifges(6,6);
t615 = t213 * Ifges(6,3);
t614 = t600 * Ifges(5,2);
t613 = t370 * Ifges(5,5);
t612 = t370 * Ifges(5,6);
t502 = qJDD(2) / 0.2e1;
t611 = mrSges(6,1) + t579;
t610 = t201 * Ifges(6,5) + t148 * Ifges(7,5) + Ifges(7,6) * t428 + t207 * Ifges(7,3) + t615 + t616;
t339 = t386 * t377;
t536 = pkin(9) * t382;
t340 = t366 + t536;
t246 = t339 * t381 - t340 * t376;
t609 = qJD(6) * t246 + t376 * t645 + t381 * t646;
t247 = t339 * t376 + t340 * t381;
t608 = -qJD(6) * t247 - t376 * t646 + t381 * t645;
t231 = t306 * t378 + t307 * t383;
t173 = t406 * t231;
t241 = t374 * t335 + t337 * t373;
t210 = -pkin(8) * t307 + t241;
t242 = t373 * t335 - t374 * t337;
t211 = pkin(8) * t306 + t242;
t165 = t210 * t378 + t211 * t383;
t158 = t382 * t165;
t263 = -pkin(3) * t306 - t355;
t408 = t383 * t306 - t307 * t378;
t159 = -pkin(4) * t408 - pkin(9) * t231 + t263;
t88 = t377 * t159 + t158;
t607 = mrSges(5,1) * t370 + mrSges(6,1) * t200 - mrSges(6,2) * t201 - mrSges(5,3) * t409;
t606 = t383 * t210 - t211 * t378;
t605 = -t238 + t639;
t604 = -t239 + t640;
t603 = t601 + t643;
t602 = t267 - t150;
t596 = -t124 + t643;
t304 = -pkin(7) * t438 + t356;
t595 = t304 * t384 + t305 * t379;
t594 = m(7) + m(6) + m(5);
t591 = t352 * t655 + t629;
t589 = -m(3) * pkin(7) + m(4) * t375 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t588 = t4 * mrSges(7,1) - t3 * mrSges(7,2);
t587 = t16 * mrSges(6,1) - t15 * mrSges(6,2);
t523 = mrSges(6,3) * t200;
t154 = -mrSges(6,2) * t213 + t523;
t522 = mrSges(6,3) * t201;
t155 = mrSges(6,1) * t213 - t522;
t61 = mrSges(6,1) * t135 - mrSges(6,3) * t96;
t586 = m(6) * ((-t377 * t76 - t382 * t75) * qJD(5) + t633) - t155 * t459 - t154 * t460 - t377 * t61;
t585 = m(3) * pkin(1) + m(4) * t355 + mrSges(2,1) - t629 - t636;
t583 = t240 * mrSges(5,1) + t75 * mrSges(6,1) + t20 * mrSges(7,1) - t76 * mrSges(6,2) - t21 * mrSges(7,2) - t124 * mrSges(5,3);
t547 = -t370 / 0.2e1;
t555 = -t213 / 0.2e1;
t559 = -t201 / 0.2e1;
t560 = -t200 / 0.2e1;
t582 = Ifges(5,1) * t552 + Ifges(5,5) * t547 + t412 * t555 + t414 * t560 + t416 * t559 - t399 - t473 / 0.2e1 - t630;
t553 = -t600 / 0.2e1;
t581 = Ifges(6,5) * t559 + Ifges(7,5) * t563 - Ifges(5,2) * t553 - Ifges(5,6) * t547 + Ifges(6,6) * t560 + Ifges(7,6) * t565 + Ifges(6,3) * t555 + Ifges(7,3) * t557 - t583;
t558 = t201 / 0.2e1;
t551 = t409 / 0.2e1;
t548 = -t289 / 0.2e1;
t544 = pkin(2) * t379;
t543 = pkin(4) * t351;
t542 = pkin(5) * t201;
t539 = pkin(7) * t384;
t533 = g(3) * t351;
t521 = Ifges(3,4) * t379;
t520 = Ifges(3,4) * t384;
t519 = Ifges(5,4) * t409;
t515 = pkin(5) * qJD(6);
t509 = t227 * mrSges(4,3);
t508 = t228 * mrSges(4,3);
t507 = t289 * Ifges(4,4);
t290 = t307 * qJD(2);
t291 = t306 * qJD(2);
t176 = qJD(4) * t408 - t290 * t378 + t291 * t383;
t500 = t176 * t382;
t493 = t231 * t377;
t492 = t231 * t382;
t482 = t364 * t380;
t481 = t364 * t385;
t480 = t365 * t380;
t479 = t365 * t385;
t476 = t377 * t380;
t475 = t377 * t385;
t474 = t380 * t382;
t472 = t382 * t385;
t254 = t352 * t482 + t479;
t255 = -t352 * t480 + t481;
t471 = -t254 * mrSges(7,1) + t255 * mrSges(7,2);
t256 = -t352 * t481 + t480;
t257 = t352 * t479 + t482;
t470 = t256 * mrSges(7,1) - t257 * mrSges(7,2);
t434 = qJD(2) * t375;
t285 = t379 * t434 + t463;
t286 = t384 * t434 - t464;
t209 = t374 * t285 + t373 * t286;
t468 = pkin(3) * t362 + t367;
t456 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t467) * t539;
t455 = Ifges(7,5) * t44 + Ifges(7,6) * t45 + Ifges(7,3) * t127;
t454 = Ifges(6,5) * t96 + Ifges(6,6) * t97 + Ifges(6,3) * t135;
t359 = pkin(2) * t465;
t443 = t231 * t459;
t253 = pkin(3) * t290 + t359;
t208 = -t285 * t373 + t374 * t286;
t185 = -pkin(8) * t291 + t208;
t186 = -pkin(8) * t290 + t209;
t83 = qJD(4) * t606 + t185 * t378 + t186 * t383;
t177 = qJD(4) * t231 + t383 * t290 + t291 * t378;
t93 = pkin(4) * t177 - pkin(9) * t176 + t253;
t435 = -t377 * t83 + t382 * t93;
t431 = -t236 * mrSges(4,1) + t237 * mrSges(4,2);
t430 = -t137 * mrSges(5,1) + t136 * mrSges(5,2);
t87 = t382 * t159 - t165 * t377;
t425 = pkin(9) * t485 - t380 * t543;
t424 = pkin(9) * t484 - t385 * t543;
t280 = -pkin(4) - t283;
t422 = mrSges(3,1) * t379 + mrSges(3,2) * t384;
t419 = mrSges(5,1) * t351 + mrSges(5,2) * t352;
t417 = -mrSges(7,1) * t364 - mrSges(7,2) * t365;
t415 = t384 * Ifges(3,2) + t521;
t413 = Ifges(3,5) * t384 - Ifges(3,6) * t379;
t73 = -pkin(5) * t408 - pkin(10) * t492 + t87;
t78 = -pkin(10) * t493 + t88;
t36 = -t376 * t78 + t381 * t73;
t37 = t376 * t73 + t381 * t78;
t411 = -t377 * t75 + t382 * t76;
t410 = t154 * t382 - t155 * t377;
t403 = t455 + t588;
t402 = pkin(1) * t422;
t278 = -t352 * t475 + t474;
t276 = t352 * t476 + t472;
t401 = t176 * t377 + t443;
t400 = t231 * t460 - t500;
t398 = t379 * (Ifges(3,1) * t384 - t521);
t22 = t159 * t459 - t165 * t460 + t377 * t93 + t382 * t83;
t84 = qJD(4) * t165 - t383 * t185 + t186 * t378;
t369 = -pkin(8) + t375;
t357 = Ifges(3,4) * t466;
t334 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t466;
t324 = -pkin(3) * t361 - t544;
t317 = pkin(1) + t468;
t303 = t385 * t324;
t302 = t380 * t324;
t298 = Ifges(3,1) * t467 + Ifges(3,5) * qJD(2) + t357;
t297 = Ifges(3,6) * qJD(2) + qJD(1) * t415;
t282 = Ifges(4,4) * t288;
t279 = t352 * t472 + t476;
t277 = -t352 * t474 + t475;
t260 = qJD(2) * mrSges(4,1) + mrSges(4,3) * t289;
t259 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t288;
t258 = t280 - t540;
t226 = -mrSges(4,1) * t288 - mrSges(4,2) * t289;
t223 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t237;
t222 = -qJDD(2) * mrSges(4,2) + mrSges(4,3) * t236;
t217 = -t289 * Ifges(4,1) + Ifges(4,5) * qJD(2) + t282;
t216 = t288 * Ifges(4,2) + Ifges(4,6) * qJD(2) - t507;
t204 = -mrSges(5,2) * t370 + mrSges(5,3) * t600;
t172 = t314 * t231;
t170 = -mrSges(5,1) * t600 + mrSges(5,2) * t409;
t167 = Ifges(5,1) * t409 + t212 + t613;
t166 = t519 + t612 + t614;
t119 = -mrSges(5,2) * t368 + mrSges(5,3) * t137;
t118 = mrSges(5,1) * t368 - mrSges(5,3) * t136;
t114 = pkin(5) * t493 - t606;
t105 = mrSges(7,1) * t207 - mrSges(7,3) * t148;
t104 = -mrSges(7,2) * t207 + mrSges(7,3) * t428;
t82 = -mrSges(7,1) * t428 + mrSges(7,2) * t148;
t62 = -mrSges(6,2) * t135 + mrSges(6,3) * t97;
t57 = t173 * t593 - t314 * t176;
t56 = -t176 * t406 - t231 * t239;
t54 = pkin(5) * t401 + t84;
t50 = -mrSges(6,1) * t97 + mrSges(6,2) * t96;
t27 = -mrSges(7,2) * t127 + mrSges(7,3) * t45;
t26 = mrSges(7,1) * t127 - mrSges(7,3) * t44;
t25 = t381 * t64 - t505;
t24 = -t376 * t64 - t503;
t23 = -qJD(5) * t88 + t435;
t18 = -pkin(10) * t401 + t22;
t17 = -pkin(10) * t500 + pkin(5) * t177 + (-t158 + (pkin(10) * t231 - t159) * t377) * qJD(5) + t435;
t13 = -mrSges(7,1) * t45 + mrSges(7,2) * t44;
t6 = -qJD(6) * t37 + t17 * t381 - t18 * t376;
t5 = qJD(6) * t36 + t17 * t376 + t18 * t381;
t1 = [m(6) * (t15 * t88 + t16 * t87 + t22 * t76 + t23 * t75) + (t610 / 0.2e1 - Ifges(5,4) * t551 + Ifges(7,3) * t556 + Ifges(6,5) * t558 - t612 / 0.2e1 + t615 / 0.2e1 + t616 / 0.2e1 - t614 / 0.2e1 - t166 / 0.2e1 + Ifges(7,5) * t562 + Ifges(7,6) * t564 + t583) * t177 + (-Ifges(7,5) * t173 - Ifges(7,6) * t172) * t567 + (-Ifges(7,4) * t173 - Ifges(7,2) * t172) * t574 + (-t172 * t3 + t173 * t4 - t20 * t56 + t21 * t57) * mrSges(7,3) + (-Ifges(7,1) * t173 - Ifges(7,4) * t172) * t575 + t33 * (mrSges(7,1) * t172 - mrSges(7,2) * t173) + (Ifges(5,1) * t551 + t613 / 0.2e1 + t212 / 0.2e1 + t167 / 0.2e1 + t440 + t630) * t176 + (Ifges(4,1) * t291 - Ifges(4,4) * t290) * t548 + (Ifges(4,5) * t291 - Ifges(4,6) * t290) * t530 + t325 * (mrSges(4,1) * t290 + mrSges(4,2) * t291) + t288 * (Ifges(4,4) * t291 - Ifges(4,2) * t290) / 0.2e1 + (Ifges(7,5) * t56 + Ifges(7,6) * t57) * t556 + t520 * t624 + t415 * t625 + (-t476 * t579 - t279 * mrSges(6,1) - t257 * mrSges(7,1) - t278 * mrSges(6,2) - t256 * mrSges(7,2) - t594 * (t385 * t317 - t380 * t369) + t589 * t380 + (-t585 + t626) * t385) * g(2) + t213 * (-Ifges(6,5) * t400 - Ifges(6,6) * t401) / 0.2e1 + (-Ifges(6,1) * t400 - Ifges(6,4) * t401) * t558 + (Ifges(7,4) * t56 + Ifges(7,2) * t57) * t564 + (t618 / 0.2e1 + t617 / 0.2e1 - mrSges(3,2) * t539 + Ifges(3,6) * t384 / 0.2e1) * qJDD(2) + (t617 + t618) * t502 + t443 * t651 + t200 * (-Ifges(6,4) * t400 - Ifges(6,2) * t401) / 0.2e1 + (t199 * mrSges(5,2) - t53 * mrSges(5,3) + Ifges(5,1) * t136 + Ifges(5,4) * t137 + Ifges(5,5) * t368 + t111 * t436 + t412 * t566 + t414 * t568 + t416 * t569 + t418 * t49) * t231 + (t322 * t539 + t595) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t595) - t336 * t501 - t40 * t493 / 0.2e1 + m(4) * (t174 * t241 + t175 * t242 + t208 * t227 + t209 * t228 - t275 * t355 + t325 * t359) - t334 * t449 - (mrSges(5,1) * t199 - mrSges(5,3) * t52 - Ifges(5,4) * t136 + Ifges(6,5) * t569 + Ifges(7,5) * t575 - Ifges(5,2) * t137 - Ifges(5,6) * t368 + Ifges(6,6) * t568 + Ifges(7,6) * t574 + Ifges(6,3) * t566 + Ifges(7,3) * t567 + t587 + t588) * t408 - (t455 + t454) * t408 / 0.2e1 + (-t277 * mrSges(6,1) - t255 * mrSges(7,1) - t276 * mrSges(6,2) - t254 * mrSges(7,2) + (t369 * t594 + t589 - t623) * t385 + (-m(7) * (-t317 - t599) - m(6) * (-t317 - t597) + m(5) * t317 + t585) * t380) * g(1) - (-m(5) * t53 + m(6) * t49 - t118 + t50) * t606 + m(7) * (t102 * t54 + t114 * t33 + t20 * t6 + t21 * t5 + t3 * t37 + t36 * t4) + t226 * t359 - t290 * t508 - t297 * t465 / 0.2e1 + m(5) * (t124 * t83 + t165 * t52 + t199 * t263 + t240 * t253) - t402 * t458 + (Ifges(7,1) * t56 + Ifges(7,4) * t57) * t562 - t291 * t509 + t263 * t430 - t355 * t431 + (-m(5) * t123 + m(6) * t116 - t607) * t84 + (t384 * (-Ifges(3,2) * t379 + t520) + t398) * t458 / 0.2e1 + Ifges(2,3) * qJDD(1) + t36 * t26 + t37 * t27 - pkin(1) * (-mrSges(3,1) * t322 + mrSges(3,2) * t323) + t275 * (-mrSges(4,1) * t306 + mrSges(4,2) * t307) + t236 * (Ifges(4,4) * t307 + Ifges(4,2) * t306) + t237 * (Ifges(4,1) * t307 + Ifges(4,4) * t306) - t290 * t216 / 0.2e1 + t291 * t217 / 0.2e1 + t209 * t259 + t208 * t260 + t241 * t223 + t242 * t222 + t253 * t170 + t116 * (mrSges(6,1) * t401 - mrSges(6,2) * t400) + t83 * t204 + t165 * t119 + (-t15 * t493 - t16 * t492 + t400 * t75 - t401 * t76) * mrSges(6,3) + (t413 * t530 - t456) * qJD(2) + (Ifges(3,4) * t624 + Ifges(3,2) * t625 + Ifges(3,6) * t502 + t298 * t530) * t384 + (Ifges(3,1) * t323 - pkin(7) * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t323) + Ifges(3,4) * t625 + 0.2e1 * Ifges(3,5) * t502) * t379 + t54 * t82 + t87 * t61 + t88 * t62 + t56 * t570 + t57 * t572 + t492 * t576 - t173 * t577 - t172 * t578 + (-t174 * t307 + t175 * t306) * mrSges(4,3) + t102 * (-mrSges(7,1) * t57 + mrSges(7,2) * t56) + t5 * t104 + t6 * t105 + t114 * t13 + t22 * t154 + t23 * t155; (t267 * t411 + t280 * t49 - t75 * t80 - t76 * t81 - g(1) * (t303 + t424) - g(2) * (t302 + t425) + t601 * t116) * m(6) + (-t123 * t601 + t124 * t602 - t240 * t252 + t283 * t53 + t284 * t52) * m(5) + t602 * t204 + t603 * t82 + (pkin(7) * t334 + t297 / 0.2e1) * t467 + (-m(4) * t367 - m(7) * (t599 + t468) - m(5) * t468 - m(6) * (t468 + t597) + t591 + t636) * g(3) - t601 * t607 + t631 * (m(4) * t544 - m(5) * t324 + mrSges(4,1) * t361 + mrSges(4,2) * t362 + t419 + t422) + t610 * t552 + t641 + (t456 + (t402 - t398 / 0.2e1) * qJD(1)) * qJD(1) + t586 * t281 + (-t227 * t234 - t228 * t235 - t325 * t358 + (t174 * t374 + t175 * t373) * pkin(2)) * m(4) + t62 * t489 - t226 * t358 + (Ifges(5,4) * t553 - t167 / 0.2e1 + t582) * t600 + t288 * t509 - t413 * t458 / 0.2e1 - t289 * t508 + t654 * t154 + t653 * t155 + t289 * (Ifges(4,1) * t288 + t507) / 0.2e1 + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + t223 * t545 + t222 * t546 + t216 * t548 + (-Ifges(5,4) * t552 + t581 + t166 / 0.2e1) * t409 - (-Ifges(3,2) * t467 + t298 + t357) * t466 / 0.2e1 - (Ifges(4,2) * t289 + t217 + t282) * t288 / 0.2e1 - t325 * (-mrSges(4,1) * t289 + mrSges(4,2) * t288) + Ifges(3,6) * t322 + Ifges(3,5) * t323 - t304 * mrSges(3,2) - t305 * mrSges(3,1) - qJD(2) * (Ifges(4,5) * t288 + Ifges(4,6) * t289) / 0.2e1 + t283 * t118 + t284 * t119 + t280 * t50 + t258 * t13 - t235 * t259 - t234 * t260 - t252 * t170 + Ifges(4,6) * t236 + Ifges(4,5) * t237 + t183 * t26 + t184 * t27 + t174 * mrSges(4,1) - t175 * mrSges(4,2) + t619 * t104 + t620 * t105 + (-g(1) * t303 - g(2) * t302 + t603 * t102 + t183 * t4 + t184 * t3 + t620 * t20 + t619 * t21 + t258 * t33) * m(7); (-t82 + t607) * t409 + t430 + t431 + t382 * t61 + t377 * t62 - t406 * t26 + (-t204 - t410) * t600 + t410 * qJD(5) + t314 * t27 - t289 * t260 - t288 * t259 + t604 * t105 + t605 * t104 + (-g(1) * t380 + g(2) * t385) * (m(4) + t594) + (-t102 * t409 + t20 * t604 + t21 * t605 + t3 * t314 - t4 * t406) * m(7) + (-t116 * t409 + t15 * t377 + t16 * t382 + t213 * t411) * m(6) + (t123 * t409 - t124 * t600 + t199) * m(5) + (-t227 * t289 - t228 * t288 + t275) * m(4); t631 * t419 + (t591 + t626) * g(3) + t641 + t608 * t105 + t609 * t104 + (t102 * t596 + t20 * t608 + t21 * t609 + t246 * t4 + t247 * t3 - t33 * t354) * m(7) + (-t519 + t610) * t552 + t596 * t82 + t586 * pkin(9) + t582 * t600 + (-pkin(4) * t49 - g(1) * t424 - g(2) * t425 - t116 * t124 - t75 * t85 - t76 * t86) * m(6) + (t212 + t167) * t553 + t607 * t124 + t166 * t551 + t62 * t536 + t581 * t409 - t354 * t13 + t246 * t26 + t247 * t27 - t123 * t204 - pkin(4) * t50 - t86 * t154 - t85 * t155; (-mrSges(6,2) * t277 + t276 * t611 - t471) * g(2) + (mrSges(6,2) * t279 - t278 * t611 - t470) * g(1) + t587 + t454 + (t523 - t154) * t75 + (t522 + t155) * t76 + (-Ifges(6,2) * t201 + t111 + t198) * t560 + (-t417 + t418 + t623) * t533 - t82 * t542 - m(7) * (t102 * t542 + t20 * t24 + t21 * t25) + (t381 * t515 - t25) * t104 + (-t376 * t515 - t24) * t105 + t403 + (t26 * t381 + t27 * t376) * pkin(5) + (Ifges(6,5) * t200 - Ifges(6,6) * t201) * t555 - (Ifges(7,4) * t563 + Ifges(7,2) * t565 + Ifges(7,6) * t557 + t573 - t628) * t148 + (Ifges(7,1) * t563 + Ifges(7,4) * t565 + Ifges(7,5) * t557 + t571 - t627) * t428 - t116 * (mrSges(6,1) * t201 + mrSges(6,2) * t200) + t110 * t558 + (Ifges(6,1) * t200 - t510) * t559 + (t3 * t376 + t381 * t4 + (-t20 * t376 + t21 * t381) * qJD(6)) * t579; -t102 * (mrSges(7,1) * t148 + mrSges(7,2) * t428) + (Ifges(7,1) * t428 - t516) * t563 + t71 * t562 + (Ifges(7,5) * t428 - Ifges(7,6) * t148) * t557 - t20 * t104 + t21 * t105 - g(1) * t470 - g(2) * t471 - t417 * t533 + (t148 * t21 + t20 * t428) * mrSges(7,3) + t403 + (-Ifges(7,2) * t148 + t141 + t72) * t565;];
tau  = t1;
