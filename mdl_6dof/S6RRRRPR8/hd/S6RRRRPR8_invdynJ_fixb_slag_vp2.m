% Calculate vector of inverse dynamics joint torques for
% S6RRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPR8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR8_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR8_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR8_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR8_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR8_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:38:52
% EndTime: 2019-03-09 22:40:15
% DurationCPUTime: 48.44s
% Computational Cost: add. (17720->978), mult. (38023->1245), div. (0->0), fcn. (26571->12), ass. (0->455)
t377 = sin(qJ(2));
t546 = Ifges(3,4) * t377;
t382 = cos(qJ(2));
t564 = t382 / 0.2e1;
t711 = Ifges(3,2) * t564 + t546 / 0.2e1;
t665 = mrSges(5,1) + mrSges(6,1);
t664 = mrSges(5,2) - mrSges(6,3);
t380 = cos(qJ(4));
t489 = qJD(4) * t380;
t448 = pkin(2) * t382 + pkin(8) * t377;
t318 = -pkin(1) - t448;
t289 = t318 * qJD(1);
t497 = qJD(1) * t382;
t363 = pkin(7) * t497;
t326 = qJD(2) * pkin(8) + t363;
t376 = sin(qJ(3));
t381 = cos(qJ(3));
t212 = t381 * t289 - t326 * t376;
t498 = qJD(1) * t377;
t471 = t381 * t498;
t300 = qJD(2) * t376 + t471;
t169 = -pkin(9) * t300 + t212;
t213 = t376 * t289 + t381 * t326;
t466 = t376 * t498;
t495 = qJD(2) * t381;
t413 = t466 - t495;
t170 = -pkin(9) * t413 + t213;
t375 = sin(qJ(4));
t523 = t375 * t170;
t92 = t380 * t169 - t523;
t717 = pkin(3) * t489 - t92;
t301 = t375 * t376 - t380 * t381;
t623 = qJD(3) + qJD(4);
t217 = t623 * t301;
t410 = t301 * t382;
t245 = qJD(1) * t410;
t716 = t217 - t245;
t302 = t375 * t381 + t376 * t380;
t218 = t623 * t302;
t244 = t302 * t497;
t694 = t218 - t244;
t345 = qJD(3) - t497;
t335 = qJD(4) + t345;
t567 = -t335 / 0.2e1;
t394 = t380 * t300 - t375 * t413;
t579 = -t394 / 0.2e1;
t659 = Ifges(6,2) + Ifges(5,3);
t660 = Ifges(6,4) + Ifges(5,5);
t715 = t567 * t659 + t579 * t660;
t321 = qJD(6) - t335;
t568 = t321 / 0.2e1;
t204 = t375 * t300 + t380 * t413;
t374 = sin(qJ(6));
t379 = cos(qJ(6));
t426 = t204 * t374 + t379 * t394;
t586 = t426 / 0.2e1;
t123 = t204 * t379 - t374 * t394;
t588 = t123 / 0.2e1;
t714 = Ifges(7,5) * t586 + Ifges(7,6) * t588 + Ifges(7,3) * t568;
t582 = -t204 / 0.2e1;
t713 = -t300 / 0.2e1;
t712 = -t345 / 0.2e1;
t668 = t413 / 0.2e1;
t662 = Ifges(5,1) + Ifges(6,1);
t490 = qJD(4) * t375;
t513 = t380 * t170;
t91 = t169 * t375 + t513;
t710 = pkin(3) * t490 - t91;
t643 = qJD(5) + t717;
t447 = pkin(2) * t377 - pkin(8) * t382;
t306 = t447 * qJD(1);
t223 = pkin(7) * t466 + t381 * t306;
t511 = t381 * t382;
t421 = pkin(3) * t377 - pkin(9) * t511;
t188 = qJD(1) * t421 + t223;
t281 = t376 * t306;
t518 = t377 * t381;
t520 = t376 * t382;
t209 = t281 + (-pkin(7) * t518 - pkin(9) * t520) * qJD(1);
t118 = t375 * t188 + t380 * t209;
t108 = qJ(5) * t498 + t118;
t384 = -pkin(9) - pkin(8);
t473 = qJD(3) * t384;
t307 = t376 * t473;
t308 = t381 * t473;
t327 = t384 * t376;
t328 = t384 * t381;
t145 = t380 * t307 + t375 * t308 + t327 * t489 + t328 * t490;
t709 = t145 - t108;
t117 = t188 * t380 - t375 * t209;
t222 = t375 * t327 - t380 * t328;
t146 = qJD(4) * t222 + t307 * t375 - t380 * t308;
t708 = -t146 - t117;
t385 = -pkin(4) - pkin(5);
t666 = pkin(10) * t394;
t149 = pkin(3) * t345 + t169;
t85 = t380 * t149 - t523;
t60 = t85 + t666;
t53 = t335 * t385 + qJD(5) - t60;
t317 = t335 * qJ(5);
t682 = pkin(10) * t204;
t86 = t375 * t149 + t513;
t61 = t86 + t682;
t57 = t317 + t61;
t21 = t374 * t53 + t379 * t57;
t558 = mrSges(7,3) * t21;
t362 = pkin(7) * t498;
t325 = -qJD(2) * pkin(2) + t362;
t227 = pkin(3) * t413 + t325;
t390 = qJ(5) * t394 - t227;
t71 = t204 * t385 + t390;
t707 = -mrSges(7,1) * t71 + t558;
t449 = m(7) * (-pkin(10) - t384) - mrSges(7,3);
t701 = -mrSges(5,3) - mrSges(6,2);
t625 = t701 * t377;
t706 = -t449 * t377 + t625;
t566 = t335 / 0.2e1;
t578 = t394 / 0.2e1;
t581 = t204 / 0.2e1;
t74 = -pkin(4) * t335 + qJD(5) - t85;
t94 = t204 * pkin(4) - t390;
t693 = t227 * mrSges(5,2) + mrSges(6,2) * t74 - mrSges(5,3) * t85 - t94 * mrSges(6,3);
t705 = Ifges(5,4) * t582 + Ifges(6,5) * t581 + t566 * t660 + t578 * t662 + t693;
t20 = -t374 * t57 + t379 * t53;
t559 = mrSges(7,3) * t20;
t704 = mrSges(7,2) * t71 - t559;
t661 = -Ifges(5,4) + Ifges(6,5);
t658 = -Ifges(5,6) + Ifges(6,6);
t478 = t385 * t377;
t700 = pkin(10) * t716 - qJD(1) * t478 - t708;
t699 = -pkin(10) * t694 - t709;
t488 = qJD(1) * qJD(2);
t310 = t382 * qJDD(1) - t377 * t488;
t311 = qJDD(1) * t377 + t382 * t488;
t531 = qJDD(1) * pkin(1);
t216 = -pkin(2) * t310 - pkin(8) * t311 - t531;
t295 = t310 * pkin(7);
t267 = qJDD(2) * pkin(8) + t295;
t491 = qJD(3) * t381;
t493 = qJD(3) * t376;
t100 = t376 * t216 + t381 * t267 + t289 * t491 - t326 * t493;
t698 = t100 * mrSges(4,2);
t101 = -qJD(3) * t213 + t381 * t216 - t267 * t376;
t697 = t101 * mrSges(4,1);
t696 = -t666 + t643;
t695 = -t682 + t710;
t472 = t376 * t497;
t633 = -t363 + (-t472 + t493) * pkin(3);
t120 = Ifges(7,4) * t123;
t692 = Ifges(7,2) * t426 - t120;
t198 = Ifges(6,5) * t394;
t111 = t335 * Ifges(6,6) + t204 * Ifges(6,3) + t198;
t541 = Ifges(5,4) * t394;
t114 = -t204 * Ifges(5,2) + t335 * Ifges(5,6) + t541;
t75 = t317 + t86;
t691 = -mrSges(6,2) * t75 - mrSges(5,3) * t86 + t227 * mrSges(5,1) + t94 * mrSges(6,1) + t111 / 0.2e1 - t114 / 0.2e1;
t373 = qJ(3) + qJ(4);
t366 = sin(t373);
t367 = cos(t373);
t383 = cos(qJ(1));
t378 = sin(qJ(1));
t514 = t378 * t382;
t256 = t366 * t514 + t367 * t383;
t509 = t383 * t366;
t257 = t367 * t514 - t509;
t425 = t256 * t374 + t257 * t379;
t627 = t256 * t379 - t257 * t374;
t690 = mrSges(7,1) * t627 - t425 * mrSges(7,2);
t423 = t366 * t374 + t367 * t379;
t424 = t366 * t379 - t367 * t374;
t441 = -mrSges(4,1) * t381 + mrSges(4,2) * t376;
t689 = m(4) * pkin(2) + t423 * mrSges(7,1) + t424 * mrSges(7,2) - t366 * t664 + t367 * t665 - t441;
t56 = Ifges(7,1) * t426 + Ifges(7,5) * t321 + t120;
t688 = Ifges(7,1) * t586 + Ifges(7,4) * t588 + Ifges(7,5) * t568 + t56 / 0.2e1;
t540 = Ifges(7,4) * t426;
t55 = Ifges(7,2) * t123 + Ifges(7,6) * t321 + t540;
t687 = Ifges(7,4) * t586 + Ifges(7,2) * t588 + Ifges(7,6) * t568 + t55 / 0.2e1;
t569 = -t321 / 0.2e1;
t587 = -t426 / 0.2e1;
t686 = (t123 * t20 + t21 * t426) * mrSges(7,3) + (Ifges(7,5) * t123 - Ifges(7,6) * t426) * t569 + (Ifges(7,1) * t123 - t540) * t587 - t71 * (mrSges(7,1) * t426 + mrSges(7,2) * t123);
t684 = t20 * mrSges(7,1) + t74 * mrSges(6,1) + t86 * mrSges(5,2) + Ifges(4,5) * t713 + Ifges(4,6) * t668 + Ifges(4,3) * t712 + t658 * t582 + Ifges(3,6) * qJD(2) / 0.2e1 + qJD(1) * t711 - t21 * mrSges(7,2) - t75 * mrSges(6,3) - t85 * mrSges(5,1) + t714 + t715;
t199 = Ifges(5,4) * t204;
t539 = Ifges(6,5) * t204;
t639 = t335 * t660 + t394 * t662 - t199 + t539;
t683 = t639 / 0.2e1;
t533 = qJ(5) * t204;
t547 = mrSges(5,3) * t394;
t175 = mrSges(5,1) * t335 - t547;
t176 = -mrSges(6,1) * t335 + mrSges(6,2) * t394;
t636 = -t175 + t176;
t681 = qJ(5) * t716 - qJD(5) * t302 + t633;
t494 = qJD(2) * t382;
t470 = t376 * t494;
t405 = t377 * t491 + t470;
t510 = t382 * t383;
t278 = -t376 * t510 + t378 * t381;
t678 = -Ifges(5,2) * t582 + Ifges(6,3) * t581 + t566 * t658 + t578 * t661 + t691;
t677 = -Ifges(5,2) * t581 + Ifges(6,3) * t582 + t567 * t658 + t579 * t661 - t691;
t676 = Ifges(5,4) * t581 + Ifges(6,5) * t582 + t567 * t660 + t579 * t662 - t693;
t401 = t413 * qJD(3);
t195 = qJDD(2) * t376 + t311 * t381 - t401;
t196 = -qJD(3) * t300 + qJDD(2) * t381 - t311 * t376;
t82 = -qJD(4) * t204 + t380 * t195 + t375 * t196;
t83 = qJD(4) * t394 + t375 * t195 - t380 * t196;
t25 = qJD(6) * t123 + t374 * t83 + t379 * t82;
t605 = t25 / 0.2e1;
t26 = -qJD(6) * t426 - t374 * t82 + t379 * t83;
t604 = t26 / 0.2e1;
t597 = t82 / 0.2e1;
t595 = t83 / 0.2e1;
t670 = -m(3) - m(4);
t669 = -m(7) - m(6);
t584 = t195 / 0.2e1;
t583 = t196 / 0.2e1;
t297 = qJDD(3) - t310;
t288 = qJDD(4) + t297;
t271 = qJDD(6) - t288;
t574 = t271 / 0.2e1;
t573 = t288 / 0.2e1;
t572 = t297 / 0.2e1;
t667 = pkin(4) * t394;
t663 = -mrSges(3,3) + mrSges(2,2);
t657 = t288 * t660 + t661 * t83 + t662 * t82;
t221 = -t380 * t327 - t328 * t375;
t182 = -pkin(10) * t302 + t221;
t183 = pkin(10) * t301 + t222;
t105 = t182 * t374 + t183 * t379;
t656 = -qJD(6) * t105 + t374 * t699 + t379 * t700;
t104 = t182 * t379 - t183 * t374;
t655 = qJD(6) * t104 + t374 * t700 - t379 * t699;
t653 = mrSges(3,2) * t382;
t648 = t385 * t694 - t681;
t555 = pkin(3) * t380;
t357 = -pkin(4) - t555;
t351 = -pkin(5) + t357;
t556 = pkin(3) * t375;
t352 = qJ(5) + t556;
t252 = t351 * t379 - t352 * t374;
t647 = qJD(6) * t252 + t374 * t695 + t379 * t696;
t253 = t351 * t374 + t352 * t379;
t646 = -qJD(6) * t253 - t374 * t696 + t379 * t695;
t264 = t302 * t377;
t640 = t385 * t394;
t638 = pkin(4) * t694 + t681;
t258 = -t378 * t367 + t382 * t509;
t259 = t366 * t378 + t367 * t510;
t171 = t258 * t379 - t259 * t374;
t172 = t258 * t374 + t259 * t379;
t507 = -mrSges(7,1) * t171 + mrSges(7,2) * t172;
t173 = -mrSges(6,2) * t204 + mrSges(6,3) * t335;
t548 = mrSges(5,3) * t204;
t174 = -mrSges(5,2) * t335 - t548;
t637 = t173 + t174;
t635 = (mrSges(7,1) * t424 - mrSges(7,2) * t423) * t377;
t293 = Ifges(4,4) * t413;
t187 = t300 * Ifges(4,1) + t345 * Ifges(4,5) - t293;
t361 = Ifges(3,4) * t497;
t634 = Ifges(3,1) * t498 + Ifges(3,5) * qJD(2) + t381 * t187 + t361;
t632 = qJD(2) * mrSges(3,1) - mrSges(4,1) * t413 - t300 * mrSges(4,2) - mrSges(3,3) * t498;
t370 = t382 * pkin(4);
t532 = qJ(5) * t382;
t631 = t366 * t532 + t367 * t370;
t524 = t367 * t377;
t630 = -mrSges(6,3) * t524 + t635;
t349 = pkin(7) * t511;
t235 = t376 * t318 + t349;
t492 = qJD(3) * t377;
t406 = -t376 * t492 + t381 * t494;
t629 = t288 * t659 + t658 * t83 + t660 * t82;
t296 = t311 * pkin(7);
t628 = t295 * t382 + t296 * t377;
t626 = t100 * t381 - t101 * t376;
t443 = t382 * mrSges(3,1) - t377 * mrSges(3,2);
t622 = t377 * mrSges(4,3) + mrSges(2,1) + t443;
t620 = t258 * t665 + t259 * t664 - t507;
t619 = t256 * t665 + t257 * t664 + t690;
t618 = m(7) * pkin(5) + t665;
t62 = pkin(3) * t297 - pkin(9) * t195 + t101;
t72 = pkin(9) * t196 + t100;
t16 = t149 * t489 - t170 * t490 + t375 * t62 + t380 * t72;
t13 = t288 * qJ(5) + t335 * qJD(5) + t16;
t446 = qJDD(2) * pkin(2) - t296;
t408 = pkin(3) * t196 + t446;
t392 = qJ(5) * t82 + qJD(5) * t394 + t408;
t22 = pkin(4) * t83 - t392;
t596 = -t83 / 0.2e1;
t610 = -mrSges(5,1) * t408 + mrSges(6,1) * t22 - mrSges(6,2) * t13 - mrSges(5,3) * t16 + 0.2e1 * Ifges(6,3) * t595 - t82 * Ifges(5,4) / 0.2e1 - t288 * Ifges(5,6) / 0.2e1 + (t661 + Ifges(6,5)) * t597 + (t658 + Ifges(6,6)) * t573 + (-t596 + t595) * Ifges(5,2);
t608 = Ifges(7,4) * t605 + Ifges(7,2) * t604 + Ifges(7,6) * t574;
t607 = Ifges(7,1) * t605 + Ifges(7,4) * t604 + Ifges(7,5) * t574;
t606 = m(5) * pkin(3);
t601 = -t55 / 0.2e1;
t599 = -t56 / 0.2e1;
t594 = Ifges(4,1) * t584 + Ifges(4,4) * t583 + Ifges(4,5) * t572;
t589 = -t123 / 0.2e1;
t571 = t300 / 0.2e1;
t557 = pkin(3) * t300;
t554 = pkin(5) * t367;
t368 = t377 * pkin(7);
t551 = -qJD(1) / 0.2e1;
t550 = qJD(3) / 0.2e1;
t545 = Ifges(3,4) * t382;
t544 = Ifges(4,4) * t300;
t543 = Ifges(4,4) * t376;
t542 = Ifges(4,4) * t381;
t538 = t213 * mrSges(4,3);
t537 = t367 * mrSges(5,2);
t525 = t366 * t377;
t522 = t376 * t377;
t521 = t376 * t378;
t519 = t376 * t383;
t516 = t377 * t384;
t358 = pkin(3) * t381 + pkin(2);
t331 = t382 * t358;
t299 = t381 * t318;
t211 = -pkin(9) * t518 + t299 + (-pkin(7) * t376 - pkin(3)) * t382;
t220 = -pkin(9) * t522 + t235;
t136 = t375 * t211 + t380 * t220;
t309 = t447 * qJD(2);
t496 = qJD(2) * t377;
t480 = pkin(7) * t496;
t501 = t381 * t309 + t376 * t480;
t329 = qJ(5) * t524;
t346 = pkin(3) * t522;
t500 = t329 - t346;
t312 = t368 + t346;
t499 = t383 * pkin(1) + t378 * pkin(7);
t486 = Ifges(7,5) * t25 + Ifges(7,6) * t26 + Ifges(7,3) * t271;
t364 = pkin(7) * t494;
t479 = m(4) * pkin(8) + mrSges(4,3);
t476 = t383 * t516;
t475 = Ifges(4,5) * t195 + Ifges(4,6) * t196 + Ifges(4,3) * t297;
t371 = t383 * pkin(7);
t474 = pkin(3) * t519 + t378 * t516 + t371;
t230 = pkin(3) * t405 + t364;
t186 = -Ifges(4,2) * t413 + Ifges(4,6) * t345 + t544;
t465 = -t376 * t186 / 0.2e1;
t66 = -t288 * mrSges(6,1) + t82 * mrSges(6,2);
t454 = -t256 * pkin(4) + qJ(5) * t257;
t453 = -t258 * pkin(4) + qJ(5) * t259;
t135 = t211 * t380 - t375 * t220;
t451 = t331 - t516;
t450 = pkin(3) * t521 + t358 * t510 + t499;
t128 = t136 - t532;
t265 = t301 * t377;
t444 = -qJ(5) * t265 - t312;
t129 = -t135 + t370;
t442 = mrSges(3,1) * t377 + t653;
t440 = mrSges(4,1) * t376 + mrSges(4,2) * t381;
t437 = Ifges(4,1) * t381 - t543;
t436 = Ifges(4,1) * t376 + t542;
t434 = -Ifges(4,2) * t376 + t542;
t433 = Ifges(4,2) * t381 + t543;
t432 = Ifges(3,5) * t382 - Ifges(3,6) * t377;
t431 = Ifges(4,5) * t381 - t376 * Ifges(4,6);
t430 = Ifges(4,5) * t376 + Ifges(4,6) * t381;
t429 = -t557 - t533;
t93 = pkin(5) * t382 + pkin(10) * t265 + t129;
t95 = pkin(10) * t264 + t128;
t48 = -t374 * t95 + t379 * t93;
t49 = t374 * t93 + t379 * t95;
t315 = t379 * qJ(5) + t374 * t385;
t314 = -t374 * qJ(5) + t379 * t385;
t102 = -mrSges(7,2) * t321 + mrSges(7,3) * t123;
t103 = mrSges(7,1) * t321 - mrSges(7,3) * t426;
t427 = t102 * t379 - t103 * t374;
t177 = t264 * t379 + t265 * t374;
t178 = t264 * t374 - t265 * t379;
t207 = t301 * t379 - t302 * t374;
t208 = t301 * t374 + t302 * t379;
t422 = t278 * pkin(3);
t17 = -t149 * t490 - t170 * t489 - t375 * t72 + t380 * t62;
t420 = qJ(5) * t302 + t358;
t134 = t421 * qJD(2) + (-t349 + (pkin(9) * t377 - t318) * t376) * qJD(3) + t501;
t167 = t376 * t309 + t318 * t491 + (-t377 * t495 - t382 * t493) * pkin(7);
t139 = -pkin(9) * t405 + t167;
t47 = t134 * t380 - t375 * t139 - t211 * t490 - t220 * t489;
t10 = pkin(10) * t83 + t13;
t411 = qJDD(5) - t17;
t9 = -pkin(10) * t82 + t288 * t385 + t411;
t2 = qJD(6) * t20 + t10 * t379 + t374 * t9;
t3 = -qJD(6) * t21 - t10 * t374 + t379 * t9;
t419 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + t486;
t418 = -pkin(4) * t367 - qJ(5) * t366 - t358;
t416 = pkin(1) * t442;
t276 = t376 * t514 + t381 * t383;
t415 = t325 * t440;
t414 = t377 * (Ifges(3,1) * t382 - t546);
t46 = t375 * t134 + t380 * t139 + t211 * t489 - t220 * t490;
t407 = t276 * pkin(3);
t404 = t259 * pkin(4) + qJ(5) * t258 + t450;
t403 = t413 * mrSges(4,3);
t402 = -g(1) * t258 - g(2) * t256 - g(3) * t525;
t140 = -qJD(2) * t410 - t264 * t623;
t400 = qJ(5) * t140 - qJD(5) * t265 - t230;
t398 = t422 + t453;
t397 = Ifges(4,5) * t377 + t382 * t437;
t396 = Ifges(4,6) * t377 + t382 * t434;
t395 = Ifges(4,3) * t377 + t382 * t431;
t44 = qJ(5) * t496 - qJD(5) * t382 + t46;
t393 = -t407 + t454;
t14 = -pkin(4) * t288 + t411;
t388 = t17 * mrSges(5,1) - t14 * mrSges(6,1) - t16 * mrSges(5,2) + t13 * mrSges(6,3) - t419 + t629;
t323 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t497;
t290 = t440 * t377;
t279 = t381 * t510 + t521;
t277 = -t378 * t511 + t519;
t261 = -t374 * qJD(5) - qJD(6) * t315;
t260 = t379 * qJD(5) + qJD(6) * t314;
t249 = t258 * pkin(5);
t246 = t256 * pkin(5);
t234 = -pkin(7) * t520 + t299;
t226 = mrSges(4,1) * t345 - mrSges(4,3) * t300;
t225 = -t345 * mrSges(4,2) - t403;
t224 = -pkin(7) * t471 + t281;
t200 = pkin(4) * t301 - t420;
t168 = -qJD(3) * t235 + t501;
t162 = t244 * t374 - t245 * t379;
t161 = t244 * t379 + t245 * t374;
t160 = t301 * t385 + t420;
t153 = pkin(4) * t264 - t444;
t152 = -mrSges(4,2) * t297 + mrSges(4,3) * t196;
t151 = mrSges(4,1) * t297 - mrSges(4,3) * t195;
t141 = -t490 * t522 + (t518 * t623 + t470) * t380 + t406 * t375;
t133 = mrSges(5,1) * t204 + mrSges(5,2) * t394;
t132 = mrSges(6,1) * t204 - mrSges(6,3) * t394;
t131 = t533 + t667;
t130 = t264 * t385 + t444;
t119 = -mrSges(4,1) * t196 + mrSges(4,2) * t195;
t109 = -pkin(4) * t498 - t117;
t107 = -t429 + t667;
t98 = t195 * Ifges(4,4) + t196 * Ifges(4,2) + t297 * Ifges(4,6);
t90 = -qJD(6) * t208 + t217 * t374 + t218 * t379;
t89 = qJD(6) * t207 - t217 * t379 + t218 * t374;
t88 = -t533 + t640;
t73 = t429 + t640;
t68 = -mrSges(6,2) * t83 + mrSges(6,3) * t288;
t67 = -mrSges(5,2) * t288 - mrSges(5,3) * t83;
t65 = mrSges(5,1) * t288 - mrSges(5,3) * t82;
t58 = -mrSges(7,1) * t123 + mrSges(7,2) * t426;
t52 = pkin(4) * t141 - t400;
t45 = -pkin(4) * t496 - t47;
t41 = t141 * t385 + t400;
t40 = mrSges(5,1) * t83 + mrSges(5,2) * t82;
t39 = mrSges(6,1) * t83 - mrSges(6,3) * t82;
t30 = pkin(10) * t141 + t44;
t29 = -pkin(10) * t140 + qJD(2) * t478 - t47;
t28 = t374 * t61 + t379 * t60;
t27 = -t374 * t60 + t379 * t61;
t19 = -mrSges(7,2) * t271 + mrSges(7,3) * t26;
t18 = mrSges(7,1) * t271 - mrSges(7,3) * t25;
t11 = t385 * t83 + t392;
t8 = -mrSges(7,1) * t26 + mrSges(7,2) * t25;
t5 = -qJD(6) * t49 + t29 * t379 - t30 * t374;
t4 = qJD(6) * t48 + t29 * t374 + t30 * t379;
t1 = [-t382 * t697 + t382 * t698 + (t212 * mrSges(4,1) - t213 * mrSges(4,2) + Ifges(5,6) * t582 + Ifges(6,6) * t581 + t659 * t566 + t660 * t578 - t684 - t714) * t496 + (-Ifges(6,5) * t265 - Ifges(6,6) * t382) * t595 + (-Ifges(5,4) * t265 - Ifges(5,6) * t382) * t596 + t17 * (-mrSges(5,1) * t382 + mrSges(5,3) * t265) + t14 * (mrSges(6,1) * t382 - mrSges(6,2) * t265) + t678 * t141 + t610 * t264 + (pkin(7) * t382 * mrSges(3,3) + pkin(1) * mrSges(3,1) + 0.2e1 * t711) * t310 + m(5) * (t135 * t17 + t136 * t16 + t227 * t230 - t312 * t408 + t46 * t86 + t47 * t85) + m(4) * (t100 * t235 + t101 * t234 + t167 * t213 + t168 * t212 + (t325 * t494 - t377 * t446) * pkin(7)) - t446 * t290 + t443 * t531 + (-m(5) * t474 - t277 * mrSges(4,1) + t425 * mrSges(7,1) - t276 * mrSges(4,2) + t627 * mrSges(7,2) + t669 * (-t257 * pkin(4) - qJ(5) * t256 + t474) + t663 * t383 + t670 * t371 + t618 * t257 - t664 * t256 + (m(3) * pkin(1) - m(4) * t318 + (-m(7) * pkin(10) - mrSges(7,3)) * t377 + (-m(5) + t669) * (-pkin(1) - t331) + t622 - t625) * t378) * g(1) + t312 * t40 - t413 * (qJD(2) * t396 - t433 * t492) / 0.2e1 + t345 * (qJD(2) * t395 - t430 * t492) / 0.2e1 + (t377 * Ifges(3,1) + t545 / 0.2e1 - pkin(1) * mrSges(3,2) + mrSges(3,3) * t368 + Ifges(3,4) * t564) * t311 + (t16 * t382 + t265 * t408) * mrSges(5,2) - t98 * t522 / 0.2e1 - t416 * t488 + (-Ifges(4,6) * t382 + t377 * t434) * t583 + (-Ifges(4,5) * t382 + t377 * t437) * t584 + (qJD(2) * t397 - t436 * t492) * t571 + (-Ifges(4,3) * t382 + t377 * t431) * t572 + (Ifges(7,5) * t178 + Ifges(7,6) * t177 + Ifges(7,3) * t382) * t574 + t486 * t564 + (-t100 * t522 - t101 * t518 - t212 * t406 - t213 * t405) * mrSges(4,3) + (-t13 * t382 + t22 * t265) * mrSges(6,3) + t119 * t368 + (t688 + t704) * (qJD(6) * t177 + t140 * t379 + t141 * t374) + (t683 + t705) * t140 + (-m(5) * (t450 - t476) - m(6) * (t404 - t476) - t279 * mrSges(4,1) - t278 * mrSges(4,2) - m(7) * t404 - t172 * mrSges(7,1) - t171 * mrSges(7,2) + t670 * t499 + t663 * t378 - t618 * t259 + t664 * t258 + (-m(4) * t448 - t622 + t706) * t383) * g(2) + (t687 + t707) * (-qJD(6) * t178 - t140 * t374 + t141 * t379) - t323 * t480 + (t382 * (-Ifges(3,2) * t377 + t545) + t414) * t488 / 0.2e1 - (t186 * t381 + t187 * t376) * t492 / 0.2e1 + t234 * t151 + t235 * t152 + t167 * t225 + t168 * t226 + t230 * t133 + qJD(2) ^ 2 * t432 / 0.2e1 + Ifges(2,3) * qJDD(1) + t45 * t176 + t11 * (-mrSges(7,1) * t177 + mrSges(7,2) * t178) + t44 * t173 + t46 * t174 + t47 * t175 + t153 * t39 + t135 * t65 + t136 * t67 + t128 * t68 + t129 * t66 + t130 * t8 + t52 * t132 + t4 * t102 + t5 * t103 + t465 * t494 + t41 * t58 + t48 * t18 + t49 * t19 + m(6) * (t128 * t13 + t129 * t14 + t153 * t22 + t44 * t75 + t45 * t74 + t52 * t94) + m(7) * (t11 * t130 + t2 * t49 + t20 * t5 + t21 * t4 + t3 * t48 + t41 * t71) + t325 * (mrSges(4,1) * t405 + mrSges(4,2) * t406) + (-mrSges(3,1) * t368 + Ifges(3,5) * t377 + 0.2e1 * Ifges(3,6) * t564 - pkin(7) * t653) * qJDD(2) + t518 * t594 + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t628) + t628 * mrSges(3,3) - (t475 + t629) * t382 / 0.2e1 - t632 * t364 + t634 * t494 / 0.2e1 + (Ifges(7,4) * t178 + Ifges(7,2) * t177 + Ifges(7,6) * t382) * t604 + (Ifges(7,1) * t178 + Ifges(7,4) * t177 + Ifges(7,5) * t382) * t605 + t178 * t607 + t177 * t608 + t2 * (-mrSges(7,2) * t382 + mrSges(7,3) * t177) + t3 * (mrSges(7,1) * t382 - mrSges(7,3) * t178) - t657 * t265 / 0.2e1 + (-t265 * t660 - t382 * t659) * t573 + (-t265 * t662 - t382 * t660) * t597; (-t212 * (mrSges(4,1) * t377 - mrSges(4,3) * t511) - t213 * (-mrSges(4,2) * t377 - mrSges(4,3) * t520) + t396 * t668 + (-t414 / 0.2e1 + t416) * qJD(1)) * qJD(1) + (Ifges(7,1) * t162 + Ifges(7,4) * t161) * t587 + (t395 * t551 + t431 * t550) * t345 + t677 * t244 + t678 * t218 + (t684 + Ifges(3,2) * t497 / 0.2e1 + Ifges(5,6) * t581 + Ifges(6,6) * t582 - Ifges(7,5) * t587 - Ifges(7,6) * t589 - Ifges(7,3) * t569 + t715) * t498 + (t415 + t465) * qJD(3) + (g(1) * t383 + g(2) * t378) * (t442 + (-t449 - t479 + (m(5) + m(6)) * t384 + t701) * t382 + (m(5) * t358 - m(6) * t418 - m(7) * (t418 - t554) + t689) * t377) + t610 * t301 + t639 * (t245 / 0.2e1 - t217 / 0.2e1) + ((-t162 + t89) * mrSges(7,2) + (t161 - t90) * mrSges(7,1)) * t71 + (t657 / 0.2e1 - mrSges(5,2) * t408 + mrSges(6,2) * t14 - mrSges(5,3) * t17 - mrSges(6,3) * t22 + Ifges(5,4) * t596 + Ifges(6,5) * t595 + t573 * t660 + t597 * t662) * t302 + (pkin(2) * t446 - t212 * t223 - t213 * t224 - t325 * t363) * m(4) - t446 * t441 - t415 * t497 + Ifges(3,6) * t310 + Ifges(3,5) * t311 + (t68 + t67) * t222 + (t66 - t65) * t221 + (Ifges(7,5) * t162 + Ifges(7,6) * t161) * t569 + t187 * t491 / 0.2e1 + (t397 * t551 + t437 * t550) * t300 - t676 * t245 - (t361 + t634) * t497 / 0.2e1 - t432 * t488 / 0.2e1 + (t558 + t687) * t90 + (-t559 + t688) * t89 + (t408 * t358 + t16 * t222 - t17 * t221 + (t145 - t118) * t86 + t708 * t85 + t633 * t227) * m(5) + (t13 * t222 + t14 * t221 + t200 * t22 + t638 * t94 + t709 * t75 + (t146 - t109) * t74) * m(6) + t433 * t583 + t436 * t584 + t430 * t572 + (Ifges(7,5) * t208 + Ifges(7,6) * t207) * t574 + (-t161 * t21 + t162 * t20 + t2 * t207 - t208 * t3) * mrSges(7,3) - t493 * t538 + t186 * t472 / 0.2e1 - t705 * t217 + (-t443 - t377 * t479 - m(7) * (t331 + t631) - m(5) * t451 - m(6) * (t451 + t631) + (-m(7) * t554 - t689) * t382 + t706) * g(3) - t224 * t225 - t223 * t226 - t434 * t401 / 0.2e1 + t11 * (-mrSges(7,1) * t207 + mrSges(7,2) * t208) + t200 * t39 - t108 * t173 - t118 * t174 - t117 * t175 - t109 * t176 + t160 * t8 - t295 * mrSges(3,2) - t296 * mrSges(3,1) - pkin(2) * t119 + t104 * t18 + t105 * t19 + Ifges(3,3) * qJDD(2) + t376 * t594 + t162 * t599 + t161 * t601 + (Ifges(7,4) * t162 + Ifges(7,2) * t161) * t589 - t358 * t40 + (m(4) * ((-t212 * t381 - t213 * t376) * qJD(3) + t626) - t225 * t493 - t226 * t491 - t376 * t151 + t381 * t152) * pkin(8) + (-t212 * t491 + t626) * mrSges(4,3) + t323 * t362 + t632 * t363 + t633 * t133 + t636 * t146 + t637 * t145 + t638 * t132 + (Ifges(7,4) * t208 + Ifges(7,2) * t207) * t604 + (Ifges(7,1) * t208 + Ifges(7,4) * t207) * t605 + t208 * t607 + t207 * t608 + t648 * t58 + t381 * t98 / 0.2e1 + t655 * t102 + t656 * t103 + (t104 * t3 + t105 * t2 + t11 * t160 + t20 * t656 + t21 * t655 + t648 * t71) * m(7); (-t225 - t403) * t212 + (-Ifges(4,5) * t413 - Ifges(4,6) * t300) * t712 + (-Ifges(4,1) * t413 - t544) * t713 + t677 * t394 + (-t107 * t94 + t13 * t352 + t14 * t357 + t643 * t75) * m(6) + (-Ifges(4,2) * t300 + t187 - t293) * t668 + t717 * t174 + t475 - t698 + t186 * t571 + t65 * t555 + t67 * t556 - t133 * t557 - m(5) * (t227 * t557 - t85 * t91 + t86 * t92) + (m(6) * t74 + t636) * t710 + t300 * t538 - (Ifges(7,1) * t587 + Ifges(7,4) * t589 + Ifges(7,5) * t569 + t599 - t704) * t123 + (Ifges(7,4) * t587 + Ifges(7,2) * t589 + Ifges(7,6) * t569 + t601 - t707) * t426 + t697 + t252 * t18 + t253 * t19 + t213 * t226 - t107 * t132 - t73 * t58 - t325 * (t300 * mrSges(4,1) - mrSges(4,2) * t413) + (t683 - t676) * t204 + (-m(6) * t393 + m(5) * t407 + mrSges(4,1) * t276 - mrSges(4,2) * t277 - m(7) * (-t246 + t393) + t619) * g(2) + (-m(6) * t398 - m(5) * t422 - mrSges(4,1) * t278 + mrSges(4,2) * t279 - m(7) * (-t249 + t398) + t620) * g(1) + t352 * t68 + t357 * t66 + t388 + (t290 - m(6) * (-pkin(4) * t525 + t500) + mrSges(6,1) * t525 - (-t366 * mrSges(5,1) - t376 * t606 - t537) * t377 - m(7) * (t366 * t478 + t500) + t630) * g(3) + (t16 * t375 + t17 * t380 + (-t375 * t85 + t380 * t86) * qJD(4)) * t606 + t643 * t173 + t646 * t103 + t647 * t102 + (t2 * t253 + t20 * t646 + t21 * t647 + t252 * t3 - t71 * t73) * m(7); -t686 + t692 * t589 + (t204 * t74 + t394 * t75) * mrSges(6,2) + (Ifges(6,3) * t394 - t539) * t582 - t227 * (mrSges(5,1) * t394 - mrSges(5,2) * t204) - t94 * (mrSges(6,1) * t394 + mrSges(6,3) * t204) - t123 * t599 + t314 * t18 + t315 * t19 + t114 * t578 - m(6) * (t131 * t94 + t74 * t86 + t75 * t85) - m(7) * (t20 * t27 + t21 * t28 + t71 * t88) + m(6) * (-pkin(4) * t14 + qJ(5) * t13 + qJD(5) * t75) + m(7) * (t2 * t315 + t20 * t261 + t21 * t260 + t3 * t314) + t426 * t601 + qJD(5) * t173 - t131 * t132 - t88 * t58 + qJ(5) * t68 - pkin(4) * t66 + (t669 * t329 + (t537 + (m(6) * pkin(4) - m(7) * t385 + t665) * t366) * t377 + t630) * g(3) + (-m(6) * t454 - m(7) * (-t246 + t454) + t619) * g(2) + (-m(7) * (-t249 + t453) - m(6) * t453 + t620) * g(1) + t388 + (t547 - t636) * t86 + (-t548 - t637) * t85 + (-Ifges(5,2) * t394 - t199 + t639) * t581 + (t261 - t27) * t103 + (-t204 * t660 + t394 * t658) * t567 + (-t204 * t662 + t111 + t198 - t541) * t579 + (t260 - t28) * t102; t379 * t18 + t374 * t19 + (t132 - t58) * t394 + t427 * qJD(6) + (-t173 - t427) * t335 + t66 + (t2 * t374 - t394 * t71 + t3 * t379 + t402 + t321 * (-t20 * t374 + t21 * t379)) * m(7) + (-t335 * t75 + t394 * t94 + t14 + t402) * m(6); t55 * t586 - t20 * t102 + t21 * t103 + g(1) * t507 - g(2) * t690 - g(3) * t635 + t419 + (t56 - t692) * t589 + t686;];
tau  = t1;
