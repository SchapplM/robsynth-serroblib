% Calculate vector of inverse dynamics joint torques for
% S6RRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR6_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR6_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR6_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR6_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR6_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:18:42
% EndTime: 2019-03-09 22:20:19
% DurationCPUTime: 57.40s
% Computational Cost: add. (31291->1059), mult. (70182->1410), div. (0->0), fcn. (51800->18), ass. (0->469)
t438 = sin(qJ(2));
t443 = cos(qJ(2));
t487 = pkin(2) * t438 - pkin(8) * t443;
t375 = t487 * qJD(1);
t442 = cos(qJ(3));
t437 = sin(qJ(3));
t526 = qJD(1) * t438;
t506 = t437 * t526;
t288 = pkin(7) * t506 + t442 * t375;
t536 = t442 * t443;
t466 = pkin(3) * t438 - pkin(9) * t536;
t445 = -pkin(9) - pkin(8);
t507 = qJD(3) * t445;
t720 = -qJD(1) * t466 + t442 * t507 - t288;
t346 = t437 * t375;
t539 = t438 * t442;
t541 = t437 * t443;
t719 = t346 + (-pkin(7) * t539 - pkin(9) * t541) * qJD(1) - t437 * t507;
t436 = sin(qJ(4));
t441 = cos(qJ(4));
t370 = t436 * t442 + t437 * t441;
t656 = qJD(3) + qJD(4);
t275 = t656 * t370;
t460 = t370 * t443;
t308 = qJD(1) * t460;
t718 = t275 - t308;
t392 = t445 * t437;
t393 = t445 * t442;
t517 = qJD(4) * t441;
t518 = qJD(4) * t436;
t665 = t392 * t517 + t393 * t518 + t720 * t436 - t719 * t441;
t283 = t436 * t392 - t441 * t393;
t664 = -qJD(4) * t283 + t719 * t436 + t720 * t441;
t467 = t436 * t437 - t441 * t442;
t274 = t656 * t467;
t461 = t443 * t467;
t309 = qJD(1) * t461;
t717 = -pkin(4) * t526 - qJD(5) * t370 + t664 + (t274 - t309) * qJ(5);
t716 = -t718 * qJ(5) - qJD(5) * t467 + t665;
t433 = sin(pkin(11));
t434 = cos(pkin(11));
t202 = t274 * t433 - t275 * t434;
t221 = -t308 * t434 + t309 * t433;
t715 = t202 - t221;
t523 = qJD(2) * t442;
t367 = -t506 + t523;
t504 = t442 * t526;
t368 = qJD(2) * t437 + t504;
t266 = t367 * t436 + t368 * t441;
t489 = t441 * t367 - t368 * t436;
t196 = t266 * t434 + t433 * t489;
t435 = sin(qJ(6));
t440 = cos(qJ(6));
t687 = -t266 * t433 + t434 * t489;
t115 = t196 * t440 + t435 * t687;
t418 = pkin(7) * t526;
t390 = -qJD(2) * pkin(2) + t418;
t296 = -pkin(3) * t367 + t390;
t210 = -pkin(4) * t489 + qJD(5) + t296;
t130 = -pkin(5) * t687 + t210;
t688 = -t196 * t435 + t440 * t687;
t516 = qJD(1) * qJD(2);
t380 = qJDD(1) * t438 + t443 * t516;
t254 = qJD(3) * t367 + qJDD(2) * t437 + t380 * t442;
t255 = -qJD(3) * t368 + qJDD(2) * t442 - t380 * t437;
t144 = qJD(4) * t489 + t254 * t441 + t255 * t436;
t145 = -qJD(4) * t266 - t254 * t436 + t255 * t441;
t86 = -t144 * t433 + t145 * t434;
t87 = t144 * t434 + t145 * t433;
t29 = qJD(6) * t688 + t435 * t86 + t440 * t87;
t30 = -qJD(6) * t115 - t435 * t87 + t440 * t86;
t379 = qJDD(1) * t443 - t438 * t516;
t364 = qJDD(3) - t379;
t353 = qJDD(4) + t364;
t337 = qJDD(6) + t353;
t513 = Ifges(7,5) * t29 + Ifges(7,6) * t30 + Ifges(7,3) * t337;
t525 = qJD(1) * t443;
t404 = qJD(3) - t525;
t394 = qJD(4) + t404;
t488 = pkin(2) * t443 + pkin(8) * t438;
t384 = -pkin(1) - t488;
t354 = t384 * qJD(1);
t419 = pkin(7) * t525;
t391 = qJD(2) * pkin(8) + t419;
t269 = t442 * t354 - t391 * t437;
t228 = -pkin(9) * t368 + t269;
t216 = pkin(3) * t404 + t228;
t270 = t354 * t437 + t391 * t442;
t229 = pkin(9) * t367 + t270;
t223 = t436 * t229;
t153 = t441 * t216 - t223;
t690 = qJ(5) * t266;
t124 = t153 - t690;
t116 = pkin(4) * t394 + t124;
t225 = t441 * t229;
t154 = t216 * t436 + t225;
t671 = qJ(5) * t489;
t125 = t154 + t671;
t118 = t433 * t125;
t61 = t434 * t116 - t118;
t691 = pkin(10) * t196;
t46 = pkin(5) * t394 + t61 - t691;
t545 = t434 * t125;
t62 = t433 * t116 + t545;
t684 = pkin(10) * t687;
t50 = t62 + t684;
t16 = -t435 * t50 + t440 * t46;
t552 = qJDD(1) * pkin(1);
t273 = -pkin(2) * t379 - pkin(8) * t380 - t552;
t362 = t379 * pkin(7);
t333 = qJDD(2) * pkin(8) + t362;
t177 = -qJD(3) * t270 + t442 * t273 - t333 * t437;
t126 = pkin(3) * t364 - pkin(9) * t254 + t177;
t519 = qJD(3) * t442;
t521 = qJD(3) * t437;
t176 = t437 * t273 + t442 * t333 + t354 * t519 - t391 * t521;
t138 = pkin(9) * t255 + t176;
t49 = -qJD(4) * t154 + t441 * t126 - t138 * t436;
t35 = pkin(4) * t353 - qJ(5) * t144 - qJD(5) * t266 + t49;
t48 = t436 * t126 + t441 * t138 + t216 * t517 - t229 * t518;
t39 = qJ(5) * t145 + qJD(5) * t489 + t48;
t12 = t434 * t35 - t39 * t433;
t6 = pkin(5) * t353 - pkin(10) * t87 + t12;
t13 = t433 * t35 + t434 * t39;
t7 = pkin(10) * t86 + t13;
t2 = qJD(6) * t16 + t435 * t6 + t440 * t7;
t17 = t435 * t46 + t440 * t50;
t3 = -qJD(6) * t17 - t435 * t7 + t440 * t6;
t696 = -t3 * mrSges(7,1) + t2 * mrSges(7,2);
t465 = t513 - t696;
t583 = mrSges(7,3) * t17;
t584 = mrSges(7,3) * t16;
t589 = -t394 / 0.2e1;
t385 = qJD(6) + t394;
t591 = -t385 / 0.2e1;
t604 = -t196 / 0.2e1;
t606 = -t687 / 0.2e1;
t615 = -t115 / 0.2e1;
t617 = -t688 / 0.2e1;
t106 = Ifges(6,1) * t196 + Ifges(6,4) * t687 + Ifges(6,5) * t394;
t619 = -t106 / 0.2e1;
t105 = Ifges(6,4) * t196 + Ifges(6,2) * t687 + Ifges(6,6) * t394;
t621 = -t105 / 0.2e1;
t107 = Ifges(7,4) * t688;
t59 = Ifges(7,1) * t115 + Ifges(7,5) * t385 + t107;
t627 = -t59 / 0.2e1;
t561 = Ifges(7,4) * t115;
t58 = Ifges(7,2) * t688 + Ifges(7,6) * t385 + t561;
t629 = -t58 / 0.2e1;
t649 = mrSges(6,2) * t210 - mrSges(6,3) * t61;
t651 = -mrSges(6,1) * t210 + mrSges(6,3) * t62;
t683 = Ifges(5,3) + Ifges(6,3);
t659 = Ifges(5,5) * t144 + Ifges(6,5) * t87 + Ifges(5,6) * t145 + Ifges(6,6) * t86 + t353 * t683;
t693 = -t49 * mrSges(5,1) - t12 * mrSges(6,1) + t48 * mrSges(5,2) + t13 * mrSges(6,2);
t714 = -(Ifges(6,4) * t604 + Ifges(6,2) * t606 + Ifges(6,6) * t589 + t621 - t651) * t196 + (Ifges(6,1) * t604 + Ifges(6,4) * t606 + Ifges(6,5) * t589 + t619 - t649) * t687 + t465 + t659 - t693 + (-mrSges(7,2) * t130 + Ifges(7,1) * t615 + Ifges(7,4) * t617 + Ifges(7,5) * t591 + t584 + t627) * t688 - (mrSges(7,1) * t130 + Ifges(7,4) * t615 + Ifges(7,2) * t617 + Ifges(7,6) * t591 - t583 + t629) * t115;
t713 = pkin(5) * t196;
t680 = -t716 * t433 + t434 * t717;
t679 = t433 * t717 + t716 * t434;
t162 = -t228 * t436 - t225;
t131 = t162 - t671;
t163 = t441 * t228 - t223;
t132 = t163 - t690;
t544 = t434 * t436;
t560 = pkin(3) * qJD(4);
t674 = -t434 * t131 + t132 * t433 + (-t433 * t441 - t544) * t560;
t546 = t433 * t436;
t673 = -t433 * t131 - t434 * t132 + (t434 * t441 - t546) * t560;
t710 = mrSges(4,3) + mrSges(5,3) + mrSges(6,3) + mrSges(7,3);
t708 = m(4) * pkin(7);
t707 = t380 / 0.2e1;
t203 = -t274 * t434 - t275 * t433;
t222 = -t308 * t433 - t309 * t434;
t706 = -pkin(5) * t526 + t680 + (-t203 + t222) * pkin(10);
t705 = pkin(10) * t715 + t679;
t704 = t684 + t674;
t703 = t673 + t691;
t505 = t437 * t525;
t660 = -t419 + (-t505 + t521) * pkin(3);
t484 = mrSges(3,1) * t438 + mrSges(3,2) * t443;
t566 = Ifges(3,4) * t438;
t700 = -pkin(1) * t484 + t438 * (Ifges(3,1) * t443 - t566) / 0.2e1;
t432 = qJ(3) + qJ(4);
t424 = cos(t432);
t411 = pkin(4) * t424;
t428 = t442 * pkin(3);
t383 = t411 + t428;
t422 = pkin(11) + t432;
t408 = cos(t422);
t400 = pkin(5) * t408;
t321 = t400 + t383;
t314 = pkin(2) + t321;
t374 = pkin(2) + t383;
t412 = qJ(6) + t422;
t401 = sin(t412);
t402 = cos(t412);
t407 = sin(t422);
t414 = t428 + pkin(2);
t423 = sin(t432);
t483 = -mrSges(4,1) * t442 + mrSges(4,2) * t437;
t698 = m(4) * pkin(2) + m(5) * t414 + m(6) * t374 + m(7) * t314 + mrSges(5,1) * t424 + mrSges(6,1) * t408 + mrSges(7,1) * t402 - mrSges(5,2) * t423 - mrSges(6,2) * t407 - mrSges(7,2) * t401 - t483;
t431 = -qJ(5) + t445;
t425 = -pkin(10) + t431;
t697 = -m(4) * pkin(8) + m(5) * t445 + m(6) * t431 + m(7) * t425 - t710;
t695 = -t177 * mrSges(4,1) + t176 * mrSges(4,2);
t694 = -t269 * mrSges(4,1) + t270 * mrSges(4,2);
t476 = t443 * Ifges(3,2) + t566;
t692 = t154 * mrSges(5,2) + t17 * mrSges(7,2) + t62 * mrSges(6,2) + Ifges(3,6) * qJD(2) / 0.2e1 + qJD(1) * t476 / 0.2e1 - t153 * mrSges(5,1) - t16 * mrSges(7,1) - t61 * mrSges(6,1);
t579 = pkin(4) * t266;
t663 = t718 * pkin(4) + t660;
t522 = qJD(2) * t443;
t454 = t437 * t522 + t438 * t519;
t439 = sin(qJ(1));
t444 = cos(qJ(1));
t689 = g(1) * t444 + g(2) * t439;
t635 = m(5) * pkin(3);
t633 = t29 / 0.2e1;
t632 = t30 / 0.2e1;
t623 = t86 / 0.2e1;
t622 = t87 / 0.2e1;
t686 = m(6) + m(7);
t613 = t144 / 0.2e1;
t612 = t145 / 0.2e1;
t602 = t254 / 0.2e1;
t601 = t255 / 0.2e1;
t596 = t337 / 0.2e1;
t595 = t353 / 0.2e1;
t594 = t364 / 0.2e1;
t282 = t441 * t392 + t393 * t436;
t241 = -qJ(5) * t370 + t282;
t242 = -qJ(5) * t467 + t283;
t172 = t434 * t241 - t242 * t433;
t262 = t370 * t434 - t433 * t467;
t140 = -pkin(10) * t262 + t172;
t173 = t433 * t241 + t434 * t242;
t261 = -t370 * t433 - t434 * t467;
t141 = pkin(10) * t261 + t173;
t78 = t140 * t440 - t141 * t435;
t682 = qJD(6) * t78 + t435 * t706 + t440 * t705;
t79 = t140 * t435 + t141 * t440;
t681 = -qJD(6) * t79 - t435 * t705 + t440 * t706;
t580 = pkin(3) * t441;
t413 = pkin(4) + t580;
t329 = -pkin(3) * t546 + t434 * t413;
t322 = pkin(5) + t329;
t331 = pkin(3) * t544 + t413 * t433;
t237 = t322 * t440 - t331 * t435;
t678 = qJD(6) * t237 + t435 * t704 + t440 * t703;
t238 = t322 * t435 + t331 * t440;
t677 = -qJD(6) * t238 - t435 * t703 + t440 * t704;
t576 = pkin(4) * t434;
t409 = pkin(5) + t576;
t577 = pkin(4) * t433;
t330 = t409 * t440 - t435 * t577;
t67 = -t124 * t433 - t545;
t51 = t67 - t684;
t68 = t434 * t124 - t118;
t52 = t68 - t691;
t676 = t330 * qJD(6) - t435 * t51 - t440 * t52;
t332 = t409 * t435 + t440 * t577;
t675 = -t332 * qJD(6) + t435 * t52 - t440 * t51;
t672 = t635 + mrSges(4,1);
t328 = t467 * t438;
t666 = -pkin(5) * t715 + t663;
t366 = t442 * t384;
t268 = -pkin(9) * t539 + t366 + (-pkin(7) * t437 - pkin(3)) * t443;
t406 = pkin(7) * t536;
t303 = t437 * t384 + t406;
t543 = t437 * t438;
t277 = -pkin(9) * t543 + t303;
t205 = t436 * t268 + t441 * t277;
t360 = Ifges(4,4) * t367;
t245 = t368 * Ifges(4,1) + t404 * Ifges(4,5) + t360;
t417 = Ifges(3,4) * t525;
t662 = Ifges(3,1) * t526 + Ifges(3,5) * qJD(2) + t442 * t245 + t417;
t661 = qJD(2) * mrSges(3,1) + mrSges(4,1) * t367 - mrSges(4,2) * t368 - mrSges(3,3) * t526;
t363 = t380 * pkin(7);
t658 = t362 * t443 + t363 * t438;
t657 = t176 * t442 - t177 * t437;
t479 = -mrSges(7,1) * t401 - mrSges(7,2) * t402;
t655 = mrSges(5,1) * t423 + mrSges(6,1) * t407 + mrSges(5,2) * t424 + mrSges(6,2) * t408 - t479;
t535 = t443 * t444;
t300 = -t407 * t535 + t408 * t439;
t301 = t407 * t439 + t408 * t535;
t318 = -t423 * t535 + t424 * t439;
t319 = t423 * t439 + t424 * t535;
t286 = -t401 * t535 + t402 * t439;
t287 = t401 * t439 + t402 * t535;
t533 = t286 * mrSges(7,1) - t287 * mrSges(7,2);
t654 = -t318 * mrSges(5,1) - t300 * mrSges(6,1) + t319 * mrSges(5,2) + t301 * mrSges(6,2) - t533;
t538 = t439 * t443;
t298 = t407 * t538 + t408 * t444;
t299 = t407 * t444 - t408 * t538;
t316 = t423 * t538 + t424 * t444;
t317 = t423 * t444 - t424 * t538;
t284 = t401 * t538 + t402 * t444;
t285 = t401 * t444 - t402 * t538;
t534 = -t284 * mrSges(7,1) + t285 * mrSges(7,2);
t653 = t316 * mrSges(5,1) + t298 * mrSges(6,1) - t317 * mrSges(5,2) - t299 * mrSges(6,2) - t534;
t652 = -mrSges(5,1) * t296 + mrSges(5,3) * t154;
t650 = mrSges(5,2) * t296 - mrSges(5,3) * t153;
t648 = -m(3) - m(5) - m(4) - t686;
t647 = t368 * Ifges(4,5) + t266 * Ifges(5,5) + t196 * Ifges(6,5) + t115 * Ifges(7,5) + t367 * Ifges(4,6) + Ifges(5,6) * t489 + Ifges(6,6) * t687 + t688 * Ifges(7,6) + t404 * Ifges(4,3) + t385 * Ifges(7,3) + t394 * t683;
t578 = pkin(4) * t423;
t358 = -pkin(5) * t407 - t578;
t581 = pkin(3) * t437;
t320 = -t358 + t581;
t382 = t578 + t581;
t646 = -m(6) * t382 - m(7) * t320;
t485 = mrSges(3,1) * t443 - mrSges(3,2) * t438;
t644 = t438 * t710 + mrSges(2,1) + t485;
t643 = mrSges(2,2) - mrSges(3,3) + t646;
t637 = Ifges(7,4) * t633 + Ifges(7,2) * t632 + Ifges(7,6) * t596;
t636 = Ifges(7,1) * t633 + Ifges(7,4) * t632 + Ifges(7,5) * t596;
t634 = m(6) * pkin(4);
t631 = Ifges(6,4) * t622 + Ifges(6,2) * t623 + Ifges(6,6) * t595;
t630 = Ifges(6,1) * t622 + Ifges(6,4) * t623 + Ifges(6,5) * t595;
t628 = t58 / 0.2e1;
t626 = t59 / 0.2e1;
t625 = Ifges(5,4) * t613 + Ifges(5,2) * t612 + Ifges(5,6) * t595;
t624 = Ifges(5,1) * t613 + Ifges(5,4) * t612 + Ifges(5,5) * t595;
t620 = t105 / 0.2e1;
t618 = t106 / 0.2e1;
t616 = t688 / 0.2e1;
t614 = t115 / 0.2e1;
t611 = Ifges(4,1) * t602 + Ifges(4,4) * t601 + Ifges(4,5) * t594;
t562 = Ifges(5,4) * t266;
t184 = Ifges(5,2) * t489 + Ifges(5,6) * t394 + t562;
t610 = -t184 / 0.2e1;
t609 = t184 / 0.2e1;
t258 = Ifges(5,4) * t489;
t185 = Ifges(5,1) * t266 + Ifges(5,5) * t394 + t258;
t608 = -t185 / 0.2e1;
t607 = t185 / 0.2e1;
t605 = t687 / 0.2e1;
t603 = t196 / 0.2e1;
t600 = -t489 / 0.2e1;
t599 = t489 / 0.2e1;
t598 = -t266 / 0.2e1;
t597 = t266 / 0.2e1;
t592 = t368 / 0.2e1;
t590 = t385 / 0.2e1;
t588 = t394 / 0.2e1;
t582 = pkin(3) * t368;
t573 = g(3) * t438;
t426 = t438 * pkin(7);
t208 = -qJD(2) * t461 - t275 * t438;
t524 = qJD(2) * t438;
t378 = t487 * qJD(2);
t511 = pkin(7) * t524;
t528 = t442 * t378 + t437 * t511;
t201 = t466 * qJD(2) + (-t406 + (pkin(9) * t438 - t384) * t437) * qJD(3) + t528;
t226 = t437 * t378 + t384 * t519 + (-t438 * t523 - t443 * t521) * pkin(7);
t207 = -pkin(9) * t454 + t226;
t98 = -qJD(4) * t205 + t441 * t201 - t207 * t436;
t65 = pkin(4) * t524 - qJ(5) * t208 + qJD(5) * t328 + t98;
t209 = -qJD(2) * t460 + t328 * t656;
t327 = t370 * t438;
t97 = t436 * t201 + t441 * t207 + t268 * t517 - t277 * t518;
t69 = qJ(5) * t209 - qJD(5) * t327 + t97;
t32 = t433 * t65 + t434 * t69;
t568 = mrSges(5,3) * t489;
t567 = mrSges(5,3) * t266;
t565 = Ifges(3,4) * t443;
t564 = Ifges(4,4) * t437;
t563 = Ifges(4,4) * t442;
t559 = t269 * mrSges(4,3);
t558 = t270 * mrSges(4,3);
t557 = t368 * Ifges(4,4);
t542 = t437 * t439;
t540 = t437 * t444;
t204 = t441 * t268 - t277 * t436;
t168 = -pkin(4) * t443 + qJ(5) * t328 + t204;
t179 = -qJ(5) * t327 + t205;
t100 = t433 * t168 + t434 * t179;
t381 = pkin(3) * t543 + t426;
t520 = qJD(3) * t438;
t421 = pkin(7) * t522;
t508 = Ifges(4,5) * t254 + Ifges(4,6) * t255 + Ifges(4,3) * t364;
t297 = pkin(3) * t454 + t421;
t244 = t367 * Ifges(4,2) + t404 * Ifges(4,6) + t557;
t501 = -t437 * t244 / 0.2e1;
t42 = -t86 * mrSges(6,1) + t87 * mrSges(6,2);
t10 = -t30 * mrSges(7,1) + t29 * mrSges(7,2);
t31 = -t433 * t69 + t434 * t65;
t99 = t434 * t168 - t179 * t433;
t312 = pkin(4) * t467 - t414;
t271 = pkin(4) * t327 + t381;
t220 = t582 + t579;
t334 = -qJDD(2) * pkin(2) + t363;
t482 = mrSges(4,1) * t437 + mrSges(4,2) * t442;
t478 = Ifges(4,1) * t442 - t564;
t477 = Ifges(4,1) * t437 + t563;
t475 = -Ifges(4,2) * t437 + t563;
t474 = Ifges(4,2) * t442 + t564;
t473 = Ifges(3,5) * t443 - Ifges(3,6) * t438;
t472 = Ifges(4,5) * t442 - Ifges(4,6) * t437;
t471 = Ifges(4,5) * t437 + Ifges(4,6) * t442;
t233 = -t327 * t433 - t328 * t434;
t85 = -pkin(5) * t443 - pkin(10) * t233 + t99;
t232 = -t327 * t434 + t328 * t433;
t89 = pkin(10) * t232 + t100;
t43 = -t435 * t89 + t440 * t85;
t44 = t435 * t85 + t440 * t89;
t164 = t232 * t440 - t233 * t435;
t165 = t232 * t435 + t233 * t440;
t193 = t261 * t440 - t262 * t435;
t197 = t261 * t435 + t262 * t440;
t470 = t314 * t443 - t425 * t438;
t469 = t374 * t443 - t431 * t438;
t468 = t443 * t414 - t438 * t445;
t182 = -pkin(4) * t209 + t297;
t343 = -t437 * t535 + t439 * t442;
t341 = t437 * t538 + t442 * t444;
t463 = t390 * t482;
t217 = -pkin(3) * t255 + t334;
t455 = -t437 * t520 + t442 * t522;
t453 = Ifges(4,5) * t438 + t443 * t478;
t452 = Ifges(4,6) * t438 + t443 * t475;
t451 = Ifges(4,3) * t438 + t443 * t472;
t101 = -pkin(4) * t145 + qJDD(5) + t217;
t387 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t525;
t359 = t400 + t411;
t355 = t482 * t438;
t344 = t442 * t535 + t542;
t342 = -t439 * t536 + t540;
t302 = -pkin(7) * t541 + t366;
t295 = mrSges(4,1) * t404 - mrSges(4,3) * t368;
t294 = -mrSges(4,2) * t404 + mrSges(4,3) * t367;
t289 = -pkin(7) * t504 + t346;
t231 = mrSges(5,1) * t394 - t567;
t230 = -mrSges(5,2) * t394 + t568;
t227 = -qJD(3) * t303 + t528;
t219 = -mrSges(4,2) * t364 + mrSges(4,3) * t255;
t218 = mrSges(4,1) * t364 - mrSges(4,3) * t254;
t215 = -pkin(5) * t261 + t312;
t200 = -mrSges(5,1) * t489 + mrSges(5,2) * t266;
t189 = -mrSges(4,1) * t255 + mrSges(4,2) * t254;
t188 = -pkin(5) * t232 + t271;
t181 = mrSges(6,1) * t394 - mrSges(6,3) * t196;
t180 = -mrSges(6,2) * t394 + mrSges(6,3) * t687;
t169 = t254 * Ifges(4,4) + t255 * Ifges(4,2) + t364 * Ifges(4,6);
t156 = t221 * t435 + t222 * t440;
t155 = t221 * t440 - t222 * t435;
t146 = t579 + t713;
t139 = t220 + t713;
t136 = t208 * t434 + t209 * t433;
t135 = -t208 * t433 + t209 * t434;
t134 = -mrSges(5,2) * t353 + mrSges(5,3) * t145;
t133 = mrSges(5,1) * t353 - mrSges(5,3) * t144;
t117 = -mrSges(6,1) * t687 + mrSges(6,2) * t196;
t103 = mrSges(7,1) * t385 - mrSges(7,3) * t115;
t102 = -mrSges(7,2) * t385 + mrSges(7,3) * t688;
t94 = -pkin(5) * t135 + t182;
t91 = -qJD(6) * t197 + t202 * t440 - t203 * t435;
t90 = qJD(6) * t193 + t202 * t435 + t203 * t440;
t88 = -mrSges(5,1) * t145 + mrSges(5,2) * t144;
t77 = mrSges(6,1) * t353 - mrSges(6,3) * t87;
t76 = -mrSges(6,2) * t353 + mrSges(6,3) * t86;
t60 = -mrSges(7,1) * t688 + mrSges(7,2) * t115;
t54 = -qJD(6) * t165 + t135 * t440 - t136 * t435;
t53 = qJD(6) * t164 + t135 * t435 + t136 * t440;
t45 = -pkin(5) * t86 + t101;
t25 = -mrSges(7,2) * t337 + mrSges(7,3) * t30;
t24 = mrSges(7,1) * t337 - mrSges(7,3) * t29;
t21 = pkin(10) * t135 + t32;
t20 = pkin(5) * t524 - pkin(10) * t136 + t31;
t5 = -qJD(6) * t44 + t20 * t440 - t21 * t435;
t4 = qJD(6) * t43 + t20 * t435 + t21 * t440;
t1 = [(-t16 * t53 + t164 * t2 - t165 * t3 + t17 * t54) * mrSges(7,3) + (-qJDD(2) * mrSges(3,1) + mrSges(3,3) * t380 + t189) * t426 + m(4) * (t176 * t303 + t177 * t302 + t226 * t270 + t227 * t269) + (Ifges(6,1) * t136 + Ifges(6,4) * t135) * t603 + (Ifges(6,1) * t233 + Ifges(6,4) * t232) * t622 + (-Ifges(5,4) * t328 - Ifges(5,2) * t327) * t612 + (Ifges(5,4) * t208 + Ifges(5,2) * t209) * t599 + (Ifges(7,4) * t165 + Ifges(7,2) * t164) * t632 + (Ifges(7,4) * t53 + Ifges(7,2) * t54) * t616 + (Ifges(5,5) * t208 + Ifges(6,5) * t136 + Ifges(5,6) * t209 + Ifges(6,6) * t135) * t588 + (Ifges(6,4) * t136 + Ifges(6,2) * t135) * t605 + (Ifges(6,4) * t233 + Ifges(6,2) * t232) * t623 + t658 * mrSges(3,3) + (-t176 * t543 - t177 * t539 - t269 * t455 - t270 * t454) * mrSges(4,3) + (t566 + t476) * t379 / 0.2e1 + (-t12 * t233 + t13 * t232 + t135 * t62 - t136 * t61) * mrSges(6,3) + t700 * t516 + qJD(2) ^ 2 * t473 / 0.2e1 - t387 * t511 + t565 * t707 - t661 * t421 + m(6) * (t100 * t13 + t101 * t271 + t12 * t99 + t182 * t210 + t31 * t61 + t32 * t62) + m(7) * (t130 * t94 + t16 * t5 + t17 * t4 + t188 * t45 + t2 * t44 + t3 * t43) + m(5) * (t153 * t98 + t154 * t97 + t204 * t49 + t205 * t48 + t217 * t381 + t296 * t297) + (Ifges(7,5) * t53 + Ifges(7,6) * t54) * t590 + (Ifges(7,5) * t165 + Ifges(7,6) * t164) * t596 + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t658) + t485 * t552 + t390 * (mrSges(4,1) * t454 + mrSges(4,2) * t455) + (Ifges(7,1) * t165 + Ifges(7,4) * t164) * t633 + (Ifges(7,1) * t53 + Ifges(7,4) * t54) * t614 - t169 * t543 / 0.2e1 + t217 * (mrSges(5,1) * t327 - mrSges(5,2) * t328) + (-Ifges(5,1) * t328 - Ifges(5,4) * t327) * t613 + (Ifges(5,1) * t208 + Ifges(5,4) * t209) * t597 + (-t153 * t208 + t154 * t209 - t327 * t48 + t328 * t49) * mrSges(5,3) + (Ifges(7,6) * t616 - t692 + t683 * t588 + Ifges(7,3) * t590 + Ifges(5,5) * t597 + Ifges(5,6) * t599 + Ifges(6,5) * t603 + Ifges(6,6) * t605 + Ifges(7,5) * t614 + t647 / 0.2e1 - t694) * t524 + (-Ifges(5,5) * t328 + Ifges(6,5) * t233 - Ifges(5,6) * t327 + Ifges(6,6) * t232) * t595 + t165 * t636 + t164 * t637 - t328 * t624 - t327 * t625 + t53 * t626 + t54 * t628 + t233 * t630 + t232 * t631 + Ifges(2,3) * qJDD(1) - pkin(1) * (-mrSges(3,1) * t379 + mrSges(3,2) * t380) + t381 * t88 + (-t540 * t635 - t342 * mrSges(4,1) - t317 * mrSges(5,1) - t299 * mrSges(6,1) - t285 * mrSges(7,1) - t341 * mrSges(4,2) - t316 * mrSges(5,2) - t298 * mrSges(6,2) - t284 * mrSges(7,2) + (-m(7) * (-pkin(1) - t470) - m(6) * (-pkin(1) - t469) - m(5) * (-pkin(1) - t468) - m(4) * t384 + m(3) * pkin(1) + t644) * t439 + (pkin(7) * t648 + t643) * t444) * g(1) + t334 * t355 + t297 * t200 + t302 * t218 + t303 * t219 + t226 * t294 + t227 * t295 + t296 * (-mrSges(5,1) * t209 + mrSges(5,2) * t208) + t367 * (qJD(2) * t452 - t474 * t520) / 0.2e1 + t404 * (qJD(2) * t451 - t471 * t520) / 0.2e1 + t271 * t42 + t101 * (-mrSges(6,1) * t232 + mrSges(6,2) * t233) + t97 * t230 + t98 * t231 + t210 * (-mrSges(6,1) * t135 + mrSges(6,2) * t136) + t204 * t133 + t205 * t134 + t188 * t10 + t32 * t180 + t31 * t181 + t182 * t117 + t45 * (-mrSges(7,1) * t164 + mrSges(7,2) * t165) + t130 * (-mrSges(7,1) * t54 + mrSges(7,2) * t53) + t43 * t24 + t44 * t25 + (-t542 * t635 - t344 * mrSges(4,1) - t319 * mrSges(5,1) - t301 * mrSges(6,1) - t287 * mrSges(7,1) - t343 * mrSges(4,2) - t318 * mrSges(5,2) - t300 * mrSges(6,2) - t286 * mrSges(7,2) + t648 * (t444 * pkin(1) + t439 * pkin(7)) + t643 * t439 + (-m(4) * t488 - m(5) * t468 - m(6) * t469 - m(7) * t470 - t644) * t444) * g(2) - (t442 * t244 + t437 * t245) * t520 / 0.2e1 + (qJD(2) * t453 - t477 * t520) * t592 + t208 * t607 + t209 * t609 + t539 * t611 + t136 * t618 + t135 * t620 + ((pkin(7) * mrSges(3,3) + Ifges(3,2) / 0.2e1) * t379 - Ifges(5,6) * t612 - Ifges(5,5) * t613 + t693 + t696 - t659 / 0.2e1 - t513 / 0.2e1 - t508 / 0.2e1 - Ifges(7,5) * t633 - Ifges(7,6) * t632 + (-mrSges(3,2) * pkin(7) + Ifges(3,6)) * qJDD(2) - t683 * t595 + (-Ifges(3,2) * t438 + t565) * t516 / 0.2e1 - Ifges(4,3) * t594 - Ifges(7,3) * t596 - Ifges(4,6) * t601 - Ifges(4,5) * t602 - Ifges(6,5) * t622 - Ifges(6,6) * t623 + Ifges(3,4) * t707 + t695) * t443 + (t662 / 0.2e1 + t390 * t708 + t501) * t522 + (Ifges(3,1) * t380 + Ifges(3,5) * qJDD(2) + t334 * t708 + t472 * t594 + t475 * t601 + t478 * t602) * t438 + t94 * t60 + t99 * t77 + t100 * t76 + t4 * t102 + t5 * t103; -t485 * g(3) + t689 * t484 + (Ifges(7,4) * t156 + Ifges(7,2) * t155) * t617 + (Ifges(6,4) * t222 + Ifges(6,2) * t221) * t606 + (-pkin(2) * t334 - t269 * t288 - t270 * t289 - t390 * t419) * m(4) + (Ifges(7,4) * t614 + Ifges(7,2) * t616 + Ifges(7,6) * t590 + t583 + t628) * t91 + (Ifges(7,1) * t156 + Ifges(7,4) * t155) * t615 + (Ifges(6,1) * t222 + Ifges(6,4) * t221) * t604 + (Ifges(5,5) * t598 + Ifges(6,5) * t604 + Ifges(7,5) * t615 + Ifges(5,6) * t600 + Ifges(6,6) * t606 + Ifges(7,6) * t617 + Ifges(7,3) * t591 + t683 * t589 + t692) * t526 + ((t269 * t536 + t270 * t541) * qJD(1) + t657) * mrSges(4,3) + (-Ifges(5,4) * t309 - Ifges(5,2) * t308) * t600 + (-Ifges(5,1) * t309 - Ifges(5,4) * t308) * t598 + (-t153 * t309 + t154 * t308 - t370 * t49 - t467 * t48) * mrSges(5,3) + (-Ifges(5,5) * t309 + Ifges(6,5) * t222 - Ifges(5,6) * t308 + Ifges(6,6) * t221) * t589 - t521 * t558 + t244 * t505 / 0.2e1 - t700 * qJD(1) ^ 2 + (-t559 + t245 / 0.2e1) * t519 + (Ifges(6,1) * t603 + Ifges(6,4) * t605 + Ifges(6,5) * t588 + t618 + t649) * t203 - (Ifges(5,1) * t597 + Ifges(5,4) * t599 + Ifges(5,5) * t588 + t607 + t650) * t274 + (Ifges(6,4) * t603 + Ifges(6,2) * t605 + Ifges(6,6) * t588 + t620 + t651) * t202 - (Ifges(5,4) * t597 + Ifges(5,2) * t599 + Ifges(5,6) * t588 + t609 + t652) * t275 - t463 * t525 + t660 * t200 + t661 * t419 - (-Ifges(3,2) * t526 + t417 + t662) * t525 / 0.2e1 + t663 * t117 + t664 * t231 + (t153 * t664 + t154 * t665 - t217 * t414 + t282 * t49 + t283 * t48 + t296 * t660) * m(5) + t665 * t230 + t666 * t60 + (Ifges(7,1) * t614 + Ifges(7,4) * t616 + Ifges(7,5) * t590 - t584 + t626) * t90 + (g(3) * t697 + qJD(1) * t694 + t689 * t698) * t438 + (-g(3) * t698 + t689 * t697) * t443 + t387 * t418 + (m(4) * ((-t269 * t442 - t270 * t437) * qJD(3) + t657) - t295 * t519 - t294 * t521 - t437 * t218 + t442 * t219) * pkin(8) + (Ifges(5,5) * t370 + Ifges(6,5) * t262 - Ifges(5,6) * t467 + Ifges(6,6) * t261) * t595 + t217 * (mrSges(5,1) * t467 + mrSges(5,2) * t370) + (Ifges(5,4) * t370 - Ifges(5,2) * t467) * t612 + (Ifges(5,1) * t370 - Ifges(5,4) * t467) * t613 - t467 * t625 + (t463 + t501) * qJD(3) - t473 * t516 / 0.2e1 - t647 * t526 / 0.2e1 - t296 * (mrSges(5,1) * t308 - mrSges(5,2) * t309) + t334 * t483 + (-t155 * t17 + t156 * t16 + t193 * t2 - t197 * t3) * mrSges(7,3) + (Ifges(7,4) * t197 + Ifges(7,2) * t193) * t632 + (Ifges(7,1) * t197 + Ifges(7,4) * t193) * t633 + t197 * t636 + t193 * t637 + t370 * t624 + t156 * t627 + t155 * t629 + t262 * t630 + t261 * t631 + t442 * t169 / 0.2e1 - t414 * t88 + Ifges(3,5) * t380 + Ifges(3,6) * t379 - t362 * mrSges(3,2) - t363 * mrSges(3,1) + Ifges(3,3) * qJDD(2) + t312 * t42 - t289 * t294 - t288 * t295 + t282 * t133 + t283 * t134 + ((-t156 + t90) * mrSges(7,2) + (t155 - t91) * mrSges(7,1)) * t130 + t101 * (-mrSges(6,1) * t261 + mrSges(6,2) * t262) - t210 * (-mrSges(6,1) * t221 + mrSges(6,2) * t222) + t215 * t10 + t45 * (-mrSges(7,1) * t193 + mrSges(7,2) * t197) - pkin(2) * t189 + t172 * t77 + t173 * t76 + (Ifges(7,5) * t156 + Ifges(7,6) * t155) * t591 + t679 * t180 + t680 * t181 + (t101 * t312 + t12 * t172 + t13 * t173 + t210 * t663 + t61 * t680 + t62 * t679) * m(6) + t681 * t103 + t682 * t102 + (t130 * t666 + t16 * t681 + t17 * t682 + t2 * t79 + t215 * t45 + t3 * t78) * m(7) + (t367 * t475 + t368 * t478 + t404 * t472) * qJD(3) / 0.2e1 - (t367 * t452 + t368 * t453 + t404 * t451) * qJD(1) / 0.2e1 + t471 * t594 + (Ifges(7,5) * t197 + Ifges(7,6) * t193) * t596 + t474 * t601 + t477 * t602 - t309 * t608 - t308 * t610 + t437 * t611 + t222 * t619 + t221 * t621 + (Ifges(6,1) * t262 + Ifges(6,4) * t261) * t622 + (Ifges(6,4) * t262 + Ifges(6,2) * t261) * t623 + t78 * t24 + t79 * t25 + (-t12 * t262 + t13 * t261 - t221 * t62 + t222 * t61) * mrSges(6,3); (Ifges(5,1) * t598 + Ifges(5,4) * t600 + Ifges(5,5) * t589 + t608 - t650) * t489 + (t134 * t436 + t230 * t517 - t231 * t518) * pkin(3) + (-mrSges(4,2) * t342 - m(6) * (-t382 * t538 - t383 * t444) - m(7) * (-t320 * t538 - t321 * t444) + t672 * t341 + t653) * g(2) + (-m(7) * (-t320 * t535 + t321 * t439) + mrSges(4,2) * t344 - m(6) * (-t382 * t535 + t383 * t439) - t672 * t343 + t654) * g(1) + t673 * t180 + (t12 * t329 + t13 * t331 - t210 * t220 + t61 * t674 + t62 * t673) * m(6) + t674 * t181 + t677 * t103 + (-t130 * t139 + t16 * t677 + t17 * t678 + t2 * t238 + t237 * t3) * m(7) + t678 * t102 - (Ifges(5,4) * t598 + Ifges(5,2) * t600 + Ifges(5,6) * t589 + t610 - t652) * t266 + (m(5) * t581 - t646 + t655) * t573 - t368 * (Ifges(4,1) * t367 - t557) / 0.2e1 - t200 * t582 - m(5) * (t153 * t162 + t154 * t163 + t296 * t582) + t508 + t714 + (t436 * t48 + t441 * t49 + (-t153 * t436 + t154 * t441) * qJD(4)) * t635 - t404 * (Ifges(4,5) * t367 - Ifges(4,6) * t368) / 0.2e1 - t390 * (mrSges(4,1) * t368 + mrSges(4,2) * t367) + g(3) * t355 + t329 * t77 + t331 * t76 - t269 * t294 + t270 * t295 + t237 * t24 + t238 * t25 - t163 * t230 - t162 * t231 - t220 * t117 - t139 * t60 - t695 - (-Ifges(4,2) * t368 + t245 + t360) * t367 / 0.2e1 + t368 * t558 + t367 * t559 + t133 * t580 + t244 * t592; (-t230 + t568) * t153 + t675 * t103 + (-t130 * t146 + t16 * t675 + t17 * t676 + t2 * t332 + t3 * t330) * m(7) + t676 * t102 + (-Ifges(5,2) * t266 + t185 + t258) * t600 + (-m(7) * (t358 * t538 - t359 * t444) + t316 * t634 + t653) * g(2) + (-m(7) * (t358 * t535 + t359 * t439) - t318 * t634 + t654) * g(1) + (m(6) * t578 - m(7) * t358 + t655) * t573 - t117 * t579 - m(6) * (t210 * t579 + t61 * t67 + t62 * t68) + t714 + (t231 + t567) * t154 + (t12 * t434 + t13 * t433) * t634 + t330 * t24 + t332 * t25 - t67 * t181 - t68 * t180 - t146 * t60 - t296 * (mrSges(5,1) * t266 + mrSges(5,2) * t489) + (Ifges(5,5) * t489 - Ifges(5,6) * t266) * t589 + (Ifges(5,1) * t489 - t562) * t598 + t77 * t576 + t76 * t577 + t184 * t597; -t688 * t102 + t115 * t103 - t687 * t180 + t196 * t181 + t10 + t42 + (t443 * g(3) - t438 * t689) * t686 + (t115 * t16 - t17 * t688 + t45) * m(7) + (t196 * t61 - t62 * t687 + t101) * m(6); -t130 * (mrSges(7,1) * t115 + mrSges(7,2) * t688) + (Ifges(7,1) * t688 - t561) * t615 + t58 * t614 + (Ifges(7,5) * t688 - Ifges(7,6) * t115) * t591 - t16 * t102 + t17 * t103 - g(1) * t533 - g(2) * t534 - t479 * t573 + (t115 * t17 + t16 * t688) * mrSges(7,3) + t465 + (-Ifges(7,2) * t115 + t107 + t59) * t617;];
tau  = t1;
