% Calculate vector of inverse dynamics joint torques for
% S6PRRRRR2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 00:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRRR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:42:28
% EndTime: 2019-03-09 00:43:15
% DurationCPUTime: 26.87s
% Computational Cost: add. (16513->839), mult. (36339->1144), div. (0->0), fcn. (27987->18), ass. (0->403)
t385 = sin(qJ(2));
t379 = sin(pkin(6));
t487 = qJD(1) * t379;
t453 = t385 * t487;
t384 = sin(qJ(3));
t481 = qJD(3) * t384;
t606 = pkin(3) * t481 - t453;
t646 = mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t383 = sin(qJ(4));
t388 = cos(qJ(4));
t389 = cos(qJ(3));
t325 = t383 * t384 - t388 * t389;
t404 = t325 * qJD(4);
t257 = -qJD(3) * t325 - t404;
t327 = t383 * t389 + t384 * t388;
t375 = qJD(3) + qJD(4);
t258 = t375 * t327;
t663 = pkin(4) * t258 - pkin(10) * t257 + t606;
t392 = -pkin(9) - pkin(8);
t455 = qJD(3) * t392;
t334 = t384 * t455;
t335 = t389 * t455;
t390 = cos(qJ(2));
t452 = t390 * t487;
t348 = t392 * t384;
t350 = t392 * t389;
t608 = t388 * t348 + t350 * t383;
t610 = qJD(4) * t608 + t325 * t452 + t334 * t388 + t335 * t383;
t376 = qJ(5) + qJ(6);
t369 = sin(t376);
t371 = cos(t376);
t382 = sin(qJ(5));
t387 = cos(qJ(5));
t425 = -mrSges(6,1) * t387 + mrSges(6,2) * t382;
t658 = -mrSges(7,1) * t371 + mrSges(7,2) * t369 - mrSges(5,1) + t425;
t366 = pkin(3) * t389 + pkin(2);
t244 = pkin(4) * t325 - pkin(10) * t327 - t366;
t276 = t348 * t383 - t350 * t388;
t476 = qJD(5) * t387;
t477 = qJD(5) * t382;
t617 = t244 * t476 - t276 * t477 + t382 * t663 + t387 * t610;
t478 = qJD(4) * t388;
t466 = pkin(3) * t478;
t338 = qJD(2) * pkin(8) + t453;
t434 = pkin(9) * qJD(2) + t338;
t380 = cos(pkin(6));
t486 = qJD(1) * t380;
t451 = t384 * t486;
t264 = t389 * t434 + t451;
t251 = t383 * t264;
t359 = t389 * t486;
t263 = -t384 * t434 + t359;
t149 = t263 * t388 - t251;
t482 = qJD(2) * t389;
t484 = qJD(2) * t384;
t318 = -t383 * t484 + t388 * t482;
t319 = t327 * qJD(2);
t242 = pkin(4) * t319 - pkin(10) * t318;
t469 = pkin(3) * t484;
t201 = t242 + t469;
t90 = t387 * t149 + t382 * t201;
t660 = t387 * t466 - t90;
t89 = -t149 * t382 + t387 * t201;
t659 = -t382 * t466 - t89;
t657 = -t382 * t610 + t387 * t663;
t364 = pkin(5) * t387 + pkin(4);
t377 = qJ(3) + qJ(4);
t370 = sin(t377);
t372 = cos(t377);
t391 = -pkin(11) - pkin(10);
t593 = m(6) * pkin(10);
t594 = m(6) * pkin(4);
t656 = (-m(7) * t364 - t594 + t658) * t372 + (m(7) * t391 - t593 + t646) * t370;
t262 = t387 * t276;
t525 = t257 * t387;
t655 = -pkin(11) * t525 + pkin(5) * t258 + (-t262 + (pkin(11) * t327 - t244) * t382) * qJD(5) + t657;
t448 = t327 * t476;
t408 = t257 * t382 + t448;
t654 = pkin(11) * t408 - t617;
t561 = pkin(3) * t383;
t363 = pkin(10) + t561;
t552 = -pkin(11) - t363;
t439 = qJD(5) * t552;
t522 = t318 * t382;
t471 = pkin(11) * t522;
t653 = t382 * t439 + t471 + t660;
t521 = t318 * t387;
t428 = t319 * pkin(5) - pkin(11) * t521;
t652 = t387 * t439 - t428 + t659;
t454 = qJD(5) * t391;
t253 = qJD(3) * pkin(3) + t263;
t144 = t253 * t388 - t251;
t95 = t387 * t144 + t382 * t242;
t651 = t382 * t454 + t471 - t95;
t94 = -t144 * t382 + t387 * t242;
t650 = t387 * t454 - t428 - t94;
t269 = -t319 * t382 + t375 * t387;
t270 = t319 * t387 + t375 * t382;
t547 = mrSges(5,3) * t319;
t609 = mrSges(5,1) * t375 + mrSges(6,1) * t269 - mrSges(6,2) * t270 - t547;
t381 = sin(qJ(6));
t386 = cos(qJ(6));
t162 = t269 * t381 + t270 * t386;
t432 = t386 * t269 - t270 * t381;
t88 = -mrSges(7,1) * t432 + mrSges(7,2) * t162;
t649 = t88 - t609;
t592 = m(7) * pkin(5);
t138 = -pkin(4) * t375 - t144;
t99 = -pkin(5) * t269 + t138;
t644 = mrSges(7,2) * t99;
t252 = t388 * t264;
t145 = t253 * t383 + t252;
t139 = pkin(10) * t375 + t145;
t307 = -qJD(2) * t366 - t452;
t183 = -pkin(4) * t318 - pkin(10) * t319 + t307;
t82 = -t139 * t382 + t387 * t183;
t643 = t82 * mrSges(6,1);
t83 = t139 * t387 + t183 * t382;
t642 = t83 * mrSges(6,2);
t309 = qJD(5) - t318;
t301 = qJD(6) + t309;
t627 = t309 * Ifges(6,3);
t628 = t269 * Ifges(6,6);
t641 = t270 * Ifges(6,5) + t162 * Ifges(7,5) + Ifges(7,6) * t432 + t301 * Ifges(7,3) + t627 + t628;
t148 = t263 * t383 + t252;
t479 = qJD(4) * t383;
t640 = -pkin(3) * t479 + t148;
t424 = t382 * mrSges(6,1) + mrSges(6,2) * t387;
t639 = -t369 * mrSges(7,1) - t371 * mrSges(7,2) - t382 * t592 - mrSges(5,3) - t424;
t513 = t379 * t385;
t314 = t380 * t389 - t384 * t513;
t528 = cos(pkin(12));
t436 = t528 * t390;
t378 = sin(pkin(12));
t515 = t378 * t385;
t313 = -t380 * t515 + t436;
t512 = t379 * t389;
t638 = -t313 * t384 + t378 * t512;
t374 = qJDD(3) + qJDD(4);
t282 = t338 * t389 + t451;
t485 = qJD(2) * t379;
t445 = qJD(1) * t485;
t352 = t390 * t445;
t473 = qJDD(1) * t379;
t303 = t385 * t473 + t352;
t294 = qJDD(2) * pkin(8) + t303;
t472 = qJDD(1) * t380;
t168 = -qJD(3) * t282 - t294 * t384 + t389 * t472;
t475 = qJD(2) * qJD(3);
t337 = qJDD(2) * t384 + t389 * t475;
t135 = qJDD(3) * pkin(3) - pkin(9) * t337 + t168;
t167 = qJD(3) * t359 + t389 * t294 - t338 * t481 + t384 * t472;
t336 = qJDD(2) * t389 - t384 * t475;
t140 = pkin(9) * t336 + t167;
t49 = t383 * t135 + t388 * t140 + t253 * t478 - t264 * t479;
t46 = pkin(10) * t374 + t49;
t187 = -qJD(2) * t404 + t336 * t383 + t337 * t388;
t188 = -qJD(4) * t319 + t336 * t388 - t383 * t337;
t351 = t385 * t445;
t302 = t390 * t473 - t351;
t293 = -qJDD(2) * pkin(2) - t302;
t243 = -pkin(3) * t336 + t293;
t80 = -pkin(4) * t188 - pkin(10) * t187 + t243;
t14 = -t139 * t477 + t183 * t476 + t382 * t80 + t387 * t46;
t15 = -qJD(5) * t83 - t382 * t46 + t387 * t80;
t637 = t14 * t387 - t15 * t382;
t437 = t528 * t385;
t514 = t378 * t390;
t311 = t380 * t437 + t514;
t438 = t379 * t528;
t245 = -t311 * t370 - t372 * t438;
t246 = t311 * t372 - t370 * t438;
t636 = t245 * t658 + t246 * t646;
t516 = t378 * t379;
t247 = -t313 * t370 + t372 * t516;
t248 = t313 * t372 + t370 * t516;
t635 = t247 * t658 + t248 * t646;
t295 = -t370 * t513 + t372 * t380;
t296 = t370 * t380 + t372 * t513;
t634 = t295 * t658 + t296 * t646;
t68 = pkin(11) * t269 + t83;
t529 = t386 * t68;
t67 = -pkin(11) * t270 + t82;
t58 = pkin(5) * t309 + t67;
t22 = t381 * t58 + t529;
t633 = -mrSges(7,1) * t99 + t22 * mrSges(7,3);
t186 = qJDD(5) - t188;
t182 = qJDD(6) + t186;
t576 = t182 / 0.2e1;
t113 = qJD(5) * t269 + t187 * t387 + t374 * t382;
t114 = -qJD(5) * t270 - t187 * t382 + t374 * t387;
t39 = -qJD(6) * t162 - t113 * t381 + t114 * t386;
t588 = t39 / 0.2e1;
t38 = qJD(6) * t432 + t113 * t386 + t114 * t381;
t589 = t38 / 0.2e1;
t590 = Ifges(7,1) * t589 + Ifges(7,4) * t588 + Ifges(7,5) * t576;
t591 = Ifges(7,4) * t589 + Ifges(7,2) * t588 + Ifges(7,6) * t576;
t582 = t113 / 0.2e1;
t581 = t114 / 0.2e1;
t575 = t186 / 0.2e1;
t632 = t336 / 0.2e1;
t631 = t337 / 0.2e1;
t146 = t387 * t244 - t276 * t382;
t519 = t327 * t387;
t102 = pkin(5) * t325 - pkin(11) * t519 + t146;
t147 = t382 * t244 + t262;
t520 = t327 * t382;
t122 = -pkin(11) * t520 + t147;
t56 = t102 * t381 + t122 * t386;
t630 = -qJD(6) * t56 + t654 * t381 + t386 * t655;
t55 = t102 * t386 - t122 * t381;
t629 = qJD(6) * t55 + t381 * t655 - t654 * t386;
t626 = t375 * Ifges(5,5);
t625 = t375 * Ifges(5,6);
t624 = t389 * Ifges(4,2);
t623 = -mrSges(6,1) - t592;
t322 = t552 * t382;
t373 = t387 * pkin(11);
t517 = t363 * t387;
t323 = t373 + t517;
t229 = t322 * t381 + t323 * t386;
t621 = -qJD(6) * t229 - t381 * t653 + t386 * t652;
t347 = t391 * t382;
t349 = pkin(10) * t387 + t373;
t273 = t347 * t386 - t349 * t381;
t620 = qJD(6) * t273 + t381 * t650 + t386 * t651;
t275 = t347 * t381 + t349 * t386;
t619 = -qJD(6) * t275 - t381 * t651 + t386 * t650;
t618 = -qJD(5) * t147 + t657;
t228 = t322 * t386 - t323 * t381;
t616 = qJD(6) * t228 + t381 * t652 + t386 * t653;
t172 = mrSges(5,1) * t374 - mrSges(5,3) * t187;
t57 = -mrSges(6,1) * t114 + mrSges(6,2) * t113;
t615 = t57 - t172;
t614 = t138 * t424;
t413 = t381 * t382 - t386 * t387;
t217 = t413 * t327;
t297 = pkin(5) * t522;
t465 = pkin(5) * t477;
t607 = t465 - t297 - t640;
t605 = t167 * t389 - t168 * t384;
t604 = m(6) + m(7) + m(5);
t603 = qJD(5) + qJD(6);
t531 = t381 * t68;
t21 = t386 * t58 - t531;
t7 = pkin(5) * t186 - pkin(11) * t113 + t15;
t8 = pkin(11) * t114 + t14;
t3 = qJD(6) * t21 + t381 * t7 + t386 * t8;
t4 = -qJD(6) * t22 - t381 * t8 + t386 * t7;
t599 = t4 * mrSges(7,1) - t3 * mrSges(7,2);
t598 = t15 * mrSges(6,1) - t14 * mrSges(6,2);
t464 = m(4) * pkin(8) + mrSges(4,3);
t597 = mrSges(3,2) - t464 + t639;
t427 = -t389 * mrSges(4,1) + t384 * mrSges(4,2);
t403 = m(4) * pkin(2) - t427;
t596 = mrSges(3,1) + t403 - t656;
t533 = t270 * Ifges(6,4);
t131 = t269 * Ifges(6,2) + t309 * Ifges(6,6) + t533;
t326 = t381 * t387 + t382 * t386;
t193 = t326 * t318;
t194 = t413 * t318;
t532 = t319 * Ifges(5,4);
t202 = t318 * Ifges(5,2) + t532 + t625;
t306 = Ifges(5,4) * t318;
t203 = t319 * Ifges(5,1) + t306 + t626;
t50 = t135 * t388 - t383 * t140 - t253 * t479 - t264 * t478;
t47 = -pkin(4) * t374 - t50;
t25 = -pkin(5) * t114 + t47;
t255 = t603 * t413;
t256 = t603 * t326;
t417 = Ifges(6,5) * t387 - Ifges(6,6) * t382;
t541 = Ifges(6,4) * t387;
t419 = -Ifges(6,2) * t382 + t541;
t542 = Ifges(6,4) * t382;
t421 = Ifges(6,1) * t387 - t542;
t43 = t113 * Ifges(6,4) + t114 * Ifges(6,2) + t186 * Ifges(6,6);
t443 = -t477 / 0.2e1;
t265 = Ifges(6,4) * t269;
t132 = t270 * Ifges(6,1) + t309 * Ifges(6,5) + t265;
t506 = t387 * t132;
t447 = t506 / 0.2e1;
t509 = t382 * t131;
t537 = t144 * mrSges(5,3);
t564 = t319 / 0.2e1;
t566 = t318 / 0.2e1;
t568 = -t309 / 0.2e1;
t569 = t301 / 0.2e1;
t570 = -t301 / 0.2e1;
t572 = -t270 / 0.2e1;
t573 = -t269 / 0.2e1;
t577 = t162 / 0.2e1;
t578 = -t162 / 0.2e1;
t579 = t432 / 0.2e1;
t580 = -t432 / 0.2e1;
t157 = Ifges(7,4) * t432;
t75 = t162 * Ifges(7,1) + t301 * Ifges(7,5) + t157;
t583 = t75 / 0.2e1;
t584 = -t75 / 0.2e1;
t540 = Ifges(7,4) * t162;
t74 = Ifges(7,2) * t432 + Ifges(7,6) * t301 + t540;
t585 = t74 / 0.2e1;
t586 = -t74 / 0.2e1;
t587 = Ifges(6,1) * t582 + Ifges(6,4) * t581 + Ifges(6,5) * t575;
t595 = -t22 * (-mrSges(7,2) * t319 - mrSges(7,3) * t193) - t99 * (mrSges(7,1) * t193 - mrSges(7,2) * t194) + (-Ifges(7,1) * t194 - Ifges(7,4) * t193 + Ifges(7,5) * t319) * t578 + (-Ifges(7,4) * t194 - Ifges(7,2) * t193 + Ifges(7,6) * t319) * t580 + (-Ifges(7,5) * t194 - Ifges(7,6) * t193 + Ifges(7,3) * t319) * t570 - t307 * (mrSges(5,1) * t319 + mrSges(5,2) * t318) + (t521 * t82 + t522 * t83 + t637) * mrSges(6,3) + (t269 * t419 + t270 * t421 + t309 * t417) * qJD(5) / 0.2e1 - (-Ifges(5,2) * t319 + t203 + t306 + t506) * t318 / 0.2e1 - (Ifges(5,1) * t318 - t532 + t641) * t319 / 0.2e1 + (t614 + t447) * qJD(5) + (mrSges(7,1) * t25 - mrSges(7,3) * t3 - 0.2e1 * t591) * t413 + (mrSges(7,2) * t25 + 0.2e1 * t590) * t326 + (t21 * t255 - t326 * t4) * mrSges(7,3) + t319 * t642 - t21 * (mrSges(7,1) * t319 + mrSges(7,3) * t194) + t131 * t443 + (Ifges(6,5) * t382 + Ifges(6,6) * t387) * t575 + (Ifges(6,2) * t387 + t542) * t581 + (Ifges(6,1) * t382 + t541) * t582 - t194 * t584 - t193 * t586 + t382 * t587 + t202 * t564 + t509 * t566 + (Ifges(6,3) * t319 + t318 * t417) * t568 + (Ifges(6,5) * t319 + t318 * t421) * t572 + (Ifges(6,6) * t319 + t318 * t419) * t573 + t318 * t537 - (Ifges(7,4) * t577 + Ifges(7,2) * t579 + Ifges(7,6) * t569 + t585 + t633) * t256 + t387 * t43 / 0.2e1 - t375 * (Ifges(5,5) * t318 - Ifges(5,6) * t319) / 0.2e1 + Ifges(5,3) * t374 + Ifges(5,5) * t187 + Ifges(5,6) * t188 - t49 * mrSges(5,2) + t50 * mrSges(5,1) - t319 * t643 - (Ifges(7,1) * t577 + Ifges(7,4) * t579 + Ifges(7,5) * t569 + t583 + t644) * t255 - t318 * t614 + t47 * t425;
t393 = qJD(2) ^ 2;
t571 = t270 / 0.2e1;
t562 = mrSges(6,3) * t82;
t560 = pkin(3) * t388;
t559 = pkin(5) * t270;
t557 = g(3) * t379;
t554 = t83 * mrSges(6,3);
t546 = mrSges(6,3) * t269;
t545 = mrSges(6,3) * t270;
t544 = Ifges(4,4) * t384;
t543 = Ifges(4,4) * t389;
t539 = pkin(5) * qJD(6);
t536 = t145 * mrSges(5,3);
t76 = mrSges(6,1) * t186 - mrSges(6,3) * t113;
t530 = t382 * t76;
t511 = t379 * t390;
t310 = -t380 * t436 + t515;
t505 = (-t246 * t369 + t310 * t371) * mrSges(7,1) + (-t246 * t371 - t310 * t369) * mrSges(7,2);
t312 = t380 * t514 + t437;
t504 = (-t248 * t369 + t312 * t371) * mrSges(7,1) + (-t248 * t371 - t312 * t369) * mrSges(7,2);
t502 = t245 * t364 - t246 * t391;
t501 = t247 * t364 - t248 * t391;
t498 = (-t296 * t369 - t371 * t511) * mrSges(7,1) + (-t296 * t371 + t369 * t511) * mrSges(7,2);
t493 = t295 * t364 - t296 * t391;
t483 = qJD(2) * t385;
t480 = qJD(3) * t389;
t470 = Ifges(7,5) * t38 + Ifges(7,6) * t39 + Ifges(7,3) * t182;
t462 = mrSges(4,3) * t484;
t461 = mrSges(4,3) * t482;
t459 = t382 * t511;
t457 = t387 * t511;
t456 = Ifges(6,5) * t113 + Ifges(6,6) * t114 + Ifges(6,3) * t186;
t450 = t379 * t483;
t449 = t390 * t485;
t442 = t245 * pkin(4) + t246 * pkin(10);
t441 = t247 * pkin(4) + pkin(10) * t248;
t440 = t295 * pkin(4) + pkin(10) * t296;
t435 = t475 / 0.2e1;
t429 = t638 * pkin(3);
t420 = t544 + t624;
t418 = Ifges(4,5) * t389 - Ifges(4,6) * t384;
t315 = t380 * t384 + t385 * t512;
t210 = t314 * t383 + t315 * t388;
t179 = -t210 * t382 - t457;
t410 = -t210 * t387 + t459;
t91 = t179 * t386 + t381 * t410;
t92 = t179 * t381 - t386 * t410;
t414 = t388 * t314 - t315 * t383;
t412 = t314 * pkin(3);
t411 = t470 + t599;
t407 = t327 * t477 - t525;
t339 = -qJD(2) * pkin(2) - t452;
t406 = t339 * (mrSges(4,1) * t384 + mrSges(4,2) * t389);
t405 = t384 * (Ifges(4,1) * t389 - t544);
t400 = -t311 * t384 - t389 * t438;
t398 = t400 * pkin(3);
t166 = qJD(4) * t276 + t334 * t383 - t388 * t335;
t397 = (-t382 * t83 - t387 * t82) * qJD(5) + t637;
t367 = Ifges(4,4) * t482;
t365 = -pkin(4) - t560;
t344 = -qJD(3) * mrSges(4,2) + t461;
t343 = qJD(3) * mrSges(4,1) - t462;
t342 = -t364 - t560;
t329 = t427 * qJD(2);
t317 = Ifges(4,1) * t484 + Ifges(4,5) * qJD(3) + t367;
t316 = Ifges(4,6) * qJD(3) + qJD(2) * t420;
t305 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t337;
t304 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t336;
t283 = -mrSges(5,2) * t375 + mrSges(5,3) * t318;
t281 = -t338 * t384 + t359;
t266 = -mrSges(4,1) * t336 + mrSges(4,2) * t337;
t261 = qJD(3) * t314 + t389 * t449;
t260 = -qJD(3) * t315 - t384 * t449;
t239 = -mrSges(5,1) * t318 + mrSges(5,2) * t319;
t216 = t326 * t327;
t195 = pkin(5) * t520 - t608;
t190 = mrSges(6,1) * t309 - t545;
t189 = -mrSges(6,2) * t309 + t546;
t173 = -mrSges(5,2) * t374 + mrSges(5,3) * t188;
t125 = mrSges(7,1) * t301 - mrSges(7,3) * t162;
t124 = -mrSges(7,2) * t301 + mrSges(7,3) * t432;
t115 = t145 + t297;
t98 = -mrSges(5,1) * t188 + mrSges(5,2) * t187;
t97 = qJD(4) * t210 - t388 * t260 + t261 * t383;
t96 = qJD(4) * t414 + t260 * t383 + t261 * t388;
t93 = pkin(5) * t408 + t166;
t77 = -mrSges(6,2) * t186 + mrSges(6,3) * t114;
t70 = t217 * t603 - t326 * t257;
t69 = -t256 * t327 - t257 * t413;
t62 = qJD(5) * t410 - t382 * t96 + t387 * t450;
t61 = qJD(5) * t179 + t382 * t450 + t387 * t96;
t27 = -mrSges(7,2) * t182 + mrSges(7,3) * t39;
t26 = mrSges(7,1) * t182 - mrSges(7,3) * t38;
t24 = t386 * t67 - t531;
t23 = -t381 * t67 - t529;
t18 = -qJD(6) * t92 - t381 * t61 + t386 * t62;
t17 = qJD(6) * t91 + t381 * t62 + t386 * t61;
t16 = -mrSges(7,1) * t39 + mrSges(7,2) * t38;
t1 = [m(2) * qJDD(1) + t17 * t124 + t18 * t125 + t210 * t173 + t179 * t76 - t410 * t77 + t61 * t189 + t62 * t190 + t91 * t26 + t260 * t343 + t261 * t344 + t92 * t27 + t96 * t283 + t315 * t304 + t314 * t305 + t649 * t97 - (t16 + t615) * t414 + (-m(2) - m(3) - m(4) - t604) * g(3) + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t393 - t266 - t98) * t390 + (-mrSges(3,1) * t393 - mrSges(3,2) * qJDD(2) + (t239 + t329) * qJD(2)) * t385) * t379 + m(4) * (t167 * t315 + t168 * t314 + t260 * t281 + t261 * t282 + (-t293 * t390 + t339 * t483) * t379) + m(5) * (-t144 * t97 + t145 * t96 + t414 * t50 + t210 * t49 + (-t243 * t390 + t307 * t483) * t379) + m(3) * (qJDD(1) * t380 ^ 2 + (t302 * t390 + t303 * t385) * t379) + m(6) * (t138 * t97 - t14 * t410 + t15 * t179 - t414 * t47 + t61 * t83 + t62 * t82) + m(7) * (t17 * t22 + t18 * t21 - t25 * t414 + t3 * t92 + t4 * t91 + t97 * t99); t543 * t631 + t420 * t632 + (-t390 * t557 + t302 + t351) * mrSges(3,1) + (t385 * t557 - t303 + t352) * mrSges(3,2) + (t470 + t456) * t325 / 0.2e1 + (-pkin(2) * t293 - (t339 * t385 + (-t281 * t384 + t282 * t389) * t390) * t487) * m(4) + (Ifges(7,5) * t69 + Ifges(7,6) * t70) * t569 + t309 * (-Ifges(6,5) * t407 - Ifges(6,6) * t408) / 0.2e1 + (Ifges(7,4) * t69 + Ifges(7,2) * t70) * t579 + t269 * (-Ifges(6,4) * t407 - Ifges(6,2) * t408) / 0.2e1 + (t641 / 0.2e1 - t536 + Ifges(7,5) * t577 + Ifges(7,6) * t579 - Ifges(5,4) * t564 - Ifges(5,2) * t566 + Ifges(7,3) * t569 + Ifges(6,5) * t571 - t625 / 0.2e1 + t307 * mrSges(5,1) - t202 / 0.2e1 + t21 * mrSges(7,1) - t22 * mrSges(7,2) + t643 - t642 + t628 / 0.2e1 + t627 / 0.2e1) * t258 + (Ifges(7,1) * t69 + Ifges(7,4) * t70) * t577 + (-Ifges(6,1) * t407 - Ifges(6,4) * t408) * t571 + (-Ifges(7,5) * t217 - Ifges(7,6) * t216) * t576 + (-Ifges(7,4) * t217 - Ifges(7,2) * t216) * t588 + (-Ifges(7,1) * t217 - Ifges(7,4) * t216) * t589 + (-t21 * t69 - t216 * t3 + t217 * t4 + t22 * t70) * mrSges(7,3) + t25 * (mrSges(7,1) * t216 - mrSges(7,2) * t217) + t405 * t435 - t217 * t590 - t216 * t591 + t69 * t583 + t70 * t585 + t519 * t587 + (-t14 * t520 - t15 * t519 + t407 * t82 - t408 * t83) * mrSges(6,3) + (-t144 * t166 + t145 * t610 - t243 * t366 + t276 * t49 + t307 * t606 + t50 * t608) * m(5) + (t138 * t166 + t14 * t147 + t146 * t15 - t47 * t608 + t617 * t83 + t618 * t82) * m(6) - t615 * t608 - t43 * t520 / 0.2e1 + Ifges(3,3) * qJDD(2) + t317 * t480 / 0.2e1 - t316 * t481 / 0.2e1 + qJDD(3) * (Ifges(4,5) * t384 + Ifges(4,6) * t389) + (t406 + t418 * qJD(3) / 0.2e1) * qJD(3) - t329 * t453 - t366 * t98 - t131 * t448 / 0.2e1 + t276 * t173 - pkin(2) * t266 + ((t656 * t390 + (t392 * t604 + t639) * t385) * t379 - t604 * t366 * t511) * g(3) + t195 * t16 + t146 * t76 + t147 * t77 + t99 * (-mrSges(7,1) * t70 + mrSges(7,2) * t69) + t93 * t88 + t55 * t26 + t56 * t27 - (t385 * t464 + t390 * t403) * t557 + (t243 * mrSges(5,2) - t50 * mrSges(5,3) + Ifges(5,1) * t187 + Ifges(5,4) * t188 + Ifges(5,5) * t374 + t132 * t443 + t417 * t575 + t419 * t581 + t421 * t582 + t424 * t47 + (m(5) * t144 - m(6) * t138 - m(7) * t99 - t649) * t452) * t327 + (mrSges(5,1) * t243 - mrSges(5,3) * t49 - Ifges(5,4) * t187 + Ifges(6,5) * t582 + Ifges(7,5) * t589 - Ifges(5,2) * t188 - Ifges(5,6) * t374 + Ifges(6,6) * t581 + Ifges(7,6) * t588 + Ifges(6,3) * t575 + Ifges(7,3) * t576 + t598 + t599) * t325 + t138 * (mrSges(6,1) * t408 - mrSges(6,2) * t407) + (-t604 * (-t312 * t366 - t313 * t392) + t597 * t313 + t596 * t312) * g(1) + (-t604 * (-t310 * t366 - t311 * t392) + t597 * t311 + t596 * t310) * g(2) + (m(4) * ((-t281 * t389 - t282 * t384) * qJD(3) + t605) + t389 * t304 - t384 * t305 - t343 * t480 - t344 * t481) * pkin(8) + (-t281 * t480 - t282 * t481 + t605) * mrSges(4,3) + t606 * t239 - t609 * t166 + t610 * t283 + t617 * t189 + t618 * t190 + (-t537 + t447 + Ifges(5,1) * t564 + Ifges(5,4) * t566 - t509 / 0.2e1 + t626 / 0.2e1 + t307 * mrSges(5,2) + t203 / 0.2e1) * t257 + t293 * t427 + t629 * t124 + t630 * t125 + (t195 * t25 + t21 * t630 + t22 * t629 + t3 * t56 + t4 * t55 + t93 * t99) * m(7) + (Ifges(4,4) * t631 + Ifges(4,2) * t632 - t344 * t452 + t543 * t435) * t389 + (Ifges(4,1) * t337 + Ifges(4,4) * t632 + t343 * t452 - t435 * t624) * t384; -(-Ifges(4,2) * t484 + t317 + t367) * t482 / 0.2e1 + (-m(6) * (t412 + t440) - m(7) * (t412 + t493) - m(5) * t412 - mrSges(4,1) * t314 + mrSges(4,2) * t315 + t634) * g(3) + (-m(6) * (t398 + t442) - m(7) * (t398 + t502) - m(5) * t398 - t400 * mrSges(4,1) - (-t311 * t389 + t384 * t438) * mrSges(4,2) + t636) * g(2) + (-m(6) * (t429 + t441) - m(7) * (t429 + t501) - t638 * mrSges(4,1) - (-t313 * t389 - t384 * t516) * mrSges(4,2) - m(5) * t429 + t635) * g(1) + t640 * t609 + (t365 * t47 + (t138 * t383 + (-t382 * t82 + t387 * t83) * t388) * qJD(4) * pkin(3) + t397 * t363 - t138 * t148 - t82 * t89 - t83 * t90) * m(6) + ((t383 * t49 + t388 * t50 + (-t144 * t383 + t145 * t388) * qJD(4)) * pkin(3) + t144 * t148 - t145 * t149 - t307 * t469) * m(5) + t595 - t363 * t530 + t172 * t560 + t173 * t561 + t319 * t536 + t77 * t517 + (-t149 + t466) * t283 + (-t344 + t461) * t281 + (t343 + t462) * t282 + t316 * t484 / 0.2e1 + Ifges(4,3) * qJDD(3) - t418 * t475 / 0.2e1 - t239 * t469 + t365 * t57 + Ifges(4,6) * t336 + Ifges(4,5) * t337 + t342 * t16 + (-t363 * t477 + t660) * t189 + (-t363 * t476 + t659) * t190 + t228 * t26 + t229 * t27 - t167 * mrSges(4,2) + t168 * mrSges(4,1) - t477 * t554 - t476 * t562 - t393 * t405 / 0.2e1 - qJD(2) * t406 + t607 * t88 + t616 * t124 + t621 * t125 + (t621 * t21 + t616 * t22 + t228 * t4 + t229 * t3 + t25 * t342 + t607 * t99) * m(7); (t609 + t547) * t145 + t595 - m(6) * (t138 * t145 + t82 * t94 + t83 * t95) + (-t501 * g(1) - t502 * g(2) - t493 * g(3) - t25 * t364 + t273 * t4 + t275 * t3 + (-t115 + t465) * t99 + t620 * t22 + t619 * t21) * m(7) + t620 * t124 + (-m(6) * t442 + t636) * g(2) + (-m(6) * t441 + t635) * g(1) + (-m(6) * t440 + t634) * g(3) + t619 * t125 + (t387 * t77 - t530) * pkin(10) + t397 * t593 - t364 * t16 - t144 * t283 + t273 * t26 + t275 * t27 - t95 * t189 - t94 * t190 - t115 * t88 - pkin(4) * t57 + ((-pkin(10) * t190 - t562) * t387 + (pkin(5) * t88 - pkin(10) * t189 - t554) * t382) * qJD(5) - t47 * t594; (-(-t296 * t387 + t459) * mrSges(6,2) - t498 + t623 * (-t296 * t382 - t457)) * g(3) + (-t504 - (-t248 * t387 - t312 * t382) * mrSges(6,2) + t623 * (-t248 * t382 + t312 * t387)) * g(1) + (-t505 - (-t246 * t387 - t310 * t382) * mrSges(6,2) + t623 * (-t246 * t382 + t310 * t387)) * g(2) - (Ifges(7,4) * t578 + Ifges(7,2) * t580 + Ifges(7,6) * t570 + t586 - t633) * t162 + (t21 * mrSges(7,3) + Ifges(7,1) * t578 + Ifges(7,4) * t580 + Ifges(7,5) * t570 + t584 - t644) * t432 + (-t381 * t539 - t23) * t125 + t456 + (t3 * t381 + t386 * t4 + (-t21 * t381 + t22 * t386) * qJD(6)) * t592 + (Ifges(6,5) * t269 - Ifges(6,6) * t270) * t568 + t131 * t571 + (Ifges(6,1) * t269 - t533) * t572 + (t26 * t386 + t27 * t381) * pkin(5) + t598 + t411 - t138 * (mrSges(6,1) * t270 + mrSges(6,2) * t269) + (t545 + t190) * t83 + (t546 - t189) * t82 + (t386 * t539 - t24) * t124 + (-Ifges(6,2) * t270 + t132 + t265) * t573 - t88 * t559 - m(7) * (t21 * t23 + t22 * t24 + t559 * t99); -t99 * (mrSges(7,1) * t162 + mrSges(7,2) * t432) - t21 * t124 + t22 * t125 + (Ifges(7,1) * t432 - t540) * t578 + t74 * t577 + (Ifges(7,5) * t432 - Ifges(7,6) * t162) * t570 - g(1) * t504 - g(2) * t505 - g(3) * t498 + (t162 * t22 + t21 * t432) * mrSges(7,3) + t411 + (-Ifges(7,2) * t162 + t157 + t75) * t580;];
tau  = t1;
