% Calculate vector of inverse dynamics joint torques for
% S6RRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRP3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:39:09
% EndTime: 2019-03-09 16:39:51
% DurationCPUTime: 24.69s
% Computational Cost: add. (16775->798), mult. (36569->1030), div. (0->0), fcn. (26768->14), ass. (0->355)
t606 = Ifges(6,1) + Ifges(7,1);
t604 = Ifges(7,4) + Ifges(6,5);
t634 = -mrSges(6,3) - mrSges(7,2);
t638 = mrSges(5,3) - t634;
t608 = mrSges(6,1) + mrSges(7,1);
t569 = m(7) * pkin(5) + t608;
t390 = sin(qJ(3));
t391 = sin(qJ(2));
t394 = cos(qJ(3));
t395 = cos(qJ(2));
t319 = t390 * t395 + t391 * t394;
t301 = t319 * qJD(1);
t384 = qJD(2) + qJD(3);
t386 = sin(pkin(10));
t387 = cos(pkin(10));
t272 = t301 * t387 + t384 * t386;
t389 = sin(qJ(5));
t393 = cos(qJ(5));
t421 = t301 * t386 - t384 * t387;
t187 = t389 * t272 + t393 * t421;
t465 = qJD(1) * qJD(2);
t327 = qJDD(1) * t395 - t391 * t465;
t328 = qJDD(1) * t391 + t395 * t465;
t318 = t390 * t391 - t394 * t395;
t410 = t318 * qJD(3);
t215 = -qJD(1) * t410 + t327 * t390 + t328 * t394;
t382 = qJDD(2) + qJDD(3);
t192 = -t215 * t386 + t382 * t387;
t193 = t215 * t387 + t382 * t386;
t73 = -qJD(5) * t187 + t389 * t192 + t393 * t193;
t560 = t73 / 0.2e1;
t618 = t393 * t272 - t389 * t421;
t74 = qJD(5) * t618 - t393 * t192 + t389 * t193;
t559 = -t74 / 0.2e1;
t411 = t319 * qJD(3);
t216 = qJD(1) * t411 - t394 * t327 + t328 * t390;
t214 = qJDD(5) + t216;
t547 = t214 / 0.2e1;
t629 = Ifges(6,4) - Ifges(7,5);
t603 = Ifges(7,2) + Ifges(6,3);
t628 = Ifges(6,6) - Ifges(7,6);
t472 = qJD(1) * t395;
t473 = qJD(1) * t391;
t300 = t390 * t473 - t394 * t472;
t294 = qJD(5) + t300;
t541 = t294 / 0.2e1;
t550 = t618 / 0.2e1;
t552 = t187 / 0.2e1;
t553 = -t187 / 0.2e1;
t397 = -pkin(8) - pkin(7);
t345 = t397 * t395;
t322 = qJD(1) * t345;
t302 = t390 * t322;
t344 = t397 * t391;
t321 = qJD(1) * t344;
t309 = qJD(2) * pkin(2) + t321;
t256 = t394 * t309 + t302;
t231 = -t384 * pkin(3) + qJD(4) - t256;
t173 = pkin(4) * t421 + t231;
t381 = t395 * pkin(2);
t372 = t381 + pkin(1);
t342 = t372 * qJD(1);
t227 = pkin(3) * t300 - qJ(4) * t301 - t342;
t480 = t394 * t322;
t257 = t390 * t309 - t480;
t237 = qJ(4) * t384 + t257;
t148 = t387 * t227 - t237 * t386;
t105 = pkin(4) * t300 - pkin(9) * t272 + t148;
t149 = t386 * t227 + t387 * t237;
t115 = -pkin(9) * t421 + t149;
t33 = t105 * t389 + t115 * t393;
t30 = qJ(6) * t294 + t33;
t59 = t187 * pkin(5) - qJ(6) * t618 + t173;
t183 = Ifges(7,5) * t618;
t95 = Ifges(7,6) * t294 + Ifges(7,3) * t187 + t183;
t508 = Ifges(6,4) * t618;
t98 = -Ifges(6,2) * t187 + Ifges(6,6) * t294 + t508;
t632 = -mrSges(6,3) * t33 + mrSges(6,1) * t173 + mrSges(7,1) * t59 - mrSges(7,2) * t30 - t98 / 0.2e1 + t95 / 0.2e1;
t636 = -Ifges(6,2) * t553 + Ifges(7,3) * t552 - t541 * t628 - t550 * t629 + t632;
t503 = qJDD(1) * pkin(1);
t290 = -pkin(2) * t327 - t503;
t106 = pkin(3) * t216 - qJ(4) * t215 - qJD(4) * t301 + t290;
t314 = t328 * pkin(7);
t267 = qJDD(2) * pkin(2) - pkin(8) * t328 - t314;
t313 = t327 * pkin(7);
t274 = pkin(8) * t327 + t313;
t469 = qJD(3) * t394;
t470 = qJD(3) * t390;
t133 = t390 * t267 + t394 * t274 + t309 * t469 + t322 * t470;
t125 = qJ(4) * t382 + qJD(4) * t384 + t133;
t35 = t387 * t106 - t125 * t386;
t24 = pkin(4) * t216 - pkin(9) * t193 + t35;
t36 = t386 * t106 + t387 * t125;
t27 = pkin(9) * t192 + t36;
t467 = qJD(5) * t393;
t468 = qJD(5) * t389;
t5 = t105 * t467 - t115 * t468 + t389 * t24 + t393 * t27;
t1 = qJ(6) * t214 + qJD(6) * t294 + t5;
t134 = t267 * t394 - t390 * t274 - t309 * t470 + t322 * t469;
t126 = -pkin(3) * t382 + qJDD(4) - t134;
t85 = -pkin(4) * t192 + t126;
t10 = pkin(5) * t74 - qJ(6) * t73 - qJD(6) * t618 + t85;
t558 = t74 / 0.2e1;
t635 = mrSges(6,1) * t85 + mrSges(7,1) * t10 - mrSges(7,2) * t1 - Ifges(6,4) * t73 / 0.2e1 - Ifges(6,6) * t214 / 0.2e1 + 0.2e1 * Ifges(7,3) * t558 - mrSges(6,3) * t5 + (-t559 + t558) * Ifges(6,2) + (-t629 + Ifges(7,5)) * t560 + (-t628 + Ifges(7,6)) * t547;
t607 = mrSges(6,2) - mrSges(7,3);
t184 = Ifges(6,4) * t187;
t507 = Ifges(7,5) * t187;
t598 = t604 * t294 + t606 * t618 - t184 + t507;
t385 = qJ(2) + qJ(3);
t378 = sin(t385);
t383 = pkin(10) + qJ(5);
t376 = sin(t383);
t493 = t376 * t378;
t519 = mrSges(5,2) * t386;
t633 = -mrSges(6,2) * t493 - t378 * t519;
t6 = -qJD(5) * t33 + t24 * t393 - t27 * t389;
t3 = -pkin(5) * t214 + qJDD(6) - t6;
t630 = mrSges(6,2) * t85 + mrSges(7,2) * t3 - mrSges(6,3) * t6 - mrSges(7,3) * t10 + Ifges(7,5) * t558 + 0.2e1 * t547 * t604 + 0.2e1 * t560 * t606 + (Ifges(6,4) + t629) * t559;
t625 = t148 * mrSges(5,1);
t624 = t149 * mrSges(5,2);
t250 = pkin(3) * t301 + qJ(4) * t300;
t456 = pkin(2) * t473;
t232 = t250 + t456;
t260 = t321 * t394 + t302;
t163 = t387 * t232 - t260 * t386;
t499 = t300 * t387;
t438 = t301 * pkin(4) + pkin(9) * t499;
t118 = t163 + t438;
t164 = t386 * t232 + t387 * t260;
t500 = t300 * t386;
t463 = pkin(9) * t500;
t135 = t463 + t164;
t316 = t386 * t389 - t393 * t387;
t453 = pkin(2) * t469;
t357 = qJD(4) + t453;
t531 = pkin(2) * t390;
t367 = qJ(4) + t531;
t306 = (-pkin(9) - t367) * t386;
t380 = t387 * pkin(9);
t307 = t367 * t387 + t380;
t420 = t393 * t306 - t307 * t389;
t596 = qJD(5) * t420 - t389 * t118 - t393 * t135 - t316 * t357;
t248 = t306 * t389 + t307 * t393;
t317 = t386 * t393 + t387 * t389;
t595 = -qJD(5) * t248 - t118 * t393 + t135 * t389 - t317 * t357;
t165 = t387 * t250 - t256 * t386;
t122 = t165 + t438;
t166 = t386 * t250 + t387 * t256;
t141 = t463 + t166;
t388 = -pkin(9) - qJ(4);
t337 = t388 * t386;
t338 = qJ(4) * t387 + t380;
t419 = t393 * t337 - t338 * t389;
t591 = -qJD(4) * t316 + qJD(5) * t419 - t389 * t122 - t393 * t141;
t270 = t337 * t389 + t338 * t393;
t590 = -qJD(4) * t317 - qJD(5) * t270 - t122 * t393 + t141 * t389;
t225 = t317 * t300;
t226 = t316 * t300;
t284 = pkin(4) * t500;
t296 = t316 * qJD(5);
t297 = t317 * qJD(5);
t623 = -qJD(6) * t317 + t284 + (t226 + t296) * qJ(6) + (t225 + t297) * pkin(5);
t425 = -t35 * t386 + t36 * t387;
t521 = mrSges(5,1) * t387;
t433 = t519 - t521;
t621 = Ifges(5,5) * t272 - t421 * Ifges(5,6) + Ifges(5,3) * t300 - t187 * t628 + t294 * t603 + t604 * t618;
t620 = t638 * t378;
t368 = pkin(4) * t387 + pkin(3);
t377 = cos(t383);
t619 = (t521 - m(7) * (-qJ(6) * t376 - t368) + t376 * mrSges(7,3) + t569 * t377) * t378;
t613 = -m(6) - m(7);
t549 = t192 / 0.2e1;
t548 = t193 / 0.2e1;
t546 = t216 / 0.2e1;
t610 = t327 / 0.2e1;
t532 = t395 / 0.2e1;
t37 = -mrSges(7,2) * t74 + mrSges(7,3) * t214;
t40 = -mrSges(6,2) * t214 - mrSges(6,3) * t74;
t601 = t37 + t40;
t38 = mrSges(6,1) * t214 - mrSges(6,3) * t73;
t39 = -t214 * mrSges(7,1) + t73 * mrSges(7,2);
t600 = t39 - t38;
t599 = t395 * Ifges(3,2);
t289 = t301 * qJ(6);
t597 = -t289 + t596;
t527 = pkin(5) * t301;
t594 = t527 - t595;
t593 = -t257 + t623;
t592 = -t289 + t591;
t589 = t527 - t590;
t432 = mrSges(5,1) * t386 + mrSges(5,2) * t387;
t588 = t231 * t432;
t255 = pkin(3) * t318 - qJ(4) * t319 - t372;
t276 = t344 * t390 - t345 * t394;
t176 = t387 * t255 - t276 * t386;
t495 = t319 * t387;
t140 = pkin(4) * t318 - pkin(9) * t495 + t176;
t177 = t386 * t255 + t387 * t276;
t496 = t319 * t386;
t155 = -pkin(9) * t496 + t177;
t587 = t389 * t140 + t393 * t155;
t518 = mrSges(7,2) * t187;
t159 = mrSges(7,3) * t294 - t518;
t515 = mrSges(6,3) * t187;
t160 = -mrSges(6,2) * t294 - t515;
t586 = -t159 - t160;
t514 = mrSges(6,3) * t618;
t161 = mrSges(6,1) * t294 - t514;
t517 = mrSges(7,2) * t618;
t162 = -mrSges(7,1) * t294 + t517;
t585 = t161 - t162;
t259 = t321 * t390 - t480;
t454 = pkin(2) * t470;
t584 = t454 - t259 + t623;
t506 = t301 * mrSges(4,3);
t583 = mrSges(4,1) * t384 - mrSges(5,1) * t421 - t272 * mrSges(5,2) - t506;
t582 = t394 * t344 + t345 * t390;
t379 = cos(t385);
t581 = t379 * mrSges(4,1) - mrSges(4,2) * t378;
t580 = t214 * t603 + t604 * t73 - t628 * t74;
t525 = pkin(7) * t395;
t526 = pkin(7) * t391;
t579 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t473) * t525 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t472) * t526;
t219 = -t300 * mrSges(5,2) - mrSges(5,3) * t421;
t220 = mrSges(5,1) * t300 - mrSges(5,3) * t272;
t578 = t387 * t219 - t386 * t220;
t123 = -mrSges(5,2) * t216 + mrSges(5,3) * t192;
t124 = mrSges(5,1) * t216 - mrSges(5,3) * t193;
t577 = t387 * t123 - t386 * t124;
t576 = t313 * t395 + t314 * t391;
t392 = sin(qJ(1));
t396 = cos(qJ(1));
t575 = g(1) * t396 + g(2) * t392;
t487 = t379 * t392;
t574 = t392 * t633 - t487 * t638;
t486 = t379 * t396;
t573 = t396 * t633 - t486 * t638;
t488 = t379 * t388;
t418 = -t368 * t378 - t488;
t530 = pkin(2) * t391;
t572 = -m(7) * (-t488 - t530) - m(6) * (t418 - t530) - m(5) * (-pkin(3) * t378 - t530) + t619;
t571 = -m(6) * t418 + m(7) * t488 + t619;
t570 = -m(3) * pkin(7) + m(5) * t397 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - t432;
t264 = qJD(2) * t319 + t411;
t263 = -qJD(2) * t318 - t410;
t501 = t263 * t387;
t471 = qJD(2) * t391;
t455 = pkin(2) * t471;
t156 = pkin(3) * t264 - qJ(4) * t263 - qJD(4) * t319 + t455;
t448 = qJD(2) * t397;
t323 = t391 * t448;
t324 = t395 * t448;
t194 = qJD(3) * t582 + t323 * t394 + t324 * t390;
t92 = t387 * t156 - t194 * t386;
t55 = pkin(4) * t264 - pkin(9) * t501 + t92;
t502 = t263 * t386;
t93 = t386 * t156 + t387 * t194;
t76 = -pkin(9) * t502 + t93;
t12 = -qJD(5) * t587 - t389 * t76 + t393 * t55;
t491 = t377 * t379;
t492 = t376 * t379;
t568 = t379 * t433 - t491 * t608 + t607 * t492 - t581 - t620;
t567 = m(7) * qJ(6) - t607;
t566 = m(5) * t231 - t583;
t341 = -mrSges(3,1) * t395 + mrSges(3,2) * t391;
t565 = -(m(5) * pkin(3) - t433) * t379 - mrSges(2,1) - m(3) * pkin(1) + t341 - t581;
t564 = t6 * mrSges(6,1) - t3 * mrSges(7,1) - t5 * mrSges(6,2) + t1 * mrSges(7,3);
t557 = Ifges(5,1) * t548 + Ifges(5,4) * t549 + Ifges(5,5) * t546;
t551 = -t618 / 0.2e1;
t542 = -t294 / 0.2e1;
t539 = -t300 / 0.2e1;
t538 = t300 / 0.2e1;
t536 = t301 / 0.2e1;
t533 = t387 / 0.2e1;
t529 = pkin(2) * t394;
t528 = pkin(4) * t386;
t522 = g(3) * t378;
t516 = mrSges(4,3) * t300;
t513 = Ifges(3,4) * t391;
t512 = Ifges(3,4) * t395;
t511 = Ifges(4,4) * t301;
t510 = Ifges(5,4) * t386;
t509 = Ifges(5,4) * t387;
t361 = t378 * qJ(4);
t490 = t378 * t392;
t489 = t378 * t396;
t330 = t379 * t368;
t481 = t392 * t377;
t479 = t396 * t376;
t474 = t379 * pkin(3) + t361;
t447 = -(Ifges(5,4) * t272 - Ifges(5,2) * t421 + Ifges(5,6) * t300) * t386 / 0.2e1;
t446 = (Ifges(5,1) * t272 - Ifges(5,4) * t421 + Ifges(5,5) * t300) * t533;
t22 = t74 * mrSges(6,1) + t73 * mrSges(6,2);
t21 = t74 * mrSges(7,1) - t73 * mrSges(7,3);
t443 = t465 / 0.2e1;
t112 = -t192 * mrSges(5,1) + t193 * mrSges(5,2);
t441 = -t378 * t388 + t330;
t348 = t396 * t372;
t440 = -t392 * t397 + t348;
t229 = pkin(4) * t496 - t582;
t436 = mrSges(3,1) * t391 + mrSges(3,2) * t395;
t434 = mrSges(4,1) * t378 + mrSges(4,2) * t379;
t431 = Ifges(5,1) * t387 - t510;
t430 = t513 + t599;
t429 = -Ifges(5,2) * t386 + t509;
t428 = Ifges(3,5) * t395 - Ifges(3,6) * t391;
t427 = Ifges(5,5) * t387 - Ifges(5,6) * t386;
t32 = t105 * t393 - t115 * t389;
t67 = t140 * t393 - t155 * t389;
t422 = -t148 * t386 + t149 * t387;
t415 = pkin(5) * t491 + qJ(6) * t492 + t441;
t413 = pkin(1) * t436;
t412 = t391 * (Ifges(3,1) * t395 - t513);
t11 = t140 * t467 - t155 * t468 + t389 * t55 + t393 * t76;
t249 = pkin(5) * t316 - qJ(6) * t317 - t368;
t195 = qJD(3) * t276 + t323 * t390 - t394 * t324;
t144 = pkin(4) * t502 + t195;
t235 = -Ifges(4,2) * t300 + Ifges(4,6) * t384 + t511;
t292 = Ifges(4,4) * t300;
t236 = Ifges(4,1) * t301 + Ifges(4,5) * t384 - t292;
t29 = -pkin(5) * t294 + qJD(6) - t32;
t79 = Ifges(5,4) * t193 + Ifges(5,2) * t192 + Ifges(5,6) * t216;
t398 = (Ifges(6,4) * t226 + Ifges(6,6) * t301) * t552 + (Ifges(6,2) * t552 - Ifges(7,3) * t553 + t542 * t628 + t551 * t629 + t632) * t225 + (Ifges(7,5) * t226 + Ifges(7,6) * t301) * t553 + t300 * t588 + t635 * t316 + t636 * t297 + t598 * (-t296 / 0.2e1 - t226 / 0.2e1) + (t226 * t59 - t30 * t301) * mrSges(7,3) + t630 * t317 + (-Ifges(4,2) * t301 + t236 - t292) * t538 + t300 * t446 + t300 * t447 - t32 * (mrSges(6,1) * t301 - mrSges(6,3) * t226) - t29 * (-mrSges(7,1) * t301 + mrSges(7,2) * t226) + t342 * (mrSges(4,1) * t301 - mrSges(4,2) * t300) + t301 * t624 - Ifges(4,6) * t216 + Ifges(4,5) * t215 - t133 * mrSges(4,2) + t134 * mrSges(4,1) + (t226 * t604 + t301 * t603) * t542 + (t226 * t606 + t301 * t604) * t551 + (-t173 * mrSges(6,2) - t29 * mrSges(7,2) + t32 * mrSges(6,3) + t59 * mrSges(7,3) - Ifges(6,4) * t553 - Ifges(7,5) * t552 - t604 * t541 - t606 * t550) * t296 - t301 * t625 + (-t148 * t499 - t149 * t500 + t425) * mrSges(5,3) - (-Ifges(4,1) * t300 - t511 + t621) * t301 / 0.2e1 + Ifges(4,3) * t382 - t384 * (-Ifges(4,5) * t300 - Ifges(4,6) * t301) / 0.2e1 + (-t173 * t226 + t301 * t33) * mrSges(6,2) + t421 * (Ifges(5,6) * t301 - t300 * t429) / 0.2e1 - t272 * (Ifges(5,5) * t301 - t300 * t431) / 0.2e1 + t126 * t433 + t257 * t506 - t256 * t516 + t79 * t533 + t235 * t536 + (Ifges(5,3) * t301 - t300 * t427) * t539 + (Ifges(5,5) * t386 + Ifges(5,6) * t387) * t546 + (Ifges(5,1) * t386 + t509) * t548 + (Ifges(5,2) * t387 + t510) * t549 + t386 * t557;
t374 = Ifges(3,4) * t472;
t371 = -pkin(3) - t529;
t347 = qJ(4) * t486;
t346 = qJ(4) * t487;
t336 = -t368 - t529;
t299 = Ifges(3,1) * t473 + Ifges(3,5) * qJD(2) + t374;
t298 = Ifges(3,6) * qJD(2) + qJD(1) * t430;
t282 = t376 * t392 + t377 * t486;
t281 = t379 * t479 - t481;
t280 = t379 * t481 - t479;
t279 = t376 * t487 + t377 * t396;
t277 = -mrSges(4,2) * t384 - t516;
t252 = mrSges(4,1) * t300 + mrSges(4,2) * t301;
t241 = t316 * t319;
t240 = t317 * t319;
t234 = t249 - t529;
t211 = t259 - t284;
t204 = -t284 + t257;
t201 = -mrSges(4,2) * t382 - mrSges(4,3) * t216;
t200 = mrSges(4,1) * t382 - mrSges(4,3) * t215;
t121 = t263 * t317 + t467 * t495 - t468 * t496;
t120 = -t263 * t316 - t297 * t319;
t113 = pkin(5) * t240 + qJ(6) * t241 + t229;
t111 = mrSges(6,1) * t187 + mrSges(6,2) * t618;
t110 = mrSges(7,1) * t187 - mrSges(7,3) * t618;
t109 = pkin(5) * t618 + qJ(6) * t187;
t57 = -pkin(5) * t318 - t67;
t56 = qJ(6) * t318 + t587;
t25 = pkin(5) * t121 - qJ(6) * t120 + qJD(6) * t241 + t144;
t9 = -pkin(5) * t264 - t12;
t7 = qJ(6) * t264 + qJD(6) * t318 + t11;
t2 = [(Ifges(7,5) * t120 + Ifges(7,6) * t264) * t552 + (t290 * mrSges(4,2) - t134 * mrSges(4,3) + Ifges(4,1) * t215 - Ifges(4,4) * t216 + Ifges(4,5) * t382 + t126 * t432 + t427 * t546 + t429 * t549 + t431 * t548) * t319 + (-t148 * t501 - t149 * t502 - t35 * t495 - t36 * t496) * mrSges(5,3) + (-t256 * t263 - t257 * t264) * mrSges(4,3) + (t395 * t512 + t412) * t443 - (-m(4) * t134 + m(5) * t126 + t112 - t200) * t582 + (-m(4) * t440 - m(5) * t348 + t634 * t489 + t613 * (t368 * t486 - t388 * t489 + t392 * t528 + t440) - t569 * t282 - t567 * t281 + (-(m(5) * qJ(4) + mrSges(5,3)) * t378 + t565) * t396 + t570 * t392) * g(2) + t263 * t588 + t635 * t240 + t636 * t121 + m(7) * (t1 * t56 + t10 * t113 + t25 * t59 + t29 * t9 + t3 * t57 + t30 * t7) - t630 * t241 - t342 * (mrSges(4,1) * t264 + mrSges(4,2) * t263) + m(5) * (t148 * t92 + t149 * t93 + t176 * t35 + t177 * t36) + (Ifges(3,4) * t328 + Ifges(3,2) * t327) * t532 + (t327 * t525 + t328 * t526 + t576) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t576) + (t299 * t532 + t428 * qJD(2) / 0.2e1 - t579) * qJD(2) + (Ifges(5,5) * t193 + Ifges(5,6) * t192 + Ifges(5,3) * t216 + t580) * t318 / 0.2e1 + t430 * t610 + t263 * t446 + t263 * t447 + (t290 * mrSges(4,1) + t35 * mrSges(5,1) - t36 * mrSges(5,2) - t133 * mrSges(4,3) - Ifges(4,4) * t215 + Ifges(5,5) * t548 + Ifges(4,2) * t216 - Ifges(4,6) * t382 + Ifges(5,6) * t549 + Ifges(6,6) * t559 + Ifges(7,6) * t558 + Ifges(5,3) * t546 + t547 * t603 + t560 * t604 + t564) * t318 + (t120 * t604 + t264 * t603) * t541 + (t120 * t606 + t264 * t604) * t550 - t372 * (mrSges(4,1) * t216 + mrSges(4,2) * t215) + (-m(4) * t256 + t566) * t195 + (Ifges(6,4) * t120 + Ifges(6,6) * t264) * t553 + (-mrSges(3,1) * t526 - mrSges(3,2) * t525 + 0.2e1 * Ifges(3,6) * t532) * qJDD(2) + (Ifges(3,1) * t328 + Ifges(3,4) * t610 + Ifges(3,5) * qJDD(2) - t443 * t599) * t391 + t264 * t625 + t252 * t455 + m(6) * (t11 * t33 + t12 * t32 + t144 * t173 + t229 * t85 + t5 * t587 + t6 * t67) + t587 * t40 + m(4) * (t133 * t276 + t194 * t257 - t290 * t372 - t342 * t455) + t328 * t512 / 0.2e1 + t598 * t120 / 0.2e1 + (t120 * t173 - t264 * t33) * mrSges(6,2) - pkin(1) * (-mrSges(3,1) * t327 + mrSges(3,2) * t328) + Ifges(2,3) * qJDD(1) + t194 * t277 + t276 * t201 + t29 * (-mrSges(7,1) * t264 + mrSges(7,2) * t120) - t264 * t235 / 0.2e1 + t32 * (mrSges(6,1) * t264 - mrSges(6,3) * t120) + t263 * t236 / 0.2e1 + t229 * t22 - t341 * t503 + t93 * t219 + t92 * t220 - t79 * t496 / 0.2e1 + t176 * t124 + t177 * t123 + t11 * t160 + t12 * t161 + t9 * t162 + t7 * t159 + t144 * t111 - t298 * t471 / 0.2e1 - t413 * t465 + t25 * t110 + t113 * t21 - t264 * t624 + t67 * t38 + t56 * t37 + t57 * t39 + (t613 * t388 * t490 + t569 * t280 + t567 * t279 + (-m(5) * (-t372 - t361) + m(4) * t372 + t613 * (-t372 - t330) - t565 + t620) * t392 + (t613 * (-t397 + t528) + m(4) * t397 + t570) * t396) * g(1) + t621 * t264 / 0.2e1 + t384 * (Ifges(4,5) * t263 - Ifges(4,6) * t264) / 0.2e1 + (-t120 * t59 + t264 * t30) * mrSges(7,3) - t421 * (Ifges(5,6) * t264 + t263 * t429) / 0.2e1 + t272 * (Ifges(5,5) * t264 + t263 * t431) / 0.2e1 + (Ifges(4,1) * t263 - Ifges(4,4) * t264) * t536 + (Ifges(5,3) * t264 + t263 * t427) * t538 + (Ifges(4,4) * t263 - Ifges(4,2) * t264) * t539 + t495 * t557; (t579 + (t413 - t412 / 0.2e1) * qJD(1)) * qJD(1) + (t126 * t371 - t148 * t163 - t149 * t164 - t231 * t259 + t357 * t422 + t367 * t425) * m(5) + t371 * t112 + t577 * t367 + t578 * t357 + t583 * t259 + t584 * t110 + t398 + t601 * t248 + (t453 - t260) * t277 - (-Ifges(3,2) * t473 + t299 + t374) * t472 / 0.2e1 + (m(6) * t173 + t111 + t566) * t454 - t600 * t420 + (-t173 * t211 + t248 * t5 + t595 * t32 + t596 * t33 + t336 * t85 + t420 * t6) * m(6) + (t1 * t248 + t10 * t234 + t594 * t29 - t3 * t420 + t597 * t30 + t584 * t59) * m(7) + ((t133 * t390 + t134 * t394 + (-t256 * t390 + t257 * t394) * qJD(3)) * pkin(2) + t256 * t259 - t257 * t260 + t342 * t456) * m(4) + t594 * t162 + t595 * t161 + t596 * t160 + t597 * t159 + (t341 - m(4) * t381 - m(5) * (t381 + t474) - m(7) * (t381 + t415) - m(6) * (t381 + t441) + t568) * g(3) + Ifges(3,5) * t328 + t336 * t22 + Ifges(3,6) * t327 - t313 * mrSges(3,2) - t314 * mrSges(3,1) + t234 * t21 - t164 * t219 - t163 * t220 - t211 * t111 + t298 * t473 / 0.2e1 - t428 * t465 / 0.2e1 - t252 * t456 + Ifges(3,3) * qJDD(2) + (m(4) * t530 + t434 + t436) * t575 + (-m(5) * t347 + t396 * t572 + t573) * g(1) + (-m(5) * t346 + t392 * t572 + t574) * g(2) + t200 * t529 + t201 * t531; t575 * t434 - t368 * t22 + t577 * qJ(4) + t578 * qJD(4) + t583 * t257 + t589 * t162 + t590 * t161 + t398 + t601 * t270 + (-pkin(3) * t126 + qJ(4) * t425 + qJD(4) * t422 - t148 * t165 - t149 * t166 - t231 * t257) * m(5) - t600 * t419 + (-t173 * t204 + t270 * t5 + t32 * t590 + t33 * t591 - t368 * t85 + t419 * t6) * m(6) + (t1 * t270 + t10 * t249 + t29 * t589 - t3 * t419 + t30 * t592 + t59 * t593) * m(7) + t591 * t160 + t592 * t159 + t593 * t110 + (-m(5) * t474 - m(6) * t441 - m(7) * t415 + t568) * g(3) - t256 * t277 + t249 * t21 - t166 * t219 - t165 * t220 - t204 * t111 - pkin(3) * t112 + (-m(5) * (-pkin(3) * t489 + t347) + t571 * t396 + t573) * g(1) + (-m(5) * (-pkin(3) * t490 + t346) + t571 * t392 + t574) * g(2); t421 * t219 + t272 * t220 - t586 * t187 + t585 * t618 + t112 + t21 + t22 + (t187 * t30 - t29 * t618 + t10) * m(7) + (t187 * t33 + t32 * t618 + t85) * m(6) + (t148 * t272 + t149 * t421 + t126) * m(5) + (t379 * g(3) - t378 * t575) * (m(5) - t613); t580 + t564 + (-t187 * t604 - t618 * t628) * t542 + (-Ifges(6,2) * t618 - t184 + t598) * t552 + (-m(7) * t30 - t515 + t586) * t32 + (t608 * t281 + t607 * t282) * g(1) + (-m(7) * t29 + t514 + t585) * t33 + (t608 * t279 + t607 * t280) * g(2) + (t376 * t608 + t377 * t607) * t522 + (-t109 * t59 - pkin(5) * t3 + qJ(6) * t1 + qJD(6) * t30 - (-pkin(5) * t376 + qJ(6) * t377) * t522 - g(2) * (-pkin(5) * t279 + qJ(6) * t280) - g(1) * (-pkin(5) * t281 + qJ(6) * t282)) * m(7) + (-t187 * t606 + t183 - t508 + t95) * t551 - t59 * (mrSges(7,1) * t618 + mrSges(7,3) * t187) - t173 * (mrSges(6,1) * t618 - mrSges(6,2) * t187) + qJD(6) * t159 - t109 * t110 + qJ(6) * t37 - pkin(5) * t39 + t30 * t517 + t29 * t518 + t98 * t550 + (Ifges(7,3) * t618 - t507) * t553; t618 * t110 - t294 * t159 + (-g(1) * t281 - g(2) * t279 - g(3) * t493 - t30 * t294 + t59 * t618 + t3) * m(7) + t39;];
tau  = t2;
