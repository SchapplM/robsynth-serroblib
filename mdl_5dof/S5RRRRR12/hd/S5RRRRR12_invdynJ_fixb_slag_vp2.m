% Calculate vector of inverse dynamics joint torques for
% S5RRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR12_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR12_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR12_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR12_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR12_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR12_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR12_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:46:37
% EndTime: 2019-12-31 22:48:07
% DurationCPUTime: 45.61s
% Computational Cost: add. (24714->957), mult. (69486->1346), div. (0->0), fcn. (58170->14), ass. (0->462)
t366 = sin(pkin(5));
t372 = sin(qJ(2));
t376 = cos(qJ(2));
t483 = qJD(1) * qJD(2);
t313 = (qJDD(1) * t376 - t372 * t483) * t366;
t368 = cos(pkin(5));
t482 = qJDD(1) * t368;
t351 = qJDD(2) + t482;
t365 = sin(pkin(6));
t367 = cos(pkin(6));
t241 = -t313 * t365 + t351 * t367 + qJDD(3);
t661 = t241 / 0.2e1;
t692 = 0.2e1 * t661;
t375 = cos(qJ(3));
t507 = t372 * t375;
t371 = sin(qJ(3));
t508 = t371 * t376;
t406 = t367 * t508 + t507;
t352 = qJD(1) * t368 + qJD(2);
t521 = t365 * t371;
t476 = t352 * t521;
t493 = qJD(1) * t366;
t314 = (qJDD(1) * t372 + t376 * t483) * t366;
t411 = t313 * t367 + t351 * t365;
t622 = -t314 * t371 + t375 * t411;
t142 = (-t406 * t493 - t476) * qJD(3) + t622;
t662 = t142 / 0.2e1;
t691 = 0.2e1 * t662;
t505 = t375 * t376;
t509 = t371 * t372;
t408 = t367 * t505 - t509;
t518 = t365 * t375;
t228 = t352 * t518 + t408 * t493;
t141 = qJD(3) * t228 + t314 * t375 + t371 * t411;
t663 = t141 / 0.2e1;
t690 = 0.2e1 * t663;
t572 = pkin(1) * t368;
t361 = t376 * t572;
t349 = qJD(1) * t361;
t570 = pkin(9) * t367;
t469 = pkin(8) + t570;
t442 = t366 * t469;
t410 = t372 * t442;
t258 = -qJD(1) * t410 + t349;
t360 = t372 * t572;
t388 = -t376 * t442 - t360;
t259 = t388 * qJD(1);
t571 = pkin(9) * t365;
t398 = t366 * (pkin(2) * t372 - t376 * t571);
t294 = qJD(1) * t398;
t353 = pkin(9) * t521;
t513 = t367 * t375;
t326 = pkin(2) * t513 - t353;
t514 = t367 * t371;
t651 = t326 * qJD(3) - t258 * t375 - t259 * t514 - t294 * t521;
t186 = -t259 * t365 + t294 * t367;
t407 = t367 * t507 + t508;
t394 = t407 * t366;
t277 = qJD(1) * t394;
t405 = -t367 * t509 + t505;
t291 = t405 * t366;
t278 = qJD(1) * t291;
t492 = qJD(3) * t365;
t689 = -pkin(3) * t277 + pkin(10) * t278 - t186 + (pkin(3) * t371 - pkin(10) * t375) * t492;
t464 = t372 * t493;
t445 = t365 * t464;
t688 = pkin(10) * t445 - t651;
t393 = t406 * t366;
t229 = qJD(1) * t393 + t476;
t138 = qJD(3) * t229 + qJDD(4) - t622;
t135 = Ifges(5,3) * t138;
t463 = t376 * t493;
t280 = t352 * t367 - t365 * t463 + qJD(3);
t370 = sin(qJ(4));
t374 = cos(qJ(4));
t176 = t229 * t374 + t280 * t370;
t81 = -qJD(4) * t176 - t141 * t370 + t241 * t374;
t77 = Ifges(5,6) * t81;
t175 = -t229 * t370 + t280 * t374;
t80 = qJD(4) * t175 + t141 * t374 + t241 * t370;
t78 = Ifges(5,5) * t80;
t22 = t78 + t77 + t135;
t516 = t366 * t376;
t222 = t352 * t571 + (t469 * t516 + t360) * qJD(1);
t227 = pkin(2) * t352 + t258;
t520 = t365 * t372;
t273 = (-pkin(2) * t376 - pkin(9) * t520 - pkin(1)) * t493;
t416 = t227 * t367 + t273 * t365;
t117 = t222 * t375 + t371 * t416;
t103 = pkin(10) * t280 + t117;
t357 = pkin(8) * t516;
t480 = qJD(2) * t572;
t450 = qJD(1) * t480;
t517 = t366 * t372;
t462 = qJD(2) * t517;
t451 = pkin(8) * t462;
t477 = pkin(1) * t482;
t225 = -qJD(1) * t451 + qJDD(1) * t357 + t372 * t477 + t376 * t450;
t171 = pkin(9) * t411 + t225;
t226 = -pkin(8) * t314 - t372 * t450 + t376 * t477;
t174 = pkin(2) * t351 - t314 * t570 + t226;
t573 = pkin(1) * t366;
t205 = -pkin(2) * t313 - qJDD(1) * t573 - t314 * t571;
t490 = qJD(3) * t375;
t458 = t367 * t490;
t460 = t365 * t490;
t491 = qJD(3) * t371;
t45 = t171 * t375 + t174 * t514 + t205 * t521 - t222 * t491 + t227 * t458 + t273 * t460;
t40 = pkin(10) * t241 + t45;
t113 = -t174 * t365 + t205 * t367;
t47 = -pkin(3) * t142 - pkin(10) * t141 + t113;
t488 = qJD(4) * t374;
t489 = qJD(4) * t370;
t170 = -t227 * t365 + t273 * t367;
t99 = -pkin(3) * t228 - pkin(10) * t229 + t170;
t10 = -t103 * t489 + t370 * t47 + t374 * t40 + t488 * t99;
t51 = t103 * t374 + t370 * t99;
t11 = -qJD(4) * t51 - t370 * t40 + t374 * t47;
t440 = -mrSges(5,1) * t11 + mrSges(5,2) * t10;
t687 = Ifges(4,4) * t690 + Ifges(4,2) * t691 + Ifges(4,6) * t692 + t440 + t45 * mrSges(4,3) - t22 / 0.2e1;
t216 = t278 * t370 - t374 * t445;
t519 = t365 * t374;
t321 = t367 * t370 + t371 * t519;
t265 = qJD(4) * t321 + t370 * t460;
t686 = t216 - t265;
t217 = t278 * t374 + t370 * t445;
t320 = -t367 * t374 + t370 * t521;
t264 = -qJD(4) * t320 + t374 * t460;
t685 = t217 - t264;
t328 = pkin(2) * t514 + pkin(9) * t518;
t684 = -qJD(3) * t328 + t258 * t371;
t461 = t365 * t491;
t683 = t277 - t461;
t600 = t80 / 0.2e1;
t682 = Ifges(5,4) * t600;
t300 = pkin(10) * t367 + t328;
t301 = (-pkin(3) * t375 - pkin(10) * t371 - pkin(2)) * t365;
t654 = -t300 * t489 + t301 * t488 + t370 * t689 - t374 * t688;
t653 = t259 * t513 - (-pkin(3) * t464 - t294 * t375) * t365 - t684;
t224 = -t228 + qJD(4);
t599 = t81 / 0.2e1;
t589 = t138 / 0.2e1;
t681 = -pkin(11) * t683 + t654;
t680 = -pkin(4) * t686 + pkin(11) * t685 + t653;
t459 = t367 * t491;
t46 = t375 * (t174 * t367 + t205 * t365) - t371 * t171 - t222 * t490 - t227 * t459 - t273 * t461;
t441 = t46 * mrSges(4,1) - t45 * mrSges(4,2);
t136 = Ifges(4,6) * t142;
t137 = Ifges(4,5) * t141;
t235 = Ifges(4,3) * t241;
t68 = t137 + t136 + t235;
t678 = t441 + t68 / 0.2e1;
t50 = -t103 * t370 + t374 * t99;
t677 = t170 * mrSges(4,1) + t50 * mrSges(5,1) - t51 * mrSges(5,2) - t117 * mrSges(4,3);
t645 = Ifges(5,1) * t600 + Ifges(5,5) * t589;
t611 = Ifges(5,4) * t599 + t645;
t676 = -mrSges(5,3) * t11 + 0.2e1 * t611;
t675 = -t46 * mrSges(4,3) + Ifges(4,1) * t690 + Ifges(4,4) * t691 + Ifges(4,5) * t692;
t630 = t300 * t374 + t301 * t370;
t655 = -qJD(4) * t630 + t370 * t688 + t374 * t689;
t577 = -t280 / 0.2e1;
t580 = -t228 / 0.2e1;
t582 = -t224 / 0.2e1;
t584 = -t176 / 0.2e1;
t586 = -t175 / 0.2e1;
t669 = Ifges(5,5) * t584 - Ifges(4,2) * t580 - Ifges(4,6) * t577 + Ifges(5,6) * t586 + Ifges(5,3) * t582 - t677;
t579 = -t229 / 0.2e1;
t116 = -t222 * t371 + t375 * t416;
t623 = t170 * mrSges(4,2) - t116 * mrSges(4,3);
t668 = Ifges(4,1) * t579 + Ifges(4,5) * t577 - t623;
t369 = sin(qJ(5));
t373 = cos(qJ(5));
t123 = -t176 * t369 + t224 * t373;
t124 = t176 * t373 + t224 * t369;
t173 = qJD(5) - t175;
t52 = Ifges(6,5) * t124 + Ifges(6,6) * t123 + Ifges(6,3) * t173;
t555 = Ifges(5,4) * t176;
t94 = Ifges(5,2) * t175 + Ifges(5,6) * t224 + t555;
t648 = -t94 / 0.2e1 + t52 / 0.2e1;
t161 = pkin(3) * t229 - pkin(10) * t228;
t85 = t116 * t374 + t161 * t370;
t667 = pkin(10) * t489 + pkin(11) * t229 + t85;
t666 = -pkin(10) * qJD(5) * t374 - t117 + t224 * (pkin(4) * t370 - pkin(11) * t374);
t79 = qJDD(5) - t81;
t601 = t79 / 0.2e1;
t37 = -qJD(5) * t124 + t138 * t373 - t369 * t80;
t609 = t37 / 0.2e1;
t36 = qJD(5) * t123 + t138 * t369 + t373 * t80;
t610 = t36 / 0.2e1;
t41 = -pkin(3) * t241 - t46;
t13 = -pkin(4) * t81 - pkin(11) * t80 + t41;
t43 = pkin(11) * t224 + t51;
t102 = -pkin(3) * t280 - t116;
t59 = -pkin(4) * t175 - pkin(11) * t176 + t102;
t20 = -t369 * t43 + t373 * t59;
t5 = pkin(11) * t138 + t10;
t1 = qJD(5) * t20 + t13 * t369 + t373 * t5;
t21 = t369 * t59 + t373 * t43;
t2 = -qJD(5) * t21 + t13 * t373 - t369 * t5;
t619 = -t2 * mrSges(6,1) + t1 * mrSges(6,2);
t7 = Ifges(6,5) * t36 + Ifges(6,6) * t37 + Ifges(6,3) * t79;
t665 = t619 - Ifges(6,5) * t610 + 0.2e1 * Ifges(5,2) * t599 + 0.2e1 * Ifges(5,6) * t589 - Ifges(6,6) * t609 - Ifges(6,3) * t601 + t682 - t7 / 0.2e1;
t575 = t376 / 0.2e1;
t299 = t353 + (-pkin(2) * t375 - pkin(3)) * t367;
t194 = pkin(4) * t320 - pkin(11) * t321 + t299;
t197 = -pkin(11) * t518 + t630;
t122 = t194 * t369 + t197 * t373;
t660 = -qJD(5) * t122 - t369 * t681 + t373 * t680;
t121 = t194 * t373 - t197 * t369;
t659 = qJD(5) * t121 + t369 * t680 + t373 * t681;
t658 = mrSges(5,1) * t102;
t657 = t102 * mrSges(5,2);
t656 = pkin(4) * t683 - t655;
t652 = -(t259 * t367 + t294 * t365) * t375 + t684;
t574 = sin(qJ(1));
t465 = t574 * t376;
t377 = cos(qJ(1));
t503 = t377 * t372;
t397 = t368 * t465 + t503;
t390 = t397 * t365;
t467 = t367 * t574;
t382 = t366 * t467 + t390;
t511 = t369 * t374;
t157 = -t228 * t511 + t229 * t373;
t650 = t369 * t488 + t157;
t466 = t574 * t372;
t504 = t376 * t377;
t322 = -t368 * t504 + t466;
t515 = t366 * t377;
t262 = -t322 * t365 + t367 * t515;
t433 = -mrSges(6,1) * t373 + mrSges(6,2) * t369;
t400 = m(6) * pkin(4) - t433;
t435 = -mrSges(5,1) * t374 + mrSges(5,2) * t370;
t478 = m(6) * pkin(11) + mrSges(6,3);
t629 = t370 * t478 + t374 * t400 - t435;
t637 = m(5) + m(6);
t647 = pkin(3) * t637 + mrSges(4,1) + t629;
t618 = mrSges(6,1) * t20 - mrSges(6,2) * t21;
t323 = t368 * t503 + t465;
t473 = t365 * t515;
t201 = t322 * t514 - t323 * t375 + t371 * t473;
t644 = t201 * t374 + t262 * t370;
t643 = t201 * t370 - t262 * t374;
t587 = t173 / 0.2e1;
t592 = t124 / 0.2e1;
t594 = t123 / 0.2e1;
t642 = Ifges(6,5) * t592 + Ifges(6,6) * t594 + Ifges(6,3) * t587;
t641 = mrSges(5,3) * t51 - t618 - t658;
t638 = -mrSges(5,3) * t10 - t665 - t682;
t636 = mrSges(5,1) + t400;
t452 = mrSges(5,2) - t478;
t340 = -pkin(4) * t374 - pkin(11) * t370 - pkin(3);
t487 = qJD(5) * t369;
t635 = -t340 * t487 + t369 * t667 + t373 * t666;
t485 = qJD(5) * t373;
t634 = t340 * t485 + t369 * t666 - t373 * t667;
t257 = pkin(2) * t368 + t361 - t410;
t474 = t365 * t517;
t496 = pkin(2) * t516 + pkin(9) * t474;
t284 = -t496 - t573;
t183 = -t257 * t365 + t284 * t367;
t252 = -t366 * t408 - t368 * t518;
t251 = t252 * pkin(3);
t253 = t368 * t521 + t393;
t453 = pkin(10) * t253 - t251;
t114 = t183 - t453;
t495 = t357 + t360;
t244 = (t365 * t368 + t367 * t516) * pkin(9) + t495;
t140 = t244 * t375 + t257 * t514 + t284 * t521;
t319 = -t365 * t516 + t367 * t368;
t119 = pkin(10) * t319 + t140;
t633 = t114 * t370 + t119 * t374;
t632 = t370 * t485 + t650;
t506 = t373 * t374;
t158 = t228 * t506 + t229 * t369;
t486 = qJD(5) * t370;
t631 = t369 * t486 - t373 * t488 + t158;
t628 = t10 * t374 - t11 * t370;
t627 = t1 * t373 - t2 * t369;
t139 = -t371 * t244 + t375 * (t257 * t367 + t284 * t365);
t444 = t365 * t462;
t350 = t376 * t480;
t260 = -qJD(2) * t410 + t350;
t261 = t388 * qJD(2);
t295 = qJD(2) * t398;
t88 = -t244 * t491 + t257 * t458 + t260 * t375 + t261 * t514 + t284 * t460 + t295 * t521;
t82 = pkin(10) * t444 + t88;
t184 = t368 * t461 + (qJD(2) * t407 + qJD(3) * t406) * t366;
t185 = t368 * t460 + (qJD(2) * t405 + qJD(3) * t408) * t366;
t187 = -t261 * t365 + t295 * t367;
t97 = pkin(3) * t184 - pkin(10) * t185 + t187;
t19 = -qJD(4) * t633 - t370 * t82 + t374 * t97;
t432 = mrSges(6,1) * t369 + mrSges(6,2) * t373;
t399 = m(6) * pkin(10) + t432;
t621 = -m(5) * pkin(10) + mrSges(4,2) - t399;
t620 = t226 * mrSges(3,1) - t225 * mrSges(3,2) + Ifges(3,5) * t314 + Ifges(3,6) * t313 + Ifges(3,3) * t351;
t306 = -pkin(8) * t464 + t349;
t308 = t495 * qJD(1);
t348 = Ifges(3,4) * t463;
t529 = t376 * Ifges(3,2);
t557 = Ifges(3,4) * t372;
t617 = -(t306 * t376 + t308 * t372) * mrSges(3,3) - (-mrSges(4,1) * t116 + mrSges(4,2) * t117 - Ifges(4,5) * t229 - Ifges(4,6) * t228 - Ifges(4,3) * t280) * t520 - (t529 + t557) * t464 / 0.2e1 + (Ifges(3,1) * t464 + t348) * t575 + (0.2e1 * Ifges(3,5) * t575 - Ifges(3,6) * t372) * t352;
t588 = -t173 / 0.2e1;
t593 = -t124 / 0.2e1;
t595 = -t123 / 0.2e1;
t616 = Ifges(6,5) * t593 + Ifges(6,6) * t595 + Ifges(6,3) * t588 - t618;
t8 = Ifges(6,4) * t36 + Ifges(6,2) * t37 + Ifges(6,6) * t79;
t614 = t8 / 0.2e1;
t613 = Ifges(6,1) * t610 + Ifges(6,4) * t609 + Ifges(6,5) * t601;
t552 = Ifges(6,4) * t124;
t53 = Ifges(6,2) * t123 + Ifges(6,6) * t173 + t552;
t607 = -t53 / 0.2e1;
t606 = t53 / 0.2e1;
t120 = Ifges(6,4) * t123;
t54 = Ifges(6,1) * t124 + Ifges(6,5) * t173 + t120;
t605 = -t54 / 0.2e1;
t604 = t54 / 0.2e1;
t539 = t224 * Ifges(5,3);
t540 = t176 * Ifges(5,5);
t541 = t175 * Ifges(5,6);
t93 = t539 + t540 + t541;
t598 = t93 / 0.2e1;
t531 = t280 * Ifges(4,6);
t538 = t228 * Ifges(4,2);
t556 = Ifges(4,4) * t229;
t133 = t531 + t538 + t556;
t591 = -t133 / 0.2e1;
t223 = Ifges(4,4) * t228;
t532 = t280 * Ifges(4,5);
t536 = t229 * Ifges(4,1);
t134 = t223 + t532 + t536;
t590 = t134 / 0.2e1;
t585 = t175 / 0.2e1;
t583 = t176 / 0.2e1;
t581 = t224 / 0.2e1;
t578 = t229 / 0.2e1;
t576 = t372 / 0.2e1;
t6 = -pkin(4) * t138 - t11;
t566 = t370 * t6;
t559 = mrSges(5,3) * t175;
t558 = mrSges(5,3) * t176;
t554 = Ifges(5,4) * t370;
t553 = Ifges(5,4) * t374;
t551 = Ifges(6,4) * t369;
t550 = Ifges(6,4) * t373;
t528 = t175 * t369;
t527 = t175 * t373;
t526 = t228 * t370;
t525 = t228 * t374;
t522 = t365 * t370;
t512 = t369 * t370;
t510 = t370 * t373;
t266 = -t321 * t369 - t373 * t518;
t166 = qJD(5) * t266 + t264 * t373 + t369 * t461;
t169 = t217 * t373 + t277 * t369;
t502 = t166 - t169;
t409 = -t321 * t373 + t369 * t518;
t167 = qJD(5) * t409 - t264 * t369 + t373 * t461;
t168 = -t217 * t369 + t277 * t373;
t501 = t167 - t168;
t468 = t366 * t574;
t494 = pkin(1) * t377 + pkin(8) * t468;
t479 = pkin(10) * t488;
t455 = t488 / 0.2e1;
t454 = -t486 / 0.2e1;
t448 = -t244 * t490 - t257 * t459 - t260 * t371 - t284 * t461;
t447 = t365 * t468;
t443 = -pkin(1) * t574 + pkin(8) * t515;
t437 = mrSges(4,1) * t252 + mrSges(4,2) * t253;
t192 = t253 * t370 - t319 * t374;
t193 = t253 * t374 + t319 * t370;
t436 = mrSges(5,1) * t192 + mrSges(5,2) * t193;
t147 = -t193 * t369 + t252 * t373;
t148 = t193 * t373 + t252 * t369;
t434 = mrSges(6,1) * t147 - mrSges(6,2) * t148;
t431 = Ifges(5,1) * t374 - t554;
t430 = Ifges(6,1) * t373 - t551;
t429 = Ifges(6,1) * t369 + t550;
t428 = -Ifges(5,2) * t370 + t553;
t427 = -Ifges(6,2) * t369 + t550;
t426 = Ifges(6,2) * t373 + t551;
t425 = Ifges(5,5) * t374 - Ifges(5,6) * t370;
t424 = Ifges(6,5) * t373 - Ifges(6,6) * t369;
t423 = Ifges(6,5) * t369 + Ifges(6,6) * t373;
t61 = pkin(11) * t252 + t633;
t118 = -pkin(3) * t319 - t139;
t72 = pkin(4) * t192 - pkin(11) * t193 + t118;
t26 = t369 * t72 + t373 * t61;
t25 = -t369 * t61 + t373 * t72;
t420 = -pkin(10) * t637 + mrSges(4,2) - mrSges(5,3);
t62 = t114 * t374 - t119 * t370;
t84 = -t116 * t370 + t161 * t374;
t208 = -t300 * t370 + t301 * t374;
t42 = -pkin(4) * t224 - t50;
t404 = t42 * t432;
t403 = t102 * (mrSges(5,1) * t370 + mrSges(5,2) * t374);
t18 = t114 * t488 - t119 * t489 + t370 * t97 + t374 * t82;
t391 = mrSges(3,2) + (-mrSges(4,3) + (-m(4) - t637) * pkin(9)) * t365;
t389 = t397 * t375;
t387 = -pkin(2) * t323 + pkin(9) * t262 + t443;
t386 = t420 - t432;
t198 = t322 * t513 + t323 * t371 + t375 * t473;
t324 = -t368 * t466 + t504;
t381 = pkin(2) * t324 + pkin(9) * t382 + t494;
t380 = -pkin(1) * (mrSges(3,1) * t372 + mrSges(3,2) * t376) + (Ifges(3,1) * t376 - t557) * t576;
t203 = t324 * t375 + (-t367 * t397 + t447) * t371;
t379 = pkin(3) * t203 + t381;
t83 = -t261 * t513 + (-pkin(3) * t462 - t295 * t375) * t365 - t448;
t327 = -pkin(8) * t517 + t361;
t325 = (-mrSges(3,1) * t376 + mrSges(3,2) * t372) * t366;
t317 = t397 * pkin(2);
t315 = t322 * pkin(2);
t312 = t495 * qJD(2);
t311 = t350 - t451;
t305 = -mrSges(3,2) * t352 + mrSges(3,3) * t463;
t304 = mrSges(3,1) * t352 - mrSges(3,3) * t464;
t297 = pkin(10) * t506 + t340 * t369;
t296 = -pkin(10) * t511 + t340 * t373;
t221 = -t324 * t514 - t389;
t219 = -t322 * t375 - t323 * t514;
t202 = t324 * t371 + t367 * t389 - t375 * t447;
t196 = pkin(4) * t518 - t208;
t182 = mrSges(4,1) * t280 - mrSges(4,3) * t229;
t181 = -mrSges(4,2) * t280 + mrSges(4,3) * t228;
t172 = Ifges(5,4) * t175;
t160 = -mrSges(4,1) * t228 + mrSges(4,2) * t229;
t156 = t203 * t374 + t370 * t382;
t155 = t203 * t370 - t374 * t382;
t127 = mrSges(5,1) * t224 - t558;
t126 = -mrSges(5,2) * t224 + t559;
t109 = -qJD(4) * t192 + t185 * t374 + t370 * t444;
t108 = qJD(4) * t193 + t185 * t370 - t374 * t444;
t107 = -mrSges(4,2) * t241 + mrSges(4,3) * t142;
t106 = mrSges(4,1) * t241 - mrSges(4,3) * t141;
t105 = pkin(4) * t176 - pkin(11) * t175;
t104 = -mrSges(5,1) * t175 + mrSges(5,2) * t176;
t101 = t156 * t373 + t202 * t369;
t100 = -t156 * t369 + t202 * t373;
t95 = Ifges(5,1) * t176 + Ifges(5,5) * t224 + t172;
t92 = mrSges(6,1) * t173 - mrSges(6,3) * t124;
t91 = -mrSges(6,2) * t173 + mrSges(6,3) * t123;
t89 = (t261 * t367 + t295 * t365) * t375 + t448;
t87 = -mrSges(4,1) * t142 + mrSges(4,2) * t141;
t71 = -mrSges(6,1) * t123 + mrSges(6,2) * t124;
t66 = -pkin(4) * t229 - t84;
t60 = -pkin(4) * t252 - t62;
t58 = qJD(5) * t147 + t109 * t373 + t184 * t369;
t57 = -qJD(5) * t148 - t109 * t369 + t184 * t373;
t49 = -mrSges(5,2) * t138 + mrSges(5,3) * t81;
t48 = mrSges(5,1) * t138 - mrSges(5,3) * t80;
t38 = -mrSges(5,1) * t81 + mrSges(5,2) * t80;
t31 = t105 * t369 + t373 * t50;
t30 = t105 * t373 - t369 * t50;
t27 = pkin(4) * t108 - pkin(11) * t109 + t83;
t17 = -mrSges(6,2) * t79 + mrSges(6,3) * t37;
t16 = mrSges(6,1) * t79 - mrSges(6,3) * t36;
t15 = -pkin(4) * t184 - t19;
t14 = pkin(11) * t184 + t18;
t12 = -mrSges(6,1) * t37 + mrSges(6,2) * t36;
t4 = -qJD(5) * t26 - t14 * t369 + t27 * t373;
t3 = qJD(5) * t25 + t14 * t373 + t27 * t369;
t9 = [(t102 * t109 - t184 * t51) * mrSges(5,2) + (-t116 * t185 - t117 * t184) * mrSges(4,3) + (Ifges(5,5) * t600 + Ifges(5,6) * t599 + Ifges(5,3) * t589 - t687) * t252 + (Ifges(5,1) * t109 + Ifges(5,5) * t184) * t583 + (Ifges(6,5) * t58 + Ifges(6,6) * t57) * t587 + (Ifges(6,5) * t148 + Ifges(6,6) * t147) * t601 + m(3) * (t225 * t495 + t226 * t327 - t306 * t312 + t308 * t311) + (-m(3) * t494 - t324 * mrSges(3,1) + t397 * mrSges(3,2) - m(4) * t381 - t203 * mrSges(4,1) - mrSges(4,3) * t390 - m(6) * (t156 * pkin(4) + t379) - t101 * mrSges(6,1) - t100 * mrSges(6,2) - m(5) * t379 - t156 * mrSges(5,1) - t377 * mrSges(2,1) + t574 * mrSges(2,2) + t420 * t202 + t452 * t155) * g(2) + (-m(3) * t443 + t323 * mrSges(3,1) - t322 * mrSges(3,2) + t574 * mrSges(2,1) + t377 * mrSges(2,2) - m(4) * t387 - t201 * mrSges(4,1) - t262 * mrSges(4,3) - t636 * t644 - t386 * t198 + t452 * t643 + t637 * (-pkin(3) * t201 - t387)) * g(1) + (-Ifges(5,4) * t583 - Ifges(5,2) * t585 - Ifges(5,6) * t581 - t641 + t642 + t648) * t108 + t620 * t368 + (Ifges(6,4) * t58 + Ifges(6,2) * t57) * t594 + (Ifges(6,4) * t148 + Ifges(6,2) * t147) * t609 - t6 * t434 + t495 * (-mrSges(3,2) * t351 + mrSges(3,3) * t313) + t228 * (Ifges(4,4) * t185 - Ifges(4,2) * t184) / 0.2e1 + m(5) * (t10 * t633 + t102 * t83 + t11 * t62 + t118 * t41 + t18 * t51 + t19 * t50) + t633 * t49 + (-g(2) * mrSges(4,3) * t467 + (-g(1) * t377 - g(2) * t574) * mrSges(3,3) + (mrSges(3,1) * t313 - mrSges(3,2) * t314 + (m(3) * t573 - t325) * qJDD(1)) * pkin(1) + (mrSges(3,3) * t225 + Ifges(3,4) * t314 + Ifges(3,2) * t313 + Ifges(3,6) * t351) * t376 + (-mrSges(3,3) * t226 + Ifges(3,1) * t314 + Ifges(3,4) * t313 + Ifges(3,5) * t351) * t372 + (((Ifges(3,4) * t376 - Ifges(3,2) * t372) * t575 + t380) * t493 + t617) * qJD(2)) * t366 + (Ifges(6,1) * t58 + Ifges(6,4) * t57) * t592 + (Ifges(6,1) * t148 + Ifges(6,4) * t147) * t610 + t638 * t192 + (Ifges(5,4) * t109 + Ifges(5,6) * t184) * t585 + m(6) * (t1 * t26 + t15 * t42 + t2 * t25 + t20 * t4 + t21 * t3 + t6 * t60) + m(4) * (t113 * t183 + t116 * t89 + t117 * t88 + t139 * t46 + t140 * t45 + t170 * t187) + Ifges(2,3) * qJDD(1) + t327 * (mrSges(3,1) * t351 - mrSges(3,3) * t314) + t311 * t305 - t312 * t304 + t25 * t16 + t26 * t17 + (t1 * t147 - t148 * t2 - t20 * t58 + t21 * t57) * mrSges(6,3) + t42 * (-mrSges(6,1) * t57 + mrSges(6,2) * t58) + t60 * t12 + t62 * t48 + t15 * t71 + t3 * t91 + t4 * t92 + t83 * t104 + t109 * t95 / 0.2e1 + t118 * t38 + t18 * t126 + t19 * t127 + t139 * t106 + t140 * t107 + (Ifges(5,5) * t109 + Ifges(5,3) * t184) * t581 + t88 * t181 + t89 * t182 + t183 * t87 + t50 * (mrSges(5,1) * t184 - mrSges(5,3) * t109) + t170 * (mrSges(4,1) * t184 + mrSges(4,2) * t185) + t187 * t160 + t280 * (Ifges(4,5) * t185 - Ifges(4,6) * t184) / 0.2e1 + t41 * t436 + t113 * t437 + t675 * t253 + t676 * t193 + (Ifges(4,5) * t663 + Ifges(4,6) * t662 + Ifges(4,3) * t661 + t678) * t319 + (Ifges(4,1) * t185 - Ifges(4,4) * t184) * t578 + t185 * t590 + t184 * t591 + t184 * t598 + t58 * t604 + t57 * t606 + t148 * t613 + t147 * t614; (-t134 / 0.2e1 + Ifges(4,4) * t580 + t668) * t278 + (-t93 / 0.2e1 + t133 / 0.2e1 - Ifges(4,4) * t579 + t669) * t277 + (t50 * mrSges(5,3) - t657) * t685 + t641 * t686 + (-pkin(2) * t87 + (t113 * mrSges(4,2) + t675) * t371 + (-t113 * mrSges(4,1) - t135 / 0.2e1 - t78 / 0.2e1 - t77 / 0.2e1 + t687) * t375) * t365 + (t94 - t52) * (t216 / 0.2e1 - t265 / 0.2e1) + t620 + (Ifges(5,5) * t217 - Ifges(5,6) * t216) * t582 + (t397 * mrSges(3,1) - t221 * mrSges(4,1) + t391 * t324 - t636 * (t221 * t374 + t324 * t522) + t386 * (t324 * t513 - t371 * t397) + t452 * (t221 * t370 - t324 * t519) - t637 * (pkin(3) * t221 - t317)) * g(1) + (-t291 * mrSges(4,1) - mrSges(4,3) * t474 + t325 - t636 * (t291 * t374 + t370 * t474) + t386 * t394 + t452 * (t291 * t370 - t374 * t474) - t637 * (pkin(3) * t291 + t496)) * g(3) + (-t219 * mrSges(4,1) + mrSges(3,1) * t322 + t391 * t323 - t636 * (t219 * t374 + t323 * t522) + t386 * (-t322 * t371 + t323 * t513) + t452 * (t219 * t370 - t323 * t519) - t637 * (pkin(3) * t219 - t315)) * g(2) + t654 * t126 + (t10 * t630 + t102 * t653 + t11 * t208 + t299 * t41 + t50 * t655 + t51 * t654) * m(5) + t655 * t127 + t656 * t71 + t659 * t91 + (t1 * t122 + t121 * t2 + t196 * t6 + t20 * t660 + t21 * t659 + t42 * t656) * m(6) + t660 * t92 + (-t217 / 0.2e1 + t264 / 0.2e1) * t95 + (Ifges(5,4) * t217 - Ifges(5,2) * t216) * t586 + (Ifges(5,1) * t217 - Ifges(5,4) * t216) * t584 + (t167 / 0.2e1 - t168 / 0.2e1) * t53 + t651 * t181 + t652 * t182 + (-pkin(2) * t113 * t365 + g(1) * t317 + g(2) * t315 - g(3) * t496 + t116 * t652 + t117 * t651 - t170 * t186 + t326 * t46 + t328 * t45) * m(4) + t653 * t104 + t630 * t49 + (-Ifges(6,4) * t409 + Ifges(6,2) * t266) * t609 + (t1 * t266 + t2 * t409 - t20 * t502 + t21 * t501) * mrSges(6,3) + (-Ifges(6,1) * t409 + Ifges(6,4) * t266) * t610 + t6 * (-mrSges(6,1) * t266 - mrSges(6,2) * t409) + (-Ifges(6,5) * t409 + Ifges(6,6) * t266) * t601 - t409 * t613 + (-mrSges(6,1) * t501 + mrSges(6,2) * t502) * t42 + (-t376 * t348 / 0.2e1 + (t529 * t576 - t380) * t493 - t617) * t493 + (mrSges(5,1) * t41 + t638) * t320 + (t166 / 0.2e1 - t169 / 0.2e1) * t54 + t326 * t106 + t328 * t107 + t121 * t16 + t122 * t17 - t186 * t160 + (mrSges(5,2) * t41 + t676) * t321 + ((t223 / 0.2e1 + t536 / 0.2e1 + t532 / 0.2e1 + t590 + t623) * t375 + (-t538 / 0.2e1 - t556 / 0.2e1 - t531 / 0.2e1 + t541 / 0.2e1 + t540 / 0.2e1 + t539 / 0.2e1 + t591 + t598 + t677) * t371) * t492 + (t137 / 0.2e1 + t136 / 0.2e1 + t235 / 0.2e1 + t678) * t367 + t196 * t12 + t208 * t48 + t299 * t38 - t306 * t305 + t308 * t304 + (Ifges(5,5) * t264 - Ifges(5,6) * t265) * t581 + (Ifges(5,1) * t264 - Ifges(5,4) * t265) * t583 + (Ifges(5,4) * t264 - Ifges(5,2) * t265) * t585 + (Ifges(6,5) * t166 + Ifges(6,6) * t167 + Ifges(6,3) * t265) * t587 + (Ifges(6,5) * t169 + Ifges(6,6) * t168 + Ifges(6,3) * t216) * t588 + (Ifges(6,1) * t166 + Ifges(6,4) * t167 + Ifges(6,5) * t265) * t592 + (Ifges(6,1) * t169 + Ifges(6,4) * t168 + Ifges(6,5) * t216) * t593 + (Ifges(6,4) * t166 + Ifges(6,2) * t167 + Ifges(6,6) * t265) * t594 + (Ifges(6,4) * t169 + Ifges(6,2) * t168 + Ifges(6,6) * t216) * t595 + t266 * t614; (t425 * t582 + t428 * t586 + t431 * t584 - t403 + t668) * t228 + t669 * t229 + (pkin(10) * t49 + (t424 * t587 + t427 * t594 + t430 * t592) * qJD(4) + t665) * t374 + (-t102 * t117 - t50 * t84 - t51 * t85 - pkin(3) * t41 + ((-t370 * t51 - t374 * t50) * qJD(4) + t628) * pkin(10)) * m(5) + (t616 - t648) * t526 + (-t479 - t84) * t127 + (-t104 + t182) * t117 + t554 * t599 + t553 * t600 + (t223 + t134) * t580 + (-t48 + t12) * pkin(10) * t370 + (-t556 + t93) * t579 + (t202 * t647 + t203 * t621) * g(1) + (t198 * t647 - t201 * t621) * g(2) + (-pkin(10) * t126 + t618 + t648) * t489 + (Ifges(6,5) * t158 + Ifges(6,6) * t157) * t588 + t441 - t8 * t512 / 0.2e1 + (Ifges(6,1) * t158 + Ifges(6,4) * t157) * t593 + t68 + t650 * t607 + (t175 * t428 + t176 * t431 + t224 * t425) * qJD(4) / 0.2e1 + (-t423 * t587 - t426 * t594 - t429 * t592) * t486 + (Ifges(6,4) * t158 + Ifges(6,2) * t157) * t595 + (-g(1) * t203 + g(2) * t201 - g(3) * t253 + (-t489 + t526) * t51 + (-t488 + t525) * t50 + t628) * mrSges(5,3) + (qJD(4) * t642 + t424 * t601 + t427 * t609 + t430 * t610 + t611 + t645) * t370 + (-m(5) * t453 + m(6) * t251 + t252 * t629 - t253 * t399 + t437) * g(3) + (mrSges(6,1) * t632 - mrSges(6,2) * t631) * t42 + (-t1 * t512 - t2 * t510 + t20 * t631 - t21 * t632) * mrSges(6,3) + t369 * t54 * t454 + (t479 - t66) * t71 - pkin(3) * t38 + t634 * t91 + t635 * t92 + (-t42 * t66 + t1 * t297 + t2 * t296 + (t42 * t488 + t566) * pkin(10) + t634 * t21 + t635 * t20) * m(6) + qJD(4) * t403 - t85 * t126 - t116 * t181 + t41 * t435 + (t454 * t53 + t455 * t54) * t373 + (-t525 / 0.2e1 + t455) * t95 + t296 * t16 + t297 * t17 + t432 * t566 + t133 * t578 + t158 * t605 + t510 * t613; (-t452 * t644 - t636 * t643) * g(2) + (-pkin(4) * t6 - t20 * t30 - t21 * t31) * m(6) + ((-t487 + t528) * t21 + (-t485 + t527) * t20 + t627) * mrSges(6,3) + (m(6) * ((-t20 * t373 - t21 * t369) * qJD(5) + t627) - t92 * t485 - t91 * t487 - t369 * t16 + t373 * t17) * pkin(11) + (Ifges(5,1) * t584 + Ifges(5,5) * t582 + t424 * t588 + t427 * t595 + t430 * t593 - t404 - t657) * t175 + (-Ifges(5,2) * t586 - Ifges(5,6) * t582 + t616 - t658) * t176 - t440 + t6 * t433 + (t192 * t400 - t193 * t478 + t436) * g(3) + t22 + (t123 * t427 + t124 * t430 + t173 * t424) * qJD(5) / 0.2e1 + (-t555 + t52) * t584 + (-m(6) * t42 + t127 + t558 - t71) * t51 + (-t126 + t559) * t50 + (t155 * t636 + t156 * t452) * g(1) + (t172 + t95) * t586 - t31 * t91 - t30 * t92 + qJD(5) * t404 - pkin(4) * t12 + t94 * t583 + t423 * t601 + t485 * t604 + t527 * t605 + t528 * t606 + t487 * t607 + t426 * t609 + t429 * t610 + t369 * t613 + t373 * t614; -t42 * (mrSges(6,1) * t124 + mrSges(6,2) * t123) + (Ifges(6,1) * t123 - t552) * t593 + t53 * t592 + (Ifges(6,5) * t123 - Ifges(6,6) * t124) * t588 - t20 * t91 + t21 * t92 - g(1) * (mrSges(6,1) * t100 - mrSges(6,2) * t101) - g(2) * ((t198 * t373 + t369 * t644) * mrSges(6,1) + (-t198 * t369 + t373 * t644) * mrSges(6,2)) - g(3) * t434 + (t123 * t20 + t124 * t21) * mrSges(6,3) + t7 + (-Ifges(6,2) * t124 + t120 + t54) * t595 - t619;];
tau = t9;
