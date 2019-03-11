% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 09:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPPRR8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR8_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR8_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR8_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR8_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:23:02
% EndTime: 2019-03-09 09:23:25
% DurationCPUTime: 13.43s
% Computational Cost: add. (23706->723), mult. (49183->992), div. (0->0), fcn. (52111->8), ass. (0->361)
t420 = sin(pkin(10));
t421 = cos(pkin(10));
t423 = sin(qJ(5));
t598 = cos(qJ(5));
t365 = t598 * t420 - t423 * t421;
t422 = sin(qJ(6));
t425 = cos(qJ(6));
t444 = t423 * t420 + t421 * t598;
t453 = t365 * t425 - t422 * t444;
t480 = -t365 * t422 - t425 * t444;
t507 = Ifges(7,5) * t480 - Ifges(7,6) * t453;
t587 = -pkin(8) + qJ(3);
t382 = t587 * t420;
t384 = t587 * t421;
t294 = t423 * t382 + t384 * t598;
t243 = -pkin(9) * t444 + t294;
t293 = t598 * t382 - t384 * t423;
t452 = -pkin(9) * t365 + t293;
t136 = t243 * t425 + t422 * t452;
t661 = -t243 * t422 + t425 * t452;
t717 = -t136 * mrSges(7,1) - t661 * mrSges(7,2);
t26 = t507 + t717;
t725 = t26 * qJD(6);
t424 = sin(qJ(2));
t338 = t365 * t424;
t339 = t444 * t424;
t237 = -t425 * t338 + t339 * t422;
t426 = cos(qJ(2));
t234 = -t338 * t422 - t339 * t425;
t580 = Ifges(7,4) * t234;
t116 = -Ifges(7,2) * t237 + t426 * Ifges(7,6) - t580;
t192 = -mrSges(7,2) * t426 - mrSges(7,3) * t237;
t194 = mrSges(7,1) * t426 + mrSges(7,3) * t234;
t413 = t424 * qJ(4);
t395 = t421 * t413;
t662 = pkin(7) + (pkin(3) + pkin(4)) * t420;
t295 = -t424 * t662 + t395;
t221 = -pkin(5) * t338 + t295;
t379 = -t421 * pkin(3) - t420 * qJ(4) - pkin(2);
t356 = t421 * pkin(4) - t379;
t287 = pkin(5) * t444 + t356;
t645 = -t136 / 0.2e1;
t675 = t237 * mrSges(7,2);
t696 = t234 * mrSges(7,1);
t702 = -t696 - t675;
t676 = t480 * mrSges(7,2);
t695 = t453 * mrSges(7,1);
t703 = t695 + t676;
t705 = -Ifges(7,1) * t237 + t580;
t722 = t661 / 0.2e1;
t724 = (t705 / 0.4e1 - t116 / 0.4e1) * t453 + t192 * t722 + t194 * t645 + t221 * t703 / 0.2e1 + t287 * t702 / 0.2e1;
t723 = 0.2e1 * mrSges(7,2);
t381 = -pkin(2) * t426 - t424 * qJ(3) - pkin(1);
t517 = t420 * t426;
t398 = pkin(7) * t517;
t417 = t426 * pkin(3);
t255 = pkin(4) * t426 + t398 + t417 + (-pkin(8) * t424 - t381) * t421;
t514 = t421 * t426;
t310 = pkin(7) * t514 + t420 * t381;
t297 = -qJ(4) * t426 + t310;
t518 = t420 * t424;
t281 = pkin(8) * t518 + t297;
t145 = t598 * t255 - t423 * t281;
t113 = -t339 * pkin(9) + t145;
t110 = pkin(5) * t426 + t113;
t146 = t423 * t255 + t281 * t598;
t114 = t338 * pkin(9) + t146;
t539 = t114 * t422;
t56 = t110 * t425 - t539;
t65 = t113 * t425 - t539;
t721 = -t56 + t65;
t720 = t221 * t702;
t719 = t287 * t703;
t446 = t422 * t598 + t425 * t423;
t342 = t446 * t426;
t445 = t422 * t423 - t425 * t598;
t343 = t445 * t426;
t478 = -t426 * t598 / 0.2e1;
t652 = m(7) * pkin(5);
t504 = -t652 / 0.2e1;
t650 = -mrSges(7,2) / 0.2e1;
t712 = mrSges(7,1) / 0.2e1;
t669 = t342 * t712 + t343 * t650;
t718 = (t342 * t425 + t343 * t422) * t504 + mrSges(6,2) * t478 - t669;
t557 = t453 * Ifges(7,4);
t154 = Ifges(7,2) * t480 + t557;
t706 = Ifges(7,1) * t480 - t557;
t716 = -t154 / 0.4e1 + t706 / 0.4e1;
t713 = 0.2e1 * mrSges(7,1);
t640 = t237 / 0.2e1;
t711 = -t480 / 0.2e1;
t710 = -t696 / 0.2e1;
t709 = m(7) * t287;
t538 = t114 * t425;
t57 = t110 * t422 + t538;
t64 = -t113 * t422 - t538;
t708 = t57 + t64;
t448 = -t675 / 0.2e1 + t710;
t511 = t237 * t650 + t710;
t70 = -t448 + t511;
t677 = Ifges(7,5) * t237;
t698 = Ifges(7,6) * t234;
t510 = -t677 + t698;
t228 = Ifges(7,4) * t237;
t118 = -Ifges(7,1) * t234 + t426 * Ifges(7,5) - t228;
t690 = Ifges(7,2) * t234 - t228;
t704 = t690 + t118;
t484 = -mrSges(7,1) * t446 + mrSges(7,2) * t445;
t701 = qJD(6) * t484;
t687 = t453 / 0.2e1;
t77 = -t695 / 0.2e1 + mrSges(7,1) * t687 - t480 * t650 - t676 / 0.2e1;
t486 = -t698 / 0.2e1 + t677 / 0.2e1;
t699 = -t453 / 0.2e1;
t76 = t699 * t713 + t711 * t723;
t641 = t234 / 0.2e1;
t603 = -t424 / 0.2e1;
t683 = -t237 / 0.2e1;
t340 = t426 * t365;
t341 = t444 * t426;
t240 = t340 * t422 + t341 * t425;
t635 = t240 / 0.2e1;
t236 = t340 * t425 - t341 * t422;
t638 = t236 / 0.2e1;
t450 = Ifges(7,5) * t635 + Ifges(7,6) * t638;
t541 = qJ(3) * t426;
t385 = t424 * pkin(2) - t541;
t492 = -pkin(7) * t420 - pkin(3);
t266 = (-pkin(8) * t426 - t385) * t421 + (-pkin(4) + t492) * t424;
t515 = t421 * t424;
t323 = -pkin(7) * t515 + t420 * t385;
t299 = t323 + t413;
t282 = pkin(8) * t517 + t299;
t149 = t598 * t266 - t282 * t423;
t111 = -pkin(5) * t424 - pkin(9) * t341 + t149;
t150 = t423 * t266 + t598 * t282;
t120 = pkin(9) * t340 + t150;
t62 = t111 * t425 - t120 * t422;
t63 = t111 * t422 + t120 * t425;
t476 = Ifges(7,3) * t603 + t62 * t712 + t63 * t650 + t450;
t270 = Ifges(7,4) * t480;
t689 = -Ifges(7,2) * t453 + t270;
t682 = Ifges(5,4) + Ifges(4,5);
t681 = Ifges(4,6) - Ifges(5,6);
t674 = t445 * t683;
t129 = mrSges(7,1) * t237 - mrSges(7,2) * t234;
t244 = -mrSges(6,1) * t338 + mrSges(6,2) * t339;
t673 = t129 + t244;
t152 = -mrSges(7,1) * t480 + mrSges(7,2) * t453;
t284 = mrSges(6,1) * t444 + mrSges(6,2) * t365;
t672 = t152 + t284;
t329 = Ifges(6,4) * t338;
t231 = Ifges(6,1) * t339 + t426 * Ifges(6,5) + t329;
t671 = -Ifges(6,2) * t339 + t231 + t329;
t360 = Ifges(6,4) * t444;
t286 = t365 * Ifges(6,1) - t360;
t670 = -Ifges(6,2) * t365 + t286 - t360;
t668 = -Ifges(6,5) * t444 - Ifges(6,6) * t365 + t507;
t666 = Ifges(6,5) * t338 - Ifges(6,6) * t339 + t510;
t516 = t421 * t385;
t322 = pkin(7) * t518 + t516;
t665 = -t322 * t420 + t323 * t421;
t304 = t424 * t492 - t516;
t664 = t299 * t421 + t304 * t420;
t663 = mrSges(4,2) / 0.2e1 - mrSges(5,3) / 0.2e1;
t659 = m(4) / 0.2e1;
t658 = m(5) / 0.2e1;
t657 = -m(6) / 0.2e1;
t656 = m(6) / 0.2e1;
t655 = -m(7) / 0.2e1;
t654 = m(7) / 0.2e1;
t651 = mrSges(6,1) / 0.2e1;
t646 = -t661 / 0.2e1;
t156 = Ifges(7,1) * t453 + t270;
t643 = t156 / 0.2e1;
t636 = -t234 / 0.2e1;
t632 = t480 / 0.2e1;
t300 = -mrSges(6,2) * t426 + t338 * mrSges(6,3);
t625 = -t300 / 0.2e1;
t624 = t338 / 0.2e1;
t622 = -t339 / 0.2e1;
t621 = t339 / 0.2e1;
t620 = t339 / 0.4e1;
t619 = t340 / 0.2e1;
t618 = t341 / 0.2e1;
t617 = t444 / 0.4e1;
t616 = -t444 / 0.2e1;
t615 = -t365 / 0.2e1;
t614 = t365 / 0.2e1;
t613 = t365 / 0.4e1;
t611 = -t445 / 0.2e1;
t610 = t446 / 0.2e1;
t609 = -t446 / 0.2e1;
t608 = -t420 / 0.2e1;
t607 = t420 / 0.2e1;
t606 = t421 / 0.2e1;
t605 = -t422 / 0.2e1;
t604 = -t423 / 0.2e1;
t601 = t425 / 0.2e1;
t599 = t426 / 0.2e1;
t595 = pkin(7) * t426;
t594 = t339 * pkin(5);
t593 = t56 * mrSges(7,2);
t592 = t57 * mrSges(7,1);
t589 = t64 * mrSges(7,1);
t588 = t65 * mrSges(7,2);
t30 = (-(t699 + t687) * t446 - (t711 + t632) * t445) * mrSges(7,3);
t586 = t30 * qJD(5);
t585 = mrSges(6,3) * t444;
t584 = Ifges(4,4) * t420;
t583 = Ifges(4,4) * t421;
t582 = Ifges(6,4) * t339;
t581 = Ifges(6,4) * t365;
t579 = Ifges(5,5) * t420;
t578 = Ifges(5,5) * t421;
t574 = pkin(5) * qJD(5);
t567 = t236 * mrSges(7,1);
t564 = t240 * mrSges(7,2);
t561 = t480 * mrSges(7,3);
t558 = t453 * mrSges(7,3);
t117 = Ifges(7,4) * t240 + Ifges(7,2) * t236 - Ifges(7,6) * t424;
t119 = Ifges(7,1) * t240 + Ifges(7,4) * t236 - Ifges(7,5) * t424;
t130 = t564 - t567;
t193 = mrSges(7,2) * t424 + mrSges(7,3) * t236;
t195 = -mrSges(7,1) * t424 - mrSges(7,3) * t240;
t396 = qJ(4) * t514;
t296 = t426 * t662 - t396;
t222 = t340 * pkin(5) + t296;
t229 = Ifges(6,2) * t338 + t426 * Ifges(6,6) + t582;
t230 = Ifges(6,4) * t341 + Ifges(6,2) * t340 - Ifges(6,6) * t424;
t232 = Ifges(6,1) * t341 + Ifges(6,4) * t340 - Ifges(6,5) * t424;
t552 = t341 * mrSges(6,2);
t553 = t340 * mrSges(6,1);
t245 = t552 - t553;
t309 = t421 * t381 - t398;
t298 = -t309 + t417;
t546 = t426 * mrSges(6,1);
t554 = t339 * mrSges(6,3);
t301 = t546 - t554;
t302 = mrSges(6,2) * t424 + mrSges(6,3) * t340;
t303 = -mrSges(6,1) * t424 - mrSges(6,3) * t341;
t330 = Ifges(5,6) * t424 + (t420 * Ifges(5,3) + t578) * t426;
t331 = Ifges(4,6) * t424 + (-t420 * Ifges(4,2) + t583) * t426;
t332 = Ifges(5,4) * t424 + (Ifges(5,1) * t421 + t579) * t426;
t333 = Ifges(4,5) * t424 + (Ifges(4,1) * t421 - t584) * t426;
t491 = pkin(3) * t420 + pkin(7);
t334 = t424 * t491 - t395;
t335 = t426 * t491 - t396;
t475 = t420 * mrSges(5,1) - t421 * mrSges(5,3);
t357 = t475 * t424;
t358 = t475 * t426;
t359 = (t420 * mrSges(4,1) + t421 * mrSges(4,2)) * t426;
t370 = mrSges(4,2) * t426 - mrSges(4,3) * t518;
t371 = -t424 * mrSges(4,2) - mrSges(4,3) * t517;
t372 = -t426 * mrSges(4,1) - mrSges(4,3) * t515;
t373 = t426 * mrSges(5,1) + mrSges(5,2) * t515;
t374 = t424 * mrSges(4,1) - mrSges(4,3) * t514;
t397 = mrSges(5,2) * t514;
t547 = t424 * mrSges(5,1);
t375 = t397 - t547;
t376 = -mrSges(5,2) * t517 + t424 * mrSges(5,3);
t377 = -mrSges(5,2) * t518 - mrSges(5,3) * t426;
t451 = Ifges(6,5) * t618 + Ifges(6,6) * t619;
t3 = (-pkin(1) * mrSges(3,2) + (m(4) * pkin(7) ^ 2 + Ifges(3,1) - Ifges(3,2) - Ifges(4,3) - Ifges(5,2) - Ifges(6,3) - Ifges(7,3) + (pkin(7) * mrSges(4,2) + (Ifges(4,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t421) * t421 + (pkin(7) * mrSges(4,1) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t420 + (-Ifges(4,4) + Ifges(5,5)) * t421) * t420) * t424 + t450 + t451 + (t420 * t681 - t421 * t682 + Ifges(3,4)) * t426) * t426 + (Ifges(7,5) * t641 + Ifges(7,6) * t640 + Ifges(6,5) * t622 - Ifges(6,6) * t338 / 0.2e1 + pkin(7) * t359 - pkin(1) * mrSges(3,1) - Ifges(3,4) * t424 + (t333 / 0.2e1 + t332 / 0.2e1 + (Ifges(4,5) / 0.2e1 + Ifges(5,4) / 0.2e1) * t424) * t421 + (t330 / 0.2e1 - t331 / 0.2e1 + (-Ifges(4,6) / 0.2e1 + Ifges(5,6) / 0.2e1) * t424) * t420) * t424 + m(5) * (t297 * t299 + t298 * t304 + t334 * t335) + m(6) * (t145 * t149 + t146 * t150 - t295 * t296) + m(7) * (-t221 * t222 + t56 * t62 + t57 * t63) + t117 * t683 + m(4) * (t309 * t322 + t310 * t323) + t323 * t370 + t310 * t371 + t322 * t372 + t304 * t373 + t309 * t374 + t298 * t375 + t297 * t376 + t299 * t377 + t335 * t357 + t334 * t358 + t145 * t303 + t295 * t245 - t296 * t244 + t150 * t300 + t149 * t301 + t146 * t302 - t222 * t129 + t221 * t130 + t63 * t192 + t57 * t193 + t62 * t194 + t56 * t195 + t231 * t618 + t229 * t619 + t232 * t621 + t230 * t624 + t118 * t635 + t119 * t636 + t116 * t638;
t556 = t3 * qJD(1);
t555 = t339 * mrSges(6,1);
t549 = t365 * mrSges(6,1);
t548 = t365 * mrSges(6,3);
t545 = t56 * t480;
t544 = t57 * t453;
t470 = Ifges(6,1) * t338 - t582;
t474 = t338 * mrSges(6,2) + t555;
t6 = m(7) * (t221 * t594 + t56 * t64 + t57 * t65) + t65 * t192 + t64 * t194 + t129 * t594 + t720 + t705 * t636 + t116 * t641 + t470 * t621 + t229 * t622 + t295 * t474 - t146 * t301 + t145 * t300 + (t234 * t57 + t237 * t56) * mrSges(7,3) + (-t145 * t338 - t146 * t339) * mrSges(6,3) + t671 * t624 + t666 * t599 + t704 * t683;
t543 = t6 * qJD(1);
t7 = t56 * t192 - t57 * t194 + t720 + t510 * t599 - (-t57 * mrSges(7,3) + t705 / 0.2e1 - t116 / 0.2e1) * t234 - (-t56 * mrSges(7,3) + t118 / 0.2e1 + t690 / 0.2e1) * t237;
t542 = t7 * qJD(1);
t19 = t237 * t192 - t234 * t194 - t338 * t300 + t339 * t301 + m(7) * (-t234 * t56 + t237 * t57) + m(6) * (t145 * t339 - t146 * t338) + ((-t372 + t373) * t421 + (-t370 - t377) * t420 + m(5) * (-t297 * t420 + t298 * t421) + m(4) * (-t309 * t421 - t310 * t420)) * t424;
t540 = qJD(1) * t19;
t537 = t136 * t234;
t447 = t192 * t611 + t194 * t609;
t433 = (-t234 * t609 + t674) * mrSges(7,3) + t447;
t21 = t433 - t669;
t536 = t21 * qJD(1);
t490 = t598 * t300;
t512 = t423 * t426;
t22 = m(7) * (t342 * t56 + t343 * t57) + t343 * t192 + t342 * t194 + t301 * t512 + (-t490 + m(6) * (t145 * t423 - t146 * t598) - m(5) * t297 - t377) * t426 + (-m(5) * t334 + m(6) * t295 + m(7) * t221 - t357 + t673) * t515;
t535 = t22 * qJD(1);
t534 = t234 * t425;
t533 = t237 * t422;
t530 = t453 * t425;
t529 = t480 * t422;
t528 = t293 * t338;
t523 = t339 * t152;
t522 = t339 * t365;
t520 = t365 * t129;
t503 = t652 / 0.2e1;
t502 = t655 + t657;
t501 = t598 / 0.2e1;
t498 = mrSges(7,3) * t646;
t497 = t561 / 0.2e1;
t496 = -t558 / 0.2e1;
t495 = -t554 / 0.2e1;
t479 = t504 - mrSges(6,1) / 0.2e1;
t477 = t338 * t501;
t473 = -mrSges(6,2) * t444 + t549;
t469 = -Ifges(6,1) * t444 - t581;
t23 = t719 + (t706 / 0.2e1 - t154 / 0.2e1) * t453 + (t643 + t689 / 0.2e1) * t480;
t430 = (t118 / 0.4e1 + t690 / 0.4e1) * t480 - (t498 + t156 / 0.4e1 + t689 / 0.4e1) * t237 - (mrSges(7,3) * t645 + t716) * t234 + t426 * t507 / 0.4e1 + t724;
t5 = t430 - t476;
t462 = t5 * qJD(1) + t23 * qJD(2);
t31 = (-t453 ^ 2 - t480 ^ 2) * mrSges(7,3) + (-t365 ^ 2 - t444 ^ 2) * mrSges(6,3) + m(7) * (-t136 * t480 + t453 * t661) + m(6) * (t293 * t365 + t294 * t444) + (mrSges(4,3) + mrSges(5,2) + 0.4e1 * (m(5) / 0.4e1 + m(4) / 0.4e1) * qJ(3)) * (t420 ^ 2 + t421 ^ 2);
t429 = (t373 / 0.2e1 - t372 / 0.2e1) * t420 + (t377 / 0.2e1 + t370 / 0.2e1) * t421 + (-t234 * t699 + t237 * t632) * mrSges(7,3) + (-t338 * t616 - t522 / 0.2e1) * mrSges(6,3) + (-t309 * t420 + t310 * t421) * t659 + (t297 * t421 + t298 * t420) * t658 + (t145 * t365 + t146 * t444 + t293 * t339 - t294 * t338) * t656 + (t136 * t237 - t234 * t661 + t453 * t56 - t480 * t57) * t654 + t194 * t687 + t192 * t711 - t444 * t625 + t301 * t614;
t434 = -m(5) * t335 / 0.2e1 + t296 * t657 + t222 * t655 - t567 / 0.2e1 + t564 / 0.2e1 - t553 / 0.2e1 + t552 / 0.2e1;
t9 = t429 + (-m(4) * pkin(7) / 0.2e1 - t663 * t421 + (-mrSges(4,1) / 0.2e1 - mrSges(5,1) / 0.2e1) * t420) * t426 + t434;
t461 = -qJD(1) * t9 - qJD(2) * t31;
t460 = t502 * t515;
t383 = -mrSges(5,1) * t421 - mrSges(5,3) * t420;
t428 = (-t420 * t334 + (-t379 * t424 - t541) * t421) * t658 + (t420 * t295 + (t293 * t423 - t294 * t598) * t426) * t656 + (t136 * t343 + t221 * t420 + t342 * t661) * t654 + t357 * t608 + t342 * t496 + t343 * t497 - t512 * t548 / 0.2e1 - t478 * t585 + t673 * t607 + (t356 * t656 + t287 * t654 - t383 / 0.2e1 + t672 / 0.2e1) * t515;
t431 = t304 * t658 + (t149 * t598 + t423 * t150) * t656 + (-t445 * t62 + t446 * t63) * t654 + t195 * t611 + t193 * t610 + t423 * t302 / 0.2e1 - t547 / 0.2e1 + t303 * t501;
t11 = t428 - t397 - t431;
t67 = (-m(5) * t379 + m(6) * t356 - t383 + t672 + t709) * t420;
t459 = -qJD(1) * t11 - qJD(2) * t67;
t432 = (t234 * t610 + t674) * mrSges(7,3) + (-t708 * t445 + t446 * t721) * t654 + t490 / 0.2e1 + t447;
t16 = -mrSges(6,3) * t477 + (-t301 / 0.2e1 + t495 - t546 / 0.2e1) * t423 + t432 + t718;
t458 = -t16 * qJD(1) - t30 * qJD(2);
t41 = t640 * t723 + t641 * t713 + (t622 + t534 / 0.2e1 - t533 / 0.2e1) * t652 - t474;
t43 = (t615 - t530 / 0.2e1 + t529 / 0.2e1) * t652 - t473 + t76;
t457 = qJD(1) * t41 + qJD(2) * t43;
t438 = (-t423 * t338 + t339 * t598) * t656 + (t234 * t445 + t237 * t446) * t654;
t68 = -m(5) * t515 - t438 + t460;
t437 = (t365 * t598 + t423 * t444) * t656 + (-t445 * t453 - t446 * t480) * t654;
t74 = (m(5) - t502) * t420 + t437;
t456 = qJD(1) * t68 - qJD(2) * t74;
t69 = t448 + t511;
t455 = qJD(1) * t69 - qJD(2) * t76;
t454 = t221 * t365 + t287 * t339;
t285 = -Ifges(6,2) * t444 + t581;
t10 = t154 * t699 + t285 * t615 + t719 + t356 * t473 + t706 * t687 + t469 * t614 + t480 * t643 + t670 * t616 + t632 * t689 + (t152 + t709) * pkin(5) * t365;
t427 = t293 * t625 + t294 * t301 / 0.2e1 - t295 * t473 / 0.2e1 + t285 * t620 - t356 * t474 / 0.2e1 + t229 * t613 - t670 * t338 / 0.4e1 + t671 * t617 - t668 * t426 / 0.4e1 - t704 * t480 / 0.4e1 + (t689 + t156) * t237 / 0.4e1 + t716 * t234 - t724;
t435 = Ifges(6,3) * t603 + t149 * t651 - t150 * mrSges(6,2) / 0.2e1 + t451 + t476;
t439 = t195 * t601 + (t422 * t63 + t425 * t62) * t654 + t422 * t193 / 0.2e1;
t442 = t136 * t721 + t708 * t661;
t2 = t435 + t427 + Ifges(6,4) * t522 / 0.2e1 + t442 * t655 + (-t520 / 0.2e1 - t523 / 0.2e1 + t454 * t655 + t439) * pkin(5) + (t294 * t621 + t528 / 0.2e1) * mrSges(6,3) + (t544 / 0.2e1 + t65 * t711 + t545 / 0.2e1 + t64 * t687 - t537 / 0.2e1 + t661 * t683) * mrSges(7,3) + (-t338 * t613 + t339 * t617) * Ifges(6,1);
t443 = -t2 * qJD(1) + t10 * qJD(2) + t30 * qJD(4);
t436 = (t192 * t601 + t194 * t605 + (-t234 * t605 + t425 * t640) * mrSges(7,3)) * pkin(5) - t486;
t14 = (-t56 / 0.2e1 + t65 / 0.2e1) * mrSges(7,2) + (-t57 / 0.2e1 - t64 / 0.2e1) * mrSges(7,1) + t436 + t486;
t27 = (t646 + t722) * mrSges(7,2) + (t645 + t136 / 0.2e1) * mrSges(7,1);
t378 = (mrSges(7,1) * t422 + mrSges(7,2) * t425) * pkin(5);
t441 = -t14 * qJD(1) - t27 * qJD(2) + t378 * qJD(5);
t366 = t378 * qJD(6);
t75 = t437 + (m(6) + m(7)) * t608;
t71 = t460 + t438;
t44 = (-t529 + t530) * t503 + t549 / 0.2e1 + t479 * t365 + t77;
t40 = (t533 - t534) * t503 + t555 / 0.2e1 + t479 * t339 + t70;
t20 = t433 + t669;
t15 = t301 * t604 + t512 * t651 + (t339 * t604 - t477) * mrSges(6,3) + t432 - t718;
t13 = -t592 / 0.2e1 - t593 / 0.2e1 + t589 / 0.2e1 - t588 / 0.2e1 + t436 - t486;
t12 = t428 + t431;
t8 = t595 * t659 + t429 - t434 + t663 * t514 + (mrSges(4,1) + mrSges(5,1)) * t517 / 0.2e1;
t4 = t430 + t476;
t1 = t435 - t427 + t294 * t495 - mrSges(6,3) * t528 / 0.2e1 + t65 * t497 - mrSges(7,3) * t545 / 0.2e1 + t64 * t496 - t237 * t498 + t469 * t620 + t470 * t613 + (pkin(5) * t454 + t442) * t654 + t439 * pkin(5) + (t520 + t523) * pkin(5) / 0.2e1 + (-t544 + t537) * mrSges(7,3) / 0.2e1;
t17 = [qJD(2) * t3 + qJD(3) * t19 + qJD(4) * t22 + qJD(5) * t6 + qJD(6) * t7, t556 + (t664 * mrSges(5,2) + t665 * mrSges(4,3) + (t333 + t332) * t607 + (pkin(7) * mrSges(3,2) - Ifges(3,6)) * t424 + (t682 * t420 + t681 * t421) * t424 / 0.2e1 - t149 * t548 - t62 * t558 + t119 * t687 + (Ifges(6,5) * t365 + Ifges(7,5) * t453 - Ifges(6,6) * t444 + Ifges(7,6) * t480) * t603 + t661 * t195 - t421 * t330 / 0.2e1 + t335 * t383 + t379 * t358 + t356 * t245 - pkin(2) * t359 + t293 * t303 - t296 * t284 + t294 * t302 + t287 * t130 - t222 * t152 + t136 * t193 + t63 * t561 - t150 * t585 + t331 * t606 + t232 * t614 + t230 * t616 + t286 * t618 + t285 * t619 + t117 * t632 + t156 * t635 + t154 * t638) * qJD(2) + t8 * qJD(3) + t12 * qJD(4) + t1 * qJD(5) + t4 * qJD(6) + ((t371 + t376) * t421 + (-t374 + t375) * t420) * qJD(2) * qJ(3) + ((Ifges(4,2) * t421 + t584) * t608 + (-Ifges(5,3) * t421 + t579) * t607 + Ifges(3,5) + (-mrSges(4,1) * t421 + mrSges(4,2) * t420 - mrSges(3,1)) * pkin(7) + (-t578 + t583 + (Ifges(4,1) + Ifges(5,1)) * t420) * t606) * qJD(2) * t426 + 0.2e1 * ((-pkin(2) * t595 + qJ(3) * t665) * t659 + (qJ(3) * t664 + t335 * t379) * t658 + (t149 * t293 + t150 * t294 - t296 * t356) * t656 + (t136 * t63 - t222 * t287 + t62 * t661) * t654) * qJD(2), qJD(2) * t8 + qJD(4) * t71 + qJD(5) * t40 + qJD(6) * t70 + t540, t535 + t12 * qJD(2) + t71 * qJD(3) + m(7) * (-t342 * t445 + t343 * t446) * qJD(4) + t15 * qJD(5) + t20 * qJD(6), t543 + t1 * qJD(2) + t40 * qJD(3) + t15 * qJD(4) + (-t146 * mrSges(6,1) - t145 * mrSges(6,2) - t588 + t589 + t666) * qJD(5) + t13 * qJD(6) + (m(7) * (t422 * t65 + t425 * t64) + (t234 * t422 + t237 * t425) * mrSges(7,3)) * t574, t542 + t4 * qJD(2) + t70 * qJD(3) + t20 * qJD(4) + t13 * qJD(5) + (t510 - t592 - t593) * qJD(6); qJD(3) * t9 + qJD(4) * t11 - qJD(5) * t2 + qJD(6) * t5 - t556, qJD(3) * t31 + qJD(4) * t67 + qJD(5) * t10 + qJD(6) * t23, qJD(4) * t75 + qJD(5) * t44 + qJD(6) * t77 - t461, qJD(3) * t75 - t459 + t586, t44 * qJD(3) + (-t294 * mrSges(6,1) - t293 * mrSges(6,2) + t668 + t717) * qJD(5) + t725 + (m(7) * (-t136 * t425 + t422 * t661) + (-t422 * t453 - t425 * t480) * mrSges(7,3)) * t574 + t443, t77 * qJD(3) + t26 * qJD(5) + t462 + t725; -qJD(2) * t9 + qJD(4) * t68 + qJD(5) * t41 - qJD(6) * t69 - t540, -qJD(4) * t74 + qJD(5) * t43 + qJD(6) * t76 + t461, 0, t456, t457, -t455; -qJD(2) * t11 - qJD(3) * t68 + qJD(5) * t16 + qJD(6) * t21 - t535, qJD(3) * t74 + t459 + t586, -t456, 0 (-t423 * mrSges(6,1) - t598 * mrSges(6,2) + t484 + (-t422 * t445 - t425 * t446) * t652) * qJD(5) + t701 - t458, qJD(5) * t484 + t536 + t701; qJD(2) * t2 - qJD(3) * t41 - qJD(4) * t16 + qJD(6) * t14 - t543, -qJD(3) * t43 + qJD(6) * t27 - t443, -t457, t458, -t366, -t366 - t441; -qJD(2) * t5 + qJD(3) * t69 - qJD(4) * t21 - qJD(5) * t14 - t542, -qJD(3) * t76 - qJD(5) * t27 - t462, t455, -t536, t441, 0;];
Cq  = t17;
