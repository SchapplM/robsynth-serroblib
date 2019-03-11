% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 05:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRRP1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:55:42
% EndTime: 2019-03-09 05:55:59
% DurationCPUTime: 10.64s
% Computational Cost: add. (18006->580), mult. (34930->743), div. (0->0), fcn. (35356->8), ass. (0->347)
t397 = sin(qJ(5));
t394 = t397 ^ 2;
t398 = cos(qJ(5));
t395 = t398 ^ 2;
t687 = t394 + t395;
t614 = -t397 / 0.2e1;
t612 = -t398 / 0.2e1;
t679 = Ifges(7,4) + Ifges(6,5);
t686 = Ifges(6,6) - Ifges(7,6);
t609 = cos(qJ(4));
t685 = t687 * t609;
t613 = t397 / 0.2e1;
t654 = mrSges(7,2) + mrSges(6,3);
t432 = t654 * (t395 / 0.2e1 + t394 / 0.2e1);
t389 = Ifges(7,5) * t397;
t590 = Ifges(7,3) * t398;
t366 = t389 - t590;
t581 = t398 * Ifges(7,5);
t370 = Ifges(7,1) * t397 - t581;
t392 = Ifges(6,4) * t398;
t372 = Ifges(6,1) * t397 + t392;
t684 = t366 * t614 + (t370 + t372) * t612;
t607 = sin(qJ(4));
t608 = sin(qJ(3));
t610 = cos(qJ(3));
t353 = t607 * t608 - t609 * t610;
t354 = -t607 * t610 - t608 * t609;
t367 = Ifges(7,3) * t397 + t581;
t171 = -Ifges(7,6) * t354 - t353 * t367;
t466 = Ifges(6,2) * t397 - t392;
t173 = -Ifges(6,6) * t354 + t353 * t466;
t371 = Ifges(7,1) * t398 + t389;
t175 = -Ifges(7,4) * t354 - t353 * t371;
t593 = Ifges(6,4) * t397;
t373 = Ifges(6,1) * t398 - t593;
t177 = -Ifges(6,5) * t354 - t353 * t373;
t611 = t398 / 0.2e1;
t620 = -t354 / 0.2e1;
t368 = Ifges(6,2) * t398 + t593;
t543 = t397 * t368;
t659 = t543 / 0.2e1 + t684;
t665 = t679 * t397 + t686 * t398;
t429 = Ifges(5,6) * t354 + t171 * t612 + t173 * t611 + t665 * t620 + (t175 + t177) * t613 + (-Ifges(5,5) + t659) * t353;
t361 = -t398 * mrSges(7,1) - t397 * mrSges(7,3);
t559 = t353 * t397;
t330 = pkin(5) * t559;
t558 = t353 * t398;
t469 = -qJ(6) * t558 + t330;
t470 = sin(pkin(10)) * pkin(1) + pkin(7);
t444 = t610 * t470;
t340 = pkin(8) * t610 + t444;
t443 = t608 * t470;
t431 = -pkin(8) * t608 - t443;
t645 = t609 * t340 + t607 * t431;
t648 = -t469 + t645;
t667 = t648 * t361;
t362 = -t398 * mrSges(6,1) + t397 * mrSges(6,2);
t668 = t645 * t362;
t676 = t645 * mrSges(5,1);
t237 = t340 * t607 - t609 * t431;
t677 = t237 * mrSges(5,2);
t683 = t429 + t667 + t668 + t677 - t676;
t682 = -t667 / 0.2e1 - t668 / 0.2e1 + t676 / 0.2e1 - t677 / 0.2e1;
t259 = mrSges(7,2) * t559 - mrSges(7,3) * t354;
t681 = -t259 / 0.2e1;
t549 = t397 * qJ(6);
t463 = -t398 * pkin(5) - t549;
t244 = t463 * t354;
t605 = m(7) * t244;
t680 = pkin(4) * t645;
t678 = Ifges(6,3) + Ifges(7,2);
t363 = t397 * pkin(5) - qJ(6) * t398;
t113 = -t354 * t363 + t237;
t675 = t113 * t648;
t674 = t237 * t607;
t565 = t237 * t645;
t359 = -pkin(4) + t463;
t516 = t609 * pkin(3);
t338 = -t516 + t359;
t673 = t338 * t648;
t672 = t359 * t648;
t385 = -t516 - pkin(4);
t671 = t385 * t645;
t670 = t397 * t237;
t669 = t398 * t237;
t556 = t354 * t397;
t513 = mrSges(6,3) * t556;
t253 = -t353 * mrSges(6,2) + t513;
t585 = t353 * mrSges(7,3);
t258 = mrSges(7,2) * t556 + t585;
t649 = t253 + t258;
t664 = t685 * pkin(3);
t604 = m(7) * t363;
t364 = t397 * mrSges(7,1) - mrSges(7,3) * t398;
t618 = t364 / 0.2e1;
t663 = t604 / 0.2e1 + t618;
t662 = t687 * t353;
t540 = t398 * t645;
t498 = -cos(pkin(10)) * pkin(1) - pkin(2);
t358 = -pkin(3) * t610 + t498;
t219 = t353 * pkin(4) + t354 * pkin(9) + t358;
t547 = t397 * t219;
t96 = t540 + t547;
t578 = t398 * t96;
t95 = t219 * t398 - t397 * t645;
t661 = t397 * t95 - t578 + t645;
t562 = t353 * qJ(6);
t71 = t96 + t562;
t73 = -t353 * pkin(5) - t95;
t462 = -t397 * t73 - t398 * t71;
t660 = t648 + t462;
t245 = t353 * t364;
t596 = mrSges(6,2) * t398;
t365 = t397 * mrSges(6,1) + t596;
t246 = t353 * t365;
t247 = t364 * t354;
t248 = t365 * t354;
t252 = mrSges(6,2) * t354 + mrSges(6,3) * t559;
t254 = -t354 * mrSges(6,1) + mrSges(6,3) * t558;
t508 = mrSges(7,2) * t558;
t584 = t354 * mrSges(7,1);
t255 = -t508 + t584;
t260 = -t354 * mrSges(5,1) - t353 * mrSges(5,2);
t658 = -t113 * t245 - t237 * t246 - t648 * t247 - t645 * t248 + t96 * t252 + t95 * t254 + t73 * t255 + t71 * t259 + t358 * t260;
t657 = (t371 + t373) * t613 + (t367 + t466) * t612;
t606 = m(7) * qJ(6);
t656 = mrSges(6,1) + mrSges(7,1);
t655 = -mrSges(6,2) + mrSges(7,3);
t653 = t71 - t96;
t652 = t73 + t95;
t505 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t651 = t505 * t353;
t650 = t252 + t259;
t530 = t338 + t359;
t647 = t362 + t361;
t646 = Ifges(7,6) * t397 + t398 * t679;
t487 = t361 * t611;
t644 = (t487 - mrSges(7,1) / 0.2e1) * t354 + t508;
t642 = qJ(6) * t681 + pkin(5) * t255 / 0.2e1;
t514 = t607 * pkin(3);
t384 = t514 + pkin(9);
t552 = t384 * t398;
t641 = t650 * t552;
t555 = t354 * t398;
t256 = t353 * mrSges(6,1) + mrSges(6,3) * t555;
t257 = -t353 * mrSges(7,1) - mrSges(7,2) * t555;
t635 = -m(7) / 0.2e1;
t408 = t256 * t612 + t257 * t611 + (t653 * t397 - t652 * t398) * t635 + t649 * t614;
t492 = t558 / 0.2e1;
t493 = -t558 / 0.2e1;
t640 = t584 / 0.2e1 + t354 * t487 + (t493 + t492) * mrSges(7,2);
t346 = -t543 / 0.2e1;
t639 = t346 + t657 - t684;
t638 = 2 * qJD(4);
t637 = m(6) / 0.2e1;
t636 = m(6) / 0.4e1;
t634 = m(7) / 0.2e1;
t633 = m(7) / 0.4e1;
t632 = -mrSges(6,1) / 0.2e1;
t631 = mrSges(6,1) / 0.2e1;
t630 = mrSges(7,1) / 0.2e1;
t629 = mrSges(6,2) / 0.2e1;
t628 = -mrSges(7,3) / 0.2e1;
t627 = Ifges(7,6) / 0.2e1;
t626 = t96 / 0.2e1;
t598 = t354 * pkin(4);
t599 = t353 * pkin(9);
t261 = -t598 + t599;
t515 = t608 * pkin(3);
t249 = t515 + t261;
t105 = t249 * t398 + t670;
t597 = t354 * pkin(5);
t81 = -t105 + t597;
t625 = m(7) * t81;
t109 = t261 * t398 + t670;
t90 = -t109 + t597;
t624 = m(7) * t90;
t623 = -t338 / 0.2e1;
t622 = -t353 / 0.2e1;
t619 = -t359 / 0.2e1;
t617 = t365 / 0.2e1;
t603 = pkin(4) * t246;
t602 = pkin(4) * t365;
t600 = pkin(9) * t398;
t591 = Ifges(6,6) * t353;
t583 = t397 * t81;
t582 = t397 * t90;
t387 = t398 * mrSges(7,2);
t106 = t397 * t249 - t669;
t339 = t354 * qJ(6);
t80 = -t339 + t106;
t580 = t398 * t80;
t110 = t397 * t261 - t669;
t88 = -t339 + t110;
t579 = t398 * t88;
t33 = m(7) * (t113 * t555 + t71 * t353) + t353 * t258 - t247 * t555;
t576 = qJD(1) * t33;
t575 = t105 * t397;
t574 = t106 * t398;
t573 = t109 * t397;
t424 = t248 / 0.2e1 + t247 / 0.2e1 + (t681 - t252 / 0.2e1) * t398 + (t254 / 0.2e1 - t255 / 0.2e1) * t397;
t482 = -t258 / 0.2e1 - t253 / 0.2e1;
t483 = t256 / 0.2e1 - t257 / 0.2e1;
t425 = -t246 / 0.2e1 - t245 / 0.2e1 + t482 * t398 + t483 * t397;
t451 = t574 - t575;
t461 = t580 + t583;
t11 = t424 * t354 + t425 * t353 + 0.2e1 * ((-t237 - t451) * t636 + (-t113 - t461) * t633) * t354 + 0.2e1 * (t633 * t660 + t636 * t661) * t353;
t572 = t11 * qJD(1);
t571 = t110 * t398;
t450 = t571 - t573;
t460 = t579 + t582;
t12 = ((-t237 - t450) * t637 + (-t113 - t460) * t634 + t424) * t354 + (t634 * t660 + t637 * t661 + t425) * t353;
t568 = t12 * qJD(1);
t13 = (t354 * t432 + t408) * t354 + (t647 * t354 + t605) * t622;
t567 = t13 * qJD(1);
t563 = t338 * t245;
t263 = t338 * t354;
t557 = t354 * t359;
t289 = t354 * t361;
t290 = t354 * t362;
t554 = t359 * t245;
t37 = 0.4e1 * (t636 + t633) * (-0.1e1 + t687) * t354 * t353;
t553 = t37 * qJD(2);
t551 = t385 * t246;
t550 = t385 * t365;
t548 = t397 * t113;
t546 = t397 * t254;
t545 = t397 * t255;
t471 = (t631 + t630) * t397;
t52 = t330 * t634 + (-t365 / 0.2e1 - t364 / 0.2e1 - t604 / 0.2e1 + t471 + (t628 + t629 - t606 / 0.2e1) * t398) * t353;
t535 = t52 * qJD(2);
t532 = t662 * t384;
t531 = t662 * pkin(9);
t526 = qJD(5) * t354;
t525 = qJD(3) + qJD(4);
t524 = mrSges(7,2) * t580;
t523 = mrSges(7,2) * t583;
t522 = mrSges(7,2) * t579;
t521 = mrSges(7,2) * t582;
t240 = m(7) * t559;
t53 = mrSges(6,2) * t492 + mrSges(7,3) * t493 + t469 * t634 + (t471 + t617 + t663) * t353;
t520 = t53 * qJD(5) - t240 * qJD(6) + t553;
t519 = -t52 * qJD(5) - t553;
t518 = m(6) * t607;
t517 = m(7) * t607;
t512 = mrSges(6,3) * t575;
t511 = mrSges(6,3) * t574;
t510 = mrSges(6,3) * t573;
t509 = mrSges(6,3) * t571;
t504 = Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t503 = t627 - Ifges(6,6) / 0.2e1;
t502 = t95 / 0.2e1 + t73 / 0.2e1;
t501 = t626 - t71 / 0.2e1;
t500 = t384 * t546;
t499 = t384 * t545;
t495 = t397 * t609;
t494 = t398 * t609;
t490 = t247 * t614;
t176 = Ifges(7,4) * t353 - t354 * t371;
t178 = Ifges(6,5) * t353 - t354 * t373;
t484 = t178 / 0.4e1 + t176 / 0.4e1;
t480 = t364 + t604;
t333 = t354 * t516;
t331 = t353 * t514;
t477 = -t516 / 0.2e1;
t476 = t516 / 0.2e1;
t328 = Ifges(7,5) * t555;
t172 = Ifges(7,6) * t353 - Ifges(7,3) * t556 - t328;
t174 = t354 * t466 + t591;
t468 = t172 / 0.4e1 - t174 / 0.4e1 - t328 / 0.4e1;
t458 = t397 * t477;
t457 = t398 * t476;
t456 = -t333 * t687 + t331 - t532;
t406 = -Ifges(5,4) * t354 + (-t177 / 0.2e1 - t175 / 0.2e1 + t505 * t354) * t398 + (t173 / 0.2e1 - t171 / 0.2e1 + t503 * t354) * t397 + (-Ifges(5,2) + Ifges(5,1) - t678) * t353;
t410 = Ifges(5,4) * t353 + (-t178 / 0.2e1 - t176 / 0.2e1 - t651) * t398 + (t174 / 0.2e1 - t172 / 0.2e1 - t503 * t353) * t397;
t442 = mrSges(4,1) * t608 + mrSges(4,2) * t610;
t3 = (mrSges(5,1) * t515 + t410) * t353 + m(7) * (t71 * t80 + t73 * t81 + t675) + m(6) * (t105 * t95 + t106 * t96 + t565) + m(5) * t358 * t515 + (-mrSges(5,2) * t515 + t406) * t354 + t106 * t253 + t105 * t256 + t81 * t257 + t80 * t258 + t498 * t442 + (-Ifges(4,2) + Ifges(4,1)) * t610 * t608 + (-t608 ^ 2 + t610 ^ 2) * Ifges(4,4) + t658;
t455 = t3 * qJD(1) + t11 * qJD(2);
t4 = t410 * t353 + m(7) * (t71 * t88 + t73 * t90 + t675) + m(6) * (t109 * t95 + t110 * t96 + t565) + t406 * t354 + t110 * t253 + t109 * t256 + t90 * t257 + t88 * t258 + t658;
t454 = t4 * qJD(1) + t12 * qJD(2);
t295 = t363 * t361;
t435 = t295 + t639;
t44 = t338 * t480 + t435 + t550;
t418 = (Ifges(7,1) / 0.4e1 + Ifges(6,1) / 0.4e1) * t397 + (Ifges(6,4) / 0.2e1 - Ifges(7,5) / 0.4e1) * t398 - t367 / 0.4e1 - t466 / 0.4e1 + t370 / 0.4e1 + t372 / 0.4e1;
t430 = (Ifges(6,2) / 0.4e1 + Ifges(7,3) / 0.4e1) * t398 - t366 / 0.4e1 + t368 / 0.4e1 - t371 / 0.4e1 - t373 / 0.4e1;
t401 = t432 * t384 + (mrSges(7,1) * t623 + t385 * t632 + t430) * t398 + (mrSges(7,3) * t623 + t385 * t629 + t418) * t397;
t411 = t113 * t618 + t237 * t617 + t244 * t361 / 0.2e1 - t363 * t247 / 0.2e1 + t646 * t353 / 0.4e1;
t407 = t411 + t642;
t415 = (-pkin(5) * t81 + qJ(6) * t80) * t635 + t105 * t632 + t106 * t629 + t80 * t628 + t81 * t630;
t422 = (-0.3e1 / 0.4e1 * Ifges(6,6) + t627) * t353 + t501 * mrSges(7,2) + t468;
t427 = mrSges(7,2) * t502 + t484 + t651;
t440 = -t653 * t634 + t482;
t441 = t652 * t634 - t483;
t100 = t363 * t113;
t447 = (t338 * t244 + t100) * t634;
t6 = t407 + (t401 + t504) * t354 + t447 + (t384 * t441 + t427) * t398 + (t384 * t440 + t422) * t397 + t415;
t453 = t6 * qJD(1) + t44 * qJD(3);
t10 = -t237 * t290 + t244 * t247 + (t172 * t611 - mrSges(6,3) * t578 + (t346 + t372 * t611 + (-t590 + t371) * t613) * t354 + t462 * mrSges(7,2) + t665 * t622 + (t176 + t178) * t614 + (t174 + t328) * t612) * t354 + (-t289 - t605) * t113 + (-m(7) * t73 + t256 - t257) * t96 + (-m(7) * t71 + t513 - t649) * t95;
t452 = -t10 * qJD(1) - t13 * qJD(2);
t251 = (m(7) * t338 + t361) * t397;
t419 = (-t548 + (t353 * t384 + t263) * t398) * t634 - t490;
t27 = -t625 / 0.2e1 + t419 + t644;
t449 = -qJD(1) * t27 + qJD(3) * t251;
t386 = mrSges(7,3) + t606;
t42 = t585 + 0.2e1 * (t562 / 0.2e1 + t547 / 0.4e1 + t540 / 0.4e1 - t96 / 0.4e1) * m(7);
t448 = qJD(1) * t42 + qJD(5) * t386;
t446 = (t359 * t244 + t100) * t634;
t304 = t385 * t354;
t409 = (-t304 + t456) * t637 + (-t263 + t456) * t634;
t426 = (-t531 + t598) * t637 + (-t531 - t557) * t634;
t25 = t409 - t426;
t399 = (t671 + t450 * t384 + (t494 * t96 - t495 * t95 + t674) * pkin(3)) * t637 + (t673 + t460 * t384 + (t113 * t607 + t494 * t71 + t495 * t73) * pkin(3)) * t634 - t563 / 0.2e1 - t551 / 0.2e1 - t510 / 0.2e1 + t509 / 0.2e1 - t500 / 0.2e1 + t499 / 0.2e1 + t522 / 0.2e1 + t521 / 0.2e1 + t256 * t458 + t397 * t257 * t476 + t641 / 0.2e1 + (-t247 - t248) * t514 / 0.2e1 + t649 * t457 - t682;
t400 = -m(6) * (pkin(9) * t451 - t680) / 0.2e1 + (pkin(9) * t461 + t672) * t635 - t603 / 0.2e1 + t554 / 0.2e1 + t512 / 0.2e1 - t511 / 0.2e1 + pkin(9) * t546 / 0.2e1 - pkin(9) * t545 / 0.2e1 - t524 / 0.2e1 - t523 / 0.2e1 - t650 * t600 / 0.2e1 + t682;
t5 = t400 + t399;
t405 = -t609 * mrSges(5,2) + t654 * t685 + (-mrSges(5,1) + t647) * t607;
t56 = (t338 * t517 + t385 * t518 + t405) * pkin(3) + (m(7) + m(6)) * t664 * t384;
t438 = t5 * qJD(1) + t25 * qJD(2) + t56 * qJD(3);
t404 = (-pkin(5) * t495 + qJ(6) * t494) * pkin(3) * t634 + t477 * t596 + mrSges(7,3) * t457 + t656 * t458;
t413 = t295 - t602 / 0.2e1 + t550 / 0.2e1 + t663 * t530;
t31 = t404 - t413 - t657 + t659;
t54 = t359 * t480 + t435 - t602;
t403 = t432 * pkin(9) + (mrSges(7,1) * t619 + pkin(4) * t631 + t430) * t398 + (-pkin(4) * mrSges(6,2) / 0.2e1 + mrSges(7,3) * t619 + t418) * t397;
t414 = (-pkin(5) * t90 + qJ(6) * t88) * t635 + t109 * t632 + t110 * t629 + t88 * t628 + t90 * t630;
t8 = t407 + (t403 + t504) * t354 + t446 + (pkin(9) * t441 + t427) * t398 + (pkin(9) * t440 + t422) * t397 + t414;
t437 = t8 * qJD(1) - t31 * qJD(3) + t54 * qJD(4);
t133 = (t361 + (t476 + t338 / 0.2e1 + t359 / 0.2e1) * m(7)) * t397;
t262 = (m(7) * t359 + t361) * t397;
t421 = (-t548 + (t557 + t599) * t398) * t634 - t490;
t30 = -t624 / 0.2e1 + t421 + t644;
t436 = -qJD(1) * t30 + qJD(3) * t133 + qJD(4) * t262;
t434 = -t654 * t662 - t260 - t289 - t290;
t420 = (-mrSges(7,2) * t549 - Ifges(6,6) * t397 - pkin(5) * t387 + t646) * qJD(5) - t535;
t416 = qJD(5) * (m(7) * t463 + t647);
t402 = t484 * t398 + (t397 * t501 + t398 * t502) * mrSges(7,2) + (-t591 / 0.4e1 + t468) * t397 + t411 - t642 + t678 * t620 + (-Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t559 + t679 * t493;
t360 = m(7) * t600 + t387;
t337 = m(7) * t552 + t387;
t134 = m(7) * t530 * t614 + (m(7) * t476 - t361) * t397;
t38 = (t96 + 0.2e1 * t562) * t634 + m(7) * t626 + t258;
t32 = t404 + t413 + t639;
t29 = t624 / 0.2e1 + t421 + t640;
t26 = t625 / 0.2e1 + t419 + t640;
t17 = t409 + t426 + t434;
t9 = pkin(9) * t408 + t354 * t403 + t402 - t414 + t446;
t7 = t354 * t401 + t384 * t408 + t402 - t415 + t447;
t2 = -t400 + t399 + t429;
t1 = qJD(3) * t11 + qJD(4) * t12 - qJD(5) * t13;
t14 = [qJD(3) * t3 + qJD(4) * t4 - qJD(5) * t10 + qJD(6) * t33, t1 (t683 + t641 + t524 + t499 - t563 + t523 - t500 - Ifges(4,6) * t608 + Ifges(4,5) * t610 + t511 + mrSges(4,2) * t443 - t512 + m(7) * (t384 * t461 + t673) + m(6) * (t384 * t451 + t671) + m(5) * (-t609 * t645 - t674) * pkin(3) - t551 + (t353 * t516 + t354 * t514) * mrSges(5,3) - mrSges(4,1) * t444) * qJD(3) + t2 * qJD(4) + t7 * qJD(5) + t26 * qJD(6) + t455, t2 * qJD(3) + t9 * qJD(5) + t29 * qJD(6) + ((pkin(9) * t450 - t680) * t637 + (pkin(9) * t460 + t672) * t634) * t638 + t454 + (t509 - t510 + t521 + t522 - t554 + t603 + (t650 * t398 + (-t254 + t255) * t397) * pkin(9) + t683) * qJD(4), t7 * qJD(3) + t9 * qJD(4) + ((-m(7) * pkin(5) - t656) * t96 + (t655 + t606) * t95) * qJD(5) + t38 * qJD(6) + ((mrSges(7,2) * qJ(6) + t686) * t398 + (-mrSges(7,2) * pkin(5) + t679) * t397) * t526 + t452, qJD(3) * t26 + qJD(4) * t29 + qJD(5) * t38 + t576; t1, t525 * t37, t17 * qJD(4) + t520 + t572 + (t434 - t442 + 0.2e1 * (-t304 - t532) * t637 + 0.2e1 * (-t263 - t532) * t634 + m(5) * (-t331 + t333)) * qJD(3), t17 * qJD(3) + qJD(4) * t434 + t426 * t638 + t520 + t568, -t567 - qJD(5) * t605 + t525 * t53 + (t655 * qJD(5) * t397 + (-m(7) * qJD(6) + t656 * qJD(5)) * t398) * t354, -m(7) * t398 * t526 - t240 * t525; qJD(4) * t5 + qJD(5) * t6 + qJD(6) * t27 - t455, qJD(4) * t25 + t519 - t572, qJD(4) * t56 + qJD(5) * t44 - qJD(6) * t251, t32 * qJD(5) + t134 * qJD(6) + (-pkin(4) * t518 + t359 * t517 + t405) * qJD(4) * pkin(3) + t438 + (t637 + t634) * t638 * t664 * pkin(9), t32 * qJD(4) + t337 * qJD(6) + t384 * t416 + t420 + t453, qJD(4) * t134 + qJD(5) * t337 - t449; -qJD(3) * t5 + qJD(5) * t8 + qJD(6) * t30 - t454, -qJD(3) * t25 + t519 - t568, -qJD(5) * t31 - qJD(6) * t133 - t438, qJD(5) * t54 - qJD(6) * t262, pkin(9) * t416 + t360 * qJD(6) + t420 + t437, qJD(5) * t360 - t436; -qJD(3) * t6 - qJD(4) * t8 + qJD(6) * t42 - t452, t52 * t525 + t567, qJD(4) * t31 - t453 + t535, -t437 + t535, t386 * qJD(6), t448; -qJD(3) * t27 - qJD(4) * t30 - qJD(5) * t42 - t576, 0, qJD(4) * t133 + t449, t436, -t448, 0;];
Cq  = t14;
