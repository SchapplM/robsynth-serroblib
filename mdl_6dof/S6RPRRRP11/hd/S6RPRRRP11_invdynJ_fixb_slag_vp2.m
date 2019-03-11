% Calculate vector of inverse dynamics joint torques for
% S6RPRRRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRP11_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP11_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP11_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP11_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP11_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP11_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP11_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:34:31
% EndTime: 2019-03-09 06:35:44
% DurationCPUTime: 44.90s
% Computational Cost: add. (30363->930), mult. (96730->1251), div. (0->0), fcn. (83741->14), ass. (0->418)
t347 = sin(pkin(6));
t350 = cos(pkin(6));
t345 = sin(pkin(12));
t349 = cos(pkin(7));
t354 = sin(qJ(3));
t348 = cos(pkin(12));
t554 = cos(qJ(3));
t470 = t554 * t348;
t381 = -t345 * t354 + t349 * t470;
t346 = sin(pkin(7));
t477 = t346 * t554;
t364 = t347 * t381 + t350 * t477;
t255 = t364 * qJD(1);
t250 = -t255 + qJD(4);
t553 = sin(qJ(1));
t469 = t553 * t348;
t555 = cos(qJ(1));
t473 = t555 * t345;
t301 = t350 * t473 + t469;
t341 = t553 * t345;
t472 = t555 * t348;
t396 = -t350 * t472 + t341;
t476 = t347 * t555;
t682 = t346 * t476 + t349 * t396;
t229 = -t301 * t554 + t354 * t682;
t353 = sin(qJ(4));
t356 = cos(qJ(4));
t366 = t396 * t346 - t349 * t476;
t187 = t229 * t356 - t353 * t366;
t226 = t301 * t354 + t554 * t682;
t352 = sin(qJ(5));
t355 = cos(qJ(5));
t689 = t187 * t352 + t226 * t355;
t688 = t187 * t355 - t226 * t352;
t654 = Ifges(6,4) + Ifges(7,4);
t506 = t349 * t354;
t374 = t347 * (-t345 * t506 + t470);
t282 = qJD(1) * t374;
t457 = qJD(3) * t554;
t431 = t346 * t457;
t683 = t431 - t282;
t655 = Ifges(6,1) + Ifges(7,1);
t653 = -Ifges(7,5) - Ifges(6,5);
t652 = Ifges(6,2) + Ifges(7,2);
t651 = Ifges(6,6) + Ifges(7,6);
t471 = t554 * t345;
t383 = t348 * t506 + t471;
t489 = qJD(1) * qJD(3);
t510 = t346 * t350;
t202 = (qJD(1) * t457 + qJDD(1) * t354) * t510 + (qJDD(1) * t383 + t381 * t489) * t347;
t509 = t346 * t354;
t268 = t347 * t383 + t350 * t509;
t258 = t268 * qJD(1);
t508 = t347 * t348;
t392 = t346 * t508 - t349 * t350;
t290 = -qJD(1) * t392 + qJD(3);
t215 = -t353 * t258 + t290 * t356;
t289 = -qJDD(1) * t392 + qJDD(3);
t127 = qJD(4) * t215 + t202 * t356 + t353 * t289;
t580 = t127 / 0.2e1;
t687 = Ifges(5,4) * t580;
t498 = qJD(1) * t347;
t467 = t348 * t498;
t552 = pkin(1) * t350;
t486 = qJD(1) * t552;
t295 = qJ(2) * t467 + t345 * t486;
t507 = t347 * t349;
t380 = (t348 * t507 + t510) * pkin(9);
t245 = qJD(1) * t380 + t295;
t549 = pkin(9) * t345;
t284 = (-pkin(2) * t348 - t346 * t549 - pkin(1)) * t347;
t274 = qJD(1) * t284 + qJD(2);
t332 = t348 * t486;
t512 = t345 * t347;
t551 = pkin(2) * t350;
t375 = t551 + (-pkin(9) * t349 - qJ(2)) * t512;
t251 = qJD(1) * t375 + t332;
t516 = t251 * t349;
t163 = t554 * t245 + (t274 * t346 + t516) * t354;
t684 = -t163 + t250 * (pkin(4) * t353 - pkin(11) * t356);
t216 = t258 * t356 + t353 * t290;
t128 = -qJD(4) * t216 - t353 * t202 + t289 * t356;
t123 = qJDD(5) - t128;
t171 = t216 * t355 + t250 * t352;
t312 = (qJ(2) * qJDD(1) + qJD(1) * qJD(2)) * t347;
t487 = qJDD(1) * t350;
t481 = pkin(1) * t487;
t278 = -t312 * t345 + t348 * t481;
t238 = (-t507 * t549 + t551) * qJDD(1) + t278;
t270 = qJDD(1) * t284 + qJDD(2);
t192 = -t238 * t346 + t349 * t270;
t203 = t347 * (qJDD(1) * t381 - t383 * t489) - (-qJDD(1) * t554 + t354 * t489) * t510;
t107 = -pkin(3) * t203 - pkin(10) * t202 + t192;
t209 = -t251 * t346 + t349 * t274;
t139 = -pkin(3) * t255 - pkin(10) * t258 + t209;
t142 = t290 * pkin(10) + t163;
t494 = qJD(4) * t356;
t495 = qJD(4) * t353;
t279 = t348 * t312 + t345 * t481;
t235 = qJDD(1) * t380 + t279;
t474 = t349 * t554;
t439 = t251 * t474;
t496 = qJD(3) * t354;
t91 = qJD(3) * t439 + t554 * t235 + t238 * t506 - t245 * t496 + t270 * t509 + t274 * t431;
t86 = pkin(10) * t289 + t91;
t19 = t353 * t107 + t139 * t494 - t142 * t495 + t356 * t86;
t200 = qJDD(4) - t203;
t11 = pkin(11) * t200 + t19;
t465 = t346 * t496;
t92 = -t354 * t235 + t238 * t474 - t245 * t457 + t270 * t477 - t274 * t465 - t496 * t516;
t87 = -t289 * pkin(3) - t92;
t31 = -t128 * pkin(4) - t127 * pkin(11) + t87;
t74 = t353 * t139 + t142 * t356;
t67 = t250 * pkin(11) + t74;
t162 = -t354 * t245 + t274 * t477 + t439;
t141 = -t290 * pkin(3) - t162;
t90 = -t215 * pkin(4) - t216 * pkin(11) + t141;
t33 = t352 * t90 + t355 * t67;
t4 = -qJD(5) * t33 - t11 * t352 + t355 * t31;
t170 = -t216 * t352 + t250 * t355;
t63 = qJD(5) * t170 + t127 * t355 + t200 * t352;
t1 = pkin(5) * t123 - qJ(6) * t63 - qJD(6) * t171 + t4;
t569 = t200 / 0.2e1;
t579 = t128 / 0.2e1;
t582 = t123 / 0.2e1;
t64 = -qJD(5) * t171 - t127 * t352 + t200 * t355;
t585 = t64 / 0.2e1;
t586 = t63 / 0.2e1;
t650 = -Ifges(7,3) - Ifges(6,3);
t491 = qJD(5) * t355;
t493 = qJD(5) * t352;
t3 = t355 * t11 + t352 * t31 + t90 * t491 - t493 * t67;
t2 = qJ(6) * t64 + qJD(6) * t170 + t3;
t662 = t4 * mrSges(6,1) - t3 * mrSges(6,2) - t2 * mrSges(7,2);
t681 = t1 * mrSges(7,1) - 0.2e1 * Ifges(5,2) * t579 - 0.2e1 * Ifges(5,6) * t569 - t650 * t582 + t651 * t585 - t653 * t586 + t662 - t687;
t680 = t229 * t353 + t356 * t366;
t679 = t654 * t585 + t655 * t586;
t343 = pkin(5) * t355 + pkin(4);
t678 = m(7) * t343;
t677 = t654 * t170;
t504 = t352 * t356;
t193 = -t255 * t504 + t355 * t258;
t464 = t352 * t494;
t626 = t353 * t491 + t193 + t464;
t307 = -t356 * t349 + t353 * t509;
t468 = t345 * t498;
t433 = t346 * t468;
t622 = -qJD(4) * t307 - t353 * t433 + t356 * t683;
t373 = t347 * (t348 * t354 + t349 * t471);
t281 = qJD(1) * t373;
t676 = t465 - t281;
t520 = t215 * t352;
t675 = t493 - t520;
t674 = t654 * t171;
t351 = -qJ(6) - pkin(11);
t673 = -m(7) * t351 + mrSges(6,3) + mrSges(7,3);
t214 = qJD(5) - t215;
t645 = t170 * t651 - t171 * t653 - t214 * t650;
t672 = t645 / 0.2e1;
t671 = -t653 * t123 / 0.2e1 + t679;
t648 = t123 * t651 + t63 * t654 + t64 * t652;
t670 = t648 / 0.2e1;
t644 = t170 * t652 + t214 * t651 + t674;
t643 = t171 * t655 - t653 * t214 + t677;
t483 = pkin(10) * t494;
t204 = pkin(3) * t258 - pkin(10) * t255;
t112 = -t353 * t162 + t204 * t356;
t99 = -t258 * pkin(4) - t112;
t669 = t483 - t99;
t113 = t356 * t162 + t353 * t204;
t100 = pkin(11) * t258 + t113;
t484 = pkin(10) * t495;
t668 = t684 * t355 + (t100 + t484) * t352;
t425 = -pkin(4) * t356 - t353 * pkin(11);
t323 = -pkin(3) + t425;
t667 = -t355 * t100 + t323 * t491 + t352 * t684;
t541 = mrSges(4,3) * t258;
t666 = -mrSges(4,1) * t290 - mrSges(5,1) * t215 + mrSges(5,2) * t216 + t541;
t308 = t353 * t349 + t356 * t509;
t621 = qJD(4) * t308 + t353 * t683 + t356 * t433;
t376 = t350 * t469 + t473;
t475 = t347 * t553;
t362 = t376 * t346 + t349 * t475;
t665 = t654 * t355;
t664 = t654 * t352;
t663 = Ifges(5,1) * t580 + Ifges(5,5) * t569;
t587 = Ifges(5,4) * t579 + t663;
t568 = -t214 / 0.2e1;
t576 = -t171 / 0.2e1;
t578 = -t170 / 0.2e1;
t661 = -t650 * t568 - t653 * t576 + t651 * t578;
t588 = m(7) * pkin(5);
t73 = t139 * t356 - t353 * t142;
t659 = t73 * mrSges(5,1);
t658 = t74 * mrSges(5,2);
t657 = -mrSges(6,1) - mrSges(7,1);
t656 = mrSges(6,2) + mrSges(7,2);
t649 = -t123 * t650 - t63 * t653 + t64 * t651;
t24 = -mrSges(6,1) * t64 + mrSges(6,2) * t63;
t95 = mrSges(5,1) * t200 - mrSges(5,3) * t127;
t646 = t24 - t95;
t642 = Ifges(4,5) * t202;
t641 = Ifges(4,6) * t203;
t640 = Ifges(4,3) * t289;
t502 = t355 * t356;
t194 = t255 * t502 + t352 * t258;
t340 = pkin(10) * t502;
t490 = qJD(6) * t355;
t515 = t255 * t353;
t639 = -pkin(5) * t515 + qJ(6) * t194 - t353 * t490 + (pkin(5) * t353 - qJ(6) * t502) * qJD(4) + (-t340 + (qJ(6) * t353 - t323) * t352) * qJD(5) + t668;
t503 = t353 * t355;
t638 = -qJ(6) * t193 + (-pkin(10) * qJD(4) - qJ(6) * qJD(5)) * t503 + (-qJD(6) * t353 + (-pkin(10) * qJD(5) - qJ(6) * qJD(4)) * t356) * t352 + t667;
t637 = (-t355 * t495 - t356 * t493) * pkin(10) + t667;
t292 = t352 * t323 + t340;
t636 = -qJD(5) * t292 + t668;
t635 = pkin(5) * t626 + t669;
t443 = qJD(5) * t351;
t152 = pkin(4) * t216 - pkin(11) * t215;
t53 = t352 * t152 + t355 * t73;
t634 = qJ(6) * t520 + t352 * t443 + t490 - t53;
t519 = t215 * t355;
t52 = t355 * t152 - t352 * t73;
t633 = -pkin(5) * t216 + qJ(6) * t519 - qJD(6) * t352 + t355 * t443 - t52;
t632 = pkin(5) * t675 - t74;
t631 = t588 + mrSges(7,1);
t304 = qJ(2) * t508 + t345 * t552;
t259 = t380 + t304;
t339 = t348 * t552;
t269 = t339 + t375;
t174 = -t354 * t259 + t269 * t474 + t284 * t477;
t160 = pkin(3) * t392 - t174;
t224 = t268 * t353 + t356 * t392;
t225 = t268 * t356 - t353 * t392;
t104 = t224 * pkin(4) - t225 * pkin(11) + t160;
t217 = -t269 * t346 + t349 * t284;
t444 = pkin(3) * t364 + t268 * pkin(10);
t154 = t217 - t444;
t248 = t554 * t259;
t175 = t269 * t506 + t284 * t509 + t248;
t161 = -pkin(10) * t392 + t175;
t94 = t353 * t154 + t356 * t161;
t77 = -pkin(11) * t364 + t94;
t41 = t352 * t104 + t355 * t77;
t109 = -mrSges(6,1) * t170 + mrSges(6,2) * t171;
t539 = mrSges(5,3) * t216;
t173 = mrSges(5,1) * t250 - t539;
t627 = t173 - t109;
t492 = qJD(5) * t353;
t625 = t352 * t492 - t355 * t494 + t194;
t620 = -t346 * t475 + t376 * t349;
t415 = -mrSges(7,1) * t355 + mrSges(7,2) * t352;
t418 = -mrSges(6,1) * t355 + mrSges(6,2) * t352;
t619 = m(6) * pkin(4) - t415 - t418 + t678;
t414 = mrSges(7,1) * t352 + mrSges(7,2) * t355;
t417 = mrSges(6,1) * t352 + mrSges(6,2) * t355;
t66 = -t250 * pkin(4) - t73;
t51 = -t170 * pkin(5) + qJD(6) + t66;
t618 = t51 * t414 + t66 * t417;
t617 = -t352 * t653 + t355 * t651;
t616 = -t352 * t651 - t355 * t653;
t615 = t355 * t652 + t664;
t614 = -t352 * t652 + t665;
t613 = t352 * t655 + t665;
t612 = t355 * t655 - t664;
t611 = -m(6) * pkin(11) - t673;
t609 = -t491 + t519;
t607 = t495 - t515;
t20 = t107 * t356 - t139 * t495 - t142 * t494 - t353 * t86;
t606 = t19 * t356 - t20 * t353;
t605 = t3 * t355 - t352 * t4;
t604 = -m(7) - m(5) - m(6);
t603 = -mrSges(6,1) - t631;
t602 = mrSges(5,1) + t619;
t589 = mrSges(5,2) + t611;
t420 = -mrSges(5,1) * t356 + t353 * mrSges(5,2);
t600 = -m(6) * t425 + t353 * t673 + t356 * t678 - t420;
t597 = -g(1) * t475 + g(2) * t476 - g(3) * t350;
t596 = mrSges(4,1) + t600;
t478 = pkin(5) * t352 + pkin(10);
t595 = -m(7) * t478 + mrSges(4,2) - mrSges(5,3);
t594 = -m(5) * t141 - t666;
t593 = -t92 * mrSges(4,1) + t91 * mrSges(4,2);
t592 = t20 * mrSges(5,1) - t19 * mrSges(5,2);
t534 = Ifges(5,4) * t216;
t125 = t215 * Ifges(5,2) + t250 * Ifges(5,6) + t534;
t581 = -t125 / 0.2e1;
t577 = t170 / 0.2e1;
t575 = t171 / 0.2e1;
t567 = t214 / 0.2e1;
t566 = -t215 / 0.2e1;
t565 = -t216 / 0.2e1;
t564 = t216 / 0.2e1;
t562 = -t250 / 0.2e1;
t559 = t258 / 0.2e1;
t550 = pkin(5) * t171;
t182 = -t225 * t352 - t355 * t364;
t546 = t182 * pkin(5);
t540 = mrSges(5,3) * t215;
t538 = mrSges(6,3) * t170;
t537 = mrSges(6,3) * t171;
t536 = mrSges(7,3) * t170;
t535 = mrSges(7,3) * t171;
t533 = Ifges(5,4) * t353;
t532 = Ifges(5,4) * t356;
t12 = -t200 * pkin(4) - t20;
t527 = t12 * t353;
t526 = t162 * mrSges(4,3);
t521 = t258 * Ifges(4,4);
t518 = t229 * t352;
t302 = -t341 * t350 + t472;
t231 = t302 * t554 - t354 * t620;
t517 = t231 * t352;
t514 = t255 * t356;
t513 = t268 * t352;
t505 = t352 * t353;
t499 = t555 * pkin(1) + qJ(2) * t475;
t497 = qJD(2) * t347;
t488 = qJDD(1) * t347;
t480 = Ifges(5,5) * t127 + Ifges(5,6) * t128 + Ifges(5,3) * t200;
t479 = t640 + t641 + t642;
t466 = t345 * t497;
t23 = -t64 * mrSges(7,1) + t63 * mrSges(7,2);
t456 = t345 * t488;
t455 = t348 * t488;
t453 = t494 / 0.2e1;
t448 = -t492 / 0.2e1;
t32 = -t352 * t67 + t355 * t90;
t40 = t355 * t104 - t352 * t77;
t93 = t154 * t356 - t353 * t161;
t432 = t346 * t466;
t426 = -pkin(1) * t553 + qJ(2) * t476;
t422 = -mrSges(4,1) * t364 + mrSges(4,2) * t268;
t421 = mrSges(5,1) * t224 + mrSges(5,2) * t225;
t183 = t225 * t355 - t352 * t364;
t419 = mrSges(6,1) * t182 - mrSges(6,2) * t183;
t416 = -t182 * mrSges(7,1) + t183 * mrSges(7,2);
t413 = Ifges(5,1) * t356 - t533;
t408 = -Ifges(5,2) * t353 + t532;
t403 = Ifges(5,5) * t356 - Ifges(5,6) * t353;
t189 = t231 * t356 + t353 * t362;
t230 = t302 * t354 + t554 * t620;
t134 = -t189 * t352 + t230 * t355;
t397 = -(-qJ(2) * t468 + t332) * t345 + t295 * t348;
t147 = qJD(2) * t374 + qJD(3) * t174;
t256 = t364 * qJD(3);
t257 = t268 * qJD(3);
t191 = pkin(3) * t257 - pkin(10) * t256 + t432;
t48 = -t353 * t147 - t154 * t495 - t161 * t494 + t191 * t356;
t395 = -mrSges(3,1) * t455 + mrSges(3,2) * t456;
t394 = mrSges(3,1) * t350 - mrSges(3,3) * t512;
t393 = -mrSges(3,2) * t350 + mrSges(3,3) * t508;
t76 = pkin(4) * t364 - t93;
t25 = -qJ(6) * t171 + t32;
t47 = t356 * t147 + t154 * t494 - t161 * t495 + t353 * t191;
t44 = pkin(11) * t257 + t47;
t148 = qJD(2) * t373 + (t248 + (t269 * t349 + t284 * t346) * t354) * qJD(3);
t180 = qJD(4) * t225 + t256 * t353;
t181 = -qJD(4) * t224 + t256 * t356;
t70 = t180 * pkin(4) - t181 * pkin(11) + t148;
t7 = t104 * t491 + t352 * t70 + t355 * t44 - t493 * t77;
t275 = -t352 * t308 - t355 * t477;
t382 = -t355 * t308 + t352 * t477;
t45 = -t257 * pkin(4) - t48;
t8 = -qJD(5) * t41 - t352 * t44 + t355 * t70;
t363 = -t301 * pkin(2) - pkin(9) * t366 + t426;
t361 = t302 * pkin(2) + pkin(9) * t362 + t499;
t360 = t229 * pkin(3) + t363;
t359 = t231 * pkin(3) + t361;
t358 = -pkin(10) * t226 + t360;
t357 = t230 * pkin(10) + t359;
t333 = -pkin(1) * t488 + qJDD(2);
t325 = t351 * t355;
t324 = t351 * t352;
t320 = t478 * t353;
t315 = t355 * t323;
t306 = t393 * qJD(1);
t305 = t394 * qJD(1);
t303 = -qJ(2) * t512 + t339;
t291 = -pkin(10) * t504 + t315;
t277 = -qJ(6) * t505 + t292;
t266 = -qJ(6) * t503 + t315 + (-pkin(10) * t352 - pkin(5)) * t356;
t249 = Ifges(4,4) * t255;
t218 = -mrSges(4,2) * t290 + mrSges(4,3) * t255;
t213 = Ifges(5,4) * t215;
t201 = -mrSges(4,1) * t255 + mrSges(4,2) * t258;
t188 = t231 * t353 - t356 * t362;
t179 = t258 * Ifges(4,1) + t290 * Ifges(4,5) + t249;
t178 = t255 * Ifges(4,2) + t290 * Ifges(4,6) + t521;
t177 = -mrSges(4,2) * t289 + mrSges(4,3) * t203;
t176 = mrSges(4,1) * t289 - mrSges(4,3) * t202;
t172 = -mrSges(5,2) * t250 + t540;
t135 = t189 * t355 + t230 * t352;
t129 = -mrSges(4,1) * t203 + mrSges(4,2) * t202;
t126 = t216 * Ifges(5,1) + t250 * Ifges(5,5) + t213;
t124 = t216 * Ifges(5,5) + t215 * Ifges(5,6) + t250 * Ifges(5,3);
t118 = mrSges(6,1) * t214 - t537;
t117 = mrSges(7,1) * t214 - t535;
t116 = -mrSges(6,2) * t214 + t538;
t115 = -mrSges(7,2) * t214 + t536;
t108 = -mrSges(7,1) * t170 + mrSges(7,2) * t171;
t106 = qJD(5) * t182 + t181 * t355 + t257 * t352;
t105 = -qJD(5) * t183 - t181 * t352 + t257 * t355;
t96 = -mrSges(5,2) * t200 + mrSges(5,3) * t128;
t65 = -mrSges(5,1) * t128 + mrSges(5,2) * t127;
t56 = t76 - t546;
t39 = -mrSges(6,2) * t123 + mrSges(6,3) * t64;
t38 = -mrSges(7,2) * t123 + mrSges(7,3) * t64;
t37 = mrSges(6,1) * t123 - mrSges(6,3) * t63;
t36 = mrSges(7,1) * t123 - mrSges(7,3) * t63;
t30 = qJ(6) * t182 + t41;
t27 = pkin(5) * t224 - qJ(6) * t183 + t40;
t26 = qJ(6) * t170 + t33;
t22 = -t105 * pkin(5) + t45;
t21 = pkin(5) * t214 + t25;
t9 = -t64 * pkin(5) + qJDD(6) + t12;
t6 = qJ(6) * t105 + qJD(6) * t182 + t7;
t5 = pkin(5) * t180 - qJ(6) * t106 - qJD(6) * t183 + t8;
t10 = [m(4) * (t147 * t163 + t174 * t92 + t175 * t91 + t192 * t217 + t209 * t432) - (Ifges(5,6) * t579 + Ifges(5,5) * t580 - Ifges(4,6) * t289 - t91 * mrSges(4,3) - Ifges(4,2) * t203 - Ifges(4,4) * t202 + Ifges(5,3) * t569 + t480 / 0.2e1 + t592) * t364 + t180 * t581 + m(5) * (t160 * t87 + t19 * t94 + t20 * t93 + t47 * t74 + t48 * t73) + t278 * t394 + t279 * t393 + t257 * t659 + m(3) * (t278 * t303 + t279 * t304) + (-pkin(1) * t395 + t333 * (-mrSges(3,1) * t348 + mrSges(3,2) * t345) + m(3) * (-pkin(1) * t333 + qJD(2) * t397) + (Ifges(3,4) * t345 + Ifges(3,2) * t348) * t455 + (Ifges(3,1) * t345 + Ifges(3,4) * t348) * t456 + (Ifges(3,5) * t345 + Ifges(3,6) * t348) * t487) * t347 + (-mrSges(4,3) * t92 + Ifges(4,1) * t202 + Ifges(4,4) * t203 + Ifges(4,5) * t289) * t268 + t180 * t672 + t201 * t432 + t643 * t106 / 0.2e1 + t644 * t105 / 0.2e1 + (-mrSges(5,3) * t20 + 0.2e1 * t587) * t225 + (-m(4) * t162 - t594) * t148 + (-t640 / 0.2e1 - t641 / 0.2e1 - t642 / 0.2e1 - t479 / 0.2e1 + t593) * t392 + t87 * t421 + t192 * t422 + t9 * t416 - t12 * t419 + (Ifges(3,5) * t456 + Ifges(3,6) * t455 + Ifges(3,3) * t487) * t350 + m(7) * (t1 * t27 + t2 * t30 + t21 * t5 + t22 * t51 + t26 * t6 + t56 * t9) + m(6) * (t12 * t76 + t3 * t41 + t32 * t8 + t33 * t7 + t4 * t40 + t45 * t66) + (-m(6) * (t187 * pkin(4) + t358) - m(5) * t358 - t187 * mrSges(5,1) - m(4) * t363 - t229 * mrSges(4,1) + mrSges(4,3) * t366 + mrSges(2,1) * t553 + mrSges(2,2) * t555 - m(3) * t426 + t301 * mrSges(3,1) - mrSges(3,2) * t396 - mrSges(3,3) * t476 - m(7) * (t187 * t343 + t360) - t595 * t226 + t657 * t688 + t656 * t689 + t589 * t680) * g(1) + (t105 * t651 - t106 * t653 - t180 * t650) * t567 + (t105 * t652 + t106 * t654 + t180 * t651) * t577 + (t105 * t654 + t106 * t655 - t180 * t653) * t575 + (-t180 * t74 - t181 * t73) * mrSges(5,3) + (Ifges(4,1) * t256 - Ifges(4,4) * t257) * t559 + (Ifges(5,1) * t181 - Ifges(5,4) * t180 + Ifges(5,5) * t257) * t564 + (t649 / 0.2e1 - t687 - t19 * mrSges(5,3) + t681) * t224 + (t303 * t394 + t304 * t393 + Ifges(2,3)) * qJDD(1) + (-t4 * mrSges(6,3) - t1 * mrSges(7,3) - t653 * t582 + t671 + t679) * t183 - t256 * t526 + t290 * (Ifges(4,5) * t256 - Ifges(4,6) * t257) / 0.2e1 + t256 * t179 / 0.2e1 + t250 * (Ifges(5,5) * t181 - Ifges(5,6) * t180 + Ifges(5,3) * t257) / 0.2e1 + t215 * (Ifges(5,4) * t181 - Ifges(5,2) * t180 + Ifges(5,6) * t257) / 0.2e1 + t257 * t124 / 0.2e1 + t255 * (Ifges(4,4) * t256 - Ifges(4,2) * t257) / 0.2e1 + t209 * (mrSges(4,1) * t257 + mrSges(4,2) * t256) - t257 * t178 / 0.2e1 + t217 * t129 + t147 * t218 + t141 * (mrSges(5,1) * t180 + mrSges(5,2) * t181) + t181 * t126 / 0.2e1 + t33 * (-mrSges(6,2) * t180 + mrSges(6,3) * t105) - t305 * t466 + t21 * (mrSges(7,1) * t180 - mrSges(7,3) * t106) + t32 * (mrSges(6,1) * t180 - mrSges(6,3) * t106) + t26 * (-mrSges(7,2) * t180 + mrSges(7,3) * t105) + t47 * t172 + t48 * t173 + t174 * t176 + t175 * t177 + t160 * t65 + t5 * t117 + t8 * t118 + t6 * t115 + t7 * t116 + t22 * t108 + t45 * t109 + t51 * (-mrSges(7,1) * t105 + mrSges(7,2) * t106) + t66 * (-mrSges(6,1) * t105 + mrSges(6,2) * t106) + t93 * t95 + t94 * t96 + (-m(6) * (t189 * pkin(4) + t357) - m(5) * t357 - t189 * mrSges(5,1) - m(4) * t361 - t231 * mrSges(4,1) - mrSges(4,3) * t362 - mrSges(2,1) * t555 + mrSges(2,2) * t553 - m(3) * t499 - t302 * mrSges(3,1) + mrSges(3,2) * t376 - mrSges(3,3) * t475 - m(7) * (t189 * t343 + t359) + t595 * t230 + t657 * t135 - t656 * t134 + t589 * t188) * g(2) - t257 * t658 - t163 * t257 * mrSges(4,3) + t27 * t36 + t30 * t38 + t40 * t37 + t41 * t39 + (t3 * mrSges(6,3) + t2 * mrSges(7,3) + t651 * t582 + t652 * t585 + t654 * t586 + t670) * t182 + t348 * t306 * t497 + t56 * t23 + t76 * t24; (t36 + t37) * t275 - t201 * t433 + (-t141 * t281 + t19 * t308 - t20 * t307 + (t141 * t496 - t554 * t87) * t346 + t622 * t74 - t621 * t73 + t597) * m(5) + t622 * t172 + (t23 + t646) * t307 + (t176 - t65) * t477 + t683 * t218 - (t38 + t39) * t382 + t395 + (t162 * t281 - t163 * t282 - t209 * t433 + t192 * t349 + (t554 * t92 + t354 * t91 + (-t162 * t354 + t163 * t554) * qJD(3)) * t346 + t597) * m(4) + (-t397 * t498 + t333 + t597) * m(3) + (t12 * t307 + t275 * t4 - t3 * t382 + t621 * t66 + t597) * m(6) + (t1 * t275 - t2 * t382 + t307 * t9 + t51 * t621 + t597) * m(7) + t305 * t468 + t349 * t129 + t308 * t96 - t306 * t467 + t177 * t509 + (m(6) * t32 + m(7) * t21 + t117 + t118) * (qJD(5) * t382 - t352 * t622 + t355 * t676) + (m(6) * t33 + m(7) * t26 + t115 + t116) * (qJD(5) * t275 + t352 * t676 + t355 * t622) + t666 * t676 + t621 * (t108 - t627); (-g(1) * t231 + g(2) * t229 - g(3) * t268 - t607 * t74 + (-t494 + t514) * t73 + t606) * mrSges(5,3) + (t518 * t588 - mrSges(4,2) * t229 + t604 * (-t226 * pkin(3) - pkin(10) * t229) + t657 * (-t226 * t502 - t518) - t656 * (t226 * t504 - t229 * t355) + t596 * t226) * g(2) + (-t513 * t588 + t422 + t604 * t444 + t657 * (t364 * t502 + t513) - t656 * (t268 * t355 - t364 * t504) - t600 * t364) * g(3) + (-t484 - t113) * t172 + (t193 * t651 - t194 * t653) * t568 + (t193 * t652 + t194 * t654) * t578 + (t193 * t654 + t194 * t655) * t576 - t593 + (t9 * t414 + t582 * t616 + t585 * t614 + t586 * t612 + t587 + t663) * t353 + t646 * pkin(10) * t353 + (-t1 * t503 - t2 * t505) * mrSges(7,3) + (t215 * t408 + t216 * t413 + t250 * t403) * qJD(4) / 0.2e1 - (Ifges(4,1) * t255 + t124 - t521) * t258 / 0.2e1 - (-Ifges(4,2) * t258 + t179 + t249) * t255 / 0.2e1 + t533 * t579 + (-t112 * t73 - t113 * t74 - t87 * pkin(3) + ((-t353 * t74 - t356 * t73) * qJD(4) + t606) * pkin(10)) * m(5) + t250 * t141 * (mrSges(5,1) * t353 + mrSges(5,2) * t356) + (-t483 - t112) * t173 + t258 * t658 + t503 * t671 - t648 * t505 / 0.2e1 - t649 * t356 / 0.2e1 + t636 * t118 + t637 * t116 + (-t66 * t99 + t4 * t291 + t3 * t292 + (t494 * t66 + t527) * pkin(10) + t637 * t33 + t636 * t32) * m(6) + t638 * t115 + (t1 * t266 + t2 * t277 + t21 * t639 + t26 * t638 + t320 * t9 + t51 * t635) * m(7) + t639 * t117 + (-t645 / 0.2e1 + t125 / 0.2e1 + t661) * t515 + t532 * t580 + (t541 + t594) * t163 + t87 * t420 + t635 * t108 + (-t514 / 0.2e1 + t453) * t126 + (-t617 * t492 + (-t353 * t650 + t356 * t616) * qJD(4)) * t567 + (-t615 * t492 + (t353 * t651 + t356 * t614) * qJD(4)) * t577 + (-t613 * t492 + (-t353 * t653 + t356 * t612) * qJD(4)) * t575 + (mrSges(6,1) * t607 + mrSges(6,3) * t625) * t32 + (mrSges(7,1) * t607 + mrSges(7,3) * t625) * t21 + (mrSges(7,1) * t626 - mrSges(7,2) * t625) * t51 + (t626 * mrSges(6,1) - t625 * mrSges(6,2)) * t66 + (-mrSges(7,2) * t607 - mrSges(7,3) * t626) * t26 + (-mrSges(6,2) * t607 - mrSges(6,3) * t626) * t33 + t178 * t559 + (Ifges(5,3) * t258 + t255 * t403) * t562 + (Ifges(5,5) * t258 + t255 * t413) * t565 + (Ifges(5,6) * t258 + t255 * t408) * t566 + (pkin(10) * t96 - t681) * t356 + t320 * t23 - t290 * (Ifges(4,5) * t255 - Ifges(4,6) * t258) / 0.2e1 + t291 * t37 + t292 * t39 + t277 * t38 + t266 * t36 - t209 * (mrSges(4,1) * t258 + mrSges(4,2) * t255) - t162 * t218 - t258 * t659 + t669 * t109 + (-t517 * t588 + mrSges(4,2) * t231 + t604 * (-t230 * pkin(3) + t231 * pkin(10)) + t657 * (-t230 * t502 + t517) - t656 * (t230 * t504 + t231 * t355) + t596 * t230) * g(1) + t479 + t255 * t526 + t417 * t527 + t643 * (t352 * t448 + t355 * t453 - t194 / 0.2e1) + t644 * (t355 * t448 - t464 / 0.2e1 - t193 / 0.2e1) + (t581 + t672) * t495 - pkin(3) * t65 + (-t3 * t505 - t4 * t503) * mrSges(6,3); t613 * t586 + t615 * t585 + t617 * t582 + t592 + (-t141 * mrSges(5,2) + Ifges(5,1) * t565 + Ifges(5,5) * t562 + t568 * t616 + t576 * t612 + t578 * t614 - t618) * t215 + t618 * qJD(5) + (t224 * t619 + t225 * t611 + t421) * g(3) + (-t118 * t491 - t116 * t493 + m(6) * ((-t32 * t355 - t33 * t352) * qJD(5) + t605) + t355 * t39 - t352 * t37) * pkin(11) + (-t187 * t589 - t602 * t680) * g(2) + t355 * t670 + t352 * t671 + (-t534 + t645) * t565 + (t188 * t602 + t189 * t589) * g(1) + (-t141 * mrSges(5,1) - t32 * mrSges(6,1) - t21 * mrSges(7,1) + t33 * mrSges(6,2) + t26 * mrSges(7,2) - Ifges(5,2) * t566 - Ifges(5,6) * t562 + t661) * t216 + (t170 * t614 + t171 * t612 + t214 * t616) * qJD(5) / 0.2e1 + (-t493 / 0.2e1 + t520 / 0.2e1) * t644 + (t213 + t126) * t566 + t9 * t415 + t12 * t418 + t632 * t108 + t633 * t117 + t634 * t115 + (t1 * t324 - t2 * t325 + t21 * t633 + t26 * t634 - t343 * t9 + t51 * t632) * m(7) + (t491 / 0.2e1 - t519 / 0.2e1) * t643 + (-m(6) * t66 + t539 + t627) * t74 + t125 * t564 + (t32 * t609 - t33 * t675 + t605) * mrSges(6,3) + (-t1 * t352 + t2 * t355 + t21 * t609 - t26 * t675) * mrSges(7,3) + (-pkin(4) * t12 - t32 * t52 - t33 * t53) * m(6) - t343 * t23 + t324 * t36 - t325 * t38 + (t540 - t172) * t73 - t52 * t118 - t53 * t116 - pkin(4) * t24 + t480; (-m(7) * t550 - mrSges(7,1) * t171 - mrSges(7,2) * t170) * t51 + (t537 + t118) * t33 + (t538 - t116) * t32 + (t535 - m(7) * (-t21 + t25) + t117) * t26 + (t603 * t689 - t656 * t688) * g(2) + t644 * t575 + t662 + t631 * t1 + (-m(7) * t546 + t416 - t419) * g(3) + t21 * t536 - t108 * t550 + (-t171 * t652 + t643 + t677) * t578 + (t134 * t603 + t135 * t656) * g(1) + (-t170 * t653 - t171 * t651) * t568 + (t170 * t655 - t674) * t576 - t66 * (mrSges(6,1) * t171 + mrSges(6,2) * t170) - t25 * t115 + pkin(5) * t36 + t649; -t170 * t115 + t171 * t117 + (-g(1) * t188 + g(2) * t680 - g(3) * t224 - t170 * t26 + t171 * t21 + t9) * m(7) + t23;];
tau  = t10;
