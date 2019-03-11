% Calculate vector of inverse dynamics joint torques for
% S6RPRRRP12
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
% Datum: 2019-03-09 06:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRP12_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP12_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP12_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP12_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP12_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP12_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP12_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:44:05
% EndTime: 2019-03-09 06:45:25
% DurationCPUTime: 47.63s
% Computational Cost: add. (30590->978), mult. (97313->1307), div. (0->0), fcn. (84256->14), ass. (0->427)
t359 = sin(pkin(6));
t362 = cos(pkin(6));
t357 = sin(pkin(12));
t361 = cos(pkin(7));
t365 = sin(qJ(3));
t360 = cos(pkin(12));
t572 = cos(qJ(3));
t485 = t572 * t360;
t393 = -t357 * t365 + t361 * t485;
t358 = sin(pkin(7));
t492 = t358 * t572;
t377 = t359 * t393 + t362 * t492;
t267 = t377 * qJD(1);
t372 = -t267 + qJD(4);
t571 = sin(qJ(1));
t484 = t571 * t360;
t573 = cos(qJ(1));
t488 = t573 * t357;
t317 = t362 * t488 + t484;
t354 = t571 * t357;
t487 = t573 * t360;
t407 = -t362 * t487 + t354;
t491 = t359 * t573;
t694 = t358 * t491 + t361 * t407;
t237 = -t317 * t572 + t365 * t694;
t284 = t407 * t358 - t361 * t491;
t364 = sin(qJ(4));
t367 = cos(qJ(4));
t185 = t237 * t367 - t284 * t364;
t234 = t317 * t365 + t572 * t694;
t363 = sin(qJ(5));
t366 = cos(qJ(5));
t699 = t185 * t363 + t234 * t366;
t698 = t185 * t366 - t234 * t363;
t522 = t361 * t365;
t385 = t359 * (-t357 * t522 + t485);
t299 = qJD(1) * t385;
t475 = qJD(3) * t572;
t446 = t358 * t475;
t695 = t446 - t299;
t666 = Ifges(6,1) + Ifges(7,1);
t664 = Ifges(6,5) + Ifges(7,4);
t513 = qJD(1) * t359;
t482 = t360 * t513;
t570 = pkin(1) * t362;
t499 = qJD(1) * t570;
t312 = qJ(2) * t482 + t357 * t499;
t523 = t359 * t361;
t526 = t358 * t362;
t392 = (t360 * t523 + t526) * pkin(9);
t256 = qJD(1) * t392 + t312;
t567 = pkin(9) * t357;
t301 = (-pkin(2) * t360 - t358 * t567 - pkin(1)) * t359;
t289 = qJD(1) * t301 + qJD(2);
t346 = t360 * t499;
t528 = t357 * t359;
t569 = pkin(2) * t362;
t386 = t569 + (-pkin(9) * t361 - qJ(2)) * t528;
t263 = qJD(1) * t386 + t346;
t531 = t263 * t361;
t157 = t572 * t256 + (t289 * t358 + t531) * t365;
t693 = -pkin(10) * qJD(5) * t367 - t157 + t372 * (pkin(4) * t364 - pkin(11) * t367);
t692 = t237 * t364 + t284 * t367;
t486 = t572 * t357;
t395 = t360 * t522 + t486;
t502 = qJD(1) * qJD(3);
t200 = (qJD(1) * t475 + qJDD(1) * t365) * t526 + (qJDD(1) * t395 + t393 * t502) * t359;
t525 = t358 * t365;
t282 = t359 * t395 + t362 * t525;
t270 = t282 * qJD(1);
t524 = t359 * t360;
t403 = t358 * t524 - t361 * t362;
t306 = -qJD(1) * t403 + qJD(3);
t213 = t270 * t367 + t364 * t306;
t305 = -qJDD(1) * t403 + qJDD(3);
t123 = -qJD(4) * t213 - t364 * t200 + t305 * t367;
t118 = qJDD(5) - t123;
t596 = t118 / 0.2e1;
t212 = -t364 * t270 + t306 * t367;
t122 = qJD(4) * t212 + t200 * t367 + t364 * t305;
t201 = t359 * (qJDD(1) * t393 - t395 * t502) - (-qJDD(1) * t572 + t365 * t502) * t526;
t198 = qJDD(4) - t201;
t504 = t366 * qJD(4);
t521 = t363 * t213;
t62 = t366 * t122 + t363 * t198 + (-t267 * t366 + t504 - t521) * qJD(5);
t604 = t62 / 0.2e1;
t691 = t664 * t596 + t666 * t604;
t165 = t366 * t213 + t363 * t372;
t63 = qJD(5) * t165 + t363 * t122 - t366 * t198;
t603 = -t63 / 0.2e1;
t602 = t63 / 0.2e1;
t665 = -Ifges(6,4) + Ifges(7,5);
t690 = Ifges(6,6) - Ifges(7,6);
t662 = Ifges(6,3) + Ifges(7,2);
t323 = -t361 * t367 + t364 * t525;
t483 = t357 * t513;
t448 = t358 * t483;
t637 = -qJD(4) * t323 - t364 * t448 + t367 * t695;
t384 = t359 * (t360 * t365 + t361 * t486);
t298 = qJD(1) * t384;
t511 = qJD(3) * t365;
t480 = t358 * t511;
t689 = t480 - t298;
t505 = t364 * qJD(4);
t518 = t364 * t267;
t688 = t505 - t518;
t211 = qJD(5) - t212;
t327 = (qJ(2) * qJDD(1) + qJD(1) * qJD(2)) * t359;
t500 = qJDD(1) * t362;
t495 = pkin(1) * t500;
t293 = -t327 * t357 + t360 * t495;
t248 = (-t523 * t567 + t569) * qJDD(1) + t293;
t285 = qJDD(1) * t301 + qJDD(2);
t190 = -t248 * t358 + t361 * t285;
t102 = -pkin(3) * t201 - pkin(10) * t200 + t190;
t207 = -t263 * t358 + t361 * t289;
t134 = -pkin(3) * t267 - pkin(10) * t270 + t207;
t137 = t306 * pkin(10) + t157;
t503 = t367 * qJD(4);
t294 = t360 * t327 + t357 * t495;
t244 = qJDD(1) * t392 + t294;
t489 = t361 * t572;
t455 = t263 * t489;
t87 = qJD(3) * t455 + t572 * t244 + t248 * t522 - t256 * t511 + t285 * t525 + t289 * t446;
t83 = pkin(10) * t305 + t87;
t20 = t364 * t102 + t134 * t503 - t137 * t505 + t367 * t83;
t12 = pkin(11) * t198 + t20;
t88 = -t365 * t244 + t248 * t489 - t256 * t475 + t285 * t492 - t289 * t480 - t511 * t531;
t84 = -t305 * pkin(3) - t88;
t27 = -t123 * pkin(4) - t122 * pkin(11) + t84;
t71 = t364 * t134 + t367 * t137;
t66 = pkin(11) * t372 + t71;
t156 = -t365 * t256 + t289 * t492 + t455;
t136 = -t306 * pkin(3) - t156;
t86 = -t212 * pkin(4) - t213 * pkin(11) + t136;
t29 = t363 * t86 + t366 * t66;
t4 = -qJD(5) * t29 - t12 * t363 + t27 * t366;
t2 = -pkin(5) * t118 + qJDD(6) - t4;
t507 = qJD(5) * t366;
t509 = qJD(5) * t363;
t3 = t366 * t12 + t363 * t27 + t86 * t507 - t509 * t66;
t674 = t3 * mrSges(6,2);
t1 = qJ(6) * t118 + qJD(6) * t211 + t3;
t675 = t1 * mrSges(7,3);
t687 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t674 + t675;
t686 = t665 * t602 + t691;
t164 = -t366 * t372 + t521;
t163 = Ifges(6,4) * t164;
t545 = Ifges(7,5) * t164;
t656 = t165 * t666 + t664 * t211 - t163 + t545;
t262 = Ifges(4,4) * t267;
t685 = Ifges(4,2) * t267;
t684 = Ifges(4,6) * t306;
t683 = Ifges(5,6) * t212;
t682 = t306 * Ifges(4,5);
t568 = pkin(4) * t367;
t339 = -t364 * pkin(11) - pkin(3) - t568;
t202 = pkin(3) * t270 - pkin(10) * t267;
t108 = t367 * t156 + t364 * t202;
t96 = pkin(11) * t270 + t108;
t648 = -t364 * t504 * pkin(10) + t339 * t507 + t363 * t693 - t366 * t96;
t647 = -t339 * t509 + (pkin(10) * t505 + t96) * t363 + t693 * t366;
t70 = t367 * t134 - t364 * t137;
t65 = -pkin(4) * t372 - t70;
t110 = -mrSges(7,2) * t164 + mrSges(7,3) * t211;
t553 = mrSges(6,3) * t164;
t111 = -mrSges(6,2) * t211 - t553;
t681 = -t110 - t111;
t552 = mrSges(6,3) * t165;
t112 = mrSges(6,1) * t211 - t552;
t113 = -mrSges(7,1) * t211 + mrSges(7,2) * t165;
t680 = t112 - t113;
t554 = mrSges(4,3) * t270;
t679 = -mrSges(4,1) * t306 - mrSges(5,1) * t212 + mrSges(5,2) * t213 + t554;
t324 = t364 * t361 + t367 * t525;
t636 = qJD(4) * t324 + t364 * t695 + t367 * t448;
t387 = t362 * t484 + t488;
t490 = t359 * t571;
t375 = t387 * t358 + t361 * t490;
t678 = Ifges(6,4) * t604 + Ifges(6,6) * t596 - Ifges(7,5) * t62 / 0.2e1 - Ifges(7,6) * t118 / 0.2e1 + (Ifges(6,2) + Ifges(7,3)) * t603;
t585 = t198 / 0.2e1;
t593 = t123 / 0.2e1;
t594 = t122 / 0.2e1;
t605 = Ifges(5,1) * t594 + Ifges(5,4) * t593 + Ifges(5,5) * t585;
t676 = -m(6) - m(7);
t673 = -t372 / 0.2e1;
t672 = t372 / 0.2e1;
t670 = t70 * mrSges(5,1);
t669 = t71 * mrSges(5,2);
t668 = mrSges(4,2) - mrSges(5,3);
t667 = mrSges(6,3) + mrSges(7,2);
t661 = t118 * t662 + t62 * t664 - t63 * t690;
t33 = -mrSges(7,2) * t63 + mrSges(7,3) * t118;
t36 = -mrSges(6,2) * t118 - mrSges(6,3) * t63;
t659 = t33 + t36;
t34 = mrSges(6,1) * t118 - mrSges(6,3) * t62;
t35 = -t118 * mrSges(7,1) + t62 * mrSges(7,2);
t658 = t35 - t34;
t657 = -t164 * t690 + t165 * t664 + t211 * t662;
t23 = mrSges(6,1) * t63 + mrSges(6,2) * t62;
t91 = mrSges(5,1) * t198 - mrSges(5,3) * t122;
t655 = -t91 + t23;
t654 = Ifges(4,5) * t200;
t653 = Ifges(4,6) * t201;
t652 = Ifges(4,3) * t305;
t651 = qJ(6) * t688 - qJD(6) * t367 + t648;
t650 = -pkin(5) * t688 - t647;
t519 = t363 * t367;
t191 = t267 * t519 - t366 * t270;
t515 = t366 * t367;
t192 = t267 * t515 + t363 * t270;
t417 = pkin(5) * t363 - qJ(6) * t366;
t408 = pkin(10) + t417;
t418 = pkin(5) * t366 + qJ(6) * t363;
t107 = -t364 * t156 + t202 * t367;
t95 = -t270 * pkin(4) - t107;
t649 = -t191 * pkin(5) + t192 * qJ(6) + (qJD(5) * t418 - qJD(6) * t366) * t364 + t408 * t503 - t95;
t646 = -qJD(6) * t363 + t211 * t417 - t71;
t438 = -mrSges(5,1) * t367 + t364 * mrSges(5,2);
t645 = -t438 + mrSges(4,1);
t105 = mrSges(6,1) * t164 + mrSges(6,2) * t165;
t536 = t213 * mrSges(5,3);
t167 = mrSges(5,1) * t372 - t536;
t641 = t167 - t105;
t320 = qJ(2) * t524 + t357 * t570;
t275 = t392 + t320;
t352 = t360 * t570;
t283 = t352 + t386;
t168 = -t365 * t275 + t283 * t489 + t301 * t492;
t28 = -t363 * t66 + t366 * t86;
t638 = -t28 + qJD(6);
t435 = -t366 * mrSges(7,1) - t363 * mrSges(7,3);
t437 = mrSges(6,1) * t366 - mrSges(6,2) * t363;
t635 = m(7) * t418 - t435 + t437;
t634 = -t358 * t490 + t387 * t361;
t32 = t164 * pkin(5) - t165 * qJ(6) + t65;
t434 = mrSges(7,1) * t363 - mrSges(7,3) * t366;
t436 = mrSges(6,1) * t363 + mrSges(6,2) * t366;
t633 = t32 * t434 + t65 * t436;
t632 = t363 * t664 + t366 * t690;
t631 = -t363 * t690 + t366 * t664;
t543 = Ifges(7,5) * t366;
t546 = Ifges(6,4) * t366;
t630 = t363 * t666 - t543 + t546;
t544 = Ifges(7,5) * t363;
t547 = Ifges(6,4) * t363;
t629 = t366 * t666 + t544 - t547;
t534 = t212 * t366;
t627 = t507 - t534;
t535 = t212 * t363;
t626 = -t509 + t535;
t21 = t102 * t367 - t134 * t505 - t137 * t503 - t364 * t83;
t625 = t20 * t367 - t21 * t364;
t624 = t3 * t366 - t363 * t4;
t623 = t1 * t366 + t2 * t363;
t540 = t165 * Ifges(6,4);
t80 = -t164 * Ifges(6,2) + t211 * Ifges(6,6) + t540;
t599 = -t80 / 0.2e1;
t162 = Ifges(7,5) * t165;
t77 = t211 * Ifges(7,6) + t164 * Ifges(7,3) + t162;
t622 = t77 / 0.2e1 + t599;
t621 = -t77 / 0.2e1 + t80 / 0.2e1;
t620 = mrSges(5,2) - t667;
t619 = mrSges(5,1) + t635;
t616 = -g(1) * t490 + g(2) * t491 - g(3) * t362;
t269 = t282 * qJD(3);
t141 = qJD(2) * t385 + qJD(3) * t168;
t214 = -t283 * t358 + t361 * t301;
t464 = pkin(3) * t377 + pkin(10) * t282;
t148 = t214 - t464;
t261 = t572 * t275;
t169 = t283 * t522 + t301 * t525 + t261;
t155 = -pkin(10) * t403 + t169;
t268 = t377 * qJD(3);
t512 = qJD(2) * t359;
t481 = t357 * t512;
t447 = t358 * t481;
t189 = pkin(3) * t269 - pkin(10) * t268 + t447;
t46 = t367 * t141 + t148 * t503 - t155 * t505 + t364 * t189;
t40 = pkin(11) * t269 + t46;
t90 = t364 * t148 + t367 * t155;
t74 = -pkin(11) * t377 + t90;
t154 = pkin(3) * t403 - t168;
t232 = t282 * t364 + t367 * t403;
t233 = t282 * t367 - t364 * t403;
t461 = -t232 * pkin(4) + t233 * pkin(11);
t99 = t154 - t461;
t555 = t363 * t99 + t366 * t74;
t142 = qJD(2) * t384 + (t261 + (t283 * t361 + t301 * t358) * t365) * qJD(3);
t174 = qJD(4) * t233 + t268 * t364;
t175 = -qJD(4) * t232 + t268 * t367;
t68 = t174 * pkin(4) - t175 * pkin(11) + t142;
t9 = -qJD(5) * t555 - t363 * t40 + t366 * t68;
t460 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t457 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t615 = -m(5) * t136 - t679;
t614 = -t88 * mrSges(4,1) + t87 * mrSges(4,2);
t613 = t21 * mrSges(5,1) - t20 * mrSges(5,2);
t611 = -mrSges(7,2) * t1 - mrSges(6,3) * t3 - t678;
t318 = -t354 * t362 + t487;
t239 = t318 * t572 - t365 * t634;
t187 = t239 * t367 + t364 * t375;
t610 = -g(1) * t187 + g(2) * t185 - g(3) * t233;
t25 = qJ(6) * t211 + t29;
t609 = -mrSges(6,1) * t65 - mrSges(7,1) * t32 + mrSges(7,2) * t25 + mrSges(6,3) * t29;
t24 = -pkin(5) * t211 + t638;
t584 = -t211 / 0.2e1;
t590 = -t165 / 0.2e1;
t591 = t164 / 0.2e1;
t592 = -t164 / 0.2e1;
t608 = -t28 * mrSges(6,1) + t24 * mrSges(7,1) + t29 * mrSges(6,2) - t25 * mrSges(7,3) + Ifges(6,6) * t591 + Ifges(7,6) * t592 + t584 * t662 + t590 * t664;
t550 = Ifges(5,4) * t213;
t120 = Ifges(5,2) * t212 + t372 * Ifges(5,6) + t550;
t595 = -t120 / 0.2e1;
t589 = t165 / 0.2e1;
t583 = t211 / 0.2e1;
t582 = -t212 / 0.2e1;
t581 = -t213 / 0.2e1;
t580 = t213 / 0.2e1;
t576 = t270 / 0.2e1;
t560 = t377 * pkin(4);
t551 = Ifges(4,4) * t270;
t549 = Ifges(5,4) * t364;
t548 = Ifges(5,4) * t367;
t13 = -t198 * pkin(4) - t21;
t542 = t13 * t364;
t541 = t156 * mrSges(4,3);
t537 = t212 * mrSges(5,3);
t146 = pkin(4) * t213 - pkin(11) * t212;
t52 = t363 * t146 + t366 * t70;
t533 = t234 * t364;
t238 = t318 * t365 + t572 * t634;
t532 = t238 * t364;
t530 = t267 * t367;
t529 = t377 * t364;
t517 = t364 * t366;
t516 = t366 * t339;
t308 = pkin(10) * t515 + t363 * t339;
t514 = t573 * pkin(1) + qJ(2) * t490;
t508 = qJD(5) * t364;
t501 = qJDD(1) * t359;
t498 = pkin(10) * t503;
t497 = pkin(11) * t509;
t496 = pkin(11) * t507;
t494 = Ifges(5,5) * t122 + Ifges(5,6) * t123 + Ifges(5,3) * t198;
t493 = t652 + t653 + t654;
t479 = t363 * t503;
t474 = t357 * t501;
t473 = t360 * t501;
t470 = -t508 / 0.2e1;
t469 = t507 / 0.2e1;
t467 = t503 / 0.2e1;
t466 = -t234 * pkin(3) - pkin(10) * t237;
t465 = -t238 * pkin(3) + pkin(10) * t239;
t89 = t148 * t367 - t364 * t155;
t453 = t363 * t492;
t443 = -pkin(1) * t571 + qJ(2) * t491;
t440 = -mrSges(4,1) * t377 + mrSges(4,2) * t282;
t439 = mrSges(5,1) * t232 + mrSges(5,2) * t233;
t433 = Ifges(5,1) * t367 - t549;
t428 = -Ifges(5,2) * t364 + t548;
t427 = -Ifges(6,2) * t363 + t546;
t426 = Ifges(6,2) * t366 + t547;
t423 = Ifges(5,5) * t367 - Ifges(5,6) * t364;
t420 = Ifges(7,3) * t363 + t543;
t419 = -Ifges(7,3) * t366 + t544;
t37 = -t363 * t74 + t366 * t99;
t51 = t146 * t366 - t363 * t70;
t177 = t233 * t366 - t363 * t377;
t176 = t233 * t363 + t366 * t377;
t412 = -(-qJ(2) * t483 + t346) * t357 + t312 * t360;
t47 = -t364 * t141 - t148 * t505 - t155 * t503 + t189 * t367;
t406 = -mrSges(3,1) * t473 + mrSges(3,2) * t474;
t405 = mrSges(3,1) * t362 - mrSges(3,3) * t528;
t404 = -mrSges(3,2) * t362 + mrSges(3,3) * t524;
t73 = -t89 + t560;
t8 = t363 * t68 + t366 * t40 + t99 * t507 - t509 * t74;
t394 = -t363 * t324 - t366 * t492;
t391 = -t363 * t508 + t366 * t503;
t390 = t364 * t507 + t479;
t41 = -t269 * pkin(4) - t47;
t376 = -t317 * pkin(2) - pkin(9) * t284 + t443;
t374 = t318 * pkin(2) + pkin(9) * t375 + t514;
t371 = t237 * pkin(3) - pkin(10) * t234 + t376;
t370 = t239 * pkin(3) + t238 * pkin(10) + t374;
t347 = -pkin(1) * t501 + qJDD(2);
t334 = -pkin(4) - t418;
t322 = t404 * qJD(1);
t321 = t405 * qJD(1);
t319 = -qJ(2) * t528 + t352;
t310 = t408 * t364;
t307 = -pkin(10) * t519 + t516;
t297 = -t516 + (pkin(10) * t363 + pkin(5)) * t367;
t296 = -qJ(6) * t367 + t308;
t291 = t366 * t324 - t453;
t215 = -mrSges(4,2) * t306 + mrSges(4,3) * t267;
t210 = Ifges(5,4) * t212;
t199 = -mrSges(4,1) * t267 + mrSges(4,2) * t270;
t186 = t239 * t364 - t367 * t375;
t173 = t270 * Ifges(4,1) + t262 + t682;
t172 = t551 + t684 + t685;
t171 = -mrSges(4,2) * t305 + mrSges(4,3) * t201;
t170 = mrSges(4,1) * t305 - mrSges(4,3) * t200;
t166 = -mrSges(5,2) * t372 + t537;
t130 = t187 * t366 + t238 * t363;
t129 = t187 * t363 - t238 * t366;
t124 = -mrSges(4,1) * t201 + mrSges(4,2) * t200;
t121 = Ifges(5,1) * t213 + t372 * Ifges(5,5) + t210;
t119 = Ifges(5,5) * t213 + Ifges(5,3) * t372 + t683;
t104 = mrSges(7,1) * t164 - mrSges(7,3) * t165;
t103 = pkin(5) * t165 + qJ(6) * t164;
t101 = -qJD(5) * t176 + t175 * t366 + t269 * t363;
t100 = qJD(5) * t177 + t175 * t363 - t269 * t366;
t92 = -mrSges(5,2) * t198 + mrSges(5,3) * t123;
t64 = -mrSges(5,1) * t123 + mrSges(5,2) * t122;
t54 = Ifges(5,4) * t122 + Ifges(5,2) * t123 + Ifges(5,6) * t198;
t48 = t176 * pkin(5) - t177 * qJ(6) + t73;
t43 = -pkin(5) * t213 - t51;
t42 = qJ(6) * t213 + t52;
t31 = -pkin(5) * t232 - t37;
t30 = qJ(6) * t232 + t555;
t22 = mrSges(7,1) * t63 - mrSges(7,3) * t62;
t10 = t100 * pkin(5) - t101 * qJ(6) - t177 * qJD(6) + t41;
t7 = -pkin(5) * t174 - t9;
t6 = qJ(6) * t174 + qJD(6) * t232 + t8;
t5 = t63 * pkin(5) - t62 * qJ(6) - t165 * qJD(6) + t13;
t11 = [(Ifges(6,4) * t101 + Ifges(6,6) * t174) * t592 + (-t493 / 0.2e1 - t654 / 0.2e1 - t652 / 0.2e1 - t653 / 0.2e1 + t614) * t403 + t656 * t101 / 0.2e1 + t657 * t174 / 0.2e1 + t555 * t36 + (Ifges(7,5) * t101 + Ifges(7,6) * t174) * t591 + (-t101 * t32 + t174 * t25) * mrSges(7,3) + (t101 * t664 + t174 * t662) * t583 + (t101 * t666 + t174 * t664) * t589 + t190 * t440 + t212 * (Ifges(5,4) * t175 - Ifges(5,2) * t174) / 0.2e1 + (Ifges(5,5) * t175 - Ifges(5,6) * t174) * t672 + (t319 * t405 + t320 * t404 + Ifges(2,3)) * qJDD(1) + (mrSges(6,1) * t13 + mrSges(7,1) * t5 - Ifges(6,2) * t603 + Ifges(7,3) * t602 - t596 * t690 + t604 * t665 + t611) * t176 + (-Ifges(6,2) * t592 + Ifges(7,3) * t591 - t583 * t690 + t589 * t665 - t609 + t622) * t100 - t321 * t481 + m(6) * (t13 * t73 + t28 * t9 + t29 * t8 + t3 * t555 + t37 * t4 + t41 * t65) + (Ifges(5,1) * t175 - Ifges(5,4) * t174) * t580 + (-t174 * t71 - t175 * t70) * mrSges(5,3) + (-m(3) * t443 - m(4) * t376 - m(5) * t371 + mrSges(2,1) * t571 + t317 * mrSges(3,1) - t237 * mrSges(4,1) - t185 * mrSges(5,1) + mrSges(2,2) * t573 - mrSges(3,2) * t407 - mrSges(3,3) * t491 + mrSges(4,3) * t284 + t676 * (t185 * pkin(4) + pkin(11) * t692 + t371) - t668 * t234 - t460 * t698 + t457 * t699 + t620 * t692) * g(1) + m(5) * (t154 * t84 + t20 * t90 + t21 * t89 + t46 * t71 + t47 * t70) + t199 * t447 + (-mrSges(5,3) * t21 + 0.2e1 * t605) * t233 + (t13 * mrSges(6,2) + t2 * mrSges(7,2) - t4 * mrSges(6,3) - t5 * mrSges(7,3) + Ifges(6,4) * t603 + Ifges(7,5) * t602 + t686 + t691) * t177 + (-m(4) * t156 - t615) * t142 + t84 * t439 + (-m(3) * t514 - m(4) * t374 - m(5) * t370 - mrSges(2,1) * t573 - t318 * mrSges(3,1) - t239 * mrSges(4,1) - t187 * mrSges(5,1) + mrSges(2,2) * t571 + mrSges(3,2) * t387 - mrSges(3,3) * t490 - mrSges(4,3) * t375 + t676 * (t187 * pkin(4) + t186 * pkin(11) + t370) + t668 * t238 - t460 * t130 + t457 * t129 + t620 * t186) * g(2) + (-t669 + Ifges(5,3) * t672 - t684 / 0.2e1 + t683 / 0.2e1 + t119 / 0.2e1 - t685 / 0.2e1 + t207 * mrSges(4,1) - t172 / 0.2e1 - Ifges(4,4) * t576 + Ifges(5,5) * t580 - t157 * mrSges(4,3) + t670) * t269 + (-t541 + t682 / 0.2e1 + t262 / 0.2e1 + t207 * mrSges(4,2) + t173 / 0.2e1 + Ifges(4,1) * t576) * t268 + m(7) * (t1 * t30 + t10 * t32 + t2 * t31 + t24 * t7 + t25 * t6 + t48 * t5) + (-mrSges(4,3) * t88 + Ifges(4,1) * t200 + Ifges(4,4) * t201 + Ifges(4,5) * t305) * t282 + t294 * t404 + t293 * t405 + t360 * t322 * t512 + (Ifges(3,5) * t474 + Ifges(3,6) * t473 + Ifges(3,3) * t500) * t362 + m(3) * (t293 * t319 + t294 * t320) + (m(3) * (-pkin(1) * t347 + qJD(2) * t412) + t347 * (-mrSges(3,1) * t360 + mrSges(3,2) * t357) + (Ifges(3,4) * t357 + Ifges(3,2) * t360) * t473 + (Ifges(3,1) * t357 + Ifges(3,4) * t360) * t474 - pkin(1) * t406 + (Ifges(3,5) * t357 + Ifges(3,6) * t360) * t500) * t359 + (t101 * t65 - t174 * t29) * mrSges(6,2) + (Ifges(6,6) * t603 - t20 * mrSges(5,3) - Ifges(5,4) * t594 - Ifges(5,2) * t593 + Ifges(7,6) * t602 + t662 * t596 + t664 * t604 - Ifges(5,6) * t585 + t661 / 0.2e1 - t54 / 0.2e1 + t687) * t232 + t214 * t124 + t141 * t215 + t24 * (-mrSges(7,1) * t174 + mrSges(7,2) * t101) + t28 * (mrSges(6,1) * t174 - mrSges(6,3) * t101) + t136 * (mrSges(5,1) * t174 + mrSges(5,2) * t175) + t175 * t121 / 0.2e1 + t46 * t166 + t47 * t167 + t168 * t170 + t169 * t171 + m(4) * (t141 * t157 + t168 * t88 + t169 * t87 + t190 * t214 + t207 * t447) + t154 * t64 + t174 * t595 + t8 * t111 + t9 * t112 + t7 * t113 + t6 * t110 + t10 * t104 + t41 * t105 - (t494 / 0.2e1 - Ifges(4,4) * t200 - Ifges(4,6) * t305 - Ifges(4,2) * t201 - t87 * mrSges(4,3) + Ifges(5,6) * t593 + Ifges(5,5) * t594 + Ifges(5,3) * t585 + t613) * t377 + t30 * t33 + t31 * t35 + t37 * t34 + t48 * t22 + t73 * t23 + t89 * t91 + t90 * t92; (t22 + t655) * t323 - t658 * t394 + t659 * t291 + (t156 * t298 - t157 * t299 - t207 * t448 + t190 * t361 + (t572 * t88 + t365 * t87 + (-t156 * t365 + t157 * t572) * qJD(3)) * t358 + t616) * m(4) + (-t412 * t513 + t347 + t616) * m(3) - t199 * t448 + t695 * t215 + (-t64 + t170) * t492 - t322 * t482 + t321 * t483 + t406 + (t20 * t324 - t21 * t323 + (t136 * t511 - t572 * t84) * t358 - t136 * t298 + t637 * t71 - t636 * t70 + t616) * m(5) + t637 * t166 + (t13 * t323 + t291 * t3 + t394 * t4 + t636 * t65 + t616) * m(6) + (t1 * t291 - t2 * t394 + t32 * t636 + t323 * t5 + t616) * m(7) + t361 * t124 + t324 * t92 + t171 * t525 + (m(6) * t28 - m(7) * t24 + t680) * (qJD(5) * t453 - t324 * t507 - t363 * t637 + t366 * t689) + (m(6) * t29 + m(7) * t25 - t681) * (qJD(5) * t394 + t363 * t689 + t366 * t637) + t679 * t689 - t636 * (-t104 + t641); t4 * (-mrSges(6,1) * t367 - mrSges(6,3) * t517) + t2 * (mrSges(7,1) * t367 + mrSges(7,2) * t517) + (t363 * t467 + t364 * t469) * t77 + t656 * (t363 * t470 + t366 * t467 - t192 / 0.2e1) + t29 * (-mrSges(6,2) * t505 - mrSges(6,3) * t390) + t24 * (-mrSges(7,1) * t505 + mrSges(7,2) * t391) + t25 * (-mrSges(7,2) * t390 + mrSges(7,3) * t505) + t28 * (mrSges(6,1) * t505 - mrSges(6,3) * t391) - t614 + (-t530 / 0.2e1 + t467) * t121 + (t364 * t629 - t367 * t664) * t604 + (-t65 * mrSges(6,2) - t24 * mrSges(7,2) + t28 * mrSges(6,3) + t32 * mrSges(7,3) + Ifges(6,4) * t591 + Ifges(7,5) * t592 + t584 * t664 + t590 * t666) * t192 - t270 * t670 + (t498 - t95) * t105 + t267 * t541 + t647 * t112 + t648 * t111 + t649 * t104 + t650 * t113 + (t1 * t296 + t2 * t297 + t24 * t650 + t25 * t651 + t310 * t5 + t32 * t649) * m(7) + t651 * t110 + (-Ifges(6,2) * t591 + Ifges(7,3) * t592 - t584 * t690 + t590 * t665 + t609 + t621) * t191 + (t212 * t428 + t213 * t433 + t372 * t423) * qJD(4) / 0.2e1 - (Ifges(4,1) * t267 + t119 - t551) * t270 / 0.2e1 - (-Ifges(4,2) * t270 + t173 + t262) * t267 / 0.2e1 + t611 * t363 * t364 + t517 * t686 + t366 * t80 * t470 + (t505 / 0.2e1 - t518 / 0.2e1) * t657 + (t554 + t615) * t157 + (t120 / 0.2e1 + t608) * t518 + t493 + t84 * t438 + (-t84 * pkin(3) - t107 * t70 - t108 * t71) * m(5) + (t647 * t28 + t648 * t29 + t3 * t308 + t4 * t307 - t65 * t95) * m(6) + (((-t364 * t71 - t367 * t70) * qJD(4) + t625) * m(5) + (t503 * t65 + t542) * m(6) + t655 * t364 - t166 * t505 + t367 * t92) * pkin(10) - t367 * t675 + (-m(5) * t465 + mrSges(4,2) * t239 + t667 * t532 + t676 * (-pkin(11) * t532 - t238 * t568 + t465) + t645 * t238 - t460 * (-t238 * t515 + t239 * t363) + t457 * (-t238 * t519 - t239 * t366)) * g(1) + t65 * (mrSges(6,1) * t390 + mrSges(6,2) * t391) + t32 * (mrSges(7,1) * t390 - mrSges(7,3) * t391) + t372 * t136 * (mrSges(5,1) * t364 + mrSges(5,2) * t367) + t5 * t434 * t364 + (-t498 - t107) * t167 - t661 * t367 / 0.2e1 + (-t632 * t508 + (t364 * t662 + t367 * t631) * qJD(4)) * t583 + (t364 * t631 - t367 * t662) * t596 + (-t630 * t508 + (t364 * t664 + t367 * t629) * qJD(4)) * t589 + (Ifges(5,3) * t270 + t267 * t423) * t673 + t367 * t674 + t367 * t54 / 0.2e1 + t310 * t22 - t306 * (Ifges(4,5) * t267 - Ifges(4,6) * t270) / 0.2e1 + t307 * t34 + t308 * t36 + t296 * t33 + t297 * t35 - t207 * (mrSges(4,1) * t270 + mrSges(4,2) * t267) - t156 * t215 - t108 * t166 + t479 * t599 + (-Ifges(7,6) * t367 + t364 * t420) * t602 + (-Ifges(6,6) * t367 + t364 * t427) * t603 + t364 * t605 + (-t426 * t508 + (Ifges(6,6) * t364 + t367 * t427) * qJD(4)) * t592 + (Ifges(5,2) * t367 + t549) * t593 + (Ifges(5,1) * t364 + t548) * t594 + t505 * t595 + (-t419 * t508 + (Ifges(7,6) * t364 + t367 * t420) * qJD(4)) * t591 + (Ifges(5,6) * t270 + t267 * t428) * t582 + (Ifges(5,5) * t364 + Ifges(5,6) * t367) * t585 + t172 * t576 + (Ifges(5,5) * t270 + t267 * t433) * t581 + (-m(5) * t466 - mrSges(4,2) * t237 + t667 * t533 + t676 * (-pkin(11) * t533 - t234 * t568 + t466) + t645 * t234 - t460 * (-t234 * t515 - t237 * t363) + t457 * (-t234 * t519 + t237 * t366)) * g(2) + (-g(1) * t239 + g(2) * t237 - g(3) * t282 - t688 * t71 + (-t503 + t530) * t70 + t625) * mrSges(5,3) + (-m(5) * t464 + t377 * t438 + t440 - t667 * t529 + t676 * (pkin(11) * t529 + t367 * t560 + t464) - t460 * (t282 * t363 + t377 * t515) + t457 * (-t282 * t366 + t377 * t519)) * g(3) + t436 * t542 - pkin(3) * t64 + t270 * t669; (-t497 - t52) * t111 + (t24 * t627 + t25 * t626 + t610 + t623) * mrSges(7,2) + (-t28 * t627 + t29 * t626 + t610 + t624) * mrSges(6,3) + t613 + (-t550 + t657) * t581 + (-t427 / 0.2e1 + t420 / 0.2e1) * qJD(5) * t164 + t630 * t604 + t632 * t596 + t633 * qJD(5) + (-mrSges(5,2) * t185 + t676 * (pkin(4) * t692 - pkin(11) * t185) - t619 * t692) * g(2) + t646 * t104 + (t334 * t5 + ((t24 * t366 - t25 * t363) * qJD(5) + t623) * pkin(11) - t24 * t43 - t25 * t42 + t646 * t32) * m(7) + (t165 * t629 + t211 * t631) * qJD(5) / 0.2e1 + t678 * t366 + (-t496 - t51) * t112 + (t496 - t43) * t113 + t363 * t686 - pkin(4) * t23 + (t210 + t121) * t582 + t494 + t5 * t435 - t13 * t437 + (-m(6) * t65 + t536 + t641) * t71 + (t635 * t232 + t461 * t676 + t439) * g(3) + (mrSges(5,2) * t187 + t676 * (-t186 * pkin(4) + pkin(11) * t187) + t619 * t186) * g(1) + (-t136 * mrSges(5,1) - Ifges(5,2) * t582 + Ifges(5,6) * t672 + t608) * t213 + (-t136 * mrSges(5,2) + Ifges(5,1) * t581 + Ifges(5,5) * t673 + t420 * t592 + t427 * t591 + t584 * t631 + t590 * t629 - t633) * t212 + (-t534 / 0.2e1 + t469) * t656 + (-t497 - t42) * t110 + (t537 - t166) * t70 + t621 * t535 + t622 * t509 + (-pkin(4) * t13 + ((-t28 * t366 - t29 * t363) * qJD(5) + t624) * pkin(11) - t28 * t51 - t29 * t52) * m(6) + t658 * pkin(11) * t363 + t659 * pkin(11) * t366 + t334 * t22 + t419 * t602 + t426 * t603 + t120 * t580; (-Ifges(6,2) * t165 - t163 + t656) * t591 + (-t164 * t664 - t165 * t690) * t584 + (-t164 * t666 + t162 - t540 + t77) * t590 + (t552 + t680) * t29 + (-t553 + t681) * t28 + (-t457 * t698 - t460 * t699) * g(2) + (t129 * t460 + t130 * t457) * g(1) + (t176 * t460 + t177 * t457) * g(3) + t661 + (-pkin(5) * t2 + qJ(6) * t1 - t103 * t32 - t24 * t29 + t25 * t638) * m(7) + t687 + (t164 * t24 + t165 * t25) * mrSges(7,2) - t32 * (mrSges(7,1) * t165 + mrSges(7,3) * t164) - t65 * (mrSges(6,1) * t165 - mrSges(6,2) * t164) + (Ifges(7,3) * t165 - t545) * t592 + t80 * t589 + qJD(6) * t110 - t103 * t104 + qJ(6) * t33 - pkin(5) * t35; t165 * t104 - t211 * t110 + (-g(1) * t129 + g(2) * t699 - g(3) * t176 + t165 * t32 - t211 * t25 + t2) * m(7) + t35;];
tau  = t11;
