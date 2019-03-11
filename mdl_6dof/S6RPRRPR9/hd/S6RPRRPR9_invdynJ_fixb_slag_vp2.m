% Calculate vector of inverse dynamics joint torques for
% S6RPRRPR9
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPR9_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR9_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR9_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR9_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_invdynJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR9_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR9_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR9_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:27:02
% EndTime: 2019-03-09 05:28:15
% DurationCPUTime: 46.47s
% Computational Cost: add. (38662->950), mult. (123333->1301), div. (0->0), fcn. (107303->18), ass. (0->427)
t377 = sin(pkin(12));
t379 = sin(pkin(6));
t381 = cos(pkin(12));
t390 = cos(qJ(3));
t382 = cos(pkin(7));
t387 = sin(qJ(3));
t499 = t382 * t387;
t402 = (-t377 * t499 + t381 * t390) * t379;
t316 = qJD(1) * t402;
t378 = sin(pkin(7));
t486 = qJD(3) * t390;
t460 = t378 * t486;
t675 = t316 - t460;
t502 = t379 * t381;
t362 = qJ(2) * t502;
t383 = cos(pkin(6));
t553 = pkin(1) * t383;
t478 = qJD(1) * t553;
t328 = qJD(1) * t362 + t377 * t478;
t501 = t379 * t382;
t505 = t378 * t383;
t404 = (t381 * t501 + t505) * pkin(9);
t270 = qJD(1) * t404 + t328;
t548 = pkin(9) * t377;
t318 = (-pkin(2) * t381 - t378 * t548 - pkin(1)) * t379;
t309 = qJD(1) * t318 + qJD(2);
t360 = t381 * t478;
t507 = t377 * t379;
t552 = pkin(2) * t383;
t400 = t552 + (-pkin(9) * t382 - qJ(2)) * t507;
t277 = qJD(1) * t400 + t360;
t520 = t277 * t382;
t174 = -t387 * t270 + t390 * (t309 * t378 + t520);
t498 = t382 * t390;
t417 = -t377 * t387 + t381 * t498;
t503 = t378 * t390;
t397 = t379 * t417 + t383 * t503;
t281 = t397 * qJD(1);
t416 = t377 * t390 + t381 * t499;
t504 = t378 * t387;
t300 = t379 * t416 + t383 * t504;
t284 = t300 * qJD(1);
t215 = pkin(3) * t284 - pkin(10) * t281;
t386 = sin(qJ(4));
t389 = cos(qJ(4));
t121 = -t174 * t386 + t389 * t215;
t384 = -qJ(5) - pkin(10);
t456 = qJD(4) * t384;
t518 = t281 * t389;
t674 = -pkin(4) * t284 + qJ(5) * t518 - qJD(5) * t386 + t389 * t456 - t121;
t122 = t389 * t174 + t386 * t215;
t519 = t281 * t386;
t673 = -qJ(5) * t519 - qJD(5) * t389 - t386 * t456 + t122;
t342 = t382 * t389 - t386 * t504;
t488 = qJD(1) * t379;
t463 = t377 * t488;
t447 = t378 * t463;
t644 = qJD(4) * t342 - t386 * t447 - t389 * t675;
t343 = t382 * t386 + t389 * t504;
t643 = -qJD(4) * t343 + t386 * t675 - t389 * t447;
t555 = sin(qJ(1));
t464 = t555 * t381;
t391 = cos(qJ(1));
t496 = t391 * t377;
t334 = t383 * t496 + t464;
t367 = t555 * t377;
t495 = t391 * t381;
t441 = t383 * t495 - t367;
t500 = t379 * t391;
t667 = -t378 * t500 + t382 * t441;
t243 = t334 * t390 + t387 * t667;
t302 = t441 * t378 + t382 * t500;
t515 = t302 * t386;
t672 = -t243 * t389 + t515;
t276 = qJD(4) - t281;
t566 = -t276 / 0.2e1;
t333 = -t378 * t502 + t382 * t383;
t323 = qJD(1) * t333 + qJD(3);
t224 = -t284 * t386 + t323 * t389;
t225 = t284 * t389 + t323 * t386;
t376 = sin(pkin(13));
t380 = cos(pkin(13));
t426 = t224 * t376 + t380 * t225;
t572 = -t426 / 0.2e1;
t453 = t380 * t224 - t225 * t376;
t574 = -t453 / 0.2e1;
t671 = -Ifges(6,4) * t572 - Ifges(6,2) * t574 - Ifges(6,6) * t566;
t151 = qJD(6) - t453;
t575 = t151 / 0.2e1;
t385 = sin(qJ(6));
t388 = cos(qJ(6));
t128 = t276 * t385 + t388 * t426;
t581 = t128 / 0.2e1;
t127 = t276 * t388 - t385 * t426;
t583 = t127 / 0.2e1;
t670 = Ifges(7,5) * t581 + Ifges(7,6) * t583 + Ifges(7,3) * t575;
t655 = Ifges(6,3) + Ifges(5,3);
t624 = t376 * t674 - t673 * t380;
t646 = t376 * t643 + t380 * t644;
t403 = (t377 * t498 + t381 * t387) * t379;
t315 = qJD(1) * t403;
t487 = qJD(3) * t387;
t461 = t378 * t487;
t668 = -t315 + t461;
t175 = t390 * t270 + t277 * t499 + t309 * t504;
t485 = qJD(4) * t386;
t617 = -t175 + (t485 - t519) * pkin(4);
t149 = -pkin(3) * t323 - t174;
t120 = -pkin(4) * t224 + qJD(5) + t149;
t216 = -t277 * t378 + t382 * t309;
t148 = -pkin(3) * t281 - pkin(10) * t284 + t216;
t150 = pkin(10) * t323 + t175;
t100 = t389 * t148 - t150 * t386;
t78 = -qJ(5) * t225 + t100;
t75 = pkin(4) * t276 + t78;
t101 = t148 * t386 + t150 * t389;
t79 = qJ(5) * t224 + t101;
t76 = t380 * t79;
t39 = t376 * t75 + t76;
t35 = pkin(11) * t276 + t39;
t65 = -pkin(5) * t453 - pkin(11) * t426 + t120;
t16 = -t35 * t385 + t388 * t65;
t17 = t35 * t388 + t385 * t65;
t666 = mrSges(6,1) * t120 + mrSges(7,1) * t16 - mrSges(7,2) * t17 - mrSges(6,3) * t39 + t670 - t671;
t375 = qJ(4) + pkin(13);
t372 = sin(t375);
t373 = cos(t375);
t190 = t243 * t373 - t302 * t372;
t665 = -t243 * t372 - t302 * t373;
t664 = t243 * t386 + t302 * t389;
t436 = -t388 * mrSges(7,1) + t385 * mrSges(7,2);
t630 = m(7) * pkin(5) + mrSges(6,1) - t436;
t451 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t663 = -pkin(11) * t284 + t624;
t346 = t376 * t389 + t380 * t386;
t202 = t346 * t281;
t345 = t376 * t386 - t380 * t389;
t203 = t345 * t281;
t336 = t346 * qJD(4);
t338 = t345 * qJD(4);
t662 = t617 + (-t203 + t338) * pkin(11) + (-t202 + t336) * pkin(5);
t344 = (qJ(2) * qJDD(1) + qJD(1) * qJD(2)) * t379;
t474 = qJDD(1) * t553;
t310 = -t344 * t377 + t381 * t474;
t255 = (-t501 * t548 + t552) * qJDD(1) + t310;
t305 = qJDD(1) * t318 + qJDD(2);
t311 = t381 * t344 + t377 * t474;
t645 = qJD(3) * t520 + qJDD(1) * t404 + t311;
t108 = t255 * t499 - t270 * t487 + t305 * t504 + t309 * t460 + t390 * t645;
t109 = t390 * (t255 * t382 + t305 * t378) - t270 * t486 - t309 * t461 - t645 * t387;
t481 = qJD(1) * qJD(3);
t213 = (qJDD(1) * t387 + t390 * t481) * t505 + (qJDD(1) * t416 + t417 * t481) * t379;
t214 = t379 * (qJDD(1) * t417 - t416 * t481) - (-qJDD(1) * t390 + t387 * t481) * t505;
t322 = qJDD(1) * t333 + qJDD(3);
t661 = t109 * mrSges(4,1) - t108 * mrSges(4,2) + Ifges(4,5) * t213 + Ifges(4,6) * t214 + Ifges(4,3) * t322;
t211 = qJDD(4) - t214;
t570 = t211 / 0.2e1;
t136 = qJD(4) * t224 + t213 * t389 + t322 * t386;
t137 = -qJD(4) * t225 - t213 * t386 + t322 * t389;
t89 = t136 * t380 + t137 * t376;
t589 = t89 / 0.2e1;
t88 = -t136 * t376 + t137 * t380;
t590 = t88 / 0.2e1;
t86 = qJDD(6) - t88;
t591 = t86 / 0.2e1;
t55 = -qJD(6) * t128 + t211 * t388 - t385 * t89;
t599 = t55 / 0.2e1;
t54 = qJD(6) * t127 + t211 * t385 + t388 * t89;
t600 = t54 / 0.2e1;
t104 = -pkin(3) * t322 - t109;
t66 = -pkin(4) * t137 + qJDD(5) + t104;
t19 = -pkin(5) * t88 - pkin(11) * t89 + t66;
t103 = pkin(10) * t322 + t108;
t204 = -t255 * t378 + t382 * t305;
t118 = -pkin(3) * t214 - pkin(10) * t213 + t204;
t33 = -qJD(4) * t101 - t103 * t386 + t389 * t118;
t21 = pkin(4) * t211 - qJ(5) * t136 - qJD(5) * t225 + t33;
t484 = qJD(4) * t389;
t32 = t389 * t103 + t386 * t118 + t148 * t484 - t150 * t485;
t27 = qJ(5) * t137 + qJD(5) * t224 + t32;
t8 = t376 * t21 + t380 * t27;
t6 = pkin(11) * t211 + t8;
t1 = qJD(6) * t16 + t19 * t385 + t388 * t6;
t2 = -qJD(6) * t17 + t19 * t388 - t385 * t6;
t610 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t9 = Ifges(7,5) * t54 + Ifges(7,6) * t55 + Ifges(7,3) * t86;
t659 = t610 + mrSges(6,1) * t66 - mrSges(6,3) * t8 + Ifges(7,5) * t600 + Ifges(7,6) * t599 + Ifges(7,3) * t591 + t9 / 0.2e1 + (-t570 - t211 / 0.2e1) * Ifges(6,6) + (-t590 - t88 / 0.2e1) * Ifges(6,2) + (-t589 - t89 / 0.2e1) * Ifges(6,4);
t565 = t276 / 0.2e1;
t571 = t426 / 0.2e1;
t573 = t453 / 0.2e1;
t657 = -Ifges(6,4) * t571 - Ifges(6,2) * t573 - Ifges(6,6) * t565 + t666 + t670;
t7 = t21 * t380 - t27 * t376;
t656 = t33 * mrSges(5,1) + t7 * mrSges(6,1) - t32 * mrSges(5,2) - t8 * mrSges(6,2);
t240 = -t300 * t386 + t333 * t389;
t443 = pkin(4) * t240;
t526 = t376 * t79;
t38 = t380 * t75 - t526;
t545 = t38 * mrSges(6,3);
t275 = Ifges(4,4) * t281;
t654 = t216 * mrSges(4,2);
t653 = t224 * Ifges(5,6);
t652 = t281 * Ifges(4,2);
t651 = t225 * Ifges(5,5) + Ifges(6,5) * t426 + Ifges(6,6) * t453 + t276 * t655 + t653;
t626 = t673 * t376 + t380 * t674;
t130 = mrSges(6,1) * t276 - mrSges(6,3) * t426;
t80 = -mrSges(7,1) * t127 + mrSges(7,2) * t128;
t650 = t80 - t130;
t268 = t342 * t376 + t343 * t380;
t418 = -t268 * t388 + t385 * t503;
t649 = qJD(6) * t418 - t385 * t646 + t388 * t668;
t248 = -t268 * t385 - t388 * t503;
t648 = qJD(6) * t248 + t385 * t668 + t388 * t646;
t647 = t376 * t644 - t380 * t643;
t339 = t377 * t553 + t362;
t291 = t404 + t339;
t366 = t381 * t553;
t301 = t366 + t400;
t423 = t301 * t382 + t318 * t378;
t180 = -t387 * t291 + t423 * t390;
t405 = t383 * t464 + t496;
t465 = t379 * t555;
t304 = t405 * t378 + t382 * t465;
t642 = Ifges(5,5) * t136 + Ifges(6,5) * t89 + Ifges(5,6) * t137 + Ifges(6,6) * t88 + t211 * t655;
t335 = -t367 * t383 + t495;
t618 = -t378 * t465 + t405 * t382;
t247 = t335 * t390 - t387 * t618;
t197 = -t247 * t386 + t304 * t389;
t576 = -t151 / 0.2e1;
t582 = -t128 / 0.2e1;
t584 = -t127 / 0.2e1;
t640 = Ifges(7,5) * t582 + Ifges(7,6) * t584 + Ifges(7,3) * t576 - t666 + t671;
t244 = -t334 * t387 + t390 * t667;
t99 = Ifges(6,1) * t426 + Ifges(6,4) * t453 + t276 * Ifges(6,5);
t637 = Ifges(6,1) * t571 + Ifges(6,4) * t573 + Ifges(6,5) * t565 + t99 / 0.2e1;
t636 = -mrSges(6,3) * t7 + 0.2e1 * Ifges(6,1) * t589 + 0.2e1 * Ifges(6,4) * t590 + 0.2e1 * Ifges(6,5) * t570;
t558 = -t323 / 0.2e1;
t635 = -t216 * mrSges(4,1) - t100 * mrSges(5,1) - t38 * mrSges(6,1) + t101 * mrSges(5,2) - t558 * Ifges(4,6);
t632 = -m(7) - m(6);
t578 = t136 / 0.2e1;
t577 = t137 / 0.2e1;
t18 = -mrSges(7,1) * t55 + mrSges(7,2) * t54;
t68 = mrSges(6,1) * t211 - mrSges(6,3) * t89;
t631 = t18 - t68;
t371 = pkin(4) * t389 + pkin(3);
t285 = pkin(5) * t345 - pkin(11) * t346 - t371;
t353 = t384 * t389;
t458 = t384 * t386;
t313 = -t380 * t353 + t376 * t458;
t217 = t285 * t388 - t313 * t385;
t628 = qJD(6) * t217 + t385 * t662 + t388 * t663;
t218 = t285 * t385 + t313 * t388;
t627 = -qJD(6) * t218 - t385 * t663 + t388 * t662;
t625 = pkin(5) * t284 - t626;
t166 = t203 * t385 + t284 * t388;
t482 = qJD(6) * t388;
t414 = -t385 * t338 + t346 * t482;
t620 = t166 + t414;
t167 = -t203 * t388 + t284 * t385;
t483 = qJD(6) * t385;
t413 = t388 * t338 + t346 * t483;
t619 = t167 + t413;
t226 = -t301 * t378 + t382 * t318;
t165 = -pkin(3) * t397 - pkin(10) * t300 + t226;
t274 = t390 * t291;
t181 = t301 * t499 + t318 * t504 + t274;
t173 = pkin(10) * t333 + t181;
t111 = t386 * t165 + t389 * t173;
t539 = mrSges(4,3) * t284;
t491 = -mrSges(4,1) * t323 - mrSges(5,1) * t224 + mrSges(5,2) * t225 + t539;
t241 = t300 * t389 + t333 * t386;
t177 = t240 * t376 + t241 * t380;
t290 = t397 * t388;
t143 = -t177 * t385 - t290;
t435 = mrSges(7,1) * t385 + mrSges(7,2) * t388;
t440 = -m(5) * pkin(10) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t612 = -t435 + t440;
t616 = t32 * t389 - t33 * t386;
t615 = t1 * t388 - t2 * t385;
t614 = t120 * mrSges(6,2) - t545;
t613 = t377 * ((-mrSges(4,1) * t281 + mrSges(4,2) * t284) * t378 - (mrSges(3,1) * t383 - mrSges(3,3) * t507) * qJD(1)) + (-mrSges(3,2) * t383 + mrSges(3,3) * t502) * qJD(1) * t381;
t438 = -mrSges(5,1) * t389 + mrSges(5,2) * t386;
t611 = m(5) * pkin(3) - t451 * t372 + t373 * t630 + mrSges(4,1) - t438;
t607 = Ifges(6,1) * t572 + Ifges(6,4) * t574 + Ifges(6,5) * t566 - t99 / 0.2e1;
t605 = 0.2e1 * t383;
t603 = Ifges(7,1) * t600 + Ifges(7,4) * t599 + Ifges(7,5) * t591;
t536 = Ifges(7,4) * t128;
t61 = Ifges(7,2) * t127 + Ifges(7,6) * t151 + t536;
t595 = t61 / 0.2e1;
t125 = Ifges(7,4) * t127;
t62 = Ifges(7,1) * t128 + Ifges(7,5) * t151 + t125;
t594 = -t62 / 0.2e1;
t593 = Ifges(5,4) * t578 + Ifges(5,2) * t577 + Ifges(5,6) * t570;
t592 = Ifges(5,1) * t578 + Ifges(5,4) * t577 + Ifges(5,5) * t570;
t531 = t225 * Ifges(5,4);
t134 = t224 * Ifges(5,2) + t276 * Ifges(5,6) + t531;
t580 = t134 / 0.2e1;
t223 = Ifges(5,4) * t224;
t135 = t225 * Ifges(5,1) + t276 * Ifges(5,5) + t223;
t579 = t135 / 0.2e1;
t569 = -t224 / 0.2e1;
t568 = -t225 / 0.2e1;
t567 = t225 / 0.2e1;
t564 = -t281 / 0.2e1;
t560 = t284 / 0.2e1;
t556 = t388 / 0.2e1;
t551 = pkin(4) * t225;
t550 = pkin(4) * t376;
t549 = pkin(4) * t380;
t282 = t397 * qJD(3);
t188 = qJD(4) * t240 + t282 * t389;
t283 = t300 * qJD(3);
t159 = qJD(2) * t402 + qJD(3) * t180;
t446 = qJD(2) * t378 * t507;
t201 = pkin(3) * t283 - pkin(10) * t282 + t446;
t64 = -qJD(4) * t111 - t159 * t386 + t389 * t201;
t46 = pkin(4) * t283 - qJ(5) * t188 - qJD(5) * t241 + t64;
t187 = -qJD(4) * t241 - t282 * t386;
t63 = t389 * t159 + t165 * t484 - t173 * t485 + t386 * t201;
t50 = qJ(5) * t187 + qJD(5) * t240 + t63;
t15 = t376 * t46 + t380 * t50;
t110 = t389 * t165 - t173 * t386;
t87 = -pkin(4) * t397 - qJ(5) * t241 + t110;
t96 = qJ(5) * t240 + t111;
t48 = t376 * t87 + t380 * t96;
t542 = mrSges(3,1) * t381;
t541 = mrSges(3,2) * t377;
t540 = mrSges(4,3) * t281;
t538 = Ifges(5,4) * t386;
t537 = Ifges(5,4) * t389;
t535 = Ifges(7,4) * t385;
t534 = Ifges(7,4) * t388;
t533 = t100 * mrSges(5,3);
t532 = t101 * mrSges(5,3);
t530 = t284 * Ifges(4,4);
t528 = t323 * Ifges(4,5);
t525 = t453 * t385;
t524 = t453 * t388;
t517 = t397 * t385;
t513 = t304 * t386;
t509 = t346 * t385;
t508 = t346 * t388;
t489 = t391 * pkin(1) + qJ(2) * t465;
t480 = qJDD(1) * t379;
t469 = t62 * t556;
t107 = -mrSges(6,1) * t453 + mrSges(6,2) * t426;
t467 = -t107 - t491;
t44 = -t88 * mrSges(6,1) + t89 * mrSges(6,2);
t457 = -t483 / 0.2e1;
t445 = t664 * pkin(4);
t444 = t197 * pkin(4);
t442 = -pkin(1) * t555 + qJ(2) * t500;
t439 = mrSges(5,1) * t240 - mrSges(5,2) * t241;
t434 = Ifges(5,1) * t389 - t538;
t433 = Ifges(7,1) * t388 - t535;
t432 = -Ifges(5,2) * t386 + t537;
t431 = -Ifges(7,2) * t385 + t534;
t430 = Ifges(5,5) * t389 - Ifges(5,6) * t386;
t429 = Ifges(7,5) * t388 - Ifges(7,6) * t385;
t14 = -t376 * t50 + t380 * t46;
t47 = -t376 * t96 + t380 * t87;
t43 = -pkin(11) * t397 + t48;
t172 = -pkin(3) * t333 - t180;
t126 = t172 - t443;
t176 = -t380 * t240 + t241 * t376;
t71 = pkin(5) * t176 - pkin(11) * t177 + t126;
t23 = t385 * t71 + t388 * t43;
t22 = -t385 * t43 + t388 * t71;
t93 = -mrSges(7,2) * t151 + mrSges(7,3) * t127;
t94 = mrSges(7,1) * t151 - mrSges(7,3) * t128;
t427 = -t385 * t94 + t388 * t93;
t144 = t177 * t388 - t517;
t422 = -(-qJ(2) * t463 + t360) * t377 + t328 * t381;
t34 = -pkin(5) * t276 - t38;
t415 = t34 * t435;
t412 = t149 * (mrSges(5,1) * t386 + mrSges(5,2) * t389);
t246 = t335 * t387 + t390 * t618;
t407 = -g(1) * t246 + g(2) * t244 + g(3) * t397;
t396 = -t334 * pkin(2) + pkin(9) * t302 + t442;
t394 = t335 * pkin(2) + pkin(9) * t304 + t489;
t392 = pkin(4) * t513 - t246 * t384 + t247 * t371 + t394;
t160 = qJD(2) * t403 + (t387 * t423 + t274) * qJD(3);
t119 = -pkin(4) * t187 + t160;
t370 = -pkin(5) - t549;
t361 = -pkin(1) * t480 + qJDD(2);
t354 = t480 * t541;
t337 = -qJ(2) * t507 + t366;
t312 = -t353 * t376 - t380 * t458;
t267 = -t380 * t342 + t343 * t376;
t234 = t300 * t373 + t333 * t372;
t231 = -mrSges(4,2) * t323 + t540;
t198 = t247 * t389 + t513;
t194 = t247 * t373 + t304 * t372;
t193 = t247 * t372 - t304 * t373;
t186 = t284 * Ifges(4,1) + t275 + t528;
t185 = t323 * Ifges(4,6) + t530 + t652;
t183 = -mrSges(4,2) * t322 + mrSges(4,3) * t214;
t182 = mrSges(4,1) * t322 - mrSges(4,3) * t213;
t179 = mrSges(5,1) * t276 - mrSges(5,3) * t225;
t178 = -mrSges(5,2) * t276 + mrSges(5,3) * t224;
t145 = -mrSges(4,1) * t214 + mrSges(4,2) * t213;
t142 = t194 * t388 + t246 * t385;
t141 = -t194 * t385 + t246 * t388;
t129 = -mrSges(6,2) * t276 + mrSges(6,3) * t453;
t124 = t187 * t376 + t188 * t380;
t123 = -t380 * t187 + t188 * t376;
t113 = -mrSges(5,2) * t211 + mrSges(5,3) * t137;
t112 = mrSges(5,1) * t211 - mrSges(5,3) * t136;
t91 = pkin(5) * t426 - pkin(11) * t453 + t551;
t90 = -mrSges(5,1) * t137 + mrSges(5,2) * t136;
t74 = -qJD(6) * t144 - t124 * t385 + t283 * t388;
t73 = qJD(6) * t143 + t124 * t388 + t283 * t385;
t67 = -mrSges(6,2) * t211 + mrSges(6,3) * t88;
t51 = pkin(5) * t123 - pkin(11) * t124 + t119;
t42 = pkin(5) * t397 - t47;
t41 = t380 * t78 - t526;
t40 = t376 * t78 + t76;
t29 = -mrSges(7,2) * t86 + mrSges(7,3) * t55;
t28 = mrSges(7,1) * t86 - mrSges(7,3) * t54;
t25 = t385 * t91 + t388 * t41;
t24 = -t385 * t41 + t388 * t91;
t13 = pkin(11) * t283 + t15;
t12 = -pkin(5) * t283 - t14;
t10 = t54 * Ifges(7,4) + t55 * Ifges(7,2) + t86 * Ifges(7,6);
t5 = -pkin(5) * t211 - t7;
t4 = -qJD(6) * t23 - t13 * t385 + t388 * t51;
t3 = qJD(6) * t22 + t13 * t388 + t385 * t51;
t11 = [t224 * (Ifges(5,4) * t188 + Ifges(5,2) * t187) / 0.2e1 + (Ifges(5,4) * t241 + Ifges(5,2) * t240) * t577 + (Ifges(7,5) * t73 + Ifges(7,6) * t74) * t575 + (Ifges(7,5) * t144 + Ifges(7,6) * t143) * t591 - t104 * t439 + t636 * t177 + (mrSges(4,2) * t204 - mrSges(4,3) * t109 + Ifges(4,1) * t213 + Ifges(4,4) * t214 + Ifges(4,5) * t322) * t300 + (t120 * t124 + t177 * t66) * mrSges(6,2) + (Ifges(2,3) + (mrSges(3,1) * t337 - mrSges(3,2) * t339 + Ifges(3,3) * t383) * t383 + ((-mrSges(3,3) * t337 + Ifges(3,1) * t507 + Ifges(3,5) * t605) * t377 + (t339 * mrSges(3,3) + Ifges(3,6) * t605 + (mrSges(3,1) * pkin(1) + 0.2e1 * Ifges(3,4) * t377 + Ifges(3,2) * t381) * t379) * t381) * t379) * qJDD(1) + (-m(4) * t396 + t243 * mrSges(4,1) - t302 * mrSges(4,3) - m(3) * t442 + t334 * mrSges(3,1) + t441 * mrSges(3,2) - mrSges(3,3) * t500 - m(5) * (-pkin(3) * t243 + t396) - t672 * mrSges(5,1) - t664 * mrSges(5,2) + t555 * mrSges(2,1) + t391 * mrSges(2,2) + t451 * t665 + t630 * t190 + t612 * t244 + t632 * (pkin(4) * t515 - t243 * t371 - t244 * t384 + t396)) * g(1) + (-pkin(1) * t354 + t361 * (t541 - t542) + (-t310 * t377 + t311 * t381) * mrSges(3,3) + t613 * qJD(2)) * t379 + t491 * t160 + m(3) * (t310 * t337 + t311 * t339 + (-pkin(1) * t361 + qJD(2) * t422) * t379) + (Ifges(7,1) * t73 + Ifges(7,4) * t74) * t581 + (Ifges(7,1) * t144 + Ifges(7,4) * t143) * t600 + m(6) * (t119 * t120 + t126 * t66 + t14 * t38 + t15 * t39 + t47 * t7 + t48 * t8) + m(5) * (t100 * t64 + t101 * t63 + t104 * t172 + t110 * t33 + t111 * t32 + t149 * t160) + m(7) * (t1 * t23 + t12 * t34 + t16 * t4 + t17 * t3 + t2 * t22 + t42 * t5) + (t637 - t545) * t124 + (Ifges(5,5) * t188 + Ifges(5,6) * t187) * t565 + (Ifges(5,5) * t241 + Ifges(5,6) * t240) * t570 + (-m(5) * (t247 * pkin(3) + t394) - t198 * mrSges(5,1) - t197 * mrSges(5,2) - m(3) * t489 - t335 * mrSges(3,1) + t405 * mrSges(3,2) - mrSges(3,3) * t465 - t391 * mrSges(2,1) + t555 * mrSges(2,2) - m(4) * t394 - t247 * mrSges(4,1) - t304 * mrSges(4,3) - m(6) * t392 - t194 * mrSges(6,1) - m(7) * (t194 * pkin(5) + t392) - t142 * mrSges(7,1) - t141 * mrSges(7,2) + t451 * t193 + t440 * t246) * g(2) + (t528 + t186) * t282 / 0.2e1 + (mrSges(3,1) * t310 - mrSges(3,2) * t311) * t383 + (Ifges(5,1) * t188 + Ifges(5,4) * t187) * t567 + (Ifges(5,1) * t241 + Ifges(5,4) * t240) * t578 - (mrSges(4,1) * t204 - mrSges(4,3) * t108 - Ifges(4,4) * t213 + Ifges(5,5) * t578 + Ifges(6,5) * t589 - Ifges(4,2) * t214 - Ifges(4,6) * t322 + Ifges(5,6) * t577 + Ifges(6,6) * t590 + t570 * t655 + t642 / 0.2e1 + t656) * t397 + (-t174 * mrSges(4,3) + t654 + t275 / 0.2e1 + Ifges(4,1) * t560) * t282 + (t1 * t143 - t144 * t2 - t16 * t73 + t17 * t74) * mrSges(7,3) + (Ifges(7,4) * t73 + Ifges(7,2) * t74) * t583 + (Ifges(7,4) * t144 + Ifges(7,2) * t143) * t599 + (-t100 * t188 + t101 * t187 + t240 * t32 - t241 * t33) * mrSges(5,3) + m(4) * (t108 * t181 + t109 * t180 + t159 * t175 - t160 * t174 + t204 * t226 + t216 * t446) + t226 * t145 + t159 * t231 + t149 * (-mrSges(5,1) * t187 + mrSges(5,2) * t188) + t64 * t179 + t180 * t182 + t181 * t183 + t63 * t178 + t172 * t90 + t5 * (-mrSges(7,1) * t143 + mrSges(7,2) * t144) + t143 * t10 / 0.2e1 + t15 * t129 + t14 * t130 + t126 * t44 + t119 * t107 + t110 * t112 + t111 * t113 + t3 * t93 + t4 * t94 + t12 * t80 + t34 * (-mrSges(7,1) * t74 + mrSges(7,2) * t73) + t73 * t62 / 0.2e1 + t48 * t67 + t47 * t68 + t42 * t18 + t657 * t123 + t659 * t176 + t661 * t333 + t188 * t579 + t187 * t580 + t241 * t592 + t240 * t593 + t74 * t595 + t144 * t603 + (-t39 * mrSges(6,2) + Ifges(6,6) * t573 - t175 * mrSges(4,3) + Ifges(6,5) * t571 + t653 / 0.2e1 - t652 / 0.2e1 - t185 / 0.2e1 - Ifges(4,4) * t560 + Ifges(5,5) * t567 + t655 * t565 - t635 + t651 / 0.2e1) * t283 + t22 * t28 + t23 * t29; (-qJD(1) * t613 - qJDD(1) * t542) * t379 + t467 * t315 + (t387 * t183 + (t182 - t44 - t90) * t390 + (t231 * t390 - t387 * t467) * qJD(3)) * t378 + t631 * t267 + t644 * t178 + t648 * t93 + t643 * t179 - t316 * t231 + t268 * t67 + t248 * t28 - t418 * t29 + t382 * t145 + t342 * t112 + t343 * t113 + t649 * t94 + t354 + t646 * t129 + t650 * t647 + ((-g(1) * t555 + g(2) * t391) * t379 - t383 * g(3)) * (m(3) + m(4) + m(5) - t632) + (-t1 * t418 + t16 * t649 + t17 * t648 + t2 * t248 + t267 * t5 + t34 * t647) * m(7) + (-t267 * t7 + t268 * t8 + (t120 * t487 - t390 * t66) * t378 - t120 * t315 + t646 * t39 - t647 * t38) * m(6) + (t32 * t343 + t33 * t342 + (-t104 * t390 + t149 * t487) * t378 - t149 * t315 + t644 * t101 + t643 * t100) * m(5) + (t174 * t315 - t175 * t316 - t216 * t447 + t204 * t382 + (t108 * t387 + t109 * t390 + (-t174 * t387 + t175 * t390) * qJD(3)) * t378) * m(4) + (-t422 * t488 + t361) * m(3); t104 * t438 + (-m(5) * t149 - t491 + t539) * t175 + (mrSges(7,1) * t620 - mrSges(7,2) * t619) * t34 - t620 * t61 / 0.2e1 + (-t1 * t509 + t16 * t619 - t17 * t620 - t2 * t508) * mrSges(7,3) + (m(5) * ((-t100 * t389 - t101 * t386) * qJD(4) + t616) + t389 * t113 - t178 * t485 - t386 * t112 - t179 * t484) * pkin(10) + (t100 * t518 + t101 * t519 + t616) * mrSges(5,3) + t617 * t107 + (mrSges(6,2) * t66 + t429 * t591 + t431 * t599 + t433 * t600 + t435 * t5 + t457 * t62 + t636) * t346 + (-pkin(3) * t104 - t100 * t121 - t101 * t122) * m(5) + (-Ifges(7,5) * t413 - Ifges(7,6) * t414) * t575 + t624 * t129 + t625 * t80 + (t120 * t617 - t312 * t7 + t313 * t8 - t371 * t66 + t38 * t626 + t39 * t624) * m(6) + t626 * t130 + t627 * t94 + t628 * t93 + (t1 * t218 + t16 * t627 + t17 * t628 + t2 * t217 + t312 * t5 + t34 * t625) * m(7) - t10 * t509 / 0.2e1 - (t469 + t614 + t637) * t338 + t631 * t312 - (t607 + t545) * t203 + (t224 * t432 + t225 * t434 + t276 * t430) * qJD(4) / 0.2e1 + (t275 + t186) * t564 + (t540 - t231) * t174 + (Ifges(4,5) * t558 + t430 * t566 + t432 * t569 + t434 * t568 - t412 - t654) * t281 + qJD(4) * t412 + (Ifges(5,5) * t568 + Ifges(6,5) * t572 - Ifges(4,2) * t564 + Ifges(5,6) * t569 + Ifges(6,6) * t574 + t566 * t655 + t635) * t284 + (Ifges(7,5) * t167 + Ifges(7,6) * t166) * t576 - t135 * t518 / 0.2e1 + t640 * t202 + (t632 * (-t243 * t384 + t244 * t371) + t612 * t243 - t611 * t244) * g(2) + (t632 * (-t300 * t384 + t371 * t397) + t612 * t300 - t611 * t397) * g(3) + (t120 * t203 + t284 * t39) * mrSges(6,2) - (Ifges(4,1) * t281 - t530 + t651) * t284 / 0.2e1 + (-Ifges(7,1) * t413 - Ifges(7,4) * t414) * t581 + (Ifges(7,1) * t167 + Ifges(7,4) * t166) * t582 + (-Ifges(7,4) * t413 - Ifges(7,2) * t414) * t583 + (Ifges(7,4) * t167 + Ifges(7,2) * t166) * t584 + t661 + (-t533 + t579) * t484 + (-t134 / 0.2e1 - t532) * t485 + t313 * t67 + t217 * t28 + t218 * t29 - t121 * t179 - t122 * t178 - pkin(3) * t90 + (t632 * (-t246 * t371 - t247 * t384) + t612 * t247 + t611 * t246) * g(1) + t657 * t336 + t659 * t345 + t185 * t560 + (Ifges(5,5) * t386 + Ifges(5,6) * t389) * t570 + (Ifges(5,2) * t389 + t538) * t577 + (Ifges(5,1) * t386 + t537) * t578 + t519 * t580 + t386 * t592 + t389 * t593 + t167 * t594 + t508 * t603 - t371 * t44; (-Ifges(5,2) * t225 + t135 + t223) * t569 + t5 * t436 + (-g(1) * t194 - g(2) * t190 - g(3) * t234 + (-t483 + t525) * t17 + (-t482 + t524) * t16 + t615) * mrSges(7,3) + (m(7) * ((-t16 * t388 - t17 * t385) * qJD(6) + t615) + t388 * t29 - t94 * t482 - t93 * t483 - t385 * t28) * (pkin(11) + t550) + (t429 * t576 + t431 * t584 + t433 * t582 - t415 + t607 - t614) * t453 + (-m(6) * t444 + t194 * mrSges(6,2) - m(7) * (pkin(11) * t194 + t444) - mrSges(5,1) * t197 + mrSges(5,2) * t198 + t630 * t193) * g(1) + (-t439 - m(6) * t443 + t234 * mrSges(6,2) - m(7) * (pkin(11) * t234 + t443) - t630 * (-t300 * t372 + t333 * t373)) * g(3) + (t415 + t469) * qJD(6) - t107 * t551 + (t127 * t431 + t128 * t433 + t151 * t429) * qJD(6) / 0.2e1 + t642 + (-t16 * t24 - t17 * t25 - t34 * t40 + t370 * t5) * m(7) + (t664 * mrSges(5,1) - t672 * mrSges(5,2) + m(6) * t445 + t190 * mrSges(6,2) - m(7) * (pkin(11) * t190 - t445) - t630 * t665) * g(2) + t640 * t426 + (-t120 * t551 + t38 * t40 - t39 * t41 + (t376 * t8 + t380 * t7) * pkin(4)) * m(6) + t225 * t532 + t224 * t533 - t650 * t40 + t656 - t149 * (mrSges(5,1) * t225 + mrSges(5,2) * t224) + t101 * t179 - t100 * t178 - t41 * t129 - t25 * t93 - t24 * t94 + t68 * t549 + t67 * t550 + t61 * t457 + t10 * t556 + (Ifges(5,5) * t224 - Ifges(5,6) * t225) * t566 + t134 * t567 + (Ifges(5,1) * t224 - t531) * t568 + (Ifges(7,5) * t385 + Ifges(7,6) * t388) * t591 + t524 * t594 + t525 * t595 + (Ifges(7,2) * t388 + t535) * t599 + (Ifges(7,1) * t385 + t534) * t600 + t385 * t603 + t370 * t18; t388 * t28 + t385 * t29 - t650 * t426 + t427 * qJD(6) + (-t129 - t427) * t453 + t44 + (t1 * t385 - t426 * t34 + t2 * t388 + t407 + t151 * (-t16 * t385 + t17 * t388)) * m(7) + (t38 * t426 - t39 * t453 + t407 + t66) * m(6); -t16 * t93 + t17 * t94 - t34 * (mrSges(7,1) * t128 + mrSges(7,2) * t127) + (Ifges(7,1) * t127 - t536) * t582 + t61 * t581 + (Ifges(7,5) * t127 - Ifges(7,6) * t128) * t576 - g(1) * (mrSges(7,1) * t141 - mrSges(7,2) * t142) - g(2) * ((-t190 * t385 - t244 * t388) * mrSges(7,1) + (-t190 * t388 + t244 * t385) * mrSges(7,2)) - g(3) * ((-t234 * t385 - t290) * mrSges(7,1) + (-t234 * t388 + t517) * mrSges(7,2)) + (t127 * t16 + t128 * t17) * mrSges(7,3) + t9 + (-Ifges(7,2) * t128 + t125 + t62) * t584 + t610;];
tau  = t11;
