% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PPRPRR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRPRR1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRPRR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_invdynB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRPRR1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRPRR1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRPRR1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:52:55
% EndTime: 2019-05-04 19:53:15
% DurationCPUTime: 14.91s
% Computational Cost: add. (268730->252), mult. (463290->332), div. (0->0), fcn. (380510->16), ass. (0->122)
t516 = sin(pkin(11));
t521 = cos(pkin(11));
t509 = -g(1) * t521 - g(2) * t516;
t515 = sin(pkin(12));
t520 = cos(pkin(12));
t508 = g(1) * t516 - g(2) * t521;
t513 = -g(3) + qJDD(1);
t518 = sin(pkin(6));
t523 = cos(pkin(6));
t536 = t508 * t523 + t513 * t518;
t480 = -t515 * t509 + t520 * t536;
t481 = t520 * t509 + t515 * t536;
t495 = -t508 * t518 + t513 * t523 + qJDD(2);
t526 = sin(qJ(3));
t522 = cos(pkin(7));
t529 = cos(qJ(3));
t549 = t522 * t529;
t517 = sin(pkin(7));
t551 = t517 * t529;
t471 = t480 * t549 - t481 * t526 + t495 * t551;
t469 = qJDD(3) * pkin(3) + t471;
t550 = t522 * t526;
t552 = t517 * t526;
t472 = t480 * t550 + t481 * t529 + t495 * t552;
t531 = qJD(3) ^ 2;
t470 = -pkin(3) * t531 + t472;
t514 = sin(pkin(13));
t519 = cos(pkin(13));
t465 = t469 * t514 + t470 * t519;
t463 = -pkin(4) * t531 + qJDD(3) * pkin(9) + t465;
t476 = -t480 * t517 + t495 * t522;
t475 = qJDD(4) + t476;
t525 = sin(qJ(5));
t528 = cos(qJ(5));
t459 = t463 * t528 + t475 * t525;
t504 = (-mrSges(6,1) * t528 + mrSges(6,2) * t525) * qJD(3);
t544 = qJD(3) * qJD(5);
t542 = t525 * t544;
t507 = qJDD(3) * t528 - t542;
t546 = qJD(3) * t525;
t510 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t546;
t505 = (-pkin(5) * t528 - pkin(10) * t525) * qJD(3);
t530 = qJD(5) ^ 2;
t545 = t528 * qJD(3);
t457 = -pkin(5) * t530 + qJDD(5) * pkin(10) + t505 * t545 + t459;
t464 = t519 * t469 - t470 * t514;
t462 = -qJDD(3) * pkin(4) - t531 * pkin(9) - t464;
t541 = t528 * t544;
t506 = qJDD(3) * t525 + t541;
t460 = (-t506 - t541) * pkin(10) + (-t507 + t542) * pkin(5) + t462;
t524 = sin(qJ(6));
t527 = cos(qJ(6));
t454 = -t457 * t524 + t460 * t527;
t502 = qJD(5) * t527 - t524 * t546;
t488 = qJD(6) * t502 + qJDD(5) * t524 + t506 * t527;
t503 = qJD(5) * t524 + t527 * t546;
t489 = -mrSges(7,1) * t502 + mrSges(7,2) * t503;
t512 = qJD(6) - t545;
t493 = -mrSges(7,2) * t512 + mrSges(7,3) * t502;
t501 = qJDD(6) - t507;
t452 = m(7) * t454 + mrSges(7,1) * t501 - mrSges(7,3) * t488 - t489 * t503 + t493 * t512;
t455 = t457 * t527 + t460 * t524;
t487 = -qJD(6) * t503 + qJDD(5) * t527 - t506 * t524;
t494 = mrSges(7,1) * t512 - mrSges(7,3) * t503;
t453 = m(7) * t455 - mrSges(7,2) * t501 + mrSges(7,3) * t487 + t489 * t502 - t494 * t512;
t537 = -t452 * t524 + t453 * t527;
t445 = m(6) * t459 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t507 - qJD(5) * t510 + t504 * t545 + t537;
t548 = t528 * t475;
t458 = -t463 * t525 + t548;
t511 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t545;
t456 = -qJDD(5) * pkin(5) - t530 * pkin(10) - t548 + (qJD(3) * t505 + t463) * t525;
t533 = -m(7) * t456 + mrSges(7,1) * t487 - mrSges(7,2) * t488 + t493 * t502 - t494 * t503;
t450 = m(6) * t458 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t506 + qJD(5) * t511 - t504 * t546 + t533;
t538 = t445 * t528 - t450 * t525;
t437 = m(5) * t465 - mrSges(5,1) * t531 - qJDD(3) * mrSges(5,2) + t538;
t446 = t452 * t527 + t453 * t524;
t532 = -m(6) * t462 + mrSges(6,1) * t507 - mrSges(6,2) * t506 - t510 * t546 + t511 * t545 - t446;
t442 = m(5) * t464 + qJDD(3) * mrSges(5,1) - mrSges(5,2) * t531 + t532;
t432 = t437 * t514 + t442 * t519;
t430 = m(4) * t471 + qJDD(3) * mrSges(4,1) - mrSges(4,2) * t531 + t432;
t539 = t437 * t519 - t442 * t514;
t431 = m(4) * t472 - mrSges(4,1) * t531 - qJDD(3) * mrSges(4,2) + t539;
t440 = t445 * t525 + t450 * t528;
t543 = m(5) * t475 + t440;
t439 = m(4) * t476 + t543;
t417 = t430 * t549 + t431 * t550 - t439 * t517;
t413 = m(3) * t480 + t417;
t423 = -t430 * t526 + t431 * t529;
t422 = m(3) * t481 + t423;
t555 = t413 * t520 + t422 * t515;
t416 = t430 * t551 + t431 * t552 + t439 * t522;
t415 = m(3) * t495 + t416;
t403 = -t518 * t415 + t523 * t555;
t401 = m(2) * t508 + t403;
t409 = -t413 * t515 + t422 * t520;
t408 = m(2) * t509 + t409;
t547 = t401 * t521 + t408 * t516;
t402 = t523 * t415 + t518 * t555;
t540 = -t401 * t516 + t408 * t521;
t482 = Ifges(7,5) * t503 + Ifges(7,6) * t502 + Ifges(7,3) * t512;
t484 = Ifges(7,1) * t503 + Ifges(7,4) * t502 + Ifges(7,5) * t512;
t447 = -mrSges(7,1) * t456 + mrSges(7,3) * t455 + Ifges(7,4) * t488 + Ifges(7,2) * t487 + Ifges(7,6) * t501 - t482 * t503 + t484 * t512;
t483 = Ifges(7,4) * t503 + Ifges(7,2) * t502 + Ifges(7,6) * t512;
t448 = mrSges(7,2) * t456 - mrSges(7,3) * t454 + Ifges(7,1) * t488 + Ifges(7,4) * t487 + Ifges(7,5) * t501 + t482 * t502 - t483 * t512;
t496 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t525 + Ifges(6,6) * t528) * qJD(3);
t497 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t525 + Ifges(6,2) * t528) * qJD(3);
t433 = mrSges(6,2) * t462 - mrSges(6,3) * t458 + Ifges(6,1) * t506 + Ifges(6,4) * t507 + Ifges(6,5) * qJDD(5) - pkin(10) * t446 - qJD(5) * t497 - t447 * t524 + t448 * t527 + t496 * t545;
t498 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t525 + Ifges(6,4) * t528) * qJD(3);
t434 = -mrSges(6,1) * t462 - mrSges(7,1) * t454 + mrSges(7,2) * t455 + mrSges(6,3) * t459 + Ifges(6,4) * t506 - Ifges(7,5) * t488 + Ifges(6,2) * t507 + Ifges(6,6) * qJDD(5) - Ifges(7,6) * t487 - Ifges(7,3) * t501 - pkin(5) * t446 + qJD(5) * t498 - t483 * t503 + t484 * t502 - t496 * t546;
t418 = mrSges(5,2) * t475 - mrSges(5,3) * t464 + Ifges(5,5) * qJDD(3) - Ifges(5,6) * t531 - pkin(9) * t440 + t433 * t528 - t434 * t525;
t424 = Ifges(5,6) * qJDD(3) + t531 * Ifges(5,5) - mrSges(5,1) * t475 + mrSges(5,3) * t465 - Ifges(6,5) * t506 - Ifges(6,6) * t507 - Ifges(6,3) * qJDD(5) - mrSges(6,1) * t458 + mrSges(6,2) * t459 - t524 * t448 - t527 * t447 - pkin(5) * t533 - pkin(10) * t537 - pkin(4) * t440 + (-t497 * t525 + t498 * t528) * qJD(3);
t404 = -mrSges(4,1) * t476 + mrSges(4,3) * t472 + t531 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t543 + qJ(4) * t539 + t514 * t418 + t519 * t424;
t405 = mrSges(4,2) * t476 - mrSges(4,3) * t471 + Ifges(4,5) * qJDD(3) - Ifges(4,6) * t531 - qJ(4) * t432 + t418 * t519 - t424 * t514;
t535 = pkin(8) * t423 + t404 * t529 + t405 * t526;
t410 = mrSges(4,1) * t471 - mrSges(4,2) * t472 + mrSges(5,1) * t464 - mrSges(5,2) * t465 + t525 * t433 + t528 * t434 + pkin(4) * t532 + pkin(9) * t538 + pkin(3) * t432 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3);
t398 = -mrSges(3,1) * t495 + mrSges(3,3) * t481 - pkin(2) * t416 - t517 * t410 + t522 * t535;
t399 = mrSges(3,2) * t495 - mrSges(3,3) * t480 - t526 * t404 + t529 * t405 + (-t416 * t517 - t417 * t522) * pkin(8);
t534 = qJ(2) * t409 + t398 * t520 + t399 * t515;
t397 = mrSges(3,1) * t480 - mrSges(3,2) * t481 + pkin(2) * t417 + t522 * t410 + t517 * t535;
t396 = mrSges(2,2) * t513 - mrSges(2,3) * t508 - t515 * t398 + t520 * t399 + (-t402 * t518 - t403 * t523) * qJ(2);
t395 = -mrSges(2,1) * t513 + mrSges(2,3) * t509 - pkin(1) * t402 - t518 * t397 + t523 * t534;
t1 = [-m(1) * g(1) + t540; -m(1) * g(2) + t547; -m(1) * g(3) + m(2) * t513 + t402; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t547 - t395 * t516 + t396 * t521; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t540 + t521 * t395 + t516 * t396; -mrSges(1,1) * g(2) + mrSges(2,1) * t508 + mrSges(1,2) * g(1) - mrSges(2,2) * t509 + pkin(1) * t403 + t523 * t397 + t518 * t534;];
tauB  = t1;
