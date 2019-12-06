% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRPR5
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRPR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR5_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR5_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR5_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR5_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:25:59
% EndTime: 2019-12-05 16:26:12
% DurationCPUTime: 6.44s
% Computational Cost: add. (69986->268), mult. (147916->353), div. (0->0), fcn. (101914->12), ass. (0->115)
t512 = sin(pkin(10));
t515 = cos(pkin(10));
t519 = sin(qJ(3));
t522 = cos(qJ(3));
t490 = (t512 * t519 - t515 * t522) * qJD(2);
t513 = sin(pkin(9));
t516 = cos(pkin(9));
t504 = t513 * g(1) - t516 * g(2);
t505 = -t516 * g(1) - t513 * g(2);
t511 = -g(3) + qJDD(1);
t514 = sin(pkin(5));
t517 = cos(pkin(5));
t520 = sin(qJ(2));
t523 = cos(qJ(2));
t472 = -t520 * t505 + (t504 * t517 + t511 * t514) * t523;
t546 = 2 * qJD(4);
t525 = qJD(2) ^ 2;
t528 = -qJDD(2) * pkin(2) - t472;
t467 = -t525 * pkin(7) + t528;
t539 = qJD(2) * qJD(3);
t538 = t522 * t539;
t502 = t519 * qJDD(2) + t538;
t503 = t522 * qJDD(2) - t519 * t539;
t541 = qJD(2) * t519;
t507 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t541;
t540 = qJD(2) * t522;
t508 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t540;
t543 = t517 * t520;
t544 = t514 * t520;
t473 = t504 * t543 + t523 * t505 + t511 * t544;
t468 = -t525 * pkin(2) + qJDD(2) * pkin(7) + t473;
t487 = -t514 * t504 + t517 * t511;
t454 = -t519 * t468 + t522 * t487;
t451 = (-t502 + t538) * qJ(4) + (t519 * t522 * t525 + qJDD(3)) * pkin(3) + t454;
t455 = t522 * t468 + t519 * t487;
t506 = qJD(3) * pkin(3) - qJ(4) * t541;
t510 = t522 ^ 2;
t452 = -t510 * t525 * pkin(3) + t503 * qJ(4) - qJD(3) * t506 + t455;
t447 = t512 * t451 + t515 * t452 - t490 * t546;
t491 = (t512 * t522 + t515 * t519) * qJD(2);
t476 = t490 * pkin(4) - t491 * pkin(8);
t524 = qJD(3) ^ 2;
t445 = -t524 * pkin(4) + qJDD(3) * pkin(8) - t490 * t476 + t447;
t453 = -t503 * pkin(3) + qJDD(4) + t506 * t541 + (-qJ(4) * t510 - pkin(7)) * t525 + t528;
t479 = -t512 * t502 + t515 * t503;
t480 = t515 * t502 + t512 * t503;
t448 = (qJD(3) * t490 - t480) * pkin(8) + (qJD(3) * t491 - t479) * pkin(4) + t453;
t518 = sin(qJ(5));
t521 = cos(qJ(5));
t442 = -t518 * t445 + t521 * t448;
t481 = t521 * qJD(3) - t518 * t491;
t462 = t481 * qJD(5) + t518 * qJDD(3) + t521 * t480;
t482 = t518 * qJD(3) + t521 * t491;
t463 = -t481 * mrSges(6,1) + t482 * mrSges(6,2);
t489 = qJD(5) + t490;
t465 = -t489 * mrSges(6,2) + t481 * mrSges(6,3);
t478 = qJDD(5) - t479;
t440 = m(6) * t442 + t478 * mrSges(6,1) - t462 * mrSges(6,3) - t482 * t463 + t489 * t465;
t443 = t521 * t445 + t518 * t448;
t461 = -t482 * qJD(5) + t521 * qJDD(3) - t518 * t480;
t466 = t489 * mrSges(6,1) - t482 * mrSges(6,3);
t441 = m(6) * t443 - t478 * mrSges(6,2) + t461 * mrSges(6,3) + t481 * t463 - t489 * t466;
t432 = t521 * t440 + t518 * t441;
t485 = -qJD(3) * mrSges(5,2) - t490 * mrSges(5,3);
t486 = qJD(3) * mrSges(5,1) - t491 * mrSges(5,3);
t527 = m(5) * t453 - t479 * mrSges(5,1) + t480 * mrSges(5,2) + t490 * t485 + t491 * t486 + t432;
t526 = -m(4) * t467 + t503 * mrSges(4,1) - t502 * mrSges(4,2) - t507 * t541 + t508 * t540 - t527;
t428 = m(3) * t472 + qJDD(2) * mrSges(3,1) - t525 * mrSges(3,2) + t526;
t545 = t428 * t523;
t475 = t490 * mrSges(5,1) + t491 * mrSges(5,2);
t534 = -t518 * t440 + t521 * t441;
t431 = m(5) * t447 - qJDD(3) * mrSges(5,2) + t479 * mrSges(5,3) - qJD(3) * t486 - t490 * t475 + t534;
t533 = -t515 * t451 + t512 * t452;
t446 = -0.2e1 * qJD(4) * t491 - t533;
t444 = -qJDD(3) * pkin(4) - t524 * pkin(8) + (t546 + t476) * t491 + t533;
t529 = -m(6) * t444 + t461 * mrSges(6,1) - t462 * mrSges(6,2) + t481 * t465 - t482 * t466;
t436 = m(5) * t446 + qJDD(3) * mrSges(5,1) - t480 * mrSges(5,3) + qJD(3) * t485 - t491 * t475 + t529;
t425 = t512 * t431 + t515 * t436;
t501 = (-mrSges(4,1) * t522 + mrSges(4,2) * t519) * qJD(2);
t423 = m(4) * t454 + qJDD(3) * mrSges(4,1) - t502 * mrSges(4,3) + qJD(3) * t508 - t501 * t541 + t425;
t535 = t515 * t431 - t512 * t436;
t424 = m(4) * t455 - qJDD(3) * mrSges(4,2) + t503 * mrSges(4,3) - qJD(3) * t507 + t501 * t540 + t535;
t536 = -t519 * t423 + t522 * t424;
t414 = m(3) * t473 - t525 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t536;
t417 = t522 * t423 + t519 * t424;
t416 = m(3) * t487 + t417;
t404 = t414 * t543 - t514 * t416 + t517 * t545;
t402 = m(2) * t504 + t404;
t410 = t523 * t414 - t520 * t428;
t409 = m(2) * t505 + t410;
t542 = t516 * t402 + t513 * t409;
t403 = t414 * t544 + t517 * t416 + t514 * t545;
t537 = -t513 * t402 + t516 * t409;
t456 = Ifges(6,5) * t482 + Ifges(6,6) * t481 + Ifges(6,3) * t489;
t458 = Ifges(6,1) * t482 + Ifges(6,4) * t481 + Ifges(6,5) * t489;
t433 = -mrSges(6,1) * t444 + mrSges(6,3) * t443 + Ifges(6,4) * t462 + Ifges(6,2) * t461 + Ifges(6,6) * t478 - t482 * t456 + t489 * t458;
t457 = Ifges(6,4) * t482 + Ifges(6,2) * t481 + Ifges(6,6) * t489;
t434 = mrSges(6,2) * t444 - mrSges(6,3) * t442 + Ifges(6,1) * t462 + Ifges(6,4) * t461 + Ifges(6,5) * t478 + t481 * t456 - t489 * t457;
t469 = Ifges(5,5) * t491 - Ifges(5,6) * t490 + Ifges(5,3) * qJD(3);
t470 = Ifges(5,4) * t491 - Ifges(5,2) * t490 + Ifges(5,6) * qJD(3);
t418 = mrSges(5,2) * t453 - mrSges(5,3) * t446 + Ifges(5,1) * t480 + Ifges(5,4) * t479 + Ifges(5,5) * qJDD(3) - pkin(8) * t432 - qJD(3) * t470 - t518 * t433 + t521 * t434 - t490 * t469;
t471 = Ifges(5,1) * t491 - Ifges(5,4) * t490 + Ifges(5,5) * qJD(3);
t419 = -mrSges(5,1) * t453 - mrSges(6,1) * t442 + mrSges(6,2) * t443 + mrSges(5,3) * t447 + Ifges(5,4) * t480 - Ifges(6,5) * t462 + Ifges(5,2) * t479 + Ifges(5,6) * qJDD(3) - Ifges(6,6) * t461 - Ifges(6,3) * t478 - pkin(4) * t432 + qJD(3) * t471 - t482 * t457 + t481 * t458 - t491 * t469;
t493 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t519 + Ifges(4,6) * t522) * qJD(2);
t495 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t519 + Ifges(4,4) * t522) * qJD(2);
t405 = -mrSges(4,1) * t467 + mrSges(4,3) * t455 + Ifges(4,4) * t502 + Ifges(4,2) * t503 + Ifges(4,6) * qJDD(3) - pkin(3) * t527 + qJ(4) * t535 + qJD(3) * t495 + t512 * t418 + t515 * t419 - t493 * t541;
t494 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t519 + Ifges(4,2) * t522) * qJD(2);
t406 = mrSges(4,2) * t467 - mrSges(4,3) * t454 + Ifges(4,1) * t502 + Ifges(4,4) * t503 + Ifges(4,5) * qJDD(3) - qJ(4) * t425 - qJD(3) * t494 + t515 * t418 - t512 * t419 + t493 * t540;
t399 = mrSges(3,2) * t487 - mrSges(3,3) * t472 + Ifges(3,5) * qJDD(2) - t525 * Ifges(3,6) - pkin(7) * t417 - t519 * t405 + t522 * t406;
t400 = -pkin(2) * t417 + mrSges(3,3) * t473 - mrSges(3,1) * t487 + Ifges(3,6) * qJDD(2) - pkin(3) * t425 - Ifges(4,5) * t502 - Ifges(4,6) * t503 - mrSges(4,1) * t454 + mrSges(4,2) * t455 - t521 * t433 - pkin(4) * t529 - pkin(8) * t534 - Ifges(5,5) * t480 - Ifges(5,6) * t479 - mrSges(5,1) * t446 + mrSges(5,2) * t447 - t518 * t434 - t491 * t470 - t490 * t471 + t525 * Ifges(3,5) + (-Ifges(4,3) - Ifges(5,3)) * qJDD(3) + (-t519 * t494 + t522 * t495) * qJD(2);
t530 = pkin(6) * t410 + t399 * t520 + t400 * t523;
t398 = mrSges(3,1) * t472 - mrSges(3,2) * t473 + Ifges(3,3) * qJDD(2) + pkin(2) * t526 + pkin(7) * t536 + t522 * t405 + t519 * t406;
t397 = mrSges(2,2) * t511 - mrSges(2,3) * t504 + t523 * t399 - t520 * t400 + (-t403 * t514 - t404 * t517) * pkin(6);
t396 = -mrSges(2,1) * t511 + mrSges(2,3) * t505 - pkin(1) * t403 - t514 * t398 + t530 * t517;
t1 = [-m(1) * g(1) + t537; -m(1) * g(2) + t542; -m(1) * g(3) + m(2) * t511 + t403; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t542 - t513 * t396 + t516 * t397; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t537 + t516 * t396 + t513 * t397; -mrSges(1,1) * g(2) + mrSges(2,1) * t504 + mrSges(1,2) * g(1) - mrSges(2,2) * t505 + pkin(1) * t404 + t517 * t398 + t530 * t514;];
tauB = t1;
