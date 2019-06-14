% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 22:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRPR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:39:02
% EndTime: 2019-05-05 22:39:07
% DurationCPUTime: 4.56s
% Computational Cost: add. (38732->299), mult. (94630->362), div. (0->0), fcn. (73331->10), ass. (0->131)
t568 = Ifges(5,4) + Ifges(6,6);
t578 = -Ifges(5,2) - Ifges(6,3);
t575 = Ifges(5,6) - Ifges(6,5);
t577 = Ifges(5,1) + Ifges(6,2);
t576 = Ifges(5,5) - Ifges(6,4);
t574 = Ifges(5,3) + Ifges(6,1);
t527 = sin(pkin(10));
t528 = cos(pkin(10));
t531 = sin(qJ(3));
t534 = cos(qJ(3));
t548 = -t527 * t531 + t528 * t534;
t511 = t548 * qJD(1);
t549 = t527 * t534 + t528 * t531;
t512 = t549 * qJD(1);
t530 = sin(qJ(4));
t571 = cos(qJ(4));
t494 = -t511 * t571 + t512 * t530;
t495 = t530 * t511 + t512 * t571;
t526 = qJD(3) + qJD(4);
t573 = t578 * t494 + t568 * t495 + t575 * t526;
t536 = qJD(1) ^ 2;
t572 = -2 * qJD(5);
t570 = pkin(2) * t536;
t569 = mrSges(5,1) - mrSges(6,2);
t567 = pkin(7) * qJDD(1);
t566 = t494 * t526;
t532 = sin(qJ(1));
t535 = cos(qJ(1));
t552 = -g(1) * t535 - g(2) * t532;
t513 = -pkin(1) * t536 + qJDD(1) * qJ(2) + t552;
t558 = qJD(1) * qJD(2);
t556 = -g(3) * t528 - 0.2e1 * t527 * t558;
t485 = (t528 * t570 - t513 - t567) * t527 + t556;
t503 = -g(3) * t527 + (t513 + 0.2e1 * t558) * t528;
t525 = t528 ^ 2;
t487 = -t525 * t570 + t528 * t567 + t503;
t459 = t534 * t485 - t487 * t531;
t559 = qJD(3) * t511;
t501 = qJDD(1) * t549 + t559;
t428 = (-t501 + t559) * pkin(8) + (t511 * t512 + qJDD(3)) * pkin(3) + t459;
t460 = t531 * t485 + t534 * t487;
t500 = -qJD(3) * t512 + qJDD(1) * t548;
t506 = qJD(3) * pkin(3) - pkin(8) * t512;
t510 = t511 ^ 2;
t436 = -pkin(3) * t510 + pkin(8) * t500 - qJD(3) * t506 + t460;
t424 = t428 * t571 - t530 * t436;
t457 = -t494 * qJD(4) + t530 * t500 + t501 * t571;
t472 = mrSges(5,1) * t494 + mrSges(5,2) * t495;
t523 = qJDD(3) + qJDD(4);
t471 = pkin(4) * t494 - qJ(5) * t495;
t522 = t526 ^ 2;
t420 = -t523 * pkin(4) - t522 * qJ(5) + t495 * t471 + qJDD(5) - t424;
t414 = (t494 * t495 - t523) * pkin(9) + (t457 + t566) * pkin(5) + t420;
t456 = qJD(4) * t495 - t500 * t571 + t501 * t530;
t483 = pkin(5) * t495 - pkin(9) * t526;
t490 = t494 ^ 2;
t557 = t532 * g(1) - t535 * g(2);
t551 = qJDD(2) - t557;
t561 = -t527 ^ 2 - t525;
t499 = (-pkin(2) * t528 - pkin(1)) * qJDD(1) + (pkin(7) * t561 - qJ(2)) * t536 + t551;
t444 = -t500 * pkin(3) - t510 * pkin(8) + t512 * t506 + t499;
t537 = (-t457 + t566) * qJ(5) + t444 + (t526 * pkin(4) + t572) * t495;
t415 = t537 + (pkin(4) + pkin(9)) * t456 - t490 * pkin(5) - t495 * t483;
t529 = sin(qJ(6));
t533 = cos(qJ(6));
t412 = t414 * t533 - t415 * t529;
t475 = t494 * t533 - t526 * t529;
t433 = qJD(6) * t475 + t456 * t529 + t523 * t533;
t455 = qJDD(6) + t457;
t476 = t494 * t529 + t526 * t533;
t458 = -mrSges(7,1) * t475 + mrSges(7,2) * t476;
t489 = qJD(6) + t495;
t461 = -mrSges(7,2) * t489 + mrSges(7,3) * t475;
t409 = m(7) * t412 + mrSges(7,1) * t455 - mrSges(7,3) * t433 - t458 * t476 + t461 * t489;
t413 = t414 * t529 + t415 * t533;
t432 = -qJD(6) * t476 + t456 * t533 - t523 * t529;
t462 = mrSges(7,1) * t489 - mrSges(7,3) * t476;
t410 = m(7) * t413 - mrSges(7,2) * t455 + mrSges(7,3) * t432 + t458 * t475 - t462 * t489;
t400 = t409 * t533 + t410 * t529;
t473 = -mrSges(6,2) * t494 - mrSges(6,3) * t495;
t544 = -m(6) * t420 - t457 * mrSges(6,1) - t495 * t473 - t400;
t481 = mrSges(6,1) * t494 - mrSges(6,3) * t526;
t562 = -mrSges(5,2) * t526 - mrSges(5,3) * t494 - t481;
t396 = m(5) * t424 - mrSges(5,3) * t457 - t472 * t495 + t523 * t569 + t526 * t562 + t544;
t425 = t530 * t428 + t571 * t436;
t480 = mrSges(5,1) * t526 - mrSges(5,3) * t495;
t543 = -pkin(4) * t522 + qJ(5) * t523 - t471 * t494 + t425;
t418 = t526 * t572 - t543;
t482 = mrSges(6,1) * t495 + mrSges(6,2) * t526;
t417 = -pkin(5) * t456 - pkin(9) * t490 + ((2 * qJD(5)) + t483) * t526 + t543;
t545 = -m(7) * t417 + mrSges(7,1) * t432 - t433 * mrSges(7,2) + t461 * t475 - t476 * t462;
t541 = -m(6) * t418 + t523 * mrSges(6,3) + t526 * t482 - t545;
t406 = m(5) * t425 - mrSges(5,2) * t523 - t480 * t526 + (-t472 - t473) * t494 + (-mrSges(5,3) - mrSges(6,1)) * t456 + t541;
t393 = t571 * t396 + t530 * t406;
t498 = -mrSges(4,1) * t511 + mrSges(4,2) * t512;
t504 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t511;
t391 = m(4) * t459 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t501 + qJD(3) * t504 - t498 * t512 + t393;
t505 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t512;
t554 = -t396 * t530 + t571 * t406;
t392 = m(4) * t460 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t500 - qJD(3) * t505 + t498 * t511 + t554;
t565 = t534 * t391 + t531 * t392;
t564 = t575 * t494 - t576 * t495 - t574 * t526;
t563 = -t568 * t494 + t577 * t495 + t576 * t526;
t555 = -t531 * t391 + t534 * t392;
t553 = -t529 * t409 + t533 * t410;
t550 = -mrSges(3,1) * t528 + mrSges(3,2) * t527;
t547 = mrSges(3,3) * qJDD(1) + t536 * t550;
t422 = t456 * pkin(4) + t537;
t546 = m(6) * t422 - t457 * mrSges(6,3) - t495 * t482 + t553;
t438 = Ifges(7,4) * t476 + Ifges(7,2) * t475 + Ifges(7,6) * t489;
t439 = Ifges(7,1) * t476 + Ifges(7,4) * t475 + Ifges(7,5) * t489;
t542 = mrSges(7,1) * t412 - mrSges(7,2) * t413 + Ifges(7,5) * t433 + Ifges(7,6) * t432 + Ifges(7,3) * t455 + t476 * t438 - t475 * t439;
t540 = m(5) * t444 + t457 * mrSges(5,2) + t456 * t569 + t495 * t480 + t494 * t562 + t546;
t399 = mrSges(6,2) * t523 + t481 * t526 - t544;
t437 = Ifges(7,5) * t476 + Ifges(7,6) * t475 + Ifges(7,3) * t489;
t402 = -mrSges(7,1) * t417 + mrSges(7,3) * t413 + Ifges(7,4) * t433 + Ifges(7,2) * t432 + Ifges(7,6) * t455 - t437 * t476 + t439 * t489;
t403 = mrSges(7,2) * t417 - mrSges(7,3) * t412 + Ifges(7,1) * t433 + Ifges(7,4) * t432 + Ifges(7,5) * t455 + t437 * t475 - t438 * t489;
t539 = -mrSges(5,2) * t425 - mrSges(6,3) * t418 - pkin(4) * t399 - pkin(9) * t400 - t529 * t402 + t533 * t403 + t494 * t563 + qJ(5) * (-t473 * t494 + t541) + mrSges(6,2) * t420 + mrSges(5,1) * t424 + t574 * t523 + t573 * t495 + t576 * t457 + (-qJ(5) * mrSges(6,1) - t575) * t456;
t538 = m(4) * t499 - t500 * mrSges(4,1) + t501 * mrSges(4,2) - t511 * t504 + t512 * t505 + t540;
t509 = -qJDD(1) * pkin(1) - qJ(2) * t536 + t551;
t502 = -t513 * t527 + t556;
t493 = Ifges(4,1) * t512 + Ifges(4,4) * t511 + Ifges(4,5) * qJD(3);
t492 = Ifges(4,4) * t512 + Ifges(4,2) * t511 + Ifges(4,6) * qJD(3);
t491 = Ifges(4,5) * t512 + Ifges(4,6) * t511 + Ifges(4,3) * qJD(3);
t397 = -t456 * mrSges(6,2) - t494 * t481 + t546;
t394 = mrSges(3,3) * t536 * t561 + m(3) * t509 + qJDD(1) * t550 + t538;
t387 = mrSges(6,1) * t420 + mrSges(5,2) * t444 - mrSges(5,3) * t424 - mrSges(6,3) * t422 + pkin(5) * t400 - qJ(5) * t397 - t568 * t456 + t577 * t457 + t564 * t494 + t576 * t523 - t573 * t526 + t542;
t386 = -mrSges(5,1) * t444 - mrSges(6,1) * t418 + mrSges(6,2) * t422 + mrSges(5,3) * t425 - pkin(4) * t397 - pkin(5) * t545 - pkin(9) * t553 - t533 * t402 - t529 * t403 + t578 * t456 + t568 * t457 + t564 * t495 + t575 * t523 + t563 * t526;
t385 = mrSges(4,2) * t499 - mrSges(4,3) * t459 + Ifges(4,1) * t501 + Ifges(4,4) * t500 + Ifges(4,5) * qJDD(3) - pkin(8) * t393 - qJD(3) * t492 - t530 * t386 + t387 * t571 + t511 * t491;
t384 = -mrSges(4,1) * t499 + mrSges(4,3) * t460 + Ifges(4,4) * t501 + Ifges(4,2) * t500 + Ifges(4,6) * qJDD(3) - pkin(3) * t540 + pkin(8) * t554 + qJD(3) * t493 + t386 * t571 + t530 * t387 - t512 * t491;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t557 - mrSges(2,2) * t552 + t527 * (mrSges(3,2) * t509 - mrSges(3,3) * t502 + t534 * t385 - t531 * t384 - pkin(7) * t565 + (Ifges(3,1) * t527 + Ifges(3,4) * t528) * qJDD(1)) + t528 * (-mrSges(3,1) * t509 + mrSges(3,3) * t503 + t531 * t385 + t534 * t384 - pkin(2) * t538 + pkin(7) * t555 + (Ifges(3,4) * t527 + Ifges(3,2) * t528) * qJDD(1)) - pkin(1) * t394 + qJ(2) * ((m(3) * t503 + t528 * t547 + t555) * t528 + (-m(3) * t502 + t527 * t547 - t565) * t527); t394; mrSges(4,1) * t459 - mrSges(4,2) * t460 + Ifges(4,5) * t501 + Ifges(4,6) * t500 + Ifges(4,3) * qJDD(3) + pkin(3) * t393 + t512 * t492 - t511 * t493 + t539; t539; t399; t542;];
tauJ  = t1;
