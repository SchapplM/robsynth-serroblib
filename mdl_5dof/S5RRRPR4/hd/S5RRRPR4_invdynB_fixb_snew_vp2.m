% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRPR4
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRPR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR4_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR4_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR4_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR4_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR4_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:11
% EndTime: 2019-12-31 21:11:14
% DurationCPUTime: 2.02s
% Computational Cost: add. (24586->253), mult. (31242->310), div. (0->0), fcn. (16014->8), ass. (0->102)
t554 = Ifges(4,1) + Ifges(5,1);
t548 = Ifges(4,4) - Ifges(5,5);
t547 = Ifges(4,5) + Ifges(5,4);
t553 = Ifges(4,2) + Ifges(5,3);
t546 = Ifges(4,6) - Ifges(5,6);
t552 = Ifges(4,3) + Ifges(5,2);
t551 = 2 * qJD(4);
t513 = qJD(1) + qJD(2);
t509 = t513 ^ 2;
t550 = t509 * pkin(7);
t549 = mrSges(4,3) + mrSges(5,2);
t519 = sin(qJ(3));
t545 = t513 * t519;
t523 = cos(qJ(3));
t544 = t513 * t523;
t521 = sin(qJ(1));
t525 = cos(qJ(1));
t503 = t521 * g(1) - t525 * g(2);
t496 = qJDD(1) * pkin(1) + t503;
t504 = -t525 * g(1) - t521 * g(2);
t527 = qJD(1) ^ 2;
t497 = -t527 * pkin(1) + t504;
t520 = sin(qJ(2));
t524 = cos(qJ(2));
t466 = t520 * t496 + t524 * t497;
t511 = qJDD(1) + qJDD(2);
t464 = -t509 * pkin(2) + t511 * pkin(7) + t466;
t454 = -t519 * g(3) + t523 * t464;
t487 = (-mrSges(4,1) * t523 + mrSges(4,2) * t519) * t513;
t489 = -qJD(3) * t545 + t523 * t511;
t498 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t545;
t485 = (-pkin(3) * t523 - qJ(4) * t519) * t513;
t526 = qJD(3) ^ 2;
t447 = -t526 * pkin(3) + qJDD(3) * qJ(4) + qJD(3) * t551 + t485 * t544 + t454;
t486 = (-mrSges(5,1) * t523 - mrSges(5,3) * t519) * t513;
t499 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t545;
t502 = -qJD(3) * pkin(4) - pkin(8) * t545;
t517 = t523 ^ 2;
t443 = -t517 * t509 * pkin(4) - t489 * pkin(8) + qJD(3) * t502 + t447;
t453 = -t523 * g(3) - t519 * t464;
t448 = -qJDD(3) * pkin(3) - t526 * qJ(4) + t485 * t545 + qJDD(4) - t453;
t539 = qJD(3) * t523;
t538 = t513 * t539;
t488 = t519 * t511 + t538;
t444 = (-t488 + t538) * pkin(8) + (-t509 * t519 * t523 - qJDD(3)) * pkin(4) + t448;
t518 = sin(qJ(5));
t522 = cos(qJ(5));
t439 = -t518 * t443 + t522 * t444;
t478 = (-t518 * t519 - t522 * t523) * t513;
t452 = t478 * qJD(5) + t522 * t488 - t518 * t489;
t479 = (-t518 * t523 + t519 * t522) * t513;
t462 = -t478 * mrSges(6,1) + t479 * mrSges(6,2);
t512 = -qJD(3) + qJD(5);
t467 = -t512 * mrSges(6,2) + t478 * mrSges(6,3);
t510 = -qJDD(3) + qJDD(5);
t437 = m(6) * t439 + t510 * mrSges(6,1) - t452 * mrSges(6,3) - t479 * t462 + t512 * t467;
t440 = t522 * t443 + t518 * t444;
t451 = -t479 * qJD(5) - t518 * t488 - t522 * t489;
t468 = t512 * mrSges(6,1) - t479 * mrSges(6,3);
t438 = m(6) * t440 - t510 * mrSges(6,2) + t451 * mrSges(6,3) + t478 * t462 - t512 * t468;
t534 = -t518 * t437 + t522 * t438;
t530 = m(5) * t447 + qJDD(3) * mrSges(5,3) + qJD(3) * t499 + t486 * t544 + t534;
t428 = m(4) * t454 - qJDD(3) * mrSges(4,2) - qJD(3) * t498 + t487 * t544 + t549 * t489 + t530;
t500 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t544;
t430 = t522 * t437 + t518 * t438;
t501 = mrSges(5,2) * t544 + qJD(3) * mrSges(5,3);
t529 = -m(5) * t448 + qJDD(3) * mrSges(5,1) + qJD(3) * t501 - t430;
t429 = m(4) * t453 + qJDD(3) * mrSges(4,1) + qJD(3) * t500 + (-t486 - t487) * t545 - t549 * t488 + t529;
t535 = t523 * t428 - t519 * t429;
t423 = m(3) * t466 - t509 * mrSges(3,1) - t511 * mrSges(3,2) + t535;
t465 = t524 * t496 - t520 * t497;
t533 = t511 * pkin(2) + t465;
t531 = -t488 * qJ(4) - t533;
t445 = -t489 * pkin(3) - t550 + (-0.2e1 * qJD(4) * t519 + (pkin(3) * t519 - qJ(4) * t523) * qJD(3)) * t513 + t531;
t442 = (-pkin(8) * t517 + pkin(7)) * t509 + (pkin(3) + pkin(4)) * t489 + (qJ(4) * t539 + (-pkin(3) * qJD(3) + t502 + t551) * t519) * t513 - t531;
t532 = -m(6) * t442 + t451 * mrSges(6,1) - t452 * mrSges(6,2) + t478 * t467 - t479 * t468;
t435 = m(5) * t445 - t489 * mrSges(5,1) - t488 * mrSges(5,3) - t499 * t545 - t501 * t544 + t532;
t463 = -t533 - t550;
t528 = -m(4) * t463 + t489 * mrSges(4,1) - t488 * mrSges(4,2) - t498 * t545 + t500 * t544 - t435;
t432 = m(3) * t465 + t511 * mrSges(3,1) - t509 * mrSges(3,2) + t528;
t420 = t520 * t423 + t524 * t432;
t418 = m(2) * t503 + qJDD(1) * mrSges(2,1) - t527 * mrSges(2,2) + t420;
t536 = t524 * t423 - t520 * t432;
t419 = m(2) * t504 - t527 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t536;
t543 = t525 * t418 + t521 * t419;
t424 = t519 * t428 + t523 * t429;
t542 = (t547 * t519 + t546 * t523) * t513 + t552 * qJD(3);
t541 = (-t548 * t519 - t553 * t523) * t513 - t546 * qJD(3);
t540 = (t554 * t519 + t548 * t523) * t513 + t547 * qJD(3);
t537 = -t521 * t418 + t525 * t419;
t457 = Ifges(6,1) * t479 + Ifges(6,4) * t478 + Ifges(6,5) * t512;
t456 = Ifges(6,4) * t479 + Ifges(6,2) * t478 + Ifges(6,6) * t512;
t455 = Ifges(6,5) * t479 + Ifges(6,6) * t478 + Ifges(6,3) * t512;
t434 = mrSges(6,2) * t442 - mrSges(6,3) * t439 + Ifges(6,1) * t452 + Ifges(6,4) * t451 + Ifges(6,5) * t510 + t478 * t455 - t512 * t456;
t433 = -mrSges(6,1) * t442 + mrSges(6,3) * t440 + Ifges(6,4) * t452 + Ifges(6,2) * t451 + Ifges(6,6) * t510 - t479 * t455 + t512 * t457;
t414 = mrSges(4,2) * t463 + mrSges(5,2) * t448 - mrSges(4,3) * t453 - mrSges(5,3) * t445 - pkin(8) * t430 - qJ(4) * t435 + t541 * qJD(3) + t547 * qJDD(3) - t518 * t433 + t522 * t434 + t554 * t488 + t548 * t489 + t542 * t544;
t413 = -mrSges(4,1) * t463 - mrSges(5,1) * t445 + mrSges(5,2) * t447 + mrSges(4,3) * t454 - pkin(3) * t435 - pkin(4) * t532 - pkin(8) * t534 + t540 * qJD(3) + t546 * qJDD(3) - t522 * t433 - t518 * t434 + t548 * t488 + t553 * t489 - t542 * t545;
t412 = -pkin(3) * t529 + mrSges(3,1) * g(3) - qJ(4) * t530 + Ifges(6,3) * t510 + Ifges(3,6) * t511 + t509 * Ifges(3,5) - t478 * t457 + t479 * t456 - mrSges(4,1) * t453 + mrSges(4,2) * t454 + mrSges(3,3) * t466 - mrSges(5,3) * t447 + mrSges(5,1) * t448 + Ifges(6,6) * t451 + Ifges(6,5) * t452 + mrSges(6,1) * t439 - mrSges(6,2) * t440 + pkin(4) * t430 - pkin(2) * t424 + (-qJ(4) * mrSges(5,2) - t546) * t489 + (pkin(3) * mrSges(5,2) - t547) * t488 - t552 * qJDD(3) + (t540 * t523 + (pkin(3) * t486 + t541) * t519) * t513;
t411 = -mrSges(3,2) * g(3) - mrSges(3,3) * t465 + Ifges(3,5) * t511 - t509 * Ifges(3,6) - pkin(7) * t424 - t519 * t413 + t523 * t414;
t410 = -mrSges(2,2) * g(3) - mrSges(2,3) * t503 + Ifges(2,5) * qJDD(1) - t527 * Ifges(2,6) - pkin(6) * t420 + t524 * t411 - t520 * t412;
t409 = Ifges(2,6) * qJDD(1) + t527 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t504 + t520 * t411 + t524 * t412 - pkin(1) * (-m(3) * g(3) + t424) + pkin(6) * t536;
t1 = [-m(1) * g(1) + t537; -m(1) * g(2) + t543; (-m(1) - m(2) - m(3)) * g(3) + t424; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t543 - t521 * t409 + t525 * t410; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t537 + t525 * t409 + t521 * t410; -mrSges(1,1) * g(2) + mrSges(2,1) * t503 + mrSges(3,1) * t465 + mrSges(1,2) * g(1) - mrSges(2,2) * t504 - mrSges(3,2) * t466 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * t511 + pkin(1) * t420 + pkin(2) * t528 + pkin(7) * t535 + t523 * t413 + t519 * t414;];
tauB = t1;
