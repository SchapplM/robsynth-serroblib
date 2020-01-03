% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPRP1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP1_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP1_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP1_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:25:25
% EndTime: 2020-01-03 11:25:28
% DurationCPUTime: 2.12s
% Computational Cost: add. (15273->226), mult. (31695->283), div. (0->0), fcn. (17540->8), ass. (0->101)
t491 = sin(qJ(1));
t493 = cos(qJ(1));
t472 = -t493 * g(2) - t491 * g(3);
t468 = qJDD(1) * pkin(1) + t472;
t471 = -t491 * g(2) + t493 * g(3);
t494 = qJD(1) ^ 2;
t469 = -t494 * pkin(1) + t471;
t487 = sin(pkin(7));
t489 = cos(pkin(7));
t441 = t487 * t468 + t489 * t469;
t539 = -t494 * pkin(2) + qJDD(1) * qJ(3) + (2 * qJD(1) * qJD(3)) + t441;
t538 = Ifges(5,1) + Ifges(6,1);
t531 = Ifges(5,4) + Ifges(6,4);
t530 = Ifges(5,5) + Ifges(6,5);
t537 = Ifges(5,2) + Ifges(6,2);
t536 = Ifges(5,6) + Ifges(6,6);
t535 = -Ifges(5,3) - Ifges(6,3);
t485 = -g(1) + qJDD(2);
t486 = sin(pkin(8));
t488 = cos(pkin(8));
t431 = t488 * t485 - t539 * t486;
t501 = -pkin(3) * t488 - pkin(6) * t486;
t467 = t501 * qJD(1);
t519 = t486 * qJD(1);
t429 = t467 * t519 - t431;
t490 = sin(qJ(4));
t492 = cos(qJ(4));
t516 = qJD(1) * qJD(4);
t459 = (-qJDD(1) * t490 - t492 * t516) * t486;
t460 = (qJDD(1) * t492 - t490 * t516) * t486;
t518 = t488 * qJD(1);
t474 = qJD(4) - t518;
t510 = t492 * t519;
t454 = t474 * pkin(4) - qJ(5) * t510;
t527 = t486 ^ 2 * t494;
t514 = t490 ^ 2 * t527;
t427 = -t459 * pkin(4) - qJ(5) * t514 + t454 * t510 + qJDD(5) + t429;
t508 = m(6) * t427 - t459 * mrSges(6,1);
t532 = -mrSges(5,2) - mrSges(6,2);
t534 = -m(5) * t429 + t459 * mrSges(5,1) + t532 * t460 - t508;
t533 = t488 ^ 2;
t528 = mrSges(4,2) * t486;
t432 = t486 * t485 + t539 * t488;
t465 = (-mrSges(4,1) * t488 + t528) * qJD(1);
t430 = t467 * t518 + t432;
t440 = t489 * t468 - t487 * t469;
t497 = -t494 * qJ(3) + qJDD(3) - t440;
t435 = (-pkin(2) + t501) * qJDD(1) + t497;
t434 = t492 * t435;
t425 = -t490 * t430 + t434;
t511 = t490 * t519;
t453 = -t474 * mrSges(5,2) - mrSges(5,3) * t511;
t515 = t488 * qJDD(1);
t473 = qJDD(4) - t515;
t457 = (mrSges(6,1) * t490 + mrSges(6,2) * t492) * t519;
t500 = (-t457 - (mrSges(5,1) * t490 + mrSges(5,2) * t492) * t519) * t519;
t502 = -0.2e1 * qJD(5) * t519;
t422 = t492 * t502 + t473 * pkin(4) - t460 * qJ(5) + t434 + (-pkin(4) * t492 * t527 - qJ(5) * t474 * t519 - t430) * t490;
t452 = -t474 * mrSges(6,2) - mrSges(6,3) * t511;
t513 = m(6) * t422 + t473 * mrSges(6,1) + t474 * t452;
t416 = m(5) * t425 + t473 * mrSges(5,1) + t474 * t453 + (-mrSges(5,3) - mrSges(6,3)) * t460 + t492 * t500 + t513;
t426 = t492 * t430 + t490 * t435;
t455 = t474 * mrSges(6,1) - mrSges(6,3) * t510;
t521 = -t474 * mrSges(5,1) + mrSges(5,3) * t510 - t455;
t424 = -pkin(4) * t514 + t459 * qJ(5) - t474 * t454 + t490 * t502 + t426;
t525 = m(6) * t424 + t459 * mrSges(6,3);
t417 = m(5) * t426 + t459 * mrSges(5,3) + t532 * t473 + t521 * t474 + t490 * t500 + t525;
t504 = -t490 * t416 + t492 * t417;
t520 = qJDD(1) * mrSges(4,3);
t412 = m(4) * t432 + (qJD(1) * t465 + t520) * t488 + t504;
t496 = t521 * t492 + (-t452 - t453) * t490;
t419 = m(4) * t431 + (-t520 + (-t465 + t496) * qJD(1)) * t486 + t534;
t505 = t488 * t412 - t486 * t419;
t405 = m(3) * t441 - t494 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t505;
t414 = t492 * t416 + t490 * t417;
t438 = -qJDD(1) * pkin(2) + t497;
t495 = -m(4) * t438 + mrSges(4,1) * t515 - t414 + (t494 * t533 + t527) * mrSges(4,3);
t409 = m(3) * t440 - t494 * mrSges(3,2) + (mrSges(3,1) - t528) * qJDD(1) + t495;
t506 = t489 * t405 - t487 * t409;
t399 = m(2) * t471 - t494 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t506;
t401 = t487 * t405 + t489 * t409;
t400 = m(2) * t472 + qJDD(1) * mrSges(2,1) - t494 * mrSges(2,2) + t401;
t526 = t491 * t399 + t493 * t400;
t406 = t486 * t412 + t488 * t419;
t524 = (t536 * t490 - t530 * t492) * t519 + t535 * t474;
t523 = (t537 * t490 - t531 * t492) * t519 - t536 * t474;
t522 = (t531 * t490 - t538 * t492) * t519 - t530 * t474;
t512 = m(3) * t485 + t406;
t507 = -t493 * t399 + t491 * t400;
t499 = Ifges(4,5) * t486 + Ifges(4,6) * t488;
t466 = t499 * qJD(1);
t420 = -t460 * mrSges(6,3) - t457 * t510 + t513;
t413 = mrSges(5,2) * t429 + mrSges(6,2) * t427 - mrSges(5,3) * t425 - mrSges(6,3) * t422 - qJ(5) * t420 + t531 * t459 + t538 * t460 + t530 * t473 + t523 * t474 + t524 * t511;
t407 = -mrSges(5,1) * t429 + mrSges(5,3) * t426 - mrSges(6,1) * t427 + mrSges(6,3) * t424 - pkin(4) * t508 + qJ(5) * t525 + (-qJ(5) * t455 - t522) * t474 + (-qJ(5) * mrSges(6,2) + t536) * t473 + (-pkin(4) * mrSges(6,2) + t531) * t460 + t537 * t459 + ((-pkin(4) * t452 - qJ(5) * t457) * t490 + (-pkin(4) * t455 + t524) * t492) * t519;
t402 = Ifges(4,2) * t515 - mrSges(4,1) * t438 - mrSges(5,1) * t425 - mrSges(6,1) * t422 + mrSges(5,2) * t426 + mrSges(6,2) * t424 + mrSges(4,3) * t432 - pkin(3) * t414 - pkin(4) * t420 + t535 * t473 - t530 * t460 - t536 * t459 + (Ifges(4,4) * qJDD(1) + (t522 * t490 + t523 * t492 - t466) * qJD(1)) * t486;
t395 = t466 * t518 + mrSges(4,2) * t438 - mrSges(4,3) * t431 - pkin(6) * t414 - t490 * t407 + t492 * t413 + (Ifges(4,1) * t486 + Ifges(4,4) * t488) * qJDD(1);
t394 = t494 * Ifges(3,5) - mrSges(3,1) * t485 + mrSges(3,3) * t441 - mrSges(4,1) * t431 + mrSges(4,2) * t432 - t490 * t413 - t492 * t407 - pkin(3) * t534 - pkin(6) * t504 - pkin(2) * t406 + (Ifges(3,6) - t499) * qJDD(1) + (Ifges(4,4) * t533 * qJD(1) + (-pkin(3) * t496 + (-Ifges(4,4) * t486 + (Ifges(4,1) - Ifges(4,2)) * t488) * qJD(1)) * t486) * qJD(1);
t393 = mrSges(3,2) * t485 - mrSges(3,3) * t440 + Ifges(3,5) * qJDD(1) - t494 * Ifges(3,6) - qJ(3) * t406 + t488 * t395 - t486 * t402;
t392 = -mrSges(2,2) * g(1) - mrSges(2,3) * t472 + Ifges(2,5) * qJDD(1) - t494 * Ifges(2,6) - qJ(2) * t401 + t489 * t393 - t487 * t394;
t391 = mrSges(2,1) * g(1) + mrSges(2,3) * t471 + t494 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t512 + qJ(2) * t506 + t487 * t393 + t489 * t394;
t1 = [(-m(1) - m(2)) * g(1) + t512; -m(1) * g(2) + t526; -m(1) * g(3) + t507; pkin(1) * t401 + t486 * t395 + t488 * t402 + pkin(2) * t495 + qJ(3) * t505 + mrSges(3,1) * t440 - mrSges(3,2) * t441 + mrSges(2,1) * t472 - mrSges(2,2) * t471 - mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + (-pkin(2) * t528 + Ifges(2,3) + Ifges(3,3)) * qJDD(1); mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t507 + t493 * t391 + t491 * t392; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t526 + t491 * t391 - t493 * t392;];
tauB = t1;
