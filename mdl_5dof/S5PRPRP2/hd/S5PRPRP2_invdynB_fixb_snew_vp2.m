% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRPRP2
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRPRP2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:30:44
% EndTime: 2019-12-05 15:30:48
% DurationCPUTime: 2.00s
% Computational Cost: add. (13627->215), mult. (29219->274), div. (0->0), fcn. (17540->8), ass. (0->99)
t469 = sin(pkin(7));
t471 = cos(pkin(7));
t453 = t469 * g(1) - t471 * g(2);
t454 = -t471 * g(1) - t469 * g(2);
t473 = sin(qJ(2));
t475 = cos(qJ(2));
t432 = t473 * t453 + t475 * t454;
t476 = qJD(2) ^ 2;
t521 = -t476 * pkin(2) + qJDD(2) * qJ(3) + (2 * qJD(2) * qJD(3)) + t432;
t520 = Ifges(5,1) + Ifges(6,1);
t513 = Ifges(5,4) + Ifges(6,4);
t512 = Ifges(5,5) + Ifges(6,5);
t519 = Ifges(5,2) + Ifges(6,2);
t518 = Ifges(5,6) + Ifges(6,6);
t517 = -Ifges(5,3) - Ifges(6,3);
t467 = -g(3) + qJDD(1);
t468 = sin(pkin(8));
t470 = cos(pkin(8));
t416 = t470 * t467 - t521 * t468;
t483 = -pkin(3) * t470 - pkin(6) * t468;
t452 = t483 * qJD(2);
t501 = t468 * qJD(2);
t414 = t452 * t501 - t416;
t472 = sin(qJ(4));
t474 = cos(qJ(4));
t498 = qJD(2) * qJD(4);
t444 = (-qJDD(2) * t472 - t474 * t498) * t468;
t445 = (qJDD(2) * t474 - t472 * t498) * t468;
t500 = t470 * qJD(2);
t457 = qJD(4) - t500;
t492 = t474 * t501;
t439 = t457 * pkin(4) - qJ(5) * t492;
t509 = t468 ^ 2 * t476;
t496 = t472 ^ 2 * t509;
t410 = -t444 * pkin(4) - qJ(5) * t496 + t439 * t492 + qJDD(5) + t414;
t490 = m(6) * t410 - t444 * mrSges(6,1);
t514 = -mrSges(5,2) - mrSges(6,2);
t516 = -m(5) * t414 + t444 * mrSges(5,1) + t514 * t445 - t490;
t515 = t470 ^ 2;
t510 = mrSges(4,2) * t468;
t417 = t468 * t467 + t521 * t470;
t450 = (-mrSges(4,1) * t470 + t510) * qJD(2);
t415 = t452 * t500 + t417;
t431 = t475 * t453 - t473 * t454;
t479 = -t476 * qJ(3) + qJDD(3) - t431;
t420 = (-pkin(2) + t483) * qJDD(2) + t479;
t419 = t474 * t420;
t411 = -t472 * t415 + t419;
t493 = t472 * t501;
t438 = -t457 * mrSges(5,2) - mrSges(5,3) * t493;
t497 = t470 * qJDD(2);
t456 = qJDD(4) - t497;
t442 = (mrSges(6,1) * t472 + mrSges(6,2) * t474) * t501;
t482 = (-t442 - (mrSges(5,1) * t472 + mrSges(5,2) * t474) * t501) * t501;
t484 = -0.2e1 * qJD(5) * t501;
t407 = t474 * t484 + t456 * pkin(4) - t445 * qJ(5) + t419 + (-pkin(4) * t474 * t509 - qJ(5) * t457 * t501 - t415) * t472;
t437 = -t457 * mrSges(6,2) - mrSges(6,3) * t493;
t495 = m(6) * t407 + t456 * mrSges(6,1) + t457 * t437;
t401 = m(5) * t411 + t456 * mrSges(5,1) + t457 * t438 + (-mrSges(5,3) - mrSges(6,3)) * t445 + t474 * t482 + t495;
t412 = t474 * t415 + t472 * t420;
t440 = t457 * mrSges(6,1) - mrSges(6,3) * t492;
t503 = -t457 * mrSges(5,1) + mrSges(5,3) * t492 - t440;
t409 = -pkin(4) * t496 + t444 * qJ(5) - t457 * t439 + t472 * t484 + t412;
t507 = m(6) * t409 + t444 * mrSges(6,3);
t402 = m(5) * t412 + t444 * mrSges(5,3) + t514 * t456 + t503 * t457 + t472 * t482 + t507;
t486 = -t472 * t401 + t474 * t402;
t502 = qJDD(2) * mrSges(4,3);
t397 = m(4) * t417 + (qJD(2) * t450 + t502) * t470 + t486;
t478 = t503 * t474 + (-t437 - t438) * t472;
t404 = m(4) * t416 + (-t502 + (-t450 + t478) * qJD(2)) * t468 + t516;
t487 = t470 * t397 - t468 * t404;
t390 = m(3) * t432 - t476 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t487;
t399 = t474 * t401 + t472 * t402;
t423 = -qJDD(2) * pkin(2) + t479;
t477 = -m(4) * t423 + mrSges(4,1) * t497 - t399 + (t476 * t515 + t509) * mrSges(4,3);
t394 = m(3) * t431 - t476 * mrSges(3,2) + (mrSges(3,1) - t510) * qJDD(2) + t477;
t386 = t473 * t390 + t475 * t394;
t384 = m(2) * t453 + t386;
t488 = t475 * t390 - t473 * t394;
t385 = m(2) * t454 + t488;
t508 = t471 * t384 + t469 * t385;
t391 = t468 * t397 + t470 * t404;
t506 = (t518 * t472 - t512 * t474) * t501 + t517 * t457;
t505 = (t519 * t472 - t513 * t474) * t501 - t518 * t457;
t504 = (t513 * t472 - t520 * t474) * t501 - t512 * t457;
t494 = m(3) * t467 + t391;
t489 = -t469 * t384 + t471 * t385;
t481 = Ifges(4,5) * t468 + Ifges(4,6) * t470;
t451 = t481 * qJD(2);
t405 = -t445 * mrSges(6,3) - t442 * t492 + t495;
t398 = mrSges(5,2) * t414 + mrSges(6,2) * t410 - mrSges(5,3) * t411 - mrSges(6,3) * t407 - qJ(5) * t405 + t513 * t444 + t520 * t445 + t512 * t456 + t505 * t457 + t506 * t493;
t392 = -mrSges(5,1) * t414 + mrSges(5,3) * t412 - mrSges(6,1) * t410 + mrSges(6,3) * t409 - pkin(4) * t490 + qJ(5) * t507 + (-qJ(5) * t440 - t504) * t457 + (-qJ(5) * mrSges(6,2) + t518) * t456 + (-pkin(4) * mrSges(6,2) + t513) * t445 + t519 * t444 + ((-pkin(4) * t437 - qJ(5) * t442) * t472 + (-pkin(4) * t440 + t506) * t474) * t501;
t387 = Ifges(4,2) * t497 - mrSges(4,1) * t423 - mrSges(5,1) * t411 - mrSges(6,1) * t407 + mrSges(5,2) * t412 + mrSges(6,2) * t409 + mrSges(4,3) * t417 - pkin(3) * t399 - pkin(4) * t405 + t517 * t456 - t512 * t445 - t518 * t444 + (Ifges(4,4) * qJDD(2) + (t504 * t472 + t505 * t474 - t451) * qJD(2)) * t468;
t380 = t451 * t500 + mrSges(4,2) * t423 - mrSges(4,3) * t416 - pkin(6) * t399 - t472 * t392 + t474 * t398 + (Ifges(4,1) * t468 + Ifges(4,4) * t470) * qJDD(2);
t379 = t476 * Ifges(3,5) - mrSges(3,1) * t467 + mrSges(3,3) * t432 - mrSges(4,1) * t416 + mrSges(4,2) * t417 - t472 * t398 - t474 * t392 - pkin(3) * t516 - pkin(6) * t486 - pkin(2) * t391 + (Ifges(3,6) - t481) * qJDD(2) + (Ifges(4,4) * t515 * qJD(2) + (-pkin(3) * t478 + (-Ifges(4,4) * t468 + (Ifges(4,1) - Ifges(4,2)) * t470) * qJD(2)) * t468) * qJD(2);
t378 = mrSges(3,2) * t467 - mrSges(3,3) * t431 + Ifges(3,5) * qJDD(2) - t476 * Ifges(3,6) - qJ(3) * t391 + t470 * t380 - t468 * t387;
t377 = mrSges(2,2) * t467 - mrSges(2,3) * t453 - pkin(5) * t386 + t475 * t378 - t473 * t379;
t376 = -mrSges(2,1) * t467 + mrSges(2,3) * t454 - pkin(1) * t494 + pkin(5) * t488 + t473 * t378 + t475 * t379;
t1 = [-m(1) * g(1) + t489; -m(1) * g(2) + t508; -m(1) * g(3) + m(2) * t467 + t494; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t508 - t469 * t376 + t471 * t377; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t489 + t471 * t376 + t469 * t377; pkin(1) * t386 + mrSges(2,1) * t453 - mrSges(2,2) * t454 + t470 * t387 + pkin(2) * (-qJDD(2) * t510 + t477) + qJ(3) * t487 + mrSges(3,1) * t431 - mrSges(3,2) * t432 + t468 * t380 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(3,3) * qJDD(2);];
tauB = t1;
