% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRRR12
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRRR12_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR12_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR12_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR12_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR12_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR12_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR12_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:12:41
% EndTime: 2019-12-31 19:12:44
% DurationCPUTime: 2.41s
% Computational Cost: add. (24471->265), mult. (47907->327), div. (0->0), fcn. (29818->8), ass. (0->104)
t494 = sin(qJ(1));
t498 = cos(qJ(1));
t479 = -t498 * g(1) - t494 * g(2);
t506 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t479;
t492 = sin(qJ(4));
t493 = sin(qJ(3));
t496 = cos(qJ(4));
t497 = cos(qJ(3));
t466 = (t492 * t497 + t493 * t496) * qJD(1);
t522 = -pkin(1) - pkin(6);
t521 = mrSges(2,1) - mrSges(3,2);
t520 = Ifges(2,5) - Ifges(3,4);
t519 = (-Ifges(2,6) + Ifges(3,5));
t478 = t494 * g(1) - t498 * g(2);
t499 = qJD(1) ^ 2;
t505 = -t499 * qJ(2) + qJDD(2) - t478;
t459 = t522 * qJDD(1) + t505;
t449 = t493 * g(3) + t497 * t459;
t515 = qJD(1) * qJD(3);
t513 = t493 * t515;
t474 = t497 * qJDD(1) - t513;
t432 = (-t474 - t513) * pkin(7) + (-t493 * t497 * t499 + qJDD(3)) * pkin(3) + t449;
t450 = -t497 * g(3) + t493 * t459;
t473 = -t493 * qJDD(1) - t497 * t515;
t516 = qJD(1) * t497;
t477 = qJD(3) * pkin(3) - pkin(7) * t516;
t488 = t493 ^ 2;
t433 = -t488 * t499 * pkin(3) + t473 * pkin(7) - qJD(3) * t477 + t450;
t422 = t492 * t432 + t496 * t433;
t467 = (-t492 * t493 + t496 * t497) * qJD(1);
t439 = -t467 * qJD(4) + t496 * t473 - t492 * t474;
t447 = t466 * mrSges(5,1) + t467 * mrSges(5,2);
t485 = qJD(3) + qJD(4);
t457 = t485 * mrSges(5,1) - t467 * mrSges(5,3);
t484 = qJDD(3) + qJDD(4);
t436 = -t473 * pkin(3) + t477 * t516 + (-pkin(7) * t488 + t522) * t499 + t506;
t440 = -t466 * qJD(4) + t492 * t473 + t496 * t474;
t418 = (t466 * t485 - t440) * pkin(8) + (t467 * t485 - t439) * pkin(4) + t436;
t448 = t466 * pkin(4) - t467 * pkin(8);
t483 = t485 ^ 2;
t420 = -t483 * pkin(4) + t484 * pkin(8) - t466 * t448 + t422;
t491 = sin(qJ(5));
t495 = cos(qJ(5));
t416 = t495 * t418 - t491 * t420;
t451 = -t491 * t467 + t495 * t485;
t425 = t451 * qJD(5) + t495 * t440 + t491 * t484;
t452 = t495 * t467 + t491 * t485;
t434 = -t451 * mrSges(6,1) + t452 * mrSges(6,2);
t438 = qJDD(5) - t439;
t462 = qJD(5) + t466;
t441 = -t462 * mrSges(6,2) + t451 * mrSges(6,3);
t414 = m(6) * t416 + t438 * mrSges(6,1) - t425 * mrSges(6,3) - t452 * t434 + t462 * t441;
t417 = t491 * t418 + t495 * t420;
t424 = -t452 * qJD(5) - t491 * t440 + t495 * t484;
t442 = t462 * mrSges(6,1) - t452 * mrSges(6,3);
t415 = m(6) * t417 - t438 * mrSges(6,2) + t424 * mrSges(6,3) + t451 * t434 - t462 * t442;
t509 = -t491 * t414 + t495 * t415;
t405 = m(5) * t422 - t484 * mrSges(5,2) + t439 * mrSges(5,3) - t466 * t447 - t485 * t457 + t509;
t421 = t496 * t432 - t492 * t433;
t456 = -t485 * mrSges(5,2) - t466 * mrSges(5,3);
t419 = -t484 * pkin(4) - t483 * pkin(8) + t467 * t448 - t421;
t503 = -m(6) * t419 + t424 * mrSges(6,1) - t425 * mrSges(6,2) + t451 * t441 - t452 * t442;
t410 = m(5) * t421 + t484 * mrSges(5,1) - t440 * mrSges(5,3) - t467 * t447 + t485 * t456 + t503;
t399 = t492 * t405 + t496 * t410;
t472 = (mrSges(4,1) * t493 + mrSges(4,2) * t497) * qJD(1);
t517 = qJD(1) * t493;
t475 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t517;
t397 = m(4) * t449 + qJDD(3) * mrSges(4,1) - t474 * mrSges(4,3) + qJD(3) * t475 - t472 * t516 + t399;
t476 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t516;
t510 = t496 * t405 - t492 * t410;
t398 = m(4) * t450 - qJDD(3) * mrSges(4,2) + t473 * mrSges(4,3) - qJD(3) * t476 - t472 * t517 + t510;
t393 = t497 * t397 + t493 * t398;
t461 = -qJDD(1) * pkin(1) + t505;
t504 = -m(3) * t461 + (t499 * mrSges(3,3)) - t393;
t391 = m(2) * t478 - (t499 * mrSges(2,2)) + t521 * qJDD(1) + t504;
t460 = t499 * pkin(1) - t506;
t458 = t522 * t499 + t506;
t406 = t495 * t414 + t491 * t415;
t502 = m(5) * t436 - t439 * mrSges(5,1) + t440 * mrSges(5,2) + t466 * t456 + t467 * t457 + t406;
t501 = -m(4) * t458 + t473 * mrSges(4,1) - t474 * mrSges(4,2) - t475 * t517 - t476 * t516 - t502;
t500 = -m(3) * t460 + (t499 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t501;
t402 = m(2) * t479 - (t499 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t500;
t518 = t498 * t391 + t494 * t402;
t512 = -t494 * t391 + t498 * t402;
t511 = -t493 * t397 + t497 * t398;
t465 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t497 - Ifges(4,4) * t493) * qJD(1);
t464 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t497 - Ifges(4,2) * t493) * qJD(1);
t463 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t497 - Ifges(4,6) * t493) * qJD(1);
t445 = Ifges(5,1) * t467 - Ifges(5,4) * t466 + Ifges(5,5) * t485;
t444 = Ifges(5,4) * t467 - Ifges(5,2) * t466 + Ifges(5,6) * t485;
t443 = Ifges(5,5) * t467 - Ifges(5,6) * t466 + Ifges(5,3) * t485;
t428 = Ifges(6,1) * t452 + Ifges(6,4) * t451 + Ifges(6,5) * t462;
t427 = Ifges(6,4) * t452 + Ifges(6,2) * t451 + Ifges(6,6) * t462;
t426 = Ifges(6,5) * t452 + Ifges(6,6) * t451 + Ifges(6,3) * t462;
t408 = mrSges(6,2) * t419 - mrSges(6,3) * t416 + Ifges(6,1) * t425 + Ifges(6,4) * t424 + Ifges(6,5) * t438 + t451 * t426 - t462 * t427;
t407 = -mrSges(6,1) * t419 + mrSges(6,3) * t417 + Ifges(6,4) * t425 + Ifges(6,2) * t424 + Ifges(6,6) * t438 - t452 * t426 + t462 * t428;
t395 = -mrSges(5,1) * t436 - mrSges(6,1) * t416 + mrSges(6,2) * t417 + mrSges(5,3) * t422 + Ifges(5,4) * t440 - Ifges(6,5) * t425 + Ifges(5,2) * t439 + Ifges(5,6) * t484 - Ifges(6,6) * t424 - Ifges(6,3) * t438 - pkin(4) * t406 - t452 * t427 + t451 * t428 - t467 * t443 + t485 * t445;
t394 = mrSges(5,2) * t436 - mrSges(5,3) * t421 + Ifges(5,1) * t440 + Ifges(5,4) * t439 + Ifges(5,5) * t484 - pkin(8) * t406 - t491 * t407 + t495 * t408 - t466 * t443 - t485 * t444;
t392 = -m(3) * g(3) + t511;
t389 = mrSges(4,2) * t458 - mrSges(4,3) * t449 + Ifges(4,1) * t474 + Ifges(4,4) * t473 + Ifges(4,5) * qJDD(3) - pkin(7) * t399 - qJD(3) * t464 + t496 * t394 - t492 * t395 - t463 * t517;
t388 = -mrSges(4,1) * t458 + mrSges(4,3) * t450 + Ifges(4,4) * t474 + Ifges(4,2) * t473 + Ifges(4,6) * qJDD(3) - pkin(3) * t502 + pkin(7) * t510 + qJD(3) * t465 + t492 * t394 + t496 * t395 - t463 * t516;
t387 = (t497 * t464 + t493 * t465) * qJD(1) + pkin(8) * t509 + pkin(4) * t503 + t491 * t408 + t495 * t407 - mrSges(2,3) * t478 + Ifges(5,3) * t484 + t466 * t445 + t467 * t444 + Ifges(4,6) * t473 + Ifges(4,5) * t474 + mrSges(3,1) * t461 + mrSges(4,1) * t449 - mrSges(4,2) * t450 + Ifges(5,6) * t439 + Ifges(5,5) * t440 + mrSges(5,1) * t421 - mrSges(5,2) * t422 + pkin(3) * t399 + Ifges(4,3) * qJDD(3) + pkin(2) * t393 - qJ(2) * t392 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (t519 * t499) + t520 * qJDD(1);
t386 = -mrSges(3,1) * t460 + mrSges(2,3) * t479 - pkin(1) * t392 - pkin(2) * t501 - pkin(6) * t511 + t521 * g(3) - t519 * qJDD(1) - t497 * t388 - t493 * t389 + t520 * t499;
t1 = [-m(1) * g(1) + t512; -m(1) * g(2) + t518; (-m(1) - m(2) - m(3)) * g(3) + t511; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t518 - t494 * t386 + t498 * t387; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t512 + t498 * t386 + t494 * t387; pkin(1) * t504 + qJ(2) * t500 + t497 * t389 - t493 * t388 - pkin(6) * t393 + mrSges(2,1) * t478 - mrSges(2,2) * t479 + mrSges(3,2) * t461 - mrSges(3,3) * t460 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB = t1;
