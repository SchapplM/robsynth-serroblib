% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRPR14
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRPR14_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR14_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR14_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR14_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR14_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR14_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:34:42
% EndTime: 2019-12-31 18:34:46
% DurationCPUTime: 2.24s
% Computational Cost: add. (21311->265), mult. (45597->327), div. (0->0), fcn. (27730->8), ass. (0->104)
t489 = sin(qJ(1));
t492 = cos(qJ(1));
t475 = -t492 * g(1) - t489 * g(2);
t501 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t475;
t485 = sin(pkin(8));
t486 = cos(pkin(8));
t488 = sin(qJ(3));
t491 = cos(qJ(3));
t459 = (t485 * t491 + t486 * t488) * qJD(1);
t519 = 2 * qJD(4);
t518 = -pkin(1) - pkin(6);
t517 = mrSges(2,1) - mrSges(3,2);
t516 = Ifges(2,5) - Ifges(3,4);
t515 = (-Ifges(2,6) + Ifges(3,5));
t474 = t489 * g(1) - t492 * g(2);
t494 = qJD(1) ^ 2;
t500 = -t494 * qJ(2) + qJDD(2) - t474;
t454 = t518 * qJDD(1) + t500;
t444 = t488 * g(3) + t491 * t454;
t511 = qJD(1) * qJD(3);
t509 = t488 * t511;
t470 = t491 * qJDD(1) - t509;
t427 = (-t470 - t509) * qJ(4) + (-t488 * t491 * t494 + qJDD(3)) * pkin(3) + t444;
t445 = -t491 * g(3) + t488 * t454;
t469 = -t488 * qJDD(1) - t491 * t511;
t512 = qJD(1) * t491;
t472 = qJD(3) * pkin(3) - qJ(4) * t512;
t482 = t488 ^ 2;
t428 = -t482 * t494 * pkin(3) + t469 * qJ(4) - qJD(3) * t472 + t445;
t417 = t485 * t427 + t486 * t428 - t459 * t519;
t513 = qJD(1) * t488;
t460 = -t485 * t513 + t486 * t512;
t438 = t459 * mrSges(5,1) + t460 * mrSges(5,2);
t442 = t486 * t469 - t485 * t470;
t453 = qJD(3) * mrSges(5,1) - t460 * mrSges(5,3);
t439 = t459 * pkin(4) - t460 * pkin(7);
t493 = qJD(3) ^ 2;
t414 = -t493 * pkin(4) + qJDD(3) * pkin(7) - t459 * t439 + t417;
t430 = -t469 * pkin(3) + qJDD(4) + t472 * t512 + (-qJ(4) * t482 + t518) * t494 + t501;
t443 = t485 * t469 + t486 * t470;
t415 = (qJD(3) * t459 - t443) * pkin(7) + (qJD(3) * t460 - t442) * pkin(4) + t430;
t487 = sin(qJ(5));
t490 = cos(qJ(5));
t411 = -t487 * t414 + t490 * t415;
t446 = t490 * qJD(3) - t487 * t460;
t424 = t446 * qJD(5) + t487 * qJDD(3) + t490 * t443;
t447 = t487 * qJD(3) + t490 * t460;
t431 = -t446 * mrSges(6,1) + t447 * mrSges(6,2);
t457 = qJD(5) + t459;
t432 = -t457 * mrSges(6,2) + t446 * mrSges(6,3);
t441 = qJDD(5) - t442;
t409 = m(6) * t411 + t441 * mrSges(6,1) - t424 * mrSges(6,3) - t447 * t431 + t457 * t432;
t412 = t490 * t414 + t487 * t415;
t423 = -t447 * qJD(5) + t490 * qJDD(3) - t487 * t443;
t433 = t457 * mrSges(6,1) - t447 * mrSges(6,3);
t410 = m(6) * t412 - t441 * mrSges(6,2) + t423 * mrSges(6,3) + t446 * t431 - t457 * t433;
t505 = -t487 * t409 + t490 * t410;
t400 = m(5) * t417 - qJDD(3) * mrSges(5,2) + t442 * mrSges(5,3) - qJD(3) * t453 - t459 * t438 + t505;
t504 = -t486 * t427 + t485 * t428;
t416 = -0.2e1 * qJD(4) * t460 - t504;
t452 = -qJD(3) * mrSges(5,2) - t459 * mrSges(5,3);
t413 = -qJDD(3) * pkin(4) - t493 * pkin(7) + (t519 + t439) * t460 + t504;
t498 = -m(6) * t413 + t423 * mrSges(6,1) - t424 * mrSges(6,2) + t446 * t432 - t447 * t433;
t405 = m(5) * t416 + qJDD(3) * mrSges(5,1) - t443 * mrSges(5,3) + qJD(3) * t452 - t460 * t438 + t498;
t394 = t485 * t400 + t486 * t405;
t468 = (mrSges(4,1) * t488 + mrSges(4,2) * t491) * qJD(1);
t471 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t513;
t392 = m(4) * t444 + qJDD(3) * mrSges(4,1) - t470 * mrSges(4,3) + qJD(3) * t471 - t468 * t512 + t394;
t473 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t512;
t506 = t486 * t400 - t485 * t405;
t393 = m(4) * t445 - qJDD(3) * mrSges(4,2) + t469 * mrSges(4,3) - qJD(3) * t473 - t468 * t513 + t506;
t388 = t491 * t392 + t488 * t393;
t458 = -qJDD(1) * pkin(1) + t500;
t499 = -m(3) * t458 + (t494 * mrSges(3,3)) - t388;
t386 = m(2) * t474 - (t494 * mrSges(2,2)) + t517 * qJDD(1) + t499;
t455 = t494 * pkin(1) - t501;
t451 = t518 * t494 + t501;
t401 = t490 * t409 + t487 * t410;
t497 = m(5) * t430 - t442 * mrSges(5,1) + t443 * mrSges(5,2) + t459 * t452 + t460 * t453 + t401;
t496 = -m(4) * t451 + t469 * mrSges(4,1) - t470 * mrSges(4,2) - t471 * t513 - t473 * t512 - t497;
t495 = -m(3) * t455 + (t494 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t496;
t397 = m(2) * t475 - (t494 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t495;
t514 = t492 * t386 + t489 * t397;
t508 = -t489 * t386 + t492 * t397;
t507 = -t488 * t392 + t491 * t393;
t463 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t491 - Ifges(4,4) * t488) * qJD(1);
t462 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t491 - Ifges(4,2) * t488) * qJD(1);
t461 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t491 - Ifges(4,6) * t488) * qJD(1);
t436 = Ifges(5,1) * t460 - Ifges(5,4) * t459 + (Ifges(5,5) * qJD(3));
t435 = Ifges(5,4) * t460 - Ifges(5,2) * t459 + (Ifges(5,6) * qJD(3));
t434 = Ifges(5,5) * t460 - Ifges(5,6) * t459 + (Ifges(5,3) * qJD(3));
t420 = Ifges(6,1) * t447 + Ifges(6,4) * t446 + Ifges(6,5) * t457;
t419 = Ifges(6,4) * t447 + Ifges(6,2) * t446 + Ifges(6,6) * t457;
t418 = Ifges(6,5) * t447 + Ifges(6,6) * t446 + Ifges(6,3) * t457;
t403 = mrSges(6,2) * t413 - mrSges(6,3) * t411 + Ifges(6,1) * t424 + Ifges(6,4) * t423 + Ifges(6,5) * t441 + t446 * t418 - t457 * t419;
t402 = -mrSges(6,1) * t413 + mrSges(6,3) * t412 + Ifges(6,4) * t424 + Ifges(6,2) * t423 + Ifges(6,6) * t441 - t447 * t418 + t457 * t420;
t390 = -mrSges(5,1) * t430 - mrSges(6,1) * t411 + mrSges(6,2) * t412 + mrSges(5,3) * t417 + Ifges(5,4) * t443 - Ifges(6,5) * t424 + Ifges(5,2) * t442 + Ifges(5,6) * qJDD(3) - Ifges(6,6) * t423 - Ifges(6,3) * t441 - pkin(4) * t401 + qJD(3) * t436 - t447 * t419 + t446 * t420 - t460 * t434;
t389 = mrSges(5,2) * t430 - mrSges(5,3) * t416 + Ifges(5,1) * t443 + Ifges(5,4) * t442 + Ifges(5,5) * qJDD(3) - pkin(7) * t401 - qJD(3) * t435 - t487 * t402 + t490 * t403 - t459 * t434;
t387 = -m(3) * g(3) + t507;
t384 = mrSges(4,2) * t451 - mrSges(4,3) * t444 + Ifges(4,1) * t470 + Ifges(4,4) * t469 + Ifges(4,5) * qJDD(3) - qJ(4) * t394 - qJD(3) * t462 + t486 * t389 - t485 * t390 - t461 * t513;
t383 = -mrSges(4,1) * t451 + mrSges(4,3) * t445 + Ifges(4,4) * t470 + Ifges(4,2) * t469 + Ifges(4,6) * qJDD(3) - pkin(3) * t497 + qJ(4) * t506 + qJD(3) * t463 + t485 * t389 + t486 * t390 - t461 * t512;
t382 = t490 * t402 + pkin(7) * t505 + t487 * t403 - mrSges(2,3) * t474 + t460 * t435 + Ifges(4,6) * t469 + Ifges(4,5) * t470 + mrSges(3,1) * t458 + t459 * t436 + pkin(4) * t498 + Ifges(5,6) * t442 + Ifges(5,5) * t443 + mrSges(4,1) * t444 - mrSges(4,2) * t445 + mrSges(5,1) * t416 - mrSges(5,2) * t417 + pkin(3) * t394 + pkin(2) * t388 - qJ(2) * t387 + (t515 * t494) + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + t516 * qJDD(1) + (t491 * t462 + t488 * t463) * qJD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t381 = -mrSges(3,1) * t455 + mrSges(2,3) * t475 - pkin(1) * t387 - pkin(2) * t496 - pkin(6) * t507 + t517 * g(3) - t515 * qJDD(1) - t491 * t383 - t488 * t384 + t516 * t494;
t1 = [-m(1) * g(1) + t508; -m(1) * g(2) + t514; (-m(1) - m(2) - m(3)) * g(3) + t507; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t514 - t489 * t381 + t492 * t382; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t508 + t492 * t381 + t489 * t382; pkin(1) * t499 + qJ(2) * t495 + t491 * t384 - t488 * t383 - pkin(6) * t388 + mrSges(2,1) * t474 - mrSges(2,2) * t475 + mrSges(3,2) * t458 - mrSges(3,3) * t455 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB = t1;
