% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PPRRP3
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PPRRP3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP3_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP3_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP3_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP3_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:10:42
% EndTime: 2019-12-05 15:10:43
% DurationCPUTime: 1.19s
% Computational Cost: add. (9428->182), mult. (16179->228), div. (0->0), fcn. (9142->8), ass. (0->82)
t507 = Ifges(5,1) + Ifges(6,1);
t500 = Ifges(5,4) - Ifges(6,5);
t499 = -Ifges(5,5) - Ifges(6,4);
t506 = Ifges(5,2) + Ifges(6,3);
t498 = Ifges(5,6) - Ifges(6,6);
t505 = Ifges(5,3) + Ifges(6,2);
t475 = sin(qJ(4));
t477 = cos(qJ(4));
t452 = (-mrSges(6,1) * t477 - mrSges(6,3) * t475) * qJD(3);
t489 = qJD(3) * qJD(4);
t454 = t475 * qJDD(3) + t477 * t489;
t472 = sin(pkin(7));
t474 = cos(pkin(7));
t460 = -t474 * g(1) - t472 * g(2);
t470 = -g(3) + qJDD(1);
t471 = sin(pkin(8));
t473 = cos(pkin(8));
t437 = t473 * t460 + t471 * t470;
t459 = t472 * g(1) - t474 * g(2);
t458 = qJDD(2) - t459;
t476 = sin(qJ(3));
t478 = cos(qJ(3));
t432 = t478 * t437 + t476 * t458;
t480 = qJD(3) ^ 2;
t430 = -t480 * pkin(3) + qJDD(3) * pkin(6) + t432;
t451 = (-pkin(4) * t477 - qJ(5) * t475) * qJD(3);
t479 = qJD(4) ^ 2;
t436 = t471 * t460 - t473 * t470;
t496 = t477 * t436;
t424 = -qJDD(4) * pkin(4) - t479 * qJ(5) - t496 + qJDD(5) + (qJD(3) * t451 + t430) * t475;
t490 = qJD(3) * t477;
t464 = mrSges(6,2) * t490 + qJD(4) * mrSges(6,3);
t483 = -m(6) * t424 + qJDD(4) * mrSges(6,1) + qJD(4) * t464;
t491 = qJD(3) * t475;
t421 = t454 * mrSges(6,2) + t452 * t491 - t483;
t427 = t477 * t430 + t475 * t436;
t423 = -t479 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t451 * t490 + t427;
t426 = -t475 * t430 + t496;
t455 = t477 * qJDD(3) - t475 * t489;
t462 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t491;
t484 = m(6) * t423 + qJDD(4) * mrSges(6,3) + qJD(4) * t462 + t452 * t490;
t492 = -t499 * qJD(4) + (t507 * t475 + t500 * t477) * qJD(3);
t494 = -t498 * qJD(4) + (-t500 * t475 - t506 * t477) * qJD(3);
t504 = -(t494 * t475 + t492 * t477) * qJD(3) + t505 * qJDD(4) - t499 * t454 + t498 * t455 + mrSges(5,1) * t426 - mrSges(6,1) * t424 - mrSges(5,2) * t427 + mrSges(6,3) * t423 - pkin(4) * t421 + qJ(5) * (t455 * mrSges(6,2) + t484);
t501 = mrSges(5,3) + mrSges(6,2);
t453 = (-mrSges(5,1) * t477 + mrSges(5,2) * t475) * qJD(3);
t461 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t491;
t416 = m(5) * t427 - qJDD(4) * mrSges(5,2) - qJD(4) * t461 + t453 * t490 + t501 * t455 + t484;
t463 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t490;
t417 = m(5) * t426 + qJDD(4) * mrSges(5,1) + qJD(4) * t463 - t501 * t454 + (-t452 - t453) * t491 + t483;
t413 = t477 * t416 - t475 * t417;
t410 = m(4) * t432 - t480 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t413;
t431 = -t476 * t437 + t478 * t458;
t429 = -qJDD(3) * pkin(3) - t480 * pkin(6) - t431;
t425 = -t455 * pkin(4) - t454 * qJ(5) + (-0.2e1 * qJD(5) * t475 + (pkin(4) * t475 - qJ(5) * t477) * qJD(4)) * qJD(3) + t429;
t419 = m(6) * t425 - t455 * mrSges(6,1) - t454 * mrSges(6,3) - t462 * t491 - t464 * t490;
t418 = -m(5) * t429 + t455 * mrSges(5,1) - t454 * mrSges(5,2) - t461 * t491 + t463 * t490 - t419;
t414 = m(4) * t431 + qJDD(3) * mrSges(4,1) - t480 * mrSges(4,2) + t418;
t485 = t478 * t410 - t476 * t414;
t403 = m(3) * t437 + t485;
t412 = t475 * t416 + t477 * t417;
t411 = (-m(3) - m(4)) * t436 - t412;
t486 = t473 * t403 - t471 * t411;
t396 = m(2) * t460 + t486;
t405 = t476 * t410 + t478 * t414;
t404 = m(3) * t458 + t405;
t402 = m(2) * t459 - t404;
t495 = t472 * t396 + t474 * t402;
t397 = t471 * t403 + t473 * t411;
t493 = t505 * qJD(4) + (-t499 * t475 + t498 * t477) * qJD(3);
t488 = m(2) * t470 + t397;
t487 = t474 * t396 - t472 * t402;
t406 = -mrSges(5,1) * t429 - mrSges(6,1) * t425 + mrSges(6,2) * t423 + mrSges(5,3) * t427 - pkin(4) * t419 + t492 * qJD(4) + t498 * qJDD(4) + t500 * t454 + t506 * t455 - t493 * t491;
t407 = mrSges(5,2) * t429 + mrSges(6,2) * t424 - mrSges(5,3) * t426 - mrSges(6,3) * t425 - qJ(5) * t419 + t494 * qJD(4) - t499 * qJDD(4) + t507 * t454 + t500 * t455 + t493 * t490;
t481 = mrSges(4,1) * t431 - mrSges(4,2) * t432 + Ifges(4,3) * qJDD(3) + pkin(3) * t418 + pkin(6) * t413 + t477 * t406 + t475 * t407;
t398 = -mrSges(4,1) * t436 + mrSges(4,3) * t432 + t480 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t412 - t504;
t393 = mrSges(4,2) * t436 - mrSges(4,3) * t431 + Ifges(4,5) * qJDD(3) - t480 * Ifges(4,6) - pkin(6) * t412 - t475 * t406 + t477 * t407;
t392 = -mrSges(3,1) * t458 + mrSges(3,3) * t437 - pkin(2) * t405 - t481;
t391 = mrSges(3,2) * t458 + mrSges(3,3) * t436 - pkin(5) * t405 + t478 * t393 - t476 * t398;
t390 = -mrSges(2,1) * t470 + mrSges(2,3) * t460 + mrSges(3,1) * t436 + mrSges(3,2) * t437 - t476 * t393 - t478 * t398 - pkin(2) * (-m(4) * t436 - t412) - pkin(5) * t485 - pkin(1) * t397;
t389 = mrSges(2,2) * t470 - mrSges(2,3) * t459 - qJ(2) * t397 + t473 * t391 - t471 * t392;
t1 = [-m(1) * g(1) + t487; -m(1) * g(2) + t495; -m(1) * g(3) + t488; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t495 + t474 * t389 - t472 * t390; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t487 + t472 * t389 + t474 * t390; -mrSges(1,1) * g(2) + mrSges(2,1) * t459 + mrSges(1,2) * g(1) - mrSges(2,2) * t460 - pkin(1) * t404 + qJ(2) * t486 + t471 * t391 + t473 * t392; t488; t404; t481; t504; t421;];
tauJB = t1;
