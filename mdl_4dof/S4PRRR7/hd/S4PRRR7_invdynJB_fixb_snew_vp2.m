% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+4)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4PRRR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR7_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR7_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR7_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR7_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR7_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:14
% EndTime: 2019-12-31 16:36:16
% DurationCPUTime: 1.64s
% Computational Cost: add. (17609->194), mult. (33057->258), div. (0->0), fcn. (21435->10), ass. (0->91)
t468 = sin(pkin(8));
t470 = cos(pkin(8));
t460 = t468 * g(1) - t470 * g(2);
t461 = -t470 * g(1) - t468 * g(2);
t467 = -g(3) + qJDD(1);
t469 = sin(pkin(4));
t471 = cos(pkin(4));
t474 = sin(qJ(2));
t477 = cos(qJ(2));
t431 = -t474 * t461 + (t460 * t471 + t467 * t469) * t477;
t496 = t471 * t474;
t497 = t469 * t474;
t432 = t460 * t496 + t477 * t461 + t467 * t497;
t479 = qJD(2) ^ 2;
t430 = -t479 * pkin(2) + qJDD(2) * pkin(6) + t432;
t444 = -t469 * t460 + t471 * t467;
t473 = sin(qJ(3));
t476 = cos(qJ(3));
t427 = t476 * t430 + t473 * t444;
t457 = (-pkin(3) * t476 - pkin(7) * t473) * qJD(2);
t478 = qJD(3) ^ 2;
t492 = t476 * qJD(2);
t424 = -t478 * pkin(3) + qJDD(3) * pkin(7) + t457 * t492 + t427;
t429 = -qJDD(2) * pkin(2) - t479 * pkin(6) - t431;
t491 = qJD(2) * qJD(3);
t489 = t476 * t491;
t458 = t473 * qJDD(2) + t489;
t490 = t473 * t491;
t459 = t476 * qJDD(2) - t490;
t425 = (-t458 - t489) * pkin(7) + (-t459 + t490) * pkin(3) + t429;
t472 = sin(qJ(4));
t475 = cos(qJ(4));
t421 = -t472 * t424 + t475 * t425;
t493 = qJD(2) * t473;
t454 = t475 * qJD(3) - t472 * t493;
t439 = t454 * qJD(4) + t472 * qJDD(3) + t475 * t458;
t455 = t472 * qJD(3) + t475 * t493;
t440 = -t454 * mrSges(5,1) + t455 * mrSges(5,2);
t465 = qJD(4) - t492;
t442 = -t465 * mrSges(5,2) + t454 * mrSges(5,3);
t451 = qJDD(4) - t459;
t418 = m(5) * t421 + t451 * mrSges(5,1) - t439 * mrSges(5,3) - t455 * t440 + t465 * t442;
t422 = t475 * t424 + t472 * t425;
t438 = -t455 * qJD(4) + t475 * qJDD(3) - t472 * t458;
t443 = t465 * mrSges(5,1) - t455 * mrSges(5,3);
t419 = m(5) * t422 - t451 * mrSges(5,2) + t438 * mrSges(5,3) + t454 * t440 - t465 * t443;
t412 = -t472 * t418 + t475 * t419;
t495 = t476 * t444;
t423 = -qJDD(3) * pkin(3) - t478 * pkin(7) - t495 + (qJD(2) * t457 + t430) * t473;
t433 = Ifges(5,5) * t455 + Ifges(5,6) * t454 + Ifges(5,3) * t465;
t435 = Ifges(5,1) * t455 + Ifges(5,4) * t454 + Ifges(5,5) * t465;
t413 = -mrSges(5,1) * t423 + mrSges(5,3) * t422 + Ifges(5,4) * t439 + Ifges(5,2) * t438 + Ifges(5,6) * t451 - t455 * t433 + t465 * t435;
t434 = Ifges(5,4) * t455 + Ifges(5,2) * t454 + Ifges(5,6) * t465;
t414 = mrSges(5,2) * t423 - mrSges(5,3) * t421 + Ifges(5,1) * t439 + Ifges(5,4) * t438 + Ifges(5,5) * t451 + t454 * t433 - t465 * t434;
t420 = -m(5) * t423 + t438 * mrSges(5,1) - t439 * mrSges(5,2) + t454 * t442 - t455 * t443;
t426 = -t473 * t430 + t495;
t447 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t473 + Ifges(4,2) * t476) * qJD(2);
t448 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t473 + Ifges(4,4) * t476) * qJD(2);
t499 = mrSges(4,1) * t426 - mrSges(4,2) * t427 + Ifges(4,5) * t458 + Ifges(4,6) * t459 + Ifges(4,3) * qJDD(3) + pkin(3) * t420 + pkin(7) * t412 + t475 * t413 + t472 * t414 + (t473 * t447 - t476 * t448) * qJD(2);
t411 = t475 * t418 + t472 * t419;
t462 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t493;
t463 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t492;
t482 = -m(4) * t429 + t459 * mrSges(4,1) - t458 * mrSges(4,2) - t462 * t493 + t463 * t492 - t411;
t407 = m(3) * t431 + qJDD(2) * mrSges(3,1) - t479 * mrSges(3,2) + t482;
t498 = t407 * t477;
t456 = (-mrSges(4,1) * t476 + mrSges(4,2) * t473) * qJD(2);
t410 = m(4) * t427 - qJDD(3) * mrSges(4,2) + t459 * mrSges(4,3) - qJD(3) * t462 + t456 * t492 + t412;
t416 = m(4) * t426 + qJDD(3) * mrSges(4,1) - t458 * mrSges(4,3) + qJD(3) * t463 - t456 * t493 + t420;
t487 = t476 * t410 - t473 * t416;
t401 = m(3) * t432 - t479 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t487;
t404 = t473 * t410 + t476 * t416;
t403 = m(3) * t444 + t404;
t391 = t401 * t496 - t469 * t403 + t471 * t498;
t389 = m(2) * t460 + t391;
t395 = t477 * t401 - t474 * t407;
t394 = m(2) * t461 + t395;
t494 = t470 * t389 + t468 * t394;
t390 = t401 * t497 + t471 * t403 + t469 * t498;
t488 = -t468 * t389 + t470 * t394;
t486 = m(2) * t467 + t390;
t446 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t473 + Ifges(4,6) * t476) * qJD(2);
t396 = mrSges(4,2) * t429 - mrSges(4,3) * t426 + Ifges(4,1) * t458 + Ifges(4,4) * t459 + Ifges(4,5) * qJDD(3) - pkin(7) * t411 - qJD(3) * t447 - t472 * t413 + t475 * t414 + t446 * t492;
t481 = mrSges(5,1) * t421 - mrSges(5,2) * t422 + Ifges(5,5) * t439 + Ifges(5,6) * t438 + Ifges(5,3) * t451 + t455 * t434 - t454 * t435;
t397 = -mrSges(4,1) * t429 + mrSges(4,3) * t427 + Ifges(4,4) * t458 + Ifges(4,2) * t459 + Ifges(4,6) * qJDD(3) - pkin(3) * t411 + qJD(3) * t448 - t446 * t493 - t481;
t386 = mrSges(3,2) * t444 - mrSges(3,3) * t431 + Ifges(3,5) * qJDD(2) - t479 * Ifges(3,6) - pkin(6) * t404 + t476 * t396 - t473 * t397;
t387 = -mrSges(3,1) * t444 + mrSges(3,3) * t432 + t479 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t404 - t499;
t483 = pkin(5) * t395 + t386 * t474 + t387 * t477;
t385 = mrSges(3,1) * t431 - mrSges(3,2) * t432 + Ifges(3,3) * qJDD(2) + pkin(2) * t482 + pkin(6) * t487 + t473 * t396 + t476 * t397;
t384 = mrSges(2,2) * t467 - mrSges(2,3) * t460 + t477 * t386 - t474 * t387 + (-t390 * t469 - t391 * t471) * pkin(5);
t383 = -mrSges(2,1) * t467 + mrSges(2,3) * t461 - pkin(1) * t390 - t469 * t385 + t483 * t471;
t1 = [-m(1) * g(1) + t488; -m(1) * g(2) + t494; -m(1) * g(3) + t486; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t494 - t468 * t383 + t470 * t384; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t488 + t470 * t383 + t468 * t384; -mrSges(1,1) * g(2) + mrSges(2,1) * t460 + mrSges(1,2) * g(1) - mrSges(2,2) * t461 + pkin(1) * t391 + t471 * t385 + t483 * t469; t486; t385; t499; t481;];
tauJB = t1;
