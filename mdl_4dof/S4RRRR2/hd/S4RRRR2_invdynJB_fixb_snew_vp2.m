% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RRRR2
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RRRR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR2_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR2_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR2_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR2_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR2_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:14
% EndTime: 2019-12-31 17:23:16
% DurationCPUTime: 1.43s
% Computational Cost: add. (19210->200), mult. (24674->256), div. (0->0), fcn. (13904->8), ass. (0->86)
t469 = qJDD(1) + qJDD(2);
t475 = sin(qJ(3));
t479 = cos(qJ(3));
t471 = qJD(1) + qJD(2);
t495 = qJD(3) * t471;
t451 = t475 * t469 + t479 * t495;
t477 = sin(qJ(1));
t481 = cos(qJ(1));
t462 = t477 * g(1) - t481 * g(2);
t456 = qJDD(1) * pkin(1) + t462;
t463 = -t481 * g(1) - t477 * g(2);
t482 = qJD(1) ^ 2;
t457 = -t482 * pkin(1) + t463;
t476 = sin(qJ(2));
t480 = cos(qJ(2));
t440 = t476 * t456 + t480 * t457;
t467 = t471 ^ 2;
t437 = -t467 * pkin(2) + t469 * pkin(6) + t440;
t497 = t475 * t437;
t500 = pkin(3) * t467;
t419 = qJDD(3) * pkin(3) - t451 * pkin(7) - t497 + (pkin(7) * t495 + t475 * t500 - g(3)) * t479;
t429 = -t475 * g(3) + t479 * t437;
t452 = t479 * t469 - t475 * t495;
t499 = t471 * t475;
t460 = qJD(3) * pkin(3) - pkin(7) * t499;
t473 = t479 ^ 2;
t420 = t452 * pkin(7) - qJD(3) * t460 - t473 * t500 + t429;
t474 = sin(qJ(4));
t478 = cos(qJ(4));
t417 = t478 * t419 - t474 * t420;
t446 = (-t474 * t475 + t478 * t479) * t471;
t427 = t446 * qJD(4) + t478 * t451 + t474 * t452;
t447 = (t474 * t479 + t475 * t478) * t471;
t435 = -t446 * mrSges(5,1) + t447 * mrSges(5,2);
t470 = qJD(3) + qJD(4);
t441 = -t470 * mrSges(5,2) + t446 * mrSges(5,3);
t468 = qJDD(3) + qJDD(4);
t414 = m(5) * t417 + t468 * mrSges(5,1) - t427 * mrSges(5,3) - t447 * t435 + t470 * t441;
t418 = t474 * t419 + t478 * t420;
t426 = -t447 * qJD(4) - t474 * t451 + t478 * t452;
t442 = t470 * mrSges(5,1) - t447 * mrSges(5,3);
t415 = m(5) * t418 - t468 * mrSges(5,2) + t426 * mrSges(5,3) + t446 * t435 - t470 * t442;
t405 = t478 * t414 + t474 * t415;
t428 = -t479 * g(3) - t497;
t444 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t475 + Ifges(4,2) * t479) * t471;
t445 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t475 + Ifges(4,4) * t479) * t471;
t431 = Ifges(5,4) * t447 + Ifges(5,2) * t446 + Ifges(5,6) * t470;
t432 = Ifges(5,1) * t447 + Ifges(5,4) * t446 + Ifges(5,5) * t470;
t486 = -mrSges(5,1) * t417 + mrSges(5,2) * t418 - Ifges(5,5) * t427 - Ifges(5,6) * t426 - Ifges(5,3) * t468 - t447 * t431 + t446 * t432;
t501 = mrSges(4,1) * t428 - mrSges(4,2) * t429 + Ifges(4,5) * t451 + Ifges(4,6) * t452 + Ifges(4,3) * qJDD(3) + pkin(3) * t405 + (t475 * t444 - t479 * t445) * t471 - t486;
t498 = t471 * t479;
t450 = (-mrSges(4,1) * t479 + mrSges(4,2) * t475) * t471;
t459 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t498;
t403 = m(4) * t428 + qJDD(3) * mrSges(4,1) - t451 * mrSges(4,3) + qJD(3) * t459 - t450 * t499 + t405;
t458 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t499;
t491 = -t474 * t414 + t478 * t415;
t404 = m(4) * t429 - qJDD(3) * mrSges(4,2) + t452 * mrSges(4,3) - qJD(3) * t458 + t450 * t498 + t491;
t492 = -t475 * t403 + t479 * t404;
t397 = m(3) * t440 - t467 * mrSges(3,1) - t469 * mrSges(3,2) + t492;
t439 = t480 * t456 - t476 * t457;
t489 = -t469 * pkin(2) - t439;
t436 = -t467 * pkin(6) + t489;
t421 = t460 * t499 - t452 * pkin(3) + (-pkin(7) * t473 - pkin(6)) * t467 + t489;
t487 = m(5) * t421 - t426 * mrSges(5,1) + t427 * mrSges(5,2) - t446 * t441 + t447 * t442;
t484 = -m(4) * t436 + t452 * mrSges(4,1) - t451 * mrSges(4,2) - t458 * t499 + t459 * t498 - t487;
t409 = m(3) * t439 + t469 * mrSges(3,1) - t467 * mrSges(3,2) + t484;
t392 = t476 * t397 + t480 * t409;
t389 = m(2) * t462 + qJDD(1) * mrSges(2,1) - t482 * mrSges(2,2) + t392;
t493 = t480 * t397 - t476 * t409;
t390 = m(2) * t463 - t482 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t493;
t496 = t481 * t389 + t477 * t390;
t399 = t479 * t403 + t475 * t404;
t494 = -t477 * t389 + t481 * t390;
t430 = Ifges(5,5) * t447 + Ifges(5,6) * t446 + Ifges(5,3) * t470;
t406 = -mrSges(5,1) * t421 + mrSges(5,3) * t418 + Ifges(5,4) * t427 + Ifges(5,2) * t426 + Ifges(5,6) * t468 - t447 * t430 + t470 * t432;
t407 = mrSges(5,2) * t421 - mrSges(5,3) * t417 + Ifges(5,1) * t427 + Ifges(5,4) * t426 + Ifges(5,5) * t468 + t446 * t430 - t470 * t431;
t443 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t475 + Ifges(4,6) * t479) * t471;
t385 = -mrSges(4,1) * t436 + mrSges(4,3) * t429 + Ifges(4,4) * t451 + Ifges(4,2) * t452 + Ifges(4,6) * qJDD(3) - pkin(3) * t487 + pkin(7) * t491 + qJD(3) * t445 + t478 * t406 + t474 * t407 - t443 * t499;
t394 = mrSges(4,2) * t436 - mrSges(4,3) * t428 + Ifges(4,1) * t451 + Ifges(4,4) * t452 + Ifges(4,5) * qJDD(3) - pkin(7) * t405 - qJD(3) * t444 - t474 * t406 + t478 * t407 + t443 * t498;
t488 = mrSges(3,1) * t439 - mrSges(3,2) * t440 + Ifges(3,3) * t469 + pkin(2) * t484 + pkin(6) * t492 + t479 * t385 + t475 * t394;
t485 = mrSges(2,1) * t462 - mrSges(2,2) * t463 + Ifges(2,3) * qJDD(1) + pkin(1) * t392 + t488;
t383 = mrSges(3,1) * g(3) + mrSges(3,3) * t440 + t467 * Ifges(3,5) + Ifges(3,6) * t469 - pkin(2) * t399 - t501;
t382 = -mrSges(3,2) * g(3) - mrSges(3,3) * t439 + Ifges(3,5) * t469 - t467 * Ifges(3,6) - pkin(6) * t399 - t475 * t385 + t479 * t394;
t381 = -mrSges(2,2) * g(3) - mrSges(2,3) * t462 + Ifges(2,5) * qJDD(1) - t482 * Ifges(2,6) - pkin(5) * t392 + t480 * t382 - t476 * t383;
t380 = Ifges(2,6) * qJDD(1) + t482 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t463 + t476 * t382 + t480 * t383 - pkin(1) * (-m(3) * g(3) + t399) + pkin(5) * t493;
t1 = [-m(1) * g(1) + t494; -m(1) * g(2) + t496; (-m(1) - m(2) - m(3)) * g(3) + t399; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t496 - t477 * t380 + t481 * t381; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t494 + t481 * t380 + t477 * t381; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t485; t485; t488; t501; -t486;];
tauJB = t1;
