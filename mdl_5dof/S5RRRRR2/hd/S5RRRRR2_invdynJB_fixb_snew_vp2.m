% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRRR2
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-12-05 18:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRRR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_invdynJB_fixb_snew_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:53:11
% EndTime: 2019-12-05 18:53:14
% DurationCPUTime: 1.98s
% Computational Cost: add. (24804->238), mult. (33632->306), div. (0->0), fcn. (23058->10), ass. (0->95)
t479 = qJD(1) + qJD(2);
t482 = sin(qJ(4));
t483 = sin(qJ(3));
t487 = cos(qJ(4));
t488 = cos(qJ(3));
t460 = (t482 * t483 - t487 * t488) * t479;
t485 = sin(qJ(1));
t490 = cos(qJ(1));
t471 = -t490 * g(1) - t485 * g(2);
t491 = qJD(1) ^ 2;
t466 = -t491 * pkin(1) + t471;
t484 = sin(qJ(2));
t489 = cos(qJ(2));
t470 = t485 * g(1) - t490 * g(2);
t499 = qJDD(1) * pkin(1) + t470;
t451 = t489 * t466 + t484 * t499;
t447 = -t483 * g(3) + t488 * t451;
t475 = t479 ^ 2;
t441 = (-t475 * t488 ^ 2 - qJD(3) ^ 2) * pkin(2) + t447;
t446 = -t488 * g(3) - t483 * t451;
t497 = (t475 * t483 * t488 + qJDD(3)) * pkin(2) + t446;
t424 = t487 * t441 + t482 * t497;
t450 = t484 * t466 - t489 * t499;
t477 = qJDD(1) + qJDD(2);
t503 = qJD(3) * t479;
t502 = t483 * t503;
t464 = t488 * t477 - t502;
t437 = (-t464 + t502) * pkin(2) + t450;
t481 = sin(qJ(5));
t486 = cos(qJ(5));
t420 = -t481 * t424 + t486 * t437;
t463 = t483 * t477 + t488 * t503;
t436 = -t460 * qJD(4) + t487 * t463 + t482 * t464;
t461 = (t482 * t488 + t483 * t487) * t479;
t478 = qJD(3) + qJD(4);
t452 = -t481 * t461 + t486 * t478;
t476 = qJDD(3) + qJDD(4);
t426 = t452 * qJD(5) + t486 * t436 + t481 * t476;
t453 = t486 * t461 + t481 * t478;
t431 = -t452 * mrSges(6,1) + t453 * mrSges(6,2);
t435 = -t461 * qJD(4) - t482 * t463 + t487 * t464;
t434 = qJDD(5) - t435;
t456 = qJD(5) + t460;
t439 = -t456 * mrSges(6,2) + t452 * mrSges(6,3);
t418 = m(6) * t420 + t434 * mrSges(6,1) - t426 * mrSges(6,3) - t453 * t431 + t456 * t439;
t421 = t486 * t424 + t481 * t437;
t425 = -t453 * qJD(5) - t481 * t436 + t486 * t476;
t440 = t456 * mrSges(6,1) - t453 * mrSges(6,3);
t419 = m(6) * t421 - t434 * mrSges(6,2) + t425 * mrSges(6,3) + t452 * t431 - t456 * t440;
t445 = t460 * mrSges(5,1) + t461 * mrSges(5,2);
t455 = t478 * mrSges(5,1) - t461 * mrSges(5,3);
t410 = m(5) * t424 - t476 * mrSges(5,2) + t435 * mrSges(5,3) - t481 * t418 + t486 * t419 - t460 * t445 - t478 * t455;
t423 = t482 * t441 - t487 * t497;
t454 = -t478 * mrSges(5,2) - t460 * mrSges(5,3);
t417 = t476 * mrSges(5,1) + t425 * mrSges(6,1) - t426 * mrSges(6,2) - t436 * mrSges(5,3) + t452 * t439 - t453 * t440 - t461 * t445 + t478 * t454 + (-m(5) - m(6)) * t423;
t405 = t482 * t410 + t487 * t417;
t458 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t483 + Ifges(4,2) * t488) * t479;
t459 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t483 + Ifges(4,4) * t488) * t479;
t427 = Ifges(6,5) * t453 + Ifges(6,6) * t452 + Ifges(6,3) * t456;
t429 = Ifges(6,1) * t453 + Ifges(6,4) * t452 + Ifges(6,5) * t456;
t414 = -mrSges(6,1) * t423 + mrSges(6,3) * t421 + Ifges(6,4) * t426 + Ifges(6,2) * t425 + Ifges(6,6) * t434 - t453 * t427 + t456 * t429;
t428 = Ifges(6,4) * t453 + Ifges(6,2) * t452 + Ifges(6,6) * t456;
t415 = mrSges(6,2) * t423 - mrSges(6,3) * t420 + Ifges(6,1) * t426 + Ifges(6,4) * t425 + Ifges(6,5) * t434 + t452 * t427 - t456 * t428;
t443 = Ifges(5,4) * t461 - Ifges(5,2) * t460 + Ifges(5,6) * t478;
t444 = Ifges(5,1) * t461 - Ifges(5,4) * t460 + Ifges(5,5) * t478;
t495 = mrSges(5,1) * t423 + mrSges(5,2) * t424 - Ifges(5,5) * t436 - Ifges(5,6) * t435 - Ifges(5,3) * t476 - t486 * t414 - t481 * t415 - t461 * t443 - t460 * t444;
t508 = mrSges(4,1) * t446 - mrSges(4,2) * t447 + Ifges(4,5) * t463 + Ifges(4,6) * t464 + Ifges(4,3) * qJDD(3) + pkin(2) * t405 + (t483 * t458 - t488 * t459) * t479 - t495;
t507 = t479 * t483;
t506 = t479 * t488;
t462 = (-mrSges(4,1) * t488 + mrSges(4,2) * t483) * t479;
t468 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t506;
t403 = m(4) * t446 + qJDD(3) * mrSges(4,1) - t463 * mrSges(4,3) + qJD(3) * t468 - t462 * t507 + t405;
t467 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t507;
t404 = m(4) * t447 - qJDD(3) * mrSges(4,2) + t464 * mrSges(4,3) - qJD(3) * t467 + t487 * t410 - t482 * t417 + t462 * t506;
t397 = m(3) * t451 - t475 * mrSges(3,1) - t477 * mrSges(3,2) - t483 * t403 + t488 * t404;
t493 = m(5) * t437 - t435 * mrSges(5,1) + t436 * mrSges(5,2) + t486 * t418 + t481 * t419 + t460 * t454 + t461 * t455;
t408 = t477 * mrSges(3,1) + t464 * mrSges(4,1) - t475 * mrSges(3,2) - t463 * mrSges(4,2) + (-t467 * t483 + t468 * t488) * t479 + (-m(3) - m(4)) * t450 - t493;
t505 = t484 * t397 + t489 * t408;
t504 = t488 * t403 + t483 * t404;
t442 = Ifges(5,5) * t461 - Ifges(5,6) * t460 + Ifges(5,3) * t478;
t406 = mrSges(5,2) * t437 + mrSges(5,3) * t423 + Ifges(5,1) * t436 + Ifges(5,4) * t435 + Ifges(5,5) * t476 - t481 * t414 + t486 * t415 - t460 * t442 - t478 * t443;
t494 = mrSges(6,1) * t420 - mrSges(6,2) * t421 + Ifges(6,5) * t426 + Ifges(6,6) * t425 + Ifges(6,3) * t434 + t453 * t428 - t452 * t429;
t411 = -mrSges(5,1) * t437 + mrSges(5,3) * t424 + Ifges(5,4) * t436 + Ifges(5,2) * t435 + Ifges(5,6) * t476 - t461 * t442 + t478 * t444 - t494;
t457 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t483 + Ifges(4,6) * t488) * t479;
t395 = -mrSges(4,1) * t450 + mrSges(4,3) * t447 + Ifges(4,4) * t463 + Ifges(4,2) * t464 + Ifges(4,6) * qJDD(3) - pkin(2) * t493 + qJD(3) * t459 + t482 * t406 + t487 * t411 - t457 * t507;
t400 = mrSges(4,2) * t450 - mrSges(4,3) * t446 + Ifges(4,1) * t463 + Ifges(4,4) * t464 + Ifges(4,5) * qJDD(3) - qJD(3) * t458 + t487 * t406 - t482 * t411 + t457 * t506;
t498 = -mrSges(3,1) * t450 - mrSges(3,2) * t451 + Ifges(3,3) * t477 + t488 * t395 + t483 * t400;
t496 = mrSges(2,1) * t470 - mrSges(2,2) * t471 + Ifges(2,3) * qJDD(1) + pkin(1) * t505 + t498;
t398 = mrSges(3,1) * g(3) + mrSges(3,3) * t451 + t475 * Ifges(3,5) + Ifges(3,6) * t477 - t508;
t392 = m(2) * t471 - t491 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t489 * t397 - t484 * t408;
t391 = m(2) * t470 + qJDD(1) * mrSges(2,1) - t491 * mrSges(2,2) + t505;
t390 = -mrSges(3,2) * g(3) + mrSges(3,3) * t450 + Ifges(3,5) * t477 - t475 * Ifges(3,6) - t483 * t395 + t488 * t400;
t389 = -mrSges(2,2) * g(3) - mrSges(2,3) * t470 + Ifges(2,5) * qJDD(1) - t491 * Ifges(2,6) + t489 * t390 - t484 * t398;
t388 = Ifges(2,6) * qJDD(1) + t491 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t471 + t484 * t390 + t489 * t398 - pkin(1) * (-m(3) * g(3) + t504);
t1 = [-m(1) * g(1) - t485 * t391 + t490 * t392; -m(1) * g(2) + t490 * t391 + t485 * t392; (-m(1) - m(2) - m(3)) * g(3) + t504; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - t485 * t388 + t490 * t389; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t490 * t388 + t485 * t389; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t496; t496; t498; t508; -t495; t494;];
tauJB = t1;
