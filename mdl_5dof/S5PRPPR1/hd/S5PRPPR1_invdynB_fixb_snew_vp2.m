% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRPPR1
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRPPR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR1_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR1_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR1_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR1_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR1_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:22:01
% EndTime: 2019-12-05 15:22:04
% DurationCPUTime: 3.34s
% Computational Cost: add. (27901->225), mult. (67654->319), div. (0->0), fcn. (44132->10), ass. (0->111)
t489 = sin(pkin(7));
t492 = cos(pkin(7));
t473 = t489 * g(1) - t492 * g(2);
t474 = -t492 * g(1) - t489 * g(2);
t494 = sin(qJ(2));
t496 = cos(qJ(2));
t454 = t494 * t473 + t496 * t474;
t497 = qJD(2) ^ 2;
t536 = -t497 * pkin(2) + qJDD(2) * qJ(3) + (2 * qJD(2) * qJD(3)) + t454;
t453 = t496 * t473 - t494 * t474;
t503 = -t497 * qJ(3) + qJDD(3) - t453;
t488 = sin(pkin(8));
t491 = cos(pkin(8));
t509 = -pkin(3) * t491 - qJ(4) * t488;
t527 = t488 * qJD(2);
t535 = (-pkin(2) + t509) * qJDD(2) + t503 - 0.2e1 * qJD(4) * t527;
t486 = -g(3) + qJDD(1);
t440 = t491 * t486 - t536 * t488;
t534 = mrSges(4,2) * t488;
t533 = Ifges(4,6) * t491;
t485 = t488 ^ 2;
t532 = t485 * t497;
t487 = sin(pkin(9));
t531 = t487 * t488;
t490 = cos(pkin(9));
t530 = t488 * t490;
t441 = t488 * t486 + t536 * t491;
t469 = (-mrSges(4,1) * t491 + t534) * qJD(2);
t468 = t509 * qJD(2);
t526 = t491 * qJD(2);
t433 = t468 * t526 + t441;
t508 = -pkin(4) * t491 - pkin(6) * t530;
t528 = t535 * t490;
t426 = t508 * qJDD(2) + (-t433 + (-pkin(4) * t485 * t490 + pkin(6) * t488 * t491) * t497) * t487 + t528;
t429 = t490 * t433 + t535 * t487;
t464 = t508 * qJD(2);
t522 = t487 ^ 2 * t532;
t524 = qJDD(2) * t488;
t427 = -t487 * pkin(6) * t524 - pkin(4) * t522 + t464 * t526 + t429;
t493 = sin(qJ(5));
t495 = cos(qJ(5));
t424 = t495 * t426 - t493 * t427;
t505 = (-t487 * t495 - t490 * t493) * t488;
t458 = qJD(2) * t505;
t504 = (-t487 * t493 + t490 * t495) * t488;
t459 = qJD(2) * t504;
t442 = -t458 * mrSges(6,1) + t459 * mrSges(6,2);
t446 = t458 * qJD(5) + qJDD(2) * t504;
t476 = qJD(5) - t526;
t451 = -t476 * mrSges(6,2) + t458 * mrSges(6,3);
t523 = t491 * qJDD(2);
t475 = qJDD(5) - t523;
t422 = m(6) * t424 + t475 * mrSges(6,1) - t446 * mrSges(6,3) - t459 * t442 + t476 * t451;
t425 = t493 * t426 + t495 * t427;
t445 = -t459 * qJD(5) + qJDD(2) * t505;
t452 = t476 * mrSges(6,1) - t459 * mrSges(6,3);
t423 = m(6) * t425 - t475 * mrSges(6,2) + t445 * mrSges(6,3) + t458 * t442 - t476 * t452;
t414 = t495 * t422 + t493 * t423;
t428 = -t487 * t433 + t528;
t512 = mrSges(5,1) * t487 + mrSges(5,2) * t490;
t460 = t512 * t527;
t506 = mrSges(5,2) * t491 - mrSges(5,3) * t531;
t462 = t506 * qJD(2);
t507 = -mrSges(5,1) * t491 - mrSges(5,3) * t530;
t412 = m(5) * t428 + t507 * qJDD(2) + (-t460 * t530 - t462 * t491) * qJD(2) + t414;
t463 = t507 * qJD(2);
t515 = -t493 * t422 + t495 * t423;
t413 = m(5) * t429 + t506 * qJDD(2) + (-t460 * t531 + t463 * t491) * qJD(2) + t515;
t516 = -t487 * t412 + t490 * t413;
t409 = m(4) * t441 + (qJDD(2) * mrSges(4,3) + qJD(2) * t469) * t491 + t516;
t432 = t468 * t527 + qJDD(4) - t440;
t430 = -pkin(6) * t522 + (pkin(4) * qJDD(2) * t487 + qJD(2) * t464 * t490) * t488 + t432;
t502 = m(6) * t430 - t445 * mrSges(6,1) + t446 * mrSges(6,2) - t458 * t451 + t459 * t452;
t498 = -m(5) * t432 - t502;
t418 = m(4) * t440 + ((-mrSges(4,3) - t512) * qJDD(2) + (-t462 * t487 - t463 * t490 - t469) * qJD(2)) * t488 + t498;
t517 = t491 * t409 - t488 * t418;
t401 = m(3) * t454 - t497 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t517;
t410 = t490 * t412 + t487 * t413;
t449 = -qJDD(2) * pkin(2) + t503;
t499 = -m(4) * t449 + mrSges(4,1) * t523 - t410 + (t491 ^ 2 * t497 + t532) * mrSges(4,3);
t406 = m(3) * t453 - t497 * mrSges(3,2) + (mrSges(3,1) - t534) * qJDD(2) + t499;
t397 = t494 * t401 + t496 * t406;
t395 = m(2) * t473 + t397;
t518 = t496 * t401 - t494 * t406;
t396 = m(2) * t474 + t518;
t529 = t492 * t395 + t489 * t396;
t402 = t488 * t409 + t491 * t418;
t521 = m(3) * t486 + t402;
t519 = -t489 * t395 + t492 * t396;
t511 = Ifges(4,1) * t488 + Ifges(4,4) * t491;
t510 = Ifges(5,5) * t490 - Ifges(5,6) * t487;
t501 = -Ifges(5,5) * t491 + (Ifges(5,1) * t490 - Ifges(5,4) * t487) * t488;
t500 = -Ifges(5,6) * t491 + (Ifges(5,4) * t490 - Ifges(5,2) * t487) * t488;
t470 = (Ifges(4,5) * t488 + t533) * qJD(2);
t457 = t501 * qJD(2);
t456 = t500 * qJD(2);
t455 = (-Ifges(5,3) * t491 + t510 * t488) * qJD(2);
t436 = Ifges(6,1) * t459 + Ifges(6,4) * t458 + Ifges(6,5) * t476;
t435 = Ifges(6,4) * t459 + Ifges(6,2) * t458 + Ifges(6,6) * t476;
t434 = Ifges(6,5) * t459 + Ifges(6,6) * t458 + Ifges(6,3) * t476;
t416 = mrSges(6,2) * t430 - mrSges(6,3) * t424 + Ifges(6,1) * t446 + Ifges(6,4) * t445 + Ifges(6,5) * t475 + t458 * t434 - t476 * t435;
t415 = -mrSges(6,1) * t430 + mrSges(6,3) * t425 + Ifges(6,4) * t446 + Ifges(6,2) * t445 + Ifges(6,6) * t475 - t459 * t434 + t476 * t436;
t404 = mrSges(5,2) * t432 - mrSges(5,3) * t428 - pkin(6) * t414 - t493 * t415 + t495 * t416 + (-t455 * t531 + t456 * t491) * qJD(2) + t501 * qJDD(2);
t403 = -mrSges(5,1) * t432 + mrSges(5,3) * t429 + t493 * t416 + t495 * t415 - pkin(4) * t502 + pkin(6) * t515 + (-t455 * t530 - t491 * t457) * qJD(2) + t500 * qJDD(2);
t398 = -mrSges(4,1) * t449 - mrSges(5,1) * t428 - mrSges(6,1) * t424 + mrSges(5,2) * t429 + mrSges(6,2) * t425 + mrSges(4,3) * t441 - Ifges(6,5) * t446 - Ifges(6,6) * t445 - Ifges(6,3) * t475 - pkin(3) * t410 - pkin(4) * t414 - t459 * t435 + t458 * t436 + (Ifges(4,2) + Ifges(5,3)) * t523 + ((Ifges(4,4) - t510) * qJDD(2) + (-t456 * t490 - t457 * t487 - t470) * qJD(2)) * t488;
t391 = mrSges(4,2) * t449 - mrSges(4,3) * t440 - qJ(4) * t410 + t511 * qJDD(2) - t487 * t403 + t490 * t404 + t470 * t526;
t390 = t497 * Ifges(3,5) - mrSges(3,1) * t486 + mrSges(3,3) * t454 - mrSges(4,1) * t440 + mrSges(4,2) * t441 - t487 * t404 - t490 * t403 - pkin(3) * t498 - qJ(4) * t516 - pkin(2) * t402 + (-t533 + Ifges(3,6) + (pkin(3) * t512 - Ifges(4,5)) * t488) * qJDD(2) + (-pkin(3) * (-t462 * t531 - t463 * t530) + (-t488 * (Ifges(4,4) * t488 + Ifges(4,2) * t491) + t491 * t511) * qJD(2)) * qJD(2);
t389 = mrSges(3,2) * t486 - mrSges(3,3) * t453 + Ifges(3,5) * qJDD(2) - t497 * Ifges(3,6) - qJ(3) * t402 + t491 * t391 - t488 * t398;
t388 = mrSges(2,2) * t486 - mrSges(2,3) * t473 - pkin(5) * t397 + t496 * t389 - t494 * t390;
t387 = -mrSges(2,1) * t486 + mrSges(2,3) * t474 - pkin(1) * t521 + pkin(5) * t518 + t494 * t389 + t496 * t390;
t1 = [-m(1) * g(1) + t519; -m(1) * g(2) + t529; -m(1) * g(3) + m(2) * t486 + t521; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t529 - t489 * t387 + t492 * t388; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t519 + t492 * t387 + t489 * t388; pkin(1) * t397 + mrSges(2,1) * t473 - mrSges(2,2) * t474 + t491 * t398 + pkin(2) * (-mrSges(4,2) * t524 + t499) + qJ(3) * t517 + t488 * t391 + mrSges(3,1) * t453 - mrSges(3,2) * t454 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(3,3) * qJDD(2);];
tauB = t1;
