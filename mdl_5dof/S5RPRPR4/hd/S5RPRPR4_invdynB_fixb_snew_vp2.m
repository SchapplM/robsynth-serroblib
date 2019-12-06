% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRPR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRPR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR4_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR4_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR4_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:53:32
% EndTime: 2019-12-05 17:53:38
% DurationCPUTime: 4.54s
% Computational Cost: add. (53872->269), mult. (115951->344), div. (0->0), fcn. (73622->10), ass. (0->106)
t516 = sin(qJ(1));
t519 = cos(qJ(1));
t499 = t519 * g(2) + t516 * g(3);
t490 = qJDD(1) * pkin(1) + t499;
t498 = t516 * g(2) - t519 * g(3);
t520 = qJD(1) ^ 2;
t492 = -t520 * pkin(1) + t498;
t511 = sin(pkin(8));
t513 = cos(pkin(8));
t470 = t511 * t490 + t513 * t492;
t465 = -t520 * pkin(2) + qJDD(1) * pkin(6) + t470;
t509 = -g(1) + qJDD(2);
t515 = sin(qJ(3));
t518 = cos(qJ(3));
t454 = -t515 * t465 + t518 * t509;
t533 = qJD(1) * qJD(3);
t531 = t518 * t533;
t493 = t515 * qJDD(1) + t531;
t451 = (-t493 + t531) * qJ(4) + (t515 * t518 * t520 + qJDD(3)) * pkin(3) + t454;
t455 = t518 * t465 + t515 * t509;
t494 = t518 * qJDD(1) - t515 * t533;
t535 = qJD(1) * t515;
t495 = qJD(3) * pkin(3) - qJ(4) * t535;
t508 = t518 ^ 2;
t452 = -t508 * t520 * pkin(3) + t494 * qJ(4) - qJD(3) * t495 + t455;
t510 = sin(pkin(9));
t512 = cos(pkin(9));
t480 = (t510 * t518 + t512 * t515) * qJD(1);
t434 = -0.2e1 * qJD(4) * t480 + t512 * t451 - t510 * t452;
t472 = t512 * t493 + t510 * t494;
t479 = (-t510 * t515 + t512 * t518) * qJD(1);
t432 = (qJD(3) * t479 - t472) * pkin(7) + (t479 * t480 + qJDD(3)) * pkin(4) + t434;
t435 = 0.2e1 * qJD(4) * t479 + t510 * t451 + t512 * t452;
t471 = -t510 * t493 + t512 * t494;
t475 = qJD(3) * pkin(4) - t480 * pkin(7);
t478 = t479 ^ 2;
t433 = -t478 * pkin(4) + t471 * pkin(7) - qJD(3) * t475 + t435;
t514 = sin(qJ(5));
t517 = cos(qJ(5));
t430 = t517 * t432 - t514 * t433;
t462 = t517 * t479 - t514 * t480;
t441 = t462 * qJD(5) + t514 * t471 + t517 * t472;
t463 = t514 * t479 + t517 * t480;
t450 = -t462 * mrSges(6,1) + t463 * mrSges(6,2);
t507 = qJD(3) + qJD(5);
t456 = -t507 * mrSges(6,2) + t462 * mrSges(6,3);
t506 = qJDD(3) + qJDD(5);
t428 = m(6) * t430 + t506 * mrSges(6,1) - t441 * mrSges(6,3) - t463 * t450 + t507 * t456;
t431 = t514 * t432 + t517 * t433;
t440 = -t463 * qJD(5) + t517 * t471 - t514 * t472;
t457 = t507 * mrSges(6,1) - t463 * mrSges(6,3);
t429 = m(6) * t431 - t506 * mrSges(6,2) + t440 * mrSges(6,3) + t462 * t450 - t507 * t457;
t420 = t517 * t428 + t514 * t429;
t467 = -t479 * mrSges(5,1) + t480 * mrSges(5,2);
t473 = -qJD(3) * mrSges(5,2) + t479 * mrSges(5,3);
t418 = m(5) * t434 + qJDD(3) * mrSges(5,1) - t472 * mrSges(5,3) + qJD(3) * t473 - t480 * t467 + t420;
t474 = qJD(3) * mrSges(5,1) - t480 * mrSges(5,3);
t526 = -t514 * t428 + t517 * t429;
t419 = m(5) * t435 - qJDD(3) * mrSges(5,2) + t471 * mrSges(5,3) - qJD(3) * t474 + t479 * t467 + t526;
t414 = t512 * t418 + t510 * t419;
t491 = (-mrSges(4,1) * t518 + mrSges(4,2) * t515) * qJD(1);
t534 = qJD(1) * t518;
t497 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t534;
t412 = m(4) * t454 + qJDD(3) * mrSges(4,1) - t493 * mrSges(4,3) + qJD(3) * t497 - t491 * t535 + t414;
t496 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t535;
t527 = -t510 * t418 + t512 * t419;
t413 = m(4) * t455 - qJDD(3) * mrSges(4,2) + t494 * mrSges(4,3) - qJD(3) * t496 + t491 * t534 + t527;
t528 = -t515 * t412 + t518 * t413;
t405 = m(3) * t470 - t520 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t528;
t469 = t513 * t490 - t511 * t492;
t523 = -qJDD(1) * pkin(2) - t469;
t464 = -t520 * pkin(6) + t523;
t453 = -t494 * pkin(3) + qJDD(4) + t495 * t535 + (-qJ(4) * t508 - pkin(6)) * t520 + t523;
t437 = -t471 * pkin(4) - t478 * pkin(7) + t480 * t475 + t453;
t525 = m(6) * t437 - t440 * mrSges(6,1) + t441 * mrSges(6,2) - t462 * t456 + t463 * t457;
t522 = m(5) * t453 - t471 * mrSges(5,1) + t472 * mrSges(5,2) - t479 * t473 + t480 * t474 + t525;
t521 = -m(4) * t464 + t494 * mrSges(4,1) - t493 * mrSges(4,2) - t496 * t535 + t497 * t534 - t522;
t424 = m(3) * t469 + qJDD(1) * mrSges(3,1) - t520 * mrSges(3,2) + t521;
t402 = t511 * t405 + t513 * t424;
t406 = t518 * t412 + t515 * t413;
t532 = m(3) * t509 + t406;
t529 = t513 * t405 - t511 * t424;
t400 = m(2) * t498 - t520 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t529;
t401 = m(2) * t499 + qJDD(1) * mrSges(2,1) - t520 * mrSges(2,2) + t402;
t530 = t519 * t400 - t516 * t401;
t524 = -t516 * t400 - t519 * t401;
t486 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t515 + Ifges(4,4) * t518) * qJD(1);
t485 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t515 + Ifges(4,2) * t518) * qJD(1);
t484 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t515 + Ifges(4,6) * t518) * qJD(1);
t461 = Ifges(5,1) * t480 + Ifges(5,4) * t479 + Ifges(5,5) * qJD(3);
t460 = Ifges(5,4) * t480 + Ifges(5,2) * t479 + Ifges(5,6) * qJD(3);
t459 = Ifges(5,5) * t480 + Ifges(5,6) * t479 + Ifges(5,3) * qJD(3);
t444 = Ifges(6,1) * t463 + Ifges(6,4) * t462 + Ifges(6,5) * t507;
t443 = Ifges(6,4) * t463 + Ifges(6,2) * t462 + Ifges(6,6) * t507;
t442 = Ifges(6,5) * t463 + Ifges(6,6) * t462 + Ifges(6,3) * t507;
t422 = mrSges(6,2) * t437 - mrSges(6,3) * t430 + Ifges(6,1) * t441 + Ifges(6,4) * t440 + Ifges(6,5) * t506 + t462 * t442 - t507 * t443;
t421 = -mrSges(6,1) * t437 + mrSges(6,3) * t431 + Ifges(6,4) * t441 + Ifges(6,2) * t440 + Ifges(6,6) * t506 - t463 * t442 + t507 * t444;
t408 = mrSges(5,2) * t453 - mrSges(5,3) * t434 + Ifges(5,1) * t472 + Ifges(5,4) * t471 + Ifges(5,5) * qJDD(3) - pkin(7) * t420 - qJD(3) * t460 - t514 * t421 + t517 * t422 + t479 * t459;
t407 = -mrSges(5,1) * t453 + mrSges(5,3) * t435 + Ifges(5,4) * t472 + Ifges(5,2) * t471 + Ifges(5,6) * qJDD(3) - pkin(4) * t525 + pkin(7) * t526 + qJD(3) * t461 + t517 * t421 + t514 * t422 - t480 * t459;
t398 = mrSges(4,2) * t464 - mrSges(4,3) * t454 + Ifges(4,1) * t493 + Ifges(4,4) * t494 + Ifges(4,5) * qJDD(3) - qJ(4) * t414 - qJD(3) * t485 - t510 * t407 + t512 * t408 + t484 * t534;
t397 = -mrSges(4,1) * t464 + mrSges(4,3) * t455 + Ifges(4,4) * t493 + Ifges(4,2) * t494 + Ifges(4,6) * qJDD(3) - pkin(3) * t522 + qJ(4) * t527 + qJD(3) * t486 + t512 * t407 + t510 * t408 - t484 * t535;
t396 = mrSges(6,2) * t431 - mrSges(5,1) * t434 + mrSges(5,2) * t435 - mrSges(6,1) * t430 - pkin(4) * t420 + (-Ifges(4,3) - Ifges(5,3)) * qJDD(3) - pkin(3) * t414 + (-t515 * t485 + t518 * t486) * qJD(1) + t520 * Ifges(3,5) - mrSges(3,1) * t509 - Ifges(6,3) * t506 - t480 * t460 - Ifges(4,5) * t493 - Ifges(4,6) * t494 - Ifges(5,6) * t471 - Ifges(5,5) * t472 + t479 * t461 + t462 * t444 - t463 * t443 + mrSges(3,3) * t470 - mrSges(4,1) * t454 + mrSges(4,2) * t455 - Ifges(6,6) * t440 - Ifges(6,5) * t441 - pkin(2) * t406 + Ifges(3,6) * qJDD(1);
t395 = mrSges(3,2) * t509 - mrSges(3,3) * t469 + Ifges(3,5) * qJDD(1) - t520 * Ifges(3,6) - pkin(6) * t406 - t515 * t397 + t518 * t398;
t394 = -mrSges(2,2) * g(1) - mrSges(2,3) * t499 + Ifges(2,5) * qJDD(1) - t520 * Ifges(2,6) - qJ(2) * t402 + t513 * t395 - t511 * t396;
t393 = mrSges(2,1) * g(1) + mrSges(2,3) * t498 + t520 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t532 + qJ(2) * t529 + t511 * t395 + t513 * t396;
t1 = [(-m(1) - m(2)) * g(1) + t532; -m(1) * g(2) + t524; -m(1) * g(3) + t530; pkin(1) * t402 + pkin(6) * t528 + t515 * t398 + t518 * t397 + pkin(2) * t521 + mrSges(3,1) * t469 - mrSges(3,2) * t470 + mrSges(2,1) * t499 - mrSges(2,2) * t498 - mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1); mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t530 - t519 * t393 - t516 * t394; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t524 - t516 * t393 + t519 * t394;];
tauB = t1;
