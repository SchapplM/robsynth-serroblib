% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRRR9
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRRR9_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR9_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR9_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR9_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR9_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR9_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR9_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:19:42
% EndTime: 2019-12-05 17:19:51
% DurationCPUTime: 6.60s
% Computational Cost: add. (83105->267), mult. (159827->347), div. (0->0), fcn. (110866->12), ass. (0->115)
t504 = sin(pkin(10));
t506 = cos(pkin(10));
t495 = t504 * g(1) - t506 * g(2);
t496 = -t506 * g(1) - t504 * g(2);
t503 = -g(3) + qJDD(1);
t505 = sin(pkin(5));
t507 = cos(pkin(5));
t511 = sin(qJ(2));
t515 = cos(qJ(2));
t460 = -t511 * t496 + (t495 * t507 + t503 * t505) * t515;
t517 = qJD(2) ^ 2;
t532 = t507 * t511;
t533 = t505 * t511;
t461 = t495 * t532 + t515 * t496 + t503 * t533;
t457 = -t517 * pkin(2) + qJDD(2) * pkin(7) + t461;
t476 = -t505 * t495 + t507 * t503;
t510 = sin(qJ(3));
t514 = cos(qJ(3));
t452 = t514 * t457 + t510 * t476;
t492 = (-pkin(3) * t514 - pkin(8) * t510) * qJD(2);
t516 = qJD(3) ^ 2;
t529 = t514 * qJD(2);
t443 = -t516 * pkin(3) + qJDD(3) * pkin(8) + t492 * t529 + t452;
t456 = -qJDD(2) * pkin(2) - t517 * pkin(7) - t460;
t528 = qJD(2) * qJD(3);
t527 = t514 * t528;
t493 = t510 * qJDD(2) + t527;
t502 = t510 * t528;
t494 = t514 * qJDD(2) - t502;
t446 = (-t493 - t527) * pkin(8) + (-t494 + t502) * pkin(3) + t456;
t509 = sin(qJ(4));
t513 = cos(qJ(4));
t435 = -t509 * t443 + t513 * t446;
t530 = qJD(2) * t510;
t489 = t513 * qJD(3) - t509 * t530;
t468 = t489 * qJD(4) + t509 * qJDD(3) + t513 * t493;
t486 = qJDD(4) - t494;
t490 = t509 * qJD(3) + t513 * t530;
t501 = qJD(4) - t529;
t433 = (t489 * t501 - t468) * pkin(9) + (t489 * t490 + t486) * pkin(4) + t435;
t436 = t513 * t443 + t509 * t446;
t467 = -t490 * qJD(4) + t513 * qJDD(3) - t509 * t493;
t475 = t501 * pkin(4) - t490 * pkin(9);
t485 = t489 ^ 2;
t434 = -t485 * pkin(4) + t467 * pkin(9) - t501 * t475 + t436;
t508 = sin(qJ(5));
t512 = cos(qJ(5));
t431 = t512 * t433 - t508 * t434;
t469 = t512 * t489 - t508 * t490;
t440 = t469 * qJD(5) + t508 * t467 + t512 * t468;
t470 = t508 * t489 + t512 * t490;
t453 = -t469 * mrSges(6,1) + t470 * mrSges(6,2);
t500 = qJD(5) + t501;
t458 = -t500 * mrSges(6,2) + t469 * mrSges(6,3);
t482 = qJDD(5) + t486;
t429 = m(6) * t431 + t482 * mrSges(6,1) - t440 * mrSges(6,3) - t470 * t453 + t500 * t458;
t432 = t508 * t433 + t512 * t434;
t439 = -t470 * qJD(5) + t512 * t467 - t508 * t468;
t459 = t500 * mrSges(6,1) - t470 * mrSges(6,3);
t430 = m(6) * t432 - t482 * mrSges(6,2) + t439 * mrSges(6,3) + t469 * t453 - t500 * t459;
t421 = t512 * t429 + t508 * t430;
t471 = -t489 * mrSges(5,1) + t490 * mrSges(5,2);
t473 = -t501 * mrSges(5,2) + t489 * mrSges(5,3);
t419 = m(5) * t435 + t486 * mrSges(5,1) - t468 * mrSges(5,3) - t490 * t471 + t501 * t473 + t421;
t474 = t501 * mrSges(5,1) - t490 * mrSges(5,3);
t523 = -t508 * t429 + t512 * t430;
t420 = m(5) * t436 - t486 * mrSges(5,2) + t467 * mrSges(5,3) + t489 * t471 - t501 * t474 + t523;
t417 = t513 * t419 + t509 * t420;
t497 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t530;
t498 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t529;
t519 = -m(4) * t456 + t494 * mrSges(4,1) - t493 * mrSges(4,2) - t497 * t530 + t498 * t529 - t417;
t413 = m(3) * t460 + qJDD(2) * mrSges(3,1) - t517 * mrSges(3,2) + t519;
t534 = t413 * t515;
t491 = (-mrSges(4,1) * t514 + mrSges(4,2) * t510) * qJD(2);
t524 = -t509 * t419 + t513 * t420;
t416 = m(4) * t452 - qJDD(3) * mrSges(4,2) + t494 * mrSges(4,3) - qJD(3) * t497 + t491 * t529 + t524;
t451 = -t510 * t457 + t514 * t476;
t442 = -qJDD(3) * pkin(3) - t516 * pkin(8) + t492 * t530 - t451;
t437 = -t467 * pkin(4) - t485 * pkin(9) + t490 * t475 + t442;
t520 = m(6) * t437 - t439 * mrSges(6,1) + t440 * mrSges(6,2) - t469 * t458 + t470 * t459;
t518 = -m(5) * t442 + t467 * mrSges(5,1) - t468 * mrSges(5,2) + t489 * t473 - t490 * t474 - t520;
t425 = m(4) * t451 + qJDD(3) * mrSges(4,1) - t493 * mrSges(4,3) + qJD(3) * t498 - t491 * t530 + t518;
t525 = t514 * t416 - t510 * t425;
t405 = m(3) * t461 - t517 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t525;
t409 = t510 * t416 + t514 * t425;
t408 = m(3) * t476 + t409;
t395 = t405 * t532 - t505 * t408 + t507 * t534;
t393 = m(2) * t495 + t395;
t400 = t515 * t405 - t511 * t413;
t399 = m(2) * t496 + t400;
t531 = t506 * t393 + t504 * t399;
t394 = t405 * t533 + t507 * t408 + t505 * t534;
t526 = -t504 * t393 + t506 * t399;
t447 = Ifges(6,5) * t470 + Ifges(6,6) * t469 + Ifges(6,3) * t500;
t449 = Ifges(6,1) * t470 + Ifges(6,4) * t469 + Ifges(6,5) * t500;
t422 = -mrSges(6,1) * t437 + mrSges(6,3) * t432 + Ifges(6,4) * t440 + Ifges(6,2) * t439 + Ifges(6,6) * t482 - t470 * t447 + t500 * t449;
t448 = Ifges(6,4) * t470 + Ifges(6,2) * t469 + Ifges(6,6) * t500;
t423 = mrSges(6,2) * t437 - mrSges(6,3) * t431 + Ifges(6,1) * t440 + Ifges(6,4) * t439 + Ifges(6,5) * t482 + t469 * t447 - t500 * t448;
t462 = Ifges(5,5) * t490 + Ifges(5,6) * t489 + Ifges(5,3) * t501;
t464 = Ifges(5,1) * t490 + Ifges(5,4) * t489 + Ifges(5,5) * t501;
t406 = -mrSges(5,1) * t442 + mrSges(5,3) * t436 + Ifges(5,4) * t468 + Ifges(5,2) * t467 + Ifges(5,6) * t486 - pkin(4) * t520 + pkin(9) * t523 + t512 * t422 + t508 * t423 - t490 * t462 + t501 * t464;
t463 = Ifges(5,4) * t490 + Ifges(5,2) * t489 + Ifges(5,6) * t501;
t410 = mrSges(5,2) * t442 - mrSges(5,3) * t435 + Ifges(5,1) * t468 + Ifges(5,4) * t467 + Ifges(5,5) * t486 - pkin(9) * t421 - t508 * t422 + t512 * t423 + t489 * t462 - t501 * t463;
t479 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t510 + Ifges(4,6) * t514) * qJD(2);
t480 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t510 + Ifges(4,2) * t514) * qJD(2);
t396 = mrSges(4,2) * t456 - mrSges(4,3) * t451 + Ifges(4,1) * t493 + Ifges(4,4) * t494 + Ifges(4,5) * qJDD(3) - pkin(8) * t417 - qJD(3) * t480 - t509 * t406 + t513 * t410 + t479 * t529;
t481 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t510 + Ifges(4,4) * t514) * qJD(2);
t401 = Ifges(4,4) * t493 + Ifges(4,2) * t494 + Ifges(4,6) * qJDD(3) - t479 * t530 + qJD(3) * t481 - mrSges(4,1) * t456 + mrSges(4,3) * t452 - Ifges(5,5) * t468 - Ifges(5,6) * t467 - Ifges(5,3) * t486 - t490 * t463 + t489 * t464 - mrSges(5,1) * t435 + mrSges(5,2) * t436 - Ifges(6,5) * t440 - Ifges(6,6) * t439 - Ifges(6,3) * t482 - t470 * t448 + t469 * t449 - mrSges(6,1) * t431 + mrSges(6,2) * t432 - pkin(4) * t421 - pkin(3) * t417;
t390 = mrSges(3,2) * t476 - mrSges(3,3) * t460 + Ifges(3,5) * qJDD(2) - t517 * Ifges(3,6) - pkin(7) * t409 + t514 * t396 - t510 * t401;
t391 = Ifges(3,6) * qJDD(2) + t517 * Ifges(3,5) - mrSges(3,1) * t476 + mrSges(3,3) * t461 - Ifges(4,5) * t493 - Ifges(4,6) * t494 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t451 + mrSges(4,2) * t452 - t509 * t410 - t513 * t406 - pkin(3) * t518 - pkin(8) * t524 - pkin(2) * t409 + (-t510 * t480 + t514 * t481) * qJD(2);
t521 = pkin(6) * t400 + t390 * t511 + t391 * t515;
t389 = mrSges(3,1) * t460 - mrSges(3,2) * t461 + Ifges(3,3) * qJDD(2) + pkin(2) * t519 + pkin(7) * t525 + t510 * t396 + t514 * t401;
t388 = mrSges(2,2) * t503 - mrSges(2,3) * t495 + t515 * t390 - t511 * t391 + (-t394 * t505 - t395 * t507) * pkin(6);
t387 = -mrSges(2,1) * t503 + mrSges(2,3) * t496 - pkin(1) * t394 - t505 * t389 + t521 * t507;
t1 = [-m(1) * g(1) + t526; -m(1) * g(2) + t531; -m(1) * g(3) + m(2) * t503 + t394; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t531 - t504 * t387 + t506 * t388; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t526 + t506 * t387 + t504 * t388; -mrSges(1,1) * g(2) + mrSges(2,1) * t495 + mrSges(1,2) * g(1) - mrSges(2,2) * t496 + pkin(1) * t395 + t507 * t389 + t521 * t505;];
tauB = t1;
