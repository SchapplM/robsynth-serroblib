% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% m [6x1]
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
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPRP1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP1_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP1_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP1_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:22
% EndTime: 2022-01-23 09:12:24
% DurationCPUTime: 2.03s
% Computational Cost: add. (15273->226), mult. (31695->283), div. (0->0), fcn. (17540->8), ass. (0->101)
t492 = sin(qJ(1));
t494 = cos(qJ(1));
t472 = t492 * g(1) - t494 * g(2);
t469 = qJDD(1) * pkin(1) + t472;
t473 = -t494 * g(1) - t492 * g(2);
t495 = qJD(1) ^ 2;
t470 = -t495 * pkin(1) + t473;
t488 = sin(pkin(7));
t490 = cos(pkin(7));
t442 = t488 * t469 + t490 * t470;
t540 = -t495 * pkin(2) + qJDD(1) * qJ(3) + (2 * qJD(1) * qJD(3)) + t442;
t539 = Ifges(5,1) + Ifges(6,1);
t532 = Ifges(5,4) + Ifges(6,4);
t531 = Ifges(5,5) + Ifges(6,5);
t538 = Ifges(5,2) + Ifges(6,2);
t537 = Ifges(5,6) + Ifges(6,6);
t536 = -Ifges(5,3) - Ifges(6,3);
t486 = -g(3) + qJDD(2);
t487 = sin(pkin(8));
t489 = cos(pkin(8));
t432 = t489 * t486 - t540 * t487;
t502 = -pkin(3) * t489 - pkin(6) * t487;
t468 = t502 * qJD(1);
t520 = t487 * qJD(1);
t430 = t468 * t520 - t432;
t491 = sin(qJ(4));
t493 = cos(qJ(4));
t517 = qJD(1) * qJD(4);
t460 = (-qJDD(1) * t491 - t493 * t517) * t487;
t461 = (qJDD(1) * t493 - t491 * t517) * t487;
t519 = t489 * qJD(1);
t475 = qJD(4) - t519;
t511 = t493 * t520;
t455 = t475 * pkin(4) - qJ(5) * t511;
t528 = t487 ^ 2 * t495;
t515 = t491 ^ 2 * t528;
t428 = -t460 * pkin(4) - qJ(5) * t515 + t455 * t511 + qJDD(5) + t430;
t509 = m(6) * t428 - t460 * mrSges(6,1);
t533 = -mrSges(5,2) - mrSges(6,2);
t535 = -m(5) * t430 + t460 * mrSges(5,1) + t533 * t461 - t509;
t534 = t489 ^ 2;
t529 = mrSges(4,2) * t487;
t433 = t487 * t486 + t540 * t489;
t466 = (-mrSges(4,1) * t489 + t529) * qJD(1);
t431 = t468 * t519 + t433;
t441 = t490 * t469 - t488 * t470;
t498 = -t495 * qJ(3) + qJDD(3) - t441;
t436 = (-pkin(2) + t502) * qJDD(1) + t498;
t435 = t493 * t436;
t426 = -t491 * t431 + t435;
t512 = t491 * t520;
t454 = -t475 * mrSges(5,2) - mrSges(5,3) * t512;
t516 = t489 * qJDD(1);
t474 = qJDD(4) - t516;
t458 = (mrSges(6,1) * t491 + mrSges(6,2) * t493) * t520;
t501 = (-t458 - (mrSges(5,1) * t491 + mrSges(5,2) * t493) * t520) * t520;
t503 = -0.2e1 * qJD(5) * t520;
t423 = t493 * t503 + t474 * pkin(4) - t461 * qJ(5) + t435 + (-pkin(4) * t493 * t528 - qJ(5) * t475 * t520 - t431) * t491;
t453 = -t475 * mrSges(6,2) - mrSges(6,3) * t512;
t514 = m(6) * t423 + t474 * mrSges(6,1) + t475 * t453;
t417 = m(5) * t426 + t474 * mrSges(5,1) + t475 * t454 + (-mrSges(5,3) - mrSges(6,3)) * t461 + t493 * t501 + t514;
t427 = t493 * t431 + t491 * t436;
t456 = t475 * mrSges(6,1) - mrSges(6,3) * t511;
t522 = -t475 * mrSges(5,1) + mrSges(5,3) * t511 - t456;
t425 = -pkin(4) * t515 + t460 * qJ(5) - t475 * t455 + t491 * t503 + t427;
t526 = m(6) * t425 + t460 * mrSges(6,3);
t418 = m(5) * t427 + t460 * mrSges(5,3) + t474 * t533 + t475 * t522 + t491 * t501 + t526;
t505 = -t491 * t417 + t493 * t418;
t521 = qJDD(1) * mrSges(4,3);
t413 = m(4) * t433 + (qJD(1) * t466 + t521) * t489 + t505;
t497 = t522 * t493 + (-t453 - t454) * t491;
t420 = m(4) * t432 + (-t521 + (-t466 + t497) * qJD(1)) * t487 + t535;
t506 = t489 * t413 - t487 * t420;
t406 = m(3) * t442 - t495 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t506;
t415 = t493 * t417 + t491 * t418;
t439 = -qJDD(1) * pkin(2) + t498;
t496 = -m(4) * t439 + mrSges(4,1) * t516 - t415 + (t495 * t534 + t528) * mrSges(4,3);
t410 = m(3) * t441 - t495 * mrSges(3,2) + (mrSges(3,1) - t529) * qJDD(1) + t496;
t402 = t488 * t406 + t490 * t410;
t400 = m(2) * t472 + qJDD(1) * mrSges(2,1) - t495 * mrSges(2,2) + t402;
t507 = t490 * t406 - t488 * t410;
t401 = m(2) * t473 - t495 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t507;
t527 = t494 * t400 + t492 * t401;
t407 = t487 * t413 + t489 * t420;
t525 = (t537 * t491 - t531 * t493) * t520 + t536 * t475;
t524 = (t538 * t491 - t532 * t493) * t520 - t537 * t475;
t523 = (t532 * t491 - t539 * t493) * t520 - t531 * t475;
t513 = m(3) * t486 + t407;
t508 = -t492 * t400 + t494 * t401;
t500 = Ifges(4,5) * t487 + Ifges(4,6) * t489;
t467 = t500 * qJD(1);
t421 = -t461 * mrSges(6,3) - t458 * t511 + t514;
t414 = mrSges(5,2) * t430 + mrSges(6,2) * t428 - mrSges(5,3) * t426 - mrSges(6,3) * t423 - qJ(5) * t421 + t532 * t460 + t539 * t461 + t531 * t474 + t524 * t475 + t525 * t512;
t408 = -mrSges(5,1) * t430 + mrSges(5,3) * t427 - mrSges(6,1) * t428 + mrSges(6,3) * t425 - pkin(4) * t509 + qJ(5) * t526 + (-qJ(5) * t456 - t523) * t475 + (-qJ(5) * mrSges(6,2) + t537) * t474 + (-pkin(4) * mrSges(6,2) + t532) * t461 + t538 * t460 + ((-pkin(4) * t453 - qJ(5) * t458) * t491 + (-pkin(4) * t456 + t525) * t493) * t520;
t403 = Ifges(4,2) * t516 - mrSges(4,1) * t439 - mrSges(5,1) * t426 - mrSges(6,1) * t423 + mrSges(5,2) * t427 + mrSges(6,2) * t425 + mrSges(4,3) * t433 - pkin(3) * t415 - pkin(4) * t421 + t536 * t474 - t531 * t461 - t537 * t460 + (Ifges(4,4) * qJDD(1) + (t491 * t523 + t493 * t524 - t467) * qJD(1)) * t487;
t396 = t467 * t519 + mrSges(4,2) * t439 - mrSges(4,3) * t432 - pkin(6) * t415 - t491 * t408 + t493 * t414 + (Ifges(4,1) * t487 + Ifges(4,4) * t489) * qJDD(1);
t395 = t495 * Ifges(3,5) - mrSges(3,1) * t486 + mrSges(3,3) * t442 - mrSges(4,1) * t432 + mrSges(4,2) * t433 - t491 * t414 - t493 * t408 - pkin(3) * t535 - pkin(6) * t505 - pkin(2) * t407 + (Ifges(3,6) - t500) * qJDD(1) + (Ifges(4,4) * t534 * qJD(1) + (-pkin(3) * t497 + (-Ifges(4,4) * t487 + (Ifges(4,1) - Ifges(4,2)) * t489) * qJD(1)) * t487) * qJD(1);
t394 = mrSges(3,2) * t486 - mrSges(3,3) * t441 + Ifges(3,5) * qJDD(1) - t495 * Ifges(3,6) - qJ(3) * t407 + t489 * t396 - t487 * t403;
t393 = -mrSges(2,2) * g(3) - mrSges(2,3) * t472 + Ifges(2,5) * qJDD(1) - t495 * Ifges(2,6) - qJ(2) * t402 + t490 * t394 - t488 * t395;
t392 = mrSges(2,1) * g(3) + mrSges(2,3) * t473 + t495 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t513 + qJ(2) * t507 + t488 * t394 + t490 * t395;
t1 = [-m(1) * g(1) + t508; -m(1) * g(2) + t527; (-m(1) - m(2)) * g(3) + t513; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t527 - t492 * t392 + t494 * t393; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t508 + t494 * t392 + t492 * t393; pkin(1) * t402 + mrSges(2,1) * t472 - mrSges(2,2) * t473 + qJ(3) * t506 + t487 * t396 + t489 * t403 + pkin(2) * t496 + mrSges(3,1) * t441 - mrSges(3,2) * t442 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(2) * t529 + Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB = t1;
