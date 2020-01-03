% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPRR12
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPRR12_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR12_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR12_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR12_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR12_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR12_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:07:07
% EndTime: 2019-12-31 18:07:09
% DurationCPUTime: 2.05s
% Computational Cost: add. (19016->240), mult. (42079->292), div. (0->0), fcn. (27150->8), ass. (0->104)
t495 = sin(qJ(1));
t498 = cos(qJ(1));
t475 = t495 * g(1) - t498 * g(2);
t500 = qJD(1) ^ 2;
t507 = -t500 * qJ(2) + qJDD(2) - t475;
t529 = -pkin(1) - qJ(3);
t535 = -(2 * qJD(1) * qJD(3)) + t529 * qJDD(1) + t507;
t491 = sin(pkin(8));
t485 = t491 ^ 2;
t492 = cos(pkin(8));
t526 = t492 ^ 2 + t485;
t519 = t526 * mrSges(4,3);
t476 = -t498 * g(1) - t495 * g(2);
t534 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t476;
t494 = sin(qJ(4));
t497 = cos(qJ(4));
t511 = t491 * t497 + t492 * t494;
t471 = t511 * qJD(1);
t510 = -t491 * t494 + t492 * t497;
t472 = t510 * qJD(1);
t523 = t472 * qJD(4);
t457 = -t511 * qJDD(1) - t523;
t533 = pkin(3) * t500;
t532 = mrSges(2,1) - mrSges(3,2);
t531 = -Ifges(3,4) + Ifges(2,5);
t530 = -Ifges(2,6) + Ifges(3,5);
t528 = mrSges(4,2) * t492;
t454 = t491 * g(3) + t535 * t492;
t441 = (-pkin(6) * qJDD(1) - t491 * t533) * t492 + t454;
t455 = -t492 * g(3) + t535 * t491;
t521 = qJDD(1) * t491;
t442 = -pkin(6) * t521 - t485 * t533 + t455;
t430 = t494 * t441 + t497 * t442;
t451 = t471 * mrSges(5,1) + t472 * mrSges(5,2);
t466 = qJD(4) * mrSges(5,1) - t472 * mrSges(5,3);
t456 = t471 * pkin(4) - t472 * pkin(7);
t499 = qJD(4) ^ 2;
t427 = -t499 * pkin(4) + qJDD(4) * pkin(7) - t471 * t456 + t430;
t506 = qJDD(3) + t534;
t446 = pkin(3) * t521 + (-t526 * pkin(6) + t529) * t500 + t506;
t524 = t471 * qJD(4);
t458 = t510 * qJDD(1) - t524;
t428 = (-t458 + t524) * pkin(7) + (-t457 + t523) * pkin(4) + t446;
t493 = sin(qJ(5));
t496 = cos(qJ(5));
t424 = -t493 * t427 + t496 * t428;
t460 = t496 * qJD(4) - t493 * t472;
t437 = t460 * qJD(5) + t493 * qJDD(4) + t496 * t458;
t461 = t493 * qJD(4) + t496 * t472;
t439 = -t460 * mrSges(6,1) + t461 * mrSges(6,2);
t469 = qJD(5) + t471;
t443 = -t469 * mrSges(6,2) + t460 * mrSges(6,3);
t453 = qJDD(5) - t457;
t422 = m(6) * t424 + t453 * mrSges(6,1) - t437 * mrSges(6,3) - t461 * t439 + t469 * t443;
t425 = t496 * t427 + t493 * t428;
t436 = -t461 * qJD(5) + t496 * qJDD(4) - t493 * t458;
t444 = t469 * mrSges(6,1) - t461 * mrSges(6,3);
t423 = m(6) * t425 - t453 * mrSges(6,2) + t436 * mrSges(6,3) + t460 * t439 - t469 * t444;
t515 = -t493 * t422 + t496 * t423;
t413 = m(5) * t430 - qJDD(4) * mrSges(5,2) + t457 * mrSges(5,3) - qJD(4) * t466 - t471 * t451 + t515;
t429 = t497 * t441 - t494 * t442;
t465 = -qJD(4) * mrSges(5,2) - t471 * mrSges(5,3);
t426 = -qJDD(4) * pkin(4) - t499 * pkin(7) + t472 * t456 - t429;
t504 = -m(6) * t426 + t436 * mrSges(6,1) - t437 * mrSges(6,2) + t460 * t443 - t461 * t444;
t418 = m(5) * t429 + qJDD(4) * mrSges(5,1) - t458 * mrSges(5,3) + qJD(4) * t465 - t472 * t451 + t504;
t407 = t494 * t413 + t497 * t418;
t509 = -qJDD(1) * mrSges(4,3) - t500 * (mrSges(4,1) * t491 + t528);
t405 = m(4) * t454 + t509 * t492 + t407;
t516 = t497 * t413 - t494 * t418;
t406 = m(4) * t455 + t509 * t491 + t516;
t401 = t492 * t405 + t491 * t406;
t470 = -qJDD(1) * pkin(1) + t507;
t505 = -m(3) * t470 + t500 * mrSges(3,3) - t401;
t399 = m(2) * t475 - t500 * mrSges(2,2) + t532 * qJDD(1) + t505;
t468 = t500 * pkin(1) - t534;
t464 = t529 * t500 + t506;
t414 = t496 * t422 + t493 * t423;
t503 = m(5) * t446 - t457 * mrSges(5,1) + t458 * mrSges(5,2) + t471 * t465 + t472 * t466 + t414;
t502 = -m(4) * t464 - mrSges(4,1) * t521 - qJDD(1) * t528 - t503;
t501 = -m(3) * t468 + t500 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t502;
t410 = t501 + (-mrSges(2,1) - t519) * t500 - qJDD(1) * mrSges(2,2) + m(2) * t476;
t527 = t498 * t399 + t495 * t410;
t512 = Ifges(4,5) * t492 - Ifges(4,6) * t491;
t525 = t500 * t512;
t518 = -t495 * t399 + t498 * t410;
t517 = -t491 * t405 + t492 * t406;
t514 = Ifges(4,1) * t492 - Ifges(4,4) * t491;
t513 = Ifges(4,4) * t492 - Ifges(4,2) * t491;
t449 = Ifges(5,1) * t472 - Ifges(5,4) * t471 + Ifges(5,5) * qJD(4);
t448 = Ifges(5,4) * t472 - Ifges(5,2) * t471 + Ifges(5,6) * qJD(4);
t447 = Ifges(5,5) * t472 - Ifges(5,6) * t471 + Ifges(5,3) * qJD(4);
t433 = Ifges(6,1) * t461 + Ifges(6,4) * t460 + Ifges(6,5) * t469;
t432 = Ifges(6,4) * t461 + Ifges(6,2) * t460 + Ifges(6,6) * t469;
t431 = Ifges(6,5) * t461 + Ifges(6,6) * t460 + Ifges(6,3) * t469;
t416 = mrSges(6,2) * t426 - mrSges(6,3) * t424 + Ifges(6,1) * t437 + Ifges(6,4) * t436 + Ifges(6,5) * t453 + t460 * t431 - t469 * t432;
t415 = -mrSges(6,1) * t426 + mrSges(6,3) * t425 + Ifges(6,4) * t437 + Ifges(6,2) * t436 + Ifges(6,6) * t453 - t461 * t431 + t469 * t433;
t403 = -mrSges(5,1) * t446 - mrSges(6,1) * t424 + mrSges(6,2) * t425 + mrSges(5,3) * t430 + Ifges(5,4) * t458 - Ifges(6,5) * t437 + Ifges(5,2) * t457 + Ifges(5,6) * qJDD(4) - Ifges(6,6) * t436 - Ifges(6,3) * t453 - pkin(4) * t414 + qJD(4) * t449 - t461 * t432 + t460 * t433 - t472 * t447;
t402 = mrSges(5,2) * t446 - mrSges(5,3) * t429 + Ifges(5,1) * t458 + Ifges(5,4) * t457 + Ifges(5,5) * qJDD(4) - pkin(7) * t414 - qJD(4) * t448 - t493 * t415 + t496 * t416 - t471 * t447;
t400 = -m(3) * g(3) + t517;
t397 = mrSges(4,2) * t464 - mrSges(4,3) * t454 - pkin(6) * t407 + t514 * qJDD(1) + t497 * t402 - t494 * t403 - t491 * t525;
t396 = -mrSges(4,1) * t464 + mrSges(4,3) * t455 - pkin(3) * t503 + pkin(6) * t516 + t513 * qJDD(1) + t494 * t402 + t497 * t403 - t492 * t525;
t395 = Ifges(5,3) * qJDD(4) + pkin(7) * t515 + t493 * t416 + t496 * t415 + t471 * t449 + t472 * t448 - mrSges(2,3) * t475 + pkin(4) * t504 + mrSges(3,1) * t470 + mrSges(4,1) * t454 - mrSges(4,2) * t455 + Ifges(5,6) * t457 + Ifges(5,5) * t458 + mrSges(5,1) * t429 - mrSges(5,2) * t430 + pkin(3) * t407 + pkin(2) * t401 - qJ(2) * t400 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (t512 + t531) * qJDD(1) + (t491 * t514 + t492 * t513 + t530) * t500;
t394 = mrSges(2,3) * t476 - mrSges(3,1) * t468 - t491 * t397 - t492 * t396 - pkin(2) * t502 - qJ(3) * t517 - pkin(1) * t400 - t530 * qJDD(1) + t532 * g(3) + (-pkin(2) * t519 + t531) * t500;
t1 = [-m(1) * g(1) + t518; -m(1) * g(2) + t527; (-m(1) - m(2) - m(3)) * g(3) + t517; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t527 - t495 * t394 + t498 * t395; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t518 + t498 * t394 + t495 * t395; pkin(1) * t505 + qJ(2) * (-t500 * t519 + t501) + t492 * t397 - t491 * t396 - qJ(3) * t401 + mrSges(2,1) * t475 - mrSges(2,2) * t476 + mrSges(3,2) * t470 - mrSges(3,3) * t468 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB = t1;
