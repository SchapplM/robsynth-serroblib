% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRPP3
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRPP3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP3_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP3_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP3_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP3_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:15
% EndTime: 2019-12-31 20:53:17
% DurationCPUTime: 1.34s
% Computational Cost: add. (11952->229), mult. (15037->266), div. (0->0), fcn. (6460->6), ass. (0->90)
t538 = Ifges(4,1) + Ifges(5,2) + Ifges(6,3);
t526 = Ifges(4,4) + Ifges(5,6) - Ifges(6,6);
t525 = Ifges(4,5) - Ifges(5,4) + Ifges(6,5);
t537 = -Ifges(4,2) - Ifges(6,2) - Ifges(5,3);
t524 = Ifges(4,6) - Ifges(5,5) - Ifges(6,4);
t536 = Ifges(4,3) + Ifges(5,1) + Ifges(6,1);
t499 = sin(qJ(3));
t494 = qJD(1) + qJD(2);
t527 = qJD(3) * t494;
t520 = t499 * t527;
t531 = t494 * t499;
t534 = -2 * qJD(4);
t535 = pkin(3) * t520 + t531 * t534;
t533 = -2 * qJD(5);
t532 = -mrSges(6,1) - mrSges(4,3);
t502 = cos(qJ(3));
t530 = t494 * t502;
t501 = sin(qJ(1));
t504 = cos(qJ(1));
t485 = t501 * g(1) - t504 * g(2);
t476 = qJDD(1) * pkin(1) + t485;
t486 = -t504 * g(1) - t501 * g(2);
t506 = qJD(1) ^ 2;
t477 = -t506 * pkin(1) + t486;
t500 = sin(qJ(2));
t503 = cos(qJ(2));
t442 = t500 * t476 + t503 * t477;
t492 = t494 ^ 2;
t493 = qJDD(1) + qJDD(2);
t440 = -t492 * pkin(2) + t493 * pkin(7) + t442;
t436 = -t499 * g(3) + t502 * t440;
t466 = (-mrSges(4,1) * t502 + mrSges(4,2) * t499) * t494;
t469 = t502 * t493 - t520;
t478 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t531;
t464 = (-pkin(3) * t502 - qJ(4) * t499) * t494;
t505 = qJD(3) ^ 2;
t508 = -t505 * pkin(3) + qJDD(3) * qJ(4) + t464 * t530 + t436;
t433 = qJD(3) * t534 - t508;
t465 = (mrSges(5,2) * t502 - mrSges(5,3) * t499) * t494;
t484 = mrSges(5,1) * t531 + qJD(3) * mrSges(5,2);
t480 = pkin(4) * t531 - qJD(3) * qJ(5);
t498 = t502 ^ 2;
t431 = -t498 * t492 * qJ(5) + t469 * pkin(4) + qJDD(5) + ((2 * qJD(4)) + t480) * qJD(3) + t508;
t467 = (-mrSges(6,2) * t499 - mrSges(6,3) * t502) * t494;
t481 = mrSges(6,1) * t531 - qJD(3) * mrSges(6,3);
t515 = -m(6) * t431 - qJDD(3) * mrSges(6,2) - qJD(3) * t481 - t467 * t530;
t509 = -m(5) * t433 + qJDD(3) * mrSges(5,3) + qJD(3) * t484 + t465 * t530 - t515;
t423 = t466 * t530 + m(4) * t436 - qJDD(3) * mrSges(4,2) - qJD(3) * t478 + (mrSges(5,1) - t532) * t469 + t509;
t435 = -t502 * g(3) - t499 * t440;
t519 = t502 * t527;
t468 = t499 * t493 + t519;
t479 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t530;
t482 = -mrSges(5,1) * t530 - qJD(3) * mrSges(5,3);
t434 = -qJDD(3) * pkin(3) - t505 * qJ(4) + t464 * t531 + qJDD(4) - t435;
t430 = qJD(3) * t533 + (-t492 * t499 * t502 - qJDD(3)) * qJ(5) + (t468 - t519) * pkin(4) + t434;
t483 = mrSges(6,1) * t530 + qJD(3) * mrSges(6,2);
t514 = m(6) * t430 - qJDD(3) * mrSges(6,3) - qJD(3) * t483;
t511 = -m(5) * t434 - t468 * mrSges(5,1) - t514;
t528 = -t465 - t467;
t424 = m(4) * t435 + t532 * t468 + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t479 - t482) * qJD(3) + (-t466 + t528) * t531 + t511;
t516 = t502 * t423 - t499 * t424;
t416 = m(3) * t442 - t492 * mrSges(3,1) - t493 * mrSges(3,2) + t516;
t441 = t503 * t476 - t500 * t477;
t512 = -t493 * pkin(2) - t441;
t439 = -t492 * pkin(7) + t512;
t432 = -t469 * pkin(3) + (-t468 - t519) * qJ(4) + t439 + t535;
t428 = -t468 * qJ(4) + (-pkin(4) * t498 - pkin(7)) * t492 + (-pkin(3) - qJ(5)) * t469 + (-t480 * t499 + (-qJ(4) * qJD(3) + t533) * t502) * t494 + t512 + t535;
t513 = m(6) * t428 - t468 * mrSges(6,2) - t469 * mrSges(6,3) - t481 * t531 - t483 * t530;
t510 = -m(5) * t432 - t469 * mrSges(5,2) + t484 * t531 - t513;
t507 = -m(4) * t439 + t479 * t530 + t469 * mrSges(4,1) + (-t478 * t499 - t482 * t502) * t494 + (-mrSges(4,2) + mrSges(5,3)) * t468 + t510;
t419 = m(3) * t441 + t493 * mrSges(3,1) - t492 * mrSges(3,2) + t507;
t412 = t500 * t416 + t503 * t419;
t410 = m(2) * t485 + qJDD(1) * mrSges(2,1) - t506 * mrSges(2,2) + t412;
t517 = t503 * t416 - t500 * t419;
t411 = m(2) * t486 - t506 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t517;
t529 = t504 * t410 + t501 * t411;
t417 = t499 * t423 + t502 * t424;
t523 = (t525 * t499 + t524 * t502) * t494 + t536 * qJD(3);
t522 = (-t526 * t499 + t537 * t502) * t494 - t524 * qJD(3);
t521 = (t538 * t499 + t526 * t502) * t494 + t525 * qJD(3);
t518 = -t501 * t410 + t504 * t411;
t426 = t468 * mrSges(6,1) + t467 * t531 + t514;
t425 = -t468 * mrSges(5,3) + t482 * t530 - t510;
t413 = mrSges(5,1) * t434 + mrSges(6,1) * t430 + mrSges(4,2) * t439 - mrSges(6,2) * t428 - mrSges(4,3) * t435 - mrSges(5,3) * t432 + pkin(4) * t426 - qJ(4) * t425 + t522 * qJD(3) + t525 * qJDD(3) + t538 * t468 + t526 * t469 + t523 * t530;
t406 = -mrSges(4,1) * t439 + mrSges(4,3) * t436 - mrSges(5,1) * t433 + mrSges(5,2) * t432 + mrSges(6,1) * t431 - mrSges(6,3) * t428 - pkin(4) * t515 - qJ(5) * t513 - pkin(3) * t425 - t523 * t531 + (pkin(4) * mrSges(6,1) - t537) * t469 + t526 * t468 + t524 * qJDD(3) + t521 * qJD(3);
t405 = mrSges(3,1) * g(3) - pkin(2) * t417 + qJ(5) * t426 + mrSges(6,3) * t430 - mrSges(6,2) * t431 + mrSges(5,3) * t433 - mrSges(5,2) * t434 - mrSges(4,1) * t435 + mrSges(4,2) * t436 + mrSges(3,3) * t442 + t492 * Ifges(3,5) + Ifges(3,6) * t493 - qJ(4) * t509 - pkin(3) * (-qJD(3) * t482 + t511) + (pkin(3) * mrSges(6,1) - t525) * t468 + (pkin(3) * mrSges(5,2) - t536) * qJDD(3) + (-qJ(4) * (mrSges(5,1) + mrSges(6,1)) - t524) * t469 + (t521 * t502 + (-pkin(3) * t528 + t522) * t499) * t494;
t404 = -mrSges(3,2) * g(3) - mrSges(3,3) * t441 + Ifges(3,5) * t493 - t492 * Ifges(3,6) - pkin(7) * t417 - t499 * t406 + t502 * t413;
t403 = -mrSges(2,2) * g(3) - mrSges(2,3) * t485 + Ifges(2,5) * qJDD(1) - t506 * Ifges(2,6) - pkin(6) * t412 + t503 * t404 - t500 * t405;
t402 = Ifges(2,6) * qJDD(1) + t506 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t486 + t500 * t404 + t503 * t405 - pkin(1) * (-m(3) * g(3) + t417) + pkin(6) * t517;
t1 = [-m(1) * g(1) + t518; -m(1) * g(2) + t529; (-m(1) - m(2) - m(3)) * g(3) + t417; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t529 - t501 * t402 + t504 * t403; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t518 + t504 * t402 + t501 * t403; -mrSges(1,1) * g(2) + mrSges(2,1) * t485 + mrSges(3,1) * t441 + mrSges(1,2) * g(1) - mrSges(2,2) * t486 - mrSges(3,2) * t442 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * t493 + pkin(1) * t412 + pkin(2) * t507 + pkin(7) * t516 + t502 * t406 + t499 * t413;];
tauB = t1;
