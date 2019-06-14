% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 18:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRPP2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:01:40
% EndTime: 2019-05-07 18:01:45
% DurationCPUTime: 3.02s
% Computational Cost: add. (25829->306), mult. (51701->361), div. (0->0), fcn. (35728->8), ass. (0->121)
t585 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t563 = Ifges(5,4) - Ifges(6,5) - Ifges(7,4);
t562 = Ifges(5,5) + Ifges(6,4) - Ifges(7,5);
t584 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t561 = -Ifges(6,6) + Ifges(7,6) + Ifges(5,6);
t583 = Ifges(5,3) + Ifges(6,2) + Ifges(7,3);
t536 = sin(qJ(3));
t537 = sin(qJ(2));
t539 = cos(qJ(2));
t577 = cos(qJ(3));
t516 = (t536 * t539 + t577 * t537) * qJD(1);
t565 = qJD(1) * qJD(2);
t522 = qJDD(1) * t537 + t539 * t565;
t564 = qJDD(1) * t539;
t547 = -t537 * t565 + t564;
t486 = -qJD(3) * t516 - t522 * t536 + t577 * t547;
t566 = qJD(1) * t539;
t567 = qJD(1) * t537;
t515 = -t536 * t567 + t577 * t566;
t487 = t515 * qJD(3) + t577 * t522 + t536 * t547;
t525 = qJD(2) * pkin(2) - pkin(8) * t567;
t534 = t539 ^ 2;
t541 = qJD(1) ^ 2;
t538 = sin(qJ(1));
t540 = cos(qJ(1));
t555 = t538 * g(1) - t540 * g(2);
t549 = -qJDD(1) * pkin(1) - t555;
t488 = -t547 * pkin(2) + t525 * t567 + (-pkin(8) * t534 - pkin(7)) * t541 + t549;
t533 = qJD(2) + qJD(3);
t433 = (-t515 * t533 - t487) * pkin(9) + (t516 * t533 - t486) * pkin(3) + t488;
t550 = -g(1) * t540 - g(2) * t538;
t518 = -pkin(1) * t541 + qJDD(1) * pkin(7) + t550;
t569 = t518 * t537;
t574 = pkin(2) * t541;
t477 = qJDD(2) * pkin(2) - pkin(8) * t522 - t569 + (pkin(8) * t565 + t537 * t574 - g(3)) * t539;
t505 = -t537 * g(3) + t539 * t518;
t478 = t547 * pkin(8) - qJD(2) * t525 - t534 * t574 + t505;
t440 = t536 * t477 + t577 * t478;
t500 = -pkin(3) * t515 - pkin(9) * t516;
t531 = t533 ^ 2;
t532 = qJDD(2) + qJDD(3);
t437 = -pkin(3) * t531 + pkin(9) * t532 + t500 * t515 + t440;
t535 = sin(qJ(4));
t576 = cos(qJ(4));
t430 = t576 * t433 - t535 * t437;
t502 = t516 * t535 - t576 * t533;
t503 = t576 * t516 + t535 * t533;
t472 = pkin(4) * t502 - qJ(5) * t503;
t485 = qJDD(4) - t486;
t511 = qJD(4) - t515;
t510 = t511 ^ 2;
t428 = -t485 * pkin(4) - t510 * qJ(5) + t503 * t472 + qJDD(5) - t430;
t489 = -mrSges(6,2) * t502 + mrSges(6,3) * t511;
t582 = -m(6) * t428 + t485 * mrSges(6,1) + t511 * t489;
t448 = -t502 * qJD(4) + t576 * t487 + t535 * t532;
t439 = t577 * t477 - t536 * t478;
t546 = pkin(3) * t532 + pkin(9) * t531 - t516 * t500 + t439;
t570 = t502 * t511;
t581 = (-t448 + t570) * qJ(5) - t546;
t490 = mrSges(7,2) * t511 + mrSges(7,3) * t502;
t579 = -0.2e1 * t503;
t421 = qJD(6) * t579 + (-t448 - t570) * qJ(6) + (t502 * t503 - t485) * pkin(5) + t428;
t474 = -mrSges(7,1) * t502 + mrSges(7,2) * t503;
t551 = -m(7) * t421 + t448 * mrSges(7,3) + t503 * t474;
t419 = -mrSges(7,1) * t485 - t490 * t511 - t551;
t473 = mrSges(6,1) * t502 - mrSges(6,3) * t503;
t416 = t448 * mrSges(6,2) + t473 * t503 + t419 - t582;
t431 = t535 * t433 + t576 * t437;
t578 = 2 * qJD(5);
t427 = -pkin(4) * t510 + t485 * qJ(5) - t502 * t472 + t511 * t578 + t431;
t447 = qJD(4) * t503 + t487 * t535 - t576 * t532;
t492 = -pkin(5) * t511 - qJ(6) * t503;
t501 = t502 ^ 2;
t423 = -pkin(5) * t501 + qJ(6) * t447 + 0.2e1 * qJD(6) * t502 + t492 * t511 + t427;
t493 = -mrSges(7,1) * t511 - mrSges(7,3) * t503;
t495 = -mrSges(6,1) * t511 + mrSges(6,2) * t503;
t559 = m(7) * t423 + t447 * mrSges(7,3) + t502 * t474;
t548 = m(6) * t427 + t485 * mrSges(6,3) + t511 * t495 + t559;
t556 = -t563 * t502 + t585 * t503 + t562 * t511;
t557 = t584 * t502 + t563 * t503 + t561 * t511;
t580 = -t561 * t447 + t562 * t448 + t583 * t485 + t556 * t502 + t557 * t503 + mrSges(5,1) * t430 - mrSges(6,1) * t428 - mrSges(7,1) * t421 - mrSges(5,2) * t431 + mrSges(7,2) * t423 + mrSges(6,3) * t427 - pkin(4) * t416 - pkin(5) * t419 + qJ(5) * (-t447 * mrSges(6,2) + mrSges(7,2) * t485 - t473 * t502 + t493 * t511 + t548);
t573 = -mrSges(5,3) - mrSges(6,2);
t572 = Ifges(3,6) * qJD(2);
t571 = qJD(2) * mrSges(3,1);
t499 = -mrSges(4,1) * t515 + mrSges(4,2) * t516;
t507 = mrSges(4,1) * t533 - mrSges(4,3) * t516;
t491 = -mrSges(5,2) * t511 - mrSges(5,3) * t502;
t568 = -mrSges(5,1) * t502 - mrSges(5,2) * t503 - t473;
t409 = m(5) * t430 + (t490 + t491) * t511 + t568 * t503 + (mrSges(5,1) + mrSges(7,1)) * t485 + t573 * t448 + t551 + t582;
t494 = mrSges(5,1) * t511 - mrSges(5,3) * t503;
t414 = m(5) * t431 + (t493 - t494) * t511 + t568 * t502 + (-mrSges(5,2) + mrSges(7,2)) * t485 + t573 * t447 + t548;
t553 = -t409 * t535 + t576 * t414;
t405 = m(4) * t440 - mrSges(4,2) * t532 + mrSges(4,3) * t486 + t499 * t515 - t507 * t533 + t553;
t506 = -mrSges(4,2) * t533 + mrSges(4,3) * t515;
t425 = -qJ(6) * t501 + qJDD(6) + (-pkin(4) - pkin(5)) * t447 + (-pkin(4) * t511 + t492 + t578) * t503 - t581;
t420 = m(7) * t425 - t447 * mrSges(7,1) + t448 * mrSges(7,2) - t502 * t490 + t503 * t493;
t429 = qJD(5) * t579 + (t503 * t511 + t447) * pkin(4) + t581;
t418 = m(6) * t429 + mrSges(6,1) * t447 - t448 * mrSges(6,3) + t489 * t502 - t503 * t495 - t420;
t542 = m(5) * t546 - t447 * mrSges(5,1) - mrSges(5,2) * t448 - t502 * t491 - t494 * t503 - t418;
t411 = m(4) * t439 + mrSges(4,1) * t532 - mrSges(4,3) * t487 - t499 * t516 + t506 * t533 + t542;
t400 = t536 * t405 + t577 * t411;
t407 = t576 * t409 + t535 * t414;
t558 = t561 * t502 - t562 * t503 - t583 * t511;
t554 = t577 * t405 - t536 * t411;
t399 = mrSges(5,1) * t546 + mrSges(5,3) * t431 - mrSges(6,1) * t429 + mrSges(6,2) * t427 + mrSges(7,1) * t425 - mrSges(7,3) * t423 + pkin(5) * t420 - qJ(6) * t559 - pkin(4) * t418 + (-qJ(6) * t493 + t556) * t511 + t558 * t503 + (-mrSges(7,2) * qJ(6) + t561) * t485 + t563 * t448 + t584 * t447;
t402 = -mrSges(5,2) * t546 + mrSges(6,2) * t428 + mrSges(7,2) * t425 - mrSges(5,3) * t430 - mrSges(6,3) * t429 - mrSges(7,3) * t421 - qJ(5) * t418 - qJ(6) * t419 - t563 * t447 + t585 * t448 + t562 * t485 + t558 * t502 - t557 * t511;
t497 = Ifges(4,4) * t516 + Ifges(4,2) * t515 + Ifges(4,6) * t533;
t498 = Ifges(4,1) * t516 + Ifges(4,4) * t515 + Ifges(4,5) * t533;
t545 = mrSges(4,1) * t439 - mrSges(4,2) * t440 + Ifges(4,5) * t487 + Ifges(4,6) * t486 + Ifges(4,3) * t532 + pkin(3) * t542 + pkin(9) * t553 + t576 * t399 + t535 * t402 + t516 * t497 - t515 * t498;
t544 = m(4) * t488 - t486 * mrSges(4,1) + t487 * mrSges(4,2) - t515 * t506 + t516 * t507 + t407;
t524 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t566;
t523 = -mrSges(3,3) * t567 + t571;
t521 = (-mrSges(3,1) * t539 + mrSges(3,2) * t537) * qJD(1);
t517 = -t541 * pkin(7) + t549;
t514 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t537 + Ifges(3,4) * t539) * qJD(1);
t513 = t572 + (Ifges(3,4) * t537 + Ifges(3,2) * t539) * qJD(1);
t504 = -g(3) * t539 - t569;
t496 = Ifges(4,5) * t516 + Ifges(4,6) * t515 + Ifges(4,3) * t533;
t397 = -mrSges(4,1) * t488 + mrSges(4,3) * t440 + Ifges(4,4) * t487 + Ifges(4,2) * t486 + Ifges(4,6) * t532 - pkin(3) * t407 - t516 * t496 + t533 * t498 - t580;
t396 = mrSges(4,2) * t488 - mrSges(4,3) * t439 + Ifges(4,1) * t487 + Ifges(4,4) * t486 + Ifges(4,5) * t532 - pkin(9) * t407 - t535 * t399 + t576 * t402 + t515 * t496 - t533 * t497;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t555 - mrSges(2,2) * t550 + t537 * (mrSges(3,2) * t517 - mrSges(3,3) * t504 + Ifges(3,1) * t522 + Ifges(3,4) * t547 + Ifges(3,5) * qJDD(2) - pkin(8) * t400 - qJD(2) * t513 + t577 * t396 - t536 * t397) + t539 * (-mrSges(3,1) * t517 + mrSges(3,3) * t505 + Ifges(3,4) * t522 + Ifges(3,2) * t547 + Ifges(3,6) * qJDD(2) - pkin(2) * t544 + pkin(8) * t554 + qJD(2) * t514 + t536 * t396 + t577 * t397) + pkin(1) * (mrSges(3,1) * t564 - m(3) * t517 - t522 * mrSges(3,2) + (t539 * t524 + (-t523 - t571) * t537) * qJD(1) - t544) + pkin(7) * (t539 * (m(3) * t505 - qJDD(2) * mrSges(3,2) + t547 * mrSges(3,3) - qJD(2) * t523 + t521 * t566 + t554) - t537 * (m(3) * t504 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t522 + qJD(2) * t524 - t521 * t567 + t400)); Ifges(3,6) * t564 + Ifges(3,3) * qJDD(2) + pkin(2) * t400 + t545 + (-t539 * t514 + (t513 - t572) * t537) * qJD(1) + mrSges(3,1) * t504 - mrSges(3,2) * t505 + Ifges(3,5) * t522; t545; t580; t416; t420;];
tauJ  = t1;
