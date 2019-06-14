% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
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
% Datum: 2019-05-05 16:43
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPPR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:41:39
% EndTime: 2019-05-05 16:41:42
% DurationCPUTime: 1.90s
% Computational Cost: add. (4523->258), mult. (8998->297), div. (0->0), fcn. (4292->8), ass. (0->110)
t521 = sin(qJ(3));
t524 = cos(qJ(3));
t550 = Ifges(4,4) - Ifges(5,5) + Ifges(6,4);
t575 = t521 * (Ifges(4,1) + Ifges(5,1) + Ifges(6,2)) + t524 * t550;
t574 = -t524 * (Ifges(4,2) + Ifges(5,3) + Ifges(6,1)) - t521 * t550;
t549 = Ifges(4,5) + Ifges(5,4) + Ifges(6,6);
t548 = Ifges(4,6) - Ifges(5,6) + Ifges(6,5);
t569 = t521 * (t574 * qJD(1) - t548 * qJD(3));
t522 = sin(qJ(1));
t525 = cos(qJ(1));
t542 = t522 * g(1) - g(2) * t525;
t482 = qJDD(1) * pkin(1) + t542;
t527 = qJD(1) ^ 2;
t540 = -g(1) * t525 - g(2) * t522;
t487 = -pkin(1) * t527 + t540;
t517 = sin(pkin(9));
t518 = cos(pkin(9));
t444 = t517 * t482 + t518 * t487;
t433 = -pkin(2) * t527 + qJDD(1) * pkin(7) + t444;
t515 = -g(3) + qJDD(2);
t429 = t524 * t433 + t521 * t515;
t483 = (-t524 * pkin(3) - t521 * qJ(4)) * qJD(1);
t553 = qJD(1) * t524;
t568 = qJDD(3) * qJ(4) + t483 * t553 + t429;
t552 = qJD(1) * qJD(3);
t544 = t524 * t552;
t489 = qJDD(1) * t521 + t544;
t551 = qJD(1) * qJD(5);
t567 = -0.2e1 * t521 * t551 + (-t489 + t544) * qJ(5);
t526 = qJD(3) ^ 2;
t428 = -t521 * t433 + t515 * t524;
t554 = qJD(1) * t521;
t538 = t483 * t554 + qJDD(4) - t428;
t426 = -qJDD(3) * pkin(3) - qJ(4) * t526 + t538;
t500 = mrSges(5,2) * t553 + qJD(3) * mrSges(5,3);
t566 = m(5) * t426 - qJDD(3) * mrSges(5,1) - qJD(3) * t500;
t565 = -2 * qJD(4);
t564 = 2 * qJD(4);
t563 = -pkin(3) - pkin(8);
t562 = pkin(4) + pkin(8);
t560 = pkin(3) * t526;
t559 = mrSges(4,3) + mrSges(5,2);
t557 = t524 ^ 2 * t527;
t556 = t524 * t527;
t543 = t521 * t552;
t490 = qJDD(1) * t524 - t543;
t494 = -qJD(3) * pkin(4) - qJ(5) * t554;
t443 = t518 * t482 - t517 * t487;
t432 = -qJDD(1) * pkin(2) - t527 * pkin(7) - t443;
t537 = -t490 * pkin(3) + t432 + (-t489 - t544) * qJ(4);
t530 = -qJ(5) * t557 + qJDD(5) - t537 + (t494 + t564) * t554;
t415 = t562 * t490 + pkin(5) * t489 + (pkin(5) * t524 + t563 * t521) * t552 + t530;
t488 = (t521 * pkin(5) + t524 * pkin(8)) * qJD(1);
t418 = (-pkin(5) - qJ(4)) * t526 + (-pkin(4) * t556 - qJD(1) * t488) * t521 + (-pkin(3) - t562) * qJDD(3) + t538 + t567;
t520 = sin(qJ(6));
t523 = cos(qJ(6));
t413 = t415 * t523 - t418 * t520;
t480 = -qJD(3) * t523 + t520 * t553;
t442 = qJD(6) * t480 - qJDD(3) * t520 - t490 * t523;
t481 = -qJD(3) * t520 - t523 * t553;
t445 = -mrSges(7,1) * t480 + mrSges(7,2) * t481;
t503 = qJD(6) + t554;
t446 = -mrSges(7,2) * t503 + mrSges(7,3) * t480;
t478 = qJDD(6) + t489;
t410 = m(7) * t413 + mrSges(7,1) * t478 - mrSges(7,3) * t442 - t445 * t481 + t446 * t503;
t414 = t415 * t520 + t418 * t523;
t441 = -qJD(6) * t481 - qJDD(3) * t523 + t490 * t520;
t447 = mrSges(7,1) * t503 - mrSges(7,3) * t481;
t411 = m(7) * t414 - mrSges(7,2) * t478 + mrSges(7,3) * t441 + t445 * t480 - t447 * t503;
t401 = t523 * t410 + t520 * t411;
t555 = -t520 * t410 + t523 * t411;
t546 = t575 * qJD(1) + t549 * qJD(3);
t545 = pkin(1) * t518 + pkin(2);
t484 = (-t524 * mrSges(5,1) - t521 * mrSges(5,3)) * qJD(1);
t485 = (-t524 * mrSges(4,1) + t521 * mrSges(4,2)) * qJD(1);
t498 = -qJD(3) * mrSges(6,1) + mrSges(6,3) * t553;
t499 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t553;
t422 = (-t521 * t556 - qJDD(3)) * pkin(4) + t426 + t567;
t486 = (t521 * mrSges(6,1) - t524 * mrSges(6,2)) * qJD(1);
t539 = -m(6) * t422 + t486 * t554 - t555;
t396 = m(4) * t428 + (mrSges(4,1) - mrSges(6,2)) * qJDD(3) + (-t498 + t499) * qJD(3) + (-t484 - t485) * t554 + (mrSges(6,3) - t559) * t489 + t539 - t566;
t496 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t554;
t508 = qJD(3) * t564;
t425 = t508 - t560 + t568;
t497 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t554;
t535 = pkin(4) * t557 + qJ(5) * t490 - t568;
t421 = 0.2e1 * t524 * t551 + t560 + (t565 - t494) * qJD(3) + t535;
t495 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t554;
t417 = qJDD(3) * pkin(5) + qJD(3) * t494 + t508 + t563 * t526 + (-0.2e1 * qJD(5) - t488) * t553 - t535;
t536 = -m(7) * t417 + mrSges(7,1) * t441 - t442 * mrSges(7,2) + t446 * t480 - t481 * t447;
t531 = -m(6) * t421 + qJDD(3) * mrSges(6,1) - t490 * mrSges(6,3) + qJD(3) * t495 - t536;
t529 = m(5) * t425 + qJDD(3) * mrSges(5,3) + qJD(3) * t497 + t484 * t553 + t531;
t403 = t559 * t490 + m(4) * t429 - qJD(3) * t496 - qJDD(3) * mrSges(4,2) + (t485 - t486) * t553 + t529;
t541 = -t396 * t521 + t524 * t403;
t420 = -pkin(3) * t543 + pkin(4) * t490 + t530;
t399 = m(6) * t420 + t489 * mrSges(6,1) - t490 * mrSges(6,2) + t495 * t554 - t498 * t553 + t401;
t435 = Ifges(7,4) * t481 + Ifges(7,2) * t480 + Ifges(7,6) * t503;
t436 = Ifges(7,1) * t481 + Ifges(7,4) * t480 + Ifges(7,5) * t503;
t534 = mrSges(7,1) * t413 - mrSges(7,2) * t414 + Ifges(7,5) * t442 + Ifges(7,6) * t441 + Ifges(7,3) * t478 + t481 * t435 - t480 * t436;
t533 = qJDD(3) * mrSges(6,2) + qJD(3) * t498 - t539;
t423 = (pkin(3) * qJD(3) + t565) * t554 + t537;
t532 = m(5) * t423 - t489 * mrSges(5,3) - t497 * t554 - t500 * t553 - t399;
t528 = -m(4) * t432 + t490 * mrSges(4,1) - t496 * t554 + t499 * t553 - t532;
t434 = Ifges(7,5) * t481 + Ifges(7,6) * t480 + Ifges(7,3) * t503;
t405 = mrSges(7,2) * t417 - mrSges(7,3) * t413 + Ifges(7,1) * t442 + Ifges(7,4) * t441 + Ifges(7,5) * t478 + t434 * t480 - t435 * t503;
t404 = -mrSges(7,1) * t417 + mrSges(7,3) * t414 + Ifges(7,4) * t442 + Ifges(7,2) * t441 + Ifges(7,6) * t478 - t434 * t481 + t436 * t503;
t400 = -mrSges(6,3) * t489 + t533;
t398 = t484 * t554 + (mrSges(5,2) - mrSges(6,3)) * t489 + t533 + t566;
t397 = -mrSges(5,1) * t490 + t532;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t542 - mrSges(2,2) * t540 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t443 - mrSges(3,2) * t444 + t521 * (mrSges(6,1) * t420 + mrSges(4,2) * t432 + mrSges(5,2) * t426 - mrSges(4,3) * t428 - mrSges(5,3) * t423 - mrSges(6,3) * t422 + pkin(5) * t401 - qJ(4) * t397 - qJ(5) * t400 + t534) + t524 * (pkin(8) * t401 - pkin(3) * t397 + pkin(4) * t399 - mrSges(6,2) * t420 + mrSges(6,3) * t421 - mrSges(5,1) * t423 + mrSges(5,2) * t425 + mrSges(4,3) * t429 - mrSges(4,1) * t432 + t520 * t404 - t523 * t405 - qJ(5) * (-t486 * t553 + t531)) + pkin(2) * t528 + pkin(7) * t541 + pkin(1) * (t517 * (m(3) * t444 - mrSges(3,1) * t527 - qJDD(1) * mrSges(3,2) + t541) + t518 * (m(3) * t443 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t527 + t528)) + (t545 * mrSges(5,1) - t574) * t490 + (-t545 * mrSges(4,2) + t575) * t489 + (t521 * t549 + t524 * t548) * qJDD(3) + (t524 * t546 + t569) * qJD(3); m(3) * t515 + t396 * t524 + t403 * t521; -pkin(3) * t398 - pkin(4) * t400 - pkin(8) * t555 - mrSges(6,1) * t421 + mrSges(6,2) * t422 + mrSges(5,3) * t425 - mrSges(5,1) * t426 + mrSges(4,1) * t428 - mrSges(4,2) * t429 - pkin(5) * t536 - t520 * t405 - t523 * t404 + qJ(4) * t529 + (qJ(4) * mrSges(5,2) + t548) * t490 + t549 * t489 + (Ifges(4,3) + Ifges(5,2) + Ifges(6,3)) * qJDD(3) + (-t569 + (-qJ(4) * t486 - t546) * t524) * qJD(1); t398; t399; t534;];
tauJ  = t1;
