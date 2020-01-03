% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRPP2
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
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRPP2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:36
% EndTime: 2019-12-31 20:51:38
% DurationCPUTime: 1.26s
% Computational Cost: add. (12003->223), mult. (15096->269), div. (0->0), fcn. (6511->6), ass. (0->86)
t548 = Ifges(4,1) + Ifges(5,1) + Ifges(6,1);
t538 = Ifges(4,4) - Ifges(5,5) - Ifges(6,4);
t537 = Ifges(4,5) + Ifges(5,4) - Ifges(6,5);
t547 = Ifges(4,2) + Ifges(5,3) + Ifges(6,2);
t536 = Ifges(4,6) - Ifges(5,6) + Ifges(6,6);
t546 = -Ifges(4,3) - Ifges(5,2) - Ifges(6,3);
t545 = 2 * qJD(4);
t506 = qJD(1) + qJD(2);
t504 = t506 ^ 2;
t544 = t504 * pkin(7);
t543 = mrSges(4,3) + mrSges(5,2);
t513 = sin(qJ(3));
t542 = t506 * t513;
t516 = cos(qJ(3));
t541 = t506 * t516;
t515 = sin(qJ(1));
t518 = cos(qJ(1));
t498 = t515 * g(1) - t518 * g(2);
t489 = qJDD(1) * pkin(1) + t498;
t499 = -t518 * g(1) - t515 * g(2);
t520 = qJD(1) ^ 2;
t490 = -t520 * pkin(1) + t499;
t514 = sin(qJ(2));
t517 = cos(qJ(2));
t450 = t514 * t489 + t517 * t490;
t505 = qJDD(1) + qJDD(2);
t448 = -t504 * pkin(2) + t505 * pkin(7) + t450;
t444 = -t513 * g(3) + t516 * t448;
t477 = (mrSges(6,1) * t516 + mrSges(6,2) * t513) * t506;
t478 = (-mrSges(4,1) * t516 + mrSges(4,2) * t513) * t506;
t480 = -qJD(3) * t542 + t516 * t505;
t493 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t542;
t475 = (-pkin(3) * t516 - qJ(4) * t513) * t506;
t519 = qJD(3) ^ 2;
t441 = -t519 * pkin(3) + qJDD(3) * qJ(4) + qJD(3) * t545 + t475 * t541 + t444;
t476 = (-mrSges(5,1) * t516 - mrSges(5,3) * t513) * t506;
t494 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t542;
t491 = -qJD(3) * pkin(4) - qJ(5) * t542;
t512 = t516 ^ 2;
t535 = -0.2e1 * qJD(5) * t506;
t437 = -t512 * t504 * pkin(4) - t480 * qJ(5) + qJD(3) * t491 + t516 * t535 + t441;
t492 = -qJD(3) * mrSges(6,1) - mrSges(6,3) * t542;
t527 = m(6) * t437 + qJDD(3) * mrSges(6,2) - t480 * mrSges(6,3) + qJD(3) * t492;
t523 = m(5) * t441 + qJDD(3) * mrSges(5,3) + qJD(3) * t494 + t476 * t541 + t527;
t430 = m(4) * t444 - qJDD(3) * mrSges(4,2) - qJD(3) * t493 + (-t477 + t478) * t541 + t543 * t480 + t523;
t443 = -t516 * g(3) - t513 * t448;
t539 = qJD(3) * t516;
t531 = t506 * t539;
t479 = t513 * t505 + t531;
t496 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t541;
t442 = -qJDD(3) * pkin(3) - t519 * qJ(4) + t475 * t542 + qJDD(4) - t443;
t438 = t513 * t535 + (-t479 + t531) * qJ(5) + (-t504 * t513 * t516 - qJDD(3)) * pkin(4) + t442;
t495 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t541;
t433 = m(6) * t438 - qJDD(3) * mrSges(6,1) - t479 * mrSges(6,3) - qJD(3) * t495 - t477 * t542;
t497 = mrSges(5,2) * t541 + qJD(3) * mrSges(5,3);
t522 = -m(5) * t442 + qJDD(3) * mrSges(5,1) + qJD(3) * t497 - t433;
t431 = m(4) * t443 + qJDD(3) * mrSges(4,1) + qJD(3) * t496 + (-t476 - t478) * t542 - t543 * t479 + t522;
t528 = t516 * t430 - t513 * t431;
t423 = m(3) * t450 - t504 * mrSges(3,1) - t505 * mrSges(3,2) + t528;
t449 = t517 * t489 - t514 * t490;
t526 = t505 * pkin(2) + t449;
t524 = -t479 * qJ(4) - t526;
t439 = -t480 * pkin(3) - t544 + (-0.2e1 * qJD(4) * t513 + (pkin(3) * t513 - qJ(4) * t516) * qJD(3)) * t506 + t524;
t435 = qJDD(5) + (-qJ(5) * t512 + pkin(7)) * t504 + (pkin(3) + pkin(4)) * t480 + (qJ(4) * t539 + (-pkin(3) * qJD(3) + t491 + t545) * t513) * t506 - t524;
t525 = -m(6) * t435 - t480 * mrSges(6,1) - t479 * mrSges(6,2) - t492 * t542 - t495 * t541;
t432 = m(5) * t439 - t480 * mrSges(5,1) - t479 * mrSges(5,3) - t494 * t542 - t497 * t541 + t525;
t447 = -t526 - t544;
t521 = -m(4) * t447 + t480 * mrSges(4,1) - t479 * mrSges(4,2) - t493 * t542 + t496 * t541 - t432;
t426 = m(3) * t449 + t505 * mrSges(3,1) - t504 * mrSges(3,2) + t521;
t419 = t514 * t423 + t517 * t426;
t417 = m(2) * t498 + qJDD(1) * mrSges(2,1) - t520 * mrSges(2,2) + t419;
t529 = t517 * t423 - t514 * t426;
t418 = m(2) * t499 - t520 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t529;
t540 = t518 * t417 + t515 * t418;
t424 = t513 * t430 + t516 * t431;
t534 = (-t537 * t513 - t536 * t516) * t506 + t546 * qJD(3);
t533 = (-t538 * t513 - t547 * t516) * t506 - t536 * qJD(3);
t532 = (t548 * t513 + t538 * t516) * t506 + t537 * qJD(3);
t530 = -t515 * t417 + t518 * t418;
t420 = mrSges(4,2) * t447 + mrSges(5,2) * t442 + mrSges(6,2) * t435 - mrSges(4,3) * t443 - mrSges(5,3) * t439 - mrSges(6,3) * t438 - qJ(4) * t432 - qJ(5) * t433 + t533 * qJD(3) + t537 * qJDD(3) + t548 * t479 + t538 * t480 - t534 * t541;
t413 = -mrSges(4,1) * t447 + mrSges(4,3) * t444 - mrSges(5,1) * t439 + mrSges(5,2) * t441 + mrSges(6,1) * t435 - mrSges(6,3) * t437 - pkin(4) * t525 - qJ(5) * t527 - pkin(3) * t432 + (qJ(5) * t477 * t516 + t534 * t513) * t506 + t547 * t480 + t538 * t479 + t536 * qJDD(3) + t532 * qJD(3);
t412 = -pkin(3) * t522 - qJ(4) * t523 + t504 * Ifges(3,5) + Ifges(3,6) * t505 + mrSges(3,3) * t450 + mrSges(6,1) * t438 - mrSges(5,3) * t441 + mrSges(5,1) * t442 - mrSges(4,1) * t443 + mrSges(4,2) * t444 + pkin(4) * t433 - mrSges(6,2) * t437 - pkin(2) * t424 + mrSges(3,1) * g(3) + (-qJ(4) * mrSges(5,2) - t536) * t480 + (pkin(3) * mrSges(5,2) - t537) * t479 + t546 * qJDD(3) + ((qJ(4) * t477 + t532) * t516 + (pkin(3) * t476 + t533) * t513) * t506;
t411 = -mrSges(3,2) * g(3) - mrSges(3,3) * t449 + Ifges(3,5) * t505 - t504 * Ifges(3,6) - pkin(7) * t424 - t513 * t413 + t516 * t420;
t410 = -mrSges(2,2) * g(3) - mrSges(2,3) * t498 + Ifges(2,5) * qJDD(1) - t520 * Ifges(2,6) - pkin(6) * t419 + t517 * t411 - t514 * t412;
t409 = Ifges(2,6) * qJDD(1) + t520 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t499 + t514 * t411 + t517 * t412 - pkin(1) * (-m(3) * g(3) + t424) + pkin(6) * t529;
t1 = [-m(1) * g(1) + t530; -m(1) * g(2) + t540; (-m(1) - m(2) - m(3)) * g(3) + t424; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t540 - t515 * t409 + t518 * t410; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t530 + t518 * t409 + t515 * t410; -mrSges(1,1) * g(2) + mrSges(2,1) * t498 + mrSges(3,1) * t449 + mrSges(1,2) * g(1) - mrSges(2,2) * t499 - mrSges(3,2) * t450 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * t505 + pkin(1) * t419 + pkin(2) * t521 + pkin(7) * t528 + t516 * t413 + t513 * t420;];
tauB = t1;
