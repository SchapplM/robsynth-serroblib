% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRPP5
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
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
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRPP5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP5_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP5_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP5_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP5_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:14
% EndTime: 2019-12-31 18:16:15
% DurationCPUTime: 1.03s
% Computational Cost: add. (3960->225), mult. (7598->255), div. (0->0), fcn. (2878->4), ass. (0->87)
t544 = -2 * qJD(1);
t543 = Ifges(4,1) + Ifges(5,1) + Ifges(6,1);
t517 = -Ifges(4,4) + Ifges(5,5) + Ifges(6,4);
t516 = Ifges(4,5) + Ifges(5,4) - Ifges(6,5);
t542 = Ifges(4,2) + Ifges(5,3) + Ifges(6,2);
t515 = Ifges(4,6) - Ifges(5,6) + Ifges(6,6);
t541 = (Ifges(4,3) + Ifges(5,2) + Ifges(6,3));
t492 = sin(qJ(1));
t494 = cos(qJ(1));
t473 = -t494 * g(1) - t492 * g(2);
t496 = qJD(1) ^ 2;
t432 = (t496 * pkin(1)) - qJDD(1) * qJ(2) + (qJD(2) * t544) - t473;
t430 = -(t496 * pkin(6)) - t432;
t491 = sin(qJ(3));
t493 = cos(qJ(3));
t519 = qJD(1) * qJD(3);
t459 = t491 * qJDD(1) + t493 * t519;
t460 = t493 * qJDD(1) - t491 * t519;
t526 = t460 * qJ(4);
t527 = qJ(4) * t491;
t422 = t459 * pkin(3) - t526 + (-0.2e1 * qJD(4) * t493 + (pkin(3) * t493 + t527) * qJD(3)) * qJD(1) + t430;
t521 = qJD(1) * t491;
t471 = -mrSges(5,2) * t521 + (qJD(3) * mrSges(5,3));
t520 = qJD(1) * t493;
t467 = -(qJD(3) * pkin(4)) - qJ(5) * t520;
t487 = t491 ^ 2;
t535 = -pkin(3) - pkin(4);
t536 = 0.2e1 * qJD(4);
t417 = t526 + qJDD(5) + (-qJ(5) * t487 + pkin(6)) * t496 + t535 * t459 + (-qJD(3) * t527 + (-pkin(3) * qJD(3) + t467 + t536) * t493) * qJD(1) + t432;
t465 = (qJD(3) * mrSges(6,2)) + mrSges(6,3) * t521;
t506 = -m(6) * t417 + t459 * mrSges(6,1) + t465 * t521;
t540 = -m(5) * t422 - t459 * mrSges(5,1) + (mrSges(6,2) + mrSges(5,3)) * t460 - t471 * t521 - t506;
t468 = -qJD(3) * mrSges(6,1) - mrSges(6,3) * t520;
t470 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t520;
t508 = qJD(1) * (t468 + t470);
t538 = -m(3) * t432 + (t496 * mrSges(3,2)) + qJDD(1) * mrSges(3,3);
t466 = -(qJD(3) * mrSges(4,2)) - mrSges(4,3) * t521;
t469 = (qJD(3) * mrSges(4,1)) - mrSges(4,3) * t520;
t537 = -m(4) * t430 - t459 * mrSges(4,1) - t460 * mrSges(4,2) - t466 * t521 - t469 * t520;
t533 = pkin(4) * t496;
t532 = mrSges(2,1) - mrSges(3,2);
t531 = -mrSges(4,3) - mrSges(5,2);
t530 = Ifges(2,5) - Ifges(3,4);
t529 = (-Ifges(2,6) + Ifges(3,5));
t472 = t492 * g(1) - t494 * g(2);
t502 = -t496 * qJ(2) + qJDD(2) - t472;
t431 = (-pkin(1) - pkin(6)) * qJDD(1) + t502;
t427 = -t493 * g(3) + t491 * t431;
t455 = (pkin(3) * t491 - qJ(4) * t493) * qJD(1);
t495 = qJD(3) ^ 2;
t500 = -t495 * pkin(3) + qJDD(3) * qJ(4) + qJD(3) * t536 + t427;
t424 = -t455 * t521 + t500;
t419 = -t487 * t533 + t459 * qJ(5) + qJD(3) * t467 + ((2 * qJD(5)) - t455) * t521 + t500;
t457 = (-mrSges(6,1) * t491 + mrSges(6,2) * t493) * qJD(1);
t504 = m(6) * t419 + qJDD(3) * mrSges(6,2) + t459 * mrSges(6,3) + qJD(3) * t468 + t457 * t521;
t499 = m(5) * t424 + qJDD(3) * mrSges(5,3) + qJD(3) * t470 + t504;
t456 = (mrSges(5,1) * t491 - mrSges(5,3) * t493) * qJD(1);
t509 = qJD(1) * (-t456 - (mrSges(4,1) * t491 + mrSges(4,2) * t493) * qJD(1));
t413 = m(4) * t427 - qJDD(3) * mrSges(4,2) - qJD(3) * t469 + t531 * t459 + t491 * t509 + t499;
t525 = t491 * t413;
t426 = t491 * g(3) + t493 * t431;
t505 = -t495 * qJ(4) + t455 * t520 + qJDD(4);
t420 = -t460 * qJ(5) + ((qJD(5) * t544) - t431) * t493 + t535 * qJDD(3) + (-qJ(5) * t519 + t493 * t533 - g(3)) * t491 + t505;
t416 = m(6) * t420 - qJDD(3) * mrSges(6,1) - t460 * mrSges(6,3) - qJD(3) * t465 - t457 * t520;
t425 = -qJDD(3) * pkin(3) - t426 + t505;
t498 = -m(5) * t425 + qJDD(3) * mrSges(5,1) + qJD(3) * t471 - t416;
t414 = m(4) * t426 + qJDD(3) * mrSges(4,1) + qJD(3) * t466 + t531 * t460 + t493 * t509 + t498;
t524 = t493 * t414;
t433 = -qJDD(1) * pkin(1) + t502;
t503 = -m(3) * t433 + t496 * mrSges(3,3) - t525;
t406 = m(2) * t472 - (t496 * mrSges(2,2)) + t532 * qJDD(1) + t503 - t524;
t415 = -t493 * t508 - t540;
t411 = m(2) * t473 - (t496 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t415 - t537 + t538;
t523 = t494 * t406 + t492 * t411;
t514 = -(t541 * qJD(3)) + (t515 * t491 - t516 * t493) * qJD(1);
t513 = -t515 * qJD(3) + (t542 * t491 + t517 * t493) * qJD(1);
t512 = t516 * qJD(3) + (t517 * t491 + t543 * t493) * qJD(1);
t511 = -t492 * t406 + t494 * t411;
t510 = t493 * t413 - t491 * t414;
t497 = t537 + t540;
t408 = t524 + t525;
t407 = -m(3) * g(3) + t510;
t404 = mrSges(4,2) * t430 + mrSges(5,2) * t425 + mrSges(6,2) * t417 - mrSges(4,3) * t426 - mrSges(5,3) * t422 - mrSges(6,3) * t420 - qJ(4) * t415 - qJ(5) * t416 + t513 * qJD(3) + t516 * qJDD(3) + t517 * t459 + t543 * t460 + t514 * t521;
t403 = -mrSges(4,1) * t430 + mrSges(4,3) * t427 - mrSges(5,1) * t422 + mrSges(5,2) * t424 + mrSges(6,1) * t417 - mrSges(6,3) * t419 - pkin(4) * t506 - qJ(5) * t504 - pkin(3) * t415 + ((pkin(4) * mrSges(6,2)) - t517) * t460 - t542 * t459 + t515 * qJDD(3) + t512 * qJD(3) + (pkin(4) * t468 + t514) * t520;
t402 = mrSges(4,1) * t426 - mrSges(4,2) * t427 + mrSges(3,1) * t433 - pkin(4) * t416 + mrSges(6,2) * t419 - mrSges(6,1) * t420 + mrSges(5,3) * t424 - mrSges(5,1) * t425 + pkin(2) * t408 - qJ(2) * t407 + qJ(4) * t499 + pkin(3) * t498 - mrSges(2,3) * t472 + (t529 * t496) + t530 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (-pkin(3) * mrSges(5,2) + t516) * t460 + (-qJ(4) * mrSges(5,2) - t515) * t459 + t541 * qJDD(3) + ((-pkin(3) * t456 - t513) * t493 + (-qJ(4) * t456 + t512) * t491) * qJD(1);
t401 = mrSges(2,3) * t473 - mrSges(3,1) * t432 - t491 * t404 - pkin(2) * t497 - pkin(6) * t510 - pkin(1) * t407 + t530 * t496 + (-pkin(2) * t508 - t403) * t493 - t529 * qJDD(1) + t532 * g(3);
t1 = [-m(1) * g(1) + t511; -m(1) * g(2) + t523; (-m(1) - m(2) - m(3)) * g(3) + t510; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t523 - t492 * t401 + t494 * t402; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t511 + t494 * t401 + t492 * t402; qJ(2) * (-t497 + t538) + mrSges(2,1) * t472 - mrSges(2,2) * t473 + pkin(1) * t503 + mrSges(3,2) * t433 - mrSges(3,3) * t432 - t491 * t403 - pkin(6) * t408 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * t414 - qJ(2) * t508 + t404) * t493 + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB = t1;
