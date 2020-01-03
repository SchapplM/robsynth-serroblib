% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPRP3
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPRP3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP3_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP3_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP3_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP3_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:50:49
% EndTime: 2019-12-31 17:50:50
% DurationCPUTime: 0.94s
% Computational Cost: add. (6147->201), mult. (10956->232), div. (0->0), fcn. (4655->6), ass. (0->82)
t459 = Ifges(5,1) + Ifges(6,1);
t451 = Ifges(5,4) + Ifges(6,4);
t449 = Ifges(5,5) + Ifges(6,5);
t458 = -Ifges(5,2) - Ifges(6,2);
t447 = Ifges(5,6) + Ifges(6,6);
t457 = (Ifges(5,3) + Ifges(6,3));
t420 = sin(qJ(1));
t422 = cos(qJ(1));
t403 = t420 * g(1) - t422 * g(2);
t391 = qJDD(1) * pkin(1) + t403;
t404 = -t422 * g(1) - t420 * g(2);
t423 = qJD(1) ^ 2;
t394 = -t423 * pkin(1) + t404;
t417 = sin(pkin(7));
t418 = cos(pkin(7));
t374 = t417 * t391 + t418 * t394;
t429 = qJDD(1) * qJ(3) + (2 * qJD(3) * qJD(1)) + t374;
t455 = -pkin(2) - pkin(6);
t369 = t455 * t423 + t429;
t419 = sin(qJ(4));
t421 = cos(qJ(4));
t439 = qJD(1) * qJD(4);
t395 = -t419 * qJDD(1) - t421 * t439;
t396 = t421 * qJDD(1) - t419 * t439;
t441 = qJD(1) * t419;
t399 = -(qJD(4) * mrSges(5,2)) - mrSges(5,3) * t441;
t440 = qJD(1) * t421;
t402 = (qJD(4) * mrSges(5,1)) - mrSges(5,3) * t440;
t400 = (qJD(4) * pkin(4)) - qJ(5) * t440;
t413 = t419 ^ 2;
t363 = t400 * t440 - t395 * pkin(4) + qJDD(5) + (-qJ(5) * t413 + t455) * t423 + t429;
t398 = -(qJD(4) * mrSges(6,2)) - mrSges(6,3) * t441;
t401 = (qJD(4) * mrSges(6,1)) - mrSges(6,3) * t440;
t430 = m(6) * t363 + t396 * mrSges(6,2) + t398 * t441 + t401 * t440;
t456 = -m(5) * t369 - t396 * mrSges(5,2) + (mrSges(5,1) + mrSges(6,1)) * t395 - t399 * t441 - t402 * t440 - t430;
t454 = pkin(4) * t423;
t453 = mrSges(3,1) - mrSges(4,2);
t450 = Ifges(3,5) - Ifges(4,4);
t448 = Ifges(3,6) - Ifges(4,5);
t373 = t418 * t391 - t417 * t394;
t427 = -t423 * qJ(3) + qJDD(3) - t373;
t370 = t455 * qJDD(1) + t427;
t367 = t421 * t370;
t414 = -g(3) + qJDD(2);
t364 = -t419 * t414 + t367;
t392 = (mrSges(6,1) * t419 + mrSges(6,2) * t421) * qJD(1);
t431 = qJD(1) * (-t392 - (mrSges(5,1) * t419 + mrSges(5,2) * t421) * qJD(1));
t437 = -2 * qJD(1) * qJD(5);
t360 = t421 * t437 + (qJDD(4) * pkin(4)) - t396 * qJ(5) + t367 + (-qJ(5) * t439 - t421 * t454 - t414) * t419;
t436 = m(6) * t360 + (qJDD(4) * mrSges(6,1)) + qJD(4) * t398;
t355 = m(5) * t364 + (qJDD(4) * mrSges(5,1)) + qJD(4) * t399 + (-mrSges(5,3) - mrSges(6,3)) * t396 + t421 * t431 + t436;
t365 = t419 * t370 + t421 * t414;
t361 = t395 * qJ(5) - qJD(4) * t400 - t413 * t454 + t419 * t437 + t365;
t445 = m(6) * t361 + t395 * mrSges(6,3);
t356 = m(5) * t365 + t395 * mrSges(5,3) + ((-mrSges(5,2) - mrSges(6,2)) * qJDD(4)) + (-t401 - t402) * qJD(4) + t419 * t431 + t445;
t350 = t421 * t355 + t419 * t356;
t372 = -qJDD(1) * pkin(2) + t427;
t425 = -m(4) * t372 + t423 * mrSges(4,3) - t350;
t347 = m(3) * t373 - t423 * mrSges(3,2) + t453 * qJDD(1) + t425;
t371 = t423 * pkin(2) - t429;
t424 = -m(4) * t371 + t423 * mrSges(4,2) + qJDD(1) * mrSges(4,3) - t456;
t353 = m(3) * t374 - t423 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t424;
t344 = t418 * t347 + t417 * t353;
t342 = m(2) * t403 + qJDD(1) * mrSges(2,1) - t423 * mrSges(2,2) + t344;
t433 = -t417 * t347 + t418 * t353;
t343 = m(2) * t404 - t423 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t433;
t446 = t422 * t342 + t420 * t343;
t444 = -(t457 * qJD(4)) + (t447 * t419 - t449 * t421) * qJD(1);
t443 = t447 * qJD(4) + (t458 * t419 + t451 * t421) * qJD(1);
t442 = t449 * qJD(4) + (-t451 * t419 + t459 * t421) * qJD(1);
t434 = -t420 * t342 + t422 * t343;
t432 = -t419 * t355 + t421 * t356;
t349 = m(4) * t414 + t432;
t428 = m(3) * t414 + t349;
t357 = -t396 * mrSges(6,3) - t392 * t440 + t436;
t348 = mrSges(5,2) * t369 + mrSges(6,2) * t363 - mrSges(5,3) * t364 - mrSges(6,3) * t360 - qJ(5) * t357 - t443 * qJD(4) + t449 * qJDD(4) + t451 * t395 + t459 * t396 + t444 * t441;
t345 = -mrSges(5,1) * t369 + mrSges(5,3) * t365 - mrSges(6,1) * t363 + mrSges(6,3) * t361 - pkin(4) * t430 + qJ(5) * t445 + t451 * t396 + ((pkin(4) * mrSges(6,1)) - t458) * t395 + (-qJ(5) * mrSges(6,2) + t447) * qJDD(4) + (-qJ(5) * t401 + t442) * qJD(4) + (-qJ(5) * t392 * t419 + t444 * t421) * qJD(1);
t338 = mrSges(4,1) * t372 + mrSges(5,1) * t364 + mrSges(6,1) * t360 - mrSges(5,2) * t365 - mrSges(6,2) * t361 - mrSges(3,3) * t373 + pkin(3) * t350 + pkin(4) * t357 - qJ(3) * t349 - t448 * t423 + (mrSges(3,2) - mrSges(4,3)) * t414 + t449 * t396 + t447 * t395 + (t457 * qJDD(4)) + t450 * qJDD(1) + (t442 * t419 + t443 * t421) * qJD(1);
t337 = -mrSges(4,1) * t371 + mrSges(3,3) * t374 - pkin(2) * t349 - pkin(3) * t456 - pkin(6) * t432 + t448 * qJDD(1) - t421 * t345 - t419 * t348 - t453 * t414 + t450 * t423;
t336 = -mrSges(2,2) * g(3) - mrSges(2,3) * t403 + Ifges(2,5) * qJDD(1) - t423 * Ifges(2,6) - qJ(2) * t344 - t417 * t337 + t418 * t338;
t335 = mrSges(2,1) * g(3) + mrSges(2,3) * t404 + t423 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t428 + qJ(2) * t433 + t418 * t337 + t417 * t338;
t1 = [-m(1) * g(1) + t434; -m(1) * g(2) + t446; (-m(1) - m(2)) * g(3) + t428; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t446 - t420 * t335 + t422 * t336; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t434 + t422 * t335 + t420 * t336; pkin(1) * t344 + mrSges(2,1) * t403 - mrSges(2,2) * t404 + pkin(2) * t425 + qJ(3) * t424 + mrSges(3,1) * t373 - mrSges(3,2) * t374 + t421 * t348 - t419 * t345 - pkin(6) * t350 + mrSges(4,2) * t372 - mrSges(4,3) * t371 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(2) * mrSges(4,2) + Ifges(4,1) + Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB = t1;
