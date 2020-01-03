% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRPPR5
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
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
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRPPR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR5_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR5_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR5_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR5_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:08
% EndTime: 2019-12-31 17:38:10
% DurationCPUTime: 0.94s
% Computational Cost: add. (8558->170), mult. (13304->206), div. (0->0), fcn. (5978->8), ass. (0->76)
t428 = qJD(2) ^ 2;
t422 = sin(pkin(7));
t444 = cos(pkin(7));
t409 = -t444 * g(1) - t422 * g(2);
t419 = -g(3) + qJDD(1);
t425 = sin(qJ(2));
t427 = cos(qJ(2));
t395 = t427 * t409 + t425 * t419;
t434 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t395;
t448 = -pkin(2) - pkin(3);
t389 = t448 * t428 + t434;
t394 = -t425 * t409 + t427 * t419;
t431 = -t428 * qJ(3) + qJDD(3) - t394;
t391 = t448 * qJDD(2) + t431;
t421 = sin(pkin(8));
t423 = cos(pkin(8));
t386 = t423 * t389 + t421 * t391;
t384 = -t428 * pkin(4) - qJDD(2) * pkin(6) + t386;
t408 = t422 * g(1) - t444 * g(2);
t407 = qJDD(4) + t408;
t424 = sin(qJ(5));
t426 = cos(qJ(5));
t381 = -t424 * t384 + t426 * t407;
t382 = t426 * t384 + t424 * t407;
t397 = Ifges(6,6) * qJD(5) + (-Ifges(6,4) * t424 - Ifges(6,2) * t426) * qJD(2);
t398 = Ifges(6,5) * qJD(5) + (-Ifges(6,1) * t424 - Ifges(6,4) * t426) * qJD(2);
t440 = qJD(2) * qJD(5);
t405 = -t424 * qJDD(2) - t426 * t440;
t406 = -t426 * qJDD(2) + t424 * t440;
t451 = mrSges(6,1) * t381 - mrSges(6,2) * t382 + Ifges(6,5) * t405 + Ifges(6,6) * t406 + Ifges(6,3) * qJDD(5) - (t397 * t424 - t398 * t426) * qJD(2);
t404 = (mrSges(6,1) * t426 - mrSges(6,2) * t424) * qJD(2);
t441 = qJD(2) * t426;
t411 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t441;
t442 = qJD(2) * t424;
t378 = m(6) * t381 + qJDD(5) * mrSges(6,1) - t405 * mrSges(6,3) + qJD(5) * t411 + t404 * t442;
t410 = qJD(5) * mrSges(6,1) + mrSges(6,3) * t442;
t379 = m(6) * t382 - qJDD(5) * mrSges(6,2) + t406 * mrSges(6,3) - qJD(5) * t410 - t404 * t441;
t371 = -t424 * t378 + t426 * t379;
t365 = m(5) * t386 - t428 * mrSges(5,1) + qJDD(2) * mrSges(5,2) + t371;
t385 = -t421 * t389 + t423 * t391;
t383 = qJDD(2) * pkin(4) - t428 * pkin(6) - t385;
t380 = -m(6) * t383 + t406 * mrSges(6,1) - t405 * mrSges(6,2) + t410 * t442 - t411 * t441;
t374 = m(5) * t385 - qJDD(2) * mrSges(5,1) - t428 * mrSges(5,2) + t380;
t363 = t421 * t365 + t423 * t374;
t393 = -qJDD(2) * pkin(2) + t431;
t362 = m(4) * t393 - qJDD(2) * mrSges(4,1) - t428 * mrSges(4,3) + t363;
t396 = Ifges(6,3) * qJD(5) + (-Ifges(6,5) * t424 - Ifges(6,6) * t426) * qJD(2);
t372 = -mrSges(6,1) * t383 + mrSges(6,3) * t382 + Ifges(6,4) * t405 + Ifges(6,2) * t406 + Ifges(6,6) * qJDD(5) + qJD(5) * t398 + t396 * t442;
t373 = mrSges(6,2) * t383 - mrSges(6,3) * t381 + Ifges(6,1) * t405 + Ifges(6,4) * t406 + Ifges(6,5) * qJDD(5) - qJD(5) * t397 - t396 * t441;
t392 = -t428 * pkin(2) + t434;
t435 = t423 * t365 - t421 * t374;
t432 = m(4) * t392 + qJDD(2) * mrSges(4,3) + t435;
t450 = (Ifges(3,3) + Ifges(4,2) + Ifges(5,3)) * qJDD(2) + mrSges(3,1) * t394 - mrSges(4,1) * t393 - mrSges(5,1) * t385 - mrSges(3,2) * t395 + mrSges(5,2) * t386 + mrSges(4,3) * t392 - pkin(2) * t362 - pkin(3) * t363 - pkin(4) * t380 - pkin(6) * t371 + qJ(3) * (-t428 * mrSges(4,1) + t432) - t426 * t372 - t424 * t373;
t449 = m(3) + m(4);
t447 = mrSges(3,1) + mrSges(4,1);
t446 = Ifges(4,4) + Ifges(3,5);
t445 = Ifges(3,6) - Ifges(4,6);
t358 = m(3) * t395 - qJDD(2) * mrSges(3,2) - t447 * t428 + t432;
t359 = m(3) * t394 + qJDD(2) * mrSges(3,1) - t428 * mrSges(3,2) - t362;
t436 = t427 * t358 - t425 * t359;
t352 = m(2) * t409 + t436;
t370 = t426 * t378 + t424 * t379;
t369 = m(5) * t407 + t370;
t367 = (m(2) + t449) * t408 + t369;
t443 = t422 * t352 + t444 * t367;
t353 = t425 * t358 + t427 * t359;
t438 = m(2) * t419 + t353;
t437 = t444 * t352 - t422 * t367;
t368 = -m(4) * t408 - t369;
t360 = -mrSges(5,1) * t407 + mrSges(5,3) * t386 + t428 * Ifges(5,5) - Ifges(5,6) * qJDD(2) - pkin(4) * t370 - t451;
t354 = mrSges(5,2) * t407 - mrSges(5,3) * t385 - Ifges(5,5) * qJDD(2) - t428 * Ifges(5,6) - pkin(6) * t370 - t424 * t372 + t426 * t373;
t349 = mrSges(4,2) * t393 - mrSges(3,3) * t394 - qJ(3) * t368 - qJ(4) * t363 + t423 * t354 - t421 * t360 - t445 * t428 + (-mrSges(3,2) + mrSges(4,3)) * t408 + t446 * qJDD(2);
t348 = mrSges(4,2) * t392 + mrSges(3,3) * t395 - pkin(2) * t368 + pkin(3) * t369 - qJ(4) * t435 + t445 * qJDD(2) - t421 * t354 - t423 * t360 + t447 * t408 + t446 * t428;
t347 = -mrSges(2,1) * t419 + mrSges(2,3) * t409 - pkin(1) * t353 - t450;
t346 = mrSges(2,2) * t419 - mrSges(2,3) * t408 - pkin(5) * t353 - t425 * t348 + t427 * t349;
t1 = [-m(1) * g(1) + t437; -m(1) * g(2) + t443; -m(1) * g(3) + t438; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t443 + t444 * t346 - t422 * t347; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t437 + t422 * t346 + t444 * t347; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - mrSges(2,2) * t409 + t425 * t349 + t427 * t348 + pkin(1) * t369 + pkin(5) * t436 + (pkin(1) * t449 + mrSges(2,1)) * t408; t438; t450; t362; t369; t451;];
tauJB = t1;
