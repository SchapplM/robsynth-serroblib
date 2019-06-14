% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPPRRP5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2019-05-05 14:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPPRRP5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:57:31
% EndTime: 2019-05-05 14:57:33
% DurationCPUTime: 1.03s
% Computational Cost: add. (4058->209), mult. (7509->243), div. (0->0), fcn. (3836->6), ass. (0->82)
t448 = Ifges(6,1) + Ifges(7,1);
t438 = Ifges(6,4) + Ifges(7,4);
t437 = Ifges(6,5) + Ifges(7,5);
t447 = Ifges(6,2) + Ifges(7,2);
t436 = Ifges(6,6) + Ifges(7,6);
t403 = sin(qJ(5));
t406 = cos(qJ(5));
t407 = cos(qJ(4));
t428 = qJD(1) * t407;
t385 = qJD(4) * t406 - t403 * t428;
t404 = sin(qJ(4));
t427 = qJD(1) * qJD(4);
t422 = t404 * t427;
t390 = qJDD(1) * t407 - t422;
t359 = qJD(5) * t385 + qJDD(4) * t403 + t390 * t406;
t386 = qJD(4) * t403 + t406 * t428;
t364 = -mrSges(7,1) * t385 + mrSges(7,2) * t386;
t410 = qJD(1) ^ 2;
t405 = sin(qJ(1));
t408 = cos(qJ(1));
t430 = t405 * g(1) - t408 * g(2);
t377 = -qJDD(1) * pkin(1) - t410 * qJ(2) + qJDD(2) - t430;
t369 = -qJDD(1) * qJ(3) - 0.2e1 * qJD(3) * qJD(1) + t377;
t366 = -t410 * pkin(7) - t369;
t421 = t407 * t427;
t389 = -qJDD(1) * t404 - t421;
t342 = (-t390 + t422) * pkin(8) + (-t389 + t421) * pkin(4) + t366;
t418 = -t408 * g(1) - t405 * g(2);
t440 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t418;
t370 = qJDD(3) + (-pkin(1) - qJ(3)) * t410 + t440;
t367 = -qJDD(1) * pkin(7) + t370;
t361 = -g(3) * t407 + t404 * t367;
t388 = (pkin(4) * t404 - pkin(8) * t407) * qJD(1);
t409 = qJD(4) ^ 2;
t429 = qJD(1) * t404;
t345 = -pkin(4) * t409 + qJDD(4) * pkin(8) - t388 * t429 + t361;
t338 = t406 * t342 - t403 * t345;
t384 = qJDD(5) - t389;
t393 = qJD(5) + t429;
t334 = -0.2e1 * qJD(6) * t386 + (t385 * t393 - t359) * qJ(6) + (t385 * t386 + t384) * pkin(5) + t338;
t371 = -mrSges(7,2) * t393 + mrSges(7,3) * t385;
t424 = m(7) * t334 + t384 * mrSges(7,1) + t393 * t371;
t331 = -t359 * mrSges(7,3) - t386 * t364 + t424;
t339 = t403 * t342 + t406 * t345;
t358 = -qJD(5) * t386 + qJDD(4) * t406 - t390 * t403;
t373 = pkin(5) * t393 - qJ(6) * t386;
t383 = t385 ^ 2;
t336 = -pkin(5) * t383 + qJ(6) * t358 + 0.2e1 * qJD(6) * t385 - t373 * t393 + t339;
t433 = -t447 * t385 - t438 * t386 - t436 * t393;
t442 = t438 * t385 + t448 * t386 + t437 * t393;
t444 = Ifges(6,3) + Ifges(7,3);
t446 = mrSges(6,1) * t338 + mrSges(7,1) * t334 - mrSges(6,2) * t339 - mrSges(7,2) * t336 + pkin(5) * t331 + t436 * t358 + t437 * t359 + t444 * t384 - t442 * t385 - t433 * t386;
t365 = -mrSges(6,1) * t385 + mrSges(6,2) * t386;
t372 = -mrSges(6,2) * t393 + mrSges(6,3) * t385;
t326 = m(6) * t338 + t384 * mrSges(6,1) + t393 * t372 + (-t364 - t365) * t386 + (-mrSges(6,3) - mrSges(7,3)) * t359 + t424;
t423 = m(7) * t336 + t358 * mrSges(7,3) + t385 * t364;
t374 = mrSges(7,1) * t393 - mrSges(7,3) * t386;
t431 = -mrSges(6,1) * t393 + mrSges(6,3) * t386 - t374;
t439 = -mrSges(6,2) - mrSges(7,2);
t329 = m(6) * t339 + t358 * mrSges(6,3) + t385 * t365 + t439 * t384 + t431 * t393 + t423;
t324 = t406 * t326 + t403 * t329;
t391 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t429;
t392 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t428;
t441 = m(4) * t369 - m(5) * t366 + t389 * mrSges(5,1) - t390 * mrSges(5,2) - t324 + (-t391 * t404 - t392 * t407) * qJD(1);
t360 = g(3) * t404 + t367 * t407;
t387 = (mrSges(5,1) * t404 + mrSges(5,2) * t407) * qJD(1);
t344 = -qJDD(4) * pkin(4) - pkin(8) * t409 + t388 * t428 - t360;
t337 = -pkin(5) * t358 - qJ(6) * t383 + t373 * t386 + qJDD(6) + t344;
t419 = -m(7) * t337 + t358 * mrSges(7,1) + t385 * t371;
t411 = -m(6) * t344 + t358 * mrSges(6,1) + t439 * t359 + t385 * t372 + t431 * t386 + t419;
t420 = -t326 * t403 + t406 * t329;
t435 = t404 * (m(5) * t361 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t389 - qJD(4) * t392 - t387 * t429 + t420) + t407 * (m(5) * t360 + qJDD(4) * mrSges(5,1) - t390 * mrSges(5,3) + qJD(4) * t391 - t387 * t428 + t411);
t434 = -t436 * t385 - t437 * t386 - t444 * t393;
t415 = m(4) * t370 + qJDD(1) * mrSges(4,2) - mrSges(4,3) * t410 + t435;
t381 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t407 - Ifges(5,4) * t404) * qJD(1);
t380 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t407 - Ifges(5,2) * t404) * qJD(1);
t376 = pkin(1) * t410 - t440;
t332 = t359 * mrSges(7,2) + t386 * t374 - t419;
t323 = mrSges(6,2) * t344 + mrSges(7,2) * t337 - mrSges(6,3) * t338 - mrSges(7,3) * t334 - qJ(6) * t331 + t438 * t358 + t448 * t359 + t437 * t384 - t434 * t385 + t433 * t393;
t321 = m(3) * t377 + (-mrSges(4,2) - mrSges(3,3)) * t410 + (mrSges(3,2) - mrSges(4,3)) * qJDD(1) + t441;
t320 = -mrSges(6,1) * t344 + mrSges(6,3) * t339 - mrSges(7,1) * t337 + mrSges(7,3) * t336 - pkin(5) * t332 + qJ(6) * t423 + (-qJ(6) * t374 + t442) * t393 + t434 * t386 + (-mrSges(7,2) * qJ(6) + t436) * t384 + t438 * t359 + t447 * t358;
t1 = [-pkin(1) * t321 + mrSges(2,1) * t430 - mrSges(2,2) * t418 + qJ(2) * (-m(3) * t376 + mrSges(3,2) * t410 + t415) + mrSges(3,2) * t377 - mrSges(3,3) * t376 - t404 * (-mrSges(5,1) * t366 + mrSges(5,3) * t361 + Ifges(5,4) * t390 + Ifges(5,2) * t389 + Ifges(5,6) * qJDD(4) - pkin(4) * t324 + qJD(4) * t381 - t446) - pkin(7) * t435 + mrSges(4,2) * t370 - mrSges(4,3) * t369 + t407 * (mrSges(5,2) * t366 - mrSges(5,3) * t360 + Ifges(5,1) * t390 + Ifges(5,4) * t389 + Ifges(5,5) * qJDD(4) - pkin(8) * t324 - qJD(4) * t380 - t320 * t403 + t323 * t406) + (mrSges(3,3) * qJ(2) + Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1) + (t410 * mrSges(4,2) + qJDD(1) * mrSges(4,3) - t441) * qJ(3); t321; t415; Ifges(5,5) * t390 + Ifges(5,6) * t389 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t360 - mrSges(5,2) * t361 + t403 * t323 + t406 * t320 + pkin(4) * t411 + pkin(8) * t420 + (t380 * t407 + t381 * t404) * qJD(1); t446; t332;];
tauJ  = t1;
