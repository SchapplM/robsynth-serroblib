% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPPRRP6
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
% Datum: 2019-05-05 15:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPPRRP6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP6_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP6_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP6_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP6_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP6_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:00:18
% EndTime: 2019-05-05 15:00:20
% DurationCPUTime: 1.10s
% Computational Cost: add. (3993->207), mult. (7311->243), div. (0->0), fcn. (3687->6), ass. (0->81)
t447 = Ifges(6,1) + Ifges(7,1);
t435 = Ifges(6,4) - Ifges(7,5);
t443 = -Ifges(6,5) - Ifges(7,4);
t446 = Ifges(6,2) + Ifges(7,3);
t433 = Ifges(6,6) - Ifges(7,6);
t402 = sin(qJ(5));
t405 = cos(qJ(4));
t425 = t405 * qJD(1);
t437 = cos(qJ(5));
t381 = -t437 * qJD(4) + t402 * t425;
t403 = sin(qJ(4));
t424 = qJD(1) * qJD(4);
t420 = t403 * t424;
t386 = t405 * qJDD(1) - t420;
t354 = -t381 * qJD(5) + t402 * qJDD(4) + t437 * t386;
t382 = t402 * qJD(4) + t437 * t425;
t361 = t381 * mrSges(7,1) - t382 * mrSges(7,3);
t408 = qJD(1) ^ 2;
t404 = sin(qJ(1));
t406 = cos(qJ(1));
t427 = t404 * g(1) - t406 * g(2);
t373 = -qJDD(1) * pkin(1) - t408 * qJ(2) + qJDD(2) - t427;
t366 = -qJDD(1) * qJ(3) - 0.2e1 * qJD(3) * qJD(1) + t373;
t363 = -t408 * pkin(7) - t366;
t419 = t405 * t424;
t385 = -t403 * qJDD(1) - t419;
t339 = (-t386 + t420) * pkin(8) + (-t385 + t419) * pkin(4) + t363;
t416 = -t406 * g(1) - t404 * g(2);
t438 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t416;
t367 = qJDD(3) + (-pkin(1) - qJ(3)) * t408 + t438;
t364 = -qJDD(1) * pkin(7) + t367;
t356 = -t405 * g(3) + t403 * t364;
t384 = (t403 * pkin(4) - t405 * pkin(8)) * qJD(1);
t407 = qJD(4) ^ 2;
t426 = t403 * qJD(1);
t342 = -t407 * pkin(4) + qJDD(4) * pkin(8) - t384 * t426 + t356;
t336 = t437 * t339 - t402 * t342;
t360 = t381 * pkin(5) - t382 * qJ(6);
t380 = qJDD(5) - t385;
t390 = qJD(5) + t426;
t389 = t390 ^ 2;
t334 = -t380 * pkin(5) - t389 * qJ(6) + t382 * t360 + qJDD(6) - t336;
t371 = -t381 * mrSges(7,2) + t390 * mrSges(7,3);
t417 = -m(7) * t334 + t380 * mrSges(7,1) + t390 * t371;
t331 = t354 * mrSges(7,2) + t382 * t361 - t417;
t337 = t402 * t339 + t437 * t342;
t333 = -t389 * pkin(5) + t380 * qJ(6) + 0.2e1 * qJD(6) * t390 - t381 * t360 + t337;
t353 = t382 * qJD(5) - t437 * qJDD(4) + t402 * t386;
t370 = -t390 * mrSges(7,1) + t382 * mrSges(7,2);
t421 = m(7) * t333 + t380 * mrSges(7,3) + t390 * t370;
t430 = t446 * t381 - t435 * t382 - t433 * t390;
t440 = t435 * t381 - t447 * t382 + t443 * t390;
t442 = -Ifges(6,3) - Ifges(7,2);
t445 = -t443 * t354 - t440 * t381 - t433 * t353 - t442 * t380 + mrSges(6,1) * t336 - mrSges(7,1) * t334 - mrSges(6,2) * t337 + mrSges(7,3) * t333 - pkin(5) * t331 + qJ(6) * (-t353 * mrSges(7,2) - t381 * t361 + t421) - t430 * t382;
t369 = t390 * mrSges(6,1) - t382 * mrSges(6,3);
t428 = -t381 * mrSges(6,1) - t382 * mrSges(6,2) - t361;
t436 = -mrSges(6,3) - mrSges(7,2);
t326 = m(6) * t337 - t380 * mrSges(6,2) + t436 * t353 - t390 * t369 + t428 * t381 + t421;
t368 = -t390 * mrSges(6,2) - t381 * mrSges(6,3);
t328 = m(6) * t336 + t380 * mrSges(6,1) + t436 * t354 + t390 * t368 + t428 * t382 + t417;
t322 = t402 * t326 + t437 * t328;
t387 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t426;
t388 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t425;
t439 = m(4) * t366 - m(5) * t363 + t385 * mrSges(5,1) - t386 * mrSges(5,2) - t322 + (-t387 * t403 - t388 * t405) * qJD(1);
t355 = t403 * g(3) + t405 * t364;
t383 = (t403 * mrSges(5,1) + t405 * mrSges(5,2)) * qJD(1);
t341 = -qJDD(4) * pkin(4) - t407 * pkin(8) + t384 * t425 - t355;
t335 = -0.2e1 * qJD(6) * t382 + (t381 * t390 - t354) * qJ(6) + (t382 * t390 + t353) * pkin(5) + t341;
t329 = m(7) * t335 + t353 * mrSges(7,1) - t354 * mrSges(7,3) - t382 * t370 + t381 * t371;
t409 = -m(6) * t341 - t353 * mrSges(6,1) - t354 * mrSges(6,2) - t381 * t368 - t382 * t369 - t329;
t418 = t437 * t326 - t402 * t328;
t432 = t403 * (m(5) * t356 - qJDD(4) * mrSges(5,2) + t385 * mrSges(5,3) - qJD(4) * t388 - t383 * t426 + t418) + t405 * (m(5) * t355 + qJDD(4) * mrSges(5,1) - t386 * mrSges(5,3) + qJD(4) * t387 - t383 * t425 + t409);
t431 = t433 * t381 + t443 * t382 + t442 * t390;
t413 = m(4) * t367 + qJDD(1) * mrSges(4,2) - t408 * mrSges(4,3) + t432;
t377 = Ifges(5,5) * qJD(4) + (t405 * Ifges(5,1) - t403 * Ifges(5,4)) * qJD(1);
t376 = Ifges(5,6) * qJD(4) + (t405 * Ifges(5,4) - t403 * Ifges(5,2)) * qJD(1);
t372 = t408 * pkin(1) - t438;
t320 = mrSges(6,2) * t341 + mrSges(7,2) * t334 - mrSges(6,3) * t336 - mrSges(7,3) * t335 - qJ(6) * t329 - t435 * t353 + t447 * t354 - t443 * t380 + t431 * t381 + t430 * t390;
t319 = -mrSges(6,1) * t341 - mrSges(7,1) * t335 + mrSges(7,2) * t333 + mrSges(6,3) * t337 - pkin(5) * t329 - t446 * t353 + t435 * t354 + t433 * t380 + t431 * t382 - t440 * t390;
t318 = m(3) * t373 + (-mrSges(4,2) - mrSges(3,3)) * t408 + (mrSges(3,2) - mrSges(4,3)) * qJDD(1) + t439;
t1 = [qJ(2) * (-m(3) * t372 + t408 * mrSges(3,2) + t413) - pkin(1) * t318 + mrSges(2,1) * t427 - mrSges(2,2) * t416 + mrSges(3,2) * t373 - mrSges(3,3) * t372 + t405 * (mrSges(5,2) * t363 - mrSges(5,3) * t355 + Ifges(5,1) * t386 + Ifges(5,4) * t385 + Ifges(5,5) * qJDD(4) - pkin(8) * t322 - qJD(4) * t376 - t402 * t319 + t437 * t320) - t403 * (-mrSges(5,1) * t363 + mrSges(5,3) * t356 + Ifges(5,4) * t386 + Ifges(5,2) * t385 + Ifges(5,6) * qJDD(4) - pkin(4) * t322 + qJD(4) * t377 - t445) - pkin(7) * t432 + mrSges(4,2) * t367 - mrSges(4,3) * t366 + (qJ(2) * mrSges(3,3) + Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1) + (t408 * mrSges(4,2) + qJDD(1) * mrSges(4,3) - t439) * qJ(3); t318; t413; Ifges(5,5) * t386 + Ifges(5,6) * t385 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t355 - mrSges(5,2) * t356 + t402 * t320 + t437 * t319 + pkin(4) * t409 + pkin(8) * t418 + (t405 * t376 + t403 * t377) * qJD(1); t445; t331;];
tauJ  = t1;
