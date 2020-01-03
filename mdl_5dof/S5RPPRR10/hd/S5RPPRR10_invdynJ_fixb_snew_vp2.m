% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPPRR10
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPRR10_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR10_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR10_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR10_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR10_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR10_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR10_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:04:09
% EndTime: 2019-12-31 18:04:10
% DurationCPUTime: 1.20s
% Computational Cost: add. (5571->203), mult. (13438->260), div. (0->0), fcn. (8986->8), ass. (0->90)
t419 = qJD(1) ^ 2;
t415 = sin(qJ(1));
t418 = cos(qJ(1));
t429 = -t418 * g(1) - t415 * g(2);
t450 = -t419 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t429;
t449 = mrSges(4,2) + mrSges(3,3);
t412 = cos(pkin(8));
t407 = t412 ^ 2;
t411 = sin(pkin(8));
t433 = (-t411 ^ 2 - t407) * t419;
t442 = t415 * g(1) - t418 * g(2);
t385 = -qJDD(1) * pkin(1) - t419 * qJ(2) + qJDD(2) - t442;
t436 = qJDD(1) * t412;
t440 = qJD(1) * t411;
t446 = t411 * qJ(3);
t373 = -pkin(2) * t436 - 0.2e1 * qJD(3) * t440 - qJDD(1) * t446 + t385;
t376 = -t412 * g(3) - t450 * t411;
t448 = pkin(3) * t419;
t447 = Ifges(3,4) - Ifges(4,5);
t390 = (-t412 * pkin(2) - t446) * qJD(1);
t361 = t390 * t440 + qJDD(3) - t376;
t357 = (-pkin(6) * qJDD(1) - t412 * t448) * t411 + t361;
t377 = -t411 * g(3) + t450 * t412;
t439 = qJD(1) * t412;
t362 = t390 * t439 + t377;
t359 = -pkin(6) * t436 - t407 * t448 + t362;
t414 = sin(qJ(4));
t417 = cos(qJ(4));
t340 = t417 * t357 - t414 * t359;
t426 = t411 * t417 - t412 * t414;
t425 = -t411 * t414 - t412 * t417;
t387 = t425 * qJD(1);
t438 = t387 * qJD(4);
t375 = t426 * qJDD(1) + t438;
t388 = t426 * qJD(1);
t336 = (-t375 + t438) * pkin(7) + (t387 * t388 + qJDD(4)) * pkin(4) + t340;
t341 = t414 * t357 + t417 * t359;
t374 = -t388 * qJD(4) + t425 * qJDD(1);
t380 = qJD(4) * pkin(4) - t388 * pkin(7);
t386 = t387 ^ 2;
t337 = -t386 * pkin(4) + t374 * pkin(7) - qJD(4) * t380 + t341;
t413 = sin(qJ(5));
t416 = cos(qJ(5));
t334 = t416 * t336 - t413 * t337;
t368 = t416 * t387 - t413 * t388;
t348 = t368 * qJD(5) + t413 * t374 + t416 * t375;
t369 = t413 * t387 + t416 * t388;
t354 = -t368 * mrSges(6,1) + t369 * mrSges(6,2);
t408 = qJD(4) + qJD(5);
t363 = -t408 * mrSges(6,2) + t368 * mrSges(6,3);
t405 = qJDD(4) + qJDD(5);
t331 = m(6) * t334 + t405 * mrSges(6,1) - t348 * mrSges(6,3) - t369 * t354 + t408 * t363;
t335 = t413 * t336 + t416 * t337;
t347 = -t369 * qJD(5) + t416 * t374 - t413 * t375;
t364 = t408 * mrSges(6,1) - t369 * mrSges(6,3);
t332 = m(6) * t335 - t405 * mrSges(6,2) + t347 * mrSges(6,3) + t368 * t354 - t408 * t364;
t323 = t416 * t331 + t413 * t332;
t371 = -t387 * mrSges(5,1) + t388 * mrSges(5,2);
t378 = -qJD(4) * mrSges(5,2) + t387 * mrSges(5,3);
t321 = m(5) * t340 + qJDD(4) * mrSges(5,1) - t375 * mrSges(5,3) + qJD(4) * t378 - t388 * t371 + t323;
t379 = qJD(4) * mrSges(5,1) - t388 * mrSges(5,3);
t431 = -t413 * t331 + t416 * t332;
t322 = m(5) * t341 - qJDD(4) * mrSges(5,2) + t374 * mrSges(5,3) - qJD(4) * t379 + t387 * t371 + t431;
t445 = t417 * t321 + t414 * t322;
t443 = ((Ifges(3,6) - Ifges(4,6)) * t412 + (Ifges(4,4) + Ifges(3,5)) * t411) * qJD(1);
t432 = -t414 * t321 + t417 * t322;
t430 = m(4) * t361 + t445;
t428 = -t412 * mrSges(4,1) - t411 * mrSges(4,3);
t360 = pkin(3) * t436 + pkin(6) * t433 - t373;
t339 = -t374 * pkin(4) - t386 * pkin(7) + t388 * t380 + t360;
t427 = m(6) * t339 - t347 * mrSges(6,1) + t348 * mrSges(6,2) - t368 * t363 + t369 * t364;
t391 = t428 * qJD(1);
t424 = qJDD(1) * mrSges(4,2) + qJD(1) * t391;
t350 = Ifges(6,4) * t369 + Ifges(6,2) * t368 + Ifges(6,6) * t408;
t351 = Ifges(6,1) * t369 + Ifges(6,4) * t368 + Ifges(6,5) * t408;
t422 = mrSges(6,1) * t334 - mrSges(6,2) * t335 + Ifges(6,5) * t348 + Ifges(6,6) * t347 + Ifges(6,3) * t405 + t369 * t350 - t368 * t351;
t421 = -m(5) * t360 + t374 * mrSges(5,1) - t375 * mrSges(5,2) + t387 * t378 - t388 * t379 - t427;
t420 = m(4) * t373 + t421;
t392 = (-t412 * mrSges(3,1) + t411 * mrSges(3,2)) * qJD(1);
t367 = Ifges(5,1) * t388 + Ifges(5,4) * t387 + Ifges(5,5) * qJD(4);
t366 = Ifges(5,4) * t388 + Ifges(5,2) * t387 + Ifges(5,6) * qJD(4);
t365 = Ifges(5,5) * t388 + Ifges(5,6) * t387 + Ifges(5,3) * qJD(4);
t349 = Ifges(6,5) * t369 + Ifges(6,6) * t368 + Ifges(6,3) * t408;
t327 = mrSges(4,2) * t433 + t428 * qJDD(1) + t420;
t326 = ((-mrSges(3,1) - mrSges(4,1)) * t412 + (mrSges(3,2) - mrSges(4,3)) * t411) * qJDD(1) + t420 + m(3) * t385 + t449 * t433;
t325 = mrSges(6,2) * t339 - mrSges(6,3) * t334 + Ifges(6,1) * t348 + Ifges(6,4) * t347 + Ifges(6,5) * t405 + t368 * t349 - t408 * t350;
t324 = -mrSges(6,1) * t339 + mrSges(6,3) * t335 + Ifges(6,4) * t348 + Ifges(6,2) * t347 + Ifges(6,6) * t405 - t369 * t349 + t408 * t351;
t317 = mrSges(5,2) * t360 - mrSges(5,3) * t340 + Ifges(5,1) * t375 + Ifges(5,4) * t374 + Ifges(5,5) * qJDD(4) - pkin(7) * t323 - qJD(4) * t366 - t413 * t324 + t416 * t325 + t387 * t365;
t316 = -mrSges(5,1) * t360 + mrSges(5,3) * t341 + Ifges(5,4) * t375 + Ifges(5,2) * t374 + Ifges(5,6) * qJDD(4) - pkin(4) * t427 + pkin(7) * t431 + qJD(4) * t367 + t416 * t324 + t413 * t325 - t388 * t365;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t442 - mrSges(2,2) * t429 + t411 * (mrSges(3,2) * t385 - mrSges(3,3) * t376 + mrSges(4,2) * t361 - mrSges(4,3) * t373 + t417 * t317 - t414 * t316 - pkin(6) * t445 - qJ(3) * t327 + t443 * t439 + (t412 * t447 + (Ifges(3,1) + Ifges(4,1)) * t411) * qJDD(1)) + t412 * (-mrSges(3,1) * t385 + mrSges(3,3) * t377 - mrSges(4,1) * t373 + mrSges(4,2) * t362 - t414 * t317 - t417 * t316 - pkin(3) * t421 - pkin(6) * t432 - pkin(2) * t327 - t443 * t440 + ((Ifges(3,2) + Ifges(4,3)) * t412 + t411 * t447) * qJDD(1)) - pkin(1) * t326 + qJ(2) * (t412 * (m(3) * t377 + m(4) * t362 + (t449 * qJDD(1) + (t391 + t392) * qJD(1)) * t412 + t432) + (-m(3) * t376 + (qJDD(1) * mrSges(3,3) + qJD(1) * t392 + t424) * t411 + t430) * t411); t326; t424 * t411 + t430; mrSges(5,1) * t340 - mrSges(5,2) * t341 + Ifges(5,5) * t375 + Ifges(5,6) * t374 + Ifges(5,3) * qJDD(4) + pkin(4) * t323 + t388 * t366 - t387 * t367 + t422; t422;];
tauJ = t1;
