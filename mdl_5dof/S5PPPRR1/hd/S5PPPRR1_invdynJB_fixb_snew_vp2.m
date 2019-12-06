% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PPPRR1
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
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
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PPPRR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPPRR1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPPRR1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPPRR1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:58:01
% EndTime: 2019-12-05 14:58:02
% DurationCPUTime: 1.18s
% Computational Cost: add. (12902->151), mult. (19761->200), div. (0->0), fcn. (14386->10), ass. (0->74)
t368 = sin(pkin(7));
t371 = cos(pkin(7));
t361 = -t371 * g(1) - t368 * g(2);
t365 = -g(3) + qJDD(1);
t367 = sin(pkin(8));
t370 = cos(pkin(8));
t350 = t370 * t361 + t367 * t365;
t360 = t368 * g(1) - t371 * g(2);
t359 = qJDD(2) - t360;
t366 = sin(pkin(9));
t369 = cos(pkin(9));
t346 = -t366 * t350 + t369 * t359;
t347 = t369 * t350 + t366 * t359;
t373 = sin(qJ(4));
t375 = cos(qJ(4));
t343 = t373 * t346 + t375 * t347;
t376 = qJD(4) ^ 2;
t341 = -t376 * pkin(4) + qJDD(4) * pkin(6) + t343;
t349 = -t367 * t361 + t370 * t365;
t348 = qJDD(3) - t349;
t372 = sin(qJ(5));
t374 = cos(qJ(5));
t338 = -t372 * t341 + t374 * t348;
t356 = (-mrSges(6,1) * t374 + mrSges(6,2) * t372) * qJD(4);
t385 = qJD(4) * qJD(5);
t357 = t372 * qJDD(4) + t374 * t385;
t386 = qJD(4) * t374;
t363 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t386;
t387 = qJD(4) * t372;
t335 = m(6) * t338 + qJDD(5) * mrSges(6,1) - t357 * mrSges(6,3) + qJD(5) * t363 - t356 * t387;
t339 = t374 * t341 + t372 * t348;
t358 = t374 * qJDD(4) - t372 * t385;
t362 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t387;
t336 = m(6) * t339 - qJDD(5) * mrSges(6,2) + t358 * mrSges(6,3) - qJD(5) * t362 + t356 * t386;
t326 = t374 * t335 + t372 * t336;
t325 = (m(4) + m(5)) * t348 + t326;
t352 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t372 + Ifges(6,2) * t374) * qJD(4);
t353 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t372 + Ifges(6,4) * t374) * qJD(4);
t390 = mrSges(6,1) * t338 - mrSges(6,2) * t339 + Ifges(6,5) * t357 + Ifges(6,6) * t358 + Ifges(6,3) * qJDD(5) + (t352 * t372 - t353 * t374) * qJD(4);
t327 = -t372 * t335 + t374 * t336;
t322 = m(5) * t343 - t376 * mrSges(5,1) - qJDD(4) * mrSges(5,2) + t327;
t342 = t375 * t346 - t373 * t347;
t340 = -qJDD(4) * pkin(4) - t376 * pkin(6) - t342;
t337 = -m(6) * t340 + t358 * mrSges(6,1) - t357 * mrSges(6,2) - t362 * t387 + t363 * t386;
t331 = m(5) * t342 + qJDD(4) * mrSges(5,1) - t376 * mrSges(5,2) + t337;
t319 = t373 * t322 + t375 * t331;
t317 = m(4) * t346 + t319;
t380 = t375 * t322 - t373 * t331;
t318 = m(4) * t347 + t380;
t381 = -t366 * t317 + t369 * t318;
t311 = m(3) * t350 + t381;
t324 = m(3) * t349 - t325;
t382 = t370 * t311 - t367 * t324;
t305 = m(2) * t361 + t382;
t313 = t369 * t317 + t366 * t318;
t312 = m(3) * t359 + t313;
t310 = m(2) * t360 - t312;
t388 = t368 * t305 + t371 * t310;
t306 = t367 * t311 + t370 * t324;
t384 = m(2) * t365 + t306;
t383 = t371 * t305 - t368 * t310;
t351 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t372 + Ifges(6,6) * t374) * qJD(4);
t328 = -mrSges(6,1) * t340 + mrSges(6,3) * t339 + Ifges(6,4) * t357 + Ifges(6,2) * t358 + Ifges(6,6) * qJDD(5) + qJD(5) * t353 - t351 * t387;
t329 = mrSges(6,2) * t340 - mrSges(6,3) * t338 + Ifges(6,1) * t357 + Ifges(6,4) * t358 + Ifges(6,5) * qJDD(5) - qJD(5) * t352 + t351 * t386;
t377 = mrSges(5,1) * t342 - mrSges(5,2) * t343 + Ifges(5,3) * qJDD(4) + pkin(4) * t337 + pkin(6) * t327 + t374 * t328 + t372 * t329;
t315 = -mrSges(5,1) * t348 + mrSges(5,3) * t343 + t376 * Ifges(5,5) + Ifges(5,6) * qJDD(4) - pkin(4) * t326 - t390;
t314 = mrSges(5,2) * t348 - mrSges(5,3) * t342 + Ifges(5,5) * qJDD(4) - t376 * Ifges(5,6) - pkin(6) * t326 - t372 * t328 + t374 * t329;
t302 = mrSges(4,2) * t348 - mrSges(4,3) * t346 - pkin(5) * t319 + t375 * t314 - t373 * t315;
t301 = -mrSges(4,1) * t348 + mrSges(4,3) * t347 + t373 * t314 + t375 * t315 - pkin(3) * (m(5) * t348 + t326) + pkin(5) * t380;
t300 = -mrSges(3,1) * t359 - mrSges(4,1) * t346 + mrSges(4,2) * t347 + mrSges(3,3) * t350 - pkin(2) * t313 - pkin(3) * t319 - t377;
t299 = mrSges(3,2) * t359 - mrSges(3,3) * t349 - qJ(3) * t313 - t366 * t301 + t369 * t302;
t298 = -mrSges(2,1) * t365 - mrSges(3,1) * t349 + mrSges(3,2) * t350 + mrSges(2,3) * t361 - pkin(1) * t306 + pkin(2) * t325 - qJ(3) * t381 - t369 * t301 - t366 * t302;
t297 = mrSges(2,2) * t365 - mrSges(2,3) * t360 - qJ(2) * t306 + t370 * t299 - t367 * t300;
t1 = [-m(1) * g(1) + t383; -m(1) * g(2) + t388; -m(1) * g(3) + t384; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t388 + t371 * t297 - t368 * t298; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t383 + t368 * t297 + t371 * t298; -mrSges(1,1) * g(2) + mrSges(2,1) * t360 + mrSges(1,2) * g(1) - mrSges(2,2) * t361 - pkin(1) * t312 + qJ(2) * t382 + t367 * t299 + t370 * t300; t384; t312; t325; t377; t390;];
tauJB = t1;
