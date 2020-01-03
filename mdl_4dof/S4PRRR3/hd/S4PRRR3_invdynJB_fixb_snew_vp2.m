% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4PRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+4)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4PRRR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR3_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR3_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR3_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR3_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR3_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:37
% EndTime: 2019-12-31 16:31:38
% DurationCPUTime: 0.71s
% Computational Cost: add. (8032->142), mult. (11131->187), div. (0->0), fcn. (6720->8), ass. (0->68)
t362 = sin(pkin(7));
t363 = cos(pkin(7));
t351 = t362 * g(1) - t363 * g(2);
t352 = -t363 * g(1) - t362 * g(2);
t366 = sin(qJ(2));
t369 = cos(qJ(2));
t336 = t369 * t351 - t366 * t352;
t333 = qJDD(2) * pkin(2) + t336;
t337 = t366 * t351 + t369 * t352;
t370 = qJD(2) ^ 2;
t334 = -t370 * pkin(2) + t337;
t365 = sin(qJ(3));
t368 = cos(qJ(3));
t330 = t365 * t333 + t368 * t334;
t359 = qJD(2) + qJD(3);
t357 = t359 ^ 2;
t358 = qJDD(2) + qJDD(3);
t327 = -t357 * pkin(3) + t358 * pkin(6) + t330;
t361 = -g(3) + qJDD(1);
t364 = sin(qJ(4));
t367 = cos(qJ(4));
t324 = -t364 * t327 + t367 * t361;
t325 = t367 * t327 + t364 * t361;
t339 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t364 + Ifges(5,2) * t367) * t359;
t340 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t364 + Ifges(5,4) * t367) * t359;
t383 = qJD(4) * t359;
t344 = t364 * t358 + t367 * t383;
t345 = t367 * t358 - t364 * t383;
t387 = mrSges(5,1) * t324 - mrSges(5,2) * t325 + Ifges(5,5) * t344 + Ifges(5,6) * t345 + Ifges(5,3) * qJDD(4) + (t339 * t364 - t340 * t367) * t359;
t386 = t359 * t364;
t385 = t359 * t367;
t343 = (-mrSges(5,1) * t367 + mrSges(5,2) * t364) * t359;
t350 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t385;
t322 = m(5) * t324 + qJDD(4) * mrSges(5,1) - t344 * mrSges(5,3) + qJD(4) * t350 - t343 * t386;
t349 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t386;
t323 = m(5) * t325 - qJDD(4) * mrSges(5,2) + t345 * mrSges(5,3) - qJD(4) * t349 + t343 * t385;
t378 = -t364 * t322 + t367 * t323;
t309 = m(4) * t330 - t357 * mrSges(4,1) - t358 * mrSges(4,2) + t378;
t329 = t368 * t333 - t365 * t334;
t326 = -t358 * pkin(3) - t357 * pkin(6) - t329;
t373 = -m(5) * t326 + t345 * mrSges(5,1) - t344 * mrSges(5,2) - t349 * t386 + t350 * t385;
t317 = m(4) * t329 + t358 * mrSges(4,1) - t357 * mrSges(4,2) + t373;
t306 = t365 * t309 + t368 * t317;
t303 = m(3) * t336 + qJDD(2) * mrSges(3,1) - t370 * mrSges(3,2) + t306;
t379 = t368 * t309 - t365 * t317;
t304 = m(3) * t337 - t370 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t379;
t297 = t369 * t303 + t366 * t304;
t295 = m(2) * t351 + t297;
t380 = -t366 * t303 + t369 * t304;
t296 = m(2) * t352 + t380;
t384 = t363 * t295 + t362 * t296;
t311 = t367 * t322 + t364 * t323;
t382 = m(4) * t361 + t311;
t381 = -t362 * t295 + t363 * t296;
t377 = m(3) * t361 + t382;
t376 = m(2) * t361 + t377;
t338 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t364 + Ifges(5,6) * t367) * t359;
t314 = -mrSges(5,1) * t326 + mrSges(5,3) * t325 + Ifges(5,4) * t344 + Ifges(5,2) * t345 + Ifges(5,6) * qJDD(4) + qJD(4) * t340 - t338 * t386;
t315 = mrSges(5,2) * t326 - mrSges(5,3) * t324 + Ifges(5,1) * t344 + Ifges(5,4) * t345 + Ifges(5,5) * qJDD(4) - qJD(4) * t339 + t338 * t385;
t374 = mrSges(4,1) * t329 - mrSges(4,2) * t330 + Ifges(4,3) * t358 + pkin(3) * t373 + pkin(6) * t378 + t367 * t314 + t364 * t315;
t371 = mrSges(3,1) * t336 - mrSges(3,2) * t337 + Ifges(3,3) * qJDD(2) + pkin(2) * t306 + t374;
t299 = -mrSges(4,1) * t361 + mrSges(4,3) * t330 + t357 * Ifges(4,5) + Ifges(4,6) * t358 - pkin(3) * t311 - t387;
t298 = mrSges(4,2) * t361 - mrSges(4,3) * t329 + Ifges(4,5) * t358 - t357 * Ifges(4,6) - pkin(6) * t311 - t364 * t314 + t367 * t315;
t291 = mrSges(3,2) * t361 - mrSges(3,3) * t336 + Ifges(3,5) * qJDD(2) - t370 * Ifges(3,6) - pkin(5) * t306 + t368 * t298 - t365 * t299;
t290 = -mrSges(3,1) * t361 + mrSges(3,3) * t337 + t370 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t382 + pkin(5) * t379 + t365 * t298 + t368 * t299;
t289 = mrSges(2,2) * t361 - mrSges(2,3) * t351 - pkin(4) * t297 - t366 * t290 + t369 * t291;
t288 = -mrSges(2,1) * t361 + mrSges(2,3) * t352 - pkin(1) * t377 + pkin(4) * t380 + t369 * t290 + t366 * t291;
t1 = [-m(1) * g(1) + t381; -m(1) * g(2) + t384; -m(1) * g(3) + t376; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t384 - t362 * t288 + t363 * t289; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t381 + t363 * t288 + t362 * t289; -mrSges(1,1) * g(2) + mrSges(2,1) * t351 + mrSges(1,2) * g(1) - mrSges(2,2) * t352 + pkin(1) * t297 + t371; t376; t371; t374; t387;];
tauJB = t1;
