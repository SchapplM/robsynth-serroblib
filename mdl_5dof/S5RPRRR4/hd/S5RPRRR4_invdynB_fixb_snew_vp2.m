% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRRR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2020-01-03 11:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRRR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR4_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR4_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR4_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:52:13
% EndTime: 2020-01-03 11:52:16
% DurationCPUTime: 2.71s
% Computational Cost: add. (36781->185), mult. (47778->233), div. (0->0), fcn. (25560->10), ass. (0->82)
t366 = qJD(1) + qJD(3);
t359 = qJD(4) + t366;
t370 = sin(qJ(5));
t391 = t359 * t370;
t374 = cos(qJ(5));
t390 = t359 * t374;
t373 = sin(qJ(1));
t377 = cos(qJ(1));
t355 = -t373 * g(2) + t377 * g(3);
t378 = qJD(1) ^ 2;
t356 = -t377 * g(2) - t373 * g(3);
t353 = qJDD(1) * pkin(1) + t356;
t354 = -t378 * pkin(1) + t355;
t368 = sin(pkin(9));
t369 = cos(pkin(9));
t338 = t369 * t353 - t368 * t354;
t336 = qJDD(1) * pkin(2) + t338;
t339 = t368 * t353 + t369 * t354;
t337 = -t378 * pkin(2) + t339;
t372 = sin(qJ(3));
t376 = cos(qJ(3));
t331 = t376 * t336 - t372 * t337;
t365 = qJDD(1) + qJDD(3);
t329 = t365 * pkin(3) + t331;
t332 = t372 * t336 + t376 * t337;
t364 = t366 ^ 2;
t330 = -t364 * pkin(3) + t332;
t371 = sin(qJ(4));
t375 = cos(qJ(4));
t326 = t371 * t329 + t375 * t330;
t357 = t359 ^ 2;
t358 = qJDD(4) + t365;
t324 = -t357 * pkin(4) + t358 * pkin(8) + t326;
t367 = -g(1) + qJDD(2);
t321 = -t370 * t324 + t374 * t367;
t345 = (-mrSges(6,1) * t374 + mrSges(6,2) * t370) * t359;
t388 = qJD(5) * t359;
t346 = t370 * t358 + t374 * t388;
t352 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t390;
t319 = m(6) * t321 + qJDD(5) * mrSges(6,1) - t346 * mrSges(6,3) + qJD(5) * t352 - t345 * t391;
t322 = t374 * t324 + t370 * t367;
t347 = t374 * t358 - t370 * t388;
t351 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t391;
t320 = m(6) * t322 - qJDD(5) * mrSges(6,2) + t347 * mrSges(6,3) - qJD(5) * t351 + t345 * t390;
t382 = -t370 * t319 + t374 * t320;
t310 = m(5) * t326 - t357 * mrSges(5,1) - t358 * mrSges(5,2) + t382;
t325 = t375 * t329 - t371 * t330;
t323 = -t358 * pkin(4) - t357 * pkin(8) - t325;
t379 = -m(6) * t323 + t347 * mrSges(6,1) - t346 * mrSges(6,2) - t351 * t391 + t352 * t390;
t315 = m(5) * t325 + t358 * mrSges(5,1) - t357 * mrSges(5,2) + t379;
t307 = t371 * t310 + t375 * t315;
t304 = m(4) * t331 + t365 * mrSges(4,1) - t364 * mrSges(4,2) + t307;
t383 = t375 * t310 - t371 * t315;
t305 = m(4) * t332 - t364 * mrSges(4,1) - t365 * mrSges(4,2) + t383;
t299 = t376 * t304 + t372 * t305;
t297 = m(3) * t338 + qJDD(1) * mrSges(3,1) - t378 * mrSges(3,2) + t299;
t384 = -t372 * t304 + t376 * t305;
t298 = m(3) * t339 - t378 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t384;
t385 = -t368 * t297 + t369 * t298;
t289 = m(2) * t355 - t378 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t385;
t291 = t369 * t297 + t368 * t298;
t290 = m(2) * t356 + qJDD(1) * mrSges(2,1) - t378 * mrSges(2,2) + t291;
t389 = t373 * t289 + t377 * t290;
t311 = t374 * t319 + t370 * t320;
t387 = m(5) * t367 + t311;
t386 = -t377 * t289 + t373 * t290;
t381 = m(4) * t367 + t387;
t380 = m(3) * t367 + t381;
t342 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t370 + Ifges(6,4) * t374) * t359;
t341 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t370 + Ifges(6,2) * t374) * t359;
t340 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t370 + Ifges(6,6) * t374) * t359;
t313 = mrSges(6,2) * t323 - mrSges(6,3) * t321 + Ifges(6,1) * t346 + Ifges(6,4) * t347 + Ifges(6,5) * qJDD(5) - qJD(5) * t341 + t340 * t390;
t312 = -mrSges(6,1) * t323 + mrSges(6,3) * t322 + Ifges(6,4) * t346 + Ifges(6,2) * t347 + Ifges(6,6) * qJDD(5) + qJD(5) * t342 - t340 * t391;
t306 = -mrSges(5,1) * t367 - mrSges(6,1) * t321 + mrSges(6,2) * t322 + mrSges(5,3) * t326 + t357 * Ifges(5,5) - Ifges(6,5) * t346 + Ifges(5,6) * t358 - Ifges(6,6) * t347 - Ifges(6,3) * qJDD(5) - pkin(4) * t311 + (-t341 * t370 + t342 * t374) * t359;
t300 = mrSges(5,2) * t367 - mrSges(5,3) * t325 + Ifges(5,5) * t358 - t357 * Ifges(5,6) - pkin(8) * t311 - t370 * t312 + t374 * t313;
t293 = mrSges(4,2) * t367 - mrSges(4,3) * t331 + Ifges(4,5) * t365 - t364 * Ifges(4,6) - pkin(7) * t307 + t375 * t300 - t371 * t306;
t292 = -mrSges(4,1) * t367 + mrSges(4,3) * t332 + t364 * Ifges(4,5) + Ifges(4,6) * t365 - pkin(3) * t387 + pkin(7) * t383 + t371 * t300 + t375 * t306;
t285 = mrSges(3,2) * t367 - mrSges(3,3) * t338 + Ifges(3,5) * qJDD(1) - t378 * Ifges(3,6) - pkin(6) * t299 - t372 * t292 + t376 * t293;
t284 = -mrSges(3,1) * t367 + mrSges(3,3) * t339 + t378 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t381 + pkin(6) * t384 + t376 * t292 + t372 * t293;
t283 = -mrSges(2,2) * g(1) - mrSges(2,3) * t356 + Ifges(2,5) * qJDD(1) - t378 * Ifges(2,6) - qJ(2) * t291 - t368 * t284 + t369 * t285;
t282 = mrSges(2,1) * g(1) + mrSges(2,3) * t355 + t378 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t380 + qJ(2) * t385 + t369 * t284 + t368 * t285;
t1 = [(-m(1) - m(2)) * g(1) + t380; -m(1) * g(2) + t389; -m(1) * g(3) + t386; pkin(1) * t291 + pkin(2) * t299 + mrSges(3,1) * t338 - mrSges(3,2) * t339 + pkin(3) * t307 + mrSges(4,1) * t331 - mrSges(4,2) * t332 + t374 * t312 + pkin(4) * t379 + pkin(8) * t382 + mrSges(5,1) * t325 - mrSges(5,2) * t326 + t370 * t313 + mrSges(2,1) * t356 - mrSges(2,2) * t355 + Ifges(4,3) * t365 + Ifges(5,3) * t358 - mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1); mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t386 + t377 * t282 + t373 * t283; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t389 + t373 * t282 - t377 * t283;];
tauB = t1;
