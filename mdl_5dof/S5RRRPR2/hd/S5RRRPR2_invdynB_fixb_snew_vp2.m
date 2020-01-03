% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRPR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:07:31
% EndTime: 2020-01-03 12:07:33
% DurationCPUTime: 2.17s
% Computational Cost: add. (41161->187), mult. (47778->233), div. (0->0), fcn. (25560->10), ass. (0->81)
t386 = -m(3) - m(4);
t362 = qJD(1) + qJD(2);
t357 = qJD(3) + t362;
t366 = sin(qJ(5));
t385 = t357 * t366;
t370 = cos(qJ(5));
t384 = t357 * t370;
t369 = sin(qJ(1));
t373 = cos(qJ(1));
t353 = -t369 * g(2) + t373 * g(3);
t374 = qJD(1) ^ 2;
t354 = -t373 * g(2) - t369 * g(3);
t351 = qJDD(1) * pkin(1) + t354;
t352 = -t374 * pkin(1) + t353;
t368 = sin(qJ(2));
t372 = cos(qJ(2));
t336 = t372 * t351 - t368 * t352;
t361 = qJDD(1) + qJDD(2);
t334 = t361 * pkin(2) + t336;
t337 = t368 * t351 + t372 * t352;
t360 = t362 ^ 2;
t335 = -t360 * pkin(2) + t337;
t367 = sin(qJ(3));
t371 = cos(qJ(3));
t329 = t371 * t334 - t367 * t335;
t356 = qJDD(3) + t361;
t327 = t356 * pkin(3) + t329;
t330 = t367 * t334 + t371 * t335;
t355 = t357 ^ 2;
t328 = -t355 * pkin(3) + t330;
t364 = sin(pkin(9));
t365 = cos(pkin(9));
t324 = t364 * t327 + t365 * t328;
t322 = -t355 * pkin(4) + t356 * pkin(8) + t324;
t363 = -g(1) + qJDD(4);
t319 = -t366 * t322 + t370 * t363;
t343 = (-mrSges(6,1) * t370 + mrSges(6,2) * t366) * t357;
t382 = qJD(5) * t357;
t344 = t366 * t356 + t370 * t382;
t350 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t384;
t317 = m(6) * t319 + qJDD(5) * mrSges(6,1) - t344 * mrSges(6,3) + qJD(5) * t350 - t343 * t385;
t320 = t370 * t322 + t366 * t363;
t345 = t370 * t356 - t366 * t382;
t349 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t385;
t318 = m(6) * t320 - qJDD(5) * mrSges(6,2) + t345 * mrSges(6,3) - qJD(5) * t349 + t343 * t384;
t376 = -t366 * t317 + t370 * t318;
t308 = m(5) * t324 - t355 * mrSges(5,1) - t356 * mrSges(5,2) + t376;
t323 = t365 * t327 - t364 * t328;
t321 = -t356 * pkin(4) - t355 * pkin(8) - t323;
t375 = -m(6) * t321 + t345 * mrSges(6,1) - t344 * mrSges(6,2) - t349 * t385 + t350 * t384;
t313 = m(5) * t323 + t356 * mrSges(5,1) - t355 * mrSges(5,2) + t375;
t305 = t364 * t308 + t365 * t313;
t302 = m(4) * t329 + t356 * mrSges(4,1) - t355 * mrSges(4,2) + t305;
t377 = t365 * t308 - t364 * t313;
t303 = m(4) * t330 - t355 * mrSges(4,1) - t356 * mrSges(4,2) + t377;
t297 = t371 * t302 + t367 * t303;
t295 = m(3) * t336 + t361 * mrSges(3,1) - t360 * mrSges(3,2) + t297;
t378 = -t367 * t302 + t371 * t303;
t296 = m(3) * t337 - t360 * mrSges(3,1) - t361 * mrSges(3,2) + t378;
t379 = -t368 * t295 + t372 * t296;
t287 = m(2) * t353 - t374 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t379;
t289 = t372 * t295 + t368 * t296;
t288 = m(2) * t354 + qJDD(1) * mrSges(2,1) - t374 * mrSges(2,2) + t289;
t383 = t369 * t287 + t373 * t288;
t309 = t370 * t317 + t366 * t318;
t381 = m(5) * t363 + t309;
t380 = -t373 * t287 + t369 * t288;
t340 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t366 + Ifges(6,4) * t370) * t357;
t339 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t366 + Ifges(6,2) * t370) * t357;
t338 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t366 + Ifges(6,6) * t370) * t357;
t311 = mrSges(6,2) * t321 - mrSges(6,3) * t319 + Ifges(6,1) * t344 + Ifges(6,4) * t345 + Ifges(6,5) * qJDD(5) - qJD(5) * t339 + t338 * t384;
t310 = -mrSges(6,1) * t321 + mrSges(6,3) * t320 + Ifges(6,4) * t344 + Ifges(6,2) * t345 + Ifges(6,6) * qJDD(5) + qJD(5) * t340 - t338 * t385;
t304 = -mrSges(5,1) * t363 - mrSges(6,1) * t319 + mrSges(6,2) * t320 + mrSges(5,3) * t324 + t355 * Ifges(5,5) - Ifges(6,5) * t344 + Ifges(5,6) * t356 - Ifges(6,6) * t345 - Ifges(6,3) * qJDD(5) - pkin(4) * t309 + (-t339 * t366 + t340 * t370) * t357;
t298 = mrSges(5,2) * t363 - mrSges(5,3) * t323 + Ifges(5,5) * t356 - t355 * Ifges(5,6) - pkin(8) * t309 - t366 * t310 + t370 * t311;
t291 = -mrSges(4,2) * g(1) - mrSges(4,3) * t329 + Ifges(4,5) * t356 - t355 * Ifges(4,6) - qJ(4) * t305 + t365 * t298 - t364 * t304;
t290 = mrSges(4,1) * g(1) + mrSges(4,3) * t330 + t355 * Ifges(4,5) + Ifges(4,6) * t356 - pkin(3) * t381 + qJ(4) * t377 + t364 * t298 + t365 * t304;
t283 = -mrSges(3,2) * g(1) - mrSges(3,3) * t336 + Ifges(3,5) * t361 - t360 * Ifges(3,6) - pkin(7) * t297 - t367 * t290 + t371 * t291;
t282 = Ifges(3,6) * t361 + t360 * Ifges(3,5) + mrSges(3,1) * g(1) + mrSges(3,3) * t337 + t367 * t291 + t371 * t290 - pkin(2) * (-m(4) * g(1) + t381) + pkin(7) * t378;
t281 = -mrSges(2,2) * g(1) - mrSges(2,3) * t354 + Ifges(2,5) * qJDD(1) - t374 * Ifges(2,6) - pkin(6) * t289 - t368 * t282 + t372 * t283;
t280 = Ifges(2,6) * qJDD(1) + t374 * Ifges(2,5) + mrSges(2,3) * t353 + t368 * t283 + t372 * t282 - pkin(1) * t381 + pkin(6) * t379 + (-pkin(1) * t386 + mrSges(2,1)) * g(1);
t1 = [(-m(1) - m(2) + t386) * g(1) + t381; -m(1) * g(2) + t383; -m(1) * g(3) + t380; pkin(1) * t289 + pkin(2) * t297 - mrSges(3,2) * t337 + mrSges(3,1) * t336 + pkin(3) * t305 + mrSges(4,1) * t329 - mrSges(4,2) * t330 + pkin(4) * t375 + pkin(8) * t376 + mrSges(5,1) * t323 - mrSges(5,2) * t324 + t366 * t311 + t370 * t310 + mrSges(2,1) * t354 - mrSges(2,2) * t353 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * t361 - mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + (Ifges(5,3) + Ifges(4,3)) * t356; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t380 + t373 * t280 + t369 * t281; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t383 + t369 * t280 - t373 * t281;];
tauB = t1;
