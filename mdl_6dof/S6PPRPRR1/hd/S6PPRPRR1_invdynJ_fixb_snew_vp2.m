% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PPRPRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
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
% Datum: 2019-05-04 19:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PPRPRR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRPRR1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRPRR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_invdynJ_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRPRR1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRPRR1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRPRR1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:52:54
% EndTime: 2019-05-04 19:52:55
% DurationCPUTime: 1.29s
% Computational Cost: add. (11944->173), mult. (20606->230), div. (0->0), fcn. (16927->16), ass. (0->91)
t360 = sin(pkin(11));
t365 = cos(pkin(11));
t353 = -t365 * g(1) - t360 * g(2);
t359 = sin(pkin(12));
t364 = cos(pkin(12));
t352 = t360 * g(1) - t365 * g(2);
t357 = -g(3) + qJDD(1);
t362 = sin(pkin(6));
t367 = cos(pkin(6));
t379 = t352 * t367 + t357 * t362;
t324 = -t359 * t353 + t364 * t379;
t339 = -t362 * t352 + t367 * t357 + qJDD(2);
t361 = sin(pkin(7));
t366 = cos(pkin(7));
t394 = t324 * t366 + t339 * t361;
t381 = -t361 * t324 + t366 * t339;
t320 = qJDD(4) + t381;
t372 = cos(qJ(5));
t391 = t372 * t320;
t325 = t364 * t353 + t359 * t379;
t370 = sin(qJ(3));
t373 = cos(qJ(3));
t316 = -t370 * t325 + t373 * t394;
t314 = qJDD(3) * pkin(3) + t316;
t317 = t373 * t325 + t370 * t394;
t375 = qJD(3) ^ 2;
t315 = -t375 * pkin(3) + t317;
t358 = sin(pkin(13));
t363 = cos(pkin(13));
t310 = t358 * t314 + t363 * t315;
t308 = -t375 * pkin(4) + qJDD(3) * pkin(9) + t310;
t369 = sin(qJ(5));
t304 = t372 * t308 + t369 * t320;
t348 = (-t372 * mrSges(6,1) + t369 * mrSges(6,2)) * qJD(3);
t387 = qJD(3) * qJD(5);
t385 = t369 * t387;
t351 = t372 * qJDD(3) - t385;
t389 = qJD(3) * t369;
t354 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t389;
t349 = (-t372 * pkin(5) - t369 * pkin(10)) * qJD(3);
t374 = qJD(5) ^ 2;
t388 = t372 * qJD(3);
t302 = -t374 * pkin(5) + qJDD(5) * pkin(10) + t349 * t388 + t304;
t309 = t363 * t314 - t358 * t315;
t307 = -qJDD(3) * pkin(4) - t375 * pkin(9) - t309;
t384 = t372 * t387;
t350 = t369 * qJDD(3) + t384;
t305 = (-t350 - t384) * pkin(10) + (-t351 + t385) * pkin(5) + t307;
t368 = sin(qJ(6));
t371 = cos(qJ(6));
t299 = -t368 * t302 + t371 * t305;
t346 = t371 * qJD(5) - t368 * t389;
t332 = t346 * qJD(6) + t368 * qJDD(5) + t371 * t350;
t347 = t368 * qJD(5) + t371 * t389;
t333 = -t346 * mrSges(7,1) + t347 * mrSges(7,2);
t356 = qJD(6) - t388;
t337 = -t356 * mrSges(7,2) + t346 * mrSges(7,3);
t345 = qJDD(6) - t351;
t297 = m(7) * t299 + t345 * mrSges(7,1) - t332 * mrSges(7,3) - t347 * t333 + t356 * t337;
t300 = t371 * t302 + t368 * t305;
t331 = -t347 * qJD(6) + t371 * qJDD(5) - t368 * t350;
t338 = t356 * mrSges(7,1) - t347 * mrSges(7,3);
t298 = m(7) * t300 - t345 * mrSges(7,2) + t331 * mrSges(7,3) + t346 * t333 - t356 * t338;
t382 = -t368 * t297 + t371 * t298;
t290 = m(6) * t304 - qJDD(5) * mrSges(6,2) + t351 * mrSges(6,3) - qJD(5) * t354 + t348 * t388 + t382;
t303 = -t369 * t308 + t391;
t355 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t388;
t301 = -qJDD(5) * pkin(5) - t374 * pkin(10) - t391 + (qJD(3) * t349 + t308) * t369;
t378 = -m(7) * t301 + t331 * mrSges(7,1) - t332 * mrSges(7,2) + t346 * t337 - t347 * t338;
t295 = m(6) * t303 + qJDD(5) * mrSges(6,1) - t350 * mrSges(6,3) + qJD(5) * t355 - t348 * t389 + t378;
t383 = t372 * t290 - t369 * t295;
t284 = m(5) * t310 - t375 * mrSges(5,1) - qJDD(3) * mrSges(5,2) + t383;
t291 = t371 * t297 + t368 * t298;
t377 = -m(6) * t307 + t351 * mrSges(6,1) - t350 * mrSges(6,2) - t354 * t389 + t355 * t388 - t291;
t287 = m(5) * t309 + qJDD(3) * mrSges(5,1) - t375 * mrSges(5,2) + t377;
t390 = t358 * t284 + t363 * t287;
t386 = m(5) * t320 + t369 * t290 + t372 * t295;
t281 = m(4) * t316 + qJDD(3) * mrSges(4,1) - t375 * mrSges(4,2) + t390;
t282 = m(4) * t317 - t375 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t363 * t284 - t358 * t287;
t380 = t281 * t373 + t282 * t370;
t327 = Ifges(7,4) * t347 + Ifges(7,2) * t346 + Ifges(7,6) * t356;
t328 = Ifges(7,1) * t347 + Ifges(7,4) * t346 + Ifges(7,5) * t356;
t376 = mrSges(7,1) * t299 - mrSges(7,2) * t300 + Ifges(7,5) * t332 + Ifges(7,6) * t331 + Ifges(7,3) * t345 + t347 * t327 - t346 * t328;
t342 = Ifges(6,5) * qJD(5) + (t369 * Ifges(6,1) + t372 * Ifges(6,4)) * qJD(3);
t341 = Ifges(6,6) * qJD(5) + (t369 * Ifges(6,4) + t372 * Ifges(6,2)) * qJD(3);
t326 = Ifges(7,5) * t347 + Ifges(7,6) * t346 + Ifges(7,3) * t356;
t293 = mrSges(7,2) * t301 - mrSges(7,3) * t299 + Ifges(7,1) * t332 + Ifges(7,4) * t331 + Ifges(7,5) * t345 + t346 * t326 - t356 * t327;
t292 = -mrSges(7,1) * t301 + mrSges(7,3) * t300 + Ifges(7,4) * t332 + Ifges(7,2) * t331 + Ifges(7,6) * t345 - t347 * t326 + t356 * t328;
t285 = m(4) * t381 + t386;
t280 = m(3) * t339 + t366 * t285 + t380 * t361;
t1 = [m(2) * t357 + t367 * t280 + (t359 * (m(3) * t325 - t370 * t281 + t373 * t282) + t364 * (m(3) * t324 - t361 * t285 + t380 * t366)) * t362; t280; mrSges(4,1) * t316 - mrSges(4,2) * t317 + mrSges(5,1) * t309 - mrSges(5,2) * t310 + t369 * (mrSges(6,2) * t307 - mrSges(6,3) * t303 + Ifges(6,1) * t350 + Ifges(6,4) * t351 + Ifges(6,5) * qJDD(5) - pkin(10) * t291 - qJD(5) * t341 - t368 * t292 + t371 * t293) + t372 * (-mrSges(6,1) * t307 + mrSges(6,3) * t304 + Ifges(6,4) * t350 + Ifges(6,2) * t351 + Ifges(6,6) * qJDD(5) - pkin(5) * t291 + qJD(5) * t342 - t376) + pkin(4) * t377 + pkin(9) * t383 + pkin(3) * t390 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3); t386; Ifges(6,5) * t350 + Ifges(6,6) * t351 + Ifges(6,3) * qJDD(5) + mrSges(6,1) * t303 - mrSges(6,2) * t304 + t368 * t293 + t371 * t292 + pkin(5) * t378 + pkin(10) * t382 + (t369 * t341 - t372 * t342) * qJD(3); t376;];
tauJ  = t1;
