% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% m [6x1]
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
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPRP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP1_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP1_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:21
% EndTime: 2022-01-23 09:12:22
% DurationCPUTime: 1.31s
% Computational Cost: add. (2323->164), mult. (4820->209), div. (0->0), fcn. (2671->8), ass. (0->80)
t346 = sin(qJ(1));
t348 = cos(qJ(1));
t360 = t346 * g(1) - g(2) * t348;
t326 = qJDD(1) * pkin(1) + t360;
t349 = qJD(1) ^ 2;
t354 = -g(1) * t348 - g(2) * t346;
t327 = -pkin(1) * t349 + t354;
t342 = sin(pkin(7));
t344 = cos(pkin(7));
t301 = t342 * t326 + t344 * t327;
t391 = -pkin(2) * t349 + qJDD(1) * qJ(3) + (2 * qJD(1) * qJD(3)) + t301;
t345 = sin(qJ(4));
t347 = cos(qJ(4));
t381 = Ifges(5,4) + Ifges(6,4);
t390 = t347 * (Ifges(5,1) + Ifges(6,1)) - t345 * t381;
t389 = t347 * t381 - t345 * (Ifges(5,2) + Ifges(6,2));
t380 = Ifges(5,5) + Ifges(6,5);
t379 = Ifges(5,6) + Ifges(6,6);
t343 = cos(pkin(8));
t369 = t343 * qJD(1);
t330 = qJD(4) - t369;
t341 = sin(pkin(8));
t370 = t341 * qJD(1);
t359 = (t379 * t330 + t389 * t370) * t347;
t373 = t380 * t330 + t390 * t370;
t387 = t373 * t345 + t359;
t340 = -g(3) + qJDD(2);
t291 = t340 * t343 - t391 * t341;
t315 = (t345 * mrSges(6,1) + t347 * mrSges(6,2)) * t370;
t367 = qJD(1) * qJD(4);
t318 = (qJDD(1) * t347 - t345 * t367) * t341;
t362 = t347 * t370;
t292 = t341 * t340 + t391 * t343;
t355 = -t343 * pkin(3) - t341 * pkin(6);
t325 = t355 * qJD(1);
t290 = t325 * t369 + t292;
t300 = t326 * t344 - t342 * t327;
t352 = -qJ(3) * t349 + qJDD(3) - t300;
t295 = (-pkin(2) + t355) * qJDD(1) + t352;
t294 = t347 * t295;
t366 = t343 * qJDD(1);
t329 = qJDD(4) - t366;
t356 = -0.2e1 * qJD(5) * t370;
t376 = t341 ^ 2 * t349;
t282 = t347 * t356 + pkin(4) * t329 - qJ(5) * t318 + t294 + (-pkin(4) * t347 * t376 - qJ(5) * t330 * t370 - t290) * t345;
t363 = t345 * t370;
t310 = -mrSges(6,2) * t330 - mrSges(6,3) * t363;
t364 = m(6) * t282 + t329 * mrSges(6,1) + t330 * t310;
t279 = -mrSges(6,3) * t318 - t315 * t362 + t364;
t286 = t347 * t290 + t345 * t295;
t312 = pkin(4) * t330 - qJ(5) * t362;
t317 = (-qJDD(1) * t345 - t347 * t367) * t341;
t365 = t345 ^ 2 * t376;
t284 = -pkin(4) * t365 + qJ(5) * t317 - t312 * t330 + t345 * t356 + t286;
t285 = -t290 * t345 + t294;
t384 = mrSges(5,1) * t285 + mrSges(6,1) * t282 - mrSges(5,2) * t286 - mrSges(6,2) * t284 + pkin(4) * t279 + t379 * t317 + t380 * t318 + (Ifges(5,3) + Ifges(6,3)) * t329;
t289 = t325 * t370 - t291;
t287 = -pkin(4) * t317 - qJ(5) * t365 + t312 * t362 + qJDD(5) + t289;
t383 = m(6) * t287;
t382 = -mrSges(5,2) - mrSges(6,2);
t377 = mrSges(4,2) * t341;
t375 = m(6) * t284 + t317 * mrSges(6,3);
t313 = mrSges(6,1) * t330 - mrSges(6,3) * t362;
t372 = -mrSges(5,1) * t330 + mrSges(5,3) * t362 - t313;
t371 = qJDD(1) * mrSges(4,3);
t311 = -mrSges(5,2) * t330 - mrSges(5,3) * t363;
t353 = (-t315 - (t345 * mrSges(5,1) + t347 * mrSges(5,2)) * t370) * t370;
t276 = m(5) * t285 + mrSges(5,1) * t329 + t311 * t330 + (-mrSges(5,3) - mrSges(6,3)) * t318 + t347 * t353 + t364;
t277 = m(5) * t286 + mrSges(5,3) * t317 + t382 * t329 + t372 * t330 + t345 * t353 + t375;
t323 = (-mrSges(4,1) * t343 + t377) * qJD(1);
t273 = m(4) * t292 - t276 * t345 + t277 * t347 + (qJD(1) * t323 + t371) * t343;
t278 = m(4) * t291 - m(5) * t289 - t383 + t382 * t318 + (mrSges(5,1) + mrSges(6,1)) * t317 + (-t371 + (-t323 + t372 * t347 + (-t310 - t311) * t345) * qJD(1)) * t341;
t358 = t343 * t273 - t278 * t341;
t275 = t276 * t347 + t277 * t345;
t298 = -qJDD(1) * pkin(2) + t352;
t351 = -m(4) * t298 + mrSges(4,1) * t366 - t275 + (t343 ^ 2 * t349 + t376) * mrSges(4,3);
t324 = (Ifges(4,5) * t341 + Ifges(4,6) * t343) * qJD(1);
t280 = t383 - mrSges(6,1) * t317 + mrSges(6,2) * t318 + (t310 * t345 + t313 * t347) * t370;
t274 = qJDD(1) * t377 - t351;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t360 - mrSges(2,2) * t354 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t300 - mrSges(3,2) * t301 + t341 * (t324 * t369 + mrSges(4,2) * t298 - mrSges(4,3) * t291 + t347 * (mrSges(5,2) * t289 + mrSges(6,2) * t287 - mrSges(5,3) * t285 - mrSges(6,3) * t282 - qJ(5) * t279) - t345 * (-mrSges(5,1) * t289 + mrSges(5,3) * t286 - mrSges(6,1) * t287 + mrSges(6,3) * t284 - pkin(4) * t280 + qJ(5) * (-t315 * t363 + t375)) - pkin(6) * t275 + (-t359 - t345 * (-qJ(5) * t313 + t373)) * t330 + (t347 * t380 - t345 * (-qJ(5) * mrSges(6,2) + t379)) * t329 + t390 * t318 + t389 * t317 + (Ifges(4,1) * t341 + Ifges(4,4) * t343) * qJDD(1)) + t343 * (Ifges(4,2) * t366 - mrSges(4,1) * t298 + mrSges(4,3) * t292 - pkin(3) * t275 + (Ifges(4,4) * qJDD(1) + (-t324 - t387) * qJD(1)) * t341 - t384) - pkin(2) * t274 + qJ(3) * t358 + pkin(1) * (t342 * (m(3) * t301 - mrSges(3,1) * t349 - qJDD(1) * mrSges(3,2) + t358) + t344 * (m(3) * t300 - mrSges(3,2) * t349 + (mrSges(3,1) - t377) * qJDD(1) + t351)); m(3) * t340 + t273 * t341 + t278 * t343; t274; t387 * t370 + t384; t280;];
tauJ = t1;
