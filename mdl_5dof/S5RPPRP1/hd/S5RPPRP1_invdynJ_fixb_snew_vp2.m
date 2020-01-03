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
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:25:21
% EndTime: 2020-01-03 11:25:25
% DurationCPUTime: 1.38s
% Computational Cost: add. (2323->164), mult. (4820->209), div. (0->0), fcn. (2671->8), ass. (0->80)
t345 = sin(qJ(1));
t347 = cos(qJ(1));
t353 = -g(2) * t347 - g(3) * t345;
t325 = qJDD(1) * pkin(1) + t353;
t348 = qJD(1) ^ 2;
t359 = -g(2) * t345 + t347 * g(3);
t326 = -pkin(1) * t348 + t359;
t341 = sin(pkin(7));
t343 = cos(pkin(7));
t300 = t341 * t325 + t343 * t326;
t390 = -pkin(2) * t348 + qJDD(1) * qJ(3) + (2 * qJD(1) * qJD(3)) + t300;
t344 = sin(qJ(4));
t346 = cos(qJ(4));
t380 = Ifges(5,4) + Ifges(6,4);
t389 = t346 * (Ifges(5,1) + Ifges(6,1)) - t344 * t380;
t388 = t346 * t380 - t344 * (Ifges(5,2) + Ifges(6,2));
t379 = Ifges(5,5) + Ifges(6,5);
t378 = Ifges(5,6) + Ifges(6,6);
t342 = cos(pkin(8));
t369 = qJD(1) * t342;
t329 = qJD(4) - t369;
t340 = sin(pkin(8));
t370 = qJD(1) * t340;
t358 = (t378 * t329 + t388 * t370) * t346;
t372 = t379 * t329 + t389 * t370;
t386 = t372 * t344 + t358;
t339 = -g(1) + qJDD(2);
t290 = t339 * t342 - t390 * t340;
t314 = (t344 * mrSges(6,1) + t346 * mrSges(6,2)) * t370;
t366 = qJD(1) * qJD(4);
t317 = (qJDD(1) * t346 - t344 * t366) * t340;
t361 = t346 * t370;
t291 = t340 * t339 + t390 * t342;
t354 = -t342 * pkin(3) - t340 * pkin(6);
t324 = t354 * qJD(1);
t289 = t324 * t369 + t291;
t299 = t325 * t343 - t341 * t326;
t351 = -qJ(3) * t348 + qJDD(3) - t299;
t294 = (-pkin(2) + t354) * qJDD(1) + t351;
t293 = t346 * t294;
t365 = qJDD(1) * t342;
t328 = qJDD(4) - t365;
t355 = -0.2e1 * qJD(5) * t370;
t375 = t340 ^ 2 * t348;
t281 = t346 * t355 + pkin(4) * t328 - qJ(5) * t317 + t293 + (-pkin(4) * t346 * t375 - qJ(5) * t329 * t370 - t289) * t344;
t362 = t344 * t370;
t309 = -mrSges(6,2) * t329 - mrSges(6,3) * t362;
t363 = m(6) * t281 + t328 * mrSges(6,1) + t329 * t309;
t278 = -mrSges(6,3) * t317 - t314 * t361 + t363;
t285 = t346 * t289 + t344 * t294;
t311 = pkin(4) * t329 - qJ(5) * t361;
t316 = (-qJDD(1) * t344 - t346 * t366) * t340;
t364 = t344 ^ 2 * t375;
t283 = -pkin(4) * t364 + qJ(5) * t316 - t311 * t329 + t344 * t355 + t285;
t284 = -t289 * t344 + t293;
t383 = mrSges(5,1) * t284 + mrSges(6,1) * t281 - mrSges(5,2) * t285 - mrSges(6,2) * t283 + pkin(4) * t278 + t378 * t316 + t379 * t317 + (Ifges(5,3) + Ifges(6,3)) * t328;
t288 = t324 * t370 - t290;
t286 = -pkin(4) * t316 - qJ(5) * t364 + t311 * t361 + qJDD(5) + t288;
t382 = m(6) * t286;
t381 = -mrSges(5,2) - mrSges(6,2);
t376 = t340 * mrSges(4,2);
t374 = m(6) * t283 + t316 * mrSges(6,3);
t312 = mrSges(6,1) * t329 - mrSges(6,3) * t361;
t371 = -mrSges(5,1) * t329 + mrSges(5,3) * t361 - t312;
t368 = qJDD(1) * mrSges(4,3);
t310 = -mrSges(5,2) * t329 - mrSges(5,3) * t362;
t352 = (-t314 - (t344 * mrSges(5,1) + t346 * mrSges(5,2)) * t370) * t370;
t275 = m(5) * t284 + mrSges(5,1) * t328 + t310 * t329 + (-mrSges(5,3) - mrSges(6,3)) * t317 + t346 * t352 + t363;
t276 = m(5) * t285 + mrSges(5,3) * t316 + t381 * t328 + t371 * t329 + t344 * t352 + t374;
t322 = (-t342 * mrSges(4,1) + t376) * qJD(1);
t272 = m(4) * t291 - t275 * t344 + t276 * t346 + (qJD(1) * t322 + t368) * t342;
t277 = m(4) * t290 - m(5) * t288 - t382 + t381 * t317 + (mrSges(5,1) + mrSges(6,1)) * t316 + (-t368 + (-t322 + t371 * t346 + (-t309 - t310) * t344) * qJD(1)) * t340;
t357 = t342 * t272 - t277 * t340;
t274 = t275 * t346 + t276 * t344;
t297 = -qJDD(1) * pkin(2) + t351;
t350 = -m(4) * t297 + mrSges(4,1) * t365 - t274 + (t342 ^ 2 * t348 + t375) * mrSges(4,3);
t323 = (Ifges(4,5) * t340 + Ifges(4,6) * t342) * qJD(1);
t279 = t382 - mrSges(6,1) * t316 + mrSges(6,2) * t317 + (t309 * t344 + t312 * t346) * t370;
t273 = qJDD(1) * t376 - t350;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t353 - mrSges(2,2) * t359 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t299 - mrSges(3,2) * t300 + t340 * (t323 * t369 + mrSges(4,2) * t297 - mrSges(4,3) * t290 + t346 * (mrSges(5,2) * t288 + mrSges(6,2) * t286 - mrSges(5,3) * t284 - mrSges(6,3) * t281 - qJ(5) * t278) - t344 * (-mrSges(5,1) * t288 + mrSges(5,3) * t285 - mrSges(6,1) * t286 + mrSges(6,3) * t283 - pkin(4) * t279 + qJ(5) * (-t314 * t362 + t374)) - pkin(6) * t274 + (-t358 - t344 * (-qJ(5) * t312 + t372)) * t329 + (t346 * t379 - t344 * (-qJ(5) * mrSges(6,2) + t378)) * t328 + t389 * t317 + t388 * t316 + (Ifges(4,1) * t340 + Ifges(4,4) * t342) * qJDD(1)) + t342 * (Ifges(4,2) * t365 - mrSges(4,1) * t297 + mrSges(4,3) * t291 - pkin(3) * t274 + (Ifges(4,4) * qJDD(1) + (-t323 - t386) * qJD(1)) * t340 - t383) - pkin(2) * t273 + qJ(3) * t357 + pkin(1) * (t341 * (m(3) * t300 - mrSges(3,1) * t348 - qJDD(1) * mrSges(3,2) + t357) + t343 * (m(3) * t299 - mrSges(3,2) * t348 + (mrSges(3,1) - t376) * qJDD(1) + t350)); m(3) * t339 + t272 * t340 + t277 * t342; t273; t370 * t386 + t383; t279;];
tauJ = t1;
