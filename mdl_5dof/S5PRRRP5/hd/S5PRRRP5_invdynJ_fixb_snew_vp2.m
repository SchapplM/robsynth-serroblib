% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRRP5
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRRP5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP5_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP5_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP5_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP5_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:47:58
% EndTime: 2019-12-05 16:48:00
% DurationCPUTime: 0.93s
% Computational Cost: add. (3391->190), mult. (6868->234), div. (0->0), fcn. (4293->8), ass. (0->78)
t378 = Ifges(5,4) + Ifges(6,4);
t384 = Ifges(5,2) + Ifges(6,2);
t381 = Ifges(5,6) + Ifges(6,6);
t383 = Ifges(5,1) + Ifges(6,1);
t382 = Ifges(5,5) + Ifges(6,5);
t380 = Ifges(5,3) + Ifges(6,3);
t357 = sin(qJ(4));
t358 = sin(qJ(3));
t360 = cos(qJ(4));
t361 = cos(qJ(3));
t330 = (-t358 * t357 + t361 * t360) * qJD(2);
t331 = (t361 * t357 + t358 * t360) * qJD(2);
t352 = qJD(3) + qJD(4);
t379 = t384 * t330 + t378 * t331 + t381 * t352;
t355 = sin(pkin(8));
t356 = cos(pkin(8));
t343 = -t356 * g(1) - t355 * g(2);
t354 = -g(3) + qJDD(1);
t359 = sin(qJ(2));
t362 = cos(qJ(2));
t325 = t362 * t343 + t359 * t354;
t363 = qJD(2) ^ 2;
t318 = -t363 * pkin(2) + qJDD(2) * pkin(6) + t325;
t342 = -t355 * g(1) + t356 * g(2);
t302 = -t358 * t318 + t361 * t342;
t373 = qJD(2) * qJD(3);
t370 = t361 * t373;
t340 = t358 * qJDD(2) + t370;
t286 = (-t340 + t370) * pkin(7) + (t358 * t361 * t363 + qJDD(3)) * pkin(3) + t302;
t303 = t361 * t318 + t358 * t342;
t341 = t361 * qJDD(2) - t358 * t373;
t374 = t358 * qJD(2);
t346 = qJD(3) * pkin(3) - pkin(7) * t374;
t353 = t361 ^ 2;
t287 = -t353 * t363 * pkin(3) + t341 * pkin(7) - qJD(3) * t346 + t303;
t281 = t360 * t286 - t357 * t287;
t301 = t330 * qJD(4) + t360 * t340 + t357 * t341;
t313 = -t330 * mrSges(6,1) + t331 * mrSges(6,2);
t314 = -t330 * mrSges(5,1) + t331 * mrSges(5,2);
t320 = -t352 * mrSges(5,2) + t330 * mrSges(5,3);
t351 = qJDD(3) + qJDD(4);
t275 = -0.2e1 * qJD(5) * t331 + (t330 * t352 - t301) * qJ(5) + (t330 * t331 + t351) * pkin(4) + t281;
t319 = -t352 * mrSges(6,2) + t330 * mrSges(6,3);
t372 = m(6) * t275 + t351 * mrSges(6,1) + t352 * t319;
t266 = m(5) * t281 + t351 * mrSges(5,1) + t352 * t320 + (-t313 - t314) * t331 + (-mrSges(5,3) - mrSges(6,3)) * t301 + t372;
t282 = t357 * t286 + t360 * t287;
t300 = -t331 * qJD(4) - t357 * t340 + t360 * t341;
t322 = t352 * mrSges(6,1) - t331 * mrSges(6,3);
t323 = t352 * mrSges(5,1) - t331 * mrSges(5,3);
t321 = t352 * pkin(4) - t331 * qJ(5);
t326 = t330 ^ 2;
t277 = -t326 * pkin(4) + t300 * qJ(5) + 0.2e1 * qJD(5) * t330 - t352 * t321 + t282;
t371 = m(6) * t277 + t300 * mrSges(6,3) + t330 * t313;
t269 = m(5) * t282 + t300 * mrSges(5,3) + t330 * t314 + (-t322 - t323) * t352 + (-mrSges(5,2) - mrSges(6,2)) * t351 + t371;
t264 = t360 * t266 + t357 * t269;
t377 = -t381 * t330 - t382 * t331 - t380 * t352;
t376 = -t378 * t330 - t383 * t331 - t382 * t352;
t375 = qJD(2) * t361;
t339 = (-t361 * mrSges(4,1) + t358 * mrSges(4,2)) * qJD(2);
t344 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t374;
t345 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t375;
t368 = -t357 * t266 + t360 * t269;
t369 = -t358 * (m(4) * t302 + qJDD(3) * mrSges(4,1) - t340 * mrSges(4,3) + qJD(3) * t345 - t339 * t374 + t264) + t361 * (m(4) * t303 - qJDD(3) * mrSges(4,2) + t341 * mrSges(4,3) - qJD(3) * t344 + t339 * t375 + t368);
t324 = -t359 * t343 + t362 * t354;
t367 = -qJDD(2) * pkin(2) - t324;
t288 = -t341 * pkin(3) + t346 * t374 + (-pkin(7) * t353 - pkin(6)) * t363 + t367;
t279 = -t300 * pkin(4) - t326 * qJ(5) + t331 * t321 + qJDD(5) + t288;
t272 = m(6) * t279 - t300 * mrSges(6,1) + t301 * mrSges(6,2) - t330 * t319 + t331 * t322;
t366 = m(5) * t288 - t300 * mrSges(5,1) + t301 * mrSges(5,2) - t330 * t320 + t331 * t323 + t272;
t271 = -t301 * mrSges(6,3) - t331 * t313 + t372;
t365 = mrSges(5,1) * t281 + mrSges(6,1) * t275 - mrSges(5,2) * t282 - mrSges(6,2) * t277 + pkin(4) * t271 + t381 * t300 + t382 * t301 + t376 * t330 + t379 * t331 + t380 * t351;
t317 = -t363 * pkin(6) + t367;
t364 = -m(4) * t317 + t341 * mrSges(4,1) - t340 * mrSges(4,2) - t344 * t374 + t345 * t375 - t366;
t329 = Ifges(4,5) * qJD(3) + (t358 * Ifges(4,1) + t361 * Ifges(4,4)) * qJD(2);
t328 = Ifges(4,6) * qJD(3) + (t358 * Ifges(4,4) + t361 * Ifges(4,2)) * qJD(2);
t262 = mrSges(5,2) * t288 + mrSges(6,2) * t279 - mrSges(5,3) * t281 - mrSges(6,3) * t275 - qJ(5) * t271 + t378 * t300 + t383 * t301 - t377 * t330 + t382 * t351 - t379 * t352;
t260 = -mrSges(5,1) * t288 + mrSges(5,3) * t282 - mrSges(6,1) * t279 + mrSges(6,3) * t277 - pkin(4) * t272 + qJ(5) * t371 + (-qJ(5) * t322 - t376) * t352 + (-qJ(5) * mrSges(6,2) + t381) * t351 + t377 * t331 + t378 * t301 + t384 * t300;
t1 = [m(2) * t354 + t359 * (m(3) * t325 - t363 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t369) + t362 * (m(3) * t324 + qJDD(2) * mrSges(3,1) - t363 * mrSges(3,2) + t364); Ifges(3,3) * qJDD(2) + mrSges(3,1) * t324 - mrSges(3,2) * t325 + t358 * (mrSges(4,2) * t317 - mrSges(4,3) * t302 + Ifges(4,1) * t340 + Ifges(4,4) * t341 + Ifges(4,5) * qJDD(3) - pkin(7) * t264 - qJD(3) * t328 - t357 * t260 + t360 * t262) + t361 * (-mrSges(4,1) * t317 + mrSges(4,3) * t303 + Ifges(4,4) * t340 + Ifges(4,2) * t341 + Ifges(4,6) * qJDD(3) - pkin(3) * t366 + pkin(7) * t368 + qJD(3) * t329 + t360 * t260 + t357 * t262) + pkin(2) * t364 + pkin(6) * t369; (t358 * t328 - t361 * t329) * qJD(2) + Ifges(4,5) * t340 + Ifges(4,6) * t341 + mrSges(4,1) * t302 - mrSges(4,2) * t303 + pkin(3) * t264 + t365 + Ifges(4,3) * qJDD(3); t365; t272;];
tauJ = t1;
