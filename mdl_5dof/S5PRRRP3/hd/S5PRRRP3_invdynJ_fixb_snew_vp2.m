% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRRP3
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
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRRP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP3_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP3_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP3_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP3_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:43:44
% EndTime: 2019-12-05 16:43:46
% DurationCPUTime: 0.96s
% Computational Cost: add. (3287->190), mult. (6692->231), div. (0->0), fcn. (4199->8), ass. (0->80)
t379 = Ifges(5,4) + Ifges(6,4);
t386 = Ifges(5,2) + Ifges(6,2);
t383 = Ifges(5,6) + Ifges(6,6);
t385 = Ifges(5,1) + Ifges(6,1);
t384 = Ifges(5,5) + Ifges(6,5);
t382 = Ifges(5,3) + Ifges(6,3);
t355 = sin(qJ(4));
t356 = sin(qJ(3));
t358 = cos(qJ(4));
t359 = cos(qJ(3));
t330 = (-t355 * t356 + t358 * t359) * qJD(2);
t331 = (t355 * t359 + t356 * t358) * qJD(2);
t350 = qJD(3) + qJD(4);
t381 = t386 * t330 + t379 * t331 + t383 * t350;
t371 = qJD(2) * qJD(3);
t367 = t359 * t371;
t337 = qJDD(2) * t356 + t367;
t338 = qJDD(2) * t359 - t356 * t371;
t304 = -qJD(4) * t331 - t337 * t355 + t338 * t358;
t321 = -mrSges(6,2) * t350 + mrSges(6,3) * t330;
t380 = -mrSges(6,1) * t304 - t321 * t330;
t361 = qJD(2) ^ 2;
t353 = sin(pkin(8));
t354 = cos(pkin(8));
t339 = g(1) * t353 - g(2) * t354;
t340 = -g(1) * t354 - g(2) * t353;
t357 = sin(qJ(2));
t360 = cos(qJ(2));
t374 = t357 * t339 + t360 * t340;
t319 = -pkin(2) * t361 + qJDD(2) * pkin(6) + t374;
t352 = -g(3) + qJDD(1);
t306 = -t319 * t356 + t359 * t352;
t292 = (-t337 + t367) * pkin(7) + (t356 * t359 * t361 + qJDD(3)) * pkin(3) + t306;
t307 = t359 * t319 + t356 * t352;
t373 = qJD(2) * t356;
t343 = qJD(3) * pkin(3) - pkin(7) * t373;
t351 = t359 ^ 2;
t293 = -pkin(3) * t351 * t361 + pkin(7) * t338 - qJD(3) * t343 + t307;
t287 = t358 * t292 - t293 * t355;
t305 = qJD(4) * t330 + t337 * t358 + t338 * t355;
t316 = -mrSges(6,1) * t330 + mrSges(6,2) * t331;
t317 = -mrSges(5,1) * t330 + mrSges(5,2) * t331;
t322 = -mrSges(5,2) * t350 + mrSges(5,3) * t330;
t349 = qJDD(3) + qJDD(4);
t281 = -0.2e1 * qJD(5) * t331 + (t330 * t350 - t305) * qJ(5) + (t330 * t331 + t349) * pkin(4) + t287;
t370 = m(6) * t281 + t349 * mrSges(6,1) + t350 * t321;
t272 = m(5) * t287 + mrSges(5,1) * t349 + t322 * t350 + (-t316 - t317) * t331 + (-mrSges(5,3) - mrSges(6,3)) * t305 + t370;
t288 = t355 * t292 + t358 * t293;
t324 = mrSges(6,1) * t350 - mrSges(6,3) * t331;
t325 = mrSges(5,1) * t350 - mrSges(5,3) * t331;
t323 = pkin(4) * t350 - qJ(5) * t331;
t326 = t330 ^ 2;
t283 = -pkin(4) * t326 + qJ(5) * t304 + 0.2e1 * qJD(5) * t330 - t323 * t350 + t288;
t369 = m(6) * t283 + t304 * mrSges(6,3) + t330 * t316;
t275 = m(5) * t288 + mrSges(5,3) * t304 + t317 * t330 + (-t324 - t325) * t350 + (-mrSges(5,2) - mrSges(6,2)) * t349 + t369;
t270 = t358 * t272 + t355 * t275;
t376 = -t383 * t330 - t384 * t331 - t382 * t350;
t375 = -t379 * t330 - t385 * t331 - t384 * t350;
t372 = qJD(2) * t359;
t365 = t339 * t360 - t357 * t340;
t364 = -qJDD(2) * pkin(2) - t365;
t294 = -pkin(3) * t338 + t343 * t373 + (-pkin(7) * t351 - pkin(6)) * t361 + t364;
t285 = -pkin(4) * t304 - qJ(5) * t326 + t323 * t331 + qJDD(5) + t294;
t368 = m(6) * t285 + t305 * mrSges(6,2) + t331 * t324;
t366 = -t272 * t355 + t358 * t275;
t363 = m(5) * t294 + mrSges(5,2) * t305 + t325 * t331 + t368;
t277 = -mrSges(6,3) * t305 - t316 * t331 + t370;
t362 = mrSges(5,1) * t287 + mrSges(6,1) * t281 - mrSges(5,2) * t288 - mrSges(6,2) * t283 + pkin(4) * t277 + t383 * t304 + t384 * t305 + t375 * t330 + t381 * t331 + t382 * t349;
t342 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t372;
t341 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t373;
t336 = (-mrSges(4,1) * t359 + mrSges(4,2) * t356) * qJD(2);
t329 = Ifges(4,5) * qJD(3) + (t356 * Ifges(4,1) + t359 * Ifges(4,4)) * qJD(2);
t328 = Ifges(4,6) * qJD(3) + (t356 * Ifges(4,4) + t359 * Ifges(4,2)) * qJD(2);
t318 = -pkin(6) * t361 + t364;
t278 = t368 + t380;
t269 = m(4) * t307 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t338 - qJD(3) * t341 + t336 * t372 + t366;
t268 = m(4) * t306 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t337 + qJD(3) * t342 - t336 * t373 + t270;
t267 = mrSges(5,2) * t294 + mrSges(6,2) * t285 - mrSges(5,3) * t287 - mrSges(6,3) * t281 - qJ(5) * t277 + t379 * t304 + t385 * t305 - t376 * t330 + t384 * t349 - t381 * t350;
t266 = -mrSges(5,1) * t294 + mrSges(5,3) * t288 - mrSges(6,1) * t285 + mrSges(6,3) * t283 - pkin(4) * t278 + qJ(5) * t369 + (-qJ(5) * t324 - t375) * t350 + (-qJ(5) * mrSges(6,2) + t383) * t349 + t376 * t331 + t379 * t305 + t386 * t304;
t1 = [t268 * t359 + t269 * t356 + (m(2) + m(3)) * t352; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t365 - mrSges(3,2) * t374 + t356 * (mrSges(4,2) * t318 - mrSges(4,3) * t306 + Ifges(4,1) * t337 + Ifges(4,4) * t338 + Ifges(4,5) * qJDD(3) - pkin(7) * t270 - qJD(3) * t328 - t266 * t355 + t267 * t358) + t359 * (Ifges(4,4) * t337 + Ifges(4,2) * t338 + Ifges(4,6) * qJDD(3) + qJD(3) * t329 - mrSges(4,1) * t318 + mrSges(4,3) * t307 + t355 * t267 + t358 * t266 - pkin(3) * ((-t321 - t322) * t330 + (-mrSges(5,1) - mrSges(6,1)) * t304 + t363) + pkin(7) * t366) + pkin(6) * (-t268 * t356 + t359 * t269) + (-m(4) * t318 + mrSges(4,1) * t338 + mrSges(5,1) * t304 - mrSges(4,2) * t337 + t322 * t330 - t363 + (-t341 * t356 + t342 * t359) * qJD(2) - t380) * pkin(2); Ifges(4,5) * t337 + Ifges(4,6) * t338 + mrSges(4,1) * t306 - mrSges(4,2) * t307 + pkin(3) * t270 + t362 + Ifges(4,3) * qJDD(3) + (t356 * t328 - t359 * t329) * qJD(2); t362; t278;];
tauJ = t1;
