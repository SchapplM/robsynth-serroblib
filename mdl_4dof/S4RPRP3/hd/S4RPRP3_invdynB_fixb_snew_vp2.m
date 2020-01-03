% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RPRP3
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RPRP3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP3_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP3_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP3_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP3_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:36
% EndTime: 2019-12-31 16:42:37
% DurationCPUTime: 0.81s
% Computational Cost: add. (4328->177), mult. (8185->213), div. (0->0), fcn. (3674->6), ass. (0->75)
t385 = Ifges(4,1) + Ifges(5,1);
t378 = Ifges(4,4) + Ifges(5,4);
t377 = Ifges(4,5) + Ifges(5,5);
t384 = Ifges(4,2) + Ifges(5,2);
t383 = Ifges(4,6) + Ifges(5,6);
t382 = Ifges(4,3) + Ifges(5,3);
t356 = qJD(1) ^ 2;
t353 = sin(qJ(1));
t355 = cos(qJ(1));
t341 = t353 * g(1) - t355 * g(2);
t329 = qJDD(1) * pkin(1) + t341;
t342 = -t355 * g(1) - t353 * g(2);
t332 = -t356 * pkin(1) + t342;
t350 = sin(pkin(6));
t351 = cos(pkin(6));
t312 = t351 * t329 - t350 * t332;
t358 = -qJDD(1) * pkin(2) - t312;
t310 = -t356 * pkin(5) + t358;
t352 = sin(qJ(3));
t354 = cos(qJ(3));
t368 = qJD(1) * qJD(3);
t363 = t354 * t368;
t333 = t352 * qJDD(1) + t363;
t334 = t354 * qJDD(1) - t352 * t368;
t369 = qJD(1) * t354;
t340 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t369;
t370 = qJD(1) * t352;
t336 = qJD(3) * pkin(3) - qJ(4) * t370;
t348 = t354 ^ 2;
t306 = t336 * t370 - t334 * pkin(3) + qJDD(4) + (-qJ(4) * t348 - pkin(5)) * t356 + t358;
t339 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t369;
t359 = m(5) * t306 - t334 * mrSges(5,1) - t339 * t369;
t379 = -mrSges(4,2) - mrSges(5,2);
t381 = -m(4) * t310 + t334 * mrSges(4,1) + t379 * t333 + t340 * t369 - t359;
t380 = pkin(3) * t356;
t313 = t350 * t329 + t351 * t332;
t311 = -t356 * pkin(2) + qJDD(1) * pkin(5) + t313;
t349 = -g(3) + qJDD(2);
t308 = t354 * t311 + t352 * t349;
t331 = (-mrSges(4,1) * t354 + mrSges(4,2) * t352) * qJD(1);
t367 = qJD(1) * qJD(4);
t305 = t334 * qJ(4) - qJD(3) * t336 - t348 * t380 + 0.2e1 * t354 * t367 + t308;
t330 = (-mrSges(5,1) * t354 + mrSges(5,2) * t352) * qJD(1);
t365 = m(5) * t305 + t334 * mrSges(5,3) + t330 * t369;
t337 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t370;
t371 = -qJD(3) * mrSges(4,1) + mrSges(4,3) * t370 - t337;
t300 = m(4) * t308 + t334 * mrSges(4,3) + t371 * qJD(3) + t379 * qJDD(3) + t331 * t369 + t365;
t298 = t354 * t300;
t344 = t354 * t349;
t307 = -t352 * t311 + t344;
t304 = qJDD(3) * pkin(3) + t344 + (-t333 + t363) * qJ(4) + (t354 * t380 - t311 - 0.2e1 * t367) * t352;
t366 = m(5) * t304 + qJDD(3) * mrSges(5,1) + qJD(3) * t339;
t299 = m(4) * t307 + qJDD(3) * mrSges(4,1) + qJD(3) * t340 + (-mrSges(4,3) - mrSges(5,3)) * t333 + (-t330 - t331) * t370 + t366;
t292 = m(3) * t313 - t356 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t352 * t299 + t298;
t360 = qJD(1) * t371;
t295 = m(3) * t312 + qJDD(1) * mrSges(3,1) - t356 * mrSges(3,2) + t352 * t360 + t381;
t287 = t350 * t292 + t351 * t295;
t285 = m(2) * t341 + qJDD(1) * mrSges(2,1) - t356 * mrSges(2,2) + t287;
t361 = t351 * t292 - t350 * t295;
t286 = m(2) * t342 - t356 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t361;
t375 = t355 * t285 + t353 * t286;
t293 = t354 * t299 + t352 * t300;
t374 = t382 * qJD(3) + (t377 * t352 + t383 * t354) * qJD(1);
t373 = -t383 * qJD(3) + (-t378 * t352 - t384 * t354) * qJD(1);
t372 = t377 * qJD(3) + (t385 * t352 + t378 * t354) * qJD(1);
t364 = m(3) * t349 + t293;
t362 = -t353 * t285 + t355 * t286;
t301 = -t333 * mrSges(5,3) - t330 * t370 + t366;
t289 = mrSges(4,2) * t310 + mrSges(5,2) * t306 - mrSges(4,3) * t307 - mrSges(5,3) * t304 - qJ(4) * t301 + t373 * qJD(3) + t377 * qJDD(3) + t385 * t333 + t378 * t334 + t374 * t369;
t288 = -mrSges(4,1) * t310 + mrSges(4,3) * t308 - mrSges(5,1) * t306 + mrSges(5,3) * t305 - pkin(3) * t359 + qJ(4) * t365 + t384 * t334 + (-pkin(3) * mrSges(5,2) + t378) * t333 + (-qJ(4) * mrSges(5,2) + t383) * qJDD(3) + (-qJ(4) * t337 + t372) * qJD(3) + (-pkin(3) * t337 - t374) * t370;
t281 = -mrSges(3,1) * t349 - mrSges(4,1) * t307 - mrSges(5,1) * t304 + mrSges(4,2) * t308 + mrSges(5,2) * t305 + mrSges(3,3) * t313 + t356 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t293 - pkin(3) * t301 - t383 * t334 - t377 * t333 - t382 * qJDD(3) + (t373 * t352 + t372 * t354) * qJD(1);
t280 = mrSges(3,2) * t349 - mrSges(3,3) * t312 + Ifges(3,5) * qJDD(1) - t356 * Ifges(3,6) - pkin(5) * t293 - t352 * t288 + t354 * t289;
t279 = -mrSges(2,2) * g(3) - mrSges(2,3) * t341 + Ifges(2,5) * qJDD(1) - t356 * Ifges(2,6) - qJ(2) * t287 + t351 * t280 - t350 * t281;
t278 = mrSges(2,1) * g(3) + mrSges(2,3) * t342 + t356 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t364 + qJ(2) * t361 + t350 * t280 + t351 * t281;
t1 = [-m(1) * g(1) + t362; -m(1) * g(2) + t375; (-m(1) - m(2)) * g(3) + t364; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t375 - t353 * t278 + t355 * t279; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t362 + t355 * t278 + t353 * t279; pkin(1) * t287 + mrSges(2,1) * t341 - mrSges(2,2) * t342 + t354 * t288 + pkin(2) * t381 + pkin(5) * t298 + mrSges(3,1) * t312 - mrSges(3,2) * t313 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (pkin(2) * t360 - pkin(5) * t299 + t289) * t352 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB = t1;
