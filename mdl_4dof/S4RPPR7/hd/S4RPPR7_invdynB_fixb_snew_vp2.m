% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RPPR7
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
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RPPR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR7_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR7_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR7_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR7_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:38
% EndTime: 2019-12-31 16:41:39
% DurationCPUTime: 0.72s
% Computational Cost: add. (4234->166), mult. (8911->204), div. (0->0), fcn. (4946->6), ass. (0->74)
t350 = sin(qJ(1));
t352 = cos(qJ(1));
t332 = g(1) * t350 - t352 * g(2);
t353 = qJD(1) ^ 2;
t359 = -qJ(2) * t353 + qJDD(2) - t332;
t378 = -pkin(1) - qJ(3);
t384 = -(2 * qJD(1) * qJD(3)) + t378 * qJDD(1) + t359;
t347 = sin(pkin(6));
t341 = t347 ^ 2;
t348 = cos(pkin(6));
t375 = t348 ^ 2 + t341;
t370 = t375 * mrSges(4,3);
t333 = -g(1) * t352 - g(2) * t350;
t383 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t333;
t382 = pkin(3) * t353;
t381 = mrSges(2,1) - mrSges(3,2);
t380 = -Ifges(3,4) + Ifges(2,5);
t379 = -Ifges(2,6) + Ifges(3,5);
t377 = mrSges(4,2) * t348;
t315 = t347 * g(3) + t384 * t348;
t305 = (-pkin(5) * qJDD(1) - t347 * t382) * t348 + t315;
t316 = -g(3) * t348 + t384 * t347;
t372 = qJDD(1) * t347;
t306 = -pkin(5) * t372 - t341 * t382 + t316;
t349 = sin(qJ(4));
t351 = cos(qJ(4));
t303 = t305 * t351 - t306 * t349;
t363 = -t347 * t351 - t348 * t349;
t328 = t363 * qJD(1);
t362 = -t347 * t349 + t348 * t351;
t329 = t362 * qJD(1);
t313 = -mrSges(5,1) * t328 + mrSges(5,2) * t329;
t318 = qJD(4) * t328 + t362 * qJDD(1);
t323 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t328;
t301 = m(5) * t303 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t318 + qJD(4) * t323 - t313 * t329;
t304 = t305 * t349 + t306 * t351;
t317 = -qJD(4) * t329 + t363 * qJDD(1);
t324 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t329;
t302 = m(5) * t304 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t317 - qJD(4) * t324 + t313 * t328;
t292 = t351 * t301 + t349 * t302;
t361 = -qJDD(1) * mrSges(4,3) - t353 * (mrSges(4,1) * t347 + t377);
t290 = m(4) * t315 + t361 * t348 + t292;
t367 = -t301 * t349 + t351 * t302;
t291 = m(4) * t316 + t361 * t347 + t367;
t288 = t290 * t348 + t291 * t347;
t327 = -qJDD(1) * pkin(1) + t359;
t356 = -m(3) * t327 + t353 * mrSges(3,3) - t288;
t286 = m(2) * t332 - mrSges(2,2) * t353 + t381 * qJDD(1) + t356;
t326 = pkin(1) * t353 - t383;
t358 = qJDD(3) + t383;
t322 = t378 * t353 + t358;
t308 = pkin(3) * t372 + (-t375 * pkin(5) + t378) * t353 + t358;
t357 = m(5) * t308 - mrSges(5,1) * t317 + t318 * mrSges(5,2) - t323 * t328 + t329 * t324;
t355 = -m(4) * t322 - mrSges(4,1) * t372 - qJDD(1) * t377 - t357;
t354 = -m(3) * t326 + t353 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t355;
t297 = t354 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - t370) * t353 + m(2) * t333;
t376 = t352 * t286 + t350 * t297;
t364 = Ifges(4,5) * t348 - Ifges(4,6) * t347;
t374 = t353 * t364;
t369 = -t286 * t350 + t352 * t297;
t368 = -t347 * t290 + t348 * t291;
t366 = Ifges(4,1) * t348 - Ifges(4,4) * t347;
t365 = Ifges(4,4) * t348 - Ifges(4,2) * t347;
t311 = Ifges(5,1) * t329 + Ifges(5,4) * t328 + Ifges(5,5) * qJD(4);
t310 = Ifges(5,4) * t329 + Ifges(5,2) * t328 + Ifges(5,6) * qJD(4);
t309 = Ifges(5,5) * t329 + Ifges(5,6) * t328 + Ifges(5,3) * qJD(4);
t294 = mrSges(5,2) * t308 - mrSges(5,3) * t303 + Ifges(5,1) * t318 + Ifges(5,4) * t317 + Ifges(5,5) * qJDD(4) - qJD(4) * t310 + t309 * t328;
t293 = -mrSges(5,1) * t308 + mrSges(5,3) * t304 + Ifges(5,4) * t318 + Ifges(5,2) * t317 + Ifges(5,6) * qJDD(4) + qJD(4) * t311 - t309 * t329;
t287 = -m(3) * g(3) + t368;
t284 = mrSges(4,2) * t322 - mrSges(4,3) * t315 - pkin(5) * t292 + t366 * qJDD(1) - t293 * t349 + t294 * t351 - t347 * t374;
t283 = -mrSges(4,1) * t322 + mrSges(4,3) * t316 - pkin(3) * t357 + pkin(5) * t367 + t365 * qJDD(1) + t351 * t293 + t349 * t294 - t348 * t374;
t282 = mrSges(3,1) * t327 + mrSges(4,1) * t315 + mrSges(5,1) * t303 - mrSges(4,2) * t316 - mrSges(5,2) * t304 - mrSges(2,3) * t332 + Ifges(5,5) * t318 + Ifges(5,6) * t317 + Ifges(5,3) * qJDD(4) + pkin(2) * t288 + pkin(3) * t292 - qJ(2) * t287 + t329 * t310 - t328 * t311 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (t364 + t380) * qJDD(1) + (t347 * t366 + t348 * t365 + t379) * t353;
t281 = mrSges(2,3) * t333 - mrSges(3,1) * t326 - t347 * t284 - t348 * t283 - pkin(2) * t355 - qJ(3) * t368 - pkin(1) * t287 - t379 * qJDD(1) + t381 * g(3) + (-pkin(2) * t370 + t380) * t353;
t1 = [-m(1) * g(1) + t369; -m(1) * g(2) + t376; (-m(1) - m(2) - m(3)) * g(3) + t368; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t376 - t350 * t281 + t352 * t282; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t369 + t352 * t281 + t350 * t282; pkin(1) * t356 + qJ(2) * (-t353 * t370 + t354) + t348 * t284 - t347 * t283 - qJ(3) * t288 + mrSges(2,1) * t332 - mrSges(2,2) * t333 + mrSges(3,2) * t327 - mrSges(3,3) * t326 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB = t1;
