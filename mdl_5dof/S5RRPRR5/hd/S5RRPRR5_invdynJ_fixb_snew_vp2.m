% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPRR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR5_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR5_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:03:45
% EndTime: 2020-01-03 12:03:46
% DurationCPUTime: 1.32s
% Computational Cost: add. (17827->191), mult. (24620->246), div. (0->0), fcn. (16678->10), ass. (0->90)
t365 = qJD(1) + qJD(2);
t359 = t365 ^ 2;
t367 = cos(pkin(9));
t398 = pkin(3) * t367;
t366 = sin(pkin(9));
t397 = mrSges(4,2) * t366;
t363 = t367 ^ 2;
t395 = t359 * t363;
t361 = qJDD(1) + qJDD(2);
t394 = t361 * t367;
t371 = sin(qJ(1));
t375 = cos(qJ(1));
t385 = -g(2) * t375 - t371 * g(3);
t349 = qJDD(1) * pkin(1) + t385;
t389 = -g(2) * t371 + t375 * g(3);
t350 = -qJD(1) ^ 2 * pkin(1) + t389;
t370 = sin(qJ(2));
t374 = cos(qJ(2));
t337 = t370 * t349 + t374 * t350;
t334 = -pkin(2) * t359 + qJ(3) * t361 + t337;
t392 = qJD(3) * t365;
t390 = -t367 * g(1) - 0.2e1 * t366 * t392;
t315 = (-pkin(7) * t361 + t359 * t398 - t334) * t366 + t390;
t319 = -g(1) * t366 + (t334 + 0.2e1 * t392) * t367;
t316 = -pkin(3) * t395 + pkin(7) * t394 + t319;
t369 = sin(qJ(4));
t373 = cos(qJ(4));
t297 = t373 * t315 - t316 * t369;
t381 = t366 * t373 + t367 * t369;
t380 = -t366 * t369 + t367 * t373;
t342 = t380 * t365;
t391 = t342 * qJD(4);
t333 = t381 * t361 + t391;
t343 = t381 * t365;
t293 = (-t333 + t391) * pkin(8) + (t342 * t343 + qJDD(4)) * pkin(4) + t297;
t298 = t369 * t315 + t373 * t316;
t332 = -t343 * qJD(4) + t380 * t361;
t340 = qJD(4) * pkin(4) - pkin(8) * t343;
t341 = t342 ^ 2;
t294 = -pkin(4) * t341 + pkin(8) * t332 - qJD(4) * t340 + t298;
t368 = sin(qJ(5));
t372 = cos(qJ(5));
t291 = t293 * t372 - t294 * t368;
t325 = t342 * t372 - t343 * t368;
t305 = qJD(5) * t325 + t332 * t368 + t333 * t372;
t326 = t342 * t368 + t343 * t372;
t311 = -mrSges(6,1) * t325 + mrSges(6,2) * t326;
t364 = qJD(4) + qJD(5);
t320 = -mrSges(6,2) * t364 + mrSges(6,3) * t325;
t360 = qJDD(4) + qJDD(5);
t288 = m(6) * t291 + mrSges(6,1) * t360 - mrSges(6,3) * t305 - t311 * t326 + t320 * t364;
t292 = t293 * t368 + t294 * t372;
t304 = -qJD(5) * t326 + t332 * t372 - t333 * t368;
t321 = mrSges(6,1) * t364 - mrSges(6,3) * t326;
t289 = m(6) * t292 - mrSges(6,2) * t360 + mrSges(6,3) * t304 + t311 * t325 - t321 * t364;
t280 = t372 * t288 + t368 * t289;
t330 = -mrSges(5,1) * t342 + mrSges(5,2) * t343;
t338 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t342;
t278 = m(5) * t297 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t333 + qJD(4) * t338 - t330 * t343 + t280;
t339 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t343;
t386 = -t288 * t368 + t372 * t289;
t279 = m(5) * t298 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t332 - qJD(4) * t339 + t330 * t342 + t386;
t393 = t373 * t278 + t369 * t279;
t318 = -t366 * t334 + t390;
t382 = mrSges(4,3) * t361 + (-mrSges(4,1) * t367 + t397) * t359;
t387 = -t278 * t369 + t373 * t279;
t388 = -(m(4) * t318 - t382 * t366 + t393) * t366 + t367 * (m(4) * t319 + t382 * t367 + t387);
t336 = t349 * t374 - t370 * t350;
t384 = qJDD(3) - t336;
t362 = t366 ^ 2;
t317 = (-pkin(2) - t398) * t361 + (-qJ(3) + (-t362 - t363) * pkin(7)) * t359 + t384;
t296 = -pkin(4) * t332 - pkin(8) * t341 + t340 * t343 + t317;
t383 = m(6) * t296 - t304 * mrSges(6,1) + t305 * mrSges(6,2) - t325 * t320 + t326 * t321;
t306 = Ifges(6,5) * t326 + Ifges(6,6) * t325 + Ifges(6,3) * t364;
t308 = Ifges(6,1) * t326 + Ifges(6,4) * t325 + Ifges(6,5) * t364;
t281 = -mrSges(6,1) * t296 + mrSges(6,3) * t292 + Ifges(6,4) * t305 + Ifges(6,2) * t304 + Ifges(6,6) * t360 - t306 * t326 + t308 * t364;
t307 = Ifges(6,4) * t326 + Ifges(6,2) * t325 + Ifges(6,6) * t364;
t282 = mrSges(6,2) * t296 - mrSges(6,3) * t291 + Ifges(6,1) * t305 + Ifges(6,4) * t304 + Ifges(6,5) * t360 + t306 * t325 - t307 * t364;
t322 = Ifges(5,5) * t343 + Ifges(5,6) * t342 + Ifges(5,3) * qJD(4);
t324 = Ifges(5,1) * t343 + Ifges(5,4) * t342 + Ifges(5,5) * qJD(4);
t271 = -mrSges(5,1) * t317 + mrSges(5,3) * t298 + Ifges(5,4) * t333 + Ifges(5,2) * t332 + Ifges(5,6) * qJDD(4) - pkin(4) * t383 + pkin(8) * t386 + qJD(4) * t324 + t372 * t281 + t368 * t282 - t343 * t322;
t323 = Ifges(5,4) * t343 + Ifges(5,2) * t342 + Ifges(5,6) * qJD(4);
t272 = mrSges(5,2) * t317 - mrSges(5,3) * t297 + Ifges(5,1) * t333 + Ifges(5,4) * t332 + Ifges(5,5) * qJDD(4) - pkin(8) * t280 - qJD(4) * t323 - t281 * t368 + t282 * t372 + t322 * t342;
t331 = -pkin(2) * t361 - qJ(3) * t359 + t384;
t377 = m(5) * t317 - t332 * mrSges(5,1) + mrSges(5,2) * t333 - t342 * t338 + t339 * t343 + t383;
t376 = -m(4) * t331 + mrSges(4,1) * t394 - t377 + (t359 * t362 + t395) * mrSges(4,3);
t284 = t361 * t397 - t376;
t379 = -mrSges(3,2) * t337 + t367 * (-mrSges(4,1) * t331 + mrSges(4,3) * t319 + t369 * t272 + t373 * t271 - pkin(3) * t377 + pkin(7) * t387 + (Ifges(4,4) * t366 + Ifges(4,2) * t367) * t361) + t366 * (mrSges(4,2) * t331 - mrSges(4,3) * t318 + t373 * t272 - t369 * t271 - pkin(7) * t393 + (Ifges(4,1) * t366 + Ifges(4,4) * t367) * t361) + qJ(3) * t388 - pkin(2) * t284 + mrSges(3,1) * t336 + Ifges(3,3) * t361;
t378 = mrSges(6,1) * t291 - mrSges(6,2) * t292 + Ifges(6,5) * t305 + Ifges(6,6) * t304 + Ifges(6,3) * t360 + t326 * t307 - t308 * t325;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t385 - mrSges(2,2) * t389 + pkin(1) * (t370 * (m(3) * t337 - mrSges(3,1) * t359 - mrSges(3,2) * t361 + t388) + t374 * (t376 + (mrSges(3,1) - t397) * t361 + m(3) * t336 - mrSges(3,2) * t359)) + t379; t379; t284; mrSges(5,1) * t297 - mrSges(5,2) * t298 + Ifges(5,5) * t333 + Ifges(5,6) * t332 + Ifges(5,3) * qJDD(4) + pkin(4) * t280 + t323 * t343 - t324 * t342 + t378; t378;];
tauJ = t1;
