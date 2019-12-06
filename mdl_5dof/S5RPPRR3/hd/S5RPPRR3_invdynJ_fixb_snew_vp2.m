% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPPRR3
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-05 17:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPRR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR3_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR3_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:41:54
% EndTime: 2019-12-05 17:41:55
% DurationCPUTime: 1.15s
% Computational Cost: add. (7485->192), mult. (16642->248), div. (0->0), fcn. (11244->10), ass. (0->90)
t377 = qJD(1) ^ 2;
t369 = cos(pkin(9));
t399 = t369 * pkin(3);
t367 = sin(pkin(9));
t398 = mrSges(4,2) * t367;
t364 = t369 ^ 2;
t397 = t364 * t377;
t373 = sin(qJ(1));
t376 = cos(qJ(1));
t394 = t376 * g(2) + t373 * g(3);
t349 = qJDD(1) * pkin(1) + t394;
t389 = t373 * g(2) - g(3) * t376;
t350 = -pkin(1) * t377 + t389;
t368 = sin(pkin(8));
t370 = cos(pkin(8));
t337 = t368 * t349 + t370 * t350;
t330 = -pkin(2) * t377 + qJDD(1) * qJ(3) + t337;
t366 = -g(1) + qJDD(2);
t391 = qJD(1) * qJD(3);
t395 = t369 * t366 - 0.2e1 * t367 * t391;
t316 = (-pkin(6) * qJDD(1) + t377 * t399 - t330) * t367 + t395;
t320 = t367 * t366 + (t330 + 0.2e1 * t391) * t369;
t390 = t369 * qJDD(1);
t317 = -pkin(3) * t397 + pkin(6) * t390 + t320;
t372 = sin(qJ(4));
t375 = cos(qJ(4));
t298 = t375 * t316 - t317 * t372;
t383 = t367 * t375 + t369 * t372;
t382 = -t367 * t372 + t369 * t375;
t342 = t382 * qJD(1);
t392 = qJD(4) * t342;
t335 = t383 * qJDD(1) + t392;
t343 = t383 * qJD(1);
t294 = (-t335 + t392) * pkin(7) + (t342 * t343 + qJDD(4)) * pkin(4) + t298;
t299 = t372 * t316 + t375 * t317;
t334 = -qJD(4) * t343 + t382 * qJDD(1);
t340 = qJD(4) * pkin(4) - pkin(7) * t343;
t341 = t342 ^ 2;
t295 = -pkin(4) * t341 + pkin(7) * t334 - qJD(4) * t340 + t299;
t371 = sin(qJ(5));
t374 = cos(qJ(5));
t292 = t294 * t374 - t295 * t371;
t328 = t342 * t374 - t343 * t371;
t306 = qJD(5) * t328 + t334 * t371 + t335 * t374;
t329 = t342 * t371 + t343 * t374;
t312 = -mrSges(6,1) * t328 + mrSges(6,2) * t329;
t365 = qJD(4) + qJD(5);
t321 = -mrSges(6,2) * t365 + mrSges(6,3) * t328;
t362 = qJDD(4) + qJDD(5);
t289 = m(6) * t292 + mrSges(6,1) * t362 - mrSges(6,3) * t306 - t312 * t329 + t321 * t365;
t293 = t294 * t371 + t295 * t374;
t305 = -qJD(5) * t329 + t334 * t374 - t335 * t371;
t322 = mrSges(6,1) * t365 - mrSges(6,3) * t329;
t290 = m(6) * t293 - mrSges(6,2) * t362 + mrSges(6,3) * t305 + t312 * t328 - t322 * t365;
t282 = t374 * t289 + t371 * t290;
t332 = -mrSges(5,1) * t342 + mrSges(5,2) * t343;
t338 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t342;
t280 = m(5) * t298 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t335 + qJD(4) * t338 - t332 * t343 + t282;
t339 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t343;
t386 = -t289 * t371 + t374 * t290;
t281 = m(5) * t299 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t334 - qJD(4) * t339 + t332 * t342 + t386;
t396 = t375 * t280 + t372 * t281;
t319 = -t330 * t367 + t395;
t381 = mrSges(4,3) * qJDD(1) + t377 * (-t369 * mrSges(4,1) + t398);
t275 = m(4) * t319 - t381 * t367 + t396;
t387 = -t280 * t372 + t375 * t281;
t276 = m(4) * t320 + t381 * t369 + t387;
t388 = -t275 * t367 + t369 * t276;
t336 = t349 * t370 - t368 * t350;
t385 = qJDD(3) - t336;
t363 = t367 ^ 2;
t318 = (-pkin(2) - t399) * qJDD(1) + (-qJ(3) + (-t363 - t364) * pkin(6)) * t377 + t385;
t297 = -pkin(4) * t334 - pkin(7) * t341 + t340 * t343 + t318;
t384 = m(6) * t297 - t305 * mrSges(6,1) + t306 * mrSges(6,2) - t328 * t321 + t329 * t322;
t308 = Ifges(6,4) * t329 + Ifges(6,2) * t328 + Ifges(6,6) * t365;
t309 = Ifges(6,1) * t329 + Ifges(6,4) * t328 + Ifges(6,5) * t365;
t380 = mrSges(6,1) * t292 - mrSges(6,2) * t293 + Ifges(6,5) * t306 + Ifges(6,6) * t305 + Ifges(6,3) * t362 + t329 * t308 - t309 * t328;
t379 = m(5) * t318 - t334 * mrSges(5,1) + mrSges(5,2) * t335 - t342 * t338 + t339 * t343 + t384;
t324 = -qJDD(1) * pkin(2) - qJ(3) * t377 + t385;
t378 = -m(4) * t324 + mrSges(4,1) * t390 - t379 + (t363 * t377 + t397) * mrSges(4,3);
t327 = Ifges(5,1) * t343 + Ifges(5,4) * t342 + Ifges(5,5) * qJD(4);
t326 = Ifges(5,4) * t343 + Ifges(5,2) * t342 + Ifges(5,6) * qJD(4);
t325 = Ifges(5,5) * t343 + Ifges(5,6) * t342 + Ifges(5,3) * qJD(4);
t307 = Ifges(6,5) * t329 + Ifges(6,6) * t328 + Ifges(6,3) * t365;
t285 = qJDD(1) * t398 - t378;
t284 = mrSges(6,2) * t297 - mrSges(6,3) * t292 + Ifges(6,1) * t306 + Ifges(6,4) * t305 + Ifges(6,5) * t362 + t307 * t328 - t308 * t365;
t283 = -mrSges(6,1) * t297 + mrSges(6,3) * t293 + Ifges(6,4) * t306 + Ifges(6,2) * t305 + Ifges(6,6) * t362 - t307 * t329 + t309 * t365;
t273 = mrSges(5,2) * t318 - mrSges(5,3) * t298 + Ifges(5,1) * t335 + Ifges(5,4) * t334 + Ifges(5,5) * qJDD(4) - pkin(7) * t282 - qJD(4) * t326 - t283 * t371 + t284 * t374 + t325 * t342;
t272 = -mrSges(5,1) * t318 + mrSges(5,3) * t299 + Ifges(5,4) * t335 + Ifges(5,2) * t334 + Ifges(5,6) * qJDD(4) - pkin(4) * t384 + pkin(7) * t386 + qJD(4) * t327 + t374 * t283 + t371 * t284 - t343 * t325;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t394 - mrSges(2,2) * t389 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t336 - mrSges(3,2) * t337 + t367 * (mrSges(4,2) * t324 - mrSges(4,3) * t319 + t375 * t273 - t372 * t272 - pkin(6) * t396 + (Ifges(4,1) * t367 + Ifges(4,4) * t369) * qJDD(1)) + t369 * (-mrSges(4,1) * t324 + mrSges(4,3) * t320 + t372 * t273 + t375 * t272 - pkin(3) * t379 + pkin(6) * t387 + (Ifges(4,4) * t367 + Ifges(4,2) * t369) * qJDD(1)) - pkin(2) * t285 + qJ(3) * t388 + pkin(1) * (t368 * (m(3) * t337 - mrSges(3,1) * t377 - qJDD(1) * mrSges(3,2) + t388) + t370 * (t378 + (mrSges(3,1) - t398) * qJDD(1) + m(3) * t336 - mrSges(3,2) * t377)); m(3) * t366 + t275 * t369 + t276 * t367; t285; mrSges(5,1) * t298 - mrSges(5,2) * t299 + Ifges(5,5) * t335 + Ifges(5,6) * t334 + Ifges(5,3) * qJDD(4) + pkin(4) * t282 + t326 * t343 - t327 * t342 + t380; t380;];
tauJ = t1;
