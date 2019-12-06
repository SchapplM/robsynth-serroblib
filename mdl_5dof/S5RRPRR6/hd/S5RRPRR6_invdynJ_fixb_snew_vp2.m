% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPRR6
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
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPRR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR6_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR6_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:35:55
% EndTime: 2019-12-05 18:35:56
% DurationCPUTime: 1.33s
% Computational Cost: add. (12557->186), mult. (17171->250), div. (0->0), fcn. (10524->10), ass. (0->89)
t369 = sin(qJ(1));
t373 = cos(qJ(1));
t391 = t373 * g(2) + t369 * g(3);
t343 = qJDD(1) * pkin(1) + t391;
t384 = t369 * g(2) - g(3) * t373;
t344 = -qJD(1) ^ 2 * pkin(1) + t384;
t368 = sin(qJ(2));
t372 = cos(qJ(2));
t326 = t368 * t343 + t372 * t344;
t362 = qJD(1) + qJD(2);
t359 = t362 ^ 2;
t360 = qJDD(1) + qJDD(2);
t399 = -pkin(2) * t359 + qJ(3) * t360 + 0.2e1 * qJD(3) * t362 + t326;
t364 = sin(pkin(9));
t365 = cos(pkin(9));
t312 = -t365 * g(1) - t399 * t364;
t398 = mrSges(4,2) * t364;
t397 = mrSges(4,3) * t360;
t396 = t359 * t364 ^ 2;
t395 = t360 * t365;
t394 = t362 * t364;
t393 = t362 * t365;
t313 = -g(1) * t364 + t399 * t365;
t381 = -pkin(3) * t365 - pkin(7) * t364;
t338 = t381 * t362;
t301 = t338 * t393 + t313;
t325 = t343 * t372 - t368 * t344;
t379 = -qJ(3) * t359 + qJDD(3) - t325;
t311 = (-pkin(2) + t381) * t360 + t379;
t371 = cos(qJ(4));
t310 = t371 * t311;
t367 = sin(qJ(4));
t389 = qJD(4) * t362;
t335 = (t360 * t371 - t367 * t389) * t364;
t347 = qJDD(4) - t395;
t348 = qJD(4) - t393;
t293 = pkin(4) * t347 - pkin(8) * t335 + t310 + (-pkin(4) * t371 * t396 - pkin(8) * t348 * t394 - t301) * t367;
t296 = t371 * t301 + t367 * t311;
t386 = t371 * t394;
t333 = pkin(4) * t348 - pkin(8) * t386;
t334 = (-t360 * t367 - t371 * t389) * t364;
t388 = t367 ^ 2 * t396;
t294 = -pkin(4) * t388 + pkin(8) * t334 - t333 * t348 + t296;
t366 = sin(qJ(5));
t370 = cos(qJ(5));
t291 = t293 * t370 - t294 * t366;
t327 = (-t366 * t371 - t367 * t370) * t394;
t305 = qJD(5) * t327 + t334 * t366 + t335 * t370;
t328 = (-t367 * t366 + t371 * t370) * t394;
t314 = -mrSges(6,1) * t327 + mrSges(6,2) * t328;
t346 = qJD(5) + t348;
t319 = -mrSges(6,2) * t346 + mrSges(6,3) * t327;
t345 = qJDD(5) + t347;
t288 = m(6) * t291 + mrSges(6,1) * t345 - mrSges(6,3) * t305 - t314 * t328 + t319 * t346;
t292 = t293 * t366 + t294 * t370;
t304 = -qJD(5) * t328 + t334 * t370 - t335 * t366;
t320 = mrSges(6,1) * t346 - mrSges(6,3) * t328;
t289 = m(6) * t292 - mrSges(6,2) * t345 + mrSges(6,3) * t304 + t314 * t327 - t320 * t346;
t281 = t370 * t288 + t366 * t289;
t387 = t367 * t394;
t295 = -t301 * t367 + t310;
t330 = -mrSges(5,2) * t348 - mrSges(5,3) * t387;
t332 = (t367 * mrSges(5,1) + t371 * mrSges(5,2)) * t394;
t279 = m(5) * t295 + mrSges(5,1) * t347 - mrSges(5,3) * t335 + t330 * t348 - t332 * t386 + t281;
t331 = mrSges(5,1) * t348 - mrSges(5,3) * t386;
t382 = -t288 * t366 + t370 * t289;
t280 = m(5) * t296 - mrSges(5,2) * t347 + mrSges(5,3) * t334 - t331 * t348 - t332 * t387 + t382;
t300 = t338 * t394 - t312;
t336 = (-mrSges(4,1) * t365 + t398) * t362;
t297 = -pkin(4) * t334 - pkin(8) * t388 + t333 * t386 + t300;
t375 = m(6) * t297 - t304 * mrSges(6,1) + t305 * mrSges(6,2) - t319 * t327 + t328 * t320;
t383 = -(m(4) * t312 - m(5) * t300 + mrSges(5,1) * t334 - mrSges(5,2) * t335 + (-t397 + (-t330 * t367 - t331 * t371 - t336) * t362) * t364 - t375) * t364 + t365 * (m(4) * t313 - t279 * t367 + t280 * t371 + (t336 * t362 + t397) * t365);
t278 = t279 * t371 + t280 * t367;
t322 = Ifges(5,6) * t348 + (t371 * Ifges(5,4) - t367 * Ifges(5,2)) * t394;
t323 = Ifges(5,5) * t348 + (t371 * Ifges(5,1) - t367 * Ifges(5,4)) * t394;
t380 = t371 * t322 + t367 * t323;
t317 = -pkin(2) * t360 + t379;
t376 = -m(4) * t317 + mrSges(4,1) * t395 - t278 + (t359 * t365 ^ 2 + t396) * mrSges(4,3);
t277 = t360 * t398 - t376;
t306 = Ifges(6,5) * t328 + Ifges(6,6) * t327 + Ifges(6,3) * t346;
t308 = Ifges(6,1) * t328 + Ifges(6,4) * t327 + Ifges(6,5) * t346;
t282 = -mrSges(6,1) * t297 + mrSges(6,3) * t292 + Ifges(6,4) * t305 + Ifges(6,2) * t304 + Ifges(6,6) * t345 - t306 * t328 + t308 * t346;
t307 = Ifges(6,4) * t328 + Ifges(6,2) * t327 + Ifges(6,6) * t346;
t283 = mrSges(6,2) * t297 - mrSges(6,3) * t291 + Ifges(6,1) * t305 + Ifges(6,4) * t304 + Ifges(6,5) * t345 + t306 * t327 - t307 * t346;
t337 = (Ifges(4,5) * t364 + Ifges(4,6) * t365) * t362;
t377 = -mrSges(6,1) * t291 + mrSges(6,2) * t292 - Ifges(6,5) * t305 - Ifges(6,6) * t304 - Ifges(6,3) * t345 - t328 * t307 + t327 * t308;
t374 = mrSges(5,1) * t295 - mrSges(5,2) * t296 + Ifges(5,5) * t335 + Ifges(5,6) * t334 + Ifges(5,3) * t347 + pkin(4) * t281 - t377;
t378 = -mrSges(3,2) * t326 + t364 * (t337 * t393 + mrSges(4,2) * t317 - mrSges(4,3) * t312 + t371 * (mrSges(5,2) * t300 - mrSges(5,3) * t295 + Ifges(5,1) * t335 + Ifges(5,4) * t334 + Ifges(5,5) * t347 - pkin(8) * t281 - t282 * t366 + t283 * t370 - t322 * t348) - t367 * (-mrSges(5,1) * t300 + mrSges(5,3) * t296 + Ifges(5,4) * t335 + Ifges(5,2) * t334 + Ifges(5,6) * t347 - pkin(4) * t375 + pkin(8) * t382 + t370 * t282 + t366 * t283 + t348 * t323) - pkin(7) * t278 + (Ifges(4,1) * t364 + Ifges(4,4) * t365) * t360) + t365 * ((Ifges(4,4) * t360 + (-t337 - t380) * t362) * t364 - mrSges(4,1) * t317 + mrSges(4,3) * t313 - pkin(3) * t278 + Ifges(4,2) * t395 - t374) + qJ(3) * t383 - pkin(2) * t277 + mrSges(3,1) * t325 + Ifges(3,3) * t360;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t391 - mrSges(2,2) * t384 + pkin(1) * (t368 * (m(3) * t326 - mrSges(3,1) * t359 - mrSges(3,2) * t360 + t383) + t372 * (m(3) * t325 - mrSges(3,2) * t359 + (mrSges(3,1) - t398) * t360 + t376)) + t378; t378; t277; t380 * t394 + t374; -t377;];
tauJ = t1;
