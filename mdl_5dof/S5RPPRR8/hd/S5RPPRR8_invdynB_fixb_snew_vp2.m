% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPRR8
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPRR8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR8_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR8_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR8_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR8_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR8_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:09
% EndTime: 2019-12-31 18:01:10
% DurationCPUTime: 1.02s
% Computational Cost: add. (13720->182), mult. (18963->217), div. (0->0), fcn. (6806->8), ass. (0->76)
t360 = -pkin(1) - pkin(2);
t359 = mrSges(2,1) + mrSges(3,1);
t358 = Ifges(3,4) + Ifges(2,5);
t357 = Ifges(2,6) - Ifges(3,6);
t328 = -qJD(1) + qJD(4);
t336 = sin(qJ(5));
t356 = t328 * t336;
t339 = cos(qJ(5));
t355 = t328 * t339;
t338 = sin(qJ(1));
t341 = cos(qJ(1));
t323 = -t341 * g(1) - t338 * g(2);
t342 = qJD(1) ^ 2;
t346 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t323;
t313 = -t342 * pkin(1) + t346;
t307 = t360 * t342 + t346;
t322 = t338 * g(1) - t341 * g(2);
t345 = -t342 * qJ(2) + qJDD(2) - t322;
t309 = t360 * qJDD(1) + t345;
t334 = sin(pkin(8));
t335 = cos(pkin(8));
t302 = -t334 * t307 + t335 * t309;
t300 = -qJDD(1) * pkin(3) + t302;
t303 = t335 * t307 + t334 * t309;
t301 = -t342 * pkin(3) + t303;
t337 = sin(qJ(4));
t340 = cos(qJ(4));
t297 = t337 * t300 + t340 * t301;
t326 = t328 ^ 2;
t327 = -qJDD(1) + qJDD(4);
t295 = -t326 * pkin(4) + t327 * pkin(7) + t297;
t332 = g(3) + qJDD(3);
t292 = -t336 * t295 + t339 * t332;
t317 = (-mrSges(6,1) * t339 + mrSges(6,2) * t336) * t328;
t353 = qJD(5) * t328;
t318 = t336 * t327 + t339 * t353;
t321 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t355;
t290 = m(6) * t292 + qJDD(5) * mrSges(6,1) - t318 * mrSges(6,3) + qJD(5) * t321 - t317 * t356;
t293 = t339 * t295 + t336 * t332;
t319 = t339 * t327 - t336 * t353;
t320 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t356;
t291 = m(6) * t293 - qJDD(5) * mrSges(6,2) + t319 * mrSges(6,3) - qJD(5) * t320 + t317 * t355;
t349 = -t336 * t290 + t339 * t291;
t280 = m(5) * t297 - t326 * mrSges(5,1) - t327 * mrSges(5,2) + t349;
t296 = t340 * t300 - t337 * t301;
t294 = -t327 * pkin(4) - t326 * pkin(7) - t296;
t343 = -m(6) * t294 + t319 * mrSges(6,1) - t318 * mrSges(6,2) - t320 * t356 + t321 * t355;
t286 = m(5) * t296 + t327 * mrSges(5,1) - t326 * mrSges(5,2) + t343;
t277 = t337 * t280 + t340 * t286;
t274 = m(4) * t302 - qJDD(1) * mrSges(4,1) - t342 * mrSges(4,2) + t277;
t350 = t340 * t280 - t337 * t286;
t275 = m(4) * t303 - t342 * mrSges(4,1) + qJDD(1) * mrSges(4,2) + t350;
t351 = -t334 * t274 + t335 * t275;
t347 = m(3) * t313 + qJDD(1) * mrSges(3,3) + t351;
t269 = m(2) * t323 - qJDD(1) * mrSges(2,2) - t359 * t342 + t347;
t271 = t335 * t274 + t334 * t275;
t316 = -qJDD(1) * pkin(1) + t345;
t344 = -m(3) * t316 + qJDD(1) * mrSges(3,1) + t342 * mrSges(3,3) - t271;
t270 = m(2) * t322 + qJDD(1) * mrSges(2,1) - t342 * mrSges(2,2) + t344;
t354 = t338 * t269 + t341 * t270;
t282 = t339 * t290 + t336 * t291;
t352 = t341 * t269 - t338 * t270;
t348 = (-m(4) - m(5)) * t332 - t282;
t312 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t336 + Ifges(6,4) * t339) * t328;
t311 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t336 + Ifges(6,2) * t339) * t328;
t310 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t336 + Ifges(6,6) * t339) * t328;
t284 = mrSges(6,2) * t294 - mrSges(6,3) * t292 + Ifges(6,1) * t318 + Ifges(6,4) * t319 + Ifges(6,5) * qJDD(5) - qJD(5) * t311 + t310 * t355;
t283 = -mrSges(6,1) * t294 + mrSges(6,3) * t293 + Ifges(6,4) * t318 + Ifges(6,2) * t319 + Ifges(6,6) * qJDD(5) + qJD(5) * t312 - t310 * t356;
t281 = -m(3) * g(3) + t348;
t276 = -mrSges(5,1) * t332 - mrSges(6,1) * t292 + mrSges(6,2) * t293 + mrSges(5,3) * t297 + t326 * Ifges(5,5) - Ifges(6,5) * t318 + Ifges(5,6) * t327 - Ifges(6,6) * t319 - Ifges(6,3) * qJDD(5) - pkin(4) * t282 + (-t311 * t336 + t312 * t339) * t328;
t272 = mrSges(5,2) * t332 - mrSges(5,3) * t296 + Ifges(5,5) * t327 - t326 * Ifges(5,6) - pkin(7) * t282 - t336 * t283 + t339 * t284;
t265 = mrSges(4,2) * t332 - mrSges(4,3) * t302 - Ifges(4,5) * qJDD(1) - t342 * Ifges(4,6) - pkin(6) * t277 + t340 * t272 - t337 * t276;
t264 = -Ifges(4,6) * qJDD(1) + t342 * Ifges(4,5) - mrSges(4,1) * t332 + mrSges(4,3) * t303 + t337 * t272 + t340 * t276 - pkin(3) * (m(5) * t332 + t282) + pkin(6) * t350;
t263 = mrSges(3,2) * t316 - mrSges(2,3) * t322 - qJ(2) * t281 - qJ(3) * t271 - t334 * t264 + t335 * t265 - t357 * t342 + t358 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t262 = mrSges(3,2) * t313 + mrSges(2,3) * t323 - pkin(1) * t281 - pkin(2) * t348 + t359 * g(3) - qJ(3) * t351 + t357 * qJDD(1) - t335 * t264 - t334 * t265 + t358 * t342;
t1 = [-m(1) * g(1) + t352; -m(1) * g(2) + t354; (-m(1) - m(2) - m(3)) * g(3) + t348; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t354 - t338 * t262 + t341 * t263; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t352 + t341 * t262 + t338 * t263; qJ(2) * (-t342 * mrSges(3,1) + t347) + mrSges(2,1) * t322 - mrSges(2,2) * t323 + pkin(1) * t344 - pkin(2) * t271 + mrSges(3,3) * t313 - mrSges(3,1) * t316 - pkin(3) * t277 - mrSges(4,1) * t302 + mrSges(4,2) * t303 - t336 * t284 - t339 * t283 - pkin(4) * t343 - pkin(7) * t349 - mrSges(5,1) * t296 + mrSges(5,2) * t297 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - Ifges(5,3) * t327 + (Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1);];
tauB = t1;
