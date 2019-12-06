% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-05 18:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPRP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP1_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP1_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:22:08
% EndTime: 2019-12-05 18:22:09
% DurationCPUTime: 0.59s
% Computational Cost: add. (3017->146), mult. (3889->177), div. (0->0), fcn. (1891->8), ass. (0->71)
t351 = Ifges(5,1) + Ifges(6,1);
t346 = Ifges(5,4) + Ifges(6,4);
t345 = Ifges(5,5) + Ifges(6,5);
t350 = Ifges(5,2) + Ifges(6,2);
t344 = Ifges(5,6) + Ifges(6,6);
t349 = Ifges(5,3) + Ifges(6,3);
t313 = qJD(1) + qJD(2);
t311 = t313 ^ 2;
t348 = pkin(4) * t311;
t347 = -mrSges(5,2) - mrSges(6,2);
t319 = sin(qJ(4));
t343 = t313 * t319;
t322 = cos(qJ(4));
t342 = t313 * t322;
t321 = sin(qJ(1));
t324 = cos(qJ(1));
t336 = t324 * g(2) + t321 * g(3);
t297 = qJDD(1) * pkin(1) + t336;
t330 = t321 * g(2) - t324 * g(3);
t298 = -qJD(1) ^ 2 * pkin(1) + t330;
t320 = sin(qJ(2));
t323 = cos(qJ(2));
t275 = t323 * t297 - t320 * t298;
t312 = qJDD(1) + qJDD(2);
t272 = t312 * pkin(2) + t275;
t276 = t320 * t297 + t323 * t298;
t273 = -t311 * pkin(2) + t276;
t317 = sin(pkin(8));
t318 = cos(pkin(8));
t268 = t317 * t272 + t318 * t273;
t265 = -t311 * pkin(3) + t312 * pkin(7) + t268;
t316 = -g(1) + qJDD(3);
t305 = t322 * t316;
t261 = -t319 * t265 + t305;
t289 = (-mrSges(6,1) * t322 + mrSges(6,2) * t319) * t313;
t290 = (-mrSges(5,1) * t322 + mrSges(5,2) * t319) * t313;
t335 = qJD(4) * t313;
t331 = t322 * t335;
t291 = t319 * t312 + t331;
t303 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t342;
t334 = qJD(5) * t313;
t258 = qJDD(4) * pkin(4) + t305 + (-t291 + t331) * qJ(5) + (t322 * t348 - t265 - 0.2e1 * t334) * t319;
t302 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t342;
t333 = m(6) * t258 + qJDD(4) * mrSges(6,1) + qJD(4) * t302;
t252 = m(5) * t261 + qJDD(4) * mrSges(5,1) + qJD(4) * t303 + (-t289 - t290) * t343 + (-mrSges(5,3) - mrSges(6,3)) * t291 + t333;
t262 = t322 * t265 + t319 * t316;
t292 = t322 * t312 - t319 * t335;
t299 = qJD(4) * pkin(4) - qJ(5) * t343;
t315 = t322 ^ 2;
t259 = t292 * qJ(5) - qJD(4) * t299 - t315 * t348 + 0.2e1 * t322 * t334 + t262;
t332 = m(6) * t259 + t292 * mrSges(6,3) + t289 * t342;
t300 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t343;
t337 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t343 - t300;
t253 = m(5) * t262 + t292 * mrSges(5,3) + t337 * qJD(4) + t347 * qJDD(4) + t290 * t342 + t332;
t329 = -t319 * t252 + t322 * t253;
t246 = m(4) * t268 - t311 * mrSges(4,1) - t312 * mrSges(4,2) + t329;
t267 = t318 * t272 - t317 * t273;
t327 = -t312 * pkin(3) - t267;
t264 = -t311 * pkin(7) + t327;
t260 = t299 * t343 - t292 * pkin(4) + qJDD(5) + (-qJ(5) * t315 - pkin(7)) * t311 + t327;
t328 = -m(6) * t260 + t292 * mrSges(6,1) + t302 * t342;
t325 = -m(5) * t264 + t292 * mrSges(5,1) + t347 * t291 + t303 * t342 + t337 * t343 + t328;
t249 = m(4) * t267 + t312 * mrSges(4,1) - t311 * mrSges(4,2) + t325;
t341 = t317 * t246 + t318 * t249;
t340 = (t319 * t345 + t344 * t322) * t313 + t349 * qJD(4);
t339 = (t319 * t346 + t350 * t322) * t313 + t344 * qJD(4);
t338 = (-t351 * t319 - t322 * t346) * t313 - t345 * qJD(4);
t254 = t291 * mrSges(6,2) + t300 * t343 - t328;
t255 = -t291 * mrSges(6,3) - t289 * t343 + t333;
t326 = -mrSges(3,2) * t276 - mrSges(4,2) * t268 + pkin(2) * t341 + t322 * (-mrSges(5,1) * t264 + mrSges(5,3) * t262 - mrSges(6,1) * t260 + mrSges(6,3) * t259 - pkin(4) * t254 + qJ(5) * t332 - t340 * t343 + t350 * t292 + t346 * t291 + (-qJ(5) * mrSges(6,2) + t344) * qJDD(4) + (-qJ(5) * t300 - t338) * qJD(4)) + t319 * (mrSges(5,2) * t264 + mrSges(6,2) * t260 - mrSges(5,3) * t261 - mrSges(6,3) * t258 - qJ(5) * t255 - t339 * qJD(4) + t345 * qJDD(4) + t351 * t291 + t346 * t292 + t340 * t342) + pkin(7) * t329 + pkin(3) * t325 + mrSges(4,1) * t267 + mrSges(3,1) * t275 + (Ifges(4,3) + Ifges(3,3)) * t312;
t1 = [t326 + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t336 + pkin(1) * (t320 * (m(3) * t276 - t311 * mrSges(3,1) - t312 * mrSges(3,2) + t318 * t246 - t317 * t249) + t323 * (m(3) * t275 + t312 * mrSges(3,1) - t311 * mrSges(3,2) + t341)) - mrSges(2,2) * t330; t326; m(4) * t316 + t322 * t252 + t319 * t253; mrSges(5,1) * t261 + mrSges(6,1) * t258 - mrSges(5,2) * t262 - mrSges(6,2) * t259 + pkin(4) * t255 + t344 * t292 + t345 * t291 + t349 * qJDD(4) + (t339 * t319 + t338 * t322) * t313; t254;];
tauJ = t1;
