% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RPPR5
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
% tauJB [(6+4)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RPPR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR5_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR5_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR5_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR5_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR5_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:45
% EndTime: 2019-12-31 16:39:46
% DurationCPUTime: 0.51s
% Computational Cost: add. (3809->149), mult. (6306->180), div. (0->0), fcn. (2254->6), ass. (0->63)
t337 = qJD(1) ^ 2;
t334 = sin(qJ(1));
t336 = cos(qJ(1));
t320 = -t336 * g(1) - t334 * g(2);
t341 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t320;
t353 = -pkin(1) - pkin(2);
t302 = t353 * t337 + t341;
t319 = t334 * g(1) - t336 * g(2);
t340 = -t337 * qJ(2) + qJDD(2) - t319;
t305 = t353 * qJDD(1) + t340;
t331 = sin(pkin(6));
t332 = cos(pkin(6));
t299 = t332 * t302 + t331 * t305;
t296 = -t337 * pkin(3) - qJDD(1) * pkin(5) + t299;
t329 = g(3) + qJDD(3);
t333 = sin(qJ(4));
t335 = cos(qJ(4));
t293 = -t333 * t296 + t335 * t329;
t294 = t335 * t296 + t333 * t329;
t309 = Ifges(5,6) * qJD(4) + (-Ifges(5,4) * t333 - Ifges(5,2) * t335) * qJD(1);
t310 = Ifges(5,5) * qJD(4) + (-Ifges(5,1) * t333 - Ifges(5,4) * t335) * qJD(1);
t346 = qJD(1) * qJD(4);
t314 = -t333 * qJDD(1) - t335 * t346;
t315 = -t335 * qJDD(1) + t333 * t346;
t354 = mrSges(5,1) * t293 - mrSges(5,2) * t294 + Ifges(5,5) * t314 + Ifges(5,6) * t315 + Ifges(5,3) * qJDD(4) - (t309 * t333 - t310 * t335) * qJD(1);
t352 = mrSges(2,1) + mrSges(3,1);
t351 = Ifges(3,4) + Ifges(2,5);
t350 = Ifges(2,6) - Ifges(3,6);
t306 = -t337 * pkin(1) + t341;
t313 = (mrSges(5,1) * t335 - mrSges(5,2) * t333) * qJD(1);
t347 = qJD(1) * t335;
t318 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t347;
t348 = qJD(1) * t333;
t290 = m(5) * t293 + qJDD(4) * mrSges(5,1) - t314 * mrSges(5,3) + qJD(4) * t318 + t313 * t348;
t317 = qJD(4) * mrSges(5,1) + mrSges(5,3) * t348;
t291 = m(5) * t294 - qJDD(4) * mrSges(5,2) + t315 * mrSges(5,3) - qJD(4) * t317 - t313 * t347;
t285 = -t333 * t290 + t335 * t291;
t281 = m(4) * t299 - t337 * mrSges(4,1) + qJDD(1) * mrSges(4,2) + t285;
t298 = -t331 * t302 + t332 * t305;
t295 = qJDD(1) * pkin(3) - t337 * pkin(5) - t298;
t292 = -m(5) * t295 + t315 * mrSges(5,1) - t314 * mrSges(5,2) + t317 * t348 - t318 * t347;
t288 = m(4) * t298 - qJDD(1) * mrSges(4,1) - t337 * mrSges(4,2) + t292;
t344 = t332 * t281 - t331 * t288;
t342 = m(3) * t306 + qJDD(1) * mrSges(3,3) + t344;
t273 = m(2) * t320 - qJDD(1) * mrSges(2,2) - t352 * t337 + t342;
t279 = t331 * t281 + t332 * t288;
t307 = -qJDD(1) * pkin(1) + t340;
t278 = m(3) * t307 - qJDD(1) * mrSges(3,1) - t337 * mrSges(3,3) + t279;
t274 = m(2) * t319 + qJDD(1) * mrSges(2,1) - t337 * mrSges(2,2) - t278;
t349 = t334 * t273 + t336 * t274;
t345 = t336 * t273 - t334 * t274;
t284 = t335 * t290 + t333 * t291;
t283 = m(4) * t329 + t284;
t308 = Ifges(5,3) * qJD(4) + (-Ifges(5,5) * t333 - Ifges(5,6) * t335) * qJD(1);
t286 = -mrSges(5,1) * t295 + mrSges(5,3) * t294 + Ifges(5,4) * t314 + Ifges(5,2) * t315 + Ifges(5,6) * qJDD(4) + qJD(4) * t310 + t308 * t348;
t287 = mrSges(5,2) * t295 - mrSges(5,3) * t293 + Ifges(5,1) * t314 + Ifges(5,4) * t315 + Ifges(5,5) * qJDD(4) - qJD(4) * t309 - t308 * t347;
t338 = -mrSges(3,1) * t307 - mrSges(4,1) * t298 - mrSges(2,2) * t320 - pkin(2) * t279 - pkin(3) * t292 - pkin(5) * t285 - t335 * t286 - t333 * t287 + qJ(2) * (-t337 * mrSges(3,1) + t342) - pkin(1) * t278 + mrSges(4,2) * t299 + mrSges(3,3) * t306 + mrSges(2,1) * t319 + (Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1);
t282 = -m(3) * g(3) - t283;
t275 = -mrSges(4,1) * t329 + mrSges(4,3) * t299 + t337 * Ifges(4,5) - Ifges(4,6) * qJDD(1) - pkin(3) * t284 - t354;
t269 = mrSges(4,2) * t329 - mrSges(4,3) * t298 - Ifges(4,5) * qJDD(1) - t337 * Ifges(4,6) - pkin(5) * t284 - t333 * t286 + t335 * t287;
t268 = mrSges(3,2) * t307 - mrSges(2,3) * t319 - qJ(2) * t282 - qJ(3) * t279 + t332 * t269 - t331 * t275 - t350 * t337 + t351 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t267 = mrSges(3,2) * t306 + mrSges(2,3) * t320 - pkin(1) * t282 + pkin(2) * t283 + t352 * g(3) - qJ(3) * t344 + t350 * qJDD(1) - t331 * t269 - t332 * t275 + t351 * t337;
t1 = [-m(1) * g(1) + t345; -m(1) * g(2) + t349; (-m(1) - m(2) - m(3)) * g(3) - t283; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t349 - t334 * t267 + t336 * t268; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t345 + t336 * t267 + t334 * t268; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t338; t338; t278; t283; t354;];
tauJB = t1;
