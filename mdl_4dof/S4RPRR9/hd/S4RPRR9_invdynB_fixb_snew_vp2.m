% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RPRR9
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
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RPRR9_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR9_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR9_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR9_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR9_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR9_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR9_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:12
% EndTime: 2019-12-31 16:56:13
% DurationCPUTime: 0.71s
% Computational Cost: add. (4707->191), mult. (8774->233), div. (0->0), fcn. (4459->6), ass. (0->78)
t335 = sin(qJ(1));
t338 = cos(qJ(1));
t325 = -t338 * g(1) - t335 * g(2);
t362 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t325;
t361 = -pkin(1) - pkin(5);
t334 = sin(qJ(3));
t360 = t334 * g(3);
t359 = mrSges(2,1) - mrSges(3,2);
t358 = -Ifges(3,4) + Ifges(2,5);
t357 = (Ifges(3,5) - Ifges(2,6));
t324 = t335 * g(1) - t338 * g(2);
t340 = qJD(1) ^ 2;
t345 = -t340 * qJ(2) + qJDD(2) - t324;
t306 = t361 * qJDD(1) + t345;
t337 = cos(qJ(3));
t300 = -t337 * g(3) + t334 * t306;
t318 = (mrSges(4,1) * t334 + mrSges(4,2) * t337) * qJD(1);
t353 = qJD(1) * qJD(3);
t350 = t337 * t353;
t320 = -t334 * qJDD(1) - t350;
t355 = qJD(1) * t337;
t323 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t355;
t305 = t361 * t340 - t362;
t351 = t334 * t353;
t321 = t337 * qJDD(1) - t351;
t288 = (-t321 + t351) * pkin(6) + (-t320 + t350) * pkin(3) + t305;
t319 = (pkin(3) * t334 - pkin(6) * t337) * qJD(1);
t339 = qJD(3) ^ 2;
t354 = t334 * qJD(1);
t290 = -t339 * pkin(3) + qJDD(3) * pkin(6) - t319 * t354 + t300;
t333 = sin(qJ(4));
t336 = cos(qJ(4));
t286 = t336 * t288 - t333 * t290;
t316 = t336 * qJD(3) - t333 * t355;
t297 = t316 * qJD(4) + t333 * qJDD(3) + t336 * t321;
t317 = t333 * qJD(3) + t336 * t355;
t298 = -t316 * mrSges(5,1) + t317 * mrSges(5,2);
t326 = qJD(4) + t354;
t301 = -t326 * mrSges(5,2) + t316 * mrSges(5,3);
t315 = qJDD(4) - t320;
t284 = m(5) * t286 + t315 * mrSges(5,1) - t297 * mrSges(5,3) - t317 * t298 + t326 * t301;
t287 = t333 * t288 + t336 * t290;
t296 = -t317 * qJD(4) + t336 * qJDD(3) - t333 * t321;
t302 = t326 * mrSges(5,1) - t317 * mrSges(5,3);
t285 = m(5) * t287 - t315 * mrSges(5,2) + t296 * mrSges(5,3) + t316 * t298 - t326 * t302;
t347 = -t333 * t284 + t336 * t285;
t276 = m(4) * t300 - qJDD(3) * mrSges(4,2) + t320 * mrSges(4,3) - qJD(3) * t323 - t318 * t354 + t347;
t299 = t337 * t306 + t360;
t322 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t354;
t289 = -qJDD(3) * pkin(3) - t339 * pkin(6) - t360 + (qJD(1) * t319 - t306) * t337;
t342 = -m(5) * t289 + t296 * mrSges(5,1) - t297 * mrSges(5,2) + t316 * t301 - t317 * t302;
t280 = m(4) * t299 + qJDD(3) * mrSges(4,1) - t321 * mrSges(4,3) + qJD(3) * t322 - t318 * t355 + t342;
t271 = t334 * t276 + t337 * t280;
t308 = -qJDD(1) * pkin(1) + t345;
t344 = -m(3) * t308 + (t340 * mrSges(3,3)) - t271;
t269 = m(2) * t324 - (t340 * mrSges(2,2)) + t359 * qJDD(1) + t344;
t307 = t340 * pkin(1) + t362;
t277 = t336 * t284 + t333 * t285;
t343 = -m(4) * t305 + t320 * mrSges(4,1) - t321 * mrSges(4,2) - t322 * t354 - t323 * t355 - t277;
t341 = -m(3) * t307 + (t340 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t343;
t274 = m(2) * t325 - (t340 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t341;
t356 = t338 * t269 + t335 * t274;
t349 = -t335 * t269 + t338 * t274;
t348 = t337 * t276 - t334 * t280;
t311 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t337 - Ifges(4,4) * t334) * qJD(1);
t310 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t337 - Ifges(4,2) * t334) * qJD(1);
t309 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t337 - Ifges(4,6) * t334) * qJD(1);
t293 = Ifges(5,1) * t317 + Ifges(5,4) * t316 + Ifges(5,5) * t326;
t292 = Ifges(5,4) * t317 + Ifges(5,2) * t316 + Ifges(5,6) * t326;
t291 = Ifges(5,5) * t317 + Ifges(5,6) * t316 + Ifges(5,3) * t326;
t279 = mrSges(5,2) * t289 - mrSges(5,3) * t286 + Ifges(5,1) * t297 + Ifges(5,4) * t296 + Ifges(5,5) * t315 + t316 * t291 - t326 * t292;
t278 = -mrSges(5,1) * t289 + mrSges(5,3) * t287 + Ifges(5,4) * t297 + Ifges(5,2) * t296 + Ifges(5,6) * t315 - t317 * t291 + t326 * t293;
t270 = -m(3) * g(3) + t348;
t267 = -mrSges(4,1) * t305 - mrSges(5,1) * t286 + mrSges(5,2) * t287 + mrSges(4,3) * t300 + Ifges(4,4) * t321 - Ifges(5,5) * t297 + Ifges(4,2) * t320 + Ifges(4,6) * qJDD(3) - Ifges(5,6) * t296 - Ifges(5,3) * t315 - pkin(3) * t277 + qJD(3) * t311 - t317 * t292 + t316 * t293 - t309 * t355;
t266 = mrSges(4,2) * t305 - mrSges(4,3) * t299 + Ifges(4,1) * t321 + Ifges(4,4) * t320 + Ifges(4,5) * qJDD(3) - pkin(6) * t277 - qJD(3) * t310 - t333 * t278 + t336 * t279 - t309 * t354;
t265 = -qJ(2) * t270 - mrSges(2,3) * t324 + pkin(2) * t271 + mrSges(3,1) * t308 + t336 * t278 + pkin(3) * t342 + pkin(6) * t347 + Ifges(4,5) * t321 + Ifges(4,6) * t320 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t299 - mrSges(4,2) * t300 + t333 * t279 + (t357 * t340) + t358 * qJDD(1) + (t337 * t310 + t334 * t311) * qJD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t264 = -mrSges(3,1) * t307 + mrSges(2,3) * t325 - pkin(1) * t270 - pkin(2) * t343 - pkin(5) * t348 + t359 * g(3) - t357 * qJDD(1) - t334 * t266 - t337 * t267 + t358 * t340;
t1 = [-m(1) * g(1) + t349; -m(1) * g(2) + t356; (-m(1) - m(2) - m(3)) * g(3) + t348; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t356 - t335 * t264 + t338 * t265; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t349 + t338 * t264 + t335 * t265; mrSges(2,1) * t324 - mrSges(2,2) * t325 + pkin(1) * t344 + qJ(2) * t341 + mrSges(3,2) * t308 - mrSges(3,3) * t307 + t337 * t266 - t334 * t267 - pkin(5) * t271 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB = t1;
