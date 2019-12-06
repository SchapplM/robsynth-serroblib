% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPPRR1
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPRR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR1_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR1_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR1_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR1_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR1_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:38:05
% EndTime: 2019-12-05 17:38:07
% DurationCPUTime: 0.60s
% Computational Cost: add. (1733->158), mult. (3301->199), div. (0->0), fcn. (1606->6), ass. (0->62)
t350 = 2 * qJD(1);
t329 = qJD(1) ^ 2;
t325 = sin(qJ(1));
t328 = cos(qJ(1));
t346 = t325 * g(1) - t328 * g(2);
t297 = -qJDD(1) * pkin(1) - (t329 * qJ(2)) + qJDD(2) - t346;
t333 = qJDD(1) * qJ(3) + (qJD(3) * t350) - t297;
t289 = -(t329 * pkin(6)) + t333;
t324 = sin(qJ(4));
t327 = cos(qJ(4));
t343 = qJD(1) * qJD(4);
t304 = -t324 * qJDD(1) - t327 * t343;
t340 = t324 * t343;
t305 = t327 * qJDD(1) - t340;
t345 = qJD(1) * t324;
t306 = -(qJD(4) * mrSges(5,2)) - mrSges(5,3) * t345;
t344 = t327 * qJD(1);
t307 = (qJD(4) * mrSges(5,1)) - mrSges(5,3) * t344;
t308 = (qJD(4) * pkin(4)) - pkin(7) * t344;
t321 = t324 ^ 2;
t273 = t308 * t344 - t304 * pkin(4) + (-pkin(7) * t321 - pkin(6)) * t329 + t333;
t323 = sin(qJ(5));
t326 = cos(qJ(5));
t302 = (-t324 * t323 + t327 * t326) * qJD(1);
t278 = -t302 * qJD(5) + t326 * t304 - t323 * t305;
t301 = (-t327 * t323 - t324 * t326) * qJD(1);
t279 = t301 * qJD(5) + t323 * t304 + t326 * t305;
t316 = qJD(4) + qJD(5);
t294 = -t316 * mrSges(6,2) + t301 * mrSges(6,3);
t295 = t316 * mrSges(6,1) - t302 * mrSges(6,3);
t332 = -m(6) * t273 + t278 * mrSges(6,1) - t279 * mrSges(6,2) + t301 * t294 - t302 * t295;
t349 = -m(4) * t333 - m(5) * t289 + t304 * mrSges(5,1) - t305 * mrSges(5,2) + t332 + (-t306 * t324 - t307 * t327) * qJD(1);
t338 = -t328 * g(1) - t325 * g(2);
t348 = qJDD(1) * qJ(2) + (qJD(2) * t350) + t338;
t293 = qJDD(3) + (-pkin(1) - qJ(3)) * t329 + t348;
t290 = -qJDD(1) * pkin(6) + t293;
t284 = t324 * g(3) + t327 * t290;
t270 = (-t305 - t340) * pkin(7) + (-t324 * t327 * t329 + qJDD(4)) * pkin(4) + t284;
t285 = -t327 * g(3) + t324 * t290;
t271 = -t321 * t329 * pkin(4) + t304 * pkin(7) - qJD(4) * t308 + t285;
t268 = t326 * t270 - t323 * t271;
t286 = -t301 * mrSges(6,1) + t302 * mrSges(6,2);
t315 = qJDD(4) + qJDD(5);
t265 = m(6) * t268 + t315 * mrSges(6,1) - t279 * mrSges(6,3) - t302 * t286 + t316 * t294;
t269 = t323 * t270 + t326 * t271;
t266 = m(6) * t269 - t315 * mrSges(6,2) + t278 * mrSges(6,3) + t301 * t286 - t316 * t295;
t258 = t326 * t265 + t323 * t266;
t303 = (t324 * mrSges(5,1) + t327 * mrSges(5,2)) * qJD(1);
t339 = -t323 * t265 + t326 * t266;
t347 = t324 * (m(5) * t285 - qJDD(4) * mrSges(5,2) + t304 * mrSges(5,3) - qJD(4) * t307 - t303 * t345 + t339) + t327 * (m(5) * t284 + qJDD(4) * mrSges(5,1) - t305 * mrSges(5,3) + qJD(4) * t306 - t303 * t344 + t258);
t335 = m(4) * t293 + qJDD(1) * mrSges(4,2) - (t329 * mrSges(4,3)) + t347;
t281 = Ifges(6,4) * t302 + Ifges(6,2) * t301 + Ifges(6,6) * t316;
t282 = Ifges(6,1) * t302 + Ifges(6,4) * t301 + Ifges(6,5) * t316;
t331 = mrSges(6,1) * t268 - mrSges(6,2) * t269 + Ifges(6,5) * t279 + Ifges(6,6) * t278 + Ifges(6,3) * t315 + t302 * t281 - t301 * t282;
t300 = (Ifges(5,5) * qJD(4)) + (t327 * Ifges(5,1) - t324 * Ifges(5,4)) * qJD(1);
t299 = (Ifges(5,6) * qJD(4)) + (t327 * Ifges(5,4) - t324 * Ifges(5,2)) * qJD(1);
t296 = t329 * pkin(1) - t348;
t280 = Ifges(6,5) * t302 + Ifges(6,6) * t301 + Ifges(6,3) * t316;
t261 = m(3) * t297 + (-mrSges(4,2) - mrSges(3,3)) * t329 + (mrSges(3,2) - mrSges(4,3)) * qJDD(1) + t349;
t260 = mrSges(6,2) * t273 - mrSges(6,3) * t268 + Ifges(6,1) * t279 + Ifges(6,4) * t278 + Ifges(6,5) * t315 + t301 * t280 - t316 * t281;
t259 = -mrSges(6,1) * t273 + mrSges(6,3) * t269 + Ifges(6,4) * t279 + Ifges(6,2) * t278 + Ifges(6,6) * t315 - t302 * t280 + t316 * t282;
t1 = [qJ(2) * (-m(3) * t296 + (t329 * mrSges(3,2)) + t335) - pkin(1) * t261 + mrSges(2,1) * t346 - mrSges(2,2) * t338 + mrSges(3,2) * t297 - mrSges(3,3) * t296 + t327 * (mrSges(5,2) * t289 - mrSges(5,3) * t284 + Ifges(5,1) * t305 + Ifges(5,4) * t304 + Ifges(5,5) * qJDD(4) - pkin(7) * t258 - qJD(4) * t299 - t323 * t259 + t326 * t260) - t324 * (-mrSges(5,1) * t289 + mrSges(5,3) * t285 + Ifges(5,4) * t305 + Ifges(5,2) * t304 + Ifges(5,6) * qJDD(4) + pkin(4) * t332 + pkin(7) * t339 + qJD(4) * t300 + t326 * t259 + t323 * t260) - pkin(6) * t347 + mrSges(4,2) * t293 + mrSges(4,3) * t333 + (qJ(2) * mrSges(3,3) + Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1) + (t329 * mrSges(4,2) + qJDD(1) * mrSges(4,3) - t349) * qJ(3); t261; t335; mrSges(5,1) * t284 - mrSges(5,2) * t285 + Ifges(5,5) * t305 + Ifges(5,6) * t304 + Ifges(5,3) * qJDD(4) + pkin(4) * t258 + (t327 * t299 + t324 * t300) * qJD(1) + t331; t331;];
tauJ = t1;
