% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRRR6
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRRR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR6_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR6_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR6_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR6_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR6_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:09:46
% EndTime: 2019-12-05 17:09:47
% DurationCPUTime: 0.55s
% Computational Cost: add. (5057->160), mult. (6402->213), div. (0->0), fcn. (3939->10), ass. (0->74)
t304 = qJD(2) + qJD(3);
t310 = sin(qJ(4));
t329 = t304 * t310;
t314 = cos(qJ(4));
t328 = t304 * t314;
t307 = sin(pkin(9));
t308 = cos(pkin(9));
t295 = -t308 * g(1) - t307 * g(2);
t306 = -g(3) + qJDD(1);
t312 = sin(qJ(2));
t316 = cos(qJ(2));
t276 = -t312 * t295 + t316 * t306;
t274 = qJDD(2) * pkin(2) + t276;
t277 = t316 * t295 + t312 * t306;
t317 = qJD(2) ^ 2;
t275 = -t317 * pkin(2) + t277;
t311 = sin(qJ(3));
t315 = cos(qJ(3));
t263 = t311 * t274 + t315 * t275;
t300 = t304 ^ 2;
t302 = qJDD(2) + qJDD(3);
t255 = -t300 * pkin(3) + t302 * pkin(7) + t263;
t294 = -t307 * g(1) + t308 * g(2);
t250 = -t310 * t255 + t314 * t294;
t326 = qJD(4) * t304;
t325 = t314 * t326;
t286 = t310 * t302 + t325;
t247 = (-t286 + t325) * pkin(8) + (t300 * t310 * t314 + qJDD(4)) * pkin(4) + t250;
t251 = t314 * t255 + t310 * t294;
t287 = t314 * t302 - t310 * t326;
t293 = qJD(4) * pkin(4) - pkin(8) * t329;
t305 = t314 ^ 2;
t248 = -t305 * t300 * pkin(4) + t287 * pkin(8) - qJD(4) * t293 + t251;
t309 = sin(qJ(5));
t313 = cos(qJ(5));
t245 = t313 * t247 - t309 * t248;
t281 = (-t309 * t310 + t313 * t314) * t304;
t260 = t281 * qJD(5) + t313 * t286 + t309 * t287;
t282 = (t309 * t314 + t310 * t313) * t304;
t268 = -t281 * mrSges(6,1) + t282 * mrSges(6,2);
t303 = qJD(4) + qJD(5);
t269 = -t303 * mrSges(6,2) + t281 * mrSges(6,3);
t301 = qJDD(4) + qJDD(5);
t242 = m(6) * t245 + t301 * mrSges(6,1) - t260 * mrSges(6,3) - t282 * t268 + t303 * t269;
t246 = t309 * t247 + t313 * t248;
t259 = -t282 * qJD(5) - t309 * t286 + t313 * t287;
t270 = t303 * mrSges(6,1) - t282 * mrSges(6,3);
t243 = m(6) * t246 - t301 * mrSges(6,2) + t259 * mrSges(6,3) + t281 * t268 - t303 * t270;
t233 = t313 * t242 + t309 * t243;
t285 = (-mrSges(5,1) * t314 + mrSges(5,2) * t310) * t304;
t291 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t329;
t292 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t328;
t323 = -t309 * t242 + t313 * t243;
t324 = -t310 * (m(5) * t250 + qJDD(4) * mrSges(5,1) - t286 * mrSges(5,3) + qJD(4) * t292 - t285 * t329 + t233) + t314 * (m(5) * t251 - qJDD(4) * mrSges(5,2) + t287 * mrSges(5,3) - qJD(4) * t291 + t285 * t328 + t323);
t229 = m(4) * t263 - t300 * mrSges(4,1) - t302 * mrSges(4,2) + t324;
t262 = t315 * t274 - t311 * t275;
t322 = -t302 * pkin(3) - t262;
t254 = -t300 * pkin(7) + t322;
t249 = t293 * t329 - t287 * pkin(4) + (-pkin(8) * t305 - pkin(7)) * t300 + t322;
t320 = m(6) * t249 - t259 * mrSges(6,1) + t260 * mrSges(6,2) - t281 * t269 + t282 * t270;
t318 = -m(5) * t254 + t287 * mrSges(5,1) - t286 * mrSges(5,2) - t291 * t329 + t292 * t328 - t320;
t237 = m(4) * t262 + t302 * mrSges(4,1) - t300 * mrSges(4,2) + t318;
t327 = t311 * t229 + t315 * t237;
t264 = Ifges(6,5) * t282 + Ifges(6,6) * t281 + Ifges(6,3) * t303;
t266 = Ifges(6,1) * t282 + Ifges(6,4) * t281 + Ifges(6,5) * t303;
t234 = -mrSges(6,1) * t249 + mrSges(6,3) * t246 + Ifges(6,4) * t260 + Ifges(6,2) * t259 + Ifges(6,6) * t301 - t282 * t264 + t303 * t266;
t265 = Ifges(6,4) * t282 + Ifges(6,2) * t281 + Ifges(6,6) * t303;
t235 = mrSges(6,2) * t249 - mrSges(6,3) * t245 + Ifges(6,1) * t260 + Ifges(6,4) * t259 + Ifges(6,5) * t301 + t281 * t264 - t303 * t265;
t278 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t310 + Ifges(5,6) * t314) * t304;
t279 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t310 + Ifges(5,2) * t314) * t304;
t280 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t310 + Ifges(5,4) * t314) * t304;
t321 = -mrSges(4,2) * t263 + t314 * (-mrSges(5,1) * t254 + mrSges(5,3) * t251 + Ifges(5,4) * t286 + Ifges(5,2) * t287 + Ifges(5,6) * qJDD(4) - pkin(4) * t320 + pkin(8) * t323 + qJD(4) * t280 + t313 * t234 + t309 * t235 - t278 * t329) + t310 * (mrSges(5,2) * t254 - mrSges(5,3) * t250 + Ifges(5,1) * t286 + Ifges(5,4) * t287 + Ifges(5,5) * qJDD(4) - pkin(8) * t233 - qJD(4) * t279 - t309 * t234 + t313 * t235 + t278 * t328) + pkin(7) * t324 + pkin(3) * t318 + mrSges(4,1) * t262 + Ifges(4,3) * t302;
t319 = mrSges(6,1) * t245 - mrSges(6,2) * t246 + Ifges(6,5) * t260 + Ifges(6,6) * t259 + Ifges(6,3) * t301 + t282 * t265 - t281 * t266;
t1 = [m(2) * t306 + t312 * (m(3) * t277 - t317 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t315 * t229 - t311 * t237) + t316 * (m(3) * t276 + qJDD(2) * mrSges(3,1) - t317 * mrSges(3,2) + t327); mrSges(3,1) * t276 - mrSges(3,2) * t277 + Ifges(3,3) * qJDD(2) + pkin(2) * t327 + t321; t321; mrSges(5,1) * t250 - mrSges(5,2) * t251 + Ifges(5,5) * t286 + Ifges(5,6) * t287 + Ifges(5,3) * qJDD(4) + pkin(4) * t233 + (t279 * t310 - t280 * t314) * t304 + t319; t319;];
tauJ = t1;
