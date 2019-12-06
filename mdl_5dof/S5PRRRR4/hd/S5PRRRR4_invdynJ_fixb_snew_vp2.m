% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRRR4
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
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRRR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR4_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR4_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR4_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR4_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR4_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:51
% EndTime: 2019-12-05 17:07:52
% DurationCPUTime: 0.48s
% Computational Cost: add. (4415->155), mult. (5878->205), div. (0->0), fcn. (3729->10), ass. (0->72)
t298 = qJD(2) + qJD(3);
t304 = sin(qJ(4));
t323 = t298 * t304;
t308 = cos(qJ(4));
t322 = t298 * t308;
t301 = sin(pkin(9));
t302 = cos(pkin(9));
t288 = t301 * g(1) - t302 * g(2);
t289 = -t302 * g(1) - t301 * g(2);
t306 = sin(qJ(2));
t310 = cos(qJ(2));
t316 = t310 * t288 - t306 * t289;
t268 = qJDD(2) * pkin(2) + t316;
t321 = t306 * t288 + t310 * t289;
t269 = -qJD(2) ^ 2 * pkin(2) + t321;
t305 = sin(qJ(3));
t309 = cos(qJ(3));
t254 = t305 * t268 + t309 * t269;
t294 = t298 ^ 2;
t296 = qJDD(2) + qJDD(3);
t250 = -t294 * pkin(3) + t296 * pkin(7) + t254;
t300 = -g(3) + qJDD(1);
t246 = -t304 * t250 + t308 * t300;
t320 = qJD(4) * t298;
t319 = t308 * t320;
t280 = t304 * t296 + t319;
t243 = (-t280 + t319) * pkin(8) + (t294 * t304 * t308 + qJDD(4)) * pkin(4) + t246;
t247 = t308 * t250 + t304 * t300;
t281 = t308 * t296 - t304 * t320;
t287 = qJD(4) * pkin(4) - pkin(8) * t323;
t299 = t308 ^ 2;
t244 = -t299 * t294 * pkin(4) + t281 * pkin(8) - qJD(4) * t287 + t247;
t303 = sin(qJ(5));
t307 = cos(qJ(5));
t241 = t307 * t243 - t303 * t244;
t275 = (-t303 * t304 + t307 * t308) * t298;
t259 = t275 * qJD(5) + t307 * t280 + t303 * t281;
t276 = (t303 * t308 + t304 * t307) * t298;
t264 = -t275 * mrSges(6,1) + t276 * mrSges(6,2);
t297 = qJD(4) + qJD(5);
t270 = -t297 * mrSges(6,2) + t275 * mrSges(6,3);
t295 = qJDD(4) + qJDD(5);
t238 = m(6) * t241 + t295 * mrSges(6,1) - t259 * mrSges(6,3) - t276 * t264 + t297 * t270;
t242 = t303 * t243 + t307 * t244;
t258 = -t276 * qJD(5) - t303 * t280 + t307 * t281;
t271 = t297 * mrSges(6,1) - t276 * mrSges(6,3);
t239 = m(6) * t242 - t295 * mrSges(6,2) + t258 * mrSges(6,3) + t275 * t264 - t297 * t271;
t231 = t307 * t238 + t303 * t239;
t279 = (-mrSges(5,1) * t308 + mrSges(5,2) * t304) * t298;
t286 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t322;
t229 = m(5) * t246 + qJDD(4) * mrSges(5,1) - t280 * mrSges(5,3) + qJD(4) * t286 - t279 * t323 + t231;
t285 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t323;
t317 = -t303 * t238 + t307 * t239;
t230 = m(5) * t247 - qJDD(4) * mrSges(5,2) + t281 * mrSges(5,3) - qJD(4) * t285 + t279 * t322 + t317;
t318 = -t304 * t229 + t308 * t230;
t253 = t309 * t268 - t305 * t269;
t315 = -t296 * pkin(3) - t253;
t245 = t287 * t323 - t281 * pkin(4) + (-pkin(8) * t299 - pkin(7)) * t294 + t315;
t260 = Ifges(6,5) * t276 + Ifges(6,6) * t275 + Ifges(6,3) * t297;
t262 = Ifges(6,1) * t276 + Ifges(6,4) * t275 + Ifges(6,5) * t297;
t232 = -mrSges(6,1) * t245 + mrSges(6,3) * t242 + Ifges(6,4) * t259 + Ifges(6,2) * t258 + Ifges(6,6) * t295 - t276 * t260 + t297 * t262;
t261 = Ifges(6,4) * t276 + Ifges(6,2) * t275 + Ifges(6,6) * t297;
t233 = mrSges(6,2) * t245 - mrSges(6,3) * t241 + Ifges(6,1) * t259 + Ifges(6,4) * t258 + Ifges(6,5) * t295 + t275 * t260 - t297 * t261;
t249 = -t294 * pkin(7) + t315;
t272 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t304 + Ifges(5,6) * t308) * t298;
t273 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t304 + Ifges(5,2) * t308) * t298;
t274 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t304 + Ifges(5,4) * t308) * t298;
t313 = m(6) * t245 - t258 * mrSges(6,1) + t259 * mrSges(6,2) - t275 * t270 + t276 * t271;
t311 = -m(5) * t249 + t281 * mrSges(5,1) - t280 * mrSges(5,2) - t285 * t323 + t286 * t322 - t313;
t314 = -mrSges(4,2) * t254 + t308 * (-mrSges(5,1) * t249 + mrSges(5,3) * t247 + Ifges(5,4) * t280 + Ifges(5,2) * t281 + Ifges(5,6) * qJDD(4) - pkin(4) * t313 + pkin(8) * t317 + qJD(4) * t274 + t307 * t232 + t303 * t233 - t272 * t323) + t304 * (mrSges(5,2) * t249 - mrSges(5,3) * t246 + Ifges(5,1) * t280 + Ifges(5,4) * t281 + Ifges(5,5) * qJDD(4) - pkin(8) * t231 - qJD(4) * t273 - t303 * t232 + t307 * t233 + t272 * t322) + pkin(7) * t318 + pkin(3) * t311 + mrSges(4,1) * t253 + Ifges(4,3) * t296;
t312 = mrSges(6,1) * t241 - mrSges(6,2) * t242 + Ifges(6,5) * t259 + Ifges(6,6) * t258 + Ifges(6,3) * t295 + t276 * t261 - t275 * t262;
t1 = [t308 * t229 + t304 * t230 + (m(2) + m(3) + m(4)) * t300; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t316 - mrSges(3,2) * t321 + pkin(2) * (t305 * (m(4) * t254 - t294 * mrSges(4,1) - t296 * mrSges(4,2) + t318) + t309 * (m(4) * t253 + t296 * mrSges(4,1) - t294 * mrSges(4,2) + t311)) + t314; t314; mrSges(5,1) * t246 - mrSges(5,2) * t247 + Ifges(5,5) * t280 + Ifges(5,6) * t281 + Ifges(5,3) * qJDD(4) + pkin(4) * t231 + (t273 * t304 - t274 * t308) * t298 + t312; t312;];
tauJ = t1;
