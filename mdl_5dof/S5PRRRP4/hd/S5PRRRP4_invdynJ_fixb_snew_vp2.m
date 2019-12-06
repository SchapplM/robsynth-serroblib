% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRRP4
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRRP4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP4_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP4_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP4_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP4_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP4_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:45:45
% EndTime: 2019-12-05 16:45:45
% DurationCPUTime: 0.46s
% Computational Cost: add. (2073->136), mult. (2559->172), div. (0->0), fcn. (1338->8), ass. (0->63)
t321 = Ifges(5,1) + Ifges(6,1);
t317 = Ifges(5,4) - Ifges(6,5);
t316 = Ifges(5,5) + Ifges(6,4);
t320 = Ifges(5,2) + Ifges(6,3);
t315 = Ifges(5,6) - Ifges(6,6);
t319 = Ifges(5,3) + Ifges(6,2);
t318 = mrSges(5,3) + mrSges(6,2);
t288 = qJD(2) + qJD(3);
t294 = sin(qJ(4));
t314 = t288 * t294;
t297 = cos(qJ(4));
t313 = t288 * t297;
t292 = sin(pkin(8));
t293 = cos(pkin(8));
t280 = -t292 * g(1) + t293 * g(2);
t312 = t297 * t280;
t281 = -t293 * g(1) - t292 * g(2);
t291 = -g(3) + qJDD(1);
t296 = sin(qJ(2));
t299 = cos(qJ(2));
t254 = -t296 * t281 + t299 * t291;
t251 = qJDD(2) * pkin(2) + t254;
t255 = t299 * t281 + t296 * t291;
t301 = qJD(2) ^ 2;
t252 = -t301 * pkin(2) + t255;
t295 = sin(qJ(3));
t298 = cos(qJ(3));
t247 = t295 * t251 + t298 * t252;
t286 = t288 ^ 2;
t287 = qJDD(2) + qJDD(3);
t244 = -t286 * pkin(3) + t287 * pkin(7) + t247;
t240 = -t294 * t244 + t312;
t241 = t297 * t244 + t294 * t280;
t268 = (-mrSges(6,1) * t297 - mrSges(6,3) * t294) * t288;
t269 = (-mrSges(5,1) * t297 + mrSges(5,2) * t294) * t288;
t307 = qJD(4) * t288;
t270 = t294 * t287 + t297 * t307;
t271 = t297 * t287 - t294 * t307;
t276 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t314;
t278 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t313;
t267 = (-pkin(4) * t297 - qJ(5) * t294) * t288;
t300 = qJD(4) ^ 2;
t239 = -qJDD(4) * pkin(4) - t300 * qJ(5) - t312 + qJDD(5) + (t267 * t288 + t244) * t294;
t279 = mrSges(6,2) * t313 + qJD(4) * mrSges(6,3);
t304 = -m(6) * t239 + qJDD(4) * mrSges(6,1) + qJD(4) * t279;
t238 = -t300 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t267 * t313 + t241;
t277 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t314;
t305 = m(6) * t238 + qJDD(4) * mrSges(6,3) + qJD(4) * t277 + t268 * t313;
t306 = -t294 * (m(5) * t240 + qJDD(4) * mrSges(5,1) + qJD(4) * t278 + (-t268 - t269) * t314 - t318 * t270 + t304) + t297 * (m(5) * t241 - qJDD(4) * mrSges(5,2) - qJD(4) * t276 + t269 * t313 + t318 * t271 + t305);
t227 = m(4) * t247 - t286 * mrSges(4,1) - t287 * mrSges(4,2) + t306;
t246 = t298 * t251 - t295 * t252;
t243 = -t287 * pkin(3) - t286 * pkin(7) - t246;
t236 = -t271 * pkin(4) - t270 * qJ(5) + (-0.2e1 * qJD(5) * t294 + (pkin(4) * t294 - qJ(5) * t297) * qJD(4)) * t288 + t243;
t234 = m(6) * t236 - t271 * mrSges(6,1) - t270 * mrSges(6,3) - t277 * t314 - t279 * t313;
t302 = -m(5) * t243 + t271 * mrSges(5,1) - t270 * mrSges(5,2) - t276 * t314 + t278 * t313 - t234;
t230 = m(4) * t246 + t287 * mrSges(4,1) - t286 * mrSges(4,2) + t302;
t311 = t295 * t227 + t298 * t230;
t310 = (-t294 * t317 - t320 * t297) * t288 - t315 * qJD(4);
t309 = (t294 * t316 + t297 * t315) * t288 + t319 * qJD(4);
t308 = (t321 * t294 + t297 * t317) * t288 + t316 * qJD(4);
t303 = -mrSges(4,2) * t247 + t297 * (-mrSges(5,1) * t243 - mrSges(6,1) * t236 + mrSges(6,2) * t238 + mrSges(5,3) * t241 - pkin(4) * t234 + t308 * qJD(4) + t315 * qJDD(4) + t317 * t270 + t320 * t271 - t309 * t314) + t294 * (mrSges(5,2) * t243 + mrSges(6,2) * t239 - mrSges(5,3) * t240 - mrSges(6,3) * t236 - qJ(5) * t234 + t310 * qJD(4) + t316 * qJDD(4) + t321 * t270 + t317 * t271 + t309 * t313) + pkin(7) * t306 + pkin(3) * t302 + mrSges(4,1) * t246 + Ifges(4,3) * t287;
t235 = t270 * mrSges(6,2) + t268 * t314 - t304;
t1 = [m(2) * t291 + t296 * (m(3) * t255 - t301 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t298 * t227 - t295 * t230) + t299 * (m(3) * t254 + qJDD(2) * mrSges(3,1) - t301 * mrSges(3,2) + t311); mrSges(3,1) * t254 - mrSges(3,2) * t255 + Ifges(3,3) * qJDD(2) + pkin(2) * t311 + t303; t303; mrSges(5,1) * t240 - mrSges(5,2) * t241 - mrSges(6,1) * t239 + mrSges(6,3) * t238 - pkin(4) * t235 + qJ(5) * t305 + (qJ(5) * mrSges(6,2) + t315) * t271 + t316 * t270 + t319 * qJDD(4) + (-t310 * t294 - t308 * t297) * t288; t235;];
tauJ = t1;
