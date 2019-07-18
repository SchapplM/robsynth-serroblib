% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRRR1
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-07-18 13:26
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRRR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(1,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR1_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_invdynB_fixb_snew_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR1_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR1_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR1_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:25:01
% EndTime: 2019-07-18 13:25:03
% DurationCPUTime: 0.76s
% Computational Cost: add. (5162->219), mult. (9148->273), div. (0->0), fcn. (6375->8), ass. (0->77)
t307 = mrSges(2,1) + mrSges(3,1);
t306 = -mrSges(2,2) + mrSges(3,3);
t305 = Ifges(3,4) + Ifges(2,5);
t304 = Ifges(2,6) - Ifges(3,6);
t291 = sin(qJ(4));
t295 = cos(qJ(4));
t292 = sin(qJ(3));
t302 = qJD(1) * t292;
t279 = t295 * qJD(3) - t291 * t302;
t296 = cos(qJ(3));
t300 = qJD(1) * qJD(3);
t282 = t292 * qJDD(1) + t296 * t300;
t262 = t279 * qJD(4) + t291 * qJDD(3) + t295 * t282;
t280 = t291 * qJD(3) + t295 * t302;
t301 = t296 * qJD(1);
t288 = qJD(4) - t301;
t290 = sin(qJ(5));
t294 = cos(qJ(5));
t264 = -t290 * t280 + t294 * t288;
t283 = t296 * qJDD(1) - t292 * t300;
t278 = qJDD(4) - t283;
t246 = t264 * qJD(5) + t294 * t262 + t290 * t278;
t293 = sin(qJ(1));
t297 = cos(qJ(1));
t287 = -t297 * g(1) - t293 * g(2);
t272 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t287;
t267 = -t292 * g(3) + t296 * t272;
t286 = t293 * g(1) - t297 * g(2);
t298 = qJD(1) ^ 2;
t277 = -t298 * qJ(2) + qJDD(2) - t286;
t254 = t295 * t267 + t291 * t277;
t266 = t296 * g(3) + t292 * t272;
t247 = -t290 * t254 + t294 * t266;
t265 = t294 * t280 + t290 * t288;
t252 = -t264 * mrSges(6,1) + t265 * mrSges(6,2);
t276 = qJD(5) - t279;
t255 = -t276 * mrSges(6,2) + t264 * mrSges(6,3);
t261 = -t280 * qJD(4) + t295 * qJDD(3) - t291 * t282;
t260 = qJDD(5) - t261;
t243 = m(6) * t247 + t260 * mrSges(6,1) - t246 * mrSges(6,3) - t265 * t252 + t276 * t255;
t245 = -t265 * qJD(5) - t290 * t262 + t294 * t278;
t248 = t294 * t254 + t290 * t266;
t256 = t276 * mrSges(6,1) - t265 * mrSges(6,3);
t244 = m(6) * t248 - t260 * mrSges(6,2) + t245 * mrSges(6,3) + t264 * t252 - t276 * t256;
t263 = -t279 * mrSges(5,1) + t280 * mrSges(5,2);
t269 = t288 * mrSges(5,1) - t280 * mrSges(5,3);
t238 = m(5) * t254 - t278 * mrSges(5,2) + t261 * mrSges(5,3) - t290 * t243 + t294 * t244 + t279 * t263 - t288 * t269;
t253 = t291 * t267 - t295 * t277;
t268 = -t288 * mrSges(5,2) + t279 * mrSges(5,3);
t240 = t278 * mrSges(5,1) + t245 * mrSges(6,1) - t246 * mrSges(6,2) - t262 * mrSges(5,3) + t264 * t255 - t265 * t256 - t280 * t263 + t288 * t268 + (-m(5) - m(6)) * t253;
t281 = (-mrSges(4,1) * t296 + mrSges(4,2) * t292) * qJD(1);
t284 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t302;
t233 = m(4) * t267 - qJDD(3) * mrSges(4,2) + t283 * mrSges(4,3) - qJD(3) * t284 + t295 * t238 - t291 * t240 + t281 * t301;
t285 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t301;
t237 = -t281 * t302 + qJDD(3) * mrSges(4,1) + t261 * mrSges(5,1) - t262 * mrSges(5,2) - t282 * mrSges(4,3) + qJD(3) * t285 - t294 * t243 - t290 * t244 + t279 * t268 - t280 * t269 + (-m(4) - m(5)) * t266;
t303 = t292 * t233 + t296 * t237;
t299 = m(3) * t272 + qJDD(1) * mrSges(3,3) + t296 * t233 - t292 * t237;
t275 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t292 + Ifges(4,4) * t296) * qJD(1);
t274 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t292 + Ifges(4,2) * t296) * qJD(1);
t273 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t292 + Ifges(4,6) * t296) * qJD(1);
t259 = Ifges(5,1) * t280 + Ifges(5,4) * t279 + Ifges(5,5) * t288;
t258 = Ifges(5,4) * t280 + Ifges(5,2) * t279 + Ifges(5,6) * t288;
t257 = Ifges(5,5) * t280 + Ifges(5,6) * t279 + Ifges(5,3) * t288;
t251 = Ifges(6,1) * t265 + Ifges(6,4) * t264 + Ifges(6,5) * t276;
t250 = Ifges(6,4) * t265 + Ifges(6,2) * t264 + Ifges(6,6) * t276;
t249 = Ifges(6,5) * t265 + Ifges(6,6) * t264 + Ifges(6,3) * t276;
t242 = mrSges(6,2) * t253 - mrSges(6,3) * t247 + Ifges(6,1) * t246 + Ifges(6,4) * t245 + Ifges(6,5) * t260 + t264 * t249 - t276 * t250;
t241 = -mrSges(6,1) * t253 + mrSges(6,3) * t248 + Ifges(6,4) * t246 + Ifges(6,2) * t245 + Ifges(6,6) * t260 - t265 * t249 + t276 * t251;
t239 = -mrSges(5,1) * t266 - mrSges(6,1) * t247 + mrSges(6,2) * t248 + mrSges(5,3) * t254 + Ifges(5,4) * t262 - Ifges(6,5) * t246 + Ifges(5,2) * t261 + Ifges(5,6) * t278 - Ifges(6,6) * t245 - Ifges(6,3) * t260 - t265 * t250 + t264 * t251 - t280 * t257 + t288 * t259;
t235 = mrSges(5,2) * t266 + mrSges(5,3) * t253 + Ifges(5,1) * t262 + Ifges(5,4) * t261 + Ifges(5,5) * t278 - t290 * t241 + t294 * t242 + t279 * t257 - t288 * t258;
t234 = Ifges(4,4) * t282 + Ifges(4,2) * t283 + Ifges(4,6) * qJDD(3) - t273 * t302 + qJD(3) * t275 - mrSges(4,1) * t277 + mrSges(4,3) * t267 - Ifges(5,5) * t262 - Ifges(5,6) * t261 - Ifges(5,3) * t278 - t280 * t258 + t279 * t259 + mrSges(5,1) * t253 + mrSges(5,2) * t254 - t290 * t242 - t294 * t241;
t230 = m(2) * t286 + t283 * mrSges(4,1) - t282 * mrSges(4,2) - t291 * t238 - t295 * t240 + t306 * t298 + (-m(3) - m(4)) * t277 + t307 * qJDD(1) + (-t284 * t292 + t285 * t296) * qJD(1);
t229 = mrSges(4,2) * t277 + mrSges(4,3) * t266 + Ifges(4,1) * t282 + Ifges(4,4) * t283 + Ifges(4,5) * qJDD(3) - qJD(3) * t274 + t295 * t235 - t291 * t239 + t273 * t301;
t228 = mrSges(4,1) * t266 + mrSges(3,2) * t272 + mrSges(4,2) * t267 + mrSges(2,3) * t287 - Ifges(4,5) * t282 - Ifges(4,6) * t283 - Ifges(4,3) * qJDD(3) - t291 * t235 - t295 * t239 + t305 * t298 + t304 * qJDD(1) + (-t292 * t274 + t296 * t275) * qJD(1) + t307 * g(3);
t227 = m(2) * t287 - qJDD(1) * mrSges(2,2) - t307 * t298 + t299;
t226 = -mrSges(2,3) * t286 + mrSges(3,2) * t277 + t296 * t229 - t292 * t234 - qJ(2) * t303 - t304 * t298 + t305 * qJDD(1) + (qJ(2) * m(3) + t306) * g(3);
t1 = [-m(1) * g(1) + t297 * t227 - t293 * t230; -m(1) * g(2) + t293 * t227 + t297 * t230; (-m(1) - m(2) - m(3)) * g(3) + t303; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t297 * t226 - t293 * t228; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t293 * t226 + t297 * t228; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t286 - mrSges(2,2) * t287 - mrSges(3,1) * t277 + mrSges(3,3) * t272 + t292 * t229 + t296 * t234 + qJ(2) * (-t298 * mrSges(3,1) + t299) + (Ifges(3,2) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
