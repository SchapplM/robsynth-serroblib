% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RRRR5
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% tauJ [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RRRR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR5_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR5_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR5_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR5_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR5_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:27:32
% EndTime: 2019-12-31 17:27:33
% DurationCPUTime: 0.81s
% Computational Cost: add. (6132->201), mult. (12216->259), div. (0->0), fcn. (7782->8), ass. (0->82)
t300 = qJD(1) ^ 2;
t294 = sin(qJ(1));
t298 = cos(qJ(1));
t308 = t294 * g(1) - t298 * g(2);
t273 = -qJDD(1) * pkin(1) - t300 * pkin(5) - t308;
t293 = sin(qJ(2));
t297 = cos(qJ(2));
t310 = qJD(1) * qJD(2);
t309 = t297 * t310;
t282 = t293 * qJDD(1) + t309;
t288 = t293 * t310;
t283 = t297 * qJDD(1) - t288;
t245 = (-t282 - t309) * pkin(6) + (-t283 + t288) * pkin(2) + t273;
t305 = -t298 * g(1) - t294 * g(2);
t274 = -t300 * pkin(1) + qJDD(1) * pkin(5) + t305;
t264 = -t293 * g(3) + t297 * t274;
t281 = (-t297 * pkin(2) - t293 * pkin(6)) * qJD(1);
t299 = qJD(2) ^ 2;
t311 = t297 * qJD(1);
t248 = -t299 * pkin(2) + qJDD(2) * pkin(6) + t281 * t311 + t264;
t292 = sin(qJ(3));
t296 = cos(qJ(3));
t236 = t296 * t245 - t292 * t248;
t312 = qJD(1) * t293;
t278 = t296 * qJD(2) - t292 * t312;
t257 = t278 * qJD(3) + t292 * qJDD(2) + t296 * t282;
t277 = qJDD(3) - t283;
t279 = t292 * qJD(2) + t296 * t312;
t287 = qJD(3) - t311;
t227 = (t278 * t287 - t257) * pkin(7) + (t278 * t279 + t277) * pkin(3) + t236;
t237 = t292 * t245 + t296 * t248;
t256 = -t279 * qJD(3) + t296 * qJDD(2) - t292 * t282;
t265 = t287 * pkin(3) - t279 * pkin(7);
t276 = t278 ^ 2;
t228 = -t276 * pkin(3) + t256 * pkin(7) - t287 * t265 + t237;
t291 = sin(qJ(4));
t295 = cos(qJ(4));
t225 = t295 * t227 - t291 * t228;
t258 = t295 * t278 - t291 * t279;
t235 = t258 * qJD(4) + t291 * t256 + t295 * t257;
t259 = t291 * t278 + t295 * t279;
t242 = -t258 * mrSges(5,1) + t259 * mrSges(5,2);
t286 = qJD(4) + t287;
t249 = -t286 * mrSges(5,2) + t258 * mrSges(5,3);
t275 = qJDD(4) + t277;
t222 = m(5) * t225 + t275 * mrSges(5,1) - t235 * mrSges(5,3) - t259 * t242 + t286 * t249;
t226 = t291 * t227 + t295 * t228;
t234 = -t259 * qJD(4) + t295 * t256 - t291 * t257;
t250 = t286 * mrSges(5,1) - t259 * mrSges(5,3);
t223 = m(5) * t226 - t275 * mrSges(5,2) + t234 * mrSges(5,3) + t258 * t242 - t286 * t250;
t216 = t295 * t222 + t291 * t223;
t263 = -t297 * g(3) - t293 * t274;
t260 = -t278 * mrSges(4,1) + t279 * mrSges(4,2);
t261 = -t287 * mrSges(4,2) + t278 * mrSges(4,3);
t214 = m(4) * t236 + t277 * mrSges(4,1) - t257 * mrSges(4,3) - t279 * t260 + t287 * t261 + t216;
t262 = t287 * mrSges(4,1) - t279 * mrSges(4,3);
t306 = -t291 * t222 + t295 * t223;
t215 = m(4) * t237 - t277 * mrSges(4,2) + t256 * mrSges(4,3) + t278 * t260 - t287 * t262 + t306;
t307 = -t292 * t214 + t296 * t215;
t212 = t296 * t214 + t292 * t215;
t247 = -qJDD(2) * pkin(2) - t299 * pkin(6) + t281 * t312 - t263;
t229 = -t256 * pkin(3) - t276 * pkin(7) + t279 * t265 + t247;
t304 = m(5) * t229 - t234 * mrSges(5,1) + t235 * mrSges(5,2) - t258 * t249 + t259 * t250;
t239 = Ifges(5,4) * t259 + Ifges(5,2) * t258 + Ifges(5,6) * t286;
t240 = Ifges(5,1) * t259 + Ifges(5,4) * t258 + Ifges(5,5) * t286;
t303 = -mrSges(5,1) * t225 + mrSges(5,2) * t226 - Ifges(5,5) * t235 - Ifges(5,6) * t234 - Ifges(5,3) * t275 - t259 * t239 + t258 * t240;
t302 = -m(4) * t247 + t256 * mrSges(4,1) - t257 * mrSges(4,2) + t278 * t261 - t279 * t262 - t304;
t252 = Ifges(4,4) * t279 + Ifges(4,2) * t278 + Ifges(4,6) * t287;
t253 = Ifges(4,1) * t279 + Ifges(4,4) * t278 + Ifges(4,5) * t287;
t301 = mrSges(4,1) * t236 - mrSges(4,2) * t237 + Ifges(4,5) * t257 + Ifges(4,6) * t256 + Ifges(4,3) * t277 + pkin(3) * t216 + t279 * t252 - t278 * t253 - t303;
t285 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t311;
t284 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t312;
t280 = (-t297 * mrSges(3,1) + t293 * mrSges(3,2)) * qJD(1);
t272 = Ifges(3,5) * qJD(2) + (t293 * Ifges(3,1) + t297 * Ifges(3,4)) * qJD(1);
t271 = Ifges(3,6) * qJD(2) + (t293 * Ifges(3,4) + t297 * Ifges(3,2)) * qJD(1);
t251 = Ifges(4,5) * t279 + Ifges(4,6) * t278 + Ifges(4,3) * t287;
t238 = Ifges(5,5) * t259 + Ifges(5,6) * t258 + Ifges(5,3) * t286;
t218 = mrSges(5,2) * t229 - mrSges(5,3) * t225 + Ifges(5,1) * t235 + Ifges(5,4) * t234 + Ifges(5,5) * t275 + t258 * t238 - t286 * t239;
t217 = -mrSges(5,1) * t229 + mrSges(5,3) * t226 + Ifges(5,4) * t235 + Ifges(5,2) * t234 + Ifges(5,6) * t275 - t259 * t238 + t286 * t240;
t211 = mrSges(4,2) * t247 - mrSges(4,3) * t236 + Ifges(4,1) * t257 + Ifges(4,4) * t256 + Ifges(4,5) * t277 - pkin(7) * t216 - t291 * t217 + t295 * t218 + t278 * t251 - t287 * t252;
t210 = -mrSges(4,1) * t247 + mrSges(4,3) * t237 + Ifges(4,4) * t257 + Ifges(4,2) * t256 + Ifges(4,6) * t277 - pkin(3) * t304 + pkin(7) * t306 + t295 * t217 + t291 * t218 - t279 * t251 + t287 * t253;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t308 - mrSges(2,2) * t305 + t293 * (mrSges(3,2) * t273 - mrSges(3,3) * t263 + Ifges(3,1) * t282 + Ifges(3,4) * t283 + Ifges(3,5) * qJDD(2) - pkin(6) * t212 - qJD(2) * t271 - t292 * t210 + t296 * t211) + t297 * (-mrSges(3,1) * t273 + mrSges(3,3) * t264 + Ifges(3,4) * t282 + Ifges(3,2) * t283 + Ifges(3,6) * qJDD(2) - pkin(2) * t212 + qJD(2) * t272 - t301) + pkin(1) * (-m(3) * t273 + t283 * mrSges(3,1) - t282 * mrSges(3,2) + (-t284 * t293 + t285 * t297) * qJD(1) - t212) + pkin(5) * (t297 * (m(3) * t264 - qJDD(2) * mrSges(3,2) + t283 * mrSges(3,3) - qJD(2) * t284 + t280 * t311 + t307) - t293 * (m(3) * t263 + qJDD(2) * mrSges(3,1) - t282 * mrSges(3,3) + qJD(2) * t285 - t280 * t312 + t302)); Ifges(3,5) * t282 + Ifges(3,6) * t283 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t263 - mrSges(3,2) * t264 + t292 * t211 + t296 * t210 + pkin(2) * t302 + pkin(6) * t307 + (t293 * t271 - t297 * t272) * qJD(1); t301; -t303;];
tauJ = t1;
