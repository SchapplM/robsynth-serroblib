% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RRPR9
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RRPR9_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR9_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR9_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR9_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR9_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR9_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR9_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:09:30
% EndTime: 2019-12-31 17:09:31
% DurationCPUTime: 0.81s
% Computational Cost: add. (5133->199), mult. (11030->260), div. (0->0), fcn. (6922->8), ass. (0->79)
t296 = qJD(1) ^ 2;
t291 = sin(qJ(1));
t294 = cos(qJ(1));
t302 = t291 * g(1) - t294 * g(2);
t271 = -qJDD(1) * pkin(1) - t296 * pkin(5) - t302;
t290 = sin(qJ(2));
t293 = cos(qJ(2));
t304 = qJD(1) * qJD(2);
t303 = t293 * t304;
t279 = t290 * qJDD(1) + t303;
t284 = t290 * t304;
t280 = t293 * qJDD(1) - t284;
t242 = (-t279 - t303) * qJ(3) + (-t280 + t284) * pkin(2) + t271;
t299 = -t294 * g(1) - t291 * g(2);
t272 = -t296 * pkin(1) + qJDD(1) * pkin(5) + t299;
t256 = -t290 * g(3) + t293 * t272;
t277 = (-t293 * pkin(2) - t290 * qJ(3)) * qJD(1);
t295 = qJD(2) ^ 2;
t305 = t293 * qJD(1);
t245 = -t295 * pkin(2) + qJDD(2) * qJ(3) + t277 * t305 + t256;
t287 = sin(pkin(7));
t288 = cos(pkin(7));
t306 = qJD(1) * t290;
t275 = t287 * qJD(2) + t288 * t306;
t229 = -0.2e1 * qJD(3) * t275 + t288 * t242 - t287 * t245;
t261 = t287 * qJDD(2) + t288 * t279;
t274 = t288 * qJD(2) - t287 * t306;
t227 = (-t274 * t305 - t261) * pkin(6) + (t274 * t275 - t280) * pkin(3) + t229;
t230 = 0.2e1 * qJD(3) * t274 + t287 * t242 + t288 * t245;
t260 = t288 * qJDD(2) - t287 * t279;
t262 = -pkin(3) * t305 - t275 * pkin(6);
t273 = t274 ^ 2;
t228 = -t273 * pkin(3) + t260 * pkin(6) + t262 * t305 + t230;
t289 = sin(qJ(4));
t292 = cos(qJ(4));
t225 = t292 * t227 - t289 * t228;
t252 = t292 * t274 - t289 * t275;
t234 = t252 * qJD(4) + t289 * t260 + t292 * t261;
t253 = t289 * t274 + t292 * t275;
t239 = -t252 * mrSges(5,1) + t253 * mrSges(5,2);
t283 = qJD(4) - t305;
t246 = -t283 * mrSges(5,2) + t252 * mrSges(5,3);
t276 = qJDD(4) - t280;
t222 = m(5) * t225 + t276 * mrSges(5,1) - t234 * mrSges(5,3) - t253 * t239 + t283 * t246;
t226 = t289 * t227 + t292 * t228;
t233 = -t253 * qJD(4) + t292 * t260 - t289 * t261;
t247 = t283 * mrSges(5,1) - t253 * mrSges(5,3);
t223 = m(5) * t226 - t276 * mrSges(5,2) + t233 * mrSges(5,3) + t252 * t239 - t283 * t247;
t216 = t292 * t222 + t289 * t223;
t255 = -t293 * g(3) - t290 * t272;
t254 = -t274 * mrSges(4,1) + t275 * mrSges(4,2);
t258 = mrSges(4,2) * t305 + t274 * mrSges(4,3);
t214 = m(4) * t229 - t280 * mrSges(4,1) - t261 * mrSges(4,3) - t275 * t254 - t258 * t305 + t216;
t259 = -mrSges(4,1) * t305 - t275 * mrSges(4,3);
t300 = -t289 * t222 + t292 * t223;
t215 = m(4) * t230 + t280 * mrSges(4,2) + t260 * mrSges(4,3) + t274 * t254 + t259 * t305 + t300;
t301 = -t287 * t214 + t288 * t215;
t212 = t288 * t214 + t287 * t215;
t244 = -qJDD(2) * pkin(2) - t295 * qJ(3) + t277 * t306 + qJDD(3) - t255;
t231 = -t260 * pkin(3) - t273 * pkin(6) + t275 * t262 + t244;
t298 = m(5) * t231 - t233 * mrSges(5,1) + t234 * mrSges(5,2) - t252 * t246 + t253 * t247;
t236 = Ifges(5,4) * t253 + Ifges(5,2) * t252 + Ifges(5,6) * t283;
t237 = Ifges(5,1) * t253 + Ifges(5,4) * t252 + Ifges(5,5) * t283;
t297 = mrSges(5,1) * t225 - mrSges(5,2) * t226 + Ifges(5,5) * t234 + Ifges(5,6) * t233 + Ifges(5,3) * t276 + t253 * t236 - t252 * t237;
t224 = m(4) * t244 - t260 * mrSges(4,1) + t261 * mrSges(4,2) - t274 * t258 + t275 * t259 + t298;
t282 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t305;
t281 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t306;
t278 = (-t293 * mrSges(3,1) + t290 * mrSges(3,2)) * qJD(1);
t270 = Ifges(3,5) * qJD(2) + (t290 * Ifges(3,1) + t293 * Ifges(3,4)) * qJD(1);
t269 = Ifges(3,6) * qJD(2) + (t290 * Ifges(3,4) + Ifges(3,2) * t293) * qJD(1);
t250 = Ifges(4,1) * t275 + Ifges(4,4) * t274 - Ifges(4,5) * t305;
t249 = Ifges(4,4) * t275 + Ifges(4,2) * t274 - Ifges(4,6) * t305;
t248 = Ifges(4,5) * t275 + Ifges(4,6) * t274 - Ifges(4,3) * t305;
t235 = Ifges(5,5) * t253 + Ifges(5,6) * t252 + Ifges(5,3) * t283;
t218 = mrSges(5,2) * t231 - mrSges(5,3) * t225 + Ifges(5,1) * t234 + Ifges(5,4) * t233 + Ifges(5,5) * t276 + t252 * t235 - t283 * t236;
t217 = -mrSges(5,1) * t231 + mrSges(5,3) * t226 + Ifges(5,4) * t234 + Ifges(5,2) * t233 + Ifges(5,6) * t276 - t253 * t235 + t283 * t237;
t211 = mrSges(4,2) * t244 - mrSges(4,3) * t229 + Ifges(4,1) * t261 + Ifges(4,4) * t260 - Ifges(4,5) * t280 - pkin(6) * t216 - t289 * t217 + t292 * t218 + t274 * t248 + t249 * t305;
t210 = -mrSges(4,1) * t244 + mrSges(4,3) * t230 + Ifges(4,4) * t261 + Ifges(4,2) * t260 - Ifges(4,6) * t280 - pkin(3) * t298 + pkin(6) * t300 + t292 * t217 + t289 * t218 - t275 * t248 - t250 * t305;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t302 - mrSges(2,2) * t299 + t290 * (mrSges(3,2) * t271 - mrSges(3,3) * t255 + Ifges(3,1) * t279 + Ifges(3,4) * t280 + Ifges(3,5) * qJDD(2) - qJ(3) * t212 - qJD(2) * t269 - t287 * t210 + t288 * t211) + t293 * (-mrSges(3,1) * t271 - mrSges(4,1) * t229 + mrSges(4,2) * t230 + mrSges(3,3) * t256 + Ifges(3,4) * t279 - Ifges(4,5) * t261 + Ifges(3,6) * qJDD(2) - Ifges(4,6) * t260 - pkin(2) * t212 - pkin(3) * t216 + qJD(2) * t270 - t275 * t249 + t274 * t250 - t297 + (Ifges(3,2) + Ifges(4,3)) * t280) + pkin(1) * (-m(3) * t271 + t280 * mrSges(3,1) - t279 * mrSges(3,2) + (-t281 * t290 + t282 * t293) * qJD(1) - t212) + pkin(5) * (t293 * (m(3) * t256 - qJDD(2) * mrSges(3,2) + t280 * mrSges(3,3) - qJD(2) * t281 + t278 * t305 + t301) - t290 * (m(3) * t255 + qJDD(2) * mrSges(3,1) - t279 * mrSges(3,3) + qJD(2) * t282 - t278 * t306 - t224)); Ifges(3,5) * t279 + Ifges(3,6) * t280 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t255 - mrSges(3,2) * t256 + t287 * t211 + t288 * t210 - pkin(2) * t224 + qJ(3) * t301 + (t290 * t269 - t293 * t270) * qJD(1); t224; t297;];
tauJ = t1;
