% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PPRRR2
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PPRRR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR2_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR2_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR2_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR2_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR2_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:22
% EndTime: 2019-12-05 15:14:22
% DurationCPUTime: 0.41s
% Computational Cost: add. (2239->148), mult. (4081->199), div. (0->0), fcn. (2686->10), ass. (0->68)
t278 = sin(pkin(8));
t280 = cos(pkin(8));
t268 = -g(1) * t280 - g(2) * t278;
t276 = -g(3) + qJDD(1);
t277 = sin(pkin(9));
t279 = cos(pkin(9));
t253 = -t268 * t277 + t276 * t279;
t254 = t268 * t279 + t276 * t277;
t283 = sin(qJ(3));
t286 = cos(qJ(3));
t240 = t283 * t253 + t286 * t254;
t287 = qJD(3) ^ 2;
t235 = -pkin(3) * t287 + qJDD(3) * pkin(6) + t240;
t267 = -g(1) * t278 + g(2) * t280 + qJDD(2);
t282 = sin(qJ(4));
t285 = cos(qJ(4));
t230 = -t235 * t282 + t285 * t267;
t295 = qJD(3) * qJD(4);
t294 = t285 * t295;
t265 = qJDD(3) * t282 + t294;
t227 = (-t265 + t294) * pkin(7) + (t282 * t285 * t287 + qJDD(4)) * pkin(4) + t230;
t231 = t285 * t235 + t282 * t267;
t266 = qJDD(3) * t285 - t282 * t295;
t296 = t282 * qJD(3);
t271 = qJD(4) * pkin(4) - pkin(7) * t296;
t275 = t285 ^ 2;
t228 = -pkin(4) * t275 * t287 + pkin(7) * t266 - qJD(4) * t271 + t231;
t281 = sin(qJ(5));
t284 = cos(qJ(5));
t225 = t227 * t284 - t228 * t281;
t258 = (-t282 * t281 + t285 * t284) * qJD(3);
t242 = qJD(5) * t258 + t265 * t284 + t266 * t281;
t259 = (t285 * t281 + t282 * t284) * qJD(3);
t247 = -mrSges(6,1) * t258 + mrSges(6,2) * t259;
t274 = qJD(4) + qJD(5);
t251 = -mrSges(6,2) * t274 + mrSges(6,3) * t258;
t273 = qJDD(4) + qJDD(5);
t222 = m(6) * t225 + mrSges(6,1) * t273 - mrSges(6,3) * t242 - t247 * t259 + t251 * t274;
t226 = t227 * t281 + t228 * t284;
t241 = -qJD(5) * t259 - t265 * t281 + t266 * t284;
t252 = mrSges(6,1) * t274 - mrSges(6,3) * t259;
t223 = m(6) * t226 - mrSges(6,2) * t273 + mrSges(6,3) * t241 + t247 * t258 - t252 * t274;
t215 = t284 * t222 + t281 * t223;
t297 = qJD(3) * t285;
t264 = (-t285 * mrSges(5,1) + t282 * mrSges(5,2)) * qJD(3);
t270 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t297;
t213 = m(5) * t230 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t265 + qJD(4) * t270 - t264 * t296 + t215;
t269 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t296;
t292 = -t222 * t281 + t284 * t223;
t214 = m(5) * t231 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t266 - qJD(4) * t269 + t264 * t297 + t292;
t293 = -t213 * t282 + t285 * t214;
t239 = t253 * t286 - t283 * t254;
t291 = -qJDD(3) * pkin(3) - t239;
t229 = t271 * t296 - pkin(4) * t266 + (-pkin(7) * t275 - pkin(6)) * t287 + t291;
t290 = m(6) * t229 - t241 * mrSges(6,1) + mrSges(6,2) * t242 - t258 * t251 + t252 * t259;
t244 = Ifges(6,4) * t259 + Ifges(6,2) * t258 + Ifges(6,6) * t274;
t245 = Ifges(6,1) * t259 + Ifges(6,4) * t258 + Ifges(6,5) * t274;
t289 = mrSges(6,1) * t225 - mrSges(6,2) * t226 + Ifges(6,5) * t242 + Ifges(6,6) * t241 + Ifges(6,3) * t273 + t259 * t244 - t245 * t258;
t234 = -pkin(6) * t287 + t291;
t288 = -m(5) * t234 + t266 * mrSges(5,1) - mrSges(5,2) * t265 - t269 * t296 + t270 * t297 - t290;
t257 = Ifges(5,5) * qJD(4) + (t282 * Ifges(5,1) + t285 * Ifges(5,4)) * qJD(3);
t256 = Ifges(5,6) * qJD(4) + (t282 * Ifges(5,4) + t285 * Ifges(5,2)) * qJD(3);
t243 = Ifges(6,5) * t259 + Ifges(6,6) * t258 + Ifges(6,3) * t274;
t218 = m(4) * t239 + qJDD(3) * mrSges(4,1) - mrSges(4,2) * t287 + t288;
t217 = mrSges(6,2) * t229 - mrSges(6,3) * t225 + Ifges(6,1) * t242 + Ifges(6,4) * t241 + Ifges(6,5) * t273 + t243 * t258 - t244 * t274;
t216 = -mrSges(6,1) * t229 + mrSges(6,3) * t226 + Ifges(6,4) * t242 + Ifges(6,2) * t241 + Ifges(6,6) * t273 - t243 * t259 + t245 * t274;
t211 = m(4) * t240 - mrSges(4,1) * t287 - qJDD(3) * mrSges(4,2) + t293;
t1 = [m(2) * t276 + t277 * (m(3) * t254 + t211 * t286 - t218 * t283) + t279 * (m(3) * t253 + t211 * t283 + t218 * t286); t213 * t285 + t214 * t282 + (m(3) + m(4)) * t267; Ifges(4,3) * qJDD(3) + mrSges(4,1) * t239 - mrSges(4,2) * t240 + t282 * (mrSges(5,2) * t234 - mrSges(5,3) * t230 + Ifges(5,1) * t265 + Ifges(5,4) * t266 + Ifges(5,5) * qJDD(4) - pkin(7) * t215 - qJD(4) * t256 - t216 * t281 + t217 * t284) + t285 * (-mrSges(5,1) * t234 + mrSges(5,3) * t231 + Ifges(5,4) * t265 + Ifges(5,2) * t266 + Ifges(5,6) * qJDD(4) - pkin(4) * t290 + pkin(7) * t292 + qJD(4) * t257 + t284 * t216 + t281 * t217) + pkin(3) * t288 + pkin(6) * t293; mrSges(5,1) * t230 - mrSges(5,2) * t231 + Ifges(5,5) * t265 + Ifges(5,6) * t266 + Ifges(5,3) * qJDD(4) + pkin(4) * t215 + (t282 * t256 - t285 * t257) * qJD(3) + t289; t289;];
tauJ = t1;
