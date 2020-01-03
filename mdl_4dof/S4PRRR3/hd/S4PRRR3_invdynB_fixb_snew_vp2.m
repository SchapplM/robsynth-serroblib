% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4PRRR3
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4PRRR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR3_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR3_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR3_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR3_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR3_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:37
% EndTime: 2019-12-31 16:31:37
% DurationCPUTime: 0.61s
% Computational Cost: add. (7374->142), mult. (10224->187), div. (0->0), fcn. (6170->8), ass. (0->64)
t266 = qJD(2) + qJD(3);
t270 = sin(qJ(4));
t287 = t266 * t270;
t273 = cos(qJ(4));
t286 = t266 * t273;
t268 = sin(pkin(7));
t269 = cos(pkin(7));
t260 = t268 * g(1) - t269 * g(2);
t261 = -t269 * g(1) - t268 * g(2);
t272 = sin(qJ(2));
t275 = cos(qJ(2));
t245 = t275 * t260 - t272 * t261;
t243 = qJDD(2) * pkin(2) + t245;
t246 = t272 * t260 + t275 * t261;
t276 = qJD(2) ^ 2;
t244 = -t276 * pkin(2) + t246;
t271 = sin(qJ(3));
t274 = cos(qJ(3));
t240 = t271 * t243 + t274 * t244;
t264 = t266 ^ 2;
t265 = qJDD(2) + qJDD(3);
t238 = -t264 * pkin(3) + t265 * pkin(6) + t240;
t267 = -g(3) + qJDD(1);
t235 = -t270 * t238 + t273 * t267;
t252 = (-mrSges(5,1) * t273 + mrSges(5,2) * t270) * t266;
t284 = qJD(4) * t266;
t253 = t270 * t265 + t273 * t284;
t259 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t286;
t233 = m(5) * t235 + qJDD(4) * mrSges(5,1) - t253 * mrSges(5,3) + qJD(4) * t259 - t252 * t287;
t236 = t273 * t238 + t270 * t267;
t254 = t273 * t265 - t270 * t284;
t258 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t287;
t234 = m(5) * t236 - qJDD(4) * mrSges(5,2) + t254 * mrSges(5,3) - qJD(4) * t258 + t252 * t286;
t279 = -t270 * t233 + t273 * t234;
t224 = m(4) * t240 - t264 * mrSges(4,1) - t265 * mrSges(4,2) + t279;
t239 = t274 * t243 - t271 * t244;
t237 = -t265 * pkin(3) - t264 * pkin(6) - t239;
t277 = -m(5) * t237 + t254 * mrSges(5,1) - t253 * mrSges(5,2) - t258 * t287 + t259 * t286;
t229 = m(4) * t239 + t265 * mrSges(4,1) - t264 * mrSges(4,2) + t277;
t221 = t271 * t224 + t274 * t229;
t219 = m(3) * t245 + qJDD(2) * mrSges(3,1) - t276 * mrSges(3,2) + t221;
t280 = t274 * t224 - t271 * t229;
t220 = m(3) * t246 - t276 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t280;
t213 = t275 * t219 + t272 * t220;
t211 = m(2) * t260 + t213;
t281 = -t272 * t219 + t275 * t220;
t212 = m(2) * t261 + t281;
t285 = t269 * t211 + t268 * t212;
t225 = t273 * t233 + t270 * t234;
t283 = m(4) * t267 + t225;
t282 = -t268 * t211 + t269 * t212;
t278 = m(3) * t267 + t283;
t249 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t270 + Ifges(5,4) * t273) * t266;
t248 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t270 + Ifges(5,2) * t273) * t266;
t247 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t270 + Ifges(5,6) * t273) * t266;
t227 = mrSges(5,2) * t237 - mrSges(5,3) * t235 + Ifges(5,1) * t253 + Ifges(5,4) * t254 + Ifges(5,5) * qJDD(4) - qJD(4) * t248 + t247 * t286;
t226 = -mrSges(5,1) * t237 + mrSges(5,3) * t236 + Ifges(5,4) * t253 + Ifges(5,2) * t254 + Ifges(5,6) * qJDD(4) + qJD(4) * t249 - t247 * t287;
t215 = -mrSges(4,1) * t267 - mrSges(5,1) * t235 + mrSges(5,2) * t236 + mrSges(4,3) * t240 + t264 * Ifges(4,5) - Ifges(5,5) * t253 + Ifges(4,6) * t265 - Ifges(5,6) * t254 - Ifges(5,3) * qJDD(4) - pkin(3) * t225 + (-t248 * t270 + t249 * t273) * t266;
t214 = mrSges(4,2) * t267 - mrSges(4,3) * t239 + Ifges(4,5) * t265 - t264 * Ifges(4,6) - pkin(6) * t225 - t270 * t226 + t273 * t227;
t207 = mrSges(3,2) * t267 - mrSges(3,3) * t245 + Ifges(3,5) * qJDD(2) - t276 * Ifges(3,6) - pkin(5) * t221 + t274 * t214 - t271 * t215;
t206 = -mrSges(3,1) * t267 + mrSges(3,3) * t246 + t276 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t283 + pkin(5) * t280 + t271 * t214 + t274 * t215;
t205 = mrSges(2,2) * t267 - mrSges(2,3) * t260 - pkin(4) * t213 - t272 * t206 + t275 * t207;
t204 = -mrSges(2,1) * t267 + mrSges(2,3) * t261 - pkin(1) * t278 + pkin(4) * t281 + t275 * t206 + t272 * t207;
t1 = [-m(1) * g(1) + t282; -m(1) * g(2) + t285; -m(1) * g(3) + m(2) * t267 + t278; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t285 - t268 * t204 + t269 * t205; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t282 + t269 * t204 + t268 * t205; pkin(1) * t213 + mrSges(2,1) * t260 - mrSges(2,2) * t261 + pkin(2) * t221 + mrSges(3,1) * t245 - mrSges(3,2) * t246 + t270 * t227 + t273 * t226 + pkin(3) * t277 + pkin(6) * t279 + mrSges(4,1) * t239 - mrSges(4,2) * t240 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(4,3) * t265 + Ifges(3,3) * qJDD(2);];
tauB = t1;
