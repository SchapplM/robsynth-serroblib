% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4PRRR5
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
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4PRRR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR5_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR5_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR5_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR5_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:38
% EndTime: 2019-12-31 16:33:39
% DurationCPUTime: 0.56s
% Computational Cost: add. (6191->144), mult. (7850->187), div. (0->0), fcn. (4448->8), ass. (0->63)
t291 = m(3) + m(4);
t290 = cos(pkin(7));
t271 = qJD(2) + qJD(3);
t274 = sin(qJ(4));
t289 = t271 * t274;
t277 = cos(qJ(4));
t288 = t271 * t277;
t273 = sin(pkin(7));
t265 = -t290 * g(1) - t273 * g(2);
t272 = -g(3) + qJDD(1);
t276 = sin(qJ(2));
t279 = cos(qJ(2));
t251 = -t276 * t265 + t279 * t272;
t249 = qJDD(2) * pkin(2) + t251;
t252 = t279 * t265 + t276 * t272;
t280 = qJD(2) ^ 2;
t250 = -t280 * pkin(2) + t252;
t275 = sin(qJ(3));
t278 = cos(qJ(3));
t246 = t275 * t249 + t278 * t250;
t269 = t271 ^ 2;
t270 = qJDD(2) + qJDD(3);
t244 = -t269 * pkin(3) + t270 * pkin(6) + t246;
t264 = t273 * g(1) - t290 * g(2);
t241 = -t274 * t244 - t277 * t264;
t258 = (-mrSges(5,1) * t277 + mrSges(5,2) * t274) * t271;
t286 = qJD(4) * t271;
t259 = t274 * t270 + t277 * t286;
t263 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t288;
t239 = m(5) * t241 + qJDD(4) * mrSges(5,1) - t259 * mrSges(5,3) + qJD(4) * t263 - t258 * t289;
t242 = t277 * t244 - t274 * t264;
t260 = t277 * t270 - t274 * t286;
t262 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t289;
t240 = m(5) * t242 - qJDD(4) * mrSges(5,2) + t260 * mrSges(5,3) - qJD(4) * t262 + t258 * t288;
t282 = -t274 * t239 + t277 * t240;
t228 = m(4) * t246 - t269 * mrSges(4,1) - t270 * mrSges(4,2) + t282;
t245 = t278 * t249 - t275 * t250;
t243 = -t270 * pkin(3) - t269 * pkin(6) - t245;
t281 = -m(5) * t243 + t260 * mrSges(5,1) - t259 * mrSges(5,2) - t262 * t289 + t263 * t288;
t235 = m(4) * t245 + t270 * mrSges(4,1) - t269 * mrSges(4,2) + t281;
t225 = t275 * t228 + t278 * t235;
t223 = m(3) * t251 + qJDD(2) * mrSges(3,1) - t280 * mrSges(3,2) + t225;
t283 = t278 * t228 - t275 * t235;
t224 = m(3) * t252 - t280 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t283;
t284 = -t276 * t223 + t279 * t224;
t216 = m(2) * t265 + t284;
t231 = t277 * t239 + t274 * t240;
t230 = (m(2) + t291) * t264 - t231;
t287 = t273 * t216 + t290 * t230;
t217 = t279 * t223 + t276 * t224;
t285 = t290 * t216 - t273 * t230;
t255 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t274 + Ifges(5,4) * t277) * t271;
t254 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t274 + Ifges(5,2) * t277) * t271;
t253 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t274 + Ifges(5,6) * t277) * t271;
t233 = mrSges(5,2) * t243 - mrSges(5,3) * t241 + Ifges(5,1) * t259 + Ifges(5,4) * t260 + Ifges(5,5) * qJDD(4) - qJD(4) * t254 + t253 * t288;
t232 = -mrSges(5,1) * t243 + mrSges(5,3) * t242 + Ifges(5,4) * t259 + Ifges(5,2) * t260 + Ifges(5,6) * qJDD(4) + qJD(4) * t255 - t253 * t289;
t219 = mrSges(4,1) * t264 - mrSges(5,1) * t241 + mrSges(5,2) * t242 + mrSges(4,3) * t246 + t269 * Ifges(4,5) - Ifges(5,5) * t259 + Ifges(4,6) * t270 - Ifges(5,6) * t260 - Ifges(5,3) * qJDD(4) - pkin(3) * t231 + (-t254 * t274 + t255 * t277) * t271;
t218 = -mrSges(4,2) * t264 - mrSges(4,3) * t245 + Ifges(4,5) * t270 - t269 * Ifges(4,6) - pkin(6) * t231 - t274 * t232 + t277 * t233;
t213 = -mrSges(3,2) * t264 - mrSges(3,3) * t251 + Ifges(3,5) * qJDD(2) - t280 * Ifges(3,6) - pkin(5) * t225 + t278 * t218 - t275 * t219;
t212 = Ifges(3,6) * qJDD(2) + t280 * Ifges(3,5) + mrSges(3,1) * t264 + mrSges(3,3) * t252 + t275 * t218 + t278 * t219 - pkin(2) * (-m(4) * t264 + t231) + pkin(5) * t283;
t211 = -mrSges(2,1) * t272 - mrSges(3,1) * t251 - mrSges(4,1) * t245 + mrSges(3,2) * t252 + mrSges(4,2) * t246 + mrSges(2,3) * t265 - Ifges(3,3) * qJDD(2) - Ifges(4,3) * t270 - pkin(1) * t217 - pkin(2) * t225 - pkin(3) * t281 - pkin(6) * t282 - t277 * t232 - t274 * t233;
t210 = mrSges(2,2) * t272 - mrSges(2,3) * t264 - pkin(4) * t217 - t276 * t212 + t279 * t213;
t1 = [-m(1) * g(1) + t285; -m(1) * g(2) + t287; -m(1) * g(3) + m(2) * t272 + t217; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t287 + t290 * t210 - t273 * t211; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t285 + t273 * t210 + t290 * t211; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - mrSges(2,2) * t265 + t276 * t213 + t279 * t212 - pkin(1) * t231 + pkin(4) * t284 + (pkin(1) * t291 + mrSges(2,1)) * t264;];
tauB = t1;
