% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RRPR5
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RRPR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR5_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR5_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR5_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR5_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR5_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:30
% EndTime: 2019-12-31 17:03:31
% DurationCPUTime: 0.43s
% Computational Cost: add. (4080->149), mult. (4891->180), div. (0->0), fcn. (2134->6), ass. (0->63)
t291 = -m(3) - m(4);
t290 = -pkin(2) - pkin(6);
t289 = mrSges(3,1) - mrSges(4,2);
t288 = -Ifges(4,4) + Ifges(3,5);
t287 = Ifges(4,5) - Ifges(3,6);
t267 = qJD(1) + qJD(2);
t268 = sin(qJ(4));
t286 = t267 * t268;
t271 = cos(qJ(4));
t285 = t267 * t271;
t270 = sin(qJ(1));
t273 = cos(qJ(1));
t259 = t270 * g(1) - t273 * g(2);
t255 = qJDD(1) * pkin(1) + t259;
t260 = -t273 * g(1) - t270 * g(2);
t274 = qJD(1) ^ 2;
t256 = -t274 * pkin(1) + t260;
t269 = sin(qJ(2));
t272 = cos(qJ(2));
t241 = t272 * t255 - t269 * t256;
t265 = t267 ^ 2;
t266 = qJDD(1) + qJDD(2);
t277 = -t265 * qJ(3) + qJDD(3) - t241;
t238 = t290 * t266 + t277;
t234 = t268 * g(3) + t271 * t238;
t249 = (mrSges(5,1) * t268 + mrSges(5,2) * t271) * t267;
t283 = qJD(4) * t267;
t251 = t271 * t266 - t268 * t283;
t257 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t286;
t232 = m(5) * t234 + qJDD(4) * mrSges(5,1) - t251 * mrSges(5,3) + qJD(4) * t257 - t249 * t285;
t235 = -t271 * g(3) + t268 * t238;
t250 = -t268 * t266 - t271 * t283;
t258 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t285;
t233 = m(5) * t235 - qJDD(4) * mrSges(5,2) + t250 * mrSges(5,3) - qJD(4) * t258 - t249 * t286;
t225 = t271 * t232 + t268 * t233;
t240 = -t266 * pkin(2) + t277;
t276 = -m(4) * t240 + t265 * mrSges(4,3) - t225;
t223 = m(3) * t241 - t265 * mrSges(3,2) + t289 * t266 + t276;
t242 = t269 * t255 + t272 * t256;
t278 = t266 * qJ(3) + 0.2e1 * qJD(3) * t267 + t242;
t239 = t265 * pkin(2) - t278;
t237 = t290 * t265 + t278;
t279 = -m(5) * t237 + t250 * mrSges(5,1) - t251 * mrSges(5,2) - t257 * t286 - t258 * t285;
t275 = -m(4) * t239 + t265 * mrSges(4,2) + t266 * mrSges(4,3) - t279;
t228 = m(3) * t242 - t265 * mrSges(3,1) - t266 * mrSges(3,2) + t275;
t221 = t272 * t223 + t269 * t228;
t219 = m(2) * t259 + qJDD(1) * mrSges(2,1) - t274 * mrSges(2,2) + t221;
t281 = -t269 * t223 + t272 * t228;
t220 = m(2) * t260 - t274 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t281;
t284 = t273 * t219 + t270 * t220;
t282 = -t270 * t219 + t273 * t220;
t280 = t268 * t232 - t271 * t233;
t245 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t271 - Ifges(5,4) * t268) * t267;
t244 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t271 - Ifges(5,2) * t268) * t267;
t243 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t271 - Ifges(5,6) * t268) * t267;
t230 = mrSges(5,2) * t237 - mrSges(5,3) * t234 + Ifges(5,1) * t251 + Ifges(5,4) * t250 + Ifges(5,5) * qJDD(4) - qJD(4) * t244 - t243 * t286;
t229 = -mrSges(5,1) * t237 + mrSges(5,3) * t235 + Ifges(5,4) * t251 + Ifges(5,2) * t250 + Ifges(5,6) * qJDD(4) + qJD(4) * t245 - t243 * t285;
t224 = -m(4) * g(3) - t280;
t215 = mrSges(4,1) * t240 + mrSges(5,1) * t234 - mrSges(5,2) * t235 - mrSges(3,3) * t241 + Ifges(5,5) * t251 + Ifges(5,6) * t250 + Ifges(5,3) * qJDD(4) + pkin(3) * t225 - qJ(3) * t224 + (t271 * t244 + t268 * t245) * t267 + t288 * t266 + t287 * t265 + (-mrSges(3,2) + mrSges(4,3)) * g(3);
t214 = -mrSges(4,1) * t239 + mrSges(3,3) * t242 - pkin(2) * t224 - pkin(3) * t279 + pkin(6) * t280 + t289 * g(3) - t271 * t229 - t268 * t230 + t288 * t265 - t287 * t266;
t213 = -mrSges(2,2) * g(3) - mrSges(2,3) * t259 + Ifges(2,5) * qJDD(1) - t274 * Ifges(2,6) - pkin(5) * t221 - t269 * t214 + t272 * t215;
t212 = Ifges(2,6) * qJDD(1) + t274 * Ifges(2,5) + mrSges(2,3) * t260 + t269 * t215 + t272 * t214 + pkin(1) * t280 + pkin(5) * t281 + (-pkin(1) * t291 + mrSges(2,1)) * g(3);
t1 = [-m(1) * g(1) + t282; -m(1) * g(2) + t284; (-m(1) - m(2) + t291) * g(3) - t280; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t284 - t270 * t212 + t273 * t213; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t282 + t273 * t212 + t270 * t213; pkin(1) * t221 + mrSges(2,1) * t259 - mrSges(2,2) * t260 + qJ(3) * t275 + pkin(2) * t276 + t271 * t230 - t268 * t229 - pkin(6) * t225 + mrSges(3,1) * t241 - mrSges(3,2) * t242 + mrSges(4,2) * t240 - mrSges(4,3) * t239 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + (-pkin(2) * mrSges(4,2) + Ifges(4,1) + Ifges(3,3)) * t266;];
tauB = t1;
