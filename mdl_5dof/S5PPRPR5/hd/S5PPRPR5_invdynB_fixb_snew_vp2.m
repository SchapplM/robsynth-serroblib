% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PPRPR5
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
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
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PPRPR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR5_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR5_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR5_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR5_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR5_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:27
% EndTime: 2019-12-31 17:33:27
% DurationCPUTime: 0.44s
% Computational Cost: add. (3091->151), mult. (4873->179), div. (0->0), fcn. (2426->6), ass. (0->65)
t293 = -m(4) - m(5);
t292 = -pkin(3) - pkin(6);
t291 = mrSges(4,1) - mrSges(5,2);
t290 = -mrSges(2,2) + mrSges(3,3);
t289 = -Ifges(5,4) + Ifges(4,5);
t288 = Ifges(5,5) - Ifges(4,6);
t266 = sin(pkin(7));
t267 = cos(pkin(7));
t257 = t266 * g(1) - t267 * g(2);
t255 = qJDD(2) - t257;
t258 = -t267 * g(1) - t266 * g(2);
t269 = sin(qJ(3));
t271 = cos(qJ(3));
t241 = t271 * t255 - t269 * t258;
t272 = qJD(3) ^ 2;
t276 = -t272 * qJ(4) + qJDD(4) - t241;
t238 = t292 * qJDD(3) + t276;
t263 = g(3) - qJDD(1);
t268 = sin(qJ(5));
t270 = cos(qJ(5));
t234 = t270 * t238 - t268 * t263;
t252 = (mrSges(6,1) * t268 + mrSges(6,2) * t270) * qJD(3);
t283 = qJD(3) * qJD(5);
t254 = t270 * qJDD(3) - t268 * t283;
t285 = qJD(3) * t268;
t259 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t285;
t284 = qJD(3) * t270;
t232 = m(6) * t234 + qJDD(5) * mrSges(6,1) - t254 * mrSges(6,3) + qJD(5) * t259 - t252 * t284;
t235 = t268 * t238 + t270 * t263;
t253 = -t268 * qJDD(3) - t270 * t283;
t260 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t284;
t233 = m(6) * t235 - qJDD(5) * mrSges(6,2) + t253 * mrSges(6,3) - qJD(5) * t260 - t252 * t285;
t225 = t270 * t232 + t268 * t233;
t240 = -qJDD(3) * pkin(3) + t276;
t274 = -m(5) * t240 + t272 * mrSges(5,3) - t225;
t222 = m(4) * t241 - t272 * mrSges(4,2) + t291 * qJDD(3) + t274;
t242 = t269 * t255 + t271 * t258;
t275 = qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) + t242;
t239 = t272 * pkin(3) - t275;
t237 = t292 * t272 + t275;
t278 = -m(6) * t237 + t253 * mrSges(6,1) - t254 * mrSges(6,2) - t259 * t285 - t260 * t284;
t273 = -m(5) * t239 + t272 * mrSges(5,2) + qJDD(3) * mrSges(5,3) - t278;
t227 = m(4) * t242 - t272 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t273;
t221 = t271 * t222 + t269 * t227;
t277 = -m(3) * t255 - t221;
t219 = m(2) * t257 + t277;
t280 = -t269 * t222 + t271 * t227;
t279 = m(3) * t258 + t280;
t220 = m(2) * t258 + t279;
t287 = t267 * t219 + t266 * t220;
t286 = -t268 * t232 + t270 * t233;
t282 = -m(3) * t263 - t286;
t281 = -t266 * t219 + t267 * t220;
t245 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t270 - Ifges(6,4) * t268) * qJD(3);
t244 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t270 - Ifges(6,2) * t268) * qJD(3);
t243 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t270 - Ifges(6,6) * t268) * qJD(3);
t229 = mrSges(6,2) * t237 - mrSges(6,3) * t234 + Ifges(6,1) * t254 + Ifges(6,4) * t253 + Ifges(6,5) * qJDD(5) - qJD(5) * t244 - t243 * t285;
t228 = -mrSges(6,1) * t237 + mrSges(6,3) * t235 + Ifges(6,4) * t254 + Ifges(6,2) * t253 + Ifges(6,6) * qJDD(5) + qJD(5) * t245 - t243 * t284;
t224 = m(5) * t263 + t286;
t223 = t293 * t263 + t282;
t215 = mrSges(5,1) * t240 + mrSges(6,1) * t234 - mrSges(6,2) * t235 - mrSges(4,3) * t241 + Ifges(6,5) * t254 + Ifges(6,6) * t253 + Ifges(6,3) * qJDD(5) + pkin(4) * t225 - qJ(4) * t224 + t288 * t272 + (mrSges(4,2) - mrSges(5,3)) * t263 + t289 * qJDD(3) + (t270 * t244 + t268 * t245) * qJD(3);
t214 = -mrSges(5,1) * t239 + mrSges(4,3) * t242 - pkin(3) * t224 - pkin(4) * t278 - pkin(6) * t286 - t288 * qJDD(3) - t270 * t228 - t268 * t229 - t291 * t263 + t289 * t272;
t213 = mrSges(3,2) * t255 - mrSges(2,3) * t257 - pkin(5) * t221 - qJ(2) * t223 - t269 * t214 + t271 * t215 + t290 * t263;
t212 = -t269 * t215 - t271 * t214 + pkin(2) * t286 - pkin(5) * t280 - pkin(1) * t223 + (mrSges(3,2) + mrSges(2,3)) * t258 + (-pkin(2) * t293 + mrSges(2,1) + mrSges(3,1)) * t263;
t1 = [-m(1) * g(1) + t281; -m(1) * g(2) + t287; -m(1) * g(3) + (-m(2) + t293) * t263 + t282; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t287 - t266 * t212 + t267 * t213; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t281 + t267 * t212 + t266 * t213; pkin(1) * t277 + qJ(2) * t279 + mrSges(2,1) * t257 - pkin(2) * t221 - mrSges(3,1) * t255 - pkin(3) * t274 - qJ(4) * t273 + pkin(6) * t225 - mrSges(4,1) * t241 + mrSges(4,2) * t242 - mrSges(5,2) * t240 + mrSges(5,3) * t239 - t270 * t229 + t268 * t228 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t290 * t258 + (pkin(3) * mrSges(5,2) - Ifges(5,1) - Ifges(4,3)) * qJDD(3);];
tauB = t1;
