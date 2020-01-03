% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RPPR4
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
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RPPR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR4_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR4_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR4_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR4_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR4_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:51
% EndTime: 2019-12-31 16:38:51
% DurationCPUTime: 0.41s
% Computational Cost: add. (2907->146), mult. (4891->178), div. (0->0), fcn. (2134->6), ass. (0->61)
t279 = -pkin(2) - pkin(5);
t278 = mrSges(3,1) - mrSges(4,2);
t277 = -Ifges(4,4) + Ifges(3,5);
t276 = Ifges(4,5) - Ifges(3,6);
t259 = sin(qJ(1));
t261 = cos(qJ(1));
t247 = t259 * g(1) - g(2) * t261;
t240 = qJDD(1) * pkin(1) + t247;
t248 = -g(1) * t261 - g(2) * t259;
t262 = qJD(1) ^ 2;
t242 = -pkin(1) * t262 + t248;
t256 = sin(pkin(6));
t257 = cos(pkin(6));
t229 = t240 * t257 - t256 * t242;
t266 = -qJ(3) * t262 + qJDD(3) - t229;
t226 = t279 * qJDD(1) + t266;
t253 = -g(3) + qJDD(2);
t258 = sin(qJ(4));
t260 = cos(qJ(4));
t222 = t226 * t260 - t253 * t258;
t241 = (mrSges(5,1) * t258 + mrSges(5,2) * t260) * qJD(1);
t272 = qJD(1) * qJD(4);
t244 = qJDD(1) * t260 - t258 * t272;
t274 = qJD(1) * t258;
t245 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t274;
t273 = qJD(1) * t260;
t220 = m(5) * t222 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t244 + qJD(4) * t245 - t241 * t273;
t223 = t226 * t258 + t253 * t260;
t243 = -qJDD(1) * t258 - t260 * t272;
t246 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t273;
t221 = m(5) * t223 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t243 - qJD(4) * t246 - t241 * t274;
t213 = t220 * t260 + t221 * t258;
t228 = -qJDD(1) * pkin(2) + t266;
t264 = -m(4) * t228 + t262 * mrSges(4,3) - t213;
t211 = m(3) * t229 - mrSges(3,2) * t262 + t278 * qJDD(1) + t264;
t230 = t256 * t240 + t257 * t242;
t265 = qJDD(1) * qJ(3) + 0.2e1 * qJD(3) * qJD(1) + t230;
t227 = pkin(2) * t262 - t265;
t225 = t279 * t262 + t265;
t267 = -m(5) * t225 + mrSges(5,1) * t243 - t244 * mrSges(5,2) - t245 * t274 - t246 * t273;
t263 = -m(4) * t227 + t262 * mrSges(4,2) + qJDD(1) * mrSges(4,3) - t267;
t216 = m(3) * t230 - mrSges(3,1) * t262 - qJDD(1) * mrSges(3,2) + t263;
t209 = t257 * t211 + t256 * t216;
t207 = m(2) * t247 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t262 + t209;
t270 = -t256 * t211 + t257 * t216;
t208 = m(2) * t248 - mrSges(2,1) * t262 - qJDD(1) * mrSges(2,2) + t270;
t275 = t261 * t207 + t259 * t208;
t271 = -t207 * t259 + t261 * t208;
t269 = -t220 * t258 + t260 * t221;
t212 = m(4) * t253 + t269;
t268 = m(3) * t253 + t212;
t236 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t260 - Ifges(5,4) * t258) * qJD(1);
t235 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t260 - Ifges(5,2) * t258) * qJD(1);
t234 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t260 - Ifges(5,6) * t258) * qJD(1);
t218 = mrSges(5,2) * t225 - mrSges(5,3) * t222 + Ifges(5,1) * t244 + Ifges(5,4) * t243 + Ifges(5,5) * qJDD(4) - qJD(4) * t235 - t234 * t274;
t217 = -mrSges(5,1) * t225 + mrSges(5,3) * t223 + Ifges(5,4) * t244 + Ifges(5,2) * t243 + Ifges(5,6) * qJDD(4) + qJD(4) * t236 - t234 * t273;
t203 = mrSges(4,1) * t228 + mrSges(5,1) * t222 - mrSges(5,2) * t223 - mrSges(3,3) * t229 + Ifges(5,5) * t244 + Ifges(5,6) * t243 + Ifges(5,3) * qJDD(4) + pkin(3) * t213 - qJ(3) * t212 + t276 * t262 + (mrSges(3,2) - mrSges(4,3)) * t253 + t277 * qJDD(1) + (t260 * t235 + t258 * t236) * qJD(1);
t202 = -mrSges(4,1) * t227 + mrSges(3,3) * t230 - pkin(2) * t212 - pkin(3) * t267 - pkin(5) * t269 - t276 * qJDD(1) - t260 * t217 - t258 * t218 - t278 * t253 + t277 * t262;
t201 = -mrSges(2,2) * g(3) - mrSges(2,3) * t247 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t262 - qJ(2) * t209 - t202 * t256 + t203 * t257;
t200 = mrSges(2,1) * g(3) + mrSges(2,3) * t248 + t262 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t268 + qJ(2) * t270 + t257 * t202 + t256 * t203;
t1 = [-m(1) * g(1) + t271; -m(1) * g(2) + t275; (-m(1) - m(2)) * g(3) + t268; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t275 - t259 * t200 + t261 * t201; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t271 + t261 * t200 + t259 * t201; pkin(1) * t209 + mrSges(2,1) * t247 - mrSges(2,2) * t248 + pkin(2) * t264 + qJ(3) * t263 - t258 * t217 - pkin(5) * t213 + mrSges(3,1) * t229 - mrSges(3,2) * t230 + mrSges(4,2) * t228 - mrSges(4,3) * t227 + t260 * t218 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(2) * mrSges(4,2) + Ifges(4,1) + Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB = t1;
