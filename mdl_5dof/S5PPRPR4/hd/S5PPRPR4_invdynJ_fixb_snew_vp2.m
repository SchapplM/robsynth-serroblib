% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PPRPR4
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PPRPR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR4_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR4_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR4_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR4_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR4_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:20
% EndTime: 2019-12-31 17:32:21
% DurationCPUTime: 0.32s
% Computational Cost: add. (1176->115), mult. (2401->153), div. (0->0), fcn. (1522->8), ass. (0->59)
t256 = qJD(3) ^ 2;
t250 = cos(pkin(8));
t246 = t250 ^ 2;
t272 = 0.2e1 * t250;
t271 = pkin(4) * t250;
t248 = sin(pkin(8));
t270 = mrSges(5,2) * t248;
t269 = t246 * t256;
t249 = sin(pkin(7));
t251 = cos(pkin(7));
t236 = -g(1) * t249 + g(2) * t251 + qJDD(2);
t237 = -g(1) * t251 - g(2) * t249;
t253 = sin(qJ(3));
t255 = cos(qJ(3));
t226 = t253 * t236 + t255 * t237;
t222 = -pkin(3) * t256 + qJDD(3) * qJ(4) + t226;
t247 = g(3) - qJDD(1);
t266 = qJD(3) * qJD(4);
t267 = t250 * t247 - 0.2e1 * t248 * t266;
t209 = (-pkin(6) * qJDD(3) + t256 * t271 - t222) * t248 + t267;
t213 = t250 * t222 + t248 * t247 + t266 * t272;
t265 = qJDD(3) * t250;
t210 = -pkin(4) * t269 + pkin(6) * t265 + t213;
t252 = sin(qJ(5));
t254 = cos(qJ(5));
t207 = t209 * t254 - t210 * t252;
t260 = -t248 * t252 + t250 * t254;
t229 = t260 * qJD(3);
t261 = t248 * t254 + t250 * t252;
t230 = t261 * qJD(3);
t220 = -mrSges(6,1) * t229 + mrSges(6,2) * t230;
t224 = t229 * qJD(5) + t261 * qJDD(3);
t227 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t229;
t205 = m(6) * t207 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t224 + qJD(5) * t227 - t220 * t230;
t208 = t209 * t252 + t210 * t254;
t223 = -t230 * qJD(5) + t260 * qJDD(3);
t228 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t230;
t206 = m(6) * t208 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t223 - qJD(5) * t228 + t220 * t229;
t268 = t254 * t205 + t252 * t206;
t212 = -t222 * t248 + t267;
t259 = mrSges(5,3) * qJDD(3) + t256 * (-mrSges(5,1) * t250 + t270);
t197 = m(5) * t212 - t259 * t248 + t268;
t263 = -t205 * t252 + t254 * t206;
t198 = m(5) * t213 + t259 * t250 + t263;
t264 = -t197 * t248 + t250 * t198;
t225 = t236 * t255 - t253 * t237;
t262 = qJDD(4) - t225;
t245 = t248 ^ 2;
t211 = (-pkin(3) - t271) * qJDD(3) + (-qJ(4) + (-t245 - t246) * pkin(6)) * t256 + t262;
t258 = m(6) * t211 - t223 * mrSges(6,1) + mrSges(6,2) * t224 - t229 * t227 + t228 * t230;
t219 = -qJDD(3) * pkin(3) - qJ(4) * t256 + t262;
t257 = -m(5) * t219 + mrSges(5,1) * t265 - t258 + (t245 * t256 + t269) * mrSges(5,3);
t216 = Ifges(6,1) * t230 + Ifges(6,4) * t229 + Ifges(6,5) * qJD(5);
t215 = Ifges(6,4) * t230 + Ifges(6,2) * t229 + Ifges(6,6) * qJD(5);
t214 = Ifges(6,5) * t230 + Ifges(6,6) * t229 + Ifges(6,3) * qJD(5);
t201 = qJDD(3) * t270 - t257;
t200 = mrSges(6,2) * t211 - mrSges(6,3) * t207 + Ifges(6,1) * t224 + Ifges(6,4) * t223 + Ifges(6,5) * qJDD(5) - qJD(5) * t215 + t214 * t229;
t199 = -mrSges(6,1) * t211 + mrSges(6,3) * t208 + Ifges(6,4) * t224 + Ifges(6,2) * t223 + Ifges(6,6) * qJDD(5) + qJD(5) * t216 - t214 * t230;
t1 = [-t250 * t197 - t248 * t198 + (-m(2) - m(3) - m(4)) * t247; m(3) * t236 + t253 * (m(4) * t226 - mrSges(4,1) * t256 - qJDD(3) * mrSges(4,2) + t264) + t255 * (m(4) * t225 - mrSges(4,2) * t256 + (mrSges(4,1) - t270) * qJDD(3) + t257); mrSges(4,1) * t225 - mrSges(4,2) * t226 + t248 * (mrSges(5,2) * t219 - mrSges(5,3) * t212 - pkin(6) * t268 - t252 * t199 + t254 * t200) + t250 * (-mrSges(5,1) * t219 + mrSges(5,3) * t213 - pkin(4) * t258 + pkin(6) * t263 + t254 * t199 + t252 * t200) - pkin(3) * t201 + qJ(4) * t264 + (Ifges(5,2) * t246 + Ifges(4,3) + (Ifges(5,1) * t248 + Ifges(5,4) * t272) * t248) * qJDD(3); t201; mrSges(6,1) * t207 - mrSges(6,2) * t208 + Ifges(6,5) * t224 + Ifges(6,6) * t223 + Ifges(6,3) * qJDD(5) + t215 * t230 - t216 * t229;];
tauJ = t1;
