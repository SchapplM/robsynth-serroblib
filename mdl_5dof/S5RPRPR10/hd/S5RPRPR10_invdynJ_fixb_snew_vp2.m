% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRPR10
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRPR10_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR10_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR10_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR10_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR10_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR10_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:26:00
% EndTime: 2019-12-31 18:26:00
% DurationCPUTime: 0.37s
% Computational Cost: add. (3201->111), mult. (4228->144), div. (0->0), fcn. (1514->8), ass. (0->56)
t251 = qJD(1) ^ 2;
t247 = sin(qJ(1));
t250 = cos(qJ(1));
t256 = -t250 * g(1) - t247 * g(2);
t254 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t256;
t262 = -pkin(1) - pkin(2);
t222 = t262 * t251 + t254;
t257 = t247 * g(1) - t250 * g(2);
t253 = -t251 * qJ(2) + qJDD(2) - t257;
t223 = t262 * qJDD(1) + t253;
t246 = sin(qJ(3));
t249 = cos(qJ(3));
t217 = -t246 * t222 + t249 * t223;
t239 = -qJDD(1) + qJDD(3);
t215 = t239 * pkin(3) + t217;
t218 = t249 * t222 + t246 * t223;
t240 = -qJD(1) + qJD(3);
t238 = t240 ^ 2;
t216 = -t238 * pkin(3) + t218;
t243 = sin(pkin(8));
t244 = cos(pkin(8));
t212 = t243 * t215 + t244 * t216;
t210 = -(t238 * pkin(4)) + t239 * pkin(7) + t212;
t242 = g(3) + qJDD(4);
t245 = sin(qJ(5));
t248 = cos(qJ(5));
t207 = -t245 * t210 + t248 * t242;
t231 = (-mrSges(6,1) * t248 + mrSges(6,2) * t245) * t240;
t258 = qJD(5) * t240;
t232 = t245 * t239 + t248 * t258;
t259 = t240 * t248;
t235 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t259;
t260 = t240 * t245;
t205 = m(6) * t207 + qJDD(5) * mrSges(6,1) - t232 * mrSges(6,3) + qJD(5) * t235 - t231 * t260;
t208 = t248 * t210 + t245 * t242;
t233 = t248 * t239 - t245 * t258;
t234 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t260;
t206 = m(6) * t208 - qJDD(5) * mrSges(6,2) + t233 * mrSges(6,3) - qJD(5) * t234 + t231 * t259;
t198 = -t245 * t205 + t248 * t206;
t197 = m(5) * t212 - (t238 * mrSges(5,1)) - t239 * mrSges(5,2) + t198;
t211 = t244 * t215 - t243 * t216;
t209 = -t239 * pkin(4) - t238 * pkin(7) - t211;
t203 = -m(6) * t209 + t233 * mrSges(6,1) - t232 * mrSges(6,2) - t234 * t260 + t235 * t259;
t202 = m(5) * t211 + t239 * mrSges(5,1) - t238 * mrSges(5,2) + t203;
t195 = t243 * t197 + t244 * t202;
t224 = (Ifges(6,3) * qJD(5)) + (Ifges(6,5) * t245 + Ifges(6,6) * t248) * t240;
t225 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t245 + Ifges(6,2) * t248) * t240;
t226 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t245 + Ifges(6,4) * t248) * t240;
t263 = (Ifges(4,3) + Ifges(5,3)) * t239 + mrSges(4,1) * t217 + mrSges(5,1) * t211 - mrSges(4,2) * t218 - mrSges(5,2) * t212 + pkin(3) * t195 + pkin(4) * t203 + pkin(7) * t198 + t248 * (-mrSges(6,1) * t209 + mrSges(6,3) * t208 + Ifges(6,4) * t232 + Ifges(6,2) * t233 + Ifges(6,6) * qJDD(5) + qJD(5) * t226 - t224 * t260) + t245 * (mrSges(6,2) * t209 - mrSges(6,3) * t207 + Ifges(6,1) * t232 + Ifges(6,4) * t233 + Ifges(6,5) * qJDD(5) - qJD(5) * t225 + t224 * t259);
t193 = m(4) * t217 + t239 * mrSges(4,1) - t238 * mrSges(4,2) + t195;
t194 = m(4) * t218 - t238 * mrSges(4,1) - t239 * mrSges(4,2) + t244 * t197 - t243 * t202;
t255 = t249 * t193 + t246 * t194;
t230 = -qJDD(1) * pkin(1) + t253;
t227 = -t251 * pkin(1) + t254;
t192 = m(3) * t230 - qJDD(1) * mrSges(3,1) - t251 * mrSges(3,3) + t255;
t1 = [mrSges(3,3) * t227 - mrSges(3,1) * t230 + mrSges(2,1) * t257 - pkin(2) * t255 - mrSges(2,2) * t256 - pkin(1) * t192 + qJ(2) * (m(3) * t227 - t251 * mrSges(3,1) - t246 * t193 + t249 * t194) + (qJ(2) * mrSges(3,3) + Ifges(3,2) + Ifges(2,3)) * qJDD(1) - t263; t192; t263; m(5) * t242 + t248 * t205 + t245 * t206; mrSges(6,1) * t207 - mrSges(6,2) * t208 + Ifges(6,5) * t232 + Ifges(6,6) * t233 + Ifges(6,3) * qJDD(5) + (t225 * t245 - t226 * t248) * t240;];
tauJ = t1;
