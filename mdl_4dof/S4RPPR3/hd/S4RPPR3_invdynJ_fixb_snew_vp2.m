% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RPPR3
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
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RPPR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR3_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR3_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR3_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR3_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR3_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:53
% EndTime: 2019-12-31 16:37:53
% DurationCPUTime: 0.38s
% Computational Cost: add. (1244->118), mult. (2611->159), div. (0->0), fcn. (1524->8), ass. (0->60)
t253 = qJD(1) ^ 2;
t247 = cos(pkin(7));
t271 = t247 * pkin(3);
t245 = sin(pkin(7));
t270 = mrSges(4,2) * t245;
t243 = t247 ^ 2;
t269 = t243 * t253;
t250 = sin(qJ(1));
t252 = cos(qJ(1));
t263 = t250 * g(1) - t252 * g(2);
t232 = qJDD(1) * pkin(1) + t263;
t260 = -t252 * g(1) - t250 * g(2);
t233 = -t253 * pkin(1) + t260;
t246 = sin(pkin(6));
t248 = cos(pkin(6));
t222 = t246 * t232 + t248 * t233;
t215 = -t253 * pkin(2) + qJDD(1) * qJ(3) + t222;
t244 = -g(3) + qJDD(2);
t265 = qJD(1) * qJD(3);
t267 = t247 * t244 - 0.2e1 * t245 * t265;
t205 = (-pkin(5) * qJDD(1) + t253 * t271 - t215) * t245 + t267;
t209 = t245 * t244 + (t215 + 0.2e1 * t265) * t247;
t264 = t247 * qJDD(1);
t206 = -pkin(3) * t269 + pkin(5) * t264 + t209;
t249 = sin(qJ(4));
t251 = cos(qJ(4));
t203 = t251 * t205 - t249 * t206;
t257 = -t245 * t249 + t247 * t251;
t225 = t257 * qJD(1);
t258 = t245 * t251 + t247 * t249;
t226 = t258 * qJD(1);
t217 = -t225 * mrSges(5,1) + t226 * mrSges(5,2);
t220 = t225 * qJD(4) + t258 * qJDD(1);
t223 = -qJD(4) * mrSges(5,2) + t225 * mrSges(5,3);
t201 = m(5) * t203 + qJDD(4) * mrSges(5,1) - t220 * mrSges(5,3) + qJD(4) * t223 - t226 * t217;
t204 = t249 * t205 + t251 * t206;
t219 = -t226 * qJD(4) + t257 * qJDD(1);
t224 = qJD(4) * mrSges(5,1) - t226 * mrSges(5,3);
t202 = m(5) * t204 - qJDD(4) * mrSges(5,2) + t219 * mrSges(5,3) - qJD(4) * t224 + t225 * t217;
t268 = t251 * t201 + t249 * t202;
t208 = -t245 * t215 + t267;
t256 = mrSges(4,3) * qJDD(1) + t253 * (-t247 * mrSges(4,1) + t270);
t193 = m(4) * t208 - t256 * t245 + t268;
t261 = -t249 * t201 + t251 * t202;
t194 = m(4) * t209 + t256 * t247 + t261;
t262 = -t245 * t193 + t247 * t194;
t221 = t248 * t232 - t246 * t233;
t259 = qJDD(3) - t221;
t242 = t245 ^ 2;
t207 = (-pkin(2) - t271) * qJDD(1) + (-qJ(3) + (-t242 - t243) * pkin(5)) * t253 + t259;
t255 = m(5) * t207 - t219 * mrSges(5,1) + t220 * mrSges(5,2) - t225 * t223 + t226 * t224;
t211 = -qJDD(1) * pkin(2) - t253 * qJ(3) + t259;
t254 = -m(4) * t211 + mrSges(4,1) * t264 - t255 + (t242 * t253 + t269) * mrSges(4,3);
t214 = Ifges(5,1) * t226 + Ifges(5,4) * t225 + Ifges(5,5) * qJD(4);
t213 = Ifges(5,4) * t226 + Ifges(5,2) * t225 + Ifges(5,6) * qJD(4);
t212 = Ifges(5,5) * t226 + Ifges(5,6) * t225 + Ifges(5,3) * qJD(4);
t197 = qJDD(1) * t270 - t254;
t196 = mrSges(5,2) * t207 - mrSges(5,3) * t203 + Ifges(5,1) * t220 + Ifges(5,4) * t219 + Ifges(5,5) * qJDD(4) - qJD(4) * t213 + t225 * t212;
t195 = -mrSges(5,1) * t207 + mrSges(5,3) * t204 + Ifges(5,4) * t220 + Ifges(5,2) * t219 + Ifges(5,6) * qJDD(4) + qJD(4) * t214 - t226 * t212;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t263 - mrSges(2,2) * t260 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t221 - mrSges(3,2) * t222 + t245 * (mrSges(4,2) * t211 - mrSges(4,3) * t208 + t251 * t196 - t249 * t195 - pkin(5) * t268 + (Ifges(4,1) * t245 + Ifges(4,4) * t247) * qJDD(1)) + t247 * (-mrSges(4,1) * t211 + mrSges(4,3) * t209 + t249 * t196 + t251 * t195 - pkin(3) * t255 + pkin(5) * t261 + (Ifges(4,4) * t245 + Ifges(4,2) * t247) * qJDD(1)) - pkin(2) * t197 + qJ(3) * t262 + pkin(1) * (t246 * (m(3) * t222 - t253 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t262) + t248 * (m(3) * t221 - t253 * mrSges(3,2) + (mrSges(3,1) - t270) * qJDD(1) + t254)); m(3) * t244 + t247 * t193 + t245 * t194; t197; mrSges(5,1) * t203 - mrSges(5,2) * t204 + Ifges(5,5) * t220 + Ifges(5,6) * t219 + Ifges(5,3) * qJDD(4) + t226 * t213 - t225 * t214;];
tauJ = t1;
