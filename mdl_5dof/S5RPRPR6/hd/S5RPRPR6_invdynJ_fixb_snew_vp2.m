% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRPR6
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRPR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR6_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR6_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR6_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR6_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR6_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:46
% EndTime: 2019-12-31 18:17:46
% DurationCPUTime: 0.31s
% Computational Cost: add. (1708->110), mult. (2279->140), div. (0->0), fcn. (1112->8), ass. (0->54)
t268 = -pkin(3) - pkin(7);
t246 = qJD(1) + qJD(3);
t250 = sin(qJ(5));
t267 = t246 * t250;
t253 = cos(qJ(5));
t266 = t246 * t253;
t252 = sin(qJ(1));
t255 = cos(qJ(1));
t263 = t252 * g(1) - t255 * g(2);
t234 = qJDD(1) * pkin(1) + t263;
t256 = qJD(1) ^ 2;
t262 = -t255 * g(1) - t252 * g(2);
t235 = -t256 * pkin(1) + t262;
t248 = sin(pkin(8));
t249 = cos(pkin(8));
t220 = t249 * t234 - t248 * t235;
t218 = qJDD(1) * pkin(2) + t220;
t221 = t248 * t234 + t249 * t235;
t219 = -t256 * pkin(2) + t221;
t251 = sin(qJ(3));
t254 = cos(qJ(3));
t213 = t254 * t218 - t251 * t219;
t244 = t246 ^ 2;
t245 = qJDD(1) + qJDD(3);
t260 = -t244 * qJ(4) + qJDD(4) - t213;
t208 = t268 * t245 + t260;
t247 = -g(3) + qJDD(2);
t204 = t253 * t208 - t250 * t247;
t228 = (mrSges(6,1) * t250 + mrSges(6,2) * t253) * t246;
t264 = qJD(5) * t246;
t230 = t253 * t245 - t250 * t264;
t236 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t267;
t202 = m(6) * t204 + qJDD(5) * mrSges(6,1) - t230 * mrSges(6,3) + qJD(5) * t236 - t228 * t266;
t205 = t250 * t208 + t253 * t247;
t229 = -t250 * t245 - t253 * t264;
t237 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t266;
t203 = m(6) * t205 - qJDD(5) * mrSges(6,2) + t229 * mrSges(6,3) - qJD(5) * t237 - t228 * t267;
t198 = t253 * t202 + t250 * t203;
t211 = -t245 * pkin(3) + t260;
t259 = -m(5) * t211 + t244 * mrSges(5,3) - t198;
t193 = m(4) * t213 - t244 * mrSges(4,2) + (mrSges(4,1) - mrSges(5,2)) * t245 + t259;
t214 = t251 * t218 + t254 * t219;
t261 = t245 * qJ(4) + 0.2e1 * qJD(4) * t246 + t214;
t207 = t268 * t244 + t261;
t209 = t244 * pkin(3) - t261;
t258 = -m(5) * t209 + m(6) * t207 - t229 * mrSges(6,1) + t244 * mrSges(5,2) + t230 * mrSges(6,2) + t245 * mrSges(5,3) + t236 * t267 + t237 * t266;
t197 = m(4) * t214 - t244 * mrSges(4,1) - t245 * mrSges(4,2) + t258;
t265 = t254 * t193 + t251 * t197;
t195 = t245 * mrSges(5,2) - t259;
t222 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t253 - Ifges(6,6) * t250) * t246;
t223 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t253 - Ifges(6,2) * t250) * t246;
t224 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t253 - Ifges(6,4) * t250) * t246;
t257 = -mrSges(4,2) * t214 - mrSges(5,3) * t209 - pkin(3) * t195 - pkin(7) * t198 - t250 * (-mrSges(6,1) * t207 + mrSges(6,3) * t205 + Ifges(6,4) * t230 + Ifges(6,2) * t229 + Ifges(6,6) * qJDD(5) + qJD(5) * t224 - t222 * t266) + t253 * (mrSges(6,2) * t207 - mrSges(6,3) * t204 + Ifges(6,1) * t230 + Ifges(6,4) * t229 + Ifges(6,5) * qJDD(5) - qJD(5) * t223 - t222 * t267) + qJ(4) * t258 + mrSges(5,2) * t211 + mrSges(4,1) * t213 + (Ifges(4,3) + Ifges(5,1)) * t245;
t1 = [t257 + mrSges(2,1) * t263 - mrSges(2,2) * t262 + pkin(1) * (t248 * (m(3) * t221 - t256 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t251 * t193 + t254 * t197) + t249 * (m(3) * t220 + qJDD(1) * mrSges(3,1) - t256 * mrSges(3,2) + t265)) + mrSges(3,1) * t220 - mrSges(3,2) * t221 + pkin(2) * t265 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * qJDD(1); -t250 * t202 + t253 * t203 + (m(3) + m(4) + m(5)) * t247; t257; t195; mrSges(6,1) * t204 - mrSges(6,2) * t205 + Ifges(6,5) * t230 + Ifges(6,6) * t229 + Ifges(6,3) * qJDD(5) + (t223 * t253 + t224 * t250) * t246;];
tauJ = t1;
