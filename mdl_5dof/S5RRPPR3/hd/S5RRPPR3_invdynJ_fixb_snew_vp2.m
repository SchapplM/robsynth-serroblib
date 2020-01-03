% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPPR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR3_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR3_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR3_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR3_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR3_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:32
% EndTime: 2019-12-31 19:26:32
% DurationCPUTime: 0.32s
% Computational Cost: add. (2048->109), mult. (2559->139), div. (0->0), fcn. (1250->8), ass. (0->53)
t269 = -pkin(3) - pkin(7);
t248 = qJD(1) + qJD(2);
t252 = sin(qJ(5));
t268 = t248 * t252;
t255 = cos(qJ(5));
t267 = t248 * t255;
t254 = sin(qJ(1));
t257 = cos(qJ(1));
t264 = t254 * g(1) - t257 * g(2);
t235 = qJDD(1) * pkin(1) + t264;
t263 = -t257 * g(1) - t254 * g(2);
t236 = -qJD(1) ^ 2 * pkin(1) + t263;
t253 = sin(qJ(2));
t256 = cos(qJ(2));
t221 = t256 * t235 - t253 * t236;
t247 = qJDD(1) + qJDD(2);
t218 = t247 * pkin(2) + t221;
t222 = t253 * t235 + t256 * t236;
t246 = t248 ^ 2;
t219 = -t246 * pkin(2) + t222;
t250 = sin(pkin(8));
t251 = cos(pkin(8));
t213 = t251 * t218 - t250 * t219;
t261 = -t246 * qJ(4) + qJDD(4) - t213;
t208 = t269 * t247 + t261;
t249 = -g(3) + qJDD(3);
t204 = t255 * t208 - t252 * t249;
t229 = (mrSges(6,1) * t252 + mrSges(6,2) * t255) * t248;
t265 = qJD(5) * t248;
t231 = t255 * t247 - t252 * t265;
t237 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t268;
t202 = m(6) * t204 + qJDD(5) * mrSges(6,1) - t231 * mrSges(6,3) + qJD(5) * t237 - t229 * t267;
t205 = t252 * t208 + t255 * t249;
t230 = -t252 * t247 - t255 * t265;
t238 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t267;
t203 = m(6) * t205 - qJDD(5) * mrSges(6,2) + t230 * mrSges(6,3) - qJD(5) * t238 - t229 * t268;
t198 = t255 * t202 + t252 * t203;
t211 = -t247 * pkin(3) + t261;
t260 = -m(5) * t211 + t246 * mrSges(5,3) - t198;
t193 = m(4) * t213 - t246 * mrSges(4,2) + (mrSges(4,1) - mrSges(5,2)) * t247 + t260;
t214 = t250 * t218 + t251 * t219;
t262 = t247 * qJ(4) + 0.2e1 * qJD(4) * t248 + t214;
t207 = t269 * t246 + t262;
t209 = t246 * pkin(3) - t262;
t259 = -m(5) * t209 + m(6) * t207 - t230 * mrSges(6,1) + t246 * mrSges(5,2) + t231 * mrSges(6,2) + t247 * mrSges(5,3) + t237 * t268 + t238 * t267;
t197 = m(4) * t214 - t246 * mrSges(4,1) - t247 * mrSges(4,2) + t259;
t266 = t251 * t193 + t250 * t197;
t195 = t247 * mrSges(5,2) - t260;
t223 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t255 - Ifges(6,6) * t252) * t248;
t224 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t255 - Ifges(6,2) * t252) * t248;
t225 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t255 - Ifges(6,4) * t252) * t248;
t258 = -mrSges(3,2) * t222 - mrSges(4,2) * t214 - mrSges(5,3) * t209 + pkin(2) * t266 - pkin(3) * t195 - pkin(7) * t198 - t252 * (-mrSges(6,1) * t207 + mrSges(6,3) * t205 + Ifges(6,4) * t231 + Ifges(6,2) * t230 + Ifges(6,6) * qJDD(5) + qJD(5) * t225 - t223 * t267) + t255 * (mrSges(6,2) * t207 - mrSges(6,3) * t204 + Ifges(6,1) * t231 + Ifges(6,4) * t230 + Ifges(6,5) * qJDD(5) - qJD(5) * t224 - t223 * t268) + qJ(4) * t259 + mrSges(5,2) * t211 + mrSges(4,1) * t213 + mrSges(3,1) * t221 + (Ifges(4,3) + Ifges(3,3) + Ifges(5,1)) * t247;
t1 = [t258 + Ifges(2,3) * qJDD(1) + pkin(1) * (t253 * (m(3) * t222 - t246 * mrSges(3,1) - t247 * mrSges(3,2) - t250 * t193 + t251 * t197) + t256 * (m(3) * t221 + t247 * mrSges(3,1) - t246 * mrSges(3,2) + t266)) + mrSges(2,1) * t264 - mrSges(2,2) * t263; t258; -t252 * t202 + t255 * t203 + (m(4) + m(5)) * t249; t195; mrSges(6,1) * t204 - mrSges(6,2) * t205 + Ifges(6,5) * t231 + Ifges(6,6) * t230 + Ifges(6,3) * qJDD(5) + (t224 * t255 + t225 * t252) * t248;];
tauJ = t1;
