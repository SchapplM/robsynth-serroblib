% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RRRP3
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
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RRRP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP3_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP3_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP3_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP3_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:04
% EndTime: 2019-12-31 17:14:05
% DurationCPUTime: 0.39s
% Computational Cost: add. (1357->124), mult. (1700->157), div. (0->0), fcn. (760->6), ass. (0->54)
t270 = Ifges(4,1) + Ifges(5,1);
t265 = Ifges(4,4) - Ifges(5,5);
t264 = Ifges(4,5) + Ifges(5,4);
t269 = Ifges(4,2) + Ifges(5,3);
t263 = Ifges(4,6) - Ifges(5,6);
t268 = Ifges(4,3) + Ifges(5,2);
t246 = cos(qJ(3));
t267 = t246 * g(3);
t266 = mrSges(4,3) + mrSges(5,2);
t240 = qJD(1) + qJD(2);
t243 = sin(qJ(3));
t262 = t240 * t243;
t261 = t240 * t246;
t260 = (-t243 * t265 - t269 * t246) * t240 - t263 * qJD(3);
t259 = (t243 * t264 + t246 * t263) * t240 + t268 * qJD(3);
t258 = (t270 * t243 + t246 * t265) * t240 + t264 * qJD(3);
t245 = sin(qJ(1));
t248 = cos(qJ(1));
t256 = t245 * g(1) - t248 * g(2);
t229 = qJDD(1) * pkin(1) + t256;
t252 = -t248 * g(1) - t245 * g(2);
t230 = -qJD(1) ^ 2 * pkin(1) + t252;
t244 = sin(qJ(2));
t247 = cos(qJ(2));
t206 = t244 * t229 + t247 * t230;
t257 = qJD(3) * t240;
t238 = t240 ^ 2;
t239 = qJDD(1) + qJDD(2);
t203 = -t238 * pkin(2) + t239 * pkin(6) + t206;
t200 = -t243 * g(3) + t246 * t203;
t199 = -t243 * t203 - t267;
t220 = (-mrSges(5,1) * t246 - mrSges(5,3) * t243) * t240;
t221 = (-mrSges(4,1) * t246 + mrSges(4,2) * t243) * t240;
t222 = t243 * t239 + t246 * t257;
t223 = t246 * t239 - t243 * t257;
t231 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t262;
t233 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t261;
t219 = (-pkin(3) * t246 - qJ(4) * t243) * t240;
t249 = qJD(3) ^ 2;
t198 = -qJDD(3) * pkin(3) + t267 - t249 * qJ(4) + qJDD(4) + (t219 * t240 + t203) * t243;
t234 = mrSges(5,2) * t261 + qJD(3) * mrSges(5,3);
t253 = -m(5) * t198 + qJDD(3) * mrSges(5,1) + qJD(3) * t234;
t197 = -t249 * pkin(3) + qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) + t219 * t261 + t200;
t232 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t262;
t254 = m(5) * t197 + qJDD(3) * mrSges(5,3) + qJD(3) * t232 + t220 * t261;
t255 = -t243 * (m(4) * t199 + qJDD(3) * mrSges(4,1) + qJD(3) * t233 + (-t220 - t221) * t262 - t266 * t222 + t253) + t246 * (m(4) * t200 - qJDD(3) * mrSges(4,2) - qJD(3) * t231 + t221 * t261 + t266 * t223 + t254);
t205 = t247 * t229 - t244 * t230;
t202 = -t239 * pkin(2) - t238 * pkin(6) - t205;
t195 = -t223 * pkin(3) - t222 * qJ(4) + (-0.2e1 * qJD(4) * t243 + (pkin(3) * t243 - qJ(4) * t246) * qJD(3)) * t240 + t202;
t193 = m(5) * t195 - t223 * mrSges(5,1) - t222 * mrSges(5,3) - t232 * t262 - t234 * t261;
t250 = -m(4) * t202 + t223 * mrSges(4,1) - t222 * mrSges(4,2) - t231 * t262 + t233 * t261 - t193;
t251 = -mrSges(3,2) * t206 + t246 * (-mrSges(4,1) * t202 - mrSges(5,1) * t195 + mrSges(5,2) * t197 + mrSges(4,3) * t200 - pkin(3) * t193 + t258 * qJD(3) + t263 * qJDD(3) + t265 * t222 + t269 * t223 - t259 * t262) + t243 * (mrSges(4,2) * t202 + mrSges(5,2) * t198 - mrSges(4,3) * t199 - mrSges(5,3) * t195 - qJ(4) * t193 + t260 * qJD(3) + t264 * qJDD(3) + t270 * t222 + t265 * t223 + t259 * t261) + pkin(6) * t255 + pkin(2) * t250 + mrSges(3,1) * t205 + Ifges(3,3) * t239;
t194 = t222 * mrSges(5,2) + t220 * t262 - t253;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t256 - mrSges(2,2) * t252 + pkin(1) * (t244 * (m(3) * t206 - t238 * mrSges(3,1) - t239 * mrSges(3,2) + t255) + t247 * (m(3) * t205 + t239 * mrSges(3,1) - t238 * mrSges(3,2) + t250)) + t251; t251; mrSges(4,1) * t199 - mrSges(4,2) * t200 - mrSges(5,1) * t198 + mrSges(5,3) * t197 - pkin(3) * t194 + qJ(4) * t254 + (qJ(4) * mrSges(5,2) + t263) * t223 + t264 * t222 + t268 * qJDD(3) + (-t243 * t260 - t246 * t258) * t240; t194;];
tauJ = t1;
