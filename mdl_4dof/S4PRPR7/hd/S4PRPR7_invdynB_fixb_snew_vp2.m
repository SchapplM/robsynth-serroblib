% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4PRPR7
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
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4PRPR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR7_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR7_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR7_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR7_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR7_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:39
% EndTime: 2019-12-31 16:25:39
% DurationCPUTime: 0.35s
% Computational Cost: add. (2248->137), mult. (3710->169), div. (0->0), fcn. (1786->6), ass. (0->59)
t269 = m(3) + m(4);
t268 = -pkin(2) - pkin(5);
t267 = mrSges(3,1) - mrSges(4,2);
t266 = -Ifges(4,4) + Ifges(3,5);
t265 = Ifges(4,5) - Ifges(3,6);
t264 = cos(pkin(6));
t246 = sin(pkin(6));
t237 = -t264 * g(1) - t246 * g(2);
t243 = -g(3) + qJDD(1);
t248 = sin(qJ(2));
t250 = cos(qJ(2));
t223 = -t248 * t237 + t250 * t243;
t251 = qJD(2) ^ 2;
t255 = -t251 * qJ(3) + qJDD(3) - t223;
t220 = t268 * qJDD(2) + t255;
t236 = t246 * g(1) - t264 * g(2);
t247 = sin(qJ(4));
t249 = cos(qJ(4));
t216 = t249 * t220 + t247 * t236;
t233 = (mrSges(5,1) * t247 + mrSges(5,2) * t249) * qJD(2);
t259 = qJD(2) * qJD(4);
t235 = t249 * qJDD(2) - t247 * t259;
t261 = qJD(2) * t247;
t238 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t261;
t260 = qJD(2) * t249;
t214 = m(5) * t216 + qJDD(4) * mrSges(5,1) - t235 * mrSges(5,3) + qJD(4) * t238 - t233 * t260;
t217 = t247 * t220 - t249 * t236;
t234 = -t247 * qJDD(2) - t249 * t259;
t239 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t260;
t215 = m(5) * t217 - qJDD(4) * mrSges(5,2) + t234 * mrSges(5,3) - qJD(4) * t239 - t233 * t261;
t206 = t249 * t214 + t247 * t215;
t222 = -qJDD(2) * pkin(2) + t255;
t253 = -m(4) * t222 + t251 * mrSges(4,3) - t206;
t202 = m(3) * t223 - t251 * mrSges(3,2) + t267 * qJDD(2) + t253;
t224 = t250 * t237 + t248 * t243;
t254 = qJDD(2) * qJ(3) + 0.2e1 * qJD(3) * qJD(2) + t224;
t221 = t251 * pkin(2) - t254;
t219 = t268 * t251 + t254;
t256 = -m(5) * t219 + t234 * mrSges(5,1) - t235 * mrSges(5,2) - t238 * t261 - t239 * t260;
t252 = -m(4) * t221 + t251 * mrSges(4,2) + qJDD(2) * mrSges(4,3) - t256;
t209 = m(3) * t224 - t251 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t252;
t257 = -t248 * t202 + t250 * t209;
t199 = m(2) * t237 + t257;
t262 = -t247 * t214 + t249 * t215;
t204 = (m(2) + t269) * t236 - t262;
t263 = t246 * t199 + t264 * t204;
t200 = t250 * t202 + t248 * t209;
t258 = t264 * t199 - t246 * t204;
t227 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t249 - Ifges(5,4) * t247) * qJD(2);
t226 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t249 - Ifges(5,2) * t247) * qJD(2);
t225 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t249 - Ifges(5,6) * t247) * qJD(2);
t211 = mrSges(5,2) * t219 - mrSges(5,3) * t216 + Ifges(5,1) * t235 + Ifges(5,4) * t234 + Ifges(5,5) * qJDD(4) - qJD(4) * t226 - t225 * t261;
t210 = -mrSges(5,1) * t219 + mrSges(5,3) * t217 + Ifges(5,4) * t235 + Ifges(5,2) * t234 + Ifges(5,6) * qJDD(4) + qJD(4) * t227 - t225 * t260;
t205 = -m(4) * t236 + t262;
t196 = mrSges(4,1) * t222 + mrSges(5,1) * t216 - mrSges(5,2) * t217 - mrSges(3,3) * t223 + Ifges(5,5) * t235 + Ifges(5,6) * t234 + Ifges(5,3) * qJDD(4) + pkin(3) * t206 - qJ(3) * t205 + t265 * t251 + (-mrSges(3,2) + mrSges(4,3)) * t236 + t266 * qJDD(2) + (t249 * t226 + t247 * t227) * qJD(2);
t195 = -mrSges(4,1) * t221 + mrSges(3,3) * t224 - pkin(2) * t205 - pkin(3) * t256 - pkin(5) * t262 - t265 * qJDD(2) - t249 * t210 - t247 * t211 + t267 * t236 + t266 * t251;
t194 = -pkin(1) * t200 + mrSges(2,3) * t237 - pkin(2) * t253 - qJ(3) * t252 - mrSges(4,2) * t222 + mrSges(4,3) * t221 - t249 * t211 + t247 * t210 + pkin(5) * t206 - mrSges(3,1) * t223 + mrSges(3,2) * t224 - mrSges(2,1) * t243 + (pkin(2) * mrSges(4,2) - Ifges(4,1) - Ifges(3,3)) * qJDD(2);
t193 = mrSges(2,2) * t243 - mrSges(2,3) * t236 - pkin(4) * t200 - t248 * t195 + t250 * t196;
t1 = [-m(1) * g(1) + t258; -m(1) * g(2) + t263; -m(1) * g(3) + m(2) * t243 + t200; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t263 + t264 * t193 - t246 * t194; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t258 + t246 * t193 + t264 * t194; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - mrSges(2,2) * t237 + t248 * t196 + t250 * t195 - pkin(1) * t262 + pkin(4) * t257 + (pkin(1) * t269 + mrSges(2,1)) * t236;];
tauB = t1;
