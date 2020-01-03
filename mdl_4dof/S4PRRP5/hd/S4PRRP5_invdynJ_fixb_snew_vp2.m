% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4PRRP5
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4PRRP5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP5_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP5_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP5_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP5_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:28:54
% EndTime: 2019-12-31 16:28:55
% DurationCPUTime: 0.45s
% Computational Cost: add. (524->120), mult. (1000->147), div. (0->0), fcn. (477->6), ass. (0->56)
t237 = cos(qJ(3));
t254 = Ifges(4,4) + Ifges(5,4);
t261 = t237 * t254;
t260 = Ifges(4,1) + Ifges(5,1);
t253 = Ifges(4,5) + Ifges(5,5);
t259 = Ifges(4,2) + Ifges(5,2);
t258 = Ifges(4,6) + Ifges(5,6);
t235 = sin(qJ(3));
t257 = -t253 * qJD(3) + (-t260 * t235 - t261) * qJD(2);
t239 = qJD(2) ^ 2;
t256 = pkin(3) * t239;
t255 = -mrSges(4,2) - mrSges(5,2);
t233 = sin(pkin(6));
t234 = cos(pkin(6));
t223 = -g(1) * t234 - g(2) * t233;
t232 = -g(3) + qJDD(1);
t236 = sin(qJ(2));
t238 = cos(qJ(2));
t200 = t238 * t223 + t236 * t232;
t198 = -pkin(2) * t239 + qJDD(2) * pkin(5) + t200;
t222 = -g(1) * t233 + g(2) * t234;
t195 = t237 * t198 + t235 * t222;
t252 = t258 * qJD(3) + (t254 * t235 + t237 * t259) * qJD(2);
t249 = t235 * qJD(2);
t225 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t249;
t251 = -qJD(3) * mrSges(4,1) + mrSges(4,3) * t249 - t225;
t250 = qJD(2) * t237;
t248 = qJD(2) * qJD(3);
t247 = qJD(2) * qJD(4);
t215 = t237 * t222;
t244 = t237 * t248;
t219 = qJDD(2) * t235 + t244;
t191 = qJDD(3) * pkin(3) + t215 + (-t219 + t244) * qJ(4) + (t237 * t256 - t198 - 0.2e1 * t247) * t235;
t227 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t250;
t246 = m(5) * t191 + qJDD(3) * mrSges(5,1) + qJD(3) * t227;
t220 = qJDD(2) * t237 - t235 * t248;
t224 = qJD(3) * pkin(3) - qJ(4) * t249;
t231 = t237 ^ 2;
t192 = qJ(4) * t220 - qJD(3) * t224 - t231 * t256 + 0.2e1 * t237 * t247 + t195;
t217 = (-t237 * mrSges(5,1) + t235 * mrSges(5,2)) * qJD(2);
t245 = m(5) * t192 + t220 * mrSges(5,3) + t217 * t250;
t199 = -t236 * t223 + t232 * t238;
t243 = qJD(2) * t251;
t241 = -qJDD(2) * pkin(2) - t199;
t193 = t224 * t249 - pkin(3) * t220 + qJDD(4) + (-qJ(4) * t231 - pkin(5)) * t239 + t241;
t242 = m(5) * t193 - t220 * mrSges(5,1) - t227 * t250;
t197 = -pkin(5) * t239 + t241;
t228 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t250;
t240 = -m(4) * t197 + t220 * mrSges(4,1) + t228 * t250 - t242;
t218 = (-t237 * mrSges(4,1) + t235 * mrSges(4,2)) * qJD(2);
t194 = -t198 * t235 + t215;
t188 = mrSges(5,2) * t219 + t225 * t249 + t242;
t187 = -mrSges(5,3) * t219 - t217 * t249 + t246;
t186 = m(4) * t194 + qJDD(3) * mrSges(4,1) + qJD(3) * t228 + (-mrSges(4,3) - mrSges(5,3)) * t219 + (-t217 - t218) * t249 + t246;
t185 = t237 * (m(4) * t195 + mrSges(4,3) * t220 + t251 * qJD(3) + t255 * qJDD(3) + t218 * t250 + t245);
t1 = [m(2) * t232 + t236 * (m(3) * t200 - mrSges(3,1) * t239 - qJDD(2) * mrSges(3,2) - t186 * t235 + t185) + t238 * (m(3) * t199 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t239 + t255 * t219 + t235 * t243 + t240); Ifges(3,3) * qJDD(2) + mrSges(3,1) * t199 - mrSges(3,2) * t200 + t237 * (-mrSges(4,1) * t197 - mrSges(5,1) * t193 + mrSges(4,3) * t195 + mrSges(5,3) * t192 - pkin(3) * t188 + qJ(4) * t245 + t259 * t220 + (-qJ(4) * mrSges(5,2) + t258) * qJDD(3) + (-qJ(4) * t225 - t257) * qJD(3)) + pkin(2) * t240 + pkin(5) * t185 + (pkin(2) * t255 + t261) * t219 + (mrSges(4,2) * t197 + mrSges(5,2) * t193 - mrSges(4,3) * t194 - mrSges(5,3) * t191 + pkin(2) * t243 - pkin(5) * t186 - qJ(4) * t187 - t252 * qJD(3) + t253 * qJDD(3) + t260 * t219 + t254 * t220) * t235; mrSges(4,1) * t194 + mrSges(5,1) * t191 - mrSges(4,2) * t195 - mrSges(5,2) * t192 + pkin(3) * t187 + t258 * t220 + t253 * t219 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (t252 * t235 + t257 * t237) * qJD(2); t188;];
tauJ = t1;
