% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4PRRP3
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
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4PRRP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP3_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP3_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP3_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP3_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:46
% EndTime: 2019-12-31 16:26:47
% DurationCPUTime: 0.42s
% Computational Cost: add. (480->112), mult. (938->138), div. (0->0), fcn. (458->6), ass. (0->53)
t235 = sin(qJ(3));
t237 = cos(qJ(3));
t255 = Ifges(4,4) + Ifges(5,4);
t263 = t235 * t255 + t237 * (Ifges(4,2) + Ifges(5,2));
t262 = t235 * (Ifges(4,1) + Ifges(5,1)) + t237 * t255;
t254 = Ifges(4,5) + Ifges(5,5);
t253 = Ifges(4,6) + Ifges(5,6);
t259 = (t263 * qJD(2) + t253 * qJD(3)) * t235;
t246 = qJD(2) * qJD(3);
t219 = t237 * qJDD(2) - t235 * t246;
t247 = t235 * qJD(2);
t223 = qJD(3) * pkin(3) - qJ(4) * t247;
t231 = t237 ^ 2;
t239 = qJD(2) ^ 2;
t233 = sin(pkin(6));
t234 = cos(pkin(6));
t221 = t233 * g(1) - t234 * g(2);
t222 = -t234 * g(1) - t233 * g(2);
t236 = sin(qJ(2));
t238 = cos(qJ(2));
t241 = t238 * t221 - t236 * t222;
t240 = -qJDD(2) * pkin(2) - t241;
t199 = t223 * t247 - t219 * pkin(3) + qJDD(4) + (-qJ(4) * t231 - pkin(5)) * t239 + t240;
t258 = m(5) * t199;
t257 = pkin(3) * t239;
t256 = -mrSges(4,2) - mrSges(5,2);
t250 = t236 * t221 + t238 * t222;
t204 = -t239 * pkin(2) + qJDD(2) * pkin(5) + t250;
t232 = -g(3) + qJDD(1);
t201 = t237 * t204 + t235 * t232;
t251 = -t262 * qJD(2) - t254 * qJD(3);
t224 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t247;
t249 = -qJD(3) * mrSges(4,1) + mrSges(4,3) * t247 - t224;
t248 = qJD(2) * t237;
t245 = qJD(2) * qJD(4);
t242 = t237 * t246;
t218 = t235 * qJDD(2) + t242;
t229 = t237 * t232;
t197 = qJDD(3) * pkin(3) + t229 + (-t218 + t242) * qJ(4) + (t237 * t257 - t204 - 0.2e1 * t245) * t235;
t226 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t248;
t244 = m(5) * t197 + qJDD(3) * mrSges(5,1) + qJD(3) * t226;
t198 = t219 * qJ(4) - qJD(3) * t223 - t231 * t257 + 0.2e1 * t237 * t245 + t201;
t216 = (-t237 * mrSges(5,1) + t235 * mrSges(5,2)) * qJD(2);
t243 = m(5) * t198 + t219 * mrSges(5,3) + t216 * t248;
t227 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t248;
t217 = (-t237 * mrSges(4,1) + t235 * mrSges(4,2)) * qJD(2);
t203 = -t239 * pkin(5) + t240;
t200 = -t235 * t204 + t229;
t194 = t258 - t219 * mrSges(5,1) + t218 * mrSges(5,2) + (t224 * t235 - t226 * t237) * qJD(2);
t193 = -t218 * mrSges(5,3) - t216 * t247 + t244;
t192 = m(4) * t201 + t219 * mrSges(4,3) + t249 * qJD(3) + t256 * qJDD(3) + t217 * t248 + t243;
t191 = m(4) * t200 + qJDD(3) * mrSges(4,1) + qJD(3) * t227 + (-mrSges(4,3) - mrSges(5,3)) * t218 + (-t216 - t217) * t247 + t244;
t1 = [t237 * t191 + t235 * t192 + (m(2) + m(3)) * t232; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t241 - mrSges(3,2) * t250 + t235 * (mrSges(4,2) * t203 + mrSges(5,2) * t199 - mrSges(4,3) * t200 - mrSges(5,3) * t197 - qJ(4) * t193) + t237 * (-mrSges(4,1) * t203 - mrSges(5,1) * t199 + mrSges(4,3) * t201 + mrSges(5,3) * t198 - pkin(3) * t194 + qJ(4) * t243) + pkin(5) * (-t235 * t191 + t237 * t192) + (t235 * t254 + t237 * (-qJ(4) * mrSges(5,2) + t253)) * qJDD(3) + (-t259 + t237 * (-qJ(4) * t224 - t251)) * qJD(3) + t263 * t219 + t262 * t218 + (-m(4) * t203 - t258 + (mrSges(4,1) + mrSges(5,1)) * t219 + t256 * t218 + ((t226 + t227) * t237 + t249 * t235) * qJD(2)) * pkin(2); mrSges(4,1) * t200 + mrSges(5,1) * t197 - mrSges(4,2) * t201 - mrSges(5,2) * t198 + pkin(3) * t193 + t253 * t219 + t254 * t218 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (t251 * t237 + t259) * qJD(2); t194;];
tauJ = t1;
