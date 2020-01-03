% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RPRP4
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
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RPRP4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP4_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP4_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP4_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP4_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP4_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:40
% EndTime: 2019-12-31 16:43:41
% DurationCPUTime: 0.56s
% Computational Cost: add. (668->118), mult. (1245->151), div. (0->0), fcn. (557->6), ass. (0->51)
t237 = sin(qJ(3));
t239 = cos(qJ(3));
t259 = Ifges(4,4) - Ifges(5,5);
t263 = t237 * (Ifges(4,1) + Ifges(5,1)) + t239 * t259;
t262 = t237 * t259 + t239 * (Ifges(4,2) + Ifges(5,3));
t256 = Ifges(5,4) + Ifges(4,5);
t255 = Ifges(4,6) - Ifges(5,6);
t258 = t237 * (-qJD(1) * t262 - t255 * qJD(3)) + t239 * (t263 * qJD(1) + t256 * qJD(3));
t257 = mrSges(4,3) + mrSges(5,2);
t234 = -g(3) + qJDD(2);
t254 = t234 * t239;
t238 = sin(qJ(1));
t240 = cos(qJ(1));
t248 = t238 * g(1) - g(2) * t240;
t216 = qJDD(1) * pkin(1) + t248;
t242 = qJD(1) ^ 2;
t244 = -g(1) * t240 - g(2) * t238;
t220 = -pkin(1) * t242 + t244;
t235 = sin(pkin(6));
t236 = cos(pkin(6));
t200 = t235 * t216 + t236 * t220;
t198 = -pkin(2) * t242 + qJDD(1) * pkin(5) + t200;
t195 = t239 * t198 + t237 * t234;
t251 = qJD(1) * t237;
t250 = qJD(1) * t239;
t249 = qJD(1) * qJD(3);
t219 = (-mrSges(4,1) * t239 + mrSges(4,2) * t237) * qJD(1);
t222 = qJDD(1) * t239 - t237 * t249;
t225 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t251;
t217 = (-pkin(3) * t239 - qJ(4) * t237) * qJD(1);
t241 = qJD(3) ^ 2;
t192 = -pkin(3) * t241 + qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) + t217 * t250 + t195;
t218 = (-mrSges(5,1) * t239 - mrSges(5,3) * t237) * qJD(1);
t226 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t251;
t246 = m(5) * t192 + qJDD(3) * mrSges(5,3) + qJD(3) * t226 + t218 * t250;
t186 = m(4) * t195 - qJDD(3) * mrSges(4,2) - qJD(3) * t225 + t219 * t250 + t257 * t222 + t246;
t194 = -t198 * t237 + t254;
t221 = qJDD(1) * t237 + t239 * t249;
t227 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t250;
t193 = -qJDD(3) * pkin(3) - qJ(4) * t241 - t254 + qJDD(4) + (qJD(1) * t217 + t198) * t237;
t228 = mrSges(5,2) * t250 + qJD(3) * mrSges(5,3);
t245 = -m(5) * t193 + qJDD(3) * mrSges(5,1) + qJD(3) * t228;
t187 = m(4) * t194 + qJDD(3) * mrSges(4,1) + qJD(3) * t227 - t257 * t221 + (-t218 - t219) * t251 + t245;
t247 = t239 * t186 - t187 * t237;
t199 = t216 * t236 - t235 * t220;
t197 = -qJDD(1) * pkin(2) - pkin(5) * t242 - t199;
t190 = -pkin(3) * t222 - qJ(4) * t221 + (-0.2e1 * qJD(4) * t237 + (pkin(3) * t237 - qJ(4) * t239) * qJD(3)) * qJD(1) + t197;
t188 = m(5) * t190 - mrSges(5,1) * t222 - t221 * mrSges(5,3) - t226 * t251 - t228 * t250;
t243 = -m(4) * t197 + t222 * mrSges(4,1) - mrSges(4,2) * t221 - t225 * t251 + t227 * t250 - t188;
t189 = mrSges(5,2) * t221 + t218 * t251 - t245;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t248 - mrSges(2,2) * t244 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t199 - mrSges(3,2) * t200 + t237 * (mrSges(4,2) * t197 + mrSges(5,2) * t193 - mrSges(4,3) * t194 - mrSges(5,3) * t190 - qJ(4) * t188) + t239 * (-mrSges(4,1) * t197 - mrSges(5,1) * t190 + mrSges(5,2) * t192 + mrSges(4,3) * t195 - pkin(3) * t188) + pkin(2) * t243 + pkin(5) * t247 + pkin(1) * (t235 * (m(3) * t200 - mrSges(3,1) * t242 - qJDD(1) * mrSges(3,2) + t247) + t236 * (m(3) * t199 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t242 + t243)) + (t237 * t256 + t239 * t255) * qJDD(3) + t258 * qJD(3) + t262 * t222 + t263 * t221; m(3) * t234 + t186 * t237 + t187 * t239; mrSges(4,1) * t194 - mrSges(4,2) * t195 - mrSges(5,1) * t193 + mrSges(5,3) * t192 - pkin(3) * t189 + qJ(4) * t246 + (mrSges(5,2) * qJ(4) + t255) * t222 + t256 * t221 + (Ifges(4,3) + Ifges(5,2)) * qJDD(3) - t258 * qJD(1); t189;];
tauJ = t1;
