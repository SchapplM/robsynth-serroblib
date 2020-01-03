% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RPRP7
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RPRP7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP7_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP7_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_invdynJ_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP7_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP7_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP7_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:05
% EndTime: 2019-12-31 16:47:06
% DurationCPUTime: 0.50s
% Computational Cost: add. (545->118), mult. (998->146), div. (0->0), fcn. (382->4), ass. (0->49)
t226 = sin(qJ(3));
t228 = cos(qJ(3));
t249 = Ifges(4,4) - Ifges(5,5);
t258 = t228 * (Ifges(4,1) + Ifges(5,1)) - t226 * t249;
t257 = (Ifges(4,2) + Ifges(5,3)) * t226 - t228 * t249;
t248 = Ifges(5,4) + Ifges(4,5);
t247 = Ifges(4,6) - Ifges(5,6);
t254 = t228 * (t257 * qJD(1) - t247 * qJD(3));
t227 = sin(qJ(1));
t229 = cos(qJ(1));
t235 = -t229 * g(1) - t227 * g(2);
t253 = -qJDD(1) * qJ(2) - 0.2e1 * qJD(2) * qJD(1) - t235;
t252 = -pkin(1) - pkin(5);
t251 = t226 * g(3);
t250 = -mrSges(4,3) - mrSges(5,2);
t245 = t258 * qJD(1) + t248 * qJD(3);
t244 = qJD(1) * t226;
t243 = qJD(1) * t228;
t242 = qJD(1) * qJD(3);
t231 = qJD(1) ^ 2;
t198 = t252 * t231 - t253;
t213 = t226 * qJDD(1) + t228 * t242;
t214 = t228 * qJDD(1) - t226 * t242;
t191 = t213 * pkin(3) - t214 * qJ(4) + (-0.2e1 * qJD(4) * t228 + (pkin(3) * t228 + qJ(4) * t226) * qJD(3)) * qJD(1) + t198;
t220 = -mrSges(5,2) * t244 + qJD(3) * mrSges(5,3);
t240 = m(5) * t191 + t213 * mrSges(5,1) + t220 * t244;
t238 = t227 * g(1) - t229 * g(2);
t232 = -t231 * qJ(2) + qJDD(2) - t238;
t199 = t252 * qJDD(1) + t232;
t196 = -t228 * g(3) + t226 * t199;
t210 = (t226 * pkin(3) - t228 * qJ(4)) * qJD(1);
t230 = qJD(3) ^ 2;
t193 = -t230 * pkin(3) + qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) - t210 * t244 + t196;
t219 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t243;
t239 = m(5) * t193 + qJDD(3) * mrSges(5,3) + qJD(3) * t219;
t211 = (t226 * mrSges(5,1) - t228 * mrSges(5,3)) * qJD(1);
t237 = qJD(1) * (-t211 - (t226 * mrSges(4,1) + t228 * mrSges(4,2)) * qJD(1));
t194 = -qJDD(3) * pkin(3) - t251 - t230 * qJ(4) + qJDD(4) + (qJD(1) * t210 - t199) * t228;
t236 = -m(5) * t194 + qJDD(3) * mrSges(5,1) + qJD(3) * t220;
t195 = t228 * t199 + t251;
t217 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t244;
t218 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t243;
t234 = t226 * (m(4) * t196 - qJDD(3) * mrSges(4,2) - qJD(3) * t218 + t250 * t213 + t226 * t237 + t239) + t228 * (m(4) * t195 + qJDD(3) * mrSges(4,1) + qJD(3) * t217 + t250 * t214 + t228 * t237 + t236);
t201 = -qJDD(1) * pkin(1) + t232;
t200 = t231 * pkin(1) + t253;
t189 = t214 * mrSges(5,2) + t211 * t243 - t236;
t188 = -t214 * mrSges(5,3) - t219 * t243 + t240;
t185 = m(3) * t201 + qJDD(1) * mrSges(3,2) - t231 * mrSges(3,3) + t234;
t1 = [mrSges(2,1) * t238 - mrSges(2,2) * t235 + mrSges(3,2) * t201 - mrSges(3,3) * t200 + t228 * (mrSges(4,2) * t198 + mrSges(5,2) * t194 - mrSges(4,3) * t195 - mrSges(5,3) * t191 - qJ(4) * t188) - t226 * (-mrSges(4,1) * t198 - mrSges(5,1) * t191 + mrSges(5,2) * t193 + mrSges(4,3) * t196 - pkin(3) * t188) - pkin(5) * t234 - pkin(1) * t185 + t257 * t213 + (-t226 * t247 + t228 * t248) * qJDD(3) + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-t226 * t245 + t254) * qJD(3) + t258 * t214 + (-m(3) * t200 + m(4) * t198 + t231 * mrSges(3,2) + t240 + t213 * mrSges(4,1) + qJDD(1) * mrSges(3,3) + (mrSges(4,2) - mrSges(5,3)) * t214 + (t217 * t226 + (t218 - t219) * t228) * qJD(1)) * qJ(2); t185; mrSges(4,1) * t195 - mrSges(4,2) * t196 - mrSges(5,1) * t194 + mrSges(5,3) * t193 - pkin(3) * t189 + qJ(4) * t239 + t248 * t214 + (-qJ(4) * mrSges(5,2) - t247) * t213 + (Ifges(4,3) + Ifges(5,2)) * qJDD(3) + (-t254 + (-qJ(4) * t211 + t245) * t226) * qJD(1); t189;];
tauJ = t1;
