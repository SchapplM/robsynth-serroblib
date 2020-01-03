% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RPRR4
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RPRR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR4_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR4_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR4_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR4_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR4_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:19
% EndTime: 2019-12-31 16:50:19
% DurationCPUTime: 0.36s
% Computational Cost: add. (1569->145), mult. (2935->190), div. (0->0), fcn. (1618->8), ass. (0->65)
t229 = -g(3) + qJDD(2);
t236 = cos(qJ(3));
t252 = t236 * t229;
t234 = sin(qJ(1));
t237 = cos(qJ(1));
t246 = t234 * g(1) - t237 * g(2);
t218 = qJDD(1) * pkin(1) + t246;
t239 = qJD(1) ^ 2;
t243 = -t237 * g(1) - t234 * g(2);
t220 = -t239 * pkin(1) + t243;
t230 = sin(pkin(7));
t231 = cos(pkin(7));
t203 = t230 * t218 + t231 * t220;
t194 = -t239 * pkin(2) + qJDD(1) * pkin(5) + t203;
t233 = sin(qJ(3));
t191 = t236 * t194 + t233 * t229;
t251 = qJD(1) * t233;
t250 = t236 * qJD(1);
t249 = qJD(1) * qJD(3);
t248 = t233 * t249;
t247 = t236 * t249;
t219 = (-t236 * mrSges(4,1) + t233 * mrSges(4,2)) * qJD(1);
t223 = t236 * qJDD(1) - t248;
t224 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t251;
t202 = t231 * t218 - t230 * t220;
t193 = -qJDD(1) * pkin(2) - t239 * pkin(5) - t202;
t222 = t233 * qJDD(1) + t247;
t187 = (-t222 - t247) * pkin(6) + (-t223 + t248) * pkin(3) + t193;
t221 = (-t236 * pkin(3) - t233 * pkin(6)) * qJD(1);
t238 = qJD(3) ^ 2;
t189 = -t238 * pkin(3) + qJDD(3) * pkin(6) + t221 * t250 + t191;
t232 = sin(qJ(4));
t235 = cos(qJ(4));
t185 = t235 * t187 - t232 * t189;
t216 = t235 * qJD(3) - t232 * t251;
t201 = t216 * qJD(4) + t232 * qJDD(3) + t235 * t222;
t217 = t232 * qJD(3) + t235 * t251;
t204 = -t216 * mrSges(5,1) + t217 * mrSges(5,2);
t226 = qJD(4) - t250;
t205 = -t226 * mrSges(5,2) + t216 * mrSges(5,3);
t215 = qJDD(4) - t223;
t183 = m(5) * t185 + t215 * mrSges(5,1) - t201 * mrSges(5,3) - t217 * t204 + t226 * t205;
t186 = t232 * t187 + t235 * t189;
t200 = -t217 * qJD(4) + t235 * qJDD(3) - t232 * t222;
t206 = t226 * mrSges(5,1) - t217 * mrSges(5,3);
t184 = m(5) * t186 - t215 * mrSges(5,2) + t200 * mrSges(5,3) + t216 * t204 - t226 * t206;
t244 = -t232 * t183 + t235 * t184;
t177 = m(4) * t191 - qJDD(3) * mrSges(4,2) + t223 * mrSges(4,3) - qJD(3) * t224 + t219 * t250 + t244;
t190 = -t233 * t194 + t252;
t225 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t250;
t188 = -qJDD(3) * pkin(3) - t238 * pkin(6) - t252 + (qJD(1) * t221 + t194) * t233;
t242 = -m(5) * t188 + t200 * mrSges(5,1) - t201 * mrSges(5,2) + t216 * t205 - t217 * t206;
t181 = m(4) * t190 + qJDD(3) * mrSges(4,1) - t222 * mrSges(4,3) + qJD(3) * t225 - t219 * t251 + t242;
t245 = t236 * t177 - t233 * t181;
t178 = t235 * t183 + t232 * t184;
t241 = -m(4) * t193 + t223 * mrSges(4,1) - t222 * mrSges(4,2) - t224 * t251 + t225 * t250 - t178;
t196 = Ifges(5,4) * t217 + Ifges(5,2) * t216 + Ifges(5,6) * t226;
t197 = Ifges(5,1) * t217 + Ifges(5,4) * t216 + Ifges(5,5) * t226;
t240 = mrSges(5,1) * t185 - mrSges(5,2) * t186 + Ifges(5,5) * t201 + Ifges(5,6) * t200 + Ifges(5,3) * t215 + t217 * t196 - t216 * t197;
t212 = Ifges(4,5) * qJD(3) + (t233 * Ifges(4,1) + t236 * Ifges(4,4)) * qJD(1);
t211 = Ifges(4,6) * qJD(3) + (t233 * Ifges(4,4) + t236 * Ifges(4,2)) * qJD(1);
t195 = Ifges(5,5) * t217 + Ifges(5,6) * t216 + Ifges(5,3) * t226;
t180 = mrSges(5,2) * t188 - mrSges(5,3) * t185 + Ifges(5,1) * t201 + Ifges(5,4) * t200 + Ifges(5,5) * t215 + t216 * t195 - t226 * t196;
t179 = -mrSges(5,1) * t188 + mrSges(5,3) * t186 + Ifges(5,4) * t201 + Ifges(5,2) * t200 + Ifges(5,6) * t215 - t217 * t195 + t226 * t197;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t246 - mrSges(2,2) * t243 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t202 - mrSges(3,2) * t203 + t233 * (mrSges(4,2) * t193 - mrSges(4,3) * t190 + Ifges(4,1) * t222 + Ifges(4,4) * t223 + Ifges(4,5) * qJDD(3) - pkin(6) * t178 - qJD(3) * t211 - t232 * t179 + t235 * t180) + t236 * (-mrSges(4,1) * t193 + mrSges(4,3) * t191 + Ifges(4,4) * t222 + Ifges(4,2) * t223 + Ifges(4,6) * qJDD(3) - pkin(3) * t178 + qJD(3) * t212 - t240) + pkin(2) * t241 + pkin(5) * t245 + pkin(1) * (t230 * (m(3) * t203 - t239 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t245) + t231 * (m(3) * t202 + qJDD(1) * mrSges(3,1) - t239 * mrSges(3,2) + t241)); m(3) * t229 + t233 * t177 + t236 * t181; Ifges(4,5) * t222 + Ifges(4,6) * t223 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t190 - mrSges(4,2) * t191 + t232 * t180 + t235 * t179 + pkin(3) * t242 + pkin(6) * t244 + (t233 * t211 - t236 * t212) * qJD(1); t240;];
tauJ = t1;
