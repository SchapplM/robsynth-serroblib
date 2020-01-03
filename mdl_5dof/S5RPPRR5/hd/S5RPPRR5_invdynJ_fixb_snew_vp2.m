% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPPRR5
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPRR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR5_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR5_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR5_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR5_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR5_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:33
% EndTime: 2019-12-31 17:56:33
% DurationCPUTime: 0.30s
% Computational Cost: add. (1856->112), mult. (2597->143), div. (0->0), fcn. (1100->8), ass. (0->54)
t252 = -pkin(2) - pkin(3);
t227 = -qJD(1) + qJD(4);
t235 = sin(qJ(5));
t251 = t227 * t235;
t238 = cos(qJ(5));
t250 = t227 * t238;
t241 = qJD(1) ^ 2;
t237 = sin(qJ(1));
t240 = cos(qJ(1));
t248 = t237 * g(1) - t240 * g(2);
t219 = qJDD(1) * pkin(1) + t248;
t246 = -t240 * g(1) - t237 * g(2);
t220 = -t241 * pkin(1) + t246;
t233 = sin(pkin(8));
t234 = cos(pkin(8));
t207 = t233 * t219 + t234 * t220;
t247 = qJDD(1) * qJ(3) + (2 * qJD(3) * qJD(1)) + t207;
t201 = t252 * t241 + t247;
t206 = t234 * t219 - t233 * t220;
t244 = -t241 * qJ(3) + qJDD(3) - t206;
t203 = t252 * qJDD(1) + t244;
t236 = sin(qJ(4));
t239 = cos(qJ(4));
t198 = t239 * t201 + t236 * t203;
t249 = qJD(5) * t227;
t225 = t227 ^ 2;
t226 = -qJDD(1) + qJDD(4);
t196 = -t225 * pkin(4) + t226 * pkin(7) + t198;
t231 = g(3) - qJDD(2);
t193 = -t235 * t196 + t238 * t231;
t213 = (-mrSges(6,1) * t238 + mrSges(6,2) * t235) * t227;
t214 = t235 * t226 + t238 * t249;
t222 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t250;
t191 = m(6) * t193 + qJDD(5) * mrSges(6,1) - t214 * mrSges(6,3) + qJD(5) * t222 - t213 * t251;
t194 = t238 * t196 + t235 * t231;
t215 = t238 * t226 - t235 * t249;
t221 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t251;
t192 = m(6) * t194 - qJDD(5) * mrSges(6,2) + t215 * mrSges(6,3) - qJD(5) * t221 + t213 * t250;
t185 = -t235 * t191 + t238 * t192;
t184 = m(5) * t198 - t225 * mrSges(5,1) - t226 * mrSges(5,2) + t185;
t197 = -t236 * t201 + t239 * t203;
t195 = -t226 * pkin(4) - t225 * pkin(7) - t197;
t189 = -m(6) * t195 + t215 * mrSges(6,1) - t214 * mrSges(6,2) - t221 * t251 + t222 * t250;
t188 = m(5) * t197 + t226 * mrSges(5,1) - t225 * mrSges(5,2) + t189;
t245 = t236 * t184 + t239 * t188;
t204 = -t241 * pkin(2) + t247;
t243 = m(4) * t204 - t241 * mrSges(4,1) + qJDD(1) * mrSges(4,3) + t239 * t184 - t236 * t188;
t205 = -qJDD(1) * pkin(2) + t244;
t182 = m(4) * t205 - qJDD(1) * mrSges(4,1) - t241 * mrSges(4,3) + t245;
t208 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t235 + Ifges(6,6) * t238) * t227;
t209 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t235 + Ifges(6,2) * t238) * t227;
t210 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t235 + Ifges(6,4) * t238) * t227;
t242 = mrSges(5,1) * t197 - mrSges(5,2) * t198 + Ifges(5,3) * t226 + pkin(4) * t189 + pkin(7) * t185 + t238 * (-mrSges(6,1) * t195 + mrSges(6,3) * t194 + Ifges(6,4) * t214 + Ifges(6,2) * t215 + Ifges(6,6) * qJDD(5) + qJD(5) * t210 - t208 * t251) + t235 * (mrSges(6,2) * t195 - mrSges(6,3) * t193 + Ifges(6,1) * t214 + Ifges(6,4) * t215 + Ifges(6,5) * qJDD(5) - qJD(5) * t209 + t208 * t250);
t1 = [-t242 + pkin(1) * (t233 * (m(3) * t207 - t241 * mrSges(3,1) + t243) + t234 * (m(3) * t206 - t241 * mrSges(3,2) - t182)) + (pkin(1) * (t234 * mrSges(3,1) - t233 * mrSges(3,2)) + Ifges(2,3) + Ifges(3,3) + Ifges(4,2)) * qJDD(1) + qJ(3) * t243 - pkin(3) * t245 + mrSges(2,1) * t248 - mrSges(2,2) * t246 + mrSges(4,3) * t204 - mrSges(4,1) * t205 + mrSges(3,1) * t206 - mrSges(3,2) * t207 - pkin(2) * t182; -t238 * t191 - t235 * t192 + (-m(3) - m(4) - m(5)) * t231; t182; t242; mrSges(6,1) * t193 - mrSges(6,2) * t194 + Ifges(6,5) * t214 + Ifges(6,6) * t215 + Ifges(6,3) * qJDD(5) + (t209 * t235 - t210 * t238) * t227;];
tauJ = t1;
