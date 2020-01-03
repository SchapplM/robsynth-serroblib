% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4PPRR4
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
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4PPRR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR4_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR4_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR4_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR4_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:31
% EndTime: 2019-12-31 16:18:32
% DurationCPUTime: 0.51s
% Computational Cost: add. (4274->131), mult. (7032->175), div. (0->0), fcn. (4448->8), ass. (0->59)
t234 = sin(pkin(6));
t236 = cos(pkin(6));
t229 = -t236 * g(1) - t234 * g(2);
t232 = -g(3) + qJDD(1);
t233 = sin(pkin(7));
t235 = cos(pkin(7));
t217 = -t233 * t229 + t235 * t232;
t218 = t235 * t229 + t233 * t232;
t238 = sin(qJ(3));
t240 = cos(qJ(3));
t214 = t238 * t217 + t240 * t218;
t241 = qJD(3) ^ 2;
t212 = -t241 * pkin(3) + qJDD(3) * pkin(5) + t214;
t228 = t234 * g(1) - t236 * g(2);
t227 = qJDD(2) - t228;
t237 = sin(qJ(4));
t239 = cos(qJ(4));
t209 = -t237 * t212 + t239 * t227;
t224 = (-mrSges(5,1) * t239 + mrSges(5,2) * t237) * qJD(3);
t248 = qJD(3) * qJD(4);
t225 = t237 * qJDD(3) + t239 * t248;
t249 = qJD(3) * t239;
t231 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t249;
t250 = qJD(3) * t237;
t207 = m(5) * t209 + qJDD(4) * mrSges(5,1) - t225 * mrSges(5,3) + qJD(4) * t231 - t224 * t250;
t210 = t239 * t212 + t237 * t227;
t226 = t239 * qJDD(3) - t237 * t248;
t230 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t250;
t208 = m(5) * t210 - qJDD(4) * mrSges(5,2) + t226 * mrSges(5,3) - qJD(4) * t230 + t224 * t249;
t244 = -t237 * t207 + t239 * t208;
t196 = m(4) * t214 - t241 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t244;
t213 = t240 * t217 - t238 * t218;
t211 = -qJDD(3) * pkin(3) - t241 * pkin(5) - t213;
t242 = -m(5) * t211 + t226 * mrSges(5,1) - t225 * mrSges(5,2) - t230 * t250 + t231 * t249;
t203 = m(4) * t213 + qJDD(3) * mrSges(4,1) - t241 * mrSges(4,2) + t242;
t193 = t238 * t196 + t240 * t203;
t191 = m(3) * t217 + t193;
t245 = t240 * t196 - t238 * t203;
t192 = m(3) * t218 + t245;
t246 = -t233 * t191 + t235 * t192;
t184 = m(2) * t229 + t246;
t199 = t239 * t207 + t237 * t208;
t243 = (-m(3) - m(4)) * t227 - t199;
t198 = m(2) * t228 + t243;
t251 = t234 * t184 + t236 * t198;
t185 = t235 * t191 + t233 * t192;
t247 = t236 * t184 - t234 * t198;
t221 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t237 + Ifges(5,4) * t239) * qJD(3);
t220 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t237 + Ifges(5,2) * t239) * qJD(3);
t219 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t237 + Ifges(5,6) * t239) * qJD(3);
t201 = mrSges(5,2) * t211 - mrSges(5,3) * t209 + Ifges(5,1) * t225 + Ifges(5,4) * t226 + Ifges(5,5) * qJDD(4) - qJD(4) * t220 + t219 * t249;
t200 = -mrSges(5,1) * t211 + mrSges(5,3) * t210 + Ifges(5,4) * t225 + Ifges(5,2) * t226 + Ifges(5,6) * qJDD(4) + qJD(4) * t221 - t219 * t250;
t187 = -mrSges(4,1) * t227 - mrSges(5,1) * t209 + mrSges(5,2) * t210 + mrSges(4,3) * t214 + t241 * Ifges(4,5) - Ifges(5,5) * t225 + Ifges(4,6) * qJDD(3) - Ifges(5,6) * t226 - Ifges(5,3) * qJDD(4) - pkin(3) * t199 + (-t220 * t237 + t221 * t239) * qJD(3);
t186 = mrSges(4,2) * t227 - mrSges(4,3) * t213 + Ifges(4,5) * qJDD(3) - t241 * Ifges(4,6) - pkin(5) * t199 - t237 * t200 + t239 * t201;
t181 = mrSges(3,2) * t227 - mrSges(3,3) * t217 - pkin(4) * t193 + t240 * t186 - t238 * t187;
t180 = -mrSges(3,1) * t227 + mrSges(3,3) * t218 + t238 * t186 + t240 * t187 - pkin(2) * (m(4) * t227 + t199) + pkin(4) * t245;
t179 = -mrSges(2,1) * t232 - mrSges(3,1) * t217 - mrSges(4,1) * t213 + mrSges(3,2) * t218 + mrSges(4,2) * t214 + mrSges(2,3) * t229 - Ifges(4,3) * qJDD(3) - pkin(1) * t185 - pkin(2) * t193 - pkin(3) * t242 - pkin(5) * t244 - t239 * t200 - t237 * t201;
t178 = mrSges(2,2) * t232 - mrSges(2,3) * t228 - qJ(2) * t185 - t233 * t180 + t235 * t181;
t1 = [-m(1) * g(1) + t247; -m(1) * g(2) + t251; -m(1) * g(3) + m(2) * t232 + t185; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t251 + t236 * t178 - t234 * t179; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t247 + t234 * t178 + t236 * t179; -mrSges(1,1) * g(2) + mrSges(2,1) * t228 + mrSges(1,2) * g(1) - mrSges(2,2) * t229 + pkin(1) * t243 + qJ(2) * t246 + t235 * t180 + t233 * t181;];
tauB = t1;
