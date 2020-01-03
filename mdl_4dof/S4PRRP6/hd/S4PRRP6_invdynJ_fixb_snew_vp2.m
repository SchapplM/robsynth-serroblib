% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4PRRP6
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
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4PRRP6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP6_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP6_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP6_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP6_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:16
% EndTime: 2019-12-31 16:30:16
% DurationCPUTime: 0.41s
% Computational Cost: add. (519->114), mult. (968->145), div. (0->0), fcn. (461->6), ass. (0->49)
t224 = sin(qJ(3));
t226 = cos(qJ(3));
t244 = Ifges(4,4) - Ifges(5,5);
t251 = t224 * (Ifges(4,1) + Ifges(5,1)) + t226 * t244;
t250 = -t224 * t244 - t226 * (Ifges(4,2) + Ifges(5,3));
t243 = Ifges(5,4) + Ifges(4,5);
t248 = Ifges(5,6) - Ifges(4,6);
t246 = t224 * (t250 * qJD(2) + t248 * qJD(3)) + t226 * (t251 * qJD(2) + t243 * qJD(3));
t245 = mrSges(4,3) + mrSges(5,2);
t235 = qJD(2) * qJD(3);
t208 = t226 * qJDD(2) - t224 * t235;
t241 = t208 * mrSges(5,1);
t222 = sin(pkin(6));
t223 = cos(pkin(6));
t211 = -t222 * g(1) + t223 * g(2);
t240 = t226 * t211;
t212 = -t223 * g(1) - t222 * g(2);
t221 = -g(3) + qJDD(1);
t225 = sin(qJ(2));
t227 = cos(qJ(2));
t188 = t227 * t212 + t225 * t221;
t229 = qJD(2) ^ 2;
t186 = -t229 * pkin(2) + qJDD(2) * pkin(5) + t188;
t183 = t226 * t186 + t224 * t211;
t237 = qJD(2) * t224;
t236 = qJD(2) * t226;
t182 = -t224 * t186 + t240;
t205 = (-t226 * mrSges(5,1) - t224 * mrSges(5,3)) * qJD(2);
t206 = (-t226 * mrSges(4,1) + t224 * mrSges(4,2)) * qJD(2);
t207 = t224 * qJDD(2) + t226 * t235;
t213 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t237;
t215 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t236;
t204 = (-t226 * pkin(3) - t224 * qJ(4)) * qJD(2);
t228 = qJD(3) ^ 2;
t181 = -qJDD(3) * pkin(3) - t228 * qJ(4) - t240 + qJDD(4) + (qJD(2) * t204 + t186) * t224;
t216 = mrSges(5,2) * t236 + qJD(3) * mrSges(5,3);
t232 = -m(5) * t181 + qJDD(3) * mrSges(5,1) + qJD(3) * t216;
t180 = -t228 * pkin(3) + qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) + t204 * t236 + t183;
t214 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t237;
t233 = m(5) * t180 + qJDD(3) * mrSges(5,3) + qJD(3) * t214 + t205 * t236;
t234 = -t224 * (m(4) * t182 + qJDD(3) * mrSges(4,1) + qJD(3) * t215 - t245 * t207 + (-t205 - t206) * t237 + t232) + t226 * (m(4) * t183 - qJDD(3) * mrSges(4,2) - qJD(3) * t213 + t206 * t236 + t245 * t208 + t233);
t187 = -t225 * t212 + t227 * t221;
t185 = -qJDD(2) * pkin(2) - t229 * pkin(5) - t187;
t178 = -t208 * pkin(3) - t207 * qJ(4) + (-0.2e1 * qJD(4) * t224 + (pkin(3) * t224 - qJ(4) * t226) * qJD(3)) * qJD(2) + t185;
t231 = m(5) * t178 - t207 * mrSges(5,3) - t214 * t237 - t216 * t236;
t230 = -m(4) * t185 + t208 * mrSges(4,1) - t213 * t237 + t215 * t236 - t231;
t177 = t207 * mrSges(5,2) + t205 * t237 - t232;
t176 = t231 - t241;
t1 = [m(2) * t221 + t225 * (m(3) * t188 - t229 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t234) + t227 * (m(3) * t187 + qJDD(2) * mrSges(3,1) - t229 * mrSges(3,2) - t207 * mrSges(4,2) + t230 + t241); Ifges(3,3) * qJDD(2) + mrSges(3,1) * t187 - mrSges(3,2) * t188 + t224 * (mrSges(4,2) * t185 + mrSges(5,2) * t181 - mrSges(4,3) * t182 - mrSges(5,3) * t178 - qJ(4) * t176) + t226 * (-mrSges(4,1) * t185 - mrSges(5,1) * t178 + mrSges(5,2) * t180 + mrSges(4,3) * t183 - pkin(3) * t176) + pkin(2) * t230 + pkin(5) * t234 + (pkin(2) * mrSges(5,1) - t250) * t208 + (-pkin(2) * mrSges(4,2) + t251) * t207 + (t224 * t243 - t226 * t248) * qJDD(3) + t246 * qJD(3); mrSges(4,1) * t182 - mrSges(4,2) * t183 - mrSges(5,1) * t181 + mrSges(5,3) * t180 - pkin(3) * t177 + qJ(4) * t233 + (qJ(4) * mrSges(5,2) - t248) * t208 + t243 * t207 + (Ifges(4,3) + Ifges(5,2)) * qJDD(3) - t246 * qJD(2); t177;];
tauJ = t1;
