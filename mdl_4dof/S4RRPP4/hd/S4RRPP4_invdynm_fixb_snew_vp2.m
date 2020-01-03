% Calculate vector of cutting torques with Newton-Euler for
% S4RRPP4
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
%   pkin=[a2,a3,a4,d1,d2]';
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
% m [3x5]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RRPP4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP4_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_invdynm_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP4_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP4_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP4_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:58:54
% EndTime: 2019-12-31 16:58:56
% DurationCPUTime: 1.06s
% Computational Cost: add. (4660->235), mult. (9566->281), div. (0->0), fcn. (3820->4), ass. (0->82)
t207 = sin(qJ(1));
t209 = cos(qJ(1));
t190 = -g(1) * t209 - g(2) * t207;
t211 = qJD(1) ^ 2;
t149 = -pkin(1) * t211 + qJDD(1) * pkin(5) + t190;
t206 = sin(qJ(2));
t208 = cos(qJ(2));
t126 = -t208 * g(3) - t206 * t149;
t127 = -g(3) * t206 + t208 * t149;
t140 = Ifges(4,6) * qJD(2) + (Ifges(4,5) * t206 - Ifges(4,3) * t208) * qJD(1);
t144 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t206 + Ifges(3,2) * t208) * qJD(1);
t172 = (-mrSges(4,1) * t208 - mrSges(4,3) * t206) * qJD(1);
t231 = qJD(1) * qJD(2);
t227 = t208 * t231;
t175 = qJDD(1) * t206 + t227;
t176 = qJDD(1) * t208 - t206 * t231;
t171 = (-pkin(2) * t208 - qJ(3) * t206) * qJD(1);
t210 = qJD(2) ^ 2;
t233 = qJD(1) * t206;
t125 = -qJDD(2) * pkin(2) - qJ(3) * t210 + t171 * t233 + qJDD(3) - t126;
t230 = -0.2e1 * qJD(1) * qJD(4);
t119 = t206 * t230 + (-t175 + t227) * qJ(4) + (-t206 * t208 * t211 - qJDD(2)) * pkin(3) + t125;
t173 = (mrSges(5,1) * t208 + mrSges(5,2) * t206) * qJD(1);
t232 = qJD(1) * t208;
t186 = qJD(2) * mrSges(5,2) - mrSges(5,3) * t232;
t110 = m(5) * t119 - qJDD(2) * mrSges(5,1) - t175 * mrSges(5,3) - qJD(2) * t186 - t173 * t233;
t242 = 2 * qJD(3);
t123 = -pkin(2) * t210 + qJDD(2) * qJ(3) + qJD(2) * t242 + t171 * t232 + t127;
t182 = -qJD(2) * pkin(3) - qJ(4) * t233;
t205 = t208 ^ 2;
t118 = -pkin(3) * t205 * t211 - qJ(4) * t176 + qJD(2) * t182 + t208 * t230 + t123;
t142 = -Ifges(5,6) * qJD(2) + (Ifges(5,4) * t206 - Ifges(5,2) * t208) * qJD(1);
t145 = -Ifges(5,5) * qJD(2) + (Ifges(5,1) * t206 - Ifges(5,4) * t208) * qJD(1);
t221 = mrSges(5,1) * t119 - mrSges(5,2) * t118 + Ifges(5,5) * t175 - Ifges(5,6) * t176 - Ifges(5,3) * qJDD(2) + t142 * t233 + t145 * t232;
t215 = -mrSges(4,1) * t125 + mrSges(4,3) * t123 + Ifges(4,4) * t175 + Ifges(4,2) * qJDD(2) - Ifges(4,6) * t176 - pkin(3) * t110 - t221;
t188 = mrSges(4,2) * t232 + qJD(2) * mrSges(4,3);
t217 = -m(4) * t125 + qJDD(2) * mrSges(4,1) + qJD(2) * t188 - t110;
t185 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t233;
t183 = -qJD(2) * mrSges(5,1) - mrSges(5,3) * t233;
t225 = m(5) * t118 + qJDD(2) * mrSges(5,2) - t176 * mrSges(5,3) + qJD(2) * t183;
t219 = m(4) * t123 + qJDD(2) * mrSges(4,3) + qJD(2) * t185 + t172 * t232 + t225;
t229 = t173 * t232;
t146 = Ifges(4,4) * qJD(2) + (Ifges(4,1) * t206 - Ifges(4,5) * t208) * qJD(1);
t235 = t146 + Ifges(3,5) * qJD(2) + (Ifges(3,1) * t206 + Ifges(3,4) * t208) * qJD(1);
t245 = -((t140 - t144) * t206 + t235 * t208) * qJD(1) + mrSges(3,1) * t126 - mrSges(3,2) * t127 + Ifges(3,5) * t175 + Ifges(3,6) * t176 + Ifges(3,3) * qJDD(2) + pkin(2) * (-mrSges(4,2) * t175 - t172 * t233 + t217) + qJ(3) * (mrSges(4,2) * t176 + t219 - t229) + t215;
t139 = -Ifges(5,3) * qJD(2) + (Ifges(5,5) * t206 - Ifges(5,6) * t208) * qJD(1);
t244 = -Ifges(5,5) * qJDD(2) - t139 * t232;
t241 = pkin(5) * t211;
t240 = mrSges(3,3) + mrSges(4,2);
t239 = Ifges(4,6) - Ifges(5,6);
t238 = qJ(3) * t208;
t143 = Ifges(4,2) * qJD(2) + (Ifges(4,4) * t206 - Ifges(4,6) * t208) * qJD(1);
t237 = -t139 + t143;
t189 = t207 * g(1) - t209 * g(2);
t174 = (-mrSges(3,1) * t208 + mrSges(3,2) * t206) * qJD(1);
t184 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t233;
t103 = m(3) * t127 - qJDD(2) * mrSges(3,2) - qJD(2) * t184 + t240 * t176 + (-t173 + t174) * t232 + t219;
t187 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t232;
t104 = m(3) * t126 + qJDD(2) * mrSges(3,1) + qJD(2) * t187 - t240 * t175 + (-t172 - t174) * t233 + t217;
t226 = t208 * t103 - t104 * t206;
t224 = qJDD(1) * pkin(1) + t189;
t220 = -qJ(3) * t175 - t224;
t113 = qJDD(4) + (-qJ(4) * t205 + pkin(5)) * t211 + (pkin(2) + pkin(3)) * t176 + (qJD(2) * t238 + (-pkin(2) * qJD(2) + t182 + t242) * t206) * qJD(1) - t220;
t108 = -m(5) * t113 - t176 * mrSges(5,1) - t175 * mrSges(5,2) - t183 * t233 - t186 * t232;
t223 = mrSges(5,1) * t113 - mrSges(5,3) * t118 - Ifges(5,4) * t175 + Ifges(5,2) * t176;
t222 = mrSges(5,2) * t113 - mrSges(5,3) * t119 + Ifges(5,1) * t175 - Ifges(5,4) * t176 + qJD(2) * t142;
t120 = -pkin(2) * t176 - t241 + (-0.2e1 * qJD(3) * t206 + (pkin(2) * t206 - t238) * qJD(2)) * qJD(1) + t220;
t105 = m(4) * t120 - mrSges(4,1) * t176 - t175 * mrSges(4,3) - t185 * t233 - t188 * t232 + t108;
t148 = -t224 - t241;
t213 = -m(3) * t148 + t176 * mrSges(3,1) - mrSges(3,2) * t175 - t184 * t233 + t187 * t232 - t105;
t141 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t206 + Ifges(3,6) * t208) * qJD(1);
t216 = mrSges(4,1) * t120 - mrSges(4,2) * t123 + pkin(3) * t108 + qJ(4) * (t225 - t229) - t223;
t94 = -t216 + (t145 + t235) * qJD(2) + (Ifges(3,6) - t239) * qJDD(2) + (Ifges(3,4) - Ifges(4,5)) * t175 + (Ifges(3,2) + Ifges(4,3)) * t176 + (-t141 - t237) * t233 - mrSges(3,1) * t148 + mrSges(3,3) * t127 - pkin(2) * t105;
t214 = mrSges(4,2) * t125 - mrSges(4,3) * t120 + Ifges(4,1) * t175 + Ifges(4,4) * qJDD(2) - Ifges(4,5) * t176 - qJ(4) * t110 + qJD(2) * t140 + t143 * t232 + t222;
t96 = t214 + (-t139 + t141) * t232 + (Ifges(3,5) - Ifges(5,5)) * qJDD(2) + Ifges(3,1) * t175 + Ifges(3,4) * t176 - qJD(2) * t144 + mrSges(3,2) * t148 - mrSges(3,3) * t126 - qJ(3) * t105;
t218 = mrSges(2,1) * t189 - mrSges(2,2) * t190 + Ifges(2,3) * qJDD(1) + pkin(1) * t213 + pkin(5) * t226 + t206 * t96 + t208 * t94;
t100 = m(2) * t189 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t211 + t213;
t99 = t103 * t206 + t104 * t208;
t97 = m(2) * t190 - mrSges(2,1) * t211 - qJDD(1) * mrSges(2,2) + t226;
t92 = mrSges(2,1) * g(3) + mrSges(2,3) * t190 + t211 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t99 - t245;
t91 = -mrSges(2,2) * g(3) - mrSges(2,3) * t189 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t211 - pkin(5) * t99 - t206 * t94 + t208 * t96;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t209 * t91 - t207 * t92 - pkin(4) * (t100 * t209 + t207 * t97), t91, t96, t214 + t244, t222 + t244; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t207 * t91 + t209 * t92 + pkin(4) * (-t100 * t207 + t209 * t97), t92, t94, t215 + (-t206 * t140 - t208 * t146) * qJD(1), -Ifges(5,6) * qJDD(2) - qJD(2) * t145 - t139 * t233 - t223; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t218, t218, t245, Ifges(4,5) * t175 - Ifges(4,3) * t176 + t239 * qJDD(2) + (-t145 - t146) * qJD(2) + t237 * t233 + t216, t221;];
m_new = t1;
