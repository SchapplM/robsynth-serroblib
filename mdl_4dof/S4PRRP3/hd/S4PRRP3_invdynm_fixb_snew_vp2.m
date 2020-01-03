% Calculate vector of cutting torques with Newton-Euler for
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
% m [3x5]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PRRP3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP3_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP3_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP3_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP3_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:45
% EndTime: 2019-12-31 16:26:46
% DurationCPUTime: 0.79s
% Computational Cost: add. (5989->178), mult. (11625->224), div. (0->0), fcn. (5708->6), ass. (0->72)
t159 = sin(pkin(6));
t160 = cos(pkin(6));
t143 = t159 * g(1) - t160 * g(2);
t144 = -t160 * g(1) - t159 * g(2);
t162 = sin(qJ(2));
t164 = cos(qJ(2));
t112 = t162 * t143 + t164 * t144;
t165 = qJD(2) ^ 2;
t109 = -t165 * pkin(2) + qJDD(2) * pkin(5) + t112;
t158 = -g(3) + qJDD(1);
t163 = cos(qJ(3));
t151 = t163 * t158;
t161 = sin(qJ(3));
t105 = -t161 * t109 + t151;
t106 = t163 * t109 + t161 * t158;
t120 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t161 + Ifges(4,2) * t163) * qJD(2);
t121 = Ifges(5,5) * qJD(3) + (Ifges(5,1) * t161 + Ifges(5,4) * t163) * qJD(2);
t122 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t161 + Ifges(4,4) * t163) * qJD(2);
t181 = qJD(2) * qJD(3);
t177 = t163 * t181;
t139 = t161 * qJDD(2) + t177;
t140 = t163 * qJDD(2) - t161 * t181;
t180 = qJD(2) * qJD(4);
t188 = pkin(3) * t165;
t101 = qJDD(3) * pkin(3) + t151 + (-t139 + t177) * qJ(4) + (t163 * t188 - t109 - 0.2e1 * t180) * t161;
t183 = qJD(2) * t161;
t145 = qJD(3) * pkin(3) - qJ(4) * t183;
t157 = t163 ^ 2;
t102 = t140 * qJ(4) - qJD(3) * t145 - t157 * t188 + 0.2e1 * t163 * t180 + t106;
t119 = Ifges(5,6) * qJD(3) + (Ifges(5,4) * t161 + Ifges(5,2) * t163) * qJD(2);
t171 = -mrSges(5,1) * t101 + mrSges(5,2) * t102 - Ifges(5,5) * t139 - Ifges(5,6) * t140 - Ifges(5,3) * qJDD(3) - t119 * t183;
t137 = (-mrSges(5,1) * t163 + mrSges(5,2) * t161) * qJD(2);
t182 = qJD(2) * t163;
t148 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t182;
t179 = m(5) * t101 + qJDD(3) * mrSges(5,1) + qJD(3) * t148;
t96 = -t139 * mrSges(5,3) - t137 * t183 + t179;
t190 = mrSges(4,1) * t105 - mrSges(4,2) * t106 + Ifges(4,5) * t139 + Ifges(4,6) * t140 + Ifges(4,3) * qJDD(3) + pkin(3) * t96 - (-t161 * t120 + (t121 + t122) * t163) * qJD(2) - t171;
t187 = -mrSges(4,2) - mrSges(5,2);
t138 = (-mrSges(4,1) * t163 + mrSges(4,2) * t161) * qJD(2);
t149 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t182;
t93 = m(4) * t105 + qJDD(3) * mrSges(4,1) + qJD(3) * t149 + (-mrSges(4,3) - mrSges(5,3)) * t139 + (-t137 - t138) * t183 + t179;
t178 = m(5) * t102 + t140 * mrSges(5,3) + t137 * t182;
t146 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t183;
t184 = -qJD(3) * mrSges(4,1) + mrSges(4,3) * t183 - t146;
t94 = m(4) * t106 + t140 * mrSges(4,3) + t184 * qJD(3) + t187 * qJDD(3) + t138 * t182 + t178;
t176 = -t161 * t93 + t163 * t94;
t85 = m(3) * t112 - t165 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t176;
t111 = t164 * t143 - t162 * t144;
t173 = -qJDD(2) * pkin(2) - t111;
t108 = -t165 * pkin(5) + t173;
t104 = t145 * t183 - t140 * pkin(3) + qJDD(4) + (-qJ(4) * t157 - pkin(5)) * t165 + t173;
t174 = -m(5) * t104 + t140 * mrSges(5,1) + t148 * t182;
t167 = -m(4) * t108 + t140 * mrSges(4,1) + t187 * t139 + t149 * t182 + t184 * t183 + t174;
t89 = m(3) * t111 + qJDD(2) * mrSges(3,1) - t165 * mrSges(3,2) + t167;
t78 = t162 * t85 + t164 * t89;
t87 = t161 * t94 + t163 * t93;
t175 = -t162 * t89 + t164 * t85;
t172 = -mrSges(5,1) * t104 + mrSges(5,3) * t102 + Ifges(5,4) * t139 + Ifges(5,2) * t140 + Ifges(5,6) * qJDD(3) + qJD(3) * t121;
t117 = Ifges(5,3) * qJD(3) + (Ifges(5,5) * t161 + Ifges(5,6) * t163) * qJD(2);
t170 = mrSges(5,2) * t104 - mrSges(5,3) * t101 + Ifges(5,1) * t139 + Ifges(5,4) * t140 + Ifges(5,5) * qJDD(3) + t117 * t182;
t118 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t161 + Ifges(4,6) * t163) * qJD(2);
t80 = Ifges(4,4) * t139 + Ifges(4,2) * t140 + Ifges(4,6) * qJDD(3) + qJD(3) * t122 - mrSges(4,1) * t108 + mrSges(4,3) * t106 - pkin(3) * (t139 * mrSges(5,2) - t174) + qJ(4) * (-qJDD(3) * mrSges(5,2) - qJD(3) * t146 + t178) + (-pkin(3) * t146 - t117 - t118) * t183 + t172;
t82 = t118 * t182 + mrSges(4,2) * t108 - mrSges(4,3) * t105 + Ifges(4,1) * t139 + Ifges(4,4) * t140 + Ifges(4,5) * qJDD(3) - qJ(4) * t96 + (-t119 - t120) * qJD(3) + t170;
t169 = mrSges(3,1) * t111 - mrSges(3,2) * t112 + Ifges(3,3) * qJDD(2) + pkin(2) * t167 + pkin(5) * t176 + t161 * t82 + t163 * t80;
t168 = mrSges(2,1) * t143 - mrSges(2,2) * t144 + pkin(1) * t78 + t169;
t76 = m(2) * t144 + t175;
t75 = m(2) * t143 + t78;
t74 = -mrSges(3,1) * t158 + mrSges(3,3) * t112 + t165 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t87 - t190;
t73 = mrSges(3,2) * t158 - mrSges(3,3) * t111 + Ifges(3,5) * qJDD(2) - t165 * Ifges(3,6) - pkin(5) * t87 - t161 * t80 + t163 * t82;
t72 = mrSges(2,2) * t158 - mrSges(2,3) * t143 - pkin(4) * t78 - t162 * t74 + t164 * t73;
t71 = -mrSges(2,1) * t158 + mrSges(2,3) * t144 + t162 * t73 + t164 * t74 - pkin(1) * (m(3) * t158 + t87) + pkin(4) * t175;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t160 * t72 - t159 * t71 - qJ(1) * (t159 * t76 + t160 * t75), t72, t73, t82, -qJD(3) * t119 + t170; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t159 * t72 + t160 * t71 + qJ(1) * (-t159 * t75 + t160 * t76), t71, t74, t80, -t117 * t183 + t172; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t168, t168, t169, t190, -t121 * t182 - t171;];
m_new = t1;
