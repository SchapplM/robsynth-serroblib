% Calculate vector of cutting torques with Newton-Euler for
% S4PRPR6
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
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PRPR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR6_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR6_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR6_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR6_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:22
% EndTime: 2019-12-31 16:24:23
% DurationCPUTime: 1.19s
% Computational Cost: add. (11565->162), mult. (24722->211), div. (0->0), fcn. (15575->8), ass. (0->78)
t156 = qJD(2) ^ 2;
t150 = sin(pkin(6));
t178 = cos(pkin(6));
t137 = -t178 * g(1) - t150 * g(2);
t148 = -g(3) + qJDD(1);
t153 = sin(qJ(2));
t155 = cos(qJ(2));
t124 = t155 * t137 + t153 * t148;
t120 = -t156 * pkin(2) + qJDD(2) * qJ(3) + t124;
t149 = sin(pkin(7));
t136 = t150 * g(1) - t178 * g(2);
t151 = cos(pkin(7));
t174 = qJD(2) * qJD(3);
t176 = -t151 * t136 - 0.2e1 * t149 * t174;
t105 = -t149 * t120 + t176;
t106 = -t149 * t136 + (t120 + 0.2e1 * t174) * t151;
t180 = pkin(3) * t151;
t102 = (-pkin(5) * qJDD(2) + t156 * t180 - t120) * t149 + t176;
t173 = qJDD(2) * t151;
t146 = t151 ^ 2;
t177 = t146 * t156;
t103 = -pkin(3) * t177 + pkin(5) * t173 + t106;
t152 = sin(qJ(4));
t154 = cos(qJ(4));
t100 = t154 * t102 - t152 * t103;
t101 = t152 * t102 + t154 * t103;
t164 = -t149 * t152 + t151 * t154;
t125 = t164 * qJD(2);
t165 = t149 * t154 + t151 * t152;
t126 = t165 * qJD(2);
t109 = Ifges(5,4) * t126 + Ifges(5,2) * t125 + Ifges(5,6) * qJD(4);
t110 = Ifges(5,1) * t126 + Ifges(5,4) * t125 + Ifges(5,5) * qJD(4);
t116 = -t126 * qJD(4) + t164 * qJDD(2);
t117 = t125 * qJD(4) + t165 * qJDD(2);
t160 = -mrSges(5,1) * t100 + mrSges(5,2) * t101 - Ifges(5,5) * t117 - Ifges(5,6) * t116 - Ifges(5,3) * qJDD(4) - t126 * t109 + t125 * t110;
t169 = Ifges(4,4) * t149 + Ifges(4,2) * t151;
t170 = Ifges(4,1) * t149 + Ifges(4,4) * t151;
t112 = -t125 * mrSges(5,1) + t126 * mrSges(5,2);
t121 = -qJD(4) * mrSges(5,2) + t125 * mrSges(5,3);
t96 = m(5) * t100 + qJDD(4) * mrSges(5,1) - t117 * mrSges(5,3) + qJD(4) * t121 - t126 * t112;
t122 = qJD(4) * mrSges(5,1) - t126 * mrSges(5,3);
t97 = m(5) * t101 - qJDD(4) * mrSges(5,2) + t116 * mrSges(5,3) - qJD(4) * t122 + t125 * t112;
t89 = t152 * t97 + t154 * t96;
t181 = -mrSges(4,1) * t105 + mrSges(4,2) * t106 - pkin(3) * t89 - (t149 * t169 - t151 * t170) * t156 + t160;
t179 = mrSges(4,2) * t149;
t168 = Ifges(4,5) * t149 + Ifges(4,6) * t151;
t175 = t156 * t168;
t163 = mrSges(4,3) * qJDD(2) + t156 * (-mrSges(4,1) * t151 + t179);
t87 = m(4) * t105 - t163 * t149 + t89;
t172 = -t152 * t96 + t154 * t97;
t88 = m(4) * t106 + t163 * t151 + t172;
t85 = -t149 * t87 + t151 * t88;
t81 = m(3) * t124 - t156 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t85;
t123 = -t153 * t137 + t155 * t148;
t167 = qJDD(3) - t123;
t119 = -qJDD(2) * pkin(2) - t156 * qJ(3) + t167;
t145 = t149 ^ 2;
t107 = (-pkin(2) - t180) * qJDD(2) + (-qJ(3) + (-t145 - t146) * pkin(5)) * t156 + t167;
t161 = m(5) * t107 - t116 * mrSges(5,1) + t117 * mrSges(5,2) - t125 * t121 + t126 * t122;
t159 = -m(4) * t119 + mrSges(4,1) * t173 - t161 + (t145 * t156 + t177) * mrSges(4,3);
t92 = m(3) * t123 - t156 * mrSges(3,2) + (mrSges(3,1) - t179) * qJDD(2) + t159;
t171 = -t153 * t92 + t155 * t81;
t84 = t149 * t88 + t151 * t87;
t108 = Ifges(5,5) * t126 + Ifges(5,6) * t125 + Ifges(5,3) * qJD(4);
t90 = -mrSges(5,1) * t107 + mrSges(5,3) * t101 + Ifges(5,4) * t117 + Ifges(5,2) * t116 + Ifges(5,6) * qJDD(4) + qJD(4) * t110 - t126 * t108;
t91 = mrSges(5,2) * t107 - mrSges(5,3) * t100 + Ifges(5,1) * t117 + Ifges(5,4) * t116 + Ifges(5,5) * qJDD(4) - qJD(4) * t109 + t125 * t108;
t75 = -mrSges(4,1) * t119 + mrSges(4,3) * t106 - pkin(3) * t161 + pkin(5) * t172 + t169 * qJDD(2) - t149 * t175 + t152 * t91 + t154 * t90;
t79 = mrSges(4,2) * t119 - mrSges(4,3) * t105 - pkin(5) * t89 + t170 * qJDD(2) + t151 * t175 - t152 * t90 + t154 * t91;
t72 = -mrSges(3,2) * t136 - mrSges(3,3) * t123 + Ifges(3,5) * qJDD(2) - t156 * Ifges(3,6) - qJ(3) * t84 - t149 * t75 + t151 * t79;
t74 = (Ifges(3,6) - t168) * qJDD(2) + t156 * Ifges(3,5) + mrSges(3,1) * t136 + mrSges(3,3) * t124 - pkin(2) * t84 + t181;
t162 = -mrSges(2,2) * t137 + pkin(4) * t171 + t153 * t72 + t155 * t74 + mrSges(2,1) * t136 + pkin(1) * (m(3) * t136 - t84);
t157 = mrSges(3,1) * t123 - mrSges(3,2) * t124 + Ifges(3,3) * qJDD(2) + pkin(2) * (-qJDD(2) * t179 + t159) + qJ(3) * t85 + t149 * t79 + t151 * t75;
t82 = (m(2) + m(3)) * t136 - t84;
t78 = t153 * t81 + t155 * t92;
t76 = m(2) * t137 + t171;
t70 = -mrSges(2,1) * t148 + mrSges(2,3) * t137 - pkin(1) * t78 - t157;
t69 = mrSges(2,2) * t148 - mrSges(2,3) * t136 - pkin(4) * t78 - t153 * t74 + t155 * t72;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t178 * t69 - t150 * t70 - qJ(1) * (t150 * t76 + t178 * t82), t69, t72, t79, t91; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t150 * t69 + t178 * t70 + qJ(1) * (-t150 * t82 + t178 * t76), t70, t74, t75, t90; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t162, t162, t157, t168 * qJDD(2) - t181, -t160;];
m_new = t1;
