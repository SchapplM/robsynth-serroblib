% Calculate vector of cutting torques with Newton-Euler for
% S4PRRP5
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
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PRRP5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP5_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP5_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP5_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP5_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:28:52
% EndTime: 2019-12-31 16:28:54
% DurationCPUTime: 0.78s
% Computational Cost: add. (5695->179), mult. (10838->224), div. (0->0), fcn. (5206->6), ass. (0->72)
t155 = sin(pkin(6));
t180 = cos(pkin(6));
t141 = -t180 * g(1) - t155 * g(2);
t154 = -g(3) + qJDD(1);
t157 = sin(qJ(2));
t159 = cos(qJ(2));
t108 = t159 * t141 + t157 * t154;
t160 = qJD(2) ^ 2;
t106 = -t160 * pkin(2) + qJDD(2) * pkin(5) + t108;
t140 = t155 * g(1) - t180 * g(2);
t158 = cos(qJ(3));
t132 = t158 * t140;
t156 = sin(qJ(3));
t102 = -t156 * t106 - t132;
t103 = t158 * t106 - t156 * t140;
t116 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t156 + Ifges(4,2) * t158) * qJD(2);
t117 = Ifges(5,5) * qJD(3) + (Ifges(5,1) * t156 + Ifges(5,4) * t158) * qJD(2);
t118 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t156 + Ifges(4,4) * t158) * qJD(2);
t174 = qJD(2) * qJD(3);
t170 = t158 * t174;
t136 = t156 * qJDD(2) + t170;
t137 = t158 * qJDD(2) - t156 * t174;
t115 = Ifges(5,6) * qJD(3) + (Ifges(5,4) * t156 + Ifges(5,2) * t158) * qJD(2);
t176 = qJD(2) * t156;
t173 = qJD(2) * qJD(4);
t182 = pkin(3) * t160;
t98 = qJDD(3) * pkin(3) - t132 + (-t136 + t170) * qJ(4) + (t158 * t182 - t106 - 0.2e1 * t173) * t156;
t142 = qJD(3) * pkin(3) - qJ(4) * t176;
t153 = t158 ^ 2;
t99 = t137 * qJ(4) - qJD(3) * t142 - t153 * t182 + 0.2e1 * t158 * t173 + t103;
t165 = -mrSges(5,1) * t98 + mrSges(5,2) * t99 - Ifges(5,5) * t136 - Ifges(5,6) * t137 - Ifges(5,3) * qJDD(3) - t115 * t176;
t134 = (-mrSges(5,1) * t158 + mrSges(5,2) * t156) * qJD(2);
t175 = qJD(2) * t158;
t145 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t175;
t171 = m(5) * t98 + qJDD(3) * mrSges(5,1) + qJD(3) * t145;
t93 = -t136 * mrSges(5,3) - t134 * t176 + t171;
t184 = mrSges(4,1) * t102 - mrSges(4,2) * t103 + Ifges(4,5) * t136 + Ifges(4,6) * t137 + Ifges(4,3) * qJDD(3) + pkin(3) * t93 - (-t156 * t116 + (t117 + t118) * t158) * qJD(2) - t165;
t181 = -mrSges(4,2) - mrSges(5,2);
t143 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t176;
t177 = -qJD(3) * mrSges(4,1) + mrSges(4,3) * t176 - t143;
t172 = m(5) * t99 + t137 * mrSges(5,3) + t134 * t175;
t135 = (-mrSges(4,1) * t158 + mrSges(4,2) * t156) * qJD(2);
t146 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t175;
t90 = m(4) * t102 + qJDD(3) * mrSges(4,1) + qJD(3) * t146 + (-mrSges(4,3) - mrSges(5,3)) * t136 + (-t134 - t135) * t176 + t171;
t91 = m(4) * t103 + t137 * mrSges(4,3) + t177 * qJD(3) + t181 * qJDD(3) + t135 * t175 + t172;
t87 = -t156 * t90 + t158 * t91;
t83 = m(3) * t108 - t160 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t87;
t107 = -t157 * t141 + t159 * t154;
t167 = -qJDD(2) * pkin(2) - t107;
t105 = -t160 * pkin(5) + t167;
t101 = t142 * t176 - t137 * pkin(3) + qJDD(4) + (-qJ(4) * t153 - pkin(5)) * t160 + t167;
t168 = -m(5) * t101 + t137 * mrSges(5,1) + t145 * t175;
t92 = -m(4) * t105 + t137 * mrSges(4,1) + t181 * t136 + t146 * t175 + t177 * t176 + t168;
t88 = m(3) * t107 + qJDD(2) * mrSges(3,1) - t160 * mrSges(3,2) + t92;
t169 = -t157 * t88 + t159 * t83;
t86 = t156 * t91 + t158 * t90;
t113 = Ifges(5,3) * qJD(3) + (Ifges(5,5) * t156 + Ifges(5,6) * t158) * qJD(2);
t114 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t156 + Ifges(4,6) * t158) * qJD(2);
t164 = -mrSges(5,1) * t101 + mrSges(5,3) * t99 + Ifges(5,4) * t136 + Ifges(5,2) * t137 + Ifges(5,6) * qJDD(3) + qJD(3) * t117;
t80 = Ifges(4,4) * t136 + Ifges(4,2) * t137 + Ifges(4,6) * qJDD(3) + qJD(3) * t118 - mrSges(4,1) * t105 + mrSges(4,3) * t103 - pkin(3) * (t136 * mrSges(5,2) - t168) + qJ(4) * (-qJDD(3) * mrSges(5,2) - qJD(3) * t143 + t172) + (-pkin(3) * t143 - t113 - t114) * t176 + t164;
t163 = mrSges(5,2) * t101 - mrSges(5,3) * t98 + Ifges(5,1) * t136 + Ifges(5,4) * t137 + Ifges(5,5) * qJDD(3) + t113 * t175;
t81 = t114 * t175 + mrSges(4,2) * t105 - mrSges(4,3) * t102 + Ifges(4,1) * t136 + Ifges(4,4) * t137 + Ifges(4,5) * qJDD(3) - qJ(4) * t93 + (-t115 - t116) * qJD(3) + t163;
t74 = -mrSges(3,2) * t140 - mrSges(3,3) * t107 + Ifges(3,5) * qJDD(2) - t160 * Ifges(3,6) - pkin(5) * t86 - t156 * t80 + t158 * t81;
t76 = mrSges(3,1) * t140 + mrSges(3,3) * t108 + t160 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t86 - t184;
t166 = -mrSges(2,2) * t141 + pkin(4) * t169 + t157 * t74 + t159 * t76 + mrSges(2,1) * t140 + pkin(1) * (m(3) * t140 - t86);
t162 = mrSges(3,1) * t107 - mrSges(3,2) * t108 + Ifges(3,3) * qJDD(2) + pkin(2) * t92 + pkin(5) * t87 + t156 * t81 + t158 * t80;
t84 = (m(2) + m(3)) * t140 - t86;
t79 = t157 * t83 + t159 * t88;
t77 = m(2) * t141 + t169;
t72 = -mrSges(2,1) * t154 + mrSges(2,3) * t141 - pkin(1) * t79 - t162;
t71 = mrSges(2,2) * t154 - mrSges(2,3) * t140 - pkin(4) * t79 - t157 * t76 + t159 * t74;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t180 * t71 - t155 * t72 - qJ(1) * (t155 * t77 + t180 * t84), t71, t74, t81, -qJD(3) * t115 + t163; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t155 * t71 + t180 * t72 + qJ(1) * (-t155 * t84 + t180 * t77), t72, t76, t80, -t113 * t176 + t164; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t166, t166, t162, t184, -t117 * t175 - t165;];
m_new = t1;
