% Calculate vector of cutting torques with Newton-Euler for
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
% m [3x5]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PRRP6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP6_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP6_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP6_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP6_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:14
% EndTime: 2019-12-31 16:30:16
% DurationCPUTime: 0.80s
% Computational Cost: add. (5546->179), mult. (10303->227), div. (0->0), fcn. (4932->6), ass. (0->67)
t149 = sin(qJ(3));
t151 = cos(qJ(3));
t108 = Ifges(5,6) * qJD(3) + (Ifges(5,5) * t149 - Ifges(5,3) * t151) * qJD(2);
t111 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t149 + Ifges(4,2) * t151) * qJD(2);
t127 = (-mrSges(5,1) * t151 - mrSges(5,3) * t149) * qJD(2);
t164 = qJD(2) * qJD(3);
t129 = t149 * qJDD(2) + t151 * t164;
t130 = t151 * qJDD(2) - t149 * t164;
t126 = (-pkin(3) * t151 - qJ(4) * t149) * qJD(2);
t153 = qJD(3) ^ 2;
t165 = qJD(2) * t151;
t148 = sin(pkin(6));
t170 = cos(pkin(6));
t135 = -t170 * g(1) - t148 * g(2);
t147 = -g(3) + qJDD(1);
t150 = sin(qJ(2));
t152 = cos(qJ(2));
t104 = t152 * t135 + t150 * t147;
t154 = qJD(2) ^ 2;
t102 = -t154 * pkin(2) + qJDD(2) * pkin(5) + t104;
t134 = t148 * g(1) - t170 * g(2);
t99 = t151 * t102 - t149 * t134;
t95 = -t153 * pkin(3) + qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) + t126 * t165 + t99;
t169 = t151 * t134;
t97 = -qJDD(3) * pkin(3) - t153 * qJ(4) + t169 + qJDD(4) + (qJD(2) * t126 + t102) * t149;
t159 = -mrSges(5,1) * t97 + mrSges(5,3) * t95 + Ifges(5,4) * t129 + Ifges(5,2) * qJDD(3) - Ifges(5,6) * t130;
t139 = mrSges(5,2) * t165 + qJD(3) * mrSges(5,3);
t161 = -m(5) * t97 + qJDD(3) * mrSges(5,1) + qJD(3) * t139;
t166 = qJD(2) * t149;
t137 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t166;
t162 = m(5) * t95 + qJDD(3) * mrSges(5,3) + qJD(3) * t137 + t127 * t165;
t112 = Ifges(5,4) * qJD(3) + (Ifges(5,1) * t149 - Ifges(5,5) * t151) * qJD(2);
t167 = t112 + Ifges(4,5) * qJD(3) + (Ifges(4,1) * t149 + Ifges(4,4) * t151) * qJD(2);
t98 = -t149 * t102 - t169;
t173 = -(t167 * t151 + (t108 - t111) * t149) * qJD(2) + mrSges(4,1) * t98 - mrSges(4,2) * t99 + Ifges(4,5) * t129 + Ifges(4,6) * t130 + Ifges(4,3) * qJDD(3) + pkin(3) * (-t129 * mrSges(5,2) - t127 * t166 + t161) + qJ(4) * (t130 * mrSges(5,2) + t162) + t159;
t171 = mrSges(4,3) + mrSges(5,2);
t128 = (-mrSges(4,1) * t151 + mrSges(4,2) * t149) * qJD(2);
t136 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t166;
t87 = m(4) * t99 - qJDD(3) * mrSges(4,2) - qJD(3) * t136 + t128 * t165 + t171 * t130 + t162;
t138 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t165;
t88 = m(4) * t98 + qJDD(3) * mrSges(4,1) + qJD(3) * t138 - t171 * t129 + (-t127 - t128) * t166 + t161;
t83 = -t149 * t88 + t151 * t87;
t79 = m(3) * t104 - t154 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t83;
t103 = -t150 * t135 + t152 * t147;
t101 = -qJDD(2) * pkin(2) - t154 * pkin(5) - t103;
t92 = -t130 * pkin(3) - t129 * qJ(4) + (-0.2e1 * qJD(4) * t149 + (pkin(3) * t149 - qJ(4) * t151) * qJD(3)) * qJD(2) + t101;
t89 = m(5) * t92 - t130 * mrSges(5,1) - t129 * mrSges(5,3) - t137 * t166 - t139 * t165;
t85 = -m(4) * t101 + t130 * mrSges(4,1) - t129 * mrSges(4,2) - t136 * t166 + t138 * t165 - t89;
t84 = m(3) * t103 + qJDD(2) * mrSges(3,1) - t154 * mrSges(3,2) + t85;
t163 = -t150 * t84 + t152 * t79;
t160 = -mrSges(5,1) * t92 + mrSges(5,2) * t95;
t82 = t149 * t87 + t151 * t88;
t109 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t149 + Ifges(4,6) * t151) * qJD(2);
t110 = Ifges(5,2) * qJD(3) + (Ifges(5,4) * t149 - Ifges(5,6) * t151) * qJD(2);
t76 = -mrSges(4,1) * t101 + mrSges(4,3) * t99 - pkin(3) * t89 + (Ifges(4,2) + Ifges(5,3)) * t130 + (Ifges(4,4) - Ifges(5,5)) * t129 + (Ifges(4,6) - Ifges(5,6)) * qJDD(3) + t167 * qJD(3) + (-t109 - t110) * t166 + t160;
t157 = mrSges(5,2) * t97 - mrSges(5,3) * t92 + Ifges(5,1) * t129 + Ifges(5,4) * qJDD(3) - Ifges(5,5) * t130 + qJD(3) * t108 + t110 * t165;
t77 = mrSges(4,2) * t101 - mrSges(4,3) * t98 + Ifges(4,1) * t129 + Ifges(4,4) * t130 + Ifges(4,5) * qJDD(3) - qJ(4) * t89 - qJD(3) * t111 + t109 * t165 + t157;
t70 = -mrSges(3,2) * t134 - mrSges(3,3) * t103 + Ifges(3,5) * qJDD(2) - t154 * Ifges(3,6) - pkin(5) * t82 - t149 * t76 + t151 * t77;
t72 = mrSges(3,1) * t134 + mrSges(3,3) * t104 + t154 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t82 - t173;
t158 = -mrSges(2,2) * t135 + pkin(4) * t163 + t150 * t70 + t152 * t72 + mrSges(2,1) * t134 + pkin(1) * (m(3) * t134 - t82);
t156 = mrSges(3,1) * t103 - mrSges(3,2) * t104 + Ifges(3,3) * qJDD(2) + pkin(2) * t85 + pkin(5) * t83 + t149 * t77 + t151 * t76;
t80 = (m(2) + m(3)) * t134 - t82;
t75 = t150 * t79 + t152 * t84;
t73 = m(2) * t135 + t163;
t68 = -mrSges(2,1) * t147 + mrSges(2,3) * t135 - pkin(1) * t75 - t156;
t67 = mrSges(2,2) * t147 - mrSges(2,3) * t134 - pkin(4) * t75 - t150 * t72 + t152 * t70;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t170 * t67 - t148 * t68 - qJ(1) * (t148 * t73 + t170 * t80), t67, t70, t77, t157; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t148 * t67 + t170 * t68 + qJ(1) * (-t148 * t80 + t170 * t73), t68, t72, t76, (-t149 * t108 - t151 * t112) * qJD(2) + t159; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t158, t158, t156, t173, Ifges(5,5) * t129 + Ifges(5,6) * qJDD(3) - Ifges(5,3) * t130 - qJD(3) * t112 + t110 * t166 - t160;];
m_new = t1;
