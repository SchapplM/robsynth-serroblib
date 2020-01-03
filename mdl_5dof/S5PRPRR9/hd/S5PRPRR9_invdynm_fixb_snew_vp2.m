% Calculate vector of cutting torques with Newton-Euler for
% S5PRPRR9
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% m [3x6]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRPRR9_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR9_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR9_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR9_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR9_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR9_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR9_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:39
% EndTime: 2019-12-31 17:39:41
% DurationCPUTime: 1.29s
% Computational Cost: add. (17807->173), mult. (24403->211), div. (0->0), fcn. (11856->8), ass. (0->78)
t157 = qJD(2) ^ 2;
t149 = sin(pkin(8));
t150 = cos(pkin(8));
t133 = t149 * g(1) - t150 * g(2);
t134 = -t150 * g(1) - t149 * g(2);
t153 = sin(qJ(2));
t156 = cos(qJ(2));
t118 = t153 * t133 + t156 * t134;
t169 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t118;
t175 = -pkin(2) - pkin(3);
t109 = t175 * t157 + t169;
t117 = t156 * t133 - t153 * t134;
t166 = -t157 * qJ(3) + qJDD(3) - t117;
t112 = t175 * qJDD(2) + t166;
t152 = sin(qJ(4));
t155 = cos(qJ(4));
t106 = t155 * t109 + t152 * t112;
t139 = -qJD(2) + qJD(4);
t137 = t139 ^ 2;
t138 = -qJDD(2) + qJDD(4);
t103 = -(t137 * pkin(4)) + t138 * pkin(7) + t106;
t146 = g(3) - qJDD(1);
t151 = sin(qJ(5));
t154 = cos(qJ(5));
t100 = -t151 * t103 + t154 * t146;
t101 = t154 * t103 + t151 * t146;
t120 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t151 + Ifges(6,2) * t154) * t139;
t121 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t151 + Ifges(6,4) * t154) * t139;
t171 = qJD(5) * t139;
t125 = t151 * t138 + t154 * t171;
t126 = t154 * t138 - t151 * t171;
t176 = mrSges(6,1) * t100 - mrSges(6,2) * t101 + Ifges(6,5) * t125 + Ifges(6,6) * t126 + Ifges(6,3) * qJDD(5) + (t120 * t151 - t121 * t154) * t139;
t174 = mrSges(3,1) + mrSges(4,1);
t113 = -t157 * pkin(2) + t169;
t124 = (-mrSges(6,1) * t154 + mrSges(6,2) * t151) * t139;
t172 = t139 * t154;
t131 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t172;
t173 = t139 * t151;
t98 = m(6) * t100 + qJDD(5) * mrSges(6,1) - t125 * mrSges(6,3) + qJD(5) * t131 - t124 * t173;
t130 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t173;
t99 = m(6) * t101 - qJDD(5) * mrSges(6,2) + t126 * mrSges(6,3) - qJD(5) * t130 + t124 * t172;
t92 = -t151 * t98 + t154 * t99;
t88 = m(5) * t106 - (t137 * mrSges(5,1)) - t138 * mrSges(5,2) + t92;
t105 = -t152 * t109 + t155 * t112;
t102 = -t138 * pkin(4) - t137 * pkin(7) - t105;
t96 = -m(6) * t102 + t126 * mrSges(6,1) - t125 * mrSges(6,2) - t130 * t173 + t131 * t172;
t95 = m(5) * t105 + t138 * mrSges(5,1) - t137 * mrSges(5,2) + t96;
t86 = -t152 * t95 + t155 * t88;
t167 = m(4) * t113 + qJDD(2) * mrSges(4,3) + t86;
t80 = m(3) * t118 - qJDD(2) * mrSges(3,2) - t174 * t157 + t167;
t115 = -qJDD(2) * pkin(2) + t166;
t85 = t152 * t88 + t155 * t95;
t164 = -m(4) * t115 + qJDD(2) * mrSges(4,1) + t157 * mrSges(4,3) - t85;
t81 = m(3) * t117 + qJDD(2) * mrSges(3,1) - t157 * mrSges(3,2) + t164;
t74 = t153 * t80 + t156 * t81;
t170 = -t153 * t81 + t156 * t80;
t91 = t151 * t99 + t154 * t98;
t90 = -m(5) * t146 - t91;
t119 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t151 + Ifges(6,6) * t154) * t139;
t93 = -mrSges(6,1) * t102 + mrSges(6,3) * t101 + Ifges(6,4) * t125 + Ifges(6,2) * t126 + Ifges(6,6) * qJDD(5) + qJD(5) * t121 - t119 * t173;
t94 = mrSges(6,2) * t102 - mrSges(6,3) * t100 + Ifges(6,1) * t125 + Ifges(6,4) * t126 + Ifges(6,5) * qJDD(5) - qJD(5) * t120 + t119 * t172;
t76 = mrSges(5,2) * t146 - mrSges(5,3) * t105 + Ifges(5,5) * t138 - (t137 * Ifges(5,6)) - pkin(7) * t91 - t151 * t93 + t154 * t94;
t84 = -mrSges(5,1) * t146 + mrSges(5,3) * t106 + t137 * Ifges(5,5) + Ifges(5,6) * t138 - pkin(4) * t91 - t176;
t165 = mrSges(4,2) * t115 + Ifges(4,4) * qJDD(2) + t157 * Ifges(4,6) - pkin(6) * t85 - t152 * t84 + t155 * t76;
t163 = mrSges(4,2) * t113 - pkin(3) * t90 - pkin(6) * t86 - t152 * t76 - t155 * t84;
t161 = mrSges(5,1) * t105 - mrSges(5,2) * t106 + Ifges(5,3) * t138 + pkin(4) * t96 + pkin(7) * t92 + t151 * t94 + t154 * t93;
t160 = -mrSges(4,1) * t115 + mrSges(4,3) * t113 + Ifges(4,2) * qJDD(2) - pkin(3) * t85 - t161;
t159 = -mrSges(3,2) * t118 + mrSges(3,1) * t117 + Ifges(3,3) * qJDD(2) + t160 + qJ(3) * (-t157 * mrSges(4,1) + t167) + pkin(2) * t164;
t158 = mrSges(2,1) * t133 - mrSges(2,2) * t134 + pkin(1) * t74 + t159;
t135 = m(4) * t146;
t89 = -t135 + t90;
t72 = m(2) * t134 + t170;
t71 = m(2) * t133 + t74;
t70 = -mrSges(3,3) * t117 + Ifges(3,5) * qJDD(2) - t157 * Ifges(3,6) - qJ(3) * t89 + (-mrSges(3,2) + mrSges(4,3)) * t146 + t165;
t69 = mrSges(3,3) * t118 - pkin(2) * t89 + (Ifges(4,4) + Ifges(3,5)) * t157 + t174 * t146 + (Ifges(3,6) - Ifges(4,6)) * qJDD(2) + t163;
t68 = -mrSges(2,2) * t146 - mrSges(2,3) * t133 - pkin(5) * t74 - t153 * t69 + t156 * t70;
t67 = mrSges(2,3) * t134 + t153 * t70 + t156 * t69 - pkin(1) * (-t135 - t91) + pkin(5) * t170 + (mrSges(2,1) - pkin(1) * (-m(3) - m(5))) * t146;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t150 * t68 - t149 * t67 - qJ(1) * (t149 * t72 + t150 * t71), t68, t70, mrSges(4,3) * t146 + t165, t76, t94; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t149 * t68 + t150 * t67 + qJ(1) * (-t149 * t71 + t150 * t72), t67, t69, t160, t84, t93; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t158, t158, t159, -mrSges(4,1) * t146 - t157 * Ifges(4,4) + Ifges(4,6) * qJDD(2) - t163, t161, t176;];
m_new = t1;
