% Calculate vector of cutting torques with Newton-Euler for
% S4PRRR6
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PRRR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR6_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR6_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR6_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR6_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR6_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:34:42
% EndTime: 2019-12-31 16:34:44
% DurationCPUTime: 1.33s
% Computational Cost: add. (14413->183), mult. (28248->240), div. (0->0), fcn. (17051->8), ass. (0->79)
t153 = sin(pkin(7));
t174 = cos(pkin(7));
t142 = -t174 * g(1) - t153 * g(2);
t152 = -g(3) + qJDD(1);
t156 = sin(qJ(2));
t159 = cos(qJ(2));
t125 = t159 * t142 + t156 * t152;
t160 = qJD(2) ^ 2;
t121 = -t160 * pkin(2) + qJDD(2) * pkin(5) + t125;
t141 = t153 * g(1) - t174 * g(2);
t155 = sin(qJ(3));
t158 = cos(qJ(3));
t112 = -t155 * t121 - t158 * t141;
t113 = t158 * t121 - t155 * t141;
t127 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t155 + Ifges(4,2) * t158) * qJD(2);
t128 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t155 + Ifges(4,4) * t158) * qJD(2);
t171 = qJD(2) * qJD(3);
t170 = t158 * t171;
t138 = t155 * qJDD(2) + t170;
t139 = t158 * qJDD(2) - t155 * t171;
t103 = (-t138 + t170) * pkin(6) + (t155 * t158 * t160 + qJDD(3)) * pkin(3) + t112;
t173 = qJD(2) * t155;
t145 = qJD(3) * pkin(3) - pkin(6) * t173;
t151 = t158 ^ 2;
t104 = -t151 * t160 * pkin(3) + t139 * pkin(6) - qJD(3) * t145 + t113;
t154 = sin(qJ(4));
t157 = cos(qJ(4));
t101 = t157 * t103 - t154 * t104;
t102 = t154 * t103 + t157 * t104;
t130 = (t154 * t158 + t155 * t157) * qJD(2);
t110 = -t130 * qJD(4) - t154 * t138 + t157 * t139;
t129 = (-t154 * t155 + t157 * t158) * qJD(2);
t111 = t129 * qJD(4) + t157 * t138 + t154 * t139;
t150 = qJD(3) + qJD(4);
t115 = Ifges(5,4) * t130 + Ifges(5,2) * t129 + Ifges(5,6) * t150;
t116 = Ifges(5,1) * t130 + Ifges(5,4) * t129 + Ifges(5,5) * t150;
t149 = qJDD(3) + qJDD(4);
t163 = -mrSges(5,1) * t101 + mrSges(5,2) * t102 - Ifges(5,5) * t111 - Ifges(5,6) * t110 - Ifges(5,3) * t149 - t130 * t115 + t129 * t116;
t118 = -t129 * mrSges(5,1) + t130 * mrSges(5,2);
t122 = -t150 * mrSges(5,2) + t129 * mrSges(5,3);
t98 = m(5) * t101 + t149 * mrSges(5,1) - t111 * mrSges(5,3) - t130 * t118 + t150 * t122;
t123 = t150 * mrSges(5,1) - t130 * mrSges(5,3);
t99 = m(5) * t102 - t149 * mrSges(5,2) + t110 * mrSges(5,3) + t129 * t118 - t150 * t123;
t90 = t154 * t99 + t157 * t98;
t175 = mrSges(4,1) * t112 - mrSges(4,2) * t113 + Ifges(4,5) * t138 + Ifges(4,6) * t139 + Ifges(4,3) * qJDD(3) + pkin(3) * t90 + (t155 * t127 - t158 * t128) * qJD(2) - t163;
t172 = qJD(2) * t158;
t169 = -t154 * t98 + t157 * t99;
t137 = (-mrSges(4,1) * t158 + mrSges(4,2) * t155) * qJD(2);
t144 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t172;
t88 = m(4) * t112 + qJDD(3) * mrSges(4,1) - t138 * mrSges(4,3) + qJD(3) * t144 - t137 * t173 + t90;
t143 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t173;
t89 = m(4) * t113 - qJDD(3) * mrSges(4,2) + t139 * mrSges(4,3) - qJD(3) * t143 + t137 * t172 + t169;
t86 = -t155 * t88 + t158 * t89;
t82 = m(3) * t125 - t160 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t86;
t124 = -t156 * t142 + t159 * t152;
t166 = -qJDD(2) * pkin(2) - t124;
t120 = -t160 * pkin(5) + t166;
t105 = t145 * t173 - t139 * pkin(3) + (-pkin(6) * t151 - pkin(5)) * t160 + t166;
t164 = m(5) * t105 - t110 * mrSges(5,1) + t111 * mrSges(5,2) - t129 * t122 + t130 * t123;
t94 = -m(4) * t120 + t139 * mrSges(4,1) - t138 * mrSges(4,2) - t143 * t173 + t144 * t172 - t164;
t93 = m(3) * t124 + qJDD(2) * mrSges(3,1) - t160 * mrSges(3,2) + t94;
t168 = -t156 * t93 + t159 * t82;
t85 = t155 * t89 + t158 * t88;
t126 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t155 + Ifges(4,6) * t158) * qJD(2);
t114 = Ifges(5,5) * t130 + Ifges(5,6) * t129 + Ifges(5,3) * t150;
t91 = -mrSges(5,1) * t105 + mrSges(5,3) * t102 + Ifges(5,4) * t111 + Ifges(5,2) * t110 + Ifges(5,6) * t149 - t130 * t114 + t150 * t116;
t92 = mrSges(5,2) * t105 - mrSges(5,3) * t101 + Ifges(5,1) * t111 + Ifges(5,4) * t110 + Ifges(5,5) * t149 + t129 * t114 - t150 * t115;
t76 = -mrSges(4,1) * t120 + mrSges(4,3) * t113 + Ifges(4,4) * t138 + Ifges(4,2) * t139 + Ifges(4,6) * qJDD(3) - pkin(3) * t164 + pkin(6) * t169 + qJD(3) * t128 - t126 * t173 + t154 * t92 + t157 * t91;
t80 = mrSges(4,2) * t120 - mrSges(4,3) * t112 + Ifges(4,1) * t138 + Ifges(4,4) * t139 + Ifges(4,5) * qJDD(3) - pkin(6) * t90 - qJD(3) * t127 + t126 * t172 - t154 * t91 + t157 * t92;
t73 = -mrSges(3,2) * t141 - mrSges(3,3) * t124 + Ifges(3,5) * qJDD(2) - t160 * Ifges(3,6) - pkin(5) * t85 - t155 * t76 + t158 * t80;
t75 = mrSges(3,1) * t141 + mrSges(3,3) * t125 + t160 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t85 - t175;
t165 = -mrSges(2,2) * t142 + pkin(4) * t168 + t156 * t73 + t159 * t75 + mrSges(2,1) * t141 + pkin(1) * (m(3) * t141 - t85);
t162 = mrSges(3,1) * t124 - mrSges(3,2) * t125 + Ifges(3,3) * qJDD(2) + pkin(2) * t94 + pkin(5) * t86 + t155 * t80 + t158 * t76;
t83 = (m(2) + m(3)) * t141 - t85;
t79 = t156 * t82 + t159 * t93;
t77 = m(2) * t142 + t168;
t71 = -mrSges(2,1) * t152 + mrSges(2,3) * t142 - pkin(1) * t79 - t162;
t70 = mrSges(2,2) * t152 - mrSges(2,3) * t141 - pkin(4) * t79 - t156 * t75 + t159 * t73;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t174 * t70 - t153 * t71 - qJ(1) * (t153 * t77 + t174 * t83), t70, t73, t80, t92; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t153 * t70 + t174 * t71 + qJ(1) * (-t153 * t83 + t174 * t77), t71, t75, t76, t91; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t165, t165, t162, t175, -t163;];
m_new = t1;
