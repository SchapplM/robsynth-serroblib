% Calculate vector of cutting torques with Newton-Euler for
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
% m [3x6]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPPRR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR5_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR5_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR5_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR5_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR5_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:31
% EndTime: 2019-12-31 17:56:33
% DurationCPUTime: 1.39s
% Computational Cost: add. (20239->184), mult. (28059->222), div. (0->0), fcn. (11856->8), ass. (0->80)
t162 = qJD(1) ^ 2;
t158 = sin(qJ(1));
t161 = cos(qJ(1));
t136 = t158 * g(1) - g(2) * t161;
t131 = qJDD(1) * pkin(1) + t136;
t137 = -g(1) * t161 - g(2) * t158;
t132 = -pkin(1) * t162 + t137;
t154 = sin(pkin(8));
t155 = cos(pkin(8));
t119 = t154 * t131 + t155 * t132;
t174 = qJDD(1) * qJ(3) + (2 * qJD(3) * qJD(1)) + t119;
t180 = -pkin(2) - pkin(3);
t110 = t162 * t180 + t174;
t118 = t155 * t131 - t154 * t132;
t171 = -t162 * qJ(3) + qJDD(3) - t118;
t113 = qJDD(1) * t180 + t171;
t157 = sin(qJ(4));
t160 = cos(qJ(4));
t107 = t160 * t110 + t157 * t113;
t143 = -qJD(1) + qJD(4);
t141 = t143 ^ 2;
t142 = -qJDD(1) + qJDD(4);
t104 = -(pkin(4) * t141) + pkin(7) * t142 + t107;
t151 = g(3) - qJDD(2);
t156 = sin(qJ(5));
t159 = cos(qJ(5));
t101 = -t104 * t156 + t151 * t159;
t102 = t104 * t159 + t151 * t156;
t121 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t156 + Ifges(6,2) * t159) * t143;
t122 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t156 + Ifges(6,4) * t159) * t143;
t176 = qJD(5) * t143;
t126 = t142 * t156 + t159 * t176;
t127 = t142 * t159 - t156 * t176;
t181 = mrSges(6,1) * t101 - mrSges(6,2) * t102 + Ifges(6,5) * t126 + Ifges(6,6) * t127 + Ifges(6,3) * qJDD(5) + (t121 * t156 - t122 * t159) * t143;
t179 = mrSges(3,1) + mrSges(4,1);
t114 = -pkin(2) * t162 + t174;
t125 = (-mrSges(6,1) * t159 + mrSges(6,2) * t156) * t143;
t178 = t143 * t156;
t133 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t178;
t177 = t143 * t159;
t100 = m(6) * t102 - qJDD(5) * mrSges(6,2) + t127 * mrSges(6,3) - qJD(5) * t133 + t125 * t177;
t134 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t177;
t99 = m(6) * t101 + qJDD(5) * mrSges(6,1) - t126 * mrSges(6,3) + qJD(5) * t134 - t125 * t178;
t93 = t159 * t100 - t156 * t99;
t89 = m(5) * t107 - (mrSges(5,1) * t141) - mrSges(5,2) * t142 + t93;
t106 = -t110 * t157 + t113 * t160;
t103 = -pkin(4) * t142 - pkin(7) * t141 - t106;
t97 = -m(6) * t103 + t127 * mrSges(6,1) - t126 * mrSges(6,2) - t133 * t178 + t134 * t177;
t96 = m(5) * t106 + mrSges(5,1) * t142 - mrSges(5,2) * t141 + t97;
t87 = -t157 * t96 + t160 * t89;
t172 = m(4) * t114 + qJDD(1) * mrSges(4,3) + t87;
t81 = m(3) * t119 - qJDD(1) * mrSges(3,2) - t162 * t179 + t172;
t116 = -qJDD(1) * pkin(2) + t171;
t86 = t157 * t89 + t160 * t96;
t169 = -m(4) * t116 + qJDD(1) * mrSges(4,1) + t162 * mrSges(4,3) - t86;
t82 = m(3) * t118 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t162 + t169;
t75 = t154 * t81 + t155 * t82;
t175 = -t154 * t82 + t155 * t81;
t92 = t156 * t100 + t159 * t99;
t91 = -m(5) * t151 - t92;
t120 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t156 + Ifges(6,6) * t159) * t143;
t94 = -mrSges(6,1) * t103 + mrSges(6,3) * t102 + Ifges(6,4) * t126 + Ifges(6,2) * t127 + Ifges(6,6) * qJDD(5) + qJD(5) * t122 - t120 * t178;
t95 = mrSges(6,2) * t103 - mrSges(6,3) * t101 + Ifges(6,1) * t126 + Ifges(6,4) * t127 + Ifges(6,5) * qJDD(5) - qJD(5) * t121 + t120 * t177;
t77 = mrSges(5,2) * t151 - mrSges(5,3) * t106 + Ifges(5,5) * t142 - (Ifges(5,6) * t141) - pkin(7) * t92 - t156 * t94 + t159 * t95;
t85 = -mrSges(5,1) * t151 + mrSges(5,3) * t107 + t141 * Ifges(5,5) + Ifges(5,6) * t142 - pkin(4) * t92 - t181;
t170 = mrSges(4,2) * t116 + Ifges(4,4) * qJDD(1) + t162 * Ifges(4,6) - pkin(6) * t86 - t157 * t85 + t160 * t77;
t168 = mrSges(4,2) * t114 - pkin(3) * t91 - pkin(6) * t87 - t157 * t77 - t160 * t85;
t166 = mrSges(5,1) * t106 - mrSges(5,2) * t107 + Ifges(5,3) * t142 + pkin(4) * t97 + pkin(7) * t93 + t156 * t95 + t159 * t94;
t165 = -mrSges(4,1) * t116 + mrSges(4,3) * t114 + Ifges(4,2) * qJDD(1) - pkin(3) * t86 - t166;
t164 = -mrSges(3,2) * t119 + mrSges(3,1) * t118 + Ifges(3,3) * qJDD(1) + t165 + qJ(3) * (-mrSges(4,1) * t162 + t172) + pkin(2) * t169;
t163 = mrSges(2,1) * t136 - mrSges(2,2) * t137 + Ifges(2,3) * qJDD(1) + pkin(1) * t75 + t164;
t138 = m(4) * t151;
t90 = -t138 + t91;
t73 = m(2) * t137 - mrSges(2,1) * t162 - qJDD(1) * mrSges(2,2) + t175;
t72 = m(2) * t136 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t162 + t75;
t71 = -mrSges(3,3) * t118 + Ifges(3,5) * qJDD(1) - t162 * Ifges(3,6) - qJ(3) * t90 + (-mrSges(3,2) + mrSges(4,3)) * t151 + t170;
t70 = mrSges(3,3) * t119 - pkin(2) * t90 + (Ifges(4,4) + Ifges(3,5)) * t162 + t179 * t151 + (Ifges(3,6) - Ifges(4,6)) * qJDD(1) + t168;
t69 = -mrSges(2,2) * g(3) - mrSges(2,3) * t136 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t162 - qJ(2) * t75 - t154 * t70 + t155 * t71;
t68 = Ifges(2,6) * qJDD(1) + t162 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t137 + t154 * t71 + t155 * t70 - pkin(1) * (-t138 + (-m(3) - m(5)) * t151 - t92) + qJ(2) * t175;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t161 * t69 - t158 * t68 - pkin(5) * (t158 * t73 + t161 * t72), t69, t71, mrSges(4,3) * t151 + t170, t77, t95; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t158 * t69 + t161 * t68 + pkin(5) * (-t158 * t72 + t161 * t73), t68, t70, t165, t85, t94; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t163, t163, t164, -mrSges(4,1) * t151 - t162 * Ifges(4,4) + Ifges(4,6) * qJDD(1) - t168, t166, t181;];
m_new = t1;
