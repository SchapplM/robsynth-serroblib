% Calculate vector of cutting torques with Newton-Euler for
% S4RPRR4
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RPRR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR4_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR4_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR4_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR4_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR4_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:17
% EndTime: 2019-12-31 16:50:19
% DurationCPUTime: 1.36s
% Computational Cost: add. (14838->192), mult. (27676->247), div. (0->0), fcn. (15216->8), ass. (0->82)
t155 = sin(qJ(1));
t158 = cos(qJ(1));
t143 = t155 * g(1) - t158 * g(2);
t134 = qJDD(1) * pkin(1) + t143;
t144 = -t158 * g(1) - t155 * g(2);
t160 = qJD(1) ^ 2;
t136 = -t160 * pkin(1) + t144;
t151 = sin(pkin(7));
t152 = cos(pkin(7));
t119 = t151 * t134 + t152 * t136;
t109 = -t160 * pkin(2) + qJDD(1) * pkin(5) + t119;
t154 = sin(qJ(3));
t150 = -g(3) + qJDD(2);
t157 = cos(qJ(3));
t174 = t157 * t150;
t105 = -t154 * t109 + t174;
t106 = t157 * t109 + t154 * t150;
t127 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t154 + Ifges(4,2) * t157) * qJD(1);
t128 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t154 + Ifges(4,4) * t157) * qJD(1);
t171 = qJD(1) * qJD(3);
t169 = t157 * t171;
t138 = t154 * qJDD(1) + t169;
t170 = t154 * t171;
t139 = t157 * qJDD(1) - t170;
t153 = sin(qJ(4));
t156 = cos(qJ(4));
t118 = t152 * t134 - t151 * t136;
t108 = -qJDD(1) * pkin(2) - t160 * pkin(5) - t118;
t102 = (-t138 - t169) * pkin(6) + (-t139 + t170) * pkin(3) + t108;
t137 = (-pkin(3) * t157 - pkin(6) * t154) * qJD(1);
t159 = qJD(3) ^ 2;
t172 = t157 * qJD(1);
t104 = -t159 * pkin(3) + qJDD(3) * pkin(6) + t137 * t172 + t106;
t100 = t156 * t102 - t153 * t104;
t173 = qJD(1) * t154;
t132 = t156 * qJD(3) - t153 * t173;
t116 = t132 * qJD(4) + t153 * qJDD(3) + t156 * t138;
t133 = t153 * qJD(3) + t156 * t173;
t120 = -t132 * mrSges(5,1) + t133 * mrSges(5,2);
t145 = qJD(4) - t172;
t121 = -t145 * mrSges(5,2) + t132 * mrSges(5,3);
t131 = qJDD(4) - t139;
t97 = m(5) * t100 + t131 * mrSges(5,1) - t116 * mrSges(5,3) - t133 * t120 + t145 * t121;
t101 = t153 * t102 + t156 * t104;
t115 = -t133 * qJD(4) + t156 * qJDD(3) - t153 * t138;
t122 = t145 * mrSges(5,1) - t133 * mrSges(5,3);
t98 = m(5) * t101 - t131 * mrSges(5,2) + t115 * mrSges(5,3) + t132 * t120 - t145 * t122;
t91 = -t153 * t97 + t156 * t98;
t103 = -qJDD(3) * pkin(3) - t159 * pkin(6) - t174 + (qJD(1) * t137 + t109) * t154;
t110 = Ifges(5,5) * t133 + Ifges(5,6) * t132 + Ifges(5,3) * t145;
t112 = Ifges(5,1) * t133 + Ifges(5,4) * t132 + Ifges(5,5) * t145;
t92 = -mrSges(5,1) * t103 + mrSges(5,3) * t101 + Ifges(5,4) * t116 + Ifges(5,2) * t115 + Ifges(5,6) * t131 - t133 * t110 + t145 * t112;
t111 = Ifges(5,4) * t133 + Ifges(5,2) * t132 + Ifges(5,6) * t145;
t93 = mrSges(5,2) * t103 - mrSges(5,3) * t100 + Ifges(5,1) * t116 + Ifges(5,4) * t115 + Ifges(5,5) * t131 + t132 * t110 - t145 * t111;
t99 = -m(5) * t103 + t115 * mrSges(5,1) - t116 * mrSges(5,2) + t132 * t121 - t133 * t122;
t175 = mrSges(4,1) * t105 - mrSges(4,2) * t106 + Ifges(4,5) * t138 + Ifges(4,6) * t139 + Ifges(4,3) * qJDD(3) + pkin(3) * t99 + pkin(6) * t91 + t153 * t93 + t156 * t92 + (t154 * t127 - t157 * t128) * qJD(1);
t135 = (-mrSges(4,1) * t157 + mrSges(4,2) * t154) * qJD(1);
t141 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t173;
t89 = m(4) * t106 - qJDD(3) * mrSges(4,2) + t139 * mrSges(4,3) - qJD(3) * t141 + t135 * t172 + t91;
t142 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t172;
t95 = m(4) * t105 + qJDD(3) * mrSges(4,1) - t138 * mrSges(4,3) + qJD(3) * t142 - t135 * t173 + t99;
t167 = -t154 * t95 + t157 * t89;
t81 = m(3) * t119 - t160 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t167;
t90 = t153 * t98 + t156 * t97;
t163 = -m(4) * t108 + t139 * mrSges(4,1) - t138 * mrSges(4,2) - t141 * t173 + t142 * t172 - t90;
t85 = m(3) * t118 + qJDD(1) * mrSges(3,1) - t160 * mrSges(3,2) + t163;
t74 = t151 * t81 + t152 * t85;
t83 = t154 * t89 + t157 * t95;
t168 = -t151 * t85 + t152 * t81;
t126 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t154 + Ifges(4,6) * t157) * qJD(1);
t76 = mrSges(4,2) * t108 - mrSges(4,3) * t105 + Ifges(4,1) * t138 + Ifges(4,4) * t139 + Ifges(4,5) * qJDD(3) - pkin(6) * t90 - qJD(3) * t127 + t126 * t172 - t153 * t92 + t156 * t93;
t162 = mrSges(5,1) * t100 - mrSges(5,2) * t101 + Ifges(5,5) * t116 + Ifges(5,6) * t115 + Ifges(5,3) * t131 + t133 * t111 - t132 * t112;
t78 = -mrSges(4,1) * t108 + mrSges(4,3) * t106 + Ifges(4,4) * t138 + Ifges(4,2) * t139 + Ifges(4,6) * qJDD(3) - pkin(3) * t90 + qJD(3) * t128 - t126 * t173 - t162;
t165 = mrSges(3,1) * t118 - mrSges(3,2) * t119 + Ifges(3,3) * qJDD(1) + pkin(2) * t163 + pkin(5) * t167 + t154 * t76 + t157 * t78;
t164 = mrSges(2,1) * t143 - mrSges(2,2) * t144 + Ifges(2,3) * qJDD(1) + pkin(1) * t74 + t165;
t72 = m(2) * t144 - t160 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t168;
t71 = m(2) * t143 + qJDD(1) * mrSges(2,1) - t160 * mrSges(2,2) + t74;
t70 = -mrSges(3,1) * t150 + mrSges(3,3) * t119 + t160 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t83 - t175;
t69 = mrSges(3,2) * t150 - mrSges(3,3) * t118 + Ifges(3,5) * qJDD(1) - t160 * Ifges(3,6) - pkin(5) * t83 - t154 * t78 + t157 * t76;
t68 = -mrSges(2,2) * g(3) - mrSges(2,3) * t143 + Ifges(2,5) * qJDD(1) - t160 * Ifges(2,6) - qJ(2) * t74 - t151 * t70 + t152 * t69;
t67 = Ifges(2,6) * qJDD(1) + t160 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t144 + t151 * t69 + t152 * t70 - pkin(1) * (m(3) * t150 + t83) + qJ(2) * t168;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t158 * t68 - t155 * t67 - pkin(4) * (t155 * t72 + t158 * t71), t68, t69, t76, t93; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t155 * t68 + t158 * t67 + pkin(4) * (-t155 * t71 + t158 * t72), t67, t70, t78, t92; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t164, t164, t165, t175, t162;];
m_new = t1;
