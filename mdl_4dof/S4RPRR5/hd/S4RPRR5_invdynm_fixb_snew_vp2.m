% Calculate vector of cutting torques with Newton-Euler for
% S4RPRR5
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
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RPRR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR5_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR5_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR5_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR5_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR5_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:34
% EndTime: 2019-12-31 16:51:35
% DurationCPUTime: 0.63s
% Computational Cost: add. (7064->151), mult. (9211->186), div. (0->0), fcn. (3300->6), ass. (0->64)
t110 = -qJDD(1) + qJDD(3);
t120 = sin(qJ(4));
t123 = cos(qJ(4));
t111 = -qJD(1) + qJD(3);
t139 = qJD(4) * t111;
t100 = t120 * t110 + t123 * t139;
t101 = t123 * t110 - t120 * t139;
t109 = t111 ^ 2;
t121 = sin(qJ(3));
t124 = cos(qJ(3));
t127 = qJD(1) ^ 2;
t122 = sin(qJ(1));
t125 = cos(qJ(1));
t106 = -t125 * g(1) - t122 * g(2);
t136 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t106;
t143 = -pkin(1) - pkin(2);
t87 = t143 * t127 + t136;
t105 = t122 * g(1) - t125 * g(2);
t135 = -t127 * qJ(2) + qJDD(2) - t105;
t90 = t143 * qJDD(1) + t135;
t84 = t121 * t90 + t124 * t87;
t81 = -t109 * pkin(3) + t110 * pkin(6) + t84;
t78 = t123 * g(3) - t120 * t81;
t79 = t120 * g(3) + t123 * t81;
t92 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t120 + Ifges(5,2) * t123) * t111;
t93 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t120 + Ifges(5,4) * t123) * t111;
t144 = mrSges(5,1) * t78 - mrSges(5,2) * t79 + Ifges(5,5) * t100 + Ifges(5,6) * t101 + Ifges(5,3) * qJDD(4) + (t120 * t92 - t123 * t93) * t111;
t142 = mrSges(2,1) + mrSges(3,1);
t141 = t111 * t120;
t140 = t111 * t123;
t103 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t140;
t99 = (-mrSges(5,1) * t123 + mrSges(5,2) * t120) * t111;
t76 = m(5) * t78 + qJDD(4) * mrSges(5,1) - t100 * mrSges(5,3) + qJD(4) * t103 - t99 * t141;
t102 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t141;
t77 = m(5) * t79 - qJDD(4) * mrSges(5,2) + t101 * mrSges(5,3) - qJD(4) * t102 + t99 * t140;
t70 = -t120 * t76 + t123 * t77;
t66 = m(4) * t84 - t109 * mrSges(4,1) - t110 * mrSges(4,2) + t70;
t83 = -t121 * t87 + t124 * t90;
t80 = -t110 * pkin(3) - t109 * pkin(6) - t83;
t74 = -m(5) * t80 + t101 * mrSges(5,1) - t100 * mrSges(5,2) - t102 * t141 + t103 * t140;
t73 = m(4) * t83 + t110 * mrSges(4,1) - t109 * mrSges(4,2) + t74;
t64 = -t121 * t73 + t124 * t66;
t69 = t120 * t77 + t123 * t76;
t63 = t121 * t66 + t124 * t73;
t94 = -t127 * pkin(1) + t136;
t137 = m(3) * t94 + qJDD(1) * mrSges(3,3) + t64;
t98 = -qJDD(1) * pkin(1) + t135;
t134 = -m(3) * t98 + qJDD(1) * mrSges(3,1) + t127 * mrSges(3,3) - t63;
t91 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t120 + Ifges(5,6) * t123) * t111;
t71 = -mrSges(5,1) * t80 + mrSges(5,3) * t79 + Ifges(5,4) * t100 + Ifges(5,2) * t101 + Ifges(5,6) * qJDD(4) + qJD(4) * t93 - t91 * t141;
t72 = mrSges(5,2) * t80 - mrSges(5,3) * t78 + Ifges(5,1) * t100 + Ifges(5,4) * t101 + Ifges(5,5) * qJDD(4) - qJD(4) * t92 + t91 * t140;
t57 = mrSges(4,2) * g(3) - mrSges(4,3) * t83 + Ifges(4,5) * t110 - t109 * Ifges(4,6) - pkin(6) * t69 - t120 * t71 + t123 * t72;
t60 = -mrSges(4,1) * g(3) + mrSges(4,3) * t84 + t109 * Ifges(4,5) + Ifges(4,6) * t110 - pkin(3) * t69 - t144;
t133 = mrSges(3,2) * t98 + mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) + t127 * Ifges(3,6) - pkin(5) * t63 - t121 * t60 + t124 * t57;
t132 = mrSges(3,2) * t94 - pkin(2) * (-m(4) * g(3) - t69) - pkin(5) * t64 - t121 * t57 - t124 * t60;
t130 = mrSges(4,1) * t83 - mrSges(4,2) * t84 + Ifges(4,3) * t110 + pkin(3) * t74 + pkin(6) * t70 + t120 * t72 + t123 * t71;
t129 = -mrSges(3,1) * t98 + mrSges(3,3) * t94 + Ifges(3,2) * qJDD(1) - pkin(2) * t63 - t130;
t128 = -mrSges(2,2) * t106 + mrSges(2,1) * t105 + Ifges(2,3) * qJDD(1) + t129 + qJ(2) * (-t127 * mrSges(3,1) + t137) + pkin(1) * t134;
t67 = (-m(3) - m(4)) * g(3) - t69;
t59 = m(2) * t105 + qJDD(1) * mrSges(2,1) - t127 * mrSges(2,2) + t134;
t58 = m(2) * t106 - qJDD(1) * mrSges(2,2) - t142 * t127 + t137;
t55 = -mrSges(2,2) * g(3) - mrSges(2,3) * t105 + Ifges(2,5) * qJDD(1) - t127 * Ifges(2,6) - qJ(2) * t67 + t133;
t54 = mrSges(2,3) * t106 - pkin(1) * t67 + (Ifges(3,4) + Ifges(2,5)) * t127 + (Ifges(2,6) - Ifges(3,6)) * qJDD(1) + t142 * g(3) + t132;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t125 * t55 - t122 * t54 - pkin(4) * (t122 * t58 + t125 * t59), t55, t133, t57, t72; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t122 * t55 + t125 * t54 + pkin(4) * (-t122 * t59 + t125 * t58), t54, t129, t60, t71; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t128, t128, -mrSges(3,1) * g(3) - t127 * Ifges(3,4) + Ifges(3,6) * qJDD(1) - t132, t130, t144;];
m_new = t1;
