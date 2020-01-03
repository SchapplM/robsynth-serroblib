% Calculate vector of cutting torques with Newton-Euler for
% S4PRRR3
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
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PRRR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR3_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR3_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR3_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR3_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR3_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:36
% EndTime: 2019-12-31 16:31:37
% DurationCPUTime: 0.87s
% Computational Cost: add. (11184->138), mult. (15432->183), div. (0->0), fcn. (9296->8), ass. (0->65)
t124 = qJD(2) + qJD(3);
t129 = sin(qJ(4));
t132 = cos(qJ(4));
t105 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t129 + Ifges(5,2) * t132) * t124;
t106 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t129 + Ifges(5,4) * t132) * t124;
t123 = qJDD(2) + qJDD(3);
t146 = qJD(4) * t124;
t110 = t129 * t123 + t132 * t146;
t111 = t132 * t123 - t129 * t146;
t126 = -g(3) + qJDD(1);
t122 = t124 ^ 2;
t127 = sin(pkin(7));
t128 = cos(pkin(7));
t118 = t127 * g(1) - t128 * g(2);
t119 = -t128 * g(1) - t127 * g(2);
t131 = sin(qJ(2));
t134 = cos(qJ(2));
t103 = t131 * t118 + t134 * t119;
t135 = qJD(2) ^ 2;
t100 = -t135 * pkin(2) + t103;
t130 = sin(qJ(3));
t133 = cos(qJ(3));
t102 = t134 * t118 - t131 * t119;
t99 = qJDD(2) * pkin(2) + t102;
t96 = t133 * t100 + t130 * t99;
t93 = -t122 * pkin(3) + t123 * pkin(6) + t96;
t90 = t132 * t126 - t129 * t93;
t91 = t129 * t126 + t132 * t93;
t149 = mrSges(5,1) * t90 - mrSges(5,2) * t91 + Ifges(5,5) * t110 + Ifges(5,6) * t111 + Ifges(5,3) * qJDD(4) + (t105 * t129 - t106 * t132) * t124;
t109 = (-mrSges(5,1) * t132 + mrSges(5,2) * t129) * t124;
t147 = t124 * t132;
t116 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t147;
t148 = t124 * t129;
t88 = m(5) * t90 + qJDD(4) * mrSges(5,1) - t110 * mrSges(5,3) + qJD(4) * t116 - t109 * t148;
t115 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t148;
t89 = m(5) * t91 - qJDD(4) * mrSges(5,2) + t111 * mrSges(5,3) - qJD(4) * t115 + t109 * t147;
t144 = -t129 * t88 + t132 * t89;
t75 = m(4) * t96 - t122 * mrSges(4,1) - t123 * mrSges(4,2) + t144;
t95 = -t130 * t100 + t133 * t99;
t92 = -t123 * pkin(3) - t122 * pkin(6) - t95;
t139 = -m(5) * t92 + t111 * mrSges(5,1) - t110 * mrSges(5,2) - t115 * t148 + t116 * t147;
t83 = m(4) * t95 + t123 * mrSges(4,1) - t122 * mrSges(4,2) + t139;
t72 = t130 * t75 + t133 * t83;
t69 = m(3) * t102 + qJDD(2) * mrSges(3,1) - t135 * mrSges(3,2) + t72;
t143 = -t130 * t83 + t133 * t75;
t70 = m(3) * t103 - t135 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t143;
t63 = t131 * t70 + t134 * t69;
t77 = t129 * t89 + t132 * t88;
t145 = m(4) * t126 + t77;
t142 = -t131 * t69 + t134 * t70;
t104 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t129 + Ifges(5,6) * t132) * t124;
t80 = -mrSges(5,1) * t92 + mrSges(5,3) * t91 + Ifges(5,4) * t110 + Ifges(5,2) * t111 + Ifges(5,6) * qJDD(4) + qJD(4) * t106 - t104 * t148;
t81 = mrSges(5,2) * t92 - mrSges(5,3) * t90 + Ifges(5,1) * t110 + Ifges(5,4) * t111 + Ifges(5,5) * qJDD(4) - qJD(4) * t105 + t104 * t147;
t140 = mrSges(4,1) * t95 - mrSges(4,2) * t96 + Ifges(4,3) * t123 + pkin(3) * t139 + pkin(6) * t144 + t129 * t81 + t132 * t80;
t137 = mrSges(3,1) * t102 - mrSges(3,2) * t103 + Ifges(3,3) * qJDD(2) + pkin(2) * t72 + t140;
t136 = mrSges(2,1) * t118 - mrSges(2,2) * t119 + pkin(1) * t63 + t137;
t65 = -mrSges(4,1) * t126 + mrSges(4,3) * t96 + t122 * Ifges(4,5) + Ifges(4,6) * t123 - pkin(3) * t77 - t149;
t64 = mrSges(4,2) * t126 - mrSges(4,3) * t95 + Ifges(4,5) * t123 - t122 * Ifges(4,6) - pkin(6) * t77 - t129 * t80 + t132 * t81;
t61 = m(2) * t119 + t142;
t60 = m(2) * t118 + t63;
t59 = mrSges(3,2) * t126 - mrSges(3,3) * t102 + Ifges(3,5) * qJDD(2) - t135 * Ifges(3,6) - pkin(5) * t72 - t130 * t65 + t133 * t64;
t58 = -mrSges(3,1) * t126 + mrSges(3,3) * t103 + t135 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t145 + pkin(5) * t143 + t130 * t64 + t133 * t65;
t57 = mrSges(2,2) * t126 - mrSges(2,3) * t118 - pkin(4) * t63 - t131 * t58 + t134 * t59;
t56 = -mrSges(2,1) * t126 + mrSges(2,3) * t119 + t131 * t59 + t134 * t58 - pkin(1) * (m(3) * t126 + t145) + pkin(4) * t142;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t128 * t57 - t127 * t56 - qJ(1) * (t127 * t61 + t128 * t60), t57, t59, t64, t81; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t127 * t57 + t128 * t56 + qJ(1) * (-t127 * t60 + t128 * t61), t56, t58, t65, t80; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t136, t136, t137, t140, t149;];
m_new = t1;
