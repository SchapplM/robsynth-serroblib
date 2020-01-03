% Calculate vector of cutting torques with Newton-Euler for
% S4PPRR4
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
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PPRR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR4_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR4_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR4_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR4_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:28
% EndTime: 2019-12-31 16:18:30
% DurationCPUTime: 0.70s
% Computational Cost: add. (7040->127), mult. (11563->171), div. (0->0), fcn. (7300->8), ass. (0->61)
t117 = sin(qJ(4));
t119 = cos(qJ(4));
t132 = qJD(3) * qJD(4);
t103 = t117 * qJDD(3) + t119 * t132;
t104 = t119 * qJDD(3) - t117 * t132;
t114 = sin(pkin(6));
t116 = cos(pkin(6));
t107 = t114 * g(1) - t116 * g(2);
t106 = qJDD(2) - t107;
t121 = qJD(3) ^ 2;
t118 = sin(qJ(3));
t120 = cos(qJ(3));
t108 = -t116 * g(1) - t114 * g(2);
t112 = -g(3) + qJDD(1);
t113 = sin(pkin(7));
t115 = cos(pkin(7));
t95 = -t113 * t108 + t115 * t112;
t96 = t115 * t108 + t113 * t112;
t92 = t118 * t95 + t120 * t96;
t89 = -t121 * pkin(3) + qJDD(3) * pkin(5) + t92;
t86 = t119 * t106 - t117 * t89;
t87 = t117 * t106 + t119 * t89;
t98 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t117 + Ifges(5,2) * t119) * qJD(3);
t99 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t117 + Ifges(5,4) * t119) * qJD(3);
t135 = mrSges(5,1) * t86 - mrSges(5,2) * t87 + Ifges(5,5) * t103 + Ifges(5,6) * t104 + Ifges(5,3) * qJDD(4) + (t117 * t98 - t119 * t99) * qJD(3);
t102 = (-mrSges(5,1) * t119 + mrSges(5,2) * t117) * qJD(3);
t133 = qJD(3) * t119;
t110 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t133;
t134 = qJD(3) * t117;
t82 = m(5) * t86 + qJDD(4) * mrSges(5,1) - t103 * mrSges(5,3) + qJD(4) * t110 - t102 * t134;
t109 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t134;
t83 = m(5) * t87 - qJDD(4) * mrSges(5,2) + t104 * mrSges(5,3) - qJD(4) * t109 + t102 * t133;
t130 = -t117 * t82 + t119 * t83;
t67 = m(4) * t92 - t121 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t130;
t91 = -t118 * t96 + t120 * t95;
t88 = -qJDD(3) * pkin(3) - t121 * pkin(5) - t91;
t124 = -m(5) * t88 + t104 * mrSges(5,1) - t103 * mrSges(5,2) - t109 * t134 + t110 * t133;
t78 = m(4) * t91 + qJDD(3) * mrSges(4,1) - t121 * mrSges(4,2) + t124;
t64 = t118 * t67 + t120 * t78;
t71 = t117 * t83 + t119 * t82;
t62 = m(3) * t95 + t64;
t129 = -t118 * t78 + t120 * t67;
t63 = m(3) * t96 + t129;
t131 = -t113 * t62 + t115 * t63;
t127 = (-m(3) - m(4)) * t106 - t71;
t97 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t117 + Ifges(5,6) * t119) * qJD(3);
t75 = -mrSges(5,1) * t88 + mrSges(5,3) * t87 + Ifges(5,4) * t103 + Ifges(5,2) * t104 + Ifges(5,6) * qJDD(4) + qJD(4) * t99 - t97 * t134;
t76 = mrSges(5,2) * t88 - mrSges(5,3) * t86 + Ifges(5,1) * t103 + Ifges(5,4) * t104 + Ifges(5,5) * qJDD(4) - qJD(4) * t98 + t97 * t133;
t59 = mrSges(4,2) * t106 - mrSges(4,3) * t91 + Ifges(4,5) * qJDD(3) - t121 * Ifges(4,6) - pkin(5) * t71 - t117 * t75 + t119 * t76;
t60 = -mrSges(4,1) * t106 + mrSges(4,3) * t92 + t121 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t71 - t135;
t53 = -mrSges(3,1) * t106 + mrSges(3,3) * t96 + t118 * t59 + t120 * t60 - pkin(2) * (m(4) * t106 + t71) + pkin(4) * t129;
t55 = mrSges(3,2) * t106 - mrSges(3,3) * t95 - pkin(4) * t64 - t118 * t60 + t120 * t59;
t126 = mrSges(2,1) * t107 - mrSges(2,2) * t108 + pkin(1) * t127 + qJ(2) * t131 + t113 * t55 + t115 * t53;
t125 = mrSges(4,1) * t91 - mrSges(4,2) * t92 + Ifges(4,3) * qJDD(3) + pkin(3) * t124 + pkin(5) * t130 + t117 * t76 + t119 * t75;
t122 = mrSges(3,1) * t95 - mrSges(3,2) * t96 + pkin(2) * t64 + t125;
t68 = m(2) * t107 + t127;
t58 = t113 * t63 + t115 * t62;
t56 = m(2) * t108 + t131;
t51 = -mrSges(2,1) * t112 + mrSges(2,3) * t108 - pkin(1) * t58 - t122;
t50 = mrSges(2,2) * t112 - mrSges(2,3) * t107 - qJ(2) * t58 - t113 * t53 + t115 * t55;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t116 * t50 - t114 * t51 - qJ(1) * (t114 * t56 + t116 * t68), t50, t55, t59, t76; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t114 * t50 + t116 * t51 + qJ(1) * (-t114 * t68 + t116 * t56), t51, t53, t60, t75; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t126, t126, t122, t125, t135;];
m_new = t1;
