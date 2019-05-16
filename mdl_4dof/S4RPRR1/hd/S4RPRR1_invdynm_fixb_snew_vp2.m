% Calculate vector of cutting torques with Newton-Euler for
% S4RPRR1
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
% Datum: 2019-05-04 19:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RPRR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR1_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR1_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR1_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR1_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR1_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:16:27
% EndTime: 2019-05-04 19:16:28
% DurationCPUTime: 0.90s
% Computational Cost: add. (14321->120), mult. (21234->147), div. (0->0), fcn. (11538->8), ass. (0->58)
t118 = sin(pkin(7));
t119 = cos(pkin(7));
t126 = qJD(1) ^ 2;
t121 = sin(qJ(3));
t124 = cos(qJ(3));
t114 = qJD(1) + qJD(3);
t112 = t114 ^ 2;
t113 = qJDD(1) + qJDD(3);
t120 = sin(qJ(4));
t123 = cos(qJ(4));
t108 = qJD(4) + t114;
t106 = t108 ^ 2;
t107 = qJDD(4) + t113;
t122 = sin(qJ(1));
t125 = cos(qJ(1));
t103 = -g(1) * t125 - g(2) * t122;
t100 = -pkin(1) * t126 + t103;
t102 = g(1) * t122 - g(2) * t125;
t99 = qJDD(1) * pkin(1) + t102;
t94 = -t100 * t118 + t119 * t99;
t91 = qJDD(1) * pkin(2) + t94;
t95 = t100 * t119 + t118 * t99;
t92 = -pkin(2) * t126 + t95;
t86 = -t121 * t92 + t124 * t91;
t83 = pkin(3) * t113 + t86;
t87 = t121 * t91 + t124 * t92;
t84 = -pkin(3) * t112 + t87;
t81 = -t120 * t84 + t123 * t83;
t78 = m(5) * t81 + mrSges(5,1) * t107 - mrSges(5,2) * t106;
t82 = t120 * t83 + t123 * t84;
t79 = m(5) * t82 - mrSges(5,1) * t106 - mrSges(5,2) * t107;
t72 = t120 * t79 + t123 * t78;
t69 = m(4) * t86 + mrSges(4,1) * t113 - mrSges(4,2) * t112 + t72;
t132 = -t120 * t78 + t123 * t79;
t70 = m(4) * t87 - mrSges(4,1) * t112 - mrSges(4,2) * t113 + t132;
t63 = t121 * t70 + t124 * t69;
t60 = m(3) * t94 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t126 + t63;
t131 = -t121 * t69 + t124 * t70;
t61 = m(3) * t95 - mrSges(3,1) * t126 - qJDD(1) * mrSges(3,2) + t131;
t56 = t118 * t61 + t119 * t60;
t117 = -g(3) + qJDD(2);
t134 = (m(4) + m(5)) * t117;
t133 = -t118 * t60 + t119 * t61;
t130 = mrSges(5,1) * t81 - mrSges(5,2) * t82 + Ifges(5,3) * t107;
t129 = mrSges(4,1) * t86 - mrSges(4,2) * t87 + Ifges(4,3) * t113 + pkin(3) * t72 + t130;
t128 = mrSges(3,1) * t94 - mrSges(3,2) * t95 + Ifges(3,3) * qJDD(1) + pkin(2) * t63 + t129;
t127 = mrSges(2,1) * t102 - mrSges(2,2) * t103 + Ifges(2,3) * qJDD(1) + pkin(1) * t56 + t128;
t74 = mrSges(5,2) * t117 - mrSges(5,3) * t81 + Ifges(5,5) * t107 - Ifges(5,6) * t106;
t73 = -mrSges(5,1) * t117 + mrSges(5,3) * t82 + Ifges(5,5) * t106 + Ifges(5,6) * t107;
t65 = mrSges(4,2) * t117 - mrSges(4,3) * t86 + Ifges(4,5) * t113 - Ifges(4,6) * t112 - pkin(6) * t72 - t120 * t73 + t123 * t74;
t64 = Ifges(4,6) * t113 + t112 * Ifges(4,5) + mrSges(4,3) * t87 + t120 * t74 + t123 * t73 + pkin(6) * t132 + (-m(5) * pkin(3) - mrSges(4,1)) * t117;
t54 = m(2) * t103 - mrSges(2,1) * t126 - qJDD(1) * mrSges(2,2) + t133;
t53 = m(2) * t102 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t126 + t56;
t52 = mrSges(3,2) * t117 - mrSges(3,3) * t94 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t126 - pkin(5) * t63 - t121 * t64 + t124 * t65;
t51 = -mrSges(3,1) * t117 + mrSges(3,3) * t95 + t126 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t134 + pkin(5) * t131 + t121 * t65 + t124 * t64;
t50 = -mrSges(2,2) * g(3) - mrSges(2,3) * t102 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t126 - qJ(2) * t56 - t118 * t51 + t119 * t52;
t49 = Ifges(2,6) * qJDD(1) + t126 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t103 + t118 * t52 + t119 * t51 - pkin(1) * (m(3) * t117 + t134) + qJ(2) * t133;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t125 * t50 - t122 * t49 - pkin(4) * (t122 * t54 + t125 * t53), t50, t52, t65, t74; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t122 * t50 + t125 * t49 + pkin(4) * (-t122 * t53 + t125 * t54), t49, t51, t64, t73; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t127, t127, t128, t129, t130;];
m_new  = t1;
