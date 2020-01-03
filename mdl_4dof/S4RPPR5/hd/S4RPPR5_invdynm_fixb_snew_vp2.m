% Calculate vector of cutting torques with Newton-Euler for
% S4RPPR5
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
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RPPR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR5_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR5_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR5_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR5_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR5_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:44
% EndTime: 2019-12-31 16:39:45
% DurationCPUTime: 0.58s
% Computational Cost: add. (5576->149), mult. (9211->185), div. (0->0), fcn. (3300->6), ass. (0->63)
t118 = sin(qJ(4));
t120 = cos(qJ(4));
t113 = g(3) + qJDD(3);
t123 = qJD(1) ^ 2;
t116 = sin(pkin(6));
t117 = cos(pkin(6));
t119 = sin(qJ(1));
t121 = cos(qJ(1));
t103 = -t121 * g(1) - t119 * g(2);
t132 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t103;
t139 = -pkin(1) - pkin(2);
t84 = t139 * t123 + t132;
t102 = t119 * g(1) - t121 * g(2);
t131 = -t123 * qJ(2) + qJDD(2) - t102;
t87 = t139 * qJDD(1) + t131;
t81 = t116 * t87 + t117 * t84;
t78 = -t123 * pkin(3) - qJDD(1) * pkin(5) + t81;
t75 = t120 * t113 - t118 * t78;
t76 = t118 * t113 + t120 * t78;
t92 = Ifges(5,6) * qJD(4) + (-Ifges(5,4) * t118 - Ifges(5,2) * t120) * qJD(1);
t93 = Ifges(5,5) * qJD(4) + (-Ifges(5,1) * t118 - Ifges(5,4) * t120) * qJD(1);
t135 = qJD(1) * qJD(4);
t97 = -t118 * qJDD(1) - t120 * t135;
t98 = -t120 * qJDD(1) + t118 * t135;
t140 = mrSges(5,1) * t75 - mrSges(5,2) * t76 + Ifges(5,5) * t97 + Ifges(5,6) * t98 + Ifges(5,3) * qJDD(4) - (t118 * t92 - t120 * t93) * qJD(1);
t138 = mrSges(2,1) + mrSges(3,1);
t137 = qJD(1) * t118;
t136 = qJD(1) * t120;
t101 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t136;
t96 = (mrSges(5,1) * t120 - mrSges(5,2) * t118) * qJD(1);
t72 = m(5) * t75 + qJDD(4) * mrSges(5,1) - t97 * mrSges(5,3) + qJD(4) * t101 + t96 * t137;
t100 = qJD(4) * mrSges(5,1) + mrSges(5,3) * t137;
t73 = m(5) * t76 - qJDD(4) * mrSges(5,2) + t98 * mrSges(5,3) - qJD(4) * t100 - t96 * t136;
t67 = -t118 * t72 + t120 * t73;
t63 = m(4) * t81 - t123 * mrSges(4,1) + qJDD(1) * mrSges(4,2) + t67;
t80 = -t116 * t84 + t117 * t87;
t77 = qJDD(1) * pkin(3) - t123 * pkin(5) - t80;
t74 = -m(5) * t77 + t98 * mrSges(5,1) - t97 * mrSges(5,2) + t100 * t137 - t101 * t136;
t70 = m(4) * t80 - qJDD(1) * mrSges(4,1) - t123 * mrSges(4,2) + t74;
t61 = -t116 * t70 + t117 * t63;
t60 = t116 * t63 + t117 * t70;
t66 = t118 * t73 + t120 * t72;
t88 = -t123 * pkin(1) + t132;
t133 = m(3) * t88 + qJDD(1) * mrSges(3,3) + t61;
t65 = -m(4) * t113 - t66;
t90 = -qJDD(1) * pkin(1) + t131;
t130 = -m(3) * t90 + qJDD(1) * mrSges(3,1) + t123 * mrSges(3,3) - t60;
t91 = Ifges(5,3) * qJD(4) + (-Ifges(5,5) * t118 - Ifges(5,6) * t120) * qJD(1);
t68 = -mrSges(5,1) * t77 + mrSges(5,3) * t76 + Ifges(5,4) * t97 + Ifges(5,2) * t98 + Ifges(5,6) * qJDD(4) + qJD(4) * t93 + t91 * t137;
t69 = mrSges(5,2) * t77 - mrSges(5,3) * t75 + Ifges(5,1) * t97 + Ifges(5,4) * t98 + Ifges(5,5) * qJDD(4) - qJD(4) * t92 - t91 * t136;
t54 = mrSges(4,2) * t113 - mrSges(4,3) * t80 - Ifges(4,5) * qJDD(1) - t123 * Ifges(4,6) - pkin(5) * t66 - t118 * t68 + t120 * t69;
t57 = -mrSges(4,1) * t113 + mrSges(4,3) * t81 + t123 * Ifges(4,5) - Ifges(4,6) * qJDD(1) - pkin(3) * t66 - t140;
t129 = mrSges(3,2) * t90 + mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) + t123 * Ifges(3,6) - qJ(3) * t60 - t116 * t57 + t117 * t54;
t128 = mrSges(3,2) * t88 - pkin(2) * t65 - qJ(3) * t61 - t116 * t54 - t117 * t57;
t126 = mrSges(4,1) * t80 - mrSges(4,2) * t81 - Ifges(4,3) * qJDD(1) + pkin(3) * t74 + pkin(5) * t67 + t118 * t69 + t120 * t68;
t125 = -mrSges(3,1) * t90 + mrSges(3,3) * t88 + Ifges(3,2) * qJDD(1) - pkin(2) * t60 - t126;
t124 = -mrSges(2,2) * t103 + Ifges(2,3) * qJDD(1) + t125 + qJ(2) * (-t123 * mrSges(3,1) + t133) + pkin(1) * t130 + mrSges(2,1) * t102;
t64 = -m(3) * g(3) + t65;
t56 = m(2) * t102 + qJDD(1) * mrSges(2,1) - t123 * mrSges(2,2) + t130;
t55 = m(2) * t103 - qJDD(1) * mrSges(2,2) - t138 * t123 + t133;
t52 = -mrSges(2,2) * g(3) - mrSges(2,3) * t102 + Ifges(2,5) * qJDD(1) - t123 * Ifges(2,6) - qJ(2) * t64 + t129;
t51 = mrSges(2,3) * t103 - pkin(1) * t64 + (Ifges(3,4) + Ifges(2,5)) * t123 + (Ifges(2,6) - Ifges(3,6)) * qJDD(1) + t138 * g(3) + t128;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t121 * t52 - t119 * t51 - pkin(4) * (t119 * t55 + t121 * t56), t52, t129, t54, t69; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t119 * t52 + t121 * t51 + pkin(4) * (-t119 * t56 + t121 * t55), t51, t125, t57, t68; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t124, t124, -mrSges(3,1) * g(3) - t123 * Ifges(3,4) + Ifges(3,6) * qJDD(1) - t128, t126, t140;];
m_new = t1;
