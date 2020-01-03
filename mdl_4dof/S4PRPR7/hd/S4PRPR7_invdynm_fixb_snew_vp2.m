% Calculate vector of cutting torques with Newton-Euler for
% S4PRPR7
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
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PRPR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR7_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR7_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR7_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR7_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR7_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:38
% EndTime: 2019-12-31 16:25:39
% DurationCPUTime: 0.46s
% Computational Cost: add. (3808->139), mult. (6297->171), div. (0->0), fcn. (3026->6), ass. (0->63)
t144 = m(3) + m(4);
t143 = -pkin(2) - pkin(5);
t142 = mrSges(3,1) - mrSges(4,2);
t141 = -Ifges(4,4) + Ifges(3,5);
t140 = Ifges(4,5) - Ifges(3,6);
t121 = sin(qJ(4));
t123 = cos(qJ(4));
t104 = (mrSges(5,1) * t121 + mrSges(5,2) * t123) * qJD(2);
t136 = qJD(2) * qJD(4);
t106 = t123 * qJDD(2) - t121 * t136;
t138 = qJD(2) * t121;
t110 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t138;
t137 = qJD(2) * t123;
t120 = sin(pkin(6));
t139 = cos(pkin(6));
t108 = t120 * g(1) - t139 * g(2);
t125 = qJD(2) ^ 2;
t109 = -t139 * g(1) - t120 * g(2);
t117 = -g(3) + qJDD(1);
t122 = sin(qJ(2));
t124 = cos(qJ(2));
t90 = -t122 * t109 + t124 * t117;
t133 = -t125 * qJ(3) + qJDD(3) - t90;
t86 = t143 * qJDD(2) + t133;
t82 = t121 * t108 + t123 * t86;
t78 = m(5) * t82 + qJDD(4) * mrSges(5,1) - t106 * mrSges(5,3) + qJD(4) * t110 - t104 * t137;
t105 = -t121 * qJDD(2) - t123 * t136;
t111 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t137;
t83 = -t123 * t108 + t121 * t86;
t79 = m(5) * t83 - qJDD(4) * mrSges(5,2) + t105 * mrSges(5,3) - qJD(4) * t111 - t104 * t138;
t69 = -t121 * t78 + t123 * t79;
t91 = t124 * t109 + t122 * t117;
t68 = t121 * t79 + t123 * t78;
t89 = -qJDD(2) * pkin(2) + t133;
t130 = -m(4) * t89 + t125 * mrSges(4,3) - t68;
t63 = m(3) * t90 - t125 * mrSges(3,2) + t142 * qJDD(2) + t130;
t132 = qJDD(2) * qJ(3) + 0.2e1 * qJD(3) * qJD(2) + t91;
t85 = t143 * t125 + t132;
t80 = -m(5) * t85 + t105 * mrSges(5,1) - t106 * mrSges(5,2) - t110 * t138 - t111 * t137;
t87 = t125 * pkin(2) - t132;
t75 = -m(4) * t87 + t125 * mrSges(4,2) + qJDD(2) * mrSges(4,3) - t80;
t71 = m(3) * t91 - t125 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t75;
t135 = -t122 * t63 + t124 * t71;
t94 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t123 - Ifges(5,6) * t121) * qJD(2);
t96 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t123 - Ifges(5,4) * t121) * qJD(2);
t73 = -mrSges(5,1) * t85 + mrSges(5,3) * t83 + Ifges(5,4) * t106 + Ifges(5,2) * t105 + Ifges(5,6) * qJDD(4) + qJD(4) * t96 - t94 * t137;
t95 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t123 - Ifges(5,2) * t121) * qJD(2);
t74 = mrSges(5,2) * t85 - mrSges(5,3) * t82 + Ifges(5,1) * t106 + Ifges(5,4) * t105 + Ifges(5,5) * qJDD(4) - qJD(4) * t95 - t94 * t138;
t128 = -mrSges(4,1) * t87 - pkin(3) * t80 - pkin(5) * t69 - t121 * t74 - t123 * t73;
t67 = -m(4) * t108 + t69;
t57 = mrSges(3,3) * t91 - pkin(2) * t67 - t140 * qJDD(2) + t142 * t108 + t141 * t125 + t128;
t131 = mrSges(5,1) * t82 - mrSges(5,2) * t83 + Ifges(5,5) * t106 + Ifges(5,6) * t105 + Ifges(5,3) * qJDD(4) + t95 * t137 + t96 * t138;
t127 = mrSges(4,1) * t89 + pkin(3) * t68 + t131;
t59 = -mrSges(3,3) * t90 - qJ(3) * t67 + t140 * t125 + (-mrSges(3,2) + mrSges(4,3)) * t108 + t141 * qJDD(2) + t127;
t134 = -mrSges(2,2) * t109 + pkin(4) * t135 + t122 * t59 + t124 * t57 + mrSges(2,1) * t108 + pkin(1) * (t144 * t108 - t69);
t129 = mrSges(4,2) * t89 - mrSges(4,3) * t87 + Ifges(4,1) * qJDD(2) - pkin(5) * t68 - t121 * t73 + t123 * t74;
t126 = mrSges(3,1) * t90 - mrSges(3,2) * t91 + Ifges(3,3) * qJDD(2) + pkin(2) * (-qJDD(2) * mrSges(4,2) + t130) + qJ(3) * t75 + t129;
t65 = (m(2) + t144) * t108 - t69;
t62 = t122 * t71 + t124 * t63;
t60 = m(2) * t109 + t135;
t55 = -mrSges(2,1) * t117 + mrSges(2,3) * t109 - pkin(1) * t62 - t126;
t54 = mrSges(2,2) * t117 - mrSges(2,3) * t108 - pkin(4) * t62 - t122 * t57 + t124 * t59;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t139 * t54 - t120 * t55 - qJ(1) * (t120 * t60 + t139 * t65), t54, t59, t129, t74; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t120 * t54 + t139 * t55 + qJ(1) * (-t120 * t65 + t139 * t60), t55, t57, -mrSges(4,3) * t108 + Ifges(4,4) * qJDD(2) - t125 * Ifges(4,5) - t127, t73; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t134, t134, t126, mrSges(4,2) * t108 + t125 * Ifges(4,4) + Ifges(4,5) * qJDD(2) - t128, t131;];
m_new = t1;
