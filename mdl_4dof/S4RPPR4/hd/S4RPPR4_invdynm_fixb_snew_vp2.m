% Calculate vector of cutting torques with Newton-Euler for
% S4RPPR4
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
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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

function m_new = S4RPPR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR4_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR4_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR4_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR4_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR4_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:50
% EndTime: 2019-12-31 16:38:50
% DurationCPUTime: 0.51s
% Computational Cost: add. (4635->148), mult. (7778->182), div. (0->0), fcn. (3398->6), ass. (0->64)
t148 = -pkin(2) - pkin(5);
t147 = mrSges(3,1) - mrSges(4,2);
t146 = -Ifges(4,4) + Ifges(3,5);
t145 = Ifges(4,5) - Ifges(3,6);
t124 = sin(pkin(6));
t125 = cos(pkin(6));
t130 = qJD(1) ^ 2;
t126 = sin(qJ(4));
t128 = cos(qJ(4));
t105 = (mrSges(5,1) * t126 + mrSges(5,2) * t128) * qJD(1);
t142 = qJD(1) * qJD(4);
t108 = t128 * qJDD(1) - t126 * t142;
t144 = qJD(1) * t126;
t110 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t144;
t143 = qJD(1) * t128;
t121 = -g(3) + qJDD(2);
t127 = sin(qJ(1));
t129 = cos(qJ(1));
t112 = t127 * g(1) - t129 * g(2);
t104 = qJDD(1) * pkin(1) + t112;
t113 = -t129 * g(1) - t127 * g(2);
t106 = -t130 * pkin(1) + t113;
t89 = t125 * t104 - t124 * t106;
t140 = -t130 * qJ(3) + qJDD(3) - t89;
t84 = t148 * qJDD(1) + t140;
t80 = -t126 * t121 + t128 * t84;
t77 = m(5) * t80 + qJDD(4) * mrSges(5,1) - t108 * mrSges(5,3) + qJD(4) * t110 - t105 * t143;
t107 = -t126 * qJDD(1) - t128 * t142;
t111 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t143;
t81 = t128 * t121 + t126 * t84;
t78 = m(5) * t81 - qJDD(4) * mrSges(5,2) + t107 * mrSges(5,3) - qJD(4) * t111 - t105 * t144;
t66 = t126 * t78 + t128 * t77;
t87 = -qJDD(1) * pkin(2) + t140;
t137 = -m(4) * t87 + t130 * mrSges(4,3) - t66;
t63 = m(3) * t89 - t130 * mrSges(3,2) + t147 * qJDD(1) + t137;
t90 = t124 * t104 + t125 * t106;
t139 = qJDD(1) * qJ(3) + 0.2e1 * qJD(3) * qJD(1) + t90;
t83 = t148 * t130 + t139;
t76 = -m(5) * t83 + t107 * mrSges(5,1) - t108 * mrSges(5,2) - t110 * t144 - t111 * t143;
t85 = t130 * pkin(2) - t139;
t135 = -m(4) * t85 + t130 * mrSges(4,2) + qJDD(1) * mrSges(4,3) - t76;
t70 = m(3) * t90 - t130 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t135;
t61 = t124 * t70 + t125 * t63;
t141 = -t124 * t63 + t125 * t70;
t67 = -t126 * t77 + t128 * t78;
t65 = m(4) * t121 + t67;
t97 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t128 - Ifges(5,2) * t126) * qJD(1);
t98 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t128 - Ifges(5,4) * t126) * qJD(1);
t138 = mrSges(5,1) * t80 - mrSges(5,2) * t81 + Ifges(5,5) * t108 + Ifges(5,6) * t107 + Ifges(5,3) * qJDD(4) + t97 * t143 + t98 * t144;
t96 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t128 - Ifges(5,6) * t126) * qJD(1);
t72 = -mrSges(5,1) * t83 + mrSges(5,3) * t81 + Ifges(5,4) * t108 + Ifges(5,2) * t107 + Ifges(5,6) * qJDD(4) + qJD(4) * t98 - t96 * t143;
t73 = mrSges(5,2) * t83 - mrSges(5,3) * t80 + Ifges(5,1) * t108 + Ifges(5,4) * t107 + Ifges(5,5) * qJDD(4) - qJD(4) * t97 - t96 * t144;
t136 = mrSges(4,2) * t87 - mrSges(4,3) * t85 + Ifges(4,1) * qJDD(1) - pkin(5) * t66 - t126 * t72 + t128 * t73;
t134 = -mrSges(4,1) * t85 - pkin(3) * t76 - pkin(5) * t67 - t126 * t73 - t128 * t72;
t133 = mrSges(4,1) * t87 + pkin(3) * t66 + t138;
t132 = -mrSges(3,2) * t90 + Ifges(3,3) * qJDD(1) + t136 + pkin(2) * (-qJDD(1) * mrSges(4,2) + t137) + qJ(3) * t135 + mrSges(3,1) * t89;
t131 = mrSges(2,1) * t112 - mrSges(2,2) * t113 + Ifges(2,3) * qJDD(1) + pkin(1) * t61 + t132;
t59 = m(2) * t113 - t130 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t141;
t58 = m(2) * t112 + qJDD(1) * mrSges(2,1) - t130 * mrSges(2,2) + t61;
t57 = -mrSges(3,3) * t89 - qJ(3) * t65 + t145 * t130 + (mrSges(3,2) - mrSges(4,3)) * t121 + t146 * qJDD(1) + t133;
t56 = mrSges(3,3) * t90 - pkin(2) * t65 - t145 * qJDD(1) - t147 * t121 + t146 * t130 + t134;
t55 = -mrSges(2,2) * g(3) - mrSges(2,3) * t112 + Ifges(2,5) * qJDD(1) - t130 * Ifges(2,6) - qJ(2) * t61 - t124 * t56 + t125 * t57;
t54 = Ifges(2,6) * qJDD(1) + t130 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t113 + t124 * t57 + t125 * t56 - pkin(1) * (m(3) * t121 + t65) + qJ(2) * t141;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t129 * t55 - t127 * t54 - pkin(4) * (t127 * t59 + t129 * t58), t55, t57, t136, t73; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t127 * t55 + t129 * t54 + pkin(4) * (-t127 * t58 + t129 * t59), t54, t56, mrSges(4,3) * t121 + Ifges(4,4) * qJDD(1) - t130 * Ifges(4,5) - t133, t72; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t131, t131, t132, -mrSges(4,2) * t121 + t130 * Ifges(4,4) + Ifges(4,5) * qJDD(1) - t134, t138;];
m_new = t1;
