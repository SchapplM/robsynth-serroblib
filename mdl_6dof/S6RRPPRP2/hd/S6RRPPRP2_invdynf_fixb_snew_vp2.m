% Calculate vector of cutting forces with Newton-Euler
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x7]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 09:13
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPPRP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:08:59
% EndTime: 2019-05-06 09:09:03
% DurationCPUTime: 1.68s
% Computational Cost: add. (15119->203), mult. (34859->247), div. (0->0), fcn. (22852->8), ass. (0->92)
t151 = -2 * qJD(4);
t101 = cos(qJ(2));
t104 = qJD(1) ^ 2;
t102 = cos(qJ(1));
t99 = sin(qJ(1));
t121 = t99 * g(1) - t102 * g(2);
t114 = -qJDD(1) * pkin(1) - t121;
t98 = sin(qJ(2));
t130 = qJD(1) * t98;
t126 = qJD(1) * qJD(2);
t87 = t101 * qJDD(1) - t98 * t126;
t88 = qJD(2) * pkin(2) - qJ(3) * t130;
t95 = t101 ^ 2;
t107 = -t87 * pkin(2) + qJDD(3) + t88 * t130 + (-qJ(3) * t95 - pkin(7)) * t104 + t114;
t100 = cos(qJ(5));
t127 = qJD(1) * t101;
t131 = cos(pkin(9));
t96 = sin(pkin(9));
t79 = -t131 * t127 + t96 * t130;
t129 = qJD(2) * t79;
t86 = t98 * qJDD(1) + t101 * t126;
t62 = t131 * t86 + t96 * t87;
t80 = (t101 * t96 + t131 * t98) * qJD(1);
t105 = (-t62 + t129) * qJ(4) + t107 + (pkin(3) * qJD(2) + t151) * t80;
t103 = qJD(2) ^ 2;
t118 = -t102 * g(1) - t99 * g(2);
t83 = -t104 * pkin(1) + qJDD(1) * pkin(7) + t118;
t140 = t98 * t83;
t143 = pkin(2) * t104;
t39 = qJDD(2) * pkin(2) - t86 * qJ(3) - t140 + (qJ(3) * t126 + t98 * t143 - g(3)) * t101;
t120 = -t98 * g(3) + t101 * t83;
t40 = t87 * qJ(3) - qJD(2) * t88 - t95 * t143 + t120;
t117 = t131 * t39 - t96 * t40;
t54 = t79 * pkin(3) - t80 * qJ(4);
t24 = -qJDD(2) * pkin(3) - t103 * qJ(4) + qJDD(4) - t117 + ((2 * qJD(3)) + t54) * t80;
t18 = (t79 * t80 - qJDD(2)) * pkin(8) + (t62 + t129) * pkin(4) + t24;
t61 = -t131 * t87 + t96 * t86;
t71 = t80 * pkin(4) - qJD(2) * pkin(8);
t78 = t79 ^ 2;
t22 = -t78 * pkin(4) - t80 * t71 + (pkin(3) + pkin(8)) * t61 + t105;
t97 = sin(qJ(5));
t136 = t100 * t22 + t97 * t18;
t65 = t100 * qJD(2) + t97 * t79;
t33 = -t65 * qJD(5) - t97 * qJDD(2) + t100 * t61;
t64 = -t97 * qJD(2) + t100 * t79;
t41 = -t64 * mrSges(7,1) + t65 * mrSges(7,2);
t77 = qJD(5) + t80;
t49 = t77 * pkin(5) - t65 * qJ(6);
t63 = t64 ^ 2;
t123 = m(7) * (-t63 * pkin(5) + t33 * qJ(6) + 0.2e1 * qJD(6) * t64 - t77 * t49 + t136) + t64 * t41 + t33 * mrSges(7,3);
t42 = -t64 * mrSges(6,1) + t65 * mrSges(6,2);
t50 = t77 * mrSges(7,1) - t65 * mrSges(7,3);
t51 = t77 * mrSges(6,1) - t65 * mrSges(6,3);
t60 = qJDD(5) + t62;
t11 = m(6) * t136 + t33 * mrSges(6,3) + t64 * t42 + (-t51 - t50) * t77 + (-mrSges(6,2) - mrSges(7,2)) * t60 + t123;
t70 = t80 * mrSges(5,1) + qJD(2) * mrSges(5,2);
t119 = t100 * t18 - t97 * t22;
t34 = t64 * qJD(5) + t100 * qJDD(2) + t97 * t61;
t47 = -t77 * mrSges(7,2) + t64 * mrSges(7,3);
t124 = m(7) * (-0.2e1 * qJD(6) * t65 + (t64 * t77 - t34) * qJ(6) + (t64 * t65 + t60) * pkin(5) + t119) + t77 * t47 + t60 * mrSges(7,1);
t48 = -t77 * mrSges(6,2) + t64 * mrSges(6,3);
t9 = m(6) * t119 + t60 * mrSges(6,1) + t77 * t48 + (-t42 - t41) * t65 + (-mrSges(6,3) - mrSges(7,3)) * t34 + t124;
t113 = -t100 * t11 + t97 * t9 - m(5) * (t61 * pkin(3) + t105) + t80 * t70 + t62 * mrSges(5,3);
t69 = t79 * mrSges(5,1) - qJD(2) * mrSges(5,3);
t132 = -qJD(2) * mrSges(4,2) - t79 * mrSges(4,3) - t69;
t139 = mrSges(4,1) - mrSges(5,2);
t68 = qJD(2) * mrSges(4,1) - t80 * mrSges(4,3);
t106 = m(4) * t107 + t62 * mrSges(4,2) + t132 * t79 + t139 * t61 + t80 * t68 - t113;
t89 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t130;
t90 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t127;
t150 = t106 - (t101 * t90 - t98 * t89) * qJD(1) - t87 * mrSges(3,1) + t86 * mrSges(3,2) + m(3) * (-t104 * pkin(7) + t114);
t135 = t131 * t40 + t96 * t39;
t148 = t103 * pkin(3) - qJDD(2) * qJ(4) + qJD(2) * t151 + t79 * t54 - t135;
t128 = qJD(3) * t79;
t74 = -0.2e1 * t128;
t108 = -t61 * pkin(4) - t78 * pkin(8) + qJD(2) * t71 - t148 + t74;
t122 = m(7) * (-t33 * pkin(5) - t63 * qJ(6) + t65 * t49 + qJDD(6) + t108) + t65 * t50 + t34 * mrSges(7,2);
t109 = m(6) * t108 + t34 * mrSges(6,2) + (-t48 - t47) * t64 - (mrSges(6,1) + mrSges(7,1)) * t33 + t65 * t51 + t122;
t149 = t109 - m(5) * (0.2e1 * t128 + t148);
t111 = -m(5) * t24 - t100 * t9 - t97 * t11;
t56 = -t79 * mrSges(5,2) - t80 * mrSges(5,3);
t133 = -t79 * mrSges(4,1) - t80 * mrSges(4,2) - t56;
t137 = -mrSges(4,3) - mrSges(5,1);
t7 = m(4) * t117 + (-0.2e1 * m(4) * qJD(3) + t133) * t80 + t137 * t62 + t139 * qJDD(2) + t132 * qJD(2) + t111;
t8 = m(4) * (t74 + t135) + t133 * t79 + t137 * t61 + (-mrSges(4,2) + mrSges(5,3)) * qJDD(2) + (-t68 + t70) * qJD(2) + t149;
t85 = (-mrSges(3,1) * t101 + mrSges(3,2) * t98) * qJD(1);
t4 = m(3) * (-t101 * g(3) - t140) - t86 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t85 * t130 + qJD(2) * t90 + t96 * t8 + t131 * t7;
t5 = m(3) * t120 - qJDD(2) * mrSges(3,2) + t87 * mrSges(3,3) - qJD(2) * t89 + t85 * t127 + t131 * t8 - t96 * t7;
t146 = t101 * t4 + t98 * t5;
t6 = m(2) * t121 + qJDD(1) * mrSges(2,1) - t104 * mrSges(2,2) - t150;
t1 = m(2) * t118 - t104 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t101 * t5 - t98 * t4;
t2 = [-m(1) * g(1) + t102 * t1 - t99 * t6, t1, t5, t8, -t61 * mrSges(5,2) - t79 * t69 - t113, t11, -t60 * mrSges(7,2) - t77 * t50 + t123; -m(1) * g(2) + t99 * t1 + t102 * t6, t6, t4, t7, t61 * mrSges(5,1) - qJDD(2) * mrSges(5,3) - qJD(2) * t70 + t79 * t56 - t149, t9, -t34 * mrSges(7,3) - t65 * t41 + t124; (-m(1) - m(2)) * g(3) + t146, -m(2) * g(3) + t146, t150, t106, t62 * mrSges(5,1) + qJDD(2) * mrSges(5,2) + qJD(2) * t69 + t80 * t56 - t111, t109, -t33 * mrSges(7,1) - t64 * t47 + t122;];
f_new  = t2;
