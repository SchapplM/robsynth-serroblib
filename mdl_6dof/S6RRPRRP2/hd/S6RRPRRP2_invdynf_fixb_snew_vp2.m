% Calculate vector of cutting forces with Newton-Euler
% S6RRPRRP2
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 17:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRRP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:22:29
% EndTime: 2019-05-06 17:22:39
% DurationCPUTime: 3.82s
% Computational Cost: add. (47052->204), mult. (108957->262), div. (0->0), fcn. (79913->10), ass. (0->97)
t102 = sin(qJ(2));
t105 = cos(qJ(2));
t107 = qJD(1) ^ 2;
t103 = sin(qJ(1));
t106 = cos(qJ(1));
t121 = t103 * g(1) - t106 * g(2);
t115 = -qJDD(1) * pkin(1) - t121;
t127 = qJD(1) * t102;
t125 = qJD(1) * qJD(2);
t87 = t105 * qJDD(1) - t102 * t125;
t88 = qJD(2) * pkin(2) - qJ(3) * t127;
t97 = t105 ^ 2;
t111 = -t87 * pkin(2) + qJDD(3) + t88 * t127 + (-qJ(3) * t97 - pkin(7)) * t107 + t115;
t100 = sin(qJ(5));
t86 = t102 * qJDD(1) + t105 * t125;
t98 = sin(pkin(10));
t99 = cos(pkin(10));
t71 = -t98 * t86 + t99 * t87;
t81 = (t102 * t99 + t105 * t98) * qJD(1);
t75 = qJD(2) * pkin(3) - t81 * pkin(8);
t80 = (-t102 * t98 + t105 * t99) * qJD(1);
t79 = t80 ^ 2;
t109 = -t71 * pkin(3) - t79 * pkin(8) + t81 * t75 + t111;
t136 = cos(qJ(5));
t101 = sin(qJ(4));
t104 = cos(qJ(4));
t117 = -t106 * g(1) - t103 * g(2);
t83 = -t107 * pkin(1) + qJDD(1) * pkin(7) + t117;
t128 = t102 * t83;
t135 = pkin(2) * t107;
t55 = qJDD(2) * pkin(2) - t86 * qJ(3) - t128 + (qJ(3) * t125 + t102 * t135 - g(3)) * t105;
t120 = -t102 * g(3) + t105 * t83;
t56 = t87 * qJ(3) - qJD(2) * t88 - t97 * t135 + t120;
t118 = -0.2e1 * qJD(3) * t81 + t99 * t55 - t98 * t56;
t72 = t99 * t86 + t98 * t87;
t25 = (qJD(2) * t80 - t72) * pkin(8) + (t80 * t81 + qJDD(2)) * pkin(3) + t118;
t122 = 0.2e1 * qJD(3) * t80 + t98 * t55 + t99 * t56;
t28 = -t79 * pkin(3) + t71 * pkin(8) - qJD(2) * t75 + t122;
t131 = t101 * t25 + t104 * t28;
t65 = -t101 * t81 + t104 * t80;
t66 = t101 * t80 + t104 * t81;
t51 = -t65 * pkin(4) - t66 * pkin(9);
t96 = qJD(2) + qJD(4);
t94 = t96 ^ 2;
t95 = qJDD(2) + qJDD(4);
t21 = -t94 * pkin(4) + t95 * pkin(9) + t65 * t51 + t131;
t40 = -t66 * qJD(4) - t101 * t72 + t104 * t71;
t41 = t65 * qJD(4) + t101 * t71 + t104 * t72;
t23 = (-t65 * t96 - t41) * pkin(9) + (t66 * t96 - t40) * pkin(4) + t109;
t132 = t100 * t23 + t136 * t21;
t39 = qJDD(5) - t40;
t59 = t100 * t66 - t136 * t96;
t60 = t100 * t96 + t136 * t66;
t42 = t59 * pkin(5) - t60 * qJ(6);
t64 = qJD(5) - t65;
t48 = -t64 * mrSges(7,1) + t60 * mrSges(7,2);
t63 = t64 ^ 2;
t124 = m(7) * (-t63 * pkin(5) + t39 * qJ(6) + 0.2e1 * qJD(6) * t64 - t59 * t42 + t132) + t64 * t48 + t39 * mrSges(7,3);
t43 = t59 * mrSges(7,1) - t60 * mrSges(7,3);
t130 = -t59 * mrSges(6,1) - t60 * mrSges(6,2) - t43;
t133 = -mrSges(6,3) - mrSges(7,2);
t30 = t60 * qJD(5) + t100 * t41 - t136 * t95;
t47 = t64 * mrSges(6,1) - t60 * mrSges(6,3);
t12 = m(6) * t132 - t39 * mrSges(6,2) + t130 * t59 + t133 * t30 - t64 * t47 + t124;
t114 = -t100 * t21 + t136 * t23;
t137 = m(7) * (-t39 * pkin(5) - t63 * qJ(6) + t60 * t42 + qJDD(6) - t114);
t31 = -t59 * qJD(5) + t100 * t95 + t136 * t41;
t45 = -t59 * mrSges(7,2) + t64 * mrSges(7,3);
t46 = -t64 * mrSges(6,2) - t59 * mrSges(6,3);
t14 = m(6) * t114 - t137 + (t46 + t45) * t64 + t130 * t60 + (mrSges(6,1) + mrSges(7,1)) * t39 + t133 * t31;
t61 = -t96 * mrSges(5,2) + t65 * mrSges(5,3);
t62 = t96 * mrSges(5,1) - t66 * mrSges(5,3);
t113 = -m(5) * t109 + t40 * mrSges(5,1) - t41 * mrSges(5,2) - t100 * t12 - t136 * t14 + t65 * t61 - t66 * t62;
t73 = -qJD(2) * mrSges(4,2) + t80 * mrSges(4,3);
t74 = qJD(2) * mrSges(4,1) - t81 * mrSges(4,3);
t110 = -m(4) * t111 + t71 * mrSges(4,1) - t72 * mrSges(4,2) + t80 * t73 - t81 * t74 + t113;
t89 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t127;
t126 = qJD(1) * t105;
t90 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t126;
t140 = (t102 * t89 - t105 * t90) * qJD(1) + m(3) * (-t107 * pkin(7) + t115) - t87 * mrSges(3,1) + t86 * mrSges(3,2) - t110;
t119 = -t101 * t28 + t104 * t25;
t20 = -t95 * pkin(4) - t94 * pkin(9) + t66 * t51 - t119;
t123 = m(7) * (-0.2e1 * qJD(6) * t60 + (t59 * t64 - t31) * qJ(6) + (t60 * t64 + t30) * pkin(5) + t20) + t30 * mrSges(7,1) + t59 * t45;
t139 = m(6) * t20 + t30 * mrSges(6,1) + (t47 - t48) * t60 + (mrSges(6,2) - mrSges(7,3)) * t31 + t59 * t46 + t123;
t50 = -t65 * mrSges(5,1) + t66 * mrSges(5,2);
t10 = m(5) * t119 + t95 * mrSges(5,1) - t41 * mrSges(5,3) - t66 * t50 + t96 * t61 - t139;
t69 = -t80 * mrSges(4,1) + t81 * mrSges(4,2);
t9 = m(5) * t131 - t95 * mrSges(5,2) + t40 * mrSges(5,3) - t100 * t14 + t136 * t12 + t65 * t50 - t96 * t62;
t6 = m(4) * t118 + qJDD(2) * mrSges(4,1) - t72 * mrSges(4,3) + qJD(2) * t73 + t104 * t10 + t101 * t9 - t81 * t69;
t7 = m(4) * t122 - qJDD(2) * mrSges(4,2) + t71 * mrSges(4,3) - qJD(2) * t74 - t101 * t10 + t104 * t9 + t80 * t69;
t85 = (-mrSges(3,1) * t105 + mrSges(3,2) * t102) * qJD(1);
t4 = m(3) * (-t105 * g(3) - t128) - t86 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t85 * t127 + qJD(2) * t90 + t98 * t7 + t99 * t6;
t5 = m(3) * t120 - qJDD(2) * mrSges(3,2) + t87 * mrSges(3,3) - qJD(2) * t89 + t85 * t126 - t98 * t6 + t99 * t7;
t138 = t102 * t5 + t105 * t4;
t8 = m(2) * t121 + qJDD(1) * mrSges(2,1) - t107 * mrSges(2,2) - t140;
t1 = m(2) * t117 - t107 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t102 * t4 + t105 * t5;
t2 = [-m(1) * g(1) + t106 * t1 - t103 * t8, t1, t5, t7, t9, t12, -t30 * mrSges(7,2) - t59 * t43 + t124; -m(1) * g(2) + t103 * t1 + t106 * t8, t8, t4, t6, t10, t14, -t31 * mrSges(7,3) - t60 * t48 + t123; (-m(1) - m(2)) * g(3) + t138, -m(2) * g(3) + t138, t140, -t110, -t113, t139, -t39 * mrSges(7,1) + t31 * mrSges(7,2) + t60 * t43 - t64 * t45 + t137;];
f_new  = t2;
