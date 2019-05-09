% Calculate vector of cutting forces with Newton-Euler
% S6RRRRPR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-07 19:46
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRPR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 19:41:05
% EndTime: 2019-05-07 19:41:19
% DurationCPUTime: 6.15s
% Computational Cost: add. (111069->207), mult. (243382->276), div. (0->0), fcn. (184329->12), ass. (0->106)
t110 = sin(qJ(2));
t114 = cos(qJ(2));
t116 = qJD(1) ^ 2;
t104 = t114 ^ 2;
t111 = sin(qJ(1));
t115 = cos(qJ(1));
t130 = t111 * g(1) - t115 * g(2);
t124 = -qJDD(1) * pkin(1) - t130;
t135 = qJD(1) * t110;
t133 = qJD(1) * qJD(2);
t94 = t114 * qJDD(1) - t110 * t133;
t97 = qJD(2) * pkin(2) - pkin(8) * t135;
t121 = -t94 * pkin(2) + t97 * t135 + (-pkin(8) * t104 - pkin(7)) * t116 + t124;
t105 = sin(pkin(11));
t106 = cos(pkin(11));
t109 = sin(qJ(3));
t113 = cos(qJ(3));
t88 = (t109 * t114 + t110 * t113) * qJD(1);
t93 = t110 * qJDD(1) + t114 * t133;
t70 = -t88 * qJD(3) - t109 * t93 + t113 * t94;
t103 = qJD(2) + qJD(3);
t83 = t103 * pkin(3) - t88 * pkin(9);
t87 = (-t109 * t110 + t113 * t114) * qJD(1);
t86 = t87 ^ 2;
t118 = -t70 * pkin(3) - t86 * pkin(9) + t88 * t83 + t121;
t107 = sin(qJ(6));
t112 = cos(qJ(6));
t108 = sin(qJ(4));
t140 = cos(qJ(4));
t102 = qJDD(2) + qJDD(3);
t126 = -t115 * g(1) - t111 * g(2);
t90 = -t116 * pkin(1) + qJDD(1) * pkin(7) + t126;
t136 = t110 * t90;
t139 = pkin(2) * t116;
t60 = qJDD(2) * pkin(2) - t93 * pkin(8) - t136 + (pkin(8) * t133 + t110 * t139 - g(3)) * t114;
t131 = -t110 * g(3) + t114 * t90;
t63 = t94 * pkin(8) - qJD(2) * t97 - t104 * t139 + t131;
t129 = -t109 * t63 + t113 * t60;
t71 = t87 * qJD(3) + t109 * t94 + t113 * t93;
t32 = (t103 * t87 - t71) * pkin(9) + (t87 * t88 + t102) * pkin(3) + t129;
t137 = t109 * t60 + t113 * t63;
t36 = -t86 * pkin(3) + t70 * pkin(9) - t103 * t83 + t137;
t138 = t108 * t32 + t140 * t36;
t78 = t108 * t88 - t140 * t87;
t79 = t108 * t87 + t140 * t88;
t55 = pkin(4) * t78 - qJ(5) * t79;
t100 = qJD(4) + t103;
t98 = t100 ^ 2;
t99 = qJDD(4) + t102;
t23 = -t98 * pkin(4) + t99 * qJ(5) - t78 * t55 + t138;
t45 = t79 * qJD(4) + t108 * t71 - t140 * t70;
t46 = -t78 * qJD(4) + t108 * t70 + t140 * t71;
t26 = (t100 * t78 - t46) * qJ(5) + (t100 * t79 + t45) * pkin(4) + t118;
t67 = t105 * t100 + t106 * t79;
t128 = -0.2e1 * qJD(5) * t67 - t105 * t23 + t106 * t26;
t38 = t105 * t99 + t106 * t46;
t66 = t106 * t100 - t105 * t79;
t17 = (t66 * t78 - t38) * pkin(10) + (t66 * t67 + t45) * pkin(5) + t128;
t132 = 0.2e1 * qJD(5) * t66 + t105 * t26 + t106 * t23;
t37 = -t105 * t46 + t106 * t99;
t53 = pkin(5) * t78 - t67 * pkin(10);
t64 = t66 ^ 2;
t18 = -t64 * pkin(5) + t37 * pkin(10) - t53 * t78 + t132;
t48 = -t107 * t67 + t112 * t66;
t29 = t48 * qJD(6) + t107 * t37 + t112 * t38;
t49 = t107 * t66 + t112 * t67;
t33 = -mrSges(7,1) * t48 + mrSges(7,2) * t49;
t75 = qJD(6) + t78;
t41 = -mrSges(7,2) * t75 + t48 * mrSges(7,3);
t44 = qJDD(6) + t45;
t14 = m(7) * (-t107 * t18 + t112 * t17) - t29 * mrSges(7,3) + t44 * mrSges(7,1) - t49 * t33 + t75 * t41;
t28 = -t49 * qJD(6) - t107 * t38 + t112 * t37;
t42 = mrSges(7,1) * t75 - t49 * mrSges(7,3);
t15 = m(7) * (t107 * t17 + t112 * t18) + t28 * mrSges(7,3) - t44 * mrSges(7,2) + t48 * t33 - t75 * t42;
t50 = -mrSges(6,1) * t66 + mrSges(6,2) * t67;
t51 = -mrSges(6,2) * t78 + t66 * mrSges(6,3);
t12 = m(6) * t128 + t45 * mrSges(6,1) - t38 * mrSges(6,3) + t107 * t15 + t112 * t14 - t67 * t50 + t78 * t51;
t52 = mrSges(6,1) * t78 - t67 * mrSges(6,3);
t13 = m(6) * t132 - t45 * mrSges(6,2) + t37 * mrSges(6,3) - t107 * t14 + t112 * t15 + t66 * t50 - t78 * t52;
t73 = -t100 * mrSges(5,2) - t78 * mrSges(5,3);
t74 = t100 * mrSges(5,1) - t79 * mrSges(5,3);
t122 = m(5) * t118 + t45 * mrSges(5,1) + t46 * mrSges(5,2) + t105 * t13 + t106 * t12 + t78 * t73 + t79 * t74;
t81 = -t103 * mrSges(4,2) + t87 * mrSges(4,3);
t82 = t103 * mrSges(4,1) - t88 * mrSges(4,3);
t120 = -m(4) * t121 + t70 * mrSges(4,1) - t71 * mrSges(4,2) + t87 * t81 - t88 * t82 - t122;
t95 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t135;
t134 = qJD(1) * t114;
t96 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t134;
t142 = (t110 * t95 - t114 * t96) * qJD(1) + m(3) * (-t116 * pkin(7) + t124) - t94 * mrSges(3,1) + t93 * mrSges(3,2) - t120;
t127 = -t108 * t36 + t140 * t32;
t22 = -t99 * pkin(4) - t98 * qJ(5) + t79 * t55 + qJDD(5) - t127;
t123 = t28 * mrSges(7,1) + t48 * t41 - m(7) * (-t37 * pkin(5) - t64 * pkin(10) + t67 * t53 + t22) - t29 * mrSges(7,2) - t49 * t42;
t119 = m(6) * t22 - t37 * mrSges(6,1) + t38 * mrSges(6,2) - t66 * t51 + t67 * t52 - t123;
t56 = mrSges(5,1) * t78 + mrSges(5,2) * t79;
t16 = m(5) * t127 + t99 * mrSges(5,1) - t46 * mrSges(5,3) + t100 * t73 - t79 * t56 - t119;
t80 = -mrSges(4,1) * t87 + mrSges(4,2) * t88;
t9 = m(5) * t138 - t99 * mrSges(5,2) - t45 * mrSges(5,3) - t100 * t74 - t105 * t12 + t106 * t13 - t78 * t56;
t6 = m(4) * t129 + t102 * mrSges(4,1) - t71 * mrSges(4,3) + t103 * t81 + t108 * t9 + t140 * t16 - t88 * t80;
t7 = m(4) * t137 - t102 * mrSges(4,2) + t70 * mrSges(4,3) - t103 * t82 - t108 * t16 + t140 * t9 + t87 * t80;
t92 = (-mrSges(3,1) * t114 + mrSges(3,2) * t110) * qJD(1);
t4 = m(3) * (-t114 * g(3) - t136) - t93 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t92 * t135 + qJD(2) * t96 + t109 * t7 + t113 * t6;
t5 = m(3) * t131 - qJDD(2) * mrSges(3,2) + t94 * mrSges(3,3) - qJD(2) * t95 - t109 * t6 + t113 * t7 + t92 * t134;
t141 = t110 * t5 + t114 * t4;
t8 = m(2) * t130 + qJDD(1) * mrSges(2,1) - t116 * mrSges(2,2) - t142;
t1 = m(2) * t126 - t116 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t110 * t4 + t114 * t5;
t2 = [-m(1) * g(1) + t115 * t1 - t111 * t8, t1, t5, t7, t9, t13, t15; -m(1) * g(2) + t111 * t1 + t115 * t8, t8, t4, t6, t16, t12, t14; (-m(1) - m(2)) * g(3) + t141, -m(2) * g(3) + t141, t142, -t120, t122, t119, -t123;];
f_new  = t2;
