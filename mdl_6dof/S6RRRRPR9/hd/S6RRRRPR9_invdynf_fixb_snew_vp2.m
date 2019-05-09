% Calculate vector of cutting forces with Newton-Euler
% S6RRRRPR9
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-07 22:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRPR9_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR9_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR9_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR9_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR9_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR9_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR9_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 22:11:38
% EndTime: 2019-05-07 22:12:04
% DurationCPUTime: 9.95s
% Computational Cost: add. (192558->214), mult. (412357->295), div. (0->0), fcn. (336796->14), ass. (0->115)
t149 = cos(qJ(4));
t113 = cos(pkin(6));
t148 = t113 * g(3);
t115 = sin(qJ(4));
t111 = sin(pkin(6));
t121 = cos(qJ(2));
t140 = qJD(1) * t121;
t135 = t111 * t140;
t103 = qJD(3) - t135;
t116 = sin(qJ(3));
t120 = cos(qJ(3));
t107 = t113 * qJD(1) + qJD(2);
t105 = t107 ^ 2;
t106 = t113 * qJDD(1) + qJDD(2);
t117 = sin(qJ(2));
t123 = qJD(1) ^ 2;
t118 = sin(qJ(1));
t122 = cos(qJ(1));
t134 = t118 * g(1) - t122 * g(2);
t95 = t123 * t111 * pkin(8) + qJDD(1) * pkin(1) + t134;
t144 = t113 * t95;
t130 = -t122 * g(1) - t118 * g(2);
t138 = qJDD(1) * t111;
t96 = -t123 * pkin(1) + pkin(8) * t138 + t130;
t145 = t117 * t144 + t121 * t96;
t141 = qJD(1) * t111;
t98 = (-pkin(2) * t121 - pkin(9) * t117) * t141;
t63 = -t105 * pkin(2) + t106 * pkin(9) + (-g(3) * t117 + t98 * t140) * t111 + t145;
t136 = t117 * t141;
t100 = -qJD(2) * t136 + t121 * t138;
t99 = (qJD(2) * t140 + qJDD(1) * t117) * t111;
t64 = -t100 * pkin(2) - t99 * pkin(9) - t148 + (-t95 + (pkin(2) * t117 - pkin(9) * t121) * t107 * qJD(1)) * t111;
t133 = -t116 * t63 + t120 * t64;
t87 = t120 * t107 - t116 * t136;
t75 = t87 * qJD(3) + t116 * t106 + t120 * t99;
t88 = t116 * t107 + t120 * t136;
t92 = qJDD(3) - t100;
t33 = (t103 * t87 - t75) * pkin(10) + (t87 * t88 + t92) * pkin(3) + t133;
t146 = t116 * t64 + t120 * t63;
t74 = -t88 * qJD(3) + t120 * t106 - t116 * t99;
t82 = t103 * pkin(3) - t88 * pkin(10);
t86 = t87 ^ 2;
t37 = -t86 * pkin(3) + t74 * pkin(10) - t103 * t82 + t146;
t147 = t115 * t33 + t149 * t37;
t143 = t111 * t117;
t142 = t111 * t121;
t110 = sin(pkin(12));
t112 = cos(pkin(12));
t129 = -g(3) * t142 - t117 * t96 + t121 * t144;
t62 = -t106 * pkin(2) - t105 * pkin(9) + t98 * t136 - t129;
t126 = -t74 * pkin(3) - t86 * pkin(10) + t88 * t82 + t62;
t114 = sin(qJ(6));
t119 = cos(qJ(6));
t102 = qJD(4) + t103;
t101 = t102 ^ 2;
t77 = t115 * t88 - t149 * t87;
t78 = t115 * t87 + t149 * t88;
t57 = t77 * pkin(4) - t78 * qJ(5);
t91 = qJDD(4) + t92;
t25 = -t101 * pkin(4) + t91 * qJ(5) - t77 * t57 + t147;
t47 = t78 * qJD(4) + t115 * t75 - t149 * t74;
t48 = -t77 * qJD(4) + t115 * t74 + t149 * t75;
t28 = (t102 * t77 - t48) * qJ(5) + (t102 * t78 + t47) * pkin(4) + t126;
t69 = t110 * t102 + t112 * t78;
t132 = -0.2e1 * qJD(5) * t69 - t110 * t25 + t112 * t28;
t44 = t110 * t91 + t112 * t48;
t68 = t112 * t102 - t110 * t78;
t19 = (t68 * t77 - t44) * pkin(11) + (t68 * t69 + t47) * pkin(5) + t132;
t137 = 0.2e1 * qJD(5) * t68 + t110 * t28 + t112 * t25;
t43 = -t110 * t48 + t112 * t91;
t55 = t77 * pkin(5) - t69 * pkin(11);
t67 = t68 ^ 2;
t20 = -t67 * pkin(5) + t43 * pkin(11) - t77 * t55 + t137;
t50 = -t114 * t69 + t119 * t68;
t31 = t50 * qJD(6) + t114 * t43 + t119 * t44;
t51 = t114 * t68 + t119 * t69;
t38 = -t50 * mrSges(7,1) + t51 * mrSges(7,2);
t76 = qJD(6) + t77;
t41 = -t76 * mrSges(7,2) + t50 * mrSges(7,3);
t46 = qJDD(6) + t47;
t17 = m(7) * (-t114 * t20 + t119 * t19) - t31 * mrSges(7,3) + t46 * mrSges(7,1) - t51 * t38 + t76 * t41;
t30 = -t51 * qJD(6) - t114 * t44 + t119 * t43;
t42 = t76 * mrSges(7,1) - t51 * mrSges(7,3);
t18 = m(7) * (t114 * t19 + t119 * t20) + t30 * mrSges(7,3) - t46 * mrSges(7,2) + t50 * t38 - t76 * t42;
t52 = -t68 * mrSges(6,1) + t69 * mrSges(6,2);
t53 = -t77 * mrSges(6,2) + t68 * mrSges(6,3);
t14 = m(6) * t132 + t47 * mrSges(6,1) - t44 * mrSges(6,3) + t114 * t18 + t119 * t17 - t69 * t52 + t77 * t53;
t54 = t77 * mrSges(6,1) - t69 * mrSges(6,3);
t15 = m(6) * t137 - t47 * mrSges(6,2) + t43 * mrSges(6,3) - t114 * t17 + t119 * t18 + t68 * t52 - t77 * t54;
t70 = -t102 * mrSges(5,2) - t77 * mrSges(5,3);
t71 = t102 * mrSges(5,1) - t78 * mrSges(5,3);
t127 = m(5) * t126 + t47 * mrSges(5,1) + t48 * mrSges(5,2) + t110 * t15 + t112 * t14 + t77 * t70 + t78 * t71;
t80 = -t103 * mrSges(4,2) + t87 * mrSges(4,3);
t81 = t103 * mrSges(4,1) - t88 * mrSges(4,3);
t124 = m(4) * t62 - t74 * mrSges(4,1) + t75 * mrSges(4,2) - t87 * t80 + t88 * t81 + t127;
t94 = -t107 * mrSges(3,2) + mrSges(3,3) * t135;
t97 = (-mrSges(3,1) * t121 + mrSges(3,2) * t117) * t141;
t10 = m(3) * t129 + t106 * mrSges(3,1) - t99 * mrSges(3,3) + t107 * t94 - t97 * t136 - t124;
t58 = t77 * mrSges(5,1) + t78 * mrSges(5,2);
t11 = m(5) * t147 - t91 * mrSges(5,2) - t47 * mrSges(5,3) - t102 * t71 - t110 * t14 + t112 * t15 - t77 * t58;
t131 = -t115 * t37 + t149 * t33;
t24 = -t91 * pkin(4) - t101 * qJ(5) + t78 * t57 + qJDD(5) - t131;
t128 = t30 * mrSges(7,1) + t50 * t41 - m(7) * (-t43 * pkin(5) - t67 * pkin(11) + t69 * t55 + t24) - t31 * mrSges(7,2) - t51 * t42;
t125 = m(6) * t24 - t43 * mrSges(6,1) + t44 * mrSges(6,2) - t68 * t53 + t69 * t54 - t128;
t16 = m(5) * t131 + t91 * mrSges(5,1) - t48 * mrSges(5,3) + t102 * t70 - t78 * t58 - t125;
t79 = -t87 * mrSges(4,1) + t88 * mrSges(4,2);
t7 = m(4) * t133 + t92 * mrSges(4,1) - t75 * mrSges(4,3) + t103 * t80 + t115 * t11 + t149 * t16 - t88 * t79;
t8 = m(4) * t146 - t92 * mrSges(4,2) + t74 * mrSges(4,3) - t103 * t81 + t149 * t11 - t115 * t16 + t87 * t79;
t93 = t107 * mrSges(3,1) - mrSges(3,3) * t136;
t4 = m(3) * (-g(3) * t143 + t145) + t100 * mrSges(3,3) - t106 * mrSges(3,2) + t97 * t135 - t107 * t93 + t120 * t8 - t116 * t7;
t6 = m(3) * (-t111 * t95 - t148) + t99 * mrSges(3,2) - t100 * mrSges(3,1) + t116 * t8 + t120 * t7 + (t117 * t93 - t121 * t94) * t141;
t139 = t10 * t142 + t113 * t6 + t4 * t143;
t2 = m(2) * t130 - t123 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t117 * t10 + t121 * t4;
t1 = m(2) * t134 + qJDD(1) * mrSges(2,1) - t123 * mrSges(2,2) - t111 * t6 + (t121 * t10 + t117 * t4) * t113;
t3 = [-m(1) * g(1) - t118 * t1 + t122 * t2, t2, t4, t8, t11, t15, t18; -m(1) * g(2) + t122 * t1 + t118 * t2, t1, t10, t7, t16, t14, t17; (-m(1) - m(2)) * g(3) + t139, -m(2) * g(3) + t139, t6, t124, t127, t125, -t128;];
f_new  = t3;
