% Calculate vector of cutting forces with Newton-Euler
% S6RRRRPR7
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
% Datum: 2019-05-07 21:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRPR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR7_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR7_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR7_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR7_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR7_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 21:09:50
% EndTime: 2019-05-07 21:10:21
% DurationCPUTime: 10.94s
% Computational Cost: add. (209296->216), mult. (458903->295), div. (0->0), fcn. (375573->14), ass. (0->116)
t147 = 2 * qJD(5);
t109 = cos(pkin(6));
t146 = t109 * g(3);
t111 = sin(qJ(4));
t116 = cos(qJ(4));
t112 = sin(qJ(3));
t117 = cos(qJ(3));
t103 = t109 * qJD(1) + qJD(2);
t101 = t103 ^ 2;
t102 = t109 * qJDD(1) + qJDD(2);
t107 = sin(pkin(6));
t113 = sin(qJ(2));
t118 = cos(qJ(2));
t138 = qJD(1) * t118;
t120 = qJD(1) ^ 2;
t114 = sin(qJ(1));
t119 = cos(qJ(1));
t132 = t114 * g(1) - t119 * g(2);
t91 = t120 * t107 * pkin(8) + qJDD(1) * pkin(1) + t132;
t142 = t109 * t91;
t129 = -t119 * g(1) - t114 * g(2);
t137 = qJDD(1) * t107;
t92 = -t120 * pkin(1) + pkin(8) * t137 + t129;
t143 = t113 * t142 + t118 * t92;
t139 = qJD(1) * t107;
t94 = (-pkin(2) * t118 - pkin(9) * t113) * t139;
t65 = -t101 * pkin(2) + t102 * pkin(9) + (-g(3) * t113 + t94 * t138) * t107 + t143;
t95 = (qJD(2) * t138 + qJDD(1) * t113) * t107;
t134 = t113 * t139;
t96 = -qJD(2) * t134 + t118 * t137;
t66 = -t96 * pkin(2) - t95 * pkin(9) - t146 + (-t91 + (pkin(2) * t113 - pkin(9) * t118) * t103 * qJD(1)) * t107;
t130 = -t112 * t65 + t117 * t66;
t83 = t117 * t103 - t112 * t134;
t72 = t83 * qJD(3) + t112 * t102 + t117 * t95;
t84 = t112 * t103 + t117 * t134;
t88 = qJDD(3) - t96;
t133 = t107 * t138;
t99 = qJD(3) - t133;
t32 = (t83 * t99 - t72) * pkin(10) + (t83 * t84 + t88) * pkin(3) + t130;
t144 = t112 * t66 + t117 * t65;
t71 = -t84 * qJD(3) + t117 * t102 - t112 * t95;
t79 = t99 * pkin(3) - t84 * pkin(10);
t82 = t83 ^ 2;
t38 = -t82 * pkin(3) + t71 * pkin(10) - t99 * t79 + t144;
t145 = t111 * t32 + t116 * t38;
t141 = t107 * t113;
t140 = t107 * t118;
t127 = -g(3) * t140 - t113 * t92 + t118 * t142;
t64 = -t102 * pkin(2) - t101 * pkin(9) + t94 * t134 - t127;
t124 = -t71 * pkin(3) - t82 * pkin(10) + t84 * t79 + t64;
t110 = sin(qJ(6));
t115 = cos(qJ(6));
t75 = t111 * t83 + t116 * t84;
t48 = -t75 * qJD(4) - t111 * t72 + t116 * t71;
t98 = qJD(4) + t99;
t68 = t98 * pkin(4) - t75 * qJ(5);
t74 = -t111 * t84 + t116 * t83;
t73 = t74 ^ 2;
t122 = -t48 * pkin(4) - t73 * qJ(5) + t75 * t68 + qJDD(5) + t124;
t106 = sin(pkin(12));
t108 = cos(pkin(12));
t131 = -t111 * t38 + t116 * t32;
t49 = t74 * qJD(4) + t111 * t71 + t116 * t72;
t87 = qJDD(4) + t88;
t23 = (t74 * t98 - t49) * qJ(5) + (t74 * t75 + t87) * pkin(4) + t131;
t25 = -t73 * pkin(4) + t48 * qJ(5) - t98 * t68 + t145;
t58 = -t106 * t75 + t108 * t74;
t135 = t106 * t23 + t108 * t25 + t58 * t147;
t59 = t106 * t74 + t108 * t75;
t46 = -t58 * pkin(5) - t59 * pkin(11);
t97 = t98 ^ 2;
t20 = -t97 * pkin(5) + t87 * pkin(11) + t58 * t46 + t135;
t35 = -t106 * t49 + t108 * t48;
t36 = t106 * t48 + t108 * t49;
t21 = t122 + (-t58 * t98 - t36) * pkin(11) + (t59 * t98 - t35) * pkin(5);
t50 = -t110 * t59 + t115 * t98;
t29 = t50 * qJD(6) + t110 * t87 + t115 * t36;
t34 = qJDD(6) - t35;
t51 = t110 * t98 + t115 * t59;
t39 = -t50 * mrSges(7,1) + t51 * mrSges(7,2);
t57 = qJD(6) - t58;
t40 = -t57 * mrSges(7,2) + t50 * mrSges(7,3);
t17 = m(7) * (-t110 * t20 + t115 * t21) - t29 * mrSges(7,3) + t34 * mrSges(7,1) - t51 * t39 + t57 * t40;
t28 = -t51 * qJD(6) - t110 * t36 + t115 * t87;
t41 = t57 * mrSges(7,1) - t51 * mrSges(7,3);
t18 = m(7) * (t110 * t21 + t115 * t20) + t28 * mrSges(7,3) - t34 * mrSges(7,2) + t50 * t39 - t57 * t41;
t52 = -t98 * mrSges(6,2) + t58 * mrSges(6,3);
t53 = t98 * mrSges(6,1) - t59 * mrSges(6,3);
t126 = -m(6) * t122 + t35 * mrSges(6,1) - t36 * mrSges(6,2) - t110 * t18 - t115 * t17 + t58 * t52 - t59 * t53;
t67 = -t98 * mrSges(5,2) + t74 * mrSges(5,3);
t69 = t98 * mrSges(5,1) - t75 * mrSges(5,3);
t123 = -m(5) * t124 + t48 * mrSges(5,1) - t49 * mrSges(5,2) + t74 * t67 - t75 * t69 + t126;
t77 = -t99 * mrSges(4,2) + t83 * mrSges(4,3);
t78 = t99 * mrSges(4,1) - t84 * mrSges(4,3);
t121 = m(4) * t64 - t71 * mrSges(4,1) + t72 * mrSges(4,2) - t83 * t77 + t84 * t78 - t123;
t90 = -t103 * mrSges(3,2) + mrSges(3,3) * t133;
t93 = (-mrSges(3,1) * t118 + mrSges(3,2) * t113) * t139;
t12 = m(3) * t127 + t102 * mrSges(3,1) - t95 * mrSges(3,3) + t103 * t90 - t93 * t134 - t121;
t45 = -t58 * mrSges(6,1) + t59 * mrSges(6,2);
t13 = m(6) * t135 - t87 * mrSges(6,2) + t35 * mrSges(6,3) - t110 * t17 + t115 * t18 + t58 * t45 - t98 * t53;
t128 = t106 * t25 - t108 * t23;
t125 = m(7) * (-t87 * pkin(5) - t97 * pkin(11) + (t147 + t46) * t59 + t128) - t28 * mrSges(7,1) + t29 * mrSges(7,2) - t50 * t40 + t51 * t41;
t14 = m(6) * (-0.2e1 * qJD(5) * t59 - t128) - t36 * mrSges(6,3) + t87 * mrSges(6,1) - t59 * t45 + t98 * t52 - t125;
t60 = -t74 * mrSges(5,1) + t75 * mrSges(5,2);
t10 = m(5) * t145 - t87 * mrSges(5,2) + t48 * mrSges(5,3) - t106 * t14 + t108 * t13 + t74 * t60 - t98 * t69;
t76 = -t83 * mrSges(4,1) + t84 * mrSges(4,2);
t9 = m(5) * t131 + t87 * mrSges(5,1) - t49 * mrSges(5,3) + t106 * t13 + t108 * t14 - t75 * t60 + t98 * t67;
t7 = m(4) * t130 + t88 * mrSges(4,1) - t72 * mrSges(4,3) + t111 * t10 + t116 * t9 - t84 * t76 + t99 * t77;
t8 = m(4) * t144 - t88 * mrSges(4,2) + t71 * mrSges(4,3) + t116 * t10 - t111 * t9 + t83 * t76 - t99 * t78;
t89 = t103 * mrSges(3,1) - mrSges(3,3) * t134;
t4 = m(3) * (-g(3) * t141 + t143) + t96 * mrSges(3,3) - t102 * mrSges(3,2) + t93 * t133 - t103 * t89 + t117 * t8 - t112 * t7;
t6 = m(3) * (-t107 * t91 - t146) + t95 * mrSges(3,2) - t96 * mrSges(3,1) + t112 * t8 + t117 * t7 + (t113 * t89 - t118 * t90) * t139;
t136 = t109 * t6 + t12 * t140 + t4 * t141;
t2 = m(2) * t129 - t120 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t113 * t12 + t118 * t4;
t1 = m(2) * t132 + qJDD(1) * mrSges(2,1) - t120 * mrSges(2,2) - t107 * t6 + (t113 * t4 + t118 * t12) * t109;
t3 = [-m(1) * g(1) - t114 * t1 + t119 * t2, t2, t4, t8, t10, t13, t18; -m(1) * g(2) + t119 * t1 + t114 * t2, t1, t12, t7, t9, t14, t17; (-m(1) - m(2)) * g(3) + t136, -m(2) * g(3) + t136, t6, t121, -t123, -t126, t125;];
f_new  = t3;
