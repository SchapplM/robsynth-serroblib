% Calculate vector of cutting forces with Newton-Euler
% S6RRRRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 10:25
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRRR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 10:03:51
% EndTime: 2019-05-08 10:04:24
% DurationCPUTime: 11.43s
% Computational Cost: add. (218694->216), mult. (473985->293), div. (0->0), fcn. (390445->14), ass. (0->117)
t108 = cos(pkin(6));
t147 = t108 * g(3);
t110 = sin(qJ(5));
t116 = cos(qJ(5));
t111 = sin(qJ(4));
t117 = cos(qJ(4));
t107 = sin(pkin(6));
t119 = cos(qJ(2));
t138 = qJD(1) * t119;
t134 = t107 * t138;
t100 = qJD(3) - t134;
t112 = sin(qJ(3));
t118 = cos(qJ(3));
t104 = t108 * qJD(1) + qJD(2);
t102 = t104 ^ 2;
t103 = t108 * qJDD(1) + qJDD(2);
t113 = sin(qJ(2));
t121 = qJD(1) ^ 2;
t114 = sin(qJ(1));
t120 = cos(qJ(1));
t133 = t114 * g(1) - t120 * g(2);
t91 = t121 * t107 * pkin(8) + qJDD(1) * pkin(1) + t133;
t142 = t108 * t91;
t130 = -t120 * g(1) - t114 * g(2);
t137 = qJDD(1) * t107;
t92 = -t121 * pkin(1) + pkin(8) * t137 + t130;
t143 = t113 * t142 + t119 * t92;
t139 = qJD(1) * t107;
t94 = (-pkin(2) * t119 - pkin(9) * t113) * t139;
t64 = -t102 * pkin(2) + t103 * pkin(9) + (-g(3) * t113 + t94 * t138) * t107 + t143;
t95 = (qJD(2) * t138 + qJDD(1) * t113) * t107;
t135 = t113 * t139;
t96 = -qJD(2) * t135 + t119 * t137;
t65 = -t96 * pkin(2) - t95 * pkin(9) - t147 + (-t91 + (pkin(2) * t113 - pkin(9) * t119) * t104 * qJD(1)) * t107;
t131 = -t112 * t64 + t118 * t65;
t83 = t118 * t104 - t112 * t135;
t71 = t83 * qJD(3) + t112 * t103 + t118 * t95;
t84 = t112 * t104 + t118 * t135;
t88 = qJDD(3) - t96;
t36 = (t100 * t83 - t71) * pkin(10) + (t83 * t84 + t88) * pkin(3) + t131;
t144 = t112 * t65 + t118 * t64;
t70 = -t84 * qJD(3) + t118 * t103 - t112 * t95;
t78 = t100 * pkin(3) - t84 * pkin(10);
t82 = t83 ^ 2;
t38 = -t82 * pkin(3) + t70 * pkin(10) - t100 * t78 + t144;
t132 = -t111 * t38 + t117 * t36;
t73 = -t111 * t84 + t117 * t83;
t49 = t73 * qJD(4) + t111 * t70 + t117 * t71;
t74 = t111 * t83 + t117 * t84;
t87 = qJDD(4) + t88;
t99 = qJD(4) + t100;
t23 = (t73 * t99 - t49) * pkin(11) + (t73 * t74 + t87) * pkin(4) + t132;
t145 = t111 * t36 + t117 * t38;
t48 = -t74 * qJD(4) - t111 * t71 + t117 * t70;
t68 = t99 * pkin(4) - t74 * pkin(11);
t72 = t73 ^ 2;
t25 = -t72 * pkin(4) + t48 * pkin(11) - t99 * t68 + t145;
t146 = t110 * t23 + t116 * t25;
t141 = t107 * t113;
t140 = t107 * t119;
t128 = -g(3) * t140 - t113 * t92 + t119 * t142;
t63 = -t103 * pkin(2) - t102 * pkin(9) + t94 * t135 - t128;
t125 = -t70 * pkin(3) - t82 * pkin(10) + t84 * t78 + t63;
t109 = sin(qJ(6));
t115 = cos(qJ(6));
t123 = -t48 * pkin(4) - t72 * pkin(11) + t74 * t68 + t125;
t57 = -t110 * t74 + t116 * t73;
t58 = t110 * t73 + t116 * t74;
t46 = -t57 * pkin(5) - t58 * pkin(12);
t81 = qJDD(5) + t87;
t98 = qJD(5) + t99;
t97 = t98 ^ 2;
t20 = -t97 * pkin(5) + t81 * pkin(12) + t57 * t46 + t146;
t32 = -t58 * qJD(5) - t110 * t49 + t116 * t48;
t33 = t57 * qJD(5) + t110 * t48 + t116 * t49;
t21 = t123 + (-t57 * t98 - t33) * pkin(12) + (t58 * t98 - t32) * pkin(5);
t50 = -t109 * t58 + t115 * t98;
t27 = t50 * qJD(6) + t109 * t81 + t115 * t33;
t31 = qJDD(6) - t32;
t51 = t109 * t98 + t115 * t58;
t39 = -t50 * mrSges(7,1) + t51 * mrSges(7,2);
t56 = qJD(6) - t57;
t40 = -t56 * mrSges(7,2) + t50 * mrSges(7,3);
t17 = m(7) * (-t109 * t20 + t115 * t21) - t27 * mrSges(7,3) + t31 * mrSges(7,1) - t51 * t39 + t56 * t40;
t26 = -t51 * qJD(6) - t109 * t33 + t115 * t81;
t41 = t56 * mrSges(7,1) - t51 * mrSges(7,3);
t18 = m(7) * (t109 * t21 + t115 * t20) + t26 * mrSges(7,3) - t31 * mrSges(7,2) + t50 * t39 - t56 * t41;
t52 = -t98 * mrSges(6,2) + t57 * mrSges(6,3);
t53 = t98 * mrSges(6,1) - t58 * mrSges(6,3);
t127 = -m(6) * t123 + t32 * mrSges(6,1) - t33 * mrSges(6,2) - t109 * t18 - t115 * t17 + t57 * t52 - t58 * t53;
t66 = -t99 * mrSges(5,2) + t73 * mrSges(5,3);
t67 = t99 * mrSges(5,1) - t74 * mrSges(5,3);
t124 = -m(5) * t125 + t48 * mrSges(5,1) - t49 * mrSges(5,2) + t73 * t66 - t74 * t67 + t127;
t76 = -t100 * mrSges(4,2) + t83 * mrSges(4,3);
t77 = t100 * mrSges(4,1) - t84 * mrSges(4,3);
t122 = m(4) * t63 - t70 * mrSges(4,1) + t71 * mrSges(4,2) - t83 * t76 + t84 * t77 - t124;
t90 = -t104 * mrSges(3,2) + mrSges(3,3) * t134;
t93 = (-mrSges(3,1) * t119 + mrSges(3,2) * t113) * t139;
t12 = m(3) * t128 + t103 * mrSges(3,1) - t95 * mrSges(3,3) + t104 * t90 - t93 * t135 - t122;
t45 = -t57 * mrSges(6,1) + t58 * mrSges(6,2);
t13 = m(6) * t146 - t81 * mrSges(6,2) + t32 * mrSges(6,3) - t109 * t17 + t115 * t18 + t57 * t45 - t98 * t53;
t129 = -t110 * t25 + t116 * t23;
t126 = m(7) * (-t81 * pkin(5) - t97 * pkin(12) + t58 * t46 - t129) - t26 * mrSges(7,1) + t27 * mrSges(7,2) - t50 * t40 + t51 * t41;
t14 = m(6) * t129 + t81 * mrSges(6,1) - t33 * mrSges(6,3) - t58 * t45 + t98 * t52 - t126;
t59 = -t73 * mrSges(5,1) + t74 * mrSges(5,2);
t10 = m(5) * t145 - t87 * mrSges(5,2) + t48 * mrSges(5,3) - t110 * t14 + t116 * t13 + t73 * t59 - t99 * t67;
t75 = -t83 * mrSges(4,1) + t84 * mrSges(4,2);
t9 = m(5) * t132 + t87 * mrSges(5,1) - t49 * mrSges(5,3) + t110 * t13 + t116 * t14 - t74 * t59 + t99 * t66;
t7 = m(4) * t131 + t88 * mrSges(4,1) - t71 * mrSges(4,3) + t111 * t10 + t100 * t76 + t117 * t9 - t84 * t75;
t8 = m(4) * t144 - t88 * mrSges(4,2) + t70 * mrSges(4,3) + t117 * t10 - t100 * t77 - t111 * t9 + t83 * t75;
t89 = t104 * mrSges(3,1) - mrSges(3,3) * t135;
t4 = m(3) * (-g(3) * t141 + t143) + t96 * mrSges(3,3) - t103 * mrSges(3,2) + t93 * t134 - t104 * t89 + t118 * t8 - t112 * t7;
t6 = m(3) * (-t107 * t91 - t147) + t95 * mrSges(3,2) - t96 * mrSges(3,1) + t112 * t8 + t118 * t7 + (t113 * t89 - t119 * t90) * t139;
t136 = t108 * t6 + t12 * t140 + t4 * t141;
t2 = m(2) * t130 - t121 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t113 * t12 + t119 * t4;
t1 = m(2) * t133 + qJDD(1) * mrSges(2,1) - t121 * mrSges(2,2) - t107 * t6 + (t113 * t4 + t119 * t12) * t108;
t3 = [-m(1) * g(1) - t114 * t1 + t120 * t2, t2, t4, t8, t10, t13, t18; -m(1) * g(2) + t120 * t1 + t114 * t2, t1, t12, t7, t9, t14, t17; (-m(1) - m(2)) * g(3) + t136, -m(2) * g(3) + t136, t6, t122, -t124, -t127, t126;];
f_new  = t3;
