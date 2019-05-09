% Calculate vector of cutting forces with Newton-Euler
% S6RRRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 05:18
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRRP6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP6_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP6_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP6_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP6_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP6_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP6_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP6_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 05:08:09
% EndTime: 2019-05-08 05:08:23
% DurationCPUTime: 4.49s
% Computational Cost: add. (64164->201), mult. (130030->254), div. (0->0), fcn. (93122->10), ass. (0->100)
t100 = sin(qJ(3));
t101 = sin(qJ(2));
t104 = cos(qJ(3));
t105 = cos(qJ(2));
t103 = cos(qJ(4));
t125 = qJD(1) * qJD(2);
t120 = t105 * t125;
t108 = qJD(1) ^ 2;
t102 = sin(qJ(1));
t106 = cos(qJ(1));
t123 = t102 * g(1) - t106 * g(2);
t78 = -qJDD(1) * pkin(1) - t108 * pkin(7) - t123;
t87 = t101 * qJDD(1) + t120;
t95 = t101 * t125;
t88 = t105 * qJDD(1) - t95;
t55 = (-t87 - t120) * pkin(8) + (-t88 + t95) * pkin(2) + t78;
t107 = qJD(2) ^ 2;
t118 = -t106 * g(1) - t102 * g(2);
t79 = -t108 * pkin(1) + qJDD(1) * pkin(7) + t118;
t122 = -t101 * g(3) + t105 * t79;
t126 = t105 * qJD(1);
t86 = (-pkin(2) * t105 - pkin(8) * t101) * qJD(1);
t58 = -t107 * pkin(2) + qJDD(2) * pkin(8) + t86 * t126 + t122;
t119 = -t100 * t58 + t104 * t55;
t127 = qJD(1) * t101;
t83 = t104 * qJD(2) - t100 * t127;
t65 = t83 * qJD(3) + t100 * qJDD(2) + t104 * t87;
t84 = t100 * qJD(2) + t104 * t127;
t69 = -t83 * mrSges(4,1) + t84 * mrSges(4,2);
t134 = cos(qJ(5));
t82 = qJDD(3) - t88;
t94 = qJD(3) - t126;
t29 = (t83 * t94 - t65) * pkin(9) + (t83 * t84 + t82) * pkin(3) + t119;
t129 = t100 * t55 + t104 * t58;
t64 = -t84 * qJD(3) + t104 * qJDD(2) - t100 * t87;
t72 = t94 * pkin(3) - t84 * pkin(9);
t81 = t83 ^ 2;
t31 = -t81 * pkin(3) + t64 * pkin(9) - t94 * t72 + t129;
t99 = sin(qJ(4));
t121 = t103 * t29 - t99 * t31;
t67 = t103 * t83 - t99 * t84;
t42 = t67 * qJD(4) + t103 * t65 + t99 * t64;
t68 = t103 * t84 + t99 * t83;
t80 = qJDD(4) + t82;
t93 = qJD(4) + t94;
t18 = (t67 * t93 - t42) * pkin(10) + (t67 * t68 + t80) * pkin(4) + t121;
t131 = t103 * t31 + t99 * t29;
t41 = -t68 * qJD(4) + t103 * t64 - t99 * t65;
t61 = t93 * pkin(4) - t68 * pkin(10);
t66 = t67 ^ 2;
t20 = -t66 * pkin(4) + t41 * pkin(10) - t93 * t61 + t131;
t98 = sin(qJ(5));
t132 = t134 * t20 + t98 * t18;
t50 = -t134 * t67 + t98 * t68;
t51 = t134 * t68 + t98 * t67;
t34 = t50 * pkin(5) - t51 * qJ(6);
t90 = qJD(5) + t93;
t47 = -t90 * mrSges(7,1) + t51 * mrSges(7,2);
t77 = qJDD(5) + t80;
t89 = t90 ^ 2;
t124 = m(7) * (-t89 * pkin(5) + t77 * qJ(6) + 0.2e1 * qJD(6) * t90 - t50 * t34 + t132) + t90 * t47 + t77 * mrSges(7,3);
t35 = t50 * mrSges(7,1) - t51 * mrSges(7,3);
t130 = -t50 * mrSges(6,1) - t51 * mrSges(6,2) - t35;
t133 = -mrSges(6,3) - mrSges(7,2);
t25 = t51 * qJD(5) - t134 * t41 + t98 * t42;
t46 = t90 * mrSges(6,1) - t51 * mrSges(6,3);
t11 = m(6) * t132 - t77 * mrSges(6,2) + t130 * t50 + t133 * t25 - t90 * t46 + t124;
t116 = t134 * t18 - t98 * t20;
t135 = m(7) * (-t77 * pkin(5) - t89 * qJ(6) + t51 * t34 + qJDD(6) - t116);
t26 = -t50 * qJD(5) + t134 * t42 + t98 * t41;
t44 = -t50 * mrSges(7,2) + t90 * mrSges(7,3);
t45 = -t90 * mrSges(6,2) - t50 * mrSges(6,3);
t12 = m(6) * t116 - t135 + (t45 + t44) * t90 + (mrSges(6,1) + mrSges(7,1)) * t77 + t130 * t51 + t133 * t26;
t52 = -t67 * mrSges(5,1) + t68 * mrSges(5,2);
t59 = -t93 * mrSges(5,2) + t67 * mrSges(5,3);
t7 = m(5) * t121 + t80 * mrSges(5,1) - t42 * mrSges(5,3) + t98 * t11 + t134 * t12 - t68 * t52 + t93 * t59;
t70 = -t94 * mrSges(4,2) + t83 * mrSges(4,3);
t60 = t93 * mrSges(5,1) - t68 * mrSges(5,3);
t8 = m(5) * t131 - t80 * mrSges(5,2) + t41 * mrSges(5,3) + t134 * t11 - t98 * t12 + t67 * t52 - t93 * t60;
t5 = m(4) * t119 + t82 * mrSges(4,1) - t65 * mrSges(4,3) + t103 * t7 - t84 * t69 + t94 * t70 + t99 * t8;
t71 = t94 * mrSges(4,1) - t84 * mrSges(4,3);
t6 = m(4) * t129 - t82 * mrSges(4,2) + t64 * mrSges(4,3) + t103 * t8 + t83 * t69 - t99 * t7 - t94 * t71;
t91 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t127;
t92 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t126;
t137 = m(3) * t78 - t88 * mrSges(3,1) + t87 * mrSges(3,2) + t100 * t6 + t104 * t5 + (t101 * t91 - t105 * t92) * qJD(1);
t128 = -t105 * g(3) - t101 * t79;
t57 = -qJDD(2) * pkin(2) - t107 * pkin(8) + t86 * t127 - t128;
t113 = -t64 * pkin(3) - t81 * pkin(9) + t84 * t72 + t57;
t111 = -t41 * pkin(4) - t66 * pkin(10) + t68 * t61 + t113;
t115 = t26 * mrSges(7,3) + t51 * t47 - m(7) * (t111 + (t50 * t90 - t26) * qJ(6) + (t51 * t90 + t25) * pkin(5) - 0.2e1 * qJD(6) * t51) - t25 * mrSges(7,1) - t50 * t44;
t112 = m(6) * t111 + t25 * mrSges(6,1) + t26 * mrSges(6,2) + t50 * t45 + t51 * t46 - t115;
t110 = -m(5) * t113 + t41 * mrSges(5,1) - t42 * mrSges(5,2) + t67 * t59 - t68 * t60 - t112;
t109 = m(4) * t57 - t64 * mrSges(4,1) + t65 * mrSges(4,2) - t83 * t70 + t84 * t71 - t110;
t85 = (-mrSges(3,1) * t105 + mrSges(3,2) * t101) * qJD(1);
t10 = m(3) * t128 + qJDD(2) * mrSges(3,1) - t87 * mrSges(3,3) + qJD(2) * t92 - t85 * t127 - t109;
t4 = m(3) * t122 - qJDD(2) * mrSges(3,2) + t88 * mrSges(3,3) - qJD(2) * t91 - t100 * t5 + t104 * t6 + t85 * t126;
t136 = t105 * t10 + t101 * t4;
t2 = m(2) * t123 + qJDD(1) * mrSges(2,1) - t108 * mrSges(2,2) - t137;
t1 = m(2) * t118 - t108 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t101 * t10 + t105 * t4;
t3 = [-m(1) * g(1) + t106 * t1 - t102 * t2, t1, t4, t6, t8, t11, -t25 * mrSges(7,2) - t50 * t35 + t124; -m(1) * g(2) + t102 * t1 + t106 * t2, t2, t10, t5, t7, t12, -t115; (-m(1) - m(2)) * g(3) + t136, -m(2) * g(3) + t136, t137, t109, -t110, t112, -t77 * mrSges(7,1) + t26 * mrSges(7,2) + t51 * t35 - t90 * t44 + t135;];
f_new  = t3;
