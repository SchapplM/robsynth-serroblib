% Calculate vector of cutting forces with Newton-Euler
% S6RRRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-05-07 18:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRPP4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:14:45
% EndTime: 2019-05-07 18:14:56
% DurationCPUTime: 4.47s
% Computational Cost: add. (62851->201), mult. (129118->256), div. (0->0), fcn. (92078->10), ass. (0->99)
t100 = sin(qJ(2));
t103 = cos(qJ(3));
t104 = cos(qJ(2));
t102 = cos(qJ(4));
t125 = qJD(1) * qJD(2);
t118 = t104 * t125;
t107 = qJD(1) ^ 2;
t101 = sin(qJ(1));
t105 = cos(qJ(1));
t122 = t101 * g(1) - t105 * g(2);
t78 = -qJDD(1) * pkin(1) - t107 * pkin(7) - t122;
t87 = t100 * qJDD(1) + t118;
t94 = t100 * t125;
t88 = t104 * qJDD(1) - t94;
t56 = (-t87 - t118) * pkin(8) + (-t88 + t94) * pkin(2) + t78;
t106 = qJD(2) ^ 2;
t117 = -t105 * g(1) - t101 * g(2);
t79 = -t107 * pkin(1) + qJDD(1) * pkin(7) + t117;
t121 = -t100 * g(3) + t104 * t79;
t126 = t104 * qJD(1);
t86 = (-pkin(2) * t104 - pkin(8) * t100) * qJD(1);
t59 = -t106 * pkin(2) + qJDD(2) * pkin(8) + t86 * t126 + t121;
t99 = sin(qJ(3));
t119 = t103 * t56 - t99 * t59;
t127 = qJD(1) * t100;
t83 = t103 * qJD(2) - t99 * t127;
t66 = t83 * qJD(3) + t99 * qJDD(2) + t103 * t87;
t128 = cos(pkin(10));
t82 = qJDD(3) - t88;
t84 = t99 * qJD(2) + t103 * t127;
t93 = qJD(3) - t126;
t29 = (t83 * t93 - t66) * pkin(9) + (t83 * t84 + t82) * pkin(3) + t119;
t130 = t103 * t59 + t99 * t56;
t65 = -t84 * qJD(3) + t103 * qJDD(2) - t99 * t87;
t73 = t93 * pkin(3) - t84 * pkin(9);
t81 = t83 ^ 2;
t31 = -t81 * pkin(3) + t65 * pkin(9) - t93 * t73 + t130;
t98 = sin(qJ(4));
t120 = t102 * t29 - t98 * t31;
t68 = t102 * t83 - t98 * t84;
t42 = t68 * qJD(4) + t102 * t66 + t98 * t65;
t69 = t102 * t84 + t98 * t83;
t80 = qJDD(4) + t82;
t92 = qJD(4) + t93;
t18 = (t68 * t92 - t42) * qJ(5) + (t68 * t69 + t80) * pkin(4) + t120;
t132 = t102 * t31 + t98 * t29;
t41 = -t69 * qJD(4) + t102 * t65 - t98 * t66;
t61 = t92 * pkin(4) - t69 * qJ(5);
t67 = t68 ^ 2;
t20 = -t67 * pkin(4) + t41 * qJ(5) - t92 * t61 + t132;
t97 = sin(pkin(10));
t115 = t128 * t18 - t97 * t20;
t51 = -t128 * t68 + t97 * t69;
t52 = t128 * t69 + t97 * t68;
t35 = t51 * mrSges(7,1) - t52 * mrSges(7,3);
t131 = -t51 * mrSges(6,1) - t52 * mrSges(6,2) - t35;
t133 = -mrSges(6,3) - mrSges(7,2);
t34 = t51 * pkin(5) - t52 * qJ(6);
t91 = t92 ^ 2;
t135 = m(7) * (-t80 * pkin(5) - t91 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t34) * t52 - t115);
t136 = -2 * qJD(5);
t26 = t128 * t42 + t97 * t41;
t44 = -t51 * mrSges(7,2) + t92 * mrSges(7,3);
t45 = -t92 * mrSges(6,2) - t51 * mrSges(6,3);
t10 = m(6) * t115 - t135 + (t45 + t44) * t92 + (mrSges(6,1) + mrSges(7,1)) * t80 + (m(6) * t136 + t131) * t52 + t133 * t26;
t53 = -t68 * mrSges(5,1) + t69 * mrSges(5,2);
t60 = -t92 * mrSges(5,2) + t68 * mrSges(5,3);
t123 = t128 * t20 + t51 * t136 + t97 * t18;
t47 = -t92 * mrSges(7,1) + t52 * mrSges(7,2);
t124 = m(7) * (-t91 * pkin(5) + t80 * qJ(6) + 0.2e1 * qJD(6) * t92 - t51 * t34 + t123) + t92 * t47 + t80 * mrSges(7,3);
t25 = -t128 * t41 + t97 * t42;
t46 = t92 * mrSges(6,1) - t52 * mrSges(6,3);
t9 = m(6) * t123 - t80 * mrSges(6,2) + t131 * t51 + t133 * t25 - t92 * t46 + t124;
t7 = m(5) * t120 + t80 * mrSges(5,1) - t42 * mrSges(5,3) + t128 * t10 - t69 * t53 + t92 * t60 + t97 * t9;
t70 = -t83 * mrSges(4,1) + t84 * mrSges(4,2);
t71 = -t93 * mrSges(4,2) + t83 * mrSges(4,3);
t62 = t92 * mrSges(5,1) - t69 * mrSges(5,3);
t8 = m(5) * t132 - t80 * mrSges(5,2) + t41 * mrSges(5,3) - t97 * t10 + t128 * t9 + t68 * t53 - t92 * t62;
t5 = m(4) * t119 + t82 * mrSges(4,1) - t66 * mrSges(4,3) + t102 * t7 - t84 * t70 + t93 * t71 + t98 * t8;
t72 = t93 * mrSges(4,1) - t84 * mrSges(4,3);
t6 = m(4) * t130 - t82 * mrSges(4,2) + t65 * mrSges(4,3) + t102 * t8 - t98 * t7 + t83 * t70 - t93 * t72;
t89 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t127;
t90 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t126;
t137 = m(3) * t78 - t88 * mrSges(3,1) + t87 * mrSges(3,2) + t103 * t5 + t99 * t6 + (t100 * t89 - t104 * t90) * qJD(1);
t129 = -t104 * g(3) - t100 * t79;
t58 = -qJDD(2) * pkin(2) - t106 * pkin(8) + t86 * t127 - t129;
t112 = -t65 * pkin(3) - t81 * pkin(9) + t84 * t73 + t58;
t110 = -t41 * pkin(4) - t67 * qJ(5) + t69 * t61 + qJDD(5) + t112;
t114 = t26 * mrSges(7,3) + t52 * t47 - m(7) * (t110 + (t52 * t92 + t25) * pkin(5) + (t51 * t92 - t26) * qJ(6) - 0.2e1 * qJD(6) * t52) - t25 * mrSges(7,1) - t51 * t44;
t111 = m(6) * t110 + t25 * mrSges(6,1) + t26 * mrSges(6,2) + t51 * t45 + t52 * t46 - t114;
t109 = -m(5) * t112 + t41 * mrSges(5,1) - t42 * mrSges(5,2) + t68 * t60 - t69 * t62 - t111;
t108 = m(4) * t58 - t65 * mrSges(4,1) + t66 * mrSges(4,2) - t83 * t71 + t84 * t72 - t109;
t85 = (-mrSges(3,1) * t104 + mrSges(3,2) * t100) * qJD(1);
t12 = m(3) * t129 + qJDD(2) * mrSges(3,1) - t87 * mrSges(3,3) + qJD(2) * t90 - t85 * t127 - t108;
t4 = m(3) * t121 - qJDD(2) * mrSges(3,2) + t88 * mrSges(3,3) - qJD(2) * t89 + t103 * t6 + t85 * t126 - t99 * t5;
t134 = t100 * t4 + t104 * t12;
t2 = m(2) * t122 + qJDD(1) * mrSges(2,1) - t107 * mrSges(2,2) - t137;
t1 = m(2) * t117 - t107 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t100 * t12 + t104 * t4;
t3 = [-m(1) * g(1) + t105 * t1 - t101 * t2, t1, t4, t6, t8, t9, -t25 * mrSges(7,2) - t51 * t35 + t124; -m(1) * g(2) + t101 * t1 + t105 * t2, t2, t12, t5, t7, t10, -t114; (-m(1) - m(2)) * g(3) + t134, -m(2) * g(3) + t134, t137, t108, -t109, t111, -t80 * mrSges(7,1) + t26 * mrSges(7,2) + t52 * t35 - t92 * t44 + t135;];
f_new  = t3;
