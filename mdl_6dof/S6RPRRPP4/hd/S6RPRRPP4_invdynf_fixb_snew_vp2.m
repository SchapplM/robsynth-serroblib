% Calculate vector of cutting forces with Newton-Euler
% S6RPRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-05-05 21:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRPP4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:30:07
% EndTime: 2019-05-05 21:30:15
% DurationCPUTime: 3.24s
% Computational Cost: add. (39663->190), mult. (95060->240), div. (0->0), fcn. (70984->10), ass. (0->97)
t101 = qJD(1) ^ 2;
t93 = cos(pkin(9));
t90 = t93 ^ 2;
t92 = sin(pkin(9));
t128 = t92 ^ 2 + t90;
t137 = t128 * mrSges(3,3);
t95 = sin(qJ(3));
t98 = cos(qJ(3));
t111 = t92 * t95 - t93 * t98;
t82 = t111 * qJD(1);
t112 = t92 * t98 + t93 * t95;
t83 = t112 * qJD(1);
t124 = t83 * qJD(3);
t70 = -t111 * qJDD(1) - t124;
t136 = -2 * qJD(5);
t113 = -t93 * mrSges(3,1) + t92 * mrSges(3,2);
t110 = qJDD(1) * mrSges(3,3) + t101 * t113;
t123 = qJD(1) * qJD(2);
t119 = -t93 * g(3) - 0.2e1 * t92 * t123;
t100 = qJD(3) ^ 2;
t126 = pkin(7) * qJDD(1);
t133 = pkin(2) * t101;
t96 = sin(qJ(1));
t99 = cos(qJ(1));
t115 = -t99 * g(1) - t96 * g(2);
t84 = -t101 * pkin(1) + qJDD(1) * qJ(2) + t115;
t60 = (t93 * t133 - t126 - t84) * t92 + t119;
t116 = -t92 * g(3) + (0.2e1 * t123 + t84) * t93;
t61 = t93 * t126 - t90 * t133 + t116;
t117 = t98 * t60 - t95 * t61;
t68 = t82 * pkin(3) - t83 * pkin(8);
t26 = -qJDD(3) * pkin(3) - t100 * pkin(8) + t83 * t68 - t117;
t125 = t82 * qJD(3);
t71 = t112 * qJDD(1) - t125;
t94 = sin(qJ(4));
t97 = cos(qJ(4));
t75 = t94 * qJD(3) + t97 * t83;
t47 = -t75 * qJD(4) + t97 * qJDD(3) - t94 * t71;
t80 = qJD(4) + t82;
t56 = t80 * pkin(4) - t75 * qJ(5);
t74 = t97 * qJD(3) - t94 * t83;
t73 = t74 ^ 2;
t104 = -t47 * pkin(4) - t73 * qJ(5) + t75 * t56 + qJDD(5) + t26;
t127 = cos(pkin(10));
t48 = t74 * qJD(4) + t94 * qJDD(3) + t97 * t71;
t91 = sin(pkin(10));
t30 = -t127 * t47 + t91 * t48;
t31 = t127 * t48 + t91 * t47;
t50 = -t127 * t74 + t91 * t75;
t41 = -t50 * mrSges(7,2) + t80 * mrSges(7,3);
t51 = t127 * t75 + t91 * t74;
t44 = -t80 * mrSges(7,1) + t51 * mrSges(7,2);
t108 = t31 * mrSges(7,3) + t51 * t44 - m(7) * (-0.2e1 * qJD(6) * t51 + (t50 * t80 - t31) * qJ(6) + (t51 * t80 + t30) * pkin(5) + t104) - t30 * mrSges(7,1) - t50 * t41;
t42 = -t80 * mrSges(6,2) - t50 * mrSges(6,3);
t43 = t80 * mrSges(6,1) - t51 * mrSges(6,3);
t105 = m(6) * t104 + t30 * mrSges(6,1) + t31 * mrSges(6,2) + t50 * t42 + t51 * t43 - t108;
t55 = -t80 * mrSges(5,2) + t74 * mrSges(5,3);
t57 = t80 * mrSges(5,1) - t75 * mrSges(5,3);
t102 = m(5) * t26 - t47 * mrSges(5,1) + t48 * mrSges(5,2) - t74 * t55 + t75 * t57 + t105;
t65 = t82 * mrSges(4,1) + t83 * mrSges(4,2);
t76 = -qJD(3) * mrSges(4,2) - t82 * mrSges(4,3);
t14 = m(4) * t117 + qJDD(3) * mrSges(4,1) - t71 * mrSges(4,3) + qJD(3) * t76 - t83 * t65 - t102;
t129 = t95 * t60 + t98 * t61;
t27 = -t100 * pkin(3) + qJDD(3) * pkin(8) - t82 * t68 + t129;
t120 = t96 * g(1) - t99 * g(2);
t114 = qJDD(2) - t120;
t103 = (-pkin(2) * t93 - pkin(1)) * qJDD(1) + (-t128 * pkin(7) - qJ(2)) * t101 + t114;
t34 = (-t71 + t125) * pkin(8) + (-t70 + t124) * pkin(3) + t103;
t118 = -t94 * t27 + t97 * t34;
t67 = qJDD(4) - t70;
t20 = (t74 * t80 - t48) * qJ(5) + (t74 * t75 + t67) * pkin(4) + t118;
t131 = t97 * t27 + t94 * t34;
t22 = -t73 * pkin(4) + t47 * qJ(5) - t80 * t56 + t131;
t121 = t127 * t22 + t50 * t136 + t91 * t20;
t37 = t50 * pkin(5) - t51 * qJ(6);
t79 = t80 ^ 2;
t122 = m(7) * (-t79 * pkin(5) + t67 * qJ(6) + 0.2e1 * qJD(6) * t80 - t50 * t37 + t121) + t80 * t44 + t67 * mrSges(7,3);
t38 = t50 * mrSges(7,1) - t51 * mrSges(7,3);
t130 = -t50 * mrSges(6,1) - t51 * mrSges(6,2) - t38;
t132 = -mrSges(6,3) - mrSges(7,2);
t12 = m(6) * t121 - t67 * mrSges(6,2) + t130 * t50 + t132 * t30 - t80 * t43 + t122;
t109 = t127 * t20 - t91 * t22;
t134 = m(7) * (-t67 * pkin(5) - t79 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t37) * t51 - t109);
t13 = m(6) * t109 - t134 + (t42 + t41) * t80 + (mrSges(6,1) + mrSges(7,1)) * t67 + (m(6) * t136 + t130) * t51 + t132 * t31;
t52 = -t74 * mrSges(5,1) + t75 * mrSges(5,2);
t10 = m(5) * t118 + t67 * mrSges(5,1) - t48 * mrSges(5,3) + t91 * t12 + t127 * t13 - t75 * t52 + t80 * t55;
t11 = m(5) * t131 - t67 * mrSges(5,2) + t47 * mrSges(5,3) + t127 * t12 - t91 * t13 + t74 * t52 - t80 * t57;
t77 = qJD(3) * mrSges(4,1) - t83 * mrSges(4,3);
t7 = m(4) * t129 - qJDD(3) * mrSges(4,2) + t70 * mrSges(4,3) - qJD(3) * t77 - t94 * t10 + t97 * t11 - t82 * t65;
t4 = m(3) * t119 + t95 * t7 + t98 * t14 + (-m(3) * t84 - t110) * t92;
t5 = m(3) * t116 + t110 * t93 - t95 * t14 + t98 * t7;
t135 = t93 * t4 + t92 * t5;
t107 = -m(4) * t103 + t70 * mrSges(4,1) - t71 * mrSges(4,2) - t97 * t10 - t94 * t11 - t82 * t76 - t83 * t77;
t106 = m(3) * (-qJDD(1) * pkin(1) - t101 * qJ(2) + t114) - t107;
t6 = m(2) * t120 + (-mrSges(2,2) + t137) * t101 + (mrSges(2,1) - t113) * qJDD(1) - t106;
t1 = m(2) * t115 - t101 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t92 * t4 + t93 * t5;
t2 = [-m(1) * g(1) + t99 * t1 - t96 * t6, t1, t5, t7, t11, t12, -t30 * mrSges(7,2) - t50 * t38 + t122; -m(1) * g(2) + t96 * t1 + t99 * t6, t6, t4, t14, t10, t13, -t108; (-m(1) - m(2)) * g(3) + t135, -m(2) * g(3) + t135, t113 * qJDD(1) - t101 * t137 + t106, -t107, t102, t105, -t67 * mrSges(7,1) + t31 * mrSges(7,2) + t51 * t38 - t80 * t41 + t134;];
f_new  = t2;
