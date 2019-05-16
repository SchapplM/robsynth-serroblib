% Calculate vector of cutting forces with Newton-Euler
% S6RRRRPP1
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
% Datum: 2019-05-07 18:00
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRPP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 17:53:55
% EndTime: 2019-05-07 17:54:05
% DurationCPUTime: 3.79s
% Computational Cost: add. (52883->203), mult. (107139->262), div. (0->0), fcn. (76770->10), ass. (0->98)
t100 = sin(qJ(2));
t103 = cos(qJ(3));
t104 = cos(qJ(2));
t99 = sin(qJ(3));
t81 = (t100 * t99 - t103 * t104) * qJD(1);
t106 = qJD(1) ^ 2;
t128 = cos(pkin(10));
t137 = -2 * qJD(5);
t102 = cos(qJ(4));
t101 = sin(qJ(1));
t105 = cos(qJ(1));
t122 = t101 * g(1) - t105 * g(2);
t115 = -qJDD(1) * pkin(1) - t122;
t127 = qJD(1) * t100;
t125 = qJD(1) * qJD(2);
t88 = t104 * qJDD(1) - t100 * t125;
t91 = qJD(2) * pkin(2) - pkin(8) * t127;
t96 = t104 ^ 2;
t111 = -t88 * pkin(2) + t91 * t127 + (-pkin(8) * t96 - pkin(7)) * t106 + t115;
t82 = (t100 * t103 + t104 * t99) * qJD(1);
t87 = t100 * qJDD(1) + t104 * t125;
t63 = -t82 * qJD(3) + t103 * t88 - t99 * t87;
t64 = -t81 * qJD(3) + t103 * t87 + t99 * t88;
t95 = qJD(2) + qJD(3);
t27 = (t81 * t95 - t64) * pkin(9) + (t82 * t95 - t63) * pkin(3) + t111;
t118 = -t105 * g(1) - t101 * g(2);
t84 = -t106 * pkin(1) + qJDD(1) * pkin(7) + t118;
t129 = t100 * t84;
t134 = pkin(2) * t106;
t56 = qJDD(2) * pkin(2) - t87 * pkin(8) - t129 + (pkin(8) * t125 + t100 * t134 - g(3)) * t104;
t121 = -t100 * g(3) + t104 * t84;
t57 = t88 * pkin(8) - qJD(2) * t91 - t96 * t134 + t121;
t130 = t103 * t57 + t99 * t56;
t71 = t81 * pkin(3) - t82 * pkin(9);
t93 = t95 ^ 2;
t94 = qJDD(2) + qJDD(3);
t34 = -t93 * pkin(3) + t94 * pkin(9) - t81 * t71 + t130;
t98 = sin(qJ(4));
t120 = t102 * t27 - t98 * t34;
t74 = t102 * t95 - t98 * t82;
t43 = t74 * qJD(4) + t102 * t64 + t98 * t94;
t61 = qJDD(4) - t63;
t75 = t102 * t82 + t98 * t95;
t80 = qJD(4) + t81;
t20 = (t74 * t80 - t43) * qJ(5) + (t74 * t75 + t61) * pkin(4) + t120;
t132 = t102 * t34 + t98 * t27;
t42 = -t75 * qJD(4) + t102 * t94 - t98 * t64;
t67 = t80 * pkin(4) - t75 * qJ(5);
t73 = t74 ^ 2;
t22 = -t73 * pkin(4) + t42 * qJ(5) - t80 * t67 + t132;
t97 = sin(pkin(10));
t50 = -t128 * t74 + t97 * t75;
t123 = t128 * t22 + t50 * t137 + t97 * t20;
t51 = t128 * t75 + t97 * t74;
t37 = t50 * pkin(5) - t51 * qJ(6);
t47 = -t80 * mrSges(7,1) + t51 * mrSges(7,2);
t79 = t80 ^ 2;
t124 = m(7) * (-t79 * pkin(5) + t61 * qJ(6) + 0.2e1 * qJD(6) * t80 - t50 * t37 + t123) + t80 * t47 + t61 * mrSges(7,3);
t38 = t50 * mrSges(7,1) - t51 * mrSges(7,3);
t131 = -t50 * mrSges(6,1) - t51 * mrSges(6,2) - t38;
t133 = -mrSges(6,3) - mrSges(7,2);
t30 = -t128 * t42 + t97 * t43;
t46 = t80 * mrSges(6,1) - t51 * mrSges(6,3);
t12 = m(6) * t123 - t61 * mrSges(6,2) + t131 * t50 + t133 * t30 - t80 * t46 + t124;
t114 = t128 * t20 - t97 * t22;
t135 = m(7) * (-t61 * pkin(5) - t79 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t37) * t51 - t114);
t31 = t128 * t43 + t97 * t42;
t44 = -t50 * mrSges(7,2) + t80 * mrSges(7,3);
t45 = -t80 * mrSges(6,2) - t50 * mrSges(6,3);
t13 = m(6) * t114 - t135 + (t45 + t44) * t80 + (mrSges(6,1) + mrSges(7,1)) * t61 + (m(6) * t137 + t131) * t51 + t133 * t31;
t55 = -t74 * mrSges(5,1) + t75 * mrSges(5,2);
t66 = -t80 * mrSges(5,2) + t74 * mrSges(5,3);
t10 = m(5) * t120 + t61 * mrSges(5,1) - t43 * mrSges(5,3) + t97 * t12 + t128 * t13 - t75 * t55 + t80 * t66;
t68 = t80 * mrSges(5,1) - t75 * mrSges(5,3);
t11 = m(5) * t132 - t61 * mrSges(5,2) + t42 * mrSges(5,3) + t128 * t12 - t97 * t13 + t74 * t55 - t80 * t68;
t76 = -t95 * mrSges(4,2) - t81 * mrSges(4,3);
t77 = t95 * mrSges(4,1) - t82 * mrSges(4,3);
t112 = -m(4) * t111 + t63 * mrSges(4,1) - t64 * mrSges(4,2) - t102 * t10 - t98 * t11 - t81 * t76 - t82 * t77;
t89 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t127;
t126 = qJD(1) * t104;
t90 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t126;
t138 = (t100 * t89 - t104 * t90) * qJD(1) + m(3) * (-t106 * pkin(7) + t115) - t88 * mrSges(3,1) + t87 * mrSges(3,2) - t112;
t119 = t103 * t56 - t99 * t57;
t33 = -t94 * pkin(3) - t93 * pkin(9) + t82 * t71 - t119;
t109 = -t42 * pkin(4) - t73 * qJ(5) + t75 * t67 + qJDD(5) + t33;
t113 = t31 * mrSges(7,3) + t51 * t47 - m(7) * (-0.2e1 * qJD(6) * t51 + (t50 * t80 - t31) * qJ(6) + (t51 * t80 + t30) * pkin(5) + t109) - t30 * mrSges(7,1) - t50 * t44;
t110 = m(6) * t109 + t30 * mrSges(6,1) + t31 * mrSges(6,2) + t50 * t45 + t51 * t46 - t113;
t107 = m(5) * t33 - t42 * mrSges(5,1) + t43 * mrSges(5,2) - t74 * t66 + t75 * t68 + t110;
t70 = t81 * mrSges(4,1) + t82 * mrSges(4,2);
t14 = m(4) * t119 + t94 * mrSges(4,1) - t64 * mrSges(4,3) - t82 * t70 + t95 * t76 - t107;
t7 = m(4) * t130 - t94 * mrSges(4,2) + t63 * mrSges(4,3) - t98 * t10 + t102 * t11 - t81 * t70 - t95 * t77;
t86 = (-mrSges(3,1) * t104 + mrSges(3,2) * t100) * qJD(1);
t4 = m(3) * (-t104 * g(3) - t129) - t87 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t86 * t127 + qJD(2) * t90 + t99 * t7 + t103 * t14;
t5 = m(3) * t121 - qJDD(2) * mrSges(3,2) + t88 * mrSges(3,3) - qJD(2) * t89 + t103 * t7 + t86 * t126 - t99 * t14;
t136 = t100 * t5 + t104 * t4;
t6 = m(2) * t122 + qJDD(1) * mrSges(2,1) - t106 * mrSges(2,2) - t138;
t1 = m(2) * t118 - t106 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t100 * t4 + t104 * t5;
t2 = [-m(1) * g(1) + t105 * t1 - t101 * t6, t1, t5, t7, t11, t12, -t30 * mrSges(7,2) - t50 * t38 + t124; -m(1) * g(2) + t101 * t1 + t105 * t6, t6, t4, t14, t10, t13, -t113; (-m(1) - m(2)) * g(3) + t136, -m(2) * g(3) + t136, t138, -t112, t107, t110, -t61 * mrSges(7,1) + t31 * mrSges(7,2) + t51 * t38 - t80 * t44 + t135;];
f_new  = t2;
