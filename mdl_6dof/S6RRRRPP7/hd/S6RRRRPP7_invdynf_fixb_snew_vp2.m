% Calculate vector of cutting forces with Newton-Euler
% S6RRRRPP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-05-07 18:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRPP7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP7_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP7_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP7_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP7_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP7_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:40:11
% EndTime: 2019-05-07 18:40:24
% DurationCPUTime: 5.34s
% Computational Cost: add. (95918->210), mult. (206059->282), div. (0->0), fcn. (163368->12), ass. (0->108)
t100 = sin(pkin(6));
t104 = sin(qJ(2));
t108 = cos(qJ(2));
t127 = qJD(1) * qJD(2);
t91 = (-qJDD(1) * t108 + t104 * t127) * t100;
t143 = -2 * qJD(5);
t133 = cos(pkin(11));
t102 = sin(qJ(4));
t106 = cos(qJ(4));
t103 = sin(qJ(3));
t107 = cos(qJ(3));
t129 = qJD(1) * t108;
t101 = cos(pkin(6));
t110 = qJD(1) ^ 2;
t105 = sin(qJ(1));
t109 = cos(qJ(1));
t122 = t105 * g(1) - t109 * g(2);
t141 = pkin(8) * t100;
t86 = qJDD(1) * pkin(1) + t110 * t141 + t122;
t134 = t101 * t86;
t119 = -t109 * g(1) - t105 * g(2);
t87 = -t110 * pkin(1) + qJDD(1) * t141 + t119;
t135 = t104 * t134 + t108 * t87;
t130 = qJD(1) * t100;
t89 = (-pkin(2) * t108 - pkin(9) * t104) * t130;
t96 = t101 * qJD(1) + qJD(2);
t94 = t96 ^ 2;
t95 = t101 * qJDD(1) + qJDD(2);
t53 = -t94 * pkin(2) + t95 * pkin(9) + (-g(3) * t104 + t89 * t129) * t100 + t135;
t140 = t101 * g(3);
t90 = (qJDD(1) * t104 + t108 * t127) * t100;
t54 = t91 * pkin(2) - t90 * pkin(9) - t140 + (-t86 + (pkin(2) * t104 - pkin(9) * t108) * t96 * qJD(1)) * t100;
t136 = t103 * t54 + t107 * t53;
t124 = t104 * t130;
t79 = -t103 * t124 + t107 * t96;
t80 = t103 * t96 + t107 * t124;
t68 = -t79 * pkin(3) - t80 * pkin(10);
t83 = qJDD(3) + t91;
t123 = t100 * t129;
t93 = qJD(3) - t123;
t92 = t93 ^ 2;
t27 = -t92 * pkin(3) + t83 * pkin(10) + t79 * t68 + t136;
t131 = t100 * t108;
t118 = -g(3) * t131 - t104 * t87 + t108 * t134;
t52 = -t95 * pkin(2) - t94 * pkin(9) + t89 * t124 - t118;
t65 = -t80 * qJD(3) - t103 * t90 + t107 * t95;
t66 = t79 * qJD(3) + t103 * t95 + t107 * t90;
t30 = (-t79 * t93 - t66) * pkin(10) + (t80 * t93 - t65) * pkin(3) + t52;
t121 = -t102 * t27 + t106 * t30;
t71 = -t102 * t80 + t106 * t93;
t43 = t71 * qJD(4) + t102 * t83 + t106 * t66;
t64 = qJDD(4) - t65;
t72 = t102 * t93 + t106 * t80;
t78 = qJD(4) - t79;
t20 = (t71 * t78 - t43) * qJ(5) + (t71 * t72 + t64) * pkin(4) + t121;
t138 = t102 * t30 + t106 * t27;
t42 = -t72 * qJD(4) - t102 * t66 + t106 * t83;
t61 = t78 * pkin(4) - t72 * qJ(5);
t70 = t71 ^ 2;
t22 = -t70 * pkin(4) + t42 * qJ(5) - t78 * t61 + t138;
t99 = sin(pkin(11));
t117 = t133 * t20 - t99 * t22;
t56 = -t133 * t71 + t99 * t72;
t57 = t133 * t72 + t99 * t71;
t37 = t56 * pkin(5) - t57 * qJ(6);
t77 = t78 ^ 2;
t142 = m(7) * (-t64 * pkin(5) - t77 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t37) * t57 - t117);
t139 = -mrSges(6,3) - mrSges(7,2);
t38 = t56 * mrSges(7,1) - t57 * mrSges(7,3);
t137 = -t56 * mrSges(6,1) - t57 * mrSges(6,2) - t38;
t132 = t100 * t104;
t120 = -t103 * t53 + t107 * t54;
t26 = -t83 * pkin(3) - t92 * pkin(10) + t80 * t68 - t120;
t113 = -t42 * pkin(4) - t70 * qJ(5) + t72 * t61 + qJDD(5) + t26;
t33 = -t133 * t42 + t99 * t43;
t34 = t133 * t43 + t99 * t42;
t44 = -t56 * mrSges(7,2) + t78 * mrSges(7,3);
t47 = -t78 * mrSges(7,1) + t57 * mrSges(7,2);
t115 = t34 * mrSges(7,3) + t57 * t47 - m(7) * (-0.2e1 * qJD(6) * t57 + (t56 * t78 - t34) * qJ(6) + (t57 * t78 + t33) * pkin(5) + t113) - t33 * mrSges(7,1) - t56 * t44;
t45 = -t78 * mrSges(6,2) - t56 * mrSges(6,3);
t46 = t78 * mrSges(6,1) - t57 * mrSges(6,3);
t114 = m(6) * t113 + t33 * mrSges(6,1) + t34 * mrSges(6,2) + t56 * t45 + t57 * t46 - t115;
t60 = -t78 * mrSges(5,2) + t71 * mrSges(5,3);
t62 = t78 * mrSges(5,1) - t72 * mrSges(5,3);
t111 = m(5) * t26 - t42 * mrSges(5,1) + t43 * mrSges(5,2) - t71 * t60 + t72 * t62 + t114;
t67 = -t79 * mrSges(4,1) + t80 * mrSges(4,2);
t73 = -t93 * mrSges(4,2) + t79 * mrSges(4,3);
t14 = m(4) * t120 + t83 * mrSges(4,1) - t66 * mrSges(4,3) - t80 * t67 + t93 * t73 - t111;
t84 = t96 * mrSges(3,1) - mrSges(3,3) * t124;
t88 = (-mrSges(3,1) * t108 + mrSges(3,2) * t104) * t130;
t125 = t133 * t22 + t56 * t143 + t99 * t20;
t126 = m(7) * (-t77 * pkin(5) + t64 * qJ(6) + 0.2e1 * qJD(6) * t78 - t56 * t37 + t125) + t78 * t47 + t64 * mrSges(7,3);
t12 = m(6) * t125 - t64 * mrSges(6,2) + t137 * t56 + t139 * t33 - t78 * t46 + t126;
t13 = m(6) * t117 - t142 + (t45 + t44) * t78 + (mrSges(6,1) + mrSges(7,1)) * t64 + (m(6) * t143 + t137) * t57 + t139 * t34;
t58 = -t71 * mrSges(5,1) + t72 * mrSges(5,2);
t10 = m(5) * t121 + t64 * mrSges(5,1) - t43 * mrSges(5,3) + t99 * t12 + t133 * t13 - t72 * t58 + t78 * t60;
t11 = m(5) * t138 - t64 * mrSges(5,2) + t42 * mrSges(5,3) + t133 * t12 - t99 * t13 + t71 * t58 - t78 * t62;
t74 = t93 * mrSges(4,1) - t80 * mrSges(4,3);
t9 = m(4) * t136 - t83 * mrSges(4,2) + t65 * mrSges(4,3) - t102 * t10 + t106 * t11 + t79 * t67 - t93 * t74;
t4 = m(3) * (-g(3) * t132 + t135) - t91 * mrSges(3,3) - t95 * mrSges(3,2) + t88 * t123 - t96 * t84 + t107 * t9 - t103 * t14;
t85 = -t96 * mrSges(3,2) + mrSges(3,3) * t123;
t6 = m(3) * (-t100 * t86 - t140) + t90 * mrSges(3,2) + t91 * mrSges(3,1) + t103 * t9 + t107 * t14 + (t104 * t84 - t108 * t85) * t130;
t112 = m(4) * t52 - t65 * mrSges(4,1) + t66 * mrSges(4,2) + t106 * t10 + t102 * t11 - t79 * t73 + t80 * t74;
t8 = m(3) * t118 + t95 * mrSges(3,1) - t90 * mrSges(3,3) - t88 * t124 + t96 * t85 - t112;
t128 = t101 * t6 + t8 * t131 + t4 * t132;
t2 = m(2) * t119 - t110 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t104 * t8 + t108 * t4;
t1 = m(2) * t122 + qJDD(1) * mrSges(2,1) - t110 * mrSges(2,2) - t100 * t6 + (t104 * t4 + t108 * t8) * t101;
t3 = [-m(1) * g(1) - t105 * t1 + t109 * t2, t2, t4, t9, t11, t12, -t33 * mrSges(7,2) - t56 * t38 + t126; -m(1) * g(2) + t109 * t1 + t105 * t2, t1, t8, t14, t10, t13, -t115; (-m(1) - m(2)) * g(3) + t128, -m(2) * g(3) + t128, t6, t112, t111, t114, -t64 * mrSges(7,1) + t34 * mrSges(7,2) + t57 * t38 - t78 * t44 + t142;];
f_new  = t3;
