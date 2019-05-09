% Calculate vector of cutting forces with Newton-Euler
% S6RRRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 09:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRRR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 08:51:35
% EndTime: 2019-05-08 08:51:55
% DurationCPUTime: 6.73s
% Computational Cost: add. (131092->209), mult. (261359->273), div. (0->0), fcn. (194250->12), ass. (0->108)
t110 = sin(qJ(2));
t116 = cos(qJ(2));
t118 = qJD(1) ^ 2;
t107 = sin(qJ(5));
t113 = cos(qJ(5));
t106 = sin(qJ(6));
t112 = cos(qJ(6));
t108 = sin(qJ(4));
t114 = cos(qJ(4));
t104 = qJD(2) + qJD(3);
t105 = t116 ^ 2;
t111 = sin(qJ(1));
t117 = cos(qJ(1));
t132 = t111 * g(1) - t117 * g(2);
t126 = -qJDD(1) * pkin(1) - t132;
t136 = qJD(1) * t110;
t134 = qJD(1) * qJD(2);
t96 = t116 * qJDD(1) - t110 * t134;
t99 = qJD(2) * pkin(2) - pkin(8) * t136;
t122 = -t96 * pkin(2) + t99 * t136 + (-pkin(8) * t105 - pkin(7)) * t118 + t126;
t109 = sin(qJ(3));
t115 = cos(qJ(3));
t89 = (t109 * t116 + t110 * t115) * qJD(1);
t95 = t110 * qJDD(1) + t116 * t134;
t69 = -t89 * qJD(3) - t109 * t95 + t115 * t96;
t135 = qJD(1) * t116;
t88 = -t109 * t136 + t115 * t135;
t70 = t88 * qJD(3) + t109 * t96 + t115 * t95;
t38 = (-t104 * t88 - t70) * pkin(9) + (t104 * t89 - t69) * pkin(3) + t122;
t102 = t104 ^ 2;
t103 = qJDD(2) + qJDD(3);
t128 = -t117 * g(1) - t111 * g(2);
t91 = -t118 * pkin(1) + qJDD(1) * pkin(7) + t128;
t137 = t110 * t91;
t141 = pkin(2) * t118;
t61 = qJDD(2) * pkin(2) - t95 * pkin(8) - t137 + (pkin(8) * t134 + t110 * t141 - g(3)) * t116;
t133 = -t110 * g(3) + t116 * t91;
t62 = t96 * pkin(8) - qJD(2) * t99 - t105 * t141 + t133;
t138 = t109 * t61 + t115 * t62;
t77 = -t88 * pkin(3) - t89 * pkin(9);
t41 = -t102 * pkin(3) + t103 * pkin(9) + t88 * t77 + t138;
t130 = -t108 * t41 + t114 * t38;
t79 = t114 * t104 - t108 * t89;
t49 = t79 * qJD(4) + t108 * t103 + t114 * t70;
t67 = qJDD(4) - t69;
t80 = t108 * t104 + t114 * t89;
t87 = qJD(4) - t88;
t26 = (t79 * t87 - t49) * pkin(10) + (t79 * t80 + t67) * pkin(4) + t130;
t139 = t108 * t38 + t114 * t41;
t48 = -t80 * qJD(4) + t114 * t103 - t108 * t70;
t74 = t87 * pkin(4) - t80 * pkin(10);
t78 = t79 ^ 2;
t28 = -t78 * pkin(4) + t48 * pkin(10) - t87 * t74 + t139;
t131 = -t107 * t28 + t113 * t26;
t55 = -t107 * t80 + t113 * t79;
t35 = t55 * qJD(5) + t107 * t48 + t113 * t49;
t56 = t107 * t79 + t113 * t80;
t64 = qJDD(5) + t67;
t85 = qJD(5) + t87;
t17 = (t55 * t85 - t35) * pkin(11) + (t55 * t56 + t64) * pkin(5) + t131;
t140 = t107 * t26 + t113 * t28;
t34 = -t56 * qJD(5) - t107 * t49 + t113 * t48;
t52 = t85 * pkin(5) - t56 * pkin(11);
t54 = t55 ^ 2;
t18 = -t54 * pkin(5) + t34 * pkin(11) - t85 * t52 + t140;
t45 = -t106 * t56 + t112 * t55;
t23 = t45 * qJD(6) + t106 * t34 + t112 * t35;
t46 = t106 * t55 + t112 * t56;
t32 = -t45 * mrSges(7,1) + t46 * mrSges(7,2);
t83 = qJD(6) + t85;
t42 = -t83 * mrSges(7,2) + t45 * mrSges(7,3);
t63 = qJDD(6) + t64;
t15 = m(7) * (-t106 * t18 + t112 * t17) - t23 * mrSges(7,3) + t63 * mrSges(7,1) - t46 * t32 + t83 * t42;
t22 = -t46 * qJD(6) - t106 * t35 + t112 * t34;
t43 = t83 * mrSges(7,1) - t46 * mrSges(7,3);
t16 = m(7) * (t106 * t17 + t112 * t18) + t22 * mrSges(7,3) - t63 * mrSges(7,2) + t45 * t32 - t83 * t43;
t47 = -t55 * mrSges(6,1) + t56 * mrSges(6,2);
t50 = -t85 * mrSges(6,2) + t55 * mrSges(6,3);
t12 = m(6) * t131 + t64 * mrSges(6,1) - t35 * mrSges(6,3) + t106 * t16 + t112 * t15 - t56 * t47 + t85 * t50;
t51 = t85 * mrSges(6,1) - t56 * mrSges(6,3);
t13 = m(6) * t140 - t64 * mrSges(6,2) + t34 * mrSges(6,3) - t106 * t15 + t112 * t16 + t55 * t47 - t85 * t51;
t60 = -t79 * mrSges(5,1) + t80 * mrSges(5,2);
t72 = -t87 * mrSges(5,2) + t79 * mrSges(5,3);
t10 = m(5) * t130 + t67 * mrSges(5,1) - t49 * mrSges(5,3) + t107 * t13 + t113 * t12 - t80 * t60 + t87 * t72;
t73 = t87 * mrSges(5,1) - t80 * mrSges(5,3);
t11 = m(5) * t139 - t67 * mrSges(5,2) + t48 * mrSges(5,3) - t107 * t12 + t113 * t13 + t79 * t60 - t87 * t73;
t81 = -t104 * mrSges(4,2) + t88 * mrSges(4,3);
t82 = t104 * mrSges(4,1) - t89 * mrSges(4,3);
t124 = -m(4) * t122 + t69 * mrSges(4,1) - t70 * mrSges(4,2) - t114 * t10 - t108 * t11 + t88 * t81 - t89 * t82;
t97 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t136;
t98 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t135;
t143 = (t110 * t97 - t116 * t98) * qJD(1) + m(3) * (-t118 * pkin(7) + t126) - t96 * mrSges(3,1) + t95 * mrSges(3,2) - t124;
t129 = -t109 * t62 + t115 * t61;
t40 = -t103 * pkin(3) - t102 * pkin(9) + t89 * t77 - t129;
t123 = -t48 * pkin(4) - t78 * pkin(10) + t80 * t74 + t40;
t125 = t22 * mrSges(7,1) + t45 * t42 - m(7) * (-t34 * pkin(5) - t54 * pkin(11) + t56 * t52 + t123) - t23 * mrSges(7,2) - t46 * t43;
t121 = -m(6) * t123 + t34 * mrSges(6,1) - t35 * mrSges(6,2) + t55 * t50 - t56 * t51 + t125;
t119 = m(5) * t40 - t48 * mrSges(5,1) + t49 * mrSges(5,2) - t79 * t72 + t80 * t73 - t121;
t76 = -t88 * mrSges(4,1) + t89 * mrSges(4,2);
t14 = m(4) * t129 + t103 * mrSges(4,1) - t70 * mrSges(4,3) + t104 * t81 - t89 * t76 - t119;
t7 = m(4) * t138 - t103 * mrSges(4,2) + t69 * mrSges(4,3) - t108 * t10 - t104 * t82 + t114 * t11 + t88 * t76;
t94 = (-mrSges(3,1) * t116 + mrSges(3,2) * t110) * qJD(1);
t4 = m(3) * (-t116 * g(3) - t137) - t95 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t94 * t136 + qJD(2) * t98 + t109 * t7 + t115 * t14;
t5 = m(3) * t133 - qJDD(2) * mrSges(3,2) + t96 * mrSges(3,3) - qJD(2) * t97 - t109 * t14 + t115 * t7 + t94 * t135;
t142 = t110 * t5 + t116 * t4;
t6 = m(2) * t132 + qJDD(1) * mrSges(2,1) - t118 * mrSges(2,2) - t143;
t1 = m(2) * t128 - t118 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t110 * t4 + t116 * t5;
t2 = [-m(1) * g(1) + t117 * t1 - t111 * t6, t1, t5, t7, t11, t13, t16; -m(1) * g(2) + t111 * t1 + t117 * t6, t6, t4, t14, t10, t12, t15; (-m(1) - m(2)) * g(3) + t142, -m(2) * g(3) + t142, t143, -t124, t119, -t121, -t125;];
f_new  = t2;
