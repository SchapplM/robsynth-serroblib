% Calculate vector of cutting forces with Newton-Euler
% S6RRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-06 19:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRRR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 19:29:12
% EndTime: 2019-05-06 19:29:28
% DurationCPUTime: 7.23s
% Computational Cost: add. (114225->208), mult. (274888->276), div. (0->0), fcn. (209772->12), ass. (0->106)
t107 = sin(qJ(2));
t112 = cos(qJ(2));
t114 = qJD(1) ^ 2;
t101 = t112 ^ 2;
t108 = sin(qJ(1));
t113 = cos(qJ(1));
t130 = t108 * g(1) - t113 * g(2);
t123 = -qJDD(1) * pkin(1) - t130;
t134 = qJD(1) * t107;
t132 = qJD(1) * qJD(2);
t91 = t112 * qJDD(1) - t107 * t132;
t92 = qJD(2) * pkin(2) - qJ(3) * t134;
t120 = -t91 * pkin(2) + qJDD(3) + t92 * t134 + (-qJ(3) * t101 - pkin(7)) * t114 + t123;
t102 = sin(pkin(11));
t103 = cos(pkin(11));
t90 = t107 * qJDD(1) + t112 * t132;
t75 = -t102 * t90 + t103 * t91;
t85 = (t102 * t112 + t103 * t107) * qJD(1);
t79 = qJD(2) * pkin(3) - t85 * pkin(8);
t84 = (-t102 * t107 + t103 * t112) * qJD(1);
t83 = t84 ^ 2;
t118 = -t75 * pkin(3) - t83 * pkin(8) + t85 * t79 + t120;
t104 = sin(qJ(6));
t109 = cos(qJ(6));
t106 = sin(qJ(4));
t111 = cos(qJ(4));
t70 = t106 * t84 + t111 * t85;
t76 = t102 * t91 + t103 * t90;
t46 = -t70 * qJD(4) - t106 * t76 + t111 * t75;
t100 = qJD(2) + qJD(4);
t67 = t100 * pkin(4) - t70 * pkin(9);
t69 = -t106 * t85 + t111 * t84;
t68 = t69 ^ 2;
t116 = -t46 * pkin(4) - t68 * pkin(9) + t70 * t67 + t118;
t105 = sin(qJ(5));
t110 = cos(qJ(5));
t126 = -t113 * g(1) - t108 * g(2);
t87 = -t114 * pkin(1) + qJDD(1) * pkin(7) + t126;
t135 = t107 * t87;
t138 = pkin(2) * t114;
t61 = qJDD(2) * pkin(2) - t90 * qJ(3) - t135 + (qJ(3) * t132 + t107 * t138 - g(3)) * t112;
t129 = -t107 * g(3) + t112 * t87;
t62 = t91 * qJ(3) - qJD(2) * t92 - t101 * t138 + t129;
t127 = -0.2e1 * qJD(3) * t85 - t102 * t62 + t103 * t61;
t34 = (qJD(2) * t84 - t76) * pkin(8) + (t84 * t85 + qJDD(2)) * pkin(3) + t127;
t131 = 0.2e1 * qJD(3) * t84 + t102 * t61 + t103 * t62;
t36 = -t83 * pkin(3) + t75 * pkin(8) - qJD(2) * t79 + t131;
t128 = -t106 * t36 + t111 * t34;
t47 = t69 * qJD(4) + t106 * t75 + t111 * t76;
t99 = qJDD(2) + qJDD(4);
t21 = (t100 * t69 - t47) * pkin(9) + (t69 * t70 + t99) * pkin(4) + t128;
t136 = t106 * t34 + t111 * t36;
t23 = -t68 * pkin(4) + t46 * pkin(9) - t100 * t67 + t136;
t137 = t105 * t21 + t110 * t23;
t55 = -t105 * t70 + t110 * t69;
t56 = t105 * t69 + t110 * t70;
t42 = -t55 * pkin(5) - t56 * pkin(10);
t97 = qJD(5) + t100;
t95 = t97 ^ 2;
t96 = qJDD(5) + t99;
t18 = -t95 * pkin(5) + t96 * pkin(10) + t55 * t42 + t137;
t30 = -t56 * qJD(5) - t105 * t47 + t110 * t46;
t31 = t55 * qJD(5) + t105 * t46 + t110 * t47;
t19 = (-t55 * t97 - t31) * pkin(10) + (t56 * t97 - t30) * pkin(5) + t116;
t48 = -t104 * t56 + t109 * t97;
t25 = t48 * qJD(6) + t104 * t96 + t109 * t31;
t29 = qJDD(6) - t30;
t49 = t104 * t97 + t109 * t56;
t37 = -t48 * mrSges(7,1) + t49 * mrSges(7,2);
t52 = qJD(6) - t55;
t38 = -t52 * mrSges(7,2) + t48 * mrSges(7,3);
t15 = m(7) * (-t104 * t18 + t109 * t19) - t25 * mrSges(7,3) + t29 * mrSges(7,1) - t49 * t37 + t52 * t38;
t24 = -t49 * qJD(6) - t104 * t31 + t109 * t96;
t39 = t52 * mrSges(7,1) - t49 * mrSges(7,3);
t16 = m(7) * (t104 * t19 + t109 * t18) + t24 * mrSges(7,3) - t29 * mrSges(7,2) + t48 * t37 - t52 * t39;
t50 = -t97 * mrSges(6,2) + t55 * mrSges(6,3);
t51 = t97 * mrSges(6,1) - t56 * mrSges(6,3);
t122 = -m(6) * t116 + t30 * mrSges(6,1) - t31 * mrSges(6,2) - t104 * t16 - t109 * t15 + t55 * t50 - t56 * t51;
t65 = -t100 * mrSges(5,2) + t69 * mrSges(5,3);
t66 = t100 * mrSges(5,1) - t70 * mrSges(5,3);
t119 = -m(5) * t118 + t46 * mrSges(5,1) - t47 * mrSges(5,2) + t69 * t65 - t70 * t66 + t122;
t77 = -qJD(2) * mrSges(4,2) + t84 * mrSges(4,3);
t78 = qJD(2) * mrSges(4,1) - t85 * mrSges(4,3);
t117 = -m(4) * t120 + t75 * mrSges(4,1) - t76 * mrSges(4,2) + t84 * t77 - t85 * t78 + t119;
t93 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t134;
t133 = qJD(1) * t112;
t94 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t133;
t140 = (t107 * t93 - t112 * t94) * qJD(1) + m(3) * (-t114 * pkin(7) + t123) - t91 * mrSges(3,1) + t90 * mrSges(3,2) - t117;
t73 = -t84 * mrSges(4,1) + t85 * mrSges(4,2);
t41 = -t55 * mrSges(6,1) + t56 * mrSges(6,2);
t11 = m(6) * t137 - t96 * mrSges(6,2) + t30 * mrSges(6,3) - t104 * t15 + t109 * t16 + t55 * t41 - t97 * t51;
t125 = -t105 * t23 + t110 * t21;
t121 = m(7) * (-t96 * pkin(5) - t95 * pkin(10) + t56 * t42 - t125) - t24 * mrSges(7,1) + t25 * mrSges(7,2) - t48 * t38 + t49 * t39;
t12 = m(6) * t125 + t96 * mrSges(6,1) - t31 * mrSges(6,3) - t56 * t41 + t97 * t50 - t121;
t57 = -t69 * mrSges(5,1) + t70 * mrSges(5,2);
t8 = m(5) * t128 + t99 * mrSges(5,1) - t47 * mrSges(5,3) + t100 * t65 + t105 * t11 + t110 * t12 - t70 * t57;
t9 = m(5) * t136 - t99 * mrSges(5,2) + t46 * mrSges(5,3) - t100 * t66 - t105 * t12 + t110 * t11 + t69 * t57;
t6 = m(4) * t127 + qJDD(2) * mrSges(4,1) - t76 * mrSges(4,3) + qJD(2) * t77 + t106 * t9 + t111 * t8 - t85 * t73;
t7 = m(4) * t131 - qJDD(2) * mrSges(4,2) + t75 * mrSges(4,3) - qJD(2) * t78 - t106 * t8 + t111 * t9 + t84 * t73;
t89 = (-mrSges(3,1) * t112 + mrSges(3,2) * t107) * qJD(1);
t4 = m(3) * (-t112 * g(3) - t135) - t90 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t89 * t134 + qJD(2) * t94 + t102 * t7 + t103 * t6;
t5 = m(3) * t129 - qJDD(2) * mrSges(3,2) + t91 * mrSges(3,3) - qJD(2) * t93 - t102 * t6 + t103 * t7 + t89 * t133;
t139 = t107 * t5 + t112 * t4;
t10 = m(2) * t130 + qJDD(1) * mrSges(2,1) - t114 * mrSges(2,2) - t140;
t1 = m(2) * t126 - t114 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t107 * t4 + t112 * t5;
t2 = [-m(1) * g(1) + t113 * t1 - t108 * t10, t1, t5, t7, t9, t11, t16; -m(1) * g(2) + t108 * t1 + t113 * t10, t10, t4, t6, t8, t12, t15; (-m(1) - m(2)) * g(3) + t139, -m(2) * g(3) + t139, t140, -t117, -t119, -t122, t121;];
f_new  = t2;
