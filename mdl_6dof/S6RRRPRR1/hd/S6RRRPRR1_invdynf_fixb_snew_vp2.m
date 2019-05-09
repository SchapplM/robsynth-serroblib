% Calculate vector of cutting forces with Newton-Euler
% S6RRRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 09:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPRR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:39:56
% EndTime: 2019-05-07 09:40:14
% DurationCPUTime: 7.16s
% Computational Cost: add. (119849->208), mult. (279724->276), div. (0->0), fcn. (212222->12), ass. (0->106)
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
t94 = qJD(2) * pkin(2) - pkin(8) * t134;
t120 = -t91 * pkin(2) + t94 * t134 + (-pkin(8) * t101 - pkin(7)) * t114 + t123;
t106 = sin(qJ(3));
t111 = cos(qJ(3));
t85 = (t106 * t112 + t107 * t111) * qJD(1);
t90 = t107 * qJDD(1) + t112 * t132;
t65 = -t85 * qJD(3) - t106 * t90 + t111 * t91;
t100 = qJD(2) + qJD(3);
t80 = t100 * pkin(3) - t85 * qJ(4);
t84 = (-t106 * t107 + t111 * t112) * qJD(1);
t83 = t84 ^ 2;
t118 = -t65 * pkin(3) - t83 * qJ(4) + t85 * t80 + qJDD(4) + t120;
t104 = sin(qJ(6));
t109 = cos(qJ(6));
t102 = sin(pkin(11));
t103 = cos(pkin(11));
t66 = t84 * qJD(3) + t106 * t91 + t111 * t90;
t46 = -t102 * t66 + t103 * t65;
t77 = t102 * t84 + t103 * t85;
t70 = t100 * pkin(4) - t77 * pkin(9);
t76 = -t102 * t85 + t103 * t84;
t73 = t76 ^ 2;
t116 = -t46 * pkin(4) - t73 * pkin(9) + t77 * t70 + t118;
t105 = sin(qJ(5));
t110 = cos(qJ(5));
t126 = -t113 * g(1) - t108 * g(2);
t87 = -t114 * pkin(1) + qJDD(1) * pkin(7) + t126;
t135 = t107 * t87;
t138 = pkin(2) * t114;
t61 = qJDD(2) * pkin(2) - t90 * pkin(8) - t135 + (pkin(8) * t132 + t107 * t138 - g(3)) * t112;
t129 = -t107 * g(3) + t112 * t87;
t62 = t91 * pkin(8) - qJD(2) * t94 - t101 * t138 + t129;
t128 = -t106 * t62 + t111 * t61;
t99 = qJDD(2) + qJDD(3);
t34 = (t100 * t84 - t66) * qJ(4) + (t84 * t85 + t99) * pkin(3) + t128;
t136 = t106 * t61 + t111 * t62;
t36 = -t83 * pkin(3) + t65 * qJ(4) - t100 * t80 + t136;
t127 = -0.2e1 * qJD(4) * t77 - t102 * t36 + t103 * t34;
t47 = t102 * t65 + t103 * t66;
t21 = (t100 * t76 - t47) * pkin(9) + (t76 * t77 + t99) * pkin(4) + t127;
t131 = 0.2e1 * qJD(4) * t76 + t102 * t34 + t103 * t36;
t23 = -t73 * pkin(4) + t46 * pkin(9) - t100 * t70 + t131;
t137 = t105 * t21 + t110 * t23;
t55 = -t105 * t77 + t110 * t76;
t56 = t105 * t76 + t110 * t77;
t42 = -t55 * pkin(5) - t56 * pkin(10);
t97 = qJD(5) + t100;
t95 = t97 ^ 2;
t96 = qJDD(5) + t99;
t18 = -t95 * pkin(5) + t96 * pkin(10) + t55 * t42 + t137;
t30 = -t56 * qJD(5) - t105 * t47 + t110 * t46;
t31 = t55 * qJD(5) + t105 * t46 + t110 * t47;
t19 = t116 + (-t55 * t97 - t31) * pkin(10) + (t56 * t97 - t30) * pkin(5);
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
t68 = -t100 * mrSges(5,2) + t76 * mrSges(5,3);
t69 = t100 * mrSges(5,1) - t77 * mrSges(5,3);
t119 = -m(5) * t118 + t46 * mrSges(5,1) - t47 * mrSges(5,2) + t76 * t68 - t77 * t69 + t122;
t79 = -t100 * mrSges(4,2) + t84 * mrSges(4,3);
t81 = t100 * mrSges(4,1) - t85 * mrSges(4,3);
t117 = -m(4) * t120 + t65 * mrSges(4,1) - t66 * mrSges(4,2) + t84 * t79 - t85 * t81 + t119;
t92 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t134;
t133 = qJD(1) * t112;
t93 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t133;
t140 = (t107 * t92 - t112 * t93) * qJD(1) + m(3) * (-t114 * pkin(7) + t123) - t91 * mrSges(3,1) + t90 * mrSges(3,2) - t117;
t78 = -t84 * mrSges(4,1) + t85 * mrSges(4,2);
t41 = -t55 * mrSges(6,1) + t56 * mrSges(6,2);
t11 = m(6) * t137 - t96 * mrSges(6,2) + t30 * mrSges(6,3) - t104 * t15 + t109 * t16 + t55 * t41 - t97 * t51;
t125 = -t105 * t23 + t110 * t21;
t121 = m(7) * (-t96 * pkin(5) - t95 * pkin(10) + t56 * t42 - t125) - t24 * mrSges(7,1) + t25 * mrSges(7,2) - t48 * t38 + t49 * t39;
t12 = m(6) * t125 + t96 * mrSges(6,1) - t31 * mrSges(6,3) - t56 * t41 + t97 * t50 - t121;
t57 = -t76 * mrSges(5,1) + t77 * mrSges(5,2);
t8 = m(5) * t127 + t99 * mrSges(5,1) - t47 * mrSges(5,3) + t100 * t68 + t105 * t11 + t110 * t12 - t77 * t57;
t9 = m(5) * t131 - t99 * mrSges(5,2) + t46 * mrSges(5,3) - t100 * t69 - t105 * t12 + t110 * t11 + t76 * t57;
t6 = m(4) * t128 + t99 * mrSges(4,1) - t66 * mrSges(4,3) + t100 * t79 + t102 * t9 + t103 * t8 - t85 * t78;
t7 = m(4) * t136 - t99 * mrSges(4,2) + t65 * mrSges(4,3) - t100 * t81 - t102 * t8 + t103 * t9 + t84 * t78;
t89 = (-mrSges(3,1) * t112 + mrSges(3,2) * t107) * qJD(1);
t4 = m(3) * (-t112 * g(3) - t135) - t90 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t89 * t134 + qJD(2) * t93 + t106 * t7 + t111 * t6;
t5 = m(3) * t129 - qJDD(2) * mrSges(3,2) + t91 * mrSges(3,3) - qJD(2) * t92 - t106 * t6 + t111 * t7 + t133 * t89;
t139 = t107 * t5 + t112 * t4;
t10 = m(2) * t130 + qJDD(1) * mrSges(2,1) - t114 * mrSges(2,2) - t140;
t1 = m(2) * t126 - t114 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t107 * t4 + t112 * t5;
t2 = [-m(1) * g(1) + t113 * t1 - t108 * t10, t1, t5, t7, t9, t11, t16; -m(1) * g(2) + t108 * t1 + t113 * t10, t10, t4, t6, t8, t12, t15; (-m(1) - m(2)) * g(3) + t139, -m(2) * g(3) + t139, t140, -t117, -t119, -t122, t121;];
f_new  = t2;
