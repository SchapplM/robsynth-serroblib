% Calculate vector of cutting forces with Newton-Euler
% S6RRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-05-06 12:52
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRPR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:47:21
% EndTime: 2019-05-06 12:47:34
% DurationCPUTime: 5.66s
% Computational Cost: add. (95052->206), mult. (223328->278), div. (0->0), fcn. (167537->12), ass. (0->104)
t111 = sin(qJ(2));
t114 = cos(qJ(2));
t116 = qJD(1) ^ 2;
t104 = t114 ^ 2;
t112 = sin(qJ(1));
t115 = cos(qJ(1));
t130 = g(1) * t112 - t115 * g(2);
t124 = -qJDD(1) * pkin(1) - t130;
t136 = qJD(1) * t111;
t134 = qJD(1) * qJD(2);
t96 = qJDD(1) * t114 - t111 * t134;
t97 = qJD(2) * pkin(2) - qJ(3) * t136;
t121 = -pkin(2) * t96 + qJDD(3) + t97 * t136 + (-qJ(3) * t104 - pkin(7)) * t116 + t124;
t105 = sin(pkin(11));
t107 = cos(pkin(11));
t106 = sin(pkin(10));
t108 = cos(pkin(10));
t95 = qJDD(1) * t111 + t114 * t134;
t79 = -t106 * t95 + t108 * t96;
t90 = (t106 * t114 + t108 * t111) * qJD(1);
t83 = qJD(2) * pkin(3) - pkin(8) * t90;
t89 = (-t106 * t111 + t108 * t114) * qJD(1);
t88 = t89 ^ 2;
t118 = -pkin(3) * t79 - pkin(8) * t88 + t90 * t83 + t121;
t109 = sin(qJ(6));
t113 = cos(qJ(6));
t103 = qJD(2) + qJD(4);
t101 = t103 ^ 2;
t102 = qJDD(2) + qJDD(4);
t110 = sin(qJ(4));
t140 = cos(qJ(4));
t126 = -g(1) * t115 - g(2) * t112;
t92 = -pkin(1) * t116 + qJDD(1) * pkin(7) + t126;
t137 = t111 * t92;
t139 = pkin(2) * t116;
t60 = qJDD(2) * pkin(2) - qJ(3) * t95 - t137 + (qJ(3) * t134 + t111 * t139 - g(3)) * t114;
t131 = -g(3) * t111 + t114 * t92;
t61 = qJ(3) * t96 - qJD(2) * t97 - t104 * t139 + t131;
t128 = -0.2e1 * qJD(3) * t90 - t106 * t61 + t108 * t60;
t80 = t106 * t96 + t108 * t95;
t32 = (qJD(2) * t89 - t80) * pkin(8) + (t89 * t90 + qJDD(2)) * pkin(3) + t128;
t132 = 0.2e1 * qJD(3) * t89 + t106 * t60 + t108 * t61;
t36 = -pkin(3) * t88 + pkin(8) * t79 - qJD(2) * t83 + t132;
t138 = t110 * t32 + t140 * t36;
t72 = t110 * t90 - t140 * t89;
t73 = t110 * t89 + t140 * t90;
t55 = pkin(4) * t72 - qJ(5) * t73;
t23 = -pkin(4) * t101 + qJ(5) * t102 - t55 * t72 + t138;
t46 = qJD(4) * t73 + t110 * t80 - t140 * t79;
t47 = -t72 * qJD(4) + t110 * t79 + t140 * t80;
t26 = (t103 * t72 - t47) * qJ(5) + (t103 * t73 + t46) * pkin(4) + t118;
t67 = t103 * t105 + t107 * t73;
t129 = -0.2e1 * qJD(5) * t67 - t105 * t23 + t107 * t26;
t42 = t102 * t105 + t107 * t47;
t66 = t103 * t107 - t105 * t73;
t17 = (t66 * t72 - t42) * pkin(9) + (t66 * t67 + t46) * pkin(5) + t129;
t133 = 0.2e1 * qJD(5) * t66 + t105 * t26 + t107 * t23;
t41 = t102 * t107 - t105 * t47;
t53 = pkin(5) * t72 - t67 * pkin(9);
t65 = t66 ^ 2;
t18 = -t65 * pkin(5) + t41 * pkin(9) - t53 * t72 + t133;
t48 = -t109 * t67 + t113 * t66;
t29 = t48 * qJD(6) + t109 * t41 + t113 * t42;
t49 = t109 * t66 + t113 * t67;
t33 = -mrSges(7,1) * t48 + mrSges(7,2) * t49;
t71 = qJD(6) + t72;
t37 = -mrSges(7,2) * t71 + t48 * mrSges(7,3);
t45 = qJDD(6) + t46;
t15 = m(7) * (-t109 * t18 + t113 * t17) - t29 * mrSges(7,3) + t45 * mrSges(7,1) - t49 * t33 + t71 * t37;
t28 = -t49 * qJD(6) - t109 * t42 + t113 * t41;
t38 = mrSges(7,1) * t71 - t49 * mrSges(7,3);
t16 = m(7) * (t109 * t17 + t113 * t18) + t28 * mrSges(7,3) - t45 * mrSges(7,2) + t48 * t33 - t71 * t38;
t50 = -mrSges(6,1) * t66 + mrSges(6,2) * t67;
t51 = -mrSges(6,2) * t72 + t66 * mrSges(6,3);
t12 = m(6) * t129 + t46 * mrSges(6,1) - t42 * mrSges(6,3) + t109 * t16 + t113 * t15 - t67 * t50 + t72 * t51;
t52 = mrSges(6,1) * t72 - t67 * mrSges(6,3);
t13 = m(6) * t133 - t46 * mrSges(6,2) + t41 * mrSges(6,3) - t109 * t15 + t113 * t16 + t66 * t50 - t72 * t52;
t69 = -mrSges(5,2) * t103 - mrSges(5,3) * t72;
t70 = mrSges(5,1) * t103 - mrSges(5,3) * t73;
t122 = m(5) * t118 + t46 * mrSges(5,1) + t47 * mrSges(5,2) + t105 * t13 + t107 * t12 + t72 * t69 + t73 * t70;
t81 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t89;
t82 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t90;
t120 = -m(4) * t121 + t79 * mrSges(4,1) - t80 * mrSges(4,2) + t89 * t81 - t90 * t82 - t122;
t98 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t136;
t135 = qJD(1) * t114;
t99 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t135;
t142 = (t111 * t98 - t114 * t99) * qJD(1) + m(3) * (-pkin(7) * t116 + t124) - t96 * mrSges(3,1) + t95 * mrSges(3,2) - t120;
t127 = -t110 * t36 + t140 * t32;
t22 = -t102 * pkin(4) - t101 * qJ(5) + t73 * t55 + qJDD(5) - t127;
t123 = t28 * mrSges(7,1) + t48 * t37 - m(7) * (-t41 * pkin(5) - t65 * pkin(9) + t67 * t53 + t22) - t29 * mrSges(7,2) - t49 * t38;
t119 = m(6) * t22 - t41 * mrSges(6,1) + t42 * mrSges(6,2) - t66 * t51 + t67 * t52 - t123;
t56 = mrSges(5,1) * t72 + mrSges(5,2) * t73;
t14 = m(5) * t127 + t102 * mrSges(5,1) - t47 * mrSges(5,3) + t103 * t69 - t73 * t56 - t119;
t76 = -mrSges(4,1) * t89 + mrSges(4,2) * t90;
t9 = m(5) * t138 - t102 * mrSges(5,2) - t46 * mrSges(5,3) - t103 * t70 - t105 * t12 + t107 * t13 - t72 * t56;
t6 = m(4) * t128 + qJDD(2) * mrSges(4,1) - t80 * mrSges(4,3) + qJD(2) * t81 + t110 * t9 + t140 * t14 - t90 * t76;
t7 = m(4) * t132 - qJDD(2) * mrSges(4,2) + t79 * mrSges(4,3) - qJD(2) * t82 - t110 * t14 + t140 * t9 + t89 * t76;
t94 = (-mrSges(3,1) * t114 + mrSges(3,2) * t111) * qJD(1);
t4 = m(3) * (-g(3) * t114 - t137) - t95 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t94 * t136 + qJD(2) * t99 + t106 * t7 + t108 * t6;
t5 = m(3) * t131 - qJDD(2) * mrSges(3,2) + t96 * mrSges(3,3) - qJD(2) * t98 - t106 * t6 + t108 * t7 + t94 * t135;
t141 = t111 * t5 + t114 * t4;
t8 = m(2) * t130 + qJDD(1) * mrSges(2,1) - t116 * mrSges(2,2) - t142;
t1 = m(2) * t126 - t116 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t111 * t4 + t114 * t5;
t2 = [-m(1) * g(1) + t1 * t115 - t112 * t8, t1, t5, t7, t9, t13, t16; -m(1) * g(2) + t1 * t112 + t115 * t8, t8, t4, t6, t14, t12, t15; (-m(1) - m(2)) * g(3) + t141, -m(2) * g(3) + t141, t142, -t120, t122, t119, -t123;];
f_new  = t2;
