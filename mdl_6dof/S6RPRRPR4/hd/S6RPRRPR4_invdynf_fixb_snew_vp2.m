% Calculate vector of cutting forces with Newton-Euler
% S6RPRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-05-05 22:31
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRPR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:26:45
% EndTime: 2019-05-05 22:26:56
% DurationCPUTime: 5.31s
% Computational Cost: add. (88678->194), mult. (216093->255), div. (0->0), fcn. (171451->12), ass. (0->104)
t111 = qJD(1) ^ 2;
t101 = sin(pkin(10));
t103 = cos(pkin(10));
t98 = t103 ^ 2;
t135 = t101 ^ 2 + t98;
t141 = t135 * mrSges(3,3);
t106 = sin(qJ(3));
t109 = cos(qJ(3));
t122 = -t103 * mrSges(3,1) + t101 * mrSges(3,2);
t119 = qJDD(1) * mrSges(3,3) + t111 * t122;
t132 = qJD(1) * qJD(2);
t129 = -t103 * g(3) - 0.2e1 * t101 * t132;
t105 = sin(qJ(4));
t134 = pkin(7) * qJDD(1);
t138 = pkin(2) * t111;
t107 = sin(qJ(1));
t110 = cos(qJ(1));
t123 = -t110 * g(1) - t107 * g(2);
t90 = -t111 * pkin(1) + qJDD(1) * qJ(2) + t123;
t67 = (t103 * t138 - t134 - t90) * t101 + t129;
t127 = -t101 * g(3) + (0.2e1 * t132 + t90) * t103;
t68 = t103 * t134 - t98 * t138 + t127;
t128 = -t106 * t68 + t109 * t67;
t139 = cos(qJ(4));
t120 = -t101 * t106 + t103 * t109;
t88 = t120 * qJD(1);
t133 = t88 * qJD(3);
t121 = t101 * t109 + t103 * t106;
t80 = t121 * qJDD(1) + t133;
t89 = t121 * qJD(1);
t33 = (-t80 + t133) * pkin(8) + (t88 * t89 + qJDD(3)) * pkin(3) + t128;
t136 = t106 * t67 + t109 * t68;
t79 = -t89 * qJD(3) + t120 * qJDD(1);
t83 = qJD(3) * pkin(3) - t89 * pkin(8);
t87 = t88 ^ 2;
t36 = -t87 * pkin(3) + t79 * pkin(8) - qJD(3) * t83 + t136;
t125 = -t105 * t36 + t139 * t33;
t70 = t105 * t89 - t139 * t88;
t71 = t105 * t88 + t139 * t89;
t55 = t70 * pkin(4) - t71 * qJ(5);
t99 = qJD(3) + qJD(4);
t95 = t99 ^ 2;
t96 = qJDD(3) + qJDD(4);
t22 = -t96 * pkin(4) - t95 * qJ(5) + t71 * t55 + qJDD(5) - t125;
t104 = sin(qJ(6));
t108 = cos(qJ(6));
t100 = sin(pkin(11));
t102 = cos(pkin(11));
t47 = -t70 * qJD(4) + t105 * t79 + t139 * t80;
t39 = -t100 * t47 + t102 * t96;
t40 = t100 * t96 + t102 * t47;
t60 = -t100 * t71 + t102 * t99;
t61 = t100 * t99 + t102 * t71;
t49 = t104 * t60 + t108 * t61;
t28 = -t49 * qJD(6) - t104 * t40 + t108 * t39;
t48 = -t104 * t61 + t108 * t60;
t29 = t48 * qJD(6) + t104 * t39 + t108 * t40;
t69 = qJD(6) + t70;
t37 = -t69 * mrSges(7,2) + t48 * mrSges(7,3);
t38 = t69 * mrSges(7,1) - t49 * mrSges(7,3);
t53 = t70 * pkin(5) - t61 * pkin(9);
t59 = t60 ^ 2;
t118 = t28 * mrSges(7,1) + t48 * t37 - m(7) * (-t39 * pkin(5) - t59 * pkin(9) + t61 * t53 + t22) - t29 * mrSges(7,2) - t49 * t38;
t51 = -t70 * mrSges(6,2) + t60 * mrSges(6,3);
t52 = t70 * mrSges(6,1) - t61 * mrSges(6,3);
t113 = m(6) * t22 - t39 * mrSges(6,1) + t40 * mrSges(6,2) - t60 * t51 + t61 * t52 - t118;
t56 = t70 * mrSges(5,1) + t71 * mrSges(5,2);
t64 = -t99 * mrSges(5,2) - t70 * mrSges(5,3);
t14 = m(5) * t125 + t96 * mrSges(5,1) - t47 * mrSges(5,3) - t71 * t56 + t99 * t64 - t113;
t76 = -t88 * mrSges(4,1) + t89 * mrSges(4,2);
t81 = -qJD(3) * mrSges(4,2) + t88 * mrSges(4,3);
t137 = t105 * t33 + t139 * t36;
t23 = -t95 * pkin(4) + t96 * qJ(5) - t70 * t55 + t137;
t130 = t107 * g(1) - t110 * g(2);
t124 = qJDD(2) - t130;
t115 = (-pkin(2) * t103 - pkin(1)) * qJDD(1) + (-t135 * pkin(7) - qJ(2)) * t111 + t124;
t112 = -t79 * pkin(3) - t87 * pkin(8) + t89 * t83 + t115;
t46 = t71 * qJD(4) + t105 * t80 - t139 * t79;
t26 = (t70 * t99 - t47) * qJ(5) + (t71 * t99 + t46) * pkin(4) + t112;
t126 = -0.2e1 * qJD(5) * t61 - t100 * t23 + t102 * t26;
t17 = (t60 * t70 - t40) * pkin(9) + (t60 * t61 + t46) * pkin(5) + t126;
t131 = 0.2e1 * qJD(5) * t60 + t100 * t26 + t102 * t23;
t18 = -t59 * pkin(5) + t39 * pkin(9) - t70 * t53 + t131;
t31 = -t48 * mrSges(7,1) + t49 * mrSges(7,2);
t45 = qJDD(6) + t46;
t15 = m(7) * (-t104 * t18 + t108 * t17) - t29 * mrSges(7,3) + t45 * mrSges(7,1) - t49 * t31 + t69 * t37;
t16 = m(7) * (t104 * t17 + t108 * t18) + t28 * mrSges(7,3) - t45 * mrSges(7,2) + t48 * t31 - t69 * t38;
t50 = -t60 * mrSges(6,1) + t61 * mrSges(6,2);
t12 = m(6) * t126 + t46 * mrSges(6,1) - t40 * mrSges(6,3) + t104 * t16 + t108 * t15 - t61 * t50 + t70 * t51;
t13 = m(6) * t131 - t46 * mrSges(6,2) + t39 * mrSges(6,3) - t104 * t15 + t108 * t16 + t60 * t50 - t70 * t52;
t65 = t99 * mrSges(5,1) - t71 * mrSges(5,3);
t9 = m(5) * t137 - t96 * mrSges(5,2) - t46 * mrSges(5,3) - t100 * t12 + t102 * t13 - t70 * t56 - t99 * t65;
t6 = m(4) * t128 + qJDD(3) * mrSges(4,1) - t80 * mrSges(4,3) + qJD(3) * t81 + t105 * t9 + t139 * t14 - t89 * t76;
t82 = qJD(3) * mrSges(4,1) - t89 * mrSges(4,3);
t7 = m(4) * t136 - qJDD(3) * mrSges(4,2) + t79 * mrSges(4,3) - qJD(3) * t82 - t105 * t14 + t139 * t9 + t88 * t76;
t4 = m(3) * t129 + t106 * t7 + t109 * t6 + (-m(3) * t90 - t119) * t101;
t5 = m(3) * t127 + t119 * t103 - t106 * t6 + t109 * t7;
t140 = t101 * t5 + t103 * t4;
t117 = m(5) * t112 + t46 * mrSges(5,1) + t47 * mrSges(5,2) + t100 * t13 + t102 * t12 + t70 * t64 + t71 * t65;
t116 = -m(4) * t115 + t79 * mrSges(4,1) - t80 * mrSges(4,2) + t88 * t81 - t89 * t82 - t117;
t114 = m(3) * (-qJDD(1) * pkin(1) - t111 * qJ(2) + t124) - t116;
t8 = -t114 + (-mrSges(2,2) + t141) * t111 + (mrSges(2,1) - t122) * qJDD(1) + m(2) * t130;
t1 = m(2) * t123 - t111 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t101 * t4 + t103 * t5;
t2 = [-m(1) * g(1) + t110 * t1 - t107 * t8, t1, t5, t7, t9, t13, t16; -m(1) * g(2) + t107 * t1 + t110 * t8, t8, t4, t6, t14, t12, t15; (-m(1) - m(2)) * g(3) + t140, -m(2) * g(3) + t140, t122 * qJDD(1) - t111 * t141 + t114, -t116, t117, t113, -t118;];
f_new  = t2;
