% Calculate vector of cutting forces with Newton-Euler
% S6RRRPRR4
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
% Datum: 2019-05-07 10:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPRR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 10:24:31
% EndTime: 2019-05-07 10:24:48
% DurationCPUTime: 6.44s
% Computational Cost: add. (121174->207), mult. (250292->275), div. (0->0), fcn. (184056->12), ass. (0->106)
t110 = sin(qJ(2));
t115 = cos(qJ(2));
t117 = qJD(1) ^ 2;
t108 = sin(qJ(5));
t113 = cos(qJ(5));
t107 = sin(qJ(6));
t112 = cos(qJ(6));
t105 = sin(pkin(11));
t106 = cos(pkin(11));
t103 = qJD(2) + qJD(3);
t104 = t115 ^ 2;
t111 = sin(qJ(1));
t116 = cos(qJ(1));
t131 = t111 * g(1) - t116 * g(2);
t125 = -qJDD(1) * pkin(1) - t131;
t136 = qJD(1) * t110;
t134 = qJD(1) * qJD(2);
t95 = t115 * qJDD(1) - t110 * t134;
t98 = qJD(2) * pkin(2) - pkin(8) * t136;
t122 = -t95 * pkin(2) + t98 * t136 + (-pkin(8) * t104 - pkin(7)) * t117 + t125;
t109 = sin(qJ(3));
t114 = cos(qJ(3));
t88 = (t109 * t115 + t110 * t114) * qJD(1);
t94 = t110 * qJDD(1) + t115 * t134;
t68 = t88 * qJD(3) + t109 * t94 - t114 * t95;
t135 = qJD(1) * t115;
t87 = t109 * t136 - t114 * t135;
t69 = -t87 * qJD(3) + t109 * t95 + t114 * t94;
t38 = (t103 * t87 - t69) * qJ(4) + (t103 * t88 + t68) * pkin(3) + t122;
t101 = t103 ^ 2;
t102 = qJDD(2) + qJDD(3);
t127 = -t116 * g(1) - t111 * g(2);
t90 = -t117 * pkin(1) + qJDD(1) * pkin(7) + t127;
t137 = t110 * t90;
t140 = pkin(2) * t117;
t61 = qJDD(2) * pkin(2) - t94 * pkin(8) - t137 + (pkin(8) * t134 + t110 * t140 - g(3)) * t115;
t132 = -t110 * g(3) + t115 * t90;
t62 = t95 * pkin(8) - qJD(2) * t98 - t104 * t140 + t132;
t138 = t109 * t61 + t114 * t62;
t75 = t87 * pkin(3) - t88 * qJ(4);
t41 = -t101 * pkin(3) + t102 * qJ(4) - t87 * t75 + t138;
t81 = t105 * t103 + t106 * t88;
t128 = -0.2e1 * qJD(4) * t81 - t105 * t41 + t106 * t38;
t56 = t105 * t102 + t106 * t69;
t80 = t106 * t103 - t105 * t88;
t23 = (t80 * t87 - t56) * pkin(9) + (t80 * t81 + t68) * pkin(4) + t128;
t133 = 0.2e1 * qJD(4) * t80 + t105 * t38 + t106 * t41;
t55 = t106 * t102 - t105 * t69;
t73 = t87 * pkin(4) - t81 * pkin(9);
t79 = t80 ^ 2;
t25 = -t79 * pkin(4) + t55 * pkin(9) - t87 * t73 + t133;
t130 = -t108 * t25 + t113 * t23;
t53 = -t108 * t81 + t113 * t80;
t37 = t53 * qJD(5) + t108 * t55 + t113 * t56;
t54 = t108 * t80 + t113 * t81;
t66 = qJDD(5) + t68;
t86 = qJD(5) + t87;
t17 = (t53 * t86 - t37) * pkin(10) + (t53 * t54 + t66) * pkin(5) + t130;
t139 = t108 * t23 + t113 * t25;
t36 = -t54 * qJD(5) - t108 * t56 + t113 * t55;
t50 = t86 * pkin(5) - t54 * pkin(10);
t52 = t53 ^ 2;
t18 = -t52 * pkin(5) + t36 * pkin(10) - t86 * t50 + t139;
t45 = -t107 * t54 + t112 * t53;
t28 = t45 * qJD(6) + t107 * t36 + t112 * t37;
t46 = t107 * t53 + t112 * t54;
t32 = -t45 * mrSges(7,1) + t46 * mrSges(7,2);
t85 = qJD(6) + t86;
t42 = -t85 * mrSges(7,2) + t45 * mrSges(7,3);
t63 = qJDD(6) + t66;
t14 = m(7) * (-t107 * t18 + t112 * t17) - t28 * mrSges(7,3) + t63 * mrSges(7,1) - t46 * t32 + t85 * t42;
t27 = -t46 * qJD(6) - t107 * t37 + t112 * t36;
t43 = t85 * mrSges(7,1) - t46 * mrSges(7,3);
t15 = m(7) * (t107 * t17 + t112 * t18) + t27 * mrSges(7,3) - t63 * mrSges(7,2) + t45 * t32 - t85 * t43;
t47 = -t53 * mrSges(6,1) + t54 * mrSges(6,2);
t48 = -t86 * mrSges(6,2) + t53 * mrSges(6,3);
t12 = m(6) * t130 + t66 * mrSges(6,1) - t37 * mrSges(6,3) + t107 * t15 + t112 * t14 - t54 * t47 + t86 * t48;
t49 = t86 * mrSges(6,1) - t54 * mrSges(6,3);
t13 = m(6) * t139 - t66 * mrSges(6,2) + t36 * mrSges(6,3) - t107 * t14 + t112 * t15 + t53 * t47 - t86 * t49;
t58 = -t80 * mrSges(5,1) + t81 * mrSges(5,2);
t71 = -t87 * mrSges(5,2) + t80 * mrSges(5,3);
t10 = m(5) * t128 + t68 * mrSges(5,1) - t56 * mrSges(5,3) + t108 * t13 + t113 * t12 - t81 * t58 + t87 * t71;
t72 = t87 * mrSges(5,1) - t81 * mrSges(5,3);
t11 = m(5) * t133 - t68 * mrSges(5,2) + t55 * mrSges(5,3) - t108 * t12 + t113 * t13 + t80 * t58 - t87 * t72;
t82 = -t103 * mrSges(4,2) - t87 * mrSges(4,3);
t83 = t103 * mrSges(4,1) - t88 * mrSges(4,3);
t123 = m(4) * t122 + t68 * mrSges(4,1) + t69 * mrSges(4,2) + t106 * t10 + t105 * t11 + t87 * t82 + t88 * t83;
t96 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t136;
t97 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t135;
t142 = (t110 * t96 - t115 * t97) * qJD(1) + m(3) * (-t117 * pkin(7) + t125) - t95 * mrSges(3,1) + t94 * mrSges(3,2) + t123;
t129 = -t109 * t62 + t114 * t61;
t40 = -t102 * pkin(3) - t101 * qJ(4) + t88 * t75 + qJDD(4) - t129;
t119 = -t55 * pkin(4) - t79 * pkin(9) + t81 * t73 + t40;
t124 = t27 * mrSges(7,1) + t45 * t42 - m(7) * (-t36 * pkin(5) - t52 * pkin(10) + t54 * t50 + t119) - t28 * mrSges(7,2) - t46 * t43;
t121 = -m(6) * t119 + t36 * mrSges(6,1) - t37 * mrSges(6,2) + t53 * t48 - t54 * t49 + t124;
t118 = m(5) * t40 - t55 * mrSges(5,1) + t56 * mrSges(5,2) - t80 * t71 + t81 * t72 - t121;
t76 = t87 * mrSges(4,1) + t88 * mrSges(4,2);
t16 = m(4) * t129 + t102 * mrSges(4,1) - t69 * mrSges(4,3) + t103 * t82 - t88 * t76 - t118;
t7 = m(4) * t138 - t102 * mrSges(4,2) - t68 * mrSges(4,3) - t105 * t10 - t103 * t83 + t106 * t11 - t87 * t76;
t93 = (-mrSges(3,1) * t115 + mrSges(3,2) * t110) * qJD(1);
t4 = m(3) * (-t115 * g(3) - t137) - t94 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t93 * t136 + qJD(2) * t97 + t109 * t7 + t114 * t16;
t5 = m(3) * t132 - qJDD(2) * mrSges(3,2) + t95 * mrSges(3,3) - qJD(2) * t96 - t109 * t16 + t114 * t7 + t135 * t93;
t141 = t110 * t5 + t115 * t4;
t6 = m(2) * t131 + qJDD(1) * mrSges(2,1) - t117 * mrSges(2,2) - t142;
t1 = m(2) * t127 - t117 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t110 * t4 + t115 * t5;
t2 = [-m(1) * g(1) + t116 * t1 - t111 * t6, t1, t5, t7, t11, t13, t15; -m(1) * g(2) + t111 * t1 + t116 * t6, t6, t4, t16, t10, t12, t14; (-m(1) - m(2)) * g(3) + t141, -m(2) * g(3) + t141, t142, t123, t118, -t121, -t124;];
f_new  = t2;
