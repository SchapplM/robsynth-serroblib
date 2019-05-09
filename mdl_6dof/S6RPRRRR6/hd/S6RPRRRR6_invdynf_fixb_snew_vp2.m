% Calculate vector of cutting forces with Newton-Euler
% S6RPRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-05-06 03:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRRR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR6_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR6_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR6_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR6_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR6_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 03:46:28
% EndTime: 2019-05-06 03:46:40
% DurationCPUTime: 5.72s
% Computational Cost: add. (98140->197), mult. (233256->254), div. (0->0), fcn. (180688->12), ass. (0->107)
t113 = qJD(1) ^ 2;
t100 = sin(pkin(11));
t101 = cos(pkin(11));
t99 = t101 ^ 2;
t136 = t100 ^ 2 + t99;
t142 = t136 * mrSges(3,3);
t105 = sin(qJ(3));
t110 = cos(qJ(3));
t123 = -t101 * mrSges(3,1) + t100 * mrSges(3,2);
t121 = qJDD(1) * mrSges(3,3) + t113 * t123;
t133 = qJD(1) * qJD(2);
t130 = -t101 * g(3) - 0.2e1 * t100 * t133;
t112 = qJD(3) ^ 2;
t140 = pkin(2) * t113;
t106 = sin(qJ(1));
t111 = cos(qJ(1));
t124 = -t111 * g(1) - t106 * g(2);
t91 = -t113 * pkin(1) + qJDD(1) * qJ(2) + t124;
t65 = (-pkin(7) * qJDD(1) + t101 * t140 - t91) * t100 + t130;
t126 = -t100 * g(3) + (0.2e1 * t133 + t91) * t101;
t132 = qJDD(1) * t101;
t66 = pkin(7) * t132 - t99 * t140 + t126;
t127 = -t105 * t66 + t110 * t65;
t135 = t100 * t105;
t89 = (t101 * t110 - t135) * qJD(1);
t122 = t100 * t110 + t101 * t105;
t90 = t122 * qJD(1);
t74 = -t89 * pkin(3) - t90 * pkin(8);
t37 = -qJDD(3) * pkin(3) - t112 * pkin(8) + t90 * t74 - t127;
t104 = sin(qJ(4));
t109 = cos(qJ(4));
t134 = t89 * qJD(3);
t77 = t122 * qJDD(1) + t134;
t80 = t104 * qJD(3) + t109 * t90;
t52 = -t80 * qJD(4) + t109 * qJDD(3) - t104 * t77;
t87 = qJD(4) - t89;
t64 = t87 * pkin(4) - t80 * pkin(9);
t79 = t109 * qJD(3) - t104 * t90;
t78 = t79 ^ 2;
t117 = -t52 * pkin(4) - t78 * pkin(9) + t80 * t64 + t37;
t102 = sin(qJ(6));
t107 = cos(qJ(6));
t103 = sin(qJ(5));
t108 = cos(qJ(5));
t53 = t79 * qJD(4) + t104 * qJDD(3) + t109 * t77;
t56 = t103 * t79 + t108 * t80;
t34 = -t56 * qJD(5) - t103 * t53 + t108 * t52;
t55 = -t103 * t80 + t108 * t79;
t35 = t55 * qJD(5) + t103 * t52 + t108 * t53;
t46 = t102 * t55 + t107 * t56;
t22 = -t46 * qJD(6) - t102 * t35 + t107 * t34;
t45 = -t102 * t56 + t107 * t55;
t23 = t45 * qJD(6) + t102 * t34 + t107 * t35;
t85 = qJD(5) + t87;
t83 = qJD(6) + t85;
t41 = -t83 * mrSges(7,2) + t45 * mrSges(7,3);
t42 = t83 * mrSges(7,1) - t46 * mrSges(7,3);
t50 = t85 * pkin(5) - t56 * pkin(10);
t54 = t55 ^ 2;
t120 = t22 * mrSges(7,1) + t45 * t41 - m(7) * (-t34 * pkin(5) - t54 * pkin(10) + t56 * t50 + t117) - t23 * mrSges(7,2) - t46 * t42;
t48 = -t85 * mrSges(6,2) + t55 * mrSges(6,3);
t49 = t85 * mrSges(6,1) - t56 * mrSges(6,3);
t116 = -m(6) * t117 + t34 * mrSges(6,1) - t35 * mrSges(6,2) + t55 * t48 - t56 * t49 + t120;
t60 = -t87 * mrSges(5,2) + t79 * mrSges(5,3);
t61 = t87 * mrSges(5,1) - t80 * mrSges(5,3);
t114 = m(5) * t37 - t52 * mrSges(5,1) + t53 * mrSges(5,2) - t79 * t60 + t80 * t61 - t116;
t71 = -t89 * mrSges(4,1) + t90 * mrSges(4,2);
t81 = -qJD(3) * mrSges(4,2) + t89 * mrSges(4,3);
t14 = m(4) * t127 + qJDD(3) * mrSges(4,1) - t77 * mrSges(4,3) + qJD(3) * t81 - t90 * t71 - t114;
t137 = t105 * t65 + t110 * t66;
t38 = -t112 * pkin(3) + qJDD(3) * pkin(8) + t89 * t74 + t137;
t131 = t106 * g(1) - t111 * g(2);
t125 = qJDD(2) - t131;
t115 = (-pkin(2) * t101 - pkin(1)) * qJDD(1) + (-t136 * pkin(7) - qJ(2)) * t113 + t125;
t86 = t90 * qJD(3);
t76 = -qJDD(1) * t135 + t110 * t132 - t86;
t43 = (-t77 - t134) * pkin(8) + (-t76 + t86) * pkin(3) + t115;
t128 = -t104 * t38 + t109 * t43;
t73 = qJDD(4) - t76;
t26 = (t79 * t87 - t53) * pkin(9) + (t79 * t80 + t73) * pkin(4) + t128;
t138 = t104 * t43 + t109 * t38;
t28 = -t78 * pkin(4) + t52 * pkin(9) - t87 * t64 + t138;
t129 = -t103 * t28 + t108 * t26;
t69 = qJDD(5) + t73;
t17 = (t55 * t85 - t35) * pkin(10) + (t55 * t56 + t69) * pkin(5) + t129;
t139 = t103 * t26 + t108 * t28;
t18 = -t54 * pkin(5) + t34 * pkin(10) - t85 * t50 + t139;
t32 = -t45 * mrSges(7,1) + t46 * mrSges(7,2);
t67 = qJDD(6) + t69;
t15 = m(7) * (-t102 * t18 + t107 * t17) - t23 * mrSges(7,3) + t67 * mrSges(7,1) - t46 * t32 + t83 * t41;
t16 = m(7) * (t102 * t17 + t107 * t18) + t22 * mrSges(7,3) - t67 * mrSges(7,2) + t45 * t32 - t83 * t42;
t47 = -t55 * mrSges(6,1) + t56 * mrSges(6,2);
t12 = m(6) * t129 + t69 * mrSges(6,1) - t35 * mrSges(6,3) + t102 * t16 + t107 * t15 - t56 * t47 + t85 * t48;
t13 = m(6) * t139 - t69 * mrSges(6,2) + t34 * mrSges(6,3) - t102 * t15 + t107 * t16 + t55 * t47 - t85 * t49;
t57 = -t79 * mrSges(5,1) + t80 * mrSges(5,2);
t10 = m(5) * t128 + t73 * mrSges(5,1) - t53 * mrSges(5,3) + t103 * t13 + t108 * t12 - t80 * t57 + t87 * t60;
t11 = m(5) * t138 - t73 * mrSges(5,2) + t52 * mrSges(5,3) - t103 * t12 + t108 * t13 + t79 * t57 - t87 * t61;
t82 = qJD(3) * mrSges(4,1) - t90 * mrSges(4,3);
t7 = m(4) * t137 - qJDD(3) * mrSges(4,2) + t76 * mrSges(4,3) - qJD(3) * t82 - t104 * t10 + t109 * t11 + t89 * t71;
t4 = m(3) * t130 + t105 * t7 + t110 * t14 + (-m(3) * t91 - t121) * t100;
t5 = m(3) * t126 + t121 * t101 - t105 * t14 + t110 * t7;
t141 = t100 * t5 + t101 * t4;
t119 = -m(4) * t115 + t76 * mrSges(4,1) - t77 * mrSges(4,2) - t109 * t10 - t104 * t11 + t89 * t81 - t90 * t82;
t118 = m(3) * (-qJDD(1) * pkin(1) - t113 * qJ(2) + t125) - t119;
t6 = m(2) * t131 + (-mrSges(2,2) + t142) * t113 + (mrSges(2,1) - t123) * qJDD(1) - t118;
t1 = m(2) * t124 - t113 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t100 * t4 + t101 * t5;
t2 = [-m(1) * g(1) + t111 * t1 - t106 * t6, t1, t5, t7, t11, t13, t16; -m(1) * g(2) + t106 * t1 + t111 * t6, t6, t4, t14, t10, t12, t15; (-m(1) - m(2)) * g(3) + t141, -m(2) * g(3) + t141, t123 * qJDD(1) - t113 * t142 + t118, -t119, t114, -t116, -t120;];
f_new  = t2;
