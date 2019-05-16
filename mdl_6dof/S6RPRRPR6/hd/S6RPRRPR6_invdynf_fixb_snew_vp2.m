% Calculate vector of cutting forces with Newton-Euler
% S6RPRRPR6
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
% Datum: 2019-05-05 22:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRPR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR6_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR6_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR6_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR6_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR6_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:48:08
% EndTime: 2019-05-05 22:48:20
% DurationCPUTime: 5.54s
% Computational Cost: add. (94197->196), mult. (227260->256), div. (0->0), fcn. (174724->12), ass. (0->105)
t113 = qJD(1) ^ 2;
t101 = sin(pkin(10));
t103 = cos(pkin(10));
t99 = t103 ^ 2;
t137 = t101 ^ 2 + t99;
t142 = t137 * mrSges(3,3);
t106 = sin(qJ(3));
t110 = cos(qJ(3));
t123 = -t103 * mrSges(3,1) + t101 * mrSges(3,2);
t121 = qJDD(1) * mrSges(3,3) + t113 * t123;
t134 = qJD(1) * qJD(2);
t130 = -t103 * g(3) - 0.2e1 * t101 * t134;
t112 = qJD(3) ^ 2;
t140 = pkin(2) * t113;
t107 = sin(qJ(1));
t111 = cos(qJ(1));
t124 = -t111 * g(1) - t107 * g(2);
t91 = -t113 * pkin(1) + qJDD(1) * qJ(2) + t124;
t67 = (-pkin(7) * qJDD(1) + t103 * t140 - t91) * t101 + t130;
t127 = -t101 * g(3) + (0.2e1 * t134 + t91) * t103;
t133 = qJDD(1) * t103;
t68 = pkin(7) * t133 - t99 * t140 + t127;
t128 = -t106 * t68 + t110 * t67;
t136 = t101 * t106;
t89 = (t103 * t110 - t136) * qJD(1);
t122 = t101 * t110 + t103 * t106;
t90 = t122 * qJD(1);
t75 = -t89 * pkin(3) - t90 * pkin(8);
t34 = -qJDD(3) * pkin(3) - t112 * pkin(8) + t90 * t75 - t128;
t105 = sin(qJ(4));
t109 = cos(qJ(4));
t135 = t89 * qJD(3);
t78 = t122 * qJDD(1) + t135;
t81 = t105 * qJD(3) + t109 * t90;
t52 = -t81 * qJD(4) + t109 * qJDD(3) - t105 * t78;
t87 = qJD(4) - t89;
t63 = t87 * pkin(4) - t81 * qJ(5);
t80 = t109 * qJD(3) - t105 * t90;
t79 = t80 ^ 2;
t116 = -t52 * pkin(4) - t79 * qJ(5) + t81 * t63 + qJDD(5) + t34;
t104 = sin(qJ(6));
t108 = cos(qJ(6));
t100 = sin(pkin(11));
t102 = cos(pkin(11));
t53 = t80 * qJD(4) + t105 * qJDD(3) + t109 * t78;
t37 = -t100 * t53 + t102 * t52;
t38 = t100 * t52 + t102 * t53;
t57 = -t100 * t81 + t102 * t80;
t58 = t100 * t80 + t102 * t81;
t46 = t104 * t57 + t108 * t58;
t27 = -t46 * qJD(6) - t104 * t38 + t108 * t37;
t45 = -t104 * t58 + t108 * t57;
t28 = t45 * qJD(6) + t104 * t37 + t108 * t38;
t85 = qJD(6) + t87;
t41 = -t85 * mrSges(7,2) + t45 * mrSges(7,3);
t42 = t85 * mrSges(7,1) - t46 * mrSges(7,3);
t50 = t87 * pkin(5) - t58 * pkin(9);
t56 = t57 ^ 2;
t120 = t27 * mrSges(7,1) + t45 * t41 - m(7) * (-t37 * pkin(5) - t56 * pkin(9) + t58 * t50 + t116) - t28 * mrSges(7,2) - t46 * t42;
t48 = -t87 * mrSges(6,2) + t57 * mrSges(6,3);
t49 = t87 * mrSges(6,1) - t58 * mrSges(6,3);
t117 = -m(6) * t116 + t37 * mrSges(6,1) - t38 * mrSges(6,2) + t57 * t48 - t58 * t49 + t120;
t62 = -t87 * mrSges(5,2) + t80 * mrSges(5,3);
t64 = t87 * mrSges(5,1) - t81 * mrSges(5,3);
t114 = m(5) * t34 - t52 * mrSges(5,1) + t53 * mrSges(5,2) - t80 * t62 + t81 * t64 - t117;
t72 = -t89 * mrSges(4,1) + t90 * mrSges(4,2);
t82 = -qJD(3) * mrSges(4,2) + t89 * mrSges(4,3);
t16 = m(4) * t128 + qJDD(3) * mrSges(4,1) - t78 * mrSges(4,3) + qJD(3) * t82 - t90 * t72 - t114;
t138 = t106 * t67 + t110 * t68;
t35 = -t112 * pkin(3) + qJDD(3) * pkin(8) + t89 * t75 + t138;
t131 = t107 * g(1) - t111 * g(2);
t125 = qJDD(2) - t131;
t115 = (-pkin(2) * t103 - pkin(1)) * qJDD(1) + (-t137 * pkin(7) - qJ(2)) * t113 + t125;
t86 = t90 * qJD(3);
t77 = -qJDD(1) * t136 + t110 * t133 - t86;
t43 = (-t78 - t135) * pkin(8) + (-t77 + t86) * pkin(3) + t115;
t129 = -t105 * t35 + t109 * t43;
t74 = qJDD(4) - t77;
t23 = (t80 * t87 - t53) * qJ(5) + (t80 * t81 + t74) * pkin(4) + t129;
t139 = t105 * t43 + t109 * t35;
t25 = -t79 * pkin(4) + t52 * qJ(5) - t87 * t63 + t139;
t126 = -0.2e1 * qJD(5) * t58 - t100 * t25 + t102 * t23;
t17 = (t57 * t87 - t38) * pkin(9) + (t57 * t58 + t74) * pkin(5) + t126;
t132 = 0.2e1 * qJD(5) * t57 + t100 * t23 + t102 * t25;
t18 = -t56 * pkin(5) + t37 * pkin(9) - t87 * t50 + t132;
t32 = -t45 * mrSges(7,1) + t46 * mrSges(7,2);
t70 = qJDD(6) + t74;
t14 = m(7) * (-t104 * t18 + t108 * t17) - t28 * mrSges(7,3) + t70 * mrSges(7,1) - t46 * t32 + t85 * t41;
t15 = m(7) * (t104 * t17 + t108 * t18) + t27 * mrSges(7,3) - t70 * mrSges(7,2) + t45 * t32 - t85 * t42;
t47 = -t57 * mrSges(6,1) + t58 * mrSges(6,2);
t12 = m(6) * t126 + t74 * mrSges(6,1) - t38 * mrSges(6,3) + t104 * t15 + t108 * t14 - t58 * t47 + t87 * t48;
t13 = m(6) * t132 - t74 * mrSges(6,2) + t37 * mrSges(6,3) - t104 * t14 + t108 * t15 + t57 * t47 - t87 * t49;
t59 = -t80 * mrSges(5,1) + t81 * mrSges(5,2);
t10 = m(5) * t129 + t74 * mrSges(5,1) - t53 * mrSges(5,3) + t100 * t13 + t102 * t12 - t81 * t59 + t87 * t62;
t11 = m(5) * t139 - t74 * mrSges(5,2) + t52 * mrSges(5,3) - t100 * t12 + t102 * t13 + t80 * t59 - t87 * t64;
t83 = qJD(3) * mrSges(4,1) - t90 * mrSges(4,3);
t7 = m(4) * t138 - qJDD(3) * mrSges(4,2) + t77 * mrSges(4,3) - qJD(3) * t83 - t105 * t10 + t109 * t11 + t89 * t72;
t4 = m(3) * t130 + t106 * t7 + t110 * t16 + (-m(3) * t91 - t121) * t101;
t5 = m(3) * t127 + t121 * t103 - t106 * t16 + t110 * t7;
t141 = t101 * t5 + t103 * t4;
t119 = -m(4) * t115 + t77 * mrSges(4,1) - t78 * mrSges(4,2) - t109 * t10 - t105 * t11 + t89 * t82 - t90 * t83;
t118 = m(3) * (-qJDD(1) * pkin(1) - t113 * qJ(2) + t125) - t119;
t6 = m(2) * t131 + (-mrSges(2,2) + t142) * t113 + (mrSges(2,1) - t123) * qJDD(1) - t118;
t1 = m(2) * t124 - t113 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t101 * t4 + t103 * t5;
t2 = [-m(1) * g(1) + t111 * t1 - t107 * t6, t1, t5, t7, t11, t13, t15; -m(1) * g(2) + t107 * t1 + t111 * t6, t6, t4, t16, t10, t12, t14; (-m(1) - m(2)) * g(3) + t141, -m(2) * g(3) + t141, t123 * qJDD(1) - t113 * t142 + t118, -t119, t114, -t117, -t120;];
f_new  = t2;
