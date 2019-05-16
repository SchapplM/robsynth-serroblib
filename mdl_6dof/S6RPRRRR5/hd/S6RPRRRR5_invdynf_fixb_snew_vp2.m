% Calculate vector of cutting forces with Newton-Euler
% S6RPRRRR5
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
% Datum: 2019-05-06 03:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRRR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 03:29:48
% EndTime: 2019-05-06 03:30:00
% DurationCPUTime: 5.52s
% Computational Cost: add. (92796->196), mult. (221950->253), div. (0->0), fcn. (177175->12), ass. (0->106)
t113 = qJD(1) ^ 2;
t101 = sin(pkin(11));
t102 = cos(pkin(11));
t99 = t102 ^ 2;
t136 = t101 ^ 2 + t99;
t142 = t136 * mrSges(3,3);
t106 = sin(qJ(3));
t111 = cos(qJ(3));
t124 = -t102 * mrSges(3,1) + t101 * mrSges(3,2);
t121 = qJDD(1) * mrSges(3,3) + t113 * t124;
t133 = qJD(1) * qJD(2);
t131 = -t102 * g(3) - 0.2e1 * t101 * t133;
t105 = sin(qJ(4));
t110 = cos(qJ(4));
t135 = pkin(7) * qJDD(1);
t140 = pkin(2) * t113;
t107 = sin(qJ(1));
t112 = cos(qJ(1));
t125 = -t112 * g(1) - t107 * g(2);
t91 = -t113 * pkin(1) + qJDD(1) * qJ(2) + t125;
t66 = (t102 * t140 - t135 - t91) * t101 + t131;
t127 = -t101 * g(3) + (0.2e1 * t133 + t91) * t102;
t67 = t102 * t135 - t99 * t140 + t127;
t129 = -t106 * t67 + t111 * t66;
t100 = qJD(3) + qJD(4);
t122 = -t101 * t106 + t102 * t111;
t89 = t122 * qJD(1);
t134 = t89 * qJD(3);
t123 = t101 * t111 + t102 * t106;
t81 = t123 * qJDD(1) + t134;
t90 = t123 * qJD(1);
t33 = (-t81 + t134) * pkin(8) + (t89 * t90 + qJDD(3)) * pkin(3) + t129;
t137 = t106 * t66 + t111 * t67;
t80 = -t90 * qJD(3) + t122 * qJDD(1);
t84 = qJD(3) * pkin(3) - t90 * pkin(8);
t88 = t89 ^ 2;
t38 = -t88 * pkin(3) + t80 * pkin(8) - qJD(3) * t84 + t137;
t128 = -t105 * t38 + t110 * t33;
t71 = -t105 * t90 + t110 * t89;
t72 = t105 * t89 + t110 * t90;
t57 = -t71 * pkin(4) - t72 * pkin(9);
t96 = t100 ^ 2;
t97 = qJDD(3) + qJDD(4);
t22 = -t97 * pkin(4) - t96 * pkin(9) + t72 * t57 - t128;
t103 = sin(qJ(6));
t108 = cos(qJ(6));
t104 = sin(qJ(5));
t109 = cos(qJ(5));
t48 = t71 * qJD(4) + t105 * t80 + t110 * t81;
t60 = t104 * t100 + t109 * t72;
t34 = -t60 * qJD(5) - t104 * t48 + t109 * t97;
t59 = t109 * t100 - t104 * t72;
t35 = t59 * qJD(5) + t104 * t97 + t109 * t48;
t50 = t103 * t59 + t108 * t60;
t25 = -t50 * qJD(6) - t103 * t35 + t108 * t34;
t49 = -t103 * t60 + t108 * t59;
t26 = t49 * qJD(6) + t103 * t34 + t108 * t35;
t70 = qJD(5) - t71;
t68 = qJD(6) + t70;
t39 = -t68 * mrSges(7,2) + t49 * mrSges(7,3);
t40 = t68 * mrSges(7,1) - t50 * mrSges(7,3);
t54 = t70 * pkin(5) - t60 * pkin(10);
t58 = t59 ^ 2;
t120 = t25 * mrSges(7,1) + t49 * t39 - m(7) * (-t34 * pkin(5) - t58 * pkin(10) + t60 * t54 + t22) - t26 * mrSges(7,2) - t50 * t40;
t52 = -t70 * mrSges(6,2) + t59 * mrSges(6,3);
t53 = t70 * mrSges(6,1) - t60 * mrSges(6,3);
t115 = m(6) * t22 - t34 * mrSges(6,1) + t35 * mrSges(6,2) - t59 * t52 + t60 * t53 - t120;
t56 = -t71 * mrSges(5,1) + t72 * mrSges(5,2);
t63 = -t100 * mrSges(5,2) + t71 * mrSges(5,3);
t14 = m(5) * t128 + t97 * mrSges(5,1) - t48 * mrSges(5,3) + t100 * t63 - t72 * t56 - t115;
t77 = -t89 * mrSges(4,1) + t90 * mrSges(4,2);
t82 = -qJD(3) * mrSges(4,2) + t89 * mrSges(4,3);
t138 = t105 * t33 + t110 * t38;
t23 = -t96 * pkin(4) + t97 * pkin(9) + t71 * t57 + t138;
t132 = t107 * g(1) - t112 * g(2);
t126 = qJDD(2) - t132;
t117 = (-pkin(2) * t102 - pkin(1)) * qJDD(1) + (-t136 * pkin(7) - qJ(2)) * t113 + t126;
t114 = -t80 * pkin(3) - t88 * pkin(8) + t90 * t84 + t117;
t47 = -t72 * qJD(4) - t105 * t81 + t110 * t80;
t29 = (-t100 * t71 - t48) * pkin(9) + (t100 * t72 - t47) * pkin(4) + t114;
t130 = -t104 * t23 + t109 * t29;
t46 = qJDD(5) - t47;
t17 = (t59 * t70 - t35) * pkin(10) + (t59 * t60 + t46) * pkin(5) + t130;
t139 = t104 * t29 + t109 * t23;
t18 = -t58 * pkin(5) + t34 * pkin(10) - t70 * t54 + t139;
t31 = -t49 * mrSges(7,1) + t50 * mrSges(7,2);
t44 = qJDD(6) + t46;
t15 = m(7) * (-t103 * t18 + t108 * t17) - t26 * mrSges(7,3) + t44 * mrSges(7,1) - t50 * t31 + t68 * t39;
t16 = m(7) * (t103 * t17 + t108 * t18) + t25 * mrSges(7,3) - t44 * mrSges(7,2) + t49 * t31 - t68 * t40;
t51 = -t59 * mrSges(6,1) + t60 * mrSges(6,2);
t12 = m(6) * t130 + t46 * mrSges(6,1) - t35 * mrSges(6,3) + t103 * t16 + t108 * t15 - t60 * t51 + t70 * t52;
t13 = m(6) * t139 - t46 * mrSges(6,2) + t34 * mrSges(6,3) - t103 * t15 + t108 * t16 + t59 * t51 - t70 * t53;
t64 = t100 * mrSges(5,1) - t72 * mrSges(5,3);
t9 = m(5) * t138 - t97 * mrSges(5,2) + t47 * mrSges(5,3) - t100 * t64 - t104 * t12 + t109 * t13 + t71 * t56;
t6 = m(4) * t129 + qJDD(3) * mrSges(4,1) - t81 * mrSges(4,3) + qJD(3) * t82 + t105 * t9 + t110 * t14 - t90 * t77;
t83 = qJD(3) * mrSges(4,1) - t90 * mrSges(4,3);
t7 = m(4) * t137 - qJDD(3) * mrSges(4,2) + t80 * mrSges(4,3) - qJD(3) * t83 - t105 * t14 + t110 * t9 + t89 * t77;
t4 = m(3) * t131 + t106 * t7 + t111 * t6 + (-m(3) * t91 - t121) * t101;
t5 = m(3) * t127 + t121 * t102 - t106 * t6 + t111 * t7;
t141 = t101 * t5 + t102 * t4;
t119 = -m(5) * t114 + t47 * mrSges(5,1) - t48 * mrSges(5,2) - t104 * t13 - t109 * t12 + t71 * t63 - t72 * t64;
t118 = -m(4) * t117 + t80 * mrSges(4,1) - t81 * mrSges(4,2) + t89 * t82 - t90 * t83 + t119;
t116 = m(3) * (-qJDD(1) * pkin(1) - t113 * qJ(2) + t126) - t118;
t8 = (-mrSges(2,2) + t142) * t113 + (mrSges(2,1) - t124) * qJDD(1) - t116 + m(2) * t132;
t1 = m(2) * t125 - t113 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t101 * t4 + t102 * t5;
t2 = [-m(1) * g(1) + t112 * t1 - t107 * t8, t1, t5, t7, t9, t13, t16; -m(1) * g(2) + t107 * t1 + t112 * t8, t8, t4, t6, t14, t12, t15; (-m(1) - m(2)) * g(3) + t141, -m(2) * g(3) + t141, t124 * qJDD(1) - t113 * t142 + t116, -t118, -t119, t115, -t120;];
f_new  = t2;
