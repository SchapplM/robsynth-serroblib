% Calculate vector of cutting forces with Newton-Euler
% S6RPRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-05-05 19:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPRR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR6_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR6_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR6_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR6_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR6_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 19:01:19
% EndTime: 2019-05-05 19:01:31
% DurationCPUTime: 5.44s
% Computational Cost: add. (89921->195), mult. (222189->256), div. (0->0), fcn. (170494->12), ass. (0->105)
t112 = qJD(1) ^ 2;
t100 = sin(pkin(10));
t102 = cos(pkin(10));
t98 = t102 ^ 2;
t137 = t100 ^ 2 + t98;
t142 = t137 * mrSges(3,3);
t105 = sin(qJ(3));
t109 = cos(qJ(3));
t122 = -t102 * mrSges(3,1) + t100 * mrSges(3,2);
t120 = qJDD(1) * mrSges(3,3) + t112 * t122;
t133 = qJD(1) * qJD(2);
t129 = -t102 * g(3) - 0.2e1 * t100 * t133;
t111 = qJD(3) ^ 2;
t140 = pkin(2) * t112;
t106 = sin(qJ(1));
t110 = cos(qJ(1));
t123 = -t110 * g(1) - t106 * g(2);
t90 = -t112 * pkin(1) + qJDD(1) * qJ(2) + t123;
t60 = (-pkin(7) * qJDD(1) + t102 * t140 - t90) * t100 + t129;
t126 = -t100 * g(3) + (0.2e1 * t133 + t90) * t102;
t132 = qJDD(1) * t102;
t66 = pkin(7) * t132 - t98 * t140 + t126;
t127 = -t105 * t66 + t109 * t60;
t136 = t100 * t105;
t88 = (-t102 * t109 + t136) * qJD(1);
t121 = t100 * t109 + t102 * t105;
t89 = t121 * qJD(1);
t70 = t88 * pkin(3) - t89 * qJ(4);
t35 = -qJDD(3) * pkin(3) - t111 * qJ(4) + t89 * t70 + qJDD(4) - t127;
t101 = cos(pkin(11));
t99 = sin(pkin(11));
t81 = t99 * qJD(3) + t101 * t89;
t63 = t88 * pkin(4) - t81 * pkin(8);
t135 = t88 * qJD(3);
t76 = t121 * qJDD(1) - t135;
t64 = t101 * qJDD(3) - t99 * t76;
t80 = t101 * qJD(3) - t99 * t89;
t79 = t80 ^ 2;
t115 = -t64 * pkin(4) - t79 * pkin(8) + t81 * t63 + t35;
t103 = sin(qJ(6));
t107 = cos(qJ(6));
t104 = sin(qJ(5));
t108 = cos(qJ(5));
t54 = t104 * t80 + t108 * t81;
t65 = t99 * qJDD(3) + t101 * t76;
t36 = -t54 * qJD(5) - t104 * t65 + t108 * t64;
t53 = -t104 * t81 + t108 * t80;
t37 = t53 * qJD(5) + t104 * t64 + t108 * t65;
t46 = t103 * t53 + t107 * t54;
t27 = -t46 * qJD(6) - t103 * t37 + t107 * t36;
t45 = -t103 * t54 + t107 * t53;
t28 = t45 * qJD(6) + t103 * t36 + t107 * t37;
t86 = qJD(5) + t88;
t85 = qJD(6) + t86;
t41 = -t85 * mrSges(7,2) + t45 * mrSges(7,3);
t42 = t85 * mrSges(7,1) - t46 * mrSges(7,3);
t50 = t86 * pkin(5) - t54 * pkin(9);
t52 = t53 ^ 2;
t119 = t27 * mrSges(7,1) + t45 * t41 - m(7) * (-t36 * pkin(5) - t52 * pkin(9) + t54 * t50 + t115) - t28 * mrSges(7,2) - t46 * t42;
t48 = -t86 * mrSges(6,2) + t53 * mrSges(6,3);
t49 = t86 * mrSges(6,1) - t54 * mrSges(6,3);
t116 = -m(6) * t115 + t36 * mrSges(6,1) - t37 * mrSges(6,2) + t53 * t48 - t54 * t49 + t119;
t61 = -t88 * mrSges(5,2) + t80 * mrSges(5,3);
t62 = t88 * mrSges(5,1) - t81 * mrSges(5,3);
t113 = m(5) * t35 - t64 * mrSges(5,1) + t65 * mrSges(5,2) - t80 * t61 + t81 * t62 - t116;
t71 = t88 * mrSges(4,1) + t89 * mrSges(4,2);
t82 = -qJD(3) * mrSges(4,2) - t88 * mrSges(4,3);
t16 = m(4) * t127 + qJDD(3) * mrSges(4,1) - t76 * mrSges(4,3) + qJD(3) * t82 - t89 * t71 - t113;
t138 = t105 * t60 + t109 * t66;
t38 = -t111 * pkin(3) + qJDD(3) * qJ(4) - t88 * t70 + t138;
t130 = t106 * g(1) - t110 * g(2);
t124 = qJDD(2) - t130;
t114 = (-pkin(2) * t102 - pkin(1)) * qJDD(1) + (-t137 * pkin(7) - qJ(2)) * t112 + t124;
t134 = t89 * qJD(3);
t75 = qJDD(1) * t136 - t109 * t132 + t134;
t43 = (-t76 + t135) * qJ(4) + (t75 + t134) * pkin(3) + t114;
t125 = -0.2e1 * qJD(4) * t81 + t101 * t43 - t99 * t38;
t23 = (t80 * t88 - t65) * pkin(8) + (t80 * t81 + t75) * pkin(4) + t125;
t131 = 0.2e1 * qJD(4) * t80 + t101 * t38 + t99 * t43;
t25 = -t79 * pkin(4) + t64 * pkin(8) - t88 * t63 + t131;
t128 = -t104 * t25 + t108 * t23;
t73 = qJDD(5) + t75;
t17 = (t53 * t86 - t37) * pkin(9) + (t53 * t54 + t73) * pkin(5) + t128;
t139 = t104 * t23 + t108 * t25;
t18 = -t52 * pkin(5) + t36 * pkin(9) - t86 * t50 + t139;
t32 = -t45 * mrSges(7,1) + t46 * mrSges(7,2);
t68 = qJDD(6) + t73;
t14 = m(7) * (-t103 * t18 + t107 * t17) - t28 * mrSges(7,3) + t68 * mrSges(7,1) - t46 * t32 + t85 * t41;
t15 = m(7) * (t103 * t17 + t107 * t18) + t27 * mrSges(7,3) - t68 * mrSges(7,2) + t45 * t32 - t85 * t42;
t47 = -t53 * mrSges(6,1) + t54 * mrSges(6,2);
t12 = m(6) * t128 + t73 * mrSges(6,1) - t37 * mrSges(6,3) + t103 * t15 + t107 * t14 - t54 * t47 + t86 * t48;
t13 = m(6) * t139 - t73 * mrSges(6,2) + t36 * mrSges(6,3) - t103 * t14 + t107 * t15 + t53 * t47 - t86 * t49;
t55 = -t80 * mrSges(5,1) + t81 * mrSges(5,2);
t10 = m(5) * t125 + t75 * mrSges(5,1) - t65 * mrSges(5,3) + t104 * t13 + t108 * t12 - t81 * t55 + t88 * t61;
t11 = m(5) * t131 - t75 * mrSges(5,2) + t64 * mrSges(5,3) - t104 * t12 + t108 * t13 + t80 * t55 - t88 * t62;
t83 = qJD(3) * mrSges(4,1) - t89 * mrSges(4,3);
t7 = m(4) * t138 - qJDD(3) * mrSges(4,2) - t75 * mrSges(4,3) - qJD(3) * t83 - t99 * t10 + t101 * t11 - t88 * t71;
t4 = m(3) * t129 + t105 * t7 + t109 * t16 + (-m(3) * t90 - t120) * t100;
t5 = m(3) * t126 + t120 * t102 - t105 * t16 + t109 * t7;
t141 = t100 * t5 + t102 * t4;
t118 = m(4) * t114 + t75 * mrSges(4,1) + t76 * mrSges(4,2) + t101 * t10 + t99 * t11 + t88 * t82 + t89 * t83;
t117 = m(3) * (-qJDD(1) * pkin(1) - t112 * qJ(2) + t124) + t118;
t6 = m(2) * t130 + (-mrSges(2,2) + t142) * t112 + (mrSges(2,1) - t122) * qJDD(1) - t117;
t1 = m(2) * t123 - t112 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t100 * t4 + t102 * t5;
t2 = [-m(1) * g(1) + t110 * t1 - t106 * t6, t1, t5, t7, t11, t13, t15; -m(1) * g(2) + t106 * t1 + t110 * t6, t6, t4, t16, t10, t12, t14; (-m(1) - m(2)) * g(3) + t141, -m(2) * g(3) + t141, t122 * qJDD(1) - t112 * t142 + t117, t118, t113, -t116, -t119;];
f_new  = t2;
