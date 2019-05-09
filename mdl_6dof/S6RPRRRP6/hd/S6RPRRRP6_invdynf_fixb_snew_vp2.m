% Calculate vector of cutting forces with Newton-Euler
% S6RPRRRP6
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 01:38
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRRP6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP6_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP6_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP6_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP6_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP6_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:34:04
% EndTime: 2019-05-06 01:34:13
% DurationCPUTime: 3.54s
% Computational Cost: add. (42260->192), mult. (99888->240), div. (0->0), fcn. (75548->10), ass. (0->97)
t106 = qJD(3) ^ 2;
t100 = sin(qJ(3));
t104 = cos(qJ(3));
t130 = qJD(1) * qJD(2);
t96 = sin(pkin(10));
t97 = cos(pkin(10));
t126 = -t97 * g(3) - 0.2e1 * t96 * t130;
t107 = qJD(1) ^ 2;
t140 = pkin(2) * t107;
t101 = sin(qJ(1));
t105 = cos(qJ(1));
t116 = -t105 * g(1) - t101 * g(2);
t87 = -t107 * pkin(1) + qJDD(1) * qJ(2) + t116;
t62 = (-pkin(7) * qJDD(1) + t97 * t140 - t87) * t96 + t126;
t119 = -t96 * g(3) + (0.2e1 * t130 + t87) * t97;
t131 = qJDD(1) * t97;
t95 = t97 ^ 2;
t63 = pkin(7) * t131 - t95 * t140 + t119;
t120 = -t100 * t63 + t104 * t62;
t133 = t100 * t96;
t85 = (t104 * t97 - t133) * qJD(1);
t115 = t100 * t97 + t104 * t96;
t86 = t115 * qJD(1);
t71 = -t85 * pkin(3) - t86 * pkin(8);
t33 = -qJDD(3) * pkin(3) - t106 * pkin(8) + t86 * t71 - t120;
t103 = cos(qJ(4));
t132 = t85 * qJD(3);
t74 = t115 * qJDD(1) + t132;
t99 = sin(qJ(4));
t77 = t99 * qJD(3) + t103 * t86;
t49 = -t77 * qJD(4) + t103 * qJDD(3) - t99 * t74;
t83 = qJD(4) - t85;
t61 = t83 * pkin(4) - t77 * pkin(9);
t76 = t103 * qJD(3) - t99 * t86;
t75 = t76 ^ 2;
t110 = -t49 * pkin(4) - t75 * pkin(9) + t77 * t61 + t33;
t102 = cos(qJ(5));
t50 = t76 * qJD(4) + t99 * qJDD(3) + t103 * t74;
t98 = sin(qJ(5));
t53 = t102 * t77 + t98 * t76;
t29 = -t53 * qJD(5) + t102 * t49 - t98 * t50;
t52 = t102 * t76 - t98 * t77;
t30 = t52 * qJD(5) + t102 * t50 + t98 * t49;
t81 = qJD(5) + t83;
t45 = t81 * pkin(5) - t53 * qJ(6);
t46 = t81 * mrSges(7,1) - t53 * mrSges(7,3);
t51 = t52 ^ 2;
t127 = m(7) * (-t29 * pkin(5) - t51 * qJ(6) + t53 * t45 + qJDD(6) + t110) + t30 * mrSges(7,2) + t53 * t46;
t43 = -t81 * mrSges(7,2) + t52 * mrSges(7,3);
t44 = -t81 * mrSges(6,2) + t52 * mrSges(6,3);
t47 = t81 * mrSges(6,1) - t53 * mrSges(6,3);
t145 = m(6) * t110 + t30 * mrSges(6,2) + t53 * t47 + t127 - (t44 + t43) * t52 - (mrSges(6,1) + mrSges(7,1)) * t29;
t57 = -t83 * mrSges(5,2) + t76 * mrSges(5,3);
t58 = t83 * mrSges(5,1) - t77 * mrSges(5,3);
t144 = m(5) * t33 - t49 * mrSges(5,1) + t50 * mrSges(5,2) - t76 * t57 + t77 * t58 + t145;
t134 = t96 ^ 2 + t95;
t143 = t134 * mrSges(3,3);
t117 = -t97 * mrSges(3,1) + t96 * mrSges(3,2);
t113 = qJDD(1) * mrSges(3,3) + t107 * t117;
t68 = -t85 * mrSges(4,1) + t86 * mrSges(4,2);
t78 = -qJD(3) * mrSges(4,2) + t85 * mrSges(4,3);
t14 = m(4) * t120 + qJDD(3) * mrSges(4,1) - t74 * mrSges(4,3) + qJD(3) * t78 - t86 * t68 - t144;
t135 = t100 * t62 + t104 * t63;
t34 = -t106 * pkin(3) + qJDD(3) * pkin(8) + t85 * t71 + t135;
t123 = t101 * g(1) - t105 * g(2);
t118 = qJDD(2) - t123;
t109 = (-pkin(2) * t97 - pkin(1)) * qJDD(1) + (-t134 * pkin(7) - qJ(2)) * t107 + t118;
t82 = t86 * qJD(3);
t73 = -qJDD(1) * t133 + t104 * t131 - t82;
t37 = (-t74 - t132) * pkin(8) + (-t73 + t82) * pkin(3) + t109;
t121 = t103 * t37 - t99 * t34;
t70 = qJDD(4) - t73;
t21 = (t76 * t83 - t50) * pkin(9) + (t76 * t77 + t70) * pkin(4) + t121;
t137 = t103 * t34 + t99 * t37;
t23 = -t75 * pkin(4) + t49 * pkin(9) - t83 * t61 + t137;
t122 = t102 * t21 - t98 * t23;
t66 = qJDD(5) + t70;
t129 = m(7) * (-0.2e1 * qJD(6) * t53 + (t52 * t81 - t30) * qJ(6) + (t52 * t53 + t66) * pkin(5) + t122) + t81 * t43 + t66 * mrSges(7,1);
t40 = -t52 * mrSges(7,1) + t53 * mrSges(7,2);
t41 = -t52 * mrSges(6,1) + t53 * mrSges(6,2);
t12 = m(6) * t122 + t66 * mrSges(6,1) + t81 * t44 + (-t41 - t40) * t53 + (-mrSges(6,3) - mrSges(7,3)) * t30 + t129;
t138 = t102 * t23 + t98 * t21;
t128 = m(7) * (-t51 * pkin(5) + t29 * qJ(6) + 0.2e1 * qJD(6) * t52 - t81 * t45 + t138) + t29 * mrSges(7,3) + t52 * t40;
t13 = m(6) * t138 + t29 * mrSges(6,3) + t52 * t41 + (-t47 - t46) * t81 + (-mrSges(6,2) - mrSges(7,2)) * t66 + t128;
t54 = -t76 * mrSges(5,1) + t77 * mrSges(5,2);
t10 = m(5) * t121 + t70 * mrSges(5,1) - t50 * mrSges(5,3) + t102 * t12 + t98 * t13 - t77 * t54 + t83 * t57;
t11 = m(5) * t137 - t70 * mrSges(5,2) + t49 * mrSges(5,3) + t102 * t13 - t98 * t12 + t76 * t54 - t83 * t58;
t79 = qJD(3) * mrSges(4,1) - t86 * mrSges(4,3);
t7 = m(4) * t135 - qJDD(3) * mrSges(4,2) + t73 * mrSges(4,3) - qJD(3) * t79 - t99 * t10 + t103 * t11 + t85 * t68;
t4 = m(3) * t126 + t100 * t7 + t104 * t14 + (-m(3) * t87 - t113) * t96;
t5 = m(3) * t119 - t100 * t14 + t104 * t7 + t113 * t97;
t141 = t97 * t4 + t96 * t5;
t112 = -m(4) * t109 + t73 * mrSges(4,1) - t74 * mrSges(4,2) - t103 * t10 - t99 * t11 + t85 * t78 - t86 * t79;
t111 = m(3) * (-qJDD(1) * pkin(1) - t107 * qJ(2) + t118) - t112;
t6 = m(2) * t123 + (-mrSges(2,2) + t143) * t107 + (mrSges(2,1) - t117) * qJDD(1) - t111;
t1 = m(2) * t116 - t107 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t96 * t4 + t97 * t5;
t2 = [-m(1) * g(1) + t105 * t1 - t101 * t6, t1, t5, t7, t11, t13, -t66 * mrSges(7,2) - t81 * t46 + t128; -m(1) * g(2) + t101 * t1 + t105 * t6, t6, t4, t14, t10, t12, -t30 * mrSges(7,3) - t53 * t40 + t129; (-m(1) - m(2)) * g(3) + t141, -m(2) * g(3) + t141, t117 * qJDD(1) - t107 * t143 + t111, -t112, t144, t145, -t29 * mrSges(7,1) - t52 * t43 + t127;];
f_new  = t2;
