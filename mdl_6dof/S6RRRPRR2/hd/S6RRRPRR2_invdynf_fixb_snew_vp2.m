% Calculate vector of cutting forces with Newton-Euler
% S6RRRPRR2
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
% Datum: 2019-05-07 10:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPRR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:56:23
% EndTime: 2019-05-07 09:56:36
% DurationCPUTime: 5.93s
% Computational Cost: add. (103963->208), mult. (233590->276), div. (0->0), fcn. (175413->12), ass. (0->106)
t105 = sin(pkin(11));
t106 = cos(pkin(11));
t102 = qJDD(2) + qJDD(3);
t103 = qJD(2) + qJD(3);
t109 = sin(qJ(3));
t114 = cos(qJ(3));
t110 = sin(qJ(2));
t115 = cos(qJ(2));
t134 = qJD(1) * qJD(2);
t117 = qJD(1) ^ 2;
t111 = sin(qJ(1));
t116 = cos(qJ(1));
t127 = -t116 * g(1) - t111 * g(2);
t92 = -t117 * pkin(1) + qJDD(1) * pkin(7) + t127;
t138 = t110 * t92;
t141 = pkin(2) * t117;
t95 = t110 * qJDD(1) + t115 * t134;
t61 = qJDD(2) * pkin(2) - t95 * pkin(8) - t138 + (pkin(8) * t134 + t110 * t141 - g(3)) * t115;
t104 = t115 ^ 2;
t132 = -t110 * g(3) + t115 * t92;
t96 = t115 * qJDD(1) - t110 * t134;
t136 = qJD(1) * t110;
t99 = qJD(2) * pkin(2) - pkin(8) * t136;
t62 = t96 * pkin(8) - qJD(2) * t99 - t104 * t141 + t132;
t128 = -t109 * t62 + t114 * t61;
t89 = (-t109 * t110 + t114 * t115) * qJD(1);
t70 = t89 * qJD(3) + t109 * t96 + t114 * t95;
t90 = (t109 * t115 + t110 * t114) * qJD(1);
t32 = (t103 * t89 - t70) * qJ(4) + (t89 * t90 + t102) * pkin(3) + t128;
t139 = t109 * t61 + t114 * t62;
t69 = -t90 * qJD(3) - t109 * t95 + t114 * t96;
t84 = t103 * pkin(3) - t90 * qJ(4);
t88 = t89 ^ 2;
t36 = -t88 * pkin(3) + t69 * qJ(4) - t103 * t84 + t139;
t81 = t105 * t89 + t106 * t90;
t144 = -0.2e1 * qJD(4) * t81 - t105 * t36 + t106 * t32;
t131 = t111 * g(1) - t116 * g(2);
t125 = -qJDD(1) * pkin(1) - t131;
t122 = -t96 * pkin(2) + t99 * t136 + (-pkin(8) * t104 - pkin(7)) * t117 + t125;
t108 = sin(qJ(5));
t113 = cos(qJ(5));
t119 = -t69 * pkin(3) - t88 * qJ(4) + t90 * t84 + qJDD(4) + t122;
t107 = sin(qJ(6));
t112 = cos(qJ(6));
t101 = t103 ^ 2;
t80 = -t105 * t90 + t106 * t89;
t133 = 0.2e1 * qJD(4) * t80 + t105 * t32 + t106 * t36;
t57 = -t80 * pkin(4) - t81 * pkin(9);
t23 = -t101 * pkin(4) + t102 * pkin(9) + t80 * t57 + t133;
t49 = -t105 * t70 + t106 * t69;
t50 = t105 * t69 + t106 * t70;
t26 = (-t103 * t80 - t50) * pkin(9) + (t103 * t81 - t49) * pkin(4) + t119;
t129 = -t108 * t23 + t113 * t26;
t66 = t113 * t103 - t108 * t81;
t38 = t66 * qJD(5) + t108 * t102 + t113 * t50;
t48 = qJDD(5) - t49;
t67 = t108 * t103 + t113 * t81;
t77 = qJD(5) - t80;
t17 = (t66 * t77 - t38) * pkin(10) + (t66 * t67 + t48) * pkin(5) + t129;
t140 = t108 * t26 + t113 * t23;
t37 = -t67 * qJD(5) + t113 * t102 - t108 * t50;
t54 = t77 * pkin(5) - t67 * pkin(10);
t64 = t66 ^ 2;
t18 = -t64 * pkin(5) + t37 * pkin(10) - t77 * t54 + t140;
t45 = -t107 * t67 + t112 * t66;
t29 = t45 * qJD(6) + t107 * t37 + t112 * t38;
t46 = t107 * t66 + t112 * t67;
t33 = -t45 * mrSges(7,1) + t46 * mrSges(7,2);
t74 = qJD(6) + t77;
t41 = -t74 * mrSges(7,2) + t45 * mrSges(7,3);
t44 = qJDD(6) + t48;
t15 = m(7) * (-t107 * t18 + t112 * t17) - t29 * mrSges(7,3) + t44 * mrSges(7,1) - t46 * t33 + t74 * t41;
t28 = -t46 * qJD(6) - t107 * t38 + t112 * t37;
t42 = t74 * mrSges(7,1) - t46 * mrSges(7,3);
t16 = m(7) * (t107 * t17 + t112 * t18) + t28 * mrSges(7,3) - t44 * mrSges(7,2) + t45 * t33 - t74 * t42;
t51 = -t66 * mrSges(6,1) + t67 * mrSges(6,2);
t52 = -t77 * mrSges(6,2) + t66 * mrSges(6,3);
t12 = m(6) * t129 + t48 * mrSges(6,1) - t38 * mrSges(6,3) + t107 * t16 + t112 * t15 - t67 * t51 + t77 * t52;
t53 = t77 * mrSges(6,1) - t67 * mrSges(6,3);
t13 = m(6) * t140 - t48 * mrSges(6,2) + t37 * mrSges(6,3) - t107 * t15 + t112 * t16 + t66 * t51 - t77 * t53;
t72 = -t103 * mrSges(5,2) + t80 * mrSges(5,3);
t73 = t103 * mrSges(5,1) - t81 * mrSges(5,3);
t123 = -m(5) * t119 + t49 * mrSges(5,1) - t50 * mrSges(5,2) - t108 * t13 - t113 * t12 + t80 * t72 - t81 * t73;
t83 = -t103 * mrSges(4,2) + t89 * mrSges(4,3);
t85 = t103 * mrSges(4,1) - t90 * mrSges(4,3);
t121 = -m(4) * t122 + t69 * mrSges(4,1) - t70 * mrSges(4,2) + t89 * t83 - t90 * t85 + t123;
t97 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t136;
t135 = qJD(1) * t115;
t98 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t135;
t143 = (t110 * t97 - t115 * t98) * qJD(1) + m(3) * (-t117 * pkin(7) + t125) - t96 * mrSges(3,1) + t95 * mrSges(3,2) - t121;
t22 = -t102 * pkin(4) - t101 * pkin(9) + t81 * t57 - t144;
t124 = t28 * mrSges(7,1) + t45 * t41 - m(7) * (-t37 * pkin(5) - t64 * pkin(10) + t67 * t54 + t22) - t29 * mrSges(7,2) - t46 * t42;
t120 = m(6) * t22 - t37 * mrSges(6,1) + t38 * mrSges(6,2) - t66 * t52 + t67 * t53 - t124;
t56 = -t80 * mrSges(5,1) + t81 * mrSges(5,2);
t14 = m(5) * t144 + t102 * mrSges(5,1) - t50 * mrSges(5,3) + t103 * t72 - t81 * t56 - t120;
t82 = -t89 * mrSges(4,1) + t90 * mrSges(4,2);
t9 = m(5) * t133 - t102 * mrSges(5,2) + t49 * mrSges(5,3) - t103 * t73 - t108 * t12 + t113 * t13 + t80 * t56;
t6 = m(4) * t128 + t102 * mrSges(4,1) - t70 * mrSges(4,3) + t103 * t83 + t105 * t9 + t106 * t14 - t90 * t82;
t7 = m(4) * t139 - t102 * mrSges(4,2) + t69 * mrSges(4,3) - t103 * t85 - t105 * t14 + t106 * t9 + t89 * t82;
t94 = (-mrSges(3,1) * t115 + mrSges(3,2) * t110) * qJD(1);
t4 = m(3) * (-t115 * g(3) - t138) - t95 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t94 * t136 + qJD(2) * t98 + t109 * t7 + t114 * t6;
t5 = m(3) * t132 - qJDD(2) * mrSges(3,2) + t96 * mrSges(3,3) - qJD(2) * t97 - t109 * t6 + t114 * t7 + t94 * t135;
t142 = t110 * t5 + t115 * t4;
t8 = m(2) * t131 + qJDD(1) * mrSges(2,1) - t117 * mrSges(2,2) - t143;
t1 = m(2) * t127 - t117 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t110 * t4 + t115 * t5;
t2 = [-m(1) * g(1) + t116 * t1 - t111 * t8, t1, t5, t7, t9, t13, t16; -m(1) * g(2) + t111 * t1 + t116 * t8, t8, t4, t6, t14, t12, t15; (-m(1) - m(2)) * g(3) + t142, -m(2) * g(3) + t142, t143, -t121, -t123, t120, -t124;];
f_new  = t2;
