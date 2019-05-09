% Calculate vector of cutting forces with Newton-Euler
% S6RRPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-06 19:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRRR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 19:51:50
% EndTime: 2019-05-06 19:52:08
% DurationCPUTime: 5.88s
% Computational Cost: add. (99170->208), mult. (229185->276), div. (0->0), fcn. (173261->12), ass. (0->106)
t116 = cos(qJ(2));
t136 = qJD(1) * t116;
t100 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t136;
t111 = sin(qJ(2));
t118 = qJD(1) ^ 2;
t105 = t116 ^ 2;
t112 = sin(qJ(1));
t117 = cos(qJ(1));
t132 = t112 * g(1) - t117 * g(2);
t126 = -qJDD(1) * pkin(1) - t132;
t137 = qJD(1) * t111;
t135 = qJD(1) * qJD(2);
t97 = t116 * qJDD(1) - t111 * t135;
t98 = qJD(2) * pkin(2) - qJ(3) * t137;
t123 = -t97 * pkin(2) + qJDD(3) + t98 * t137 + (-qJ(3) * t105 - pkin(7)) * t118 + t126;
t109 = sin(qJ(5));
t114 = cos(qJ(5));
t108 = sin(qJ(6));
t113 = cos(qJ(6));
t104 = qJD(2) + qJD(4);
t102 = t104 ^ 2;
t103 = qJDD(2) + qJDD(4);
t110 = sin(qJ(4));
t115 = cos(qJ(4));
t106 = sin(pkin(11));
t107 = cos(pkin(11));
t128 = -t117 * g(1) - t112 * g(2);
t93 = -t118 * pkin(1) + qJDD(1) * pkin(7) + t128;
t138 = t111 * t93;
t141 = pkin(2) * t118;
t96 = t111 * qJDD(1) + t116 * t135;
t61 = qJDD(2) * pkin(2) - t96 * qJ(3) - t138 + (qJ(3) * t135 + t111 * t141 - g(3)) * t116;
t133 = -t111 * g(3) + t116 * t93;
t62 = t97 * qJ(3) - qJD(2) * t98 - t105 * t141 + t133;
t91 = (t106 * t116 + t107 * t111) * qJD(1);
t129 = -0.2e1 * qJD(3) * t91 - t106 * t62 + t107 * t61;
t81 = t106 * t97 + t107 * t96;
t90 = (-t106 * t111 + t107 * t116) * qJD(1);
t32 = (qJD(2) * t90 - t81) * pkin(8) + (t90 * t91 + qJDD(2)) * pkin(3) + t129;
t134 = 0.2e1 * qJD(3) * t90 + t106 * t61 + t107 * t62;
t80 = -t106 * t96 + t107 * t97;
t84 = qJD(2) * pkin(3) - t91 * pkin(8);
t89 = t90 ^ 2;
t36 = -t89 * pkin(3) + t80 * pkin(8) - qJD(2) * t84 + t134;
t139 = t110 * t32 + t115 * t36;
t73 = -t110 * t91 + t115 * t90;
t74 = t110 * t90 + t115 * t91;
t57 = -t73 * pkin(4) - t74 * pkin(9);
t23 = -t102 * pkin(4) + t103 * pkin(9) + t73 * t57 + t139;
t120 = -t80 * pkin(3) - t89 * pkin(8) + t91 * t84 + t123;
t47 = -t74 * qJD(4) - t110 * t81 + t115 * t80;
t48 = t73 * qJD(4) + t110 * t80 + t115 * t81;
t29 = (-t104 * t73 - t48) * pkin(9) + (t104 * t74 - t47) * pkin(4) + t120;
t131 = -t109 * t23 + t114 * t29;
t66 = t114 * t104 - t109 * t74;
t38 = t66 * qJD(5) + t109 * t103 + t114 * t48;
t46 = qJDD(5) - t47;
t67 = t109 * t104 + t114 * t74;
t72 = qJD(5) - t73;
t17 = (t66 * t72 - t38) * pkin(10) + (t66 * t67 + t46) * pkin(5) + t131;
t140 = t109 * t29 + t114 * t23;
t37 = -t67 * qJD(5) + t114 * t103 - t109 * t48;
t54 = t72 * pkin(5) - t67 * pkin(10);
t64 = t66 ^ 2;
t18 = -t64 * pkin(5) + t37 * pkin(10) - t72 * t54 + t140;
t49 = -t108 * t67 + t113 * t66;
t26 = t49 * qJD(6) + t108 * t37 + t113 * t38;
t50 = t108 * t66 + t113 * t67;
t33 = -t49 * mrSges(7,1) + t50 * mrSges(7,2);
t70 = qJD(6) + t72;
t39 = -t70 * mrSges(7,2) + t49 * mrSges(7,3);
t44 = qJDD(6) + t46;
t15 = m(7) * (-t108 * t18 + t113 * t17) - t26 * mrSges(7,3) + t44 * mrSges(7,1) - t50 * t33 + t70 * t39;
t25 = -t50 * qJD(6) - t108 * t38 + t113 * t37;
t40 = t70 * mrSges(7,1) - t50 * mrSges(7,3);
t16 = m(7) * (t108 * t17 + t113 * t18) + t25 * mrSges(7,3) - t44 * mrSges(7,2) + t49 * t33 - t70 * t40;
t51 = -t66 * mrSges(6,1) + t67 * mrSges(6,2);
t52 = -t72 * mrSges(6,2) + t66 * mrSges(6,3);
t12 = m(6) * t131 + t46 * mrSges(6,1) - t38 * mrSges(6,3) + t108 * t16 + t113 * t15 - t67 * t51 + t72 * t52;
t53 = t72 * mrSges(6,1) - t67 * mrSges(6,3);
t13 = m(6) * t140 - t46 * mrSges(6,2) + t37 * mrSges(6,3) - t108 * t15 + t113 * t16 + t66 * t51 - t72 * t53;
t68 = -t104 * mrSges(5,2) + t73 * mrSges(5,3);
t69 = t104 * mrSges(5,1) - t74 * mrSges(5,3);
t124 = -m(5) * t120 + t47 * mrSges(5,1) - t48 * mrSges(5,2) - t109 * t13 - t114 * t12 + t73 * t68 - t74 * t69;
t82 = -qJD(2) * mrSges(4,2) + t90 * mrSges(4,3);
t83 = qJD(2) * mrSges(4,1) - t91 * mrSges(4,3);
t122 = -m(4) * t123 + t80 * mrSges(4,1) - t81 * mrSges(4,2) + t90 * t82 - t91 * t83 + t124;
t99 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t137;
t143 = -(t116 * t100 - t111 * t99) * qJD(1) + m(3) * (-t118 * pkin(7) + t126) - t97 * mrSges(3,1) + t96 * mrSges(3,2) - t122;
t130 = -t110 * t36 + t115 * t32;
t22 = -t103 * pkin(4) - t102 * pkin(9) + t74 * t57 - t130;
t125 = t25 * mrSges(7,1) + t49 * t39 - m(7) * (-t37 * pkin(5) - t64 * pkin(10) + t67 * t54 + t22) - t26 * mrSges(7,2) - t50 * t40;
t121 = m(6) * t22 - t37 * mrSges(6,1) + t38 * mrSges(6,2) - t66 * t52 + t67 * t53 - t125;
t56 = -t73 * mrSges(5,1) + t74 * mrSges(5,2);
t14 = m(5) * t130 + t103 * mrSges(5,1) - t48 * mrSges(5,3) + t104 * t68 - t74 * t56 - t121;
t77 = -t90 * mrSges(4,1) + t91 * mrSges(4,2);
t9 = m(5) * t139 - t103 * mrSges(5,2) + t47 * mrSges(5,3) - t104 * t69 - t109 * t12 + t114 * t13 + t73 * t56;
t6 = m(4) * t129 + qJDD(2) * mrSges(4,1) - t81 * mrSges(4,3) + qJD(2) * t82 + t110 * t9 + t115 * t14 - t91 * t77;
t7 = m(4) * t134 - qJDD(2) * mrSges(4,2) + t80 * mrSges(4,3) - qJD(2) * t83 - t110 * t14 + t115 * t9 + t90 * t77;
t95 = (-mrSges(3,1) * t116 + mrSges(3,2) * t111) * qJD(1);
t4 = m(3) * (-t116 * g(3) - t138) - t96 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t95 * t137 + qJD(2) * t100 + t106 * t7 + t107 * t6;
t5 = m(3) * t133 - qJDD(2) * mrSges(3,2) + t97 * mrSges(3,3) - qJD(2) * t99 - t106 * t6 + t107 * t7 + t95 * t136;
t142 = t111 * t5 + t116 * t4;
t8 = m(2) * t132 + qJDD(1) * mrSges(2,1) - t118 * mrSges(2,2) - t143;
t1 = m(2) * t128 - t118 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t111 * t4 + t116 * t5;
t2 = [-m(1) * g(1) + t117 * t1 - t112 * t8, t1, t5, t7, t9, t13, t16; -m(1) * g(2) + t112 * t1 + t117 * t8, t8, t4, t6, t14, t12, t15; (-m(1) - m(2)) * g(3) + t142, -m(2) * g(3) + t142, t143, -t122, -t124, t121, -t125;];
f_new  = t2;
