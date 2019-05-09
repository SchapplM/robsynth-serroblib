% Calculate vector of cutting forces with Newton-Euler
% S6RRRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2019-05-07 04:11
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPPR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:06:42
% EndTime: 2019-05-07 04:06:55
% DurationCPUTime: 5.91s
% Computational Cost: add. (99845->206), mult. (227733->277), div. (0->0), fcn. (169689->12), ass. (0->105)
t145 = -2 * qJD(4);
t106 = sin(pkin(10));
t138 = cos(pkin(10));
t102 = qJDD(2) + qJDD(3);
t103 = qJD(2) + qJD(3);
t109 = sin(qJ(3));
t113 = cos(qJ(3));
t110 = sin(qJ(2));
t114 = cos(qJ(2));
t134 = qJD(1) * qJD(2);
t116 = qJD(1) ^ 2;
t111 = sin(qJ(1));
t115 = cos(qJ(1));
t127 = -t115 * g(1) - t111 * g(2);
t92 = -t116 * pkin(1) + qJDD(1) * pkin(7) + t127;
t139 = t110 * t92;
t141 = pkin(2) * t116;
t95 = t110 * qJDD(1) + t114 * t134;
t60 = qJDD(2) * pkin(2) - t95 * pkin(8) - t139 + (pkin(8) * t134 + t110 * t141 - g(3)) * t114;
t104 = t114 ^ 2;
t131 = -t110 * g(3) + t114 * t92;
t96 = t114 * qJDD(1) - t110 * t134;
t136 = qJD(1) * t110;
t99 = qJD(2) * pkin(2) - pkin(8) * t136;
t61 = t96 * pkin(8) - qJD(2) * t99 - t104 * t141 + t131;
t129 = -t109 * t61 + t113 * t60;
t89 = (-t109 * t110 + t113 * t114) * qJD(1);
t71 = t89 * qJD(3) + t109 * t96 + t113 * t95;
t90 = (t109 * t114 + t110 * t113) * qJD(1);
t32 = (t103 * t89 - t71) * qJ(4) + (t89 * t90 + t102) * pkin(3) + t129;
t140 = t109 * t60 + t113 * t61;
t70 = -t90 * qJD(3) - t109 * t95 + t113 * t96;
t84 = t103 * pkin(3) - t90 * qJ(4);
t88 = t89 ^ 2;
t36 = -t88 * pkin(3) + t70 * qJ(4) - t103 * t84 + t140;
t81 = t106 * t89 + t138 * t90;
t144 = -t106 * t36 + t138 * t32 + t145 * t81;
t130 = t111 * g(1) - t115 * g(2);
t124 = -qJDD(1) * pkin(1) - t130;
t121 = -t96 * pkin(2) + t99 * t136 + (-pkin(8) * t104 - pkin(7)) * t116 + t124;
t105 = sin(pkin(11));
t107 = cos(pkin(11));
t118 = -t70 * pkin(3) - t88 * qJ(4) + t90 * t84 + qJDD(4) + t121;
t108 = sin(qJ(6));
t112 = cos(qJ(6));
t101 = t103 ^ 2;
t80 = t106 * t90 - t138 * t89;
t132 = t106 * t32 + t138 * t36 + t145 * t80;
t55 = t80 * pkin(4) - t81 * qJ(5);
t23 = -t101 * pkin(4) + t102 * qJ(5) - t80 * t55 + t132;
t49 = t106 * t71 - t138 * t70;
t50 = t106 * t70 + t138 * t71;
t26 = (t103 * t80 - t50) * qJ(5) + (t103 * t81 + t49) * pkin(4) + t118;
t67 = t105 * t103 + t107 * t81;
t128 = -0.2e1 * qJD(5) * t67 - t105 * t23 + t107 * t26;
t43 = t105 * t102 + t107 * t50;
t66 = t107 * t103 - t105 * t81;
t17 = (t66 * t80 - t43) * pkin(9) + (t66 * t67 + t49) * pkin(5) + t128;
t133 = 0.2e1 * qJD(5) * t66 + t105 * t26 + t107 * t23;
t42 = t107 * t102 - t105 * t50;
t53 = t80 * pkin(5) - t67 * pkin(9);
t65 = t66 ^ 2;
t18 = -t65 * pkin(5) + t42 * pkin(9) - t80 * t53 + t133;
t44 = -t108 * t67 + t112 * t66;
t29 = t44 * qJD(6) + t108 * t42 + t112 * t43;
t45 = t108 * t66 + t112 * t67;
t33 = -t44 * mrSges(7,1) + t45 * mrSges(7,2);
t77 = qJD(6) + t80;
t39 = -t77 * mrSges(7,2) + t44 * mrSges(7,3);
t47 = qJDD(6) + t49;
t15 = m(7) * (-t108 * t18 + t112 * t17) - t29 * mrSges(7,3) + t47 * mrSges(7,1) - t45 * t33 + t77 * t39;
t28 = -t45 * qJD(6) - t108 * t43 + t112 * t42;
t40 = t77 * mrSges(7,1) - t45 * mrSges(7,3);
t16 = m(7) * (t108 * t17 + t112 * t18) + t28 * mrSges(7,3) - t47 * mrSges(7,2) + t44 * t33 - t77 * t40;
t48 = -t66 * mrSges(6,1) + t67 * mrSges(6,2);
t51 = -t80 * mrSges(6,2) + t66 * mrSges(6,3);
t12 = m(6) * t128 + t49 * mrSges(6,1) - t43 * mrSges(6,3) + t108 * t16 + t112 * t15 - t67 * t48 + t80 * t51;
t52 = t80 * mrSges(6,1) - t67 * mrSges(6,3);
t13 = m(6) * t133 - t49 * mrSges(6,2) + t42 * mrSges(6,3) - t108 * t15 + t112 * t16 + t66 * t48 - t80 * t52;
t73 = -t103 * mrSges(5,2) - t80 * mrSges(5,3);
t74 = t103 * mrSges(5,1) - t81 * mrSges(5,3);
t122 = m(5) * t118 + t49 * mrSges(5,1) + t50 * mrSges(5,2) + t105 * t13 + t107 * t12 + t80 * t73 + t81 * t74;
t83 = -t103 * mrSges(4,2) + t89 * mrSges(4,3);
t85 = t103 * mrSges(4,1) - t90 * mrSges(4,3);
t120 = -m(4) * t121 + t70 * mrSges(4,1) - t71 * mrSges(4,2) + t89 * t83 - t90 * t85 - t122;
t97 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t136;
t135 = qJD(1) * t114;
t98 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t135;
t143 = (t110 * t97 - t114 * t98) * qJD(1) + m(3) * (-t116 * pkin(7) + t124) - t96 * mrSges(3,1) + t95 * mrSges(3,2) - t120;
t22 = -t102 * pkin(4) - t101 * qJ(5) + t81 * t55 + qJDD(5) - t144;
t123 = t28 * mrSges(7,1) + t44 * t39 - m(7) * (-t42 * pkin(5) - t65 * pkin(9) + t67 * t53 + t22) - t29 * mrSges(7,2) - t45 * t40;
t119 = m(6) * t22 - t42 * mrSges(6,1) + t43 * mrSges(6,2) - t66 * t51 + t67 * t52 - t123;
t56 = t80 * mrSges(5,1) + t81 * mrSges(5,2);
t14 = m(5) * t144 + t102 * mrSges(5,1) - t50 * mrSges(5,3) + t103 * t73 - t81 * t56 - t119;
t82 = -t89 * mrSges(4,1) + t90 * mrSges(4,2);
t9 = m(5) * t132 - t102 * mrSges(5,2) - t49 * mrSges(5,3) - t103 * t74 - t105 * t12 + t107 * t13 - t80 * t56;
t6 = m(4) * t129 + t102 * mrSges(4,1) - t71 * mrSges(4,3) + t103 * t83 + t106 * t9 + t138 * t14 - t90 * t82;
t7 = m(4) * t140 - t102 * mrSges(4,2) + t70 * mrSges(4,3) - t103 * t85 - t106 * t14 + t138 * t9 + t89 * t82;
t94 = (-mrSges(3,1) * t114 + mrSges(3,2) * t110) * qJD(1);
t4 = m(3) * (-t114 * g(3) - t139) - t95 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t94 * t136 + qJD(2) * t98 + t109 * t7 + t113 * t6;
t5 = m(3) * t131 - qJDD(2) * mrSges(3,2) + t96 * mrSges(3,3) - qJD(2) * t97 - t109 * t6 + t113 * t7 + t94 * t135;
t142 = t110 * t5 + t114 * t4;
t8 = m(2) * t130 + qJDD(1) * mrSges(2,1) - t116 * mrSges(2,2) - t143;
t1 = m(2) * t127 - t116 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t110 * t4 + t114 * t5;
t2 = [-m(1) * g(1) + t115 * t1 - t111 * t8, t1, t5, t7, t9, t13, t16; -m(1) * g(2) + t111 * t1 + t115 * t8, t8, t4, t6, t14, t12, t15; (-m(1) - m(2)) * g(3) + t142, -m(2) * g(3) + t142, t143, -t120, t122, t119, -t123;];
f_new  = t2;
