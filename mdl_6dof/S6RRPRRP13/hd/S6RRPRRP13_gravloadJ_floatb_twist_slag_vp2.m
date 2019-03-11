% Calculate Gravitation load on the joints for
% S6RRPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP13_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP13_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP13_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP13_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:55:36
% EndTime: 2019-03-09 12:55:38
% DurationCPUTime: 1.14s
% Computational Cost: add. (554->124), mult. (1321->173), div. (0->0), fcn. (1529->10), ass. (0->67)
t119 = -mrSges(5,2) + mrSges(7,3);
t82 = m(5) + m(6) + m(7);
t78 = m(4) + t82;
t118 = qJ(3) * t78;
t117 = mrSges(6,1) + mrSges(7,1);
t112 = mrSges(6,2) + mrSges(7,2);
t49 = sin(qJ(4));
t53 = cos(qJ(4));
t116 = m(6) * (pkin(4) * t49 - pkin(10) * t53) - t53 * mrSges(6,3);
t113 = mrSges(3,2) - mrSges(4,3);
t52 = cos(qJ(5));
t44 = pkin(5) * t52 + pkin(4);
t47 = -qJ(6) - pkin(10);
t115 = -m(7) * (t44 * t49 + t47 * t53) - t49 * mrSges(5,1) + t113 + t119 * t53;
t103 = m(7) * pkin(5);
t111 = mrSges(3,1) - mrSges(4,2) + mrSges(5,3);
t48 = sin(qJ(5));
t110 = m(6) * pkin(4) + m(7) * t44 - t112 * t48 + t117 * t52 + mrSges(5,1);
t59 = -m(6) * pkin(10) + m(7) * t47 - mrSges(6,3) - t119;
t109 = -t103 - t117;
t62 = -t48 * t103 - t111;
t108 = pkin(9) * t82 - t62;
t107 = -t113 + t118;
t50 = sin(qJ(2));
t51 = sin(qJ(1));
t54 = cos(qJ(2));
t55 = cos(qJ(1));
t83 = cos(pkin(6));
t71 = t55 * t83;
t31 = t50 * t71 + t51 * t54;
t30 = t50 * t51 - t54 * t71;
t46 = sin(pkin(6));
t95 = t46 * t55;
t63 = -t30 * t49 + t53 * t95;
t106 = t31 * t52 + t48 * t63;
t105 = -t31 * t48 + t52 * t63;
t104 = t115 - t116 - t118;
t102 = t30 * t48;
t72 = t51 * t83;
t32 = t55 * t50 + t54 * t72;
t99 = t32 * t48;
t98 = t46 * t50;
t97 = t46 * t51;
t96 = t46 * t54;
t94 = t48 * t49;
t93 = t48 * t50;
t92 = t49 * t52;
t91 = t50 * t52;
t85 = pkin(2) * t96 + qJ(3) * t98;
t84 = t55 * pkin(1) + pkin(8) * t97;
t33 = -t50 * t72 + t54 * t55;
t80 = t33 * pkin(2) + t84;
t75 = -pkin(1) * t51 + pkin(8) * t95;
t70 = pkin(3) * t97 + t80;
t69 = -t31 * pkin(2) + t75;
t67 = mrSges(2,2) + (-mrSges(4,1) - mrSges(3,3)) * t46;
t14 = t32 * t49 + t53 * t97;
t1 = -t14 * t48 + t33 * t52;
t64 = pkin(3) * t95 + t69;
t15 = t30 * t53 + t49 * t95;
t29 = -t49 * t96 + t53 * t83;
t28 = t49 * t83 + t53 * t96;
t26 = t32 * pkin(2);
t24 = t30 * pkin(2);
t13 = -t32 * t53 + t49 * t97;
t2 = t14 * t52 + t33 * t48;
t3 = [(-t55 * mrSges(2,1) - m(3) * t84 - m(4) * t80 - m(5) * t70 - t14 * mrSges(5,1) - m(6) * (pkin(4) * t14 + t70) - m(7) * (t14 * t44 + t70) - t107 * t32 - t117 * t2 - t112 * t1 + t67 * t51 + t59 * t13 - t108 * t33) * g(2) + (t51 * mrSges(2,1) - m(3) * t75 - m(4) * t69 - m(5) * t64 - t63 * mrSges(5,1) - m(6) * (pkin(4) * t63 + t64) - m(7) * (t44 * t63 + t64) - t117 * t105 + t107 * t30 + t112 * t106 + t67 * t55 + t59 * t15 + t108 * t31) * g(1) (t102 * t103 + m(4) * t24 - t117 * (t31 * t92 - t102) - t82 * (-pkin(9) * t30 - t24) - t112 * (-t30 * t52 - t31 * t94) + t111 * t30 + t104 * t31) * g(2) + (t99 * t103 + m(4) * t26 - t112 * (-t32 * t52 - t33 * t94) - t82 * (-pkin(9) * t32 - t26) - t117 * (t33 * t92 - t99) + t111 * t32 + t104 * t33) * g(1) + (-m(4) * t85 - t116 * t98 - t82 * (pkin(9) * t96 + t85) + (-t117 * (t48 * t54 + t49 * t91) - t112 * (-t49 * t93 + t52 * t54) + t62 * t54 + t115 * t50) * t46) * g(3) (-g(1) * t32 - g(2) * t30 + g(3) * t96) * t78 (t110 * t28 + t59 * t29) * g(3) + (-t110 * t15 - t59 * t63) * g(2) + (t110 * t13 + t59 * t14) * g(1) (-t112 * (-t29 * t52 - t46 * t93) + t109 * (-t29 * t48 + t46 * t91)) * g(3) + (-t112 * t105 + t109 * t106) * g(2) + (t109 * t1 + t112 * t2) * g(1) (-g(1) * t13 + g(2) * t15 - g(3) * t28) * m(7)];
taug  = t3(:);
