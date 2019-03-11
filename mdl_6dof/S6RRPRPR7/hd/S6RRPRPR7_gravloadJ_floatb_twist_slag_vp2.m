% Calculate Gravitation load on the joints for
% S6RRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 10:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:45:20
% EndTime: 2019-03-09 10:45:23
% DurationCPUTime: 1.06s
% Computational Cost: add. (447->145), mult. (730->171), div. (0->0), fcn. (733->10), ass. (0->81)
t42 = sin(qJ(6));
t46 = cos(qJ(6));
t119 = m(7) * pkin(5) + t46 * mrSges(7,1) - t42 * mrSges(7,2) + mrSges(6,1);
t118 = mrSges(6,2) - mrSges(7,3);
t45 = sin(qJ(1));
t48 = cos(qJ(2));
t86 = qJ(3) * t48;
t28 = t45 * t86;
t44 = sin(qJ(2));
t47 = cos(qJ(4));
t35 = pkin(4) * t47 + pkin(3);
t98 = -pkin(2) - t35;
t78 = t44 * t98;
t117 = t45 * t78 + t28;
t116 = -m(4) - m(5);
t115 = (mrSges(3,1) + mrSges(4,1)) * t48 + (-mrSges(3,2) + mrSges(4,3)) * t44;
t114 = -mrSges(2,1) - t115;
t113 = pkin(9) * m(7) - t118;
t112 = pkin(8) * m(5) + mrSges(2,2) - mrSges(4,2) - mrSges(3,3) + mrSges(5,3) + mrSges(6,3);
t111 = m(6) + m(7);
t49 = cos(qJ(1));
t39 = t49 * pkin(7);
t41 = -qJ(5) - pkin(8);
t43 = sin(qJ(4));
t110 = (-pkin(1) + t98 * t48 + (-pkin(4) * t43 - qJ(3)) * t44) * t45 + t49 * t41 + t39;
t85 = qJ(4) + pkin(10);
t73 = sin(t85);
t74 = cos(t85);
t12 = t44 * t74 - t48 * t73;
t109 = g(1) * t49 + g(2) * t45;
t91 = t47 * t48;
t96 = t43 * t44;
t58 = t91 + t96;
t10 = t58 * t49;
t54 = t49 * t73;
t55 = t49 * t74;
t5 = -t44 * t54 - t48 * t55;
t6 = -t44 * t55 + t48 * t54;
t92 = t44 * t49;
t83 = t47 * t92;
t90 = t48 * t49;
t84 = t43 * t90;
t9 = -t83 + t84;
t108 = t9 * mrSges(5,1) + t10 * mrSges(5,2) - t118 * t5 + t119 * t6;
t3 = t12 * t45;
t11 = t44 * t73 + t48 * t74;
t4 = t11 * t45;
t93 = t44 * t47;
t95 = t43 * t48;
t59 = -t93 + t95;
t7 = t59 * t45;
t8 = t58 * t45;
t107 = t7 * mrSges(5,1) + t8 * mrSges(5,2) + t118 * t4 - t119 * t3;
t106 = t58 * mrSges(5,1) - t59 * mrSges(5,2) + t119 * t11 + t118 * t12;
t105 = -pkin(2) - pkin(3);
t104 = pkin(4) * t45;
t103 = pkin(9) * t12;
t38 = t48 * pkin(2);
t36 = t44 * qJ(3);
t88 = t38 + t36;
t87 = t49 * pkin(1) + t45 * pkin(7);
t31 = pkin(4) * t96;
t82 = t44 * t105;
t21 = t95 * t104;
t80 = pkin(9) * t4 - t21;
t24 = pkin(4) * t84;
t79 = -t5 * pkin(9) - t24;
t77 = -pkin(1) - t36;
t76 = t48 * t35 + t31 + t88;
t75 = pkin(2) * t90 + t49 * t36 + t87;
t72 = -pkin(4) * t91 - t31;
t66 = -t4 * t46 - t42 * t49;
t65 = t4 * t42 - t46 * t49;
t53 = t49 * t31 + t35 * t90 + t45 * t41 + t75;
t51 = t48 * mrSges(4,3) + (-m(4) * pkin(2) - mrSges(4,1)) * t44;
t29 = t49 * t86;
t23 = pkin(4) * t83;
t20 = t93 * t104;
t2 = -t42 * t45 - t46 * t5;
t1 = t42 * t5 - t45 * t46;
t13 = [(-m(3) * t87 - m(4) * t75 - m(5) * (pkin(3) * t90 + t75) - t10 * mrSges(5,1) + t9 * mrSges(5,2) - m(6) * t53 + t5 * mrSges(6,1) - m(7) * (-pkin(5) * t5 + t53) - t2 * mrSges(7,1) - t1 * mrSges(7,2) - t113 * t6 + t114 * t49 + t112 * t45) * g(2) + (t8 * mrSges(5,1) - t7 * mrSges(5,2) - t110 * m(6) + t4 * mrSges(6,1) - t66 * mrSges(7,1) - t65 * mrSges(7,2) - (-pkin(5) * t4 + t110) * m(7) - t113 * t3 + (-m(3) + t116) * t39 + (m(3) * pkin(1) - m(4) * (t77 - t38) - m(5) * (t105 * t48 + t77) - t114) * t45 + t112 * t49) * g(1), t109 * (mrSges(3,1) * t44 + mrSges(3,2) * t48) + (-m(4) * t28 - t51 * t45 - m(5) * (t45 * t82 + t28) - m(6) * (t117 + t21) - m(7) * (t117 - t80) - t107) * g(2) + (-m(4) * t29 - t51 * t49 - m(5) * (t49 * t82 + t29) - m(6) * (t49 * t78 + t24 + t29) - m(7) * (t98 * t92 + t29 - t79) - t108) * g(1) + (-m(4) * t88 - m(5) * (pkin(3) * t48 + t88) - m(6) * t76 - m(7) * (t76 - t103) - t106 - t115) * g(3) (t48 * g(3) - t109 * t44) * (t111 - t116) (-m(6) * t72 - m(7) * (t72 + t103) + t106) * g(3) + (-m(6) * (-t21 + t20) - m(7) * (t20 + t80) + t107) * g(2) + (-m(6) * (-t24 + t23) - m(7) * (t23 + t79) + t108) * g(1), t111 * (g(1) * t45 - g(2) * t49) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * (-t65 * mrSges(7,1) + t66 * mrSges(7,2)) - g(3) * (-mrSges(7,1) * t42 - mrSges(7,2) * t46) * t12];
taug  = t13(:);
