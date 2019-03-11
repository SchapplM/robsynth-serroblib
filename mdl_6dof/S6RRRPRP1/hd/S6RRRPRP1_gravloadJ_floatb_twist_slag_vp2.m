% Calculate Gravitation load on the joints for
% S6RRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:30:44
% EndTime: 2019-03-09 16:30:47
% DurationCPUTime: 0.89s
% Computational Cost: add. (532->111), mult. (505->123), div. (0->0), fcn. (433->10), ass. (0->64)
t122 = mrSges(6,1) + mrSges(7,1);
t120 = -mrSges(6,3) - mrSges(7,3);
t44 = cos(qJ(5));
t121 = t122 * t44;
t114 = mrSges(6,2) + mrSges(7,2);
t39 = qJ(2) + qJ(3);
t34 = pkin(10) + t39;
t29 = sin(t34);
t119 = t29 * t114;
t30 = cos(t34);
t35 = sin(t39);
t36 = cos(t39);
t118 = mrSges(4,1) * t35 + mrSges(5,1) * t29 + mrSges(4,2) * t36 + mrSges(5,2) * t30;
t117 = t121 * t29;
t43 = sin(qJ(1));
t46 = cos(qJ(1));
t109 = g(1) * t46 + g(2) * t43;
t116 = t36 * mrSges(4,1) + t30 * mrSges(5,1) - t35 * mrSges(4,2) + (-mrSges(5,2) - t120) * t29;
t101 = m(7) * pkin(5);
t32 = pkin(5) * t44 + pkin(4);
t40 = -qJ(6) - pkin(9);
t59 = -t29 * t40 + t30 * t32;
t65 = t30 * pkin(4) + t29 * pkin(9);
t58 = -t29 * t32 - t30 * t40;
t98 = pkin(4) * t29;
t99 = pkin(3) * t35;
t112 = -m(7) * (t58 - t99) - m(6) * (-t98 - t99) + t117;
t41 = sin(qJ(5));
t84 = t41 * t46;
t87 = t30 * t46;
t111 = -t84 * t119 + t120 * t87;
t85 = t41 * t43;
t88 = t30 * t43;
t110 = -t85 * t119 + t120 * t88;
t108 = m(5) + m(6) + m(7);
t106 = t101 + t122;
t47 = -pkin(8) - pkin(7);
t105 = -m(3) * pkin(7) + m(4) * t47 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t104 = -t116 + (t114 * t41 - t121) * t30;
t103 = m(6) * t98 - m(7) * t58 + t117;
t45 = cos(qJ(2));
t37 = t45 * pkin(2);
t42 = sin(qJ(2));
t64 = t45 * mrSges(3,1) - t42 * mrSges(3,2);
t102 = mrSges(2,1) + m(4) * (t37 + pkin(1)) + m(3) * pkin(1) + t64 + t116;
t100 = pkin(2) * t42;
t31 = pkin(3) * t36;
t83 = t43 * t44;
t81 = t44 * t46;
t78 = t31 + t37;
t72 = t31 + t65;
t66 = t31 + t59;
t3 = -t30 * t84 + t83;
t1 = t30 * t85 + t81;
t38 = -qJ(4) + t47;
t22 = pkin(9) * t87;
t21 = pkin(9) * t88;
t14 = -t99 - t100;
t13 = pkin(1) + t78;
t7 = t46 * t14;
t6 = t43 * t14;
t4 = t30 * t81 + t85;
t2 = -t30 * t83 + t84;
t5 = [(-t85 * t101 - t108 * (t46 * t13 - t38 * t43) - t122 * t4 - t114 * t3 + t105 * t43 + (-m(6) * t65 - m(7) * t59 - t102) * t46) * g(2) + (-t122 * t2 - t114 * t1 + (-t41 * t101 + t108 * t38 + t105) * t46 + (m(5) * t13 - m(6) * (-t13 - t65) - m(7) * (-t13 - t59) + t102) * t43) * g(1) (-m(6) * (t21 + t6) - m(7) * t6 + t103 * t43 + t110) * g(2) + (-m(6) * (t22 + t7) - m(7) * t7 + t103 * t46 + t111) * g(1) + (-t64 - m(4) * t37 - m(5) * t78 - m(6) * (t37 + t72) - m(7) * (t37 + t66) + t104) * g(3) + t109 * (m(4) * t100 - m(5) * t14 + mrSges(3,1) * t42 + mrSges(3,2) * t45 + t118) (-m(6) * t21 + t112 * t43 + t110) * g(2) + (-m(6) * t22 + t112 * t46 + t111) * g(1) + (-m(5) * t31 - m(6) * t72 - m(7) * t66 + t104) * g(3) + (m(5) * t99 + t118) * t109 (-g(1) * t43 + g(2) * t46) * t108 (t106 * t41 + t114 * t44) * g(3) * t29 + (t106 * t1 - t114 * t2) * g(2) + (-t106 * t3 + t114 * t4) * g(1) (g(3) * t30 - t109 * t29) * m(7)];
taug  = t5(:);
