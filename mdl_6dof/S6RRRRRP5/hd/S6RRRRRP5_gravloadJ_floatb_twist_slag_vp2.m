% Calculate Gravitation load on the joints for
% S6RRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:18:11
% EndTime: 2019-03-10 01:18:13
% DurationCPUTime: 1.07s
% Computational Cost: add. (649->138), mult. (702->162), div. (0->0), fcn. (657->10), ass. (0->72)
t118 = mrSges(6,1) + mrSges(7,1);
t117 = mrSges(6,2) + mrSges(7,2);
t119 = mrSges(4,3) + mrSges(5,3) + mrSges(6,3) + mrSges(7,3);
t47 = qJ(3) + qJ(4);
t39 = cos(t47);
t34 = pkin(4) * t39;
t51 = cos(qJ(3));
t43 = t51 * pkin(3);
t32 = t34 + t43;
t41 = qJ(5) + t47;
t36 = cos(t41);
t33 = pkin(5) * t36;
t23 = t33 + t32;
t17 = pkin(2) + t23;
t30 = pkin(2) + t32;
t35 = sin(t41);
t37 = t43 + pkin(2);
t38 = sin(t47);
t48 = sin(qJ(3));
t116 = -m(4) * pkin(2) - m(5) * t37 - m(6) * t30 - m(7) * t17 - t51 * mrSges(4,1) - t39 * mrSges(5,1) + t48 * mrSges(4,2) + t38 * mrSges(5,2) + t117 * t35 - t118 * t36;
t54 = -pkin(9) - pkin(8);
t46 = -pkin(10) + t54;
t40 = -qJ(6) + t46;
t115 = -m(4) * pkin(8) + m(5) * t54 + m(6) * t46 + m(7) * t40 - t119;
t114 = mrSges(6,1) * t35 + t117 * t36;
t50 = sin(qJ(1));
t53 = cos(qJ(1));
t113 = g(1) * t53 + g(2) * t50;
t99 = m(5) * pkin(3);
t112 = mrSges(4,1) + t99;
t52 = cos(qJ(2));
t79 = t53 * t52;
t11 = -t35 * t79 + t36 * t50;
t12 = t35 * t50 + t36 * t79;
t109 = -t118 * t11 + t117 * t12;
t80 = t50 * t52;
t10 = t35 * t53 - t36 * t80;
t9 = t35 * t80 + t36 * t53;
t108 = -t117 * t10 + t118 * t9;
t107 = mrSges(5,1) * t38 + mrSges(7,1) * t35 + mrSges(5,2) * t39 + t114;
t20 = -t38 * t79 + t39 * t50;
t21 = t38 * t50 + t39 * t79;
t106 = -t20 * mrSges(5,1) + t21 * mrSges(5,2) + t109;
t18 = t38 * t80 + t39 * t53;
t19 = t38 * t53 - t39 * t80;
t105 = t18 * mrSges(5,1) - t19 * mrSges(5,2) + t108;
t104 = -m(3) - m(4) - m(5) - m(6) - m(7);
t91 = pkin(4) * t38;
t28 = -pkin(5) * t35 - t91;
t92 = pkin(3) * t48;
t22 = -t28 + t92;
t31 = t91 + t92;
t102 = -m(6) * t31 - m(7) * t22;
t49 = sin(qJ(2));
t70 = t52 * mrSges(3,1) - t49 * mrSges(3,2);
t101 = t119 * t49 + mrSges(2,1) + t70;
t100 = mrSges(2,2) - mrSges(3,3) + t102;
t98 = m(6) * pkin(4);
t97 = m(7) * pkin(5);
t88 = g(3) * t49;
t86 = t48 * t50;
t85 = t48 * t53;
t71 = pkin(2) * t52 + pkin(8) * t49;
t65 = t17 * t52 - t40 * t49;
t64 = t30 * t52 - t46 * t49;
t63 = t52 * t37 - t49 * t54;
t26 = -t48 * t79 + t50 * t51;
t24 = t48 * t80 + t51 * t53;
t29 = t33 + t34;
t27 = t51 * t79 + t86;
t25 = -t51 * t80 + t85;
t1 = [(-t86 * t99 - t27 * mrSges(4,1) - t21 * mrSges(5,1) - t26 * mrSges(4,2) - t20 * mrSges(5,2) + t104 * (t53 * pkin(1) + t50 * pkin(7)) + t100 * t50 - t118 * t12 - t117 * t11 + (-m(4) * t71 - m(5) * t63 - m(6) * t64 - m(7) * t65 - t101) * t53) * g(2) + (-t85 * t99 - t25 * mrSges(4,1) - t19 * mrSges(5,1) - t24 * mrSges(4,2) - t18 * mrSges(5,2) - t117 * t9 - t118 * t10 + (m(3) * pkin(1) - m(4) * (-pkin(1) - t71) - m(5) * (-pkin(1) - t63) - m(6) * (-pkin(1) - t64) - m(7) * (-pkin(1) - t65) + t101) * t50 + (t104 * pkin(7) + t100) * t53) * g(1), -g(3) * t70 + (t116 * g(3) + t113 * (mrSges(3,2) + t115)) * t52 + (t115 * g(3) + t113 * (mrSges(3,1) - t116)) * t49 (m(5) * t92 + mrSges(4,1) * t48 + mrSges(4,2) * t51 - t102 + t107) * t88 + (-t25 * mrSges(4,2) - m(6) * (-t31 * t80 - t32 * t53) - m(7) * (-t22 * t80 - t23 * t53) + t112 * t24 + t105) * g(2) + (t27 * mrSges(4,2) - m(6) * (-t31 * t79 + t32 * t50) - m(7) * (-t22 * t79 + t23 * t50) - t112 * t26 + t106) * g(1) (m(6) * t91 - m(7) * t28 + t107) * t88 + (t18 * t98 - m(7) * (t28 * t80 - t29 * t53) + t105) * g(2) + (-t20 * t98 - m(7) * (t28 * t79 + t29 * t50) + t106) * g(1) (-(-mrSges(7,1) - t97) * t35 + t114) * t88 + (t9 * t97 + t108) * g(2) + (-t11 * t97 + t109) * g(1) (g(3) * t52 - t113 * t49) * m(7)];
taug  = t1(:);
