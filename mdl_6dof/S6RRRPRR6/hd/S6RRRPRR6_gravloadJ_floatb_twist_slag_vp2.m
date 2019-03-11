% Calculate Gravitation load on the joints for
% S6RRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
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
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:25:09
% EndTime: 2019-03-09 18:25:11
% DurationCPUTime: 0.93s
% Computational Cost: add. (614->142), mult. (637->170), div. (0->0), fcn. (590->12), ass. (0->73)
t110 = mrSges(4,3) + mrSges(5,3) + mrSges(6,3) + mrSges(7,3);
t48 = qJ(3) + pkin(11);
t40 = cos(t48);
t53 = cos(qJ(3));
t44 = t53 * pkin(3);
t30 = pkin(4) * t40 + t44;
t41 = qJ(5) + t48;
t36 = cos(t41);
t23 = pkin(5) * t36 + t30;
t17 = pkin(2) + t23;
t28 = pkin(2) + t30;
t37 = qJ(6) + t41;
t32 = sin(t37);
t33 = cos(t37);
t35 = sin(t41);
t38 = t44 + pkin(2);
t39 = sin(t48);
t50 = sin(qJ(3));
t109 = -m(4) * pkin(2) - m(5) * t38 - m(6) * t28 - m(7) * t17 - t53 * mrSges(4,1) - t40 * mrSges(5,1) - t36 * mrSges(6,1) - t33 * mrSges(7,1) + t50 * mrSges(4,2) + t39 * mrSges(5,2) + t35 * mrSges(6,2) + t32 * mrSges(7,2);
t49 = -qJ(4) - pkin(8);
t47 = -pkin(9) + t49;
t42 = -pkin(10) + t47;
t108 = -m(4) * pkin(8) + m(5) * t49 + m(6) * t47 + m(7) * t42 - t110;
t52 = sin(qJ(1));
t55 = cos(qJ(1));
t107 = g(1) * t55 + g(2) * t52;
t96 = m(5) * pkin(3);
t106 = mrSges(4,1) + t96;
t67 = -mrSges(7,1) * t32 - mrSges(7,2) * t33;
t105 = mrSges(6,1) * t35 + mrSges(6,2) * t36 - t67;
t54 = cos(qJ(2));
t79 = t54 * t55;
t15 = -t35 * t79 + t52 * t36;
t16 = t52 * t35 + t36 * t79;
t7 = -t32 * t79 + t52 * t33;
t8 = t52 * t32 + t33 * t79;
t93 = t7 * mrSges(7,1) - t8 * mrSges(7,2);
t104 = -t15 * mrSges(6,1) + t16 * mrSges(6,2) - t93;
t80 = t52 * t54;
t13 = t35 * t80 + t36 * t55;
t14 = t35 * t55 - t36 * t80;
t5 = t32 * t80 + t33 * t55;
t6 = t32 * t55 - t33 * t80;
t94 = -t5 * mrSges(7,1) + t6 * mrSges(7,2);
t103 = t13 * mrSges(6,1) - t14 * mrSges(6,2) - t94;
t102 = m(5) + m(6) + m(7);
t101 = -m(3) - m(4) - t102;
t92 = pkin(3) * t50;
t29 = pkin(4) * t39 + t92;
t91 = pkin(5) * t35;
t22 = t29 + t91;
t99 = -m(6) * t29 - m(7) * t22;
t51 = sin(qJ(2));
t70 = t54 * mrSges(3,1) - t51 * mrSges(3,2);
t98 = t110 * t51 + mrSges(2,1) + t70;
t97 = mrSges(2,2) - mrSges(3,3) + t99;
t95 = m(7) * pkin(5);
t88 = g(3) * t51;
t86 = t50 * t55;
t81 = t52 * t50;
t71 = t54 * pkin(2) + t51 * pkin(8);
t66 = t54 * t17 - t51 * t42;
t65 = t54 * t28 - t51 * t47;
t64 = t54 * t38 - t51 * t49;
t26 = -t50 * t79 + t52 * t53;
t24 = t50 * t80 + t53 * t55;
t27 = t53 * t79 + t81;
t25 = -t53 * t80 + t86;
t21 = t52 * t39 + t40 * t79;
t20 = -t39 * t79 + t52 * t40;
t19 = t39 * t55 - t40 * t80;
t18 = t39 * t80 + t40 * t55;
t1 = [(-t81 * t96 - t27 * mrSges(4,1) - t21 * mrSges(5,1) - t16 * mrSges(6,1) - t8 * mrSges(7,1) - t26 * mrSges(4,2) - t20 * mrSges(5,2) - t15 * mrSges(6,2) - t7 * mrSges(7,2) + t101 * (t55 * pkin(1) + t52 * pkin(7)) + t97 * t52 + (-m(4) * t71 - m(5) * t64 - m(6) * t65 - m(7) * t66 - t98) * t55) * g(2) + (-t86 * t96 - t25 * mrSges(4,1) - t19 * mrSges(5,1) - t14 * mrSges(6,1) - t6 * mrSges(7,1) - t24 * mrSges(4,2) - t18 * mrSges(5,2) - t13 * mrSges(6,2) - t5 * mrSges(7,2) + (m(3) * pkin(1) - m(4) * (-pkin(1) - t71) - m(5) * (-pkin(1) - t64) - m(6) * (-pkin(1) - t65) - m(7) * (-pkin(1) - t66) + t98) * t52 + (t101 * pkin(7) + t97) * t55) * g(1), -g(3) * t70 + (t109 * g(3) + t107 * (mrSges(3,2) + t108)) * t54 + (t108 * g(3) + t107 * (mrSges(3,1) - t109)) * t51 (m(5) * t92 + mrSges(4,1) * t50 + mrSges(5,1) * t39 + mrSges(4,2) * t53 + mrSges(5,2) * t40 + t105 - t99) * t88 + (-t25 * mrSges(4,2) + t18 * mrSges(5,1) - t19 * mrSges(5,2) - m(6) * (-t29 * t80 - t55 * t30) - m(7) * (-t22 * t80 - t55 * t23) + t106 * t24 + t103) * g(2) + (t27 * mrSges(4,2) - t20 * mrSges(5,1) + t21 * mrSges(5,2) - m(6) * (-t29 * t79 + t52 * t30) - m(7) * (-t22 * t79 + t52 * t23) - t106 * t26 + t104) * g(1) (g(3) * t54 - t107 * t51) * t102 (m(7) * t91 + t105) * t88 + (t13 * t95 + t103) * g(2) + (-t15 * t95 + t104) * g(1), -g(1) * t93 - g(2) * t94 - t67 * t88];
taug  = t1(:);
