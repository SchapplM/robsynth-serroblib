% Calculate Gravitation load on the joints for
% S5PRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR8_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:14:49
% EndTime: 2019-12-05 17:14:52
% DurationCPUTime: 0.69s
% Computational Cost: add. (404->88), mult. (706->129), div. (0->0), fcn. (805->12), ass. (0->50)
t42 = sin(qJ(5));
t45 = cos(qJ(5));
t101 = -t45 * mrSges(6,1) + t42 * mrSges(6,2) - mrSges(5,1);
t100 = mrSges(5,2) - mrSges(6,3);
t43 = sin(qJ(3));
t46 = cos(qJ(3));
t73 = cos(pkin(5));
t41 = sin(pkin(5));
t44 = sin(qJ(2));
t81 = t41 * t44;
t99 = -t43 * t81 + t73 * t46;
t47 = cos(qJ(2));
t40 = sin(pkin(10));
t64 = t40 * t73;
t72 = cos(pkin(10));
t29 = -t44 * t64 + t72 * t47;
t80 = t41 * t46;
t98 = -t29 * t43 + t40 * t80;
t87 = -m(4) * pkin(7) - t42 * mrSges(6,1) - t45 * mrSges(6,2) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t39 = qJ(3) + qJ(4);
t37 = sin(t39);
t38 = cos(t39);
t96 = -m(4) * pkin(2) - t46 * mrSges(4,1) + t43 * mrSges(4,2) - mrSges(3,1) + (-m(6) * pkin(4) + t101) * t38 + (-m(6) * pkin(9) + t100) * t37;
t95 = m(5) + m(6);
t22 = -t37 * t81 + t73 * t38;
t23 = t73 * t37 + t38 * t81;
t91 = t100 * t23 + t101 * t22;
t82 = t40 * t41;
t13 = -t29 * t37 + t38 * t82;
t14 = t29 * t38 + t37 * t82;
t90 = t100 * t14 + t101 * t13;
t55 = t73 * t72;
t27 = t40 * t47 + t44 * t55;
t63 = t41 * t72;
t11 = -t27 * t37 - t38 * t63;
t12 = t27 * t38 - t37 * t63;
t89 = t100 * t12 + t101 * t11;
t79 = t41 * t47;
t67 = t11 * pkin(4) + t12 * pkin(9);
t66 = t13 * pkin(4) + t14 * pkin(9);
t65 = t22 * pkin(4) + t23 * pkin(9);
t61 = t98 * pkin(3);
t56 = t99 * pkin(3);
t51 = -t27 * t43 - t46 * t63;
t50 = t51 * pkin(3);
t48 = -pkin(8) - pkin(7);
t36 = t46 * pkin(3) + pkin(2);
t28 = t72 * t44 + t47 * t64;
t26 = t40 * t44 - t47 * t55;
t1 = [(-m(2) - m(3) - m(4) - t95) * g(3), (-t95 * (-t26 * t36 - t27 * t48) + t87 * t27 - t96 * t26) * g(2) + (-t95 * (-t28 * t36 - t29 * t48) + t87 * t29 - t96 * t28) * g(1) + (-t95 * t36 * t79 + (t96 * t47 + (t95 * t48 + t87) * t44) * t41) * g(3), (-t99 * mrSges(4,1) - (-t73 * t43 - t44 * t80) * mrSges(4,2) - m(5) * t56 - m(6) * (t56 + t65) + t91) * g(3) + (-t51 * mrSges(4,1) - (-t27 * t46 + t43 * t63) * mrSges(4,2) - m(5) * t50 - m(6) * (t50 + t67) + t89) * g(2) + (-t98 * mrSges(4,1) - (-t29 * t46 - t43 * t82) * mrSges(4,2) - m(5) * t61 - m(6) * (t61 + t66) + t90) * g(1), (-m(6) * t65 + t91) * g(3) + (-m(6) * t67 + t89) * g(2) + (-m(6) * t66 + t90) * g(1), -g(1) * ((-t14 * t42 + t28 * t45) * mrSges(6,1) + (-t14 * t45 - t28 * t42) * mrSges(6,2)) - g(2) * ((-t12 * t42 + t26 * t45) * mrSges(6,1) + (-t12 * t45 - t26 * t42) * mrSges(6,2)) - g(3) * ((-t23 * t42 - t45 * t79) * mrSges(6,1) + (-t23 * t45 + t42 * t79) * mrSges(6,2))];
taug = t1(:);
