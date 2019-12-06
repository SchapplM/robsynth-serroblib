% Calculate Gravitation load on the joints for
% S5PRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRP6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP6_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:50:47
% EndTime: 2019-12-05 16:50:49
% DurationCPUTime: 0.49s
% Computational Cost: add. (237->72), mult. (352->93), div. (0->0), fcn. (329->8), ass. (0->38)
t72 = mrSges(5,1) + mrSges(6,1);
t71 = mrSges(5,2) - mrSges(6,3);
t73 = m(5) + m(6);
t24 = qJ(3) + qJ(4);
t22 = sin(t24);
t23 = cos(t24);
t27 = sin(qJ(3));
t29 = cos(qJ(3));
t70 = m(4) * pkin(2) + t29 * mrSges(4,1) - t27 * mrSges(4,2) - t71 * t22 + t72 * t23 + mrSges(3,1);
t69 = -m(4) * pkin(6) + mrSges(3,2) - mrSges(6,2) - mrSges(4,3) - mrSges(5,3);
t25 = sin(pkin(8));
t26 = cos(pkin(8));
t30 = cos(qJ(2));
t52 = t27 * t30;
t68 = t25 * t29 - t26 * t52;
t50 = t30 * t22;
t13 = -t25 * t23 + t26 * t50;
t49 = t30 * t23;
t14 = t25 * t22 + t26 * t49;
t66 = t72 * t13 + t71 * t14;
t11 = t26 * t23 + t25 * t50;
t12 = -t26 * t22 + t25 * t49;
t65 = t72 * t11 + t71 * t12;
t59 = pkin(3) * t27;
t28 = sin(qJ(2));
t56 = g(3) * t28;
t55 = mrSges(5,2) * t23;
t51 = t29 * t30;
t46 = (-m(6) * qJ(5) - mrSges(6,3)) * t23 * t28;
t45 = -t11 * pkin(4) + t12 * qJ(5);
t43 = -t13 * pkin(4) + t14 * qJ(5);
t39 = pkin(4) * t23 + qJ(5) * t22;
t38 = t68 * pkin(3);
t37 = -t25 * t52 - t26 * t29;
t35 = t37 * pkin(3);
t31 = -pkin(7) - pkin(6);
t21 = t29 * pkin(3) + pkin(2);
t1 = [(-m(2) - m(3) - m(4) - t73) * g(3), (-t73 * (t30 * t21 - t28 * t31) + (-m(6) * t39 - t70) * t30 + t69 * t28) * g(3) + (g(1) * t26 + g(2) * t25) * ((t73 * t31 + t69) * t30 + (m(5) * t21 - m(6) * (-t21 - t39) + t70) * t28), -g(3) * ((m(6) * (-pkin(4) * t22 - t59) - t22 * mrSges(6,1)) * t28 - t46) + (m(5) * t59 + mrSges(4,1) * t27 + mrSges(5,1) * t22 + mrSges(4,2) * t29 + t55) * t56 + (-t37 * mrSges(4,1) - (-t25 * t51 + t26 * t27) * mrSges(4,2) - m(5) * t35 - m(6) * (t35 + t45) + t65) * g(2) + (-t68 * mrSges(4,1) - (-t25 * t27 - t26 * t51) * mrSges(4,2) - m(5) * t38 - m(6) * (t38 + t43) + t66) * g(1), ((t55 + (m(6) * pkin(4) + t72) * t22) * t28 + t46) * g(3) + (-m(6) * t45 + t65) * g(2) + (-m(6) * t43 + t66) * g(1), (-g(1) * t13 - g(2) * t11 - t22 * t56) * m(6)];
taug = t1(:);
