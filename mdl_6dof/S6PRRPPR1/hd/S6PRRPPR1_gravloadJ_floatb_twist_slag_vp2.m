% Calculate Gravitation load on the joints for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:08
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRPPR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:07:56
% EndTime: 2018-11-23 15:07:57
% DurationCPUTime: 0.93s
% Computational Cost: add. (1053->96), mult. (1080->129), div. (0->0), fcn. (1023->18), ass. (0->50)
t38 = pkin(12) + qJ(6);
t34 = sin(t38);
t36 = cos(t38);
t40 = sin(pkin(12));
t43 = cos(pkin(12));
t92 = mrSges(5,1) + m(7) * (pkin(5) * t43 + pkin(4)) + t36 * mrSges(7,1) - t34 * mrSges(7,2) + m(6) * pkin(4) + t43 * mrSges(6,1) - t40 * mrSges(6,2);
t89 = mrSges(5,2) - m(6) * qJ(5) - mrSges(6,3) + m(7) * (-pkin(9) - qJ(5)) - mrSges(7,3);
t39 = qJ(3) + pkin(11);
t35 = sin(t39);
t37 = cos(t39);
t46 = sin(qJ(3));
t48 = cos(qJ(3));
t91 = m(4) * pkin(2) + t48 * mrSges(4,1) - t46 * mrSges(4,2) - t35 * t89 + t37 * t92 + mrSges(3,1);
t98 = m(6) + m(7);
t94 = m(5) + t98;
t78 = pkin(6) + qJ(2);
t66 = cos(t78) / 0.2e1;
t79 = pkin(6) - qJ(2);
t71 = cos(t79);
t23 = t66 - t71 / 0.2e1;
t81 = cos(pkin(6));
t100 = t23 * t46 + t48 * t81;
t65 = sin(t78) / 0.2e1;
t70 = sin(t79);
t22 = t65 - t70 / 0.2e1;
t41 = sin(pkin(10));
t49 = cos(qJ(2));
t80 = cos(pkin(10));
t60 = -t22 * t41 + t49 * t80;
t42 = sin(pkin(6));
t83 = t41 * t42;
t99 = -t46 * t60 + t48 * t83;
t61 = t22 * t80 + t41 * t49;
t73 = t42 * t80;
t96 = -t46 * t61 - t48 * t73;
t90 = -m(4) * pkin(8) - t34 * mrSges(7,1) - t43 * mrSges(6,2) - t36 * mrSges(7,2) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) + (-m(7) * pkin(5) - mrSges(6,1)) * t40;
t52 = t71 / 0.2e1 + t66;
t47 = sin(qJ(2));
t44 = -qJ(4) - pkin(8);
t33 = pkin(3) * t48 + pkin(2);
t21 = t65 + t70 / 0.2e1;
t14 = t41 * t52 + t47 * t80;
t11 = t41 * t47 - t52 * t80;
t8 = -t23 * t37 + t35 * t81;
t7 = -t23 * t35 - t37 * t81;
t4 = t35 * t83 + t37 * t60;
t3 = t35 * t60 - t37 * t83;
t2 = -t35 * t73 + t37 * t61;
t1 = t35 * t61 + t37 * t73;
t5 = [(-m(2) - m(3) - m(4) - t94) * g(3) (-t94 * (t21 * t33 + t23 * t44) - t90 * t23 - t91 * t21) * g(3) + (-t94 * (-t11 * t33 - t44 * t61) + t90 * t61 + t91 * t11) * g(2) + (-t94 * (-t14 * t33 - t44 * t60) + t90 * t60 + t91 * t14) * g(1) (-t100 * mrSges(4,1) - (t23 * t48 - t46 * t81) * mrSges(4,2) + t89 * t8 + t92 * t7) * g(3) + (-(t46 * t73 - t48 * t61) * mrSges(4,2) - mrSges(4,1) * t96 + t89 * t2 + t92 * t1) * g(2) + (-t99 * mrSges(4,1) - (-t46 * t83 - t48 * t60) * mrSges(4,2) + t89 * t4 + t92 * t3) * g(1) + (-g(1) * t99 - t96 * g(2) - g(3) * t100) * t94 * pkin(3), t94 * (-g(1) * t14 - g(2) * t11 + g(3) * t21) t98 * (-g(1) * t3 - g(2) * t1 - g(3) * t7) -g(1) * ((t14 * t36 - t34 * t4) * mrSges(7,1) + (-t14 * t34 - t36 * t4) * mrSges(7,2)) - g(2) * ((t11 * t36 - t2 * t34) * mrSges(7,1) + (-t11 * t34 - t2 * t36) * mrSges(7,2)) - g(3) * ((-t21 * t36 - t34 * t8) * mrSges(7,1) + (t21 * t34 - t36 * t8) * mrSges(7,2))];
taug  = t5(:);
