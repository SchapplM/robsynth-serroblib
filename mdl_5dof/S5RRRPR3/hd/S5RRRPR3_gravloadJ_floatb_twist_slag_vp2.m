% Calculate Gravitation load on the joints for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% m [6x1]
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
% Datum: 2022-01-20 11:44
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:42:28
% EndTime: 2022-01-20 11:42:29
% DurationCPUTime: 0.28s
% Computational Cost: add. (289->60), mult. (227->62), div. (0->0), fcn. (173->10), ass. (0->34)
t35 = qJ(3) + pkin(9);
t29 = qJ(5) + t35;
t21 = sin(t29);
t22 = cos(t29);
t67 = t22 * mrSges(6,1) - t21 * mrSges(6,2);
t27 = sin(t35);
t28 = cos(t35);
t38 = sin(qJ(3));
t66 = -t28 * mrSges(5,1) + t38 * mrSges(4,2) + t27 * mrSges(5,2) - t67;
t40 = cos(qJ(3));
t65 = -t40 * mrSges(4,1) - mrSges(3,1) + t66;
t63 = mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t62 = m(5) + m(6);
t36 = qJ(1) + qJ(2);
t30 = sin(t36);
t31 = cos(t36);
t60 = g(1) * t31 + g(2) * t30;
t37 = -qJ(4) - pkin(7);
t32 = t40 * pkin(3);
t52 = pkin(4) * t28 + t32;
t51 = t31 * pkin(2) + t30 * pkin(7);
t50 = m(5) * pkin(3) + mrSges(4,1);
t34 = pkin(8) - t37;
t4 = pkin(2) + t52;
t49 = t30 * t34 + t31 * t4;
t26 = t32 + pkin(2);
t48 = t31 * t26 - t30 * t37;
t46 = mrSges(6,1) * t21 + mrSges(6,2) * t22;
t44 = t63 * t30 + t65 * t31;
t42 = (m(4) * pkin(2) + m(5) * t26 + m(6) * t4 - t65) * t30 + (-m(4) * pkin(7) + m(5) * t37 - m(6) * t34 + t63) * t31;
t41 = cos(qJ(1));
t39 = sin(qJ(1));
t33 = t41 * pkin(1);
t1 = [(t39 * mrSges(2,2) - m(4) * (t33 + t51) - m(5) * (t33 + t48) - m(6) * (t33 + t49) + (-m(3) * pkin(1) - mrSges(2,1)) * t41 + t44) * g(2) + (t41 * mrSges(2,2) + (mrSges(2,1) + (m(3) + m(4) + t62) * pkin(1)) * t39 + t42) * g(1), (-m(4) * t51 - m(5) * t48 - m(6) * t49 + t44) * g(2) + t42 * g(1), (-m(6) * t52 - t50 * t40 + t66) * g(3) + t60 * (-m(6) * (-t38 * pkin(3) - pkin(4) * t27) + mrSges(5,1) * t27 + mrSges(4,2) * t40 + mrSges(5,2) * t28 + t50 * t38 + t46), t62 * (-g(1) * t30 + g(2) * t31), -g(3) * t67 + t60 * t46];
taug = t1(:);
