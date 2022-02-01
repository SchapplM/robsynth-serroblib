% Calculate Gravitation load on the joints for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2022-01-20 12:09
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:07:40
% EndTime: 2022-01-20 12:07:41
% DurationCPUTime: 0.31s
% Computational Cost: add. (327->63), mult. (253->64), div. (0->0), fcn. (194->10), ass. (0->36)
t36 = qJ(3) + qJ(4);
t32 = qJ(5) + t36;
t25 = sin(t32);
t26 = cos(t32);
t69 = t26 * mrSges(6,1) - t25 * mrSges(6,2);
t28 = sin(t36);
t30 = cos(t36);
t68 = -t30 * mrSges(5,1) + t28 * mrSges(5,2) - t69;
t38 = sin(qJ(3));
t67 = t38 * mrSges(4,2) + t68;
t40 = cos(qJ(3));
t66 = -t40 * mrSges(4,1) - mrSges(3,1) + t67;
t64 = mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t37 = qJ(1) + qJ(2);
t29 = sin(t37);
t31 = cos(t37);
t62 = g(1) * t31 + g(2) * t29;
t42 = -pkin(8) - pkin(7);
t24 = pkin(4) * t30;
t55 = t31 * pkin(2) + t29 * pkin(7);
t33 = t40 * pkin(3);
t54 = t24 + t33;
t53 = m(5) * pkin(3) + mrSges(4,1);
t35 = pkin(9) - t42;
t4 = pkin(2) + t54;
t52 = t29 * t35 + t31 * t4;
t27 = t33 + pkin(2);
t51 = t31 * t27 - t29 * t42;
t50 = mrSges(6,1) * t25 + mrSges(6,2) * t26;
t48 = mrSges(5,2) * t30 + t50;
t45 = t64 * t29 + t66 * t31;
t43 = (m(4) * pkin(2) + m(5) * t27 + m(6) * t4 - t66) * t29 + (-m(4) * pkin(7) + m(5) * t42 - m(6) * t35 + t64) * t31;
t41 = cos(qJ(1));
t39 = sin(qJ(1));
t34 = t41 * pkin(1);
t1 = [(t39 * mrSges(2,2) - m(4) * (t34 + t55) - m(5) * (t34 + t51) - m(6) * (t34 + t52) + (-m(3) * pkin(1) - mrSges(2,1)) * t41 + t45) * g(2) + (t41 * mrSges(2,2) + (mrSges(2,1) + (m(3) + m(4) + m(5) + m(6)) * pkin(1)) * t39 + t43) * g(1), (-m(4) * t55 - m(5) * t51 - m(6) * t52 + t45) * g(2) + t43 * g(1), (-m(6) * t54 - t53 * t40 + t67) * g(3) + t62 * (-m(6) * (-t38 * pkin(3) - pkin(4) * t28) + mrSges(5,1) * t28 + mrSges(4,2) * t40 + t53 * t38 + t48), (-m(6) * t24 + t68) * g(3) + t62 * ((m(6) * pkin(4) + mrSges(5,1)) * t28 + t48), -g(3) * t69 + t62 * t50];
taug = t1(:);
