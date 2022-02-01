% Calculate Gravitation load on the joints for
% S5RPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2022-01-20 09:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:48:34
% EndTime: 2022-01-20 09:48:35
% DurationCPUTime: 0.25s
% Computational Cost: add. (252->50), mult. (169->52), div. (0->0), fcn. (125->10), ass. (0->30)
t26 = qJ(4) + qJ(5);
t22 = sin(t26);
t23 = cos(t26);
t53 = t23 * mrSges(6,1) - mrSges(6,2) * t22;
t27 = sin(qJ(4));
t52 = t27 * mrSges(5,2) - t53;
t29 = cos(qJ(4));
t51 = -t29 * mrSges(5,1) - mrSges(4,1) + t52;
t49 = mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t25 = qJ(1) + pkin(9);
t21 = qJ(3) + t25;
t16 = sin(t21);
t17 = cos(t21);
t48 = g(1) * t17 + g(2) * t16;
t42 = t17 * pkin(3) + t16 * pkin(7);
t20 = cos(t25);
t30 = cos(qJ(1));
t41 = t30 * pkin(1) + pkin(2) * t20;
t40 = m(4) + m(5) + m(6);
t39 = m(3) + t40;
t38 = m(6) * pkin(4) + mrSges(5,1);
t18 = pkin(4) * t29 + pkin(3);
t31 = -pkin(8) - pkin(7);
t37 = -t16 * t31 + t17 * t18;
t36 = mrSges(6,1) * t22 + mrSges(6,2) * t23;
t33 = t49 * t16 + t51 * t17;
t32 = (m(5) * pkin(3) + m(6) * t18 - t51) * t16 + (-m(5) * pkin(7) + m(6) * t31 + t49) * t17;
t28 = sin(qJ(1));
t19 = sin(t25);
t1 = [(mrSges(2,2) * t28 - mrSges(3,1) * t20 + mrSges(3,2) * t19 - m(4) * t41 - m(5) * (t41 + t42) - m(6) * (t37 + t41) + (-m(3) * pkin(1) - mrSges(2,1)) * t30 + t33) * g(2) + (mrSges(2,2) * t30 + mrSges(3,2) * t20 + (t40 * pkin(2) + mrSges(3,1)) * t19 + (t39 * pkin(1) + mrSges(2,1)) * t28 + t32) * g(1), -t39 * g(3), (-m(5) * t42 - m(6) * t37 + t33) * g(2) + t32 * g(1), (-t38 * t29 + t52) * g(3) + t48 * (mrSges(5,2) * t29 + t38 * t27 + t36), -g(3) * t53 + t48 * t36];
taug = t1(:);
