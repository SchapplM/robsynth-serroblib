% Calculate potential energy for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPRR4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR4_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:16:11
% EndTime: 2022-01-23 09:16:11
% DurationCPUTime: 0.59s
% Computational Cost: add. (170->80), mult. (182->69), div. (0->0), fcn. (158->10), ass. (0->34)
t23 = qJ(3) + pkin(6);
t55 = -mrSges(4,3) - mrSges(5,3) - m(4) * qJ(3) - m(5) * t23 + mrSges(3,2) + m(6) * (-pkin(7) - t23) - mrSges(6,3);
t17 = pkin(9) + qJ(4);
t10 = cos(t17);
t11 = qJ(5) + t17;
t6 = sin(t11);
t7 = cos(t11);
t21 = cos(pkin(9));
t8 = t21 * pkin(3) + pkin(2);
t54 = -m(4) * pkin(2) - m(5) * t8 - mrSges(3,1) - m(6) * (pkin(4) * t10 + t8) - t7 * mrSges(6,1) + t6 * mrSges(6,2);
t52 = -m(3) - m(6);
t20 = sin(pkin(8));
t22 = cos(pkin(8));
t51 = -mrSges(2,1) + (-m(4) - m(5)) * pkin(1) + t54 * t22 + t55 * t20;
t50 = -m(2) - m(5);
t48 = -m(6) + t50;
t19 = sin(pkin(9));
t43 = t19 * pkin(3);
t45 = m(5) * (qJ(2) + t43) - mrSges(2,2) + mrSges(3,3) + t6 * mrSges(6,1) + t7 * mrSges(6,2);
t24 = sin(qJ(1));
t9 = sin(t17);
t42 = t24 * t9;
t25 = cos(qJ(1));
t41 = t25 * t9;
t40 = t24 * t10;
t39 = t24 * t19;
t37 = t24 * t21;
t36 = t25 * t10;
t35 = t25 * t19;
t33 = t25 * t21;
t30 = -t25 * qJ(2) + r_base(2);
t14 = t24 * pkin(1);
t3 = pkin(4) * t9 + t43;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + (-m(3) - m(4) + t48) * (pkin(5) + r_base(3)) - t55 * t22 + (-t21 * mrSges(4,1) - t10 * mrSges(5,1) + t19 * mrSges(4,2) + t9 * mrSges(5,2) + t54) * t20) * g(3) + (-mrSges(1,2) - m(3) * (t14 + t30) - m(4) * t30 - (t22 * t37 - t35) * mrSges(4,1) - (-t22 * t39 - t33) * mrSges(4,2) - (t22 * t40 - t41) * mrSges(5,1) - (-t22 * t42 - t36) * mrSges(5,2) - m(6) * t14 + (-m(1) + t48) * r_base(2) + (-m(6) * (-qJ(2) - t3) + t45) * t25 + t51 * t24) * g(2) + (-mrSges(1,1) - (t22 * t33 + t39) * mrSges(4,1) - (-t22 * t35 + t37) * mrSges(4,2) - (t22 * t36 + t42) * mrSges(5,1) - (-t22 * t41 + t40) * mrSges(5,2) + (-m(1) + t50) * r_base(1) + (-m(6) * t3 - t45) * t24 + (-m(4) + t52) * (t24 * qJ(2) + r_base(1)) + (t52 * pkin(1) + t51) * t25) * g(1);
U = t1;
