% Calculate potential energy for
% S6RPRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRR3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPRR3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR3_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR3_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:40:37
% EndTime: 2019-03-09 03:40:37
% DurationCPUTime: 0.43s
% Computational Cost: add. (266->68), mult. (212->54), div. (0->0), fcn. (182->12), ass. (0->28)
t19 = pkin(11) + qJ(5);
t11 = cos(t19);
t21 = sin(pkin(11));
t22 = cos(pkin(11));
t13 = qJ(6) + t19;
t6 = sin(t13);
t7 = cos(t13);
t8 = t22 * pkin(4) + pkin(3);
t9 = sin(t19);
t52 = -m(5) * pkin(3) - m(6) * t8 - mrSges(5,1) * t22 - mrSges(6,1) * t11 + mrSges(5,2) * t21 + mrSges(6,2) * t9 - mrSges(4,1) - m(7) * (pkin(5) * t11 + t8) - mrSges(7,1) * t7 + mrSges(7,2) * t6;
t23 = -pkin(8) - qJ(4);
t51 = -m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3) + m(7) * (-pkin(9) + t23) - mrSges(7,3) + m(6) * t23 - mrSges(6,3);
t50 = -m(1) - m(2);
t49 = m(4) + m(5) + m(6) + m(7);
t24 = sin(qJ(3));
t26 = cos(qJ(3));
t45 = t51 * t24 + t52 * t26 - mrSges(3,1);
t43 = pkin(4) * t21;
t44 = m(6) * t43 + m(7) * (pkin(5) * t9 + t43) - mrSges(3,2) + mrSges(4,3) + t21 * mrSges(5,1) + t22 * mrSges(5,2) + t9 * mrSges(6,1) + t11 * mrSges(6,2) + t6 * mrSges(7,1) + t7 * mrSges(7,2);
t42 = pkin(6) + r_base(3);
t25 = sin(qJ(1));
t41 = t25 * pkin(1) + r_base(2);
t27 = cos(qJ(1));
t40 = t27 * pkin(1) + r_base(1);
t20 = qJ(1) + pkin(10);
t12 = cos(t20);
t10 = sin(t20);
t1 = (-m(1) * r_base(3) - m(2) * t42 - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + (-m(3) - t49) * (qJ(2) + t42) - t51 * t26 + t52 * t24) * g(3) + (-m(3) * t41 - t25 * mrSges(2,1) - t27 * mrSges(2,2) - mrSges(1,2) + t50 * r_base(2) - t49 * (t10 * pkin(2) + t41) + (t49 * pkin(7) + t44) * t12 + t45 * t10) * g(2) + (-m(3) * t40 - t27 * mrSges(2,1) + t25 * mrSges(2,2) - mrSges(1,1) + t50 * r_base(1) - t49 * (t12 * pkin(2) + t10 * pkin(7) + t40) + t45 * t12 - t44 * t10) * g(1);
U  = t1;
