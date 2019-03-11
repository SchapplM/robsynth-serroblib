% Calculate potential energy for
% S6RRRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRR6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRR6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR6_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR6_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:24:59
% EndTime: 2019-03-09 18:25:00
% DurationCPUTime: 0.48s
% Computational Cost: add. (266->68), mult. (254->54), div. (0->0), fcn. (228->12), ass. (0->28)
t23 = qJ(3) + pkin(11);
t15 = qJ(5) + t23;
t10 = cos(t15);
t29 = cos(qJ(3));
t12 = t29 * pkin(3) + pkin(2);
t13 = sin(t23);
t14 = cos(t23);
t26 = sin(qJ(3));
t3 = pkin(4) * t14 + t12;
t11 = qJ(6) + t15;
t5 = sin(t11);
t6 = cos(t11);
t9 = sin(t15);
t56 = -m(4) * pkin(2) - m(5) * t12 - m(6) * t3 - t29 * mrSges(4,1) - t14 * mrSges(5,1) - t10 * mrSges(6,1) + t26 * mrSges(4,2) + t13 * mrSges(5,2) + t9 * mrSges(6,2) - mrSges(3,1) - m(7) * (pkin(5) * t10 + t3) - t6 * mrSges(7,1) + t5 * mrSges(7,2);
t25 = -qJ(4) - pkin(8);
t22 = -pkin(9) + t25;
t55 = mrSges(3,2) + m(7) * (-pkin(10) + t22) - mrSges(7,3) + m(6) * t22 - mrSges(6,3) + m(5) * t25 - mrSges(5,3) - m(4) * pkin(8) - mrSges(4,3);
t54 = -m(1) - m(2);
t49 = m(3) + m(4) + m(5) + m(6) + m(7);
t27 = sin(qJ(2));
t30 = cos(qJ(2));
t48 = t55 * t27 + t56 * t30 - mrSges(2,1);
t17 = t26 * pkin(3);
t4 = pkin(4) * t13 + t17;
t47 = m(5) * t17 + m(6) * t4 + m(7) * (pkin(5) * t9 + t4) - mrSges(2,2) + mrSges(3,3) + t13 * mrSges(5,1) + t14 * mrSges(5,2) + t26 * mrSges(4,1) + t29 * mrSges(4,2) + t9 * mrSges(6,1) + t10 * mrSges(6,2) + t5 * mrSges(7,1) + t6 * mrSges(7,2);
t31 = cos(qJ(1));
t28 = sin(qJ(1));
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + (-m(2) - t49) * (pkin(6) + r_base(3)) - t55 * t30 + t56 * t27) * g(3) + (-mrSges(1,2) + t54 * r_base(2) - t49 * (t28 * pkin(1) + r_base(2)) + (t49 * pkin(7) + t47) * t31 + t48 * t28) * g(2) + (-mrSges(1,1) + t54 * r_base(1) - t49 * (t31 * pkin(1) + t28 * pkin(7) + r_base(1)) + t48 * t31 - t47 * t28) * g(1);
U  = t1;
