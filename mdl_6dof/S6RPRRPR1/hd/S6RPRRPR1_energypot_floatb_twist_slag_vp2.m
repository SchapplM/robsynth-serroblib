% Calculate potential energy for
% S6RPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 04:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPR1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRPR1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR1_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR1_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:57:17
% EndTime: 2019-03-09 04:57:17
% DurationCPUTime: 0.37s
% Computational Cost: add. (244->70), mult. (162->54), div. (0->0), fcn. (124->12), ass. (0->33)
t22 = sin(qJ(6));
t25 = cos(qJ(6));
t50 = -m(7) * pkin(5) - t25 * mrSges(7,1) + t22 * mrSges(7,2) - mrSges(6,1);
t49 = -m(7) * pkin(9) + mrSges(6,2) - mrSges(7,3);
t48 = -m(1) - m(2);
t47 = -m(3) - m(4);
t46 = m(6) + m(7);
t45 = -m(5) + t47;
t21 = qJ(3) + qJ(4);
t13 = sin(t21);
t14 = cos(t21);
t23 = sin(qJ(3));
t26 = cos(qJ(3));
t11 = pkin(11) + t21;
t5 = sin(t11);
t6 = cos(t11);
t8 = t26 * pkin(3) + pkin(2);
t43 = -m(4) * pkin(2) - m(5) * t8 - t26 * mrSges(4,1) - t14 * mrSges(5,1) + t23 * mrSges(4,2) + t13 * mrSges(5,2) + t49 * t5 + t50 * t6 - mrSges(3,1);
t28 = -pkin(8) - pkin(7);
t42 = m(4) * pkin(7) - m(5) * t28 + t22 * mrSges(7,1) + t25 * mrSges(7,2) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t41 = pkin(6) + r_base(3);
t24 = sin(qJ(1));
t40 = t24 * pkin(1) + r_base(2);
t27 = cos(qJ(1));
t39 = t27 * pkin(1) + r_base(1);
t12 = qJ(2) + t41;
t36 = t23 * pkin(3) + t12;
t20 = qJ(1) + pkin(10);
t19 = -qJ(5) + t28;
t10 = cos(t20);
t9 = sin(t20);
t3 = pkin(4) * t14 + t8;
t1 = (-m(1) * r_base(3) - m(2) * t41 - m(5) * t36 - t23 * mrSges(4,1) - t13 * mrSges(5,1) - t26 * mrSges(4,2) - t14 * mrSges(5,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - t46 * (pkin(4) * t13 + t36) - t49 * t6 + t50 * t5 + t47 * t12) * g(3) + (-t24 * mrSges(2,1) - t27 * mrSges(2,2) - mrSges(1,2) + t48 * r_base(2) - t46 * (t10 * t19 + t9 * t3 + t40) + t45 * t40 + t43 * t9 + t42 * t10) * g(2) + (-t27 * mrSges(2,1) + t24 * mrSges(2,2) - mrSges(1,1) + t48 * r_base(1) - t46 * (t10 * t3 + t39) + t45 * t39 + (t46 * t19 - t42) * t9 + t43 * t10) * g(1);
U  = t1;
