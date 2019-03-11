% Calculate potential energy for
% S6RPRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRRP6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP6_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP6_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:14:08
% EndTime: 2019-03-09 06:14:09
% DurationCPUTime: 0.48s
% Computational Cost: add. (247->74), mult. (226->62), div. (0->0), fcn. (196->10), ass. (0->34)
t27 = cos(qJ(4));
t11 = t27 * pkin(4) + pkin(3);
t21 = qJ(4) + qJ(5);
t15 = cos(t21);
t25 = sin(qJ(4));
t56 = -m(6) * t11 - m(7) * (pkin(5) * t15 + t11) - mrSges(4,1) - m(5) * pkin(3) - t27 * mrSges(5,1) + t25 * mrSges(5,2);
t29 = -pkin(9) - pkin(8);
t55 = m(6) * t29 + m(7) * (-qJ(6) + t29) + mrSges(4,2) - mrSges(6,3) - mrSges(7,3) - m(5) * pkin(8) - mrSges(5,3);
t54 = -m(2) - m(3);
t53 = -mrSges(6,1) - mrSges(7,1);
t52 = mrSges(6,2) + mrSges(7,2);
t51 = -m(1) + t54;
t50 = -m(4) - m(6) - m(7);
t49 = -m(5) + t50;
t19 = pkin(10) + qJ(3);
t12 = sin(t19);
t13 = cos(t19);
t22 = sin(pkin(10));
t23 = cos(pkin(10));
t47 = -m(3) * pkin(1) - t23 * mrSges(3,1) + t22 * mrSges(3,2) + t55 * t12 + t56 * t13 - mrSges(2,1);
t14 = sin(t21);
t45 = pkin(4) * t25;
t46 = m(3) * qJ(2) + m(6) * t45 + m(7) * (pkin(5) * t14 + t45) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + t25 * mrSges(5,1) + t27 * mrSges(5,2);
t26 = sin(qJ(1));
t44 = t14 * t26;
t28 = cos(qJ(1));
t43 = t14 * t28;
t42 = t15 * t26;
t41 = t15 * t28;
t20 = pkin(6) + r_base(3);
t9 = pkin(2) * t23 + pkin(1);
t40 = t28 * t9 + r_base(1);
t24 = -pkin(7) - qJ(2);
t1 = (-m(1) * r_base(3) - t22 * mrSges(3,1) - t23 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + t54 * t20 + t49 * (t22 * pkin(2) + t20) - t55 * t13 + (t14 * t52 + t15 * t53 + t56) * t12) * g(3) + (-mrSges(1,2) + t53 * (t13 * t42 - t43) - t52 * (-t13 * t44 - t41) + t51 * r_base(2) + t49 * (t28 * t24 + t26 * t9 + r_base(2)) + t46 * t28 + t47 * t26) * g(2) + (-m(5) * t40 - mrSges(1,1) + t53 * (t13 * t41 + t44) - t52 * (-t13 * t43 + t42) + t51 * r_base(1) + t50 * (-t26 * t24 + t40) + (m(5) * t24 - t46) * t26 + t47 * t28) * g(1);
U  = t1;
