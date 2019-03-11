% Calculate potential energy for
% S6RRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 11:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRP3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP3_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP3_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:47:13
% EndTime: 2019-03-09 11:47:14
% DurationCPUTime: 0.48s
% Computational Cost: add. (247->74), mult. (226->62), div. (0->0), fcn. (196->10), ass. (0->34)
t26 = cos(qJ(4));
t10 = t26 * pkin(4) + pkin(3);
t21 = qJ(4) + qJ(5);
t15 = cos(t21);
t23 = sin(qJ(4));
t56 = -m(6) * t10 - m(7) * (pkin(5) * t15 + t10) - mrSges(4,1) - m(5) * pkin(3) - t26 * mrSges(5,1) + t23 * mrSges(5,2);
t29 = -pkin(9) - pkin(8);
t55 = m(6) * t29 + m(7) * (-qJ(6) + t29) + mrSges(4,2) - mrSges(6,3) - mrSges(7,3) - m(5) * pkin(8) - mrSges(5,3);
t54 = -m(2) - m(3);
t53 = -mrSges(6,1) - mrSges(7,1);
t52 = mrSges(6,2) + mrSges(7,2);
t51 = -m(1) + t54;
t50 = -m(4) - m(6) - m(7);
t49 = -m(5) + t50;
t19 = qJ(2) + pkin(10);
t12 = sin(t19);
t13 = cos(t19);
t24 = sin(qJ(2));
t27 = cos(qJ(2));
t47 = -m(3) * pkin(1) - t27 * mrSges(3,1) + t24 * mrSges(3,2) + t55 * t12 + t56 * t13 - mrSges(2,1);
t14 = sin(t21);
t45 = pkin(4) * t23;
t46 = m(3) * pkin(7) + m(6) * t45 + m(7) * (pkin(5) * t14 + t45) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + t23 * mrSges(5,1) + t26 * mrSges(5,2);
t25 = sin(qJ(1));
t44 = t14 * t25;
t28 = cos(qJ(1));
t43 = t14 * t28;
t42 = t15 * t25;
t41 = t15 * t28;
t20 = pkin(6) + r_base(3);
t11 = pkin(2) * t27 + pkin(1);
t40 = t28 * t11 + r_base(1);
t22 = -qJ(3) - pkin(7);
t1 = (-m(1) * r_base(3) - t24 * mrSges(3,1) - t27 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + t54 * t20 + t49 * (t24 * pkin(2) + t20) - t55 * t13 + (t14 * t52 + t15 * t53 + t56) * t12) * g(3) + (-mrSges(1,2) + t53 * (t13 * t42 - t43) - t52 * (-t13 * t44 - t41) + t51 * r_base(2) + t49 * (t25 * t11 + t28 * t22 + r_base(2)) + t46 * t28 + t47 * t25) * g(2) + (-m(5) * t40 - mrSges(1,1) + t53 * (t13 * t41 + t44) - t52 * (-t13 * t43 + t42) + t51 * r_base(1) + t50 * (-t25 * t22 + t40) + (m(5) * t22 - t46) * t25 + t47 * t28) * g(1);
U  = t1;
