% Calculate potential energy for
% S6RRPRRP1
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
% Datum: 2019-03-09 11:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRP1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP1_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP1_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:38:54
% EndTime: 2019-03-09 11:38:54
% DurationCPUTime: 0.46s
% Computational Cost: add. (246->72), mult. (200->61), div. (0->0), fcn. (166->10), ass. (0->33)
t27 = cos(qJ(5));
t52 = -m(6) * pkin(4) - m(7) * (t27 * pkin(5) + pkin(4)) - mrSges(5,1);
t51 = m(6) * pkin(9) - m(7) * (-qJ(6) - pkin(9)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t50 = m(7) * pkin(5);
t49 = -m(2) - m(3);
t48 = -mrSges(6,1) - mrSges(7,1);
t47 = mrSges(6,2) + mrSges(7,2);
t46 = -m(5) - m(6) - m(7);
t45 = -m(1) - m(4) + t49;
t20 = qJ(2) + pkin(10);
t16 = qJ(4) + t20;
t10 = sin(t16);
t11 = cos(t16);
t28 = cos(qJ(2));
t13 = t28 * pkin(2) + pkin(1);
t14 = sin(t20);
t15 = cos(t20);
t25 = sin(qJ(2));
t44 = -m(3) * pkin(1) - m(4) * t13 - t28 * mrSges(3,1) - t15 * mrSges(4,1) + t25 * mrSges(3,2) + t14 * mrSges(4,2) - t51 * t10 + t52 * t11 - mrSges(2,1);
t23 = -qJ(3) - pkin(7);
t43 = m(3) * pkin(7) - m(4) * t23 - mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3);
t24 = sin(qJ(5));
t26 = sin(qJ(1));
t42 = t26 * t24;
t41 = t26 * t27;
t29 = cos(qJ(1));
t40 = t29 * t24;
t39 = t29 * t27;
t21 = pkin(6) + r_base(3);
t37 = t25 * pkin(2) + t21;
t19 = -pkin(8) + t23;
t7 = pkin(3) * t15 + t13;
t1 = (-m(1) * r_base(3) - m(4) * t37 - t25 * mrSges(3,1) - t14 * mrSges(4,1) - t28 * mrSges(3,2) - t15 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) + t49 * t21 + t46 * (pkin(3) * t14 + t37) + t51 * t11 + (t47 * t24 + t48 * t27 + t52) * t10) * g(3) + (t40 * t50 - mrSges(1,2) + t46 * (t29 * t19 + t26 * t7 + r_base(2)) + t48 * (t11 * t41 - t40) - t47 * (-t11 * t42 - t39) + t45 * r_base(2) + t43 * t29 + t44 * t26) * g(2) + (-t42 * t50 - mrSges(1,1) + t48 * (t11 * t39 + t42) + t46 * (-t26 * t19 + t29 * t7 + r_base(1)) - t47 * (-t11 * t40 + t41) + t45 * r_base(1) - t43 * t26 + t44 * t29) * g(1);
U  = t1;
