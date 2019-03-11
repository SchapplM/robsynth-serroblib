% Calculate potential energy for
% S6RPPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-03-09 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRP1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPRRP1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP1_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP1_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:57:32
% EndTime: 2019-03-09 01:57:32
% DurationCPUTime: 0.46s
% Computational Cost: add. (245->72), mult. (186->61), div. (0->0), fcn. (152->10), ass. (0->33)
t26 = cos(qJ(5));
t51 = -m(6) * pkin(4) - m(7) * (t26 * pkin(5) + pkin(4)) - mrSges(5,1);
t50 = m(6) * pkin(8) - m(7) * (-qJ(6) - pkin(8)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t49 = m(7) * pkin(5);
t48 = -m(1) - m(2);
t47 = -m(3) - m(4);
t46 = -mrSges(6,1) - mrSges(7,1);
t45 = mrSges(6,2) + mrSges(7,2);
t44 = -m(5) - m(6) - m(7);
t18 = pkin(10) + qJ(4);
t10 = sin(t18);
t12 = cos(t18);
t20 = sin(pkin(10));
t21 = cos(pkin(10));
t43 = -m(4) * pkin(2) - t21 * mrSges(4,1) + t20 * mrSges(4,2) - t50 * t10 + t51 * t12 - mrSges(3,1);
t42 = m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t19 = qJ(1) + pkin(9);
t11 = sin(t19);
t24 = sin(qJ(5));
t41 = t11 * t24;
t40 = t11 * t26;
t13 = cos(t19);
t39 = t13 * t24;
t38 = t13 * t26;
t37 = pkin(6) + r_base(3);
t25 = sin(qJ(1));
t36 = t25 * pkin(1) + r_base(2);
t27 = cos(qJ(1));
t35 = t27 * pkin(1) + r_base(1);
t14 = qJ(2) + t37;
t23 = -pkin(7) - qJ(3);
t8 = t21 * pkin(3) + pkin(2);
t1 = (-m(1) * r_base(3) - m(2) * t37 - t20 * mrSges(4,1) - t21 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t47 * t14 + t44 * (t20 * pkin(3) + t14) + t50 * t12 + (t45 * t24 + t46 * t26 + t51) * t10) * g(3) + (t39 * t49 - t25 * mrSges(2,1) - t27 * mrSges(2,2) - mrSges(1,2) + t48 * r_base(2) + t47 * t36 + t44 * (t11 * t8 + t13 * t23 + t36) + t46 * (t12 * t40 - t39) - t45 * (-t12 * t41 - t38) + t42 * t13 + t43 * t11) * g(2) + (-t41 * t49 - t27 * mrSges(2,1) + t25 * mrSges(2,2) - mrSges(1,1) + t48 * r_base(1) + t46 * (t12 * t38 + t41) + t47 * t35 + t44 * (-t11 * t23 + t13 * t8 + t35) - t45 * (-t12 * t39 + t40) - t42 * t11 + t43 * t13) * g(1);
U  = t1;
