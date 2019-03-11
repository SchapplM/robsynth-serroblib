% Calculate potential energy for
% S6RPPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRP7_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPRRP7_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP7_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP7_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP7_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:12:40
% EndTime: 2019-03-09 02:12:41
% DurationCPUTime: 0.48s
% Computational Cost: add. (184->78), mult. (194->66), div. (0->0), fcn. (160->8), ass. (0->34)
t49 = mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t48 = m(7) * pkin(5);
t47 = -m(1) - m(2);
t46 = m(3) + m(4);
t45 = -mrSges(6,1) - mrSges(7,1);
t44 = mrSges(6,2) + mrSges(7,2);
t43 = -m(5) - m(6) - m(7);
t15 = sin(pkin(9));
t16 = cos(pkin(9));
t13 = pkin(9) + qJ(4);
t7 = sin(t13);
t8 = cos(t13);
t42 = t15 * mrSges(4,1) + t7 * mrSges(5,1) + t16 * mrSges(4,2) + t49 * t8 - mrSges(2,2) + mrSges(3,3);
t41 = -m(4) * qJ(3) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t40 = pkin(3) * t15;
t19 = sin(qJ(5));
t22 = cos(qJ(1));
t37 = t19 * t22;
t20 = sin(qJ(1));
t36 = t20 * t19;
t21 = cos(qJ(5));
t35 = t20 * t21;
t34 = t21 * t22;
t14 = pkin(6) + r_base(3);
t33 = t20 * pkin(1) + r_base(2);
t32 = pkin(2) + t14;
t31 = t22 * pkin(1) + t20 * qJ(2) + r_base(1);
t30 = -qJ(2) - t40;
t28 = pkin(4) * t7 - pkin(8) * t8;
t17 = -qJ(6) - pkin(8);
t6 = pkin(5) * t21 + pkin(4);
t27 = t17 * t8 + t6 * t7;
t18 = -pkin(7) - qJ(3);
t1 = (-m(1) * r_base(3) - m(4) * t32 - t16 * mrSges(4,1) + t15 * mrSges(4,2) - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * t14 + (-m(6) * pkin(4) - m(7) * t6 + t44 * t19 + t45 * t21 - mrSges(5,1)) * t8 + t43 * (t16 * pkin(3) + t32) + (-m(6) * pkin(8) + m(7) * t17 + t49) * t7) * g(3) + (-t36 * t48 - mrSges(1,2) + t47 * r_base(2) + t45 * (-t34 * t7 + t36) - t46 * t33 - t44 * (t37 * t7 + t35) + t43 * (-t20 * t18 + t33) + t41 * t20 + (-m(5) * t30 - m(6) * (-t28 + t30) - m(7) * (-t27 + t30) + t46 * qJ(2) + t42) * t22) * g(2) + (-t37 * t48 - mrSges(1,1) + t47 * r_base(1) - t46 * t31 + t43 * (-t18 * t22 + t20 * t40 + t31) + t45 * (t35 * t7 + t37) - t44 * (-t36 * t7 + t34) + t41 * t22 + (-m(6) * t28 - m(7) * t27 - t42) * t20) * g(1);
U  = t1;
