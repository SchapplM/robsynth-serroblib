% Calculate potential energy for
% S6RPRRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP8_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRRP8_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP8_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP8_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP8_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:22:45
% EndTime: 2019-03-09 06:22:46
% DurationCPUTime: 0.46s
% Computational Cost: add. (192->77), mult. (207->66), div. (0->0), fcn. (177->8), ass. (0->34)
t53 = -m(1) - m(2);
t52 = m(3) + m(4);
t51 = -m(6) - m(7);
t50 = mrSges(6,3) + mrSges(7,2);
t18 = qJ(3) + qJ(4);
t11 = sin(t18);
t12 = cos(t18);
t20 = sin(qJ(3));
t23 = cos(qJ(3));
t49 = -t20 * mrSges(4,1) - t11 * mrSges(5,1) - t23 * mrSges(4,2) - t12 * mrSges(5,2) + mrSges(2,2) - mrSges(3,3);
t48 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t47 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t46 = -m(4) * pkin(7) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t45 = pkin(3) * t20;
t44 = pkin(4) * t11;
t21 = sin(qJ(1));
t43 = t12 * t21;
t19 = sin(qJ(5));
t42 = t19 * t21;
t24 = cos(qJ(1));
t41 = t19 * t24;
t22 = cos(qJ(5));
t40 = t21 * t22;
t39 = t24 * t22;
t17 = pkin(6) + r_base(3);
t38 = t21 * pkin(1) + r_base(2);
t37 = pkin(2) + t17;
t36 = -qJ(2) - t45;
t35 = t24 * pkin(1) + t21 * qJ(2) + r_base(1);
t34 = t23 * pkin(3) + t37;
t25 = -pkin(8) - pkin(7);
t31 = -t21 * t25 + t38;
t27 = t21 * t45 - t24 * t25 + t35;
t1 = (-m(1) * r_base(3) - m(4) * t37 - m(5) * t34 - t23 * mrSges(4,1) + t20 * mrSges(4,2) - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + t51 * (t12 * pkin(4) + t11 * pkin(9) + t34) + (-m(2) - m(3)) * t17 + (-t47 * t19 + t48 * t22 - mrSges(5,1)) * t12 + (mrSges(5,2) - t50) * t11) * g(3) + (-m(5) * t31 - mrSges(1,2) + t53 * r_base(2) + t51 * (t24 * t12 * pkin(9) + t31) + t48 * (-t11 * t39 + t42) - t52 * t38 + t47 * (t11 * t41 + t40) + t46 * t21 + (-m(5) * t36 + t51 * (t36 - t44) - t50 * t12 + t52 * qJ(2) - t49) * t24) * g(2) + (-m(5) * t27 - mrSges(1,1) + t53 * r_base(1) + t50 * t43 - t52 * t35 + t51 * (-pkin(9) * t43 + t21 * t44 + t27) + t48 * (t11 * t40 + t41) - t47 * (t11 * t42 - t39) + t46 * t24 + t49 * t21) * g(1);
U  = t1;
