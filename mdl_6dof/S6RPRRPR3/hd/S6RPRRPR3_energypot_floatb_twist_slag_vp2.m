% Calculate potential energy for
% S6RPRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPR3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRPR3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR3_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR3_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:04:30
% EndTime: 2019-03-09 05:04:31
% DurationCPUTime: 0.55s
% Computational Cost: add. (271->82), mult. (274->78), div. (0->0), fcn. (264->10), ass. (0->40)
t57 = -m(1) - m(2);
t56 = -m(6) - m(7);
t26 = sin(qJ(3));
t30 = cos(qJ(3));
t55 = -mrSges(4,1) * t30 + mrSges(4,2) * t26 - mrSges(3,1);
t54 = mrSges(3,2) - mrSges(4,3);
t25 = sin(qJ(4));
t29 = cos(qJ(4));
t53 = (pkin(4) * t29 + qJ(5) * t25) * t26;
t52 = mrSges(6,2) + mrSges(5,3) - mrSges(7,3);
t51 = m(7) * pkin(9) - t52;
t24 = sin(qJ(6));
t28 = cos(qJ(6));
t50 = -t24 * mrSges(7,1) - t28 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t49 = -m(7) * pkin(5) - t28 * mrSges(7,1) + t24 * mrSges(7,2) - mrSges(5,1) - mrSges(6,1);
t48 = pkin(3) * t30;
t23 = qJ(1) + pkin(10);
t17 = sin(t23);
t47 = t17 * t26;
t18 = cos(t23);
t46 = t18 * t26;
t45 = t25 * t30;
t44 = t29 * t30;
t43 = pkin(6) + r_base(3);
t27 = sin(qJ(1));
t42 = t27 * pkin(1) + r_base(2);
t31 = cos(qJ(1));
t41 = t31 * pkin(1) + r_base(1);
t19 = qJ(2) + t43;
t40 = t26 * pkin(3) + t19;
t39 = t18 * pkin(2) + t17 * pkin(7) + t41;
t37 = t17 * pkin(2) - pkin(7) * t18 + t42;
t36 = pkin(8) * t46 + t18 * t48 + t39;
t35 = -pkin(8) * t30 + t40;
t34 = pkin(8) * t47 + t17 * t48 + t37;
t6 = t17 * t25 + t18 * t44;
t5 = -t17 * t29 + t18 * t45;
t4 = t17 * t44 - t18 * t25;
t3 = t17 * t45 + t18 * t29;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t43 - mrSges(2,3) - mrSges(3,3) - m(5) * t35 - m(6) * (t35 + t53) - m(7) * (t40 + t53) + (-m(3) - m(4)) * t19 + (-mrSges(4,2) - m(7) * (-pkin(8) + pkin(9)) + t52) * t30 + (t50 * t25 + t49 * t29 - mrSges(4,1)) * t26) * g(3) + (-m(3) * t42 - m(4) * t37 - m(5) * t34 - t27 * mrSges(2,1) - t31 * mrSges(2,2) - mrSges(1,2) + t57 * r_base(2) + t56 * (t4 * pkin(4) + qJ(5) * t3 + t34) + t49 * t4 + t50 * t3 - t54 * t18 + t55 * t17 + t51 * t47) * g(2) + (-m(3) * t41 - m(4) * t39 - m(5) * t36 - t31 * mrSges(2,1) + t27 * mrSges(2,2) - mrSges(1,1) + t57 * r_base(1) + t56 * (t6 * pkin(4) + qJ(5) * t5 + t36) + t49 * t6 + t50 * t5 + t55 * t18 + t54 * t17 + t51 * t46) * g(1);
U  = t1;
