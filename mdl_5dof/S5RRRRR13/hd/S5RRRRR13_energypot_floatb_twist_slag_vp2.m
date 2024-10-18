% Calculate potential energy for
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR13_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRRR13_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR13_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR13_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR13_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:30:17
% EndTime: 2024-09-27 17:30:17
% DurationCPUTime: 0.08s
% Computational Cost: add. (213->66), mult. (140->55), div. (0->0), fcn. (106->16), ass. (0->30)
t26 = sin(qJ(4));
t28 = cos(qJ(4));
t51 = -t28 * mrSges(5,2) + (-m(6) * pkin(4) - mrSges(5,1)) * t26;
t49 = m(6) * (pkin(9) + pkin(10)) + mrSges(5,3) + mrSges(6,3) + m(5) * pkin(9);
t46 = -m(1) - m(2);
t45 = -m(4) - m(5) - m(6);
t22 = qJ(4) + qJ(5);
t44 = -m(5) * pkin(3) - m(6) * (t28 * pkin(4) + pkin(3)) - t28 * mrSges(5,1) - cos(t22) * mrSges(6,1) + t26 * mrSges(5,2) + sin(t22) * mrSges(6,2) - mrSges(4,1);
t24 = sin(pkin(5));
t25 = cos(pkin(5));
t13 = pkin(5) + t22;
t4 = sin(t13) / 0.2e1;
t14 = pkin(5) - t22;
t5 = cos(t14) / 0.2e1;
t6 = sin(t14);
t7 = cos(t13);
t43 = -(t4 - t6 / 0.2e1) * mrSges(6,1) - (t5 + t7 / 0.2e1) * mrSges(6,2) - mrSges(4,2) + t49 * t24 + t51 * t25;
t23 = qJ(1) + qJ(2);
t37 = pkin(6) + r_base(3);
t27 = sin(qJ(1));
t36 = t27 * pkin(1) + r_base(2);
t29 = cos(qJ(1));
t35 = t29 * pkin(1) + r_base(1);
t34 = pkin(7) + t37;
t19 = qJ(3) + t23;
t18 = cos(t23);
t16 = sin(t23);
t11 = cos(t19);
t10 = sin(t19);
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t37 - mrSges(2,3) - m(3) * t34 - mrSges(3,3) - mrSges(4,3) - (t5 - t7 / 0.2e1) * mrSges(6,1) - (t4 + t6 / 0.2e1) * mrSges(6,2) + t51 * t24 + t45 * (pkin(8) + t34) - t49 * t25) * g(3) + (-m(3) * t36 - t27 * mrSges(2,1) - t16 * mrSges(3,1) - t29 * mrSges(2,2) - t18 * mrSges(3,2) - mrSges(1,2) + t46 * r_base(2) + t45 * (pkin(2) * t16 + t36) + t44 * t10 + t43 * t11) * g(2) + (-m(3) * t35 - t29 * mrSges(2,1) - t18 * mrSges(3,1) + t27 * mrSges(2,2) + t16 * mrSges(3,2) - mrSges(1,1) + t46 * r_base(1) + t45 * (pkin(2) * t18 + t35) + t44 * t11 - t43 * t10) * g(1);
U = t1;
