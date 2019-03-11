% Calculate potential energy for
% S6RPPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
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
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPPRR1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPPRR1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR1_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR1_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:29:47
% EndTime: 2019-03-09 01:29:47
% DurationCPUTime: 0.34s
% Computational Cost: add. (180->67), mult. (141->47), div. (0->0), fcn. (103->8), ass. (0->27)
t13 = sin(qJ(6));
t16 = cos(qJ(6));
t39 = -m(7) * pkin(5) - t16 * mrSges(7,1) + t13 * mrSges(7,2) - mrSges(6,1);
t38 = -m(7) * pkin(8) + mrSges(6,2) - mrSges(7,3);
t37 = -m(1) - m(2);
t36 = m(6) + m(7);
t14 = sin(qJ(5));
t17 = cos(qJ(5));
t34 = t39 * t14 - t38 * t17 - mrSges(3,1) + mrSges(4,2) - mrSges(5,3);
t33 = t13 * mrSges(7,1) + t16 * mrSges(7,2) + mrSges(3,2) - mrSges(5,2) - mrSges(4,3) + mrSges(6,3);
t32 = pkin(6) + r_base(3);
t15 = sin(qJ(1));
t31 = t15 * pkin(1) + r_base(2);
t18 = cos(qJ(1));
t30 = t18 * pkin(1) + r_base(1);
t12 = qJ(1) + pkin(9);
t7 = sin(t12);
t29 = t7 * pkin(2) + t31;
t9 = qJ(2) + t32;
t8 = cos(t12);
t28 = t8 * pkin(2) + t7 * qJ(3) + t30;
t27 = pkin(3) + t9;
t22 = -qJ(3) * t8 + t29;
t1 = t7 * qJ(4);
t21 = t1 + t22;
t5 = t8 * pkin(7);
t2 = (-m(1) * r_base(3) - m(2) * t32 - m(5) * t27 - mrSges(4,1) - mrSges(5,1) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + (-m(3) - m(4)) * t9 - t36 * (pkin(4) + t27) + t39 * t17 + t38 * t14) * g(3) + (-mrSges(1,2) - t15 * mrSges(2,1) - mrSges(2,2) * t18 - m(3) * t31 - m(4) * t22 - m(5) * t21 - m(6) * (t21 + t5) - m(7) * (t1 + t5 + t29) + t37 * r_base(2) + (m(7) * qJ(3) - t33) * t8 + t34 * t7) * g(2) + (-m(3) * t30 - m(4) * t28 - mrSges(2,1) * t18 + t15 * mrSges(2,2) - mrSges(1,1) + t37 * r_base(1) + (-m(5) - t36) * (t8 * qJ(4) + t28) + t34 * t8 + (t36 * pkin(7) + t33) * t7) * g(1);
U  = t2;
