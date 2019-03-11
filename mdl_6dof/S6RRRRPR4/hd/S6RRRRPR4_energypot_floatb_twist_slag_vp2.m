% Calculate potential energy for
% S6RRRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPR4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRPR4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR4_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR4_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:06:00
% EndTime: 2019-03-09 22:06:00
% DurationCPUTime: 0.45s
% Computational Cost: add. (257->68), mult. (226->54), div. (0->0), fcn. (196->12), ass. (0->29)
t18 = qJ(4) + pkin(11);
t10 = sin(t18);
t11 = cos(t18);
t22 = sin(qJ(4));
t25 = cos(qJ(4));
t12 = qJ(6) + t18;
t5 = sin(t12);
t6 = cos(t12);
t7 = t25 * pkin(4) + pkin(3);
t53 = -m(5) * pkin(3) - m(6) * t7 - t25 * mrSges(5,1) - t11 * mrSges(6,1) + t22 * mrSges(5,2) + t10 * mrSges(6,2) - mrSges(4,1) - m(7) * (pkin(5) * t11 + t7) - t6 * mrSges(7,1) + t5 * mrSges(7,2);
t21 = -qJ(5) - pkin(9);
t52 = mrSges(4,2) + m(7) * (-pkin(10) + t21) - mrSges(7,3) + m(6) * t21 - mrSges(6,3) - m(5) * pkin(9) - mrSges(5,3);
t51 = -m(2) - m(3);
t50 = -m(1) + t51;
t49 = m(4) + m(5) + m(6) + m(7);
t20 = qJ(2) + qJ(3);
t13 = sin(t20);
t14 = cos(t20);
t23 = sin(qJ(2));
t26 = cos(qJ(2));
t45 = -m(3) * pkin(1) - t26 * mrSges(3,1) + t23 * mrSges(3,2) + t52 * t13 + t53 * t14 - mrSges(2,1);
t43 = pkin(4) * t22;
t44 = m(3) * pkin(7) + m(6) * t43 + m(7) * (pkin(5) * t10 + t43) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + t10 * mrSges(6,1) + t11 * mrSges(6,2) + t22 * mrSges(5,1) + t25 * mrSges(5,2) + t5 * mrSges(7,1) + t6 * mrSges(7,2);
t19 = pkin(6) + r_base(3);
t28 = -pkin(8) - pkin(7);
t27 = cos(qJ(1));
t24 = sin(qJ(1));
t8 = pkin(2) * t26 + pkin(1);
t1 = (-m(1) * r_base(3) - t23 * mrSges(3,1) - t26 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + t51 * t19 - t49 * (t23 * pkin(2) + t19) - t52 * t14 + t53 * t13) * g(3) + (-mrSges(1,2) + t50 * r_base(2) - t49 * (t24 * t8 + t27 * t28 + r_base(2)) + t44 * t27 + t45 * t24) * g(2) + (-mrSges(1,1) + t50 * r_base(1) - t49 * (t27 * t8 + r_base(1)) + t45 * t27 + (t49 * t28 - t44) * t24) * g(1);
U  = t1;
