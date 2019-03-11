% Calculate potential energy for
% S6RRPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRR2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR2_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR2_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:16:34
% EndTime: 2019-03-09 13:16:34
% DurationCPUTime: 0.43s
% Computational Cost: add. (256->68), mult. (200->53), div. (0->0), fcn. (166->12), ass. (0->29)
t20 = qJ(5) + qJ(6);
t13 = sin(t20);
t14 = cos(t20);
t22 = sin(qJ(5));
t25 = cos(qJ(5));
t51 = -mrSges(5,1) - m(6) * pkin(4) - mrSges(6,1) * t25 + mrSges(6,2) * t22 - m(7) * (t25 * pkin(5) + pkin(4)) - mrSges(7,1) * t14 + mrSges(7,2) * t13;
t50 = mrSges(5,2) + m(7) * (-pkin(10) - pkin(9)) - mrSges(7,3) - m(6) * pkin(9) - mrSges(6,3);
t49 = -m(2) - m(3);
t48 = m(5) + m(6) + m(7);
t47 = -m(1) - m(4) + t49;
t18 = qJ(2) + pkin(11);
t10 = sin(t18);
t11 = cos(t18);
t23 = sin(qJ(2));
t26 = cos(qJ(2));
t12 = qJ(4) + t18;
t6 = sin(t12);
t7 = cos(t12);
t9 = t26 * pkin(2) + pkin(1);
t44 = -m(3) * pkin(1) - m(4) * t9 - t26 * mrSges(3,1) - t11 * mrSges(4,1) + t23 * mrSges(3,2) + t10 * mrSges(4,2) + t50 * t6 + t51 * t7 - mrSges(2,1);
t21 = -qJ(3) - pkin(7);
t43 = m(3) * pkin(7) - m(4) * t21 + t13 * mrSges(7,1) + t25 * mrSges(6,2) + t14 * mrSges(7,2) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3) + (m(7) * pkin(5) + mrSges(6,1)) * t22;
t19 = pkin(6) + r_base(3);
t39 = t23 * pkin(2) + t19;
t27 = cos(qJ(1));
t24 = sin(qJ(1));
t17 = -pkin(8) + t21;
t3 = pkin(3) * t11 + t9;
t1 = (-m(1) * r_base(3) - m(4) * t39 - t23 * mrSges(3,1) - t10 * mrSges(4,1) - t26 * mrSges(3,2) - t11 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) + t49 * t19 - t48 * (pkin(3) * t10 + t39) - t50 * t7 + t51 * t6) * g(3) + (-mrSges(1,2) - t48 * (t27 * t17 + t24 * t3 + r_base(2)) + t47 * r_base(2) + t43 * t27 + t44 * t24) * g(2) + (-mrSges(1,1) - t48 * (t27 * t3 + r_base(1)) + t47 * r_base(1) + t44 * t27 + (t17 * t48 - t43) * t24) * g(1);
U  = t1;
