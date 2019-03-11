% Calculate potential energy for
% S6RPPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRR4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPRRR4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR4_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR4_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:25:17
% EndTime: 2019-03-09 02:25:18
% DurationCPUTime: 0.46s
% Computational Cost: add. (204->66), mult. (287->53), div. (0->0), fcn. (305->10), ass. (0->27)
t52 = m(5) + m(6) + m(7);
t18 = qJ(5) + qJ(6);
t10 = sin(t18);
t11 = cos(t18);
t19 = sin(qJ(5));
t21 = cos(qJ(5));
t51 = mrSges(5,1) + m(6) * pkin(4) + mrSges(6,1) * t21 - mrSges(6,2) * t19 + m(7) * (pkin(5) * t21 + pkin(4)) + mrSges(7,1) * t11 - mrSges(7,2) * t10;
t50 = mrSges(5,2) + m(7) * (-pkin(9) - pkin(8)) - mrSges(7,3) - m(6) * pkin(8) - mrSges(6,3);
t49 = -m(1) - m(2);
t48 = -mrSges(2,1) - mrSges(3,1);
t47 = mrSges(2,2) - mrSges(3,3);
t20 = sin(qJ(4));
t22 = cos(qJ(4));
t43 = -t50 * t20 + t51 * t22 + mrSges(4,1);
t42 = t10 * mrSges(7,1) + t21 * mrSges(6,2) + t11 * mrSges(7,2) - mrSges(4,2) + mrSges(5,3) + (m(7) * pkin(5) + mrSges(6,1)) * t19 + t52 * pkin(7);
t41 = cos(qJ(1));
t40 = sin(qJ(1));
t39 = cos(pkin(10));
t38 = sin(pkin(10));
t17 = pkin(6) + r_base(3);
t36 = t41 * pkin(1) + t40 * qJ(2) + r_base(1);
t35 = t41 * pkin(2) + t36;
t31 = t40 * pkin(1) - t41 * qJ(2) + r_base(2);
t27 = t40 * pkin(2) + t31;
t4 = t41 * t38 - t40 * t39;
t3 = -t40 * t38 - t41 * t39;
t1 = (-m(1) * r_base(3) - mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + mrSges(4,3) + (-m(2) - m(3)) * t17 + (-m(4) - t52) * (-qJ(3) + t17) + t50 * t22 + t51 * t20) * g(3) + (-m(3) * t31 - m(4) * t27 - mrSges(1,2) + t49 * r_base(2) - t47 * t41 + t48 * t40 - t52 * (-t4 * pkin(3) + t27) + t43 * t4 + t42 * t3) * g(2) + (-m(3) * t36 - m(4) * t35 - mrSges(1,1) + t49 * r_base(1) + t48 * t41 + t47 * t40 - t52 * (-t3 * pkin(3) + t35) - t42 * t4 + t43 * t3) * g(1);
U  = t1;
