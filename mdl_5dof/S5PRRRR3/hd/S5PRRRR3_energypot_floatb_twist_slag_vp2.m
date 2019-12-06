% Calculate potential energy for
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% m_mdh [6x1]
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
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRR3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRRR3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR3_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR3_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:01
% EndTime: 2019-12-05 17:06:02
% DurationCPUTime: 0.23s
% Computational Cost: add. (152->57), mult. (86->44), div. (0->0), fcn. (50->10), ass. (0->25)
t34 = -m(1) - m(2);
t33 = -m(5) - m(6);
t19 = sin(qJ(5));
t20 = cos(qJ(5));
t32 = -m(6) * pkin(4) - mrSges(6,1) * t20 + mrSges(6,2) * t19 - mrSges(5,1);
t31 = m(6) * pkin(8) - mrSges(5,2) + mrSges(6,3);
t16 = pkin(9) + qJ(2);
t17 = sin(pkin(9));
t30 = t17 * pkin(1) + r_base(2);
t18 = cos(pkin(9));
t29 = t18 * pkin(1) + r_base(1);
t28 = qJ(1) + r_base(3);
t11 = sin(t16);
t27 = pkin(2) * t11 + t30;
t12 = cos(t16);
t26 = pkin(2) * t12 + t29;
t25 = pkin(5) + t28;
t13 = qJ(3) + t16;
t22 = pkin(6) + t25;
t10 = qJ(4) + t13;
t9 = cos(t13);
t8 = sin(t13);
t4 = cos(t10);
t3 = sin(t10);
t1 = (-m(1) * r_base(3) - m(2) * t28 - m(3) * t25 - m(4) * t22 - t19 * mrSges(6,1) - t20 * mrSges(6,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) + t33 * (pkin(7) + t22)) * g(3) + (-m(3) * t30 - m(4) * t27 - mrSges(2,1) * t17 - t11 * mrSges(3,1) - t8 * mrSges(4,1) - mrSges(2,2) * t18 - t12 * mrSges(3,2) - t9 * mrSges(4,2) - mrSges(1,2) + t34 * r_base(2) + t33 * (pkin(3) * t8 + t27) + t31 * t4 + t32 * t3) * g(2) + (-m(3) * t29 - m(4) * t26 - mrSges(2,1) * t18 - t12 * mrSges(3,1) - t9 * mrSges(4,1) + mrSges(2,2) * t17 + t11 * mrSges(3,2) + t8 * mrSges(4,2) - mrSges(1,1) + t34 * r_base(1) + t32 * t4 + t33 * (pkin(3) * t9 + t26) - t31 * t3) * g(1);
U = t1;
