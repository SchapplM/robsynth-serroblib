% Calculate potential energy for
% S5RRRRR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-05 19:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRRR6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR6_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:59:59
% EndTime: 2019-12-05 19:00:00
% DurationCPUTime: 0.30s
% Computational Cost: add. (153->58), mult. (109->44), div. (0->0), fcn. (73->10), ass. (0->24)
t33 = -m(1) - m(2);
t32 = -m(3) - m(4);
t31 = -m(5) - m(6) + t32;
t17 = sin(qJ(3));
t19 = cos(qJ(3));
t15 = qJ(3) + qJ(4);
t9 = qJ(5) + t15;
t2 = sin(t9);
t3 = cos(t9);
t4 = t19 * pkin(3) + pkin(2);
t5 = sin(t15);
t7 = cos(t15);
t30 = m(4) * pkin(2) + t19 * mrSges(4,1) - t17 * mrSges(4,2) + mrSges(3,1) + m(6) * (pkin(4) * t7 + t4) + t3 * mrSges(6,1) - t2 * mrSges(6,2) + m(5) * t4 + t7 * mrSges(5,1) - t5 * mrSges(5,2);
t21 = -pkin(8) - pkin(7);
t29 = -m(4) * pkin(7) + m(5) * t21 + m(6) * (-pkin(9) + t21) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t28 = pkin(5) + r_base(1);
t10 = pkin(6) + t28;
t25 = t17 * pkin(3) + t10;
t20 = cos(qJ(1));
t18 = sin(qJ(1));
t16 = qJ(1) + qJ(2);
t8 = cos(t16);
t6 = sin(t16);
t1 = (-mrSges(2,1) * t20 + mrSges(2,2) * t18 - mrSges(1,3) + t33 * r_base(3) - t30 * t8 + t31 * (t20 * pkin(1) + r_base(3)) + t29 * t6) * g(3) + (mrSges(2,1) * t18 + mrSges(2,2) * t20 - mrSges(1,2) + t33 * r_base(2) + t31 * (-t18 * pkin(1) + r_base(2)) + t29 * t8 + t30 * t6) * g(2) + (-m(1) * r_base(1) - mrSges(1,1) - m(2) * t28 - mrSges(2,3) - mrSges(3,3) - mrSges(4,1) * t17 - mrSges(4,2) * t19 - m(5) * t25 - t5 * mrSges(5,1) - t7 * mrSges(5,2) - m(6) * (pkin(4) * t5 + t25) - t2 * mrSges(6,1) - t3 * mrSges(6,2) + t32 * t10) * g(1);
U = t1;
