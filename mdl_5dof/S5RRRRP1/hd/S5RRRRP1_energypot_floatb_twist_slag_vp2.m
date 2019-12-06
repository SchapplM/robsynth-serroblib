% Calculate potential energy for
% S5RRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRP1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRRP1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP1_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP1_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:44:40
% EndTime: 2019-12-05 18:44:41
% DurationCPUTime: 0.28s
% Computational Cost: add. (150->52), mult. (122->37), div. (0->0), fcn. (86->8), ass. (0->23)
t31 = -m(5) - m(6);
t34 = mrSges(5,2) + mrSges(6,2);
t33 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t32 = -m(2) - m(3);
t17 = sin(qJ(2));
t19 = cos(qJ(2));
t16 = qJ(2) + qJ(3);
t11 = qJ(4) + t16;
t5 = sin(t11);
t6 = cos(t11);
t7 = t19 * pkin(2) + pkin(1);
t8 = sin(t16);
t9 = cos(t16);
t30 = -m(3) * pkin(1) - m(4) * t7 - t19 * mrSges(3,1) - t9 * mrSges(4,1) + t17 * mrSges(3,2) + t8 * mrSges(4,2) - mrSges(2,1) + t31 * (pkin(3) * t9 + t7) + t33 * t6 + t34 * t5;
t29 = -m(1) - m(4) + t31 + t32;
t21 = -pkin(7) - pkin(6);
t15 = -pkin(8) + t21;
t28 = m(3) * pkin(6) - m(4) * t21 - m(5) * t15 - m(6) * (-qJ(5) + t15) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t14 = pkin(5) + r_base(3);
t27 = t17 * pkin(2) + t14;
t20 = cos(qJ(1));
t18 = sin(qJ(1));
t1 = (-m(1) * r_base(3) - m(4) * t27 - mrSges(3,1) * t17 - t8 * mrSges(4,1) - mrSges(3,2) * t19 - t9 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - t34 * t6 + t31 * (pkin(3) * t8 + t27) + t33 * t5 + t32 * t14) * g(3) + (t18 * t30 + t20 * t28 + t29 * r_base(2) - mrSges(1,2)) * g(2) + (-t18 * t28 + t20 * t30 + t29 * r_base(1) - mrSges(1,1)) * g(1);
U = t1;
