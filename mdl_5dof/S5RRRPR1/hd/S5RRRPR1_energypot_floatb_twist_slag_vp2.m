% Calculate potential energy for
% S5RRRPR1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPR1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRPR1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR1_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR1_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:37:43
% EndTime: 2019-12-05 18:37:44
% DurationCPUTime: 0.27s
% Computational Cost: add. (156->57), mult. (122->44), div. (0->0), fcn. (86->10), ass. (0->25)
t34 = -m(2) - m(3);
t22 = cos(qJ(2));
t10 = t22 * pkin(2) + pkin(1);
t19 = qJ(2) + qJ(3);
t12 = sin(t19);
t13 = cos(t19);
t2 = pkin(3) * t13 + t10;
t20 = sin(qJ(2));
t11 = pkin(9) + t19;
t9 = qJ(5) + t11;
t3 = sin(t9);
t4 = cos(t9);
t5 = sin(t11);
t6 = cos(t11);
t33 = -m(3) * pkin(1) - m(4) * t10 - t22 * mrSges(3,1) - t13 * mrSges(4,1) + t20 * mrSges(3,2) + t12 * mrSges(4,2) - mrSges(2,1) - m(6) * (pkin(4) * t6 + t2) - t4 * mrSges(6,1) + t3 * mrSges(6,2) - m(5) * t2 - t6 * mrSges(5,1) + t5 * mrSges(5,2);
t32 = -m(1) - m(4) - m(5) - m(6) + t34;
t24 = -pkin(7) - pkin(6);
t17 = -qJ(4) + t24;
t31 = m(3) * pkin(6) - m(4) * t24 - m(5) * t17 - m(6) * (-pkin(8) + t17) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t18 = pkin(5) + r_base(3);
t30 = t20 * pkin(2) + t18;
t29 = pkin(3) * t12 + t30;
t23 = cos(qJ(1));
t21 = sin(qJ(1));
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - mrSges(3,1) * t20 - mrSges(3,2) * t22 - m(4) * t30 - t12 * mrSges(4,1) - t13 * mrSges(4,2) - m(5) * t29 - t5 * mrSges(5,1) - t6 * mrSges(5,2) - m(6) * (pkin(4) * t5 + t29) - t3 * mrSges(6,1) - t4 * mrSges(6,2) + t34 * t18) * g(3) + (t33 * t21 + t31 * t23 + t32 * r_base(2) - mrSges(1,2)) * g(2) + (-t31 * t21 + t33 * t23 + t32 * r_base(1) - mrSges(1,1)) * g(1);
U = t1;
