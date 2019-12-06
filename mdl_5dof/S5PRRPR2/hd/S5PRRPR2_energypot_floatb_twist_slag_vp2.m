% Calculate potential energy for
% S5PRRPR2
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRPR2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRPR2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR2_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR2_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:15
% EndTime: 2019-12-05 16:17:15
% DurationCPUTime: 0.32s
% Computational Cost: add. (170->57), mult. (117->45), div. (0->0), fcn. (85->10), ass. (0->25)
t19 = sin(qJ(5));
t20 = cos(qJ(5));
t39 = -m(6) * pkin(4) - t20 * mrSges(6,1) + t19 * mrSges(6,2) - mrSges(5,1);
t38 = -m(6) * pkin(7) + mrSges(5,2) - mrSges(6,3);
t37 = -m(1) - m(2);
t36 = m(5) + m(6);
t15 = sin(pkin(9));
t17 = cos(pkin(9));
t35 = t38 * t15 + t17 * t39 - mrSges(4,1);
t34 = -t19 * mrSges(6,1) - t20 * mrSges(6,2) + mrSges(4,2) - mrSges(5,3);
t14 = pkin(8) + qJ(2);
t16 = sin(pkin(8));
t32 = t16 * pkin(1) + r_base(2);
t18 = cos(pkin(8));
t31 = t18 * pkin(1) + r_base(1);
t30 = qJ(1) + r_base(3);
t9 = sin(t14);
t29 = pkin(2) * t9 + t32;
t10 = cos(t14);
t28 = pkin(2) * t10 + t31;
t27 = pkin(5) + t30;
t11 = qJ(3) + t14;
t7 = cos(t11);
t6 = sin(t11);
t1 = (-m(1) * r_base(3) - m(2) * t30 - m(3) * t27 - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) + (-m(4) - t36) * (pkin(6) + t27) - t38 * t17 + t39 * t15) * g(3) + (-m(3) * t32 - m(4) * t29 - mrSges(2,1) * t16 - t9 * mrSges(3,1) - mrSges(2,2) * t18 - t10 * mrSges(3,2) - mrSges(1,2) + t37 * r_base(2) - t36 * (t6 * pkin(3) + t29) + (t36 * qJ(4) - t34) * t7 + t35 * t6) * g(2) + (-m(3) * t31 - m(4) * t28 - mrSges(2,1) * t18 - t10 * mrSges(3,1) + mrSges(2,2) * t16 + t9 * mrSges(3,2) - mrSges(1,1) + t37 * r_base(1) - t36 * (t7 * pkin(3) + t6 * qJ(4) + t28) + t35 * t7 + t34 * t6) * g(1);
U = t1;
