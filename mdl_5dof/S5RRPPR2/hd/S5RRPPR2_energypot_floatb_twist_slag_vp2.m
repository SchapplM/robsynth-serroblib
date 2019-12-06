% Calculate potential energy for
% S5RRPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPPR2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPPR2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR2_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:20:00
% EndTime: 2019-12-05 18:20:00
% DurationCPUTime: 0.37s
% Computational Cost: add. (170->61), mult. (117->49), div. (0->0), fcn. (85->10), ass. (0->25)
t15 = sin(qJ(5));
t17 = cos(qJ(5));
t35 = -m(6) * pkin(4) - t17 * mrSges(6,1) + t15 * mrSges(6,2) - mrSges(5,1);
t34 = -m(1) - m(2);
t33 = -m(5) - m(6);
t32 = -t15 * mrSges(6,1) - t17 * mrSges(6,2) + mrSges(4,2) - mrSges(5,3);
t31 = m(6) * pkin(7) + mrSges(6,3);
t13 = sin(pkin(9));
t14 = cos(pkin(9));
t30 = t13 * mrSges(5,2) + t35 * t14 - mrSges(4,1);
t12 = qJ(1) + qJ(2);
t29 = pkin(5) + r_base(1);
t18 = cos(qJ(1));
t28 = t18 * pkin(1) + r_base(3);
t27 = pkin(6) + t29;
t10 = cos(t12);
t26 = pkin(2) * t10 + t28;
t16 = sin(qJ(1));
t25 = -pkin(1) * t16 + r_base(2);
t9 = sin(t12);
t21 = -pkin(2) * t9 + t25;
t8 = pkin(8) + t12;
t5 = cos(t8);
t4 = sin(t8);
t1 = (-m(3) * t28 - m(4) * t26 - mrSges(2,1) * t18 - t10 * mrSges(3,1) + t16 * mrSges(2,2) + t9 * mrSges(3,2) - mrSges(1,3) + t34 * r_base(3) + t33 * (t5 * pkin(3) + t4 * qJ(4) + t26) + (-t31 * t13 + t30) * t5 + t32 * t4) * g(3) + (-m(3) * t25 - m(4) * t21 + t16 * mrSges(2,1) + t9 * mrSges(3,1) + mrSges(2,2) * t18 + t10 * mrSges(3,2) - mrSges(1,2) + t34 * r_base(2) + t33 * (t5 * qJ(4) + t21) + t32 * t5 + (m(5) * pkin(3) - m(6) * (-pkin(7) * t13 - pkin(3)) + t13 * mrSges(6,3) - t30) * t4) * g(2) + (-m(1) * r_base(1) - m(2) * t29 - m(3) * t27 - mrSges(1,1) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) + (-m(4) + t33) * (qJ(3) + t27) + (-mrSges(5,2) + t31) * t14 + t35 * t13) * g(1);
U = t1;
