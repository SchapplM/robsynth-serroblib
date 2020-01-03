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
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:57:17
% EndTime: 2020-01-03 11:57:18
% DurationCPUTime: 0.40s
% Computational Cost: add. (170->59), mult. (117->47), div. (0->0), fcn. (85->10), ass. (0->25)
t34 = m(5) + m(6);
t13 = sin(qJ(5));
t15 = cos(qJ(5));
t33 = -m(6) * pkin(4) - t15 * mrSges(6,1) + t13 * mrSges(6,2) - mrSges(5,1);
t32 = -m(1) - m(2);
t30 = -m(4) - t34;
t29 = m(6) * pkin(7) + mrSges(6,3);
t11 = sin(pkin(9));
t12 = cos(pkin(9));
t28 = t11 * mrSges(5,2) + t33 * t12 - mrSges(4,1);
t27 = t13 * mrSges(6,1) + t15 * mrSges(6,2) + t34 * qJ(4) - mrSges(4,2) + mrSges(5,3);
t10 = qJ(1) + qJ(2);
t26 = pkin(5) + r_base(1);
t14 = sin(qJ(1));
t25 = t14 * pkin(1) + r_base(2);
t24 = pkin(6) + t26;
t7 = sin(t10);
t23 = pkin(2) * t7 + t25;
t16 = cos(qJ(1));
t21 = -pkin(1) * t16 + r_base(3);
t8 = cos(t10);
t6 = pkin(8) + t10;
t3 = cos(t6);
t2 = sin(t6);
t1 = (-m(3) * t21 + mrSges(2,1) * t16 + t8 * mrSges(3,1) - t14 * mrSges(2,2) - t7 * mrSges(3,2) - mrSges(1,3) + t32 * r_base(3) + (m(5) * pkin(3) - m(6) * (-pkin(7) * t11 - pkin(3)) + t11 * mrSges(6,3) - t28) * t3 + t30 * (-pkin(2) * t8 + t21) + t27 * t2) * g(3) + (-m(3) * t25 - m(4) * t23 - t14 * mrSges(2,1) - t7 * mrSges(3,1) - mrSges(2,2) * t16 - t8 * mrSges(3,2) - mrSges(1,2) + t32 * r_base(2) - t34 * (t2 * pkin(3) + t23) + t27 * t3 + (-t29 * t11 + t28) * t2) * g(2) + (-m(1) * r_base(1) - m(2) * t26 - m(3) * t24 - mrSges(1,1) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) + t30 * (qJ(3) + t24) + (-mrSges(5,2) + t29) * t12 + t33 * t11) * g(1);
U = t1;
