% Calculate potential energy for
% S5PRRPR8
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
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRPR8_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRPR8_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR8_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR8_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR8_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:14
% EndTime: 2019-12-31 17:42:15
% DurationCPUTime: 0.33s
% Computational Cost: add. (168->57), mult. (142->45), div. (0->0), fcn. (110->10), ass. (0->26)
t19 = sin(qJ(5));
t21 = cos(qJ(5));
t41 = -m(6) * pkin(4) - t21 * mrSges(6,1) + t19 * mrSges(6,2) - mrSges(5,1);
t40 = -m(6) * pkin(7) + mrSges(5,2) - mrSges(6,3);
t39 = -m(2) - m(3);
t38 = m(5) + m(6);
t37 = -m(1) - m(4) + t39;
t16 = qJ(2) + qJ(3);
t10 = sin(t16);
t11 = cos(t16);
t20 = sin(qJ(2));
t22 = cos(qJ(2));
t9 = pkin(9) + t16;
t5 = sin(t9);
t6 = cos(t9);
t8 = t22 * pkin(2) + pkin(1);
t35 = -m(3) * pkin(1) - m(4) * t8 - t22 * mrSges(3,1) - t11 * mrSges(4,1) + t20 * mrSges(3,2) + t10 * mrSges(4,2) + t40 * t5 + t41 * t6 - mrSges(2,1);
t23 = -pkin(6) - pkin(5);
t34 = m(3) * pkin(5) - m(4) * t23 + t19 * mrSges(6,1) + t21 * mrSges(6,2) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3);
t14 = qJ(1) + r_base(3);
t31 = t20 * pkin(2) + t14;
t18 = cos(pkin(8));
t17 = sin(pkin(8));
t15 = -qJ(4) + t23;
t3 = pkin(3) * t11 + t8;
t1 = (-m(1) * r_base(3) - m(4) * t31 - t20 * mrSges(3,1) - t10 * mrSges(4,1) - t22 * mrSges(3,2) - t11 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - t38 * (pkin(3) * t10 + t31) - t40 * t6 + t41 * t5 + t39 * t14) * g(3) + (-mrSges(1,2) - t38 * (t18 * t15 + t17 * t3 + r_base(2)) + t37 * r_base(2) + t34 * t18 + t35 * t17) * g(2) + (-mrSges(1,1) - t38 * (t18 * t3 + r_base(1)) + t37 * r_base(1) + t35 * t18 + (t38 * t15 - t34) * t17) * g(1);
U = t1;
