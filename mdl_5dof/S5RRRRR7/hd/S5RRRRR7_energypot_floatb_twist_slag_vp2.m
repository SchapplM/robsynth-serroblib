% Calculate potential energy for
% S5RRRRR7
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
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR7_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRRR7_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR7_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR7_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR7_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:21:18
% EndTime: 2019-12-31 22:21:18
% DurationCPUTime: 0.33s
% Computational Cost: add. (168->57), mult. (142->45), div. (0->0), fcn. (110->10), ass. (0->26)
t17 = sin(qJ(5));
t20 = cos(qJ(5));
t41 = -m(6) * pkin(4) - t20 * mrSges(6,1) + t17 * mrSges(6,2) - mrSges(5,1);
t40 = -m(6) * pkin(9) + mrSges(5,2) - mrSges(6,3);
t39 = -m(2) - m(3);
t38 = m(5) + m(6);
t37 = -m(1) - m(4) + t39;
t16 = qJ(2) + qJ(3);
t10 = cos(t16);
t18 = sin(qJ(2));
t21 = cos(qJ(2));
t11 = qJ(4) + t16;
t6 = sin(t11);
t7 = cos(t11);
t8 = t21 * pkin(2) + pkin(1);
t9 = sin(t16);
t35 = -m(3) * pkin(1) - m(4) * t8 - t21 * mrSges(3,1) - t10 * mrSges(4,1) + t18 * mrSges(3,2) + t9 * mrSges(4,2) + t40 * t6 + t41 * t7 - mrSges(2,1);
t23 = -pkin(7) - pkin(6);
t34 = m(3) * pkin(6) - m(4) * t23 + t17 * mrSges(6,1) + t20 * mrSges(6,2) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3);
t14 = pkin(5) + r_base(3);
t31 = t18 * pkin(2) + t14;
t22 = cos(qJ(1));
t19 = sin(qJ(1));
t15 = -pkin(8) + t23;
t3 = pkin(3) * t10 + t8;
t1 = (-m(1) * r_base(3) - m(4) * t31 - t18 * mrSges(3,1) - t9 * mrSges(4,1) - t21 * mrSges(3,2) - t10 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - t38 * (pkin(3) * t9 + t31) - t40 * t7 + t41 * t6 + t39 * t14) * g(3) + (-mrSges(1,2) - t38 * (t22 * t15 + t19 * t3 + r_base(2)) + t37 * r_base(2) + t34 * t22 + t35 * t19) * g(2) + (-mrSges(1,1) - t38 * (t22 * t3 + r_base(1)) + t37 * r_base(1) + t35 * t22 + (t38 * t15 - t34) * t19) * g(1);
U = t1;
