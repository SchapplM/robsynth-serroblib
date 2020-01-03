% Calculate potential energy for
% S5RPRRR6
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRR6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRRR6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR6_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR6_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:00:53
% EndTime: 2019-12-31 19:00:53
% DurationCPUTime: 0.33s
% Computational Cost: add. (167->57), mult. (129->45), div. (0->0), fcn. (97->10), ass. (0->26)
t15 = sin(qJ(5));
t18 = cos(qJ(5));
t40 = -m(6) * pkin(4) - t18 * mrSges(6,1) + t15 * mrSges(6,2) - mrSges(5,1);
t39 = -m(6) * pkin(8) + mrSges(5,2) - mrSges(6,3);
t38 = -m(1) - m(2);
t37 = -m(3) - m(4);
t36 = m(5) + m(6);
t16 = sin(qJ(3));
t19 = cos(qJ(3));
t14 = qJ(3) + qJ(4);
t8 = sin(t14);
t9 = cos(t14);
t34 = -m(4) * pkin(2) - t19 * mrSges(4,1) + t16 * mrSges(4,2) + t39 * t8 + t40 * t9 - mrSges(3,1);
t33 = m(4) * pkin(6) + t15 * mrSges(6,1) + t18 * mrSges(6,2) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t32 = pkin(5) + r_base(3);
t17 = sin(qJ(1));
t31 = t17 * pkin(1) + r_base(2);
t20 = cos(qJ(1));
t30 = t20 * pkin(1) + r_base(1);
t7 = qJ(2) + t32;
t21 = -pkin(7) - pkin(6);
t13 = qJ(1) + pkin(9);
t6 = cos(t13);
t5 = sin(t13);
t4 = pkin(3) * t19 + pkin(2);
t1 = (-m(1) * r_base(3) - m(2) * t32 - mrSges(4,1) * t16 - t19 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - t36 * (t16 * pkin(3) + t7) - t39 * t9 + t40 * t8 + t37 * t7) * g(3) + (-t17 * mrSges(2,1) - t20 * mrSges(2,2) - mrSges(1,2) + t38 * r_base(2) + t37 * t31 - t36 * (t6 * t21 + t5 * t4 + t31) + t33 * t6 + t34 * t5) * g(2) + (-mrSges(2,1) * t20 + t17 * mrSges(2,2) - mrSges(1,1) + t38 * r_base(1) + t37 * t30 - t36 * (t6 * t4 + t30) + t34 * t6 + (t36 * t21 - t33) * t5) * g(1);
U = t1;
