% Calculate potential energy for
% S5RPRRR10
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
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRR10_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRRR10_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR10_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR10_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR10_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:09:38
% EndTime: 2019-12-31 19:09:39
% DurationCPUTime: 0.36s
% Computational Cost: add. (169->56), mult. (166->44), div. (0->0), fcn. (138->10), ass. (0->23)
t17 = sin(qJ(4));
t19 = cos(qJ(4));
t13 = qJ(4) + qJ(5);
t8 = sin(t13);
t9 = cos(t13);
t42 = -m(5) * pkin(3) - t19 * mrSges(5,1) + t17 * mrSges(5,2) - mrSges(4,1) - m(6) * (t19 * pkin(4) + pkin(3)) - t9 * mrSges(6,1) + t8 * mrSges(6,2);
t41 = mrSges(4,2) + m(6) * (-pkin(8) - pkin(7)) - mrSges(6,3) - m(5) * pkin(7) - mrSges(5,3);
t40 = -m(2) - m(3);
t39 = -m(1) + t40;
t38 = m(4) + m(5) + m(6);
t14 = sin(pkin(9));
t15 = cos(pkin(9));
t11 = pkin(9) + qJ(3);
t6 = sin(t11);
t7 = cos(t11);
t35 = -m(3) * pkin(1) - t15 * mrSges(3,1) + t14 * mrSges(3,2) + t41 * t6 + t42 * t7 - mrSges(2,1);
t34 = m(3) * qJ(2) + t8 * mrSges(6,1) + t19 * mrSges(5,2) + t9 * mrSges(6,2) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + (m(6) * pkin(4) + mrSges(5,1)) * t17;
t12 = pkin(5) + r_base(3);
t20 = cos(qJ(1));
t18 = sin(qJ(1));
t16 = -pkin(6) - qJ(2);
t3 = t15 * pkin(2) + pkin(1);
t1 = (-m(1) * r_base(3) - t14 * mrSges(3,1) - t15 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + t40 * t12 - t38 * (t14 * pkin(2) + t12) - t41 * t7 + t42 * t6) * g(3) + (-mrSges(1,2) + t39 * r_base(2) - t38 * (t20 * t16 + t18 * t3 + r_base(2)) + t34 * t20 + t35 * t18) * g(2) + (-mrSges(1,1) + t39 * r_base(1) - t38 * (t20 * t3 + r_base(1)) + t35 * t20 + (t38 * t16 - t34) * t18) * g(1);
U = t1;
