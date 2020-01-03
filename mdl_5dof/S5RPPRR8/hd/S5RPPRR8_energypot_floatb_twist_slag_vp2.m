% Calculate potential energy for
% S5RPPRR8
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR8_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPRR8_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR8_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR8_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR8_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:01
% EndTime: 2019-12-31 18:01:01
% DurationCPUTime: 0.31s
% Computational Cost: add. (138->63), mult. (138->52), div. (0->0), fcn. (120->8), ass. (0->29)
t40 = -m(1) - m(2);
t39 = -m(3) - m(4);
t38 = -m(5) - m(6);
t19 = sin(qJ(5));
t21 = cos(qJ(5));
t37 = m(6) * pkin(4) + t21 * mrSges(6,1) - t19 * mrSges(6,2) + mrSges(5,1);
t36 = mrSges(2,2) - mrSges(3,3);
t35 = -m(4) * pkin(2) - mrSges(2,1) - mrSges(3,1);
t34 = m(6) * pkin(7) - mrSges(5,2) + mrSges(6,3);
t17 = sin(pkin(8));
t20 = sin(qJ(1));
t33 = t20 * t17;
t22 = cos(qJ(1));
t32 = t22 * t17;
t16 = pkin(5) + r_base(3);
t31 = pkin(8) + qJ(4);
t30 = t20 * pkin(1) + r_base(2);
t29 = -qJ(3) + t16;
t28 = t22 * pkin(1) + t20 * qJ(2) + r_base(1);
t27 = cos(t31);
t26 = sin(t31);
t18 = cos(pkin(8));
t11 = pkin(3) * t18 + pkin(2);
t5 = t20 * t11;
t4 = t20 * t18 - t32;
t3 = -t22 * t18 - t33;
t2 = -t20 * t27 + t22 * t26;
t1 = -t20 * t26 - t22 * t27;
t6 = (-m(1) * r_base(3) - m(4) * t29 + t19 * mrSges(6,1) + t21 * mrSges(6,2) - mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + mrSges(4,3) + mrSges(5,3) + (-m(2) - m(3)) * t16 + t38 * (-pkin(6) + t29)) * g(3) + (-mrSges(1,2) - t4 * mrSges(4,1) - t3 * mrSges(4,2) - m(5) * (t5 + t30) - m(6) * (-pkin(3) * t32 + t5) + t40 * r_base(2) + (-m(6) + t39) * (-t22 * qJ(2) + t30) + (-m(5) * (-pkin(3) * t17 - qJ(2)) - t36) * t22 + t35 * t20 + t37 * t2 + t34 * t1) * g(2) + (t3 * mrSges(4,1) - t4 * mrSges(4,2) - mrSges(1,1) + t40 * r_base(1) + t39 * t28 + t38 * (pkin(3) * t33 + t22 * t11 + t28) + t35 * t22 + t36 * t20 - t34 * t2 + t37 * t1) * g(1);
U = t6;
