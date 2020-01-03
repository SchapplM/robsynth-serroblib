% Calculate potential energy for
% S5RRPRP6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRP6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPRP6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP6_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP6_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:56:41
% EndTime: 2019-12-31 19:56:41
% DurationCPUTime: 0.39s
% Computational Cost: add. (159->60), mult. (166->52), div. (0->0), fcn. (138->8), ass. (0->27)
t20 = cos(qJ(4));
t43 = -m(5) * pkin(3) - m(6) * (t20 * pkin(4) + pkin(3)) - mrSges(4,1);
t42 = m(5) * pkin(7) - m(6) * (-qJ(5) - pkin(7)) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t41 = m(6) * pkin(4);
t40 = -m(2) - m(3);
t39 = -mrSges(5,1) - mrSges(6,1);
t38 = mrSges(5,2) + mrSges(6,2);
t37 = -m(1) + t40;
t36 = -m(4) - m(5) - m(6);
t13 = qJ(2) + pkin(8);
t10 = sin(t13);
t11 = cos(t13);
t18 = sin(qJ(2));
t21 = cos(qJ(2));
t35 = -m(3) * pkin(1) - t21 * mrSges(3,1) + t18 * mrSges(3,2) - t42 * t10 + t43 * t11 - mrSges(2,1);
t34 = m(3) * pkin(6) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3);
t17 = sin(qJ(4));
t19 = sin(qJ(1));
t33 = t19 * t17;
t32 = t19 * t20;
t22 = cos(qJ(1));
t31 = t22 * t17;
t30 = t22 * t20;
t14 = pkin(5) + r_base(3);
t16 = -qJ(3) - pkin(6);
t9 = t21 * pkin(2) + pkin(1);
t1 = (-m(1) * r_base(3) - t18 * mrSges(3,1) - t21 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + t40 * t14 + t36 * (t18 * pkin(2) + t14) + t42 * t11 + (t38 * t17 + t39 * t20 + t43) * t10) * g(3) + (t31 * t41 - mrSges(1,2) + t36 * (t22 * t16 + t19 * t9 + r_base(2)) + t39 * (t11 * t32 - t31) - t38 * (-t11 * t33 - t30) + t37 * r_base(2) + t34 * t22 + t35 * t19) * g(2) + (-t33 * t41 - mrSges(1,1) + t39 * (t11 * t30 + t33) - t38 * (-t11 * t31 + t32) + t36 * (-t19 * t16 + t22 * t9 + r_base(1)) + t37 * r_base(1) - t34 * t19 + t35 * t22) * g(1);
U = t1;
