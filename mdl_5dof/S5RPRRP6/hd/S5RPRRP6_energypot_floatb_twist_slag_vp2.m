% Calculate potential energy for
% S5RPRRP6
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRRP6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP6_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP6_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:42:05
% EndTime: 2019-12-31 18:42:05
% DurationCPUTime: 0.38s
% Computational Cost: add. (165->60), mult. (153->54), div. (0->0), fcn. (125->8), ass. (0->26)
t19 = cos(qJ(4));
t42 = -m(5) * pkin(3) - m(6) * (t19 * pkin(4) + pkin(3)) - mrSges(4,1);
t41 = m(5) * pkin(7) - m(6) * (-qJ(5) - pkin(7)) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t40 = m(6) * pkin(4);
t39 = -m(1) - m(2);
t38 = -mrSges(5,1) - mrSges(6,1);
t37 = mrSges(3,2) - mrSges(4,3);
t36 = mrSges(5,2) + mrSges(6,2);
t35 = -m(4) - m(5) - m(6);
t17 = sin(qJ(3));
t20 = cos(qJ(3));
t34 = -t41 * t17 + t42 * t20 - mrSges(3,1);
t16 = sin(qJ(4));
t14 = qJ(1) + pkin(8);
t9 = sin(t14);
t33 = t9 * t16;
t10 = cos(t14);
t32 = t10 * t16;
t31 = t16 * t20;
t30 = t19 * t20;
t29 = pkin(5) + r_base(3);
t18 = sin(qJ(1));
t28 = t18 * pkin(1) + r_base(2);
t21 = cos(qJ(1));
t27 = t21 * pkin(1) + r_base(1);
t1 = (-m(1) * r_base(3) - m(2) * t29 - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + (-m(3) + t35) * (qJ(2) + t29) + t41 * t20 + (t36 * t16 + t38 * t19 + t42) * t17) * g(3) + (t32 * t40 - m(3) * t28 - t18 * mrSges(2,1) - t21 * mrSges(2,2) - mrSges(1,2) + t39 * r_base(2) + t35 * (t9 * pkin(2) - t10 * pkin(6) + t28) + t38 * (t30 * t9 - t32) - t37 * t10 - t36 * (-t10 * t19 - t31 * t9) + t34 * t9) * g(2) + (-t33 * t40 - m(3) * t27 - t21 * mrSges(2,1) + t18 * mrSges(2,2) - mrSges(1,1) + t39 * r_base(1) + t37 * t9 + t38 * (t10 * t30 + t33) - t36 * (-t10 * t31 + t9 * t19) + t35 * (t10 * pkin(2) + t9 * pkin(6) + t27) + t34 * t10) * g(1);
U = t1;
