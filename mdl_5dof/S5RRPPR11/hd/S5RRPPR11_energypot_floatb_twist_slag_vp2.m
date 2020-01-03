% Calculate potential energy for
% S5RRPPR11
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPPR11_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPPR11_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR11_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR11_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR11_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:46:19
% EndTime: 2019-12-31 19:46:20
% DurationCPUTime: 0.49s
% Computational Cost: add. (134->71), mult. (184->63), div. (0->0), fcn. (156->8), ass. (0->30)
t14 = sin(pkin(8));
t12 = pkin(8) + qJ(5);
t6 = sin(t12);
t7 = cos(t12);
t51 = -m(6) * pkin(4) * t14 - t6 * mrSges(6,1) - t7 * mrSges(6,2) + mrSges(3,2) - mrSges(4,3);
t50 = -mrSges(3,1) + mrSges(4,2) + m(6) * (-pkin(7) - qJ(4)) - mrSges(6,3) - m(5) * qJ(4);
t48 = -m(1) - m(2);
t47 = -m(5) - m(6);
t18 = sin(qJ(1));
t17 = sin(qJ(2));
t33 = qJ(3) * t17;
t19 = cos(qJ(2));
t37 = t18 * t19;
t46 = pkin(2) * t37 + t18 * t33;
t45 = m(4) - t47;
t42 = t7 * mrSges(6,1) - t6 * mrSges(6,2) + mrSges(4,1) - mrSges(2,2) + mrSges(3,3);
t41 = t17 * t51 + t19 * t50 - mrSges(2,1);
t39 = t18 * t14;
t15 = cos(pkin(8));
t38 = t18 * t15;
t20 = cos(qJ(1));
t36 = t20 * t14;
t35 = t20 * t15;
t34 = t20 * t19;
t13 = pkin(5) + r_base(3);
t31 = pkin(1) * t18 + r_base(2);
t29 = pkin(1) * t20 + pkin(6) * t18 + r_base(1);
t22 = -t20 * pkin(6) + t31;
t5 = pkin(4) * t15 + pkin(3);
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * t13 - t45 * (pkin(2) * t17 + t13) + (t14 * mrSges(5,1) + t15 * mrSges(5,2) + qJ(3) * t45 - t51) * t19 + (-mrSges(5,3) + t50) * t17) * g(3) + (-mrSges(1,2) - m(3) * t22 - m(4) * (t22 + t46) - (t17 * t39 - t35) * mrSges(5,1) - (t17 * t38 + t36) * mrSges(5,2) - mrSges(5,3) * t37 + t48 * r_base(2) + t47 * (t31 + t46) + (-m(5) * (-pkin(3) - pkin(6)) - m(6) * (-pkin(6) - t5) + t42) * t20 + t41 * t18) * g(2) + (-mrSges(1,1) - m(3) * t29 - (t17 * t36 + t38) * mrSges(5,1) - (t17 * t35 - t39) * mrSges(5,2) - mrSges(5,3) * t34 + t48 * r_base(1) - t45 * (pkin(2) * t34 + t20 * t33 + t29) + t41 * t20 + (-m(5) * pkin(3) - m(6) * t5 - t42) * t18) * g(1);
U = t1;
