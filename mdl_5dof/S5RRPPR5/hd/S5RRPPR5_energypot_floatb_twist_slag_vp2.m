% Calculate potential energy for
% S5RRPPR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPPR5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPPR5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR5_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR5_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:28:48
% EndTime: 2019-12-31 19:28:49
% DurationCPUTime: 0.39s
% Computational Cost: add. (155->62), mult. (156->50), div. (0->0), fcn. (126->8), ass. (0->27)
t45 = mrSges(4,2) - mrSges(5,3);
t13 = qJ(2) + pkin(8);
t10 = sin(t13);
t11 = cos(t13);
t44 = pkin(3) * t11 + qJ(4) * t10;
t43 = -m(6) * pkin(4) - mrSges(4,1) - mrSges(5,1);
t41 = -m(2) - m(3);
t40 = -m(5) - m(6);
t21 = cos(qJ(1));
t39 = t44 * t21;
t38 = -m(1) + t41;
t17 = sin(qJ(2));
t20 = cos(qJ(2));
t16 = sin(qJ(5));
t19 = cos(qJ(5));
t26 = t10 * t16 + t11 * t19;
t27 = t10 * t19 - t11 * t16;
t37 = -m(3) * pkin(1) - t20 * mrSges(3,1) - mrSges(6,1) * t26 + t17 * mrSges(3,2) - mrSges(6,2) * t27 + t45 * t10 + t43 * t11 - mrSges(2,1);
t36 = -m(3) * pkin(6) + mrSges(2,2) - mrSges(5,2) - mrSges(3,3) - mrSges(4,3) + mrSges(6,3);
t14 = pkin(5) + r_base(3);
t9 = t20 * pkin(2) + pkin(1);
t33 = t21 * t9 + r_base(1);
t31 = t17 * pkin(2) + t14;
t15 = -qJ(3) - pkin(6);
t18 = sin(qJ(1));
t25 = -t18 * t15 + t33;
t1 = (-m(1) * r_base(3) - m(4) * t31 - t17 * mrSges(3,1) - t27 * mrSges(6,1) - t20 * mrSges(3,2) + t26 * mrSges(6,2) - mrSges(1,3) - mrSges(2,3) + t40 * (t10 * pkin(3) - qJ(4) * t11 + t31) + t41 * t14 - t45 * t11 + t43 * t10) * g(3) + (-mrSges(1,2) + t38 * r_base(2) + (-m(6) * pkin(7) - t36) * t21 + (-m(4) + t40) * (t21 * t15 + t18 * t9 + r_base(2)) + (t40 * t44 + t37) * t18) * g(2) + (-mrSges(1,1) - m(4) * t25 - m(5) * (t25 + t39) - m(6) * (t33 + t39) + t38 * r_base(1) + t37 * t21 + (-m(6) * (-pkin(7) - t15) + t36) * t18) * g(1);
U = t1;
