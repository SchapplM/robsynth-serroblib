% Calculate potential energy for
% S5RPRRP4
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
% Datum: 2019-12-05 18:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRRP4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP4_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:05:24
% EndTime: 2019-12-05 18:05:25
% DurationCPUTime: 0.51s
% Computational Cost: add. (160->70), mult. (192->65), div. (0->0), fcn. (168->8), ass. (0->28)
t19 = sin(qJ(3));
t21 = cos(qJ(3));
t44 = -m(4) * pkin(2) - t21 * mrSges(4,1) + t19 * mrSges(4,2) - mrSges(3,1);
t43 = mrSges(3,2) - mrSges(5,3) - mrSges(6,3);
t42 = -m(1) - m(2);
t41 = -mrSges(5,1) - mrSges(6,1);
t40 = mrSges(5,2) + mrSges(6,2);
t39 = -m(3) - m(4) - m(5) - m(6);
t38 = m(4) * pkin(6) + mrSges(4,3);
t17 = sin(pkin(8));
t18 = cos(pkin(8));
t37 = t43 * t17 + t44 * t18 - mrSges(2,1);
t35 = t19 * pkin(3);
t16 = qJ(3) + qJ(4);
t8 = sin(t16);
t36 = -m(5) * t35 - m(6) * (pkin(4) * t8 + t35) + mrSges(2,2) - mrSges(3,3) - t19 * mrSges(4,1) - t21 * mrSges(4,2);
t23 = -pkin(7) - pkin(6);
t7 = t21 * pkin(3) + pkin(2);
t20 = sin(qJ(1));
t32 = t18 * t20;
t22 = cos(qJ(1));
t31 = t18 * t22;
t14 = -qJ(5) + t23;
t9 = cos(t16);
t5 = pkin(4) * t9 + t7;
t28 = -t14 * t17 + t18 * t5;
t27 = -t17 * t23 + t18 * t7;
t1 = (-mrSges(1,3) + t42 * r_base(3) + t41 * (t20 * t8 + t9 * t31) - t40 * (t20 * t9 - t8 * t31) + t39 * (t22 * pkin(1) + t20 * qJ(2) + r_base(3)) + t36 * t20 + (-m(5) * t27 - m(6) * t28 - t38 * t17 + t37) * t22) * g(3) + (-mrSges(1,2) + t42 * r_base(2) + t41 * (t22 * t8 - t9 * t32) - t40 * (t22 * t9 + t8 * t32) + t39 * (t22 * qJ(2) + r_base(2)) + t36 * t22 + (m(3) * pkin(1) - m(4) * (-pkin(6) * t17 - pkin(1)) + t17 * mrSges(4,3) - m(5) * (-pkin(1) - t27) - m(6) * (-pkin(1) - t28) - t37) * t20) * g(2) + (-m(1) * r_base(1) - mrSges(1,1) - mrSges(2,3) + (-m(2) + t39) * (pkin(5) + r_base(1)) + (-m(5) * t23 - m(6) * t14 + t38 - t43) * t18 + (-m(5) * t7 - m(6) * t5 + t40 * t8 + t41 * t9 + t44) * t17) * g(1);
U = t1;
