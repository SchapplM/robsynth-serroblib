% Calculate potential energy for
% S5RPRRP3
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
% m [6x1]
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
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRRP3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP3_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:29:29
% EndTime: 2022-01-23 09:29:29
% DurationCPUTime: 0.29s
% Computational Cost: add. (147->53), mult. (109->37), div. (0->0), fcn. (73->8), ass. (0->22)
t30 = -m(5) - m(6);
t34 = mrSges(5,2) + mrSges(6,2);
t33 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t32 = -m(1) - m(2);
t31 = -m(3) - m(4);
t29 = t30 + t31;
t15 = sin(qJ(3));
t17 = cos(qJ(3));
t14 = qJ(3) + qJ(4);
t6 = sin(t14);
t7 = cos(t14);
t28 = -m(4) * pkin(2) - t17 * mrSges(4,1) + t15 * mrSges(4,2) - mrSges(3,1) + t30 * (t17 * pkin(3) + pkin(2)) + t33 * t7 + t34 * t6;
t19 = -pkin(7) - pkin(6);
t27 = m(4) * pkin(6) - m(5) * t19 - m(6) * (-qJ(5) + t19) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t26 = pkin(5) + r_base(3);
t5 = qJ(2) + t26;
t18 = cos(qJ(1));
t16 = sin(qJ(1));
t13 = qJ(1) + pkin(8);
t4 = cos(t13);
t3 = sin(t13);
t1 = (-m(1) * r_base(3) - m(2) * t26 - mrSges(4,1) * t15 - mrSges(4,2) * t17 - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - t34 * t7 + t30 * (t15 * pkin(3) + t5) + t33 * t6 + t31 * t5) * g(3) + (-t16 * mrSges(2,1) - mrSges(2,2) * t18 - mrSges(1,2) + t32 * r_base(2) + t29 * (t16 * pkin(1) + r_base(2)) + t27 * t4 + t28 * t3) * g(2) + (-mrSges(2,1) * t18 + t16 * mrSges(2,2) - mrSges(1,1) + t32 * r_base(1) + t28 * t4 + t29 * (t18 * pkin(1) + r_base(1)) - t27 * t3) * g(1);
U = t1;
